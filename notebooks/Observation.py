# %%
from astropy.table import Table
import warnings
warnings.filterwarnings("ignore")
from ipywidgets import Button, Layout, jslink, IntText, IntSlider, interactive, interact, HBox, Layout, VBox
from astropy.modeling.functional_models import Gaussian2D, Gaussian1D
import os
from IPython.display import display, clear_output
import ipywidgets as widgets
import matplotlib.pyplot as plt
import inspect
from astropy.io import fits
import numpy as np
from functools import wraps
import inspect
from scipy.sparse import dia_matrix
from scipy.interpolate import interpn
from scipy.special import erf
from scipy import special
from scipy.ndimage import gaussian_filter1d, gaussian_filter
import pandas as pd
import functools
import re
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy.visualization import quantity_support
# from astropy.visualization import ZScaleInterval
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from astropy.wcs import WCS
import numpy as np
from astropy.cosmology import Planck18 as cosmo

np.seterr(invalid='ignore')
 

from scipy.ndimage import map_coordinates





def download(url, file=""):
    """Download a file
    """
    from tqdm import tqdm  # , tqdm_gui
    import requests

    try:
        response = requests.get(url, stream=True)
    except requests.exceptions.RequestException as e:
        verboseprint(e)
        return False
    else:
        total_size_in_bytes = int(response.headers.get("content-length", 0))
        block_size = 1024  # 1 Kibibyte
        progress_bar = tqdm(
            total=0.95 * total_size_in_bytes, unit="iB", unit_scale=True
        )
        with open(file, "wb") as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
        # progress_bar.close()
        # tqdm_gui.close(progress_bar)
        # progress_bar.display()
        # progress_bar.plt.close(progress_bar.fig)
        # plt.show(block=False)
        # plt.close('all')
        # plt.close(progress_bar.fig)
        if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
            verboseprint("ERROR, something went wrong")
            return False
        else:
            return True


def create_wcs_header(wave_range_nm, spatial_extent_arcsec, output_shape):
    """Crée un header WCS avec affichage correct du slicer wavelength dans DS9."""
    from astropy.wcs import WCS

    wcs = WCS(naxis=3)
    wcs.wcs.crpix = [output_shape[1] // 2 + 1, output_shape[2] // 2 + 1, 1]
    wcs.wcs.cdelt = [spatial_extent_arcsec / 3600, spatial_extent_arcsec / 3600, (wave_range_nm[1] - wave_range_nm[0]) / output_shape[0]]
    wcs.wcs.crval = [0, 0, wave_range_nm[0]]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    wcs.wcs.cunit = ["deg", "deg", "nm"]
    wcs.wcs.specsys = "TOPOCENT"
    wcs.wcs.restfrq = 0  # Ajout pour lecture correcte en spectre
    wcs.wcs.restwav = (wave_range_nm[0] + wave_range_nm[1]) / 2  # Onde centrale
    
    header = wcs.to_header()
    header['SPECSYS'] = 'TOPOCENT'
    header['WCSNAME'] = 'SIMULATED'
    header['BUNIT'] = 'erg/cm2/s/arcsec2/A'
    header['OBJECT'] = 'Simulated Cube'
    
    return header


def convert_fits_cube(input_fits, output_fits, output_shape, wave_range_nm, spatial_extent_arcsec,redshift=0.7,phys=False):
    """Convertit un cube FITS avec WCS basé sur wave_range (nm) et spatial_extent (arcsec)."""
    input_cube = fits.open(input_fits)[0]
    
    resampled_cube = resample_cube(input_cube, wave_range_nm, spatial_extent_arcsec, output_shape,phys=phys)
    resampled_cube = convert_LU2ergs(resampled_cube, np.mean(np.array(wave_range_nm)))

    header = create_wcs_header(wave_range_nm, spatial_extent_arcsec, output_shape)

    hdu = fits.PrimaryHDU(data=resampled_cube, header=header)
    hdu.writeto(output_fits, overwrite=True)
    return output_fits


def convert_LU2ergs(LU,wave_nm): #TODO here it should not be 200 nm but 1216 so we need indeed to account for the Redshift!!!
    wave =wave_nm * 1e-7 #/ (1+Redshift) It gives the good result from mature paper : at z =fwhm_sigma_ratio, 1 LU/˚A is 1 ph cm-2 s-1 sr-1˚A-1 = 1.2 ˆ 10-22 erg cm-2 s-1˚A-1 
    Energy = 6.62e-27 * 3e10 / wave
    angle =  np.pi / (180 * 3600)
    flux_ergs = LU * Energy * angle * angle
    return flux_ergs

def convert_ergs2LU(flux_ergs,wave_nm):
    wave =wave_nm * 1e-7 #/ (1+Redshift)
    Energy = 6.62e-27 * 3e10 / wave
    angle =    np.pi / (180 * 3600) 
    LU = flux_ergs/ (Energy  * angle * angle)
    flux_ergs = LU * Energy * angle * angle
    return LU

    # solid_angle = (np.pi / (180 * 3600))**2  # Convert arcsec² to steradians    
    # LU = flux_ergs / (Energy * solid_angle)  # Convert ergs to LU


def resample_cube(input_cube, wave_range_nm, spatial_extent_arcsec, output_shape, phys=False):
    """
    Resample un datacube avec interpolation/extrapolation selon longueur d'onde et spatial.
    """
    header = input_cube.header
    input_cube=input_cube.data
    wave_start_out, wave_end_out = wave_range_nm
    input_nz, input_nx, input_ny = input_cube.data.shape
    # print(input_cube.data.shape)
    cunit3 = header.get('CUNIT3', 'nm').strip().lower()
    # Définition du facteur de conversion
    unit_conversion = {
        'm': 1e9,       # mètres → nanomètres
        'cm': 1e7,      # centimètres → nanomètres
        'mm': 1e6,      # millimètres → nanomètres
        'um': 1e3,      # micromètres → nanomètres
        'nm': 1,        # nanomètres (aucun changement)
        'angstrom': 0.1 # Ångströms → nanomètres
    }
    factor = unit_conversion.get(cunit3, 1)  # Valeur par défaut = 1 (nm)
    # Calcul correct des longueurs d'onde en tenant compte de l'unité
    wave_start_in = (header['CRVAL3'] + (1 - header['CRPIX3']) * header['CDELT3']) * factor
    wave_end_in = (header['CRVAL3'] + (input_nz - header['CRPIX3']) * header['CDELT3']) * factor
    spatial_radius_in_x = header['CDELT1'] * (input_nx - 1) / 2.0 * 3600  # Convert from degrees to arcseconds
    spatial_radius_in_y = header['CDELT2'] * (input_ny - 1) / 2.0 * 3600  # Convert from degrees to arcseconds
    spatial_radius_out = spatial_extent_arcsec / 2.0


    if phys==False:
        # print("phys",phys,False)
        input_nx, input_ny, input_nz = input_cube.shape
        x_out = np.linspace(0, input_nx - 1, output_shape[0])
        y_out = np.linspace(0, input_ny - 1, output_shape[1])
        z_out = np.linspace(0, input_nz - 1, output_shape[2])

        x_grid, y_grid, z_grid = np.meshgrid(x_out, y_out, z_out, indexing='ij')
        resampled_cube = map_coordinates(input_cube, [x_grid, y_grid, z_grid], order=1, mode='nearest')
        return resampled_cube
    else:
        # print("phys",phys,True)
        # if we are out of limit, we change the reshift of the cube to end in the limit:
        # wave_start_in, wave_end_in = 203,205
        # wave_start_out, wave_end_out = 400, 500
        input_wave_center = (wave_start_in+ wave_end_in)/2
        inst_wave_center = (wave_start_out+ wave_end_out)/2
        if 1>0 : #(input_wave_center < wave_start_out ) | (input_wave_center> wave_end_out):
            ratio = inst_wave_center/input_wave_center
            # print(ratio)
            # = 121.5
            redshift = inst_wave_center/(input_wave_center/(1+0.7))-1
            # print(redshift)
            wave_start_in, wave_end_in = wave_start_in*ratio, wave_end_in*ratio
        input_wave = np.linspace(wave_start_in, wave_end_in, input_nz)
        # print(input_wave)
        input_x = np.linspace(-abs(spatial_radius_in_x), abs(spatial_radius_in_x), input_nx)
        input_y = np.linspace(-abs(spatial_radius_in_y),abs( spatial_radius_in_y), input_ny)

        output_wave = np.linspace(wave_start_out, wave_end_out, output_shape[0])
        # print(wave_start_in, wave_end_in, wave_start_out, wave_end_out)
        output_x = np.linspace(-spatial_radius_out, spatial_radius_out, output_shape[2])
        output_y = np.linspace(-spatial_radius_out, spatial_radius_out, output_shape[1])        

        interpolator = RegularGridInterpolator(
            (input_wave, input_x, input_y),
            input_cube,
            method='linear',
            bounds_error=False,
            fill_value=None
        )
        W, X, Y = np.meshgrid(output_wave, output_x, output_y, indexing='ij')
        resampled_cube = interpolator((W, X, Y))
        if np.nanmin(resampled_cube)<0:
            resampled_cube -= np.nanmin(resampled_cube)
        scale_factor = (
            (output_x[1] - output_x[0]) / (input_x[1] - input_x[0]) *
            (output_y[1] - output_y[0]) / (input_y[1] - input_y[0]) *
            (output_wave[1] - output_wave[0]) / (input_wave[1] - input_wave[0])
        )
        # print(np.max(resampled_cube),scale_factor)
        return resampled_cube * scale_factor
        # TODO est ce des pixels ou autre chose ?? si c'est des pixels il faut multiplier par le rapport de la taille des pixels



def to_float_or_nan(arr):
    def convert(val):
        try:
            return float(val)  # Try to convert to float
        except (ValueError, TypeError):  # If conversion fails, return NaN
            return np.nan

    # Apply conversion to entire array using list comprehension
    return np.array([convert(val) for val in arr])







def generate_gaussian_galaxy(amplitude, redshift, platescale, PSF_RMS, size_kpc):

    # Define the final image size
    final_image_size = [500, 100]  # Number of pixels per side of the output image
    # Convert the physical size of the galaxy to angular size in arcseconds
    angular_size_arcsec = size_kpc / cosmo.kpc_proper_per_arcmin(redshift).to(u.kpc / u.arcsec).value
    # Create a grid of coordinates
    x = np.linspace(-platescale * final_image_size[0] / 2, platescale * final_image_size[0] / 2, final_image_size[0])
    y = np.linspace(-platescale * final_image_size[1] / 2, platescale * final_image_size[1] / 2, final_image_size[1])
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2 + yy**2)  # Radial distance from the center
    # Create a Gaussian profile for the galaxy
    galaxy = np.exp(-r**2 / (2 * (angular_size_arcsec /2.355)**2))  # Gaussian profile
    # Normalize the galaxy to conserve total energy
    galaxy_sum_initial = np.max(galaxy) * 100  # Approximate normalization factor
    galaxy *= amplitude / galaxy_sum_initial  # Normalize so that the peak is equal to the amplitude
    galaxy *= size_kpc * size_kpc  # Scale to maintain the total energy constant
    # Convolve the galaxy with the PSF
    PSF_RMS_pix = np.sqrt(PSF_RMS**2 + 0.05**2) / platescale  # Convert PSF RMS to pixels
    galaxy_convolved = gaussian_filter(galaxy, PSF_RMS_pix)  # Apply Gaussian smoothing
    return galaxy_convolved

def float_to_latex(mumber):
    try:
        return "$"+ ("%.1E"%(mumber)).replace("E"," 10^{")+"}$"
    except TypeError as e:
        print(e,mumber)
    # return ("%.1E"%(mumber)).replace("E"," 10$^{")+"}$"

def rsetattr(obj, attr, val):
    pre, _, post = attr.rpartition('.')
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)


def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)
    return functools.reduce(_getattr, [obj] + attr.split('.'))



def generate_spiral_galaxy(amplitude,redshift, pixel_scale, PSF_RMS,size_kpc,kpc=True,fwhm_arcsec=None, final_image_size = [500,100]):
    """
    Generates a simulated galaxy image of defined size and PSF.
    
    Parameters:
    pixel_scale (float): Redshift of the galaxy.
    platescale (float): Size of one pixel in arcsec/pix.
    PSF_RMS (float): RMS size of the PSF in arcsec.
    
    Returns:
    np.ndarray: Simulated and convolved 2D galaxy image.
    """
    # size_kpc = 10.0  # Physical size of the galaxy in kpc
    intensity = 0.4  # Number of pixels per side of the output image
    if kpc:
        angular_size_arcsec = size_kpc / cosmo.kpc_proper_per_arcmin(redshift).to(u.kpc / u.arcsec).value
        angular_size_pix = angular_size_arcsec / pixel_scale
    else:
        angular_size_pix = fwhm_arcsec / pixel_scale  # Conversion FWHM -> sigma en pixels

    # field_of_view_arcsec = 
    x = np.linspace(-pixel_scale * final_image_size[0] / 2, pixel_scale * final_image_size[0] / 2, final_image_size[0])
    y = np.linspace(-pixel_scale * final_image_size[1] / 2, pixel_scale * final_image_size[1] / 2, final_image_size[1])
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2 + yy**2)
    theta = np.arctan2(yy, xx)
    
    core = np.exp(-r**2 / (2 * (angular_size_pix / 5 /2.355)**2))
    spiral = np.exp(-r**2 / (2 * (angular_size_pix /2.355)**2)) * (1 + intensity * np.sin(2 * theta + r * 20 / angular_size_pix))
    galaxy = core + spiral

    # Normalize the galaxy to conserve total energy
    galaxy_sum_initial = np.max(galaxy)*100
    galaxy *= amplitude /galaxy_sum_initial  # Normalize so that the sum is 1
    if kpc:
        galaxy *= size_kpc * size_kpc  # Scale to maintain the total energy constant

    PSF_RMS_pix = np.sqrt(PSF_RMS**2 + 0.05**2) / pixel_scale
    galaxy_convolved = gaussian_filter(galaxy, PSF_RMS_pix)

    return galaxy_convolved

def generate_filament(amplitude, redshift, pixel_scale, PSF_RMS, size_kpc, size_arcsec=None,final_image_size=[500, 100]):
    """
    Generates a simulated filament image that runs diagonally (x = y).
    """
    if size_arcsec is None:
        angular_size_arcsec = size_kpc / cosmo.kpc_proper_per_arcmin(redshift).to(u.kpc / u.arcsec).value
    else:
        angular_size_arcsec = size_arcsec
    x = np.linspace(-pixel_scale * final_image_size[0] / 2, pixel_scale * final_image_size[0] / 2, final_image_size[0])
    y = np.linspace(-pixel_scale * final_image_size[1] / 2, pixel_scale * final_image_size[1] / 2, final_image_size[1])
    xx, yy = np.meshgrid(x, y)
    diagonal = (xx - yy) / np.sqrt(2)
    filament = np.exp(-diagonal**2 / (2 * (angular_size_arcsec /2.35)**2))
    filament_sum_initial = np.max(filament) * 100
    filament *= amplitude / filament_sum_initial
    if size_kpc is not None:
        filament *= size_kpc * size_kpc
    PSF_RMS_pix = np.sqrt(PSF_RMS**2 + 0.05**2) / pixel_scale
    filament_convolved = gaussian_filter(filament, PSF_RMS_pix)
    return filament_convolved


def create_fits_cube(output_path="/tmp/test.fits", cube_size=(100, 100, 500), pixel_scale=2, wave_range=(195, 215), continuum_flux=None, continuum_fwhm=2, line_flux=1e-19, line_center=206, line_fwhm=10, line_spatial_fwhm=7, galaxy_shape=False, kpc_size=None, redshift=None,filament=False):
    wavelengths = np.linspace(wave_range[0]/10, wave_range[1]/10, cube_size[2])
    dw = np.mean(np.diff(wavelengths))
    cube = np.zeros(cube_size, dtype=np.float32)
    
    if galaxy_shape:
        # if kpc_size is None:
        continuum_profile = generate_spiral_galaxy(amplitude=np.max(continuum_flux), redshift=redshift, pixel_scale=pixel_scale, PSF_RMS=0.1,fwhm_arcsec=None, size_kpc=kpc_size, kpc=(kpc_size is not None), final_image_size=[cube_size[0], cube_size[1]])
        # fitswrite(continuum_profile,"/tmp/test_image.fits")
        # else:
        #     continuum_profile = generate_spiral_galaxy(amplitude=np.max(continuum_flux), redshift=redshift, pixel_scale=pixel_scale, PSF_RMS=0.1, size_kpc=kpc_size, kpc=True, final_image_size=[cube_size[0], cube_size[1]])
    else:
        sigma_continuum = continuum_fwhm / pixel_scale
        continuum_profile = np.zeros((cube_size[0], cube_size[1]))
        continuum_profile[cube_size[0]//2, cube_size[1]//2] = 1
        continuum_profile = gaussian_filter(continuum_profile, sigma_continuum)
        continuum_profile /= np.max(continuum_profile)
    cube += continuum_flux[None, None, :] * continuum_profile[:, :, None]
    # fitswrite(cube,"/tmp/testtest.fits")
    sigma_line = line_spatial_fwhm / (2.355 * pixel_scale)
    line_profile = np.zeros((cube_size[0], cube_size[1]))
    line_profile[cube_size[0]//2, cube_size[1]//2] = 1
    line_profile = gaussian_filter(line_profile, sigma_line)
    line_profile /= np.max(line_profile)
    sigma_wave = line_fwhm / (2.355 * 10)
    line_spectrum = np.exp(-0.5 * ((wavelengths - line_center) / sigma_wave) ** 2)
    line_spectrum /= np.max(line_spectrum)
    cube += line_flux * line_profile[:, :, None] * line_spectrum[None, None, :]
    if filament & galaxy_shape:
        filament_profile = generate_filament(amplitude=np.max(line_spectrum), redshift=redshift, pixel_scale=pixel_scale, PSF_RMS=0.1, size_kpc=kpc_size/2,  size_arcsec=line_spatial_fwhm/2, final_image_size=[cube_size[0], cube_size[1]])
        cube += line_flux * filament_profile[:, :, None] * line_spectrum[None, None, :]
    # fitswrite(cube,"/tmp/testtest2.fits")
    wcs = WCS(naxis=3)
    wcs.wcs.crpix = [cube_size[1] // 2, cube_size[0] // 2, 1]
    wcs.wcs.cdelt = [pixel_scale / 3600, pixel_scale / 3600, dw]
    wcs.wcs.crval = [0, 0, wave_range[0]]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    hdu = fits.PrimaryHDU(cube.T, header=wcs.to_header())
    hdu.writeto(output_path, overwrite=True)
    return cube

def initializer(func):
    """
    Automatically assigns the parameters.

    >>> class process:
    ...     @initializer
    ...     def __init__(self, cmd, reachable=False, user='root'):
    ...         pass
    >>> p = process('halt', True)
    >>> p.cmd, p.reachable, p.user
    ('halt', True, 'root')
    """
    # names, varargs, keywords, defaults = inspect.getargspec(func)
    names, varargs, keywords, defaults,_,_,_ = inspect.getfullargspec(func)

    @wraps(func)
    def wrapper(self, *args, **kargs):
        for name, arg in list(zip(names[1:], args)) + list(kargs.items()):
            setattr(self, name, arg)

        for name, default in zip(reversed(names), reversed(defaults)):
            if not hasattr(self, name):
                setattr(self, name, default)

        func(self, *args, **kargs)

    return wrapper


# Initialization of the thresholding functions. So that files are not read several times
n=10
type_="" #"new_" #""
#new is for when we don't use fraction and use RN (false I think), "" is with fraction true positives and RN/gain, seems better 
path=""
# path = "/Users/Vincent/Github/fireball2-etc/notebooks/"
table_threshold = fits.open(path+"../data/Instruments/EMCCD/%sthreshold_%s.fits"%(type_,n))[0].data
table_snr = fits.open(path+"../data/Instruments/EMCCD/%ssnr_max_%s.fits"%(type_,n))[0].data
table_fraction_rn = fits.open(path+"../data/Instruments/EMCCD/%sfraction_rn_%s.fits"%(type_,n))[0].data
table_fraction_flux = fits.open(path+"../data/Instruments/EMCCD/%sfraction_flux_%s.fits"%(type_,n))[0].data




def variable_smearing_kernels(image, smearing=1.5, SmearExpDecrement=50000):
    """Creates variable smearing kernels for inversion
    """
    import numpy as np
    
    smearing_length = smearing * np.exp(-image / SmearExpDecrement)
    smearing_kernels = np.exp(-np.arange(6)[:, np.newaxis, np.newaxis] / smearing_length)
    smearing_kernels /= smearing_kernels.sum(axis=0)
    return smearing_kernels   

#temperature=-100,
class Observation:
    @initializer
    # def __init__(self, instrument="FIREBall-2 2023", Atmosphere=0.5, Throughput=0.13*0.9, exposure_time=50, counting_mode=False, Signal=1e-16, EM_gain=1400, RN=109, CIC_charge=0.005, Dark_current=0.08, Sky=10000, readout_time=1.5, extra_background = 0,acquisition_time = 2,smearing=0,i=25,plot_=False,temperature=-100,n=n,PSF_RMS_mask=5, PSF_RMS_det=8, QE = 0.45,cosmic_ray_loss_per_sec=0.005,Size_source=16,lambda_stack=1,Slitwidth=5,Bandwidth=200,Collecting_area=1,Δx=0,Δλ=0,pixel_scale=np.nan, Spectral_resolution=np.nan, dispersion=np.nan,Line_width=np.nan,wavelength=np.nan, pixel_size=np.nan,len_xaxis=50):#,photon_kept=0.7#, flight_background_damping = 0.9
    def __init__(self, instruments=None, instrument="FIREBall-2 2025", Atmosphere=None, Throughput=None, exposure_time=None, counting_mode=False, Signal=None, EM_gain=None, RN=None, CIC_charge=None, Dark_current=None, Sky=None, readout_time=None, extra_background = None,acquisition_time = None,smearing=None,i=33,plot_=False,n=n,PSF_RMS_mask=None, PSF_RMS_det=None, QE = None,cosmic_ray_loss_per_sec=None,Size_source=None,lambda_stack=1,Slitwidth=None,Bandwidth=None,Collecting_area=None,Δx=None,Δλ=None,pixel_scale=None, Spectral_resolution=None, dispersion=None,Line_width=None,wavelength=None, pixel_size=None,len_xaxis=50,Slitlength=None,IFS=None, Redshift=None, Throughput_FWHM=None, SNR_res="per pix",spectrograph=True,test=True):#,photon_kept=0.7#, flight_background_damping = 0.9
    # def __init__(self, instrument="FIREBall-2 2023", Atmosphere=0.5, Throughput=0.13, exposure_time=50, counting_mode=False, Signal=1e-17, EM_gain=1500, RN=40, CIC_charge=0.005, Dark_current=1, Sky=2e-18, readout_time=5, extra_background = 0.5,acquisition_time = 2,smearing=1.50,i=33,plot_=False,n=n,PSF_RMS_mask=2.5, PSF_RMS_det=3, QE = 0.4,cosmic_ray_loss_per_sec=0.005,Size_source=16,lambda_stack=0.21,Slitwidth=6,Bandwidth=160,Collecting_area=0.707,Δx=0,Δλ=0,pixel_scale=1.1, Spectral_resolution=1300, dispersion=0.21,Line_width=15,wavelength=200, pixel_size=13,len_xaxis=50,Slitlength=10):#,photon_kept=0.7#, flight_background_damping = 0.9
        """
        ETC calculator: computes the noise budget at the detector level based on instrument/detector parameters
        This is currently optimized for slit spectrographs and EMCCD but could be pretty easily generalized to other instrument type if needed
        """
        self.instruments = instruments
        self.instruments_dict = {name: {key: val for key, val in zip(instruments["Charact."][:], instruments[name][:]) if not isinstance(key, np.ma.core.MaskedConstant) and not isinstance(val, np.ma.core.MaskedConstant)} for name in instruments.colnames[3:]}

        for key, val in zip(instruments["Charact."][:],instruments[instrument][:]):
            if hasattr(self,key):
                if getattr(self,key) is None:
                    setattr(self, key, float(val))

        self.initilize(IFS=IFS)
    
    # TODO delete all the *= that actually multiply again the signal
    def initilize(self,IFS):
        self.precise = False
        # self.Signal = Gaussian2D(amplitude=self.Signal,x_mean=0,y_mean=0,x_stddev=self.Size_source,y_stddev=self.Line_width,theta=0)(self.Δx,self.Δλ)
        # print("\ni",self.i,"\nAtmosphere",self.Atmosphere, "\nThroughput=",self.Throughput,"\nSky=",self.Sky, "\nacquisition_time=",self.acquisition_time,"\ncounting_mode=",self.counting_mode,"\nSignal=",self.Signal,"\nEM_gain=",self.EM_gain,"RN=",self.RN,"CIC_charge=",self.CIC_charge,"Dark_current=",self.Dark_current,"\nreadout_time=",self.readout_time,"\nsmearing=",self.smearing,"\nextra_background=",self.extra_background,"\nPSF_RMS_mask=",self.PSF_RMS_mask,"\nPSF_RMS_det=",self.PSF_RMS_det,"\nQE=",self.QE,"\ncosmic_ray_loss_per_sec=",self.cosmic_ray_loss_per_sec,"\nlambda_stack",self.lambda_stack,"\nSlitwidth",self.Slitwidth, "\nBandwidth",self.Bandwidth,"\nSize_source",self.Size_source,"\nCollecting_area",self.Collecting_area)
        # print("\Collecting_area",self.Collecting_area, "\nΔx=",self.Δx,"\nΔλ=",self.Δλ, "\napixel_scale=",self.pixel_scale,"\nSpectral_resolution=",self.Spectral_resolution,"\ndispersion=",self.dispersion,"\nLine_width=",self.Line_width,"wavelength=",self.wavelength,"pixel_size=",self.pixel_size)
      
      
        self.spectro = self.spectrograph#False if np.isnan(self.instruments_dict[self.instrument]["dispersion"]) else True
        # Simple hack to me able to use UV magnitudes (not used for the ETC)
        if np.max([self.Signal])>1:
            self.Signal = 10**(-(self.Signal-20.08)/2.5)*2.06*1E-16
        #TODO be sure we account for potentialfwhm_sigma_ratio ratio here
        #convolve input flux by instrument PSF
        if type(self.Slitlength) == np.float64:
            if (self.Slitlength==self.Slitwidth):
                self.Signal *= np.pi/4 # ratio between fiber disk and square slit
        
        if self.precise: # TODO are we sure we should do that here?
            self.Signal *= (erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) )
            #convolve input flux by spectral resolution
            self.spectro_resolution_A = 10*self.wavelength/self.Spectral_resolution
            if self.spectro:
                self.Signal *= (erf(self.Line_width / (2 * np.sqrt(2) * self.spectro_resolution_A  )) )
            # print("Factor spatial and spectral",  (erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) ),   (erf(self.Line_width / (2 * np.sqrt(2) * 10*self.wavelength/self.Spectral_resolution)) ))
        # print(self.IFS)
        if ~np.isnan(self.Slitwidth).all(): #& (self.precise) # & (self.SNR_res!="per Source")
            # assess flux fraction going through slit
            # self.flux_fraction = ((1+erf(self.Slitwidth/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)     *    ((1+erf(self.Slitlength/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)
            self.flux_fraction =   ((1+erf(self.Slitlength/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)
            if    (~self.IFS): 
                self.flux_fraction *=((1+erf(self.Slitwidth/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)  
        else:
            self.flux_fraction = 1
        self.flux_fraction_slit_applied = self.flux_fraction



        # if self.smearing>0:
        # self.Signal *= 1 - np.exp(-1/(self.smearing+1e-15)) - np.exp(-2/(self.smearing+1e-15))  - np.exp(-3/(self.smearing+1e-15))
        
        
        self.PSF_lambda_pix = 10*self.wavelength / self.Spectral_resolution / self.dispersion
        fwhm_sigma_ratio =1.0 # 2.355
        if self.spectro:
            source_spatial_pixels = np.maximum(1,np.minimum(np.sqrt(self.Size_source**2+self.PSF_RMS_mask**2) * fwhm_sigma_ratio / self.pixel_scale, self.Slitwidth / self.pixel_scale))
            # source_spatial_pixels = np.minimum(self.Size_source * fwhm_sigma_ratio / self.pixel_scale, self.Slitwidth / self.pixel_scale)
            source_spectral_pixels = np.maximum(1, np.sqrt(self.PSF_lambda_pix**2 + (np.minimum(self.Line_width, self.Bandwidth) / self.dispersion)**2))
            self.source_size = np.maximum(np.minimum(self.Size_source * fwhm_sigma_ratio, self.Slitlength) / self.pixel_scale,1)    * np.sqrt(source_spatial_pixels**2 + source_spectral_pixels**2)
            self.pixels_total_source =  self.source_size  * ( np.ceil(np.sqrt(self.Size_source**2+self.PSF_RMS_mask**2)*fwhm_sigma_ratio / self.Slitwidth) if self.IFS else 1)
            # source_spectral_pixels = np.sqrt(    (np.minimum(self.Line_width, self.Bandwidth) / self.dispersion)**2    + (self.pixel_scale)**2)
            # self.source_size = np.minimum(self.Size_source * fwhm_sigma_ratio, self.Slitlength) / self.pixel_scale * np.sqrt(source_spatial_pixels**2 + source_spectral_pixels**2)
            # self.pixels_total_source = self.source_size * (np.ceil(self.Size_source * fwhm_sigma_ratio / self.Slitwidth) if self.IFS else 1)
# np.minimum(self.PSF_RMS_mask *fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale     *  np.sqrt(np.minimum(self.PSF_RMS_det/self.pixel_scale,self.PSF_lambda_pix)**2+self.PSF_lambda_pix**2)
            # TODO this one gives the same result for tot and not tot but I think we should replace PSF_RMS_mask by PSF_RMS_det
            self.elem_size =  np.ceil(np.minimum(self.PSF_RMS_mask *fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale)      *  np.ceil(np.sqrt(np.minimum(self.PSF_RMS_det/self.pixel_scale,self.PSF_lambda_pix)**2+np.minimum(self.Line_width/self.dispersion,self.PSF_lambda_pix)**2))   * ( np.ceil(self.PSF_RMS_mask*fwhm_sigma_ratio/ self.Slitwidth) if self.IFS else 1)
            # self.elem_size =  np.maximum(np.minimum(self.PSF_RMS_mask *fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale,1)      *  np.maximum(np.sqrt(np.minimum(self.PSF_RMS_det/self.pixel_scale,self.PSF_lambda_pix)**2+self.PSF_lambda_pix**2),1)   * ( np.ceil(self.PSF_RMS_mask*fwhm_sigma_ratio/ self.Slitwidth) if self.IFS else 1)
            # self.elem_size =  np.ceil(np.minimum(self.PSF_RMS_det *fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale)         *    np.ceil(np.sqrt((self.PSF_RMS_mask/self.pixel_scale)**2+self.PSF_lambda_pix**2))                                                                                    * ( np.ceil(self.PSF_RMS_mask*fwhm_sigma_ratio/ self.Slitwidth) if self.IFS else 1)
            # print( np.ceil(np.minimum(self.PSF_RMS_det *fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale), np.ceil(np.sqrt((self.PSF_RMS_mask/self.pixel_scale)**2+self.PSF_lambda_pix**2))   , ( np.ceil(self.PSF_RMS_mask*fwhm_sigma_ratio/ self.Slitwidth) if self.IFS else 1))
            # print(self.PSF_RMS_mask *fwhm_sigma_ratio,self.Slitlength ,self.pixel_scale)
        else:
            self.source_size =  (self.Size_source *fwhm_sigma_ratio /self.pixel_scale) **2
            self.elem_size = (self.PSF_RMS_det *fwhm_sigma_ratio /self.pixel_scale) **2

        if self.SNR_res=="per pix":
          self.number_pixels_used = 1
        elif self.SNR_res=="per Res elem": # is that true ? when not IFS, the SNR won't get bigger than the slit , the rest will be cut
            self.number_pixels_used = np.ceil(self.elem_size)
        elif self.SNR_res=="per Source":
            self.number_pixels_used = np.ceil(self.pixels_total_source)

        # print(self.IFS,self.elem_size,self.source_size, self.pixels_total_source)
        # print(np.minimum(self.PSF_RMS_mask *fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale, np.sqrt(np.minimum(self.PSF_RMS_det/self.pixel_scale,self.PSF_lambda_pix)**2+np.minimum(self.Line_width/self.dispersion,self.PSF_lambda_pix)**2),( np.ceil(self.PSF_RMS_mask*fwhm_sigma_ratio/ self.Slitwidth) if self.IFS else 1),     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength) /self.pixel_scale      , np.sqrt(np.minimum(self.Slitwidth/self.pixel_scale,self.Bandwidth/self.dispersion)**2+(np.minimum(self.Line_width,self.Bandwidth)/self.dispersion)**2)  , (np.ceil(self.Size_source*fwhm_sigma_ratio / self.Slitwidth) if self.IFS else 1) , )

        red, blue, violet, yellow, green, pink, grey  = '#E24A33','#348ABD','#988ED5','#FBC15E','#8EBA42','#FFB5B8','#777777'
        self.colors= [red, violet, yellow  ,blue, green, pink, grey ]


        self.ENF = 1 if (self.counting_mode | ((self.EM_gain*np.ones(self.len_xaxis))[self.i]<2)) else 2 # Excess Noise Factor 

        self.CIC_noise = np.sqrt(self.CIC_charge * self.ENF) 
        self.Dark_current_f = self.Dark_current * self.exposure_time / 3600 # e/pix/frame
        self.Dark_current_noise =  np.sqrt(self.Dark_current_f * self.ENF)
        self.effective_area =  self.QE * self.Throughput * self.Atmosphere  *    (self.Collecting_area * 100 * 100)
        # For now we put the regular QE without taking into account the photon kept fracton, because then infinite loop. 
        # Two methods to compute it: interpolate_optimal_threshold & compute_optimal_threshold
        # self.pixel_scale  = (self.pixel_scale*np.pi/180/3600) #go from arcsec/pix to str/pix 
        self.arcsec2str = (np.pi/180/3600)**2
        self.Sky_CU = convert_ergs2LU(self.Sky, self.wavelength)  
        # print(self.Sky_CU,self.Sky, self.wavelength, self.Redshift)
        if (self.counting_mode) : #& (self.EM_gain>=1)  Normaly if counting mode is on EM_gain is >1
            if self.spectro:
                if  (self.SNR_res!="per Source"):
                    self.factor_CU2el =  self.effective_area * np.minimum(self.Slitwidth,self.Size_source)  * self.arcsec2str  * self.dispersion
                else:
                    self.factor_CU2el =  self.effective_area  * self.Size_source  * self.arcsec2str  * self.dispersion
            else:
                self.factor_CU2el =  self.effective_area

            self.sky = self.Sky_CU*self.factor_CU2el*self.exposure_time  # el/pix/frame
            self.Sky_noise_pre_thresholding = np.sqrt(self.sky * self.ENF) 
            self.signal_pre_thresholding = self.Signal*self.factor_CU2el*self.exposure_time  # el/pix/frame
            self.n_threshold, self.Photon_fraction_kept, self.RN_fraction_kept, self.gain_thresholding = self.interpolate_optimal_threshold(plot_=self.plot_, i=self.i)#,flux=self.signal_pre_thresholding)
            # self.n_threshold, self.Photon_fraction_kept, self.RN_fraction_kept, self.gain_thresholding = self.compute_optimal_threshold(plot_=plot_, i=i,flux=self.signal_pre_thresholding)
        else:
            self.n_threshold, self.Photon_fraction_kept, self.RN_fraction_kept, self.gain_thresholding = np.zeros(self.len_xaxis),np.ones(self.len_xaxis),np.ones(self.len_xaxis), np.zeros(self.len_xaxis)
                
        # The faction of detector lost by cosmic ray masking (taking into account ~5-10 impact per seconds and around 2000 pixels loss per impact (0.01%))
        self.cosmic_ray_loss = np.minimum(self.cosmic_ray_loss_per_sec*(self.exposure_time+self.readout_time/2),1)
        self.QE_efficiency = self.Photon_fraction_kept * self.QE
        self.effective_area =  self.QE_efficiency * self.Throughput * self.Atmosphere  *    (self.Collecting_area * 100 * 100)
        # TODO verify that indeed it should not depend on self.pixel_scale**2 !! We still see some dependancy, why that??

        self.nfibers = 1
        # TODO need to verify that the evolution the sky and signal makes sense with nfibers... Maybe this has been solved
        self.source_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength)   *          (self.Size_source   if self.IFS else   np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitwidth)   )
        self.slit_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength)   *          (np.maximum(self.Size_source,self.Slitwidth)   if self.IFS else   self.Slitwidth   )


        if self.spectro: # previously was multiplying by self.nfibers *
            # mat's solution provides a local optimum in dispersion that I don't get with my solution!!!
            self.factor_CU2el_tot = 1*self.effective_area * self.arcsec2str * np.minimum(self.Line_width,self.Bandwidth)   *  self.source_size_arcsec_after_slit  / self.pixels_total_source  
            self.factor_CU2el_sky_tot = 1*self.effective_area * self.arcsec2str *  np.minimum(self.Line_width,self.Bandwidth) *  self.slit_size_arcsec_after_slit / self.pixels_total_source  #it works only for a line emission and we take the total sky flux over the same pixels

            # working almost for dispersion! But actually we also need to have an optimal with pixel_scale!!!
            # self.factor_CU2el = self.effective_area * self.arcsec2str  * np.minimum(self.Slitwidth, self.Size_source)  * np.minimum( self.dispersion, np.minimum(self.Line_width, self.Bandwidth)/2.35) * self.pixel_scale
            # # il faut ici que je rajoute une line width
            # # quand la ligne width devient importante 
            # sky_spectral_coverage =  np.minimum( self.dispersion, np.minimum(self.Line_width, self.Bandwidth)/2.35)#np.minimum(self.dispersion, self.Bandwidth)
            # self.factor_CU2el_sky = self.effective_area * self.arcsec2str  * self.Slitwidth * sky_spectral_coverage * self.pixel_scale


            # working for dispersion but adding pixels_scale and slit_width!!!
            source_spatial_arcsec = np.minimum(self.Slitwidth, self.Size_source)
            source_spectral_angstrom = np.minimum(self.Line_width, self.Bandwidth)#/2.35
            spatial_per_pix = source_spatial_arcsec / self.pixel_scale
            spectral_per_pix = source_spectral_angstrom / self.dispersion
            # Combine spatiale + spectrale en quadrature :
            effective_pix_size = np.sqrt(spatial_per_pix**2 + spectral_per_pix**2)
            # Puis :
            self.factor_CU2el_average = self.effective_area * self.arcsec2str * source_spatial_arcsec * source_spectral_angstrom / effective_pix_size
            sky_spatial_arcsec = self.Slitwidth
            sky_spectral_angstrom =  self.Bandwidth  #/2.35            #np.minimum(self.Line_width, self.Bandwidth)
            spatial_per_pix_sky = sky_spatial_arcsec / self.pixel_scale
            spectral_per_pix_sky = sky_spectral_angstrom / self.dispersion
            # Combine spatiale + spectrale en quadrature :
            effective_pix_size_sky = np.sqrt(spatial_per_pix_sky**2 + spectral_per_pix_sky**2)
            # Puis :
            self.factor_CU2el_sky_average = self.effective_area * self.arcsec2str * sky_spatial_arcsec * sky_spectral_angstrom / effective_pix_size_sky
            
            difference = np.abs(self.factor_CU2el_tot - self.factor_CU2el_average)
            ratio = difference / np.abs(self.factor_CU2el_tot)
            # if (difference > 0.1) | (ratio > 0.1):
            #     print("Warning: difference or ratio between the two methods to compute the factor is too high: ", difference, ratio)
            #     print("factor_CU2el_tot", self.factor_CU2el_tot, "factor_CU2el_average", self.factor_CU2el_average)            
            if  self.test:
                self.factor_CU2el = self.factor_CU2el_tot
                self.factor_CU2el_sky = self.factor_CU2el_sky_tot
            else:
                self.factor_CU2el = self.factor_CU2el_average
                self.factor_CU2el_sky = self.factor_CU2el_sky_average

        else: 
            # TODO for imager that already have some throughput, integrate over the throughput curve.
            self.factor_CU2el =   self.pixel_scale**2 * self.Throughput_FWHM 
            self.factor_CU2el_sky = self.pixel_scale**2 * self.Throughput_FWHM  


        self.N_images = self.acquisition_time*3600/(self.exposure_time + self.readout_time)
        self.N_images_true = self.N_images * (1-self.cosmic_ray_loss)

        #TODO Here should be sure about the calculation. There is another way of taking the entire flux if it is a line 
        self.sky = self.Sky_CU*self.factor_CU2el_sky*self.exposure_time  # el/pix/frame
        self.Sky_noise = np.sqrt(self.sky * self.ENF) 
        self.Signal_LU = convert_ergs2LU(self.Signal,self.wavelength)
        # if 1==0: # if line is totally resolved (for cosmic web for instance)
        #     self.Signal_el =  self.Signal_LU*self.factor_CU2el*self.exposure_time * self.flux_fraction_slit_applied  / self.spectral_resolution_pixel # el/pix/frame#     Signal * (sky / Sky_)  #el/pix
        # else: # if line is unresolved for QSO for instance
        self.Signal_el =  self.Signal_LU * self.factor_CU2el * self.exposure_time * self.flux_fraction_slit_applied   # el/pix/frame#     Signal * (sky / Sky_)  #el/pix






        # TODO in counting mode, Photon_fraction_kept should also be used for CIC
        self.RN_final = self.RN  * self.RN_fraction_kept / self.EM_gain 
        self.Additional_background = self.extra_background/3600 * self.exposure_time# e/pix/exp
        self.Additional_background_noise = np.sqrt(self.Additional_background * self.ENF)
        
        # number of images taken during one field acquisition (~2h)
        # print(self.flux_fraction_slit_applied)

        self.signal_noise = np.sqrt(self.Signal_el * self.ENF)     #el / resol/ N frame

        if self.spectro:
            self.lambda_stack = 1 
        self.N_resol_element_A = self.lambda_stack 
        self.factor =   np.sqrt(self.number_pixels_used) * np.sqrt(self.N_resol_element_A) * np.sqrt(self.N_images_true)
        self.Signal_resolution = self.Signal_el * self.factor**2# el/N exposure/resol
        self.signal_noise_nframe = self.signal_noise * self.factor
        self.Total_noise_final = self.factor*np.sqrt(self.signal_noise**2 + self.Dark_current_noise**2  + self.Additional_background_noise**2 + self.Sky_noise**2 + self.CIC_noise**2 + self.RN_final**2   ) #e/  pix/frame
        self.SNR = self.Signal_resolution / self.Total_noise_final
        self.snrs_per_pixel = self.Signal_el * self.N_images_true /   (self.Total_noise_final / np.sqrt(self.number_pixels_used) / np.sqrt(self.N_resol_element_A)) 
        
        if type(self.Total_noise_final + self.Signal_resolution) == np.float64:
            n=0
        else:
            n =len(self.Total_noise_final + self.Signal_resolution) 
        if n>1:
            for name in ["signal_noise","Dark_current_noise", "Additional_background_noise","Sky_noise", "CIC_noise", "RN_final","Signal_resolution","Signal_el","sky","CIC_charge","Dark_current_f","RN","Additional_background"]:
                setattr(self, name, getattr(self,name)*np.ones(n))
        self.factor = self.factor*np.ones(n) if type(self.factor)== np.float64 else self.factor
        if type(self.number_pixels_used) == np.float64 : 
            self.noises_per_exp = np.sqrt(self.number_pixels_used) * np.array([self.signal_noise,  self.Dark_current_noise,  self.Sky_noise, self.RN_final, self.CIC_noise, self.Additional_background_noise, np.sqrt(self.Signal_el)]).T
        else:
            self.noises_per_exp = (np.sqrt(self.number_pixels_used) * np.array([self.signal_noise,  self.Dark_current_noise,  self.Sky_noise, self.RN_final, self.CIC_noise, self.Additional_background_noise, np.sqrt(self.Signal_el)])).T
        self.noises = np.array([self.signal_noise*self.factor,  self.Dark_current_noise*self.factor,  self.Sky_noise*self.factor, self.RN_final*self.factor, self.CIC_noise*self.factor, self.Additional_background_noise*self.factor, self.Signal_resolution]).T
        self.electrons_per_pix =  np.array([self.Signal_el,  self.Dark_current_f,  self.sky,  0*self.RN_final, self.CIC_charge, self.Additional_background]).T
        self.names = ["Signal","Dark current", "Sky", "Read noise","CIC", "Extra background"]
        # self.snrs = self.Signal_resolution /self.Total_noise_final

        if np.ndim(self.noises)==2:
            self.percents =  100* np.array(self.noises).T[:-1,:]**2/self.Total_noise_final**2
        else:
            self.percents =  100* np.array(self.noises).T[:-1]**2/self.Total_noise_final**2            
        
        self.el_per_pix = self.Signal_el + self.sky + self.CIC_charge +  self.Dark_current_f
        n_sigma = 5
        self.signal_nsig_e_resol_nframe = (n_sigma**2 * self.ENF + n_sigma * np.sqrt(4*self.Total_noise_final**2 - 4*self.signal_noise_nframe**2 + self.ENF**2*n_sigma**2))/2
        # self.signal_nsig_e_resol_nframe = 457*np.ones(self.len_xaxis)
        self.eresolnframe2lu = self.Signal_LU/self.Signal_resolution #TBV
        self.signal_nsig_LU = self.signal_nsig_e_resol_nframe * self.eresolnframe2lu #TBV
        self.signal_nsig_ergs = convert_LU2ergs(self.signal_nsig_LU, self.wavelength) 
        self.extended_source_5s = self.signal_nsig_ergs * (self.PSF_RMS_det*2.35)**2

        self.SB_lim_per_pix = self.signal_nsig_ergs
        self.SB_lim_per_res = self.signal_nsig_ergs / self.elem_size
        self.SB_lim_per_source = self.signal_nsig_ergs / self.source_size

        #TODO change this ratio of 1.30e57
        self.point_source_5s = self.extended_source_5s * 1.30e57
        self.time2reach_n_sigma_SNR = self.acquisition_time *  np.square(n_sigma / self.SNR)


        # print("Effective_area = ",np.unique(self.effective_area), "\n",
        # "Exposure time = ", np.unique(self.exposure_time[self.i]),"\n",
        # "Pixels_involved = %i * %i x %i"%( np.sqrt(source_spatial_pixels**2 + source_spectral_pixels**2)  ,   np.maximum(np.minimum(self.Size_source * fwhm_sigma_ratio, self.Slitlength) / self.pixel_scale,1)    ,    ( np.ceil(np.sqrt(self.Size_source**2+self.PSF_RMS_mask**2)*fwhm_sigma_ratio / self.Slitwidth) if self.IFS else 1)), np.unique(self.pixels_total_source),"\n",
        # "elem size pix = ", np.unique(self.elem_size),"\n",
        # "Max_exposures = ",np.unique(self.N_images[self.i]), np.unique(self.N_images_true[self.i]),"\n",
        # "Sky_CU = ",np.unique(self.Sky_CU),"\n",
        # "Signal_CU = ",np.unique(self.Signal_LU)," reduction (%i) \n"%(100*self.flux_fraction),
        # "Signal_e-_res_tot = ",np.unique(self.Signal_resolution[self.i]),"\n",
        # "self.Signal_el = ",np.unique(self.Signal_el[self.i]),"\n",
        # "self.Signal_resolution = ",np.unique(self.Signal_resolution[self.i]),"\n",
        # "self.contributions_resolution [signal, dark, sky, rn, cic, add] = ",self.electrons_per_pix[self.i]* self.factor[self.i]**2,"good - \n",
        # "noises_per_exp = ",self.noises_per_exp[self.i],"\n",
        # "noises = ",self.noises[self.i],"\n",
        # "verifications noises = ", (self.factor*np.array([self.signal_noise, self.Sky_noise, self.RN_final, self.Dark_current_noise, self.CIC_noise, self.Additional_background_noise])   )[:,self.i],"\n",
        # "total noise = ", (self.factor*np.sqrt(self.signal_noise**2 + self.Dark_current_noise**2  + self.Additional_background_noise**2 + self.Sky_noise**2 + self.CIC_noise**2 + self.RN_final**2   ))[self.i],"\n",
        # "ENF = ",self.ENF,"\n",
        # # "Source_e-_source_tot = ",,"\n",
        # "SNR = ",np.unique(self.SNR[self.i]))


        # self.atm_trans=np.nan
        # self.final_sky=np.nan
        # self.Throughput_curve=np.nan
        # print(self.acquisition_time, self.exposure_time[self.i] , self.readout_time)
        # print(self.N_images_true[self.i], self.N_images[self.i] , self.cosmic_ray_loss[self.i])
        # print("E2E troughput",int((100*self.QE_efficiency * self.Throughput * self.Atmosphere)[self.i]) , "\nFrame number=",self.N_images_true[self.i],"\nResolElem=",self.number_pixels_used, "\nSignal=",self.Signal_resolution[self.i])
        # print("Sigma=5")
        # print("Flux (e/rsol/Nframe), σ==5 :",self.signal_nsig_e_resol_nframe[self.i])
        # print("Flux LU, σ==5 : %0.1E"%(self.signal_nsig_LU[self.i]))
        # print("Flux	erg/cm2/s/''2/Å :",self.signal_nsig_ergs[self.i])
        # print("Flux	overPSF :",self.extended_source_5s[self.i])
        # print("Flux	point source :",self.point_source_5s[self.i])
        # print("factor=",self.factor[self.i])
        # print("N_images_true=",np.sqrt(self.N_images_true)[self.i] )
        # print("resolution_element=", self.number_pixels_used)
        # print("N_resol_element_A=",np.sqrt(self.N_resol_element_A))
        # print("lambda_stack=",self.lambda_stack)
        # print("dispersion=",self.dispersion)
        # print("cosmic_ray_loss=",np.sqrt(self.cosmic_ray_loss)[self.i])
        # print("N_images=",np.sqrt(self.N_images)[self.i])
        # from astropy.cosmology import Planck15 as cosmo
        # 4*np.pi* (cosmo.luminosity_distance(z=0.7).to("cm").value)**2 = 2.30e57


       

    def PlotNoise(self,title='',x='exposure_time', lw=8):
        """
        Generate a plot of the evolution of the noise budget with one parameter:
        exposure_time, Sky_CU, acquisition_time, Signal, EM_gain, RN, CIC_charge, Dark_current, readout_time, smearing, temperature, PSF_RMS_det, PSF_RMS_mask, QE, extra_background, cosmic_ray_loss_per_sec
        """
        fontsize = 8
        fig, axes= plt.subplots(4, 1, figsize=(12, 8), sharex=True) # fig, (ax1, ax2,ax3) = plt.subplots(3, 1, figsize=(12, 7), sharex=True) #figsize=(9, 5.5)
        ax1, ax2,ax3, ax4  = axes # TODO maybe add a sqrt here
        labels = ['%s: %0.3f (%0.1f%%)'%(name, self.number_pixels_used *self.electrons_per_pix[self.i,j],100*self.electrons_per_pix[self.i,j]/np.nansum(self.electrons_per_pix[self.i,:])) for j,name in enumerate(self.names)]



        # ax1 
        for i,(name,c) in enumerate(zip(self.names,self.colors)):
            # ax1.plot(getattr(self,x), self.noises[:,i]/self.factor,label='%s: %0.2f (%0.1f%%)'%(name,self.noises[self.i,i]/self.factor[self.i],self.percents[i,self.i]),lw=lw,alpha=0.8,c=c)
            ax1.plot(getattr(self,x), self.noises_per_exp[:,i],label='%s: %0.2f (%0.1f%%)'%(name,self.noises_per_exp[self.i,i],self.percents[i,self.i]),lw=lw,alpha=0.8,c=c)
        # ax1.plot(getattr(self,x), np.nansum(self.noises[:,:-1],axis=1)/self.factor,label='%s: %0.2f (%0.1f%%)'%("Total",np.nansum(self.noises[self.i,-1])/self.factor[self.i],np.nansum(self.percents[:,self.i])),lw=lw,alpha=0.4,c="k")
        ax1.plot(getattr(self,x), np.sqrt(np.nansum(np.multiply(self.noises_per_exp[:,:-1],self.noises_per_exp[:,:-1]),axis=1)) ,label='%s: %0.2f'%("Quadratic sum",np.sqrt(np.nansum(np.multiply(self.noises_per_exp[self.i,:-1],self.noises_per_exp[self.i,:-1])))   ) ,lw=lw,alpha=0.4,c="k") #np.nansum(self.percents[:,self.i])

        ax1.legend(loc='upper left', fontsize=fontsize)



        # ax2 
        ax2.grid(False)
        # self.stackplot1 = ax2.stackplot(getattr(self,x),  np.array(self.electrons_per_pix).T[:,:],alpha=0.7,colors=self.colors,labels=labels)
        self.stackplot1 = ax2.stackplot(getattr(self,x),  self.number_pixels_used * np.array(self.electrons_per_pix).T[:,:],alpha=0.7,colors=self.colors,labels=labels)
        ax2.legend(loc='upper left',title="Overall background: %0.3f (%0.1f%%)"%( self.number_pixels_used * np.nansum(self.electrons_per_pix[self.i,1:]),100*np.nansum(self.electrons_per_pix[self.i,1:])/np.nansum(self.electrons_per_pix[self.i,:])), fontsize=fontsize)
        ax2.set_xlim((getattr(self,x).min(),getattr(self,x).max()))

        # ax3
        ax3.grid(False)
        self.stackplot2 = ax3.stackplot(getattr(self,x), self.SNR * np.array(self.noises).T[:-1,:]**2/self.Total_noise_final**2,alpha=0.7,colors=self.colors)
        ax3.set_ylim((0,np.nanmax(self.SNR)))

        if self.SNR_res:
            ax1.set_ylabel('Noise (e-/Res/exp)')
            ax2.set_ylabel('e-/Res/frame')
            ax3.set_ylabel('SNR (Res, N frames)')        
        else:
            ax1.set_ylabel('Noise (e-/pix/exp)')
            ax2.set_ylabel('e-/pix/frame')
            ax3.set_ylabel('SNR (pix, N frames)')        

# self.extended_source_5s * 1.30e57
        ax4.plot(getattr(self,x), np.log10(self.SB_lim_per_pix),"-",lw=lw-1,label="SNR=5 limiting SB/power per pixel (%0.2f-%0.2f)"%(np.log10(1.30e57*self.SB_lim_per_pix[self.i]),np.nanmin(np.log10(1.30e57*self.SB_lim_per_pix))),c="k")
        ax4.plot(getattr(self,x), np.log10(self.SB_lim_per_res),"-",lw=lw-1,label="SNR=5 limiting SB/power per elem resolution (%0.2f-%0.2f)"%(np.log10(1.30e57*self.SB_lim_per_res[self.i]),np.nanmin(np.log10(1.30e57*self.SB_lim_per_res))),c="k",alpha=0.5)
        ax4.plot(getattr(self,x), np.log10(self.SB_lim_per_source),"-",lw=lw-1,label="SNR=5 limiting SB/power per source (%0.2f-%0.2f)"%(np.log10(1.30e57*self.SB_lim_per_source[self.i]),np.nanmin(np.log10(1.30e57*self.SB_lim_per_source))),c="k",alpha=0.3)
        # ax4.plot(getattr(self,x), np.log10(self.extended_source_5s/np.sqrt(2)),"-",lw=lw-1,label="Two elem resolution (%0.2f-%0.2f)"%(np.log10(self.point_source_5s[self.i]/np.sqrt(2)),np.nanmin(np.log10(self.point_source_5s/np.sqrt(2)))),c="grey")

        T2 =  lambda x:np.log10(10**x/1.30e57)
        self.pow_2018 = 42.95
        self.pow_best = 41.74
        ax4b = ax4.secondary_yaxis("right", functions=(lambda x:np.log10(10**x * 1.30e57),T2))
        if ("FIREBall" in self.instrument) & (1==0):
            ax4.plot([getattr(self,x).min(),getattr(self,x).min(),np.nan,getattr(self,x).max(),getattr(self,x).max()],[T2(self.pow_2018),T2(self.pow_best),np.nan,T2(self.pow_2018),T2(self.pow_best)],lw=lw,label="2018 flight (%0.1f) - most optimistic case (%0.1f)"%(self.pow_2018,self.pow_best),c="r",alpha=0.5)
        self.T2=T2
        self.ax4b = ax4b
        ax4.legend(loc="upper right", fontsize=fontsize,title=r"$\mathbf{Left}$: Ext. Source surface brightness, $\mathbf{Right}$: Point source power",title_fontsize=fontsize )
        ax4.set_ylabel(r"Log(erg/cm$^2$/s/asec$^2$)")
        ax4b.set_ylabel(r" Log(erg/s)")

        axes[-1].set_xlabel(x)
        ax1.tick_params(labelright=True,right=True)
        ax2.tick_params(labelright=True,right=True)
        ax3.tick_params(labelright=True,right=True)
        fig.tight_layout(h_pad=0.01)
        return fig 

    
    def compute_optimal_threshold(self,flux = 0.1,dark_cic_sky_noise=None,plot_=False,title='',i=0,axes=None,size= (int(1e3),int(1e3)),size_bin=25, threshold=-1000):
        """ 
        Create a ADU value histogram and defin the threshold so that it gives the optimal SNR based on RN, smearing, noise, flux, gain
        Function is pretty slow so output of this function has been saved and can then directly be used with interpolation (see function interpolate_optimal_threshold)
        """
        #self.Signal_el if np.isscalar(self.Signal_el) else 0.3
        EM_gain = self.EM_gain if np.isscalar(self.EM_gain) else self.EM_gain[i]#1000
        RN = self.RN if np.isscalar(self.RN) else self.RN[i]#80
        CIC_noise = self.CIC_noise if np.isscalar(self.CIC_noise) else self.CIC_noise[i]
        dark_noise = self.Dark_current_noise if np.isscalar(self.Dark_current_noise) else self.Dark_current_noise[i]
        try:
            Sky_noise = self.Sky_noise_pre_thresholding if np.isscalar(self.Sky_noise_pre_thresholding) else self.Sky_noise_pre_thresholding[i]
        except AttributeError:
            raise AttributeError('You must use counting_mode=True to use compute_optimal_threshold method.')


        dark = dark_noise**2
        CIC = CIC_noise**2
        sky = Sky_noise**2
        im = np.random.poisson(flux+dark+CIC+sky, size=size)
        values,bins = np.histogram(im,bins=[-0.5,0.5,1.5,2.5])
        ConversionGain=1#/4.5
        imaADU = np.random.gamma(im, EM_gain) *ConversionGain
        bins = np.arange(np.min(imaADU)-5*RN*ConversionGain,np.max(imaADU)+5*RN*ConversionGain,25)
        # bins = np.linspace(-500,10000,400)
        #imaADU = (np.random.gamma(im, EM_gain) + np.random.normal(0, RN, size=size))*ConversionGain
        imaADU_copy = imaADU.copy()
        imaADU_copy += np.random.normal(0, RN, size=size)*ConversionGain

        if plot_:
            if axes is None:
                fig, (ax1, ax2) = plt.subplots(2,1,sharex=True,figsize=(12, 7))#,figsize=(9,5))
            else:
                fig=0
                ax1, ax2 = axes
                ax1.clear()
                ax2.clear()
            val0,_,l0 = ax1.hist(imaADU_copy[im==0],bins=bins,alpha=0.5,log=True,histtype='step',lw=0.7,color='k',label='Before ampl & smearing')
            val1,_,l1 = ax1.hist(imaADU_copy[im==1],bins=bins,alpha=0.5,log=True,histtype='step',lw=0.7,color='k')
            val2,_,l2 = ax1.hist(imaADU_copy[im==2],bins=bins,alpha=0.5,log=True,histtype='step',lw=0.5,color='k')
            # _,_,_ = ax1.hist(imaADU_copy.flatten(),bins=bins,alpha=0.5,log=True,histtype='step',lw=0.5,color='k')
            # val3,_,l3 = ax1.hist(imaADU[im==3],bins=bins,alpha=0.5,log=True,histtype='step',lw=0.5,color='k')
            # val4,_,l4 = ax1.hist(imaADU[im==4],bins=bins,alpha=0.5,log=True,histtype='step',lw=0.5,color='k')
            # val5,_,l5 = ax1.hist(imaADU[im==5],bins=bins,alpha=0.5,log=True,histtype='step',lw=0.5,color='k')


        if self.smearing > 0:
            # print(SmearExpDecrement)
            smearing_kernels = variable_smearing_kernels(
                imaADU, self.smearing, SmearExpDecrement=5e4)
            offsets = np.arange(6)
            A = dia_matrix(
                (smearing_kernels.reshape((6, -1)), offsets),
                shape=(imaADU.size, imaADU.size))

            imaADU = A.dot(imaADU.ravel()).reshape(imaADU.shape)
        imaADU += np.random.normal(0, RN, size=size)*ConversionGain

        if 1==1:
            b = (bins[:-1]+bins[1:])/2
            rn_frac = np.array([np.sum(val0[b>bi]) for bi in b])/np.sum(val0) 
            rn_noise = (RN/(EM_gain * ConversionGain)) * rn_frac #/(EM_gain*ConversionGain)#/(EM_gain*ConversionGain)
            # rn_noise = RN * np.array([np.sum(val0[b>bi]) for bi in b])/np.sum(val0) #/(EM_gain*ConversionGain)#/(EM_gain*ConversionGain)
            signal12 = flux * np.array([np.sum(val1[b>bi])+np.sum(val2[b>bi]) for bi in b])/(np.sum(val1)+np.sum(val2))
            kept = np.array([np.sum(val1[b>bi]) for bi in b])/np.sum(val1)
            signal1 = flux * kept
            pc = np.ones(len(b))#                 # ([np.sum(val1[b>bi])for bi in b]/(np.array([np.sum(val1[b>bi])for bi in b])+np.array([np.sum(val0[b>bi]) for bi in b])))
            pc =  ([np.sum(val1[b>bi])for bi in b]/(np.array([np.sum(val1[b>bi])for bi in b])+np.array([np.sum(val0[b>bi]) for bi in b])))

            if dark_cic_sky_noise is None:
                noise = CIC_noise**2+dark_noise**2+Sky_noise**2
            else:
                noise = dark_cic_sky_noise
            # SNR1 = signal1/np.sqrt(signal1+noise+np.array(rn_noise)**2)#
            SNR1 = pc*signal1/np.sqrt((flux+dark+CIC+sky)*kept+  ((1-(flux+dark+CIC+sky))*rn_frac)**2   )#

            SNR12 = pc*signal12/ np.sqrt(signal12+noise+np.array(rn_noise)**2)
            SNR_analogic = flux/np.sqrt(2*flux+2*noise+(RN/(EM_gain * ConversionGain))**2)
            threshold_55 = 5.5*RN*ConversionGain
            id_55 =  np.argmin(abs(threshold_55 - b))
            # b = (bins[:-1]+bins[1:])/2
            # rn_frac = np.array([np.sum(val0[b>bi]) for bi in b])/np.sum(val0) 
            # rn_noise = (RN/(EM_gain * ConversionGain)) * rn_frac #/(EM_gain*ConversionGain)#/(EM_gain*ConversionGain)
            # # rn_noise = RN * np.array([np.sum(val0[b>bi]) for bi in b])/np.sum(val0) #/(EM_gain*ConversionGain)#/(EM_gain*ConversionGain)
            # signal12 = flux * np.array([np.sum(val1[b>bi])+np.sum(val2[b>bi]) for bi in b])/(np.sum(val1)+np.sum(val2))
            # kept = np.array([np.sum(val1[b>bi]) for bi in b])/np.sum(val1)
            # signal1 = flux * kept

            # pc = np.ones(len(b))# 
            #     # ([np.sum(val1[b>bi])for bi in b]/(np.array([np.sum(val1[b>bi])for bi in b])+np.array([np.sum(val0[b>bi]) for bi in b])))
            # pc =  ([np.sum(val1[b>bi])for bi in b]/(np.array([np.sum(val1[b>bi])for bi in b])+np.array([np.sum(val0[b>bi]) for bi in b])))

            # if dark_cic_sky_noise is None:
            #     noise = np.sqrt(CIC_noise**2+dark_noise**2+Sky_noise**2)
            #     # noise = np.sqrt(CIC_noise+dark_noise+Sky_noise)
            # else:
            #     noise = dark_cic_sky_noise
            # # noise = np.sqrt(noise)
            # # print('noises = ',noise)
            # print("signal + dark + cic + Sky + rn = 1      :", "%0.2f + %0.2f+ %0.2f+%0.2f+%0.2f=%0.2f"%(flux,dark,CIC,sky,np.sum(val0)/np.sum(val0+val1+val2),  flux+dark+CIC+sky+np.sum(val0)/np.sum(val0+val1+val2)    ))
            # SNR1 = signal1/np.sqrt(signal1+noise+np.array(rn_noise)**2)#
            # # SNR1 = signal1/np.sqrt((flux+dark+CIC+sky)*kept+  ((1-(flux+dark+CIC+sky))*rn_frac)   )#
            # SNR12 = pc*signal12/ np.sqrt(signal12+noise+np.array(rn_noise)**2)
            # SNR_analogic = flux/np.sqrt(2*flux+2*noise+(RN/(EM_gain * ConversionGain))**2)
            # print('SNR_analogic = ',SNR_analogic)

            lw=3
            if plot_:
                ax2.plot(b,rn_frac,lw=lw,ls=":",c="C0")#,label='RN(RN>T):  %0.2f%% ➛ %0.2f%%'%(100*rn_frac[id_55],100*rn_frac[id_t])
                ax2.plot(b,signal1/flux,lw=lw,ls=":",c="C1")#,label='Signal(Signal>T):  %0.1f%% ➛ %0.1f%%'%(100*signal1[id_55]/flux,100*signal1[id_t]/flux)
                # ax2.plot(b,np.array(rn_noise)**2,label='(RN(RN>T)/EM_gain)**2',lw=lw)
                ax2.plot(b,pc,lw=lw,ls=":",c="C2")#,label='Fraction(T) of true positive: %0.1f%% ➛ %0.1f%%'%(100*pc[id_55],100*pc[id_t])
                #ax2.plot(b,SNR1/pc,label='SNR without fraction')

                # ax2.plot(b,SNR1/np.nanmax(SNR1),lw=lw,c="C4") # ,label='SNR1: %0.2f%% ➛ %0.2f%%'%(SNR1[id_55],SNR1[id_t])#'%(100*np.sum(val0[id_t:])/np.sum(val0),100*np.sum(val1[id_t:])/np.sum(val1)),lw=lw)
                # ax2.plot(b,SNR12,':',label='SNR12, [N1+N2]/[N0] = %0.2f, frac(N1+N2)=%i%%'%((val1[np.nanargmax(SNR12)]+val2[np.nanargmax(SNR12)])/val0[np.nanargmax(SNR12)],100*np.sum(val1[np.nanargmax(SNR12):]+val2[np.nanargmax(SNR12):])/(np.sum(val1)+np.sum(val2))),lw=lw)
                # TODO check if we should have pc or not!!!
                ax2.plot(b,SNR1/SNR_analogic,lw=lw,ls=":",c="C3")#,label='SNR1 PC / SNR analogic: %0.2f ➛ %0.2f'%(SNR1[id_55]/SNR_analogic,SNR1[id_t]/SNR_analogic)

        if plot_:
            val0,_,l0 = ax1.hist(imaADU[im==0],bins=bins,alpha=0.5,label='0',log=True)
            val1,_,l1 = ax1.hist(imaADU[im==1],bins=bins,alpha=0.5,label='1',log=True)
            val2,_,l2 = ax1.hist(imaADU[im==2],bins=bins,alpha=0.5,label='2',log=True)
            # val3,_,l3 = ax1.hist(imaADU[im==3],bins=bins,alpha=0.5,label='3',log=True)
            # val4,_,l4 = ax1.hist(imaADU[im==4],bins=bins,alpha=0.5,label='4',log=True)
            # val5,_,l5 = ax1.hist(imaADU[im==5],bins=bins,alpha=0.5,label='5',log=True)
            ax1.hist(imaADU.flatten(),bins=bins,label='Total histogram',log=True,histtype='step',lw=1,color='k')
            # ax1.fill_between([bins[np.argmin(val0>val1)],bins[np.argmax(val0>val1)]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C0")
            # a0 = np.where(val0>val1)[0].max()
            # a1 = np.where((val1>val2)&(val1>val0))[0].max()
            # a2 = np.where((val2>val3)&(val2>val1))[0].max()
            # a3 = np.where((val3>val4)&(val3>val2))[0].max()
            # a4 = np.where((val4>val5)&(val4>val3))[0].max()
            # ax1.fill_between([ bins[0],bins[a0]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C0")
            # ax1.fill_between([bins[a0],bins[a1]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C1")
            # ax1.fill_between([bins[a1],bins[a2]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C2")
            # ax1.fill_between([bins[a2],bins[a3]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C3")
            # ax1.fill_between([bins[a3],bins[a4]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C4")
            # ax1.fill_between([bins[a4],bins[-1]],[val0.max(),val0.max()],[1.2*val0.max(),1.2*val0.max()],alpha=0.3,color="C5")
        else:
            val0,_ = np.histogram(imaADU[im==0],bins=bins)#,alpha=0.5,label='0',log=True)
            val1,_ = np.histogram(imaADU[im==1],bins=bins)#,alpha=0.5,label='1',log=True)
            val2,_ = np.histogram(imaADU[im==2],bins=bins)#,alpha=0.5,label='2',log=True)
            val3,_ = np.histogram(imaADU[im==3],bins=bins)
            val4,_ = np.histogram(imaADU[im==4],bins=bins)
            val5,_ = np.histogram(imaADU[im==5],bins=bins)



        b = (bins[:-1]+bins[1:])/2
        rn_frac = np.array([np.sum(val0[b>bi]) for bi in b])/np.sum(val0) 
        rn_noise = (RN/(EM_gain * ConversionGain)) * rn_frac #/(EM_gain*ConversionGain)#/(EM_gain*ConversionGain)
        # rn_noise = RN * np.array([np.sum(val0[b>bi]) for bi in b])/np.sum(val0) #/(EM_gain*ConversionGain)#/(EM_gain*ConversionGain)
        signal12 = flux * np.array([np.sum(val1[b>bi])+np.sum(val2[b>bi]) for bi in b])/(np.sum(val1)+np.sum(val2))
        kept = np.array([np.sum(val1[b>bi]) for bi in b])/np.sum(val1)
        signal1 = flux * kept

        pc = np.ones(len(b))# 
              # ([np.sum(val1[b>bi])for bi in b]/(np.array([np.sum(val1[b>bi])for bi in b])+np.array([np.sum(val0[b>bi]) for bi in b])))
        pc =  ([np.sum(val1[b>bi])for bi in b]/(np.array([np.sum(val1[b>bi])for bi in b])+np.array([np.sum(val0[b>bi]) for bi in b])))

        if dark_cic_sky_noise is None:
            noise = np.sqrt(CIC_noise**2+dark_noise**2+Sky_noise**2)
            # noise = np.sqrt(CIC_noise+dark_noise+Sky_noise)
        else:
            noise = dark_cic_sky_noise
        # noise = np.sqrt(noise)
        # print('noises = ',noise)
        print("signal + dark + cic + Sky + rn = 1      :", "%0.2f + %0.2f+ %0.2f+%0.2f+%0.2f=%0.2f"%(flux,dark,CIC,sky,np.sum(val0)/np.sum(val0+val1+val2),  flux+dark+CIC+sky+np.sum(val0)/np.sum(val0+val1+val2)    ))
        # SNR1 = signal1 / np.sqrt(signal1+noise+np.array(rn_noise)**2)#
        SNR1 = pc*signal1/np.sqrt((flux+dark+CIC+sky)*kept+  ((1-(flux+dark+CIC+sky))*rn_frac)   )  /2#
        SNR12 = pc*signal12/ np.sqrt(signal12+noise+np.array(rn_noise)**2)
        SNR_analogic = flux/np.sqrt(2*flux+2*noise+(RN/(EM_gain * ConversionGain))**2)
        # print('SNR_analogic = ',SNR_analogic)
        threshold_55 = 5.5 * RN * ConversionGain
        id_55 =  np.nanargmin(abs(threshold_55 - b))
        if threshold<-5:
            id_t = np.nanargmax(SNR1)
            threshold = b[id_t]
        else:
            threshold *= RN*ConversionGain
            id_t = np.nanargmin(abs(threshold - b))
        # print(threshold)
        fraction_signal = np.nansum(val1[id_t:])/np.nansum(val1)
        fraction_rn = np.nansum(val0[id_t:])/np.nansum(val0)
        lw=3
        if plot_:
            ax2.plot(b,rn_frac,label='RN(RN>T):  %0.2f%% ➛ %0.2f%%'%(100*rn_frac[id_55],100*rn_frac[id_t]),lw=lw,c="C0")
            ax2.plot(b,signal1/flux,label='Signal(Signal>T):  %0.1f%% ➛ %0.1f%%'%(100*signal1[id_55]/flux,100*signal1[id_t]/flux),lw=lw,c="C1")
            # ax2.plot(b,np.array(rn_noise)**2,label='(RN(RN>T)/EM_gain)**2',lw=lw)
            ax2.plot(b,pc,label='Fraction(T) of true positive: %0.1f%% ➛ %0.1f%%'%(100*pc[id_55],100*pc[id_t]),lw=lw,c="C2")
            #ax2.plot(b,SNR1/pc,label='SNR without fraction')
            # print(SNR1)
            # print(SNR1/SNR1.max())
            # ax2.plot(b,SNR1/np.nanmax(SNR1),label='SNR1: %0.2f%% ➛ %0.2f%%'%(SNR1[id_55],SNR1[id_t]),lw=lw,c="C4") #'%(100*np.sum(val0[id_t:])/np.sum(val0),100*np.sum(val1[id_t:])/np.sum(val1)),lw=lw)
            # ax2.plot(b,SNR12,':',label='SNR12, [N1+N2]/[N0] = %0.2f, frac(N1+N2)=%i%%'%((val1[np.nanargmax(SNR12)]+val2[np.nanargmax(SNR12)])/val0[np.nanargmax(SNR12)],100*np.sum(val1[np.nanargmax(SNR12):]+val2[np.nanargmax(SNR12):])/(np.sum(val1)+np.sum(val2))),lw=lw)



            ax2.plot(b,SNR1/SNR_analogic,label='SNR1 PC / SNR analogic: %0.2f ➛ %0.2f'%(SNR1[id_55]/SNR_analogic,SNR1[id_t]/SNR_analogic),lw=lw,c="C3")
            # ax2.plot(b,SNR12/SNR_analogic,':',label='SNR12 PC / SNR analogic',lw=lw)
            # ax2.set_yscale('log')
            ax2.set_ylim(ymin=1e-5)
            
            # ax2.plot(b,SNR1,label='[N1]/[N0] = %0.2f, frac(N1)=%i%%'%(val1[id_t]/val0[id_t],100*np.sum(val1[id_t:])/np.sum(val1)))
            # ax2.plot(b,SNR12,label='[N1+N2]/[N0] = %0.2f, frac(N1+N2)=%i%%'%((val1[np.nanargmax(SNR12)]+val2[np.nanargmax(SNR12)])/val0[np.nanargmax(SNR12)],100*np.sum(val1[np.nanargmax(SNR12):]+val2[np.nanargmax(SNR12):])/(np.sum(val1)+np.sum(val2))))

            ax2.legend(title = "T = 5.5σ ➛ %0.1fσ "%(threshold/(RN*ConversionGain)), fontsize=10)
            ax2.legend(title = "T = 5.5σ ➛ %0.1fσ "%(threshold/(RN*ConversionGain)), fontsize=13)
            ax2.set_xlabel('ADU',fontsize=13)
            ax1.set_ylabel('Occurence',fontsize=13)
            ax2.set_ylabel('SNR',fontsize=13)
            ax1.plot([threshold,threshold],[0,np.max(val0)],':',c='k',label=r"SNR optimal threshold")
            ax2.plot([threshold,threshold],[0,1],':',c='k')
            ax1.plot([threshold_55,threshold_55],[0,np.max(val0)],'-.',c='k',label=r"5.5 $\sigma_{RN}$ threshold")
            ax2.plot([threshold_55,threshold_55],[0,1],'-.',c='k')
            L = ax1.legend(fontsize=10,loc="upper right")
            L = ax1.legend(fontsize=13,loc="upper right")
            # L.get_texts()[1].set_text('0 e- : %i%%, fraction kept: %0.2f%%'%(100*values[0]/(size[0]*size[1]),100*np.sum(val0[id_t:])/np.sum(val0)))
            # L.get_texts()[2].set_text('1 e- : %i%%, fraction kept: %0.2f%%'%(100*values[1]/(size[0]*size[1]),100*np.sum(val1[id_t:])/np.sum(val1)))
            # L.get_texts()[3].set_text('2 e- : %i%%, fraction kept: %0.2f%%'%(100*values[2]/(size[0]*size[1]),100*np.sum(val2[id_t:])/np.sum(val2)))

            L.get_texts()[1].set_text('0 e$^-$ : %0.1f%% of pixels'%(100*values[0]/(size[0]*size[1])))
            L.get_texts()[2].set_text('1 e$^-$ : %0.1f%% of pixels'%(100*values[1]/(size[0]*size[1])))
            L.get_texts()[3].set_text('2 e$^-$ : %0.1f%% of pixels'%(100*values[2]/(size[0]*size[1])))

            ax1.tick_params(axis='both',labelsize=12)
            ax2.tick_params(axis='both',labelsize=12)


            ax1.set_title(title+'Gain = %i, RN = %i, flux = %0.2f, smearing=%0.1f, Threshold = %i = %0.2f$σ$'%(EM_gain,RN,flux,self.smearing, threshold,threshold/(RN*ConversionGain)))
            ax1.set_xlim(xmin=bins.min(),xmax=5000)#bins.max())
            # ax1.set_xlim(xmin=bins.min(),xmax=bins.max()/2)
            if axes is None:
                fig.tight_layout()
            fig.savefig("/tmp/histogram_analysis.svg",bbox_inches='tight')
            plt.show()
            # return fig
        # sys.exit()
        return threshold/(RN*ConversionGain), fraction_signal, fraction_rn, np.nanmax(SNR1/SNR_analogic)
 


    def interpolate_optimal_threshold(self,flux = 0.1,dark_cic_sky_noise=None,plot_=False,title='',i=0):
        """
        Return the threshold optimizing the SNR
        """
        #self.Signal_el if np.isscalar(self.Signal_el) else 0.3
        EM_gain = self.EM_gain #if np.isscalar(self.EM_gain) else self.EM_gain[i]
        RN= self.RN #if np.isscalar(self.RN) else self.RN[i]#80
        CIC_noise = self.CIC_noise #if np.isscalar(self.CIC_noise) else self.CIC_noise[i]
        dark_noise = self.Dark_current_noise #if np.isscalar(self.Dark_current_noise) else self.Dark_current_noise[i]
         
        try:
            Sky_noise = self.Sky_noise_pre_thresholding #if np.isscalar(self.Sky_noise_pre_thresholding) else self.Sky_noise_pre_thresholding[i]
        except AttributeError:
            raise AttributeError('You must use counting_mode=True to use interpolate_optimal_threshold method.')

        noise_value = CIC_noise**2+dark_noise**2+Sky_noise**2
        
        gains=np.linspace(800,2500,n)#self.len_xaxis)
        rons=np.linspace(30,120,n)#self.len_xaxis)
        fluxes=np.linspace(0.01,0.7,n)#self.len_xaxis)
        smearings=np.linspace(0,2,n)#self.len_xaxis)
        noise=np.linspace(0.002,0.05,n)#self.len_xaxis)
        if (n==6)|(n==10):
            coords = (gains, rons, fluxes, smearings)
            point = (EM_gain, RN, flux, self.smearing)            
        elif n==5:
            coords = (gains, rons, fluxes, smearings,noise)
            point = (EM_gain, RN, flux, self.smearing,noise_value)
        else:
            print(n,EM_gain, RN, flux, self.smearing,noise_value)
            
        if ~np.isscalar(noise_value) |  ~np.isscalar(self.smearing) | ~np.isscalar(EM_gain) | ~np.isscalar(RN):
            point = np.repeat(np.zeros((4,1)), self.len_xaxis, axis=1).T
            point[:,0] =  self.EM_gain
            point[:,1] = self.RN
            point[:,2] = flux
            point[:,3] = self.smearing
        fraction_rn =interpn(coords, table_fraction_rn, point,bounds_error=False,fill_value=None)
        fraction_signal =interpn(coords, table_fraction_flux, point,bounds_error=False,fill_value=None)
        threshold = interpn(coords, table_threshold, point,bounds_error=False,fill_value=None)
        snr_ratio = interpn(coords, table_snr, point,bounds_error=False,fill_value=None)

        if type(self.smearing)==float:
            if self.smearing == 0:
                a = Table.read("fraction_flux.csv")
                threshold = 5.5
                fraction_signal = np.interp(self.EM_gain/self.RN,a["G/RN"],a["fractionflux"])
            # fraction_rn = f(flux=0.1,EM_gain=self.EM_gain, RN=self.RN)
            # fraction_signal = f2(flux=0.1,EM_gain=self.EM_gain, RN=self.RN)
            # snr_ratio = f3(flux=0.1,EM_gain=self.EM_gain, RN=self.RN)

        return threshold, fraction_signal, fraction_rn, snr_ratio#np.nanmax(SNR1/SNR_analogic)
 



    def SimulateFIREBallemCCDImage(self,  Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Slit", size=[100, 100], OSregions=[0, 100], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=True,atmlambda=True,fraction_lya=0.05, Full_well=60, conversion_gain=1, Throughput_FWHM=200, Altitude=35,sky_lines=True,Redshift=0,IFS=False, source_image=False):
        # self.EM_gain=1500; Bias=0; self.RN=80; self.CIC_charge=1; p_sCIC=0; self.Dark_current=1/3600; self.smearing=1; SmearExpDecrement=50000; self.exposure_time=50; flux=1; self.Sky=4; source="Spectra m=17"; Rx=8; Ry=8;  size=[100, 100]; OSregions=[0, 120]; name="Auto"; spectra="Spectra m=17"; cube="-"; n_registers=604; save=False;self.readout_time=5;stack=100;self.QE=0.5
        from astropy.modeling.functional_models import Gaussian2D, Gaussian1D
        from scipy.sparse import dia_matrix
        from scipy.interpolate import interp1d
        self.fast = True
        self.Redshift=Redshift
        for key in list(self.instruments["Charact."]) + ["Signal_el","N_images_true","Dark_current_f","sky"]:
            if hasattr(self,key ):
                if (type(getattr(self,key)) != float) & (type(getattr(self,key)) != int) &  (type(getattr(self,key)) != np.float64) &  (type(getattr(self,key)) != bool):
                    setattr(self, key,getattr(self,key)[self.i])
                    # print(getattr(self,key))


        self.point_source_spectral_resolution = (10*self.wavelength)/self.Spectral_resolution
        self.diffuse_spectral_resolution = np.sqrt(self.point_source_spectral_resolution**2+(self.Slitwidth*self.dispersion/self.pixel_scale)**2)
        conv_gain=conversion_gain
        OS1, OS2 = OSregions
        ConversionGain = conv_gain
        Bias=0
        image = np.zeros((size[1], size[0]), dtype="float64")
        image_without_source = np.zeros((size[1], size[0]), dtype="float64")
        image_only_source = np.zeros((size[1], size[0]), dtype="float64")
        image_stack = np.zeros((size[1], size[0]), dtype="float64")
        image_stack_without_source = np.zeros((size[1], size[0]), dtype="float64")
        image_stack_only_source = np.zeros((size[1], size[0]), dtype="float64")

        source_im = 0 * image[:, OSregions[0] : OSregions[1]]
        self.sky_im = np.ones(image[:, OSregions[0] : OSregions[1]].shape) 
        source_im_wo_atm = 0 * image[:, OSregions[0] : OSregions[1]]
        lx, ly = source_im.shape
        y = np.linspace(0, lx - 1, lx)
        x = np.linspace(0, ly - 1, ly)
        x, y = np.meshgrid(x, y)

        stack = int(self.N_images_true)
        flux = (self.Signal_el /self.exposure_time)
        Rx = self.PSF_RMS_det/self.pixel_scale
        PSF_x = np.sqrt((np.nanmin([self.Size_source/self.pixel_scale,self.Slitlength/self.pixel_scale]))**2 + (Rx)**2)
        PSF_λ = np.sqrt((self.diffuse_spectral_resolution/self.dispersion)**2 + (self.Line_width/self.dispersion)**2)

        wave_min, wave_max = 10*self.wavelength - (size[0]/2) * self.dispersion , 10*self.wavelength + (size[0]/2) * self.dispersion
        nsize2, nsize = size
        wavelengths = np.linspace(wave_min,wave_max,nsize2)
        if os.path.exists("../data/Instruments/%s/Throughput.csv"%(self.instrument.upper().replace(" ","_"))):
            QE = Table.read("../data/Instruments/%s/Throughput.csv"%(self.instrument.upper().replace(" ","_")))
            QE = interp1d(QE[QE.colnames[0]]*10,QE[QE.colnames[1]])#
            self.Throughput_curve = QE(wavelengths)/np.nanmax(QE(wavelengths))  if QElambda else Gaussian1D.evaluate(wavelengths,  1,  self.wavelength*10, Throughput_FWHM )
        
        else:
            self.Throughput_curve = Gaussian1D.evaluate(wavelengths,  1,  self.wavelength*10, Throughput_FWHM )  if QElambda else 1  #self.QE
        if os.path.exists("../data/Instruments/%s/Atmosphere_transmission.csv"%(self.instrument.upper().replace(" ","_"))):
            trans = Table.read("../data/Instruments/%s/Atmosphere_transmission.csv"%(self.instrument.upper().replace(" ","_")))
            resolution_atm = self.diffuse_spectral_resolution/(10*(wavelengths[2]-wavelengths[1]))
            trans["trans_conv"] = gaussian_filter1d(trans[trans.colnames[0]], resolution_atm/2.35)
            self.atm_trans_before_convolution =  interp1d(list(trans[trans.colnames[0]]*10),list(trans[trans.colnames[1]]))(wavelengths)
            self.atm_trans = gaussian_filter1d(self.atm_trans_before_convolution, resolution_atm/2.35)
            self.atm_trans = self.atm_trans/np.nanmax(self.atm_trans)   if (atmlambda & (Altitude<100) ) else 1#self.Atmosphere

        elif Altitude<10: # only for ground instruments (based on altitude column)
            trans = Table.read("../data/Atm_transmission/pwv_atm_combined_ground.csv")
            self.atm_trans_before_convolution =  interp1d(list(trans["wave_microns"]*1000), list(trans["transmission"]))(wavelengths)
            resolution_atm = self.diffuse_spectral_resolution/(wavelengths[1]-wavelengths[0])
            self.atm_trans = gaussian_filter1d(self.atm_trans_before_convolution/np.nanmax(self.atm_trans_before_convolution),resolution_atm/2.35)       if atmlambda else 1 #self.Atmosphere
            # plt.figure();plt.plot(wavelengths,QE);plt.title("throughtput_fwhm: %f"%(Throughput_FWHM));plt.xlabel("Angstrom");plt.show()
        else:
            self.atm_trans_before_convolution = 1#self.Atmosphere
            self.atm_trans = 1#self.Atmosphere
            # self.Throughput_curve = Gaussian1D.evaluate(wavelengths,  self.QE,  self.wavelength*10, Throughput_FWHM )  if QElambda else self.QE
            # self.Throughput_curve = Gaussian1D.evaluate(wavelengths,  1,  self.wavelength*10, Throughput_FWHM )  if QElambda else 1  #self.QE
        # TODO should we really devide by atmosphere?
        atm_qe_normalized_shape =  np.ones(nsize2) * self.atm_trans * self.Throughput_curve #/ (self.QE*self.Atmosphere) 

        if (Altitude<100)  & (wavelengths.max()<10425) & (sky_lines):   #& (wavelengths.min()>3141)
            if os.path.exists("../data/Instruments/%s/Sky_emission_lines.csv"%(self.instrument.replace(" ","_"))):
                sky_lines = Table.read("../data/Instruments/%s/Sky_emission_lines.csv"%(self.instrument.replace(" ","_")))
            else:
                sky_lines = Table.read("../data/Sky_emission_lines/spectra_0.2A.csv")
                mask = (10*sky_lines[sky_lines.colnames[0]]>wavelengths.min()-10*self.dispersion) & (10*sky_lines[sky_lines.colnames[0]]<wavelengths.max()+10*self.dispersion)
                sky_lines = sky_lines[mask]

            sky = interp1d(10*sky_lines[sky_lines.colnames[0]],sky_lines[sky_lines.colnames[1]])(wavelengths)
            self.final_sky_before_convolution = (self.sky/self.exposure_time) * sky * (self.Sky/1e-16) #/np.mean(self.sky)   # 
            sky_model_interval = 10 * (sky_lines[sky_lines.colnames[0]][1]-sky_lines[sky_lines.colnames[0]][0])
            self.final_sky = gaussian_filter1d(self.final_sky_before_convolution, self.diffuse_spectral_resolution/2.35/sky_model_interval)

        else:
            self.final_sky_before_convolution = self.sky/self.exposure_time*np.ones(nsize2)
            self.final_sky = self.sky/self.exposure_time*np.ones(nsize2)

        length = self.Slitlength/2/self.pixel_scale
        a_ = special.erf((length - (np.linspace(0,nsize,nsize) - nsize/2)) / np.sqrt(2 * Rx ** 2))
        b_ = special.erf((length + (np.linspace(0,nsize,nsize) - nsize/2)) / np.sqrt(2 * Rx ** 2))
        slit_profile = (a_ + b_) / np.ptp(a_ + b_)  # Shape: (100,)
        slit_profile = (slit_profile - np.nanmin(slit_profile) )/ np.nanmax(slit_profile - np.nanmin(slit_profile) )

        if ("baseline" in source.lower()) | ("Spectra" in source) | ("Salvato" in source) | ("COSMOS" in source)| ("Blackbody" in source)| ("kpc spiral galaxy" in source)| ("Gaussian" in source)| ("cube" in source):
            # TODO should I replace PSF_x by PSF_x**2+Rx**2???? Issue with the normilization maybe... Need to compare to Bruno's code
            spatial_profile = Gaussian1D.evaluate(np.arange(size[1]),  1,  size[1]/2, PSF_x)
            if self.spectro:
                if ("baseline" in source.lower()) | (("UVSpectra=" in source) & (self.wavelength>300)): 
                    if ("CIV" in source)| ("CIII" in source)| ("OVI" in source)| ("Lyα" in source):
                        w_nm = (float(re.search(r'\d+', source).group())/10 ) *(1+self.Redshift)
                    else:
                        w_nm=self.wavelength
    
                    
                    spectra = flux* Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2 + (w_nm-self.wavelength)*10/self.dispersion, PSF_λ)#
                    # print(w_nm,self.wavelength,self.dispersion, size[0]/2 + (w_nm-self.wavelength)*10*self.dispersion,np.min(spectra),np.max(spectra))
                    if spectra.sum()>0:
                        spectra /= spectra.sum()
                    # print(spectra)
                    # else:
                    #     spectra = np.ones(size[0])
                    # / Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2 + (w_nm-self.wavelength)*10*self.dispersion, self.PSF_lambda_pix**2/(PSF_λ**2 + self.PSF_lambda_pix**2)).sum()
                    spectra *= atm_qe_normalized_shape   
                    source_im =  np.outer(spectra,spatial_profile*slit_profile ).T /Gaussian1D.evaluate(np.arange(size[1]),  1,  50, Rx**2/(PSF_x**2+Rx**2)).sum()
                elif "blackbody" in source.lower():
                    temperature = int(re.search(r'\d+', source).group()) * u.K
                    blackbody_spectra = BlackBody(temperature=temperature/(1+self.Redshift))(wavelengths * u.nm)
                    # Convertir le spectre de corps noir en unités désirées: erg / (cm^2 * s * Å * arcsec^2)
                    flux_in_erg = blackbody_spectra.to(
                        u.erg / (u.cm**2 * u.s * u.AA * u.arcsec**2),
                        equivalencies=u.spectral_density(wavelengths * u.nm))
                    # Ajuster le spectre du corps noir pour correspondre au flux moyen donné sur le détecteur
                    spectra = atm_qe_normalized_shape  * (flux_in_erg.value/np.mean(flux_in_erg.value)) * flux    #* self.QE*self.Atmosphere
                    source_im =  np.outer(spectra,spatial_profile*slit_profile ).T /Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, Rx**2/(PSF_x**2+Rx**2)).sum()
                elif ("cube" in source):
                    n_wave=size[0]
                    if os.path.isfile("../data/Emission_cube/"+ source.split("-")[0] + ".fits"):
                        # print(source)
                        phys = True if source.split("-")[1]=="resampled_phys" else False
                        # if source.split("-")[0] == "lya_cube_merged_with_artificial_source_CU_1pc":
                        name_ = "_resampled_phys.fits" if phys else "_resampled.fits"
                        new_cube = convert_fits_cube("../data/Emission_cube/"+ source.split("-")[0] + ".fits","../data/Emission_cube/"+ source.split("-")[0] + name_, output_shape=(nsize2,100,100),wave_range_nm=(wave_min/10,wave_max/10), spatial_extent_arcsec=100*self.pixel_scale,redshift=Redshift,phys=phys)                         
                        cube_detector =  fits.open(new_cube)[0].data
                        # print(1, np.max(cube_detector), np.min(cube_detector))
                        cube_detector -= np.nanmedian(cube_detector)
                        cube_detector[cube_detector<np.nanmedian(cube_detector)]= np.nanmedian(cube_detector)
                        # print(2, np.max(cube_detector), np.min(cube_detector))
                        # else:
                        #     cube_detector = fits.open("../data/Emission_cube/"+ source + ".fits")[0].data
                        # if np.nanmin(cube_detector) < 0:
                        # cube_detector[cube_detector==np.nanmin(cube_detector)] =np.nanmedian(cube_detector)
                        # print(np.nanmin(cube_detector))
                        # if remap:
                        # else:
                        # cube_detector = self.Signal_el * cube_detector/ np.nanmax(cube_detector)
                        fitswrite(cube_detector,"/tmp/gal_simu.fits")
                        cube_detector = np.transpose(cube_detector, (1, 2, 0))
                        # print(3, np.nanmin(cube_detector), np.nanmin(cube_detector))
                    else:
                        cube_detector = create_fits_cube(cube_size=(size[1], size[1], n_wave), pixel_scale=self.pixel_scale, wave_range=(wave_min,wave_max),continuum_flux=np.ones(n_wave)*self.extra_background * int(self.exposure_time)/3600 ,continuum_fwhm=None, line_flux=self.Signal_el, line_center=self.wavelength, line_fwhm=self.Line_width, line_spatial_fwhm=self.Size_source,galaxy_shape=True, kpc_size=int(source.split(" ")[1]), redshift=self.Redshift,filament=True)
                        fitswrite(cube_detector,"/tmp/gal_me.fits")
                    # print(cube_detector.min(),cube_detector.max(),np.argmax(cube_detector))
                    if IFS:

                        if source_image=="Source" :
                            return self.Signal * cube_detector/ np.percentile(cube_detector,99.999) , self.Signal_el * cube_detector/ np.percentile(cube_detector,99.999)
                        cube_detector = self.Signal_el * cube_detector/ np.percentile(cube_detector,99.999)              #np.nanmax(cube_detector)
                        cube_spatial = np.array([gaussian_filter(cube_detector[:, :, i], sigma=Rx) for i in range(cube_detector.shape[2])])
                        cube_spatial *= np.nanmax(cube_detector) /  np.nanmax(cube_spatial)
                        #TODO add 
                        fitswrite(cube_spatial, "/tmp/cube_spatial.fits")
                        cube_spatial = np.transpose(cube_spatial, (1, 2, 0))
                        # cube_spectral = np.array([gaussian_filter(cube_spatial[i, j, :], sigma=PSF_λ) for i in range(cube_spatial.shape[0]) for j in range(cube_spatial.shape[1])])
                        cube_detector = cube_spatial#.reshape(cube_spatial.shape)
            

                        cube_detector += self.sky + self.Dark_current_f  
                        Nx, Ny, Nλ = cube_detector.shape
                        reduction = int(self.Slitwidth/self.pixel_scale)
                        reduction = self.Slitwidth/self.pixel_scale
                        # new_Nx = Nx // reduction  # Nouvelle taille en X
                        new_Nx = int(Nx * self.pixel_scale/ self.Slitwidth)
                        # cube_reduced = np.array([np.nanmean(cube_detector[i * reduction:(i + 1) * reduction, :, :], axis=0) for i in range(new_Nx)])
                        cube_reduced = np.array([np.nanmean(cube_detector[int(i * reduction):int((i + 1) * reduction), :, :], axis=0) for i in range(new_Nx)])
                        cube_detector = np.array(cube_reduced)  # Remettre en array numpy
                        if source_image=="Convolved source" :
                            return cube_detector, cube_detector

                        if source_image=="SNR" : #need to be sure we are good for the SNR
                            SNR_all = (cube_reduced - np.nanmin(cube_reduced)) / np.ptp(cube_reduced - np.nanmin(cube_reduced)) * self.snrs_per_pixel[self.i] #/ np.sqrt(source_im + source_background + self.CIC_charge + self.RN**2)
                            SNR_single = SNR_all / np.sqrt(np.ones(50) * self.N_images_true)[self.i]
                            return SNR_single, SNR_all

                        cube_detector_stack = np.ones(cube_detector.shape)
                        if (self.EM_gain>1) : # TODO does not work!!!
                            if self.counting_mode : 
                                # print(2)
                                list_im =[(np.random.gamma(np.random.poisson(cube_detector)  + np.array(np.random.rand(*cube_detector.shape)<self.CIC_charge,dtype=int) , self.EM_gain)) + np.random.normal(0, self.RN, cube_detector.shape) for i in range(int(stack))]
                                # print(list_im)
                                cube_detector_stack = np.mean([np.array(l>5.5*self.RN) for l in list_im],axis=0)   
                                # print(cube_detector_stack)
                            else:
                                CIC =  np.mean([  np.array(np.random.rand(*cube_detector.shape)<self.CIC_charge,dtype=int)   for i in range(int(stack)) ],axis=0)
                                if int(stack)>0 :
                                    cube_detector_stack[:,:,:] =  np.random.gamma( np.random.poisson(cube_detector * stack) / int(stack) + CIC  , self.EM_gain) + np.random.normal(0, self.RN/np.sqrt(int(stack)), cube_detector.shape)                        
                                else :
                                    cube_detector_stack[:,:,:] =  0
                            cube_detector[:,:,:] += np.random.gamma( np.random.poisson(cube_detector) + np.array(np.random.rand(*cube_detector.shape)<self.CIC_charge,dtype=int) , self.EM_gain)
                        else:
                            # print(np.max(cube_detector * stack),np.nanmax(cube_detector * stack))
                            # cube_detector[:,:,:]+= np.random.poisson(cube_detector)
                            N_poisson = np.random.poisson(cube_detector * stack)           # Somme des Poisson sur la stack
                            N_CIC = np.random.binomial(stack, self.CIC_charge, cube_detector.shape) # Somme des CIC sur la stack
                            cube_detector_stack = np.random.gamma(N_poisson + N_CIC, self.EM_gain / stack)  + np.random.normal(0, self.RN/np.sqrt(int(stack)), cube_detector.shape)      if stack>0 else cube_detector*0                   
                        # fitswrite(np.transpose(cube_reduced, (2, 1, 0)), "/tmp/cube_reduced.fits")
                        cube_detector = np.random.poisson(cube_detector) + np.random.normal(0, self.RN, cube_detector.shape)
                        if self.counting_mode : 
                            counting_image = np.zeros(cube_detector.shape)
                            counting_image[cube_detector>5.5*self.RN]=1
                            cube_detector=counting_image
                        return cube_detector, cube_detector_stack
                        # return cube_detector, cube_detector

                        # return np.transpose(cube_detector, (1, 0, 2)), np.transpose(cube_detector_stack, (1, 0, ))  , 
                    else:
                        cube_detector = self.Signal_el * cube_detector/ np.percentile(cube_detector,99.999)              #np.nanmax(cube_detector)

                        pre_im = np.mean(cube_detector[:,np.max([0,int(50-self.Slitwidth/self.pixel_scale/2+self.Δx)]):int(50+self.Slitwidth/self.pixel_scale/2 +self.Δx+1),:], axis=1)
                        source_im =  atm_qe_normalized_shape.T*(pre_im.T * slit_profile).T/int(self.exposure_time)
                        if "Source" in source_image:
                            return source_im  * int(self.exposure_time) , source_im  * int(self.exposure_time)
                        # TODO add slit size!!!
                        image_test = gaussian_filter(source_im, sigma=(Rx,self.diffuse_spectral_resolution/self.dispersion))
                        # image_test *= np.nanmax(source_im) / np.nanmax(image_test)
                        source_im = image_test
                        if source_image=="Convolved source" :
                            return source_im  * int(self.exposure_time), source_im  * int(self.exposure_time)


                        # pre_im = np.mean(cube_detector[:,np.max([0,int(50-self.Slitwidth/self.pixel_scale/2+self.Δx)]):int(50+self.Slitwidth/self.pixel_scale/2 +self.Δx+1),:], axis=1)
                        # pre_im = atm_qe_normalized_shape.T * pre_im.T / int(self.exposure_time)
                        # if source_image=="Source" :
                        #     return pre_im, pre_im
                        # slit_profile = (slit_profile - np.nanmin(slit_profile) )/ np.nanmax(slit_profile - np.nanmin(slit_profile) )
                        # # TODO add slit size!!!
                        # image_test = gaussian_filter(pre_im, sigma=(Rx,self.diffuse_spectral_resolution/self.dispersion))
                        # image_test *= np.nanmax(pre_im) / np.nanmax(image_test)
                        # source_im = (image_test.T * slit_profile).T
                        # if source_image=="Convolved source" :
                        #     return source_im, source_im


                else:
                    if "_" not in source:
                        flux_name,wave_name ="FLUX", "WAVELENGTH"
                        fname = "h_%sfos_spc.fits"%(source.split(" ")[2])
                        # print(fname)
                        try:
                            a = Table.read("../data/Spectra/"+fname)
                        except FileNotFoundError: 
                            a = Table.read("/Users/Vincent/Github/notebooks/Spectra/" + fname)
                        a["photons"] = a[flux_name]/9.93E-12   
                        # TODO why there is no QE in the formula??
                        # a["e_pix_sec"]  = a["photons"] * self.Throughput * self.Atmosphere  * self.Collecting_area*100*100 *self.dispersion
                        a["e_pix_sec"]  = flux * a["photons"] / np.nanmax(a["photons"])
                    elif "COSMOS" in source:
                        a = Table.read("../data/Spectra/GAL_COSMOS_SED/%s.txt"%(source.split(" ")[2]),format="ascii")
                        wave_name,flux_name ="col1", "col2"
                        a[wave_name] = a[wave_name]*(1+self.Redshift)
                        mask = (a[wave_name]>wave_min - 100) & (a[wave_name]<wave_max+100)
                        a = a[mask]
                        a["e_pix_sec"] = flux * a[flux_name] / np.nanmax(a[flux_name])
                    elif "Salvato" in source:
                        a = Table.read("../data/Spectra/QSO_SALVATO2015/%s.txt"%(source.split(" ")[2]),format="ascii")
                        wave_name,flux_name ="col1", "col2"
                        a[wave_name] = a[wave_name]*(1+self.Redshift)
                        mask = (a[wave_name]>wave_min - 100) & (a[wave_name]<wave_max+100)
                        a = a[mask]
                        a["e_pix_sec"] = flux * a[flux_name] / np.nanmax(a[flux_name])
                    min_interval = np.nanmin(wavelengths[1:] - wavelengths[:-1])
                    mask = (a[wave_name]>wave_min) & (a[wave_name]<wave_max)
                    slits = None 
                    source_background=np.zeros((nsize,nsize2))
                    f = interp1d(a[wave_name],a["e_pix_sec"])
                    
                    spectra = gaussian_filter1d(f(wavelengths),  self.diffuse_spectral_resolution/min_interval/2.35) * atm_qe_normalized_shape      #* self.QE*self.Atmosphere
                    spatial_profile =  Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, PSF_x) #/Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, PSF_x).sum()
                    subim = np.zeros((nsize2,nsize))
                    source_im =  np.outer(spectra,spatial_profile*slit_profile ).T /Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, Rx**2/(PSF_x**2+Rx**2)).sum()
    

                if np.isfinite(length) & (np.ptp(a_ + b_)>0):
                    if self.Slitlength/self.pixel_scale<nsize:
                        self.sky_im =   np.outer(self.final_sky * self.Throughput_curve /self.QE, slit_profile ).T
                    else:
                        self.sky_im =   np.outer(self.final_sky * self.Throughput_curve /self.QE,  np.ones(nsize) / nsize ).T
                else:
                    self.sky_im =   np.outer(self.final_sky * self.Throughput_curve /self.QE, np.ones(size[1])   ).T                
            
            
                if IFS : 
                    n_wave=size[0]
                    cube_detector = create_fits_cube(cube_size=(size[1], size[1], n_wave), pixel_scale=self.pixel_scale, wave_range=(wave_min,wave_max),continuum_flux=spectra* int(self.exposure_time),continuum_fwhm=self.Size_source, line_flux=self.Signal_el*0, line_center=self.wavelength, line_fwhm=self.Line_width, line_spatial_fwhm=self.Size_source)
                    cube_detector += self.sky + self.Dark_current_f  + self.extra_background * int(self.exposure_time)/3600 
                    if "Source" in source_image:
                        return self.Signal * cube_detector/ np.percentile(cube_detector,99.999) ,cube_detector
                    cube_spatial = np.array([gaussian_filter(cube_detector[:, :, i], sigma=Rx) for i in range(cube_detector.shape[2])])
                    fitswrite(cube_spatial, "/tmp/cube_spatial.fits")
                    cube_spatial = np.transpose(cube_spatial, (1, 2, 0))
                    # cube_spectral = np.array([gaussian_filter(cube_spatial[i, j, :], sigma=PSF_λ) for i in range(cube_spatial.shape[0]) for j in range(cube_spatial.shape[1])])
                    cube_spatial *= np.nanmax(cube_detector) /  np.nanmax(cube_spatial)
                    cube_detector = cube_spatial#.reshape(cube_spatial.shape)
                    # if source_image:
                        

                    Nx, Ny, Nλ = cube_detector.shape
                    # n3 = int(np.sqrt(60*60*self.instruments_dict[self.instrument]["FOV_size"])/self.Slitwidth)
                    reduction = self.Slitwidth/self.pixel_scale
                    new_Nx = Nx // reduction  # Nouvelle taille en X
                    new_Nx = int(Nx * self.pixel_scale/ self.Slitwidth)
                    cube_reduced = np.array([np.nanmean(cube_detector[int(i * reduction):int((i + 1) * reduction), :, :], axis=0) for i in range(new_Nx)])
                    cube_detector = np.array(cube_reduced)  # Remettre en array numpy
                    cube_detector_stack = cube_detector.copy()
                    if source_image=="Convolved source" :
                        return cube_detector, cube_detector
            
                    if source_image=="SNR" : # TODO need to be sure we are good for the SNR which sqrt(50)
                        SNR_all = (cube_detector - np.nanmin(cube_detector)) / np.ptp(cube_detector - np.nanmin(cube_detector)) * self.snrs_per_pixel[self.i] #/ np.sqrt(source_im + source_background + self.CIC_charge + self.RN**2)
                        SNR_single = SNR_all / np.sqrt(np.ones(50) * self.N_images_true)[self.i]
                        return SNR_single, SNR_all
                    if (self.EM_gain>1) :
                        if self.counting_mode : 
                            list_im =[(np.random.gamma(np.random.poisson(cube_detector)  + np.array(np.random.rand(*cube_detector.shape)<self.CIC_charge,dtype=int) , self.EM_gain)) + np.random.normal(0, self.RN, cube_detector.shape) for i in range(int(stack))]
                            cube_detector_stack = np.mean([np.array(l>5.5*self.RN) for l in list_im],axis=0)   
                        else:
                            if self.fast :
                                N_poisson = np.random.poisson(cube_detector * stack)           # Somme des Poisson sur la stack
                                N_CIC = np.random.binomial(stack, self.CIC_charge, cube_detector.shape) # Somme des CIC sur la stack
                                cube_detector_stack = np.random.gamma(N_poisson + N_CIC, self.EM_gain / stack)  + np.random.normal(0, self.RN/np.sqrt(int(stack)), cube_detector.shape)     if stack>0 else cube_detector*0                    
                            else:
                                cube_detector_stack[:,:,:] = np.mean([(np.random.gamma(np.random.poisson(cube_detector)  + np.array(np.random.rand(*cube_detector.shape)<self.CIC_charge,dtype=int) , self.EM_gain)) + np.random.normal(0, self.RN, cube_detector.shape) for i in range(int(stack))],axis=0)                        
                        cube_detector[:,:,:] += np.random.gamma( np.random.poisson(cube_detector) + np.array(np.random.rand(*cube_detector.shape)<self.CIC_charge,dtype=int) , self.EM_gain)
                    else: # TODO weird that cube_detector_stack does not depend on exposure time. it only depends on stack
                        cube_detector_stack[:,:,:] = np.random.poisson(cube_detector*int(stack)) / int(stack)  + np.random.normal(0, self.RN/np.sqrt(int(stack)), cube_detector.shape) if int(stack)>0 else 0
                        cube_detector[:,:,:]+= np.random.poisson(cube_detector)
                        # cube_detector_stack[:,:,:] = np.mean(np.random.poisson(np.repeat(cube_detector[:, :,:, np.newaxis], int(stack), axis=3)),axis=3)
                        # cube_detector_stack += np.random.normal(0, self.RN/np.sqrt(int(stack)), cube_detector.shape)
                    # add dark
                    cube_detector += np.random.normal(0, self.RN, cube_detector.shape)
                    if self.counting_mode : 
                        counting_image = np.zeros(cube_detector.shape)
                        counting_image[cube_detector>5.5*self.RN]=1
                        cube_detector=counting_image
                    fitswrite(np.transpose(cube_reduced, (2, 1, 0)), "/tmp/cube_reduced.fits")
                    return cube_detector, cube_detector_stack           
            
            
            
            else:
                if ("Gaussian" in source):
                    self.sky_im *= self.sky/self.exposure_time
                    # TODO integrer le flux sur la fenetre (la correlation entre la bandpass et la raie d'emission )
                    gaussian = Gaussian2D(amplitude=flux, x_mean=nsize2/2, y_mean=nsize/2, x_stddev=PSF_x, y_stddev=PSF_x, theta=0)
                    X, Y = np.meshgrid(np.arange(nsize2), np.arange(nsize))
                    source_im = gaussian(X, Y)

                elif "kpc spiral galaxy" in source:
                    source_im = generate_spiral_galaxy(amplitude=flux, redshift=self.Redshift, pixel_scale=self.pixel_scale, PSF_RMS=self.PSF_RMS_det, size_kpc=int(source.split(" ")[0]))
                    CGM = generate_gaussian_galaxy(amplitude=flux*0.0001, redshift=self.Redshift, platescale=self.pixel_scale, PSF_RMS=self.PSF_RMS_det, size_kpc=300)
                    if "CGM" in source:
                        source_im += CGM
                    if "filament" in source:
                        source_im =  np.tile( CGM[:,250][:, np.newaxis], (1, 500))

                # TODO add star + planet
                # just IGM filament / cosmic web
                # hide the size of the PSF is size is fixed
                # reflechir pour les flux par 
                
        # if source_image :

        source_im_only_source =  source_im  * int(self.exposure_time)
        if ("Source" in source_image) | (source_image=="Convolved source") :
            # print(source_im_only_source.shape)
            return source_im_only_source, source_im_only_source
        source_background = self.sky_im  * int(self.exposure_time) + self.Dark_current_f  + self.extra_background * int(self.exposure_time)/3600 
        source_im =  source_background +  source_im_only_source

        if source_image=="SNR" :
            SNR_all = (source_im - np.nanmin(source_im)) / np.ptp(source_im - np.nanmin(source_im)) * self.snrs_per_pixel[self.i] #/ np.sqrt(source_im + source_background + self.CIC_charge + self.RN**2)
            SNR_single = SNR_all / np.sqrt(np.ones(50) * self.N_images_true)[self.i]
            return SNR_single, SNR_all

        # source_im_wo_atm = self.Dark_current_f + self.extra_background * int(self.exposure_time)/3600 +  source_im_wo_atm * int(self.exposure_time)
        y_pix=1000
        if (self.readout_time/self.exposure_time > 0.2) & (self.fast is False):
            cube = np.array([(self.readout_time/self.exposure_time/y_pix)*np.vstack((np.zeros((i,len(source_im))),source_im[::-1,:][:-i,:]))[::-1,:] for i in np.arange(1,len(source_im))],dtype=float)
            source_im = source_im+np.sum(cube,axis=0)
        # if self.cosmic_ray_loss_per_sec is None:
        #     self.cosmic_ray_loss_per_sec = np.minimum(0.005*(self.exposure_time+self.readout_time/2),1)#+self.readout_time/2
        stack = int(self.N_images_true)
        cube_stack = -np.ones((stack,size[1], size[0]), dtype="int32")


        n_smearing=6
        if (self.EM_gain>1) & (self.CIC_charge>0):
            image[:, OSregions[0] : OSregions[1]] += np.random.gamma( np.random.poisson(source_im) + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)
            # if self.fast is False:
            #     image_without_source[:, OSregions[0] : OSregions[1]] +=  np.random.gamma( np.random.poisson(source_background) + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)
            #     image_only_source[:, OSregions[0] : OSregions[1]] +=  np.random.gamma( np.random.poisson(source_im_only_source) , self.EM_gain)
        else:
            image[:, OSregions[0] : OSregions[1]] += np.random.poisson(source_im)
            # if self.fast is False:
            #     image_without_source[:, OSregions[0] : OSregions[1]] +=  np.random.poisson(source_background)
            #     image_only_source[:, OSregions[0] : OSregions[1]] +=  np.random.poisson(source_im_only_source)
        if self.fast :
            image_without_source, image_only_source, image_stack_only_source, image_stack_without_source, source_im_wo_atm, imaADU_stack_only_source, imaADU_without_source, imaADU_stack_without_source, imaADU_source = 0, 0, 0, 0, 0, 0, 0, 0,0


        # take into acount CR losses #18%
        # image_stack[:, OSregions[0] : OSregions[1]] = np.nanmean([np.where(np.random.rand(size[1], OSregions[1]-OSregions[0]) < self.cosmic_ray_loss_per_sec/n_smearing,np.nan,1) * (np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)) for i in range(int(stack))],axis=0)
        if self.EM_gain>1:
            if self.counting_mode : 
                # if self.fast is True:
                list_im =[(np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)) for i in range(int(stack))]
                image_stack[:, OSregions[0] : OSregions[1]] = np.mean([np.array(l>5.5*self.RN) for l in list_im],axis=0)   
                # else:
                #     # Génération vectorisée de toutes les images en une seule étape
                #     shape = (int(stack), *source_im.shape)  # Dimensions (stack, H, W)
                #     # 1. Génération Poisson + CIC pour toute la stack
                #     poisson_all = np.random.poisson(source_im, shape)
                #     cic_all = np.random.binomial(1, self.CIC_charge, shape)
                #     # 2. Application Gamma et RN en une seule opération
                #     gamma_all = np.random.gamma(poisson_all + cic_all, self.EM_gain)
                #     noisy_imgs = gamma_all + np.random.normal(0, self.RN, shape)
                #     # 3. Seuillage et moyenne vectorisés
                #     image_stack[:, OSregions[0] : OSregions[1]] = np.mean(noisy_imgs > 5.5*self.RN, axis=0)
            else:
                # image_stack[:, OSregions[0] : OSregions[1]] = np.mean([(np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)) for i in range(int(stack))],axis=0)
                N_poisson = np.random.poisson(source_im * stack)           # Somme des Poisson sur la stack
                N_CIC = np.random.binomial(stack, self.CIC_charge, (size[1], OSregions[1]-OSregions[0])) # Somme des CIC sur la stack
                image_stack[:, OSregions[0] : OSregions[1]]  = np.random.gamma(N_poisson + N_CIC, self.EM_gain / stack)  if stack>0 else 0#+ np.random.normal(0, self.RN/np.sqrt(int(stack)), cube_detector.shape)                        

            
            # if self.fast is False:
            #     image_stack_only_source[:, OSregions[0] : OSregions[1]] = 0#np.mean([(np.random.gamma(np.random.poisson(source_im_only_source)  , self.EM_gain)) for i in range(int(stack))],axis=0)
            #     image_stack_without_source[:, OSregions[0] : OSregions[1]] = 0#np.mean([(np.random.gamma(np.random.poisson(source_background)  , self.EM_gain)) for i in range(int(stack))],axis=0)
        else:
            # image_stack[:, OSregions[0] : OSregions[1]] = np.mean(np.random.poisson(np.repeat(source_im[:, :, np.newaxis], int(stack), axis=2)),axis=2)
            image_stack[:, OSregions[0] : OSregions[1]] = np.random.poisson(source_im[:, :]*int(stack)) / int(stack)  if int(stack)>0 else 0
            # if self.fast is False:
            #     image_stack_only_source[:, OSregions[0] : OSregions[1]] = np.mean(np.random.poisson(np.repeat(source_im_only_source[:, :, np.newaxis], int(stack), axis=2)),axis=2)
            #     image_stack_without_source[:, OSregions[0] : OSregions[1]] = np.mean(np.random.poisson(np.repeat(source_background[:, :, np.newaxis], int(stack), axis=2)),axis=2)

        # a = (np.where(np.random.rand(int(stack), size[1],OSregions[1]-OSregions[0]) < self.cosmic_ray_loss_per_sec/n_smearing,np.nan,1) * np.array([ (np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand( OSregions[1]-OSregions[0],size[1]).T<self.CIC_charge,dtype=int) , self.EM_gain))  for i in range(int(stack))]))
        # Addition of the phyical image on the 2 overscan regions
        #image += source_im2
        #TODO add this to background too
        if (p_sCIC>0) & (self.fast is False):
            image +=  np.random.gamma( np.array(np.random.rand(size[1], size[0])<p_sCIC,dtype=int) , np.random.randint(1, n_registers, size=image.shape))
            #30%
            image_stack += np.random.gamma( np.array(np.random.rand(size[1], size[0])<int(stack)*p_sCIC,dtype=int) , np.random.randint(1, n_registers, size=image.shape))

        # if self.counting_mode:
        #     a = np.array([ (np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand( OSregions[1]-OSregions[0],size[1]).T<self.CIC_charge,dtype="int32") , self.EM_gain))  for i in range(int(stack))])
        #     cube_stack[:,:, OSregions[0] : OSregions[1]] = a
        #     cube_stack += np.random.gamma( np.array(np.random.rand(int(stack),size[1], size[0])<int(stack)*p_sCIC,dtype=int) , np.random.randint(1, n_registers, size=image.shape)).astype("int32")
            # print(cube_stack.shape)
        #         # addition of pCIC (stil need to add sCIC before EM registers)
        #         prob_pCIC = np.random.rand(size[1], size[0])  # Draw a number prob in [0,1]
        #         image[prob_pCIC < self.CIC_charge] += 1
        #         source_im2_stack[prob_pCIC < p_pCIC*stack] += 1

        #         # EM amp (of source + self.Dark_current + pCIC)
        #         id_nnul = image != 0
        #         image[id_nnul] = np.random.gamma(image[id_nnul], self.EM_gain)
                # Addition of sCIC inside EM registers (ie partially amplified)
        #         prob_sCIC = np.random.rand(size[1], size[0])  # Draw a number prob in [0,1]
        #         id_scic = prob_sCIC < p_sCIC  # sCIC positions
        #         # partial amplification of sCIC
        #         register = np.random.randint(1, n_registers, size=id_scic.sum())  # Draw at which stage of the EM register the electoself.RN is created
        #         image[id_scic] += np.random.exponential(np.power(self.EM_gain, register / n_registers))
            # semaring post EM amp (sgest noise reduction)
            #TODO must add self.smearing for cube! But would take too much time...
        # print(np.ptp(source_im), np.ptp(source_im_only_source))

        if (self.fast is False) & (self.smearing > 0):
            # self.smearing dependant on flux
            #2%
            smearing_kernels = variable_smearing_kernels(image, self.smearing, SmearExpDecrement)
            offsets = np.arange(n_smearing)
            A = dia_matrix((smearing_kernels.reshape((n_smearing, -1)), offsets), shape=(image.size, image.size))

            image = A.dot(image.ravel()).reshape(image.shape)
            image_stack = A.dot(image_stack.ravel()).reshape(image_stack.shape)
        #     if self.readout_time > 0:
        #         # self.smearing dependant on flux
        #         self.smearing_kernels = variable_smearing.smearing_keself.RNels(image.T, self.readout_time, SmearExpDecrement)#.swapaxes(1,2)
        #         offsets = np.arange(n_smearing)
        #         A = dia_matrix((self.smearing_kernels.reshape((n_smearing, -1)), offsets), shape=(image.size, image.size))#.swapaxes(0,1)
        #         image = A.dot(image.ravel()).reshape(image.shape)#.T
        #         image_stack = A.dot(image_stack.ravel()).reshape(image_stack.shape)#.T
        type_ = "int32"
        type_ = "float64"
        readout = np.random.normal(Bias, self.RN, (size[1], size[0]))
        readout_stack = np.random.normal(Bias, self.RN/np.sqrt(int(stack)), (size[1], size[0]))
        if self.counting_mode:
            readout_cube = np.random.normal(Bias, self.RN, (int(stack),size[1], size[0])).astype("int32")
        cr_loss = self.cosmic_ray_loss[self.i]  if hasattr(self.cosmic_ray_loss,"__len__") else self.cosmic_ray_loss
        if (0==1) & (cr_loss>0) : 
            from scipy.ndimage import binary_dilation
            cr = np.zeros(readout.shape)
            cr[np.random.rand(*readout.shape) < cr_loss/8/20]=np.nan
            dilated_mask = binary_dilation(cr, structure=np.ones((8, 20)))
            readout[dilated_mask==1] = np.nan #=dilated_mask

        # TODO Maybe add the readnoise not in this function but in the ETC one
        imaADU_wo_RN = (image * ConversionGain).round().astype(type_)
        imaADU_RN = (readout * ConversionGain).round().astype(type_)
        imaADU = ((image + 1*readout) * ConversionGain).round().astype(type_)
        imaADU_without_source = ((image_without_source + 1*readout) * ConversionGain).round().astype(type_)
        imaADU_source = ((image_only_source + 0*readout) * ConversionGain).round().astype(type_)
        if self.fast is False:
            imaADU_stack_only_source = ((image_stack_only_source + 0*readout_stack ) * ConversionGain).astype(type_)
            imaADU_stack_without_source = ((image_stack_without_source + 1*readout_stack ) * ConversionGain).astype(type_) # same
        if self.counting_mode:
            # imaADU_cube = ((cube_stack + 1*readout_cube) * ConversionGain).round().astype("int32")
            imaADU_stack = image_stack
            imaADU[imaADU<5.5*self.RN] = 0
            imaADU[(imaADU>5.5*self.RN)&(imaADU<Full_well*1000)] = 1
            imaADU[(imaADU>Full_well*1000)] = np.nan
        else:
            imaADU_stack = ((image_stack + 1*readout_stack) * ConversionGain).astype(type_)
            imaADU_cube = imaADU_stack
        imaADU_cube = imaADU_stack
        imaADU[imaADU>Full_well*1000] = np.nan

        return imaADU, imaADU_stack #, imaADU_cube, source_im, source_im_wo_atm, imaADU_stack_only_source, imaADU_without_source, imaADU_stack_without_source, imaADU_source#imaADU_wo_RN, imaADU_RN




def fitswrite(fitsimage, filename, verbose=True, header=None):
    """Write fits image function with different tests
    """

    if type(fitsimage) == np.ndarray:
        try:
            fitsimage = fits.HDUList([fits.PrimaryHDU(fitsimage, header=header)])[0]
        except KeyError as e:
            print(fitsimage)
            print("discarding header due to error: ", e)
            fitsimage = fits.HDUList([fits.PrimaryHDU(fitsimage)])[0]

    if len(filename) == 0:
        filename = tmp_image
        fitsimage.header["NAXIS3"], fitsimage.header["NAXIS1"] = (
            fitsimage.header["NAXIS1"],
            fitsimage.header["NAXIS3"],
        )
        fitsimage.writeto(filename, overwrite=True)
    if hasattr(fitsimage, "header"):
        if "NAXIS3" in fitsimage.header:
            # verboseprint("2D array: Removing NAXIS3 from header...")
            fitsimage.header.remove("NAXIS3")
        if "SKEW" in fitsimage.header:
            fitsimage.header.remove("SKEW")
    elif hasattr(fitsimage[0], "header"):
        if "NAXIS3" in fitsimage[0].header:
            # verboseprint("2D array: Removing NAXIS3 from header...")
            fitsimage[0].header.remove("NAXIS3")
        if "SKEW" in fitsimage[0].header:
            fitsimage[0].header.remove("SKEW")

    elif not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    if not os.path.exists(os.path.dirname(filename)):
        # verboseprint("%s not existing: creating folde..." % (os.path.dirname(filename)))
        os.makedirs(os.path.dirname(filename))
    try:
        fitsimage.writeto(filename, overwrite=True)
    except IOError:
        # verboseprint("Can not write in this repository : " + filename)
        filename = "%s/%s" % (
            os.path.dirname(os.path.dirname(filename)),
            os.path.basename(filename),
        )
        # filename = "/tmp/" + os.path.basename(filename)
        # verboseprint("Instead writing new file in : " + filename) 
        fitsimage.writeto(filename, overwrite=True)
    # verboseprint("Image saved: %s" % (filename))
    return filename
##%%

if __name__ == "__main__":


    instruments = Table.read("../data/Instruments/Instruments.csv")
    instruments.write("../data/Instruments/instruments.csv",overwrite=True)
    self = Observation( instruments=instruments,instrument="SCWI",Signal=0*1e-27,Sky=1e-15,Line_width=600,Slitlength=6, IFS=True)
    # imaADU, imaADU_stack, imaADU_cube, source_im, source_im_wo_atm, imaADU_stack_only_source, imaADU_without_source, imaADU_stack_without_source, imaADU_source = self.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Baseline Spectra", size=[500, 100], OSregions=[0, 500], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=False,atmlambda=False,fraction_lya=0.05,sky_lines=True)#,Altitude=3
    imaADU, imaADU_stack, imaADU_cube, source_im, source_im_wo_atm, imaADU_stack_only_source, imaADU_without_source, imaADU_stack_without_source, imaADU_source = self.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="cube", size=[500, 100], OSregions=[0, 500], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=False,atmlambda=False,fraction_lya=0.05,sky_lines=True)#,Altitude=3
    if 1==1:
        arrays = {
        'imaADU': imaADU,
        'imaADU_stack': imaADU_stack,
        'imaADU_cube': imaADU_cube,
        'source_im': source_im,
        'source_im_wo_atm': source_im_wo_atm,
        'imaADU_stack_only_source': imaADU_stack_only_source,
        'imaADU_without_source': imaADU_without_source,
        'imaADU_stack_without_source': imaADU_stack_without_source,
        'imaADU_source': imaADU_source
        }
        fig, axes = plt.subplots(3, 3, figsize=(15, 8),sharex=True, sharey=True)
        fig.suptitle('Visualisation des Arrays')

        # Parcours de chaque tableau et affichage
        for ax, (title, data) in zip(axes.flatten(), arrays.items()):
            im = ax.imshow(data, cmap='viridis', aspect='auto')  # Choix du colormap
            ax.set_title(title)
            # ax.axis('off')  # Enlève les axes pour une présentation plus claire
            fig.colorbar(im, ax=ax)  # Ajoute une colorbar à chaque subplot

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Ajuste le layout
        plt.show()


# %load_ext line_profiler
# %lprun -f Observation  Observation()#.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Slit", size=[100, 100], OSregions=[0, 100], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=True,atmlambda=True,fraction_lya=0.05)

## %%
# %lprun -u 1e-1 -T /tmp/initilize.py -s -r -f  Observation.initilize  Observation(exposure_time=np.linspace(50,1500,50))
# %lprun -u 1e-1 -T /tmp/interpolate_optimal_threshold.py -s -r -f  Observation.interpolate_optimal_threshold  Observation(exposure_time=np.linspace(50,1500,50),counting_mode=True,plot_=False).interpolate_optimal_threshold()
# %lprun -u 1e-1 -T /tmp/PlotNoise.py -s -r -f  Observation.PlotNoise  Observation(exposure_time=np.linspace(50,1500,50)).PlotNoise()
# %lprun -u 1e-1 -T /tmp/SimulateFIREBallemCCDImage.py -s -r -f  Observation.SimulateFIREBallemCCDImage  Observation(exposure_time=np.linspace(50,1500,50)).SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Slit", size=[100, 100], OSregions=[0, 100], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=True,atmlambda=True,fraction_lya=0.05)
## %%


# %%
