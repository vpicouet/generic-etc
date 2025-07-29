# %%
from astropy.table import Table
import warnings
warnings.filterwarnings("ignore")
from ipywidgets import Button, Layout, jslink, IntText, IntSlider, interactive, interact, HBox, Layout, VBox
from astropy.modeling.functional_models import Gaussian2D, Gaussian1D
import os, glob
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
import cmocean
from scipy.stats import norm
from matplotlib.offsetbox import AnchoredText
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable

np.seterr(invalid='ignore')
 

from scipy.ndimage import map_coordinates


gaus = lambda x, a, xo, sigma, offset: a ** 2  * np.exp(-(x - xo)**2 / (2 * sigma**2)) + offset
n2,n1=100,500


def download(url, file=""):
    """Download a file
    """
    from tqdm import tqdm  # , tqdm_gui
    import requests

    try:
        response = requests.get(url, stream=True)
    except requests.exceptions.RequestException as e:
        print(e)
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
            print("ERROR, something went wrong")
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


# def resample_cube(input_cube, wave_range_nm, spatial_extent_arcsec, output_shape, phys=False):
#     """
#     Resample un datacube avec interpolation/extrapolation selon longueur d'onde et spatial.
#     """
#     header = input_cube.header
#     input_cube=input_cube.data
#     wave_start_out, wave_end_out = wave_range_nm
#     input_nz, input_nx, input_ny = input_cube.data.shape
#     # print(input_cube.data.shape)
#     cunit3 = header.get('CUNIT3', 'nm').strip().lower()
#     # Définition du facteur de conversion
#     unit_conversion = {
#         'm': 1e9,       # mètres → nanomètres
#         'cm': 1e7,      # centimètres → nanomètres
#         'mm': 1e6,      # millimètres → nanomètres
#         'um': 1e3,      # micromètres → nanomètres
#         'nm': 1,        # nanomètres (aucun changement)
#         'angstrom': 0.1 # Ångströms → nanomètres
#     }
#     factor = unit_conversion.get(cunit3, 1)  # Valeur par défaut = 1 (nm)
#     # Calcul correct des longueurs d'onde en tenant compte de l'unité
#     wave_start_in = (header['CRVAL3'] + (1 - header['CRPIX3']) * header['CDELT3']) * factor
#     wave_end_in = (header['CRVAL3'] + (input_nz - header['CRPIX3']) * header['CDELT3']) * factor

#     print(wave_start_in,wave_end_in)

def resample_cube(input_hdu, wave_range_nm, spatial_extent_arcsec, output_shape, phys=False):
    """
    Resample a datacube with interpolation/extrapolation in wavelength and spatial directions.
    """
    header = input_hdu.header
    input_cube=input_hdu.data
    wave_start_out, wave_end_out = wave_range_nm
    input_nz, input_nx, input_ny = input_cube.shape

    # Get wavelength axis information from header
    crval = header['CRVAL3']
    cdelt = header['CDELT3']
    crpix = header['CRPIX3']
    cunit = header.get('CUNIT3', 'm')  # Default to nm if missing

    try:
        unit = u.Unit(cunit)
    except ValueError:
        raise ValueError(f"Unknown unit in CUNIT3: {cunit}")

    pix = np.arange(input_nz)
    wave_axis = (crval + (pix + 1 - crpix) * cdelt) * unit
    # print(wave_axis)
    wave_axis_nm = wave_axis.to(u.nm).value

    wave_start_in = wave_axis_nm[0]
    wave_end_in = wave_axis_nm[-1]
    # print(wave_start_in,wave_end_in)
    spatial_radius_in_x = header['CDELT1'] * (input_nx - 1) / 2.0 * 3600  # Convert from degrees to arcseconds
    spatial_radius_in_y = header['CDELT2'] * (input_ny - 1) / 2.0 * 3600  # Convert from degrees to arcseconds
    spatial_radius_out = spatial_extent_arcsec / 2.0

    if phys==False:
        # in this case I just keep the exact same cube and act as if it was exactly the size of the FOV
        input_nx, input_ny, input_nz = input_cube.shape
        x_out = np.linspace(0, input_nx - 1, output_shape[0])
        y_out = np.linspace(0, input_ny - 1, output_shape[1])
        z_out = np.linspace(0, input_nz - 1, output_shape[2])

        x_grid, y_grid, z_grid = np.meshgrid(x_out, y_out, z_out, indexing='ij')
        resampled_cube = map_coordinates(input_cube, [x_grid, y_grid, z_grid], order=1, mode='nearest')
        return resampled_cube
    else:
        # in this case I just keep the exact same cube and act as if it was exactly the size of the FOV
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



def load_instruments(sheet_id= "1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU", sheet_name= "instruments.csv", database="Online DB"):
    """
    Load instruments data from a Google Sheet or local database.

    Parameters:
    - sheet_id (str): The Google Sheet ID.
    - sheet_name (str): The name of the sheet to load.
    - database (str): Either "Online DB" or "Local DB".

    Returns:
    - instruments (Table): The loaded instruments data.
    - database (str): The database source used.
    """
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    if database == "Online DB":
        try:
            instruments =   Table.from_pandas(pd.read_csv(url)).filled(np.nan)
            database = "Online DB"
        except Exception as e:
            print(e)
            instruments = Table.from_pandas(pd.read_excel("../instruments.xlsx"))
            instruments = instruments[instruments.colnames]
            for col_name in instruments.colnames[3:]:
                instruments[col_name] = to_float_or_nan(instruments[col_name])
            database = "Local DB"
    else:
        instruments = Table.from_pandas(pd.read_excel("../instruments.xlsx"))
        instruments = instruments[instruments.colnames]
        for col_name in instruments.colnames[3:]:
            instruments[col_name] = to_float_or_nan(instruments[col_name])
        database = "Local DB"
    return instruments, database

# instruments, database = load_instruments(sheet_id= "1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU", sheet_name= "instruments.csv")
# instruments_SNR, database_SNR = load_instruments(sheet_name="SNR")


def dark_plateau(T=-80,a=0.23362821,plateau=0.44384484):
    T = T+273.15
    E = 1.1557 - (7.021*1e-4*T**2)/(1108+T)
    k = 8.62e-5
    return 2.55*1e15*a*3600*0.0013**2*T**(3/2)*np.exp(-E/(2*k*T))+plateau



def mostFrequent(arr):
    n=len(arr)
    # Insert all elements in Hash.
    Hash = dict()
    for i in range(n):
        if arr[i] in Hash.keys():
            Hash[arr[i]] += 1
        else:
            Hash[arr[i]] = 1
    # find the max frequency
    max_count = 0
    res = -1
    for i in Hash:
        if (max_count < Hash[i]):
            res = i
            max_count = Hash[i]
         
    return res





class ExposureTimeCalulator(widgets.HBox):
    @initializer
    def __init__(self, instruments=None,database=None, instrument="FIREBall-2 2025",x_axis='exposure_time', time_max = 2005,SNR_res="per Res elem" ,spectrograph=True, **kwargs):#, Atmosphere=0.5, Throughput=0.13*0.9, follow_temp=False, acquisition_time=1, Sky=4, Signal=24, EM_gain=1400,RN=109,CIC_charge=0.005, Dark_current=0.08,readout_time=1.5,counting_mode=False,smearing=0.7,extra_background=0,temperature=-100,PSF_RMS_mask=2.5,PSF_RMS_det = 3.5,QE=0.45,cosmic_ray_loss_per_sec=0.005,
        """
        Generate an ETC app containing multiple widghet that allow to change the ETC parameters
        as well as plotting the result (e- and noise budget, limiting flux, SNR) in terms of the different parameters.
        """
        super().__init__()
        self.instruments=instruments
        args, _, _, locals_ = inspect.getargvalues(inspect.currentframe())
        self.Redshift=0
        for i in range(len(instruments)):
            setattr(self, instruments["Charact."][i], instruments[instrument][i] ) if  isinstance(type(instruments[instrument][i]), (int, float, complex,np.float64))  else setattr(self, instruments["Charact."][i], float(instruments[instrument][i]) )#.replace("%","")
        exposure_time=np.logspace(0,np.log10(time_max))
        self.instruments_dict = {name: {key: val for key, val in zip(instruments["Charact."][:], instruments[name][:]) if not isinstance(key, np.ma.core.MaskedConstant) and not isinstance(val, np.ma.core.MaskedConstant)} for name in instruments.colnames[3:]}

        self.output = widgets.Output()
        time=exposure_time
        self.time =  float(instruments[instrument][instruments["Charact."]=="exposure_time"][0])
        i = np.argmin(abs(self.time - exposure_time))
        self.arg=i
        # print(self.arg,i,time[i],time[i],self.time)
        IFS = False if self.dimensions==2 else True

        self.new = Observation(instruments=instruments, instrument=instrument,Redshift=self.Redshift,Throughput_FWHM=self.Throughput_FWHM, smearing=self.smearing,counting_mode=False,exposure_time=exposure_time,Sky=self.Sky,acquisition_time=self.acquisition_time,Signal=self.Signal,EM_gain=self.EM_gain,RN=self.RN, CIC_charge=self.CIC_charge, Dark_current=self.Dark_current,PSF_RMS_mask=self.PSF_RMS_mask,PSF_RMS_det=self.PSF_RMS_det,QE=self.QE, extra_background=self.extra_background,
            Collecting_area=self.Collecting_area, pixel_scale=self.pixel_scale, Throughput=self.Throughput, Spectral_resolution=self.Spectral_resolution, Slitwidth=self.Slitwidth, dispersion=self.dispersion,
            Size_source=self.Size_source,Line_width=self.Line_width,wavelength=self.wavelength, Atmosphere=self.Atmosphere, pixel_size=self.pixel_size,cosmic_ray_loss_per_sec=self.cosmic_ray_loss_per_sec,readout_time=self.readout_time, Slitlength=self.Slitlength,i=self.arg,Δλ=0,Δx=0,IFS=IFS,SNR_res=SNR_res,spectrograph=spectrograph)


        self.x = exposure_time
        
        
        style={}
        width = '400px'
        width = '500px'
        small = '247px'
        small = '230px'
        psmall = '186px'
        vsmall = '147px'
        c_update=False
        # ALL THE DIFFERENT VARIATION WITH TEMPERATURE
        self.smearing_poly = np.poly1d([-0.0306087, -2.2226087])#np.poly1d([-0.0453913, -3.5573913])
        # self.dark_poly = np.poly1d([2.13640462e-05, 7.83596239e-03, 9.57682651e-01, 3.86154296e+01])#with plateau
        # self.dark_poly = np.poly1d([0.07127906, 6.83562573]) #without plateau# does to low because only calibrated down to -100 but then it kinda saturated. maybe because of CIC?
        self.dark_poly = dark_plateau
        # self.CIC_poly = np.poly1d([1.57925408e-05, 2.80396270e-03, 1.34276224e-01]) #without plateau# does to low because only calibrated down to -100 but then it kinda saturated. maybe because of CIC?


        # self.other_options =  ['exposure_time','acquisition_time',"Signal","RN","Dark_current" ,'Sky',"readout_time","PSF_RMS_det","QE","cosmic_ray_loss_per_sec","Throughput","Atmosphere","lambda_stack","Size_source",'Collecting_area',"Δx","Δλ"]#[:-3]
        self.other_options =  [
                               "------DECTECTOR PERFORMANCE","QE","RN","Dark_current" ,"cosmic_ray_loss_per_sec","EM_gain","extra_background","CIC_charge", #,"pixel_size"
                               "-------------OBSERVED SOURCE","Signal",'Sky',"Size_source","Line_width", #,"Redshift"
                               "--------OBSERVATION STRATEGY", "Atmosphere",'exposure_time','acquisition_time',"readout_time","lambda_stack",  "wavelength", #,"Δx","Δλ"
                               "-----------INSTRUMENT DESIGN",'Collecting_area',"pixel_scale","Throughput","PSF_RMS_mask","PSF_RMS_det",
                               "---------SPECTROGRAPH DESIGN","Spectral_resolution","Slitwidth","Slitlength","dispersion"] #,"Bandwidth"
        self.other_options_imager =  [
                               "------DECTECTOR PERFORMANCE","QE","RN","Dark_current" ,"cosmic_ray_loss_per_sec", #,"pixel_size"
                               "-------------OBSERVED SOURCE","Signal",'Sky',"Size_source",#,"Redshift"
                               "--------OBSERVATION STRATEGY", "Atmosphere",'exposure_time','acquisition_time',"readout_time",
                               "-----------INSTRUMENT DESIGN",'Collecting_area',"pixel_scale","Throughput","PSF_RMS_det"]
        self.fb_options_no_temp = self.other_options #+ ["----------AMPLIFIED DECTECTOR","smearing"]
        self.fb_options = self.fb_options_no_temp + ["temperature"]

        self.instrument = widgets.Dropdown(options=instruments.colnames[3:],value=instrument,description="Instrument", layout=Layout(width=small),description_tooltip="Instrument characteristics",continuous_update=c_update)
        # print(self.instruments_dict[self.instrument.value]["dispersion"],type(self.instruments_dict[self.instrument.value]["dispersion"]))
        self.spectro = False if np.isnan(self.instruments_dict[self.instrument.value]["dispersion"]) else True

        self.ylog = widgets.Checkbox(value=False,description='ylog',disabled=False,style =dict(description_width='initial'), layout=Layout(width='147px'),description_tooltip="Use this check box to use a log scale for the y axis",continuous_update=c_update)
        self.yscale="linear"
        self.xlog = widgets.Checkbox(value=True,description='xlog',disabled=False,style =dict(description_width='initial'), layout=Layout(width='147px'),description_tooltip="Use this check box to use a log scale for the x axis",continuous_update=c_update)
        # self.SNR_res = widgets.Checkbox(value=SNR_res,description='SNR(Res)',disabled=False, style =dict(description_width='initial'),layout=Layout(width='147px'),description_tooltip="Use this check box to plot the SNR per pixel or per element resolution",continuous_update=c_update)
        
        self.SNR_res = widgets.Dropdown(value=SNR_res,description='SNR',options=["per pix","per Res elem","per Source"],disabled=False, style =dict(description_width='initial'),layout=Layout(width='147px'),description_tooltip="Use this check box to plot the SNR per pixel or per element resolution",continuous_update=c_update)
        #,"λPix/xRes","xPix/λRes"
        
        self.IFS = widgets.Checkbox(value=IFS,description='IFS',disabled=False, layout=Layout(width='147px'),description_tooltip="Check this box if the instrument is an integral field spectrograph (not just a single slit of fiber)",continuous_update=c_update)
        # self.source_im = widgets.Checkbox(value=False,description='Source',disabled=False, layout=Layout(width='179px'),description_tooltip="Check this box to image the source",continuous_update=c_update)
        self.source_im = widgets.Dropdown(options=["Sim image","Source","Convolved source","SNR"],value="Sim image",description='Type',disabled=False, layout=Layout(width='179px'),description_tooltip="Type of the images. Either the source, the simulated images or the SNR",continuous_update=c_update)
        #,"SNR"
        self.spectrograph = widgets.Checkbox(value=self.spectro,description='Spectro',disabled=False, visible=False,layout=Layout(width='167px'),description_tooltip="Check this box if the instrument is not just an imager (= has a dispersive element)",continuous_update=c_update)
        # self.IFS_value = ("IFS" if IFS else "Spectro") is self.spectro else "Imager"
        # self.IFS_value = "IFS" if IFS else ("Spectro" if self.spectro else "Imager")
        # self.IFS = widgets.Dropdown(value="Spectro",options=["Imager","Spectro","IFS"],description='Type',disabled=False, layout=Layout(width='147px'),description_tooltip="IFS is instrument is an integral field spectrograph (not just a single slit of fiber)",continuous_update=True)
       
        self.Signal = widgets.FloatLogSlider(min=-21, max=-9 ,value=self.Signal,description='Source brightness', style =dict(description_width='initial'), layout=Layout(width=width),description_tooltip="Flux of the diffuse source in ergs/cm2/s/arcsec2/Å.",continuous_update=c_update)
        self.Sky = widgets.FloatLogSlider( min=-23, max=-15,value=self.Sky,base=10, style =dict(description_width='initial'), layout=Layout(width=width),description='Sky    brightness',description_tooltip="Level of sky background illumination (zodiacal and galactic) in ergs/cm2/s/arcsec2/Å ",continuous_update=c_update)
        self.acquisition_time = widgets.FloatLogSlider( min=-1, max=3,value=self.acquisition_time,style =style ,base=10,layout=Layout(width=width),description='Taq (h)',description_tooltip="Total acquisition time [hours]",continuous_update=c_update)
        self.exposure = widgets.FloatRangeSlider( min=0, max=time_max,value=(self.readout_time,self.time),style = style, layout=Layout(width=width),description='Rd/Exp time',step=0.1,readout_format='.0f',description_tooltip="Readout time and exposure time [seconds]",continuous_update=c_update)

        self.fwhm = widgets.FloatRangeSlider( min=0.01, max=6,value=(self.PSF_RMS_mask,self.PSF_RMS_det),style = style, layout=Layout(width=width),description='Mask/det σ',step=0.01,readout_format='.2f',description_tooltip="Spatial resolution in arcseconds respectively at the mask and detector level. To be multiplied by 2.35 to have the FWHM.",continuous_update=c_update)
        self.RN = widgets.FloatSlider( min=0.01, max=120,value=self.RN, style = style, step=0.1, layout=Layout(width=width),description='Read noise',description_tooltip="Detector readout noise in electrons/pixel",continuous_update=c_update)
        self.QE = widgets.FloatSlider( min=0.01, max=1,value=self.QE,style = style, layout=Layout(width=width),description='QE',step=0.01,readout_format='.2f',description_tooltip="Detector quantum efficiency",continuous_update=c_update)
        self.Dark_current = widgets.FloatSlider( min=0, max=50,value=self.Dark_current, style = style, layout=Layout(width=width),description='Dark current',step=0.0011,readout_format='.2f',description_tooltip="Detector dark current [e-/pix/hour]",continuous_update=c_update)

        self.extra_background = widgets.FloatSlider( min=0, max=200,value=self.extra_background,style = style, layout=Layout(width=width),description='Extra bckgnd',step=0.2,readout_format='.1f',description_tooltip="Additional background on the detector [e-/pix/hour]",continuous_update=c_update)
        self.EM_gain = widgets.IntSlider( min=1, max=3500,value=self.EM_gain, style = style, layout=Layout(width=width),description='EM gain',description_tooltip="EMCCD amplification gain in e-/e-",continuous_update=c_update)
        self.CIC_charge = widgets.FloatSlider( min=0, max=0.07,value=self.CIC_charge,style = style, layout=Layout(width=width),description='CIC charge',step=0.001,readout_format='.3f',description_tooltip="EMCCD spurious charges due to amplification in electrons [e-/pix]",continuous_update=c_update)
        self.follow_temp = widgets.Checkbox(value=False,description='Temp',disabled=False, layout=Layout(width=vsmall),description_tooltip="Check this box to force charge transfer efficiency and dark current levels to be fixed by the temperature widget. Interesting feature to optimize EMCCD temperature.",continuous_update=c_update)
        self.counting_mode = widgets.Checkbox(value=False,description='γ-Threshold',disabled=False, layout=Layout(width=psmall),description_tooltip="Check this box to apply thresholding photon counting processing. The efficiency of this process is determined by the gain, read noise, smearing, flux.",continuous_update=c_update)
        self.temperature = widgets.FloatSlider( min=-120, max=-85,value=-115, style = style,description=r'Temp (C)',step=0.1, layout=Layout(width=width),description_tooltip="EMCCD's Temperature in Celcius degrees: determines its charge transfer efficiency and dark current rate.",continuous_update=c_update)
        self.smearing = widgets.FloatSlider( min=0, max=self.smearing_poly(-120),value=self.smearing, layout=Layout(width=width),description='Smearing',step=0.01,description_tooltip="Smearing length of the EMCCD (exponential length in pixels). This length, representing the charge transfer efficiency is fixed by the temperature when the Temp checkbox is checked.",continuous_update=c_update)
       
        self.Collecting_area = widgets.FloatLogSlider( min=-2, max=3,value=self.Collecting_area, style =style,base=10, layout=Layout(width=width),description='Area',description_tooltip="Collecting area of the instrument in square meter",continuous_update=c_update)
        self.pixel_scale = widgets.FloatSlider( min=0.01, max=5,value=self.pixel_scale,base=10, style =style, layout=Layout(width=width),description='Pixel scale',description_tooltip="Pixel plate scale in  ''/pix",continuous_update=c_update)
        self.Throughput = widgets.FloatSlider( min=0.01, max=1,value=self.Throughput,base=10, style =style, layout=Layout(width=width),description='Throughput',description_tooltip="Instrument throughput at effective wavelength (not accounting for detector quantum efficiency and atmospheric transmission)",continuous_update=c_update)

        self.database = widgets.Dropdown(options=["Online DB","Local DB"],value=database,description="", layout=Layout(width='90px'),description_tooltip="Instrument characteristics",continuous_update=c_update, style =dict(description_width='initial'))
        self.interpolation = widgets.Dropdown(options=["None", 'gaussian', 'none', 'nearest', 'bilinear', 'bicubic', 'spline16','spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric','catrom', 'bessel', 'mitchell', 'sinc', 'lanczos'],value="None",description="Interpolation", layout=Layout(width='350px'),description_tooltip="Interpolation method in the the imshow method",continuous_update=c_update, style =dict(description_width='initial'))


        self.dependencies = widgets.SelectMultiple(options=['dispersion=10*wavelength/Spectral_resolution/2','pixel_scale=2.35*PSF_RMS_det/2', 'CIC_charge=EM_gain/2'],value=[],rows=2,description='Dep.',disabled=False)

        
        self.Atmosphere = widgets.FloatSlider( min=0.1, max=1,value=self.Atmosphere,base=10, style =style, layout=Layout(width=width),description='Atmosphere',description_tooltip="Atmospheric transmission",continuous_update=c_update)
        self.pixel_size = widgets.FloatSlider( min=2, max=40,value=self.pixel_size,base=10, style =style, layout=Layout(width=width),description='Pix size',description_tooltip="Pixel size in microns",continuous_update=c_update)
        self.pixel_size.layout.visibility = 'hidden'  

        self.Size_source = widgets.FloatSlider( min=0.01, max=100,value=self.Size_source,base=10, style =style, layout=Layout(width=width),description='σ Source',description_tooltip="Spatial extension of the source in arcseconds",continuous_update=c_update)
        self.Line_width = widgets.FloatSlider( min=0.01, max=100,value=self.Line_width,base=10, style =style, layout=Layout(width=width),description='Eq width (Å)',description_tooltip="Spectral extension of the source/emission line in Å ",continuous_update=c_update)


        self.wavelength = widgets.FloatSlider( min=50, max=1000,value=self.wavelength,base=10, style =style, step=0.1, layout=Layout(width=width),description='Observed λ',description_tooltip="Oberved λ in Å (only used for conversions)",continuous_update=c_update)
        self.Δλ = widgets.FloatSlider( min=-250, max=250,value=-10,base=10, style =style, layout=Layout(width='400px'),description='Δλ',description_tooltip="Distance to the emission line being analyzed in pixels",continuous_update=c_update)
        self.Δx = widgets.FloatSlider( min=-50, max=50,value=0,base=10, style =style, layout=Layout(width='400px'),description='Δx',description_tooltip="Distance to the source being analyzed in pixels",continuous_update=c_update)
        self.Throughput_FWHM = widgets.FloatLogSlider( min=1, max=5,value=np.log10(self.Throughput_FWHM),base=10, style =style, layout=Layout(width='400px'),description='TH FWHM',description_tooltip="Instrument throughput FWHM in  Å. If an instrument-specific λ-dependent throughput is added to the repository, this csv file will be used.")#5.57e-18
        self.Redshift = widgets.FloatSlider( min=0.01, max=10,value=self.Redshift,base=10, step=0.01, style =style, layout=Layout(width='400px'),description='Redshift',description_tooltip="Redshift only considered to shift blackbody spectra or QSO/Star spectra. Flux will stay based on source observed surface brightness. Changes also the kpc size",continuous_update=c_update)
        # TODO verify that this stack is well taken into account into the CNR computation
        self.lambda_stack = widgets.IntSlider( min=1, max=200,value=1, layout=Layout(width=width),description='λ Stack (pix)',step=0.1,description_tooltip="Number of spectral slices used to stack cube (spectral pix)",continuous_update=c_update)

        self.Spectral_resolution = widgets.IntSlider( min=90, max=10000,value=self.Spectral_resolution,base=10, style =style, layout=Layout(width=width),description='R (λ/dλ) ',description_tooltip="Instrument spectral resolution λ/dλ",continuous_update=c_update)
        # self.Slitwidth = widgets.FloatSlider( min=0.1, max=600,value=self.Slitwidth,base=10, style =style, layout=Layout(width=width),description='Slit ["]',description_tooltip="Width of the slit [''] ")
        self.SlitDims = widgets.FloatRangeSlider( min=0.001, max=600,value=(self.Slitwidth,self.Slitlength),base=10, step=0.001, style =style, layout=Layout(width=width),description='Slit dims["]',description_tooltip="Width and length of the slit [''] (put same value for fibers) ")
        
        self.minmax = widgets.FloatRangeSlider( min=0, max=1,value=(0,1),base=10, step=0.001, style =style, layout=Layout(width=width),description='Vmin/Vmax',description_tooltip="Color map min/max",continuous_update=c_update)
        
        self.dispersion = widgets.FloatSlider( min=0.01, max=5,value=self.dispersion,base=10, style =style, step=0.001, readout_format='.2f', layout=Layout(width=width),description='Dispersion',description_tooltip="Dispersion at the detector Å/pix ",continuous_update=c_update)

        self.cosmic_ray_loss_per_sec = widgets.FloatSlider( min=0, max=0.1,value=self.cosmic_ray_loss_per_sec,base=10, style =style,step=0.0001,readout_format='.4f', layout=Layout(width=width),description='CR loss',description_tooltip="Cosmic ray loss per second. eg. 0.01 would mean that 1 sec image looses 1% pixels due to cosmic rays",continuous_update=c_update)

        self.gals = ["Rest-frame: COSMOS " + os.path.basename(f).replace(".txt","") for f in glob.glob("../data/Spectra/GAL_COSMOS_SED/*.txt")]
        self.QSOs = ["Rest-frame: Salvato " + os.path.basename(f).replace(".txt","") for f in glob.glob("../data/Spectra/QSO_SALVATO2015/*.txt")]
        self.spectra_options = ["Observed-frame: Baseline Spectra", "Rest-frame: Baseline Lyα (1216Å)", "Rest-frame: Baseline CIV (1549Å)", "Rest-frame: Baseline OVI (1033Å)", "Rest-frame: Baseline CIII (1908Å)","galaxy_disk_cube-resampled_phys","galaxy_disk_cube-resampled","galaxy_and_cgm_cube-resampled_phys","galaxy_and_cgm_cube-resampled","CGM_cube-resampled_phys","CGM_cube-resampled",]  + self.gals + self.QSOs + ["Rest-frame: Blackbody 5900 K (09V)","Rest-frame: Blackbody 1500 K (BOV)","Rest-frame: Blackbody 9000 K (B3V)","Rest-frame: Blackbody 480 K (AOV)","Rest-frame: Blackbody 8810 K (A2V)","Rest-frame: Blackbody 8160 K (A5V)","Rest-frame: Blackbody 7020 K (FOV)","Rest-frame: Blackbody 6750 K (F2V)","Rest-frame: Blackbody 6530 K (F5V)","Rest-frame: Blackbody 930 K (GOV)","Rest-frame: Blackbody 5830 K (G2V)","Rest-frame: Blackbody 5560 K (G5V)","Rest-frame: Blackbody 240 K (KOV)","Rest-frame: Blackbody 5010 K (K2V)","Rest-frame: Blackbody 4560 K (K4V)","Rest-frame: Blackbody 4340 K (K5V)","Rest-frame: Blackbody 4040 K (K7V)","Rest-frame: Blackbody 3800 K (MOV)","Rest-frame: Blackbody 3530 K (M2V)","Rest-frame: Blackbody 3380 K (M3V)","Rest-frame: Blackbody 3180 K (M4V)","Rest-frame: Blackbody 3030 K (M5V)","Rest-frame: Blackbody 2850 K (M6V)"] + ["Observed-frame: UVSpectra 1538p477 NUV~16.6","Observed-frame: UVSpectra 1821p643 NUV~14",'Observed-frame: UVSpectra 0044p030 NUV~16.5',"Observed-frame: UVSpectra mrk509","Observed-frame: UVSpectra 2344p092","Observed-frame: UVSpectra 1637p574","Observed-frame: UVSpectra 1115p080","Observed-frame: UVSpectra 0414m060","Observed-frame: UVSpectra 0115p027","Observed-frame: UVSpectra 2251p113","Observed-frame: UVSpectra 2201p315","Observed-frame: UVSpectra 1928p738","Observed-frame: UVSpectra 1700p518","cube 10 kpc galaxy + Lya em CGM+Filament","cube 30 kpc galaxy + Lya em CGM+Filament","cube 100 kpc galaxy + Lya em CGM+Filament","lya_cube_merged_with_artificial_source_CU_1pc-resampled_phys","lya_cube_merged_with_artificial_source_CU_1pc-resampled","cube_01-resampled_phys","cube_01-resampled"]#,"lya_cube_merged_with_artificial_source_CU_1pc_map","lya_cube_merged_with_artificial_source_CU_1pc_resampled"] 
        # self.spectra     = widgets.Dropdown(options=self.spectra_options, layout=Layout(width='350px'),description='Spectra',value="galaxy_and_cgm_cube-resampled_phys",continuous_update=c_update)#Observed-frame: Baseline Spectra   "lya_cube_merged_with_artificial_source_CU_1pc-remap"
        self.spectra     = widgets.Dropdown(options=self.spectra_options, layout=Layout(width='350px'),description='Spectra',value="Observed-frame: Baseline Spectra",continuous_update=c_update)#Observed-frame: Baseline Spectra   "lya_cube_merged_with_artificial_source_CU_1pc-remap"
        self.units       = widgets.Dropdown(options=["ADU/frame","e-/frame","photons/frame","e-/hour","photons/hour","e-/second","photons/second"], layout=Layout(width='350px'),description='Units',value="ADU/frame")# TODO add ergs/cm2/... "amplified e-/frame","amplified e-/hour",

        # widgets.Dropdown(options=["S2: 0.053 ADU/e-, FW=5.6 KADU","S2_hdr: 0.97 ADU/e-, FW=52 KADU","1: 0.02 ADU/e-, FW=2.1 KADU","1': 0.4 ADU/e-, FW=40 KADU","2: 0.04 ADU/e-, FW=4.7 KADU","2018: 0.5 ADU/e-, FW=56 KADU","2022: 0.2 ADU/e-, FW=22 KADU","2023_noOS: 0.04 ADU/e-, FW=39 KADU"], layout=Layout(width='350px'),description='RO seq',value="S2_hdr: 0.97 ADU/e-, FW=52 KADU",continuous_update=c_update)
        self.QElambda = widgets.Checkbox(value=True,description='Throughput(λ)',disabled=False,tooltip="Check this box to apply λ QE dependancy",layout=Layout(width="217px"))
        self.atmlambda = widgets.Checkbox(value=True,description='atm(λ)',disabled=False,tooltip="Check this box to apply λ atm transmission dependancy",layout=Layout(width="217px"))
        self.sky_lines = widgets.Checkbox(value=True,description='Sky lines',disabled=False,tooltip="Check this box to add sky emission lines",layout=Layout(width="217px"))
        self.test = widgets.Checkbox(value=True,description='F',disabled=False,tooltip="Check this box to use the method where we compute all flux and then divide by pixels or compute directly the flux per pixel",layout=Layout(width="100px"))
        # self.test.layout.visibility = 'hidden'  


        self.fraction_lya = widgets.FloatSlider( min=0, max=0.2,value=0.05,style = style, layout=Layout(width='217px'),description='Lya fraction',step=0.001,readout_format='.2f',tooltip="Fraction of E(Lya)/E(NUV)")
        self.fraction_lya.layout.visibility = 'hidden'  

        self.save_plot_button = widgets.Button(description="Save Plot", layout=Layout(width='auto'))
        self.save_data_button = widgets.Button(description="Save Data", layout=Layout(width='auto'))


        self.reset = widgets.Button(value=False,description='↺',disabled=False,button_style='', layout=Layout(width="30px")) 


        self.change = widgets.Checkbox(value=True,description='change',disabled=False, layout=Layout(width='147px'),description_tooltip="Use this check box to use a log scale for the x axis")
        self.change.layout.visibility = 'hidden'  
        self.spectrograph.layout.visibility = 'hidden'  
        # self.lambda_stack = widgets.FloatSlider( min=10*self.wavelength/self.Spectral_resolution, max=self.Bandwidth*10,value=self.Bandwidth, layout=Layout(width=width),description='λ width [Å]',step=0.1,description_tooltip="Wavelength range used to stack signal. Min = 1 resolution element = 5pixels = 1Å, Max = total spectra =bandwidth = 200Å")   
       
        if "FIREBall" in instrument:
            options = self.fb_options if self.follow_temp.value else self.fb_options_no_temp
        else:
            options = self.other_options
        self.x_axis=widgets.Dropdown(options=options,value=self.x_axis,description='X axis', layout=Layout(width=small),description_tooltip="Variable used to analyze the evolution of the SNR.")

        self.follow_temp.layout.visibility = 'hidden'  
        # if ("FIREBall" not in instrument) & ("SCWI" not in instrument) :
        if self.EM_gain.value == 1 :
            self.counting_mode.layout.visibility = 'hidden'  

        self.smearing.layout.visibility = 'hidden'
        self.temperature.layout.visibility = 'hidden'
        if self.follow_temp.value:
            self.Dark_current.value = self.dark_poly(self.temperature.value)#10**
            self.smearing.value = self.smearing_poly(self.temperature.value)
            # self.CIC_charge.value = self.CIC_poly(self.temperature.value)
  

        # print(self.new.QE,self.new.N_images_true)
        # self.Signal_el = self.new.Signal_el

        wids = widgets.interactive(self.update,x_axis=self.x_axis,xlog=self.xlog,log=self.ylog, SNR_res=self.SNR_res,smearing=self.smearing,counting_mode=self.counting_mode,exposure=self.exposure,Sky=self.Sky,acquisition_time=self.acquisition_time,Signal=self.Signal,EM_gain=self.EM_gain,RN=self.RN, CIC_charge=self.CIC_charge, Dark_current=self.Dark_current,temperature=self.temperature,follow_temp=self.follow_temp,fwhm = self.fwhm,QE=self.QE, extra_background=self.extra_background,
            Collecting_area=self.Collecting_area, pixel_scale=self.pixel_scale, Throughput=self.Throughput, Spectral_resolution=self.Spectral_resolution, SlitDims=self.SlitDims, dispersion=self.dispersion,
            Size_source=self.Size_source,Line_width=self.Line_width,wavelength=self.wavelength,Δλ=self.Δλ,Δx=self.Δx, Atmosphere=self.Atmosphere, pixel_size=self.pixel_size,cosmic_ray_loss_per_sec=self.cosmic_ray_loss_per_sec,lambda_stack=self.lambda_stack,#change=self.change,
            spectra=self.spectra,units=self.units,Throughput_FWHM=self.Throughput_FWHM, QElambda=self.QElambda, atmlambda=self.atmlambda, fraction_lya=self.fraction_lya,sky_lines=self.sky_lines,Redshift=self.Redshift, IFS=self.IFS, spectrograph=self.spectrograph,interpolation=self.interpolation,source_im=self.source_im,minmax=self.minmax, test=self.test,dependencies=self.dependencies)
        
        wids2 = widgets.interactive(self.update_instrument,instrument=self.instrument)
        # Question: why do we need to reset it to True here?
        self.change.value=True

        def reset(_):
            self.update_instrument(self.instrument.value)
            self.Signal.value = self.instruments_dict[self.instrument.value]["Signal"]
        self.reset.on_click(reset)

        def save_plot(_):
            self.fig.savefig("/tmp/fig1.png", dpi=100, bbox_inches="tight")
            self.fig2.savefig("/tmp/fig2.png", dpi=100, bbox_inches="tight")
            self.fig3.savefig("/tmp/fig3.png", dpi=100, bbox_inches="tight")
        self.save_plot_button.on_click(save_plot)

        def save_data(_):
            self.f = lambda x: self.wavelength.value + (self.dispersion.value/10) * (x - n1/2)
            fitswrite(self.ifs_cube.T,"/tmp/ifs_cube.fits")
            fitswrite(np.transpose(self.ifs_cube_stack,(1,2,0)),"/tmp/ifs_cube_stack.fits")
            # fitswrite(self.imaADU_without_source,"/tmp/imaADU_without_source.fits")
            # fitswrite(self.imaADU_source,"/tmp/imaADU_source.fits")
            # fitswrite(self.imaADU_stack_without_source,"/tmp/imaADU_stack_without_source.fits")
            # fitswrite(self.imaADU_stack_only_source,"/tmp/imaADU_stack_only_source.fits")
            np.savetxt("/tmp/spectra.csv", np.asarray([ self.f(np.arange(n1)), self.ifs_spectra[0].get_ydata(), self.ifs_spectra_stack[0].get_ydata(), self.ifs_spectra_background[0].get_ydata(), self.ifs_spectra_background_stack[0].get_ydata()  ]), delimiter=",",header="wavelength,ifs_spectra, ifs_spectra_stack, ifs_spectra_background, ifs_spectra_background_stack")
            return
        self.save_data_button.on_click(save_data)


        self.obs_tab = VBox([HBox([self.Signal,self.Size_source ]), HBox([self.Sky,self.Line_width])]) 
        self.strat_tab = VBox([HBox([self.Atmosphere,self.acquisition_time ]), HBox([self.exposure,self.wavelength])  , HBox([self.lambda_stack])     ])
        self.inst_tab = VBox([HBox([self.Collecting_area,self.pixel_scale, self.spectrograph ]), HBox([self.Throughput,self.fwhm ])])
        self.spectro_tab = VBox([HBox([self.Spectral_resolution,self.SlitDims ]), HBox([self.dispersion,self.IFS ])]) 
        self.det_tab = VBox([HBox([self.QE,self.RN ]), HBox([self.Dark_current,self.extra_background]),          HBox([self.cosmic_ray_loss_per_sec,self.EM_gain ]), HBox([self.CIC_charge,self.counting_mode])])#,self.pixel_size
        # self.amp_tab = HBox([self.EM_gain,self.CIC_charge])
        self.fb_tab = VBox([HBox([self.follow_temp ]),HBox([self.temperature,self.smearing]),HBox([self.change]) ])
        self.im_tab =         VBox([HBox([self.spectra,self.Redshift,self.Throughput_FWHM ]), HBox([self.units,self.Δx, self.Δλ])   , HBox([self.interpolation,self.QElambda,self.atmlambda,self.sky_lines,self.source_im]) , HBox([self.minmax])   ]) #self.fraction_lya
        
        # if ("FIREBall-2" in instrument):
        #     tab_contents = [ "Observed Source", "Observation strategy" , "Instrument Design","Detector performance",  "Spectrograph Design","FIREBall specific"]#, ''] #,"Amplified Detector"
        #     children = [ self.obs_tab, self.strat_tab,  self.inst_tab,  self.det_tab,self.spectro_tab, self.fb_tab]#, self.]
        # else:
        tab_contents = [ "Observed Source", "Observation strategy" , "Instrument Design","Detector performance",  "Spectrograph Design"]#,"Imaging"] #,"Amplified Detector"
        children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab, self.spectro_tab]#, self.im_tab] 




        self.controls = widgets.Tab()# Accordion
        self.controls.children = children
        for i, name in enumerate(tab_contents):
            self.controls.set_title(i, name)

        self.out1 = widgets.Output()
        self.out2 = widgets.Output()
        self.out3 = widgets.Output()
        # self.output_tabs = widgets.Tab(children = [self.out2, self.out1]);self.output_tabs.set_title(0, 'Image');self.output_tabs.set_title(1, 'SNR')
        if self.IFS.value is False:
            self.output_tabs = widgets.Tab(children = [self.out1, self.out2]); self.output_tabs.set_title(0, 'SNR');self.output_tabs.set_title(1, 'Spectral image')
            #;self.output_tabs.set_title(2, 'IFS Image')
            self.plot_shown = False#  # Reset the plot flag when hiding tab 3    # self.output_tabs.selected_index = 0  # Sélectionner un autre onglet par défaut
            self.first_plot = True
        else:    
            self.output_tabs = widgets.Tab(children = [self.out1, self.out2, self.out3]); self.output_tabs.set_title(0, 'SNR');self.output_tabs.set_title(1, 'Spectral image');self.output_tabs.set_title(2, 'IFS Image')
            self.plot_shown = True#  # Reset the plot flag when hiding tab 3    # self.output_tabs.selected_index = 0  # Sélectionner un autre onglet par défaut
            self.first_plot = False

        def on_tab_change(change): 
            self.on_instrument_change()



        self.output_tabs.observe(on_tab_change, names='selected_index')
        
        def database_change(change):
            sheet_id = "1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU"
            sheet_name = "instruments.csv"
            url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
            with self.out1:
                inst = self.instrument.value
                if self.database.value == "Online DB":
                    self.instruments = Table.from_pandas(pd.read_csv(url))
                elif self.database.value == "Local DB":
                    # instruments = Table.read("../data/Instruments/Instruments.csv")
                    # instruments = Table.read("../Instruments.csv")
                    instruments = Table.from_pandas(pd.read_excel("../instruments.xlsx"))

                    self.instruments = instruments[instruments.colnames]
                    for col_name in self.instruments.colnames[3:]:
                        self.instruments[col_name] = to_float_or_nan(self.instruments[col_name])
                self.instruments_dict = {name: {key: val for key, val in zip(self.instruments["Charact."][:], self.instruments[name][:]) if not isinstance(key, np.ma.core.MaskedConstant) and not isinstance(val, np.ma.core.MaskedConstant)} for name in self.instruments.colnames[3:]}
                self.instrument.options = self.instruments.colnames[3:]
                if inst in self.instruments.colnames[3:]:
                    # self.instrument.value = inst
                    self.update_instrument(instrument)
                else:
                    self.instrument.value = self.instruments.options[0]
            return
            
        self.database.observe(database_change, names='value')
        
        
        new = VBox([ HBox([self.instrument, self.database, self.x_axis,self.SNR_res,self.xlog, self.ylog,self.reset,self.save_plot_button,self.save_data_button ,self.test,self.dependencies  ])  ,self.controls , self.output_tabs ]) #


        with self.out1: # before because then the arrays are transformed into numbers
            # print(self.IFS.value)
            self.fig = self.new.PlotNoise(x=x_axis)
            args, _, _, locals_ = inspect.getargvalues(inspect.currentframe())
            self.v=[]
            for j, ax in enumerate(self.fig.axes):
                if j==2:
                    label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f [%0.2f]'%(x_axis,float_to_latex(self.time),float_to_latex(exposure_time[np.nanargmax(self.new.SNR)]),self.new.SNR[self.arg],np.nanmax(self.new.SNR))#, self.new.gain_thresholding[arg])
                    self.v.append(ax.axvline(self.time,ls=':',c='k',label=label))
                    ax.legend(loc='upper right')
                else:
                    self.v.append(ax.axvline(self.time,ls=':',c='k'))
            self.ax0 =  self.fig.axes[0]
            self.ax1 =  self.fig.axes[1]
            self.ax2 =  self.fig.axes[2]
            self.ax3 =  self.fig.axes[3]
            self.ax3.set_ylim((-19,-12.5))
            self.ax0.set_xscale('log')
            self.fig.canvas.toolbar_position = 'bottom'
            title = 'Instrument=%s, FOV=%samin$^2$, λ=%inm, Throughput=%i%%, Atm=%i%%, Platescale=%.1f, area=%0.1fm$^2$'%(instrument, self.instruments[instrument][self.instruments["Charact."]=="FOV_size"][0], self.instruments[instrument][self.instruments["Charact."]=="wavelength"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Throughput"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Atmosphere"][0], self.instruments[instrument][self.instruments["Charact."]=="pixel_scale"][0], self.instruments[instrument][self.instruments["Charact."]=="Collecting_area"][0])
            self.ax0.set_title(title,y=0.97,fontsize=10)

            plt.show(self.fig)  
        self.im,self.im_stack = self.new.SimulateFIREBallemCCDImage( Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Baseline Spectra",size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=self.QElambda.value,atmlambda=self.QElambda.value,fraction_lya=self.fraction_lya.value, Full_well=self.Full_well, Altitude=self.Altitude, conversion_gain=self.conversion_gain, Throughput_FWHM=self.Throughput_FWHM.value,sky_lines=self.sky_lines.value, Redshift=self.Redshift.value,source_image=self.source_im.value)
        # self.im,self.im_stack, self.cube_stack, self.im0, source_im_wo_atm, self.imaADU_stack_only_source, self.imaADU_without_source, self.imaADU_stack_without_source, self.imaADU_source = self.new.SimulateFIREBallemCCDImage( Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Baseline Spectra",size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=self.QElambda.value,atmlambda=self.QElambda.value,fraction_lya=self.fraction_lya.value, Full_well=self.Full_well, Altitude=self.Altitude, conversion_gain=self.conversion_gain, Throughput_FWHM=self.Throughput_FWHM.value,sky_lines=self.sky_lines.value, Redshift=self.Redshift.value)

        center = n1/2
        f = lambda x: self.wavelength.value + self.dispersion.value * (x - center)
        g = lambda x: (x - self.wavelength.value) / self.dispersion.value + center
     
        with self.out2:
            self.current_cmap = cmocean.cm.deep# current_cmap = cmocean.cm.solar# self.current_cmap = cmocean.cm.thermal
            

            self.current_cmap.set_bad(color='black')
            self.bins=np.arange(-100,4000,100)
            self.bins=np.linspace(-100,np.nanmax(self.im),100)

            self.mod = mostFrequent(self.im_stack[:20,:].flatten())

            self.limit = self.mod+self.new.n_threshold * self.RN.value

            self.fig2 = plt.figure(figsize=(12, 8),)
            gs = self.fig2.add_gridspec(2, 2, height_ratios=[0.5, 1])
            self.nax = self.fig2.add_subplot(gs[0, 0])      # Top-left
            self.nax0 = self.fig2.add_subplot(gs[0, 1], sharex = self.nax, sharey = self.nax)     # Top-right
            gs_nax2 = gs[1, 0].subgridspec(2, 1, height_ratios=[1, 1])
            self.nax2_1 = self.fig2.add_subplot(gs_nax2[0, 0], sharex = self.nax)  # Bottom-left, top half
            self.nax2 = self.fig2.add_subplot(gs_nax2[1, 0])    # Bottom-left, bottom half
            self.nax1 = self.fig2.add_subplot(gs[1, 1], sharex = self.nax) 
            
            self.nax1_secondary = self.nax1.secondary_xaxis("top", functions=(f,g))
            self.nax2_1_secondary = self.nax2_1.secondary_xaxis("top", functions=(f,g))
            self.nax2_1.set_title('Wavelength (nm)',fontsize=10)
            self.nax2_1.set_yscale(self.yscale)

            # if self.counting_mode.value:
            #     stacked_image = np.nansum(self.cube_stack>self.limit,axis=0)
            #     im0 = self.nax0.imshow(stacked_image, aspect="auto",cmap=self.current_cmap)
            # else:
            im0 = self.nax0.imshow(self.im_stack, aspect="auto",cmap=self.current_cmap)#,interpolation=interpolation)
            im = self.nax.imshow(self.im, aspect="auto",cmap=self.current_cmap)#,interpolation=interpolation)
            labels =  ['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]
            self.nax.plot(0,0,".",label="\n".join(labels))
            self.nax.legend(loc="upper left",handlelength=0, handletextpad=0, fancybox=True,markerscale=0,fontsize=6)
            # self.nax.text(0,0,"\n".join(labels))
            im2=self.nax0.imshow(self.im_stack, aspect="auto",cmap=self.current_cmap)#,interpolation=interpolation)
            self.nax0.get_xaxis().set_ticks([])
            self.nax0.get_yaxis().set_ticks([])
            self.nax.get_xaxis().set_ticks([])
            self.nax.get_yaxis().set_ticks([])
            self.nax.set_title('Single image: FOV = %i" × %iÅ, λ~%iÅ'%(100*self.pixel_scale.value,500*self.dispersion.value,10*self.wavelength.value))
            self.nax0.set_title('Stacked image: Pixel size = %0.2f" × %0.2fÅ'%(self.pixel_scale.value,self.dispersion.value))
            self.nax1.set_title('Wavelength (nm)',fontsize=10)


            self.nax2.set_xlabel("Pixels' values [Unit]")
            self.nax1.set_xlabel('Pixels')
            # FIXME besure to only use correctly the factor 2.35
            self.l1 = self.nax1.plot(self.im[:,int(n1/2-self.Line_width.value/self.pixel_scale.value):int(n1/2+self.Line_width.value/self.pixel_scale.value)].mean(axis=1),ls='-',lw=3,label='Single exp profiles',alpha=0.2)
            self.l2 = self.nax1.plot(np.linspace(0,n1,n1),self.im[int(n2/2 - self.Size_source.value/self.pixel_scale.value/2/2.35):int(n2/2 + self.Size_source.value/self.pixel_scale.value/2/2.35),:].mean(axis=0),ls='-',lw=3,alpha=0.2,c="k")



            # self.final_sky_before_convolution =  self.nax2_1.plot(self.new.final_sky_before_convolution*np.ones(n1),ls='-',c="k",label="final_sky_before_convolution")#,lw=3,alpha=0.2)
            # self.atm_trans_before_convolution =  self.nax2_1.plot(self.new.atm_trans_before_convolution*np.ones(n1),ls=':',c="k",label="atm_trans_before_convolution")#,lw=3,alpha=0.2)

            self.absorption         =  self.nax2_1.plot(self.new.atm_trans*np.ones(n1),ls=':',c="k",label=r"Atm$_{Trans}$(λ) for Signal")#,lw=3,alpha=0.2)
            self.emission_lines     =  self.nax2_1.plot(self.new.final_sky*np.ones(n1),ls='-',c="k",label=r"Sky(λ) Em lines")#,lw=3,alpha=0.2)
            self.Throughput_curve     =  self.nax2_1.plot(self.new.Throughput_curve*np.ones(n1),ls='-.',c="k",label="[QExThroughtput](λ)")#,lw=3,alpha=0.2)
            self.nax2_1.legend(loc='lower right',fontsize=8,title="Transmission curves")
            # self.nax2_1.set_ylim((-0.05,1.05))
            self.nax2_1.set_ylim((1e-2,1.05))


            # self.nax1bis = self.nax1.twinx()
            self.hw, self.hl =  self.Slitwidth/2/self.pixel_scale.value ,  self.Slitlength/2/self.pixel_scale.value
            self.nax.plot([n1/2 - self.hw,n1/2 + self.hw,n1/2 + self.hw,n1/2 - self.hw,n1/2 - self.hw],[n2/2 - self.hl,n2/2 - self.hl,n2/2 + self.hl,n2/2 + self.hl,n2/2 - self.hl],"-k")
            self.nax0.plot([n1/2 - self.hw,n1/2 + self.hw,n1/2 + self.hw,n1/2 - self.hw,n1/2 - self.hw],[n2/2 - self.hl,n2/2 - self.hl,n2/2 + self.hl,n2/2 + self.hl,n2/2 - self.hl],"-k", label="Slit=%0.1f'' × %0.1f''"%(self.Slitwidth,self.Slitlength))
            self.slit_text = AnchoredText("Slit=%0.1f'' × %0.1f''"%(self.Slitwidth,self.Slitlength), frameon=True, loc=4, pad=0.5,prop={'fontsize': 8})
            plt.setp(self.slit_text.patch, facecolor='white', alpha=0.5)
            self.nax.add_artist(self.slit_text)
            
            
            self.nax0.legend(loc='upper right',fontsize=8)

            self.profile = np.mean(im0.get_array().data[:,int(n1/2-self.Line_width.value/self.pixel_scale.value):int(n1/2+self.Line_width.value/self.pixel_scale.value)],axis=1)
            spatial_profile = self.im[:,:].mean(axis=1)
            # what means this convolution
            self.nax1.lines[0].set_ydata(spatial_profile)#np.convolve(spatial_profile,3,mode="same"))
            self.profile = np.mean(im0.get_array().data[:,:],axis=1)
            self.l1_s = self.nax1.plot(self.profile,label='Stack. spatial prof',c=self.l1[0].get_color())
            self.l2_s = self.nax1.plot(self.im_stack[int(n2/2 - self.Size_source.value/self.pixel_scale.value/2/2.35):int(n2/2 + self.Size_source.value/self.pixel_scale.value/2/2.35),:].mean(axis=0),label='Stack. spectral prof',c=self.l2[0].get_color())
            
            self.nax1.get_xaxis().set_ticks_position('bottom')

            self.nax1.tick_params(axis='x', labelbottom=True)

            
            try:
                self.popt, self.pcov = curve_fit(gaus,np.arange(len(self.profile)),self.profile,p0=[np.ptp(self.profile), 50, 5, self.profile.min()])
            except RuntimeError:
                self.popt = [0,0,0,0]
            # self.fit = PlotFit1D(x= np.arange(len(self.profile)),y=self.profile,deg="gaus", plot_=False,ax=self.nax1bis,c="k",ls=":",P0=[np.ptp(self.profile), 50, 5, self.profile.min()])
            # self.nax1bis.plot( np.arange(len(self.profile)),gaus( np.arange(len(self.profile)),*self.popt),":k",label="SNR=%0.1f/%0.1f=%0.1f"%(self.popt[0]**2,self.profile[:20].std(),self.popt[0]**2/self.profile[:20].std()))
            self.l3_s = self.nax1.plot( np.arange(len(self.profile)),gaus( np.arange(len(self.profile)),*self.popt),":k",label="SNR=%0.1f/%0.1f=%0.1f"%(self.popt[0]**2,self.profile[:20].std(),self.popt[0]**2/self.profile[:20].std()))

            self.nax1.set_xlim((0,n1))
            self.nax1.legend(loc="upper right",fontsize=8,title="Averaged profiles")
            # self.nax1bis.legend(loc='upper right',fontsize=8,title="Averaged profiles")
            _,_,self.bars1 = self.nax2.hist(self.im.flatten(),bins=self.bins,alpha=0.3,color=self.l1[0].get_color(),label='Single image')
            _,_,self.bars2 = self.nax2.hist(self.im_stack.flatten(),bins=self.bins,alpha=0.3,color=self.l2[0].get_color(),label='Averaged stack')
            # TODO change the 40 in the next formula
            title = 'Signal kept=%i%%, RN kept=%i%%, Signal/tot=%i%%'%(100*self.new.Photon_fraction_kept[0], 100*self.new.RN_fraction_kept[0],100*(np.mean(self.im_stack[40:-40,:])-np.mean(self.im_stack[:20,:]))/np.mean(self.im_stack[40:-40,:]))
            self.nax2.plot([self.mod,self.mod],[0,100],c="k",ls=":",label="Bias %0.3f, PC limit %0.3f (%s):\n%s"%(self.mod,self.limit[i], self.new.counting_mode, title))
            self.nax2.plot([self.limit,self.limit],[0,100],c="k",ls=":")#,label="PC limit %i: %s"%(self.limit, title))
                        
            self.nax2.legend(loc='upper right',fontsize=8,title="Histogram - Pixels'values")
            self.nax2.set_xlim(xmin=-10, xmax=np.nanmax(self.bins))
            self.cax = make_axes_locatable(self.nax).append_axes('bottom', size='15%', pad=0.05)
            self.cax0 = make_axes_locatable(self.nax0).append_axes('bottom', size='15%', pad=0.05)
            self.cbar1 = self.fig2.colorbar(im, cax=self.cax, orientation='horizontal')
            self.cbar2 = self.fig2.colorbar(im2, cax=self.cax0, orientation='horizontal')
            self.cbar1.formatter.set_powerlimits((0, 0))
            self.cbar2.formatter.set_powerlimits((0, 0))
            self.cbar1.formatter.set_useMathText(True)
            self.cbar2.formatter.set_useMathText(True)


            self.fig2.canvas.toolbar_position = 'bottom'
            self.fig2.tight_layout()
            plt.show(self.fig2)


        with self.out3:
            self.out3.clear_output(wait=True)
            self.current_cmap = cmocean.cm.deep# current_cmap = cmocean.cm.solar# self.current_cmap = cmocean.cm.thermal
            self.current_cmap.set_bad(color='black')
            n3 = int(np.sqrt(60*60*self.FOV_size)/self.Slitwidth)

            self.fig3 = plt.figure(figsize=(12, 8))
            gs = self.fig3.add_gridspec(2,2,height_ratios=[2,0.5])
            self.nax20 = self.fig3.add_subplot(gs[0,0])
            self.nax21 = self.fig3.add_subplot(gs[0,1], sharex = self.nax20, sharey = self.nax20)
            self.only_spectra=False
            if self.only_spectra:
                self.nax2s = self.fig3.add_subplot(gs[1,:]) # HACK vincent subplot change this
            else:
                self.nax2s = self.fig3.add_subplot(gs[1,0]) # HACK vincent subplot change this
                self.nax2_rp = self.fig3.add_subplot(gs[1,1]) # HACK vincent subplot change this
                self.nax2_rp2 = self.nax2_rp.twinx()

            self.nax2s_secondary = self.nax2s.secondary_xaxis("top", functions=(f,g))





            self.nax20.set_xlabel("Spatial pixel")
            self.nax20.set_ylabel("Spatial pixel")
            self.nax21.set_xlabel("Spatial pixel")
            self.nax21.set_ylabel("Spatial pixel")
            self.nax20.set_title("Single cube: %ipix × %ipix"%(n2,n3),fontsize=8)
            self.nax21.set_title("Stacked Cube: %0.1f' × %0.1f', full FOV= %0.1f' × %0.1f'"%(n2/self.pixel_scale.value/60,n3/self.pixel_scale.value/60,np.sqrt(self.FOV_size),np.sqrt(self.FOV_size)),fontsize=8)


            self.nax2s.set_title("Wavelength (nm)",fontsize=8)
            self.nax2s.set_xlabel("Spectral pixel")
            self.nax2s.set_ylabel("Flux")
            if ~self.only_spectra:
                self.nax2_rp.set_title("Radial profile",fontsize=8)
                self.nax2_rp.set_xlabel("Spatial pixel")
                self.nax2_rp.set_ylabel("Flux")
                self.rp2 = self.nax2_rp.plot([np.nan],[np.nan],"k-", alpha=0.2,lw=3,label="Single")
                self.rp = self.nax2_rp.plot([np.nan],[np.nan],"k-",label="Stack")
                self.rp_fit = self.nax2_rp.plot([np.nan],[np.nan],"k:",label="Fit")
                self.res1 =  self.nax2_rp.axvline(x=np.nan, color='red', linestyle='--',label="Half resolution elements",lw=0.5)
                self.res2 =  self.nax2_rp.axvline(x=np.nan, color='red', linestyle='--',lw=0.5)
                self.nax2_rp.legend(fontsize=7,loc="upper right",title="Radial profiles",title_fontsize=8)
                self.diff = self.nax2_rp2.plot([np.nan],[np.nan],"k-",label="diff",alpha=0.2)


            self.ifs_cube = np.zeros((n2,n1,n3))

            self.wavelength_line1 = self.nax2s.axvline(int(n1/2)+0.5,ls='--',c='k',lw=0.7)
            self.wavelength_line2 = self.nax2s.axvline(int(n1/2)+0.5,ls='--',c='k',lw=0.7)
            if IFS:

                self.ifs_spectra = self.nax2s.plot(self.im_stack[int(n2/2),:],"k-",lw=3,alpha=0.1,label="Spectra: On source")
                self.ifs_spectra_stack = self.nax2s.plot(self.im[int(n2/2),:],"k--",alpha=0.5,label="Stacked spectra: On source")
            
                self.ifs_integ_spectra_stack = self.nax2s.plot(       np.nanmean(self.im[int(n2/2-self.Slitwidth/self.pixel_scale.value):int(n2/2+self.Slitwidth/self.pixel_scale.value),:],axis=0)      ,"k-",label="Integrated stacked spectra: On source")
                x1,x2 = n3/2-(n2/n2)*self.Size_source.value/self.pixel_scale.value,n3/2+(n2/n2)*self.Size_source.value/self.pixel_scale.value
                y1, y2 = n2/2-self.Size_source.value/self.pixel_scale.value,n2/2+self.Size_source.value/self.pixel_scale.value
                self.stack_square = self.nax21.plot([x1,x2,x2,x1,x1],[y2,y2,y1,y1,y2],"k:")
                
                # self.ifs_spectra_background  = self.nax2s.plot(self.imaADU_without_source[int(n2/2),:],"k-",alpha=0)
                # self.ifs_spectra_background_stack = self.nax2s.plot((self.im_stack-self.imaADU_stack_only_source)[int(n2/2),:],"k:",label="Spectra: Field edge (background)")


                gaussian = norm.pdf(np.arange(-n3,n3,2), loc=0, scale=self.Size_source.value/2)
                ratios = (gaussian - gaussian.min()) / (gaussian.max() - gaussian.min())
                ratios_reshaped = ratios[np.newaxis, np.newaxis, :]
                indices = np.array([np.random.permutation(n3) for _ in range(n1 * n2)])
                indices = indices.reshape(n2, n1, n3)                                

                self.ifs_cube = np.repeat(self.im[:, :, np.newaxis], n3, axis=2)
                self.ifs_cube = np.take_along_axis(self.ifs_cube, indices, axis=1)
                # self.ifs_cube += self.im[:, :, np.newaxis] * ratios_reshaped
                self.ifs_cube_stack = np.ones((n2,n1,n3))#np.repeat( self.im_stack[:, :, np.newaxis], n3, axis=2)
                # self.ifs_cube_stack = np.take_along_axis(self.im_stack, indices, axis=1)
                # self.ifs_cube_stack += self.im_source[:, :, np.newaxis] * ratios_reshaped

                # self.ifs_cube = np.repeat(self.imaADU_without_source[:, :, np.newaxis], n3, axis=2)
                # self.ifs_cube = np.take_along_axis(self.ifs_cube, indices, axis=1)
                # self.ifs_cube += self.imaADU_source[:, :, np.newaxis] * ratios_reshaped
                # self.ifs_cube_stack = np.repeat( self.imaADU_source[:, :, np.newaxis], n3, axis=2)
                # self.ifs_cube_stack = np.take_along_axis(self.ifs_cube_stack, indices, axis=1)
                # self.ifs_cube_stack += self.imaADU_source[:, :, np.newaxis] * ratios_reshaped
            else:
                self.stack_square = self.nax21.plot([np.nan],[np.nan],"k:")
                self.ifs_slice = self.nax20.imshow(np.nan*self.ifs_cube[:,0,:], aspect="auto",cmap=self.current_cmap)#,interpolation=interpolation)
                self.ifs_slice_stack = self.nax21.imshow(np.nan*self.ifs_cube[:,0,:], aspect="auto",cmap=self.current_cmap)#,interpolation=interpolation)

                self.ifs_spectra = self.nax2s.plot([0,n1],[np.nan,np.nan],"k-",lw=3,alpha=0.1,label="Spectra: On source")
                self.ifs_spectra_stack = self.nax2s.plot([0,n1],[np.nan,np.nan],"k--",label="Stacked spectra: On source")
                self.ifs_integ_spectra_stack = self.nax2s.plot([0,n1],[np.nan,np.nan],"k-",label="Integrated stacked spectra: On source")
                self.ifs_spectra_background       = self.nax2s.plot([0,n1],[np.nan,np.nan],"k-",alpha=0)
                self.ifs_spectra_background_stack = self.nax2s.plot([0,n1],[np.nan,np.nan],"k:",label="Spectra: Field edge (background)")
                # self.ifs_slice = 0  # actually you need to initiate these variable in any case
            self.nax2s.legend(fontsize=7,loc="upper right")
            self.cax_slicer = make_axes_locatable(self.nax20).append_axes('right', size='5%', pad=0.05)
            self.cax_slicer0 = make_axes_locatable(self.nax21).append_axes('right', size='5%', pad=0.05)
            self.cbar_slicer1 = self.fig2.colorbar(im, cax=self.cax_slicer, orientation='vertical')
            self.cbar_slicer2 = self.fig2.colorbar(im2, cax=self.cax_slicer0, orientation='vertical')
            self.cbar_slicer1.formatter.set_powerlimits((0, 0))
            self.cbar_slicer2.formatter.set_powerlimits((0, 0))
            self.cbar_slicer2.formatter.set_useMathText(True)
            self.cbar_slicer1.formatter.set_useMathText(True)


            self.position1 = self.nax20.plot(int(n3/2),int(n2/2),"ro")
            self.position2 = self.nax21.plot(int(n3/2),int(n2/2),"ro")
            self.fig3.tight_layout()
            self.fig3.canvas.toolbar_position = 'bottom'
            plt.show(self.fig3)
            # if self.IFS.value:
            #     plt.show(self.fig3)
            # else:
            #     # self.fig3.clf()
            #     plt.close(self.fig3)
        display(HBox([self.output,new]))



    def hide_tab3(self):
        self.out3.layout.display = 'none'
        self.output_tabs.children = [self.out1, self.out2]  # Réinitialiser les onglets sans out3
        self.plot_shown = False  # Reset the plot flag when hiding tab 3    # self.output_tabs.selected_index = 0  # Sélectionner un autre onglet par défaut

    def show_tab3(self):
        self.out3.layout.display = 'block'
        self.output_tabs.children = [self.out1, self.out2, self.out3]  # Ajouter out3 à nouveau
        self.output_tabs.set_title(2, 'IFS Image')
        if not self.plot_shown:
            with self.out3:
                if self.first_plot:
                    self.fig3.canvas.draw()
                    self.first_plot = False
                else:
                    plt.show(self.fig3)
                self.plot_shown = True



    def on_instrument_change(self):
        # with self.out1:
        #     print(self.instruments_dict[self.instrument.value]["dispersion"])
        self.spectro = False if np.isnan(self.instruments_dict[self.instrument.value]["dispersion"]) else True
        self.spectrograph.value = False if np.isnan(self.instruments_dict[self.instrument.value]["dispersion"]) else True
        self.Δλ.layout.visibility = 'hidden' if self.output_tabs.selected_index==1 else 'visible'  

        self.change.value=True
        selected_index = self.output_tabs.selected_index
        if selected_index == 0:
            self.x_axis.layout.visibility = 'visible'  
            self.SNR_res.layout.visibility = 'visible'  
            self.xlog.layout.visibility = 'visible'  
            # self.ylog.layout.visibility = 'visible'  
            with self.out1:
                self.Size_source.value-=0.01
            # if ("FIREBall-2" in self.instrument.value):
            #     self.spectro_tab.layout.display = 'block'
            #     self.fb_tab.layout.display = 'block'
            #     self.controls.children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab, self.spectro_tab, self.fb_tab]
            #     self.controls.set_title(5, 'FIREBall specific')
            # else:
            if self.spectro:
                self.spectro_tab.layout.display = 'block'
                self.fb_tab.layout.display = 'none'
                self.controls.children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab, self.spectro_tab]                    
                self.controls.set_title(4, 'Spectrograph design')
            else:
                self.fb_tab.layout.display = 'none'
                self.spectro_tab.layout.display = 'none'
                self.controls.children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab]                    
        else:
            self.x_axis.layout.visibility = 'hidden'  
            self.SNR_res.layout.visibility = 'hidden'  
            self.xlog.layout.visibility = 'hidden'  
            # self.ylog.layout.visibility = 'hidden'  

            # if ("FIREBall-2" in self.instrument.value):
            #     self.spectro_tab.layout.display = 'block'
            #     self.fb_tab.layout.display = 'block'
            #     self.controls.children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab, self.spectro_tab, self.im_tab, self.fb_tab]
            #     self.controls.set_title(4, 'Spectrograph design')
            #     self.controls.set_title(5, 'Imaging')
            #     self.controls.set_title(6, 'FIREBall specific')
            # else:
            if self.spectro:
                self.spectro_tab.layout.display = 'block'
                self.fb_tab.layout.display = 'none'
                self.controls.children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab, self.spectro_tab, self.im_tab]                    
                self.controls.set_title(4, 'Spectrograph design')
                self.controls.set_title(5, 'Imaging')
            else:
                self.fb_tab.layout.display = 'none'
                self.spectro_tab.layout.display = 'none'
                self.controls.children = [ self.obs_tab, self.strat_tab,  self.inst_tab, self.det_tab, self.im_tab]                    
                self.controls.set_title(4, 'Imaging')

            if selected_index == 1:
                with self.out2:
                    self.Size_source.value-=0.01
            elif selected_index == 2:
                with self.out3:
                    self.Size_source.value+=0.01
        
        return
    


    def update(self,x_axis, counting_mode,Sky,acquisition_time,Signal,EM_gain,RN,CIC_charge,Dark_current,exposure,smearing,temperature,follow_temp,fwhm,QE,extra_background, log,xlog,SNR_res,
    Collecting_area, pixel_scale, Throughput, Spectral_resolution, SlitDims, dispersion,
    Size_source,Line_width,wavelength,Δλ,Δx, Atmosphere, pixel_size,cosmic_ray_loss_per_sec,lambda_stack, #change, 
    spectra,units,Throughput_FWHM, QElambda, atmlambda, fraction_lya,sky_lines, Redshift, IFS, spectrograph,interpolation, source_im,minmax ,test, dependencies
    ):
        """
        Update values in the ETC plot
        """
        self.show_tab3() if self.IFS.value else self.hide_tab3()

        with self.out1:



            # print(counting_mode,Sky,acquisition_time,Signal,EM_gain,RN,CIC_charge,Dark_current,exposure,smearing,temperature,follow_temp,fwhm,QE,extra_background, log,xlog,SNR_res,    Collecting_area, pixel_scale, Throughput, Spectral_resolution, SlitDims, dispersion,    Size_source,Line_width,wavelength,Δλ,Δx, Atmosphere, pixel_size,cosmic_ray_loss_per_sec,lambda_stack, change,     spectra,units,Throughput_FWHM, QElambda, atmlambda, fraction_lya,sky_lines, Redshift, IFS)
            # self.yscale="symlog" if log else "linear"
            self.yscale="log" if log else "linear"

            if self.change.value:
                PSF_RMS_mask=fwhm[0]
                PSF_RMS_det=fwhm[1]
                # self.Slitwidth=SlitDims[0]
                # self.Slitlength=SlitDims[1]
                Slitwidth=SlitDims[0]
                Slitlength=SlitDims[1]
                self.Slitwidth=SlitDims[0]
                self.Slitlength=SlitDims[1]
                readout_time=exposure[0]
                exposure_time=exposure[1]
                if follow_temp:
                    self.follow_temp.value=follow_temp
                    self.Dark_current.value = self.dark_poly(temperature)
                    self.smearing.value = self.smearing_poly(temperature)
                    # self.CIC_charge.value = self.CIC_poly(self.temperature.value)
                if 1==1: # HACK here is just to be able to plot x_axis for fireball. Issue is that when change x_axis options, it comes back to the first choice
                    if "FIREBall" in self.instrument.value:
                        options = [x_axis] + self.fb_options if self.follow_temp.value else [x_axis] + self.fb_options_no_temp
                    else:
                        options = [x_axis] +  self.other_options if self.spectrograph.value else [x_axis] +  self.other_options_imager
                    self.x_axis.options = options
                    if x_axis in options:
                        self.x_axis.value=x_axis
                    if  "OBSERVED SOURCE" in self.x_axis.value:
                        self.x_axis.value, x_axis = "Signal", "Signal"
                    elif  "OBSERVATION STRATEGY" in self.x_axis.value:
                        self.x_axis.value, x_axis = "Atmosphere","Atmosphere"
                    elif  "INSTRUMENT DESIGN" in self.x_axis.value:
                        self.x_axis.value, x_axis = "Collecting_area","Collecting_area"
                    elif  "SPECTROGRAPH DESIGN" in self.x_axis.value:
                        self.x_axis.value, x_axis = "Spectral_resolution","Spectral_resolution"
                    elif  "DECTECTOR PERFORMANCE" in self.x_axis.value:
                        self.x_axis.value, x_axis = "exposure_time", "exposure_time"
                    elif  "AMPLIFIED" in self.x_axis.value:
                        self.x_axis.value, x_axis = "EM_gain","EM_gain"



                args, _, _, locals_ = inspect.getargvalues(inspect.currentframe())
                if x_axis in locals_:
                    value = locals_[x_axis]
                else:
                    value = getattr(self,x_axis)
                    if (type(value) != float) & (type(value) != int):
                        value = self.instruments_dict[self.instrument.value][x_axis]
                    # error here! getattr(self,x_axis) can be 

                names = ["Signal","Dark current","Sky", "CIC", "Read noise","Extra Background"]
 

                self.smearing.layout.visibility = 'visible' if ("FIREBall-2" in self.instrument.value) & (self.counting_mode.value)    else 'hidden'
                self.temperature.layout.visibility = 'visible' if ("FIREBall-2" in self.instrument.value) &  (self.follow_temp.value)  else 'hidden'
                if self.spectrograph.value:
                    self.Throughput_FWHM.layout.visibility = 'visible' if ((self.QElambda.value) &  ~(os.path.exists("../data/Instruments/%s/Throughput.csv"%(self.instrument.value.upper().replace(" ","_"))) )      )  else 'hidden'
                else:
                    self.Throughput_FWHM.layout.visibility = 'visible'

                if x_axis == 'temperature':
                    temperature=np.linspace(self.temperature.min, self.temperature.max)
                    Dark_current = 10**self.dark_poly(temperature)
                    # smearing = np.poly1d([-0.0306087, -2.2226087])(temperature)
                    smearing = self.smearing_poly(temperature)
                # d = {name:np.linspace(rgetattr(self, '%s.min'%(name)), rgetattr(self, '%s.max'%(name))   name for self.fb_options_no_temp}

                self.len_xaxis = 50
                def space(a, b):
                    if (self.xlog.value) & (a>=0):
                        if a==0:
                            y = np.logspace(np.log10(np.max([a,0.0001])),np.log10(b),self.len_xaxis) 
                        else:
                            y = np.logspace(np.log10(a),np.log10(b),self.len_xaxis) 
    
                    else:
                        y = np.linspace(a,b,self.len_xaxis)
                    return y

                if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out1):  

                    if x_axis == 'exposure_time':
                        exposure_time=space(1,self.time_max)
                    elif x_axis == 'Sky':
                        # Sky=np.logspace(-19,-15)
                        Sky = space(10**self.Sky.min,10**self.Sky.max)
                    elif x_axis == 'Signal':
                        Signal=space(10**self.Signal.min,10**self.Signal.max)
                    elif x_axis == 'EM_gain':
                        EM_gain=space(self.EM_gain.min,self.EM_gain.max)
                    elif x_axis == 'acquisition_time':
                        acquisition_time=space(self.acquisition_time.min,self.acquisition_time.max)
                    elif x_axis == 'RN':
                        RN=space(self.RN.min,self.RN.max)
                    elif x_axis == 'CIC_charge':
                        CIC_charge=space(0.001,self.CIC_charge.max)
                    elif x_axis == 'Dark_current':
                        Dark_current=space(self.Dark_current.min,self.Dark_current.max)
                    elif x_axis == 'readout_time':
                        readout_time=space(self.exposure.min,10*exposure_time)
                    elif x_axis == 'smearing':
                        smearing=space(self.smearing.min,self.smearing.max)
                    elif x_axis == 'Redshift':
                        Redshift=space(self.Redshift.min,self.Redshift.max)
                    elif x_axis == 'temperature':
                        temperature=space(self.temperature.min,self.temperature.max)
                    elif x_axis == 'QE':
                        QE=space(self.QE.min,self.QE.max)
                    elif x_axis == 'PSF_RMS_mask':
                        PSF_RMS_mask=space(self.fwhm.min,self.fwhm.value[1])
                    elif x_axis == 'PSF_RMS_det':
                        PSF_RMS_det=space(self.fwhm.min,self.fwhm.max)
                    elif x_axis == 'extra_background':
                        extra_background=space(self.extra_background.min,self.extra_background.max)
                    elif x_axis == 'Δx':
                        Δx=space(self.Δx.min,self.Δx.max)#pixels
                    elif x_axis == 'Δλ':
                        Δλ=space(self.Δλ.min,self.Δλ.max)#pixels
                    elif x_axis == 'Throughput':
                        Throughput=space(self.Throughput.min,self.Throughput.max)
                    elif x_axis == 'Atmosphere':
                        Atmosphere=space(self.Atmosphere.min,self.Atmosphere.max)
                    elif x_axis == 'Line_width':
                        Line_width=space(self.Line_width.min,self.Line_width.max)
                    elif x_axis == 'Size_source':
                        Size_source =  space(self.Size_source.min,self.Size_source.max)
                    elif x_axis == 'pixel_scale':
                        pixel_scale = space(self.pixel_scale.min,self.pixel_scale.max)
                    elif x_axis == 'pixel_size':
                        pixel_size = space(self.pixel_size.min,self.pixel_size.max)
                    elif x_axis == 'wavelength':
                        wavelength =space(self.wavelength.min,self.wavelength.max)
                    elif x_axis == 'Slitwidth':
                        Slitwidth =space(self.SlitDims.min,self.SlitDims.max)
                    elif x_axis == 'Slitlength':
                        Slitlength =space(self.SlitDims.min,self.SlitDims.max)
                    elif x_axis == 'Spectral_resolution':
                        Spectral_resolution = space(self.Spectral_resolution.min,self.Spectral_resolution.max)
                    elif x_axis == 'dispersion':
                        dispersion = space(self.dispersion.min,self.dispersion.max)
                    elif x_axis == 'lambda_stack':
                        lambda_stack = space(self.dispersion.value,self.Bandwidth)
                    elif x_axis == 'Collecting_area':
                        Collecting_area = space(10**self.Collecting_area.min,10**self.Collecting_area.max)
                    elif x_axis == "cosmic_ray_loss_per_sec":
                        cosmic_ray_loss_per_sec=space(0,1/exposure[1])
                    elif x_axis == "cosmic_ray_loss_per_sec":
                        cosmic_ray_loss_per_sec=space(0,1/exposure[1])
                    elif x_axis == "Bandwidth":
                        self.Bandwidth=space(0,10*self.Bandwidth)



                    # if x_axis == 'exposure_time':
                    #     dict_val["exposure_time"] = space(1, self.time_max)
                    # elif x_axis == 'Sky':
                    #     dict_val["Sky"] = space(10**self.Sky.min, 10**self.Sky.max)
                    # elif x_axis == 'Signal':
                    #     dict_val["Signal"] = space(10**self.Signal.min, 10**self.Signal.max)
                    # elif x_axis == 'EM_gain':
                    #     dict_val["EM_gain"] = space(self.EM_gain.min, self.EM_gain.max)
                    # elif x_axis == 'acquisition_time':
                    #     dict_val["acquisition_time"] = space(self.acquisition_time.min, self.acquisition_time.max)
                    # elif x_axis == 'RN':
                    #     dict_val["RN"] = space(self.RN.min, self.RN.max)
                    # elif x_axis == 'CIC_charge':
                    #     dict_val["CIC_charge"] = space(0.001, self.CIC_charge.max)
                    # elif x_axis == 'Dark_current':
                    #     dict_val["Dark_current"] = space(self.Dark_current.min, self.Dark_current.max)
                    # elif x_axis == 'readout_time':
                    #     dict_val["readout_time"] = space(self.exposure.min, 10 * dict_val.get("exposure_time", 1))  # fallback si non défini
                    # elif x_axis == 'smearing':
                    #     dict_val["smearing"] = space(self.smearing.min, self.smearing.max)
                    # elif x_axis == 'Redshift':
                    #     dict_val["Redshift"] = space(self.Redshift.min, self.Redshift.max)
                    # elif x_axis == 'temperature':
                    #     dict_val["temperature"] = space(self.temperature.min, self.temperature.max)
                    # elif x_axis == 'QE':
                    #     dict_val["QE"] = space(self.QE.min, self.QE.max)
                    # elif x_axis == 'PSF_RMS_mask':
                    #     dict_val["PSF_RMS_mask"] = space(self.fwhm.min, self.fwhm.value[1])
                    # elif x_axis == 'PSF_RMS_det':
                    #     dict_val["PSF_RMS_det"] = space(self.fwhm.min, self.fwhm.max)
                    # elif x_axis == 'extra_background':
                    #     dict_val["extra_background"] = space(self.extra_background.min, self.extra_background.max)
                    # elif x_axis == 'Δx':
                    #     dict_val["Δx"] = space(self.Δx.min, self.Δx.max)
                    # elif x_axis == 'Δλ':
                    #     dict_val["Δλ"] = space(self.Δλ.min, self.Δλ.max)
                    # elif x_axis == 'Throughput':
                    #     dict_val["Throughput"] = space(self.Throughput.min, self.Throughput.max)
                    # elif x_axis == 'Atmosphere':
                    #     dict_val["Atmosphere"] = space(self.Atmosphere.min, self.Atmosphere.max)
                    # elif x_axis == 'Line_width':
                    #     dict_val["Line_width"] = space(self.Line_width.min, self.Line_width.max)
                    # elif x_axis == 'Size_source':
                    #     dict_val["Size_source"] = space(self.Size_source.min, self.Size_source.max)
                    # elif x_axis == 'pixel_scale':
                    #     dict_val["pixel_scale"] = space(self.pixel_scale.min, self.pixel_scale.max)
                    # elif x_axis == 'pixel_size':
                    #     dict_val["pixel_size"] = space(self.pixel_size.min, self.pixel_size.max)
                    # elif x_axis == 'wavelength':
                    #     dict_val["wavelength"] = space(self.wavelength.min, self.wavelength.max)
                    # elif x_axis == 'Slitwidth':
                    #     dict_val["Slitwidth"] = space(self.SlitDims.min, self.SlitDims.max)
                    # elif x_axis == 'Slitlength':
                    #     dict_val["Slitlength"] = space(self.SlitDims.min, self.SlitDims.max)
                    # elif x_axis == 'Spectral_resolution':
                    #     dict_val["Spectral_resolution"] = space(self.Spectral_resolution.min, self.Spectral_resolution.max)
                    # elif x_axis == 'dispersion':
                    #     dict_val["dispersion"] = space(self.dispersion.min, self.dispersion.max)
                    # elif x_axis == 'lambda_stack':
                    #     dict_val["lambda_stack"] = space(self.dispersion.value, self.Bandwidth)
                    # elif x_axis == 'Collecting_area':
                    #     dict_val["Collecting_area"] = space(10**self.Collecting_area.min, 10**self.Collecting_area.max)
                    # elif x_axis == "cosmic_ray_loss_per_sec":
                    #     dict_val["cosmic_ray_loss_per_sec"] = space(0, 1 / dict_val.get("exposure_time", 1))
                    # elif x_axis == "Bandwidth":
                    #     self.Bandwidth = space(0, 10 * self.Bandwidth)
                    #     dict_val["Bandwidth"] = self.Bandwidth


                dict_val = {
                    "exposure_time": exposure[1],
                    "Sky": Sky,
                    "Signal": Signal,
                    "EM_gain": EM_gain,
                    "acquisition_time": acquisition_time,
                    "RN": RN,
                    "CIC_charge": CIC_charge,
                    "Dark_current": Dark_current,
                    "readout_time": exposure[0],
                    "smearing": smearing,
                    "Redshift": Redshift,
                    "temperature": temperature,
                    "QE": QE,
                    "PSF_RMS_mask": fwhm[0],
                    "PSF_RMS_det": fwhm[1],
                    "extra_background": extra_background,
                    "Δx": Δx,
                    "Δλ": Δλ,
                    "Throughput": Throughput,
                    "Atmosphere": Atmosphere,
                    "Line_width": Line_width,
                    "Size_source": Size_source,
                    "pixel_scale": pixel_scale,
                    "pixel_size": pixel_size,
                    "wavelength": wavelength,
                    "Slitwidth": SlitDims[0],
                    "Slitlength": SlitDims[1],
                    "Spectral_resolution": Spectral_resolution,
                    "dispersion": dispersion,
                    "lambda_stack": lambda_stack,
                    "Collecting_area": Collecting_area,
                    "cosmic_ray_loss_per_sec": cosmic_ray_loss_per_sec,
                    "Bandwidth": self.Bandwidth,
                }

                for vname in ["CIC_charge", "dispersion", "pixel_scale"]:
                    if vname not in [ d.split('=')[0].strip() for d in dependencies]:
                        rsetattr(self, '%s.layout.visibility'%(vname), 'visible'  )  
                if len(dependencies)>0: #TODO besure that we remnove it when one is seleted but not the others 

                    for dep in dependencies:
                        exec(dep, dict_val)
                        var_name = dep.split('=')[0].strip()
                        if var_name == "CIC_charge":
                            CIC_charge = dict_val[var_name]
                            rsetattr(self, '%s.layout.visibility'%(var_name), 'hidden'  ) 
                        if var_name == "dispersion":
                            dispersion = dict_val[var_name]
                            rsetattr(self, '%s.layout.visibility'%(var_name), 'hidden'  )  
                        if var_name == "pixel_scale":
                            pixel_scale = dict_val[var_name]
                            rsetattr(self, '%s.layout.visibility'%(var_name), 'hidden'  )  

                        # print(EM_gain, CIC_charge)              
                # else:
                #     rsetattr(self, '%s.layout.visibility'%("CIC_charge"), 'visible'  )                
                #     rsetattr(self, '%s.layout.visibility'%("dispersion"), 'visible'  )                
                #     rsetattr(self, '%s.layout.visibility'%("pixel_scale"), 'visible'  )                

                        

                

                args, _, _, locals_ = inspect.getargvalues(inspect.currentframe())
                try:
                    new_value = locals_[x_axis]
                except KeyError:
                    new_value = getattr(self,x_axis)
                arg = np.argmin(abs(new_value - value))
                self.new = Observation(instruments=self.instruments, instrument=self.instrument.value,  Throughput_FWHM=Throughput_FWHM, Redshift=Redshift, exposure_time=exposure_time,Sky=Sky, acquisition_time=acquisition_time,counting_mode=counting_mode,Signal=Signal,EM_gain=EM_gain,RN=RN,CIC_charge=CIC_charge,Dark_current=Dark_current,readout_time=readout_time,smearing=smearing,extra_background=extra_background,i=arg,PSF_RMS_mask=PSF_RMS_mask,PSF_RMS_det=PSF_RMS_det,QE=QE,cosmic_ray_loss_per_sec=cosmic_ray_loss_per_sec, Throughput=Throughput, Atmosphere=Atmosphere,lambda_stack=lambda_stack,Slitwidth=Slitwidth, Bandwidth=self.Bandwidth,Size_source=Size_source,Collecting_area=Collecting_area,Δx=Δx,Δλ=Δλ,
                pixel_scale=pixel_scale, Spectral_resolution=Spectral_resolution,  dispersion=dispersion,Line_width=Line_width,wavelength=wavelength,  pixel_size=pixel_size,len_xaxis=self.len_xaxis, Slitlength=Slitlength,IFS=IFS,SNR_res=SNR_res,spectrograph=spectrograph,test=test)
                # self.new = Observation(
                #     instruments=self.instruments,
                #     instrument=self.instrument.value,
                #     Throughput_FWHM=Throughput_FWHM,
                #     Redshift=dict_val["Redshift"],
                #     exposure_time=dict_val["exposure_time"],
                #     Sky=dict_val["Sky"],
                #     acquisition_time=dict_val["acquisition_time"],
                #     counting_mode=counting_mode,
                #     Signal=dict_val["Signal"],
                #     EM_gain=dict_val["EM_gain"],
                #     RN=dict_val["RN"],
                #     CIC_charge=dict_val["CIC_charge"],
                #     Dark_current=dict_val["Dark_current"],
                #     readout_time=dict_val["readout_time"],
                #     smearing=dict_val["smearing"],
                #     extra_background=dict_val["extra_background"],
                #     i=arg,
                #     PSF_RMS_mask=dict_val["PSF_RMS_mask"],
                #     PSF_RMS_det=dict_val["PSF_RMS_det"],
                #     QE=dict_val["QE"],
                #     cosmic_ray_loss_per_sec=dict_val["cosmic_ray_loss_per_sec"],
                #     Throughput=dict_val["Throughput"],
                #     Atmosphere=dict_val["Atmosphere"],
                #     lambda_stack=dict_val["lambda_stack"],
                #     Slitwidth=dict_val["Slitwidth"],
                #     Bandwidth=dict_val["Bandwidth"],
                #     Size_source=dict_val["Size_source"],
                #     Collecting_area=dict_val["Collecting_area"],
                #     Δx=dict_val["Δx"],
                #     Δλ=dict_val["Δλ"],
                #     pixel_scale=dict_val["pixel_scale"],
                #     Spectral_resolution=dict_val["Spectral_resolution"],
                #     dispersion=dict_val["dispersion"],
                #     Line_width=dict_val["Line_width"],
                #     wavelength=dict_val["wavelength"],
                #     pixel_size=dict_val["pixel_size"],
                #     len_xaxis=self.len_xaxis,
                #     Slitlength=dict_val["Slitlength"],
                #     IFS=IFS,
                #     SNR_res=SNR_res,
                #     spectrograph=spectrograph,
                #     test=test
                # )
                self.colors=self.new.colors

                
                if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out1): 
                    ft=8
                    
                    try:
                        arg = np.argmin(abs(getattr(self.new,x_axis) - value))
                    except AttributeError:
                        arg= np.argmin(temperature - value)
                    try:
                        label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f, SNR=%0.2f'%(self.x_axis.value,float_to_latex(value),float_to_latex(new_value[np.nanargmax(self.new.SNR)]),self.new.SNR[arg],np.nanmax(self.new.SNR))
                    except (TypeError,ValueError) as e:
                        if ("FIREBall" in self.instrument.value) & (self.counting_mode.value):
                            label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f, SNR=%0.2f\nT=%0.1f sigma\nSignal kept=%i%%, RN kept=%i%%'%(self.x_axis.value,float_to_latex(value),float_to_latex(new_value[np.nanargmax(self.new.SNR)]),self.new.SNR[arg],np.nanmax(self.new.SNR),self.new.n_threshold[arg], 100*self.new.Photon_fraction_kept[arg], 100*self.new.RN_fraction_kept[arg])#, self.new.gain_thresholding[arg])
                        else:
                            label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f, SNR=%0.2f'%(self.x_axis.value,float_to_latex(value),float_to_latex(new_value[np.nanargmax(self.new.SNR)]),self.new.SNR[arg],np.nanmax(self.new.SNR))#, self.new.gain_thresholding[arg])

                    max_,min_=[],[]

                    for i,name in enumerate(self.new.names): 
                        self.ax0.lines[i].set_xdata(new_value)
                        self.ax0.lines[i].set_ydata(self.new.noises_per_exp[:,i] )
                        if self.new.percents[i,self.new.i] ==np.max(self.new.percents[:,self.new.i]):
                            self.ax0.lines[i].set_label(r"$\bf{➛%s}$: %0.2f (%0.1f%%)"%(name,self.new.noises_per_exp[self.new.i,i],self.new.percents[i,self.new.i]))
                        else:
                            self.ax0.lines[i].set_label('%s: %0.2f (%0.1f%%)'%(name,self.new.noises_per_exp[self.new.i,i],self.new.percents[i,self.new.i]))
                        max_.append(np.nanmax(self.new.noises_per_exp[:,i]))
                        min_.append(np.nanmin(self.new.noises_per_exp[:,i]))
                    self.ax0.lines[i+1].set_xdata(new_value)
                    self.ax0.lines[i+1].set_ydata(   np.sqrt(np.nansum(np.multiply(self.new.noises_per_exp[:,:-1],self.new.noises_per_exp[:,:-1]),axis=1))   )


                    
                    self.ax0.lines[i+1].set_label('%s: %0.2f'%("Quadratic sum",np.sqrt(np.nansum(np.multiply(self.new.noises_per_exp[self.new.i,:-1],self.new.noises_per_exp[self.new.i,:-1])))   ))     
                    self.ax0.legend(loc='upper left', fontsize=ft,title_fontsize=ft )
                    if x_axis in ["exposure_time","readout_time","PSF_RMS_mask","PSF_RMS_det","Slitwidth","Slitlength"]:
                        self.ax3.set_xlabel(x_axis.replace("_"," "))
                    else:
                        try:
                            self.ax3.set_xlabel(rgetattr(self, '%s.description_tooltip'%(x_axis)) )
                        except AttributeError:
                            self.ax3.set_xlabel(x_axis + "  [%s]"%(self.instruments["Unit"][self.instruments["Charact."]==x_axis][0]))

                    self.ax3.lines[0].set_data(new_value,  np.log10(self.new.SB_lim_per_pix))
                    self.ax3.lines[1].set_data(new_value,  np.log10(self.new.SB_lim_per_res))
                    self.ax3.lines[2].set_data(new_value,  np.log10(self.new.SB_lim_per_source))

                    self.ax3.lines[0].set_label("SNR=5 limiting SB/power per pixel (%0.2fergs/s~%0.1ELU)"%(np.log10(1.30e57*self.new.SB_lim_per_pix[self.new.i]),   convert_ergs2LU(self.new.SB_lim_per_pix[self.new.i],self.wavelength.value)   ))
                    self.ax3.lines[1].set_label("SNR=5 limiting SB/power per elem resolution (%0.2fergs/s~%0.1ELU)"%(np.log10(1.30e57*self.new.SB_lim_per_res[self.new.i]),   convert_ergs2LU(self.new.SB_lim_per_res[self.new.i],self.wavelength.value)   ))
                    self.ax3.lines[2].set_label("SNR=5 limiting SB/power per source (%0.2fergs/s~%0.1ELU)"%(np.log10(1.30e57*self.new.SB_lim_per_source[self.new.i]),   convert_ergs2LU(self.new.SB_lim_per_source[self.new.i],self.wavelength.value)   ))
                    self.ax3.set_ylim(np.nanmin([np.nanmin(np.log10(0.1*self.new.SB_lim_per_source)),np.nanmin(0.1*np.log10(self.new.SB_lim_per_res))]),np.nanmax(np.log10(self.new.SB_lim_per_pix[np.isfinite(self.new.SB_lim_per_pix)])))


                    for v in self.v:
                        v.set_xdata([value,value])
                    self.v[-2].set_label(label)


                    for artist in self.ax2.collections+self.ax1.collections:
                        artist.remove()

                    self.stackplot2 = self.ax2.stackplot(new_value,self.new.SNR * np.array(self.new.noises).T[:-1,:]**2/self.new.Total_noise_final**2,alpha=0.7,colors=self.colors)
                    # self.stackplot2 = ax3.stackplot(getattr(self,x), self.SNR * np.array(self.noises).T[:-1,:]**2/self.Total_noise_final**2,alpha=0.7,colors=self.colors)
      
                                
                    self.ax0.set_ylabel('Noise (e-/%s/exp)'%(SNR_res.split(" ")[1]))     
                    self.ax1.set_ylabel('Contrib (e-/%s/exp)'%(SNR_res.split(" ")[1]))     
                    self.ax2.set_ylabel('SNR (%s, N frames)'%(SNR_res.split(" ")[1]))     


                    # labels =  ['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]
                    if (type(self.new.number_pixels_used)==np.float64) | (type(self.new.number_pixels_used) is int)| (type(self.new.number_pixels_used) is float): 
                        self.number_pixels_used = self.new.number_pixels_used
                    else:
                        self.number_pixels_used = self.new.number_pixels_used[arg]
                    contributions = self.number_pixels_used * getattr(self.new, "electrons_per_pix")[self.new.i, :]

                    percentages = 100 * contributions / np.sum(contributions)
                    # Find the index of the maximum contribution
                    max_index = np.argmax(percentages)
                    # Generate the labels, making the largest contribution bold
                    labels = [
                        r'%s: %0.3f (%0.1f%%)' % (name, contributions[j], percentages[j])
                        if j != max_index
                        else r"➛$\mathbf{%s}$: %0.3f (%0.1f%%)" % (name, contributions[j], percentages[j])
                        for j, name in enumerate(self.new.names)
                    ]                    
                    self.stackplot1 = self.ax1.stackplot(new_value, self.number_pixels_used * np.array(self.new.electrons_per_pix).T,alpha=0.7,colors=self.colors,labels=labels)

                    self.ax1.legend(loc='upper left',title="Overall background: %0.3f (%0.1f%%)"%(np.nansum(self.number_pixels_used *self.new.electrons_per_pix[self.new.i,1:]),100*np.nansum(self.number_pixels_used * self.new.electrons_per_pix[self.new.i,1:])/np.nansum(self.number_pixels_used * self.new.electrons_per_pix[self.new.i,:])),fontsize=ft,title_fontsize=ft)
                    self.ax2.legend(loc='upper right', fontsize=ft,title_fontsize=ft )
                    self.ax3.legend(loc='upper right',title=r"$\mathbf{Left}$: Ext. Source surface brightness, $\mathbf{Right}$: Point source power", fontsize=ft,title_fontsize=ft )
                    self.ax2.set_xlim((np.max([np.min(new_value),1e-6]),np.max(new_value)))
                    self.ax2.set_xlim((np.min(new_value),np.max(new_value)))


                    self.ax0.set_yscale(self.yscale)
                    self.ax1.set_yscale(self.yscale)
                    self.ax2.set_yscale(self.yscale)
                    if log:
                        # self.ax0.set_ylim(ymin=np.nanmin(self.new.noises[:,:-1]/self.new.factor[:,None]),ymax=np.nanmax(np.nansum(self.new.noises[:,:-1],axis=1)/self.new.factor))
                        self.ax0.set_ylim(ymin=0,ymax=   np.nanmax(  np.sqrt(np.nansum(np.multiply(self.new.noises_per_exp[:,:-1],self.new.noises_per_exp[:,:-1]),axis=1))  )   )
                        self.ax1.set_ylim(ymin=np.nanmin(self.number_pixels_used * np.array(self.new.electrons_per_pix[:,0])),ymax=np.max(self.number_pixels_used *np.sum(getattr(self.new,'electrons_per_pix'),axis=1)))
                        self.ax2.set_ylim(ymin=np.nanmin(np.array( self.new.SNR * np.array(self.new.noises).T[:-1,:]**2/self.new.Total_noise_final**2)[:,0]),ymax=np.nanmax(getattr(self.new,'SNR')))

                        # if SNR:
                        # else:
                        #     self.ax2.set_ylim(ymin=np.nanmin(np.array( self.new.snrs_per_pixel * np.array(self.new.noises).T[:-1,:]**2/self.new.Total_noise_final**2)[:,0]),ymax=np.nanmax(getattr(self.new,'snrs_per_pixel')))

                    else:
                        self.ax0.set_ylim((-0.1,np.nanmax(  np.sqrt(np.nansum(np.multiply(self.new.noises_per_exp[:,:-1],self.new.noises_per_exp[:,:-1]),axis=1))  )   ))
                        self.ax1.set_ylim((0,  self.number_pixels_used *  np.max(np.sum(getattr(self.new,'electrons_per_pix'),axis=1))))
                        self.ax2.set_ylim((0,np.nanmax(self.new.SNR)))
                        # if SNR:
                        #     self.ax2.set_ylim((0,np.nanmax(self.new.SNR)))
                        # else:
                        #     self.ax2.set_ylim((0,np.nanmax(self.new.snrs_per_pixel)))

                    
                    if xlog:
                        try:
                            if (np.nanmin(new_value)<=0) & (np.nanmin(new_value)  < value < np.nanmax(new_value)  ):
                            # if (rgetattr(self,"%s.min"%(self.x_axis.value))<=0) & ( rgetattr(self,"%s.min"%(self.x_axis.value))  <  rgetattr(self,"%s.value"%(self.x_axis.value)) <  rgetattr(self,"%s.max"%(self.x_axis.value))):
                                self.ax0.set_xscale("symlog")  
                            else:
                                self.ax0.set_xscale("log")
                        except AttributeError:
                            self.ax0.set_xscale("log")

                    else:
                        self.ax0.set_xscale("linear")
                    #TODO add that if dispersion is an array we must 
                    disp = dispersion[arg] if type(dispersion)==np.ndarray else dispersion

                    title = '%s : FOV=%0.1famin$^2$, λ=%inm, Total Throughput=%i%%, Effective area=%0.1fcm$^2$, Platescale=%.1f,  PSF$_{x,λ}$=%0.1f, %0.1f pix, Npix = %i '%(self.instrument.value,self.FOV_size, self.wavelength.value, 100*self.Throughput.value*self.QE.value*self.Atmosphere.value, 100*100*self.Throughput.value*self.QE.value*self.Atmosphere.value*self.Collecting_area.value,  pixel_scale, 2.35*self.PSF_RMS_det/pixel_scale, 10*self.wavelength.value/self.Spectral_resolution.value/disp, self.number_pixels_used    )
                    # title = 'Instrument=%s, FOV=%0.1famin$^2$, λ=%inm, Throughput=%i%%, Atm=%i%%, Platescale=%.1f, area=%0.1fm$^2$'%(instrument,self.instruments[instrument][self.instruments["Charact."]=="FOV_size"][0], self.instruments[instrument][self.instruments["Charact."]=="wavelength"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Throughput"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Atmosphere"][0], self.instruments[instrument][self.instruments["Charact."]=="pixel_scale"][0], self.instruments[instrument][self.instruments["Charact."]=="Collecting_area"][0])
                    self.ax0.set_title(title,y=0.97,fontsize=10)

                    self.fig.canvas.draw()
                    # self.fig2.canvas.draw()

                else:


                    if self.spectrograph.value:
                        self.Redshift.disabled = False if "Rest-frame" in self.spectra.value else True
                    else:
                        self.Redshift.disabled = True if "Gaussian" in self.spectra.value else False
                    
                    self.Size_source.disabled = True if "cube" in self.spectra.value.lower() else False
                    self.Line_width.disabled = False if "baseline" in self.spectra.value.lower() else True

                    # if ("CIV" in spectra)| ("CIII" in spectra)| ("OVI" in spectra)| ("Lyα" in spectra):
                    #     self.Redshift.value =   wavelength*10 /float(re.search(r'\d+', spectra).group()) -1
                    if "Spectra mNUV=" in x_axis:
                        Signal = float(spectra.split("=")[-1])
                        Size_source = 0.1
                        self.fraction_lya.layout.visibility = 'visible'
                    else: 
                        Signal = 20
                        Size_source = 4
                        self.fraction_lya.layout.visibility = 'hidden'
                    Sky = (self.new.Sky/exposure_time)#[arg]
                    flux = (self.new.Signal_el/exposure_time)#[arg]
                    IFS = True if self.output_tabs.get_state()["selected_index"]==2 else False
                    # try:
                        # self.im,self.im_stack, self.cube_stack, self.im0, source_im_wo_atm, self.imaADU_stack_only_source, self.imaADU_without_source, self.imaADU_stack_without_source, self.imaADU_source = self.new.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source=spectra,size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=QElambda,atmlambda=atmlambda,fraction_lya= fraction_lya, Full_well=self.Full_well, conversion_gain=self.conversion_gain, Altitude=self.Altitude,Throughput_FWHM=self.Throughput_FWHM.value,sky_lines=self.sky_lines.value, Redshift=self.Redshift.value,IFS=IFS)
                    self.im,self.im_stack = self.new.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source=spectra,size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=QElambda,atmlambda=atmlambda,fraction_lya= fraction_lya, Full_well=self.Full_well, conversion_gain=self.conversion_gain, Altitude=self.Altitude,Throughput_FWHM=self.Throughput_FWHM.value,sky_lines=self.sky_lines.value, Redshift=self.Redshift.value,IFS=IFS,source_image=self.source_im.value)
                    self.cube_detector, self.cube_detector_stack = self.im,self.im_stack
                    # except ValueError as e:
                        # self.cube_detector, self.cube_detector_stack = self.new.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source=spectra,size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=QElambda,atmlambda=atmlambda,fraction_lya= fraction_lya, Full_well=self.Full_well, conversion_gain=self.conversion_gain, Altitude=self.Altitude,Throughput_FWHM=self.Throughput_FWHM.value,sky_lines=self.sky_lines.value, Redshift=self.Redshift.value,IFS=IFS,source_image=self.source_im.value)


                    self.bins=np.linspace(np.nanmin(self.im),np.nanmax(self.im),100)

                    self.current_cmap.set_bad('red',1.)
                    gain_unit = 1 if self.counting_mode else EM_gain
                    if units=="ADU/frame": #ADU/frame: ok basic
                        factor=1
                    elif units=="e-/frame": #e-/frame: divide by conversion gain and amplification gain
                        factor=1/self.conversion_gain/gain_unit
                    elif units=="photons/frame": #photons/frame: account for QE
                        factor=1/self.conversion_gain/gain_unit/QE
                    elif units=="e-/hour": #e-/hour: divide by exptime
                        factor=3600/self.conversion_gain/gain_unit/exposure_time
                    elif units=="photons/hour": #photons/hour: divide by exptime, account for QE
                        factor=3600/self.conversion_gain/gain_unit/QE/exposure_time
                    elif units=="e-/second": #e-/hour: divide by exptime
                        factor=self.conversion_gain/gain_unit/exposure_time
                    elif units=="photons/second": #photons/hour: divide by exptime, account for QE
                        factor=self.conversion_gain/gain_unit/QE/exposure_time

                if self.source_im.value == "Sim image":
                    t1, t2 = 'Single image', 'Stacked image'
                    t12, t22 = 'Single cube', 'Stacked cube'
                elif self.source_im.value == "Source":
                    t1, t2 = r'$F_S$  (e-/pix/exp)', r'$F_S$ (e-/pix/exp)'
                    t12, t22 = r'$F_S$ (erg/s/cm²/arcsec²/Å)', r'$F_S$ (e-/pix/exp)'
                elif self.source_im.value == "Convolved source":
                    t1, t2 =  r'$F_S$ (e-/pix/exp)', r'$F_S$ (e-/pix/exp)'
                    t12, t22 =  r'$F_S$ (e-/pix/exp)', r'$F_S$ (e-/pix/exp)'
                elif self.source_im.value == "SNR":
                    t1, t2 = 'SNR (single)', 'SNR (stack)'
                    t12, t22 = 'SNR (single)', 'SNR (stack)'                
                self.f = lambda x: wavelength + (dispersion/10) * (x - n1/2) #if spectro else lambda x:x
                self.g = lambda x: (x - wavelength) / (dispersion/10) + n1/2 #if spectro else lambda x:x
                self.vmin, self.vmax = self.minmax.value
                if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out2): 
                    with self.out2:
                        self.nax1_secondary.remove()
                        self.nax1_secondary = self.nax1.secondary_xaxis("top", functions=(self.f,self.g)) if self.spectrograph.value else self.nax1.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
                        if self.spectrograph.value:
                            self.slit_text.txt.set_text("Slit=%0.1f'' × %0.1f''\n    = %0.2f × %0.2fkpc$^2$"%(self.Slitwidth,self.Slitlength,self.Slitwidth/cosmo.arcsec_per_kpc_proper(Redshift+0.001).value,self.Slitlength/cosmo.arcsec_per_kpc_proper(Redshift+0.001).value))
                            self.nax1.set_title('Wavelength (nm) - R = %i (λ/dλ) = %0.1fÅ = %0.1f km/s'%(Spectral_resolution, 10*wavelength/Spectral_resolution,299792.458/Spectral_resolution),fontsize=10)
                        else:
                            self.slit_text.txt.set_text("")
                            self.nax1.set_title('',fontsize=10)#Wavelength (nm) - R = %i (λ/dλ) = %0.1fÅ = %0.1f km/s'%(Spectral_resolution, 10*wavelength/Spectral_resolution,299792.458/Spectral_resolution))
                        if (spectra =="Observed-frame: Baseline Spectra") & self.spectrograph.value:
                            self.l2_s[0].set_label("Stack. spectral prof \nFWHM$_v$ = %0.1f km/s"%(Line_width * 299792.458/(10*wavelength )))
                        else:
                            self.l2_s[0].set_label("Stack. spectral prof")

                        self.nan_if_imager = 1 if self.spectrograph.value else np.nan
                        self.absorption[0].set_ydata(  self.Atmosphere.value*   self.new.atm_trans*np.ones(n1))             
                        self.emission_lines[0].set_ydata(  (self.new.final_sky*np.ones(n1)   - np.min(self.new.final_sky)) / np.ptp(self.new.final_sky*np.ones(n1)   - np.min(self.new.final_sky))) 
                        self.Throughput_curve[0].set_ydata(  self.new.Throughput_curve*np.ones(n1)  )#  / np.max(self.new.Throughput_curve))   #
                        self.nax2_1_secondary.remove()
                        self.nax2_1_secondary = self.nax2_1.secondary_xaxis("top", functions=(self.f,self.g)) 
                        self.nax2_1.set_yscale(self.yscale)
                        norm = "linear" if log==False else SymLogNorm(linthresh=1*factor) #
                                
                        self.mod = mostFrequent(self.im_stack[:20,:].flatten())
                        # self.limit = self.mod+threshold*RN
                        self.limit = self.mod+self.new.n_threshold * RN

                        if log==False:
                            im = self.nax.imshow(self.im*factor, aspect="auto",cmap=self.current_cmap,interpolation=interpolation ,vmin=np.nanmin(self.im*factor) + self.vmin*np.ptp(self.im*factor),vmax= np.nanmax(self.im*factor)-  (1-self.vmax)*np.ptp(self.im*factor))
                            im0 = self.nax0.imshow(self.im_stack*factor, aspect="auto",cmap=self.current_cmap,interpolation=interpolation ,vmin=np.nanmin(self.im_stack*factor) + self.vmin*np.ptp(self.im_stack*factor),vmax=np.nanmax(self.im_stack*factor)-  (1-self.vmax)*np.ptp(self.im_stack*factor))
                        else:
                            im = self.nax.imshow(self.im*factor, aspect="auto",cmap=self.current_cmap,interpolation=interpolation ,norm=SymLogNorm(linthresh=1*factor,vmin=np.nanmin(self.im*factor) + self.vmin*np.ptp(self.im*factor),vmax= np.nanmax(self.im*factor)-  (1-self.vmax)*np.ptp(self.im*factor)))
                            im0 = self.nax0.imshow(self.im_stack*factor, aspect="auto",cmap=self.current_cmap,interpolation=interpolation,norm=SymLogNorm(linthresh=1*factor ,vmin=np.nanmin(self.im_stack*factor) + self.vmin*np.ptp(self.im_stack*factor),vmax=np.nanmax(self.im_stack*factor)-  (1-self.vmax)*np.ptp(self.im_stack*factor)))

                        self.cbar1 = self.fig2.colorbar(im, cax=self.cax, orientation='horizontal')
                        self.cbar2 = self.fig2.colorbar(im0, cax=self.cax0, orientation='horizontal')
                        if log==False:
                            self.cbar1.formatter.set_powerlimits((0, 0))
                            self.cbar2.formatter.set_powerlimits((0, 0))
                            self.cbar1.formatter.set_useMathText(True)
                            self.cbar2.formatter.set_useMathText(True)
                        # else:
                        #     self.cbar1.set_yscale('log')
                        #     self.cbar2.set_yscale('log')

                        labels =  ['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]

                        stacked_profile = np.nanmean(im0.get_array().data[:,int(n1/2 - Line_width/dispersion/2):int(n1/2 + Line_width/dispersion/2)],axis=1)
                        spatial_profile = np.nanmean(im.get_array().data[:,int(n1/2 - Line_width/dispersion/2):int(n1/2 + Line_width/dispersion/2)],axis=1)
                        self.nax1.set_yscale(self.yscale)
                        self.nax1.lines[0].set_ydata(spatial_profile)
                        self.nax1.lines[1].set_ydata(np.nanmean(im.get_array().data[int(n2/2 - Size_source/pixel_scale/2):int(n2/2 + Size_source/pixel_scale/2),:],axis=0)) 
                        self.nax.lines[0].set_label("  \n".join(labels)) 
                        self.nax.legend(loc="upper left",handlelength=0, handletextpad=0, fancybox=True,markerscale=0,fontsize=8)
                
                        try:
                            self.popt, self.pcov = curve_fit(gaus,np.arange(len(stacked_profile)),stacked_profile,p0=[np.ptp(stacked_profile), 50, 5, stacked_profile.min()])
                        except (RuntimeError, ValueError) as e:
                            self.popt = [0,0,0,0]

                        self.noise_res_element = self.im_stack[int(n2/2 - Size_source/pixel_scale/2/2.35):int(n2/2 + Size_source/pixel_scale/2/2.35),int(n1/2 - Line_width/dispersion/2/2.35):int(n1/2 + Line_width/dispersion/2/2.35)].std()/np.sqrt(self.new.number_pixels_used)
                        self.measured_SNR = self.popt[0]**2 / self.noise_res_element
                        
                        self.Flux_ADU =  np.sum(gaus( np.arange(len(stacked_profile)),*self.popt)-self.popt[-1])/factor 
                        # self.Flux_ADU_counting =  np.sum(-np.log(1-( self.fit["function"]( np.arange(len(stacked_profile)),*self.fit["popt"])-self.fit["popt"][-1] )/(np.exp(-threshold*RN/EM_gain))))
                        self.e_s_pix = self.Flux_ADU * self.new.dispersion / exposure_time / self.new.N_images_true/self.conversion_gain  if counting_mode else  self.Flux_ADU * self.new.dispersion / EM_gain / exposure_time/self.conversion_gain
                        self.flux = self.e_s_pix / self.new.Throughput/ self.new.Atmosphere / QE / self.new.Collecting_area
                        photon_energy_erg = 9.93e-12
                        self.mag = -2.5*np.log10(self.flux*photon_energy_erg/(2.06*1E-16))+20.08
                        # Power 2 in SNR as it is in the definition of the gaussian fit.
                        # TODO we do no account for the poisson noise of the source, neither the fact that we use one resolution element



                        self.l1_s[0].set_ydata(stacked_profile) 
                        self.l2_s[0].set_ydata(np.nanmean(im0.get_array().data[int(n2/2 - Size_source/pixel_scale/2):int(n2/2 + Size_source/pixel_scale/2),:],axis=0))
                        self.l3_s[0].set_ydata(gaus( np.arange(len(stacked_profile)),*self.popt))
                        if self.spectrograph.value:
                            self.l3_s[0].set_label("σ=%0.1f'', SNR=%0.2f/%0.2f=%0.2f, mag=%0.1f"%(self.popt[2]*pixel_scale,self.popt[0]**2,self.noise_res_element,self.measured_SNR ,self.mag))
                        else:
                            self.l3_s[0].set_label("σ=%0.1f'', SNR=%0.2f/%0.2f=%0.2f"%(self.popt[2]*pixel_scale,self.popt[0]**2,self.noise_res_element,self.measured_SNR))



                        if self.spectrograph.value:
                            self.nax.set_title(t1+': FOV = %i" × %iÅ, Slit=%0.1f" × %0.1f"'%(100*pixel_scale,500*dispersion, self.Slitwidth,self.Slitlength)) 
                            self.nax0.set_title(t2+': Pixel size = %0.2f" × %0.2fÅ'%(pixel_scale,dispersion))
                        else:
                            self.nax.set_title(t1+': FOV = %i" × %i" '%(100*pixel_scale,500*pixel_scale)) 
                            self.nax0.set_title(t2+': Pixel size = %0.2f" × %0.2f"'%(pixel_scale,pixel_scale))
                        self.nax1.legend(loc="upper right",fontsize=8,title="Averaged profiles")
                        self.nax1.relim()
                        self.nax1.autoscale_view()

                        [b.remove() for b in self.bars2]
                        [b.remove() for b in self.bars1]
                        n_,_,self.bars1 = self.nax2.hist(self.im[3:-3,3:-3].flatten() * factor,bins=self.bins* factor,alpha=0.3,color=self.l1[0].get_color(),label='Single image') #,log=True
                        _,_,self.bars2 = self.nax2.hist(self.im_stack[3:-3,3:-3].flatten()* factor,bins=self.bins* factor,alpha=0.3,color=self.l2[0].get_color(),label='Averaged stack')#log=True,
                        self.nax2.set_yscale("log")                       
                        self.nax2.set_xlim(np.nanmin(self.bins)* factor,np.nanmax(self.bins)* factor)
                        self.nax2.set_ylim(1,np.nanmax(n_))
                        # self.nax2.relim()
                        # self.nax2.autoscale_view()


                        self.hw, self.hl = Slitwidth/2/pixel_scale ,  self.Slitlength/2/pixel_scale
                        if (np.round(self.Slitlength,3)!=np.round(Slitwidth,3)):
                            self.nax.lines[1].set_data(self.nan_if_imager *np.array([n1/2 - self.hw,n1/2 + self.hw,n1/2 + self.hw,n1/2 - self.hw,n1/2 - self.hw]),[n2/2 - self.hl,n2/2 - self.hl,n2/2 + self.hl,n2/2 + self.hl,n2/2 - self.hl])
                            self.nax0.lines[0].set_data(self.nan_if_imager *np.array([n1/2 - self.hw,n1/2 + self.hw,n1/2 + self.hw,n1/2 - self.hw,n1/2 - self.hw]),[n2/2 - self.hl,n2/2 - self.hl,n2/2 + self.hl,n2/2 + self.hl,n2/2 - self.hl])
                        else:
                            self.nax.lines[1].set_data(n1/2 + Slitwidth/2/pixel_scale * self.nan_if_imager * np.cos( np.linspace( 0 , 2 * np.pi , 50 ) ),n2/2+Slitwidth/2/pixel_scale * np.sin( np.linspace( 0 , 2 * np.pi , 50 ) ))
                            self.nax0.lines[0].set_data(n1/2 + Slitwidth/2/pixel_scale * self.nan_if_imager *np.cos( np.linspace( 0 , 2 * np.pi , 50 ) ),n2/2+Slitwidth/2/pixel_scale * np.sin( np.linspace( 0 , 2 * np.pi , 50 ) ))

                        # labels = "  \n".join(['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)])
                        # self.nax0.lines[0].set_label()
                        # self.nax0.legend(loc="upper left", handlelength=0, handletextpad=0, fancybox=True, markerscale=0, fontsize=8, title="Signal contribution")
                        contributions = getattr(self.new, "electrons_per_pix")[self.new.i, :]
                        max_index = np.argmax(100 * contributions / np.sum(contributions))
                        labels = "  \n".join([
                            r'%s: %0.3f (%0.1f%%)' % (name, contributions[j], 100 * contributions[j] / np.sum(contributions))
                            if j != max_index
                            else r"➛$\mathbf{%s}$: %0.3f (%0.1f%%)" % (name, contributions[j], 100 * contributions[j] / np.sum(contributions))
                            for j, name in enumerate(self.new.names)])

                        # Update the label and legend for the plot
                        self.nax0.lines[0].set_label(r'{}'.format(labels))
                        self.nax0.legend(loc="upper left", handlelength=0, handletextpad=0, fancybox=True, markerscale=0, fontsize=8, title="Signal contribution")
                        labels = [     r'%s: %0.2f (%0.1f%%)' % (name, self.new.noises[self.new.i, i] / self.new.factor[self.new.i], self.new.percents[i, self.new.i])     if self.new.percents[i, self.new.i] < np.max(self.new.percents[:, self.new.i])     else r"➛$\mathbf{%s}$: %0.2f (%0.1f%%)" % (name, self.new.noises[self.new.i, i] / self.new.factor[self.new.i], self.new.percents[i, self.new.i])     for i, name in enumerate(self.new.names)]
                        self.nax.lines[0].set_label(r'{}'.format("\n".join(labels)))
                        self.nax.legend(loc="upper left", handlelength=0, handletextpad=0, fancybox=True, markerscale=0, fontsize=8, title="%s Noise"%(self.instrument.value))

                        self.nax2.lines[0].set_xdata([self.mod,self.mod])
                        self.nax2.lines[1].set_xdata([self.limit[arg],self.limit[arg]])
                        if "FIREBall" in  self.instrument.value:
                            try:
                                title = 'Signal kept=%i%%, RN kept=%i%%, Signal/tot=%i%%'%(100*self.new.Photon_fraction_kept[0], 100*self.new.RN_fraction_kept[0],100*(np.mean(self.im_stack[int(n2/2 - self.Size_source.value/pixel_scale/2/2.35):int(n2/2 + self.Size_source.value/pixel_scale/2/2.35),:])-np.mean(self.im_stack[:20,:]))/np.mean(self.im_stack[int(n2/2 - self.Size_source.value/pixel_scale/2/2.35):int(n2/2 + self.Size_source.value/pixel_scale/2/2.35),:]))
                            except IndexError as e:
                                title = 'Signal kept=%i%%, RN kept=%i%%, Signal/tot=%i%%'%(100*self.new.Photon_fraction_kept, 100*self.new.RN_fraction_kept[0],100*(np.mean(self.im_stack[int(n2/2 - self.Size_source.value/pixel_scale/2/2.35):int(n2/2 + self.Size_source.value/pixel_scale/2/2.35),:])-np.mean(self.im_stack[:20,:]))/np.mean(self.im_stack[int(n2/2 - self.Size_source.value/pixel_scale/2/2.35):int(n2/2 + self.Size_source.value/pixel_scale/2/2.35),:]))
                            self.nax2.lines[0].set_label("Bias %0.3f, PC limit %0.3f (%s):\n%s "%(self.mod,self.limit[arg], counting_mode, title))
                        else:
                            self.nax2.lines[0].set_label(" ")
                        # self.l1[0].set_ydata(factor*self.im[:,int(n1/2-self.Line_width.value/self.pixel_scale.value):int(n1/2+self.Line_width.value/self.pixel_scale.value)].mean(axis=1)) 
                        # self.l2[0].set_ydata(factor*self.im[int(n2/2 - self.Size_source.value/self.pixel_scale.value/2/2.35):int(n2/2 + self.Size_source.value/self.pixel_scale.value/2/2.35),:].mean(axis=0)) 

                        self.nax2.legend(loc='upper right',fontsize=8,title="Histogram - Pixels'values")
                        self.fig2.tight_layout()
                        self.fig2.canvas.draw()


                if len(self.output_tabs.children)>2: 
                    if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out3): 
                        with self.out3:
                            if IFS:
                                # TODO change
                                # n3 = int(np.sqrt(60*60*self.instruments_dict[self.instrument.value]["FOV_size"])/self.Slitwidth)
                                n3 = (n2 *pixel_scale) / Slitwidth
                                if n3>n2 : 
                                    n3=n2
                                center = n1/2
                                # print(n1,n2,n3)
                                self.f = lambda x: wavelength + (dispersion/10) * (x - center)
                                self.g = lambda x: (x - wavelength) / (dispersion/10) + center
                                self.nax2s_secondary.remove()
                                self.nax2s_secondary = self.nax2s.secondary_xaxis("top", functions=(self.f,self.g))
                                if hasattr(self, 'cube_detector'):
                                    self.cube_detector  *= factor
                                    # print(self.cube_detector.shape)
                                    self.cube_detector_stack  *= factor
                                    self.im = self.cube_detector 
                                    self.im_stack =  self.cube_detector_stack
                                    # n1,n2,n3 = self.im.shape
                                    # try:
                                    ax, ay, az = self.im.shape

                                    self.ifs_cube = np. transpose (self.cube_detector, (1, 2, 0))  
                                    self.ifs_cube_stack =np.transpose (self.cube_detector_stack, (1, 2, 0))  
                                    test = self.ifs_cube.copy()
                                    # self.ifs_cube_stack = subtract_continuum(self.ifs_cube_stack,self.ifs_cube_stack)
                                    # self.ifs_cube = subtract_continuum(test,test)


                                    self.ifs_spectra[0].set_ydata(self.im[int(ax/2),int(ay/2)+int(Δx),:]) 
                                    self.ifs_spectra_stack[0].set_ydata(self.im_stack[int(ax/2),int(ay/2)+int(Δx),:])
                        
                                    # ifs_integ_spectra_stack =       np.nanmean(self.ifs_cube_stack[int(n3/2-(n3/n2)*Slitwidth/pixel_scale):int(n3/2+(n3/n2)*Slitwidth/pixel_scale),int(n2/2-Slitwidth/pixel_scale):int(n2/2+Slitwidth/pixel_scale),:],axis=(0,1)) 
                                    # self.ifs_integ_spectra_stack[0].set_ydata( ifs_integ_spectra_stack)


                                if self.Slitlength>Slitwidth:
                                    self.nax20.set_title(t12 + ": %i pixels × %i spaxels"%(n2,n3),fontsize=9)
                                    # self.nax21.set_title("%0.1f' × %0.1f', full FOV= %0.1f' × %0.1f'"%(n2/pixel_scale/60,n3/pixel_scale/60,np.sqrt(self.FOV_size),np.sqrt(self.FOV_size)))
                                    self.nax21.set_title(t22+ ": %0.1f' × %0.1f', full FOV= %0.1f' × %0.1f'"%(n2*pixel_scale/60,n3*Slitwidth/60, np.sqrt(self.FOV_size),np.sqrt(self.FOV_size)),fontsize=9)
                                    self.nax20.set_ylabel("Real pixels")
                                    self.nax20.set_xlabel("Spatial pixels/slices (%0.1f'')"%(Slitwidth))
                                    self.nax21.set_xlabel("Spatial pixels/slices (%0.1f'')"%(Slitwidth))
                                else:
                                    self.nax20.set_xlabel("Spatial pixels/fibers (%0.1f''Ø)"%(Slitwidth))
                                    self.nax20.set_ylabel("Spatial pixels/fibers (%0.1f''Ø)"%(Slitwidth))
                                    self.nax21.set_xlabel("Spatial pixels/fibers (%0.1f''Ø)"%(Slitwidth))
                                    self.nax20.set_title(r"%i × %i spaxels$^2$"%(n3,n3),fontsize=9)
                                    self.nax21.set_title("%0.1f' × %0.1f', full FOV= %0.1f' × %0.1f'"%(n3*pixel_scale/60,n3/pixel_scale/60,np.sqrt(self.FOV_size),np.sqrt(self.FOV_size)),fontsize=9)

                                # labels =  ['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]
                                labels = [r"➛$\mathbf{%s}$: %0.3f (%0.1f%%)" % (name, getattr(self.new, "electrons_per_pix")[self.new.i, j], 100 * getattr(self.new, "electrons_per_pix")[self.new.i, j] / np.sum(getattr(self.new, "electrons_per_pix")[self.new.i, :])) 
                                        if getattr(self.new, "electrons_per_pix")[self.new.i, j] == np.max(getattr(self.new, "electrons_per_pix")[self.new.i, :]) 
                                        else '%s: %0.3f (%0.1f%%)' % (name, getattr(self.new, "electrons_per_pix")[self.new.i, j], 100 * getattr(self.new, "electrons_per_pix")[self.new.i, j] / np.sum(getattr(self.new, "electrons_per_pix")[self.new.i, :])) 
                                        for j, name in enumerate(self.new.names)]
                                self.position2[0].set_label("  \n".join(labels)) 
                                self.nax21.legend(loc="upper left",handlelength=0, handletextpad=0, fancybox=True,markerscale=0,fontsize=8,title="Signal contribution (/pix)")#Overall background: %0.3f (%0.1f%%)"%(np.nansum(self.new.electrons_per_pix[self.new.i,1:]),100*np.nansum(self.new.electrons_per_pix[self.new.i,1:])/np.nansum(self.new.electrons_per_pix[self.new.i,:])))
                                
                                labels = [     r'%s: %0.2f (%0.1f%%)' % (name, self.new.noises[self.new.i, i] / self.new.factor[self.new.i], self.new.percents[i, self.new.i])     if self.new.percents[i, self.new.i] < np.max(self.new.percents[:, self.new.i])     else r"➛$\mathbf{%s}$: %0.2f (%0.1f%%)" % (name, self.new.noises[self.new.i, i] / self.new.factor[self.new.i], self.new.percents[i, self.new.i])     for i, name in enumerate(self.new.names)]
                                self.position1[0].set_label(r'{}'.format("\n".join(labels)))
                                self.nax20.legend(loc="upper left", handlelength=0, handletextpad=0, fancybox=True, markerscale=0, fontsize=8, title="%s Noise (/pix)"%(self.instrument.value))

                                # x1,x2 = (self.im.shape[0]-1)/2-(self.im.shape[0]/n2)*self.Size_source.value/pixel_scale/2.35 ,  (self.im.shape[0]-1)/2+(self.im.shape[0]/n2)*self.Size_source.value/pixel_scale/2.35
                                # y1, y2 = Δx + (n2-1)/2-self.Size_source.value/pixel_scale/2.35 , Δx + (n2-1)/2+self.Size_source.value/pixel_scale/2.35
                                x1,x2 = (self.im.shape[0]-1)/2-(self.im.shape[0]/n2)*self.Size_source.value/pixel_scale ,  (self.im.shape[0]-1)/2+(self.im.shape[0]/n2)*self.Size_source.value/pixel_scale
                                y1, y2 = Δx + (n2-1)/2-self.Size_source.value/pixel_scale , Δx + (n2-1)/2+self.Size_source.value/pixel_scale

                                self.position1[0].set_data([ (x1+x2)/2],[(y1+y2)/2]) 
                                self.position2[0].set_data([(x1+x2)/2],[(y1+y2)/2]) 

                                self.stack_square[0].set_data([x1,x2,x2,x1,x1],[y2,y2,y1,y1,y2])


                                
                                ifs_integ_spectra_stack =        np.nanmean(self.ifs_cube_stack[int(y1):int(y2),:,int(x1):int(x2)],axis=(0,2)) 
                                self.ifs_integ_spectra_stack[0].set_ydata( ifs_integ_spectra_stack)


                                im1 = np.nanmean(self.ifs_cube_stack[:,min(max(0,int(n1/2+0.5-lambda_stack/2)+int(Δλ)),n1-2)   :min(max(1,int(n1/2+0.5+lambda_stack/2)+int(Δλ)),n1-1),:],axis=1)
                                im2 = np.nanmean(self.ifs_cube[:,min(max(0,int(n1/2+0.5-lambda_stack/2)+int(Δλ)),n1-2)   :min(max(1,int(n1/2+0.5+lambda_stack/2)+int(Δλ)),n1-1),:],axis=1)
                                
                                
                                if log==False:
                                    self.ifs_slice_stack = self.nax21.imshow(im1, aspect="auto",cmap=self.current_cmap,interpolation=interpolation, vmin=np.nanmin(im1)+ self.vmin*np.ptp(im1),vmax=np.nanmax(im1)-  (1-self.vmax)*np.ptp(im1))#,vmin=np.nanmin(self.ifs_cube_stack),vmax=np.nanmax(self.ifs_cube_stack))
                                    self.ifs_slice       = self.nax20.imshow(im2, aspect="auto",cmap=self.current_cmap,interpolation=interpolation, vmin=np.nanmin(im1)+ self.vmin*np.ptp(im2),vmax=np.nanmax(im2)-  (1-self.vmax)*np.ptp(im2))
                                else:
                                    self.ifs_slice_stack = self.nax21.imshow(im1, aspect="auto",cmap=self.current_cmap,interpolation=interpolation,norm=SymLogNorm(linthresh=np.nanmax(im1)/10000, vmin=np.nanmin(im1)+ self.vmin*np.ptp(im1),vmax=np.nanmax(im1)-  (1-self.vmax)*np.ptp(im1)))
                                    self.ifs_slice       = self.nax20.imshow(im2, aspect="auto",cmap=self.current_cmap,interpolation=interpolation,norm=SymLogNorm(linthresh=np.nanmax(im2)/10000, vmin=np.nanmin(im1)+ self.vmin*np.ptp(im2),vmax=np.nanmax(im2)-  (1-self.vmax)*np.ptp(im2)))
                                    # self.ifs_slice_stack = self.nax21.imshow(im1, aspect="auto",cmap=self.current_cmap,interpolation=interpolation,norm=LogNorm( vmin=np.nanmin(im1)+ self.vmin*np.ptp(im1),vmax=np.nanmax(im1)-  (1-self.vmax)*np.ptp(im1)))
                                    # self.ifs_slice       = self.nax20.imshow(im2, aspect="auto",cmap=self.current_cmap,interpolation=interpolation,norm=LogNorm( vmin=np.nanmin(im1)+ self.vmin*np.ptp(im2),vmax=np.nanmax(im2)-  (1-self.vmax)*np.ptp(im2)))
                                self.cbar_slicer1 = self.fig3.colorbar(self.ifs_slice, cax=self.cax_slicer, orientation='vertical')
                                self.cbar_slicer2 = self.fig3.colorbar(self.ifs_slice_stack, cax=self.cax_slicer0, orientation='vertical')
                                if log==False:
                                    self.cbar_slicer1.formatter.set_powerlimits((0, 0))
                                    self.cbar_slicer2.formatter.set_powerlimits((0, 0))
                                    self.cbar_slicer1.formatter.set_useMathText(True)
                                    self.cbar_slicer2.formatter.set_useMathText(True)
                                    self.nax2_rp.set_yscale('linear')
                                    self.nax2s.set_yscale('linear')
                                else:
                                    self.nax2_rp.set_yscale('log')
                                    self.nax2_rp.set_ylim(ymin=0.001)
                                    self.nax2s.set_yscale('log')
                                #     self.cbar_slicer1.set_yscale('log')
                                #     self.cbar_slicer2.set_yscale('log')
                                # self.ifs_slice_stack = self.nax21.imshow(np.nanmean(self.ifs_cube_stack[:,0:2,:],axis=1), aspect="auto",cmap=self.current_cmap)
                                # print(min(max(0,int(n1/2+0.5-lambda_stack/2)+int(Δλ/dispersion)),n1-1)  ,min(max(0,int(n1/2+0.5+lambda_stack/2)+int(Δλ/dispersion)),n1-1))
                                if ~self.only_spectra:
                                    ratio = int(im2.shape[0]/im2.shape[1])+1
                                    # print(im2.shape[0],im2.shape[1],ratio)
                                    radial = np.hstack([im2[:, i:i+1].repeat(ratio, axis=1) for i in range(im2.shape[1])])
                                    # print(radial.shape)
                                    self.res1.set_xdata([2.355*PSF_RMS_det/2,2.355*PSF_RMS_det/2])
                                    self.res2.set_xdata([Slitwidth/2,Slitwidth/2])
                                    try : 
                                        rsurf, rmean, profile, EE, NewCenter, stddev = radial_profile_normalized(radial,center=(ratio*(x1+x2)/2, (y1+y2)/2),anisotrope=False,n=n,center_type="maximum",size=100)
                                        rsurf2, rmean2, profile2, EE, NewCenter, stddev = radial_profile_normalized( np.hstack([im1[:, i:i+1].repeat(ratio, axis=1) for i in range(im2.shape[1])]),center=NewCenter,anisotrope=False,n=n,center_type="User",size=100)
                                        self.rp[0].set_data(rmean2[:40],profile2[:40])
                                        if self.source_im.value != "Source":
                                            self.rp2[0].set_data(rmean[:40],profile[:40])
                                        else:
                                            self.rp2[0].set_data([np.nan,np.nan],[np.nan,np.nan])
                                        # gaus_profile = lambda  x, a, sigma : a ** 2 * np.exp(-np.square(x / sigma) / 2)
                                        gaus_profile = lambda x, a, sigma: a**2 * np.exp(- (x)**2 / (2 * sigma**2))
                                        popt, pcov = curve_fit(gaus_profile, rmean2[:40],profile2[:40], p0=[np.ptp(profile2[:40]), 2])
                                        self.rp_fit[0].set_data(np.linspace(0,np.max(rmean[:40]),100),gaus_profile(np.linspace(0,np.max(rmean[:40]),100),*popt))
                                        self.rp_fit[0].set_label("Gaussian fit: σ={:.1f}''".format(popt[1]*pixel_scale))#
                                        # self.diff[0].set_data(rmean2[:40],np.cumsum((profile2[:40]-gaus_profile(rmean2[:40],*popt))/ np.cumsum(gaus_profile(rmean2[:40],*popt))))
                                        
                                        self.diff[0].set_data(rmean2[:40],(np.cumsum(profile2[:40]*rmean2[:40]-rmean2[:40]*gaus_profile(rmean2[:40],*popt))/np.cumsum(rmean2[:40]*gaus_profile(rmean2[:40],*popt))))
                                    except ValueError:
                                        self.rp_fit[0].set_data([np.nan,np.nan],[np.nan,np.nan])
                                        self.diff[0].set_data([np.nan,np.nan],[np.nan,np.nan])
                                    

                                    self.nax2_rp.legend(fontsize=7,loc="upper right",title="Radial profiles",title_fontsize=8)
                                    # self.nax2_rp.plot(rmean2, profile2,":")
                                    # self.nax2_rp.set_ylim((np.nanmin(profile),np.nanmax(profile)))
                                    # self.nax2_rp.set_xlim((np.nanmin(rmean),np.nanmax(rmean)))
                                    self.nax2_rp.relim()                 # Recalculate limits based on current data
                                    self.nax2_rp.autoscale(enable=True, axis='both')
                                    self.nax2_rp2.relim()                 # Recalculate limits based on current data
                                    self.nax2_rp2.autoscale(enable=True, axis='both')
                                
 

                                self.wavelength_line1.set_xdata([int(n1/2)+int(Δλ)+0.5+(lambda_stack-1)/2,int(n1/2)+int(Δλ)+0.5+(lambda_stack-1)/2])
                                self.wavelength_line2.set_xdata([int(n1/2)+int(Δλ)+0.5-(lambda_stack-1)/2,int(n1/2)+int(Δλ)+0.5-(lambda_stack-1)/2])
                                self.wavelength_line1.set_label("λ = %0.1f-%0.1f nm"%(self.f(int(n1/2)+int(Δλ)+0.5-(lambda_stack-1)/2), self.f(int(n1/2)+int(Δλ)+0.5+(lambda_stack-1)/2)))
                                self.nax2s.legend(fontsize=7,loc="upper right")
                                self.ifs_slice.set_clim(vmin=np.nanmin(self.ifs_cube), vmax=np.nanmax(self.ifs_cube))
                                # self.nax2s.set_ylim((np.nanmin(list(self.im_stack[int(ax/2),int(n2/2),:]) +list(ifs_integ_spectra_stack)),np.nanmax(list(self.im_stack[int(ax/2),int(n2/2),:]) + list(ifs_integ_spectra_stack) )))
                                # self.nax2s.set_ylim((np.nanmin(list(self.im_stack[int(ax/2),int(n2/2),:]) ),np.nanmax(list(self.im_stack[int(ax/2),int(n2/2),:])  )))

                                title = '%s : FOV=%0.1famin$^2$, λ=%inm, Total Throughput=%i%%, Effective area=%0.1fcm$^2$, Platescale=%.1f,  PSF$_{x,λ}$=%0.1f, %0.1f pix, Npix = %i, t=%ih, F=%0.1E'%(self.instrument.value,self.FOV_size, self.wavelength.value, 100*self.Throughput.value*self.QE.value*self.Atmosphere.value, 100*100*self.Throughput.value*self.QE.value*self.Atmosphere.value*self.Collecting_area.value,  pixel_scale, 2.35*self.PSF_RMS_det/pixel_scale, 10*self.wavelength.value/self.Spectral_resolution.value/dispersion, self.number_pixels_used, self.acquisition_time.value, self.Signal.value)   
                                # self.fig3.suptitle("Spectral image of %s at %0.1f nm"%(self.instrument.value, wavelength), fontsize=10)
                                self.fig3.suptitle(title,y=0.97,fontsize=10)

                    # title = 'Instrument=%s, FOV=%0.1famin$^2$, λ=%inm, Throughput=%i%%, Atm=%i%%, Platescale=%.1f, area=%0.1fm$^2$'%(instrument,self.instruments[instrument][self.instruments["Charact."]=="FOV_size"][0], self.instruments[instrument][self.instruments["Charact."]=="wavelength"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Throughput"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Atmosphere"][0], self.instruments[instrument][self.instruments["Charact."]=="pixel_scale"][0], self.instruments[instrument][self.instruments["Charact."]=="Collecting_area"][0])


                                self.nax2s.relim()                 # Recalculate limits based on current data
                                self.nax2s.autoscale(enable=True, axis='both')


                                self.nax2s.set_xlim((200,300))
                                # self.nax2s.set_xlim((240,260))
                                self.fig3.tight_layout()
                                self.fig3.canvas.draw()
                            self.change.value=True

         
                    
    def update_instrument(self, instrument):
        # with self.output:
        with self.out1:
            # print(1)
            #HACK this was added here just because there is another issue with imagers
            self.change.value=False
            # print(1.2)
            self.spectro = False if np.isnan(self.instruments_dict[self.instrument.value]["dispersion"]) else True
            # print(1.3)
            self.spectrograph.value = self.spectro
            # print(1.4)
            self.spectra.options =self.spectra_options if self.spectrograph.value else ["Gaussian","10 kpc spiral galaxy","10 kpc spiral galaxy + CGM","30 kpc spiral galaxy","100 kpc spiral galaxy"]#,"10 kpc spiral galaxy + CGM + 0.0001 filament","10 kpc spiral galaxy + CGM + 0.001 filament",
            # print(2)
            if self.spectrograph.value:
                self.output_tabs.set_title(1, 'Spectral image')
                self.Redshift.disabled = False if "Rest-frame" in self.spectra.value else True
            else:
                self.output_tabs.set_title(1, 'Image')
                self.Redshift.disabled = True if "Gaussian" in self.spectra.value else False
            # print(3)

            self.follow_temp.layout.visibility = 'visible' if ("FIREBall-2" in instrument) else 'hidden'
            self.Line_width.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'
            self.lambda_stack.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'
            self.wavelength.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'

            # print(4)
            self.QElambda.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'
            # self.atmlambda.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'
            # self.sky_lines.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'
            self.Δλ.layout.visibility = 'visible' if self.spectrograph.value else 'hidden'

            # print(5)

            # self.spectra.layout.visibility = 'visible' if self.spectro else 'hidden'
            # self.spectra.value ="10 kpc spiral galaxy"
            keys = list(self.instruments_dict[instrument].keys())
            keys.remove("Signal")
            print(5)
            for key in keys:
                if ~np.isnan(self.instruments_dict[instrument][key]):
                    # print(key, self.instruments_dict[instrument][key])
                    # rsetattr(self, '%s.value'%(key), self.instruments_dict[instrument][key]) 
                    try:
                        rsetattr(self, '%s.value'%(key), self.instruments_dict[instrument][key]) 
                        # print(key, self.instruments_dict[instrument][key])
                    except AttributeError as e:
                        setattr( self, key, self.instruments_dict[instrument][key])
                        print(e)
                    # print(key, self.instruments_dict[instrument][key])
            self.IFS.value = True if ((self.dimensions==3) & self.spectrograph.value) else False
            self.counting_mode.layout.visibility = 'visible' if self.EM_gain.value > 1 else 'hidden'

            # print(6)
            self.atmlambda.disabled = False if ((self.Altitude<100) & self.spectrograph.value) else True
            self.sky_lines.disabled = False if (self.spectrograph.value &((self.Altitude<10)   |   (os.path.exists("../data/Instruments/%s/Sky_emission_lines.csv"%(instrument.upper().replace(" ","_"))) )    )) else True  
            if self.spectrograph.value:
                self.Throughput_FWHM.layout.visibility = 'visible' if ((self.QElambda.value) &   ~(os.path.exists("../data/Instruments/%s/Throughput.csv"%(instrument.upper().replace(" ","_"))) )      )  else 'hidden'
            else:
                self.Throughput_FWHM.layout.visibility = 'visible'
            if ~self.spectrograph.value:
                self.QElambda.value = True
                self.atmlambda.value = True if (self.Altitude<100) else False
                self.sky_lines.value = True if (self.Altitude<10) else False

            # print(7)
            self.fwhm.value = (self.PSF_RMS_mask, self.PSF_RMS_det )
            self.SlitDims.value = (self.Slitwidth, self.Slitlength) 
            self.exposure.value = (self.readout_time, self.instruments[instrument][self.instruments["Charact."]=="exposure_time"][0] )
            self.smearing.layout.visibility = 'visible' if (("FIREBall-2" in instrument)|("SCWI" in instrument)) & (self.counting_mode.value)    else 'hidden'
            self.temperature.layout.visibility = 'visible' if ("FIREBall-2" in instrument) &  (self.follow_temp.value)  else 'hidden'
            # title = 'Instrument=%s, FOV=%0.1famin$^2$, λ=%inm, Throughput=%i%%, Atm=%i%%, Platescale=%.1f, area=%0.1fm$^2$, PSF_{x,λ}=%0.1,%0.1pix '%(instrument,self.instruments[instrument][self.instruments["Charact."]=="FOV_size"][0], self.instruments[instrument][self.instruments["Charact."]=="wavelength"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Throughput"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Atmosphere"][0], self.instruments[instrument][self.instruments["Charact."]=="pixel_scale"][0], self.instruments[instrument][self.instruments["Charact."]=="Collecting_area"][0],    )
            title = '%s : FOV=%0.1famin$^2$, λ=%inm, Total Throughput=%i%%, Effective area=%0.1fcm$^2$, Platescale=%.1f,  PSF$_{x,λ}$=%0.1f, %0.1f pix, Npix = %i '%(self.instrument.value,self.FOV_size, self.wavelength.value, 100*self.Throughput.value*self.QE.value*self.Atmosphere.value,100*100* self.Throughput.value*self.QE.value*self.Atmosphere.value*self.Collecting_area.value,  self.pixel_scale.value, 2.35*self.PSF_RMS_det/self.pixel_scale.value, 10*self.wavelength.value/self.Spectral_resolution.value/self.dispersion.value, self.number_pixels_used    )
            # title = 'Instrument=%s, FOV=%0.1famin$^2$, λ=%inm, Throughput=%i%%, Atm=%i%%, Platescale=%.1f, area=%0.1fm$^2$'%(instrument,self.instruments[instrument][self.instruments["Charact."]=="FOV_size"][0], self.instruments[instrument][self.instruments["Charact."]=="wavelength"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Throughput"][0], 100*self.instruments[instrument][self.instruments["Charact."]=="Atmosphere"][0], self.instruments[instrument][self.instruments["Charact."]=="pixel_scale"][0], self.instruments[instrument][self.instruments["Charact."]=="Collecting_area"][0])
            self.ax0.set_title(title,y=0.97,fontsize=10)
            # print(8)
            self.change.value=True
            self.fig.tight_layout

            self.on_instrument_change()
            self.show_tab3() if self.IFS.value else self.hide_tab3()
            # print(9)
            self.update(x_axis=self.x_axis.value, counting_mode=self.counting_mode.value,Sky=self.Sky.value,acquisition_time=self.acquisition_time.value,Signal=self.Signal.value,EM_gain=self.EM_gain.value,RN=self.RN.value,CIC_charge=self.CIC_charge.value,Dark_current=self.Dark_current.value,exposure=self.exposure.value,smearing=self.smearing.value,temperature=self.temperature.value,follow_temp=self.follow_temp.value,fwhm=self.fwhm.value,QE=self.QE.value,extra_background=self.extra_background.value,log=self.ylog.value,SNR_res=self.SNR_res.value,
            xlog=self.xlog.value , Collecting_area=self.Collecting_area.value , pixel_scale=self.pixel_scale.value , Throughput=self.Throughput.value , Spectral_resolution=self.Spectral_resolution.value , SlitDims=self.SlitDims.value , dispersion=self.dispersion.value , Size_source=self.Size_source.value , Line_width=self.Line_width.value , wavelength=self.wavelength.value , Δλ=self.Δλ.value , Δx=self.Δx.value , Atmosphere=self.Atmosphere.value , pixel_size=self.pixel_size.value , cosmic_ray_loss_per_sec=self.cosmic_ray_loss_per_sec.value, lambda_stack = self.lambda_stack.value,change=self.change,         test=self.test,
            spectra=self.spectra.value,units=self.units.value,Throughput_FWHM=self.Throughput_FWHM.value, QElambda=self.QElambda.value, atmlambda=self.atmlambda.value, fraction_lya=self.fraction_lya.value, IFS=self.IFS.value,sky_lines=self.sky_lines.value, Redshift=self.Redshift.value, Bandwidth=self.Bandwidth)
            # print(10)
        return
    

def subtract_continuum(im,continuum_im):
    """
    Subtracts a linear continuum interpolated along the spectral axis (axis=1)
    from the datacubes self.ifs_cube_stack and self.ifs_cube.
    """
    # Get shape (ny, nwave, nx)
    ny, nwave, nx = im.shape

    # Step 1: compute average continuum images on the spectral edges
    im1 = np.nanmean(continuum_im[:, :10, :], axis=1)  # shape (ny, nx)
    im2 = np.nanmean(continuum_im[:, -10:, :], axis=1)  # shape (ny, nx)

    # Step 2: create linear interpolation along axis=1 (wavelength axis)
    # Generate normalized interpolation weights (from 0 to 1)
    weights = np.linspace(0, 1, nwave)[np.newaxis, :, np.newaxis]  # shape (1, nwave, 1)

    # Broadcast im1 and im2 for interpolation
    continuum = (1 - weights) * im1[:, np.newaxis, :] + weights * im2[:, np.newaxis, :]  # shape (ny, nwave, nx)

    # Step 3: subtract continuum from cubes
    ifs_im_contsub = im - continuum
    return ifs_im_contsub




class Observation:
    @initializer
    # def __init__(self, instrument="FIREBall-2 2023", Atmosphere=0.5, Throughput=0.13*0.9, exposure_time=50, counting_mode=False, Signal=1e-16, EM_gain=1400, RN=109, CIC_charge=0.005, Dark_current=0.08, Sky=10000, readout_time=1.5, extra_background = 0,acquisition_time = 2,smearing=0,i=25,plot_=False,temperature=-100,n=n,PSF_RMS_mask=5, PSF_RMS_det=8, QE = 0.45,cosmic_ray_loss_per_sec=0.005,Size_source=16,lambda_stack=1,Slitwidth=5,Bandwidth=200,Collecting_area=1,Δx=0,Δλ=0,pixel_scale=np.nan, Spectral_resolution=np.nan, dispersion=np.nan,Line_width=np.nan,wavelength=np.nan, pixel_size=np.nan,len_xaxis=50):#,photon_kept=0.7#, flight_background_damping = 0.9
    def __init__(self, instruments=None, instrument="FIREBall-2 2025", Atmosphere=None, Throughput=None, exposure_time=None, counting_mode=False, Signal=None, EM_gain=None, RN=None, CIC_charge=None, Dark_current=None, Sky=None, readout_time=None, extra_background = None,acquisition_time = None,smearing=None,i=33,plot_=False,n=n,PSF_RMS_mask=None, PSF_RMS_det=None, QE = None,cosmic_ray_loss_per_sec=None,Size_source=None,lambda_stack=1,Slitwidth=None,Bandwidth=None,Collecting_area=None,Δx=None,Δλ=None,pixel_scale=None, Spectral_resolution=None, dispersion=None,Line_width=None,wavelength=None, pixel_size=None,len_xaxis=50,Slitlength=None,IFS=None, Redshift=None, Throughput_FWHM=None, SNR_res="per pix",spectrograph=True,test=True):#,photon_kept=0.7#, flight_background_damping = 0.9
    # def __init__(self, instrument="FIREBall-2 2023", Atmosphere=0.5, Throughput=0.13, exposure_time=50, counting_mode=False, Signal=1e-17, EM_gain=1500, RN=40, CIC_charge=0.005, Dark_current=1, Sky=2e-18, readout_time=5, extra_background = 0.5,acquisition_time = 2,smearing=1.50,i=33,plot_=False,n=n,PSF_RMS_mask=2.5, PSF_RMS_det=3, QE = 0.4,cosmic_ray_loss_per_sec=0.005,Size_source=16,lambda_stack=0.21,Slitwidth=6,Bandwidth=160,Collecting_area=0.707,Δx=0,Δλ=0,pixel_scale=1.1, Spectral_resolution=1300, dispersion=0.21,Line_width=15,wavelength=200, pixel_size=13,len_xaxis=50,Slitlength=10):#,photon_kept=0.7#, flight_background_damping = 0.9
        """
        ETC calculator: computes the noise budget at the detector level based on instrument/detector parameters
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
        #####################################
        # Any of the parameter below can either be a float or an array allowing to check the evolution of the SNR 
        #####################################
                
        # self.Signal = Gaussian2D(amplitude=self.Signal,x_mean=0,y_mean=0,x_stddev=self.Size_source,y_stddev=self.Line_width,theta=0)(self.Δx,self.Δλ)
        # print("\ni",self.i,"\nAtmosphere",self.Atmosphere, "\nThroughput=",self.Throughput,"\nSky=",self.Sky, "\nacquisition_time=",self.acquisition_time,"\ncounting_mode=",self.counting_mode,"\nSignal=",self.Signal,"\nEM_gain=",self.EM_gain,"RN=",self.RN,"CIC_charge=",self.CIC_charge,"Dark_current=",self.Dark_current,"\nreadout_time=",self.readout_time,"\nsmearing=",self.smearing,"\nextra_background=",self.extra_background,"\nPSF_RMS_mask=",self.PSF_RMS_mask,"\nPSF_RMS_det=",self.PSF_RMS_det,"\nQE=",self.QE,"\ncosmic_ray_loss_per_sec=",self.cosmic_ray_loss_per_sec,"\nlambda_stack",self.lambda_stack,"\nSlitwidth",self.Slitwidth, "\nBandwidth",self.Bandwidth,"\nSize_source",self.Size_source,"\nCollecting_area",self.Collecting_area)
        # print("\Collecting_area",self.Collecting_area, "\nΔx=",self.Δx,"\nΔλ=",self.Δλ, "\napixel_scale=",self.pixel_scale,"\nSpectral_resolution=",self.Spectral_resolution,"\ndispersion=",self.dispersion,"\nLine_width=",self.Line_width,"wavelength=",self.wavelength,"pixel_size=",self.pixel_size)
      
      
        self.spectro = self.spectrograph#False if np.isnan(self.instruments_dict[self.instrument]["dispersion"]) else True
        # Simple hack to me able to use UV magnitudes (not used for the ETC)
        if np.max([self.Signal])>1:
            self.Signal = 10**(-(self.Signal-20.08)/2.5)*2.06*1E-16
        #TODO be sure we account for potentialfwhm_sigma_ratio ratio here

        #####################################
        # Adjust signal for circular aperture (fiber) geometry
        #####################################
        if type(self.Slitlength) == np.float64:
            if (self.Slitlength==self.Slitwidth):
                self.Signal *= np.pi/4 # ratio between fiber disk and square slit
        #####################################
        # If precision mode is on (currently not), apply spatial and spectral resolution dimming
        #####################################
        if self.precise: # TODO are we sure we should do that here?
            self.Signal *= (erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) )
            #convolve input flux by spectral resolution
            self.spectro_resolution_A = 10*self.wavelength/self.Spectral_resolution
            if self.spectro:
                self.Signal *= (erf(self.Line_width / (2 * np.sqrt(2) * self.spectro_resolution_A  )) )
            # print("Factor spatial and spectral",  (erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) ),   (erf(self.Line_width / (2 * np.sqrt(2) * 10*self.wavelength/self.Spectral_resolution)) ))
        #####################################
        # Compute fraction of signal passing through the slit (if slit spectro)
        #####################################
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
        
        # @mat do you think that's right?
        self.PSF_lambda_pix = 10*self.wavelength / self.Spectral_resolution / self.dispersion
        fwhm_sigma_ratio =2.35#1.0 # 2.355
        if self.spectro:
            source_spatial_pixels = np.maximum(1,np.minimum(np.sqrt(self.Size_source**2+self.PSF_RMS_mask**2) * fwhm_sigma_ratio / self.pixel_scale, self.Slitlength / self.pixel_scale))
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
        #####################################
        # Determine number of pixels used in SNR estimation
        #####################################
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

        #####################################
        # Compute noise sources: CIC, dark current, and effective area
        #####################################
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
        # self.source_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio/2,self.Slitlength)   *          (self.Size_source   if self.IFS else   np.minimum(self.Size_source*fwhm_sigma_ratio/2,self.Slitwidth)   )
        # self.slit_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio/2,self.Slitlength)   *          (np.maximum(self.Size_source/2,self.Slitwidth)   if self.IFS else   self.Slitwidth   )

        self.source_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength)   *   (self.Size_source   if self.IFS else   np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitwidth)   )
        self.slit_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength)   *   (np.maximum(self.Size_source,self.Slitwidth)   if self.IFS else   self.Slitwidth   )

        #####################################
        # Convert sky background to LU and then to electrons
        #####################################
        # @mat do you think that's right?
        if self.spectro: # previously was multiplying by self.nfibers *
            # mat's solution provides a local optimum in dispersion that I don't get with my solution!!!
            self.factor_CU2el_tot =     1*self.effective_area * self.arcsec2str * np.minimum(self.Line_width,self.Bandwidth) *  self.source_size_arcsec_after_slit  / self.pixels_total_source  
            # self.factor_CU2el_sky_tot = 1*self.effective_area * self.arcsec2str * np.minimum(self.Line_width,self.Bandwidth) *  self.slit_size_arcsec_after_slit    / self.pixels_total_source  #it works only for a line emission and we take the total sky flux over the same pixels
            self.factor_CU2el_sky_tot = 1*self.effective_area * self.arcsec2str * np.maximum(np.minimum(self.Line_width,self.Bandwidth),self.dispersion) *  self.slit_size_arcsec_after_slit    / self.pixels_total_source  #it works only for a line emission and we take the total sky flux over the same pixels
            self.factor_CU2el = self.factor_CU2el_tot
            self.factor_CU2el_sky = self.factor_CU2el_sky_tot

            # working almost for dispersion! But actually we also need to have an optimal with pixel_scale!!!
            # self.factor_CU2el = self.effective_area * self.arcsec2str  * np.minimum(self.Slitwidth, self.Size_source)  * np.minimum( self.dispersion, np.minimum(self.Line_width, self.Bandwidth)/2.35) * self.pixel_scale
            # # il faut ici que je rajoute une line width
            # # quand la ligne width devient importante 
            # sky_spectral_coverage =  np.minimum( self.dispersion, np.minimum(self.Line_width, self.Bandwidth)/2.35)#np.minimum(self.dispersion, self.Bandwidth)
            # self.factor_CU2el_sky = self.effective_area * self.arcsec2str  * self.Slitwidth * sky_spectral_coverage * self.pixel_scale


            # working for dispersion but adding pixels_scale and slit_width!!!
            # source_spatial_arcsec = self.source_size_arcsec_after_slit#
            # # source_spatial_arcsec = np.minimum(self.Slitwidth, self.Size_source)
            # source_spectral_angstrom = np.minimum(self.Line_width, self.Bandwidth)
            # spatial_per_pix = source_spatial_arcsec / self.pixel_scale
            # spectral_per_pix = source_spectral_angstrom / self.dispersion
            # # Combine spatiale + spectrale en quadrature :
            # effective_pix_size = np.sqrt(spatial_per_pix**2 + spectral_per_pix**2)
            # # Puis :
            # self.factor_CU2el_average = self.effective_area * self.arcsec2str * source_spatial_arcsec * source_spectral_angstrom / self.pixels_total_source #effective_pix_size
            # sky_spatial_arcsec = self.Slitwidth
            # sky_spatial_arcsec = self.slit_size_arcsec_after_slit
            # sky_spectral_angstrom =  self.Bandwidth  #/2.35            #np.minimum(self.Line_width, self.Bandwidth)
            # sky_spectral_angstrom =  np.minimum(self.Line_width, self.Bandwidth)
            # spatial_per_pix_sky = sky_spatial_arcsec / self.pixel_scale
            # spectral_per_pix_sky = sky_spectral_angstrom / self.dispersion
            # # Combine spatiale + spectrale en quadrature :
            # effective_pix_size_sky = np.sqrt(spatial_per_pix_sky**2 + spectral_per_pix_sky**2)
            # # Puis :
            # self.factor_CU2el_sky_average = self.effective_area * self.arcsec2str * sky_spatial_arcsec * sky_spectral_angstrom / self.pixels_total_source #effective_pix_size_sky
                        
            # difference = np.abs(self.factor_CU2el_tot - self.factor_CU2el_average)
            # ratio = difference / np.abs(self.factor_CU2el_tot)
            # # print(ratio, difference)
            # # if (difference > 0.1) | (ratio > 0.1):
            # #     print("Warning: difference or ratio between the two methods to compute the factor is too high: ", difference, ratio)
            # #     print("factor_CU2el_tot", self.factor_CU2el_tot, "factor_CU2el_average", self.factor_CU2el_average)            
            # if  self.test:
            #     self.factor_CU2el = self.factor_CU2el_tot
            #     self.factor_CU2el_sky = self.factor_CU2el_sky_tot
            # else:
            #     self.factor_CU2el = self.factor_CU2el_average
            #     self.factor_CU2el_sky = self.factor_CU2el_sky_average

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




        #####################################
        # Compute other noise contributions (read noise, extra background)
        #####################################


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

    def revew_for_mat(self,IFS):



        #####################################
        # Any of the parameter below can either be a float or an array allowing to check the evolution of the SNR 
        #####################################
                
        self.precise = False
        self.spectro = self.spectrograph#False if np.isnan(self.instruments_dict[self.instrument]["dispersion"]) else True
        #TODO be sure we account for potentialfwhm_sigma_ratio ratio here
        #convolve input flux by instrument PSF
        
        if self.precise: # TODO are we sure we should do that here?
            self.Signal *= (erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) )
            #convolve input flux by spectral resolution
            self.spectro_resolution_A = 10*self.wavelength/self.Spectral_resolution
            if self.spectro:
                self.Signal *= (erf(self.Line_width / (2 * np.sqrt(2) * self.spectro_resolution_A  )) )
            # print("Factor spatial and spectral",  (erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) ),   (erf(self.Line_width / (2 * np.sqrt(2) * 10*self.wavelength/self.Spectral_resolution)) ))



        #####################################
        # Adjust signal for circular aperture (fiber) geometry
        #####################################
        
        if type(self.Slitlength) == np.float64:
            if (self.Slitlength==self.Slitwidth):
                self.Signal *= np.pi/4 # ratio between fiber disk and square slit

        #####################################
        # If precision mode is on (currently not), apply spatial and spectral resolution dimming
        #####################################
        if self.precise:
            self.Signal *= erf(self.Size_source / (2 * np.sqrt(2) * self.PSF_RMS_det))
            self.spectro_resolution_A = 10 * self.wavelength / self.Spectral_resolution
            if self.spectro:
                self.Signal *= erf(self.Line_width / (2 * np.sqrt(2) * self.spectro_resolution_A))


        #####################################
        # Compute fraction of signal passing through the slit (if slit spectro)
        #####################################
        
        if ~np.isnan(self.Slitwidth).all(): #& (self.precise) # & (self.SNR_res!="per Source")
            # assess flux fraction going through slit
            # self.flux_fraction = ((1+erf(self.Slitwidth/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)     *    ((1+erf(self.Slitlength/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)
            self.flux_fraction =   ((1+erf(self.Slitlength/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)
            if    (~self.IFS): 
                self.flux_fraction *=((1+erf(self.Slitwidth/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1)  
        else:
            self.flux_fraction = 1
        self.flux_fraction_slit_applied = self.flux_fraction


        # compute size of the spectral PSF in pixels
        self.PSF_lambda_pix = 10*self.wavelength / self.Spectral_resolution / self.dispersion
        fwhm_sigma_ratio =2.35#1.0 # 2.355
        # compute the size of the source/res elem in pixels
        if self.spectro:
            source_spatial_pixels = np.maximum(1,np.minimum(np.sqrt(self.Size_source**2+self.PSF_RMS_mask**2) * fwhm_sigma_ratio / self.pixel_scale, self.Slitlength / self.pixel_scale))
            source_spectral_pixels = np.maximum(1, np.sqrt((self.Slitwidth/self.pixel_scale)**2 +self.PSF_lambda_pix**2 + (np.minimum(self.Line_width, self.Bandwidth) / self.dispersion)**2))
            self.source_size = np.maximum(np.minimum(self.Size_source * fwhm_sigma_ratio, self.Slitlength) / self.pixel_scale,1)    * np.sqrt(source_spatial_pixels**2 + source_spectral_pixels**2)
            self.pixels_total_source =  self.source_size  * ( np.ceil(np.sqrt(self.Size_source**2+self.PSF_RMS_mask**2)*fwhm_sigma_ratio / self.Slitwidth) if self.IFS else 1)
        else:
            self.source_size =  (self.Size_source *fwhm_sigma_ratio /self.pixel_scale) **2
            self.elem_size = (self.PSF_RMS_det *fwhm_sigma_ratio /self.pixel_scale) **2

        #####################################
        # Determine number of pixels used in SNR estimation
        #####################################

        if self.SNR_res=="per pix":
          self.number_pixels_used = 1
        elif self.SNR_res=="per Res elem": # is that true ? when not IFS, the SNR won't get bigger than the slit , the rest will be cut
            self.number_pixels_used = np.ceil(self.elem_size)
        elif self.SNR_res=="per Source":
            self.number_pixels_used = np.ceil(self.pixels_total_source)

        #####################################
        # Compute noise sources: CIC, dark current, and effective area
        #####################################
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
                
        # The faction of detector lost by cosmic ray masking (taking into account ~5-10 impact per seconds and around 2000 pixels loss per impact (0.01%))
        self.cosmic_ray_loss = np.minimum(self.cosmic_ray_loss_per_sec*(self.exposure_time+self.readout_time/2),1)
        self.QE_efficiency = self.Photon_fraction_kept * self.QE
        self.effective_area =  self.QE_efficiency * self.Throughput * self.Atmosphere  *    (self.Collecting_area * 100 * 100)
        # TODO need to verify that the evolution the sky and signal makes sense with nfibers... Maybe this has been solved

        self.source_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength)   *   (self.Size_source   if self.IFS else   np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitwidth)   )
        self.slit_size_arcsec_after_slit =     np.minimum(self.Size_source*fwhm_sigma_ratio,self.Slitlength)   *   (np.maximum(self.Size_source,self.Slitwidth)   if self.IFS else   self.Slitwidth   )

        #####################################
        # Convert sky background to LU and then to electrons
        #####################################

        if self.spectro: # previously was multiplying by self.nfibers *
            # mat's solution provides a local optimum in dispersion that I don't get with my solution!!!
            self.factor_CU2el_tot =     1*self.effective_area * self.arcsec2str * np.minimum(self.Line_width,self.Bandwidth) *  self.source_size_arcsec_after_slit  / self.pixels_total_source  
            self.factor_CU2el_sky_tot = 1*self.effective_area * self.arcsec2str * np.maximum(np.minimum(self.Line_width,self.Bandwidth),self.dispersion) *  self.slit_size_arcsec_after_slit    / self.pixels_total_source  #it works only for a line emission and we take the total sky flux over the same pixels
            # TODO maybe missing a factor 2.35 here...

            # print(ratio, difference)
            # if (difference > 0.1) | (ratio > 0.1):
            #     print("Warning: difference or ratio between the two methods to compute the factor is too high: ", difference, ratio)
            #     print("factor_CU2el_tot", self.factor_CU2el_tot, "factor_CU2el_average", self.factor_CU2el_average)            
            self.factor_CU2el = self.factor_CU2el_tot
            self.factor_CU2el_sky = self.factor_CU2el_sky_tot


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




        #####################################
        # Compute other noise contributions (read noise, extra background)
        #####################################

        # TODO in counting mode, Photon_fraction_kept should also be used for CIC
        self.RN_final = self.RN  * self.RN_fraction_kept / self.EM_gain 
        self.Additional_background = self.extra_background/3600 * self.exposure_time# e/pix/exp
        self.Additional_background_noise = np.sqrt(self.Additional_background * self.ENF)
        
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

        if np.ndim(self.noises)==2:
            self.percents =  100* np.array(self.noises).T[:-1,:]**2/self.Total_noise_final**2
        else:
            self.percents =  100* np.array(self.noises).T[:-1]**2/self.Total_noise_final**2            
        
        self.el_per_pix = self.Signal_el + self.sky + self.CIC_charge +  self.Dark_current_f
        
        
        # NO NEED TO REVIEW AFTER THAT
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
        self.point_source_5s = self.extended_source_5s * 1.30e57
        self.time2reach_n_sigma_SNR = self.acquisition_time *  np.square(n_sigma / self.SNR)
       

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
                fraction_signal = interp(self.EM_gain/self.RN,a["G/RN"],a["fractionflux"])
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
    
                    
                    spectra =  flux *Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2 + (w_nm-self.wavelength)*10/self.dispersion, PSF_λ)   / Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2, self.PSF_lambda_pix**2/(PSF_λ**2 + self.PSF_lambda_pix**2)).sum()
                    # print(w_nm,self.wavelength,self.dispersion, size[0]/2 + (w_nm-self.wavelength)*10*self.dispersion,np.min(spectra),np.max(spectra))
                    # if spectra.sum()>0:
                    #     spectra *= flux/spectra.sum()
                    # print(spectra)
                    # else:
                    #     spectra = np.ones(size[0])
                    # / Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2 + (w_nm-self.wavelength)*10*self.dispersion, self.PSF_lambda_pix**2/(PSF_λ**2 + self.PSF_lambda_pix**2)).sum()
                    spectra *= atm_qe_normalized_shape   
                    source_im =  np.outer(spectra,spatial_profile*slit_profile ).T /Gaussian1D.evaluate(np.arange(size[1]),  1,  50, Rx**2/(PSF_x**2+Rx**2)).sum()
                    # source_im = np.outer(spectra, spatial_profile).T
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
                        # print("../data/Emission_cube/"+ source.split("-")[0] + ".fits")
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
                    # cube_detector = create_fits_cube(cube_size=(size[1], size[1], n_wave), pixel_scale=self.pixel_scale, wave_range=(wave_min,wave_max),continuum_flux=spectra* int(self.exposure_time),continuum_fwhm=PSF_x, line_flux=self.Signal_el*0, line_center=self.wavelength, line_fwhm=self.Line_width, line_spatial_fwhm=self.Size_source)
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
            # print("2D array: Removing NAXIS3 from header...")
            fitsimage.header.remove("NAXIS3")
        if "SKEW" in fitsimage.header:
            fitsimage.header.remove("SKEW")
    elif hasattr(fitsimage[0], "header"):
        if "NAXIS3" in fitsimage[0].header:
            # print("2D array: Removing NAXIS3 from header...")
            fitsimage[0].header.remove("NAXIS3")
        if "SKEW" in fitsimage[0].header:
            fitsimage[0].header.remove("SKEW")

    elif not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    if not os.path.exists(os.path.dirname(filename)):
        # print("%s not existing: creating folde..." % (os.path.dirname(filename)))
        os.makedirs(os.path.dirname(filename))
    try:
        fitsimage.writeto(filename, overwrite=True)
    except IOError:
        # print("Can not write in this repository : " + filename)
        filename = "%s/%s" % (
            os.path.dirname(os.path.dirname(filename)),
            os.path.basename(filename),
        )
        # filename = "/tmp/" + os.path.basename(filename)
        # print("Instead writing new file in : " + filename) 
        fitsimage.writeto(filename, overwrite=True)
    # print("Image saved: %s" % (filename))
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


def estimate_background(data, center, radius=30, n=1.8):
    """Function that estimate the Background behing a source given an inner
    radius and a factor n to the outter radius such as the background is
    computed on the area which is on C2(n*radius)/C1(radius)
    """
    import numpy as np

    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    r = r.astype(int)
    mask = (r >= radius) & (r <= n * radius)
    fond = np.nanmean(data[mask])
    return fond


def radial_profile_normalized(
    data,
    center,
    anisotrope=False,
    angle=30,
    radius=40,
    n=1.5,
    center_type="barycentre",
    radius_bg=70,
    n1=20,
    size=70,
):
    """Function that returns the radial profile of a spot
    given an input image + center.
    """
    from scipy import ndimage
    import numpy as np

    y, x = np.indices((data.shape))
    # print("center_type = %s" % (center_type))
    n1 = np.nanmin([n1, int(center[1]), int(center[0])])
    image = data[
        int(center[1]) - n1 : int(center[1]) + n1,
        int(center[0]) - n1 : int(center[0]) + n1,
    ]
    if center_type.lower() == "maximum":
        barycentre = np.array(
            [
                np.where(image == image.max())[0][0],
                np.where(image == image.max())[1][0],
            ]
        )
    if center_type.lower() == "barycentre":
        background = estimate_background(data, center, radius, 1.8)
        new_image = image - background
        index = new_image > 0.5 * np.nanmax(new_image)
        new_image[~index] = 0
        barycentre = ndimage.measurements.center_of_mass(new_image)
    if center_type.lower() == "user":
        barycentre = [n1, n1]
    else:
        # print("Center type not understood, taking barycenter one")
        background = estimate_background(data, center, radius, 1.8)
        new_image = image - background
        index = new_image > 0.5 * np.nanmax(new_image)  # .max()
        new_image[~index] = 0
        barycentre = ndimage.measurements.center_of_mass(
            new_image
        )  # background#np.nanmin(image)
    new_center = np.array(center) + barycentre[::-1] - n1
    # print(
    #     "new_center = {}, defined with center type: {}".format(new_center, center_type)
    # )
    if radius_bg:
        fond = estimate_background(data, new_center, radius, n)
    else:
        fond = 0
    image = data - fond  # (data - fond).astype(int)

    r = np.sqrt((x - new_center[0]) ** 2 + (y - new_center[1]) ** 2)
    #    r = np.around(r)-1
    rint = r.astype(int)

    image_normalized = image
    # / np.nansum(image[r<radius])
    if anisotrope == True:
        theta = abs(180 * np.arctan((y - new_center[1]) / (x - new_center[0])) / np.pi)
        tbin_spectral = np.bincount(
            r[theta < angle].ravel(), image_normalized[theta < angle].ravel()
        )
        tbin_spatial = np.bincount(
            r[theta > 90 - angle].ravel(), image_normalized[theta > 90 - angle].ravel(),
        )
        nr_spectral = np.bincount(r[theta < angle].ravel())
        nr_spatial = np.bincount(r[theta > 90 - angle].ravel())
        EE_spatial = (
            100
            * np.nancumsum(tbin_spatial)
            / np.nanmax(np.nancumsum(tbin_spatial)[:100] + 1e-5)
        )
        EE_spectral = (
            100
            * np.nancumsum(tbin_spectral)
            / np.nanmax(np.nancumsum(tbin_spectral)[:100] + 1e-5)
        )
        return (
            tbin_spectral / nr_spectral,
            tbin_spatial / nr_spatial,
            EE_spectral,
            EE_spatial,
        )
    else:
        tbin = np.bincount(rint.ravel(), image_normalized.ravel())
        nr = np.bincount(rint.ravel())
        rsurf = np.sqrt(np.nancumsum(nr) / np.pi)
        rmean = np.bincount(rint.ravel(), r.ravel()) / nr
        dist = np.array(rint[rint < radius].ravel(), dtype=int)
        data = image[rint < radius].ravel()
        stdd = [
            np.nanstd(data[dist == distance]) / np.sqrt(len(data[dist == distance]))
            for distance in np.arange(size)
        ]
        radialprofile = tbin / nr
        EE = np.nancumsum(tbin) * 100 / np.nanmax(np.nancumsum(tbin)[:radius] + 1e-5)
        return (
            rsurf[:size],
            rmean[:size],
            radialprofile[:size],
            EE[:size],
            new_center[:size],
            stdd[:size],
        )

# %load_ext line_profiler
# %lprun -f Observation  Observation()#.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Slit", size=[100, 100], OSregions=[0, 100], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=True,atmlambda=True,fraction_lya=0.05)

## %%
# %lprun -u 1e-1 -T /tmp/initilize.py -s -r -f  Observation.initilize  Observation(exposure_time=np.linspace(50,1500,50))
# %lprun -u 1e-1 -T /tmp/interpolate_optimal_threshold.py -s -r -f  Observation.interpolate_optimal_threshold  Observation(exposure_time=np.linspace(50,1500,50),counting_mode=True,plot_=False).interpolate_optimal_threshold()
# %lprun -u 1e-1 -T /tmp/PlotNoise.py -s -r -f  Observation.PlotNoise  Observation(exposure_time=np.linspace(50,1500,50)).PlotNoise()
# %lprun -u 1e-1 -T /tmp/SimulateFIREBallemCCDImage.py -s -r -f  Observation.SimulateFIREBallemCCDImage  Observation(exposure_time=np.linspace(50,1500,50)).SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Slit", size=[100, 100], OSregions=[0, 100], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=True,atmlambda=True,fraction_lya=0.05)
## %%


# %%
