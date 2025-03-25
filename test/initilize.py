Timer unit: 0.1 s

Total time: 0.005715 s
File: /Users/Vincent/Github/fireball2-etc/notebooks/Observation.py
Function: initilize at line 140

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   140                                               def initilize(self):
   141         1          0.0      0.0      0.1          self.precise = True
   142         1          0.0      0.0     50.0          self.Signal = Gaussian2D(amplitude=self.Signal,x_mean=0,y_mean=0,x_stddev=self.PSF_source,y_stddev=4,theta=0)(self.Δx,self.Δλ)
   143                                                   # print("\nAtmosphere",self.Atmosphere, "\nThroughput=",self.Throughput,"\nSky=",Sky, "\nacquisition_time=",acquisition_time,"\ncounting_mode=",counting_mode,"\nSignal=",Signal,"\nEM_gain=",EM_gain,"RN=",RN,"CIC_charge=",CIC_charge,"Dard_current=",Dard_current,"\nreadout_time=",readout_time,"\n_smearing=",smearing,"\nextra_background=",extra_background,"\ntemperature=",temperature,"\nPSF_RMS_mask=",PSF_RMS_mask,"\nPSF_RMS_det=",PSF_RMS_det,"\nQE=",QE,"\ncosmic_ray_loss_per_sec=",self.cosmic_ray_loss_per_sec,"\nlambda_stack",self.lambda_stack,"\nSlitwidth",self.Slitwidth, "\nBandwidth",self.Bandwidth,"\nPSF_source",self.PSF_source,"\nCollecting_area",self.Collecting_area)
   144                                                   # print("\Collecting_area",self.Collecting_area, "\nΔx=",self.Δx,"\nΔλ=",Δλ, "\napixel_scale=",pixel_scale,"\nSpectral_resolution=",Spectral_resolution,"\ndispersion=",dispersion,"\nLine_width=",Line_width,"wavelength=",wavelength,"pixel_size=",pixel_size)
   145                                                   
   146                                                   # Simple hack to me able to use UV magnitudes (not used for the ETC)
   147         1          0.0      0.0      0.8          if np.max([self.Signal])>1:
   148                                                       self.Signal = 10**(-(Signal-20.08)/2.5)*2.06*1E-16
   149                                                   #TODO be sure we account for potential 2.35 ratio here
   150                                                   #convolve input flux by instrument PSF
   151         1          0.0      0.0      0.1          if self.precise:
   152         1          0.0      0.0     19.8              self.Signal *= (erf(self.PSF_source / (2 * np.sqrt(2) * self.PSF_RMS_det)) )
   153                                                       #convolve input flux by spectral resolution
   154                                                       # self.spectro_resolution_A = self.wavelength * self.spectral
   155         1          0.0      0.0      0.2              self.Signal *= (erf(self.Line_width / (2 * np.sqrt(2) * 10*self.wavelength/self.Spectral_resolution)) )
   156                                           
   157                                           
   158         1          0.0      0.0      3.7          if ~np.isnan(self.Slitwidth).all() & self.precise:
   159                                                       # assess flux fraction going through slit
   160         1          0.0      0.0      0.1              self.flux_fraction_slit = (1+erf(self.Slitwidth/(2*np.sqrt(2)*self.PSF_RMS_mask)))-1
   161                                                   else:
   162                                                       self.flux_fraction_slit = 1
   163                                                   
   164         1          0.0      0.0      0.1          self.resolution_element= self.PSF_RMS_det * 2.35 /self.pixel_scale  # in pix (before it was in arcseconds)
   165         1          0.0      0.0      0.0          self.PSF_lambda_pix = self.wavelength / self.Spectral_resolution / self.dispersion
   166                                           
   167         1          0.0      0.0      0.0          red, blue, violet, yellow, green, pink, grey  = '#E24A33','#348ABD','#988ED5','#FBC15E','#8EBA42','#FFB5B8','#777777'
   168                                                   # self.colors= ['#E24A33','#348ABD','#988ED5','#FBC15E','#FFB5B8','#8EBA42','#777777']
   169                                                   # self.colors= ['#E24A33','#348ABD','#988ED5','#FBC15E','#8EBA42','#FFB5B8','#777777']
   170         1          0.0      0.0      0.1          self.colors= [red, violet, yellow  ,blue, green, pink, grey ]
   171                                                   # self.Sky_CU =  convert_ergs2LU(self.Sky_,self.wavelength,self.pixel_scale)
   172                                                   # self.Sky_ = self.Sky_CU*self.lu2ergs# ergs/cm2/s/arcsec^2 
   173                                           
   174         1          0.0      0.0      0.0          self.ENF = 1 if self.counting_mode else 2 # Excess Noise Factor 
   175         1          0.0      0.0      0.1          self.CIC_noise = np.sqrt(self.CIC_charge * self.ENF) 
   176         1          0.0      0.0      0.2          self.Dark_current_f = self.Dard_current * self.exposure_time / 3600 # e/pix/frame
   177         1          0.0      0.0      0.1          self.Dark_current_noise =  np.sqrt(self.Dark_current_f * self.ENF)
   178                                                   
   179                                                   # For now we put the regular QE without taking into account the photon kept fracton, because then infinite loop. 
   180                                                   # Two methods to compute it: interpolate_optimal_threshold & compute_optimal_threshold
   181         1          0.0      0.0      0.1          self.pixel_size_arcsec = self.pixel_scale
   182                                                   # self.pixel_scale  = (self.pixel_scale*np.pi/180/3600) #go from arcsec/pix to str/pix 
   183         1          0.0      0.0      0.1          self.arcsec2str = (np.pi/180/3600)**2
   184         1          0.0      0.0      0.2          self.Sky_CU = convert_ergs2LU(self.Sky, self.wavelength,self.pixel_size_arcsec) 
   185                                                   # self.Sky_ = convert_LU2ergs(self.Sky_CU, self.wavelength,self.pixel_size_arcsec) 
   186                                                   # self.Collecting_area *= 100 * 100#m2 to cm2
   187                                                   # TODO use astropy.unit
   188         1          0.0      0.0      0.0          if self.counting_mode:
   189                                                       self.factor_CU2el =  self.QE * self.Throughput * self.Atmosphere  *    (self.Collecting_area * 100 * 100)  * self.Slitwidth * self.arcsec2str  * self.dispersion
   190                                                       self.sky = self.Sky_CU*self.factor_CU2el*self.exposure_time  # el/pix/frame
   191                                                       self.Sky_noise_pre_thresholding = np.sqrt(self.sky * self.ENF) 
   192                                                       self.signal_pre_thresholding = self.Signal*self.factor_CU2el*self.exposure_time  # el/pix/frame
   193                                                       self.n_threshold, self.Photon_fraction_kept, self.RN_fraction_kept, self.gain_thresholding = self.interpolate_optimal_threshold(plot_=self.plot_, i=self.i)#,flux=self.signal_pre_thresholding)
   194                                                       # self.n_threshold, self.Photon_fraction_kept, self.RN_fraction_kept, self.gain_thresholding = self.compute_optimal_threshold(plot_=plot_, i=i,flux=self.signal_pre_thresholding)
   195                                                   else:
   196         1          0.0      0.0      5.2              self.n_threshold, self.Photon_fraction_kept, self.RN_fraction_kept, self.gain_thresholding = np.zeros(self.len_xaxis),np.ones(self.len_xaxis),np.ones(self.len_xaxis), np.zeros(self.len_xaxis)
   197                                                   # The faction of detector lost by cosmic ray masking (taking into account ~5-10 impact per seconds and around 2000 pixels loss per impact (0.01%))
   198         1          0.0      0.0      3.3          self.cosmic_ray_loss = np.minimum(self.cosmic_ray_loss_per_sec*(self.exposure_time+self.readout_time/2),1)
   199         1          0.0      0.0      0.1          self.QE_efficiency = self.Photon_fraction_kept * self.QE
   200                                                   # TODO verify that indeed it should not depend on self.pixel_scale**2 
   201                                                   # Compute ratio to convert CU to el/pix 
   202         1          0.0      0.0      0.1          if np.isnan(self.Slitwidth).all():
   203                                                       # If instrument is not a spectro?
   204                                                       self.factor_CU2el = self.QE_efficiency * self.Throughput * self.Atmosphere  *    (self.Collecting_area * 100 * 100)   * self.Bandwidth  * self.arcsec2str *self.pixel_scale**2 #but here it's total number of electrons we don't know if it is per A or not and so if we need to devide by dispersion: 1LU/A = .. /A. OK So we need to know if sky is LU or LU/A            
   205                                                   else:
   206         1          0.0      0.0      0.4              self.factor_CU2el = self.QE_efficiency * self.Throughput * self.Atmosphere  *    (self.Collecting_area * 100 * 100)  * self.Slitwidth * self.arcsec2str  * self.dispersion *self.pixel_scale**2
   207                                                   
   208                                           
   209         1          0.0      0.0      0.1          self.sky = self.Sky_CU*self.factor_CU2el*self.exposure_time  # el/pix/frame
   210         1          0.0      0.0      0.1          self.Sky_noise = np.sqrt(self.sky * self.ENF) 
   211                                                       
   212                                                   # TODO in counting mode, Photon_fraction_kept should also be used for CIC
   213         1          0.0      0.0      0.1          self.RN_final = self.RN  * self.RN_fraction_kept / self.EM_gain 
   214         1          0.0      0.0      0.1          self.Additional_background = self.extra_background/3600 * self.exposure_time# e/pix/exp
   215         1          0.0      0.0      0.1          self.Additional_background_noise = np.sqrt(self.Additional_background * self.ENF)
   216                                                   
   217                                                   # number of images taken during one field acquisition (~2h)
   218         1          0.0      0.0      0.3          self.N_images = self.acquisition_time*3600/(self.exposure_time + self.readout_time)
   219         1          0.0      0.0      0.1          self.N_images_true = self.N_images * (1-self.cosmic_ray_loss)
   220                                           
   221         1          0.0      0.0      0.1          self.Signal_LU = convert_ergs2LU(self.Signal,self.wavelength,self.pixel_size_arcsec)
   222                                                   # if 1==0: # if line is totally resolved (for cosmic web for instance)
   223                                                   #     self.Signal_el =  self.Signal_LU*self.factor_CU2el*self.exposure_time * self.flux_fraction_slit  / self.spectral_resolution_pixel # el/pix/frame#     Signal * (sky / Sky_)  #el/pix
   224                                                   # else: # if line is unresolved for QSO for instance
   225         1          0.0      0.0      0.8          self.Signal_el =  self.Signal_LU * self.factor_CU2el * self.exposure_time * self.flux_fraction_slit   # el/pix/frame#     Signal * (sky / Sky_)  #el/pix
   226                                                   # print(self.flux_fraction_slit)
   227                                           
   228         1          0.0      0.0      0.1          self.signal_noise = np.sqrt(self.Signal_el * self.ENF)     #el / resol/ N frame
   229                                           
   230         1          0.0      0.0      0.0          self.N_resol_element_A = self.lambda_stack / self.dispersion# / (1/self.dispersion)#/ (10*self.wavelength/self.Spectral_resolution) # should work even when no spectral resolution
   231         1          0.0      0.0      0.1          self.factor = np.sqrt(self.N_images_true) * self.resolution_element * np.sqrt(self.N_resol_element_A)
   232         1          0.0      0.0      0.1          self.Signal_resolution = self.Signal_el * self.factor**2# el/N exposure/resol
   233         1          0.0      0.0      0.1          self.signal_noise_nframe = self.signal_noise * self.factor
   234         1          0.0      0.0      0.2          self.Total_noise_final = self.factor*np.sqrt(self.signal_noise**2 + self.Dark_current_noise**2  + self.Additional_background_noise**2 + self.Sky_noise**2 + self.CIC_noise**2 + self.RN_final**2   ) #e/  pix/frame
   235         1          0.0      0.0      0.1          self.SNR = self.Signal_resolution / self.Total_noise_final
   236                                                   
   237         1          0.0      0.0      0.1          if type(self.Total_noise_final + self.Signal_resolution) == np.float64:
   238                                                       n=0
   239                                                   else:
   240         1          0.0      0.0      0.1              n =len(self.Total_noise_final + self.Signal_resolution) 
   241         1          0.0      0.0      0.0          if n>1:
   242        14          0.0      0.0      0.5              for name in ["signal_noise","Dark_current_noise", "Additional_background_noise","Sky_noise", "CIC_noise", "RN_final","Signal_resolution","Signal_el","sky","CIC_charge","Dark_current_f","RN","Additional_background"]:
   243        13          0.0      0.0      1.5                  setattr(self, name, getattr(self,name)*np.ones(n))
   244         1          0.0      0.0      0.1          self.factor = self.factor*np.ones(n) if type(self.factor)== np.float64 else self.factor
   245         1          0.0      0.0      2.7          self.noises = np.array([self.signal_noise*self.factor,  self.Dark_current_noise*self.factor,  self.Sky_noise*self.factor, self.RN_final*self.factor, self.CIC_noise*self.factor, self.Additional_background_noise*self.factor, self.Signal_resolution]).T
   246         1          0.0      0.0      0.2          self.electrons_per_pix =  np.array([self.Signal_el,  self.Dark_current_f,  self.sky,  self.RN_final, self.CIC_charge, self.Additional_background]).T
   247         1          0.0      0.0      0.1          self.names = ["Signal","Dark current", "Sky", "Read noise","CIC", "Extra background"]
   248         1          0.0      0.0      0.1          self.snrs=self.Signal_resolution /self.Total_noise_final
   249                                           
   250         1          0.0      0.0      0.1          if np.ndim(self.noises)==2:
   251         1          0.0      0.0      1.1              self.percents =  100* np.array(self.noises).T[:-1,:]**2/self.Total_noise_final**2
   252                                                   else:
   253                                                       self.percents =  100* np.array(self.noises).T[:-1]**2/self.Total_noise_final**2            
   254                                                   
   255         1          0.0      0.0      0.1          self.el_per_pix = self.Signal_el + self.sky + self.CIC_charge +  self.Dark_current_f
   256         1          0.0      0.0      0.1          n_sigma = 5
   257         1          0.0      0.0      0.3          self.signal_nsig_e_resol_nframe = (n_sigma**2 * self.ENF + n_sigma**2 * np.sqrt(4*self.Total_noise_final**2 - 4*self.signal_noise_nframe**2 + self.ENF**2*n_sigma**2))/2
   258         1          0.0      0.0      5.2          self.eresolnframe2lu = self.Signal_LU/self.Signal_resolution
   259         1          0.0      0.0      0.1          self.signal_nsig_LU = self.signal_nsig_e_resol_nframe * self.eresolnframe2lu
   260         1          0.0      0.0      0.2          self.signal_nsig_ergs = convert_LU2ergs(self.signal_nsig_LU, self.wavelength,self.pixel_size_arcsec) # self.signal_nsig_LU * self.lu2ergs
   261         1          0.0      0.0      0.1          self.extended_source_5s = self.signal_nsig_ergs * (self.pixel_scale*self.PSF_RMS_det)**2
   262         1          0.0      0.0      0.1          self.point_source_5s = self.extended_source_5s * 1.30e57
   263         1          0.0      0.0      0.3          self.time2reach_n_sigma_SNR = self.acquisition_time *  np.square(n_sigma / self.snrs)
   264                                                   # print("factor=",self.factor[self.i])
   265                                                   # print("N_images_true=",np.sqrt(self.N_images_true)[self.i] )
   266                                                   # print("resolution_element=", self.resolution_element)
   267                                                   # print("N_resol_element_A=",np.sqrt(self.N_resol_element_A))
   268                                                   # print("lambda_stack=",self.lambda_stack)
   269                                                   # print("dispersion=",self.dispersion)
   270                                                   # print("cosmic_ray_loss=",np.sqrt(self.cosmic_ray_loss)[self.i])
   271                                                   # print("N_images=",np.sqrt(self.N_images)[self.i])
   272                                           
   273                                                   #TODO change this ratio of 1.30e57
   274                                                   # from astropy.cosmology import Planck15 as cosmo
   275                                                   # 4*np.pi* (cosmo.luminosity_distance(z=0.7).to("cm").value)**2 = 2.30e57