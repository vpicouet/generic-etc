
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


        # compute size of the spectral PSF in pixels # add without 
        self.PSF_lambda_pix = 10*self.wavelength / self.Spectral_resolution / self.dispersion # usually based on FWHM
        fwhm_sigma_ratio =2.35#1.0 # 2.355
        # compute the size of the source/res elem in pixels
        if self.spectro:
            source_spatial_pixels = np.maximum(1,np.minimum(np.sqrt(self.Size_source**2+self.PSF_RMS_det**2) * fwhm_sigma_ratio / self.pixel_scale, self.Slitlength / self.pixel_scale))
            if self.Slitwidth>>self.PSF_RMS_mask:
                source_spectral_pixels = np.maximum(1, np.sqrt((self.Slitwidth/self.pixel_scale/4)**2 +self.PSF_lambda_pix**2 + (np.minimum(self.Line_width, self.Bandwidth) / self.dispersion)**2))
            else:
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