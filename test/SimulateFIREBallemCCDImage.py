Timer unit: 0.1 s

Total time: 0.018031 s
File: /Users/Vincent/Github/fireball2-etc/notebooks/Observation.py
Function: SimulateFIREBallemCCDImage at line 539

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   539                                               def SimulateFIREBallemCCDImage(self,  Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source="Slit", size=[100, 100], OSregions=[0, 100], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=True,atmlambda=True,fraction_lya=0.05, Full_well=60, conversion_gain=1, Throughput_FWHM=20, Altitude=35):
   540                                                   # self.EM_gain=1500; Bias=0; self.RN=80; self.CIC_charge=1; p_sCIC=0; self.Dard_current=1/3600; self.smearing=1; SmearExpDecrement=50000; self.exposure_time=50; flux=1; self.Sky=4; source="Spectra m=17"; Rx=8; Ry=8;  size=[100, 100]; OSregions=[0, 120]; name="Auto"; spectra="Spectra m=17"; cube="-"; n_registers=604; save=False;self.readout_time=5;stack=100;self.QE=0.5
   541         1          0.0      0.0      0.0          from astropy.modeling.functional_models import Gaussian2D, Gaussian1D
   542         1          0.0      0.0      0.1          from scipy.sparse import dia_matrix
   543         1          0.0      0.0      0.0          from scipy.interpolate import interp1d
   544        45          0.0      0.0      1.3          for key in list(instruments["Charact."]) + ["Signal_el","N_images_true","Dark_current_f","sky"]:
   545        44          0.0      0.0      1.1              if hasattr(self,key ):
   546        35          0.0      0.0      1.0                  if (type(getattr(self,key)) != float) & (type(getattr(self,key)) != int) &  (type(getattr(self,key)) != np.float64):
   547         7          0.0      0.0      0.3                      setattr(self, key,getattr(self,key)[self.i])
   548                                                               # print(getattr(self,key))
   549                                                               # print(key, self.i)
   550                                                               # print(getattr(self,key)[self.i])
   551                                                               # a = getattr(self,key)[self.i]
   552                                                               # print(a,self.setattr(key),)
   553                                           
   554                                           
   555         1          0.0      0.0      0.0          conv_gain=conversion_gain
   556         1          0.0      0.0      0.0          OS1, OS2 = OSregions
   557                                                   # ConversionGain=1
   558         1          0.0      0.0      0.0          ConversionGain = conv_gain
   559         1          0.0      0.0      0.0          Bias=0
   560         1          0.0      0.0      0.7          image = np.zeros((size[1], size[0]), dtype="float64")
   561         1          0.0      0.0      0.4          image_stack = np.zeros((size[1], size[0]), dtype="float64")
   562                                           
   563                                                   # self.Dard_current & flux
   564         1          0.0      0.0      0.3          source_im = 0 * image[:, OSregions[0] : OSregions[1]]
   565         1          0.0      0.0      0.3          source_im_wo_atm = 0 * image[:, OSregions[0] : OSregions[1]]
   566         1          0.0      0.0      0.0          lx, ly = source_im.shape
   567         1          0.0      0.0      0.8          y = np.linspace(0, lx - 1, lx)
   568         1          0.0      0.0      0.3          x = np.linspace(0, ly - 1, ly)
   569         1          0.0      0.0      0.7          x, y = np.meshgrid(x, y)
   570                                           
   571                                                   # Source definition. For now the flux is not normalized at all, need to fix this
   572                                                   # Cubes still needs to be implememted, link to detector model or putting it here?
   573                                                   # if os.path.isfile(cube):
   574                                                   # throughput = self.Throughput.value#0.13*0.9
   575                                                   # atm = self.Atmosphere.value#0.45
   576                                                   # area = self.Collecting_area.value/100/100#7854
   577                                                   # dispersion = 1/self.dispersion.value#46.6/10
   578                                                   # wavelength=self.wavelength.value/10 #2000
   579         1          0.0      0.0      0.0          stack = int(self.N_images_true)
   580         1          0.0      0.0      0.0          flux = (self.Signal_el /self.exposure_time)
   581         1          0.0      0.0      0.0          Rx = self.PSF_RMS_det/self.pixel_scale
   582         1          0.0      0.0      0.4          PSF_x = np.sqrt((np.nanmin([self.PSF_source/self.pixel_scale,self.Slitlength/self.pixel_scale]))**2 + (Rx)**2)
   583         1          0.0      0.0      0.0          PSF_λ = np.sqrt(self.PSF_lambda_pix**2 + (self.Line_width/self.dispersion)**2)
   584                                                               
   585                                                   # nsize,nsize2 = size[1],size[0]
   586         1          0.0      0.0      0.0          wave_min, wave_max = 10*self.wavelength - (size[0]/2) * self.dispersion , 10*self.wavelength + (size[0]/2) * self.dispersion
   587                                                   # wavelengths = np.linspace(lmax-nsize2/2*self.dispersion,lmax+nsize2/2*self.dispersion,nsize2)
   588         1          0.0      0.0      0.0          nsize2, nsize = size
   589                                                   # nsize,nsize2 = 100,500
   590         1          0.0      0.0      0.3          wavelengths = np.linspace(wave_min,wave_max,nsize2)
   591         1          0.0      0.0      0.0          if "FIREBall" in self.instrument:
   592         1          0.1      0.1     33.0              trans = Table.read("interpolate/transmission_pix_resolution.csv")
   593         1          0.0      0.0     26.2              QE = Table.read("interpolate/QE_2022.csv")
   594         1          0.0      0.0      1.3              QE = interp1d(QE["wave"]*10,QE["QE_corr"])#
   595                                                       # print(trans["col2"],trans)
   596                                                       # trans["trans_conv"] = np.convolve(trans["col2"],np.ones(int(self.PSF_lambda_pix))/int(self.PSF_lambda_pix),mode="same")
   597         1          0.0      0.0      1.0              trans["trans_conv"] = np.convolve(trans["col2"],np.ones(int(5))/int(5),mode="same")
   598                                                       # trans = trans[:-5]
   599         1          0.0      0.0      6.3              atm_trans =  interp1d(list(trans["col1"]*10),list(trans["trans_conv"]))#
   600                                                       # print(wavelengths.min(),wavelengths.max(), trans["col1"].min(),trans["col1"].max())
   601         1          0.0      0.0      0.4              QE = QE(wavelengths)  if QElambda else self.QE
   602         1          0.0      0.0      0.2              atm_trans = atm_trans(wavelengths)   if (atmlambda & (Altitude<100) ) else self.Atmosphere
   603                                                   else:
   604                                                       trans = Table.read("interpolate/transmission_ground.csv")
   605                                                       atm_trans =  interp1d(list(trans["wave_microns"]*1000), list(trans["transmission"]))#
   606                                                       # print(wavelengths.min(),wavelengths.max(),(trans["wave_microns"]/1000).min(),(trans["wave_microns"]/1000).max())
   607                                                       atm_trans = atm_trans(wavelengths)  if (atmlambda & (Altitude<100) ) else self.Atmosphere
   608                                                       QE = Gaussian1D.evaluate(wavelengths,  self.QE,  self.wavelength*10, Throughput_FWHM )  if QElambda else self.QE
   609                                                       # print(QE)
   610         1          0.0      0.0      0.1          atm_qe =  atm_trans * QE / (self.QE*self.Atmosphere) 
   611                                           
   612         1          0.0      0.0      0.0          length = self.Slitlength/2/self.pixel_scale
   613         1          0.0      0.0      0.4          a = special.erf((length - (np.linspace(0,100,100) - 50)) / np.sqrt(2 * Rx ** 2))
   614         1          0.0      0.0      0.3          b = special.erf((length + (np.linspace(0,100,100) - 50)) / np.sqrt(2 * Rx ** 2))
   615                                           
   616                                           
   617         1          0.0      0.0      0.0          if ("Spectra" in source) | ("Salvato" in source) | ("COSMOS" in source):
   618                                                       if ("baseline" in source.lower()) | (("UVSpectra=" in source) & (self.wavelength>300  )):
   619                                                           # print(PSF_x,PSF_λ)
   620                                                           # print(self.PSF_source,self.pixel_scale,self.PSF_RMS_det,self.pixel_scale)
   621                                                           with_line = flux* Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2, PSF_λ)/ Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2, self.PSF_lambda_pix**2/(PSF_λ**2 + self.PSF_lambda_pix**2)).sum()
   622                                                           # print("QE",QE)
   623                                                           # print("atm_trans",atm_trans)
   624                                                           with_line *= atm_qe
   625                                                           # print(self.Signal_el ,self.exposure_time,flux,with_line)
   626                                                           # source_im[50:55,:] += elec_pix #Gaussian2D.evaluate(x, y, flux, ly / 2, lx / 2, 100 * Ry, Rx, 0)
   627                                                           spatial_profile = Gaussian1D.evaluate(np.arange(size[1]),  1,  size[1]/2, PSF_x)
   628                                                           # spatial_profile += (self.sky/self.exposure_time) * (a + b) / (a + b).ptp()  # 4 * l
   629                                                           # print( length, a, b, Rx )
   630                                                           # print(PSF_x,self.sky,self.exposure_time,length, np.isfinite(length))
   631                                                           profile =  np.outer(with_line,spatial_profile ) /Gaussian1D.evaluate(np.arange(size[1]),  1,  50, Rx**2/(PSF_x**2+Rx**2)).sum()
   632                                                           if np.isfinite(length) & ( (a + b).ptp()>0):
   633                                                               # profile += (self.sky/self.exposure_time) * (a + b) / (a + b).ptp()  * atm_qe
   634                                                               profile +=   np.outer(atm_qe, (self.sky/self.exposure_time) * (a + b) / (a + b).ptp() )
   635                                                           else:
   636                                                               profile +=   np.outer(atm_qe, np.ones(size[1]) *  (self.sky/self.exposure_time)  )  
   637                                                           # print(with_line,spatial_profile ,profile)
   638                                                           # print(self.PSF_source,self.pixel_scale,PSF_x)
   639                                                           # print(flux,with_line ,PSF_x)
   640                                                           source_im = source_im.T
   641                                                           source_im[:,:] += profile
   642                                                           source_im = source_im.T 
   643                                                           # print(source_im,self.Slitlength,profile,atm_trans , QE,self.QE,self.Atmosphere)
   644                                                           #TODO take into account the PSF: += Gaussian1D.evaluate(np.arange(size[1]),  1,  50, PSF_x) with special.erf
   645                                                           #     source_im[50-int(length):50+int(length),:] += self.sky/self.exposure_time  
   646                                                           # else:
   647                                                           #     source_im += self.sky/self.exposure_time  
   648                                                           # print(source_im, self.sky,self.exposure_time  )
   649                                           
   650                                                       # source_im[50:55,:] += elec_pix #Gaussian2D.evaluate(x, y, flux, ly / 2, lx / 2, 100 * Ry, Rx, 0)
   651                                           
   652                                           
   653                                           
   654                                                       elif "mNUV=" in source:
   655                                                           #%%
   656                                                           mag=float(source.split("mNUV=")[-1])
   657                                                           factor_lya = fraction_lya
   658                                                           flux = 10**(-(mag-20.08)/2.5)*2.06*1E-16/((6.62E-34*300000000/(self.wavelength*0.0000000001)/0.0000001))
   659                                                           elec_pix = flux * self.Throughput  * self.Collecting_area*100*100 *self.dispersion  * trans * QE# * self.Atmosphere * self.QE # should not be multiplied by self.exposure_time time here
   660                                                           with_line = elec_pix*(1-factor_lya) + factor_lya * (3700/1)*elec_pix* Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2,PSF_λ)/ Gaussian1D.evaluate(np.arange(size[0]),  1,  size[0]/2, PSF_λ).sum()
   661                                                           # source_im[50:55,:] += elec_pix #Gaussian2D.evaluate(x, y, flux, ly / 2, lx / 2, 100 * Ry, Rx, 0)
   662                                                           profile =  np.outer(with_line,Gaussian1D.evaluate(np.arange(size[1]),  1,  size[1]/2, PSF_x) /Gaussian1D.evaluate(np.arange(size[1]),  1,  size[1]/2, Rx).sum())
   663                                                           source_im = source_im.T
   664                                                           source_im[:,:] += profile
   665                                                           # source_im = source_im.T
   666                                                           # a = Table(data=([np.linspace(1500,2500,nsize2),np.zeros(nsize2)]),names=("WAVELENGTH","e_pix_sec"))
   667                                                           # a["e_pix_sec"] = elec_pix*(1-factor_lya) + factor_lya * (3700/1)*elec_pix* Gaussian1D.evaluate(a["WAVELENGTH"],  1,  line["wave"], 8) 
   668                                                           # f = interp1d(a["WAVELENGTH"],a["e_pix_sec"])
   669                                                           # profile =   Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, Rx) /Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, Rx).sum()
   670                                                           # subim = np.zeros((nsize2,nsize))
   671                                                           # wavelengths = np.linspace(2060-yi*dispersion,2060+(1000-yi)*dispersion,nsize2)
   672                                                           # source_im[int(xi-nsize/2):int(xi+nsize/2), OSregions[0] : OSregions[1]] +=  (subim+profile).T*f(wavelengths) * atm_trans(wavelengths) * self.QE(wavelengths)
   673                                                           # source_im_wo_atm[int(xi-nsize/2):int(xi+nsize/2), OSregions[0] : OSregions[1]] +=  (subim+profile).T*f(wavelengths) #* atm_trans(wavelengths)
   674                                                       else:
   675                                                           # for file in glob.glob("/Users/Vincent/Downloads/FOS_spectra/FOS_spectra_for_FB/CIV/*.fits"):
   676                                           
   677                                                           # print(wave_min, wave_max)
   678                                                           if "_" not in source:
   679                                                               try:
   680                                                                   a = Table.read("Spectra/h_%sfos_spc.fits"%(source.split(" ")[1]))
   681                                                                   flux_name,wave_name ="FLUX", "WAVELENGTH"
   682                                                               except FileNotFoundError: 
   683                                                                   a = Table.read("/Users/Vincent/Github/notebooks/Spectra/h_%sfos_spc.fits"%(source.split(" ")[-1]))
   684                                                                   # slits = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Targets/2022/" + field).to_pandas()
   685                                                                   # trans = Table.read("/Users/Vincent/Github/FIREBall_IMO/Python Package/FireBallIMO-1.0/FireBallIMO/transmission_pix_resolution.csv")
   686                                                                   # self.QE = Table.read("interpolate/QE_2022.csv")
   687                                                               a["photons"] = a[flux_name]/9.93E-12   
   688                                                               a["e_pix_sec"]  = a["photons"] * self.Throughput * self.Atmosphere  * self.Collecting_area*100*100 *self.dispersion
   689                                                           elif "COSMOS" in source:
   690                                                               a = Table.read("Spectra/GAL_COSMOS_SED/%s.txt"%(source.split(" ")[1]),format="ascii")
   691                                                               wave_name,flux_name ="col1", "col2"
   692                                                               mask = (a[wave_name]>wave_min - 100) & (a[wave_name]<wave_max+100)
   693                                                               a = a[mask]
   694                                                               a["e_pix_sec"] = a[flux_name] * flux / np.nanmax(a[flux_name])
   695                                                           elif "Salvato" in source:
   696                                                               a = Table.read("Spectra/Salvato/%s.txt"%(source.split(" ")[1]),format="ascii")
   697                                                               wave_name,flux_name ="col1", "col2"
   698                                                               mask = (a[wave_name]>wave_min - 100) & (a[wave_name]<wave_max+100)
   699                                                               a = a[mask]
   700                                                               a["e_pix_sec"] = a[flux_name] * flux / np.nanmax(a[flux_name])
   701                                                           mask = (a[wave_name]>wave_min) & (a[wave_name]<wave_max)
   702                                                           slits = None #Table.read("Targets/2022/" + field).to_pandas()
   703                                                           source_im=np.zeros((nsize,nsize2))
   704                                                           source_im_wo_atm=np.zeros((nsize2,nsize))
   705                                                           # mask = (a[wave_name]>1960) & (a[wave_name]<2280)
   706                                                           # lmax = a[wave_name][mask][np.argmax( a["e_pix_sec"][mask])]
   707                                                           # plt.plot( a["WAVELENGTH"],a["e_pix_sec"])
   708                                                           # plt.plot( a["WAVELENGTH"][mask],a["e_pix_sec"][mask])
   709                                           
   710                                                           f = interp1d(a[wave_name],a["e_pix_sec"])#
   711                                                           profile =   Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, PSF_x) /Gaussian1D.evaluate(np.arange(nsize),  1,  nsize/2, PSF_x).sum()
   712                                                           subim = np.zeros((nsize2,nsize))
   713                                                           source_im[:,:] +=  (subim+profile).T*f(wavelengths) * atm_trans * QE
   714                                           
   715                                                           if np.isfinite(length) & ( (a + b).ptp()>0):
   716                                                               # profile += (self.sky/self.exposure_time) * (a + b) / (a + b).ptp()  * atm_qe
   717                                                               profile +=   np.outer(atm_qe, (self.sky/self.exposure_time) * (a + b) / (a + b).ptp() )
   718                                                           else:
   719                                                               profile +=   np.outer(atm_qe, np.ones(size[1]) *  (self.sky/self.exposure_time)  )  
   720                                           
   721                                                           # source_im_wo_atm[:,:] +=  (subim+profile).T*f(wavelengths) #* atm_trans(wavelengths)
   722                                                           if 1==0:
   723                                                               fig,(ax0,ax1,ax2) = plt.subplots(3,1)
   724                                                               ax0.fill_between(wavelengths, profile.max()*f(wavelengths),profile.max()* f(wavelengths) * atm_trans(wavelengths),label="Atmosphere impact",alpha=0.3)
   725                                                               ax0.fill_between(wavelengths, profile.max()*f(wavelengths)* atm_trans(wavelengths)*QE(wavelengths),profile.max()* f(wavelengths) * atm_trans(wavelengths),label="self.QE impact",alpha=0.3)
   726                                                               ax1.plot(wavelengths,f(wavelengths)/f(wavelengths).ptp(),label="Spectra")
   727                                                               ax1.plot(wavelengths, f(wavelengths)* atm_trans(wavelengths)/(f(wavelengths)* atm_trans(wavelengths)).ptp(),label="Spectra * Atm")
   728                                                               ax1.plot(wavelengths, f(wavelengths)* atm_trans(wavelengths)*QE/( f(wavelengths)* atm_trans*QE).ptp(),label="Spectra * Atm * self.QE")
   729                                                               ax2.plot(wavelengths,atm_trans(wavelengths) ,label="Atmosphere")
   730                                                               ax2.plot(wavelengths,QE ,label="self.QE")
   731                                                               ax0.legend()
   732                                                               ax1.legend()
   733                                                               ax2.legend()
   734                                                               ax0.set_ylabel("e/pix/sec")
   735                                                               ax1.set_ylabel("Moself.RNalized prof")
   736                                                               ax2.set_ylabel("%")
   737                                                               ax2.set_xlabel("wavelength")
   738                                                               ax0.set_title(source.split(" ")[-1])
   739                                                               fig.savefig("/Users/Vincent/Github/notebooks/Spectra/h_%sfos_spc.png"%(source.split(" ")[-1]))
   740                                                               plt.show()
   741         1          0.0      0.0      1.2          source_im = self.Dark_current_f  + self.extra_background * int(self.exposure_time)/3600 +  source_im  * int(self.exposure_time)
   742                                                   # print(self.Dark_current_f, self.extra_background , int(self.exposure_time)/3600 ,  source_im  , int(self.exposure_time))
   743         1          0.0      0.0      0.3          source_im_wo_atm = self.Dark_current_f + self.extra_background * int(self.exposure_time)/3600 +  source_im_wo_atm * int(self.exposure_time)
   744         1          0.0      0.0      0.0          y_pix=1000
   745                                                   # print(len(source_im),source_im.shape)
   746         1          0.0      0.0      0.0          self.long = False
   747         1          0.0      0.0      0.1          if (self.readout_time/self.exposure_time > 0.2) & (self.long):
   748                                                       # print(source_im)
   749                                                       cube = np.array([(self.readout_time/self.exposure_time/y_pix)*np.vstack((np.zeros((i,len(source_im))),source_im[::-1,:][:-i,:]))[::-1,:] for i in np.arange(1,len(source_im))],dtype=float)
   750                                                       source_im = source_im+np.sum(cube,axis=0)
   751         1          0.0      0.0      0.0          if self.cosmic_ray_loss_per_sec is None:
   752                                                       self.cosmic_ray_loss_per_sec = np.minimum(0.005*(self.exposure_time+self.readout_time/2),1)#+self.readout_time/2
   753                                                   # stack = np.max([int(stack * (1-self.cosmic_ray_loss_per_sec)),1])
   754         1          0.0      0.0      0.0          stack = int(self.N_images_true)
   755         1          0.0      0.0      0.1          cube_stack = -np.ones((stack,size[1], size[0]), dtype="int32")
   756                                           
   757                                                   # print(self.cosmic_ray_loss_per_sec)
   758         1          0.0      0.0      0.0          n_smearing=6
   759                                                   # image[:, OSregions[0] : OSregions[1]] += source_im
   760                                                   # print(image[:, OSregions[0] : OSregions[1]].shape,source_im.shape)
   761         1          0.0      0.0      0.0          if (self.EM_gain>1) & (self.CIC_charge>0):
   762         1          0.0      0.0      4.9              image[:, OSregions[0] : OSregions[1]] += np.random.gamma( np.random.poisson(source_im) + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)
   763                                                   else:
   764                                                       # print(source_im)
   765                                                       image[:, OSregions[0] : OSregions[1]] += np.random.poisson(source_im)
   766                                                   # take into acount CR losses
   767                                                   #18%
   768                                                   # image_stack[:, OSregions[0] : OSregions[1]] = np.nanmean([np.where(np.random.rand(size[1], OSregions[1]-OSregions[0]) < self.cosmic_ray_loss_per_sec/n_smearing,np.nan,1) * (np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)) for i in range(int(stack))],axis=0)
   769         1          0.0      0.0      0.4          image_stack[:, OSregions[0] : OSregions[1]] = np.mean([(np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand(size[1], OSregions[1]-OSregions[0])<self.CIC_charge,dtype=int) , self.EM_gain)) for i in range(int(stack))],axis=0)
   770                                                   
   771                                                   # a = (np.where(np.random.rand(int(stack), size[1],OSregions[1]-OSregions[0]) < self.cosmic_ray_loss_per_sec/n_smearing,np.nan,1) * np.array([ (np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand( OSregions[1]-OSregions[0],size[1]).T<self.CIC_charge,dtype=int) , self.EM_gain))  for i in range(int(stack))]))
   772                                                   # Addition of the phyical image on the 2 overscan regions
   773                                                   #image += source_im2
   774         1          0.0      0.0      0.0          if p_sCIC>0:
   775                                                       image +=  np.random.gamma( np.array(np.random.rand(size[1], size[0])<p_sCIC,dtype=int) , np.random.randint(1, n_registers, size=image.shape))
   776                                                       #30%
   777                                                       image_stack += np.random.gamma( np.array(np.random.rand(size[1], size[0])<int(stack)*p_sCIC,dtype=int) , np.random.randint(1, n_registers, size=image.shape))
   778         1          0.0      0.0      0.0          if self.counting_mode:
   779                                                       a = np.array([ (np.random.gamma(np.random.poisson(source_im)  + np.array(np.random.rand( OSregions[1]-OSregions[0],size[1]).T<self.CIC_charge,dtype="int32") , self.EM_gain))  for i in range(int(stack))])
   780                                                       cube_stack[:,:, OSregions[0] : OSregions[1]] = a
   781                                                       cube_stack += np.random.gamma( np.array(np.random.rand(int(stack),size[1], size[0])<int(stack)*p_sCIC,dtype=int) , np.random.randint(1, n_registers, size=image.shape)).astype("int32")
   782                                                       # print(cube_stack.shape)
   783                                                   #         # addition of pCIC (stil need to add sCIC before EM registers)
   784                                                   #         prob_pCIC = np.random.rand(size[1], size[0])  # Draw a number prob in [0,1]
   785                                                   #         image[prob_pCIC < self.CIC_charge] += 1
   786                                                   #         source_im2_stack[prob_pCIC < p_pCIC*stack] += 1
   787                                           
   788                                                   #         # EM amp (of source + self.Dard_current + pCIC)
   789                                                   #         id_nnul = image != 0
   790                                                   #         image[id_nnul] = np.random.gamma(image[id_nnul], self.EM_gain)
   791                                                           # Addition of sCIC inside EM registers (ie partially amplified)
   792                                                   #         prob_sCIC = np.random.rand(size[1], size[0])  # Draw a number prob in [0,1]
   793                                                   #         id_scic = prob_sCIC < p_sCIC  # sCIC positions
   794                                                   #         # partial amplification of sCIC
   795                                                   #         register = np.random.randint(1, n_registers, size=id_scic.sum())  # Draw at which stage of the EM register the electoself.RN is created
   796                                                   #         image[id_scic] += np.random.exponential(np.power(self.EM_gain, register / n_registers))
   797                                                       # semaring post EM amp (sgest noise reduction)
   798                                                       #TODO must add self.smearing for cube!
   799         1          0.0      0.0      0.0          if self.smearing > 0:
   800                                                       # self.smearing dependant on flux
   801                                                       #2%
   802         1          0.0      0.0      8.0              smearing_kernels = variable_smearing_kernels(image, self.smearing, SmearExpDecrement)
   803         1          0.0      0.0      0.0              offsets = np.arange(n_smearing)
   804         1          0.0      0.0      0.6              A = dia_matrix((smearing_kernels.reshape((n_smearing, -1)), offsets), shape=(image.size, image.size))
   805                                           
   806         1          0.0      0.0      0.4              image = A.dot(image.ravel()).reshape(image.shape)
   807         1          0.0      0.0      1.6              image_stack = A.dot(image_stack.ravel()).reshape(image_stack.shape)
   808                                                   #     if self.readout_time > 0:
   809                                                   #         # self.smearing dependant on flux
   810                                                   #         self.smearing_kernels = variable_smearing.smearing_keself.RNels(image.T, self.readout_time, SmearExpDecrement)#.swapaxes(1,2)
   811                                                   #         offsets = np.arange(n_smearing)
   812                                                   #         A = dia_matrix((self.smearing_kernels.reshape((n_smearing, -1)), offsets), shape=(image.size, image.size))#.swapaxes(0,1)
   813                                                   #         image = A.dot(image.ravel()).reshape(image.shape)#.T
   814                                                   #         image_stack = A.dot(image_stack.ravel()).reshape(image_stack.shape)#.T
   815         1          0.0      0.0      0.0          type_ = "int32"
   816         1          0.0      0.0      0.0          type_ = "float64"
   817         1          0.0      0.0      1.4          readout = np.random.normal(Bias, self.RN, (size[1], size[0]))
   818         1          0.0      0.0      1.5          readout_stack = np.random.normal(Bias, self.RN/np.sqrt(int(stack)), (size[1], size[0]))
   819         1          0.0      0.0      0.0          if self.counting_mode:
   820                                                       readout_cube = np.random.normal(Bias, self.RN, (int(stack),size[1], size[0])).astype("int32")
   821                                                       # print((np.random.rand(source_im.shape[0], source_im.shape[1]) < self.cosmic_ray_loss_per_sec).mean())
   822                                                       #TOKEEP  for cosmic ray masking readout[np.random.rand(source_im.shape[0], source_im.shape[1]) < self.cosmic_ray_loss_per_sec]=np.nan
   823                                                       #print(np.max(((image + readout) * ConversionGain).round()))
   824                                                   #     if np.max(((image + readout) * ConversionGain).round()) > 2 ** 15:
   825         1          0.0      0.0      0.1          imaADU_wo_RN = (image * ConversionGain).round().astype(type_)
   826         1          0.0      0.0      0.3          imaADU_RN = (readout * ConversionGain).round().astype(type_)
   827         1          0.0      0.0      1.1          imaADU = ((image + 1*readout) * ConversionGain).round().astype(type_)
   828         1          0.0      0.0      0.1          imaADU_stack = ((image_stack + 1*readout_stack) * ConversionGain).round().astype(type_)
   829                                                   # print(image_stack,readout_stack)
   830         1          0.0      0.0      0.0          if self.counting_mode:
   831                                                       imaADU_cube = ((cube_stack + 1*readout_cube) * ConversionGain).round().astype("int32")
   832                                                   else:
   833         1          0.0      0.0      0.0              imaADU_cube = imaADU_stack
   834                                                   # print(imaADU_cube.shape)
   835         1          0.0      0.0      0.1          imaADU[imaADU>Full_well*1000] = np.nan
   836         1          0.0      0.0      0.0          return imaADU, imaADU_stack, imaADU_cube, source_im, source_im_wo_atm#imaADU_wo_RN, imaADU_RN