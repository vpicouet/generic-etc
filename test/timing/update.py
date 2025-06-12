Timer unit: 0.1 s

Total time: 1.38002 s
File: /var/folders/m8/f6l41h_51qxdzrz8p1xqr3f80000gp/T/ipykernel_57347/2605358845.py
Function: update at line 429

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   429                                               def update(self,x_axis, counting_mode,Sky,acquisition_time,Signal,EM_gain,RN,CIC_charge,Dard_current,exposure,smearing,temperature,follow_temp,fwhm,QE,extra_background, log,xlog,
   430                                               Collecting_area, pixel_scale, Throughput, Spectral_resolution, Slitwidth, dispersion,
   431                                               PSF_source,Line_width,wavelength,Δλ,Δx, Atmosphere, pixel_size,cosmic_ray_loss_per_sec,lambda_stack, change, 
   432                                               spectra,units,Throughput_FWHM, QElambda, atmlambda, fraction_lya
   433                                               ):
   434                                                   """
   435                                                   Update values in the ETC plot
   436                                                   """
   437         1          0.0      0.0      0.0          if self.change.value:
   438         1          0.0      0.0      0.0              PSF_RMS_mask=fwhm[0]
   439         1          0.0      0.0      0.0              PSF_RMS_det=fwhm[1]
   440         1          0.0      0.0      0.0              readout_time=exposure[0]
   441         1          0.0      0.0      0.0              exposure_time=exposure[1]
   442                                           
   443         1          0.0      0.0      0.2              with self.out1:
   444         1          0.0      0.0      0.0                  if  "OBSERVED SOURCE" in self.x_axis.value:
   445                                                               self.x_axis.value, x_axis = "Signal", "Signal"
   446         1          0.0      0.0      0.0                  elif  "OBSERVATION STRATEGY" in self.x_axis.value:
   447                                                               self.x_axis.value, x_axis = "Atmosphere","Atmosphere"
   448         1          0.0      0.0      0.0                  elif  "INSTRUMENT DESIGN" in self.x_axis.value:
   449                                                               self.x_axis.value, x_axis = "Collecting_area","Collecting_area"
   450         1          0.0      0.0      0.0                  elif  "SPECTROGRAPH DESIGN" in self.x_axis.value:
   451                                                               self.x_axis.value, x_axis = "Spectral_resolution","Spectral_resolution"
   452         1          0.0      0.0      0.0                  elif  "DECTECTOR PERFORMANCE" in self.x_axis.value:
   453                                                               self.x_axis.value, x_axis = "RN", "RN"
   454         1          0.0      0.0      0.0                  elif  "AMPLIFIED" in self.x_axis.value:
   455                                                               self.x_axis.value, x_axis = "EM_gain","EM_gain"
   456                                           
   457                                           
   458                                           
   459         1          0.0      0.0      0.0                  args, _, _, locals_ = inspect.getargvalues(inspect.currentframe())
   460         1          0.0      0.0      0.0                  if x_axis in locals_:
   461         1          0.0      0.0      0.0                      value = locals_[x_axis]
   462                                                           else:
   463                                                               value = getattr(self,x_axis)
   464                                                               if (type(value) != float) & (type(value) != int):
   465                                                                   value = instruments_dict[self.instrument.value][x_axis]
   466                                                               # error here! getattr(self,x_axis) can be 
   467                                           
   468         1          0.0      0.0      0.0                  names = ["Signal","Dark current","Sky", "CIC", "Read noise","Extra Background"]
   469         1          0.0      0.0      0.0                  if follow_temp:
   470                                                               self.Dard_current.value = 10**self.dark_poly(temperature)
   471                                                               self.smearing.value = self.smearing_poly(temperature)
   472                                                               self.CIC_charge.value = self.CIC_poly(self.temperature.value)
   473                                           
   474         1          0.0      0.0      0.0                  self.smearing.layout.visibility = 'visible' if ("FIREBall-2" in self.instrument.value) & (self.counting_mode.value)    else 'hidden'
   475         1          0.0      0.0      0.0                  self.temperature.layout.visibility = 'visible' if ("FIREBall-2" in self.instrument.value) &  (self.follow_temp.value)  else 'hidden'
   476                                           
   477         1          0.0      0.0      0.0                  if x_axis == 'temperature':
   478                                                               temperature=np.linspace(self.temperature.min, self.temperature.max)
   479                                                               Dard_current = 10**self.dark_poly(temperature)
   480                                                               # smearing = np.poly1d([-0.0306087, -2.2226087])(temperature)
   481                                                               smearing = self.smearing_poly(temperature)
   482                                                           # d = {name:np.linspace(rgetattr(self, '%s.min'%(name)), rgetattr(self, '%s.max'%(name))   name for self.fb_options_no_temp}
   483                                           
   484         1          0.0      0.0      0.0                  self.len_xaxis = 50
   485         1          0.0      0.0      0.0                  def space(a, b):
   486                                                               if (self.xlog.value) & (a>=0):
   487                                                                   if a==0:
   488                                                                       y = np.logspace(np.log10(np.max([a,0.0001])),np.log10(b),self.len_xaxis) 
   489                                                                   else:
   490                                                                       y = np.logspace(np.log10(a),np.log10(b),self.len_xaxis) 
   491                                               
   492                                                               else:
   493                                                                   y = np.linspace(a,b,self.len_xaxis)
   494                                                               return y
   495                                           
   496                                                           
   497                                                           # globals()[x_axis] = np.linspace(rgetattr(self, '%s.min'%(x_axis)), rgetattr(self, '%s.max'%(x_axis)) )
   498                                                           # print(eval("x_axis"))
   499                                                           # exec("%s = np.linspace(self.%s.min,self.%s.max)"%(x_axis,x_axis,x_axis))
   500                                                           # print(eval("x_axis"))
   501         1          0.0      0.0      0.0                  if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out1):  #1==1: #
   502                                           
   503                                                               if x_axis == 'exposure_time':
   504                                                                   exposure_time=space(0.1,self.time_max)
   505                                                               elif x_axis == 'Sky':
   506                                                                   # Sky=np.logspace(-19,-15)
   507                                                                   Sky = space(10**self.Sky.min,10**self.Sky.max)
   508                                                               elif x_axis == 'Signal':
   509                                                                   # Signal=np.logspace(-19,-13)
   510                                                                   Signal=space(10**self.Signal.min,10**self.Signal.max)
   511                                                               elif x_axis == 'EM_gain':
   512                                                                   # EM_gain=np.linspace(self.EM_gain.min,self.EM_gain.max)
   513                                                                   EM_gain=space(self.EM_gain.min,self.EM_gain.max)
   514                                                               elif x_axis == 'acquisition_time':
   515                                                                   acquisition_time=space(self.acquisition_time.min,self.acquisition_time.max)
   516                                                               elif x_axis == 'RN':
   517                                                                   RN=space(self.RN.min,self.RN.max)
   518                                                                   # print(self.RN.min,self.RN.max)
   519                                                                   # print(RN)
   520                                                               elif x_axis == 'CIC_charge':
   521                                                                   CIC_charge=space(0.001,self.CIC_charge.max)
   522                                                               elif x_axis == 'Dard_current':
   523                                                                   Dard_current=space(self.Dard_current.min,self.Dard_current.max)
   524                                                                   # Dard_current = np.logspace(-1,np.log10(self.Dard_current.max))
   525                                                               elif x_axis == 'readout_time':
   526                                                                   readout_time=space(self.exposure.min,exposure_time)
   527                                                               elif x_axis == 'smearing':
   528                                                                   smearing=space(self.smearing.min,self.smearing.max)
   529                                                               elif x_axis == 'temperature':
   530                                                                   temperature=space(self.temperature.min,self.temperature.max)
   531                                                               elif x_axis == 'QE':
   532                                                                   QE=space(self.QE.min,self.QE.max)
   533                                                               elif x_axis == 'PSF_RMS_mask':
   534                                                                   PSF_RMS_mask=space(self.fwhm.min,self.fwhm.value[1])
   535                                                               elif x_axis == 'PSF_RMS_det':
   536                                                                   PSF_RMS_det=space(self.fwhm.min,self.fwhm.max)
   537                                                               elif x_axis == 'extra_background':
   538                                                                   extra_background=space(self.extra_background.min,self.extra_background.max)
   539                                                               elif x_axis == 'Δx':
   540                                                                   Δx=space(self.Δx.min,self.Δx.max)#arcseconds
   541                                                               elif x_axis == 'Δλ':
   542                                                                   Δλ=space(self.Δλ.min,self.Δλ.max)#angstrom
   543                                                               elif x_axis == 'Throughput':
   544                                                                   Throughput=space(self.Throughput.min,self.Throughput.max)
   545                                                               elif x_axis == 'Atmosphere':
   546                                                                   Atmosphere=space(self.Atmosphere.min,self.Atmosphere.max)
   547                                                               elif x_axis == 'Line_width':
   548                                                                   Line_width=space(self.Line_width.min,self.Line_width.max)
   549                                                               elif x_axis == 'PSF_source':
   550                                                                   PSF_source =  space(self.PSF_source.min,self.PSF_source.max)
   551                                                               elif x_axis == 'pixel_scale':
   552                                                                   pixel_scale = space(self.pixel_scale.min,self.pixel_scale.max)
   553                                                               elif x_axis == 'pixel_size':
   554                                                                   pixel_size = space(self.pixel_size.min,self.pixel_size.max)
   555                                                               elif x_axis == 'wavelength':
   556                                                                   wavelength =space(self.wavelength.min,self.wavelength.max)
   557                                                               elif x_axis == 'Slitwidth':
   558                                                                   Slitwidth =space(self.Slitwidth.min,self.Slitwidth.max)
   559                                                               elif x_axis == 'Spectral_resolution':
   560                                                                   Spectral_resolution = space(self.Spectral_resolution.min,self.Spectral_resolution.max)
   561                                                               elif x_axis == 'dispersion':
   562                                                                   dispersion = space(self.dispersion.min,self.dispersion.max)
   563                                                               elif x_axis == 'lambda_stack':
   564                                                                   # lambda_stack = np.linspace(np.log10(self.dispersion.value),np.log10(self.Bandwidth))
   565                                                                   lambda_stack = space(self.dispersion.value,self.Bandwidth)
   566                                                               elif x_axis == 'Collecting_area':
   567                                                                   Collecting_area = space(self.Collecting_area.min,self.Collecting_area.max)
   568                                                                   # Collecting_area = np.linspace(0.01, 10)
   569                                                               elif x_axis == "cosmic_ray_loss_per_sec":
   570                                                                   cosmic_ray_loss_per_sec=space(0,1/exposure[1])
   571                                           
   572                                           
   573                                           
   574         1          0.0      0.0      0.0                  args, _, _, locals_ = inspect.getargvalues(inspect.currentframe())
   575         1          0.0      0.0      0.0                  try:
   576         1          0.0      0.0      0.0                      new_value = locals_[x_axis]
   577                                                           except KeyError:
   578                                                               new_value = getattr(self,x_axis)
   579         1          0.0      0.0      0.0                  arg = np.argmin(abs(new_value - value))
   580                                                           # print("instrument=",self.instrument.value," \nAtmosphere",self.Atmosphere, "\nThroughput=",self.Throughput,"\nSky=",Sky, "\nacquisition_time=",acquisition_time,"\ncounting_mode=",counting_mode,"\nSignal=",Signal,"\nEM_gain=",EM_gain,"RN=",RN,"CIC_charge=",CIC_charge,"Dard_current=",Dard_current,"\nreadout_time=",readout_time,"\nsmearing=",smearing,"\nextra_background=",extra_background,"\ntemperature=",temperature,"\nPSF_RMS_mask=",PSF_RMS_mask,"\nPSF_RMS_det=",PSF_RMS_det,"\nQE=",QE,"\ncosmic_ray_loss_per_sec=",self.cosmic_ray_loss_per_sec,"\nlambda_stack",self.lambda_stack,"\nSlitwidth",self.Slitwidth, "\nBandwidth",self.Bandwidth,"\nPSF_source",self.PSF_source,"\nCollecting_area",self.Collecting_area)
   581         2          0.0      0.0      0.2                  self.new = Observation(instrument=self.instrument.value,  exposure_time=exposure_time,Sky=Sky, acquisition_time=acquisition_time,counting_mode=counting_mode,Signal=Signal,EM_gain=EM_gain,RN=RN,CIC_charge=CIC_charge,Dard_current=Dard_current,readout_time=readout_time,smearing=smearing,extra_background=extra_background,i=arg,temperature=temperature,PSF_RMS_mask=PSF_RMS_mask,PSF_RMS_det=PSF_RMS_det,QE=QE,cosmic_ray_loss_per_sec=cosmic_ray_loss_per_sec, Throughput=Throughput, Atmosphere=Atmosphere,lambda_stack=lambda_stack,Slitwidth=Slitwidth, Bandwidth=self.Bandwidth,PSF_source=PSF_source,Collecting_area=Collecting_area,Δx=Δx,Δλ=Δλ,
   582         1          0.0      0.0      0.0                  pixel_scale=pixel_scale, Spectral_resolution=Spectral_resolution,  dispersion=dispersion,
   583         1          0.0      0.0      0.0                  Line_width=Line_width,wavelength=wavelength,  pixel_size=pixel_size,len_xaxis=self.len_xaxis, Slitlength=self.Slitlength)
   584                                                           # self.Signal_el = self.new.Signal_el
   585                                                           
   586                                                       #         self.Total_noise_final = self.factor*np.sqrt(self.signal_noise**2 + self.Dark_current_noise**2  + self.Additional_background_noise**2 + self.Sky_noise**2 + self.CIC_noise**2 + self.RN_final**2   ) #e/  pix/frame
   587                                                       # self.SNR = self.Signal_resolution / self.Total_noise_final
   588                                           
   589         1          0.0      0.0      0.0                  self.colors=self.new.colors
   590                                                           
   591         1          0.0      0.0      0.1                  if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out1): #1==1: #
   592                                                               arg = np.argmin(abs(getattr(self.new,x_axis) - value))
   593                                                               try:
   594                                                                   label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f, SNR=%0.2f'%(self.x_axis.value,float_to_latex(value),float_to_latex(new_value[np.nanargmax(self.new.SNR)]),self.new.SNR[arg],np.nanmax(self.new.SNR))#, self.new.gain_thresholding[arg])
   595                                                               except (TypeError,ValueError) as e:
   596                                                                   # print(e)
   597                                                                   if ("FIREBall" in self.instrument.value) & (self.counting_mode.value):
   598                                                                       label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f, SNR=%0.2f\nT=%0.1f sigma\nSignal kept=%i%%, RN kept=%i%%'%(self.x_axis.value,float_to_latex(value),float_to_latex(new_value[np.nanargmax(self.new.SNR)]),self.new.SNR[arg],np.nanmax(self.new.SNR),self.new.n_threshold[arg], 100*self.new.Photon_fraction_kept[arg], 100*self.new.RN_fraction_kept[arg])#, self.new.gain_thresholding[arg])
   599                                                                   else:
   600                                                                       label = '%s [Best]=%s [%s]\nSNR [Best]=%0.2f, SNR=%0.2f'%(self.x_axis.value,float_to_latex(value),float_to_latex(new_value[np.nanargmax(self.new.SNR)]),self.new.SNR[arg],np.nanmax(self.new.SNR))#, self.new.gain_thresholding[arg])
   601                                           
   602                                                               max_,min_=[],[]
   603                                           
   604                                                               for i,name in enumerate(self.new.names): 
   605                                                                   self.ax0.lines[i].set_xdata(new_value)
   606                                                                   self.ax0.lines[i].set_ydata(self.new.noises[:,i]/self.new.factor)
   607                                                                   # print(self.new.factor,self.new.i)
   608                                                                   # print(self.new.factor[self.new.i])
   609                                                                   if self.new.percents[i,self.new.i] ==np.max(self.new.percents[:,self.new.i]):
   610                                                                       self.ax0.lines[i].set_label(r"$\bf{%s}$: %0.2f (%0.1f%%)"%(name,self.new.noises[self.new.i,i]/self.new.factor[self.new.i],self.new.percents[i,self.new.i]))
   611                                                                   else:
   612                                                                       self.ax0.lines[i].set_label('%s: %0.2f (%0.1f%%)'%(name,self.new.noises[self.new.i,i]/self.new.factor[self.new.i],self.new.percents[i,self.new.i]))
   613                                                                   max_.append(np.nanmax(self.new.noises[:,i]/self.new.factor))
   614                                                                   min_.append(np.nanmin(self.new.noises[:,i]/self.new.factor))
   615                                                               self.ax0.lines[i+1].set_xdata(new_value)
   616                                                               self.ax0.lines[i+1].set_ydata(np.nansum(self.new.noises[:,:-1],axis=1)/self.new.factor)
   617                                                               self.ax0.lines[i+1].set_label('%s: %0.2f (%0.1f%%)'%("Total",np.nansum(self.new.noises[self.new.i,:-1])/self.new.factor[self.new.i],np.nansum(self.new.percents[:,self.new.i])))
   618                                           
   619                                                               self.ax0.legend(loc='upper right')
   620                                                               # self.ax3.set_xlabel(x_axis)
   621                                                               if x_axis in ["exposure_time","readout_time","PSF_RMS_mask","PSF_RMS_det"]:
   622                                                                   self.ax3.set_xlabel(x_axis)
   623                                                               else:
   624                                                                   try:
   625                                                                       self.ax3.set_xlabel(rgetattr(self, '%s.description_tooltip'%(x_axis)) )
   626                                                                   except AttributeError:
   627                                                                       self.ax3.set_xlabel(x_axis + "  [%s]"%(instruments["Unit"][instruments["Charact."]==x_axis][0]))
   628                                           
   629                                                               self.ax3.lines[0].set_data(new_value,  np.log10(self.new.extended_source_5s))
   630                                                               self.ax3.lines[0].set_label("SNR=5 Flux/Pow on one elem resolution (%0.2f-%0.2f)"%(np.log10(self.new.point_source_5s[arg]),np.nanmin(np.log10(self.new.point_source_5s))))
   631                                           
   632                                                               # print(3)
   633                                           
   634                                                               if "FIREBall" in self.instrument.value:
   635                                                                   self.ax3.lines[1].set_data(new_value,  np.log10(self.new.extended_source_5s/np.sqrt(2)))
   636                                                                   self.ax3.lines[1].set_label("Two elem resolution  (%0.2f-%0.2f)"%(np.log10(self.new.point_source_5s[arg]/np.sqrt(2)),np.nanmin(np.log10(self.new.point_source_5s/np.sqrt(2)))))
   637                                                                   if x_axis == 'Sky':
   638                                                                       for line, text, text2, sky,name in zip(self.test,self.text,self.text2,[2e-16,2e-17,2e-18],["2018/2","Nominal","Best"]):
   639                                                                           line.set_xdata(sky)
   640                                                                           text.set_position((sky,1))
   641                                                                           text.set_text(name)
   642                                                                           text2.set_position((sky,self.new.SNR[np.argmin(abs(new_value-sky))]))
   643                                                                           text2.set_text( "SNR=%0.1f"%(self.new.SNR[np.argmin(abs(new_value-sky))]))
   644                                                                   else:
   645                                                                       for line, text, text2, sky,name in zip(self.test,self.text,self.text2,[2e-16,2e-17,2e-18],["2018/2","Nominal","Best"]):
   646                                                                           line.set_xdata(np.nan)
   647                                                                           text.set_position((np.nan,1))
   648                                                                           text.set_text("")
   649                                                                           text2.set_position((np.nan,np.nan))
   650                                                                           text2.set_text( "")
   651                                                                       # self.ax2.clear()
   652                                                                       # self.ax2.axvline(4e-17,c="k",alpha=0.5,ls=":")
   653                                                                       # self.ax2.text(4e-17,1, "Nominal")
   654                                                                       # self.ax2.axvline(1.1e-17,c="k",alpha=0.5,ls=":")
   655                                                                       # self.ax2.text(1.1e-17,1, "Best")
   656                                                                       # self.ax2.axvline(1.1e-17,c="k",alpha=0.5,ls=":")
   657                                                                       # self.ax2.text(1.1e-17,1, "Best")
   658                                                                   # elif x_axis == 'Signal':  
   659                                                                   #     self.ax2.axvline(1.1e-19,c="k",alpha=0.5,ls=":")
   660                                                                   #     self.ax2.text(1.1e-19,1, "Chen 2021 limit")
   661                                                                   #     self.ax2.axvline(9e-14,c="k",alpha=0.5,ls=":")
   662                                                                   #     self.ax2.text(3e-14,1, "IRAS QSO")
   663                                                                   #     self.ax2.axvline(1e-14,c="k",alpha=0.5,ls=":")
   664                                                                   #     self.ax2.text(1e-14,1, "PG QSOs")
   665                                                                   #     self.ax2.axvline(1e-16,c="k",alpha=0.5,ls=":")
   666                                                                   #     self.ax2.text(1e-16,1, "~Other targeted QSOs")
   667                                                               # ax2.axhline(0.1,ls=":",color="k")
   668                                                               else:
   669                                                                   self.ax3.lines[1].set_data(new_value,  np.log10(self.new.extended_source_5s/np.sqrt(2)))
   670                                           
   671                                           
   672                                                               for v in self.v:
   673                                                                   v.set_xdata([value,value])
   674                                                               self.v[-2].set_label(label)
   675                                           
   676                                           
   677                                                               for artist in self.ax2.collections+self.ax1.collections:
   678                                                                   artist.remove()
   679                                                               # print(4)
   680                                                               self.stackplot2 = self.ax2.stackplot(new_value,self.new.snrs * np.array(self.new.noises).T[:-1,:]**2/self.new.Total_noise_final**2,alpha=0.7,colors=self.colors)
   681                                                               labels =  ['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]
   682                                                               self.stackplot1 = self.ax1.stackplot(new_value,  np.array(self.new.electrons_per_pix).T,alpha=0.7,colors=self.colors,labels=labels)
   683                                                               self.ax1.legend(loc='upper right',title="Overall background: %0.3f (%0.1f%%)"%(np.nansum(self.new.electrons_per_pix[self.new.i,1:]),100*np.nansum(self.new.electrons_per_pix[self.new.i,1:])/np.nansum(self.new.electrons_per_pix[self.new.i,:])))
   684                                                               self.ax2.legend(loc='upper right')
   685                                                               self.ax3.legend(loc='upper right',fontsize=8)
   686                                                               self.ax2.set_xlim((np.max([np.min(new_value),1e-6]),np.max(new_value)))
   687                                                               self.ax2.set_xlim((np.min(new_value),np.max(new_value)))
   688                                           
   689                                           
   690                                                               if log:
   691                                                                   self.ax0.set_yscale("log")
   692                                                                   self.ax1.set_yscale("log")
   693                                                                   self.ax2.set_yscale("log")
   694                                                                   # self.ax0.set_ylim(ymin=np.nanmin(self.new.noises[:,:-1]/self.new.factor[:,None]),ymax=np.nanmax(np.nansum(self.new.noises[:,:-1],axis=1)/self.new.factor))
   695                                                                   self.ax0.set_ylim(ymin=0,ymax=np.nanmax(np.nansum(self.new.noises[:,:-1],axis=1)/self.new.factor))
   696                                                                   self.ax1.set_ylim(ymin=np.nanmin( np.array(self.new.electrons_per_pix[:,0])),ymax=np.max(np.sum(getattr(self.new,'electrons_per_pix'),axis=1)))
   697                                                                   self.ax2.set_ylim(ymin=np.nanmin(np.array( self.new.snrs * np.array(self.new.noises).T[:-1,:]**2/self.new.Total_noise_final**2)[:,0]),ymax=np.nanmax(getattr(self.new,'SNR')))
   698                                                               else:
   699                                                                   self.ax0.set_yscale("linear")
   700                                                                   self.ax1.set_yscale("linear")
   701                                                                   self.ax2.set_yscale("linear")
   702                                                                   self.ax0.set_ylim((-0.1,np.nanmax(np.nansum(self.new.noises[:,:-1],axis=1)/self.new.factor)))
   703                                                                   self.ax1.set_ylim((0,np.max(np.sum(getattr(self.new,'electrons_per_pix'),axis=1))))
   704                                                                   self.ax2.set_ylim((0,np.nanmax(getattr(self.new,'SNR'))))
   705                                                               
   706                                                               if xlog:
   707                                                                   try:
   708                                                                       if (rgetattr(self,"%s.min"%(x_axis))<=0) & ( rgetattr(self,"%s.min"%(x_axis))  <  rgetattr(self,"%s.value"%(x_axis)) <  rgetattr(self,"%s.max"%(x_axis))):
   709                                                                           self.ax0.set_xscale("symlog")
   710                                                                       else:
   711                                                                           self.ax0.set_xscale("log")
   712                                                                   except AttributeError:
   713                                                                       self.ax0.set_xscale("log")
   714                                           
   715                                                               else:
   716                                                                   self.ax0.set_xscale("linear")
   717                                           
   718                                                               self.fig.canvas.draw()
   719         1          0.0      0.0      0.0              if self.output_tabs.get_state()["selected_index"]==self.output_tabs.children.index(self.out2): #1==1: #
   720         1          0.0      0.0      0.1                  with self.out2:
   721         1          0.0      0.0      0.0                      if "Spectra mNUV=" in x_axis:
   722                                                                   Signal = float(spectra.split("=")[-1])
   723                                                                   psf_source = 0.1
   724                                                                   self.fraction_lya.layout.visibility = 'visible'
   725                                                               else: 
   726         1          0.0      0.0      0.0                          Signal = 20
   727         1          0.0      0.0      0.0                          psf_source = 4
   728         1          0.0      0.0      0.0                          self.fraction_lya.layout.visibility = 'hidden'
   729                                                               # if counting_mode:
   730                                                               #     self.threshold.layout.visibility = 'visible'
   731                                                               # else:
   732                                                               #     self.threshold.layout.visibility = 'hidden'
   733                                                               # self.new = Observation(exposure_time=exposure_time,Sky=Sky, acquisition_time=acquisition_time,counting_mode=counting_mode,Signal=Signal,EM_gain=EM_gain,RN=RN,CIC_charge=CIC_charge,Dard_current=Dard_current,readout_time=readout_time,smearing=smearing,extra_background=extra_background,QE=QE)
   734         1          0.0      0.0      0.0                      Sky = (self.new.Sky/exposure_time)#[arg]
   735         1          0.0      0.0      0.0                      flux = (self.new.Signal_el/exposure_time)#[arg]
   736                                                               # ros_number, conv_gain, full_well = np.array(re.findall(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?', Throughput_FWHM),dtype=float)
   737                                                               # conv_gain = 1/conv_gain
   738                                                               # print(EM_gain,RN, CIC_charge, (Dard_current+extra_background)/3600, smearing, exposure_time[arg], flux,Sky, x_axis, int(self.new.N_images_true[arg]),readout_time,QE,QElambda,atmlambda,conv_gain,fraction_lya)
   739                                                               # self.im,self.im_stack, self.cube_stack,self.im0,_  =  self.SimulateFIREBallemCCDImage(EmGain=EM_gain, Bias=0, RN=RN, p_pCIC=CIC_charge, p_sCIC=0, Dark=(Dard_current+extra_background)/3600, Smearing=smearing, SmearExpDecrement=50000, exposure=exposure_time, flux=flux, Sky=Sky, source=spectra, Rx=fwhm[-1], Ry=fwhm[1], size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False,stack=int(self.new.N_images_true),readout_time=readout_time,QE=QE,QElambda=QElambda,atmlambda=atmlambda,conv_gain=conv_gain,fraction_lya=fraction_lya)
   740         1          3.9      3.9     28.5                      self.im,self.im_stack, self.cube_stack, self.im0, _ = self.new.SimulateFIREBallemCCDImage(Bias="Auto",  p_sCIC=0,  SmearExpDecrement=50000,  source=spectra,size=[n1, n2], OSregions=[0, max(n2,n1)], name="Auto", spectra="-", cube="-", n_registers=604, save=False, field="targets_F2.csv",QElambda=QElambda,atmlambda=atmlambda,fraction_lya= fraction_lya, Full_well=self.Full_well, conversion_gain=self.conversion_gain, Altitude=self.Altitude,Throughput_FWHM=self.Throughput_FWHM.value)
   741                                           
   742         1          0.0      0.0      0.0                      self.bins=np.linspace(np.nanmin(self.im),np.nanmax(self.im),100)
   743                                           
   744         1          0.0      0.0      0.0                      self.current_cmap.set_bad('red',1.)
   745         1          0.0      0.0      0.0                      if units=="ADU/frame": #ADU/frame: ok basic
   746         1          0.0      0.0      0.0                          factor=1
   747                                                               elif units=="amplified e-/frame": #e-/frame: divide by conversion gain and amplification gain
   748                                                                   factor=1/self.conversion_gain
   749                                                               elif units=="e-/frame": #e-/frame: divide by conversion gain and amplification gain
   750                                                                   factor=1/self.conversion_gain/EM_gain
   751                                                               elif units=="photons/frame": #photons/frame: account for QE
   752                                                                   factor=1/self.conversion_gain/EM_gain/QE
   753                                                               elif units=="amplified e-/hour": #e-/hour: divide by exptime
   754                                                                   factor=1/self.conversion_gain/exposure_time/3600
   755                                                               elif units=="e-/hour": #e-/hour: divide by exptime
   756                                                                   factor=1/self.conversion_gain/EM_gain/exposure_time/3600
   757                                                               elif units=="photons/hour": #photons/hour: divide by exptime, account for QE
   758                                                                   factor=1/self.conversion_gain/EM_gain/QE/exposure_time/3600
   759                                           
   760         1          0.0      0.0      0.1                      im = self.nax.imshow(self.im*factor, aspect="auto",cmap=self.current_cmap)
   761         1          0.1      0.1      0.5                      self.mod = mostFrequent(self.im_stack[:20,:].flatten())
   762                                                               # self.limit = self.mod+threshold*RN
   763         1          0.0      0.0      0.0                      self.limit = self.mod+self.new.n_threshold * RN
   764                                           
   765         1          0.0      0.0      0.0                      if counting_mode:
   766                                                                   stacked_image = np.nansum(self.cube_stack>self.limit,axis=0)
   767                                                                   im0 = self.nax0.imshow(stacked_image, aspect="auto",cmap=self.current_cmap)
   768                                                               else:
   769         1          0.0      0.0      0.1                          im0 = self.nax0.imshow(self.im_stack*factor, aspect="auto",cmap=self.current_cmap)
   770         1          0.1      0.1      0.8                      self.cbar1 = self.fig2.colorbar(im, cax=self.cax, orientation='horizontal')
   771         1          0.1      0.1      0.7                      self.cbar2 = self.fig2.colorbar(im0, cax=self.cax0, orientation='horizontal')
   772         1          0.0      0.0      0.0                      self.cbar1.formatter.set_powerlimits((0, 0))
   773         1          0.0      0.0      0.0                      self.cbar2.formatter.set_powerlimits((0, 0))
   774                                           
   775         1          0.0      0.0      0.0                      labels =  ['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]
   776                                           
   777         1          0.0      0.0      0.0                      if "Spectra m=" not in x_axis:
   778         1          0.0      0.0      0.0                          self.nax1.lines[0].set_ydata(self.im[:,45:-45].mean(axis=1))
   779         1          0.0      0.0      0.0                          stacked_profile = np.mean(im0.get_array().data[:,45:-45],axis=1)
   780         1          0.0      0.0      0.0                          self.nax1.lines[1].set_ydata(self.im[40:-40,:].mean(axis=0)) # spectral direction, only if not spectra
   781         1          0.0      0.0      0.0                          self.nax1bis.lines[1].set_ydata(np.mean(im0.get_array().data[40:-40,:],axis=0))
   782         1          0.0      0.0      0.0                          self.nax.lines[0].set_label(" ")                 
   783                                                               else:
   784                                                                   # self.profile = np.mean(im0.get_array().data[:,:],axis=1)
   785                                                                   stacked_profile = np.mean(im0.get_array().data[:,:],axis=1)
   786                                                                   # print(x_axis)
   787                                                                   spatial_profile = self.im[:,:].mean(axis=1)
   788                                                                   self.nax1.lines[0].set_ydata(np.convolve(spatial_profile,3,mode="same"))
   789                                                                   self.nax1bis.lines[1].set_ydata(np.mean(im0.get_array().data[40:-40,:],axis=0))
   790                                                                   self.nax1.lines[1].set_ydata(self.im[40:-40,:].mean(axis=0)) 
   791                                                                   self.nax.lines[0].set_label("  \n".join(labels)) 
   792         1          0.0      0.0      0.1                      self.nax.legend(loc="upper left",handlelength=0, handletextpad=0, fancybox=True,markerscale=0,fontsize=8)
   793                                                       
   794         1          0.0      0.0      0.0                      try:
   795         1          0.0      0.0      0.1                          self.popt, self.pcov = curve_fit(gaus,np.arange(len(stacked_profile)),stacked_profile,p0=[stacked_profile.ptp(), 50, 5, stacked_profile.min()])
   796                                                               except RuntimeError:
   797                                                                   self.popt = [0,0,0,0]
   798                                                               # self.fit = PlotFit1D(x= np.arange(len(stacked_profile)),y=stacked_profile,deg="gaus", plot_=False,ax=self.nax1bis,c="k",ls=":",P0=[stacked_profile.ptp(), 50, 5, stacked_profile.min()])
   799         1          0.0      0.0      0.0                      self.nax1bis.lines[0].set_ydata(stacked_profile)
   800         1          0.0      0.0      0.0                      self.SNR = self.popt[0]**2/stacked_profile[:20].std()
   801         1          0.0      0.0      0.0                      self.Flux_ADU =  np.sum(gaus( np.arange(len(stacked_profile)),*self.popt)-self.popt[-1]) 
   802                                                               # self.Flux_ADU_counting =  np.sum(-np.log(1-( self.fit["function"]( np.arange(len(stacked_profile)),*self.fit["popt"])-self.fit["popt"][-1] )/(np.exp(-threshold*RN/EM_gain))))
   803         1          0.0      0.0      0.0                      self.e_s_pix = self.Flux_ADU * self.new.dispersion / exposure_time / self.new.N_images_true/self.conversion_gain  if counting_mode else  self.Flux_ADU * self.new.dispersion / EM_gain / exposure_time/self.conversion_gain
   804         1          0.0      0.0      0.0                      self.flux = self.e_s_pix / self.new.Throughput/ self.new.Atmosphere / QE / self.new.Collecting_area
   805         1          0.0      0.0      0.0                      photon_energy_erg = 9.93e-12
   806         1          0.0      0.0      0.0                      self.mag = -2.5*np.log10(self.flux*photon_energy_erg/(2.06*1E-16))+20.08
   807                                                               # self.nax1bis.lines[2].set_label("SNR=%s/%s=%s, mag=%s"%(self.popt[0]**2,stacked_profile[:20].std(),self.SNR ,self.mag[arg]))
   808         1          0.0      0.0      0.0                      self.nax1bis.lines[2].set_label("SNR=%0.1f/%0.1f=%0.1f, mag=%0.1f"%(self.popt[0]**2,stacked_profile[:20].std(),self.SNR ,self.mag))
   809         1          0.0      0.0      0.0                      self.nax1bis.lines[2].set_ydata(gaus( np.arange(len(stacked_profile)),*self.popt))
   810         1          0.0      0.0      0.0                      self.nax.set_title('Single image: FOV = %i" × %iÅ, λ~%iÅ'%(100*self.pixel_scale.value,500*self.dispersion.value, 10*self.wavelength.value))
   811         1          0.0      0.0      0.0                      self.nax0.set_title('Stacked image: Pixel size = %0.2f" × %0.2fÅ'%(self.pixel_scale.value,self.dispersion.value))
   812                                           
   813         1          0.0      0.0      0.3                      self.nax1bis.legend(loc="upper right",fontsize=8)
   814                                           
   815         1          0.0      0.0      0.0                      self.nax1.relim()
   816         1          0.0      0.0      0.0                      self.nax1.autoscale_view()
   817         1          0.0      0.0      0.0                      self.nax1bis.relim()
   818         1          0.0      0.0      0.0                      self.nax1bis.autoscale_view()
   819                                                               # for patch in self.nax2.patches:
   820                                                               #     patch.clear()
   821                                                               # self.nax2.cla()
   822         1          0.0      0.0      0.1                      [b.remove() for b in self.bars2]
   823         1          0.0      0.0      0.0                      [b.remove() for b in self.bars1]
   824         1          0.8      0.8      5.9                      _,_,self.bars1 = self.nax2.hist(self.im.flatten(),bins=self.bins,log=True,alpha=0.3,color=self.l1[0].get_color(),label='Single image')
   825         1          0.8      0.8      5.9                      _,_,self.bars2 = self.nax2.hist(self.im_stack.flatten(),bins=self.bins,log=True,alpha=0.3,color=self.l2[0].get_color(),label='Averaged stack')
   826         1          0.0      0.0      0.1                      self.nax2.set_xlim(self.bins.min(),self.bins.max())
   827                                           
   828                                                               # self.fig2.axes[3].patches.clear()
   829                                                               # self.fig2.axes[3].hist(self.im.flatten(),bins=self.bins,log=True,alpha=0.3,color=self.l1[0].get_color(),label='Single image')
   830                                                               # self.fig2.axes[3].hist(self.im_stack.flatten(),bins=self.bins,log=True,alpha=0.3,color=self.l2[0].get_color(),label='Averaged stack')
   831                                           
   832         1          0.0      0.0      0.0                      self.hw, self.hl =  self.Slitwidth.value/2/self.pixel_scale.value ,  self.Slitlength/2/self.pixel_scale.value
   833         1          0.0      0.0      0.0                      try:
   834         1          0.0      0.0      0.0                          self.nax.lines[1].set_data([250 - self.hw,250 + self.hw,250 + self.hw,250 - self.hw,250 - self.hw],[50 - self.hl,50 - self.hl,50 + self.hl,50 + self.hl,50 - self.hl])
   835         1          0.0      0.0      0.0                          self.nax0.lines[0].set_data([250 - self.hw,250 + self.hw,250 + self.hw,250 - self.hw,250 - self.hw],[50 - self.hl,50 - self.hl,50 + self.hl,50 + self.hl,50 - self.hl])
   836                                                               except ValueError:
   837                                                                   self.nax.lines[1].set_data(np.nan,np.nan)
   838                                                                   self.nax0.lines[0].set_data(np.nan,np.nan)
   839         1          0.0      0.0      0.0                      try:
   840         1          0.0      0.0      0.0                          self.nax0.lines[0].set_label('Slit=%i" × %i"'%(self.Slitwidth.value,self.Slitlength))
   841                                                               except ValueError:
   842                                                                   self.nax0.lines[0].set_label(" ")
   843         1          0.0      0.0      0.0                      self.nax.lines[0].set_label("  \n".join(['%s: %0.3f (%0.1f%%)'%(name,getattr(self.new,"electrons_per_pix")[self.new.i,j],100*getattr(self.new,"electrons_per_pix")[self.new.i,j]/np.sum(getattr(self.new,'electrons_per_pix')[self.new.i,:])) for j,name in enumerate(self.new.names)]))
   844                                           
   845         1          0.0      0.0      0.1                      self.nax.legend(loc="upper left",handlelength=0, handletextpad=0, fancybox=True,markerscale=0,fontsize=8)
   846         1          0.0      0.0      0.1                      self.nax0.legend(loc='upper right',fontsize=8)
   847         1          0.0      0.0      0.2                      self.nax2.legend(loc='upper right',fontsize=8)
   848         1          0.0      0.0      0.0                      self.nax2.lines[0].set_xdata([self.mod,self.mod])
   849         1          0.0      0.0      0.0                      self.nax2.lines[1].set_xdata([self.limit[arg],self.limit[arg]])
   850                                                               # title = 'Signal kept=%i%%, RN kept=%i%%\n'%(100*self.new.Photon_fraction_kept[0], 100*self.new.RN_fraction_kept[0])
   851         1          0.0      0.0      0.0                      if "FIREBall" in  self.instrument.value:
   852         1          0.0      0.0      0.0                          try:
   853         1          0.0      0.0      0.0                              title = 'Signal kept=%i%%, RN kept=%i%%, Signal/tot=%i%%'%(100*self.new.Photon_fraction_kept[0], 100*self.new.RN_fraction_kept[0],100*(np.mean(self.im_stack[40:-40,:])-np.mean(self.im_stack[:20,:]))/np.mean(self.im_stack[40:-40,:]))
   854                                                                   except IndexError as e:
   855                                                                       title = 'Signal kept=%i%%, RN kept=%i%%, Signal/tot=%i%%'%(100*self.new.Photon_fraction_kept, 100*self.new.RN_fraction_kept[0],100*(np.mean(self.im_stack[40:-40,:])-np.mean(self.im_stack[:20,:]))/np.mean(self.im_stack[40:-40,:]))
   856         1          0.0      0.0      0.0                          self.nax2.lines[0].set_label("Bias %s, PC limit %s (%s):\n%s "%(self.mod,self.limit[arg], counting_mode, title))
   857                                                               else:
   858                                                                   self.nax2.lines[0].set_label(" ")
   859                                           
   860                                                               # self.nax2.set_xlim(xmax=)       
   861         1          0.0      0.0      0.2                      self.nax2.legend(fontsize=8,loc="upper right")
   862                                           
   863                                                               # title = '%i stacks, CIC=%0.1f, Dark=%0.1f, RN=%s, Flux=%0.2E, Sky=%0.2E'%(self.new.N_images_true,CIC_charge,Dard_current/3600,RN,flux,Sky)
   864                                                               # self.fig2.suptitle(title)
   865         1          1.5      1.5     11.0                      self.fig2.tight_layout()
   866         1          6.1      6.1     44.3                      self.fig2.canvas.draw()
   867                                                               # print(self.new.Signal_el)