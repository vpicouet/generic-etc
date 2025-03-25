Timer unit: 0.1 s

Total time: 0.009539 s
File: /Users/Vincent/Github/fireball2-etc/notebooks/Observation.py
Function: interpolate_optimal_threshold at line 483

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   483                                               def interpolate_optimal_threshold(self,flux = 0.1,dark_cic_sky_noise=None,plot_=False,title='',i=0):
   484                                                   """
   485                                                   Return the threshold optimizing the SNR
   486                                                   """
   487                                                   #self.Signal_el if np.isscalar(self.Signal_el) else 0.3
   488         2          0.0      0.0      0.1          EM_gain = self.EM_gain #if np.isscalar(self.EM_gain) else self.EM_gain[i]
   489         2          0.0      0.0      0.0          RN= self.RN #if np.isscalar(self.RN) else self.RN[i]#80
   490         2          0.0      0.0      0.0          CIC_noise = self.CIC_noise #if np.isscalar(self.CIC_noise) else self.CIC_noise[i]
   491         2          0.0      0.0      0.0          dark_noise = self.Dark_current_noise #if np.isscalar(self.Dark_current_noise) else self.Dark_current_noise[i]
   492                                                    
   493         2          0.0      0.0      0.0          try:
   494         2          0.0      0.0      0.0              Sky_noise = self.Sky_noise_pre_thresholding #if np.isscalar(self.Sky_noise_pre_thresholding) else self.Sky_noise_pre_thresholding[i]
   495                                                   except AttributeError:
   496                                                       raise AttributeError('You must use counting_mode=True to use compute_optimal_threshold method.')
   497                                           
   498         2          0.0      0.0      0.1          noise_value = CIC_noise**2+dark_noise**2+Sky_noise**2
   499                                                   
   500         2          0.0      0.0      1.1          gains=np.linspace(800,2500,n)#self.len_xaxis)
   501         2          0.0      0.0      1.0          rons=np.linspace(30,120,n)#self.len_xaxis)
   502         2          0.0      0.0      0.9          fluxes=np.linspace(0.01,0.7,n)#self.len_xaxis)
   503         2          0.0      0.0      0.9          smearings=np.linspace(0,2,n)#self.len_xaxis)
   504         2          0.0      0.0      0.8          noise=np.linspace(0.002,0.05,n)#self.len_xaxis)
   505         2          0.0      0.0      0.0          if (n==6)|(n==10):
   506         2          0.0      0.0      0.0              coords = (gains, rons, fluxes, smearings)
   507         2          0.0      0.0      0.0              point = (EM_gain, RN, flux, self.smearing)            
   508                                                   elif n==5:
   509                                                       coords = (gains, rons, fluxes, smearings,noise)
   510                                                       point = (EM_gain, RN, flux, self.smearing,noise_value)
   511                                                   else:
   512                                                       print(n,EM_gain, RN, flux, self.smearing,noise_value)
   513                                                       
   514         2          0.0      0.0      0.2          if ~np.isscalar(noise_value) |  ~np.isscalar(self.smearing) | ~np.isscalar(EM_gain) | ~np.isscalar(RN):
   515         2          0.0      0.0      9.3              point = np.repeat(np.zeros((4,1)), self.len_xaxis, axis=1).T
   516         2          0.0      0.0      0.1              point[:,0] =  self.EM_gain
   517         2          0.0      0.0      0.0              point[:,1] = self.RN
   518         2          0.0      0.0      0.0              point[:,2] = flux
   519         2          0.0      0.0      0.0              point[:,3] = self.smearing
   520         2          0.0      0.0     26.0          fraction_rn =interpn(coords, table_fraction_rn, point,bounds_error=False,fill_value=None)
   521         2          0.0      0.0     19.5          fraction_signal =interpn(coords, table_fraction_flux, point,bounds_error=False,fill_value=None)
   522         2          0.0      0.0     20.9          threshold = interpn(coords, table_threshold, point,bounds_error=False,fill_value=None)
   523         2          0.0      0.0     18.8          snr_ratio = interpn(coords, table_snr, point,bounds_error=False,fill_value=None)
   524                                           
   525         2          0.0      0.0      0.0          if type(self.smearing)==float:
   526         2          0.0      0.0      0.0              if self.smearing == 0:
   527                                                           a = Table.read("fraction_flux.csv")
   528                                                           threshold = 5.5
   529                                                           fraction_signal = np.interp(self.EM_gain/self.RN,a["G/RN"],a["fractionflux"])
   530                                                       # fraction_rn = f(flux=0.1,EM_gain=self.EM_gain, RN=self.RN)
   531                                                       # fraction_signal = f2(flux=0.1,EM_gain=self.EM_gain, RN=self.RN)
   532                                                       # snr_ratio = f3(flux=0.1,EM_gain=self.EM_gain, RN=self.RN)
   533                                           
   534         2          0.0      0.0      0.0          return threshold, fraction_signal, fraction_rn, snr_ratio#np.nanmax(SNR1/SNR_analogic)