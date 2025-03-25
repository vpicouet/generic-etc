Timer unit: 0.1 s

Total time: 0.341224 s
File: /Users/Vincent/Github/fireball2-etc/notebooks/Observation.py
Function: PlotNoise at line 280

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   280                                               def PlotNoise(self,title='',x='exposure_time', lw=8):
   281                                                   """
   282                                                   Generate a plot of the evolution of the noise budget with one parameter:
   283                                                   exposure_time, Sky_CU, acquisition_time, Signal, EM_gain, RN, CIC_charge, Dard_current, readout_time, smearing, temperature, PSF_RMS_det, PSF_RMS_mask, QE, extra_background, cosmic_ray_loss_per_sec
   284                                                   """
   285         1          0.9      0.9     25.5          fig, axes= plt.subplots(4, 1, figsize=(12, 8), sharex=True) # fig, (ax1, ax2,ax3) = plt.subplots(3, 1, figsize=(12, 7), sharex=True) #figsize=(9, 5.5)
   286         1          0.0      0.0      0.0          ax1, ax2,ax3, ax4  = axes
   287         1          0.0      0.0      0.0          labels = ['%s: %0.3f (%0.1f%%)'%(name,self.electrons_per_pix[self.i,j],100*self.electrons_per_pix[self.i,j]/np.nansum(self.electrons_per_pix[self.i,:])) for j,name in enumerate(self.names)]
   288                                           
   289                                                   # ax1 
   290         7          0.0      0.0      0.0          for i,(name,c) in enumerate(zip(self.names,self.colors)):
   291         6          0.0      0.0      1.1              ax1.plot(getattr(self,x), self.noises[:,i]/self.factor,label='%s: %0.2f (%0.1f%%)'%(name,self.noises[self.i,i]/self.factor[self.i],self.percents[i,self.i]),lw=lw,alpha=0.8,c=c)
   292         1          0.0      0.0      0.2          ax1.plot(getattr(self,x), np.nansum(self.noises[:,:-1],axis=1)/self.factor,label='%s: %0.2f (%0.1f%%)'%("Total",np.nansum(self.noises[self.i,-1])/self.factor[self.i],np.nansum(self.percents[:,self.i])),lw=lw,alpha=0.4,c="k")
   293         1          0.1      0.1      2.1          ax1.legend(loc='upper right')
   294         1          0.0      0.0      0.0          ax1.set_ylabel('Noise (e-/pix/exp)')
   295                                           
   296                                                   # ax1b = ax1.secondary_yaxis("right", functions=( lambda x:  x * self.factor[self.i], lambda x:x / self.factor[self.i] ))
   297                                                   # self.ax1b = ax1b
   298                                                   # ax1b.set_ylabel("Noise (e-/res/N frames)")#r"%0.1f,%0.1f,%0.1f"%(self.factor[self.i],self.resolution_element , np.sqrt(self.N_resol_element_A)))
   299                                           
   300                                           
   301                                                   # ax2 
   302         1          0.0      0.0      0.2          ax2.grid(False)
   303         1          0.1      0.1      2.3          self.stackplot1 = ax2.stackplot(getattr(self,x),  np.array(self.electrons_per_pix).T[:,:],alpha=0.7,colors=self.colors,labels=labels)
   304         1          0.0      0.0      0.0          ax2.set_ylabel('e-/pix/frame')
   305         1          0.0      0.0      1.2          ax2.legend(loc='upper right',title="Overall background: %0.3f (%0.1f%%)"%(np.nansum(self.electrons_per_pix[self.i,1:]),100*np.nansum(self.electrons_per_pix[self.i,1:])/np.nansum(self.electrons_per_pix[self.i,:])))
   306         1          0.0      0.0      0.5          ax2.set_xlim((getattr(self,x).min(),getattr(self,x).max()))
   307                                           
   308                                           
   309                                                   # ax2b = ax2.secondary_yaxis("right", functions=( lambda x:  x * self.factor[self.i]**2, lambda x:x / self.factor[self.i]**2 ))
   310                                                   # self.ax2b = ax2b
   311                                                   # ax2b.set_ylabel(r"%0.1f,%0.1f,%0.1f"%(self.factor[self.i],self.resolution_element , np.sqrt(self.N_resol_element_A)))
   312                                           
   313                                           
   314                                           
   315                                                   # ax3
   316         1          0.0      0.0      0.2          ax3.grid(False)
   317         1          0.1      0.1      2.2          self.stackplot2 = ax3.stackplot(getattr(self,x), self.snrs * np.array(self.noises).T[:-1,:]**2/self.Total_noise_final**2,alpha=0.7,colors=self.colors)
   318         1          0.0      0.0      0.1          ax3.set_ylim((0,np.nanmax(self.SNR)))
   319         1          0.0      0.0      0.0          ax3.set_ylabel('SNR (res, N frames)')        
   320                                           
   321                                                   # ax3b = ax3.secondary_yaxis("right", functions=( lambda x: x / self.factor[self.i]**2, lambda x: x * self.factor[self.i]**2))
   322                                                   # self.ax3b = ax3b
   323                                                   # ax3b.set_ylabel(r" SNR(res/N frames")
   324                                           
   325                                           
   326                                           
   327                                                   # ax4
   328         1          0.0      0.0      0.2          ax4.plot(getattr(self,x), np.log10(self.extended_source_5s),"-",lw=lw-1,label="SNR=5 Flux/Pow on one elem resolution (%0.2f-%0.2f)"%(np.log10(self.point_source_5s[self.i]),np.nanmin(np.log10(self.point_source_5s))),c="k")
   329                                                   # if self.instrument==FIREBall:
   330         1          0.0      0.0      0.0          if "FIREBall" in self.instrument:
   331                                           
   332         1          0.0      0.0      0.2              ax4.plot(getattr(self,x), np.log10(self.extended_source_5s/np.sqrt(2)),"-",lw=lw-1,label="Two elem resolution (%0.2f-%0.2f)"%(np.log10(self.point_source_5s[self.i]/np.sqrt(2)),np.nanmin(np.log10(self.point_source_5s/np.sqrt(2)))),c="grey")
   333                                                       # ax4.plot(getattr(self,x), np.log10(self.extended_source_5s/np.sqrt(40)),"-",lw=lw-1,label="20 sources stacked on 2 res elem. (%0.2f-%0.2f)"%(np.log10(self.point_source_5s[self.i]/np.sqrt(40)),np.nanmin(np.log10(self.point_source_5s/np.sqrt(40)))),c="lightgrey")
   334                                                       # ax4.plot(getattr(self,x), np.log10(self.extended_source_5s/np.sqrt(2)/30),"-",lw=lw-1,label="Sources transported to high z: (%0.2f-%0.2f) \ngain of factor 22-50 depending on line resolution"%(np.log10(self.point_source_5s[self.i]/np.sqrt(2)/30),np.nanmin(np.log10(self.point_source_5s/np.sqrt(2)/30))),c="whitesmoke")
   335         1          0.0      0.0      0.0          T2 =  lambda x:np.log10(10**x/1.30e57)
   336         1          0.0      0.0      0.0          self.pow_2018 = 42.95
   337         1          0.0      0.0      0.0          self.pow_best = 41.74
   338         1          0.2      0.2      5.5          ax4b = ax4.secondary_yaxis("right", functions=(lambda x:np.log10(10**x * 1.30e57),T2))
   339         1          0.0      0.0      0.0          if ("FIREBall" in self.instrument) & (1==0):
   340                                                       ax4.plot([getattr(self,x).min(),getattr(self,x).min(),np.nan,getattr(self,x).max(),getattr(self,x).max()],[T2(self.pow_2018),T2(self.pow_best),np.nan,T2(self.pow_2018),T2(self.pow_best)],lw=lw,label="2018 flight (%0.1f) - most optimistic case (%0.1f)"%(self.pow_2018,self.pow_best),c="r",alpha=0.5)
   341         1          0.0      0.0      0.0          self.T2=T2
   342         1          0.0      0.0      0.0          self.ax4b = ax4b
   343         1          0.0      0.0      0.8          ax4.legend(loc="upper right", fontsize=8,title="Left: Extend. source F, Right: Point source power" )
   344         1          0.0      0.0      0.0          ax4.set_ylabel(r"Log(erg/cm$^2$/s/asec$^2$)")
   345         1          0.0      0.0      0.0          ax4b.set_ylabel(r" Log(erg/s)")
   346                                           
   347         1          0.0      0.0      0.0          axes[-1].set_xlabel(x)
   348         1          0.0      0.0      0.2          ax1.tick_params(labelright=True,right=True)
   349         1          0.0      0.0      0.2          ax2.tick_params(labelright=True,right=True)
   350         1          0.0      0.0      0.2          ax3.tick_params(labelright=True,right=True)
   351         1          1.9      1.9     57.1          fig.tight_layout(h_pad=0.01)
   352         1          0.0      0.0      0.0          return fig