from ipywidgets import interactive, interact, interact_manual

N = 512; L = 24
s = columnsfmri.simulation(N,L)
gwn = sim.gwnoise()
rho,deltaRelative = 0.5,0.5
fwhm = 2
beta = 0.05
w = 1

def f(rho,deltaRelative,fwhm,w):
    columnPattern = sim.columnPattern(rho,deltaRelative,gwn)
    boldPattern,_,_ = sim.bold(fwhm,beta,columnPattern)
    mriPattern = sim.mri(w,1+boldPattern)
    sim.plotColumnsBoldMRI(columnPattern,boldPattern,mriPattern)
    
interact(f,rho=[0.1,0.2,0.4,0.8,1.6,3.2],
      deltaRelative=[0.01, 0.25, 0.5, 0.75, 1],
      fwhm=[0.01, 0.5, 1, 1.5, 2,2.5, 3, 3.5],
      w=[0.25,0.5,1,1.5,2,3]);
      
# Show column pattern, bold response and fMRI image next to each other.     
sim.plotColumnsBoldMRI(columnPattern,boldPattern,mriPattern)

s# Calculate CNR and plot as a function of voxel width.
wRange = parameters['wRange']
cnr = cr/noiseOfW
cnrMean = np.mean(cnr,axis=0)
cnrStd = np.std(cnr,axis=0)
cnrErrorPlus = cnrMean+2*cnrStd
cnrErrorMinus = cnrMean-2*cnrStd

plt.plot(wRange, cnrMean)
plt.fill_between(wRange, cnrErrorMinus, cnrErrorPlus,alpha=0.5)
plt.xlabel('voxel width [mm]')
plt.ylabel('CNR')
plt.show()


import importlib
importlib.reload(my_module)

    def plotColumnsBoldMRI(self,columns,bold,mri):
        fig = plt.figure()
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
        minx = min(self.x)
        maxx = max(self.x)
        extent = [minx,maxx,minx,maxx]
        im1 = ax1.imshow(columns,'gray',
                       extent=extent,
                       interpolation='bilinear')
        im2 = ax2.imshow(bold,'gray',
                       extent=extent,
                       interpolation='bilinear')
        im3 = ax3.imshow(mri,'gray',
                       extent=extent,
                       interpolation='none')
        plt.show()
        
        
    def plotPatterns(self,plotList):
        minx = min(self.x)
        maxx = max(self.x)
        nPlots = len(plotPlist)
        fig,axes = plt.subplots(1,nPlots,
                                sharex='col',sharey='row',
                                figsize=(nPlots*3+1,4))
        fig.subplots_adjust(wspace=0.3)
        for ax in axes:
            # CONTINUE HERE:
            im = ax.imshow(y,cmap,
                           extent=[minx,maxx,minx,maxx],
                           interpolation='none')
            ax.set_title(title)
            ax.set_xlabel('position [mm]')
        fig.colorbar(im, ax=ax)
        plt.show()
        
        
        
def printResults(results):
    print('optimized quantity                 | optimal value | optimal voxel width')
    print('-----------------------------------+---------------+--------------------')
    print('univariate detection probability   | %.2f          | %.2f mm            ' %
          (results['pDetectUnivariate_max'],results['pDetectUnivariate_optW']))
    print('multivariate detection probability | %.2f          | %.2f mm            ' %
          (results['pDetectMultivariate_max'],results['pDetectMultivariate_optW']))
    print('decoding probability - 2 classes   | %.2f          | %.2f mm            ' %
          (results['accuracyDecode2_max'],results['pDecode2_optW']))
    print('decoding accuracy    - 2 classes   | %.2f          | %.2f mm            ' %
          (results['pDecode2_max'],results['accuracyDecode2_optW']))
    print('decoding probability - 4 classes   | %.2f          | %.2f mm            ' %
          (results['pDecode4_max'],results['pDecode4_optW']))
    print('decoding accuracy    - 4 classes   | %.2f          | %.2f mm            ' %
          (results['accuracyDecode4_max'],results['accuracyDecode4_optW']))
    print('decoding probability - 8 classes   | %.2f          | %.2f mm            ' %
          (results['pDecode8_max'],results['pDecode8_optW']))
    print('decoding accuracy    - 8 classes   | %.2f          | %.2f mm            ' %
          (results['accuracyDecode8_max'],results['accuracyDecode8_optW']))
    print('pattern correlation                | %.2f          | %.2f mm            ' %
          (results['patternCorrelation_max'],results['patternCorrelation_optW']))