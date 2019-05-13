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