import numpy as np
from scipy.stats import norm, chi2, gamma, weibull_min
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def displayFigureA(results):
    """
    figure A - optimization of all 3 objectives
    """
    fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharex='col',figsize=(10,3))
    fig.subplots_adjust(wspace=0.3)
    
    ax1.plot(results['w'],results['pDetectUnivariate'],'C0')
    ax1.plot([results['pDetectUnivariate_optW'],results['pDetectUnivariate_optW']],
             [0, 1],'C0--',)
    ax1.axis([min(results['w']),max(results['w']),0,1])
    ax1.set_xlabel('voxel width [mm]')
    ax1.set_ylabel('probability')
    ax1.set_title('Univariate detection')
        
    ax2.plot(results['w'],results['accuracyDecode2'],'C1')
    ax2.plot([results['accuracyDecode2_optW'],results['accuracyDecode2_optW']],
             [0.5, 1],'C1--',)
    ax2.axis([min(results['w']),max(results['w']),0.5,1])
    ax2.set_xlabel('voxel width [mm]')
    ax2.set_ylabel('accuracy')
    ax2.set_title('Multivariate decoding')
    
    ax3.plot(results['w'],results['patternCorrelation'],'C2')
    ax3.plot([results['patternCorrelation_optW'],results['patternCorrelation_optW']],
             [0, 1],'C2--',)
    ax3.axis([min(results['w']),max(results['w']),0,1])   
    ax3.set_xlabel('voxel width [mm]')
    ax3.set_ylabel('correlation')
    ax3.set_title('Reconstruction')

    plt.show()
    
def printResults(results):
    data = [['univariate detection probability',
             results['pDetectUnivariate_max'],
             results['pDetectUnivariate_optW']],
            ['multivariate detection probability',
             results['pDetectMultivariate_max'],
             results['pDetectUnivariate_optW']],
            ['decoding probability - 2 classes',
             results['pDecode2_max'],
             results['pDecode2_optW']],
            ['decoding accuracy    - 2 classes',
             results['accuracyDecode2_max'],
             results['accuracyDecode2_optW']],
            ['decoding probability - 4 classes',
             results['pDecode4_max'],
             results['pDecode4_optW']],
            ['decoding accuracy    - 4 classes',
             results['accuracyDecode4_max'],
             results['accuracyDecode4_optW']],
            ['decoding probability - 8 classes',
             results['pDecode8_max'],
             results['pDecode8_optW']],
            ['decoding accuracy    - 8 classes',
             results['accuracyDecode8_max'],
             results['accuracyDecode8_optW']],
            ['pattern correlation',
             results['patternCorrelation_max'],
             results['patternCorrelation_optW']]]

    table = pd.DataFrame(data,columns=['optimized quantity',
                                       'optimal value',
                                       'optimal voxel width'])
    return table

def cdfModelWeibull(a,b,y0):
    f = lambda x: y0+(weibull_min.cdf(x,b,scale=a)-weibull_min.cdf(0,b,scale=a))*(1-y0)
    return f

def simulatefMRIOfColumnPatterns(parameters):
    # set parameters
    nTrials        = parameters['nTrials']
    Nsim           = parameters['N']
    L              = parameters['L']
    rho            = parameters['rho']
    deltaRelative  = parameters['deltaRelative']
    fwhm           = parameters['fwhm']
    beta           = parameters['beta']
    wRange         = parameters['wRange']
    sliceThickness = parameters['sliceThickness']
    TR             = parameters['TR']
    nT             = parameters['nT']
    noiseType      = parameters['noiseType']
    AFlat          = parameters['AFlat']
    
    nWRange = len(wRange)
    
    # setup simulation
    sim = simulation(Nsim,L);

    # voxel volume as a function of voxel width
    voxelVOfW = sliceThickness*wRange**2;

    # number of voxels as a function of voxel width (assuming single slice)
    nVoxelsOfW = AFlat/(wRange**2);

    # noise level as a function of voxel volume
    differentialFlag = True;
    noiseOfW = noiseModel(voxelVOfW,TR,nT,differentialFlag,noiseType=noiseType);
    
    # initialize contrast range and pattern correlation result arrays
    cr = np.zeros((nTrials,nWRange));
    cor = np.zeros((nTrials,nWRange));
                   
    for zTrial in range(nTrials):
        np.random.seed(parameters['randomNumberSeed']+zTrial)
    
        # initialize noise pattern for simulation of columns
        gwn = sim.gwnoise()
        # simulate column pattern
        columnPattern,_ =  sim.columnPattern(rho,deltaRelative,gwn) 
        # simulate BOLD response pattern
        boldPattern,_,_ = sim.bold(fwhm,beta,columnPattern);
    
        # for each tested voxel width
        for zW in range(nWRange):
            w = wRange[zW];
            # simulate MRI voxel sampling
            mriPattern = sim.mri(w,1+boldPattern)
            # calculate contrast range
            cr[zTrial,zW] = np.std(mriPattern);
            # add noise
            mriPlusNoisePattern = mriPattern + noiseOfW[zW] * np.random.randn(*mriPattern.shape)
            # calculate correlation to original (neuronal) pattern        
            cor[zTrial,zW] = sim.patternCorrelation(columnPattern,mriPlusNoisePattern);
    
    # average over trials    
    cr = np.mean(cr,axis=0)
    cor = np.mean(cor,axis=0)

    # calculate CNR and mvCNR
    CNR = cr/noiseOfW
    a = 0.26
    mvCNR = CNR * nVoxelsOfW**a
    
    # calculate univariate detection probability
    pUniOfW = detectionProbability(CNR,1)
    
    # calculate multivariate detection probability
    fPMultiDetect = cdfModelWeibull(1.9038,3.5990,0.05)
    pMultiOfW = fPMultiDetect(mvCNR)
    
    # calculate decoding probability
    fPDecode2 = cdfModelWeibull(2.0675,3.4101,0.05)
    fPDecode4 = cdfModelWeibull(2.8451,3.4770,0.05)
    fPDecode8 = cdfModelWeibull(3.8086,3.2733,0.05)

    pDecodeOfW2 = fPDecode2(mvCNR)
    pDecodeOfW4 = fPDecode4(mvCNR)
    pDecodeOfW8 = fPDecode8(mvCNR)

    # calculate decoding accuracy
    fMeanClassPerf2 = cdfModelWeibull(4.7018,1.9954,1/2)
    fMeanClassPerf4 = cdfModelWeibull(8.0900,2.1585,1/4)
    fMeanClassPerf8 = cdfModelWeibull(15.8022,1.9954,1/8)

    meanClassPerfOfW2 = fMeanClassPerf2(mvCNR)
    meanClassPerfOfW4 = fMeanClassPerf4(mvCNR)
    meanClassPerfOfW8 = fMeanClassPerf8(mvCNR)
    
    # find optima
    idxMaxPUni = np.argmax(pUniOfW)
    maxPUni = pUniOfW[idxMaxPUni]
    optWPUni = wRange[idxMaxPUni]
    
    idxMaxMVCNR = np.argmax(mvCNR)
    maxMVCNR = mvCNR[idxMaxMVCNR]    
    optWMVCNR = wRange[idxMaxMVCNR]

    maxPMulti = np.amax(pMultiOfW)
    assert maxPMulti==pMultiOfW[idxMaxMVCNR]
    
    maxPDecode2 = np.max(pDecodeOfW2)
    maxPDecode4 = np.max(pDecodeOfW4)
    maxPDecode8 = np.max(pDecodeOfW8)
    assert maxPDecode2==pDecodeOfW2[idxMaxMVCNR]
    assert maxPDecode4==pDecodeOfW4[idxMaxMVCNR]
    assert maxPDecode8==pDecodeOfW8[idxMaxMVCNR]
    
    maxMeanClassPerf2 = np.max(meanClassPerfOfW2)
    maxMeanClassPerf4 = np.max(meanClassPerfOfW4)
    maxMeanClassPerf8 = np.max(meanClassPerfOfW8)
    assert maxMeanClassPerf2==meanClassPerfOfW2[idxMaxMVCNR]
    assert maxMeanClassPerf4==meanClassPerfOfW4[idxMaxMVCNR]
    assert maxMeanClassPerf8==meanClassPerfOfW8[idxMaxMVCNR]
    
    idxMaxCor = np.argmax(cor)
    maxCor = cor[idxMaxCor]
    optWCor = wRange[idxMaxCor]
    
    # assign results
    results = dict()
    
    results['w']                        = wRange
    results['cr']                       = cr
    results['noiseLevel']               = noiseOfW
    results['SNR']                      = 1/noiseOfW
    results['CNR']                      = CNR

    results['pDetectUnivariate']        = pUniOfW
    results['pDetectUnivariate_max']    = maxPUni
    results['pDetectUnivariate_optW']   = optWPUni

    results['mvCNR']                    = mvCNR
    results['mvCNR_max']                = maxMVCNR
    results['mvCNR_opt']                = optWMVCNR

    results['pDetectMultivariate']      = pMultiOfW
    results['pDetectMultivariate_max']  = maxPMulti
    results['pDetectMultivariate_optW'] = optWMVCNR

    results['pDecode2']                  = pDecodeOfW2
    results['pDecode2_max']              = maxPDecode2
    results['pDecode2_optW']             = optWMVCNR

    results['accuracyDecode2']           = meanClassPerfOfW2
    results['accuracyDecode2_max']       = maxMeanClassPerf2
    results['accuracyDecode2_optW']      = optWMVCNR

    results['pDecode4']                  = pDecodeOfW4
    results['pDecode4_max']              = maxPDecode4
    results['pDecode4_optW']             = optWMVCNR

    results['accuracyDecode4']           = meanClassPerfOfW4
    results['accuracyDecode4_max']       = maxMeanClassPerf4
    results['accuracyDecode4_optW']      = optWMVCNR

    results['pDecode8']                  = pDecodeOfW8
    results['pDecode8_max']              = maxPDecode8
    results['pDecode8_optW']             = optWMVCNR

    results['accuracyDecode8']           = meanClassPerfOfW8
    results['accuracyDecode8_max']       = maxMeanClassPerf8
    results['accuracyDecode8_optW']      = optWMVCNR

    results['patternCorrelation']       = cor
    results['patternCorrelation_max']   = maxCor
    results['patternCorrelation_optW']  = optWCor
    
    # visualize results
    # printResults(results)
    # displayFigureA(results)
    # displayFigureB(results)
    # displayFigureC(results)

    return results
    
def setParameters(*args):
    """
    Sets parameters for simulatefMRIOfColumnPatterns.
    
    parameters = setParameters() returns default parameters.

    parameters = setParameters(s1,...) sets parameters according to
    predefined scenarios.
    s1,... are strings that select one of multiple scanner/pulse
    sequence and/or pattern irregularity scenarios:
    '3TGE','7TGE','7TSE','regular','irregular 
 
    parameters is a dictionary consisting of the following entries:
 
    randomNumberSeed    random number seed
 
    nTrials             number of simulation trials
    N                   simulation grid points along one dimension
    L                   simulation size along one dimension [mm]
    wRange              list of MRI voxel widths
                        (need to be divisors of L)                          
    rho                 main pattern frequency 
                        (~1/(2*column spacing) [1/mm]
    deltaRelative       (relative) irregularity (= bandpass FWHM/rho)
 
    fwhm                BOLD PSF FWHM [mm]
    beta                BOLD PSF amplitude [relative signal ch.]
 
    sliceThickness      slice thickness [mm]
    AFlat               FOV area, determines number of voxels [mm^2]
    TR                  TR (repetition time)
    nT                  number of volumes (measurements), 
                        (nT/2 per one of two conditions)
    noiseType           noise type, either '3T' or '7T',
                        can either be a numeric vector
                        [k, l, T1] consisting of noise model
                        parameters kappa and lambda 
                        (see Triantafyllou et al. 2005) and
                        longitudinal relaxation time T1
 
    see also simulatefMRIOfColumnPatterns
    """
    parameters = dict()
    
    # seed for random number generator
    parameters['randomNumberSeed'] = 23 

    # simulation parameters
    parameters['nTrials'] = 32  # number of simulation trials
    parameters['N'] = 512  # simulation grid points
    parameters['L'] = 24  # mm, for simulation/calc. of contrast range

    # range of voxel widths
    parameters['wRange'] = parameters['L']/(2*np.hstack((np.arange(3,29),[30,32,34,37,40,44,48,53,60,69,80,96,120,160,240])))
    
    # pattern parameters
    parameters['rho'] = 1/1.6  # main pattern frequency = 1/(2*column spacing)
    parameters['deltaRelative'] = 0.5  # (moderate) irregularity (bandpass FWHM/rho)

    # PSF parameters
    parameters['fwhm'] = 1.02   # PSF FWHM in mm
    parameters['beta'] = 0.035  # amplitute

    # further MRI parameters affecting SNR calculation
    parameters['sliceThickness'] = 2.5  # mm
    parameters['AFlat'] = 87  # mm^2 area for calculating nVoxels (from odc rois)
    parameters['TR'] = 2  # s
    parameters['nT'] = 1000  # number of volumes (sum of both conditions)
    parameters['noiseType'] = '7T'  # alternatively: '3T'

    if '3TGE' in args:
        assert not('7TGE' in args) and not('7TSE' in args),\
        'Only one field strength/pulse sequence scenario allowed!'
        parameters['fwhm'] = 2.8  # PSF FWHM in mm
        parameters['beta'] = 0.03 # amplitute
        parameters['noiseType'] = '3T' # alternatively: '3T'

    if '7TGE' in args:
        assert not('3TGE' in args) and not('7TSE' in args),\
        'Only one field strength/pulse sequence scenario allowed!'      
        parameters['fwhm'] = 1.02   # PSF FWHM in mm
        parameters['beta'] = 0.035  # amplitute
        parameters['noiseType'] = '7T'  # alternatively: '3T'

    if '7TSE' in args:
        assert not('3TGE' in args) and not('7TGE' in args),\
        'Only one field strength/pulse sequence scenario allowed!'
        parameters['fwhm'] = 0.82   # PSF FWHM in mm
        parameters['beta'] = 0.025  # amplitute
        parameters['noiseType'] = '7T'  # alternatively: '3T'

    if 'irregular' in args:
        assert not('regular' in args),\
        'Only one regularity scenario allowed!'        
        parameters['deltaRelative'] = 1 

    if 'regular' in args:
        assert not('irregular' in args),\
        'Only one regularity scenario allowed!'        
        parameters['deltaRelative'] = 0 
        
    return parameters

def meanpower(s):
    return np.mean(np.abs(s**2))

def noiseModel(V,TR,nT,differentialFlag,*args,**kwargs):
    TR0 = 5.4
    
    noiseType = kwargs.get('noiseType', None)
    k = kwargs.get('k', None)
    l =  kwargs.get('l', None)
    T1 = kwargs.get('T1', None)
    
    if noiseType == None:
        if l == None or k == None or T1 == None:
            raise ValueError('k,l or T1 not specified!')
    else:
        if l != None or k != None or T1 != None:
            raise ValueError('specify either noiseType or (k,l and T1), not both!')
    if noiseType == '3T':
        k = 6.6567 
        l = 0.0129 
        T1 = 1.607 
    if noiseType =='7T':
        k = 9.9632 
        l = 0.0113 
        T1 = 1.939 
    if noiseType =='Thermal':
        k = 9.9632 
        l = 0 
        T1 = 1.939 
    if noiseType == 'Physiological':
        k = np.Inf  # TEST THIS!
        l = 0.0113 
        T1 = 1.939 
    
    if not(differentialFlag) and nT != 1:
        raise ValueError('for multiple measurements only differential implemented!')
    elif nT == 1:
        F = np.sqrt(np.tanh(TR/(2*T1))/np.tanh(TR0/(2*T1)))
        k = k * F
        sigma = np.sqrt(1+l**2*k**2*V**2)/(k*V)
        return sigma
    
    s = 0
    assert nT%2==0 
    for t1 in range(1,int(nT/2)+1):
        for t2 in range(1,int(nT/2)+1):
            s = s + np.exp((-TR*abs(t1-t2))/15)
    
    F = np.sqrt(np.tanh(TR/(2*T1))/np.tanh(TR0/(2*T1)))
    k = k * F
    sigma = np.sqrt((4/(k**2*V**2*nT)) + ((2*l**2)/(nT/2)**2)*s)
    return sigma

def detectionProbability(cnr,N):
    if np.size(cnr)>1:
        p = np.zeros(cnr.shape)
        if np.size(N)==1:
            N = np.ones(cnr.shape) * N
        for z in range(np.size(cnr)):
            p[z] = detectionProbability(cnr[z],N[z])
    else:
        alpha = 0.05
        x_crit = chi2.ppf(1-alpha, df=N)
        a = N/2
        b = 2*(1+cnr**2)
        p = 1-gamma.cdf(x_crit,a,scale=b)
    return p

class simulation:
    def __init__(self,N,L):
        self.N = N
        self.L = L
        self.dx = L/self.N
        self.dk = 1/self.L
        self.Fs = 1/self.dx
        
        self.x = np.fft.fftshift(
            np.concatenate(
                (np.arange(0,self.N/2),np.arange(-self.N/2,0))) * self.dx)
        self.k = np.fft.fftshift(
            np.concatenate(
                (np.arange(0,self.N/2),np.arange(-self.N/2,0))) * self.dk)
        self.x1, self.x2 = np.meshgrid(self.x, self.x)
        self.k1, self.k2 = np.meshgrid(self.k, self.k)
        
    def gwnoise(self):
        return np.random.randn(self.N,self.N) + 1j* np.random.randn(self.N,self.N)
    
    def ft2(self,y):
        return (self.L**2/self.N**2)*np.fft.fftshift(
            np.fft.fft2(np.fft.ifftshift(y)))
    
    def ift2(self,fy):
        return (self.N**2 * self.dk**2)*np.fft.ifftshift(
            np.fft.ifft2(np.fft.fftshift(fy)))
    
    def columnPattern(self,rho,deltaRelative,gwnoise):
        fwhmfactor = 2*np.sqrt(2*np.log(2))
        r = np.sqrt(self.k1**2+self.k2**2)
        if deltaRelative==0:
            FORIENTNotNormalized = np.double(abs(r - rho)<self.dk/2)
        else:
            FORIENTNotNormalized = \
            norm.pdf(r,loc= rho,scale=(rho*deltaRelative)/fwhmfactor) + \
            norm.pdf(r,loc=-rho,scale=(rho*deltaRelative)/fwhmfactor)
        C = (np.sqrt(meanpower(FORIENTNotNormalized)))*np.sqrt(np.pi/8)
        FORIENT = FORIENTNotNormalized/C
        noiseF = self.ft2(gwnoise)
        gamma = self.ift2(FORIENT*noiseF)
        neuronal = np.real(gamma)
        preferredOrientation = np.angle(gamma)
        return neuronal, preferredOrientation
    
    def bold(self,fwhm,beta,y):
        if fwhm==0:
            by = beta * y
            psf = None
            MTF = np.ones(np.shape(y))
        else:
            fwhmfactor = 2*np.sqrt(2*np.log(2))
            psf = beta * \
            norm.pdf(self.x1,loc=0,scale=fwhm/fwhmfactor) * \
            norm.pdf(self.x2,loc=0,scale=fwhm/fwhmfactor)
            MTF = beta * np.exp(-(self.k1**2+self.k2**2) * 
                                2*(fwhm/fwhmfactor)**2*np.pi**2) 
            by = np.real(self.ift2(MTF*self.ft2(y)))
        return by,psf,MTF
    
    def mri(self,w,y):
        nSamplesHalf = self.L/(2*w)
        if nSamplesHalf % 1 == 0:
            nSamplesHalf = int(nSamplesHalf)
            yk = self.ft2(y)
            centerIdx = int(self.N/2)
            downSampledY = yk[centerIdx-nSamplesHalf:centerIdx+nSamplesHalf,
                              centerIdx-nSamplesHalf:centerIdx+nSamplesHalf]
            my = abs(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(downSampledY))))/w**2
        else:
            my = None
        return my
    
    def upsample(self,y):
        Fy = np.fft.fft2(y)
        nx, ny = Fy.shape # size of original matrix
        nxAdd = self.N - nx # number of zeros to add
        nyAdd = self.N - ny
        # quadrants of the FFT, starting in the upper left
        q1 = Fy[:int(nx/2),:int(ny/2)]
        q2 = Fy[:int(nx/2),int(ny/2):]
        q3 = Fy[int(nx/2):,int(ny/2):]
        q4 = Fy[int(nx/2):,:int(ny/2)]
        zeroPaddRows = np.zeros((nxAdd,int(ny/2)))
        zeroPaddColumns = np.zeros((nxAdd+nx,nyAdd))
        zeroPaddFy = np.hstack(
            (np.vstack((q1,zeroPaddRows,q4)),
             zeroPaddColumns,
             np.vstack((q2,zeroPaddRows,q3))))
        upPattern = np.real(np.fft.ifft2(zeroPaddFy)) * self.N**2/(nx*ny)
        return upPattern
    
    def patternCorrelation(self,orig,mri):
        mriUp = self.upsample(mri)
        c = np.corrcoef(orig.flatten(),mriUp.flatten())
        r = c[0,1]
        return r
        
    def plotPattern(self,y,cmap='gray',title=None, ax=None):
        if not(ax):
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        minx = min(self.x)
        maxx = max(self.x)
        im = ax.imshow(y,cmap,
                       extent=[minx,maxx,minx,maxx],
                       interpolation='none')
        ax.set_title(title)
        ax.set_xlabel('position [mm]')
        fig.colorbar(im,ax=ax)
        #plt.show()

            
    