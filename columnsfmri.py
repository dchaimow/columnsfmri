import numpy as np
from scipy.stats import norm, chi2, gamma
import matplotlib.pyplot as plt
import seaborn as sns

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
        k = 6.6567;
        l = 0.0129;
        T1 = 1.607;
    if noiseType =='7T':
        k = 9.9632;
        l = 0.0113;
        T1 = 1.939;
    if noiseType =='Thermal':
        k = 9.9632;
        l = 0;
        T1 = 1.939;
    if noiseType == 'Physiological':
        k = np.Inf; # TEST THIS!
        l = 0.0113;
        T1 = 1.939;
    
    if not(differentialFlag) and nT != 1:
        raise ValueError('for multiple measurements only differential implemented!')
    elif nT == 1:
        F = np.sqrt(np.tanh(TR/(2*T1))/np.tanh(TR0/(2*T1)))
        k = k * F
        sigma = np.sqrt(1+l**2*k**2*V**2)/(k*V)
        return sigma
    
    s = 0
    assert(nT%2==0)
    for t1 in range(1,int(nT/2)+1):
        for t2 in range(1,int(nT/2)+1):
            s = s + np.exp((-TR*abs(t1-t2))/15)
    
    F = np.sqrt(np.tanh(TR/(2*T1))/np.tanh(TR0/(2*T1)))
    k = k * F
    sigma = np.sqrt((4/(k**2*V**2*nT)) + ((2*l**2)/(nT/2)**2)*s)
    return sigma

def detectionProbability(cnr,N):
    if np.size(cnr)>1:
        p = np.zeros(shape(cnr))
        if np.size(N)==1:
            N = np.ones(shape(cnr)) * N
        for z in range(size(cnr)):
            p[z] = detectionProbability(cnr[z],N[z])
    else:
        alpha = 0.05
        x_crit = chi2.ppf(1-alpha, df=N)
        a = N/2
        b = 2*(1+cnr**2)
        p = 1-gamma.cdf(x_crit,a,scale=b)
    return p

class sim:
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
        return np.random.randn(self.N,self.N)
    
    def ft2(self,y):
        return (self.L**2/self.N**2)*np.fft.fftshift(
            np.fft.fft2(np.fft.ifftshift(y)))
    
    def ift2(self,fy):
        return (self.N**2 * self.dk**2)*np.fft.fftshift(
            np.fft.ifft2(np.fft.ifftshift(fy)))
    
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
        neuronal = np.real(self.ift2(FORIENT*noiseF))
        return neuronal
    
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
                                2*(fwhm/fwhmfactor)**2*np.pi**2);
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
        
    def plotPattern(self,y):
        fig, ax = plt.subplots()
        minx = min(self.x)
        maxx = max(self.x)
        im = ax.imshow(y,'gray',
                       extent=[minx,maxx,minx,maxx],
                       interpolation='bilinear')
        fig.colorbar(im, ax=ax)
        plt.show()
    
    def plotVoxels(self,y):
        fig, ax = plt.subplots()
        minx = min(self.x)
        maxx = max(self.x)
        im = ax.imshow(y,'gray',
                       extent=[minx,maxx,minx,maxx],
                       interpolation='none')
        fig.colorbar(im, ax=ax)
        plt.show()
    
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
    