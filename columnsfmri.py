import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

def meanpower(s):
    return np.mean(np.abs(s**2))

def noiseModel(V,TR,nT,diffFlag,*args,**kwargs):
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
        neuronal = self.ift2(FORIENT*noiseF)
        return neuronal
    
    def bold(self,fwhm,beta,y):
        if fwhm==0:
            by = beta * y
            psf = None
            MTF = np.ones(np.size(y))
        else:
            fwhmfactor = 2*np.sqrt(2*np.log(2))
            psf = beta * \
            norm.pdf(self.x1,loc=0,scale=fwhm/fwhmfactor) * \
            norm.pdf(self.x2,loc=0,scale=fwhm/fwhmfactor)
            MTF = beta * np.exp(-(self.k1**2+self.k2**2) * 
                                2*(fwhm/fwhmfactor)**2*np.pi**2);
            by = self.ift2(MTF*self.ft2(y))
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
        
    def plotPattern(self,y):
        fig, ax = plt.subplots()
        im = ax.imshow(np.real(y),
                       extent=[min(self.x),max(self.x),min(self.x),max(self.x)],
                       interpolation='bilinear')
        fig.colorbar(im, ax=ax)
        plt.show()
    
    def plotVoxels(self,y):
        fig, ax = plt.subplots()
        im = ax.imshow(np.real(y),
                       extent=[min(self.x),max(self.x),min(self.x),max(self.x)],
                       interpolation='none')
        fig.colorbar(im, ax=ax)
        plt.show()