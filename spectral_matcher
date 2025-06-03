import numpy as np
from scipy.fft import fft,fftshift,ifft,ifftshift

"""
Spectral matching algorithm used to generate statistically equivalent fluctuations to input experimental data.
"""
class SpectralMatcher:
    """
    A spectral matching algorithm class.
    Contains methods to generate the psd of a signal, introduce randomness, and recover a
    statistically equivalent signal.
    Inputs:
        exp_data - experimental data sequence of melt pool measurements (width, depth, or cap) uniformly spaced by pw. assumed to be in microns (um)
        pw - pixel width (or resolution) of measurement sequence spacing. assumed to be in microns (um)
    """
    def __init__(self,exp_data,pw):
        self.exp_data,self.pw = exp_data,pw
    
    def generate_surface(self,return_all=False):
        # generates a 1D surface with the same power spectral density as inputted experimental data
        # this is the key function in this object class
        m = len(self.exp_data)
        if m%2!=0:
            self.exp_data = self.exp_data[:-1]
            m += -1

        # take fft and shift
        B = fftshift(fft(self.exp_data))

        # estimate power
        C = (self.pw/(2*np.pi*m))*np.power((np.abs(B)),2)



        # generate frequency vector for plotting
        q = np.zeros(m)
        for k in range(len(q)):
            q[k]=(2*np.pi/m)*k
        q_2 = fftshift(q)
        q_3 = np.unwrap(q_2-2*np.pi)
        q = q_3/self.pw    


        # sqrt to get fourier magnitude
        B = np.sqrt(C/(self.pw/(2*np.pi*m)))

        # conjugate symmetry operation
        Bq = B
        Bq[0] = 0
        Bq[m//2] = 0
        Bq[1:m//2] = np.flip(Bq[m//2+1:])


        # generate random phases
        phi = np.random.uniform(size=len(Bq))*2*np.pi-np.pi

        # apply conjugate symmetry
        phi[0] = 0
        phi[m//2] = 0
        phi[1:m//2] = -np.flip(phi[m//2+1:])

        # helper function polar->cartesian
        def pol2cart(rho, phi):
            x = rho * np.cos(phi)
            y = rho * np.sin(phi)
            return [x,y]

        # convert to complex coords
        [a,b] = pol2cart(Bq,phi)
        Hm = a+1j*b
        z = ifft(ifftshift(Hm))

        # take posterior PSD for verification
        B = fftshift(fft(z))

        # estimate power
        Cf = (self.pw/(2*np.pi*(m)))*np.power((np.abs(B)),2)
        
        if return_all:
            # returns initial and final PSDs, frequency plotting vector q, and equivalent signal z
            return C,Cf,q,z
        else:
            return z
    
    
    
    def gen_equivalent_fluctuation(self,return_all=False):
        # returns equivalent fluctuation with identical pixel width as input.
        if return_all:
            C,Cf,q,z_prime = self.generate_surface(return_all=return_all)
            return C,Cf,q,np.real(z_prime)+np.mean(self.exp_data)
        else: 
            z = self.generate_surface(self.exp_data,self.pw)
            return np.real(z)+np.mean(self.exp_data)