import numpy as np
from scipy.fft import fft,fftshift,ifft,ifftshift
import matplotlib.pyplot as plt

"""
TEST CASES


1. Synthetic data consisting of prescribed superposition of cosines with different amplitude, frequency and phase

2. Synthetic data consisting of (fixed seed) weakly stationary random variable

3. Data used in the analysis of stochastic defects in IN718

Outputs will show PSD of the input sequence and generated sequence to ensure equivalence, and an example of a generated track.
Note that each time the test cases are run, different output sequences will be computed.
"""

# import SpectralMatcher class
from spectral_matcher import SpectralMatcher

# TEST CASE 1: SYNTHETIC DATA
L = 1000 # length of time series
# synthetic fluctuation sequences follow
a1,t1,phi1 = 2,100,np.pi/4
a2,t2,phi2 = 4,500,3*np.pi/4
a3,t3,phi3 = 0.5,750,7*np.pi/4

sig1 = a1*np.cos(np.arange(L)/t1*2*np.pi +phi1)
sig2 = a2*np.cos(np.arange(L)/t2*2*np.pi +phi2)
sig3 = a3*np.cos(np.arange(L)/t3*2*np.pi +phi3)

synthetic_data = sig1+sig2+sig3

fig1,ax1 = plt.subplots(2,2,figsize=(10,5))

ax1[0,0].plot(sig1)
ax1[0,0].plot(sig2)
ax1[0,0].plot(sig3)
ax1[0,0].set_xlabel('t')
ax1[0,0].set_title('Synthetic data components')

ax1[0,1].plot(synthetic_data,c='r')
ax1[0,1].set_xlabel('t')
ax1[0,1].set_title('Synthetic data superposition')


# assume uniform sampling/spacing of 1 unit, i.e. pw=1
synthetic_pw = 1
sm = SpectralMatcher(synthetic_data,synthetic_pw)

# generating statistically equivalent fluctuation
C_syn,Cf_syn,q_syn,z_syn = sm.gen_equivalent_fluctuation(return_all=True)

ax1[1,0].loglog(q_syn,C_syn,label='synthetic input PSD',c='r')
ax1[1,0].loglog(q_syn,Cf_syn,label='generated sequence PSD',c='b',ls='--')
ax1[1,0].set_ylim(1e-7,1e4)
ax1[1,0].legend()
ax1[1,0].set_xlabel('q')
ax1[1,0].set_ylabel('C')

ax1[1,1].plot(synthetic_data,label='synthetic input',c='r')
ax1[1,1].plot(z_syn,label='generated sequence',c='b',ls='--')
ax1[1,1].legend()
ax1[1,1].set_xlabel('t')


plt.tight_layout()
plt.show()

# TEST CASE 2: SYNTHETIC RANDOM VARIABLE
N = 50 # number of random cosines to superposition
np.random.seed(271) # set seed for reproduction
# generating amplitudes, periods, phases for random cosines
a_i = np.random.uniform(low=0.1,high=5,size=N)
t_i = np.random.uniform(low=50,high=950,size=N)
phi_i = np.random.uniform(low=0,high=2*np.pi,size=N)
sigs = [ai*np.cos(2*np.pi*np.arange(L)/ti + phii) for ai,ti,phii in zip(a_i,t_i,phi_i)]
synthetic_random = np.sum(np.array(sigs),axis=0)

sm = SpectralMatcher(synthetic_random,synthetic_pw)
C_syn,Cf_syn,q_syn,z_syn = sm.gen_equivalent_fluctuation(return_all=True)
fig2,ax2 = plt.subplots(1,2,figsize=(10,3))
ax2[0].loglog(q_syn,C_syn,c='r',label='synthetic RV PSD')
ax2[0].loglog(q_syn,Cf_syn,c='b',ls='--',label='generated sequence PSD')
ax2[0].legend()
ax2[0].set_xlabel('q')
ax2[0].set_ylabel('C')
ax2[0].set_ylim(1e-6,1e5)

ax2[1].plot(synthetic_random,label='synthetic RV')
ax2[1].plot(z_syn,label='generated sequence')
ax2[1].set_xlabel('t')
ax2[1].legend()

plt.tight_layout()
plt.show()

# TEST CASE 3: EXTRACTED WIDTH, DEPTH DATA
width_data = np.loadtxt('widths_weaver2022.csv')
width_pw = 3.03 # resolution of width measurements
x_w = np.arange(len(width_data))*width_pw # distance along width scan in microns (um)
depth_data = np.loadtxt('depths_nadammal2021.csv')
depth_pw = 0.338 # resolution of depth measurements
x_d = np.arange(len(depth_data))*depth_pw



sm_w = SpectralMatcher(width_data,width_pw)
sm_d = SpectralMatcher(depth_data,depth_pw)
C_w,Cf_w,q_w,z_w = sm_w.gen_equivalent_fluctuation(return_all=True)
C_d,Cf_d,q_d,z_d = sm_d.gen_equivalent_fluctuation(return_all=True)


fig3,ax3 = plt.subplots(2,2,figsize=(10,5))
ax3[0,0].loglog(q_w,C_w,label='input width PSD',c='r')
ax3[0,0].loglog(q_w,Cf_w,label='generated width PSD',c='b',ls='--')
ax3[0,0].set_xlabel('q')
ax3[0,0].set_xlabel('c')
ax3[0,0].set_ylim(1e-6,1e5)
ax3[0,0].legend()


ax3[1,0].loglog(q_d,C_d,label='input depth PSD',c='r')
ax3[1,0].loglog(q_d,Cf_d,label='generated depth PSD',c='b',ls='--')
ax3[1,0].set_xlabel('q')
ax3[1,0].set_xlabel('c')
ax3[1,0].set_ylim(1e-8,1e8)
ax3[1,0].legend()

ax3[0,1].plot(x_w,width_data,label='input width',c='r')
ax3[0,1].plot(x_w[:-1],z_w,label='generated width',c='b')
ax3[0,1].set_xlabel('x')
ax3[0,1].legend()

ax3[1,1].plot(x_d,depth_data,label='input depth',c='r')
ax3[1,1].plot(x_d[:-1],z_d,label='generated depth',c='b')
ax3[1,1].set_xlabel('x')
ax3[1,1].legend()

plt.tight_layout()
plt.show()







