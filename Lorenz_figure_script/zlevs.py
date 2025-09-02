import numpy as np
from netCDF4 import Dataset
import os
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from itertools import chain
from scipy.ndimage import gaussian_filter
import scipy.io as sio

def csf(sc, theta_s,theta_b):
    one64 = np.float64(1)
    
    if theta_s > 0.:
        csrf = ((one64 - np.cosh(theta_s * sc))/ (np.cosh(theta_s) - one64))
    else:
        csrf = -sc ** 2
    sc1 = csrf + one64

    if theta_b > 0.:
        Cs = ((np.exp(theta_b * sc1) - one64)/ (np.exp(theta_b) - one64) - one64)
    else:
        Cs = csrf
    return Cs

def zlevs(h,zeta,theta_s,theta_b,hc,N,ptype,vtransform):

    [M,L]=h.shape
    Ndim=N-1
    sc_r=np.zeros((N))
    Cs_r=np.zeros((N))
    sc_w=np.zeros((N+1))
    Cs_w=np.zeros((N+1))

    if (vtransform == 2):
        ds=1/N
     
        if 'w' in ptype:
            sc_w[0] = -1.0
            sc_w[N] =  0
            Cs_w[0] = -1.0
            Cs_w[N] =  0
          
            sc_w[1:-1]=(ds*(np.linspace(1,N-1,N-1)-N))
            Cs_w=csf(sc_w, theta_s,theta_b)
            N=N+1
        elif 'r' in ptype:
            sc= ds*(np.linspace(1,N,N)-N-0.5)
            Cs_r=csf(sc, theta_s,theta_b)
            sc_r=sc


    elif (vtransform==1):
        
        cff1=1/np.sinh(theta_s)
        cff2=0.5/np.tanh(0.5*theta_s)

        if 'w' in ptype:
            sc= (np.linspace(0,N,N+1)-N)/N
            N=N+1
        elif 'r' in ptype:
            sc=(np.linspace(1,N,N)-N-0.5)/N


        Cs=(1-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)


    Dcrit=0.2
    h[h==0]=1.e-14
    zeta[zeta<(Dcrit-h)]=Dcrit-h[zeta<(Dcrit-h)]
    hinv=1/h
    z  = np.empty((int(N),) + h.shape, dtype=np.float64)
    
    if (vtransform == 2):
        if 'w' in ptype:
            cff1=Cs_w
            cff2=sc_w+1
            sc=sc_w
        elif 'r' in ptype:
            cff1=Cs_r
            cff2=sc_r+1
            sc=sc_r

        hinv = 1. / (h + hc)
        cff = hc * sc

        for k in np.arange(N, dtype=int):
            z[k] = zeta + (zeta + h) * (cff[k] + cff1[k] * h) * hinv


    elif (vtransform==1):
        cff1=Cs
        cff2=sc+1
        cff=hc*(sc-Cs)
        cff2=sc+1

        for k in np.arange(N, dtype=int):
            z0=cff[k]+cff1[k]*h;
            z[k,:,:]=z0+zeta*(1.+z0*hinv);

    return z

