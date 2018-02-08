import numpy as np
from scipy.integrate import quad, dblquad, nquad
import scipy.interpolate as itp
from scipy import LowLevelCallable
from itertools import combinations_with_replacement as cwr
import itertools as itr
from _BFisher import lib, ffi
import ctypes
from ctypes import *

kk, PP = np.loadtxt('test_matterpower.dat',unpack=True)
Piso = itp.interp1d(kk, PP, 'cubic')

libb = ctypes.cdll.LoadLibrary('./BFisherutils.so')
libb.ICovB.restype = c_double
libb.ICovB.argtypes = (c_int, POINTER(c_double), c_void_p)
libb.ICrossPB.restype = c_double
libb.ICrossPB.argtypes = (c_int, POINTER(c_double), c_void_p) 

def CovP(mu,k,l1,l2,parc,pars):
    var = (k, mu)
    Pk1 = Piso(k)
    cov = lib.CovP(Pk1,var,parc,pars)*lib.Legandre(l1,mu)*lib.Legandre(l2,mu)
    return cov

def Pcov(k1,l1,k2,l2,parc,pars):
    cov = quad(CovP,0,1,args=(k1,l1,l2,parc,pars),epsrel=1e-3)[0]
    return cov

def Bcov(k1,k2,k3,l1,k4,k5,k6,l2,parc,pars):

    Pk1 = Piso(k1)
    Pk2 = Piso(k2)
    Pk3 = Piso(k3)
    apar, aper, f, b1, b2 = parc
    nave, Vs = pars
    udata = (c_double*12)(Pk1,Pk2,Pk3,apar,aper,f,b1,b2,nave,Vs)
    user_data = cast(pointer(udata),c_void_p)

    Ifunct = LowLevelCallable(libb.ICovB,user_data)
    option = {'epsrel' : 1e-3}
    arguments = (c_double*3)(k1,k2,k3)
    cov = nquad(Ifunct,[[0,1],[0,2*np.pi]],args=arguments,opts=option)[0]
    return cov

def PBcross(k1,l1,k2,k3,k4,l2,parc,pars):
    Pk1 = Piso(k1)
    Pk2 = Piso(k2)
    Pk3 = Piso(k3)
    apar, aper, f, b1, b2 = parc
    nave, Vs = pars
    udata = (c_double*12)(Pk1,Pk2,Pk3,apar,aper,f,b1,b2,nave,Vs)
    user_data = cast(pointer(udata),c_void_p)

    Ifunct = LowLevelCallable(libb.ICrossPB,user_data)
    option = {'epsrel' : 1e-3}
    arguments = (c_double*3)(k1,k2,k3)
    cov = nquad(Ifunct,[[0,1],[0,2*np.pi]],args=arguments,opts=option)[0]
    return cov

kmax = 0.2
Nk = 20
kedges = np.linspace(0,kmax,Nk+1)
klow = kedges[:-1]
khigh = kedges[1:]
kmid = (klow + khigh)/2

n = 1e-4
V = 1000

pars = (n, V)

apar = 1
aper = 1
f = 0.7
b1 = 2
b2 = 0.3

parc = (apar, aper, f, b1, b2)

multipoles = (0,2,4)
Pk_sequence = list(itr.product(multipoles,kmid))
Bk_sequence = list(itr.product(multipoles,cwr(kmid,3)))
print('NBk', len(Bk_sequence))
Bk_sequence = [k for k in Bk_sequence if k[1][0] + k[1][1] > k[1][2]]
NP = len(Pk_sequence)
NB = len(Bk_sequence)
print('NB',NB)
NPB = NP + NB
print('NPB',NPB)
PBcov = np.zeros((NPB,NPB))

print('PP covariance')
for i, pi in enumerate(Pk_sequence):
    l1, k1 = pi
    for j, pj in enumerate(Pk_sequence):
        l2, k2 = pj
        if k1 == k2:
            PBcov[i,j] = Pcov(k1,l1,k2,l2,parc,pars)
            PBcov[j,i] = PBcov[j,i]

print('BB covariance')
for i, bi in enumerate(Bk_sequence):
    print(i)
    l1, (k1, k2, k3) = bi
    for j, bj in enumerate(Bk_sequence):
        l2, (k4, k5, k6) = bj
        if k1 == k4 and k2 == k5 and k3 == k5:
            PBcov[NP+i,NP+j] = Bcov(k1,k2,k3,l1,k4,k5,k6,l2,parc,pars)
            PBcov[NP+j,NP+i] = PBcov[NP+i,NP+j]

print('PB cross-covariance')
for i, pi in enumerate(Pk_sequence):
    print(i)
    l1, k1 = pi
    for j, bj in enumerate(Bk_sequence):
        l2, (k2, k3, k4) = bj
        if k1 == k2 or k1 == k3 or k1 == k4:
            PBcov[i,NP+j] = PBcross(k1,l1,k2,k3,k4,l2,parc,pars)
            PBcov[NP+j,i] = PBcov[i,NP+j]

print('Invert cov to get fisher')
PBfisher = np.linalg.inv(PBcov)
'''
dPB = np.zeros((np.size(parc),NPB))

for i, pi in enumerate(Pk_sequence):
    dPB[:,i] = dP(var,parc)

for i, bi in enumerate(Bk_sequence):
    dPB[:,NP+i] = dB(var,parc,pars)
'''
