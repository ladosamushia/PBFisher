import scipy.interpolate as itp
import numpy as np

from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef("double Pk(double, double *, double *);")
ffibuilder.cdef("double Bk(double, double, double, double *, double *, double *);")

ffibuilder.set_source("_BFisher",
r""" 
#include <math.h>

 double Pk(double Pk1, double *var, double *parc) {
   double k = var[0];
   double mu = var[1];
   double apar = parc[0];
   double aper = parc[1];
   double f = parc[2];
   double b1 = parc[3];

   k = k*sqrt(aper*aper+mu*mu*(apar*apar-aper*aper));
   mu = mu*apar/sqrt(aper*aper+mu*mu*(apar*apar-aper*aper));

   return (b1 + f*mu)*(b1 + f*mu)*Pk1;
 }

 double Bk(double Pk1, double Pk2, double Pk3, double *var, double *parc,
           double *pars) {
   double k1  = var[0];
   double k2  = var[1];
   double k3  = var[2];
   double mu1  = var[3];
   double phi12 = var[4];

   double apar = parc[0];
   double aper = parc[1];
   double f = parc[2];
   double b1 = parc[3];
   double b2 = parc[4];

   double nave = pars[0];

   double mu12 = (k3*k3 - k1*k1 - k2*k2)/2/k1/k2;
   double mu2 = mu1*mu12 - sqrt(1 - mu1*mu1)*sqrt(1 - mu12*mu12)*cos(phi12);
   double mu3 = -(mu1*k1 + mu2*k2)/k3;

   k1 = k1*sqrt(aper*aper+mu1*mu1*(apar*apar-aper*aper));
   k2 = k2*sqrt(aper*aper+mu2*mu2*(apar*apar-aper*aper));
   k3 = k3*sqrt(aper*aper+mu3*mu3*(apar*apar-aper*aper));
   mu1 = mu1*apar/sqrt(aper*aper+mu1*mu1*(apar*apar-aper*aper));
   mu2 = mu2*apar/sqrt(aper*aper+mu2*mu2*(apar*apar-aper*aper));
   mu3 = mu3*apar/sqrt(aper*aper+mu3*mu3*(apar*apar-aper*aper));

   mu12 = (k3*k3 - k1*k1 - k2*k2)/2/k1/k2;
   double mu31 = -(k1 + k2*mu12)/k3;
   double mu23 = -(k1*mu12 + k2)/k3;

   double k12 = sqrt(k1*k1 + k2*k2 + 2*k1*k2*mu12);
   double k23 = sqrt(k2*k2 + k3*k3 + 2*k2*k3*mu23);
   double k31 = sqrt(k3*k3 + k1*k1 + 2*k3*k1*mu31);

   double Z1k1 = b1 + f*mu1*mu1;
   double Z1k2 = b1 + f*mu2*mu2;
   double Z1k3 = b1 + f*mu3*mu3;

   double F12 = 5./7. + mu12/2*(k1/k2 + k2/k1) + 2./7.*mu12*mu12;
   double F23 = 5./7. + mu23/2*(k2/k3 + k3/k2) + 2./7.*mu23*mu23;
   double F31 = 5./7. + mu31/2*(k3/k1 + k1/k3) + 2./7.*mu31*mu31;

   double G12 = 3./7. + mu12/2*(k1/k2 + k2/k1) + 4./7.*mu12*mu12;
   double G23 = 3./7. + mu23/2*(k2/k3 + k3/k2) + 4./7.*mu23*mu23;
   double G31 = 3./7. + mu31/2*(k3/k1 + k1/k3) + 4./7.*mu31*mu31;

   double mu1p2 = (mu1*k1 + mu2*k2)/k12;
   double mu2p3 = (mu2*k2 + mu3*k3)/k23;
   double mu3p1 = (mu3*k3 + mu1*k1)/k31;

   double Z2k12 = b2/2. + b1*F12 + f*mu1p2*mu1p2*G12;
   Z2k12 += f*mu1p2*k12/2.*(mu1/k1*Z1k2 + mu2/k2*Z1k1);
   double Z2k23 = b2/2. + b1*F23 + f*mu2p3*mu2p3*G23;
   Z2k23 += f*mu2p3*k23/2.*(mu2/k2*Z1k3 + mu3/k3*Z1k2);
   double Z2k31 = b2/2. + b1*F31 + f*mu3p1*mu3p1*G31;
   Z2k31 += f*mu3p1*k31/2.*(mu3/k3*Z1k1 + mu1/k1*Z1k3);

   double Bi = 2*Z2k12*Z1k1*Z1k2*Pk1*Pk2;
   Bi += 2*Z2k23*Z1k2*Z1k3*Pk2*Pk3;
   Bi += 2*Z2k31*Z1k3*Z1k1*Pk3*Pk1;
   
   double var1[2] = {k1,mu1};
   double var2[2] = {k2,mu2};
   double var3[2] = {k3,mu3};

   Bi += (Pk(Pk1,var1,parc) + Pk(Pk2,var2,parc) + Pk(Pk3,var3,parc))/nave;
   
   return Bi;
 }

 double CovP(double Pk1, double* var,double* parc, double* pars){
    double navg = pars[0];
    double C = (Pk(Pk1,var,parc) + 1/navg);
    C *= C;
    return C;
 }

 double CovB(double Pk1, double Pk2, double Pk3, double* var, double* parc, double* pars){
 
    double k1 = var[0];
    double k2 = var[1];
    double k3 = var[2];
    double mu1 = var[3];
    double phi12 = var[4];

    double navg = pars[0];
    double Vs = pars[1];

    double mu12 = (k3*k3 - k1*k1 - k2*k2)/2/k1/k2;
    double mu2 = mu1*mu12 - sqrt(1 - mu1*mu1)*sqrt(1 - mu12*mu12)*cos(phi12);
    double mu3 = -(mu1*k1 + mu2*k2)/k3;
    
    double pvar1[2] = {k1,mu1};
    double pvar2[2] = {k2,mu2};
    double pvar3[2] = {k3,mu3};

    double C = Pk(Pk1,pvar1,parc) + 1/navg;
    C *= Pk(Pk2,pvar2,parc) + 1/navg;
    C *= Pk(Pk3,pvar3,parc) + 1/navg;
    C *= Vs;

    return C;
 }

 double CrossPB(double Pk1, double Pk2, double Pk3, double* var1,double* var2,
                double* parc,double* pars){

    double k1 = var2[0];
    double k2 = var2[1];
    double k3 = var2[2];
    double mu1 = var2[3];
    double phi12 = var2[4];

    double navg = pars[0];

    double mu12 = (k3*k3 - k1*k1 - k2*k2)/2/k1/k2;
    double mu2 = mu1*mu12 - sqrt(1 - mu1*mu1)*sqrt(1 - mu12*mu12)*cos(phi12);
    double mu3 = -(mu1*k1 + mu2*k2)/k3;
    
    double pvar1[2] = {k1,mu1};
    double pvar2[2] = {k2,mu2};
    double pvar3[2] = {k3,mu3};

    double C = Bk(Pk1,Pk2,Pk3,var2,parc,pars);
    C += Pk(Pk1,pvar1,parc)/navg;
    C += Pk(Pk2,pvar2,parc)/navg; 
    C += Pk(Pk3,pvar3,parc)/navg;  
    C += 1/navg/navg;

    return C;
 }

""")

if __name__ == "__main__":
   ffibuilder.compile(verbose=True) 

"""

def dP(var,parc):
    '''
    Derivative of Power spectrum with respect to cosmological parameters - parc.
    var are the 5D variables k1,k2,k3,mu1,phi12 k1, k2, k3, mu1, phi12 = var
    '''
    eps = 1e-6
    derP = np.zeros(np.shape(parc))
    P0 = Pk(var,parc)
    for i in range(np.shape(parc)[0]):
        parc1 = np.copy(parc)
        parc1[i,:] += eps
        P1 = Pk(var,parc1)
        derP[i,:] = (P1 - P0)/eps
    return derP
    
def dB(var,parc,pars):
    '''
    Derivative of Bispectrum with respect to cosmological parameters - parc. var are
    the 5D variables k1,k2,k3,mu1,phi12
    k1, k2, k3, mu1, phi12 = var
    '''
    eps = 1e-6
    derB = np.zeros(np.shape(parc))
    B0 = Bisp(var,parc,pars)
    for i in range(np.shape(parc)[0]):
        parc1 = np.copy(parc)
        parc1[i,:] += eps
        B1 = Bisp(var,parc1,pars)
        derB[i,:] = (B1 - B0)/eps
    return derB


def Fishz(pars,parc):
    '''
    Compute Bispectrum Fisher Matrix at a fixed redshift for survey parameters -
    pars (that are known/fixed), and cosmological parameters - parc (that are
    unknown).
    pars = [nbar, Volume]
    parc = [Omega_m, sigma_8, n_s, w_0, w_a, Omega_b, h0]
    '''
    # Bispectrum in redshift space is a 5D function. Integrate over 5D with MC.
    RR = random.rand(NMC,5)
    for i in range(NMC):
        eta1 = etamax*RR[i,0]
        eta2 = etamax*RR[i,1]
        eta3 = etamax*RR[i,2]
        k1 = np.sqrt(2*eta1)
        k2 = np.sqrt(2*eta2)
        k3 = np.sqrt(2*eta3)
        mu1 = 2*mumax*RR[i,3] - 1
        phi12 = 2*np.pi*RR[i,4]
        var = [k1, k2, k3, mu1, phi12]
        CB = CovB()
        derB = dB()
        FM += np.outer(dB,dB)/CB

    return FM
"""
