#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:33:56 2019

@author: cg411
"""

# GW data analysis function 
from matplotlib import pyplot as plt
import pttools.gw_code.ssmtools as ssm
import pttools.bubble.bubble as b
import numpy as np

h0 = 67.4*1e3 *(1/(1e6 * 3.086e16)) #( 2.1840570317563187e-18 planck  2018)  


#LISA Sensitivity curve 

def LISA_noise_curve_om(f):
    # https://www.cosmos.esa.int/web/lisa/lisa-documents : science requirements section 3.2 
    #update coming soon (noted:18/9/19)
    # note: compare with what i was already using 
    
    
    f1 = 0.4 *1e-3   #0.4mHZ
    f2 = 25 * 1e-3   #25mHZ
    R = 1 + (f/f2)**2 
    
    Si = 5.76 * 1e-48 * (1 + (f1/f)**2) #( 1/s**4Hz????)
    Sii = 3.6 * 1e-41 #( 1/Hz)
    S = 1/2 * 20/3 * ( Si/(2*np.pi*f)**4 + Sii)* R
    
    return S_to_om(S,f)


def get_LISA_setup(freqs, n_yrs):
        
    
    yr=365.25*86400.
    t_obs = 3*yr
    Om_n_ana  = LISA_noise_curve_om(freqs)
    
    return Om_n_ana ,t_obs



def calc_pl(freqs, om_noise_curve,SNR, Tobs ): 
    
    """
    Calculates the power law sensitvity for  a given SNR and Tob. 
    The LISA noise curve  is eneterd as om_noise_curve
    Check factor route 2 
    """
    
    beta = np.arange(-7,7,0.2)
    f_ref = 1e-2
    
    om_beta = np.zeros(len(beta))

    PLs = np.zeros([len(beta),len(freqs)])
    
    
    for i, B in enumerate(beta):

        y = (freqs/f_ref)**(2*B) /om_noise_curve**2 
        integration = np.trapz(y,freqs)
        om_beta[i]= SNR * 1/(np.sqrt(Tobs*integration))
        PLs[i,:] = om_beta[i] *  (freqs/f_ref)**(B)
       
    PL_curve = np.amax(PLs, axis=0)
    return PL_curve 

def CalcSNR(freqs, GWSignalSpectrum, NoiseSpectrum,Tobs ):
    
    
    
    y = GWSignalSpectrum**2 /NoiseSpectrum**2 
    integration = np.trapz(y,freqs)
    return np.sqrt(Tobs*integration)


# Models for the GW power spectrum for PTs

def gw_spec_ssm_zb(z,zb,Om_peak,rb):
        """
        GW power spectrum from the SSM (double broken power law)
        p and b correspond to peak and break respectively.
        s =~ k/kp ~f/fp
        
        arXiv:1909.10040
        """
        bp = 1/rb
        zp = bp * zb
        s = z/zp
        m = (9 + bp**4)/(1 + bp**4) 
        M = s**9*((1 + bp**4)/(1 + (bp*s)**4))**2 * (5/(5 - m + m*s**2 ))**(5/2)
        GW = Om_peak*M
        
        return GW 
        

def gw_spec_ssm(s,Om_peak,rp):
        """
        GW power spectrum from the SSM (double broken power law)
        p and b correspond to peak and break respectively.
        s =~ k/kp ~f/fp
        
        arXiv:1909.10040
        """
        bp = 1/rp
        m = (9 + bp**4)/(1 + bp**4) 
        M = s**9*((1 + bp**4)/(1 + (bp*s)**4))**2 * (5/(5 - m + m*s**2 ))**(5/2)
        GW = Om_peak*M
        
        return GW
    
    
def M(s,rp):
    
        bp = 1/rp
        m = (9 + bp**4)/(1 + bp**4) 
        M = s**9*((1 + bp**4)/(1 + (bp*s)**4))**2 * (5/(5 - m + m*s**2 ))**(5/2)
        return M 
    
    
    

    
    
    
def gw_spec_cwg(s,om_peak):
        """
        GW power spectrum (broken power law), Used for PTs in the LISA Cwg
        """
        
        return om_peak * s**3 * (7./(4. + 3.*s**2))**3.5
def M_cwg(s):
    return s**3 * (7./(4. + 3.*s**2))**3.5
    



# Conversion functions 
def S_to_om(S,freq):
    
    """
    Converts strain spectral sensitivity into energy density 
    """
     
    hn_ana_squared = freq * S
    
    return 2*np.pi **2 /(3 * h0**2) * freq**2 *hn_ana_squared




def om_to_S(om, freq):
    """
    converts energy density to noise
    """
    
    return 3 * h0**2 * om /(2 * np.pi**3 * freq **3 )
    
    
def S_to_h(S,freq):
    """
    Converts noise function into charateristic strain 
    """
    return np.sqrt(freq * S)



def h_to_om(freq,h):
    """ converts characteristic strai into energy density 
    """
    return 2*np.pi **2 /(3 * h0**2) * freq**2 *h**2
    

# Foreground models 
#    
    
def om_galactic_binary(freqs):
    
    """
    Analytic prediction for the SGWB / foreground from galactic binaries ( white dwarfs) 
    
    # parameters taken from Cornish and Robson 2017 for 4yr integration
    """
    alpha = 0.138
    beta = -221
    kappa = 521
    gamma =1680
    fk = 0.00113
    
    X = -freqs**alpha +beta*freqs*np.sin(kappa*freqs)
    Y = 1 + np.tanh(gamma*(fk - freqs))
    A = 1.80e-44
    S = A*freqs**(-7/3)*np.exp(X)*Y
    S = np.nan_to_num(S)
    
    return  S_to_om(S,freqs)



def om_compact_binary_ligo(freqs):
    """
    Analytic estimation for foreground for unresolved BBHs and BNSs
    arXiv:1809.10360v2
    LIGO band continuation
    """
    A = 1.8e-9/(25**(2/3)) # + 2.7 - 1.3
    return A*freqs**(2/3)
    



def om_to_om_tilde(om_p,vw,alpha, rb):
        
        k= b.get_ke_frac(vw,alpha) # efficientcy factor 

        mu = 4.78 - 6.27 * rb + 3.34 * rb**2
        
        P =3*(k)**2/mu   

        return om_p/P
    
def alpha_vw_to_rb(alpha, vw):
    """
    Relationship between rb and vw and alpha, obtain by fitting double broken power law to the scaled power spectra frpom ssmtools
    """
    
    a0 =  -32.892* alpha + 6.18 #11.4
    a1 = 173.9*alpha -26.7
    a2 = -287.75*alpha +	37.4
    a3 = 150.3* alpha -16.58
    
    return a0 + a1 * vw + a2 * vw**2 + a3 * vw**3

def get_plot_with_sensitivity(power_law=False):
    """
    power_law=False
    Returns axis for a plot, already containing 5yr sensitivity,
    and optionally (``power_law=True``) power law sensitivity.
    """
    
    freqs = np.logspace(-7,1,num=100,base =10)
    
    yr=365.25*86400.
    Tobs = 3*yr


    PLS = calc_pl(freqs, LISA_noise_curve_om(freqs),10, Tobs )
        
        
    plt.figure()
    plt.loglog(freqs, LISA_noise_curve_om(freqs),linewidth = 2.2,color = 'k', label = 'LISA SR sensitivity curve')
#    #plt.loglog(freqs, OmEff, ls=':', label=r"$\Omega_\mathrm{eff}$ sensitivity")
#    plt.loglog(freqs, np.sqrt(bin_vars), label=r"$\Omega_\mathrm{eff}$ 5yr sensitivity")
    if power_law:
        plt.loglog(freqs, PLS, ls=':', color='k',linewidth = 2.2,
                   label="PL sensitivity curve, SNR=10, Tobs = 3 years ")
    plt.xlabel(r"$f$ [Hz]")
    plt.ylabel(r"$\Omega(f)$")
    plt.axis([1e-5, 1, 1e-14, 1e-6 ])
    plt.grid(True, which='both',alpha=0.5)
#    plt.legend()
    plt.tight_layout()
    return plt.gca()




#
#
#def gw_spec_ssm_z2(z,zp,Om_peak,zb):
#        """
#        GW power spectrum from the SSM (double broken power law)
#        p and b correspond to peak and break respectively.
#        s =~ k/kp ~f/fp
#        
#        arXiv:1909.10040
#        """
#        s = z/zp
#        bp = zp/zb
#        m = (9 + bp**4)/(1 + bp**4) 
#        M = s**9*((1 + bp**4)/(1 + (bp*s)**4))**2 * (5/(5 - m + m*s**2 ))**(5/2)
#        GW = Om_peak*M
#        
#        return GW
#    
#def gw_spec_ssm_z3(z,zb,Om_peak,bp):
#        """
#        GW power spectrum from the SSM (double broken power law)
#        p and b correspond to peak and break respectively.
#        s =~ k/kp ~f/fp
#        
#        arXiv:1909.10040
#        """
#        zp = bp * zb
#        s = z/zp
#        m = (9 + bp**4)/(1 + bp**4) 
#        M = s**9*((1 + bp**4)/(1 + (bp*s)**4))**2 * (5/(5 - m + m*s**2 ))**(5/2)
#        GW = Om_peak*M
#        return GW    

#
#def gw_spec_ssm_rescale(z,zb,Om_peak,rb):
#        """
#        GW power spectrum from the SSM (double broken power law)
#        p and b correspond to peak and break respectively.
#        s =~ k/kp ~f/fp
#        
#        arXiv:1909.10040
#        """
#        bp = 1/rb
#        zp = bp * zb
#        s = z/zp
#        m = (9 + bp**4)/(1 + bp**4) 
#        M = (bp * s)**9 * 4/(1 + (bp*s)**4)**2 * ((5 - m + m/bp**2)/(5 - m + m*s**2 ))**(5/2)
#        GW = Om_peak*M
#        
#        return GW 
#
#def gw_spec_ssm_const_zb(z,Om_peak,rb):
#        """
#        GW power spectrum from the SSM (double broken power law)
#        p and b correspond to peak and break respectively.
#        s =~ k/kp ~f/fp
#        
#        arXiv:1909.10040
#        """
#        zp = 1/rb
#        bp = 1/rb
#        s = z/zp
#        m = (9 + bp**4)/(1 + bp**4) 
#        M = s**9*((1 + bp**4)/(1 + (bp*s)**4))**2 * (5/(5 - m + m*s**2 ))**(5/2)
#        GW = Om_peak*M
#        
#        return GW    
#
        

