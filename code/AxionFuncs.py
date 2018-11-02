import numpy as np
from numpy import pi, sqrt, exp, zeros, size, shape, sinc, linspace, trapz, loadtxt, interp
from scipy.integrate import cumtrapz, quad

#================================AxionFuncs.py=================================#
# Contains:
#==============================================================================#


#==============================================================================#
def AxionFlux_Primakoff(gag,E):
    norm = 6.02e10*(gag/1e-10)**2.0
    return norm*((E**2.481)/exp(E/1.205))

def AxionFlux_Compton(gae,E):
    norm = 13.314e6*(gae/1e-13)**2.0
    return norm*((E**2.987)/exp(E*0.776))

def AxionFlux_Brem(gae,E):
    norm = 26.311e8*(gae/1e-13)**2.0
    return norm*E*exp(-0.77*E)/(1+0.667*E**1.278)

def AxionFlux_AxioRecomb(gae,E):
    #Solar differential Axion Flux from the axion electron coupling
    #column 1 = Energy [keV]
    #column 2 = Axion Flux 1/[10^19 keV cm^2 day]
    #I calculated the flux using gae = 0.511*10^-10 
    #for other values of gae use:
    #FLUX = Table*[gae/(0.511*10^-10)]^2
    data = loadtxt('gaeflux.txt')
    E1 = data[:,0]
    F1 = data[:,1]
    norm = 1e19*(gae/(0.511e-10))**2.0/(3600*24)
    Flux = interp(E,E1,F1)*norm
    return Flux
#==============================================================================#





#==============================================================================#
def PhotonNumber_Primakoff(E,m_a):
    norm = 220893
    return norm*((E**2.481)/exp(E/1.205))*(sinc(25380.710659898483/pi*m_a**2.0/E))**2.0

def PhotonNumber_Electron(Flux,E,m_a):
    norm = 220893/(6.02e10)
    return norm*Flux*(sinc(25380.710659898483/pi*m_a**2.0/E))**2.0

def smear(dN,E,E_res):
    n = size(dN)
    Norm = 1.0/sqrt(2*pi*E_res**2.0)
    dN_smeared = zeros(shape=n)
    for i in range(0,n):
        K = Norm*exp(-(E-E[i])**2.0/(2*E_res**2.0))
        dN_smeared[i] = trapz(K*dN,E)
    return dN_smeared

def smearFast(dN,E,E_res):
    n = size(dN)
    dE = E[1]-E[0]
    irange = int(3*E_res/dE)
    Norm = 1.0/sqrt(2*pi*E_res**2.0)
    dN_smeared = zeros(shape=n)
    for i in range(0,n):
        i1 = max(0,i-irange)
        i2 = min(n-1,i+irange)
        Eint = E[i1:i2]
        K = Norm*exp(-(Eint-E[i])**2.0/(2*E_res**2.0))
        dN_smeared[i] = trapz(K*dN[i1:i2],Eint)
    return dN_smeared


def BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=100,res_on=False): 
    E_bin_edges = linspace(E_min,E_max,nE_bins+1)
    E_bw = (E_max-E_min)/(nE_bins+1.0)
    E_bins = (E_bin_edges[1:]+E_bin_edges[:-1])/2

    # TRAPZ
    nm = size(m_vals)
    R1_tab = zeros(shape=(nE_bins,nm))
    #Ei = linspace(E_min,E_max,nfine*nE_bins)
    Ei = zeros(shape=(nE_bins*nfine))
    for i in range(0,nE_bins):
        Ei[i*nfine:(i+1)*nfine] = linspace(E_bin_edges[i],E_bin_edges[i+1]-E_bw/nfine,nfine)
        
    for j in range(0,nm):
        dN = PhotonNumber_Primakoff(Ei,m_vals[j])
        if res_on:
            dN = smear(dN,Ei,E_min)
        for i in range(0,nE_bins):
            Ebin = Ei[i*nfine:(i+1)*nfine]
            dNbin = dN[i*nfine:(i+1)*nfine]
            R1_tab[i,j] = sum(0.5*(Ebin[1:]-Ebin[0:-1])*(dNbin[1:]+dNbin[0:-1]))
    
    # Get m = 0 rate
    R0 = zeros(shape=(nE_bins))
    dN = PhotonNumber_Primakoff(Ei,0.0)
    if res_on:
        dN = smear(dN,Ei,E_min)
    for i in range(0,nE_bins):
        Ebin = Ei[i*nfine:(i+1)*nfine]
        dNbin = dN[i*nfine:(i+1)*nfine]
        R0[i] = sum(0.5*(Ebin[1:]-Ebin[0:-1])*(dNbin[1:]+dNbin[0:-1]))

    return E_bins,R1_tab,R0

def BinnedPhotonNumberTable_Electron(m_vals,E_min,E_max,nE_bins,nfine=100,res_on=False): 
    E_bin_edges = linspace(E_min,E_max,nE_bins+1)
    E_bw = (E_max-E_min)/(nE_bins+1.0)
    E_bins = (E_bin_edges[1:]+E_bin_edges[:-1])/2

    # TRAPZ
    nm = size(m_vals)
    R1_tab = zeros(shape=(nE_bins,nm))
    #Ei = linspace(E_min,E_max,nfine*nE_bins)
    Ei = zeros(shape=(nE_bins*nfine))
    for i in range(0,nE_bins):
        Ei[i*nfine:(i+1)*nfine] = linspace(E_bin_edges[i],E_bin_edges[i+1]-E_bw/nfine,nfine)
       
    Flux = AxionFlux_AxioRecomb(1e-10,Ei)
    for j in range(0,nm):
        dN = PhotonNumber_Electron(Flux,Ei,m_vals[j])
        if res_on:
            dN = smearFast(dN,Ei,E_min)
        for i in range(0,nE_bins):
            Ebin = Ei[i*nfine:(i+1)*nfine]
            dNbin = dN[i*nfine:(i+1)*nfine]
            R1_tab[i,j] = sum(0.5*(Ebin[1:]-Ebin[0:-1])*(dNbin[1:]+dNbin[0:-1]))
        
    # Get m = 0 rate
    R0 = zeros(shape=(nE_bins))
    dN =  PhotonNumber_Electron(Flux,Ei,0.0)
    if res_on:
        dN = smear(dN,Ei,E_min)
    for i in range(0,nE_bins):
        Ebin = Ei[i*nfine:(i+1)*nfine]
        dNbin = dN[i*nfine:(i+1)*nfine]
        R0[i] = sum(0.5*(Ebin[1:]-Ebin[0:-1])*(dNbin[1:]+dNbin[0:-1]))

    return E_bins,R1_tab,R0
