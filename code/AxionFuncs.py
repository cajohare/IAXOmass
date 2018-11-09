#================================AxionFuncs.py=================================#
# Written by C. O'Hare
# Contains:
# Functions for calculating Solar Axion fluxes for photon & electron coupling
# Functions for calculating X-ray spectra in a haloscope
# Functions to smear X-ray spectra by an angular resolution
# Script to generate binned X-ray data for given flux+experiment
#==============================================================================#


from numpy import pi, sqrt, exp, zeros, size, shape
from numpy import sinc, linspace, trapz, loadtxt, interp
from scipy.integrate import cumtrapz, quad


#==============================================================================#
def AxionFlux_Primakoff(gag,E):
    # Parameterised differential Axion Flux in [cm^-1 s^-1 keV^-1]
    # gag = Axion-photon coupling in GeV^-1
    # E = Axion/X-ray energy in keV
    norm = 6.02e10*(gag/1e-10)**2.0
    return norm*((E**2.481)/exp(E/1.205))

def AxionFlux_Axioelectron(gae,E):
    # Differential Axion Flux from the axion electron coupling
    # Flux = AxionRecomb+Compton+Bremsstrahlung
    # column 1 = Energy [keV]
    # column 2 = Axion Flux 1/[10^19 keV cm^2 day]
    # Output: flux in cm^-1 s^-1 keV^-1
    # gae = Axion-electron coupling in GeV^-1
    # E = Axion/Xray energy in keV
    data = loadtxt('gaeflux.txt')
    E1 = data[:,0]
    F1 = data[:,1]
    norm = 1e19*(gae/(0.511e-10))**2.0/(3600*24)
    Flux = interp(E,E1,F1)*norm
    return Flux

def AxionFlux_Compton(gae,E):
    # Parameterised Compton axion flux (unused in paper)
    norm = 13.314e6*(gae/1e-13)**2.0
    return norm*((E**2.987)/exp(E*0.776))

def AxionFlux_Brem(gae,E):
    # Parameterised Bremsstrahlung axion flux (unused in paper)
    norm = 26.311e8*(gae/1e-13)**2.0
    return norm*E*exp(-0.77*E)/(1+0.667*E**1.278)

#==============================================================================#



#==============================================================================#
def PhotonNumber_Primakoff(Flux_scale,E,m_a,\
                           Bfield=2.5,Exposure=1.5,Length=20.0,\
                           N_bores=8,BoreDiameter=60.0,eps_D=0.7,eps_T=0.8):
    # differential Xray count dN/dE (in keV^-1) for axionphoton flux
    # (Optional) Flux_scale = scaling for normalisation (set to 1 for units used in paper)
    # E = Xray energy (keV)
    # m_a = axion mass (eV)
    norm,normq = NgammaNorm(Bfield,Exposure,Length,N_bores,BoreDiameter,eps_D,eps_T)
    norm = Flux_scale*norm
    return norm*((E**2.481)/exp(E/1.205))*(sinc(normq/pi*m_a**2.0/E))**2.0  # keV^-1

def PhotonNumber_Electron(Flux,E,m_a,\
                           Bfield=2.5,Exposure=1.5,Length=20.0,\
                           N_bores=8,BoreDiameter=60.0,eps_D=0.7,eps_T=0.8):
    # differential Xray count dN/dE (in keV^-1) for axionelectron flux
    # Flux_scale = scaling for normalisation (set to 1 for units used in paper)
    # E = Xray energy (keV)
    # m_a = axion mass (eV)
    norm,normq = NgammaNorm(Bfield,Exposure,Length,N_bores,BoreDiameter,eps_D,eps_T)
    norm = norm/(6.02e10)
    return norm*Flux*(sinc(normq/pi*m_a**2.0/E))**2.0 # keV^-1

def NgammaNorm(Bfield,Exposure,Length,N_bores,BoreDiameter,eps_D,eps_T):
    # Nnorm = normalisation of overall photon number to get it in keV^-1 and constant that enters into t
    S_cm = N_bores*pi*(BoreDiameter/2.0)**2.0 # cm^2
    L_eV = Length/1.97e-7 # eV^-1
    t_secs = Exposure*3600*24*365 # s
    B = Bfield*(1e-19*195)
    norm = 6.02e10*t_secs*S_cm*eps_D*eps_T*(B*L_eV/2.0)**2.0
    normq = L_eV/(4*1000)
    return norm,normq
#==============================================================================#


#==============================================================================#
def smear(dN,E,E_res):
    # Smear spectrum dN(E) by energy resolution Eres
    # dN = spectrum (arbitrary units)
    # E = Energies defining dN
    # E_res = Energy resolution to smear by
    n = size(dN)
    Norm = 1.0/sqrt(2*pi*E_res**2.0)
    dN_smeared = zeros(shape=n)
    for i in range(0,n):
        # Each new energy is the full spectrum convolved by a gaussian
        K = Norm*exp(-(E-E[i])**2.0/(2*E_res**2.0))
        dN_smeared[i] = trapz(K*dN,E)
    return dN_smeared

def smearFast(dN,E,E_res):
    # Does the same as 'smear' but is faster and less accurate for E_res>100 eV
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
#==============================================================================#




#==============================================================================#
def EnergyBins(E_min,E_max,nfine,nE_bins):
    # Define energy array for doing the trapz integration below
    # E_min = energy threshold
    # E_max = max energy
    # nfine = number of energies within one bin to integrate over
    # nE_bins = number of energy bins between E_min and E_max
    E_bin_edges = linspace(E_min,E_max,nE_bins+1)
    E_bw = (E_max-E_min)/(nE_bins+1.0)
    E_bins = (E_bin_edges[1:]+E_bin_edges[:-1])/2
    
    Ei = zeros(shape=(nE_bins*nfine))
    for i in range(0,nE_bins):
        Ei[i*nfine:(i+1)*nfine] = linspace(E_bin_edges[i],E_bin_edges[i+1]-E_bw/nfine,nfine)
        
    return Ei,E_bins

def BinnedPhotonNumberTable(m_vals,E_min,E_max,nE_bins,coupling='Photon',\
                            nfine=100,res_on=False,\
                           Bfield=2.5,Exposure=1.5,Length=20.0,\
                           N_bores=8,BoreDiameter=60.0,eps_D=0.7,eps_T=0.8): 
    # Generate tabulated values of data for a range of axion masses 
    # OUTPUT: R1_tab = Tabulated values of the binned Xray counts (columns) vs axion mass (rows)
    # R0 = massless data
    # E_bins = centers of energy bins
    # INPUT: m_vals = masses to add to the tabulation 
    # E_min = threshold energy (also resolution if res_on=True)
    # E_max = maximum energy
    # nE_bins = number of energy bins
    # coupling = 'Photon' or 'Electron' for g_ag or g_ae
    # nfine = number of points to integrate over within one bin (controls accuracy)
    # res_on = True/False, whether to do energy resolution integral or not
    nm = size(m_vals)
    R1_tab = zeros(shape=(nE_bins,nm))
    Ei,E_bins = EnergyBins(E_min,E_max,nfine,nE_bins)
    
    if coupling=='Electron':
        Flux = AxionFlux_Axioelectron(1e-10,Ei)
        dN_func = PhotonNumber_Electron
    else:
        Flux = 1.0
        dN_func = PhotonNumber_Primakoff
        
    # Tabulate m != 0 rates    
    for j in range(0,nm):
        dN = dN_func(Flux,Ei,m_vals[j],\
                     Bfield,Exposure,Length,\
                     N_bores,BoreDiameter,eps_D,eps_T)  
        if res_on:
            dN = smear(dN,Ei,E_min)
        for i in range(0,nE_bins):
            Ebin = Ei[i*nfine:(i+1)*nfine]
            dNbin = dN[i*nfine:(i+1)*nfine]
            R1_tab[i,j] = sum(0.5*(Ebin[1:]-Ebin[0:-1])*(dNbin[1:]+dNbin[0:-1]))
        
    # Get m = 0 rate
    R0 = zeros(shape=(nE_bins))
    dN =  dN_func(Flux,Ei,0.0)
    if res_on:
        dN = smear(dN,Ei,E_min)
    for i in range(0,nE_bins):
        Ebin = Ei[i*nfine:(i+1)*nfine]
        dNbin = dN[i*nfine:(i+1)*nfine]
        R0[i] = sum(0.5*(Ebin[1:]-Ebin[0:-1])*(dNbin[1:]+dNbin[0:-1]))

    return E_bins,R1_tab,R0
#==============================================================================#
