from numpy import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz, quad
from iminuit import minimize
from scipy.special import gammaln
import AxionFuncs
import Like

# Generate Tabulated spectra
nm = 10000
m_vals = logspace(-4.0,2e0,nm)
E_max = 11.0
nE_bins = 500
n_per_bin = 10

n_DL = 10000
m_DL_vals = logspace(log10(1e-3),log10(2e-1),n_DL)


E_min = 1.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res001eV.txt',data)
print E_min


E_min = 10.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res010eV.txt',data)
print E_min


E_min = 50.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res050eV.txt',data)
print 1


E_min = 100.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res100eV.txt',data)
print 1



E_min = 200.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res200eV.txt',data)
print 1



E_min = 300.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res300eV.txt',data)
print 1



E_min = 400.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res400eV.txt',data)
print 1



E_min = 500.0e-3
E_bins,R1_tab,R0 = AxionFuncs.BinnedPhotonNumberTable_Primakoff(m_vals,E_min,E_max,nE_bins,nfine=n_per_bin,res_on=True)
data = zeros(shape=(nE_bins+1,nm+1))
data[1:,0] = E_bins
data[0,1:] = m_vals
data[1:,1:] = R1_tab
savetxt('XrayTab_gag_500bins_res500eV.txt',data)
print 1
