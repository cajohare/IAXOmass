# IAXOmass
Code for reproducing results from arXiv:[1811.?????] on how IAXO can use the energy dependence of the X-ray spectrum from Solar axion conversion to measure the axion mass over a small range when the discrimination is not destroyed by oscillations.

The code is written in python, and all the plots are generated using jupyter notebooks. This repository contains all the code needed to reproduce the paper.

---


## Fig. 1
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/AxionLimits.png" width="420" height="400">

[Click here for the notebook](https://github.com/cajohare/IAXOmass/blob/master/code/plot_FinalLimit.ipynb)

The axion parameter space g-m_a. The main result of this paper is shown as a blue region: the median range of masses and couplings for which IAXO can determine the axion mass to be non-zero with 3sigma significance. The QCD axion band is shaded in orange and the common benchmark KSVZ and DFSZ models are drawn as straight brown lines. In various shades of green are axions excluded by astrophysical arguments. The best experimental results shown in shades of red. All the limit data can be found in the 'limit_data' folder

---

## Fig. 2
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/SolarAxionFlux.png" width="420" height="400">

[Click here for the notebook](https://github.com/cajohare/IAXOmass/blob/master/code/plot_MassDiscovery_Electron.ipynb)

The Solar axion fluxes expected on Earth and their components due to g_ae and g_ag. The g_ag flux is analytic and can be found in AxionFuncs.py. The g_ae flux is tabulated in gaeflux.txt also found in the code folder.

---

## Fig. 3
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/XraySpectra_Photon.png" width="840" height="400">
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/XraySpectra_Electron.png" width="840" height="400">

[Click here for the notebook (upper panel)](https://github.com/cajohare/IAXOmass/blob/master/code/plot_XraySpectra_Photon.ipynb)

[Click here for the notebook (lower panel)](https://github.com/cajohare/IAXOmass/blob/master/code/plot_XraySpectra_Electron.ipynb)

X-ray spectrum from Solar axion conversion in 2.5 T magnet with a projected length of L=20 m. Displaying spectra for different values of the axion mass as well as for both the Solar axioelectron (upper panels) and Primakoff (lower panels) fluxes. The left hand panel in both cases show the underlying idealised spectra, whereas the right hand panels show the spectra after being convolved with a Gaussian of width E_0 = 100 eV. For comparison, we have normalised all spectra to one.

---

## Fig. 4
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/Like0.png" width="420" height="400">
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/Like1.png" width="420" height="400">

<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/Like2.png" width="420" height="400">
<img src="https://github.com/cajohare/IAXOmass/blob/master/plots/plots_png/Like3.png" width="420" height="400">

[Click here for the notebook](https://github.com/cajohare/IAXOmass/blob/master/code/plot_Like.ipynb)
