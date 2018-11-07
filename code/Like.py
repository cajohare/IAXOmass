#================================Like.py=======================================#
# Written by C. O'Hare
# Contains:
# InterpExpectedEvents: convert tabulated data to data for arbitrary mass
# lnPF: logarithm of poisson pdf
# llhood2: -1*log-likelihood for params = (m_a,g)
# llhood2_marg: profile likelihood for params = (m_a)
# llhood2_marg0: value of llhood2_marg at m_a=0
# llhood1: massless model likelihood params = (g)
# ConstantObsNumberLine: constant event number lines
# MassDiscoveryLimit_Simple: analytic median mass discovery limit
# MassDiscoveryLimit_Minuit: numerical median mass discovery limit
# MassEstimationDiscoveryLimit: numerical median mass estimation limit
#==============================================================================#


from numpy import pi, sqrt, exp, zeros, size, shape, sinc, linspace, logspace
from numpy import log10, floor, log, interp, append, flipud
from scipy.integrate import cumtrapz, quad
from iminuit import minimize
from scipy.special import gammaln

#==============================================================================#
# Important function that interpolates N_exp events at an arbitrary mass (m)
# using the tabulated data stored in R1_tab, then rescales by g
def InterpExpectedEvents(g,m,m_vals,R1_tab):
    nm = size(m_vals)
    m1 = m_vals[0]
    m2 = m_vals[-1]
    i1 = int(floor((log10(m)-log10(m_vals[0]))*(nm-1.0)/(log10(m_vals[-1])-log10(m_vals[0]))+1.0))-1
    i2 = i1+1
    N_exp = 1e40*(g**4.0)*((R1_tab[:,i1]*(m_vals[i2]-m)+R1_tab[:,i2]*(m - m_vals[i1]))/(m_vals[i2]- m_vals[i1]))
    return N_exp
#==============================================================================#




#==============================================================================#
# Likelihoods:

# Log of poisson pdf
def lnPF(Nob,Nex):
    return sum(Nob*log(Nex) - Nex)# - gammaln(Nob+1.0)) #factorial removed for speed

# 2D likelihood for (g,m)
def llhood2(X,N_obs,m_vals,R1_tab):
    m = X[1]
    g = 10.0**X[0]
    N_exp = InterpExpectedEvents(g,m,m_vals,R1_tab)
    LL = -1.0*lnPF(N_obs,N_exp)
    return LL
    
# Profile likelihood for (m)
def llhood2_marg(m,N_obs,m_vals,R1_tab):
    N_exp_10 = InterpExpectedEvents(1e-10,m,m_vals,R1_tab)
    g0 = ((sum(N_obs)/sum(N_exp_10))**0.25)*1e-10
    LL = llhood2([log10(g0),m],N_obs,m_vals,R1_tab)
    return LL

# Profile likelihood for (m=0)
def llhood2_marg0(N_obs,R0):
    N_exp = (sum(N_obs)/sum(R0))*R0
    LL = -1.0*lnPF(N_obs,N_exp)
    return LL

# Profile likelihood for (g,m=0)
def llhood1(X,N_obs,R0):
    g = 10.0**X[0]
    N_exp = 1e40*(g**4.0)*R0
    LL = -1.0*lnPF(N_obs,N_exp)
    return LL
#==============================================================================#




#==============================================================================#
def ConstantObsNumberLine(N_ob,mi,m_vals,R1_tab):
    ni = size(mi)
    g = zeros(shape=ni)
    for i in range(0,ni):
        N_exp_10 = InterpExpectedEvents(1e-10,mi[i],m_vals,R1_tab)
        g[i] = 1e-10*(N_ob/sum(N_exp_10))**0.25
    return g
#==============================================================================#


#==============================================================================#
# Simple analytic calculation of mass discovery limit using formulae derived in
# the paper.
def MassDiscoveryLimit_Simple(m_vals,R1_tab,R0,m_DL_vals):
    nm = size(m_vals)
    n_DL = size(m_DL_vals)
    DL = zeros(shape=n_DL)
    for im in range(0,n_DL):
        m0 = m_DL_vals[im]
        i0 = int(floor((log10(m0)-log10(m_vals[0]))*(nm-1.0)/(log10(m_vals[-1])-log10(m_vals[0]))+1.0))-1
        N = R1_tab[:,i0]
        N0 = R0
        D = sum(R1_tab[:,i0])/sum(R0)
        DL[im] = 1e-10*(9.0/sum(2*N*log(N/(D*N0))))**0.25
    return DL

# Algorithmic calculation of the discovery limit using Minuit likelihood 
# minimisation (currently unused but will be needed for future studies)
def MassDiscoveryLimit_Minuit(m_vals,R1_tab,R0,m_DL_vals,gmin=1e-12,gmax=1e-7,ng=100):
    nm = size(m_vals)
    n_DL = size(m_DL_vals)
    DL = zeros(shape=(n_DL))
    g_vals = logspace(log10(gmin),log10(gmax),ng)
    for im in range(0,n_DL):
        for j in range(0,ng):
            g = g_vals[j]
            m0 = m_DL_vals[im]
            N_obs = InterpExpectedEvents(g,m0,m_vals,R1_tab)

            #print log10(g),log10(m0),sum(N_obs)
            D12_prev = 0.0
            g_prev = g
            if sum(N_obs)>3:
                # ----- Massive case -------- #
                L2 = -1.0*lnPF(N_obs,N_obs)

                #------ Massless case ------#
                X_in1 = [log10(g)]
                res = minimize(llhood1, X_in1, args=(N_obs,R0))
                L1 = res.fun

                # Test statistic
                D12 = -2.0*(L2-L1) # significance for measuring mass
                if D12>9.0: # Median 3sigma detection -> D = 9
                    DL[im] = 10.0**(interp(9.0,[D12_prev,D12],[log10(g_prev),log10(g)]))
                    break
                g_prev = g # Reset for interpolation
                D12_prev = D12
    return DL
#==============================================================================#


#==============================================================================#
# Function for calculating the mass estimation discovery limit 
# (Fig.6 of the paper)
# Can work in two modes:
#
# Mode 1 when gmin_vals = scalar: The scan over coupling values begins with 
# the value that produced N_obs=2 events
#
# Mode 2 (when gmin_vals = array of size (n_DL)): The scan begins at the values 
# specified
#
# The two modes are needed because the first scan over the parameter space is needed
# to eliminate all of the spurious regions at low and high masses that give rise 
# to likelihoods within the confidence interval band 
# (i.e. likelihood ratio>-4 for example) but are no where near the correct mass 
# (i.e. within err away from m0)
def MassEstimationDiscoveryLimit(err,m_vals,R0,R1_tab,m_DL_vals,sigmas=2,gmin_vals=1.0,gmax=1e-9,ng=500,nL=1000):
    itot = nL
    nm = size(m_vals)
    n_DL = size(m_DL_vals)
    DLmerr = zeros(shape=(n_DL))
    for im in range(0,n_DL):
        m0 = m_DL_vals[im]
        i0 = int(floor((log10(m0)-log10(m_vals[0]))\
                       *(nm-1.0)/(log10(m_vals[-1])-log10(m_vals[0]))+1.0))-1
        iupper = int(itot*i0/nm)
        ilower = itot-iupper
        mlower = flipud(logspace(log10(m_vals[1]),log10(m0*(1-err/2)),ilower))
        mupper = logspace(log10(m0*(1+err/2)),log10(m_vals[-2]),iupper)
        if size(gmin_vals)>1:
            gmin = gmin_vals[im]
        else:
            N = R1_tab[:,i0]
            N0 = R0
            D = sum(R1_tab[:,i0])/sum(R0)
            gmin = 1e-10*(1.0/sum(2*N*log(N/(D*N0))))**0.25
            
        g_vals = logspace(log10(gmin),log10(gmax),ng)
        for j in range(0,ng):
            g0 = g_vals[j]
            N_obs = R1_tab[:,i0]*(g0/1e-10)**4.0
            Lmax = llhood2_marg(m0,N_obs,m_vals,R1_tab)
            Lprof0 = -2*(llhood2_marg0(N_obs,R0) - Lmax) 
            Lprofend = -2*(llhood2_marg(m_vals[-2],N_obs,m_vals,R1_tab) - Lmax)
            if (Lprof0<-sigmas**2.0)&(Lprofend<-sigmas**2.0):
                c = True
                for ii in range(0,ilower):
                    Lprof = -2*(llhood2_marg(mlower[ii],N_obs,m_vals,R1_tab) - Lmax)
                    if Lprof>-sigmas**2.0:
                        c = False
                        break
                for ii in range(0,iupper):
                    Lprof = -2*(llhood2_marg(mupper[ii],N_obs,m_vals,R1_tab) - Lmax)
                    if Lprof>-sigmas**2.0:
                        c = False
                        break    
                if c:
                    DLmerr[im] = g0
                    break
                    
    return DLmerr
#==============================================================================#








#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#


# OLD CODE:

#==============================================================================#
# Given an input mass and coupling (m0,g0), calculates the confidence interval
# around the best fit mass to a given significance level (sigmas)
# also requires err_l and err_u which is roughly the fractional size over which
# to search (needs to be very large ~100% for most values of m0)
# This function currently not quite extensive enough to give accurate results
# Unless one requires err<0.01
def MassMeasurement(m0,g0,err_l,err_u,sigmas,m_vals,R1_tab,nL=100): 
    nm = size(m_vals)
    mlim = []
    N_obs = InterpExpectedEvents(g0,m0,m_vals,R1_tab)
    mi = linspace((1-err_l)*m0,(1+err_u)*m0,nL)
    LL = zeros(shape=nL)
    maxL = llhood2([log10(g0),m0],N_obs,m_vals,R1_tab) 
    up = True
    for ii in range(0,nL):
        m = mi[ii]
        dL = -2*(llhood2_marg(m,N_obs,m_vals,R1_tab)-maxL)
        if up:
            if dL>-sigmas**2.0:
                mlim = append(mlim,m)
                up = False
        else:
            if dL<-sigmas**2.0:
                mlim = append(mlim,m)
                up = True
    return mlim
#==============================================================================#


#==============================================================================#
# Old way of calculating the Mass estimate discovery limit, currently it fails 
# for most axion masses because the likelihood is so badly behaved
def MassErrorDiscoveryLimit_Old(err,m_vals,R0,R1_tab,m_DL_vals,gmin=1e-12,gmax=1e-7,ng=100,nL=100):
    nLike = nL
    nm = size(m_vals)
    n_DL = size(m_DL_vals)
    DLmerr = zeros(shape=(n_DL))
    g_vals = logspace(log10(gmin),log10(gmax),ng)
    for im in range(0,n_DL):
        m0 = m_DL_vals[im]
        for j in range(0,ng):
            g0 = g_vals[j]
            N_obs = InterpExpectedEvents(g0,m0,m_vals,R1_tab)
            Lmax = llhood2_marg(m0,N_obs,m_vals,R1_tab)
            D0 = -2*(llhood2_marg0(N_obs,R0)-Lmax)
            if D0<-1:
                X_in2 = [log10(g0),m0]
                maxL = llhood2(X_in2,N_obs,m_vals,R1_tab)
                Lu = -2*(llhood2_marg(m0*(1+0.5),N_obs,m_vals,R1_tab)-maxL)
                Ll = -2*(llhood2_marg(m0*(1-0.5),N_obs,m_vals,R1_tab)-maxL)
                if (Lu<-1)&(Ll<-1):
                    mlim = MassMeasurement(m0,g0,err,err,1,m_vals,R1_tab,nL=nLike)
                    if size(mlim)<=2:
                        rel_err = (max(mlim)-min(mlim))/(2*m0)
                        #rel_err = (max(log10(mlim))-min(log10(mlim)))/(2*log10(m0))                         
                        if rel_err<err:
                            DLmerr[im] = g0
                            break
                              
    return DLmerr
#==============================================================================#