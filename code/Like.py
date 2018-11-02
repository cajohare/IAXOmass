from numpy import pi, sqrt, exp, zeros, size, shape, sinc, linspace, logspace
from numpy import log10, floor, log, interp, append, flipud
from scipy.integrate import cumtrapz, quad
from iminuit import minimize
from scipy.special import gammaln

#================================Like.py=======================================#
# Contains:
#==============================================================================#


def lnPF(Nob,Nex):
    return sum(Nob*log(Nex) - Nex)# - gammaln(Nob+1.0)) #factorial removed for speed

def InterpExpectedEvents(g,m,m_vals,R1_tab):
    nm = size(m_vals)
    m1 = m_vals[0]
    m2 = m_vals[-1]
    i1 = int(floor((log10(m)-log10(m_vals[0]))*(nm-1.0)/(log10(m_vals[-1])-log10(m_vals[0]))+1.0))-1
    i2 = i1+1
    N_exp = 1e40*(g**4.0)*((R1_tab[:,i1]*(m_vals[i2]-m)+R1_tab[:,i2]*(m - m_vals[i1]))/(m_vals[i2]- m_vals[i1]))
    return N_exp

def llhood2(X,N_obs,m_vals,R1_tab):
    m = X[1]
    g = 10.0**X[0]
    N_exp = InterpExpectedEvents(g,m,m_vals,R1_tab)
    LL = -1.0*lnPF(N_obs,N_exp)
    return LL
    
def llhood2_marg(m,N_obs,m_vals,R1_tab):
    N_exp_10 = InterpExpectedEvents(1e-10,m,m_vals,R1_tab)
    g0 = ((sum(N_obs)/sum(N_exp_10))**0.25)*1e-10
    LL = llhood2([log10(g0),m],N_obs,m_vals,R1_tab)
    return LL

def llhood2_marg0(N_obs,R0):
    N_exp = (sum(N_obs)/sum(R0))*R0
    LL = -1.0*lnPF(N_obs,N_exp)
    return LL
    
def llhood1(X,N_obs,R0):
    g = 10.0**X[0]
    N_exp = 1e40*(g**4.0)*R0
    LL = -1.0*lnPF(N_obs,N_exp)
    return LL

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

def MassDiscoveryLimit(m_vals,R1_tab,R0,m_DL_vals,gmin=1e-12,gmax=1e-7,ng=100):
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
                L1 = llhood2_marg0(N_obs,R0)
                #res = minimize(llhood1, X_in1, args=(N_obs,R0))
                #L1 = res.fun

                # Test statistic
                D12 = -2.0*(L2-L1) # significance for measuring mass
                if D12>9.0: # Median 3sigma detection -> D = 9
                    DL[im] = 10.0**(interp(9.0,[D12_prev,D12],[log10(g_prev),log10(g)]))
                    break
                g_prev = g # Reset for interpolation
                D12_prev = D12
    return DL


def CouplingDiscoveryLimit(m_vals,R1_tab,m_DL_vals,gmin=1e-12,gmax=1e-7,ng=100):
    nm = size(m_vals)
    n_DL = size(m_DL_vals)
    DL0 = zeros(shape=(n_DL))
    g_vals = logspace(log10(gmin),log10(gmax),ng)
    for im in range(0,n_DL):
        for j in range(0,ng):
            g = g_vals[j]
            m0 = m_DL_vals[im]
            N_obs = InterpExpectedEvents(g,m0,m_vals,R1_tab)
            if sum(N_obs)>4.9:
                DL0[im] = g
                break
    return DL0


def MassMeasurement(m0,g0,err_l,err_u,sig,m_vals,R1_tab,nL=100): 
    nm = size(m_vals)
    if sig==1:
        crit = -1
    elif sig==2:
        crit = -4
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
            if dL>crit:
                mlim = append(mlim,m)
                up = False
        else:
            if dL<crit:
                 mlim = append(mlim,m)
                 up = True
    return mlim


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

def MassErrorDiscoveryLimit(err,m_vals,R0,R1_tab,m_DL_vals,gmax=1e-9,ng=1000,nL=1000):
    itot = nL
    nm = size(m_vals)
    n_DL = size(m_DL_vals)
    DLmerr = zeros(shape=(n_DL))
    for im in range(0,n_DL):
        m0 = m_DL_vals[im]
        i0 = int(floor((log10(m0)-log10(m_vals[0]))*(nm-1.0)/(log10(m_vals[-1])-log10(m_vals[0]))+1.0))-1
        iupper = int(itot*i0/nm)
        ilower = itot-iupper
        mlower = flipud(logspace(log10(m_vals[1]),log10(m0*(1-err/2)),ilower))
        mupper = logspace(log10(m0*(1+err/2)),log10(m_vals[-2]),iupper)
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
            if (Lprof0<-1)&(Lprofend<-1):
                c = True
                for ii in range(0,ilower):
                    Lprof = -2*(llhood2_marg(mlower[ii],N_obs,m_vals,R1_tab) - Lmax)
                    if Lprof>-1:
                        c = False
                        break
                for ii in range(0,iupper):
                    Lprof = -2*(llhood2_marg(mupper[ii],N_obs,m_vals,R1_tab) - Lmax)
                    if Lprof>-1:
                        c = False
                        break    
                if c:
                    DLmerr[im] = g0
                    break
                    
    return DLmerr
                
def ConstantObsNumberLine(N_ob,mi,m_vals,R1_tab):
    ni = size(mi)
    g = zeros(shape=ni)
    for i in range(0,ni):
        N_exp_10 = InterpExpectedEvents(1e-10,mi[i],m_vals,R1_tab)
        g[i] = 1e-10*(N_ob/sum(N_exp_10))**0.25
    return g