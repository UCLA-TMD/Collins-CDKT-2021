import lhapdf
import numpy as np
from scipy.special import gamma,digamma,beta
from scipy.integrate import quad,fixed_quad
from mpmath import hyp2f1
from mpmath import psi as polygamma

from TMDs.Evolve.kernel_q import kernel_q, alphas
from TMDs.default_params import defaparams
from TMDs.Numerical.FBT import FBT

import warnings
warnings.filterwarnings("ignore")

class DISTRIB:

    def __init__(self,collFF,collPDF,params):
      self.ffpip = lhapdf.mkPDFs(collFF)
      self.pdfp  = lhapdf.mkPDFs(collPDF)

      self.ffpip = lhapdf.mkPDFs(collFF)
      self.pdfp  = lhapdf.mkPDFs(collPDF)

      self.fbt0 = FBT(0)
      self.fbt1 = FBT(1)

      #-- assign parameters
      #- collinear PDF
      self.Au = defaparams['f1']['MSTW2008']['NLO']['Au']
      self.n1 = defaparams['f1']['MSTW2008']['NLO']['n1']
      self.n2 = defaparams['f1']['MSTW2008']['NLO']['n2']
      self.eu = defaparams['f1']['MSTW2008']['NLO']['eu']
      self.gu = defaparams['f1']['MSTW2008']['NLO']['gu']
      self.Ad = defaparams['f1']['MSTW2008']['NLO']['Ad']
      self.n3 = defaparams['f1']['MSTW2008']['NLO']['n3']
      self.n4 = defaparams['f1']['MSTW2008']['NLO']['n4']
      self.ed = defaparams['f1']['MSTW2008']['NLO']['ed']
      self.gd = defaparams['f1']['MSTW2008']['NLO']['gd']
      self.AS = defaparams['f1']['MSTW2008']['NLO']['AS']
      self.dS = defaparams['f1']['MSTW2008']['NLO']['dS']
      self.nS = defaparams['f1']['MSTW2008']['NLO']['nS']
      self.eS = defaparams['f1']['MSTW2008']['NLO']['eS']
      self.gS = defaparams['f1']['MSTW2008']['NLO']['gS']
      self.AD = defaparams['f1']['MSTW2008']['NLO']['AD']
      self.nD = defaparams['f1']['MSTW2008']['NLO']['nD']
      self.nS = defaparams['f1']['MSTW2008']['NLO']['nS']
      self.gD = defaparams['f1']['MSTW2008']['NLO']['gD']
      self.dD = defaparams['f1']['MSTW2008']['NLO']['dD']
      self.Ap = defaparams['f1']['MSTW2008']['NLO']['Ap']
      self.dS = defaparams['f1']['MSTW2008']['NLO']['dS']
      self.npp= defaparams['f1']['MSTW2008']['NLO']['npp']
      self.eS = defaparams['f1']['MSTW2008']['NLO']['eS']
      self.gS = defaparams['f1']['MSTW2008']['NLO']['gS']
      self.Am = defaparams['f1']['MSTW2008']['NLO']['Am']
      self.dm = defaparams['f1']['MSTW2008']['NLO']['dm']
      self.nm = defaparams['f1']['MSTW2008']['NLO']['nm']
      self.x0 = defaparams['f1']['MSTW2008']['NLO']['x0']

      #- collinear fragmentation
      self.Nup = defaparams['g1']['DSSV']['N']['u+']
      self.aup = defaparams['g1']['DSSV']['a']['u+']
      self.bup = defaparams['g1']['DSSV']['b']['u+']
      self.gup = defaparams['g1']['DSSV']['g']['u+']
      self.nup = defaparams['g1']['DSSV']['e']['u+']
      self.Ndp = defaparams['g1']['DSSV']['N']['d+']
      self.adp = defaparams['g1']['DSSV']['a']['d+']
      self.bdp = defaparams['g1']['DSSV']['b']['d+']
      self.gdp = defaparams['g1']['DSSV']['g']['d+']
      self.ndp = defaparams['g1']['DSSV']['e']['d+']
      self.Nub = defaparams['g1']['DSSV']['N']['ub']
      self.aub = defaparams['g1']['DSSV']['a']['ub']
      self.bub = defaparams['g1']['DSSV']['b']['ub']
      self.gub = defaparams['g1']['DSSV']['g']['ub']
      self.nub = defaparams['g1']['DSSV']['e']['ub']
      self.Ndb = defaparams['g1']['DSSV']['N']['db']
      self.adb = defaparams['g1']['DSSV']['a']['db']
      self.bdb = defaparams['g1']['DSSV']['b']['db']
      self.gdb = defaparams['g1']['DSSV']['g']['db']
      self.ndb = defaparams['g1']['DSSV']['e']['db']
      self.Nss = defaparams['g1']['DSSV']['N']['ss']
      self.ass = defaparams['g1']['DSSV']['a']['ss']
      self.bss = defaparams['g1']['DSSV']['b']['ss']
      self.gss = defaparams['g1']['DSSV']['g']['ss']
      self.nss = defaparams['g1']['DSSV']['e']['ss']

      #- transversity
      self.NuuT = params[0]#defaparams['h1']['KPSY']['N']['uu']
      self.auuT = params[1] #defaparams['h1']['KPSY']['a']['uu']
      self.buuT = params[2] #defaparams['h1']['KPSY']['b']['uu']
      self.NddT = params[3] #defaparams['h1']['KPSY']['N']['uu']
      self.addT = params[4] #defaparams['h1']['KPSY']['a']['dd']
      self.bddT = params[5] #defaparams['h1']['KPSY']['b']['dd']
      self.NssT = params[6] #defaparams['h1']['KPSY']['N']['ss']
      self.assT = params[7] #defaparams['h1']['KPSY']['a']['ss']
      self.bssT = params[8] #defaparams['h1']['KPSY']['b']['ss']
      self.NsbT = params[9] #defaparams['h1']['KPSY']['N']['sb']
      self.asbT = params[10] #defaparams['h1']['KPSY']['a']['sb']
      self.bsbT = params[11] #defaparams['h1']['KPSY']['b']['sb']
      self.NdbT = params[12] #defaparams['h1']['KPSY']['N']['db']
      self.adbT = params[13] #defaparams['h1']['KPSY']['a']['db']
      self.bdbT = params[14] #defaparams['h1']['KPSY']['b']['db']
      self.NubT = params[15] #defaparams['h1']['KPSY']['N']['ub']
      self.aubT = params[16] #defaparams['h1']['KPSY']['a']['ub']
      self.bubT = params[17] #defaparams['h1']['KPSY']['b']['ub']

      self.Nfav =     params[18] #defaparams['H3']['KPSY']['N']['fav']
      self.Nufv =     params[19] #defaparams['H3']['KPSY']['N']['unf']
      self.alphafav = params[20] #defaparams['H3']['KPSY']['a']['fav']
      self.alphaufv = params[21] #defaparams['H3']['KPSY']['a']['unf']
      self.betafav =  params[22] #defaparams['H3']['KPSY']['b']['fav']
      self.betaufv =  params[23] #defaparams['H3']['KPSY']['N']['unf']

      self.NupDSS = defaparams['D1']['DSS']['N']['u+']
      self.NdpDSS = defaparams['D1']['DSS']['N']['d+']
      self.NubDSS = defaparams['D1']['DSS']['N']['ub']
      self.NddDSS = defaparams['D1']['DSS']['N']['dd']

      self.aupDSS = defaparams['D1']['DSS']['a']['u+']
      self.adpDSS = defaparams['D1']['DSS']['a']['d+']
      self.aubDSS = defaparams['D1']['DSS']['a']['ub']
      self.addDSS = defaparams['D1']['DSS']['a']['dd']

      self.alphaupDSS = defaparams['D1']['DSS']['a']['u+']
      self.alphadpDSS = defaparams['D1']['DSS']['a']['d+']
      self.alphaubDSS = defaparams['D1']['DSS']['a']['ub']
      self.alphaddDSS = defaparams['D1']['DSS']['a']['dd']

      self.betaupDSS = defaparams['D1']['DSS']['b']['u+']
      self.betadpDSS = defaparams['D1']['DSS']['b']['d+']
      self.betaubDSS = defaparams['D1']['DSS']['b']['ub']
      self.betaddDSS = defaparams['D1']['DSS']['b']['dd']

      self.gammaupDSS = defaparams['D1']['DSS']['g']['u+']
      self.gammadpDSS = defaparams['D1']['DSS']['g']['d+']
      self.gammaubDSS = defaparams['D1']['DSS']['g']['ub']
      self.gammaddDSS = defaparams['D1']['DSS']['g']['dd']

      self.deltaupDSS = defaparams['D1']['DSS']['d']['u+']
      self.deltadpDSS = defaparams['D1']['DSS']['d']['d+']
      self.deltaubDSS = defaparams['D1']['DSS']['d']['ub']
      self.deltaddDSS = defaparams['D1']['DSS']['d']['dd']


    #-- define collinear fragmentation function (FF)
    def CxFF(self,z,Q):
        #if Q<1.1:
        #    Q = 1.1
        D = self.ffpip[0].xfxQ( 1,z,Q)/z
        U = self.ffpip[0].xfxQ( 2,z,Q)/z
        S = self.ffpip[0].xfxQ( 3,z,Q)/z
        SB= self.ffpip[0].xfxQ(-3,z,Q)/z
        UB= self.ffpip[0].xfxQ(-2,z,Q)/z
        DB= self.ffpip[0].xfxQ(-1,z,Q)/z
        return 0,D,U,S,SB,UB,DB

    #-- define collinear parton distribution function (PDF)
    def CxPDF(self,x,Q):
        #if Q<1.1:
        #    Q = 1.1
        D = self.pdfp[0].xfxQ( 1,x,Q)/x
        U = self.pdfp[0].xfxQ( 2,x,Q)/x
        S = self.pdfp[0].xfxQ( 3,x,Q)/x
        SB= self.pdfp[0].xfxQ(-3,x,Q)/x
        UB= self.pdfp[0].xfxQ(-2,x,Q)/x
        DB= self.pdfp[0].xfxQ(-1,x,Q)/x
        return 0,D,U,S,SB,UB,DB

    #-- define unpolarized TMD PDF
    def TMDPDF(self,x,Q,b):
        bmax = 1.5
        bstar = b/np.sqrt(1.+(b/bmax)**2.)
        c0 = 1.122919
        Q0 = c0/bstar
        Revo = kernel_q(bstar,Q0,Q,Q0,Q,3)
        ggi,DDi,UUi,SSi,SBi,UBi,DBi = self.CxPDF(x,Q0)
        ktw = 0.424
        kappa2 = 0.84
        Qini = np.sqrt(2.4)
        kappa1 = ktw/4.
        FNP = np.exp( -(kappa1)*b*b-kappa2/2.*np.log(b/bstar)*np.log(Q/Qini) )
        DD= DDi*FNP*Revo
        UU= UUi*FNP*Revo
        SS= SSi*FNP*Revo
        SB= SBi*FNP*Revo
        UB= UBi*FNP*Revo
        DB= DBi*FNP*Revo
        return DD,UU,SS,SB,UB,DB

    def TMDFF(self,z,Q,b):
        bmax = 1.5
        bstar = b/np.sqrt(1.+(b/bmax)**2.)
        c0 = 1.122919
        Q0 = c0/bstar
        Revo = kernel_q(bstar,Q0,Q,Q0,Q,3)
        ggi,DDi,UUi,SSi,SBi,UBi,DBi = self.CxFF(z,Q0)
        ktw = 0.168
        kappa2 = 0.84
        Qini = np.sqrt(2.4)
        kappa1 = ktw/4./z/z
        FNP = np.exp( -(kappa1)*b*b-kappa2/2.*np.log(b/bstar)*np.log(Q/Qini) )
        DD= DDi*FNP*Revo/z/z
        UU= UUi*FNP*Revo/z/z
        SS= SSi*FNP*Revo/z/z
        SB= SBi*FNP*Revo/z/z
        UB= UBi*FNP*Revo/z/z
        DB= DBi*FNP*Revo/z/z
        return DD,UU,SS,SB,UB,DB

    def TransvuuQ0(self,n):
        pref = (self.auuT+self.buuT)**(self.auuT+self.buuT)/self.auuT**self.auuT/self.buuT**self.buuT
        val =  \
        -(self.Nub*self.NuuT*gamma(-1 + n + self.aub + self.auuT)*gamma(1 + self.bub + self.buuT))/(2.*gamma(n + self.aub + self.auuT + self.bub + self.buuT)) - \
        (self.Nub*self.NuuT*self.gub*gamma(-0.5 + n + self.aub + self.auuT)*gamma(1 + self.bub + self.buuT))/(2.*gamma(0.5 + n + self.aub + self.auuT + self.bub + self.buuT)) - \
        (self.Nub*self.NuuT*self.nub*gamma(n + self.aub + self.auuT)*gamma(1 + self.bub + self.buuT))/(2.*gamma(1 + n + self.aub + self.auuT + self.bub + self.buuT)) + \
        (self.Nup*self.NuuT*gamma(-1 + n + self.aup + self.auuT)*gamma(1 + self.bup + self.buuT))/(2.*gamma(n + self.aup + self.auuT + self.bup + self.buuT)) + \
        (self.Nup*self.NuuT*self.gup*gamma(-0.5 + n + self.aup + self.auuT)*gamma(1 + self.bup + self.buuT))/(2.*gamma(0.5 + n + self.aup + self.auuT + self.bup + self.buuT)) + \
        (self.Nup*self.NuuT*self.nup*gamma(n + self.aup + self.auuT)*gamma(1 + self.bup + self.buuT))/(2.*gamma(1 + n + self.aup + self.auuT + self.bup + self.buuT)) + \
        (self.Au*self.NuuT*gamma(-1 + n + self.auuT + self.n1)*gamma(1 + self.buuT + self.n2))/(2.*gamma(n + self.auuT + self.buuT + self.n1 + self.n2)) + \
        (self.Au*self.NuuT*self.eu*gamma(-0.5 + n + self.auuT + self.n1)*gamma(1 + self.buuT + self.n2))/(2.*gamma(0.5 + n + self.auuT + self.buuT + self.n1 + self.n2)) + \
        (self.Au*self.NuuT*self.gu*gamma(n + self.auuT + self.n1)*gamma(1 + self.buuT + self.n2))/(2.*gamma(1 + n + self.auuT + self.buuT + self.n1 + self.n2)) - \
        (self.Ap*self.NuuT*gamma(-1 + n + self.auuT + self.dS)*gamma(1 + self.buuT + self.npp))/(8.*gamma(n + self.auuT + self.buuT + self.dS + self.npp)) - \
        (self.Ap*self.NuuT*self.eS*gamma(-0.5 + n + self.auuT + self.dS)*gamma(1 + self.buuT + self.npp))/(8.*gamma(0.5 + n + self.auuT + self.buuT + self.dS + self.npp)) - \
        (self.Ap*self.NuuT*self.gS*gamma(n + self.auuT + self.dS)*gamma(1 + self.buuT + self.npp))/(8.*gamma(1 + n + self.auuT + self.buuT + self.dS + self.npp)) + \
        (self.AS*self.NuuT*gamma(-1 + n + self.auuT + self.dS)*gamma(1 + self.buuT + self.nS))/(8.*gamma(n + self.auuT + self.buuT + self.dS + self.nS)) + \
        (self.AS*self.NuuT*self.eS*gamma(-0.5 + n + self.auuT + self.dS)*gamma(1 + self.buuT + self.nS))/(8.*gamma(0.5 + n + self.auuT + self.buuT + self.dS + self.nS)) + \
        (self.AS*self.NuuT*self.gS*gamma(n + self.auuT + self.dS)*gamma(1 + self.buuT + self.nS))/(8.*gamma(1 + n + self.auuT + self.buuT + self.dS + self.nS)) - \
        (self.AD*self.NuuT*gamma(1 + self.buuT + self.nS)*gamma(-1 + n + self.auuT + self.nD))/(4.*gamma(n + self.auuT + self.buuT + self.nS + self.nD)) - \
        (self.AD*self.NuuT*self.gD*gamma(1 + self.buuT + self.nS)*gamma(-0.5 + n + self.auuT + self.nD))/(4.*gamma(0.5 + n + self.auuT + self.buuT + self.nS + self.nD)) + \
        (self.AD*self.NuuT*gamma(1 + self.buuT + self.nS)*gamma(n + self.auuT + self.nD))/(2.*gamma(1 + n + self.auuT + self.buuT + self.nS + self.nD)) - \
        (self.AD*self.NuuT*self.dD*gamma(1 + self.buuT + self.nS)*gamma(n + self.auuT + self.nD))/(4.*gamma(1 + n + self.auuT + self.buuT + self.nS + self.nD)) + \
        (self.AD*self.NuuT*self.gD*gamma(1 + self.buuT + self.nS)*gamma(0.5 + n + self.auuT + self.nD))/(2.*gamma(1.5 + n + self.auuT + self.buuT + self.nS + self.nD)) - \
        (self.AD*self.NuuT*gamma(1 + self.buuT + self.nS)*gamma(1 + n + self.auuT + self.nD))/(4.*gamma(2 + n + self.auuT + self.buuT + self.nS + self.nD)) + \
        (self.AD*self.NuuT*self.dD*gamma(1 + self.buuT + self.nS)*gamma(1 + n + self.auuT + self.nD))/(2.*gamma(2 + n + self.auuT + self.buuT + self.nS + self.nD)) - \
        (self.AD*self.NuuT*self.gD*gamma(1 + self.buuT + self.nS)*gamma(1.5 + n + self.auuT + self.nD))/(4.*gamma(2.5 + n + self.auuT + self.buuT + self.nS + self.nD)) - \
        (self.AD*self.NuuT*self.dD*gamma(1 + self.buuT + self.nS)*gamma(2 + n + self.auuT + self.nD))/(4.*gamma(3 + n + self.auuT + self.buuT + self.nS + self.nD))
        return pref*val

    def TransvddQ0(self,n):
        pref = (self.addT+self.bddT)**(self.addT+self.bddT)/self.addT**self.addT/self.bddT**self.bddT
        val = \
        -(self.Ndb*self.NddT*gamma(-1 + n + self.adb + self.addT)*gamma(1 + self.bdb + self.bddT))/(2.*gamma(n + self.adb + self.addT + self.bdb + self.bddT)) - \
        (self.Ndb*self.NddT*self.gdb*gamma(-0.5 + n + self.adb + self.addT)*gamma(1 + self.bdb + self.bddT))/(2.*gamma(0.5 + n + self.adb + self.addT + self.bdb + self.bddT)) - \
        (self.Ndb*self.NddT*self.Ndb*gamma(n + self.adb + self.addT)*gamma(1 + self.bdb + self.bddT))/(2.*gamma(1 + n + self.adb + self.addT + self.bdb + self.bddT)) + \
        (self.NddT*self.Ndp*gamma(-1 + n + self.addT + self.adp)*gamma(1 + self.bddT + self.bdp))/(2.*gamma(n + self.addT + self.adp + self.bddT + self.bdp)) + \
        (self.NddT*self.Ndp*self.gdp*gamma(-0.5 + n + self.addT + self.adp)*gamma(1 + self.bddT + self.bdp))/(2.*gamma(0.5 + n + self.addT + self.adp + self.bddT + self.bdp)) + \
        (self.NddT*self.Ndp*self.Ndp*gamma(n + self.addT + self.adp)*gamma(1 + self.bddT + self.bdp))/(2.*gamma(1 + n + self.addT + self.adp + self.bddT + self.bdp)) + \
        (self.Ad*self.NddT*gamma(-1 + n + self.addT + self.n3)*gamma(1 + self.bddT + self.n4))/(2.*gamma(n + self.addT + self.bddT + self.n3 + self.n4)) + \
        (self.Ad*self.NddT*self.ed*gamma(-0.5 + n + self.addT + self.n3)*gamma(1 + self.bddT + self.n4))/(2.*gamma(0.5 + n + self.addT + self.bddT + self.n3 + self.n4)) + \
        (self.Ad*self.NddT*self.gd*gamma(n + self.addT + self.n3)*gamma(1 + self.bddT + self.n4))/(2.*gamma(1 + n + self.addT + self.bddT + self.n3 + self.n4)) - \
        (self.Ap*self.NddT*gamma(-1 + n + self.addT + self.dS)*gamma(1 + self.bddT + self.npp))/(8.*gamma(n + self.addT + self.bddT + self.dS + self.npp)) - \
        (self.Ap*self.NddT*self.eS*gamma(-0.5 + n + self.addT + self.dS)*gamma(1 + self.bddT + self.npp))/(8.*gamma(0.5 + n + self.addT + self.bddT + self.dS + self.npp)) - \
        (self.Ap*self.NddT*self.gS*gamma(n + self.addT + self.dS)*gamma(1 + self.bddT + self.npp))/(8.*gamma(1 + n + self.addT + self.bddT + self.dS + self.npp)) + \
        (self.AS*self.NddT*gamma(-1 + n + self.addT + self.dS)*gamma(1 + self.bddT + self.nS))/(8.*gamma(n + self.addT + self.bddT + self.dS + self.nS)) + \
        (self.AS*self.NddT*self.eS*gamma(-0.5 + n + self.addT + self.dS)*gamma(1 + self.bddT + self.nS))/(8.*gamma(0.5 + n + self.addT + self.bddT + self.dS + self.nS)) + \
        (self.AS*self.NddT*self.gS*gamma(n + self.addT + self.dS)*gamma(1 + self.bddT + self.nS))/(8.*gamma(1 + n + self.addT + self.bddT + self.dS + self.nS)) + \
        (self.AD*self.NddT*gamma(1 + self.bddT + self.nS)*gamma(-1 + n + self.addT + self.nD))/(4.*gamma(n + self.addT + self.bddT + self.nS + self.nD)) + \
        (self.AD*self.NddT*self.gD*gamma(1 + self.bddT + self.nS)*gamma(-0.5 + n + self.addT + self.nD))/(4.*gamma(0.5 + n + self.addT + self.bddT + self.nS + self.nD)) - \
        (self.AD*self.NddT*gamma(1 + self.bddT + self.nS)*gamma(n + self.addT + self.nD))/(2.*gamma(1 + n + self.addT + self.bddT + self.nS + self.nD)) + \
        (self.AD*self.NddT*self.dD*gamma(1 + self.bddT + self.nS)*gamma(n + self.addT + self.nD))/(4.*gamma(1 + n + self.addT + self.bddT + self.nS + self.nD)) - \
        (self.AD*self.NddT*self.gD*gamma(1 + self.bddT + self.nS)*gamma(0.5 + n + self.addT + self.nD))/(2.*gamma(1.5 + n + self.addT + self.bddT + self.nS + self.nD)) + \
        (self.AD*self.NddT*gamma(1 + self.bddT + self.nS)*gamma(1 + n + self.addT + self.nD))/(4.*gamma(2 + n + self.addT + self.bddT + self.nS + self.nD)) - \
        (self.AD*self.NddT*self.dD*gamma(1 + self.bddT + self.nS)*gamma(1 + n + self.addT + self.nD))/(2.*gamma(2 + n + self.addT + self.bddT + self.nS + self.nD)) + \
        (self.AD*self.NddT*self.gD*gamma(1 + self.bddT + self.nS)*gamma(1.5 + n + self.addT + self.nD))/(4.*gamma(2.5 + n + self.addT + self.bddT + self.nS + self.nD)) + \
        (self.AD*self.NddT*self.dD*gamma(1 + self.bddT + self.nS)*gamma(2 + n + self.addT + self.nD))/(4.*gamma(3 + n + self.addT + self.bddT + self.nS + self.nD))
        return pref*val

    def Pqqh1(self,n):
        return 4./3.*(-2.*(digamma(n)+1./n+np.euler_gamma)+1.5)

    def DGLAP(self,n,Q):
        Q0 = 1.
        mc = 1.275
        mb = 4.5
        if Q<= mc:
            nf = 3.
            b0 = 11.-2./3.*nf
            evo  = (alphas(Q)/alphas(Q0))**(-self.Pqqh1(n)/b0)
        elif Q<= mb:
            nf = 3.
            b0 = 11.-2./3.*nf
            evoc = (alphas(mc)/alphas(Q0))**(-self.Pqqh1(n)/b0)
            nf = 4.
            b0 = 11.-2./3.*nf
            evob = (alphas(Q )/alphas(mc))**(-self.Pqqh1(n)/b0)
            evo  = evoc*evob
        elif Q > mb:
            nf = 3.
            b0 = 11.-2./3.*nf
            evoc = (alphas(mc)/alphas(Q0))**(-self.Pqqh1(n)/b0)
            nf = 4.
            b0 = 11.-2./3.*nf
            evob = (alphas(mb)/alphas(mc))**(-self.Pqqh1(n)/b0)
            nf = 5.
            b0 = 11.-2./3.*nf
            evoQ = (alphas(Q )/alphas(mb))**(-self.Pqqh1(n)/b0)
            evo  = evoc*evob*evoQ
        return evo

    #-- evolved Transversity
    def TransvuuQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        return self.TransvuuQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))

    def TransvddQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        return self.TransvddQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))

    def TransvssQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        return self.TransvssQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))

    def TransvubQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        return self.TransvubQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))

    def TransvdbQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        return self.TransvdbQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))

    def TransvsbQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        return self.TransvsbQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))

    def TransvQ(self,fl,n,Q):

        if fl == 'uu':
            TransvQ0 = self.TransvuuQ0(n)
        elif fl == 'dd':
            TransvQ0 = self.TransvddQ0(n)
        elif fl == 'ss':
            TransvQ0 = self.TransvssQ0(n)
        elif fl == 'ub':
            TransvQ0 = self.TransvubQ0(n)
        elif fl == 'db':
            TransvQ0 = self.TransvdbQ0(n)
        elif fl == 'sb':
            TransvQ0 = self.TransvsbQ0(n)

        evo = self.DGLAP(n,Q)

        return self.TransvuuQ0(n)*evo*(1.-2.*4./3./np.pi*alphas(Q))


    #--- transversity in x and Q
    def TransvuuxQ(self,x,Q):
        phi = 3.*np.pi/4.
        c = 2.
        integrand = lambda z: 1./np.pi*np.imag(np.exp(1j*phi)*x**(-c-z*np.exp(1j*phi))*self.TransvuuQ(c+z*np.exp(1j*phi),Q))
        return quad(integrand,0.,20.,epsabs = 0.,epsrel = 0.05)[0]

    def TransvddxQ(self,x,Q):
        phi = 3.*np.pi/4.
        c = 2.
        integrand = lambda z: 1./np.pi*np.imag(np.exp(1j*phi)*x**(-c-z*np.exp(1j*phi))*self.TransvddQ(c+z*np.exp(1j*phi),Q))
        return quad(integrand,0.,20.,epsabs = 0.,epsrel = 0.05)[0]


    def h1buu(self,x,Q,b):
        bmax = 1.5
        bstar = b/np.sqrt(1.+(b/bmax)**2.)
        c0 = 1.122919
        Q0 = c0/bstar
        Revo = kernel_q(bstar,Q0,Q,Q0,Q,3)
        UUi = self.TransvuuxQ(x,Q0)
        ktw = 0.424
        kappa2 = 0.84
        Qini = np.sqrt(2.4)
        kappa1 = ktw/4.
        FNP = np.exp( -(kappa1)*b*b-kappa2/2.*np.log(b/bstar)*np.log(Q/Qini) )
        UU= UUi*FNP*Revo
        return UU

    def h1bdd(self,x,Q,b):
        bmax = 1.5
        bstar = b/np.sqrt(1.+(b/bmax)**2.)
        c0 = 1.122919
        Q0 = c0/bstar
        Revo = kernel_q(bstar,Q0,Q,Q0,Q,3)
        DDi = self.TransvddxQ(x,Q0)
        ktw = 0.424
        kappa2 = 0.84
        Qini = np.sqrt(2.4)
        kappa1 = ktw/4.
        FNP = np.exp( -(kappa1)*b*b-kappa2/2.*np.log(b/bstar)*np.log(Q/Qini) )
        DD= DDi*FNP*Revo
        return DD

    def ColnsCDKTddQ0(self,n):
        return \
        (self.NddDSS*self.Nufv*gamma(self.addDSS + n + self.alphaufv)*(self.gammaddDSS*gamma(1 + self.addDSS + n + self.alphaufv + self.betaddDSS + self.betaufv)*gamma(1 + self.betaddDSS + self.betaufv + self.deltaddDSS) + \
        gamma(1 + self.betaddDSS + self.betaufv)*gamma(1 + self.addDSS + n + self.alphaufv + self.betaddDSS + self.betaufv + self.deltaddDSS)))/ \
        ((beta(2 + self.alphaddDSS,1 + self.betaddDSS) + self.gammaddDSS*beta(2 + self.alphaddDSS,1 + self.betaddDSS + self.deltaddDSS))* \
        gamma(1 + self.addDSS + n + self.alphaufv + self.betaddDSS + self.betaufv)*gamma(1 + self.addDSS + n + self.alphaufv + self.betaddDSS + self.betaufv + self.deltaddDSS))

    def ColnsCDKTuuQ0(self,n):
        return \
        self.Nfav*gamma(1 + self.betafav)*(-((self.NubDSS*gamma(self.aubDSS + n + self.alphafav)*(hyp2f1(self.aubDSS + n + self.alphafav,-self.betaubDSS,1 + self.aubDSS + n + self.alphafav + self.betafav,1) + \
        self.gammaubDSS*hyp2f1(self.aubDSS + n + self.alphafav,-self.betaubDSS - self.deltaubDSS,1 + self.aubDSS + n + self.alphafav + self.betafav,1)))/ \
        ((beta(2 + self.alphaubDSS,1 + self.betaubDSS) + self.gammaubDSS*beta(2 + self.alphaubDSS,1 + self.betaubDSS + self.deltaubDSS))*gamma(1 + self.aubDSS + n + self.alphafav + self.betafav))) + \
        (self.NupDSS*gamma(self.aupDSS + n + self.alphafav)*(hyp2f1(self.aupDSS + n + self.alphafav,-self.betaupDSS,1 + self.aupDSS + n + self.alphafav + self.betafav,1) + \
        self.gammaupDSS*hyp2f1(self.aupDSS + n + self.alphafav,-self.betaupDSS - self.deltaupDSS,1 + self.aupDSS + n + self.alphafav + self.betafav,1)))/ \
        ((beta(2 + self.alphaupDSS,1 + self.betaupDSS) + self.gammaupDSS*beta(2 + self.alphaupDSS,1 + self.betaupDSS + self.deltaupDSS))*gamma(1 + self.aupDSS + n + self.alphafav + self.betafav)))


    def H1TuuQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        CF =  4./3.
        As = alphas(Q)
        return self.ColnsCDKTuuQ0(n)*evo*(1.+As/np.pi*(-2*CF*polygamma(1,n+1))-2.*4./3./np.pi*As)

    def H1TddQ(self,n,Q):
        evo = self.DGLAP(n,Q)
        CF =  4./3.
        As = alphas(Q)
        return self.ColnsCDKTddQ0(n)*evo*(1.+As/np.pi*(-2*CF*polygamma(1,n+1))-2.*4./3./np.pi*As)

    def H1TuuxQ(self,x,Q):
        phi = 3.*np.pi/4.
        c = 3.
        integrand = lambda z: 1./np.pi*np.imag(np.exp(1j*phi)*x**(-c-z*np.exp(1j*phi))*self.H1TuuQ(c+z*np.exp(1j*phi),Q))
        return quad(integrand,0.,20.,epsabs = 0.,epsrel = 0.05)[0]

    def H1TddxQ(self,x,Q):
        phi = 3.*np.pi/4.
        c = 3.
        integrand = lambda z: 1./np.pi*np.imag(np.exp(1j*phi)*x**(-c-z*np.exp(1j*phi))*self.H1TddQ(c+z*np.exp(1j*phi),Q))
        return quad(integrand,0.,20.,epsabs = 0.,epsrel = 0.05)[0]

    def H1Tuub(self,z,Q,b):
        bmax = 1.5
        bstar = b/np.sqrt(1.+(b/bmax)**2.)
        c0 = 1.122919
        Q0 = c0/bstar
        Revo = kernel_q(bstar,Q0,Q,Q0,Q,3)
        UUi = self.H1TuuxQ(z,Q0)
        ptw = (0.042-0.0236)*4
        kappa2 = 0.84
        Qini = np.sqrt(2.4)
        kappa1 = ptw/4.
        FNP = np.exp( -(kappa1)*b*b/z/z-kappa2/2.*np.log(b/bstar)*np.log(Q/Qini) )
        UU= UUi*FNP*Revo/z/z
        return UU

    def H1Tddb(self,z,Q,b):
        bmax = 1.5
        bstar = b/np.sqrt(1.+(b/bmax)**2.)
        c0 = 1.122919
        Q0 = c0/bstar
        Revo = kernel_q(bstar,Q0,Q,Q0,Q,3)
        DDi = self.H1TddxQ(z,Q0)
        ptw = (0.042-0.0236)*4
        kappa2 = 0.84
        Qini = np.sqrt(2.4)
        kappa1 = ptw/4.
        FNP = np.exp( -(kappa1)*b*b/z/z-kappa2/2.*np.log(b/bstar)*np.log(Q/Qini) )
        DD= DDi*FNP*Revo/z/z
        return DD

    #-- SIDIS structure functions


    def FUU_b_SIDIS(self,b,x,z,Q):
        FDD,FUU,FSS,FSB,FUB,FDB = self.TMDPDF(x,Q,b)
        DDD,DUU,DSS,DSB,DUB,DDB = self.TMDFF(z,Q,b)
        als = alphas(Q)
        CF = 4./3.
        HQ = 1#+als*CF/2./np.pi*(-8.)
        eu2 = 4./9.
        ed2 = 1./9.
        return HQ*(eu2*(FUU*DUU+FUB*DUB)+ed2*(FDD*DDD+FDB*DDB+FSS*DSS+FSB*DSB))
    def FUU_p_SIDIS(self,PhT,x,z,Q):
        FUUb = np.vectorize(lambda b: b*self.FUU_b_SIDIS(b,x,z,Q))
        return self.fbt0.fbt(FUUb, PhT/z, 20,Q)

    def FUT_b_SIDIS(self,b,x,z,Q):
        FDD,FUU = self.h1bdd(x,Q,b) ,self.h1buu(x,Q,b)
        DDD,DUU = self.H1Tddb(z,Q,b),self.H1Tddb(z,Q,b)
        als = alphas(Q)
        CF = 4./3.
        HQ = 1#+als*CF/2./np.pi*(-8.)
        eu2 = 4./9.
        ed2 = 1./9.
        return HQ*(eu2*(FUU*DUU)+ed2*(FDD*DDD))
    def FUT_p_SIDIS(self,PhT,x,z,Q):
        FUTb = np.vectorize(lambda b: b*b*self.FUT_b_SIDIS(b,x,z,Q))
        return self.fbt1.fbt(FUTb, PhT/z, 20,Q)

    def Collins_Asym(self,PhT,x,z,Q):
        FUU = self.FUU_p_SIDIS(PhT,x,z,Q)
        FUT = self.FUT_p_SIDIS(PhT,x,z,Q)

        return FUT/FUU
