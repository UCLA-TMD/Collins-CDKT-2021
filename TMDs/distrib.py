import lhapdf
import numpy as np
from Evolve.kernel_q import kernel_q, alphas

class DISTRIB:

    def __init__(self,collFF,collPDF):
      self.ffpip = lhapdf.mkPDFs(collFF)
      self.pdfp  = lhapdf.mkPDFs(collPDF)

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
