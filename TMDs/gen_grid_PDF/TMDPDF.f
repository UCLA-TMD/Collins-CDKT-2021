c--------------------------------------------------------------
c     TMDPDF
c--------------------------------------------------------------
      subroutine TMD_PDF(x,Q,b,loops,logs,DD,UU,SS,SB,UB,DB)
Cf2py intent(in)  x
Cf2py intent(in)  Q
Cf2py intent(in)  b
Cf2py intent(in)  loops
Cf2py intent(in)  logs
Cf2py intent(out) DD
Cf2py intent(out) UU
Cf2py intent(out) SS
Cf2py intent(out) SB
Cf2py intent(out) UB
Cf2py intent(out) DB
      implicit none
      real*8 x,Q,b
      integer loops,logs
      real*8 DD,UU,SS,SB,UB,DB
      real*8 bmax,bstar,Q0,revo,Qini,ktw,kappa2,kappa1,c0
      real*8 UUi,DDi,UBi,DBi,SSi,SBi
      real*8 FNP
      integer init_flag
      common /in_flag/ init_flag

*-----b* prescription
      bmax = 1.5d0
      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      c0 = 1.122919d0
      Q0 = c0/bstar

*-----Perturbative evolution
      call kernel_q(bstar,Q0,Q,Q0,Q,logs,Revo) !logs =  1/2/3/4 for LL/NLL/NNLL/NNNLL

*-----OPE
      if(init_flag.ne.10) then
      call InitPDFsetByNameM(1,"HERAPDF20_NLO_VAR")
      call InitPDFM(1,0)
      init_flag = 10
      endif
      call CxPDF(x,Q0,loops,DDi,UUi,SSi,SBi,UBi,DBi) !0 = LO, 1 = NLO

*-----NP parameterization from 1406.3073v2
      ktw = 0.424d0
      kappa2 = 0.84d0
      Qini = dsqrt(2.4d0)
      kappa1 = ktw/4
      FNP = dexp( -(kappa1)*b*b-kappa2/2.*dlog(b/bstar)*dlog(Q/Qini) )

*-----TMDPDFs
      DD= DDi*FNP*REVO
      UU= UUi*FNP*REVO
      SS= SSi*FNP*REVO
      SB= SBi*FNP*REVO
      UB= UBi*FNP*REVO
      DB= DBi*FNP*REVO

      return
      end

c--------------------------------------------------------------
c     OPE
c--------------------------------------------------------------

      subroutine CxPDF(x,Q,loops,DD,UU,SS,SB,UB,DB)
Cf2py intent(in)  x
Cf2py intent(in)  Q
Cf2py intent(in)  loops
Cf2py intent(out) DD
Cf2py intent(out) UU
Cf2py intent(out) SS
Cf2py intent(out) SB
Cf2py intent(out) UB
Cf2py intent(out) DB
      implicit none
      real*8 x,Q
      integer loops
      real*8 DD,UU,SS,SB,UB,DB
      real*8 pdf(-6:6)
      real*8 DDn,UUn,SSn,SBn,UBn,DBn
      real*8 as,alphas,pref,pi
      data pi /3.14159265359d0/

      if (Q.le.1d0) then
          Q=1d0
      endif
     
      call evolvePDFM(1,x,Q,pdf)

      DD=pdf( 1)/x
      UU=pdf( 2)/x
      SS=pdf( 3)/x
      SB=pdf(-3)/x
      UB=pdf(-2)/x
      DB=pdf(-1)/x

      if (loops.eq.1) then
          as = alphas(Q)
          pref = as*4d0/3d0/4d0/pi*(-pi*pi/6d0)
          call PDF_NLO(x,Q,DDn,UUn,SSn,SBn,UBn,DBn)
          DD=DD*(1d0+pref)+DDn
          UU=UU*(1d0+pref)+UUn
          SS=SS*(1d0+pref)+SSn
          SB=SB*(1d0+pref)+SBn
          UB=UB*(1d0+pref)+UBn
          DB=DB*(1d0+pref)+DBn
      end if
 
      return
      end

      subroutine PDF_P(x,Q,DD,UU,SS,SB,UB,DB,GL)
Cf2py intent(in)  x
Cf2py intent(in)  Q
Cf2py intent(out) DD
Cf2py intent(out) UU
Cf2py intent(out) SS
Cf2py intent(out) SB
Cf2py intent(out) UB
Cf2py intent(out) DB
Cf2py intent(out) GL
      implicit none
      real*8 x,Q
      real*8 DD,UU,SS,SB,UB,DB,GL
      real*8 pdf(-6:6)

      if (Q.le.1d0) then
          Q=1d0
      endif
     
      call evolvePDFM(1,x,Q,pdf)

      DD=pdf( 1)/x
      UU=pdf( 2)/x
      SS=pdf( 3)/x
      SB=pdf(-3)/x
      UB=pdf(-2)/x
      DB=pdf(-1)/x
      GL=pdf( 0)/x
 
      return
      end

c-----------------------------------------------------------------------
c PDF at NLO
c-----------------------------------------------------------------------

      subroutine PDF_NLO(x,Q,DDN,UUN,SSN,SBN,UBN,DBN)
Cf2py intent(in)  x
Cf2py intent(in)  Q
Cf2py intent(out) DDN
Cf2py intent(out) UUN
Cf2py intent(out) SSN
Cf2py intent(out) SBN
Cf2py intent(out) UBN
Cf2py intent(out) DBN
      implicit none
      real*8 x,Q
      real*8 DDN,UUN,SSN,SBN,UBN,DBN
      real*8 cdd,cuu,css,csb,cub,cdb

      call conv_qgauss_PDF(50,x,Q,cdd,cuu,css,csb,cub,cdb)
      
      DDN=cdd
      UUN=cuu
      SSN=css
      SBN=csb
      UBN=cub
      DBN=cdb
      return
      
      end

c-----------------------------------------------------------------------
c NLO Coefficient at mu = mub and xi = mub**2
c-----------------------------------------------------------------------
      real*8 function qq(x,Q)
      implicit none
      real*8 x,Nc,pi,euler,Q,alfas,Cf
      real*8, external :: alphas
      
      pi=3.141592653589793d0
      
      Cf=1.3333333d0
      qq=alphas(Q)*Cf/2d0/pi*(1d0-x)/x
      
      end function
      
      real*8 function gq(x,Q)
      implicit none
      real*8 x,Nc,pi,euler,Q,alfas,TR
      real*8, external :: alphas
      
      pi=3.141592653589793d0
      
      TR=0.5d0
      gq=alphas(Q)*TR/pi*x*(1d0-x)/x
      
      end function

c-----------------------------------------------------------------------
c     Int_x_i^x_f func(x) dx using gaussian integration (for convolving)
c-----------------------------------------------------------------------
      subroutine conv_qgauss_PDF(n,x,Q,cdd,cuu,css,csb,cub,cdb)
Cf2py intent(in)  n
Cf2py intent(in)  x
Cf2py intent(in)  Q
Cf2py intent(out) cdd
Cf2py intent(out) cuu
Cf2py intent(out) css
Cf2py intent(out) csb
Cf2py intent(out) cub
Cf2py intent(out) cdb
      implicit none
      integer n
      real*8 x,Q
      real*8 cdd,cuu,css,csb,cub,cdb
      real*8 xi,xf,xn,x1,x2
      real*8 cddi,cuui,cssi,csbi,cubi,cdbi
      integer i


      cdd=0d0
      cuu=0d0
      css=0d0
      csb=0d0
      cub=0d0
      cdb=0d0

      if(n.le.1) then           ! same as n=1
         x1=x
         x2=1d0
         call conv_gauss_PDF(x1,x2,x,Q,cddi,cuui,cssi,csbi,cubi,cdbi)
         cdd=cddi
         cuu=cuui
         css=cssi
         csb=csbi
         cub=cubi
         cdb=cdbi
         return
      endif
	
      xi=x
      xf=1d0
      xn=(xf-xi)/float(n)
      x2=x
      Do 100 i=1,n
         x1=x2
         x2=x1+xn
         call conv_gauss_PDF(x1,x2,x,Q,cddi,cuui,cssi,csbi,cubi,cdbi)
         cdd=cdd+cddi
         cuu=cuu+cuui
         css=css+cssi
         csb=csb+csbi
         cub=cub+cubi
         cdb=cdb+cdbi
 100  continue
      return
      end

c-------------------------------------------------
      subroutine conv_gauss_PDF(xi,xf,x,Q,cdd,cuu,css,csb,cub,cdb)
Cf2py intent(in)  xi
Cf2py intent(in)  xf
Cf2py intent(in)  x
Cf2py intent(in)  Q
Cf2py intent(out) cdd
Cf2py intent(out) cuu
Cf2py intent(out) css
Cf2py intent(out) csb
Cf2py intent(out) cub
Cf2py intent(out) cdb
      implicit none
      real*8 xi,xf,x,Q
      real*8 cdd,cuu,css,csb,cub,cdb
      real*8 xm,xr,d,xd(8),w(8),eps
      real*8 DDp,UUp,SSp,SBp,UBp,DBp,GLp
      real*8 DDm,UUm,SSm,SBm,UBm,DBm,GLm
      real*8 gluon
      real*8 qq,gq
      integer j
      data eps /1.0d-25/
      data w
     1   / 0.02715 24594 11754 09485 17805 725D0,
     2     0.06225 35239 38647 89286 28438 370D0,
     3     0.09515 85116 82492 78480 99251 076D0,
     4     0.12462 89712 55533 87205 24762 822D0,
     5     0.14959 59888 16576 73208 15017 305D0,
     6     0.16915 65193 95002 53818 93120 790D0,
     7     0.18260 34150 44923 58886 67636 680D0,
     8     0.18945 06104 55068 49628 53967 232D0 /
      DATA XD
     1   / 0.98940 09349 91649 93259 61541 735D0,
     2     0.94457 50230 73232 57607 79884 155D0,
     3     0.86563 12023 87831 74388 04678 977D0,
     4     0.75540 44083 55003 03389 51011 948D0,
     5     0.61787 62444 02643 74844 66717 640D0,
     6     0.45801 67776 57227 38634 24194 430D0,
     7     0.28160 35507 79258 91323 04605 015D0,
     8     0.09501 25098 37637 44018 53193 354D0 /
      
      xm=0.5d0*(xf+xi)
      xr=0.5d0*(xf-xi)
      if (abs(xr).lt.eps) print *,
     >     'WARNING: Too high accuracy required for QGAUSS!'
 
      cdd=0d0
      cuu=0d0
      css=0d0
      csb=0d0
      cub=0d0
      cdb=0d0

      Do 100 j=1,8
         d=xr*xd(j)
         call PDF_P(x/(xm+d),Q,DDp,UUp,SSp,SBp,UBp,DBp,GLp)
         call PDF_P(x/(xm-d),Q,DDm,UUm,SSm,SBm,UBm,DBm,GLm)
         gluon=gq(xm+d,Q)*GLp+gq(xm-d,Q)*GLm
         cdd=cdd+w(j)*(qq(xm+d,Q)*DDp+qq(xm-d,Q)*DDm+gluon)
         cuu=cuu+w(j)*(qq(xm+d,Q)*UUp+qq(xm-d,Q)*UUm+gluon)
         css=css+w(j)*(qq(xm+d,Q)*SSp+qq(xm-d,Q)*SSm+gluon)
         csb=csb+w(j)*(qq(xm+d,Q)*SBp+qq(xm-d,Q)*SBm+gluon)
         cub=cub+w(j)*(qq(xm+d,Q)*UBp+qq(xm-d,Q)*UBm+gluon)
         cdb=cdb+w(j)*(qq(xm+d,Q)*DBp+qq(xm-d,Q)*DBm+gluon)
 100  continue
	
      cdd=xr*cdd
      cuu=xr*cuu
      css=xr*css
      csb=xr*csb
      cub=xr*cub
      cdb=xr*cdb
      return
      end

c--------------------------------------------------------------
c     alphas(q) -- NLO strong coupling constant (from CTEQ)
c--------------------------------------------------------------      
      function alphas(q)
      implicit none
      real*8 q,q2,lambda,lambda2,alphas,b0,b1,tt,pi,mb
      integer nf

      pi = atan(1d0)*4d0
      mb = 4.5d0
      nf=4
      b0=11d0-2d0/3d0*nf
      b1=102d0-38d0/3d0*nf
      lambda=0.227506d0

      if (q.le.0.23d0) then
        q=0.23d0
      endif
      q2=q*q
      lambda2=lambda*lambda
      tt=dlog(q2/lambda2)
      alphas=4d0*pi/(b0*tt)*(1d0-b1/(b0*b0)*dlog(tt)/tt)
      pi = atan(1d0)*4d0
      mb = 4.5d0
      if (q.le.mb) then
         nf=4
         lambda=0.326d0
      else
         nf=5
         lambda=0.226d0
      endif
      b0=11d0-2d0/3d0*nf
      b1=102d0-38d0/3d0*nf
*
      q2=q*q
      lambda2=lambda*lambda
      tt=dlog(q2/lambda2)
      alphas=4d0*pi/(b0*tt)*(1d0-b1/(b0*b0)*dlog(tt)/tt)
      return
      end

c--------------------------------------------------------------
c     alphas/(4pi) in terms of Lambda_QCD at NNLL
c--------------------------------------------------------------
      function aspi(Q,order)
      implicit none
      real*8 aspi,Q,x,xlog
      integer order
      real*8 b0,b1,b2,b1b0,b2b0,alphas
      real*8 c0,c02,pi,CF,CA,mc,mb,lambda
      integer ord,flag
      real*8 lambda3(3)
      real*8 lambda4(3)
      real*8 lambda5(3)
      real*8 Beta3(3)  
      real*8 Beta4(3)  
      real*8 Beta5(3)  
      data c0/1.122919d0/c02/1.260947d0/ 
      data pi/3.141592653589793d0/CF/1.333333d0/CA/3d0/
      data mc/1.275d0/mb/4.5d0/
      data lambda3/0.266d0, 0.266d0, 0.266d0/
      data lambda4/0.246d0, 0.246d0, 0.246d0/
      data lambda5/0.224d0, 0.224d0, 0.224d0/
      data Beta3/9.000000000000000d0, 64.000000000000000d0,
     >     643.8333333333334d0/
      data Beta4/8.333333333333333d0, 51.333333333333336d0,
     >     406.35185185185185d0/
      data Beta5/7.666666666666667d0, 38.666666666666664d0,
     >     180.90740740740742d0/

      flag = 2

      if(flag.eq.1) then

      ord = order
      ord = 3 

      if(Q.lt.mc) then
        lambda=lambda3(ord)
        b0=Beta3(1)
        b1=Beta3(2)
        b2=Beta3(3)
      elseif(Q.ge.mc.and.Q.lt.mb) then
        lambda=lambda4(ord)
        b0=Beta4(1)
        b1=Beta4(2)
        b2=Beta4(3)
      elseif(Q.ge.mb) then
        lambda=lambda5(ord)
        b0=Beta5(1)
        b1=Beta5(2)
        b2=Beta5(3)
      endif

      b1b0 = b1/(b0**2d0)
      b2b0 = b2/(b0*b0*b0)

      x = 2d0*dlog(Q/lambda)
      xlog = dlog(x)

      if(ord.eq.1) then 
        aspi = 1d0/b0/x
      elseif(ord.eq.2) then 
        aspi = 1d0/b0/x - b1b0/b0*(xlog/x/x)
      elseif(ord.eq.3) then 
        aspi = 1d0/b0/x - b1b0/b0*(xlog/x/x)
     >       + b1b0*b1b0/b0*(xlog*xlog-xlog-1d0)/(x*x*x)
     >       + b2b0/b0/(x*x*x)
      elseif(ord.gt.3) then
        write(*,*) 'Order in aspi larger than 3',ord
        stop
      endif


      elseif(flag.eq.2) then
c     calling the subroutine initialized in the main file
c     alphas(mu_r) returns the coupling constant (see alphas.f)
      aspi=alphas(Q)/(4d0*pi)

      endif

      return
      end function aspi



c---------------------------------------------------------
c     aspi_int(mui,muf,power,order) =
c     = \int_{mui}^{muf}dln(mu) aspi(mu,order)^power
c---------------------------------------------------------
      function aspi_int(mui,muf,power,order)
      real*8 aspi_int,mui,muf,aspi_int_aux,GAUSSP3
      integer power, pow, order, ord, steps
      common /order_aspi/ ord,pow

      ord = order
      pow = power
      steps = 4
      aspi_int = GAUSSP3(mui,muf,0.01d0)

      return
      end function aspi_int

c--------------------------------------------------------------------------
c               / B
c               |
c       value = | dx F(x), ACC is the percentage accuracy
c               |
c              / A
c--------------------------------------------------------------------------

      
      FUNCTION GAUSSP3(A,B,ACC)
      IMPLICIT NONE
      real*8 aspi_int_aux
      REAL*8 GAUSSP3,A,B,ACC,W(12),X(12),CONST,DELTA,AA,BB,Y
      REAL*8 C1,C2,S8,S16,U,SPLIT
      INTEGER I
      EXTERNAL F
      DATA CONST,SPLIT/1D-5,2.D0/
      DATA W/0.1012285,.2223810,.3137067,.3623838,.0271525,
     &     .0622535,0.0951585,.1246290,.1495960,.1691565,
     &     .1826034,.1894506/
      DATA X/0.9602899,.7966665,.5255324,.1834346,.9894009,
     &     .9445750,0.8656312,.7554044,.6178762,.4580168,
     &     .2816036,.0950125/

      DELTA=CONST*ABS(A-B)
      GAUSSP3=0.0D0
      AA=A
 5    Y=B-AA
      IF(ABS(Y).LE.DELTA) RETURN
 2    BB=AA+Y
      C1=0.5D0*(AA+BB)
      C2=C1-AA
      S8=0.0D0
      S16=0.0D0
      DO 1 I=1,4
      U=X(I)*C2
 1    S8=S8+W(I)*(aspi_int_aux(C1+U)+aspi_int_aux(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
 3    S16=S16+W(I)*(aspi_int_aux(C1+U)+aspi_int_aux(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(ABS(S16-S8).GT.ACC*ABS(S16)) GOTO 4
 8    GAUSSP3=GAUSSP3+S16
      AA=BB
      GOTO 5
 4    Y=Y/SPLIT
      
      IF(ABS(Y).GT.DELTA) THEN
         GOTO 2
      ELSE
         GOTO 8
      ENDIF
      
      RETURN
      END

c---------------------------------------------------------
c     Here I define the integrand to be input of qgauss
      function aspi_int_aux(mu)
      real*8 aspi_int_aux,aspi,mu
      integer ord,pow
      common /order_aspi/ ord,pow

      aspi_int_aux = (1d0/mu)
     &               *(aspi(mu,ord))**(float(pow))

      return
      end function aspi_int_aux
c---------------------------------------------------------


c---------------------------------------------------------
c     aspi_log_int(mui,muf,power,order) =
c     = \int_{mui}^{muf}dln(mu) aspi(mu,order)^power ln(muf^2/mu^2)
c---------------------------------------------------------
      function aspi_log_int(mui,muf,Q,power,order)
      real*8 aspi_log_int,mui,muf,aspi_log_int_aux,GAUSSP4,muff
      real*8 Q,QQ
      integer power, pow, order, ord, steps
      common /order_aspi2/ ord,pow,muff,QQ

      ord = order
      pow = power
      muff=muf
      QQ=Q
      steps = 4
      aspi_log_int = GAUSSP4(mui,muf,0.01d0)

      return
      end function aspi_log_int

c--------------------------------------------------------------------------
c               / B
c               |
c       value = | dx F(x), ACC is the percentage accuracy
c               |
c              / A
c--------------------------------------------------------------------------

      
      FUNCTION GAUSSP4(A,B,ACC)
      IMPLICIT NONE
      real*8 aspi_log_int_aux
      REAL*8 GAUSSP4,A,B,ACC,W(12),X(12),CONST,DELTA,AA,BB,Y
      REAL*8 C1,C2,S8,S16,U,SPLIT
      INTEGER I
      EXTERNAL F
      DATA CONST,SPLIT/1D-5,2.D0/
      DATA W/0.1012285,.2223810,.3137067,.3623838,.0271525,
     &     .0622535,0.0951585,.1246290,.1495960,.1691565,
     &     .1826034,.1894506/
      DATA X/0.9602899,.7966665,.5255324,.1834346,.9894009,
     &     .9445750,0.8656312,.7554044,.6178762,.4580168,
     &     .2816036,.0950125/

      DELTA=CONST*ABS(A-B)
      GAUSSP4=0.0D0
      AA=A
 5    Y=B-AA
      IF(ABS(Y).LE.DELTA) RETURN
 2    BB=AA+Y
      C1=0.5D0*(AA+BB)
      C2=C1-AA
      S8=0.0D0
      S16=0.0D0
      DO 1 I=1,4
      U=X(I)*C2
 1    S8=S8+W(I)*(aspi_log_int_aux(C1+U)+
     >            aspi_log_int_aux(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
 3    S16=S16+W(I)*(aspi_log_int_aux(C1+U)+
     >              aspi_log_int_aux(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(ABS(S16-S8).GT.ACC*ABS(S16)) GOTO 4
 8    GAUSSP4=GAUSSP4+S16
      AA=BB
      GOTO 5
 4    Y=Y/SPLIT
      
      IF(ABS(Y).GT.DELTA) THEN
         GOTO 2
      ELSE
         GOTO 8
      ENDIF
      
      RETURN
      END

c---------------------------------------------------------
c     Here I define the integrand to be input of qgauss
      function aspi_log_int_aux(mu)
      real*8 aspi_log_int_aux,aspi,mu,muff,QQ
      integer ord,pow
      common /order_aspi2/ ord,pow,muff,QQ
      real*8 factor
      common /fac2/ factor

      aspi_log_int_aux = (1d0/mu)
     &                   *(aspi(mu,ord))**(float(pow))
     &                   *2d0*dlog(QQ/factor/mu)

      return
      end function aspi_log_int_aux
c---------------------------------------------------------


c---------------------------------------------------------
c     FOR QUARK TMDs
c     \exp[ - \int_{mui}^{muf} dln(mu) ( A \ln(Q^2/mu^2) + B ) ]
c     Inputs: mui,muf
c     The integrand is taken consistently at LL,NLL,NNLL and NNNLL
c     depending on the value of the variable order: 1,2,3,4
c---------------------------------------------------------
      function ExpGamma_q(muii,muff,order)
      real*8 ExpGamma_q,mui,muf,muii,muff,value
      integer order
      real*8 aspi_int, aspi_log_int
      integer i
      real*8 c0,c02
      real*8 pi,CF,CA,mc,mb
      real*8 CuspQ3(4)      
      real*8 CuspQ4(4)      
      real*8 CuspQ5(4)      
      real*8 NCuspQ3(3)     
      real*8 NCuspQ4(3)     
      real*8 NCuspQ5(3)     
      common /NCuspQ/NCuspQ3,NCuspQ4,NCuspQ5
      data c0/1.122919d0/c02/1.260947d0/ 
      data pi/3.141592653589793d0/CF/1.333333d0/CA/3d0/
      data mc/1.275d0/mb/4.5d0/
      data CuspQ3/5.333333333333333d0, 48.69544319419009d0,
     >     618.2248693918798d0, 7848.82d0/
      data CuspQ4/5.333333333333333d0, 42.76951726826417d0,
     >     429.5065747522099d0, 4313.26d0/
      data CuspQ5/5.333333333333333d0, 36.84359134233824d0,
     >     239.20803319895987d0, 1553.06d0/

      if(muff.gt.muii) then
        muf=muff
        mui=muii
      elseif(muff.lt.muii) then
        muf=muii
        mui=muff
      else
        ExpGamma_q=1.0d0
        return
      end if

      if (mui.lt.mc.and.muf.lt.mc) then
        if(order.gt.0) then
          value = CuspQ3(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ3(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspQ3(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mc .and. muf.ge.mc
     &       .and. mui.lt.mb .and. muf.lt.mb) then
        if(order.gt.0) then
          value = CuspQ4(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ4(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspQ4(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mb .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspQ5(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ5(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspQ5(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.lt.mc .and. muf.ge.mc .and. muf.lt.mb) then
        if(order.gt.0) then
          value = CuspQ3(1)*aspi_log_int(mui,mc,muff,1,order)
     &            +CuspQ4(1)*aspi_log_int(mc,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ3(i)*aspi_log_int(mui,mc,muff,i,order)
     &          +CuspQ4(i)*aspi_log_int(mc,muf,muff,i,order)
     &          +NCuspQ3(i-1)*aspi_int(mui,mc,i-1,order)
     &          +NCuspQ4(i-1)*aspi_int(mc,muf,i-1,order)
          end do
        endif
      elseif(mui.lt.mc .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspQ3(1)*aspi_log_int(mui,mc,muff,1,order)
     &            +CuspQ4(1)*aspi_log_int(mc,mb,muff,1,order)
     &            +CuspQ5(1)*aspi_log_int(mb,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ3(i)*aspi_log_int(mui,mc,muff,i,order)
     &          +CuspQ4(i)*aspi_log_int(mc,mb,muff,i,order)
     &          +CuspQ5(i)*aspi_log_int(mb,muf,muff,i,order)
     &          +NCuspQ3(i-1)*aspi_int(mui,mc,i-1,order)
     &          +NCuspQ4(i-1)*aspi_int(mc,mb,i-1,order)
     &          +NCuspQ5(i-1)*aspi_int(mb,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mc .and. mui.lt.mb .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspQ4(1)*aspi_log_int(mui,mb,muff,1,order)
     &            +CuspQ5(1)*aspi_log_int(mb,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ4(i)*aspi_log_int(mui,mb,muff,i,order)
     &          +CuspQ5(i)*aspi_log_int(mb,muf,muff,i,order)
     &          +NCuspQ4(i-1)*aspi_int(mui,mb,i-1,order)
     &          +NCuspQ5(i-1)*aspi_int(mb,muf,i-1,order)
          end do
        endif
      end if

      if(muff.gt.muii) then
        ExpGamma_q=dexp(-1d0*value) 
      elseif(muff.lt.muii) then
        ExpGamma_q=dexp(+1.0d0*value) 
      else
        ExpGamma_q=1.0d0
        return
      end if

      return
      end function ExpGamma_q
c---------------------------------------------------------

c--------------------------------------------------------------
c     D function for quarks at fixed order,
c     with order=1,2,3,4 for LL,NLL,NNLL,NNNLL
c--------------------------------------------------------------
      subroutine Dterm_q(bt,mu,order,res)
Cf2py intent(in)  bt
Cf2py intent(in)  mu
Cf2py intent(in)  order
Cf2py intent(out) res
      implicit none
      real*8 bt,mu,res
      integer order,nf
      real*8 lperp,alf,alf2,alf3,aspi
      real*8 DtermQ11(3:5)
      real*8 DtermQ10(3:5)
      real*8 DtermQ22(3:5)
      real*8 DtermQ21(3:5)
      real*8 DtermQ20(3:5)
      real*8 DtermQ33(3:5)
      real*8 DtermQ32(3:5)
      real*8 DtermQ31(3:5)
      real*8 DtermQ30(3:5)
      common /DtermQ/DtermQ11,DtermQ10,DtermQ22,DtermQ21,DtermQ20,
     &               DtermQ33,DtermQ32,DtermQ31,DtermQ30
      real*8 c0,c02,pi,CF,CA,mc,mb
      data c0/1.122919d0/c02/1.260947d0/
      data pi/3.141592653589793d0/CF/1.333333d0/CA/3d0/
      data mc/1.275d0/mb/4.5d0/

      if(mu.lt.mc) then
        nf=3
      elseif(mu.ge.mc.and.mu.lt.mb) then
        nf=4
      else
        nf=5
      endif

      lperp = dlog(mu*mu*bt*bt/c02)
      alf  = aspi(mu,order)
      alf2 = alf*alf
      alf3 = alf2*alf

      res = 0d0
c        write(*,*) order,alf,res

      if(order.ge.2) then
        res = alf*(DtermQ11(nf)*lperp)
      end if

      if(order.ge.3) then
        res = res +
     &            alf2*(DtermQ22(nf)*lperp*lperp
     &                  + DtermQ21(nf)*lperp
     &                  + DtermQ20(nf))
      end if


      if(order.eq.4) then
        res = res +
     &            alf3*(DtermQ33(nf)*lperp*lperp*lperp
     &                  + DtermQ32(nf)*lperp*lperp
     &                  + DtermQ31(nf)*lperp + DtermQ30(nf))
      end if
c        write(*,*) 'D32=',nf,DtermQ32(nf)

      return
      end subroutine
c--------------------------------------------------------------


c---------------------------------------------------------
c     KERNEL FOR QUARK TMDs
c     The kernel is taken consistently at LL,NLL,NNLL and NNNLL
c     depending on the value of the variable order: 1,2,3,4
c---------------------------------------------------------
      subroutine kernel_q(b,mui,muf,Qi,Qf,order,res)
Cf2py intent(in)  b
Cf2py intent(in)  mui
Cf2py intent(in)  muf
Cf2py intent(in)  Qi
Cf2py intent(in)  Qf
Cf2py intent(in)  order
Cf2py intent(out) res
      real*8 b,mui,muf,Qi,Qf,res
      integer order
      real*8 ExpGamma_q,res1,res2
      real*8 factor,ffactor
      common /fac2/ ffactor

      factor = 1d0

      ffactor = factor

      res1 = ExpGamma_q(mui,muf,order)
      call Dterm_q(b,mui,order,res2)

      res = res1*dexp(-res2*2d0*dlog(Qf/Qi))

      return
      end subroutine
c---------------------------------------------------------







