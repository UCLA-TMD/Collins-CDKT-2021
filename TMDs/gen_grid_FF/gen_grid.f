      program master
      implicit none
      CHARACTER*72 fdss_pip_DIR, fdss_pip_GRIDS
      LOGICAL DIRECTORY_EXISTS
      INTEGER I_PIP
      CHARACTER*72 MKDIR
      integer NX, NQ
      PARAMETER (NX=47, NQ=24)
      real*8 zh
      real*8 Q2(NQ), QQ(NQ)
      real*8 XB(NX),Q
      integer I,J
      INTEGER PARTICLE_ID(-3:21)
      real*8 dd,uu,ss,sb,ub,db
      integer loops
      integer FINI14
      COMMON / FRAGINI14 / FINI14

C-----CHOOSE GRIDS TO EXPORT
      I_PIP = 1

      loops = 1

      FINI14 = 0

C.....Q**2 VALUES OF THE GRID.
       DATA Q2 / 1.d0, 1.25D0, 1.5D0, 2.5D0,
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
C......Z VALUES OF THE GRID.
       DATA XB /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

C-----PARTICLE ID NUMBERS: Based on Standard PDG ID Numbers:
C------https://pdg.lbl.gov/2006/reviews/pdf-files/montecarlo-web.pdf
      PARTICLE_ID( 1) =  1  ! d
      PARTICLE_ID(-1) = -1  ! dbar
      PARTICLE_ID( 2) =  2  ! u
      PARTICLE_ID(-2) = -2  ! ubar
      PARTICLE_ID( 3) =  3  ! s
      PARTICLE_ID(-3) = -3  ! sbar
      PARTICLE_ID(21) = 21  ! g

*-----Set Grid Directory Names
      fdss_pip_DIR = 'FDSS_PIP/'

*-----GRID FILE NAMES:
      fdss_pip_GRIDS = 'FDSS_PIP_0000.dat'

C-------HE

      IF (I_PIP.eq.1) THEN

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO
      INQUIRE(FILE=fdss_pip_DIR, EXIST=DIRECTORY_EXISTS)
      IF (.not.DIRECTORY_EXISTS) THEN
          CALL SYSTEM(MKDIR//fdss_pip_DIR)
      ENDIF
*-------Generate LHAPDF Grids
      OPEN(UNIT = 4, FILE = trim(fdss_pip_DIR)//fdss_pip_GRIDS)
      WRITE(4,*) 'PdfType: central'
      WRITE(4,*) 'Format: lhagrid1'
      WRITE(4,*) '---'
      WRITE(4,100) XB
      WRITE(4,101) QQ
      WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

      DO I = 1, NX
        zh = XB(I)
        DO J = 1, NQ
          Q = QQ(J)
          print *, zh,Q
          call CxFF(zh,Q,1,DD,UU,SS,SB,UB,DB)
          WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*DD,zh*UU,zh*SS,0d0
        ENDDO
      ENDDO


      WRITE(4,*) '---'
      CLOSE(4)
      ENDIF

*------FORMATTING

  100 FORMAT(47(E12.6E2,' '))
  101 FORMAT(51(E12.6E2,' '))
  102 FORMAT(' ',3(I2,' '),3(I1,' '),I2)
  103 FORMAT(7(E12.6E2,'  '))

      return
      end
