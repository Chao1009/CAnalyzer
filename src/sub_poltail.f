!     These subroutines are based on Seonho's version of POLRAD (elastic part)
!     Hacked it slightly to handle both He-3 and proton. Also updated He-3 Form factors.
!-----------------------------------------------------------------------------------
!     Downloaded on 12/07/06 from Karl Slifer's website by Jaideep Singh
!     changed all comment line markers to !
      BLOCK DATA
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)

       !DATA AMT/2.8094D0/,     TARA/3D0/, TARZ/2D0/
       !DATA AMT/0.93827231D0/, TARA/1D0/, TARZ/1D0/                 ! CPROT


      DATA AMM/2.7928456D0/
      DATA AMN/-1.913148D0/
      DATA CHBAR/.197328D0/
      DATA BARN/.389379D6/
      DATA AML/.511000D-3/
      DATA AML2/.261112D-6/
      DATA AL2/.522240D-6/
      DATA FERMOM/.164D0/     !  - ? FERMI MOMENTUM IN 3HE
      DATA ISF20/4/
      DATA PI/3.1415926D0/,PI2/9.869604D0/,ALFA/.729735D-2/
      DATA AMC2/1.151857D0/
      DATA AMP/.938272D0/,AMH/.938272D0/
      DATA I2/1,1,1,2,3,3,1,2/,I1/3,3,4,4,3,3,3,3/
      END

!-----------------------------------------------------------------------------------


      !-----------------------------------------------------------------
      SUBROUTINE POLSIG_EL(E0,EP,TH_RAD,IPOL,SIGMA_BORN,SIGMA_RAD)
!
!     CALCULATION OF THE ELASTIC TAIL
!
!     INPUT
!
!       E0 = INCIDENT ELECTRON ENERGY IN MEV
!       EP = SCATTERED ELECTRON ENERGY IN MEV
!       TH_RAD = SCATTERING ANGLE IN RADIANS
!
!       IPOL = 0 : PARALLEL, UNPOLARIZED
!              1 : PARALLEL, POLARIZED
!              2 : PERPENDICULAR, UNPOLARIZED
!              3 : PERPENDICULAR, POLARIZED
!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 E0,EP,TH_RAD,SIGMA_BORN,SIGMA_RAD
      INTEGER IPOL

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &           FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/CHOI/      W2_ELAS
      COMMON/SOFTPHOTON/EXTAI2,EXTAI3

      COMMON/SXY/S   ,X   ,SX  ,SXP ,Y   ,YM,W2   ,ALS  ,ALX,ALM,ALY,
     &           SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS ,YS ,TPL,TMI

      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)

      COMMON/TAIL/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH
      !------------------------------------------------------------------
      Q2  = 4.0*E0*EP*SIN(0.5*TH_RAD)**2 * 1.0E-6
      XNU = (E0-EP)                      * 1.0E-3

!$$$  W2  = 0.938272*(0.938272 + 2.0*XNU) - Q2
      W2  = AMT*(AMT + 2.0*XNU) - Q2

      IF (W2.LE.0.0) THEN
         SIGMA_BORN = 0.0
         SIGMA_RAD  = 0.0
         RETURN
      ENDIF

!$$$  W       = SQRT(W2)
!$$$  XJACOB2 = 2.0*EP*W/(2.0*0.938272*EP*1.0E-3+Q2)

      XS      =              Q2 / (2.0*AMT*XNU)
      YS      =              XNU/ (E0*1.0E-3)
      XJACOB  = 2.0D0*PI*AMH*YS / (1.0-YS)

      XJACOB  = 1.0D-3/XJACOB                    ! 1.0D-3 CONVERTS FROM NB/GEV TO NB/MEV



!$$$  XJACOB  = 1.0                              ! TO COMPARE WITH THE ORIGINAL POLRAD20

      Y_ELAS  = 1.0/(1.0 + AMT/(2.0* E0*1.0E-3 * SIN(0.5*TH_RAD)**2) )


!$$$  W2_ELAS = AMH**2 - 2.0*(AMT - AMH)*Y_ELAS*E0*1.0E-3 ! USE PROTON MASS FOR W2
      W2_ELAS = AMT**2                                    ! USE HE3 MASS FOR W2

      SNUC    = 2.0*AMH*SQRT((E0*1.0E-3)**2+AML2)
      Y       = SNUC*XS*YS*(AMT/AMH)                      ! Q^2

      YMA     = 1D0/(1D0+AMT**2*XS/SNUC) ! MAXIMUM Y FOR A GIVEN E AND X. IT OCCURRS
                                         ! AT BACK SCATTERING, THETA = 180 DEGREES

      !AMC2    = (MP + M_PI)**2                    ! PION PRODUCTION THRESHOLD

      YMI     = (AMC2-AMP**2) / (SNUC*(1D0-XS))    ! MINIMUM Y FOR PION PRODUCTION

!$$$  IF (YMI.LT.0.0) YMI = 1.0

      IF(YS.GT.YMA .OR. YS.LT.Y_ELAS .OR. XS.LT.0D0 .OR. XS.GT.2.0) THEN
         SIGMA_BORN = 0.0
         SIGMA_RAD  = 0.0
         RETURN
      ENDIF

      CALL CONKIN(SNUC,AMT,IPOL)

!
!-----DELTA IS FACTORIZING PART OF VIRTUAL AND REAL LEPTONIC BREMSSTRAHLUNG
!
      IF (YS.GE.Y_ELAS) THEN
         CALL DELTAS(DELTA,DELINF,TR,FACTOR1,DEL_SUB)
      ELSE
         DELTA   = 0.0
         DELINF  = 0.0
         TR      = 0.0
         FACTOR1 = 0.0
         DEL_SUB = 0.0
      ENDIF

      ISF3 = 1

      IF (IPOL.EQ.0 .OR. IPOL.EQ.2) THEN  ! CROSS SECTION FOR UNPOLARIZED HADRON TARGET
         UN = 1.0               ! UN = 1 CALCULATES F1 AND F2 STRUCTURE FUNCTION
         PL = 1.0               ! LEPTON POLARIZATION
         PN = 0.0               ! UNPOLARIZED TARGET, G1 = G2 = 0
         QN = 0.0               ! QN = 0 FOR THE HADRONIC TENSOR FOR SPIN 1/2 PARTICLE
         ISF1 = 1
         ISF2 = 2
      ELSE IF (IPOL.EQ.1.OR.IPOL.EQ.3) THEN ! DIFFERENCE BETWEEN TWO HADRON
                                            ! POLARIZATION DIRECTIONS
         UN = 0.0               ! UN = 0 MEANS F1 = F2 = 0
         PL = 1.0               ! LEPTON POLARIZATION
         PN = 1.0               ! PN DEFINES HADRON POLARIZATION, G1 AND G2 NON-ZERO
         QN = 0.0
         ISF1 = 3
         ISF2 = 4
      ELSE
         STOP "IPOL Unknown"
      ENDIF

      EXTAI2 = 1.0D0
      EXTAI3 = 1.0D0

      CALL BORNIN(SIB)



      EXTAI2 = ( (SX-Y/TARA)**2 /S /(S-Y/TARA) )**TR
      EXTAI3 = (      (SX-Y)**2 /S /(S-Y)      )**TR    ! NOT USED NOW

      TAIL       = AN  * ALFA/PI * TAIL_INTEG_EL(TAMIN,TAMAX)
      SIGMA_BORN = SIB * XJACOB


      IF(TARZ.EQ.1. .AND. TARA.EQ.1.) THEN ! CPROT
        EXTAI1     = EXP(ALFA/PI*DELINF)
        SIGMA_RAD  = (SIB * EXTAI1 * (1.+ ALFA/PI * (DELTA-DELINF) ) + TAIL) * XJACOB ! C$$$
      ELSEIF(TARZ.EQ.2. .AND. TARA.EQ.3. ) THEN
        SIGMA_RAD  = (SIB * FACTOR1* (1.+ ALFA/PI *         DEL_SUB) + TAIL) * XJACOB
      ELSE
        STOP "Problem 1"
      ENDIF

      IF (IPOL.EQ.1 .OR. IPOL.EQ.3) THEN
         SIGMA_BORN = -SIGMA_BORN
         SIGMA_RAD  = -SIGMA_RAD
      ENDIF

      RETURN
      END

!-----------------------------------------------------------------------------------
      REAL*8 FUNCTION TAIL_INT_OVER_R_EL(TALN)

      REAL*8   TALN,TAU,TAIL_INTEGRAND,TM(8,6)
      EXTERNAL          TAIL_INTEGRAND

      REAL*8 R_EL,SUM,TEMP
      !REAL*8 R_QE,DR_E_1SIG,DR_Q_1SIG,RCUT(100)
      REAL*8 TAU_PASS

      COMMON/TAIL_INTEGRAL/TAU_PASS,TM

      REAL*8 AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,AMT,TARA,
     &     TARZ,FERMOM,AMM,AMN,CHBAR,BARN,S,X,SX,SXP,Y,YM,W2,ALS,
     &     ALX,ALM,ALY,SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,
     &     TMI

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     &     SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      DATA NENTRY/0/

      NENTRY = NENTRY + 1
      TAU    = DEXP(TALN) - XS

      CALL TAILS(TAU,TM)

!
!     ELASTIC R
!
!$$$  R_EL = ((AMT/AMH)*SX-Y)/(AMT/AMH + TAU)
      R_EL = (SX-Y)/(1.0 + TAU)

      TAU_PASS = TAU   ! PASS TAU ARGUMENT TO TAIL_INTEGRAND
      SUM      = 0.0D0
!
! ADD CONTRIBUTION FROM THE ELASTIC PEAK
!
      TEMP = TAIL_INTEGRAND(R_EL) * 2.0 * AMT/(1.0+TAU)
      SUM  = SUM + TEMP

      TAIL_INT_OVER_R_EL = SUM * (XS + TAU)

      RETURN
      END
!-----------------------------------------------------------------------------------
      REAL*8 FUNCTION TAIL_INTEG_EL(TAU_MIN,TAU_MAX)
      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     &     FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     &     SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/TAIL_INTEGRAL/TAU_PASS,TM

      REAL*8   TAU_MIN,TAU_MAX
      REAL*8   TAIL_INT_OVER_R_EL,TEMP
      EXTERNAL TAIL_INT_OVER_R_EL
      REAL*8   TCUT(100)

      DIMENSION NGAUSSPT(7)

      REAL*8 AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,AMT,TARA,
     &     TARZ,FERMOM,AMM,AMN,CHBAR,BARN,S,X,SX,SXP,Y,YM,W2,ALS,
     &     ALX,ALM,ALY,SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,
     &     TMI


      REAL*8 TAU_PASS,TM(8,6) !,R,TAIL_INTEGRAND

      REAL*8 DGAUSS ! DGQUAD                 ! CSL

      !REAL TEMP1,TEMP2,TEMP3,WEIGHT
      DATA      NGAUSSPT/8,16,16,8,16,16,8/
      !---------------------------------------------------------------

      TCUT(1) = DLOG(TAU_MIN+XS)
      TCUT(2) = DLOG(XS)
      TCUT(3) = DLOG(TAU_MAX+XS)

      TAIL_INTEG_EL = 0.0D0

      DO I = 1,2
         TEMP  =  DGAUSS( TAIL_INT_OVER_R_EL, TCUT(I), TCUT(I+1), 1.0D-6 )
         TAIL_INTEG_EL = TAIL_INTEG_EL + TEMP
      ENDDO

      RETURN
      END
!------------------------------------------------------------------------------------
      SUBROUTINE CONKIN(SNUC,AMTAR,IPOL)
!     SET OF KINEMATICAL CONSTANTS

      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/POL/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS

      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)
!
!     POLARIZATION SPECIFICATION
!
!     IPOL = 0 : PARALLEL, UNPOLARIZED
!     IPOL = 1 : PARALLEL, POLARIZED
!     IPOL = 2 : PERPENDICULAR, UNPOLARIZED
!     IPOL = 3 : PERPENDICULAR, POLARIZED
!

      AMP  = AMTAR              ! CSL AMP is redefined here
      AP   = 2.*AMP
      AMP2 = AMP**2
      AP2  = 2.*AMP**2          ! 2.0*MP**2

      S    = SNUC * AMP/AMH     ! S = S*(M_TARGET/M_HADRON)
      X    = S    * (1.-YS)     ! X = (1-Y)S
      SX   = S - X              ! SX
      SXP  = S + X              ! SP
      YM   = Y + AL2            ! Q^2 + 2M^2 = Q_M^2

      TPL  = S**2 + X**2
      TMI  = S**2 - X**2


!SL   W2   = AMP2 + S - Y - X            ! USE PROTON MASS FOR W2
      W2   = AMT * (AMT + (S-X)/AMH) - Y ! USE HE3 MASS FOR W2
!SL   Makes no difference.



      ALS  = S*S - AL2 * AP2             ! S^2 - 2M^2*M^2 = LAMBDA_S = S^2 IF M=0
      ALX  = X*X - AL2 * AP2
      ALM  = Y*Y  + 4.*AML2*Y
      ALY  = SX**2+ 4.*AMP2*Y            !  LAMBDA_Q

      SQLS = DSQRT(ALS)                  ! SQRT(LAMBDA_S) = S
      SQLX = DSQRT(ALX)
      SQLY = DSQRT(ALY)                  ! SQRT(LAMBDA_Q)
      SQLM = DSQRT(ALM)

      ALLM  = DLOG( (SQLM+Y)/(SQLM-Y) )/SQLM
      AXY   = PI*(S-X)
      AN    = 2. * ALFA**2/SQLS * AXY *BARN * AMH/AMP

!     TAMIN = (SX-SQLY)/AP2
      TAMAX = (SX+SQLY)/AP2               ! TAU_MAX
      TAMIN =  -Y /AMP2 /TAMAX            ! TAU_MIN
      AS    = S /2. /AML /SQLS
      BS    = 0.
      CS    = -AML/SQLS

      IF (IPOL/2.EQ.0) THEN ! PARALLEL CONFIGURATION
        AE  = AMP/SQLS
        BE  = 0.
        CE  = -S/AP/SQLS
      ELSE                  ! PERPENDICULAR CONFIGURATION
        SQN = DSQRT( S*X*Y - ALY*AML2 - AMP2*Y*Y)
        AE  = (-S*X+AP2*YM)/SQLS/SQN/2.
        BE  = SQLS/SQN/2.
        CE  = -(S*Y+AL2*SX)/SQLS/SQN/2.
      ENDIF

      APQ   =           -Y * (AE-BE) + CE*SX  ! Q*ETA
      APN   =  (Y+4.*AML2) * (AE+BE) + CE*SXP ! (K1+K2)*ETA

      DK2KS = AS*YM + AL2*BS + CS*X           ! K_2*KSI
      DKSP1 = AS*S  +  BS*X  + CS*AP2         ! KSI*P

      DAPKS = 2. * (AL2* (AS*AE+BS*BE) + AP2*CS*CE + YM*(AS*BE+BS*AE)
     .           +    S* (AS*CE+CS*AE)             +  X*(BS*CE+CS*BE)) ! KSI*ETA

      RETURN
      END

!DECK  ID>, BORNIN.
!***************** BORNIN *************************************

      SUBROUTINE BORNIN(SIBOR)
!
!     SIBOR IS BORN CROSS SECTION WITH POLARIZED INITIAL
!     LEPTON AND POLARIZED TARGET
!     SIAMM IS CONTRIBUTION OF ANOMALOUS MAGNETIC MOMENT.
!
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/POL/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS

      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)

      COMMON/TAIL/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH

      COMMON/PRINT/IPRI1
      DIMENSION SFM0(8),TM(8),SFM(8)
!
!     FIRST DETERMINE REACTION REGION (ELASTIC, QUASI-ELASTIC OR INELASTIC)
!
      IPRI1 = 1
      CALL STRF(0D0,0D0,SFM0,SFM)

      IPRI1 = 0

      TM(1) =  -(2.*AML2-Y)
      TM(2) = (-(AMP2*Y-S*X)) / (2.*AMP2)
      TM(3) = (2.* (APQ*DK2KS - DAPKS*Y) *AML) /AMP
      TM(4) = APQ/AMP * ( -(DK2KS*SX - 2.*DKSP1*Y) * AML) /AMP2
      TM(7) = (-(4.*AML2 + 3.*APN**2 - 3.*APQ**2+Y))      /2.
      TM(8) = APQ/AMP * (-3.*(APN*SXP- APQ*SX))           /(2.*AP)

      EK    = (3.*APQ**2-Y)/AMP2
      TM(5) = -EK*TM(1)
      TM(6) = -EK*TM(2)
      SSUM  = 0.
      DO 1 ISF=ISF1,ISF2,ISF3
         PPOL = 1.
         IF (ISF.EQ.3 .OR. ISF.EQ.4) PPOL = -PN
         IF (ISF.GE.5)               PPOL = QN/6

         SSUM = SSUM + TM(ISF)*SFM0(ISF)*PPOL
    1 CONTINUE

      SIBOR = SSUM*2.0*AN /Y**2.   ! 2.0*AN/Y**2 = 4*PI*ALPHA**2/LAMBDA_S*(S*SX)/Q^4
      RETURN
      END

!DECK  ID>, DELTAS.
!***************** DELTAS *************************************

      SUBROUTINE DELTAS(DELTA,DELINF,TR,FACTOR1,DEL_SUB)
!
! DELTA IS FACTORIZING PART OF VIRTUAL AND REAL LEPTONIC BREMSSTRAHLUNG
!
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/CHOI/W2_ELAS
      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)

      DEL1  = -YM*(ALM*ALLM**2/2.+2.*FSPEN(2D0*SQLM/(Y+SQLM))-PI2/2.)/SQLM
      DEL2  = (3.*Y/2.+4.*AML2)*ALLM-2.

      SUM   = VACPOL(Y)

      AJ0   = 2.*(YM*ALLM-1.)
      DELTAI= AJ0 * DLOG( DABS(W2-W2_ELAS) /AML /DSQRT(W2) )

      SS    = X + Y
      XX    = S - Y

      ALSS  = SS**2 - 2.*W2*AL2
      ALXX  = XX**2 - 2.*W2*AL2

      SQLSS = DSQRT(ALSS)
      SQLXX = DSQRT(ALXX)

      ALLSS = DLOG( (SQLSS+SS) / (-SQLSS+SS) ) /SQLSS
      ALLXX = DLOG( (SQLXX+XX) / (-SQLXX+XX) ) /SQLXX

      DLM   = DLOG(Y/AML2)
      SFPR  = DLM**2/2. - DLM * DLOG( SS*XX/(AML2*W2) )
     &                  - (DLOG( SS/XX ))**2 /2.
     &                  + FSPEN( ( S*X - Y*AMP2) / ( SS*XX )) - PI2/3.


      DELTA0  = ( SS*ALLSS + XX*ALLXX ) /2. + SFPR

      DELTA   = DELTAI + DELTA0 + DEL1 + DEL2 + SUM

      TR      = ALFA/PI * (DLM-1.)
!SL   TR      = TR + 1.33333*0.0052341

! JSHERE - 032607 - mod b/c it kept stopping here
!      IF (TARZ.EQ.1. .AND. TARA.EQ.1. ) THEN  ! CPROT
       IF(TARZ.EQ.1.D0) THEN
         DELINF  = (DLM-1.) * DLOG( (W2-AMC2   )**2 /(SS*XX) ) ! C$$$
         FACTOR1 = ( (W2-W2_ELAS)**2 / (SS*XX) )**( ALFA/PI*0.5*AJ0 )
!      ELSEIF(TARZ.EQ.2. .AND. TARA.EQ.3.) THEN
       ELSEIF(TARZ.EQ.2.D0) THEN
         DELINF  = (DLM-1.) * DLOG( (W2-W2_ELAS)**2 /(SS*XX) )
         FACTOR1 = ( (W2-W2_ELAS)**2 / (SS*XX) )**TR
      ELSE
         WRITE(*,*) 'JSHERE sub_poltail.f 479'
         WRITE(6,*) TARA,TARZ
         STOP "Problem 3"
      ENDIF



      DEL_SUB = AJ0 * DLOG( DSQRT(SS*XX) /AML /DSQRT(W2) )
     &          + DELTA0 + DEL1 + DEL2 + SUM
      RETURN
      END
!-----------------------------------------------------------------------------------
      REAL*8 FUNCTION VACPOL(T)

      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)

      DIMENSION AM2(3)                   ! AM2 : SQUARED MASSES OF CHARGE LEPTONS
      DATA      AM2/.26110D-6,.111637D-1,3.18301D0/

      SUML = 0.
      DO 10 I=1,3
         A2    = 2.*AM2(I)
         SQLMI = DSQRT( T*T + 2.*A2*T)
         ALLMI =  DLOG( (SQLMI+T) / (SQLMI-T) ) /SQLMI
  10  SUML = SUML + 2.*(T+A2)*ALLMI/3. - 10./9. + 4.*A2*(1.-A2*ALLMI) /3. /T

      IF     (T.LT.1.D0) THEN
        AAA = -1.345D-9
        BBB = -2.302D-3
        CCC =  4.091
      ELSEIF (T.LT.64D0) THEN
        AAA = -1.512D-3
        BBB = -2.822D-3
        CCC =  1.218
      ELSE
        AAA = -1.1344D-3
        BBB = -3.0680D-3
        CCC =  9.9992D-1
      ENDIF

      SUMH  = -( AAA + BBB * LOG(1.+CCC*T) ) *2*PI/ALFA

      VACPOL= SUML+SUMH

      END

!DECK  ID>, FSPENS.
!***************** FSPENS *************************************

      REAL*8 FUNCTION FSPENS(X)
!
!     SPENCE FUNCTION
!
      IMPLICIT REAL*8(A-H,O-Z)

      F   = 0.D0
      A   = 1.D0
      AN  = 0.D0
      TCH = 1.D-16
  1   AN  = AN + 1.D0
      A   = A*X
      B   = A/AN**2
      F   = F+B

      IF (B-TCH) 2,2,1
  2   FSPENS = F

      RETURN
      END
!DECK  ID>, FSPEN.
!***************** FSPEN **************************************

      REAL*8 FUNCTION FSPEN(X)
      IMPLICIT REAL*8(A-H,O-Z)

      DATA F1/1.644934D0/

      IF (X     ) 8,1,1
  1   IF (X-.5D0) 2,2,3

    2 FSPEN = FSPENS(X)
      RETURN

    3 IF (X-1D0 ) 4,4,5
    4 FSPEN = F1 -   DLOG(X) * DLOG( 1D0-X+1D-10  ) - FSPENS(1D0-X)
      RETURN

    5 IF (X-2D0 ) 6,6,7
    6 FSPEN = F1 -.5*DLOG(X) * DLOG( (X-1D0)**2/X ) + FSPENS(1D0-1D0/X)
      RETURN

    7 FSPEN = 2D0*F1 - .5D0  * DLOG(X)**2           - FSPENS(1D0/X)
      RETURN

    8 IF (X+1D0) 10,9,9
   9  FSPEN =        - .5D0  * DLOG( 1D0-X )**2     - FSPENS( X/(X-1D0) )
      RETURN

  10  FSPEN = -.5 * DLOG(1.-X) * DLOG(X**2/(1D0-X)) - F1 + FSPENS(1D0/(1D0-X))

      RETURN
      END


!DECK  ID>, TAILS.
!***************** TAILS **************************************

       SUBROUTINE TAILS(TA,TM)

       IMPLICIT REAL*8(A-H,O-Z)

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/POL/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS

      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)

      COMMON/BSEO/OIS,OIR,OI12,EEIS,EEIR,EEI12,
     . EEI1I2,EB,EEB,TM3(6,4,3)

      DIMENSION TM(8,6),AJM2(2),AJM3(3),II(8)
      DATA II/1,2,3,4,1,2,5,6/

      B2   = ( -ALY*TA + SXP*SX*TA + 2.*SXP*Y) /2. ! B_2
      B1   = ( -ALY*TA - SXP*SX*TA - 2.*SXP*Y) /2. ! B_1

      C1   = -(4.* (AMP2* TA**2 - SX*TA - Y ) * AML2- (S*TA+Y)**2 ) ! C_1
      C2   = -(4.* (AMP2* TA**2 - SX*TA - Y ) * AML2- (TA*X-Y)**2 ) ! C_2

      BB   = 1./SQLY                                 ! F
      SC1  = DSQRT(C1)
      SC2  = DSQRT(C2)

      BI12  = (SXP*(SX*TA+2.*Y))/(SC1*SC2*(SC1+SC2)) ! F_D
      BI1PI2= 1./SC2+1./SC1                          ! F_1+

      BIS   = -B1/SC1/C1+B2/SC2/C2                   ! F_2+
      BIR   =  B2/SC2/C2+B1/SC1/C1                   ! F_2-

      B1I   = -B1/ALY/SQLY                           ! F_I
      B11I  = (3.*B1**2 - ALY*C1)/2. /ALY**2 /SQLY   ! F_II
      SPS   = AS + BS                                ! S_KSI
      SPE   = AE + BE                                ! S_ETA
      CCPE  = (AE-BE) * TA + 2.*CE                   ! R_ETA
      CCPS  = (AS-BS) * TA + 2.*CS                   ! R_KSI

      SIS   = (2.*BI1PI2*SPS    + BIR*SPS*TA + BIS*CCPS)   /2. ! F_{2+}^\KSI
      SIR   = ( (2.*BI12*SPS*TA + BIR*CCPS   + BIS*SPS*TA))/2. ! F_{2-}^\KSI

      SI12  = (BI12*CCPS + BI1PI2*SPS)/2.                      ! F_D^\KSI
      EIS   = (2.*BI1PI2*SPE + BIR*SPE*TA  + BIS*CCPE)   /2.   ! F_{2+}^\ETA
      EIR   = ( (2.*BI12*SPE*TA + BIR*CCPE + BIS*SPE*TA))/2.   ! F_{2-}^\ETA
      EI12  = (BI12*CCPE + BI1PI2*SPE)/2.                      ! F_D^\ETA

      OIS = ((2.*BI1PI2 + BIR*TA)*(CCPE*SPS+CCPS*SPE)+
     .       (CCPE*CCPS + SPE*SPS*TA**2)*BIS         +
     .    8.*BB*SPE*SPS + 4.*BI12*SPE*SPS*TA**2)/4.            ! F_{2+}^{\KSI\ETA}

      OIR = ( ((2.*BI12 + BIS)*(CCPE*SPS      + CCPS*SPE)*TA +
     .      ( CCPE*CCPS + SPE*SPS*TA**2)*BIR  +
     .     4.*BI1PI2*SPE*SPS*TA))/4.                           ! F_{2-}^{\KSI\ETA}

      OI12= ((CCPE*CCPS + SPE*SPS*TA**2)*BI12 +
     .      (CCPE*SPS   + CCPS*SPE)*BI1PI2 + 4.*BB*SPE*SPS)/4. ! F_D^{\KSI\ETA}

      EEIS = ((CCPE**2  + SPE**2*TA**2)*BIS   + 8.*BB*SPE**2 +
     .     4.*BI12*SPE**2*TA**2 + 4.*BI1PI2*CCPE*SPE +
     .     2.*BIR*CCPE*SPE*TA)/4.                              ! F_{1+}^{\ETA\ETA}

      EEIR  = ( ((CCPE**2   + SPE**2*TA**2)*BIR +
     .      4.*BI12*CCPE*SPE*TA + 4.*BI1PI2*SPE**2*TA +
     .      2.*BIS*CCPE*SPE*TA))/4.

      EEI12 = ((CCPE**2 + SPE**2*TA**2)*BI12 + 4.*BB*SPE**2
     .    + 2.*BI1PI2*CCPE*SPE)/4.

      EI1PI2 = (4.*BB*SPE + BI12*SPE*TA**2       + BI1PI2*CCPE)/2.
      EEI1I2 = ( (CCPE**2 + SPE**2*TA**2)*BI1PI2 + 4.*(2.*CCPE-SPE*TA)
     .       * BB*SPE + 8.*B1I*SPE**2 + 2.*BI12*CCPE*SPE*TA**2)/4.

      EB     = ((CCPE - SPE*TA)*BB + 2.*B1I*SPE)/2.

      EEB    = ((CCPE - SPE*TA)**2*BB + 4.*(CCPE - SPE*TA)*B1I*SPE + 4.*B11I
     .       *SPE**2)/4.

      CALL FFU(1,BB ,BIS ,BIR ,BI12 ,BI1PI2,SIR,SIS,SI12,EIS ,EIR ,EI12 ,EI1PI2,TA)
      CALL FFU(2,EB ,EIS ,EIR ,EI12 ,EI1PI2,OIR,OIS,OI12,EEIS,EEIR,EEI12,EEI1I2,TA)
      CALL FFU(3,EEB,EEIS,EEIR,EEI12,EEI1I2,0D0,0D0,0D0 ,0D0 ,0D0 ,0D0  ,0D0   ,TA)

      AJM2(1) =            APQ/AMP
      AJM2(2) =            -1./AMP
      AJM3(1) = (Y-3.*APQ**2) /AMP2
      AJM3(2) =         6.*APQ/AMP2
      AJM3(3) =            -3./AMP2

      DO 15 I = 1,8
        DO 13 L=1,6
   13   TM(I,L)=0
        DO 10 K = 1,I2(I)
          AJK = 1.

          IF (I.EQ.4 .OR. I.EQ.8 ) AJK = AJM2(K)
          IF (I.EQ.5 .OR. I.EQ.6 ) AJK = AJM3(K)

          DO 10 J = K,I1(I)+K-1
            TM(I,J) = TM(I,J) + TM3(II(I),J-K+1,K)* AJK

            IF( (I.EQ.5 .OR. I.EQ.6) .AND. K.EQ.2 )
     .      TM(I,J) = TM(I,J) + TM3(II(I),J-K+1,1)* TA/AMP2
  10    CONTINUE
  15  CONTINUE
      RETURN
      END

!DECK  ID>, FFU.
!***************** FFU ****************************************

       SUBROUTINE FFU(N,BB,BIS,BIR,BI12,BI1PI2,SIR,SIS,SI12
     .        ,EIS,EIR,EI12,EI1PI2,TA)
       IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/POL/AS,BS,CS,AE,BE,CE,APN,APQ,DK2KS,DKSP1,DAPKS

      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)
      COMMON/BSEO/OIS,OIR,OI12,EEIS,EEIR,EEI12,EEI1I2,EB,EEB,TM3(6,4,3)

      HI2  = AML2*BIS - YM*BI12
      SHI2 = AML2*SIS - YM*SI12
      EHI2 = AML2*EIS - YM*EI12
      OHI2 = AML2*OIS - YM*OI12

      GOTO(10,20,30)N

  10  CONTINUE
      TM3(3,1,N) = (8.*(APQ*DK2KS-DAPKS*Y)*AML*HI2)/AMP

      TM3(3,2,N) = (-2.*((2.*(BI12*DK2KS*TA-2.*SHI2)*APQ+(2.*SHI2-
     .       SIR*Y+SIS*YM)*APN+4.*DAPKS*HI2*TA)-4.*((2.*EI12-EIS)*
     .       DK2KS-(SI12-SIS)*APN)*AML2)*AML)/AMP

      TM3(3,3,N) = (2.*(((2.*SI12+SIR-SIS)*APN*TA-2.*DK2KS*EI12*TA
     .             -6.*OHI2-OIR*Y+OIS*YM)-4.*AML2*OI12)*AML)/AMP

      TM3(3,4,N) = (2.*(2.*OI12-OIR+OIS)*AML*TA)/AMP

      TM3(5,1,N) = -2.*(4.*AML2+3.*APN**2-3.*APQ**2+Y)*HI2

      TM3(5,2,N) = -2.*(6.*AML2*APN*EIR-3.*APN**2*BI12*TA+3.*APN*
     .              APQ*BI1PI2+6.*APQ*EHI2+HI2*TA)

      TM3(5,3,N) = -(24.*AML2*EEI12-6.*APN*EI1PI2-6.*APQ*EI12*TA-
     .             2.*BB-BI12*TA**2)

  20  CONTINUE
      TM3(4,1,N) = (-4.*(DK2KS*SX-2.*DKSP1*Y)*AML*HI2)/AMP2

      TM3(4,2,N) = (((2.*(SXP-2.*SX)*SHI2+2.*BI12*DK2KS*SX*TA+8.*
     .  DKSP1*HI2*TA-SIR*SXP*Y+SIS*SXP*YM)-4.*(2.*BI12*DK2KS-BIS*
     .  DK2KS-SI12*SXP+SIS*SXP)*AML2)*AML)/AMP2

      TM3(4,3,N) = ((((SXP*TA-YM)*SIS-(SXP*TA-Y)*SIR+2.*BI12*DK2KS
     .       *TA+6.*SHI2-2.*SI12*SXP*TA)+4.*AML2*SI12)*AML)/AMP2

      TM3(4,4,N) = (-(2.*SI12-SIR+SIS)*AML*TA)/AMP2

      TM3(6,1,N) = (-3.*(APN*SXP-APQ*SX)*HI2)/AMP

      TM3(6,2,N) = (-3.*(2.*(APN*BIR+EIR*SXP)*AML2-(2.*BI12*SXP*TA
     .   -BI1PI2*SX)*APN+(BI1PI2*SXP+2.*HI2)*APQ+2.*EHI2*SX))/(2.*AMP)

      TM3(6,3,N) = (-3.*(8.*AML2*EI12-APN*BI1PI2-APQ*BI12*TA-EI12*
     .    SX*TA-EI1PI2*SXP))/(2.*AMP)

  30  CONTINUE

      TM3(1,1,N) = -4.*(2.*AML2-Y)*HI2

      TM3(1,2,N) = 4.*HI2*TA

      TM3(1,3,N) = -2.*(2.*BB+BI12*TA**2)

      TM3(2,1,N) = (((SXP**2-SX**2)-4.*AMP2*Y)*HI2)/(2.*AMP2)

      TM3(2,2,N) = (2.*AML2*BIR*SXP-4.*AMP2*HI2*TA-BI12*SXP**2*TA+
     .             BI1PI2*SXP*SX+2.*HI2*SX)/(2.*AMP2)

      TM3(2,3,N) = (2.*(2.*BB+BI12*TA**2)*AMP2+4.*AML2*BI12-BI12*
     .             SX*TA-BI1PI2*SXP)/(2.*AMP2)
       RETURN
       END

!-----------------------------------------------------------------------------------
!DECK  ID>, STRF.
      SUBROUTINE STRF(TA,RR,SFM,SFM0)
!
!     THE PROGRAMM CALCULATES DEEP INELASTIC (ITA=1),
!     ELASTIC (ITA=2), QUASIELASTIC (ITA=3) STRUCTURE FUNCTIONS
!     IN KINEMATICAL POINT (TA,RR).
!          RR=SX-TT,
!          TA=(T-Y)/RR,
!    WHERE TT=T+AMF2-AMP2, AMF2 IS INVARINT MASS OF FINAL HADRONS
!
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/SOFTPHOTON/EXTAI2,EXTAI3

      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI

      COMMON/TAIL/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH

      COMMON/PRINT/IPRI1

      DIMENSION SFM(8),SFM0(8) !,FIGI(4,2)

      COMMON/DEBUG/IDEBUG1
!
! INITIALIZE FORM FACTORS
!
      DO I = 1,8
         SFM(I) = 0.0D0
         SFM0(I) = 0.0D0
      ENDDO

      T   = Y  + RR*TA            ! SAME AS Q**2 WHEN THERE IS NO REAL PHOTON EMITTED
      TT  = SX - RR               ! SAME AS 2*M*NU WHEN THERE IS NO REAL PHOTON
      AMF2= TT - T + AMP2         ! SAME AS W**2 WHEN THERE IS NO REAL PHOTON
      AKS = T / TT                ! SAME AS X_BJORKEN WHEN THERE IS NO REAL PHOTON
      ANU = TT/ AP                ! SAME AS NU WHEN THERE IS NO REAL PHOTON

      B1  = 0.D0
      B2  = 0.D0
      B3  = 0.D0
      B4  = 0.D0
!
!     IF T<0, THEN INELASTIC (OR QUASI-ELASTIC) CHANNEL IS NOT POSSIBLE
!
      IF (T.LE.0.0   ) RETURN
      IF (AMF2.LT.0.0) RETURN
!
!     CONTRIBUTION FROM THE ELASTIC CHANNEL
!
      EPSI = AP2/T*(AMT/AMP)**2
      TAU  = T/4./AMT**2
      TAU1 = 1.+TAU

      IF    (TARZ.EQ.1. .AND. TARA.EQ.1.) THEN ! CPROT
        GE = (1.0D0 + T/0.71D0)**(-2)          ! C$$$  PROTON FORM FACTOR USED BY MO AND TSAI
        GM = 2.793*GE                          ! C$$$
      ELSEIF(TARZ.EQ.2. .AND. TARA.EQ.3.) THEN
        CALL FFHE3(T,GE,GM)
      ELSE
        STOP "Problem 4"
      ENDIF


      XNU_EL = 0.5*T/AMT

      IF (ABS(ANU-XNU_EL).GT.1.0D-3) THEN
         SPREAD = 0.0
      ELSE
         SPREAD = 1.0
      ENDIF

      F1 =     AMT*TAU *       GM**2            * SPREAD
      F2 = 2.0*AMT*TAU * (GE**2+TAU*GM**2)/TAU1 * SPREAD
      G1 =     AMT*TAU *    GM*(GE+TAU*GM)/TAU1 * SPREAD
      G2 =     AMT*TAU**2 * GM*(GE-GM    )/TAU1 * SPREAD

      FACTOR  = AMP/AMT
      FACTOR2 = FACTOR**2
      FACTOR3 = FACTOR2*FACTOR

      SFM(1) = SFM(1) + EXTAI2 * FACTOR            * (UN*F1+QN/6.*B1)
      SFM(2) = SFM(2) + EXTAI2 * FACTOR  * EPSI    * (UN*F2+QN/6.*B2)
      SFM(3) = SFM(3) + EXTAI2 * FACTOR2 * EPSI    * (G1+G2)
      SFM(4) = SFM(4) + EXTAI2 * FACTOR3 * EPSI**2 * G2
      SFM(5) = SFM(5) + EXTAI2 * FACTOR3 * EPSI**2 * B1
      SFM(6) = SFM(6) + EXTAI2 * FACTOR3 * EPSI**3 * (B2/3.+B3+B4)
      SFM(7) = SFM(7) + EXTAI2 * FACTOR  * EPSI    * (B2/3.-B3)
      SFM(8) = SFM(8) + EXTAI2 * FACTOR2 * EPSI**2 * (B2/3.-B4)
      RETURN
      END
!-----------------------------------------------------------------------------------
!SL      SUBROUTINE FFHE3(T,GE,GM)
!SL      IMPLICIT REAL*8(A-H,O-Z)
!SL      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
!SL     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
!        COMMON/POLTARG/AMT,TARA,TARZ
!SL      TF=T/CHBAR**2
!SL      IF (TF.GE.0.0) THEN
!SL         QF=SQRT(TF)
!SL      ELSE
!SL         QF=-SQRT(-TF)
!SL      ENDIF
!SL      A=.675
!SL      B=.366
!SL      C=.836
!SL      AM=.654
!SL      BM=.456
!SL      CM=.821
!SL      D=-6.78D-3
!SL      P=.9
!SL      Q0=3.98
!SL      F0=DDEXP(-A**2*TF) - B**2*QF**2*DDEXP(-C**2*TF)
!SL      FM=DDEXP(-AM**2*TF) - BM**2*QF**2*DDEXP(-CM**2*TF)
!SL      DF=D*DDEXP(-((QF-Q0)/P)**2)
!SL      GE=(F0+DF)*TARZ
!SL      GM=FM*TARA      * (-2.13)
!SL      END

      SUBROUTINE FFHE3(T,GE,GM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      QSQP=ABS( T/CHBAR**2 )  ! fm^-2

      FCH  = FORMC(QSQP)
      FMAG = FORMM(QSQP)

      GE=(FCH)*TARZ
      GM= FMAG*TARA      * (-2.13)
      END


!-----------------------------------------------------------------------------------
      REAL*8 FUNCTION TAIL_INTEGRAND(R)
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 SFM0(8),SFM(8),TM(8,6)
      COMMON/CMP/AMP,AMP2,AP,AP2,AML,AML2,AL2,AMC2,AMH,
     .FERMOM,AMM,AMN,CHBAR,BARN,ISF20
      COMMON/POLTARG/AMT,TARA,TARZ
      COMMON/SXY/S,X,SX,SXP,Y,YM,W2,ALS,ALX,ALM,ALY,
     .SQLS,SQLX,SQLY,SQLM,ALLM,AN,TAMIN,TAMAX,XS,YS,TPL,TMI
      COMMON/TAIL/UN,PL,PN,QN,ITA,ISF1,ISF2,ISF3,IRE,ICH
      COMMON/P/PI,PI2,ALFA,I1(8),I2(8)
      COMMON/TAIL_INTEGRAL/TAU,TM

      CALL STRF(TAU,R,SFM,SFM0)
      SUM = 0.0D0
      DO ISF = ISF1,ISF2,ISF3
         IF (ISF.EQ.3.OR.ISF.EQ.4) THEN
            PPOL = -PN
         ELSE IF (ISF.EQ.5) THEN
            PPOL = QN/6.0D0
         ELSE
            PPOL = 1.0D0
         ENDIF
         DO IRR = 1,I1(ISF)+I2(ISF)-1
            IF (IRR.EQ.1) THEN
               TEMP = -0.5*TM(ISF,IRR)*R**(IRR-2)*PPOL*(SFM(ISF)/(Y+
     &              R*TAU)**2-SFM0(ISF)/Y**2)
            ELSE
               TEMP = -0.5*TM(ISF,IRR)*R**(IRR-2)*PPOL*SFM(ISF)/(Y+
     &              R*TAU)**2
            ENDIF
            SUM = SUM + TEMP
         ENDDO
      ENDDO
      TAIL_INTEGRAND = SUM
      RETURN
      END

      REAL*8 FUNCTION DDEXP(X)
      IMPLICIT REAL*8(A-H,O-Z)
        DDEXP=0.
        IF(X.GT.-50.)DDEXP=EXP(X)
      RETURN
      END


      REAL*8 FUNCTION GET_W(E0,TH,XNU)
      IMPLICIT REAL*8(A-H,O-Z)
      XM  = 938.272
      XMT = XM
      PI  = ACOS(-1.)
      THR = TH*PI/180.0
      S2  = ( SIN(THR/2.0) )**2

      EP = E0-XNU
      Q2 = 4.0*E0*EP*S2
      W2 = XM**2 - Q2 + 2.0*XM*XNU
      GET_W  = SQRT(W2)
      RETURN
      END
