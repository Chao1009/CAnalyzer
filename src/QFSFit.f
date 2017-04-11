!=======================================================================
! An interface to get the cross section
! Accepts unit in MeV and radian, output in ub/MeV/sr
      SUBROUTINE GETXS(Z, A, Ei, Ef, theta, xs)
     &bind(C, name = "QFS_xs")
!-----------------------------------------------------------------------
      USE, intrinsic :: ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(IN), value :: Z, A, Ei, Ef, theta
      real(C_DOUBLE), intent(OUT) :: xs
      real*8 sig, nu, E, TH, PF/130.0/, EPS/10.0/, EPSD/0.0/
      real*8 NORM/1D7/ ! convert the cross section to ub/MeV/sr
      REAL*8 PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      DATA PM/938.9187D0/                ! avg. proton/neutron mass in MeV
      DATA DM/1232.D0/                   ! Delta mass in MeV
      DATA ALPHA/7.297352568D-3/          ! 1.0D0/137.03599911D0
      DATA HBARC/197.326968D0/           ! hbar*c in MeV*fermi
      DATA PI/3.1415926535897932D0/      ! 17 digits
      DATA EMASS/0.510998918D0/          ! electron mass in MeV
      DATA UP/2.792847351D0/             ! proton  m.m. in nuclear magnetons
      DATA UN/-1.91304273D0/             ! neutron m.m. in nuclear magnetons
      DATA PIMASS/139.57018D0/           ! charged pion mass in MeV
      LOGICAL Qdep
      REAL*8 XPAR0,XPAR1
      COMMON/QDEPENDENCE/XPAR0,XPAR1,Qdep
      DATA Qdep/.FALSE./
      DATA XPAR0/0.686704D0/
      DATA XPAR1/0.0453828D0/
      REAL*8 SIGQFS, SIG2N, SIGDEL, SIGR1, SIGR2, SIGX

      nu = (Ei - Ef)
      TH = theta/PI*180D0
      sig = 0D0
      sig = sig + SIGQFS(Ei,TH,nu,Z,A,EPS,PF)               ! Q.E. peak
      sig = sig + SIG2N (Ei,TH,nu,Z,A,PF)                   ! Dip region
      sig = sig + SIGDEL(Ei,TH,nu,A,EPSD,PF)                ! Delta
      sig = sig + SIGR1 (Ei,TH,nu,A,PF)                     ! Resonance 1500 MeV
      sig = sig + SIGR2 (Ei,TH,nu,A,PF)                     ! Resonance 1700 MeV
      sig = sig + SIGX  (Ei,TH,nu,A)                        ! DIS


      xs = sig*NORM
      END SUBROUTINE GETXS
!=======================================================================

!-------------------------------------------------------------------------------
! QFS FORTRAN 77 Code
! by: J.W. Lightbody, Jr. & J.S. O'Connell
! October 1987
! based on the paper:
!           "Modeling single arm electron scattering and nucleon production
!            from nuclei by GeV electrons"
! published in: Computers in Physics, May/June 1988, pp. 57-64
!
!  CALCULATE QUASIELASTIC, TWO-NUCLEON, RESONANCE, AND X-SCALING
!  ELECTRON SCATTERING CROSS SECTION FOR ARBITRARY NUCLEI
!
!  WRITTEN BY J.W. LIGHTBODY JR. AND J.S. O'CONNELL
!  NATIONAL BUREAU OF STANDARDS, GAITHERSBURG, MD 20899
!  OCTOBER 1987
!-------------------------------------------------------------------------------
! input  = "input_qfs.txt" see INPUT FILE FORMAT info below
! output = two files into the './output/' directory
!          "[TAG]_sig_qfs.out": includes info from input file and lists in order
!                               nu, W, radiated xs, and born xs
!          "[TAG]_sig_qfs_all.out": first line is energy, angle, and scale.
!                                   then each row lists in order: nu , W ,
!                                   QE xs , dip region, delta xs, W=1500 xs,
!                                   W=1700 xs, dis xs, tot born xs, radiated xs
!-------------------------------------------------------------------------------
!  09/16/02 - K. Slifer, Temple U.
!     Modified subroutine "RADIATE" to include external radiation using
!     formalism from S. Stein et al. Phys. Rev. D 12 1884-1919 (1975)
!     In the resonance region, the new code gives identical results
!     to original code if external radiation is turned off.
!     At larger nu, the new code predicts a much larger tail than the original.
!
!     Included option 'Qdep' to introduce a Q^2 dependent correction to the
!     Q.E. peak, based on a comparison with Saclay and SLAC carbon data.
!     The correction is only to the peak height.  A more sophistacated
!     correction would adjust the gaussian width as well.
!
!     The most recent F.F. are provided as an option, but they are commented
!     out at present. See functions GEP,GEN for example
!
!     If the radiated cross section exhibits glitches adjust the parameters
!     "PREC" or "DEL" in subroutine RADIATE
!-------------------------------------------------------------------------------
! 03/19/07 - downloaded by Jaideep Singh (UVa) from:
!            http://www.jlab.org/~slifer/codes/qfs/nqfs.tar.gz
!         (0) started with "nqfs.f"
!         (1) cleaned up general formating
!         (2) "COMMON"ed and updated physical constants using CODATA 2002
!         (3) programs terminates after one iteration
!         (4) removed "not quite debugged" rates calculation
!         (5) removed stuff specific to E94010, i.e. "TARGETCELL"
!         (6) removed interactive mode; inputs are now read from a file
!         (7) cleaned up if/then/else statements
!         (8) commented out vestigial calcs in the main loop, see "qsfv.f.vax"
!         (9) better descriptions of variables, functions, subroutines
!        (10) converted all constants to REAL*8: "1.0"->"1.0D0" & "1.E2"->"1.D2"
!        (11) SIGDEL has some ambiguity for where the Delta is located
!        (12) I believe that all energy units are MeV
!        (13) labelled equation numbers in QFS paper as (LON),
!                  where "LO" = Lightbody/O'Connell and N = equation #
!             labelled equation numbers in STEIN as (AN),
!                  where "A" = stein appendix A and N = equation #
!        (14) implemented my own recursive spence function "JSPENCE.F" in radcor
!             results in radiated xs that are systematically smaller by 0.5%
!        (15) cleaned up and organized radcor notation, still STEIN formalism
!        (16) user is now forced to input XI, the collistional thickness
!             instead of using old and usually wrong calculation
!        (17) moved many program options into an external control file
!        (18) there is a "kink" near the pion threshold, this is due to the
!             awkward nature of the exponential Delta "turn on" at this xW
!             exponential turn on = 1-exp(-x/x0) , x0 = 5 MeV, x = xW-xW_thr
!        (19) the variable W is really nu because they used to call it omega
!             therefore xW = missing mass
!        (20) works with g77 on a linux machines, but no other guarnatees
!        (21) uses extensions that are not strictly in 77 ANSI standard:
!             http://www.fortran.com/fortran/F77_std/rjcnf0001.html
!        (22) compile with "g77 -ffixed-line-length-none", otherwise errors
!        (23) made some modifications to the radiative corrections, see RADIATE
!-------------------------------------------------------------------------------
! REFERENCES: MOTS69, Radiative Corrections to Elastic and Inelastic ep and
!                     mu-p Scattering by: L.W. Mo and Y.S. Tsai
!                     Reviews of Modern Physics, 41, p205-235, Jan 1969
!          ***http://link.aps.org/abstract/RMP/v41/p205
!
!             TSAI71, Radiative Corrections to Electron Scattering,
!                     by: Yung-Su Tsai, SLAC PUB 848, Jan 1971
!          ***http://www.slac.stanford.edu/pubs/slacpubs/0750/slac-pub-0848.pdf
!
!             STEIN, Electron scattering at 4° with energies of 4.5-20 GeV
!                    by: S. Stein, W. B. Atwood, E. D. Bloom, R. L. A. Cottrell,
!                        H. DeStaebler, C. L. Jordan§, H. G. Piel,
!                        C. Y. Prescott, R. Siemann, and R. E. Taylor
!                    Phys. Rev. D 12, 1884 - 1919 (October 1975)
!          ***http://link.aps.org/abstract/PRD/v12/p1884
!             Miller72, Inelastic Electron-Proton Scattering at Large Momentum
!                       Transfers and the Inelastic Structure Functions
!                       of the Proton, by: G. Miller, E. D. Bloom, G. Buschhorn,
!                       D. H. Coward, H. DeStaebler, J. Drees, C. L. Jordan,
!                       L. W. Mo, and R. E. Taylor, J. I. Friedman,
!                       G. C. Hartmanna, H. W. Kendall, and R. Verdier
!                       Phys. Rev. D 5, 528 - 544 (Feb. 1972)
!          ***http://link.aps.org/abstract/PRD/v5/p528
!-------------------------------------------------------------------------------
! INPUT FILE FORMAT:
!
! TAG
! Z , A
! E , TH
! units , DELTAW , WMIN , WMAX
! PF , EPS , EPSD
! Tb , Ta
! xib , xia
! PREC , DEL
! RAD , extrad , Qdep
!-------------------------------------------------------------------------------
! SAMPLE INPUT FILE:
!
! HE_2135_06
! 2.0 , 3.0
! 2135.0 , 6.0
! nb , 5.0 , 55.0 , 1225.0
! 130.0 , 10.0 , 15.0
! 2.022E-03 , 6.389E-02
! 0.13405E-01 , 0.42356E+00
! 0.00001 , 10.0
! T , T , F
!-------------------------------------------------------------------------------
! DESCRIPTION OF VARIABLES:
!
! TAG    = 10 character string, base of output filenames, must be exactly 10
!          characters long, otherwise bad/weird things happen to good people
! Z      = double, charge of nucleus (for He3, Z = 2)
! A      = double, relative atomic mass (for He3, A = 3)
! E      = double, beam energy in MeV
! TH     = double, scattering angle in degrees
! UNITS  = 2 character string, chooses units for cross section
!          (cm,millibarn,microbarn,nanobarn,picobarn): cm , mb , ub , nb , pb
! DELTAW = double, output energy (nu) bin size
! WMIN   = double,  low end of output energy (nu) range in MeV
! WMAX   = double, high end of output energy (nu) range in MeV
! PF     = double, FERMI MOMENTUM [MEV/C].         (Affects width of QE),
!          typ = 130.0 (MeV)
! EPS    = double, NUCLEON SEPARATION ENERGY [MEV].(Shifts QE central value)
!          typ =  10.0 (MeV)
! EPSD   = double, DELTA SEPARATION ENERGY [MEV].  (Shifts Delta central value)
!          typ =   0.0 (MeV)
! Tb     = double, total ext. radiation length Before scattering, unitless
! Ta     = double, total ext. radiation length After  scattering, unitless
! xib    = double, total "collisional" length Before scattering, MeV
! xia    = double, total "collisional" length After  scattering, MeV
! PREC   = double, relative precision desired for integration
!          this parameter helps tell the integration routine when to stop
!          start with 0.01, if there are "glitches" in the radiated spectrum,
!                           then divide PREC by 10 and re-run program
!          repeat this process of dividing by 10 until "glitches" are gone
! DEL    = double, size of low energy corner, MeV, typ = 5 or 10 MeV
!          should satisfy: DEL > 10*(xia+xib) according to TSAI71
! RAD    = boolean, T if radcor is to be calculated
! extrad = boolean, T if external radcor is to be calculated
! Qdep   = boolean, T if Q.E. peak is to be corrected to agree with world data
!          as of 03/24/07, Karl was not happy with this, so leave F
! note for radiative corrections:   no radcor            , RAD = F
!                                   internal only        , RAD = T & extrad = F
!                                   internal and external, RAD = T & extrad = T
!-------------------------------------------------------------------------------

! GETSCALE Start ---------------------------------------------------------------
! converts from fm^2/sr/MeV to "units"/sr/MeV
! 1 barn = 10^{-24} cm^2 = 100 fm^2
      SUBROUTINE GETSCALE(units,SCALE,UNSCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER units*2

!-----DETERMINE DESIRED UNITS OF CROSS SECTION
      if      (units.EQ.'cm') then ! cm^2/(MeV sr)
          SCALE=1.D-26
          UNSCALE=1.0D0
      else if (units.EQ.'mb') then ! mb/(MeV sr)
          SCALE=1.D+01
          UNSCALE=1.D-27
      else if (units.EQ.'ub')then  ! ub/(MeV sr)
          SCALE=1.D+04
          UNSCALE=1.D-30
      else if (units.EQ.'nb')then  ! nb/(MeV sr)
          SCALE=1.D+07
          UNSCALE=1.D-33
      else if (units.EQ.'pb') then ! pb/(MeV sr)
          SCALE=1.D+10
          UNSCALE=1.D-36
      else
          WRITE(6,*) 'Unknown scale'
          STOP
      endif
!      write(8,'(A,A)') '#',units
      RETURN
      END
! GETSCALE End #################################################################

! FD Start ---------------------------------------------------------------------
! FD = the dipole "form", (LO5), unitless
! QMS and A^2 must have the same units
      REAL*8 FUNCTION FD(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      FD=1.D0/(1.D0+QMS/A**2)**2
      RETURN
      END
! FD End #######################################################################

! FM Start ---------------------------------------------------------------------
! not used? not even in original code? i mean it's just the sqrt of FD, unitless
! QMS and A^2 must have the same units
      REAL*8 FUNCTION FM(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      FM=1.D0/(1.D0+QMS/A**2)
      RETURN
      END
! FM End #######################################################################

! FYUKAWA Start ----------------------------------------------------------------
! no idea what this is? not even used!
! QMS and A^2 must have the same units
      REAL*8 FUNCTION FYUKAWA(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (QMS.LT.1.D-5.OR.A.LT.1.D-5) THEN
          FYUKAWA=0.D0
      ELSE
          ARG=SQRT(QMS/2.D0)/A
          FYUKAWA=ATAN(ARG)/ARG            ! unitless
      ENDIF
      RETURN
      END
! FYUKAWA End ##################################################################

! SIGMOT Start -----------------------------------------------------------------
! SIGMOT = Mott cross section without the elastic recoil factor, fm^2/sr
! see first part of (LO7)
! E in MeV
! THR in radians
      REAL*8 FUNCTION SIGMOT(E,THR)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      SIGMOT=(ALPHA*HBARC*COS(THR/2.D0)/2.D0/E/SIN(THR/2.D0)**2)**2
      RETURN
      END
! SIGMOT End ###################################################################

! RECOIL Start -----------------------------------------------------------------
! RECOIL = recoil factor for elastic scattering, unitless
! see second part of (LO7)
! E and TM must have same units
! THR is radians
      REAL*8 FUNCTION RECOIL(E,THR,TM)
      IMPLICIT REAL*8 (A-H,O-Z)
      RECOIL=1.D0/(1.D0+2.D0*E*SIN(THR/2.D0)**2/TM)
      RETURN
      END
! RECOIL End ###################################################################

! SIGQFS Start -----------------------------------------------------------------
! SIGQFS = quasielastic piece, fm^2/sr/MeV
! incoherent sum of elastic p & n xs w/ gaussian "smearing"
! all in MeV except TH (deg) and Z,A (unitless)
      REAL*8 FUNCTION SIGQFS(E,TH,W,Z,A,EPS,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL Qdep
      COMMON/QDEPENDENCE/XPAR0,XPAR1,Qdep
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      NA   = INT(A)
!---- defining reference kinematics
      GAMR = 120.D0
      PFR  = 230.D0
      QMSRQ= 4.D0*730.D0*(730.D0-115.D0)*SIN(37.1D0*PI/180.D0/2.D0)**2
      QVSRQ= QMSRQ+115.D0**2
!---- parameter AP0 is a function of NA
      AP0  = 840.D0     ! MeV
      AP1  = 750.D0     ! MeV
      IF (NA.EQ.1) THEN
          AP=AP0
      ELSEIF (NA.LT.4) THEN
          AP=AP0+(A-1.D0)*(AP1-AP0)/3.D0
      ELSE
          AP=AP1
      ENDIF
!---- this kinematics
      THR = TH*PI/180.D0
      QMS = 4.D0*E*(E-W)*SIN(THR/2.D0)**2
      QVS = QMS+W**2
!---- START QFS SECTION
      SIGNS  = SIGMOT(E,THR)*RECOIL(E,THR,PM)                          ! (LO7)
!---- proton elastic
      SIGEP = GEP(QMS,AP)**2 + TAU(QMS) * GMP(QMS,AP)**2
      SIGEP = SIGEP/(1.0D0+TAU(QMS) )
      SIGEP = SIGEP+2.0D0*TAU(QMS)*GMP(QMS,AP)**2 * (TAN(THR/2.))**2
      SIGEP = SIGNS*SIGEP                                              ! (LO6)
!---- neutron elastic
      SIGEN = GEN(QMS,AP)**2 + TAU(QMS) * GMN(QMS,AP)**2
      SIGEN = SIGEN/(1.0D0+TAU(QMS) )
      SIGEN = SIGEN+2.0D0*TAU(QMS)*GMN(QMS,AP)**2 * (TAN(THR/2.D0))**2
      SIGEN = SIGNS*SIGEN
!---- final energy of electron
      EPQ    = 4.D0*E**2*SIN(THR/2.D0)**2/2.D0/PM
      EPQ    = EPQ/(1.D0+2.D0*E*SIN(THR/2.D0)**2/PM)+EPS
      EPQ    = E-EPQ

!DEBUG   QSEPQ=4.*E*EPQ*SIN(THR/2.)**2  ! Q^2 at Q.E. Peak
!---- gaussian parameters
      IF (INT(A).EQ.1) THEN
          ARG = (E-W-EPQ)/SQRT(2.D0)/1.D0
          DEN = 2.51D0                                   ! MeV
      ELSE
          GAMQ = GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
          ARG  = (E-W-EPQ)/1.20D0/(GAMQ/2.D0)
          DEN  = 2.13D0*(GAMQ/2.D0)                      ! MeV
      ENDIF
      NQ = INT(ARG)
      IF (ABS(NQ).GT.10) THEN   ! sets the xs to zero if very far from peak
          SIGQ = 0.D0
      ELSE
          SIGQ = (Z*SIGEP+(A-Z)*SIGEN)*EXP(-ARG**2)/DEN
      ENDIF
      SIGQFS = SIGQ

!DEBUG Q2 DEPENDENCE BY KS
      if (Qdep) then
          QQ    = QMS/1.D6
          XCOR  = XPAR1/QQ + XPAR0
          SIGQFS= SIGQFS/XCOR
      endif

      RETURN
      END
! SIGQFS End ###################################################################

! GETQDEP Start ----------------------------------------------------------------
! Introduce Q-dependent correction to Q.E. Peak height
! BY KS
      SUBROUTINE GETQDEP(A,Qdep,XPAR0,XPAR1)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL Qdep


      XPAR1 = 0.0453828D0 ! GeV^2  ! 12/13/02 All available CARBON(only) FIT
      XPAR0 = 0.686704D0  ! GeV^2


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!DEBUG     XPAR1 = 0.0683579 ! GeV^2  ! He3 Fit. Bates data, T.S. Ueng Thesis
!DEBUG     XPAR0 = 0.757269  ! GeV^2
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
! GETQDEP End ##################################################################

! SIG2N Start ------------------------------------------------------------------
! SIG2N = "dip" region between QE and delta, units assumed to be fm^2/sr/MeV
! same form as (LO8), (LO9), (LO10)
! A2 = dipole scale factor
! GAMR = natural FWHM
! GAMQ = FWHM due to fermi motion, not used
! GAMQFR = FWHM due to fermi motion for reference kinematics, not used
! GAM2N = threshold scale parameter in exponential
! all in MeV except TH (deg) and Z,A (unitless)
      REAL*8 FUNCTION SIG2N(E,TH,W,Z,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
!      PM=940.D0           ! js-i am assuming that this is supposed to be 939.D0
!---- defining shape parameters
      A2=550.D0                                         ! MeV
      GAMREF=300.D0                                     ! MeV
      GAMR=GAMREF                                       ! "natural" FWHM
      GAM2N=20.D0                                       ! MeV
!---- defining reference kinematics
      PFR=60.D0
      SIGREF=0.20D-7                                    ! reference xs, units?
      QMSR=4.D0*596.8D0*(596.8D0-380.D0)*SIN(60.D0*PI/180.D0/2.D0)**2
      QVSR=QMSR+380.D0**2
      SIGKIN=0.5D0*SIGMOT(596.8D0,60.D0*PI/180.D0)
      SIGKIN=SIGKIN*(QMSR/2.D0/QVSR+TAN(60.D0*PI/180.D0/2.D0)**2)
      SIGKIN=SIGKIN*QVSR*FD(QMSR,A2)**2
      SIGKIN=SIGKIN*GAMR/GAMREF
      SIGCON=SIGREF/SIGKIN
!---- defining our kinematics
      THR=TH*PI/180.D0
      QMS=4.D0*E*(E-W)*SIN(THR/2.D0)**2
      QVS=QMS+W**2
      EFFMASS=(PM+DM)/2.D0
      SIG=(Z*(A-Z)/A)*SIGMOT(E,THR)                     ! 1st part of (LO8)
      SIG=SIG*(QMS/2.D0/QVS+TAN(THR/2.D0)**2)           ! 2nd part of (LO8)
!---- adding each part of the response function (LO9)
      SIG=SIG*QVS*FD(QMS,A2)**2
      EKAPPA=W-QMS/2.D0/PM
      CMTOT2=PM**2+2.D0*PM*EKAPPA                       ! (LO10)
!      GAMQFR=40.D0               ! FWHM due to fermi motion, not used
!      GAMQF=GAMQFR*(PF/PFR)*(SQRT(QVS)/SQRT(QVSR))
!      GAM=SQRT(GAMR**2+GAMQF**2)
      GAM=GAMR                    ! tot FWHM is just "natural width"
      SIG=SIG*CMTOT2*GAM**2
      SIG=SIG/((CMTOT2-EFFMASS**2)**2+CMTOT2*GAM**2)
      SIG=SIG*(GAMR/GAM)*SIGCON
      SIG2N=SIG
!---- xs is 0 below threshold
      WTHRESH=QMS/4.D0/PM
      IF (W.GT.WTHRESH) THEN
          THRESH=1.D0-EXP(-(W-WTHRESH)/GAM2N)           ! unitless
      ELSE
          THRESH=0.D0                                   ! unitless
      ENDIF
      SIG2N=SIG2N*THRESH                                ! units of SIGREF
      RETURN
      END
! SIG2N End ####################################################################

! SIGDEL Start -----------------------------------------------------------------
! SIGDEL = delta region, units assumed to be fm^2/sr/MeV
! originally used PM = 1219 and EPSD = 15 to give PM+EPSD = 1234 MeV
! changed it such that PM = 1232 and EPSD is chosen accordingly
! same form as (LO8), (LO9), (LO10)
! AD = dipole "scale" parameter in dipole "form" FD
! GAMPI = threshold scale factor in exponential
! DM+EPSD = peak of lorentzian
! GAM = FWHM of lorentzian, quadrature sum of:
!       GAMDP (natural width), GAMQ (fermi motion), GSPRDA (nuclear medium)
! all in MeV except TH (deg) and A (unitless)
      REAL*8 FUNCTION SIGDEL(E,TH,W,A,EPSD,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
!      DM     = 1219.D0   ! js-i am asuuming that this is supposed to be 1232.D0
!---- defining shape parameters
      AD1    = 700.D0                                  ! MeV
      AD0    = 774.D0                                  ! MeV
      GAMDP  = 110.D0                                  ! MeV
      GAMSPRD= 140.D0                                  ! MeV
      GAMR   = 120.D0
      GAMPI  = 5.D0                                    ! MeV
!---- defining reference kinematics
!---- QMSR  and related quantities give the reference for the shape/strength
!---- QMSRQ and related quantities give the reference for the fermi motion
      QFDP   = 1.02D-7                                 ! reference xs, units?
      PFR    = 230.D0
      QMSR   = 4.D0*730.D0*(730.D0-390.D0)*SIN(37.1D0*PI/180.D0/2.D0)**2
      QVSR   = QMSR+390.D0**2
      QMSRQ  = 4.D0*730.D0*(730.D0-115.D0)*SIN(37.1D0*PI/180.D0/2.D0)**2
      QVSRQ  = QMSRQ+115.D0**2
!---- AD and FWHM due to nuclear medium effects depend on A
      NA=INT(A)
      IF (NA.EQ.1) THEN
          QFD=QFDP
          GSPRDA=0.D0
          AD=AD0
      ELSEIF (NA.LT.4) THEN
          QFD=QFDP
          GSPRDA=(A-1.D0)*GAMSPRD/3.D0
          AD=AD0+(A-1.D0)*(AD1-AD0)/3.D0
      ELSE
          AD=AD1
          GSPRDA=GAMSPRD
          QFD=QFDP
      ENDIF
!---- defining our kinematics
      THR = TH*PI/180.D0
      QMS = 4.D0*E*(E-W)*SIN(THR/2.D0)**2
      QVS = QMS+W**2
      EKAPPA = W-QMS/2.D0/PM
      CMTOT2 = PM**2+2.D0*PM*EKAPPA
!  BEGIN DELTA CALCULATION
      IF (NA.GT.1) THEN                              ! FWHM due to fermi motion
          GAMQ=GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
      ELSE
          GAMQ=0.D0
      ENDIF
!---- calc kinematics for peak
      EPD = E-(DM-PM)*(DM+PM)/2.D0/PM
      EPD = EPD/(1.D0+2.D0*E*SIN(THR/2.D0)**2/PM)
      EPD = EPD-EPSD
      WD  = E-EPD
      QMSPK = 4.D0*E*EPD*SIN(THR/2.D0)**2
      QVSPK = QMSPK+WD**2
!
! NOTE WIDTH INCLUDES E-DEPENDENCE,FERMI BROADENING,& SPREADING
!
      WTHRESH = 4.D0*E**2*SIN(THR/2.D0)**2+PIMASS**2+2.D0*PIMASS*PM
      WTHRESH = WTHRESH/2.D0/PM
      THRESHD = 1.D0+PF/PM+PF**2/2.D0/PM**2+2.D0*E*SIN(THR/2.D0)**2/PM
      WTHRESH = WTHRESH/THRESHD
!---- xs is 0 below pion threshold
      IF (W.GT.WTHRESH) THEN
          THRESH=1.D0-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
          THRESH=0.D0
      ENDIF
!---- (032607- js) this form of the behaviour
!---- near the pion threshold gives a wierd kink in the cross sections
!---- this same form of the threshold beahviour is not a problem for the
!---- "dip" region of the 1500 and 1700 resonances
      SIGD = THRESH
!      SIGD = 1.D0
      GAMD = GAMDP                            ! natural FWHM
      GAM  = SQRT(GAMD**2+GAMQ**2+GSPRDA**2)  ! tot FWHM
!---- cross section is formed by taking ratios with reference quantities
      SIGD = SIGD*QFDP*(GAMDP/GAM)
      SIGD = SIGD*CMTOT2*GAM**2
      SIGD = SIGD/((CMTOT2-(DM+EPSD)**2)**2+CMTOT2*GAM**2)
      SIGD = SIGD*FD(QMS,AD)**2/FD(QMSR,AD)**2
      SIGD = SIGD*QVS/QVSR
      SIGD = SIGD*(QMS/2.D0/QVS+TAN(THR/2.D0)**2)
      SIGD = SIGD/(QMSR/2.D0/QVSR+TAN(37.1D0*PI/180.D0/2.D0)**2)
      SIGD = SIGD*SIGMOT(E,THR)/SIGMOT(730.D0,37.1D0*PI/180.D0)
      SIGD = SIGD*A
      SIGDEL = SIGD                                        ! units of QFDP

!      DM = 1232.0D0 ! js-shifting back to 1232.D0 for all other subrouts/funcs
      RETURN
      END
! SIGDEL End ###################################################################

! SIGR1 Start ------------------------------------------------------------------
! SIGR1 = resonance at 1500 MeV, units assumed to be fm^2/sr/MeV
! same form as Delta: (LO8), (LO9), (LO10)
! all in MeV except TH (deg) and A (unitless)
      REAL*8 FUNCTION SIGR1(E,TH,W,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      RM=1500.D0                                      ! MeV
!---- defining shape parameters
      AR0=1000.D0                                     ! MeV
      AR1=1000.D0                                     ! MeV
      GAMR=110.D0                                     ! MeV
      GAMSPRD=140.D0                                  ! MeV
      GAMQFR=120.D0
      GAMPI=5.D0                                      ! MeV
!---- defining reference kinematics
!---- QMSRR  and related quantities give the reference for the shape/strength
!---- QMSQFR and related quantities give the reference for the fermi motion
      PFR=230.D0
      EPSR=0.D0
      QFRP=1.20D-7                                    ! reference xs, units?
      QMSQFR=4.D0*730.D0*(730.D0-115.D0)*SIN(37.1D0*PI/180.D0/2.D0)**2
      QVSQFR=QMSQFR+115.D0**2
      QMSRR=4.D0*10000.D0*(10000.D0-1240.D0)*SIN(6.D0*PI/180.D0/2.D0)**2
      QVSRR=QMSRR+1240.D0**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2.D0/QVSRR+TAN(6.D0*PI/180.D0/2.D0)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.D0*PI/180.D0)
!---- AR and FWHM due to nuclear medium effects depend on A
      NA=INT(A)
      IF (NA.EQ.1) THEN
          QFR=QFRP
          GSPRDA=0.D0
          AR=AR0
      ELSEIF (NA.LT.4) THEN
          QFR=QFRP
          GSPRDA=(A-1.D0)*GAMSPRD/3.D0
          AR=AR0+(A-1.D0)*(AR1-AR0)/3.D0
      ELSE
          AR=AR1
          GSPRDA=GAMSPRD
          QFR=QFRP
      ENDIF
! defining our kinematics
      THR=TH*PI/180.D0
      QMS=4.D0*E*(E-W)*SIN(THR/2.D0)**2
      QVS=QMS+W**2
      IF (NA.GT.1) THEN                              ! FWHM due to fermi motion
          GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
          GAMQ=0.D0
      ENDIF
      CMTOT2=PM**2+2.D0*PM*W-QMS
!---- xs is 0 if below threshold
      WTHRESH=4.D0*E**2*SIN(THR/2.D0)**2+PIMASS**2+2.D0*PIMASS*PM
      WTHRESH=WTHRESH/2.D0/PM
      THRESHD=1.D0+PF/PM+PF**2/2.D0/PM**2+2.D0*E*SIN(THR/2.D0)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF (W.GT.WTHRESH) THEN
          THRESH=1.D0-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
          THRESH=0.D0
      ENDIF
!---- calculating kinematics at peak
      EPR=E-(RM-PM)*(RM+PM)/2.D0/PM
      EPR=EPR/(1.D0+2.D0*E*SIN(THR/2.D0)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)            ! tot FWHM
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2.D0/QVS+TAN(THR/2.D0)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR1=A*THRESH*SIGR
      SIGR1=SIGR1                                     ! units of QFRP
      RETURN
      END
! SIGR1 End ####################################################################

! SIGR2 Start ------------------------------------------------------------------
! SIGR1 = resonance at 1500 MeV, units assumed to be fm^2/sr/MeV
! same form as Delta: (LO8), (LO9), (LO10)
! identical to SIGR1 except for:
!                                RM   =    1500 <--> 1700     (resonance mass)
!                                AR   =    1000 <--> 1200     (dipole scale)
!                                QFRP = 1.20D-7 <--> 0.68D-7  (ref. xs, units?)
!                                kinematics for QMSRR
! all in MeV except TH (deg) and A (unitless)
      REAL*8 FUNCTION SIGR2(E,TH,W,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      RM=1700.D0
!---- defining shape parameters
      AR0=1200.D0
      AR1=1200.D0
      GAMR=110.D0
      GAMSPRD=140.D0
      GAMQFR=120.D0
      GAMPI=5.D0
!---- defining reference kinematics
!---- QMSRR  and related quantities give the reference for the shape/strength
!---- QMSQFR and related quantities give the reference for the fermi motion
      PFR=230.D0
      EPSR=0.D0
      QFRP=0.68D-7                                    ! reference xs, units?
      QMSQFR=4.D0*730.D0*(730.D0-115.D0)*SIN(37.1D0*PI/180.D0/2.D0)**2
      QVSQFR=QMSQFR+115.D0**2
      QMSRR=4.D0*10000.D0*(10000.D0-1520.D0)*SIN(6.D0*PI/180.D0/2.D0)**2
      QVSRR=QMSRR+1520.D0**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2.D0/QVSRR+TAN(6.D0*PI/180.D0/2.D0)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.D0*PI/180.D0)
! AR and FWHM dues to nuclear medium effects depend on A
      NA=INT(A)
      IF (NA.EQ.1) THEN
          QFR=QFRP
          GSPRDA=0.D0
          AR=AR0
      ELSEIF (NA.LT.4) THEN
          QFR=QFRP
          GSPRDA=(A-1.D0)*GAMSPRD/3.D0
          AR=AR0+(A-1.D0)*(AR1-AR0)/3.D0
      ELSE
          AR=AR1
          GSPRDA=GAMSPRD
          QFR=QFRP
      ENDIF
!---- defining our kinematics
      THR=TH*PI/180.D0
      QMS=4.D0*E*(E-W)*SIN(THR/2.D0)**2
      QVS=QMS+W**2
      IF (NA.GT.1) THEN                               ! FWHM due to fermi motion
          GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
          GAMQ=0.D0
      ENDIF
      CMTOT2=PM**2+2.D0*PM*W-QMS
!---- xs is 0 if below threshold
      WTHRESH=4.D0*E**2*SIN(THR/2.D0)**2+PIMASS**2+2.D0*PIMASS*PM
      WTHRESH=WTHRESH/2.D0/PM
      THRESHD=1.D0+PF/PM+PF**2/2.D0/PM**2+2.D0*E*SIN(THR/2.D0)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF (W.GT.WTHRESH) THEN
          THRESH=1.D0-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
          THRESH=0.D0
      ENDIF
!---- calculating kinematics at peak
      EPR=E-(RM-PM)*(RM+PM)/2.D0/PM
      EPR=EPR/(1.D0+2.D0*E*SIN(THR/2.D0)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)             ! tot FWHM
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2.D0/QVS+TAN(THR/2.D0)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR2=A*THRESH*SIGR
      SIGR2=SIGR2                                     ! units of QFRP
      RETURN
      END
! SIGR2 End ####################################################################

! SIGX Start -------------------------------------------------------------------
! SIGX = cross section in DIS, (LO11), fm^2/sr/MeV
! E and W must be MeV
! TH must be degrees, A is unitless
      REAL*8 FUNCTION SIGX(E,TH,W,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
!---- defining real photon parameters
!---- in the paper, these are given by 100 ub and 54x10^3 ub*MeV
!---- shouldn't they be per steradian as well? we'll assume this is the case...
!     SIG0=111.*1.D-30
      SIG0=100.D-4                                         ! fm^2/sr??
!     SIG1=60.*1.D-27
      SIG1=54.D-1                                          ! fm^2*MeV/sr??
!     GAM0=550.D0
      GAM0=650.D0                                          ! MeV
!     R=0.10D0
!      AQ=250.D0 ! not used even in original code!
!---- defining our kinematics
      THR=TH*PI/180.D0
      IF (W.LT.1.D-5) GO TO 4
      QMS=4.D0*E*(E-W)*SIN(THR/2.D0)**2
!---- shape of the real photon cross section
!---- parameterization in this code is more sophisticated than (LO15)
      ARG0=W-QMS/2.D0/PM-PIMASS-PIMASS**2/2.D0/PM          ! w - w_pi (LO16)
      ARG1=ARG0/GAM0
      ARG=ARG1**2/2.D0
      IF (ARG1.GT.8.D0) THEN
          SHAPE=1.D0+SIG1/SIG0/ARG0
      ELSEIF (ARG1.LT.1.D-5) THEN
          SHAPE=0.D0
      ELSEIF (ARG1.LT.0.1D0) THEN
          SHAPE=SIG1*ARG0/2.D0/GAM0**2/SIG0
      ELSE
          SHAPE=(1.D0-EXP(-ARG))*(1.D0+SIG1/SIG0/ARG0)
      ENDIF
      SIGGAM=SIG0*SHAPE                                    ! fm^2
!---- calculating non-real photon xs parts of electron xs
      EKAPPA=W-QMS/2.D0/PM                                 ! (LO12), MeV
      QS=QMS+W**2
      EPS=1.D0/(1.D0+2.D0*QS*TAN(THR/2.D0)**2/QMS)         ! (LO13)
      FLUX=ALPHA*EKAPPA*(E-W)/2.D0/PI**2/QMS/E/(1.D0-EPS)   ! (LO12), 1/MeV
      IF (FLUX.LT.1.D-20) FLUX=0.D0
      SIGEE=FLUX*SIGGAM*FPHENOM(QMS)**2                    ! fm^2/MeV
!     SIGEE=FLUX*SIGGAM
      R=0.56D0*1.D6/(QMS+PM**2)                            ! (LO14), unitless
      FACTOR1=1.D0+EPS*R
      SIGEE=SIGEE*FACTOR1                                  ! (LO11)
 4    SIGX=A*SIGEE                                         ! fm^2/MeV
!---- added to insure there is no funny stuff going on
      IF (W.LE.0.D0) THEN
           SIGX = 0.D0
      ENDIF
      RETURN
      END
! SIGX End #####################################################################

! FPHENOM Start ----------------------------------------------------------------
! FPHENOM = F_x^2, (LO17)
! phenomenological ``form factor'' for virtual photon scattering, unitless
! QMS must be MeV^2
! in the paper, the constants are labelled a1 to a7
      REAL*8 FUNCTION FPHENOM(QMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      A1=.55D0                                   ! a1
      A2=20.D0/1.D6                              ! a2, 1/MeV^2
      B1=.45D0                                   ! a3
      B2=.45D0/1.D6                              ! a4, 1/MeV^2
      C1=0.03D0                                  ! a5
      C2=0.2D0/1.D12                             ! a6, 1/MeV^4
      C3=4.5D6                                   ! a7, MeV^2
      FPHENOM=A1*EXP(-A2*QMS)+B1*EXP(-B2*QMS)
      FPHENOM=FPHENOM+C1*EXP(-C2*(QMS-C3)**2)
      FPHENOM=SQRT(FPHENOM)
      RETURN
      END
! FPHENOM End ##################################################################

! TAU Start --------------------------------------------------------------------
! TAU = kinematic factor found in Rosenbluth xs, unitless
! QMS must be in MeV^2
! proton and neutron mass is assumed to be equal their average
! this causes an "error" of 0.1% for tau
      REAL*8 FUNCTION TAU(QMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      TAU  = QMS/4.0D0/PM**2
      RETURN
      END
! TAU End ######################################################################

!------------------------------ Elastic Form Factors ---------------------------
!     following by KS
!     These are the most recent form factors available.
!     Do not use option 'Qdep' with these FF.
!      REAL*8 FUNCTION GEP(QMS,AP)
!      IMPLICIT REAL*8 (A-H,O-Z)
!      GEP=1./(1.+QMS/AP**2)**2
!      RETURN
!      END
!
!      REAL*8 FUNCTION GEN(QMS,AP)
!      IMPLICIT REAL*8 (A-H,O-Z)
!C     H.ZHU et al. PRL 87, Number 8, 081801-1 (2001)
!
!      PM   = 939.
!      UN   = -1.91304184
!      xp   = 4.36   ! +- 1.11
!      xa   = 0.895  ! +- 0.039
!
!      GEN = -UN * xa * TAU(QMS)/( 1.0+xp*TAU(QMS) )
!      GEN = GEN * GEP(QMS,AP)
!      RETURN
!      END
!
!      REAL*8 FUNCTION GMP(QMS,AP)
!      IMPLICIT REAL*8 (A-H,O-Z)
!C     Gayou et al. PRL 88,number 9, 092301-1 (2002)
!      UP   =  2.7928456
!      GMP  =  UP * GEP(QMS,AP)
!      GMP  =  GMP/(1.0 - 0.13*(QMS-0.04) )
!      RETURN
!      END
!
!      REAL*8 FUNCTION GMN(QMS,AP)
!      IMPLICIT REAL*8 (A-H,O-Z)
!      UN   = -1.91304184
!      GMN  = UN * GEP(QMS,AP)
!      RETURN
!      END
!

! js-The following are the "classical" form factors as listed in (LO5)

! GEP Start --------------------------------------------------------------------
! proton electric form factor, normarlized to proton charge +1
! QMS and AP^2 must have the same units
! simple dipole form
      REAL*8 FUNCTION GEP(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      GEP=1.D0/(1.D0+QMS/AP**2)**2
      RETURN
      END
! GEP End ######################################################################

! GEN Start --------------------------------------------------------------------
! neutron electric form factor
! modified dipole form
      REAL*8 FUNCTION GEN(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      GEN = -UN
      GEN = GEN * TAU(QMS)/( 1.0D0+5.6D0*TAU(QMS) )
      GEN = GEN * GEP(QMS,AP)
      RETURN
      END
! GEN End ######################################################################

! GMP Start --------------------------------------------------------------------
! proton magnetic form factor, normalized to proton m.m. in nuclear magnetons
! simple dipole form
      REAL*8 FUNCTION GMP(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      GMP  =  UP * GEP(QMS,AP)
      RETURN
      END
! GMP End ######################################################################

! GMN Start --------------------------------------------------------------------
! neutron magnetic form factor, normalized to neutron m.m. in nuclear magnetons
! simple dipole form
      REAL*8 FUNCTION GMN(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/QFS_CONST/PM,DM,ALPHA,HBARC,PI,EMASS,UP,UN,PIMASS
      GMN  = UN * GEP(QMS,AP)
      RETURN
      END
! GMN End ######################################################################

! ROM Start --------------------------------------------------------------------
! ROMBERG METHOD OF INTEGRATION - see Numerical Recipes by Press et al
! (0) performs integration using trapezoid rule by iteration
! (1) initial integration bin size is the full range
! (2) integration bin size is halved after each subsequent iteration
! (3) after each iteration, the "area vs. step size" curve is extrapolated
!     to zero step size using Neville's Algorithm for Polynomial Int./Ext.
! (4) loop is terminated when desired presicion is achieved
! A     = lower limit of integration
! B     = upper limit of integration
! EPS   = desired relative precision (must be positive!)
! ANS   = the value of the integral
! K     = number of iterations
! IFUNC = integrand label, 1 for integral over dEs' ; 2 for integral over dEp'
      SUBROUTINE ROM(A,B,EPS,ANS,K,IFUNC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(50,50)
      H=B-A
      K=0
      CALL VALY(A,FA,IFUNC)
      CALL VALY(B,FB,IFUNC)
      W(1,1)=(FA+FB)*H/2.D0
    4 K=K+1
      IF (K.GE.49) GO TO 5
      H=H/2.D0
      SIG=0.D0
      M=2**(K-1)
      DO 1 J=1,M
      J1=2*J-1
      X=A+FLOAT(J1)*H
      CALL VALY(X,F,IFUNC)
!      write(6,*) "DEBUG: ",k,IFUNC,f
    1 SIG=SIG+F
      W(K+1,1)=W(K,1)/2.D0+H*SIG
      DO 2 L=1,K
      IU=K+1-L
      IV=L+1
    2 W(IU,IV)=(4.D0**(IV-1)*W(IU+1,IV-1)-W(IU,IV-1)) ! js-originally one line
     &         /(4.D0**(IV-1)-1.D0)                   ! had to split it
      E=(W(IU,IV)-W(IU,IV-1))/W(IU,IV)
      IF (ABS(E)-EPS) 3,3,4
    3 ANS=W(1,IV)
      RETURN
    5 PRINT 100
  100 FORMAT(' K OVERFLOW')
      CALL EXIT(0)
      END
! ROM End ######################################################################
