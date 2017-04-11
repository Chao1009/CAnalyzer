! *****output is ub/sr-GeV********
! This code is to provide an interface to C/C++ code for generating elastic 3He
! tail with radiation effect
! Target and some of the boolean values are fixed for this specific task

!   a little history;
!      the program was changed to read in the radiation lenghts of the
!      pre window, pre target, post window, post target. The change in the code
!      should have been correct based on previous code.  ie before (find cttt)
!            ttb= (.5 *target thickness / target angle ) + pre window
!      so ttb was assumed to be the total target lenght + window before
!      however examine extern  and the rad lenght is the sum of ttb + twb
!      so the twb in no longer added to pre target to get tta.
!      now the program aggrees with sylvie and tta and ttb are just the target
!      thickness values.
!
!        other features added in the attemp to understand.  the form factors
!        for the proton are now calculated by donal nform.  the subr nform
!        inserted in newtail and a common section included to pass some
!        arguments.
!        fbar is now used for the proton calculation
! 12/07/06 - Jaideep Singh
!
! notes and modifications:
!
!     Implementation of collisional loss, i assume that the incident
!     energy and nu range given by the user refer to the beam energy and
!     spectrometer setting.  therfore i use these values and user supplied
!     coll. energy losses before and after scattering to find the "true"
!     kinematics of the elastic scattering that generates the tail.  then
!     the elastic tail is mapped back to the nu corresponding to the
!     spectrometer setting.  note that the beam energy is the energy before
!     coll. loss and the spectrometer momentum is the final energy after
!     both coll. losses have occurred.  in this way, interpolation is not
!     needed which is good b/c i couldn't get the old version to make sense.
!
!
!     ORIGINALLY SLDTAIL.FOR (IN OFLNE5.RTAIL) MODIFIED BY KLG
!     AND RMS, THICKNESSES BEFORE, TARGET AND AFTER NOW IN RADIATION
!     LENGTHS, READ IN FROM LOGICAL UNIT 5   MARCH 15,1988
!     SUBROUTINE TARGETS REMOVED - PRINTS OUT THICKNESSS
!     REMOVED MAIN PROGRAM COMPUTATION OF THICKNESSES (COMMENTED OUT
!     BY CTTT)
!     ------------------------------------------------------------------
!     THIS PROGRAM COMPUTES RADIATIVE TAIL CORRECTIONS TO ELECTRON
!     SCATTERING CROSS SECTIONS.
!
!     ------------------------------------------------------------------
!     |                                                                |
!     |      ***  DEVELOPED IN UNIU. OF VIRGINIA  ***                  |
!     |                                                                |
!     | 1980 >>> : WRITTEN BY ROSEMARY ALTULMUS AND JOHN WISE          :
!     | 1980 SEP : REVISED BY MIKE JOHNSON                             |
!     | 1983 MAR : REVISED BY T. S. UENG, CHAIRMAN FOR THE REVOLUTION  |
!     | 1986 JUN : DBD
!     |                                                                |
!     ------------------------------------------------------------------
!     -----------------------------------------------------------------
!     THIS UPDATE VERSION (MARCH 1983) IS A LITTLE DIFFERENT FROM THE
!     ORIGINAL VERSION. ALL THE PROGRAMS HAVE BEEN CONVERTED INTO DOU-
!     BLE PRECESSION IN ORDER TO WORK IN THE VAX/11 OR OTHER MEDIUM
!     SIZE COMPUTER. THE INPUTS FORMAT HAVE BEEN CHANGED. THE ORDER OF
!     THE INPUT DATA IS:
!
!           1. EPS: CONVERGENCE FOR SUBROUTINE SIMPSN, ( <0.00001 )
!           2. ITARGT, AT, AWB, AWA, ZT, ZWB, ZWA
!           3. TTPRE, TTPOST, TWB, TWA
!           4. (READ IN FORM FACTOR IF NECESSARY)
!           5. -1, 0., 0., 0. (IF FORM FACTOR IS READED, USED THIS CARD
!                              TO TERMINATE IT)
!           6. (KEYWORDS)
!           7. ANGLE, EINC, TANGLE
!           8. EFIN (OR EMAX,EMIN,EDEL)
!           9. ANY NEGATIVE NUMBER TO TERMINATE IT (OPTION)
!          10. ........ REPEAT FROM 7.
!              (END)
!
!     -----------------------------------------------------------------
!     THE TOTLE ROUTINE NEED FOR GASTAIL ARE:
!            1. STAILS
!            2. EXTERNL
!            3. FBAR
!            4. TARGET
!            5. GTARGET
!            6. XSECT
!            7. XSECTP
!            8. FMFAC
!            9. RADLENG
!           10. SIGBAR
!           11. SPENCE
!           12. SIMPSN
!           13. TERP
!
!     -----------------------------------------------------------------
!     THE POSSIBLE TAIL CORRECTIONS (ALL OPTIONAL) ARE:
!     -----------------------------------------------------------------
!     1. THE EXTERNAL BREMSSTRAHLUNG (I.E., TARGET) CORRECTION.
!        THIS ACCOUNTS FOR THE FINITE TARGET THICKNESS AND THE CONTAIN-
!        ER (i.e., THE WINDOWS). (STEIN, et. al. PHYS REV D, 12, 1884
!        (1975).)
!     2. THE INTERNAL BREMSSTRAHLUNG CORRECTION.
!        THIS ACCOUNTS FOR BREMSSTRAHLUNG AT THE TARGET NUCLEUS ITSELF.
!        THIS COMPUTATION CAN BE CARRIED OUT ONE OF TWO WAYS:
!        A. USING THE EXACT CALCULATION FROM MO AND TSAI, REV MOD PHYS,
!           41, 205 (1969), EQ. B.5.  MO AND TSAI'S VERSION IS HEREIN
!           MODIFIED SLIGHTLY EXCEPT WHEN THE BUILT-IN PROTON FORM FAC-
!           TORS ARE USED (SEE BELOW), THE FORM FACTORS ARE MULTIPLIED
!           BY A FACTOR FBAR (SEE STEIN, EQ.A42 TO A44, OR TSAI,
!           SLAC-PUB-848, JAN. 1971).
!        B. USING THE PEAKING APPROXIMATION FROM STEIN, EQ. A56.  THIS
!           AND THE EXACT CALCULATION CANNOT BE CARRIED OUT DURING THE
!           SAME RUN.
!
!     OTHER PROGRAM OPTIONS INCLUDE:
!
!     3. THE MULTIPLE-PHOTON CORRECTION TERM (FROM STEIN) CAN BE INCLUD-
!        ED IN THE ABOVE TAILS OR NOT.
!     4. THE FORM FACTORS CAN BE READ IN FROM DATA FILES, OR THE PROGRAM
!        CAN USE ONE OF THE BUILT-IN ANALYTICAL EXPRESSIONS.  CURRENTLY,
!        THESE ARE INCLUDED FOR THE PROTON (TO SEE IF THE PROGRAM AGREES
!        WITH M.T.), D2, HE3 AND HE4.
!     5. ANY NUMBER OF INCIDENT ELECTRON ENERGIES AND ANGLES CAN BE
!        USED CORRESPONDING TO A GIVEN PAIR, UP TO 400 FINAL ELECTRON
!        ENERGIES ARE ALLOWED, OR YOU CAN CHANGE THE DIMENSION TO IN-
!        CREASE THE NUMBER.
!
!     -----------------------------------------------------------------
!     CONVENTIONS:
!     -----------------------------------------------------------------
!     1. THE PROGRAM GENERALLY USES UNITS OF (MEV) AND (CM). EXCEPTIONS
!        (IN FMFAC) ARE NOTED.
!     2. OPTIONS ARE CHOSEN AND DATA ENTERED USING KEYWORDS (EXPLAINED
!        BELLOW). THESE CAN BE IN ANY ORDER, EXCEPT THAT KEYWORD 5 (RUN)
!        MUST FOLLOW ALL OTHERS.
!     3. THIS VERSION OF THE PROGRAM USES THE SOLID TARGET GEOMETRY.
!
!     -----------------------------------------------------------------
!     MEANING OF THE INPUT VARIABLES:
!     -----------------------------------------------------------------
!     ITARGT CHOOSES THE TARGET MATERIAL (BY SELECTING THE FORM FACTOR
!     TO BE USED).
!     -1 = PROGRAM WILL READ IN ELASTIC AND MAGNETIC FORM FACTORS FROM
!          CARDS.
!      2 = PROTON.  PROGRAM USES DIPOLE FORM FACTORS FROM MO AND TSAI.
!      3 = HELIUM 3. PROGRAM USES FORM FACTORS FROM MCCARTHY, SICK AND
!                    WHITNEY, PHYS REV C, 15, 1396 (1977).
!      4 = HELIUM 4. PROGRAM USES FORM FACTORS FROM MCCARTHY, SICK AND
!          WHITNEY.
!      5 = DEUTERIUM.  PROGRAM USES FORM FACTORS FROM STEIN.
!     ------------------------------------------------------------------
! CTARG  = 4 character string, target material: Prot , He-3 , He-4 , Deut ,
!          Nitr , Carb , Bery , Alum , Copp , Gold
!     EPS      : THE RELATIVE CONVERGENCE FACTOR FOR FUNCTION SIMPSN
!     AT       : THE ATOMIC WEIGHT OF THE SCATTERING NUCLEUS.
!     ZT       : Z FOR THE SCATTERING NUCLEUS.
!     ANGLE    : SCATTERING ANGLE (DEGREE)
!     EINC_0   : Double precision, beam energy in MeV *before* coll. loss
!     EINC     : INCIDENT ENERGY OF THE ELECTRON (MEV)
!     TANGLE   : TARGET ANGLE (DEGREE)
!     EFIN     : FINAL ENERGY OF THE ELECTRON (MEV)
!     EFMAX    : CHOOSED MAXIMUM FINAL ENERGY (MEV)
!     EFMIN    : CHOOSED MINIMUM FINAL ENERGY (MEV)
!     EDEL     : CHOOSED FINAL ENERGY SEPARATION (MEV)
!     -----------------------------------------------------------------
!     THE SOLID TARGET USED HERE IS THE REGULAR THIN PLATE TARGET.
!     ----------------------------------------------------------------
!
!     MEANING OF THE KEYWORDS:
!     -----------------------------------------------------------------
!     INTERN = .TRUE. IF THE INTERNAL BREMSSTRAHLUNG IS TO BE CALCULAT-
!                     ED EXACTLY (FROM MO AND TSAI).
!     PEAK   = .TRUE. IF THE PEAKING APPROXIMATION (FROM STEIN) IS USED
!                     FOR  THE INTERNAL TAIL.
!     EXTERN = .TRUE. IF THE EXTERNAL TAIL IS TO BE CALCULATED (STEIN).
!     MULTIP = .TRUE. IF THE MULTIPLE-PHOTON CORRECTION IS TO BE
!                     INCLUDED IN THE TOTAL TAIL CORRECTIONS.
!     AMROUN_OR_MSW = .TRUE. IF the AMROUN 3He paramterization should be
!                            used instead of McCarthy, Sick, and Whitney.
!     COLLISION_LOSS = boolean, TRUE if collisional loss is to be calculated,
!                      changess scattering kinematics
!     IN_DELTAS = double precision, *most probable* collisonal loss
!                 before scattering in MeV
!     IN_DELTAP = double precision, *most probable* collisonal loss
!                 after scattering in MeV
!     USERXI = boolean, TRUE if user supplied values of XI should be used
!     USER_XIB = double precision, XI factor before scattering in MeV
!     USER_XIA = double precision, XI factor after scattering in MeV
!     POLARIZED  = boolean, TRUE if polarized radcor is to be calculated,
!                  only for He-3 or Prot
!     STHETA     = double precision, angle of target spin in
!                  degrees (parallel,perpendicular): 180.0 , 270.0
!        *** BOTH INTERN AND PEAK CAN NOT BE TURE AT THE SAME TIME.
!
!     THE ROUTINE READS IN A NUMBER WHICH SELECTS ONE OF THE ABOVE OP-
!     TIONS.  THE NUMBERS AND THE CORRESPONDING KEYWORDS ARE:
!
!     1 : INTERN > INCLUDING THIS CARD TELLS THE PROGRAM TO COMPUTE THE
!                  INTERNAL BREMSTRAHLUNG CORRECTION.
!     2 : PEAK   > INCLUDING THIS CARD TELLS THE PROGRAM TO COMPUTE THE
!                  PEAKING APPROXIMATION TO THE INTERNAL TAIL. THIS AND
!     3 : EXTERN > INCLUDING THIS CARD TELLS THE PROGRAM TO COMPUTE THE
!                  EXTERNAL (i.e., TARGET EFFECT) CORRECTION.
!     4 : MULTIP > INCLUDING THIS CARD TELLS THE PROGRAM TO USE THE
!                  MULTIPLE-PHOTON CORRECTION.
!     5 : RUN    > THIS CARD IS ALWAYS REQUIRED, AND MUST FOLLOW ALL
!                  OTHER KEY WORDS. IT IS FOLLOWED BY THE (1)SCATTERING
!                  ANGLE, INCIDENT ENERGY, AND (2) FINAL ENERGY DATA.
!     -----------------------------------------------------------------

!-----------------------------------------------------------------------
! A interface to C/C++, initialize the rtail program with user defined
! XI related settings and polarization settings
! It also initializes all the COMMON BLOCKS that will be used by the subroutines
! throughout the whole program, the target is fixed to Helium-3 for our case
      SUBROUTINE RTAILS_INIT(xi_bool, xib, xia, pol_bool, pol_theta)
     &BIND(C, name = "rtails_init")
!-----------------------------------------------------------------------
      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT REAL*8 (A-H, O-Z)
      REAL(C_DOUBLE), INTENT(IN), VALUE :: xib, xia, pol_theta
      LOGICAL(C_BOOL), INTENT(IN), VALUE :: xi_bool, pol_bool
      LOGICAL USERXI,POLARIZED,AMROUN_OR_MSW
      REAL*8 USER_XIB,USER_XIA,EINC,EFIN,THETA,M,MT,MTF
      INTEGER IPOL,ITARGET,NKLG,IG

      COMMON /MASS/M,MT,MTF
      COMMON /RUN/EINC,EFIN,THETA,STHETA
      COMMON /TARGET/TWB,TWA,TTB,TTA,ZWB,ZWA,ZT
      COMMON /FMFACCONSTANT/ITARGET
      COMMON /klgcom/rmod,imod,nklg,ig
      COMMON /POLTARG/AMT,TARA,TARZ
      COMMON /CONSTANTS/FINES,ALPHA,PI
      COMMON /XISTUFF/USER_XIB,USER_XIA,USERXI
      COMMON /HE3FF/AMROUN_OR_MSW
      COMMON /LTYPE/POLARIZED
      COMMON /POLINP/IPOL

      DATA  FINES/137.03599911D0/         ! 1/alpha
      DATA  ALPHA/7.297352568D-3/         ! 1.0D0/137.03599911D0
      DATA     PI/3.14159265359D0/        ! tasty
      DATA      M/0.510998918D0/          ! electron mass, MeV

      USERXI = xi_bool
      USER_XIB = xib
      USER_XIA = xia
      POLARIZED = pol_bool
      STHETA = pol_theta

      IF (POLARIZED) THEN
          IF (STHETA.EQ.180) THEN
               IPOL=1
          ELSE IF (STHETA.EQ.270) THEN
               IPOL=3
          ELSE
               IPOL=-1
               PRINT *, 'bad spin angle: ' , STHETA
               PRINT *, '180 for para and 270 for perp,'
               PRINT *, 'assumed unpolarized!'
          ENDIF
      ELSE
          IPOL=-1
      ENDIF

      ! needed by POLRAD
      TARA = 3.D0
      TARZ = 2.D0
      AMT = 3.0160293097D0*931.494043D0/1D3 ! He-3 + electrons only

      ! needed by external
      MT = 3.0160293097D0*931.494043D0
      MTF = MT

      ! needed by radiative corrections
      ! target thickness, included in TW
      TTB = 0D0
      TTA = 0D0
      ! for helium-3
      ZWB = 13D0
      ZWA = 13D0
      ZT = 2D0

      ITARGET = 3
      nklg = 100
      ig = 13

      AMROUN_OR_MSW = .TRUE.

      END SUBROUTINE RTAILS_INIT



!-----------------------------------------------------------------------
! A interface to C/C++, generate the cross section with radiation effects
! It accepts Es Ep in MeV, theta in rad, and radiation length before and after
! It returns sigrad in ub/MeV/sr
      SUBROUTINE RTAILS_RAD_CXSN(es, ep, th, rlin, rlout, sigrad)
     &BIND(C, name = "rtails_rad_cxsn")
!-----------------------------------------------------------------------
      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT REAL*8 (A-H, O-Z)
      REAL(C_DOUBLE), INTENT(IN), VALUE :: es, ep, th, rlin, rlout
      REAL(C_DOUBLE), INTENT(OUT) :: sigrad
      LOGICAL INCREM,POLARIZED
      INTEGER IPOL
      EXTERNAL XSECT
      COMMON /TARGET/TWB,TWA,TTB,TTA,ZWB,ZWA,ZT
      COMMON /RUN/EINC,EFIN,THETA,STHETA
      COMMON /POLTARG/AMT,TARA,TARZ
      COMMON /POLINP/ IPOL
      COMMON /LTYPE/POLARIZED
      COMMON /SINGUL/INCREM
      COMMON /CONSTANTS/FINES,ALPHA,PI
      DATA HBARC/197.326968D-13/

      HCSQ = HBARC**2

      EINC = es
      EFIN = ep
      THETA = th
      TWB = rlin
      TWA = rlout

      CALL FBAR(THETA)
      CALL EXTERNL(XEXTB, XEXTA, FSOFT, XINT)

      RCEXT = HCSQ*FSOFT*(XEXTB + XEXTA)*1.0D20
      IF (POLARIZED) THEN
          CALL POLSIG_EL(EINC, EFIN, THETA, IPOL, SIGMA_BORN, SIGMA_RAD)
          RCINT = SIGMA_RAD
      ELSE
          CALL XSECTP(EINC, EFIN, THETA)
          INCREM = .TRUE.
          RCINT = (ALPHA**3/(4.D0*PI**2))*(EFIN/(EINC*AMT*1.D3))
     &          * SIMPSN(XSECT, -0.999999999D0, 0.999999999D0, 0.00001D0)
     &          * FSOFT * HCSQ * 1.0D30
      ENDIF

      sigrad = RCEXT + RCINT
      END SUBROUTINE RTAILS_RAD_CXSN



!************************************************************************
!   this a version of the subroutines newtail with TWO mods.
!        1. calls  donals nform subroutinne which is inserted at the end of
!           newtailto get the proton form factors.  this will be used
!            to test the difference between various form factors on tail calc.
!                  new common/klgcom/  has variables used by nform
!                      that are read in in rosetail see sub fmfac
!                  insert call to nform in subr fmfac the correct Q**2 units
!                      ect are all in the same section in fmfac
!
!        2.  also change the proton calculation (itargt=2) so that
!             fbar is used (find cklg). previous version had not used fbar so
!             that the hydrogen tail could be compared with tsai
!                    this only required that I comment out an "if" statement

      SUBROUTINE EXTERNL(XEXTB,XEXTA,FSOFT,XINT)

!     -----------------------------------------------------------------
!     COMPUTES THE RADIATIVE CORRECTIONS FOR REAL BREMSSTRAHLUNG AND
!     IONIZATION LOSS IN THE TARGET AND THE TARGET WINDOWS (i.e., THE
!     CONTAINER). THIS IS ALSO KNOWN AS THE EXTERNAL CORRECTION. (STEIN,
!     et.al., PHYS. REV. D VOL.12, NO.7, PP. 1915-1916.) ALSO COMPUTES
!     THE MULTIPLE CORRECTION TERM FSOFT (EQ. A58 IN STEIN) AND CAN, IF
!     DESIRED, COMPUTE THE PEAKING APPROXIMATION TO THE INTERNAL TAIL.
!     THE ROUTINE ALWAYS COMPUTES FSOFT.  IT COMPUTES THE EXTERNAL TAIL
!     ONLY IF EXTERN IS .TRUE., AND COMPUTES THE PEAKING APPROXIMATION
!     TO THE INTERNAL TAIL ONLY IF PEAK IS .TRUE.
!     THE EXTERNAL TAIL IS BROKEN INTO TWO PARTS, XEXTB AND XEXTA, WHICH
!     ARE THE CONTRIBUTIONS FROM, RESPECTIVELY, BEFORE AND AFTER SCAT-
!     TERING.
!     RETURNS CROSS SECTIONS IN UNITS OF HBARC**2/(STR-MEV**3), SO MUL-
!     TIPLYING BY HBARC**2 = (197.3E-13 MEV-CM)**2 GIVES THE OUTPUT
!     CROSS SECTIONS IN UNITS OF CM**2/(STR-MEV).
!-----------------there is a factor of 1E10, -xiaochao
!     THIS ROUTINE GENERALLY FOLLOWS STEINS NOTATION.
!     SUBSCRIPTS S AND P REFER, RESPECTIVELY, TO THE INCIDENT AND OUTGO-
!     ING ELECTRON (e.g., ES AND EP ARE THE ENERGIES OF, RESPECTIVELY,
!     THE INCIDENT AND OUTGOING ELECTRON, IN MEV).
!     SUBSCRIPT T GENERALLY REFERS TO THE TARGET (e.g., MT IS THE MASS
!     OF THE TARGET NUCLEUS IN MEV AND ZT IS THE Z OF THE TARGET).
!     SUBSCRIPT W REFERS TO THE WINDOW.
!     SUBSCRIPTS B OR A (USUALLY APPENDED TO SUBSCRIPTS T OR W) REFER,
!     RESPECTIVELY, TO BEFORE AND AFTER SCATTERING (E.G., TWB IS THE
!     THICKNESS OF THE WINDOW BEFORE SCATTERING, IN RADIATION LENGTHS).
!     THETA IS THE SCATTERING ANGLE (IN RADIANS).
!     M IS THE MASS OF THE ELECTRON IN MEV.
!     TWB, TTB, TTA, AND TWA ARE ALL LENGTHS (IN RADIATION LENGTHS).
!     ZT, ZWB, AND ZWA ARE THE VALUES OF Z.
!     NUMBERS IN THE FAR RIGHT COLUMNS INDICATE THE EQUATION NUMBERS
!     (IN STEIN FROM) WHICH THE MARKED LINES CAME.
!
!     0101115 - vs - implemented from jrosetail (032807 - js):
!     an option for a user specified xi
!     -----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MT,M,MTF,NUS,NUP
      COMMON /MASS/M,MT,MTF
      COMMON /RUN/ES,EP,THETA,STHETA
      COMMON /TARGET/TWB,TWA,TTB,TTA,ZWB,ZWA,ZT
      COMMON /FBARCONSTANT/FBAR1,FBAR2,FBAR3,BTTOT
      COMMON /CONSTANTS/FINES,ALPHA,PI
      COMMON /XISTUFF/USER_XIB,USER_XIA,USERXI
      COMMON /LTYPE/POLARIZED
      LOGICAL POLARIZED,USERXI

!     -----------------------------------------------------------------
!     DEFINE FUNCTIONS USED BY STEIN AND NECESSARY CONSTANTS.
!     -----------------------------------------------------------------
      TERM(Z) = DLOG(184.15D0/(Z**(1.D0/3.D0)))
     @     - 1.2D0*(Z*1/137.03604D0)**2
      ETA(Z,TERMZ) = DLOG(1194.D0/(Z**(2.D0/3.D0)))/TERMZ ! A46
      B(Z,ETAZ,TERMZ) = (4.D0/3.D0)*(1.D0+(Z+1.D0)/(12.D0*(Z+ETAZ)
     @     *TERMZ))             ! A45
      PHI(X) = 1.D0 - X + .75D0*X*X ! A54

!     STEIN'S FORMULA
!     TERM(Z) = DLOG(183.D0/(Z**(1.D0/3.D0)))
!     ETA(Z,TERMZ) = DLOG(1440.D0/(Z**(2.D0/3.D0)))/TERMZ      ! A46
!     B(Z,ETAZ,TERMZ) = (4.D0/3.D0)*(1.D0+(Z+1.D0)/(9.D0*(Z+ETAZ)*
!     @        TERMZ))                                             ! A45


!     -----------------------------------------------------------------
!     COMPUTE TERMS NEEDED FOR ALL CALCULATIONS.
!     -----------------------------------------------------------------
      SINSQ = DSIN(THETA/2.D0)**2
      OMEGAS = ES - EP/(1.D0-2.D0*EP*SINSQ/MT)                    ! A50
      OMEGAP = ES/(1.D0+2.D0*ES*SINSQ/MT) - EP                    ! A51

      NUS = OMEGAS/ES                                             ! A52
      NUP = OMEGAP/(EP+OMEGAP)                                    ! A53

!     -----------------------------------------------------------------
!     STEINS FORMULATION IS HERE EXTENDED TO ACCOUNT FOR THE WINDOW
!     THICKNESSES BEFORE AND AFTER SCATTERING, IN ADDITION TO THE TARGET
!     THICKNESSES BEFORE AND AFTER.
!     -----------------------------------------------------------------
      TERMWB = TERM(ZWB)
      TERMT = TERM(ZT)
      TERMWA = TERM(ZWA)

      ETAWB = ETA(ZWB,TERMWB)
      ETAT = ETA(ZT,TERMT)
      ETAWA = ETA(ZWA,TERMWA)

      BWB = B(ZWB,ETAWB,TERMWB)
      BT = B(ZT,ETAT,TERMT)
      BWA = B(ZWA,ETAWA,TERMWA)

      BTB = BWB*TWB + BT*TTB
      BTA = BT*TTA + BWA*TWA                                      ! A49
      BTTOT = BTB + BTA

!     -----------------------------------------------------------------
!     COMPUTE FSOFT, THE MULTIPLE-PHOTON CORRECTION DEFINED BY STEIN,
!     EQ. A58.
!     -----------------------------------------------------------------
      QSQ = 4.D0*ES*EP*SINSQ
      BTR = (ALPHA/PI)*(DLOG(QSQ/(M*M)) - 1.D0)                   ! A57
      FSOFT = (NUS**(BTB+BTR))*(NUP**(BTA+BTR))                   ! A58

!     -----------------------------------------------------------------
!     COMPUTE TERMS NEEDED FOR BOTH THE EXTERNAL CORRECTION AND THE
!     PEAKING APPROXIMATION.
!     -----------------------------------------------------------------
!      IF ((.NOT.EXTERN) .AND. (.NOT.PEAK)) RETURN

      IF (POLARIZED) THEN
              E2 = (ES-OMEGAS)/1000.0
              E3 = (ES)       /1000.0
              ASYM2 = ASYM_CALC3( E2,THETA, STHETA)
              ASYM3 = ASYM_CALC3( E3,THETA, STHETA)
          ELSE
              ASYM2=1.0D0
              ASYM3=1.0D0
          ENDIF

      TERM1 = (MT + 2.D0*(ES-OMEGAS)*SINSQ)/(MT - 2.D0*EP*SINSQ) !zxc A49
      TERM2 = SIGBAR(ES-OMEGAS)*ASYM2
      TERM3 = SIGBAR(ES)*ASYM3
!      IF(.NOT. EXTERN) GO TO 10

!---- calculating xi - old and sometimes wrong calculation of XI
!---- (1) the radiation logrithms used in STEIN are inaccurate and don't
!         include the Coulomb correction. the effect is at few % level
!         for more info, see: Tsai, Rev. Mod. Phys. 46, 815-51 (1974)
!            and the ERRATUM: Tsai, Rev. Mod. Phys. 49, 421    (1977)
!---- (2) for a composite material and/or "layers" of materials, the total
!         collisional thickness "XI" does not have a simple relationship to
!         the total radiation thickness, T=Tb_Ta
!---- (3) the before and after collisional thickness are different by an order
!         of magnitude for typical He3 experiments!  STEIN justs uses an avg
!---- (4) using the formulas below can result in collisional thicknesses that
!         are too large by factors upto 30! it all depends on how close the
!         target Z is to the effective Z of the materials that contribute most
!         to the radiation thickness. the proper way to calculate the
!         collisional thickness in general is outlined in:
!                 http://www.jlab.org/~singhj/docs/radlength_sagdhv[latest].ps
!---- (5) the !! formulas use more accurate rad logs and distinguish between
!         before and after collisional thicknesses. however, they still assume
!         a simple relationship with the radiation thickness!

      IF (USERXI) THEN
         XIB = USER_XIB
         XIA = USER_XIA

!--------"before" and "after parts of (A49)
         XEXTB = BTB*PHI(NUS)/OMEGAS
         XEXTB = XEXTB + XIB/( OMEGAS**2 )
         XEXTB = XEXTB*TERM2*TERM1
         XEXTB = XEXTB*1.D10

         XEXTA = BTA*PHI(NUP)/OMEGAP
         XEXTA = XEXTA + XIA/( OMEGAP**2 )
         XEXTA = XEXTA*TERM3
         XEXTA = XEXTA*1.D10
      ELSE
!     -----------------------------------------------------------------
!     COMPUTE THE EXTERNAL CORRECTION (EQ. A49).
!     XI HERE IS GENERALIZED TO INCLUDE THE WINDOW.
!     -----------------------------------------------------------------
         XIWB = TWB/((ZWB+ETAWB)*TERMWB)
         XIT = (TTB+TTA)/((ZT+ETAT)*TERMT)
         XIWA = TWA/((ZWA+ETAWA)*TERMWA)
         XI = (PI*M/(2.D0*ALPHA))*(XIWB+XIT+XIWA) ! A52

!     -----------------------------------------------------------------
!     THE TOTAL EXTERNAL TAIL IS XEXT = XEXTA + XEXTB
!     -----------------------------------------------------------------
         XEXTB = 1.D10*(TERM1*TERM2*(BTB*PHI(NUS)/OMEGAS +
     .        XI/(2.D0*OMEGAS*OMEGAS))) ! A49
         XEXTA = 1.D10*(TERM3*(BTA*PHI(NUP)/OMEGAP +
     .        XI/(2.D0*OMEGAP*OMEGAP)))
      ENDIF

!     -----------------------------------------------------------------
!     COMPUTE THE PEAKING APPROX. TO THE INTERNAL CORRECTION (EQ. A56).
!     -----------------------------------------------------------------
!10    IF (.NOT.PEAK) RETURN
!      XINT =1.D10*(TERM1*TERM2*(BTR*PHI(NUS)/OMEGAS) +
!     .	TERM3*(BTR*PHI(NUP)/OMEGAP))


      RETURN
      END

!     *****************************************************************


      SUBROUTINE FBAR(THETA)

!     ------------------------------------------------------------------
!     COMPUTES TERMS USED TO CALCULATE THE FUNCTION FBAR (SEE STEIN,
!     EQ. A43).  THIS SUBROUTINE MUST BE CALLED EACH TIME THETA CHANGES.
!     ------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,M-Z)
      COMMON /FBARCONSTANT/FBAR1,FBAR2,FBAR3,BTTOT
      COMMON /CONSTANTS/FINES,ALPHA,PI

      AP = ALPHA/PI

      SINSQ = DSIN(THETA/2.D0)**2
      COSSQ = DCOS(THETA/2.D0)**2
      PHI = SPENCE(COSSQ,SINSQ)

      FBAR1 = 2.D0*AP*(-14.D0/9.D0) + AP*(PI*PI/6.D0 - PHI)
      FBAR2 = 2.D0*AP*(13.D0/12.D0)
      FBAR3 = - .5D0*AP

      RETURN
      END


!     *****************************************************************

      FUNCTION XSECT(X)

!     -----------------------------------------------------------------
!     CALCULATES THE INTEGRAND FOR THE MO AND TSAI INTERNAL
!     BREMSSTRAHLUNG  CORRECTION (MT, EQ. B.5).
!
!     THIS IS DIFFERENT FROM MO AND TSAI IN ONE WAY: THE FORM FACTORS
!     FJQ AND GJU USED IN MT ARE MULTIPLIED BY A FUNCTION FBAR (AS IN
!     STEIN EQ. A42-A44, OR TSAI 71).  THIS IS NOT DONE IF THE TARGET
!     IS A PROTON (WHICH CASE IS INCLUDED IN THIS PROGRAM ONLY TO TEST
!     MT, WHICH DOES NOT USE FBAR).  FOR A GIVEN CHOICE OF ES, EP AND
!     THETA, SUBROUTINE XSECTP MUST BE CALLED BEFORE USING THIS FUNCTION
!     (XSECTP DEFINES THE VALUES CONTAINED IN COMMON /XSECT/ IN AN
!     ATTEMPT TO MAKE THIS ROUTINE MORE EFFICIENT). IF INCREM IS .TRUE.,
!     THE ROUTINE WILL INCREMENT THETAK TO AVOID SINGULARITIES IN THE
!     INTEGRAND (THIS IS VEVY NECESSARY TO AVOID NUMERICAL ERROR).
!
!     THE NOTATION USED HERE GENERALLY FOLLOWS MT.
!     THE RELEVANT 4-VECTORS ARE:
!       S  = INCIDENT ELECTRON = (ES,SVEC)
!       P  = OUTGOING ELECTRON = (EP,PVEC)
!       P1 = TARGET = (MTR,0)
!       K  = REAL PHOTON = (OMEGA,KVEC)
!       U  = TOTAL MOMENTUM TRANSFER FROM ELECTRON = S + P1 - P
!          = (U0,UVEC)
!       Q  = MOMENTUM TRANSFER TO TARGET = S - P - K.
!     THE METRIC USED IS (A0,AVEC)(B0,BVEC) = A0*B0 - AVEC*BVEC.
!     THE NOTATIONS A3V AND A4V REFER, RESPECTIVELY, TO THE MAGNITUDE OF
!     THE VECTOR AVEC AND THE 4-VECTOR A.
!     -----------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ME,MTR,MTRF,NU,NUPRM,NU1,NU2,NU3,NUPRM1,MESQ,INC
      COMMON /MASS/ME,MTR,MTRF
      COMMON /RUN/ES,EP,THETA,STHETA
      COMMON /SECT/U0,U3V,MESQ,SP,OMEGA1,QSQ1,NU1,NU2,NU3,NUPRM1,A1,
     .             APRM1,B1,BPRM1
      COMMON /TABLE/TERM1,TERM2,TERM4,TERM5,TERM6,TERM7,TERM9,TERM10,NU,
     .              NUPRM,QSQ,OMEGA,FJQ,GJQ,TERMI
      COMMON /FBARCONSTANT/FBAR1,FBAR2,FBAR3,BT
      COMMON /FMFACCONSTANT/ITARGT
      COMMON /LTYPE/POLARIZED
      COMMON /SINGUL/INCREM
      LOGICAL INCREM,POLARIZED
      INTEGER HI,CYCLES
      COMMON /CONSTANTS/FINES,ALPHA,PI   ! js - 121206 - constants are "common" now
!      DATA PI/3.1415926536D0/
!     ----------------------------------------------------------------

      INC = .0003D0
      KOUNT = 0

!     -----------------------------------------------------------------
!     THETAK IS THE ANGLE BETWEEN K3V AND U3V.   X = COS(THETAK).
!     OMEGA IS THE PHOTON ENERGY (MT, P.230).
!     -----------------------------------------------------------------
      THETAK = DACOS(X)
      IF(THETAK .EQ. 0.D0) GO TO 170

100   COSTHK = DCOS(THETAK)
      SINTHK = DSIN(THETAK)
105   OMEDON = U0-U3V*COSTHK
      OMEGA  = OMEGA1/OMEDON

!     ----------------------------------------------------------------
!     QSQ IS Q4V**2.
!     (MT, P.230 IS WRONG. SEE STEIN, A30, OR TSAI 71, P.40).
!     ----------------------------------------------------------------
      QSQ = 2.D0*(QSQ1 - OMEGA*(ES-EP) + OMEGA*U3V*COSTHK)

!     ----------------------------------------------------------------
!     COMPUTE TERMS FROM MT, P. 231.  SUBSCRIPT PRM MEANS PRIMED.
!     ----------------------------------------------------------------
      DENOM = OMEGA*(NU2 + NU3*COSTHK)
      NU = NU1/DENOM
      NUPRM = NUPRM1/DENOM
      A = OMEGA*(EP-A1*COSTHK)
      APRM = OMEGA*(ES-APRM1*COSTHK)
      B = - OMEGA*B1*SINTHK
      BPRM = - OMEGA*BPRM1*SINTHK

!     -----------------------------------------------------------------
!     THE INTEGRAND BECOMES INDEFINITE WHEN A*BPRM = APRM*B (SEE MT,
!     APP. D). WHEN THIS CONDITION BECOMES APPROXIMATELY TRUE, WE
!     INCREASE THETAK TO AVOID A SMALL REGION AROUND THE SINGULARITY
!     (IF INCREMENT IS .TRUE.)
!     KOUNT IS A FLAG THAT INSURES THE INCREMENTS ARE NOT TOO SMALL.
!     -----------------------------------------------------------------
      IF(DABS(A) .LE. DABS(B) .OR. DABS(APRM) .LE. DABS(BPRM)) GO TO 170
      IF (.NOT.INCREM) GO TO 110
      ABPRM = A*BPRM
      APRMB = APRM*B
      DAB   = ABPRM -APRMB
      IF(DAB .EQ. 0.D0) GO TO 1000
110   CONTINUE

      AB  = DSQRT(A*A - B*B)
      ABP = DSQRT(APRM*APRM - BPRM*BPRM)

!     ----------------------------------------------------------------
!     COMPUTE FORM FACTORS FJQ AND GJQ (BEST DEFINED BY MT, EQ. B.3).
!     MULTIPLY THE FORM FACTORS BY FBAR, AS IN STEIN, EQ. A42-A44
!     (UNLESS THE TARGET IS A PROTON).
!     ----------------------------------------------------------------
      CALL FMFAC(QSQ,FJQ,GJQ)
!      write(23,*) -QSQ*1.0D-6, FJQ, GJQ
      IF (ITARGT.EQ.2) GO TO 115
      FBAR = 1.D0 + .5772D0*BT + FBAR1 + FBAR2*DLOG(-QSQ/MESQ) +
     .       FBAR3*DLOG(ES/EP)**2
!      print *,'in fmfac:  ', QSQ,FJQ,GJQ,FBAR
      FJQ  = FJQ*FBAR
      GJQ  = GJQ*FBAR
115   CONTINUE

!     ----------------------------------------------------------------
!     COMPUTE AND SUM TERMS 1 THROUGH 6 FOR THE ELECTRIC FORM FACTOR
!     CONTRIBUTION   STEIN (A24)
!     ----------------------------------------------------------------
      FTERM = 0.D0
      IF(FJQ .LT. 1.D-10) GO TO 120

      TERM1 = (-2.D0*PI*A*MESQ/(AB**3))*(2.D0*ES*(EP+OMEGA)+.5D0*QSQ)
      TERM2 = (-2.D0*PI*APRM*MESQ/(ABP**3))*(2.D0*EP*(ES-OMEGA) +
     .        .5D0*QSQ)
      TERM3 = -4.D0*PI
      TERM4 = 4.D0*PI*(NU/AB - NUPRM/ABP) * (MESQ*(SP-OMEGA*OMEGA)
     .        + SP*(2.D0*ES*EP-SP+OMEGA*(ES-EP)))
      TERM5 = (2.D0*PI/AB)*(2.D0*(ES*EP+ES*OMEGA+EP*EP) +
     .        .5D0*QSQ-SP-MESQ)
      TERM6 = (-2.D0*PI/ABP)*(2.D0*(ES*EP-EP*OMEGA+ES*ES) +
     .        .5D0*QSQ-SP-MESQ)
      FTERM = MTR*MTR*FJQ*(TERM1+TERM2+TERM3+TERM4+TERM5+TERM6)

!     ------------------------------------------------------------------
!     COMPUTE AND SUM TERMS 7 THROUGH 10 FOR THE MAGNETIC FORM FACTOR
!     CONTRIBUITON.
!     ------------------------------------------------------------------

120   GTERM = 0.D0
      IF(GJQ .LT. 1.D-10) GO TO 130
      TERM7 = 2.D0*PI*(A/(AB**3) + APRM/(ABP**3))*MESQ*(2.D0*MESQ+QSQ)
      TERM8 = 8.D0*PI
      TERM9 = 8.D0*PI*(NU/AB - NUPRM/ABP)*SP*(SP-2.D0*MESQ)
      TERM10= 2.D0*PI*(1.D0/AB - 1.D0/ABP)*(2.D0*(SP+MESQ)-QSQ)
      GTERM = GJQ*(TERM7+TERM8+TERM9+TERM10)

!     ------------------------------------------------------------------
!     COMPUTE THE INTEGRAND
!     ------------------------------------------------------------------
130   TERMI = OMEGA/(2.D0*QSQ*QSQ*(U0-U3V*COSTHK))
      XSECT = TERMI*(FTERM+GTERM)
!      PRINT *, 'f,g:' ,FTERM, GTERM
1000  RETURN

!     ------------------------------------------------------------------
!     THE INTEGRAND IS SINGULAR.  INCREMENT THETAK.
!     ------------------------------------------------------------------
140   WRITE(6,2)
2     FORMAT(/' APRMB = ABPRM ENCOUNTERED. ')
      KOUNT = KOUNT + 1
      IF (KOUNT.LE.10) GO TO 180
      KOUNT = 0
      INC = INC*10.D0
      GO TO 180

170   IF(INCREM) GO TO 180
      XSECT = 0.D0
      RETURN

180   THETAK = THETAK + INC
      GO TO 100

      END

!     ****************************************************************

      SUBROUTINE XSECTP(ES,EP,THETA)
!     ----------------------------------------------------------------
!     THE PREPARATION SUBROUTINE FOR XSECT.  CALL THIS ROUTINE, EACH
!     TIME THE CALLING PARAMETERS CHANGE, BEFORE USING XSECT.  THIS
!     DEFINES SOME VARIABLE NEEDED BY XSECT WHICH ARE NOT DEPENDENT ON
!     THETAK (AND THEREFORE SHOULD NOT BE REDEFINED EACH TIME XSECT IS
!     CALLED).
!
!     -----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ME,MTR,MTRF,MESQ
      REAL*8 NU1,NU2,NU3,NUPRM1
      COMMON /MASS/ME,MTR,MTRF
      COMMON /SECT/U0,U3V,MESQ,SP,OMEGA1,QSQ1,NU1,NU2,NU3,NUPRM1,A1,
     .             APRM1,B1,BPRM1
!     -----------------------------------------------------------------
      MESQ = ME*ME
      COSTH = DCOS(THETA)
      SINTH = DSIN(THETA)
      U0  = ES + MTR - EP
      S3V = DSQRT(ES*ES - MESQ)
      P3V = DSQRT(EP*EP - MESQ)
      U3V = DSQRT(P3V*P3V + S3V*S3V - 2.D0*P3V*S3V*COSTH)
      U4V = DSQRT(U0*U0 - U3V*U3V)

      THETAS = DASIN((P3V/U3V)*SINTH)
      THETAP = THETA + THETAS
      SINTHS = DSIN(THETAS)
      COSTHS = DCOS(THETAS)
      SINTHP = DSIN(THETAP)
      COSTHP = DCOS(THETAP)

      OMEGA1 = .5D0*(U4V*U4V - MTRF*MTRF)
      QSQ1 = MESQ - ES*EP + S3V*P3V*COSTH

      NU1 = -P3V*SINTHP
      NU2 = EP*S3V*SINTHS - ES*P3V*SINTHP
      NU3 = S3V*P3V*SINTH
      NUPRM1 = -S3V*SINTHS

      A1 = P3V*COSTHP
      APRM1 = S3V*COSTHS

      B1 = P3V*SINTHP
      BPRM1 = S3V*SINTHS
      SP = ES*EP - S3V*P3V*COSTH
      RETURN
      END

!     ***************************************************************

      SUBROUTINE FMFAC(QSQ,FJQ,GJQ)

!     -----------------------------------------------------------------
!     COMPUTES THE FORM FACTORS FJQ AND GJQ NEEDED BY MO AND TSAI'S
!     INTERNAL BREMSSTRAHLUNG CALCULATION.
!
!     THE FORM FACTORS ARE BEST DEFINED BY EQUATION B.3 IN M & T, WHICH
!     SHOWS THAT THE CROSS-SECTION IN TERMS OF THEM IS:
!
!        (MOTT/4Z**2)(RECOIL FACTOR)(FJQ+(2/M**2)(TAN(THETA/2)**2)GJQ)
!
!     THE RELATIONSHIP BETWEEN THESE FORM FACTORS AND THE ELECTRIC AND
!     MAGNETIC FORM FACTORS GE AND GM (IN MT) IS (EQ. III.2, III.3)_
!
!         FJQ = 4(GE**2 + TAU*GM**2)/(1+TAU)
!         GJQ = -(Q**2)(GM**2)
!     WHERE
!         TAU = -(Q**2)/(4MT**2)
!         MT = MASS OF TARGET (INITIAL MASS)
!         Q = 4-MOMENTUM (MT HAS Q**2 NEGATIVE, Q**2 = E**2 - QVEC**2).
!
!     THE GE AND GM USED THIS WAY IN MT ARE EXACTLY THE SAME AS THE
!     NORMAL ELECTRIC AND MAGNETIC FORM FACTORS ( FRAUENFELDER AND
!     HENLEY EQ.6.22 AND 6.43).  EXCEPT, GE AND GM ARE NORMALIZED TO Z
!     (i.e., GE(0)=GM(0)=Z).
!     -----------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES FORM FACTORS ONE OF FIVE WAYS: BY USING
!     INPUT FORM FACTORS (AND INTERPOLATING), OR BY USING APPROXIMATE
!     ANALYTIC EXPRESSIONS FOR PROTONS, HE3, HE4, OR DEUTERIUM. THE
!     DESIRED FORM FACTOR IS CHOSEN BY SPECIFYING ITARGT (ITARGT = -1
!     TO 5 CORRESPONDS TO THE ABOVE 5 CHOICES IN ORDER).
!
!     THE CONVENTIONS FOR THE CHOICES ARE:
!     1. CARDS  : THE ROUTINE ASSUMES THAT ARRAYS YH AND YH2 CONTAIN
!                 FJQ AND GJQ.
!     2. PROTON : THE ROUTINE COMPUTES GE AND GM USING EQ. III.4 OF MT.
!                 (Z=1, SO THERE IS NO REASON TO THINK ABOUT NORMALIZA-
!                 TION).
!     3. HE3    : THE ROUTINE COMPUTES THE ELECTRIC AND MAGNETIC FORM
!                 FACTORS FCH AND FMAG USING THE ANALYTICAL FORM IN
!                 AMROUN ET AL. NUC. PHYS A579, 596-626 (1994)  (KS Nov.2002)
!
!                 FCH AND FMAG ARE EXACTLY THE
!                 SAME AS GE AND GM, EXCEPT FOR NORMALIZATION:
!                       GE = Z*FCH
!                       GM = Z*(1+K)*FMAG
!                 WHERE K IS THE ANOMALOUS MAGNETIC MOMENT OF HE3 RE-
!                 FERRED TO HE3
!                 (i.e., MAG. MOM. = (1+K)*Z*E*HBAR/(2*MHE3*C)).
!                 K = -4.20 FOR HE3.
!     4. HE4    : EXACTLY AS HE3, EXCEPT FMAG = 0.
!     5. D2     : THE ROUTINE COMPUTES THE FORM FACTORS W1 AND W2 USING
!                 THE ANALYTICAL EXPRESSION IN STEIN, EQS. A9-A12.  W1
!                 AND W2 ARE RELATED TO FJQ AND GJQ BY_
!                       FJQ = 4*W2
!                       GJQ = 4*MT*MT*W1
!                 (COMPARE A3 OF STEIN WITH B3 OF MT.)
!     -----------------------------------------------------------------
!     THIS ROUTINE CAN BE CALLED WITH QSQ = Q**2 EITHER POSITIVE OR
!     NEGATIVE. THE FORM USED IN EACH OF THE FOLLOWING SECTIONS IS THAT
!     USED IN THE REFERENCE PAPER.
!     -----------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ME,MT,MTF
      REAL*8 K,KN,KP,MP
      COMMON /MASS/ME,MT,MTF
      COMMON /FMFACCONSTANT/ITARGT
      common/klgcom/rmod,imod,nklg,ig
      COMMON /HE3FF/AMROUN_OR_MSW
      LOGICAL AMROUN_OR_MSW

      DIMENSION H(6)

!     -----------------------------------------------------------------
!     THE CONSTANTS BELOW ARE USED TO COMPUTE THE HE3 AND HE4 FORM
!     FACTORS K IS THE ANOMALOUS MAGNETIC MOMENT OF HE3.
!     -----------------------------------------------------------------
      DATA K/-4.2D0/

!     -----------------------------------------------------------------
!     SUBSCRIPTS 1, 2, AND 3 REFER, RESPECTIVELY, TO HE3 ELECTRIC, HE3
!     MAGNETIC, AND HE4 ELECTRIC.
!     -----------------------------------------------------------------
      DATA A1,B1,C1,D1,P1,Q01/.675D0,.366D0,.836D0,-6.78D-3,
     .                        .90D0,3.98D0/
      DATA A2,B2,C2/.654D0,.456D0,.821D0/
      DATA A3,B3/.316D0,.675D0/

!     ----------------------------------------------------------------
!     THE FOLLOWING CONSTANTS ARE USED IN THE DEUTERIUM FORM FACTOR
!     COMPUTATION KN AND KP ARE THE ANOMALOUS MAGNETIC MOMENTS OF THE
!     NEUTRON AND PROTON MP IS THE MASS OF THE PROTON (IN GEV).
!     THE CONSTANTS IN ARRAY H ARE USED TO COMPUTE P(Q**2) ON THE WAY
!     TO FIND THE DEUTERIUM FORM FACTORS.
!     ----------------------------------------------------------------
      DATA KN,KP/-1.91348D0,1.7927D0/
      DATA MP/.938256D0/
      DATA H/1.0007D0,1.01807D0,1.05584D0,0.836380D0,0.6864584D0,
     .       0.672830D0/

!     -----------------------------------------------------------------

!      print *,'fmfac begin'

      IF(QSQ.le.0) THEN
        QSQP = DABS(-QSQ)
      ELSE
        QSQP = DABS(QSQ)
      ENDIF
      QSQN = -QSQP

      GO TO (10,20,30,40,80), ITARGT

!     -----------------------------------------------------------------
!     FORM FACTORS WERE READ IN FROM CARDS.  THEY HAVE BEEN TURNED INTO
!     FJQ AND GJQ ALREADY.
!     -----------------------------------------------------------------
10    FJQ = TERP(QSQP,1)
      GJQ = TERP(QSQP,2)
      RETURN

!     -----------------------------------------------------------------
!     PROTON FORM FACTORS FROM MT, EQ. III.4.
!     -----------------------------------------------------------------
20    CONTINUE
      GE = (1.D0 - (QSQN/.71D6))**(-2)
!ccklg
      GM = 2.793D0*GE
!cc
!ccccccccc   use donal form factor routine
	qsqd=qsqp*2.568162e-05  ! (1/197.38)**2  converts Mev to fm-2

 	call nform(ig,rmod,imod,QSQd,gepd,gend,gmpd,gmnd)

        write(23,676) QSQd, GE, GM, gepd,gend,gmpd,gmnd
676     FORMAT(1X,7(1X,1PD12.5,1X))

	if(nklg.eq.1)then
	write(6,'(''  day form factor sub '')')
	write(6,'(''   qsq in fermi and in energy='',2e20.5)')qsqd,qsqp
	write(6,'(''   gepdon, dipole='',2e20.5)')gepd,ge
	endif

	ge=gepd
	gm=gmpd
      GO TO 60

!     ------------------------------------------------------------------
!     HE3 FORM FACTORS FROM AMROUN ET AL. NUC. PHYS A579, 596-626 (1994)
!     NOTICE THAT THIS SECTION REQUIRES Q**2 TO BE POSITIVE AND IN FERMI
!     **(-2) INSTEAD OF MEV**2
!     ------------------------------------------------------------------
30    CONTINUE
      IF (AMROUN_OR_MSW) THEN
         Q = DSQRT(QSQP)/197.32D0 ! fm**-1
         QSQP = Q*Q             ! fm**-2
         FCH = FORMC(QSQP)
         FMAG = FORMM(QSQP)
      ELSE
!     ------------------------------------------------------------------
!     HE3 FORM FACTORS FROM MSW, P. 1403.
!     NOTICE THAT THIS SECTION REQUIRES Q**2 TO BE POSITIVE AND IN FERMI
!     **(-2) INSTEAD OF MEV**2
!     ------------------------------------------------------------------
         Q = DSQRT(QSQP)/197.32D0
         QSQP = Q*Q

         F0 = DEXP(-A1*A1*QSQP) - B1*B1*QSQP*DEXP(-C1*C1*QSQP)
         DF = D1*DEXP(-((Q-Q01)/P1)**2)
         FCH = F0 + DF

         FMAG = DEXP(-A2*A2*QSQP) - B2*B2*QSQP*DEXP(-C2*C2*QSQP)
      ENDIF

!     WRITE(26,*) QSQP*(197.32D0)**2*1.D-6,FCH,FMAG

      GO TO 50

!     ------------------------------------------------------------------
!     HE4 FORM FACTORS FROM MSW, P. 1403.
!     AS ABOVE, REQUIRE Q**2 POSITIVE, IN FERMI**(-2).
!     ------------------------------------------------------------------
40    CONTINUE
      QSQP = QSQP/(197.32D0**2)
      FCH = (1.D0 - (A3*A3*QSQP)**6)*DEXP(-B3*B3*QSQP)
      FMAG = 0.D0

!     ------------------------------------------------------------------
!     FOR HE3 AND HE4, NORMALIZE BY Z
!     ------------------------------------------------------------------
50    CONTINUE
      Z = 2.D0
      GE = Z*FCH
      GM = Z*(1.D0+K)*FMAG

!     ------------------------------------------------------------------
!     FOR PROTON, HE3 AND HE4, SQUARE THE FORM FACTORS
!     ------------------------------------------------------------------
60    CONTINUE
      GESQR = GE*GE
      GMSQR = GM*GM

!     ------------------------------------------------------------------
!     COMPUTE FJQ AND GJQ FROM MT (EQ. III.2 AND III.3, MT).
!     ------------------------------------------------------------------
70    TAU = -QSQN/(4.D0*MT*MT)
      FJQ = 4.D0*(GESQR + TAU*GMSQR)/(1.D0+TAU)
      GJQ = -QSQN*GMSQR

      RETURN

!     ------------------------------------------------------------------
!     DEUTERIUM FORM FACTORS FROM STEIN, EQ. A10.
!     INITIALLY, Q MUST BE IN GEV.
!     ------------------------------------------------------------------
80    CONTINUE
      Q = DSQRT(QSQP)*.001D0
      TAU = Q*Q/(4.D0*MP*MP)

!     ------------------------------------------------------------------
!     COMPUTE THE PROTON TERMS GEP AND GMP (EQ. A7).
!     FIRST FIND P(Q**2) (EQ. A8).
!     ------------------------------------------------------------------
      P = 0.D0
      DO 100 II = 1,6
      I = II - 1
      PROD = 1.
      DO 90 JJ = 1,6
      J = JJ - 1
      IF(J .EQ. I) GO TO 90
      PROD = PROD*(Q-J)/(I-J)
90    CONTINUE
      P = P + H(II)*PROD
100   CONTINUE

      GEP = P/((1.D0+Q*Q/.71D0)**2)
      GMP = (1.D0+KP)*GEP

!     ------------------------------------------------------------------
!     COMPUTE TERMS LEADING TO W1 AND W2 (EQ. A12)
!     ------------------------------------------------------------------
      GMN = KN*GEP
      F1N = TAU*GMN/(1.D0+TAU)
      F1P = (GEP+TAU*GMP)/(1.D0+TAU)
      F2N = GMN/(KN*(1.D0+TAU))
      F2P = (GMP-GEP)/(KP*(1.D0+TAU))
      GP = F1N + F1P
      GS = GP + KN*F2N + KP*F2P

!     ------------------------------------------------------------------
!     COMPUTE FD (EQ. A11).  THIS REQUIRES Q IN INVERSE FERMIS.
!     ------------------------------------------------------------------
      Q = DSQRT(QSQP)/197.3D0
      IF(Q .GT. 1.D-10)
     .FD = (1.580D0/Q)*(DATAN(Q/.930D0) - 2.D0*DATAN(Q/3.19D0)
     .      + DATAN(Q/5.45D0))
      IF(Q .LT. 1.D-10) FD=1.580D0*(1.D0/.930D0-2.D0/3.19D0+1.D0/5.45D0)

!     ------------------------------------------------------------------
!     COMPUTE W1 AND W2 (EQ. A10) AND FJQ GJQ
!     ------------------------------------------------------------------
      W1 = FD*FD*(2.D0/3.D0)*TAU*GS*GS
      W2 = FD*FD*GP*GP + W1
      FJQ = 4.D0*W2
      GJQ = 4.D0*MT*MT*W1

      RETURN
      END

!     *****************************************************************

! FORMC Start ------------------------------------------------------------------
!*****from Amroun et al data****************
!     AMROUN ET AL. NUC. PHYS A579, 596-626 (1994)
!     A. Deur
      REAL*8 FUNCTION FORMC(Q2)
      IMPLICIT REAL*8 (A-H,O-Z)
!     Function to compute a SOG form factor
!     Q2 : momentum transfer squared (fm-2)
!     NR : number of Gaussians
!     GA : Gaussians rms (usually 0.8)
!     RR : Gaussians central positions
!     QQ : Gaussians amplitudes
      real*8 RR(12),QQ(12)
      real*8  Q2,GA,A,B
      DATA NR/12/,GA/0.65D0/
      DATA RR/0.1D0,0.5D0,0.9D0,1.3D0,1.6D0,2.0D0,2.4D0,
     &        2.9D0,3.4D0,4.0D0,4.6D0,5.2D0/
      DATA QQ/.027614D0,.170847D0,.219805D0,.170486D0,
     &        .134453D0,.100953D0,.074310D0,.053970D0,
     +        .023689D0,.017502D0,.002034D0,.004338D0/
      Q=DSQRT(Q2)
      G2=GA*GA
      A=DEXP(-Q2*G2/4.D0)
      S=0.D0
      DO I=1,NR
          B=2.D0*RR(I)**2/G2
          QR=Q*RR(I)
          IF (QR.EQ.0.) THEN
              SS=1.D0+B
          ELSE
              SS=DCOS(QR)+B*DSIN(QR)/QR
          END IF
          SS=QQ(I)/(1.D0+B)*SS
          S=S+SS
      END DO
      FORMC=A*S
      RETURN
      END
! FORMC End ####################################################################

! FORMM Start ------------------------------------------------------------------
!     AMROUN ET AL. NUC. PHYS A579, 596-626 (1994)
!     A. Deur
      REAL*8 FUNCTION FORMM(Q2)
      IMPLICIT REAL*8 (A-H,O-Z)
!     Function to compute a SOG form factor
!     Q2 : momentum transfer squared (fm-2)
!     NR : number of Gaussians
!     GA : Gaussians rms (usually 0.8)
!     RR : Gaussians central positions
!     QQ : Gaussians amplitudes
      real*8 RR(12),QQ(12)
      real*8  Q2,GA,A,B
      DATA NR/12/,GA/0.65D0/
      DATA RR/0.1D0,0.5D0,0.9D0,1.3D0,1.6D0,2.0D0,2.4D0,2.9D0,
     &        3.4D0,4.0D0,4.6D0,5.2D0/
      DATA QQ/.059785D0,.138368D0,.281326D0,.000037D0,
     &        .289808D0,.019056D0,.114825D0,.042296D0,
     +        .028345D0,.018312D0,.007843D0,.000000D0/
      Q=DSQRT(Q2)
      G2=GA*GA
      A=DEXP(-Q2*G2/4.D0)
      S=0.D0
      DO I=1,NR
          B=2.D0*RR(I)**2/G2
          QR=Q*RR(I)
          IF (QR.EQ.0.) THEN
              SS=1.D0+B
          ELSE
              SS=DCOS(QR)+B*DSIN(QR)/QR
          END IF
          SS=QQ(I)/(1.D0+B)*SS
          S=S+SS
      END DO
      FORMM=A*S
      RETURN
      END
! FORMM End ####################################################################

! ASYM_CALC3 Start -------------------------------------------------------------
	REAL*8 FUNCTION ASYM_CALC3(E,A,ASPIN)

        IMPLICIT NONE
	REAL*8 A,ANUM,ASPIN,ASTAR
	REAL*8 DENOM,E,EF,FC,FM,GE,GM,MT,NU,PHISTAR,Q,Q3
	REAL*8 Q3SQ,QFM,T,T2,TAU,THETAQ,VL,VT,VTLP,VTP
        REAL*8 PI,C

!CCC        REAL*8 FORMC_1994,FORMM_1994
        REAL*8 FORMC,FORMM

        PI  = ACOS(-1.)
!SL	MT  = 2.80793                ! 3HE MASS = 2.80793 GEV/C**2
        MT  = 2.8094                 ! Changed. 06/09/03 from 2.8084
        EF  = E/(1+2.*E*(SIN(ABS(A)/2.))**2/MT)
        Q   = (E-EF)*2.*MT
        NU  = E - EF
        TAU = NU/(2.*MT)
        QFM = Q / 0.197328**2.       ! H-BAR C = 0.197328 GEV-FM

!CCC        FC = FORMC_1994(QFM)
!CCC        FM = FORMM_1994(QFM)
        FC = FORMC(QFM)
        FM = FORMM(QFM)

        GE = 2. * FC
        GM = -6.3663 * FM

!	CALCULATE ANGLE BETWEEN THE SPIN DIRECTION AND MOMENTUM TRANSFER.
        Q3SQ   = Q + NU**2
        Q3     = SQRT(Q3SQ)
        THETAQ = -ASIN(EF*SIN(A)/Q3)
!        WRITE(6,*) "Ef,Q3, THETAQ: ",EF,Q3,THETAQ

        ASTAR = ABS( ASPIN*PI/180. - THETAQ )
        IF (ASTAR.GE.0.0.AND.ASTAR.LE.PI) THEN
           PHISTAR = 0.0
        ELSE
           PHISTAR = PI
        ENDIF

        C=180.0/PI

!        WRITE(6,'(A,3F10.3)')
!     &        "THQ,PHIS,THS:",THETAQ*C,PHISTAR*C,ASTAR*C

! NOW FIX PHISTAR AND ASTAR TO E99117 ELASTIC CONDITIONS
!
        PHISTAR = PI
        IF (ASTAR.GE.0.0.AND.ASTAR.LE.PI/2) THEN
           ASTAR = PI-ASTAR
        ENDIF
        PHISTAR = 0
        IF (ASTAR.GE.PI/2.AND.ASTAR.LE.PI) THEN
           ASTAR=PI-ASTAR
        ENDIF
!        WRITE(6,'(A,3F10.3)')
!     &        "THQ,PHIS,THS:",THETAQ*C,PHISTAR*C,ASTAR*C

        T    = TAN( ABS(A)/2. )
        T2   = T**2

        VL   = (Q/Q3SQ)**2
        VT   = Q/(2.0*Q3SQ) + T2
        VTP  = T * SQRT(T2 + Q/Q3SQ)
        VTLP = Q * T/(1.4142*Q3SQ)

        ANUM = 2.*TAU*GM**2*COS(ASTAR)*VTP + 2.0*SQRT(2.0*TAU*
     1         (1.+TAU))*GM*GE*SIN(ASTAR)*COS(PHISTAR)*VTLP
        DENOM = VL*(1.+TAU)*GE**2 + VT*2*TAU*GM**2

        ASYM_CALC3 = -(ANUM/DENOM)
        END
! ASYM_CALC3 End ###############################################################

!     *****************************************************************
      FUNCTION RADLENG(Z,DENS,A)

!     -----------------------------------------------------------------
!     COMPUTES THE RADIATION LENGTH (CM) AS DEFINED BY G. MILLER,
!     (SLAC REP. NO.129, P.97)  EXCEPT, R0 THERE IS CHANGED TO MAKE
!     THE UNITS CORRECT (SEE JACKSON, 2D ED., EQ. 15.49).
!
!     THE MASS DENSITY IS
!        DENS = (ATOMIC DENSITY)(ATOMIC MASS IN G.)
!             = (N)(A/A0)
!     A0   : AVOGADROS NUMBER.
!     MASS : MASS OF THE INCIDENT PARTICLE.
!
!     -----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MASS,N
      DATA HBAR,C,A0/1.0545D-27,3.D+10,6.024D+23/
      DATA MASS/9.11D-28/

!     -----------------------------------------------------------------
      ALPHA = 1.D0/137.03604D0
      N = DENS*A0/A
      R0 = ALPHA*HBAR/(C*MASS)
      TERM = DLOG(191.D0/Z**(1.D0/3.D0)) - 1.2D0*(ALPHA*Z)**2
      ZETA = DLOG(1440.D0/Z**(2.D0/3.D0))/TERM
      RADLENG = 1.D0/(4.D0*N*ALPHA*R0*R0*Z*(Z+ZETA)*TERM)
      PRINT *,'RADLEN called'
      PRINT *,Z,A,DENS,RADLENG
      RETURN
      END

!     *****************************************************************

      FUNCTION SIGBAR(E)

!     -----------------------------------------------------------------
!     COMPUTES THE ELASTIC CROSS SECTION (AS A FUNCTION OF E = INCIDENT
!     ELECTRON ENERGY) TIMES A FUNCTION FBAR(QSQ) DEFINED IN STEIN, EQ.
!     A44. THE CROSS-SECTION IS RETURNED IN UNITS OF
!             HBARC**2/(STR*MEV**2),
!     SO THAT MULTIPLICATION BY HBARC**2 = (197.3E-13 MEV-CM)**2 GIVES
!     THE CROSS SECTION IN CM**2/STR.
!     -----------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ME,MT,MTF
      COMMON /MASS/ME,MT,MTF
      COMMON /RUN/ES,EP,THETA,STHETA
      COMMON /FBARCONSTANT/FBAR1,FBAR2,FBAR3,BT
      COMMON /LTYPE/POLARIZED
      COMMON /CONSTANTS/FINES,ALPHA,PI ! js - 121206 - common
      LOGICAL POLARIZED

!     -----------------------------------------------------------------
!     THE ELASTIC CROSS-SECTION COMPUTED HERE USES MO AND TSAIS FORM
!     FACTORS FJQ AND GJQ (SEE EQ. B3 IN MT).
!     THESE ARE RELATED TO W1 AND W2 OF STEIN (EQ. A3) BY_
!           FJQ = 4*W2
!           GJQ = 4*MT*MT*W1
!     (COMPARE B3 OF MT WITH A3 OF STEIN).
!     -----------------------------------------------------------------

      SINSQ = DSIN(THETA/2.D0)**2
      RECOIL = 1.D0/(1.D0 + 2.D0*E*SINSQ/MT)
      EPRM = E*RECOIL
      QSQ = 4.D0*E*EPRM*SINSQ
      CALL FMFAC(QSQ,FJQ,GJQ)

!     ----------------------------------------------------------------
!     XMOTT HERE IS REALLY THE MOTT CROSS-SECTION DIVIDED BY 4Z**2.
!     ----------------------------------------------------------------
      XMOTT = ALPHA*EPRM*DCOS(THETA/2.D0)/QSQ
      XMOTT = XMOTT*XMOTT
      XELAST = XMOTT*RECOIL*(FJQ + 2.D0*GJQ*(DTAN(THETA/2.D0)/MT)**2)

!   ----------------------------------------------------------------
!     COMPUTE FBAR (STEIN, EQ. A44).
!     ----------------------------------------------------------------
      FBAR = 1.D0 + .5772D0*BT + FBAR2*DLOG(QSQ/(ME*ME)) +
     .       FBAR3*DLOG(ES/EP)**2 + FBAR1


!     ----------------------------------------------------------------
!     COMPUTE SIGMA BAR_
!     ----------------------------------------------------------------
      SIGBAR = XELAST*FBAR
!       print *,'elastic:', eprm, qsq, xmott, fjq, gjq, xelast, sigbar

      RETURN
      END


!     *****************************************************************

      FUNCTION SPENCE(X,Y)

!     ------------------------------------------------------------------
!     COMPUTES THE SPENCE FUNCTION
!
!           PHI(X) = - INTEGRAL FROM 0 TO X OF (LOG(1-T))/T
!
!     AS DEFINED IN STEIN, EQ. A48.
!     Y EQUALS 1-X.  IT IS DEFINED OUTSIDE THE ROUTINE IN CASE IT CAN BE
!     COMPUTED MORE EXACTLY (NECESSARY IF X IS NEAR 1).
!
!     THE SPENCE INTEGRAL FOR N=2 (THE DILOGARITHM) IS DEFINED IN
!     ABRAMOWITZ ANDSTEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, AS (EQ.
!     27.7.1)
!
!           F(X) = - INTEGRAL FROM 1 TO X OF (LOG(T))/(T-1).
!
!     COMPARISON SHOWS THAT
!
!           PHI(X) = F(1-X)    FOR 0^X^1.
!
!     SO I CAN ADAPT THE SERIES EXPANSION OF F(X) (AS, EQ. 27.7.2)
!           F(X) = SUM FROM 1 TO INFINITY OF ((1-X)**K)/(K**2))
!     VALID FOR 0 < X < 2.
!
!     THE SERIES EXPANSION FOR PHI(X), VALID FOR 0 < X < 1, IS THEREFORE
!
!           PHI(X) = SUM FROM 1 TO INFINITY OF (X**K)/K**2.
!
!     THIS ROUTINE ALSO USES IDENTITY 27.7.3 FROM AS TO INCREASE THE
!     EFFICIENCY OF COMPUTATION FOR LARGE X (I.E., X NEAR 1.)
!     ------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      DATA MAXIT,EPS/10000,1.D-10/
      DATA PI/3.14159265359D0/

!     ------------------------------------------------------------------
      IF((X .LT. 0.D0) .OR. (X .GT. 1.D0)) STOP 'SPENCE'

      Z = X
      IF(X .GT. .5D0) Z = Y

      SPENCE = Z
      IF(Z .EQ. 0.D0) GO TO 30

      DO 10 K = 2,MAXIT
      ADD = (Z**K)/(K*K)
      SPENCE = SPENCE + ADD
      CHANGE = DABS(ADD/SPENCE)
      IF(CHANGE .LT. EPS) GO TO 20
10    CONTINUE
20    CONTINUE
      IF(X .LE. .5D0) RETURN

!     -----------------------------------------------------------------
!     USE IDENTITY 27.7.3 FROM AS FOR LARGE X (.5{X^1).
!     -----------------------------------------------------------------
      SPENCE = -SPENCE - DLOG(X)*DLOG(Y) + PI*PI/6.D0
      RETURN

!     -----------------------------------------------------------------
!     SPECIAL CASES_  X=0 OR X=1.
!     -----------------------------------------------------------------
30    IF(X .EQ. 1.D0) SPENCE = PI*PI/6.D0

      RETURN
      END

!     *****************************************************************


      FUNCTION SIMPSN(ARG,Y1,Y2,FERR)

!     -----------------------------------------------------------------
!     SIMPSN INTEGRATION ROUTINE WRITTEN AS FORTRAN IV FUNCTION J.SMITH
!     THIS ROUTINE BREAKS THE INTERVAL UP INTO 3 SUBINTERVALS,WHICH
!     ARE AGAIN BROKEN UP INTO 3 INTERVALS, ETC. THE INTERVALS ARE DONE
!     SEPARATLY SO AS TO REDUCE THE NUMBER OF ITERATIONS. GOOD WHERE
!     FUNCTIONAL DEPENDENCE HAS PEAK LIKE STRUCTURE. REDUCES NUMBER OF
!     ITERATIONS BY ABOUT A FACTOR OF 5.
!     -----------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F2T(20),FMT(20),F3T(20),F4T(20),FBT(20),
     1DXT(20),X1T(20),X2T(20),ART(20),EPST(20),ES2T(20),
     2ES3T(20),LEG(20),SUM1(20),SUM2(20)

!     -----------------------------------------------------------------
!     INITIAL SET-UP
!     -----------------------------------------------------------------
      A=Y1
      EPS=FERR
      B=Y2
      DA=B-A
      FA=ARG(A)
      FM=4.D0*ARG((A+B)*.5D0)
      FB=ARG(B)
      AREA=1.0D0
      EST=1.0D0
      L=1

!     ----------------------------------------------------------------
!     BEGIN SIMPSON
!     ----------------------------------------------------------------
    1 DX=DA/3.D0
      X1=A+DX
      X2=X1+DX
      F1=4.D0*ARG(A+.5D0*DX)
      F2=ARG(X1)
      F3=ARG(X2)
      F4=4.D0*ARG(A+2.5D0*DX)
      DX6=DX/6.D0
      EST1=(FA+F1+F2)*DX6
      EST2=(F2+FM+F3)*DX6
      EST3=(F3+F4+FB)*DX6
      AREA=AREA-DABS(EST)+DABS(EST1)+DABS(EST2)+DABS(EST3)
      SUM=EST1+EST2+EST3

!     ----------------------------------------------------------------
!     TEST FOR CONVERGENCE
!     ----------------------------------------------------------------
      IF(DABS(EST-SUM)-EPS*AREA)2,2,3
    2 IF(EST-1.0)6,3,6
    3 IF(L-20)5,6,6
    5 L=L+1
      LEG(L)=3

!     ---------------------------------------------------------------
!     STORE PARAMETERS FOR SIMPSON II AND III
!     ---------------------------------------------------------------
      F2T(L)=F2
      FMT(L)=FM
      F3T(L)=F3
      F4T(L)=F4
      FBT(L)=FB
      DXT(L)=DX
      X1T(L)=X1
      X2T(L)=X2
      ART(L)=AREA
      EPST(L)=EPS/1.7D0
      ES2T(L)=EST2
      ES3T(L)=EST3

!     ---------------------------------------------------------------
!     RETURN TO SIMPSON I
!     ---------------------------------------------------------------
      DA=DX
      FM=F1
      FB=F2
      EST=EST1
      EPS=EPST(L)
      GO TO 1
    6 IF(LEG(L)-2) 9,8,7
    7 SUM1(L)=SUM
      LEG(L)=2

!     ---------------------------------------------------------------
!     RETURN TO SIMPSON II
!     ---------------------------------------------------------------
      A=X1T(L)
      DA=DXT(L)
      FA=F2T(L)
      FM=FMT(L)
      FB=F3T(L)
      AREA=ART(L)
      EST=ES2T(L)
      EPS=EPST(L)
      GO TO 1
    8 SUM2(L)=SUM
      LEG(L)=1

!     ---------------------------------------------------------------
!     RETURN TO SIMPSON III
!     ---------------------------------------------------------------
      A=X2T(L)
      DA=DXT(L)
      FA=F3T(L)
      FM=F4T(L)
      FB=FBT(L)
      AREA=ART(L)
      EST=ES3T(L)
      EPS=EPST(L)
      GO TO 1
    9 SUM=SUM1(L)+SUM2(L)+SUM
      L=L-1
      IF(L-1)11,11,6
   11 SIMPSN = SUM
      RETURN
      END

!     *****************************************************************

      FUNCTION TERP(XIN,ICHANN)

!     ----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /XEPL/XH(400),YH(400),XH2(400),YH2(400)
      COMMON /CARDS/NPTS
      DIMENSION X(400),Y(400),DELTA(10),A(10)

!     -----------------------------------------------------------------
!     NTERMS-POINT INTERPOLATION ROUTINE FROM BEVINGTON (PROGRAM 13-2),
!     SLIGHTLY MODIFIED.
!     PLACE PROPER CHANNEL INTO X AND Y.
!     -----------------------------------------------------------------

      NTERMS = 3
      IF(ICHANN .NE. 1) GO TO 20
      DO 10 I = 1,NPTS
      X(I) = XH(I)
      Y(I) = YH(I)
10    CONTINUE
      GO TO 40
20    DO 30 I = 1,NPTS
      X(I) = XH2(I)
      Y(I) = YH2(I)
30    CONTINUE
40    CONTINUE

!     ---------------------------------------------------------------
!     CHECK LIMITS
!     ---------------------------------------------------------------
      IF((XIN .GE. X(1)) .AND. (XIN .LE. X(NPTS))) GO TO 60
      TERP = 0.D0
      RETURN

!     ---------------------------------------------------------------
!     SEARCH FOR APPROPRIATE VALUE OF X1.  AS SOON AS X(I) IS GREATER
!     THAN OR EQUAL TO XIN, SET I1 AND EXIT LOOP.
!     ---------------------------------------------------------------
60    CONTINUE
      DO 90 I = 1,NPTS
      IF(XIN-X(I)) 70,80,90

!     ---------------------------------------------------------------
!     XIN IS LESS THAN X(I).
!     ---------------------------------------------------------------
70    I1 = I - NTERMS/2
      IF(I1 .LE. 0) I1 = 1
      GO TO 100

!     ---------------------------------------------------------------
!     XIN EQUALS X(I)
!     ---------------------------------------------------------------
80    TERP = Y(I)
      RETURN

!     ---------------------------------------------------------------
!     XIN IS GREATER THAN X(I).
!     IF EXIT LOOP WITHOUT SUCCESS, USE THE LAST NTERMS POINTS IN THE
!     ARRAY.  (WITH THE INITIAL TEST, THIS WILL NOT HAPPEN).
!     ---------------------------------------------------------------
90    CONTINUE
      I1 = NPTS - NTERMS + 1
100   I2 = I1 + NTERMS - 1
      IF(I2 .LE. NPTS) GO TO 110
      I2 = NPTS
      I1 = I2 - NTERMS + 1
      IF(I1 .GE. 1) GO TO 110
      I1 = 1
      NTERMS = I2 - I1 + 1

!     ---------------------------------------------------------------
!     EVALUATE DEVIATIONS DELTA
!     ---------------------------------------------------------------
110   DENOM = X(I1+1) - X(I1)
      DELTAX = (XIN-X(I1))/DENOM
      DO 120 I = 1,NTERMS
      IX = I1 + I -1
      DELTA(I) = (X(IX)-X(I1))/DENOM
120   CONTINUE

!     ---------------------------------------------------------------
!     ACCUMULATE COEFFICIENTS A
!     ---------------------------------------------------------------
      A(1) = Y(I1)
      DO 140 K = 2,NTERMS
      PROD = 1.D0
      SUM = 0.D0
      IMAX = K-1
      IXMAX = I1 + K - 1
      DO 130 I = 1,IMAX
      J = K - I
      PROD = PROD*(DELTA(K) - DELTA(J))
      SUM = SUM - A(J)/PROD
130   CONTINUE
      A(K) = SUM + Y(IXMAX)/PROD
140   CONTINUE

!     ------------------------------------------------------------------
!     ACCUMULATE SUM OF EXPANSION
!     ------------------------------------------------------------------
      SUM = A(1)
      DO 160 J = 2,NTERMS
      PROD = 1.D0
      IMAX = J - 1
      DO 150 I = 1,IMAX
150   PROD = PROD*(DELTAX - DELTA(I))
      SUM = SUM + A(J)*PROD
160   CONTINUE

      TERP = SUM
      RETURN
	end


! this version of nform was copied from [adam.model] feb 24 1988
! by klg

      SUBROUTINE NFORM(IG,RMOD,IMOD,QQ,GEP,GEN,GMP,GMN)
!
!  ----------------------------------------------------------
!
!   CALCULATE NUCLEON FORM FACTORS
!
!
!   QQ = INPUT Q SQUARED (FM**-2)
!	RMOD IS MODIFIED NUCLEON RADIUS TO BE USED WITH IG=11
!	IMOD IS FLAG WHETHER TO MODIFY THE FORM FACTOR
!	IMOD = 0 >> NO MODIFICATON
!	IMOD = 1 >> MODIFICATION BY THE TERM FACTOR BELOW
!  ------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
!
      COMMON/IO/INPUT,IOUT
!
!
      DIMENSION TITLE(15)
      CHARACTER*40 TITLE
      DATA TITLE(1)/'EMPERICAL DIPOLE WITH GEN = 0.0         '/,
     *     TITLE(2)/'IJL 5 PARAMETER MODEL DIPOLE FIT        '/,
     *     TITLE(3)/'BEST FIT=GEP+GMP=IJL; GMN=GD;GEN=GALSTER'/,
     *     TITLE(4)/'BEST FIT EXCEPT GEN= 0.0                '/,
     *     TITLE(5)/'POINT NUCLEONS                          '/,
     *     TITLE(6)/'HOEHLER,NPB114,505                      '/,
     *     TITLE(7)/'BLATNIK + ZOVKO VDM FIT                 '/,
     *     TITLE(8)/'JANSSENS 1966 STANDARD FIT              '/,
     *     TITLE(9)/'HOEHLER FIT NO. 5.1                     '/,
     *     TITLE(10)/'HOEHLER FIT NO. 5.3                    '/,
     *	   title(11)/'Av. of BZ and Dipole for Gmn; BZ else  '/,
     *	   title(12)/'Simon GEP,GMN NPA333(1980)381; BZ else '/,
     *	   title(13)/'Simon-Gep,Gmp;Gen BZ;Gmn Av. BZ+dipole '/,
     *     TITLE(14)/'HOEHLER*NUCLEON SIZE CHANGE *GE ONLY*  '/,
     *     TITLE(15)/'DIPOLE WITH MODIFIED MAGNETIC MOMENT   '/
!_______________________________________________________________
!   IJL PARAMETERS FOR 5 PARAMETER DIPOLE FIT (IN GEV UNITS)
!   PHYS LETT. 43B, 191(1973)
      DATA GAM   ,BR    ,BW    ,BF    ,AF
     +   /0.25  ,0.672 ,1.102 ,0.112 ,-.052 /
      DATA RMN2  ,RMW2  ,RMF2  ,RMR2  ,GAMR  ,PI
     +   /0.8817,0.6146,1.0384,0.5852,0.112 ,3.14159/
      DATA RMPI  ,RMPI2
     +   /.139  ,.019321/
!     BZ
      DATA TRO   ,TROP  ,TROPP ,TFI   ,TOM   ,TOMP
     +   /0.585 ,1.30  ,2.10  ,1.039 ,0.614 ,1.40  /
      DATA RMUS  ,RMUV  ,BS    ,BV
     +   /-.060 ,1.853 ,-.91  ,-1.10 /
      DATA RMUP  ,RMUN
     +   /2.79278,-1.91315/
!	data for simon parametrization
!  Les masses sont en fm-2
      DATA A1S/0.312/,A2S/1.312/,A3S/-.709/,A4S/0.085/
      DATA A1V/0.694/,A2V/0.719/,A3V/-.418/,A4V/0.005/
      DATA AM1/6./,AM2/15.02/,AM3/44.08/,AM4/154.2/
      DATA AM5/8.5/,AM6/355.4/,AMUP/2.793/
!	end data for simon
!
!     QQ IN FM-2,QQG IN GEV/C**2
!
!      print *,'nform: begin'

	IF (IG .EQ. 15)THEN
	RMUP = 3.631
	RMUN = -2.487
	ENDIF
      QQG=QQ*0.197328**2
!      print *, QQG, 'GeV/c**2', QQ,'fm-2'
      	IF (IMOD .EQ. 0)THEN
	FACTOR = 1.0
	ELSE
	IF (IMOD .EQ. 1) THEN
	FACTOR = (1.+QQ*0.80**2/12.)**2/(1.+QQ*RMOD**2/12.)**2
	ELSE
	STOP 'Meaningless value of IMOD stopped in NFORM'
	ENDIF
	ENDIF
      GOTO(110,120,120,120,150,160,170,180,190,200,110,
     >170,110,210,220),IG
!     DIPOLE
  110 AAA=0.71
      GEP=1./(1.+QQG/AAA)**2
      GEN=0.
      GMP=RMUP*GEP
      GMN=RMUN*GEP
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
	if(ig .eq. 13)goto 230
	if(ig .eq. 11)goto 230
      GOTO 900
!     IJL
  120 TAU=QQG/(4.*RMN2)
      GT=0.5/(1.+GAM*QQG)**2
      T1=SQRT(QQG+4.*RMPI2)
      T2=SQRT(QQG)
      ALPH=2.*T1*LOG((T1+T2)/(2.*RMPI))/(T2*PI)
      TOP=RMR2+8.*GAMR*RMPI/PI
      BOT=RMR2+QQG+(4.*RMPI2+QQG)*GAMR*ALPH/RMPI
      RHO=TOP/BOT
      F1S=GT*((1.-BW-BF)+BW/(1.+QQG/RMW2)+BF/(1.+QQG/RMF2))
      F1V=GT*((1.-BR)+BR*RHO)
      F2S=GT*((-0.12-AF)/(1.+QQG/RMW2)+AF/(1.+QQG/RMF2))
      F2V=GT*(3.706*RHO)
      GEP=F1V+F1S-TAU*(F2V+F2S)
      GEN=F1S-F1V-TAU*(F2S-F2V)
      GMP=F1V+F1S+F2V+F2S
      GMN=F1S-F1V+F2S-F2V
      IF(IG.EQ.2) GOTO 900
      GD=1./(1.+QQG/.71)**2
      GMN=RMUN*GD
      GEN=-RMUN*TAU*GD/(1.+5.6*TAU)
      IF(IG.EQ.3) GOTO 900
      GEN=0.
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
  150 CONTINUE
      GMN=RMUN
      GEP=1.
      GEN=0.
      GMP=RMUP
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
  160 CONTINUE
!     HOEHLER FIT 8.2
      TAU=QQG/(4.*RMN2)
      F1S= 0.71/(0.613+QQG)-0.64/(1.040+QQG)-0.13/(3.24+QQG)
      F2S= -0.11/(0.613+QQG)+0.13/(1.040+QQG)-0.02/(3.24+QQG)
      F1V=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536)
      F2V=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)
      F1V=F1V+0.05/(1.464+QQG)-0.52/(6.003+QQG)+0.28/(8.703+QQG)
      F2V=F2V-1.99/(1.464+QQG)+0.20/(6.003+QQG)+0.19/(8.703+QQG)
      GEP=F1V+F1S-TAU*(F2V+F2S)
      GEN=F1S-F1V-TAU*(F2S-F2V)
      GMP=(F1V+F1S+F2V+F2S)
      GMN=(F1S-F1V+F2S-F2V)
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
!     BZ
  170 TAU=QQG/(4.*RMN2)
      RS=1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))
      RV=1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP))
      F1E=(0.5-TAU*(RMUS+2.*RMN2*BS))*RS
      F2E=(0.5-TAU*(RMUV+2.*RMN2*BV))*RV
      F1M=(0.5+RMUS-0.5*BS*QQG)*RS
      F2M=(0.5+RMUV-0.5*BV*QQG)*RV
      GEP=F1E+F2E
      GMP=F1M+F2M
      GEN=F1E-F2E
      GMN=F1M-F2M
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
	if (ig .eq. 12)goto 240 ! if simon use gmn,gen from B+Z
      GOTO 900
!     JANSENS
  180 F1=1.+QQ/15.7
      F2=1.+QQ/26.7
      F3=1.+QQ/8.19
      GES=0.5*(2.5/F1-1.6/F2+0.1)
      GMS=0.44*(3.33/F1-2.77/F2+0.44)
      GEV=0.5*(1.16/F3-0.16)
      GMV=2.353*(1.11/F3-0.11)
      GEP=GES+GEV
      GMP=GMS+GMV
      GEN=GES-GEV
      GMN=GMS-GMV
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GO TO 900
!     HOELER FIT 5.1
  190 TAU=QQG/(4.*RMN2)
      F1=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536)
      F2=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)
      F1P=F1+0.63/(0.613+QQG)-0.43/(1.124+QQG)-0.45/(2.723+QQG)
      F2P=F2+0.02/(0.613+QQG)-1.89/(1.323+QQG)+0.27/(7.952+QQG)
      GEP=F1P-TAU*F2P
      GEN=0.0
      GMP=F1P+F2P
      GMN=0.0
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
!     HOEHLER FIT 5.3
  200 TAU=QQG/(4.*RMN2)
      F1=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536)
      F2=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)
      F1P=F1+0.67/(0.613+QQG)-0.39/(0.922+QQG)-0.54/( 2.756+QQG)
      F2P=F2+0.04/(0.613+QQG)-1.88/(1.300+QQG)+0.24/(10.176+QQG)
      GEP=F1P-TAU*F2P
      GEN=0.0
      GMP=F1P+F2P
      GMN=0.0
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
!     HOEHLER FIT 8.2
!	Here we modify only the electric form factors leave the
!	magnetic unchanged so this is ig=6 with ge modified
210   TAU=QQG/(4.*RMN2)
      F1S= 0.71/(0.613+QQG)-0.64/(1.040+QQG)-0.13/(3.24+QQG)
      F2S= -0.11/(0.613+QQG)+0.13/(1.040+QQG)-0.02/(3.24+QQG)
      F1V=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536)
      F2V=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)
      F1V=F1V+0.05/(1.464+QQG)-0.52/(6.003+QQG)+0.28/(8.703+QQG)
      F2V=F2V-1.99/(1.464+QQG)+0.20/(6.003+QQG)+0.19/(8.703+QQG)
      GEP=F1V+F1S-TAU*(F2V+F2S)
      GEN=F1S-F1V-TAU*(F2S-F2V)
      GMP=(F1V+F1S+F2V+F2S)
      GMN=(F1S-F1V+F2S-F2V)
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
!     GMP=FACTOR*GMP
!     GMN=FACTOR*GMN
      GOTO 900
!	here all is modified by we also change the magnetic moment
  220 AAA=0.71
!      print *,'dipole used for proton'
      GEP=1./(1.+QQG/AAA)**2
      GEN=0.
      GMP=RMUP*GEP
      GMN=RMUN*GEP
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
  230 TAU=QQG/(4.*RMN2)
      RS=1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))
      RV=1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP))
      F1E=(0.5-TAU*(RMUS+2.*RMN2*BS))*RS
      F2E=(0.5-TAU*(RMUV+2.*RMN2*BV))*RV
      F1M=(0.5+RMUS-0.5*BS*QQG)*RS
      F2M=(0.5+RMUV-0.5*BV*QQG)*RV
      GEP=F1E+F2E
      GMP=F1M+F2M
      GEN=F1E-F2E
      GMN=(F1M-F2m + gmn)/2.    ! using gmn from dipole to av with BZ gmn
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
        print *, QQG,GEP,GEN,GMP,GMN

	if(ig .eq. 13) goto 240
      GOTO 900
240	continue
!	form factor of simon
!	fit of simon NPA 333 (1980) 381
!	parametrizations in terms of 4 poles
!	use gmn,gen from b&z above
	q2 = qq	!momentum transfer in fm-2
      GEP=A1S/(1.D0+Q2/AM1)+A2S/(1.D0+Q2/AM2)+A3S/(1.D0+Q2/AM3)
     ++A4S/(1.D0+Q2/AM4)
      GMP=A1V/(1.D0+Q2/AM5)+A2V/(1.D0+Q2/AM2)+A3V/(1.D0+Q2/AM3)
     ++A4V/(1.D0+Q2/AM6)
      GMP=GMP*AMUp
	goto 900

!
!
!
!      ENTRY NFORMI(IG,RMOD,IMOD)
!  ----------------------------------------------------------
!
!   ENTRY POINT TO WRITE FORM FACTOR LABEL IN THE OUTPUT
!
!  ----------------------------------------------------------
!
	iout=6
        print *, QQG,GEP,GEN,GMP,GMN

      print *,'nformi: end'

      WRITE(iout,9000)IG,TITLE(IG)
 9000 FORMAT(1H0,'IG = ',I5,4X,'NUCLEON FORM FACTORS USED = ',A40)

!	IF (IMOD .EQ. 1)WRITE(iout,9001)RMOD
!9001	FORMAT(2X,' These form factors have been modified by a ',
!     >	'change in the nucleon size: RMOD = ',F5.2)
!
        print *, QQ,GEP,GEN,GMP,GMN

      print *,'nformi: end'

      print *,'nform: end'
 900  RETURN
      END


