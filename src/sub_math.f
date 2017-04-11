!-----These subroutines from Seonho Choi. 
!     Downloaded on 12/07/06 from Karl Slifer's website by Jaideep Singh
!     changed all comment line markers to !

      REAL*8 FUNCTION DGAUSS(F,A,B,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGAUSS')
      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X
     1        /0.96028 98564 97536 23168 35608 68569 47D0,
     2         0.79666 64774 13626 73959 15539 36475 83D0,
     3         0.52553 24099 16328 98581 77390 49189 25D0,
     4         0.18343 46424 95649 80493 94761 42360 18D0,
     5         0.98940 09349 91649 93259 61541 73450 33D0,
     6         0.94457 50230 73232 57607 79884 15534 61D0,
     7         0.86563 12023 87831 74388 04678 97712 39D0,
     8         0.75540 44083 55003 03389 51011 94847 44D0,
     9         0.61787 62444 02643 74844 66717 64048 79D0,
     A         0.45801 67776 57227 38634 24194 42983 58D0,
     B         0.28160 35507 79258 91323 04605 01460 50D0,
     C         0.95012 50983 76374 40185 31933 54249 58D-1/

      DATA W
     1        /0.10122 85362 90376 25915 25313 54309 96D0,
     2         0.22238 10344 53374 47054 43559 94426 24D0,
     3         0.31370 66458 77887 28733 79622 01986 60D0,
     4         0.36268 37833 78361 98296 51504 49277 20D0,
     5         0.27152 45941 17540 94851 78057 24560 18D-1,
     6         0.62253 52393 86478 92862 84383 69943 78D-1,
     7         0.95158 51168 24927 84809 92510 76022 46D-1,
     8         0.12462 89712 55533 87205 24762 82192 02D0,
     9         0.14959 59888 16576 73208 15017 30547 48D0,
     A         0.16915 65193 95002 53818 93120 79030 36D0,
     B         0.18260 34150 44923 58886 67636 67969 22D0,
     C         0.18945 06104 55068 49628 53967 23208 28D0/

      H = 0
      IF(B .EQ. A) GO TO 99
      CONST = CST/ABS(B-A)

      BB = A
    1 AA = BB
      BB = B

    2 C1 = HF*(BB+AA)
      C2 = HF*(BB-AA)

      S8=0
      DO 3 I = 1,4
        U  = C2*X(I)
    3 S8 = S8 + W(I)*( F(C1+U) + F(C1-U) )

      S16=0
      DO 4 I = 5,12
        U =C2*X(I)
    4 S16 = S16+W(I)*(F(C1+U)+F(C1-U))
      S16 = C2*S16

      IF ( ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H = H+S16
       IF (BB .NE. B) GO TO 1
      ELSE
       BB = C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
       H = 0
!       CALL MTLPRT(NAME,'D103.1','TOO HIGH ACCURACY REQUIRED')
       GO TO 99
      END IF

   99 DGAUSS=H
      RETURN
      END

      !------------------------------------------------------------------
      SUBROUTINE INTERPOL_PARA(XP,YP,EYP,X,Y,EY,ITYPE_ERR)
      IMPLICIT NONE
      REAL*8 XP,YP,EYP,X,Y,EY
      DIMENSION XP(*),YP(*),EYP(*)
      INTEGER ITYPE_ERR

      REAL*8 A,B,C,D1,D2,D3,D12,D13,D23

      D1 = X - XP(1)
      D2 = X - XP(2)
      D3 = X - XP(3)

      D12 = XP(1) - XP(2)
      D13 = XP(1) - XP(3)
      D23 = XP(2) - XP(3)

      A = D1*D2/(D13*D23)
      B = D2*D3/(D12*D13)
      C = D3*D1/(-D23*D12)

      Y = A*YP(3)+B*YP(1)+C*YP(2)
      IF (ITYPE_ERR.EQ.-1) THEN
         EY = -1.0
      ELSE IF (ITYPE_ERR.EQ.0) THEN
         EY = (A*EYP(3))**2 + (B*EYP(1))**2 + (C*EYP(2))**2
         EY = SQRT(EY)
      ELSE IF (ITYPE_ERR.EQ.1) THEN
         EY = A*EYP(3) + B*EYP(1) + C*EYP(2)
      ENDIF
      RETURN
      END


!------------------------------------------------------------------

!-----From Seonho Choi.  06/07/03
      SUBROUTINE INTERPOL_1DIM(N,XP,YP,EYP,X,Y,EY,IFLAG,ITYPE_ERR,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (DELTA=1.0E-3)
      REAL*8 XP,YP,EYP,X,Y,EY
      REAL*8 Y1,EY1,Y2,EY2,W1,W2,SW,A,B
      INTEGER N,ITYPE_ERR,IFLAG,JFLAG,IERR,I,J,K,L,NDEG,NPT_MIN
      DIMENSION XP(*),YP(*),EYP(*),NPT_MIN(20)
      LOGICAL L_LEFT,L_RIGHT
      DATA NPT_MIN/2,3,3,4,4,5,6,7,8,11*0/
!
! IFLAG = 1 : linear interpolation
! IFLAG = 2 : quadratic interpolation (a x b c)
! IFLAG = 3 : quadratic interpolation (a b x c)
! IFLAG = 4 : quadratic interpolation (average if IFLAG=2 and IFLAG=3)
! IFLAG = 5 : cubic interpolation
! IFLAG = 6 : 4th order
! IFLAG = 7 : 5th order
! IFLAG = 8 : 6th order
! IFLAG = 9 : 7th order
!
! IFLAG > 10 : does the extrapolation, too
!
! ITYPE_ERR = -1 : no error calculation
! ITYPE_ERR =  0 : all the errors are independent (statistical error)
! ITYPE_ERR =  1 : all the errors are fully correlated (systematic error)
!
! IERR = 0  : normal exit (interpolation)
! IERR = 1  : normal exit (left extrapolation)
! IERR = 2  : normal exit (right extrapolation)
! IERR = -1 : left out of bound
! IERR = -2 : right out of bound
! IERR = -3 : only 1 data point
!
!
! Check if we have sufficient number of points for a given method
!
      JFLAG = MOD(IFLAG,10)
      IF (N.LT.NPT_MIN(JFLAG)) THEN
         IF (N.LT.2) THEN ! nothing is possible, return
            IERR = -3
            RETURN
         ELSE IF (N.EQ.2) THEN ! only linear interpolation is possible
            JFLAG = 1
         ELSE IF (N.EQ.3.AND.JFLAG.GT.2) THEN  ! only quadratic interpolation is possible
            JFLAG = 2
         ELSE IF (N.EQ.4.AND.JFLAG.GT.5) THEN
            JFLAG = 5
         ELSE IF (N.GE.5.AND.JFLAG.GT.6) THEN
            JFLAG = N + 1
         ENDIF
      ENDIF

      L_LEFT = .FALSE.
      L_RIGHT = .FALSE.
      IF ((X-XP(1))/X.LT.-DELTA) L_LEFT = .TRUE.
      IF ((X-XP(N))/X.GT.DELTA) L_RIGHT = .TRUE.
!
! No extrapolation when JFLAG < 10
!
      IF (IFLAG.LT.10) THEN
         IF (L_LEFT) THEN
            IERR = -1
            Y = YP(1)
            EY = EYP(1)
            RETURN
         ELSE IF (L_RIGHT) THEN
            IERR = -2
            Y = YP(N)
            EY = EYP(N)
            RETURN
         ENDIF
      ENDIF
!
! When one of XP(I) is close enough to X, no need to go further
!
      DO I = 1,N
         IF (ABS((X-XP(I))/X).LT.DELTA) THEN
            Y = YP(I)
            EY = EYP(I)
            IERR = 0
            RETURN
         ENDIF
      ENDDO

!
! Find the corresponding interval
!      
      J = 0
      IF (L_LEFT) THEN
         J = 0
      ELSE IF (L_RIGHT) THEN
         J = N
      ELSE
         J = -1
         DO I = 1,N-1
            IF (X.GT.XP(I).AND.X.LT.XP(I+1)) THEN
               J = I
            ENDIF
         ENDDO
      ENDIF
!
! If everything is OK, you should not have J = -1, but just in case
!
      IF (J.EQ.-1) THEN
         IERR = -4
         RETURN
      ENDIF
!
! At boundaries, force other methods for quadratic interpolation
!
      IF (J.EQ.0.OR.J.EQ.1) THEN
         IF (JFLAG.EQ.3.OR.JFLAG.EQ.4) JFLAG = 2
      ELSE IF (J.EQ.N-1.OR.J.EQ.N) THEN
         IF (JFLAG.EQ.2.OR.JFLAG.EQ.4) JFLAG = 3
      ENDIF

      IF (JFLAG.EQ.1) THEN
         K = J
         L = K + 1
      ELSE IF (JFLAG.EQ.2) THEN
         K = J
         L = K + 2
      ELSE IF (JFLAG.EQ.3) THEN
         K = J - 1
         L = J + 1
      ELSE IF (JFLAG.EQ.4.OR.JFLAG.EQ.5) THEN
         K = J - 1
         L = J + 2
      ELSE IF (JFLAG.EQ.6) THEN
         K = J - 1
         L = J + 3
      ELSE IF (JFLAG.EQ.7) THEN
         K = J - 2
         L = J + 3
      ELSE IF (JFLAG.EQ.8) THEN
         K = J - 2
         L = J + 4
      ELSE IF (JFLAG.EQ.9) THEN
         K = J - 3
         L = J + 4
      ENDIF
!
! Consider boundary conditions
!
      IF (K.LT.1) THEN
         L = L + (1 - K)
         K = 1
      ENDIF
      IF (L.GT.N) THEN
         K = K - (L - N)
         L = N
      ENDIF

      IF (JFLAG.EQ.1) THEN ! linear interpolation
         A = (XP(K+1)-X)/(XP(K+1)-XP(K))
         B = (X-XP(K))/(XP(K+1)-XP(K))
         Y = A*YP(K)+B*YP(K+1)
         IF (ITYPE_ERR.EQ.-1) THEN
            EY = -1.0
         ELSE IF (ITYPE_ERR.EQ.0) THEN
            EY = (A*EYP(K))**2+(B*EYP(K+1))**2
            EY = SQRT(EY)
         ELSE IF (ITYPE_ERR.EQ.1) THEN
            EY = A*EYP(K)+B*EYP(K+1)
         ENDIF
      ELSE IF (JFLAG.EQ.2.OR.JFLAG.EQ.3) THEN
         CALL INTERPOL_PARA(XP(K),YP(K),EYP(K),X,Y,EY,ITYPE_ERR)
      ELSE IF (JFLAG.EQ.4) THEN
         CALL INTERPOL_PARA(XP(K),YP(K),EYP(K),X,Y1,EY1,ITYPE_ERR)
         CALL INTERPOL_PARA(XP(K+1),YP(K+1),EYP(K+1),X,Y2,EY2,ITYPE_ERR)
         IF (EY1.GT.0.0) THEN
            W1 = 1.0/EY1**2
         ELSE
            W1 = 1.0
         ENDIF
         IF (EY2.GT.0.0) THEN
            W2 = 1.0/EY2**2
         ELSE
            W2 = 0.0
         ENDIF
         SW = 1.0/(W1 + W2)
         Y = (W1*Y1+W2*Y2)*SW
         EY = SQRT(SW)
      ELSE IF (JFLAG.EQ.5) THEN
         CALL INTERPOL_CUBIC(XP(K),YP(K),EYP(K),X,Y,EY,ITYPE_ERR)
      ELSE
         NDEG = JFLAG - 2
         CALL INTERPOL_NDEG(NDEG,XP(K),YP(K),EYP(K),X,Y,EY,ITYPE_ERR)
      ENDIF
      IERR = 0
      RETURN
      END


      SUBROUTINE INTERPOL_CUBIC(XP,YP,EYP,X,Y,EY,ITYPE_ERR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XP,YP,EYP,X,Y,EY
      DIMENSION XP(*),YP(*),EYP(*),C(4)
      INTEGER ITYPE_ERR

      D1 = X - XP(1)
      D2 = X - XP(2)
      D3 = X - XP(3)
      D4 = X - XP(4)

      D12 = XP(1) - XP(2)
      D13 = XP(1) - XP(3)
      D14 = XP(1) - XP(4)
      D23 = XP(2) - XP(3)
      D24 = XP(2) - XP(4)
      D34 = XP(3) - XP(4)

      C(1) = (D2*D3*D4)/(D12*D13*D14)
      C(2) = (D1*D3*D4)/(-D12*D23*D24)
      C(3) = (D1*D2*D4)/(D13*D23*D34)
      C(4) = (D1*D2*D3)/(-D14*D24*D34)

      Y = 0.0
      EY = 0.0
      DO I = 1,4
         Y = Y + C(I)*YP(I)
         IF (ITYPE_ERR.EQ.-1) THEN
            EY = -1.0
         ELSE IF (ITYPE_ERR.EQ.0) THEN
            EY = EY + (C(I)*EYP(I))**2
         ELSE IF (ITYPE_ERR.EQ.1) THEN
            EY = EY + C(I)*EYP(I)
         ENDIF
      ENDDO
      IF (ITYPE_ERR.EQ.0) THEN
         EY = SQRT(EY)
      ENDIF
      RETURN
      END

      SUBROUTINE INTERPOL_NDEG(NDEG,XP,YP,EYP,X,Y,EY,ITYPE_ERR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XP,YP,EYP,X,Y,EY
      INTEGER NDEG,ITYPE_ERR
      DIMENSION XP(*),YP(*),EYP(*),D(100),DIJ(100,100),C(100)

      NP = NDEG + 1

      DO I = 1,NP
         D(I) = X - XP(I)
         DO J = I+1,NP
            DIJ(I,J) = XP(I) - XP(J)
            DIJ(J,I) = -DIJ(I,J)
         ENDDO
      ENDDO

      DO I = 1,NP
         C(I) = 1.0
         DO J = 1,NP
            IF (J.NE.I) C(I) = C(I)*D(J)/DIJ(I,J)
         ENDDO
      ENDDO

      Y = 0.0
      EY = 0.0
      DO I = 1,NP
         Y = Y + C(I)*YP(I)
         IF (ITYPE_ERR.EQ.-1) THEN
            EY = -1.0
         ELSE IF (ITYPE_ERR.EQ.0) THEN
            EY = EY + (C(I)*EYP(I))**2
         ELSE IF (ITYPE_ERR.EQ.1) THEN
            EY = EY + C(I)*EYP(I)
         ENDIF
      ENDDO
      IF (ITYPE_ERR.EQ.0) THEN
         EY = SQRT(EY)
      ENDIF
      RETURN
      END
