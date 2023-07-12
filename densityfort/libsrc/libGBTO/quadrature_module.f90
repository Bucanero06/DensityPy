! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
!> Quadrature module
!> =================
!> Contains routines for integration of a b-spline with an arbitrary user-defined function.
!> The function is defined using the bound_user_function class.
!> See function_integration module for details on how this is used.
module quadrature_module
 use precisn
 use general_quadrature, only: bound_user_function
 use utils, only: xermsg
 use const, only: cfp_dummy, wp_dummy, ep_dummy

 private

 interface cfp_bsgq8
    module procedure wp_bsgq8, ep_bsgq8
 end interface

 interface cfp_bssgq8
    module procedure wp_bssgq8, ep_bssgq8
 end interface

 interface cfp_polint
    module procedure wp_polint, ep_polint
 end interface

 interface cfp_arth
    module procedure wp_arth, ep_arth
 end interface

 interface cfp_trapzd
    module procedure wp_trapzd, ep_trapzd
 end interface

 public cfp_bfqad, cfp_bsqad, cfp_bfqro

 integer, parameter :: npar_arth=16,npar2_arth=8

contains

!>Purpose:
!>============
!> \verbatim
!>         Compute the integral of a product of a function and a
!>         derivative of a K-th order B-spline:
!> \endverbatim
!>         \f[
!>          \int_{x_{1}}^{x_{2}} dr B(r)f(r) 
!>         \f]
!> \verbatim
!>AUTHOR  Amos, D. E., (SNLA)
!>DESCRIPTION
!>
!>     Abstract    **** a double precision routine ****
!>
!>         cfp_bfqad computes the integral on (X1,X2) of a product of a
!>         function F and the ID-th derivative of a K-th order B-spline,
!>         using the B-representation (T,BCOEF,N,K).  (X1,X2) must be a
!>         subinterval of T(K) .LE. X .LE. T(N+1).  An integration rou-
!>         tine, DBSGQ8 (a modification of GAUS8), integrates the product
!>         on subintervals of (X1,X2) formed by included (distinct) knots
!>
!>         The maximum number of significant digits obtainable in
!>         DBSQAD is the smaller of 18 and the number of digits
!>         carried in double precision arithmetic.
!>
!>     Description of Arguments
!>         Input      F,T,BCOEF,X1,X2,TOL are double precision
!> \endverbatim
!>
!> \param[in] F     
!> \verbatim
!>            Function of one argument for the integrand BF(X)=F(X)*BVALU(T,BCOEF,N,K,ID,X,INBV,WORK)
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>            Knot array of length N+K
!> \endverbatim
!>
!> \param[in] BCOEF
!> \verbatim
!>            Coefficient array of length N
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>            Length of coefficient array
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>            Order of B-spline, K .GE. 1
!> \endverbatim
!>
!> \param[in] ID
!> \verbatim
!>            Order of the spline derivative, 0 .LE. ID .LE. K-1. ID=0 gives the spline function
!> \endverbatim
!>
!> \param[in] X1, X2
!> \verbatim  
!>            End points of quadrature interval in  T(K) .LE. X .LE. T(N+1)
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>            Desired accuracy for the quadrature, suggest 10.*DTOL .LT. TOL .LE. .1 where DTOL is the maximum  of 1.0D-18 and double precision unit roundoff for the machine = F1MACH(4)
!> \endverbatim
!>
!> \param[out] QUAD
!>             \f$ \int_{x_{1}}^{x_{2}} dr B(r)f(r) \f$
!> \verbatim
!>             Integral of BF(X) on (X1,X2)
!> \endverbatim
!>
!> \param[out] IERR
!> \verbatim
!>             A status code
!>                    IERR=1  normal return
!>                         2  some quadrature on (X1,X2) does not meet the requested tolerance.
!> \endverbatim
!>
!> \param[in,out] WORK
!> \verbatim
!>             Work vector of length 3*K
!> \endverbatim
!> \verbatim
!>
!>     Error Conditions
!>         Improper input is a fatal error
!>         Some quadrature fails to meet the requested tolerance
!>
!>***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
!>                 B-splines, Report SAND79-1825, Sandia Laboratories,
!>                 December 1979.
!>***ROUTINES CALLED  F1MACH, DBSGQ8, INTRV, XERMSG
!> \endverbatim
!
      SUBROUTINE cfp_bfqad (F, T, BCOEF, N, K, ID, X1, X2, TOL, QUAD, IERR, WORK)
      use bspline_base
      
      IMPLICIT NONE
      INTEGER ID, IERR, IFLG, ILO, IL1, IL2, K, LEFT, MFLAG, N, NPK, NP1, INBV
      REAL(kind=cfp) A,AA,ANS,B,BB,Q,QUAD,TA,TB,TOL,WTOL, X1, X2, WORK(:)
      REAL(kind=cfp), INTENT(IN) :: BCOEF(:),T(:)
      class(bound_user_function) :: F
! ***FIRST EXECUTABLE STATEMENT  cfp_bfqad
      IERR = 1
      QUAD = 0.0_cfp
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(ID.LT.0 .OR. ID.GE.K) GO TO 110
      WTOL = F1MACH(4,cfp_dummy)
      if (cfp .eq. wp) then
         WTOL = MAX(WTOL,1.0E-18_cfp)
      else
         WTOL = MAX(WTOL,1.0E-35_cfp)
      endif
      IF (TOL.LT.WTOL .OR. TOL.GT.0.1_cfp) GO TO 30
      AA = MIN(X1,X2)
      BB = MAX(X1,X2)
      IF (AA.LT.T(K)) GO TO 20
      NP1 = N + 1
      IF (BB.GT.T(NP1)) GO TO 20
      IF (AA.EQ.BB) RETURN
      NPK = N + K
! 
      ILO = 1
      CALL INTRV(T, NPK, AA, ILO, IL1, MFLAG)
      CALL INTRV(T, NPK, BB, ILO, IL2, MFLAG)
      IF (IL2.GE.NP1) IL2 = N
      INBV = 1
      Q = 0.0_cfp
      DO 10 LEFT=IL1,IL2
        TA = T(LEFT)
        TB = T(LEFT+1)
        IF (TA.EQ.TB) GO TO 10
        A = MAX(AA,TA)
        B = MIN(BB,TB)
        CALL CFP_BSGQ8(F,T,BCOEF,N,K,ID,A,B,INBV,TOL,ANS,IFLG,WORK)
!        CALL wp_bsgq8(F,T,BCOEF,N,K,ID,A,B,INBV,TOL,ANS,IFLG,WORK)
        IF (IFLG.GT.1) IERR = 2
        Q = Q + ANS
   10 CONTINUE
      IF (X1.GT.X2) Q = -Q
      QUAD = Q
      RETURN
! 
! 
   20 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bfqad', 'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)', 2, 1)
      RETURN
   30 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bfqad', 'TOL IS LESS DTOL OR GREATER THAN 0.1', 2, 1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bfqad', 'K DOES NOT SATISFY K.GE.1', 2, 1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bfqad', 'N DOES NOT SATISFY N.GE.K', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bfqad', 'ID DOES NOT SATISFY 0.LE.ID.LT.K', 2, 1)
      RETURN
      END SUBROUTINE

!>Purpose:
!>============
!> \verbatim
!>***BEGIN PROLOGUE  wp_bsgq8
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to DBFQAD
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=wp) (BSGQ8-S, wp_bsgq8-D)
!>***AUTHOR  Jones, R. E., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract    **** A REAL(kind=wp) routine ****
!>
!>        wp_bsgq8, a modification of GAUS8, integrates the
!>        product of FUN(X) by the ID-th derivative of a spline
!>        BVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
!>
!>     Description of Arguments
!>
!>        INPUT-- FUN,XT,BC,A,B,ERR are REAL(kind=wp)
!>        FUN - Name of function which multiplies BVALU.
!>        XT  - Knot array for BVALU
!>        BC  - B-coefficient array for BVALU
!>        N   - Number of B-coefficients for BVALU
!>        KK  - Order of the spline, KK.GE.1
!>        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
!>        A   - Lower limit of integral
!>        B   - Upper limit of integral (may be less than A)
!>        INBV- Initialization parameter for BVALU
!>        ERR - Is a requested pseudorelative error tolerance.  Normally
!>              pick a value of ABS(ERR).LT.1D-3.  ANS will normally
!>              have no more error than ABS(ERR) times the integral of
!>              the absolute value of FUN(X)*BVALU(XT,BC,N,KK,X,ID,
!>              INBV,WORK).
!>
!>
!>        OUTPUT-- ERR,ANS,WORK are REAL(kind=wp)
!>        ERR - Will be an estimate of the absolute error in ANS if the
!>              input value of ERR was negative.  (ERR is unchanged if
!>              the input value of ERR was nonnegative.)  The estimated
!>              error is solely for information to the user and should
!>              not be used as a correction to the computed integral.
!>        ANS - Computed value of integral
!>        IERR- A status code
!>            --Normal Codes
!>               1 ANS most likely meets requested error tolerance,
!>                 or A=B.
!>              -1 A and B are too nearly equal to allow normal
!>                 integration.  ANS is set to zero.
!>            --Abnormal Code
!>               2 ANS probably does not meet requested error tolerance.
!>        WORK- Work vector of length 3*K for BVALU
!>
!>***SEE ALSO  DBFQAD
!>***ROUTINES CALLED  F1MACH, BVALU, I1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!>***END PROLOGUE  wp_bsgq8
!> \endverbatim
      SUBROUTINE wp_bsgq8 (FUN, XT, BC, N, KK, ID, A, B, INBV, ERR, ANS, IERR, WORK)
      use bspline_base

      IMPLICIT NONE
      INTEGER ID, IERR, INBV, K, KK, KML, KMX, L, LMN, LMX, LR(60), MXL, N, NBITS, NIB, NLMN, NLMX
      REAL(kind=wp) A, AA(60), AE, ANIB, ANS, AREA, B, C, CE, EE, EF, EPS, ERR, EST, GL, GLR, GR(60), HH(60), SQ2, TOL, VL(60)
      REAL(kind=wp), INTENT(IN) :: BC(:),XT(:)
      REAL(kind=wp) VR, W1, W2, W3, W4, X1, X2, X3, X4, X, H, WORK(:), G8
      class(bound_user_function) :: FUN
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
      DATA X1, X2, X3, X4/1.83434642495649805D-01,5.25532409916328986D-01,7.96666477413626740D-01,9.60289856497536232D-01/
      DATA W1, W2, W3, W4/3.62683783378361983D-01,3.13706645877887287D-01,2.22381034453374471D-01,1.01228536290376259D-01/
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
! Gauss rule transformed from the interval [-1;1]
      G8(X,H) = H * &
        ((W1*(FUN%eval(X-X1*H)*BVALU(XT,BC,N,KK,ID,X-X1*H,INBV,WORK) + FUN%eval(X+X1*H)*BVALU(XT,BC,N,KK,ID,X+X1*H,INBV,WORK))  &
        + W2*(FUN%eval(X-X2*H)*BVALU(XT,BC,N,KK,ID,X-X2*H,INBV,WORK) + FUN%eval(X+X2*H)*BVALU(XT,BC,N,KK,ID,X+X2*H,INBV,WORK))) &
       + (W3*(FUN%eval(X-X3*H)*BVALU(XT,BC,N,KK,ID,X-X3*H,INBV,WORK) + FUN%eval(X+X3*H)*BVALU(XT,BC,N,KK,ID,X+X3*H,INBV,WORK))  &
        + W4*(FUN%eval(X-X4*H)*BVALU(XT,BC,N,KK,ID,X-X4*H,INBV,WORK) + FUN%eval(X+X4*H)*BVALU(XT,BC,N,KK,ID,X+X4*H,INBV,WORK))))
! 
!      INITIALIZE
! 
! ***FIRST EXECUTABLE STATEMENT  wp_bsgq8
      K = I1MACH(14)
      ANIB = F1MACH(5,wp_dummy)*K/0.30102000D0
      NBITS = INT(ANIB)
      NLMX = MIN((NBITS*5)/8,60)
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A.EQ.B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B.EQ.0.0D0) GO TO 10
      IF (SIGN(1.0D0,B)*A.LE.0.0D0) GO TO 10
      C = ABS(1.0D0-A/B)
      IF (C.GT.0.1D0) GO TO 10
      IF (C.LE.0.0D0) GO TO 140
      ANIB = 0.5D0 - LOG(C)/0.69314718D0
      NIB = INT(ANIB)
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX.LT.1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0
      IF (ERR.EQ.0.0D0) TOL = SQRT(F1MACH(4,wp_dummy))
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5D0
      MXL = 0
! 
!      COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
! 
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE <= AE) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L) <= 0) THEN
         GO TO 60
      ELSE
         GO TO 80
      END IF
! 
!      CONSIDER THE LEFT HALF OF THIS LEVEL
! 
   50 IF (K.GT.KMX) LMX = KML
      IF (L.GE.LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
! 
!      PROCEED TO RIGHT HALF AT THIS LEVEL
! 
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
! 
!      RETURN ONE LEVEL
! 
   80 VR = GLR
   90 IF (L.LE.1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L) <= 0) THEN
  100    VL(L) = VL(L+1) + VR
         GO TO 70
      ELSE
  110    VR = VL(L+1) + VR
         GO TO 90
      END IF
! 
!       EXIT
! 
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
      PRINT *,'ANS=',VR
      !CALL XERMSG ('quadrature_module', 'wp_bsgq8', 'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 1) !This was here originally. We want to replace it by warning only:
      CALL XERMSG ('quadrature_module', 'wp_bsgq8', 'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 0) !issue a warning, but don't terminate the program
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('quadrature_module', 'wp_bsgq8', &
                   'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. ANS IS SET TO ZERO AND IERR TO -1.', 1, 0)
  140 CONTINUE
      IF (ERR.LT.0.0D0) ERR = CE
      RETURN
      END SUBROUTINE

!> \verbatim
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to cfp_bfqad
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=ep1) (BSGQ8-S, ep_bsgq8-D)
!>***AUTHOR  Jones, R. E., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract    **** A REAL(kind=ep1) routine ****
!>
!>        ep_bsgq8, a modification of GAUS8, integrates the
!>        product of FUN(X) by the ID-th derivative of a spline
!>        BVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
!>
!>     Description of Arguments
!>
!>        INPUT-- FUN,XT,BC,A,B,ERR are REAL(kind=ep1)
!>        FUN - Name of external function of one argument which
!>              multiplies BVALU.
!>        XT  - Knot array for BVALU
!>        BC  - B-coefficient array for BVALU
!>        N   - Number of B-coefficients for BVALU
!>        KK  - Order of the spline, KK.GE.1
!>        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
!>        A   - Lower limit of integral
!>        B   - Upper limit of integral (may be less than A)
!>        INBV- Initialization parameter for BVALU
!>        ERR - Is a requested pseudorelative error tolerance.  Normally
!>              pick a value of ABS(ERR).LT.1D-3.  ANS will normally
!>              have no more error than ABS(ERR) times the integral of
!>              the absolute value of FUN(X)*BVALU(XT,BC,N,KK,X,ID,
!>              INBV,WORK).
!>
!>
!>        OUTPUT-- ERR,ANS,WORK are REAL(kind=ep1)
!>        ERR - Will be an estimate of the absolute error in ANS if the
!>              input value of ERR was negative.  (ERR is unchanged if
!>              the input value of ERR was nonnegative.)  The estimated
!>              error is solely for information to the user and should
!>              not be used as a correction to the computed integral.
!>        ANS - Computed value of integral
!>        IERR- A status code
!>            --Normal Codes
!>               1 ANS most likely meets requested error tolerance,
!>                 or A=B.
!>              -1 A and B are too nearly equal to allow normal
!>                 integration.  ANS is set to zero.
!>            --Abnormal Code
!>               2 ANS probably does not meet requested error tolerance.
!>        WORK- Work vector of length 3*K for BVALU
!>
!>***SEE ALSO  DBFQAD
!>***ROUTINES CALLED  Q1MACH, BVALU, I1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!> \endverbatim
!
      SUBROUTINE ep_bsgq8 (FUN, XT, BC, N, KK, ID, A, B, INBV, ERR, ANS, IERR, WORK)
!
!
      use bspline_base
      IMPLICIT NONE
      INTEGER ID, IERR, INBV, K, KK, KML, KMX, L, LMN, LMX, LR, MXL, N, NBITS, NIB, NLMN, NLMX
      REAL(kind=ep1) :: XT(:),BC(:),A,B,ERR,ANS,WORK(:)
      REAL(kind=ep1) :: AA,AE,ANIB,AREA,C,CE,EE,EF,EPS,EST,GL,GLR,GR,HH,SQ2,TOL,VL,VR,W1, W2, W3, W4, X1, X2, X3, X4, X, H
      REAL(kind=ep1) :: G8
      class(bound_user_function) :: FUN
      DIMENSION AA(120), HH(120), LR(120), VL(120), GR(120)
! THE WEIGHTS AND ABSCISSAS NEED TO BE IN QUADRUPLE PRECISION AS WELL!!!
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
      DATA X1, X2, X3, X4/&
! quad precision values (Mathematica):
     &     0.18343464249564980493947614236018400_ep1, &
     &     0.52553240991632898581773904918924600_ep1,&
     &     0.79666647741362673959155393647583000_ep1, &
     &     0.96028985649753623168356086856947300_ep1/
! double precision values:
!     1     1.83434642495649805D-01,     5.25532409916328986D-01,
!     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/&
     &     0.36268378337836198296515044927719560_ep1, &
     &     0.31370664587788728733796220198660130_ep1,&
     &     0.22238103445337447054435599442624090_ep1,&
     &     0.10122853629037625915253135430996220_ep1/
! double precision values:
!     1     3.62683783378361983D-01,     3.13706645877887287D-01,
!     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA SQ2/1.41421356237309504880168872420969808_ep1/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H) = H * &
        ((W1*(FUN%eval(X-X1*H)*BVALU(XT,BC,N,KK,ID,X-X1*H,INBV,WORK) + FUN%eval(X+X1*H)*BVALU(XT,BC,N,KK,ID,X+X1*H,INBV,WORK))  &
        + W2*(FUN%eval(X-X2*H)*BVALU(XT,BC,N,KK,ID,X-X2*H,INBV,WORK) + FUN%eval(X+X2*H)*BVALU(XT,BC,N,KK,ID,X+X2*H,INBV,WORK))) &
        +(W3*(FUN%eval(X-X3*H)*BVALU(XT,BC,N,KK,ID,X-X3*H,INBV,WORK) + FUN%eval(X+X3*H)*BVALU(XT,BC,N,KK,ID,X+X3*H,INBV,WORK))  &
        + W4*(FUN%eval(X-X4*H)*BVALU(XT,BC,N,KK,ID,X-X4*H,INBV,WORK) + FUN%eval(X+X4*H)*BVALU(XT,BC,N,KK,ID,X+X4*H,INBV,WORK))))
!
!     INITIALIZE
!
!***FIRST EXECUTABLE STATEMENT  ep_bsgq8
      K = I1MACH(17)
      ANIB = F1MACH(5,ep_dummy)*K/0.30102000_ep1
      NBITS = INT(ANIB)
      NLMX = MIN((NBITS*5)/8,120)
      ANS = 0.0_ep1
      IERR = 1
      CE = 0.0_ep1
      IF (A.EQ.B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B.EQ.0.0_ep1) GO TO 10
      IF (SIGN(1.0_ep1,B)*A.LE.0.0_ep1) GO TO 10
      C = ABS(1.0_ep1-A/B)
      IF (C.GT.0.1_ep1) GO TO 10
      IF (C.LE.0.0_ep1) GO TO 140
      ANIB = 0.5_ep1 - LOG(C)/0.69314718_ep1
      NIB = INT(ANIB)
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX.LT.1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0_ep1**(5-NBITS))/2.0_ep1
      IF (ERR.EQ.0.0_ep1) TOL = SQRT(F1MACH(4,ep_dummy))
      EPS = TOL
      HH(1) = (B-A)/4.0_ep1
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0_ep1*HH(L),2.0_ep1*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5_ep1
      MXL = 0
!
!     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
!
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0_ep1*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE <= AE) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L) <= 0) THEN
         GO TO 60
      ELSE
         GO TO 80
      END IF
!
!     CONSIDER THE LEFT HALF OF THIS LEVEL
!
   50 IF (K.GT.KMX) LMX = KML
      IF (L.GE.LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5_ep1
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5_ep1
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
!
!     PROCEED TO RIGHT HALF AT THIS LEVEL
!
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0_ep1*HH(L)
      GO TO 20
!
!     RETURN ONE LEVEL
!
   80 VR = GLR
   90 IF (L.LE.1) GO TO 120
      L = L - 1
      EPS = EPS*2.0_ep1
      EF = EF*SQ2
      IF (LR(L) <= 0) THEN
  100    VL(L) = VL(L+1) + VR
         GO TO 70
      ELSE
  110    VR = VL(L+1) + VR
         GO TO 90
      END IF
!
!      EXIT
!
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0_ep1*TOL*AREA)) GO TO 140
      IERR = 2
      PRINT *,'ANS=',VR
      PRINT *,'MXL,EST',MXL,ABS(CE),2.0_ep1*TOL*AREA
      CALL XERMSG ('quadrature_module', 'ep_bsgq8', 'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 0)
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('quadrature_module', 'ep_bsgq8', &
                   'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. ANS IS SET TO ZERO AND IERR TO -1.', 1, 0)
  140 CONTINUE
      IF (ERR.LT.0.0_ep1) ERR = CE
      RETURN
      END SUBROUTINE

!>Purpose:
!>============
!> \verbatim
!>***BEGIN PROLOGUE  cfp_bsqad
!>***PURPOSE  Compute the integral of a product of a function and a square of the
!>            derivative of a K-th order B-spline.
!>***LIBRARY   SLATEC
!>***CATEGORY  H2A2A1, E3, K6
!>***TYPE      REAL(kind=cfp) (BFQAD-S, cfp_bsqad-D)
!>***KEYWORDS  INTEGRAL OF B-SPLINE, QUADRATURE
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract    **** a double precision routine ****
!>
!>         cfp_bsqad computes the integral on (X1,X2) of a product of a
!>         function F and the ID-th derivative of a K-th order B-spline,
!>         using the B-representation (T,BCOEF,N,K).  (X1,X2) must be a
!>         subinterval of T(K) .LE. X .LE. T(N+1).  An integration rou-
!>         tine, DBSGQ8 (a modification of GAUS8), integrates the product
!>         on subintervals of (X1,X2) formed by included (distinct) knots
!>
!>         The maximum number of significant digits obtainable in
!>         DBSQAD is the smaller of 18 and the number of digits
!>         carried in double precision arithmetic.
!>
!>     Description of Arguments
!>         Input      F,T,BCOEF,X1,X2,TOL are double precision
!>           F      - function of one argument for the integrand BF(X)=F(X)*BVALU(T,BCOEF,N,K,ID,X,INBV,WORK)
!>           T      - knot array of length N+K
!>           BCOEF  - coefficient array of length N
!>           N      - length of coefficient array
!>           K      - order of B-spline, K .GE. 1
!>           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1
!>                    ID=0 gives the spline function
!>           X1,X2  - end points of quadrature interval in
!>                    T(K) .LE. X .LE. T(N+1)
!>           TOL    - desired accuracy for the quadrature, suggest
!>                    10.*DTOL .LT. TOL .LE. .1 where DTOL is the maximum
!>                    of 1.0D-18 and double precision unit roundoff for
!>                    the machine = F1MACH(4)
!>
!>         Output     QUAD,WORK are double precision
!>           QUAD   - integral of BF(X) on (X1,X2)
!>           IERR   - a status code
!>                    IERR=1  normal return
!>                         2  some quadrature on (X1,X2) does not meet
!>                            the requested tolerance.
!>           WORK   - work vector of length 3*K
!>
!>     Error Conditions
!>         Improper input is a fatal error
!>         Some quadrature fails to meet the requested tolerance
!>
!>***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
!>                 B-splines, Report SAND79-1825, Sandia Laboratories,
!>                 December 1979.
!>***ROUTINES CALLED  F1MACH, DBSGQ8, INTRV, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!>***END PROLOGUE  cfp_bsqad
!>
!> \endverbatim
      SUBROUTINE cfp_bsqad (F, T, BCOEF, N, K, ID, X1, X2, TOL, QUAD, IERR, WORK)
      use bspline_base
      IMPLICIT NONE

      INTEGER ID, IERR, IFLG, ILO, IL1, IL2, K, LEFT, MFLAG, N, NPK, NP1, INBV
      REAL(kind=cfp) A,AA,ANS,B,BB,Q,QUAD,TA,TB,TOL,WTOL, X1, X2, WORK(:)
      REAL(kind=cfp), INTENT(IN) :: BCOEF(:),T(:)
      class(bound_user_function) :: F
! ***FIRST EXECUTABLE STATEMENT  cfp_bsqad
      IERR = 1
      QUAD = 0.0_cfp
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(ID.LT.0 .OR. ID.GE.K) GO TO 110
      WTOL = F1MACH(4,cfp_dummy)
      if (cfp .eq. wp) then
         WTOL = MAX(WTOL,1.E-18_cfp)
      else
         WTOL = MAX(WTOL,1.E-35_cfp)
      endif
      IF (TOL.LT.WTOL .OR. TOL.GT.0.1_cfp) GO TO 30
      AA = MIN(X1,X2)
      BB = MAX(X1,X2)
      IF (AA.LT.T(K)) GO TO 20
      NP1 = N + 1
!ZM: relaxed the test from BB.GT.T(NP1) since this is not good enough in cases
!there are negligible differences between BB and T(NP1).
      IF (BB-10*WTOL > T(NP1)) GO TO 20
      IF (AA.EQ.BB) RETURN
      NPK = N + K
! 
      ILO = 1
      CALL INTRV(T, NPK, AA, ILO, IL1, MFLAG)
      CALL INTRV(T, NPK, BB, ILO, IL2, MFLAG)
      IF (IL2.GE.NP1) IL2 = N
      INBV = 1
      Q = 0.0_cfp
      DO 10 LEFT=IL1,IL2
        TA = T(LEFT)
        TB = T(LEFT+1)
        IF (TA.EQ.TB) GO TO 10
        A = MAX(AA,TA)
        B = MIN(BB,TB)
        CALL CFP_BSSGQ8(F,T,BCOEF,N,K,ID,A,B,INBV,TOL,ANS,IFLG,WORK)
        !CALL wp_bssgq8(F,T,BCOEF,N,K,ID,A,B,INBV,TOL,ANS,IFLG,WORK)
        IF (IFLG.GT.1) IERR = 2
        Q = Q + ANS
   10 CONTINUE
      IF (X1.GT.X2) Q = -Q
      QUAD = Q
      RETURN
! 
! 
   20 CONTINUE
      print *,AA,T(K)
      print *,BB,T(NP1)
      CALL XERMSG ('quadrature_module', 'cfp_bsqad', 'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)', 2, 1)
      RETURN
   30 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bsqad', 'TOL IS LESS DTOL OR GREATER THAN 0.1', 2, 1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bsqad', 'K DOES NOT SATISFY K.GE.1', 2, 1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bsqad', 'N DOES NOT SATISFY N.GE.K', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('quadrature_module', 'cfp_bsqad', 'ID DOES NOT SATISFY 0.LE.ID.LT.K', 2, 1)
      RETURN
      END SUBROUTINE

!>Purpose:
!>============
!> \verbatim
!>***BEGIN PROLOGUE  wp_bssgq8
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to DBFQAD
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=wp) (BSGQ8-S, wp_bssgq8-D)
!>***AUTHOR  Jones, R. E., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract    **** A REAL(kind=wp) routine ****
!>
!>        wp_bssgq8, a modification of GAUS8, integrates the
!>        product of FUN(X) by the square of the ID-th derivative of a spline
!>        BVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
!>
!>     Description of Arguments
!>
!>        INPUT-- FUN,XT,BC,A,B,ERR are REAL(kind=wp)
!>        FUN - Name of function which multiplies BVALU^2.
!>        XT  - Knot array for BVALU
!>        BC  - B-coefficient array for BVALU
!>        N   - Number of B-coefficients for BVALU
!>        KK  - Order of the spline, KK.GE.1
!>        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
!>        A   - Lower limit of integral
!>        B   - Upper limit of integral (may be less than A)
!>        INBV- Initialization parameter for BVALU
!>        ERR - Is a requested pseudorelative error tolerance.  Normally
!>              pick a value of ABS(ERR).LT.1D-3.  ANS will normally
!>              have no more error than ABS(ERR) times the integral of
!>              the absolute value of FUN(X)*BVALU(XT,BC,N,KK,X,ID,
!>              INBV,WORK).
!>
!>
!>        OUTPUT-- ERR,ANS,WORK are REAL(kind=wp)
!>        ERR - Will be an estimate of the absolute error in ANS if the
!>              input value of ERR was negative.  (ERR is unchanged if
!>              the input value of ERR was nonnegative.)  The estimated
!>              error is solely for information to the user and should
!>              not be used as a correction to the computed integral.
!>        ANS - Computed value of integral
!>        IERR- A status code
!>            --Normal Codes
!>               1 ANS most likely meets requested error tolerance,
!>                 or A=B.
!>              -1 A and B are too nearly equal to allow normal
!>                 integration.  ANS is set to zero.
!>            --Abnormal Code
!>               2 ANS probably does not meet requested error tolerance.
!>        WORK- Work vector of length 3*K for BVALU
!>
!>***SEE ALSO  DBFQAD
!>***ROUTINES CALLED  F1MACH, BVALU, I1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!>***END PROLOGUE  wp_bssgq8
!> \endverbatim
      SUBROUTINE wp_bssgq8 (FUN, XT, BC, N, KK, ID, A, B, INBV, ERR, ANS, IERR, WORK)
      use bspline_base
      IMPLICIT NONE

      INTEGER ID, IERR, INBV, K, KK, KML, KMX, L, LMN, LMX, LR(60), MXL, N, NBITS, NIB, NLMN, NLMX
      REAL(kind=wp) A, AA(60), AE, ANIB, ANS, AREA, B, C, CE, EE, EF, EPS, ERR, EST, GL, GLR, GR(60), HH(60), SQ2, TOL, VL(60)
      REAL(kind=wp), INTENT(IN) :: BC(:),XT(:)
      REAL(kind=wp) VR, W1, W2, W3, W4, X1, X2, X3, X4, X, H, WORK(:), G8
      class(bound_user_function) :: FUN
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
      DATA X1, X2, X3, X4/1.83434642495649805D-01,5.25532409916328986D-01,7.96666477413626740D-01,9.60289856497536232D-01/
      DATA W1, W2, W3, W4/3.62683783378361983D-01,3.13706645877887287D-01,2.22381034453374471D-01,1.01228536290376259D-01/
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
! Gauss rule transformed from the interval [-1;1]
      G8(X,H) = H * &
        ((W1*(FUN%eval(X-X1*H)*BVALU(XT,BC,N,KK,ID,X-X1*H,INBV,WORK)**2.0_wp  &
            + FUN%eval(X+X1*H)*BVALU(XT,BC,N,KK,ID,X+X1*H,INBV,WORK)**2.0_wp) &
        + W2*(FUN%eval(X-X2*H)*BVALU(XT,BC,N,KK,ID,X-X2*H,INBV,WORK)**2.0_wp  &
            + FUN%eval(X+X2*H)*BVALU(XT,BC,N,KK,ID,X+X2*H,INBV,WORK)**2.0_wp))&
        +(W3*(FUN%eval(X-X3*H)*BVALU(XT,BC,N,KK,ID,X-X3*H,INBV,WORK)**2.0_wp  &
            + FUN%eval(X+X3*H)*BVALU(XT,BC,N,KK,ID,X+X3*H,INBV,WORK)**2.0_wp) &
        + W4*(FUN%eval(X-X4*H)*BVALU(XT,BC,N,KK,ID,X-X4*H,INBV,WORK)**2.0_wp  &
            + FUN%eval(X+X4*H)*BVALU(XT,BC,N,KK,ID,X+X4*H,INBV,WORK)**2.0_wp)))
! 
!      INITIALIZE
! 
! ***FIRST EXECUTABLE STATEMENT  wp_bssgq8
      K = I1MACH(14)
      ANIB = F1MACH(5,wp_dummy)*K/0.30102000D0
      NBITS = INT(ANIB)
      NLMX = MIN((NBITS*5)/8,60)
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A.EQ.B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B.EQ.0.0D0) GO TO 10
      IF (SIGN(1.0D0,B)*A.LE.0.0D0) GO TO 10
      C = ABS(1.0D0-A/B)
      IF (C.GT.0.1D0) GO TO 10
      IF (C.LE.0.0D0) GO TO 140
      ANIB = 0.5D0 - LOG(C)/0.69314718D0
      NIB = INT(ANIB)
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX.LT.1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0
      IF (ERR.EQ.0.0D0) TOL = SQRT(F1MACH(4,wp_dummy))
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5D0
      MXL = 0
!>
!>     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
!>
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE <= AE) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L) <= 0) THEN
         GO TO 60
      ELSE
         GO TO 80
      END IF
! 
!      CONSIDER THE LEFT HALF OF THIS LEVEL
! 
   50 IF (K.GT.KMX) LMX = KML
      IF (L.GE.LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
! 
!      PROCEED TO RIGHT HALF AT THIS LEVEL
! 
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
! 
!      RETURN ONE LEVEL
! 
   80 VR = GLR
   90 IF (L.LE.1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L) <= 0) THEN
  100    VL(L) = VL(L+1) + VR
         GO TO 70
      ELSE
  110    VR = VL(L+1) + VR
         GO TO 90
      END IF
! 
!       EXIT
! 
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
      PRINT *,'ANS=',VR
      CALL XERMSG ('quadrature_module', 'wp_bssgq8', 'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 0)
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('quadrature_module', 'wp_bssgq8', &
                   'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. ANS IS SET TO ZERO AND IERR TO -1.', 1, 0)
  140 CONTINUE
      IF (ERR.LT.0.0D0) ERR = CE
      RETURN
      END SUBROUTINE
!> \verbatim
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to DBFQAD
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=ep1) (BSGQ8-S, ep_bssgq8-D)
!>***AUTHOR  Jones, R. E., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract    **** A REAL(kind=ep1) routine ****
!>
!>        ep_bssgq8, a modification of GAUS8, integrates the
!>        product of FUN(X) and the square of ID-th derivative of a spline
!>        QBVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
!>
!>     Description of Arguments
!>
!>        INPUT-- FUN,XT,BC,A,B,ERR are REAL(kind=ep1)
!>        FUN - Name of external function of one argument which
!>              multiplies QBVALU.
!>        XT  - Knot array for QBVALU
!>        BC  - B-coefficient array for QBVALU
!>        N   - Number of B-coefficients for QBVALU
!>        KK  - Order of the spline, KK.GE.1
!>        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
!>        A   - Lower limit of integral
!>        B   - Upper limit of integral (may be less than A)
!>        INBV- Initialization parameter for QBVALU
!>        ERR - Is a requested pseudorelative error tolerance.  Normally
!>              pick a value of ABS(ERR).LT.1D-3.  ANS will normally
!>              have no more error than ABS(ERR) times the integral of
!>              the absolute value of FUN(X)*QBVALU(XT,BC,N,KK,X,ID,
!>              INBV,WORK).
!>
!>
!>        OUTPUT-- ERR,ANS,WORK are REAL(kind=ep1)
!>        ERR - Will be an estimate of the absolute error in ANS if the
!>              input value of ERR was negative.  (ERR is unchanged if
!>              the input value of ERR was nonnegative.)  The estimated
!>              error is solely for information to the user and should
!>              not be used as a correction to the computed integral.
!>        ANS - Computed value of integral
!>        IERR- A status code
!>            --Normal Codes
!>               1 ANS most likely meets requested error tolerance,
!>                 or A=B.
!>              -1 A and B are too nearly equal to allow normal
!>                 integration.  ANS is set to zero.
!>            --Abnormal Code
!>               2 ANS probably does not meet requested error tolerance.
!>        WORK- Work vector of length 3*K for QBVALU
!>
!>***SEE ALSO  DBFQAD
!>***ROUTINES CALLED  Q1MACH, QBVALU, I1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!> \endverbatim
!
      SUBROUTINE ep_bssgq8 (FUN, XT, BC, N, KK, ID, A, B, INBV, ERR, ANS, IERR, WORK)
!
!
      use precisn
      use bspline_base
      IMPLICIT NONE
      INTEGER ID, IERR, INBV, K, KK, KML, KMX, L, LMN, LMX, LR, MXL, N, NBITS, NIB, NLMN, NLMX
      REAL(kind=ep1) A, AA, AE, ANIB, ANS, AREA, B, BC(:), C, CE, EE, EF, EPS, ERR, EST, GL, GLR, GR, HH, SQ2, TOL, VL, VR
      REAL(kind=ep1) :: WORK(:), W1, W2, W3, W4, XT(:), X1, X2, X3, X4, X, H, G8
      class(bound_user_function) :: FUN
      DIMENSION AA(120), HH(120), LR(120), VL(120), GR(120)
! THE WEIGHTS AND ABSCISSAS NEED TO BE IN QUADRUPLE PRECISION AS WELL!!!
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
      DATA X1, X2, X3, X4/ &
! quad precision values (Mathematica):
     &     0.18343464249564980493947614236018400_ep1, &
     &     0.52553240991632898581773904918924600_ep1, &
     &     0.79666647741362673959155393647583000_ep1, &
     &     0.96028985649753623168356086856947300_ep1/
! double precision values:
!     1     1.83434642495649805D-01,     5.25532409916328986D-01,
!     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/ &
     &     0.36268378337836198296515044927719560_ep1, &
     &     0.31370664587788728733796220198660130_ep1, &
     &     0.22238103445337447054435599442624090_ep1, &
     &     0.10122853629037625915253135430996220_ep1/
! double precision values:
!     1     3.62683783378361983D-01,     3.13706645877887287D-01,
!     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA SQ2/1.414213560_ep1/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)= H * &
         (W1*(FUN%eval(X-X1*H)*(BVALU(XT,BC,N,KK,ID,X-X1*H,INBV,WORK))**2.0_ep1  &
            + FUN%eval(X+X1*H)*(BVALU(XT,BC,N,KK,ID,X+X1*H,INBV,WORK))**2.0_ep1) &
        + W2*(FUN%eval(X-X2*H)*(BVALU(XT,BC,N,KK,ID,X-X2*H,INBV,WORK))**2.0_ep1  &
            + FUN%eval(X+X2*H)*(BVALU(XT,BC,N,KK,ID,X+X2*H,INBV,WORK))**2.0_ep1) &
        + W3*(FUN%eval(X-X3*H)*(BVALU(XT,BC,N,KK,ID,X-X3*H,INBV,WORK))**2.0_ep1  &
            + FUN%eval(X+X3*H)*(BVALU(XT,BC,N,KK,ID,X+X3*H,INBV,WORK))**2.0_ep1) &
        + W4*(FUN%eval(X-X4*H)*(BVALU(XT,BC,N,KK,ID,X-X4*H,INBV,WORK))**2.0_ep1  &
            + FUN%eval(X+X4*H)*(BVALU(XT,BC,N,KK,ID,X+X4*H,INBV,WORK))**2.0_ep1))
!
!     INITIALIZE
!
!***FIRST EXECUTABLE STATEMENT  ep_bssgq8
      K = I1MACH(17)
      ANIB = F1MACH(5,ep_dummy)*K/0.30102000_ep1
      NBITS = INT(ANIB)
      NLMX = MIN((NBITS*5)/8,120)
      ANS = 0.0_ep1
      IERR = 1
      CE = 0.0_ep1
      IF (A.EQ.B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B.EQ.0.0_ep1) GO TO 10
      IF (SIGN(1.0_ep1,B)*A.LE.0.0_ep1) GO TO 10
      C = ABS(1.0_ep1-A/B)
      IF (C.GT.0.1_ep1) GO TO 10
      IF (C.LE.0.0_ep1) GO TO 140
      ANIB = 0.5_ep1 - LOG(C)/0.69314718_ep1
      NIB = INT(ANIB)
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX.LT.1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0_ep1**(5-NBITS))/2.0_ep1
      IF (ERR.EQ.0.0_ep1) TOL = SQRT(F1MACH(4,ep_dummy))
      EPS = TOL
      HH(1) = (B-A)/4.0_ep1
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0_ep1*HH(L),2.0_ep1*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5_ep1
      MXL = 0
!
!     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
!
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0_ep1*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE <= AE) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L) <= 0) THEN
         GO TO 60
      ELSE
         GO TO 80
      END IF
!
!     CONSIDER THE LEFT HALF OF THIS LEVEL
!
   50 IF (K.GT.KMX) LMX = KML
      IF (L.GE.LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5_ep1
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5_ep1
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
!
!     PROCEED TO RIGHT HALF AT THIS LEVEL
!
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0_ep1*HH(L)
      GO TO 20
!
!     RETURN ONE LEVEL
!
   80 VR = GLR
   90 IF (L.LE.1) GO TO 120
      L = L - 1
      EPS = EPS*2.0_ep1
      EF = EF*SQ2
      IF (LR(L) <= 0) THEN
  100    VL(L) = VL(L+1) + VR
         GO TO 70
      ELSE
  110    VR = VL(L+1) + VR
         GO TO 90
      END IF
!
!      EXIT
!
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0_ep1*TOL*AREA)) GO TO 140
      IERR = 2
      PRINT *,'ANS=',VR
      CALL XERMSG ('quadrature_module', 'ep_bssgq8', 'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 1)
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('quadrature_module', 'ep_bssgq8', &
                   'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. ANS IS SET TO ZERO AND IERR TO -1.', 1, 0)
  140 CONTINUE
      IF (ERR.LT.0.0_ep1) ERR = CE
      RETURN
      END SUBROUTINE

      SUBROUTINE cfp_bfqro (F, T, BCOEF, N, KK, ID, X1, X2, TOL, QUAD, IERR, WORK)
!
      IMPLICIT NONE
      INTEGER :: ID, IERR, KK, N, INBV
      REAL(kind=cfp) :: QUAD, TOL, X1, X2
      REAL(kind=cfp), INTENT(IN) :: BCOEF(:), T(:)
      REAL(kind=cfp), INTENT(INOUT) :: WORK(:)
      class(bound_user_function) :: F
!
      INTEGER, PARAMETER :: JMAX=400,JMAXP=JMAX+1,K=5,KM=K-1 !20
      REAL(kind=cfp) :: EPS
      REAL(kind=cfp), DIMENSION(JMAXP) :: h,s
      REAL(kind=cfp) :: dqromb
      INTEGER :: j
         IERR = 0
         EPS = F1MACH(4,cfp_dummy)
         EPS = MAX(EPS,TOL)
         h(1)=1.0_cfp
         INBV = 1
         do j=1,JMAX
            call cfp_trapzd(F,X1,X2,s(j),j,T,BCOEF,N,KK,ID,WORK,INBV)
            if (j >= K) then
               call cfp_polint(h(j-KM:j),s(j-KM:j),0.0_cfp,QUAD,dqromb)
               if (abs(dqromb) <= EPS*abs(QUAD)) RETURN
            end if
            s(j+1)=s(j)
            h(j+1)=0.25_cfp*h(j)
         end do
         IERR = 1
         call xermsg ('quadrature_module', 'cfp_bfqro', 'too many steps.',ierr,1)
      END SUBROUTINE cfp_bfqro
   
      SUBROUTINE wp_polint(xa,ya,x,y,dy)
      IMPLICIT NONE
      REAL(kind=wp), DIMENSION(:), INTENT(IN) :: xa,ya
      REAL(kind=wp), INTENT(IN) :: x
      REAL(kind=wp), INTENT(OUT) :: y,dy
      INTEGER :: m,n,ns
      REAL(kind=wp), DIMENSION(size(xa)) :: c,d,den,ho
      INTEGER, DIMENSION(1) :: imin
         if (size(xa) .eq. size(ya)) then
            n = size(xa)
         else
            call xermsg ('quadrature_module', 'wp_polint', 'sizes of the input arrays xa, ya must be the same.',1,1)
         endif
         c=ya
         d=ya
         ho=xa-x
         !ns=iminloc(abs(x-xa))
         imin=minloc(abs(x-xa))
         ns=imin(1)
         y=ya(ns)
         ns=ns-1
         do m=1,n-1
            den(1:n-m)=ho(1:n-m)-ho(1+m:n)
            if (any(den(1:n-m) == 0.0_wp)) call xermsg ('quadrature_module', 'wp_polint', 'Calculation failure.',1,1)
            den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
            d(1:n-m)=ho(1+m:n)*den(1:n-m)
            c(1:n-m)=ho(1:n-m)*den(1:n-m)
            if (2*ns < n-m) then
               dy=c(ns+1)
            else
               dy=d(ns)
               ns=ns-1
            end if
            y=y+dy
         end do
      END SUBROUTINE wp_polint

      SUBROUTINE ep_polint(xa,ya,x,y,dy)
      IMPLICIT NONE
      REAL(kind=ep1), DIMENSION(:), INTENT(IN) :: xa,ya
      REAL(kind=ep1), INTENT(IN) :: x
      REAL(kind=ep1), INTENT(OUT) :: y,dy
      INTEGER :: m,n,ns
      REAL(kind=ep1), DIMENSION(size(xa)) :: c,d,den,ho
      INTEGER, DIMENSION(1) :: imin
         if (size(xa) .eq. size(ya)) then
            n = size(xa)
         else
            call xermsg ('quadrature_module', 'ep_polint', 'sizes of the input arrays xa, ya must be the same.',1,1)
         endif
         c=ya
         d=ya
         ho=xa-x
         !ns=iminloc(abs(x-xa))
         imin=minloc(abs(x-xa))
         ns=imin(1)
         y=ya(ns)
         ns=ns-1
         do m=1,n-1
            den(1:n-m)=ho(1:n-m)-ho(1+m:n)
            if (any(den(1:n-m) == 0.0_ep1)) call xermsg ('quadrature_module', 'ep_polint', 'Calculation failure.',1,1)
            den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
            d(1:n-m)=ho(1+m:n)*den(1:n-m)
            c(1:n-m)=ho(1:n-m)*den(1:n-m)
            if (2*ns < n-m) then
               dy=c(ns+1)
            else
               dy=d(ns)
               ns=ns-1
            end if
            y=y+dy
         end do
      END SUBROUTINE ep_polint
   
      SUBROUTINE wp_trapzd(F,a,b,s,n,T,BCOEF,NB,K,ID,WORK,INBV)
      use bspline_base
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ID, K, NB, INBV
      REAL(kind=wp), INTENT(INOUT) :: WORK(:)
      REAL(kind=wp), INTENT(IN) :: BCOEF(:), T(:)
      class(bound_user_function) :: F
      REAL(kind=wp), INTENT(IN) :: a,b
      REAL(kind=wp), INTENT(INOUT) :: s
      INTEGER, INTENT(IN) :: n
      REAL(kind=wp) :: del,fsum,x(max(1,2**(n-2)))
      INTEGER :: it, j
         if (n == 1) then
            s = 0.5_wp*(b-a)*(F%eval(a)*BVALU(T,BCOEF,NB,K,ID,a,INBV,WORK) + F%eval(b)*BVALU(T,BCOEF,NB,K,ID,b,INBV,WORK))
            !s=0.5_wp*(b-a)*sum(func( (/ a,b /) ))
         else
            it=2**(n-2)
            del=(b-a)/it
            x = cfp_arth(a+0.5_wp*del,del,it)
            fsum = 0.0_wp
            do j=1,n
               fsum = fsum + F%eval(x(j))*BVALU(T,BCOEF,NB,K,ID,x(j),INBV,WORK)
            enddo
            !fsum=sum(func(cfp_arth(a+0.5_dp*del,del,it)))
            s=0.5_wp*(s+del*fsum)
         end if
      END SUBROUTINE wp_trapzd

      SUBROUTINE ep_trapzd(F,a,b,s,n,T,BCOEF,NB,K,ID,WORK,INBV)
      use bspline_base
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ID, K, NB, INBV
      REAL(kind=ep1), INTENT(INOUT) :: WORK(:)
      REAL(kind=ep1), INTENT(IN) :: BCOEF(:), T(:)
      class(bound_user_function) :: F
      REAL(kind=ep1), INTENT(IN) :: a,b
      REAL(kind=ep1), INTENT(INOUT) :: s
      INTEGER, INTENT(IN) :: n
      REAL(kind=ep1) :: del,fsum,x(max(1,2**(n-2)))
      INTEGER :: it, j
         if (n == 1) then
            s = 0.5_ep1*(b-a)*(F%eval(a)*BVALU(T,BCOEF,NB,K,ID,a,INBV,WORK) + F%eval(b)*BVALU(T,BCOEF,NB,K,ID,b,INBV,WORK))
            !s=0.5_ep1*(b-a)*sum(func( (/ a,b /) ))
         else
            it=2**(n-2)
            del=(b-a)/it
            x = cfp_arth(a+0.5_ep1*del,del,it)
            fsum = 0.0_ep1
            do j=1,it
               fsum = fsum + F%eval(x(j))*BVALU(T,BCOEF,NB,K,ID,x(j),INBV,WORK)
            enddo
            !fsum=sum(func(cfp_arth(a+0.5_wp*del,del,it)))
            s=0.5_ep1*(s+del*fsum)
         end if
      END SUBROUTINE ep_trapzd

      FUNCTION wp_arth(first,increment,n)
      IMPLICIT NONE
      REAL(kind=wp), INTENT(IN) :: first,increment
      INTEGER, INTENT(IN) :: n
      REAL(kind=wp), DIMENSION(n) :: wp_arth
      INTEGER :: k,k2
      REAL(kind=wp) :: temp
         if (n > 0) wp_arth(1)=first
         if (n <= NPAR_ARTH) then
            do k=2,n
               wp_arth(k)=wp_arth(k-1)+increment
            end do
         else
            do k=2,NPAR2_ARTH
               wp_arth(k)=wp_arth(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
               if (k >= n) exit
               k2=k+k
               wp_arth(k+1:min(k2,n))=temp+wp_arth(1:min(k,n-k))
               temp=temp+temp
               k=k2
            end do
         end if
      END FUNCTION wp_arth

      FUNCTION ep_arth(first,increment,n)
      IMPLICIT NONE
      REAL(kind=ep1), INTENT(IN) :: first,increment
      INTEGER, INTENT(IN) :: n
      REAL(kind=ep1), DIMENSION(n) :: ep_arth
      INTEGER :: k,k2
      REAL(kind=ep1) :: temp
         if (n > 0) ep_arth(1)=first
         if (n <= NPAR_ARTH) then
            do k=2,n
               ep_arth(k)=ep_arth(k-1)+increment
            end do
         else
            do k=2,NPAR2_ARTH
               ep_arth(k)=ep_arth(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
               if (k >= n) exit
               k2=k+k
               ep_arth(k+1:min(k2,n))=temp+ep_arth(1:min(k,n-k))
               temp=temp+temp
               k=k2
            end do
         end if
      END FUNCTION ep_arth

end module quadrature_module
