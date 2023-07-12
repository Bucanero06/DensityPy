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
module bspline_base
use utils
use precisn, only: cfp_dummy, f1mach

   private

   interface intrv
      module procedure dintrv, qintrv
   end interface

   interface bvalu
      module procedure dbvalu, qbvalu
   end interface

   public intrv, bvalu, map_knots_to_grid

contains 

  !> Find the indices mapping the start and end of each B-spline to the (quadrature) points r.
  !> For a B-spline with index ind which doesn't overlap with the quadrature grid the values in bspline_start_end_r are: bspline_start_end_r(1,ind) = 0, bspline_start_end_r(2,ind) = -1.
  subroutine map_knots_to_grid(knots,order,last_bspline,r,bspline_start_end_r)
     use precisn
     implicit none
     integer, intent(in) :: last_bspline, order
     real(kind=cfp), intent(in) :: knots(:), r(:)
     !OUTPUT:
     integer, allocatable :: bspline_start_end_r(:,:)

     integer :: i, j, err, ind, n_total_points

        n_total_points = size(r)

        if (allocated(bspline_start_end_r)) deallocate(bspline_start_end_r)

        allocate(bspline_start_end_r(2,last_bspline),stat=err)
        if (err .ne. 0) call xermsg('bspline_base','map_knots_to_grid','Memory allocation failed.',err,1)
        do ind=1,last_bspline

           bspline_start_end_r(1,ind) = 0

           do i=1,n_total_points
              if (r(i) > knots(ind)) then !the first r quadrature point beyond the B-spline starting point
                 bspline_start_end_r(1,ind) = i
                 exit
              endif
           enddo !i

           if (bspline_start_end_r(1,ind) .eq. 0) then !this B-spline doesn't overlap with the quadrature grid.
              bspline_start_end_r(2,ind) = -1
           else
              bspline_start_end_r(2,ind) = n_total_points
              do i=bspline_start_end_r(1,ind),n_total_points
                 if (r(i) > knots(ind+order)) then !the last r quadrature point beyond the B-spline end point was the previous one
                    bspline_start_end_r(2,ind) = i-1
                    exit
                 endif
              enddo !i
           endif

       enddo !ind

  end subroutine map_knots_to_grid

!> \verbatim
!>***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
!>            such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!>            the X interval.
!>***LIBRARY   SLATEC
!>***CATEGORY  E3, K6
!>***TYPE      real(kind=wp) (INTRV-S, DINTRV-D)
!>***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>     Written by Carl de Boor and modified by D. E. Amos
!>
!>     Abstract    **** a double precision routine ****
!>         DINTRV is the INTERV routine of the reference.
!>
!>         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
!>         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!>         the X interval.  Precisely,
!>
!>                      X .LT. XT(1)                1         -1
!>         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
!>           XT(LXT) .LE. X                         LXT        1,
!>
!>         That is, when multiplicities are present in the break point
!>         to the left of X, the largest index is taken for ILEFT.
!>
!>     Description of Arguments
!>
!>         Input      XT,X are double precision
!>          XT      - XT is a knot or break point vector of length LXT
!>          LXT     - length of the XT vector
!>          X       - argument
!>          ILO     - an initialization parameter which must be set
!>                    to 1 the first time the spline array XT is
!>                    processed by DINTRV.
!>
!>         Output
!>          ILO     - ILO contains information for efficient process-
!>                    ing after the initial call and ILO must not be
!>                    changed by the user.  Distinct splines require
!>                    distinct ILO parameters.
!>          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
!>          MFLAG   - signals when X lies out of bounds
!>
!>     Error Conditions
!>         None
!>
!>***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!>                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!>                 pp. 441-472.
!>***ROUTINES CALLED  (NONE)
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
      SUBROUTINE DINTRV (XT, LXT, X, ILO, ILEFT, MFLAG)
      use precisn
      IMPLICIT NONE
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      real(kind=wp) X, XT
      DIMENSION XT(*)
! ***FIRST EXECUTABLE STATEMENT  DINTRV
      ISTEP = 0
      MIDDLE = 0
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
! 
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
! 
!  *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
!  *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
! 
!  *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
!      NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
!  *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END SUBROUTINE

!>\verbatim
!>***PURPOSE  Evaluate the B-representation of a B-spline at X for the
!>            function value or any of its derivatives.
!>***LIBRARY   SLATEC
!>***CATEGORY  E3, K6
!>***TYPE      real(kind=cfp) (BVALU-S, DBVALU-D)
!>***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>     Written by Carl de Boor and modified by D. E. Amos
!>
!>     Abstract   **** a double precision routine ****
!>         DBVALU is the BVALUE function of the reference.
!>
!>         DBVALU evaluates the B-representation (T,A,N,K) of a B-spline
!>         at X for the function value on IDERIV=0 or any of its
!>         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!>         (right derivatives) are returned except at the right end
!>         point X=T(N+1) where left limiting values are computed.  The
!>         spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
!>         a fatal error message when X is outside of this interval.
!>
!>         To compute left derivatives or left limiting values at a
!>         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!>
!>         DBVALU calls DINTRV
!>
!>     Description of Arguments
!>
!>         Input      T,A,X are double precision
!>          T       - knot vector of length N+K
!>          A       - B-spline coefficient vector of length N
!>          N       - number of B-spline coefficients
!>                    N = sum of knot multiplicities-K
!>          K       - order of the B-spline, K .GE. 1
!>          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!>                    IDERIV = 0 returns the B-spline value
!>          X       - argument, T(K) .LE. X .LE. T(N+1)
!>          INBV    - an initialization parameter which must be set
!>                    to 1 the first time DBVALU is called.
!>
!>         Output     WORK,DBVALU are double precision
!>          INBV    - INBV contains information for efficient process-
!>                    ing after the initial call and INBV must not
!>                    be changed by the user.  Distinct splines require
!>                    distinct INBV parameters.
!>          WORK    - work vector of length 3*K.
!>          DBVALU  - value of the IDERIV-th derivative at X
!>
!>     Error Conditions
!>         An improper input is a fatal error
!>
!>***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!>                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!>                 pp. 441-472.
!>***ROUTINES CALLED  DINTRV, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
      real(kind=wp) FUNCTION DBVALU (T, A, N, K, IDERIV, X, INBV, WORK)
      use precisn
      IMPLICIT NONE
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ, IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      real(kind=wp) A, FKMJ, T, WORK, X, WTOL
      DIMENSION T(*), A(*), WORK(*)
! ***FIRST EXECUTABLE STATEMENT  DBVALU
      DBVALU = 0.0_wp
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
! 
!  *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!      (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
!ZM: relaxed the test to accommodate cases of negligible differences between X and T(I).
      WTOL = F1MACH(4,wp_dummy)
      IF (X-10*WTOL > T(I)) GO TO 130
!      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
! 
!  *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
!      WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
! 
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
! 
!  *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
!      GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)*WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 DBVALU = WORK(1)
      RETURN
! 
! 
  101 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'N DOES NOT SATISFY N.GE.K', 2, 1)
      RETURN
  102 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'K DOES NOT SATISFY K.GE.1', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
      RETURN
  130 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
      RETURN
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END FUNCTION

!> \verbatim
!>***PURPOSE  Calculate the value of the IDERIV-th derivative of the
!>            B-spline from the PP-representation.
!>***LIBRARY   SLATEC
!>***CATEGORY  E3, K6
!>***TYPE      real(kind=cfp) (PPVAL-S, DPPVAL-D)
!>***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>     Written by Carl de Boor and modified by D. E. Amos
!>
!>     Abstract    **** a double precision routine ****
!>         DPPVAL is the PPVALU function of the reference.
!>
!>         DPPVAL calculates (at X) the value of the IDERIV-th
!>         derivative of the B-spline from the PP-representation
!>         (C,XI,LXI,K).  The Taylor expansion about XI(J) for X in
!>         the interval XI(J) .LE. X .LT. XI(J+1) is evaluated, J=1,LXI.
!>         Right limiting values at X=XI(J) are obtained.  DPPVAL will
!>         extrapolate beyond XI(1) and XI(LXI+1).
!>
!>         To obtain left limiting values (left derivatives) at XI(J)
!>         replace LXI by J-1 and set X=XI(J),J=2,LXI+1.
!>
!>     Description of Arguments
!>
!>         Input      C,XI,X are double precision
!>          LDC     - leading dimension of C matrix, LDC .GE. K
!>          C       - matrix of dimension at least (K,LXI) containing
!>                    right derivatives at break points XI(*).
!>          XI      - break point vector of length LXI+1
!>          LXI     - number of polynomial pieces
!>          K       - order of B-spline, K .GE. 1
!>          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!>                    IDERIV=0 gives the B-spline value
!>          X       - argument, XI(1) .LE. X .LE. XI(LXI+1)
!>          INPPV   - an initialization parameter which must be set
!>                    to 1 the first time DPPVAL is called.
!>
!>         Output     DPPVAL is double precision
!>          INPPV   - INPPV contains information for efficient process-
!>                    ing after the initial call and INPPV must not
!>                    be changed by the user.  Distinct splines require
!>                    distinct INPPV parameters.
!>          DPPVAL  - value of the IDERIV-th derivative at X
!>
!>     Error Conditions
!>         Improper input is a fatal error
!>
!>***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!>                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!>                 pp. 441-472.
!>***ROUTINES CALLED  DINTRV, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
      real(kind=cfp) FUNCTION DPPVAL (LDC, C, XI, LXI, K, IDERIV, X, INPPV)
      use precisn
      IMPLICIT NONE
      INTEGER I, IDERIV, INPPV, J, K, LDC, LXI, NDUMMY, KK
      real(kind=wp) C, DX, X, XI
      DIMENSION XI(*), C(LDC,*)
! ***FIRST EXECUTABLE STATEMENT  DPPVAL
      DPPVAL = 0.0_wp
      IF(K.LT.1) GO TO 90
      IF(LDC.LT.K) GO TO 80
      IF(LXI.LT.1) GO TO 85
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 95
      I = K - IDERIV
      KK = I
      CALL DINTRV(XI, LXI, X, INPPV, I, NDUMMY)
      DX = X - XI(I)
      J = K
   10 DPPVAL = (DPPVAL/KK)*DX + C(J,I)
      J = J - 1
      KK = KK - 1
      IF (KK.GT.0) GO TO 10
      RETURN
! 
! 
   80 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'LDC DOES NOT SATISFY LDC.GE.K', 2, 1)
      RETURN
   85 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'LXI DOES NOT SATISFY LXI.GE.1', 2, 1)
      RETURN
   90 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'K DOES NOT SATISFY K.GE.1', 2, 1)
      RETURN
   95 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
      END FUNCTION
!> \verbatim
!>***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
!>            such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!>            the X interval.
!>***LIBRARY   SLATEC
!>***CATEGORY  E3, K6
!>***TYPE      real(kind=cfp) (INTRV-S, QINTRV-D)
!>***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>     Written by Carl de Boor and modified by D. E. Amos
!>
!>     Abstract    **** a double precision routine ****
!>         QINTRV is the INTERV routine of the reference.
!>
!>         QINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
!>         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!>         the X interval.  Precisely,
!>
!>                      X .LT. XT(1)                1         -1
!>         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
!>           XT(LXT) .LE. X                         LXT        1,
!>
!>         That is, when multiplicities are present in the break point
!>         to the left of X, the largest index is taken for ILEFT.
!>
!>     Description of Arguments
!>
!>         Input      XT,X are double precision
!>          XT      - XT is a knot or break point vector of length LXT
!>          LXT     - length of the XT vector
!>          X       - argument
!>          ILO     - an initialization parameter which must be set
!>                    to 1 the first time the spline array XT is
!>                    processed by QINTRV.
!>
!>         Output
!>          ILO     - ILO contains information for efficient process-
!>                    ing after the initial call and ILO must not be
!>                    changed by the user.  Distinct splines require
!>                    distinct ILO parameters.
!>          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
!>          MFLAG   - signals when X lies out of bounds
!>
!>     Error Conditions
!>         None
!>
!>***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!>                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!>                 pp. 441-472.
!>***ROUTINES CALLED  (NONE)
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
!
      SUBROUTINE QINTRV (XT, LXT, X, ILO, ILEFT, MFLAG)
!
!
      use precisn
      IMPLICIT NONE
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      REAL(kind=ep1) :: X, XT
      DIMENSION XT(*)
!***FIRST EXECUTABLE STATEMENT  QINTRV
      ISTEP = 0
      MIDDLE = 0
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
!
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
!
! *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
! *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
!
! *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
!     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
! *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END SUBROUTINE
!> \verbatim
!>***PURPOSE  Evaluate the B-representation of a B-spline at X for the
!>            function value or any of its derivatives.
!>***LIBRARY   SLATEC
!>***CATEGORY  E3, K6
!>***TYPE      REAL(kind=ep1) (BVALU-S, QBVALU-D)
!>***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>     Written by Carl de Boor and modified by D. E. Amos
!>
!>     Abstract   **** a double precision routine ****
!>         QBVALU is the BVALUE function of the reference.
!>
!>         QBVALU evaluates the B-representation (T,A,N,K) of a B-spline
!>         at X for the function value on IDERIV=0 or any of its
!>         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!>         (right derivatives) are returned except at the right end
!>         point X=T(N+1) where left limiting values are computed.  The
!>         spline is defined on T(K) .LE. X .LE. T(N+1).  QBVALU returns
!>         a fatal error message when X is outside of this interval.
!>
!>         To compute left derivatives or left limiting values at a
!>         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!>
!>         QBVALU calls QINTRV
!>
!>     Description of Arguments
!>
!>         Input      T,A,X are double precision
!>          T       - knot vector of length N+K
!>          A       - B-spline coefficient vector of length N
!>          N       - number of B-spline coefficients
!>                    N = sum of knot multiplicities-K
!>          K       - order of the B-spline, K .GE. 1
!>          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!>                    IDERIV = 0 returns the B-spline value
!>          X       - argument, T(K) .LE. X .LE. T(N+1)
!>          INBV    - an initialization parameter which must be set
!>                    to 1 the first time QBVALU is called.
!>
!>         Output     WORK,QBVALU are double precision
!>          INBV    - INBV contains information for efficient process-
!>                    ing after the initial call and INBV must not
!>                    be changed by the user.  Distinct splines require
!>                    distinct INBV parameters.
!>          WORK    - work vector of length 3*K.
!>          QBVALU  - value of the IDERIV-th derivative at X
!>
!>     Error Conditions
!>         An improper input is a fatal error
!>
!>***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!>                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!>                 pp. 441-472.
!>***ROUTINES CALLED  QINTRV, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800901  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
!
      REAL(kind=ep1) FUNCTION QBVALU (T, A, N, K, IDERIV, X, INBV, WORK)
!
!
      use precisn
      IMPLICIT NONE
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ, IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      REAL(kind=ep1) :: A, FKMJ, T, WORK, X
      DIMENSION T(*), A(*), WORK(*)
!***FIRST EXECUTABLE STATEMENT  QBVALU
      QBVALU = 0.0_ep1
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
!
! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL QINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
!
! *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
!     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
!
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
!
! *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
!     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)*WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 QBVALU = WORK(1)
      RETURN
!
!
  101 CONTINUE
      CALL XERMSG ('SLATEC', 'QBVALU', 'N DOES NOT SATISFY N.GE.K', 2, 1)
      RETURN
  102 CONTINUE
      CALL XERMSG ('SLATEC', 'QBVALU', 'K DOES NOT SATISFY K.GE.1', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'QBVALU', 'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'QBVALU', 'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
      RETURN
  130 CONTINUE
      CALL XERMSG ('SLATEC', 'QBVALU', 'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
      RETURN
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'QBVALU', 'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END FUNCTION

end module bspline_base
