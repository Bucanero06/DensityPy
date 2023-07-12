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
module general_quadrature
use precisn
use utils, only: xermsg
use const, only: limit

 private

 public bound_user_function, function_2d, n_7, x_7, w_7, n_10, x_10, w_10, gl_expand_A_B
 
 !DQAGS is the routine that performs the numerical quadrature
 public DQAGS

 public quad2d, gl2d, dqelg

 !> \class <bound_user_function>
 !> This is a class that defines an abstract function of one variable, whose specific implementation is deferred.
 !> The purpose is to use this abstract function in some bspline-related routines which require a user-defined function as a parameters.
 type, abstract :: bound_user_function
   !no data components
 contains
      !> This must be used in all routines using this object to evaluate the function. This symbol is resolved into one of wp_eval,
      !> ep_eval depending on the floating point type of X on input.
      generic, public :: eval => wp_eval, ep_eval
      !> \memberof bound_user_function
      procedure(wp_user_function_interface), deferred :: wp_eval
      !> \memberof bound_user_function
      procedure(ep_user_function_interface), deferred :: ep_eval
 end type bound_user_function
 abstract interface
      real(wp) function wp_user_function_interface(data,x)
         import :: bound_user_function, wp
         class(bound_user_function) :: data
         real(kind=wp), intent(in) :: x
      end function wp_user_function_interface
 end interface
 abstract interface
      real(ep1) function ep_user_function_interface(data,x)
         import :: bound_user_function, ep1
         class(bound_user_function) :: data
         real(kind=ep1), intent(in) :: x
      end function ep_user_function_interface
 end interface

 !> \class <function_2d>
 !> This is a class that defines an abstract function of two variables, whose specific implementation is deferred.
 !> The purpose is to use this abstract function in some bspline-related routines which require a user-defined function as a parameters.
 type, abstract :: function_2d
    !> Number of function evaulation, number of subdivisions (applicable for adaptive quadratures).
    integer :: neval = 0, ndiv = 0
    !> Maximum number of sub-divisions of the area to integrate over: currently this is not used in quad2d.
    integer :: max_div = 2*limit
 contains
    procedure(fn2d_eval_interface), deferred :: eval
 end type function_2d
 abstract interface
      real(cfp) function fn2d_eval_interface(this,x,y)
         import :: function_2d, cfp
         class(function_2d) :: this
         real(kind=cfp), intent(in) :: x, y
      end function fn2d_eval_interface
 end interface

 !> Order of the Gauss-Legendre quadrature to which the x_7 and w_7 arrays correspond.
 integer, parameter :: n_7 = 7

 !> Weights for the Gauss-Legendre quadrature of order 7 on interval [0,1].
 real(kind=cfp), parameter :: w_7(2*n_7+1) = (/0.015376620998058634177314196788602209_cfp,&
                                           &0.035183023744054062354633708225333669_cfp,&
                                           &0.05357961023358596750593477334293465_cfp,&
                                           &0.06978533896307715722390239725551416_cfp,&
                                           &0.08313460290849696677660043024060441_cfp,&
                                           &0.09308050000778110551340028093321141_cfp,&
                                           &0.09921574266355578822805916322191966_cfp,&
                                           &0.10128912096278063644031009998375966_cfp,&
                                           &0.09921574266355578822805916322191966_cfp,&
                                           &0.09308050000778110551340028093321141_cfp,&
                                           &0.08313460290849696677660043024060441_cfp,&
                                           &0.06978533896307715722390239725551416_cfp,&
                                           &0.05357961023358596750593477334293465_cfp,&
                                           &0.035183023744054062354633708225333669_cfp,&
                                           &0.015376620998058634177314196788602209_cfp/)

 !> Abscissas for the Gauss-Legendre quadrature of order 7 on interval [0,1].
 real(kind=cfp), parameter :: x_7(2*n_7+1) = (/0.0060037409897572857552171407066937094_cfp,&
                                          &0.031363303799647047846120526144895264_cfp,&
                                          &0.075896708294786391899675839612891574_cfp,&
                                          &0.13779113431991497629190697269303100_cfp,&
                                          &0.21451391369573057623138663137304468_cfp,&
                                          &0.30292432646121831505139631450947727_cfp,&
                                          &0.39940295300128273884968584830270190_cfp,&
                                          &0.50000000000000000000000000000000000_cfp,&
                                          &0.60059704699871726115031415169729810_cfp,&
                                          &0.69707567353878168494860368549052273_cfp,&
                                          &0.78548608630426942376861336862695532_cfp,&
                                          &0.86220886568008502370809302730696900_cfp,&
                                          &0.92410329170521360810032416038710843_cfp,&
                                          &0.96863669620035295215387947385510474_cfp,&
                                          &0.99399625901024271424478285929330629_cfp/)

 !> Order of the Gauss-Legendre quadrature to which the x_10 and w_10 arrays correspond.
 integer, parameter :: n_10 = 10

 !> Abscissas for the Gauss-Legendre quadrature of order 10 on interval [0,1].
 real(kind=cfp), parameter :: x_10(2*n_10+1) = (/0.0031239146898052498698789820310295354_cfp,&
                                        &0.016386580716846852841688892546152419_cfp,&
                                        &0.039950332924799585604906433142515553_cfp,&
                                        &0.073318317708341358176374680706216165_cfp,&
                                        &0.11578001826216104569206107434688598_cfp,&
                                        &0.16643059790129384034701666500483042_cfp,&
                                        &0.22419058205639009647049060163784336_cfp,&
                                        &0.28782893989628060821316555572810597_cfp,&
                                        &0.35598934159879945169960374196769984_cfp,&
                                        &0.42721907291955245453148450883065683_cfp,&
                                        &0.50000000000000000000000000000000000_cfp,&
                                        &0.57278092708044754546851549116934317_cfp,&
                                        &0.64401065840120054830039625803230016_cfp,&
                                        &0.71217106010371939178683444427189403_cfp,&
                                        &0.77580941794360990352950939836215664_cfp,&
                                        &0.83356940209870615965298333499516958_cfp,&
                                        &0.88421998173783895430793892565311402_cfp,&
                                        &0.92668168229165864182362531929378384_cfp,&
                                        &0.96004966707520041439509356685748445_cfp,&
                                        &0.98361341928315314715831110745384758_cfp,&
                                        &0.99687608531019475013012101796897046_cfp/)

 !> Weights for the Gauss-Legendre quadrature of order 10 on interval [0,1].
 real(kind=cfp), parameter :: w_10(2*n_10+1) =  (/0.008008614128887166662112308429235508_cfp,&
                                           &0.018476894885426246899975334149664833_cfp,&
                                           &0.028567212713428604141817913236223979_cfp,&
                                           &0.038050056814189651008525826650091590_cfp,&
                                           &0.046722211728016930776644870556966044_cfp,&
                                           &0.05439864958357418883173728903505282_cfp,&
                                           &0.06091570802686426709768358856286680_cfp,&
                                           &0.06613446931666873089052628724838780_cfp,&
                                           &0.06994369739553657736106671193379156_cfp,&
                                           &0.07226220199498502953191358327687627_cfp,&
                                           &0.07304056682484521359599257384168559_cfp,&
                                           &0.07226220199498502953191358327687627_cfp,&
                                           &0.06994369739553657736106671193379156_cfp,&
                                           &0.06613446931666873089052628724838780_cfp,&
                                           &0.06091570802686426709768358856286680_cfp,&
                                           &0.05439864958357418883173728903505282_cfp,&
                                           &0.046722211728016930776644870556966044_cfp,&
                                           &0.038050056814189651008525826650091590_cfp,&
                                           &0.028567212713428604141817913236223979_cfp,&
                                           &0.018476894885426246899975334149664833_cfp,&
                                           &0.008008614128887166662112308429235508_cfp/)

contains

  !> Takes the Gauss-Legendre rule for the interval [0,1] and expands it for the given interval [A,B].
  subroutine gl_expand_A_B(x,w,n,x_AB,w_AB,A,B)
     implicit none
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: A, B
     real(kind=cfp), intent(in) :: x(2*n+1), w(2*n+1)
     real(kind=cfp), intent(out) :: x_AB(2*n+1), w_AB(2*n+1)

     integer :: i
     real(kind=cfp) :: delta

        delta = B-A
        do i=1,2*n+1
           x_AB(i) = x(i)*delta + A
           w_AB(i) = w(i)*delta
        enddo !i

  end subroutine gl_expand_A_B

 !> Adaptive 2D quadrature on rectangle based on Gauss-Kronrod rule. The algorithm is based on that of Romanowski published in Int. J. Q. Chem.
 !> \param [in] f The 2D function to be integrated.
 !> \param [in] Ax Rectangle X-coordinate start.
 !> \param [in] Bx Rectangle X-coordinate end.
 !> \param [in] Ay Rectangle Y-coordinate start.
 !> \param [in] By Rectangle Y-coordinate end.
 !> \param [in] eps Required relative precision for the integral.
 !> \param [in] Qest Optional: an estimate of the integral over the specified rectangle as obtained by a call to gl2d. 
 !>                  Note that the estimate must be obtained using gl2d for the algorithm to proceed correctly since it relies
 !>                  on the fixed-point G-K quadrature to calculate integrals on the sub-rectangles.
 recursive function quad2d(f,Ax,Bx,Ay,By,eps,Qest) result(I)
     implicit none
     class(function_2d) :: f
     real(kind=cfp), intent(in) :: Ax,Bx,Ay,By,eps
     real(kind=cfp), optional, intent(in) :: Qest

     real(kind=cfp) :: Xhalf, Yhalf, Lx, Ly, QP1, QP2, QP3, QP4, Q, I

        !Quit if too many sub-divisions have been performed: todo this should trigger an error/warning message!!
!        if (f%ndiv > f%max_div) return

        Xhalf = (Ax+Bx)*0.5_cfp
        Yhalf = (Ay+By)*0.5_cfp

        Lx = Bx-Ax
        Ly = By-Ay

        !If the estimate of the quadrature on this rectangle using a previous call to gl2d is given then use it otherwise calculate it.
        if (present(Qest)) then
           Q = Qest
        else
           Q = gl2d(f,Ax,Bx,Ay,By)
        endif

        !Calculate quadratures on the four sub-rectangles
        QP1 = gl2d(f,Ax,Xhalf,Ay,Yhalf)
        QP2 = gl2d(f,Ax,Xhalf,Yhalf,By)
        QP3 = gl2d(f,Xhalf,Bx,Yhalf,By)
        QP4 = gl2d(f,Xhalf,Bx,Ay,Yhalf)

        I = QP1+QP2+QP3+QP4

        !Continue recursively on each rectangle if the desired precision has not been reached
        !todo instead of dividing all four try dividing the largest one, then the second, etc. and see if we converge without improving all at the same time.
        if (abs((Q-I)/I) > eps) then
           Q = quad2d(f,Ax,Xhalf,Ay,Yhalf,eps,QP1)
           I = Q
           Q = quad2d(f,Ax,Xhalf,Yhalf,By,eps,QP2)
           I = I + Q
           Q = quad2d(f,Xhalf,Bx,Yhalf,By,eps,QP3)
           I = I + Q
           Q = quad2d(f,Xhalf,Bx,Ay,Yhalf,eps,QP4)
           I = I + Q
        endif

        f%ndiv = f%ndiv + 1

 end function quad2d

 !> 2D Quadrature on rectangle using the Gauss-Kronrod rule of order 8. The meaning of the input variables is identical to quad2d input parameters.
 function gl2d(f,Ax,Bx,Ay,By)
     implicit none
     class(function_2d) :: f
     real(kind=cfp), intent(in) :: Ax,Bx,Ay,By
     real(kind=cfp) :: gl2d

     !Abscissae and weights for the Gauss-Kronrod rule for interval [0,1] obtained using the Mathematica command:
     !n2d=8;
     !crule=NIntegrate`CartesianRuleData[{{"GaussKronrodRule","GaussPoints"->n2d},{"GaussKronrodRule","GaussPoints"->n2d}},36]
     integer, parameter :: n = 8
     real(kind=cfp), parameter :: x(2*n+1) = (/ &
        0.003310062059141922032055965490164602_cfp, &
        0.019855071751231884158219565715263505_cfp, &
        0.052939546576271789025819491230874338_cfp, &
        0.101666761293186630204223031762084782_cfp, &
        0.163822964527420661421844630953584475_cfp, &
        0.237233795041835507091130475405376825_cfp, &
        0.319649451035934021403725688851554280_cfp, &
        0.408282678752175097530261928819908010_cfp, &
        0.500000000000000000000000000000000000_cfp, &
        0.591717321247824902469738071180091990_cfp, &
        0.680350548964065978596274311148445720_cfp, &
        0.762766204958164492908869524594623175_cfp, &
        0.836177035472579338578155369046415525_cfp, &
        0.898333238706813369795776968237915218_cfp, &
        0.947060453423728210974180508769125662_cfp, &
        0.980144928248768115841780434284736495_cfp, &
        0.996689937940858077967944034509835398_cfp  &
     /)
     real(kind=cfp), parameter :: w(2*n+1) = (/ &
        0.00891119166035517757639348060137489490_cfp, &
        0.024719697501069654250181984723498447_cfp,   &
        0.0412411494656791653443125967228039477_cfp,  &
        0.055823185413419806611054079466970742_cfp,   &
        0.0681315546275861076311693726272531016_cfp,  &
        0.078326303084094200245124044243484369_cfp,   &
        0.0860353042776056559286474401019285433_cfp,  &
        0.09070001253401732153087426258627522_cfp,    &
        0.0922232028723458217644854778528214649_cfp,  &
        0.09070001253401732153087426258627522_cfp,    &
        0.0860353042776056559286474401019285433_cfp,  &
        0.078326303084094200245124044243484369_cfp,   &
        0.0681315546275861076311693726272531016_cfp,  &
        0.055823185413419806611054079466970742_cfp,   &
        0.0412411494656791653443125967228039477_cfp,  &
        0.024719697501069654250181984723498447_cfp,   &
        0.00891119166035517757639348060137489490_cfp  &
     /)

     integer :: i, j
     real(kind=cfp) :: Lx, Ly

        Lx = Bx-Ax
        Ly = By-Ay

        gl2d = 0.0_cfp
        do i=1,2*n+1
           do j=1,2*n+1
              gl2d = gl2d + w(i)*w(j)*f%eval(x(i)*Lx+Ax,x(j)*Ly+Ay)
           enddo !j
        enddo !i
        gl2d = gl2d*Lx*Ly

 end function gl2d

!>***BEGIN PROLOGUE  DQAGS
!>***PURPOSE  The routine calculates an approximation result to a given
!>            Definite integral  I = Integral of F over (A,B),
!>            Hopefully satisfying following claim for accuracy
!>            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A1
!>***TYPE      real(kind=cfp) (QAGS-S, DQAGS-D)
!>***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!>             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
!>             QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>        Computation of a definite integral
!>        Standard fortran subroutine
!>        Double precision version
!>
!>
!>        PARAMETERS
!>         ON ENTRY
!>            F      - class(bound_user_function)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>            A      - Double precision
!>                     Lower limit of integration
!>
!>            B      - Double precision
!>                     Upper limit of integration
!>
!>            EPSABS - Double precision
!>                     Absolute accuracy requested
!>            EPSREL - Double precision
!>                     Relative accuracy requested
!>                     If  EPSABS.LE.0
!>                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
!>                     The routine will end with IER = 6.
!>
!>         ON RETURN
!>            RESULT - Double precision
!>                     Approximation to the integral
!>
!>            ABSERR - Double precision
!>                     Estimate of the modulus of the absolute error,
!>                     which should equal or exceed ABS(I-RESULT)
!>
!>            NEVAL  - Integer
!>                     Number of integrand evaluations
!>
!>            IER    - Integer
!>                     IER = 0 Normal and reliable termination of the
!>                             routine. It is assumed that the requested
!>                             accuracy has been achieved.
!>                     IER.GT.0 Abnormal termination of the routine
!>                             The estimates for integral and error are
!>                             less reliable. It is assumed that the
!>                             requested accuracy has not been achieved.
!>            ERROR MESSAGES
!>                     IER = 1 Maximum number of subdivisions allowed
!>                             has been achieved. One can allow more sub-
!>                             divisions by increasing the value of LIMIT
!>                             (and taking the according dimension
!>                             adjustments into account. However, if
!>                             this yields no improvement it is advised
!>                             to analyze the integrand in order to
!>                             determine the integration difficulties. If
!>                             the position of a local difficulty can be
!>                             determined (E.G. SINGULARITY,
!>                             DISCONTINUITY WITHIN THE INTERVAL) one
!>                             will probably gain from splitting up the
!>                             interval at this point and calling the
!>                             integrator on the subranges. If possible,
!>                             an appropriate special-purpose integrator
!>                             should be used, which is designed for
!>                             handling the type of difficulty involved.
!>                         = 2 The occurrence of roundoff error is detec-
!>                             ted, which prevents the requested
!>                             tolerance from being achieved.
!>                             The error may be under-estimated.
!>                         = 3 Extremely bad integrand behaviour
!>                             occurs at some points of the integration
!>                             interval.
!>                         = 4 The algorithm does not converge.
!>                             Roundoff error is detected in the
!>                             Extrapolation table. It is presumed that
!>                             the requested tolerance cannot be
!>                             achieved, and that the returned result is
!>                             the best which can be obtained.
!>                         = 5 The integral is probably divergent, or
!>                             slowly convergent. It must be noted that
!>                             divergence can occur with any other value
!>                             of IER.
!>                         = 6 The input is invalid, because
!>                             (EPSABS.LE.0 AND
!>                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)
!>                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
!>                             RESULT, ABSERR, NEVAL, LAST are set to
!>                             zero.  Except when LIMIT or LENW is
!>                             invalid, IWORK(1), WORK(LIMIT*2+1) and
!>                             WORK(LIMIT*3+1) are set to zero, WORK(1)
!>                             is set to A and WORK(LIMIT+1) TO B.
!>
!>         DIMENSIONING PARAMETERS
!>            LIMIT - Integer
!>                    DIMENSIONING PARAMETER FOR IWORK
!>                    LIMIT determines the maximum number of subintervals
!>                    in the partition of the given integration interval
!>                    (A,B), LIMIT.GE.1.
!>                    IF LIMIT.LT.1, the routine will end with IER = 6.
!>
!>            LENW  - Integer
!>                    DIMENSIONING PARAMETER FOR WORK
!>                    LENW must be at least LIMIT*4.
!>                    If LENW.LT.LIMIT*4, the routine will end
!>                    with IER = 6.
!>
!>            LAST  - Integer
!>                    On return, LAST equals the number of subintervals
!>                    produced in the subdivision process, determines the
!>                    number of significant elements actually in the WORK
!>                    Arrays.
!>
!>         WORK ARRAYS
!>            IWORK - Integer
!>                    Vector of dimension at least LIMIT, the first K
!>                    elements of which contain pointers
!>                    to the error estimates over the subintervals
!>                    such that WORK(LIMIT*3+IWORK(1)),... ,
!>                    WORK(LIMIT*3+IWORK(K)) form a decreasing
!>                    sequence, with K = LAST IF LAST.LE.(LIMIT/2+2),
!>                    and K = LIMIT+1-LAST otherwise
!>
!>            WORK  - Double precision
!>                    Vector of dimension at least LENW
!>                    on return
!>                    WORK(1), ..., WORK(LAST) contain the left
!>                     end-points of the subintervals in the
!>                     partition of (A,B),
!>                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
!>                     the right end-points,
!>                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
!>                     the integral approximations over the subintervals,
!>                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
!>                     contain the error estimates.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  DQAGSE, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>***END PROLOGUE  DQAGS
!>
!>
      SUBROUTINE DQAGS (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      class(bound_user_function) :: F
      real(kind=cfp) A,ABSERR,B,EPSABS,EPSREL,RESULT,WORK!,F
      INTEGER IER,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
!
      DIMENSION IWORK(*),WORK(*)
!
      !DECLARE F AS: real(kind=cfp), EXTERNAL :: F
      !EXTERNAL F
!
!         CHECK VALIDITY OF LIMIT AND LENW.
!
!***FIRST EXECUTABLE STATEMENT  DQAGS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0_cfp
      ABSERR = 0.0_cfp
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
!
!         PREPARE CALL FOR DQAGSE.
!
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
!
      CALL DQAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
!
!         CALL ERROR HANDLER IF NECESSARY.
!
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) THEN
         PRINT *,RESULT
         CALL XERMSG ('SLATEC', 'DQAGS', 'ABNORMAL RETURN', IER, LVL)
      ENDIF
      RETURN
      END SUBROUTINE DQAGS

!>***BEGIN PROLOGUE  DQAGSE
!>***PURPOSE  The routine calculates an approximation result to a given
!>            definite integral I = Integral of F over (A,B),
!>            hopefully satisfying following claim for accuracy
!>            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A1
!>***TYPE      real(kind=cfp) (QAGSE-S, DQAGSE-D)
!>***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!>             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
!>             QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>        Computation of a definite integral
!>        Standard fortran subroutine
!>        Double precision version
!>
!>        PARAMETERS
!>         ON ENTRY
!>            F      - class(bound_user_function)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>            A      - Double precision
!>                     Lower limit of integration
!>
!>            B      - Double precision
!>                     Upper limit of integration
!>
!>            EPSABS - Double precision
!>                     Absolute accuracy requested
!>            EPSREL - Double precision
!>                     Relative accuracy requested
!>                     If  EPSABS.LE.0
!>                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
!>                     the routine will end with IER = 6.
!>
!>            LIMIT  - Integer
!>                     Gives an upper bound on the number of subintervals
!>                     in the partition of (A,B)
!>
!>         ON RETURN
!>            RESULT - Double precision
!>                     Approximation to the integral
!>
!>            ABSERR - Double precision
!>                     Estimate of the modulus of the absolute error,
!>                     which should equal or exceed ABS(I-RESULT)
!>
!>            NEVAL  - Integer
!>                     Number of integrand evaluations
!>
!>            IER    - Integer
!>                     IER = 0 Normal and reliable termination of the
!>                             routine. It is assumed that the requested
!>                             accuracy has been achieved.
!>                     IER.GT.0 Abnormal termination of the routine
!>                             the estimates for integral and error are
!>                             less reliable. It is assumed that the
!>                             requested accuracy has not been achieved.
!>            ERROR MESSAGES
!>                         = 1 Maximum number of subdivisions allowed
!>                             has been achieved. One can allow more sub-
!>                             divisions by increasing the value of LIMIT
!>                             (and taking the according dimension
!>                             adjustments into account). However, if
!>                             this yields no improvement it is advised
!>                             to analyze the integrand in order to
!>                             determine the integration difficulties. If
!>                             the position of a local difficulty can be
!>                             determined (e.g. singularity,
!>                             discontinuity within the interval) one
!>                             will probably gain from splitting up the
!>                             interval at this point and calling the
!>                             integrator on the subranges. If possible,
!>                             an appropriate special-purpose integrator
!>                             should be used, which is designed for
!>                             handling the type of difficulty involved.
!>                         = 2 The occurrence of roundoff error is detec-
!>                             ted, which prevents the requested
!>                             tolerance from being achieved.
!>                             The error may be under-estimated.
!>                         = 3 Extremely bad integrand behaviour
!>                             occurs at some points of the integration
!>                             interval.
!>                         = 4 The algorithm does not converge.
!>                             Roundoff error is detected in the
!>                             extrapolation table.
!>                             It is presumed that the requested
!>                             tolerance cannot be achieved, and that the
!>                             returned result is the best which can be
!>                             obtained.
!>                         = 5 The integral is probably divergent, or
!>                             slowly convergent. It must be noted that
!>                             divergence can occur with any other value
!>                             of IER.
!>                         = 6 The input is invalid, because
!>                             EPSABS.LE.0 and
!>                             EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28).
!>                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
!>                             IORD(1) and ELIST(1) are set to zero.
!>                             ALIST(1) and BLIST(1) are set to A and B
!>                             respectively.
!>
!>            ALIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the left end points
!>                     of the subintervals in the partition of the
!>                     given integration range (A,B)
!>
!>            BLIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the right end points
!>                     of the subintervals in the partition of the given
!>                     integration range (A,B)
!>
!>            RLIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the integral
!>                     approximations on the subintervals
!>
!>            ELIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the moduli of the
!>                     absolute error estimates on the subintervals
!>
!>            IORD   - Integer
!>                     Vector of dimension at least LIMIT, the first K
!>                     elements of which are pointers to the
!>                     error estimates over the subintervals,
!>                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
!>                     form a decreasing sequence, with K = LAST
!>                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
!>                     otherwise
!>
!>            LAST   - Integer
!>                     Number of subintervals actually produced in the
!>                     subdivision process
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, DQELG, DQK21, DQPSRT
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  DQAGSE
!>
!>            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
!>            LIMEXP IN SUBROUTINE DQELG (RLIST2 SHOULD BE OF DIMENSION
!>            (LIMEXP+2) AT LEAST).
!>
!>            LIST OF MAJOR VARIABLES
!>            -----------------------
!>
!>           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
!>                       CONSIDERED UP TO NOW
!>           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
!>                       CONSIDERED UP TO NOW
!>           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
!>                       (ALIST(I),BLIST(I))
!>           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2 CONTAINING
!>                       THE PART OF THE EPSILON TABLE WHICH IS STILL
!>                       NEEDED FOR FURTHER COMPUTATIONS
!>           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!>           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
!>                       ESTIMATE
!>           ERRMAX    - ELIST(MAXERR)
!>           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
!>                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
!>           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
!>           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
!>           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
!>                       ABS(RESULT))
!>           *****1    - VARIABLE FOR THE LEFT INTERVAL
!>           *****2    - VARIABLE FOR THE RIGHT INTERVAL
!>           LAST      - INDEX FOR SUBDIVISION
!>           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
!>           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
!>                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
!>                       INTEGRAL HAS BEEN OBTAINED IT IS PUT IN
!>                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
!>                       BY ONE.
!>           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
!>                       TO NOW, MULTIPLIED BY 1.5
!>           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
!>                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
!>           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
!>                       ATTEMPTING TO PERFORM EXTRAPOLATION I.E. BEFORE
!>                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
!>                       DECREASE THE VALUE OF ERLARG.
!>           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
!>                       IS NO LONGER ALLOWED (TRUE VALUE)
!>
!>            MACHINE DEPENDENT CONSTANTS
!>            ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!>           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
      SUBROUTINE DQAGSE (F, A, B, EPSABS, EPSREL, LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
      use const, only: max_epstab
      class(bound_user_function) :: F
      real(kind=cfp) A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,BLIST,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2, &
     &  DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,ERRMAX, &
     &  ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,OFLOW,RESABS,RESEPS,RESULT,&
     &  RES3LA,RLIST,RLIST2,SMALL,UFLOW
      real(kind=cfp) cfp_dummy
      INTEGER ID,IER,IERRO,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,KSGN,KTMIN,LAST,LIMIT,MAXERR,NEVAL,NRES,NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
!
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),RES3LA(3),RLIST(*),RLIST2(max_epstab)
!
!***FIRST EXECUTABLE STATEMENT  DQAGSE
      EPMACH = F1MACH(4,cfp_dummy)
!
!            TEST ON VALIDITY OF PARAMETERS
!            ------------------------------
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0_cfp
      ABSERR = 0.0_cfp
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0_cfp
      ELIST(1) = 0.0_cfp
      IF(EPSABS.LE.0.0_cfp.AND.EPSREL.LT.MAX(0.5E+02_cfp*EPMACH,0.5E-28_cfp)) IER = 6
      IF(IER.EQ.6) GO TO 999
!
!           FIRST APPROXIMATION TO THE INTEGRAL
!           -----------------------------------
!
      UFLOW = F1MACH(1,cfp_dummy)
      OFLOW = F1MACH(2,cfp_dummy)
      IERRO = 0
      CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
!
!           TEST ON ACCURACY.
!
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      IF(ABSERR.LE.1.0E+02_cfp*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.ABSERR.EQ.0.0_cfp) GO TO 140
!
!           INITIALIZATION
!           --------------
!
      RLIST2(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      ABSERR = OFLOW
      NRMAX = 1
      NRES = 0
      NUMRL2 = 2
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1E+01_cfp-0.5E+02_cfp*EPMACH)*DEFABS) KSGN = 1
!
!           MAIN DO-LOOP
!           ------------
!
      DO 90 LAST = 2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
!           ESTIMATE.
!
        A1 = ALIST(MAXERR)
        B1 = 0.5_cfp*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 15
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1E-04_cfp*ABS(AREA12).OR.ERRO12.LT.0.99_cfp*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
!
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
!
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
!           EQUALS LIMIT.
!
        IF(LAST.EQ.LIMIT) IER = 1
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT A POINT OF THE INTEGRATION RANGE.
!
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1E+01_cfp+0.1E+03_cfp*EPMACH)*(ABS(A2)+0.1E+04_cfp*UFLOW)) IER = 4
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
!
!           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
!           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
!
   30   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 115
! ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(ABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
!
!           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
!           SMALLEST INTERVAL.
!
        IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
!
!           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
!           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
!
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
! ***JUMP OUT OF DO-LOOP
          IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
!
!           PERFORM EXTRAPOLATION.
!
   60   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1E-02_cfp*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
! ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LE.ERTEST) GO TO 100
!
!           PREPARE BISECTION OF THE SMALLEST INTERVAL.
!
   70   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5_cfp
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = ABS(B-A)*0.375_cfp
        ERLARG = ERRSUM
        ERTEST = ERRBND
        RLIST2(2) = AREA
   90 CONTINUE
!
!           SET FINAL RESULT AND ERROR ESTIMATE.
!           ------------------------------------
!
  100 IF(ABSERR.EQ.OFLOW) GO TO 115
      IF(IER+IERRO.EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0_cfp.AND.AREA.NE.0.0_cfp) GO TO 105
      IF(ABSERR.GT.ERRSUM) GO TO 115
      IF(AREA.EQ.0.0_cfp) GO TO 130
      GO TO 110
  105 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA)) GO TO 115
!
!           TEST ON DIVERGENCE.
!
  110 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.DEFABS*0.1E-01_cfp) GO TO 130
      IF(0.1E-01_cfp.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1E+03_cfp.OR.ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 130
!
!           COMPUTE GLOBAL INTEGRAL SUM.
!
  115 RESULT = 0.0_cfp
      DO 120 K = 1,LAST
         RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 IF(IER.GT.2) IER = IER-1
  140 NEVAL = 42*LAST-21
  999 RETURN
      END SUBROUTINE DQAGSE

!>***BEGIN PROLOGUE  DQELG
!>***SUBSIDIARY
!>***PURPOSE  The routine determines the limit of a given sequence of
!>            approximations, by means of the Epsilon algorithm of
!>            P.Wynn. An estimate of the absolute error is also given.
!>            The condensed Epsilon table is computed. Only those
!>            elements needed for the computation of the next diagonal
!>            are preserved.
!>***LIBRARY   SLATEC
!>***TYPE      real(kind=cfp) (QELG-S, DQELG-D)
!>***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Epsilon algorithm
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS
!>              N      - Integer
!>                       EPSTAB(N) contains the new element in the
!>                       first column of the epsilon table.
!>
!>              EPSTAB - Double precision
!>                       Vector of dimension 52 containing the elements
!>                       of the two lower diagonals of the triangular
!>                       epsilon table. The elements are numbered
!>                       starting at the right-hand corner of the
!>                       triangle.
!>
!>              RESULT - Double precision
!>                       Resulting approximation to the integral
!>
!>              ABSERR - Double precision
!>                       Estimate of the absolute error computed from
!>                       RESULT and the 3 previous results
!>
!>              RES3LA - Double precision
!>                       Vector of dimension 3 containing the last 3
!>                       results
!>
!>              NRES   - Integer
!>                       Number of calls to the routine
!>                       (should be zero at first call)
!>
!>***SEE ALSO  DQAGIE, DQAGOE, DQAGPE, DQAGSE
!>***ROUTINES CALLED  F1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>***END PROLOGUE  DQELG
!>
!>           LIST OF MAJOR VARIABLES
!>           -----------------------
!>
!>           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
!>           E1       ELEMENT IN THE EPSILON TABLE IS BASED
!>           E2
!>           E3                 E0
!>                        E3    E1    NEW
!>                              E2
!>           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
!>                    DIAGONAL
!>           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
!>           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
!>                    OF ERROR
!>
!>           MACHINE DEPENDENT CONSTANTS
!>           ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!>           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
!>           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
!>           DIAGONAL OF THE EPSILON TABLE IS DELETED.
!>
      SUBROUTINE DQELG (N, EPSTAB, RESULT, ABSERR, RES3LA, NRES)
      use const, only: max_epstab
      real(kind=cfp) ABSERR, DELTA1, DELTA2, DELTA3, EPMACH, EPSINF, EPSTAB, ERROR, ERR1, ERR2, ERR3, E0, E1, E1ABS
      real(kind=cfp) E2, E3, OFLOW, RES, RESULT, RES3LA, SS, TOL1, TOL2, TOL3, cfp_dummy
      INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
      DIMENSION EPSTAB(max_epstab),RES3LA(3)
!***FIRST EXECUTABLE STATEMENT  DQELG
      EPMACH = F1MACH(4,cfp_dummy)
      OFLOW = F1MACH(2,cfp_dummy)
      NRES = NRES+1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF(N.LT.3) GO TO 100
      LIMEXP = max_epstab-2
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 40 I = 1,NEWELM
        K2 = K1-1
        K3 = K1-2
        RES = EPSTAB(K1+2)
        E0 = EPSTAB(K3)
        E1 = EPSTAB(K2)
        E2 = RES
        E1ABS = ABS(E1)
        DELTA2 = E2-E1
        ERR2 = ABS(DELTA2)
        TOL2 = MAX(ABS(E2),E1ABS)*EPMACH
        DELTA3 = E1-E0
        ERR3 = ABS(DELTA3)
        TOL3 = MAX(E1ABS,ABS(E0))*EPMACH
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
!
!           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
!           ACCURACY, CONVERGENCE IS ASSUMED.
!           RESULT = E2
!           ABSERR = ABS(E1-E0)+ABS(E2-E1)
!
        RESULT = RES
        ABSERR = ERR2+ERR3
! ***JUMP OUT OF DO-LOOP
        GO TO 100
   10   E3 = EPSTAB(K1)
        EPSTAB(K1) = E1
        DELTA1 = E1-E3
        ERR1 = ABS(DELTA1)
        TOL1 = MAX(E1ABS,ABS(E3))*EPMACH
!
!           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
!           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
!
        IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
        SS = 0.1E+01_cfp/DELTA1+0.1E+01_cfp/DELTA2-0.1E+01_cfp/DELTA3
        EPSINF = ABS(SS*E1)
!
!           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
!           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
!           OF N.
!
        IF(EPSINF.GT.0.1E-03_cfp) GO TO 30
   20   N = I+I-1
! ***JUMP OUT OF DO-LOOP
        GO TO 50
!
!           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
!           THE VALUE OF RESULT.
!
   30   RES = E1+0.1E+01_cfp/SS
        EPSTAB(K1) = RES
        K1 = K1-2
        ERROR = ERR2+ABS(RES-E2)+ERR3
        IF(ERROR.GT.ABSERR) GO TO 40
        ABSERR = ERROR
        RESULT = RES
   40 CONTINUE
!
!           SHIFT THE TABLE.
!
   50 IF(N.EQ.LIMEXP) N = 2*(LIMEXP/2)-1
      IB = 1
      IF((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM+1
      DO 60 I=1,IE
        IB2 = IB+2
        EPSTAB(IB) = EPSTAB(IB2)
        IB = IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      INDX = NUM-N+1
      DO 70 I = 1,N
        EPSTAB(I)= EPSTAB(INDX)
        INDX = INDX+1
   70 CONTINUE
   80 IF(NRES.GE.4) GO TO 90
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 100
!
!           COMPUTE ERROR ESTIMATE
!
   90 ABSERR = ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2))+ABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  100 ABSERR = MAX(ABSERR,0.5E+01_cfp*EPMACH*ABS(RESULT))
      RETURN
      END SUBROUTINE DQELG

!>***BEGIN PROLOGUE  DQK21
!>***PURPOSE  To compute I = Integral of F over (A,B), with error
!>                           estimate
!>                       J = Integral of ABS(F) over (A,B)
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A2
!>***TYPE      real(kind=cfp) (QK21-S, DQK21-D)
!>***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Integration rules
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS
!>            ON ENTRY
!>            F      - class(bound_user_function)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>              A      - Double precision
!>                       Lower limit of integration
!>
!>              B      - Double precision
!>                       Upper limit of integration
!>
!>            ON RETURN
!>              RESULT - Double precision
!>                       Approximation to the integral I
!>                       RESULT is computed by applying the 21-POINT
!>                       KRONROD RULE (RESK) obtained by optimal addition
!>                       of abscissae to the 10-POINT GAUSS RULE (RESG).
!>
!>              ABSERR - Double precision
!>                       Estimate of the modulus of the absolute error,
!>                       which should not exceed ABS(I-RESULT)
!>
!>              RESABS - Double precision
!>                       Approximation to the integral J
!>
!>              RESASC - Double precision
!>                       Approximation to the integral of ABS(F-I/(B-A))
!>                       over (A,B)
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  DQK21
!>
!>
!>           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!>           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!>           CORRESPONDING WEIGHTS ARE GIVEN.
!>
!>           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
!>                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
!>                    GAUSS RULE
!>                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!>                    ADDED TO THE 10-POINT GAUSS RULE
!>
!>           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
!>
!>           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
!>
!>
!> GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
!> AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
!> BELL LABS, NOV. 1981.
!>
!>
!>           LIST OF MAJOR VARIABLES
!>           -----------------------
!>
!>           CENTR  - MID POINT OF THE INTERVAL
!>           HLGTH  - HALF-LENGTH OF THE INTERVAL
!>           ABSC   - ABSCISSA
!>           FVAL*  - FUNCTION VALUE
!>           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
!>           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
!>           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!>                    I.E. TO I/(B-A)
!>
!>
!>           MACHINE DEPENDENT CONSTANTS
!>           ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!>
      SUBROUTINE DQK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
      class(bound_user_function) :: F
      real(kind=cfp) A, ABSC, ABSERR, B, CENTR, DHLGTH, EPMACH, FC, FSUM, FVAL1, FVAL2, FV1, FV2, HLGTH, RESABS, RESASC
      real(kind=cfp) RESG, RESK, RESKH, RESULT, UFLOW, WG, WGK, XGK, cfp_dummy
      INTEGER J,JTW,JTWM1
!
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
!
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.066671344308688137593568809893332_cfp /
      DATA WG  (  2) / 0.149451349150580593145776339657697_cfp /
      DATA WG  (  3) / 0.219086362515982043995534934228163_cfp /
      DATA WG  (  4) / 0.269266719309996355091226921569469_cfp /
      DATA WG  (  5) / 0.295524224714752870173892994651338_cfp /
!
      DATA XGK (  1) / 0.995657163025808080735527280689003_cfp /
      DATA XGK (  2) / 0.973906528517171720077964012084452_cfp /
      DATA XGK (  3) / 0.930157491355708226001207180059508_cfp /
      DATA XGK (  4) / 0.865063366688984510732096688423493_cfp /
      DATA XGK (  5) / 0.780817726586416897063717578345042_cfp /
      DATA XGK (  6) / 0.679409568299024406234327365114874_cfp /
      DATA XGK (  7) / 0.562757134668604683339000099272694_cfp /
      DATA XGK (  8) / 0.433395394129247190799265943165784_cfp /
      DATA XGK (  9) / 0.294392862701460198131126603103866_cfp /
      DATA XGK ( 10) / 0.148874338981631210884826001129720_cfp /
      DATA XGK ( 11) / 0.000000000000000000000000000000000_cfp /
!
      DATA WGK (  1) / 0.011694638867371874278064396062192_cfp /
      DATA WGK (  2) / 0.032558162307964727478818972459390_cfp /
      DATA WGK (  3) / 0.054755896574351996031381300244580_cfp /
      DATA WGK (  4) / 0.075039674810919952767043140916190_cfp /
      DATA WGK (  5) / 0.093125454583697605535065465083366_cfp /
      DATA WGK (  6) / 0.109387158802297641899210590325805_cfp /
      DATA WGK (  7) / 0.123491976262065851077958109831074_cfp /
      DATA WGK (  8) / 0.134709217311473325928054001771707_cfp /
      DATA WGK (  9) / 0.142775938577060080797094273138717_cfp /
      DATA WGK ( 10) / 0.147739104901338491374841515972068_cfp /
      DATA WGK ( 11) / 0.149445554002916905664936468389821_cfp /
!
!***FIRST EXECUTABLE STATEMENT  DQK21
      EPMACH = F1MACH(4,cfp_dummy)
      UFLOW = F1MACH(1,cfp_dummy)
!
      CENTR = 0.5_cfp*(A+B)
      HLGTH = 0.5_cfp*(B-A)
      DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
      RESG = 0.0_cfp
      FC = F%eval(CENTR)
      RESK = WGK(11)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F%eval(CENTR-ABSC)
        FVAL2 = F%eval(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,5
        JTWM1 = 2*J-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F%eval(CENTR-ABSC)
        FVAL2 = F%eval(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5_cfp
      RESASC = WGK(11)*ABS(FC-RESKH)
      DO 20 J=1,10
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0_cfp.AND.ABSERR.NE.0.0_cfp) ABSERR = RESASC*MIN(0.1E+01_cfp,(0.2E+03_cfp*ABSERR/RESASC)**1.5E+00_cfp)
      IF(RESABS.GT.UFLOW/(0.5E+02_cfp*EPMACH)) ABSERR = MAX((EPMACH*0.5E+02_cfp)*RESABS,ABSERR)
      RETURN
      END SUBROUTINE DQK21

!>***BEGIN PROLOGUE  DQPSRT
!>***SUBSIDIARY
!>***PURPOSE  This routine maintains the descending ordering in the
!>            list of the local error estimated resulting from the
!>            interval subdivision process. At each call two error
!>            estimates are inserted using the sequential search
!>            method, top-down for the largest error estimate and
!>            bottom-up for the smallest error estimate.
!>***LIBRARY   SLATEC
!>***TYPE      real(kind=cfp) (QPSRT-S, DQPSRT-D)
!>***KEYWORDS  SEQUENTIAL SORTING
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Ordering routine
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS (MEANING AT OUTPUT)
!>              LIMIT  - Integer
!>                       Maximum number of error estimates the list
!>                       can contain
!>
!>              LAST   - Integer
!>                       Number of error estimates currently in the list
!>
!>              MAXERR - Integer
!>                       MAXERR points to the NRMAX-th largest error
!>                       estimate currently in the list
!>
!>              ERMAX  - Double precision
!>                       NRMAX-th largest error estimate
!>                       ERMAX = ELIST(MAXERR)
!>
!>              ELIST  - Double precision
!>                       Vector of dimension LAST containing
!>                       the error estimates
!>
!>              IORD   - Integer
!>                       Vector of dimension LAST, the first K elements
!>                       of which contain pointers to the error
!>                       estimates, such that
!>                       ELIST(IORD(1)),...,  ELIST(IORD(K))
!>                       form a decreasing sequence, with
!>                       K = LAST if LAST.LE.(LIMIT/2+2), and
!>                       K = LIMIT+1-LAST otherwise
!>
!>              NRMAX  - Integer
!>                       MAXERR = IORD(NRMAX)
!>
!>***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE
!>***ROUTINES CALLED  (NONE)
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>***END PROLOGUE  DQPSRT
!>
      SUBROUTINE DQPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
      real(kind=cfp) ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,NRMAX
      DIMENSION ELIST(*),IORD(*)
!
!           CHECK WHETHER THE LIST CONTAINS MORE THAN
!           TWO ERROR ESTIMATES.
!
!***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
!
!           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
!           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
!           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
!           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
!
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
!
!           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
!           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
!           SUBDIVISIONS STILL ALLOWED.
!
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
!
!           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
!           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
!
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
!
!           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
!
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
!
!           SET MAXERR AND ERMAX.
!
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END SUBROUTINE DQPSRT

end module general_quadrature
