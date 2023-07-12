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
module precisn
! Definition of precision-related parameters
! from Martin Plummer, Nov 2009
  implicit none
  public

! Uncomment if needed
!  logical, parameter    :: real64 = .TRUE.
!  integer, parameter    :: ibyte = 4
  integer, parameter    :: rbyte = 8
!  integer, parameter    :: zbyte = 16

!  Set kind type parameters.
!  The argument to the function `selected_real_kind' is the 
!  requested minimum number of decimal digits of accuracy.
  integer, parameter    :: sp=selected_real_kind(6)   ! `single' precision 
  integer, parameter    :: sp_bytes = 4               ! number of bytes in float of type real(kind=sp)
  integer, parameter    :: wp=selected_real_kind(12)  ! `double' precision
  integer, parameter    :: wp_bytes = 8               ! number of bytes in float of type real(kind=wp)
  integer, parameter    :: ep1=selected_real_kind(19) ! `quad', extended, precision
  integer, parameter    :: ep1_bytes = 16             ! number of bytes in float of type real(kind=ep1)
  integer, parameter    :: ep = MAX(ep1, wp) ! extended precision if possible; double precision if not

! The use of `selected_int_kind' below is to ensure that the corresponding integer kind type 
! parameters can support integers up to at least 10^9 and 10^10, respectively.
  integer, parameter    :: shortint=selected_int_kind(9) ! i.e. a 32-bit signed integer
  integer, parameter    :: longint=selected_int_kind(10) ! i.e. a 64-bit signed integer

! Uncomment if needed
!  real(wp), parameter   :: acc8 = 2.0e-16_wp
!  real(wp), parameter   :: acc16 = 3.0e-33_ep
!  real(wp), parameter   :: fpmax = 1.0e60_wp
!  real(wp), parameter   :: fpmin = 1.0e-60_wp
!  real(wp), parameter   :: fplmax = 140.0_wp
!  real(wp), parameter   :: fplmin = -140.0_wp

  !> Floating point precision for which the library is to be compiled (Current-Float-Precision). It must be one of: wp, ep, ep1. Single precision (sp) is not supported.
  !> \warning Under no circumstances chnage or swap the definitions of sp,wp,ep,ep1 - these must correspond to the single,double and extended (quadruple) precision kind parameters.
#ifdef usequadprec
  integer, parameter :: cfp = ep
#else
  integer, parameter :: cfp = wp
#endif
  !> Number of bytes in float of type real(kind=cfp).
  !> \warning This number is NOT the same as wp,ep1 contrary to popular belief! The number of bytes used to represent the data type with the selected kind value is compiler-specific.
#ifdef usequadprec
  integer, parameter :: cfp_bytes = ep1_bytes
#else
  integer, parameter :: cfp_bytes = wp_bytes
#endif

  !> Storage unit used when writing values of kind cfp into disk.
  integer, parameter :: storage_unit_cfp = cfp_bytes

  !> Storage unit used when writing values of default integer kind into disk.
  integer, parameter :: storage_unit_int = bit_size(wp)/8

  !> Dummy variable used on input to F1MACH in various routines.
  real(kind=cfp), parameter :: cfp_dummy = 1.0_cfp

  !> Dummy variable used on input to F1MACH in various routines.
  real(kind=wp), parameter :: wp_dummy = 1.0_cfp

  !> Dummy variable used on input to F1MACH in various routines.
  real(kind=ep1), parameter :: ep_dummy = 1.0_cfp

  !> This function resolves into the particular routines depending on the type of the second dummy argument 'p' which specifies the data type for which the machine parameters are required.
  interface f1mach
     module procedure r1mach, d1mach, q1mach
  end interface

  public i1mach, print_precision_params
  private r1mach, d1mach, q1mach

  contains

  subroutine print_precision_params(stdout)
#ifdef usempi
     use mpi
#endif
     implicit none
     integer, intent(in) :: stdout 
     real(kind=cfp) :: test
#ifdef usempi
     integer(kind=kind(MPI_COMM_WORLD)) :: test_i  ! kind of this integer must match "mpiint" from mpi_mod
#endif
        write(stdout,*) 'Smallest positive magnitude:',F1MACH(1,test)
        write(stdout,*) 'Largest magnitude:',F1MACH(2,test)
        write(stdout,*) 'Smallest relative spacing:',F1MACH(3,test)
        write(stdout,*) 'Largest relative spacing:',F1MACH(4,test)
        write(stdout,*) 'Log10(Base)',F1MACH(5,test)
        write(stdout,*) 'Decimal precision for real(kind=cfp):',precision(test)
        write(stdout,*) 'Number of bits in default integer:',bit_size(cfp)
#ifdef usempi
        write(stdout,*) 'Number of bits in MPI integers:',bit_size(test_i)
#endif
  end subroutine print_precision_params

!> \verbatim
!>***PURPOSE  Return floating point machine dependent constants.
!>***LIBRARY   SLATEC
!>***CATEGORY  R1
!>***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!>***KEYWORDS  MACHINE CONSTANTS
!>***AUTHOR  Fox, P. A., (Bell Labs)
!>           Hall, A. D., (Bell Labs)
!>           Schryer, N. L., (Bell Labs)
!>***DESCRIPTION
!>
!>   D1MACH can be used to obtain machine-dependent parameters for the
!>   local machine environment.  It is a function subprogram with one
!>   (input) argument, and can be referenced as follows:
!>
!>        D = D1MACH(I)
!>
!>   where I=1,...,5.  The (output) value of D above is determined by
!>   the (input) value of I.  The results for various values of I are
!>   discussed below.
!>
!>   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!>   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!>   D1MACH( 3) = B**(-T), the smallest relative spacing.
!>   D1MACH( 4) = B**(1-T), the largest relative spacing.
!>   D1MACH( 5) = LOG10(B)
!>
!>***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!>                 a portable library, ACM Transactions on Mathematical
!>                 Software 4, 2 (June 1978), pp. 177-188.
!> \endverbatim
!
      FUNCTION D1MACH (I,p)
!
!>
!>
      INTEGER :: I
      REAL(kind=wp) :: X,B, D1MACH, p
!>
!>***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) STOP "PRECISN/D1MACH: I OUT OF BOUNDS." !CALL XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
!>
      X = 1.0_wp
      B = RADIX(X)
      SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)
        CASE (2)
          D1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          D1MACH = EPSILON(X)            ! the largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
          STOP
      END SELECT
!>
      END FUNCTION

!> \verbatim
!>***PURPOSE  Return floating point machine dependent constants.
!>***LIBRARY   SLATEC
!>***CATEGORY  R1
!>***TYPE      DOUBLE PRECISION (R1MACH-S, Q1MACH-D)
!>***KEYWORDS  MACHINE CONSTANTS
!>***AUTHOR  Fox, P. A., (Bell Labs)
!>           Hall, A. D., (Bell Labs)
!>           Schryer, N. L., (Bell Labs)
!>***DESCRIPTION
!>
!>   Q1MACH can be used to obtain machine-dependent parameters for the
!>   local machine environment.  It is a function subprogram with one
!>   (input) argument, and can be referenced as follows:
!>
!>        D = Q1MACH(I)
!>
!>   where I=1,...,5.  The (output) value of D above is determined by
!>   the (input) value of I.  The results for various values of I are
!>   discussed below.
!>
!>   Q1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!>   Q1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!>   Q1MACH( 3) = B**(-T), the smallest relative spacing.
!>   Q1MACH( 4) = B**(1-T), the largest relative spacing.
!>   Q1MACH( 5) = LOG10(B)
!>
!>***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!>                 a portable library, ACM Transactions on Mathematical
!>                 Software 4, 2 (June 1978), pp. 177-188.
!> \endverbatim
!
      FUNCTION Q1MACH (I,p)
!
!
!
      REAL(kind=ep1) :: X,B, Q1MACH, p
      INTEGER :: I
!
!***FIRST EXECUTABLE STATEMENT  Q1MACH
      IF (I .LT. 1 .OR. I .GT. 5) STOP "PRECISN/Q1MACH: I OUT OF BOUNDS." !CALL XERMSG ('SLATEC', 'Q1MACH', 'I OUT OF BOUNDS', 1, 2)
!
      X = 1.0_ep1
      B = RADIX(X)
      SELECT CASE (I)
        CASE (1)
          Q1MACH = TINY(X)
        CASE (2)
          Q1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          Q1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          Q1MACH = EPSILON(X)            ! the largest relative spacing.
        CASE (5)
          Q1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN Q1MACH - I OUT OF BOUNDS')
          STOP
      END SELECT
!
      END FUNCTION
!> \verbatim
!>***PURPOSE  Return integer machine dependent constants.
!>***LIBRARY   SLATEC
!>***CATEGORY  R1
!>***TYPE      INTEGER (I1MACH-I)
!>***KEYWORDS  MACHINE CONSTANTS
!>***AUTHOR  Fox, P. A., (Bell Labs)
!>           Hall, A. D., (Bell Labs)
!>           Schryer, N. L., (Bell Labs)
!>***DESCRIPTION
!>
!>   I1MACH can be used to obtain machine-dependent parameters for the
!>   local machine environment.  It is a function subprogram with one
!>   (input) argument and can be referenced as follows:
!>
!>        K = I1MACH(I)
!>
!>   where I=1,...,16.  The (output) value of K above is determined by
!>   the (input) value of I.  The results for various values of I are
!>   discussed below.
!>
!>   I/O unit numbers:
!>     I1MACH( 1) = the standard input unit.
!>     I1MACH( 2) = the standard output unit.
!>     I1MACH( 3) = the standard punch unit.
!>     I1MACH( 4) = the standard error message unit.
!>
!>   Words:
!>     I1MACH( 5) = the number of bits per integer storage unit.
!>     I1MACH( 6) = the number of characters per integer storage unit.
!>
!>   Integers:
!>     assume integers are represented in the S-digit, base-A form
!>
!>                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!>
!>                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!>     I1MACH( 7) = A, the base.
!>     I1MACH( 8) = S, the number of base-A digits.
!>     I1MACH( 9) = A**S - 1, the largest magnitude.
!>
!>   Floating-Point Numbers:
!>     Assume floating-point numbers are represented in the T-digit,
!>     base-B form
!>                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!>
!>                where 0 .LE. X(I) .LT. B for I=1,...,T,
!>                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!>     I1MACH(10) = B, the base.
!>
!>   Single-Precision:
!>     I1MACH(11) = T, the number of base-B digits.
!>     I1MACH(12) = EMIN, the smallest exponent E.
!>     I1MACH(13) = EMAX, the largest exponent E.
!>
!>   Double-Precision:
!>     I1MACH(14) = T, the number of base-B digits.
!>     I1MACH(15) = EMIN, the smallest exponent E.
!>     I1MACH(16) = EMAX, the largest exponent E.
!>
!>   Quad-Precision:
!>     I1MACH(17) = T, the number of base-B digits.
!>     I1MACH(18) = EMIN, the smallest exponent E.
!>     I1MACH(19) = EMAX, the largest exponent E.
!>
!>***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!>                 a portable library, ACM Transactions on Mathematical
!>                 Software 4, 2 (June 1978), pp. 177-188.
!>***ROUTINES CALLED  (NONE)
!> \endverbatim
!
      INTEGER FUNCTION I1MACH (I)
      use iso_fortran_env, only: error_unit, input_unit, output_unit, character_storage_size
!
!
      INTEGER :: I
      REAL(kind=sp) :: x_sp
      REAL(kind=wp) :: x_wp
      REAL(kind=ep1) :: x_ep
!
!***FIRST EXECUTABLE STATEMENT  I1MACH
      SELECT CASE (I)
        CASE (1)
          I1MACH = input_unit
        CASE (2)
          I1MACH = output_unit
        CASE (3)
          I1MACH = output_unit
        CASE (4)
          I1MACH = error_unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          I1MACH = character_storage_size
        CASE (7)
          I1MACH = RADIX(I)
        CASE (8)
          I1MACH = DIGITS(I)
        CASE (9)
          I1MACH = HUGE(I)
        CASE (10)
          I1MACH = RADIX(x_sp)
        CASE (11)
          I1MACH = DIGITS(x_sp)
        CASE (12)
          I1MACH = MINEXPONENT(x_sp)
        CASE (13)
          I1MACH = MAXEXPONENT(x_sp)
        CASE (14)
          I1MACH = DIGITS(x_wp)
        CASE (15)
          I1MACH = MINEXPONENT(x_wp)
        CASE (16)
          I1MACH = MAXEXPONENT(x_wp)
        CASE (17)
          I1MACH = DIGITS(x_ep)
        CASE (18)
          I1MACH = MINEXPONENT(x_ep)
        CASE (19)
          I1MACH = MAXEXPONENT(x_ep)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
          STOP
      END SELECT
      END FUNCTION
!> \verbatim
!>***PURPOSE  Return floating point machine dependent constants.
!>***LIBRARY   SLATEC
!>***CATEGORY  R1
!>***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!>***KEYWORDS  MACHINE CONSTANTS
!>***AUTHOR  Fox, P. A., (Bell Labs)
!>           Hall, A. D., (Bell Labs)
!>           Schryer, N. L., (Bell Labs)
!>***DESCRIPTION
!>
!>   R1MACH can be used to obtain machine-dependent parameters for the
!>   local machine environment.  It is a function subprogram with one
!>   (input) argument, and can be referenced as follows:
!>
!>        A = R1MACH(I)
!>
!>   where I=1,...,5.  The (output) value of A above is determined by
!>   the (input) value of I.  The results for various values of I are
!>   discussed below.
!>
!>   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!>   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!>   R1MACH(3) = B**(-T), the smallest relative spacing.
!>   R1MACH(4) = B**(1-T), the largest relative spacing.
!>   R1MACH(5) = LOG10(B)
!>
!>***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!>                 a portable library, ACM Transactions on Mathematical
!>                 Software 4, 2 (June 1978), pp. 177-188.
!> \endverbatim
!
      FUNCTION R1MACH (I,p)
!
!
!
      REAL(kind=sp) :: X,B,p, R1MACH
      INTEGER :: I
!
!
!
!***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1 .OR. I .GT. 5) STOP "PRECISN/R1MACH: I OUT OF BOUNDS." !CALL XERMSG ('SLATEC', 'R1MACH', 'I OUT OF BOUNDS', 1, 2)
!
      X = 1.0_sp
      B = RADIX(X)
      SELECT CASE (I)
        CASE (1)
          R1MACH = TINY(X)
        CASE (2)
          R1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          R1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          R1MACH = EPSILON(X)            ! the largest relative spacing.
        CASE (5)
          R1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN R1MACH - I OUT OF BOUNDS')
          STOP
      END SELECT
!
      END FUNCTION

end module precisn
