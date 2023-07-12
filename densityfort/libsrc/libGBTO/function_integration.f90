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
!> Declares functions that can be integrated (together with the B-spline) using MODDBFQAD
module function_integration
 use general_quadrature, only: bound_user_function
 use quadrature_module
 use precisn

 implicit none

 !> \f$ r^{l} \f$
 type, extends(bound_user_function) :: power_function
      real(kind=cfp) :: l !power
 contains
      !> \memberof power_function
      procedure :: wp_eval => wp_power_function_evaluation
      !> \memberof power_function
      procedure :: ep_eval => ep_power_function_evaluation
 end type 

 private wp_power_function_evaluation, ep_power_function_evaluation

 !> \f$ r^{l}\exp[ar^{2}] \f$
 type, extends(bound_user_function) :: poly_exp_function
      real(kind=cfp) :: l !power
      real(kind=cfp) :: a !exponent
 contains
      !> \memberof poly_exp_function
      procedure :: wp_eval => wp_poly_exp_evaluation
      !> \memberof poly_exp_function
      procedure :: ep_eval => ep_poly_exp_evaluation
 end type

 private wp_poly_exp_evaluation, ep_poly_exp_evaluation

 !> \f$ B_{i}(r)r^l \f$
 !> This function is used in calculation of the overlap integral between two B-spline functions.
 type, extends(bound_user_function) :: bspl_prod_pow
      integer :: l   !power of the radial coordinate
      integer :: i   !index of the radial B-spline function
      integer, private :: i_ref !index of the radial B-spline function: this is a reference value used in 'eval' to check whether the index 'i' which is accessible by the user has changed or not.
      integer :: id  !derivative of the radial B-spline function. id = 0: no derivative
      integer :: k !order of the B-splines
      integer :: n !number of B-splines in the basis
      real(kind=wp), allocatable, private :: wp_knots(:) !knot sequence (1:n+k)
      real(kind=wp), allocatable, private :: wp_bcoef(:) !coefficient array for the B-spline
      real(kind=ep1), allocatable, private :: ep_knots(:) !knot sequence (1:n+k)
      real(kind=ep1), allocatable, private :: ep_bcoef(:) !coefficient array for the B-spline
      integer, private :: inbv !helper variable for BVALU.
      logical, private :: initialized = .false. !this is set to 1 after initialization
      logical, private :: wp_input = .false.
      logical, private :: ep_input = .false.
      real(kind=wp), allocatable, private :: wp_work(:) !work array for BVALU
      real(kind=ep1), allocatable, private :: ep_work(:) !work array for BVALU
 contains
      !> \memberof bspl_prod_pow_evaluation
      !> This function takes as an input an array of knots. The size of the input array must be k+n. The values k, n must be set before calling this function.
      !> The input array of knots is transferred to the internal array knots which is then used in eval.
      procedure :: init => init_bspl_prod_pow
      !> \memberof bspl_prod_pow_evaluation
      procedure :: wp_eval => wp_bspl_prod_pow_evaluation
      !> \memberof bspl_prod_pow_evaluation
      procedure :: ep_eval => ep_bspl_prod_pow_evaluation
 end type

 private wp_bspl_prod_pow_evaluation, ep_bspl_prod_pow_evaluation

 !> A complicated radial function
 type, extends(bound_user_function) :: radial_function
      integer :: li      !angular momentum of the B-spline orbital
      integer :: l1      !angular momentum of the spherical harmonic from the translation of the spherical harmonic of the GTO
      integer :: l       !angular momentum of the partial wave in the SCE of the GTO
      real(kind=cfp) :: alpha   !exponent of the GTO
      real(kind=cfp) :: RA      !distance of the atom from the CMS
      real(kind=cfp) :: NB1_li  !Norm of the B-spline orbital in the form: NB1_li = N^(1/li), where N is the normalization factor of the B-spline orbital
      real(kind=cfp) :: NGTO_ln !Norm of the GTO in the form: NGTO_ln = ln(N), where N is the normalization factor of the GTO. N = sqrt(2*(2*alpha)^(L+3/2)/Gamma(L+3/2))
 contains
      !> \memberof radial_function
      procedure :: wp_eval => wp_radial_evaluation
      !> \memberof radial_function
      procedure :: ep_eval => ep_radial_evaluation
 end type

 private wp_radial_evaluation, ep_radial_evaluation

 !> A complicated radial function that uses a supplied array to store the values of the Bessel functions. Note that this function is identical to radial_function except of 1.5_wp in one of the exponents.
 type, extends(bound_user_function) :: radial_function_buff
      integer :: li      !angular momentum of the B-spline orbital
      integer :: l1      !angular momentum of the spherical harmonic from the translation of the spherical harmonic of the GTO
      integer :: l       !angular momentum of the partial wave in the SCE of the GTO
      real(kind=cfp) :: alpha   !exponent of the GTO
      real(kind=cfp) :: RA      !distance of the atom from the CMS
      real(kind=wp), allocatable :: wp_bes(:)  !array (1:l+1) that will be used to store the values of the bessel functions in double precision
      real(kind=ep1), allocatable :: ep_bes(:)  !array (1:l+1) that will be used to store the values of the bessel functions in double precision
      real(kind=cfp), private :: f1mach_saved !value of f1mach(4)
      logical, private :: initialized = .false. !this variable is set to 1 after initialization. It is 0 otherwise
      logical, private :: wp_input = .false. !tells me whether cfp=wp
      logical, private :: ep_input = .false. !tells me whether cfp=ep
      integer :: no_eval = 0 !total number of calls to %eval.
 contains
      !> allocates the buffer
      !> \memberof radial_function_buff
      procedure :: init => init_r_f_buff
      !> deallocates the buffer
      !> \memberof radial_function_buff
      procedure :: final => final_r_f_buff
      !> \memberof radial_function_buff
      procedure :: wp_eval => wp_radial_evaluation_buff
      !> \memberof radial_function_buff
      procedure :: ep_eval => ep_radial_evaluation_buff
 end type

 private wp_radial_evaluation_buff, ep_radial_evaluation_buff

 !> Bessel function \f$J_{l+1/2}(kr)\f$, of half integral order \f$l\f$ multiplied by \f$r^{p}\f$: \f$r^{p}J_{l+1/2}(kr)\f$.
 !> This function can be used to calculate representation of the Bessel orbital in terms of the B-spline orbitals.
 type, extends(bound_user_function) :: bessel_fn
    !> Angular momentum of the Bessel function.
    integer :: l
    !> Linear momentum.
    real(kind=cfp) :: k
    !> Power of the radial coordinate.
    integer :: p
 contains
    !> Bessel function \f$J_{l+1/2}(kr)\f$ at \f$r\f$.
    procedure :: wp_eval => wp_bessel_eval
    !> Bessel function \f$J_{l+1/2}(kr)\f$ at \f$r\f$.
    procedure :: ep_eval => ep_bessel_eval
 end type

 private wp_bessel_eval, ep_bessel_eval

contains
 
 real(wp) function wp_power_function_evaluation(data,x) result(r)
      implicit none
      class(power_function) :: data
      real(wp), intent(in) :: x

         r = x**(data%l)

 end function wp_power_function_evaluation

 real(ep1) function ep_power_function_evaluation(data,x) result(r)
      implicit none
      class(power_function) :: data
      real(ep1), intent(in) :: x

         r = x**(data%l)

 end function ep_power_function_evaluation

 real(wp) function wp_poly_exp_evaluation(data,x) result(r)
      implicit none
      class(poly_exp_function) :: data
      real(wp), intent(in) :: x

         r = x**(data%l)*exp((data%a)*x**2.0_wp)

 end function wp_poly_exp_evaluation

 real(ep1) function ep_poly_exp_evaluation(data,x) result(r)
      implicit none
      class(poly_exp_function) :: data
      real(ep1), intent(in) :: x

         r = x**(data%l)*exp((data%a)*x**2.0_ep1)

 end function ep_poly_exp_evaluation

 real(wp) function wp_radial_evaluation(data,x) result(r)
      use special_functions, only: cfp_besi
      use utils, only: xermsg
      implicit none
      class(radial_function) :: data
      real(wp), intent(in) :: x
      real(wp) :: arg
      integer, parameter :: kode = 2
      real(wp), parameter :: half = 0.5_wp
      integer :: nz, err
      real(wp), allocatable :: y(:)

         allocate(y(1:data%l+1),stat=err)
         if (err /= 0) then
            call xermsg ('function_integration', 'radial_evaluation', 'Memory allocation failed; see radial_evaluation', err, 1)
         end if
   
         arg = 2.0_wp*data%alpha*x*data%RA
         CALL cfp_besi(arg, half, kode, data%l+1, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
         
         r = (data%NB1_li*x)**(data%li)*x**(data%l1+2.0_wp)*exp(-(data%alpha)*(x-(data%RA))**2.0_wp + data%NGTO_ln)*y(data%l+1)
     
         deallocate(y)

 end function wp_radial_evaluation

 real(ep1) function ep_radial_evaluation(data,x) result(r)
      use special_functions, only: cfp_besi
      use utils, only: xermsg
      implicit none
      class(radial_function) :: data
      real(ep1), intent(in) :: x
      real(ep1) :: arg
      integer, parameter :: kode = 2
      real(ep1), parameter :: half = 0.5_ep1
      integer :: nz, err
      real(ep1), allocatable :: y(:)

         allocate(y(1:data%l+1),stat=err)
         if (err /= 0) then
            call xermsg ('function_integration', 'radial_evaluation', 'Memory allocation failed; see radial_evaluation', err, 1)
         end if

         arg = 2.0_ep1*data%alpha*x*data%RA
         CALL cfp_besi(arg, half, kode, data%l+1, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.

         r = (data%NB1_li*x)**(data%li)*x**(data%l1+2.0_ep1)*exp(-(data%alpha)*(x-(data%RA))**2.0_ep1 + data%NGTO_ln)*y(data%l+1)

         deallocate(y)

 end function ep_radial_evaluation

 integer function init_r_f_buff(data) result(err)
      use utils, only: xermsg
      implicit none
      class(radial_function_buff) :: data
      real(kind=cfp) :: cfp_dummy

         err = 0
   
         if (data%l < 0) then !\todo: add error checking for the other parameters as well
            err = 1
            call xermsg ('function_integration', 'init_r_f_buff', 'Negative value of "l"; see init_r_f_buff', err, 1)
         endif

         !Determine which buffer will keep the values of the Bessel functions
         data%wp_input = .false.; data%ep_input = .false.
         if (cfp .eq. wp) then
            data%wp_input = .true.
            data%ep_input = .false.
         elseif (cfp .eq. ep1) then
            data%wp_input = .false.
            data%ep_input = .true.
         else
            call xermsg ('function_integration', 'init_r_f_buff', 'Unsupported numberic type on input; see init_r_f_buff', 1, 1)
         endif

         if (allocated(data%wp_bes)) deallocate(data%wp_bes)
         if (allocated(data%ep_bes)) deallocate(data%ep_bes)
   
         err = -1
         if (data%wp_input) then
            allocate(data%wp_bes(1:data%l+1),stat=err)
         elseif (data%ep_input) then
            allocate(data%ep_bes(1:data%l+1),stat=err)
         endif
         if (err /= 0) call xermsg ('function_integration', 'init_r_f_buff', 'Memory allocation failed; see init_r_f_buff', err, 1)

         data%f1mach_saved = f1mach(4,cfp_dummy)

         data%no_eval = 0

         data%initialized = .true.

 end function init_r_f_buff

 integer function final_r_f_buff(data) result(err)
      use utils, only: xermsg
      implicit none
      class(radial_function_buff) :: data

         err = 0
   
         if (allocated(data%wp_bes)) deallocate(data%wp_bes,stat=err)
         if (allocated(data%ep_bes)) deallocate(data%ep_bes,stat=err)
         
         if (err /= 0) then
            call xermsg ('function_integration', 'final_r_f_buff', 'Memory deallocation failed; see final_r_f_buff', err, 1)
         end if
   
         data%no_eval = 0
         data%initialized = .false.
         data%wp_input = .false.; data%ep_input = .false. 

 end function final_r_f_buff

 real(wp) function wp_radial_evaluation_buff(data,x) result(r)
      use special_functions, only: cfp_besi
      use utils, only: xermsg
      implicit none
      class(radial_function_buff) :: data
      real(wp), intent(in) :: x
      real(wp) :: arg, exp_arg, pow
      integer, parameter :: kode = 2
      real(wp), parameter :: half = 0.5_wp
      integer :: nz, pow_exponent

         if (.not. data % initialized) then
            call xermsg ('function_integration', 'radial_evaluation_buff', 'Function not initialized.', 1, 1)
         end if

         arg = 2.0_wp*data%alpha*x*data%RA
         pow_exponent = data%l+1
         call cfp_besi(arg, half, kode, pow_exponent, data%wp_bes, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
   
         if (data % l1 + data % li + 1.5_wp < 0.0_wp .and. x <= data % f1mach_saved) then
            call xermsg ('function_integration', 'radial_evaluation_buff', &
                         'The integrand would evaluate to an inaccurate number.', 2, 1)
         end if

         exp_arg = x-data%RA
         exp_arg = exp_arg*exp_arg
         exp_arg = data%alpha*exp_arg
         !r = x**(data%l1+data%li+1.5_wp)*exp(-exp_arg)*data%bes(data%l+1)
         pow_exponent = data%l1+data%li+1
         pow = 1.0_wp
         do nz = 1,abs(pow_exponent) !this is faster than: pow=x**(data%l1+data%li+1)
            pow = pow*x
         enddo
         if (sign(1,data%l1+data%li+1) < 0) then 
            pow=sqrt(x)/pow
         else
            pow = pow*sqrt(x)
         endif
         r = pow*exp(-exp_arg)*data%wp_bes(data%l+1)

         data%no_eval = data%no_eval + 1
 
 end function wp_radial_evaluation_buff

 real(ep1) function ep_radial_evaluation_buff(data,x) result(r)
      use special_functions, only: cfp_besi !_direct
      use utils, only: xermsg
      implicit none
      class(radial_function_buff) :: data
      real(ep1), intent(in) :: x
      real(ep1) :: arg, exp_arg, pow
      integer, parameter :: kode = 2
      real(ep1), parameter :: half = 0.5_ep1
      integer :: nz, pow_exponent

         if (.not. data % initialized) then
            call xermsg ('function_integration', 'radial_evaluation_buff', 'Function not initialized.', 1, 1)
         end if

         arg = 2.0_ep1*data%alpha*x*data%RA
         pow_exponent = data%l+1
         call cfp_besi(arg, half, kode, pow_exponent, data%ep_bes, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.

         if (data % l1 + data % li + 1.5_ep1 < 0.0_ep1 .and. x <= data % f1mach_saved) then
            call xermsg ('function_integration', 'radial_evaluation_buff', &
                         'The integrand would evaluate to an inaccurate number.', 2, 1)
         end if

         exp_arg = x-data%RA
         exp_arg = exp_arg*exp_arg
         exp_arg = data%alpha*exp_arg
         !r = x**(data%l1+data%li+1.5_ep1)*exp(-exp_arg)*data%bes(data%l+1) !Direct power evaluation
         pow_exponent = data%l1+data%li+1
         pow = 1.0_ep1
         do nz = 1,abs(pow_exponent) !this is faster than: pow=x**(data%l1+data%li+1)
            pow = pow*x
         enddo
         if (sign(1,data%l1+data%li+1) < 0) then
            pow=sqrt(x)/pow
         else
            pow = pow*sqrt(x)
         endif
         r = pow*exp(-exp_arg)*data%ep_bes(data%l+1)

         data%no_eval = data%no_eval + 1
 
 end function ep_radial_evaluation_buff

 integer function init_bspl_prod_pow(data,knot_inp) result(e)
      use utils, only: xermsg
      implicit none
      class(bspl_prod_pow) :: data
      real(cfp), intent(in) :: knot_inp(:)
      integer :: err

         e = 0

         if (data%n .le. 0 .or. data%k .le. 0) then
            e = 1
            call xermsg ('function_integration', 'init_bspl_prod_pow', 'Invalid value(s) of n or k; see init_bspl_prod_pow', e, 1)
         endif

!         if (data%l < 0) then
!            e = 2
!            call xermsg ('function_integration', 'init_bspl_prod_pow', 'Invalid value of l; see init_bspl_prod_pow', e, 1)
!         endif

         if (data%id < 0) then
            e = 3
            call xermsg ('function_integration', 'init_bspl_prod_pow', 'Invalid value of id; see init_bspl_prod_pow', e, 1)
         endif

         if (data%i .le. 0 .or. data%i > data%n) then
            e = 4
            call xermsg ('function_integration', 'init_bspl_prod_pow', 'Invalid value of i; see init_bspl_prod_pow', e, 1)
         endif

         if (size(knot_inp) .eq. data%n+data%k) then
            !Here we have to allocate the auxiliary arrays of the correct type so that the %eval method always passes the appropriate arrays to the BVALU function.
            data%wp_input = .false.
            data%ep_input = .false.
            if (kind(knot_inp(1)) .eq. wp) then

               if (allocated(data%wp_knots)) deallocate(data%wp_knots)
               if (allocated(data%wp_work)) deallocate(data%wp_work)
               if (allocated(data%wp_bcoef)) deallocate(data%wp_bcoef)

               allocate(data%wp_knots(1:data%n+data%k),data%wp_work(1:3*data%k),data%wp_bcoef(1:data%n),stat=err)
               if (err /= 0) then
                  call xermsg ('function_integration', 'init_bspl_prod_pow', &
                               'Memory allocation failed; see init_bspl_prod_pow', err, 1)
               end if
   
               data%wp_knots = knot_inp !transfer the knot array
               data%wp_bcoef = 0.0_cfp
               data%wp_bcoef(data%i) = 1.0_cfp

               data%wp_input = .true.
            elseif (kind(knot_inp(1)) .eq. ep1) then

               if (allocated(data%ep_knots)) deallocate(data%ep_knots)
               if (allocated(data%ep_work)) deallocate(data%ep_work)
               if (allocated(data%ep_bcoef)) deallocate(data%ep_bcoef)

               allocate(data%ep_knots(1:data%n+data%k),data%ep_work(1:3*data%k),data%ep_bcoef(1:data%n),stat=err)
               if (err /= 0) then
                  call xermsg ('function_integration', 'init_bspl_prod_pow', &
                               'Memory allocation failed; see init_bspl_prod_pow', err, 1)
               end if
   
               data%ep_knots = knot_inp !transfer the knot array
               data%ep_bcoef = 0.0_cfp
               data%ep_bcoef(data%i) = 1.0_cfp

               data%ep_input = .true.
            else
               call xermsg ('function_integration', 'init_bspl_prod_pow', &
                            'The kind of the knot array on input can be only one of: wp, ep1.', 2, 1)
            endif

            data%inbv = 1
            data%i_ref = data%i
         else
            e = 5
            call xermsg ('function_integration', 'init_bspl_prod_pow', &
                         'Size of the knot array is not compatible with N, K set in bspl_prod_pow', 1, 1)
         endif

         data%initialized = .true.

 end function init_bspl_prod_pow

 real(wp) function wp_bspl_prod_pow_evaluation(data,x) result(r)
      use bspline_base, only: bvalu
      use utils, only: xermsg
      implicit none
      class(bspl_prod_pow) :: data
      real(wp), intent(in) :: x

         if (.not. data % initialized) then
            call xermsg ('function_integration', 'bspl_prod_pow_evaluation', 'Function not initialized.', 1, 1)
         end if

         !if the user has changed the required B-spline index then we have to delete the previous coefficient from the coefficient array and put the new index there
         if (data%i .ne. data%i_ref) then
            data%wp_bcoef(data%i_ref) = 0.0_wp
            data%wp_bcoef(data%i) = 1.0_wp
            data%i_ref = data%i
         endif

         r = BVALU(data%wp_knots,data%wp_bcoef,data%n,data%k,data%id,x,data%inbv,data%wp_work)
         r = r*x**data%l

 end function wp_bspl_prod_pow_evaluation

 real(ep1) function ep_bspl_prod_pow_evaluation(data,x) result(r)
      use bspline_base, only: bvalu
      use utils, only: xermsg
      implicit none
      class(bspl_prod_pow) :: data
      real(ep1), intent(in) :: x

         if (.not. data % initialized) then
            call xermsg ('function_integration', 'bspl_prod_pow_evaluation', 'Function not initialized.', 1, 1)
         end if

         !if the user has changed the required B-spline index then we have to delete the previous coefficient from the coefficient array and put the new index there
         if (data%i .ne. data%i_ref) then
            data%ep_bcoef(data%i_ref) = 0.0_ep1
            data%ep_bcoef(data%i) = 1.0_ep1
            data%i_ref = data%i
         endif

         r = BVALU(data%ep_knots,data%ep_bcoef,data%n,data%k,data%id,x,data%inbv,data%ep_work)
         r = r*x**data%l

 end function ep_bspl_prod_pow_evaluation

 real(wp) function wp_bessel_eval(data,x) result(r)
      use special_functions, only: cfp_besj
      use utils, only: xermsg
      implicit none
      class(bessel_fn) :: data
      real(wp), intent(in) :: x

      real(wp) :: z, y(1), alpha, pow
      real(wp), parameter :: half = 0.5_wp
      integer, parameter :: n =1
      integer :: nz, i

         if (data%k <= 0 .or. data % l < 0 .or. data % p < 0) then
            call xermsg ('function_integration', 'bessel_eval', 'Invalid input parameters: k and/or l and/or p.', 1, 1)
         end if

         z = data%k*x
         alpha = half + data%l
         call cfp_besj(z, alpha, n, y, nz)

         if (nz .ne. 0) call xermsg('function_integration','bessel_eval','Underflow in the Bessel function computation.',1,0)

         !x^p
         pow = x
         do i=2,data%p
            pow = pow*x
         enddo

         r = y(1)*pow

 end function wp_bessel_eval

 real(ep1) function ep_bessel_eval(data,x) result(r)
      use special_functions, only: cfp_besj
      use utils, only: xermsg
      implicit none
      class(bessel_fn) :: data
      real(ep1), intent(in) :: x

         r = 0.0_ep1
         
         call xermsg('function_integration','ep_bessel_eval','Quad precision not implemented yet.',1,1)

 end function ep_bessel_eval

end module function_integration
