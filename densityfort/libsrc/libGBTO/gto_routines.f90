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
!> \ingroup GTO
!> This module contains elementary routines which are used for evaluation of GTO-only integrals.
module gto_routines
use precisn
use special_functions
use coupling_obj, only: couplings_type
use utils, only: xermsg

  !> This object contains routines for evaluation of the Boys function in double precision.
  !> The object must be first initialized using the method init and its argument values imax and mmax. imax corresponds to the largest power in series expansion of the Boys function that is needed to 
  !> evaluate the Boys function in the expected range of T values. For T in the range 0,60 the value of imax has been determined to be 140.
  !> The routine eval then calculates at once the string of Boys functions for a range of m values.
  type boys_function_obj
     !> Maximum i for which the fac_terms have been calculated.
     integer, private :: imax = 0
     !> Maximum m for which the fac_terms have been calculated.
     integer, private :: mmax = 0
     !> The values \f$(2*m-1)!!/(2*m+2*i+1)!\f$ calculated for \f$i=0,\dots,imax\f$ and \f$m=0,\dots,mmax\f$.
     real(kind=cfp), allocatable, private :: fac_terms(:,:)
     !> Order of the Taylor series for the Boys function calculation.
     integer, private :: k = 0
     !> The number of points for which the Boys function has been precalculated.
     integer, private :: grid_len = 0
     !> Step for which the T values have been calculated on the grid.
     real(kind=cfp), private :: step = 0.0_cfp
     !> The Boys function \f$F_{m}(T)\f$ evaluated on the regularly spaced grid of grid_len values in the range T = 0,..,Tmax. For each T the m-values are: \f$m=0,\dots,mmax+1+k\f$.
     real(kind=cfp), private :: Tmax = 0
     !> The Boys function \f$F_{m}(T)\f$ evaluated on the regularly spaced grid of grid_len values in the range T = 0,..,boys_f_dprec_asym_thr. For each T the m-values are: \f$m=0,\dots,mmax+1+k\f$.
     !> This array holds the precalculated values of \f$F_{m}(T)\f$.
     real(kind=cfp), allocatable, private :: grid_val(:)
     !> Auxiliary array used by eval_taylor.
     real(kind=cfp), allocatable, private :: a(:), f(:)
     !> Set to .true. after initialization.
     logical, private :: initialized = .false.
  contains 
     !> Calculates the values fac_terms for the i and m values whose range is determined by the arguments imax and mmax. Also, the Boys function is precalculated for the number of values given by grid_len.
     procedure :: init => init_boys
     !> Returns the string of Boys functions \f$F_{m}(T)\f$ for \f$m=0,\dots,m_{max}\f$.
     !> \f[
     !>    F_{m}(T) = \int_{0}^{1}x^{2m}\exp[-Tx^2]dx.
     !> \f]
     !> The result is a double precision array boys_function(1:mmax+1) containing the value of the Boys function for \f$m=0,\dots,m_{max}\f$.
     !> The Boys function is calculated first for \f$m = m_{max}\f$ using the series expansion:
     !> \f[
     !>    F_{m}(T) = \exp[-T]\sum_{i=0}^{\infty}\frac{(2m-1)!!(2T)^{i}}{(2m+2i+1)!!}
     !> \f]
     !> and then for \f$m = m_{max}-1,\dots,0\f$ using the recurrence:
     !> \f[
     !>    F_{m}(T) = \frac{2TF_{m+1}(T)+\exp[-T]}{2m+1}.
     !> \f]
     !> For T > 60 we use the asymptotic formula for the Boys function:
     !> \f[
     !>    F_{m}(T) \approx \frac{1}{2T^{m+1/2}}\Gamma(m+1/2).
     !> \f]
     !> \param[in] T Real value corresponding to T in \f$F_{m}(T)\f$.
     !> \param[in] mmax Integer value corresponding to \f$m_{max}\f$ in \f$F_{m}(T)\f$, \f$m=0,\dots,m_{max}\f$.
     !> boys_function Double precision array boys_function(1:mmax+1) containing the values: \f$F_{m}(T)\f$, \f$m=0,\dots,m_{max}\f$.
     !> \warning It is vital for calculation of the two electron integrals that the Boys function is calculated as accurately as possible since it is used to start the recurrent evaluation of the integrals.
     !> It is important to understand that the limits on the maximum GTO L used are connected with the precision of calculation of the corresponding Boys functions. Given the largest GTO L in the basis the 
     !> largest required mmax is 4L. The parameters involved in the calculation of the Boys function (see the module const) have been tuned for the case L_max=6, i.e. mmax=24.
     procedure :: eval => eval_boys_function
     !> Uses the combination of the Taylor expansion and the Horner scheme for polynomial evaluation and the asymptotic expansion to obtain the value of the Boys function.
     !> This routine becomes available following the initialization which determines the order of the Taylor expansion and the density of the tabulated grid of values of the Boys function.
     procedure :: eval_taylor => eval_boys_function_taylor
     !>
     procedure, private :: eval_grid => eval_taylor_grid
     !> Deallocates the array fac_terms and resets the values of imax and mmax.
     procedure :: final => final_boys
  end type boys_function_obj

  private init_boys, eval_boys_function, eval_boys_function_taylor, final_boys, eval_taylor_grid

  !> We use the same object accross this module to calculate all angular couplings. Each routine first calls precalculate which precalculates the Gaunt coefficients (if needed) and stores them efficiently.
  !> This is done to ensure the couplings are always retrieved with maximum efficiency.
  type(couplings_type), private :: cpl

  !> This array is used by the eri_tail_shell routine to perform the transposition of the calculated integrals batch if neccessary.
  real(kind=cfp), allocatable, private :: eri_tail_tmp(:)

  !> This array is used by reorder_and_index_2el to reorder the calculated 2-electron integrals when the MPI parallelization is used.
  real(kind=cfp), allocatable, private :: ints_tmp(:)

  !$OMP THREADPRIVATE(cpl,eri_tail_tmp,ints_tmp)

contains

  function init_boys(this,imax,mmax,step,k)
     use const, only: boys_f_dprec_asym_thr, stdout
     implicit none
     class(boys_function_obj) :: this
     integer, intent(in) :: imax,mmax,k
     real(kind=cfp), intent(in) :: step
     integer :: init_boys
 
     integer :: i, m, prev, no_m
     real(kind=cfp) :: fac

        init_boys = 0
        if (this%initialized) then
           init_boys = 1
           return
        endif

        if (imax < 0 .or. mmax < 0 .or. step .le. 0.0_cfp .or. k .le. 0) then
           init_boys = 2
           return
        endif

        this%imax = imax
        this%mmax = mmax + k !we need to make mmax larger by k due to the Taylor expansion method which requires, for a given mmax, m values up to mmax+k
        this%k = k

        if (imax .eq. 0) this%imax = 1
        if (mmax .eq. 0) this%mmax = 1

        no_m = this%mmax+1 !the number of m values for which each Boys function is precalculated

        allocate(this%fac_terms(0:this%imax,0:this%mmax),this%a(1:this%k+1),this%f(1:this%k),stat=init_boys)
        if (init_boys .ne. 0) then
           init_boys = 3
           return
        endif 

        !precalculate the factorial terms
        do m=0,this%mmax
           fac = 2*m + 1.0_cfp
           prev = 2*m + 1
           this%fac_terms(0,m) = 1.0_cfp/fac

           do i=1,this%imax
              prev = prev + 2
              fac = fac*prev 
              this%fac_terms(i,m) = 1.0_cfp/fac !1/fac = (2*m-1)!!/(2*m+2*i+1)!
           enddo

        enddo

        !precalculate some coefficients for the Taylor expansion
        fac = 1.0_cfp
        do i=1,this%k
           fac = -fac*i
           this%f(i) = 1.0_cfp/fac
        enddo

        this%initialized = .true.

        !precalculate the Boys function for the Tmax=boys_f_dprec_asym_thr and the T-step=0.1.
        call this%eval_grid(boys_f_dprec_asym_thr,step)

  end function init_boys

  subroutine eval_taylor_grid(this,Tmax,step)
     use utils, only: xermsg
     use const, only: stdout
     use omp_lib
     implicit none
     class(boys_function_obj) :: this
     real(kind=cfp), intent(in) :: Tmax, step

     integer :: err, i, j, no_m, grid_len
     real(kind=cfp) :: T

        if (allocated(this%grid_val)) deallocate(this%grid_val)

        no_m = this%mmax+1 !the number of m values for which each Boys function is precalculated; we assume mmax has been set to the value appropriately larger, see init_boys

        grid_len = ceiling(Tmax/step+1)
        this%grid_len = grid_len 
        this%step = step

        allocate(this%grid_val(1:grid_len*no_m),stat=err)
        if (err .ne. 0) call xermsg('gto_routines','eval_taylor_grid','Memory allocation failed.',err,1)

        !precalculate the Boys function on the given grid
        if (.not.(omp_in_parallel())) then
           write(stdout,'("Precalculating the Boys function on the grid of ",i6," points and for mmax = ",i3,&
                        &" and k = ",i3,"; Tmax= ",e20.10)') this%grid_len, this%mmax-this%k, this%k, this%step*(grid_len-1)
        endif

        this%Tmax = Tmax
        j = 1
        do i=1,grid_len
           T = this%step*(i-1) !the grid of T values starts with T=0
           call this%eval(this%grid_val(j:j+this%mmax),no_m,T,this%mmax)
           j = j + no_m
        enddo

        if (.not.(omp_in_parallel())) then
           write(stdout,'("..done")')
        endif

  end subroutine eval_taylor_grid

  subroutine eval_boys_function_taylor(this,boys,d,T,mmax)
     use precisn
     use utils, only: xermsg
     use const, only: fit_order, wp_fit_terms, ep_fit_terms, wp_fit_terms_up, ep_fit_terms_up
     use phys_const, only: rtpi_half
     implicit none
     class(boys_function_obj) :: this
     real(kind=cfp), intent(in) :: T
     integer, intent(in) :: mmax, d
     real(kind=cfp), intent(inout) :: boys(1:d)

     !These values are used in the Taylor expansion of Fm(z) and are given by the expression T(k) = (-1)**k / (k!)
     real(kind=cfp), parameter :: T2 =  0.500000000000000000000000000000000000000_cfp
     real(kind=cfp), parameter :: T3 = -0.166666666666666666666666666666666666666_cfp
     real(kind=cfp), parameter :: T4 =  0.041666666666666666666666666666666666666_cfp
     real(kind=cfp), parameter :: T5 = -0.008333333333333333333333333333333333333_cfp
     real(kind=cfp), parameter :: T6 =  0.001388888888888888888888888888888888888_cfp
     real(kind=cfp), parameter :: T7 = -0.000198412698412698412698412698412698412_cfp
     real(kind=cfp), parameter :: T8 =  0.000024801587301587301587301587301587301_cfp
     real(kind=cfp), parameter :: T9 = -2.755731922398589065255731922398589E-06_cfp
     real(kind=cfp), parameter :: T10 = 2.755731922398589065255731922398589E-07_cfp

     !These values are used in the downward recursion for Fm(z) and are given by the expression M(k) = 1/(2*k+1)
     real(kind=cfp), parameter :: M1 = 0.3333333333333333333333333333333333333333_cfp
     real(kind=cfp), parameter :: M2 = 0.2000000000000000000000000000000000000000_cfp
     real(kind=cfp), parameter :: M3 = 0.1428571428571428571428571428571428571429_cfp
     real(kind=cfp), parameter :: M4 = 0.1111111111111111111111111111111111111111_cfp
     real(kind=cfp), parameter :: M5 = 0.0909090909090909090909090909090909090909_cfp
     real(kind=cfp), parameter :: M6 = 0.0769230769230769230769230769230769230769_cfp
     real(kind=cfp), parameter :: M7 = 0.0666666666666666666666666666666666666666_cfp

     real(kind=cfp) :: T0, delta, two_T, exp_T
     integer :: i, k, m, f, test_asym, test_upward

        if (.not.(this%initialized)) then
           call xermsg('gto_routines','eval_boys_function_taylor','The Boys function object has not been initialized.',1,1)
        endif

        if (T < 0.0_cfp .or. mmax < 0 .or. d < mmax + 1) then
           call xermsg('gto_routines','eval_boys_function_taylor','Invalid input parameters.',2,1)
        endif

        if (T .eq. 0.0_cfp) then
           !todo use M1,...,M7 to do this
           forall (m=0:mmax)
              boys(1+m) = 1.0_cfp/(2*m+1.0_cfp)
           endforall
           return
        endif

        !The estimate for use of the upward recursion differs depending on whether we require full double or quad precision result. Similarly for the asymptotic method below.
        if (cfp .eq. wp) then
           test_upward = nint(cfp_eval_poly_horner(fit_order,T,wp_fit_terms_up))
        elseif (cfp .eq. ep1) then
           test_upward = nint(cfp_eval_poly_horner(fit_order,T,ep_fit_terms_up))
        endif

        if (test_upward .ge. mmax) then !use the upward recurrence whenever possible

           boys(1) = rtpi_half/(sqrt(T)) !m=0

           two_T = 1.0_cfp/(2.0_cfp*T)
           exp_T = exp(-T)

           !What follows is an unrolled version of the loop below:
           !do m=1,mmax
           !   boys(m+1) = boys(m)*two_T*(2*m-1) - exp_T
           !enddo

           if (mmax .eq. 1) then
              boys(2) = boys(1)*two_T - exp_T
           elseif (mmax .eq. 2) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
           elseif (mmax .eq. 3) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
           elseif (mmax .eq. 4) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
              boys(5) = boys(4)*two_T*7 - exp_T
           elseif (mmax .eq. 5) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
              boys(5) = boys(4)*two_T*7 - exp_T
              boys(6) = boys(5)*two_T*9 - exp_T
           elseif (mmax .eq. 6) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
              boys(5) = boys(4)*two_T*7 - exp_T
              boys(6) = boys(5)*two_T*9 - exp_T
              boys(7) = boys(6)*two_T*11 - exp_T
           elseif (mmax .eq. 7) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
              boys(5) = boys(4)*two_T*7 - exp_T 
              boys(6) = boys(5)*two_T*9 - exp_T
              boys(7) = boys(6)*two_T*11 - exp_T
              boys(8) = boys(7)*two_T*13 - exp_T
           elseif (mmax .eq. 8) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
              boys(5) = boys(4)*two_T*7 - exp_T
              boys(6) = boys(5)*two_T*9 - exp_T
              boys(7) = boys(6)*two_T*11 - exp_T
              boys(8) = boys(7)*two_T*13 - exp_T
              boys(9) = boys(8)*two_T*15 - exp_T
           elseif (mmax > 8) then
              boys(2) = boys(1)*two_T - exp_T
              boys(3) = boys(2)*two_T*3 - exp_T
              boys(4) = boys(3)*two_T*5 - exp_T
              boys(5) = boys(4)*two_T*7 - exp_T
              boys(6) = boys(5)*two_T*9 - exp_T
              boys(7) = boys(6)*two_T*11 - exp_T
              boys(8) = boys(7)*two_T*13 - exp_T
              boys(9) = boys(8)*two_T*15 - exp_T

              !todo perform this in chunks of 8
              do m=9,mmax
                 boys(m+1) = boys(m)*two_T*(2*m-1) - exp_T
              enddo
           endif

           return
        endif

        !The value of test_mmax for which Fm(T), m=0,...,test_mmax are calculated to full double/quad precision using the asymptotic formula. If the given mmax value is smaller than test_mmax
        !then we can safely use the asymptotic formula.
        if (cfp .eq. wp) then
           test_asym = nint(cfp_eval_poly_horner(fit_order,T,wp_fit_terms))
        elseif (cfp .eq. ep1) then
           test_asym = nint(cfp_eval_poly_horner(fit_order,T,ep_fit_terms))
        endif

        if (test_asym .ge. mmax) then !get F_{mmax}(T) using the asymptotic formula

           boys(mmax+1) = cfp_gamma_fun(mmax+0.5_cfp)*0.5_cfp/(T**(mmax+0.5_cfp))

        else !use the Taylor expansion around the grid point T0 closest to T to get F_{mmax}(T)

           !recalculate the grid if the requested T value lies outside of the range for which the grid has been originally evaluated.
           if (T > this%Tmax) then
              call this%eval_grid(T,this%step)
           endif
           if (mmax > this % mmax) then
              call xermsg ('gto_routines', 'eval_boys_function_taylor', &
                           'The requested mmax is larger than the mmax for which the boys function object &
                           &has been initialized.', 3, 1)
           end if

           i = nint(T/this%step)+1 !determine the index of the T value closest to this one
           T0 = this%step*(i-1) !the T value used for the starting point for the Taylor expansion

           delta = T-T0
           m = (i-1)*(this%mmax+1)+1 !the index of F_{0}(T0)

           !Taylor expansion coefficients for F_{mmax}(T) expanded around T0.
           m = m + mmax

           !The lines below are unrolled versions of the following:
           !f = 1
           !do k=1,this%k
           !   f = -f*k
           !   this%a(k+1) = this%grid_val(m+k)/real(f,kind=cfp)
           !enddo

           if (this%k .eq. 1) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)

           elseif (this%k .eq. 2) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)

           elseif (this%k .eq. 3) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)

           elseif (this%k .eq. 4) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)

           elseif (this%k .eq. 5) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)

           elseif (this%k .eq. 6) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)
              this%a(7) = T6*this%grid_val(m+6)
           
           elseif (this%k .eq. 7) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)
              this%a(7) = T6*this%grid_val(m+6)
              this%a(8) = T7*this%grid_val(m+7)

           elseif (this%k .eq. 8) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)
              this%a(7) = T6*this%grid_val(m+6)
              this%a(8) = T7*this%grid_val(m+7)
              this%a(9) = T8*this%grid_val(m+8)

           elseif (this%k .eq. 9) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)  
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)
              this%a(7) = T6*this%grid_val(m+6)
              this%a(8) = T7*this%grid_val(m+7)
              this%a(9) = T8*this%grid_val(m+8)
              this%a(10)= T9*this%grid_val(m+9)

           elseif (this%k .eq. 10) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)
              this%a(7) = T6*this%grid_val(m+6)
              this%a(8) = T7*this%grid_val(m+7)
              this%a(9) = T8*this%grid_val(m+8)
              this%a(10)= T9*this%grid_val(m+9)
              this%a(11)=T10*this%grid_val(m+10)

           elseif (this%k > 10) then

              this%a(1) = this%grid_val(m)
              this%a(2) =-this%grid_val(m+1)
              this%a(3) = T2*this%grid_val(m+2)
              this%a(4) = T3*this%grid_val(m+3)
              this%a(5) = T4*this%grid_val(m+4)
              this%a(6) = T5*this%grid_val(m+5)
              this%a(7) = T6*this%grid_val(m+6)
              this%a(8) = T7*this%grid_val(m+7)
              this%a(9) = T8*this%grid_val(m+8)
              this%a(10)= T9*this%grid_val(m+9)
              this%a(11)=T10*this%grid_val(m+10)

              !todo perform this in chunks of 8
              f = 3628800 !The value of f is obtained as in the commented loop over 'm' above.
              do k=11,this%k
                 f = -f*k
                 this%a(k+1) = this%grid_val(m+k)/real(f,kind=cfp)
              enddo

           endif

           !get F_{mmax}(T) using the Taylor expansion
           boys(mmax+1) = cfp_eval_poly_horner(this%k,delta,this%a)

        endif
   
        !use the downward recursion to calculate Boys function for m=mmax-1,...,0
        two_T = 2.0_cfp*T
        exp_T = exp(-T)

        !The lines below are unrolled versions of the following:
        !do m=mmax-1,0,-1
        !   boys(m+1) = (boys(m+2)*two_T+exp_T)/real(2*m+1,kind=cfp)
        !enddo

        if (mmax .eq. 1) then
           boys(1) = (boys(2)*two_T+exp_T)
        elseif (mmax .eq. 2) then
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax .eq. 3) then
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax .eq. 4) then
           boys(4) = (boys(5)*two_T+exp_T)*M3 !m=3
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax .eq. 5) then
           boys(5) = (boys(6)*two_T+exp_T)*M4 !m=4
           boys(4) = (boys(5)*two_T+exp_T)*M3 !m=3
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax .eq. 6) then
           boys(6) = (boys(7)*two_T+exp_T)*M5 !m=5
           boys(5) = (boys(6)*two_T+exp_T)*M4 !m=4
           boys(4) = (boys(5)*two_T+exp_T)*M3 !m=3
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax .eq. 7) then
           boys(7) = (boys(8)*two_T+exp_T)*M6 !m=6
           boys(6) = (boys(7)*two_T+exp_T)*M5 !m=5
           boys(5) = (boys(6)*two_T+exp_T)*M4 !m=4
           boys(4) = (boys(5)*two_T+exp_T)*M3 !m=3
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax .eq. 8) then
           boys(8) = (boys(9)*two_T+exp_T)*M7 !m=7
           boys(7) = (boys(8)*two_T+exp_T)*M6 !m=6
           boys(6) = (boys(7)*two_T+exp_T)*M5 !m=5
           boys(5) = (boys(6)*two_T+exp_T)*M4 !m=4
           boys(4) = (boys(5)*two_T+exp_T)*M3 !m=3
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        elseif (mmax > 8) then
           !todo perform this in chunks of 8 explicitly unrolled
           do m=mmax-1,8,-1
              boys(m+1) = (boys(m+2)*two_T+exp_T)/real(2*m+1,kind=cfp)
           enddo
           boys(8) = (boys(9)*two_T+exp_T)*M7 !m=7
           boys(7) = (boys(8)*two_T+exp_T)*M6 !m=6
           boys(6) = (boys(7)*two_T+exp_T)*M5 !m=5
           boys(5) = (boys(6)*two_T+exp_T)*M4 !m=4
           boys(4) = (boys(5)*two_T+exp_T)*M3 !m=3
           boys(3) = (boys(4)*two_T+exp_T)*M2 !m=2
           boys(2) = (boys(3)*two_T+exp_T)*M1 !m=1
           boys(1) = (boys(2)*two_T+exp_T)    !m=0
        endif

  end subroutine eval_boys_function_taylor

  subroutine eval_boys_function(this,boys,d,T,mmax)
     use precisn
     use utils, only: xermsg
     use const, only: boys_tol, wp_fit_terms, ep_fit_terms, wp_fit_terms_up, ep_fit_terms_up
     use phys_const, only: pi
     use special_functions, only: boys_function_quad
     implicit none
     class(boys_function_obj) :: this
     real(kind=cfp), intent(in) :: T
     integer, intent(in) :: mmax, d
     real(kind=cfp), intent(inout) :: boys(1:d)

     integer :: m, i, test, test_upward
     real(kind=cfp) :: T_pow, two_T, exp_T, s, s_prev, term

        if (.not.(this%initialized)) then
           call xermsg('gto_routines','eval_boys_function','The Boys function object has not been initialized.',1,1)
        endif

        if (T < 0.0_cfp .or. mmax < 0 .or. d < mmax + 1) then
           call xermsg('gto_routines','eval_boys_function','Invalid input parameters.',2,1)
        endif

        if (mmax > this%mmax) then !in this case the non-asymptotic algorithm will proceed
           print *,mmax,this%mmax
           call xermsg ('gto_routines', 'eval_boys_function', &
                        'The factorial term is required for mmax larger than the mmax for which the boys function object &
                        &has been initialized.', 3, 1)
        endif

        two_T = 2.0_cfp*T
        exp_T = exp(-T)

        !The estimate for use of the upward recursion differs depending on whether we require full double or quad precision result. Similarly for the asymptotic method below.
        if (cfp .eq. wp) then
           test_upward = nint(cfp_eval_poly_horner(2,T,wp_fit_terms_up))
        elseif (cfp .eq. ep1) then
           test_upward = nint(cfp_eval_poly_horner(2,T,ep_fit_terms_up))
        endif

        if (test_upward .ge. mmax) then !use the upward recurrence whenever possible

           boys(1) = sqrt(pi)*0.5_cfp/(sqrt(T)) !m=0

           two_T = 1.0_cfp/two_T

           do m=1,mmax
              boys(m+1) = boys(m)*two_T*(2*m-1)-exp_T
           enddo

           return
        endif

        !The value of test_mmax for which Fm(T), m=0,...,test_mmax are calculated to full double/quad precision using the asymptotic formula. If the given mmax value is smaller than test_mmax
        !then we can safely use the asymptotic formula.
        if (cfp .eq. wp) then
           test = nint(cfp_eval_poly_horner(2,T,wp_fit_terms))
        elseif (cfp .eq. ep1) then
           test = nint(cfp_eval_poly_horner(2,T,ep_fit_terms))
        endif

        if (test .ge. mmax) then !use the asymptotic formula; 85 is the smallest exponent for which the asymptotic formula gives fully accurate results up to mmax = 24, but we can afford 60.

           boys(mmax+1) = cfp_gamma_fun(mmax+0.5_cfp)*0.5_cfp/(T**(mmax+0.5_cfp))
           
        else !use the power series for F_{m}(T)

           if (cfp .eq. wp .and. T > 65.0_cfp) then
              !Branch into quad precision to avoid numerical problems (see below). This should happen only in extreme cases.
              !In any case this routine should be only used to evaluate the grid for the Taylor expansion method so the fact we are calculating in quad precision here should not lead to any significant
              !overheads in the actual calculation.
              boys(1:mmax+1) = real(boys_function_quad(real(T,kind=ep1),mmax),kind=cfp)
              return
           else !In all other cases we should be OK but keep on mind that we can potentially run into the same numerical problems even in quad precision and this case is not covered here.

              !i=0: starting term for the sum
              T_pow = 1.0_cfp
              s = T_pow*this%fac_terms(0,mmax)
              s_prev = 0.0_cfp
      
              !sum over i=1,...,until convergence
              !note that the terms responsible for numerical problems for T > 65 are T_pow and fac; if these are evaluated in quad precision then the limitations on T and mmax can be dropped.
              i = 0
              do
                 i = i + 1
                 if (i > this%imax) then
                    print *,i,this%imax
                    call xermsg ('gto_routines', 'eval_boys_function', &
                                 'The factorial term is required for i larger than the imax for which the boys function object &
                                 &has been initialized.', 4, 1)
                 endif
                 T_pow = T_pow*two_T   !(2*T)**i
                 term = T_pow*this%fac_terms(i,mmax)      !the next term in the sum: (2*m-1)!!/(2*m+2*i+1)!!*(2*T)**i
                 s = s + term
                 if (term .le. s*boys_tol) exit !convergence criterion
                 s_prev = s
              enddo
      
              boys(mmax+1) = exp_T*s !Boys function for m=mmax
           endif

        endif

        !Boys function for m=mmax-1,...,0
        do m=mmax-1,0,-1
           boys(m+1) = (boys(m+2)*two_T+exp_T)/real(2*m+1,kind=cfp)
        enddo

  end subroutine eval_boys_function

  function final_boys(this)
     implicit none
     class(boys_function_obj) :: this
     integer :: final_boys

        final_boys = 0

        if (this%initialized) then
           deallocate(this%fac_terms,this%grid_val,this%a,this%f,stat=final_boys)
           if (final_boys .ne. 0) then
              final_boys = 1
              return
           endif

           this%imax = 0
           this%mmax = 0
           this%grid_len = 0
           this%k = 0
           this%step = 0.0_cfp
           this%Tmax = 0.0_cfp

           this%initialized = .false.
        endif

  end function final_boys

  !> \par Purpose:
  !> \verbatim
  !> Calculate the normalization factor for a GTO orbital.
  !> \endverbatim
  !> \f[
  !>   N^{GTO}_{\alpha l} = \sqrt{\frac{2(2\alpha)^{l+3/2}}{\Gamma(l+3/2)}}\sqrt{\frac{2l+1}{4\pi}}.
  !> \f]
  !> \verbatim
  !>  Gamma is the gamma function.
  !> \endverbatim
  !
  !> \param[in] l, alpha
  !> \verbatim
  !>  Integer angular momentum l of the GTO and its exponent alpha (real).
  !> \endverbatim
  pure function dngto(l,alpha)
    use precisn
    use phys_const, only: fourpi
    implicit none
    integer, intent(in) :: l
    real(kind=cfp), intent(in) :: alpha
    real(kind=cfp) :: dngto
  
    dngto = sqrt(2.0_cfp*((2*alpha)**(l+1.5_cfp))/cfp_gamma_fun(l+1.5_cfp))*sqrt((2*l+1.0_cfp)/fourpi)  
  
  end function dngto

  !> This function calculates \f$ I = \int\int\int x^{i}y^{j}z^{k}\exp[-(\alpha+\beta)r^2] = {\frac{\Gamma(i/2+1/2)\Gamma(j/2+1/2)\Gamma(k/2+1/2)}{((\alpha+\beta)^{1/2(i+j+k+3)})}} \f$. 
  !> In particular, this function can be used to calculate overlap of two primitive GTOs centered on the same nucleus.
  function cart_gto_int(alp,bet,i,j,k)
     use utils, only: xermsg
     implicit none
     real(kind=cfp) :: cart_gto_int 
     integer, intent(in) :: i, j, k
     real(kind=cfp), intent(in) :: alp, bet

     integer :: it 
     real(kind=cfp) :: sum_exp

        if (alp <= 0.0_cfp .or. bet <= 0 .or. i < 0 .or. j < 0 .or. k < 0) then
            call xermsg ('gto_routines', 'cart_gto_int', &
                         'One or more of the following values are invalid: alp, bet, i, j, k.', 1, 1)
        end if

        if (mod(i,2) .ne. 0 .or. mod(j,2) .ne. 0 .or. mod(k,2) .ne. 0) then !if one of the polynomial exponents is odd then the integral vanishes since exp(-a*r^2) is symmetrical
           cart_gto_int = 0.0_cfp
           return
        endif

        sum_exp = alp+bet
        cart_gto_int = sum_exp

        do it=2,i+j+k+3 !this loop gives (alp+bet)^(i+j+k+3); it is faster than (bet*alp)**(i+j+k+3)
           cart_gto_int = cart_gto_int*sum_exp
        enddo
        cart_gto_int = sqrt(cart_gto_int) !cart_gto_int = (alp+bet)**(0.5*(i+j+k+3))

        cart_gto_int = cfp_gamma_fun(0.5_cfp*(i+1))*cfp_gamma_fun(0.5_cfp*(j+1))*cfp_gamma_fun(0.5_cfp*(k+1))/cart_gto_int

        !norm of a primitive GTO = sqrt(sum_exp**(i+j+k+1.5_cfp)/(cfp_gamma_fun(i+0.5_cfp)*cfp_gamma_fun(j+0.5_cfp)*cfp_gamma_fun(k+0.5_cfp)))

  end function cart_gto_int

  !> Calculates the normalization factor for the primitive cartesian GTO function.
  function norm_cart_gto(alp,i,j,k)
     use utils, only: xermsg
     implicit none
     real(kind=cfp) :: norm_cart_gto
     integer, intent(in) :: i, j, k
     real(kind=cfp), intent(in) :: alp

        if (alp <= 0.0_cfp .or. i < 0 .or. j < 0 .or. k < 0) then
            call xermsg ('gto_routines', 'norm_cart_gto', 'One or more of the following values are invalid: alp, i, j, k.', 1, 1)
        end if

        norm_cart_gto = 1.0_cfp/sqrt(cart_gto_int(alp,alp,2*i,2*j,2*k))

  end function norm_cart_gto

  !> This function calculates the normalization factor for a contracted cartesian GTO. The cartesian exponents are i,j,k. The contraction coefficients and exponents of the primitives are in the arrays
  !> alp(:) and ccf(:).
  function contr_cart_gto_norm(i,j,k,alp,ccf)
     use utils, only: xermsg
     implicit none
     real(kind=cfp) :: contr_cart_gto_norm
     integer, intent(in) :: i, j, k
     real(kind=cfp), intent(in) :: alp(:), ccf(:)

     integer :: it, jt, n, two_i, two_j, two_k

        n = size(alp)
        if (n /= size(ccf)) then
            call xermsg ('gto_routines', 'contr_cart_gto_norm', &
                         'The number of exponents does not match the number of contractions.', 1, 1)
        end if

        if (i < 0 .or. j < 0 .or. k < 0) then
            call xermsg ('gto_routines', 'contr_cart_gto_norm', 'One or more of the following values are invalid: i, j, k.', 2, 1)
        end if

        !cartesian exponents that enter the integrals over the primtive GTOs.
        two_i = 2*i
        two_j = 2*j
        two_k = 2*k

        contr_cart_gto_norm = 0.0_cfp !first we calculate the self-overlap of the contracted cartesian GTO
        do it=1,n !over all primitives
           do jt=1,n !over all primitives
              contr_cart_gto_norm = contr_cart_gto_norm + ccf(it)*ccf(jt)*norm_cart_gto(alp(jt),i,j,k) &
                                *norm_cart_gto(alp(it),i,j,k)*cart_gto_int(alp(it),alp(jt),two_i,two_j,two_k)
           enddo
       enddo

       contr_cart_gto_norm = 1.0_cfp/sqrt(contr_cart_gto_norm) !norm of the contracted cartesian GTO

  end function contr_cart_gto_norm

  !> This function calculates overlap of a contracted cartesian GTO with a contracted spherical GTO (both sitting on the same nucleus). 
  !> The exponents of the contracted cartesian GTO are: i,j,k. The L,M values of the spherical GTO are: l,m. The exponents and contraction coefficients of the primitives are in the arrays alp(:), ccf(:).
  function olap_ccart_csph(i,j,k,l,m,alp,ccf)
     use utils, only: xermsg
     implicit none
     real(kind=cfp) :: olap_ccart_csph
     integer, intent(in) :: i, j, k, l, m
     real(kind=cfp), intent(in) :: alp(:), ccf(:)

     integer :: n, it, jt, kt, no_terms
     integer, allocatable :: i_car(:), j_car(:), k_car(:)
     real(kind=cfp), allocatable :: c(:)
     real(kind=cfp) :: prod, n_c, n_s

        n = size(alp)
        if (n /= size(ccf)) then
            call xermsg ('gto_routines', 'olap_ccart_csph', &
                         'The number of exponents does not match the number of contractions.', 1, 1)
        end if

        call cfp_sph_to_cart_mapping(l,m,c,i_car,j_car,k_car) !obtain the exponents and coefficients of the cartesian GTOs which build the spherical GTO.
        no_terms = size(c) !number of terms in the cartesian -> spherical mapping

        olap_ccart_csph = 0.0_cfp
        do it=1,n !over all cartesian primitives
           n_c = norm_cart_gto(alp(it),i,j,k) !norm of the cartesian primitive
           do jt=1,n !over all spherical primitives
              n_s = dngto(l,alp(jt)) !norm of the spherical primitive
              prod = ccf(it)*ccf(jt)*n_c*n_s
              do kt=1,no_terms !over all cartesian terms in the spherical primitive
                 olap_ccart_csph = olap_ccart_csph + prod*c(kt)*cart_gto_int(alp(it),alp(jt),i_car(kt)+i,j_car(kt)+j,k_car(kt)+k)
                 !cart_gto_int = overlap between a primitive cartesian contributing to the jt-th primitive spherical GTO and it-th primitive cartesian GTO contributing to the contracted cartesian GTO
              enddo
           enddo
        enddo

  end function olap_ccart_csph

  !> This routine performs the conversion of the orbital coefficients (for one shell) for the cartesian contracted GTOs to the basis of contracted spherical GTOs.
  !> The input are the data for the cartesian GTOs (l,ix,iy,iz,alp,ccf) and the corresponding orbital coefficients (cart_cf) and the overall norm of the spherical contracted GTO (sph_cgto_norm). It is
  !> assumed that the coefficients correspond to normalized contracted cartesian GTOs. The output are the orbital coefficients sph_cf for the corresponding shell of contracted spherical GTOs with the
  !> overall norm given by sph_cgto_norm.
  subroutine cart_cf_sph_cf(l,ix,iy,iz,alp,ccf,sph_cgto_norm,cart_cf,sph_cf)
     use utils, only: xermsg
     implicit none
     integer, intent(in) :: l, ix(:),iy(:),iz(:)
     real(kind=cfp), intent(in) :: alp(:), ccf(:), sph_cgto_norm, cart_cf(:)
     real(kind=cfp), intent(out) :: sph_cf(:)

     integer :: s, no_cart_gtos, no_sph_gtos, k, m
     real(kind=cfp) :: olap

        if (l < 0) call xermsg ('gto_routines', 'cart_cf_sph_cf', 'The l value on input is < 0.', 1, 1)

        no_cart_gtos = (l+1)*(l+2)/2
        no_sph_gtos = 2*l+1

        s = min(size(ix),size(iy),size(iz))
        if (s < no_cart_gtos) call xermsg ('gto_routines', 'cart_cf_sph_cf', 'The ix,iy or iz input data are incomplete.', 2, 1)

        s = size(sph_cf)
        if (s < no_sph_gtos) call xermsg ('gto_routines', 'cart_cf_sph_cf', 'The output array sph_cf is too small.', 3, 1)

        !calculate the contracted spherical GTO coefficients for all M using the (contracted cartesian GTO)x(contracted spherical GTO) overlaps
        do m=-l,l
           sph_cf(m+l+1) = 0.0_cfp
           do k=1,no_cart_gtos !loop over all contracted cartesian GTOs for this shell
              olap = olap_ccart_csph(ix(k),iy(k),iz(k),l,m,alp,ccf) !overlap between the contracted spherical and the contracted cartesian GTOs
              !multiply the overlap by the normalization factor of the contracted spherical GTO and the normalization of the contracted GTO function:
              olap = olap*contr_cart_gto_norm(ix(k),iy(k),iz(k),alp,ccf)*sph_cgto_norm
              !sph_cf(m+l+1) is the MO coefficient for the contracted spherical GTO with M=m:
              sph_cf(m+l+1) = sph_cf(m+l+1) + cart_cf(k)*olap
           enddo
        enddo !m

  end subroutine cart_cf_sph_cf

  !> This routine performs the conversion of the orbital coefficients (for one shell) for the spherical contracted GTOs to the basis of contracted cartesian GTOs.
  !> The input are the data for the cartesian GTOs (l,ix,iy,iz,alp,ccf) and the orbital coefficients for the spherical GTOs (sph_cf). 
  !> It is assumed that the coefficients correspond to contracted spherical GTOs normalized with sph_cgto_norm. The output are the orbital coefficients cart_cf for the corresponding shell of normalized 
  !> contracted cartesian GTOs.
  subroutine sph_cf_cart_cf(l,ix,iy,iz,alp,ccf,sph_cgto_norm,sph_cf,cart_cf)
     use utils, only: xermsg
     implicit none
     integer, intent(in) :: l, ix(:),iy(:),iz(:)
     real(kind=cfp), intent(in) :: alp(:), ccf(:), sph_cf(:), sph_cgto_norm
     real(kind=cfp), intent(out) :: cart_cf(:)

     integer :: s, no_cart_gtos, no_sph_gtos, k, m, no_terms, a, b, n
     real(kind=cfp) :: cart_norm, na, nb
     integer, allocatable :: i_car(:), j_car(:), k_car(:)
     real(kind=cfp), allocatable :: c(:)

        if (l < 0) call xermsg ('gto_routines', 'cart_cf_sph_cf', 'The l value on input is < 0.', 1, 1)

        no_cart_gtos = (l+1)*(l+2)/2
        no_sph_gtos = 2*l+1

        s = min(size(ix),size(iy),size(iz),size(cart_cf))
        if (s < no_cart_gtos) call xermsg ('gto_routines', 'cart_cf_sph_cf', 'The ix,iy,iz or cart_cf arrays are too small.', 2, 1)

        s = size(sph_cf)
        if (s < no_sph_gtos) call xermsg ('gto_routines', 'cart_cf_sph_cf', 'The input array sph_cf is too small.', 3, 1)

        n = size(alp)

        cart_cf = 0.0_cfp
        do m=-l,l
           call cfp_sph_to_cart_mapping(l,m,c,i_car,j_car,k_car) !obtain the cartesian exponents and coefficients of the cartesian GTOs which build the spherical GTO.
           no_terms = size(c) !number of terms in the cartesian -> spherical mapping
           do k=1,no_cart_gtos
              cart_norm = contr_cart_gto_norm(ix(k),iy(k),iz(k),alp,ccf) !norm of the contracted cartesian GTO
              do s=1,no_terms
                 !calculate contributions to the k-th cartesian coefficient from the matching cartesian GTO in the cartesian -> spherical mapping 
                 if (i_car(s) .eq. ix(k) .and. j_car(s) .eq. iy(k) .and. k_car(s) .eq. iz(k)) then
                    do a=1,n !primitives of the cartesian GTO
                       na = norm_cart_gto(alp(a),ix(k),iy(k),iz(k)) !norm of the cartesian primitive GTO
                       do b=1,n !primitives of the spherical GTO
                          nb = dngto(l,alp(b)) !norm of the spherical primitive GTO
                          cart_cf(k) = cart_cf(k) + ccf(a)*ccf(b)*c(s)*cart_norm*sph_cgto_norm*na*nb*sph_cf(m+l+1) &
                                                *cart_gto_int(alp(a),alp(b),2*ix(k),2*iy(k),2*iz(k))
                       enddo !b
                    enddo !a
                 endif
              enddo !s
           enddo !k
        enddo !m

  end subroutine sph_cf_cart_cf

  !> Computes the normalization factor for a CMS-centered contracted spherical GTO and radius given by rmat_radius. This is computed from the self-overlap of the GTO inside the sphere of radius rmat_radius.
  function cms_gto_norm(rmat_radius,l,np,alp,ccf,cnorm,norms)
     use special_functions, only: cfp_gamic, cfp_gamma_fun
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), intent(in) :: rmat_radius, alp(:), ccf(:), norms(:), cnorm
     integer, intent(in) :: l, np
     real(kind=cfp) :: cms_gto_norm

     integer :: i, j
     real(kind=cfp) :: prod, fac, r_sq, alp_ab, incomplete_gamma, arg, arg2, gam

        cms_gto_norm = 0.0_cfp

        fac = fourpi/(2*l+1.0_cfp)*cnorm*cnorm
        arg = l+1.5_cfp
        r_sq = rmat_radius*rmat_radius
        gam = cfp_gamma_fun(arg)

        do i=1,np
           do j=1,np
              prod = ccf(i)*norms(i)*ccf(j)*norms(j)*fac
              alp_ab = alp(i)+alp(j)
              arg2 = r_sq*alp_ab
              incomplete_gamma = cfp_gamic(arg,arg2)
              cms_gto_norm = cms_gto_norm + 0.5_cfp*prod*(alp_ab)**(-arg)*(gam-incomplete_gamma)
           enddo !j
        enddo !i

        cms_gto_norm = 1.0_cfp/sqrt(cms_gto_norm)

  end function cms_gto_norm

  !> Calculates the reduced boundary amplitude of a CGTO shell
  function CGTO_amplitude(r,l,number_of_primitives,norm,norms,contractions,exponents)
     use phys_const, only: fourpi
     implicit none
     integer, intent(in) :: l, number_of_primitives
     real(kind=cfp), intent(in) :: r, norm, norms(number_of_primitives), &
                                    contractions(number_of_primitives), exponents(number_of_primitives)
     real(kind=cfp) :: CGTO_amplitude

     integer :: i

         CGTO_amplitude = 0.0_cfp

         do i=1,number_of_primitives !loop over primitives
            CGTO_amplitude = CGTO_amplitude + norms(i)*contractions(i)*exp(-exponents(i)*r*r)
         enddo
         CGTO_amplitude = CGTO_amplitude*norm*sqrt(fourpi/(2*l+1.0_cfp))*r**(l+1)

  end function CGTO_amplitude

   !> Overlap and kinetic energy tail integrals for the normalized contracted spherical GTOs centered on the CMS. These integrals are given by:
   !> \f[
   !>    S_{tail} = \langle S_{LM}(\mathbf{r})\exp[-\alpha r^2] \vert S_{L'M'}(\mathbf{r})\exp[-\beta r^2] \rangle = \delta_{L'L}\delta_{M'M}\frac{4\pi}{2L+1}\frac{1}{2}(\alpha+\beta)^{-L-3/2}\Gamma[3/2+L,
   !>               a^2(\alpha+\beta)],
   !> \f]
   !> \f[
   !>    K_{tail} = \langle S_{LM}(\mathbf{r})\exp[-\alpha r^2] \vert\Delta\vert S_{L'M'}(\mathbf{r})\exp[-\beta r^2] \rangle = \delta_{L'L}\delta_{M'M}\frac{4\pi}{2L+1}\beta(\alpha+\beta)^{-\frac{5}{2}-L} 
   !>               \left(-(\alpha+\beta)(3+2 L)\Gamma\left[\frac{3}{2}+L,a^2 (\alpha+\beta)\right]+2 \beta\Gamma\left[\frac{5}{2}+L,a^2 (\alpha+\beta)\right]\right).
   !> \f]
   !> The delta terms are actually not taken into account in this routine - hence it is responsibility of the calling routine to subtract the calculated integrals only when appropriate.
   !> We assume that the mass of the particle is that of electron.
   !> The formula is for a pair of unnormalized primitives. The actual routine performs contraction over the primitives and multiplies by all the norms.
   !> The Bloch terms for the (ab) combination of contracted GTOs are also calculated. These should be added to the KE integrals following the tail subtraction.
   subroutine olap_kei_tail(l,np_a,np_b,alp_a,alp_b,ccf_a,ccf_b,cnorm_a, &
                            norms_a,cnorm_b,norms_b,rmat_radius,olap_tail,kei_tail,bloch_ab)
      use special_functions, only: cfp_gamic
      use phys_const, only: fourpi
      implicit none
      real(kind=cfp), intent(in) :: rmat_radius, alp_a(:), alp_b(:), ccf_a(:), ccf_b(:), norms_a(:), norms_b(:), cnorm_a, cnorm_b
      integer, intent(in) :: l, np_a, np_b
      real(kind=cfp), intent(out) :: olap_tail, kei_tail, bloch_ab

      real(kind=cfp) :: arg, arg2, arg3, fac, prod, incomplete_gamma, incomplete_gamma2, alp_ab, r_sq, fac_bloch
      integer :: i, j

         fac = fourpi/(2*l+1.0_cfp)
         fac_bloch = fac*rmat_radius**(2*l+1)
         arg = 1.5_cfp+l
         arg3 = 2.5_cfp+l
         r_sq = rmat_radius*rmat_radius
         olap_tail = 0.0_cfp
         kei_tail = 0.0_cfp
         bloch_ab = 0.0_cfp

         !contraction loops
         do j=1,np_b
            do i=1,np_a
               prod = ccf_a(i)*ccf_b(j)*norms_a(i)*norms_b(j)
               alp_ab = alp_a(i)+alp_b(j)
               arg2 = r_sq*alp_ab
               incomplete_gamma = cfp_gamic(arg,arg2)
               olap_tail = olap_tail + prod*fac*0.5_cfp/(alp_ab)**(arg)*incomplete_gamma
               incomplete_gamma2 = cfp_gamic(arg3,arg2)
               kei_tail = kei_tail + prod*fac*alp_b(j)*(-(3+2*l)/(alp_ab)**(arg)*incomplete_gamma &
                    + (alp_ab)**(-arg3)*2*alp_b(j)*incomplete_gamma2)

               bloch_ab = bloch_ab + fac_bloch*prod*exp(-alp_ab*r_sq)*(1+l-2*alp_b(j)*r_sq)
            enddo 
         enddo 

         !multiply the result by the norms of the contracted GTOs and by -1/2 for the kinetic energy
         olap_tail = olap_tail*cnorm_a*cnorm_b
         kei_tail = -0.5_cfp*kei_tail*cnorm_a*cnorm_b
         bloch_ab =  0.5_cfp*bloch_ab*cnorm_a*cnorm_b

   end subroutine olap_kei_tail

   !> Multipole property tail integrals for a pair of shells of contracted spherical GTOs centered on the CMS. We assume that the multipole is centered on CMS just like the two GTO shells.
   !> \f[
   !>    M_{tail} = \langle N_{a} S_{La,Ma}(\mathbf{r})\exp[-\alpha r^2]\vert S_{Lc,Mc} \vert N_{b} S_{Lb,Mb}(\mathbf{r})\exp[-\beta r^2] \rangle = 
   !>             = \frac{4\pi^{3/2}}{\sqrt{(2La+1)(2Lb+1)(2Lc+1)}}N_{a}N_{b}\frac{1}{2}(\alpha+\beta)^{-(La+Lb+Lc+3)/2}\Gamma[(La+Lb+Lc+3)/2,a^2(\alpha+\beta)],
   !> \f]
   !> where the first integration is only over \$ r \ge a \$.
   !> The formula is for a pair of normalized primitives. The actual routine performs contraction over the primitives. The calculated integrals are stored in the output 1D array which emulates the 
   !> 3D array with dimensions (2*la+1,2*lb+1,(l_c+1)**2) where la = max(l_a,l_b), lb = min(l_a,l_b). Hence the tail integral for the M angular numbers ma,mb,m is stored in (ma+la+1,mb+lb+1,l*l+m+l+1), i.e.
   !> in the element of the 1D array with index ma+la+1 + (2*la+1)*(mb+lb) + (2*la+1)*(2*lb+1)*(l*l+m+l). Note that the order of the tail integrals as output by this routine is the same as the one in 
   !> sph_mult_mom so subtraction of the tails is straightforward: sph_mult_mom(:) - prop_tail(:).
   subroutine prop_cms_tail(l_c,l_a,l_b,np_a,np_b,alp_a,alp_b,ccf_a,ccf_b,cnorm_a,norms_a,cnorm_b,norms_b,rmat_radius,prop_tail)
      use special_functions, only: cfp_gamic
      use phys_const, only: fourpi
      implicit none
      real(kind=cfp), intent(in) :: rmat_radius, alp_a(:), alp_b(:), ccf_a(:), ccf_b(:), norms_a(:), norms_b(:), cnorm_a, cnorm_b
      integer, intent(in) :: l_c, l_a, l_b, np_a, np_b
      real(kind=cfp), intent(out) :: prop_tail(:)

      real(kind=cfp) :: arg, arg2, fac, prod, incomplete_gamma, alp_ab, r_sq, rad, f
      integer :: i, j, m, ma, mb, l, ind, la, lb

         !precalculate the couplings if needed so that we retrieve them efficiently
         call cpl%prec_cgaunt(l_a+l_b) !l_a+l_b is the maximum allowed L in the real Gaunt coefficient <la ma|l m|lb mb> that we use below.

         !This ensures that the calculated integrals will be ordered in the same way as in sph_mult_mom. We don't need to reorder the input shell data since the formula for the tail integral is 
         !invariant with respect to exchange of the two shells. Therfore we need to order only the angular momenta l_a, l_b.
         la = max(l_a,l_b)
         lb = min(l_a,l_b)

         fac = sqrt(fourpi)*fourpi/sqrt((2*l_a+1.0_cfp)*(2*l_b+1.0_cfp))*cnorm_a*cnorm_b
         r_sq = rmat_radius*rmat_radius
         prop_tail = 0.0_cfp

         !contraction loops
         do j=1,np_b
            do i=1,np_a
               prod = ccf_a(i)*ccf_b(j)*norms_a(i)*norms_b(j)
               alp_ab = alp_a(i)+alp_b(j)
               arg2 = r_sq*alp_ab

               ind = 0
               do l=0,l_c
                  f = fac/sqrt(2*l+1.0_cfp)
                  arg = 0.5_cfp*(l+l_a+l_b+3)
                  incomplete_gamma = cfp_gamic(arg,arg2)
                  rad = prod*f*0.5_cfp/(alp_ab)**(arg)*incomplete_gamma !the complete radial integral including all norms and contraction coefficients
                  !loop over the angular dependencies
                  do m=-l,l
                     do mb=-lb,lb
                        do ma=-la,la
                           ind = ind + 1 != ma+la+1 + (2*la+1)*(lb+mb) + (2*la+1)*(2*lb+1)*(l*l+l+m), i.e. prop_tail is a 3D array with dimensions (2*la+1,2*lb+1,(l+1)**2).
                           prop_tail(ind) = prop_tail(ind) + rad*cpl%rgaunt(l,la,lb,m,ma,mb)
                        enddo !ma
                     enddo !mb
                  enddo !m
               enddo !l
            enddo 
         enddo 

   end subroutine prop_cms_tail

   !> Calculates nuclear attraction tail integrals for a pair of shells of spherical CGTOs centered on the CMS.
   !> \f[
   !>    NAI_{tail} = (-1)^{m_{a}+m_{b}}\frac{16\pi^2}{\sqrt{(2l_{a}+1)(2l_{b}+1)}}N_{a}N_{b}\sum_{i=1}^{n_{a}}\sum_{j=1}^{n_{j}}c_{i}c_{j}np_{i}np_{j}\sum_{l=\vert l_a-l_b\vert}^{l_a+l_b}
   !>                 \frac{R_{c}^{l}}{2l+1}\frac{1}{2}(\alpha_{i}+\beta_{j})^{\frac{l_a+l_b-l+2}{2}}\Gamma[\frac{l_a+l_b-l+2}{2},a^{2}(\alpha_{i}+\beta_{j})]
   !>                 \sum_{m=-l}^{l}\langle l_a,m_a\vert l,m\vert l_b,m_b\rangle_{R}X_{l,m}(\mathbf{R_{c}}).
   !> \f]
   !> \todo Make xc,yc,zc arrays so I can calculate all NAI tails at once saving CPU time since all terms except Rc and Xlm(Rc) are the same!
   subroutine nari_tail(xc,yc,zc,l_a,l_b,np_a,np_b,alp_a,alp_b,ccf_a,ccf_b,cnorm_a,norms_a,cnorm_b,norms_b,rmat_radius,na_tail)
      use special_functions, only: cfp_gamic, cfp_resh
      use phys_const, only: fourpi
      implicit none
      real(kind=cfp), intent(in) :: rmat_radius, alp_a(:), alp_b(:), ccf_a(:), ccf_b(:), &
                                    norms_a(:), norms_b(:), cnorm_a, cnorm_b, xc, yc, zc
      integer, intent(in) :: l_a, l_b, np_a, np_b
      real(kind=cfp), intent(out) :: na_tail(:)

      real(kind=cfp) :: fac, prod, alp_ab, s_l, r_sq, arg2, arg, SH(-(l_a+l_b):(l_a+l_b),0:l_a+l_b+1), &
                        angular, incomplete_gamma, radial, d, dist, base, rt
      integer :: i, j, l, l_min, l_max, m, ind, m_a, m_b

         l_min = abs(l_a-l_b)
         l_max = l_a+l_b

         !precalculate the couplings if needed so that we retrieve them efficiently
         call cpl%prec_cgaunt(l_max)

         dist = sqrt(xc*xc+yc*yc+zc*zc)

         if (dist .eq. 0) then !then this is almost exactly the same as the tail for the overlap integral

            if (l_a .ne. l_b) then
               na_tail(:) = 0.0_cfp
               return
            endif

            fac = fourpi/(2*l_a+1.0_cfp)*cnorm_a*cnorm_b
            arg = 1.0_cfp+l_a
            r_sq = rmat_radius*rmat_radius
            na_tail(:) = 0.0_cfp
   
            !contraction loops
            do j=1,np_b
               do i=1,np_a
                  prod = ccf_a(i)*ccf_b(j)*norms_a(i)*norms_b(j)
                  alp_ab = alp_a(i)+alp_b(j)
                  arg2 = r_sq*alp_ab
                  incomplete_gamma = cfp_gamic(arg,arg2)
                  do m_a=-l_a,l_a
                     ind = l_a+m_a+1  + (2*l_a+1)*(m_a+l_a+1 -1)
                     na_tail(ind) = na_tail(ind) + prod*fac*0.5_cfp/(alp_ab)**(arg)*incomplete_gamma
                  enddo
               enddo 
            enddo 
   
         else

            fac = fourpi*fourpi/sqrt((2*l_a+1.0_cfp)*(2*l_b+1.0_cfp))*cnorm_a*cnorm_b
            r_sq = rmat_radius*rmat_radius

            base = dist**(l_max+1)
            dist = 1.0_cfp/dist
            call cfp_resh(SH,xc,yc,zc,l_max)
   
            na_tail(:) = 0.0_cfp
   
            !contraction loops
            do j=1,np_b
               do i=1,np_a
                  prod = fac*ccf_a(i)*ccf_b(j)*norms_a(i)*norms_b(j)
                  alp_ab = alp_a(i)+alp_b(j)
                  arg2 = r_sq*alp_ab
                  rt = 1.0_cfp/sqrt(alp_ab)
   
                  d = base
                  s_l = (l_a+l_b-(l_max+1)+2)*0.5_cfp
                  arg = 0.5_cfp*alp_ab**(-s_l)
                  do l=l_max,l_min,-1 !sum starting with the smallest numbers
                     d = d*dist         != dist**l
                     s_l = s_l + 0.5_cfp != (l_a+l_b-l+2.0_cfp)*0.5_cfp
                     arg = arg*rt       != 0.5_cfp/(alp_ab**s_l)
                     incomplete_gamma = cfp_gamic(s_l,arg2)
                     radial = d/(2*l+1.0_cfp)*arg*incomplete_gamma*prod
   
                     ind = 0
                     do m_b=-l_b,l_b
                        do m_a=-l_a,l_a
                           ind = ind + 1 != l_a+m_a+1  + (2*l_a+1)*(m_b+l_b+1 -1)

                           do m=-l,l
                              angular = cpl%rgaunt(l_a,l,l_b,m_a,m,m_b)
                              if (angular .ne. 0.0_cfp) then
                                 na_tail(ind) = na_tail(ind) + radial*angular*SH(m,l)
                                 !write(*,'("contr",5i4,2e)') ind, l,m_a,m,m_b, radial*angular, na_tail(ind)
                              endif
                           enddo !m
                        enddo !m_a
                     enddo !m_b
   
                  enddo !l
   
               enddo !i
            enddo !j

         endif

   end subroutine nari_tail

   !> Calculates 2-electron tail integrals for 1 electron in the continuum for a pair of shells of spherical CGTOs centered on the CMS and a pair of shells of target CGTOs. This is the driver routine for
   !> eri_tail_shell. This routine ensures that the parameters are passed to eri_tail_shell in such an order to ensure that the tail integrals on output are ordered in the same way as the integrals
   !> calculated by eri. This allows for easy subtraction of the tails. l_a_tgt,l_b_tgt are angular momenta in the pair of shells of target CGTOs, tgt_prop and tgt_pair specify the property integrals for
   !> the pair of shells of the target CGTOs. The rest of the input parameters are related to the continuum and the R-matrix radius.
   subroutine eri_tail(tgt_prop,tgt_pair,l_a_tgt,l_b_tgt,l_a,l_b,np_a,np_b,alp_a,&
                       alp_b,ccf_a,ccf_b,norms_a,norms_b,rmat_radius,swap_ab_cd,eri_tail_int)
      implicit none
      real(kind=cfp), intent(in) :: tgt_prop(:,:), rmat_radius, alp_a(:), alp_b(:), ccf_a(:), ccf_b(:), norms_a(:), norms_b(:)
      integer, intent(in) :: l_a, l_b, np_a, np_b, l_a_tgt, l_b_tgt, tgt_pair
      real(kind=cfp), allocatable :: eri_tail_int(:)
      logical, intent(in) :: swap_ab_cd

        call eri_tail_shell(tgt_prop,tgt_pair,l_a_tgt,l_b_tgt,l_a,l_b,np_a,np_b,alp_a,&
                            alp_b,ccf_a,ccf_b,norms_a,norms_b,rmat_radius,swap_ab_cd,eri_tail_int)

   end subroutine eri_tail

   !> Calculates 2-electron tail integrals for 1 electron in the continuum for a pair of shells of spherical CGTOs centered on the CMS and a pair of shells of target CGTOs.
   !> \f[
   !>    ERI_{tail} = (-1)^{m_{a}+m_{b}}\frac{4\pi}{\sqrt{(2l_{a}+1)(2l_{b}+1)}}N_{a}N_{b}\sum_{i=1}^{n_{a}}\sum_{j=1}^{n_{j}}c_{i}c_{j}np_{i}np_{j}\sum_{l=\vert l_a-l_b\vert}^{l_a+l_b}
   !>                 \sqrt{\frac{4\pi}{2l+1}}\frac{1}{2}(\alpha_{i}+\beta_{j})^{\frac{l_a+l_b-l+2}{2}}\Gamma[\frac{l_a+l_b-l+2}{2},a^{2}(\alpha_{i}+\beta_{j})]
   !>                 \sum_{m=-l}^{l}\langle l_a,m_a\vert l,m\vert l_b,m_b\rangle_{R}P_{A,B}^{l,m}.
   !> \f]
   !> This routine calculates the tail integrals in the order (A,B,C,D) where A,B are the continuum shells. The variable swap_ab_cd controls whether the integrals on output should be ordered as (C,D,A,B).
   !> The array tgt_prop is the array of the properties calculated by sph_mult_mom at least for 
   !> L_max = l_a+l_b for the pair of shells (ab) of the target GTOs. However, here we assume that the values in tgt_prop have been calculated for a number of pairs of shells of GTOs and that these are saved
   !> in columns of tgt_prop, where the column corresponding to the current combination of target GTOs is given by the index tgt_pair.
   !> \todo The whole eri tail calculation can be improved on since all of the couplings are the same as for NARI hence the eri tails should be precalculated reverting the loop so that the properties are
   !> in the inner loop and the continuum functions the outer loop. This should speed things up considerably.
   subroutine eri_tail_shell(tgt_prop,tgt_pair,l_a_tgt,l_b_tgt,l_a,l_b,np_a,np_b,&
                             alp_a,alp_b,ccf_a,ccf_b,norms_a,norms_b,rmat_radius,swap_ab_cd,eri_tail_int)
      use utils, only: xermsg
      use special_functions, only: cfp_gamic
      use phys_const, only: fourpi
      implicit none
      real(kind=cfp), intent(in) :: tgt_prop(:,:), rmat_radius, alp_a(:), alp_b(:), ccf_a(:), ccf_b(:), norms_a(:), norms_b(:)
      integer, intent(in) :: l_a, l_b, np_a, np_b, l_a_tgt, l_b_tgt, tgt_pair
      real(kind=cfp), allocatable :: eri_tail_int(:)
      logical, intent(in) :: swap_ab_cd

      real(kind=cfp) :: fac, prod, alp_ab, s_l, r_sq, arg2, arg, angular, incomplete_gamma, radial, inv
      integer :: i, j, l, l_min, l_max, m, ind, m_a, m_b, ind_prop, m_a_tgt, m_b_tgt, sph_a_tgt, sph_b_tgt, sph_a_cont, sph_b_cont
      integer :: err

         l_min = abs(l_a-l_b)
         l_max = l_a+l_b

         !precalculate the couplings if needed so that we retrieve them efficiently
         call cpl%prec_cgaunt(l_max)

         fac = fourpi/sqrt((2*l_a+1.0_cfp)*(2*l_b+1.0_cfp))
         r_sq = rmat_radius*rmat_radius

         sph_a_tgt = 2*l_a_tgt+1
         sph_b_tgt = 2*l_b_tgt+1
         sph_a_cont = 2*l_a+1
         sph_b_cont = 2*l_b+1
         ind = sph_a_cont*sph_b_cont*sph_a_tgt*sph_b_tgt

         if (.not. allocated(eri_tail_int) .or. size(eri_tail_int) < ind) then
            if (allocated(eri_tail_int)) deallocate(eri_tail_int)
            allocate(eri_tail_int(ind),stat=err)
            if (err .ne. 0) call xermsg('gto_routines','eri_tail_shell','Memory allocation failed.',err,1)
         endif

         eri_tail_int(1:ind) = 0.0_cfp

         !contraction loops
         do j=1,np_b
            do i=1,np_a
               prod = fac*ccf_a(i)*ccf_b(j)*norms_a(i)*norms_b(j)
               alp_ab = alp_a(i)+alp_b(j)
               arg2 = r_sq*alp_ab
               inv = 1.0_cfp/alp_ab
   
               s_l = (l_a+l_b-(l_max+2)+2)*0.5_cfp
               arg = 0.5_cfp*alp_ab**(-s_l)
               do l=l_max,l_min,-2 !sum starting with the smallest contributions to the tail; we loop over l in steps of 2 since the sum L values entering rgaunt must be even and l_max+l_a+l_b = even.
                  s_l = s_l + 1.0_cfp  != (l_a+l_b-l+2.0_cfp)*0.5_cfp
                  arg = arg*inv       != 0.5_cfp/(alp_ab**s_l)
                  incomplete_gamma = cfp_gamic(s_l,arg2)
                  radial = arg*incomplete_gamma*prod*sqrt(fourpi/(2*l+1.0_cfp))

                  !We assume that the property integrals are indexed as if tgt_prop was a 3D array with dimensions (2*l_a_tgt+1,2*l_b_tgt+1,(l_max+1)**2)
                  ind_prop = sph_a_tgt*sph_b_tgt*(l*l) !now we need property integrals for L=l
                  do m=-l,l 
                     ind = 0 !the tail integrals are indexed as if eri_tail_int was a 4D array with dimensions (2*l_a+1,2*l_b+1,2*l_a_tgt+1,2*l_b_tgt+1)
                     do m_b_tgt = -l_b_tgt,l_b_tgt
                        do m_a_tgt = -l_a_tgt,l_a_tgt
                           ind_prop = ind_prop + 1
                           if (tgt_prop(ind_prop,tgt_pair) .ne. 0.0_cfp) then !we can have a non-zero contribution only if the target property is non-zero
                              ind = sph_a_cont*sph_b_cont*(m_a_tgt+l_a_tgt) + sph_a_cont*sph_b_cont*sph_a_tgt*(m_b_tgt+l_b_tgt)
                              do m_b=-l_b,l_b
                                 do m_a=-l_a,l_a

                                    ind = ind + 1

                                    angular = cpl%rgaunt(l_a,l,l_b,m_a,m,m_b)
                                    if (angular .ne. 0.0_cfp) then
                                       eri_tail_int(ind) = eri_tail_int(ind) + radial*angular*tgt_prop(ind_prop,tgt_pair)
                                    endif

                                 enddo !m_a
                              enddo !m_b
                           endif
                        enddo !m_b_tgt
                     enddo !m_a_tgt
                  enddo !m
   
               enddo !l
   
            enddo !i
         enddo !j

         !Above we calculate the tail integrals assuming the order [AB|CD], where AB are the continuum functions. If the actual order of the integrals on output should be [CD|AB] then we have to transpose 
         !the array eri_tail_int.
         if (swap_ab_cd) then

            ind = sph_a_cont*sph_b_cont*sph_a_tgt*sph_b_tgt
            i = check_real_array_size(eri_tail_tmp,ind)
            if (i .ne. 0) call xermsg('gto_routines','eri_tail_shell','Error in allocating space by check_real_array_size.',i,1)

            eri_tail_tmp(1:ind) = eri_tail_int(1:ind)

            call abcd_to_cdab(eri_tail_tmp,eri_tail_int,sph_a_cont,sph_b_cont,sph_a_tgt,sph_b_tgt)

         endif

   end subroutine eri_tail_shell

   !> This routine takes on input the batch of integrals [ab|cd] with dimensions given by na,nb,nc,nd and returns in batch_t the transposed batch [cd|ab].
   !> Note that this procedure is equivalent to transposition of a matrix M with dimensions na*nb,nc*nd.
   subroutine abcd_to_cdab(batch_in,batch_t,na,nb,nc,nd)
      use const, only: tile
      implicit none
      integer, intent(in) :: na,nb,nc,nd
      real(kind=cfp), intent(in) :: batch_in(*)
      real(kind=cfp), intent(out) :: batch_t(*)

      integer :: nrow, ncol, i, j, ii, jj, iend, jend, col

         nrow = na*nb
         ncol = nc*nd
         do jj = 1,nrow,tile
            jend = min (nrow,jj+tile-1)
            do ii = 1,ncol,tile
               iend = min (ncol,ii+tile-1)

               do j = jj,jend
                  col = ncol*(j-1)
                  do i = ii,iend
                     batch_t(i+col) = batch_in(j+nrow*(i-1))
                  enddo
               enddo

            enddo
         enddo
      
   end subroutine abcd_to_cdab

  !> Takes on input the linear array 'a' and the dimension 'd' which the array should AT LEAST have, i.e. not exactly. If the array is smaller it is reallocated (loosing its contents) to the size 'd'.
  function check_real_array_size(a,d)
     implicit none
     integer, intent(in) :: d
     real(kind=cfp), allocatable, intent(inout) :: a(:)
     integer :: check_real_array_size

         check_real_array_size = 0

         if (d .le. 0) stop "d .le. 0 in check_real_array_size"

         if (.not.allocated(a)) then
            allocate(a(d),stat=check_real_array_size)
         elseif (size(a) < d) then
            deallocate(a,stat=check_real_array_size)
            if (check_real_array_size .ne. 0) stop "error deallocating array in check_real_array_size"
            allocate(a(d),stat=check_real_array_size)
         endif

   end function check_real_array_size

   !> Calculates indices of 1-electron integrals output in the linear arrays from the integral routines given the starting indices of functions in the a and b shells.
   !> ind_a corresponds to the function with m=-la, same for ind_b. We assume that indices of the rest of the functions in the shell are sequential.
   !> The indices for each integral are reordered in such a way that the index of the a function is always .ge. index of the b function. Since the 1-electron integrals are symmetrical there is no problem.
   !> The number np denotes the number of passive indices which don't enter the ordering, i.e. we think of the integral batch as being a 3d array (2*la+1,2*lb+1,np). The third dimension useded to be used in the
   !> ordering of the multipole moment integrals where the first two dimensions run over the GTOs A,B and the third dimension over the components of the multipole moments. The requirement to index all
   !> multipole integrals separately has been dropped so this option is effectively not being used in the whole library.
   subroutine index_1el(la,lb,ind_a,ind_b,np,int_index)
      implicit none
      integer, intent(in) :: la, lb, ind_a, ind_b, np
      integer, intent(out) :: int_index(:,:) !indices of the integrals corresponding to the functions in the a=1, and b=2 shells.

      integer :: a,b,i,p,sph_shell_a

      i = 0
      sph_shell_a = 2*la+1
      do p=1,np !this is the loop over the passive indices, e.g. what used to be indices of the multipole moments
         do b=ind_b,ind_b+2*lb
            do a=ind_a,ind_a+2*la !the a shell index changes quickest, then b
               i = i + 1
               if (a .ge. b) then
                  int_index(1,i) = a
                  int_index(2,i) = b
               else
                  int_index(1,i) = b
                  int_index(2,i) = a
               endif
            enddo
         enddo
      enddo

   end subroutine index_1el

   !> Calculates indices corresponding to the 2-electron integrals (ab|cd) in the linear array sph_ints (calculated by the routine eri_shell) given the starting indices of the functions in the 
   !> a,b,c,d shells. ind_a must correspond to the function with m=-la, same for the b,c,d shells. We assume that indices of the rest of the functions in each shell are sequential.
   !> The indices of the functions are ordered, exploiting the permutational symmetry of the (ab|cd) integrals, so that:
   !> If keep_ab_cd_order .eq. .false. then the input variable swap_ab_cd is ignored and:
   !> ind_a .ge. ind_b, ind_c .ge. ind_d, ind_a+ind_b .ge. ind_c+ind_d, ind_a .ge. ind_c. This order of the functions is needed for the integral indexing scheme.
   !> If keep_ab_cd_order .eq. .true.: indices of the pairs of shells (ab|, |cd) are swapped or not depending on the value of the logical variable swap_ab_cd. Within the (ab| and |cd) shells the indices are
   !> ordered so that a.ge.b, c.ge.d.
   subroutine index_2el(la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,int_index,keep_ab_cd_order,swap_ab_cd)
      implicit none
      integer, intent(in) :: la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d
      !intent(out):
      integer, allocatable :: int_index(:,:) !indices of the integrals corresponding to the functions in the a,b,c,d shells.
      logical, intent(in) :: keep_ab_cd_order, swap_ab_cd

      integer :: a,b,c,d,i,iAB,iCD,t,ct,dt

         if (keep_ab_cd_order) then !use the swap_ab_cd logical to swap indices if necessary
            i = 0
            iCD = 0; iAB = 0
            do d=ind_d,ind_d+2*ld
               do c=ind_c,ind_c+2*lc
                  if (c .ge. d) then
                     ct = c
                     dt = d
                  else
                     ct = d
                     dt = c
                  endif
                  do b=ind_b,ind_b+2*lb
                     do a=ind_a,ind_a+2*la !the a shell index in sph_ints changes quickest, then b,c,d
                        i = i + 1
                        if (a .ge. b) then
                           int_index(1,i) = a
                           int_index(2,i) = b
                        else
                           int_index(1,i) = b
                           int_index(2,i) = a
                        endif
                        int_index(3,i) = ct
                        int_index(4,i) = dt
                        if (swap_ab_cd) then !todo this can be taken outside of the loops and the code can be unrolled into two options
                           t = int_index(1,i)
                           int_index(1,i) = int_index(3,i)
                           int_index(3,i) = t
      
                           t = int_index(2,i)
                           int_index(2,i) = int_index(4,i)
                           int_index(4,i) = t
                        endif
                        !At this point the quartet of indices has been ordered so that the order of the pairs is either (ab|cd) or (cd|ab) depending on swap_ab_cd.
                     enddo !a
                  enddo !b
               enddo !c
            enddo !d
         else !Permute the indices in the order required for full indexing of the integrals
            i = 0
            do d=ind_d,ind_d+2*ld
               do c=ind_c,ind_c+2*lc
                  if (c .ge. d) then
                     ct = c
                     dt = d
                  else
                     ct = d
                     dt = c
                  endif
                  do b=ind_b,ind_b+2*lb
                     do a=ind_a,ind_a+2*la !the a shell index in sph_ints changes quickest, then b,c,d
                        i = i + 1
                        if (a .ge. b) then
                           int_index(1,i) = a
                           int_index(2,i) = b
                        else
                           int_index(1,i) = b
                           int_index(2,i) = a
                        endif
                        int_index(3,i) = ct
                        int_index(4,i) = dt
                        !compute indices of the AB, CD pairs and see if wee have to swap; the pair index corresponds to the index of within a set of unique pairs of functions (ab) where a.ge.b.
                        iAB = int_index(1,i)*(int_index(1,i)-1)/2+int_index(2,i)
                        iCD = int_index(3,i)*(int_index(3,i)-1)/2+int_index(4,i)
                        if (iAB < iCD) then
                           t = int_index(1,i)
                           int_index(1,i) = int_index(3,i)
                           int_index(3,i) = t
      
                           t = int_index(2,i)
                           int_index(2,i) = int_index(4,i)
                           int_index(4,i) = t
                        endif
                     enddo !a
                  enddo !b
               enddo !c
            enddo !d
         endif

   end subroutine index_2el

   !> Orders the shells according to the angular momentum in the same way it is done in eri and then computes the indices for each quartet of (ab|cd) integrals. See index_2el for 
   !> description of the ordering option given by the input variable keep_ab_cd_order.
   subroutine index_2el_drv(la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,ind,keep_ab_cd_order)
      implicit none
      integer, intent(in) :: la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d
      !intent(out):
      integer, allocatable :: ind(:,:)
      logical, intent(in) :: keep_ab_cd_order

         !Order the shells according to their angular momentum: la .ge .lb, lc .ge. ld, la+lb .ge. lc+ld and compute the full integral indices
         if (la+lb .ge. lc+ld) then
            if (la .ge. lb) then
               if (lc .ge. ld) then !la,lb,lc,ld
                  call index_2el(la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,ind,keep_ab_cd_order,.false.)
               else !la,lb,ld,lc
                  call index_2el(la,lb,ld,lc,ind_a,ind_b,ind_d,ind_c,ind,keep_ab_cd_order,.false.)
               endif
            else !la < lb
               if (lc .ge. ld) then !lb,la,lc,ld
                  call index_2el(lb,la,lc,ld,ind_b,ind_a,ind_c,ind_d,ind,keep_ab_cd_order,.false.)
               else !lb,la,ld,lc
                  call index_2el(lb,la,ld,lc,ind_b,ind_a,ind_d,ind_c,ind,keep_ab_cd_order,.false.)
               endif
            endif
         else !la+lb < lc+ld
            if (la .ge. lb) then
               if (lc .ge. ld) then !lc,ld,la,lb
                  call index_2el(lc,ld,la,lb,ind_c,ind_d,ind_a,ind_b,ind,keep_ab_cd_order,.true.)
               else !ld,lc,la,lb
                  call index_2el(ld,lc,la,lb,ind_d,ind_c,ind_a,ind_b,ind,keep_ab_cd_order,.true.)
               endif
            else !la < lb
               if (lc .ge. ld) then !lc,ld,lb,la
                  call index_2el(lc,ld,lb,la,ind_c,ind_d,ind_b,ind_a,ind,keep_ab_cd_order,.true.)
               else !ld,lc,lb,la
                  call index_2el(ld,lc,lb,la,ind_d,ind_c,ind_b,ind_a,ind,keep_ab_cd_order,.true.)
               endif
            endif
         endif

   end subroutine index_2el_drv

   !> Determine the mapping of indices ind_a,ind_b,ind_c,ind_d in the array ind_orig so that ind_ap.ge.ind_bp, ind_cp.ge.ind_dp, ind_ap.ge.ind_cp.
   subroutine find_mapping(ind_orig,n,n_map,map)
      implicit none
      integer, intent(inout) :: ind_orig(4)
      integer, intent(in) :: n(4)
      integer, intent(out) :: n_map(3), map(4)

      integer :: i, j, max_AB, max_CD

         map(1:4) = (/1,2,3,4/)

         max_AB = max(ind_orig(1),ind_orig(2))
         if (max_AB > ind_orig(1)) then !swap A,B
            max_AB = ind_orig(1)
            ind_orig(1) = ind_orig(2)
            ind_orig(2) = max_AB
            map(1) = 2
            map(2) = 1
         endif

         max_CD = max(ind_orig(3),ind_orig(4))
         if (max_CD > ind_orig(3)) then !swap C,D
            max_CD = ind_orig(3)
            ind_orig(3) = ind_orig(4)
            ind_orig(4) = max_CD
            map(3) = 4
            map(4) = 3
         endif

         i = ipair(ind_orig(1))+ind_orig(2)
         j = ipair(ind_orig(3))+ind_orig(4)

         if (j > i) then !swap AB and CD
            i = map(1)
            j = map(2)
            map(1) = map(3)
            map(2) = map(4)
            map(3) = i
            map(4) = j
         endif

         n_map(1) = n(map(1))
         n_map(2) = n(map(2))*n_map(1)
         n_map(3) = n(map(3))*n_map(2)

   end subroutine find_mapping

   !> Reorder the integrals ints in columns a,b,c,d so that ap.ge.bp,cp.ge.dp,ap.ge.cp. The input values are the angular momenta and the starting indices in columns 1,2,3,4 of the ints and ind arrays.
   !> On output the indices in ind(1:4,i) are not ordered from the largest to the smallest as the original ones (a.ge.b,c.ge.d,a.ge.c) but that's OK since they're not needed to save the integrals
   !> in the MPI mode when this routine is used.
   subroutine reorder_and_index_2el(la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,column,ind,ints)
      implicit none
      integer, intent(in) :: la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,column
      !intent(out):
      integer, allocatable :: ind(:,:)
      real(kind=cfp), allocatable :: ints(:,:)

      integer :: ind_ordered, ind_orig(4), i, n(4), n_map(3), map(4), i1,i2,i3,i4, a,b,c,d, n_ints, err
      logical :: ordered

         ind_orig(1:4) = (/ind_a,ind_b,ind_c,ind_d/)
         n(1:4) = (/2*la+1,2*lb+1,2*lc+1,2*ld+1/)

         call find_mapping(ind_orig,n,n_map,map)

         !Check if the required order and the present order are the same. If they are then there is nothing to do.
         ordered = .true.
         do i=1,4
            if (map(i) .ne. i) ordered = .false.
         enddo
         if (ordered) return

         n_ints = n(1)*n(2)*n(3)*n(4)
         if (allocated(ints_tmp)) then
            if (size(ints_tmp) < n_ints) deallocate (ints_tmp)
         end if
         if (.not. allocated(ints_tmp)) then
            allocate(ints_tmp(n_ints),stat=err)
            if (err .ne. 0) call xermsg ('gto_routines', 'reorder_and_index_2el', 'Memory allocation 2 has failed.', err, 1)
         endif

         !Reorder the integrals so that 1st column corresponds to the shell with
         !the largest starting index, 2nd column corresponds to the shell with
         !the 2nd largest starting index, etc.
         ints_tmp(1:n_ints) = ints(1:n_ints,column)
         i = 0
         do d=1,n(4)
            do c=1,n(3)
               do b=1,n(2)
                  do a=1,n(1)
                     i = i + 1
                     ind_orig(1:4) = (/a,b,c,d/)
                     i1 = ind_orig(map(1))
                     i2 = ind_orig(map(2))
                     i3 = ind_orig(map(3))
                     i4 = ind_orig(map(4))
                     ind_ordered = i1 + n_map(1)*(i2-1) + n_map(2)*(i3-1) + n_map(3)*(i4-1)
                     ind(1:4,ind_ordered) = (/i1,i2,i3,i4/)
                     ints(ind_ordered,column) = ints_tmp(i)
                  enddo !a
               enddo !b
            enddo !c
         enddo !d

   end subroutine reorder_and_index_2el

   !> This is an index for two-electron integral keeping only 1p in the continuum. It works also for 2p in the continuum in which case is_CCTT must be set to .false.
   pure function index_1p_continuum(ordered_pairs,ind1,ind2,ind3,ind4,is_CCTT,last_CT,n_prec,n_TT_pairs)
      use special_functions, only: ipair
      implicit none
      integer :: index_1p_continuum
      integer, intent(in) :: ordered_pairs(:,:)
      integer, intent(in) :: ind1,ind2,ind3,ind4,last_CT,n_prec,n_TT_pairs
      logical, intent(in) :: is_CCTT

      integer :: pq,rs

         pq = ordered_pairs(1,ipair(ind1) + ind2)
         rs = ordered_pairs(1,ipair(ind3) + ind4)

         index_1p_continuum = 0

         if (is_CCTT) then
            if (pq > last_CT) then !pq is CC pair
               index_1p_continuum = n_prec + rs + n_TT_pairs*(pq-last_CT-1)
            else !rs is CC pair
               index_1p_continuum = n_prec + pq + n_TT_pairs*(rs-last_CT-1)
            endif
         else
            index_1p_continuum = ipair(max(pq,rs)) + min(pq,rs)
         endif

   end function index_1p_continuum

   subroutine normalize_cgto(number_of_primitives,l,exponents,contractions,norms,norm)
      use phys_const, only: twopi
      implicit none
      integer :: number_of_primitives,l
      real(kind=cfp), intent(in) :: exponents(number_of_primitives), contractions(number_of_primitives)
      real(kind=cfp), intent(out) :: norms(number_of_primitives), norm

      integer :: i,j,l2p1
      real(kind=cfp) :: olap, tmp, n_i, n_j

         do i=1,number_of_primitives  !calculate the norms of the primitive GTO. We normalize them to 1.
            norms(i) = dngto(l,exponents(i))
         enddo

         !Calculate the normalization factor for the contracted GTO.
         !This can of course be done more efficiently, but the loops below do the job...
         olap=0.0_cfp !self-overlap of the contracted GTO
         l2p1 = 2*l+1.0_cfp
         tmp = l+1.5_cfp
         do i=1,number_of_primitives
            n_i = norms(i) !normalization factor for the primitive GTO 
            olap = olap + contractions(i)*contractions(i)*n_i*n_i*twopi/l2p1*cfp_gamma_fun(tmp)/(2*exponents(i))**(tmp)
            do j=i+1,number_of_primitives
               n_j = norms(j) !normalization factor for the primitive GTO
               olap = olap + 2.0_cfp*contractions(i)*contractions(j)*n_i*n_j*twopi/l2p1*cfp_gamma_fun(tmp) &
                    /(exponents(i)+exponents(j))**(tmp)
            enddo
         enddo

         norm = 1.0_cfp/sqrt(olap) !norm of the contracted GTO

   end subroutine normalize_cgto

   function check_cgto_data(number_of_primitives,l,exponents,contractions,norms,number_of_functions)
      use const, only: max_contr_len
      implicit none
      integer, intent(in) :: number_of_primitives, l, number_of_functions
      real(kind=cfp), allocatable :: exponents(:), contractions(:), norms(:)
      integer :: check_cgto_data

      integer :: i,j

         check_cgto_data = 0

         if (l < 0) then
            check_cgto_data = 1
            return
         endif

         if (number_of_primitives .le. 0 .or. number_of_primitives > max_contr_len) then
            check_cgto_data = 2
            return
         endif

         if (.not.(allocated(exponents)) .or. .not.(allocated(contractions)) .or. .not.(allocated(norms))) then
            check_cgto_data = 3
            return
         endif

         if (size(exponents) /= number_of_primitives .or. &
             size(contractions) /= number_of_primitives .or. &
             size(norms) /= number_of_primitives) then
            check_cgto_data = 4
            return
         endif

         !check for duplicities in the exponents and check their values
         do i=1,number_of_primitives
            if (exponents(i) .le. 0.0_cfp) then
               check_cgto_data = 5
            endif
            do j=1,i-1
               if (exponents(i) .eq. exponents(j)) then
                  check_cgto_data = 6
                  return
               endif
            enddo
         enddo

         if (number_of_functions .ne. 2*l+1) then
            check_cgto_data = 7
            return
         endif

   end function check_cgto_data

   subroutine read_cgto(number_of_primitives,l,exponents,contractions,norms,norm,center, &
                        non_zero_at_boundary,number_of_functions,lunit,posit,pos_after_rw)
      use mpi_mod
      implicit none
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: number_of_primitives, l, pos_after_rw, number_of_functions
      real(kind=cfp), allocatable :: exponents(:), contractions(:), norms(:)
      real(kind=cfp), intent(out) :: norm, center(3)
      logical, intent(out) :: non_zero_at_boundary

      integer :: err

         if (myrank .eq. master) then
            read(lunit,pos=posit,err=10) center(1:3)
            read(lunit,err=10) l, number_of_functions
            !now data on the primitives
            read(lunit,err=10) number_of_primitives

            if (allocated(exponents)) deallocate(exponents)
            if (allocated(contractions)) deallocate(contractions)
            if (allocated(norms)) deallocate(norms)
            !as long as the write_cgto method has been used to write the data we should never have number_of_primitives .le. 0 since that is checked in write_cgto before the data are written.
            allocate(exponents(1:number_of_primitives), &
                     contractions(1:number_of_primitives),norms(1:number_of_primitives),stat=err)
            if (err .ne. 0) call xermsg ('gto_routines', 'read_cgto', 'Memory allocation 1 has failed.', err, 1)
   
            read(lunit,err=10) exponents(1:number_of_primitives), &
                               contractions(1:number_of_primitives), &
                               norms(1:number_of_primitives), norm, non_zero_at_boundary
            inquire(lunit,pos=pos_after_rw)
         endif

         !master broadcasts all its data to the other processes
         !start with the integers (i.e. dimensions, etc.)
         !todo replace by one 6-element array broadcast
         call mpi_mod_bcast(l,master)
         call mpi_mod_bcast(number_of_functions,master)
         call mpi_mod_bcast(number_of_primitives,master)
         call mpi_mod_bcast(pos_after_rw,master)
         call mpi_mod_bcast(non_zero_at_boundary,master)
         call mpi_mod_bcast(norm,master)

         if (myrank .ne. master) then
            if (allocated(exponents)) deallocate(exponents)
            if (allocated(contractions)) deallocate(contractions)
            if (allocated(norms)) deallocate(norms)
            !as long as the write_cgto method has been used to write the data we should never have number_of_primitives .le. 0 since that is checked in write_cgto before the data are written.
            allocate(exponents(1:number_of_primitives),contractions(1:number_of_primitives),norms(1:number_of_primitives),stat=err)
            if (err .ne. 0) call xermsg ('gto_routines', 'read_cgto', 'Memory allocation 2 has failed.', err, 1)
         endif

         !broadcast the array data
         call mpi_mod_bcast(center,master)
         call mpi_mod_bcast(exponents,master)
         call mpi_mod_bcast(contractions,master)
         call mpi_mod_bcast(norms,master)

         return

 10      call xermsg ('gto_routines', 'read_cgto', 'Error reading the CGTO data from the file and position given.', 1, 1)

   end subroutine read_cgto

   subroutine write_cgto(number_of_primitives,l,exponents,contractions,norms,norm,center,&
                         non_zero_at_boundary,number_of_functions,lunit,posit,pos_after_rw)
      use mpi_mod
      implicit none
      integer, intent(in) :: number_of_primitives, l, lunit, posit, number_of_functions
      real(kind=cfp), intent(in) :: exponents(number_of_primitives), contractions(number_of_primitives),&
                                    norms(number_of_primitives), norm, center(3)
      logical, intent(in) :: non_zero_at_boundary
      integer, intent(out) :: pos_after_rw

         pos_after_rw = 0
         if (myrank .eq. master) then
            write(lunit,pos=posit,err=10) center(1:3)
            write(lunit,err=10) l, number_of_functions
            !now data on the primitives
            write(lunit,err=10) number_of_primitives
            write(lunit,err=10) exponents(1:number_of_primitives), contractions(1:number_of_primitives),&
                                norms(1:number_of_primitives), norm, non_zero_at_boundary

            inquire(lunit,pos=pos_after_rw)
         endif

         !master ensures all processes know where the record ends
         call mpi_mod_bcast(pos_after_rw,master)

         return

 10      call xermsg ('gto_routines', 'write_cgto', 'Error writing the GTO data into the file and position given.', 1, 1)

   end subroutine write_cgto

   subroutine print_cgto_data(number_of_primitives,l,exponents,contractions,norms,norm,center,non_zero_at_boundary)
      use const, only: stdout
      implicit none
      integer, intent(in) :: number_of_primitives, l
      real(kind=cfp), intent(in) :: exponents(number_of_primitives), contractions(number_of_primitives), &
                                    norms(number_of_primitives), norm, center(3)
      logical, intent(in) :: non_zero_at_boundary

      integer :: n

         write(stdout,'("Contracted spherical GTO shell data:")')
         write(stdout,'("L:",i0)') l
         write(stdout,'("Center: ",3e20.10)') center
         write(stdout,'("Number of primitives:",i0)') number_of_primitives
 
         n = size(exponents)
         if (n > 0) write(stdout,'("Exponents:   ",100e20.10)') exponents

         n = size(contractions)
         if (n > 0) write(stdout,'("Contractions:",100e20.10)') contractions

         n = size(norms)
         if (n > 0) write(stdout,'("Primitive normalization factors:",100e20.10)') norms

         write(stdout,'("Contraction normalization factor:",e20.10)') norm
         write(stdout,'("Is non-zero at the boundary: ",l)') non_zero_at_boundary

   end subroutine print_cgto_data

   subroutine eval_cgto(r,n_points,number_of_primitives,l,exponents,contractions,norms,norm,center,eval_CGTO_shell)
      use special_functions, only: cfp_solh
      implicit none
      integer, intent(in) :: number_of_primitives, l, n_points
      real(kind=cfp), intent(in) :: exponents(number_of_primitives), contractions(number_of_primitives), &
                                    norms(number_of_primitives), norm, center(3)
      real(kind=cfp), intent(in) :: r(1:3,n_points)  !coordinates with respect to the center of coordinates of the points for which we want to evaluate the GTO
      real(kind=cfp) :: eval_CGTO_shell(2*l+1,n_points)

      real(kind=cfp), allocatable :: slm(:,:)
      real(kind=cfp) :: dist, rdiff(1:3), rad
      integer :: m, err, i, p

         if (l > 0) then
            allocate(slm(-l:l,0:l),stat=err)
            if (err .ne. 0) call xermsg ('gto_routines', 'eval_cgto', 'Memory allocation failed', err, 1)
         endif

         do p=1,n_points

            rdiff(1:3) = r(1:3,p)-center(1:3)
   
            !accumulate the exponential part from the primitive Gaussians mulitplied by the contraction coefficients
            dist = dot_product(rdiff,rdiff)
            rad = 0.0_cfp
            do i=1,number_of_primitives
               rad = rad + contractions(i)*norms(i)*exp(-exponents(i)*dist)
            enddo
   
            if (l > 0) then
               call cfp_solh(slm,rdiff(1),rdiff(2),rdiff(3),l) !calculate the solid harmonics at r-center up to l (we'll need only slm(m,l))
   
               do m=-l,l
                  eval_CGTO_shell(m+l+1,p) = norm*slm(m,l)*rad
               enddo !m
            else
               !solid harmonic for l=0,m=0 equals 1.0_cfp
               eval_CGTO_shell(1,p) = norm*rad
            endif

         enddo !p

   end subroutine eval_cgto

  !> Only for debugging.
  subroutine compare_print_1el_ints(tag,integrals_1,int_index_1,integrals_2,int_index_2,n_integrals,column)
     implicit none
     character(len=4), intent(in) :: tag
     real(kind=cfp), allocatable :: integrals_1(:,:), integrals_2(:,:)
     integer, allocatable :: int_index_1(:,:), int_index_2(:,:)
     integer, intent(in) :: n_integrals, column

     integer :: i, j, ind_1, ind_2
     logical :: compared

        do i=1,n_integrals
           compared = .false.
           ind_1 = int_index_1(1,i)*(int_index_1(1,i)-1)/2+int_index_1(2,i)
           do j=1,n_integrals
              ind_2 = int_index_2(1,j)*(int_index_2(1,j)-1)/2+int_index_2(2,j)
              if (ind_1 .eq. ind_2) then
                 if (integrals_1(i,column) .ne. 0.0_cfp .and. integrals_2(j,column) .ne. 0.0_cfp) then
                    write(*,'(a,i0,3e25.15)') tag,i,integrals_1(i,column),integrals_2(j,column), &
                                                abs((integrals_1(i,column)-integrals_2(j,column))/integrals_1(i,column))
                 endif
                 compared = .true.
                 exit
              endif
           enddo !j
           if (.not.(compared)) write(*,'("No matching integral found for: ",i0)') i
        enddo !i

  end subroutine compare_print_1el_ints

  !> Only for debugging.
  subroutine compare_print_2el_ints(tag,integrals_1,int_index_1,integrals_2,int_index_2,n_integrals,column)
     implicit none
     character(len=4), intent(in) :: tag
     real(kind=cfp), allocatable :: integrals_1(:,:), integrals_2(:,:)
     integer, allocatable :: int_index_1(:,:), int_index_2(:,:)
     integer, intent(in) :: n_integrals, column

     integer :: i, j, ind_1, ind_2, p,q, r,s
     logical :: compared

        do i=1,n_integrals
           compared = .false.
           p = int_index_1(1,i)*(int_index_1(1,i)-1)/2+int_index_1(2,i)
           q = int_index_1(3,i)*(int_index_1(3,i)-1)/2+int_index_1(4,i)
           r = max(p,q)
           s = min(p,q)
           ind_1 = r*(r-1)/2+s
           do j=1,n_integrals
              p = int_index_2(1,j)*(int_index_2(1,j)-1)/2+int_index_2(2,j)
              q = int_index_2(3,j)*(int_index_2(3,j)-1)/2+int_index_2(4,j)
              r = max(p,q)
              s = min(p,q)
              ind_2 = r*(r-1)/2+s
              if (ind_1 .eq. ind_2) then
                 if (integrals_1(i,column) .ne. 0.0_cfp .and. integrals_2(j,column) .ne. 0.0_cfp) then
                    write(*,'(a,i0,3e25.15)') tag,i,integrals_1(i,column),integrals_2(j,column), &
                                                abs((integrals_1(i,column)-integrals_2(j,column))/integrals_1(i,column))
                 endif
                 compared = .true.
                 exit
              endif
           enddo !j
           if (.not.(compared)) write(*,'("No matching integral found for: ",i0)') i
        enddo !i

  end subroutine compare_print_2el_ints

end module gto_routines
