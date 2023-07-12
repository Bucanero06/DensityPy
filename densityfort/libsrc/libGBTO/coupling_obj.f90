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
module coupling_obj
   use wigner_cf
   use precisn
   use utils, only: xermsg
   use common_obj, only: darray_1d
   use const, only: stdout, mib, cache_line_size
   use special_functions, only: cfp_gamma_fun

   private

   public couplings_type

   !> \class <couplings>
   !> This object calculates the coupling coefficients: Clebsch-Gordan coefficients, Gaunt coefficients for the real and complex spherical harmonics and the coefficients \f$ G(l1m1|l2m2|m3) \f$ from the
   !> addition theorem for real solid harmonics.
   !> The Gaunt coefficients for complex spherical harmonics can be precalculated using the prec_cgaunt routine. The coefficients are then retrieved by calling the routine cgaunt which decides whether
   !> to calculate the required coefficient (if it has not been precalculated) or whether to retrieve it from the array cgaunt_cf if it has been precalculated. The use of the prec_cgaunt routine is 
   !> optional, i.e. it is not required to run it before any coefficients are required. All coefficients calculated by this object (except of Clebsch-Gordan coefficients) are proportional to the Gaunt
   !> coefficients for complex spherical harmonics. Therefore in applications requiring a large number of coupling coefficients it is recommended that these are precalculated using prec_cgaunt which will
   !> guarantee the best performance. Some additional coefficients (Pochhammer symbols) are required to calculate the coefficients \f$ G(l1m1|l2m2|m3) \f$. These additional coefficients can also be 
   !> precalculated using prec_G_cf. Note that this routine automatically precalculates also the Gaunt coefficients and therefore if prec_G_cf is called prec_cgaunt doesn't have to be called at all.
   !> \todo perform checks on the calculated values of the coefficients and force to zero those that are smaller than d1mach(4).
   type couplings_type
     !> This is set to .true. if some complex Gaunt coefficients have been precalculated. The precalculated coefficients are used in cgaunt. If a coefficient has not been precalculated, 
     !> it is directly calculated.
     logical, private :: cgaunt_precalculated = .false.
     !> Set to .true. if the Pochhammer symbols from the \f$ G(l1m1|l2m2|m3) \f$ coefficients have been precalculated. These are stored in the array g_cgaunt_cf.
     logical, private :: pochham_precalculated = .false.
     !> The array of precalculated Gaunt coefficients for complex spherical harmonics. The array is controlled by precalculate_cgaunt.
     type(darray_1d), allocatable, private :: cgaunt_cf(:)
     !> If pochham_precalculated .eq. .true. then this array contains the values \f$ 2\pi\frac{(1/2){l1+1}}{(1/2)_{l2+1}(1/2)_{l1-l2+1}}, l1=0,\dots,L; l2=0,\dots,L \f$ for L=LG. 
     !> These values are used in G_real_cf to obtain \f$ G(l1m1|l2m2|m3) \f$ (see G_real_cf). This array is controlled by prec_G_cf.
     real(kind=cfp), allocatable, private :: pochham(:,:)
     !> Maximum L for which the Pochhammer symbols have been precalculated.
     integer, private :: LG = 0
     !> L from equation (18) of Pinchon (2007), i.e. the largest L for which the Gaunt coefficients have been precalculated.
     integer, private :: L = 0
     !> U(L) from equation (18) of Pinchon (2007), i.e. dimension of the array cgaunt_cf.
     integer, private :: U = 0
   contains
     !> This can be used to precalculate the Gaunt coefficients for complex spherical harmonics.
     procedure :: prec_cgaunt => precalculate_cgaunt
     !> Implements the indexing function of Pinchon, Hoggan: New index Functions for storing Gaunt coefficients. Int. J. Quant. Chem., 107, 2186-2196 (2007).
     !> It is assumed that the input parameters have been ordered as follows:
     !> l1 >= 0, [l1/2] <= l2 <= l1, l1-l2 <= l3 <= l2, l1+l2+l3 = 2*n, 0 <= m3 <= l3.
     !> This function has been made private since it doesn't perform checking of the order of the arguments as given above and in any case the array cgaunt_cf is private and there is no reason why the user
     !> should have a direct access to it. If the order of the arguments is not the one given above then the index function will return incorrect values. The correct order of the arguments prior call 
     !> to this function is always ensured in this object. The routine performs no checking for performance reasons since it is expected that it will be called many times.
     procedure, private, nopass :: cgaunt_index_f2
     !> Calculate directly the C-G coefficient for integer quantum numbers. The function can be easilly generalized to accept half-integer arguments, but we don't need it. Note that this calls my wigner_3j
     !> routine which has some problems with phases for high L values (L > 14). It should be modified to use wigner3j instead or wigner_3j must be corrected.
     procedure, nopass :: cg => calc_cg_cf
     !> Calculate the Gaunt coefficient for complex spherical harmonics. The coefficient is either calculated directly, or, if the coefficient has been precalculated, it is fetched from the memory.
     !> The indexing and fetching from memory vs. direct evaluation has been tested for all L up to L = 40 but it is completely general. The routine assumes l1 .ge. 0, l2 .ge. 0, l3 .ge. 0 and doesn't 
     !> check this. The coefficients are defined in the following way:
     !> \f[
     !>  G(l_{1},l_{2},l_{3},m_{1},m_{2},m_{3}) = \int_{0}^{2\pi}d\phi \int_{0}^{\pi}\sin\theta d\theta {\overline{Y}}_{l_{3}m_{3}}(\theta,\phi)Y_{l_{2}m_{2}}(\theta,\phi)Y_{l_{1}m_{1}}(\theta,\phi).
     !> \f]
     !> Note that this definition is equivalent to that of Homeier and Steinborn, see eq (A2) of their paper, but differs from that of Pinchon and Hoggan, see eq (6) of their paper.
     !> The coefficients are calculated using the explicit expression:
     !> \f[
     !>  G(l_{1},l_{2},l_{3},m_{1},m_{2},m_{3}) = {\sqrt{\frac{(2l_{1}+1)(2l_{2}+1)(2l_{3}+1)}{4\pi}}}W(l_{1},l_{2},l_{3};0,0,0)W(l_{1},l_{2},l_{3};m_{1},m_{2},-m_{3}),
     !> \f]
     !> where \f$ W(l_{3},l_{2},l_{1};m_{3},m_{2},m_{1}) \f$ is the standard Wigner 3j symbol.
     procedure :: cgaunt => calc_complex_gaunt_cf
     !> Calculate the Gaunt coefficient for the real spherical harmonics. These are accurate to 14 decimal digits. These coefficients are always proportional to a single Gaunt coefficient or zero.
     !> The method of calculating the real Gaunt coefficients follows the Theochem paper of Homeier and Steinborn (HS). However, they defined real spherical harmonics without the (-1)**m factor. In this code we
     !> include the (-1)**m factor in the definition of the real spherical harmonics (that is the standard convention). There is no need to modify the method of HS for our definition of
     !> the real spherical harmonics since the transformation of our real Gaunt coefficients to those of HS is achieved multiplying the HS real Gaunt coeffs. by (-1)**(m1+m2+m3). However,
     !> the sum m1+m2+m3 must be even for a non-zero real Gaunt coefficient so in fact the coefficients four our real spherical harmonics and those of HS do not differ.
     !> \todo Get rid of the calls to ulmmu and replace it by the actual numbers. This should speed up the rgaunt a bit.
     procedure :: rgaunt => calc_real_gaunt_cf
     !> Returns the bounds [l_min;l_max] on l_1 for the real Gaunt coefficient \f$ RG(l_{1},l_{2},l_{3},m_{1},m_{2},m_{3}) \f$.
     procedure, nopass :: bounds_rg => get_l_bounds_rg
     !> This routine precalculates the values needed to evaluate the the coefficients \f$ G(l1m1|l2m2|m3) \f$ for \f$ l1=0,\dots,L; l2 \le l1 \f$. This includes precalculating the Gaunt coefficients and the
     !> Pochhammer symbols.
     procedure :: prec_G_cf => precalculate_G_coeff
     !> Calculates the coefficient \f$ G(l1m1|l2m2|m3) = 2\pi\frac{(1/2){l1+1}}{(1/2)_{l2+1}(1/2)_{l1-l2+1}}<l1m1|l2m2|l1-l2 m3>_{R} \f$, where \f$ <l1m1|l2m2|l1-l2 m3>_{R} \f$ is the Gaunt coefficient 
     !> for real spherical harmonics. The coefficients \f$ G(l1m1|l2m2|m3) \f$ arise in the addition theorem for real solid harmonics:
     !> \f[
     !>    (-1)^{m1}r^{l1}_{1+2} X_{l1,m1}(r_{1+2}) = \sum_{l2=0}^{l1}\sum_{m2=-l2}^{l2}\sum_{m3=-(l1-l2)}^{l1-l2} G(l1m1|l2m2|m3) \times
     !>                                               (-1)^{m2}r^{l2}_{1} X_{l2,m2}(r_{1}) (-1)^{m3}r^{l1-l2}_{2} X_{l1-l2,m3}(r_{2}).
     !> \f]
     !> The coefficient is either calculated directly (if it has not been precalculated), or it is fetched from cgaunt_cf and pochham if prec_G_cf has been called. Since the G(l1m1|l2m2|m3) coefficient is 
     !> directly proportional to a single Gaunt coefficient for real spherical harmonics (which is in turn proportional to a single Gaunt coefficient) all that we need to calculate G efficiently is to
     !> precalculate the usual Gaunt coefficients and the Pochhammer symbols. Therefore the method of obtaining the coefficients here is optimal in the sense of minimizing FLOPS. An alternative method of
     !> storing and fetching the coefficients \f$ G(l1m1|l2m2|m3) \f$ is to note that \f$ G'(l1m1|l2m2|m3) = 2\pi\frac{1}{(1/2){l1+1}(1/2)_{l2+1}(1/2)_{l1-l2+1}}<l1m1|l2m2|l1-l2 m3> \f$ have the same 
     !> symmetry properties as the Gaunt coefficients. Therefore we can store and index \f$ G'(l1m1|l2m2|m3) \f$ in the same way as \f$ <l1m1|l2m2|l3 m3> \f$. The coefficient \f$ G(l1m1|l2m2|m3) \f$ is 
     !> then given by \f$ G(l1m1|l2m2|m3) = (1/2){l1+1}^{2}G'(l1m1|l2m2|m3) \f$. For a given maximum L there will only be L values \f$ (1/2){l+1}^{2}, l=0,\dots,L \f$ needed to obtain G from G' so this 
     !> scheme avoids calculation of \f$ L^2 \f$ Pochhammer symbols \f$\frac{1}{(1/2){l1+1}(1/2)_{l2+1}(1/2)_{l1-l2+1}}\f$. However, performance-wise this scheme will probably be indistinguishable from
     !> the one employed at the moment.
     procedure :: G_real_cf => G_coeff
     !> Writes the precalculated Gaunt coefficients to the specified file.
     procedure :: write_prec_cgaunt
     !> Reads the precalculated Gaunt coefficients from the specified file.
     procedure :: read_prec_cgaunt
   end type couplings_type

   character(len=*), parameter :: header = "CGAUNT COEFFICIENTS"

contains

   subroutine write_prec_cgaunt(this,path)
     use mpi_mod
     implicit none
     class(couplings_type) :: this
     character(len=*), intent(in) :: path

     integer :: err, lu, i, j

        if (.not. this % cgaunt_precalculated) then
            call xermsg ('coupling_obj', 'write_prec_cgaunt', &
                         'Attempt to write coefficients but these were not precalculated.', 1, 1)
        end if

        call mpi_mod_barrier(err)

        !Only master writes
        if (myrank .eq. master) then
           open(file=path,newunit=lu,access='stream',form='unformatted',status='replace',iostat=err)
           if (err .ne. 0) then
              call xermsg('coupling_obj','write_prec_cgaunt','Error opening the file for output.',err,1)
           endif

           write(lu) header
           write(lu) this%L, this%U

           do i=1,this%U
              write(lu) this%cgaunt_cf(i)%d1
              write(lu) this%cgaunt_cf(i)%a
           enddo !i
        endif

        call mpi_mod_barrier(err)

   end subroutine write_prec_cgaunt

   subroutine read_prec_cgaunt(this,path)
     use mpi_mod
     use omp_lib
     implicit none
     class(couplings_type) :: this
     character(len=*), intent(in) :: path

     integer :: err, lu, i, j
     character(len=len(header)) :: inp_header


        if (.not.(omp_in_parallel())) then
           write(stdout,'("Precalculated Gaunt coefficients will be read-in from the file: ",a)') path
        endif

        call mpi_mod_barrier(err)

        if (this%cgaunt_precalculated) then
           if (allocated(this%cgaunt_cf)) deallocate(this%cgaunt_cf)
        endif

        this%cgaunt_precalculated = .false.

        !Only master reads
        if (myrank .eq. master) then
           open(file=path,newunit=lu,access='stream',form='unformatted',status='old',iostat=err)
           if (err .ne. 0) call xermsg('coupling_obj','read_prec_cgaunt','Error opening the file for input.',err,1)

           read(lu) inp_header

           if (inp_header .ne. header) then
              print *,inp_header,header
              call xermsg('coupling_obj','read_prec_cgaunt','Header stored on the file does not match with the expected one.',1,1)
           endif
           
           read(lu) this%L, this%U

           allocate(this%cgaunt_cf(this%U),stat=err)
           if (err .ne. 0) call xermsg('coupling_obj','read_prec_cgaunt','Memory allocation failed.',err,1)

           do i=1,this%U
              read(lu) this%cgaunt_cf(i)%d1
              read(lu) this%cgaunt_cf(i)%a
           enddo !i
        endif

        call mpi_mod_barrier(err)

   end subroutine read_prec_cgaunt
   
   subroutine precalculate_cgaunt(this,l1m)
     use omp_lib
     implicit none
     class(couplings_type) :: this
     !> The limiting value for which the complex Gaunt coefficients will be precalculated. I.e. the last one will be:  <l1m l1m|l1m l1m|l1m l1m>.
     integer, intent(in) :: l1m

     type integer_list
        integer(kind=1), allocatable :: a(:) !this array should be logical however intel implements logical as integer of default size -> waste of memory!!!
     end type integer_list
     
     integer :: l1, l2, l3, m1, m2, m3, m3p, ind, err, b, tot, ld, m2d, m3pd, bd, indd, m1d, l_min, l_max, l1d,l2d,l3d
     type(integer_list), allocatable :: missing_cf(:)
     real(kind=cfp) :: start_t, end_t, gaunt_cf(0:max(1,2*l1m))

        if (l1m < 0) then
           call xermsg('coupling_obj','precalculate_cgaunt','The input L value must be .ge. 0.',1,1)
        endif

        if (.not.(this%cgaunt_precalculated) .or. l1m > this%L) then

           this%cgaunt_precalculated = .false.

           if (allocated(this%cgaunt_cf)) deallocate(this%cgaunt_cf)

           this%L = l1m
           if (mod(l1m,2) .eq. 0) then
              this%U = (l1m+2)*(l1m+4)*(3*l1m*l1m+14*l1m+24)/192
           else
              this%U = (l1m+1)*(l1m+3)*(l1m+5)*(3*l1m+5)/192
           endif

           if (.not.(omp_in_parallel())) then
              write(stdout,'(/,"--------->","couplings_type: precalculating the Gaunt coefficients for L_max = ",i0)') l1m
              write(stdout,'(/,"Size of the buffer cgaunt_prec is: ",i0)') this%U
           endif

           allocate(this%cgaunt_cf(this%U),missing_cf(this%U),stat=err)
           if (err .ne. 0) call xermsg('coupling_obj','precalculate_cgaunt','Memory allocation failed.',err,1)

           if (.not.(omp_in_parallel())) then
              write(stdout,'("Precalculating Gaunt coefficients for complex spherical harmonics...")')
           endif

           call cpu_time(start_t)
           !First allocate space for all coefficients and the counter of coefficients to calculate:
           tot = 0
           do l1=0,this%L
              do l2=ceiling(l1/2.0),l1
                 do l3=l1-l2,l2
                    do m3=0,l3

                       if (mod(l1+l2+l3,2) .eq. 0) then

                          ind = this%cgaunt_index_f2(l1,l2,l3,m3)

                          !m2, m3p are just dummy variables in this call to m2_limit since we're not interested in the optimal m2 at this stage, we just need b
                          m2 = 0
                          m3p = m3
                          b = m2_limit(l1,l2,l3,m2,m3p)

                          if (this%cgaunt_cf(ind)%d1 .eq. 0) then
                             this%cgaunt_cf(ind)%d1 = b + l2 + 1
                             allocate(this % cgaunt_cf(ind) % a(this % cgaunt_cf(ind) % d1), &
                                      missing_cf(ind) % a(this % cgaunt_cf(ind) % d1), stat = err)
                             if (err /= 0) then
                                call xermsg ('coupling_obj', 'precalculate_cgaunt', &
                                             'Memory allocation of cgaunt_prec(ind)%a failed.',err,1)
                             end if
                             this%cgaunt_cf(ind)%a(1:this%cgaunt_cf(ind)%d1) = 0.0_cfp
                             missing_cf(ind)%a(:) = 1 !1 = coefficient is missing
                          else
                             call xermsg ('coupling_obj', 'precalculate_cgaunt', &
                                          'Error in the index function cgaunt_index_f2: &
                                          &two different sets (l1,l2,l3,m3) have the same index.', 2, 1)
                          endif

                          tot = tot + this%cgaunt_cf(ind)%d1

                       endif

                    enddo !m3
                 enddo !l3
              enddo !l2
           enddo !l1

           !Now calculate all coefficients utilizing the whole string of Gaunt coefficients calculated for each combination of l1,l2,l3,m3:
           do l1=0,this%L 
              do l2=ceiling(l1/2.0),l1
                 do l3=l1-l2,l2
                    do m3=0,l3
   
                       if (mod(l1+l2+l3,2) .eq. 0) then
   
                          ind = this%cgaunt_index_f2(l1,l2,l3,m3)

                          !m2, m3p are just dummy variables in this call to m2_limit since we're not interested in the optimal m2 at this stage, we just need b
                          m2 = 0
                          m3p = m3
                          b = m2_limit(l1,l2,l3,m2,m3p)

                          do m2=-l2,b
                             !We store in this%cgaunt_cf the following quantity: (-1)**m3 * integral Y_{l1,m1}*Y_{l2,m2}*Y_{l3,m3}
                             !This is (up to the phase (-1)**m3) the Gaunt coefficient of Pinchon and Hoggan:
                             !(-1)**m3 * <l1m1|l2m2|l3m3>_Pinchon = <l1m1|l2m2|l3-m3>_Our
                             !The reasons why we don't get rid of the phase (-1)**m3 here are explained in calc_complex_gaunt_cf.
                             m1 = -m2-m3
                             if (missing_cf(ind)%a(l2+m2+1) .eq. 1) then
                                !This coefficient has not been calculated: calculate the whole string of coefficients for all allowed values of l1
                                call cgaunt_string(gaunt_cf,l_min,l_max,l2,l3,m1,m2,-m3)
                                !Save the coefficient:
                                this%cgaunt_cf(ind)%a(l2+m2+1) = gaunt_cf(l1) != this%cgaunt(l1,l2,l3,m1,m2,-m3)
                                missing_cf(ind)%a(l2+m2+1) = 0
                                !Store all the other calculated Gaunt coefficients: repeat the same procedure to calculate indices of the coefficients.
                                do ld=max(l_min,l2),min(l_max,this%L)
                                   !We need to make sure that the values of l1 (=ld) satisfy the equation (19) of the Hoggan and Pinchon paper:
                                   if (mod(ld+l2+l3,2) .ne. 0) cycle
                                   if (l2 .ge. ceiling(ld/2.0) .and. l3 .ge. ld-l2 .and. ld .ge. m1) then
                                      indd = this%cgaunt_index_f2(ld,l2,l3,m3)
                                      m2d = 0
                                      m3pd = m3
                                      bd = m2_limit(ld,l2,l3,m2d,m3pd)
                                      if (bd < m2) cycle !m2d=-l2,...,bd so we have to make sure it includes the m2 for which the Gaunt coefficient has been precalculated.
                                      if (missing_cf(indd)%a(l2+m2+1) .eq. 1) then
                                         this%cgaunt_cf(indd)%a(l2+m2+1) = gaunt_cf(ld)
                                         missing_cf(indd)%a(l2+m2+1) = 0
                                      endif
                                   endif
                                enddo !ld
                             endif
                          enddo !m2
      
                       endif
   
                    enddo !m3
                 enddo !l3
              enddo !l2
           enddo !l1
           call cpu_time(end_t)

           deallocate(missing_cf)

           this%cgaunt_precalculated = .true.

           if (.not.(omp_in_parallel())) then
              write(stdout,'("Total number of Gaunt coefficients evaluated and stored: ",i0)') tot
              write(stdout,'("Memory usage for storing the Gaunt coefficients: ",f8.3," [Mib]")') tot*8/real(mib,kind=cfp)
              write(stdout,'("Precalculating the coefficients took [s]: ",f8.3)') end_t-start_t
              write(stdout,'(/,"<---------","...finished")')
           endif

        endif

   end subroutine precalculate_cgaunt

   !> This routine determines the limit on m2 for the Gaunt coefficients which have been stored and also the smallest possible m2 compatible with the symmetries of the given Gaunt coefficient.
   !> We don't need to know the value of m1 since that is uniquely determined from: m1 = -m2-m3. This routine is the function b(l1,l2,l3,m3) from the paper of Pinchon (2007).
   function m2_limit(l1,l2,l3,m2,m3)
     implicit none
     integer, intent(in) :: l1,l2,l3
     integer, intent(inout) :: m2,m3
     integer :: m2_limit

     integer :: t

        if (l1 .eq. l2) then !m1 and m2 can be exchanged; note that we take into account the fact that <l1m1|l2m2|l3m3> = <l1-m1|l2-m2|l3-m3>
           m2_limit = -ceiling(m3/2.0)
           m2 = min(m2,-m2,-m2-m3,m2+m3) ! = min(m2,-m2,m1,-m1)
        elseif (m3 .eq. 0) then !m3 == 0, so we take into account the fact that <l1m1|l2m2|l3m3> = <l1-m1|l2-m2|l3-m3>
           m2_limit = 0
           m2 = min(m2,-m2)
        elseif (l2 .eq. l3) then !m2 and m3 can be exchanged
           m2_limit = min(l1-m3,m3)
           if (m3 < m2) then
              t = m2
              m2 = m3
              m3 = t
           endif
        else
           m2_limit = min(l1-m3,l2)
        endif

   end function m2_limit

   elemental function cgaunt_index_f2(l1,l2,l3,m3)
     implicit none
     integer, intent(in) :: l1,l2,l3,m3
     integer :: cgaunt_index_f2

        !we want to compute with integers only so we've multiplied the equation (21) by 384
        if (mod(l1,2) .eq. 0) then
           cgaunt_index_f2 = l1*l1*l1*l1*6 + l1*l1*l1*64 - l1*l1*24 -l1*64 - l2*96*(l1+1)*(l1-l2-2) + l3*l3*96 + 384*(m3 + 1)
        else
           cgaunt_index_f2 = l1*l1*l1*l1*6 + l1*l1*l1*64 + l1*l1*36 + l1*128 + 438 - l1*l2*96*(l1-l2-1) + l3*l3*96 + 384*m3
        endif

        cgaunt_index_f2 = cgaunt_index_f2/384

   end function cgaunt_index_f2

   function calc_cg_cf(j1,m1,j2,m2,J,M)
     implicit none
     real(kind=cfp) :: calc_cg_cf
     integer, intent(in) :: j1, j2, J, m1, m2, M !assumed to be positive WHOLE integers

        !\todo: some error checking needed?
   
        calc_cg_cf = C_G_cf(j1,m1,j2,m2,J,M)

   end function calc_cg_cf

   function calc_complex_gaunt_cf(this,l1,l2,l3,m1,m2,m3)
     use phys_const, only: inv_fourpi
     implicit none
     class(couplings_type), intent(in) :: this
     real(kind=cfp) :: calc_complex_gaunt_cf
     integer, intent(in) :: l1,l2,l3,m1,m2,m3

     real(kind=cfp) :: fac, w3j1, w3j2, w(l2+l3+1)
     integer :: l_min, l_max, ind, m1p, m2p, m3p, b

        calc_complex_gaunt_cf = 0.0_cfp
   
        !Sufficient conditions for the Gaunt coefficient to be zero:
        if (abs(m1) > l1 .or. abs(m2) > l2 .or. abs(m3) > l3) return
        if (l1 > l2 + l3 .or. l2 > l1 + l3 .or. l3 > l1 + l2) return
        if (mod(l1+l2+l3,2) .ne. 0) return
        if (m1 + m2 .ne. m3) return
   
        !now calculate the two Wigner 3j symbols or get the Gaunt coefficient from the buffer if we have them precalculated:
        if (this%cgaunt_precalculated .and. max(l1,l2,l3) .le. this%L) then
           !Here we map our Gaunt coefficient onto the stored Gaunt coefficient of Pinchon and Hoggan which we fetch from the memory.
           !We stored the Gaunt coefficients with a given order of l1,l2,l3 and m3 values (see documentation for cgaunt_index_f2) so we need to permute the triplets to that order.

           if (l1 .ge. l2 .and. l1 .ge. l3 .and. l2 .ge. l3) then     !l1 l2 l3
              if (l2 .ge. ceiling(l1/2.0)) then
                 m2p = m2; m3p = -m3; !Mapping between the Gaunt coefficient that we want (m2,m3) and the Gaunt coefficient of Pinchon and Hoggan (m2p,m3p). m1p = -m2p-m3p so we don't need to calculate it.
                 !The Gaunt coefficients stored all have m3p >= 0 so we may need to change signs of all m-values if m3p < 0 (<l1m1|l2m2|l3m3> = <l1-m1|l2-m2|l3-m3>).
                 if (m3p < 0) then
                    m3p = -m3p
                    m2p = -m2p
                 endif
                 b = m2_limit(l1,l2,l3,m2p,m3p)
                 if (m2p .le. b .and. abs(m2p) .le. l2) then
                    ind = this%cgaunt_index_f2(l1,l2,l3,m3p)
                    !Retrieve: '(-1)**m3 * integral Y_{l1,-m2p-m3p}*Y_{l2,m2p}*Y_{l3,m3p}' and multiply by the phase (-1)**m3 
                    !which originates in the mapping between our Gaunt coefficient and the Gaunt cf. of Pinchon and Hoggan.
                    !Note that (-1)**(m3+m3) = 1 and so we actually don't multiply any phase. This property of the mapping
                    !allows us to store in this%cgaunt_cf the values '(-1)**m3 * integral Y_{l1,-m2p-m3p}*Y_{l2,m2}*Y_{l3,m3}' instead
                    !of 'integral Y_{l1,-m2p-m3p}*Y_{l2,m2}*Y_{l3,m3}'. This is the same in all other cases.
                    calc_complex_gaunt_cf = this%cgaunt_cf(ind)%a(l2+m2p+1)
                 endif
              endif
           elseif (l2 .ge. l1 .and. l2 .ge. l3 .and. l1 .ge. l3) then !l2 l1 l3
              if (l1 .ge. ceiling(l2/2.0)) then
                 m1p = m1; m3p = -m3
                 if (m3p < 0) then
                    m1p = -m1p
                    m3p = -m3p
                 endif
                 b = m2_limit(l2,l1,l3,m1p,m3p)
                 if (m1p .le. b .and. abs(m1p) .le. l1) then
                    ind = this%cgaunt_index_f2(l2,l1,l3,m3p)
                    calc_complex_gaunt_cf = this%cgaunt_cf(ind)%a(l1+m1p+1) 
                 endif
              endif
           elseif (l3 .ge. l1 .and. l3 .ge. l2 .and. l1 .ge. l2) then !l3 l1 l2
              if (l1 .ge. ceiling(l3/2.0)) then
                 m2p = m2; m1p = m1
                 if (m2p < 0) then
                    m2p = -m2p
                    m1p = -m1p
                 endif
                 b = m2_limit(l3,l1,l2,m1p,m2p)
                 if (m1p .le. b .and. abs(m1p) .le. l1) then
                    ind = this%cgaunt_index_f2(l3,l1,l2,m2p)
                    calc_complex_gaunt_cf = this%cgaunt_cf(ind)%a(l1+m1p+1) 
                 endif
              endif
           elseif (l2 .ge. l3 .and. l2 .ge. l1 .and. l3 .ge. l1) then !l2 l3 l1
              if (l3 .ge. ceiling(l2/2.0)) then
                 m1p = m1; m3p = -m3
                 if (m1p < 0) then
                    m3p = -m3p
                    m1p = -m1p
                 endif
                 b = m2_limit(l2,l3,l1,m3p,m1p)
                 if (m3p .le. b .and. abs(m3p) .le. l3) then
                    ind = this%cgaunt_index_f2(l2,l3,l1,m1p)
                    calc_complex_gaunt_cf = this%cgaunt_cf(ind)%a(l3+m3p+1) 
                 endif
              endif
           elseif (l3 .ge. l2 .and. l3 .ge. l1 .and. l2 .ge. l1) then !l3 l2 l1
              if (l2 .ge. ceiling(l3/2.0)) then
                 m1p = m1; m2p = m2
                 if (m1p < 0) then
                    m1p = -m1p
                    m2p = -m2p
                 endif
                 b = m2_limit(l3,l2,l1,m2p,m1p)
                 if (m2p .le. b .and. abs(m2p) .le. l2) then
                    ind = this%cgaunt_index_f2(l3,l2,l1,m1p)
                    calc_complex_gaunt_cf = this%cgaunt_cf(ind)%a(l2+m2p+1) 
                 endif
              endif
           elseif (l1 .ge. l3 .and. l1 .ge. l2 .and. l3 .ge. l2) then !l1 l3 l2
              if (l3 .ge. ceiling(l1/2.0)) then
                 m2p = m2; m3p = -m3
                 if (m2p < 0) then
                    m2p = -m2p
                    m3p = -m3p
                 endif
                 b = m2_limit(l1,l3,l2,m3p,m2p)
                 if (m3p .le. b .and. abs(m3p) .le. l3) then
                    ind = this%cgaunt_index_f2(l1,l3,l2,m2p)
                    calc_complex_gaunt_cf = this%cgaunt_cf(ind)%a(l3+m3p+1) 
                 endif
              endif
           endif

        else

           call wigner3j(w,l_min,l_max,l2,l3,m1,m2,-m3)
           w3j1 = w(l1-l_min+1)
       
           if (w3j1 .eq. 0.0_cfp) return
     
           call wigner3j(w,l_min,l_max,l2,l3,0,0,0)
           w3j2 = w(l1-l_min+1)
      
           if (w3j2 .eq. 0.0_cfp) return

           fac = sqrt((2*l1+1.0_cfp)*(2*l2+1.0_cfp)*(2*l3+1.0_cfp)*inv_fourpi)
           calc_complex_gaunt_cf = fac*w3j1*w3j2

        endif

   end function calc_complex_gaunt_cf

   !> Calculates a string of (complex) Gaunt coefficients for the given values of l2,l3,m1,m2,m3. The returned values are l_min,l_max which are the range of allowed l1 values and the array gaunt_cf
   !> containing the Gaunt coefficient for each allowed l1, e.g. gaunt_cf(l1) = <l1,m1|l2,m2|l3,m3>.
   subroutine cgaunt_string(gaunt_cf,l_min,l_max,l2,l3,m1,m2,m3)
     use phys_const, only: inv_fourpi
     implicit none
     real(kind=cfp), intent(out) :: gaunt_cf(0:)
     integer, intent(out) :: l_min,l_max
     integer, intent(in) :: l2,l3,m1,m2,m3

     integer :: l1, l_min_l,l_max_l,l_min_0,l_max_0
     real(kind=cfp) :: fac, wl(l2+l3+1), w0(l2+l3+1)

        call wigner3j(wl,l_min_l,l_max_l,l2,l3,m1,m2,-m3)
        call wigner3j(w0,l_min_0,l_max_0,l2,l3,0,0,0)

        fac = sqrt((2*l2+1.0_cfp)*(2*l3+1.0_cfp)*inv_fourpi)

        l_min = max(l_min_l,l_min_0); l_max = min(l_max_l,l_max_0)
        do l1=l_min,l_max
           gaunt_cf(l1) = wl(l1-l_min_l+1)*w0(l1-l_min_0+1)*sqrt(2*l1+1.0_cfp)*fac
        enddo !l1

   end subroutine cgaunt_string

   subroutine get_l_bounds_rg(l2,l3,m2,m3,l_min,l_max)
     implicit none
     integer, intent(in) :: l2,l3,m2,m3
     integer, intent(out) :: l_min, l_max

     integer :: k

        !determine l_max and l_min (i.e. bounds on l1). Homeier, Steinborn: J.Mol. Struc. 368 (1996), 31-37.
        l_max = l2 + l3
        k = max(abs(l2-l3),min(abs(m2+m3),abs(m2-m3)))
   
        if (mod(k+l_max,2) == 0) then !k+l_max is even
           l_min = k
        else !k+l_max is odd
           l_min = k + 1
        endif

   end subroutine get_l_bounds_rg

   function calc_real_gaunt_cf(this,l1,l2,l3,m1,m2,m3)
     implicit none
     class(couplings_type), intent(in) :: this
     real(kind=cfp) :: calc_real_gaunt_cf
     integer, intent(in) :: l1,l2,l3,m1,m2,m3

     integer :: l_max, l_min, p, k

        calc_real_gaunt_cf = 0.0_cfp
   
        !determine l_max and l_min (i.e. bounds on l1). Homeier, Steinborn: J.Mol. Struc. 368 (1996), 31-37.
        l_max = l2 + l3
        k = max(abs(l2-l3),min(abs(m2+m3),abs(m2-m3)))

        if (mod(k+l_max,2) == 0) then !k+l_max is even
           l_min = k
        else !k+l_max is odd
           l_min = k + 1
        endif

   
        !sufficient conditions for the real gaunt coefficient to be zero:
        if (.not.(l1 .ge. l_min .and. l1 .le. l_max)) return
        !todo the condition on the m-values is equivalent to the requirement m1+m2+m3 = even
        if (.not.(m1 == m2+m3 .or. m1 == m2-m3 .or. m1 == -m2+m3 .or. m1 == -m2-m3)) return
        if (mod(l1+l2+l3,2) .ne. 0) return
  
        p = -2
        if (mod(m3,2) .eq. 0) p = 2 !p = 2*(-1)**m3
   
        !now determine the R-Gaunt coefficient according to cases A to C:
        !Case A:
        if (m1 .ne. 0 .and. m2 .ne. 0 .and. m3 .ne. 0) then
           if (m3 .eq. m1+m2 .or. m3 .eq. -(m1+m2)) then
              calc_real_gaunt_cf = p*this%cgaunt(l1,l2,l3,m1,m2,m1+m2)*real(conjg(ulmmu(m1+m2,m3))*ulmmu(m2,m2)*ulmmu(m1,m1))
              return
           else
              calc_real_gaunt_cf = p*this%cgaunt(l1,l2,l3,-m1,m2,-m1+m2)*real(conjg(ulmmu(-m1+m2,m3))*ulmmu(m2,m2)*ulmmu(-m1,m1))
              return
           endif
        endif
   
        !Case B: one mi vanishes
        if (m1 .eq. 0 .and. m2 .ne. 0 .and. m3 .ne. 0) then
           calc_real_gaunt_cf = p*this%cgaunt(l1,l2,l3,m1,m2,m2)*real(conjg(ulmmu(m2,m3))*ulmmu(m2,m2))
           return
        endif
        if (m2 .eq. 0 .and. m1 .ne. 0 .and. m3 .ne. 0) then
           calc_real_gaunt_cf = p*this%cgaunt(l1,l2,l3,m1,m2,m1)*real(conjg(ulmmu(m1,m3))*ulmmu(m1,m1))
           return
        endif
        if (m3 .eq. 0 .and. m1 .ne. 0 .and. m2 .ne. 0) then
           calc_real_gaunt_cf = p*this%cgaunt(l1,l2,l3,m1,-m1,m3)*real(conjg(ulmmu(m1,m1)*ulmmu(-m1,m2)))
           return
        endif
   
   
        !Case C: two or more mi vanish
        if (m1 .eq. 0 .and. m2 .eq. 0 .and. m3 .eq. 0) then
           calc_real_gaunt_cf = this%cgaunt(l1,l2,l3,0,0,0)
           return
        endif
   
   !    General formula that works for all cases. Only used for debugging.
   !     do i = -l1,l1
   !        do j = -l2,l2
   !           calc_real_gaunt_cf = calc_real_gaunt_cf + conjg(ulmmu(i+j,m3))*ulmmu(i,m1)*ulmmu(j,m2)*this%cgaunt(l1,l2,l3,i,j,i+j)
   !        enddo
   !     enddo
   !     return

   end function calc_real_gaunt_cf

   !> This is a helper function for calc_real_gaunt_cf which implements the function U_{l m}^{mu} of Homeier and Steinborn. See equation (15) in their paper.
   elemental function ulmmu(m,mu)
     use phys_const, only: imu, roneh
     implicit none
     integer, intent(in) :: m, mu
     complex(kind=cfp) :: ulmmu

        ulmmu = 0.0_cfp
   
        if (abs(m) .ne. abs(mu)) return
   
        if (mu .eq. 0) then 
           ulmmu = 1.0_cfp
           return
        endif
   
        if (mu > 0) then
           if (m > 0) then
              ulmmu = roneh
              return
           else !m < 0
              ulmmu = roneh
              if (mod(m,2) .ne. 0) ulmmu = -ulmmu !ulmmu = roneh*(-1)**m
              return
           endif
        else !mu < 0
           if (m > 0) then
              ulmmu = -imu*roneh
              return
           else !m < 0
              ulmmu = roneh*imu
              if (mod(m,2) .ne. 0) ulmmu = -ulmmu !ulmmu = roneh*imu*(-1)**m
              return
           endif
        endif
     
   end function ulmmu

   subroutine precalculate_G_coeff(this,l1m)
     use phys_const, only: twopi
     use omp_lib
     implicit none
     class(couplings_type) :: this
     integer, intent(in) :: l1m
     
     integer :: l1, l2, err
     real(kind=cfp) :: gamma_half

        if (l1m < 0) then
           call xermsg('coupling_obj','precalculate_G_coeff','The input L value must be .ge. 0.',1,1)
        endif

        if (.not.(this%pochham_precalculated) .or. l1m > this%LG) then

           if (.not.(omp_in_parallel())) then
              write(stdout,'(/,"--------->","couplings_type: precalculating the data needed for the &
                            &addition theorem of real solid harmonics for L_max = ",i0)') l1m
           endif

           this%pochham_precalculated = .false.

           this%LG = l1m

           !First precalculate the Gaunt coefficients for complex spherical harmonics which will be used to get the Gaunt coefficients for the real spherical harmonics.
           call this%prec_cgaunt(l1m)

           if (allocated(this%pochham)) deallocate(this%pochham)

           allocate(this%pochham(l1m+1,l1m+1),stat=err)
           if (err .ne. 0) call xermsg('coupling_obj','precalculate_G_coeff','Memory allocation failed.',err,1)

           gamma_half = cfp_gamma_fun(0.5_cfp)

           this%pochham(:,:) = 0.0_cfp
           do l1=0,l1m
              do l2=0,l1
                 this % pochham(l1+1,l2+1) = twopi * cfp_gamma_fun(0.5_cfp+l1+1) * gamma_half &
                                            / (cfp_gamma_fun(0.5_cfp+l2+1) * cfp_gamma_fun(0.5_cfp+(l1-l2)+1))
              enddo
           enddo

           this%pochham_precalculated = .true.

           if (.not.(omp_in_parallel())) then
              write(stdout,'(/,"<---------","...finished")')
           endif

        endif

   end subroutine precalculate_G_coeff

   function G_coeff(this,l1,l2,m1,m2,m3)
     use phys_const, only: twopi
     implicit none
     class(couplings_type), intent(in) :: this
     real(kind=cfp) :: G_coeff
     integer, intent(in) :: l1,l2,m1,m2,m3

     integer :: l3
     real(kind=cfp) :: gamma_half

        G_coeff = 0.0_cfp

        l3 = l1-l2

        !This will always be determined from the precalculated buffer if prec_G_cf has been called first. Otherwise it may or may not be taken from the buffer depending on whether the Gaunt coefficients
        !have been precalculated using a separate call to prec_cgaunt.
        G_coeff = this%rgaunt(l1,l2,l3,m1,m2,m3)

        !Calculate the Pochammer symbols directly or get them from the buffer if we have them precalculated:
        if (this%pochham_precalculated .and. max(l1,l2) .le. this%LG) then
           G_coeff = this%pochham(l1+1,l2+1)*G_coeff
        else
           gamma_half = cfp_gamma_fun(0.5_cfp)
           G_coeff = G_coeff*twopi*cfp_gamma_fun(0.5_cfp+l1+1)*gamma_half/(cfp_gamma_fun(0.5_cfp+l2+1)*cfp_gamma_fun(0.5_cfp+l3+1))
        endif

   end function G_coeff 

end module coupling_obj
