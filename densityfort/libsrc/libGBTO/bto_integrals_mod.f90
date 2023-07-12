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
!> This module contains routines to calculate the 1-electron integrals in the B-spline basis.
!> \todo Add checking that the radial quadratures have high enough order to get exact results.
!> Integrals involving radial B-splines with a non-zero derivative at r=A (i.e. at the origin of the B-spline grid) are not evaluated.
module bto_integrals_mod
   use utils
   use precisn
   use basis_data_generic_mod
   use coupling_obj
   use bspline_grid_mod
   use symmetry
!UNCOMMENT TO DEBUG
!   use cgto_pw_expansions_mod, only: dbg_cgto, init_dbg
   implicit none

   private

   public BB_initialize, BB_shell_integrals, construct_bspline_quadrature_grid

   !> Used to get various coupling coefficients, mostly the real Gaunt coefficients.
   type(couplings_type), private :: cpl

   !> Variables used to integrate the B-splines on a given quadrature grid.
   real(kind=cfp), allocatable, private :: r_points(:), weights(:), B_vals(:,:,:), temp_r(:), bspline_boundary_val(:,:)

   !> Variables used to integrate the B-splines on the quadrature grid prepared for the Legendre expansion.
   real(kind=cfp), allocatable, private :: r_points_leg(:), weights_leg(:), B_vals_leg(:,:,:), radial_prod(:)

   !> Maps start and end of each B-spline to the quadrature points of the radial grid.
   integer, allocatable, private :: bspline_start_end(:,:), bspline_start_end_leg(:,:)

   !> prop_on_grid: array containing the values sqrt(fourpi/(2*l+1.0_cfp))*r**l on the quadrature grid r_points for all l up to max_prop_l
   real(kind=cfp), allocatable, private :: prop_on_grid(:,:), leg_on_grid(:,:), re_sph_harm_center(:,:)

   !> Symmetry and nuclear data.
   type(symmetry_obj) :: symmetry_data_saved

contains

   subroutine BB_shell_integrals(shell_A,shell_B,starting_index_A,starting_index_B,&
                                     &a,max_property_l,property_center,symmetry_data,&
                                     &olap_column,kei_column,prop_column,nai_column,one_elham_column,int_index,integrals)
!UNCOMMENT TO DEBUG
!      use cgto_hgp, only: sph_nari
!      use gto_routines, only: compare_print_1el_ints
      use symmetry
      use phys_const, only: fourpi
      implicit none
      type(BTO_shell_data_obj), target, intent(in) :: shell_A, shell_B
      real(kind=cfp), intent(in) :: a, property_center(3)
      type(symmetry_obj), intent(in) :: symmetry_data
      integer, intent(in) :: starting_index_A, starting_index_B, olap_column,kei_column,nai_column,prop_column,&
                             one_elham_column,max_property_l
      !We assume that these two arrays have been allocated to the appropriate dimensions:
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)

      integer :: sph_shell_A, sph_shell_B, A_ind, B_ind, l, m, ind, i, shell_A_ind, shell_B_ind
      real(kind=cfp) :: bloch_el, radial_olap, radial_kei_1, radial_kei_2, kei, fac
      integer :: m1, m2, column, n_ang_prop, lm, terms, A_ind_leg, B_ind_leg
      real(kind=cfp) :: radial_prop, s
!UNCOMMENT TO DEBUG
!     type(CGTO_shell_data_obj) :: tmp_cgto_shell_A, tmp_cgto_shell_B
!     real(kind=cfp) :: tmp_A, tmp_B
!     real(kind=cfp), allocatable :: integrals_cgto(:,:), nari(:)
!     integer, allocatable :: int_index_cgto(:,:)
!     character(len=4) :: tag="BB"

         if (one_elham_column > 0 .and. (kei_column == 0 .or. nai_column == 0)) then
            call xermsg ('bto_gto_integrals_mod', 'BB_shell_integrals', &
                         'If one electron Hamiltonian integrals are requested then the kinetic energy and &
                         &nuclear attraction integrals must be requested too.', 1, 1)
         end if

         shell_A_ind = shell_A%bspline_index
         shell_B_ind = shell_B%bspline_index
         !todo DBG
         !if (shell_A%l .eq. 1) shell_A_ind = 10
         !if (shell_B%l .eq. 1) shell_B_ind = 10

         !Integrate over the radial interval where the B-splines overlap:
         A_ind = max(bspline_start_end(1,shell_A_ind),bspline_start_end(1,shell_B_ind))
         B_ind = min(bspline_start_end(2,shell_A_ind),bspline_start_end(2,shell_B_ind))

         sph_shell_A = 2*shell_A%l+1
         sph_shell_B = 2*shell_B%l+1
         terms = sph_shell_A*sph_shell_B

         !generate the int_index values:
         call index_1el(sph_shell_A,starting_index_A,sph_shell_B,starting_index_B,1,int_index)

         if (B_ind .le. A_ind) then
            if (olap_column > 0) integrals(1:terms,olap_column) = 0.0_cfp
            if (kei_column > 0) integrals(1:terms,kei_column) = 0.0_cfp
            if (nai_column > 0) integrals(1:terms,nai_column) = 0.0_cfp
            if (prop_column > 0) integrals(1:terms,prop_column:prop_column-1+(max_property_l+1)**2) = 0.0_cfp
            if (one_elham_column > 0) integrals(1:terms,one_elham_column) = 0.0_cfp
            return
         endif

         if (olap_column > 0) then
            integrals(1:terms,olap_column) = 0.0_cfp

            !Radial B-splines must overlap and the angular momenta of the shells A,B must be the same for the integrals to be non-zero:
            if (shell_A%l .eq. shell_B%l) then
               l = shell_A%l
   
               !Radial overlap integral:
               radial_olap = sum(B_vals(1,A_ind:B_ind,shell_A_ind)*B_vals(1,A_ind:B_ind,shell_B_ind)*weights(A_ind:B_ind))

               !Save the integrals to the appropriate place in the final array of integrals:
               !we use the fact that the overlap and kei are non-zero only for (l1,m1) .eq. (l2,m2) = (l,m)
               do m=-l,l
                  ind = m+l+1
                  ind = ind + sph_shell_A*(ind-1) !index of the (l1,m1) .eq. (l2,m2) = (l,m) pair
                  integrals(ind,olap_column) = radial_olap
               enddo !m
            endif

         endif

         if (kei_column > 0) then
            integrals(1:terms,kei_column) = 0.0_cfp
   
            !Radial B-splines must overlap and the angular momenta of the shells A,B must be the same for the integrals to be non-zero:
            if (shell_A%l .eq. shell_B%l) then
               l = shell_A%l
   
               !Two integrals that enter the evaluation of the KEI for a pair of BTOs (shell_A,shell_B) combination:
               radial_kei_1 = sum(B_vals(1,A_ind:B_ind,shell_A_ind) &
                                * B_vals(2,A_ind:B_ind,shell_B_ind) * weights(A_ind:B_ind))
               radial_kei_2 = sum(B_vals(1,A_ind:B_ind,shell_A_ind) &
                                * B_vals(1,A_ind:B_ind,shell_B_ind) * temp_r(A_ind:B_ind) * weights(A_ind:B_ind))
   
               kei = -0.5_cfp*(radial_kei_1 - l*(l+1)*radial_kei_2)
   
               !Calculate the Bloch element for the (shell_A,shell_B) combination, (most of them are zero but their calculation is so fast it doesn't matter):
               bloch_el = 0.5_cfp*(bspline_boundary_val(1,shell_A_ind)*bspline_boundary_val(2,shell_B_ind)-&
                                  &bspline_boundary_val(3,shell_A_ind)*bspline_boundary_val(4,shell_B_ind)) !we calculate the general Bloch term including the terms for both: r=r1 and r=r2.
   
               kei = kei + bloch_el
      
               !Save the integrals to the appropriate place in the final array of integrals:
               !we use the fact that the overlap and kei are non-zero only for (l1,m1) .eq. (l2,m2) = (l,m)
               do m=-l,l
                  ind = m+l+1
                  ind = ind + sph_shell_A*(ind-1) !index of the (l1,m1) .eq. (l2,m2) = (l,m) pair
                  integrals(ind,kei_column) = kei
               enddo !m
            endif

         endif

         if (prop_column > 0) then
            n_ang_prop = (max_property_l+1)**2
   
            do l=0,max_property_l
   
               !Radial property integral:
               radial_prop = sum(B_vals(1,A_ind:B_ind,shell_A_ind) &
                               * B_vals(1,A_ind:B_ind,shell_B_ind) * prop_on_grid(A_ind:B_ind,l) * weights(A_ind:B_ind))
   
               do m=-l,l
                  lm = l*l+l+m+1
   
                  column = prop_column-1+lm
                  integrals(1:terms,column) = 0.0_cfp
                  ind = 0
                  do m2=-shell_B%l,shell_B%l
                     do m1=-shell_A%l,shell_A%l
                        ind = ind + 1
                        integrals(ind,column) = radial_prop*cpl%rgaunt(l,shell_A%l,shell_B%l,m,m1,m2)
                     enddo !m1
                  enddo !m2
   
               enddo !m
   
            enddo !l
         endif

         if (nai_column > 0) then
            integrals(1:terms,nai_column) = 0.0_cfp

            !Integrate over the radial interval where the B-splines overlap:
            A_ind_leg = max(bspline_start_end_leg(1,shell_A_ind),bspline_start_end_leg(1,shell_B_ind))
            B_ind_leg = min(bspline_start_end_leg(2,shell_A_ind),bspline_start_end_leg(2,shell_B_ind))

            radial_prod(A_ind_leg:B_ind_leg) = B_vals_leg(1,A_ind_leg:B_ind_leg,shell_A_ind) &
                                             * B_vals_leg(1,A_ind_leg:B_ind_leg,shell_B_ind) * weights_leg(A_ind_leg:B_ind_leg)

!UNCOMMENT TO DEBUG
!            A_ind_leg = 1
!            B_ind_leg = maxval(bspline_start_end_leg(2,:))
!            call tmp_cgto_shell_A%make_space(1)
!            tmp_cgto_shell_A%l = shell_A%l
!            tmp_cgto_shell_A%center = 0.0_cfp
!            tmp_cgto_shell_A%number_of_primitives = 1
!            tmp_cgto_shell_A%exponents(1) = 3.0_cfp
!            tmp_cgto_shell_A%contractions(1) = 1.0_cfp
!            tmp_cgto_shell_A%non_zero_at_boundary = .false.
!            tmp_cgto_shell_A%number_of_functions = 2*tmp_cgto_shell_A%l+1
!            call tmp_cgto_shell_A%normalize
!            tmp_A = sqrt(fourpi/(2*shell_A%l+1.0_cfp))*tmp_cgto_shell_A%norm*tmp_cgto_shell_A%norms(1)
!
!            call tmp_cgto_shell_B%make_space(1)
!            tmp_cgto_shell_B%l = shell_B%l
!            tmp_cgto_shell_B%center = 0.0_cfp
!            tmp_cgto_shell_B%number_of_primitives = 1
!            tmp_cgto_shell_B%exponents(1) = 2.5_cfp
!            tmp_cgto_shell_B%contractions(1) = 1.0_cfp
!            tmp_cgto_shell_B%non_zero_at_boundary = .false.
!            tmp_cgto_shell_B%number_of_functions = 2*tmp_cgto_shell_B%l+1
!            call tmp_cgto_shell_B%normalize
!            tmp_B = sqrt(fourpi/(2*shell_B%l+1.0_cfp))*tmp_cgto_shell_B%norm*tmp_cgto_shell_B%norms(1)
!            do i=A_ind_leg,B_ind_leg
!               radial_prod(i) = tmp_A*r_points_leg(i)**tmp_cgto_shell_A%l*exp(-tmp_cgto_shell_A%exponents(1)*r_points_leg(i)**2)&
!                              &*tmp_B*r_points_leg(i)**tmp_cgto_shell_B%l*exp(-tmp_cgto_shell_B%exponents(1)*r_points_leg(i)**2)*r_points_leg(i)**2*weights_leg(i)
!            enddo !i

            ind = 0 
            do m2=-shell_B%l,shell_B%l
               do m1=-shell_A%l,shell_A%l
                  !Accumulate the final integral from the contributions of the terms from the Leg. expansion:
                  s = 0.0_cfp
                  do l=abs(shell_A%l-shell_B%l),shell_A%l+shell_B%l
                     do m=-l,l
                        lm = l*l+l+m+1
                        fac = cpl%rgaunt(l,shell_A%l,shell_B%l,m,m1,m2)
                        s = s + fac*sum(leg_on_grid(A_ind_leg:B_ind_leg,lm)*radial_prod(A_ind_leg:B_ind_leg))
                     enddo !m
                  enddo !l
                  !save in the order: (m1,m2)
                  ind = ind + 1
                  integrals(ind,nai_column) = s
               enddo !m2
            enddo !m1
!UNCOMMENT TO DEBUG
!            !loop over the nuclei.
!            ind = (2*shell_A%l+1)*(2*shell_B%l+1)
!            allocate(nari(ind),integrals_cgto(ind,nai_column),int_index_cgto(2,ind))
!    
!            integrals_cgto(1:ind,nai_column) = 0.0_cfp
!            do i=1,symmetry_data_saved%no_nuc
!    
!            call sph_nari(tmp_cgto_shell_A%number_of_primitives,tmp_cgto_shell_A%center(1),tmp_cgto_shell_A%center(2),tmp_cgto_shell_A%center(3),tmp_cgto_shell_A%norm,tmp_cgto_shell_A%norms,tmp_cgto_shell_A%l,tmp_cgto_shell_A%exponents,tmp_cgto_shell_A%contractions,starting_index_A,&
!                         &tmp_cgto_shell_B%number_of_primitives,tmp_cgto_shell_B%center(1),tmp_cgto_shell_B%center(2),tmp_cgto_shell_B%center(3),tmp_cgto_shell_B%norm,tmp_cgto_shell_B%norms,tmp_cgto_shell_B%l,tmp_cgto_shell_B%exponents,tmp_cgto_shell_B%contractions,starting_index_B,&
!                         &symmetry_data_saved%nucleus(i)%center(1),symmetry_data_saved%nucleus(i)%center(2),symmetry_data_saved%nucleus(i)%center(3), nari,int_index_cgto)
!               !accumulate the results into the nari_nuc array holding the final results and multiply the contributions by the nuclear charge.
!               integrals_cgto(1:ind,nai_column) = integrals_cgto(1:ind,nai_column) - symmetry_data_saved%nucleus(i)%charge*nari(1:ind)
!            enddo !i
!   
!            call compare_print_1el_ints(tag,integrals,int_index,integrals_cgto,int_index_cgto,ind,nai_column)

         endif

         if (one_elham_column > 0) then
            integrals(1:terms,one_elham_column) = integrals(1:terms,kei_column) + integrals(1:terms,nai_column)
         endif

   end subroutine BB_shell_integrals

   !> Calculates: r_points, weights, B_vals, temp_r, bspline_boundary_val, bspline_start_end, prop_on_grid, leg_on_grid, re_sph_harm_center.
   subroutine BB_initialize(bspline_grid,max_bspline_l,max_prop_l,symmetry_data)
      use symmetry
      use sort, only: sort_float
      use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
      use bspline_base
      use phys_const, only: fourpi
      use special_functions, only: cfp_resh
      use const, only: epsabs
      implicit none
      type(bspline_grid_obj), intent(inout) :: bspline_grid
      integer, optional, intent(in) :: max_prop_l,max_bspline_l
      type(symmetry_obj), intent(in) :: symmetry_data

      integer :: i, j, k, l, m, n_total_points, err, ind, lm, n_total_points_leg, n_centers
      real(kind=cfp) :: bto_norm, fac, R, r1, r2, val
      real(kind=cfp), allocatable :: xlm(:,:), list(:,:)

!UNCOMMENT TO DEBUG
!         call init_dbg(bspline_grid)

         symmetry_data_saved = symmetry_data

         !Precalculate the coupling coefficients: real Gaunt cfs.
         call cpl%prec_cgaunt(max(2*max_bspline_l,max_prop_l))

         !Construct the quadrature grid that will be used to integrate over all pairs of B-splines
         call construct_bspline_quadrature_grid(bspline_grid%knots,x_10,w_10,n_10,4,r_points,weights,n_total_points)

         !Find the start and end of each B-spline in the quadrature grid.
         call map_knots_to_grid(bspline_grid%knots,bspline_grid%order,bspline_grid%n,r_points,bspline_start_end)

         !Construct the quadrature grid that will be used to integrate over all pairs of B-splines when constructing the nuclear attraction integrals
         k = size(bspline_grid%knots)+symmetry_data%no_nuc
         allocate(list(k,1),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','prepare_bspline_basis','Memory allocation 0 failed.',err,1)

         do k=1,size(bspline_grid%knots)
            list(k,1) = bspline_grid%knots(k)
         enddo
         k = size(bspline_grid%knots)

         do i=1,symmetry_data%no_nuc
            k = k + 1
            list(k,1) = sqrt(dot_product(symmetry_data%nucleus(i)%center(1:3),symmetry_data%nucleus(i)%center(1:3)))
         enddo !i

         call sort_float(k,1,list)

         call construct_bspline_quadrature_grid(list(1:k,1),x_10,w_10,n_10,4,r_points_leg,weights_leg,n_total_points_leg)
   
         !Find the start and end of each B-spline in the quadrature grid.
         call map_knots_to_grid(bspline_grid%knots,bspline_grid%order,bspline_grid%n,r_points_leg,bspline_start_end_leg)
         
         if (allocated(B_vals)) deallocate(B_vals)
         if (allocated(B_vals_leg)) deallocate(B_vals_leg)
         if (allocated(temp_r)) deallocate(temp_r)
         if (allocated(bspline_boundary_val)) deallocate(bspline_boundary_val)

         allocate(B_vals(2,n_total_points,bspline_grid%n), &
                  B_vals_leg(2,n_total_points_leg,bspline_grid%n), &
                  temp_r(n_total_points), &
                  bspline_boundary_val(4,bspline_grid%n),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','prepare_bspline_basis','Memory allocation 1 failed.',err,1)
         B_vals = 0.0_cfp
         bspline_grid%bcoef = 0.0_cfp
         bspline_boundary_val = 0.0_cfp
         do ind=1,bspline_grid%n
            bspline_grid%bcoef(ind) = 1.0_cfp
            call bspline_grid%bspline_range(ind,r1,r2)
    
            !Calculate norm of the radial B-spline
            bto_norm = bspline_grid%normalize(ind)

!UNCOMMENT TO DEBUG
!            bspline_start_end(1,ind) = 1
!            bspline_start_end(2,ind) = n_total_points
!            bspline_start_end_leg(1,ind) = 1
!            bspline_start_end_leg(2,ind) = n_total_points_leg
!            bto_norm = dbg_cgto(ind)%norm*dbg_cgto(ind)%norms(1)
!            l = dbg_cgto(ind)%l
    
            !Evaluate it on the grid
            do i=bspline_start_end(1,ind),bspline_start_end(2,ind)
               !no derivative:
               B_vals(1,i,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                  bspline_grid % bcoef, &
                                                  bspline_grid % n, &
                                                  bspline_grid % order, 0, r_points(i), &
                                                  bspline_grid % inbv, &
                                                  bspline_grid % work)
               !2nd derivative:
               B_vals(2,i,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                  bspline_grid % bcoef, &
                                                  bspline_grid % n, &
                                                  bspline_grid % order, 2, r_points(i), &
                                                  bspline_grid % inbv, &
                                                  bspline_grid % work)

!UNCOMMENT TO DEBUG
!               B_vals(1,i,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*r_points(i)**l*exp(-dbg_cgto(ind)%exponents(1)*r_points(i)**2)*r_points(i)
!               B_vals(2,i,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*(exp(-dbg_cgto(ind)%exponents(1)*r_points(i)**2)*r_points(i)**(-1+l)&
!                                     &*(l*(1+l)-2*dbg_cgto(ind)%exponents(1)*(3+2*l)*r_points(i)**2+4*dbg_cgto(ind)%exponents(1)**2*r_points(i)**4))
            enddo !i

            !Evaluate it on the grid for Legendre resolution
            do i=bspline_start_end_leg(1,ind),bspline_start_end_leg(2,ind)
               !no derivative:
               B_vals_leg(1,i,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                      bspline_grid % bcoef, &
                                                      bspline_grid % n, &
                                                      bspline_grid % order, 0, r_points_leg(i), &
                                                      bspline_grid % inbv, &
                                                      bspline_grid % work)
               !2nd derivative:
               B_vals_leg(2,i,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                      bspline_grid % bcoef, &
                                                      bspline_grid % n, &
                                                      bspline_grid % order, 2, r_points_leg(i), &
                                                      bspline_grid % inbv, &
                                                      bspline_grid % work)
!UNCOMMENT TO DEBUG
!               B_vals_leg(1,i,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*r_points_leg(i)**l*exp(-dbg_cgto(ind)%exponents(1)*r_points_leg(i)**2)*r_points_leg(i)
!               B_vals_leg(2,i,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*(exp(-dbg_cgto(ind)%exponents(1)*r_points_leg(i)**2)*r_points_leg(i)**(-1+l)&
!                                     &*(l*(1+l)-2*dbg_cgto(ind)%exponents(1)*(3+2*l)*r_points_leg(i)**2+4*dbg_cgto(ind)%exponents(1)**2*r_points_leg(i)**4))
            enddo !i
   
            !Calculate the value and the first derivative of the B-spline at the boundary:
            bspline_boundary_val(1,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                           bspline_grid % bcoef, &
                                                           bspline_grid % n, &
                                                           bspline_grid % order, 0, r2, &
                                                           bspline_grid % inbv, &
                                                           bspline_grid % work)
            bspline_boundary_val(2,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                           bspline_grid % bcoef, &
                                                           bspline_grid % n, &
                                                           bspline_grid % order, 1, r2, &
                                                           bspline_grid % inbv, &
                                                           bspline_grid % work)

            !Calculate the value and the first derivative of the B-spline at the start point of the grid:
            bspline_boundary_val(3,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                           bspline_grid % bcoef, &
                                                           bspline_grid % n, &
                                                           bspline_grid % order, 0, r1, &
                                                           bspline_grid % inbv, &
                                                           bspline_grid % work)
            bspline_boundary_val(4,ind) = bto_norm * bvalu(bspline_grid % knots, &
                                                           bspline_grid % bcoef, &
                                                           bspline_grid % n, &
                                                           bspline_grid % order, 1, r1, &
                                                           bspline_grid % inbv, &
                                                           bspline_grid % work)

!UNCOMMENT TO DEBUG
!            r1 = 0.0_cfp
!            r2 = bspline_grid%B
!            bspline_boundary_val(3,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*r1**l*exp(-dbg_cgto(ind)%exponents(1)*r1**2)*r1
!            bspline_boundary_val(4,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*exp(-dbg_cgto(ind)%exponents(1)*r1**2)*r1**l&
!                                         &*(1+l-2*dbg_cgto(ind)%exponents(1)*r1**2)
!            bspline_boundary_val(1,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*r2**l*exp(-dbg_cgto(ind)%exponents(1)*r2**2)*r2
!            bspline_boundary_val(2,ind) = bto_norm*sqrt(fourpi/(2*l+1.0_cfp))*exp(-dbg_cgto(ind)%exponents(1)*r2**2)*r2**l&
!                                         &*(1+l-2*dbg_cgto(ind)%exponents(1)*r2**2)

            bspline_grid%bcoef(ind) = 0.0_cfp !clean-up for the next
         enddo !ind
   
         !Precalculate 1/r**2 on the quadrature grid.
         !The weights could be included here but to keep the code transparent we don't do it here.
         temp_r(1:n_total_points) = 1.0_cfp/(r_points(1:n_total_points))**2

         !Precalculate the radial part of the property operator on the quadrature grid.
         if (allocated(prop_on_grid)) deallocate(prop_on_grid)
         allocate(prop_on_grid(n_total_points,0:max(max_prop_l,1)),stat=err)

         !precalculate sqrt(fourpi/(2*l+1.0_cfp))*r**l on the quadrature grid
         do l=0,max_prop_l
            fac = sqrt(fourpi/(2*l+1.0_cfp))
            do i=1,n_total_points
               prop_on_grid(i,l) = fac*r_points(i)**l
            enddo !i
         enddo !l

         n_centers = symmetry_data%no_nuc

         if (allocated(leg_on_grid)) deallocate(leg_on_grid)
         if (allocated(re_sph_harm_center)) deallocate(re_sph_harm_center)
         if (allocated(radial_prod)) deallocate(radial_prod)
         allocate(leg_on_grid(n_total_points_leg,(2*max_bspline_l+1)**2), &
                  xlm(-2*max_bspline_l:max(1,2*max_bspline_l),0:max(1,2*max_bspline_l)), &
                  re_sph_harm_center(n_centers,(2*max_bspline_l+1)**2), &
                  radial_prod(n_total_points_leg), stat = err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','nuclear_attraction_ints_bto','Memory allocation 1 failed.',err,1)

         !Precalculate the values of the real spherical harmonics for each nucleus position -R times the nuclear charge times the charge of electron:
         re_sph_harm_center = 0.0_cfp
         if (max_bspline_l > 0) then
            do k=1,n_centers
               if (abs(symmetry_data%nucleus(k)%charge) .le. epsabs) cycle
               call cfp_resh(xlm, symmetry_data % nucleus(k) % center(1), &
                                  symmetry_data % nucleus(k) % center(2), &
                                  symmetry_data % nucleus(k) % center(3), 2 * max_bspline_l)
               do l=0,2*max_bspline_l
                  do m=-l,l
                     lm = l*l+l+m+1
                     re_sph_harm_center(k,lm) = -symmetry_data%nucleus(k)%charge*xlm(m,l)
                  enddo !m
               enddo !l
            enddo !k
         else !l .eq. 0
            do k=1,n_centers
               re_sph_harm_center(k,1) = -symmetry_data%nucleus(k)%charge/sqrt(fourpi)
            enddo
         endif

         leg_on_grid = 0.0_cfp
         do i=1,n_centers
            if (abs(symmetry_data%nucleus(i)%charge) .le. epsabs) cycle
            R = sqrt(dot_product(symmetry_data%nucleus(i)%center(1:3),symmetry_data%nucleus(i)%center(1:3)))
            do j=1,n_total_points_leg
               do l=0,2*max_bspline_l
                  if (r_points_leg(j) .le. R) then
                     val = ((r_points_leg(j)/R)**l)/R*fourpi/(2*l+1.0_cfp)
                  else
                     val = ((R/r_points_leg(j))**l)/r_points_leg(j)*fourpi/(2*l+1.0_cfp)
                  endif
                  do m=-l,l
                     lm = l*l+l+m+1
                     leg_on_grid(j,lm) = leg_on_grid(j,lm) + val*re_sph_harm_center(i,lm)
                  enddo !m
               enddo !l
            enddo !j
         enddo !i

   end subroutine BB_initialize

   !> Constructs a quadrature grid for the B-spline basis described by a given knot sequence. The final quadrature grid comprises 
   !> effectively a series of quadratures over subintervals generated in between each pair of distinct knots.
   subroutine construct_bspline_quadrature_grid(knots,x,w,n,n_rng_knot,r_points,weights,n_total_points)
      use general_quadrature, only: gl_expand_A_B
      implicit none
      !IN/OUT:
      integer, intent(in) :: n,n_rng_knot
      real(kind=cfp), intent(in) :: knots(:), x(2*n+1), w(2*n+1)
      real(kind=cfp), allocatable :: r_points(:), weights(:)
      integer, intent(out) :: n_total_points

      integer :: i, j, n_points, n_knots, err, cnt
      real(kind=cfp) :: delta, A, B

         n_knots = size(knots)

         if (n_knots <= 1 .or. n == 0) then
            call xermsg ('bto_integrals_mod', 'construct_quadrature_grid', &
                         'Invalid knot grid or bad Gaussian quadrature rule.', 1, 1)
         end if
         if (n_rng_knot < 1) then
            call xermsg ('bto_integrals_mod', 'construct_quadrature_grid', &
                         'On input n_rng_knot < 1 but must be at least 1.', 2, 1)
         end if

         !Calculate how many quadrature points there will be in total.
         n_points = 2*n+1 !number of points in the elementary Gaussian quadrature.
         n_total_points = 0
         do i=1,n_knots-1
            if (knots(i+1) > knots(i)) n_total_points = n_total_points + n_points*n_rng_knot
         enddo !i

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)

         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','construct_quadrature_grid','Memory allocation failed.',err,1)
         r_points = 0.0_cfp; weights = 0.0_cfp

         !Construct the full grid expanding the elementary Gaussian rule on each distinct interval of knots.
         cnt = 0
         do i=1,n_knots-1
            if (knots(i+1) > knots(i)) then
               A = knots(i)
               delta = (knots(i+1)-knots(i))/n_rng_knot
               do j=1,n_rng_knot
                  if (j .eq. n_rng_knot) then
                     B = knots(i+1)
                  else
                     B = A + delta
                  endif
                  call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),A,B)
                  cnt = cnt + n_points
                  A = B
               enddo
            endif
         enddo !i

   end subroutine construct_bspline_quadrature_grid

   subroutine index_1el(sph_shell_A,starting_index_A,sph_shell_B,starting_index_B,n_repeat,int_index)
      implicit none
      integer, intent(in) :: sph_shell_A,starting_index_A,sph_shell_B,starting_index_B,n_repeat
      integer, allocatable :: int_index(:,:)

      integer :: i,j,k, ind

         ind = 0
         do k=1,n_repeat
            do i=1,sph_shell_B
               do j=1,sph_shell_A
                  ind = ind + 1
                  int_index(1,ind) = max(i-1+starting_index_B,j-1+starting_index_A)
                  int_index(2,ind) = min(i-1+starting_index_B,j-1+starting_index_A)
               enddo !j
            enddo !i
         enddo !k

   end subroutine index_1el

end module bto_integrals_mod
