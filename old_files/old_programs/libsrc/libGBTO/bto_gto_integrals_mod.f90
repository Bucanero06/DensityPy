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
module bto_gto_integrals_mod
  use general_quadrature
  use basis_data_generic_mod
  use bspline_grid_mod
  use common_obj, only: nucleus_type
  use coupling_obj
  use cgto_pw_expansions_mod, only: CGTO_shell_pw_expansion_obj, CGTO_shell_pair_pw_expansion_obj, init_CGTO_pw_expansions_mod, &
                                    cpl, precalculate_Xlm_for_nuclei, legendre_grid_r1_r2_obj, dbg_cgto, init_dbg
  use gto_routines, only: compare_print_2el_ints, compare_print_1el_ints
  use const, only: epsrel, epsabs, max_epstab, stdout, mib, line_len, Y_lm_size_threshold
  use precisn
  use utils, only: xermsg
  implicit none

  private

  public BG_shell_integrals, BG_initialize, BG_final
  public BG_mixed_2el_initialize, BBGG_shell_integrals, BGGG_shell_integrals, BGBG_shell_integrals, max_l_BGGG, max_l_BGBG, &
         lebedev_BGGG_shell_integrals, max_l_BG


  type(legendre_grid_r1_r2_obj) :: grid_r1_r2
  type(CGTO_shell_pw_expansion_obj), allocatable :: CGTO_pw(:)
  type(CGTO_shell_pair_pw_expansion_obj) :: GG_pair_pw
  logical :: keep_ab_cd_order
  integer :: max_l_pw = -1, indexing_method = 0

  real(kind=cfp), allocatable :: epstab(:,:), res3la(:,:) !used for extrapolations
  logical, allocatable :: done(:) !used for extrapolations
  real(kind=cfp), allocatable :: couplings(:) !Auxiliary used for BBGG_shell_integrals
  real(kind=cfp), allocatable :: Y_lm_mixed_from_disk(:) !Auxiliary used for BGGG and BGBG classes if Y_lm functions are saved on disk

  !$OMP THREADPRIVATE(GG_pair_pw, epstab, res3la, done, couplings, Y_lm_mixed_from_disk)

  !> Maximum L to use in the Legendre expansion of the Coulomb potential when calculating the 1-electron integrals.
  integer :: max_l_legendre_1el = 40
  !> Maximum L to use in the Legendre expansion of the Coulomb potential when calculating the 2-electron integrals.
  integer :: max_l_legendre_2el = 40

  !> Convergence parameters specifying the requested relative precision.
  real(kind=cfp) :: overflow = -1.0_cfp, prec_goal = -1.0_cfp

  !> Maximum L that was needed to convere integrals in the BG, BGGG and BGBG class.
  integer, protected :: max_l_BGGG = -1, max_l_BGBG = -1, max_l_BG = -1

  !> Method chosen for calculation of the mixed integrals.
  integer :: mixed_ints_method = -1

  type(nucleus_type), allocatable :: nuclei(:)
  real(kind=cfp), allocatable :: Xlm_nuclei(:) !m,l,nucleus
  integer :: n_Xlm_nuclei

  logical :: module_initialized = .false.
  logical, parameter :: use_extrapolation = .false.

contains

  !> Precalculates norms of the radial B-splines.
  subroutine BG_initialize(inp_max_l_legendre_1el,inp_bspline_grid,inp_first_bspline_index,inp_max_bspline_l,inp_max_prop_l,&
                           inp_max_l_cgto,inp_nuclei,delta_r1,nai_column,inp_mixed_ints_method)
     implicit none
     type(bspline_grid_obj), intent(inout) :: inp_bspline_grid
     integer, intent(in) :: inp_max_l_legendre_1el,inp_first_bspline_index,inp_max_bspline_l,inp_max_prop_l,nai_column,&
                            inp_max_l_cgto,inp_mixed_ints_method
     type(nucleus_type), intent(in) :: inp_nuclei(:)
     real(kind=cfp), intent(in) :: delta_r1

     integer :: err
     logical, parameter :: only_on_bto_grid = .true. !We want to evaluate the CGTO pw expansion only on the radial domain of the B-spline grid.
     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied

!UNCOMMENT TO DEBUG
!        call init_dbg(inp_bspline_grid)

        if (allocated(nuclei)) deallocate(nuclei)
        allocate(nuclei,source=inp_nuclei,stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','BG_initialize','Memory allocation 1 failed.',err,1)

        mixed_ints_method = inp_mixed_ints_method

        if (mixed_ints_method .eq. 1) then
           max_l_legendre_1el = inp_max_l_legendre_1el
           write(stdout,'("Maximum L in the Leg. expansion of the Coulomb potential: ",i4)') max_l_legendre_1el

           if (max_l_legendre_1el .le. 0) then
              call xermsg('bto_gto_integrals_mod', 'BG_initialize','On input max_l_legendre_1el was .le. 0.',1,1)
           endif

           overflow = F1MACH(2,cfp_dummy)
           prec_goal = F1MACH(4,cfp_dummy)*100

           if (use_extrapolation) write(stdout,'("Requested relative precision for the extrapolated integrals: ",e25.15)') prec_goal

           if (max_l_legendre_1el < 2*inp_max_bspline_l) then
              call xermsg ('bto_gto_integrals_mod', 'BG_initialize', &
                           'On input max_l_legendre_1el was too small: it must be at least 2*max_bspline_l.', 2, 1)
           end if

           !max_bspline_l is assumed to be set to the largest angular momentum in the whole BTO basis.
           max_l_pw = max(inp_max_prop_l,inp_max_bspline_l,0)+max_l_legendre_1el

           !Real spherical harmonics for the nuclei: result in Xlm_nuclei.
           call precalculate_Xlm_for_nuclei(nuclei,max_l_legendre_1el,Xlm_nuclei,n_Xlm_nuclei)

        elseif (mixed_ints_method .eq. 2) then

           !max_bspline_l is assumed to be set to the largest angular momentum in the whole BTO basis.
           max_l_legendre_1el = 0
           max_l_pw = inp_max_bspline_l + max(inp_max_prop_l,0)

        elseif (mixed_ints_method .eq. 3) then

           !max_bspline_l is assumed to be set to the largest angular momentum in the whole BTO basis.
           max_l_legendre_1el = 0
           max_l_pw = inp_max_bspline_l + max(inp_max_prop_l,0)

        else
           call xermsg ('bto_gto_integrals_mod', 'BG_initialize', &
                        'mixed_ints_method is out of range: allowed values are 1,2,3.', mixed_ints_method, 1)
        endif

        !Initialize the module calculating the CGTO partial wave expansion on
        !the B-spline grid and on the grid needed to integrate accurately over
        !the domain of the CGTO.
        call init_CGTO_pw_expansions_mod(max_l_pw,inp_max_l_cgto)

        call grid_r1_r2%construct_r1_r2_grids(inp_bspline_grid,inp_first_bspline_index,inp_max_bspline_l,&
                        inp_max_prop_l,max(max_l_legendre_1el,0),nuclei,delta_r1,x_10,w_10,n_10,x_7,w_7,n_7,n_rng_knot)

        if (inp_mixed_ints_method .eq. 2) then
           call grid_r1_r2%eval_Xlm_on_lebedev_grid(80)
           print *,'lebedev eval'
        endif

        if (allocated(CGTO_pw)) deallocate(CGTO_pw)
        allocate(CGTO_pw(1),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_initialize','Memory allocation failed.',err,1)

        max_l_BGGG = -1
        max_l_BGBG = -1
        max_l_BG = -1

        module_initialized = .true.

        write(stdout,'(/,10X,"bto_gto_integrals_mod initialized")')

  end subroutine BG_initialize

  subroutine BG_final
     implicit none

     integer :: err, i

        call grid_r1_r2%final

        do i=1,size(CGTO_pw)
           call CGTO_pw(i)%final
        enddo !i

        deallocate(CGTO_pw,stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_final','Memory deallocation failed.',err,1)

        if (allocated(epstab)) deallocate(epstab)
        if (allocated(res3la)) deallocate(res3la)
        if (allocated(done)) deallocate(done)
        if (allocated(couplings)) deallocate(couplings)

        max_l_BGGG = -1
        max_l_BGBG = -1
        max_l_BG = -1
        module_initialized = .false.

  end subroutine BG_final

  subroutine BG_shell_integrals(cgto_shell,bto_shell,cgto_starting_index,bto_starting_index,cgto_shell_index,olap_column,&
                                kei_column,prop_column,nai_column,one_elham_column,int_index,integrals)
     use phys_const, only: fourpi
     use omp_lib
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: sph_olap_kei
!     use gto_routines, only: compare_print_1el_ints
     implicit none
     type(BTO_shell_data_obj), intent(inout) :: bto_shell
     type(CGTO_shell_data_obj), intent(inout) :: cgto_shell
     integer, intent(in) :: olap_column,kei_column,prop_column,nai_column,one_elham_column,cgto_starting_index,&
                            bto_starting_index,cgto_shell_index
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)

     real(kind=cfp) :: B, olap, kei, prop, fac, rmat_radius, bto_val(2), bloch, a0_square, &
                       bloch_full(2*bto_shell%l+1,2*cgto_shell%l+1)
     integer :: n_cgto_m, i, ind, m_ind, lm, base, l, m, n_threads, iam, lp_mp, lp, mp, j, terms, p
     integer :: first_point, last_point
     logical :: do_nai = .false.
     logical, parameter :: only_on_bto_grid = .true. !We want to evaluate the CGTO pw expansion only on the radial domain of the B-spline grid.
     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied
!UNCOMMENT TO DEBUG
!     integer :: lena, la, ind_a, lenb, lb, ind_b
!     real(kind=cfp), allocatable :: anorms(:), bnorms(:), acoefs(:), bcoefs(:), aexps(:), bexps(:), cgto_int(:,:)
!     real(kind=cfp) :: acnorm, bcnorm, xa,ya,za, xb,yb,zb, olap_tail,kei_tail_ab,bloch_ab
!     integer, allocatable :: cgto_int_index(:,:)
!     character(len=4) :: tag = "BG"

        if (.not. module_initialized) then
            call xermsg ('bto_gto_integrals_mod', 'BG_shell_integrals', &
                         'The module bto_gto_integrals_mod has not been initialized. Run BG_initialize first.', 1, 1)
        end if

        if (bto_shell%l > grid_r1_r2%max_bspline_l) then
            call xermsg ('bto_gto_integrals_mod', 'BG_shell_integrals', &
                         'Attempt to calculate integrals for BTO L larger than what the initialization was done for.', 2, 1)
        end if

        if (one_elham_column > 0 .and. (kei_column == 0 .or. nai_column == 0)) then
            call xermsg ('bto_gto_integrals_mod', 'BG_shell_integrals', &
                         'If one electron Hamiltonian integrals are requested then the kinetic energy &
                         &and nuclear attraction integrals must be requested too.', 3, 1)
        end if
!
!------ Recalculate the pw-CGTO expansion and/or the B-spline quadrature grids only if the CGTO has changed: for this reason it is most advantageous to loop over the shell pairs outside of this routine
!       in such a way that minimizes the repetition of the CGTO shells (see one_electron_integrals).
!
        if (CGTO_pw(1)%cgto_shell_index .ne. cgto_shell_index) then
           call CGTO_pw(1)%init_CGTO_shell_pw_expansion(cgto_shell,cgto_shell_index)
           call CGTO_pw(1)%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           !Evaluate the pw expansion on the r1 grid and on the grid of knots: needed for the Bloch operator.
           if (mixed_ints_method .eq. 2) then
              call CGTO_pw(1) % eval_CGTO_shell_pw_expansion (grid_r1_r2 % bspline_grid % knots, &
                                                              grid_r1_r2 % max_bspline_l, &
                                                              grid_r1_r2 % max_prop_l, 0)
              !Needed for the NAI using the angular quadrature
              call CGTO_pw(1)%eval_at_lebedev_points(grid_r1_r2)
           elseif (mixed_ints_method .eq. 1) then
              call CGTO_pw(1) % eval_CGTO_shell_pw_expansion (grid_r1_r2 % bspline_grid % knots, &
                                                              grid_r1_r2 % max_bspline_l, &
                                                              grid_r1_r2 % max_prop_l, &
                                                              grid_r1_r2 % max_l_legendre)
           elseif (mixed_ints_method .eq. 3) then
              call CGTO_pw(1) % eval_CGTO_shell_pw_expansion (grid_r1_r2 % bspline_grid % knots, &
                                                              grid_r1_r2 % max_bspline_l, &
                                                              grid_r1_r2 % max_prop_l, &
                                                              grid_r1_r2 % max_bspline_l)
              call CGTO_pw(1)%eval_NAI_X_lm_projections(grid_r1_r2,nuclei)
           endif
           !todo temp
           !call CGTO_pw(1)%write
        endif
!
!------ Recalculate the pw-CGTO expansion and/or the B-spline quadrature grids only if the CGTO has changed: for this reason it is most advantageous to loop over the shell pairs outside of this routine
!       in such a way that minimizes the repetition of the CGTO shells (see one_electron_integrals).
!
        if (nai_column > 0) do_nai = .true. !calculate the variables needed for NAI evaluation
!
!------ Calculate the radial integrals for all BTOs: at each quadrature point calculate its contribution to integrals over all BTOs non-zero at that point.
!
        first_point = grid_r1_r2%bspline_start_end_r1(1,bto_shell%bspline_index)
        last_point = grid_r1_r2%bspline_start_end_r1(2,bto_shell%bspline_index)
        l = bto_shell%l
        ind = bto_shell%bspline_index
        !todo DBG
        !if (bto_shell%l .eq. 1) ind = 10
        base = l**2 !number of preceeding (l,m) combinations up to l = bto_shell%l-1.
        rmat_radius = bto_shell%bspline_grid%B
        a0_square = rmat_radius*rmat_radius
        n_cgto_m = 2*cgto_shell%l+1
        terms = (2*l+1)*n_cgto_m

        if (first_point > last_point) then
           !This radial B-spline doesn't overlap with the radial domain of the CGTO: the integrals are negligible.
           if (olap_column > 0) then
              integrals(1:terms,olap_column) = 0.0_cfp
           endif
           if (kei_column > 0) then
              integrals(1:terms,kei_column) = 0.0_cfp
           endif
           if (nai_column > 0) then
              integrals(1:terms,nai_column) = 0.0_cfp
           endif
           if (one_elham_column > 0) then
              integrals(1:terms,one_elham_column) = 0.0_cfp
           endif
           if (prop_column > 0) then
              do lp_mp=0,(grid_r1_r2%max_prop_l+1)**2-1
                 integrals(1:terms,prop_column+lp_mp) = 0.0_cfp
              enddo
           endif
           return
        endif
!
        ! Note: This parallel section needs to use DEFAULT(SHARED) to allow work with polymorphic objects with gfortran.
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP & SHARED(CGTO_pw,n_cgto_m,bto_shell,integrals,bloch_full,grid_r1_r2,&
        !$OMP &        cgto_shell,first_point,last_point,ind,l,base,rmat_radius,a0_square,olap_column, &
        !$OMP &        kei_column,prop_column,nai_column,bloch,terms) &
        !$OMP & PRIVATE(i,j,p,iam,n_threads,olap,m_ind,lm,kei,m,lp,fac,mp,prop,lp_mp,bto_val)
        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()
!
!------ Overlap integrals
        if (olap_column > 0) then
           !For all BTO m-angular values
           do m=-l,l
              lm = base + l+m+1
              if (mod(lm,n_threads) .ne. iam) cycle !work distribution
              !For all m-values of the CGTO function
              do m_ind=1,n_cgto_m
                 p = CGTO_pw(1)%non_neg_indices_l(m_ind,lm)
                 if (p .ne. 0) then
                    olap = sum(grid_r1_r2 % bto_radial_olap(first_point:last_point,ind) &
                                * CGTO_pw(1) % angular_integrals(first_point:last_point,p))
                 else
                    olap = 0.0_cfp
                 endif
                 !The integrals are stored in the order: (CGTO m, BTO m)
                 i = m_ind + n_cgto_m*(l+m)
                 integrals(i,olap_column) = olap
              enddo !m_ind
           enddo !m
        endif
        !$OMP BARRIER
!
!------ Kinetic energy integrals
        if (kei_column > 0) then

           !$OMP SINGLE
           bloch_full = 0.0_cfp
           if (cgto_shell%is_continuum() .and. bto_shell%l .eq. cgto_shell%l) then
              !Indices in the angular_integrals_at_knots of the r_start, r_end points   
              i = ind
              j = ind + bto_shell%bspline_grid%order
!UNCOMMENT TO DEBUG
!              i = 1
!              j = size(CGTO_pw(1)%angular_integrals_at_knots,1)
   
              !first derivative of the B-spline at its end points: only the B-splines at the end points of the interval may give rise to Bloch terms
              bto_val(1) = grid_r1_r2%bto_end_points(1,ind) !1/r*Jac*BTO_derivative at r=r_start
              bto_val(2) = grid_r1_r2%bto_end_points(2,ind) !1/r*Jac*BTO_derivative at r=r_end

              if (bto_val(1) .ne. 0.0_cfp .or. bto_val(2) .ne. 0.0_cfp) then !the Bloch term is zero
                 !calculate the Bloch terms for the m_BTO=m_CGTO combinations of angular parts of the CGTO and the BTO shells
                 !do m_ind=1,n_cgto_m
                    do m=-l,l
                       lm = l*l+l+m+1
                       m_ind = l+m+1
                       !angular_integrals_at_knots = projections of the CGTO at radial points corresponding to the BTO knots on the angular parts of the BTO.
                       !Bloch term for the end-points of the radial B-spline
                       p = CGTO_pw(1)%non_neg_indices_l_at_knots(m_ind,lm)
                       if (p .ne. 0) then
                          bloch_full(m+l+1,m_ind) = 0.5_cfp*(CGTO_pw(1)%angular_integrals_at_knots(j,p)*bto_val(2) &
                                                           - CGTO_pw(1)%angular_integrals_at_knots(i,p)*bto_val(1))
                       endif
                    enddo !m
                 !enddo !m_ind
              endif
           endif
           !$OMP END SINGLE
   
           !For all m-values of the CGTO function
           do m_ind=1,n_cgto_m
              !For all BTO angular m-values
              do m=-l,l
                 lm = base+l+m+1
                 if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                 p = CGTO_pw(1)%non_neg_indices_l(m_ind,lm)
                 if (p .ne. 0) then
                    kei = -0.5_cfp * sum(grid_r1_r2 % bto_radial_kei(first_point:last_point,ind,l) &
                                            * CGTO_pw(1) % angular_integrals(first_point:last_point,p))
                 else
                    kei = 0.0_cfp
                 endif
                 !Add the Bloch terms:
                 kei = kei + bloch_full(m+l+1,m_ind)
                 !The integrals are stored in the order: (CGTO m, BTO m)
                 i = m_ind + n_cgto_m*(l+m)
                 integrals(i,kei_column) = kei
              enddo !m
              !$OMP BARRIER
           enddo !m_ind
        endif
!
!------ Property integrals: (lp,mp) = property l,m values
        if (prop_column > 0) then
           !For all property (l,m)
           do lp=0,grid_r1_r2%max_prop_l
              do mp=-lp,lp
                 lp_mp = lp*lp + lp+mp+1
   
                 fac = sqrt(fourpi/(2*lp+1.0_cfp))
                 !For all BTO angular m-values
                 do m=-l,l
                    lm = base+l+m+1
                    if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                    !For all m-values of the CGTO function
                    do m_ind=1,n_cgto_m
                       p = CGTO_pw(1)%non_neg_indices_l_lp(m_ind,lp_mp,lm)
                       if (p .ne. 0) then
                          prop = fac * sum(grid_r1_r2 % bto_radial_prop(first_point:last_point,ind,lp) &
                                            * CGTO_pw(1) % gaunt_angular_integrals(first_point:last_point,p))
                       else
                          prop = 0.0_cfp
                       endif
                       !The integrals are stored in the order: (CGTO m,BTO m)
                       i = m_ind + n_cgto_m*(l+m)
                       integrals(i,prop_column+lp_mp-1) = prop
                    enddo !m_ind
                 enddo !m
                 !$OMP BARRIER
              enddo !mp
           enddo !lp
        endif
        !$OMP END PARALLEL
!
!------ Generate the basis function indices
!
        i = 0
        do m=1,2*l+1
           do m_ind=1,n_cgto_m
              i = i + 1
              int_index(1,i) = max(m-1+bto_starting_index,m_ind-1+cgto_starting_index)
              int_index(2,i) = min(m-1+bto_starting_index,m_ind-1+cgto_starting_index)
           enddo !m_ind
        enddo !m
!
!------ Nuclear attraction integrals
        if (nai_column > 0) then
           if (mixed_ints_method .eq. 2) then
              call lebedev_BG_nai_integrals(CGTO_pw(1),bto_shell,integrals,nai_column,&
                                            int_index,cgto_starting_index,bto_starting_index)
           elseif (mixed_ints_method .eq. 1) then
              call BG_nai_integrals(CGTO_pw(1),bto_shell,max_l_legendre_1el,integrals,&
                                    nai_column,int_index,cgto_starting_index,bto_starting_index)
           elseif (mixed_ints_method .eq. 3) then
              call BG_semi_analytic_nai_integrals(CGTO_pw(1),bto_shell,integrals,nai_column,&
                                                  int_index,cgto_starting_index,bto_starting_index)
           endif
        endif
!
!------ One electron Hamiltonian integrals
        if (one_elham_column > 0) then
           integrals(1:terms,one_elham_column) = integrals(1:terms,nai_column) + integrals(1:terms,kei_column)
        endif

!UNCOMMENT TO DEBUG
!        allocate(anorms(1), bnorms(cgto_shell%number_of_primitives))
!        allocate(aexps(1), bexps(cgto_shell%number_of_primitives))
!        allocate(acoefs(1), bcoefs(cgto_shell%number_of_primitives), cgto_int(terms,max(olap_column,kei_column)))
!        allocate(cgto_int_index(2,terms))
!        lena = 1
!        la = dbg_cgto(ind)%l
!        acnorm = dbg_cgto(ind)%norm
!        anorms(1) = dbg_cgto(ind)%norms(1)
!        aexps(1) = dbg_cgto(ind)%exponents(1)
!        acoefs(1) = dbg_cgto(ind)%contractions(1)
!        xa = 0.0_cfp; ya = 0.0_cfp; za = 0.0_cfp
!        ind_a = bto_starting_index
!   
!        lenb = cgto_shell%number_of_primitives
!        lb = cgto_shell%l
!        bcnorm = cgto_shell%norm
!        bnorms(1:lenb) = cgto_shell%norms(1:lenb)
!        bexps(1:lenb) = cgto_shell%exponents(1:lenb)
!        bcoefs(1:lenb) = cgto_shell%contractions(1:lenb)
!        xb = cgto_shell%center(1); yb = cgto_shell%center(2); zb = cgto_shell%center(3)
!        ind_b = cgto_starting_index
!   
!        call sph_olap_kei(lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs,ind_a, lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs,ind_b, olap_column,kei_column,cgto_int,cgto_int_index)
!        print *,'nint olap',cgto_shell_index,i
!        call compare_print_1el_ints(tag,integrals,int_index,cgto_int,cgto_int_index,terms,olap_column)
!        print *,'nint kei',cgto_shell_index,i
!        call compare_print_1el_ints(tag,integrals,int_index,cgto_int,cgto_int_index,terms,kei_column)

  end subroutine BG_shell_integrals
 
  !> requires that integrals() have been zeroed out before.
  subroutine lebedev_BG_nai_integrals(CGTO_pw,bto_shell,integrals,nai_column,int_index,cgto_starting_index,bto_starting_index)
     use common_obj, only: nucleus_type
     use phys_const, only: fourpi
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: sph_nari
     implicit none
     type(CGTO_shell_pw_expansion_obj), intent(in) :: CGTO_pw
     type(BTO_shell_data_obj), intent(in) :: bto_shell
     integer, intent(in) :: nai_column
     real(kind=cfp), allocatable :: integrals(:,:)
     integer, allocatable :: int_index(:,:) !only for debugging
     integer :: cgto_starting_index, bto_starting_index  !only for debugging

     integer :: n_nuclei, B_start, B_end, err, ind, i, j, n, la_ma, m_ind, ma, n_integrals, base, n_angular
     real(kind=cfp), allocatable :: coulomb_bto(:,:)
     real(kind=cfp) :: r_leb(3), r(3), x, coulomb, bto_part
!UNCOMMENT TO DEBUG
!     real(kind=cfp), allocatable :: integrals_cgto(:,:), nari(:)
!     integer, allocatable :: int_index_cgto(:,:)
!     integer :: starting_index_A, starting_index_B

        if (.not. module_initialized) then
            call xermsg ('bto_gto_integrals_mod', 'BG_nai_integrals', &
                         'The module bto_gto_integrals_mod has not been initialized. Run BG_initialize first.', 1, 1)
        end if
        
        B_start = grid_r1_r2%bspline_start_end_r1(1,bto_shell%bspline_index)
        B_end = grid_r1_r2%bspline_start_end_r1(2,bto_shell%bspline_index)
        allocate(coulomb_bto(grid_r1_r2%lebedev_order,B_start:B_end),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','lebedev_BG_nai_integrals','Memory allocation failed.',err,1)

        n_nuclei = size(nuclei)
        n_angular = grid_r1_r2%lebedev_order

        do i=B_start,B_end
           x = grid_r1_r2%r1(i)
           bto_part = grid_r1_r2%B_vals_r1(i,bto_shell%bspline_index)*grid_r1_r2%r1(i)*grid_r1_r2%w1(i)
           do j=1,grid_r1_r2%lebedev_order
              r_leb(1:3) = (/x*grid_r1_r2%leb_r1(j),x*grid_r1_r2%leb_r2(j),x*grid_r1_r2%leb_r3(j)/)
              coulomb = 0.0_cfp
              do n=1,n_nuclei
                 if (abs(nuclei(n)%charge) .le. epsabs) cycle
                 r(1:3) = r_leb(1:3) - nuclei(n)%center(1:3)
                 coulomb = coulomb - nuclei(n)%charge/sqrt(dot_product(r,r))
              enddo !n
              !Multiply by the radial part corresp. to the BTO and the Lebedev
              !quadrature weights.
              coulomb_bto(j,i) = coulomb*bto_part*grid_r1_r2%leb_w(j)*fourpi
           enddo !j
        enddo !i

        ind = 0
        n_angular = grid_r1_r2%lebedev_order
        integrals(1:(2*bto_shell%l+1)*(2*CGTO_pw%cgto_shell%l+1),nai_column) = 0.0_cfp
        do ma=-bto_shell%l,bto_shell%l
           la_ma = bto_shell%l*bto_shell%l+bto_shell%l+ma+1
           base = n_angular*(la_ma-1)
           do m_ind=1,2*CGTO_pw%cgto_shell%l+1
              ind = ind + 1
              do i=B_start,B_end
                 integrals(ind,nai_column) = integrals(ind,nai_column) &
                        + sum(coulomb_bto(1:n_angular,i) &
                                * grid_r1_r2 % Xlm_Lebedev(base+1:base+n_angular) &
                                * CGTO_pw % at_lebedev_points(1:n_angular,i,m_ind))
              enddo !i
           enddo !m_ind
        enddo !ma

  end subroutine lebedev_BG_nai_integrals

  subroutine BG_semi_analytic_nai_integrals(CGTO_pw,bto_shell,integrals,nai_column,int_index,cgto_starting_index,bto_starting_index)
     use common_obj, only: nucleus_type
     use phys_const, only: fourpi
     implicit none
     type(CGTO_shell_pw_expansion_obj), intent(in) :: CGTO_pw
     type(BTO_shell_data_obj), intent(in) :: bto_shell
     integer, intent(in) :: nai_column
     real(kind=cfp), allocatable :: integrals(:,:)
     integer, allocatable :: int_index(:,:)
     integer :: cgto_starting_index, bto_starting_index

     integer :: ind,ma,la_ma,n,B_start,B_end,bto_ind,lama,m_ind

        if (.not. module_initialized) then
            call xermsg ('bto_gto_integrals_mod', 'BG_semi_analytic_nai_integrals', &
                         'The module bto_gto_integrals_mod has not been initialized. Run BG_initialize first.', 1, 1)
        end if

        bto_ind = bto_shell%bspline_index
        n = grid_r1_r2%n1_total_points
        B_start = grid_r1_r2%bspline_start_end_r1(1,bto_ind)
        B_end = grid_r1_r2%bspline_start_end_r1(2,bto_ind)
        !print *,'BTO/GTO L',bto_shell%l,CGTO_pw%cgto_shell%l

        ind = 0
        do ma=-bto_shell%l,bto_shell%l
           lama = bto_shell%l*bto_shell%l+bto_shell%l+ma+1
           do m_ind=1,2*CGTO_pw%cgto_shell%l+1
              ind = ind + 1
              integrals(ind,nai_column) = sum(CGTO_pw%NAI_X_lm_projections(B_start:B_end,m_ind,lama) &
                    *grid_r1_r2%B_vals_r1(B_start:B_end,bto_ind)*grid_r1_r2%r1(B_start:B_end)*grid_r1_r2%w1(B_start:B_end))
              !print *,'new',ind,integrals(ind,nai_column)
              !write(50,'(i10,e25.15)') ind,integrals(ind,nai_column)
           enddo
        enddo !ma

  end subroutine BG_semi_analytic_nai_integrals

  subroutine BG_nai_integrals(CGTO_pw,bto_shell,max_l_legendre_1el,integrals,nai_column,int_index, &
                                cgto_starting_index,bto_starting_index)
     use common_obj, only: nucleus_type
     use phys_const, only: fourpi
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: sph_nari
     implicit none
     type(CGTO_shell_pw_expansion_obj), intent(in) :: CGTO_pw
     type(BTO_shell_data_obj), intent(in) :: bto_shell
     integer, intent(in) :: max_l_legendre_1el, nai_column
     real(kind=cfp), allocatable :: integrals(:,:)
     integer, allocatable :: int_index(:,:) !only for debugging
     integer :: cgto_starting_index, bto_starting_index  !only for debugging

     logical :: extrapolate, finite_leg_expansion
     integer :: l1, l2, min_l, max_l, dl, n_nuclei, B_start, B_end, err, ind, i, j, lm, l, m, la_ma, m_ind, ma, n_integrals, &
                n, n_wynn, nres, not_converged, p, bto_ind, n_prec, int_ind
     real(kind=cfp), allocatable :: r_12(:), f_l_i(:), f_lm(:,:)
     real(kind=cfp) :: R, abserr, res, rel_precision_1, rel_precision_2, fac, terms(max_epstab), max_val, nai
!UNCOMMENT TO DEBUG
!     real(kind=cfp), allocatable :: integrals_cgto(:,:), nari(:)
!     integer, allocatable :: int_index_cgto(:,:)
!     integer :: starting_index_A, starting_index_B
!     character(len=4) :: tag = "BGNI"

        if (.not. module_initialized) then
            call xermsg ('bto_gto_integrals_mod', 'BG_nai_integrals', &
                         'The module bto_gto_integrals_mod has not been initialized. Run BG_initialize first.', 1, 1)
        end if

        bto_ind = bto_shell%bspline_index
!UNCOMMENT TO DEBUG
        !todo DBG
!        if (bto_shell%l .eq. 1) bto_ind = 10

        extrapolate = use_extrapolation
        finite_leg_expansion = .false.

        min_l = 0
        max_l = max_l_legendre_1el
        !If the CGTO lies on CMS then the Legendre expansion is finite and exact: no extrapolation must be done.
        if (dot_product(CGTO_pw%cgto_shell%center,CGTO_pw%cgto_shell%center) .le. epsabs) then
           extrapolate = .false.
           min_l = abs(CGTO_pw%cgto_shell%l-bto_shell%l)
           max_l = CGTO_pw%cgto_shell%l+bto_shell%l
           finite_leg_expansion = .true.
        endif

        n_integrals = (2*bto_shell%l+1)*(2*CGTO_pw%cgto_shell%l+1)
        integrals(1:n_integrals,nai_column) = 0.0_cfp

        if (.not. allocated(epstab) .or. size(epstab,2) < n_integrals) then
           if (allocated(epstab)) deallocate(epstab)
           if (allocated(res3la)) deallocate(res3la)
           allocate(epstab(max_epstab,n_integrals),res3la(3,n_integrals),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_nai_integrals','Memory allocation 2 failed.',err,1)
        endif

        B_start = grid_r1_r2%bspline_start_end_r1(1,bto_ind)
        B_end = grid_r1_r2%bspline_start_end_r1(2,bto_ind)

        !Generate the function coming from the Legendre expansion and
        !sum it over all nuclei.
        allocate(r_12(B_start:B_end),f_l_i(B_start:B_end),f_lm(B_start:B_end,(max_l+1)**2),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_nai_integrals','Memory allocation failed.',err,1)
        n_nuclei = size(nuclei)
        f_lm = 0.0_cfp
        do i=1,n_nuclei
           ind = (i-1)*n_Xlm_nuclei
           if (abs(nuclei(i)%charge) .le. epsabs) cycle

           R = sqrt(dot_product(nuclei(i)%center,nuclei(i)%center))
           if (R <= epsabs) then
              !The nucleus is sitting on CMS:
              !the integral can be represented exactly
              !as a mixed BTO/GTO overlap integral
              !but with extra 1/r term in the integrand.

              !We also have to ensure that the extrapolation is not used.
              extrapolate = .false.
              finite_leg_expansion = .true.

              int_ind = 0
              !For all BTO m-angular values
              do ma=-bto_shell%l,bto_shell%l
                 la_ma = bto_shell%l*bto_shell%l+bto_shell%l+ma+1
                 !For all m-values of the CGTO function
                 do m_ind=1,2*CGTO_pw%cgto_shell%l+1
                    int_ind = int_ind + 1
                    p = CGTO_pw % non_neg_indices_l(m_ind,la_ma)
                    if (p .ne. 0) then
                       nai = sum(grid_r1_r2 % bto_radial_olap(B_start:B_end,bto_ind) &
                                   * (-nuclei(i)%charge) / grid_r1_r2 % r1(B_start:B_end) &
                                   * CGTO_pw % angular_integrals(B_start:B_end,p))
                    else
                       nai = 0.0_cfp
                    endif
                    !The integrals are stored in the order: (CGTO m, BTO m)
                    integrals(int_ind,nai_column) = integrals(int_ind,nai_column) + nai
                 enddo !m_ind
              enddo !m

              !Go to the next nucleus
              cycle

           end if

           do j=B_start,B_end
              if (CGTO_pw%r_points(j) .le. R) then
                 r_12(j) = CGTO_pw%r_points(j)/R
              else
                 r_12(j) = R/CGTO_pw%r_points(j)
              endif
           enddo !j

           lm = 0
           do l=0,max_l
              fac = fourpi/(2*l+1.0_cfp)
              do j=B_start,B_end
                 if (CGTO_pw%r_points(j) .le. R) then
                    f_l_i(j) = r_12(j)**(l+1) !r*f_l
                 else
                    f_l_i(j) = r_12(j)**(l) !r*f_l
                 endif
              enddo !j
              do m=-l,l
                 lm = lm + 1
                 do j=B_start,B_end
                    f_lm(j,lm) = f_lm(j,lm) + fac*f_l_i(j)*(-nuclei(i)%charge)*Xlm_nuclei(ind+lm)
                 enddo !j
              enddo !m
           enddo !l
        enddo !i

        if (extrapolate) then
           max_l = min(max_l,max_epstab-3)
           if (max_l < 6) stop "max_l must be at least 6 due to extrapolation"
           epstab(1:max_epstab,1:n_integrals) = 0.0_cfp
           res3la = 0.0_cfp
        endif

        l1 = max_l
        l2 = min_l
        dl = -1

        n = 0
        max_val = 0.0_cfp
        do l=l1,l2,dl

           n = n+1

           do m=-l,l
              lm = l*l+l+m+1
              
              f_lm(B_start:B_end,lm) = f_lm(B_start:B_end,lm) &
                                        * grid_r1_r2 % B_vals_r1(B_start:B_end,bto_ind) &
                                        * grid_r1_r2 % w1(B_start:B_end)

              int_ind = 0
              do ma=-bto_shell%l,bto_shell%l
                 la_ma = bto_shell%l*bto_shell%l+bto_shell%l+ma+1
                 do m_ind=1,2*CGTO_pw%cgto_shell%l+1
                    int_ind = int_ind + 1
                    p = CGTO_pw%non_neg_indices_l_lp(m_ind,la_ma,lm)
                    if (p .eq. 0) cycle
                    if (extrapolate) then
                       epstab(n,int_ind) = epstab(n,int_ind)&
                                         + sum(f_lm(B_start:B_end,lm)*CGTO_pw%gaunt_angular_integrals(B_start:B_end,p))
                    else
                       integrals(int_ind,nai_column) = integrals(int_ind,nai_column) &
                                                + sum(f_lm(B_start:B_end,lm)*CGTO_pw%gaunt_angular_integrals(B_start:B_end,p))
                    endif
!                    write(*,'("maxval",3i4,e25.15)') l,m,int_ind,maxval(abs(CGTO_pw%gaunt_angular_integrals(B_start:B_end,m_ind,lm,la_ma)))
                 enddo !m_ind
              enddo !ma
              
           enddo !m

           if (.not.(extrapolate) .and. .not.(finite_leg_expansion)) then
              if (n .eq. 1 .or. max_val .eq. 0.0_cfp) then
                 max_val = maxval(abs(integrals(1:n_integrals,nai_column)))
                 if (max_val .ne. 0.0_cfp) print *,'max_val',max_val
              endif
           endif

        enddo !l

        max_l_BG = max(max_l_BG,l)

        if (extrapolate) then
           do i=1,n_integrals
!              integrals(i,nai_column) = sum(epstab(1:max_l+1,i))
              terms = 0.0_cfp
              nres = 0
              do n=1,max_l+1
                 terms(n) = sum(epstab(1:n,i))
!                 write(*,'("t",2i6,e25.15)') i,n,terms(n)
                 n_wynn = n
                 call dqelg(n_wynn, terms(1:max_epstab), res, abserr, res3la(1:3,i), nres)
!                 write(*,'("r",3e25.15)') integrals(i,nai_column),res,abserr
                 p = n
              enddo !n
              integrals(i,nai_column) = res
           enddo !i
        endif
!        do i=1,n_integrals
!           print *,'old',i,integrals(i,nai_column)
!        enddo !i

!UNCOMMENT TO DEBUG
!        ind = (2*bto_shell%l+1)*(2*CGTO_pw%cgto_shell%l+1)
!        allocate(nari(ind),integrals_cgto(ind,nai_column),int_index_cgto(2,ind))
!        starting_index_A = cgto_starting_index
!        starting_index_B = bto_starting_index
!
!        integrals_cgto(1:ind,nai_column) = 0.0_cfp
!        do i=1,n_nuclei
!
!        call sph_nari(CGTO_pw%cgto_shell%number_of_primitives,CGTO_pw%cgto_shell%center(1),CGTO_pw%cgto_shell%center(2),CGTO_pw%cgto_shell%center(3),CGTO_pw%cgto_shell%norm,CGTO_pw%cgto_shell%norms,CGTO_pw%cgto_shell%l,CGTO_pw%cgto_shell%exponents,CGTO_pw%cgto_shell%contractions,starting_index_A,&
!                     &dbg_cgto(bto_ind)%number_of_primitives,dbg_cgto(bto_ind)%center(1),dbg_cgto(bto_ind)%center(2),dbg_cgto(bto_ind)%center(3),dbg_cgto(bto_ind)%norm,dbg_cgto(bto_ind)%norms,dbg_cgto(bto_ind)%l,dbg_cgto(bto_ind)%exponents,dbg_cgto(bto_ind)%contractions,starting_index_B,&
!                     &nuclei(i)%center(1),nuclei(i)%center(2),nuclei(i)%center(3), nari,int_index_cgto)
!           !accumulate the results into the nari_nuc array holding the final results and multiply the contributions by the nuclear charge.
!           integrals_cgto(1:ind,nai_column) = integrals_cgto(1:ind,nai_column) - nuclei(i)%charge*nari(1:ind)
!        enddo !i
!
!        print *,'nint nai'
!        call compare_print_1el_ints(tag,integrals,int_index,integrals_cgto,int_index_cgto,ind,nai_column)

  end subroutine BG_nai_integrals

  subroutine BG_mixed_2el_initialize(inp_max_l_legendre_2el,inp_bspline_grid,inp_first_bspline_index,inp_max_bspline_l,&
                        delta_r1,cgto_shells,inp_keep_ab_cd_order,mixed_ints_method,scratch_directory,indexing_method_inp)
     use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
     use omp_lib
     implicit none
     type(bspline_grid_obj), intent(inout) :: inp_bspline_grid
     integer, intent(in) :: inp_max_l_legendre_2el, inp_first_bspline_index, inp_max_bspline_l, &
                            mixed_ints_method, indexing_method_inp
     character(len=line_len), intent(in) :: scratch_directory
     logical, intent(in) :: inp_keep_ab_cd_order
     type(CGTO_shell_data_obj), intent(in) :: cgto_shells(:)
     real(kind=cfp), intent(in) :: delta_r1

     integer :: n_gq, err, i, j, ind, no_cgto_shells, pair_index, last_tgt_shell, n1_points, n2_points, &
                k, l, B_i_start, B_i_end, B_j_start, B_j_end, max_l
     integer :: B_start, B_end, max_l_cgto
     real(kind=cfp) :: bto_norm, t1, t2
     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied
     integer, parameter :: n_rng_knot_BB = 1 !Number of intervals in between each pair of knots within which the quadrature for [BB|GG] class will be applied
     logical, parameter :: only_on_bto_grid = .true.
     logical :: save_to_scratch

        write(stdout,'("--------->BG_mixed_2el_initialize")')

!UNCOMMENT TO DEBUG
!        call init_dbg(inp_bspline_grid)

        indexing_method = indexing_method_inp
        if (indexing_method > 2 .or. indexing_method <= 0) then
            call xermsg ('bto_gto_integrals_mod', 'BG_mixed_2el_initialize', &
                         'On input indexing_method was out of range [1,2].', 1, 1)
        end if

        no_cgto_shells = size(cgto_shells)
        max_l_legendre_2el = inp_max_l_legendre_2el
        write(stdout,'("Maximum L in the Leg. expansion of the Coulomb potential: ",i4)') max_l_legendre_2el

        if (max_l_legendre_2el .le. 0) then
           call xermsg('bto_gto_integrals_mod', 'BG_mixed_2el_initialize','On input max_l_legendre_2el was .le. 0.',2,1)
        endif

        max_l_BGGG = -1
        max_l_BGBG = -1

        overflow = F1MACH(2,cfp_dummy)
        prec_goal = F1MACH(4,cfp_dummy)*100

        if (use_extrapolation) write(stdout,'("Requested relative precision for the extrapolated integrals: ",e25.15)') prec_goal

        if (mixed_ints_method .eq. 2) then
           write(stdout,'(/,"Method for calculation of the BGGG class: Lebedev")')
        elseif (mixed_ints_method .eq. 1 .or. mixed_ints_method .eq. 3) then
           write(stdout,'(/,"Method for calculation of the BGGG class: Legendre expansion")')
        else
           call xermsg ('bto_gto_integrals_mod', 'BG_mixed_2el_initialize', &
                        'mixed_ints_method is out of range: allowed values are 1,2,3.', mixed_ints_method, 1)
        endif

        !max_bspline_l is assumed to be set to the largest angular momentum in the whole BTO basis.
        max_l_pw = inp_max_bspline_l+max_l_legendre_2el
        max_l_cgto = maxval(cgto_shells(:)%l)

        call init_CGTO_pw_expansions_mod(max_l_pw,max_l_cgto)

        n_gq = 2*n_7+1 !number of points of the G-L quadrature rule
        if (2*inp_bspline_grid%order > 2*n_gq-1) then
            call xermsg ('bto_gto_integrals_mod', 'BG_mixed_2el_initialize', &
                         'The B-spline order is too large for the n_10 Gauss-Legendre quadrature rule. &
                         &Decrease the B-spline order or increase the quadrature rule order.', 3, 1)
        end if

        if (no_cgto_shells <= 0) then
            call xermsg ('bto_gto_integrals_mod', 'BG_mixed_2el_initialize', 'On input no_cgto_shells was .le. 0.', 4, 1)
        end if

        keep_ab_cd_order = inp_keep_ab_cd_order
!
!------ Evaluate PW expansions of the all CGTO shells and all pairs of CGTO shells excluding the continuum shells.
!       WE ASSUME THAT ALL TARGET CGTO SHELLS PRECEDE THE CONTINUUM SHELLS. IF
!       THIS REQ. IS NOT MET THEN THE INTEGRALS WILL BE CALCULATED INCORRECTLY.
        last_tgt_shell = no_cgto_shells
        do i=1,no_cgto_shells
           if (cgto_shells(i)%non_zero_at_boundary) then
              last_tgt_shell = i-1
              exit
           endif
        enddo

        if (scratch_directory .eq. '') then
           save_to_scratch = .false.
           write(stdout,'("All Y_lm functions will be kept in memory.")')
        else
           save_to_scratch = .true.
           write(stdout,'("Saving to disk all Y_lm functions requiring more than ",e25.15," amount of memory (MiB).")') &
                Y_lm_size_threshold
        endif

        call grid_r1_r2%construct_r1_r2_grids(inp_bspline_grid,inp_first_bspline_index,inp_max_bspline_l,0,max_l_legendre_2el,&
                        nuclei,delta_r1,x_7,w_7,n_7,x_10,w_10,n_10,n_rng_knot)

        if (mixed_ints_method .eq. 2) then
           call grid_r1_r2%eval_Xlm_on_lebedev_grid(70)
        endif

        !Calculate the Y_l function for all unique pairs of radial B-splines
        call grid_r1_r2%eval_Y_l_BTO_BTO

        !Calculate the Y_lm function for each BTO/CGTO pair
        write(stdout,'(/,10X,"Number of CGTO shells: ",i4)') no_cgto_shells
        if (allocated(CGTO_pw)) deallocate(CGTO_pw)
        allocate(CGTO_pw(no_cgto_shells),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_mixed_2el_initialize','Memory allocation failed.',err,1)

        t1 = omp_get_wtime()
        !todo loop only over the symmetrically non-redundant shells!
        do i=1,no_cgto_shells
           call CGTO_pw(i)%init_CGTO_shell_pw_expansion(cgto_shells(i),i)

           !Evaluate the BTO/CGTO Y_lm function for all radial B-splines integrating over the r2 grid.
           if (mixed_ints_method .eq. 1 .or. mixed_ints_method .eq. 3) then !don't calculate Y_lm when lebedev is used for the BGGG class.
              write(stdout,'(/,"Evaluating pw expansion and the Y_lm function for CGTO shell: ",i4)') i
              call CGTO_pw(i)%eval_BTO_CGTO_Y_lm(grid_r1_r2)

              !Save the Y_lm disk only if requested and only those Y_lm requiring more than Y_lm_size_threshold MiBs of memory:
              if (save_to_scratch .and. CGTO_pw(i)%Y_lm_size_mib > Y_lm_size_threshold) then
                 call CGTO_pw(i)%write_Y_lm_to_file(scratch_directory)
                 deallocate(CGTO_pw(i)%Y_lm_mixed)
              endif
           endif

           !Evaluate the pw expansion on the r1 grid: this is needed for the BGBG exchange integrals
           call CGTO_pw(i)%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)

           if (mixed_ints_method .eq. 2) write(stdout,'("Evaluating pw expansion on the r1 grid for CGTO shell: ",i4)') i

           !Calculate the partial wave expansion of the CGTO on the r1 grid and on the grid of knots.
           call CGTO_pw(i)%eval_CGTO_shell_pw_expansion(grid_r1_r2%bspline_grid%knots,grid_r1_r2%max_bspline_l,&
                        grid_r1_r2%max_prop_l,grid_r1_r2%max_l_legendre)

           !Only when lebedev is used for the BGGG class
           !todo this can be evaluated directly in leb_BGGG class since it will be needed only once per each CGTO shell.
           if (mixed_ints_method .eq. 2) then
              call CGTO_pw(i)%eval_at_lebedev_points(grid_r1_r2)
              print *,'cgto eval at lebedev'
           endif
        enddo !i

        !Deallocate the arrays we no longer need for the integral evaluation:
        deallocate(grid_r1_r2%B_vals_r2,grid_r1_r2%w1_w2)

        t2 = omp_get_wtime()
        write(stdout,'("Calculation of the CGTO pw expansions took [s]: ",f25.15)') t2-t1

        write(stdout,'("<---------done:BG_mixed_2el_initialize")')

  end subroutine BG_mixed_2el_initialize

  subroutine BBGG_shell_integrals(bto_shell_A,bto_shell_B,cgto_shell_C,cgto_shell_D,A,B,C,D,starting_index_A,&
                                    starting_index_B,starting_index_C,starting_index_D,two_el_column,int_index,integrals)
     use gto_routines, only: index_2el, reorder_and_index_2el
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: eri
     implicit none
     type(BTO_shell_data_obj), intent(in) :: bto_shell_A,bto_shell_B
     type(CGTO_shell_data_obj), intent(in) :: cgto_shell_C,cgto_shell_D
     integer, intent(in) :: starting_index_A, starting_index_B, starting_index_C, starting_index_D, A,B,C,D, two_el_column
     !We assume that these three arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)

     integer :: n_cgto_DC_m, lm, m, ma, mb, m_DC, n_bto_BA_m, l, err, ind, n_integrals, pair_index, l_min, l_max, lm_start, &
                lm_end, p, q, bto_ind_A, bto_ind_B
     real(kind=cfp) :: R_AB(3), R_AB_square, K_pref, prod_alp, prod_P(3)
!UNCOMMENT TO DEBUG
!     integer :: lena,lenb,lenc,lend,ind_a,ind_b,ind_c,ind_d,la,lb,lc,ld,tgt_pair
!     real(kind=cfp) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
!     real(kind=cfp), allocatable :: anorms(:), bnorms(:), cnorms(:), dnorms(:)
!     real(kind=cfp), allocatable :: aexps(:), bexps(:), cexps(:), dexps(:)
!     real(kind=cfp), allocatable :: acoefs(:), bcoefs(:), ccoefs(:), dcoefs(:)
!     real(kind=cfp), allocatable :: tgt_prop(:,:), eri_int(:,:), eri_tail_int(:)
!     logical :: ab_is_continuum, do_tails_for_this_quartet
!     integer, allocatable :: int_index_eri(:,:)
!     character(len=4) :: tag = "BBGG"

        ! -- C, D are from the outer 2el loop index, each pair being processed by a different
        !    thread. This means that it is enough to make GG_pair_pw threadprivate, and wrap this section
        !    in omp critical to cater for the dirty use of global variables during precomputing.
        if (GG_pair_pw%cgto_shell_A_index .ne. C .or. GG_pair_pw%cgto_shell_B_index .ne. D) then
           !$OMP CRITICAL
           !todo in this case the pw expansion can only be up to 2*max_bspline_l: this should save a lot of time
           call GG_pair_pw%init_CGTO_shell_pair_pw_expansion(cgto_shell_C,C,cgto_shell_D,D)
           call GG_pair_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call GG_pair_pw%eval_CGTO_shell_pair_pw_expansion
           call GG_pair_pw%eval_radial_GG_BB(grid_r1_r2) !evaluate radial_lm_BB_GG
           write(stdout,'("BBGG_shell_integrals evaluated pw expansion for CGTO pair: ",2i4)') C,D
          !$OMP END CRITICAL
        endif

        n_cgto_DC_m = (2*cgto_shell_C%l+1)*(2*cgto_shell_D%l+1)
        n_bto_BA_m = (2*bto_shell_B%l+1)*(2*bto_shell_A%l+1)
        n_integrals = n_cgto_DC_m*n_bto_BA_m

        integrals(1:n_integrals,two_el_column) = 0.0_cfp

        l_min = abs(bto_shell_A%l-bto_shell_B%l)
        l_max = bto_shell_A%l+bto_shell_B%l

        lm_start = (l_min-1 +1)**2+1
        lm_end =   (l_max   +1)**2

        if (.not. allocated(couplings) .or. size(couplings) < lm_end) then
           if (allocated(couplings)) deallocate(couplings)
           allocate(couplings(lm_end),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BBGG_shell_integrals','Memory allocation failed.',err,1)
        endif

        bto_ind_A = bto_shell_A%bspline_index
        bto_ind_B = bto_shell_B%bspline_index
        !todo DBG
        !if (bto_shell_A%l .eq. 1) bto_ind_A = 10
        !if (bto_shell_B%l .eq. 1) bto_ind_B = 10

        pair_index = max(bto_ind_A,bto_ind_B)
        pair_index = pair_index*(pair_index-1)/2 + min(bto_ind_A,bto_ind_B)

        ind = 0
        do ma=-bto_shell_A%l,bto_shell_A%l
           do mb=-bto_shell_B%l,bto_shell_B%l
              do l=l_min,l_max
                 do m=-l,l
                    lm = l*l+l+m+1
                    couplings(lm) = cpl%rgaunt(l,bto_shell_A%l,bto_shell_B%l,m,ma,mb)
                 enddo
              enddo !l
              do m_DC=1,n_cgto_DC_m
                 integrals(ind+m_DC,two_el_column) = sum(couplings(lm_start:lm_end) &
                        * GG_pair_pw%radial_lm_BB_GG(lm_start:lm_end,m_DC,pair_index))
                 !print *,ind+m_DC,integrals(ind+m_DC,two_el_column)
              enddo !m_DC
              ind = ind + n_cgto_DC_m
           enddo !mb
        enddo !ma

        !Compute indices
        if (indexing_method .eq. 2) then
           call reorder_and_index_2el(cgto_shell_D%l,cgto_shell_C%l,bto_shell_B%l,bto_shell_A%l,starting_index_D,&
                starting_index_C,starting_index_B,starting_index_A,two_el_column,int_index,integrals)
        else
           call index_2el(cgto_shell_D%l,cgto_shell_C%l,bto_shell_B%l,bto_shell_A%l,starting_index_D,starting_index_C,&
                starting_index_B,starting_index_A,int_index,keep_ab_cd_order,.false.)
        endif

!UNCOMMENT TO DEBUG
!        allocate(anorms(1), bnorms(1), cnorms(cgto_shell_C%number_of_primitives), dnorms(cgto_shell_D%number_of_primitives))
!        allocate(aexps(1),   bexps(1),  cexps(cgto_shell_C%number_of_primitives),  dexps(cgto_shell_D%number_of_primitives))
!        allocate(acoefs(1), bcoefs(1), ccoefs(cgto_shell_C%number_of_primitives), dcoefs(cgto_shell_D%number_of_primitives),eri_int(n_integrals,two_el_column))
!        allocate(int_index_eri(4,n_integrals))
!        ab_is_continuum = .false.
!        do_tails_for_this_quartet = .false.
!        tgt_pair = 0
!        !Shell A
!        lena = 1
!        la = dbg_cgto(bto_ind_A)%l
!        xa = dbg_cgto(bto_ind_A)%center(1)
!        ya = dbg_cgto(bto_ind_A)%center(2)
!        za = dbg_cgto(bto_ind_A)%center(3)
!        anorms(1) = dbg_cgto(bto_ind_A)%norms(1)*dbg_cgto(bto_ind_A)%norm
!        acoefs(1) = dbg_cgto(bto_ind_A)%contractions(1)
!        aexps(1) = dbg_cgto(bto_ind_A)%exponents(1)
!        ind_a = starting_index_A
!        !Shell B
!        lenb = 1
!        lb = dbg_cgto(bto_ind_B)%l
!        xb = dbg_cgto(bto_ind_B)%center(1)
!        yb = dbg_cgto(bto_ind_B)%center(2)
!        zb = dbg_cgto(bto_ind_B)%center(3)
!        bnorms(1) = dbg_cgto(bto_ind_B)%norms(1)*dbg_cgto(bto_ind_B)%norm
!        bcoefs(1) = dbg_cgto(bto_ind_B)%contractions(1)
!        bexps(1) = dbg_cgto(bto_ind_B)%exponents(1)
!        ind_b = starting_index_B
!        !Shell C
!        lenc = cgto_shell_C%number_of_primitives
!        lc = cgto_shell_C%l
!        xc = cgto_shell_C%center(1)
!        yc = cgto_shell_C%center(2)
!        zc = cgto_shell_C%center(3)
!        cnorms(1:lenc) = cgto_shell_C%norms(1:lenc)*cgto_shell_C%norm
!        ccoefs(1:lenc) = cgto_shell_C%contractions(1:lenc)
!        cexps(1:lenc) = cgto_shell_C%exponents(1:lenc)
!        ind_c = starting_index_C
!        !Shell D
!        lend = cgto_shell_D%number_of_primitives
!        ld = cgto_shell_D%l
!        xd = cgto_shell_D%center(1)
!        yd = cgto_shell_D%center(2)
!        zd = cgto_shell_D%center(3)
!        dnorms(1:lend) = cgto_shell_D%norms(1:lend)*cgto_shell_D%norm
!        dcoefs(1:lend) = cgto_shell_D%contractions(1:lend)
!        dexps(1:lend) = cgto_shell_D%exponents(1:lend)
!        ind_d = starting_index_D
!
!        call eri(lena,xa,ya,za,anorms,la,aexps,acoefs,ind_a, lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,ind_b, lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,ind_c, lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,ind_d,&
!                & two_el_column,int_index_eri,keep_ab_cd_order, do_tails_for_this_quartet,ab_is_continuum,tgt_prop,tgt_pair,eri_tail_int,bto_shell_A%bspline_grid%B,eri_int)
!
!        call compare_print_2el_ints(tag,integrals,int_index,eri_int,int_index_eri,n_integrals,two_el_column)

  end subroutine BBGG_shell_integrals

  subroutine lebedev_BGGG_shell_integrals(bto_shell_A,cgto_shell_B,cgto_shell_C,cgto_shell_D,A,B,C,D,&
                                starting_index_A,starting_index_B,starting_index_C,starting_index_D,&
                                       &two_el_column,int_index,integrals)
     use gto_routines, only: index_2el, reorder_and_index_2el
     use general_quadrature, only: dqelg
     use phys_const, only: fourpi
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: eri
     implicit none
     type(BTO_shell_data_obj), intent(in) :: bto_shell_A
     type(CGTO_shell_data_obj), target, intent(in) :: cgto_shell_B,cgto_shell_C,cgto_shell_D
     integer, intent(in) :: starting_index_A, starting_index_B, starting_index_C, starting_index_D, A,B,C,D, two_el_column
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)

     integer :: n_integrals, max_l, i, k, n, ma,lp,mp
     integer :: l,m,ind,err,n_shell_A,n_shell_B,n_shell_C,n_shell_D,n_shell_DC,radial,n_angular
     integer :: m_DC, lm, min_l, n1_points, cgto_m, base
!UNCOMMENT TO DEBUG
!     integer :: lena,lenb,lenc,lend,ind_a,ind_b,ind_c,ind_d,la,lb,lc,ld,tgt_pair
!     real(kind=cfp) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
!     real(kind=cfp), allocatable :: anorms(:), bnorms(:), cnorms(:), dnorms(:)
!     real(kind=cfp), allocatable :: aexps(:), bexps(:), cexps(:), dexps(:)
!     real(kind=cfp), allocatable :: acoefs(:), bcoefs(:), ccoefs(:), dcoefs(:)
!     real(kind=cfp), allocatable :: tgt_prop(:,:), eri_int(:,:), eri_tail_int(:)
!     logical :: ab_is_continuum, do_tails_for_this_quartet
!     integer, allocatable :: int_index_eri(:,:)
!     character(len=4) :: tag = "BGGG"

        if (.not. CGTO_pw(B) % initialized) then
            call xermsg ('bto_gto_integrals_mod', 'lebedev_BGGG_shell_integrals', &
                         'The requested CGTO has not been initialized.', 1, 1)
        end if

        if (CGTO_pw(B) % cgto_shell_index /= B) then
            call xermsg ('bto_gto_integrals_mod', 'lebedev_BGGG_shell_integrals', &
                         'Shell index of the CGTO does not match with the required one.', 2, 1)
        end if

        !todo switch to an equivalent of BBGG when the Gb sits on CMS!!!
        if (GG_pair_pw%cgto_shell_A_index .ne. C .or. GG_pair_pw%cgto_shell_B_index .ne. D) then
           call GG_pair_pw%init_CGTO_shell_pair_pw_expansion(cgto_shell_C,C,cgto_shell_D,D)
           call GG_pair_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
!           call GG_pair_pw%eval_CGTO_shell_pair_pw_expansion
!           call GG_pair_pw%eval_radial_GG_BB(grid_r1_r2)
           call GG_pair_pw%eval_coulomb_integrals(grid_r1_r2)
           write(stdout,'("lebedev_BGGG_shell_integrals evaluated Coulomb integrals for CGTO pair: ",2i4)') C,D
        endif

        n_shell_A = 2*bto_shell_A%l+1
        n_shell_B = 2*cgto_shell_B%l+1
        n_shell_C = 2*cgto_shell_C%l+1
        n_shell_D = 2*cgto_shell_D%l+1
        n_shell_DC = n_shell_C*n_shell_D

        n_integrals = n_shell_A*n_shell_B*n_shell_C*n_shell_D
        integrals(1:n_integrals,two_el_column) = 0.0_cfp

        ind = 0
        n_angular = grid_r1_r2%lebedev_order
        do ma=-bto_shell_A%l,bto_shell_A%l
           lm = bto_shell_A%l*bto_shell_A%l+bto_shell_A%l+ma+1
           base = n_angular*(lm-1)
           do cgto_m=1,n_shell_B
              do m_DC=1,n_shell_DC
                 ind = ind + 1
                 do radial = grid_r1_r2 % bspline_start_end_r1(1, bto_shell_A % bspline_index), &
                             grid_r1_r2 % bspline_start_end_r1(2, bto_shell_A % bspline_index)
                    integrals(ind,two_el_column) = integrals(ind,two_el_column) + &
                          grid_r1_r2 % w1(radial) &
                        * grid_r1_r2 % r1(radial) &
                        * grid_r1_r2 % B_vals_r1(radial,bto_shell_A%bspline_index) &
                        * sum(grid_r1_r2 % Xlm_Lebedev(base+1:base+n_angular) &
                              * CGTO_pw(B) % at_lebedev_points(1:n_angular,radial,cgto_m) &
                              * GG_pair_pw % coulomb_integrals(1:n_angular,radial,m_DC))
                 enddo !radial
              enddo !m_DC
           enddo !mb
        enddo !ma

        if (indexing_method .eq. 2) then
           if (GG_pair_pw%order_AB) then
              call reorder_and_index_2el (cgto_shell_C%l,cgto_shell_D%l,cgto_shell_B%l,bto_shell_A%l,starting_index_C,&
                                            starting_index_D,starting_index_B,starting_index_A,two_el_column,int_index,integrals)
           else
              call reorder_and_index_2el (cgto_shell_D%l,cgto_shell_C%l,cgto_shell_B%l,bto_shell_A%l,starting_index_D,&
                                            starting_index_C,starting_index_B,starting_index_A,two_el_column,int_index,integrals)
           endif
        else
           !compute indices for order md,mc,mb,ma
           if (GG_pair_pw%order_AB) then
              call index_2el (cgto_shell_C%l,cgto_shell_D%l,cgto_shell_B%l,bto_shell_A%l,starting_index_C,starting_index_D, &
                                starting_index_B,starting_index_A,int_index,keep_ab_cd_order,.false.)
           else
              call index_2el(cgto_shell_D%l,cgto_shell_C%l,cgto_shell_B%l,bto_shell_A%l,starting_index_D,starting_index_C, &
                                starting_index_B,starting_index_A,int_index,keep_ab_cd_order,.false.)
           endif
        endif

!UNCOMMENT TO DEBUG
!        allocate(anorms(1), bnorms(cgto_shell_B%number_of_primitives), cnorms(cgto_shell_C%number_of_primitives), dnorms(cgto_shell_D%number_of_primitives))
!        allocate(aexps(1),   bexps(cgto_shell_B%number_of_primitives),  cexps(cgto_shell_C%number_of_primitives),  dexps(cgto_shell_D%number_of_primitives))
!        allocate(acoefs(1), bcoefs(cgto_shell_B%number_of_primitives), ccoefs(cgto_shell_C%number_of_primitives), dcoefs(cgto_shell_D%number_of_primitives),eri_int(n_integrals,two_el_column))
!        allocate(int_index_eri(4,n_integrals))
!        ab_is_continuum = .false.
!        do_tails_for_this_quartet = .false.
!        tgt_pair = 0
!        !Shell A
!        lena = 1
!        la = dbg_cgto(bto_shell_A%bspline_index)%l
!        xa = dbg_cgto(bto_shell_A%bspline_index)%center(1)
!        ya = dbg_cgto(bto_shell_A%bspline_index)%center(2)
!        za = dbg_cgto(bto_shell_A%bspline_index)%center(3)
!        anorms(1) = dbg_cgto(bto_shell_A%bspline_index)%norms(1)*dbg_cgto(bto_shell_A%bspline_index)%norm
!        acoefs(1) = dbg_cgto(bto_shell_A%bspline_index)%contractions(1)
!        aexps(1) = dbg_cgto(bto_shell_A%bspline_index)%exponents(1)
!        ind_a = starting_index_A
!        !Shell B
!        lenb = cgto_shell_B%number_of_primitives
!        lb = cgto_shell_B%l
!        xb = cgto_shell_B%center(1)
!        yb = cgto_shell_B%center(2)
!        zb = cgto_shell_B%center(3)
!        bnorms(1:lenb) = cgto_shell_B%norms(1:lenb)*cgto_shell_B%norm
!        bcoefs(1:lenb) = cgto_shell_B%contractions(1:lenb)
!        bexps(1:lenb) = cgto_shell_B%exponents(1:lenb)
!        ind_b = starting_index_B
!        !Shell C
!        lenc = cgto_shell_C%number_of_primitives
!        lc = cgto_shell_C%l
!        xc = cgto_shell_C%center(1)
!        yc = cgto_shell_C%center(2)
!        zc = cgto_shell_C%center(3)
!        cnorms(1:lenc) = cgto_shell_C%norms(1:lenc)*cgto_shell_C%norm
!        ccoefs(1:lenc) = cgto_shell_C%contractions(1:lenc)
!        cexps(1:lenc) = cgto_shell_C%exponents(1:lenc)
!        ind_c = starting_index_C
!        !Shell D
!        lend = cgto_shell_D%number_of_primitives
!        ld = cgto_shell_D%l
!        xd = cgto_shell_D%center(1)
!        yd = cgto_shell_D%center(2)
!        zd = cgto_shell_D%center(3)
!        dnorms(1:lend) = cgto_shell_D%norms(1:lend)*cgto_shell_D%norm
!        dcoefs(1:lend) = cgto_shell_D%contractions(1:lend)
!        dexps(1:lend) = cgto_shell_D%exponents(1:lend)
!        ind_d = starting_index_D
!
!        call eri(lena,xa,ya,za,anorms,la,aexps,acoefs,ind_a, lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,ind_b, lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,ind_c, lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,ind_d,&
!                & two_el_column,int_index_eri,keep_ab_cd_order, do_tails_for_this_quartet,ab_is_continuum,tgt_prop,tgt_pair,eri_tail_int,bto_shell_A%bspline_grid%B,eri_int)
!
!        print *,'L=',lb,lc,ld
!        print *,'alp b',bexps
!        print *,'alp c',cexps
!        print *,'alp d',dexps
!        call compare_print_2el_ints(tag,integrals,int_index,eri_int,int_index_eri,n_integrals,two_el_column)

  end subroutine lebedev_BGGG_shell_integrals

  !> BG pair is r2, GG pair is r1.
  subroutine BGGG_shell_integrals(bto_shell_A,cgto_shell_B,cgto_shell_C,cgto_shell_D,A,B,C,D,starting_index_A,&
                    starting_index_B,starting_index_C,starting_index_D,two_el_column,int_index,integrals)
     use gto_routines, only: index_2el, reorder_and_index_2el
     use general_quadrature, only: dqelg
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: eri
     implicit none
     type(BTO_shell_data_obj), intent(in) :: bto_shell_A
     type(CGTO_shell_data_obj), target, intent(in) :: cgto_shell_B,cgto_shell_C,cgto_shell_D
     integer, intent(in) :: starting_index_A, starting_index_B, starting_index_C, starting_index_D, A,B,C,D, two_el_column
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)

     integer :: n_integrals, max_l, i, k, n, n_wynn, ma,lp,mp, l,m,ind,err,nres, &
                n_shell_A,n_shell_B,n_shell_C,n_shell_D,n_shell_DC,lama,p,bto_ind_A,l1,l2,dl
     integer :: m_DC, lm, not_converged, min_l, n1_points, cgto_m
     real(kind=cfp) :: res, abserr, rel_precision_1, rel_precision_2, cf
     logical :: extrapolate
!UNCOMMENT TO DEBUG
!     integer :: lena,lenb,lenc,lend,ind_a,ind_b,ind_c,ind_d,la,lb,lc,ld,tgt_pair
!     real(kind=cfp) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
!     real(kind=cfp), allocatable :: anorms(:), bnorms(:), cnorms(:), dnorms(:)
!     real(kind=cfp), allocatable :: aexps(:), bexps(:), cexps(:), dexps(:)
!     real(kind=cfp), allocatable :: acoefs(:), bcoefs(:), ccoefs(:), dcoefs(:)
!     real(kind=cfp), allocatable :: tgt_prop(:,:), eri_int(:,:), eri_tail_int(:)
!     logical :: ab_is_continuum, do_tails_for_this_quartet
!     integer, allocatable :: int_index_eri(:,:)
!     character(len=4) :: tag = "BGGG"

        if (.not. CGTO_pw(B) % initialized) then
            call xermsg ('bto_gto_integrals_mod', 'BGGG_shell_integrals', &
                         'The requested CGTO has not been initialized.', 1, 1)
        end if

        if (CGTO_pw(B) % cgto_shell_index /= B) then
            call xermsg ('bto_gto_integrals_mod', 'BGGG_shell_integrals', &
                         'Shell index of the CGTO does not match with the required one.', 2, 1)
        end if

        ! -- C, D are from the outer 2el loop index, each pair being processed by a different
        !    thread. This means that it is enough to make GG_pair_pw threadprivate, and wrap this section
        !    in omp critical to cater for the dirty use of global variables during precomputing.
        if (GG_pair_pw%cgto_shell_A_index .ne. C .or. GG_pair_pw%cgto_shell_B_index .ne. D) then
           !$OMP CRITICAL
           call GG_pair_pw%init_CGTO_shell_pair_pw_expansion(cgto_shell_C,C,cgto_shell_D,D)
           call GG_pair_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call GG_pair_pw%eval_CGTO_shell_pair_pw_expansion
           call GG_pair_pw%eval_radial_GG_BB(grid_r1_r2)
           write(stdout,'("BGGG_shell_integrals evaluated pw expansion for CGTO pair: ",2i4)') C,D
           !$OMP END CRITICAL
        endif

        n_shell_A = 2*bto_shell_A%l+1
        n_shell_B = 2*cgto_shell_B%l+1
        n_shell_C = 2*cgto_shell_C%l+1
        n_shell_D = 2*cgto_shell_D%l+1
        n_shell_DC = n_shell_C*n_shell_D

        n_integrals = n_shell_A*n_shell_B*n_shell_C*n_shell_D
        integrals(1:n_integrals,two_el_column) = 0.0_cfp

        bto_ind_A = bto_shell_A%bspline_index
        !todo DBG
        !if (bto_shell_A%l .eq. 1) bto_ind_A = 10

        max_l = max_l_legendre_2el
        min_l = 0
        extrapolate = use_extrapolation

        !No extrapolation must be used if the CGTO B or CGTO C and CGTO D lie on CMS since then the Legendre expansion is finite (exact).
        if (dot_product(cgto_shell_B%center,cgto_shell_B%center) .le. epsabs) then
           max_l = min(max_l,cgto_shell_B%l+bto_shell_A%l)
           min_l = max(min_l,abs(cgto_shell_B%l-bto_shell_A%l))
           extrapolate = .false.
        endif

        if (dot_product(cgto_shell_C % center, cgto_shell_C % center) <= epsabs .and. &
            dot_product(cgto_shell_D % center, cgto_shell_D % center) <= epsabs) then
           max_l = min(max_l,cgto_shell_C%l+cgto_shell_D%l)
           min_l = max(min_l,abs(cgto_shell_C%l-cgto_shell_D%l))
           extrapolate = .false.
        endif

        if (CGTO_pw(B)%Y_lm_on_disk) then
           n1_points = CGTO_pw(B)%Y_lm_dim_on_disk(1)

           if (.not. allocated(Y_lm_mixed_from_disk) .or. size(Y_lm_mixed_from_disk) < n1_points) then
              if (allocated(Y_lm_mixed_from_disk)) deallocate(Y_lm_mixed_from_disk)
              allocate(Y_lm_mixed_from_disk(n1_points),stat=err)
              if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BGGG_shell_integrals','Memory allocation 0 failed.',err,1)
           endif
        else
           n1_points = size(CGTO_pw(B)%Y_lm_mixed,1)
        endif

        if (.not. allocated(done) .or. size(done) < n_integrals) then
           if (allocated(done)) deallocate(done)
           allocate(done(n_integrals),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_nai_integrals','Memory allocation 1 failed.',err,1)
        endif
        done = .false.

        if (extrapolate .and. (.not. allocated(epstab) .or. .not. allocated(res3la) .or. size(epstab,2) < n_integrals)) then
           if (allocated(epstab)) deallocate(epstab)
           if (allocated(res3la)) deallocate(res3la)
           allocate(epstab(max_epstab,n_integrals),res3la(3,n_integrals),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BGGG_shell_integrals','Memory allocation 2 failed.',err,1)
        endif

        if (extrapolate) then
           max_l = min(max_l,max_epstab-3)
           if (max_l < 6) stop "max_l must be at least 6 due to extrapolation"
           epstab(1:max_epstab,1:n_integrals) = 0.0_cfp
           res3la = 0.0_cfp
           l1 = min_l
           l2 = max_l
           dl = 1
        else
           l1 = max_l
           l2 = min_l
           dl = -1
        endif

        n = 0
        do l=l1,l2,dl

           n = n+1

           do m=-l,l
         
              lm = l*l+l+m+1

              ind = 0
              do ma=-bto_shell_A%l,bto_shell_A%l
                 lama = bto_shell_A%l*bto_shell_A%l+bto_shell_A%l+ma+1
                 do cgto_m=1,n_shell_B
                    p = CGTO_pw(B)%Y_lm_non_neg_indices(cgto_m,lama,lm)
                    if (p .eq. 0) then
                       ind = ind + n_shell_DC
                       cycle
                    endif
                    if (CGTO_pw(B)%Y_lm_on_disk) then
                       !$OMP CRITICAL
                       call CGTO_pw(B)%read_Y_lm_from_file(p,bto_ind_A,Y_lm_mixed_from_disk)
                       !$OMP END CRITICAL

                       do m_DC=1,n_shell_DC
                          if (GG_pair_pw%neglect_m_lm(m_DC,lm) .or. done(ind+m_DC)) cycle
                          !save in order D,C,B,A
                          integrals(ind+m_DC,two_el_column) = integrals(ind+m_DC,two_el_column) &
                            + sum(GG_pair_pw%angular_integrals(1:n1_points,m_DC,lm)*grid_r1_r2%r1(1:n1_points) &
                                * Y_lm_mixed_from_disk(1:n1_points))
                       enddo !m_DC
                    else !Y_lm in memory
                       do m_DC=1,n_shell_DC
                          if (GG_pair_pw%neglect_m_lm(m_DC,lm) .or. done(ind+m_DC)) cycle
                          !save in order D,C,B,A
                          integrals(ind+m_DC,two_el_column) = integrals(ind+m_DC,two_el_column) &
                            + sum(GG_pair_pw%angular_integrals(1:n1_points,m_DC,lm)*grid_r1_r2%r1(1:n1_points) &
                                * CGTO_pw(B)%Y_lm_mixed(1:n1_points,p,bto_ind_A))
                       enddo !m_DC
                    endif
                    ind = ind + n_shell_DC
                 enddo !cgto_mb
              enddo !ma

           enddo !m

           if (extrapolate) then
              !The value of n_wynn can get changed by dqelg if numerical difficulties are encountered: 
              !typically this will happen for very small (negligible) values of the final integrals.
              rel_precision_1 = 0.0_cfp; rel_precision_2 = 0.0_cfp
              do i=1,n_integrals
                 if (done(i)) cycle
                 epstab(n,i) = integrals(i,two_el_column)
                 if (n .ge. 3) then
                    n_wynn = n
                    nres = n-3
                    call dqelg(n_wynn, epstab(1:max_epstab,i), res, abserr, res3la(1:3,i), nres)
                    if (n .ge. 6) then !the first 3 estimates of the absolute error are bogus
                       if (res .ne. 0.0_cfp) then
                          !Mark as done all integrals with projected values below the threshold
                          if (abs(res) < epsabs) then
                             done(i) = .true.
                          else
                             rel_precision_2 = abs((res3la(2,i)-res3la(3,i))/res)
                             rel_precision_1 = abs((res3la(1,i)-res3la(3,i))/res)
                             !It is important to terminate the extrapolation only
                             !if last two extrapolations are converged since in
                             !the case the pair of CGTOs are identical and lie on
                             !centers obtained by inversion through CMS then for
                             !some l-values there will be no contribution so the
                             !difference wrt previous extrapolation would be 0 but
                             !that does not have to mean convergence since the next
                             !l-value may contribute significantly.
                             if (rel_precision_1 .le. prec_goal .and. rel_precision_2 .le. prec_goal) done(i) = .true.
                          endif
                       else
                          done(i) = .true.
                       endif
                    endif
                    !todo temp
                    !if (i .eq. 1) print *,n,integrals(i,two_el_column),res
                 endif
              enddo !i
              
              !Terminate the series when all integrals have been converged to the required precision. 
              not_converged = n_integrals-count(done)
              if (not_converged .eq. 0) exit

           endif

        enddo !l

        max_l_BGGG = max(max_l_BGGG,l)

        if (extrapolate) then
           if (not_converged > 0) then
              print *,'not converged BGGG',not_converged
              print *,A,B,C,D
              do i=1,n_integrals
                 if (.not.(done(i))) then
                    rel_precision_2 = abs((res3la(2,i)-res3la(3,i))/res3la(3,i))
                    rel_precision_1 = abs((res3la(1,i)-res3la(3,i))/res3la(3,i))
                    write(*,'(i0,3e25.15)') i,res3la(3,i),rel_precision_2,rel_precision_1 !show the last 3 extrapolated value to see what the relative precision is.
                 endif
              enddo
              call xermsg ('bto_gto_integrals_mod', 'BGGG_shell_integrals', &
                           'The required relative precision has not been reached for all integrals in the quartet of shells. &
                           &Try increasing max_l_legendre_2el.', 1, 0)
           endif
           !Replace the integrals obtained by direct summation of the Legendre series with the extrapolated values.
           do i=1,n_integrals
              integrals(i,two_el_column) = res3la(3,i)
              !print *,i,integrals(i,two_el_column)
           enddo
        endif

        !Compute indices
        if (indexing_method .eq. 2) then
           call reorder_and_index_2el (cgto_shell_D%l,cgto_shell_C%l,cgto_shell_B%l,bto_shell_A%l,starting_index_D,&
                                        starting_index_C,starting_index_B,starting_index_A,two_el_column,int_index,integrals)
        else
           call index_2el (cgto_shell_D%l,cgto_shell_C%l,cgto_shell_B%l,bto_shell_A%l,starting_index_D,starting_index_C,&
                            starting_index_B,starting_index_A,int_index,keep_ab_cd_order,.false.)
        endif

!        print *,'chek',integrals(1:n_integrals,two_el_column)

!UNCOMMENT TO DEBUG
!        allocate(anorms(1), bnorms(cgto_shell_B%number_of_primitives), cnorms(cgto_shell_C%number_of_primitives), dnorms(cgto_shell_D%number_of_primitives))
!        allocate(aexps(1),   bexps(cgto_shell_B%number_of_primitives),  cexps(cgto_shell_C%number_of_primitives),  dexps(cgto_shell_D%number_of_primitives))
!        allocate(acoefs(1), bcoefs(cgto_shell_B%number_of_primitives), ccoefs(cgto_shell_C%number_of_primitives), dcoefs(cgto_shell_D%number_of_primitives),eri_int(n_integrals,two_el_column))
!        allocate(int_index_eri(4,n_integrals))
!        ab_is_continuum = .false.
!        do_tails_for_this_quartet = .false.
!        tgt_pair = 0
!        !Shell A
!        lena = 1
!        la = dbg_cgto(bto_ind_A)%l
!        xa = dbg_cgto(bto_ind_A)%center(1)
!        ya = dbg_cgto(bto_ind_A)%center(2)
!        za = dbg_cgto(bto_ind_A)%center(3)
!        anorms(1) = dbg_cgto(bto_ind_A)%norms(1)*dbg_cgto(bto_ind_A)%norm
!        acoefs(1) = dbg_cgto(bto_ind_A)%contractions(1)
!        aexps(1) = dbg_cgto(bto_ind_A)%exponents(1)
!        ind_a = starting_index_A
!        !Shell B
!        lenb = cgto_shell_B%number_of_primitives
!        lb = cgto_shell_B%l
!        xb = cgto_shell_B%center(1)
!        yb = cgto_shell_B%center(2)
!        zb = cgto_shell_B%center(3)
!        bnorms(1:lenb) = cgto_shell_B%norms(1:lenb)*cgto_shell_B%norm
!        bcoefs(1:lenb) = cgto_shell_B%contractions(1:lenb)
!        bexps(1:lenb) = cgto_shell_B%exponents(1:lenb)
!        ind_b = starting_index_B
!        !Shell C
!        lenc = cgto_shell_C%number_of_primitives
!        lc = cgto_shell_C%l
!        xc = cgto_shell_C%center(1)
!        yc = cgto_shell_C%center(2)
!        zc = cgto_shell_C%center(3)
!        cnorms(1:lenc) = cgto_shell_C%norms(1:lenc)*cgto_shell_C%norm
!        ccoefs(1:lenc) = cgto_shell_C%contractions(1:lenc)
!        cexps(1:lenc) = cgto_shell_C%exponents(1:lenc)
!        ind_c = starting_index_C
!        !Shell D
!        lend = cgto_shell_D%number_of_primitives
!        ld = cgto_shell_D%l
!        xd = cgto_shell_D%center(1)
!        yd = cgto_shell_D%center(2)
!        zd = cgto_shell_D%center(3)
!        dnorms(1:lend) = cgto_shell_D%norms(1:lend)*cgto_shell_D%norm
!        dcoefs(1:lend) = cgto_shell_D%contractions(1:lend)
!        dexps(1:lend) = cgto_shell_D%exponents(1:lend)
!        ind_d = starting_index_D
!
!        call eri(lena,xa,ya,za,anorms,la,aexps,acoefs,ind_a, lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,ind_b, lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,ind_c, lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,ind_d,&
!                & two_el_column,int_index_eri,keep_ab_cd_order, do_tails_for_this_quartet,ab_is_continuum,tgt_prop,tgt_pair,eri_tail_int,bto_shell_A%bspline_grid%B,eri_int)
!
!        call compare_print_2el_ints(tag,integrals,int_index,eri_int,int_index_eri,n_integrals,two_el_column)

  end subroutine BGGG_shell_integrals

  !> CD pair is r2, AB pair is r1.
  subroutine BGBG_shell_integrals(bto_shell_A,cgto_shell_B,bto_shell_C,cgto_shell_D,A,B,C,D,starting_index_A,&
                                starting_index_B,starting_index_C,starting_index_D,two_el_column,int_index,integrals)
     use gto_routines, only: index_2el, reorder_and_index_2el
!UNCOMMENT TO DEBUG
!     use cgto_hgp, only: eri
     use phys_const, only: fourpi
     use general_quadrature, only: dqelg
     implicit none
     type(BTO_shell_data_obj), intent(in) :: bto_shell_A,bto_shell_C
     type(CGTO_shell_data_obj), intent(in) :: cgto_shell_B,cgto_shell_D
     integer, intent(in) :: starting_index_A, starting_index_B, starting_index_C, starting_index_D, A,B,C,D, two_el_column
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)

     integer :: n_integrals, max_l, i, k, n, n_wynn, ma,lp,mp, l,m,ind,err,nres,n_shell_A,n_shell_B,n_shell_C,n_shell_D,&
                n_shell_DC,lpmp,ind_ma,lama,p,q,bto_ind_A,bto_ind_C,l1,l2,dl
     integer :: mc, md, lcmc, lm, not_converged, min_l, n1_points, cgto_m, Bc_start, Bc_end
     real(kind=cfp) :: res, abserr, rel_precision_1, rel_precision_2, cf
     logical :: extrapolate
!UNCOMMENT TO DEBUG
!     integer :: lena,lenb,lenc,lend,ind_a,ind_b,ind_c,ind_d,la,lb,lc,ld,tgt_pair
!     real(kind=cfp) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
!     real(kind=cfp), allocatable :: anorms(:), bnorms(:), cnorms(:), dnorms(:)
!     real(kind=cfp), allocatable :: aexps(:), bexps(:), cexps(:), dexps(:)
!     real(kind=cfp), allocatable :: acoefs(:), bcoefs(:), ccoefs(:), dcoefs(:)
!     real(kind=cfp), allocatable :: tgt_prop(:,:), eri_int(:,:), eri_tail_int(:)
!     logical :: ab_is_continuum, do_tails_for_this_quartet
!     integer, allocatable :: int_index_eri(:,:)
!     character(len=4) :: tag = "BGBG"

        if (CGTO_pw(B) % cgto_shell_index /= B) then
            call xermsg ('bto_gto_integrals_mod', 'BGBG_shell_integrals', &
                         'Shell index of the CGTO B does not match with the required one.', 1, 1)
        end if
        if (CGTO_pw(D) % cgto_shell_index /= D) then
            call xermsg ('bto_gto_integrals_mod', 'BGBG_shell_integrals', &
                         'Shell index of the CGTO D does not match with the required one.', 2, 1)
        end if

        ! -- The first thread that needs a particular CGTO_pw(B)%Y_lm_mixed will precompute it.
        !    Other threads will then reuse the data.
        !$OMP CRITICAL
        if (size(CGTO_pw(B)%Y_lm_mixed,1) .eq. 0) then
           write(stdout,'("Evaluating pw expansion and the Y_lm function for CGTO shell: ",i4)') B
           call CGTO_pw(B)%init_CGTO_shell_pw_expansion(cgto_shell_B,B)

           !Evaluate the BTO/CGTO Y_lm function for all radial B-splines integrating over the r2 grid.
           call CGTO_pw(B)%eval_BTO_CGTO_Y_lm(grid_r1_r2)

           !Calculate the partial wave expansion of the CGTO on the r1 grid
           call CGTO_pw(B)%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call CGTO_pw(B)%eval_CGTO_shell_pw_expansion(grid_r1_r2%bspline_grid%knots,&
                                grid_r1_r2%max_bspline_l,grid_r1_r2%max_prop_l,grid_r1_r2%max_l_legendre)
        endif
        !$OMP END CRITICAL

        if (CGTO_pw(B)%Y_lm_on_disk) then
           n1_points = CGTO_pw(B)%Y_lm_dim_on_disk(1)

           if (.not. allocated(Y_lm_mixed_from_disk) .or. size(Y_lm_mixed_from_disk) < n1_points) then
              if (allocated(Y_lm_mixed_from_disk)) deallocate(Y_lm_mixed_from_disk)
              allocate(Y_lm_mixed_from_disk(n1_points),stat=err)
              if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BGBG_shell_integrals','Memory allocation 0 failed.',err,1)
           endif
        else
           n1_points = size(CGTO_pw(B)%Y_lm_mixed,1)
        endif

        if (n1_points .ne. size(CGTO_pw(D)%gaunt_angular_integrals,1)) then
           print *,n1_points,size(CGTO_pw(D)%gaunt_angular_integrals,1)
           call xermsg('bto_gto_integrals_mod','BGBG_shell_integrals','CGTO B and D pw expansion are incompatible.',3,1)
        endif

        n_shell_A = 2*bto_shell_A%l+1
        n_shell_B = 2*cgto_shell_B%l+1
        n_shell_C = 2*bto_shell_C%l+1
        n_shell_D = 2*cgto_shell_D%l+1

        n_integrals = n_shell_A*n_shell_B*n_shell_C*n_shell_D
        integrals(1:n_integrals,two_el_column) = 0.0_cfp

        bto_ind_A = bto_shell_A%bspline_index
        bto_ind_C = bto_shell_C%bspline_index
        !todo DBG
        !if (bto_shell_A%l .eq. 1) bto_ind_A = 10
        !if (bto_shell_C%l .eq. 1) bto_ind_C = 10

        max_l = max_l_legendre_2el
        min_l = 0
        extrapolate = use_extrapolation

        !No extrapolation must be used if one of the CGTO shells lies on CMS since then the Legendre expansion is finite (exact).
        if (dot_product(cgto_shell_B%center,cgto_shell_B%center) .le. epsabs) then
           max_l = min(max_l,cgto_shell_B%l+bto_shell_A%l)
           min_l = max(min_l,abs(cgto_shell_B%l-bto_shell_A%l))
           extrapolate = .false.
        endif

        if (dot_product(cgto_shell_D%center,cgto_shell_D%center) .le. epsabs) then
           max_l = min(max_l,cgto_shell_D%l+bto_shell_C%l)
           min_l = max(min_l,abs(cgto_shell_D%l-bto_shell_C%l))
           extrapolate = .false.
        endif

        Bc_start = grid_r1_r2%bspline_start_end_r1(1,bto_ind_C)
        Bc_end   = grid_r1_r2%bspline_start_end_r1(2,bto_ind_C)

        if (.not. allocated(done) .or. size(done) < n_integrals) then
           if (allocated(done)) deallocate(done)
           allocate(done(n_integrals),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BG_nai_integrals','Memory allocation 1 failed.',err,1)
        endif
        done = .false.

        if (extrapolate .and. (.not. allocated(epstab) .or. .not. allocated(res3la) .or. size(epstab,2) < n_integrals)) then
           if (allocated(epstab)) deallocate(epstab)
           if (allocated(res3la)) deallocate(res3la)
           allocate(epstab(max_epstab,n_integrals),res3la(3,n_integrals),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BGBG_shell_integrals','Memory allocation 2 failed.',err,1)
        endif

        if (extrapolate) then
           max_l = min(max_l,max_epstab-3)
           if (max_l < 6) stop "max_l must be at least 6 due to extrapolation"
           epstab(1:max_epstab,1:n_integrals) = 0.0_cfp
           res3la = 0.0_cfp
        endif

        if (extrapolate) then
           max_l = min(max_l,max_epstab-3)
           if (max_l < 6) stop "max_l must be at least 6 due to extrapolation"
           epstab(1:max_epstab,1:n_integrals) = 0.0_cfp
           res3la = 0.0_cfp
           l1 = min_l
           l2 = max_l
           dl = 1
        else
           l1 = max_l
           l2 = min_l
           dl = -1
        endif

        n = 0
        do l=l1,l2,dl

           n = n+1

           do m=-l,l
         
              lm = l*l+l+m+1

              ind = 0
              do ma=-bto_shell_A%l,bto_shell_A%l
                 lama = bto_shell_A%l*bto_shell_A%l+bto_shell_A%l+ma+1
                 do cgto_m=1,n_shell_B
                    p = CGTO_pw(B)%Y_lm_non_neg_indices(cgto_m,lama,lm)
                    if (p .eq. 0) then
                       ind = ind + n_shell_D*n_shell_C
                       cycle
                    endif
                    if (CGTO_pw(B)%Y_lm_on_disk) then
                       !todo only the short section Bc_start:Bc_end has to be read!
                       !$OMP CRITICAL
                       call CGTO_pw(B)%read_Y_lm_from_file(p,bto_ind_A,Y_lm_mixed_from_disk)
                       !$OMP END CRITICAL

                       do mc=-bto_shell_C%l,bto_shell_C%l
                          lcmc = bto_shell_C%l*bto_shell_C%l+bto_shell_C%l+mc+1
                          do md=1,n_shell_D
                             q = CGTO_pw(D)%non_neg_indices_l_lp(md,lcmc,lm)
                             if (q .eq. 0 .or. done(ind+md)) cycle
                             !save in order D,C,B,A
                             integrals(ind+md,two_el_column) = integrals(ind+md,two_el_column) &
                                + sum(Y_lm_mixed_from_disk(Bc_start:Bc_end)&
                                    *grid_r1_r2%B_vals_r1(Bc_start:Bc_end,bto_ind_C)&
                                    *CGTO_pw(D)%gaunt_angular_integrals(Bc_start:Bc_end,q))
                          enddo !md
                          ind = ind + n_shell_D
                       enddo !mc
                    else !Y_lm in memory
                       do mc=-bto_shell_C%l,bto_shell_C%l
                          lcmc = bto_shell_C%l*bto_shell_C%l+bto_shell_C%l+mc+1
                          do md=1,n_shell_D
                             q = CGTO_pw(D)%non_neg_indices_l_lp(md,lcmc,lm)
                             if (q .eq. 0 .or. done(ind+md)) cycle
                             !save in order D,C,B,A
                             integrals(ind+md,two_el_column) = integrals(ind+md,two_el_column) &
                                + sum(CGTO_pw(B)%Y_lm_mixed(Bc_start:Bc_end,p,bto_ind_A)&
                                    *grid_r1_r2%B_vals_r1(Bc_start:Bc_end,bto_ind_C)&
                                    *CGTO_pw(D)%gaunt_angular_integrals(Bc_start:Bc_end,q))
                          enddo !md
                          ind = ind + n_shell_D
                       enddo !mc
                    endif
                 enddo !cgto_mb

              enddo !ma

           enddo !m

           if (extrapolate) then
              !The value of n_wynn can get changed by dqelg if numerical difficulties are encountered: 
              !typically this will happen for very small (negligible) values of the final integrals.
              rel_precision_1 = 0.0_cfp; rel_precision_2 = 0.0_cfp
              do i=1,n_integrals
                 if (done(i)) cycle
                 epstab(n,i) = integrals(i,two_el_column)
                 if (n .ge. 3) then
                    n_wynn = n
                    nres = n-3
                    call dqelg(n_wynn, epstab(1:max_epstab,i), res, abserr, res3la(1:3,i), nres)
                    if (n .ge. 6) then !the first 3 estimates of the absolute error are bogus
                       if (res .ne. 0.0_cfp) then
                          !Mark as done all integrals with projected values below the threshold
                          if (abs(res) < epsabs) then
                             done(i) = .true.
                          else
                             rel_precision_2 = abs((res3la(2,i)-res3la(3,i))/res)
                             rel_precision_1 = abs((res3la(1,i)-res3la(3,i))/res)
                             !It is important to terminate the extrapolation only
                             !if last two extrapolations are converged since in
                             !the case the pair of CGTOs are identical and lie on
                             !centers obtained by inversion through CMS then for
                             !some l-values there will be no contribution so the
                             !difference wrt previous extrapolation would be 0 but
                             !that does not have to mean convergence since the next
                             !l-value may contribute significantly.
                             if (rel_precision_1 .le. prec_goal .and. rel_precision_2 .le. prec_goal) done(i) = .true.
                          endif
                       else
                          done(i) = .true.
                       endif
                    endif
                    !todo temp
                    !if (i .eq. 1) print *,n,integrals(i,two_el_column),res
                 endif
              enddo !i
              
              !Terminate the series when all integrals have been converged to the required precision. 
              not_converged = n_integrals-count(done)
              if (not_converged .eq. 0) exit

           endif

        enddo !l

        max_l_BGBG = max(max_l_BGBG,l)

        if (extrapolate) then
           if (not_converged > 0) then
              print *,'not converged BGBG',not_converged
              print *,A,B,C,D
              do i=1,n_integrals
                 if (.not.(done(i))) then
                    rel_precision_2 = abs((res3la(2,i)-res3la(3,i))/res3la(3,i))
                    rel_precision_1 = abs((res3la(1,i)-res3la(3,i))/res3la(3,i))
                    write(*,'(i0,3e25.15)') i,res3la(3,i),rel_precision_2,rel_precision_1 !show the last 3 extrapolated value to see what the relative precision is.
                 endif
              enddo
              call xermsg ('bto_gto_integrals_mod', 'BGBG_shell_integrals', &
                           'The required relative precision has not been reached for all integrals in the quartet of shells. &
                           &Try increasing max_l_legendre_2el.', 1, 0)
           endif
           !Replace the integrals obtained by direct summation of the Legendre series with the extrapolated values.
           do i=1,n_integrals
              integrals(i,two_el_column) = res3la(3,i)
              !print *,i,integrals(i,two_el_column)
           enddo
        endif

        !Compute indices
        if (indexing_method .eq. 2) then
           call reorder_and_index_2el (cgto_shell_D%l,bto_shell_C%l,cgto_shell_B%l,bto_shell_A%l,starting_index_D,&
                                    starting_index_C,starting_index_B,starting_index_A,two_el_column,int_index,integrals)
        else
           call index_2el (cgto_shell_D%l,bto_shell_C%l,cgto_shell_B%l,bto_shell_A%l,starting_index_D,starting_index_C,&
                            starting_index_B,starting_index_A,int_index,keep_ab_cd_order,.false.)
        endif

!UNCOMMENT TO DEBUG
!        allocate(anorms(1), bnorms(cgto_shell_B%number_of_primitives), cnorms(1), dnorms(cgto_shell_D%number_of_primitives))
!        allocate(aexps(1),   bexps(cgto_shell_B%number_of_primitives),  cexps(1),  dexps(cgto_shell_D%number_of_primitives))
!        allocate(acoefs(1), bcoefs(cgto_shell_B%number_of_primitives), ccoefs(1), dcoefs(cgto_shell_D%number_of_primitives),eri_int(n_integrals,two_el_column))
!        allocate(int_index_eri(4,n_integrals))
!        ab_is_continuum = .false.
!        do_tails_for_this_quartet = .false.
!        tgt_pair = 0
!        !Shell A
!        lena = 1
!        la = dbg_cgto(bto_ind_A)%l
!        xa = dbg_cgto(bto_ind_A)%center(1)
!        ya = dbg_cgto(bto_ind_A)%center(2)
!        za = dbg_cgto(bto_ind_A)%center(3)
!        anorms(1) = dbg_cgto(bto_ind_A)%norms(1)*dbg_cgto(bto_ind_A)%norm
!        acoefs(1) = dbg_cgto(bto_ind_A)%contractions(1)
!        aexps(1) = dbg_cgto(bto_ind_A)%exponents(1)
!        ind_a = starting_index_A
!        !Shell C
!        lenc = dbg_cgto(bto_ind_C)%number_of_primitives
!        lc = dbg_cgto(bto_ind_C)%l
!        xc = dbg_cgto(bto_ind_C)%center(1)
!        yc = dbg_cgto(bto_ind_C)%center(2)
!        zc = dbg_cgto(bto_ind_C)%center(3)
!        cnorms(1:lenc) = dbg_cgto(bto_ind_C)%norms(1:lenc)*dbg_cgto(bto_ind_C)%norm
!        ccoefs(1:lenc) = dbg_cgto(bto_ind_C)%contractions(1:lenc)
!        cexps(1:lenc) = dbg_cgto(bto_ind_C)%exponents(1:lenc)
!        ind_c = starting_index_C
!        !Shell B
!        lenb = cgto_shell_B%number_of_primitives
!        lb = cgto_shell_B%l
!        xb = cgto_shell_B%center(1)
!        yb = cgto_shell_B%center(2)
!        zb = cgto_shell_B%center(3)
!        bnorms(1:lenb) = cgto_shell_B%norms(1:lenb)*cgto_shell_B%norm
!        bcoefs(1:lenb) = cgto_shell_B%contractions(1:lenb)
!        bexps(1:lenb) = cgto_shell_B%exponents(1:lenb)
!        ind_b = starting_index_B
!        !Shell D
!        lend = cgto_shell_D%number_of_primitives
!        ld = cgto_shell_D%l
!        xd = cgto_shell_D%center(1)
!        yd = cgto_shell_D%center(2)
!        zd = cgto_shell_D%center(3)
!        dnorms(1:lend) = cgto_shell_D%norms(1:lend)*cgto_shell_D%norm
!        dcoefs(1:lend) = cgto_shell_D%contractions(1:lend)
!        dexps(1:lend) = cgto_shell_D%exponents(1:lend)
!        ind_d = starting_index_D
!
!        call eri(lena,xa,ya,za,anorms,la,aexps,acoefs,ind_a, lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,ind_b, lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,ind_c, lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,ind_d,&
!                & two_el_column,int_index_eri,keep_ab_cd_order, do_tails_for_this_quartet,ab_is_continuum,tgt_prop,tgt_pair,eri_tail_int,bto_shell_A%bspline_grid%B,eri_int)
!
!        call compare_print_2el_ints(tag,integrals,int_index,eri_int,int_index_eri,n_integrals,two_el_column)

  end subroutine BGBG_shell_integrals

end module bto_gto_integrals_mod
