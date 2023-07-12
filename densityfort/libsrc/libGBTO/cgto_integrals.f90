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
module cgto_integrals
   use precisn
   use cgto_hgp
   use basis_data_generic_mod
   use gto_routines
   use utils

   private

   public GG_initialize, GG_shell_integrals, GGGG_shell_integrals, GGGG_initialize, GGGG_final

   real(kind=cfp), allocatable :: prop_tail(:), na_tail(:), nari(:), tgt_prop(:,:)
   integer, allocatable :: mom_index(:,:), ij_mapping(:,:)
   logical, allocatable :: is_continuum(:)
   logical :: tails_1el, tails_2el, keep_ab_cd_order_saved, two_p_continuum_saved
   real(kind=cfp) :: rmat_radius
   integer :: indexing_method

contains

   !> Determine wheter tails will be subtracted when the integrals will be calculated. If so then it calculates the norms of the continuum shells.
   !> \warning It is assumed that the CGTO shells on input ARE NOT NORMALIZED TO THE R-MATRIX SPHERE.
   subroutine GG_initialize(CGTO_shells,a)
      implicit none
      type(CGTO_shell_data_obj), intent(in) :: CGTO_shells(:)
      real(kind=cfp), intent(in) :: a

      integer :: err, i

         if (a > 0.0_cfp) then
            tails_1el = .true.
            rmat_radius = a

         else
            tails_1el = .false.
            rmat_radius = a
         endif

   end subroutine GG_initialize

   !> \warning This routine assumes that the arrays for output int_index and integrals have been allocated to the correct dimensions and that the column indices *_column
   !> are allowed, i.e. if they are greater than 0 they must be within the limits of the array integrals.
   !> \warning Not thread safe. We assume that the shells on input have been normalized to the requested R-matrix radius.
   subroutine GG_shell_integrals(shell_A,shell_B,A,B,starting_index_A,starting_index_B,&
                                     &use_spherical_cgto_alg,max_property_l,property_center,symmetry_data,&
                                     &olap_column,kei_column,prop_column,nai_column,one_elham_column,int_index,integrals)
      use cgto_hgp, only: sph_olap_kei, sph_nari, sph_mult_mom, index_1el
      use eri_sph_coord, only: olap_kei_sph
      use gto_routines, only: olap_kei_tail, nari_tail, prop_cms_tail
      use symmetry
      implicit none
      type(CGTO_shell_data_obj), target, intent(in) :: shell_A, shell_B
      real(kind=cfp), intent(in) :: property_center(3)
      type(symmetry_obj), intent(in) :: symmetry_data
      integer, intent(in) :: A,B, starting_index_A, starting_index_B, olap_column,kei_column,nai_column,prop_column,&
                             one_elham_column,max_property_l
      logical, intent(in) :: use_spherical_cgto_alg
      !We assume that these two arrays have been allocated to the appropriate dimensions:
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)

      logical :: do_tails_for_this_pair,A_is_continuum,B_is_continuum
      integer :: terms, sph_shell_a, sph_shell_b, small_l, i, j, ind_Ap, ind_Bp, no_prop, cnt, err
      real(kind=cfp) :: olap_tail,kei_tail_ab,bloch_ab,kei_tail_ba,bloch_ba
      type(CGTO_shell_data_obj), pointer :: shell_Ap, shell_Bp

         !If do_tails_for_this_pair == .true. then both shells are at CMS and non_zero_at_boundary so we have to subtract the tail integrals.
         do_tails_for_this_pair = .false.

         if (tails_1el) then
            A_is_continuum = shell_A%is_continuum()
            B_is_continuum = shell_B%is_continuum()
            if (A_is_continuum .and. B_is_continuum) do_tails_for_this_pair = .true.
         endif

         sph_shell_a = 2*shell_A%l+1
         sph_shell_b = 2*shell_B%l+1
         terms = sph_shell_a*sph_shell_b

         if (one_elham_column > 0 .and. (kei_column == 0 .or. nai_column == 0)) then
            call xermsg ('cgto_integrals', 'GG_shell_integrals', &
                         'If one electron Hamiltonian integrals are requested then the kinetic energy and &
                         &nuclear attraction integrals must be requested too.', 1, 1)
         end if

!------- Calculate integrals over all space:
         if (olap_column > 0 .or. kei_column > 0) then !overlap and/or kinetic energy integrals

            !todo the selection rule should be incorporated into sph_olap_kei:
            if (shell_A % center(1) == shell_B % center(1) .and. &
                shell_A % center(2) == shell_B % center(2) .and. &
                shell_A % center(3) == shell_B % center(3) .and. shell_A % l /= shell_B % l) then
               integrals(1:terms,olap_column) = 0.0_cfp
               integrals(1:terms,kei_column) = 0.0_cfp
               call index_1el(shell_A%l,shell_B%l,starting_index_A,starting_index_B,1,int_index)
            else

            if (use_spherical_cgto_alg) then
               call olap_kei_sph (shell_A % number_of_primitives, &
                                  shell_A % center(1), &
                                  shell_A % center(2), &
                                  shell_A % center(3), &
                                  shell_A % norm, &
                                  shell_A % norms, &
                                  shell_A % l, &
                                  shell_A % exponents, &
                                  shell_A % contractions, &
                                  starting_index_A, &
                                  shell_B % number_of_primitives, &
                                  shell_B % center(1), &
                                  shell_B % center(2), &
                                  shell_B % center(3), &
                                  shell_B % norm, &
                                  shell_B % norms, &
                                  shell_B % l, &
                                  shell_B % exponents, &
                                  shell_B % contractions, &
                                  starting_index_B,&
                                  olap_column, kei_column, integrals, int_index)
            else !cartesian-based alg.
               call sph_olap_kei (shell_A % number_of_primitives, &
                                  shell_A % center(1), &
                                  shell_A % center(2), &
                                  shell_A % center(3), &
                                  shell_A % norm, &
                                  shell_A % norms, &
                                  shell_A % l, &
                                  shell_A % exponents, &
                                  shell_A % contractions, &
                                  starting_index_A, &
                                  shell_B % number_of_primitives, &
                                  shell_B % center(1), &
                                  shell_B % center(2), &
                                  shell_B % center(3), &
                                  shell_B % norm, &
                                  shell_B % norms, &
                                  shell_B % l, &
                                  shell_B % exponents, &
                                  shell_B % contractions, &
                                  starting_index_B, &
                                  olap_column, kei_column, integrals, int_index)
            endif
            endif
         endif

         if (nai_column > 0) then !nuclear attraction integrals

            err = check_real_array_size(nari,terms)
            if (err .ne. 0) stop "nuclear_attraction_integrals: memory allocation error 2"
   
            !loop over the nuclei.
            integrals(1:terms,nai_column) = 0.0_cfp
            do i=1,symmetry_data%no_nuc
   
            !Calculate the overlap and kinetic energy integrals over the current pair of shells of contracted GTOs and for the nucleus i. The result in nari might correspond to (ab) being swapped 
            !depending on their L. The indices of the functions corresponding to the integrals are in int_index.
            call sph_nari (shell_A % number_of_primitives, &
                           shell_A % center(1), &
                           shell_A % center(2), &
                           shell_A % center(3), &
                           shell_A % norm, &
                           shell_A % norms, &
                           shell_A % l, &
                           shell_A % exponents, &
                           shell_A % contractions, &
                           starting_index_A, &
                           shell_B % number_of_primitives, &
                           shell_B % center(1), &
                           shell_B % center(2), &
                           shell_B % center(3), &
                           shell_B % norm, &
                           shell_B % norms, &
                           shell_B % l, &
                           shell_B % exponents, &
                           shell_B % contractions, &
                           starting_index_B, &
                           symmetry_data % nucleus(i) % center(1), &
                           symmetry_data % nucleus(i) % center(2), &
                           symmetry_data % nucleus(i) % center(3), nari, int_index)
               !accumulate the results into the nari_nuc array holding the final results and multiply the contributions by the nuclear charge.
               integrals(1:terms,nai_column) = integrals(1:terms,nai_column) - symmetry_data%nucleus(i)%charge*nari(1:terms)
            enddo

            if (do_tails_for_this_pair) then
               err = check_real_array_size(na_tail,terms)
               if (err .ne. 0) stop "nuclear_attraction_integrals: memory allocation error 1"
            endif
         endif

         if (prop_column > 0) then
            !Calculate the property integrals over the current pair of shells of contracted GTOs. The results in prop
            !The indices of the (ab) functions corresponding to the integrals in integrals(i,:) are in int_index(1,i), int_index(2,i). 
            !The indices of the functions corresponding to each integral are ordered so that int_index(1,i) .ge. int_index(2,i). This allows us to calculate the ordered index for each integral.
            call sph_mult_mom (shell_A % number_of_primitives, &
                               shell_A % center(1), &
                               shell_A % center(2), &
                               shell_A % center(3), &
                               shell_A % norm, &
                               shell_A % norms, &
                               shell_A % l, &
                               shell_A % exponents, &
                               shell_A % contractions, &
                               starting_index_A, & !the index of the first function in the a shell; it is assumed that the other functions's in the shell have sequential indices, same for the b shell.
                               max_property_l, property_center(1), property_center(2), property_center(3), & !the maximum L of the property and its center
                               shell_B % number_of_primitives, &
                               shell_B % center(1), &
                               shell_B % center(2), &
                               shell_B % center(3), &
                               shell_B % norm, &
                               shell_B % norms, &
                               shell_B % l, &
                               shell_B % exponents, &
                               shell_B % contractions, &
                               starting_index_B, prop_column, integrals, int_index)
   
            !Note that the integrals in prop are ordered in the opposite way than sph_shell_a, sph_shell_b if l_b > l_a.
            if (do_tails_for_this_pair) then
               no_prop = (max_property_l+1)**2
               err = check_real_array_size(prop_tail,terms*no_prop)
               if (err .ne. 0) stop "property_integrals: memory allocation error"
            endif
         endif

         !Subtract the tail integrals:
         if (do_tails_for_this_pair) then

            !Did we swap (ab) for (ba) in sph_olap_kei and sph_nari?
            !Determining la,lb below allows us to calculate indices of each of the integrals corresponding to a particular combination of functions with M_a,M_b.
            if (shell_A%l .ge. shell_B%l) then !no
               shell_Ap => shell_A
               shell_Bp => shell_B
               ind_Ap = starting_index_A
               ind_Bp = starting_index_B
            else !yes
               shell_Ap => shell_B
               shell_Bp => shell_A
               ind_Ap = starting_index_B
               ind_Bp = starting_index_A

               sph_shell_a = 2*shell_Ap%l+1
               sph_shell_b = 2*shell_Bp%l+1
            endif

            !Olap/kei tails:
            if (olap_column > 0 .or. kei_column > 0) then
               if (shell_A%l .eq. shell_B%l) then !selection rule for calculation of the olap/kei tail integrals

                  !Calculate the tail integrals and the Bloch terms for the (ab) and (ba) pairs. Since the KE integrals over the R-matrix sphere are not symmetric we have to subtract below the 
                  !KE tails corresponding to exactly the combination of functions given by the indices int_index(1,i) and int_index(2,i). These indices define which function is a and b in the 
                  !i-th -1/2*<a|Nabla|b> integral.
         
                  !(ap bp) combination
                  call olap_kei_tail (shell_Ap % l, &
                                      shell_Ap % number_of_primitives,  shell_Bp % number_of_primitives, &
                                      shell_Ap % exponents,             shell_Bp % exponents, &
                                      shell_Ap % contractions,          shell_Bp % contractions, &
                                      shell_Ap % norm, &
                                      shell_Ap % norms, &
                                      shell_Bp % norm, &
                                      shell_Bp % norms, &
                                      rmat_radius, olap_tail, kei_tail_ab,  bloch_ab)
         
                  !(bp ap) combination
                  call olap_kei_tail (shell_Bp % l, &
                                      shell_Bp % number_of_primitives,  shell_Ap % number_of_primitives, &
                                      shell_Bp % exponents,             shell_Ap % exponents, &
                                      shell_Bp % contractions,          shell_Ap % contractions, &
                                      shell_Bp % norm, &
                                      shell_Bp % norms, &
                                      shell_Ap % norm, &
                                      shell_Ap % norms, &
                                      rmat_radius, olap_tail, kei_tail_ba, bloch_ba)
         
                  small_l = min(shell_Ap%l,shell_Bp%l) !loop over the pairs of functions which have M_a=M_b; only for those pairs are the tails nonzero
                  do j=-small_l,small_l
         
                     i = (j+shell_Ap%l+1) + sph_shell_a*(j+shell_Bp%l+1 -1) !index of the integral corresponding to the pair of functions with M_a=M_b
                     if (olap_column > 0) integrals(i,olap_column) = (integrals(i,olap_column) - olap_tail)
         
                     if (kei_column > 0) then
                        !Subtract the KE tail integral depending on whether the indices define the integral as -1/2*<a|Nabla|b> or -1/2*<b|Nabla|a>.
                        !Note that we assume here sequential indexing in both shells.
                        if (int_index(1,i) .ge. ind_Ap .and. int_index(1,i) .le. ind_Ap + 2*shell_Ap%l) then !the int_index(1,i) belongs to the a shell: integrals(i,kei_column) = -1/2*<a|Nabla|b>
                           integrals(i,kei_column) = (integrals(i,kei_column) - kei_tail_ab + bloch_ab)
                        else !the int_index(1,i) belongs to the b shell: integrals(i,kei_column) = -1/2*<b|Nabla|a>
                           integrals(i,kei_column) = (integrals(i,kei_column) - kei_tail_ba + bloch_ba)
                        endif
                     endif
      
                  enddo !j 
               endif
            endif

            !Nari tails:
            if (nai_column > 0) then
               !todo the loop over nuclei should go into nari_tail since the couplings will always be the same; this is a silly way of calculating the NARI tail
               na_tail(1:terms) = 0.0_cfp
               do i=1,symmetry_data%no_nuc
                  call nari_tail (symmetry_data % nucleus(i) % center(1), &
                                  symmetry_data % nucleus(i) % center(2), &
                                  symmetry_data % nucleus(i) % center(3), &
                                  shell_Ap % l,                     shell_Bp % l, &
                                  shell_Ap % number_of_primitives,  shell_Bp % number_of_primitives, &
                                  shell_Ap % exponents,             shell_Bp % exponents, &
                                  shell_Ap % contractions,          shell_Bp % contractions, &
                                  shell_Ap % norm, &
                                  shell_Ap % norms, &
                                  shell_Bp % norm, &
                                  shell_Bp % norms, &
                                  rmat_radius, nari)
                  !Accumulate the tail integrals
                  na_tail(1:terms) = na_tail(1:terms) - symmetry_data%nucleus(i)%charge*nari(1:terms)
               enddo
               integrals(1:terms,nai_column) = (integrals(1:terms,nai_column) - na_tail(1:terms))
            endif

            !Property tails:
            if (prop_column > 0) then
               !the tails are in prop_tail and ordered in the same way as in prop.
               call prop_cms_tail (max_property_l, &
                                   shell_A % l,                     shell_B % l, &
                                   shell_A % number_of_primitives,  shell_B % number_of_primitives, &
                                   shell_A % exponents,             shell_B % exponents, &
                                   shell_A % contractions,          shell_B % contractions, &
                                   shell_A % norm, &
                                   shell_A % norms, &
                                   shell_B % norm, &
                                   shell_B % norms, &
                                   rmat_radius, prop_tail)
               cnt = 0
               do j=prop_column,prop_column+no_prop-1
                  integrals(1:terms,j) = (integrals(1:terms,j) - prop_tail(cnt+1:cnt+terms))
                  cnt = cnt + terms
               enddo !j
            endif

         endif !do_tails_for_this_pair

         !Note that this requires that the KEI and NARI are always calculated if one_elham_column > 0, i.e. nai_column > 0 .and. kei_column > 0.
         if (one_elham_column > 0) then
            integrals(1:terms,one_elham_column) = integrals(1:terms,nai_column) + integrals(1:terms,kei_column)
         end if

   end subroutine GG_shell_integrals

   !> \warning We assume that the overal norms of the contractions have been multiplied into CGTO_shells()%norms!
   subroutine GGGG_initialize(CGTO_shells,shell_starting_indices,tol,a,keep_ab_cd_order,two_p_continuum,indexing_method_inp)
      use const
      use cgto_hgp, only: sph_mult_mom
      implicit none
      type(CGTO_shell_data_obj), target, intent(in) :: CGTO_shells(:)
      integer, intent(in) :: shell_starting_indices(:), indexing_method_inp
      real(kind=cfp), intent(in) :: a, tol
      logical, intent(in) :: two_p_continuum,keep_ab_cd_order

      integer :: i, j, k, n, ij, number_of_target_shells, number_of_shells, err, max_l_tgt, max_l_cont, tgt_pairs, mom_space, &
                 l_max, no_prop, sph_shell_A, sph_shell_B, sph_shell_AB, n_prim
      real(kind=cfp), allocatable :: shell_prop(:,:)
      logical :: cgto_continuum

         !Transfer some of the input parameters into the module variables.
         keep_ab_cd_order_saved = keep_ab_cd_order
         two_p_continuum_saved = two_p_continuum
         rmat_radius = a
         indexing_method = indexing_method_inp

         if (indexing_method > 2 .or. indexing_method <= 0) then
            call xermsg ('cgto_integrals', 'GGGG_initialize', 'On input indexing_method was out of range [1,2].', 1, 1)
         end if

         number_of_shells = size(CGTO_shells)

         if (allocated(is_continuum)) deallocate(is_continuum)
         allocate(is_continuum(number_of_shells),stat=err)
         if (err .ne. 0) call xermsg('cgto_integrals', 'GGGG_initialize','Memory allocation error.',err,1)
         
         !Analyze the CGTO basis set:
         max_l_cont = -1
         max_l_tgt = -1
         number_of_target_shells = 0
         is_continuum = .false.
         n_prim = 0
         do i=1,number_of_shells
            n_prim = max(n_prim,CGTO_shells(i)%number_of_primitives)
            if (CGTO_shells(i)%is_continuum()) then
               is_continuum(i) = .true.
               max_l_cont = max(max_l_cont,CGTO_shells(i)%l)
            else
               max_l_tgt = max(max_l_tgt,CGTO_shells(i)%l)
               number_of_target_shells = number_of_target_shells + 1
            endif
         enddo !i
         tgt_pairs = number_of_target_shells*(number_of_target_shells+1)/2 !number of pairs of shells of target GTOs

         if (number_of_target_shells .eq. number_of_shells) then
            cgto_continuum = .false.
         else
            cgto_continuum = .true.
         endif

         !tell the user if the tails will be calculated based on the value of the R-matrix radius
         if (a > 0.0_cfp) then
            if (.not.(cgto_continuum)) then
               tails_2el = .false.
               write(stdout,'("a > 0 but the basis contains no CGTO continuum functions so no CGTO tails have to be subtracted.")')
            elseif (cgto_continuum .and. number_of_target_shells > 0) then !there are some functions in the basis that represent the bound electrons
               if (.not.(two_p_continuum)) then
                  write(stdout,'("Evaluation of integrals of the type [continuum,continuum|continuum,continuum] &
                                 &and [continuum,continuum|continuum,target] will be skipped.")')
                  write(stdout,'("Tail integrals for ONE ELECTRON IN THE CONTINUUM and a = ",e25.15," a.u. &
                                 &will be calculated and subtracted.")') a
                  tails_2el = .true.
               else
                  write(stdout,'("Tail integrals for TWO ELECTRONS IN THE CONTINUUM and a = ",e25.15," a.u. &
                                 &will be calculated and subtracted.")') a
                  tails_2el = .true.
               endif
               write(stdout,'("Number of shells of target functions: ",i0)') number_of_target_shells
               write(stdout,'("ASSUMPTION: all non-CMS functions fit completely inside the R-matrix sphere.")')
            elseif (cgto_continuum .and. number_of_target_shells .eq. 0) then !the basis contains only continuum GTOs
               if (.not.(two_p_continuum)) then
                  !Note that this is true only for 1p in the continuum.
                  write(stdout,'("The GTO basis set does not contain any functions representing the bound electrons. &
                                 &No 1p tails will be subtracted.")')
                  tails_2el = .false.
               else !For 2p in the continuum have to subtract the Continuum-Continuum 2p tails.
                  write(stdout,'("Tail integrals for TWO ELECTRONS IN THE CONTINUUM and a = ",e25.15," a.u. &
                                 &will be calculated and subtracted.")') a
                  tails_2el = .true.
               endif
            endif
         else !R-matrix radius .le. 0 or no continuum functions in the basis signifies that integrals over all space are required
            write(stdout,'("R-matrix radius is not > 0.0_cfp; no tail integrals will be calculated.")')
            tails_2el = .false.
         endif

         if (two_p_continuum .and. tails_2el) then
            call xermsg ('cgto_integrals', 'GGGG_initialize', 'Tails for 2p in the continuum requested but not implemented.', 4, 1)
         end if

         if (tails_2el) then

            !space for the input data for the eri_tail routines.
            !We need the property integrals for the pair of the target GTO shells from l=0 up to l=l_max.
            l_max = 2*max_l_cont !The maximum target multipole moment needed is given by the sum of the continuum GTO's L values under the assumption that these GTOs are centered on CMS.

            !Space for the property integrals for all unique pairs of shells of target GTOs.
            mom_space = (2*max_l_tgt+1)**2*(2*max_l_cont+1)**2
            if (allocated(mom_index)) deallocate(mom_index)
            if (allocated(tgt_prop)) deallocate(tgt_prop)
            if (allocated(ij_mapping)) deallocate(ij_mapping)
            allocate(shell_prop(mom_space,(l_max+1)**2),mom_index(1:2,mom_space), &
                     tgt_prop(mom_space,tgt_pairs),ij_mapping(number_of_target_shells,number_of_target_shells),stat=err)
            if (err .ne. 0) call xermsg ('two_particle_integrals_obj', 'GGGG_initialize', 'Memory allocation 3 failed', err, 1)

            mom_space  = (l_max+1)**2 !from now on mom_space is the total number of (l,m) property values calculated for each pair of target GTOs.

            !Note that the loops assume that the sequence numbers of target shells preceed the continuum ones!
            ij_mapping = 0
            ij = 0
            do i=1,number_of_target_shells
               do j=1,i

                  if (is_continuum(i) .or. is_continuum(j)) then
                    call xermsg ('two_particle_integrals_obj', 'GGGG_initialize', &
                                 'Error in GTO shell ordering: all target shells must preceed the continuum.', 5, 1)
                  end if

                  ij = ij + 1

                  ij_mapping(i,j) = ij
                  ij_mapping(j,i) = ij

                  !In order to calculate the tail integrals for 1p in the continuum we require some property integrals for the target GTO shells.
                  !todo We can precalculate also the continuum 1p potential integral terms.

                  sph_shell_a = 2*CGTO_shells(i)%l+1
                  sph_shell_b = 2*CGTO_shells(j)%l+1
                  sph_shell_ab = sph_shell_a*sph_shell_b

                  no_prop = mom_space*sph_shell_ab !Number of target properties to evaluate for this pair of shells
      
                  !The indices of the (ab) functions corresponding to the integrals in tgt_prop(i,ab) are in mom_index(1,i), mom_index(2,i). The shells are ordered in sph_mult_mom so that l_a .ge. l_b, where
                  !l_a = max(la,lb), l_b = min(la,lb) and l_a is the L value in the a-shell and l_b is the L value in the b-shell.
                  !The indices of the functions corresponding to each integral are ordered so that mom_index(1,i) .ge. mom_index(2,i) but we don't need the indices for the tail integrals calculation.
                  !Note that placing the results in the column 'ab' of tgt_prop assumes that the target shells always preceed the continuum shells.
                  !Additionally, when the angular momenta of the shells are the same the only thing that decides on the order of the shells is the order of the loops over i,j. Therefore this order must be
                  !preserved when the target moments are needed in GGGG_shell_integrals!!!
                  shell_prop = 0.0_cfp
                  call sph_mult_mom (CGTO_shells(i) % number_of_primitives, &
                                     CGTO_shells(i) % center(1), &
                                     CGTO_shells(i) % center(2), &
                                     CGTO_shells(i) % center(3), &
                                     1.0_cfp, &
                                     CGTO_shells(i) % norms, &
                                     CGTO_shells(i) % l, &
                                     CGTO_shells(i) % exponents, &
                                     CGTO_shells(i) % contractions, &
                                     shell_starting_indices(i), &
                                     l_max, 0.0_cfp, 0.0_cfp, 0.0_cfp, & !the maximum L of the property and its center
                                     CGTO_shells(j) % number_of_primitives, &
                                     CGTO_shells(j) % center(1), &
                                     CGTO_shells(j) % center(2), &
                                     CGTO_shells(j) % center(3), &
                                     1.0_cfp, &
                                     CGTO_shells(j) % norms, &
                                     CGTO_shells(j) % l, &
                                     CGTO_shells(j) % exponents, &
                                     CGTO_shells(j) % contractions, &
                                     shell_starting_indices(j), &
                                     1, shell_prop, mom_index) !save the property integrals in the column 'ab' of the tgt_prop array

                  !Transfer the properties into a single column: this is how the properties are needed in eri_tail.
                  do n=1,mom_space
                     tgt_prop((n-1)*sph_shell_ab+1:n*sph_shell_ab,ij) = shell_prop(1:sph_shell_ab,n)
                  enddo
     
                  do k=1,no_prop
                     if (abs(tgt_prop(k,ij)) < tol) tgt_prop(k,ij) = 0.0_cfp
                  enddo

               enddo !j
            enddo !i
         endif

   end subroutine GGGG_initialize

   !> \warning We assume that the overal norms of the contractions have been multiplied into CGTO_shells()%norms and that GGGG_initialize has been called before!
   subroutine GGGG_shell_integrals(shell_A,shell_B,shell_C,shell_D,A,B,C,D,&
                                   starting_index_A,starting_index_B,starting_index_C,starting_index_D,&
                                   use_spherical_cgto_alg,two_el_column,int_index,integrals)
      use cgto_hgp, only: eri
      use eri_sph_coord, only: eri_sph
      use gto_routines, only: eri_tail
      implicit none
      type(CGTO_shell_data_obj), intent(in) :: shell_A, shell_B, shell_C, shell_D
      integer, intent(in) :: starting_index_A, starting_index_B, starting_index_C, starting_index_D, A,B,C,D, two_el_column
      logical, intent(in) :: use_spherical_cgto_alg
      !We assume that these two arrays have been allocated to the appropriate dimensions:
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)

      integer :: tgt_pair, no_shell, sph_shell_a,sph_shell_b,sph_shell_c,sph_shell_d
      logical :: do_tails_for_this_quartet, ab_is_continuum, cnt(4)

         !Should the tails be subtracted for this quartet of shells?
         do_tails_for_this_quartet = .false.
         if (tails_2el) then
            !Yes: associate pointers to the target and continuum shells.
            if (is_continuum(A) .and. is_continuum(B) .and. (.not.(is_continuum(C)) .and. .not.(is_continuum(D)))) then !CC|TT
               do_tails_for_this_quartet = .true.
               ab_is_continuum = .true.
               tgt_pair = ij_mapping(C,D) !index of the target shell-pair
            elseif (.not.(is_continuum(A)) .and. .not.(is_continuum(B)) .and. (is_continuum(C) .and. is_continuum(D))) then !TT|CC
               do_tails_for_this_quartet = .true.
               ab_is_continuum = .false.
               tgt_pair = ij_mapping(A,B) !index of the target shell-pair
            else
               cnt(1:4) = (/is_continuum(A),is_continuum(B),is_continuum(C),is_continuum(D)/)
               if (count(cnt) > 2) then
                  call xermsg ('two_particle_integrals_obj', 'GGGG_shell_integrals', &
                               'Tails for 2p in the continuum have been requested but they are not implemented. &
                               &See comments in GGGG_shell_integrals.', 1, 1)
                  !Tails of the type (CC|CT) can be implemented in the current scheme only the calculation of the properties must be extended to include the CT pairs of shells.
                  !Tails of the type (CC|CC) cannot be evaluated analytically since both CC pairs overlap in the outer region. It would be easier to calculate the actual (CC|CC) integrals using the Legendre
                  !expansion (which is finite for continuum shells) integrating numerically over the radial coordinates (r1,r2) in the inner region.
               endif
            endif
         endif

         sph_shell_a = 2*shell_A%l+1
         sph_shell_b = 2*shell_B%l+1
         sph_shell_c = 2*shell_C%l+1
         sph_shell_d = 2*shell_D%l+1
         no_shell = sph_shell_a*sph_shell_b*sph_shell_c*sph_shell_d

         !calculate integrals over the (ab|cd) shells over all space; the results are in integrals(1:no_shell,two_el_column) and the ordered quartets of indices are in int_index.
         !If the variable do_tails_for_this_quartet is true, the tail integrals are calculated and subtracted. The variable ab_is_continuum determines which pair of shells is taken 
         !to be the continuum one: (ab) or (cd). Only if both (ab) or (cd) shells are continum and both (cd) or (ab) shells are target then there is a contribution to the 1p tail integral. 
         !We calculate the 1p tails using the Legendre expansion for the Coulomb interaction. The 2p tail integrals, i.e. (ab)/(cd) is continuum and (cd)/(ab) is continuum, are more difficult, 
         !involve a 2D integral that needs to be calculated numerically and we don't implement those. The target property integrals needed for the tail integral calculation are supplied 
         !in the tgt_prop array and the index of the target pair of shells is given by tgt_pair. Note that if the angular momenta of the target shells are the same the order of the pair of target 
         !shells for which the properties were calculated must be the same as the order in which they are used here: this must be ensured when calling this routine.
         if (use_spherical_cgto_alg) then
            call eri_sph (shell_A % number_of_primitives, &
                          shell_A % center(1), shell_A % center(2), shell_A % center(3), &
                          shell_A % norms, &
                          shell_A % l, &
                          shell_A % exponents, &
                          shell_A % contractions, &
                          starting_index_A, &
                          shell_B % number_of_primitives, &
                          shell_B % center(1), shell_B % center(2), shell_B % center(3), &
                          shell_B % norms, &
                          shell_B % l, &
                          shell_B % exponents, &
                          shell_B % contractions, &
                          starting_index_B, &
                          shell_C % number_of_primitives, &
                          shell_C % center(1), shell_C % center(2), shell_C % center(3), &
                          shell_C % norms, &
                          shell_C % l, &
                          shell_C % exponents, &
                          shell_C % contractions, &
                          starting_index_C, &
                          shell_D % number_of_primitives, &
                          shell_D % center(1), shell_D % center(2), shell_D % center(3), &
                          shell_D % norms, &
                          shell_D % l, &
                          shell_D % exponents, &
                          shell_D % contractions, &
                          starting_index_D, &
                          two_el_column, int_index, keep_ab_cd_order_saved, indexing_method, &
                          do_tails_for_this_quartet, ab_is_continuum, tgt_prop, tgt_pair, rmat_radius, & !extra data for the tail integral calculation and subtraction
                          integrals)
         else
            call eri (shell_A % number_of_primitives, &
                      shell_A % center(1), shell_A % center(2), shell_A % center(3), &
                      shell_A % norms, &
                      shell_A % l, &
                      shell_A % exponents, &
                      shell_A % contractions, &
                      starting_index_A, &
                      shell_B % number_of_primitives, &
                      shell_B % center(1), shell_B % center(2), shell_B % center(3), &
                      shell_B % norms, &
                      shell_B % l, &
                      shell_B % exponents, &
                      shell_B % contractions, &
                      starting_index_B, &
                      shell_C % number_of_primitives, &
                      shell_C % center(1), shell_C % center(2), shell_C % center(3), &
                      shell_C % norms, &
                      shell_C % l, &
                      shell_C % exponents, &
                      shell_C % contractions, &
                      starting_index_C, &
                      shell_D % number_of_primitives, &
                      shell_D % center(1), shell_D % center(2), shell_D % center(3), &
                      shell_D % norms, &
                      shell_D % l, &
                      shell_D % exponents, &
                      shell_D % contractions, &
                      starting_index_D, &
                      two_el_column, int_index, keep_ab_cd_order_saved, indexing_method, &
                      do_tails_for_this_quartet, ab_is_continuum, tgt_prop, tgt_pair, rmat_radius, & !extra data for the tail integral calculation and subtraction
                      integrals)
         endif

   end subroutine GGGG_shell_integrals

   !> Gets rid of all arrays local to this module related to evaluation of 2-electron integrals. Additionally, the arrays used by cgto_hgp and eri_sph are removed, too: these can be large for high ang. momentum.
   subroutine GGGG_final
      use cgto_hgp, only: cgto_hgp_final
      use eri_sph_coord, only: eri_sph_coord_final
      implicit none

         !todo deallocation of the local arrays tgt_prop, etc. gives strange erros even if debugging is used!!! Compiler error?

         call cgto_hgp_final
         call eri_sph_coord_final
   
   end subroutine GGGG_final

end module cgto_integrals
