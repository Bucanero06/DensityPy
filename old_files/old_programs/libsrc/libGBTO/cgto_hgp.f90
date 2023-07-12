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
!Replacing the array arguments with assumed size (+ Wigner-Eckart) increases the performance by ~10% for simple scattering run tests since the compiler knows that the arrays are unit-stride.
!todo note the potential loss of precision in the cartesian to spherical transformation due to very large transformation coefficients for large l values.
!Head-Gordon-Pople method for contracted 2-electron integrals (i.e. modified Obara-Saika method) and the Obara-Saika method for 1-electron integrals.
!todo all triplets of cartesian coordinates should be input as 3-member arrays to make sure the numbers are aligned in memory.
!todo define l1_cache_size and use it to manage sizes of the buffers: check if the source-target buffers fit inside the cache. If not then reallocate them (if possible) to their minimal sizes.
module cgto_hgp
 use precisn
 use phys_const, only: pi, pi54, rtwo
 use gto_routines, only: boys_function_obj, abcd_to_cdab, check_real_array_size, index_1el, &
                         index_2el, eri_tail_shell, reorder_and_index_2el

   implicit none

   !max_sum_l variable is used to hold the value of the dimension of the cart_* arrays; it must be initially set to -1 to force allocation of these arrays in the calc_can routine.
   !max_l holds the value of the largest L encountered so far. Initially it must be set to 0. It is used to control the size of the array cart_to_sph containing the coefficients for mapping of cartesian
   !shells into spherical shells.
   integer, private :: max_sum_l = -1, max_l = 0
   integer, allocatable, private :: cart_l(:), cart_m(:), cart_n(:)

   !Array with integrals (ab|cd) over the four shells of contracted spherical GTOs. Can be seen but not modified from outside this module. This linear array effectively emulates a 4D array with dimensions:
   !(sph_shell_a,sph_shell_b,sph_shell_c,sph_shell_d), where
   !   sph_shell_a = 2*la+1,
   !   sph_shell_b = 2*lb+1,
   !   sph_shell_c = 2*lc+1,
   !   sph_shell_d = 2*ld+1.
   ! Hence sph_ints(i), where i is:
   ! i = ma+la+1 + (mb+lb+1-1)*sph_shell_a + (mc+lc+1-1)*sph_shell_a*sph_shell_b + (md+ld+1-1)*sph_shell_a*sph_shell_b*sph_shell_c
   ! corresponds to the integral over the basis functions with the M values (ma,mb,mc,md).
!   real(kind=cfp), allocatable, protected :: sph_ints(:)

   real(kind=cfp), allocatable, private :: cart_to_sph(:), cart_mult_mom(:), contr_cart_mult_mom(:), cart_olap_buf(:), &
                                           cart_kei_buf(:), contr_cart_olap(:), contr_cart_kei(:)

   !indices for each L of the beginning and end of the coefficients for the cartesian to spherical mapping in cart_to_sph
   integer, allocatable, private :: sph_map_start(:), nz(:), nz_can(:), nz_can_ab(:,:)

   !buffers needed in the cartesian -> spherical harmonic transformation
   real(kind=cfp), allocatable, private :: half_sph(:), transf(:), c_nz(:), c_ab(:)

   !buffers that store at the end of the corresponding step of the algorithm the class of integrals required for the next step
   real(kind=cfp), allocatable, private :: contr_et(:), hrr1_tgt(:), hrr2_tgt(:)

   !Needed for calculation of the multipole moment integrals
   real(kind=cfp), allocatable, private :: shell_mom(:)

   !object used for evaluation of the Boys function
   type(boys_function_obj), private :: boys

   !auxiliary to keep the 2p tail integrals
   real(kind=cfp), allocatable :: eri_tail_int(:)

   !$OMP THREADPRIVATE (max_sum_l, max_l, cart_l, cart_m, cart_n, cart_to_sph, cart_mult_mom, contr_cart_mult_mom, &
   !$OMP &              cart_olap_buf, cart_kei_buf, contr_cart_olap, contr_cart_kei, sph_map_start, nz, nz_can, nz_can_ab, &
   !$OMP &              half_sph, transf, c_nz, c_ab, contr_et, hrr1_tgt, hrr2_tgt, shell_mom, boys, eri_tail_int)

contains

   subroutine cgto_hgp_final
      implicit none

         max_sum_l = -1; max_l = 0
         if (allocated(cart_l)) deallocate(cart_l)
         if (allocated(cart_m)) deallocate(cart_m)
         if (allocated(cart_n)) deallocate(cart_n) 
         if (allocated(cart_to_sph)) deallocate(cart_to_sph)
         if (allocated(cart_mult_mom)) deallocate(cart_mult_mom)
         if (allocated(contr_cart_mult_mom)) deallocate(contr_cart_mult_mom)
         if (allocated(cart_olap_buf)) deallocate(cart_olap_buf)
         if (allocated(cart_kei_buf)) deallocate(cart_kei_buf)
         if (allocated(contr_cart_olap)) deallocate(contr_cart_olap)
         if (allocated(contr_cart_kei)) deallocate(contr_cart_kei) 
         if (allocated(sph_map_start)) deallocate(sph_map_start)
         if (allocated(nz)) deallocate(nz)
         if (allocated(nz_can)) deallocate(nz_can)
         if (allocated(nz_can_ab)) deallocate(nz_can_ab) 
         if (allocated(half_sph)) deallocate(half_sph)
         if (allocated(transf)) deallocate(transf)
         if (allocated(c_nz)) deallocate(c_nz)
         if (allocated(c_ab)) deallocate(c_ab) 
         if (allocated(contr_et)) deallocate(contr_et)
         if (allocated(hrr1_tgt)) deallocate(hrr1_tgt)
         if (allocated(hrr2_tgt)) deallocate(hrr2_tgt) 
         if (allocated(shell_mom)) deallocate(shell_mom) 
         if (allocated(eri_tail_int)) deallocate(eri_tail_int) 

   end subroutine cgto_hgp_final

   !calculates the l,m,n triplets for all cartesian shells up to l=sum_l and stores the result in cart_l, cart_m, cart_n
   subroutine calc_can(sum_l)
      implicit none
      integer, intent(in) :: sum_l
      integer(kind=shortint) :: can, l, m, n, i, shell, space, err

         if (sum_l > max_sum_l) then

            !update the control variable of the largest angular momentum for which the canonical indices have been calculated
            max_sum_l = sum_l

            !space for the canonical indices
            space = ncart(max_sum_l)

            if (allocated(cart_l)) deallocate(cart_l)
            if (allocated(cart_m)) deallocate(cart_m)
            if (allocated(cart_n)) deallocate(cart_n)

            allocate(cart_l(space),cart_m(space),cart_n(space),stat=err)
            if (err .ne. 0) stop "cart_l, cart_m, cart_n allocation failed."
   
            can = 0
            cart_l = 0; cart_m = 0; cart_n = 0
            do shell=0,max_sum_l
               do i=0,shell
                  l = shell - i
                  do n=0,i
                     m = i-n
                     can = can + 1    !the full canonical index
                     cart_l(can) = l  !X power
                     cart_m(can) = m  !Y power
                     cart_n(can) = n  !Z power
                     !print *,can,l,m,n,can(shell,l,n)
                  enddo
               enddo
            enddo
        endif

   end subroutine calc_can

   subroutine allocate_space(la,lb,lc,ld,space_vrr_tgt,space_et_tgt,space_sph_ints,space_vrr_buf,space_et_buf,&
                             space_hrr1_buf,space_hrr2_buf,space_hrr1_tgt,space_hrr2_tgt)
      implicit none
      integer, intent(in) :: la, lb, lc, ld
      integer, intent(out) :: space_vrr_tgt, space_et_tgt, space_sph_ints, space_vrr_buf, space_et_buf, space_hrr1_buf, &
                              space_hrr2_buf, space_hrr1_tgt, space_hrr2_tgt

      integer :: i, j, k, l, err, sum_l, space, s_y, s_xp1, s_zp1
 
         sum_l = la+lb+lc+ld

         !Space for the final spherical integrals: note that it is actually larger than minimal which is (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1)
         space_sph_ints = (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1)
!         err = check_real_array_size(sph_ints,space_sph_ints); if (err .ne. 0) stop "sph_ints allocation failed."
 
         !space for the intermediate half spherical half cartesian integrals: (ab|cd]
         space_hrr2_tgt = (2*la+1)*(2*lb+1)*nshell(lc)*nshell(ld)
         err = check_real_array_size(hrr2_tgt,space_hrr2_tgt); if (err .ne. 0) stop "hrr2_tgt allocation failed."

         !space needed to store the target class for the VRR step
         space_vrr_tgt = ncart(sum_l)

         !compute the minimum amount of memory we will need for the VRR step
         space_vrr_buf = 0
         do j=0,sum_l

            !i = space required to store the row of auxiliaries with m=0,...,sum_l+1-s and the angular momentum 'j'.
            i = nshell(j)*(sum_l+1-j)

            if (i > space_vrr_buf) space_vrr_buf = i
         enddo
         if (space_vrr_buf .eq. 0) space_vrr_buf = 1 !this is just to make sure the automatic arrays in contr_vrr always have an allowed size

         !Compute space needed to store the target class for the ET step: all intermediates whose first angular momentum lies in the range [la,la+lb] and whose third angular momentum lies in [lc,lc+ld].
         !We always compute this since we store the target integrals from the VRR and ET steps in contr_et regardless of whether the ET steps were actually used.
         space_et_tgt = 0
         do i=lc,lc+ld
            do j=la,la+lb
               space_et_tgt = space_et_tgt + nshell(i)*nshell(j)
            enddo
         enddo

         !(re)allocate space for the target contracted ET integrals (if needed)
         err = check_real_array_size(contr_et,space_et_tgt); if (err .ne. 0) stop "contr_et allocation failed."

         !Compute space needed to store the intermediates for the ET step. We compute this only in case the ET steps are really necessary, i.e. if lc + ld > 0.
         if (lc + ld > 0) then
            space_et_buf = 0
            do i=1,lc+ld !|ys]
               k = 0
               l = i
               if (lc+ld > la) l = 0 !in this case (additional) intermediates starting with x=s are required
               do j=l,la+lb+lc+ld-i ![xs|
                  k = k + nshell(i)*nshell(j)
               enddo
               if (k > space_et_buf) space_et_buf = k !k = space required to store the current row [xs|ys] of ET intermediates
            enddo
            !update the size needed for the targets from the VRR step; in this case the VRR targets will be used also as one of the buffers required for the ET step.
            space_vrr_tgt = max(space_vrr_tgt,space_et_buf)
         else
            space_et_buf = 1 !this is just to make sure the automatic arrays in contr_vrr always have an allowed size
         endif
         !By now vrr_tgt has the size needed: 1) to store the targets from the VRR step and 2) to use them in the following ET step as one of the temporary buffers.

         !space for the HRR1 buffers
         space_hrr1_buf = 0
         if ((lb > 0 .and. ld .eq. 0) .or. (lb > 0 .and. ld > 0)) then !simulate the HRR1 outer loops to find out the maximum size of the target intermediates
            !loop over lc,lc+1,...,
            do s_y = lc, lc+ld
               !loop over p,d,...,
               do s_xp1 = 1,lb
                  !compute space for the next target group
                  space = (ncart(la+lb-s_xp1)-ncart(la-1))*nshell(s_y)*nshell(s_xp1)
                  if (space > space_hrr1_buf) space_hrr1_buf = space
               enddo
            enddo
            !space for the HRR1 target integrals
            space_hrr1_tgt = nshell(la)*nshell(lb)*(ncart(lc+ld)-ncart(lc-1))
            err = check_real_array_size(hrr1_tgt,space_hrr1_tgt); if (err .ne. 0) stop "hrr1_tgt allocation failed."
         endif

         !space for the HRR2 buffers
         space_hrr2_buf = 0
         if ((lb .eq. 0 .and. ld > 0) .or. (lb > 0 .and. ld > 0)) then !simulate the HRR2 outer loop to find out the maximum size of the target intermediates

            i = (2*la+1)*(2*lb+1)
            do s_zp1 = 1,ld
               !compute space for the target integrals in this row
               space = i*(ncart(lc+ld-s_zp1)-ncart(lc-1))*nshell(s_zp1)
               if (space > space_hrr2_buf) space_hrr2_buf = space
            enddo

            !in this case we will use hrr1_tgt as a buffer for half spherical integrals obtained from contr_et so hrr1_tgt must be large enough.
            if (space_et_tgt > space_hrr1_tgt) then
               space_hrr1_tgt = space_et_tgt
               err = check_real_array_size(hrr1_tgt,space_hrr1_tgt); if (err .ne. 0) stop "hrr1_tgt allocation 2 failed."
            endif

         endif

         !space and data for the cartesian -> spherical transform
         call allocate_space_sph_transf(la,lb,lc,ld)

   end subroutine allocate_space

   !space and data for the cartesian -> spherical transform
   subroutine allocate_space_sph_transf(la,lb,lc,ld)
      use special_functions, only: cfp_sph_shell_to_cart_lshells
      implicit none
      integer, intent(in) :: la,lb,lc,ld
      integer :: space, l, i, j, err, space_nz, ind
      real(kind=cfp), allocatable :: tmp_c_nz(:)
      integer, allocatable :: tmp_nz_can(:)

         space = (2*la+1)*(2*lb+1)*max(nshell(lc)*nshell(ld), (ncart(lc+ld)-ncart(lc-1))) !the second value is used as crt_l variable in eri_shell
         !temporary arrays used in sh_ab for the intermediate half transformed integrals.
         err = check_real_array_size(half_sph,space); if (err .ne. 0) stop "half_sph allocation failed."

         !temporary array used in sh_cd
         err = check_real_array_size(transf,space); if (err .ne. 0) stop "transf allocation failed."

         l = max(la,lb,lc,ld)
         if (l > max_l) then

            !make sure we save the already calculated coefficients: move them into tmp
            if (allocated(c_nz)) then
               deallocate(nz,nz_can,c_nz,sph_map_start)
            endif

            space_nz = (l+1)*(l+1)
            space = (1+l)*(2+l)*(3+l)*(2+3*l)/12 !space = sum_{i=0}^{l} (2*i+1)*(i+1)*(i+2)/2
            allocate(nz(space_nz),tmp_nz_can(space),tmp_c_nz(space),stat=err)
            if (err .ne. 0) stop "cart_to_sph allocation 2 failed."

            !recalculate the mapping for L=0,...,l
            call cfp_sph_shell_to_cart_lshells(l,nz,tmp_c_nz,tmp_nz_can)

            !How many non-zero coefficients there are in total? nz contains the number of non-zero coefficients for each L,m combination for L=0,...,l.
            space = sum(nz(1:space_nz))

            allocate(nz_can(space),c_nz(space),sph_map_start(space_nz),stat=err)
            if (err .ne. 0) stop "cart_to_sph allocation 3 failed."

            !transfer the calculated coefficients and the mapping data
            nz_can(1:space) = tmp_nz_can(1:space)  !in-shell canonical indices of the cartesians contributing to the cartesian->spherical mapping for each l,m
            c_nz(1:space) = tmp_c_nz(1:space)      !coefficients for the cartesian->spherical mapping

            ind = 0
            do i=0,l
               do j=-i,i
                  ind = ind + 1
                  !where in the c_nz array do the coefficients for the l,m cartesian -> spherical map start
                  sph_map_start(ind) = sum(nz(1:ind-1))+1
               enddo !j
            enddo !i

            deallocate(tmp_nz_can,tmp_c_nz)

            !finally, update the value of the maximum angular momentum encountered
            max_l = l

         endif

   end subroutine allocate_space_sph_transf

   !Calculates the full set of 2-particle integrals for the given data for the shells (ab|cd) and the starting indices of the basis functions in each of the shells. The calculated integrals are placed in
   !the module linear array sph_ints. The basis function indices are expanded and ordered and are placed in the array int_index. The array int_index must have dimensions at least (1:4,no_int),
   !where no_int is the total number of integrals corresponding to the (ab|cd) combination of shells, i.e.: no_int = (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1). The indices in int_index(1:4,i) then correspond to
   !the 2-particle integral sph_ints(i). There are two possibilities how the indices can be permuted on output. If keep_ab_cd_order .true. then each quartet of indices is permuted so it corresponds to the
   !order of the pairs of shells (ab|cd) on input, where the order within the (ab| or |cd) pair is such that a.ge.b, c.ge.d, i.e. not necessarily the same as on input. 
   !If indexing_method .eq. 1 then the integrals are not reordered and standard indexing is used. If indexing_method .eq. 2 then the columns a,b,c,d of sph_ints are reordered so that ap.ge.bp,cp.ge.dp,ap.ge.cp.
   !If keep_ab_cd_order .false. then each quartet of indices is permuted in the order required for computation of the full index of the (ab|cd) integral: 
   !ind_a .ge. ind_b, ind_c .ge. ind_d, ind_a+ind_b .ge. ind_c+ind_d, ind_a .ge. ind_c.
   !Finally, remember that the ordering of the la,lb,lc,ld shells as performed here must be matched by the ordering in the index_2el_drv routine so if the ordering here or in index_2el_drv is changed the
   !ordering in the other routine must be changed as well.
   subroutine eri(lena,xa,ya,za,anorms,la,aexps,acoefs,ind_a, &
                  lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,ind_b, &
                  lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,ind_c, &
                  lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,ind_d, &
                  two_el_column,int_index,keep_ab_cd_order,indexing_method,&
                  do_tails_for_this_quartet,ab_is_continuum,tgt_prop,tgt_pair,rmat_radius,sph_ints)
      implicit none
      integer, intent(in) :: lena, lenb, lenc, lend, la, lb, lc, ld, ind_a, ind_b, ind_c, ind_d, &
                             two_el_column, tgt_pair, indexing_method
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc, xd,yd,zd, rmat_radius
      !intent(in):
      real(kind=cfp), allocatable :: anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), allocatable :: bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), allocatable :: cnorms(:),cexps(:),ccoefs(:)
      real(kind=cfp), allocatable :: dnorms(:),dexps(:),dcoefs(:)
      real(kind=cfp), allocatable :: tgt_prop(:,:)
      logical, intent(in) :: ab_is_continuum,do_tails_for_this_quartet
      !intent(out):
      integer, allocatable :: int_index(:,:)
      logical, intent(in) :: keep_ab_cd_order
      real(kind=cfp), allocatable :: sph_ints(:,:)

      integer :: no_shell

         no_shell = (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1)

         ! Order the shells according to their angular momentum: la .ge .lb, lc .ge. ld, la+lb .ge. lc+ld, 
         ! perform the calculation (eri_shell) and then expand and order the basis function indices (index_2el)
         if (la+lb .ge. lc+ld) then
            if (la .ge. lb) then
               if (lc .ge. ld) then !la,lb,lc,ld
                  !print *,'abcd'
                  call eri_shell(lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                 lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                 lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                 lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,two_el_column,sph_ints)

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lc,ld,la,lb,lena,lenb,aexps,bexps,acoefs,bcoefs, &
                                             anorms,bnorms,rmat_radius,.false.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,la,lb,lc,ld,lenc,lend,cexps,dexps,ccoefs,dcoefs, &
                                             cnorms,dnorms,rmat_radius,.true.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(la,lb,lc,ld,ind_a,ind_b,ind_c,ind_d,int_index,keep_ab_cd_order,.false.)
                  endif
               else !la,lb,ld,lc
                  !print *,'abdc'
                  call eri_shell (lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                  lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                  lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                  lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,two_el_column,sph_ints)

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,ld,lc,la,lb,lena,lenb,aexps,bexps,acoefs,bcoefs,anorms, &
                                             bnorms,rmat_radius,.false.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,la,lb,ld,lc,lend,lenc,dexps,cexps,dcoefs,ccoefs,dnorms, &
                                             cnorms,rmat_radius,.true.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(la,lb,ld,lc,ind_a,ind_b,ind_d,ind_c,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(la,lb,ld,lc,ind_a,ind_b,ind_d,ind_c,int_index,keep_ab_cd_order,.false.)
                  endif
               endif
            else !la < lb
               if (lc .ge. ld) then !lb,la,lc,ld
                  !print *,'bacd'
                  call eri_shell (lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                  lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                  lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                  lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,two_el_column,sph_ints)

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lc,ld,lb,la,lenb,lena,bexps,aexps,bcoefs,acoefs, &
                                             bnorms,anorms,rmat_radius,.false.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lb,la,lc,ld,lenc,lend,cexps,dexps,ccoefs,dcoefs, &
                                             cnorms,dnorms,rmat_radius,.true.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(lb,la,lc,ld,ind_b,ind_a,ind_c,ind_d,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(lb,la,lc,ld,ind_b,ind_a,ind_c,ind_d,int_index,keep_ab_cd_order,.false.)
                  endif
               else !lb,la,ld,lc
                  call eri_shell (lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                  lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                  lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                  lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,two_el_column,sph_ints)
                  !print *,'badc'

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,ld,lc,lb,la,lenb,lena,bexps,aexps,bcoefs,acoefs, &
                                             bnorms,anorms,rmat_radius,.false.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lb,la,ld,lc,lend,lenc,dexps,cexps,dcoefs,ccoefs, &
                                             dnorms,cnorms,rmat_radius,.true.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(lb,la,ld,lc,ind_b,ind_a,ind_d,ind_c,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(lb,la,ld,lc,ind_b,ind_a,ind_d,ind_c,int_index,keep_ab_cd_order,.false.)
                  endif
               endif
            endif
         else !la+lb < lc+ld
            if (la .ge. lb) then
               if (lc .ge. ld) then !lc,ld,la,lb
                  call eri_shell (lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                  lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                  lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                  lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,two_el_column,sph_ints)
                  !print *,'cdab'

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lc,ld,la,lb,lena,lenb,aexps,bexps,acoefs,bcoefs, &
                                             anorms,bnorms,rmat_radius,.true.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,la,lb,lc,ld,lenc,lend,cexps,dexps,ccoefs,dcoefs, &
                                             cnorms,dnorms,rmat_radius,.false.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(lc,ld,la,lb,ind_c,ind_d,ind_a,ind_b,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(lc,ld,la,lb,ind_c,ind_d,ind_a,ind_b,int_index,keep_ab_cd_order,.true.)
                  endif
               else !ld,lc,la,lb
                  call eri_shell (lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                  lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                  lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                  lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,two_el_column,sph_ints)
                  !print *,'dcab'

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,ld,lc,la,lb,lena,lenb,aexps,bexps,acoefs,bcoefs, &
                                             anorms,bnorms,rmat_radius,.true.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,la,lb,ld,lc,lend,lenc,dexps,cexps,dcoefs,ccoefs, &
                                             dnorms,cnorms,rmat_radius,.false.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(ld,lc,la,lb,ind_d,ind_c,ind_a,ind_b,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(ld,lc,la,lb,ind_d,ind_c,ind_a,ind_b,int_index,keep_ab_cd_order,.true.)
                  endif
               endif
            else !la < lb
               if (lc .ge. ld) then !lc,ld,lb,la
                  call eri_shell (lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                  lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                  lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                  lena,xa,ya,za,anorms,la,aexps,acoefs,two_el_column,sph_ints)
                  !print *,'cdba'

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lc,ld,lb,la,lenb,lena,bexps,aexps,bcoefs,acoefs,bnorms, &
                                             anorms,rmat_radius,.true.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lb,la,lc,ld,lenc,lend,cexps,dexps,ccoefs,dcoefs,cnorms, &
                                             dnorms,rmat_radius,.false.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(lc,ld,lb,la,ind_c,ind_d,ind_b,ind_a,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(lc,ld,lb,la,ind_c,ind_d,ind_b,ind_a,int_index,keep_ab_cd_order,.true.)
                  endif
               else !ld,lc,lb,la
                  call eri_shell (lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                  lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                  lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                  lena,xa,ya,za,anorms,la,aexps,acoefs,two_el_column,sph_ints)
                  !print *,'dcba'

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,ld,lc,lb,la,lenb,lena,bexps,aexps,bcoefs,acoefs, &
                                             bnorms,anorms,rmat_radius,.true.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lb,la,ld,lc,lend,lenc,dexps,cexps,dcoefs,ccoefs, &
                                             dnorms,cnorms,rmat_radius,.false.,eri_tail_int)
                     endif
                     !Subtract the tails:
                     sph_ints(1:no_shell,two_el_column) = sph_ints(1:no_shell,two_el_column) - eri_tail_int(1:no_shell)
                  endif

                  if (indexing_method .eq. 2) then
                     call reorder_and_index_2el(ld,lc,lb,la,ind_d,ind_c,ind_b,ind_a,two_el_column,int_index,sph_ints)
                  else
                     call index_2el(ld,lc,lb,la,ind_d,ind_c,ind_b,ind_a,int_index,keep_ab_cd_order,.true.)
                  endif
               endif
            endif
         endif

   end subroutine eri

   !we assume that the shells have been already ordered so that la+lb .ge. lc+ld, la .ge. lb, lc .ge. ld
   subroutine eri_shell (lena,xa,ya,za,anorms,la,aexps,acoefs, &
                         lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                         lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                         lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,two_el_column,sph_ints)
      use const, only: mmax, imax_wp,boys_f_grid_step_wp,taylor_k_wp, imax_ep,boys_f_grid_step_ep,taylor_k_ep
      implicit none
      integer, intent(in) :: lena, lenb, lenc, lend, la, lb, lc, ld, two_el_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc, xd,yd,zd
      real(kind=cfp), allocatable :: sph_ints(:,:)
      !intent(in):
      real(kind=cfp), allocatable :: anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), allocatable :: bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), allocatable :: cnorms(:),cexps(:),ccoefs(:)
      real(kind=cfp), allocatable :: dnorms(:),dexps(:),dcoefs(:)

      integer :: i, crt_shell_c, crt_shell_d, crt_l, crt, sph_shell_a, sph_shell_b
      integer :: space_vrr_tgt, space_et_tgt, space_sph_ints, sum_l, space_vrr_buf, space_et_buf, &
                 space_hrr1_buf, space_hrr2_buf, space_hrr1_tgt, space_hrr2_tgt

         sum_l = la+lb+lc+ld
         call calc_can(sum_l)

         !allocate space for all target and intermediate buffers including the space for the final integrals
         call allocate_space (la, lb, lc, ld, space_vrr_tgt, space_et_tgt, space_sph_ints, space_vrr_buf, &
                              space_et_buf, space_hrr1_buf, space_hrr2_buf, space_hrr1_tgt, space_hrr2_tgt)

         !test for atomic case and apply the selection rule
         if (xb .eq. xa .and. yb .eq. ya .and. zb .eq. za) then !R_A = R_B
            if (xd .eq. xc .and. yd .eq. yc .and. zd .eq. zc) then !R_C = R_D
               if (xa .eq. xc .and. ya .eq. yc .and. za .eq. zc) then !R_A = R_C
                  if (abs(la-lb) > lc+ld) then !Wigner-Eckart assuming la+lb .ge. lc+ld
                     sph_ints(1:(2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1),two_el_column) = 0.0_cfp
                     return
                  endif
               endif
            endif
         endif

         ! initialize the Boys function; if i == 1 then it has been initialized already. Note that we use the appropriate parameters
         ! (found using the test program test_boys_function_obj) depending on the precision used.
         if (cfp .eq. wp) then
            i = boys%init(imax_wp,mmax,boys_f_grid_step_wp,taylor_k_wp)
         elseif (cfp .eq. ep1) then
            i = boys%init(imax_ep,mmax,boys_f_grid_step_ep,taylor_k_ep)
         else
            stop "eri_shell: unsupported numeric type"
         endif

         if (i .ne. 0 .and. i .ne. 1) then
            print *,i
            stop "eri_shell: initialization of the Boys function failed."
         endif

         !VRR and ET steps: generate all 6-index terms needed to generate the 12-index integrals; the result is in contr_et
         call contr_vrr (lena,xa,ya,za,anorms,aexps,acoefs, &
                         lenb,xb,yb,zb,bnorms,bexps,bcoefs, &
                         lenc,xc,yc,zc,cnorms,cexps,ccoefs, &
                         lend,xd,yd,zd,dnorms,dexps,dcoefs, &
                         la,lb,lc,ld,contr_et,space_et_tgt,space_vrr_tgt,space_vrr_buf,space_et_buf)

         sph_shell_a = 2*la+1
         sph_shell_b = 2*lb+1
         crt_shell_c = nshell(lc)
         crt_shell_d = nshell(ld)

         !HRR steps including the cartesian -> spherical harmonics transformation
         if (lb .eq. 0 .and. ld .eq. 0) then !no HRR steps required: the target class from the ET step are the final cartesian integrals; the target class from the ET step will be in contr_et

            if (la > 0 .and. lc .eq. 0) then  !(la s|s s)
               call sh_ab(contr_et,sph_ints(:,two_el_column),la,lb,crt_shell_c,crt_shell_d)
            elseif (la > 0 .and. lc > 0) then !(la s|lc s)
               call sh_ab(contr_et,hrr2_tgt,la,lb,crt_shell_c,crt_shell_d)
               call sh_cd(hrr2_tgt,sph_ints(:,two_el_column),sph_shell_a,sph_shell_b,lc,ld)
            elseif (la .eq. 0 .and. lc .eq. 0) then !(ss|ss)
               sph_ints(1,two_el_column) = contr_et(1)
            endif

         else !(lb .ne. 0 .or. ld .ne. 0): further HRR steps will be required: the target class from the ET step will be in contr_et buffer

            if (lb > 0 .and. ld .eq. 0) then !(la lb|lc s); no HRR2 steps required
               call hrr1(la,xa,ya,za,lb,xb,yb,zb,lc,ld,contr_et,hrr1_tgt,space_hrr1_buf)
               if (lc > 0) then
                  call sh_ab(hrr1_tgt,hrr2_tgt,la,lb,crt_shell_c,crt_shell_d)
                  call sh_cd(hrr2_tgt,sph_ints(:,two_el_column),sph_shell_a,sph_shell_b,lc,ld)
               else  !(la lb|s s)
                  call sh_ab(hrr1_tgt,sph_ints(:,two_el_column),la,lb,crt_shell_c,crt_shell_d)
               endif
            elseif (lb .eq. 0 .and. ld > 0) then !(la s|lc ld); HRR2 steps will be required
               !transform into spherical basis the ab part of all ET target batches (ab|l s), l=lc,...,lc+ld in contr_et
               crt = 1
               crt_l = (ncart(lc+ld)-ncart(lc-1))
               call sh_ab(contr_et,hrr1_tgt,la,lb,crt_l,crt)
               call hrr2(lc,xc,yc,zc,ld,xd,yd,zd,la,lb,hrr1_tgt,hrr2_tgt,space_hrr2_buf)
               call sh_cd(hrr2_tgt,sph_ints(:,two_el_column),sph_shell_a,sph_shell_b,lc,ld)
            elseif (lb > 0 .and. ld > 0) then !(la lb|lc ld); HRR2 steps will be required
               call hrr1(la,xa,ya,za,lb,xb,yb,zb,lc,ld,contr_et,hrr1_tgt,space_hrr1_buf)

               !transform into spherical basis the ab part of all ET target batches (ab|l s), l=lc,...,lc+ld in contr_et
               !note that we put the result back into contr_et; this array is guaranteed to be large enough
               crt = 1
               crt_l = (ncart(lc+ld)-ncart(lc-1))
               call sh_ab(hrr1_tgt,contr_et,la,lb,crt_l,crt)

               call hrr2(lc,xc,yc,zc,ld,xd,yd,zd,la,lb,contr_et,hrr2_tgt,space_hrr2_buf)
               call sh_cd(hrr2_tgt,sph_ints(:,two_el_column),sph_shell_a,sph_shell_b,lc,ld)
            endif

         endif

         !write(*,'(4i)') la,lb,lc,ld

         !In case some of the shells are p-type we need to permute the p-type integrals from the Cartesian order (x,y,z) into the spherical order M=-1,0,1: y,z,x.
         call reorder_p_shells(sph_ints(:,two_el_column),la,lb,lc,ld)

         !do i=1,space_sph_ints
         !   write(*,'(i,e25.15)') i, sph_ints(i)
         !enddo

   end subroutine eri_shell

   subroutine reorder_p_shells(sph_ints,la,lb,lc,ld)
      implicit none
      real(kind=cfp), intent(inout) :: sph_ints(*)
      integer, intent(in) :: la,lb,lc,ld

      integer :: sph_a, sph_b, sph_c, sph_d, i, j, err, n
      
         sph_a = 2*la+1
         sph_b = 2*lb+1
         sph_c = 2*lc+1
         sph_d = 2*ld+1

         j = sph_a*sph_b*sph_c*sph_d

         if (la .eq. 1) then
            err = check_real_array_size(transf,1)
            if (err .ne. 0) stop "transf allocation 1 failed"
            do i=1,j,3
               transf(1) = sph_ints(i)
               sph_ints(i) = sph_ints(i+1)   !y
               sph_ints(i+1) = sph_ints(i+2) !z 
               sph_ints(i+2) = transf(1)     !x
            enddo !i
         endif

         if (lb .eq. 1) then
            err = check_real_array_size(transf,sph_a)
            if (err .ne. 0) stop "transf allocation 2 failed"
            do i=1,j,sph_a*3
               transf(1:sph_a) = sph_ints(i:i+sph_a-1) !a_i,bx a_i+1,bx ... a_sph_a,bx
               sph_ints(i:i+sph_a-1) = sph_ints(i+sph_a:i+2*sph_a-1)           !y
               sph_ints(i+sph_a:i+2*sph_a-1) = sph_ints(i+2*sph_a:i+3*sph_a-1) !z
               sph_ints(i+2*sph_a:i+3*sph_a-1) = transf(1:sph_a)               !x
            enddo
         endif

         if (lc .eq. 1) then
            n = sph_a*sph_b
            err = check_real_array_size(transf,n)
            if (err .ne. 0) stop "transf allocation 2 failed"
            do i=1,j,n*3
               transf(1:n) = sph_ints(i:i+n-1)
               sph_ints(i:i+n-1) = sph_ints(i+n:i+2*n-1)           !y
               sph_ints(i+n:i+2*n-1) = sph_ints(i+2*n:i+3*n-1)     !z
               sph_ints(i+2*n:i+3*n-1) = transf(1:n)               !x
            enddo
         endif

         if (ld .eq. 1) then
            n = sph_a*sph_b*sph_c
            err = check_real_array_size(transf,n)
            if (err .ne. 0) stop "transf allocation 2 failed"
            do i=1,j,n*3
               transf(1:n) = sph_ints(i:i+n-1)
               sph_ints(i:i+n-1) = sph_ints(i+n:i+2*n-1)           !y
               sph_ints(i+n:i+2*n-1) = sph_ints(i+2*n:i+3*n-1)     !z
               sph_ints(i+2*n:i+3*n-1) = transf(1:n)               !x
            enddo
         endif

   end subroutine reorder_p_shells

   subroutine hrr1(la,xa,ya,za,lb,xb,yb,zb,lc,ld,src,tgt,space_hrr1_buf)
      implicit none
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      integer, intent(in) :: la,lb,lc,ld
      !intent(in):
      real(kind=cfp), allocatable :: src(:)
      !intent(out):
      real(kind=cfp), allocatable :: tgt(:)
      integer, intent(in) :: space_hrr1_buf

      integer :: d1_h1,d2_h1 !dimensions in the hrr1_bufA,hrr1_bufB linear arrays which are used to emulate three dimensional arrays
      real(kind=cfp) :: x, y, z, r_ab
      integer :: s_w, can_w, ind, can_x, can_xp1, s_xp1, can_wp1, s_y, off_tgt, off_src, can_y_start, can_y_end
      integer :: d1_src, d2_src, ind_wxy_base, ind_p1_base, ind_base, last
      integer :: ind_wsys_base, d1_et, before_s_xp1, before_s_wp1
      logical :: x_dir, y_dir, z_dir, hrr1_bufA_tgt
      real(kind=cfp) :: hrr1_bufA(space_hrr1_buf), hrr1_bufB(space_hrr1_buf)
 
   !HRR step 1: (w(x+1_i)|ys) = ((w+1_i)x|ys) + Rab_i * (wx|ys)

      x = xa-xb
      y = ya-yb
      z = za-zb

      hrr1_bufA_tgt = .false. !if set to true then hrr1_bufA is the array which will hold the target auxiliary integrals; otherwise it is hrr1_bufB
      last = 0 !last index in tgt: this is the number of integrals trasferred between the HRR step 1 and HRR step 2.

      d1_et = ncart(la+lb)-ncart(la-1)
      ind_wsys_base = -ncart(la-1)-d1_et*ncart(lc-1)
      
      !loop over lc,lc+1,...,
      do s_y = lc, lc+ld
  
         off_tgt = 0
         d1_h1 = nshell(s_y)
         d2_h1 = ncart(la+lb) - ncart(la-1)
         can_y_start = ncart(s_y-1)+1
         can_y_end = ncart(s_y)

         !loop over p,d,...,
         do s_xp1 = 1,lb

            !compute the offset indices for the target and the source buffers
            off_src = off_tgt
            d1_src = d1_h1
            d2_src = d2_h1

            d1_h1 = nshell(s_y)
            d2_h1 = ncart(la+lb-s_xp1)-ncart(la-1)

            off_tgt = (ncart(s_y-1)+1 + d1_h1*ncart(la-1) + d1_h1*d2_h1*ncart(s_xp1-1)) -1

            before_s_xp1 = ncart(s_xp1-2) !how many cartesians there are in all shells up to s=s_xp1-2

            !loop over la,la+1,...,
            do s_w = la,la+lb-s_xp1

               before_s_wp1 = ncart(s_w) !how many cartesians there are in all shells up to s=s_w

               !loop over the targets in the micro group
               do can_xp1 = ncart(s_xp1-1)+1,ncart(s_xp1)

                  x_dir = .false.; y_dir = .false.; z_dir = .false.
   
                  if (cart_l(can_xp1) > 0) then
                     x_dir = .true.
                     can_x = before_s_xp1 + can_shell(s_xp1-1,cart_l(can_xp1)-1,cart_n(can_xp1))
                     r_ab = x
                  elseif (cart_m(can_xp1) > 0) then
                     y_dir = .true.
                     can_x = before_s_xp1 + can_shell(s_xp1-1,cart_l(can_xp1),cart_n(can_xp1))
                     r_ab = y
                  else !z
                     z_dir = .true.
                     can_x = before_s_xp1 + can_shell(s_xp1-1,cart_l(can_xp1),cart_n(can_xp1)-1)
                     r_ab = z
                  endif

                  do can_w = ncart(s_w-1)+1,ncart(s_w)

                     if (x_dir) then
                        can_wp1 = before_s_wp1 + can_shell(s_w+1,cart_l(can_w)+1,cart_n(can_w))
                     elseif (y_dir) then
                        can_wp1 = before_s_wp1 + can_shell(s_w+1,cart_l(can_w),cart_n(can_w))
                     else !z_dir
                        can_wp1 = before_s_wp1 + can_shell(s_w+1,cart_l(can_w),cart_n(can_w)+1)
                     endif

                     ind_wxy_base = d1_src*(can_w-1)  + d1_src*d2_src*(can_x-1) - off_src
                     ind_p1_base = d1_src*(can_wp1-1)+ d1_src*d2_src*(can_x-1) - off_src
                     ind_base = d1_h1*(can_w-1)  + d1_h1*d2_h1*(can_xp1-1) - off_tgt

                     !Vectorized loop over can_y: can_y = ncart(s_y-1)+1,ncart(s_y); the result are the integrals (w(x+1_i)|ys) with w ~ l_w; x+1_i ~ l_xp1; y ~ l_y
                     !The results are stored in buf1(ind) or buf2(ind) depending on which one is the target; ind = can_y + d1_h1*(can_w-1)  + d1_h1*d2_h1*(can_xp1-1)- off_tgt :(w(x+1_i)|ys)
                     call hrr1_micro (can_y_start, can_y_end, s_xp1,d1_et, ind_wsys_base, ind_wxy_base, ind_p1_base, &
                                      ind_base, hrr1_bufA_tgt, hrr1_bufA, hrr1_bufB, src, can_w, can_wp1, r_ab)

                  enddo !can_w

               enddo !can_xp1

            enddo !s_w

            !prepare for the next row of integrals: swap the target buffer with the source one
            hrr1_bufA_tgt = .not.(hrr1_bufA_tgt)

         enddo !s_xp1

         !transfer the target integrals just calculated into the source buffer for the HRR step 2
         ind = ncart(s_y) + ind_base !index of the last element in the target buffer hrr1_bufA or hrr1_bufB
  
         !transfer the values into the tgt buffer in the normal order: la,lb,s_y; the order in hrr1_buf* is s_y,la,lb
         if (hrr1_bufA_tgt) then !hrr1_bufB contains the integrals we want to transfer - see the '.not.(hrr1_bufA_tgt)' above
            call from_hrr1_tgt_to_hrr2_src(la,lb,s_y,hrr1_bufB,tgt,last)
         else
            call from_hrr1_tgt_to_hrr2_src(la,lb,s_y,hrr1_bufA,tgt,last)
         endif

         last = last + ind !increment the number of integrals transferred to tgt
         
      enddo !s_y

   end subroutine hrr1

   subroutine hrr1_micro (can_y_start, can_y_end, s_xp1, stride, ind_wsys_base, ind_wxy_base, ind_p1_base, ind_base, &
                          buf1_tgt, buf1, buf2, et_tgt, can_w, can_wp1, r_ab)
      implicit none
      integer, intent(in) :: can_y_start,can_y_end,s_xp1,stride,ind_wxy_base,ind_p1_base,ind_base,can_w,can_wp1,ind_wsys_base
      logical, intent(in) :: buf1_tgt
      real(kind=cfp), intent(in) :: r_ab
      !intent(in):
      real(kind=cfp), allocatable :: et_tgt(:)
      real(kind=cfp), intent(inout) :: buf1(*), buf2(*)

         if (s_xp1 .eq. 1) then !the source values are in et_tgt which is the class of the target auxiliary integrals from the ET step
            if (buf1_tgt) then !buf1 contains the target auxiliary integrals
               call hrr1_micro_xp1_p(can_y_start,can_y_end,stride,ind_base,ind_wsys_base,buf1,et_tgt,can_w,can_wp1,r_ab)
            else !buf2 contains the target auxiliary integrals
               call hrr1_micro_xp1_p(can_y_start,can_y_end,stride,ind_base,ind_wsys_base,buf2,et_tgt,can_w,can_wp1,r_ab)
            endif
         else !one of buf1 or buf2 contains the source integrals and buf2 or buf1 is the target
            if (buf1_tgt) then !buf1 contains the target auxiliary integrals
               call hrr1_micro_xp1_general(can_y_start,can_y_end,ind_wxy_base,ind_p1_base,ind_base,buf1,buf2,r_ab)
            else !buf2 contains the target auxiliary integrals
               call hrr1_micro_xp1_general(can_y_start,can_y_end,ind_wxy_base,ind_p1_base,ind_base,buf2,buf1,r_ab)
            endif
         endif

   end subroutine hrr1_micro

   subroutine hrr1_micro_xp1_p(can_y_start,can_y_end,stride,ind_base,ind_wsys_base,tgt,src,can_w,can_wp1,r_ab)
      implicit none
      integer, intent(in) :: can_y_start,can_y_end,stride,ind_base,can_w,can_wp1,ind_wsys_base
      real(kind=cfp) :: r_ab
      real(kind=cfp), intent(in) :: src(*)
      real(kind=cfp), intent(out) :: tgt(*)

      integer :: can_y, src1_base, src2_base !, ind_wxy, ind_p1, ind

         src1_base = can_wp1 + ind_wsys_base
         src2_base = can_w + ind_wsys_base
         forall (can_y = can_y_start:can_y_end)                  !ncart(s_y-1)+1:ncart(s_y)
            !ind = can_y + ind_base                              ! = can_y + d1_h1*(can_w-1)  + d1_h1*d2_h1*(can_xp1-1)- off_tgt !(w(x+1_i)|ys)
            !ind_wxy = can_w + stride*(can_y-1) + ind_wsys_base  ! = (ws|ys)
            !ind_p1 = can_wp1 + stride*(can_y-1)+ ind_wsys_base  ! = ((w+1)s|ys)
            !tgt(ind) = src(ind_p1) + r_ab*src(ind_wxy)
            tgt(can_y + ind_base) = src(src1_base + stride*(can_y-1)) + r_ab*src(src2_base + stride*(can_y-1))
         endforall

   end subroutine hrr1_micro_xp1_p

   subroutine hrr1_micro_xp1_general(can_y_start,can_y_end,ind_wxy_base,ind_p1_base,ind_base,tgt,src,r_ab)
      implicit none
      integer, intent(in) :: can_y_start,can_y_end,ind_wxy_base,ind_p1_base,ind_base
      real(kind=cfp) :: r_ab
      real(kind=cfp), intent(in) :: src(*)
      real(kind=cfp), intent(out) :: tgt(*)

      integer :: can_y !, ind_wxy, ind_p1, ind

         forall (can_y = can_y_start:can_y_end)          !ncart(s_y-1)+1:ncart(s_y)
            !ind_wxy = can_y + ind_wxy_base              ! = can_y + d1_src*(can_w-1)  + d1_src*d2_src*(can_x-1)  - off_src !(wx|ys)
            !ind_p1  = can_y + ind_p1_base               ! = can_y + d1_src*(can_wp1-1)+ d1_src*d2_src*(can_x-1)  - off_src !((w+1_i)x|ys)
            !ind     = can_y + ind_base                  ! = can_y + d1_h1*(can_w-1)  + d1_h1*d2_h1*(can_xp1-1)- off_tgt !(w(x+1_i)|ys)
            !tgt(ind) = src(ind_p1) + r_ab*src(ind_wxy)
            tgt(can_y + ind_base) = src(can_y + ind_p1_base) + r_ab*src(can_y + ind_wxy_base)
         endforall !can_y

   end subroutine hrr1_micro_xp1_general

   !transfers one set of elementary targets of the type (la lb|y s) from the HRR1 step into the buffer tgt
   !which is assumed to be the source buffer for the HRR2 step, i.e. (la lb|y s) -> (y s|la lb).
   !The value last is the last used index in the buffer tgt. This routine is equivalent to transposition
   !of a matrix with dimension nshell(la)*nshell(lb),nshell(s_y).
   subroutine from_hrr1_tgt_to_hrr2_src(la,lb,s_y,src,tgt,last)
      use const, only: tile
      implicit none
      real(kind=cfp), intent(in) :: src(*)
      !intent(out):
      real(kind=cfp), allocatable :: tgt(:)
      integer, intent(in) :: last,la,lb,s_y

      integer(kind=shortint) :: i, j, ii, jj, nrow, ncol, col, iend, jend
      integer(kind=shortint) :: tile_m_1 = tile - 1

         !The target buffer is effectively a matrix with dimensions: nrow,ncol.
         nrow = nshell(s_y)
         ncol = nshell(la)*nshell(lb)
         do jj = 1,nrow,tile
            jend = min(nrow, jj + tile_m_1)
            do ii = 1,ncol,tile
               iend = min(ncol, ii + tile_m_1)

               do j = jj,jend
                  col = last+ncol*(j-1)
                  do i = ii,iend
                     tgt(i+col) = src(j+nrow*(i-1))
                     !The order in the target buffer is: a,b,y; the order in the source buffer is: y,a,b
                     !i = last + ind_a + in_shell_la*(ind_b-1) + in_shell_la*in_shell_lb*(ind_y-1) !the order in the target buffer is: a,b,y 
                     !j = ind_y + in_shell_s_y*(ind_a-1) + in_shell_s_y*in_shell_la*(ind_b-1)      !the order in the source buffer is: y,a,b
                     !tgt(i) = src(j)
                  enddo
               enddo

            enddo
         enddo

   end subroutine from_hrr1_tgt_to_hrr2_src

   subroutine hrr2(lc,xc,yc,zc,ld,xd,yd,zd,la,lb,src,tgt,space_hrr2_buf)
      implicit none
      real(kind=cfp), intent(in) :: xc,yc,zc, xd,yd,zd
      integer, intent(in) :: la,lb,lc,ld
      !intent(in):
      real(kind=cfp), allocatable :: src(:)
      !intent(out):
      real(kind=cfp), allocatable :: tgt(:)
      integer, intent(in) :: space_hrr2_buf

      integer :: d1,d2 !dimensions in the hrr1_bufA,hrr1_bufB, hrr2_bufA,hrr2_bufB linear arrays which are used to emulate three and four dimensional arrays
      real(kind=cfp) :: x, y, z, r_cd
      integer :: ind_y, s_y, ind_zp1, ind_z, ind_yp1, s_zp1, off_tgt
      integer :: ind_base, off_src1, off_src2, can_zp1, can_y
      integer :: ind_wxyz_base, ind_wxyp1z_base, in_shell_s_y, in_shell_s_yp1, in_shell_z, in_shell_zp1, in_shells_ab, in_shell_lc
      logical :: x_dir, y_dir, z_dir, hrr2_bufA_tgt
      real(kind=cfp) :: hrr2_bufA(space_hrr2_buf), hrr2_bufB(space_hrr2_buf)
 
   !HRR step 2: (wx|y(z+1_i)) = (wx|(y+1_i)z) + Rcd_i * (wx|yz)

      hrr2_bufA_tgt = .false. !if set to true then hrr2_bufA is the buffer containing the target auxiliary integrals; otherwise it is hrr2_bufB

      x = xc-xd
      y = yc-yd
      z = zc-zd

      !note that the nshell() values would apply if the ab part of the source integrals (la lb|lc s) was cartesian instead of spherical.
      d1 = 2*la+1 !nshell(la)
      d2 = 2*lb+1 !nshell(lb)
      in_shells_ab = d1*d2

      in_shell_lc = nshell(lc)

      do s_zp1 = 1,ld

         in_shell_zp1 = nshell(s_zp1)
         in_shell_z = nshell(s_zp1-1)

        !offset for the integrals in the target and source batches for the target row y=s_y 
        off_tgt = 0
        off_src1 = in_shells_ab*in_shell_lc*in_shell_z
        off_src2 = 0

        do s_y = lc,lc+ld-s_zp1

           in_shell_s_y = nshell(s_y)
           in_shell_s_yp1 = nshell(s_y+1)
  
           do ind_zp1 = 1,in_shell_zp1

              can_zp1 = ind_zp1 + ncart(s_zp1-1)

              x_dir = .false.; y_dir = .false.; z_dir = .false.
              if (cart_l(can_zp1) > 0) then
                 ind_z = can_shell(s_zp1-1,cart_l(can_zp1)-1,cart_n(can_zp1))
                 x_dir = .true.
                 r_cd = x
              elseif (cart_m(can_zp1) > 0) then
                 ind_z = can_shell(s_zp1-1,cart_l(can_zp1),cart_n(can_zp1))
                 y_dir = .true.
                 r_cd = y
              else !z
                 ind_z = can_shell(s_zp1-1,cart_l(can_zp1),cart_n(can_zp1)-1)
                 z_dir = .true.
                 r_cd = z
              endif
   
              do ind_y = 1,in_shell_s_y

                 can_y = ind_y + ncart(s_y-1)
      
                 if (x_dir) then
                    ind_yp1 = can_shell(s_y+1,cart_l(can_y)+1,cart_n(can_y))
                 elseif (y_dir) then
                    ind_yp1 = can_shell(s_y+1,cart_l(can_y),cart_n(can_y))
                 else !z_dir
                    ind_yp1 = can_shell(s_y+1,cart_l(can_y),cart_n(can_y)+1)
                 endif

                 ind_base = in_shells_ab*(ind_y-1) + in_shells_ab*in_shell_s_y*(ind_zp1-1) + off_tgt !(wx|y(z+1_i))
                 ind_wxyz_base = in_shells_ab*(ind_y-1) + in_shells_ab*in_shell_s_y*(ind_z-1) + off_src2 !(wx|yz)
                 ind_wxyp1z_base = in_shells_ab*(ind_yp1-1) + in_shells_ab*in_shell_s_yp1*(ind_z-1) + off_src1 !(wx|(y+1_i)z)

                 call hrr2_micro (ld, s_zp1, in_shells_ab, r_cd, ind_base, ind_wxyp1z_base, ind_wxyz_base, hrr2_bufA_tgt, &
                                  hrr2_bufA, hrr2_bufB, src, tgt)

              enddo !ind_y

           enddo !ind_zp1

           off_tgt = off_tgt + in_shells_ab*in_shell_s_y*in_shell_zp1
           off_src2= off_src2+ in_shells_ab*in_shell_s_y*in_shell_z
           off_src1= off_src1+ in_shells_ab*in_shell_s_yp1*in_shell_z
   
        enddo !s_y

        !prepare for the next row of integrals: swap the target buffer with the source one
        hrr2_bufA_tgt = .not.(hrr2_bufA_tgt)

      enddo !s_zp1

   end subroutine hrr2

   subroutine hrr2_micro (ld, s_zp1, in_shells_ab, r_cd, ind_base, ind_wxyp1z_base, ind_wxyz_base, &
                          hrr2_bufA_tgt, hrr2_bufA, hrr2_bufB, src, tgt)
      implicit none
      integer, intent(in) :: ld,s_zp1,in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base
      logical, intent(in) :: hrr2_bufA_tgt
      real(kind=cfp), intent(in) :: r_cd
      !intent(in):
      real(kind=cfp), allocatable :: src(:)
      real(kind=cfp), intent(inout) :: hrr2_bufA(*), hrr2_bufB(*)
      !intent(out):
      real(kind=cfp), allocatable :: tgt(:)

         if (s_zp1 .eq. 1 .and. s_zp1 .eq. ld) then
            call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,src,tgt)
            return
         endif

         if (hrr2_bufA_tgt) then
            if (s_zp1 .eq. 1 .and. s_zp1 .ne. ld) then
               call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,src,hrr2_bufA)
            elseif (s_zp1 .eq. ld) then
               call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,hrr2_bufB,tgt)
            else
               call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,hrr2_bufB,hrr2_bufA)
            endif
         else
            if (s_zp1 .eq. 1 .and. s_zp1 .ne. ld) then
               call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,src,hrr2_bufB)
            elseif (s_zp1 .eq. ld) then
               call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,hrr2_bufA,tgt)
            else
               call hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,hrr2_bufA,hrr2_bufB)
            endif
         endif

   end subroutine hrr2_micro

   !HRR2 step in the micro group (wx|y(z+1_i)): for fixed y,z+1_i and all wx.
   subroutine hrr2_micro_zp1_general(in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base,r_cd,src,tgt)
      implicit none
      integer, intent(in) :: in_shells_ab,ind_base,ind_wxyp1z_base,ind_wxyz_base
      real(kind=cfp), intent(in) :: r_cd
      real(kind=cfp), intent(in) :: src(*)
      real(kind=cfp), intent(out) :: tgt(*)

      integer :: wx

         forall (wx=1:in_shells_ab)
            tgt(wx+ind_base) = src(wx+ind_wxyp1z_base) + r_cd*src(wx+ind_wxyz_base)
         endforall

   end subroutine hrr2_micro_zp1_general

   !todo change the module not to use any *norms arrays: absorb all normalizations into the contraction coefficients. This saves compute time.
   subroutine contr_vrr (lena,xa,ya,za,anorms,aexps,acoefs, &
                         lenb,xb,yb,zb,bnorms,bexps,bcoefs, &
                         lenc,xc,yc,zc,cnorms,cexps,ccoefs, &
                         lend,xd,yd,zd,dnorms,dexps,dcoefs, &
                         la, lb, lc, ld, contr_et_tgt, size_contr_et, size_vrr_tgt, size_vrr_buff, size_et_buff)
      implicit none
      integer, intent(in) :: lena, lenb, lenc, lend, la,lb,lc,ld, size_contr_et, size_vrr_tgt, size_vrr_buff, size_et_buff
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc, xd,yd,zd
      !intent(in):
      real(kind=cfp), allocatable :: anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), allocatable :: bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), allocatable :: cnorms(:),cexps(:),ccoefs(:)
      real(kind=cfp), allocatable :: dnorms(:),dexps(:),dcoefs(:)
      !intent(out):
      real(kind=cfp), allocatable :: contr_et_tgt(:)

      real(kind=cfp) :: Fm(la+lb+lc+ld+1),vrr_buf1(size_vrr_buff), vrr_buf2(size_vrr_buff), &
                                          vrr_buf3(size_vrr_buff), vrr_tgt(size_vrr_tgt)
      real(kind=cfp) :: et_buf2(size_et_buff), et_buf3(size_et_buff), et_tgt(size_contr_et)

      integer :: i,j,k,l,ind
      real(kind=cfp) :: prod, prod_ij, prod_i, prod_ijk, rab2, rcd2

      rab2 = dist2(xa,ya,za,xb,yb,zb)
      rcd2 = dist2(xc,yc,zc,xd,yd,zd)

      !make sure only the required elements are zeroed out
      contr_et_tgt(1:size_contr_et) = 0.0_cfp

      !todo simplify these loops in case some shells are equal
      do i=1,lena
         prod_i = acoefs(i)*anorms(i)
         do j=1,lenb
            prod_ij = prod_i*bcoefs(j)*bnorms(j)
            do k=1,lenc
               prod_ijk = prod_ij*ccoefs(k)*cnorms(k)
               do l=1,lend

                  prod = prod_ijk*dcoefs(l)*dnorms(l) !product of all contraction coefficients and norms

                  !calculate the vrr and et terms for this primitive combination of shells
                  call vrr_et (xa,ya,za,aexps(i), &
                               xb,yb,zb,bexps(j), &
                               xc,yc,zc,cexps(k), &
                               xd,yd,zd,dexps(l), &
                               la,lb,lc,ld, rab2,rcd2, Fm,vrr_buf1,vrr_buf2,vrr_buf3,vrr_tgt, et_buf2,et_buf3, et_tgt)

                  !accummulate data for all shells into the contr_et_tgt; vectorized loop
                  forall (ind=1:size_contr_et)
                     contr_et_tgt(ind) = contr_et_tgt(ind) + prod*et_tgt(ind)
                  endforall

               enddo
            enddo
         enddo
      enddo

   end subroutine contr_vrr

   !VRR and ET steps of the HGP algorithm
   subroutine vrr_et (xa,ya,za,alphaa, &
                      xb,yb,zb,alphab, &
                      xc,yc,zc,alphac, &
                      xd,yd,zd,alphad, &
                      la,lb,lc,ld, rab2,rcd2, Fm,vrr_buf1,vrr_buf2,vrr_buf3,vrr_tgt, et_buf2,et_buf3, et_tgt)
      implicit none
      real(kind=cfp), parameter :: fac = rtwo*pi54
      real(kind=cfp), intent(in) :: xa,ya,za,alphaa, xb,yb,zb,alphab, xc,yc,zc,alphac, xd,yd,zd,alphad, rab2,rcd2 
      integer, intent(in) :: la,lb,lc,ld
      real(kind=cfp), intent(out) :: Fm(*), vrr_buf1(*),vrr_buf2(*),vrr_buf3(*),vrr_tgt(*), et_buf2(*),et_buf3(*), et_tgt(*)

      real(kind=cfp) :: px,py,pz,qx,qy,qz,zeta,eta,wx,wy,wz,Kcd,rpq2,T,Kab,two_zeta,two_eta,e_o_ez,z_o_ez,two_ze,&
                            deltax,deltay,deltaz,alp, wpx,wpy,wpz, pax,pay,paz, hlp
      real(kind=cfp) :: a_xa, a_ya, a_za, b_xb, b_yb, b_zb, c_xc, c_yc, c_zc, d_xd, d_yd, d_zd, inv_apb, inv_cpd, &
                            inv_zeta_p_eta, alp_zeta, zeta_inv_zeta_p_eta, eta_inv_zeta_p_eta
      integer :: s,im,m_max, no_cart_vrr

         a_xa = alphaa*xa; a_ya = alphaa*ya; a_za = alphaa*za;
         b_xb = alphab*xb; b_yb = alphab*yb; b_zb = alphab*zb;
         c_xc = alphac*xc; c_yc = alphac*yc; c_zc = alphac*zc;
         d_xd = alphad*xd; d_yd = alphad*yd; d_zd = alphad*zd;

         zeta = alphaa+alphab
         eta  = alphac+alphad

         inv_apb = 1.0_cfp/zeta
         inv_cpd = 1.0_cfp/eta
         inv_zeta_p_eta = 1.0_cfp/(zeta+eta)

         two_zeta = 0.5_cfp*inv_apb
         two_eta =  0.5_cfp*inv_cpd
         two_ze = 0.5_cfp*inv_zeta_p_eta
         e_o_ez = eta*inv_zeta_p_eta
         z_o_ez = zeta*inv_zeta_p_eta

         alp = zeta*eta*inv_zeta_p_eta
         alp_zeta = alp*inv_apb

         zeta_inv_zeta_p_eta = zeta*inv_zeta_p_eta
         eta_inv_zeta_p_eta = eta*inv_zeta_p_eta

         px = (a_xa + b_xb)*inv_apb
         py = (a_ya + b_yb)*inv_apb 
         pz = (a_za + b_zb)*inv_apb 

         qx = (c_xc + d_xd)*inv_cpd 
         qy = (c_yc + d_yd)*inv_cpd
         qz = (c_zc + d_zd)*inv_cpd

         wx = zeta_inv_zeta_p_eta*px + eta_inv_zeta_p_eta*qx !(zeta*px+eta*qx)*inv_zeta_p_eta
         wy = zeta_inv_zeta_p_eta*py + eta_inv_zeta_p_eta*qy !(zeta*py+eta*qy)*inv_zeta_p_eta
         wz = zeta_inv_zeta_p_eta*pz + eta_inv_zeta_p_eta*qz !(zeta*pz+eta*qz)*inv_zeta_p_eta

         wpx = alp_zeta*(px-qx)
         wpy = alp_zeta*(py-qy)
         wpz = alp_zeta*(pz-qz)

         pax = px-xa
         pay = py-ya
         paz = pz-za

         Kab = fac*inv_apb !*exp(-alphaa*alphab/(alphaa+alphab)*rab2)
         Kcd = fac*inv_cpd !*exp(-alphac*alphad/(alphac+alphad)*rcd2)
         rpq2 = dist2(px,py,pz,qx,qy,qz)
         T = alp*rpq2

         deltax = -(alphab*(xa-xb)+alphad*(xc-xd))*inv_cpd
         deltay = -(alphab*(ya-yb)+alphad*(yc-yd))*inv_cpd
         deltaz = -(alphab*(za-zb)+alphad*(zc-zd))*inv_cpd
 
         m_max = la+lb+lc+ld

         no_cart_vrr = ncart(m_max)

         !calculate Fm(T) for m=0,...,m_max using the routine based on the Taylor expansion
         call boys%eval_taylor(Fm,m_max+1,T,m_max)
         !Fm(1:m_max+1) = boys_function(T,m_max)

         hlp = Kab*Kcd/sqrt(zeta+eta)*exp(-alphaa*alphab*inv_apb*rab2-alphac*alphad*inv_cpd*rcd2)
         forall (im=0:m_max)
            vrr_buf1(im+1) = hlp*Fm(im+1)
         endforall

         !the VRR target integral [ss|ss]^(0)
         vrr_tgt(1) = vrr_buf1(1)

         if (m_max .eq. 0) then !la+lb+lc+ld .eq. 0: nothing else to do in this subroutine; this completes the integral calculation; the final integral is in et_tgt(1)
            et_tgt(1) = vrr_buf1(1)
            return
         endif
!
!-------- We get here only if la+lb+lc+ld > 0
!
!VRR
!
         !generate the [ps|ss]^(m) terms from the terms [ss|ss]^(m); put the target integrals into vrr_tgt
         call vrr_psss(m_max,wpx,wpy,wpz,pax,pay,paz,vrr_buf1,vrr_buf2,vrr_tgt)

         if (la .eq. 1 .and. lb+lc+ld .eq. 0) then !integral calculation has been completed; the final integrals are in vrr_tgt(2:4), so we transfer them into et_tgt which is used in contr_vrr
            et_tgt(1:3) = vrr_tgt(2:4)
            return
         endif

         !generate the triangle of [xs|ss]^(m) terms containing x=d,f,...,la+lb+lc+ld from the auxiliaries in vrr_buf1 and vrr_buf2; put the targets into vrr_tgt
         call vrr_xsss(m_max,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez,vrr_buf1,vrr_buf2,vrr_buf3,vrr_tgt)

         if (lc+ld .eq. 0) then !no ET steps will be required; the target integrals are in vrr_tgt
            s = ncart(la-1) !only the integrals in the range [xs|ss] with x = [la,la+lb] are required for the next steps of the algorithm
            et_tgt(1:(no_cart_vrr-s)) = vrr_tgt(s+1:no_cart_vrr)
            return
         endif
!
!-------- We get here only if lc+ld > 0
!
!ET
!
         !Use the target integrals from the VRR step as one of the source/buffers for the ET step:
         !generate the [xs|ys] terms from the [xs|ss] terms in vrr_tgt
         call et_xsys(m_max,la,lb,lc,ld,deltax,deltay,deltaz,zeta,eta,two_eta,vrr_tgt,et_buf2,et_buf3,et_tgt)

   end subroutine vrr_et

   ![ps|ss]^(m), m=0,m_max-1 generated from [ss|ss]^(m), m=0,m_max
   subroutine vrr_psss(m_max,wpx,wpy,wpz,pax,pay,paz,aux1,aux2,tgt)
      implicit none
      integer, intent(in) :: m_max
      real(kind=cfp), intent(in) :: wpx,wpy,wpz, pax,pay,paz
      real(kind=cfp), intent(in) :: aux1(*)
      real(kind=cfp), intent(out) :: aux2(*), tgt(*)

      integer :: m

         forall (m=0:m_max-1)
            aux2(m*3+1) = pax*aux1(m+1) - wpx*aux1(m+2)
            aux2(m*3+2) = pay*aux1(m+1) - wpy*aux1(m+2) !2 in (m*3+2) stands for the canonical index of py within the p shell. Similar for the other two rows
            aux2(m*3+3) = paz*aux1(m+1) - wpz*aux1(m+2)
         endforall

         !target integrals for the VRR step: px,py,pz
         tgt(2:4) = aux2(1:3)

   end subroutine vrr_psss

   ![xs|ss]^(m) terms starting with x=d from the auxiliaries in src; put the target integrals into tgt
   subroutine vrr_xsss(m_max,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez,aux1,aux2,aux3,tgt)
      implicit none
      integer, intent(in) :: m_max
      real(kind=cfp), intent(in) :: wpx,wpy,wpz, pax,pay,paz, two_zeta, e_o_ez
      real(kind=cfp), intent(inout) :: aux1(*), aux2(*), aux3(*)
      real(kind=cfp), intent(out) :: tgt(*)

      integer :: shell, in_shell, ind
      integer(kind=shortint) :: src1,src2

         !src1, src2 tell me which one of aux1,aux2,aux3 are the source buffers corresponding to the shells s-2, s-1; the remaining buffer corresponds to the target buffer (for shell s).
         !the buffers src1, src2 shift in each iteration
         src1 = 1
         src2 = 2
         ind = 4 !there are 4 preceeding target integrals in tgt: s,px,py,pz
         do shell=2,m_max

            in_shell = nshell(shell)
           
            !calculate the whole row of auxiliaries for the next shell 
            if (src1 .eq. 1 .and. src2 .eq. 2) then
               call xsss(m_max,shell,aux1,aux2,aux3,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez)
               tgt(ind+1:ind+in_shell) = aux3(1:in_shell) !transfer the target integrals [xs|ss]^(0) from the buffer
               src1 = 2; src2 = 3
            elseif (src1 .eq. 2 .and. src2 .eq. 3) then
               call xsss(m_max,shell,aux2,aux3,aux1,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez)
               tgt(ind+1:ind+in_shell) = aux1(1:in_shell) !transfer the target integrals [xs|ss]^(0) from the buffer
               src1 = 3; src2 = 1
            else !(src1 .eq. 3 .and. src2 .eq. 1)
               call xsss(m_max,shell,aux3,aux1,aux2,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez)
               tgt(ind+1:ind+in_shell) = aux2(1:in_shell) !transfer the target integrals [xs|ss]^(0) from the buffer
               src1 = 1; src2 = 2
            endif

            ind = ind + in_shell

         enddo

   end subroutine vrr_xsss

   !calculates one row of auxiliaries [xs|ss]^(m), m=0,shell-m_max, from the VRR step; the source buffers are aux1, aux2. The buffer of the target auxiliaries is aux3. 
   subroutine xsss(m_max,shell,aux1,aux2,aux3,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez)
      implicit none
      integer, intent(in) :: m_max, shell
      real(kind=cfp), intent(in) :: wpx,wpy,wpz, pax,pay,paz, two_zeta, e_o_ez
      real(kind=cfp), intent(inout) :: aux1(*), aux2(*), aux3(*)

      integer :: m, n, in_shell_xp1, in_shell_x, in_shell_xm1, before_shell
      integer :: can_xp1, ind_xp1, ind_x, ind_xm1
      real(kind=cfp) :: p, w

         in_shell_xp1 = nshell(shell)    !number of functions in the target shell
         in_shell_x = nshell(shell-1)    !number of functions in the target source 2 shell
         in_shell_xm1 = nshell(shell-2)  !number of functions in the target source 1 shell

         before_shell = ncart(shell-1)   !how many target integrals, i.e. [xs|ss]^(0) have been determined already

         do ind_xp1=1,in_shell_xp1 !canonical index (corresponding to [x+1_i]) within the target shell

            can_xp1 = before_shell+ind_xp1 !the full canonical index of [x+1_i]

            !determine which direction to take from [x+1_i] and determine canonical indices of the source integrals
            if (cart_l(can_xp1) > 0) then !X
               n = cart_l(can_xp1)-1
               ind_x = ind_xp1 !=can_shell(shell-1,n,cart_n(can_xp1)) ![x]
               p = pax; w = wpx

               if (n > 0) ind_xm1 = ind_xp1 !=can_shell(shell-2,n-1,cart_n(can_xp1)) ![x-1_i]
            elseif (cart_m(can_xp1) > 0) then !Y
               n = cart_m(can_xp1)-1
               ind_x = ind_xp1-shell+cart_l(can_xp1) !=can_shell(shell-1,cart_l(can_xp1),cart_n(can_xp1)) ![x]
               p = pay; w = wpy

               if (n > 0) ind_xm1 = ind_x-(shell-1)+cart_l(can_xp1) !=can_shell(shell-2,cart_l(can_xp1),cart_n(can_xp1)) ![x-1_i]
            else !Z
               n = cart_n(can_xp1)-1
               ind_x = ind_xp1-1-shell+cart_l(can_xp1) !=can_shell(shell-1,cart_l(can_xp1),n) ![x]
               p = paz; w = wpz

               if (n > 0) ind_xm1 = ind_x-1-(shell-1)+cart_l(can_xp1) !=can_shell(shell-2,cart_l(can_xp1),n-1) ![x-1_i]
            endif

            !calculate the next row of auxiliaries using the VRR depending on whether the second part of the recursion relation contributes or not
            !todo this is not very efficient due to the memory strides: switch the loops over m and x. I.e. generate the full shell of x for each m.
            if (n > 0) then
               forall (m =0:m_max-shell)
                  aux3(ind_xp1+in_shell_xp1*m) = p*aux2(ind_x+in_shell_x*m) - w*aux2(ind_x+in_shell_x*(m+1)) &
                        + n*two_zeta*(aux1(ind_xm1+in_shell_xm1*m) - e_o_ez*aux1(ind_xm1+in_shell_xm1*(m+1)))
               endforall
            else
               forall (m =0:m_max-shell)
                  aux3(ind_xp1+in_shell_xp1*m) = p*aux2(ind_x+in_shell_x*m) - w*aux2(ind_x+in_shell_x*(m+1))
               endforall
            endif

         enddo

   end subroutine xsss

   !calculates [xs|ys]. The source buffer is et_buf1 from the VRR step.
   subroutine et_xsys(m_max,la,lb,lc,ld,deltax,deltay,deltaz,zeta,eta,two_eta,et_buf1,et_buf2,et_buf3,et_tgt)
      implicit none
      integer, intent(in) :: m_max, la, lb, lc, ld
      real(kind=cfp), intent(in) :: deltax,deltay,deltaz, zeta, eta, two_eta
      real(kind=cfp), intent(inout) :: et_buf1(*), et_buf2(*), et_buf3(*), et_tgt(*)

      real(kind=cfp) :: alp_ab_cd, delta
      integer :: i, ind_yp1, s_x, before_s, stride1, stride2, stride, ind, in_shell, &
                    n_y, y, s_yp1, yp1, off_y, off_ym1, before_s_x, can_yp1
      integer :: off_yp1, ym1, before_sm1, before_s_xp1, before_s_xm1, col_y, col_yp1, col_ym1
      integer(kind=shortint) :: src1, src2
      logical :: x_dir, y_dir, z_dir, tgt_y, tgt_x, more_aux

         if (lc + ld > la) then !in this case additional auxiliaries are required for the HRR1 step
            more_aux = .true.
         else
            more_aux = .false.
         endif

         tgt_y = .false.

         alp_ab_cd = zeta/eta
         ind = 0 !last used element in the array et_tgt

         !src1, src2 tell me which one of et_buf1,et_buf2,et_buf3 are the source buffers corresponding to the shells s-2, s-1;
         !the remaining buffer corresponds to the target buffer (for shell s).
         !the buffers et_buf1, et_buf2 shift in each iteration
         src1 = 1
         src2 = 2

         do s_yp1 = 1, lc+ld

            if (s_yp1 .ge. lc) tgt_y = .true. !the |y] is in the range of the target ET auxiliaries

            before_s = ncart(s_yp1-1)
            before_sm1 = ncart(s_yp1-2)

            if (more_aux) then !in this case the HRR step requires additional auxiliaries
               i = 0
               off_yp1 = 0
               off_y = 0
               off_ym1 = 0
            else
               i = s_yp1
               off_yp1 = ncart(s_yp1   -1)
               off_y =   ncart(s_yp1-1 -1)
               off_ym1 = ncart(s_yp1-2 -1)
            endif

            stride = ncart(m_max-s_yp1)-off_yp1 !how many [xs| there are in the y+1 row
            stride1 = ncart(m_max-s_yp1+1)-off_y !how many [xs| there are in the y row
            stride2 = ncart(m_max-s_yp1+2)-off_ym1 !how many [xs| there are in the y-1 row

            do yp1 = 1,nshell(s_yp1) !|(y+1)s]
 
                x_dir = .false.; y_dir = .false.; z_dir = .false.
                can_yp1 = yp1+before_s
    
                !canonical index of [y] and [y-1]
                if (cart_l(can_yp1) > 0) then !X
                   x_dir = .true.
                   delta = deltax
                   y = can_shell(s_yp1-1,cart_l(can_yp1)-1,cart_n(can_yp1))
                   n_y = cart_l(y+before_sm1)
                   if (n_y > 0) ym1 = can_shell(s_yp1-2,cart_l(can_yp1)-2,cart_n(can_yp1))
                elseif (cart_m(can_yp1) > 0) then !Y
                   y_dir = .true.
                   delta = deltay
                   y = can_shell(s_yp1-1,cart_l(can_yp1),cart_n(can_yp1))
                   n_y = cart_m(y+before_sm1)
                   if (n_y > 0) ym1 = can_shell(s_yp1-2,cart_l(can_yp1),cart_n(can_yp1))
                else !Z
                   z_dir = .true.
                   delta = deltaz
                   y = can_shell(s_yp1-1,cart_l(can_yp1),cart_n(can_yp1)-1)
                   n_y = cart_n(y+before_sm1)
                   if (n_y > 0) ym1 = can_shell(s_yp1-2,cart_l(can_yp1),cart_n(can_yp1)-2)
                endif

                if (n_y > 0) col_ym1 = -off_ym1 + stride2*(ym1-1)
                col_y =   -off_y   + stride1*(y-1)
                col_yp1 = -off_yp1 + stride*(yp1-1)

                do s_x = i,m_max-s_yp1

                   in_shell = nshell(s_x)
 
                   if (s_x .ge. la .and. s_x .le. la+lb) then
                      tgt_x = .true.
                   else
                      tgt_x = .false.
                   endif
 
                   before_s_x = ncart(s_x-1)
                   before_s_xp1 = ncart(s_x)
                   before_s_xm1 = ncart(s_x-2)
     
                   !this routine and the line below effectively do what is commented out below
                   call et_xsys_micro (x_dir, y_dir, col_ym1, col_y, col_yp1, n_y, s_x, in_shell, before_s_xm1, before_s_x, &
                                       before_s_xp1, delta, alp_ab_cd, two_eta, src1, src2, et_buf1, et_buf2, et_buf3)
                   ind_yp1 = before_s_x+in_shell + col_yp1

!                   do x_ind = 1,in_shell ![xs|
!
!                      can_x = x_ind+before_s_x
!       
!                      !indices of the auxiliary integrals involving [x]
!                      if (n_y > 0) ind_ym1 = can_x + col_ym1
!                      ind_y   = can_x + col_y  
!                      ind_yp1 = can_x + col_yp1
! 
!                      !indices of the auxiliary integrals involving [x+1] and [x-1]
!                      n_x = 0
!                      if (x_dir) then
!                         ind_xp1 = before_s_xp1 + x_ind + col_y !before_s_xp1 + can_shell(s_x+1,cart_l(can_x)+1,cart_n(can_x)) + col_y
!                         if (cart_l(can_x) > 0) then !X-1
!                            n_x = cart_l(can_x)
!                            ind_xm1 = before_s_xm1 + x_ind + col_y !before_s_xm1 + can_shell(s_x-1,cart_l(can_x)-1,cart_n(can_x)) + col_y
!                         endif
!                      elseif (y_dir) then
!                         ind_xp1 = before_s_xp1 + x_ind + s_x - cart_l(can_x) + 1 + col_y !before_s_xp1 + can_shell(s_x+1,cart_l(can_x),cart_n(can_x)) + col_y
!                         if (cart_m(can_x) > 0) then !Y-1
!                            n_x = cart_m(can_x)
!                            ind_xm1 = before_s_xm1 + x_ind - s_x + cart_l(can_x) + col_y !before_s_xm1 + can_shell(s_x-1,cart_l(can_x),cart_n(can_x)) + col_y
!                         endif
!                      else !z
!                         ind_xp1 = before_s_xp1 + x_ind + 2 + s_x - cart_l(can_x) + col_y !before_s_xp1 + can_shell(s_x+1,cart_l(can_x),cart_n(can_x)+1) + col_y
!                         if (cart_n(can_x) > 0) then !Z-1
!                            n_x = cart_n(can_x)
!                            ind_xm1 = before_s_xm1 + x_ind - 1 -s_x + cart_l(can_x) + col_y !before_s_xm1 + can_shell(s_x-1,cart_l(can_x),cart_n(can_x)-1) + col_y
!                         endif
!                      endif
!       
!                      !ET step:
!                      ![x s|y+1 s] = delta_i*[xs|ys] - ab/cd*[(x+1_i)s|ys] + N_i(x)/(2*cd)*[(x-1_i)s|ys] + N_i(y)/(2*cd)*[xs|(y-1_i)s]
!                      if (src1 .eq. 1 .and. src2 .eq. 2) then
!                         et_buf3(ind_yp1) = delta*et_buf1(ind_y) - alp_ab_cd*et_buf1(ind_xp1)
!                         if (n_x > 0) et_buf3(ind_yp1) = et_buf3(ind_yp1) + n_x*two_eta*et_buf1(ind_xm1)
!                         if (n_y > 0) et_buf3(ind_yp1) = et_buf3(ind_yp1) + n_y*two_eta*et_buf2(ind_ym1)
!                      elseif (src1 .eq. 2 .and. src2 .eq. 3) then
!                         et_buf1(ind_yp1) = delta*et_buf2(ind_y) - alp_ab_cd*et_buf2(ind_xp1)
!                         if (n_x > 0) et_buf1(ind_yp1) = et_buf1(ind_yp1) + n_x*two_eta*et_buf2(ind_xm1)
!                         if (n_y > 0) et_buf1(ind_yp1) = et_buf1(ind_yp1) + n_y*two_eta*et_buf3(ind_ym1)
!                      else !(src1 .eq. 3 .and. src2 .eq. 1)
!                         et_buf2(ind_yp1) = delta*et_buf3(ind_y) - alp_ab_cd*et_buf3(ind_xp1)
!                         if (n_x > 0) et_buf2(ind_yp1) = et_buf2(ind_yp1) + n_x*two_eta*et_buf3(ind_xm1)
!                         if (n_y > 0) et_buf2(ind_yp1) = et_buf2(ind_yp1) + n_y*two_eta*et_buf1(ind_ym1)
!                      endif
!          
!                   enddo !x_ind

                   if (tgt_y .and. tgt_x) then !the generated set of auxiliaries [xs|ys] is one of the target auxiliaries, so we copy it into et_tgt
                      if (src1 .eq. 1) et_tgt(ind+1:ind+in_shell) = et_buf3(ind_yp1-in_shell+1:ind_yp1)
                      if (src1 .eq. 2) et_tgt(ind+1:ind+in_shell) = et_buf1(ind_yp1-in_shell+1:ind_yp1)
                      if (src1 .eq. 3) et_tgt(ind+1:ind+in_shell) = et_buf2(ind_yp1-in_shell+1:ind_yp1)
   
                      ind = ind + in_shell
                   endif

                enddo !s_x
                
             enddo !yp1

             !prepare for the next line of auxiliaries: shift the source and the target buffers
             if (src1 .eq. 1 .and. src2 .eq. 2) then
                src1 = 3; src2 = 1
             elseif (src1 .eq. 3 .and. src2 .eq. 1) then
                src1 = 2; src2 = 3
             else !(src1 .eq. 2 .and. src2 .eq. 3)
                src1 = 1; src2 = 2
             endif

          enddo !s_yp1

   end subroutine et_xsys

   !ET step for all x from one shell and a fixed (y+1):
   ![x s|y+1 s] = delta_i*[xs|ys] - ab/cd*[(x+1_i)s|ys] + N_i(x)/(2*cd)*[(x-1_i)s|ys] + N_i(y)/(2*cd)*[xs|(y-1_i)s]
   subroutine et_xsys_micro (x_dir, y_dir, col_ym1, col_y, col_yp1, n_y, s_x, in_shell, before_s_xm1, before_s_x, &
                             before_s_xp1, delta, alp_ab_cd, two_eta, src1, src2, et_buf1, et_buf2, et_buf3)
      implicit none
      logical, intent(in) :: x_dir,y_dir
      integer, intent(in) :: col_y,col_yp1,col_ym1,n_y,s_x,in_shell,before_s_xm1,before_s_x,before_s_xp1
      integer(kind=shortint), intent(in) :: src1, src2
      real(kind=cfp), intent(in) :: delta, alp_ab_cd, two_eta
      real(kind=cfp), intent(inout) :: et_buf1(*), et_buf2(*), et_buf3(*)

         if (x_dir) then
            if (src1 .eq. 1 .and. src2 .eq. 2) then
               call et_xsys_micro_X_dir(et_buf1,et_buf2,et_buf3,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            elseif (src1 .eq. 2 .and. src2 .eq. 3) then
               call et_xsys_micro_X_dir(et_buf2,et_buf3,et_buf1,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            else !(src1 .eq. 3 .and. src2 .eq. 1)
               call et_xsys_micro_X_dir(et_buf3,et_buf1,et_buf2,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            endif
         elseif (y_dir) then
            if (src1 .eq. 1 .and. src2 .eq. 2) then
               call et_xsys_micro_Y_dir(et_buf1,et_buf2,et_buf3,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            elseif (src1 .eq. 2 .and. src2 .eq. 3) then
               call et_xsys_micro_Y_dir(et_buf2,et_buf3,et_buf1,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            else !(src1 .eq. 3 .and. src2 .eq. 1)
               call et_xsys_micro_Y_dir(et_buf3,et_buf1,et_buf2,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            endif
         else !z_dir
            if (src1 .eq. 1 .and. src2 .eq. 2) then
               call et_xsys_micro_Z_dir(et_buf1,et_buf2,et_buf3,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            elseif (src1 .eq. 2 .and. src2 .eq. 3) then
               call et_xsys_micro_Z_dir(et_buf2,et_buf3,et_buf1,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            else !(src1 .eq. 3 .and. src2 .eq. 1)
               call et_xsys_micro_Z_dir(et_buf3,et_buf1,et_buf2,n_y,s_x,in_shell,before_s_xm1,&
                before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
            endif
         endif

   end subroutine et_xsys_micro

   !Vectorization of the X_dir pays off since the X direction is the one most often used.
   subroutine et_xsys_micro_X_dir(src1,src2,tgt,n_y,s_x,in_shell,before_s_xm1,&
                                    before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
      implicit none
      real(kind=cfp), intent(in) :: src1(*), src2(*)
      real(kind=cfp), intent(in) :: delta,alp_ab_cd,two_eta
      real(kind=cfp), intent(out) :: tgt(*)
      integer, intent(in) :: n_y,s_x,in_shell,before_s_xm1,before_s_x,before_s_xp1,col_ym1,col_y,col_yp1

      integer :: x_ind, n_x, ind_ym1, ind_y, ind_yp1, ind_xm1, ind_xp1, i, can, ind1_b, ind2_b, ind3_b, ind4_b, ind5_b
      real(kind=cfp) :: tx, ty


         if (n_y > 0) then
            ind_ym1 = col_ym1+before_s_x
            ind_y   = col_y  +before_s_x
            ind_yp1 = col_yp1+before_s_x
            ind_xp1 = before_s_xp1 + col_y
            ind_xm1 = before_s_xm1 + col_y

            !n_x > 0
            can = 0
            do i=0,s_x-1
               n_x = s_x-i
               tx = n_x*two_eta
               ty = n_y*two_eta
               x_ind = can+1
               !calculate the base indices for the source and target buffers below
               ind1_b = ind_yp1+x_ind
               ind2_b = ind_y+x_ind
               ind3_b = ind_xp1+x_ind
               ind4_b = ind_xm1+x_ind
               ind5_b = ind_ym1+x_ind
               if (i .eq. 0) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
               elseif (i .eq. 1) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
               elseif (i .eq. 2) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
               elseif (i .eq. 3) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
                  tgt(ind1_b+3) = delta*src1(ind2_b+3) - alp_ab_cd*src1(ind3_b+3) + tx*src1(ind4_b+3) + ty*src2(ind5_b+3)
               elseif (i .eq. 4) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
                  tgt(ind1_b+3) = delta*src1(ind2_b+3) - alp_ab_cd*src1(ind3_b+3) + tx*src1(ind4_b+3) + ty*src2(ind5_b+3)
                  tgt(ind1_b+4) = delta*src1(ind2_b+4) - alp_ab_cd*src1(ind3_b+4) + tx*src1(ind4_b+4) + ty*src2(ind5_b+4)
               elseif (i .eq. 5) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
                  tgt(ind1_b+3) = delta*src1(ind2_b+3) - alp_ab_cd*src1(ind3_b+3) + tx*src1(ind4_b+3) + ty*src2(ind5_b+3)
                  tgt(ind1_b+4) = delta*src1(ind2_b+4) - alp_ab_cd*src1(ind3_b+4) + tx*src1(ind4_b+4) + ty*src2(ind5_b+4)
                  tgt(ind1_b+5) = delta*src1(ind2_b+5) - alp_ab_cd*src1(ind3_b+5) + tx*src1(ind4_b+5) + ty*src2(ind5_b+5)
               elseif (i .eq. 6) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
                  tgt(ind1_b+3) = delta*src1(ind2_b+3) - alp_ab_cd*src1(ind3_b+3) + tx*src1(ind4_b+3) + ty*src2(ind5_b+3)
                  tgt(ind1_b+4) = delta*src1(ind2_b+4) - alp_ab_cd*src1(ind3_b+4) + tx*src1(ind4_b+4) + ty*src2(ind5_b+4)
                  tgt(ind1_b+5) = delta*src1(ind2_b+5) - alp_ab_cd*src1(ind3_b+5) + tx*src1(ind4_b+5) + ty*src2(ind5_b+5)
                  tgt(ind1_b+6) = delta*src1(ind2_b+6) - alp_ab_cd*src1(ind3_b+6) + tx*src1(ind4_b+6) + ty*src2(ind5_b+6)
               elseif (i .eq. 7) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
                  tgt(ind1_b+3) = delta*src1(ind2_b+3) - alp_ab_cd*src1(ind3_b+3) + tx*src1(ind4_b+3) + ty*src2(ind5_b+3)
                  tgt(ind1_b+4) = delta*src1(ind2_b+4) - alp_ab_cd*src1(ind3_b+4) + tx*src1(ind4_b+4) + ty*src2(ind5_b+4)
                  tgt(ind1_b+5) = delta*src1(ind2_b+5) - alp_ab_cd*src1(ind3_b+5) + tx*src1(ind4_b+5) + ty*src2(ind5_b+5)
                  tgt(ind1_b+6) = delta*src1(ind2_b+6) - alp_ab_cd*src1(ind3_b+6) + tx*src1(ind4_b+6) + ty*src2(ind5_b+6)
                  tgt(ind1_b+7) = delta*src1(ind2_b+7) - alp_ab_cd*src1(ind3_b+7) + tx*src1(ind4_b+7) + ty*src2(ind5_b+7)
               elseif (i .eq. 8) then
                  tgt(ind1_b) = delta*src1(ind2_b) - alp_ab_cd*src1(ind3_b) + tx*src1(ind4_b) + ty*src2(ind5_b)
                  tgt(ind1_b+1) = delta*src1(ind2_b+1) - alp_ab_cd*src1(ind3_b+1) + tx*src1(ind4_b+1) + ty*src2(ind5_b+1)
                  tgt(ind1_b+2) = delta*src1(ind2_b+2) - alp_ab_cd*src1(ind3_b+2) + tx*src1(ind4_b+2) + ty*src2(ind5_b+2)
                  tgt(ind1_b+3) = delta*src1(ind2_b+3) - alp_ab_cd*src1(ind3_b+3) + tx*src1(ind4_b+3) + ty*src2(ind5_b+3)
                  tgt(ind1_b+4) = delta*src1(ind2_b+4) - alp_ab_cd*src1(ind3_b+4) + tx*src1(ind4_b+4) + ty*src2(ind5_b+4)
                  tgt(ind1_b+5) = delta*src1(ind2_b+5) - alp_ab_cd*src1(ind3_b+5) + tx*src1(ind4_b+5) + ty*src2(ind5_b+5)
                  tgt(ind1_b+6) = delta*src1(ind2_b+6) - alp_ab_cd*src1(ind3_b+6) + tx*src1(ind4_b+6) + ty*src2(ind5_b+6)
                  tgt(ind1_b+7) = delta*src1(ind2_b+7) - alp_ab_cd*src1(ind3_b+7) + tx*src1(ind4_b+7) + ty*src2(ind5_b+7)
                  tgt(ind1_b+8) = delta*src1(ind2_b+8) - alp_ab_cd*src1(ind3_b+8) + tx*src1(ind4_b+8) + ty*src2(ind5_b+8)
               elseif (i > 8) then
                  forall (x_ind=can+1:can+i+1)
                     tgt(ind_yp1+x_ind) = delta*src1(ind_y+x_ind) - alp_ab_cd*src1(ind_xp1+x_ind) &
                                                + tx*src1(ind_xm1+x_ind) + ty*src2(ind_ym1+x_ind)
                  endforall
               endif
               can = can + i + 1
            enddo

            !n_x = 0: i = s_x
            ty = n_y*two_eta
            forall (x_ind=in_shell-s_x:in_shell)
               tgt(ind_yp1+x_ind) = delta*src1(ind_y+x_ind) - alp_ab_cd*src1(ind_xp1+x_ind) + ty*src2(ind_ym1+x_ind)
            endforall
         else
            ind_y   = col_y  +before_s_x
            ind_yp1 = col_yp1+before_s_x
            ind_xp1 = before_s_xp1 + col_y
            ind_xm1 = before_s_xm1 + col_y

            !n_x > 0
            can = 0
            do i=0,s_x-1
               n_x = s_x-i
               tx = n_x*two_eta
               forall (x_ind=can+1:can+i+1)
                  tgt(ind_yp1+x_ind) = delta*src1(ind_y+x_ind) - alp_ab_cd*src1(ind_xp1+x_ind) + tx*src1(ind_xm1+x_ind)
               endforall
               can = can + i + 1
            enddo

            !n_x = 0: i = s_x
            forall (x_ind=in_shell-s_x:in_shell)
               tgt(ind_yp1+x_ind) = delta*src1(ind_y+x_ind) - alp_ab_cd*src1(ind_xp1+x_ind)
            endforall
         endif

   end subroutine et_xsys_micro_X_dir

   !Note that this part could be vectorized similarly to the X_dir case. However, the loop count would be smaller and therefore it doesn't seem like it would be worth doing.
   !Also note that the contribution of the terms ind_y* can be vectorized over x_ind.
   subroutine et_xsys_micro_Y_dir(src1,src2,tgt,n_y,s_x,in_shell,before_s_xm1,before_s_x,&
                                    before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
      implicit none
      real(kind=cfp), intent(in) :: src1(*), src2(*)
      real(kind=cfp), intent(in) :: delta,alp_ab_cd,two_eta
      real(kind=cfp), intent(out) :: tgt(*) 
      integer, intent(in) :: n_y,s_x,in_shell,before_s_xm1,before_s_x,before_s_xp1,col_ym1,col_y,col_yp1

      integer :: x_ind, can_x, n_x, ind_ym1, ind_y, ind_yp1, ind_xm1, ind_xp1

         do x_ind = 1,in_shell ![xs|

            can_x = x_ind+before_s_x
       
            !indices of the auxiliary integrals involving [x]
            if (n_y > 0) ind_ym1 = can_x + col_ym1
            ind_y   = can_x + col_y  
            ind_yp1 = can_x + col_yp1
 
            !indices of the auxiliary integrals involving [x+1] and [x-1]
            n_x = 0
            ind_xp1 = before_s_xp1 + x_ind + s_x - cart_l(can_x) + 1 + col_y !before_s_xp1 + can_shell(s_x+1,cart_l(can_x),cart_n(can_x)) + col_y
            if (cart_m(can_x) > 0) then !Y-1
               n_x = cart_m(can_x)
               ind_xm1 = before_s_xm1 + x_ind - s_x + cart_l(can_x) + col_y !before_s_xm1 + can_shell(s_x-1,cart_l(can_x),cart_n(can_x)) + col_y
            endif
       
            !ET step:
            ![x s|y+1 s] = delta_i*[xs|ys] - ab/cd*[(x+1_i)s|ys] + N_i(x)/(2*cd)*[(x-1_i)s|ys] + N_i(y)/(2*cd)*[xs|(y-1_i)s]
            tgt(ind_yp1) = delta*src1(ind_y) - alp_ab_cd*src1(ind_xp1)
            if (n_x > 0) tgt(ind_yp1) = tgt(ind_yp1) + n_x*two_eta*src1(ind_xm1)
            if (n_y > 0) tgt(ind_yp1) = tgt(ind_yp1) + n_y*two_eta*src2(ind_ym1)
         
         enddo !x_ind

   end subroutine et_xsys_micro_Y_dir

   !Note that this part could be vectorized similarly to the X_dir case. However, the loop count would be even smaller than for Y_dir and therefore it doesn't seem like it would be worth doing.
   !Also note that the contribution of the terms ind_y* can be vectorized over x_ind.
   subroutine et_xsys_micro_Z_dir(src1,src2,tgt,n_y,s_x,in_shell,before_s_xm1,&
                            before_s_x,before_s_xp1,col_ym1,col_y,col_yp1,delta,alp_ab_cd,two_eta)
      implicit none
      real(kind=cfp), intent(in) :: src1(*), src2(*)
      real(kind=cfp), intent(in) :: delta,alp_ab_cd,two_eta
      real(kind=cfp), intent(out) :: tgt(*) 
      integer, intent(in) :: n_y,s_x,in_shell,before_s_xm1,before_s_x,before_s_xp1,col_ym1,col_y,col_yp1

      integer :: x_ind, can_x, n_x, ind_ym1, ind_y, ind_yp1, ind_xm1, ind_xp1

         do x_ind = 1,in_shell ![xs|

            can_x = x_ind+before_s_x
       
            !indices of the auxiliary integrals involving [x]
            if (n_y > 0) ind_ym1 = can_x + col_ym1
            ind_y   = can_x + col_y  
            ind_yp1 = can_x + col_yp1
 
            !indices of the auxiliary integrals involving [x+1] and [x-1]
            n_x = 0
            ind_xp1 = before_s_xp1 + x_ind + 2 + s_x - cart_l(can_x) + col_y !before_s_xp1 + can_shell(s_x+1,cart_l(can_x),cart_n(can_x)+1) + col_y
            if (cart_n(can_x) > 0) then !Z-1
               n_x = cart_n(can_x)
               ind_xm1 = before_s_xm1 + x_ind - 1 -s_x + cart_l(can_x) + col_y !before_s_xm1 + can_shell(s_x-1,cart_l(can_x),cart_n(can_x)-1) + col_y
            endif
       
            !ET step:
            ![x s|y+1 s] = delta_i*[xs|ys] - ab/cd*[(x+1_i)s|ys] + N_i(x)/(2*cd)*[(x-1_i)s|ys] + N_i(y)/(2*cd)*[xs|(y-1_i)s]
            tgt(ind_yp1) = delta*src1(ind_y) - alp_ab_cd*src1(ind_xp1)
            if (n_x > 0) tgt(ind_yp1) = tgt(ind_yp1) + n_x*two_eta*src1(ind_xm1)
            if (n_y > 0) tgt(ind_yp1) = tgt(ind_yp1) + n_y*two_eta*src2(ind_ym1)
         
         enddo !x_ind

   end subroutine et_xsys_micro_Z_dir

   !Transforms the 'ab' part of the integrals batch [ab|cd] into the spherical basis. la,lb are the angular momenta in the 'a' and 'b' shells to be transformed and nc,nd are the number of terms 
   !in the 'c' and 'd' shells, i.e. nc, nd are not angular momenta. This allows the routine to be used also on the [ab|cd] batch where the 'cd' part is spherical. Note that this routine can be
   !(and is used) to transform the one-electron (ab) integrals. This can be done setting nc=nd=1.
   !SVN revision 490 has a version of the cartesian->spherical transform which is more straightforward and exploits dgemm more but is slower. Especially for high L the present method is much more efficient.
   subroutine sh_ab(cart_ints,ab_sph_ints,la,lb,nc,nd)
      implicit none
      integer, intent(in) :: la, lb, nc, nd
      real(kind=cfp), intent(in) :: cart_ints(*)
      real(kind=cfp), intent(out) :: ab_sph_ints(*)

      real(kind=cfp), parameter :: alpha = 1.0_cfp
      real(kind=cfp), parameter :: beta = 0.0_cfp
      integer :: crt_shell_a, crt_shell_b, crt_shell_ab, err, ind_start, ind_end, sph_shell_a, &
                                    sph_shell_b, ind_ambm, sph_shell_ab, m, ind, i, j, k, ind_t, ind_am
      integer :: n_bcd, n_cd, space, b, ma, mb, ind_start_a, ind_end_a, ind_start_b, nz_ab, ind_a, &
                                ind_b, no_batches, no_rest, bcol1,bcol2,bcol3,bcol4,bcol5,bcol6,bcol7,bcol8, col
      real(kind=cfp) :: tr1,tr2,tr3,tr4,tr5,tr6,tr7,tr8

         if (la > max_l .or. lb > max_l) stop "sh_ab: a or b in (ab|cd) has L greater than max_l"
         if (la < 0 .or. lb < 0) stop "sh_ab: la and/or lb are < 0"
         if (nc .le. 0 .or. nd .le. 0) stop "sh_ab: nc .le. 0 .or. nd .le. 0"

         crt_shell_a = nshell(la)
         crt_shell_b = nshell(lb)
         crt_shell_ab = crt_shell_a*crt_shell_b

         sph_shell_a = 2*la+1
         sph_shell_b = 2*lb+1
         sph_shell_ab = sph_shell_a*sph_shell_b

         if (la .le. 1 .and. lb .le. 1) then
            !no need to transform a and b: transfer cart_ints to ab_sph_ints
            ab_sph_ints(1:sph_shell_ab*nc*nd) = cart_ints(1:sph_shell_ab*nc*nd)
         endif

         !Transform only a
         if (la > 1 .and. lb .le. 1) then
            n_bcd = crt_shell_b*nc*nd

            ind_t = 0
            ind_am = 0
            ind = la*la
            ind_start = sph_map_start(ind+1)
            !Loop over all M values in the 'a' shell of the [ab|cd] batch.
            do m=-la,la
               ind_am = ind_am+1 !m+l+1: index within the 'a' shell of the spherical function with M=m.
               ind = ind + 1   !nz(ind) = how many non-zero elements there are in the current cartesian-> spherical (l,m) transformation matrix.

               col = (ind_am-1)*n_bcd !column in the half_sph array corresponding to the current M value, i.e the order of the indices in half_sph is: bcda

               ind_end = ind_start + nz(ind)

               !The loop over 'i' runs over all bcd indices in the [ab|cd] batch. Therefore in each iteration the (cartesian) indices bcd are fixed (e.g. b_i c_j d_k) and we perform the 
               !cartesian -> spherical transformation on 'a' given by the dot product: (a_1 b_i c_j d_k, a_2 b_i c_j d_k, ..., a_{crt_shell_a} b_i c_j d_k).(c_1,c_2,...,c_{crt_shell_a}) = a_m b_i c_j d_k.
               !However, in order to speed up the computation we use only the non-zero cartesian->spherical 'c' coefficients (and the corresponding 'a' terms). These non-zero 'c' coefficients are stored 
               !sequentially in the c_nz array. ind_t+nz_can(.) are indices of the corresponding 'a' terms in the set (a_1 b_i c_j d_k, a_2 b_i c_j d_k, ..., a_{crt_shell_a} b_i c_j d_k). The computational
               !savings introduced by this method (as opposed to the simple dgemm multiplication of the whole cartesian 'a' shell by the cartesian -> spherical transformation) are significant especially
               !for high L: e.g. for L=10 the cartesian->spherical transformation for the whole shell contains only 17.7% non-zero elements. The value nz(ind) is the number of the non-zero 
               !cartesian->spherical terms for each l,m combination (which specified by index 'ind').
               !Unrolled loops are used to aid the computation of the dot product. The resulting values a_m b_i c_j d_k are saved in the half_sph 1d array which emulates a 2d array (n_bcd,2*la+1). Finally, 
               !this array is transposed to obtain [ab|cd].
               if (nz(ind) .eq. 1) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start))
                  enddo
               elseif (nz(ind) .eq. 2) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                     + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1))
                  enddo
               elseif (nz(ind) .eq. 3) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                        + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1)) &
                                        + c_nz(ind_start+2)*cart_ints(ind_t+nz_can(ind_start+2))
                  enddo
               elseif (nz(ind) .eq. 4) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                     + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1)) &
                                     + c_nz(ind_start+2)*cart_ints(ind_t+nz_can(ind_start+2)) &
                                     + c_nz(ind_start+3)*cart_ints(ind_t+nz_can(ind_start+3))
                  enddo
               elseif (nz(ind) .eq. 5) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                     + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1)) &
                                     + c_nz(ind_start+2)*cart_ints(ind_t+nz_can(ind_start+2)) &
                                     + c_nz(ind_start+3)*cart_ints(ind_t+nz_can(ind_start+3)) &
                                     + c_nz(ind_start+4)*cart_ints(ind_t+nz_can(ind_start+4))
                  enddo
               elseif (nz(ind) .eq. 6) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                     + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1)) &
                                     + c_nz(ind_start+2)*cart_ints(ind_t+nz_can(ind_start+2)) &
                                     + c_nz(ind_start+3)*cart_ints(ind_t+nz_can(ind_start+3)) &
                                     + c_nz(ind_start+4)*cart_ints(ind_t+nz_can(ind_start+4)) &
                                     + c_nz(ind_start+5)*cart_ints(ind_t+nz_can(ind_start+5))
                  enddo
               elseif (nz(ind) .eq. 7) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                     + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1)) &
                                     + c_nz(ind_start+2)*cart_ints(ind_t+nz_can(ind_start+2)) &
                                     + c_nz(ind_start+3)*cart_ints(ind_t+nz_can(ind_start+3)) &
                                     + c_nz(ind_start+4)*cart_ints(ind_t+nz_can(ind_start+4)) &
                                     + c_nz(ind_start+5)*cart_ints(ind_t+nz_can(ind_start+5)) &
                                     + c_nz(ind_start+6)*cart_ints(ind_t+nz_can(ind_start+6))
                  enddo
               elseif (nz(ind) .eq. 8) then
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = c_nz(ind_start)*cart_ints(ind_t+nz_can(ind_start)) &
                                     + c_nz(ind_start+1)*cart_ints(ind_t+nz_can(ind_start+1)) &
                                     + c_nz(ind_start+2)*cart_ints(ind_t+nz_can(ind_start+2)) &
                                     + c_nz(ind_start+3)*cart_ints(ind_t+nz_can(ind_start+3)) &
                                     + c_nz(ind_start+4)*cart_ints(ind_t+nz_can(ind_start+4)) &
                                     + c_nz(ind_start+5)*cart_ints(ind_t+nz_can(ind_start+5)) &
                                     + c_nz(ind_start+6)*cart_ints(ind_t+nz_can(ind_start+6)) &
                                     + c_nz(ind_start+7)*cart_ints(ind_t+nz_can(ind_start+7))
                  enddo
               elseif (nz(ind) > 8) then
                  no_batches = nz(ind)/8   !how many batches of terms (a_1 bcd,...,a_8 bcd) we have to deal with
                  no_rest = mod(nz(ind),8) !the rest of the integrals (a_9 bcd,...,a_{9+no_rest-1} bcd)
                  do i=1,n_bcd
                     ind_t = crt_shell_a*(i-1)
                     half_sph(i+col) = 0.0_cfp
                     k = ind_start
                     do j=1,no_batches     !transform batches of 8 at once
                        half_sph(i+col) = half_sph(i+col) + c_nz(k)*cart_ints(ind_t+nz_can(k)) &
                                        + c_nz(k+1)*cart_ints(ind_t+nz_can(k+1)) &
                                        + c_nz(k+2)*cart_ints(ind_t+nz_can(k+2)) &
                                        + c_nz(k+3)*cart_ints(ind_t+nz_can(k+3)) &
                                        + c_nz(k+4)*cart_ints(ind_t+nz_can(k+4)) &
                                        + c_nz(k+5)*cart_ints(ind_t+nz_can(k+5)) &
                                        + c_nz(k+6)*cart_ints(ind_t+nz_can(k+6)) &
                                        + c_nz(k+7)*cart_ints(ind_t+nz_can(k+7))
                        k = k + 8
                     enddo
                     !transform the remaining terms
                     do j=1,no_rest
                        half_sph(i+col) = half_sph(i+col) + c_nz(k+j-1)*cart_ints(ind_t+nz_can(k+j-1))
                     enddo
                  enddo !i
               endif
               ind_start = ind_start + nz(ind)

            enddo !m

            !transpose bcda (half_sph) to abcd (ab_sph_ints)
            call abcd_to_cdab(half_sph,ab_sph_ints,n_bcd,1,sph_shell_a,1)

         endif

         !Transform both: a,b
         !The principle of the transformation is the same as in the a-only case. However, the sparseness of the cartesian->spherical transformation means that from the full set of pairs of 'ab' cartesian
         !indices only a subset of those will contribute to each cartesian -> spherical transformation: [ab|cd] -> a_m b_m|cd]. Therefore the dot products within each loop (do i=1,n_cd) over the cartesian
         !indices cd are computed only from the non-zero cartesian -> spherical transformation coefficients and the corresponding cartesian integrals. Before the dot product is computed we gather the
         !non-zero coefficients for the 'ab' pairs and the corresponding indices of the 'a' and 'b' cartesians. The resulting values a_m b_m c_j d_k are saved in the half_sph 1d array which emulates a 
         !2d array (n_cd,(2*la+1)*(2lb+1)). Finally, this array is transposed to obtain [ab|cd].

         if (la > 1 .and. lb > 1) then
            n_cd = nc*nd

            ind_a = la*la
            ind_b = lb*lb
            !todo do this in the main allocation routine?
            ma = maxval(nz(ind_a+1:ind_a+sph_shell_a))
            mb = maxval(nz(ind_b+1:ind_b+sph_shell_b))
            space = ma*mb
            err = check_real_array_size(c_ab,space)
            if (err .ne. 0) stop "sh_ab: transf allocation 3 failed"
            if (allocated(nz_can_ab)) then 
               if (size(nz_can_ab,2) < space) then 
                  deallocate(nz_can_ab)
                  allocate(nz_can_ab(1:2,space),stat=err)
                  if (err .ne. 0) stop "sh_ab: transf allocation 4 failed"
               endif
            else
               allocate(nz_can_ab(1:2,space),stat=err)
               if (err .ne. 0) stop "sh_ab: transf allocation 4 failed"
            endif

            ind_b = lb*lb
            ind_t = 0
            ind_start_b = sph_map_start(ind_b+1)
            ind_ambm = 0
            do mb=-lb,lb
               ind_b = ind_b + 1
               ind_a = la*la
               ind_start_a = sph_map_start(ind_a+1)
               do ma=-la,la
                  ind_a = ind_a + 1
                  ind_ambm = ind_ambm + 1

                  nz_ab = nz(ind_a)*nz(ind_b)
                  ind_end_a = ind_start_a + nz(ind_a) - 1

                  !transform a_ma,b_mb at once
                  i = 0
                  do b=1,nz(ind_b)
                     c_ab(i+1:i+nz(ind_a)) = c_nz(ind_start_b+b-1)*c_nz(ind_start_a:ind_end_a)
                     nz_can_ab(1,i+1:i+nz(ind_a)) = nz_can(ind_start_a:ind_end_a)
                     nz_can_ab(2,i+1:i+nz(ind_a)) = nz_can(ind_start_b+b-1)
                     i = i + nz(ind_a)
                  enddo

                  col = (ind_ambm-1)*n_cd !column in the half_sph array corresponding to the current Ma,Mb values.
                  if (nz_ab .eq. 1) then
                     tr1 = c_ab(1)
                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1)
                     enddo
                  elseif (nz_ab .eq. 2) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2)
                     enddo
                  elseif (nz_ab .eq. 3) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)
                     tr3 = c_ab(3)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     bcol3 = crt_shell_a*(nz_can_ab(2,3)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2) + &
                                                              &tr3*cart_ints(ind_t+nz_can_ab(1,3)+bcol3)
                     enddo
                  elseif (nz_ab .eq. 4) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)
                     tr3 = c_ab(3)
                     tr4 = c_ab(4)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     bcol3 = crt_shell_a*(nz_can_ab(2,3)-1)
                     bcol4 = crt_shell_a*(nz_can_ab(2,4)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2) + &
                                                              &tr3*cart_ints(ind_t+nz_can_ab(1,3)+bcol3) + &
                                                              &tr4*cart_ints(ind_t+nz_can_ab(1,4)+bcol4)
                     enddo
                  elseif (nz_ab .eq. 5) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)
                     tr3 = c_ab(3)
                     tr4 = c_ab(4)
                     tr5 = c_ab(5)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     bcol3 = crt_shell_a*(nz_can_ab(2,3)-1)
                     bcol4 = crt_shell_a*(nz_can_ab(2,4)-1)
                     bcol5 = crt_shell_a*(nz_can_ab(2,5)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2) + &
                                                              &tr3*cart_ints(ind_t+nz_can_ab(1,3)+bcol3) + &
                                                              &tr4*cart_ints(ind_t+nz_can_ab(1,4)+bcol4) + &
                                                              &tr5*cart_ints(ind_t+nz_can_ab(1,5)+bcol5)
                     enddo
                  elseif (nz_ab .eq. 6) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)
                     tr3 = c_ab(3)
                     tr4 = c_ab(4)
                     tr5 = c_ab(5)
                     tr6 = c_ab(6)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     bcol3 = crt_shell_a*(nz_can_ab(2,3)-1)
                     bcol4 = crt_shell_a*(nz_can_ab(2,4)-1)
                     bcol5 = crt_shell_a*(nz_can_ab(2,5)-1)
                     bcol6 = crt_shell_a*(nz_can_ab(2,6)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2) + &
                                                              &tr3*cart_ints(ind_t+nz_can_ab(1,3)+bcol3) + &
                                                              &tr4*cart_ints(ind_t+nz_can_ab(1,4)+bcol4) + &
                                                              &tr5*cart_ints(ind_t+nz_can_ab(1,5)+bcol5) + &
                                                              &tr6*cart_ints(ind_t+nz_can_ab(1,6)+bcol6)
                     enddo
                  elseif (nz_ab .eq. 7) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)
                     tr3 = c_ab(3)
                     tr4 = c_ab(4)
                     tr5 = c_ab(5)
                     tr6 = c_ab(6)
                     tr7 = c_ab(7)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     bcol3 = crt_shell_a*(nz_can_ab(2,3)-1)
                     bcol4 = crt_shell_a*(nz_can_ab(2,4)-1)
                     bcol5 = crt_shell_a*(nz_can_ab(2,5)-1)
                     bcol6 = crt_shell_a*(nz_can_ab(2,6)-1)
                     bcol7 = crt_shell_a*(nz_can_ab(2,7)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2) + &
                                                              &tr3*cart_ints(ind_t+nz_can_ab(1,3)+bcol3) + &
                                                              &tr4*cart_ints(ind_t+nz_can_ab(1,4)+bcol4) + &
                                                              &tr5*cart_ints(ind_t+nz_can_ab(1,5)+bcol5) + &
                                                              &tr6*cart_ints(ind_t+nz_can_ab(1,6)+bcol6) + &
                                                              &tr7*cart_ints(ind_t+nz_can_ab(1,7)+bcol7)
                     enddo
                  elseif (nz_ab .eq. 8) then
                     tr1 = c_ab(1)
                     tr2 = c_ab(2)
                     tr3 = c_ab(3)
                     tr4 = c_ab(4)
                     tr5 = c_ab(5)
                     tr6 = c_ab(6)
                     tr7 = c_ab(7)
                     tr8 = c_ab(8)

                     bcol1 = crt_shell_a*(nz_can_ab(2,1)-1)
                     bcol2 = crt_shell_a*(nz_can_ab(2,2)-1)
                     bcol3 = crt_shell_a*(nz_can_ab(2,3)-1)
                     bcol4 = crt_shell_a*(nz_can_ab(2,4)-1)
                     bcol5 = crt_shell_a*(nz_can_ab(2,5)-1)
                     bcol6 = crt_shell_a*(nz_can_ab(2,6)-1)
                     bcol7 = crt_shell_a*(nz_can_ab(2,7)-1)
                     bcol8 = crt_shell_a*(nz_can_ab(2,8)-1)
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = tr1*cart_ints(ind_t+nz_can_ab(1,1)+bcol1) + &
                                                              &tr2*cart_ints(ind_t+nz_can_ab(1,2)+bcol2) + &
                                                              &tr3*cart_ints(ind_t+nz_can_ab(1,3)+bcol3) + &
                                                              &tr4*cart_ints(ind_t+nz_can_ab(1,4)+bcol4) + &
                                                              &tr5*cart_ints(ind_t+nz_can_ab(1,5)+bcol5) + &
                                                              &tr6*cart_ints(ind_t+nz_can_ab(1,6)+bcol6) + &
                                                              &tr7*cart_ints(ind_t+nz_can_ab(1,7)+bcol7) + &
                                                              &tr8*cart_ints(ind_t+nz_can_ab(1,8)+bcol8)
                     enddo
                  elseif (nz_ab > 8) then
                     no_batches = nz_ab/8   !how many batches of 8 we have to deal with
                     no_rest = mod(nz_ab,8) !the rest of the integrals
                     do i=1,n_cd
                        ind_t = crt_shell_a*crt_shell_b*(i-1)
                        half_sph(i+col) = 0.0_cfp
                        k = 0
                        do j=1,no_batches     !transform batches of 8 at once
                           tr1 = c_ab(k+1)
                           tr2 = c_ab(k+2)
                           tr3 = c_ab(k+3)
                           tr4 = c_ab(k+4)
                           tr5 = c_ab(k+5)
                           tr6 = c_ab(k+6)
                           tr7 = c_ab(k+7)
                           tr8 = c_ab(k+8)

                           bcol1 = crt_shell_a*(nz_can_ab(2,k+1)-1)
                           bcol2 = crt_shell_a*(nz_can_ab(2,k+2)-1)
                           bcol3 = crt_shell_a*(nz_can_ab(2,k+3)-1)
                           bcol4 = crt_shell_a*(nz_can_ab(2,k+4)-1)
                           bcol5 = crt_shell_a*(nz_can_ab(2,k+5)-1)
                           bcol6 = crt_shell_a*(nz_can_ab(2,k+6)-1)
                           bcol7 = crt_shell_a*(nz_can_ab(2,k+7)-1)
                           bcol8 = crt_shell_a*(nz_can_ab(2,k+8)-1)

                           half_sph(i+col) = half_sph(i+col) + &
                                                                 &tr1*cart_ints(ind_t+nz_can_ab(1,k+1)+bcol1) + &
                                                                 &tr2*cart_ints(ind_t+nz_can_ab(1,k+2)+bcol2) + &
                                                                 &tr3*cart_ints(ind_t+nz_can_ab(1,k+3)+bcol3) + &
                                                                 &tr4*cart_ints(ind_t+nz_can_ab(1,k+4)+bcol4) + &
                                                                 &tr5*cart_ints(ind_t+nz_can_ab(1,k+5)+bcol5) + &
                                                                 &tr6*cart_ints(ind_t+nz_can_ab(1,k+6)+bcol6) + &
                                                                 &tr7*cart_ints(ind_t+nz_can_ab(1,k+7)+bcol7) + &
                                                                 &tr8*cart_ints(ind_t+nz_can_ab(1,k+8)+bcol8)
                           k = k + 8
                        enddo
                        !contributions from the remaining terms
                        do j=1,no_rest
                           half_sph(i+col) = half_sph(i+col) + c_ab(k+j)*cart_ints(ind_t+nz_can_ab(1,k+j) &
                                                                +crt_shell_a*(nz_can_ab(2,k+j)-1))
                        enddo
                     enddo !i
                  endif
                  ind_start_a = ind_start_a + nz(ind_a)
               enddo !ma
               ind_start_b = ind_start_b + nz(ind_b)
            enddo !mb

            !transpose cdab (half_sph) to abcd (ab_sph_ints)
            call abcd_to_cdab(half_sph,ab_sph_ints,nc,nd,sph_shell_a,sph_shell_b)

         endif

   end subroutine sh_ab

   !Transforms all cd terms from the batch of integrals [ab|cd] into the spherical basis. na,nb are the number of terms in the 'a' and 'b' shells and lc,ld are the angular momenta in the shells 'c' and 'd'.
   !The calculation proceeds by first transposing the batch to [cd|ab]. The 'cd' part is then transformed using the sh_ab routine.
   !SVN revision 490 has a version of the cartesian->spherical transform which is more straightforward and exploits dgemm more but is generally slower. Especially for high L the present method is much more
   !efficient.
   subroutine sh_cd(ab_sph_ints,sph_ints,na,nb,lc,ld)
      implicit none
      integer, intent(in) :: na, nb, lc, ld
      real(kind=cfp), intent(inout) :: ab_sph_ints(*)
      real(kind=cfp), intent(out) :: sph_ints(*)

      real(kind=cfp), parameter :: alpha = 1.0_cfp
      real(kind=cfp), parameter :: beta = 0.0_cfp
      integer :: crt_shell_c, crt_shell_d, n_ab, sph_shell_c, sph_shell_d, sph_shell_cd

         if (lc > max_l .or. ld > max_l) stop "c or d in (ab|cd) has L greater than max_l"
         if (lc < 0 .or. ld < 0) stop "sh_cd: lc and/or ld are < 0"
         if (na .le. 0 .or. nb .le. 0) stop "sh_cd: na .le. 0 .or. nb .le. 0"

         n_ab = na*nb

         sph_shell_c = 2*lc+1
         sph_shell_d = 2*ld+1
         sph_shell_cd = sph_shell_c*sph_shell_d

         crt_shell_c = nshell(lc)
         crt_shell_d = nshell(ld)

         if (lc .le. 1 .and. ld .le. 1) then
            !no need to transform c and d
            sph_ints(1:n_ab*sph_shell_cd) = ab_sph_ints(1:n_ab*sph_shell_cd)
            return
          endif

         !transpose the batch and produce [cd|ab]: ab_sph_ints -> transf
         call abcd_to_cdab(ab_sph_ints,transf,na,nb,crt_shell_c,crt_shell_d)

         !apply the spherical transform on the [cd| part: transf -> ab_sph_ints
         call sh_ab(transf,ab_sph_ints,lc,ld,na,nb)

         !transpose the batch again and produce [ab|cd]: ab_sph_ints -> sph_ints
         call abcd_to_cdab(ab_sph_ints,sph_ints,sph_shell_c,sph_shell_d,na,nb)

   end subroutine sh_cd

   !number of cartesian functions within shell of angular momentum l
   elemental function nshell(l)
      implicit none
      integer, intent(in) :: l
      integer :: nshell

         nshell = (l+1)*(l+2)/2

   end function nshell

   !number of cartesians in all shells with angular momentum .le. l
   elemental function ncart(l)
      implicit none
      integer, intent(in) :: l
      integer :: ncart

         ncart = (l+1)*(l+2)*(l+3)/6

   end function ncart

   !calculates the full canonical index given the shell angular momentum ixyz and the x,z exponents ix, iz
   elemental function can(ixyz,ix,iz)
      implicit none
      integer :: can
      integer, intent(in) :: ixyz,ix,iz

         can = ixyz*(ixyz+1)*(ixyz+2)/6 + (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1

   end function can

   !calculates the canonical index within the shell of angular momentum ixyz
   elemental function can_shell(ixyz,ix,iz)
      implicit none
      integer :: can_shell
      integer, intent(in) :: ixyz,ix,iz

         can_shell = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1

   end function can_shell

   elemental function dist2(x1,y1,z1,x2,y2,z2)
      implicit none
      real(kind=cfp) :: dist2
      real(kind=cfp), intent(in) :: x1,y1,z1,x2,y2,z2

         dist2 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)

   end function dist2

   elemental function product_center_1D(alphaa,xa,alphab,xb)
      implicit none
      real(kind=cfp) :: product_center_1D
      real(kind=cfp), intent(in) :: alphaa,xa,alphab,xb
      
         product_center_1D = (alphaa*xa+alphab*xb)/(alphaa+alphab)

   end function product_center_1D
!
! ROUTINES FOR OVERLAP, KE AND MULTIPOLE MOMENT INTEGRALS
!
   !Calculates overlap and kinetic energy integrals for a pair of shells of spherical GTOs.
   subroutine sph_olap_kei(lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs,ind_a, &
                           lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs,ind_b, &
                           olap_column,kei_column,integrals,int_index)
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, ind_a, ind_b, olap_column, kei_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      real(kind=cfp), intent(in) :: acnorm,anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(out) :: integrals(:,:)
      integer, intent(out) :: int_index(:,:)

         !reorder the input data so that the function a has never angular momentum smaller than b: this is important for the sh_ab routine.
         if (la .ge. lb) then !la,lb
            call sph_olap_kei_shell (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                                     lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, olap_column,kei_column,integrals)
            call index_1el(la,lb,ind_a,ind_b,1,int_index)
         else !lb,la
            call sph_olap_kei_shell (lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, &
                                     lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, olap_column,kei_column,integrals)
            call index_1el(lb,la,ind_b,ind_a,1,int_index)
         endif

   end subroutine sph_olap_kei

   !Calculates the multipole moment integrals for a pair of shells of spherical GTOs and the multipole
   !centered at (xc,yc,zc) for all multipoles l .le. lc. For l .eq. 0 the multipole integrals reduce to
   !the overlap integrals between the two shells of GTOs.
   !\warning The ordering of the output integrals (i.e. la .ge. lb) is important for the calculation of the property integral
   !         tails so this needs to be consistent. See the comments in prop_cms_tail and the source for eri_tail_shell.
   subroutine sph_mult_mom (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs,ind_a, lc,xc,yc,zc, &
                            lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs,ind_b, property_column,sph_mult,int_index)
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, lc, ind_a, ind_b, property_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc
      real(kind=cfp), intent(in) :: acnorm,anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(out) :: sph_mult(:,:)
      integer, intent(out) :: int_index(:,:)

         !reorder the input data so that the function a has never angular momentum smaller than b: this is important for the sh_ab routine.
         if (la .ge. lb) then !la,lb,lc
            call sph_mult_mom_shell (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, lc,xc,yc,zc, &
                                     lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, property_column,sph_mult)
            call index_1el(la,lb,ind_a,ind_b,1,int_index)
         else !lb,la,lc
            call sph_mult_mom_shell (lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, lc,xc,yc,zc, &
                                     lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, property_column,sph_mult)
            call index_1el(lb,la,ind_b,ind_a,1,int_index)
         endif

   end subroutine sph_mult_mom

   !This routine assumes that the input data have been reorganized so that la .ge. lb.
   subroutine sph_mult_mom_shell (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, lc,xc,yc,zc, &
                                  lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, property_column,sph_mult_mom)
      use phys_const, only: pi, fourpi
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, lc, property_column
      real(kind=cfp), intent(in) :: xa,ya,za, xc,yc,zc, xb,yb,zb
      real(kind=cfp), intent(in) :: acnorm,anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(out) :: sph_mult_mom(:,:)

      integer :: i, j, m, column, i_cart_start, i_cart_end, l, err, sph_shell_a, sph_shell_b, sph_shell_c, &
                 nshell_a, nshell_b, ncart_lc, sph_shell_ab, sph_shell_mom
      real(kind=cfp) :: Rab(1:3), Rpa(1:3), Rac(1:3), alp_ab, K_ab, mu, r_sq, prod

         !calculate the canonical indices
         call calc_can(max(la,lb,lc))

         !space and data for the cartesian to spherical transform
         call allocate_space_sph_transf(la,lb,lc,0)

         sph_shell_a = 2*la+1
         sph_shell_b = 2*lb+1
         sph_shell_ab = sph_shell_a*sph_shell_b

         nshell_a = nshell(la)
         nshell_b = nshell(lb)
         ncart_lc = ncart(lc)
         sph_shell_c = 2*lc+1

         !auxiliary arrays for the cartesian intermediate integrals
         i = nshell_a*nshell_b*ncart_lc
         err = check_real_array_size(cart_mult_mom,i); if (err .ne. 0) stop "cart_mult_mom allocation failed."
         err = check_real_array_size(contr_cart_mult_mom,i); if (err .ne. 0) stop "contr_cart_mult_mom allocation failed."
         j = sph_shell_ab*sph_shell_c
         err = check_real_array_size(shell_mom,j); if (err .ne. 0) stop "shell_mom allocation failed."
         
         !we must make half_sph larger than what allocate_space_sph_transf does since we need it accommodate
         !a larger number of auxiliaries in the call to sh_ab below than what is required in eri.
         i = sph_shell_a*sph_shell_b*ncart_lc
         err = check_real_array_size(half_sph,i); if (err .ne. 0) stop "sph_mult_mom_shell: half_sph allocation failed."

         Rab(1) = xa-xb;
         Rab(2) = ya-yb;
         Rab(3) = za-zb;
         r_sq = dot_product(Rab,Rab)

         Rac(1) = xa-xc;
         Rac(2) = ya-yc;
         Rac(3) = za-zc;

         contr_cart_mult_mom = 0.0_cfp

         do i=1,lena
            do j=1,lenb

               prod = anorms(i)*acoefs(i)*bnorms(j)*bcoefs(j)*(acnorm*bcnorm)

               alp_ab = aexps(i)+bexps(j)
               mu = aexps(i)*bexps(j)/alp_ab
               K_ab = (pi/alp_ab)*sqrt(pi/alp_ab)*exp(-mu*r_sq)

               Rpa(1) = product_center_1D(aexps(i),xa,bexps(j),xb) - xa
               Rpa(2) = product_center_1D(aexps(i),ya,bexps(j),yb) - ya
               Rpa(3) = product_center_1D(aexps(i),za,bexps(j),zb) - za

               call prim_cart_mult_mom(la,lc,lb,Rab,Rpa,Rac,K_ab,alp_ab,cart_mult_mom)

               contr_cart_mult_mom = contr_cart_mult_mom + prod*cart_mult_mom

            enddo
         enddo

         if (la > 0 .or. lb > 0 .or. lc > 0) then
            !Spherical transform on the GTO A and the GTO B
            call sh_ab(contr_cart_mult_mom,cart_mult_mom,la,lb,ncart_lc,1) !the result is in cart_mult_mom
      
            if (lc > 0) then
               i_cart_start = 0
               i_cart_end = 0
               !perform the spherical transform on each shell of multipole moments and reorder the values in each shell if needed.
               column = property_column-1
               do l=0,lc
                  sph_shell_c = 2*l+1
  
                  i_cart_end = i_cart_end + sph_shell_ab*nshell(l) 

                  sph_shell_mom = 2*l+1
                  call sh_cd (cart_mult_mom(i_cart_start+1:i_cart_end), &
                              shell_mom(1:sph_shell_ab*sph_shell_mom),&
                              sph_shell_a,sph_shell_b,l,0)
   
                  !In case some of the shells are p-type we need to permute the p-type integrals from the Cartesian order (x,y,z) into the spherical order M=-1,0,1: (y,z,x).
                  call reorder_p_shells(shell_mom(1:sph_shell_ab*sph_shell_mom),la,lb,l,0)

                  !Transfer the transformed properties from the linear array to the columns of the final array:
                  i = 0
                  do m=-l,l
                     column = column+1
                     sph_mult_mom(1:sph_shell_ab,column) = shell_mom(i+1:i+sph_shell_ab)
                     i = i + sph_shell_ab
                  enddo
   
                  i_cart_start = i_cart_end
               enddo
            else !lc .eq. 0
               !In case some of the shells are p-type we need to permute the p-type integrals from the Cartesian order (x,y,z) into the spherical order M=-1,0,1: (y,z,x).
               call reorder_p_shells(cart_mult_mom,la,lb,0,0)

               sph_mult_mom(1:sph_shell_ab,property_column) = cart_mult_mom(1:sph_shell_ab)
            endif

         else !la .eq. 0, lb .eq. 0, lc .eq. 0
            sph_mult_mom(1,property_column) = contr_cart_mult_mom(1)
         endif

   end subroutine sph_mult_mom_shell

   !> Calculates the cartesian multipole moment integrals for a pair of shells of primitive GTOs and a given shell L of the multipole moment.
   subroutine prim_cart_mult_mom(la,lc,lb,Rab,Rpa,Rac,K_ab,alp_ab,cart_mom)
      use phys_const, only: fourpi
      implicit none
      integer, intent(in) :: la, lc, lb
      real(kind=cfp), intent(in) :: Rab(1:3), Rpa(1:3), Rac(1:3), K_ab, alp_ab
      real(kind=cfp), intent(out) :: cart_mom(:)

      real(kind=cfp) :: S(0:max(la+lb+lc,1),0:lb,0:lc,1:3) !the first dimension must be at least 1 since arrays with dimensions 0,0,0,: are not allowed
      integer :: i, j, k, l, ind, nshell_a, nshell_b, nshell_c, ncart_a, ncart_b, can_i, can_j, can_k

         cart_mom = 0.0_cfp

         call mult_mom_recurrence(S,Rab,Rpa,Rac,la,lc,lb,alp_ab)

         nshell_a = nshell(la)
         nshell_b = nshell(lb)

         ncart_a = ncart(la-1)
         ncart_b = ncart(lb-1)

         !Loop over the canonical indices of a,b,c and assemble the 3D integrals.
         !We calculate at once all multipole integrals for l=0,...,lc.
         ind = 0
         can_k = 0
         do l=0,lc
            nshell_c = nshell(l)
            do k=1,nshell_c
               can_k = can_k + 1 != k + ncart_c
               can_j = ncart_b
               do j=1,nshell_b
                  can_j = can_j + 1 != j + ncart_b
                  can_i = ncart_a
                  do i=1,nshell_a
                     ind = ind + 1 != i + (j-1)*nshell_a + (can_k-1)*nshell_a*nshell_b
                     can_i = can_i + 1 != i + ncart_a
                     !assemble the multipole integral from the integrals over the x,y,z coordinates, e.g. can_l(can_i) is the exponent for the x coordinate in the GTO A with the canonical index can_i.
                     cart_mom(ind) = K_ab*S(cart_l(can_i),cart_l(can_j),cart_l(can_k),1) &
                                         *S(cart_m(can_i),cart_m(can_j),cart_m(can_k),2) &
                                         *S(cart_n(can_i),cart_n(can_j),cart_n(can_k),3)
                  enddo !i
               enddo !j
            enddo !k
         enddo !l

   end subroutine prim_cart_mult_mom

   !> This routine implements the Obara-Saika recurrent relations for the GTO auxiliary overlap integrals needed for calculation of multipole moment integrals for a pair of cartesian GTOs.
   !> See Helgaker - Sections 9.3.2 for the equations.
   !> \param[out] S On output S contains the auxiliary overlap integrals. Prior call to this routine the array S(0:d1,0:d2,0:d3,1:3) must be allocated with d1 .ge. la+lb+lc, d2 .ge. lb, d3 .ge. lc.
   !> \param[in] Rab Real vector \f$\mathbf{R}_{ab}=\mathbf{R}_{a}-\mathbf{R}_{b}\f$, where \f$\mathbf{R}_{a}\f$ is the center of the GTO a and \f$\mathbf{R}_{b}\f$ is the center of the GTO b.
   !> \param[in] Rpa Real vector \f$\mathbf{R}_{pa}=\mathbf{R}_{p}-\mathbf{R}_{a}\f$, where \f$\mathbf{R}_{p}\f$ is the center of the product GTO, while \f$\mathbf{R}_{a}\f$ is the center of the GTO a.
   !> \param[in] Rac Real vector \f$\mathbf{R}_{ac}=\mathbf{R}_{a}-\mathbf{R}_{c}\f$, where \f$\mathbf{R}_{a}\f$ is the center of the GTO a and \f$\mathbf{R}_{c}\f$ is the center of the multipole moment.
   !> \param[in] la Angular momentum on the GTO a.
   !> \param[in] lc Multipole moment L.
   !> \param[in] lb Angular momentum on the GTO b.
   !> \param[in] alp_ab Sum of the exponents on the two GTOs.
   subroutine mult_mom_recurrence(S,Rab,Rpa,Rac,la,lc,lb,alp_ab)
      implicit none
      integer, intent(in) :: la, lc, lb
      real(kind=cfp), intent(out) :: S(0:max(la+lb+lc,1),0:lb,0:lc,1:3)
      real(kind=cfp), intent(in) :: Rab(1:3), Rpa(1:3), Rac(1:3), alp_ab
            
      integer :: i, j, k, l, sl
      real(kind=cfp) :: fac
            
            sl = la+lb+lc
            S = 0.0_cfp
            !starting values for the recursion
            do i=1,3
               S(0,0,0,i) = 1.0_cfp
            enddo
            !generate S(x,0)^0
            fac = 1.0_cfp/(2.0_cfp*alp_ab)
            do i=1,3
               S(1,0,0,i) = Rpa(i)*S(0,0,0,i)
               do j=2,sl
                  S(j,0,0,i) = Rpa(i)*S(j-1,0,0,i) + (j-1)*fac*S(j-2,0,0,i)
               enddo !j
            enddo !i
            !generate S(x,0)^k, k=1,...,lc
            do i=1,3
               do k=1,lc
                  do j=0,sl-k
                     S(j,0,k,i) = S(j+1,0,k-1,i) + Rac(i)*S(j,0,k-1,i)
                  enddo
               enddo
            enddo
            !Horizontal recurrence: S(x,y)^z
            !todo tiling could be used here to improve memory locality but is it worth doing?
            do i=1,3 
               do k=0,lc
                  do l=1,lb
                     do j=0,la+lb-l
                        S(j,l,k,i) = Rab(i)*S(j,l-1,k,i)+S(j+1,l-1,k,i)
                     enddo
                  enddo
               enddo
            enddo

   end subroutine mult_mom_recurrence

   !This routine assumes that the input data have been reorganized so that la .ge. lb.
   subroutine sph_olap_kei_shell (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                                  lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, olap_column,kei_column,integrals)
      use phys_const, only: pi
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, olap_column,kei_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      real(kind=cfp), intent(in) :: acnorm,anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(out) :: integrals(:,:)

      integer :: i, j, err, n
      real(kind=cfp) :: Rab(1:3), Rpa(1:3), alp_ab, K_ab, mu, r_sq, prod

         !calculate the canonical indices
         call calc_can(max(la,lb))

         !space and data for the cartesian to spherical transform
         call allocate_space_sph_transf(la,lb,0,0)

         Rab(1) = xa-xb;
         Rab(2) = ya-yb;
         Rab(3) = za-zb;
         r_sq = dot_product(Rab,Rab)

         n = nshell(la)*nshell(lb)
         err = check_real_array_size(cart_olap_buf,n); if (err .ne. 0) stop "cart_olap_buf allocation failed."
         err = check_real_array_size(cart_kei_buf,n); if (err .ne. 0) stop "cart_kei_buf allocation failed."
         err = check_real_array_size(contr_cart_olap,n); if (err .ne. 0) stop "contr_cart_olap allocation failed."
         err = check_real_array_size(contr_cart_kei,n); if (err .ne. 0) stop "contr_cart_kei allocation failed."

         contr_cart_olap = 0.0_cfp
         contr_cart_kei = 0.0_cfp

         do i=1,lena
            do j=1,lenb

               prod = anorms(i)*acoefs(i)*bnorms(j)*bcoefs(j)

               alp_ab = aexps(i)+bexps(j)
               mu = aexps(i)*bexps(j)/alp_ab
               K_ab = (pi/alp_ab)**1.5_cfp*exp(-mu*r_sq)

               Rpa(1) = product_center_1D(aexps(i),xa,bexps(j),xb) - xa
               Rpa(2) = product_center_1D(aexps(i),ya,bexps(j),yb) - ya
               Rpa(3) = product_center_1D(aexps(i),za,bexps(j),zb) - za

               call prim_cart_olap_kei(la,lb,Rab,Rpa,K_ab,aexps(i),alp_ab,cart_olap_buf,cart_kei_buf)

               contr_cart_olap = contr_cart_olap + prod*cart_olap_buf
               contr_cart_kei = contr_cart_kei + prod*cart_kei_buf

            enddo
         enddo

         n = (2*la+1)*(2*lb+1)

         if (la > 0 .or. lb > 0) then
            if (olap_column > 0) call sh_ab(contr_cart_olap,integrals(1:n,olap_column),la,lb,1,1)
            if (kei_column > 0) call sh_ab(contr_cart_kei ,integrals(1:n,kei_column) ,la,lb,1,1)

            !In case some of the shells are p-type we need to permute the p-type integrals from the Cartesian order (x,y,z) into the spherical order M=-1,0,1: y,z,x.
            if (olap_column > 0) call reorder_p_shells(integrals(1:n,olap_column),la,lb,0,0)
            if (kei_column > 0) call reorder_p_shells(integrals(1:n,kei_column) ,la,lb,0,0)
         else
            if (olap_column > 0) integrals(1:n,olap_column) = contr_cart_olap(1:n)
            if (kei_column > 0) integrals(1:n,kei_column) = contr_cart_kei(1:n)
         endif

         !multiply by the norms of the contracted GTOs
         if (olap_column > 0) integrals(1:n,olap_column) = acnorm*bcnorm*integrals(1:n,olap_column)
         if (kei_column > 0) integrals(1:n,kei_column) = acnorm*bcnorm*integrals(1:n,kei_column)
 
   end subroutine sph_olap_kei_shell

   subroutine prim_cart_olap_kei(la,lb,Rab,Rpa,K_ab,a,alp_ab,cart_olap,cart_kei)
      implicit none
      integer, intent(in) :: la, lb
      real(kind=cfp), intent(in) :: Rab(1:3), Rpa(1:3), K_ab, a, alp_ab
      real(kind=cfp), intent(out) :: cart_olap(:), cart_kei(:)

      real(kind=cfp) :: S0(0:la+lb+2,0:lb,1:3), D2(0:la,0:lb,1:3)
      integer :: i, j, ind, nshell_a, nshell_b, can_i, can_j

         cart_olap = 0.0_cfp
         cart_kei = 0.0_cfp

         call olap_ke_recurrence(S0,Rab,Rpa,la,lb,alp_ab)
         call S0_to_D2(S0,D2,a,la,lb)

         nshell_a = nshell(la)
         nshell_b = nshell(lb)

         !loop over the canonical indices of a,b
         do j=1,nshell_b
            can_j = j + ncart(lb-1)
            do i=1,nshell_a
               ind = i + (j-1)*nshell_a
               can_i = i + ncart(la-1)
               !assemble the overlap integral
               cart_olap(ind) = K_ab * S0(cart_l(can_i),cart_l(can_j),1) &
                                     * S0(cart_m(can_i),cart_m(can_j),2) &
                                     * S0(cart_n(can_i),cart_n(can_j),3)
               !assemble the kinetic energy integral
               cart_kei(ind) = -0.5_cfp * K_ab * &
               (  D2(cart_l(can_i),cart_l(can_j),1) * S0(cart_m(can_i),cart_m(can_j),2) * S0(cart_n(can_i),cart_n(can_j),3) &
                + S0(cart_l(can_i),cart_l(can_j),1) * D2(cart_m(can_i),cart_m(can_j),2) * S0(cart_n(can_i),cart_n(can_j),3) &
                + S0(cart_l(can_i),cart_l(can_j),1) * S0(cart_m(can_i),cart_m(can_j),2) * D2(cart_n(can_i),cart_n(can_j),3) )
            enddo
         enddo

   end subroutine prim_cart_olap_kei

   !> Calculates overlap integral between a pair of primitive cartesian functions. Note that this routine is used
   !> for conversion of orbital coefficients from one basis to another so it does not need to be
   !> very efficient or sophisticated. Therefore this routine is different to the sph_olap_kei which calculates the integrals
   !> over the whole pair of shells of functions and returns the KE integral as well.
   subroutine cart_olap(xa,ya,za,ix,iy,iz,aexp, xb,yb,zb,jx,jy,jz,bexp, olap)
      implicit none
      integer, intent(in) :: ix,iy,iz, jx,jy,jz
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      real(kind=cfp), intent(in) :: aexp, bexp
      real(kind=cfp), intent(out) :: olap

      integer :: la, lb

         !order the shells according to their angular momentum (la .ge. lb)
         la = ix+iy+iz
         lb = jx+jy+jz
         if (ix+iy+iz .ge. jx+jy+jz) then
            call cart_olap_pair(xa,ya,za,la,ix,iy,iz,aexp, xb,yb,zb,lb,jx,jy,jz,bexp, olap)
         else
            call cart_olap_pair(xb,yb,zb,lb,jx,jy,jz,bexp, xa,ya,za,la,ix,iy,iz,aexp, olap)
         endif

   end subroutine cart_olap

   !> Assuming la .ge. lb (la = ix+iy+iz, lb = jx+jy+jz) this routine calculates the overlap integral between a pair of primitive Gaussian functions.
   subroutine cart_olap_pair(xa,ya,za,la,ix,iy,iz,aexp, xb,yb,zb,lb,jx,jy,jz,bexp, olap)
      use phys_const, only: pi
      implicit none
      integer, intent(in) :: la,ix,iy,iz, lb,jx,jy,jz
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      real(kind=cfp), intent(in) :: aexp, bexp
      real(kind=cfp), intent(out) :: olap

      real(kind=cfp) :: Rab(1:3), Rpa(1:3), alp_ab, K_ab, mu, r_sq
      real(kind=cfp) :: S0(0:la+lb+2,0:lb,1:3)

         Rab(1) = xa-xb;
         Rab(2) = ya-yb;
         Rab(3) = za-zb;
         r_sq = dot_product(Rab,Rab)

         alp_ab = aexp+bexp
         mu = aexp*bexp/alp_ab
         K_ab = (pi/alp_ab)**1.5_cfp*exp(-mu*r_sq)

         Rpa(1) = product_center_1D(aexp,xa,bexp,xb) - xa
         Rpa(2) = product_center_1D(aexp,ya,bexp,yb) - ya
         Rpa(3) = product_center_1D(aexp,za,bexp,zb) - za

         !generate the 1D auxiliaries
         call olap_ke_recurrence(S0,Rab,Rpa,la,lb,alp_ab)

         !assemble the overlap integral
         olap = K_ab*S0(ix,jx,1)*S0(iy,jy,2)*S0(iz,jz,3)

   end subroutine cart_olap_pair
   
   !> This routine implements the Obara-Saika recurrent relations for the GTO auxiliary overlap integrals needed for calculation of overlaps and kinetic energy integrals for a pair of cartesian GTOs.
   !> See Helgaker - Sections 9.3.1, 9.3.4 for the equations.
   !> \param[out] S0 On output S0 contains the auxiliary overlap integrals. Prior call to this routine the array S0(0:d1,0:d2,1:3) must be allocated with d1 .ge. 2*max(la,lb)+2, d2 .ge. max(la,lb).
   !> \param[in] Rab Real vector \f$\mathbf{R}_{ab}=\mathbf{R}_{a}-\mathbf{R}_{b}\f$, where \f$\mathbf{R}_{a}\f$ is the center of the GTO a and \f$\mathbf{R}_{b}\f$ is the center of the GTO b.
   !> \param[in] Rpa Real vector \f$\mathbf{R}_{pa}=\mathbf{R}_{p}-\mathbf{R}_{a}\f$, where \f$\mathbf{R}_{p}\f$ is the center of the product GTO, while \f$\mathbf{R}_{a}\f$ is the center of the GTO a.
   !> \param[in] la Angular momentum on the GTO a.
   !> \param[in] lb Angular momentum on the GTO b.
   !> \param[in] alp_ab Sum of the exponents on the two GTOs.
   subroutine olap_ke_recurrence(S0,Rab,Rpa,la,lb,alp_ab)
      implicit none
      integer, intent(in) :: la, lb
      real(kind=cfp), intent(out) :: S0(0:la+lb+2,0:lb,1:3)
      real(kind=cfp), intent(in) :: Rab(1:3), Rpa(1:3), alp_ab
            
      integer :: i, j, k, sl
            
            sl = la+lb+2
            S0 = 0.0_cfp
            !starting values for the recursion
            do i=1,3
               S0(0,0,i) = 1.0_cfp
            enddo
            !generate S(x,0)
            do i=1,3
               S0(1,0,i) = Rpa(i)*S0(0,0,i)
               do j=2,sl
                  S0(j,0,i) = Rpa(i)*S0(j-1,0,i) + (j-1)/(2.0_cfp*alp_ab)*S0(j-2,0,i)
               enddo !j
            enddo !i
            !Horizontal recurrence: S(x,y)
            do i=1,3 
               do j=1,lb
                  do k=0,sl-j !the largest k needed for KEI is la+2, but we need k_max=(la+lb)+1 due to the last term in the RR below that uses S0(k+1,...)
                     S0(k,j,i) = Rab(i)*S0(k,j-1,i)+S0(k+1,j-1,i)
                  enddo !k
               enddo !j
            enddo !i

   end subroutine olap_ke_recurrence

   !Generates the D2 terms (needed to calculate the kinetic energy integral) from the S0 terms generated using the routine olap_ke_recurrence.
   subroutine S0_to_D2(S0,D2,a,la,lb)
      implicit none
      integer, intent(in) :: la,lb
      real(kind=cfp), intent(in) :: S0(0:la+lb+2,0:lb,1:3), a
      real(kind=cfp), intent(out) :: D2(0:la,0:lb,1:3)

      integer :: i,j,k

         do k=1,3
            do j=0,lb
               do i=0,la
                  D2(i,j,k) = 4.0_cfp*a*a*S0(i+2,j,k) - 2.0_cfp*a*(2*i+1)*S0(i,j,k) + i*(i-1)*S0(max(i-2,0),j,k)
               enddo
            enddo
         enddo

   end subroutine S0_to_D2
!
! ROUTINES FOR NUCLEAR ATTRACTION/REPULSION INTEGRALS
!
   subroutine sph_nari (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs,ind_a, &
                        lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs,ind_b, xc,yc,zc, sph_nari_int,int_index)
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, ind_a, ind_b
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc
      real(kind=cfp), intent(in) :: acnorm,anorms(lena),aexps(lena),acoefs(lena)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(lenb),bexps(lenb),bcoefs(lenb)
      real(kind=cfp), intent(out) :: sph_nari_int(:)
      integer, intent(out) :: int_index(:,:)

         !reorder the input data so that the function a has never angular momentum smaller than b: this is important for the sh_ab routine.
         if (la .ge. lb) then !la,lb
            call sph_nari_shell (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                                 lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, xc,yc,zc, sph_nari_int)
            call index_1el(la,lb,ind_a,ind_b,1,int_index)
         else !lb,la
            call sph_nari_shell (lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, &
                                 lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, xc,yc,zc, sph_nari_int)
            call index_1el(lb,la,ind_b,ind_a,1,int_index)
         endif

   end subroutine sph_nari

   !The calculation of the NAR integrals follows the same algorthim as that for the 2-electron integrals of the type (ab|ss) only with some modifications to the starting values for the recursion.
   !This routine assumes that the input data have been reorganized so that la .ge. lb. Note that we use the same temporary buffers like the 2-electron integral routines.
   !todo Multiply the result by Z*e, where Z and e are the charges of the nucleus and the particle.
   subroutine sph_nari_shell (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                              lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, xc,yc,zc, sph_nari_int)
      use const, only: mmax, imax_wp,boys_f_grid_step_wp,taylor_k_wp, imax_ep,boys_f_grid_step_ep,taylor_k_ep
      implicit none
      integer, intent(in) :: lena, lenb, la, lb
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc
      real(kind=cfp), intent(in) :: acnorm,anorms(lena),aexps(lena),acoefs(lena)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(lenb),bexps(lenb),bcoefs(lenb)
      real(kind=cfp), intent(out) :: sph_nari_int(:)

      integer :: sum_l, i, terms
      integer :: space_vrr_tgt, space_et_tgt, space_sph_ints, space_vrr_buf, space_et_buf, &
                 space_hrr1_buf, space_hrr2_buf, space_hrr1_tgt, space_hrr2_tgt

         !calculate the canonical indices
         sum_l = la+lb
         call calc_can(sum_l)

         !allocate space for all target and intermediate buffers including the space for the final integrals
         call allocate_space (la,lb,0,0,space_vrr_tgt,space_et_tgt,space_sph_ints,space_vrr_buf,&
                              space_et_buf,space_hrr1_buf,space_hrr2_buf,space_hrr1_tgt,space_hrr2_tgt)

         !initialize the Boys function; if i == 1 then it has been initialized already.
         if (cfp .eq. wp) then
            i = boys%init(imax_wp,mmax,boys_f_grid_step_wp,taylor_k_wp)
         elseif (cfp .eq. ep1) then
            i = boys%init(imax_ep,mmax,boys_f_grid_step_ep,taylor_k_ep)
         else 
            stop "sph_nari_shell: unsupported numeric type"
         endif

         if (i .ne. 0 .and. i .ne. 1) then
            print *,i
            stop "sph_nari_shell: initialization of the Boys function failed."
         endif

         !Calculate all (xs) auxiliaries required to generate (ab)
         call contr_vrr_nari (lena,xa,ya,za,acnorm,anorms,aexps,acoefs, &
                              lenb,xb,yb,zb,bcnorm,bnorms,bexps,bcoefs, &
                              xc,yc,zc,la,lb, contr_et,space_et_tgt,space_vrr_tgt,space_vrr_buf)

         terms = (2*la+1)*(2*lb+1)

         if (lb > 0) then !(la lb), lb .le. la: transfer angular momentum from the first function to the other
            call hrr1(la,xa,ya,za,lb,xb,yb,zb,0,0,contr_et,hrr1_tgt,space_hrr1_buf)
            call sh_ab(hrr1_tgt,sph_nari_int(1:terms),la,lb,1,1) !spherical transform
         else !(la s): no HRR steps required
            if (la > 0) then
               call sh_ab(contr_et,sph_nari_int(1:terms),la,lb,1,1) !spherical transform
            else !(ss)
               sph_nari_int(1) = contr_et(1)
            endif
         endif

         !In case some of the shells are p-type we need to permute the p-type integrals from the Cartesian order (x,y,z) into the spherical order M=-1,0,1: y,z,x.
         call reorder_p_shells(sph_nari_int(1:terms),la,lb,0,0)

   end subroutine sph_nari_shell

   subroutine contr_vrr_nari (lena,xa,ya,za,acnorm,anorms,aexps,acoefs, &
                              lenb,xb,yb,zb,bcnorm,bnorms,bexps,bcoefs, &
                              xc,yc,zc,la,lb,contr_et_tgt,size_contr_et,size_vrr_tgt,size_vrr_buff)
      implicit none
      integer, intent(in) :: lena, lenb, la,lb, size_contr_et, size_vrr_tgt, size_vrr_buff
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc, acnorm,bcnorm
      real(kind=cfp), intent(in) :: anorms(lena),aexps(lena),acoefs(lena)
      real(kind=cfp), intent(in) :: bnorms(lenb),bexps(lenb),bcoefs(lenb)
      real(kind=cfp), intent(out) :: contr_et_tgt(:)

      real(kind=cfp) :: Fm(la+lb+1),vrr_buf1(size_vrr_buff), vrr_buf2(size_vrr_buff), vrr_buf3(size_vrr_buff), vrr_tgt(size_vrr_tgt)
      real(kind=cfp) :: et_tgt(size_contr_et)

      integer :: i,j,ind
      real(kind=cfp) :: prod

      !make sure only the required elements are zeroed out
      contr_et_tgt(1:size_contr_et) = 0.0_cfp

      !todo simplify these loops in case some shells are equal
      do i=1,lena
         do j=1,lenb

            !todo this should be split into the i,j loops respectively
            prod = acnorm*acoefs(i)*anorms(i)*bcnorm*bcoefs(j)*bnorms(j) !product of all contraction coefficients and norms

            !calculate the vrr terms for this primitive combination of shells
            call vrr_nari(xa,ya,za,aexps(i),xb,yb,zb,bexps(j),xc,yc,zc,la,lb, Fm,vrr_buf1,vrr_buf2,vrr_buf3,vrr_tgt, et_tgt)

            !accummulate data for all shells into the contr_et_tgt; vectorized loop
            forall (ind=1:size_contr_et)
               contr_et_tgt(ind) = contr_et_tgt(ind) + prod*et_tgt(ind)
            endforall

         enddo
      enddo

   end subroutine contr_vrr_nari

   !VRR and ET steps of the HGP algorithm for the NAR integral.
   subroutine vrr_nari(xa,ya,za,alphaa,xb,yb,zb,alphab,xc,yc,zc,la,lb, Fm,vrr_buf1,vrr_buf2,vrr_buf3,vrr_tgt, et_tgt)
      implicit none
      real(kind=cfp), intent(in) :: xa,ya,za,alphaa, xb,yb,zb,alphab, xc,yc,zc
      integer, intent(in) :: la,lb
      real(kind=cfp), intent(out) :: Fm(:),vrr_buf1(:),vrr_buf2(:),vrr_buf3(:),vrr_tgt(:), et_tgt(:)

      real(kind=cfp) :: px,py,pz, pax,pay,paz, wpx,wpy,wpz, e_o_ez,zeta,two_zeta,alp, rab2,rwp2, T,Kab
      integer :: s,im,m_max, no_cart_vrr

         px = product_center_1D(alphaa,xa,alphab,xb)
         py = product_center_1D(alphaa,ya,alphab,yb)
         pz = product_center_1D(alphaa,za,alphab,zb)

         pax = px-xa
         pay = py-ya
         paz = pz-za

         wpx = px-xc
         wpy = py-yc
         wpz = pz-zc
         
         zeta = alphaa+alphab

         two_zeta = 1.0_cfp/(2.0_cfp*zeta)
         e_o_ez = 1.0_cfp

         alp = alphaa*alphab/zeta

         !this can be computed only once outside of the contraction loop in the contr_vrr routine
         rab2 = dist2(xa,ya,za,xb,yb,zb)

         Kab = 2.0_cfp*pi/zeta*exp(-alp*rab2)
         rwp2 = wpx*wpx+wpy*wpy+wpz*wpz
         T = zeta*rwp2

         m_max = la+lb

         no_cart_vrr = ncart(m_max)

         !calculate Fm(T) for m=0,...,m_max using the routine based on the Taylor expansion
         call boys%eval_taylor(Fm,m_max+1,T,m_max)
         !Fm(1:m_max+1) = boys_function(T,m_max)

         forall (im=0:m_max)
            vrr_buf1(im+1) = Kab*Fm(im+1)
         endforall

         !the VRR target integral [ss|ss]^(0)
         vrr_tgt(1) = vrr_buf1(1)

         if (m_max .eq. 0) then !la+lb+lc+ld .eq. 0: nothing else to do in this subroutine; this completes the integral calculation; the final integral is in et_tgt(1)
            et_tgt(1) = vrr_buf1(1)
            return
         endif
!
!-------- We get here only if la+lb > 0
!
!VRR
!
         !generate the [ps|ss]^(m) terms from the terms [ss|ss]^(m); put the target integrals into vrr_tgt
         call vrr_psss(m_max,wpx,wpy,wpz,pax,pay,paz,vrr_buf1,vrr_buf2,vrr_tgt)

         if (la .eq. 1 .and. lb .eq. 0) then !integral calculation has been completed; the final integrals are in vrr_tgt(2:4), so we transfer them into et_tgt which is used in contr_vrr
            et_tgt(1:3) = vrr_tgt(2:4)
            return
         endif

         !generate the triangle of [xs|ss]^(m) terms containing x=d,f,...,la+lb from the auxiliaries in vrr_buf1 and vrr_buf2; put the targets into vrr_tgt
         call vrr_xsss(m_max,wpx,wpy,wpz,pax,pay,paz,two_zeta,e_o_ez,vrr_buf1,vrr_buf2,vrr_buf3,vrr_tgt)

         !the target integrals are in vrr_tgt
         s = ncart(la-1) !only the integrals in the range [xs|ss] with x = [la,la+lb] are required for the next steps of the algorithm
         et_tgt(1:(no_cart_vrr-s)) = vrr_tgt(s+1:no_cart_vrr)

   end subroutine vrr_nari

end module cgto_hgp
