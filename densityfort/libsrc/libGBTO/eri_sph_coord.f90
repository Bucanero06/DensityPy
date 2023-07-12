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
!todo change the array argument definitions in the subroutines to have the 'allocatable' attribute so unit stride is assumed.
module eri_sph_coord
 use precisn
 use gto_routines, only: boys_function_obj, index_1el, index_2el, eri_tail_shell, reorder_and_index_2el
 use coupling_obj, only: couplings_type

   implicit none

   type(boys_function_obj), private :: boys
   type(couplings_type), private :: cpl

   real(kind=cfp), allocatable, private :: G_A(:), G_B(:), F_X(:), l_ijkl_L_coefficients(:)
!   real(kind=cfp), allocatable, protected :: eri_ints(:)
   real(kind=cfp), allocatable, private :: sum_AB(:), sum_CD(:), cpl_ABCD(:), contr_AB(:), cf(:), SH_ml(:,:), SH_l(:)
   real(kind=cfp), allocatable :: eri_tail_int(:)

   !> Maxium L (for the AB and CD pairs) and contraction length that the routine allocate_space has encountered so far. It is used to determine whether memory should be reallocated or not.
   integer, private :: max_l_1 = -1, max_l_2 = -2, max_len = -1

   !> Number of (n,l) values for the Lag. polynomials and a control variable which controls when should the coefficients be recalculated.
   integer, private :: n_n = -1, n_l = -1, max_l_cf = -1

   !$OMP THREADPRIVATE (boys,cpl,G_A,G_B,F_X,l_ijkl_L_coefficients,sum_AB,sum_CD,cpl_ABCD,contr_AB, &
   !$OMP &              cf,SH_ml,SH_l,eri_tail_int,max_l_1,max_l_2,max_len,n_n,n_l,max_l_cf)

contains

   subroutine eri_sph_coord_final
      implicit none

         max_l_1 = -1; max_l_2 = -1; max_len = -1
         n_n = -1; n_l = -1; max_l_cf = -1
         if (allocated(G_A)) deallocate(G_A)
         if (allocated(G_B)) deallocate(G_B)
         if (allocated(F_X)) deallocate(F_X)
         if (allocated(l_ijkl_L_coefficients)) deallocate(l_ijkl_L_coefficients)
         if (allocated(sum_AB)) deallocate(sum_AB)
         if (allocated(sum_CD)) deallocate(sum_CD)
         if (allocated(cpl_ABCD)) deallocate(cpl_ABCD)
         if (allocated(contr_AB)) deallocate(contr_AB)
         if (allocated(cf)) deallocate(cf)
         if (allocated(SH_ml)) deallocate(SH_ml)
         if (allocated(SH_l)) deallocate(SH_l)
         if (allocated(eri_tail_int)) deallocate(eri_tail_int)

   end subroutine eri_sph_coord_final

   !Calculates overlap and kinetic energy integrals for a pair of shells of spherical GTOs.
   subroutine olap_kei_sph (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs,ind_a, &
                            lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs,ind_b, olap_column,kei_column,integrals,int_index)
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, ind_a, ind_b, olap_column, kei_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      real(kind=cfp), intent(in) :: acnorm,anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(out) :: integrals(:,:)
      integer, intent(out) :: int_index(:,:)

         !reorder the input data so that the function a has never angular momentum smaller than b: this is important for the sh_ab routine.
         if (la .ge. lb) then !la,lb
            call olap_kei_shell_sph (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                                     lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, olap_column,kei_column,integrals)
            call index_1el(la,lb,ind_a,ind_b,1,int_index)
         else !lb,la
            call olap_kei_shell_sph (lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, &
                                     lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, olap_column,kei_column,integrals)
            call index_1el(lb,la,ind_b,ind_a,1,int_index)
         endif

   end subroutine olap_kei_sph

   subroutine allocate_cf_space(la,lb)
      implicit none
      integer, intent(in) :: la,lb

      integer :: sum_l, space_cf, ind, n, l, m, err
      real(kind=cfp) :: f

         if (max(la,lb) .le. max_l_cf) then 
            return
         else
            max_l_cf = max(la,lb)
            sum_l = 2*max_l_cf
         endif

         n_n = sum_l/2+2 !for overlap integrals sum_l/2+1 would be enough but for KEI +1.
         n_l = max(sum_l-2,2*max_l_cf)+1

         space_cf = n_l*n_n*(n_n+1)/2 !space for the Lag. coefficients

         !Allocate all at once - spatial locality of the arrays in memory (page faults)
         if (allocated(cf)) deallocate(cf,SH_l,SH_ml)
         allocate(cf(space_cf),SH_ml(2,(sum_l+1)**2),SH_l(sum_l+1),stat=err)
         if (err .ne. 0) stop "allocate_cf_space: memory allocation error"

         !Precalculate the coefficients for the Lag. polynomials. cf_m(l,n)
         ind = 0
         do n=0,n_n-1
            do l=0,n_l-1
               f = 1.0_cfp
               do m=0,n
                 if (m .ne. 0) f = f*m
                 ind = ind + 1 !save in order (m,l,n)
                 cf(ind) = (-1)**m/f*gen_binom(n+l+0.5_cfp,n-m) !cf(m+1)
               enddo
            enddo !l
         enddo !n

         !Precalculate the coefficients needed to generate the Solid Harmonics
         ind = 0
         do l=0,sum_l
            SH_l(l+1) = sqrt((2.0_cfp*l+1)/(2.0_cfp*l+2))
            do m=-l,l
               ind = ind + 1
               SH_ml(1,ind) = 1.0_cfp/sqrt((l+1.0_cfp+m)*(l+1.0_cfp-m))
               SH_ml(2,ind) = sqrt((l+m)*(l-m)*1.0_cfp)*SH_ml(1,ind)
            enddo
         enddo

   end subroutine allocate_cf_space

   !This routine assumes that the input data have been reorganized so that la .ge. lb.
   subroutine olap_kei_shell_sph (lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                                  lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, olap_column,kei_column,integrals)
      use phys_const, only: twopi, fourpi
      implicit none
      integer, intent(in) :: lena, lenb, la, lb, olap_column, kei_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
      real(kind=cfp), intent(in) :: acnorm,anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bcnorm,bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(out) :: integrals(:,:)

      integer :: i, j, sum_l, l, m, ind, n, range_start, range_end, ma, mb, l_min, l_max, base
      real(kind=cfp) :: adj_norms_a(lena), adj_norms_b(lenb), r_ba(3)
      real(kind=cfp) :: solid_harmonics_ab((la+lb+1)**2), contr_olap((la+lb+1)**2), contr_kei((la+lb+1)**2), sq_dist_ba
      real(kind=cfp) :: xsi, factor, factor_olap, factor_kei, prod, xsi_inv
      logical :: rA_eq_rB

         sum_l = la+lb

         !precalculate all G factors and Gaunt coefficients
         call cpl%prec_G_cf(sum_l)

         !allocate space and precalculate the coefficients for the Lag. polynomials and the solid harmonics
         call allocate_cf_space(la,lb)

         factor = sqrt(fourpi/(2*la+1.0_cfp))
         do i=1,lena
            adj_norms_a(i) = acnorm*anorms(i)/((2*aexps(i))**(la+1.5_cfp))*factor*acoefs(i)
         enddo

         factor = sqrt(fourpi/(2*lb+1.0_cfp))
         do i=1,lenb
            adj_norms_b(i) = bcnorm*bnorms(i)/((2*bexps(i))**(lb+1.5_cfp))*factor*bcoefs(i)
         enddo

         !evaluate the AB solid harmonics
         r_ba(1:3) = (/xb-xa,yb-ya,zb-za/)
         sq_dist_ba = dot_product(r_ba,r_ba)

         if (xb .eq. xa .and. yb .eq. ya .and. zb .eq. za) then
            rA_eq_rB = .true. !in this case I need solid_harmonics only for l=0, i.e. sqrt(1/fourpi)
            l_min = 0; l_max = 0
            base = 0
            integrals(1:(2*la+1)*(2*lb+1),olap_column) = 0.0_cfp
            integrals(1:(2*la+1)*(2*lb+1),kei_column) = 0.0_cfp
            if (la .ne. lb) return !atomic selection rule
            solid_harmonics_ab(1) = sqrt(1.0_cfp/fourpi) 
         else
            rA_eq_rB = .false.
            l_min = abs(la-lb); l_max = la+lb
            base = l_min*l_min
            call cfp_solh_1d(solid_harmonics_ab,r_ba(1),r_ba(2),r_ba(3),l_max)
            !we need Slm = r**l*X(l,m) but only for l_min:l_max so we skip re-normalizing the ones for l=0,...,l_min-1
            do l=l_min,l_max
               range_start = l*l+1
               range_end = range_start-1 + 2*l+1
               solid_harmonics_ab(range_start:range_end) = solid_harmonics_ab(range_start:range_end)*sqrt((2*l+1.0_cfp)/fourpi)
            enddo
         endif

         contr_olap = 0.0_cfp
         contr_kei = 0.0_cfp
         do i=1,lena
            do j=1,lenb
               xsi = aexps(i)*bexps(j)/(aexps(i)+bexps(j))
               xsi_inv = 1.0_cfp/(4.0_cfp*xsi)
               prod = (-1)**lb*(twopi)**(1.5_cfp)*exp(-xsi*sq_dist_ba)*adj_norms_a(i)*adj_norms_b(j)
               ind = 0
               do l=l_min,l_max
                  n = (la+lb-l)/2
                  factor_olap = (-1)**n * cnla(n,l,xsi_inv)*Lag_n_hlf_k(n,l,xsi*sq_dist_ba)
                  factor_kei = 0.5_cfp*(-1)**n * cnla(n+1,l,xsi_inv)*Lag_n_hlf_k(n+1,l,xsi*sq_dist_ba)
                  do m=-l,l
                     ind = ind + 1
                     contr_olap(ind) = contr_olap(ind) + prod*factor_olap*solid_harmonics_ab(base+ind)
                     contr_kei(ind) = contr_kei(ind) + prod*factor_kei*solid_harmonics_ab(base+ind)
                  enddo !m
               enddo !l
            enddo !j
         enddo !i

         if (rA_eq_rB) then
            n = (2*la+1)
            do ma=-la,la !we need to loop only over the pairs (ma,mb) where ma=mb
               i = (ma+la+1) + n*(ma+la)
               factor = cpl%rgaunt(la,la,0,ma,ma,0)
               integrals(i,olap_column) = factor*contr_olap(1)
               integrals(i,kei_column) = factor*contr_kei(1)
            enddo !ma
         else
            i = 0
            do mb=-lb,lb
              do ma=-la,la
                 i = i + 1
   
                 ind = 0
                 integrals(i,olap_column) = 0.0_cfp
                 integrals(i,kei_column) = 0.0_cfp
                 do l=l_min,l_max
                    do m=-l,l !todo there are max. 2 allowed values of m: m1 = abs(ma+mb), m2 = abs(ma-mb), m1 .ne. m2, abs(m1) .le. l, abs(m2) .le. l
                       ind = ind + 1
                       factor = cpl%rgaunt(la,lb,l,ma,mb,m)
                       integrals(i,olap_column) = integrals(i,olap_column) + factor*contr_olap(ind)
                       integrals(i,kei_column) = integrals(i,kei_column) + factor*contr_kei(ind)
                    enddo !m
                 enddo !l
   
              enddo !ma
            enddo !mb
         endif

   end subroutine olap_kei_shell_sph

   !If indexing_method .eq. 1 then the integrals are not reordered and standard indexing is used. If indexing_method .eq. 2 then the columns a,b,c,d of sph_ints are reordered so that ap.ge.bp,cp.ge.dp,ap.ge.cp.
   subroutine eri_sph (lena,xa,ya,za,anorms,la,aexps,acoefs,ind_a, &
                       lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,ind_b, &
                       lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,ind_c, &
                       lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,ind_d, &
                       two_el_column,int_index,keep_ab_cd_order,indexing_method,&
                       do_tails_for_this_quartet,ab_is_continuum,tgt_prop,tgt_pair,rmat_radius,sph_ints)
      implicit none
      integer, intent(in) :: lena, lenb, lenc, lend, la, lb, lc, ld, ind_a, ind_b, ind_c, ind_d, two_el_column, &
                             tgt_pair, indexing_method
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

         !Order the shells according to their angular momentum: la .ge .lb, lc .ge. ld, la+lb .ge. lc+ld, perform the calculation (eri_shell) and then expand and order the basis function indices (index_2el)
         if (la+lb .ge. lc+ld) then
            if (la .ge. lb) then
               if (lc .ge. ld) then !la,lb,lc,ld
                  call eri_shell_sph (lena,xa,ya,za,anorms,la,aexps,acoefs, &
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
                  call eri_shell_sph (lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                      lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                      lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                      lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,two_el_column,sph_ints)

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,ld,lc,la,lb,lena,lenb,aexps,bexps,acoefs,bcoefs, &
                                             anorms,bnorms,rmat_radius,.false.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,la,lb,ld,lc,lend,lenc,dexps,cexps,dcoefs,ccoefs, &
                                             dnorms,cnorms,rmat_radius,.true.,eri_tail_int)
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
                  call eri_shell_sph (lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
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
                  call eri_shell_sph (lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                      lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                      lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                      lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs,two_el_column,sph_ints)

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
                  call eri_shell_sph (lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                      lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                      lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                      lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,two_el_column,sph_ints)

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
                  call eri_shell_sph (lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                      lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                      lena,xa,ya,za,anorms,la,aexps,acoefs, &
                                      lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs,two_el_column,sph_ints)

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
                  call eri_shell_sph (lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                      lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                      lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                      lena,xa,ya,za,anorms,la,aexps,acoefs,two_el_column,sph_ints)

                  if (do_tails_for_this_quartet) then !Calculate and subtract the tails
                     if (ab_is_continuum) then !AB are continuum shells and CD are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lc,ld,lb,la,lenb,lena,bexps,aexps,bcoefs,acoefs, &
                                             bnorms,anorms,rmat_radius,.true.,eri_tail_int)
                     else !CD are continuum shells and AB are target shells
                        call eri_tail_shell (tgt_prop,tgt_pair,lb,la,lc,ld,lenc,lend,cexps,dexps,ccoefs,dcoefs, &
                                             cnorms,dnorms,rmat_radius,.false.,eri_tail_int)
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
                  call eri_shell_sph (lend,xd,yd,zd,dnorms,ld,dexps,dcoefs, &
                                      lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                                      lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                                      lena,xa,ya,za,anorms,la,aexps,acoefs,two_el_column,sph_ints)

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

   end subroutine eri_sph

   subroutine allocate_space(la,lb,lc,ld, lena,lenb,lenc,lend)
      implicit none
      integer, intent(in) :: la, lb, lc, ld, lena,lenb,lenc,lend

      integer :: err, sum_l, la_p_lb, n_ij, n_kl, n_ijkl, space_eri, space_contr_AB, n_l43m43, sph_c, sph_d, lc_p_ld, &
                 space_sum_CD, cont_len, l_1, l_2, sph_shell_1, sph_shell_2, ind, n, m, l
      integer :: space_cpl_AB, space_G, n_l21m21, sph_a, sph_b, space_F_X, n_contr_pair, space_sum_AB, space_cpl_ABCD, &
                 space_l_ijkl_L_coefficients, space_cf
      real(kind=cfp) :: f

         !Determine if the memory should be (re)allocated
         l_1 = max(la,lb)
         l_2 = max(lc,ld)
         cont_len = max(lena,lenb,lenc,lend)
         if (l_1 > max_l_1 .or. l_2 > max_l_2 .or. cont_len > max_len) then
            if (allocated(G_A)) deallocate(G_A,G_B,F_X,sum_AB,contr_AB,sum_CD,cpl_ABCD,l_ijkl_L_coefficients,cf,SH_ml,SH_l)
            max_l_1 = max(l_1,max_l_1)
            max_l_2 = max(l_2,max_l_2)
            max_len = max(cont_len,max_len)
         else
            return
         endif

         !The allocation is done assuming (LL|LL) shell with each contraction having max_len primitives. The commented expressions would be used considering exactly the shell parameters given on input.
         sph_shell_1 = 2*max_l_1+1
         sph_shell_2 = 2*max_l_2+1
         sum_l = 2*max_l_1+2*max_l_2   !la+lb+lc+ld
         la_p_lb = 2*max_l_1 !la+lb
         sph_a = sph_shell_1 !2*la+1
         sph_b = sph_shell_1 !2*lb+1
         lc_p_ld = 2*max_l_2 !lc+ld
         sph_c = sph_shell_2 !2*lc+1
         sph_d = sph_shell_2 !2*ld+1
         n_l21m21 = (sph_shell_1)**2 !(la+lb+1)**2
         n_l43m43 = (sph_shell_2)**2 !(lc+ld+1)**2
         n_ij = max_len*max_len !lena*lenb
         n_kl = max_len*max_len !lenc*lend
         n_ijkl = n_ij*n_kl

         space_eri = sph_a*sph_b*sph_c*sph_d

         space_cpl_AB = max(sph_a*sph_b,sph_c*sph_d)*n_l21m21
         space_G = 2*(max(sph_shell_1,sph_shell_2))**2

         n_contr_pair = max(n_ij,n_kl)
         space_F_X = n_contr_pair*max(sph_a*sph_b,sph_c*sph_d)
         space_sum_AB = sph_a*sph_b*n_contr_pair*(la_p_lb+1)*(la_p_lb+2)*(2*la_p_lb+3)/6
         space_sum_CD = sph_c*sph_d*n_contr_pair*(lc_p_ld+1)*(lc_p_ld+2)*(2*lc_p_ld+3)/6
         space_cpl_ABCD = n_l21m21*n_l43m43*(sum_l+1)*n_ijkl
         space_contr_AB = sph_a*sph_b*n_kl*n_l43m43*(lc_p_ld+1)

         space_l_ijkl_L_coefficients = (2+sum_l*(sum_l+1))/2*n_ijkl

         !space for the temporary array that keeps all intermediate results before their indices are swapped/transposed into the desired order
         n_l21m21 = max(la_p_lb+1,lc_p_ld+1)
         n_l21m21 = (n_l21m21+1)*(n_l21m21+2)*(2*n_l21m21+3)/6

         n_n = sum_l/2+1
         n_l = max(sum_l-2,2*max(max_l_1,max_l_2))+1
         space_cf = n_l*n_n*(n_n+1)/2 !space for the Lag. coefficients

         !Allocate all at once - spatial locality of the arrays in memory (page faults)
         if (allocated(cf)) deallocate(cf,SH_l,SH_ml)
         allocate(G_A(space_G),G_B(space_G),F_X(space_F_X),sum_AB(space_sum_AB),contr_AB(space_contr_AB),sum_CD(space_sum_CD), &
                  cpl_ABCD(space_cpl_ABCD),l_ijkl_L_coefficients(space_l_ijkl_L_coefficients),cf(space_cf), &
                  SH_ml(2,(sum_l+1)**2),SH_l(sum_l+1),stat=err)
         if (err .ne. 0) stop "memory allocation error"

         !Precalculate the coefficients for the Lag. polynomials. cf_m(l,n)
         ind = 0
         do n=0,n_n-1
            do l=0,n_l-1
               f = 1.0_cfp
               do m=0,n
                 if (m .ne. 0) f = f*m
                 ind = ind + 1 !save in order (m,l,n)
                 cf(ind) = (-1)**m/f*gen_binom(n+l+0.5_cfp,n-m) !cf(m+1)
               enddo
            enddo !l
         enddo !n

         !Precalculate the coefficients needed to generate the Solid Harmonics
         ind = 0
         do l=0,sum_l
            SH_l(l+1) = sqrt((2.0_cfp*l+1)/(2.0_cfp*l+2))
            do m=-l,l
               ind = ind + 1
               SH_ml(1,ind) = 1.0_cfp/sqrt((l+1.0_cfp+m)*(l+1.0_cfp-m))
               SH_ml(2,ind) = sqrt((l+m)*(l-m)*1.0_cfp)*SH_ml(1,ind)
            enddo
         enddo

   end subroutine allocate_space

   !we assume that the shells have been already ordered so that la+lb .ge. lc+ld, la .ge. lb, lc .ge. ld
   subroutine eri_shell_sph (lena,xa,ya,za,anorms,la,aexps,acoefs, &
                             lenb,xb,yb,zb,bnorms,lb,bexps,bcoefs, &
                             lenc,xc,yc,zc,cnorms,lc,cexps,ccoefs, &
                             lend,xd,yd,zd,dnorms,ld,dexps,dcoefs,two_el_column,eri_ints)
      use phys_const, only: twopi, fourpi
      use const, only: imax_wp,boys_f_grid_step_wp,taylor_k_wp, imax_ep,boys_f_grid_step_ep,taylor_k_ep,mmax
      implicit none
      real(kind=cfp), parameter :: prefactor = 32*(twopi)**(6.5_cfp)
      real(kind=cfp), allocatable :: eri_ints(:,:)
      integer, intent(in) :: lena, lenb, lenc, lend, la, lb, lc, ld, two_el_column
      real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb, xc,yc,zc, xd,yd,zd
      real(kind=cfp), intent(in) :: anorms(:),aexps(:),acoefs(:)
      real(kind=cfp), intent(in) :: bnorms(:),bexps(:),bcoefs(:)
      real(kind=cfp), intent(in) :: cnorms(:),cexps(:),ccoefs(:)
      real(kind=cfp), intent(in) :: dnorms(:),dexps(:),dcoefs(:)

      integer :: i, j, sum_l, range_start, range_end, l, number_of_Fm, ind, ij, kl, n_ij, n_kl, big_l, n, block, n_ijkl
      real(kind=cfp) :: adj_norms_a(lena), adj_norms_b(lenb), adj_norms_c(lenc), adj_norms_d(lend), r_abcd(3), r_ba(3), r_dc(3), &
                        r_alpha_beta(3,lena*lenb), r_gamma_delta(3,lenc*lend)
      real(kind=cfp) :: eta_ab(lena*lenb), eta_cd(lenc*lend), T(lena*lenb*lenc*lend), exp_T(lena*lenb*lenc*lend), &
                        eta(lena*lenb*lenc*lend), Fm_contracted((la+lb+lc+ld+1)*lena*lenb*lenc*lend)
      real(kind=cfp) :: solid_harmonics_ab((la+lb+1)**2), solid_harmonics_cd((lc+ld+1)**2), &
                        solid_harmonics_abcd((la+lb+lc+ld+1)**2*lena*lenb*lenc*lend), sq_dist_ba, sq_dist_dc, sq_dist_abcd
      real(kind=cfp) :: inverse, data_AB(1:6,lena*lenb), data_CD(1:6,lenc*lend), xsi, factor
      logical :: rC_eq_rD, rA_eq_rB, is_cd, atomic

         sum_l = la+lb+lc+ld

         !precalculate all G factors and Gaunt coefficients
         call cpl%prec_G_cf(sum_l)

         !allocate space for all intermediate arrays
         call allocate_space(la,lb,lc,ld,lena,lenb,lenc,lend)

!         write(*,'("Ls",4i)') la,lb,lc,ld
!         write(*,'("len",4i)') lena,lenb,lenc,lend
!         write(*,'("RA",3e)') xa,ya,za
!         write(*,'("RB",3e)') xb,yb,zb
!         write(*,'("RC",3e)') xc,yc,zc
!         write(*,'("RD",3e)') xd,yd,zd
!         write(*,'("a norms",100f)') anorms(1:lena)
!         write(*,'("b norms",100f)') bnorms(1:lenb)
!         write(*,'("c norms",100f)') cnorms(1:lenc)
!         write(*,'("d norms",100f)') dnorms(1:lend)
!
!         write(*,'("a exps",100f)') aexps(1:lena)
!         write(*,'("b exps",100f)') bexps(1:lenb)
!         write(*,'("c exps",100f)') cexps(1:lenc)
!         write(*,'("d exps",100f)') dexps(1:lend)
!
!         write(*,'("a coefs",100f)') acoefs(1:lena)
!         write(*,'("b coefs",100f)') bcoefs(1:lenb)
!         write(*,'("c coefs",100f)') ccoefs(1:lenc)
!         write(*,'("d coefs",100f)') dcoefs(1:lend)
!
         factor = sqrt(fourpi/(2*la+1.0_cfp))
         do i=1,lena
            adj_norms_a(i) = anorms(i)/((2*aexps(i))**(la+1.5_cfp))*factor
         enddo

         factor = sqrt(fourpi/(2*lb+1.0_cfp))
         do i=1,lenb
            adj_norms_b(i) = bnorms(i)/((2*bexps(i))**(lb+1.5_cfp))*factor
         enddo

         factor = sqrt(fourpi/(2*lc+1.0_cfp))
         do i=1,lenc
            adj_norms_c(i) = cnorms(i)/((2*cexps(i))**(lc+1.5_cfp))*factor
         enddo

         factor = sqrt(fourpi/(2*ld+1.0_cfp))
         do i=1,lend
            adj_norms_d(i) = dnorms(i)/((2*dexps(i))**(ld+1.5_cfp))*factor
         enddo

         !initialize the Boys function; if i == 1 then it has been initialized already. Note that we use the appropriate parameters (found using the test program test_boys_function_obj) depending
         !on the precision used.
         if (cfp .eq. wp) then
            i = boys%init(imax_wp,mmax,boys_f_grid_step_wp,taylor_k_wp) !mmax = sum_l
         elseif (cfp .eq. ep1) then
            i = boys%init(imax_ep,mmax,boys_f_grid_step_ep,taylor_k_ep) !mmax = sum_l
         else
            stop "eri_shell: unsupported numeric type"
         endif

         if (i .ne. 0 .and. i .ne. 1) then
            print *,i
            stop "eri_shell: initialization of the Boys function failed."
         endif

         !evaluate the AB solid harmonics
         r_ba(1:3) = (/xb-xa,yb-ya,zb-za/)
         sq_dist_ba = dot_product(r_ba,r_ba)

         if (xb .eq. xa .and. yb .eq. ya .and. zb .eq. za) then
            rA_eq_rB = .true. !in this case I need solid_harmonics only for l=0, i.e. sqrt(1/fourpi)
            solid_harmonics_ab(1) = sqrt(1.0_cfp/fourpi)
            solid_harmonics_ab(2:(la+lb+1)**2) = 0.0_cfp
         else
            rA_eq_rB = .false.
            call cfp_solh_1d(solid_harmonics_ab,r_ba(1),r_ba(2),r_ba(3),la+lb)
            !we need Slm = r**l*X(l,m)
            do l=0,la+lb
               range_start = l*l+1
               range_end = range_start-1 + 2*l+1
               solid_harmonics_ab(range_start:range_end) = solid_harmonics_ab(range_start:range_end)*sqrt((2*l+1.0_cfp)/fourpi)
            enddo
         endif

         !evaluate the CD solid harmonics
         r_dc(1:3) = (/xd-xc,yd-yc,zd-zc/)
         sq_dist_dc = dot_product(r_dc,r_dc)

         if (xd .eq. xc .and. yd .eq. yc .and. zd .eq. zc) then
            rC_eq_rD = .true. !in this case I need solid_harmonics only for l=0, i.e. sqrt(1/fourpi)
            solid_harmonics_cd(1) = sqrt(1.0_cfp/fourpi)
            solid_harmonics_cd(2:(lc+ld+1)**2) = 0.0_cfp
         else
            rC_eq_rD = .false.
            call cfp_solh_1d(solid_harmonics_cd,r_dc(1),r_dc(2),r_dc(3),lc+ld)
            !we need Slm = r**l*X(l,m)
            do l=0,lc+ld
               range_start = l*l+1
               range_end = range_start-1 + 2*l+1
               solid_harmonics_cd(range_start:range_end) = solid_harmonics_cd(range_start:range_end)*sqrt((2*l+1.0_cfp)/fourpi)
            enddo
         endif

         !Determine if this is just one-centre problem
         atomic = .false.
         if (rA_eq_rB .and. rC_eq_rD) then
            if (xa .eq. xc .and. ya .eq. yc .and. za .eq. zc) then
               atomic = .true.
               if (abs(la-lb) > lc+ld) then !Wigner-Eckart assuming la+lb .ge. lc+ld
                  eri_ints(1:(2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1),two_el_column) = 0.0_cfp
                  return
               endif
            endif
         endif

         !evaluate the AB product center and expand the product of the contraction coefficients
         n_ij = lena*lenb
         ij = 0
         do i=1,lena
            do j=1,lenb
               ij=ij+1
               eta_ab(ij) = aexps(i)+bexps(j)
               inverse = 1.0_cfp/eta_ab(ij)
   
               r_alpha_beta(1,ij) = (aexps(i)*xa+bexps(j)*xb)*inverse
               r_alpha_beta(2,ij) = (aexps(i)*ya+bexps(j)*yb)*inverse
               r_alpha_beta(3,ij) = (aexps(i)*za+bexps(j)*zb)*inverse

               xsi = aexps(i)*bexps(j)*inverse
               data_AB(1,ij) = acoefs(i)*bcoefs(j)*adj_norms_a(i)*adj_norms_b(j)
               data_AB(2,ij) = aexps(i)*inverse
               data_AB(3,ij) = bexps(j)*inverse
               data_AB(4,ij) = 1.0_cfp/(4.0_cfp*xsi)
               data_AB(5,ij) = exp(-xsi*sq_dist_ba)
               data_AB(6,ij) = xsi*sq_dist_ba
            enddo !j
         enddo !i

         !evaluate the CD product center and expand the product of the contraction coefficients
         n_kl = lenc*lend
         kl = 0
         do i=1,lenc
            do j=1,lend 
               kl=kl+1
               eta_cd(kl) = cexps(i)+dexps(j)
               inverse = 1.0_cfp/eta_cd(kl)
 
               r_gamma_delta(1,kl) = (cexps(i)*xc+dexps(j)*xd)*inverse
               r_gamma_delta(2,kl) = (cexps(i)*yc+dexps(j)*yd)*inverse
               r_gamma_delta(3,kl) = (cexps(i)*zc+dexps(j)*zd)*inverse

               xsi = cexps(i)*dexps(j)*inverse 
               data_CD(1,kl) = ccoefs(i)*dcoefs(j)*adj_norms_c(i)*adj_norms_d(j)
               data_CD(2,kl) = cexps(i)*inverse
               data_CD(3,kl) = dexps(j)*inverse
               data_CD(4,kl) = 1.0_cfp/(4.0_cfp*xsi)
               data_CD(5,kl) = exp(-xsi*sq_dist_dc)
               data_CD(6,kl) = xsi*sq_dist_dc
            enddo !j
         enddo !i

         factor = prefactor*(-1)**(lb+ld)
         n_ijkl = n_ij*n_kl
         ind = 0
         do kl=1,n_kl
            do ij=1,n_ij

               ind = ind + 1

               if (atomic) then
                  r_abcd(1:3) = 0.0_cfp
               else
                  r_abcd(1:3) = r_alpha_beta(1:3,ij) - r_gamma_delta(1:3,kl)
               endif
               eta(ind) = eta_ab(ij)*eta_cd(kl)/(eta_ab(ij)+eta_cd(kl))
               xsi = 1.0_cfp/(4.0_cfp*eta(ind))

               !calculate Fm(T) for m=0,...,mmax using the routine based on the Taylor expansion
               range_start = (ind-1)*(sum_l+1)+1
               range_end = range_start-1 + (sum_l+1)
               j = range_start !we'll need range_start in the loop over big_l below so we save it into j
               sq_dist_abcd = dot_product(r_abcd,r_abcd)
               T(ind) = eta(ind)*sq_dist_abcd != T
               exp_T(ind) = exp(-T(ind))
               number_of_Fm = range_end-range_start+1
               call boys%eval_taylor(Fm_contracted(range_start:range_end),number_of_Fm,T(ind),sum_l)

               !evaluate the solid harmonics
               range_start = (ind-1)*(sum_l+1)**2+1
               range_end = range_start-1 + (sum_l+1)**2
               call cfp_solh_1d(solid_harmonics_abcd(range_start:range_end),r_abcd(1),r_abcd(2),r_abcd(3),sum_l)
               !we need Slm = r**l*X(l,m)
               range_start = range_start + 1 !to ensure l=0 case is handled properly: 2*(l-1)+1 for l=0 gives -1
               do big_l=0,sum_l
                  range_start = range_start + 2*(big_l-1)+1
                  range_end = range_start-1 + 2*big_l+1
                  solid_harmonics_abcd(range_start:range_end) &
                    = solid_harmonics_abcd(range_start:range_end)*sqrt((2*big_l+1.0_cfp)/fourpi)

                  block = (2+(big_l-1)*big_l)/2    !the number of l values for all L < big_l; for each L there are max(1,L) l-values.
                  if (big_l .eq. 0) block = 0

                  !The coefficients l_ijkl_L_coefficients are saved in the order (l,ij,kl,L) where the size of the first dimension is given by: max(1,L)
                  i = max(1,big_l) !the number of l values corresponding to the current L
                  i = block*n_ijkl + i*(ij-1) + i*n_ij*(kl-1) !offset for the current block (l_min,ij,kl,L)
                 
                  do l=0,big_l-2
                     i = i + 1
                     n = (big_l-l)/2
                     l_ijkl_L_coefficients(i) = (-1)**n*cnla(n-1,l,xsi)*exp_T(ind)*Lag_n_hlf_k(n-1,l,T(ind))*factor
!                     if (l .eq. 0) write(*,'("true cfs",i,4i4,e)') i,l,ij,kl,big_l,l_ijkl_L_coefficients(i)
                  enddo !l

                  l = big_l
                  i = i + 1
                  n = -1
                  l_ijkl_L_coefficients(i) = cnla(n,l,xsi)*Fm_contracted(j + big_l)*factor
!                  if (l .eq. 0) write(*,'("true cfs",i,4i4,e)') i,l,ij,kl,big_l,l_ijkl_L_coefficients(i)
               enddo !big_l

            enddo !ij
         enddo !kl

         !AB pair
         !result in sum_AB((ij),(m21,l21),(ma,mb),lap_lbp)
         is_cd = .false.
         call sum_over_map_mbp_contract_i_j(la,lb,n_ij,solid_harmonics_ab,data_AB,sum_AB,is_cd,rA_eq_rB,lc+ld)
    
         !CD pair
         !result in sum_CD((kl),(m43,l43),(mc,md),lcp_ldp)
         is_cd = .true.
         call sum_over_map_mbp_contract_i_j(lc,ld,n_kl,solid_harmonics_cd,data_CD,sum_CD,is_cd,rC_eq_rD,lc+ld)
 
         !ABCD contraction coefficients
         !result in cpl_ABCD((kl),(m43,l43),(ij),(m21,l21),L)
         call calculate_ABCD_coefficients(la+lb,lc+ld,n_ij,n_kl,l_ijkl_L_coefficients,solid_harmonics_abcd,cpl_ABCD)
    
         !Final step: contract the AB,CD pairs together
         !result in eri_ints(ma,mb,mc,md)
         call contract_AB_CD(la,lb,n_ij,sum_AB,lc,ld,n_kl,sum_CD,cpl_ABCD,eri_ints,two_el_column)

   end subroutine eri_shell_sph

   subroutine contract_AB_CD(la,lb,n_ij,sum_AB,lc,ld,n_kl,sum_CD,cpl_ABCD,eri_ints,two_el_column)
      implicit none
      integer, intent(in) :: la,lb,n_ij,lc,ld,n_kl,two_el_column
      real(kind=cfp), intent(in) :: sum_AB(:), sum_CD(:), cpl_ABCD(:)
      real(kind=cfp), allocatable :: eri_ints(:,:)

      integer :: lc_p_ld, la_p_lb, L, base_sum_AB, base_cpl_ABCD, ma_mb, mc_md, stride, sph_shell_ab, sph_shell_cd, &
                 n_l21m21_preceeding, n_l21m21, n_l43m43, n_l43m43_preceeding, base_sum_CD
      integer :: ind, n_l21m21_total, base_contr_AB

         sph_shell_ab = (2*la+1)*(2*lb+1)
         stride = (lc+ld+1)**2*n_kl*(la+lb+1)**2*n_ij
         sph_shell_cd = (2*lc+1)*(2*ld+1)
         n_l21m21_total = (la+lb+1)**2*n_ij

         eri_ints(1:sph_shell_ab*sph_shell_cd,two_el_column) = 0.0_cfp

         lc_p_ld = lc+ld+1
         n_l43m43_preceeding = lc_p_ld*(lc_p_ld+1)*(2*lc_p_ld+1)/6
         base_contr_AB = n_kl*sph_shell_ab*n_l43m43_preceeding
         contr_AB(1:base_contr_AB) = 0.0_cfp

         do lc_p_ld=0,lc+ld
            n_l43m43 = (lc_p_ld+1)**2*n_kl
            n_l43m43_preceeding = lc_p_ld*(lc_p_ld+1)*(2*lc_p_ld+1)/6
            base_contr_AB = n_kl*sph_shell_ab*n_l43m43_preceeding
            do la_p_lb=0,la+lb
               L = la_p_lb+lc_p_ld
               n_l21m21 = (la_p_lb+1)**2*n_ij
               n_l21m21_preceeding = la_p_lb*(la_p_lb+1)*(2*la_p_lb+1)/6
               base_sum_AB = n_ij*sph_shell_ab*n_l21m21_preceeding
               base_cpl_ABCD = L*stride
!               ind = base_contr_AB
!               do ma_mb=1,sph_shell_ab !(ma,mb)
!                  base_cpl_ABCD = L*stride
!                  do i=1,n_l43m43 !(kl,l43m43)
!                     ind = ind + 1
!                     contr_AB(ind) = contr_AB(ind) + sum(sum_AB(base_sum_AB+1:base_sum_AB+n_l21m21)*cpl_ABCD(base_cpl_ABCD+1:base_cpl_ABCD+n_l21m21)) !save in order contr_AB(kl,l43m43,ma,mb,lc_p_ld)
!                     base_cpl_ABCD = base_cpl_ABCD + n_l21m21_total
!                  enddo !i
!                  base_sum_AB = base_sum_AB + n_l21m21
!               enddo !ma_mb
               call mat_T_mat_mul_special_blocking(cpl_ABCD,sum_AB,contr_AB,n_l43m43,sph_shell_ab,n_l21m21,n_l21m21_total,&
                                                    base_cpl_ABCD,base_sum_AB,base_contr_AB)
            enddo !la_p_lb
         enddo !lc_p_ld

         base_contr_AB = 0
         do lc_p_ld=0,lc+ld
            n_l43m43 = (lc_p_ld+1)**2*n_kl
            n_l43m43_preceeding = lc_p_ld*(lc_p_ld+1)*(2*lc_p_ld+1)/6
            base_sum_CD = n_kl*sph_shell_cd*n_l43m43_preceeding
            do ma_mb=1,sph_shell_ab
               ind = base_sum_CD
               do mc_md=1,sph_shell_cd
                  eri_ints(ma_mb+sph_shell_ab*(mc_md-1),two_el_column) = eri_ints(ma_mb+sph_shell_ab*(mc_md-1),two_el_column) &
                            + sum(contr_AB(base_contr_AB+1:base_contr_AB+n_l43m43)*sum_CD(ind+1:ind+n_l43m43))
                  ind = ind + n_l43m43
               enddo !mc_md
               base_contr_AB = base_contr_AB + n_l43m43
            enddo !ma_mb
         enddo !lc_p_ld

   end subroutine contract_AB_CD

   !> Form C := A**T*B (+ C)
   subroutine mat_T_mat_mul_special_blocking(A,B,C,m,n,k,stride_a,A_base,B_base,C_base)
      use const, only: tile
      implicit none
      integer, intent(in) :: m,n,k, A_base,B_base,C_base, stride_a
      real(kind=cfp), intent(in) :: A(:), B(:)
      real(kind=cfp), intent(inout) :: C(:)

      integer :: i,j,l, c_column, b_column, a_column, i_tile,i_end,j_tile,j_end !, l_start,l_end
      real(kind=cfp) :: temp

          do i_tile=1,m,tile
             i_end = min(i_tile+tile-1,m)

             do j_tile=1,n,tile
                j_end = min(j_tile+tile-1,n)

!                do l_start=1,k,tile
!                   l_end = min(l_start+tile-1,k)
                   do j=j_tile,j_end
                      b_column = k*(j-1) + B_base
                      c_column = m*(j-1) + C_base
                      do i=i_tile,i_end
                         a_column = stride_a*(i-1) + A_base !note that a_column has stride stride_a and not k.
                         temp = 0.0_cfp
                         do l = 1,k !l_start,l_end
                            temp = temp + a(l+a_column)*b(l+b_column) !a(l,i)*b(l,j)
                         enddo !l
                         c(i+c_column) = temp + c(i+c_column) !c(i,j)
                      enddo !i
                   enddo !j
!                enddo !l_start

             enddo !j_tile
          enddo !i_tile

   end subroutine mat_T_mat_mul_special_blocking

   !> Form C := A**T*B (+ C)
   subroutine mat_T_mat_mul_special(A,B,C,m,n,k,stride_a,A_base,B_base,C_base)
      implicit none
      integer, intent(in) :: m,n,k, A_base,B_base,C_base, stride_a
      real(kind=cfp), intent(in) :: A(:), B(:)
      real(kind=cfp), intent(inout) :: C(:)

      integer :: i,j,l, c_column, b_column, a_column
      real(kind=cfp) :: temp

          do j = 1,n
             b_column = k*(j-1) + B_base
             c_column = m*(j-1) + C_base
             do i = 1,m
                a_column = stride_a*(i-1) + A_base
                temp = 0.0_cfp
                do l = 1,k
                   temp = temp + a(l+a_column)*b(l+b_column) !a(l,i)*b(l,j)
                enddo !l
                c(i+c_column) = temp + c(i+c_column)
             enddo !i
          enddo !j

   end subroutine mat_T_mat_mul_special

   subroutine calculate_F_ij(la,lap,lb,lbp,n_ij,data_AB,solid_harmonics_ab,is_cd,rA_eq_rB,F_X)
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), parameter :: inv_fourpi_sq = 1.0_cfp/(fourpi*fourpi)
     integer, intent(in) :: la,lb,lap,lbp,n_ij
     logical, intent(in) :: is_cd,rA_eq_rB
     real(kind=cfp), intent(in) :: solid_harmonics_ab(:), data_AB(1:6,n_ij)
     real(kind=cfp), intent(out) :: F_X(:)

     real(kind=cfp) :: sum_over_m1((2*la+1)*(2*lb+1)), s, pow_AB(n_ij)
     integer :: ij,i,l1,l1_min,l1_max,n1,sgn,m1_A,m1_B,ind,n_mapp_mbpp,base,mapp,mbpp,lapp,lbpp,p

              lapp = la-lap; lbpp = lb-lbp
              n_mapp_mbpp = (2*lapp+1)*(2*lbpp+1)
              l1_min = abs(lapp-lbpp); l1_max = lapp+lbpp

              !Calculate the lap,lbp -dependent part of the F_ij terms
              if (is_cd) then !swap the order of the powers lap,lbp
                 do ij=1,n_ij
                    pow_AB(ij) = data_AB(1,ij)*(-data_AB(3,ij))**lbp*data_AB(2,ij)**lap*data_AB(5,ij)*inv_fourpi_sq
                 enddo
              else
                 do ij=1,n_ij
                    pow_AB(ij) = data_AB(1,ij)*(-data_AB(2,ij))**lap*data_AB(3,ij)**lbp*data_AB(5,ij)*inv_fourpi_sq
                 enddo
              endif

              F_X(1:n_ij*n_mapp_mbpp) = 0.0_cfp
              i = (l1_min-1 +1)**2 !skip the spherical harmonics for l1=0,...,l1_min-1
              do l1=l1_min,l1_max,2

                 n1 = (l1_max-l1)/2
                 p = (-1)**n1

                 if (rA_eq_rB .and. l1 > 0) exit !For Rba = 0 only the l1 spherical harmonic is non-zero

                 !For each mapp,mbpp combination of angular indices sum over the m1 angular indices and save the non-trivial mapp,mbpp values
                 sum_over_m1 = 0.0_cfp
                 do mbpp=-lbpp,lbpp
                    base = (2*lapp+1)*(lbpp+mbpp) !base for the full 2D index: lapp+mapp+1 + (2*lapp+1)*(lbpp+mbpp) of the element (lapp+mapp+1,lbpp+mbpp+1)
                    do mapp=-lapp,lapp

                       s = 0.0_cfp
                       !determine the m1 values that can yield a non-trivial Gaunt coefficient and sum over the corresp. terms.
                       sgn = sign(1,mbpp)*sign(1,mapp) !if one of mapp,mbpp is zero its sign is positive according to the Fortran std.
                       m1_A = sgn*abs(mapp+mbpp)
                       if (abs(m1_A) .le. l1) then 
                          s = cpl%rgaunt(lbpp,lapp,l1,mbpp,mapp,m1_A)*solid_harmonics_ab(i +m1_A+l1+1)
                       endif
                       m1_B = sgn*abs(mbpp-mapp)
                       if (abs(m1_B) .le. l1 .and. m1_B .ne. m1_A) then 
                          s = s + cpl%rgaunt(lbpp,lapp,l1,mbpp,mapp,m1_B)*solid_harmonics_ab(i +m1_B+l1+1)
                       endif

                       sum_over_m1(lapp+mapp+1 + base) = sum_over_m1(lapp+mapp+1 + base) + s*p

                    enddo !mapp
                 enddo !mbpp
                 i = i + (2*l1+1) + (2*(l1+1)+1) !skip the l1+1 harmonics (we loop over l1 in steps of 2)

                 !Sum into F_X(mapp,mbpp,ij) for each combination of mapp,mbpp indices the F_ij factors multiplied by the solid harmonics Y_{l1,m1}*<lapp,mapp|lbpp,mbpp|l1,m1>
                 ind = 0
                 do ij=1,n_ij
                    F_X(ind+1:ind+n_mapp_mbpp) = F_X(ind+1:ind+n_mapp_mbpp) &
                            + sum_over_m1(1:n_mapp_mbpp)*pow_AB(ij)*cnla(n1,l1,data_AB(4,ij))*Lag_n_hlf_k(n1,l1,data_AB(6,ij))
                    ind = ind + n_mapp_mbpp
                 enddo !ij
                 
              enddo !l1

   end subroutine calculate_F_ij

   subroutine sum_over_map_mbp_contract_i_j(la,lb,n_ij,solid_harmonics_ab,data_AB,sum_AB,is_cd,rA_eq_rB,lc_p_ld)
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), parameter :: inv_fourpi_sq = 1.0_cfp/(fourpi*fourpi)
     integer, intent(in) :: la,lb,n_ij,lc_p_ld
     logical, intent(in) :: is_cd, rA_eq_rB
     real(kind=cfp), intent(in) :: solid_harmonics_ab((la+lb+1)**2), data_AB(1:6,n_ij)
     real(kind=cfp), intent(out) :: sum_AB(:)

     integer :: ma, lap, lbp, mb, map, mbp, l21,m21, i, n_mapp, n_mbpp, lapp, lbpp, n_mapp_mbpp, l_min, l_max, ij, ind, l21m21
     integer :: base, n_l21m21, stride_A, stride_B, stride_AB, base_AB, ma_mb, base_A, base_B
     real(kind=cfp) :: coupling, G_mbp_mbpp(2)
     real(kind=cfp) :: prod_G_A_G_B_coupling(4), ij_ma_mb(n_ij*(2*la+1)*(2*lb+1))
     integer :: sph_ab, d1, d2, d3, m1_A, m1_B, sgn, sph_a, sph_b, mapp_list(2*(2*la+1)**2), mbpp_list(2*(2*lb+1)**2), &
                mapp_mbpp_list(4)

        sph_a = 2*la+1
        sph_b = 2*lb+1
        sph_ab = sph_a*sph_b
        stride_A = sph_a*2
        stride_B = sph_b*2

        base_AB = la+lb+1
        base_AB = (base_AB)*(base_AB+1)*(2*base_AB+1)/6 !total number of l21,m21 pairs for l21=0,1,...,lap+lbp-1
        base_AB = sph_ab*n_ij*base_AB !total number of preceeding combinations (m21,l21,ma,mb,ij,L) for L=0,1,...,lap+lbp-1
        sum_AB(1:base_AB) = 0.0_cfp 
        do lbp=0,lb
            lbpp = lb-lbp
            n_mbpp = 2*lbpp+1

           !Load the G coefficients G(lb,mb|lbp,mbp|mbpp) for the B centre into a linear array in the order (mb,mbpp,mbp)
           ind = 0
           do mbp=-lbp,lbp
              do mb=-lb,lb
                 G_B(ind+1:ind+2) = 0.0_cfp !initialize the coefficients corresponding to the two possible mapp values to zero.
                 mbpp_list(ind+1:ind+2) = 1

                 sgn = sign(1,mbp)*sign(1,mb) !if one of mapp,mbpp is zero its sign is positive according to the Fortran std.
                 m1_A = sgn*abs(mbp+mb) !mbpp_1
                 m1_B = sgn*abs(mbp-mb) !mbpp_2
                 if (abs(m1_A) .le. lbpp) then
                    G_B(ind+1) = cpl%G_real_cf(lb,lbp,mb,mbp,m1_A) !save in order (mapp,ma,map)
                    mbpp_list(ind+1) = m1_A+lbpp+1
                 endif
                 if (abs(m1_B) .le. lbpp .and. m1_B .ne. m1_A) then
                    G_B(ind+2) = cpl%G_real_cf(lb,lbp,mb,mbp,m1_B) !save in order (mapp,ma,map)
                    mbpp_list(ind+2) = m1_B+lbpp+1
                 endif
                 ind = ind + 2
              enddo
           enddo

           do lap=0,la
              lapp = la-lap
              n_mapp = 2*lapp+1

              !If Rba = 0 then only those lap,lbp pairs contribute which give non-zero <lbpp,mbpp|lapp,mapp|0,0> => (lapp,mapp) = (lbpp,mbpp)
              if (rA_eq_rB .and. lapp .ne. lbpp) cycle

              n_mapp_mbpp = n_mapp*n_mbpp !total number of (mapp,mbpp) pairs

              base_AB = lap+lbp
              base_AB = (base_AB)*(base_AB+1)*(2*base_AB+1)/6 !total number of l21,m21 pairs for l21=0,1,...,lap+lbp-1
              base_AB = sph_ab*n_ij*base_AB !total number of preceeding combinations (m21,l21,ma,mb,ij,L) for L=0,1,...,lap+lbp-1

!
!----- Calculate F_ij terms
!
              call calculate_F_ij(la,lap,lb,lbp,n_ij,data_AB,solid_harmonics_ab,is_cd,rA_eq_rB,F_X)
!
!-----
!
              !Load the G coefficients G(la,ma|lap,map|mapp) for the A centre into a linear array in the order (ma,mapp,map)
              !todo load them only if lbp=0 otherwise save them and save the offsets for each lap as well so I can reuse them in the following iterations.
              ind = 0
              do map=-lap,lap
                 do ma=-la,la
                    G_A(ind+1:ind+2) = 0.0_cfp !initialize the coefficients corresponding to the two possible mapp values to zero.
                    mapp_list(ind+1:ind+2) = 1

                    sgn = sign(1,map)*sign(1,ma) !if one of mapp,mbpp is zero its sign is positive according to the Fortran std.
                    m1_A = sgn*abs(map+ma) !mapp_1
                    m1_B = sgn*abs(map-ma) !mapp_2
                    if (abs(m1_A) .le. lapp) then
                       G_A(ind+1) = cpl%G_real_cf(la,lap,ma,map,m1_A) !save in order (mapp,ma,map)
                       mapp_list(ind+1) = m1_A+lapp+1
                    endif
                    if (abs(m1_B) .le. lapp .and. m1_B .ne. m1_A) then
                       G_A(ind+2) = cpl%G_real_cf(la,lap,ma,map,m1_B) !save in order (mapp,ma,map)
                       mapp_list(ind+2) = m1_B+lapp+1
                    endif
                    ind = ind + 2
                 enddo
              enddo

              !Perform for each (ma,mb,ij,m21,l21): sum (over map,mbp,mapp,mbpp) G(lb,mb|lbp,mbp,mbpp)*<lbp,mbp|l21,m21|lap,map>*G(la,ma|lap,map,mapp)*F_ij(lap,lbp)
              l_min = abs(lap-lbp)
              l_max = lap+lbp
              stride_AB = sph_a*(2*lapp+1)*sph_b*(2*lbpp+1)
              n_l21m21 = (l_max+1)**2
              d1 = (2*lapp+1)
              d2 = d1*(2*lbpp+1)
              d3 = d2*(2*la+1)
              do l21=l_min,l_max,2
                 do m21=-l21,l21

                    l21m21 = l21**2 + l21+m21+1
                    ij_ma_mb = 0.0_cfp
                    !Here we perform the double summation over the map,mbp,mapp,mbpp indices while multiplying by the F_ij factors that depend on the primitive parameters:
                    !sum (over map,mbp) sum (over mapp,mbpp) G(lb,mb|lbp,mbp,mbpp)*<lbp,mbp|l21,m21|lap,map>*G(la,ma|lap,map,mapp)*F_ij(mapp,mbpp)
                    do mbp=-lbp,lbp
                      
                       !todo For each (lap,lbp) pair it is possible to use the symmetry of <lbp,mbp|l21,m21|lap,map> to calculate at once the contribution of the pair (lbp,lap) if I also calculate 
                       !F_ij(lbp,lap). This way I'll reduce approx. by 1/2 the number of rgaunt evaluated -for (LL|LcLd) class. Generally, this method can be applied for those lap not greater than lbp.
                       !todo I can also use the Gaunt cf. generated for the A,B pair to generate sum_CD for the CD pair since the set of Gaunts required for the CD pair is a subset of those needed for the
                       !AB pair. However, if parallelization is required then maybe it is better to split the AB,CD pair coefficient generation?

                       !Determine the at most 2 possible values of map: m1_A, m1_B.
                       sgn = sign(1,mbp)*sign(1,m21) !rember that if one of mapp,mbpp is zero its sign is positive according to the Fortran std. so the tests below are correct even for mpp = 0.
                       m1_A = sgn*abs(mbp+m21) !map_1
                       m1_B = sgn*abs(mbp-m21) !map_2
                       if (abs(m1_A) .le. lap) then
                          coupling = cpl%rgaunt(lbp,l21,lap,mbp,m21,m1_A) != <lbp,mbp|l21,m21|lap,map>
                          if (coupling .ne. 0.0_cfp) then
                             base_B = stride_B*(lbp+mbp)
                             !todo doing mapp,mbpp summations here is not good - go back to the model of doing map first, then map, and then multiply by mapp,mbpp. However, I need to save the non-trivial pp.
                             ind = 0
                             do mb=1,sph_b
                                G_mbp_mbpp(1:2) = (/G_B(base_B+1)*coupling,G_B(base_B+2)*coupling/)
                                base_A = stride_A*(lap+m1_A)
                                if (G_mbp_mbpp(1) .eq. 0.0_cfp .and. G_mbp_mbpp(2) .eq. 0.0_cfp) then 
                                   base_A = base_A + 2*sph_a
                                   base_B = base_B + 2
                                   ind = ind + sph_a*n_ij
                                   cycle
                                endif
                                do ma=1,sph_a
                                   mapp_mbpp_list(1) = mapp_list(base_A+1) + d1*(mbpp_list(base_B+1)-1)
                                   mapp_mbpp_list(2) = mapp_list(base_A+1) + d1*(mbpp_list(base_B+2)-1)
                                   mapp_mbpp_list(3) = mapp_list(base_A+2) + d1*(mbpp_list(base_B+1)-1)
                                   mapp_mbpp_list(4) = mapp_list(base_A+2) + d1*(mbpp_list(base_B+2)-1)
                                   prod_G_A_G_B_coupling(1) = G_A(base_A+1)*G_mbp_mbpp(1)
                                   prod_G_A_G_B_coupling(2) = G_A(base_A+1)*G_mbp_mbpp(2)
                                   prod_G_A_G_B_coupling(3) = G_A(base_A+2)*G_mbp_mbpp(1)
                                   prod_G_A_G_B_coupling(4) = G_A(base_A+2)*G_mbp_mbpp(2)
                                   do ij=1,n_ij
                                      ij_ma_mb(ind+ij) = ij_ma_mb(ind+ij) &
                                            + prod_G_A_G_B_coupling(1)*F_X(mapp_mbpp_list(1) + n_mapp_mbpp*(ij-1))&
                                            + prod_G_A_G_B_coupling(2)*F_X(mapp_mbpp_list(2) + n_mapp_mbpp*(ij-1))&
                                            + prod_G_A_G_B_coupling(3)*F_X(mapp_mbpp_list(3) + n_mapp_mbpp*(ij-1))&
                                            + prod_G_A_G_B_coupling(4)*F_X(mapp_mbpp_list(4) + n_mapp_mbpp*(ij-1))
                                   enddo !ij
                                   ind = ind + n_ij
                                   base_A = base_A + 2
                                enddo !ma
                                base_B = base_B + 2
                             enddo !mb
                          endif
                       endif
                       if (abs(m1_B) .le. lap .and. m1_B .ne. m1_A) then 
                          coupling = cpl%rgaunt(lbp,l21,lap,mbp,m21,m1_B)
                          if (coupling .ne. 0.0_cfp) then
                             base_B = stride_B*(lbp+mbp)
                             ind = 0
                             do mb=1,sph_b
                                G_mbp_mbpp(1:2) = (/G_B(base_B+1)*coupling,G_B(base_B+2)*coupling/)
                                base_A = stride_A*(lap+m1_B)
                                if (G_mbp_mbpp(1) .eq. 0.0_cfp .and. G_mbp_mbpp(2) .eq. 0.0_cfp) then
                                   base_A = base_A + 2*sph_a
                                   base_B = base_B + 2
                                   ind = ind + sph_a*n_ij
                                   cycle
                                endif
                                do ma=1,sph_a
                                   mapp_mbpp_list(1) = mapp_list(base_A+1) + d1*(mbpp_list(base_B+1)-1)
                                   mapp_mbpp_list(2) = mapp_list(base_A+1) + d1*(mbpp_list(base_B+2)-1)
                                   mapp_mbpp_list(3) = mapp_list(base_A+2) + d1*(mbpp_list(base_B+1)-1)
                                   mapp_mbpp_list(4) = mapp_list(base_A+2) + d1*(mbpp_list(base_B+2)-1)
                                   prod_G_A_G_B_coupling(1) = G_A(base_A+1)*G_mbp_mbpp(1)
                                   prod_G_A_G_B_coupling(2) = G_A(base_A+1)*G_mbp_mbpp(2)
                                   prod_G_A_G_B_coupling(3) = G_A(base_A+2)*G_mbp_mbpp(1)
                                   prod_G_A_G_B_coupling(4) = G_A(base_A+2)*G_mbp_mbpp(2)
                                   do ij=1,n_ij
                                      ij_ma_mb(ind+ij) = ij_ma_mb(ind+ij) &
                                            + prod_G_A_G_B_coupling(1)*F_X(mapp_mbpp_list(1) + n_mapp_mbpp*(ij-1))&
                                            + prod_G_A_G_B_coupling(2)*F_X(mapp_mbpp_list(2) + n_mapp_mbpp*(ij-1))&
                                            + prod_G_A_G_B_coupling(3)*F_X(mapp_mbpp_list(3) + n_mapp_mbpp*(ij-1))&
                                            + prod_G_A_G_B_coupling(4)*F_X(mapp_mbpp_list(4) + n_mapp_mbpp*(ij-1))
                                   enddo !ij
                                   ind = ind + n_ij
                                   base_A = base_A + 2
                                enddo !ma
                                base_B = base_B + 2
                             enddo !mb
                          endif
                       endif 
                    enddo !mbp

                    !sum_AB((ij),(m21,l21),(ma,mb),L)
                    i = 0
                    do ma_mb=1,sph_ab
                       do ij=1,n_ij
                          i = i + 1
                          ind = base_AB + ij + n_ij*(l21m21-1) + n_l21m21*n_ij*(ma_mb-1) !base_AB + l21m21 + n_l21m21*(ij-1) + n_l21m21*n_ij*(ma_mb-1)
                          sum_AB(ind) = sum_AB(ind) + ij_ma_mb(i)
                       enddo !ij
                    enddo !ma_mb

                    !todo typically only less than 50% of sum_AB elements are non-zero...
                 enddo !m21
              enddo !l21
            enddo !lap
         enddo !lbp

   end subroutine sum_over_map_mbp_contract_i_j

   !> Form C := A*B (+ C)
   subroutine mat_mat_mul(A,B,C,m,n,k,add_to_C,C_base)
      implicit none
      integer, intent(in) :: m,n,k, C_base
      real(kind=cfp), intent(in) :: A(:), B(:)
      real(kind=cfp), intent(inout) :: C(:)
      logical, intent(in) :: add_to_C

      integer :: i,j,l, c_column, b_column, a_column

          do j = 1,n
             b_column = k*(j-1)
             c_column = m*(j-1) + C_base
             if (.not.(add_to_c)) then
                do i = 1,m
                   c(i + c_column) = 0.0_cfp !c(i,j)
                enddo
             end if
             do l = 1,k
                a_column = m*(l-1)
                if (b(l + b_column) .ne. 0.0_cfp) then !b(l,j)
                   do i = 1,m
                      c(i + c_column) = c(i + c_column) + b(l + b_column)*a(i + a_column) !c(i,j) = c(i,j) + b(l + b_column)*a(i,l)
                   enddo !i
                end if
             enddo !l
          enddo !j

   end subroutine mat_mat_mul

   !> Form C := A**T*B (+ C)
   subroutine mat_T_mat_mul(A,B,C,m,n,k,add_to_C,C_base)
      implicit none
      integer, intent(in) :: m,n,k, C_base
      real(kind=cfp), intent(in) :: A(:), B(:)
      real(kind=cfp), intent(inout) :: C(:)
      logical, intent(in) :: add_to_C

      integer :: i,j,l, c_column, b_column, a_column
      real(kind=cfp) :: temp

          do j = 1,n
             b_column = k*(j-1)
             c_column = m*(j-1) + C_base
             do i = 1,m
                a_column = k*(i-1)
                temp = 0.0_cfp
                do l = 1,k
                   temp = temp + a(l+a_column)*b(l+b_column) !a(l,i)*b(l,j)
                enddo !l
                if (.not.(add_to_C)) then
                   c(i+c_column) = temp !c(i,j)
                else
                   c(i+c_column) = temp + c(i+c_column)
                end if
             enddo !i
          enddo !j

   end subroutine mat_T_mat_mul

   !We assume here that la_p_lb .ge. lc_p_ld
   !todo if Rabcd = 0 (=> Rb=0) then only (l21,m21) = (l43,m43) contributes
   subroutine calculate_ABCD_coefficients(la_p_lb,lc_p_ld,n_ij,n_kl,l_ijkl_L_coefficients,solid_harmonics_abcd,cpl_ABCD)
       implicit none
       integer, intent(in) :: la_p_lb,lc_p_ld,n_ij,n_kl
       real(kind=cfp), intent(in) :: l_ijkl_L_coefficients(:), solid_harmonics_abcd(:)
       real(kind=cfp), intent(out) :: cpl_ABCD(:)

       integer :: l21, m21, l43, m43, big_l, l,m, ij, kl, base_solid_harmonics, base, ind, l_min, l_max, n_ijkl, &
                  terms_to_skip, base_l, ind_sum, ijkl, n_l, base_cf, n_l_total, i
       integer :: n_l21m21, n_l43m43, total_l, lm2, m21_min, lm_total, stride, m_A, m_B
       real(kind=cfp) :: sum_over_m((la_p_lb+lc_p_ld+1)*n_ij*n_kl), coupling
       integer :: n_harmonics, cnt, offset, l21m21_rank, l43m43_rank, l21_base, l43_base, sgn
       integer :: dAB, dBA

          n_l_total = max(1,la_p_lb+lc_p_ld)
          n_ijkl = n_ij*n_kl
          n_harmonics = (la_p_lb+lc_p_ld+1)**2
          n_l43m43 = (lc_p_ld+1)**2
          n_l21m21 = (la_p_lb+1)**2
          lm_total = n_l43m43*n_l21m21*n_ijkl !(n_l43m43+2*n_l43m43*n_l21m21-n_l43m43*n_l43m43)/2*n_ijkl !total number of (l21,m21)(ij)(l43m43)(kl) indices
          do l43=0,lc_p_ld
             l43_base = l43*l43+l43+1
             do m43=-l43,l43
                l43m43_rank = l43_base+m43 !sequence number of the pair (l43,m43)
                !loop only over the unique combinations of (l21,m21),(l43,m43) since the real Gaunt is symmetrical with respect to the 21,43 indices and therefore is cpl_ABCD.
                do l21=l43,la_p_lb
                   l21_base = l21*l21+l21+1

                   !Determine the max allowed l-value for the Gaunt coefficients <l,m|l21,m21|l43,m43>
                   l_max = l21+l43

                   !Make sure l21m21_rank .ge. l43m43_rank
                   m21_min = -l21
                   if (l21 .eq. l43) m21_min = m43

                   do m21=m21_min,l21

                      l21m21_rank = l21_base+m21 !sequence number of the pair (l21,m21)

                      !Determine the m values that can yield a non-trivial Gaunt coefficient <l,m|l21,m21|l43,m43> and sum over the Gaunt coefficients multiplied by the spherical harmonics Y(l,m)_ijkl.
                      sgn = sign(1,m43)*sign(1,m21)
                      m_A = abs(m21+m43); m_B = abs(m21-m43)
                      l_min = max(l21-l43,min(m_A,m_B)) !assuming l21.ge.l43
                      if (mod(l_max+l_min,2) .ne. 0) l_min = l_min+1
                      total_l = l_max-l_min+1 !total number of different l-values for each i,j,k,l combination (this number may change below).
                      sum_over_m(1:total_l*n_ijkl) = 0.0_cfp

                      offset = max(l_min,m_A)
                      if (mod(offset+l_max,2) .ne. 0) offset = offset+1
                      m_A = sgn*m_A
                      do l=offset,l_max,2
                         coupling = cpl%rgaunt(l43,l21,l,m43,m21,m_A)
                         cnt = l-l_min+1
                         base = l*l+l+m_A+1
                         do ijkl=1,n_ijkl
                            ind_sum = cnt + total_l*(ijkl-1) !save in order (l,ij,kl)
                            sum_over_m(ind_sum) = coupling*solid_harmonics_abcd(base + n_harmonics*(ijkl-1))
                         enddo
                      enddo

                      m = sgn*m_B
                      if (m .ne. m_A) then !second m-value (if needed)
                         offset = max(l_min,m_B)
                         if (mod(offset+l_max,2) .ne. 0) offset = offset+1
                         do l=offset,l_max,2
                            coupling = cpl%rgaunt(l43,l21,l,m43,m21,m)
                            cnt = l-l_min+1
                            base = l*l+l+m+1
                            do ijkl=1,n_ijkl
                               ind_sum = cnt + total_l*(ijkl-1) !save in order (l,ij,kl)
                               sum_over_m(ind_sum) = sum_over_m(ind_sum) &
                                            + coupling*solid_harmonics_abcd(base + n_harmonics*(ijkl-1))
                            enddo
                         enddo
                      endif

                      base_solid_harmonics = l_min**2 !there are ((l_min-1)+1)**2 combinations of the solid harmonics Y(M,L) with L=0,...,l_min-1
                      terms_to_skip = l_min*n_ijkl !number of cpl_ABCD(ij,kl,L) indices preceeding the block with L=l_min
                      total_l = l_max-l_min+1 !total number of different l-values for each i,j,k,l combination.

                      base_cf = n_ijkl*((2+(l_min-1)*l_min)/2) !number of values in l_ijkl_L_coefficients preceeding the value l_ijkl_L_coefficients(0,1,1,l_min). In parenth.: sum_L=0,...,l_min-1 max(1,L)
                      if (l_min .eq. 0) base_cf = 0

                      dAB = n_ij*(l21m21_rank-1) + n_l21m21*n_ijkl*(l43m43_rank-1)

                      dBA = n_ij*(l43m43_rank-1) + n_l21m21*n_ijkl*(l21m21_rank-1)

                      !Zero-out the terms that are zero by symmetry: terms with big_l=0,1,...,l_min-1
                      do big_l=0,l_min-1
                         stride = lm_total*big_l
                         do kl=1,n_kl
                            do ij=1,n_ij
!                            i = (2*l21m21_rank+2*ij*(1-l43m43_rank+n_l21m21)-(2 + (-3+l43m43_rank)*l43m43_rank+(-1+kl)*(-1+n_l43m43)*n_l43m43)*n_ij+2*n_l21m21*(-1+(-1+l43m43_rank+(-1+kl)*n_l43m43)*n_ij))/2
!                               ind = i + stride
!                               ind = stride + l21m21_rank + n_l21m21*(ij-1) + n_l21m21*dBA_1 + n_l21m21*n_ij*n_l43m43*(kl-1) !(l21m21,ij,l43m43,kl,L)
                               ind = stride + ij + dAB + n_l21m21*n_ij*(kl-1) !(ij,l21m21,kl,l43m43,L)
                               cpl_ABCD(ind) = 0.0_cfp
                               if (l21m21_rank .le. n_l43m43) then
!                                  ind = stride + l43m43_rank + n_l21m21*(ij-1) + n_l21m21*dAB_1 + n_l21m21*n_ij*n_l43m43*(kl-1) !(l43m43,ij,l21m21,kl,L)                               
                                   ind = stride + ij + dBA + n_l21m21*n_ij*(kl-1) !(ij,l43m43,kl,l21m21,L)
                                  cpl_ABCD(ind) = 0.0_cfp
                               endif
                            enddo !ij
                         enddo !kl
                      enddo !big_l
                     
                      !For each (ij,kl,L) we sum the terms sum (l) (sum (m)  <l,m|l21,m21|l43,m43>*Y(m,l,ij,kl))*l_ijkl_L_coefficients(l,ij,kl,L) over all (m,l) combinations
                      base = base_cf !skip the l_ijkl_L_coefficients for L=0,1,...,l_min-1
                      do big_l=l_min,l_max
                         n_l = max(1,big_l)
                         l = big_l-l_min+1
                         lm2 = big_l-2-l_min !relative (within each ijkl group) sequence number of the l-value corresponding to l=L-2
                         base_l = 0
                         stride = lm_total*big_l
                         do kl=1,n_kl
                            do ij=1,n_ij
!                            i = (2*l21m21_rank+2*ij*(1-l43m43_rank+n_l21m21)-(2 + (-3+l43m43_rank)*l43m43_rank+(-1+kl)*(-1+n_l43m43)*n_l43m43)*n_ij+2*n_l21m21*(-1+(-1+l43m43_rank+(-1+kl)*n_l43m43)*n_ij))/2
!                               ind = i + stride !save in order (l21m21,ij,l43m43,kl,L)
!                               ind = stride + l21m21_rank + n_l21m21*(ij-1) + n_l21m21*dBA_1 + n_l21m21*n_ij*n_l43m43*(kl-1) !(l21m21,ij,l43m43,kl,L)
                               ind = stride + ij + dAB + n_l21m21*n_ij*(kl-1) !(ij,l21m21,kl,l43m43,L)
                               !sum over l = l_min,...,big_l-2 in steps of 2
                               cpl_ABCD(ind) = sum(sum_over_m(base_l+1:base_l+lm2+1:2) &
                                                    * l_ijkl_L_coefficients(base+l_min+1:base+big_l-1:2))
                               !add the term l = big_l
                               cpl_ABCD(ind) = cpl_ABCD(ind) + sum_over_m(base_l+l)*l_ijkl_L_coefficients(base+n_l)
                               if (l21m21_rank .le. n_l43m43) then
!                                  i = stride + l43m43_rank + n_l21m21*(ij-1) + n_l21m21*dAB_1 + n_l21m21*n_ij*n_l43m43*(kl-1) !(l43m43,ij,l21m21,kl,L)
                                   i = stride + ij + dBA + n_l21m21*n_ij*(kl-1) !(ij,l43m43,kl,l21m21,L)
                                  cpl_ABCD(i) = cpl_ABCD(ind)
                               endif
                               base = base + n_l !shift the base by n_l
                               base_l = base_l + total_l
!                               write(*,'("cpl_ABCD",4i3,i4,i,e)') l21,m21,l43,m43,big_l,ind,cpl_ABCD(ind)
                            enddo !ij
                         enddo !kl
                      enddo !big_l

                      do big_l=l_max+1,la_p_lb+lc_p_ld
                         l = min(l_max,big_l-2)-l_min
                         base_l = 0
                         stride = lm_total*big_l
                         do kl=1,n_kl
                            do ij=1,n_ij
!                            i = (2*l21m21_rank+2*ij*(1-l43m43_rank+n_l21m21)-(2 + (-3+l43m43_rank)*l43m43_rank+(-1+kl)*(-1+n_l43m43)*n_l43m43)*n_ij+2*n_l21m21*(-1+(-1+l43m43_rank+(-1+kl)*n_l43m43)*n_ij))/2
!                               ind = i + stride !save in order (l21m21,ij,l43m43,kl,L)
!                               ind = stride + l21m21_rank + n_l21m21*(ij-1) + n_l21m21*dBA_1 + n_l21m21*n_ij*n_l43m43*(kl-1) !(l21m21,ij,l43m43,kl,L)
                               ind = stride + ij + dAB + n_l21m21*n_ij*(kl-1) !(ij,l21m21,kl,l43m43,L)
                               !sum over l = l_min,...,min(l_max,big_l-2) in steps of 2
                               cpl_ABCD(ind) = sum(sum_over_m(base_l+1:base_l+l+1:2) &
                                                    * l_ijkl_L_coefficients(base+l_min+1:base+l+l_min+1:2))
                               if (l21m21_rank .le. n_l43m43) then
!                                  i = stride + l43m43_rank + n_l21m21*(ij-1) + n_l21m21*dAB_1 + n_l21m21*n_ij*n_l43m43*(kl-1) !(l43m43,ij,l21m21,kl,L)
                                   i = stride + ij + dBA + n_l21m21*n_ij*(kl-1) !(ij,l43m43,kl,l21m21,L)
                                  cpl_ABCD(i) = cpl_ABCD(ind)
                               endif
                               base = base + big_l !shift the base by big_l
                               base_l = base_l + total_l
!                               write(*,'("cpl_ABCD",4i3,i4,i,e)') l21,m21,l43,m43,big_l,ind,cpl_ABCD(ind)
                            enddo !ij
                         enddo !kl
                      enddo !big_l

                   enddo !m21
                enddo !l21
             enddo !m43
          enddo !l43

   end subroutine calculate_ABCD_coefficients

   subroutine cfp_solh_1d(SH,x,y,z,L)
     implicit none
     integer, intent(in) :: L
     real(kind=cfp), intent(out) :: SH(:)
     real(kind=cfp), intent(in) :: x, y, z
     real(kind=cfp), parameter :: fac_l1 = sqrt(3.0_cfp/4.0_cfp)
   
     integer :: l_it, m_it, lp1, l2p1, base_l_it_p1, base_l_it, base_l_it_m1, ind
     real(kind=cfp) :: rsq, fac, z_cf
   
     SH(1:(L+1)*(L+1)) = 0.0_cfp
    
     !initialize the starting values
     SH(1) = 1.0_cfp !(0,0)

     if (x .eq. 0.0_cfp .and. y .eq. 0.0_cfp .and. z .eq. 0.0_cfp) return
    
     select case (L)
     case (1) !L=1

        SH(2) = y !(-1,1)
        SH(3) = z !(0,1)
        SH(4) = x !(1,1)

     case (2:) !L>1 

        rsq = x*x + y*y + z*z

        !l_it = 0 case
        SH(2) = y !(-1,1)
        SH(3) = z !(0,1)
        SH(4) = x !(1,1)

        !diagonal recursions
        !l_it = 1 case separately since the recursions would attempt to address SH(1,0) and SH(-1,0) which lie out of range.
        l_it = 1
        SH(9) = fac_l1*(x*SH(4)-y*SH(2))
        SH(5) = fac_l1*(y*SH(4)+x*SH(2))
        l2p1 = 3
        lp1 = 2
        !vertical recursions
        SH(6) = (l2p1*z*SH(2))/sqrt(1.0_cfp*(lp1-1)*(lp1+1)) !todo precalculate the l-dependent part
        SH(7) = (l2p1*z*SH(3)-rsq)/sqrt(1.0_cfp*lp1*lp1)
        SH(8) = (l2p1*z*SH(4))/sqrt(1.0_cfp*(lp1+1)*(lp1-1))

        !diagonal recursions (general case l_it .ge. 2)
        ind = 4
        do l_it = 2, L-1
           base_l_it_p1 = (l_it+1)*(l_it+1)+l_it+1 +1
           base_l_it = (l_it)*(l_it)+l_it +1
           base_l_it_m1 = (l_it-1)*(l_it-1)+l_it-1 +1
           fac = SH_l(l_it+1) != sqrt((2.0_cfp*l_it+1.0_cfp)/(2.0_cfp*l_it+2.0_cfp)) !use the precalculated coefficients
           SH(l_it+1+base_l_it_p1) = fac*(x*SH(l_it+base_l_it)-y*SH(-l_it+base_l_it))
           SH(-l_it-1+base_l_it_p1) = fac*(y*SH(l_it+base_l_it)+x*SH(-l_it+base_l_it))
           l2p1 = 2*l_it + 1
           lp1 = l_it + 1
           !vertical recursions (general case)
           ind = ind + l_it+1
           z_cf = l2p1*z
           !todo add compiler directive specifying no dependencies in the loop so it can be vectorized.
           do m_it = -l_it,l_it
!              SH(m_it+base_l_it_p1) = (l2p1*z*SH(m_it+base_l_it)-sqrt(1.0_cfp*(l_it+m_it)*(l_it-m_it))*rsq*SH(m_it+base_l_it_m1))/sqrt(1.0_cfp*(lp1+m_it)*(lp1-m_it))
              SH(m_it+base_l_it_p1) = z_cf*SH(m_it+base_l_it)*SH_ml(1,ind+m_it)-rsq*SH(m_it+base_l_it_m1)*SH_ml(2,ind+m_it) !use the precalculated coefficients instead
           end do
           ind = ind + l_it
        end do

     end select
  
   end subroutine cfp_solh_1d

   function cnla(n,l,alp)
      implicit none
      real(kind=cfp) :: cnla
      integer, intent(in) :: n, l
      real(kind=cfp), intent(in) :: alp

      integer :: i

         if (n .eq. -1) then
            cnla = 1.0_cfp/(2.0_cfp*alp)**(l+0.5_cfp)
            return
         endif

         cnla = 1.0_cfp

         do i=1,n
            cnla = 2.0_cfp*cnla*i
         enddo

         cnla = cnla/(2.0_cfp*alp)**(n+l+1.5_cfp)

   end function cnla

   function Lag_n_hlf_k(n,l,arg)
      use special_functions, only: cfp_eval_poly_horner
      implicit none
      integer, intent(in) :: n, l
      real(kind=cfp), intent(in) :: arg
      real(kind=cfp) :: Lag_n_hlf_k

      integer :: base
 
!      real(kind=cfp) :: cf(1:n+1), f
!      integer :: m
!
!        Lag_n_hlf_k = 0.0_cfp
!
!         f = 1.0_cfp
!         do m=0,n
!           if (m .ne. 0) f = f*m
!           cf(m+1) = (-1)**m/f*gen_binom(n+l+0.5_cfp,n-m)
!         enddo
!
!         Lag_n_hlf_k = cfp_eval_poly_horner(n,arg,cf)

         !Use the precalculated coefficients for the Lag. polynomials
         base = n_l*n*(n+1)/2 + l*(n+1) !base index for the coefficents cf_m(l,n)
         Lag_n_hlf_k = cfp_eval_poly_horner(n,arg,cf(base+1:base+n+1))

   end function Lag_n_hlf_k

   function gen_binom(x,m)
      implicit none
      real(kind=cfp) :: gen_binom
      real(kind=cfp), intent(in) :: x
      integer, intent(in) :: m

         gen_binom = gamma(x+1.0_cfp)/(gamma(m+1.0_cfp)*gamma(x-m+1.0_cfp))

   end function gen_binom

end module eri_sph_coord
