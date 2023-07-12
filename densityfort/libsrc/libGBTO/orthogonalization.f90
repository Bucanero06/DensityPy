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
module orthogonalization
   use precisn

   private

   public GS_ortho_routine, SYM_ortho_routine

contains

   !> Performs orthogonalization of several vectors (vecs(1:no_cf,a_s:a_e)) against the set of fixed vectors (vecs(1:no_cf,p_s:p_e))
   !> which are assumed to be orthogonal with respect to the overlap matrix ao_overlaps. The number of AO coefficients is no_cf 
   !> and the range of active vectors (columns of vecs) is given by p_s,p_e. The active and passive sets of vectors are assumed to be disjunct.
   !> mo2so_range gives the index of the first and the last AO which contributes to a given MO (given in the column index).
   !> Note that the use of this array ensures that symmetry adaptaion of the AO basis is handled with maximum efficiency.
   !> Note that we cannot use mo2so_range in the loops involving the orbitals to orthogonalize since coefficients of these orbitals
   !> change. However, these loops can be made more efficient if I supply the routine with the information on the range of the symmetry
   !> adapted basis set in each symmetry (once symmetry adaptation is implemented).
   subroutine GS_ortho_routine(no_cf,p_s,p_e,a_s,a_e,vecs,ao_overlaps,symmetry,mo2so_range,del_thrs)
      use const, only: thrs_gs_ortho, thrs_lin_dep_gs_ortho
      use utils, only: xermsg
      implicit none
      integer, intent(in) :: no_cf, p_s,p_e,a_s,a_e, mo2so_range(:,:), symmetry
      real(kind=cfp), intent(in) :: ao_overlaps(:,:), del_thrs(no_cf)
      real(kind=cfp), intent(inout) :: vecs(:,:)

      integer :: i,j,k,l,m
      real(kind=cfp) :: proj, norm

         do i=a_s,a_e

            do m=p_s,p_e !remove from the vector projections on all fixed (orthogonal) vectors
               proj = 0.0_cfp
               do k=1,no_cf
                  if (vecs(k,i) .ne. 0.0_cfp) then
                     do l=mo2so_range(1,m),mo2so_range(2,m)
                        proj = proj + vecs(k,i)*vecs(l,m)*ao_overlaps(l,k)
                     enddo
                  endif
               enddo
               !subtract the m-th projection on the orbital
               if (proj .ne. 0.0_cfp) then
                  do k=mo2so_range(1,m),mo2so_range(2,m)
                     vecs(k,i) = vecs(k,i) - proj*vecs(k,m)
                  enddo
               endif
            enddo !m

            !calculate norm of the i-th orbital: R(i,i)
            norm = 0.0_cfp
            do j=1,no_cf
               if (vecs(j,i) .ne. 0.0_cfp) then
                  norm = norm + vecs(j,i)*vecs(j,i)*ao_overlaps(j,j)
                  do k=j+1,no_cf
                     norm = norm + 2.0_cfp*vecs(j,i)*vecs(k,i)*ao_overlaps(j,k)
                  enddo
               endif
            enddo
      
            if (norm < thrs_lin_dep_gs_ortho) then
               print *,i,norm,thrs_lin_dep_gs_ortho
               call xermsg ('orthogonalization', 'GS_ortho_routine', &
                            'A linear dependency detected during the G-S orthogonalization.', 1, 1)
            endif
      
            norm = 1.0_cfp/sqrt(norm)
        
            do j=1,no_cf 
               vecs(j,i) = vecs(j,i)*norm !normalize the orbital: Q(:,i) = A(:,i)/R(i,i)
               if (abs(vecs(j,i)) < del_thrs(j)) vecs(j,i) = 0.0_cfp !delete coefficients smaller than a threshold
            enddo
      
            do j=i+1,a_e !remove one-by-one components of the i-th vector from the vectors i+1,...,a_e
               proj = 0.0_cfp !projection of the i-th orbital on the j-th orbital
               do k=1,no_cf
                  if (vecs(k,j) .ne. 0.0_cfp) then
                     do l=1,no_cf
                        proj = proj + vecs(l,i)*vecs(k,j)*ao_overlaps(l,k) !R(i,j) = Q(:,i)**T * A(:,j)
                     enddo
                  endif
               enddo
               !subtract the i-th projection on the j-th orbital
               if (proj .ne. 0.0_cfp) then
                  do k=1,no_cf
                     vecs(k,j) = vecs(k,j) - proj*vecs(k,i) !A(:,j) = A(:,j) - Q(:,i)*R(i,j)
                  enddo
               endif
      
            enddo !j
      
         enddo !i

   end subroutine GS_ortho_routine

   !> Performs symmetric orthogonalization of orbitals given in columns of vecs(:,:). The range of orbitals (columns)
   !> to orthogonalize is given by a_s, a_e. The overlap matrix in the basis of the orbitals
   !> being orthogonalized is given in mo_overlaps. The number of coefficients (rows of vecs) each orbital has is no_cf.
   !> The deletion threshold is given by del_thr. The logical array marking the orbitals
   !> for deletion is to_delete.
   subroutine SYM_ortho_routine(no_cf,a_s,a_e,vecs,mo_overlaps,del_thr,del_thr_ao,to_delete)
      use const, only: stdout
      use utils, only: xermsg
      use blas_lapack, only: blasint, syevr, syev
      implicit none
      integer, intent(in) :: no_cf, a_s,a_e
      real(kind=cfp), intent(in) :: mo_overlaps(:,:), del_thr, del_thr_ao(no_cf)
      real(kind=cfp), intent(inout) :: vecs(:,:)
      logical, intent(out) :: to_delete(:)

      real(kind=cfp) :: norm
      real(kind=cfp) :: cfp_dummy
      integer :: i,j

      !variables needed for syevr
      integer(blasint) :: n,lda,il,iu,m,ldz,liwork,lwork,info
      real(kind=cfp), allocatable :: work(:), evals(:), eigvec(:,:), z(:,:)
      integer(blasint), allocatable :: iwork(:), isuppz(:)
      real(kind=cfp) :: vl, vu, abstol

         n = a_e-a_s+1

         !workspace query
         allocate(work(1:1),iwork(1:1),evals(n),eigvec(n,n),z(n,n),isuppz(2*n),stat=info)
         if (info .ne. 0) call xermsg('orthogonalization','SYM_ortho_routine','Memory allocation 1 failed.',int(info),1)

         lda = n
         abstol = f1mach(1,cfp_dummy)
         ldz = n
         lwork = -1
         liwork = -1
         call syevr('V','A','U',n,eigvec,lda,vl,vu,il,iu,abstol,m,evals,z,ldz,isuppz,work,lwork,iwork,liwork,info)
         if (info .ne. 0) call xermsg('orthogonalization','SYM_ortho_routine','dsyevr: workspace query failed.',int(info),1)

         lwork = work(1)
         liwork = iwork(1)
         deallocate(work,iwork)

         allocate(work(1:lwork),iwork(1:liwork),stat=info)
         if (info .ne. 0) call xermsg('orthogonalization','SYM_ortho_routine','Memory allocation 2 failed.',int(info),1)

         !diagonalization
         eigvec(1:n,1:n) = -mo_overlaps(1:n,1:n)
         call syevr('V','A','U',n,eigvec,lda,vl,vu,il,iu,abstol,m,evals,z,ldz,isuppz,work,lwork,iwork,liwork,info)
         if (info /= 0) then
            call xermsg ('orthogonalization', 'SYM_ortho_routine', &
                         'dsyevr: diagonalization of the overlap matrix failed.', int(info), 1)
         end if
         evals(1:n) = -evals(1:n)
         eigvec = z

         write(stdout,'(/,"Eigenvalues of the orbital overlap matrix:",/)')

         to_delete(:) = .false.

         do i=1,n !count how many functions we delete and normalize the eigenvectors

            if (evals(i) > 0.0_cfp) then !make sure that the negative eigenvalues don't cause problems
               eigvec(1:n,i) = eigvec(1:n,i)/sqrt(evals(i)) !normalization of the i-th eigenvector
            endif

            do j=1,n !get rid of the coefficients that are too small - this helps to stabilize the integral transformation
               if (abs(eigvec(j,i)) < del_thr_ao(j)) eigvec(j,i) = 0.0_cfp
            enddo

            if (evals(i) .le. del_thr) then !deletion threshold check
               to_delete(a_s+i-1) = .true.
               write(stdout,'(i0,1x,e25.15," D")') i,evals(i)
            else
               write(stdout,'(i0,1x,e25.15)') i,evals(i)
            endif

         enddo

         !finally, calculate the orthogonalized continuum functions using the transformation matrix given by eigenvectors of the overlap matrix
         vecs(1:no_cf,a_s:a_e) = matmul(vecs(1:no_cf,a_s:a_e),eigvec(1:n,1:n))

   end subroutine SYM_ortho_routine

end module orthogonalization
