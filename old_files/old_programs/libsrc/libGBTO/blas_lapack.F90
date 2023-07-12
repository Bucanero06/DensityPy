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
!> Interface module for BLAS and LAPACK routine. The interface is for double and quad precision routines. The double precision routines are assumed external and the quad precision ones are here.
module blas_lapack
   use precisn
   use utils, only: xermsg

   private

   ! BLAS and LAPACK integer interface
#if defined(blas64bit)
   integer, parameter, public :: blasint = longint
#elif defined(blas32bit)
   integer, parameter, public :: blasint = shortint
#else
   integer, parameter, public :: blasint = kind(0)
#endif

   !> Interface for LAPACK dsyev routine and its quad precision implementation.
   !> \warning the Quad precision version is not an extension of the LAPACK double routine but an implementation of the NR matrix diagonalization based on the QR algorithm.
   interface syev
      module procedure dble_syev, quad_syev
   end interface

   !> Interface for LAPACK dsyevr routine and its quad precision implementation.
   !> \warning the Quad precision version is not an extension of the LAPACK double routine but an implementation of the NR matrix diagonalization based on the QR algorithm.
   interface syevr
      module procedure dble_syevr, quad_syevr
   end interface

   !> Interface for LAPACK dtrsm routine and its quad precision implementation.
   !> \warning the Quad precision version of the routine has not been implemented yet.
   interface trsm
      module procedure dble_trsm, quad_trsm
   end interface

   public syev, syevr, trsm 

contains

   subroutine dble_trsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
      !.. Scalar Arguments ..
      real(kind=wp) ALPHA
      INTEGER(blasint) LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      !..
      !.. Array Arguments ..
      real(kind=wp) A(:,:), B(:,:)  !A(LDA,*),B(LDB,*)
      !
         !Standard LAPACK double precision
         call dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

   end subroutine dble_trsm

   subroutine quad_trsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
      !.. Scalar Arguments ..
      real(kind=ep1) ALPHA
      INTEGER(blasint) LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      real(kind=ep1) A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      real(kind=ep1) temp
      integer i,info,j,k,nrowa
      logical lside,nounit,upper
!     ..
!     .. Parameters ..
      real(kind=ep1) one,zero
      parameter (one=1.0_ep1,zero=0.0_ep1)
   !     ..
   !
   !     Test the input parameters.
   !
         lside = lsame(side,'L')
         if (lside) then
             nrowa = m
         else
             nrowa = n
         end if
         nounit = lsame(diag,'N')
         upper = lsame(uplo,'U')
   !
         info = 0
         if ((.not.lside) .and. (.not.lsame(side,'R'))) then
             info = 1
         else if ((.not.upper) .and. (.not.lsame(uplo,'L'))) then
             info = 2
         else if ((.not.lsame(transa,'N')) .and. (.not.lsame(transa,'T')) .and. (.not.lsame(transa,'C'))) then
             info = 3
         else if ((.not.lsame(diag,'U')) .and. (.not.lsame(diag,'N'))) then
             info = 4
         else if (m.lt.0) then
             info = 5
         else if (n.lt.0) then
             info = 6
         else if (lda.lt.max(1,nrowa)) then
             info = 9
         else if (ldb.lt.max(1,m)) then
             info = 11
         end if
         if (info.ne.0) then
             !report the lapack error message so we can find out what the problem is.
             call xermsg ('blas_lapack', 'quad_trsm', &
                          'Error in input parameters. See LAPACK DTRSM error message list for details', info, 1)
         end if
   !
   !     Quick return if possible.
   !
         if (m.eq.0 .or. n.eq.0) return
   !
   !     And when  alpha.eq.zero.
   !
         if (alpha.eq.zero) then
             do j = 1,n
                 do i = 1,m
                     b(i,j) = zero
                 enddo
             enddo
             return
         end if
   !
   !     Start the operations.
   !
         if (lside) then
             if (lsame(transa,'N')) then
   !
   !           Form  B := alpha*inv( A )*B.
   !
                 if (upper) then
                     do j = 1,n
                         if (alpha.ne.one) then
                             do i = 1,m
                                 b(i,j) = alpha*b(i,j)
                             enddo
                         end if
                         do k = m,1,-1
                             if (b(k,j).ne.zero) then
                                 if (nounit) b(k,j) = b(k,j)/a(k,k)
                                 do i = 1,k - 1
                                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                                 enddo
                             end if
                         enddo
                     enddo
                 else
                     do j = 1,n
                         if (alpha.ne.one) then
                             do i = 1,m
                                 b(i,j) = alpha*b(i,j)
                             enddo
                         end if
                         do k = 1,m
                             if (b(k,j).ne.zero) then
                                 if (nounit) b(k,j) = b(k,j)/a(k,k)
                                 do i = k + 1,m
                                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                                 enddo
                             end if
                         enddo
                     enddo
                 end if
             else
   !
   !           Form  B := alpha*inv( A**T )*B.
   !
                 if (upper) then
                     do j = 1,n
                         do i = 1,m
                             temp = alpha*b(i,j)
                             do k = 1,i - 1
                                 temp = temp - a(k,i)*b(k,j)
                             enddo
                             if (nounit) temp = temp/a(i,i)
                             b(i,j) = temp
                         enddo
                     enddo
                 else
                     do j = 1,n
                         do i = m,1,-1
                             temp = alpha*b(i,j)
                             do k = i + 1,m
                                 temp = temp - a(k,i)*b(k,j)
                             enddo
                             if (nounit) temp = temp/a(i,i)
                             b(i,j) = temp
                         enddo
                     enddo
                 end if
             end if
         else
             if (lsame(transa,'N')) then
   !
   !           Form  B := alpha*B*inv( A ).
   !
                 if (upper) then
                     do j = 1,n
                         if (alpha.ne.one) then
                             do i = 1,m
                                 b(i,j) = alpha*b(i,j)
                             enddo
                         end if
                         do k = 1,j - 1
                             if (a(k,j).ne.zero) then
                                 do i = 1,m
                                     b(i,j) = b(i,j) - a(k,j)*b(i,k)
                                 enddo
                             end if
                         enddo
                         if (nounit) then
                             temp = one/a(j,j)
                             do i = 1,m
                                 b(i,j) = temp*b(i,j)
                             enddo
                         end if
                     enddo
                 else
                     do j = n,1,-1
                         if (alpha.ne.one) then
                             do i = 1,m
                                 b(i,j) = alpha*b(i,j)
                             enddo
                         end if
                         do k = j + 1,n
                             if (a(k,j).ne.zero) then
                                 do i = 1,m
                                     b(i,j) = b(i,j) - a(k,j)*b(i,k)
                                 enddo
                             end if
                         enddo
                         if (nounit) then
                             temp = one/a(j,j)
                             do i = 1,m
                                 b(i,j) = temp*b(i,j)
                             enddo
                         end if
                     enddo
                 end if
             else
   !
   !           Form  B := alpha*B*inv( A**T ).
   !
                 if (upper) then
                     do k = n,1,-1
                         if (nounit) then
                             temp = one/a(k,k)
                             do i = 1,m
                                 b(i,k) = temp*b(i,k)
                             enddo
                         end if
                         do j = 1,k - 1
                             if (a(j,k).ne.zero) then
                                 temp = a(j,k)
                                 do i = 1,m
                                     b(i,j) = b(i,j) - temp*b(i,k)
                                 enddo
                             end if
                         enddo
                         if (alpha.ne.one) then
                             do i = 1,m
                                 b(i,k) = alpha*b(i,k)
                             enddo
                         end if
                     enddo
                 else
                     do k = 1,n
                         if (nounit) then
                             temp = one/a(k,k)
                             do i = 1,m
                                 b(i,k) = temp*b(i,k)
                             enddo
                         end if
                         do j = k + 1,n
                             if (a(j,k).ne.zero) then
                                 temp = a(j,k)
                                 do i = 1,m
                                     b(i,j) = b(i,j) - temp*b(i,k)
                                 enddo
                             end if
                         enddo
                         if (alpha.ne.one) then
                             do i = 1,m
                                 b(i,k) = alpha*b(i,k)
                             enddo
                         end if
                     enddo
                 end if
             end if
         end if
!
   end subroutine quad_trsm

   subroutine dble_syev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      IMPLICIT NONE
      !.. Scalar Arguments ..
      CHARACTER JOBZ, UPLO
      INTEGER(blasint) INFO, LDA, LWORK, N
      !..
      !.. Array Arguments ..
      real(kind=wp) A(LDA,*), W(*), WORK(*)
      !..
         !Standard LAPACK double precision
         call dsyev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   end subroutine dble_syev

   subroutine quad_syev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      use sort, only: cfp_sort_float_int_1d
      IMPLICIT NONE
      !.. Scalar Arguments ..
      CHARACTER JOBZ, UPLO
      INTEGER(blasint) INFO, LDA, LWORK, N
      !..
      !.. Array Arguments ..
      real(kind=ep1) A(:,:), W(:), WORK(:)
      !..
      INTEGER :: i, err
      INTEGER, ALLOCATABLE :: permutation(:)
      real(kind=ep1), ALLOCATABLE :: z(:,:)
      logical :: novectors
      logical, save :: first = .true.

         if (first) then
            call xermsg ('blas_lapack', 'quad_syev', &
                         'Running matrix diagonalization algorithms tred2, tqli, NOT LAPACK ROUTINES. &
                         &This warning message will be supressed for subsequent calls to quad_syev.', 1, 0)
            first = .false.
         endif

         INFO = 0

         if (N .le. 0) call xermsg('blas_lapack','quad_syev','N .le. 0.',1,1)

         if (lwork .eq. -1) then
            work(1) = N
            return 
         endif

         if (lwork .le. 0) call xermsg('blas_lapack','quad_syev','lwork .le. 0.',2,1)

         !if (size(WORK) < N) call xermsg('blas_lapack','quad_syev','size(WORK) < N.',3,1)

         if (JOBZ .eq. 'V' .or. JOBZ .eq. 'v') then
            novectors = .false. !we use this reverse way of marking in order to be compatible with tred2 routine
         else
            novectors = .true.
         endif

         !todo using Z is not necessary - this should be improved on
         allocate(permutation(N),z(N,N),stat=err)
         if (err .ne. 0) call xermsg('blas_lapack','quad_syev','Memory allocation failed.',err,1)

         do i=1,N
            permutation(i) = i
         enddo

         call balanc(A,int(N)) !balancing of the matrix

         call tred2(A,W,WORK,novectors)

         if (.not.(novectors)) then
            z = A
            call tqli(W,WORK,z)

            !sort the eigenvalues from small to large and order the corresponding eigenvectors
            call cfp_sort_float_int_1d(int(N),W,permutation)
            do i=1,N
               A(1:N,i) = Z(1:N,permutation(i))
            enddo
         else
            call tqli(W,WORK)
            !sort the eigenvalues from small to large
            call cfp_sort_float_int_1d(int(N),W,permutation)
         endif

   end subroutine quad_syev

   subroutine dble_syevr (JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, &
                          WORK, LWORK, IWORK, LIWORK, INFO)
      IMPLICIT NONE
      !.. Scalar Arguments ..
      CHARACTER :: JOBZ, RANGE, UPLO
      INTEGER(blasint) :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
      real(kind=wp) :: ABSTOL, VL, VU
      !..
      !.. Array Arguments ..
      INTEGER(blasint) :: ISUPPZ(*), IWORK(*)
      real(kind=wp) :: A(LDA,*), W(*), WORK(*), Z(LDZ,*)
      !..

         !Standard LAPACK double precision
         call dsyevr( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

   end subroutine dble_syevr

   subroutine quad_syevr (JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, &
                          WORK, LWORK, IWORK, LIWORK, INFO)
      use sort, only: cfp_sort_float_int_1d
      IMPLICIT NONE
      !.. Scalar Arguments ..
      CHARACTER :: JOBZ, RANGE, UPLO
      INTEGER(blasint) :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
      REAL(kind=ep1) :: ABSTOL, VL, VU
      !..
      !.. Array Arguments ..
      INTEGER(blasint) :: ISUPPZ(:), IWORK(:)
      REAL(kind=ep1) :: A(:,:), W(:), WORK(:), Z(:,:)
      !..
      INTEGER :: i, err
      INTEGER, ALLOCATABLE :: permutation(:)
      LOGICAL :: novectors
      logical, save :: first = .true.

         if (first) then
            call xermsg ('blas_lapack','quad_syevr', &
                         'Running matrix diagonalization algorithms tred2, tqli, NOT LAPACK ROUTINES. &
                         &This warning message will be supressed for subsequent calls to quad_syevr.', 1, 0)
            first = .false.
         endif

         INFO = 0

         if (N .le. 0) call xermsg('blas_lapack','quad_syevr','N .le. 0.',1,1)

         if (lwork .eq. -1 .or. liwork .eq. -1) then
            work(1) = N
            iwork(1) = 1
            return
         endif

         if (lwork .le. 0 .or. liwork .le. 0) call xermsg('blas_lapack','quad_syevr','lwork .le. 0 .or. liwork .le. 0.',2,1)

         !if (size(WORK) < N) call xermsg('blas_lapack','quad_syevr','size(WORK) < N.',3,1)

         if (JOBZ .eq. 'V' .or. JOBZ .eq. 'v') then
            novectors = .false.
         else
            novectors = .true.
         endif

         allocate(permutation(N),stat=err)
         if (err .ne. 0) call xermsg('blas_lapack','quad_syev','Memory allocation failed.',err,1)

         do i=1,N
            permutation(i) = i
         enddo

         call balanc(A,int(N)) !balancing of the matrix

         call tred2(A,W,WORK,novectors) !reduction to tridiagonal form

         if (.not.(novectors)) then

            !If eigenvectors are required then A must be initialized to the matrix output by tred2.
            call tqli(W,WORK,A) !eigensystem for a tridiagonal matrix; eigenvectors will be in A

            !sort the eigenvalues from small to large and order the corresponding eigenvectors
            call cfp_sort_float_int_1d(int(N),W,permutation)
            do i=1,N
               Z(1:N,i) = A(1:N,permutation(i))
            enddo
         else
            call tqli(W,WORK)
            !sort the eigenvalues from small to large
            call cfp_sort_float_int_1d(int(N),W,permutation)
         endif

   end subroutine quad_syevr

   subroutine tqli(d,e,z)
      IMPLICIT NONE
      REAL(kind=ep1), INTENT(INOUT) :: d(:),e(:)
      REAL(kind=ep1), OPTIONAL, INTENT(INOUT) :: z(:,:)
      INTEGER :: i,iter,l,m,n,ndum
      REAL(kind=ep1) :: b,c,dd,f,g,p,r,s
      REAL(kind=ep1), DIMENSION(size(e)) :: ff

         if (size(e) .ge. size(d)) then
            n = size(d)
         else
            call xermsg('blas_lapack','tqli','Size of the input array e must be at least size(d).',1,1)
         endif

         if (present(z)) then 
            if (size(z,1) .eq. size(z,2)) then
               ndum = size(z,1)
            else
               call xermsg('blas_lapack','tqli','Dimensions of the array z must be the same.',2,1)
            endif
         endif

         e(1:n)=eoshift(e(1:n),1)
         do l=1,n
            iter=0
            iterate: do
               do m=l,n-1
                  dd=abs(d(m))+abs(d(m+1))
                  if (abs(e(m))+dd == dd) exit
               end do
               if (m == l) exit iterate
               if (iter == 30) call xermsg('blas_lapack','tqli','Too many iterations.',3,1)
               iter=iter+1
               g=(d(l+1)-d(l))/(2.0_ep1*e(l))
               r=pythag(g,1.0_ep1)
               g=d(m)-d(l)+e(l)/(g+sign(r,g))
               s=1.0_ep1
               c=1.0_ep1
               p=0.0_ep1
               do i=m-1,l,-1
                  f=s*e(i)
                  b=c*e(i)
                  r=pythag(f,g)
                  e(i+1)=r
                  if (r == 0.0_ep1) then
                     d(i+1)=d(i+1)-p
                     e(m)=0.0_ep1
                     cycle iterate
                  end if
                  s=f/r
                  c=g/r
                  g=d(i+1)-p
                  r=(d(i)-g)*s+2.0_ep1*c*b
                  p=s*r
                  d(i+1)=g+p
                  g=c*r-b
                  if (present(z)) then
                     ff(1:n)=z(1:n,i+1)
                     z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                     z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
                  end if
               end do
               d(l)=d(l)-p
               e(l)=g
               e(m)=0.0_ep1
            end do iterate
         end do

   end subroutine tqli

   subroutine tred2(a,d,e,novectors)
      IMPLICIT NONE
      REAL(kind=ep1), DIMENSION(:,:), INTENT(INOUT) :: a
      REAL(kind=ep1), DIMENSION(:), INTENT(OUT) :: d,e
      LOGICAL, OPTIONAL, INTENT(IN) :: novectors
      INTEGER :: i,j,l,n
      REAL(kind=ep1) :: f,g,h,hh,scale
      REAL(kind=ep1), DIMENSION(size(a,1)) :: gg
      LOGICAL, SAVE :: yesvec = .true.

         if (size(a,1) .eq. size(a,2) .and. size(a,1) .eq. size(d) .and. size(e) .ge. size(a,1)) then
            n = size(a,1)
         else
            call xermsg('blas_lapack','tred2','Dimensions of the input arrays a,d must be equal and size(e) .ge. size(a,1).',1,1)
         endif

         if (present(novectors)) yesvec = .not.(novectors)
         do i=n,2,-1
            l=i-1
            h=0.0_ep1
            if (l > 1) then
               scale=sum(abs(a(i,1:l)))
               if (scale == 0.0_ep1) then
                  e(i)=a(i,l)
               else
                  a(i,1:l)=a(i,1:l)/scale
                  h=sum(a(i,1:l)**2)
                  f=a(i,l)
                  g=-sign(sqrt(h),f)
                  e(i)=scale*g
                  h=h-f*g
                  a(i,l)=f-g
                  if (yesvec) a(1:l,i)=a(i,1:l)/h
                  do j=1,l
                     e(j)=(dot_product(a(j,1:j),a(i,1:j)) + dot_product(a(j+1:l,j),a(i,j+1:l)))/h
                  end do
                  f=dot_product(e(1:l),a(i,1:l))
                  hh=f/(h+h)
                  e(1:l)=e(1:l)-hh*a(i,1:l)
                  do j=1,l
                     a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
                  end do
               end if
            else
               e(i)=a(i,l)
            end if
            d(i)=h
         end do
         if (yesvec) d(1)=0.0_ep1
         e(1)=0.0_ep1
         do i=1,n
            if (yesvec) then
               l=i-1
               if (d(i) /= 0.0_ep1) then
                  gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
                  a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
               end if
               d(i)=a(i,i)
               a(i,i)=1.0_ep1
               a(i,1:l)=0.0_ep1
               a(1:l,i)=0.0_ep1
            else
               d(i)=a(i,i)
            end if
         end do

   end subroutine tred2

   function outerprod(a,b)
      IMPLICIT NONE
      REAL(kind=ep1), DIMENSION(:), INTENT(IN) :: a,b
      REAL(kind=ep1), DIMENSION(size(a),size(b)) :: outerprod
         outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
   end function outerprod

   function pythag(a,b)
      IMPLICIT NONE
      REAL(kind=ep1), INTENT(IN) :: a,b
      REAL(kind=ep1) :: pythag
      REAL(kind=ep1) :: absa,absb
         absa=abs(a)
         absb=abs(b)
         if (absa > absb) then
            pythag=absa*sqrt(1.0_ep1+(absb/absa)**2)
         else
            if (absb == 0.0_ep1) then
               pythag=0.0_ep1
            else
               pythag=absb*sqrt(1.0_ep1+(absa/absb)**2)
            end if
         end if
   end function pythag

   subroutine balanc(a,n)
      IMPLICIT NONE
      INTEGER :: n, i, j
      REAL(kind=ep1) :: a(:,:)
      REAL(kind=ep1), PARAMETER :: RDX=radix(1.0_ep1), SQRDX=RDX**2
      REAL(kind=ep1) :: last
      REAL(kind=ep1) :: c,f,g,r,s

   1     continue
         last=1.0_ep1
         do i=1,n
            c=0.0_ep1
            r=0.0_ep1
            do j=1,n
               if (j.ne.i) then
                  c=c+abs(a(j,i))
                  r=r+abs(a(i,j))
               endif
            enddo
            if (c.ne.0.0_ep1 .and. r.ne.0.0_ep1) then
               g=r/rdx
               f=1.0_ep1
               s=c+r
   2           if (c.lt.g) then
                  f=f*rdx
                  c=c*SQRDX
                  goto 2
               endif
               g=r*rdx
   3           if (c.gt.g) then
                  f=f/rdx
                  c=c/SQRDX
                  goto 3
               endif
               if ((c+r)/f.lt.0.95_ep1*s) then
                  last=0.0_ep1
                  g=1.0_ep1/f
                  do j=1,n
                     a(i,j)=a(i,j)*g
                  enddo
                  do j=1,n
                     a(j,i)=a(j,i)*f
                  enddo
               endif
            endif
         enddo

         if (last.eq.0.0_ep1) goto 1

   end subroutine balanc

   function lsame(ca,cb)
!
!  -- reference blas level1 routine (version 3.1) --
!  -- reference blas is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!     november 2011
!
!
      logical lsame
!     .. scalar arguments ..
      character ca,cb
!     ..
!
! =====================================================================
!
!     .. intrinsic functions ..
!      intrinsic ichar
!     ..
!     .. local scalars ..
      integer inta,intb,zcode
!     ..
!
!     test if the characters are equal
!
      lsame = ca .eq. cb
      if (lsame) return
!
!     now test for equivalence if both characters are alphabetic.
!
      zcode = ichar('z')
!
!     use 'z' rather than 'a' so that ascii can be detected on prime
!     machines, on which ichar returns a value with bit 8 set.
!     ichar('a') on prime machines returns 193 which is the same as
!     ichar('a') on an ebcdic machine.
!
      inta = ichar(ca)
      intb = ichar(cb)
!
      if (zcode.eq.90 .or. zcode.eq.122) then
!
!        ascii is assumed - zcode is the ascii code of either lower or
!        upper case 'z'.
!
          if (inta.ge.97 .and. inta.le.122) inta = inta - 32
          if (intb.ge.97 .and. intb.le.122) intb = intb - 32
!
      else if (zcode.eq.233 .or. zcode.eq.169) then
!
!        ebcdic is assumed - zcode is the ebcdic code of either lower or
!        upper case 'z'.
!
          if (inta.ge.129 .and. inta.le.137 .or. inta.ge.145 .and. inta.le.153 .or. inta.ge.162 .and. inta.le.169) inta = inta + 64
          if (intb.ge.129 .and. intb.le.137 .or. intb.ge.145 .and. intb.le.153 .or. intb.ge.162 .and. intb.le.169) intb = intb + 64
!
      else if (zcode.eq.218 .or. zcode.eq.250) then
!
!        ascii is assumed, on prime machines - zcode is the ascii code
!        plus 128 of either lower or upper case 'z'.
!
          if (inta.ge.225 .and. inta.le.250) inta = inta - 32
          if (intb.ge.225 .and. intb.le.250) intb = intb - 32
      end if
      lsame = inta .eq. intb
!
   end function lsame

end module blas_lapack
