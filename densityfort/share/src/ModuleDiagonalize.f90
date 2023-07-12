!! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
module ModuleDiagonalize

  use, intrinsic :: ISO_FORTRAN_ENV 
  
  use ModuleErrorHandling

  implicit none

  private

  real(kind(1d0)), parameter :: AIMAG_THRESHOLD = 1.d-14
  logical :: MODULE_DIAG_DEBUG_FLAG = .FALSE.

  private :: Short_DSYEV 
  private :: Short_ZGEEV 

  interface Short_Diag
     module procedure Short_DSYEV, Short_ZGEEV, Short_ZGEEV_LR
  end interface Short_Diag
  public :: Short_ZHEEV
  public :: Short_Diag
  public :: ModuleDiagonalize_SetDebug  
  public :: ModuleDiagonalize_UnsetDebug  

contains
  
  subroutine ModuleDiagonalize_SetDebug()
      MODULE_DIAG_DEBUG_FLAG = .TRUE.
  end subroutine ModuleDiagonalize_SetDebug

  subroutine ModuleDiagonalize_UnsetDebug()
      MODULE_DIAG_DEBUG_FLAG = .FALSE.
  end subroutine ModuleDiagonalize_UnsetDebug

  subroutine Short_DSYEV( N, A, E )
    !
    implicit none
    !
    integer        , intent(in)    :: N
    real(kind(1d0)), intent(inout) :: A(:,:)
    real(kind(1d0)), intent(out)   :: E(:)
    !
    real(kind(1d0)), allocatable :: C(:,:)
    integer :: INFO, LWORK
    real(kind(1d0)), allocatable :: WORK(:)
    character(len=8) :: ERRMSG
    !
    if(N<=0)return
    !
    allocate(C(N,N))
    C=A(1:N,1:N)
    !
    !.. Determines optimal dimension for working space
    LWORK=-1 !.. This means that the call to DSYEV is a size query
    allocate(WORK(1))!.. The optimal size is returned in WORK(1)
    call DSYEV("V","U",N,C,N,E,WORK,LWORK,INFO)
    LWORK=int(WORK(1)+0.5d0)
    deallocate(WORK)
    !
    !.. Dimension the workspace to the optimal size
    !   and performs the diagonalization for real
    allocate(WORK(LWORK))
    call DSYEV("V","U",N,C,N,E,WORK,LWORK,INFO)
    deallocate(WORK)
    if( INFO /=0 )then
       write(ERRMSG,"(i0)") INFO
       call ErrorMessage("DSYEV Failed: INFO ="//trim(ERRMSG))
       stop
    endif
    A(1:N,1:N)=C
    deallocate(C)
    !
    return
    !
  end subroutine Short_DSYEV
  

  subroutine Short_ZGEEV(n,A,E)
    
    !.. Assumes it is complex and symmetric. Not clear how useful this is

    implicit none

    integer, intent(in) :: n
    complex(kind(1d0)), intent(inout) :: A(:,:), E(:)

    real   (kind(1d0)), allocatable :: DA(:,:), DE(:)
    complex(kind(1d0))              :: DummyZMat(1,1)
    complex(kind(1d0)), allocatable :: zwork(:), C(:,:)
    real   (kind(1d0)), allocatable :: rwork(:)
    integer           , allocatable :: iperm(:)
    integer                         :: lwork, info
    complex(kind(1d0))              :: zw1
    real   (kind(1d0))              :: norm_Real_A
    real   (kind(1d0))              :: norm_Imag_A
    character(len=512) :: errmsg
    integer :: i

    allocate(zwork(1))
    allocate(rwork(2*n))
    allocate(C(n,n))
    call ZGEEV('N','V',n,A,n,E,DummyZMat,1,C,n,zwork,-1,rwork,info)
    lwork=int(abs(zwork(1))+1)
    deallocate(zwork)
    allocate(zwork(lwork))
    call ZGEEV('N','V',n,A,n,E,DummyZMat,1,C,n,zwork,lwork,rwork,info)
    if(info/=0)then
       write(errmsg,"(a,i0)") "Diagonalization failed, info=",info
       call ErrorMessage(errmsg)
       stop
    endif
    deallocate(zwork)
    deallocate(rwork)
    
    !
    !.. Rescale the right eigenvectors to obtain new right eigenvectors
    !   $U^R$ which, together with left eigenvectors $U^L = (U^R)^*$, 
    !   obey the standard orthogonality condition $(U^L)^\dagger U^R=1$.
    !..
    do i=1,n
       zw1 = (1.d0,0.d0) / sqrt(sum(C(:,i)**2))
       C(:,i) = zw1 * C(:,i)
       if(real(C(1,i))<0.d0) C(:,i)=-C(:,i)
    enddo
    !
    !.. The left eigenvectors $U^L$ of the matrix H obey the 
    !   following relations
    !
    !   1. 
    !   $$
    !       \sum_j (U^L_jk)^\dagger H_{ji} = \lambda_k (U^L_{ik})^\dagger
    !   $$
    !   
    !   2. They are normalized as
    !   $$
    !      \forall j,\qquad \sum_i ( U^L_{ij} )^2 = 1
    !   $$
    !
    !   3. The right eigenvectors, defined with the convention
    !   $$
    !        \sum_j H_{ij} U^R_{jk} = U^R_{ik} \lambda_k,
    !   $$     
    !   are related to the left eigenvectors $U^L$ by the
    !   following relation
    !   $$
    !       U^L_{ij} = ( U^R_{ij} )^*
    !   $$
    !
    !   4. The following further orthonormality relation also holds
    !   $$
    !       ( U^L )^\dagger U^R = 1
    !   $$

    !.. Order the eigenvalues in ascending order of the real part
    !..
    allocate( RWORK( n ), IPERM( n ) )
    RWORK = dble( E )
    do i = 1, n
       IPERM(i) = i
    enddo
    !
    !.. Return the permutation vector IPERM resulting from sorting
    !   RWORK in increasing order and do not sort RWORK.
    !
    !   The permutation is such that IPERM(i) is the index of the
    !   value in the original order of the RWORK array that is in
    !   the i-th location of the sorted order
    !..
    call DPSORT( RWORK, n, IPERM, 1, INFO )
    deallocate( RWORK )
    if( INFO /=0 )then
       write(errmsg,"(a,i5,a)")" Error ",INFO," in DPSORT"
       call ErrorMessage(errmsg)
       stop
    endif

    allocate(zwork(n))
    zwork=E
    do i=1,n
       E(i)=zwork(iperm(i))
       A(:,i)=C(:,iperm(i))
    enddo
    deallocate(C,zwork)
    
    return
  end subroutine Short_ZGEEV



  subroutine Short_ZHEEV(n,A,E)
    
    !.. Diagonalizes complex hermitean matrix

    implicit none

    integer           , intent(in)    :: n
    complex(kind(1d0)), intent(inout) :: A(:,:)
    real   (kind(1d0)), intent(inout) :: E(:)
    

    complex(kind(1d0)), allocatable :: zwork(:), C(:,:)
    real   (kind(1d0)), allocatable :: rwork(:)
    integer           , allocatable :: iperm(:)
    integer                         :: lwork, info, i
    character(len=512) :: errmsg
    complex(kind(1d0)) :: zw1

    allocate(zwork(1))
    allocate(rwork(3*n))
    call ZHEEV('V','U',n,A,n,E,zwork,  -1 ,rwork,info)
    lwork=int(abs(zwork(1))+1)
    deallocate(zwork)
    allocate(zwork(lwork))
    call ZHEEV('V','U',n,A,n,E,zwork,lwork,rwork,info)
    if(info/=0)then
       write(errmsg,"(a,i0)") "Diagonalization failed, info=",info
       call ErrorMessage(errmsg)
       stop
    endif
    deallocate(zwork)
    deallocate(rwork)
    
    !.. Rescale the right eigenvectors to obtain new right eigenvectors
    !   $U^R$ which, together with left eigenvectors $U^L = (U^R)^*$, 
    !   obey the standard orthogonality condition $(U^L)^\dagger U^R=1$.
    !..
    do i=1,n
       zw1 = (1.d0,0.d0) / sqrt(sum(abs(A(:,i))**2))
       A(:,i) = zw1 * A(:,i)
       if(real(A(1,i))<0.d0) A(:,i)=-A(:,i)
    enddo
    return
  end subroutine Short_ZHEEV


  subroutine Short_ZGEEV_LR(n,A,E,LC,RC)

    implicit none

    integer, intent(in) :: n
    complex(kind(1d0)), intent(inout) :: A(:,:), E(:)
    complex(kind(1d0)), intent(out)   :: LC(:,:), RC(:,:)

    complex(kind(1d0)), allocatable :: zmat(:,:)
    complex(kind(1d0)), allocatable :: zwork(:)
    real   (kind(1d0)), allocatable :: rwork(:)
    integer           , allocatable :: iperm(:)
    integer                         :: lwork, info
    complex(kind(1d0))              :: zw1
    real   (kind(1d0))              :: norm_Real_A
    real   (kind(1d0))              :: norm_Imag_A
    character(len=512) :: errmsg
    integer :: i

    allocate(zwork(1))
    allocate(rwork(2*n))
    LC=(0.d0,0.d0)
    RC=(0.d0,0.d0)
    call ZGEEV('V','V',n,A,n,E,LC,n,RC,n,zwork,-1,rwork,info)
    lwork=int(abs(zwork(1))+1)
    deallocate(zwork)
    allocate(zwork(lwork))
    call ZGEEV('V','V',n,A,n,E,LC,n,RC,n,zwork,lwork,rwork,info)
    if(info/=0)then
       write(errmsg,"(a,i0)") "Diagonalization failed, info=",info
       call ErrorMessage(errmsg)
       stop
    endif
    deallocate(zwork)
    deallocate(rwork)
    
    !
    !.. Rescale the right and left eigenvectors to obtain the standard 
    !   orthogonality condition $(U^L)^\dagger U^R=1$.
    !..
    do i=1,n
       zw1 = (1.d0,0.d0) / sqrt(sum(RC(:,i)**2))
       RC(:,i) = zw1 * RC(:,i)
       if(real(RC(1,i))<0.d0) RC(:,i)=-RC(:,i)
       zw1 = sum(conjg(LC(:,i))*RC(:,i))
       LC(:,i)=LC(:,i)/conjg(zw1)
    enddo
    !
    !.. The left eigenvectors $U^L$ of the matrix H obey the 
    !   following relations
    !
    !   1. 
    !   $$
    !       \sum_j (U^L_jk)^\dagger H_{ji} = \lambda_k (U^L_{ik})^\dagger
    !   $$
    !   
    !   2. They are normalized as
    !   $$
    !      \forall j,\qquad \sum_i ( U^L_{ij} )^2 = 1
    !   $$
    !
    !   3. The right eigenvectors, defined with the convention
    !   $$
    !        \sum_j H_{ij} U^R_{jk} = U^R_{ik} \lambda_k,
    !   $$     
    !   are related to the left eigenvectors $U^L$ by the
    !   following relation
    !   $$
    !       U^L_{ij} = ( U^R_{ij} )^*
    !   $$
    !
    !   4. The following further orthonormality relation also holds
    !   $$
    !       ( U^L )^\dagger U^R = 1
    !   $$

    !.. Order the eigenvalues in ascending order of the real part
    !..
    allocate( RWORK( n ), IPERM( n ) )
    RWORK = dble( E )
    do i = 1, n
       IPERM(i) = i
    enddo
    !
    !.. Return the permutation vector IPERM resulting from sorting
    !   RWORK in increasing order and do not sort RWORK.
    !
    !   The permutation is such that IPERM(i) is the index of the
    !   value in the original order of the RWORK array that is in
    !   the i-th location of the sorted order
    !..
    call DPSORT( RWORK, n, IPERM, 1, INFO )
    deallocate( RWORK )
    if( INFO /=0 )then
       write(errmsg,"(a,i5,a)")" Error ",INFO," in DPSORT"
       call ErrorMessage(errmsg)
       stop
    endif

    allocate(zwork(n),zmat(n,n))
    zwork=E
    do i=1,n
       E(i)=zwork(iperm(i))
    enddo
    zmat=RC
    do i=1,n
       RC(:,i)=zmat(:,iperm(i))
    enddo
    zmat=LC
    do i=1,n
       LC(:,i)=zmat(:,iperm(i))
    enddo
    deallocate(zwork,zmat)
    
    return
  end subroutine Short_ZGEEV_LR


end module ModuleDiagonalize
