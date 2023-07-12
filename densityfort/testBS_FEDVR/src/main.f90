! {{{ Detailed description

!> \mainpage Program <ProgramName> <Insert here what the program does>
!! 
!! Synopsis:
!! ---------
!!
!!     <Program Name> <mandatory run-time parameters (RTP)> [<optional RTP>]
!!
!! ___
!! Description:
!! ------------
!!
!! Input parameters:      {#Input_Parameters}
!! =================
!! [...Input](@ref ...Input) as specified in the command line.
!!
!> \file
!!
!!
! }}}
program ProgramTemplate

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO
  use ModuleConstants

  use ModuleMainInterface

  use ModuleBSpline
  use ModuleDiagonalize

  implicit none

  !.. Run-time parameters
  !..
  integer                       :: nSize
  character(len=:), allocatable :: FileName
  
  !.. Local parameters
  real(kind(1d0)), allocatable :: dMat(:,:)
  real(kind(1d0)), allocatable :: dEval(:)

  integer, parameter :: BS_ORDER = 7
  integer, parameter :: BS_NINT1 = 10
  integer, parameter :: BS_MULTB = 1
  integer, parameter :: BS_NINT2 = 20
  integer, parameter :: BS_NNODS = BS_NINT1 + BS_MULTB + BS_NINT2
  integer, parameter :: NBND     = BS_ORDER - BS_MULTB
  integer, parameter :: NINT     = BS_NINT1 + BS_MULTB - 2
  integer, parameter :: NEXT     = BS_NINT2 + BS_MULTB - 2
  real(kind(1d0))    :: BS_grid(BS_NNODS), x

  real(kind(1d0)), parameter :: BS_X0 =  0.d0 
  real(kind(1d0)), parameter :: BS_XB = 10.d0 
  real(kind(1d0)), parameter :: BS_XF = 30.d0 

  integer        , parameter :: NPLOT = 501

  type(ClassBSpline) :: BSet
  
  integer                        :: i,j,k, Bs, NBS
  real(kind(1d0))                :: normf, Efact
  real(kind(1d0))  , allocatable :: Smat(:,:)
  real(kind(1d0))  , allocatable :: Rim(:,:), Sim(:,:), Tim(:,:), Pim(:,:)
  real(kind(1d0))  , allocatable :: Sbm(:,:), Tbm(:,:)
  real(kind(1d0))  , allocatable :: Rem(:,:), Sem(:,:), Tem(:,:), Pem(:,:)
  real(kind(1d0))  , allocatable :: eSi(:), eRi(:), eSe(:), eRe(:), eSb(:)
  real(kind(1d0))  , allocatable :: Bv(:,:), H0m(:,:), H1m(:,:), eH(:)
  procedure(D2DFun), pointer     :: fptr
  
  BS_grid(1:BS_NINT1) = [ ( BS_X0 + (BS_XB-BS_X0) / dble(BS_NINT1) * dble( i ), i = 0, BS_NINT1-1 ) ]
  BS_grid(BS_NINT1+1:BS_NINT1+BS_MULTB) = BS_XB
  BS_grid(BS_NINT1+BS_MULTB+1:) = [ ( BS_XB + (BS_XF-BS_XB) / dble(BS_NINT2) * dble( i ), i = 1, BS_NINT2 ) ]
  
  call BSet%Init(BS_NNODS,BS_ORDER,BS_grid) 
  
  NBS = BSet%GetNBsplines()
  if(NBS/=NINT+NEXT+NBND+2)then
     write(*,*) "SIZE INCONSISTENCY"
     STOP
  endif
  do i = 1, NPLOT
     x = BS_X0 + ( BS_XF - BS_X0 ) / dble( NPLOT - 1 ) * dble( i - 1 )
     write(10,"(i0,*(x,e24.16))") i, x, (BSet%Eval(x,Bs),Bs=1,BSet%GetNBsplines())
  enddo

  !.. Computes overlap matrix over the whole Bspline basis, except the first and last functions
  allocate(Smat(NBS-2,NBS-2))
  Smat=0.d0
  k = BS_ORDER
  fptr=>fun1
  do i = 1, NBS-2
     do j = max(1,i-k+1),min(NBS-2,i+k-1)
        Smat(i,j) = BSet%Integral(fptr,i+1,j+1)
     enddo
  enddo

  !.. Internal metric and position
  allocate(Sim,source=Smat(1:nint,1:nint))
  allocate(eSi(nint))
  call short_Diag( nint, Sim, eSi )
  do i =1, nint
     write(*,*) "eS",i,eSi(i)
     Sim(:,i) = Sim(:,i)/sqrt(eSi(i))
  enddo
  allocate(Rim(nint,nint))
  Rim=0.d0
  fptr=>funR
  do i = 1, nint
     do j = max(1, i-k+1), min( nint,i+k-1)
        Rim(i,j) = BSet%Integral(fptr,i+1,j+1)
     enddo
  enddo
  allocate(Tim(nint,nint))
  Tim = matmul(Rim,Sim)
  Rim = matmul(transpose(Sim),Tim)
  allocate(eRi(nint))
  call short_Diag( nint, Rim, eRi )
  Tim = matmul(Sim,Rim)
  allocate(Pim,source=matmul(Tim,transpose(Tim)))
  do i = 1, nint
     if(BSet%Eval(eRi(i),Tim(:,i),Bsmin=2,Bsmax=nint+1)<0)Tim(:,i)=-Tim(:,i)
  enddo

  !.. External metric and position
  allocate(Sem,source=Smat(nint+nbnd+1:,nint+nbnd+1:))
  allocate(eSe(next))
  call short_Diag( next, Sem, eSe )
  do i =1, next
     write(*,*) "eS",i,eSe(i)
     Sem(:,i) = Sem(:,i)/sqrt(eSe(i))
  enddo
  allocate(Rem(next,next))
  Rem=0.d0
  fptr=>funR
  do i = 1, next
     do j = max(1, i-k+1), min(next,i+k-1)
        Rem(i,j) = BSet%Integral(fptr,i+nint+nbnd+1,j+nint+nbnd+1)
     enddo
  enddo
  allocate(Tem(next,next))
  Tem = matmul(Rem,Sem)
  Rem = matmul(transpose(Sem),Tem)
  allocate(eRe(next))
  call short_Diag( next, Rem, eRe )
  Tem = matmul(Sem,Rem)
  allocate(Pem,source=matmul(Tem,transpose(Tem)))
  do i = 1, next
     if(BSet%Eval(eRe(i),Tem(:,i),Bsmin=nint+nbnd+2,Bsmax=nint+nbnd+next+1)<0)Tem(:,i)=-Tem(:,i)
  enddo

!!$  do i = 1, nint
!!$     write(*,"(i0,*(x,e11.3))") i, eR(i), (BSet%Eval(eR(i),Tim(:,j),Bsmin=2,Bsmax=nint+1),j=1,nint)
!!$  enddo

  !.. Computes the boundary Bsplines
  allocate(Bv(nint+next+nbnd,nbnd))
  Bv=0.d0
  Bv(1:nint,:)=-matmul(Pim,Smat(1:nint,nint+1:nint+nbnd))
  do i=1,nbnd
     Bv(nint+i,i)=1.d0
  enddo
  Bv(nint+nbnd+1:,:)=-matmul(Pem,Smat(nint+nbnd+1:,nint+1:nint+nbnd))
  allocate(Sbm,source=matmul(transpose(Bv),matmul(Smat,Bv)))
  allocate(eSb(nbnd))
  call short_Diag(nbnd,Sbm,eSb)
  do i=1,nbnd
     Sbm(:,i)=Sbm(:,i)/sqrt(eSb(i))
  enddo
  allocate(Tbm,source=matmul(Bv,Sbm))
  
  do i = 1, NPLOT
     x = BS_X0 + ( BS_XF - BS_X0 ) / dble( NPLOT - 1 ) * dble( i - 1 )
     write(20,"(i0,*(x,e24.16))") i, x, &
          (BSet%Eval(x,Tim(:,j),Bsmin=1+1          ,Bsmax=1+nint          ),j=1,nint), &
          (BSet%Eval(x,Tbm(:,j),Bsmin=1+1          ,Bsmax=1+nint+nbnd+next),j=1,nbnd), &
          (BSet%Eval(x,Tem(:,j),Bsmin=1+nint+nbnd+1,Bsmax=1+nint+nbnd+next),j=1,next)          
  enddo

!!$  allocate(H0m(2*nint+1,2*nint+1))
!!$  H0m=0.d0
!!$  fptr=>fun1
!!$  do i=2,2*nint+2
!!$     do j=max(2,i-k+1),min(2*nint+2,i+k-1)
!!$        H0m(i-1,j-1)=BSet%Integral(fptr,i,j,1,1)*0.5d0
!!$     enddo
!!$  enddo
!!$  
!!$  !.. Hamiltonian of the particle in the box
!!$  allocate(H1m(2*nint+1,2*nint+1))
!!$
!!$  H1m(:,1:nint )=matmul(H0m(:,1:nint ),Tim)
!!$  H1m(:,1+nint )=matmul(H0m           ,Bv )
!!$  H1m(:,nint+2:)=matmul(H0m(:,nint+2:),Tim)
!!$
!!$  H0m(1:nint ,:)=matmul(transpose(Tim),H1m(1:nint ,:))
!!$  H0m(1+nint ,:)=matmul(Bv            ,H1m           )
!!$  H0m(nint+2:,:)=matmul(transpose(Tim),H1m(nint+2:,:))
!!$
!!$  Efact=0.5d0*(4.d0*atan(1.d0)/(BS_XF-BS_X0))**2
!!$  allocate(eH(2*nint+1))
!!$  call Short_Diag( 2*nint+1, H0m, eH )
!!$  H1m(1:nint ,:)=matmul(Tim,H0m(1:nint ,:))
!!$  H1m(nint+2:,:)=matmul(Tim,H0m(nint+2:,:))
!!$  do i=1,2*nint+1
!!$     write(*,"(a,x,i0,*(x,e14.6))") "eH",i,eH(i),Efact*i**2, eH(i)/(Efact*i**2)-1.d0
!!$     H1m(:,i)=H1m(:,i)+Bv*H0m(nint+1,i)
!!$  enddo
!!$  
!!$  do i = 1, NPLOT
!!$     x = BS_X0 + ( BS_XF - BS_X0 ) / dble( NPLOT - 1 ) * dble( i - 1 )
!!$     write(30,"(i0,*(x,e24.16))") i, x, (BSet%Eval(x,H1m(:,j),Bsmin=2,Bsmax=1+2*nint+1),j=1,2*nint+1)
!!$  enddo
  
  
  
  
!!$  allocate(dMat(nSize,nSize))
!!$  call random_number(dMat)
!!$  dMat = dMat + transpose(dMat)
!!$
!!$  allocate(dEval(nSize))
!!$  call Short_Diag( nSize, dMat, dEval )
!!$
!!$  call SaveVector( FileName, dEval, "formatted" )


contains
  
  Pure DoublePrecision function fun1( x, parvec ) result(res)
    DoublePrecision          , intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    res = 1.d0
  end function Fun1

  Pure DoublePrecision function funR( x, parvec ) result(res)
    DoublePrecision          , intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    res = x
  end function FunR

end program ProgramTemplate

