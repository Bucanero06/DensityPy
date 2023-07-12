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
  use ModuleDiagonalize
  use ModuleBspline
  use ModuleAngularMomentum
  use ModuleBasisJUAN
  use ModuleParametersJUAN
!  use ModuleDensityMatricesJUAN
  use bspline_grid_mod
  !.. 
  use mpi_mod 
  use ukrmol_interface
  use symmetry
  use precisn
 ! use ModuleUKRmolInterface
  !..

  implicit none

  !.. Run-time parameters
  !...

  integer*8                 :: NFTINT

  integer*8                 :: irrdim
  type(ClassBspline)             :: Bspline,Bsplinex
  
  !.. Local parameters
  integer*8        :: in8

  !.. this refers to ukrmol interface variables
  type(bspline_grid_obj) :: bspline_grid
  integer*8 :: lunit, posit
  integer*8 :: lb, bspline_index, number_of_functions, pos_after_rw
  real(kind=cfp)  :: norm
  logical*8 :: non_zero_at_boundary
  
  integer         :: ia,ib,ic,id,i,j,k,l,m,lp,mp,i1,i2,i3,i4,iflag
  real(kind(1d0)) :: ivalue,ivalue2,ivalue3,dr,r
  

  integer                        :: IOSTAT
  real(kind(1d0))                :: Integral
  character(len=12)              :: GridFile
  character(len=2)               :: PB
  integer                 :: indx
  real(kind(1d0))         :: H0d
  type(BasisElementInfo), allocatable :: ielement(:),iprim(:)

  integer                        ::  LowInBspl
  integer                        ::  LasExBspl
 ! type(ClassDensityMatrices)             :: STEX_DM

!  External Functions
  integer, external :: Kron,irr_lm
  real(kind(1d0)), external   :: GET_INTEGRALS_E,GET_INTEGRALS_Eno,GET_INTEGRALS_Xi
  real(kind(1d0)), external   :: OneBody_STEX,TwoBody_STEX
  real(kind(1d0)), external   :: OneBody_STEX_NO,TwoBody_STEX_NO
  real(kind(1d0)), external   :: OneBody_STEX_DM !,TwoBody_STEX_DM  

! MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT


 call InitFromFiles

!stop
!  call UKRmol_Orbitals%Init()
!  call UKRMol_Init_Integrals()


  ! irr_series  = nTOT
  ! irr_counter = sTOT
  ! irrdim=size(irr_series)
  ! Nukrmol=sum(irr_series)
  ! Do i=1,8
  !    write(*,*) irr_series(i)
  ! End Do   
  ! write(*,*) "Nukrmol",Nukrmol
  ! write(*,*) "ukrmol    ^ ^ ^ ^ "

!pause
  
!Here we define bsplines which correspond to an extension of the ukrmol bspline basis.
!If we want to set the same basis we must consider the following:
!The number  "no_bsplines" in scatci_integrals input corresponds to NNodes+BsplineOrder-2
! "a", the R-matrix radius is our Rmax (bspline_grid_start = 0.0)
!  bspline_order is BsplineOrder. The  bspline_indices(1,l)  and bspline_indices(2,0) are the first and last bspline elements considered. We must remove the first one for l=0 (bspline_indices(1,0)=2) and the first and second from l>0  (bspline_indices(1,l>0)=3). bspline_indices(2,l) is chossen in order to exclude the bsplines which reach the radious "a", corresponding to the last BsplineOrder-1, so it has to be set to bspline_indices(2,l)=no_bsplines-(bspline_order-1) if we want to remove them from the calculation.
!Once we do the calculation with ukrmol we have to set our bsplines in order to exactly match the bsplines up to bspline_indices(2,l). This can be set by making a grid in a larger domain where it mathc the ukrmol internal grid. There is an overlap region defined by the last bspline_order-1 orbitals considered in ukrmol.
!The bsplines which are not in the internal basis and do overlap with the last elements of the inner basis, the elment number no_bsplines-bspline_order-1, are the ones which go from no_bsplines-bspline_order to no_bsplines-2. Since this last element of the inner basis is connected to  all basis elements due to the orthonormalization processes, these lastmentioned elements are coupled with all the inner ukrmol elements.

  Rmax=a  !our new redius. If we do not want the external bsplines we set Rmax=a, ELSE we define a larger radius:
!  Rmax=30.1333333333333d0 
call InitBsplines(Bspline,Bsplinex)

call InitBasisParameters
allocate(iprim(1:NprimUkrmol))
call InitIprim(iprim)
!Here we organize the one particle basis elements
Nbasis=sum(irr_total)*2 !irr_total(1)!sum(irr_total)*2!irr_total(1)!irr_total(1)*2      ! sum(irr_total)*2   !        !Spin
write(*,*) "New dimension of the basis:",sum(irr_total)
write(*,*) "New dimension of the basis including spin :",Nbasis

  
!Here we organize the one particle basis elements
Nbasis=sum(irr_total)*2 !irr_total(1)!sum(irr_total)*2!irr_total(1)!irr_total(1)*2      ! sum(irr_total)*2   !        !Spin
write(*,*) "New dimension of the basis:",sum(irr_total)
write(*,*) "New dimension of the basis including spin :",Nbasis

!NExtOrb=0
allocate(ielement(1:Nbasis))
call InitIelement(ielement)
! .. Basis orthonormalization: Xi
! write(*,*) "inn"
! call Get_Xi(Bsplinex,ielement,iprim)
! write(*,*) "out"

! !.. Transformation Matrix 
! allocate(TMa(1:Ndim,1:Ndim))
! allocate(arr(1:Ndim))
! arr=0.d0
! TMa=0.d0
! TMa(1,1)=1.d0
! Do j=2,Ndim  !UkXi
!    Do i=2,Ndim  !UkBE
!       If(iha(i).eq.iha(j))then
!          If(ielement(iep(i))%spin.eq.ielement(iep(j))%spin)Then
!             TMa(i,j)=Xi(iep(i),iep(j))
!          End IF
!       End IF
!    End Do
! End Do

  allocate(Spq(1:Nbasis,1:Nbasis))
  Spq=0.d0
  Do j=1,Nbasis
     Do i=1,Nbasis
        Spq(j,i)=GET_INTEGRALS_E(j,i,0,0,1,Bsplinex,ielement,iprim)
     End Do
  End Do


call STEX_BASIS(ielement)

!.. Now we build the system of equations
allocate(Hirr(1:Ndim,1:Ndim),Over(1:Ndim,1:Ndim),En(1:Ndim))  
Over=0.d0
Hirr=0.d0
En=0.d0

!.. Hartree Fock Energy
NucRep=(7*7.d0/(2*0.54875d0))*(1.d0/1.8897259886)
write(*,*) "Nuclear repulsion",NucRep
!NucRep=23.6261352
ivalue2=OneBody_STEX(1,1,0,Bsplinex,ielement,iprim)
ivalue3=TwoBody_STEX(1,1,0,Bsplinex,ielement,iprim)
ivalue=ivalue2+ivalue3
write(*,*) "Hartree Fock Energy",ivalue+NucRep,ivalue2, ivalue3, NucRep
pause
!stop


!.. STEX MATRIX using ukrmol basis only
iflag=0


Do j=1,Ndim
   write(*,*) "j",j,Ndim
   Do i=1,Ndim
      Over(i,j)=1.d0*Kron(i,j)
      !      Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX_DM(j,i,1,Bsplinex,ielement,iprim)
      !      ivalue=OneBody_STEX_DM(j,i,iflag,Bsplinex,ielement,iprim)
      !      ivalue=ivalue+TwoBody_STEX_NO(j,i,0,Bsplinex,ielement,iprim)

           Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX(j,i,1,Bsplinex,ielement,iprim)     
           ivalue=OneBody_STEX(j,i,iflag,Bsplinex,ielement,iprim)
           ivalue=ivalue+TwoBody_STEX(j,i,0,Bsplinex,ielement,iprim)
      
      ! Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX_NO(j,i,1,Bsplinex,ielement,iprim)
      ! ivalue=OneBody_STEX_NO(j,i,iflag,Bsplinex,ielement,iprim)
      ! ivalue=ivalue+TwoBody_STEX_NO(j,i,0,Bsplinex,ielement,iprim)

      
      Hirr(j,i)=ivalue+Over(i,j)*NucRep
   End Do
End Do

! write(*,*) Ndim
! open(unit=10,file="matrizH.dat")
! do j=1,Ndim
! do i=1,Ndim
!    write(10,*) Hirr(j,i)
! End Do
! End Do

  
!check symmetry
H0d=0.d0
ivalue=0.d0
Do j=1,Nbasis!NprimUkrmol
   Do i=1,Nbasis!NprimUkrmol
      H0d=H0d+abs(Hirr(j,i)-Hirr(i,j))
      !      write(*,*) j,i,abs(Hirr(j,i)-Hirr(i,j))
      ivalue=ivalue+abs(Over(j,i)-Over(i,j))
   End DO
   !      write(*,*) H0d
   !     pause
End Do
write(*,*) "assymetry", H0d,ivalue*0.5/Nbasis


! Hirr=matmul(Hirr,TMa)
! Hirr=matmul(transpose(TMa),Hirr)
! Over=matmul(Over,TMa)
! Over=matmul(transpose(TMa),Over)


LWORK=Ndim*3-1
allocate(WORK(1:LWORK))
WORK=0.d0
INFO=0
call dsygv(1,'N','U',Ndim,Hirr,Ndim,Over,Ndim,En,WORK,LWORK,INFO)
write(*,*) "energies INFOpp",INFO
!pause
Do i=1,Ndim
   Write(*,*) i,En(i)
End do
open(unit=2, file="energy_ukrmol.dat")
!open(unit=2, file="energy_ukrmol_extended.dat")
Do i=1,Ndim
   write(2,*) En(i)
End Do
close(2)
deallocate(WORK)


!.. Comparison of the energies
open(unit=1, file="energy_ukrmol.dat")
open(unit=2, file="energy_ukrmol_extended.dat")
H0d=0.d0
r=0.d0
Do i=1,Ndim
   read(2,*) ivalue2
   read(1,*) ivalue
   write(*,*) ivalue-ivalue2,iha(i),iep(i),i
   If(abs(ivalue-ivalue2).gt.r)then
      r=abs(ivalue-ivalue2)
      ia=i
   End IF
   H0d=H0d+abs(ivalue-ivalue2)
End Do
close(1)
close(2)

 write(*,*) "differences eigenvallues",H0d
 write(*,*) ia,r


  

call mpi_mod_finalize


!!$  call GetRunTimeParameters( FileName, nSize )

stop

contains





  
  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( FileName, nSize )
    !
    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: FileName
    integer                      , intent(out) :: nSize
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Template Programs which diagonalizes a random matrix"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help" , "Print Command Usage" )
    call List.Add( "-o"     , "Output File" ,"eval", "optional" )
    call List.Add( "-n"     , "Matrix dim"  , 100  , "required" )

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if

    call List.Get( "-o",  strnBuf  )
    allocate(FileName,source=trim(adjustl(strnBuf)))

    call List.Get( "-n", nSize )

    call List.Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
    write(OUTPUT_UNIT,"(a)"   ) "Output File : "//FileName
    write(OUTPUT_UNIT,"(a,i0)") "Matrix size : ",nSize
    !
  end subroutine GetRunTimeParameters


  
end program ProgramTemplate

