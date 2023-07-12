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
module ModuleOrbitalBasis

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  !L0.0
  use ModuleErrorHandling
  use ModuleString
  !L0.1
  use ModuleGroups
  use ModuleBspline
  
  implicit none
  private
  
  character(len=*), public, parameter :: BASIS_SUBDIR = "basis/"

  type, private :: ClassOrbitalBasis
     character(len=:), allocatable :: StorageDir
     integer                       :: lmax
     type(ClassBSpline)            :: BsSet
     integer         , allocatable :: nMolv(:)
     integer         , allocatable :: nInAv(:)
     integer         , allocatable :: nActv(:)
     integer         , allocatable :: nHybv(:)
   contains
     procedure, public :: Init          => ClassOrbitalBasis_Init
     procedure, public :: Free          => ClassOrbitalBasis_Free
     procedure, public :: Save          => ClassOrbitalBasis_Save
     procedure, public :: Load          => ClassOrbitalBasis_Load
     procedure, public :: SetStorageDir => ClassOrbitalBasis_SetStorageDir
     procedure, public :: SetNactive    => ClassOrbitalBasis_SetNactive
     procedure, public :: Getlmax       => ClassOrbitalBasis_Getlmax
     procedure, public :: GetNOrb       => ClassOrbitalBasis_GetNOrb
     procedure, public :: GetBsOrder    => ClassOrbitalBasis_GetBsOrder
  end type ClassOrbitalBasis
  type(ClassOrbitalBasis), public :: OrbitalBasis

  public :: ReadDALTONIrreps
  
contains

  subroutine ClassOrbitalBasis_Init(self, lmax, nMolv, nHybv, BsSet )
    class(ClassOrbitalBasis), intent(inout) :: self
    integer                 , intent(in)    :: lmax
    integer, allocatable    , intent(in)    :: nMolv(:)
    integer, allocatable    , intent(in)    :: nHybv(:)
    type(ClassBspline)      , intent(in)    :: BsSet
    integer :: iIrr, iAstra
    self%lmax = lmax
    self%BsSet = BsSet
    allocate(self%nMolv,source=nMolv)
    allocate(self%nHybv,source=nHybv)
  end subroutine ClassOrbitalBasis_Init
  
  subroutine ClassOrbitalBasis_Free(self )
    class(ClassOrbitalBasis), intent(inout) :: self
    self%lmax = -1
    call self%BsSet%Free()
    if(allocated(self%nMolv))deallocate(self%nMolv)
    if(allocated(self%nHybv))deallocate(self%nHybv)
    if(allocated(self%nActv))deallocate(self%nActv)
    if(allocated(self%StorageDir)) deallocate(self%StorageDir)
  end subroutine ClassOrbitalBasis_Free
  
  subroutine ClassOrbitalBasis_SetStorageDir(self, Dir )
    class(ClassOrbitalBasis), intent(inout) :: self
    character(len=*)        , intent(in)    :: Dir
    allocate(self%StorageDir, source=Dir)
  end subroutine ClassOrbitalBasis_SetStorageDir
  
  subroutine ClassOrbitalBasis_Save(self)
    class(ClassOrbitalBasis), intent(inout) :: self
    integer :: uid
    call Execute_Command_Line("mkdir -p "//self%StorageDir//"/"//BASIS_SUBDIR)
    open(newunit=uid,file=self%StorageDir//"/"//BASIS_SUBDIR//"OrbitalBasis",&
         form="formatted",status="unknown")
    write(uid,"(*(x,i0))") self%lmax, size(self%nMolv)
    write(uid,"(*(x,i0))") self%nMolv
    write(uid,"(*(x,i0))") self%nHybv
    call self%BsSet%Save(uid)
    close(uid)
  end subroutine ClassOrbitalBasis_Save
  
  subroutine ClassOrbitalBasis_Load(self)
    class(ClassOrbitalBasis), intent(inout) :: self
    integer :: uid,nIrr
    open(newunit=uid,file=self%StorageDir//"/"//BASIS_SUBDIR//"OrbitalBasis",&
         form="formatted",status="old")
    call self%Free()
    read(uid,*) self%lmax, nIrr
    allocate(self%nMolv(nIrr),self%nHybv(nIrr))
    read(uid,*) self%nMolv
    read(uid,*) self%nHybv
    call self%BsSet%Load(uid)
    close(uid)
  end subroutine ClassOrbitalBasis_Load
  
  subroutine ClassOrbitalBasis_SetNactive(self,ninactive,nactive)
    class(ClassOrbitalBasis), intent(inout) :: self
    integer                 , intent(in)    :: ninactive(:)
    integer                 , intent(in)    :: nactive(:)
    allocate(self%ninav,source=ninactive)
    allocate(self%nactv,source=nactive)
  end subroutine ClassOrbitalBasis_SetNactive
  
  integer function ClassOrbitalBasis_GetNOrb(self, Irr, sType) result( res )
    class(ClassOrbitalBasis), intent(inout) :: self
    type(ClassIrrep)        , intent(in)    :: Irr
    character(len=*)        , intent(in)    :: sType
    integer :: iIrr
    iIrr = GlobalGroup%GetIrrepIndex(Irr)
    if    ( sType .is. "molecular" )then
       res = self%nMolv( iIrr )
    elseif( sType .is. "inactive" )then
       res = self%nInav( iIrr )
    elseif( sType .is. "active" )then
       res = self%nActv( iIrr )
    elseif( sType .is. "virtual" )then
       res = self%nMolv( iIrr ) - self%nActv( iIrr ) - self%nInav( iIrr )
    elseif( sType .is. "hybrid" )then
       res = self%nHybv( iIrr )
    elseif( sType .is. "spline" )then
       res = self%BsSet%GetNBsplines() - 2*self%BsSet%GetOrder() + 2
    endif
  end function ClassOrbitalBasis_GetNOrb

  integer function ClassOrbitalBasis_GetLMax(self) result( res )
    class(ClassOrbitalBasis), intent(inout) :: self
    res = self%lmax
  end function ClassOrbitalBasis_GetLMax

  integer function ClassOrbitalBasis_GetBsOrder(self) result( res )
    class(ClassOrbitalBasis), intent(inout) :: self
    res = self%BsSet%GetOrder()
  end function ClassOrbitalBasis_GetBsOrder


  subroutine ReadDALTONIrreps( ASTRAIrr_from_DALTONIrr, QCDir )

    integer, allocatable, intent(out) :: ASTRAIrr_from_DALTONIrr(:)
    character(len=*)    , intent(in)  :: QCDir

    character(len=*), parameter :: DALTON_OUTPUT = "DALTON.OUT"
    character(len=*), parameter :: IRREP_LINE_ID = "The irrep name for each symmetry:"

    integer :: uid, iostat, i, iIrrASTRA, iIrrDALTON, nIrr
    character(len=1000) :: line, iomsg
    character(len=3) :: IrrName
    type(ClassIrrep), pointer :: IrrepPtr

    open(newunit = uid, &
         file    = QCDir//"/"//DALTON_OUTPUT, &
         status  ="old", &
         form    ="formatted",&
         action  ="read", &
         iostat  = iostat, &
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Error opening "//DALTON_OUTPUT//" "//trim(iomsg))
       stop
    endif
    do
       read(uid,"(a)",iostat=iostat)line
       if(iostat/=0)then
          call ErrorMessage("irrep names not found in "//DALTON_OUTPUT)
          stop
       endif
       i=index(line,IRREP_LINE_ID)
       if(i>0)exit
    enddo
    line=adjustl(line(i+len(IRREP_LINE_ID):))

    nIrr = GlobalGroup%GetNIrreps()

    if(allocated(ASTRAIrr_from_DALTONIrr))deallocate(ASTRAIrr_from_DALTONIrr)
    allocate(ASTRAIrr_from_DALTONIrr(nIrr))
    ASTRAIrr_from_DALTONIrr=0
    write(*,"(a)") "n Irr DALTON,  n Irr ASTRA,   Irr Name"
    do iIrrDALTON = 1, nIrr
       i=index(line,":")
       line=adjustl(line(i+1:))
       IrrName = line(1:3)
       IrrepPtr => GlobalGroup%GetIrrep(trim(line(1:3)))
       iIrrASTRA= GlobalGroup%GetIrrepIndex(IrrepPtr)
       write(*,"(2(x,i2),2x,a3)") iIrrDALTON, iIrrASTRA, IrrName
       ASTRAIrr_from_DALTONIrr(iIrrDALTON)=iIrrASTRA
    enddo
    close(uid)


  end subroutine ReadDALTONIrreps

  
end module ModuleOrbitalBasis
