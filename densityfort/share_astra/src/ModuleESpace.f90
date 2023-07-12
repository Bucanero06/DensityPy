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
module ModuleESpace
  
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleIO
  use ModuleString
  use ModuleMatrix
  use ModuleGroups
  use ModuleParentIons
  use ModuleSymESpace
  
  implicit none
  private

  character(len=*), parameter :: CLOSE_COUPLING_DIR      = "CloseCoupling"
  character(len=*), parameter :: HamiltonianFileRootName = "H_"
  character(len=*), parameter :: KinEnergyFileRootName   = "K_"
  character(len=*), parameter :: OverlapFileRootName     = "S_"
  character(len=*), parameter :: DipoleLenFileRootName   = "DipoleLen_"
  character(len=*), parameter :: DipoleVelFileRootName   = "DipoleVel_"
  character(len=*), parameter :: QCFileExtension         = "QC"

  type, public :: ClassESpace
     private
     logical                           :: MULTIPLICITY_IS_SET = .FALSE.
     integer                           :: Multiplicity
     type(ClassGroup)                  :: Group
     integer                           :: Charge
     integer                           :: Lmax
     integer                           :: NIrreps
     integer                           :: NParentIons
     type(ClassParentIon), pointer     :: ParentIonList(:)
     integer             , allocatable :: nRefOccupied(:)
     type(ClasssymCCSpace), allocatable :: IrrepSpaceVec(:)
     character(len=:)    , allocatable :: StorageDir
     !
     ! Note: nRefOccupied is the number of orbitals, per each 
     ! irrep, that is doubly occupied in a closed-shell reference
     ! configuration. This parameter is useful only in the case
     ! of the static-exchange approximation where such reference
     ! configuration is needed to define the nature of all parent ions
   contains
     procedure, public :: Initialized      => ClassESpaceInitialized
     procedure, public :: ParseConfigFile  => ClassESpaceParseConfigFile
     procedure, public :: SetMultiplicity  => ClassESpaceSetMultiplicity
     procedure, public :: SetRootDir       => ClassESpaceSetRootDir
     procedure, public :: Show             => ClassESpaceShow
     procedure, public :: GetGroup         => ClassESpaceGetGroup
     procedure, public :: GetGroupName     => ClassESpaceGetGroupName
     procedure, public :: GetPionList      => ClassESpaceGetPionList
     procedure, public :: GetLmax          => ClassESpaceGetLmax
     procedure, public :: GetStorageDir    => ClassESpaceGetStorageDir
     procedure, public :: GetPionCharge    => ClassESpaceGetPionCharge
     procedure, public :: GetNIrreps       => ClassESpaceGetNIrreps
     procedure, public :: GetnRefOccupied  => ClassESpaceGetnRefOccupied
     generic  , public :: GetSymElectSpace => ClassESpaceGetSymElectSpace, &
          ClassESpaceGetSymElectSpaceByIrrepName
     procedure, private:: ClassESpaceGetSymElectSpace
     procedure, private:: ClassESpaceGetSymElectSpaceByIrrepName
     procedure, public :: CheckSymmetry    => ClassESpaceCheckSymmetry
  end type ClassESpace

  
  !> Class of the operators electronic space blocks: S, H, D, etc.
  type, public :: ClassESBlock
     private
     type(ClassESpace)     , pointer     :: Space
     !> Blocks correspondind to all the combinations of
     !! symmetric electronic space bras and kets, defined
     !! for the point group and the close-coupling expansion.
     type(ClassSESSESBlock), allocatable :: Block(:,:)
     logical                             :: initialized
   contains
     procedure, public :: Init => ClassESBlockInit
     procedure, public :: Load => ClassESBlockLoad
     procedure, public :: Save => ClassESBlockSave
     procedure, public :: Free => ClassESBlockFree
     final             :: ClassESBlockFinal
  end type ClassESBlock

  public :: GetCloseCouplingDir

contains

  subroutine ClassESpaceSetMultiplicity( Self, Multiplicity )
    class(ClassESpace), intent(inout) :: Self
    integer           , intent(in)    :: Multiplicity
    Self%Multiplicity = Multiplicity
    Self%MULTIPLICITY_IS_SET = .TRUE.
  end subroutine ClassESpaceSetMultiplicity

  subroutine ClassESpaceSetRootDir( self, RootDir )
    class(ClassESpace), intent(inout) :: self
    character(len=*)  , intent(in)    :: RootDir
    allocate( self%StorageDir, source = RootDir//"/"//CLOSE_COUPLING_DIR )
  end subroutine ClassESpaceSetRootDir

  function ClassESpaceGetGroup( Self ) result(Group)
    class(ClassESpace), intent(in) :: Self
    type(ClassGroup), pointer :: Group
    allocate(Group,source=Self%Group)
  end function ClassESpaceGetGroup

  function ClassESpaceGetPionList( Self ) result(IonList)
    class(ClassESpace), intent(in) :: Self
    type(ClassParentIon), pointer :: IonList(:)
    allocate(IonList,source=Self%ParentIonList)
  end function ClassESpaceGetPionList

  function ClassESpaceGetLmax( Self ) result(Lmax)
    class(ClassESpace), intent(in) :: Self
    integer :: Lmax
    Lmax=Self%LMax
  end function ClassESpaceGetLmax

  function ClassESpaceGetPionCharge( Self ) result(PionCharge)
    class(ClassESpace), intent(in) :: Self
    integer :: PionCharge
    PionCharge=Self%Charge+1
  end function ClassESpaceGetPionCharge

  function ClassESpaceGetStorageDir( Self ) result(StorageDir)
    class(ClassESpace), intent(in) :: Self
    character(len=:), allocatable :: StorageDir
    allocate(StorageDir,source=Self%StorageDir)
  end function ClassESpaceGetStorageDir

  integer function ClassESpaceGetNIrreps( Self ) result(N)
    class(ClassESpace), intent(in) :: Self
    N = Self%NIrreps
  end function ClassESpaceGetNIrreps

  !> Gets the number of doubly occupied orbitals per each irrep,
  !! in the reference configuration which is assumed to be closed shell
  function ClassESpaceGetnRefOccupied( Self ) result(nRefOccupied)
    class(ClassESpace), intent(in) :: Self
    integer, allocatable :: nRefOccupied(:)
    allocate(nRefOccupied, source = Self%nRefOccupied)
  end function ClassESpaceGetnRefOccupied

  function ClassESpaceGetGroupName( Self ) result(GName)
    class(ClassESpace), intent(in) :: Self
    character(len=:), allocatable :: GName
    allocate( GName, source = Self%Group%GetName() )
  end function ClassESpaceGetGroupName

  !> Get the symmetric electronic space giving the corresponding index in the list.
  function ClassESpaceGetSymElectSpace( Self, iIrrep ) result( SymSpace )
    class(ClassESpace), target, intent(in) :: Self
    integer                   , intent(in) :: iIrrep
    type(ClasssymCCSpace), pointer :: SymSpace
    allocate( SymSpace, source = Self%IrrepSpaceVec(iIrrep) )
  end function ClassESpaceGetSymElectSpace

  !> Getthe symmetric electronic space giving the irreducible representation name.
  function ClassESpaceGetSymElectSpaceByIrrepName( Self, IrrepName ) result( SymSpace )
    class(ClassESpace), target, intent(in) :: Self
    character(len=*)          , intent(in) :: IrrepName
    type(ClasssymCCSpace), pointer :: SymSpace
    integer :: i
    do i = 1, Self%GetNIrreps()
       if ( Self%IrrepSpaceVec(i)%GetIrrepLabel() .is. IrrepName ) then
          allocate( SymSpace, source = Self%IrrepSpaceVec(i) )
          return
       end if
    end do
  end function ClassESpaceGetSymElectSpaceByIrrepName

  subroutine ClassESpaceCheckSymmetry( Space, SymLabel )
    class(ClassESpace), target, intent(in) :: Space
    character(len=*)          , intent(in) :: SymLabel
    integer                       :: NumSymElectSpace, i
    type(ClasssymCCSpace), pointer :: SymElectSpace
    logical                       :: PresentIrrep
    NumSymElectSpace = Space%GetNIrreps()
    PresentIrrep = .false.
    do i = 1, NumSymElectSpace
       allocate( SymElectSpace, source = Space%GetSymElectSpace(i) )
       if ( trim(adjustl(SymLabel)) .is. SymElectSpace%GetIrrepLabel() ) then
          PresentIrrep = .true.
          return
       end if
       deallocate( SymElectSpace )
    end do
    if ( .not.PresentIrrep )then
       call ErrorMessage( &
         ' The irreducible representation '//SymLabel//&
         ' is not present in the close-coupling expansion for the '//&
         Space%Group%GetName()//&
         ' group.' )
       stop
    endif
  end subroutine ClassESpaceCheckSymmetry

  function GetCloseCouplingDir() result(Dir)
    character(len=:), allocatable :: Dir
    allocate( Dir, source = CLOSE_COUPLING_DIR )
  end function GetCloseCouplingDir

  subroutine ClassESpaceParseConfigFile( Self, &
       FileName, VERBOUS_ )
    class(ClassESpace), intent(inout) :: Self
    character(len=*)  , intent(in)    :: FileName
    integer, optional , intent(in)    :: VERBOUS_

    character(len=5)              :: irrepLabel
    character(len=:), allocatable :: FullText, GroupLabel, GeneralLabel, ParentIonList
    character(len=:), allocatable :: snRefOccupied
    type(ClassIrrep), pointer     :: irrepv(:)
    integer                       :: iIrrep, ichar, iflag, iSym, iIon, iIrr, iostat
    character(len=:), allocatable :: sIon, sOrb, symlabel, AuxCfgFile, PassFile
    character(len=*), parameter   :: EOLN_STRN = "<\EOLN>"
    integer :: VERBOUS, uid
    logical :: USE_FULL_BASIS
    VERBOUS = 0
    if(present(VERBOUS_))VERBOUS=VERBOUS_

    !.. Load the Configuration file
    call GetFullText( FileName, FullText, EOLN_STRN )
    call SetStringToUppercase( FullText )

    !.. Open Auxiliary Config File
    call Execute_Command_Line("mkdir -p "//self%StorageDir)
    AuxCfgFile=self%StorageDir//"/.CLSCPLNG.INP"
    open(newunit = uid       , &
         file    = AuxCfgFile, &
         form    ="formatted", &
         status  ="unknown"  , &
         action  ="write"    )
         
    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "GROUP", GroupLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GroupLabel)) call Assert("Group label missing in "//trim(FileName))
    write(*  ,"(a)") " Group               : "//GroupLabel
    write(uid,"(a)") " Group = "//GroupLabel
    call Self%Group%init( GroupLabel )

    call FetchGlobalVariable( FullText, "USE_FULL_BASIS", GeneralLabel, EOLN_STRN = EOLN_STRN )
    USE_FULL_BASIS=.FALSE.
    PassFile=FileName
    if(allocated(GeneralLabel))then
       write(*,"(a)")   " Use full basis      : "//GeneralLabel
       write(uid,"(a)") " USE_FULL_BASIS = "     //GeneralLabel
       USE_FULL_BASIS = ( GeneralLabel .is. "TRUE" )
       PassFile=AuxCfgFile
    end if

    call FetchGlobalVariable( FullText, "LMAX", GeneralLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GeneralLabel)) call Assert("lmax label missing in "//trim(FileName))
    write(*  ,"(a)") " Lmax                : "//GeneralLabel
    write(uid,"(a)") " LMAX = "//GeneralLabel
    read(GeneralLabel,*) Self%Lmax

    call FetchGlobalVariable( FullText, "CHARGE", GeneralLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GeneralLabel)) call Assert("CHARGE label missing in "//trim(FileName))
    write(*  ,"(a)") " Tot. Charge         : "//GeneralLabel
    write(uid,"(a)") " CHARGE = "//GeneralLabel
    read(GeneralLabel,*) Self%Charge

    call FetchGlobalVariable( FullText, "PARENT_ION_LIST", ParentIonList, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(ParentIonList)) call Assert("Parent-ion list missing in "//trim(FileName))
    write(*  ,"(a)") " List of parent ions : "//trim(ParentIonList)
    write(uid,"(a)") " PARENT_ION_LIST = "//trim(ParentIonList)
    Self%NParentIons = nTokens(ParentIonList)
    allocate(Self%ParentIonList(Self%NParentIons))
    do iIon = 1, Self%NParentIons
       call GetToken( ParentIonList, iIon, sIon )
       call Self%ParentIonList(iIon)%Init(Self%Group,sIon,Self%Charge)
       if(VERBOUS>0)call Self%ParentIonList(iIon)%Show()
    enddo

    call FetchGlobalVariable( FullText, "N_REF_OCCUPIED", snRefOccupied, EOLN_STRN = EOLN_STRN )
    write(uid,"(a)") " N_REF_OCCUPIED = "//trim(snRefOccupied)
    allocate(Self%nRefOccupied(Self%Group%GetNIrreps()))
    if(.not.allocated(snRefOccupied))then
       Self%nRefOccupied=0
    else
       write(*,"(a)") " HF orbital occupancy: "//trim(snRefOccupied)
       if( Self%Group%GetNIrreps() /= nTokens(snRefOccupied) )then
          write(*,"(a)")"the number of N_REF_OCCUPIED orbitals does not match the number of irreps"
          stop
       endif
       do iIrr = 1, size(Self%nRefOccupied)
          call GetToken( snRefOccupied, iIrr, sOrb )
          read(sOrb,*,iostat=iostat)Self%nRefOccupied(iIrr)
          if(iostat/=0)then
             write(*,"(a)")"syntax error in N_REF_OCCUPIED"
             stop
          endif
       enddo
    endif

    if(Self%MULTIPLICITY_IS_SET)then

       !.. Determines the number of symmetric spaces
       irrepv => Self%Group%GetIrrepList()
       Self%NIrreps = 0
       write(*,"(a)",advance="no") " List of Symmetries  :"
       do iIrrep = 1, size( irrepv )
          irrepLabel = adjustl(irrepv(iIrrep)%GetName())
          call SetStringToUppercase( irrepLabel )
          symlabel="["//AlphabeticNumber(Self%Multiplicity)//trim(irrepLabel)//"]"
          if( USE_FULL_BASIS )then
             write(*,"(a)",advance="no") " "//AlphabeticNumber(Self%Multiplicity)//trim(irrepLabel)
             write(uid,"(a)") symLabel//"{"
             Self%NIrreps = Self%NIrreps + 1
             do iIon=1,self%NParentIons
                write(uid,"(a)") self%ParentIonList(iIon)%GetLabel()//"( aiM viM hiG beS:ALL_XLM )"
             enddo
             write(uid,"(a)")"}"
          else
             ichar = index( FullText,symlabel)
             if(ichar>0)then
                Self%NIrreps = Self%NIrreps + 1
                write(*,"(a)",advance="no") " "//AlphabeticNumber(Self%Multiplicity)//trim(irrepLabel)
             endif
          endif

       enddo
       close(uid)
       write(*,*)

       !.. Parse each Symmetric Electronic Space
       allocate( Self%IrrepSpaceVec( Self%NIrreps ) )
       write(*,"(a,i0,a)") " Compute ",self%NIrreps," symmetries"
       
       iSym = 0
       do iIrrep = 1, size( irrepv )

          if(.not.USE_FULL_BASIS)then
             irrepLabel = adjustl(irrepv(iIrrep)%GetName())
             call SetStringToUppercase( irrepLabel )
             ichar = index( FullText, "["//AlphabeticNumber(Self%Multiplicity)//trim(irrepLabel)//"]" )
             if(ichar<1) cycle
          endif

          iSym = iSym + 1
          call Self%IrrepSpaceVec( iSym )%SetGroup( Self%Group )
          call Self%IrrepSpaceVec( iSym )%SetIrrep( irrepv( iIrrep ) )
          call Self%IrrepSpaceVec( iSym )%SetMultiplicity( Self%Multiplicity )
          call Self%IrrepSpaceVec( iSym )%SetRoot( Self%StorageDir )
          call Self%IrrepSpaceVec( iSym )%ParseConfigFile( PassFile, iflag )
          if(iflag /= 0 )call Assert( "Invalid syntax for "//&
               Self%IrrepSpaceVec( iSym )%GetLabel()//" in "//trim(FileName) )

       enddo
    endif
  end subroutine ClassESpaceParseConfigFile

  logical function ClassESpaceInitialized( Self ) result(res)
    class(ClassESpace), intent(in) :: Self
    res = allocated(Self%IrrepSpaceVec)
  end function ClassESpaceInitialized

  subroutine ClassESpaceShow( Self, unit )
    !> Class of the electronic space.
    class(ClassESpace), intent(in) :: Self
    integer, optional          , intent(in) :: unit
    integer :: i ,outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    write(outunit,"(a,i4)") "  Total charge       :",Self%Charge
    write(outunit,"(a,i4)") "  Space Multiplicity :",Self%Multiplicity
    write(outunit,"(a,i4)") "  Max Ang Momentum   :",Self%Lmax
    call Self%Group%Show()
    write(outunit,"(a,i4)") "  Number of Irreps   :",Self%NIrreps
    do i = 1, Self%NIrreps
       call Self%IrrepSpaceVec(i)%Show()
    end do
    write(outunit,"(a,a)" ) "  Space Storage Dir  :",Self%StorageDir
  end subroutine ClassESpaceShow


  !--------------------------------------------------
  ! Methods for ClassESBlock
  !-------------------------------------------------

  subroutine ClassESBlockFinal( Self )
    type(ClassESBlock) :: Self
    call Self%Free()
  end subroutine ClassESBlockFinal

  subroutine ClassESBlockFree( Self )
    class(ClassESBlock), intent(inout) :: Self
    Self%Space => NULL()
    if ( allocated(Self%Block) ) deallocate( Self%Block )
    self%initialized = .false.
  end subroutine ClassESBlockFree

  !> Initialize the electronic space operator matrix.
  subroutine ClassESBlockInit( Self, Space, OpLabel, IDLabel )
    class(ClassESBlock)       , intent(inout) :: Self
    class(ClassESpace), target, intent(in)    :: Space
    character(len=*)          , intent(in)    :: OpLabel
    character(len=*), optional, intent(in)    :: IDLabel
    integer :: i, j, NIrreps
    type(ClasssymCCSpace), pointer :: BraSymSpace, KetSymSpace
    Self%Space => Space
    NIrreps = Space%GetNIrreps()
    if ( allocated(Self%Block) ) deallocate( Self%Block )
    allocate( Self%Block(NIrreps,NIrreps) )
    do j = 1, NIrreps
       KetSymSpace => Space%GetSymElectSpace(j)
       do i = 1, NIrreps
          BraSymSpace => Space%GetSymElectSpace(i)
          if ( .not.ValidSymmetries( &
               OpLabel, BraSymSpace, KetSymSpace ) ) cycle
          call Self%Block(i,j)%Init( &
               BraSymSpace         , &
               KetSymSpace         , &
               OpLabel             , &
               IDLabel             )
       end do
    end do
    self%initialized = .true.
  end subroutine ClassESBlockInit

  !> Load the electronic space operator matrix.
  subroutine ClassESBlockLoad( Self )
    class(ClassESBlock), intent(inout) :: Self
    integer :: i, j, NIrreps
    NIrreps = self%Space%GetNIrreps()
    if ( allocated(Self%Block) ) deallocate( Self%Block )
    allocate( Self%Block(NIrreps,NIrreps) )
    do j = 1, NIrreps
       do i = 1, NIrreps
          if ( self%Block(i,j)%IsInitialized() ) call Self%Block(i,j)%Load()
       end do
    end do
  end subroutine ClassESBlockLoad

  !> Save the electronic space operator matrix.
  subroutine ClassESBlockSave( Self )
    class(ClassESBlock)        , intent(inout) :: Self
    integer :: i, j
    integer :: NIrreps
    NIrreps = Self%Space%GetNIrreps()
    do j = 1, NIrreps
       do i = 1, NIrreps
          call Self%Block(i,j)%Save( )
       end do
    end do
  end subroutine ClassESBlockSave

  
end module ModuleESpace
