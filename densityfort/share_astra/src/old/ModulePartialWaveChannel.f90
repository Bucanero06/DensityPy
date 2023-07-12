module ModulePartialWaveChannel

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleIO
  use ModuleString
  use ModuleConstants
  use ModuleGroups
  use ModuleSymmetryAdaptedSphericalHarmonics
  use ModuleParentIons
  use ModuleMatrix
  !use symba
  use ModuleGeneralChannel

  implicit none

  private

  character(len=*), parameter :: BsplineOrbDir = "BsplineOrbitals"


  type, public, extends( ClassGeneralChannel ) :: ClassPartialWaveChannel
     ! {{{ private attributes

     private
     type(ClassParentIon), pointer :: ParentIon
     integer                       :: TotalMultiplicity
     type(ClassIrrep),     pointer :: TotalIrrep
     type(ClassIrrep),     pointer :: OrbitalIrrep
     character(len=:), allocatable :: StorageDir
     integer                       :: Lmax
     type(ClassXlm)                :: Xlm

     ! }}}
   contains
     generic, public :: init                => ClassPartialWaveChannelInit
     generic, public :: show                => ClassPartialWaveChannelShow
     generic, public :: free                => ClassPartialWaveChannelFree
     generic, public :: GetLmax             => ClassPartialWaveChannelGetLmax
     generic, public :: GetL                => ClassPartialWaveChannelGetL
     generic, public :: GetM                => ClassPartialWaveChannelGetM
     generic, public :: GetStorageDir       => ClassPartialWaveChannelGetStorageDir
     generic, public :: GetTotIrrepName     => ClassPartialWaveChannelGetTotIrrepName
     generic, public :: GetTotIrrep         => ClassPartialWaveChannelGetTotIrrep
     generic, public :: GetOrbIrrep         => ClassPartialWaveChannelGetOrbIrrep
     generic, public :: GetTotMultiplicity  => ClassPartialWaveChannelGetTotMultiplicity
     generic, public :: GetXlmLabel         => ClassPartialWaveChannelGetXlmLabel
     generic, public :: GetSymLabel         => ClassPartialWaveChannelGetSymLabel
     procedure, public :: GetPILabel        => ClassPartialWaveChannelGetPILabel
     procedure, public :: GetPI             => ClassPartialWaveChannelGetPI
     procedure, public :: GetPINelect       => ClassPartialWaveChannelGetPINelect
     procedure, public :: GetPICharge       => ClassPartialWaveChannelGetPICharge
     procedure, public :: GetPIMultiplicity => ClassPartialWaveChannelGetPIMultiplicity
     generic, public :: GetPIEnergy         => ClassPartialWaveChannelGetPIEnergy, ClassPartialWaveChannelReadPIEnergy
     ! {{{ private procedures

     procedure, private :: ClassPartialWaveChannelInit
     procedure, private :: ClassPartialWaveChannelShow
     procedure, private :: ClassPartialWaveChannelFree
     procedure, private :: ClassPartialWaveChannelGetLMax
     procedure, private :: ClassPartialWaveChannelGetL
     procedure, private :: ClassPartialWaveChannelGetM
     procedure, private :: ClassPartialWaveChannelGetStorageDir
     procedure, private :: ClassPartialWaveChannelGetTotIrrepName
     procedure, private :: ClassPartialWaveChannelGetTotIrrep
     procedure, private :: ClassPartialWaveChannelGetOrbIrrep
     procedure, private :: ClassPartialWaveChannelGetTotMultiplicity
     procedure, private :: ClassPartialWaveChannelGetXlmLabel
     procedure, private :: ClassPartialWaveChannelGetSymLabel
     procedure, private :: ClassPartialWaveChannelGetPILabel
     procedure, private :: ClassPartialWaveChannelGetPI
     procedure, private :: ClassPartialWaveChannelGetPINelect
     procedure, private :: ClassPartialWaveChannelGetPICharge
     procedure, private :: ClassPartialWaveChannelGetPIMultiplicity
     procedure, private :: ClassPartialWaveChannelGetPIEnergy
     procedure, private :: ClassPartialWaveChannelReadPIEnergy

     ! }}}
  end type ClassPartialWaveChannel





  type, public :: ClassBsplineXlmBsplineXlmBlock

     type(ClassXlm)  , pointer     :: BraXlm
     type(ClassXlm)  , pointer     :: OperatorXlm
     character(len=:), allocatable :: OperatorLabel
     type(ClassXlm)  , pointer     :: KetXlm
     type(ClassMatrix)             :: Block

   contains

     generic, public :: init      => ClassBsplineXlmBsplineXlmBlockInit
     generic, public :: save      => ClassBsplineXlmBsplineXlmBlockSave
     generic, public :: ReadBlock => ClassBsplineXlmBsplineXlmBlockReadBlock
     generic, public :: GetFile   => ClassBsplineXlmBsplineXlmBlockGetFile
     generic, public :: Free      => ClassBsplineXlmBsplineXlmBlockFree

     procedure, private :: ClassBsplineXlmBsplineXlmBlockInit
     procedure, private :: ClassBsplineXlmBsplineXlmBlockSave
     procedure, private :: ClassBsplineXlmBsplineXlmBlockReadBlock
     procedure, private :: ClassBsplineXlmBsplineXlmBlockGetFile
     procedure, private :: ClassBsplineXlmBsplineXlmBlockFree
     
  end type ClassBsplineXlmBsplineXlmBlock



  public :: GetBsplineOrbDir



contains



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !      ClassPartialWaveChannel
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine ClassPartialWaveChannelInit( &
       Channel          , &
       ParentIon        , &
       Lmax             , &
       TotalMultiplicity, &
       TotalIrrep       , &
       Xlm              , &
       StorageDir         )
    class(ClassPartialWaveChannel) , intent(inout) :: Channel
    class(ClassParentIon), target  , intent(in)    :: ParentIon
    integer                        , intent(in)    :: Lmax
    integer                        , intent(in)    :: TotalMultiplicity
    class(ClassIrrep),     target  , intent(in)    :: TotalIrrep
    class(ClassXlm)                , intent(in)    :: Xlm
    character(len=*)               , intent(in)    :: StorageDir
    Channel.ParentIon          =>  ParentIon
    Channel.Lmax               =   Lmax
    Channel.TotalMultiplicity  =   TotalMultiplicity
    Channel.TotalIrrep         =>  TotalIrrep
    Channel.OrbitalIrrep       =>  TotalIrrep * ParentIon.GetIrrep()
    Channel.Xlm                =   Xlm
    if ( allocated(Channel.StorageDir) ) deallocate( Channel.StorageDir )
    allocate( Channel.StorageDir, source = StorageDir )
  end subroutine ClassPartialWaveChannelInit
  



  subroutine ClassPartialWaveChannelShow( Channel, unit )
    class(ClassPartialWaveChannel), intent(in) :: Channel
    integer, optional             , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    !
    write(outunit,"(a)") "Partial Wave Channels Info : "
    !
    call Channel.ParentIon.Show(unit)
    write(outunit,"(a,i4)"          ) "  Maximum Angular Momentum :",Channel.Lmax
    write(outunit,"(a,i4)"          ) "  Total Multiplicity :",Channel.TotalMultiplicity
    write(outunit,"(a)",advance="no") "  Total Symmetry     :"
    call Channel.TotalIrrep.show(unit)
    write(outunit,"(a)",advance="no") "  Orbital Symmetry     :"
    call Channel.OrbitalIrrep.show(unit)
    call Channel.Xlm.show(unit)
  end subroutine ClassPartialWaveChannelShow


  subroutine ClassPartialWaveChannelFree( Channel )
    class(ClassPartialWaveChannel), intent(inout) :: Channel
    Channel.ParentIon    => NULL()
    Channel.TotalIrrep   => NULL()
    Channel.OrbitalIrrep => NULL()
    Channel.Lmax = -1
    Channel.TotalMultiplicity = 0
    call Channel.Xlm.init(0,0)
  end subroutine ClassPartialWaveChannelFree


  integer function ClassPartialWaveChannelGetLMax( Channel ) result( LMax )
    class(ClassPartialWaveChannel), intent(in) :: Channel
    LMax = Channel.Lmax
  end function ClassPartialWaveChannelGetLMax



  integer function ClassPartialWaveChannelGetL( Channel ) result( L )
    class(ClassPartialWaveChannel), intent(in) :: Channel
    L = Channel.Xlm.GetL()
  end function ClassPartialWaveChannelGetL



  integer function ClassPartialWaveChannelGetM( Channel ) result( M )
    class(ClassPartialWaveChannel), intent(in) :: Channel
    M = Channel.Xlm.GetM()
  end function ClassPartialWaveChannelGetM


  function ClassPartialWaveChannelGetStorageDir( Self ) result( Dir )
    class(ClassPartialWaveChannel), intent(inout) :: Self
    character(len=:), allocatable :: Dir
    if ( .not.allocated(Self.StorageDir) ) call Assert( &
         'The storage directory is not allocated in PWC, impossible to fetch.' )
    allocate( Dir, source = Self.StorageDir )
  end function ClassPartialWaveChannelGetStorageDir



  function ClassPartialWaveChannelGetTotIrrepName( Self ) result(IrrepName)
    class(ClassPartialWaveChannel), intent(in) :: Self
    character(len=:), allocatable :: IrrepName
    allocate( IrrepName, source = Self.TotalIrrep.GetName() )
  end function ClassPartialWaveChannelGetTotIrrepName


  function ClassPartialWaveChannelGetTotIrrep( Self ) result(Irrep)
    class(ClassPartialWaveChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self.TotalIrrep
  end function ClassPartialWaveChannelGetTotIrrep


  function ClassPartialWaveChannelGetOrbIrrep( Self ) result(Irrep)
    class(ClassPartialWaveChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self.OrbitalIrrep
  end function ClassPartialWaveChannelGetOrbIrrep



  function ClassPartialWaveChannelGetTotMultiplicity( Self ) result(Mult)
    class(ClassPartialWaveChannel), intent(in) :: Self
    integer :: Mult
    Mult = Self.TotalMultiplicity
  end function ClassPartialWaveChannelGetTotMultiplicity




  function ClassPartialWaveChannelGetSymLabel( Self ) result(Label)
    class(ClassPartialWaveChannel), intent(in) :: Self
    character(len=:), allocatable :: Label
    allocate( Label, source = AlphabeticNumber(Self.GetTotMultiplicity())//Self.GetTotIrrepName() )
  end function ClassPartialWaveChannelGetSymLabel



  function ClassPartialWaveChannelGetXlmLabel( Channel ) result( Label )
    class(ClassPartialWaveChannel), intent(inout) :: Channel
    character(len=:), allocatable :: Label
    allocate( Label, source = Channel.Xlm.GetLabel() )
  end function ClassPartialWaveChannelGetXlmLabel




  function ClassPartialWaveChannelGetPILabel( PWC ) result( Label )
    class(ClassPartialWaveChannel), intent(inout) :: PWC
    character(len=:), allocatable :: Label
    allocate( Label, source = PWC.Parention.GetLabel() )
  end function ClassPartialWaveChannelGetPILabel



  function ClassPartialWaveChannelGetPI( PWC ) result( PI )
    class(ClassPartialWaveChannel), target, intent(in) :: PWC
    type(ClassparentIon), pointer :: PI
    allocate( PI, source = PWC.Parention )
  end function ClassPartialWaveChannelGetPI


  integer function ClassPartialWaveChannelGetPINelect( PWC ) result( PINE )
    class(ClassPartialWaveChannel), intent(in) :: PWC
    PINE = PWC.Parention.GetNelect()
  end function ClassPartialWaveChannelGetPINelect



  integer function ClassPartialWaveChannelGetPICharge( PWC ) result( PICharge )
    class(ClassPartialWaveChannel), intent(in) :: PWC
    PICharge = PWC.ParentIon.GetCharge()
  end function ClassPartialWaveChannelGetPICharge



  integer function ClassPartialWaveChannelGetPIMultiplicity( PWC ) result( PIMult )
    class(ClassPartialWaveChannel), intent(in) :: PWC
    PIMult = PWC.ParentIon.GetMultiplicity()
  end function ClassPartialWaveChannelGetPIMultiplicity



  real(kind(1d0)) function ClassPartialWaveChannelGetPIEnergy( PWC ) result( PIE )
    class(ClassPartialWaveChannel), intent(inout) :: PWC
    PIE = PWC.ParentIon.GetEnergy()
  end function ClassPartialWaveChannelGetPIEnergy



  real(kind(1d0)) function ClassPartialWaveChannelReadPIEnergy( PWC, StorageDir ) result( PIE )
    class(ClassPartialWaveChannel), intent(inout) :: PWC
    character(len=*)              , intent(in)    :: StorageDir
    PIE = PWC.ParentIon.ReadEnergy( StorageDir )
  end function ClassPartialWaveChannelReadPIEnergy


  
  subroutine ClassPartialWaveChannelFinal( Channel )
    type(ClassPartialWaveChannel) :: Channel
    call Channel.free()
  end subroutine ClassPartialWaveChannelFinal




!!$
!!$  subroutine PWCPWCDipoleXBlockBuildBasicMatElem( Dip, BasicMatrixElements, Basis, Force, Id, NumBsDropAtEnd, NumBsDropAtBegin, ConditionBs )
!!$    !
!!$    class(PWCPWCDipoleXBlock),         intent(inout) :: Dip
!!$    class( ClassBasicMatrixElements ), intent(in)    :: BasicMatrixElements
!!$    class(ClassBasis),                 intent(inout) :: Basis
!!$    logical,                           intent(in)    :: Force
!!$    character(len=*),                  intent(in)    :: Id
!!$    integer, optional,                 intent(in)    :: NumBsDropAtEnd
!!$    integer, optional,                 intent(in)    :: NumBsDropAtBegin
!!$    logical, optional,                 intent(in)    :: ConditionBs
!!$    !
!!$    type(ClassComplexMatrix) :: Overlap, CorrectMat
!!$    complex(kind(1d0)), allocatable :: OverlapArray(:,:) 
!!$    complex(kind(1d0)), allocatable :: NewArray(:,:)
!!$    complex(kind(1d0)), allocatable :: AuxArray(:,:), AuxArray2(:,:), CompNewArray(:,:)
!!$    integer :: NumBsDropEnding, NumBsDropBeginning, LBra, LKet, MBra, MKet
!!$    type(PWCPWCOverlapBlock) :: S
!!$    real(kind(1d0)) :: PIsDipole
!!$    !
!!$    !
!!$    if ( Dip.ValidBlock(Id,Force) ) return
!!$    !
!!$    if ( Hydrogen ) then
!!$       write(output_unit,*) "*** Computing Hydrogen atom ***"
!!$    elseif ( H2Plus ) then
!!$       write(output_unit,*) "*** Computing H2+ molecule ***"
!!$    end if
!!$    !
!!$    if ( present(NumBsDropAtEnd) ) then
!!$       NumBsDropEnding = NumBsDropAtEnd
!!$    else
!!$       NumBsDropEnding = 1
!!$    end if
!!$    !
!!$    if ( present(NumBsDropAtBegin) ) then
!!$       NumBsDropBeginning = NumBsDropAtBegin
!!$    else
!!$       NumBsDropBeginning = 0
!!$    end if
!!$    !
!!$    !..Read the block overlap
!!$    call S.Init( Dip.SpaceBra, Dip.SpaceKet, Dip.StorageDir, NumBsDropEnding, NumBsDropBeginning )  
!!$    !
!!$    call S.Build( BasicMatrixElements, Basis, Force, "DummyIdentifier", .true., NumBsDropEnding, NumBsDropBeginning, .false., ConditionBs ) 
!!$    !
!!$    call S.LoadBlock( Overlap )
!!$    !
!!$    call Overlap.FetchMatrix( OverlapArray )
!!$    !
!!$    !
!!$    !
!    if ( Dip.ValidBlock(Id,Force) ) then
!       return
!    else
!!$    !
!!$    PIsDipole = GetPIsDipoleX( Dip )
!!$    !
!!$    deallocate( Dip.Block )
!!$    allocate( Dip.Block(size(OverlapArray,1),size(OverlapArray,2)) )
!!$    Dip.Block = PIsDipole * OverlapArray
!!$    !
!!$    if ( (Dip.SpaceBra.GetPILabel() .is. &
!!$          Dip.SpaceKet.GetPILabel()) ) then
!!$       !
!!$       call ComputePWCPWCMonoElectronicDipoleX( Dip, BasicMatrixElements, NumBsDropBeginning, NumBsDropEnding, CorrectMat )
!!$       !
!!$       call CorrectMat.FetchMatrix( NewArray )
!!$       allocate( CompNewArray(size(NewArray,1),size(NewArray,2)) )
!!$       CompNewArray = NewArray
!!$       deallocate( NewArray )
!!$       !
!!$       LBra = Dip.SpaceBra.GetL()
!!$       MBra = Dip.SpaceBra.GetM()
!!$       LKet = Dip.SpaceKet.GetL()
!!$       MKet = Dip.SpaceKet.GetM()
!!$       !***
!!$       if ( ConditionBs ) then
!!$          call ApplyBraPreconditioning( CompNewArray, AuxArray, LBra, MBra )
!!$          call ApplyKetPreconditioning( AuxArray, AuxArray2, LKet, MKet )
!!$          deallocate( CompNewArray )
!!$          allocate( CompNewArray, source = AuxArray2 )
!!$          deallocate( AuxArray, AuxArray2 )
!!$       end if
!!$       !***
!!$       !
!!$       Dip.Block = Dip.Block + CompNewArray
!!$    end if
!!$    !
!!$    call Dip.SaveBlock(Id)
!!$
!!$    !
!!$    !
!!$  end subroutine PWCPWCDipoleXBlockBuildBasicMatElem





!!$  subroutine ComputePWCPWCMonoElectronicDipoleZ( Dip, BasicMatrixElements, NumBsDropBeginning, NumBsDropEnding, ResMat )
!!$    !
!!$    class(PWCPWCDipoleZBlock),       intent(in) :: Dip
!!$    class(ClassBasicMatrixElements), intent(in) :: BasicMatrixElements
!!$    integer,                        intent(in)  :: NumBsDropBeginning
!!$    integer,                        intent(in)  :: NumBsDropEnding
!!$    class(ClassComplexMatrix),      intent(out) :: ResMat
!!$    !
!!$    integer :: LKet, MKet, LBra, MBra
!!$    integer :: NumFunGABS, NumSkippedBs, NumGaussianGABS, NumExp
!!$    integer :: NCols, NRows
!!$    complex(kind(1d0)), allocatable :: AllBraXlmDipKetXlm(:,:), OutMat(:,:)
!!$    real(kind(1d0)), allocatable :: AuxKetLowerMat(:,:), AuxKetHigherMat(:,:)
!!$    real(kind(1d0)) :: ParHigherL, Par2_HigherL
!!$    type(ClassMatrix) :: RedDipLKetLowerMat, RedDipLKetHigherMat, CorrectRedDipLKetLowerMat, CorrectRedDipLKetHigherMat
!!$    character(len=:), allocatable :: Identifier
!!$    type(ClassComplexMatrix) :: Mat1, Mat2
!!$    logical :: IsZero, IsAsym
!!$    !
!!$    allocate( Identifier, source = "z" )
!!$    !
!!$    LKet = Dip.SpaceKet.GetL()
!!$    MKet = Dip.SpaceKet.GetM()
!!$    LBra = Dip.SpaceBra.GetL()
!!$    MBra = Dip.SpaceBra.GetM()
!!$    !
!!$    NumFunGABS = Dip.SpaceKet.GetNumGABS()
!!$    NumSkippedBs = Dip.SpaceKet.GetNumSkippedBs()
!!$    NumGaussianGABS = Dip.SpaceKet.GetNumGaussians()
!!$    NumExp = Dip.SpaceKet.GetGaussianNumExponents()
!!$    !
!!$    !
!!$    NCols = NumFunGABS - NumBsDropEnding - &
!!$         (NumGaussianGABS + NumSkippedBs + NumBsDropBeginning)
!!$    NRows = NCols
!!$    ! 
!!$    !
!!$    if ( DipoleZBlockLabel .is. LenDipoleZBlockLabel ) then
!!$       call BasicMatrixElements.FetchGeneralLengthMatrix( RedDipLKetHigherMat, LBra, LKet )  
!!$    elseif ( DipoleZBlockLabel .is. VelDipoleZBlockLabel ) then
!!$       call BasicMatrixElements.FetchGeneralVelocityMatrix(  RedDipLKetHigherMat, LBra, LKet )
!!$    end if
!!$    !
!!$    call RedDipLKetHigherMat.GetSubMatrix( &
!!$         NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$         NumFunGABS - NumBsDropEnding, &
!!$         NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$         NumFunGABS - NumBsDropEnding, &
!!$         CorrectRedDipLKetHigherMat )
!!$    call RedDipLKetHigherMat.Free()
!!$    call CorrectRedDipLKetHigherMat.FetchMatrix( AuxKetHigherMat )
!!$    call CorrectRedDipLKetHigherMat.Free()
!!$    !
!!$    !
!!$    !
!!$    allocate( AllBraXlmDipKetXlm(NRows,NCols) )
!!$    AllBraXlmDipKetXlm = Z0
!!$    !
!!$    !
!!$    !
!!$    call ComputePWCPWCXlmDipXlm( &
!!$         Identifier, &
!!$         LBra, MBra, &
!!$         LKet, MKet, &
!!$         AuxKetHigherMat, &
!!$         OutMat  )
!!$    !
!!$    AllBraXlmDipKetXlm(:,:) = OutMat(:,:) 
!!$    !
!!$    !
!!$    !
!!$    if ( allocated(OutMat) ) deallocate( OutMat )
!!$    !
!!$    !
!!$    ResMat = AllBraXlmDipKetXlm
!    !***
!    IsZero = ResMat.IsZero()
!    if ( .not.IsZero ) then
!       write(*,*) "Not Zero", LBra, MBra, LKet, MKet
!       stop
!    end if
    !***
!!$    deallocate( AllBraXlmDipKetXlm )
!!$    !
!!$  end subroutine ComputePWCPWCMonoElectronicDipoleZ


!!$
!!$  subroutine ComputePWCPWCMonoElectronicDipoleY( Dip, BasicMatrixElements, NumBsDropBeginning, NumBsDropEnding, ResMat )
!!$    !
!!$    class(PWCPWCDipoleYBlock),       intent(in) :: Dip
!!$    class(ClassBasicMatrixElements), intent(in) :: BasicMatrixElements
!!$    integer,                        intent(in)  :: NumBsDropBeginning
!!$    integer,                        intent(in)  :: NumBsDropEnding
!!$    class(ClassComplexMatrix),      intent(out) :: ResMat
!!$    !
!!$    integer :: LKet, MKet, LBra, MBra
!!$    integer :: NumFunGABS, NumSkippedBs, NumGaussianGABS, NumExp
!!$    integer :: NCols, NRows
!!$    complex(kind(1d0)), allocatable :: AllBraXlmDipKetXlm(:,:), OutMat(:,:)
!!$    real(kind(1d0)), allocatable :: AuxKetLowerMat(:,:), AuxKetHigherMat(:,:)
!!$    real(kind(1d0)) :: ParHigherL, Par2_HigherL
!!$    type(ClassMatrix) :: RedDipLKetLowerMat, RedDipLKetHigherMat, CorrectRedDipLKetLowerMat, CorrectRedDipLKetHigherMat
!!$    character(len=:), allocatable :: Identifier
!!$    type(ClassComplexMatrix) :: Mat1, Mat2
!!$    logical :: IsZero, IsAsym
!!$    !
!!$    allocate( Identifier, source = "y" )
!!$    !
!!$    LKet = Dip.SpaceKet.GetL()
!!$    MKet = Dip.SpaceKet.GetM()
!!$    LBra = Dip.SpaceBra.GetL()
!!$    MBra = Dip.SpaceBra.GetM()
!!$    !
!!$    NumFunGABS = Dip.SpaceKet.GetNumGABS()
!!$    NumSkippedBs = Dip.SpaceKet.GetNumSkippedBs()
!!$    NumGaussianGABS = Dip.SpaceKet.GetNumGaussians()
!!$    NumExp = Dip.SpaceKet.GetGaussianNumExponents()
!!$    !
!!$    !
!!$    NCols = NumFunGABS - NumBsDropEnding - &
!!$         (NumGaussianGABS + NumSkippedBs + NumBsDropBeginning)
!!$    NRows = NCols
!!$    ! 
!!$    !
!!$    if ( DipoleYBlockLabel .is. LenDipoleYBlockLabel ) then
!!$       call BasicMatrixElements.FetchGeneralLengthMatrix( RedDipLKetHigherMat, LBra, LKet )  
!!$    elseif ( DipoleYBlockLabel .is. VelDipoleYBlockLabel ) then
!!$       call BasicMatrixElements.FetchGeneralVelocityMatrix(  RedDipLKetHigherMat, LBra, LKet )
!!$    end if
!!$    !
!!$    call RedDipLKetHigherMat.GetSubMatrix( &
!!$         NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$         NumFunGABS - NumBsDropEnding, &
!!$         NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$         NumFunGABS - NumBsDropEnding, &
!!$         CorrectRedDipLKetHigherMat )
!!$    call RedDipLKetHigherMat.Free()
!!$    call CorrectRedDipLKetHigherMat.FetchMatrix( AuxKetHigherMat )
!!$    call CorrectRedDipLKetHigherMat.Free()
!!$    !
!!$    !
!!$    !
!!$    allocate( AllBraXlmDipKetXlm(NRows,NCols) )
!!$    AllBraXlmDipKetXlm = Z0
!!$    !
!!$    !
!!$    !
!!$    call ComputePWCPWCXlmDipXlm( &
!!$         Identifier, &
!!$         LBra, MBra, &
!!$         LKet, MKet, &
!!$         AuxKetHigherMat, &
!!$         OutMat  )
!!$    !
!!$    AllBraXlmDipKetXlm(:,:) = OutMat(:,:) 
!!$    !
!!$    !
!!$    !
!!$    if ( allocated(OutMat) ) deallocate( OutMat )
!!$    !
!!$    !
!!$    ResMat = AllBraXlmDipKetXlm
    !***
!    IsZero = ResMat.IsZero()
!    if ( .not.IsZero ) then
!       write(*,*) "Not Zero", LBra, MBra, LKet, MKet
!       stop
!    end if
!    !***
!!$    deallocate( AllBraXlmDipKetXlm )
!!$    !
!!$  end subroutine ComputePWCPWCMonoElectronicDipoleY



!!$  subroutine ComputePWCPWCMonoElectronicDipoleX( Dip, BasicMatrixElements, NumBsDropBeginning, NumBsDropEnding, ResMat )
!!$    !
!!$    class(PWCPWCDipoleXBlock),       intent(in) :: Dip
!!$    class(ClassBasicMatrixElements), intent(in) :: BasicMatrixElements
!!$    integer,                        intent(in)  :: NumBsDropBeginning
!!$    integer,                        intent(in)  :: NumBsDropEnding
!!$    class(ClassComplexMatrix),      intent(out) :: ResMat
!!$    !
!!$    integer :: LKet, MKet, LBra, MBra
!!$    integer :: NumFunGABS, NumSkippedBs, NumGaussianGABS, NumExp
!!$    integer :: NCols, NRows
!!$    complex(kind(1d0)), allocatable :: AllBraXlmDipKetXlm(:,:), OutMat(:,:)
!!$    real(kind(1d0)), allocatable :: AuxKetLowerMat(:,:), AuxKetHigherMat(:,:)
!!$    real(kind(1d0)) :: ParHigherL, Par2_HigherL
!!$    type(ClassMatrix) :: RedDipLKetLowerMat, RedDipLKetHigherMat, CorrectRedDipLKetLowerMat, CorrectRedDipLKetHigherMat
!!$    character(len=:), allocatable :: Identifier
!!$    type(ClassComplexMatrix) :: Mat1, Mat2
!!$    logical :: IsZero, IsAsym
!!$    !
!!$    allocate( Identifier, source = "x" )
!!$    !
!!$    LKet = Dip.SpaceKet.GetL()
!!$    MKet = Dip.SpaceKet.GetM()
!!$    LBra = Dip.SpaceBra.GetL()
!!$    MBra = Dip.SpaceBra.GetM()
!!$    !
!!$    NumFunGABS = Dip.SpaceKet.GetNumGABS()
!!$    NumSkippedBs = Dip.SpaceKet.GetNumSkippedBs()
!!$    NumGaussianGABS = Dip.SpaceKet.GetNumGaussians()
!!$    NumExp = Dip.SpaceKet.GetGaussianNumExponents()
!!$    !
!!$    !
!!$    NCols = NumFunGABS - NumBsDropEnding - &
!!$         (NumGaussianGABS + NumSkippedBs + NumBsDropBeginning)
!!$    NRows = NCols
!!$    ! 
!!$    !
!!$    if ( DipoleXBlockLabel .is. LenDipoleXBlockLabel ) then
!!$       call BasicMatrixElements.FetchGeneralLengthMatrix( RedDipLKetHigherMat, LBra, LKet )  
!!$    elseif ( DipoleXBlockLabel .is. VelDipoleXBlockLabel ) then
!!$       call BasicMatrixElements.FetchGeneralVelocityMatrix(  RedDipLKetHigherMat, LBra, LKet )
!!$    end if
!!$    !
!!$    call RedDipLKetHigherMat.GetSubMatrix( &
!!$         NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$         NumFunGABS - NumBsDropEnding, &
!!$         NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$         NumFunGABS - NumBsDropEnding, &
!!$         CorrectRedDipLKetHigherMat )
!!$    call RedDipLKetHigherMat.Free()
!!$    call CorrectRedDipLKetHigherMat.FetchMatrix( AuxKetHigherMat )
!!$    call CorrectRedDipLKetHigherMat.Free()
!!$    !
!!$    !
!!$    !
!!$    allocate( AllBraXlmDipKetXlm(NRows,NCols) )
!!$    AllBraXlmDipKetXlm = Z0
!!$    !
!!$    !
!!$    !
!!$    call ComputePWCPWCXlmDipXlm( &
!!$         Identifier, &
!!$         LBra, MBra, &
!!$         LKet, MKet, &
!!$         AuxKetHigherMat, &
!!$         OutMat  )
!!$    !
!!$    AllBraXlmDipKetXlm(:,:) = OutMat(:,:) 
!!$    !
!!$    !
!!$    !
!!$    if ( allocated(OutMat) ) deallocate( OutMat )
!!$    !
!!$    !
!!$    ResMat = AllBraXlmDipKetXlm
    !***
!    IsZero = ResMat.IsZero()
!    if ( .not.IsZero ) then
!       write(*,*) "Not Zero", LBra, MBra, LKet, MKet
!       stop
!    end if
    !***
!!$    deallocate( AllBraXlmDipKetXlm )
!!$    !
!!$  end subroutine ComputePWCPWCMonoElectronicDipoleX




!!$
!!$  subroutine ComputePWCPWCXlmDipXlm( Identifier, LBra, MBra, LKet, MKet, ReducedMatRow, OutMat )
!!$    !
!!$    character(len=*),                intent(in)  :: Identifier
!!$    integer,                         intent(in)  :: LBra
!!$    integer,                         intent(in)  :: MBra
!!$    integer,                         intent(in)  :: LKet
!!$    integer,                         intent(in)  :: MKet
!!$    real(kind(1d0)),                 intent(in)  :: ReducedMatRow(:,:)
!!$    complex(kind(1d0)), allocatable, intent(out) :: OutMat(:,:)
!!$    !
!!$    real(kind(1d0)) :: CMKet, CMBra
!!$    complex(kind(1d0)) :: KMKet, KMBra, DMKet, DMBra
!!$    complex(kind(1d0)), allocatable :: Res1(:,:), Res2(:,:), Res3(:,:), Res4(:,:)
!!$    !
!!$    if ( MKet == 0 ) then
!!$       CMKet = 1.d0
!!$       KMKet = Zi**(0.5d0*(CMKet-1.d0))/sqrt(2.d0)
!!$       DMKet = KMKet/sqrt(2.d0)
!!$    else
!!$       CMKet = dble(MKet/abs(MKet))
!!$       KMKet = Zi**(0.5d0*(CMKet-1.d0))/sqrt(2.d0)
!!$       DMKet = KMKet
!!$    end if
!!$    !
!!$    if ( MBra == 0 ) then
!!$       CMBra = 1.d0
!!$       KMBra = Zi**(0.5d0*(CMBra-1.d0))/sqrt(2.d0)
!!$       DMBra = KMBra/sqrt(2.d0)
!!$    else
!!$       CMBra = dble(MBra/abs(MBra))
!!$       KMBra = Zi**(0.5d0*(CMBra-1.d0))/sqrt(2.d0)
!!$       DMBra = KMBra
!!$    end if
!!$    !
!!$    allocate( OutMat(size(ReducedMatRow,1),size(ReducedMatRow,2)) )
!!$    OutMat = Z0
!!$    !
!!$    call ComputePWCPWCCartDipole( Identifier, LBra, LKet,  abs(MBra),  abs(MKet), ReducedMatRow, Res1 )
!!$    call ComputePWCPWCCartDipole( Identifier, LBra, LKet,  abs(MBra), -abs(MKet), ReducedMatRow, Res2 )
!!$    call ComputePWCPWCCartDipole( Identifier, LBra, LKet, -abs(MBra),  abs(MKet), ReducedMatRow, Res3 )
!!$    call ComputePWCPWCCartDipole( Identifier, LBra, LKet, -abs(MBra), -abs(MKet), ReducedMatRow, Res4 )
!!$    !
!!$    OutMat = conjg(DMBra)*DMKet * ( Res1 + (-1.d0)**MKet * CMKet * Res2  + &
!!$         (-1.d0)**MBra * CMBra * Res3 + (-1.d0)**(MKet+MBra) * CMBra*CMKet * Res4 )
!!$    !
!!$  end subroutine ComputePWCPWCXlmDipXlm



!!$  subroutine ComputePWCPWCCartDipole( Identifier, LBra, LKet, MBra, MKet, ReducedMatRow, Res )
!!$    !
!!$    character(len=*),                intent(in)  :: Identifier
!!$    integer,                         intent(in)  :: LBra
!!$    integer,                         intent(in)  :: LKet
!!$    integer,                         intent(in)  :: MBra
!!$    integer,                         intent(in)  :: MKet
!!$    real(kind(1d0)),                 intent(in)  :: ReducedMatRow(:,:)
!!$    complex(kind(1d0)), allocatable, intent(out) :: Res(:,:)
!!$    !
!!$    complex(kind(1d0)) :: Factor
!!$    !
!!$    if ( Identifier .is. "x" ) then
!!$       Factor = ( CGC(LKet,1,LBra,MKet,-1)*delta(MBra,MKet-1) - CGC(LKet,1,LBra,MKet,1)*delta(MBra,MKet+1) )/sqrt(2.d0)/sqrt(dble(2*LBra+1))
!!$    elseif ( Identifier .is. "y" ) then
!!$       Factor = ( CGC(LKet,1,LBra,MKet,-1)*delta(MBra,MKet-1) + CGC(LKet,1,LBra,MKet,1)*delta(MBra,MKet+1) )*Zi/sqrt(2.d0)/sqrt(dble(2*LBra+1))
!!$    elseif ( Identifier .is. "z" ) then
!!$       Factor = CGC(LKet,1,LBra,MKet,0)*delta(MBra,MKet)/sqrt(dble(2*LBra+1))
!!$    else
!!$       call Assert( "Invalid identifier for the dipole transition. It should be either 'x', 'y' or 'z'." )
!!$    end if
!!$    !
!!$    allocate( Res(size(ReducedMatRow,1),size(ReducedMatRow,2)) )
!!$    Res(:,:) = Factor * ReducedMatRow(:,:)
!!$    !
!!$  end subroutine ComputePWCPWCCartDipole





!!$  real(kind(1d0)) function GetPIsDipoleX( Dip ) result( Res )
!!$    !
!!$    class(PWCPWCDipoleXBlock), intent(in)     :: Dip
!!$    !
!!$    character(len=:), allocatable :: BraPILabel, KetPILabel, DipoleFile, PIStoreDir
!!$    integer :: uid, i, Counter, iostat
!!$    integer, allocatable :: MonomialIndex(:), MonomialIndexRef(:)
!!$    integer, parameter :: NMon = 3
!!$    character(len=IOMSG_LENGTH) :: iomsg
!!$    logical :: Opened
!!$    !
!!$    !*** Must eliminate any reference to absolute paths
!!$    allocate( PIStoreDir, source = Dip.StorageDir//"../"//"../"//"../"//"../ParentIons/" )
!!$    !
!!$    allocate( BraPILabel, source = Dip.SpaceBra.GetPILabel() )
!!$    allocate( KetPILabel, source = Dip.SpaceKet.GetPILabel() )
!!$    !
!!$    allocate( DipoleFile, source = PIStoreDir//"ParentIon_"//BraPILabel//"_"//KetPILabel//"_CartMultipolesL1" )
!!$    !
!!$    call OpenFile( DipoleFile, uid, "read", "formatted" )
!!$    !Skip first three lines corresponding to labels x, y, z
!!$    do i = 1, NMon
!!$       read(uid,*)
!!$    end do
!!$    !
!!$    allocate( MonomialIndex(NMon) )
!!$    MonomialIndex = -1
!!$    !
!!$    allocate( MonomialIndexRef(NMon) )
!!$    MonomialIndexRef(1) = 1
!!$    MonomialIndexRef(2) = 0
!!$    MonomialIndexRef(3) = 0
!!$    !
!!$    Counter = 0
!!$    !
!!$    do
!!$       read(uid,*) MonomialIndex(1), MonomialIndex(2), MonomialIndex(3), Res
!!$       Counter = Counter + 1
!!$       if ( (MonomialIndex(1) == MonomialIndexRef(1)) .and. &
!!$            (MonomialIndex(2) == MonomialIndexRef(2)) .and. &
!!$            (MonomialIndex(3) == MonomialIndexRef(3)) ) then
!!$          !
!!$          close( uid )
!!$          return
!!$          !
!!$       end if
!!$       if ( Counter > NMon ) then
!!$          call Assert( "There are more dipoles to compare than 3." )
!!$       end if
!!$    end do
!!$    !
!!$    INQUIRE(&
!!$         UNIT  = uid , &
!!$         OPENED= Opened  , &
!!$         IOSTAT= iostat  , &
!!$         IOMSG = iomsg     )
!!$    if(iostat/=0) call Assert(iomsg)
!!$    if ( Opened ) close( uid )
!!$    !
!!$  end function GetPIsDipoleX



!!$  real(kind(1d0)) function GetPIsDipoleY( Dip ) result( Res )
!!$    !
!!$    class(PWCPWCDipoleYBlock), intent(in)     :: Dip
!!$    !
!!$    character(len=:), allocatable :: BraPILabel, KetPILabel, DipoleFile, PIStoreDir
!!$    integer :: uid, i, Counter, iostat
!!$    integer, allocatable :: MonomialIndex(:), MonomialIndexRef(:)
!!$    integer, parameter :: NMon = 3
!!$    character(len=IOMSG_LENGTH) :: iomsg
!!$    logical :: Opened
!!$    !
!!$    allocate( PIStoreDir, source = Dip.StorageDir//"../"//"../"//"../"//"../ParentIons/" )
!!$    !
!!$    allocate( BraPILabel, source = Dip.SpaceBra.GetPILabel() )
!!$    allocate( KetPILabel, source = Dip.SpaceKet.GetPILabel() )
!!$    !
!!$    allocate( DipoleFile, source = PIStoreDir//"ParentIon_"//BraPILabel//"_"//KetPILabel//"_CartMultipolesL1" )
!!$    !
!!$    call OpenFile( DipoleFile, uid, "read", "formatted" )
!!$    !Skip first three lines corresponding to labels x, y, z
!!$    do i = 1, NMon
!!$       read(uid,*)
!!$    end do
!!$    !
!!$    allocate( MonomialIndex(NMon) )
!!$    MonomialIndex = -1
!!$    !
!!$    allocate( MonomialIndexRef(NMon) )
!!$    MonomialIndexRef(1) = 0
!!$    MonomialIndexRef(2) = 1
!!$    MonomialIndexRef(3) = 0
!!$    !
!!$    Counter = 0
!!$    !
!!$    do
!!$       read(uid,*) MonomialIndex(1), MonomialIndex(2), MonomialIndex(3), Res
!!$       Counter = Counter + 1
!!$       if ( (MonomialIndex(1) == MonomialIndexRef(1)) .and. &
!!$            (MonomialIndex(2) == MonomialIndexRef(2)) .and. &
!!$            (MonomialIndex(3) == MonomialIndexRef(3)) ) then
!!$          !
!!$          close( uid )
!!$          return
!!$          !
!!$       end if
!!$       if ( Counter > NMon ) then
!!$          call Assert( "There are more dipoles to compare than 3." )
!!$       end if
!!$    end do
!!$    !
!!$    INQUIRE(&
!!$         UNIT  = uid , &
!!$         OPENED= Opened  , &
!!$         IOSTAT= iostat  , &
!!$         IOMSG = iomsg     )
!!$    if(iostat/=0) call Assert(iomsg)
!!$    if ( Opened ) close( uid )
!!$    !
!!$  end function GetPIsDipoleY



!!$  real(kind(1d0)) function GetPIsDipoleZ( Dip ) result( Res )
!!$    !
!!$    class(PWCPWCDipoleZBlock), intent(in)     :: Dip
!!$    !
!!$    character(len=:), allocatable :: BraPILabel, KetPILabel, DipoleFile, PIStoreDir
!!$    integer :: uid, i, Counter, iostat
!!$    integer, allocatable :: MonomialIndex(:), MonomialIndexRef(:)
!!$    integer, parameter :: NMon = 3
!!$    character(len=IOMSG_LENGTH) :: iomsg
!!$    logical :: Opened
!!$    !
!!$    allocate( PIStoreDir, source = Dip.StorageDir//"../"//"../"//"../"//"../ParentIons/" )
!!$    !
!!$    allocate( BraPILabel, source = Dip.SpaceBra.GetPILabel() )
!!$    allocate( KetPILabel, source = Dip.SpaceKet.GetPILabel() )
!!$    !
!!$    allocate( DipoleFile, source = PIStoreDir//"ParentIon_"//BraPILabel//"_"//KetPILabel//"_CartMultipolesL1" )
!!$    !
!!$    call OpenFile( DipoleFile, uid, "read", "formatted" )
!!$    !Skip first three lines corresponding to labels x, y, z
!!$    do i = 1, NMon
!!$       read(uid,*)
!!$    end do
!!$    !
!!$    allocate( MonomialIndex(NMon) )
!!$    MonomialIndex = -1
!!$    !
!!$    allocate( MonomialIndexRef(NMon) )
!!$    MonomialIndexRef(1) = 0
!!$    MonomialIndexRef(2) = 0
!!$    MonomialIndexRef(3) = 1
!!$    !
!!$    Counter = 0
!!$    !
!!$    do
!!$       read(uid,*) MonomialIndex(1), MonomialIndex(2), MonomialIndex(3), Res
!!$       Counter = Counter + 1
!!$       if ( (MonomialIndex(1) == MonomialIndexRef(1)) .and. &
!!$            (MonomialIndex(2) == MonomialIndexRef(2)) .and. &
!!$            (MonomialIndex(3) == MonomialIndexRef(3)) ) then
!!$          !
!!$          close( uid )
!!$          return
!!$          !
!!$       end if
!!$       if ( Counter > NMon ) then
!!$          call Assert( "There are more dipoles to compare than 3." )
!!$       end if
!!$    end do
!!$    !
!!$    INQUIRE(&
!!$         UNIT  = uid , &
!!$         OPENED= Opened  , &
!!$         IOSTAT= iostat  , &
!!$         IOMSG = iomsg     )
!!$    if(iostat/=0) call Assert(iomsg)
!!$    if ( Opened ) close( uid )
!!$    !
!!$  end function GetPIsDipoleZ



!!$  subroutine PWCPWCDipoleYBlockBuildBasicMatElem( Dip, BasicMatrixElements, Basis, Force, Id, NumBsDropAtEnd, NumBsDropAtBegin, ConditionBs )
!!$    !
!!$    class(PWCPWCDipoleYBlock),         intent(inout) :: Dip
!!$    class( ClassBasicMatrixElements ), intent(in)    :: BasicMatrixElements
!!$    class(ClassBasis),                 intent(inout) :: Basis
!!$    logical,                           intent(in)    :: Force
!!$    character(len=*),                  intent(in)    :: Id
!!$    integer, optional,                 intent(in)    :: NumBsDropAtEnd
!!$    integer, optional,                 intent(in)    :: NumBsDropAtBegin
!!$    logical, optional,                 intent(in)    :: ConditionBs
!!$    !
!!$    type(ClassComplexMatrix) :: Overlap, CorrectMat
!!$    complex(kind(1d0)), allocatable :: OverlapArray(:,:) 
!!$    complex(kind(1d0)), allocatable :: NewArray(:,:)
!!$    complex(kind(1d0)), allocatable :: AuxArray(:,:), AuxArray2(:,:), CompNewArray(:,:)
!!$    integer :: NumBsDropEnding, NumBsDropBeginning, LBra, LKet, MBra, MKet
!!$    type(PWCPWCOverlapBlock) :: S
!!$    real(kind(1d0)) :: PIsDipole
!!$    !
!!$    !
!!$    if ( Hydrogen ) then
!!$       write(output_unit,*) "*** Computing Hydrogen atom ***"
!!$    elseif ( H2Plus ) then
!!$       write(output_unit,*) "*** Computing H2+ molecule ***"
!!$    end if
!!$    !
!!$    if ( present(NumBsDropAtEnd) ) then
!!$       NumBsDropEnding = NumBsDropAtEnd
!!$    else
!!$       NumBsDropEnding = 1
!!$    end if
!!$    !
!!$    if ( present(NumBsDropAtBegin) ) then
!!$       NumBsDropBeginning = NumBsDropAtBegin
!!$    else
!!$       NumBsDropBeginning = 0
!!$    end if
!!$    !
!!$    !..Read the block overlap
!!$    call S.Init( Dip.SpaceBra, Dip.SpaceKet, Dip.StorageDir, NumBsDropEnding, NumBsDropBeginning )  
!!$    !
!!$    call S.Build( BasicMatrixElements, Basis, Force, "DummyIdentifier", .true., NumBsDropEnding, NumBsDropBeginning, .false., ConditionBs ) 
!!$    !
!!$    call S.LoadBlock( Overlap )
!!$    !
!!$    call Overlap.FetchMatrix( OverlapArray )
!!$    !
!!$    !
!!$    !
!!$    if ( Dip.ValidBlock(Id,Force) ) then
!!$       return
!!$    else
!!$       !
!!$       PIsDipole = GetPIsDipoleY( Dip )
!!$       !
!!$       !
!!$       if ( (Dip.SpaceBra.GetPILabel() .is. &
!!$            Dip.SpaceKet.GetPILabel()) ) then
!!$          !
!!$          call ComputePWCPWCMonoElectronicDipoleY( Dip, BasicMatrixElements, NumBsDropBeginning, NumBsDropEnding, CorrectMat )
!!$          !
!!$          call CorrectMat.FetchMatrix( NewArray )
!!$          allocate( CompNewArray(size(NewArray,1),size(NewArray,2)) )
!!$          CompNewArray = NewArray
!!$          deallocate( NewArray )
!!$          !
!!$          LBra = Dip.SpaceBra.GetL()
!!$          MBra = Dip.SpaceBra.GetM()
!!$          LKet = Dip.SpaceKet.GetL()
!!$          MKet = Dip.SpaceKet.GetM()
!!$          !***
!!$          if ( ConditionBs ) then
!!$             call ApplyBraPreconditioning( CompNewArray, AuxArray, LBra, MBra )
!!$             call ApplyKetPreconditioning( AuxArray, AuxArray2, LKet, MKet )
!!$             deallocate( CompNewArray )
!!$             allocate( CompNewArray, source = AuxArray2 )
!!$             deallocate( AuxArray, AuxArray2 )
!!$             deallocate( Dip.Block )
!!$             allocate( Dip.Block(size(CompNewArray,1),size(CompNewArray,2)) )
!!$             Dip.Block = Z0
!!$          end if
!!$          !***
!!$          !
!!$          Dip.Block = PIsDipole * OverlapArray + CompNewArray
!!$             call Dip.SaveBlock(Id)
!!$          !
!!$       else
!!$          !
!!$          deallocate( Dip.Block )
!!$          allocate( Dip.Block(size(OverlapArray,1),size(OverlapArray,2)) )
!!$          Dip.Block = PIsDipole * OverlapArray
!!$          call Dip.SaveBlock(Id)
!!$          !
!!$       end if
!!$    end if
!!$    !
!!$    !
!!$    !
!!$  end subroutine PWCPWCDipoleYBlockBuildBasicMatElem



!!$  subroutine PWCPWCDipoleZBlockBuildBasicMatElem( Dip, BasicMatrixElements, Basis, Force, Id, NumBsDropAtEnd, NumBsDropAtBegin, ConditionBs )
!!$    !
!!$    class(PWCPWCDipoleZBlock),         intent(inout) :: Dip
!!$    class( ClassBasicMatrixElements ), intent(in)    :: BasicMatrixElements
!!$    class(ClassBasis),                 intent(inout) :: Basis
!!$    logical,                           intent(in)    :: Force
!!$    character(len=*),                  intent(in)    :: Id
!!$    integer, optional,                 intent(in)    :: NumBsDropAtEnd
!!$    integer, optional,                 intent(in)    :: NumBsDropAtBegin
!!$    logical, optional,                 intent(in)    :: ConditionBs
!!$    !
!!$    type(ClassComplexMatrix) :: Overlap, CorrectMat
!!$    complex(kind(1d0)), allocatable :: OverlapArray(:,:) 
!!$    complex(kind(1d0)), allocatable :: NewArray(:,:)
!!$    complex(kind(1d0)), allocatable :: AuxArray(:,:), AuxArray2(:,:), CompNewArray(:,:)
!!$    integer :: NumBsDropEnding, NumBsDropBeginning, LBra, LKet, MBra, MKet
!!$    type(PWCPWCOverlapBlock) :: S
!!$    real(kind(1d0)) :: PIsDipole
!!$    !
!!$    !
!!$    if ( Hydrogen ) then
!!$       write(output_unit,*) "*** Computing Hydrogen atom ***"
!!$    elseif ( H2Plus ) then
!!$       write(output_unit,*) "*** Computing H2+ molecule ***"
!!$    end if
!!$    !
!!$    if ( present(NumBsDropAtEnd) ) then
!!$       NumBsDropEnding = NumBsDropAtEnd
!!$    else
!!$       NumBsDropEnding = 1
!!$    end if
!!$    !
!!$    if ( present(NumBsDropAtBegin) ) then
!!$       NumBsDropBeginning = NumBsDropAtBegin
!!$    else
!!$       NumBsDropBeginning = 0
!!$    end if
!!$    !
!!$    !..Read the block overlap
!!$    call S.Init( Dip.SpaceBra, Dip.SpaceKet, Dip.StorageDir, NumBsDropEnding, NumBsDropBeginning )  
!!$    !
!!$    call S.Build( BasicMatrixElements, Basis, Force, "DummyIdentifier", .true., NumBsDropEnding, NumBsDropBeginning, .false., ConditionBs ) 
!!$    !
!!$    call S.LoadBlock( Overlap )
!!$    !
!!$    call Overlap.FetchMatrix( OverlapArray )
!!$    !
!!$    !
!!$    !
!!$    if ( Dip.ValidBlock(Id,Force) ) then
!!$       return
!!$    else
!!$       !
!!$       PIsDipole = GetPIsDipoleZ( Dip )
!!$       !
!!$       !
!!$       if ( (Dip.SpaceBra.GetPILabel() .is. &
!!$            Dip.SpaceKet.GetPILabel()) ) then
!!$          !
!!$          call ComputePWCPWCMonoElectronicDipoleZ( Dip, BasicMatrixElements, NumBsDropBeginning, NumBsDropEnding, CorrectMat )
!!$          !
!!$          call CorrectMat.FetchMatrix( NewArray )
!!$          allocate( CompNewArray(size(NewArray,1),size(NewArray,2)) )
!!$          CompNewArray = NewArray
!!$          deallocate( NewArray )
!!$          !
!!$          LBra = Dip.SpaceBra.GetL()
!!$          MBra = Dip.SpaceBra.GetM()
!!$          LKet = Dip.SpaceKet.GetL()
!!$          MKet = Dip.SpaceKet.GetM()
!!$          !***
!!$          if ( ConditionBs ) then
!!$             call ApplyBraPreconditioning( CompNewArray, AuxArray, LBra, MBra )
!!$             call ApplyKetPreconditioning( AuxArray, AuxArray2, LKet, MKet )
!!$             deallocate( CompNewArray )
!!$             allocate( CompNewArray, source = AuxArray2 )
!!$             deallocate( AuxArray, AuxArray2 )
!!$             deallocate( Dip.Block )
!!$             allocate( Dip.Block(size(CompNewArray,1),size(CompNewArray,2)) )
!!$             Dip.Block = Z0
!!$          end if
!!$          !***
!!$          !
!!$          Dip.Block = PIsDipole * OverlapArray + CompNewArray
!!$          call Dip.SaveBlock(Id)
!!$          !
!!$       else
!!$          !
!!$          deallocate( Dip.Block )
!!$          allocate( Dip.Block(size(OverlapArray,1),size(OverlapArray,2)) )
!!$          Dip.Block = PIsDipole * OverlapArray
!!$          call Dip.SaveBlock(Id)
!!$          !
!!$       end if
!!$    end if
!!$    !
!!$    !
!!$    !
!!$  end subroutine PWCPWCDipoleZBlockBuildBasicMatElem



!!$  subroutine PWCPWCCAPBlockBuild( CAP, File, Basis, Force, NumBsDropAtEnd, NumBsDropAtBegin )
!!$    !
!!$    class(PWCPWCCAPBlock),             intent(inout) :: CAP
!!$    character(len=*),                  intent(in)    :: File
!!$    type(ClassBasis),                  intent(inout) :: Basis
!!$    logical,                           intent(in)    :: Force
!!$    integer, optional,                 intent(in)    :: NumBsDropAtEnd
!!$    integer, optional,                 intent(in)    :: NumBsDropAtBegin
!!$    !
!!$    character(len=:), allocatable :: Id
!!$    integer :: NumGaussianGABS, NumSkippedBs, NumFunGABS, NumBsDropEnding, NumBsDropBeginning, L
!!$    type(ClassAbsorptionPotential) :: AbsPot
!!$    type(ClassMatrix) :: MatArray, CorrectMat 
!!$    real(kind(1d0)), allocatable :: NewArray(:,:)
!!$    !
!!$    CAP.Block = 0.d0
!!$    !
!!$    NumBsDropEnding = 1
!!$    if ( present(NumBsDropAtEnd) ) NumBsDropEnding = NumBsDropAtEnd
!!$    !
!!$    NumBsDropBeginning = 0
!!$    if ( present(NumBsDropAtBegin) ) NumBsDropBeginning = NumBsDropAtBegin
!!$    !
!!$    !.. Load the absorption potentials
!!$    call AbsPot.Parse( File )
!!$    call AbsPot.Init( Basis )
!!$    call AbsPot.Read( Basis.GetDir() )
!!$    !
!!$    L = CAP.SpaceKet.GetL()
!!$    allocate( Id, source = "CAP" )
!!$    !
!!$    if( (CAP.SpaceBra.GetPILabel() .is. CAP.SpaceKet.GetPILabel()) .and. &
!!$         CAP.SpaceBra.Xlm.Is(CAP.SpaceKet.Xlm) .and. &
!!$         .not.CAP.ValidBlock(Id,Force) )then
!!$       !
!!$       NumFunGABS = CAP.SpaceBra.Basis.GetNFun()
!!$       NumSkippedBs = CAP.SpaceBra.Basis.GetNBSplineSkip()
!!$       NumGaussianGABS = CAP.SpaceBra.Basis.GetNFun("G")
!!$       !
!!$       !*** Assumes only one CAP
!!$       call AbsPot.Fetch( 1, MatArray, L, REP_LEVEL_TOTAL )
!!$       !
!!$       call MatArray.GetSubMatrix( &
!!$            NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$            NumFunGABS - NumBsDropEnding, &
!!$            NumGaussianGABS + 1 + NumSkippedBs + NumBsDropBeginning, &
!!$            NumFunGABS - NumBsDropEnding, &
!!$            CorrectMat )
!!$       !
!!$       call CorrectMat.FetchMatrix( NewArray )
!!$       CAP.Block = NewArray
!!$       !
!!$    end if
!!$
!!$    call CAP.SaveBlock(Id)
!!$    !
!!$  end subroutine PWCPWCCAPBlockBuild




  !------------------------------------------------------
  ! Methods for ClassBsplineXlmBsplineXlmBlock
  !------------------------------------------------------



  subroutine ClassBsplineXlmBsplineXlmBlockInit( this, BraXlm, OpXlm, OpLabel, KetXlm, Block )
    class(ClassBsplineXlmBsplineXlmBlock), intent(inout) :: this
    class(ClassXlm)  , target, intent(in) :: BraXlm
    class(ClassXlm)  , target, intent(in) :: OpXlm
    character(len=*)         , intent(in) :: OpLabel
    class(ClassXlm)  , target, intent(in) :: KetXlm
    class(ClassMatrix), optional, intent(in) :: Block
    
    this.BraXlm       => BraXlm
    this.OperatorXlm  => OpXlm
    this.KetXlm       => KetXlm
    if(allocated(this.OperatorLabel))deallocate(this.OperatorLabel)
    allocate(this.OperatorLabel,source=OpLabel)
    if(present(Block)) this.Block = Block
    
  end subroutine ClassBsplineXlmBsplineXlmBlockInit



  subroutine ClassBsplineXlmBsplineXlmBlockFree( this )
    class(ClassBsplineXlmBsplineXlmBlock), intent(inout) :: this
    this.BraXlm       => NULL()
    this.OperatorXlm  => NULL()
    this.KetXlm       => NULL()
    if(allocated(this.OperatorLabel))deallocate(this.OperatorLabel)
    call this.Block.Free()
  end subroutine ClassBsplineXlmBsplineXlmBlockFree



  function ClassBsplineXlmBsplineXlmBlockGetFile( this ) result( FileName )
    class(ClassBsplineXlmBsplineXlmBlock), intent(in) :: this
    character(len=:), allocatable :: FileName
    allocate(FileName, source=&
         this.OperatorLabel//&
         AlphabeticNumber(this.OperatorXlm.Getl())//"."//AlphabeticNumber(this.OperatorXlm.Getm())//"_"//&
         AlphabeticNumber(this.BraXlm.Getl())//"."//AlphabeticNumber(this.BraXlm.Getm())//"_"//&
         AlphabeticNumber(this.KetXlm.Getl())//"."//AlphabeticNumber(this.KetXlm.Getm()) )
  end function ClassBsplineXlmBsplineXlmBlockGetFile



  subroutine ClassBsplineXlmBsplineXlmBlockSave( this, Storage )
    class(ClassBsplineXlmBsplineXlmBlock), intent(inout) :: this
    character(len=*)                     , intent(in)    :: Storage
    
    character(len=:), allocatable :: FileName
    integer                       :: uid
    
    call execute_command_line("mkdir -p "//AddSlash(Storage)//AddSlash(BsplineOrbDir) )

    allocate(FileName, source= AddSlash(Storage)//AddSlash(BsplineOrbDir)//this.GetFile())
    call OpenFile( FileName, uid, "write", "formatted" )
    call this.Block.write(uid)
    close(uid)
    
  end subroutine ClassBsplineXlmBsplineXlmBlockSave



  subroutine ClassBsplineXlmBsplineXlmBlockReadBlock( this, Storage, Matrix )
    class(ClassBsplineXlmBsplineXlmBlock), intent(in)  :: this
    character(len=*)                     , intent(in)  :: Storage
    class(ClassMatrix)                   , intent(out) :: Matrix
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid
    !
    allocate(FileName, source= AddSlash(Storage)//AddSlash(BsplineOrbDir)//this.GetFile())
    call OpenFile( FileName, uid, "read", "formatted" )
    call Matrix.Read(uid)
    close(uid)
    !
  end subroutine ClassBsplineXlmBsplineXlmBlockReadBlock



  function GetBsplineOrbDir() result(Dir)
    character(len=:), allocatable :: Dir
    allocate( Dir, source = BsplineOrbDir )
  end function GetBsplineOrbDir



end module ModulePartialWaveChannel
