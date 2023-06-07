module ModuleElementaryElectronicSpaceOperators


  use, intrinsic :: ISO_FORTRAN_ENV


  use ModuleErrorHandling
  use ModuleString
  use ModuleMatrix
  use ModuleIO
  use ModuleConstants
  use symba
  use ModuleSymmetryAdaptedSphericalHarmonics
  use ModuleSymmetricLocalizedStates
  use ModuleShortRangeChannel
  use ModulePartialWaveChannel
  use ModuleShortRangeOrbitals
  use ModuleGroups
  use ModuleParentIons
  use ModuleElectronicSpace


  !... Gauge labels
  character(len=*), public, parameter :: LengthLabel   = 'Len'
  character(len=*), public, parameter :: VelocityLabel = 'Vel'

  !..File's extensions.
  character(len=*), public, parameter :: FileExtensionLocPWC = 'LocPWC'
  character(len=*), public, parameter :: FileExtensionPWCLoc = 'PWCLoc'
  character(len=*), public, parameter :: FileExtensionSRCPWC = 'SRCPWC'
  character(len=*), public, parameter :: FileExtensionPWCSRC = 'PWCSRC'
  character(len=*), public, parameter :: FileExtensionPWCPWC = 'PWCPWC'

  !..Operator labels
  character(len=*), public, parameter :: OverlapLabel     = 'S'
  character(len=*), public, parameter :: HamiltonianLabel = 'H'
  character(len=*), public, parameter :: KinEnergyLabel   = 'K'
  character(len=*), public, parameter :: DIPOLE_LABEL     = 'Dipole'
  character(len=*), public, parameter :: CAPLabel         = 'CAP'

  !> Quenched Hamiltonian label.
  character(len=*), public, parameter :: QHLabel          = 'QH'

  !> Feshbach partitioning Hamiltonian label.
  character(len=*), public, parameter :: FeshbachHLabel   = 'FeshbachH'
  character(len=*), public, parameter :: MultipoleLabel   = 'Multipole'


  !> Class of a general operator block.
  type, private :: ClassGeneralBlock
     private
     !> Operator label
     character(len=:) , allocatable :: OpLabel
     !> Operator associated Xlm
     type(ClassXlm)   ,  pointer    :: OpXlm 
     !> Array containing the operator representation in the current space.
     type(ClassMatrix)              :: Block
     logical                        :: BoxStatesOnly       = .FALSE.
     logical                        :: BraBoxStatesOnly    = .FALSE.
     logical                        :: KetBoxStatesOnly    = .FALSE.
     logical                        :: AsympBsPermutation  = .FALSE.
     logical                        :: FormattedWrite      = .FALSE.
     character(len=:), allocatable  :: FileExtension
     integer                        :: NumBsNotOverlappingDiffuseBra = -1
     integer                        :: NumBsNotOverlappingDiffuseKet = -1
   contains
     generic, public       :: Init                    => ClassGeneralBlockInit
     generic, public       :: SetFileExtension        => ClassGeneralBlockSetFileExtension
     generic, public       :: GetBlockName            => ClassGeneralBlockGetBlockName
     generic, public       :: ValidBlock              => ClassGeneralBlockValidBlock
     generic, public       :: Save                    => ClassGeneralBlockSaveBlock
     generic, public       :: FetchBlock              => ClassGeneralBlockFetchBlock
     generic, public       :: Read                    =>  ClassGeneralBlockReadBlock, ClassGeneralBlockReadQCBlock
     generic, public       :: SetBoxOnly              =>  ClassGeneralBlockSetBoxOnly
     generic, public       :: SetBraBoxOnly           =>  ClassGeneralBlockSetBraBoxOnly
     generic, public       :: SetKetBoxOnly           =>  ClassGeneralBlockSetKetBoxOnly
     generic, public       :: UnsetBoxOnly            =>  ClassGeneralBlockUnsetBoxOnly
     generic, public       :: UnsetBraBoxOnly         =>  ClassGeneralBlockUnsetBraBoxOnly
     generic, public       :: UnsetKetBoxOnly         =>  ClassGeneralBlockUnsetKetBoxOnly
     generic, public       :: SetAsympBsPermutation   =>  ClassGeneralBlockSetAsympBsPermutation
     generic, public       :: UnsetAsympBsPermutation =>  ClassGeneralBlockUnsetAsympBsPermutation
     generic, public       :: IsBoxOnly               =>  ClassGeneralBlockIsBoxOnly
     generic, public       :: IsBraBoxOnly            =>  ClassGeneralBlockIsBraBoxOnly
     generic, public       :: IsKetBoxOnly            =>  ClassGeneralBlockIsKetBoxOnly
     generic, public       :: SetFormattedWrite       =>  ClassGeneralBlockSetFormattedWrite
     generic, public       :: GetNumBsNotOvDiffBra    =>  ClassGeneralBlockGetNumBsNotOvDiffBra
     generic, public       :: GetNumBsNotOvDiffKet    =>  ClassGeneralBlockGetNumBsNotOvDiffKet
     generic, public       :: GetNRows                =>  ClassGeneralBlockGetNRows
     generic, public       :: GetNColumns             =>  ClassGeneralBlockGetNColumns
     generic, public       :: UnsetFormattedWrite     =>  ClassGeneralBlockUnsetFormattedWrite
     generic, public       :: Free                    => ClassGeneralBlockFree
     !
     procedure, private    :: ClassGeneralBlockInit
     procedure, private    :: ClassGeneralBlockSetFileExtension
     procedure, private    :: ClassGeneralBlockGetBlockName
     procedure, private    :: ClassGeneralBlockValidBlock
     procedure, private    :: ClassGeneralBlockSaveBlock
     procedure, private    :: ClassGeneralBlockReadBlock
     procedure, private    :: ClassGeneralBlockReadQCBlock
     procedure, private    :: ClassGeneralBlockFetchBlock
     procedure, private    :: ClassGeneralBlockSetBoxOnly
     procedure, private    :: ClassGeneralBlockSetBraBoxOnly
     procedure, private    :: ClassGeneralBlockSetKetBoxOnly
     procedure, private    :: ClassGeneralBlockUnsetBoxOnly
     procedure, private    :: ClassGeneralBlockUnsetBraBoxOnly
     procedure, private    :: ClassGeneralBlockUnsetKetBoxOnly
     procedure, private    :: ClassGeneralBlockSetAsympBsPermutation
     procedure, private    :: ClassGeneralBlockUnsetAsympBsPermutation
     procedure, private    :: ClassGeneralBlockIsBoxOnly
     procedure, private    :: ClassGeneralBlockIsBraBoxOnly
     procedure, private    :: ClassGeneralBlockIsKetBoxOnly
     procedure, private    :: ClassGeneralBlockSetFormattedWrite
     procedure, private    :: ClassGeneralBlockUnsetFormattedWrite
     procedure, private    :: ClassGeneralBlockGetNumBsNotOvDiffBra
     procedure, private    :: ClassGeneralBlockGetNumBsNotOvDiffKet
     procedure, private    :: ClassGeneralBlockGetNRows
     procedure, private    :: ClassGeneralBlockGetNColumns
     procedure, private    :: ClassGeneralBlockFree
     final                 :: ClassGeneralBlockFinal
  end type ClassGeneralBlock


  type, private, extends( ClassGeneralBlock ) :: LocLocBlock
     type(ClassSymmetricLocalizedStates),  pointer :: SpaceBra
     type(ClassSymmetricLocalizedStates),  pointer :: SpaceKet
   contains
     generic, public :: Build    => LocLocBlockBuild
     generic, public :: BuildCAP => LocLocBlockBuildCAP
     generic, public :: Save     => LocLocBlockSave
     generic, public :: Load     => LocLocBlockLoad
     !
     procedure, private :: LocLocBlockBuild
     procedure, private :: LocLocBlockBuildCAP
     procedure, private :: LocLocBlockSave
     procedure, private :: LocLocBlockLoad
  end type LocLocBlock
  

  type, private, extends( ClassGeneralBlock ) :: LocSRCBlock
     type(ClassSymmetricLocalizedStates),  pointer :: SpaceBra
     type(ClassShortRangeChannel)       ,  pointer :: SpaceKet
   contains
     generic, public :: Build    => LocSRCBlockBuild
     generic, public :: BuildCAP => LocSRCBlockBuildCAP
     generic, public :: Save     => LocSRCBlockSave
     generic, public :: Load     => LocSRCBlockLoad
     !
     procedure, private :: LocSRCBlockBuild
     procedure, private :: LocSRCBlockBuildCAP
     procedure, private :: LocSRCBlockSave
     procedure, private :: LocSRCBlockLoad
  end type LocSRCBlock


  type, private, extends( ClassGeneralBlock ) :: LocPWCBlock
     type(ClassSymmetricLocalizedStates),  pointer :: SpaceBra
     type(ClassPartialWaveChannel)      ,  pointer :: SpaceKet
   contains
     generic, public :: Build      => LocPWCBlockBuild
     generic, public :: Save       => LocPWCBlockSave
     generic, public :: Load       => LocPWCBlockLoad
     generic, public :: Condition  => LocPWCBlockCondition
     !
     procedure, private :: LocPWCBlockBuild
     procedure, private :: LocPWCBlockSave
     procedure, private :: LocPWCBlockLoad
     procedure, private :: LocPWCBlockCondition
  end type LocPWCBlock


  type, private, extends( ClassGeneralBlock ) :: SRCLocBlock
     type(ClassShortRangeChannel)       ,  pointer :: SpaceBra
     type(ClassSymmetricLocalizedStates),  pointer :: SpaceKet
   contains
     generic, public :: Build    => SRCLocBlockBuild
     generic, public :: BuildCAP => SRCLocBlockBuildCAP
     generic, public :: Save     => SRCLocBlockSave
     generic, public :: Load     => SRCLocBlockLoad
     !
     procedure, private :: SRCLocBlockBuild
     procedure, private :: SRCLocBlockBuildCAP
     procedure, private :: SRCLocBlockSave
     procedure, private :: SRCLocBlockLoad
  end type SRCLocBlock


  type, private, extends( ClassGeneralBlock ) :: SRCSRCBlock
     type(ClassShortRangeChannel)       ,  pointer :: SpaceBra
     type(ClassShortRangeChannel)       ,  pointer :: SpaceKet
   contains
     generic, public :: Build    => SRCSRCBlockBuild
     generic, public :: BuildCAP => SRCSRCBlockBuildCAP
     generic, public :: Save     => SRCSRCBlockSave
     generic, public :: Load     => SRCSRCBlockLoad
     !
     procedure, private :: SRCSRCBlockBuild
     procedure, private :: SRCSRCBlockBuildCAP
     procedure, private :: SRCSRCBlockSave
     procedure, private :: SRCSRCBlockLoad
  end type SRCSRCBlock


  type, private, extends( ClassGeneralBlock ) :: SRCPWCBlock
     type(ClassShortRangeChannel)       ,  pointer :: SpaceBra
     type(ClassPartialWaveChannel)      ,  pointer :: SpaceKet
   contains
     generic, public  :: Build            => SRCPWCBlockBuild
     generic, public  :: Save             => SRCPWCBlockSave
     generic, public  :: Load             => SRCPWCBlockLoad
     generic, public  :: Condition        => SRCPWCBlockCondition
     generic, private :: BuildOverlap     => SRCPWCBlockBuildOverlap
     generic, private :: BuildHamiltonian => SRCPWCBlockBuildHamiltonian
     generic, private :: BuildKinEnergy   => SRCPWCBlockBuildKinEnergy
     generic, private :: BuildDipole      => SRCPWCBlockBuildDipole
     !
     procedure, private :: SRCPWCBlockBuild
     procedure, private :: SRCPWCBlockSave
     procedure, private :: SRCPWCBlockLoad
     procedure, private :: SRCPWCBlockCondition
     procedure, private :: SRCPWCBlockBuildOverlap
     procedure, private :: SRCPWCBlockBuildHamiltonian
     procedure, private :: SRCPWCBlockBuildKinEnergy
     procedure, private :: SRCPWCBlockBuildDipole
  end type SRCPWCBlock
  

  type, private, extends( ClassGeneralBlock ) :: PWCLocBlock
     type(ClassPartialWaveChannel)      ,  pointer :: SpaceBra
     type(ClassSymmetricLocalizedStates),  pointer :: SpaceKet
   contains
     generic, public :: Build      => PWCLocBlockBuild
     generic, public :: Save       => PWCLocBlockSave
     generic, public :: Load       => PWCLocBlockLoad
     generic, public :: Condition  => PWCLocBlockCondition
     !
     procedure, private :: PWCLocBlockBuild
     procedure, private :: PWCLocBlockSave
     procedure, private :: PWCLocBlockLoad
     procedure, private :: PWCLocBlockCondition
  end type PWCLocBlock


  type, private, extends( ClassGeneralBlock ) :: PWCSRCBlock
     type(ClassPartialWaveChannel)      ,  pointer :: SpaceBra
     type(ClassShortRangeChannel)       ,  pointer :: SpaceKet
   contains
     generic, public  :: Build              => PWCSRCBlockBuild
     generic, public  :: Save               => PWCSRCBlockSave
     generic, public  :: Load               => PWCSRCBlockLoad
     generic, public  :: Condition          => PWCSRCBlockCondition
     generic, private :: BuildOverlap       => PWCSRCBlockBuildOverlap
     generic, private :: BuildHamiltonian   => PWCSRCBlockBuildHamiltonian
     generic, private :: BuildKinEnergy     => PWCSRCBlockBuildKinEnergy
     generic, private :: BuildDipole        => PWCSRCBlockBuildDipole
     !
     procedure, private :: PWCSRCBlockBuild
     procedure, private :: PWCSRCBlockSave
     procedure, private :: PWCSRCBlockLoad
     procedure, private :: PWCSRCBlockCondition
     procedure, private :: PWCSRCBlockBuildOverlap
     procedure, private :: PWCSRCBlockBuildHamiltonian
     procedure, private :: PWCSRCBlockBuildKinEnergy
     procedure, private :: PWCSRCBlockBuildDipole
  end type PWCSRCBlock
  
  
  type, private, extends( ClassGeneralBlock ) :: PWCPWCBlock
     type(ClassPartialWaveChannel)      ,  pointer :: SpaceBra
     type(ClassPartialWaveChannel)      ,  pointer :: SpaceKet
   contains
     generic, public  :: Build            => PWCPWCBlockBuild
     generic, public  :: Save             => PWCPWCBlockSave
     generic, public  :: Load             => PWCPWCBlockLoad
     generic, public  :: Condition        => PWCPWCBlockCondition
     generic, private :: BuildOverlap     => PWCPWCBlockBuildOverlap
     generic, private :: BuildHamiltonian => PWCPWCBlockBuildHamiltonian
     generic, private :: BuildKinEnergy   => PWCPWCBlockBuildKinEnergy
     generic, private :: BuildDipole      => PWCPWCBlockBuildDipole
     generic, private :: BuildCAP         => PWCPWCBlockBuildCAP
     !
     procedure, private :: PWCPWCBlockBuild
     procedure, private :: PWCPWCBlockSave
     procedure, private :: PWCPWCBlockLoad
     procedure, private :: PWCPWCBlockCondition
     procedure, private :: PWCPWCBlockBuildOverlap
     procedure, private :: PWCPWCBlockBuildHamiltonian
     procedure, private :: PWCPWCBlockBuildKinEnergy
     procedure, private :: PWCPWCBlockBuildDipole
     procedure, private :: PWCPWCBlockBuildCAP
  end type PWCPWCBlock


  

  public :: GetOperatorXlm
  public :: GetOperatorIrrep
  public :: GetFullDipoleOperatorLabel
  public :: GetFullCAPLabel
  public :: GetFullQHLabel
  public :: GetFullFeshbachHLabel




contains


  !> Obtain a new label with the Xlm identifier added.
  function GetFullDipoleOperatorLabel( DipOperatorLabel, DipoleAxis ) result(Label)
    character(len=*), intent(in) :: DipOperatorLabel
    character(len=*), intent(in) :: DipoleAxis
    character(len=:), allocatable :: Label
    type(ClassXlm) :: OpXlm
    integer :: l, m
    character(len=:), allocatable :: Axis
    !.. Apparently extra work but gives flexibility to GetOperatorXlm function.
    allocate( Axis, source = DipoleAxis )
    OpXlm = GetOperatorXlm( DipOperatorLabel, Axis )
    l = OpXlm.GetL()
    m = OpXlm.GetM()
    allocate( Label, source = &
         DipOperatorLabel//'_'//&
         AlphabeticNumber(l)//'.'//&
         AlphabeticNumber(m) )
    !
  end function GetFullDipoleOperatorLabel


  !> Obtain the full CAP label.
  function GetFullCAPLabel( iCAP ) result(Label)
    integer, intent(in) :: iCAP
    character(len=:), allocatable :: Label
    allocate( Label, source = CAPLabel//AlphabeticNumber(iCAP)//"_" )
  end function GetFullCAPLabel


  !> Obtain the full QH label.
  function GetFullQHLabel( iCAP ) result(Label)
    integer, intent(in) :: iCAP
    character(len=:), allocatable :: Label
    allocate( Label, source = QHLabel//AlphabeticNumber(iCAP)//"_" )
  end function GetFullQHLabel


  !> Obtain the full Feshbach partitioning Hamiltonian label.
  function GetFullFeshbachHLabel( ) result(Label)
    character(len=:), allocatable :: Label
    allocate( Label, source = FeshbachHLabel//"_" )
  end function GetFullFeshbachHLabel



  !> Gets the number of B-splines that do not overlap with the diffuse orbitals.
  integer function GetNumBsNotOverlappingDiff( PWC ) result(NBs)
    class(ClassPartialWaveChannel), target, intent(in) :: PWC
    !
    type(ClassIrrep), pointer :: DiffIrrep
    type(ClassXlm) :: Xlm, OverlapXlm
    character(len=:), allocatable :: StorageDir, Axis
    type(ClassDiffuseBsplineXlmBlock)    :: DBSBlock
    type(ClassBsplineXlmBsplineXlmBlock) :: BSBSBlock
    type(ClassMatrix) :: DBSMatrix, BSBSMatrix
    !
    allocate( StorageDir, source = GetStorageDir(PWC) )
    call Xlm.Init( PWC.GetL(), PWC.GetM() )
    OverlapXlm = GetOperatorXlm(OverlapLabel,Axis)
    allocate( DiffIrrep, source = GetDiffuseIrrep(PWC,OverlapXlm) )
    !
    call DBSBlock.init( DiffIrrep, OverlapXlm, OverlapLabel, Xlm )
    call DBSBlock.ReadBlock( StorageDir, DBSMatrix )
    call DBSBlock.Free()
    !
    call BSBSBlock.init( Xlm, OverlapXlm, OverlapLabel, Xlm )
    call BSBSBlock.ReadBlock( StorageDir, BSBSMatrix )
    call BSBSBlock.Free()
    !
    NBs = GetNumBsNotOverlappingDiffuse( DBSMatrix, BSBSMatrix )
    !
  end function GetNumBsNotOverlappingDiff




  !> Giving a SRC-PWC, PWC-SRC or PWC-PWC block, put this block in the correct position, 
  !! provided that indeed the SRC part correspond to the diffuse, and not to the localized
  !! polycentric, whose corresponding sub-block is not touch at this stage (considered zero
  !! as the initialization).
  subroutine PutBlockInRigthPosition( RefMat, InBlock, OutMat )
    !> Matrix containing the correct full dimensions that serves as a reference.
    class(ClassMatrix), intent(in)  :: RefMat
    !> Block to be places in the right position.
    class(ClassMatrix), intent(in)  :: InBlock
    !> matrix with the full correct dimension containing in the right position the previous block.
    class(ClassMatrix), intent(out) :: OutMat
    !
    integer :: NTotRows, NTotCols, NParRows, NParCols
    !
    NTotRows = RefMat.NRows()
    NTotCols = RefMat.NColumns()
    NParRows = InBlock.NRows()
    NParCols = InBlock.NColumns()
    !
    call OutMat.InitFull(NTotRows,NTotCols)
    call OutMat.ComposeFromBlocks(          &
         NTotRows - NParRows + 1, NTotRows, &
         NTotCols - NParCols + 1, NTotCols, &
         InBlock                 )
    !
  end subroutine PutBlockInRigthPosition


  subroutine CoupleSpins( SpaceBra, SpaceKet, MatBlock )
    class(*), target  , intent(in)    :: SpaceBra
    class(*), target  , intent(in)    :: SpaceKet
    class(Classmatrix), intent(inout) :: MatBlock
    !
    real(kind(1d0)) :: ClebshGordonCoeffBra, ClebshGordonCoeffKet, SignFactor
    !
    ! QC phase for the diffuse part.
    SignFactor           = GetQCPhase( SpaceBra, SpaceKet )
    !
    ClebshGordonCoeffBra = GetClebshGordonCoeff(SpaceBra)  
    ClebshGordonCoeffKet = GetClebshGordonCoeff(SpaceKet)  
    !
    call MatBlock.Multiply( ClebshGordonCoeffBra * ClebshGordonCoeffKet * SignFactor )
    !
  end subroutine CoupleSpins




  !> Gets the Xlm associated with the operator, assuming that 
  !! everything which is not S, H, K or CAP is a dipole.
  function GetOperatorXlm( OpLabel, Axis ) result (OpXlm)
    character(len=*)             , intent(in)    :: OpLabel
    character(len=:), allocatable, intent(inout) :: Axis
    type(ClassXlm) :: OpXlm
    if ( (OpLabel .is. OverlapLabel) .or. &
         (OpLabel .is. HamiltonianLabel) .or. &
         (OpLabel .is. KinEnergyLabel) .or. &
         (OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel)) then
       call OpXlm.Init(0,0)
    else
       if ( allocated(Axis) ) then
          if ( Axis .is. 'X' ) then
             call OpXlm.Init(1,1)
          elseif ( Axis .is. 'Y' ) then
             call OpXlm.Init(1,-1)
          elseif ( Axis .is. 'Z' ) then
             call OpXlm.Init(1,0)
          else
             call Assert( 'Invalid orientation to get the dipole operator Xlm, it must be either "x", "y" or "z".' )
          end if
       else
          call Assert( 'Missing orientation for the dipole in order o get the operator Xlm.' )
       end if
    end if
  end function GetOperatorXlm


  function GetOperatorIrrep( Group, OpLabel, Axis ) result (OpIrrep)
    class(ClassGroup)            , intent(in)    :: Group
    character(len=*)             , intent(in)    :: OpLabel
    character(len=:), allocatable, intent(inout) :: Axis
    type(ClassIrrep) :: OpIrrep
    type(ClassXlm) :: OpXlm
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    OpIrrep = OpXlm.GetIrrep(Group)
  end function GetOperatorIrrep
  

  real(kind(1d0)) function PionElectronCGC( TotSpin, ElecProj, SpinsAddUp ) result( CGC )
    !
    real(kind(1d0)), intent(in) :: TotSpin
    real(kind(1d0)), intent(in) :: ElecProj
    logical,         intent(in) :: SpinsAddUp
    !
    if ( SpinsAddUp ) then
       if ( ElecProj > 0.d0 ) then
          CGC = 1.d0
       else
          CGC = 0.d0
       end if
    else
       if ( ElecProj > 0.d0 ) then
          CGC = -sqrt(0.5d0/(TotSpin+1.d0))
       else
          CGC = sqrt((2.d0*TotSpin+1.d0)/(2.d0*TotSpin+2.d0))
       end if
    end if
    !
  end function PionElectronCGC



  !> Given the total multiplicity, and the parent ion multiplicity, compute the total spin \f$S_{T}\f$, the total spin projection \f$\Sigma_{T}$ (as a convention in the QC program \f$\Sigma_{\alpha}=S_{T}\f$), the parent ion spin \f$S_{\alpha}\f$ and the parent ion spin projection \f$\Sigma_{\alpha}\f$ (as a convention \f$\Sigma_{\alpha}\f$ takes the highest allowed value). Then is returned the Clebsh-Gordan coefficient: 
  !! \f[
  !! $C_{S_{\alpha}\Sigma_{\alpha},\frac{1}{2}\sigma}^{S_{T}\Sigma_{T}}$
  !!
  !! \f]
  !!
  !! Where \f$\sigma\f$ stands for the electron spin projection (\f$\Sigma_{T}=\Sigma_{\alpha}+\Sigma_{\alpha}\f$) The formulas used are those in "Quantum Theory of Angular Momentum", Varshalovich et al page 271.
  subroutine ComputeClebshGordonCoeff( &
       TotMultiplicity, &
       PIMultiplicity, &
       ClebshGordonCoeff )
    !
    integer,         intent(in)  :: TotMultiplicity
    integer,         intent(in)  :: PIMultiplicity
    real(kind(1d0)), intent(out) :: ClebshGordonCoeff
    !
    real(kind(1d0)) :: TotSpin, TotProj, PISpin, PIProj, ElecProj, PISpinAux
    real(kind(1d0)), allocatable :: AllowedPIProj(:)
    integer :: i
    real(kind(1d0)), parameter :: Tol = 1.d-10
    logical :: SpinsAddUp
    !
    TotSpin = 0.5d0 * (dble(TotMultiplicity) -1.d0)
    TotProj = TotSpin
    !
    PISpin = 0.5d0 * (dble(PIMultiplicity) -1.d0)
    allocate( AllowedPIProj(PIMultiplicity) )
    !
    ! From higher to lower.
    PISpinAux = PISpin
    do i = 1, PIMultiplicity
       AllowedPIProj(i) = PISpinAux
       PISpinAux = PISpinAux - 1
    end do
    !
    ! Select the maximum posible PI spin projection
    do i = 1, PIMultiplicity
       !
       if ( abs(TotProj-(AllowedPIProj(i)+0.5d0))<Tol ) then
          PIProj   = AllowedPIProj(i)
          ElecProj = 0.5d0
          exit
       elseif ( abs(TotProj-(AllowedPIProj(i)-0.5d0))<Tol ) then
          PIProj   = AllowedPIProj(i)
          ElecProj = -0.5d0
          exit
       end if
       !
    end do
    !
    if ( abs(TotSpin-(PISpin+0.5d0))<Tol ) then
       SpinsAddUp = .true.
    elseif ( abs(TotSpin-(PISpin-0.5d0))<Tol ) then
       SpinsAddUp = .false.
    else
       call Assert( "The spin addition is not working" )
    end if
    !
    ClebshGordonCoeff = PionElectronCGC( TotSpin, ElecProj, SpinsAddUp )
    !
  end subroutine ComputeClebshGordonCoeff



  subroutine ComputeMultipoleContribution( SpaceBra, SpaceKet, MultMat )
    class(*)        , target, intent(in)    :: SpaceBra
    class(*)        , target, intent(in)    :: SpaceKet
    class(ClassMatrix)      , intent(inout) :: MultMat
    !
    integer :: l, Lmax, m, NElectPI
    real(kind(1d0)) :: Factor, FactorL
    type(ClassMultipole):: MultipolePI
    real(kind(1d0)) :: PIXlmMultipole
    type(ClassMatrix) :: MonoElectMultipole, AuxMultMat
    real(kind(1d0)), parameter :: TolMultPI = 1.d-10
    !
    NElectPI = GetNumberElectPI(SpaceBra)
    Factor = 4.d0 * PI * dble(NElectPI)
    !
    MultipolePI = GetPIMultipole(SpaceBra,SpaceKet)
    Lmax = GetLmax( SpaceBra )
    !
    call MonoElectMultipole.Free()
    !
    do l = 1, 2*Lmax
       FactorL  = 1.d0/dble(2*l+1)
       do m = -l, l
          !
          PIXlmMultipole = MultipolePI.GetMultipole(l,m)
          if ( abs(PIXlmMultipole) < TolMultPI ) cycle
          MonoElectMultipole = GetMonoElectMultipoleMat( SpaceBra, SpaceKet, l, m )           
          !
          if ( .not.MonoElectMultipole.IsInitialized() ) cycle
          !
          call MonoElectMultipole.Multiply( PIXlmMultipole * FactorL )
          !..Put the block in the right position 
          call PutBlockInRigthPosition( MultMat, MonoElectMultipole, AuxMultMat )
          call MonoElectMultipole.Free()
          !
          call MultMat.Add( AuxMultMat )
          call AuxMultMat.Free()
          !
          !
       end do
    end do
    !
    call MultMat.Multiply( Factor )
    !
  end subroutine ComputeMultipoleContribution



  function GetPIsDipole( SpaceBra, SpaceKet, OpXlm ) result(PIsDip)
    class(*)        , target, intent(in) :: SpaceBra
    class(*)        , target, intent(in) :: SpaceKet
    class(ClassXlm)         , intent(in) :: OpXlm
    real(kind(1d0)) :: PIsDip
    !
    type(ClassMultipole) :: MultipolePI
    character(len=:), allocatable :: DipoleOrientation
    !
    MultipolePI = GetPIMultipole(SpaceBra,SpaceKet)
    allocate( DipoleOrientation, source = GetDipoleOrientation(OpXlm) )
    PIsDip = MultipolePI.FetchCartesianMultipole( DipoleOrientation )
    !
  end function GetPIsDipole


  function GetPIsKinEnergy( SpaceBra, SpaceKet ) result(PIsKinEnergy)
    class(*)        , target, intent(in) :: SpaceBra
    class(*)        , target, intent(in) :: SpaceKet
    real(kind(1d0)) :: PIsKinEnergy
    !
    type(ClassKineticEnergy) :: PIKinEnerg
    type(ClassParentIon), pointer :: BraPI, KetPI
    !
    allocate( BraPI, source = GetPI(SpaceBra) )
    allocate( KetPI, source = GetPI(SpaceKet) )
    call PIKinEnerg.Init( BraPI, KetPI )
    call PIKinEnerg.Read( GetStorageDir(SpaceBra) )
    PIsKinEnergy = PIKinEnerg.GetKinEnergy()
    !
  end function GetPIsKinEnergy



  function GetDipoleOrientation( OpXlm ) result(Label)
    class(ClassXlm), intent(in) :: OpXlm
    character(len=:), allocatable :: Label
    integer :: l, m
    !
    l = OpXlm.GetL()
    m = OpXlm.GetM()
    !
    if ( l==1 .and. m==1 ) then
       allocate( Label, source = 'x' )
    elseif ( l==1 .and. m==-1 ) then
       allocate( Label, source = 'y' )
    elseif ( l==1 .and. m==0 ) then
       allocate( Label, source = 'z' )
    else
       call Assert( 'Invalid Xlm associated to the dipole operator.' )
    end if
    !
  end function GetDipoleOrientation



  subroutine GetCloseCouplingQCFileName( SpaceBra, SpaceKet, OpLabel, OpXlm, FileExtension, FileName )
    class(*)        , target      , intent(in)   :: SpaceBra
    class(*)        , target      , intent(in)   :: SpaceKet
    character(len=*)              , intent(in)   :: OpLabel
    class(ClassXlm) , target      , intent(in)   :: OpXlm
    character(len=*)              , intent(in)   :: FileExtension
    character(len=:), allocatable , intent(out)  :: FileName 
    !
    character(len=:), allocatable :: BraSymLabel, KetSymLabel
    character(len=:), allocatable :: BraPILabel , KetPILabel
    character(len=:), allocatable :: StorageDir , CCDir
    integer :: lOp, mOp
    !
    allocate( CCDir, source = GetCloseCouplingDir() )
    !
    lOp = OpXlm.GetL()
    mOp = OpXlm.GetM()
    !
    allocate( StorageDir, source = GetStorageDir(SpaceBra) )
    !
    allocate( BraSymLabel, source = GetSymLabel(SpaceBra) )
    allocate( KetSymLabel, source = GetSymLabel(SpaceKet) )
    !
    allocate( BraPILabel, source = GetPILabel(SpaceBra) )
    allocate( KetPILabel, source = GetPILabel(SpaceKet) )
    !
    allocate( FileName, source = &
         AddSlash(StorageDir)//&
         AddSlash(CCDir)//&
         BraSymLabel//'_'//AddSlash(KetSymLabel)//&
         BraPILabel//'_'//AddSlash(KetPILabel)//&
         OpLabel//'_'//AlphabeticNumber(lOp)//'.'//AlphabeticNumber(mOp)//FileExtension )
    !
  end subroutine GetCloseCouplingQCFileName



  subroutine GetDiffuseFileName( InSpace, OpLabel, OpXlm, FileName )
    class(*)        , target      , intent(in)   :: InSpace
    character(len=*)              , intent(in)   :: OpLabel
    class(ClassXlm) , target      , intent(in)   :: OpXlm
    character(len=:), allocatable , intent(out)  :: FileName 
    !
    character(len=:), allocatable :: DiffuseSymLabel
    character(len=:), allocatable :: StorageDir , DiffuseDir
    integer :: lOp, mOp, l, m
    !
    allocate( DiffuseDir, source = GetDiffuseOrbDir() )
    !
    lOp = OpXlm.GetL()
    mOp = OpXlm.GetM()
    l   = OuterElectronL( InSpace )
    m   = OuterElectronM( InSpace )
    !
    allocate( StorageDir, source = GetStorageDir(InSpace) )
    allocate( DiffuseSymLabel, source = GetDiffuseSymLabel(InSpace,OpXlm) )
    !
    allocate( FileName, source = &
         AddSlash(StorageDir)//&
         AddSlash(DiffuseDir)//&
         OpLabel//AlphabeticNumber(lOp)//'.'//AlphabeticNumber(mOp)//&
         '_'//DiffuseSymLabel//'_'//&
         AlphabeticNumber(l)//'.'//AlphabeticNumber(m) )
    !
  end subroutine GetDiffuseFileName



  subroutine GetBsplinesFileName( InSpaceBra, InSpaceKet, OpLabel, OpXlm, FileName )
    class(*)        , target      , intent(in)   :: InSpaceBra
    class(*)        , target      , intent(in)   :: InSpaceKet
    character(len=*)              , intent(in)   :: OpLabel
    class(ClassXlm) , target      , intent(in)   :: OpXlm
    character(len=:), allocatable , intent(out)  :: FileName 
    !
    character(len=:), allocatable :: StorageDir , BsDir
    integer :: lOp, mOp, lBra, mBra, lKet, mKet
    !
    allocate( BsDir, source = GetBsplineOrbDir() )
    !
    lOp = OpXlm.GetL()
    mOp = OpXlm.GetM()
    lBra   = OuterElectronL( InSpaceBra )
    mBra   = OuterElectronM( InSpaceBra )
    lKet   = OuterElectronL( InSpaceKet )
    mKet   = OuterElectronM( InSpaceKet )
    !
    allocate( StorageDir, source = GetStorageDir(InSpaceBra) )
    !
    allocate( FileName, source = &
         AddSlash(StorageDir)//&
         AddSlash(BsDir)//&
         OpLabel//AlphabeticNumber(lOp)//'.'//AlphabeticNumber(mOp)//&
         '_'//AlphabeticNumber(lBra)//'.'//AlphabeticNumber(mBra)//&
         '_'//AlphabeticNumber(lKet)//'.'//AlphabeticNumber(mKet) )
    !
  end subroutine GetBsplinesFileName


  
  function GetStorageDir( InSpace ) result(StorageDir)
    class(*), target, intent(in)  :: InSpace
    character(len=:), allocatable :: StorageDir
    type(ClassSymmetricLocalizedStates), pointer :: Loc
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassSymmetricLocalizedStates )
       allocate( Loc, source = inSpace )
       allocate( StorageDir, source = Loc.GetStorageDir() )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = inSpace )
       allocate( StorageDir, source = SRC.GetStorageDir() )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = inSpace )
       allocate( StorageDir, source = PWC.GetStorageDir() )
    class default 
       call Assert( "Invalid electronic subspace." )
    end select
       !
  end function GetStorageDir


  function GetDiffuseSymLabel( InSpace, OpXlm ) result(DiffuseSymLabel)
    class(*), target, intent(in)  :: InSpace
    class(ClassXlm) , intent(in)  :: OpXlm
    character(len=:), allocatable :: DiffuseSymLabel
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    type(ClassIrrep), pointer :: OpIrrep, OrbIrrep, ResIrrep
    type(ClassGroup), pointer :: Group
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       allocate( DiffuseSymLabel, source = SRC.GetDiffuseIrrepName() )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       OrbIrrep => PWC.GetOrbIrrep()
       Group => OrbIrrep.GetGroup()
       OpIrrep => OpXlm.GetIrrep(Group)
       ResIrrep => OpIrrep * OrbIrrep
       allocate( DiffuseSymLabel, source = ResIrrep.GetName() )
    class default 
       call Assert( "Invalid electronic subspace to get the diffuse irrep name." )
    end select
       !
  end function GetDiffuseSymLabel


  function GetDiffuseIrrep( InSpace, OpXlm ) result(DiffIrrep)
    class(*), target, intent(in)  :: InSpace
    class(ClassXlm) , intent(in)  :: OpXlm
    type(ClassIrrep), pointer :: DiffIrrep
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    type(ClassIrrep), pointer :: OpIrrep, OrbIrrep, ResIrrep
    type(ClassGroup), pointer :: Group
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       allocate( DiffIrrep, source = SRC.GetDiffuseIrrep() )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       OrbIrrep => PWC.GetOrbIrrep()
       Group => OrbIrrep.GetGroup()
       OpIrrep => OpXlm.GetIrrep(Group)
       ResIrrep => OpIrrep * OrbIrrep
       allocate( DiffIrrep, source = ResIrrep )
    class default 
       call Assert( "Invalid electronic subspace to get the diffuse irrep name." )
    end select
       !
  end function GetDiffuseIrrep


  integer function GetLMax( InSpace ) result(LMax)
    class(*), target, intent(in)  :: InSpace
    type(ClassShortRangeChannel) , pointer :: SRC
    type(ClassPartialWaveChannel), pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       LMax = SRC.GetLMax()
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       LMax = PWC.GetLMax()
    class default 
       call Assert( "Invalid electronic subspace to get the maximum angular momentum." )
    end select
    !
  end function GetLMax


  integer function GetNumPolycentricOrb( InSpace ) result(NumPoly)
    class(*), target, intent(in)  :: InSpace
    type(ClassShortRangeChannel) , pointer :: SRC
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       NumPoly = SRC.GetNumMolOrbitals()
    class default 
       call Assert( "Invalid electronic subspace to get the number of polycentric orbitals." )
    end select
    !
  end function GetNumPolycentricOrb



  integer function GetNumberElectPI( InSpace ) result(NEPI)
    class(*), target, intent(in)  :: InSpace
    type(ClassShortRangeChannel) , pointer :: SRC
    type(ClassPartialWaveChannel), pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       NEPI = SRC.GetPINelect()
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       NEPI = PWC.GetPINelect()
    class default 
       call Assert( "Invalid electronic subspace to get the number of electrons of the parent ion." )
    end select
    !
  end function GetNumberElectPI


  function GetPIMultipole( InSpaceBra, InSpaceKet ) result(PIMultipole)
    class(*), target, intent(in)  :: InSpaceBra
    class(*), target, intent(in)  :: InSpaceKet
    type(ClassMultipole) :: PIMultipole
    type(ClassParentIon), pointer :: BraPI, KetPI
    integer :: Lmax
    character(len=:), allocatable :: StorageDir
    !
    allocate( BraPI, source = GetPI(InSpaceBra) )
    allocate( KetPI, source = GetPI(InSpaceKet) )
    Lmax = GetLmax( InSpaceBra )
    !
    call PIMultipole.Init( 2*Lmax, BraPI, KetPI )
    !
    allocate( StorageDir, source = GetStorageDir(InSpaceBra) )
    call PIMultipole.Read( StorageDir )
    !
  end function GetPIMultipole


  function GetMonoElectMultipoleMat( InSpaceBra, InSpaceKet, l, m ) result(MonMat)
    class(*), target, intent(in)  :: InSpaceBra
    class(*), target, intent(in)  :: InSpaceKet
    integer         , intent(in)  :: l
    integer         , intent(in)  :: m
    type(ClassMatrix)  :: MonMat
    !
    integer :: lKet, mKet, lBra, mBra
    type(ClassShortRangeChannel) , pointer :: SRC
    type(ClassPartialWaveChannel), pointer :: PWC
    character(len=:), allocatable :: StorageDir, FileName, FullFileName
    type(ClassDiffuseBsplineXlmBlock)    :: DiffBsBlock
    type(ClassBsplineXlmBsplineXlmBlock) :: BsBsBlock
    type(ClassXlm ) :: OpXlm, BraXlm, KetXlm
    type(ClassIrrep), pointer :: BraIrrep, KetIrrep
    logical :: TransposeMat, exist
    integer :: uid
    !
    call OpXlm.Init( l, m )
    !
    select type( ptrBra => InSpaceBra )
    class is ( ClassShortRangeChannel )
       !
       allocate( SRC, source = InSpaceBra )
       BraIrrep => SRC.GetTotIrrep()
       select type( ptrKet => inSpaceKet )
       class is ( ClassShortRangeChannel )
          call Assert( 'Invalid bra and ket space to fetch monoelectronic multipole matrix.' )
       class is ( ClassPartialWaveChannel )
          lKet   = OuterElectronL( InSpaceKet )
          mKet   = OuterElectronM( InSpaceKet )
          call KetXlm.Init( lKet, mKet )
          call DiffBsBlock.Init( BraIrrep, OpXlm, MultipoleLabel, KetXlm )
          allocate( FileName, source = DiffBsBlock.GetFile() )
          TransposeMat = .false.
       class default 
          call Assert( "Invalid electronic subspace to get the parent ion class." )
       end select
       !
    class is ( ClassPartialWaveChannel )
       lBra   = OuterElectronL( InSpaceBra )
       mBra   = OuterElectronM( InSpaceBra )
       call BraXlm.Init( lBra, mBra )
       select type( ptrKet => inSpaceKet )
       class is ( ClassShortRangeChannel )
          allocate( SRC, source = InSpaceKet )
          KetIrrep => SRC.GetTotIrrep()
          call DiffBsBlock.Init( KetIrrep, OpXlm, MultipoleLabel, BraXlm )
          allocate( FileName, source = DiffBsBlock.GetFile() )
          TransposeMat = .true.
       class is ( ClassPartialWaveChannel )
          allocate( PWC, source = InSpaceKet )
          lKet   = OuterElectronL( InSpaceKet )
          mKet   = OuterElectronM( InSpaceKet )
          call KetXlm.Init( lKet, mKet )
          call BsBsBlock.Init( BraXlm, OpXlm, MultipoleLabel, KetXlm )
          allocate( FileName, source = BsBsBlock.GetFile() )
          TransposeMat = .false.
       class default 
          call Assert( "Invalid electronic subspace to get the parent ion class." )
       end select
    end select
    !
    !
    allocate( StorageDir, source = GetStorageDir(InSpaceBra) )
    allocate( FullFileName, source = AddSlash(StorageDir)//FileName )
    !
    INQUIRE( file = FullFileName, exist = exist )
    if (.not.exist ) return
    !
    call OpenFile( FullFileName, uid, 'read', 'formatted' )
    call MonMat.read( uid )
    close( uid )
    !
    if ( TransposeMat ) call MonMat.Transpose()
    !
  end function GetMonoElectMultipoleMat
  

  function GetPI ( InSpace ) result(PI)
    class(*), target, intent(in)  :: InSpace
    type(ClassParentIon), pointer :: PI
    type(ClassShortRangeChannel) , pointer :: SRC
    type(ClassPartialWaveChannel), pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       allocate( PI, source = SRC.GetPI() )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       allocate( PI, source = PWC.GetPI() )
    class default 
       call Assert( "Invalid electronic subspace to get the parent ion class." )
    end select
    !
  end function GetPI



  function OuterElectronL( InSpace ) result(l)
    class(*), target, intent(in) :: InSpace
    integer :: l
    type(ClassPartialWaveChannel), pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       l = PWC.GetL()
    class default 
       call Assert( "Invalid electronic subspace to get to asymptotic electron angular momentum." )
    end select
    !
  end function OuterElectronL


  function OuterElectronM( InSpace ) result(m)
    class(*), target, intent(in) :: InSpace
    integer :: m
    type(ClassPartialWaveChannel), pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       m = PWC.GetM()
    class default 
       call Assert( "Invalid electronic subspace to get to asymptotic electron angular momentum projection." )
    end select
    !
  end function OuterElectronM


  function GetAngularMomentumLabel( InSpace ) result(AMlabel)
    class(*), target,  intent(in) :: InSpace
    character(len=:), allocatable :: AMlabel
    type(ClassPartialWaveChannel), pointer :: PWC
    integer :: l, m
    !
    select type( ptr => InSpace )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       l = PWC.GetL()
       m = PWC.GetM()
       allocate( AMlabel, source = AlphabeticNumber(l)//'.'//AlphabeticNumber(m) )
    class default
       allocate( AMlabel, source ='' )
    end select
    !
  end function GetAngularMomentumLabel


  function GetSymLabel( InSpace ) result(SymLabel)
    class(*), target,  intent(in) :: InSpace
    character(len=:), allocatable :: SymLabel
    type(ClassSymmetricLocalizedStates), pointer :: Loc
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassSymmetricLocalizedStates )
       allocate( Loc, source = InSpace )
       allocate( SymLabel, source = Loc.GetSymLabel() )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       allocate( SymLabel, source = SRC.GetSymLabel() )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       allocate( SymLabel, source = PWC.GetSymLabel() )
    class default 
       call Assert( "Invalid electronic subspace." )
    end select
    !
  end function GetSymLabel



  function GetPILabel( InSpace ) result(PILabel)
    class(*), target,  intent(in) :: InSpace
    character(len=:), allocatable :: PILabel
    type(ClassSymmetricLocalizedStates), pointer :: Loc
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassSymmetricLocalizedStates )
       allocate( Loc, source = InSpace )
       allocate( PILabel, source = Loc.GetPILabel() )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       allocate( PILabel, source = SRC.GetPILabel() )
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       allocate( PILabel, source = PWC.GetPILabel() )
    class default 
       call Assert( "Invalid electronic subspace." )
    end select
    !
  end function GetPILabel


  function GetPIMultiplicity( InSpace ) result(PIMult)
    class(*), target,  intent(in) :: InSpace
    integer :: PIMult
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       PIMult = SRC.GetPIMultiplicity()
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       PIMult = PWC.GetPIMultiplicity()
    class default 
       call Assert( "Invalid electronic subspace." )
    end select
    !
  end function GetPIMultiplicity


  function GetQCPhase( InSpaceBra, InSpaceKet ) result(SignF)
    class(*), target,  intent(in) :: InSpaceBra
    class(*), target,  intent(in) :: InSpaceKet
    real(kind(1d0)) :: SignF
    type(ClassShortRangeChannel)       , pointer :: SRC
    integer :: PIMult
    real(kind(1d0)) :: BraSignF, KetSignF
    !
    BraSignF = 1.d0
    KetSignF = 1.d0
    !
    select type( ptr => InSpaceBra )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpaceBra )
       PIMult = SRC.GetPIMultiplicity()
       BraSignF = -1.d0 + 2.d0 * mod(PIMult,2)
    end select
    !
    select type( ptr => InSpaceKet )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpaceKet )
       PIMult = SRC.GetPIMultiplicity()
       KETSignF = -1.d0 + 2.d0 * mod(PIMult,2)
    end select
    !
    SignF = BraSignF * KETSignF
    !
  end function GetQCPhase


  real(kind(1d0)) function GetClebshGordonCoeff(InSpace) result(CGC)
    class(*), target,  intent(in) :: InSpace
    type(ClassShortRangeChannel) , pointer :: SRC
    integer :: PIMultiplicity, TotMultiplicity
    !
    CGC = 1.d0
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       TotMultiplicity = GetTotMultiplicity(SRC)
       PIMultiplicity  = GetPIMultiplicity(SRC)
       call ComputeClebshGordonCoeff( &
            TotMultiplicity, &
            PIMultiplicity , &
            CGC  )
    end select
    !
  end function GetClebshGordonCoeff


  function GetTotMultiplicity( InSpace ) result(TotMult)
    class(*), target,  intent(in) :: InSpace
    integer :: TotMult
    type(ClassShortRangeChannel)       , pointer :: SRC
    type(ClassPartialWaveChannel)      , pointer :: PWC
    !
    select type( ptr => InSpace )
    class is ( ClassShortRangeChannel )
       allocate( SRC, source = InSpace )
       TotMult = SRC.GetTotMultiplicity()
    class is ( ClassPartialWaveChannel )
       allocate( PWC, source = InSpace )
       TotMult = PWC.GetTotMultiplicity()
    class default 
       call Assert( "Invalid electronic subspace." )
    end select
    !
  end function GetTotMultiplicity



  logical function IsQCBlock( SpaceBra, SpaceKet ) result(QCB)
    class(*), target,  intent(in) :: SpaceBra
    class(*), target,  intent(in) :: SpaceKet
    !
    QCB = .true.
    !
    select type( ptr => SpaceBra )
    class is ( ClassPartialWaveChannel )
       QCB =.false.
       return
    end select
    !
    select type( ptr => SpaceKet )
    class is ( ClassPartialWaveChannel )
       QCB =.false.
       return
    end select
    !
  end function IsQCBlock


  
  !------------------------------------------------
  ! Methods for ClassGeneralBlock
  !-----------------------------------------------


  subroutine ClassGeneralBlockFinal( Self )
    type(ClassGeneralBlock) :: Self
    call Self.Free()
  end subroutine ClassGeneralBlockFinal


  subroutine ClassGeneralBlockFree( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    if ( allocated(Self.OpLabel) ) deallocate( Self.OpLabel )
    Self.OpXlm    => NULL()
    call Self.Block.Free()
    Self.BoxStatesOnly      = .false.
    Self.BraBoxStatesOnly   = .false.
    Self.KetBoxStatesOnly   = .false.
    Self.AsympBsPermutation = .false.
    Self.FormattedWrite     = .false.
    if ( allocated(Self.FileExtension) ) deallocate( Self.FileExtension )
    Self.NumBsNotOverlappingDiffuseBra = -1
    Self.NumBsNotOverlappingDiffuseKet = -1
  end subroutine ClassGeneralBlockFree


  integer function ClassGeneralBlockGetNumBsNotOvDiffBra( Self ) result( NBs )
    class(ClassGeneralBlock), intent(in) :: Self
    NBs = Self.NumBsNotOverlappingDiffuseBra
  end function ClassGeneralBlockGetNumBsNotOvDiffBra


  integer function ClassGeneralBlockGetNumBsNotOvDiffKet( Self ) result( NBs )
    class(ClassGeneralBlock), intent(in) :: Self
    NBs = Self.NumBsNotOverlappingDiffuseKet
  end function ClassGeneralBlockGetNumBsNotOvDiffKet


  integer function ClassGeneralBlockGetNRows( Self ) result( N )
    class(ClassGeneralBlock), intent(in) :: Self
    N = Self.Block.NRows()
  end function ClassGeneralBlockGetNRows


  integer function ClassGeneralBlockGetNColumns( Self ) result( N )
    class(ClassGeneralBlock), intent(in) :: Self
    N = Self.Block.NColumns()
  end function ClassGeneralBlockGetNColumns


  subroutine ClassGeneralBlockUnsetBoxOnly( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.BoxStatesOnly = .false.
  end subroutine ClassGeneralBlockUnsetBoxOnly

  subroutine ClassGeneralBlockUnsetBraBoxOnly( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.BraBoxStatesOnly = .false.
  end subroutine ClassGeneralBlockUnsetBraBoxOnly


  subroutine ClassGeneralBlockUnsetKetBoxOnly( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.KetBoxStatesOnly = .false.
  end subroutine ClassGeneralBlockUnsetKetBoxOnly


  subroutine ClassGeneralBlockSetBoxOnly( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.BoxStatesOnly = .true.
  end subroutine ClassGeneralBlockSetBoxOnly


  subroutine ClassGeneralBlockSetBraBoxOnly( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.BraBoxStatesOnly = .true.
  end subroutine ClassGeneralBlockSetBraBoxOnly


  subroutine ClassGeneralBlockSetKetBoxOnly( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.KetBoxStatesOnly = .true.
  end subroutine ClassGeneralBlockSetKetBoxOnly


  subroutine ClassGeneralBlockSetAsympBsPermutation( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.AsympBsPermutation = .true.
  end subroutine ClassGeneralBlockSetAsympBsPermutation


  subroutine ClassGeneralBlockUnsetAsympBsPermutation( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.AsympBsPermutation = .false.
  end subroutine ClassGeneralBlockUnsetAsympBsPermutation


  logical function ClassGeneralBlockIsBoxOnly( Self ) result(IsBox)
    class(ClassGeneralBlock), intent(in) :: Self
    IsBox = Self.BoxStatesOnly
  end function ClassGeneralBlockIsBoxOnly


  logical function ClassGeneralBlockIsBraBoxOnly( Self ) result(IsBox)
    class(ClassGeneralBlock), intent(in) :: Self
    IsBox = Self.BraBoxStatesOnly
  end function ClassGeneralBlockIsBraBoxOnly


  logical function ClassGeneralBlockIsKetBoxOnly( Self ) result(IsBox)
    class(ClassGeneralBlock), intent(in) :: Self
    IsBox = Self.KetBoxStatesOnly
  end function ClassGeneralBlockIsKetBoxOnly


  subroutine ClassGeneralBlockUnsetFormattedWrite( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.FormattedWrite = .false.
  end subroutine ClassGeneralBlockUnsetFormattedWrite


  subroutine ClassGeneralBlockSetFormattedWrite( Self )
    class(ClassGeneralBlock), intent(inout) :: Self
    Self.FormattedWrite = .true.
  end subroutine ClassGeneralBlockSetFormattedWrite



  subroutine ClassGeneralBlockReadBlock( Self, FileName )
    class(ClassGeneralBlock), intent(inout) :: Self
    character(len=*)        , intent(in)    :: FileName
    integer :: uid
    !
    if ( Self.FormattedWrite ) then
       call OpenFile( FileName, uid, 'read', 'formatted' )
    else
       call OpenFile( FileName, uid, 'read', 'unformatted' )
    end if
    call Self.Block.Read( uid )
    close( uid )
    !
  end subroutine ClassGeneralBlockReadBlock


  subroutine ClassGeneralBlockSaveBlock( Self, FileName )
    class(ClassGeneralBlock), intent(in) :: Self
    character(len=*)        , intent(in) :: FileName
    !
    if ( Self.FormattedWrite ) then
       call Self.Block.Write( FileName, 'formatted' )
    else
       call Self.Block.Write( FileName, 'unformatted' )
    end if
    !
  end subroutine ClassGeneralBlockSaveBlock




  function ClassGeneralBlockFetchBlock( Self ) result(NewBlock)
    class(ClassGeneralBlock), intent(in) :: Self
    type(ClassMatrix) :: NewBlock
    NewBlock = Self.Block
  end function ClassGeneralBlockFetchBlock



  logical function ClassGeneralBlockValidBlock( Self, FileName ) result(ValidB)
    class(ClassGeneralBlock), intent(in) :: Self
    character(len=*)        , intent(in) :: FileName
    !
    logical :: exist
    type(ClassMatrix) :: StoredMat
    integer :: uid
    !
    ValidB = .false.
    !
    INQUIRE( file = FileName, exist = exist )
    if ( .not.exist ) return
    !
    if ( .not.Self.Block.IsInitialized() ) return
    !
    if ( Self.FormattedWrite ) then
       call OpenFile( Filename, uid, 'read', 'formatted' )
    else
       call OpenFile( Filename, uid, 'read', 'unformatted' )
    end if
    !
    call StoredMat.Read( uid )
    close( uid )
    !
    if ( (Self.Block.NRows()==StoredMat.NRows()) .and. &
         (Self.Block.NColumns()==StoredMat.NColumns()) ) ValidB = .true.
    !
  end function ClassGeneralBlockValidBlock



  subroutine ClassGeneralBlockInit( Self, OpLabel, OpXlm )
    class(ClassGeneralBlock), intent(inout) :: Self
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    if ( allocated(Self.OpLabel) ) deallocate(Self.OpLabel)
    allocate( Self.OpLabel, source = OpLabel )
    if ( associated(Self.OpXlm) ) deallocate(Self.OpXlm)
    allocate( Self.OpXlm, source = OpXlm )
  end subroutine ClassGeneralBlockInit


  subroutine ClassGeneralBlockReadQCBlock( Self, SpaceBra, SpaceKet, OpLabel, OpXlm )
    class(ClassGeneralBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: CCFileName
    real(kind(1d0)), allocatable :: ArrayQC(:,:)
    !
    if ( allocated(Self.OpLabel) ) deallocate( Self.OpLabel )
    allocate( Self.OpLabel, source = OpLabel )
    Self.OpXlm    => OpXlm
    !
    call GetCloseCouplingQCFileName( SpaceBra, SpaceKet, OpLabel, OpXlm, Self.FileExtension, CCFileName )
    !
    call ReadMatrix( CCFileName, ArrayQC )
    Self.Block = ArrayQC
    !
  end subroutine ClassGeneralBlockReadQCBlock



  subroutine ClassGeneralBlockSetFileExtension( Self, InName )
    class(ClassGeneralBlock), intent(inout) :: Self
    character(len=*)        , intent(in)    :: InName
    if ( allocated(Self.FileExtension) ) deallocate( Self.FileExtension  )
    allocate( Self.FileExtension, source = InName )
  end subroutine ClassGeneralBlockSetFileExtension


  function ClassGeneralBlockGetBlockName( Self, SpaceBra, SpaceKet ) result(BlockName)
    class(ClassGeneralBlock), intent(in) :: Self
    class(*)  , target      , intent(in) :: SpaceBra
    class(*)  , target      , intent(in) :: SpaceKet
    character(len=:), allocatable        :: BlockName
    !
    character(len=:), allocatable :: BraSymLabel, KetSymLabel, AngMomBraLabel, AngMomKetLabel, AngMomLabel
    character(len=:), allocatable :: BraPILabel , KetPILabel
    character(len=:), allocatable :: StorageDir , CCDir
    integer :: lOp, mOp
    !
    if ( IsQCBlock(SpaceBra,SpaceKet) ) then
       call GetCloseCouplingQCFileName( SpaceBra, SpaceKet, Self.OpLabel, Self.OpXlm, Self.FileExtension, BlockName )
       return
    end if
    !
    allocate( CCDir, source = GetCloseCouplingDir() )
    !
    lOp = Self.OpXlm.GetL()
    mOp = Self.OpXlm.GetM()
    !
    allocate( StorageDir, source = GetStorageDir(SpaceBra) )
    !
    allocate( BraSymLabel, source = GetSymLabel(SpaceBra) )
    allocate( KetSymLabel, source = GetSymLabel(SpaceKet) )
    !
    allocate( BraPILabel, source = GetPILabel(SpaceBra) )
    allocate( KetPILabel, source = GetPILabel(SpaceKet) )
    !
    allocate( AngMomBraLabel, source = GetAngularMomentumLabel(SpaceBra) )
    allocate( AngMomKetLabel, source = GetAngularMomentumLabel(SpaceKet) )
    if ( len(AngMomBraLabel)>0 .and. len(AngMomKetLabel)>0 ) then
       allocate( AngMomLabel, source = AngMomBraLabel//'_'//AngMomKetLabel )
    else
       allocate( AngMomLabel, source = AngMomBraLabel//AngMomKetLabel )
    end if
    !
    allocate( BlockName, source = &
         AddSlash(StorageDir)//&
         AddSlash(CCDir)//&
         BraSymLabel//'_'//AddSlash(KetSymLabel)//&
         BraPILabel//'_'//AddSlash(KetPILabel)//&
         Self.OpLabel//'_'//AlphabeticNumber(lOp)//'.'//AlphabeticNumber(mOp)//'_'//&
         AngMomLabel//Self.FileExtension )
    !
  end function ClassGeneralBlockGetBlockName




  !--------------------------------------------------
  ! Methods for the specific blocks: 
  ! LocLoc, LocSRC, LocPWC
  ! SRCLoc, SRCSRC, SRCPWC
  ! PWCLoc, PWCSRC, PWCPWC
  !-------------------------------------------------



  !--------------------------------------------------
  ! Methods for the LocLocBlock
  !-------------------------------------------------


  subroutine LocLocBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(LocLocBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
!!$    !.. Do nothing for the blocks obtained from QC.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call Self.BuildCAP( )
!!$    end if
    !
  end subroutine LocLocBlockBuild


  subroutine LocLocBlockBuildCAP( Self )
    class(LocLocBlock)   , intent(inout) :: Self
    !
    character(len=:), allocatable :: SRCSRCQCFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: SRCSRCMat
    !
    !.. Get the dimension of the short range part.
    call GetCloseCouplingQCFileName( Self.SpaceKet, Self.SpaceKet, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NCols = SRCSRCMat.NColumns()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call GetCloseCouplingQCFileName( Self.SpaceBra, Self.SpaceBra, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NRows = SRCSRCMat.NRows()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call Self.Block.InitFull(NRows,NCols)
    !
  end subroutine LocLocBlockBuildCAP



  subroutine LocLocBlockSave( Self )
    class(LocLocBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    !.. Do nothing for the blocks obtained from QC.
    if ( Self.OpLabel(:min(len(Self.OpLabel),len(CAPLabel))) .is. CAPLabel ) then
       allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
       call Self.Save( FileName )
    end if
  end subroutine LocLocBlockSave
  

  subroutine LocLocBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm ) 
    class(LocLocBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: BlockName
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
  end subroutine LocLocBlockLoad



  !--------------------------------------------------
  ! Methods for the LocSCRBlock
  !-------------------------------------------------

  subroutine LocSRCBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(LocSRCBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
!!$    !.. Do nothing for the blocks obtained from QC.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call Self.BuildCAP( )
!!$    end if
    !
  end subroutine LocSRCBlockBuild


  subroutine LocSRCBlockBuildCAP( Self )
    class(LocSRCBlock)   , intent(inout) :: Self
    !
    character(len=:), allocatable :: SRCSRCQCFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: SRCSRCMat
    !
    !.. Get the dimension of the short range part.
    call GetCloseCouplingQCFileName( Self.SpaceKet, Self.SpaceKet, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NCols = SRCSRCMat.NColumns()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call GetCloseCouplingQCFileName( Self.SpaceBra, Self.SpaceBra, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NRows = SRCSRCMat.NRows()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call Self.Block.InitFull(NRows,NCols)
    !
  end subroutine LocSRCBlockBuildCAP


  subroutine LocSRCBlockSave( Self )
    class(LocSRCBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    !.. Do nothing for the blocks obtained from QC.
    if ( Self.OpLabel(:min(len(Self.OpLabel),len(CAPLabel))) .is. CAPLabel ) then
       allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
       call Self.Save( FileName )
    end if
  end subroutine LocSRCBlockSave


  subroutine LocSRCBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, OnlyPoly ) 
    class(LocSRCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical, optional       , intent(in)    :: OnlyPoly
    !
    character(len=:), allocatable :: BlockName
    integer :: NumRows, NumCols
    type(Classmatrix) :: AuxMat
    
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
    !
    if ( present(OnlyPoly) .and. OnlyPoly ) then
       NumRows = Self.GetNRows()
       NumCols = GetNumPolycentricOrb( Self.SpaceKet )
       call Self.Block.GetSubMatrix( 1, NumRows, 1, NumCols, AuxMat )
       call Self.Block.Free()
       Self.Block = AuxMat
       call AuxMat.Free()
    end if
    !
  end subroutine LocSRCBlockLoad


  !--------------------------------------------------
  ! Methods for the LocPWCBlock
  !-------------------------------------------------

  subroutine LocPWCBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(LocPWCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: LocLocQCFile, DiffuseFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: LocLocMat
    type(ClassMatrix) :: DiffuseMat
    type(ClassXlm) :: OverlapXLM
    !
    !
    !.. At the moment this block is zero but will be store always.
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
!!$    Self.NumBsNotOverlappingDiffuseKet = GetNumBsNotOverlappingDiff( Self.SpaceKet )
!!$    !
!!$    !.. Get the dimension of the localized part.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call GetCloseCouplingQCFileName( SpaceBra, SpaceBra, OverlapLabel, OpXlm, GetQCFileExtension(), LocLocQCFile )
!!$    else
!!$       call GetCloseCouplingQCFileName( SpaceBra, SpaceBra, OpLabel, OpXlm, GetQCFileExtension(), LocLocQCFile )
!!$    end if
!!$    call OpenFile( LocLocQCFile, uid, 'read', 'formatted' )
!!$    call LocLocMat.Read( uid )
!!$    close ( uid )
!!$    NRows = LocLocMat.NRows()
!!$    call LocLocMat.Free()
!!$    !
!!$    !.. Get the dimension of the B-splines part.
!!$    call OverlapXLM.Init(0,0)
!!$    call GetBsplinesFileName( SpaceKet, SpaceKet, OverlapLabel, OverlapXlm, DiffuseFile )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
!!$    call DiffuseMat.Read( uid )
!!$    close( uid )
!!$    NCols = DiffuseMat.NColumns()
!!$    call DiffuseMat.Free()
!!$    !
    !***
    NRows = 100
    NCols = 100
    !***
    call Self.Block.InitFull(NRows,NCols)
    !
  end subroutine LocPWCBlockBuild



  subroutine LocPWCBlockSave( Self )
    class(LocPWCBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    call Self.Save( FileName )
  end subroutine LocPWCBlockSave


  subroutine LocPWCBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm ) 
    class(LocPWCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: BlockName
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    Self.NumBsNotOverlappingDiffuseKet = GetNumBsNotOverlappingDiff( Self.SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
  end subroutine LocPWCBlockLoad


  subroutine LocPWCBlockCondition( Self, Conditioner, LocSRC )
    class(LocPWCBlock)          , intent(inout) :: Self
    class(ClassConditionerBlock), intent(inout) :: Conditioner
    class(LocSRCBlock)          , intent(in)    :: LocSRC
    !
    type(ClassIrrep), pointer :: DiffIrrep
    type(ClassXlm) :: Xlm, OverlapXlm
    type(ClassMatrix) :: ConditionMat, LocSRCMat
    real(kind(1d0)), allocatable :: Array(:,:)
    character(len=:), allocatable :: StorageDir, Axis
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    !
    call Xlm.Init( Self.SpaceKet.GetL(), Self.SpaceKet.GetM() )
    call Conditioner.SetXlm( Xlm )
    !
    OverlapXlm = GetOperatorXlm(OverlapLabel,Axis)
    allocate( DiffIrrep, source = GetDiffuseIrrep(Self.SpaceKet,OverlapXlm) )
    call Conditioner.SetIrrep( DiffIrrep )
    !
    call Conditioner.ReadBlock( StorageDir, ConditionMat )
    !
    call Self.Block.Multiply( ConditionMat, 'Right', 'N' )
    call ConditionMat.Free()
    !
    call Conditioner.ReadDiffuseBlock( StorageDir, ConditionMat )
    if ( ConditionMat.IsInitialized() ) then
       call LocSRC.Block.FetchMatrix( Array )
       LocSRCMat = Array(:,size(Array,2)-ConditionMat.NRows()+1:)
       deallocate( Array )
       call LocSRCMat.Multiply( ConditionMat, 'Right', 'N' )
       call ConditionMat.Free()
       call Self.Block.Add( LocSRCMat )
    end if
    !
  end subroutine LocPWCBlockCondition


  !--------------------------------------------------
  ! Methods for the SRCLocBlock
  !-------------------------------------------------


  subroutine SRCLocBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(SRCLocBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
!!$    !.. Do nothing for the blocks obtained from QC.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call Self.BuildCAP( )
!!$    end if
    !
  end subroutine SRCLocBlockBuild


  subroutine SRCLocBlockBuildCAP( Self )
    class(SRCLocBlock)   , intent(inout) :: Self
    !
    character(len=:), allocatable :: SRCSRCQCFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: SRCSRCMat
    !
    !.. Get the dimension of the short range part.
    call GetCloseCouplingQCFileName( Self.SpaceKet, Self.SpaceKet, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NCols = SRCSRCMat.NColumns()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call GetCloseCouplingQCFileName( Self.SpaceBra, Self.SpaceBra, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NRows = SRCSRCMat.NRows()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call Self.Block.InitFull(NRows,NCols)
    !
  end subroutine SRCLocBlockBuildCAP


  subroutine SRCLocBlockSave( Self )
    class(SRCLocBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    !.. Do nothing for the blocks obtained from QC.
!!$    if ( Self.OpLabel(:min(len(Self.OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
!!$       call Self.Save( FileName )
!!$    end if
  end subroutine SRCLocBlockSave


  subroutine SRCLocBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, OnlyPoly ) 
    class(SRCLocBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical, optional       , intent(in)    :: OnlyPoly
    !
    character(len=:), allocatable :: BlockName
    integer :: NumRows, NumCols
    type(Classmatrix) :: AuxMat
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
    !
    if ( present(OnlyPoly) .and. OnlyPoly ) then
       NumRows = GetNumPolycentricOrb( Self.SpaceBra )
       NumCols = Self.GetNColumns()
       call Self.Block.GetSubMatrix( 1, NumRows, 1, NumCols, AuxMat )
       call Self.Block.Free()
       Self.Block = AuxMat
       call AuxMat.Free()
    end if
    !
  end subroutine SRCLocBlockLoad


  !--------------------------------------------------
  ! Methods for the SRCSRCBlock
  !-------------------------------------------------



  subroutine SRCSRCBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(SRCSRCBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
!!$    !.. Do nothing for the blocks obtained from QC.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call Self.BuildCAP( )
!!$    end if
    !
  end subroutine SRCSRCBlockBuild


  subroutine SRCSRCBlockBuildCAP( Self )
    class(SRCSRCBlock)   , intent(inout) :: Self
    !
    character(len=:), allocatable :: SRCSRCQCFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: SRCSRCMat
    !
    !.. Get the dimension of the short range part.
    call GetCloseCouplingQCFileName( Self.SpaceKet, Self.SpaceKet, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NCols = SRCSRCMat.NColumns()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call GetCloseCouplingQCFileName( Self.SpaceBra, Self.SpaceBra, OverlapLabel, Self.OpXlm, GetQCFileExtension(), SRCSRCQCFile )
    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
    call SRCSRCMat.Read( uid )
    close ( uid )
    NRows = SRCSRCMat.NRows()
    call SRCSRCMat.Free()
    deallocate(SRCSRCQCFile)
    !
    call Self.Block.InitFull(NRows,NCols)
    !
  end subroutine SRCSRCBlockBuildCAP



  subroutine SRCSRCBlockSave( Self )
    class(SRCSRCBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    !.. Do nothing for the blocks obtained from QC.
    if ( Self.OpLabel(:min(len(Self.OpLabel),len(CAPLabel))) .is. CAPLabel ) then
       allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
       call Self.Save( FileName )
    end if
  end subroutine SRCSRCBlockSave
     

  subroutine SRCSRCBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, OnlyPoly ) 
    class(SRCSRCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical, optional       , intent(in)    :: OnlyPoly
    !
    character(len=:), allocatable :: BlockName
    integer :: NumRows, NumCols
    type(Classmatrix) :: AuxMat
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
    !
    if ( present(OnlyPoly) .and. OnlyPoly ) then
       NumRows = GetNumPolycentricOrb( Self.SpaceBra )
       NumCols = GetNumPolycentricOrb( Self.SpaceKet )
       call Self.Block.GetSubMatrix( 1, NumRows, 1, NumCols, AuxMat )
       call Self.Block.Free()
       Self.Block = AuxMat
       call AuxMat.Free()
    end if
    !
  end subroutine SRCSRCBlockLoad


  !--------------------------------------------------
  ! Methods for the SRCPWCBlock
  !-------------------------------------------------


  subroutine SRCPWCBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(SRCPWCBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in) :: SpaceBra
    class(*)       , target , intent(in) :: SpaceKet
    character(len=*)        , intent(in) :: OpLabel
    class(ClassXlm), target , intent(in) :: OpXlm
    logical                 , intent(in) :: Force
    !
    character(len=:), allocatable :: SRCSRCQCFile, DiffuseFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: SRCSRCMat
    type(ClassMatrix) :: DiffuseMat
    type(ClassIrrep), pointer :: OpIrrep, BraIrrep, KetIrrep
    type(ClassXlm) :: OverlapXLM
    !
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
!!$    Self.NumBsNotOverlappingDiffuseKet = GetNumBsNotOverlappingDiff( Self.SpaceKet )
!!$    !
!!$    !.. Get the dimension of the short range part.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call GetCloseCouplingQCFileName( SpaceBra, SpaceBra, OverlapLabel, OpXlm, GetQCFileExtension(), SRCSRCQCFile )
!!$    else
!!$       call GetCloseCouplingQCFileName( SpaceBra, SpaceBra, OpLabel, OpXlm, GetQCFileExtension(), SRCSRCQCFile )
!!$    end if
!!$    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
!!$    call SRCSRCMat.Read( uid )
!!$    close ( uid )
!!$    NRows = SRCSRCMat.NRows()
!!$    call SRCSRCMat.Free()
!!$    !
!!$    !.. Get the dimension of the B-splines part.
!!$    call OverlapXLM.Init(0,0)
!!$    call GetBsplinesFileName( SpaceKet, SpaceKet, OverlapLabel, OverlapXlm, DiffuseFile )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
!!$    call DiffuseMat.Read( uid )
!!$    close( uid )
!!$    NCols = DiffuseMat.NColumns()
!!$    call DiffuseMat.Free()
    !
    !***
    NRows = 100
    NCols = 100
    !***
    call Self.Block.InitFull(NRows,NCols)
    !
!!$    !.. Determines if the operator is allowed by symmetry
!!$    BraIrrep => Self.SpaceBra.GetTotIrrep()
!!$    KetIrrep => Self.SpaceKet.GetTotIrrep()
!!$    OpIrrep  => OpXlm.GetIrrep( BraIrrep.GetGroup() )
!!$    if ( .not.ValidIrreps(BraIrrep,OpIrrep,KetIrrep) ) return
!!$    !
!!$    !.. Select the operator
!!$    if ( OpLabel .is. OverlapLabel ) then
!!$       call Self.BuildOverlap( Force )
!!$    elseif ( OpLabel .is. HamiltonianLabel ) then
!!$       call Self.BuildHamiltonian( Force )
!!$    elseif ( OpLabel .is. KinEnergyLabel ) then
!!$       call Self.BuildKinEnergy( Force )
!!$    elseif ( (OpLabel .is. DIPOLE_LABEL//LengthLabel) .or. &
!!$         (OpLabel .is. DIPOLE_LABEL//VelocityLabel) ) then
!!$       call Self.BuildDipole( Force )
!!$    elseif ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       ! It is zero.
!!$       return
!!$    else
!!$       call Assert( 'Invalid operator label: '//OpLabel//'.' )
!!$    end if
!!$    !
!!$    call CoupleSpins( Self.SpaceBra, Self.SpaceKet, Self.Block )
    !
  end subroutine SRCPWCBlockBuild



  subroutine SRCPWCBlockSave( Self )
    class(SRCPWCBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
!!$    allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
!!$    call Self.Save( FileName )
  end subroutine SRCPWCBlockSave


  subroutine SRCPWCBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm ) 
    class(SRCPWCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: BlockName
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    Self.NumBsNotOverlappingDiffuseKet = GetNumBsNotOverlappingDiff( Self.SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
    !
  end subroutine SRCPWCBlockLoad



  subroutine SRCPWCBlockCondition( Self, Conditioner, SRCSRC )
    class(SRCPWCBlock)          , intent(inout) :: Self
    class(ClassConditionerBlock), intent(inout) :: Conditioner
    class(SRCSRCBlock)          , intent(in)    :: SRCSRC
    !
    type(ClassIrrep), pointer :: DiffIrrep
    type(ClassXlm) :: Xlm, OverlapXlm
    type(ClassMatrix) :: ConditionMat, SRCSRCMat
    character(len=:), allocatable :: StorageDir, Axis
    real(kind(1d0)), allocatable :: Array(:,:)
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    !
    call Xlm.Init( Self.SpaceKet.GetL(), Self.SpaceKet.GetM() )
    call Conditioner.SetXlm( Xlm )
    !
    OverlapXlm = GetOperatorXlm(OverlapLabel,Axis)
    allocate( DiffIrrep, source = GetDiffuseIrrep(Self.SpaceKet,OverlapXlm) )
    call Conditioner.SetIrrep( DiffIrrep )
    !
    call Conditioner.ReadBlock( StorageDir, ConditionMat )
    !
    call Self.Block.Multiply( ConditionMat, 'Right', 'N' )
    call ConditionMat.Free()
    !
    call Conditioner.ReadDiffuseBlock( StorageDir, ConditionMat )
    if ( ConditionMat.IsInitialized() ) then
       call SRCSRC.Block.FetchMatrix( Array )
       SRCSRCMat = Array(:,size(Array,2)-ConditionMat.NRows()+1:)
       deallocate( Array )
       call SRCSRCMat.Multiply( ConditionMat, 'Right', 'N' )
       call ConditionMat.Free()
       call Self.Block.Add( SRCSRCMat )
    end if
    !
  end subroutine SRCPWCBlockCondition



  subroutine SRCPWCBlockBuildOverlap( Self, Force )
    class(SRCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, AuxMat
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !
    !.. Get the dimension of the B-splines part.
    call GetDiffuseFileName( Self.SpaceKet, Self.OpLabel, Self.OpXlm, DiffuseFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
    call DiffuseMat.Read( uid )
    close( uid )
    !
    call PutBlockInRigthPosition( Self.Block, DiffuseMat, AuxMat )
    call DiffuseMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
  end subroutine SRCPWCBlockBuildOverlap



  subroutine SRCPWCBlockBuildHamiltonian( Self, Force )
    class(SRCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, DiffuseOverlapFile, BlockFileName, StorageDir
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, DiffuseOverlapMat, AuxMat
    real(kind(1d0)) :: PIEnergy
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    !..Multipole contribution
    call ComputeMultipoleContribution( &
         Self.SpaceBra, Self.SpaceKet, &
         Self.Block                  )
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceKet, Self.OpLabel, Self.OpXlm, DiffuseFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
    call DiffuseMat.Read( uid )
    close( uid )
    !
    !.. Get the overlap operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceKet, OverlapLabel, Self.OpXlm, DiffuseOverlapFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseOverlapFile, uid, 'read', 'formatted' )
    call DiffuseOverlapMat.Read( uid )
    close( uid )
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    PIEnergy = Self.SpaceBra.GetPIEnergy(StorageDir)
    call DiffuseOverlapMat.Multiply( PIEnergy )
    call DiffuseOverlapMat.Add( DiffuseMat )
    call DiffuseMat.Free()
    !
    call PutBlockInRigthPosition( Self.Block, DiffuseOverlapMat, AuxMat )
    call DiffuseOverlapMat.Free()
    call Self.Block.Add( AuxMat )
    call AuxMat.Free()
    !
  end subroutine SRCPWCBlockBuildHamiltonian



  subroutine SRCPWCBlockBuildDipole( Self, Force )
    class(SRCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, DiffuseOverlapFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, DiffuseOverlapMat, AuxMat
    real(kind(1d0))   :: PIsDipole
    type(ClassXlm)    :: OverlapXlm
    logical :: exist
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    call OverlapXlm.Init(0,0)
    !.. Get the overlap operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceKet, OverlapLabel, OverlapXlm, DiffuseOverlapFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseOverlapFile, uid, 'read', 'formatted' )
    call DiffuseOverlapMat.Read( uid )
    close( uid )
    !
    PIsDipole = GetPIsDipole( Self.SpaceBra, Self.SpaceKet, Self.OpXlm )
    call DiffuseOverlapMat.Multiply( PIsDipole )
    call PutBlockInRigthPosition( Self.Block, DiffuseOverlapMat, AuxMat )
    call DiffuseOverlapMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceKet, Self.OpLabel, Self.OpXlm, DiffuseFile )
    INQUIRE( File = DiffuseFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
       call DiffuseMat.Read( uid )
       close( uid )
       !
       call PutBlockInRigthPosition( Self.Block, DiffuseMat, AuxMat )
       call DiffuseMat.Free()
       call Self.Block.Add( AuxMat )
       call AuxMat.Free()
    end if
    !
  end subroutine SRCPWCBlockBuildDipole


  subroutine SRCPWCBlockBuildKinEnergy( Self, Force )
    class(SRCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, DiffuseOverlapFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, DiffuseOverlapMat, AuxMat
    real(kind(1d0))   :: PIsKinEnergy
    type(ClassXlm)    :: OverlapXlm
    logical :: exist
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    call OverlapXlm.Init(0,0)
    !.. Get the overlap operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceKet, OverlapLabel, OverlapXlm, DiffuseOverlapFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseOverlapFile, uid, 'read', 'formatted' )
    call DiffuseOverlapMat.Read( uid )
    close( uid )
    !
    PIsKinEnergy = GetPIsKinEnergy( Self.SpaceBra, Self.SpaceKet )
    call DiffuseOverlapMat.Multiply( PIsKinEnergy )
    call PutBlockInRigthPosition( Self.Block, DiffuseOverlapMat, AuxMat )
    call DiffuseOverlapMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceKet, Self.OpLabel, Self.OpXlm, DiffuseFile )
    INQUIRE( File = DiffuseFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
       call DiffuseMat.Read( uid )
       close( uid )
       !
       call PutBlockInRigthPosition( Self.Block, DiffuseMat, AuxMat )
       call DiffuseMat.Free()
       call Self.Block.Add( AuxMat )
       call AuxMat.Free()
    end if
    !
  end subroutine SRCPWCBlockBuildKinEnergy



  !--------------------------------------------------
  ! Methods for the PWCLocBlock
  !-------------------------------------------------


  subroutine PWCLocBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(PWCLocBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: LocLocQCFile, DiffuseFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: LocLocMat
    type(ClassMatrix) :: DiffuseMat
    type(ClassXlm) :: OverlapXLM
    !
    !
    !.. At the moment this block is zero but will be store always.
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
!!$    Self.NumBsNotOverlappingDiffuseBra = GetNumBsNotOverlappingDiff( Self.SpaceBra )
    !
    !.. Get the dimension of the localized part.
    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
       call GetCloseCouplingQCFileName( SpaceKet, SpaceKet, OverlapLabel, OpXlm, GetQCFileExtension(), LocLocQCFile )
    else
       call GetCloseCouplingQCFileName( SpaceKet, SpaceKet, OpLabel, OpXlm, GetQCFileExtension(), LocLocQCFile )
    end if
!!$    call OpenFile( LocLocQCFile, uid, 'read', 'formatted' )
!!$    call LocLocMat.Read( uid )
!!$    close ( uid )
!!$    NCols = LocLocMat.NColumns()
!!$    call LocLocMat.Free()
    !***
    !***
    NCols = 100
    !***
    !***
    !
    !.. Get the dimension of the B-splines part.
    call OverlapXLM.Init(0,0)
!!$    call GetBsplinesFileName( SpaceBra, SpaceBra, OverlapLabel, OverlapXlm, DiffuseFile )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
!!$    call DiffuseMat.Read( uid )
!!$    close( uid )
!!$    NRows = DiffuseMat.NColumns()
    !***
    !***
    NRows = 100
    !***
    !***
    call DiffuseMat.Free()
    !
    call Self.Block.InitFull(NRows,NCols)
    !
  end subroutine PWCLocBlockBuild


  subroutine PWCLocBlockSave( Self )
    class(PWCLocBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    call Self.Save( FileName )
  end subroutine PWCLocBlockSave


  subroutine PWCLocBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm ) 
    class(PWCLocBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: BlockName
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    Self.NumBsNotOverlappingDiffuseBra = GetNumBsNotOverlappingDiff( Self.SpaceBra )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
  end subroutine PWCLocBlockLoad



  subroutine PWCLocBlockCondition( Self, Conditioner, SRCLoc )
    class(PWCLocBlock)          , intent(inout) :: Self
    class(ClassConditionerBlock), intent(inout) :: Conditioner
    class(SRCLocBlock)          , intent(in)    :: SRCLoc
    !
    type(ClassIrrep), pointer :: DiffIrrep
    type(ClassXlm) :: Xlm, OverlapXlm
    type(ClassMatrix) :: ConditionMat, SRCLocMat
    character(len=:), allocatable :: StorageDir, Axis
    real(kind(1d0)), allocatable :: Array(:,:)
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    !
    call Xlm.Init( Self.SpaceBra.GetL(), Self.SpaceBra.GetM() )
    call Conditioner.SetXlm( Xlm )
    !
    OverlapXlm = GetOperatorXlm(OverlapLabel,Axis)
    allocate( DiffIrrep, source = GetDiffuseIrrep(Self.SpaceBra,OverlapXlm) )
    call Conditioner.SetIrrep( DiffIrrep )
    !
    call Conditioner.ReadBlock( StorageDir, ConditionMat )
    !
    call Self.Block.Multiply( ConditionMat, 'Left', 'T' )
    call ConditionMat.Free()
    !
    call Conditioner.ReadDiffuseBlock( StorageDir, ConditionMat )
    if ( ConditionMat.IsInitialized() ) then
       call SRCLoc.Block.FetchMatrix( Array )
       SRCLocMat = Array(size(Array,1)-ConditionMat.NRows()+1:,:)
       deallocate( Array )
       call SRCLocMat.Multiply( ConditionMat, 'Left', 'T' )
       call ConditionMat.Free()
       call Self.Block.Add( SRCLocMat )
    end if
    !
  end subroutine PWCLocBlockCondition



  !--------------------------------------------------
  ! Methods for the PWCSRCBlock
  !-------------------------------------------------

  subroutine PWCSRCBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(PWCSRCBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: SRCSRCQCFile, DiffuseFile
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: SRCSRCMat
    type(ClassMatrix) :: DiffuseMat
    type(ClassIrrep), pointer :: OpIrrep, BraIrrep, KetIrrep
    type(ClassXlm) :: OverlapXLM
    !
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
!!$    Self.NumBsNotOverlappingDiffuseBra = GetNumBsNotOverlappingDiff( Self.SpaceBra )
    !
    !.. Get the dimension of the short range part.
!!$    if ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call GetCloseCouplingQCFileName( SpaceKet, SpaceKet, OverlapLabel, OpXlm, GetQCFileExtension(), SRCSRCQCFile )
!!$    else
!!$       call GetCloseCouplingQCFileName( SpaceKet, SpaceKet, OpLabel, OpXlm, GetQCFileExtension(), SRCSRCQCFile )
!!$    end if
!!$    call OpenFile( SRCSRCQCFile, uid, 'read', 'formatted' )
!!$    call SRCSRCMat.Read( uid )
!!$    close ( uid )
!!$    NCols = SRCSRCMat.NColumns()
!!$    call SRCSRCMat.Free()
!!$    !
!!$    !.. Get the dimension of the B-splines part.
!!$    call OverlapXLM.Init(0,0)
!!$    call GetBsplinesFileName( SpaceBra, SpaceBra, OverlapLabel, OverlapXlm, DiffuseFile )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
!!$    call DiffuseMat.Read( uid )
!!$    close( uid )
!!$    !
!!$    call DiffuseMat.Transpose()
!!$    NRows = DiffuseMat.NRows()
!!$    call DiffuseMat.Free()
!!$    !
    !***
    !***
    NRows = 100
    NCols = 100
    !***
    !***
    call Self.Block.InitFull(NRows,NCols)
    !
!!$    !.. Determines if the operator is allowed by symmetry
!!$    BraIrrep => Self.SpaceBra.GetTotIrrep()
!!$    KetIrrep => Self.SpaceKet.GetTotIrrep()
!!$    OpIrrep => OpXlm.GetIrrep( BraIrrep.GetGroup() )
!!$    if ( .not.ValidIrreps(BraIrrep,OpIrrep,KetIrrep) ) return
!!$    !
!!$    !.. Select the operator
!!$    if ( OpLabel .is. OverlapLabel ) then
!!$       call Self.BuildOverlap( Force )
!!$    elseif ( OpLabel .is. HamiltonianLabel ) then
       call Self.BuildHamiltonian( Force )
!!$    elseif ( OpLabel .is. KinEnergyLabel ) then
!!$       call Self.BuildKinEnergy( Force )
!!$    elseif ( (OpLabel .is. DIPOLE_LABEL//LengthLabel) .or. &
!!$         (OpLabel .is. DIPOLE_LABEL//VelocityLabel) ) then
!!$       call Self.BuildDipole( Force )
!!$    elseif ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       ! It is zero.
!!$       return
!!$    else
!!$       call Assert( 'Invalid operator label: '//OpLabel//'.' )
!!$    end if
!!$    !
!!$    call CoupleSpins( Self.SpaceBra, Self.SpaceKet, Self.Block )
    !
  end subroutine PWCSRCBlockBuild



  subroutine PWCSRCBlockSave( Self )
    class(PWCSRCBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
!!$    allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
!!$    call Self.Save( FileName )
  end subroutine PWCSRCBlockSave


  subroutine PWCSRCBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm ) 
    class(PWCSRCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: BlockName
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    Self.NumBsNotOverlappingDiffuseBra = GetNumBsNotOverlappingDiff( Self.SpaceBra )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
    !
  end subroutine PWCSRCBlockLoad



  subroutine PWCSRCBlockCondition( Self, Conditioner, SRCSRC )
    class(PWCSRCBlock)          , intent(inout) :: Self
    class(ClassConditionerBlock), intent(inout) :: Conditioner
    class(SRCSRCBlock)          , intent(in)    :: SRCSRC
    !
    type(ClassIrrep), pointer :: DiffIrrep
    type(ClassXlm) :: Xlm, OverlapXlm
    type(ClassMatrix) :: ConditionMat, SRCSRCMat
    character(len=:), allocatable :: StorageDir, Axis
    real(kind(1d0)), allocatable :: Array(:,:)
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    !
    call Xlm.Init( Self.SpaceBra.GetL(), Self.SpaceBra.GetM() )
    call Conditioner.SetXlm( Xlm )
    !
    OverlapXlm = GetOperatorXlm(OverlapLabel,Axis)
    allocate( DiffIrrep, source = GetDiffuseIrrep(Self.SpaceBra,OverlapXlm) )
    call Conditioner.SetIrrep( DiffIrrep )
    !
    call Conditioner.ReadBlock( StorageDir, ConditionMat )
    !
    call Self.Block.Multiply( ConditionMat, 'Left', 'T' )
    call ConditionMat.Free()
    !
    call Conditioner.ReadDiffuseBlock( StorageDir, ConditionMat )
    if ( ConditionMat.IsInitialized() ) then
       call SRCSRC.Block.FetchMatrix( Array )
       SRCSRCMat = Array(size(Array,1)-ConditionMat.NRows()+1:,:)
       deallocate( Array )
       call SRCSRCMat.Multiply( ConditionMat, 'Left', 'T' )
       call ConditionMat.Free()
       call Self.Block.Add( SRCSRCMat )
    end if
    !
  end subroutine PWCSRCBlockCondition



  subroutine PWCSRCBlockBuildOverlap( Self, Force )
    class(PWCSRCBlock)   , intent(inout) :: Self
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, AuxMat
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !
    !.. Get the dimension of the B-splines part.
    call GetDiffuseFileName( Self.SpaceBra, Self.OpLabel, Self.OpXlm, DiffuseFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
    call DiffuseMat.Read( uid )
    close( uid )
    call DiffuseMat.Transpose()
    !
    call PutBlockInRigthPosition( Self.Block, DiffuseMat, AuxMat )
    call DiffuseMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
  end subroutine PWCSRCBlockBuildOverlap



  subroutine PWCSRCBlockBuildHamiltonian( Self, Force )
    class(PWCSRCBlock), intent(inout) :: Self
    logical           , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, DiffuseOverlapFile, BlockFileName, StorageDir
    integer           :: uid
    type(ClassMatrix) :: DiffuseMat, DiffuseOverlapMat, AuxMat
    real(kind(1d0))   :: PIEnergy
    !
!!$    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
!!$    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
!!$    !
!!$    !..Multipole contribution
!!$    call ComputeMultipoleContribution( &
!!$         Self.SpaceBra, Self.SpaceKet, &
!!$         Self.Block                  )
!!$    !
!!$    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
!!$    !
!!$    !.. Get the current operator monoelectronic block.
!!$    call GetDiffuseFileName( Self.SpaceBra, Self.OpLabel, Self.OpXlm, DiffuseFile )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
!!$    call DiffuseMat.Read( uid )
!!$    close( uid )
!!$    call DiffuseMat.Transpose()
!!$    !
!!$    !.. Get the overlap operator monoelectronic block.
!!$    call GetDiffuseFileName( Self.SpaceBra, OverlapLabel, Self.OpXlm, DiffuseOverlapFile )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( DiffuseOverlapFile, uid, 'read', 'formatted' )
!!$    call DiffuseOverlapMat.Read( uid )
!!$    close( uid )
!!$    call DiffuseOverlapMat.Transpose()
!!$    !
!!$    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
!!$    PIEnergy = Self.SpaceBra.GetPIEnergy(StorageDir)
!!$    call DiffuseOverlapMat.Multiply( PIEnergy )
!!$    call DiffuseOverlapMat.Add( DiffuseMat )
!!$    call DiffuseMat.Free()
!!$    !
!!$    call PutBlockInRigthPosition( Self.Block, DiffuseOverlapMat, AuxMat )
!!$    call DiffuseOverlapMat.Free()
!!$    call Self.Block.Add( AuxMat )
!!$    call AuxMat.Free()
    !
  end subroutine PWCSRCBlockBuildHamiltonian




  subroutine PWCSRCBlockBuildDipole( Self, Force )
    class(PWCSRCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, DiffuseOverlapFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, DiffuseOverlapMat, AuxMat
    real(kind(1d0))   :: PIsDipole
    type(ClassXlm)    :: OverlapXlm
    logical :: exist
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    call OverlapXlm.Init(0,0)
    !.. Get the overlap operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceBra, OverlapLabel, OverlapXlm, DiffuseOverlapFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseOverlapFile, uid, 'read', 'formatted' )
    call DiffuseOverlapMat.Read( uid )
    close( uid )
    call DiffuseOverlapMat.Transpose()
    !
    PIsDipole = GetPIsDipole( Self.SpaceBra, Self.SpaceKet, Self.OpXlm )
    call DiffuseOverlapMat.Multiply( PIsDipole )
    call PutBlockInRigthPosition( Self.Block, DiffuseOverlapMat, AuxMat )
    call DiffuseOverlapMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceBra, Self.OpLabel, Self.OpXlm, DiffuseFile )
    INQUIRE( File = DiffuseFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
       call DiffuseMat.Read( uid )
       close( uid )
       call DiffuseMat.Transpose( )
       !
       call PutBlockInRigthPosition( Self.Block, DiffuseMat, AuxMat )
       call DiffuseMat.Free()
       call Self.Block.Add( AuxMat )
       call AuxMat.Free()
    end if
    !
  end subroutine PWCSRCBlockBuildDipole


  subroutine PWCSRCBlockBuildKinEnergy( Self, Force )
    class(PWCSRCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: DiffuseFile, DiffuseOverlapFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: DiffuseMat, DiffuseOverlapMat, AuxMat
    real(kind(1d0))   :: PIsKinEnergy
    type(ClassXlm)    :: OverlapXlm
    logical :: exist
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    call OverlapXlm.Init(0,0)
    !.. Get the overlap operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceBra, OverlapLabel, OverlapXlm, DiffuseOverlapFile )
    !.. Assumes saved with format.
    call OpenFile( DiffuseOverlapFile, uid, 'read', 'formatted' )
    call DiffuseOverlapMat.Read( uid )
    close( uid )
    call DiffuseOverlapMat.Transpose()
    !
    PIsKinEnergy = GetPIsKinEnergy( Self.SpaceBra, Self.SpaceKet )
    call DiffuseOverlapMat.Multiply( PIsKinEnergy )
    call PutBlockInRigthPosition( Self.Block, DiffuseOverlapMat, AuxMat )
    call DiffuseOverlapMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetDiffuseFileName( Self.SpaceBra, Self.OpLabel, Self.OpXlm, DiffuseFile )
    INQUIRE( File = DiffuseFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( DiffuseFile, uid, 'read', 'formatted' )
       call DiffuseMat.Read( uid )
       close( uid )
       call DiffuseMat.Transpose( )
       !
       call PutBlockInRigthPosition( Self.Block, DiffuseMat, AuxMat )
       call DiffuseMat.Free()
       call Self.Block.Add( AuxMat )
       call AuxMat.Free()
    end if
    !
  end subroutine PWCSRCBlockBuildKinEnergy


  !--------------------------------------------------
  ! Methods for the PWCPWCBlock
  !-------------------------------------------------


  subroutine PWCPWCBlockBuild( Self, SpaceBra, SpaceKet, OpLabel, OpXlm, Force )
    class(PWCPWCBlock)   , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: BsFileBra, BsFileKet
    integer :: NRows, NCols, uid
    type(ClassMatrix) :: BsMatBra, BsMatKet
    type(ClassIrrep), pointer :: OpIrrep, BraIrrep, KetIrrep
    !
    type(ClassXlm) :: OverlapXLM
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
!!$    Self.NumBsNotOverlappingDiffuseBra = GetNumBsNotOverlappingDiff( Self.SpaceBra )
!!$    Self.NumBsNotOverlappingDiffuseKet = GetNumBsNotOverlappingDiff( Self.SpaceKet )
!!$    !
!!$    call OverlapXLM.Init(0,0)
!!$    !.. Get the dimension of the bra B-splines part.
!!$    call GetBsplinesFileName( SpaceBra, SpaceBra, OverlapLabel, OverlapXlm, BsFileBra )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( BsFileBra, uid, 'read', 'formatted' )
!!$    call BsMatBra.Read( uid )
!!$    close( uid )
!!$    NRows = BsMatBra.NRows()
!!$    call BsMatBra.Free()
!!$    !
!!$    !.. Get the dimension of the ket B-splines part.
!!$    call GetBsplinesFileName( SpaceKet, SpaceKet, OverlapLabel, OverlapXlm, BsFileKet )
!!$    !.. Assumes saved with format.
!!$    call OpenFile( BsFileKet, uid, 'read', 'formatted' )
!!$    call BsMatKet.Read( uid )
!!$    close( uid )
!!$    NCols = BsMatKet.NColumns()
!!$    call BsMatKet.Free()
    !***
    NRows = 100
    NCols = 100
    !***
    !
    call Self.Block.InitFull(NRows,NCols)
    !
!!$    !.. Determines if the operator is allowed by symmetry
!!$    BraIrrep => Self.SpaceBra.GetTotIrrep()
!!$    KetIrrep => Self.SpaceKet.GetTotIrrep()
!!$    OpIrrep => OpXlm.GetIrrep( BraIrrep.GetGroup() )
!!$    if ( .not.ValidIrreps(BraIrrep,OpIrrep,KetIrrep) ) return
!!$    !
!!$    !.. Select the operator
!!$    if ( OpLabel .is. OverlapLabel ) then
!!$       call Self.BuildOverlap( Force )
!!$    elseif ( OpLabel .is. HamiltonianLabel ) then
!!$       call Self.BuildHamiltonian( Force )
!!$    elseif ( OpLabel .is. KinEnergyLabel ) then
!!$       call Self.BuildKinEnergy( Force )
!!$    elseif ( (OpLabel .is. DIPOLE_LABEL//LengthLabel) .or. &
!!$         (OpLabel .is. DIPOLE_LABEL//VelocityLabel) ) then
!!$       call Self.BuildDipole( Force )
!!$    elseif ( OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel ) then
!!$       call Self.BuildCAP( Force )
!!$    else
!!$       call Assert( 'Invalid operator label: '//OpLabel//'.' )
!!$    end if
!!$    !
!!$    call CoupleSpins( Self.SpaceBra, Self.SpaceKet, Self.Block )
    !
  end subroutine PWCPWCBlockBuild


  subroutine PWCPWCBlockSave( Self )
    class(PWCPWCBlock), intent(in) :: Self
    character(len=:), allocatable :: FileName
    allocate( FileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    call Self.Save( FileName )
  end subroutine PWCPWCBlockSave


  subroutine PWCPWCBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, OpXlm ) 
    class(PWCPWCBlock)      , intent(inout) :: Self
    class(*)       , target , intent(in)    :: SpaceBra
    class(*)       , target , intent(in)    :: SpaceKet
    character(len=*)        , intent(in)    :: OpLabel
    class(ClassXlm), target , intent(in)    :: OpXlm
    !
    character(len=:), allocatable :: BlockName
    !
    call Self.Init( OpLabel, OpXlm )
    allocate( Self.SpaceBra, source = SpaceBra )
    allocate( Self.SpaceKet, source = SpaceKet )
    !
    Self.NumBsNotOverlappingDiffuseBra = GetNumBsNotOverlappingDiff( Self.SpaceBra )
    Self.NumBsNotOverlappingDiffuseKet = GetNumBsNotOverlappingDiff( Self.SpaceKet )
    !
    allocate( BlockName, source = Self.GetBlockName( SpaceBra, SpaceKet ) )
    call CheckFileMustBePresent( BlockName )
    call Self.Read( BlockName )
    !
  end subroutine PWCPWCBlockLoad



  subroutine PWCPWCBlockCondition( Self, Conditioner, SRCSRC, PWCSRC, SRCPWC )
    class(PWCPWCBlock)          , intent(inout) :: Self
    class(ClassConditionerBlock), intent(inout) :: Conditioner
    class(SRCSRCBlock)          , intent(in)    :: SRCSRC
    class(PWCSRCBlock)          , intent(in)    :: PWCSRC
    class(SRCPWCBlock)          , intent(in)    :: SRCPWC
    !
    type(ClassIrrep), pointer :: BraDiffIrrep, KetDiffIrrep
    type(ClassXlm) :: BraXlm, KetXlm, OverlapXlm
    type(ClassMatrix) :: BraConditionMat, KetConditionMat, BraDiffConditionMat, KetDiffConditionMat, PWCSRCMat, SRCPWCMat, SRCSRCMat
    character(len=:), allocatable :: StorageDir, Axis
    real(kind(1d0)), allocatable :: Array(:,:)
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    OverlapXlm = GetOperatorXlm(OverlapLabel,Axis)
    !
    call BraXlm.Init( Self.SpaceBra.GetL(), Self.SpaceBra.GetM() )
    call Conditioner.SetXlm( BraXlm )
    allocate( BraDiffIrrep, source = GetDiffuseIrrep(Self.SpaceBra,OverlapXlm) )
    call Conditioner.SetIrrep( BraDiffIrrep )
    call Conditioner.ReadBlock( StorageDir, BraConditionMat )
    call Conditioner.ReadDiffuseBlock( StorageDir, BraDiffConditionMat )
    !
    call KetXlm.Init( Self.SpaceKet.GetL(), Self.SpaceKet.GetM() )
    call Conditioner.SetXlm( KetXlm )
    allocate( KetDiffIrrep, source = GetDiffuseIrrep(Self.SpaceKet,OverlapXlm) )
    call Conditioner.SetIrrep( KetDiffIrrep )
    call Conditioner.ReadBlock( StorageDir, KetConditionMat )
    call Conditioner.ReadDiffuseBlock( StorageDir, KetDiffConditionMat )
    !
    call Self.Block.Multiply( BraConditionMat, 'Left', 'T' )
    call Self.Block.Multiply( KetConditionMat, 'Right', 'N' )
    !
    !
    if ( BraDiffConditionMat.IsInitialized() .and. &
         KetDiffConditionMat.IsInitialized() ) then
       !
       call SRCSRC.Block.FetchMatrix( Array )
       SRCSRCMat = Array(size(Array,1)-BraDiffConditionMat.NRows()+1:,size(Array,2)-KetDiffConditionMat.NRows()+1:)
       deallocate( Array )
       call SRCSRCMat.Multiply( BraDiffConditionMat, 'Left', 'T' )
       call SRCSRCMat.Multiply( KetDiffConditionMat, 'Right', 'N' )
       call Self.Block.Add( SRCSRCMat )
       call SRCSRCMat.Free()
       !
       call PWCSRC.Block.FetchMatrix( Array )
       PWCSRCMat = Array(:,size(Array,2)-KetDiffConditionMat.NRows()+1:)
       deallocate( Array )
       call PWCSRCMat.Multiply( BraConditionMat, 'Left', 'T' )
       call PWCSRCMat.Multiply( KetDiffConditionMat, 'Right', 'N' )
       call Self.Block.Add( PWCSRCMat )
       call PWCSRCMat.Free()
       !
       call SRCPWC.Block.FetchMatrix( Array )
       SRCPWCMat = Array(size(Array,1)-BraDiffConditionMat.NRows()+1:,:)
       deallocate( Array )
       call SRCPWCMat.Multiply( BraDiffConditionMat, 'Left', 'T' )
       call SRCPWCMat.Multiply( KetConditionMat, 'Right', 'N' )
       call Self.Block.Add( SRCPWCMat )
       call SRCPWCMat.Free()
       !
    end if
    !
  end subroutine PWCPWCBlockCondition



  subroutine PWCPWCBlockBuildOverlap( Self, Force )
    class(PWCPWCBlock)   , intent(inout) :: Self
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: BsFile, BlockFileName
    type(ClassMatrix) :: BsMat, AuxMat
    integer :: lBra, mBra, lKet, mKet, uid
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    lBra = Self.SpaceBra.GetL()
    mBra = Self.SpaceBra.GetM()
    lKet = Self.SpaceKet.GetL()
    mKet = Self.SpaceKet.GetM()
    if ( lBra/=lKet .or. mBra/=mKet ) return
    !
    !.. Get the dimension of the B-splines part.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, Self.OpLabel, Self.OpXlm, BsFile )
    !.. Assumes saved with format.
    call OpenFile( BsFile, uid, 'read', 'formatted' )
    call BsMat.Read( uid )
    close( uid )
    !
    call PutBlockInRigthPosition( Self.Block, BsMat, AuxMat )
    call BsMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
  end subroutine PWCPWCBlockBuildOverlap



  subroutine PWCPWCBlockBuildHamiltonian( Self, Force )
    class(PWCPWCBlock)   , intent(inout) :: Self
    logical                 , intent(in)    :: Force
    !
    character(len=:), allocatable :: BsFile, BsOverlapFile, BlockFileName, StorageDir
    type(ClassMatrix) :: BsMat, BsOverlapMat, AuxMat
    integer :: lBra, mBra, lKet, mKet, uid
    real(kind(1d0)) :: PIEnergy
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    !..Multipole contribution
    call ComputeMultipoleContribution( &
         Self.SpaceBra, Self.SpaceKet, &
         Self.Block                  )
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    lBra = Self.SpaceBra.GetL()
    mBra = Self.SpaceBra.GetM()
    lKet = Self.SpaceKet.GetL()
    mKet = Self.SpaceKet.GetM()
    if ( lBra/=lKet .or. mBra/=mKet ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, Self.OpLabel, Self.OpXlm, BsFile )
    !.. Assumes saved with format.
    call OpenFile( BsFile, uid, 'read', 'formatted' )
    call BsMat.Read( uid )
    close( uid )
    !
    !.. Get the overlap operator monoelectronic block. The associated Xlm is the same that for the hamiltonian
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, OverlapLabel, Self.OpXlm, BsOverlapFile )
    !.. Assumes saved with format.
    call OpenFile( BsOverlapFile, uid, 'read', 'formatted' )
    call BsOverlapMat.Read( uid )
    close( uid )
    !
    allocate( StorageDir, source = GetStorageDir(Self.SpaceBra) )
    PIEnergy = Self.SpaceBra.GetPIEnergy(StorageDir)
    call BsOverlapMat.Multiply( PIEnergy )
    call BsOverlapMat.Add( BsMat )
    call BsMat.Free()
    !
    call PutBlockInRigthPosition( Self.Block, BsOverlapMat, AuxMat )
    call BsOverlapMat.Free()
    call Self.Block.Add( AuxMat )
    call AuxMat.Free()
    !
  end subroutine PWCPWCBlockBuildHamiltonian



  subroutine PWCPWCBlockBuildDipole( Self, Force )
    class(PWCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: BsFile, BsOverlapFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: BsMat, BsOverlapMat, AuxMat
    real(kind(1d0))   :: PIsDipole
    type(ClassXlm)    :: OverlapXlm
    logical :: exist
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    call OverlapXlm.Init(0,0)
    !.. Get the overlap operator monoelectronic block.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, OverlapLabel, OverlapXlm, BsOverlapFile )
    INQUIRE( File = BsOverlapFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( BsOverlapFile, uid, 'read', 'formatted' )
       call BsOverlapMat.Read( uid )
       close( uid )
       !
       PIsDipole = GetPIsDipole( Self.SpaceBra, Self.SpaceKet, Self.OpXlm )
       call BsOverlapMat.Multiply( PIsDipole )
       call PutBlockInRigthPosition( Self.Block, BsOverlapMat, AuxMat )
       call BsOverlapMat.Free()
       Self.Block = AuxMat
       call AuxMat.Free()
    end if
    !
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, Self.OpLabel, Self.OpXlm, BsFile )
    INQUIRE( File = BsFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( BsFile, uid, 'read', 'formatted' )
       call BsMat.Read( uid )
       close( uid )
       !
       call PutBlockInRigthPosition( Self.Block, BsMat, AuxMat )
       call BsMat.Free()
       call Self.Block.Add( AuxMat )
       call AuxMat.Free()
    end if
    !
  end subroutine PWCPWCBlockBuildDipole


  subroutine PWCPWCBlockBuildKinEnergy( Self, Force )
    class(PWCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: BsFile, BsOverlapFile, BlockFileName
    integer :: uid
    type(ClassMatrix) :: BsMat, BsOverlapMat, AuxMat
    real(kind(1d0))   :: PIsKinEnergy
    type(ClassXlm)    :: OverlapXlm
    logical :: exist
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    call OverlapXlm.Init(0,0)
    !.. Get the overlap operator monoelectronic block.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, OverlapLabel, OverlapXlm, BsOverlapFile )
    INQUIRE( File = BsOverlapFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( BsOverlapFile, uid, 'read', 'formatted' )
       call BsOverlapMat.Read( uid )
       close( uid )
       !
       PIsKinEnergy = GetPIsKinEnergy( Self.SpaceBra, Self.SpaceKet )
       call BsOverlapMat.Multiply( PIsKinEnergy )
       call PutBlockInRigthPosition( Self.Block, BsOverlapMat, AuxMat )
       call BsOverlapMat.Free()
       Self.Block = AuxMat
       call AuxMat.Free()
    end if
    !
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    !.. Get the current operator monoelectronic block.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, Self.OpLabel, Self.OpXlm, BsFile )
    INQUIRE( File = BsFile, exist = exist )
    !
    if ( exist ) then
       !.. Assumes saved with format.
       call OpenFile( BsFile, uid, 'read', 'formatted' )
       call BsMat.Read( uid )
       close( uid )
       !
       call PutBlockInRigthPosition( Self.Block, BsMat, AuxMat )
       call BsMat.Free()
       call Self.Block.Add( AuxMat )
       call AuxMat.Free()
    end if
    !
  end subroutine PWCPWCBlockBuildKinEnergy


  !..Essentially is the same subroutine that the overlap, because
  !.. automatically loads the correct operator matrix.
  subroutine PWCPWCBlockBuildCAP( Self, Force )
    class(PWCPWCBlock)   , intent(inout) :: Self
    logical              , intent(in)    :: Force
    !
    character(len=:), allocatable :: BsFile, BlockFileName
    type(ClassMatrix) :: BsMat, AuxMat
    integer :: lBra, mBra, lKet, mKet, uid
    !
    allocate( BlockFileName, source = Self.GetBlockName(Self.SpaceBra,Self.SpaceKet) )
    if ( .not.Force .and. Self.ValidBlock(BlockFileName) ) return
    !
    if ( (Self.SpaceBra.GetPILabel() .isnt. Self.SpaceKet.GetPILabel()) ) return
    !
    lBra = Self.SpaceBra.GetL()
    mBra = Self.SpaceBra.GetM()
    lKet = Self.SpaceKet.GetL()
    mKet = Self.SpaceKet.GetM()
    if ( lBra/=lKet .or. mBra/=mKet ) return
    !
    !.. Get the dimension of the B-splines part.
    call GetBsplinesFileName( Self.SpaceBra, Self.SpaceKet, Self.OpLabel, Self.OpXlm, BsFile )
    !.. Assumes saved with format.
    call OpenFile( BsFile, uid, 'read', 'formatted' )
    call BsMat.Read( uid )
    close( uid )
    !
    call PutBlockInRigthPosition( Self.Block, BsMat, AuxMat )
    call BsMat.Free()
    Self.Block = AuxMat
    call AuxMat.Free()
    !
  end subroutine PWCPWCBlockBuildCAP



end module ModuleElementaryElectronicSpaceOperators
