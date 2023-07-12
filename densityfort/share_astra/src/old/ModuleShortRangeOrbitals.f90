module ModuleShortRangeOrbitals


  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleIO

  use ModuleConstants
  use ModuleGroups
  use ModuleSymmetryAdaptedSphericalHarmonics
  use ModuleMatrix

  implicit none

  private

  character(len=*), parameter :: DiffuseOrbDir          = "DiffuseOrbitals"
  character(len=*), parameter :: ConditionerDir         = "Conditioners"
  character(len=*), parameter :: DiffuseLabel           = "Diffuse"
  character(len=*), parameter :: MethodRM               = 'RM'
  character(len=*), parameter :: MethodLI               = 'LI'
  character(len=*), parameter :: MethodLIO              = 'LIO'
  character(len=*), parameter :: PerfectProjectorLabel  = "PerfectProjector"
  character(len=*), parameter :: QCConditionerLabel     = "QC"

  !**** The interface name must changed to a root filename to reflect the variable irrep
  character(len=*), parameter :: DiffuseOrbInterface = "interface"
  character(len=*), parameter :: DiffuseMonOverlapLabel = "MonOverlap"
  logical         , parameter :: PlotDiffuseOrb      = .false.
  !> Maximum number of monomial or exponential exponents
  !! intended to be read.
  integer         , parameter :: NumMaxExponents  = 10000


  type, private :: ClassDiffuseOrbitals
     ! {{{ private attributes

     private
     
     !> Maximum angular momentum the space can represent.
     integer                   :: Lmax
     !> Irrep of the orbitals 
     type(ClassIrrep), pointer :: Irrep

     !> Number of diffuse orbitals with symmetry specified by Irrep, 
     !! and computed by the Quantum Chemistry code so to be orthogonal 
     !! to the molecular orbitals of the same symmetry and to be 
     !! linearly-independent.
     integer :: NOrbitals

     !> Specifies whether the monomials corresponds to the 
     !! cartesian (.true.) or spherical functions (.false.)%
     logical :: CartesianFunc

     !> Number of symmetry-adapted cartesian (or spherical) Gaussian functions
     integer :: NCartesianGaussians

     !> Array of the cartesian monomial exponents (or spherical ) as read from the interace file.
     integer, allocatable :: CartMonExponents(:,:)

     !> Array of exponents of the exponential function.
     real(kind(1d0)), allocatable :: GaussExponents(:)

     !> Matrix of the expansion coefficients of the diffuse 
     !! orbitals (One orbital per column, NColumns = NOrbitals) 
     !! on the basis of symmetry adapted cartesian (or spherical)
     !! Gaussian functions (Number of Rows)%
     real(kind(1d0)), allocatable :: CartesianCoefficients(:,:)

     !> Full transformation matrix (all the cartesian (or spherical) monomials and Xlm multiplied by the used exponentials functions) to pass from Xlm to cartesian representation.
     complex(kind(1d0)), allocatable :: FullXlmToCartArray(:,:)
     !>
     type(ClassCartesianSymmetricSet) :: CartesianSymSet
     !
     type(ClassSpherSymmetricSet) :: SpherSymSet

     ! }}} 
   contains
     generic, public :: init                             => ClassDiffuseOrbitalsInit
     generic, public :: Load                             => ClassDiffuseOrbitalsLoad
     generic, public :: SetIrrep                         => ClassDiffuseOrbitalsSetIrrep
     generic, public :: SetLmax                          => ClassDiffuseOrbitalsSetLmax
     generic, public :: GetNOrbitals                     => ClassDiffuseOrbitalsGetNOrbitals
     generic, public :: ComputeTransfMat                 => ClassDiffuseOrbitalsComputeTransfMat
     generic, public :: GetTransfMat                     => ClassDiffuseOrbitalsGetTransfMat
     generic, public :: GetCartesianSet                  => ClassDiffuseOrbitalsGetCartesianSet
     generic, public :: GetSpherSet                      => ClassDiffuseOrbitalsGetSpherSet
     generic, public :: ComputeFullTransfMat             => ClassDiffuseOrbitalsComputeFullTransfMat
     generic, public :: ComputeFullXlmToCartTransfMat    => ClassDiffuseOrbitalsComputeFullXlmToCartTransfMat
     generic, public :: ComputeFullXlmToSpherTransfMat   => ClassDiffuseOrbitalsComputeFullXlmToSpherTransfMat
     generic, public :: FetchFullTransfMat               => ClassDiffuseOrbitalsFetchFullTransfMat
     generic, public :: FetchFullXlmToCartesianTransfMat => ClassDiffuseOrbitalsFetchFullXlmToCartesianTransfMat
     generic, public :: FetchFullXlmToSpherTransfMat     => ClassDiffuseOrbitalsFetchFullXlmToSpherTransfMat
     generic, public :: FetchTransfMatrix                => ClassDiffuseOrbitalsFetchtransfMatrix
     generic, public :: FetchXlmToCartesianMatrix        => ClassDiffuseOrbitalsFetchXlmToCartesianMatrix
     generic, public :: FetchXlmToSpherMatrix            => ClassDiffuseOrbitalsFetchXlmToSpherMatrix
     generic, public :: FetchCartMonInterface            => ClassDiffuseOrbitalsFetchCartMonInterface
     generic, public :: FetchCartMonInitialized          => ClassDiffuseOrbitalsFetchCartMonInitialized
     generic, public :: FetchXlmTransfIndexes            => ClassDiffuseOrbitalsFetchXlmTransfIndexes
     generic, public :: FetchCartesianCoeff              => ClassDiffuseOrbitalsFetchCartesianCoeff
     generic, public :: show                             => ClassDiffuseOrbitalsShow
     generic, public :: free                             => ClassDiffuseOrbitalsFree
     ! {{{ private procedures

     procedure, private :: ClassDiffuseOrbitalsInit
     procedure, private :: ClassDiffuseOrbitalsLoad
     procedure, private :: ClassDiffuseOrbitalsSetIrrep
     procedure, private :: ClassDiffuseOrbitalsSetLmax
     procedure, private :: ClassDiffuseOrbitalsGetNOrbitals
     procedure, private :: ClassDiffuseOrbitalsComputeTransfMat
     procedure, private :: ClassDiffuseOrbitalsGetTransfMat
     procedure, private :: ClassDiffuseOrbitalsGetCartesianSet
     procedure, private :: ClassDiffuseOrbitalsGetSpherSet
     procedure, private :: ClassDiffuseOrbitalsComputeFullXlmToCartTransfMat
     procedure, private :: ClassDiffuseOrbitalsComputeFullXlmToSpherTransfMat
     procedure, private :: ClassDiffuseOrbitalsComputeFullTransfMat
     procedure, private :: ClassDiffuseOrbitalsFetchFullTransfMat
     procedure, private :: ClassDiffuseOrbitalsFetchFullXlmToCartesianTransfMat
     procedure, private :: ClassDiffuseOrbitalsFetchFullXlmToSpherTransfMat
     procedure, private :: ClassDiffuseOrbitalsFetchXlmToCartesianMatrix
     procedure, private :: ClassDiffuseOrbitalsFetchXlmToSpherMatrix
     procedure, private :: ClassDiffuseOrbitalsFetchTransfMatrix
     procedure, private :: ClassDiffuseOrbitalsFetchCartMonInterface
     procedure, private :: ClassDiffuseOrbitalsFetchCartMonInitialized
     procedure, private :: ClassDiffuseOrbitalsFetchXlmTransfIndexes
     procedure, private :: ClassDiffuseOrbitalsFetchCartesianCoeff
     procedure, private :: ClassDiffuseOrbitalsShow
     procedure, private :: ClassDiffuseOrbitalsFree
     final              :: ClassDiffuseOrbitalsFinal

     ! }}}
  end type ClassDiffuseOrbitals



  type, public :: ClassShortRangeSymOrbitals
     ! {{{ private attributes

     private
     
     !> Maximum Angular Momentum the space can represent
     integer                   :: Lmax
     !> Irrep of the orbitals 
     type(ClassIrrep), pointer :: Irrep

     !> Number of localized molecular orbitals (which,
     !! by hypothesis, do not overlap with the BSplines) 
     !! with the symmetry specified by Irrep
     integer                    :: NMolecularOrbitals
     !> Diffuse orbitals class.
     type(ClassDiffuseOrbitals) :: DiffuseOrbitals

     ! }}} 
   contains
     generic, public :: init                                    => ClassShortRangeSymOrbitalsInit
     !> Read the file that contains all the information on the SR orbitals
     !! from the QC interface, for a given symmetry.
     generic, public :: Load                                    => ClassShortRangeSymOrbitalsLoad
     generic, public :: show                                    => ClassShortRangeSymOrbitalsShow
     generic, public :: GetTotNumOrbitals                       => ClassShortRangeSymOrbitalsGetTotNumOrbitals
     generic, public :: GetNumDiffuseOrbitals                   => ClassShortRangeSymOrbitalsGetNumDiffuseOrbitals
     generic, public :: GetNumMolOrbitals                       => ClassShortRangeSymOrbitalsGetNumMolOrbitals
     generic, public :: ComputeAndSaveDiffMonOverlap            => ClassSRSOComputeAndSaveDiffMonOverlap
     generic, public :: LoadDiffMonOverlap                      => ClassSRSOLoadDiffMonOverlap
     generic, public :: ComputeTransfMat                        => ClassShortRangeSymOrbitalsComputeTransfMat
     generic, public :: GetTransfMat                            => ClassShortRangeSymOrbitalsGetTransfMat 
     generic, public :: GetLmax                                 => ClassShortRangeSymOrbitalsGetLmax
     generic, public :: GetDiffuseMonOvFileName                 => ClassShortRangeSymOrbitalsGetDiffuseMonOvFileName
     generic, public :: GetCartesianSet                         => ClassShortRangeSymOrbitalsGetCartesianSet
     generic, public :: AcquireInterface                        => ClassShortRangeSymOrbitalsAcquireInterface
     generic, public :: ComputeFullTransfMat                    => ClassShortRangeSymOrbitalsComputeFullTransfMat
     generic, public :: FetchFullTransfMat                      => ClassShortRangeSymOrbitalsFetchFullTransfMat
     generic, public :: FetchXlmToCartesianMatrix               => ClassShortRangeSymOrbitalsFetchXlmToCartesianMatrix
     generic, public :: FetchCartMonInterface                   => ClassShortRangeSymOrbitalsFetchCartMonInterface
     generic, public :: FetchCartMonInitialized                 => ClassShortRangeSymOrbitalsFetchCartMonInitialized
     generic, public :: FetchXlmTransfIndexes                   => ClassShortRangeSymOrbitalsFetchXlmTransfIndexes
     generic, public :: FetchCartesianCoeff                     => ClassShortRangeSymOrbitalsFetchCartesianCoeff
     generic, public :: FetchklmMask                            => ClassShortRangeSymOrbitalsFetchklmMask
     generic, public :: free                                    => ClassShortRangeSymOrbitalsFree
     ! {{{ private procedures

     procedure, private :: ClassShortRangeSymOrbitalsInit
     procedure, private :: ClassShortRangeSymOrbitalsLoad
     procedure, private :: ClassShortRangeSymOrbitalsShow
     procedure, private :: ClassShortRangeSymOrbitalsGetTotNumOrbitals
     procedure, private :: ClassShortRangeSymOrbitalsGetNumDiffuseOrbitals
     procedure, private :: ClassShortRangeSymOrbitalsGetDiffuseMonOvFileName
     procedure, private :: ClassShortRangeSymOrbitalsGetLMax
     procedure, private :: ClassShortRangeSymOrbitalsGetNumMolOrbitals
     procedure, private :: ClassSRSOComputeAndSaveDiffMonOverlap
     procedure, private :: ClassSRSOLoadDiffMonOverlap
     procedure, private :: ClassShortRangeSymOrbitalsComputeTransfMat
     procedure, private :: ClassShortRangeSymOrbitalsGetTransfMat
     procedure, private :: ClassShortRangeSymOrbitalsAcquireInterface
     procedure, private :: ClassShortRangeSymOrbitalsGetCartesianSet
     procedure, private :: ClassShortRangeSymOrbitalsComputeFullTransfMat
     procedure, private :: ClassShortRangeSymOrbitalsFetchFullTransfMat
     procedure, private :: ClassShortRangeSymOrbitalsFetchXlmToCartesianMatrix
     procedure, private :: ClassShortRangeSymOrbitalsFetchCartMonInterface
     procedure, private :: ClassShortRangeSymOrbitalsFetchCartMonInitialized
     procedure, private :: ClassShortRangeSymOrbitalsFetchXlmTransfIndexes
     procedure, private :: ClassShortRangeSymOrbitalsFetchCartesianCoeff
     procedure, private :: ClassShortRangeSymOrbitalsFetchklmMask
     procedure, private :: ClassShortRangeSymOrbitalsFree
     final              :: ClassShortRangeSymOrbitalsFinal

     ! }}}
  end type ClassShortRangeSymOrbitals



  !> Set of the short-range orbitals for all the possible symmetries.
  !! These orbitals need to be loaded once at the beginning.
  !! Afterwards, individual symmetric components can be associated
  !! to several different parent ions in different SRC classes.
  type, public :: ClassShortRangeOrbitals 
     ! {{{ private attributes

     private
     type(ClassGroup), pointer :: Group

     !> Vectors of the symmetry-adapted short-range orbitals
     !! for all the irreps of group 
     type(ClassShortRangeSymOrbitals), allocatable :: SymSRO(:)

     !**** Must use this
     logical, allocatable :: IrrepIsPresent(:)

     ! }}}
   contains
     generic, public :: GetSymSet  => ClassShortRangeOrbitalsGetSymSet
     generic, public :: init       => ClassShortRangeOrbitalsInit
     generic, public :: Load       => ClassShortRangeOrbitalsLoad
     generic, public :: show       => ClassShortRangeOrbitalsShow
     generic, public :: free       => ClassShortRangeOrbitalsFree
     ! {{{ private procedures

     procedure, private :: ClassShortRangeOrbitalsGetSymSet
     procedure, private :: ClassShortRangeOrbitalsInit
     procedure, private :: ClassShortRangeOrbitalsLoad
     procedure, private :: ClassShortRangeOrbitalsFree
     procedure, private :: ClassShortRangeOrbitalsShow
     final :: ClassShortRangeOrbitalsFinal

     ! }}}
  end type ClassShortRangeOrbitals



  type, public :: ClassDiffuseBsplineXlmBlock

     type(ClassIrrep), pointer     :: DiffuseIrrep
     type(ClassXlm)  , pointer     :: OperatorXlm
     character(len=:), allocatable :: OperatorLabel
     type(ClassXlm)  , pointer     :: KetXlm
     type(ClassMatrix)             :: Block

   contains

     generic, public :: init      => ClassDiffuseBsplineXlmBlockInit
     generic, public :: save      => ClassDiffuseBsplineXlmBlockSave
     generic, public :: ReadBlock => ClassDiffuseBsplineXlmBlockReadBlock
     generic, public :: GetFile   => ClassDiffuseBsplineXlmBlockGetFile
     generic, public :: Free      => ClassDiffuseBsplineXlmBlockFree

     procedure, private :: ClassDiffuseBsplineXlmBlockInit
     procedure, private :: ClassDiffuseBsplineXlmBlockSave
     procedure, private :: ClassDiffuseBsplineXlmBlockReadBlock
     procedure, private :: ClassDiffuseBsplineXlmBlockGetFile
     procedure, private :: ClassDiffuseBsplineXlmBlockFree
     
  end type ClassDiffuseBsplineXlmBlock



  type, public :: ClassConditionerBlock
     
     character(len=:), allocatable :: MethodLabel
     real(kind(1d0))               :: Threshold
     integer                       :: NumBsDropBeginning
     integer                       :: NumBsDropEnding
     !> Irreducible representation of the diffuse orbitals.
     type(ClassIrrep), pointer     :: Irrep
     !> Xlm associated to the asymptotic electron.
     type(ClassXlm)  , pointer     :: KetXlm
     !> Transforms the B-splines to the new basis
     type(ClassMatrix)             :: Block
     !> Transforms the diffuse orbitals to the new basis
     type(ClassMatrix)             :: AuxBlock
     !> If FALSE, the overlap matrix among monocentric Gaussian
     !! orbitals is assumed to be the identity matrix. If TRUE
     !! the exact overlap matrix will bw taken into account.
     logical                       :: PerfectProjector

   contains

     generic, public :: init                    => ClassConditionerBlockInit
     generic, public :: save                    => ClassConditionerBlockSave
     generic, public :: ReadBlock               => ClassConditionerBlockReadBlock
     generic, public :: ReadDiffuseBlock        => ClassConditionerBlockReadDiffuseBlock
     generic, public :: GetFile                 => ClassConditionerBlockGetFile
     generic, public :: GetDiffuseFile          => ClassConditionerBlockGetDiffuseFile
     generic, public :: GetLabel                => ClassConditionerBlockGetLabel
     generic, public :: SetXlm                  => ClassConditionerBlockSetXlm
     generic, public :: SetIrrep                => ClassConditionerBlockSetIrrep
     generic, public :: ComputeConditioner      => ClassConditionerBlockComputeConditioner
     generic, public :: ComputeRMConditioner    => ClassConditionerBlockComputeRMConditioner
!!$     generic, public :: ComputeLIConditioner    => ClassConditionerBlockComputeLIConditioner
!!$     generic, public :: ComputeLIOConditioner   => ClassConditionerBlockComputeLIOConditioner
     generic, public :: Free                    => ClassConditionerBlockFree

     procedure, private :: ClassConditionerBlockInit
     procedure, private :: ClassConditionerBlockSave
     procedure, private :: ClassConditionerBlockReadBlock
     procedure, private :: ClassConditionerBlockReadDiffuseBlock
     procedure, private :: ClassConditionerBlockGetFile
     procedure, private :: ClassConditionerBlockGetDiffuseFile
     procedure, private :: ClassConditionerBlockGetLabel
     procedure, private :: ClassConditionerBlockSetXlm
     procedure, private :: ClassConditionerBlockSetIrrep
     procedure, private :: ClassConditionerBlockComputeConditioner
     procedure, private :: ClassConditionerBlockComputeRMConditioner
!!$     procedure, private :: ClassConditionerBlockComputeLIConditioner
!!$     procedure, private :: ClassConditionerBlockComputeLIOConditioner
     procedure, private :: ClassConditionerBlockFree
     
  end type ClassConditionerBlock



  public :: GetDiffuseOrbDir
  public :: GetConditionerDir
  public :: CheckThreshold
  public :: CheckMethod
  public :: CheckBsToRemove
  public :: CheckLastBsIsPresent
  public :: GetNumBsNotOverlappingDiffuse
  public :: SaveQCConditioner
  public :: ReadQCConditioner


contains


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !      ClassShortRangeOrbitals
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !> Given an irrep, returns the pointer to the corresponding symmetry-adapted short range functions
  function ClassShortRangeOrbitalsGetSymSet( ShortRangeSet, irrep ) result( SymShortRangeSet )
    class(ClassShortRangeOrbitals), target, intent(in) :: ShortRangeSet
    class(ClassIrrep)             , intent(in) :: irrep
    type(ClassShortRangeSymOrbitals), pointer  :: SymShortRangeSet
    !
    type(ClassGroup), pointer :: Group
    integer :: iIrrep
    !
    Group => Irrep%GetGroup()
    iIrrep = Group%GetIrrepIndex( irrep )
    SymShortRangeSet => ShortRangeSet%SymSRO( iIrrep )
    !
  end function ClassShortRangeOrbitalsGetSymSet





  !> Given group, initializes the short-range orbitals for all its irreps
  subroutine ClassShortRangeOrbitalsInit( ShortRangeSet, Group, Lmax, StorageDir ) 
    class(ClassShortRangeOrbitals), intent(inout) :: ShortRangeSet
    class(ClassGroup), target     , intent(in)    :: Group
    integer,                        intent(in)    :: Lmax
    character(len=*),               intent(in)    :: StorageDir
    !
    class(ClassIrrep), pointer, dimension(:) :: irrepv
    integer :: iIrrep

    ShortRangeSet%Group => Group
    irrepv => Group%GetIrrepList()
    allocate( ShortRangeSet%SymSRO( size( irrepv ) ) )
    do iIrrep = 1, size( irrepv )
       call ShortRangeSet%SymSRO( iIrrep )%init( &
            irrepv(iIrrep), &
            AddSlash(StorageDir), &
            Lmax )
    enddo
  end subroutine ClassShortRangeOrbitalsInit


  !> Given a storage directory, loads the information about the 
  !! symmetry-adapted short-range functions for all the symmetries of the group.
  subroutine ClassShortRangeOrbitalsLoad( ShortRangeSet, Storage, Lmax, UseExternalMonomials, LoadSphericalHarm ) 
    class(ClassShortRangeOrbitals), intent(inout) :: ShortRangeSet
    character(len=*)              , intent(in)    :: Storage
    integer,                        intent(in)    :: Lmax
    logical,                        intent(in)    :: UseExternalMonomials
    logical,                        intent(in)    :: LoadSphericalHarm
    !
    class(ClassIrrep), pointer, dimension(:) :: irrepv
    integer :: iIrrep
    irrepv => ShortRangeSet%Group%GetIrrepList()
    do iIrrep = 1, size( irrepv )
       call ShortRangeSet%SymSRO( iIrrep )%Load( Storage, Lmax, UseExternalMonomials, LoadSphericalHarm )
    enddo
  end subroutine ClassShortRangeOrbitalsLoad


  !> Given a storage directory, loads the information about the 
  !! symmetry-adapted short-range functions for all the symmetries of the group.
  subroutine ClassShortRangeOrbitalsShow( ShortRangeSet, unit ) 
    class(ClassShortRangeOrbitals), intent(in) :: ShortRangeSet
    integer, optional             , intent(in) :: unit
    class(ClassIrrep), pointer, dimension(:) :: irrepv
    integer :: iIrrep, outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    irrepv => ShortRangeSet%Group%GetIrrepList()
    do iIrrep = 1, size( irrepv )
       call ShortRangeSet%SymSRO( iIrrep )%show( unit )
    enddo
  end subroutine ClassShortRangeOrbitalsShow


  !> Free the short-range orbitals class
  subroutine ClassShortRangeOrbitalsFree( ShortRangeSet ) 
    class(ClassShortRangeOrbitals), intent(inout) :: ShortRangeSet
    ShortRangeSet%Group => NULL()
    if(allocated(ShortRangeSet%SymSRO))deallocate(ShortRangeSet%SymSRO)
  end subroutine ClassShortRangeOrbitalsFree


  subroutine ClassShortRangeOrbitalsFinal( ShortRangeSet ) 
    type(ClassShortRangeOrbitals) :: ShortRangeSet
    call ShortRangeSet%free()
  end subroutine ClassShortRangeOrbitalsFinal


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !      ClassShortRangeSymOrbitals
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !> Given an irrep, initializes the symmetry-adapted orbitals
  subroutine ClassShortRangeSymOrbitalsInit( SymSet, irrep, Storage, Lmax ) 
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    class(ClassIrrep), target        , intent(in)    :: irrep
    character(len=*),                  intent(in)    :: Storage
    integer,                           intent(in)    :: Lmax
    !
    character(len=:), allocatable :: FileName
    integer :: uid
    !
    SymSet%irrep => irrep
    SymSet%Lmax = Lmax
    !
    allocate( FileName, source = AddSlash(Storage)//AddSlash(DiffuseOrbDir)//DiffuseOrbInterface//SymSet%Irrep%GetName() )
    call OpenFile( FileName, uid, "read", "formatted" )
    !
    read(uid,*) 
    read(uid,*) SymSet%NMolecularOrbitals

    call SymSet%DiffuseOrbitals%init( irrep, uid, Lmax )

    close(uid)

  end subroutine ClassShortRangeSymOrbitalsInit


  !> Given a storage directory, loads the information about the 
  !! symmetry-adapted short-range functions for all the symmetries of the group.
  subroutine ClassShortRangeSymOrbitalsLoad( SymSet, Storage, Lmax, UseExternalMonomials, LoadSphericalHarm ) 
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    character(len=*)                 , intent(in)    :: Storage
    integer,                           intent(in)    :: Lmax
    logical,                           intent(in)    :: UseExternalMonomials
    logical,                           intent(in)    :: LoadSphericalHarm
    !
    character(len=:), allocatable :: FileName
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    !.. Initializes the global quantum numbers of the diffuse orbitals
    call SymSet%DiffuseOrbitals%SetIrrep( SymSet%Irrep )
    call SymSet%DiffuseOrbitals%SetLmax( SymSet%Lmax )
    !
    !.. Determines file names
    !..
    allocate( FileName, source = AddSlash(Storage)//AddSlash(DiffuseOrbDir)//DiffuseOrbInterface//SymSet%Irrep%GetName() )
    !
    !.. Opens the file
    call OpenFile( FileName, uid, "read", "formatted" )
    !
    !.. Read the part related to the molecular orbitals and global quantum numbers
    !
    !..This line that is skipped in the interface usually contains
    ! the multiplicity of the full system, but this information
    ! is provided ib other ways to the program and is not used.
    ! Maybe in the future will be better to remove this line instead.
    read(uid,*)
    read(uid,*) SymSet%NMolecularOrbitals
    !
    !.. Call a Load function for the diffuse orbitals, passing
    !   the open file unit
    !
    call SymSet%DiffuseOrbitals%Load( uid )
    !
    !
    !.. Close the file
    !
    close(uid)
    !
  end subroutine ClassShortRangeSymOrbitalsLoad


  !> Show the formatted content of the class
  subroutine ClassShortRangeSymOrbitalsShow( SymSet, unit ) 
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    integer, optional                , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    write(outunit,"(a)") "Short Range Symmetric Orbitals info :"
    call SymSet%Irrep%show(unit)
    write(outunit,"(a,i4)") "  Maximum Angular Momentum : ",SymSet%Lmax
    write(outunit,"(a,i4)") "  Number Molecular Orbitals : ",SymSet%NMolecularOrbitals
    call SymSet%DiffuseOrbitals%show( unit )
  end subroutine ClassShortRangeSymOrbitalsShow



  subroutine ClassSRSOComputeAndSaveDiffMonOverlap( SymSet, Storage )
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    character(len=*)                 , intent(in) :: Storage
    type(ClassMatrix) :: OverlapMat
    integer :: uid
    character(len=:), allocatable :: FileName
    !
    call ComputeDiffOverlap( &
         SymSet%DiffuseOrbitals%CartMonExponents, &
         SymSet%DiffuseOrbitals%GaussExponents, &
         SymSet%DiffuseOrbitals%CartesianCoefficients, &
         OverlapMat, &
         SymSet%DiffuseOrbitals%CartesianFunc  )
    !
    allocate( FileName, source = SymSet%GetDiffuseMonOvFileName(Storage) )
    !
    call OpenFile( FileName, uid, 'write', 'formatted' )
    call OverlapMat%Write( uid )
    close( uid )
    !
  end subroutine ClassSRSOComputeAndSaveDiffMonOverlap


  subroutine ClassSRSOLoadDiffMonOverlap( SymSet, Storage, OverlapMat )
    class(ClassShortRangeSymOrbitals), intent(in)  :: SymSet
    character(len=*)                 , intent(in)  :: Storage
    type(ClassMatrix)                , intent(out) :: OverlapMat
    integer :: uid
    character(len=:), allocatable :: FileName
    !
    allocate( FileName, source = SymSet%GetDiffuseMonOvFileName(Storage) )
    !
    call OpenFile( FileName, uid, 'read', 'formatted' )
    call OverlapMat%Read( uid )
    close( uid )
    !
  end subroutine ClassSRSOLoadDiffMonOverlap



  !> Gets the total number of short range orbitals: Molecular orbitals + diffuse orbitals.
  integer function ClassShortRangeSymOrbitalsGetTotNumOrbitals( SymSet ) result( NumOrb )
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    NumOrb = SymSet%NMolecularOrbitals + SymSet%DiffuseOrbitals%GetNOrbitals()
  end function ClassShortRangeSymOrbitalsGetTotNumOrbitals




  !> Gets the number of molecular orbitals (those which do not overlap with B-splines)%
  integer function ClassShortRangeSymOrbitalsGetNumMolOrbitals( SymSet ) result( NumOrb )
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    NumOrb = SymSet%NMolecularOrbitals
  end function ClassShortRangeSymOrbitalsGetNumMolOrbitals


  !> Gets the name of the file containing the overlap matrix 
  !! between the monocentric component of diffuse orbitals.
  function ClassShortRangeSymOrbitalsGetDiffuseMonOvFileName( SymSet, Storage ) result( FileName )
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    character(len=*)                 , intent(in) :: Storage
    character(len=:), allocatable :: FileName
    allocate( FileName, source = AddSlash(Storage)//AddSlash(DiffuseOrbDir)//DiffuseMonOverlapLabel//SymSet%Irrep%GetName() )
  end function ClassShortRangeSymOrbitalsGetDiffuseMonOvFileName


  !> Gets the maximum angular momentum in the expansion.
  integer function ClassShortRangeSymOrbitalsGetLmax( SymSet ) result( Lmax )
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    Lmax = SymSet%Lmax
  end function ClassShortRangeSymOrbitalsGetLmax



  !> Gets the number of diffuse orbitals.
  integer function ClassShortRangeSymOrbitalsGetNumDiffuseOrbitals( SymSet ) result( NumOrb )
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    NumOrb = SymSet%DiffuseOrbitals%GetNOrbitals()
  end function ClassShortRangeSymOrbitalsGetNumDiffuseOrbitals



  !> Computes the transformation matrices between cartesian Gaussians and symmetry adapted spherical harmonics. 
  subroutine ClassShortRangeSymOrbitalsComputeTransfMat( SymSet )
    !
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    !
    call SymSet%DiffuseOrbitals%ComputeTransfMat()
    !
  end subroutine ClassShortRangeSymOrbitalsComputeTransfMat



  subroutine ClassShortRangeSymOrbitalsGetTransfMat( SymSet, Identifier, Array )
    !
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    character(len=*),                  intent(in)    :: Identifier
    complex(kind(1d0)), allocatable,   intent(out)   :: Array(:,:)
    !
    call SymSet%DiffuseOrbitals%GetTransfMat( Identifier, Array )
    !
  end subroutine ClassShortRangeSymOrbitalsGetTransfMat




  function ClassShortRangeSymOrbitalsGetCartesianSet( SymSet ) result( CartSet )
    !
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    type(ClassCartesianSymmetricSet), pointer :: CartSet
    !
    CartSet => SymSet%DiffuseOrbitals%GetCartesianSet()
    !
  end function ClassShortRangeSymOrbitalsGetCartesianSet


  subroutine ClassShortRangeSymOrbitalsAcquireInterface( SymSet, StorageDir, InterfaceFile )
    !
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    character(len=*), intent(in) :: StorageDir
    character(len=*), intent(in) :: InterfaceFile
    !
    character(len=:), allocatable :: Dir
    !
    allocate( Dir, source = AddSlash(StorageDir)//AddSlash(DiffuseOrbDir) )
    call execute_command_line("mkdir -p "//Dir )
    call execute_command_line("cp "//InterfaceFile//" "//Dir )
    !
  end subroutine ClassShortRangeSymOrbitalsAcquireInterface


  !> Computes the full transformation matrice between cartesian Gaussians multiplied by ll the used exponentials and symmetry adapted spherical harmonics also multiplied by the exponentials, to pass from Xlm to cartesian representation. 
  subroutine ClassShortRangeSymOrbitalsComputeFullTransfMat( SymSet, NumExp )
    !
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    integer,                           intent(in)    :: NumExp
    !
    call SymSet%DiffuseOrbitals%ComputeFullTransfMat( NumExp )
    !
  end subroutine ClassShortRangeSymOrbitalsComputeFullTransfMat



  !> Fetches the transformation matrix that expresses the cartesian monomials in terms of the symmetric adapted spherical harmonics Xlm%
  subroutine ClassShortRangeSymOrbitalsFetchXlmToCartesianMatrix( SymSet, Mat )
    !
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    complex(kind(1d0)), allocatable,   intent(out)   :: Mat(:,:)
    !
    call SymSet%DiffuseOrbitals%FetchXlmToCartesianMatrix( Mat )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchXlmToCartesianMatrix



  !> Fetches the full transformation matrix that expresses the cartesian monomials in terms of the symmetric adapted spherical harmonics Xlm, where both sets are multiplied by the exponentials functions.
  subroutine ClassShortRangeSymOrbitalsFetchFullTransfMat( SymSet, Mat )
    !
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    complex(kind(1d0)), allocatable,   intent(out)   :: Mat(:,:)
    !
    call SymSet%DiffuseOrbitals%FetchFullTransfMat( Mat )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchFullTransfMat



  !> Fetches the array of cartesian monomials as appear in the interface file.
  subroutine ClassShortRangeSymOrbitalsFetchCartMonInterface( SymSet, MonArray )
    !
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    integer, allocatable,              intent(out)   :: MonArray(:,:)
    !
    call SymSet%DiffuseOrbitals%FetchCartMonInterface( MonArray )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchCartMonInterface



  !> Fetches the array of cartesian monomials as it was initialized in ClassCartesianSymSet%
  subroutine ClassShortRangeSymOrbitalsFetchCartMonInitialized( SymSet, MonArray )
    !
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    integer, allocatable,              intent(out)   :: MonArray(:,:)
    !
    call SymSet%DiffuseOrbitals%FetchCartMonInitialized( MonArray )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchCartMonInitialized



  !> Fetches the array of symmetry adapted spherical harmonics's indexes.
  subroutine ClassShortRangeSymOrbitalsFetchXlmTransfIndexes( SymSet, IndArray )
    !
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    integer, allocatable,              intent(out)   :: IndArray(:,:)
    !
    call SymSet%DiffuseOrbitals%FetchXlmTransfIndexes( IndArray )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchXlmTransfIndexes



  !> Fetches the diffuse orbitals expansion in terms of the cartesian Gaussians.
  subroutine ClassShortRangeSymOrbitalsFetchCartesianCoeff( SymSet, Array )
    !
    class(ClassShortRangeSymOrbitals), intent(in) :: SymSet
    real(kind(1d0)), allocatable,      intent(out)   :: Array(:,:)
    !
    call SymSet%DiffuseOrbitals%FetchCartesianCoeff( Array )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchCartesianCoeff


  !> Free the short-range orbitals class
  subroutine ClassShortRangeSymOrbitalsFree( SymSet ) 
    class(ClassShortRangeSymOrbitals), intent(inout) :: SymSet
    SymSet%irrep => NULL()
    SymSet%Lmax = -1
    SymSet%NMolecularOrbitals = 0
    call SymSet%DiffuseOrbitals%free()
  end subroutine ClassShortRangeSymOrbitalsFree


  subroutine ClassShortRangeSymOrbitalsFinal( SymSet ) 
    type(ClassShortRangeSymOrbitals) :: SymSet
    call SymSet%free()
  end subroutine ClassShortRangeSymOrbitalsFinal



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !      ClassDiffuseOrbitals
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



  !> Given an irrep, initializes the symmetry-adapted orbitals
  subroutine ClassDiffuseOrbitalsInit( DifSet, Irrep, uid, Lmax ) 
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    class(ClassIrrep), target  , intent(in)    :: Irrep
    integer,                     intent(in)    :: uid
    integer,                     intent(in)    :: Lmax
    !
    logical :: USE_CARTESIAN_FUNCTIONS = .FALSE.

    DifSet%Irrep         => Irrep
    DifSet%Lmax          =  Lmax
    DifSet%CartesianFunc = USE_CARTESIAN_FUNCTIONS

    call DifSet%Load( uid )
    call DifSet%ComputeTransfMat()

  end subroutine ClassDiffuseOrbitalsInit


  !> Given a storage directory, loads the information about the 
  !! symmetry-adapted short-range functions for all the symmetries of the group.
  subroutine ClassDiffuseOrbitalsLoad( DifSet, unit ) 
    use ModuleAngularMomentum
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    integer                    , intent(in)    :: unit
    !
    integer :: NumCartesians, iCartCoord, iCartMon, Lmax, i, j, NumGaussExp, Ntot, k, l
    integer :: MonomialExponents(3,NumMaxExponents)
    real(kind(1d0)), allocatable :: CartGaussNormFactor(:)
    real(kind(1d0)) :: VecExponents(NumMaxExponents)
    real(kind(1d0)) :: Prefactor
    real(kind(1d0)) :: ExtraFactor
    real(kind(1d0)), external :: DGamma
    real(kind(1d0)), allocatable :: RealVecExponents(:)
    type(ClassMatrix) :: OverlapMat, DiffPrecMat
    real(kind(1d0)) :: AuxVal
    !
    character(len=:), allocatable :: PlotDir, PlotOrbFile, PlotAllFile
    integer :: nOrb, nEval, nExp, nMon, uidOrb, uidAll
    real(kind(1d0)) :: EvalR, Val
    real(kind(1d0)), parameter :: Rmin = 1.d-4
    real(kind(1d0)), parameter :: Rmax = 400.d0
    integer, parameter :: NumPoints = 1000
    !
    !.. Read Lmax
    !
    read(unit,*) Lmax
    !
    if ( Lmax /= DifSet%Lmax ) then
       call Assert( "The maximum angular momentum read from "//&
               "the diffuse orbitals configuraion file differs from"//&
               "the value used in the diffuse orbital initialization." )
    end if
    !
    !
    if ( .not.DifSet%CartesianFunc ) then
       !.. Read set of spherical harmonics
       !
       read(unit,*) NumCartesians,( (MonomialExponents(i,j),i=1,3),j=1,NumCartesians)
       !
       do i = 1, NumCartesians
          k = MonomialExponents(1,i)
          l = MonomialExponents(2,i)
          MonomialExponents(1,i) = l + 2*k
       end do
       !
       allocate( DifSet%CartMonExponents, source = MonomialExponents(:,:NumCartesians) )
!!$       !
    else
       !.. Read set of cartesians 
       !
       read(unit,*) NumCartesians,( (MonomialExponents(i,j),i=1,3),j=1,NumCartesians)
       !
       allocate( DifSet%CartMonExponents, source = MonomialExponents(:,:NumCartesians) )
       !
    end if
    !
    !
    if ( .not.DifSet%CartesianFunc ) then
       !
       call DifSet%SpherSymSet%Init( DifSet%Lmax, &
            DifSet%Irrep, &
            DifSet%CartMonExponents )
       !
    else
       !
       call DifSet%CartesianSymSet%Init( DifSet%Lmax, &
            DifSet%Irrep, &
            DifSet%CartMonExponents )
       !
    end if
    !
    !
    ! Reads the Gaussian exponents.
    do i = 0, Lmax
       read(unit,*) l, NumGaussExp, (VecExponents(j),j=1,NumGaussExp)
    end do
    allocate( RealVecExponents, source = VecExponents(:NumGaussExp ) )
    if ( allocated(DifSet%GaussExponents) ) deallocate( DifSet%GaussExponents )
    allocate( DifSet%GaussExponents, source = RealVecExponents )
    !
    read(unit,*) DifSet%NCartesianGaussians
    read(unit,*) DifSet%NOrbitals
    !
    !.. Read coefficients
    !
    if ( allocated(DifSet%CartesianCoefficients) ) deallocate( DifSet%CartesianCoefficients )
    allocate( DifSet%CartesianCoefficients(DifSet%NCartesianGaussians,DifSet%NOrbitals) )
    !
    do j = 1, DifSet%NOrbitals
       read(unit,*) ( DifSet%CartesianCoefficients(i,j), i=1,DifSet%NCartesianGaussians )
    end do
    !
    !
    !
    ! Computes the normalization factor (NF) of the cartesian gaussian, as is expected to be carried out by the QC program: NF = ... If activated lines:
    ! Prefactor = sqrt(2.d0**(dble(Ntot))/(Pi**(1.5d0)))
    ! CartGaussNormFactor( (i-1)*NumGaussExp + j ) = &
    !           PreFactor * &
    !           exp( ( dble(Ntot) + 3.d0 / 2.d0 ) / 2.d0 * log( 2.d0*VecExponents( j ) ) )
    !
    if ( allocated(CartGaussNormFactor) ) deallocate( CartGaussNormFactor )
    allocate( CartGaussNormFactor( DifSet%NCartesianGaussians ) )
    CartGaussNormFactor = 0.d0
    !
    do i = 1, NumCartesians
       !
       if ( DifSet%CartesianFunc ) then
          !
          Ntot = sum(MonomialExponents(1:3,i))
          Prefactor = sqrt(2.d0**(dble(Ntot))/(Pi**(1.5d0)))
          !
       else
          !
          Ntot = MonomialExponents(1,i)
          PreFactor =sqrt(2.d0/DGamma( dble(Ntot) + 1.5d0 ))
          !
          !*** Possible QC normalization for Xlm%
          PreFactor = PreFactor * ((-1.d0)**MonomialExponents(3,i))
          !***
          !
       end if
       !
       do j = 1, NumGaussExp
          !
          l = l + 1
          !
          CartGaussNormFactor( (i-1)*NumGaussExp + j ) = &
               PreFactor *  &
               exp( ( dble(Ntot) + 3.d0 / 2.d0 ) / 2.d0 * log( 2.d0*VecExponents( j ) ) )
          !
       end do
       !
    end do
    !
    !
    !***
    ! Plot the diffuse orbitals part expresed exclusively by
    ! the mono-centric Gaussian in the 'z' direction.
    if (PlotDiffuseOrb ) then
       !
       write(OUTPUT_UNIT,*) "Plotting the diffuse orbitals in 'z' direction ..."
       !
       allocate( PlotDir, source = "PlotNullSpace/" )
       allocate( PlotAllFile, source = PlotDir//"OrbAll" )
       call OpenFile( PlotAllFile, uidAll, "write", "formatted" )
       !
       do nOrb = 1, DifSet%NOrbitals
          !
          write(OUTPUT_UNIT,*) "Plotting diffuse orbital "//AlphabeticNumber(nOrb)//" ..."
          allocate( PlotOrbFile, source = PlotDir//"Orb"//AlphabeticNumber(nOrb) )
          call OpenFile( PlotOrbFile, uidOrb, "write", "formatted" )
          !
          do nEval = 1, NumPoints
             !
             l = 0
             Val = 0.d0
             EvalR = Rmin + dble(nEval-1)/dble(NumPoints-1)*(Rmax-Rmin)
             !
             do nMon = 1, NumCartesians
                !
                Ntot = MonomialExponents(1,nMon)
                !
                do nExp = 1, NumGaussExp
                   !
                   l = l + 1
                   !
                   Val = Val + EvalR**(Ntot)*exp(-VecExponents(nExp)*EvalR**2)*&
                        dble(Ylm(0.d0,1.d0,MonomialExponents(2,nMon),MonomialExponents(3,nMon)))*&
                        DifSet%CartesianCoefficients(l,nOrb)
                   !
                end do
                !
             end do
             !
             write(uidOrb,*) EvalR, Val
             write(uidAll,*) EvalR, Val
             !
          end do
          !
          close( uidOrb )
          deallocate( PlotOrbFile )
          write(uidAll,*)
          !
       end do
       !
       close( uidAll )
       deallocate( PlotAllFile )
       !
       stop
       !
    end if
    !
    l = 0
    do i = 1, NumCartesians
       do j = 1, NumGaussExp
          !
          l = l + 1 
          DifSet%CartesianCoefficients(l,:) = DifSet%CartesianCoefficients(l,:) * CartGaussNormFactor(l)
       end do
    end do
    !
  end subroutine ClassDiffuseOrbitalsLoad



  !> Set the diffuse orbitals's Irrep%
  subroutine  ClassDiffuseOrbitalsSetIrrep( DifSet, Irrep )
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    class(ClassIrrep), target,   intent(in)    :: Irrep
    !
    DifSet%Irrep => Irrep
    !
  end subroutine ClassDiffuseOrbitalsSetIrrep


  !> Set the diffuse orbitals's maximum angular momentum.
  subroutine  ClassDiffuseOrbitalsSetLmax( DifSet, Lmax )
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    integer,                     intent(in)    :: Lmax
    !
    DifSet%Lmax = Lmax
    !
  end subroutine ClassDiffuseOrbitalsSetLmax


  !> Gets the number of diffuse orbitals.
  integer function ClassDiffuseOrbitalsGetNOrbitals( DifSet ) result( NumOrb )
    !
    class(ClassDiffuseOrbitals), intent(in) :: DifSet
    !
    if ( .not. allocated(DifSet%CartesianCoefficients) ) then
       call Assert( "The cartesian coefficients has not been allocated, imposible to get the number of molecular orbitals." )
    end if
    !
    NumOrb = DifSet%NOrbitals
    !
  end function ClassDiffuseOrbitalsGetNOrbitals



  !> Computes the transformation matrices between cartesian (or spherical) Gaussians and symmetry adapted spherical harmonics. 
  subroutine ClassDiffuseOrbitalsComputeTransfMat( DifSet )
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    !
    if ( DifSet%CartesianFunc ) then
       call DifSet%CartesianSymSet%ComputeTransfMat()
    else
       call DifSet%SpherSymSet%SpherComputeTransfMat()
    end if
    !
  end subroutine ClassDiffuseOrbitalsComputeTransfMat



  subroutine ClassDiffuseOrbitalsGetTransfMat( DifSet, Identifier, Array )
    !
    class(ClassDiffuseOrbitals),     intent(inout) :: DifSet
    character(len=*),                intent(in)    :: Identifier
    complex(kind(1d0)), allocatable, intent(out)   :: Array(:,:)
    !
    if ( DifSet%CartesianFunc ) then
       !
       if ( Identifier .is. "XlmToCartesian" ) then
          call DifSet%CartesianSymSet%FetchXlmToCartesianMatrix( Array )
       elseif ( Identifier .is. "CartesianToXlm" ) then
          call DifSet%CartesianSymSet%FetchCartesianToXlmMatrix( Array )
       else
          call Assert( "Invalid identifier: "//Identifier//" , it should be 'XlmToCartesian' or 'CartesianToXlm'." )
       end if
       !
    else
       !
       if ( Identifier .is. "XlmToSpher" ) then
          call DifSet%SpherSymSet%FetchXlmToSpherMatrix( Array )
       elseif ( Identifier .is. "SpherToXlm" ) then
          call DifSet%SpherSymSet%FetchSpherToXlmMatrix( Array )
       else
          call Assert( "Invalid identifier: "//Identifier//" , it should be 'XlmToSpher' or 'SpherToXlm'." )
       end if
       !
    end if
    !
  end subroutine ClassDiffuseOrbitalsGetTransfMat




  function ClassDiffuseOrbitalsGetCartesianSet( DifSet ) result( CartSet )
    !
    class(ClassDiffuseOrbitals)      , intent(inout) :: DifSet
    type (ClassCartesianSymmetricSet), pointer       :: CartSet
    type (ClassCartesianSymmetricSet), pointer       :: CopyCartSet
    !
    CartSet => NULL()
    allocate(CopyCartSet)
    CopyCartSet = DifSet%CartesianSymSet
    CartSet => CopyCartSet
    CopyCartSet => NULL()
    !
  end function ClassDiffuseOrbitalsGetCartesianSet



  function ClassDiffuseOrbitalsGetSpherSet( DifSet ) result( SpherSet )
    !
    class(ClassDiffuseOrbitals)  , intent(inout) :: DifSet
    type (ClassSpherSymmetricSet), pointer       :: SpherSet
    type (ClassSpherSymmetricSet), pointer       :: CopySpherSet
    !
    SpherSet => NULL()
    allocate(CopySpherSet) 
    CopySpherSet = DifSet%SpherSymSet
    SpherSet => CopySpherSet
    CopySpherSet => NULL()
    !
  end function ClassDiffuseOrbitalsGetSpherSet




  !> Computes the full transformation matrix between cartesian Gaussians multiplied by the used exponentials and symmetry adapted spherical harmonics also multiplied by the exponentials, to pass from Xlm to cartesian representation.
  subroutine ClassDiffuseOrbitalsComputeFullXlmToCartTransfMat( DifSet, NumExp )
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    integer,                     intent(in)    :: NumExp
    !
    complex(kind(1d0)), allocatable :: XlmToCartArray(:,:)
    integer :: NumOrigRows, NewNumRows, iOrig, jOrig, iExp
    !
    if ( .not.DifSet%CartesianFunc ) then
       call Assert( "Impossible to compute the Xlm to Cartesian full matrix because the Gaussian functions are not cartesian." )
    end if
    !
    call DifSet%CartesianSymSet%FetchXlmToCartesianMatrix( XlmToCartArray )
    !
    NumOrigRows = size(XlmToCartArray,1)
    NewNumRows  = NumExp * NumOrigRows
    !
    if ( allocated(DifSet%FullXlmToCartArray) ) deallocate( DifSet%FullXlmToCartArray )
    allocate( DifSet%FullXlmToCartArray(NewNumRows,NewNumRows) )
    DifSet%FullXlmToCartArray = Z0
    !
    if ( .true. ) then !*** Incredibly this way is avoided the fatal compilation error.
       do iOrig = 1, NumOrigRows
          do jOrig = 1, NumOrigRows
             do iExp = 1, NumExp
                !
                DifSet%FullXlmToCartArray( &
                     NumExp*(iOrig-1) + iExp, &
                     NumExp*(jOrig-1) + iExp ) = XlmToCartArray(iOrig,jOrig)
                !
             enddo
          enddo
       enddo
    end if
    !
  end subroutine ClassDiffuseOrbitalsComputeFullXlmToCartTransfMat



  !> Computes the full transformation matrix between spherical Gaussians multiplied by the used exponentials and symmetry adapted spherical harmonics also multiplied by the exponentials, to pass from Xlm to cartesian representation.
  subroutine ClassDiffuseOrbitalsComputeFullXlmToSpherTransfMat( DifSet, NumExp )
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    integer,                     intent(in)    :: NumExp
    !
    complex(kind(1d0)), allocatable :: XlmToSpherArray(:,:)
    integer :: NumOrigRows, NewNumRows, iOrig, jOrig, iExp, i, j
    !
    if ( DifSet%CartesianFunc ) then
       call Assert( "Impossible to compute the Xlm to Spherical full matrix because the Gaussian functions are not spherical." )
    end if
    !
    call DifSet%SpherSymSet%FetchXlmToSpherMatrix( XlmToSpherArray )
    !
    NumOrigRows = size(XlmToSpherArray,1)
    NewNumRows  = NumExp * NumOrigRows
    !
    if ( allocated(DifSet%FullXlmToCartArray) ) deallocate( DifSet%FullXlmToCartArray )
    allocate( DifSet%FullXlmToCartArray(NewNumRows,NewNumRows) )
    DifSet%FullXlmToCartArray = Z0
    !
    do iOrig = 1, NumOrigRows
       do jOrig = 1, NumOrigRows
          do iExp = 1, NumExp
             !
             DifSet%FullXlmToCartArray( &
                  NumExp*(iOrig-1) + iExp, &
                  NumExp*(jOrig-1) + iExp ) = XlmToSpherArray(iOrig,jOrig)
             !
          enddo
       enddo
    enddo
    !
  end subroutine ClassDiffuseOrbitalsComputeFullXlmToSpherTransfMat



  subroutine ClassDiffuseOrbitalsComputeFullTransfMat( DifSet, NumExp )
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    integer,                     intent(in)    :: NumExp
    !
    if ( DifSet%CartesianFunc ) then
       call DifSet%ComputeFullXlmToCartTransfMat( NumExp )
    else
       call DifSet%ComputeFullXlmToSpherTransfMat( NumExp )
    end if
    !
  end subroutine ClassDiffuseOrbitalsComputeFullTransfMat


  !> Fetches the transformation matrix that expresses the cartesian monomials in terms of the symmetric adapted spherical harmonics Xlm%
  subroutine ClassDiffuseOrbitalsFetchXlmToCartesianMatrix( DifSet, Mat )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    complex(kind(1d0)), allocatable, intent(out)   :: Mat(:,:)
    !
    call DifSet%CartesianSymSet%FetchXlmToCartesianMatrix( Mat )
    !
  end subroutine ClassDiffuseOrbitalsFetchXlmToCartesianMatrix


  !> Fetches the transformation matrix that expresses the spherical monomials in terms of the symmetric adapted spherical harmonics Xlm%
  subroutine ClassDiffuseOrbitalsFetchXlmToSpherMatrix( DifSet, Mat )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    complex(kind(1d0)), allocatable, intent(out)   :: Mat(:,:)
    !
    call DifSet%SpherSymSet%FetchXlmToSpherMatrix( Mat )
    !
  end subroutine ClassDiffuseOrbitalsFetchXlmToSpherMatrix



  !> Fetches the transformation matrix that expresses the spherical/cartesian monomials in terms of the symmetric adapted spherical harmonics Xlm%
  subroutine ClassDiffuseOrbitalsFetchTransfMatrix( DifSet, Mat )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    complex(kind(1d0)), allocatable, intent(out)   :: Mat(:,:)
    !
    if ( DifSet%CartesianFunc ) then
       call DifSet%FetchXlmToCartesianMatrix( Mat )
    else
       call DifSet%FetchXlmToSpherMatrix( Mat )
    end if
    !
  end subroutine ClassDiffuseOrbitalsFetchTransfMatrix



  !> Fetches the full transformation matrix that expresses the cartesian monomials in terms of the symmetric adapted spherical harmonics Xlm, where both sets are multiplied by the exponentials functions.
  subroutine ClassDiffuseOrbitalsFetchFullXlmToCartesianTransfMat( DifSet, Mat )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    complex(kind(1d0)), allocatable, intent(out)   :: Mat(:,:)
    !
    if ( .not. allocated(DifSet%FullXlmToCartArray) ) then
       call Assert( "Imposible to fetch the full transformation matrix because it has not been computed yet." )
    end if
    !
    allocate ( Mat, source = DifSet%FullXlmToCartArray )
    !
  end subroutine ClassDiffuseOrbitalsFetchFullXlmToCartesianTransfMat



  !> Fetches the full transformation matrix that expresses the spherical monomials in terms of the symmetric adapted spherical harmonics Xlm, where both sets are multiplied by the exponentials functions.
  subroutine ClassDiffuseOrbitalsFetchFullXlmToSpherTransfMat( DifSet, Mat )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    complex(kind(1d0)), allocatable, intent(out)   :: Mat(:,:)
    !
    if ( .not. allocated(DifSet%FullXlmToCartArray) ) then
       call Assert( "Imposible to fetch the full transformation matrix because it has not been computed yet." )
    end if
    !
    allocate ( Mat, source = DifSet%FullXlmToCartArray )
    !
  end subroutine ClassDiffuseOrbitalsFetchFullXlmToSpherTransfMat



  subroutine ClassDiffuseOrbitalsFetchFullTransfMat( DifSet, Mat )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    complex(kind(1d0)), allocatable, intent(out)   :: Mat(:,:)
    !
    if ( DifSet%CartesianFunc ) then
       call DifSet%FetchFullXlmToCartesianTransfMat( Mat )
    else
       call DifSet%FetchFullXlmToSpherTransfMat( Mat )
    end if
    !
  end subroutine ClassDiffuseOrbitalsFetchFullTransfMat



  !> Fetches the array of monomial orbitals as appear in the interface file.
  subroutine ClassDiffuseOrbitalsFetchCartMonInterface( DifSet, MonArray )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    integer, allocatable,            intent(out)   :: MonArray(:,:)
    !
    if ( .not. allocated(DifSet%CartMonExponents) ) then
       call Assert( "Imposible to fetch the interface file's array of cartesian monomials because it has not been allocated." )
    end if
    !
    allocate( MonArray, source = DifSet%CartMonExponents )
    !
  end subroutine ClassDiffuseOrbitalsFetchCartMonInterface



  !> Fetches the array of monomial orbitals as was initialized in ClassCartesianSymSet%
  subroutine ClassDiffuseOrbitalsFetchCartMonInitialized( DifSet, MonArray )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    integer, allocatable,            intent(out)   :: MonArray(:,:)
    !
    if ( DifSet%CartesianFunc ) then
       !
       if ( .not. DifSet%CartesianSymSet%Initialized() ) then
          call Assert( "Imposible to fetch the initialized array of cartesian monomials because it has not been allocated." )
       end if
       !
       call DifSet%CartesianSymSet%FetchMonomialExponents( MonArray )
       !
    else
       !
       if ( .not. DifSet%SpherSymSet%Initialized() ) then
          call Assert( "Imposible to fetch the initialized array of spherical monomials because it has not been allocated." )
       end if
       !
       call DifSet%SpherSymSet%FetchMonomialExponents( MonArray )
       !
    end if
    !
  end subroutine ClassDiffuseOrbitalsFetchCartMonInitialized




  !> Fetches the array of symmetry adapted spherical harmonics's indexes.
  subroutine ClassDiffuseOrbitalsFetchXlmTransfIndexes( DifSet, IndArray )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    integer, allocatable,            intent(out)   :: IndArray(:,:)
    !
    if ( DifSet%CartesianFunc ) then
       !
       if ( .not. DifSet%CartesianSymSet%Initialized() ) then
          call Assert( "Imposible to fetch the symmetry adapted spherical harmonics's indexes because it has not been allocated." )
       end if
       !
       call DifSet%CartesianSymSet%FetchXlmTransfIndexes( IndArray )
       !
    else
       !
       if ( .not. DifSet%SpherSymSet%Initialized() ) then
          call Assert( "Imposible to fetch the symmetry adapted spherical harmonics's indexes because it has not been allocated." )
       end if
       !
       call DifSet%SpherSymSet%SpherFetchXlmTransfIndexes( IndArray )
       !
    end if
    !
  end subroutine ClassDiffuseOrbitalsFetchXlmTransfIndexes



  !> Fetches the diffuse orbitals expansion in terms of the cartesian Gaussians.
  subroutine ClassDiffuseOrbitalsFetchCartesianCoeff( DifSet, Array )
    !
    class(ClassDiffuseOrbitals),     intent(in) :: DifSet
    real(kind(1d0)), allocatable,    intent(out)   :: Array(:,:)
    !
    if ( .not. allocated(DifSet%CartesianCoefficients) ) then
       call Assert( "Imposible to fetch the diffuse orbitals coefficients because they have not been allocated." )
    end if
    !
    allocate( Array, source = DifSet%CartesianCoefficients )
    !
  end subroutine ClassDiffuseOrbitalsFetchCartesianCoeff


  !> Show the formatted content of the class
  subroutine ClassDiffuseOrbitalsShow( DifSet, unit ) 
    !
    class(ClassDiffuseOrbitals), intent(in) :: DifSet
    integer, optional          , intent(in) :: unit
    !
    integer :: outunit
    !
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    !
    write(outunit,"(a)") "Diffuse Orbitals Info : "
    !
    call DifSet%Irrep%Show()
    !
    write(outunit,"(a,i4)") "  Maximum Angular Momentum : ",DifSet%Lmax
    !
    if ( DifSet%CartesianFunc ) then
       write(outunit,"(a,a)") "  Type of Gaussian : ", "Cartesian"
    else
       write(outunit,"(a,a)") "  Type of Gaussian : ", "Spherical"
    end if
    !
    write(outunit,"(a,i4)") "  Number Diffuse Orbitals : ",DifSet%NOrbitals
    write(outunit,"(a,i4)") "  Number Cartesian Gaussians : ",DifSet%NCartesianGaussians
    write(outunit,"(a)") " Cartesian Monomial Exponents : "
    write(outunit,*) DifSet%CartMonExponents
    !
    if ( DifSet%CartesianFunc ) then
       call DifSet%CartesianSymSet%Show()
    else
       call DifSet%SpherSymSet%Show()
    end if
    !
  end subroutine ClassDiffuseOrbitalsShow


  !> Free the short-range orbitals class
  subroutine ClassDiffuseOrbitalsFree( DifSet ) 
    !
    class(ClassDiffuseOrbitals), intent(inout) :: DifSet
    !
    DifSet%Irrep => NULL()
    call DifSet%CartesianSymSet%Free()
    !
    DifSet%Lmax = -1
    DifSet%NOrbitals = 0
    DifSet%NCartesianGaussians = 0
    !
    if ( allocated(DifSet%CartMonExponents) ) deallocate( DifSet%CartMonExponents )
    if ( allocated(DifSet%GaussExponents) ) deallocate( DifSet%GaussExponents )
    if ( allocated(DifSet%CartesianCoefficients) ) deallocate( DifSet%CartesianCoefficients )
    if ( allocated(DifSet%FullXlmToCartArray) ) deallocate( DifSet%FullXlmToCartArray )
    !
  end subroutine ClassDiffuseOrbitalsFree


  subroutine ClassDiffuseOrbitalsFinal( DifSet ) 
    type(ClassDiffuseOrbitals) :: DifSet
    call DifSet%free()
  end subroutine ClassDiffuseOrbitalsFinal



  subroutine ComputeDiffOverlap( CartMonExponents, VecExponents, CoeffMat, Overlap, CartesianGaussian )
    !
    integer,                      intent(in)  :: CartMonExponents(:,:)
    real(kind(1d0)),              intent(in)  :: VecExponents(:)
    real(kind(1d0)),              intent(in)  :: CoeffMat(:,:)
    class(ClassMatrix),           intent(out) :: Overlap
    logical,                      intent(in)  :: CartesianGaussian
    !
    integer :: i, j, NumDiff, uid
    real(kind(1d0)), allocatable :: OverlapMat(:,:)
    character(len=:), allocatable :: File
    logical, parameter :: PrintMat =.true.
    type(ClassMatrix) :: CopyCoeffMat, PrimitiveMat
    !
    CopyCoeffMat = CoeffMat
    Overlap      =  CoeffMat
    call Overlap%Transpose()
    !
    call ComputePrimitiveGaussianOverlap( CartMonExponents, VecExponents, CartesianGaussian, PrimitiveMat )
    !
    call Overlap%Multiply( PrimitiveMat, 'Right', 'N' )
    call PrimitiveMat%Free()
    call Overlap%Multiply( CopyCoeffMat, 'Right', 'N' )
    call CopyCoeffMat%Free()

!!$    NumDiff = size(CoeffMat,2)
!!$    !
!!$    allocate( OverlapMat(NumDiff,NumDiff) )
!!$    OverlapMat = 0.d0
!!$    call Overlap%InitFull( NumDiff, NumDiff )
!!$    !
!!$    do j = 1, NumDiff
!!$       do i = 1, j
!!$          !
!!$          OverlapMat(i,j) = IntegralDiffDiff( CoeffMat(:,i), CoeffMat(:,j), CartMonExponents, VecExponents, CartesianGaussian )
!!$          OverlapMat(j,i) = OverlapMat(i,j)
!!$          !
!!$       end do
!!$    end do
!!$    !
!!$    Overlap = OverlapMat
!!$    deallocate( OverlapMat )
    !
  end subroutine ComputeDiffOverlap



  subroutine ComputePrimitiveGaussianOverlap( CartMonExponents, VecExponents, CartesianGaussian, OutMat )
    !
    integer,            intent(in)  :: CartMonExponents(:,:)
    real(kind(1d0)),    intent(in)  :: VecExponents(:)
    logical,            intent(in)  :: CartesianGaussian
    class(ClassMatrix), intent(out) :: OutMat
    !
    integer :: i, j, k, l, m, n, NumCartesian, NumMonomial, NumExp
    real(kind(1d0)) :: Val
    !
    NumMonomial = size(CartMonExponents,2)
    NumExp = size(VecExponents)
    NumCartesian = NumMonomial * NumExp
    !
    call OutMat%InitFull( NumCartesian, NumCartesian )
    !
    k = 0
    do i = 1, NumMonomial
       do j = 1, NumExp
          k = k + 1
          !
          n = 0
          do l = 1, NumMonomial
             do m = 1, NumExp
                n = n + 1
                !
                if ( n <= k ) then
                   !
                   if ( CartesianGaussian ) then
                      !
                      Val = IntegralCartCart(     &
                           CartMonExponents(:,l), &
                           VecExponents(m)      , &
                           CartMonExponents(:,i), &
                           VecExponents(j)  )
                      !
                   else
                      !
                      Val = IntegralSpherSpher(   &
                           CartMonExponents(:,l), &
                           VecExponents(m)      , &
                           CartMonExponents(:,i), &
                           VecExponents(j)    )
                      !
                   end if
                   !
                   call OutMat%SetElement( n, k, Val )
                   call OutMat%SetElement( k, n, Val )
                   !
                end if
                !
             end do
          end do
          !
       end do
    end do
    !
  end subroutine ComputePrimitiveGaussianOverlap



  real(kind(1d0)) function IntegralDiffDiff( BraVec, KetVec, CartMonExponents, VecExponents, CartesianGaussian ) result( Integral )
    !
    real(kind(1d0)), intent(in) :: BraVec(:)
    real(kind(1d0)), intent(in) :: KetVec(:)
    integer,         intent(in) :: CartMonExponents(:,:)
    real(kind(1d0)), intent(in) :: VecExponents(:)
    logical,         intent(in) :: CartesianGaussian
    !
    integer :: i, j, k, l, m, n, NumCartesian, NumMonomial, NumExp
    !
    NumCartesian = size(BraVec)
    NumMonomial = size(CartMonExponents,2)
    NumExp = size(VecExponents)
    !
    if ( NumCartesian /= NumMonomial*NumExp ) call Assert( &
         'The total number of cartesian Gaussian do not conform.' )
    !
    Integral = 0.d0
    k = 0
    do i = 1, NumMonomial
       do j = 1, NumExp
          k = k + 1
          !
          n = 0
          do l = 1, NumMonomial
             do m = 1, NumExp
                n = n + 1
                !
                if ( CartesianGaussian ) then
                   !
                   Integral = Integral + BraVec(n)*KetVec(k)*&
                        IntegralCartCart( &
                        CartMonExponents(:,l), &
                        VecExponents(m), &
                        CartMonExponents(:,i), &
                        VecExponents(j) )
                   !
                else
                   !
                   Integral = Integral + BraVec(n)*KetVec(k)*&
                        IntegralSpherSpher( &
                        CartMonExponents(:,l), &
                        VecExponents(m), &
                        CartMonExponents(:,i), &
                        VecExponents(j) )
                   !
                end if
                !
             end do
          end do
          !
       end do
    end do
    !
  end function IntegralDiffDiff




  real(kind(1d0)) function IntegralCartCart( BraMonExp, BraVecExp, KetMonExp, KetVecExp ) result( Integral )
    !
    integer,         intent(in) :: BraMonExp(:)
    real(kind(1d0)), intent(in) :: BraVecExp
    integer,         intent(in) :: KetMonExp(:)
    real(kind(1d0)), intent(in) :: KetVecExp
    !
    integer :: L, i, DF
    !
    Integral = 0.d0
    !
    do i = 1, 3
       if ( mod(BraMonExp(i)+KetMonExp(i),2) /= 0 ) return
    end do
    !
    L = 0
    DF = 1
    do i = 1, 3
       L = L + (BraMonExp(i)+KetMonExp(i))/2
       DF = DF * DoubleFact(BraMonExp(i)+KetMonExp(i)-1)
    end do
    !
    Integral = PI**(1.5d0)/((BraVecExp+KetVecExp)**(dble(L)+1.5d0)*2.d0**dble(L)) * dble(DF)
    !
  end function IntegralCartCart





  real(kind(1d0)) function IntegralSpherSpher( BraMonExp, BraVecExp, KetMonExp, KetVecExp ) result( Integral )
    !
    integer,         intent(in) :: BraMonExp(:)
    real(kind(1d0)), intent(in) :: BraVecExp
    integer,         intent(in) :: KetMonExp(:)
    real(kind(1d0)), intent(in) :: KetVecExp
    !
    real(kind(1d0)), external :: DGamma
    integer :: N, i
    real(kind(1d0)) :: Alpha
    !
    Integral = 0.d0
    !
    if ( BraMonExp(2) /= KetMonExp(2) ) return ! l
    if ( BraMonExp(3) /= KetMonExp(3) ) return ! m
    !
    N = BraMonExp(1) + KetMonExp(1) + 2
    Alpha = BraVecExp + KetVecExp
    !
    Integral = DGamma((dble(N) + 1.d0)*0.5d0)/&
         (2.d0*(Alpha**((dble(N) + 1.d0)*0.5d0)))
    !***
    if ( Integral > 1d200 ) then
       write(*,*) "Integral", Integral
       write(*,*) "N, Alpha, Gamma, Down", N, Alpha, DGamma((dble(N) + 1.d0)*0.5d0), 2.d0*(Alpha**((dble(N) + 1.d0)*0.5d0))
    end if
    !***
    !
  end function IntegralSpherSpher



  integer function DoubleFact(n) result( DF )
    !
    integer, intent(in) :: n
    !
    integer :: i
    !
    DF = 1
    !
    if ( n <= 0 ) return
    !
    i = n
    do
       DF = DF * i
       i = n -2 
       if ( i<=1 ) exit
    end do
    !
  end function DoubleFact



  subroutine ClassShortRangeSymOrbitalsFetchklmMask( SRSOrb, NewXlmToCartMat, L, M, NumExponents )
    !
    class(ClassShortRangeSymOrbitals), intent(in)  :: SRSOrb
    complex(kind(1d0)), allocatable,   intent(out) :: NewXlmToCartMat(:,:)
    integer,                           intent(in)  :: L
    integer,                           intent(in)  :: M
    integer,                           intent(in)  :: NumExponents
    !
    !> As appear in the interface file.
    integer, allocatable :: CartMonInterface(:,:)
    !
    !> As built in the initialization, might coincide with CartMonInterface, but not in general.
    integer, allocatable :: CartMonInitialized(:,:)
    complex(kind(1d0)), allocatable :: OldXlmToCartMat(:,:)
    complex(kind(1d0)), allocatable :: IntermediateXlmToCartMat(:,:)
    !
    integer :: OldNumRows, OldNumColumns, NewNumRows, NewNumColumns, i, j, ix, iy, iz, ChanL, ChanM, Counter
    integer, allocatable :: XlmTransfIndexes(:,:), AuxVec(:)
    !
    call SRSOrb%FetchCartMonInterface( CartMonInterface )
    call SRSOrb%FetchCartMonInitialized( CartMonInitialized )
    !
    call SRSOrb%FetchFullTransfMat( OldXlmToCartMat )
    !
    OldNumRows    = size(OldXlmToCartMat,1)
    OldNumColumns = size(OldXlmToCartMat,2)
    !
    allocate( IntermediateXlmToCartMat(OldNumRows,OldNumColumns) )
    IntermediateXlmToCartMat = Z0
    !
    !.. Looks which monomials of the interface file are present in the initialized set, 
    !   and arranges the order of the transformation matrix rows if it were neccesary 
    !   in order to maintain the compatibility with the matrix of coefficients that 
    !   expresses the diffuses orbitals in terms of the cartesians, and is defined 
    !   also in the interface file.
    !.. 
    do i = 1, size(CartMonInterface,2)
       do j = 1, size(CartMonInitialized,2)
          !
          if( sum( abs(CartMonInitialized(:,j)-CartMonInterface(:,i)) )/=0 )cycle
          !
          IntermediateXlmToCartMat( (i-1) * NumExponents + 1 : i * NumExponents, : ) = &
               OldXlmToCartMat( (j-1) * NumExponents + 1 : j * NumExponents, : )
          !
          exit
          !
       end do
    end do
    !
    !.. The 1st, 2nd and 3rd rows in XlmTransfIndexes array corresponds to n, l and m in r^nXlm% 
    !   Each columns is a diferent function.
    !..
    call SRSOrb%FetchXlmTransfIndexes( XlmTransfIndexes )
    !
    allocate( AuxVec(size(XlmTransfIndexes,2)) )
    AuxVec = 0
    !
    Counter = 0
    do i = 1, size(AuxVec)
       if ( (L == XlmTransfIndexes(2,i)) .and. &
            (M == XlmTransfIndexes(3,i)) ) then
          !
          Counter = Counter + 1
          AuxVec(Counter) = i
          !
       end if
    end do
    !
    !
    NewNumRows    = size(CartMonInterface,2) * NumExponents
    NewNumColumns = Counter * NumExponents
    !
    if ( Counter == 0 ) then
       NewNumColumns = NumExponents
       allocate( NewXlmToCartMat(NewNumRows,NewNumColumns) )
       NewXlmToCartMat = Z0
       return
    end if
    !
    allocate( NewXlmToCartMat(NewNumRows,NewNumColumns) )
    NewXlmToCartMat = Z0
    !
    do i = 1, Counter
       NewXlmToCartMat(:,(i-1)*NumExponents+1:i*NumExponents) = &
               IntermediateXlmToCartMat(:NewNumRows,(AuxVec(i)-1)*&
               NumExponents+1:AuxVec(i)*NumExponents)
    end do
    !
    deallocate( IntermediateXlmToCartMat )
    deallocate( OldXlmToCartMat )
    !
  end subroutine ClassShortRangeSymOrbitalsFetchklmMask



  !---------------------------------------------
  ! ClassDiffuseBsplineXlmBlock methods
  !---------------------------------------------


  subroutine ClassDiffuseBsplineXlmBlockInit( this, irrepBra, OpXlm, OpLabel, KetXlm, Block )
    class(ClassDiffuseBsplineXlmBlock), intent(inout) :: this
    class(ClassIrrep), target, intent(in) :: irrepBra
    class(ClassXlm)  , target, intent(in) :: OpXlm
    character(len=*)         , intent(in) :: OpLabel
    class(ClassXlm)  , target, intent(in) :: KetXlm
    class(ClassMatrix), optional, intent(in) :: Block
    
    this%DiffuseIrrep => irrepBra
    this%OperatorXlm  => OpXlm
    this%KetXlm       => KetXlm
    if(allocated(this%OperatorLabel))deallocate(this%OperatorLabel)
    allocate(this%OperatorLabel,source=OpLabel)
    if(present(Block)) this%Block = Block
    
  end subroutine ClassDiffuseBsplineXlmBlockInit



  subroutine ClassDiffuseBsplineXlmBlockFree( this )
    class(ClassDiffuseBsplineXlmBlock), intent(inout) :: this
    this%DiffuseIrrep => NULL()
    this%OperatorXlm  => NULL()
    this%KetXlm       => NULL()
    if(allocated(this%OperatorLabel))deallocate(this%OperatorLabel)
    call this%Block%Free()
  end subroutine ClassDiffuseBsplineXlmBlockFree




  function ClassDiffuseBsplineXlmBlockGetFile( this ) result( FileName )
    class(ClassDiffuseBsplineXlmBlock), intent(in) :: this
    character(len=:), allocatable :: FileName
    allocate(FileName, source=&
         this%OperatorLabel//&
         AlphabeticNumber(this%OperatorXlm%Getl())//"."//AlphabeticNumber(this%OperatorXlm%Getm())//"_"//&
         this%DiffuseIrrep%GetName()//"_"//&
         AlphabeticNumber(this%KetXlm%Getl())//"."//AlphabeticNumber(this%KetXlm%Getm()) )
  end function ClassDiffuseBsplineXlmBlockGetFile


  subroutine ClassDiffuseBsplineXlmBlockSave( this, Storage, WithFormat )
    class(ClassDiffuseBsplineXlmBlock), intent(inout) :: this
    character(len=*)                  , intent(in)    :: Storage
    logical, optional                 , intent(in)    :: WithFormat
    
    character(len=:), allocatable :: FileName
    integer                       :: uid
    logical :: Formatted
    
    if ( present(WithFormat) .and. WithFormat ) then
       Formatted = .true.
    else
       Formatted = .false.
    end if

    allocate(FileName, source= AddSlash(Storage)//AddSlash(DiffuseOrbDir)//this%GetFile())
    if ( Formatted ) then
       call OpenFile( FileName, uid, "write", "formatted" )
    else
       call OpenFile( FileName, uid, "write", "unformatted" )
    end if
    call this%Block%write(uid)
    close(uid)
    
  end subroutine ClassDiffuseBsplineXlmBlockSave



  subroutine ClassDiffuseBsplineXlmBlockReadBlock( this, Storage, Matrix )
    class(ClassDiffuseBsplineXlmBlock), intent(in)  :: this
    character(len=*)                  , intent(in)  :: Storage
    class(ClassMatrix)                , intent(out) :: Matrix
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid, iostat
    character(len=32) :: form
    character(len=IOMSG_LENGTH) :: iomsg
    logical :: exist
    !
    allocate(FileName, source= AddSlash(Storage)//AddSlash(DiffuseOrbDir)//this%GetFile())
    INQUIRE( file = FileName, exist = exist )
    if ( .not.exist ) then
       call ErrorMessage( 'The file '//FileName//' does not exist.' )
       return
    end if
    !
    call OpenFile( FileName, uid, "read", "formatted" )
    call Matrix%Read(uid)
    close(uid)
    !
  end subroutine ClassDiffuseBsplineXlmBlockReadBlock



  !---------------------------------------------
  ! ClassConditionerBlock methods
  !---------------------------------------------


  subroutine ClassConditionerBlockInit( this, MethodLabel, Threshold, &
                  NumBsDropBeginning, NumBsDropEnding, Irrep, KetXlm, &
                  PerfectProjector, Block )
    class(ClassConditionerBlock), intent(inout) :: this
    character(len=*)            , intent(in)    :: MethodLabel
    real(kind(1d0))             , intent(in)    :: Threshold
    integer                     , intent(in)    :: NumBsDropBeginning
    integer                     , intent(in)    :: NumBsDropEnding
    class(ClassIrrep), target   , intent(in)    :: Irrep
    class(ClassXlm)  , target   , intent(in)    :: KetXlm
    logical                     , intent(in)    :: PerfectProjector
    class(ClassMatrix), optional, intent(in)    :: Block
    !
    if(allocated(this%MethodLabel))deallocate(this%MethodLabel)
    allocate(this%MethodLabel,source = MethodLabel)
    this%Threshold          =  Threshold
    this%NumBsDropBeginning = NumBsDropBeginning
    this%NumBsDropEnding    = NumBsDropEnding
    if ( associated(this%Irrep) ) deallocate(this%Irrep)
    allocate( this%Irrep, source = Irrep )
    if ( associated(this%KetXlm) ) deallocate(this%KetXlm)
    allocate( this%KetXlm, source = KetXlm )
    if(present(Block)) this%Block = Block
    this%PerfectProjector = PerfectProjector
    !
  end subroutine ClassConditionerBlockInit



  subroutine ClassConditionerBlockFree( this )
    class(ClassConditionerBlock), intent(inout) :: this
    this%KetXlm    => NULL()
    this%Threshold =  -1.d0
    this%NumBsDropBeginning = -1
    this%NumBsDropEnding    = -1
    if(allocated(this%MethodLabel))deallocate(this%MethodLabel)
    call this%Block%Free()
    call this%AuxBlock%Free()
    if ( associated(this%Irrep) ) deallocate(this%Irrep)
    if ( associated(this%KetXlm) ) deallocate(this%KetXlm)
    this%PerfectProjector = .false.
  end subroutine ClassConditionerBlockFree



  subroutine ClassConditionerBlockSetXlm( this, Xlm )
    class(ClassConditionerBlock), intent(inout) :: this
    class(ClassXlm)             , intent(in)    :: Xlm
    if( associated(this%KetXlm) ) deallocate( this%KetXlm )
    allocate( this%KetXlm, source = Xlm )
  end subroutine ClassConditionerBlockSetXlm


  subroutine ClassConditionerBlockSetIrrep( this, Irrep )
    class(ClassConditionerBlock), intent(inout) :: this
    class(ClassIrrep)           , intent(in)    :: Irrep
    if( associated(this%Irrep) ) deallocate( this%Irrep )
    allocate( this%Irrep, source = Irrep )
  end subroutine ClassConditionerBlockSetIrrep



  function ClassConditionerBlockGetFile( this ) result( FileName )
    class(ClassConditionerBlock), intent(in) :: this
    character(len=:), allocatable :: FileName
    character(len=64) :: Strn
    character(len=:), allocatable :: FileNameRoot, ConditionerLabel
    !
    allocate( ConditionerLabel, source = this%GetLabel() )
    !
    allocate( FileName, source = &
         ConditionerLabel//'_'//&
         this%Irrep%GetName()//'_'//&
         AlphabeticNumber(this%KetXlm%Getl())//"."//&
         AlphabeticNumber(this%KetXlm%Getm()) )
    !
  end function ClassConditionerBlockGetFile


  function ClassConditionerBlockGetDiffuseFile( this ) result( FileName )
    class(ClassConditionerBlock), intent(in) :: this
    character(len=:), allocatable :: FileName
    !
    allocate( FileName, source = this%GetFile()//'_'//DiffuseLabel )
    !
  end function ClassConditionerBlockGetDiffuseFile



  function ClassConditionerBlockGetLabel( this ) result(Label)
    class(ClassConditionerBlock), intent(in) :: this
    character(len=:), allocatable :: Label
    !
    character(len=64) :: Strn
    character(len=:), allocatable :: RootLabel
    !
    if ( this%MethodLabel .is. MethodRM ) then
       allocate( RootLabel, source = &
            'BsInit'//AlphabeticNumber(this%NumBsDropBeginning)//&
            'BsEnd'//AlphabeticNumber(this%NumBsDropEnding)//&
            '_'//this%MethodLabel )
    elseif ( this%MethodLabel .is. MethodLI ) then
       write(Strn,*) this%Threshold
       allocate( RootLabel, source = trim(adjustl(Strn))//'_'//this%MethodLabel )
    elseif ( this%MethodLabel .is. MethodLIO ) then
       write(Strn,*) this%Threshold
       allocate( RootLabel, source = trim(adjustl(Strn))//'_'//this%MethodLabel )
    else
       call Assert( "Invalid conditioning method" )
    end if
    !
    if ( this%PerfectProjector ) then
       if ( (this%MethodLabel .is. MethodLI) .or. &
            (this%MethodLabel .is. MethodLIO) ) then
          allocate( Label, source = RootLabel//'_'//PerfectProjectorLabel )
          return
       end if
    end if
    !
    allocate( Label, source = RootLabel )
    !
  end function ClassConditionerBlockGetLabel



  subroutine ClassConditionerBlockSave( this, Storage )
    class(ClassConditionerBlock), intent(inout) :: this
    character(len=*)            , intent(in)    :: Storage
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid
    !
    call execute_command_line("mkdir -p "//AddSlash(Storage)//AddSlash(ConditionerDir) )
    !
    allocate(FileName, source= AddSlash(Storage)//AddSlash(ConditionerDir)//this%GetFile())
    call OpenFile( FileName, uid, "write", "formatted" )
    call this%Block%write(uid)
    close(uid)
    deallocate( FileName )
    !
    if ( this%AuxBlock%IsInitialized() ) then
       allocate(FileName, source= AddSlash(Storage)//AddSlash(ConditionerDir)//this%GetDiffuseFile())
       call OpenFile( FileName, uid, "write", "formatted" )
       call this%AuxBlock%write(uid)
       close(uid)
       deallocate( FileName )
    end if
    !
  end subroutine ClassConditionerBlockSave



  subroutine ClassConditionerBlockReadBlock( this, Storage, Matrix )
    class(ClassConditionerBlock), intent(in)  :: this
    character(len=*)            , intent(in)  :: Storage
    class(ClassMatrix)          , intent(out) :: Matrix
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid
    logical :: exist
    !
    allocate(FileName, source= AddSlash(Storage)//AddSlash(ConditionerDir)//this%GetFile())
    INQUIRE( file = FileName, exist = exist )
    if ( .not.exist ) then
       call ErrorMessage( 'The conditioner stored in '//FileName//' does not exist.' )
       return
    end if
    call OpenFile( FileName, uid, "read", "formatted" )
    call Matrix%Read(uid)
    close(uid)
    !
  end subroutine ClassConditionerBlockReadBlock


  subroutine ClassConditionerBlockReadDiffuseBlock( this, Storage, Matrix )
    class(ClassConditionerBlock), intent(in)  :: this
    character(len=*)            , intent(in)  :: Storage
    class(ClassMatrix)          , intent(out) :: Matrix
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid
    logical :: exist
    !
    allocate(FileName, source= AddSlash(Storage)//AddSlash(ConditionerDir)//this%GetDiffuseFile())
    INQUIRE( file = FileName, exist = exist )
    if ( .not.exist ) then
       call ErrorMessage( 'The conditioner stored in '//FileName//' does not exist.' )
       return
    end if
    call OpenFile( FileName, uid, "read", "formatted" )
    call Matrix%Read(uid)
    close(uid)
    !
  end subroutine ClassConditionerBlockReadDiffuseBlock



  subroutine ClassConditionerBlockComputeConditioner( this, &
       Storage, SRSOrbSet                                 , &
       DBSMatrix, BSBSMatrix                        )
    !
    class(ClassConditionerBlock)     , intent(inout) :: this
    character(len=*)                 , intent(in)    :: Storage
    Class(ClassShortRangeSymOrbitals), intent(in)    :: SRSOrbSet
    class(ClassMatrix)               , intent(in)    :: DBSMatrix
    class(ClassMatrix)               , intent(in)    :: BSBSMatrix
    !
    if ( this%MethodLabel .is. MethodRM ) then
       call this%ComputeRMConditioner( DBSMatrix%NColumns(), this%NumBsDropBeginning, this%NumBsDropEnding )
!!$    elseif ( this%MethodLabel .is. MethodLI ) then
!!$       call this%ComputeLIConditioner( DBSMatrix, BSBSMatrix, Storage, SRSOrbSet )
!!$    elseif ( this%MethodLabel .is. MethodLIO ) then
!!$       call this%ComputeLIOConditioner( DBSMatrix, BSBSMatrix, Storage, SRSOrbSet )
    else
       call Assert( 'The conditioning method: '//this%MethodLabel//' is not valid.' )
    end if
    !
  end subroutine ClassConditionerBlockComputeConditioner




  subroutine ClassConditionerBlockComputeRMConditioner( this, NRows, NumBsDropBeginning, NumBsDropEnding )
    !
    class(ClassConditionerBlock), intent(inout) :: this
    integer                     , intent(in)    :: NRows
    integer                     , intent(in)    :: NumBsDropBeginning
    integer                     , intent(in)    :: NumBsDropEnding
    !
    integer :: NCols, i, j
    !
    NCols = NRows - NumBsDropBeginning - NumBsDropEnding
    !
    call this%Block%InitFull( NRows, NCols )
    !
    do j = 1, NCols
       call this%Block%SetElement( NumBsDropBeginning+j, j, 1.d0 )
    end do
    !
  end subroutine ClassConditionerBlockComputeRMConditioner



!!$  subroutine ClassConditionerBlockComputeLIConditioner( this, DBSMatrix, BSBSMatrix, Storage, SRSOrbSet )
!!$    !
!!$    class(ClassConditionerBlock)     , intent(inout) :: this
!!$    class(ClassMatrix)               , intent(in)    :: DBSMatrix
!!$    class(ClassMatrix)               , intent(in)    :: BSBSMatrix
!!$    character(len=*)                 , intent(in)    :: Storage
!!$    Class(ClassShortRangeSymOrbitals), intent(in)    :: SRSOrbSet
!!$    !
!!$    real(kind(1d0)), parameter :: LOWER_PROJECTOR_TOLERANCE = -1.d-5
!!$    real(kind(1d0)), parameter :: UPPER_PROJECTOR_TOLERANCE =  1.d-5
!!$    real(kind(1d0)), parameter :: OverlapThreshold = 1.d-32
!!$    !
!!$    integer :: NRows, NCols, i, j
!!$    integer :: NumDiff, NumBs, NumOverBs, NPrecon
!!$    real(kind(1d0)), allocatable :: CopyArray(:,:), FinalProjector(:,:), Eval(:)
!!$    real(kind(1d0)) :: Overlap
!!$    type(Classmatrix) :: Projector, Projector2, GaussGauss, GaussGaussInv, CorrectBsBs
!!$    type(ClassSpectralResolution) :: SpecRes
!!$    !
!!$    NumDiff = DBSMatrix%NRows()
!!$    NumBs   = DBSMatrix%NColumns()
!!$    !
!!$    call DBSMatrix%FetchMatrix( CopyArray )
!!$    NumOverBs = NumBs
!!$    OverlapCycle : do 
!!$       Overlap = dot_product(CopyArray(:,NumOverBs),CopyArray(:,NumOverBs))
!!$       if( abs(Overlap) > OverlapThreshold )exit OverlapCycle
!!$       NumOverBs = NumOverBs - 1
!!$       if(NumOverBs == 0) exit
!!$    enddo OverlapCycle
!!$    !
!!$    if( NumOverBs == 0 )then
!!$       call ErrorMessage("All BSplines are orthogonal to the Gaussian basis set.")
!!$       call this%Block%InitFull( NumBs, NumBs ) 
!!$       do i = 1, NumOverBs
!!$          call this%Block%SetElement( i, i, 1.d0 )
!!$       end do
!!$       return
!!$    endif
!!$    write(OUTPUT_UNIT,"(a,i5)") " Number of overlapping BSplines : ", NumOverBs
!!$    !
!!$    Projector = CopyArray(:,1:NumOverBs)
!!$    call Projector.Transpose( Projector2 )
!!$    !
!!$    if ( this%PerfectProjector ) then
!!$       call SRSOrbSet%LoadDiffMonOverlap( Storage, GaussGauss )
!!$       write(*,*) "GG is Symmetric", GaussGauss.IsSymmetric(1.d-10)
!!$       call GaussGauss.Inverse( GaussGaussInv )
!!$       write(*,*) "GGinv is Symmetric", GaussGaussInv.IsSymmetric(1.d-10)
!!$       call GaussGauss.Free()
!!$       call Projector.Multiply( GaussGaussInv, "Left", "N" )
!!$       call GaussGaussInv.Free()
!!$       call Projector.Multiply( Projector2, "Left", "N" )
!!$    else
!!$       call Projector.Multiply( Projector2, "Left", "N" )
!!$    end if
!!$    !
!!$    call Projector2.Free()
!!$    call BSBSMatrix%GetSubMatrix( &
!!$         1, NumOverBs, &
!!$         1, NumOverBs, &
!!$         CorrectBsBs )
!!$    call Projector.Diagonalize( CorrectBsBs, SpecRes, EPSILON(1d0) )
!!$    call CorrectBsBs.Free()
!!$    call SpecRes.Fetch( FinalProjector )
!!$    call SpecRes.Fetch( Eval )
!!$    !
!!$    !
!!$    NPrecon = 0
!!$    do i = 1, NumOverBs
!!$       !
!!$       if( Eval(i) < LOWER_PROJECTOR_TOLERANCE )then
!!$          write(output_unit,*) Eval(i)
!!$          call Assert("Eigenvalue below lower projector tolerance")
!!$       endif
!!$       !
!!$       if( (Eval(i) < 1.d0 - this%Threshold) ) then
!!$          NPrecon = NPrecon + 1
!!$          FinalProjector(:,NPrecon) = FinalProjector(:,i)
!!$       else
!!$          write(OUTPUT_UNIT,*) "Removed eigenvalue: ", Eval(i)
!!$       endif
!!$       !
!!$    end do
!!$    !
!!$    !
!!$    write(OUTPUT_UNIT,"(a,i5)") " Number of B-splines' combinations removed : ", NumOverBs-NPrecon
!!$    !
!!$    call this%Block%InitFull( NumBs, NumBs-(NumOverBs-NPrecon) ) 
!!$    do j = 1, NPrecon
!!$       do i = 1, NumOverBs
!!$          call this%Block%SetElement( i, j, FinalProjector(i,j) )
!!$       end do
!!$    end do
!!$    !
!!$    do i = NumOverBs+1, NumBs 
!!$       call this%Block%SetElement( i, i-(NumOverBs-NPrecon), 1.d0 )
!!$    end do
!!$    !
!!$  end subroutine ClassConditionerBlockComputeLIConditioner



!!$  subroutine ClassConditionerBlockComputeLIOConditioner( this, DBSMatrix, BSBSMatrix, Storage, SRSOrbSet )
!!$    !
!!$    class(ClassConditionerBlock)     , intent(inout) :: this
!!$    class(ClassMatrix)               , intent(in)    :: DBSMatrix
!!$    class(ClassMatrix)               , intent(in)    :: BSBSMatrix
!!$    character(len=*)                 , intent(in)    :: Storage
!!$    Class(ClassShortRangeSymOrbitals), intent(in)    :: SRSOrbSet
!!$    !
!!$    real(kind(1d0)), parameter :: LOWER_PROJECTOR_TOLERANCE = -1.d-5
!!$    real(kind(1d0)), parameter :: UPPER_PROJECTOR_TOLERANCE =  1.d-5
!!$    real(kind(1d0)), parameter :: OverlapThreshold = 1.d-32
!!$    !
!!$    integer :: NRows, NCols, i, j
!!$    integer :: NumDiff, NumBs, NumOverBs, NPrecon
!!$    real(kind(1d0)), allocatable :: CopyArray(:,:), FinalProjector(:,:), Eval(:)
!!$    real(kind(1d0)) :: Overlap
!!$    type(Classmatrix) :: Projector, Projector2, GaussGauss, GaussGaussInv, CorrectBsBs, AuxMat, BsExpanMat
!!$    type(ClassSpectralResolution) :: SpecRes
!!$    !
!!$    NumDiff = DBSMatrix%NRows()
!!$    NumBs   = DBSMatrix%NColumns()
!!$    !
!!$    call DBSMatrix%FetchMatrix( CopyArray )
!!$    NumOverBs = NumBs
!!$    OverlapCycle : do 
!!$       Overlap = dot_product(CopyArray(:,NumOverBs),CopyArray(:,NumOverBs))
!!$       if( abs(Overlap) > OverlapThreshold )exit OverlapCycle
!!$       NumOverBs = NumOverBs - 1
!!$       if(NumOverBs == 0) exit
!!$    enddo OverlapCycle
!!$    !
!!$    if( NumOverBs == 0 )then
!!$       call ErrorMessage("All BSplines are orthogonal to the Gaussian basis set.")
!!$       call this%Block%InitFull( NumBs, NumBs ) 
!!$       do i = 1, NumOverBs
!!$          call this%Block%SetElement( i, i, 1.d0 )
!!$       end do
!!$       call this%AuxBlock%InitFull( NumDiff, NumBs ) 
!!$       return
!!$    endif
!!$    write(OUTPUT_UNIT,"(a,i5)") " Number of overlapping BSplines : ", NumOverBs
!!$    !
!!$    Projector = CopyArray(:,1:NumOverBs)
!!$    call Projector.Transpose( Projector2 )
!!$    !
!!$    if ( this%PerfectProjector ) then
!!$       call SRSOrbSet%LoadDiffMonOverlap( Storage, GaussGauss )
!!$       write(*,*) "GG is Symmetric", GaussGauss.IsSymmetric(1.d-10)
!!$       call GaussGauss.Inverse( GaussGaussInv )
!!$       write(*,*) "GGinv is Symmetric", GaussGaussInv.IsSymmetric(1.d-10)
!!$       call AuxMat%InitFull( GaussGauss.NRows(),GaussGauss.NColumns() )
!!$       do i = 1, GaussGauss.NRows()
!!$          call AuxMat%SetElement( i, i, 1.d0 )
!!$       end do
!!$       call AuxMat%Multiply( 2.d0 )
!!$       call GaussGauss.Multiply( -1.d0 )
!!$       call AuxMat%Add( GaussGauss )
!!$       call GaussGauss.Free()
!!$       call Projector.Multiply( AuxMat, "Left", "N" )
!!$    end if
!!$    !
!!$    call Projector.Multiply( Projector2, "Left", "N" )
!!$    call Projector2.Free()
!!$    call BSBSMatrix%GetSubMatrix( &
!!$         1, NumOverBs, &
!!$         1, NumOverBs, &
!!$         CorrectBsBs )
!!$    !
!!$    call Projector.Multiply( -1.d0 )
!!$    call Projector.Add( CorrectBsBs )
!!$    call CorrectBsBs.Free()
!!$    call Projector.ConditionMatrix( this%Threshold, BsExpanMat )
!!$    !
!!$    NPrecon = BsExpanMat%NColumns()
!!$    call BsExpanMat%FetchMatrix( FinalProjector )
!!$    call BsExpanMat%Free()
!!$    !
!!$    write(OUTPUT_UNIT,"(a,i5)") " Number of B-splines' combinations removed : ", NumOverBs-NPrecon
!!$    !
!!$    call this%Block%InitFull( NumBs, NumBs-(NumOverBs-NPrecon) ) 
!!$    do j = 1, NPrecon
!!$       do i = 1, NumOverBs
!!$          call this%Block%SetElement( i, j, FinalProjector(i,j) )
!!$       end do
!!$    end do
!!$    !
!!$    do i = NumOverBs+1, NumBs 
!!$       call this%Block%SetElement( i, i-(NumOverBs-NPrecon), 1.d0 )
!!$    end do
!!$    !
!!$    this%AuxBlock = DBSMatrix
!!$    call this%AuxBlock%Multiply( -1.d0 )
!!$    call this%AuxBlock%Multiply( this%Block, 'Right', 'N' )
!!$    if ( this%PerfectProjector ) then
!!$       call this%AuxBlock%Multiply( GaussGaussInv, 'Left', 'N' )
!!$    end if
!!$    !
!!$  end subroutine ClassConditionerBlockComputeLIOConditioner



  !---------------------------------------------
  ! General purpose methods
  !---------------------------------------------

  function GetDiffuseOrbDir() result(Dir)
    character(len=:), allocatable :: Dir
    allocate( Dir, source = DiffuseOrbDir )
  end function GetDiffuseOrbDir


  function GetConditionerDir() result(Dir)
    character(len=:), allocatable :: Dir
    allocate( Dir, source = ConditionerDir )
  end function GetConditionerDir


  subroutine GetConditioningMethodList( MethodList )
    character(len=*), allocatable, intent(inout) :: MethodList(:)
    integer :: NumMethods
    NumMethods = 3
    allocate( MethodList(NumMethods) )
    MethodList(1) = MethodRM
    MethodList(2) = MethodLI
    MethodList(3) = MethodLIO
  end subroutine GetConditioningMethodList


  subroutine CheckMethod( Method )
    character(len=*), intent(in) :: Method
    character(len=64), allocatable :: MethodList(:)
    integer :: i, NumMethods
    logical :: MethodFound
    !
    MethodFound = .false.
    call GetConditioningMethodList( MethodList )
    NumMethods = size(MethodList)
    !
    do i = 1, NumMethods
       if ( Method .is. trim(adjustl(MethodList(i))) ) then
          MethodFound = .true.
          exit
       end if
    end do
    !
    if ( .not.MethodFound ) then
       write( output_unit,*) 'Invalid conditioning method. It must be either: '   , &
            (trim(adjustl(MethodList(i)))//', ',i=1,NumMethods-1), &
            trim(adjustl(MethodList(NumMethods)))//'.'
       stop
    end if
    !
  end subroutine CheckMethod


  subroutine CheckThreshold( Thr )
    real(kind(1d0)), intent(in) :: Thr
    if ( Thr < 0 ) call Assert( 'The conditioning threshold cannot be negative.' )
  end subroutine CheckThreshold


  subroutine CheckBsToRemove( NumBs )
    integer, intent(in) :: NumBs
    if ( NumBs < 0 ) call Assert( 'The number of B-splines to remove cannot be negative.' )
  end subroutine CheckBsToRemove


  subroutine CheckLastBsIsPresent( Method, NumBsDropEnding )
    character(len=*), intent(in) :: Method
    integer         , intent(in) :: NumBsDropEnding
    if ( (Method .is. MethodRM) .and. &
         (NumBsDropEnding /= 0) ) call Assert( &
         'The last B-spline is not present.' )
  end subroutine CheckLastBsIsPresent


  !> Gets the number of B-splines that do not overlap with the diffuse orbitals.
  integer function GetNumBsNotOverlappingDiffuse( DBSMatrix, BSBSMatrix ) result(NBs)
    !
    class(ClassMatrix)               , intent(in)    :: DBSMatrix
    class(ClassMatrix)               , intent(in)    :: BSBSMatrix
    !
    real(kind(1d0)), parameter :: OverlapThreshold = 1.d-32
    !
    integer :: NumDiff, NumBs, NumOverBs
    real(kind(1d0)), allocatable :: CopyArray(:,:)
    real(kind(1d0)) :: Overlap
    !
    NumDiff = DBSMatrix%NRows()
    NumBs   = DBSMatrix%NColumns()
    !
    call DBSMatrix%FetchMatrix( CopyArray )
    NumOverBs = NumBs
    OverlapCycle : do 
       Overlap = dot_product(CopyArray(:,NumOverBs),CopyArray(:,NumOverBs))
       if( abs(Overlap) > OverlapThreshold )exit OverlapCycle
       NumOverBs = NumOverBs - 1
       if(NumOverBs == 0) exit
    enddo OverlapCycle
    !
    NBs = NumBs - NumOverBs
    if ( NBs <= 0 ) call Assert( &
         "There aren't B-splines not overlapping with diffuse Gaussian." )
    !
  end function GetNumBsNotOverlappingDiffuse



  subroutine SaveQCConditioner( StorageDir, Irrep, QCConditionerMat )
    character(len=*)  , intent(in) :: StorageDir
    class(ClassIrrep) , intent(in) :: Irrep
    class(ClassMatrix), intent(in) :: QCConditionerMat
    !
    character(len=:), allocatable :: FileName
    integer :: uid
    !
    allocate( FileName, source = BuildQCConditionerFileName(StorageDir,Irrep) )
    call OpenFile( FileName, uid, "write", "formatted" )
    call QCConditionerMat%Write(uid)
    close(uid)
    !
  end subroutine SaveQCConditioner

  

  subroutine ReadQCConditioner( StorageDir, Irrep, QCConditionerMat )
    character(len=*)  , intent(in)  :: StorageDir
    class(ClassIrrep) , intent(in)  :: Irrep
    class(ClassMatrix), intent(out) :: QCConditionerMat
    !
    character(len=:), allocatable :: FileName
    integer :: uid
    !
    allocate( FileName, source = BuildQCConditionerFileName(StorageDir,Irrep) )
    call OpenFile( FileName, uid, "read", "formatted" )
    call QCConditionerMat%Read(uid)
    close(uid)
    !
  end subroutine ReadQCConditioner



  function BuildQCConditionerFileName( StorageDir, Irrep ) result(FileName)
    character(len=*)  , intent(in) :: StorageDir
    class(ClassIrrep) , intent(in) :: Irrep
    character(len=:) , allocatable :: FileName
    allocate( FileName, source = &
         AddSlash(StorageDir)//AddSlash(ConditionerDir)//&
         QCConditionerLabel//'_'//Irrep%GetName() )
  end function BuildQCConditionerFileName




end module ModuleShortRangeOrbitals
