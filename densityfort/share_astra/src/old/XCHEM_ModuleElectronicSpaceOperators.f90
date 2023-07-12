module ModuleElectronicSpaceOperators


  use, intrinsic :: ISO_FORTRAN_ENV


  use ModuleErrorHandling
  use ModuleMatrix
  use ModuleElementaryElectronicSpaceOperators
  use ModuleSymmetricElectronicSpaceOperators
  use ModuleSymmetricElectronicSpace
  use ModuleGroups
  use ModuleElectronicSpace




  !> Class of the operators electronic space blocks: S, H, D, etc.
  type, public :: ClassESBlock
     private
     type(ClassElectronicSpace), pointer :: Space
     !> Blocks correspondind to all the combinations of
     !! symmetric electronic space bras and kets, defined
     !! for the point group and the close-coupling expansion.
     type(ClassSESSESBlock), allocatable :: Block(:,:)
     logical                                :: BoxStatesOnly  = .FALSE.
     logical                                :: FormattedWrite = .FALSE.
   contains
     generic, public :: Init                => ClassESBlockInit
     generic, public :: Load                => ClassESBlockLoad
     generic, public :: Save                => ClassESBlockSave
     generic, public :: Assemble            => ClassESBlockAssemble
     generic, public :: SetBoxOnly          => ClassESBlockSetBoxOnly
     generic, public :: UnsetBoxOnly        => ClassESBlockUnsetBoxOnly
     generic, public :: SetFormattedWrite   => ClassESBlockSetFormattedWrite
     generic, public :: UnsetFormattedWrite => ClassESBlockUnsetFormattedWrite
     generic, public :: Free                => ClassESBlockFree
     ! {{{ private procedures

     procedure, private :: ClassESBlockInit
     procedure, private :: ClassESBlockLoad
     procedure, private :: ClassESBlockSave
     procedure, private :: ClassESBlockAssemble
     procedure, private :: ClassESBlockSetBoxOnly
     procedure, private :: ClassESBlockUnsetBoxOnly
     procedure, private :: ClassESBlockSetFormattedWrite
     procedure, private :: ClassESBlockUnsetFormattedWrite
     procedure, private :: ClassESBlockFree
     final              :: ClassESBlockFinal

     ! }}}
  end type ClassESBlock


  public :: GetSymmetricElectSpaces


contains


  !.. From the label of the bra space and the ket space
  !   returns a pointer to a symmetric space for both the
  !   bra and the ket, and tells whether they are compatible
  !   with the operator (it is a bit contrieved ...)
  !..
  subroutine GetSymmetricElectSpaces( Space, &
       BraSymLabel  , &
       KetSymLabel  , &
       OperatorLabel, &
       DipoleAxis   , &
       BraSymSpace  , &
       KetSymSpace  )
    class(ClassElectronicSpace), target , intent(in)    :: Space
    character(len=:), allocatable       , intent(inout) :: BraSymLabel
    character(len=:), allocatable       , intent(inout) :: KetSymLabel
    character(len=*)                    , intent(in)    :: OperatorLabel
    character(len=:), allocatable       , intent(inout) :: DipoleAxis
    !
    type(ClassSymmetricElectronicSpace), pointer, intent(out) :: BraSymSpace
    type(ClassSymmetricElectronicSpace), pointer, intent(out) :: KetSymSpace
    !
    type(ClassGroup), pointer :: Group
    type(ClassIrrep), pointer :: BraIrrep, KetIrrep, AuxIrrep
    type(ClassIrrep)          :: OpIrrep
    character(len=:), allocatable :: BraIrrepName, KetIrrepName
    !
    Group => Space.GetGroup()
    !
    if ( allocated(BraSymLabel) ) then
       call Space.CheckSymmetry( BraSymLabel )
       BraIrrep => Group.GetIrrep( BraSymLabel )
    end if
    !
    if ( allocated(KetSymLabel) ) then
       call Space.CheckSymmetry( KetSymLabel )
       KetIrrep => Group.GetIrrep( KetSymLabel )
    end if
    !
    OpIrrep = GetOperatorIrrep( Group, OperatorLabel, DipoleAxis )
    if ( allocated(BraSymLabel) .and. allocated(KetSymLabel) ) then
       AuxIrrep => OpIrrep * KetIrrep
       if ( .not.BraIrrep.NameIs( AuxIrrep.GetName() ) ) call Assert( &
            'The bra, ket and operator irreducible representations do not match.' )
    end if
    !
    if ( allocated(BraSymLabel) ) then
       KetIrrep => OpIrrep * BraIrrep
    end if
    if ( allocated(KetSymLabel) ) then
       BraIrrep => OpIrrep * KetIrrep
    end if
    !
    allocate( BraIrrepName, source = BraIrrep.GetName() )
    allocate( KetIrrepName, source = KetIrrep.GetName() )
    !
    call Space.CheckSymmetry( BraIrrepName )
    call Space.CheckSymmetry( KetIrrepName )
    !
    allocate( BraSymSpace, source = Space.GetSymElectSpace(BraIrrepName) )
    allocate( KetSymSpace, source = Space.GetSymElectSpace(KetIrrepName) )
    !
  end subroutine GetSymmetricElectSpaces


  !--------------------------------------------------
  ! Methods for ClassESBlock
  !-------------------------------------------------

  subroutine ClassESBlockFinal( Self )
    type(ClassESBlock) :: Self
    call Self.Free()
  end subroutine ClassESBlockFinal


  subroutine ClassESBlockFree( Self )
    class(ClassESBlock), intent(inout) :: Self
    Self.Space => NULL()
    if ( allocated(Self.Block) ) deallocate( Self.Block )
    Self.BoxStatesOnly  = .false.
    Self.FormattedWrite = .false.
  end subroutine ClassESBlockFree


  !> Initialize the electronic space operator matrix.
  subroutine ClassESBlockInit( Self, Space, OpLabel, Axis, Force )
    !
    !> Class of the electronic space operator matrix.
    class(ClassESBlock)                , intent(inout) :: Self
    !> Class of the electronic space.
    class(ClassElectronicSpace), target, intent(in)    :: Space
    !> Operator label.
    character(len=*)                   , intent(in)    :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable      , intent(inout) :: Axis
    logical                            , intent(in)    :: Force
    !
    integer :: i, j
    integer :: NIrreps
    type(ClassSymmetricElectronicSpace), pointer :: BraSymSpace, KetSymSpace
    !
    if ( .not.Space.Initialized() ) call Assert( "To initialize the electronic space operator block the space must be initialized." )
    !
    Self.Space => Space
    NIrreps = Space.GetNIrreps()
    if ( allocated(Self.Block) ) deallocate( Self.Block )
    allocate( Self.Block(NIrreps,NIrreps) )
    !
    do j = 1, NIrreps
       KetSymSpace => Space.GetSymElectSpace(j)
       do i = 1, NIrreps
          if ( Self.BoxStatesOnly ) call Self.Block(i,j).SetBoxOnly()
          if ( Self.FormattedWrite ) call Self.Block(i,j).SetFormattedWrite()
          BraSymSpace => Space.GetSymElectSpace(i)
          if ( .not.ValidSymmetries( OpLabel, BraSymSpace, KetSymSpace, Axis ) ) cycle
          call Self.Block(i,j).Init( &
               BraSymSpace         , &
               KetSymSpace         , &
               OpLabel             , &
               Axis                , &
               Force               )
       end do
    end do
    !
  end subroutine ClassESBlockInit



  !> Load the electronic space operator matrix.
  subroutine ClassESBlockLoad( Self, Space, OpLabel, Axis )
    !
    !> Class of the electronic space operator matrix.
    class(ClassESBlock)                , intent(inout) :: Self
    !> Class of the electronic space.
    class(ClassElectronicSpace), target, intent(in)    :: Space
    !> Operator label.
    character(len=*)                   , intent(in)    :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable      , intent(inout) :: Axis
    !
    integer :: i, j
    integer :: NIrreps
    type(ClassSymmetricElectronicSpace), pointer :: BraSymSpace, KetSymSpace
    !
    if ( .not.Space.Initialized() ) call Assert( "To load the electronic space operator block the space must be initialized." )
    !
    Self.Space => Space
    !
    NIrreps = Space.GetNIrreps()
    if ( allocated(Self.Block) ) deallocate( Self.Block )
    allocate( Self.Block(NIrreps,NIrreps) )
    !
    do j = 1, NIrreps
       KetSymSpace => Space.GetSymElectSpace(j)
       do i = 1, NIrreps
          if ( Self.BoxStatesOnly ) call Self.Block(i,j).SetBoxOnly()
          if ( Self.FormattedWrite ) call Self.Block(i,j).SetFormattedWrite()
          BraSymSpace => Space.GetSymElectSpace(i)
          if ( .not.ValidSymmetries( OpLabel, BraSymSpace, KetSymSpace, Axis ) ) cycle
          call Self.Block(i,j).Load( &
               BraSymSpace         , &
               KetSymSpace         , &
               OpLabel             , &
               Axis                )
       end do
    end do
    !
  end subroutine ClassESBlockLoad


  !> Save the electronic space operator matrix.
  subroutine ClassESBlockSave( Self )
    !
    !> Class of the electronic space operator matrix.
    class(ClassESBlock)        , intent(inout) :: Self
    !
    integer :: i, j
    integer :: NIrreps
    !
    if ( .not.Self.Space.Initialized() ) call Assert( "To load the electronic space operator block the space must be initialized." )
    !
    NIrreps = Self.Space.GetNIrreps()
    if ( .not.allocated(Self.Block) ) call Assert( 'The electronic space block must be allocated to save it.' )
    !
    do j = 1, NIrreps
       do i = 1, NIrreps
          call Self.Block(i,j).Save( )
       end do
    end do
    !
  end subroutine ClassESBlockSave



  !> Load the electronic space operator matrix.
  subroutine ClassESBlockAssemble( Self, LoadLocStates, OpMat )
    !
    !> Class of the electronic space operator matrix.
    class(ClassESBlock), intent(inout) :: Self
    logical            , intent(in)    :: LoadLocStates
    !> Full space class matrix representation of the operator.
    class(ClassMatrix) , intent(out)   :: OpMat
    !
    integer :: i, j
    integer, allocatable :: NumRows(:), NumCols(:)
    integer :: NIrreps
    type(ClassMatrix), allocatable :: SymOpMat(:,:), LoadedSymOpMat(:,:)
    !
    if ( .not.Self.Space.Initialized() ) call Assert( "To load the electronic space operator block the space must be initialized." )
    !
    NIrreps = Self.Space.GetNIrreps()
    if ( .not.allocated(Self.Block) ) call Assert( &
         'The electronic space block must be initialized before assemble.' )
    !
    allocate( SymOpMat(NIrreps,NIrreps) )
    allocate( LoadedSymOpMat(NIrreps,NIrreps) )
    allocate( NumRows(NIrreps), NumCols(NIrreps) )
    !
    NumRows = 0
    NumCols = 0
    do j = 1, NIrreps
       do i = 1, NIrreps
          if ( Self.BoxStatesOnly ) call Self.Block(i,j).SetBoxOnly()
          call Self.Block(i,j).Assemble( LoadLocStates, SymOpMat(i,j) )
          if ( Self.Block(i,j).BlockIsAvailable() ) then
             NumRows(i) = SymOpMat(i,j).NRows()
             NumCols(j) = SymOpMat(i,j).NColumns()
          end if
       end do
    end do

    do j = 1, NIrreps
       do i = 1, NIrreps
          if ( Self.Block(i,j).BlockIsAvailable() ) then
             LoadedSymOpMat(i,j) = SymOpMat(i,j)
             call SymOpMat(i,j).Free()
          else
             call LoadedSymOpMat(i,j).InitFull(NumRows(i),NumCols(j))
          end if
       end do
    end do

    !
    deallocate( SymOpMat )
    call OpMat.BuildUpMatrix( LoadedSymOpMat )
    deallocate( LoadedSymOpMat )
    !
  end subroutine ClassESBlockAssemble



  !> When loaded, the last B-spline in the operator matrix elements is removed.
  subroutine ClassESBlockSetBoxOnly( Self )
    !> Class of the electronic space operator.
    class(ClassESBlock), intent(inout) :: Self
    Self.BoxStatesOnly = .true.
  end subroutine ClassESBlockSetBoxOnly


  !> When loaded, the last B-spline in the operator matrix elements is included.
  subroutine ClassESBlockUnsetBoxOnly( Self )
    !> Class of the electronic space operator.
    class(ClassESBlock), intent(inout) :: Self
    Self.BoxStatesOnly = .false.
  end subroutine ClassESBlockUnsetBoxOnly


  !> The blocks will be write and read with format.
  subroutine ClassESBlockSetFormattedWrite( Self )
    !> Class of the electronic space operator.
    class(ClassESBlock), intent(inout) :: Self
    Self.FormattedWrite = .true.
  end subroutine ClassESBlockSetFormattedWrite



  !> The blocks will be write and read with format.
  subroutine ClassESBlockUnsetFormattedWrite( Self )
    !> Class of the electronic space operator.
    class(ClassESBlock), intent(inout) :: Self
    Self.FormattedWrite = .false.
  end subroutine ClassESBlockUnsetFormattedWrite




end module ModuleElectronicSpaceOperators
