module ModuleSymmetricElectronicSpaceOperators


  use, intrinsic :: ISO_FORTRAN_ENV


  use ModuleErrorHandling
  use ModuleMatrix
  use ModuleSymmetryAdaptedSphericalHarmonics
  use ModuleSymmetricLocalizedStates
  use ModuleShortRangeChannel
  use ModulePartialWaveChannel
  use ModuleSymmetricElectronicSpace
  use ModuleElectronicSpace
  use ModuleShortRangeOrbitals
  use ModuleGroups
  use ModuleElementaryElectronicSpaceOperators




  character(len=*), parameter :: RowsLabel = 'ROWS'
  character(len=*), parameter :: ColsLabel = 'COLUMNS'
  logical         , parameter :: DoBsPermutation = .true.
  
  
  !> Class of the symmetric electronic space blocks, devoted to manage the matrices: S, H, D, ... defined in this space.
  type, public :: ClassSESSESBlock
     ! {{{ private attributes
     private
     !
     !> Space of the Bra functions.
     type(ClassSymmetricElectronicSpace), pointer :: SpaceBra
     !
     !> Space of the Ket functions.
     type(ClassSymmetricElectronicSpace), pointer :: SpaceKet

     !.. The Q space comprises only states obtained by augmenting ions with electrons
     !   in active molecular orbitals (i.e., present in the CI expansion of the parent ions).
     !   In this program, this Q space corresponds to the so-called localized component of
     !   a channel.
     !
     !   The P spaces has three components:
     !
     !   P0: ions times native FEDVR that do not overalp with the ions at all
     !       and hence which retain well defined lm and the FEDVR radial properties.
     !
     !   P1: ions times linear combinations of FEDVR with lm well defined,
     !       which do preside the same radial region as the MOs, but which are 
     !       nonetheless orthogonal to the MOs
     !
     !   P0 and P1, which have well-defined Xlm angular dependence, formed the so-called 
     !   partial-wave component of the channel.
     !
     !   P2: This space contains hybrid FEDVR-MO functions. All the elements are
     !       at the moment treated on equal footing, but under the hood the space
     !       has in fact three distinct components 
     !       !
     !       P2a) ions times linear combinations of pure FEDVR withOUT lm well 
     !            defined, in the same region as the MOs, and orthogonal to them
     !       !
     !       P2b) ions times linear combinatins of FEDVR and MOs (real hybrid basis)
     !       !
     !       P2c) ions times virtual MOs, i.e., orbitals that are expressed exclusively 
     !            in terms of Gaussian functions but which are not represented in the 
     !            expansion of the parent ions
     !
     !   The P2a-c spaces are indicated here as the short-range component of the channel
     !..

     !> Q - Q block
     type(LocLocBlock),allocatable                :: LocLocBlock(:,:)

     !> Q - P2 block
     type(LocSRCBlock),allocatable                :: LocSRCBlock(:,:)

     !> Q - P0/1
     type(LocPWCBlock),allocatable                :: LocPWCBlock(:,:)
     
     !> P2 - Q
     type(SRCLocBlock),allocatable                :: SRCLocBlock(:,:)

     !> P2 - P2
     type(SRCSRCBlock),allocatable                :: SRCSRCBlock(:,:)

     !> P2 - P0/1
     type(SRCPWCBlock),allocatable                :: SRCPWCBlock(:,:)

     !> P0/1 - Q
     type(PWCLocBlock),allocatable                :: PWCLocBlock(:,:)

     !> P0/1 - P2
     type(PWCSRCBlock),allocatable                :: PWCSRCBlock(:,:)

     !> P0/1 - P0/1
     type(PWCPWCBlock),allocatable                :: PWCPWCBlock(:,:)

     logical                                      :: BoxStatesOnly      = .FALSE.
     logical                                      :: BraBoxStatesOnly   = .FALSE.
     logical                                      :: KetBoxStatesOnly   = .FALSE.
     logical                                      :: AsympBsPermutation = .FALSE.
     logical                                      :: FormattedWrite     = .FALSE.
     logical                                      :: AvailableBlocks    = .FALSE.
     ! }}}
     !
   contains
     !
     generic, public :: Init                    => ClassSESSESBlockInit
     generic, public :: InitAndSave             => ClassSESSESBlockInitAndSave
     generic, public :: Load                    => ClassSESSESBlockLoad
     generic, public :: LoadQC                  => ClassSESSESBlockLoadQC
     generic, public :: LoadLoc                 => ClassSESSESBlockLoadLoc
     generic, public :: LoadPoly                => ClassSESSESBlockLoadPoly
     generic, public :: Save                    => ClassSESSESBlockSave
     generic, public :: Assemble                => ClassSESSESBlockAssemble
     generic, public :: AssembleQC              => ClassSESSESBlockAssembleQC
     generic, public :: AssembleLoc             => ClassSESSESBlockAssembleLoc
     generic, public :: AssemblePoly             => ClassSESSESBlockAssemblePoly
     generic, public :: Condition               => ClassSESSESBlockCondition
     generic, public :: SetBoxOnly              => ClassSESSESBlockSetBoxOnly
     generic, public :: SetBraBoxOnly           => ClassSESSESBlockSetBraBoxOnly
     generic, public :: SetKetBoxOnly           => ClassSESSESBlockSetKetBoxOnly
     generic, public :: UnsetBoxOnly            => ClassSESSESBlockUnsetBoxOnly
     generic, public :: UnsetBraBoxOnly         => ClassSESSESBlockUnsetBraBoxOnly
     generic, public :: UnsetKetBoxOnly         => ClassSESSESBlockUnsetKetBoxOnly
     generic, public :: SetAsympBsPermutation   => ClassSESSESBlockSetAsympBsPermutation
     generic, public :: UnsetAsympBsPermutation => ClassSESSESBlockUnsetAsympBsPermutation
     generic, public :: SetFormattedWrite       => ClassSESSESBlockSetFormattedWrite
     generic, public :: UnsetFormattedWrite     => ClassSESSESBlockUnsetFormattedWrite
     generic, public :: IsBoxOnly               => ClassSESSESBlockIsBoxOnly
     generic, public :: IsBraBoxOnly            => ClassSESSESBlockIsBraBoxOnly
     generic, public :: IsKetBoxOnly            => ClassSESSESBlockIsKetBoxOnly
     generic, public :: BlockIsAvailable        => ClassSESSESBlockIsAvailable
     generic, public :: GetStorageDir           => ClassSESSESBlockGetStorageDir
     generic, public :: GetNumFunQC             => ClassSESSESBlockGetNumFunQC
     generic, public :: GetMinNumNotOvBs        => ClassSESSESBlockGetMinNumNotOvBs
     generic, public :: Free                    => ClassSESSESBlockFree
     ! {{{ private procedures

     procedure, private :: ClassSESSESBlockInit
     procedure, private :: ClassSESSESBlockInitAndSave
     procedure, private :: ClassSESSESBlockAssemble
     procedure, private :: ClassSESSESBlockAssembleQC
     procedure, private :: ClassSESSESBlockAssembleLoc
     procedure, private :: ClassSESSESBlockAssemblePoly
     procedure, private :: ClassSESSESBlockCondition
     procedure, private :: ClassSESSESBlockLoad
     procedure, private :: ClassSESSESBlockLoadQC
     procedure, private :: ClassSESSESBlockLoadLoc
     procedure, private :: ClassSESSESBlockLoadPoly
     procedure, private :: ClassSESSESBlockSave
     procedure, private :: ClassSESSESBlockSetBoxOnly
     procedure, private :: ClassSESSESBlockSetBraBoxOnly
     procedure, private :: ClassSESSESBlockSetKetBoxOnly
     procedure, private :: ClassSESSESBlockUnsetBoxOnly
     procedure, private :: ClassSESSESBlockUnsetBraBoxOnly
     procedure, private :: ClassSESSESBlockUnsetKetBoxOnly
     procedure, private :: ClassSESSESBlockSetAsympBsPermutation
     procedure, private :: ClassSESSESBlockUnsetAsympBsPermutation
     procedure, private :: ClassSESSESBlockSetFormattedWrite
     procedure, private :: ClassSESSESBlockUnsetFormattedWrite
     procedure, private :: ClassSESSESBlockIsBoxOnly
     procedure, private :: ClassSESSESBlockIsBraBoxOnly
     procedure, private :: ClassSESSESBlockIsKetBoxOnly
     procedure, private :: ClassSESSESBlockIsAvailable
     procedure, private :: ClassSESSESBlockGetStorageDir
     procedure, private :: ClassSESSESBlockGetNumFunQC
     procedure, private :: ClassSESSESBlockGetMinNumNotOvBs
     procedure, private :: ClassSESSESBlockFree
     final              :: ClassSESSESBlockFinal

     ! }}}
  end type ClassSESSESBlock

  

  public :: ValidSymmetries
  public :: FormatOperatorLabel
  public :: FormatGaugeLabel
  public :: FormatAxisLabel
  public :: GetBlockStorageDir
  public :: SetLabel




contains



  function GetBlockStorageDir( SpaceBra, SpaceKet ) result(Dir)
    class(ClassSymmetricElectronicSpace), intent(in) :: SpaceBra
    class(ClassSymmetricElectronicSpace), intent(in) :: SpaceKet
    character(len=:), allocatable :: Dir
    !
    character(len=:), allocatable :: NucConfDir, CCDir, BraSymLabel, KetSymLabel
    !
    allocate( NucConfDir , source = SpaceBra.GetStorageDir() )
    allocate( CCDir      , source = GetCloseCouplingDir()    )
    allocate( BraSymLabel, source = SpaceBra.GetLabel()      )
    allocate( KetSymLabel, source = SpaceKet.GetLabel()      )
    !
    allocate( Dir, source = &
         AddSlash(NucConfDir)//&
         AddSlash(CCDir)//&
         BraSymLabel//'_'//AddSlash(KetSymLabel) )
    !
  end function GetBlockStorageDir




  subroutine FormatOperatorLabel( OperatorLabel )
    character(len=:), allocatable, intent(inout) :: OperatorLabel
    if ( OperatorLabel .is. OverlapLabel ) then
       call SetLabel( OverlapLabel, OperatorLabel )
    elseif ( OperatorLabel .is. HamiltonianLabel ) then
       call SetLabel( HamiltonianLabel, OperatorLabel )
    elseif ( OperatorLabel .is. KinEnergyLabel ) then
       call SetLabel( KinEnergyLabel, OperatorLabel )
    elseif ( OperatorLabel .is. DIPOLE_LABEL ) then
       call SetLabel( DIPOLE_LABEL, OperatorLabel )
    elseif ( OperatorLabel .is. CAPLabel ) then
       call SetLabel( CAPLabel, OperatorLabel )
    else
       call Assert( 'Invalid operator label, it must be one of this: '//&
            OverlapLabel//', '//&
            HamiltonianLabel//', '//&
            KinEnergyLabel//', '//&
            DIPOLE_LABEL//' or '//&
            CAPLabel//'.' )
    end if
  end subroutine FormatOperatorLabel


  subroutine FormatGaugeLabel( Gauge )
    character(len=:), allocatable, intent(inout) :: Gauge
    if ( Gauge .is. 'L' ) then
       call SetLabel( LengthLabel, Gauge )
    elseif ( Gauge .is. 'V' ) then
       call SetLabel( VelocityLabel, Gauge )
    else
       call Assert( 'Invalid gauge label, it must be "l" or "v".' )
    end if
  end subroutine FormatGaugeLabel



  subroutine FormatAxisLabel( DipoleAxis )
    character(len=*), intent(in) :: DipoleAxis
    if ( (DipoleAxis .isnt. 'X') .and. &
         (DipoleAxis .isnt. 'Y') .and. &
         (DipoleAxis .isnt. 'Z') ) &
         call Assert( 'Invalid dipole orientation, it must be either "x", "y" or "z".' )
  end subroutine FormatAxisLabel



  subroutine SetLabel( NewLabel, OldLabel )
    character(len=*),              intent(in)    :: NewLabel
    character(len=:), allocatable, intent(inout) :: OldLabel
    if ( allocated(OldLabel) ) deallocate( OldLabel )
    allocate( OldLabel, source = NewLabel )
  end subroutine SetLabel



  logical function ValidSymmetries( OpLabel, BraSymSpace, KetSymSpace, Axis ) result(ValidSym)
    character(len=*)                            , intent(in) :: OpLabel
    class(ClassSymmetricElectronicSpace), target, intent(in) :: BraSymSpace
    class(ClassSymmetricElectronicSpace), target, intent(in) :: KetSymSpace    
    character(len=:) , allocatable              , intent(inout) :: Axis
    type(ClassIrrep), pointer :: BraIrrep, KetIrrep
    type(ClassIrrep) :: OpIrrep
    type(ClassGroup), pointer :: Group
    !
    ValidSym = .false.
    !
    BraIrrep => BraSymSpace.GetIrrep()
    KetIrrep => KetSymSpace.GetIrrep()
    Group    => BraIrrep.GetGroup()
    OpIrrep  =  GetOperatorIrrep( Group, OpLabel, Axis )
    !
    if ( ValidIrreps(BraIrrep,OpIrrep,KetIrrep) ) ValidSym = .true.
    !
  end function ValidSymmetries


  !> Rearrange the B-splines in such a way that those that do not
  !! overlap with the diffuse Gaussian are reorder in increasing
  !! number of B-splines' index, for all the channels. Let's suppose
  !! that from a set of 100 B-splines tha last 50 do not overlap with
  !! the diffuse orbitals, and that we have 3 PWC. If we identify the
  !! B-splines as \f$Bs_{i}^{j}\f$, where \f$i\f$ is the function
  !! index and \f$j\f$ corresponds to the channel, originally, the
  !! B-splines are organized as follows:
  !!
  !! \f$Bs_{1}^{1}, Bs_{2}^{1}, ..., Bs_{100}^{1}, 
  !! Bs_{1}^{2}, Bs_{2}^{2}, ..., Bs_{100}^{2},
  !! Bs_{1}^{3}, Bs_{2}^{3}, ..., Bs_{100}^{3}, ...\f$.
  !!
  !! After transformation the order is the following:
  !!
  !! \f$Bs_{1}^{1}, Bs_{2}^{1}, ..., Bs_{50}^{1},  
  !! Bs_{1}^{2}, Bs_{2}^{2}, ..., Bs_{50}^{2},
  !! Bs_{1}^{3}, Bs_{2}^{3}, ..., Bs_{50}^{3},
  !! Bs_{51}^{1}, Bs_{51}^{2}, ..., Bs_{51}^{3}, 
  !! Bs_{52}^{1}, Bs_{52}^{2}, ..., Bs_{52}^{3}, ...
  !! Bs_{100}^{1}, Bs_{100}^{2}, ..., Bs_{100}^{3}\f$ 
  !!
  subroutine PermuteBsplines( FullMat, SubSpacesMat, SubSpacesBlock, Which )
    !> ClassMatrix build from the SubSpacesMat array.
    class(ClassMatrix), intent(inout) :: FullMat
    !> Array of ClassMatrix containing all the blocks for a definite
    !! kind of sub-spaces, i.e.: all the < loc | O | pwc Xlm >, or
    !! < pwc Xlm | O | pwc Xlm >, etc.
    class(ClassMatrix), intent(in)    :: SubSpacesMat(:,:)
    !> Array of a general sub-spaces block i.e.: 
    !! LocPWCBlock, PWCPWCBlock, etc.
    class(*), target  , intent(in)    :: SubSpacesBlock(:,:)
    !> Spacifies whether the rows or the columns will be swapped.
    character(len=*)  , intent(in)    :: Which
    !
    integer :: MinNumNotOverlapingBs
    integer :: NumOverlapingBs, NTotOvBs
    integer :: i, N, j
    integer :: OldIndexBeg, NewIndexBeg, OldIndexEnd, NewIndexEnd
    type(ClassMatrix) :: CopyFullMat
    real(kind(1d0)), allocatable :: NewArray(:,:), OldArray(:,:)
    !
    if ( (Which .isnt. RowsLabel) .and. &
         (Which .isnt. ColsLabel) ) call Assert( &
         "Invalid parameter for permutation. "//&
         "It must be "//RowsLabel//" or "//ColsLabel//". " )
    !
    if ( (size(SubSpacesMat,1) /= size(SubSpacesBlock,1)) .or. &
         (size(SubSpacesMat,2) /= size(SubSpacesBlock,2)) ) call Assert( &
         'The dimensions of the matrix array and the block array must be the same.' )
    !
    CopyFullMat = FullMat
    call FullMat.Free()
    call CopyFullMat.FetchMatrix(OldArray)
    allocate( NewArray(CopyFullMat.NRows(),CopyFullMat.NColumns()) )
    NewArray = 0.d0
    !
    MinNumNotOverlapingBs = GetMinNumNotOverlapingBs( SubSpacesBlock, Which )
    !
    !
    if ( Which .is. RowsLabel ) then
       N = size(SubSpacesMat,1)
    else
       N = size(SubSpacesMat,2)
    end if
    
    
    !.. Permutation of the overlapping B-splines
    OldIndexBeg = 1
    NewIndexBeg = 1
    !
    NTotOvBs = 0
    do i = 1, N
       !
       if ( Which .is. RowsLabel ) then
          NumOverlapingBs = SubSpacesMat(i,1).NRows() - MinNumNotOverlapingBs
       else
          NumOverlapingBs = SubSpacesMat(1,i).NColumns() - MinNumNotOverlapingBs
       end if
       !
       NTotOvBs = NTotOvBs + NumOverlapingBs
       OldIndexEnd = OldIndexBeg + NumOverlapingBs - 1
       NewIndexEnd = NewIndexBeg + NumOverlapingBs - 1
       !
       if ( Which .is. RowsLabel ) then
          NewArray(NewIndexBeg:NewIndexEnd,:) = OldArray(OldIndexBeg:OldIndexEnd,:)
          OldIndexBeg = OldIndexBeg + SubSpacesMat(i,1).NRows()
       else
          NewArray(:,NewIndexBeg:NewIndexEnd) = OldArray(:,OldIndexBeg:OldIndexEnd)
          OldIndexBeg = OldIndexBeg + SubSpacesMat(1,i).NColumns()
       end if
       !
       NewIndexBeg = NewIndexEnd + 1
       !
    end do
    
    
    !.. 1st step in the permutation of not overlapping B-splines
    if ( Which .is. RowsLabel ) then
       OldIndexBeg = SubSpacesMat(1,1).NRows() - MinNumNotOverlapingBs + 1
    else
       OldIndexBeg = SubSpacesMat(1,1).NColumns() - MinNumNotOverlapingBs + 1
    end if
    !
    NewIndexBeg = NTotOvBs + 1
    !
    do i = 1, N
       !
       OldIndexEnd = OldIndexBeg + MinNumNotOverlapingBs - 1
       NewIndexEnd = NewIndexBeg + MinNumNotOverlapingBs - 1
       !
       if ( Which .is. RowsLabel ) then
          NewArray(NewIndexBeg:NewIndexEnd,:) = OldArray(OldIndexBeg:OldIndexEnd,:)
          if ( i < N ) then
             OldIndexBeg = OldIndexEnd + SubSpacesMat(i+1,1).NRows() - MinNumNotOverlapingBs + 1
          else
             exit
          end if
       else
          NewArray(:,NewIndexBeg:NewIndexEnd) = OldArray(:,OldIndexBeg:OldIndexEnd)
          if ( i < N ) then
             OldIndexBeg = OldIndexEnd + SubSpacesMat(1,i+1).NColumns() - MinNumNotOverlapingBs + 1
          else
             exit
          end if
       end if
       !
       NewIndexBeg = NewIndexEnd + 1
       !
    end do
    
    
    !.. 2nd step in the permutation of not overlapping B-splines
    deallocate( OldArray )
    allocate( OldArray, source = NewArray )
    !
    OldIndexBeg = NTotOvBs + 1
    NewIndexBeg = NTotOvBs + 1
    !
    do i = 1, MinNumNotOverlapingBs
       do j = 1, N
          if ( Which .is. RowsLabel ) then
             NewArray(NewIndexBeg + j-1,:) = OldArray(OldIndexBeg + (j-1)*MinNumNotOverlapingBs,:)
          else
             NewArray(:,NewIndexBeg + j-1) = OldArray(:,OldIndexBeg + (j-1)*MinNumNotOverlapingBs)
          end if
       end do
       OldIndexBeg = OldIndexBeg + 1
       NewIndexBeg = NewIndexBeg + N
    end do
    
    FullMat = NewArray
    deallocate( NewArray )
    !
  end subroutine PermuteBsplines
  
  

  !> Gets the minimum number of not overlapping B-splines among all channels.
  integer function GetMinNumNotOverlapingBs( SubSpacesBlock, Which ) result(MinNBs)
    !> Array of a general sub-spaces block i.e.: 
    !! LocPWCBlock, PWCPWCBlock, etc.
    class(*), target, intent(in) :: SubSpacesBlock(:,:)
    !> Spacifies whether the rows or the columns will be swapped.
    character(len=*), intent(in) :: Which
    !
    type(LocPWCBlock),pointer :: LocPWC
    type(SRCPWCBlock),pointer :: SRCPWC
    type(PWCLocBlock),pointer :: PWCLoc
    type(PWCSRCBlock),pointer :: PWCSRC
    type(PWCPWCBlock),pointer :: PWCPWC
    integer :: i, N, ierror
    integer, allocatable :: Nvec(:)
    !
    if ( Which .is. RowsLabel ) then
       !
       N = size(SubSpacesBlock,1)
       allocate( Nvec(N) )
       !
       do i = 1, N
          select type( ptr => SubSpacesBlock(i,1) )
          class is ( PWCLocBlock )
             allocate( PWCLoc, source = SubSpacesBlock(i,1) )
             Nvec(i) = PWCLoc.GetNumBsNotOvDiffBra()
             NULLIFY( PWCLoc )
          class is ( PWCSRCBlock )
             allocate( PWCSRC, source = SubSpacesBlock(i,1) )
             Nvec(i) = PWCSRC.GetNumBsNotOvDiffBra()
             NULLIFY( PWCSRC )
          class is ( PWCPWCBlock )
             allocate( PWCPWC, source = SubSpacesBlock(i,1) )
             Nvec(i) = PWCPWC.GetNumBsNotOvDiffBra()
             NULLIFY( PWCPWC )
          class default 
             call Assert( "Invalid electronic subspace block." )
          end select
       end do
       !
    elseif ( Which .is. ColsLabel ) then
       !
       N = size(SubSpacesBlock,2)
       allocate( Nvec(N) )
       !
       do i = 1, N
          select type( ptr => SubSpacesBlock(1,i) )
          class is ( LocPWCBlock )
             allocate( LocPWC, source = SubSpacesBlock(1,i) )
             Nvec(i) = LocPWC.GetNumBsNotOvDiffKet()
             NULLIFY( LocPWC )
          class is ( SRCPWCBlock )
             allocate( SRCPWC, source = SubSpacesBlock(1,i) )
             Nvec(i) = SRCPWC.GetNumBsNotOvDiffKet()
             NULLIFY( SRCPWC )
          class is ( PWCPWCBlock )
             allocate( PWCPWC, source = SubSpacesBlock(1,i), STAT=ierror )
             if ( ierror /= 0) call Assert( 'Error associating PWCPWC pointer.' )
             Nvec(i) = PWCPWC.GetNumBsNotOvDiffKet()
             NULLIFY( PWCPWC )
          class default 
             call Assert( "Invalid electronic subspace block." )
          end select
       end do
       !
    end if
    !
    MinNBs = MINVAL( Nvec )
    if ( MinNBs <= 0 ) call Assert( &
         'The minimum number of B-splines not overlapping with '//&
         'diffuse Gaussian must be higher than 0.' )
    !
  end function GetMinNumNotOverlapingBs



  !--------------------------------------------------
  ! Methods for ClassSESSESBlock
  !-------------------------------------------------


  subroutine ClassSESSESBlockFinal( Self )
    type(ClassSESSESBlock) :: Self
    call Self.Free()
  end subroutine ClassSESSESBlockFinal


  subroutine ClassSESSESBlockFree( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    integer :: i, j
    Self.SpaceBra => NULL()
    Self.SpaceKet => NULL()
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    if ( allocated(Self.LocSRCBlock) ) deallocate( Self.LocSRCBlock )
    if ( allocated(Self.LocPWCBlock) ) deallocate( Self.LocPWCBlock )
    if ( allocated(Self.SRCLocBlock) ) deallocate( Self.SRCLocBlock )
    if ( allocated(Self.SRCSRCBlock) ) deallocate( Self.SRCSRCBlock )
    if ( allocated(Self.SRCPWCBlock) ) deallocate( Self.SRCPWCBlock )
    if ( allocated(Self.PWCLocBlock) ) deallocate( Self.PWCLocBlock )
    if ( allocated(Self.PWCSRCBlock) ) deallocate( Self.PWCSRCBlock )
    if ( allocated(Self.PWCPWCBlock) ) deallocate( Self.PWCPWCBlock )
    Self.BoxStatesOnly      = .false.
    Self.BraBoxStatesOnly   = .false.
    Self.KetBoxStatesOnly   = .false.
    Self.AsympBsPermutation = .false.
    Self.FormattedWrite     = .false.
    Self.AvailableBlocks    = .false.
  end subroutine ClassSESSESBlockFree



  logical function ClassSESSESBlockIsBoxOnly( Self ) result( IsBoxOnly )
    class(ClassSESSESBlock), intent(in) :: Self
    IsBoxOnly = Self.BoxStatesOnly
  end function ClassSESSESBlockIsBoxOnly


  logical function ClassSESSESBlockIsBraBoxOnly( Self ) result( IsBoxOnly )
    class(ClassSESSESBlock), intent(in) :: Self
    IsBoxOnly = Self.BraBoxStatesOnly
  end function ClassSESSESBlockIsBraBoxOnly


  logical function ClassSESSESBlockIsKetBoxOnly( Self ) result( IsBoxOnly )
    class(ClassSESSESBlock), intent(in) :: Self
    IsBoxOnly = Self.KetBoxStatesOnly
  end function ClassSESSESBlockIsKetBoxOnly


  logical function ClassSESSESBlockIsAvailable( Self ) result( BlockAv )
    class(ClassSESSESBlock), intent(in) :: Self
    BlockAv = Self.AvailableBlocks
  end function ClassSESSESBlockIsAvailable



  subroutine ClassSESSESBlockUnsetFormattedWrite( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.FormattedWrite = .false.
  end subroutine ClassSESSESBlockUnsetFormattedWrite


  subroutine ClassSESSESBlockSetFormattedWrite( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.FormattedWrite = .true.
  end subroutine ClassSESSESBlockSetFormattedWrite



  subroutine ClassSESSESBlockUnsetBoxOnly( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.BoxStatesOnly = .false.
  end subroutine ClassSESSESBlockUnsetBoxOnly


  subroutine ClassSESSESBlockUnsetBraBoxOnly( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.BraBoxStatesOnly = .false.
  end subroutine ClassSESSESBlockUnsetBraBoxOnly


  subroutine ClassSESSESBlockUnsetKetBoxOnly( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.KetBoxStatesOnly = .false.
  end subroutine ClassSESSESBlockUnsetKetBoxOnly

  
  subroutine ClassSESSESBlockSetBoxOnly( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.BoxStatesOnly = .true.
  end subroutine ClassSESSESBlockSetBoxOnly


  subroutine ClassSESSESBlockSetBraBoxOnly( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.BraBoxStatesOnly = .true.
  end subroutine ClassSESSESBlockSetBraBoxOnly


  subroutine ClassSESSESBlockSetKetBoxOnly( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.KetBoxStatesOnly = .true.
  end subroutine ClassSESSESBlockSetKetBoxOnly


  subroutine ClassSESSESBlockSetAsympBsPermutation( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.AsympBsPermutation = .true.
  end subroutine ClassSESSESBlockSetAsympBsPermutation


  subroutine ClassSESSESBlockUnsetAsympBsPermutation( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self.AsympBsPermutation = .false.
  end subroutine ClassSESSESBlockUnsetAsympBsPermutation




  !> Initializes the ClassSESSESBlock class.
  subroutine ClassSESSESBlockInit( Self, SpaceBra, SpaceKet, OpLabel, Axis, Force )
    !> Class to be initialized.
    class(ClassSESSESBlock),                      intent(inout) :: Self
    !> Space of the Bra functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)     :: SpaceBra
    !> Space of the Ket functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)     :: SpaceKet
    !> Operator Label.
    character(len=*),                            intent(in)     :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable              , intent(inout)  :: Axis
    logical,                                     intent(in)     :: Force
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumPWCBra, NumPWCKet, NumLocBra, NumLocKet
    type(ClassPartialWaveChannel)       , pointer :: BraPWC, KetPWC
    type(ClassShortRangeChannel)        , pointer :: BraSRC, KetSRC
    type(ClassSymmetricLocalizedStates) , pointer :: BraLoc, KetLoc
    type(ClassXlm) :: OpXlm
    !
    Self.SpaceBra => SpaceBra
    Self.SpaceKet => SpaceKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    NumPWCBra = Self.SpaceBra.GetNumPWC()
    NumPWCKet = Self.SpaceKet.GetNumPWC()
    !
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    if ( allocated(Self.LocSRCBlock) ) deallocate( Self.LocSRCBlock )
    if ( allocated(Self.LocPWCBlock) ) deallocate( Self.LocPWCBlock )
    if ( allocated(Self.SRCLocBlock) ) deallocate( Self.SRCLocBlock )
    if ( allocated(Self.SRCSRCBlock) ) deallocate( Self.SRCSRCBlock )
    if ( allocated(Self.SRCPWCBlock) ) deallocate( Self.SRCPWCBlock )
    if ( allocated(Self.PWCLocBlock) ) deallocate( Self.PWCLocBlock )
    if ( allocated(Self.PWCSRCBlock) ) deallocate( Self.PWCSRCBlock )
    if ( allocated(Self.PWCPWCBlock) ) deallocate( Self.PWCPWCBlock )
    !
    allocate( Self.LocLocBlock(NumLocBra,NumLocKet) )
    allocate( Self.LocSRCBlock(NumLocBra,NumSRCKet) )
    allocate( Self.LocPWCBlock(NumLocBra,NumPWCKet) )
    allocate( Self.SRCLocBlock(NumSRCBra,NumLocKet) )
    allocate( Self.SRCSRCBlock(NumSRCBra,NumSRCKet) )
    allocate( Self.SRCPWCBlock(NumSRCBra,NumPWCKet) )
    allocate( Self.PWCLocBlock(NumPWCBra,NumLocKet) )
    allocate( Self.PWCSRCBlock(NumPWCBra,NumSRCKet) )
    allocate( Self.PWCPWCBlock(NumPWCBra,NumPWCKet) )
    !
    !
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    !
    !.. < * | O | loc >
    do j = 1, NumLocKet

       KetLoc  => Self.SpaceKet.GetChannelLS()

       !.. < loc | O | loc >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocLocBlock(i,j).SetFormattedWrite()
          call Self.LocLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocLocBlock(i,j).Build( BraLoc, KetLoc, OpLabel, OpXlm, Force )
       end do
       
       !.. Sum over the parent ions
       do i = 1, NumSRCBra

          !.. < src | O | loc >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCLocBlock(i,j).SetFormattedWrite()
          call Self.SRCLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCLocBlock(i,j).Build( BraSRC, KetLoc, OpLabel, OpXlm, Force )
          
          !.. < pwc | O | loc >
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             
             !.. < pwc Xlm | O | loc >
             BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
             CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
             if ( Self.BoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetKetBoxOnly()
             if ( Self.FormattedWrite ) call Self.PWCLocBlock(CounterBra,j).SetFormattedWrite()
             call Self.PWCLocBlock(CounterBra,j).SetFileExtension( FileExtensionPWCLoc  ) 
             call Self.PWCLocBlock(CounterBra,j).Build( BraPWC, KetLoc, OpLabel, OpXlm, Force )
             
          end do

       end do

    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       
       KetSRC => Self.SpaceKet.GetChannelSRC(j)

       !.. < loc | O | src >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocSRCBlock(i,j).SetFormattedWrite()
          call Self.LocSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocSRCBlock(i,j).Build( BraLoc, KetSRC, OpLabel, OpXlm, Force )
       end do

       do i = 1, NumSRCBra

          !.. < src | O | src >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCSRCBlock(i,j).SetFormattedWrite()
          call Self.SRCSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCSRCBlock(i,j).Build( BraSRC, KetSRC, OpLabel, OpXlm, Force )
          
          do n = 1, Self.SpaceBra.GetNumPWC(i)

             !.. < pwc Xlm | O | src >
             BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
             CounterBra =  Self.SpaceBra.PWCChannelIndex(i,n)
             if ( Self.BoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetKetBoxOnly()
             if ( Self.FormattedWrite ) call Self.PWCSRCBlock(CounterBra,j).SetFormattedWrite()
             call Self.PWCSRCBlock(CounterBra,j).SetFileExtension( FileExtensionPWCSRC  ) 
             call Self.PWCSRCBlock(CounterBra,j).Build( BraPWC, KetSRC, OpLabel, OpXlm, Force )
             
          end do
          
       end do

    end do


    !.. < * | O | pwc >
    do j = 1, NumSRCKet
       do k = 1, Self.SpaceKet.GetNumPWC(j)
          KetPWC => Self.SpaceKet.GetChannelPWC(j,k)
          CounterKet = Self.SpaceKet.PWCChannelIndex(j,k)

          !.. < loc | O | pwc Xlm >
          do i = 1, NumLocBra
             BraLoc  => Self.SpaceBra.GetChannelLS()
             if ( Self.BoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetKetBoxOnly()
             if ( Self.FormattedWrite ) call Self.LocPWCBlock(i,CounterKet).SetFormattedWrite()
             call Self.LocPWCBlock(i,CounterKet).SetFileExtension( FileExtensionLocPWC  ) 
             call Self.LocPWCBlock(i,CounterKet).Build( BraLoc, KetPWC, OpLabel, OpXlm, Force )
          end do
          
          !.. < src | O | pwc Xlm >
          do i = 1, NumSRCBra
             BraSRC => Self.SpaceBra.GetChannelSRC(i)
             if ( Self.BoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetKetBoxOnly()
             if ( Self.FormattedWrite ) call Self.SRCPWCBlock(i,CounterKet).SetFormattedWrite()
             call Self.SRCPWCBlock(i,CounterKet).SetFileExtension( FileExtensionSRCPWC  ) 
             call Self.SRCPWCBlock(i,CounterKet).Build( BraSRC, KetPWC, OpLabel, OpXlm, Force )
             
             do n = 1, Self.SpaceBra.GetNumPWC(i)

                !.. < pwc Xlm | O | pwc Xlm >
                BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                if ( Self.BoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetBoxOnly()
                if ( Self.BraBoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetBraBoxOnly()
                if ( Self.KetBoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetKetBoxOnly()
                if ( Self.FormattedWrite ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetFormattedWrite()
                call Self.PWCPWCBlock(CounterBra,CounterKet).SetFileExtension( FileExtensionPWCPWC  ) 
                call Self.PWCPWCBlock(CounterBra,CounterKet).Build( BraPWC, KetPWC, OpLabel, OpXlm, Force )
                
             end do
          end do

       end do
    end do


    !
  end subroutine ClassSESSESBlockInit




  !> Initializes the ClassSESSESBlock class.
  subroutine ClassSESSESBlockInitAndSave( Self, &
    SpaceBra, &
    SpaceKet, &
    OpLabel , &
    Axis    , &
    Force   , &
    ExtraLabel )
    !> Class to be initialized.
    class(ClassSESSESBlock),                      intent(inout) :: Self
    !> Space of the Bra functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)     :: SpaceBra
    !> Space of the Ket functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)     :: SpaceKet
    !> Operator Label.
    character(len=*),                            intent(in)     :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable              , intent(inout)  :: Axis
    logical,                                     intent(in)     :: Force
    character(len=*), optional                 , intent(in)    :: ExtraLabel
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumPWCBra, NumPWCKet, NumLocBra, NumLocKet
    type(ClassPartialWaveChannel)       , pointer :: BraPWC, KetPWC
    type(ClassShortRangeChannel)        , pointer :: BraSRC, KetSRC
    type(ClassSymmetricLocalizedStates) , pointer :: BraLoc, KetLoc
    type(ClassXlm) :: OpXlm
    !
    Self.SpaceBra => SpaceBra
    Self.SpaceKet => SpaceKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    NumPWCBra = Self.SpaceBra.GetNumPWC()
    NumPWCKet = Self.SpaceKet.GetNumPWC()
    !
    write(*,*) NumLocBra, NumSRCBra, NumPWCBra
    write(*,*) NumLocKet, NumSRCKet, NumPWCKet
    !
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    if ( allocated(Self.LocSRCBlock) ) deallocate( Self.LocSRCBlock )
    if ( allocated(Self.LocPWCBlock) ) deallocate( Self.LocPWCBlock )
    if ( allocated(Self.SRCLocBlock) ) deallocate( Self.SRCLocBlock )
    if ( allocated(Self.SRCSRCBlock) ) deallocate( Self.SRCSRCBlock )
    if ( allocated(Self.SRCPWCBlock) ) deallocate( Self.SRCPWCBlock )
    if ( allocated(Self.PWCLocBlock) ) deallocate( Self.PWCLocBlock )
    if ( allocated(Self.PWCSRCBlock) ) deallocate( Self.PWCSRCBlock )
    if ( allocated(Self.PWCPWCBlock) ) deallocate( Self.PWCPWCBlock )
    !
    allocate( Self.LocLocBlock(NumLocBra,NumLocKet) )
    allocate( Self.LocSRCBlock(NumLocBra,NumSRCKet) )
    allocate( Self.LocPWCBlock(NumLocBra,NumPWCKet) )
    allocate( Self.SRCLocBlock(NumSRCBra,NumLocKet) )
    allocate( Self.SRCSRCBlock(NumSRCBra,NumSRCKet) )
    allocate( Self.SRCPWCBlock(NumSRCBra,NumPWCKet) )
    allocate( Self.PWCLocBlock(NumPWCBra,NumLocKet) )
    allocate( Self.PWCSRCBlock(NumPWCBra,NumSRCKet) )
    allocate( Self.PWCPWCBlock(NumPWCBra,NumPWCKet) )
    !
    !
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    !
    !> A generic one-body operator between CC SI states has the form
    !! \f[
    !!     \langle \Psi_{A\alpha} | O | \Psi_{B\beta} \rangle = 
    !!       \delta_{\alpha\beta} \sum_{nm} o_{nm}\rho^{BA}_mn 
    !!     + \delta_{AB} o_{\alpha\beta}
    !!     - \sum_\delta \rho^{BA}_{\alpha\delta}o_{\delta\beta}
    !!     - \sum_\delta o_{\alpha\delta}\rho^{BA}_{\delta\beta}
    !! \f]
    !!
    !> The e-e two-body operator between CC SI states has the form
    !! \f[
    !!     \langle \Psi_{A\alpha} | G | \Psi_{B\beta} \rangle = 
    !!       \delta_{\alpha\beta} < A | G | B >
    !!     +   \sum_{nm}    \rho^{AB}_{nm}    { [ n m | \beta \alpha ] - [ n \alpha | m \beta ] } 
    !!     - 2 \sum_{nop}  \pi^{BA}_{n\beta,op} [ n o | p \alpha ] 
    !!     - 2 \sum_{nmo}  \pi^{BA}_{nm,o\alpha} [ n o | m \beta ] 
    !! \f]
    !!
    !
    !.. < * | O | loc >
    do j = 1, NumLocKet

       KetLoc  => Self.SpaceKet.GetChannelLS()

       !.. < loc | O | loc >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly    ) call Self.LocLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite   ) call Self.LocLocBlock(i,j).SetFormattedWrite()
          call Self.LocLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocLocBlock(i,j).Build( BraLoc, KetLoc, OpLabel, OpXlm, Force )
          call Self.LocLocBlock(i,j).Save( )
          call Self.LocLocBlock(i,j).Free( )
       end do
       
       !.. Sum over the parent ions
       do i = 1, NumSRCBra

          !.. < src | O | loc >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly    ) call Self.SRCLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite   ) call Self.SRCLocBlock(i,j).SetFormattedWrite()
          call Self.SRCLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCLocBlock(i,j).Build( BraSRC, KetLoc, OpLabel, OpXlm, Force )
          call Self.SRCLocBlock(i,j).Save( )
          call Self.SRCLocBlock(i,j).Free( )
          
          !.. < pwc | O | loc >
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             
             !.. < pwc Xlm | O | loc >
             BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
             CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
             if ( Self.BoxStatesOnly    ) call Self.PWCLocBlock(CounterBra,j).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetKetBoxOnly()
             if ( Self.FormattedWrite   ) call Self.PWCLocBlock(CounterBra,j).SetFormattedWrite()
             call Self.PWCLocBlock(CounterBra,j).SetFileExtension( FileExtensionPWCLoc  ) 
             call Self.PWCLocBlock(CounterBra,j).Build( BraPWC, KetLoc, OpLabel, OpXlm, Force )
             if ( present(ExtraLabel) ) then
                call Self.PWCLocBlock(CounterBra,j).SetFileExtension( FileExtensionPWCLoc//ExtraLabel )
             end if
             call Self.PWCLocBlock(CounterBra,j).Save( )
             call Self.PWCLocBlock(CounterBra,j).Free( )
             
          end do

       end do

    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       
       KetSRC => Self.SpaceKet.GetChannelSRC(j)

       !.. < loc | O | src >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly    ) call Self.LocSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite   ) call Self.LocSRCBlock(i,j).SetFormattedWrite()
          call Self.LocSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocSRCBlock(i,j).Build( BraLoc, KetSRC, OpLabel, OpXlm, Force )
          call Self.LocSRCBlock(i,j).Save( )
          call Self.LocSRCBlock(i,j).Free( )
       end do

       do i = 1, NumSRCBra

          !.. < src | O | src >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly    ) call Self.SRCSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite   ) call Self.SRCSRCBlock(i,j).SetFormattedWrite()
          call Self.SRCSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCSRCBlock(i,j).Build( BraSRC, KetSRC, OpLabel, OpXlm, Force )
          call Self.SRCSRCBlock(i,j).Save( )
          call Self.SRCSRCBlock(i,j).Free( )
          
          do n = 1, Self.SpaceBra.GetNumPWC(i)

             !.. < pwc Xlm | O | src >
             BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
             CounterBra =  Self.SpaceBra.PWCChannelIndex(i,n)
             if ( Self.BoxStatesOnly    ) call Self.PWCSRCBlock(CounterBra,j).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetKetBoxOnly()
             if ( Self.FormattedWrite   ) call Self.PWCSRCBlock(CounterBra,j).SetFormattedWrite()
             call Self.PWCSRCBlock(CounterBra,j).SetFileExtension( FileExtensionPWCSRC  ) 
             call Self.PWCSRCBlock(CounterBra,j).Build( BraPWC, KetSRC, OpLabel, OpXlm, Force )
             if ( present(ExtraLabel) ) then
                call Self.PWCSRCBlock(CounterBra,j).SetFileExtension( FileExtensionPWCSRC//ExtraLabel  ) 
             end if
             call Self.PWCSRCBlock(CounterBra,j).Save( )
             call Self.PWCSRCBlock(CounterBra,j).Free( )
             
          end do
          
       end do

    end do


    !.. < * | O | pwc >
    do j = 1, NumSRCKet
       do k = 1, Self.SpaceKet.GetNumPWC(j)
          KetPWC => Self.SpaceKet.GetChannelPWC(j,k)
          CounterKet = Self.SpaceKet.PWCChannelIndex(j,k)

          !.. < loc | O | pwc Xlm >
          do i = 1, NumLocBra
             BraLoc  => Self.SpaceBra.GetChannelLS()
             if ( Self.BoxStatesOnly    ) call Self.LocPWCBlock(i,CounterKet).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetKetBoxOnly()
             if ( Self.FormattedWrite   ) call Self.LocPWCBlock(i,CounterKet).SetFormattedWrite()
             call Self.LocPWCBlock(i,CounterKet).SetFileExtension( FileExtensionLocPWC  ) 
             call Self.LocPWCBlock(i,CounterKet).Build( BraLoc, KetPWC, OpLabel, OpXlm, Force )
             if ( present(ExtraLabel) ) then
                call Self.LocPWCBlock(i,CounterKet).SetFileExtension( FileExtensionLocPWC//ExtraLabel ) 
             end if
             call Self.LocPWCBlock(i,CounterKet).Save( )
             call Self.LocPWCBlock(i,CounterKet).Free( )
          end do
          
          !.. < src | O | pwc Xlm >
          do i = 1, NumSRCBra
             BraSRC => Self.SpaceBra.GetChannelSRC(i)
             if ( Self.BoxStatesOnly    ) call Self.SRCPWCBlock(i,CounterKet).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetKetBoxOnly()
             if ( Self.FormattedWrite   ) call Self.SRCPWCBlock(i,CounterKet).SetFormattedWrite()
             call Self.SRCPWCBlock(i,CounterKet).SetFileExtension( FileExtensionSRCPWC  ) 
             call Self.SRCPWCBlock(i,CounterKet).Build( BraSRC, KetPWC, OpLabel, OpXlm, Force )
             if ( present(ExtraLabel) ) then
                call Self.SRCPWCBlock(i,CounterKet).SetFileExtension( FileExtensionSRCPWC//ExtraLabel ) 
             end if
             call Self.SRCPWCBlock(i,CounterKet).Save( )
             call Self.SRCPWCBlock(i,CounterKet).Free( )
             
             do n = 1, Self.SpaceBra.GetNumPWC(i)

                !.. < pwc Xlm | O | pwc Xlm >
                BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                if ( Self.BoxStatesOnly    ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetBoxOnly()
                if ( Self.BraBoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetBraBoxOnly()
                if ( Self.KetBoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetKetBoxOnly()
                if ( Self.FormattedWrite   ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetFormattedWrite()
                call Self.PWCPWCBlock(CounterBra,CounterKet).SetFileExtension( FileExtensionPWCPWC  ) 
                call Self.PWCPWCBlock(CounterBra,CounterKet).Build( BraPWC, KetPWC, OpLabel, OpXlm, Force )
                if ( present(ExtraLabel) ) then
                   call Self.PWCPWCBlock(CounterBra,CounterKet).SetFileExtension( FileExtensionPWCPWC//ExtraLabel ) 
                end if
                call Self.PWCPWCBlock(CounterBra,CounterKet).Save( )
                call Self.PWCPWCBlock(CounterBra,CounterKet).Free( )
                
             end do
          end do

       end do
    end do
    !
  end subroutine ClassSESSESBlockInitAndSave






  !> Reads the ClassSESSESBlock class.
  subroutine ClassSESSESBlockLoad( Self, SpaceBra, SpaceKet, OpLabel, Axis, ExtraLabel )
    !> Class to be initialized.
    class(ClassSESSESBlock),                  intent(inout) :: Self
    !> Space of the Bra functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceBra
    !> Space of the Ket functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceKet
    !> Operator Label.
    character(len=*),                            intent(in)    :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable              , intent(inout) :: Axis
    !> Extra label to form the block name
    character(len=*), optional                 , intent(in)    :: ExtraLabel
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumPWCBra, NumPWCKet, NumLocBra, NumLocKet
    type(ClassPartialWaveChannel)       , pointer :: BraPWC, KetPWC
    type(ClassShortRangeChannel)        , pointer :: BraSRC, KetSRC
    type(ClassSymmetricLocalizedStates) , pointer :: BraLoc, KetLoc
    type(ClassXlm) :: OpXlm
    !
    Self.SpaceBra => SpaceBra
    Self.SpaceKet => SpaceKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    NumPWCBra = Self.SpaceBra.GetNumPWC()
    NumPWCKet = Self.SpaceKet.GetNumPWC()
    !
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    if ( allocated(Self.LocSRCBlock) ) deallocate( Self.LocSRCBlock )
    if ( allocated(Self.LocPWCBlock) ) deallocate( Self.LocPWCBlock )
    if ( allocated(Self.SRCLocBlock) ) deallocate( Self.SRCLocBlock )
    if ( allocated(Self.SRCSRCBlock) ) deallocate( Self.SRCSRCBlock )
    if ( allocated(Self.SRCPWCBlock) ) deallocate( Self.SRCPWCBlock )
    if ( allocated(Self.PWCLocBlock) ) deallocate( Self.PWCLocBlock )
    if ( allocated(Self.PWCSRCBlock) ) deallocate( Self.PWCSRCBlock )
    if ( allocated(Self.PWCPWCBlock) ) deallocate( Self.PWCPWCBlock )
    !
    allocate( Self.LocLocBlock(NumLocBra,NumLocKet) )
    allocate( Self.LocSRCBlock(NumLocBra,NumSRCKet) )
    allocate( Self.LocPWCBlock(NumLocBra,NumPWCKet) )
    allocate( Self.SRCLocBlock(NumSRCBra,NumLocKet) )
    allocate( Self.SRCSRCBlock(NumSRCBra,NumSRCKet) )
    allocate( Self.SRCPWCBlock(NumSRCBra,NumPWCKet) )
    allocate( Self.PWCLocBlock(NumPWCBra,NumLocKet) )
    allocate( Self.PWCSRCBlock(NumPWCBra,NumSRCKet) )
    allocate( Self.PWCPWCBlock(NumPWCBra,NumPWCKet) )
    !
    !
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    !
    !.. < * | O | loc >
    do j = 1, NumLocKet

       KetLoc  => Self.SpaceKet.GetChannelLS()

       !.. < loc | O | loc >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocLocBlock(i,j).SetFormattedWrite()
          call Self.LocLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocLocBlock(i,j).Load( BraLoc, KetLoc, OpLabel, OpXlm )
       end do
       
       !.. Sum over the parent ions
       do i = 1, NumSRCBra

          !.. < src | O | loc >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCLocBlock(i,j).SetFormattedWrite()
          call Self.SRCLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCLocBlock(i,j).Load( BraSRC, KetLoc, OpLabel, OpXlm )
          
          !.. < pwc | O | loc >
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             
             !.. < pwc Xlm | O | loc >
             BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
             CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
             if ( Self.BoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.PWCLocBlock(CounterBra,j).SetKetBoxOnly()
             if ( Self.AsympBsPermutation ) call Self.PWCLocBlock(CounterBra,j).SetAsympBsPermutation()
             if ( Self.FormattedWrite ) call Self.PWCLocBlock(CounterBra,j).SetFormattedWrite()
             if ( present(ExtraLabel) ) then
                call Self.PWCLocBlock(CounterBra,j).SetFileExtension( FileExtensionPWCLoc//ExtraLabel )
             else
                call Self.PWCLocBlock(CounterBra,j).SetFileExtension( FileExtensionPWCLoc  )
             end if
             call Self.PWCLocBlock(CounterBra,j).Load( BraPWC, KetLoc, OpLabel, OpXlm )
             
          end do

       end do

    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       
       KetSRC => Self.SpaceKet.GetChannelSRC(j)

       !.. < loc | O | src >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocSRCBlock(i,j).SetFormattedWrite()
          call Self.LocSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocSRCBlock(i,j).Load( BraLoc, KetSRC, OpLabel, OpXlm )
       end do

       do i = 1, NumSRCBra

          !.. < src | O | src >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCSRCBlock(i,j).SetFormattedWrite()
          call Self.SRCSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCSRCBlock(i,j).Load( BraSRC, KetSRC, OpLabel, OpXlm )
          
          do n = 1, Self.SpaceBra.GetNumPWC(i)

             !.. < pwc Xlm | O | src >
             BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
             CounterBra =  Self.SpaceBra.PWCChannelIndex(i,n)
             if ( Self.BoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.PWCSRCBlock(CounterBra,j).SetKetBoxOnly()
             if ( Self.AsympBsPermutation ) call Self.PWCSRCBlock(CounterBra,j).SetAsympBsPermutation()
             if ( Self.FormattedWrite ) call Self.PWCSRCBlock(CounterBra,j).SetFormattedWrite()
             if ( present(ExtraLabel) ) then
                call Self.PWCSRCBlock(CounterBra,j).SetFileExtension( FileExtensionPWCSRC//ExtraLabel  ) 
             else
                call Self.PWCSRCBlock(CounterBra,j).SetFileExtension( FileExtensionPWCSRC  ) 
             end if
             call Self.PWCSRCBlock(CounterBra,j).Load( BraPWC, KetSRC, OpLabel, OpXlm )
             
          end do
          
       end do

    end do


    !.. < * | O | pwc >
    do j = 1, NumSRCKet
       do k = 1, Self.SpaceKet.GetNumPWC(j)
          KetPWC => Self.SpaceKet.GetChannelPWC(j,k)
          CounterKet = Self.SpaceKet.PWCChannelIndex(j,k)

          !.. < loc | O | pwc Xlm >
          do i = 1, NumLocBra
             BraLoc  => Self.SpaceBra.GetChannelLS()
             if ( Self.BoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.LocPWCBlock(i,CounterKet).SetKetBoxOnly()
             if ( Self.AsympBsPermutation ) call Self.LocPWCBlock(i,CounterKet).SetAsympBsPermutation()
             if ( Self.FormattedWrite ) call Self.LocPWCBlock(i,CounterKet).SetFormattedWrite()
             if ( present(ExtraLabel) ) then
                call Self.LocPWCBlock(i,CounterKet).SetFileExtension( FileExtensionLocPWC//ExtraLabel ) 
             else
                call Self.LocPWCBlock(i,CounterKet).SetFileExtension( FileExtensionLocPWC ) 
             end if
             call Self.LocPWCBlock(i,CounterKet).Load( BraLoc, KetPWC, OpLabel, OpXlm )
          end do
          
          !.. < src | O | pwc Xlm >
          do i = 1, NumSRCBra
             BraSRC => Self.SpaceBra.GetChannelSRC(i)
             if ( Self.BoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetBoxOnly()
             if ( Self.BraBoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetBraBoxOnly()
             if ( Self.KetBoxStatesOnly ) call Self.SRCPWCBlock(i,CounterKet).SetKetBoxOnly()
             if ( Self.AsympBsPermutation ) call Self.SRCPWCBlock(i,CounterKet).SetAsympBsPermutation()
             if ( Self.FormattedWrite ) call Self.SRCPWCBlock(i,CounterKet).SetFormattedWrite()
             if ( present(ExtraLabel) ) then
                call Self.SRCPWCBlock(i,CounterKet).SetFileExtension( FileExtensionSRCPWC//ExtraLabel ) 
             else
                call Self.SRCPWCBlock(i,CounterKet).SetFileExtension( FileExtensionSRCPWC ) 
             end if
             call Self.SRCPWCBlock(i,CounterKet).Load( BraSRC, KetPWC, OpLabel, OpXlm )
             
             do n = 1, Self.SpaceBra.GetNumPWC(i)

                !.. < pwc Xlm | O | pwc Xlm >
                BraPWC => Self.SpaceBra.GetChannelPWC(i,n)
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                if ( Self.BoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetBoxOnly()
                if ( Self.BraBoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetBraBoxOnly()
                if ( Self.KetBoxStatesOnly ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetKetBoxOnly()
                if ( Self.AsympBsPermutation ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetAsympBsPermutation()
                if ( Self.FormattedWrite ) call Self.PWCPWCBlock(CounterBra,CounterKet).SetFormattedWrite()
                if ( present(ExtraLabel) ) then
                   call Self.PWCPWCBlock(CounterBra,CounterKet).SetFileExtension( FileExtensionPWCPWC//ExtraLabel ) 
                else
                   call Self.PWCPWCBlock(CounterBra,CounterKet).SetFileExtension( FileExtensionPWCPWC ) 
                end if
                call Self.PWCPWCBlock(CounterBra,CounterKet).Load( BraPWC, KetPWC, OpLabel, OpXlm )
                
             end do
          end do

       end do
    end do


    !
  end subroutine ClassSESSESBlockLoad



  !> Reads the QC components in ClassSESSESBlock class.
  subroutine ClassSESSESBlockLoadQC( Self, SpaceBra, SpaceKet, OpLabel, Axis, ExtraLabel )
    !> Class to be initialized.
    class(ClassSESSESBlock),                  intent(inout) :: Self
    !> Space of the Bra functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceBra
    !> Space of the Ket functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceKet
    !> Operator Label.
    character(len=*),                            intent(in)    :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable              , intent(inout) :: Axis
    !> Extra label to form the block name
    character(len=*), optional                 , intent(in)    :: ExtraLabel
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumLocBra, NumLocKet
    type(ClassShortRangeChannel)        , pointer :: BraSRC, KetSRC
    type(ClassSymmetricLocalizedStates) , pointer :: BraLoc, KetLoc
    type(ClassXlm) :: OpXlm
    !
    Self.SpaceBra => SpaceBra
    Self.SpaceKet => SpaceKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    if ( allocated(Self.LocSRCBlock) ) deallocate( Self.LocSRCBlock )
    if ( allocated(Self.SRCLocBlock) ) deallocate( Self.SRCLocBlock )
    if ( allocated(Self.SRCSRCBlock) ) deallocate( Self.SRCSRCBlock )
    !
    allocate( Self.LocLocBlock(NumLocBra,NumLocKet) )
    allocate( Self.LocSRCBlock(NumLocBra,NumSRCKet) )
    allocate( Self.SRCLocBlock(NumSRCBra,NumLocKet) )
    allocate( Self.SRCSRCBlock(NumSRCBra,NumSRCKet) )
    !
    !
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    !
    !.. < * | O | loc >
    do j = 1, NumLocKet

       KetLoc  => Self.SpaceKet.GetChannelLS()

       !.. < loc | O | loc >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocLocBlock(i,j).SetFormattedWrite()
          call Self.LocLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocLocBlock(i,j).Load( BraLoc, KetLoc, OpLabel, OpXlm )
       end do
       
       !.. Sum over the parent ions
       do i = 1, NumSRCBra

          !.. < src | O | loc >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCLocBlock(i,j).SetFormattedWrite()
          call Self.SRCLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCLocBlock(i,j).Load( BraSRC, KetLoc, OpLabel, OpXlm )
          
       end do

    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       
       KetSRC => Self.SpaceKet.GetChannelSRC(j)

       !.. < loc | O | src >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocSRCBlock(i,j).SetFormattedWrite()
          call Self.LocSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocSRCBlock(i,j).Load( BraLoc, KetSRC, OpLabel, OpXlm )
       end do

       do i = 1, NumSRCBra

          !.. < src | O | src >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCSRCBlock(i,j).SetFormattedWrite()
          call Self.SRCSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCSRCBlock(i,j).Load( BraSRC, KetSRC, OpLabel, OpXlm )
          
       end do

    end do


    !
  end subroutine ClassSESSESBlockLoadQC



  !> Reads the polycentric Gaussian components in ClassSESSESBlock class.
  subroutine ClassSESSESBlockLoadPoly( Self, SpaceBra, SpaceKet, OpLabel, Axis, ExtraLabel )
    !> Class to be initialized.
    class(ClassSESSESBlock),                  intent(inout) :: Self
    !> Space of the Bra functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceBra
    !> Space of the Ket functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceKet
    !> Operator Label.
    character(len=*),                            intent(in)    :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable              , intent(inout) :: Axis
    !> Extra label to form the block name
    character(len=*), optional                 , intent(in)    :: ExtraLabel
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumLocBra, NumLocKet
    type(ClassShortRangeChannel)        , pointer :: BraSRC, KetSRC
    type(ClassSymmetricLocalizedStates) , pointer :: BraLoc, KetLoc
    type(ClassXlm) :: OpXlm
    logical, parameter :: OnlyPoly = .true.
    !
    Self.SpaceBra => SpaceBra
    Self.SpaceKet => SpaceKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    if ( allocated(Self.LocSRCBlock) ) deallocate( Self.LocSRCBlock )
    if ( allocated(Self.SRCLocBlock) ) deallocate( Self.SRCLocBlock )
    if ( allocated(Self.SRCSRCBlock) ) deallocate( Self.SRCSRCBlock )
    !
    allocate( Self.LocLocBlock(NumLocBra,NumLocKet) )
    allocate( Self.LocSRCBlock(NumLocBra,NumSRCKet) )
    allocate( Self.SRCLocBlock(NumSRCBra,NumLocKet) )
    allocate( Self.SRCSRCBlock(NumSRCBra,NumSRCKet) )
    !
    !
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    !
    !.. < * | O | loc >
    do j = 1, NumLocKet

       KetLoc  => Self.SpaceKet.GetChannelLS()

       !.. < loc | O | loc >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocLocBlock(i,j).SetFormattedWrite()
          call Self.LocLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocLocBlock(i,j).Load( BraLoc, KetLoc, OpLabel, OpXlm )
       end do
       
       !.. Sum over the parent ions
       do i = 1, NumSRCBra

          !.. < src | O | loc >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCLocBlock(i,j).SetFormattedWrite()
          call Self.SRCLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCLocBlock(i,j).Load( BraSRC, KetLoc, OpLabel, OpXlm, OnlyPoly )
          
       end do

    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       
       KetSRC => Self.SpaceKet.GetChannelSRC(j)

       !.. < loc | O | src >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocSRCBlock(i,j).SetFormattedWrite()
          call Self.LocSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocSRCBlock(i,j).Load( BraLoc, KetSRC, OpLabel, OpXlm, OnlyPoly )
       end do

       do i = 1, NumSRCBra

          !.. < src | O | src >
          BraSRC => Self.SpaceBra.GetChannelSRC(i)
          if ( Self.BoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.SRCSRCBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.SRCSRCBlock(i,j).SetFormattedWrite()
          call Self.SRCSRCBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.SRCSRCBlock(i,j).Load( BraSRC, KetSRC, OpLabel, OpXlm, OnlyPoly )
          
       end do

    end do


    !
  end subroutine ClassSESSESBlockLoadPoly



  !> Reads the localized components in ClassSESSESBlock class.
  !! In general just the polycentric submatrix.
  subroutine ClassSESSESBlockLoadLoc( Self, SpaceBra, SpaceKet, OpLabel, Axis, ExtraLabel )
    !> Class to be initialized.
    class(ClassSESSESBlock),                  intent(inout) :: Self
    !> Space of the Bra functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceBra
    !> Space of the Ket functions.
    type(ClassSymmetricElectronicSpace), target, intent(in)    :: SpaceKet
    !> Operator Label.
    character(len=*),                            intent(in)    :: OpLabel
    !> Dipole orientation.
    character(len=:), allocatable              , intent(inout) :: Axis
    !> Extra label to form the block name
    character(len=*), optional                 , intent(in)    :: ExtraLabel
    !
    integer :: i, j, k, n
    integer :: NumLocBra, NumLocKet
    type(ClassSymmetricLocalizedStates) , pointer :: BraLoc, KetLoc
    type(ClassXlm) :: OpXlm
    !
    Self.SpaceBra => SpaceBra
    Self.SpaceKet => SpaceKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !
    if ( allocated(Self.LocLocBlock) ) deallocate( Self.LocLocBlock )
    !
    allocate( Self.LocLocBlock(NumLocBra,NumLocKet) )
    !
    !
    OpXlm = GetOperatorXlm(OpLabel,Axis)
    !
    !.. < * | O | loc >
    do j = 1, NumLocKet
       
       KetLoc  => Self.SpaceKet.GetChannelLS()

       !.. < loc | O | loc >
       do i = 1, NumLocBra
          BraLoc  => Self.SpaceBra.GetChannelLS()
          if ( Self.BoxStatesOnly ) call Self.LocLocBlock(i,j).SetBoxOnly()
          if ( Self.BraBoxStatesOnly ) call Self.LocLocBlock(i,j).SetBraBoxOnly()
          if ( Self.KetBoxStatesOnly ) call Self.LocLocBlock(i,j).SetKetBoxOnly()
          if ( Self.FormattedWrite ) call Self.LocLocBlock(i,j).SetFormattedWrite()
          call Self.LocLocBlock(i,j).SetFileExtension( GetQCFileExtension() ) 
          call Self.LocLocBlock(i,j).Load( BraLoc, KetLoc, OpLabel, OpXlm )
       end do
       
    end do
    
    !
  end subroutine ClassSESSESBlockLoadLoc




  !> Save the ClassSESSESBlock class.
  subroutine ClassSESSESBlockSave( Self, ExtraLabel )
    !> Class to be initialized.
    class(ClassSESSESBlock)   , intent(inout) :: Self
    character(len=*), optional, intent(in)    :: ExtraLabel
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumPWCBra, NumPWCKet, NumLocBra, NumLocKet
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    NumPWCBra = Self.SpaceBra.GetNumPWC()
    NumPWCKet = Self.SpaceKet.GetNumPWC()
    !
    if ( .not.allocated(Self.LocLocBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.LocSRCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.LocPWCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.SRCLocBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.SRCSRCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.SRCPWCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.PWCLocBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.PWCSRCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.PWCPWCBlock) ) call Assert( 'Sub-block not allocated.' )
    !

    !.. < * | O | loc >
    do j = 1, NumLocKet
       !.. < loc | O | loc >
       do i = 1, NumLocBra
          call Self.LocLocBlock(i,j).Save( )
       end do
       !.. Sum over the parent ions
       do i = 1, NumSRCBra
          !.. < src | O | loc >
          call Self.SRCLocBlock(i,j).Save( )
          !.. < pwc | O | loc >
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             !.. < pwc Xlm | O | loc >
             CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
             if ( present(ExtraLabel) ) then
                call Self.PWCLocBlock(CounterBra,j).SetFileExtension( FileExtensionPWCLoc//ExtraLabel )
             end if
             call Self.PWCLocBlock(CounterBra,j).Save( )
          end do
       end do
    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       !.. < loc | O | src >
       do i = 1, NumLocBra
          call Self.LocSRCBlock(i,j).Save( )
       end do
       do i = 1, NumSRCBra
          !.. < src | O | src >
          call Self.SRCSRCBlock(i,j).Save( )
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             !.. < pwc Xlm | O | src >
             CounterBra =  Self.SpaceBra.PWCChannelIndex(i,n)
             if ( present(ExtraLabel) ) then
                call Self.PWCSRCBlock(CounterBra,j).SetFileExtension( FileExtensionPWCSRC//ExtraLabel  ) 
             end if
             call Self.PWCSRCBlock(CounterBra,j).Save( )
          end do
       end do
    end do


    !.. < * | O | pwc >
    do j = 1, NumSRCKet
       do k = 1, Self.SpaceKet.GetNumPWC(j)
          CounterKet = Self.SpaceKet.PWCChannelIndex(j,k)
          !.. < loc | O | pwc Xlm >
          do i = 1, NumLocBra
             if ( present(ExtraLabel) ) then
                call Self.LocPWCBlock(i,CounterKet).SetFileExtension( FileExtensionLocPWC//ExtraLabel ) 
             end if
             call Self.LocPWCBlock(i,CounterKet).Save( )
          end do
          !.. < src | O | pwc Xlm >
          do i = 1, NumSRCBra
             if ( present(ExtraLabel) ) then
                call Self.SRCPWCBlock(i,CounterKet).SetFileExtension( FileExtensionSRCPWC//ExtraLabel ) 
             end if
             call Self.SRCPWCBlock(i,CounterKet).Save( )
             do n = 1, Self.SpaceBra.GetNumPWC(i)
                !.. < pwc Xlm | O | pwc Xlm >
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                if ( present(ExtraLabel) ) then
                   call Self.PWCPWCBlock(CounterBra,CounterKet).SetFileExtension( FileExtensionPWCPWC//ExtraLabel ) 
                end if
                call Self.PWCPWCBlock(CounterBra,CounterKet).Save( )
             end do
          end do
       end do
    end do


    !
  end subroutine ClassSESSESBlockSave



  !> Gets the number of rows and columns of the full quantum chemistry
  !! operator matrices.
  subroutine  ClassSESSESBlockGetNumFunQC( Self, LoadLocStates, NBra, NKet ) 
    !> Class to be initialized.
    class(ClassSESSESBlock)   , intent(in)  :: Self
    logical                   , intent(in)  :: LoadLocStates
    integer                   , intent(out) :: NBra
    integer                   , intent(out) :: NKet
    !
    integer :: i, j
    integer :: NumSRCBra, NumSRCKet, NumLocBra, NumLocKet
    !
    NBra = 0
    NKet = 0
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()

    if ( LoadLocStates ) then
       !.. < * | O | loc >
       do j = 1, NumLocKet
          !.. < loc | O | loc >
          do i = 1, NumLocBra
             if ( j == 1 ) NBra = NBra + Self.LocLocBlock(i,j).GetNRows()
          end do
          NKet = NKet + Self.LocLocBlock(1,j).GetNColumns()
          !.. Sum over the parent ions
          do i = 1, NumSRCBra
             !.. < src | O | loc >
             if ( j == 1 ) NBra = NBra + Self.SRCLocBlock(i,j).GetNRows()
          end do
          NKet = NKet + Self.SRCLocBlock(1,j).GetNColumns()
       end do
    end if
    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       if ( LoadLocStates ) then
          !.. < loc | O | src >
          do i = 1, NumLocBra
             if ( j == 1 ) NBra = NBra + Self.LocSRCBlock(i,j).GetNRows()
          end do
          NKet = NKet + Self.LocSRCBlock(1,j).GetNColumns()
       end if
       do i = 1, NumSRCBra
          !.. < src | O | src >
          if ( j == 1 ) NBra = NBra + Self.SRCSRCBlock(i,j).GetNRows()
       end do
       NKet = NKet + Self.SRCSRCBlock(1,j).GetNColumns()
    end do
    !
  end subroutine ClassSESSESBlockGetNumFunQC



  !> Condition the blocks.
  subroutine ClassSESSESBlockCondition( Self, Conditioner )
    !> Class to be initialized.
    class(ClassSESSESBlock)     , intent(inout) :: Self
    class(ClassConditionerBlock), intent(inout) :: Conditioner
    !
    integer :: i, j, k, n, CounterKet, CounterBra
    integer :: NumSRCBra, NumSRCKet, NumPWCBra, NumPWCKet, NumLocBra, NumLocKet
    type(LocSRCBlock), allocatable :: OriginalLocSRCBlock(:,:)
    type(SRCLocBlock), allocatable :: OriginalSRCLocBlock(:,:)
    type(SRCSRCBlock), allocatable :: OriginalSRCSRCBlock(:,:)
    type(SRCPWCBlock), allocatable :: OriginalSRCPWCBlock(:,:)
    type(PWCSRCBlock), allocatable :: OriginalPWCSRCBlock(:,:)
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    NumPWCBra = Self.SpaceBra.GetNumPWC()
    NumPWCKet = Self.SpaceKet.GetNumPWC()
    !
    if ( .not.allocated(Self.LocPWCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.SRCPWCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.PWCLocBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.PWCSRCBlock) ) call Assert( 'Sub-block not allocated.' )
    if ( .not.allocated(Self.PWCPWCBlock) ) call Assert( 'Sub-block not allocated.' )
    !

    allocate( OriginalLocSRCBlock, source = Self.LocSRCBlock(:,:) )
    allocate( OriginalSRCLocBlock, source = Self.SRCLocBlock(:,:) )
    allocate( OriginalSRCSRCBlock, source = Self.SRCSRCBlock(:,:) )
    allocate( OriginalSRCPWCBlock, source = Self.SRCPWCBlock(:,:) )
    allocate( OriginalPWCSRCBlock, source = Self.PWCSRCBlock(:,:) )


    !.. < * | O | loc >
    do j = 1, NumLocKet
       !.. Sum over the parent ions
       do i = 1, NumSRCBra
          !.. < pwc | O | loc >
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             !.. < pwc Xlm | O | loc >
             CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
             call Self.PWCLocBlock(CounterBra,j).Condition( Conditioner, OriginalSRCLocBlock(i,j) )
          end do
       end do
    end do

    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       do i = 1, NumSRCBra
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             !.. < pwc Xlm | O | src >
             CounterBra =  Self.SpaceBra.PWCChannelIndex(i,n)
             call Self.PWCSRCBlock(CounterBra,j).Condition( Conditioner, OriginalSRCSRCBlock(i,j) )
          end do
       end do
    end do


    !.. < * | O | pwc >
    do j = 1, NumSRCKet
       do k = 1, Self.SpaceKet.GetNumPWC(j)
          CounterKet = Self.SpaceKet.PWCChannelIndex(j,k)
          !.. < loc | O | pwc Xlm >
          do i = 1, NumLocBra
             call Self.LocPWCBlock(i,CounterKet).Condition( Conditioner, OriginalLocSRCBlock(i,j) )
          end do
          !.. < src | O | pwc Xlm >
          do i = 1, NumSRCBra
             call Self.SRCPWCBlock(i,CounterKet).Condition( Conditioner, OriginalSRCSRCBlock(i,j) )
             do n = 1, Self.SpaceBra.GetNumPWC(i)
                !.. < pwc Xlm | O | pwc Xlm >
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                call Self.PWCPWCBlock(CounterBra,CounterKet).Condition( Conditioner, OriginalSRCSRCBlock(i,j), OriginalPWCSRCBlock(CounterBra,j), OriginalSRCPWCBlock(i,CounterKet) )
             end do
          end do
       end do
    end do


    !
  end subroutine ClassSESSESBlockCondition



  !> Gets the minimum number of not overlapping B-splines among all channels.
  integer function ClassSESSESBlockGetMinNumNotOvBs( Self ) result(NBs)
    !> Class to be initialized.
    class(ClassSESSESBlock)   , intent(in) :: Self
    !
    if (.not.allocated(Self.PWCSRCBlock) ) call Assert( &
         'The PWCSRC block is not allocated, impossible to get '//&
         'the minimum number of B-splines that do not overlap '//&
         'with the diffuse orbitals.' )
    !.. For every symmetric electronic space only one localized channel is present.
    NBs = GetMinNumNotOverlapingBs(Self.PWCSRCBlock,'Rows')
    !
  end function ClassSESSESBlockGetMinNumNotOvBs




  !> Assemble in one matrix all the blocks in ClassSESSESBlock class.
  subroutine ClassSESSESBlockAssemble( Self, LoadLocStates, OpMat )
    !> Class to be assemble.
    class(ClassSESSESBlock)   , intent(inout) :: Self
    logical                   , intent(in)    :: LoadLocStates
    !> Operator matrix in the symmetric electronic space.
    class(ClassMatrix)        , intent(out)   :: OpMat
    !
    integer :: i, j, k, n, CounterBra, CounterKet, LocIndex
    integer :: NumSRCBra, NumSRCKet, NumPWCBra, NumPWCKet, NumLocBra, NumLocKet
    type(ClassMatrix), allocatable :: LocLocMat(:,:), LocSRCMat(:,:), LocPWCMat(:,:)
    type(ClassMatrix), allocatable :: SRCLocMat(:,:), SRCSRCMat(:,:), SRCPWCMat(:,:)
    type(Classmatrix), allocatable :: PWCLocMat(:,:), PWCSRCMat(:,:), PWCPWCMat(:,:)
    type(ClassMatrix) :: FullLocLoc, FullLocSRC, FullLocPWC 
    type(ClassMatrix) :: FullSRCLoc, FullSRCSRC, FullSRCPWC 
    type(ClassMatrix) :: FullPWCLoc, FullPWCSRC, FullPWCPWC 
    type(ClassMatrix), allocatable :: ArrayMat(:,:)
    !
    if ( allocated(Self.LocLocBlock) ) then
       Self.AvailableBlocks = .true.
    else
       Self.AvailableBlocks = .false.
       return
    end if
    !
    !.. It's better for consistencies always compute the permutation
    ! of the decouple B-splines at the assemble step.
    if ( DoBsPermutation ) then
       call Self.SetAsympBsPermutation()
    end if
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !
    NumPWCBra = Self.SpaceBra.GetNumPWC()
    NumPWCKet = Self.SpaceKet.GetNumPWC()

    if ( LoadLocStates ) then
       allocate( LocLocMat(NumLocBra,NumLocKet) )
       allocate( LocSRCMat(NumLocBra,NumSRCKet) )
       allocate( LocPWCMat(NumLocBra,NumPWCKet) )
       allocate( SRCLocMat(NumSRCBra,NumLocKet) )
       allocate( PWCLocMat(NumPWCBra,NumLocKet) )
    end if
    allocate( SRCSRCMat(NumSRCBra,NumSRCKet) )
    allocate( SRCPWCMat(NumSRCBra,NumPWCKet) )
    allocate( PWCSRCMat(NumPWCBra,NumSRCKet) )
    allocate( PWCPWCMat(NumPWCBra,NumPWCKet) )

    if ( LoadLocStates ) then
       !.. < * | O | loc >
       do j = 1, NumLocKet
          !.. < loc | O | loc >
          do i = 1, NumLocBra
             LocLocMat(i,j) = Self.LocLocBlock(i,j).FetchBlock( )
          end do
          !.. Sum over the parent ions
          do i = 1, NumSRCBra
             !.. < src | O | loc >
             SRCLocMat(i,j) = Self.SRCLocBlock(i,j).FetchBlock( )
             !.. < pwc | O | loc >
             do n = 1, Self.SpaceBra.GetNumPWC(i)
                !.. < pwc Xlm | O | loc >
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                PWCLocMat(CounterBra,j) = Self.PWCLocBlock(CounterBra,j).FetchBlock( )
                if ( Self.PWCLocBlock(CounterBra,j).IsBoxOnly() ) call PWCLocMat(CounterBra,j).RemoveRows(1,'END')
                if ( Self.PWCLocBlock(CounterBra,j).IsBraBoxOnly() ) call PWCLocMat(CounterBra,j).RemoveRows(1,'END')
             end do
          end do
       end do
    end if
    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       if ( LoadLocStates ) then
          !.. < loc | O | src >
          do i = 1, NumLocBra
             LocSRCMat(i,j) = Self.LocSRCBlock(i,j).FetchBlock( )
          end do
       end if
       do i = 1, NumSRCBra
          !.. < src | O | src >
          SRCSRCMat(i,j) = Self.SRCSRCBlock(i,j).FetchBlock( )
          do n = 1, Self.SpaceBra.GetNumPWC(i)
             !.. < pwc Xlm | O | src >
             CounterBra =  Self.SpaceBra.PWCChannelIndex(i,n)
             PWCSRCMat(CounterBra,j) = Self.PWCSRCBlock(CounterBra,j).FetchBlock( )
             if ( Self.PWCSRCBlock(CounterBra,j).IsBoxOnly() ) call PWCSRCMat(CounterBra,j).RemoveRows(1,'END')
             if ( Self.PWCSRCBlock(CounterBra,j).IsBraBoxOnly() ) call PWCSRCMat(CounterBra,j).RemoveRows(1,'END')
          end do
       end do
    end do


    !.. < * | O | pwc >
    do j = 1, NumSRCKet
       do k = 1, Self.SpaceKet.GetNumPWC(j)
          CounterKet = Self.SpaceKet.PWCChannelIndex(j,k)
          if ( LoadLocStates ) then
             !.. < loc | O | pwc Xlm >
             do i = 1, NumLocBra
                LocPWCMat(i,CounterKet) = Self.LocPWCBlock(i,CounterKet).FetchBlock( )
                if ( Self.LocPWCBlock(i,CounterKet).IsBoxOnly() ) call LocPWCMat(i,CounterKet).RemoveColumns(1,'END')
                if ( Self.LocPWCBlock(i,CounterKet).IsKetBoxOnly() ) call LocPWCMat(i,CounterKet).RemoveColumns(1,'END')
             end do
          end if
          !.. < src | O | pwc Xlm >
          do i = 1, NumSRCBra
             SRCPWCMat(i,CounterKet) = Self.SRCPWCBlock(i,CounterKet).FetchBlock( )
             if ( Self.SRCPWCBlock(i,CounterKet).IsBoxOnly() ) call SRCPWCMat(i,CounterKet).RemoveColumns(1,'END')
             if ( Self.SRCPWCBlock(i,CounterKet).IsKetBoxOnly() ) call SRCPWCMat(i,CounterKet).RemoveColumns(1,'END')
             do n = 1, Self.SpaceBra.GetNumPWC(i)
                !.. < pwc Xlm | O | pwc Xlm >
                CounterBra = Self.SpaceBra.PWCChannelIndex(i,n)
                PWCPWCMat(CounterBra,CounterKet) = Self.PWCPWCBlock(CounterBra,CounterKet).FetchBlock( )
                if ( Self.PWCPWCBlock(CounterBra,CounterKet).IsBoxOnly() ) call PWCPWCMat(CounterBra,CounterKet).RemoveColumns(1,'END')
                if ( Self.PWCPWCBlock(CounterBra,CounterKet).IsKetBoxOnly() ) call PWCPWCMat(CounterBra,CounterKet).RemoveColumns(1,'END')
                if ( Self.PWCPWCBlock(CounterBra,CounterKet).IsBoxOnly() ) call PWCPWCMat(CounterBra,CounterKet).RemoveRows(1,'END')
                if ( Self.PWCPWCBlock(CounterBra,CounterKet).IsBraBoxOnly() ) call PWCPWCMat(CounterBra,CounterKet).RemoveRows(1,'END')
             end do
          end do
       end do
    end do


    if ( LoadLocStates ) then
       !
       call FullLocLoc.BuildUpMatrix( LocLocMat )
       deallocate( LocLocMat )
       call FullLocSRC.BuildUpMatrix( LocSRCMat )
       deallocate( LocSRCMat )
       call FullLocPWC.BuildUpMatrix( LocPWCMat )
       if ( Self.AsympBsPermutation ) then
          write(output_unit,*) "Permuting LocPWCMat's columns"
          call PermuteBsplines( FullLocPWC, LocPWCMat, Self.LocPWCBlock, 'Columns' )
       end if
       deallocate( LocPWCMat )
       call FullSRCLoc.BuildUpMatrix( SRCLocMat )
       deallocate( SRCLocMat )
       call FullPWCLoc.BuildUpMatrix( PWCLocMat )
       if ( Self.AsympBsPermutation ) then
          write(output_unit,*) "Permuting PWCLocMat's rows"
          call PermuteBsplines( FullPWCLoc, PWCLocMat, Self.PWCLocBlock, 'Rows' )
       end if
       !
       deallocate( PWCLocMat )
    end if
    !
    !
    call FullSRCSRC.BuildUpMatrix( SRCSRCMat )
    deallocate( SRCSRCMat )
    call FullSRCPWC.BuildUpMatrix( SRCPWCMat )
    if ( Self.AsympBsPermutation ) then
       write(output_unit,*) "Permuting SRCPWCMat's columns"
       call PermuteBsplines( FullSRCPWC, SRCPWCMat, Self.SRCPWCBlock, 'Columns' )
    end if
    deallocate( SRCPWCMat )
    call FullPWCSRC.BuildUpMatrix( PWCSRCMat )
    if ( Self.AsympBsPermutation ) then
       write(output_unit,*) "Permuting PWCSRCMat's rows"
       call PermuteBsplines( FullPWCSRC, PWCSRCMat, Self.PWCSRCBlock, 'Rows' )
    end if
    deallocate( PWCSRCMat )
    call FullPWCPWC.BuildUpMatrix( PWCPWCMat )
    if ( Self.AsympBsPermutation ) then
       write(output_unit,*) "Permuting PWCPWCMat's rows"
       call PermuteBsplines( FullPWCPWC, PWCPWCMat, Self.PWCPWCBlock, 'Rows' )
       write(output_unit,*) "Permuting PWCPWCMat's columns"
       call PermuteBsplines( FullPWCPWC, PWCPWCMat, Self.PWCPWCBlock, 'Columns' )
    end if
    !
    deallocate( PWCPWCMat )

    LocIndex = 1
    if ( LoadLocStates ) LocIndex = 0
    
    allocate( ArrayMat(3-LocIndex,3-LocIndex) )
    !
    if ( LoadLocStates ) then
       ArrayMat(1-LocIndex,1-LocIndex) = FullLocLoc
       call FullLocLoc.Free()
       ArrayMat(1-LocIndex,2-LocIndex) = FullLocSRC
       call FullLocSRC.Free()
       ArrayMat(1-LocIndex,3-LocIndex) = FullLocPWC
       call FullLocPWC.Free()
       !
       ArrayMat(2-LocIndex,1-LocIndex) = FullSRCLoc
       call FullSRCLoc.Free()
       !
       ArrayMat(3-LocIndex,1-LocIndex) = FullPWCLoc
       call FullPWCLoc.Free()
    end if
    !
    ArrayMat(2-LocIndex,2-LocIndex) = FullSRCSRC
    call FullSRCSRC.Free()
    ArrayMat(2-LocIndex,3-LocIndex) = FullSRCPWC
    call FullSRCPWC.Free()
    ArrayMat(3-LocIndex,2-LocIndex) = FullPWCSRC
    call FullPWCSRC.Free()
    ArrayMat(3-LocIndex,3-LocIndex) = FullPWCPWC
    call FullPWCPWC.Free()
    
    call OpMat.BuildUpMatrix( ArrayMat )
    deallocate( ArrayMat )

    !
  end subroutine ClassSESSESBlockAssemble





  !> Assemble in one matrix all the QC blocks in ClassSESSESBlock class.
  subroutine ClassSESSESBlockAssembleQC( Self, LoadLocStates, OpMat )
    !> Class to be initialized.
    class(ClassSESSESBlock)   , intent(inout) :: Self
    logical                   , intent(in)    :: LoadLocStates
    !> Space of the Bra functions.
    class(ClassMatrix)        , intent(out)   :: OpMat
    !
    integer :: i, j, LocIndex
    integer :: NumSRCBra, NumSRCKet, NumLocBra, NumLocKet
    type(ClassMatrix), allocatable :: LocLocMat(:,:), LocSRCMat(:,:), SRCLocMat(:,:), SRCSRCMat(:,:)
    type(ClassMatrix) :: FullLocLoc, FullLocSRC, FullSRCLoc, FullSRCSRC
    type(ClassMatrix), allocatable :: ArrayMat(:,:)
    !
    if ( allocated(Self.LocLocBlock) ) then
       Self.AvailableBlocks = .true.
    else
       Self.AvailableBlocks = .false.
       return
    end if
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !

    if ( LoadLocStates ) then
       allocate( LocLocMat(NumLocBra,NumLocKet) )
       allocate( LocSRCMat(NumLocBra,NumSRCKet) )
       allocate( SRCLocMat(NumSRCBra,NumLocKet) )
    end if
    allocate( SRCSRCMat(NumSRCBra,NumSRCKet) )

    if ( LoadLocStates ) then
       !.. < * | O | loc >
       do j = 1, NumLocKet
          !.. < loc | O | loc >
          do i = 1, NumLocBra
             LocLocMat(i,j) = Self.LocLocBlock(i,j).FetchBlock( )
          end do
          !.. Sum over the parent ions
          do i = 1, NumSRCBra
             !.. < src | O | loc >
             SRCLocMat(i,j) = Self.SRCLocBlock(i,j).FetchBlock( )
          end do
       end do
    end if
    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       if ( LoadLocStates ) then
          !.. < loc | O | src >
          do i = 1, NumLocBra
             LocSRCMat(i,j) = Self.LocSRCBlock(i,j).FetchBlock( )
          end do
       end if
       do i = 1, NumSRCBra
          !.. < src | O | src >
          SRCSRCMat(i,j) = Self.SRCSRCBlock(i,j).FetchBlock( )
       end do
    end do


    if ( LoadLocStates ) then
       !
       call FullLocLoc.BuildUpMatrix( LocLocMat )
       deallocate( LocLocMat )
       call FullLocSRC.BuildUpMatrix( LocSRCMat )
       deallocate( LocSRCMat )
       call FullSRCLoc.BuildUpMatrix( SRCLocMat )
       deallocate( SRCLocMat )
       !
    end if
    !
    call FullSRCSRC.BuildUpMatrix( SRCSRCMat )
    deallocate( SRCSRCMat )

    LocIndex = 1
    if ( LoadLocStates ) LocIndex = 0
    
    allocate( ArrayMat(2-LocIndex,2-LocIndex) )
    !
    if ( LoadLocStates ) then
       ArrayMat(1-LocIndex,1-LocIndex) = FullLocLoc
       call FullLocLoc.Free()
       ArrayMat(1-LocIndex,2-LocIndex) = FullLocSRC
       call FullLocSRC.Free()
       ArrayMat(2-LocIndex,1-LocIndex) = FullSRCLoc
       call FullSRCLoc.Free()
    end if
    !
    ArrayMat(2-LocIndex,2-LocIndex) = FullSRCSRC
    call FullSRCSRC.Free()
    
    call OpMat.BuildUpMatrix( ArrayMat )
    deallocate( ArrayMat )

    !
  end subroutine ClassSESSESBlockAssembleQC


  !> Assemble in one matrix all the polycentric Gaussian blocks in ClassSESSESBlock class.
  subroutine ClassSESSESBlockAssemblePoly( Self, LoadLocStates, OpMat )
    !> Class to be initialized.
    class(ClassSESSESBlock)   , intent(inout) :: Self
    logical                   , intent(in)    :: LoadLocStates
    !> Space of the Bra functions.
    class(ClassMatrix)        , intent(out)   :: OpMat
    !
    integer :: i, j, LocIndex
    integer :: NumSRCBra, NumSRCKet, NumLocBra, NumLocKet
    type(ClassMatrix), allocatable :: LocLocMat(:,:), LocSRCMat(:,:), SRCLocMat(:,:), SRCSRCMat(:,:)
    type(ClassMatrix) :: FullLocLoc, FullLocSRC, FullSRCLoc, FullSRCSRC
    type(ClassMatrix), allocatable :: ArrayMat(:,:)
    !
    if ( allocated(Self.LocLocBlock) ) then
       Self.AvailableBlocks = .true.
    else
       Self.AvailableBlocks = .false.
       return
    end if
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1
    !.. For every symmetric electronic channel only one short range channel is present.
    NumSRCBra = Self.SpaceBra.GetNumChannels()
    NumSRCKet = Self.SpaceKet.GetNumChannels()
    !

    if ( LoadLocStates ) then
       allocate( LocLocMat(NumLocBra,NumLocKet) )
       allocate( LocSRCMat(NumLocBra,NumSRCKet) )
       allocate( SRCLocMat(NumSRCBra,NumLocKet) )
    end if
    allocate( SRCSRCMat(NumSRCBra,NumSRCKet) )

    if ( LoadLocStates ) then
       !.. < * | O | loc >
       do j = 1, NumLocKet
          !.. < loc | O | loc >
          do i = 1, NumLocBra
             LocLocMat(i,j) = Self.LocLocBlock(i,j).FetchBlock( )
          end do
          !.. Sum over the parent ions
          do i = 1, NumSRCBra
             !.. < src | O | loc >
             SRCLocMat(i,j) = Self.SRCLocBlock(i,j).FetchBlock( )
          end do
       end do
    end if
    
    !.. < * | O | src >
    do j = 1, NumSRCKet
       if ( LoadLocStates ) then
          !.. < loc | O | src >
          do i = 1, NumLocBra
             LocSRCMat(i,j) = Self.LocSRCBlock(i,j).FetchBlock( )
          end do
       end if
       do i = 1, NumSRCBra
          !.. < src | O | src >
          SRCSRCMat(i,j) = Self.SRCSRCBlock(i,j).FetchBlock( )
       end do
    end do


    if ( LoadLocStates ) then
       !
       call FullLocLoc.BuildUpMatrix( LocLocMat )
       deallocate( LocLocMat )
       call FullLocSRC.BuildUpMatrix( LocSRCMat )
       deallocate( LocSRCMat )
       call FullSRCLoc.BuildUpMatrix( SRCLocMat )
       deallocate( SRCLocMat )
       !
    end if
    !
    call FullSRCSRC.BuildUpMatrix( SRCSRCMat )
    deallocate( SRCSRCMat )

    LocIndex = 1
    if ( LoadLocStates ) LocIndex = 0
    
    allocate( ArrayMat(2-LocIndex,2-LocIndex) )
    !
    if ( LoadLocStates ) then
       ArrayMat(1-LocIndex,1-LocIndex) = FullLocLoc
       call FullLocLoc.Free()
       ArrayMat(1-LocIndex,2-LocIndex) = FullLocSRC
       call FullLocSRC.Free()
       ArrayMat(2-LocIndex,1-LocIndex) = FullSRCLoc
       call FullSRCLoc.Free()
    end if
    !
    ArrayMat(2-LocIndex,2-LocIndex) = FullSRCSRC
    call FullSRCSRC.Free()
    
    call OpMat.BuildUpMatrix( ArrayMat )
    deallocate( ArrayMat )

    !
  end subroutine ClassSESSESBlockAssemblePoly


  !> Assemble in one matrix all the localized block in ClassSESSESBlock class. In general just the polycentric submatrix.
  subroutine ClassSESSESBlockAssembleLoc( Self, LoadLocStates, OpMat )
    !> Class to be initialized.
    class(ClassSESSESBlock)   , intent(inout) :: Self
    logical                   , intent(in)    :: LoadLocStates
    !> Space of the Bra functions.
    class(ClassMatrix)        , intent(out)   :: OpMat
    !
    integer :: i, j, LocIndex
    integer :: NumLocBra, NumLocKet
    type(ClassMatrix), allocatable :: LocLocMat(:,:)
    type(ClassMatrix) :: FullLocLoc
    type(ClassMatrix), allocatable :: ArrayMat(:,:)
    !
    if ( allocated(Self.LocLocBlock) ) then
       Self.AvailableBlocks = .true.
    else
       Self.AvailableBlocks = .false.
       return
    end if
    !
    !.. For every symmetric electronic space only one localized channel is present.
    NumLocBra = 1
    NumLocKet = 1

    if ( LoadLocStates ) then
       allocate( LocLocMat(NumLocBra,NumLocKet) )
    else
       call Assert( 'If it is requested the assembling of just '//&
            'the localized orbitals submatrix, then the '//&
            'variable "LoadLocStates" in the program '//&
            'configuration input file, must to be set .true.')
    end if

    !.. < * | O | loc >
    do j = 1, NumLocKet
       !.. < loc | O | loc >
       do i = 1, NumLocBra
          LocLocMat(i,j) = Self.LocLocBlock(i,j).FetchBlock( )
       end do
       !.. Sum over the parent ions
    end do
    !
    
    call FullLocLoc.BuildUpMatrix( LocLocMat )
    deallocate( LocLocMat )
    !

    allocate( ArrayMat(1,1) )
    !
    ArrayMat(1,1) = FullLocLoc
    call FullLocLoc.Free()
    !
    call OpMat.BuildUpMatrix( ArrayMat )
    deallocate( ArrayMat )
    !
    !
  end subroutine ClassSESSESBlockAssembleLoc



  function ClassSESSESBlockGetStorageDir( Self ) result(Dir)
    class(ClassSESSESBlock), intent(in) :: Self
    character(len=:), allocatable :: Dir
    !
    character(len=:), allocatable :: NucConfDir, CCDir, BraSymLabel, KetSymLabel
    !
    allocate( NucConfDir, source = Self.SpaceBra.GetStorageDir() )
    allocate( CCDir,      source = GetCloseCouplingDir() )
    allocate( BraSymLabel, source = Self.SpaceBra.GetLabel() )
    allocate( KetSymLabel, source = Self.SpaceKet.GetLabel() )
    !
    allocate( Dir, source = &
         AddSlash(NucConfDir)//&
         AddSlash(CCDir)//&
         BraSymLabel//'_'//AddSlash(KetSymLabel) )
    !
  end function ClassSESSESBlockGetStorageDir



end module ModuleSymmetricElectronicSpaceOperators
