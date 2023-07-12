module ModuleCCOSymSpace
  
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV
  use omp_lib

  use ModuleSystemUtils
  use ModuleErrorHandling
  use ModuleString
  use ModuleIO
  use ModuleCSFCCSymSpace
  use ModuleParentIons
  use ModuleParameterList
    
  implicit none

  private

  character(len=*), parameter :: lstrnUC="SPDFGHIKLMNOQRTUVWXYZ"
  character(len=*), parameter :: lstrn="spdfghiklmnoqrtuvwxyz"
  character(len=*), parameter :: pstrn="eo"

  character(len=*), parameter :: AUG_SUBDIR = "aug/"
  character(len=*), parameter :: CCO_SUBDIR = "cco/"
  character(len=*), parameter :: ION_SUBDIR = "ion/"
  character(len=*), parameter :: ONE_SUBDIR = "one/"

  character       , parameter :: CCO_CH_SEPARATOR  = ASCII_NEW_LINE
  character(len=*), parameter :: DEFAULT_SELECTION = "DEFAULT"
  
  integer         , parameter ::   LW_LOG_OFFSET =   0
  integer         , parameter ::    N_LOG_OFFSET =  LW_LOG_OFFSET + 2
  integer         , parameter ::  PAR_LOG_OFFSET =   N_LOG_OFFSET + 2
  integer         , parameter ::    L_LOG_OFFSET = PAR_LOG_OFFSET + 2
  integer         , parameter :: MULT_LOG_OFFSET =   L_LOG_OFFSET + 2

  type, public :: ClassCCOSymSpace
     
     private

     character(len=:), allocatable :: Name
     real(kind(1d0))               :: TotalCharge
     character(len=:), allocatable :: Selection

     !> Storage directory
     character(len=:), allocatable :: StoreDir
     character(len=3) :: term

     integer                       :: nCCPWC
     !
     !*** This set of data is a good candidate to branch out
     !    as a new CCPWC object
     !***
     character(len=8)  , allocatable :: CCPWC_Name(:)
     character(len=3)  , allocatable :: CCPWC_Term_Pion(:)
     integer           , allocatable :: CCPWC_iCIpion(:)
     real(kind(1d0))   , allocatable :: CCPWC_pionEval(:)
     integer           , allocatable :: CCPWC_lw(:)
     integer           , allocatable :: CCPWC_dim(:)
     integer           , allocatable :: CCPWC_idim(:)
     integer                         :: CCPWC_totdim
     character(len=256), allocatable :: CCPWC_Attr(:)
     
     integer                         :: CCPWC_NMAX_LOC
     real(kind(1d0))   , allocatable :: CCPWC_CAPSF(:)
     real(kind(1d0))   , allocatable :: CCPWC_ESHFT(:) ! Shift of the parent-ion energy

     !.. Auxiliary
     type(ClassCSFCCSymSpace), pointer     :: CSFCCSymSpace
     integer                 , allocatable :: CCPWC_listCSFPWC(:,:)

   contains
     
     generic, public :: Init          => ClassCCOSymSpaceInit, &
                                         ClassCCOSymSpaceInitFromConfigFile
     generic, public :: SetCSFSymSpace =>ClassCCOSymSpaceSetCSFSymSpace
     generic, public :: Show          => ClassCCOSymSpaceShow
     generic, public :: Get_TotalSize => ClassCCOSymSpaceGet_TotalSize
     generic, public :: GetChIonEnergy => ClassCCOSymSpaceGetChIonEnergy
     generic, public :: GetChOrbAngMom => ClassCCOSymSpaceGetChOrbAngMom
     generic, public :: GetIonCharge  => ClassCCOSymSpaceGetIonCharge
     generic, public :: GetTerm       => ClassCCOSymSpaceGetTerm
     generic, public :: GetChName     => ClassCCOSymSpaceGetChName
     generic, public :: GetNCh        => ClassCCOSymSpaceGetNCh
     generic, public :: GetChIndexes  => ClassCCOSymSpaceGetChIndexes
     generic, public :: GetCAPSF      => ClassCCOSymSpaceGetCAPSF
     generic, public :: GetESHFT      => ClassCCOSymSpaceGetESHFT
     generic, public :: GetL          => ClassCCOSymSpaceGetL
     generic, public :: GetDir        => ClassCCOSymSpaceGetDir
     generic, public :: GetBlock      => ClassCCOSymSpaceGetBlock
     generic, public :: GetSideBlock   => ClassCCOSymSpaceGetSideBlock
     generic, public :: SetBlock      => ClassCCOSymSpaceSetBlock
     generic, public :: GetOpenChannels => ClassCCOSymSpaceGetOpenChannels
     generic, public :: GetThresholdSequence => ClassCCOSymSpaceGetThresholdSequence
     generic, public :: SaveAllChProp => ClassCCOSymSpaceSaveAllChProp
     generic, public :: LoadAllChProp => ClassCCOSymSpaceLoadAllChProp
     generic, public :: SaveChProp    => ClassCCOSymSpaceSaveChProp
     generic, public :: LoadChProp    => ClassCCOSymSpaceLoadChProp
     generic, public :: ParseConfigFile => ClassCCOSymSpaceParseConfigFile
     generic, public :: ParseSection   => ClassCCOSymSpaceParseSection
     generic, public :: ParseSelection => ClassCCOSymSpaceParseSelection

     generic, public :: Transform_BraKet_CSF_to_CCO => ClassCCOSymSpaceTransform_BraKet_CSF_to_CCO
     generic, public :: Transform_Ket_CSF_to_CCO => ClassCCOSymSpaceTransform_Ket_CSF_to_CCO
     generic, public :: Transform_Bra_CSF_to_CCO => ClassCCOSymSpaceTransform_Bra_CSF_to_CCO

     generic, public :: EvalWeights   => ClassCCOSymSpaceEvalWeights

     generic, public :: CoreBlockName => ClassCCOSymSpaceCoreBlockName
     generic, public :: SaveMatrixInBlocks => ClassCCOSymSpaceSaveMatrixInBlocks
     generic, public :: LoadMatrixInBlocks => ClassCCOSymSpaceLoadMatrixInBlocks
     generic, public :: SaveMatrixInStripeBlocks => ClassCCOSymSpaceSaveMatrixInStripeBlocks
     
     !.. Private procedures

     procedure, private :: ClassCCOSymSpaceSetCSFSymSpace
     procedure, private :: ClassCCOSymSpaceGetThresholdSequence
     procedure, private :: ClassCCOSymSpaceGetL
     procedure, private :: ClassCCOSymSpaceGetOpenChannels
     procedure, private :: ClassCCOSymSpaceGetChIonEnergy
     procedure, private :: ClassCCOSymSpaceInit
     procedure, private :: ClassCCOSymSpaceInitFromConfigFile
     procedure, private :: ClassCCOSymSpaceGetChOrbAngMom
     procedure, private :: ClassCCOSymSpaceParseConfigFile
     procedure, private :: ClassCCOSymSpaceParseSection
     procedure, private :: ClassCCOSymSpaceParseSelection
     procedure, private :: ClassCCOSymSpaceShow
     procedure, private :: ClassCCOSymSpaceGet_TotalSize
     procedure, private :: ClassCCOSymSpaceGetTerm
     procedure, private :: ClassCCOSymSpaceGetIonCharge
     procedure, private :: ClassCCOSymSpaceGetChName
     procedure, private :: ClassCCOSymSpaceGetNCh
     procedure, private :: ClassCCOSymSpaceGetChIndexes
     procedure, private :: ClassCCOSymSpaceGetCAPSF
     procedure, private :: ClassCCOSymSpaceGetESHFT
     procedure, private :: ClassCCOSymSpaceTransform_Ket_CSF_to_CCO
     procedure, private :: ClassCCOSymSpaceTransform_Bra_CSF_to_CCO
     procedure, private :: ClassCCOSymSpaceTransform_BraKet_CSF_to_CCO
     procedure, private :: ClassCCOSymSpaceEvalWeights
     procedure, private :: ClassCCOSymSpaceCoreBlockName
     procedure, private :: ClassCCOSymSpaceGetDir
     procedure, private :: ClassCCOSymSpaceGetBlock
     procedure, private :: ClassCCOSymSpaceGetSideBlock
     procedure, private :: ClassCCOSymSpaceSetBlock
     procedure, private :: ClassCCOSymSpaceSaveMatrixInBlocks
     procedure, private :: ClassCCOSymSpaceSaveMatrixInStripeBlocks
     procedure, private :: ClassCCOSymSpaceLoadMatrixInBlocks
     procedure, private :: ClassCCOSymSpaceLoadAllChProp
     procedure, private :: ClassCCOSymSpaceSaveAllChProp
     procedure, private :: ClassCCOSymSpaceSaveChProp
     procedure, private :: ClassCCOSymSpaceLoadChProp

     ! generic  , private :: SetupDirTree => ClassCCOSymSpaceSetupDirTree
     ! generic  , private :: IonTermDir   => ClassCCOSymSpaceIonTermDir
     ! generic  , private :: AugTermDir   => ClassCCOSymSpaceAugTermDir


     ! procedure, private :: ClassCCOSymSpaceParseConfigFile

     ! procedure, private :: ClassCCOSymSpaceSetupDirTree
     ! procedure, private :: ClassCCOSymSpaceIonTermDir
     ! procedure, private :: ClassCCOSymSpaceAugTermDir

  end type ClassCCOSymSpace

  !.. Binary function (depend on two spaces)
  public :: Transform_SymBraKet_CSF_to_CCO
  public :: SymCCO_SaveMatrixInBlocks
  public :: SymCCO_LoadMatrixInBlocks
  public :: ExtractCCOSymBlock
  public :: ExtractCCOSymvec
  public :: SetCCOSymBlock
  public :: GetCCOpairDir

contains


  !.. Defines the attributes of the localized channel
  subroutine CCLOC_ParseAttrLst( AttrStrn, NMAX_LOC )
    character(len=*), intent(in)  :: AttrStrn
    type(ClassParameterList)      :: list
    integer, intent(out) :: NMAX_LOC
    call list%Add( "NMAX", 1000000, "optional" )
    call list%ParseStrn( AttrStrn )
    call list%Get( "NMAX", NMAX_LOC )
  end subroutine CCLOC_ParseAttrLst


  !.. Defines the attributes of the localized channel
  subroutine CCPWC_ParseAttrLst( AttrStrn, CAPSF, ESHFT )
    character(len=*), intent(in)  :: AttrStrn
    type(ClassParameterList)      :: list
    real(kind(1d0)),  intent(out) :: CAPSF 
    real(kind(1d0)),  intent(out) :: ESHFT !Energy shift
    call list%Add( "CAPSF", 0.d0, "optional" )
    call list%Add( "ESHFT", 0.d0, "optional" )
    call list%ParseStrn( AttrStrn )
    call list%Get( "CAPSF", CAPSF )
    call list%Get( "ESHFT", ESHFT )
    if(list%Present("CAPSF")) write(*,*) "CAPSF =",CAPSF
    if(list%Present("ESHFT")) write(*,*) "ESHFT =",ESHFT
  end subroutine CCPWC_ParseAttrLst


  subroutine ExtractCCOSymvec( KetSpace, A, ik, blk )
    type(ClassCCOSymSpace), intent(in) :: KetSpace
    complex(kind(1d0))    , intent(in) :: A(:)
    integer               , intent(in) :: ik
    complex(kind(1d0)), allocatable , intent(inout) :: blk(:)
    integer :: jmi, jma, nj
    jmi = KetSpace.CCPWC_idim( ik - 1 ) + 1
    jma = KetSpace.CCPWC_idim( ik )
    nj  = jma - jmi + 1
    if(allocated(blk))deallocate(blk)
    allocate(blk(nj))
    blk=A(jmi:jma)
  end subroutine ExtractCCOSymvec

  subroutine ExtractCCOSymBlock( BraSpace, KetSpace, A, ib, ik, blk )
    type(ClassCCOSymSpace), intent(in) :: BraSpace
    type(ClassCCOSymSpace), intent(in) :: KetSpace
    complex(kind(1d0))    , intent(in) :: A(:,:)
    integer               , intent(in) :: ib, ik
    complex(kind(1d0)), allocatable , intent(inout) :: blk(:,:)
    integer :: imi, ima, ni, jmi, jma, nj
    imi = BraSpace.CCPWC_idim( ib - 1 ) + 1
    ima = BraSpace.CCPWC_idim( ib )
    ni  = ima - imi + 1
    jmi = KetSpace.CCPWC_idim( ik - 1 ) + 1
    jma = KetSpace.CCPWC_idim( ik )
    nj  = jma - jmi + 1
    if(allocated(blk))deallocate(blk)
    allocate(blk(ni,nj))
    blk=A(imi:ima,jmi:jma)
  end subroutine ExtractCCOSymBlock

  subroutine SetCCOSymBlock( BraSpace, KetSpace, A, ib, ik, blk )
    type(ClassCCOSymSpace), intent(in) :: BraSpace
    type(ClassCCOSymSpace), intent(in) :: KetSpace
    complex(kind(1d0))    , intent(inout) :: A(:,:)
    integer               , intent(in) :: ib, ik
    complex(kind(1d0)), allocatable , intent(inout) :: blk(:,:)
    integer :: imi, ima, ni, jmi, jma, nj, ibk, jbk
    !
    imi = BraSpace.CCPWC_idim( ib - 1 ) + 1
    ima = BraSpace.CCPWC_idim( ib )
    ni  = ima - imi + 1
    !
    jmi = KetSpace.CCPWC_idim( ik - 1 ) + 1
    jma = KetSpace.CCPWC_idim( ik )
    nj  = jma - jmi + 1
    !
    ibk=size(blk,1)
    jbk=size(blk,2)
    !
    if( ibk /= ni .or. jbk /= nj )call Assert("Size mismatch in SetCCOSymBlock")
    !
    A(imi:ima,jmi:jma)=blk
    !
  end subroutine SetCCOSymBlock

  subroutine ClassCCOSymSpaceGetBlock( self, A, ib, ik, blk )
    class(ClassCCOSymSpace), intent(in) :: self
    complex(kind(1d0))     , intent(in) :: A(:,:)
    integer                , intent(in) :: ib, ik
    complex(kind(1d0)), allocatable , intent(inout) :: blk(:,:)
    integer :: imi, ima, ni, jmi, jma, nj
    imi = self.CCPWC_idim( ib - 1 ) + 1
    ima = self.CCPWC_idim( ib )
    ni  = ima - imi + 1
    jmi = self.CCPWC_idim( ik - 1 ) + 1
    jma = self.CCPWC_idim( ik )
    nj  = jma - jmi + 1
    if(allocated(blk))deallocate(blk)
    allocate(blk(ni,nj))
    blk=A(imi:ima,jmi:jma)
  end subroutine ClassCCOSymSpaceGetBlock


  subroutine ClassCCOSymSpaceGetSideBlock( self, A, dim, ich, blk )
    class(ClassCCOSymSpace), intent(in) :: self
    complex(kind(1d0))     , intent(in) :: A(:,:)
    integer                , intent(in) :: ich, dim
    complex(kind(1d0)), allocatable , intent(inout) :: blk(:,:)
    integer :: ni, jmi, jma, nj
    !
    jmi = self.CCPWC_idim( ich - 1 ) + 1
    jma = self.CCPWC_idim( ich )
    nj  = jma - jmi + 1
    !
    select case( dim )
    case( 1 ) ! BRA CCO PWC STRIPE
       ni  = size(A,2)
       call realloc( blk, nj, ni )
       blk = A(jmi:jma,:)
    case( 2 ) ! KET CCO PWC STRIPE
       ni  = size( A, 1 )
       call realloc( blk, ni, nj )
       blk = A(:,jmi:jma)
    end select
    !
  end subroutine ClassCCOSymSpaceGetSideBlock


  subroutine ClassCCOSymSpaceSetBlock( self, A, ib, ik, blk )
    class(ClassCCOSymSpace), intent(in)    :: self
    complex(kind(1d0))     , intent(inout) :: A(:,:)
    integer                , intent(in)    :: ib, ik
    complex(kind(1d0))     , intent(in)    :: blk(:,:)
    integer :: imi, ima, ni, jmi, jma, nj, ibl, jbl
    imi = self.CCPWC_idim( ib - 1 ) + 1
    ima = self.CCPWC_idim( ib )
    ni  = ima - imi + 1
    jmi = self.CCPWC_idim( ik - 1 ) + 1
    jma = self.CCPWC_idim( ik )
    nj  = jma - jmi + 1
    ibl = size(blk,1)
    jbl = size(blk,2)
    if( ibl /= ni .or. jbl /= nj )then
       write(ERROR_UNIT,*)"Size mismatch in Set Block ",ib,ik
       write(ERROR_UNIT,*)"Loaded block sizes:",ibl, jbl
       write(ERROR_UNIT,*)"Expected block sizes:",ni,nj
       stop
    endif
    A(imi:ima,jmi:jma)=blk
  end subroutine ClassCCOSymSpaceSetBlock

  real(kind(1d0)) function ClassCCOSymSpaceGetChIonEnergy( self, ich ) result( Energy )
    class(ClassCCOSymSpace), intent(in) :: self
    integer                , intent(in) :: ich
    if(ich<=0)then
       Energy=0.d0
    else
       Energy = self.CCPWC_pionEval(ich) + self.CCPWC_ESHFT(ich) 
    endif
  end function ClassCCOSymSpaceGetChIonEnergy
  
  subroutine ClassCCOSymSpaceGetOpenChannels( self, Energy, nOpen, listOpen )
    class(ClassCCOSymSpace), intent(in)   :: self
    real(kind(1d0))        , intent(in)   :: Energy
    integer                , intent(out)  :: nOpen
    integer, allocatable   , intent(inout):: listOpen(:)

    real(kind(1d0)), allocatable :: ThresholdList(:)
    integer        , allocatable :: nThresholdChannels(:), ThresholdChannelList(:,:)
    integer                      :: nThresholds, iOpen, ich, iThr 

    call Self.GetThresholdSequence( &
         nThresholds, &
         ThresholdList, &
         nThresholdChannels, &
         ThresholdChannelList )
    !
    !.. Determines the number of open channels
    !..
    nOpen=0
    do iThr=1,nThresholds
       if(ThresholdList(iThr) >= Energy )exit
       do ich = 1, nThresholdChannels(iThr)
          nOpen=nOpen+1
       enddo
    enddo
    call realloc(listOpen,nOpen)
    !
    !.. List the open channels in order of threshold 
    !   energy following the intra-threshold order provided
    !   by GetThresholdSequence
    !..
    iOpen=0
    do iThr=1,nThresholds
       if(ThresholdList(iThr) >= Energy )exit
       do ich = 1, nThresholdChannels(iThr)
          iOpen=iOpen+1
          listOpen(iOpen) = ThresholdChannelList(ich,iThr)
       enddo
    enddo

    if(allocated( ThresholdList        )) deallocate( ThresholdList )
    if(allocated( nThresholdChannels   )) deallocate( nThresholdChannels )
    if(allocated( ThresholdChannelList )) deallocate( ThresholdChannelList )

  end subroutine ClassCCOSymSpaceGetOpenChannels
  
!!$  subroutine ClassCCOSymSpaceGetOpenChannels( self, Energy, nOpen, listOpen )
!!$    class(ClassCCOSymSpace), intent(in) :: self
!!$    real(kind(1d0))        , intent(in) :: Energy
!!$    integer                , intent(out):: nOpen
!!$    integer, allocatable   , intent(inout):: listOpen(:)
!!$    integer :: ich,iOpen
!!$    nOpen=0
!!$    do ich=1,self.GetNCh()
!!$       if(self.GetChIonEnergy(ich)<=Energy)nOpen=nOpen+1
!!$    enddo
!!$    call realloc(listOpen,nOpen)
!!$    iOpen=0
!!$    do ich=1,self.GetNCh()
!!$       if(self.GetChIonEnergy(ich)<=Energy)then
!!$          iOpen=iOpen+1
!!$          listOpen(iOpen)=ich
!!$       endif
!!$    enddo
!!$  end subroutine ClassCCOSymSpaceGetOpenChannels
  
  subroutine ClassCCOSymSpaceGetThresholdSequence( self, nThr, ThrList, nThrCh, ThrChList )
    class(ClassCCOSymSpace)     , intent(in)  :: self
    integer                     , intent(out) :: nThr
    real(kind(1d0)), allocatable, intent(out) :: ThrList(:)
    integer        , allocatable, intent(out) :: nThrCh(:)
    integer        , allocatable, intent(out) :: ThrChList(:,:)
    !
    integer, parameter   :: SORT_IN_INCREASING_ORDER=2
    integer                      :: nchan,ich,INFO, i1,i2,ndeg,ndegmax
    real(kind(1d0)), allocatable :: Thresholds(:)
    integer        , allocatable :: Channels(:)
    integer        , allocatable :: lwvec(:), lwvec2(:), iperm(:), Channels2(:)
    !
    nchan = self.GetNCh()
    call realloc( Channels  , nchan )
    call realloc( Thresholds, nchan )
    call realloc( lwvec2    , nchan )
    do ich=1,nchan
       Channels  ( ich ) = ich
       Thresholds( ich ) = self.GetChIonEnergy( ich )
       lwvec2    ( ich ) = self.GetChOrbAngMom( ich )
    enddo
    call dpsort( Thresholds, nchan, Channels, SORT_IN_INCREASING_ORDER, INFO )
    call realloc( lwvec     , nchan )
    do ich=1,nchan
       lwvec(ich)=lwvec2(Channels(ich))
    enddo
    deallocate(lwvec2)
    allocate(iperm,source=[(ich,ich=1,nchan)])
    allocate( Channels2, source=Channels )
    i1=1
    nThr=0
    ndegmax=0
    outer: do 
       nThr=nThr+1
       if( i1 == nchan )exit
       i2=i1+1
       inner: do 
          if( Thresholds(i2) > Thresholds(i1) )then
             i2=i2-1
             exit inner
          endif
          if(i2==nchan)exit inner
          i2=i2+1
       enddo inner
       ndeg=i2-i1+1
       ndegmax=max(ndeg,ndegmax)
       if( ndeg > 1 )then
          iperm(1:ndeg)=[(ich,ich=1,ndeg)]
          call ipsort( lwvec(i1:i2), ndeg, iperm(1:ndeg), SORT_IN_INCREASING_ORDER, INFO )
          do ich=i1,i2
             Channels(ich)=Channels2(i1-1+iperm(ich-i1+1))
          enddo
       endif
       if( i2 == nchan )exit outer
       i1=i2+1
    enddo outer
    deallocate(lwvec,iperm,Channels2)
    allocate(ThrList(nThr))
    allocate(nThrCh(nThr))
    allocate(ThrChList(ndegmax,nThr))
    ThrChList=0
    i1=1
    nThr=0
    do 
       nThr=nThr+1
       ThrList(nThr)=Thresholds(i1)
       if( i1 == nchan )then
          nThrCh(nThr)=1
          ThrChList(1,nThr)=Channels(i1)
          exit
       endif
       i2=i1+1
       do 
          if( Thresholds(i2) > Thresholds(i1) )then
             i2=i2-1
             exit 
          endif
          if(i2==nchan)exit 
          i2=i2+1
       enddo
       nThrCh( nThr ) = i2 - i1 + 1
       ThrChList(1:i2-i1+1,nThr)=Channels(i1:i2)
       if( i2 == nchan )exit 
       i1=i2+1
    enddo
    deallocate(Thresholds,Channels)
  end subroutine ClassCCOSymSpaceGetThresholdSequence
  
  real(kind(1d0)) function ClassCCOSymSpaceGetIonCharge( self ) result( Charge )
    class(ClassCCOSymSpace), intent(in) :: self
    Charge = self.totalCharge + 1.d0
  end function ClassCCOSymSpaceGetIonCharge
  
  function ClassCCOSymSpaceGetDir( self ) result(strn)
    class(ClassCCOSymSpace), intent(in) :: self
    character(len=:), allocatable :: strn
    allocate(strn,source=self.StoreDir//AUG_SUBDIR//self.term//"/")
  end function ClassCCOSymSpaceGetDir

  function GetCCOPairDir( BraSpace, KetSpace ) result(strn)
    type(ClassCCOSymSpace), intent(in) :: BraSpace, KetSpace
    character(len=:), allocatable :: strn
    allocate(strn,source=BraSpace.StoreDir//AUG_SUBDIR//BraSpace.term//"_"//KetSpace.term//"/")
  end function GetCCOPairDir

  function FormatPairCoreBlockName( BraSpace, KetSpace, i, j ) result(coreBlockName)
    type(ClassCCOSymSpace), intent(in) :: BraSpace, KetSpace
    integer               , intent(in) :: i,j
    character(len=:), allocatable :: coreBlockName
    allocate(coreBlockName, source = &
         BraSpace.GetTerm() // "_" // BraSpace.GetChName( i )//"_"//   &
         KetSpace.GetTerm() // "_" // KetSpace.GetChName( j ) )
  end function FormatPairCoreBlockName

  function ClassCCOSymSpaceCoreBlockName( self, i, j ) result(coreBlockName)
    class(ClassCCOSymSpace), intent(in) :: self
    integer                , intent(in) :: i,j
    character(len=:), allocatable :: coreBlockName
    allocate(coreBlockName, source = &
         self.GetTerm()//"_"//self.GetChName( i )//"_"//self.GetChName( j ) )
  end function ClassCCOSymSpaceCoreBlockName


  !> Save a matrix between CCO spaces in separate blocks
  subroutine SymCCO_SaveMatrixInBlocks( BraSpace, KetSpace, A, subdir_, prefix, suffix, form )
    type(ClassCCOSymSpace), intent(in) :: BraSpace
    type(ClassCCOSymSpace), intent(in) :: KetSpace
    complex(kind(1d0))    , intent(in) :: A(:,:)
    character(len=*)      , intent(in) :: subdir_
    character(len=*)      , intent(in) :: prefix
    character(len=*)      , intent(in) :: suffix
    character(len=*)      , intent(in) :: form
    !
    character(len=:)  , allocatable :: FileName
    character(len=:)  , allocatable :: CoreBlockName
    character(len=:)  , allocatable :: pairdir, subdir
    complex(kind(1d0)), allocatable :: blk(:,:)
    !
    integer :: ibch, ikch
    !
    pairdir = GetCCOPairDir( BraSpace, KetSpace )
    subdir = FormatAsDir( subdir_ )
    call system("mkdir -p "//pairdir)
    call system("mkdir -p "//pairdir//subdir)
    !
    do ibch = 0, BraSpace.GetNCh()
       do ikch = 0, KetSpace.GetNCh()
          !
          call ExtractCCOSymBlock( BraSpace, KetSpace, A, ibch, ikch, blk )
          !
          CoreBlockName = FormatPairCoreBlockName( BraSpace, KetSpace, ibch, ikch )
          allocate( FileName, source = pairdir // subdir // prefix // CoreBlockName // suffix )
          call SaveMatrix( FileName, blk, form )
          deallocate( FileName )
          !
       enddo
    enddo
    !
  end subroutine SymCCO_SaveMatrixInBlocks

  !> Load a matrix between CCO spaces in separate blocks
  subroutine SymCCO_LoadMatrixInBlocks( BraSpace, KetSpace, A, subdir_, prefix, suffix, form )
    type(ClassCCOSymSpace), intent(in) :: BraSpace
    type(ClassCCOSymSpace), intent(in) :: KetSpace
    complex(kind(1d0)), allocatable, intent(out) :: A(:,:)
    character(len=*)      , intent(in) :: subdir_
    character(len=*)      , intent(in) :: prefix
    character(len=*)      , intent(in) :: suffix
    character(len=*)      , intent(in) :: form
    !
    character(len=:)  , allocatable :: FileName
    character(len=:)  , allocatable :: CoreBlockName
    character(len=:)  , allocatable :: pairdir, subdir
    complex(kind(1d0)), allocatable :: blk(:,:)
    !
    integer :: ibch, ikch, ni, nj
    !
    pairdir = GetCCOPairDir( BraSpace, KetSpace )
    subdir = FormatAsDir( subdir_ )
    !
    ni = BraSpace.Get_TotalSize()
    nj = KetSpace.Get_TotalSize()
    if(allocated(A))deallocate(A)
    allocate(A(ni,nj))
    A=(0.d0,0.d0)
    do ibch = 0, BraSpace.GetNCh()
       do ikch = 0, KetSpace.GetNCh()
          !
          CoreBlockName = FormatPairCoreBlockName( BraSpace, KetSpace, ibch, ikch )
          allocate( FileName, source = pairdir // subdir // prefix // CoreBlockName // suffix )
          call LoadMatrix( FileName, blk, form )
          deallocate( FileName )
          
          call SetCCOSymBlock( BraSpace, KetSpace, A, ibch, ikch, blk )
       enddo
    enddo
    !
  end subroutine SymCCO_LoadMatrixInBlocks


  !> Save a matrix between CCO spaces in separate blocks
  subroutine ClassCCOSymSpaceSaveMatrixInBlocks( self, A, subdir_, prefix, suffix, form )
    class(ClassCCOSymSpace), intent(in) :: self
    complex(kind(1d0))    , intent(in) :: A(:,:)
    character(len=*)      , intent(in) :: subdir_
    character(len=*)      , intent(in) :: prefix
    character(len=*)      , intent(in) :: suffix
    character(len=*)      , intent(in) :: form
    !
    character(len=:)  , allocatable :: FileName
    character(len=:)  , allocatable :: CoreBlockName
    character(len=:)  , allocatable :: symdir, subdir
    complex(kind(1d0)), allocatable :: blk(:,:)
    !
    integer :: ibch, ikch
    !
    symdir = self.GetDir()
    subdir = FormatAsDir( subdir_ )
    call system("mkdir -p "//symdir)
    call system("mkdir -p "//symdir//subdir)
    !
    do ibch = 0, self.GetNCh()
       do ikch = 0, self.GetNCh()
          !
          call self.GetBlock( A, ibch, ikch, blk )
          CoreBlockName = self.CoreBlockName( ibch, ikch )
          allocate( FileName, source = symdir // subdir // prefix // CoreBlockName // suffix )
          call SaveMatrix( FileName, blk, form )
          deallocate( FileName )
          !
       enddo
    enddo
    !
  end subroutine ClassCCOSymSpaceSaveMatrixInBlocks

  !> Save a matrix between CCO spaces in separate blocks
  subroutine ClassCCOSymSpaceSaveMatrixInStripeBlocks( self, dim, dir, prefix, suffix, A, form )
    class(ClassCCOSymSpace), intent(in) :: self
    integer                , intent(in) :: dim
    character(len=*)       , intent(in) :: dir
    character(len=*)       , intent(in) :: prefix
    character(len=*)       , intent(in) :: suffix
    complex(kind(1d0))     , intent(in) :: A(:,:)
    character(len=*)       , intent(in) :: form
    !
    character(len=:)  , allocatable :: FileName
    character(len=:)  , allocatable :: ChName
    character(len=:)  , allocatable :: symdir
    complex(kind(1d0)), allocatable :: blk(:,:)
    !
    integer :: ich
    !
    symdir = self.GetDir()
    call system("mkdir -p "//symdir)
    call system("mkdir -p "//dir)
    !
    do ich = 0, self.GetNCh()
       !
       call self.GetSideBlock( A, dim, ich, blk )
       ChName = self.GetChName( ich )
       allocate( FileName, source = dir // prefix // ChName // suffix )
       call SaveMatrix( FileName, blk, form )
       deallocate( FileName )
       !
    enddo
    deallocate( blk )
    !
  end subroutine ClassCCOSymSpaceSaveMatrixInStripeBlocks


  !> Load a matrix within a CCO space from separate blocks
  subroutine ClassCCOSymSpaceLoadMatrixInBlocks( self, A, subdir_, prefix, suffix, form )
    class(ClassCCOSymSpace)        , intent(inout) :: self
    complex(kind(1d0)), allocatable, intent(inout) :: A(:,:)
    character(len=*)               , intent(in)    :: subdir_
    character(len=*)               , intent(in)    :: prefix
    character(len=*)               , intent(in)    :: suffix
    character(len=*)               , intent(in)    :: form
    !
    character(len=:)  , allocatable :: FileName
    character(len=:)  , allocatable :: CoreBlockName
    character(len=:)  , allocatable :: symdir, subdir
    complex(kind(1d0)), allocatable :: blk(:,:)
    !
    integer :: ibch, ikch, lda
    !
    symdir = self.GetDir()
    subdir = FormatAsDir( subdir_ )
    !
    if(allocated(A))deallocate(A)
    lda=self.Get_TotalSize()
    allocate(A(lda,lda))
    A=(0.d0,0.d0)
    do ibch = 0, self.GetNCh()
       do ikch = 0, self.GetNCh()
          !
          CoreBlockName = self.CoreBlockName( ibch, ikch )
          allocate( FileName, source = symdir // subdir // prefix // CoreBlockName // suffix )
          call LoadMatrix( FileName, blk, form )
          deallocate( FileName )
          call self.SetBlock( A, ibch, ikch, blk )
          !
       enddo
    enddo
    !
  end subroutine ClassCCOSymSpaceLoadMatrixInBlocks


  !> Save a matrix between CCO spaces in separate blocks
  subroutine ClassCCOSymSpaceSaveAllChProp( self )
    class(ClassCCOSymSpace), intent(in) :: self
    integer :: ich
    do ich = 0, self.GetNCh()
       call self.SaveChProp( ich )
    enddo
  end subroutine ClassCCOSymSpaceSaveAllChProp


  !> Save a matrix between CCO spaces in separate blocks
  subroutine ClassCCOSymSpaceSaveChProp( self, ich )
    class(ClassCCOSymSpace), intent(in) :: self
    integer                , intent(in) :: ich
    !
    character(len=:), allocatable :: symdir, FileName, Term_pion
    integer                       :: uid, ipion, lw, dim, icsf, ncsf_pion
    real(kind(1d0))               :: Ener_Pion
    real(kind(1d0)),  allocatable :: civec(:)
    
    symdir = self.GetDir()
    call system("mkdir -p "//symdir)
    call system("mkdir -p "//symdir//CCO_SUBDIR)
    !
    allocate( FileName, source = symdir // CCO_SUBDIR // self.CCPWC_Name(ich) )
    call OpenFile( FileName, uid, "write", "formatted" )
    if(ich==0)then
       write(uid,"(i0)") self.CCPWC_dim(0)
    else
       iPion      = self.CCPWC_iCIPion(   ich )
       Term_Pion = self.CCPWC_Term_Pion( ich )
       lw         = self.CCPWC_lw(ich)
       dim        = self.CCPWC_dim(ich)
       Ener_Pion  = ParentIons.Get_CIEnergy(iPion,Term_Pion)
       call ParentIons.Get_CIVector( iPion, Term_Pion, nCSF_pion, civec )
       write(uid,"(i0,x,a,x,i0,x,d24.16,x,i0)") iPion, Term_Pion, lw, Ener_Pion, dim
       write(uid,"(i0)") ncsf_pion
       write(uid,"(*(x,i0))") ( self.CCPWC_listCSFPWC( icsf, ich ), icsf = 1, ncsf_pion )
       write(uid,"(*(x,d24.16))") ( civec( icsf ) , icsf = 1, ncsf_pion )
       deallocate(civec)
    endif
    close(uid)
    deallocate( FileName )
    
  end subroutine ClassCCOSymSpaceSaveChProp
  

  !> Save a matrix between CCO spaces in separate blocks
  subroutine ClassCCOSymSpaceLoadAllChProp( self )
    class(ClassCCOSymSpace), intent(inout) :: self
    integer :: ich
    do ich = 0, self.GetNCh()
       call self.LoadChProp( ich )
    enddo
  end subroutine ClassCCOSymSpaceLoadAllChProp


  !> 
  subroutine ClassCCOSymSpaceLoadChProp( self, ich )
    class(ClassCCOSymSpace), intent(inout) :: self
    integer                , intent(in)    :: ich
    !
    character(len=:), allocatable :: symdir, FileName
    character(len=3)              :: Term_pion
    integer                       :: uid, ipion, lw, dim, icsf, ncsf_pion
    real(kind(1d0))               :: Ener_Pion
    !
    symdir = self.GetDir()
    !
    allocate( FileName, source = symdir // CCO_SUBDIR // self.CCPWC_Name(ich) )
    
    call OpenFile( FileName, uid, "read", "formatted" )
    if(ich==0)then
       read(uid,*) self.CCPWC_dim(0)
    else
       read(uid,*) iPion, Term_Pion, lw, Ener_Pion, dim
       self.CCPWC_iCIPion(   ich ) = iPion 
       self.CCPWC_Term_Pion( ich ) = Term_Pion
       self.CCPWC_lw(ich)          = lw
       self.CCPWC_dim(ich)         = dim
       self.CCPWC_pionEval(ich )   = Ener_Pion
       read(uid,*) ncsf_pion
       read(uid,*) (self.CCPWC_listCSFPWC( icsf, ich ), icsf = 1, ncsf_pion )
    endif
    close(uid)
    deallocate( FileName )
    !
  end subroutine ClassCCOSymSpaceLoadChProp
  

  subroutine ClassCCOSymSpaceEvalWeights( Self, zvec, nlc, wlc, npwc, wpwc )
    !
    class(ClassCCOSymSpace)     , intent(in)    :: Self
    complex(kind(1d0))          , intent(in)    :: zvec(:)
    integer                     , intent(out)   :: nlc, npwc
    real(kind(1d0)), allocatable, intent(inout) :: wlc(:), wpwc(:)
    integer :: ilc, ipwc, imi, ima

    nlc  = Self.CCPWC_dim(0)
    call realloc( wlc, nlc )
    do ilc=1,nlc
       wlc(ilc)=abs(zvec(ilc)**2)
    enddo

    npwc = Self.nCCPWC
    call realloc( wpwc, npwc )
    do ipwc=1,npwc
       imi = self.ccpwc_idim(ipwc-1) + 1
       ima = self.ccpwc_idim(ipwc)
       wpwc(ipwc)=abs(sum(zvec(imi:ima)**2)) ! <-- notice it is not the L2 norm
    enddo

  end subroutine ClassCCOSymSpaceEvalWeights

  subroutine ClassCCOSymSpaceTransform_BraKet_CSF_to_CCO( &
       self    , &
       zMat_CSF_CSF, &
       zMat_CCO_CCO  )
    !
    class(ClassCCOSymSpace)        , intent(in)    :: self
    complex(kind(1d0))             , intent(in)    :: zMat_CSF_CSF(:,:)
    complex(kind(1d0)), allocatable, intent(inout) :: zMat_CCO_CCO(:,:)
    complex(kind(1d0)), allocatable :: zMat_CSF_CCO(:,:)
    call self.Transform_Ket_CSF_to_CCO( zMat_CSF_CSF, zMat_CSF_CCO )
    call self.Transform_Bra_CSF_to_CCO( zMat_CSF_CCO, zMat_CCO_CCO )
    deallocate( zMat_CSF_CCO )
    !
  end subroutine ClassCCOSymSpaceTransform_BraKet_CSF_to_CCO

  subroutine Transform_SymBraKet_CSF_to_CCO( &
       BraSpace    , &
       KetSpace    , &
       zMat_CSF_CSF, &
       zMat_CCO_CCO  )
    !
    type(ClassCCOSymSpace)         , intent(in)    :: BraSpace
    type(ClassCCOSymSpace)         , intent(in)    :: KetSpace
    complex(kind(1d0))             , intent(in)    :: zMat_CSF_CSF(:,:)
    complex(kind(1d0)), allocatable, intent(inout) :: zMat_CCO_CCO(:,:)
    complex(kind(1d0)), allocatable :: zMat_CSF_CCO(:,:)
    call KetSpace.Transform_Ket_CSF_to_CCO( zMat_CSF_CSF, zMat_CSF_CCO )
    call BraSpace.Transform_Bra_CSF_to_CCO( zMat_CSF_CCO, zMat_CCO_CCO )
    deallocate( zMat_CSF_CCO )
    !
  end subroutine Transform_SymBraKet_CSF_to_CCO

  subroutine ClassCCOSymSpaceTransform_Ket_CSF_to_CCO( Self, zMat_XXX_CSF, zMat_XXX_CCO )
    !
    class(ClassCCOSymSpace)        , intent(in)    :: Self
    complex(kind(1d0))             , intent(in)    :: zMat_XXX_CSF(:,:)
    complex(kind(1d0)), allocatable, intent(inout) :: zMat_XXX_CCO(:,:)
    !
    character(len=*), parameter :: PROCEDURE_NAME="ClassCCOSymSpaceTransform_Ket_CSF_to_CCO"
    !
    integer :: nRows, nCols
    character(len=3) :: Term_Pion
    integer :: iCIPion, iCCPWC, jccmi, jccma
    integer :: iCSFPion, iCSFPWC, jcsfmi, jcsfma
    integer                      :: ncsf, nCCO
    real(kind(1d0)), allocatable :: civec(:)

    nRows=size(zMat_XXX_CSF,1)
    nCols=size(zMat_XXX_CSF,2)
    if( nCols /= self.CSFCCSymSpace.Get_TotalSize() )then
       write(ERROR_UNIT,"(a)"   ) "column mismatch in "//PROCEDURE_NAME
       write(ERROR_UNIT,"(a,i0)") "nCols        = ",nCols
       write(ERROR_UNIT,"(a,i0)") "Ket CSF Size = ",self.CSFCCSymSpace.Get_TotalSize()
       stop
    endif

    if(allocated(zMat_XXX_CCO))deallocate(zMat_XXX_CCO)
    nCCO = self.Get_TotalSize()
    allocate(zMat_XXX_CCO(nRows,nCCO))
    zMat_XXX_CCO=(0.d0,0.d0)

    !.. Copies the *:LOC block
    zMat_XXX_CCO(:,1:self.CCPWC_dim(0)) = zMat_XXX_CSF(:,1:self.CCPWC_dim(0))

    do iCCPWC = 1, self.nCCPWC
       
       Term_Pion = self.CCPWC_Term_Pion( iCCPWC )
       iCIPion   = self.CCPWC_ICIPion(   iCCPWC )
       
       jccmi     = self.CCPWC_idim( iCCPWC - 1 ) + 1
       jccma     = self.CCPWC_idim( iCCPWC )

       call ParentIons.Get_CIVector( iCIPion, Term_Pion, ncsf, civec )
       
       do iCSFPion = 1, ncsf
          
          iCSFPWC = self.CCPWC_listCSFPWC( iCSFPion, iCCPWC )
          
          jcsfmi= self.CSFCCSymSpace.Get_PWCidim( iCSFPWC - 1 ) + 1
          jcsfma= self.CSFCCSymSpace.Get_PWCidim( iCSFPWC )
          
          zMat_XXX_CCO( :, jccmi : jccma ) = zMat_XXX_CCO( :, jccmi : jccma ) + &
               zMat_XXX_CSF( :, jcsfmi:jcsfma ) * civec( iCSFPion )
       enddo

    enddo

  end subroutine ClassCCOSymSpaceTransform_Ket_CSF_to_CCO


  subroutine ClassCCOSymSpaceTransform_Bra_CSF_to_CCO( Self, zMat_CSF_XXX, zMat_CCO_XXX )
    !
    class(ClassCCOSymSpace)        , intent(in)    :: Self
    complex(kind(1d0))             , intent(in)    :: zMat_CSF_XXX(:,:)
    complex(kind(1d0)), allocatable, intent(inout) :: zMat_CCO_XXX(:,:)
    !
    character(len=*), parameter :: PROCEDURE_NAME="ClassCCOSymSpaceTransform_Bra_CSF_to_CCO"
    integer :: nRows, nCols
    complex(kind(1d0)), allocatable :: zmat1(:,:), zmat2(:,:)

    nRows=size(zMat_CSF_XXX,1)
    nCols=size(zMat_CSF_XXX,2)
    if( nRows /= self.CSFCCSymSpace.Get_TotalSize() )then
       write(ERROR_UNIT,"(a)"   ) "row mismatch in "//PROCEDURE_NAME
       write(ERROR_UNIT,"(a,i0)") "nRows        = ",nRows
       write(ERROR_UNIT,"(a,i0)") "Bra CSF Size = ",self.CSFCCSymSpace.Get_TotalSize()
       stop
    endif

    if(allocated(zMat_CCO_XXX))deallocate(zMat_CCO_XXX)
    allocate(zmat1(nCols,nRows))
    zmat1=transpose(zMat_CSF_XXX)
    
    call self.Transform_Ket_CSF_to_CCO( zmat1, zmat2 )
    deallocate(zmat1)
    nCols=size(zmat2,1)
    nRows=size(zmat2,2)
    allocate(zMat_CCO_XXX(nRows,nCols))
    zMat_CCO_XXX=transpose(zmat2)
    deallocate(zmat2)
    
  end subroutine ClassCCOSymSpaceTransform_Bra_CSF_to_CCO
    
  integer function ClassCCOSymSpaceGetNch( Self ) result(n)
    class(ClassCCOSymSpace) , intent(in) :: Self
    n = Self.NCCPWC
  end function ClassCCOSymSpaceGetNch

  subroutine ClassCCOSymSpaceGetChIndexes( Self, ich, i1, i2 ) 
    class(ClassCCOSymSpace), intent(in) :: Self
    integer                , intent(in) :: ich
    integer                , intent(out):: i1, i2
    i1 = Self.CCPWC_idim( ich - 1 ) + 1 
    i2 = Self.CCPWC_idim( ich )
  end subroutine ClassCCOSymSpaceGetChIndexes

  real(kind(1d0)) function ClassCCOSymSpaceGetCAPSF( Self, ich ) result( CAPSF )
    class(ClassCCOSymSpace) , intent(in) :: Self
    integer                 , intent(in) :: ich
    CAPSF=0.d0
    if(ich<1.or.ich>Self.NCCPWC)return
    CAPSF=self.CCPWC_CAPSF(ich)
  end function ClassCCOSymSpaceGetCAPSF

  real(kind(1d0)) function ClassCCOSymSpaceGetESHFT( Self, ich ) result( ESHFT )
    class(ClassCCOSymSpace) , intent(in) :: Self
    integer                 , intent(in) :: ich
    ESHFT=0.d0
    if(ich<1.or.ich>Self.NCCPWC)return
    ESHFT=self.CCPWC_ESHFT(ich)
  end function ClassCCOSymSpaceGetESHFT

  integer function ClassCCOSymSpaceGet_TotalSize( Self ) result(n)
    class(ClassCCOSymSpace) , intent(in) :: Self
    n = Self.CCPWC_totdim
  end function ClassCCOSymSpaceGet_TotalSize

  function ClassCCOSymSpaceGetTerm( Self ) result( term )
    class(ClassCCOSymSpace) , intent(in) :: Self
    character(len=:), allocatable :: term
    allocate(term,source=trim(adjustl(Self.term)))
  end function ClassCCOSymSpaceGetTerm

  integer function ClassCCOSymSpaceGetL( Self ) result( L )
    class(ClassCCOSymSpace) , intent(in) :: Self
    L=index(lstrnUC,Self.term(2:2))-1
  end function ClassCCOSymSpaceGetL

  integer function ClassCCOSymSpaceGetChOrbAngMom( Self, iCh ) result( lw )
    class(ClassCCOSymSpace) , intent(in) :: Self
    integer                 , intent(in) :: iCh
    lw=Self.CCPWC_lw(iCh)
  end function ClassCCOSymSpaceGetChOrbAngMom

  function ClassCCOSymSpaceGetChName( Self, iCh ) result( ChName )
    class(ClassCCOSymSpace) , intent(in) :: Self
    integer                 , intent(in) :: iCh
    character(len=:), allocatable :: ChName
    allocate(ChName,source=trim(adjustl(Self.CCPWC_Name(iCh))))
  end function ClassCCOSymSpaceGetChName

  subroutine ClassCCOSymSpaceSetCSFSymSpace( Self, CSFCCSymSpace )
    class(ClassCCOSymSpace) , intent(inout) :: Self
    type(ClassCSFCCSymSpace), target, intent(inout) :: CSFCCSymSpace
    Self.CSFCCSymSpace => CSFCCSymSpace
  end subroutine ClassCCOSymSpaceSetCSFSymSpace

  subroutine ClassCCOSymSpaceInit( Self, store, term, CSFCCSymSpace )

    class(ClassCCOSymSpace) , intent(inout) :: Self
    character(len=*)        , intent(in)    :: store
    character(len=*)        , intent(in)    :: term
    type(ClassCSFCCSymSpace), target, intent(inout) :: CSFCCSymSpace

    integer :: iTerm_pion, iCIPion, iCCPWC, iCSFPion, iCSFPWC, ilw, iBuf
    integer :: nterms_pion !auxiliary
    integer :: iPion
    character(len=3) :: Term_Pion
    integer, allocatable :: nlw_per_pion_term(:)
    integer, allocatable :: lwvec_per_pion_term(:,:)

    Self.StoreDir = FormatAsDir( store )
    Self.Term = adjustl(term)

    Self.CSFCCSymSpace => CSFCCSymSpace

    nterms_pion = ParentIons.GetNterms()
    call Self.CSFCCSymSpace.Get_lw_ppit( nlw_per_pion_term, lwvec_per_pion_term )

    !.. Counts the total number of CC PWC channels
    !..
    self.nCCPWC = 0
    do iTerm_pion = 1, ParentIons.GetNterms()
       do iCIPion = 1, ParentIons.Get_nCIperTerm( iTerm_Pion ) 
          self.nCCPWC = self.nCCPWC + nlw_per_pion_term( iTerm_Pion ) 
       enddo
    enddo
    !
    ! allocate( self.Ccpwc_iterm_pion( self.nCCPWC ) )
    allocate( self.Ccpwc_term_pion(  self.nCCPWC ) )
    allocate( self.CCPWC_iCIPion(    self.nCCPWC ) )
    allocate( self.CCPWC_lw(         self.nCCPWC ) )
    allocate( self.CCPWC_dim(      0:self.nCCPWC ) )
    !
    iBuf=ParentIons.GetMaxNCSF()
    allocate( self.CCPWC_listCSFPWC( iBuf, self.nCCPWC ) )
    !
    iCCPWC = 0
    self.CCPWC_dim(0) = self.CSFCCSymSpace.GetnRef()
    !
    do iTerm_pion = 1, nterms_pion
       !
       do iCIPion = 1, ParentIons.Get_nCIperTerm( iTerm_Pion ) 
          !
          do ilw = 1, nlw_per_pion_term( iTerm_Pion )
             !
             iCCPWC = iCCPWC + 1
             !
             ! self.CCPWC_iTerm_Pion( iCCPWC ) = iTerm_pion
             self.CCPWC_Term_Pion(  iCCPWC ) = ParentIons.GetTermName( iTerm_pion )
             self.CCPWC_iCIPion(    iCCPWC ) = iCIPion
             self.CCPWC_lw(         iCCPWC ) = lwvec_per_pion_term( ilw, iTerm_pion )
             !
             !.. Cycle over the CSFPWC to ascertain whether they correspond 
             !   to the same term and photoelectron angular momentum.
             !..
             iCSFPion=0
             do iCSFPWC = 1, self.CSFCCSymSpace.Get_nPWC()
                !
                if( self.CSFCCSymSpace.Get_iTermOfPWCpion( iCSFPWC ) /= iTerm_pion              ) cycle
                if( self.CSFCCSymSpace.Get_lwOfPWC(        iCSFPWC ) /= self.CCPWC_lw( iCCPWC ) ) cycle
                !
                iCSFPion = iCSFPion + 1
                !
                self.CCPWC_listCSFPWC( iCSFPion, iCCPWC ) = iCSFPWC
                self.CCPWC_dim( iCCPWC ) = self.CSFCCSymSpace.Get_PWCdim( iCSFPWC )
                !
             enddo
             
          enddo
       enddo
    enddo
    !
    allocate( self.CCPWC_idim( -1:self.nCCPWC ) )
    self.CCPWC_idim=0
    do iCCPWC = 0, self.nCCPWC
       self.CCPWC_idim( iCCPWC ) = self.CCPWC_idim( iCCPWC-1 ) + self.CCPWC_dim( iCCPWC )
    enddo
    self.CCPWC_totdim = self.CCPWC_idim( self.nCCPWC )

    allocate( self.CCPWC_Name( 0:self.nCCPWC ) )
    self.CCPWC_Name(0)="LOC"
    
    do iCCPWC = 1, self.nCCPWC
       iPion     = self.CCPWC_iCIPion(   iCCPWC )
       Term_Pion = self.CCPWC_Term_Pion( iCCPWC )
       write(self.CCPWC_Name(iCCPWC),"(i0,a,i0)") iPion, &
            "."//Term_Pion, self.CCPWC_lw( iCCPWC )
    enddo

    allocate( self.CCPWC_CAPSF( 0:self.nCCPWC ) )
    self.CCPWC_CAPSF=0.d0

    allocate( self.CCPWC_ESHFT( 0:self.nCCPWC ) )
    self.CCPWC_ESHFT=0.d0

    deallocate( nlw_per_pion_term, lwvec_per_pion_term )

  end subroutine ClassCCOSymSpaceInit


  subroutine ClassCCOSymSpaceInitFromConfigFile( Self, store, ConfigFile, term, selection_ )

    class(ClassCCOSymSpace) , intent(inout) :: Self
    character(len=*)        , intent(in)    :: store
    character(len=*)        , intent(in)    :: ConfigFile
    character(len=*)        , intent(in)    :: term
    character(len=*), optional, intent(in)  :: selection_

    character(len=:), allocatable :: selection

    if( present(selection_) )then
       allocate(selection,source=trim(adjustl(selection_)))
    else
       allocate(selection,source=trim(adjustl(DEFAULT_SELECTION)))
    endif
    
    Self.StoreDir = FormatAsDir( store )

    call Self.ParseConfigFile( ConfigFile, term, selection )
    
  end subroutine ClassCCOSymSpaceInitFromConfigFile



  subroutine ClassCCOSymSpaceParseConfigFile( Self, &
       FileName, &
       term, selection )
    class(ClassCCOSymSpace), intent(inout) :: Self
    character(len=*)       , intent(in)    :: FileName
    character(len=*)       , intent(in)    :: term, selection
    
    character(len=:), allocatable :: FullText, strnBuf
    integer :: i0, i1, i2

    !.. Load the Configuration file
    !call GetFullText( FileName, FullText )
    call GetFullText( FileName, FullText, ASCII_NEW_LINE )

    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "NAME", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("NAME missing in "//trim(FileName))
    allocate( Self.name, source=trim(adjustl(strnBuf)) )

    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "TOTAL_CHARGE", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("TOTAL_CHARGE missing in "//trim(FileName))
    read(strnBuf,*) Self.TotalCharge

    !.. Finds the section for the current term
    i0=1
    do

       i1=index(FullText(i0:), "[")
       if(i1<=0)then
          call ErrorMessage(Self.term//" not found in "//trim(FileName))
          stop
       endif
       i1=i0-1+i1
       i2=index(FullText(i1+1:),"]")
       if(i2<=0)then
          call ErrorMessage("Syntax error in "//trim(FileName)//" (unbalanced [] parentheses)")
          stop
       endif
       i2=i1+i2


       !.. Find the first character after the { that delimits the current term domain
       if( trim(adjustl(term)) /= trim( adjustl( FullText(i1+1:i2-1)  ) ) )then
          i0=i2+1
          cycle
       endif
       i1=index( FullText(i2+1:), "{" )
       i1=i2+i1+1

       !.. Find the beginning of the next term domain
       i2=index( FullText(i1:)  , "[" )
       if(i2<=0)then
          !.. the present term is also the last one in the configuration file
          i2=len_trim(FullText)
       else
          !.. reset the index to its absolute value
          i2=i1-1+i2-1
       endif
       !.. Find the end of the current term domain
       i2=index( FullText(:i2), "}", back=.TRUE. )
       i2=i2-1

       Self.Term = adjustl(term)
       exit

    enddo

    call Self.ParseSection( FullText(i1:i2), selection )
    
  end subroutine ClassCCOSymSpaceParseConfigFile


  subroutine ClassCCOSymSpaceParseSection( Self, &
       Text, selection )
    class(ClassCCOSymSpace), intent(inout) :: Self
    character(len=*)       , intent(in)    :: Text
    character(len=*)       , intent(in)    :: selection
    
    integer :: i1, i2
    
    !.. Finds the beginning of the current selection domain
    i1=index(Text,"@"//selection)
    if(i1<=0)then
       call ErrorMessage("Selection "//selection//" missing: abort")
       stop
    endif
    i2=index(Text(i1+1:),"{")
    i1=i1+i2+1

    !.. Find the beginning of the next selection domain
    i2=index( Text(i1:)  , "@" )
    if(i2<=0)then
       !.. the present term is also the last one in the configuration file
       i2=len_trim(Text)
    else
       !.. reset the index to its absolute value
       i2=i1-1+i2-1
    endif
    !.. Find the end of the current term domain
    i2=index( Text(i1:i2), "}" )+i1-2

    allocate(Self.Selection,source=selection)

    call Self.ParseSelection(Text(i1:i2))
    
  end subroutine ClassCCOSymSpaceParseSection

  function GetIonTermfromCCOChName( chstrn ) result( term )
    character(len=*), intent(in)  :: chstrn
    character(len=:), allocatable :: term
    integer :: i1
    i1=index(chstrn,".")
    allocate(term,source=chstrn(i1+1:i1+3))
  end function GetIonTermfromCCOChName

  subroutine GetCCOChQNfromStrn( chstrn, mult, L, par, n, lw )
    character(len=*), intent(in)  :: chstrn
    integer         , intent(out) :: mult, L, par, n, lw
    integer :: i1
    character(len=16) :: strnBuf
    i1=index(chstrn,".")
    strnBuf=chstrn(1:i1-1)
    read(strnBuf,*)n
    strnBuf=chstrn(i1+1:i1+1)
    read(strnBuf,*)mult
    L=index(lstrnUC,chstrn(i1+2:i1+2))-1
    par=index(pstrn  ,chstrn(i1+3:i1+3))-1
    strnBuf=chstrn(i1+4:)
    read(strnBuf,*)lw
  end subroutine GetCCOChQNfromStrn

  integer function GetCCOChId( chstrn ) result( chid )
    character(len=*), intent(in) :: chstrn
    integer :: n, lw, mult, L, par
    character(len=3) :: strn3
    strn3=trim(adjustl(chstrn))
    if(strn3.is."LOC")then
       chid=0
       return
    else
       call GetCCOChQNfromStrn( chstrn, mult, L, par, n, lw )
       chid = &
            lw   * 10**  LW_LOG_OFFSET + &
            n    * 10**   N_LOG_OFFSET + &
            par  * 10** PAR_LOG_OFFSET + &
            L    * 10**   L_LOG_OFFSET + &
            mult * 10**MULT_LOG_OFFSET
    endif
  end function GetCCOChId



  subroutine ClassCCOSymSpaceParseSelection( Self, &
       Text )
    class(ClassCCOSymSpace), intent(inout) :: Self
    character(len=*)       , intent(in)    :: Text
    
    character(len=:)  , allocatable :: ch_list, token
    character(len=256), allocatable :: ch_namevec(:)
    character(len=256), allocatable :: ch_namevec2(:)
    integer           , allocatable :: ch_idvec(:), iperm(:)
    integer :: nch, ich, ier, npwc, ncsfmax, isep, ipwc
    integer, parameter :: INCREASING_ORDER_AND_PERMUTE = 2
    character(len=256) :: NameBuf, AttrBuf

    ch_list = ExpandChannelList( Text )

    !.. Tokenize the channel list
    nch = nTokens( ch_list, CCO_CH_SEPARATOR )
    allocate( ch_namevec( nch ) )
    allocate( ch_idvec( nch ) ) 
    do ich = 1, nch
       call GetToken( ch_list, ich, token, CCO_CH_SEPARATOR )
       ch_namevec( ich ) = adjustl(token)
       ch_idvec ( ich ) = GetCCOChId( token )
       write(*,"(x,i3,x,i,x,a)") ich, ch_idvec( ich ), trim(ch_namevec(ich))
    enddo

    !.. Order the channels by their id
    allocate( iperm, source = [(ich,ich=1,nch)] )
    call IPSORT( ch_idvec, nch, iperm, INCREASING_ORDER_AND_PERMUTE, ier )
    if(ier/=0)call Assert("Error in channel id sorting")

    allocate(ch_namevec2,source=ch_namevec)
    do ich=1,nch
       ch_namevec(ich)=ch_namevec2(iperm(ich))
    enddo
    deallocate(ch_namevec2)
    
    !.. Check for duplicates
    do ich=1,nch-1
       if(ch_idvec(ich)==ch_idvec(ich+1))then
          write(*,*) "Duplicate in CCO config file "//Self.term//" "//Self.selection
          stop
       endif
    enddo
    
    !*** MUST CHECK COMPATIBILITY OF QUANTUM NUMBERS
    ! do ich=1,nch
    !    if( .not. Self.CheckQN(ch_namevec(ich)) )then
    !       call ErrorMessage("Channel "//ch_namevec(ich)//" incompatible with "//Self.term)
    !       stop
    !    endif
    ! enddo

    !.. Now it can initialize the CCO sym space

    !*** At the moment, assumes that LOC is always present
    npwc=nch-1
    self.nCCPWC=npwc

    allocate(self.CCPWC_Name(    0:npwc ))
    allocate(self.CCPWC_term_pion( npwc ))
    allocate(self.CCPWC_iCIpion(   npwc ))
    allocate(self.CCPWC_pionEval(  npwc ))
    allocate(self.CCPWC_lw(        npwc ))
    allocate(self.CCPWC_dim(     0:npwc ))
    nCSFmax=ParentIons.GetMaxNCSF()
    allocate(self.CCPWC_listCSFPWC(ncsfmax,npwc ))
    self.CCPWC_listCSFPWC=0.d0
    !
    allocate(self.CCPWC_Attr( 0:npwc ))
    allocate(self.CCPWC_CAPSF(  npwc ))
    allocate(self.CCPWC_ESHFT(  npwc ))

    do ich = 1, nch
       NameBuf = ch_namevec(ich)
       AttrBuf=" "
       isep=index(NameBuf," ")
       if(isep>0)then
          AttrBuf=NameBuf(isep+1:)
          NameBuf(isep:)=" "
       endif
       self.CCPWC_Name( ich-1 ) = NameBuf
       self.CCPWC_Attr( ich-1 ) = adjustl(AttrBuf)
       write(*,"(i4,x,a)") ich-1, self.CCPWC_Name( ich-1 )//": "//trim(self.CCPWC_Attr( ich-1 ))
    enddo

    
    !.. Parse the attribute list of the LOC channel
    call CCLOC_ParseAttrLst( self.CCPWC_Attr(0), self.CCPWC_NMAX_LOC )
    
    !.. Parse the attribute list of the PWC channels
    do ipwc=1,npwc
       call CCPWC_ParseAttrLst(    &
            self.CCPWC_Attr(ipwc), &
            self.CCPWC_CAPSF(ipwc),& 
            self.CCPWC_ESHFT(ipwc) )
    enddo

    call self.LoadAllChProp()

    allocate(self.CCPWC_idim(   -1:npwc ))
    self.CCPWC_idim(-1)=0
    do ich=0,npwc
       self.CCPWC_idim(ich)=self.CCPWC_idim(ich-1)+self.CCPWC_dim(ich)
    enddo
    self.CCPWC_totdim=self.CCPWC_idim(npwc)

  end subroutine ClassCCOSymSpaceParseSelection


  function ExpandChannelList( Text ) result ( list )
    character(len=*), intent(in)  :: Text
    character(len=:), allocatable :: list
    character(len=10000) :: strn
    character(len=:), allocatable :: token
    integer :: iToken, icur
    strn = " "
    do iToken=1, nTokens(Text, CCO_CH_SEPARATOR )
       call GetToken(Text,iToken,token, CCO_CH_SEPARATOR )
       icur=len_trim(strn)
       if(icur==0)then
          icur=1
       else
          strn(icur+1:icur+1)=CCO_CH_SEPARATOR
          icur=icur+2
       endif
       strn(icur:)=ExpandChannelToken( token )
    enddo
    allocate(list,source=trim(adjustl(strn)))
  end function ExpandChannelList

  function ExpandChannelToken( token ) result ( list )
    character(len=*), intent(in)  :: token
    character(len=:), allocatable :: list
    character(len=10000) :: strn
    character(len=256)   :: chstrn
    integer :: i1, i2, nmi, nma, n, icur
    !
    i1=index(token,"(")
    if(i1<=0) then
       allocate(list,source=trim(adjustl(token)))
       return
    endif
    !
    i2=index(token,":")
    read(token(i1+1:i2-1),*) nmi
    i1=i2
    i2=index(token,")")
    read(token(i1+1:i2-1),*) nma
    !
    strn=" "
    do n = nmi, nma
       write(chstrn,"(i0,a)")n,"."//trim(adjustl(token(i2+1:)))
       icur=len_trim(strn)
       if(icur==0)then
          icur=1
       else
          strn(icur+1:icur+1)=CCO_CH_SEPARATOR
          icur=icur+2
       endif
       strn(icur:)=trim(adjustl(chstrn))
    enddo
    !
    allocate(list,source=trim(adjustl(strn)))
    !
  end function ExpandChannelToken

  subroutine ClassCCOSymSpaceShow( Self, unit )
    !> Class of the electronic space.
    class(ClassCCOSymSpace), intent(in) :: Self
    integer, optional     , intent(in) :: unit
    integer :: outunit
    integer :: iCCPWC, iCSFPion, iPion
    character(len=3) :: Term_Pion
    !
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit

    !.. List the properties of the CCPWC basis
    !..
    write(*,*) 
    write(*,*) "Properties of the CC Partial Wave Channels in "//self.Term//" symmetry"
    write(*,*) "======================================================================"
    write(*,"(a,i0)")" Number of CC PWC:  ", self.nCCPWC
    write(*,"(a,i0)")" Number of CC REF:  ", self.CCPWC_dim(0)
    write(*,"(a,i0)")" MAX NLOC        :  ", self.CCPWC_NMAX_LOC
    write(*,"(a,i0)")" Total dimension :  ", self.CCPWC_totdim
    write(*,"(a)") " CCPWC  name  CIion    Ion Energy     l   dim  idim   CAPSF   Ener. Shift"//&
         "     CSF PWC composing the CCO PWC"
    do iCCPWC = 1, self.nCCPWC
       iPion      = self.CCPWC_iCIPion(   iCCPWC )
       Term_Pion  = self.CCPWC_Term_Pion( iCCPWC )
       write(*,"(i4,2x,a,3x,i0,a,f12.6,a,x,i2,x,i5,x,i6,x,e8.2,x,f12.6,x,*(x,i3))") iCCPWC, &
            trim(self.CCPWC_Name( iCCPWC )),&
            iPion,"."//Term_Pion//&
            " [E=",ParentIons.Get_CIEnergy(iPion,Term_Pion),"]",&
            self.CCPWC_lw(   iCCPWC ), &
            self.CCPWC_dim(  iCCPWC ), &
            self.CCPWC_idim( iCCPWC ), &
            self.CCPWC_CAPSF(iCCPWC ), &
            self.CCPWC_ESHFT(iCCPWC ), &
            ( self.CCPWC_listCSFPWC( iCSFPion, iCCPWC ), iCSFPion = 1, ParentIons.Get_nCSFperTerm( Term_Pion ) )
    enddo

  end subroutine ClassCCOSymSpaceShow



  ! function ClassCCOSymSpaceIonTermDir( Self, iTerm ) result(strnBuf)
  !   class(ClassCCOSymSpace), intent(in)  :: Self
  !   integer               , intent(in)  :: iTerm
  !   character(len=:)      , allocatable :: strnBuf
  !   allocate(strnBuf,source=Self.StoreDir//ION_SUBDIR//Self.ionTerms(iTerm)//"/")
  ! end function ClassCCOSymSpaceIonTermDir

  ! function ClassCCOSymSpaceAugTermDir( Self, iTerm ) result(strnBuf)
  !   class(ClassCCOSymSpace), intent(in)  :: Self
  !   integer               , intent(in)  :: iTerm
  !   character(len=:)      , allocatable :: strnBuf
  !   allocate(strnBuf,source=Self.StoreDir//AUG_SUBDIR//Self.AugTerms(iTerm)//"/")
  ! end function ClassCCOSymSpaceAugTermDir

  ! subroutine ClassCCOSymSpaceSetupDirTree( Self )
  !   class(ClassCCOSymSpace), intent(in) :: Self
  !   integer :: iTerm
  !   call system("mkdir -p "//Self.StoreDir)
  !   call system("mkdir -p "//Self.StoreDir//AUG_SUBDIR)
  !   call system("mkdir -p "//Self.StoreDir//ION_SUBDIR)
  !   call system("mkdir -p "//Self.StoreDir//ONE_SUBDIR)
  !   do iTerm = 1, Self.nIonTerms
  !      call system("mkdir -p "//Self.IonTermDir( iTerm ) )
  !   enddo
  !   do iTerm = 1, Self.nAugTerms
  !      call system("mkdir -p "//Self.AugTermDir( iTerm ) )
  !   enddo
  ! end subroutine ClassCCOSymSpaceSetupDirTree


end module ModuleCCOSymSpace
