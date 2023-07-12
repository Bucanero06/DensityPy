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
module ModuleSymESpace

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling        
  use ModuleSystemUtils
  use ModuleString   
  use ModuleAngularMomentum   
  use ModuleConstants
  use ModuleIO       
  use ModuleMatrix   
  use ModuleGroups               
  use ModuleXlm                  
  use ModuleOrbitalBasis
  use ModuleParentIons           
  use ModuleCloseCouplingChannels
  use ModuleDensityMatrices
  use ModuleIntegrals

  implicit none
  private

  character(len=*), public, parameter :: OverlapLabel     ='S'
  character(len=*), public, parameter :: HamiltonianLabel ='H'
  character(len=*), public, parameter :: KinEnergyLabel   ='K'
  character(len=*), public, parameter :: R_LABEL          ='R'
  character(len=*), public, parameter :: X_LABEL          ='X'
  character(len=*), public, parameter :: Y_LABEL          ='Y'
  character(len=*), public, parameter :: Z_LABEL          ='Z'
  character(len=*), public, parameter :: DX_LABEL         ='DX'
  character(len=*), public, parameter :: DY_LABEL         ='DY'
  character(len=*), public, parameter :: DZ_LABEL         ='DZ'
  character(len=*), public, parameter :: CAPLabel         ='CAP'
  integer        ,  public, parameter :: N_SESSES_ID      = 4
  integer        ,  public, parameter :: SESSES_LOAD      = 1
  integer        ,  public, parameter :: SESSES_BUILD     = 2
  integer        ,  public, parameter :: SESSES_SAVE      = 3
  integer        ,  public, parameter :: SESSES_FREE      = 4

  logical        , private, parameter :: CCBLOCK_FORMATTED_WRITE = .TRUE. 
  
  !> Class of a general operator block.
  type, private :: ClassCCBlock
     private
     character(len=:)          , allocatable   :: Root
     character(len=:)          , allocatable   :: OpLabel
     character(len=:)          , allocatable   :: IDLabel
     type(ClassXlm)            , pointer       :: OpXlm 
     class(ClassCloseCouplingChannel), pointer :: IonChBra 
     class(ClassCloseCouplingChannel), pointer :: IonChKet 
     integer                                   :: nOrbBra
     integer                                   :: nOrbKet
     type(ClassMatrix)                         :: Block
     logical                                   :: initialized
     logical                                   :: FormattedWrite
   contains
     procedure, public :: Init                 => ClassCCBlockInit
     procedure, public :: Load                 => ClassCCBlockLoad
     procedure, public :: Build                => ClassCCBlockBuild
     procedure, public :: Save                 => ClassCCBlockSave
     procedure, public :: Scale                => ClassCCBlockScale
     procedure, public :: Free                 => ClassCCBlockFree
     procedure, public :: IsInitialized        => ClassCCBlockIsInitialized
     procedure, public :: GetLabel             => ClassCCBlockGetLabel
     procedure, public :: GetIonBraLabel       => ClassCCBlockGetIonBraLabel
     procedure, public :: GetIonKetLabel       => ClassCCBlockGetIonKetLabel
     procedure, public :: IsType               => ClassCCBlockIsType
     procedure, public :: GetStorageDir        => ClassCCBlockGetStorageDir
     procedure, public :: GetBlockName         => ClassCCBlockGetBlockName
     procedure, public :: GetSize              => ClassCCBlockGetSize
     procedure, public :: FetchBlock           => ClassCCBlockFetchBlock
     procedure, public :: GetNRows             => ClassCCBlockGetNRows
     procedure, public :: GetNColumns          => ClassCCBlockGetNColumns
     procedure, public :: Transpose            => ClassCCBlockTranspose
     procedure, public :: SetFormattedWrite    => ClassCCBlockSetFormattedWrite
     procedure, public :: UnsetFormattedWrite  => ClassCCBlockUnsetFormattedWrite
     procedure, public :: GetIntegralLabel1B   => ClassCCBlockGetIntegralLabel1B
     procedure, public :: GetIntegralLabel2B   => ClassCCBlockGetIntegralLabel2B
     generic  , public :: FetchMatrix          => ClassCCBlockFetchMatrix, &
          ClassCCBlockFetchMatrixCM
     generic  , public :: assignment(=) => ClassCCBlockCopy
     procedure, private:: ClassCCBlockCopy
     procedure, private:: ClassCCBlockFetchMatrix
     procedure, private:: ClassCCBlockFetchMatrixCM
     !
     final                 :: ClassCCBlockFinal
  end type ClassCCBlock


  !> Operators between sym spaces
  type, public :: ClassSESSESBlock
     private
     class(ClasssymESpace), pointer    :: SpaceBra
     class(ClasssymESpace), pointer    :: SpaceKet
     character(len=:)    , allocatable :: OpLabel
     character(len=:)    , allocatable :: IDLabel
     type(ClassCCBlock)  , allocatable :: mBlock(:,:)
     character(len=:)    , allocatable :: StorageDir
     logical                           :: INITIALIZED
     logical                           :: AvailableBlocks = .FALSE.
   contains
     procedure, public :: Init          => ClassSESSESBlock_Init
     procedure, public :: Driver        => ClassSESSESBlock_Driver
     procedure, public :: IsInitialized => ClassSESSESBlock_IsInitialized
     procedure, public :: Load          => ClassSESSESBlock_Load
     procedure, public :: Save          => ClassSESSESBlock_Save
     procedure, public :: Free          => ClassSESSESBlock_Free
     procedure, public :: GetSize       => ClassSESSESBlock_GetSize
     procedure, public :: Assemble      => ClassSESSESBlock_Assemble
     procedure, public :: GetBlock      => ClassSESSESBlock_GetBlock
     final             :: ClassSESSESBlock_Final
  end type ClassSESSESBlock


  !.. Symmetric Electronic Space
  type, public :: ClassSymESpace
     private
     integer                          :: Charge
     integer                          :: Multiplicity
     type(ClassGroup), pointer        :: Group
     type(ClassIrrep), pointer        :: Irrep
     integer                          :: Lmax
     character(len=:)   , allocatable :: StorageDir
     character(len=:)   , allocatable :: Root
   contains
     procedure, public :: SetGroup         => ClassSymESpace_SetGroup
     procedure, public :: SetIrrep         => ClassSymESpace_SetIrrep
     procedure, public :: SetMultiplicity  => ClassSymESpace_SetMultiplicity
     procedure, public :: SetRoot          => ClassSymESpace_SetRoot
     procedure, public :: GetIrrep         => ClassSymESpace_GetIrrep
     procedure, public :: GetIrrepLabel    => ClassSymESpace_GetIrrepLabel
     !.. Overwritten
     procedure, public :: free             => ClassSymESpace_Free
     procedure, public :: show             => ClassSymESpace_Show
     procedure, public :: GetNumChannels   => ClassSymESpace_GetNumChannels
     procedure, public :: GetChannel       => ClassSymESpace_GetChannel
     procedure, public :: GetLabel         => ClassSymESpace_GetLabel
     procedure, public :: GetRoot          => ClassSymESpace_GetRoot
     !..
     final             :: ClassSymESpace_Final
  end type ClassSymESpace


  character(len=*), private,parameter :: USE_ALL_XLMSYM   = "ALL_XLM"
  character(len=*)        , parameter :: RowsLabel        = 'ROWS'
  character(len=*)        , parameter :: ColsLabel        = 'COLUMNS'
  
  logical, parameter :: AIC_VIC_M_SWITCH = .TRUE.
  
  
  type, private :: ChanContainer
     class(ClassCloseCouplingChannel), pointer :: Ch
   contains
     final :: ChanContainerFinal
  end type ChanContainer

  type, private, extends( ClassSymESpace ) :: ClassSymIonCh
     private
     type(ClassParentIon)          :: ParentIon
     integer                       :: NaiMOs ! N Act. Int. MOs
     integer                       :: NVirtIntCh
     integer                       :: NVirtExtCh
     integer                       :: NChans
     type(ChanContainer), allocatable :: Chanv(:)
   contains
     procedure, public :: free           => ClassSymIonCh_Free
     procedure, public :: show           => ClassSymIonCh_Show
     procedure, public :: GetNumChannels => ClassSymIonCh_GetNumChannels
     procedure, public :: GetChannel     => ClassSymIonCh_GetChannel
     !
     procedure, public :: init           => ClassSymIonCh_Init
     procedure, public :: GetNChans      => ClassSymIonCh_GetNChans
     procedure, public :: ParseConfLine  => ClassSymIonCh_ParseConfLine
     procedure, public :: GetPILabelFun  => ClassSymIonCh_GetPILabelFun
     procedure, public :: GetPartialSize => ClassSymIonCh_GetPartialSize
     procedure, public :: Condition      => ClassSymIonCh_Condition
  end type ClassSymIonCh

  
  !.. Symmetric Close-Coupling Space
  type, public, extends(ClassSymESpace) :: ClassSymCCSpace
     private
     integer                          :: NIons
     type(ClassSymIonCh), allocatable :: IonicSpaceV(:)
     !.. AIC Conditioning (involves multiple ions)
     integer                          :: NLinIndepAIC
     integer                          :: NLinDepAIC
     type(ClassMatrix)                :: T_aiM_aiMc
     real(kind(1d0)), allocatable     :: hEval_aiMc(:)
   contains
     procedure, public :: free             => ClassSymCCSpace_Free
     procedure, public :: show             => ClassSymCCSpace_Show
     procedure, public :: GetNumChannels   => ClassSymCCSpace_GetNumChannels
     procedure, public :: GetChannel       => ClassSymCCSpace_GetChannel
     !
     procedure, public :: parseConfigFile  => ClassSymCCSpace_ParseConfigFile
     procedure, public :: GetSize          => ClassSymCCSpace_GetSize
     procedure, public :: Condition        => ClassSymCCSpace_Condition
     final             :: ClassSymCCSpace_Final
  end type ClassSymCCSpace


  logical, private :: SESSES_build_switch(3,3,3) = .TRUE.

  public :: GetOpXlm
  public :: GetOpIrrep
  public :: ValidSymmetries

  public :: Set_SESSES_switch
  
contains

  subroutine Set_SESSES_switch( type, switch )
    character(len=*), intent(in) :: type
    logical         , intent(in) :: switch(3,3,3)
    if(type .is. "build")then
       SESSES_build_switch = switch
    endif
  end subroutine Set_SESSES_switch
  
  !-----------------------------------------------------------------------
  ! Symmetric Electronic Space
  subroutine ClassSymESpace_Show( Self, unit )
    class(ClassSymESpace), intent(in) :: Self
    integer, optional                   , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    write(outunit,"(a)")              "  Symmetric El. Space info "
    write(outunit,"(a,i4)"          ) "  Total Charge             :",Self%Charge
    write(outunit,"(a,i4)"          ) "  Space Multiplicity       :",Self%Multiplicity
    write(outunit,"(a,i4)"          ) "  Maximum Angular Momentum :",Self%Lmax
    write(outunit,"(a,a)"           ) "  Space Storage Dir        :",Self%StorageDir
    write(outunit,"(a)",advance="no") "  Space Symmetry           :"
    call Self%Irrep%show( outunit )
  end subroutine ClassSymESpace_Show
  subroutine ClassSymESpace_Free( Self )
    class(ClassSymESpace), intent(inout) :: Self
    Self%Lmax   = -1
    Self%Charge = 0
    Self%Multiplicity = 0
    Self%Group => NULL()
    Self%Irrep => NULL()
    if ( allocated(Self%Root      ) ) deallocate( Self%Root       )
    if ( allocated(Self%StorageDir) ) deallocate( Self%StorageDir )
  end subroutine ClassSymESpace_Free
  subroutine ClassSymESpace_Final( Self )
    type(ClassSymESpace) :: Self
    call Self%free()
  end subroutine ClassSymESpace_Final
  subroutine ClassSymESpace_SetGroup( Self, Group )
    class(ClassSymESpace), intent(inout) :: Self
    type(ClassGroup), target,             intent(in)    :: Group
    if ( associated(Self%Group) ) Self%Group => NULL()
    Self%Group => Group
  end subroutine ClassSymESpace_SetGroup
  subroutine ClassSymESpace_SetIrrep( Self, Irrep )
    class(ClassSymESpace), intent(inout) :: Self
    type(ClassIrrep), target,             intent(in)    :: Irrep
    if ( associated(Self%Irrep) ) Self%Irrep => NULL()
    Self%Irrep => Irrep
  end subroutine ClassSymESpace_SetIrrep
  subroutine ClassSymESpace_SetMultiplicity( Self, MUltiplicity )
    class(ClassSymESpace), intent(inout) :: Self
    integer,                              intent(in)    :: Multiplicity
    Self%Multiplicity = Multiplicity
  end subroutine ClassSymESpace_SetMultiplicity
  subroutine ClassSymESpace_SetRoot( Self, Root )
    class(ClassSymESpace), intent(inout) :: Self
    character(len=*),                     intent(in)    :: Root
    self%Root = Root
  end subroutine ClassSymESpace_SetRoot
  function ClassSymESpace_GetIrrep( Self ) result( SymIrrep )
    class(ClassSymESpace), target, intent(in) :: Self
    type(ClassIrrep), pointer :: SymIrrep
    allocate( SymIrrep, source = Self%Irrep )
  end function ClassSymESpace_GetIrrep
  function ClassSymESpace_GetIrrepLabel( Self ) result( Label )
    class(ClassSymESpace), intent(in) :: Self
    character(len=:), allocatable :: Label
    allocate( Label, source = Self%Irrep%GetName() )
  end function ClassSymESpace_GetIrrepLabel
  integer function ClassSymESpace_GetNumChannels( Self ) result(Nchan)
    class(ClassSymESpace), intent(in) :: Self
    nchan = 0
  end function ClassSymESpace_GetNumChannels
  function ClassSymESpace_GetChannel( Self, iChan ) result(Chan)
    class(ClassSymESpace), intent(in) :: Self
    integer              , intent(in) :: iChan
    class(ClassCloseCouplingChannel), pointer :: Chan
    chan => NULL()
  end function ClassSymESpace_GetChannel
  function ClassSymESpace_GetLabel( Self ) result( Label )
    class(ClassSymESpace), intent(in) :: Self
    character(len=:), allocatable :: Label
    allocate( Label, source = AlphabeticNumber(Self%Multiplicity)//Self%Irrep%GetName() )
  end function ClassSymESpace_GetLabel
  function ClassSymESpace_GetRoot( Self ) result( string )
    class(ClassSymESpace), intent(in) :: Self
    character(len=:), allocatable :: string
    allocate( string, source = self%Root ) 
  end function ClassSymESpace_GetRoot



  !-----------------------------------------------------------------------
  ! Symmetric CloseCoupling Space

  subroutine ClassSymCCSpace_ParseConfigFile( Self, &
       FileName, iflag )
    use moduleParameterList
    use moduleString
    use moduleIO
    class(ClassSymCCSpace), intent(inout) :: Self
    character(len=*)                    , intent(in)    :: FileName
    integer                             , intent(out)   :: iflag ! 0 => OK, n => some error
    !
    integer :: j, ichar, ichar2, dichar
    character(len=1024)      :: IonStrn
    character(len=:), allocatable :: FullText
    character(len=:), allocatable :: GroupLabel, GeneralLabel, SymStrn
    integer :: NumChannels
    iflag = 0 

    !.. Fetch GlobalVariables from Configuration file
    !..

    !.. 1. Fetch the full text of the configuration file
    !.. 
    call GetFullText( FileName, FullText )
    call SetStringToUppercase( FullText )

    !.. 3. Extract global variables
    !..
    call FetchGlobalVariable( FullText, "GROUP"     , GroupLabel  )
    if(.not.allocated(GroupLabel)) call Assert("Group label missing in "//trim(FileName))

    call FetchGlobalVariable( FullText, "LMAX"      , GeneralLabel )
    if(.not.allocated(GeneralLabel)) call Assert("lmax label missing in "//trim(FileName))
    read(GeneralLabel,*) Self%Lmax

    call FetchGlobalVariable( FullText, "CHARGE", GeneralLabel )
    if(.not.allocated(GeneralLabel)) call Assert("CHARGE label missing in "//trim(FileName))
    read(GeneralLabel,*) Self%Charge

    !.. Define the name of the storage directory
    SymStrn=Self%GetLabel()
    Self%StorageDir=self%Root//"/"//SymStrn

    !.. 4. Isolate the irrep in the file 
    !..
    call SetStringToUppercase( SymStrn )
    ichar = index(FullText,"["//SymStrn//"]")
    if(ichar<1)then
       iflag = 1 
       return
    endif
    !
    ichar = ichar+index(FullText(ichar+1:),"{")
    FullText=adjustl(FullText(ichar+1:))
    !
    ichar = index(FullText,"}")
    FullText=adjustl(FullText(:ichar-1))

    !.. Count the number of Ionic Channels
    !..
    NumChannels = 0
    ichar=0
    IonCycle : do 
       dichar = index(FullText(ichar+1:),")")
       if(dichar<1)exit IonCycle
       ichar = ichar + dichar
       NumChannels = NumChannels+1
    enddo IonCycle
    Self%NIons = NumChannels
    allocate( Self%IonicSpaceV(NumChannels) )

    !.. Parse the description of each ionic channel
    !..
    ichar = 0
    do j = 1, Self%NIons

       ichar2 = ichar
       ichar = ichar+index(FullText(ichar+1:),")")
       IonStrn=adjustl(FullText(ichar2+1:ichar))
!write(*,*) trim(IonStrn)
!write(*,*) "---------",ichar,ichar2
       call Self%IonicSpaceV(j)%Init( &
            Self%Lmax        , & 
            Self%Multiplicity, &
            Self%Irrep       , &
            Self%Charge      , &
            Self%StorageDir  )
       call Self%IonicSpaceV(j)%ParseConfLine( trim(IonStrn) )

    enddo

  end subroutine ClassSymCCSpace_ParseConfigFile

  subroutine ClassSymCCSpace_Show( Self, unit )
    class(ClassSymCCSpace), intent(in) :: Self
    integer, optional     , intent(in) :: unit
    integer :: outunit, iChan
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    call self%ClassSymESpace%show( unit )
    do iChan = 1, Self%NIons
       call Self%IonicSpaceV(iChan)%show( outunit )
    enddo
  end subroutine ClassSymCCSpace_Show

  subroutine ClassSymCCSpace_Free( Self )
    class(ClassSymCCSpace), intent(inout) :: Self
    integer :: iIon
    do iIon=1,self%NIons
       call Self%IonicSpaceV(iIon)%Free()
    enddo
    Self%NIons = 0
    if ( allocated(Self%IonicSpaceV) ) deallocate( Self%IonicSpaceV )
    call self%ClassSymESpace%free()
  end subroutine ClassSymCCSpace_Free

  subroutine ClassSymCCSpace_Final( Self )
    type(ClassSymCCSpace) :: Self
    call Self%free()
  end subroutine ClassSymCCSpace_Final


  !.. Return the total size for the whole space, or for a subspace 
  !   in which the kind of channel and/or the ion is as specified
  integer function ClassSymCCSpace_GetSize( Self, &
       Required_channel, Required_Ion_Label ) result( isize )
    class(ClassSymCCSpace), intent(in) :: Self
    class(ClassCloseCouplingChannel), optional, intent(in) :: Required_channel
    character(len=*)          , optional, intent(in) :: Required_Ion_Label
    character(len=:), allocatable     :: Ion_Label
    integer                           :: iIon
    type(ClassActiveInternalChannel)  :: aic
    type(ClassVirtualInternalChannel) :: vic
    type(ClassVirtualExternalChannel) :: vec
    isize = 0
    do iIon = 1, Self%NIons
       if( present( Required_Ion_Label ) )then
          Ion_Label = Self%IonicSpaceV(iIon)%GetPILabelFun()
          if( Ion_Label .isnt. Required_Ion_Label ) cycle
       endif
       if( present( Required_channel ) )then
          isize = isize + Self%IonicSpaceV(iIon)%GetPartialSize(Required_channel)
       else
          isize = isize + Self%IonicSpaceV(iIon)%GetPartialSize(aic)
          isize = isize + Self%IonicSpaceV(iIon)%GetPartialSize(vic)
          isize = isize + Self%IonicSpaceV(iIon)%GetPartialSize(vec)
       end if
    end do
  end function ClassSymCCSpace_GetSize

  integer function ClassSymCCSpace_GetNumChannels( Self ) result(Nchan)
    class(ClassSymCCSpace), intent(in) :: Self
    integer :: i, nactintch, nvirtintch,nVirtExtCh, nchan_part
    nchan = 0
    do i = 1, Self%NIons
       call Self%IonicSpaceV(i)%GetNChans( nactintch, nvirtintch, nVirtExtCh, nchan_part )
       nchan = nchan + nchan_part
    end do
  end function ClassSymCCSpace_GetNumChannels

  function ClassSymCCSpace_GetChannel( Self, iChan ) result(Chan)
    class(ClassSymCCSpace)          , intent(in) :: Self
    integer                         , intent(in) :: iChan
    class(ClassCloseCouplingChannel), pointer    :: Chan
    integer :: iIon, nactintch, nvirtintch, nVirtExtCh, nchan_part,nchan, iChanRel
    nchan=0
    do iIon = 1, Self%NIons
       call Self%IonicSpaceV(iIon)%GetNChans( nactintch, nvirtintch, nVirtExtCh, nchan_part )
       if( iChan <= nchan + nchan_part )then
          iChanRel = iChan - nchan
          Chan => Self%IonicSpaceV(iIon)%GetChannel(iChanRel)
          exit
       endif
       nchan = nchan + nchan_part
    end do

  end function ClassSymCCSpace_GetChannel

  subroutine ClassSymCCSpace_Condition( self, threshold )
    class(ClassSymCCSpace), intent(inout) :: self
    real(kind(1d0))      , intent(in)    :: threshold

    real(kind(1d0)), parameter    :: NEG_THR_DEF =-1.d-12

    type(ClassSESSESBlock)        :: SBlk, HBlk
    type(ClassMatrix)             :: SMat, HMat, sEvec, hEvec
    real(kind(1d0)), allocatable  :: sEval(:), hEval(:), sfact(:)
    logical                       :: lTasks(N_SESSES_ID)
    integer                       :: iIon

    call SBlk%Free()
    call HBlk%Free()
    call SBlk%Init(self,self,"S",self%Root)
    call HBlk%Init(self,self,"H",self%Root)

    !.. First, it conditions the aic space, which requires
    !   gathering blocks from different ions
    BLOCK 
      type(ClassActiveInternalChannel) :: aic
      integer                          :: NLinDep, i, nLinInd
      !.. Load and diagonalize the aic overlap 
      call Execute_Command_Line("mkdir -p "//self%StorageDir//"/aiM")
      lTasks=.FALSE.; lTasks(SESSES_LOAD)=.TRUE.
      call SBlk%Driver(lTasks,aic,aic)
      call SBlk%Assemble(SMat,aic,aic)
      call SMat%Write(self%StorageDir//"/aiM/S_aiM_aiM","formatted")
      call SMat%Diagonalize(sEval,sEvec)
      call sEvec%Write(self%StorageDir//"/aiM/S_aiM_aiMc_evec","formatted")
      call SaveVector(self%StorageDir//"/aiM/S_aiM_eval",sEval,"formatted")
      !.. Drop the solutions below threshold
      do i=1,size(sEval)
         if(sEval(i)<NEG_THR_DEF) call Assert("S eval < 0")
         if(sEval(i)>threshold)exit
      enddo
      NLinDep = i-1
      NLinInd = size(sEval) - NLinDep
      write(*,"(a)") "  n lin dep in aic space : ", NLinDep
      call sEvec%Drop("columns",1,NLinDep)
      do i = 1, size(sEval) - NLinDep
         sEval(i) = sEval(i+NLinDep)
      enddo
      allocate(sfact,source=1.d0/sqrt(sEval(1:NLinInd)))
      call sEvec%Multiply(sfact,"Right")
      call sEvec%Write(self%StorageDir//"/aiM/T_aiM_aiMc","formatted")
      !.. Load, transform, and diagonalize the aic hamiltonian
      !*** MUST UNIT-TEST
      call HBlk%Driver(lTasks,aic,aic)
      call HBlk%Assemble(HMat,aic,aic)
      call HMat%Write(self%StorageDir//"/aiM/H_aiM_aiM","formatted")
      call HMat%Multiply(sEvec,"Right","N")
      call HMat%Multiply(sEvec,"Left","T")
      call HMat%Write(self%StorageDir//"/aiM/H_aiMc_aiMc","formatted")
      call HMat%Diagonalize(hEval,hEvec)
      call hEvec%Write(self%StorageDir//"/aiM/H_aiMc_aiMc_evec","formatted")
      call SaveVector(self%StorageDir//"/aiM/H_aiMc_eval",hEval,"formatted")
      call hEvec%Multiply(sEvec,"Left","N")
      call hEvec%Write(self%StorageDir//"/aiM/H_aiM_aiMc_evec","formatted")
      deallocate(sEval,sfact)
      call sEvec%Free()
    END BLOCK

    !.. Condition the individual ionic channels
    do iIon = 1, self%NIons
       call self%IonicSpaceV(iIon)%Condition()
    end do

  end subroutine ClassSymCCSpace_Condition

  !-------------------------------------------------------------------------------
  ! General close-coupling Block
  !> Gets the Xlm associated with the operator, assuming that 
  !! everything which is not S, H, K or CAP is a dipole.
  function GetOpXlm( OpLabel ) result (OpXlm)
    character(len=*)             , intent(in)    :: OpLabel
    type(ClassXlm) :: OpXlm
    if ( (OpLabel .is. OverlapLabel) .or. &
         (OpLabel .is. HamiltonianLabel) .or. &
         (OpLabel .is. KinEnergyLabel) .or. &
         (OpLabel(:min(len(OpLabel),len(CAPLabel))) .is. CAPLabel)) then
       call OpXlm%Init(0,0)
    else
       if     ( ( OpLabel .is. 'X' ) .or. ( OpLabel .is. 'DX' ) ) then
          call OpXlm%Init(1,1)
       elseif ( ( OpLabel .is. 'Y' ) .or. ( OpLabel .is. 'DY' ) ) then
          call OpXlm%Init(1,-1)
       elseif ( ( OpLabel .is. 'Z' ) .or. ( OpLabel .is. 'DZ' ) ) then
          call OpXlm%Init(1,0)
       else
          call Assert( 'Invalid orientation to get the dipole operator Xlm'//&
               ', it must be either "x", "y" or "z".' )
       endif
    end if
  end function GetOpXlm

  function GetOpIrrep( Group, OpLabel ) result (OpIrrep)
    class(ClassGroup)            , intent(in)    :: Group
    character(len=*)             , intent(in)    :: OpLabel
    type(ClassIrrep), pointer :: OpIrrep
    type(ClassXlm)   :: OpXlm
    OpXlm = GetOpXlm(OpLabel)
    allocate( OpIrrep, source = OpXlm%GetIrrep(Group) )
  end function GetOpIrrep

  !------------------------------------------------
  ! Methods for ClassCCBlock
  !-----------------------------------------------

  subroutine ClassCCBlockCopy( self, from )
    class(ClassCCBlock), intent(inout) :: self
    type (ClassCCBlock), intent(in)    :: from
    call self%Free()
    call self%Init(from%IonChBra,from%IonChKet,from%OpLabel, from%Root, from%IDLabel )
    self%Block = from%Block
  end subroutine ClassCCBlockCopy

  subroutine ClassCCBlockScale( self, x )
    class(ClassCCBlock), intent(inout) :: self
    real(kind(1d0))         , intent(in)    :: x
    call self%Block%Multiply(x)
  end subroutine ClassCCBlockScale


  subroutine ClassCCBlockInit( self, IonChBra, IonChKet, OpLabel, Root, IDLabel )
    class(ClassCCBlock)            , intent(inout) :: self
    class(ClassCloseCouplingChannel), target  , intent(in)    :: IonChBra
    class(ClassCloseCouplingChannel), target  , intent(in)    :: IonChKet
    character(len=*)                    , intent(in)    :: OpLabel
    character(len=*)                    , intent(in)    :: Root
    character(len=*)          , optional, intent(in)    :: IDLabel

    call self%free()
    self%FormattedWrite = CCBLOCK_FORMATTED_WRITE
    allocate( self%OpLabel , source = OpLabel )
    if(present(IDLabel)) allocate( self%IDLabel , source = IDLabel )
    allocate( self%OpXlm   , source = GetOpXlm( OpLabel ) )
    allocate( self%Root    , source = Root )
!!$    allocate( self%StorageDir, source = AddSlash(Root)//IonChBra%GetLabel()//"_"//IonChKet%GetLabel() )
    self%IonChBra => IonChBra
    self%IonChKet => IonChKet
    select type ( ptr => self%IonChBra )
    type is ( ClassActiveInternalChannel )
       self%nOrbBra = OrbitalBasis%GetNOrb(self%IonChBra%GetOrbIrrep(),"active")
    type is ( ClassVirtualInternalChannel )
       if(ptr%IsMolecular())then
          self%nOrbBra = OrbitalBasis%GetNOrb(self%IonChBra%GetOrbIrrep(),"virtual")
       else
          self%nOrbBra = OrbitalBasis%GetNOrb(self%IonChBra%GetOrbIrrep(),"hybrid")
       endif
    type is ( ClassVirtualExternalChannel )
       self%nOrbBra = OrbitalBasis%GetNOrb(self%IonChBra%GetOrbIrrep(),"spline")
    end select
    select type ( ptr => self%IonChKet )
    type is ( ClassActiveInternalChannel )
       self%nOrbKet = OrbitalBasis%GetNOrb(self%IonChKet%GetOrbIrrep(),"active")
    type is ( ClassVirtualInternalChannel )
       if(ptr%IsMolecular())then
          self%nOrbKet = OrbitalBasis%GetNOrb(self%IonChKet%GetOrbIrrep(),"virtual")
       else
          self%nOrbKet = OrbitalBasis%GetNOrb(self%IonChKet%GetOrbIrrep(),"hybrid")
       endif
    type is ( ClassVirtualExternalChannel )
       self%nOrbKet = OrbitalBasis%GetNOrb(self%IonChKet%GetOrbIrrep(),"spline")
    end select

    self%initialized = .true.
  end subroutine ClassCCBlockInit


  logical function ClassCCBlockIsInitialized( self ) result( IsInitialized )
    class(ClassCCBlock), intent(in) :: self
    IsInitialized = self%initialized
  end function ClassCCBlockIsInitialized


  subroutine ClassCCBlockFinal( self )
    type(ClassCCBlock) :: self
    call self%Free()
  end subroutine ClassCCBlockFinal


  subroutine ClassCCBlockFree( self )
    class(ClassCCBlock), intent(inout) :: self
    if ( allocated(self%Root      ) ) deallocate( self%Root    )
    if ( allocated(self%OpLabel   ) ) deallocate( self%OpLabel )
    if ( allocated(self%IDLabel   ) ) deallocate( self%IDLabel )
    if ( associated(self%OpXlm    ) ) deallocate( self%OpXlm   )
    self%OpXlm          => NULL()
    self%IonChBra       => NULL()
    self%IonChKet       => NULL()
    self%FormattedWrite = .false.
    self%initialized    = .false.
    call self%Block%Free()
  end subroutine ClassCCBlockFree

  subroutine ClassCCBlockLoad( self, bstat ) 
    class(ClassCCBlock), intent(inout) :: self
    integer                 , intent(out)   :: bstat
    character(len=:)        , allocatable   :: BlockName
    integer :: uid
    bstat=1
    if(.not.self%initialized) call Assert("Block not initialized")
    allocate( BlockName, source = self%GetBlockName() )
    call CheckFileMustBePresent( BlockName )
    if ( self%FormattedWrite ) then
       call OpenFile( BlockName, uid, 'read', 'formatted' )
    else
       call OpenFile( BlockName, uid, 'read', 'unformatted' )
    end if
    call self%Block%Read( uid )
    close( uid )
    bstat=0
  end subroutine ClassCCBlockLoad

  subroutine ClassCCBlockSave( self )
    class(ClassCCBlock), intent(in) :: self
    !
    character(len=:), allocatable :: FileName
    allocate( FileName, source = self%GetBlockName() )
    if ( self%FormattedWrite ) then
       call self%Block%Write( FileName, 'formatted' )
    else
       call self%Block%Write( FileName, 'unformatted' )
    end if
  end subroutine ClassCCBlockSave

  function ClassCCBlockGetLabel( self, side ) result(ChLabel)
    class(ClassCCBlock), intent(in)  :: self
    character(len=*)   , intent(in)  :: side  
    character(len=:), allocatable :: ChLabel
    if(side.is."bra")then
       ChLabel = self%IonChBra%GetPILabel()//GetccExt(self%IonChBra)
    elseif(side.is."ket")then
       ChLabel = self%IonChKet%GetPILabel()//GetccExt(self%IonChKet)
    endif
  end function ClassCCBlockGetLabel

  function ClassCCBlockGetIonBraLabel( self ) result(IonLabel)
    class(ClassCCBlock)     , intent(in)  :: self
    character(len=:), allocatable :: IonLabel
    IonLabel = self%IonChBra%GetPILabel() 
  end function ClassCCBlockGetIonBraLabel

  function ClassCCBlockGetIonKetLabel( self ) result(IonLabel)
    class(ClassCCBlock)     , intent(in)  :: self
    character(len=:), allocatable :: IonLabel
    IonLabel = self%IonChKet%GetPILabel() 
  end function ClassCCBlockGetIonKetLabel

  subroutine ClassCCBlockGetSize( self, nr, nc )
    class(ClassCCBlock), intent(in)  :: self
    integer                 , intent(out) :: nr, nc
    nr = self%nOrbBra
    nc = self%nOrbKet
  end subroutine ClassCCBlockGetSize

  logical function ClassCCBlockIsType( self, bChan, kChan ) result(isType)
    class(ClassCCBlock)            , intent(inout) :: self
    class(ClassCloseCouplingChannel), optional, intent(in)  :: bChan, kChan
    isType=.TRUE.
    if(present(bChan))then
       isType=isType.and.same_type_as(self%IonChBra,bChan)
    endif
    if(present(kChan))then
       isType=isType.and.same_type_as(self%IonChKet,kChan)
    endif
  end function ClassCCBlockIsType

  function ClassCCBlockGetStorageDir( self ) result(Dir)
    class(ClassCCBlock)            , intent(in) :: self
    character(len=:), allocatable :: Dir
    if ( .not.allocated(Self%root) ) call Assert( &
         'Root directory not defined.' )
    allocate( Dir, source = AddSlash(self%Root) // &
         self%IonChBra%GetLabel() // "_" // self%IonChKet%GetLabel() )
  end function ClassCCBlockGetStorageDir


  subroutine ClassCCBlockTranspose( self ) 
    class(ClassCCBlock), intent(inout) :: self
    class(ClassCloseCouplingChannel), pointer :: cptr
    integer :: n
    cptr          => self%IonChBra
    self%IonChBra => self%IonChKet
    self%IonChKet => cptr
    n             =  self%nOrbBra
    self%nOrbBra  =  self%nOrbKet
    self%nOrbKet  =  n
    call self%Block%Transpose()
  end subroutine ClassCCBlockTranspose


  function ConvertLabel_SES_to_INT( SES_OP_label ) result( INT_OP_label )
    character(len=*), intent(in)  :: SES_OP_label
    character(len=:), allocatable :: INT_OP_LABEL
   
    !.. Determines the operator label in module Integrals
    select case( SES_OP_label )
    case ( OverlapLabel )
       INT_OP_LABEL = I1B_ID_LIST( I1B_OVERLAP )
    case ( KinEnergyLabel )
       INT_OP_LABEL = I1B_ID_LIST( I1B_KINETIC )
    case ( HamiltonianLabel )
       INT_OP_LABEL = I1B_ID_LIST( I1B_HAMILTO )
    case ( X_Label )
       INT_OP_LABEL = I1B_ID_LIST( I1B_COORD_X )
    case ( Y_Label )
       INT_OP_LABEL = I1B_ID_LIST( I1B_COORD_Y )
    case ( Z_Label )
       INT_OP_LABEL = I1B_ID_LIST( I1B_COORD_Z )
    case ( DX_Label )
       INT_OP_LABEL = I1B_ID_LIST( I1B_NABLA_X )
    case ( DY_Label )
       INT_OP_LABEL = I1B_ID_LIST( I1B_NABLA_Y )
    case ( DZ_Label )
       INT_OP_LABEL = I1B_ID_LIST( I1B_NABLA_Z )
    case DEFAULT
       INT_OP_LABEL = "___"
    end select
    
  end function ConvertLabel_SES_to_INT


  subroutine ClassCCBlockBuild( self, bstat, switch_, timeMatrix )
    class(ClassCCBlock)      , intent(inout)  :: self
    integer                  , intent(out)    :: bstat
    logical        , optional, intent(in)     :: switch_(3,3,3) !{aic,vic,vec} x {aic,vic,vec} x {S,M,B}
    real(kind(1d0)), optional, intent(inout)  :: timeMatrix(3,3)
    !
    logical                       :: IsOverlap, IsHamiltonian, IsOneBody, IsTwoBody
    character(len=:), allocatable :: IonBraLabel, IonKetLabel, IntOpLabel
    type(ClassIrrep)              :: IrrOrbBra, IrrOrbKet, IrrIonBra, IrrIonKet, IrrAB, OpIrrep
    type(ClassIrrep), allocatable :: vIrreps(:)
    integer                       :: SA2, SB2, S2, iIrrOrbBra, iIrrOrbKet, nInctvOrbBra, nInctvOrbKet
    integer         , allocatable :: nactive(:), ninactive(:)
    real(kind(1d0))               :: etime
    logical                       :: switch(3,3,3)
    !
    if( self%nOrbBra * self%nOrbKet == 0 )then
       call self%block%InitFull(self%nOrbBra,self%nOrbKet)
       bstat=0
       return
    endif

    if(present(switch_))then
       switch=switch_
    else
       switch=.TRUE.
    endif
    !*** Thes gave a wrong behaviour of swich when self%Build was called (autoreferentially)
    !switch=.TRUE.   
    !write(*,*) present(switch_),"\|/",switch_
    !if(present(switch_))switch=switch_


    !******************************************************
    !******************************************************
    !******************************************************
    !******************************************************
    !
    !  WARNING: MUST MAKE SURE THAT THE OPERATORS WE ARE
    !  COMPUTING FOR THE ION FOLLOWS THE SAME CONVENTIONS
    !  AS THE OPERATOR FOR THE PHOTOELECTRON. IN PARTICULAR,
    !  THE DIPOLE MOMENTS AND MULTIPOLES CONTAIN THE CHARGE
    !  OF THE PARTICLES IN THEM, AND HAVE A NUCLEAR (POSITIVE
    !  CHARGE) AND AN ELECTRONIC (NEGATIVE CHARGE) COMPONENT
    !  THAT MUST BE TREATED CONSISTENTLY.
    !
    !******************************************************
    !******************************************************
    !******************************************************
    !******************************************************

    bstat = 1

    IsOverlap     = ( self%OpLabel == OverlapLabel ) 
    IsHamiltonian = ( self%OpLabel == HamiltonianLabel ) 
    IsOneBody     = .not. IsOverlap
    IsTwoBody     = IsHamiltonian

    call TDM_Manager%GetNActiveEls( nInactive, nActive )

    vIrreps      = GlobalGroup%GetIrrepList()
    IonBraLabel  = self%IonChBra%GetPILabel()
    IonKetLabel  = self%IonChKet%GetPILabel()
    IrrOrbBra    = self%IonChBra%GetOrbIrrep()
    IrrOrbKet    = self%IonChKet%GetOrbIrrep()
    iIrrOrbBra   = GlobalGroup%GetIrrepIndex(IrrOrbBra)
    iIrrOrbKet   = GlobalGroup%GetIrrepIndex(IrrOrbKet)
    nInctvOrbBra = nInactive( GlobalGroup%GetIrrepIndex( IrrOrbBra ) )
    nInctvOrbKet = nInactive( GlobalGroup%GetIrrepIndex( IrrOrbKet ) )
    IrrIonBra    = self%IonChBra%GetIonIrrep()
    IrrIonKet    = self%IonChKet%GetIonIrrep()
    IrrAB        = IrrIonBra * IrrIonKet
    OpIrrep      = self%OpXlm%GetIrrep(GlobalGroup)
    SA2          = self%IonChBra%GetPIMultiplicity()-1
    SB2          = self%IonChKet%GetPIMultiplicity()-1
    S2           = self%IonChKet%GetTotMultiplicity()-1

    IntOpLabel   = ConvertLabel_SES_to_INT( self%OpLabel )

    !*** delete this line
    !write(*,*) self%IonChBra%GetLabel()
    
    select type( pB => self%IonChBra )
    type is ( ClassActiveInternalChannel )

       select type( pK => self%IonChKet )
       type is ( ClassActiveInternalChannel )

          call ElapsedTime()

          !.. AIC-AIC (CASE 1)

          if( IsOverlap .and. switch (1,1,1) )then
             !   OVERLAP (0B)
             !   \langle A,p; S\Sigma|B,q; S\Sigma\rangle =\delta_{AB}S_{pq}+B^{BA}_{pq}
             !
             !
             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing ....................................V
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             aic_aic_S : BLOCK
               integer :: nc, nr
               call self%GetSize(nr,nc)
               call self%block%InitFull(nr,nc)
               self%Block = TDM_Manager%GetB( IonBraLabel, IonKetLabel, IrrOrbBra )
               if(IonBraLabel .is. IonKetLabel) call self%Block%AddDiagonal( 1.d0 )
             end BLOCK aic_aic_S
          endif

          if( IsOneBody .and. .not. IsHamiltonian .and. switch (1,1,2) )then
             !
             ! 1.  Write code ............................................V
             ! 2.  Compile in debug ......................................V
             ! 3.  Runs w/o crashing .....................................V
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   1-BODY generic
             !   \langle A,p; S\Sigma | \hat{O} |B,q; S\Sigma\rangle =
             !   \delta_{S_A S_B}S_{pq}\Pi_{S_A}^{-1}\langle A \| \hat{O}\|B\rangle +
             !   \delta_{AB} o_{pq}+{\sum_r}'\left(B^{BA}_{pr}o_{rq} + o_{pr}B^{BA}_{rq}\right)+
             !   2\left(\sum_xo_{xx}\right)B^{BA}_{pq} + {\sum_{rs}}'o_{rs}P^{BA}_{ps,qr}.
             !
             !.. Note: in the AIC-AIC case, the Hamiltonian has already been computed by Lucia,
             !   and here it is only reconstructed with the right geometrical coefficients.
             !   For this reason, there is no Hamiltonian subsection in this one-body block.
             aic_aic_M : BLOCK
               type(ClassMatrix)             :: B, C, Op
               type(ClassMatrix4D), pointer  :: Q
               type(ClassIrrep)              :: IrrOrbMid, Irr2
               integer                       :: nr, nc, ix, iIrrX
               integer                       :: nInctvOrbMid, iIrrOrbMid
               integer                       :: iIrr1, iIrr2, ivecmi(4), ivecma(4)
               real(kind(1d0))               :: O_AB, O_XX
               
               call self%GetSize(nr,nc)
               
               call C%InitFull(nr,nc)
               C=0.d0
               

               !.. \delta_{S_A S_B}S_{pq} \langle A | \hat{O} | B \rangle
               if(SA2 == SB2)then
                  O_AB=TDM_Manager%GetMatEl(IntOpLabel, IonBraLabel, IonKetLabel)
                  call C%SetDiagonal(O_AB)
                  call self%Block%Add(C)
               endif

               !.. \delta_{AB} o_{pq}
               if(IonBraLabel .is. IonKetLabel)then
                  Op= GlobalIntegral%Get1B( IntOplabel,"MO_MO", iIrrOrbBra, iIrrOrbKet )
                  call Op%Drop("rows"   ,1,nInctvOrbBra) 
                  call Op%Drop("columns",1,nInctvOrbKet) 
                  call self%Block%Add(Op)
               endif

               !.. B^{BA}_{pr}o_{rq} 
               IrrOrbMid  = OpIrrep * IrrOrbKet
               iIrrOrbMid = GlobalGroup%GetIrrepIndex( IrrOrbMid )
               Op= GlobalIntegral%Get1B( IntOpLabel,"MO_MO", iIrrOrbMid, iIrrOrbKet )
               nInctvOrbMid = nInactive( GlobalGroup%GetIrrepIndex( IrrOrbMid ) )
               call Op%Drop("rows"   ,1,nInctvOrbMid) 
               call Op%Drop("columns",1,nInctvOrbKet) 
               B = TDM_Manager%GetB( IonBraLabel, IonKetLabel, IrrOrbBra )
               call B%Multiply(Op,"Right","N")
               call self%block%Add(B)

               !.. o_{pr}B^{BA}_{rq}
               IrrOrbMid  = IrrOrbBra * OpIrrep
               iIrrOrbMid = GlobalGroup%GetIrrepIndex( IrrOrbMid )
               Op = GlobalIntegral%Get1B( IntOpLabel,"MO_MO", iIrrOrbBra, iIrrOrbMid )
               nInctvOrbMid = nInactive( GlobalGroup%GetIrrepIndex( IrrOrbMid ) )
               call Op%Drop("rows"   ,1,nInctvOrbBra) 
               call Op%Drop("columns",1,nInctvOrbMid) 
               B = TDM_Manager%GetB( IonBraLabel, IonKetLabel, IrrOrbMid )
               call B%Multiply(Op,"Left","N")
               call self%block%Add(B)
         
               !.. ( \sum_x o_{xx} ) B^{BA}_{pq}
               O_XX=0.d0
               do iIrrX = 1, size(vIrreps)
                  if( nInactive( iIrrX ) == 0 ) cycle
                  Op = GlobalIntegral%Get1B( IntOpLabel,"MO_MO", iIrrX, iIrrX )
                  do ix = 1, nInactive( iIrrX )
                     O_XX = O_XX + Op%Element(ix,ix)
                  enddo
               enddo
               B = TDM_Manager%GetB( IonBraLabel, IonKetLabel, IrrOrbBra )
               call B%Multiply(2.d0*O_XX)
               call self%block%Add(B)

               ! o_{sr} Q^{BA}_{sr,pq}
               ivecmi=1
               ivecma(3)=nr
               ivecma(4)=nc
               do iIrr1 = 1, size(vIrreps)
                  Irr2  = vIrreps(iIrr1) * OpIrrep
                  iIrr2 = GlobalGroup%GetIrrepIndex(Irr2)
                  ivecma(1)=nActive(iIrr1)
                  ivecma(2)=nActive(iIrr2)
                  Op = GlobalIntegral%Get1B( IntOpLabel,"MO_MO", iIrrOrbBra, iIrrOrbMid )
                  call Op%Drop("rows"   ,1,nInactive(iIrr1)) 
                  call Op%Drop("columns",1,nInactive(iIrr2)) 
                  Q  =>TDM_Manager%GetQ( IonBraLabel, IonKetLabel, iIrr1, iIrr2, iIrrOrbBra, S2 )
                  C  = Op%Contract4D12(Q,ivecmi,ivecma)
                  call self%Block%Add(C)
               enddo

               call C%Free()
               call B%Free()
               call Op%Free()
               
             end BLOCK aic_aic_M
          endif

          if( IsHamiltonian .and. switch (1,1,3) )then

             !
             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing ...?
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   HAMILTONIAN (1+2B)
             !   \langle A,p; S\Sigma | \hat{H} |B,q; S\Sigma\rangle = 
             !   =(-1)^{S_B+1/2+S}\sum_T \mathsf{H}^{AB}_{[p,q]_T} \Pi_{T}
             !    \sjs{T}{1/2}{1/2}{S}{S_B}{S_A}
             !
             aic_aic_H : BLOCK
               logical                       :: l0, l1
               real(kind(1d0))               :: c0, c1, sgn
               type(ClassMatrix)             :: H0, H1
               !
               !1. Fetch the matrix elements H
               call TDM_Manager%GetHAIC( IonBraLabel, IonKetLabel, self%IonChBra%GetOrbIrrep(), H0, H1 )
               l0=H0%IsInitialized()
               l1=H1%IsInitialized()
               if(.not.(l0.or.l1))then
                  call self%Block%InitFull(self%nOrbBra,self%nOrbKet)
                  self%Block=0.d0
                  exit aic_aic_H
               endif
               !
               sgn = 1.d0 - 2.d0 * mod((SB2+S2+2)/2,2)
               !
               if(l0)then
                  c0 = sgn * SixJSymbol_HalfHalf( &
                       0.d0, dble(SB2)/2.d0, dble(SA2)/2.d0, dble(S2)/2.d0 )
                  call H0%Multiply(c0)
                  self%Block = H0
               endif
               !
               if(l1)then
                  c1 = sgn * sqrt(3.d0) * SixJSymbol_HalfHalf( &
                       1.d0, dble(SB2)/2.d0, dble(SA2)/2.d0, dble(S2)/2.d0 )
                  call H1%Multiply(c1)
                  if(.not.l0)then
                     self%Block = H1
                  else
                     call self%Block%Add( H1 )
                  endif
               endif
               !
               call H0%Free()
               call H1%Free()
             end BLOCK aic_aic_H

          end if

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(1,1)=TimeMatrix(1,1)+etime

       type is ( ClassVirtualInternalChannel )

          call ElapsedTime()

          !.. AIC-VIC  (CASE 2)
          if( IsOverlap .and. switch (1,2,1) )then
             !
             ! 1.  Write code ...............................V
             ! 2.  Compile in debug .........................V
             ! 3.  Runs w/o crashing ........................V
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   OVERLAP (0B) : identically zero
             call set_block_to_zero()

          endif

          if( IsOneBody .and. switch (1,2,2) )then
             !
             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing .................................... ?   (Hamiltonian not tested)
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   1-BODY Generic
             !   \langle A,p; S\Sigma | \hat{O} |B,i; S\Sigma\rangle = 
             !   {\sum_r}'B^{BA}_{pr}o_{ri} + \delta_{AB} o_{pi}
             !   Hamiltonian (1B) is just a special case: o_{pi} = \tilde{h}_{pi}
             aic_vic_M : BLOCK
               type(ClassMatrix) :: B, Op
               type(ClassIrrep)  :: IrrOrbMid !Irrep of the middle orbital
               integer           :: nInctvOrbMid, iIrrOrbMid
               !.. Get N inactive for the three relevant orbital symmetries
               IrrOrbMid = OpIrrep * IrrOrbKet
               iIrrOrbMid = GlobalGroup%GetIrrepIndex( IrrOrbMid )
               !nInctvOrbBra = nInactive( GlobalGroup%GetIrrepIndex( IrrOrbBra ) )
               nInctvOrbMid = nInactive( GlobalGroup%GetIrrepIndex( IrrOrbMid ) )
               B = TDM_Manager%GetB( IonBraLabel, IonKetLabel, IrrOrbBra )
               Op= GlobalIntegral%Get1B( IntOpLabel,"MO_HY", iIrrOrbMid, iIrrOrbKet )
               call Op%Drop("rows",1,nInctvOrbMid) 
               call B%Multiply(Op,"Right","N")
               self%block=B
               if(IonBraLabel .is. IonKetLabel)then
                  call self%Block%Add(Op)
               endif
             end BLOCK aic_vic_M

          endif

          if( IsTwoBody .and. switch (1,2,3) )then
             !
             ! 1.  Write code ............................................V
             ! 2.  Compile in debug ......................................V
             ! 3.  Runs w/o crashing 
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   HAMILTONIAN (2B)
             !   \langle A,p; S\Sigma | \hat{H} |B,i; S\Sigma\rangle =
             !   {\sum_{rs}}' \left( \delta_{S_AS_B} [pi|rs] A^{BA}_{sr} +
             !                                       [ps|ri] B^{BA}_{sr} \right) +
             !   {\sum_{rst}}' [it|rs] Q^{BA}_{sp,rt}
             aic_vic_H : BLOCK
               !.. Here, we must ADD the 2B result to the one computed in the
               !   monoelectronic section.
               type(ClassMatrix)  , pointer  :: A, B
               type(ClassMatrix)  , pointer  :: C
               type(ClassMatrix4D), pointer  :: Biel, Q
               type(ClassIrrep)              :: Irr1
               character(len=:), allocatable :: IntLabel
               integer                       :: ivecmi(4), ivecma(4), i, nr, nc
               integer                       :: iIrr1, iIrr2, iIrr3, iIrr4
               integer                       :: i4miB, i4maB

               call xic_vic_B_Case12()

               call self%GetSize(nr,nc)

               !.. PART C:  Q^{BA}_{sr,pt} [sr|ti]
               ivecmi(3) = nInactive( iIrrOrbBra )+1
               ivecma(3) = ivecmi(3) - 1 + nActive( iIrrOrbBra )
               IntLabel = self%GetIntegralLabel2B("12")
               if(IntLabel(5:5).is."L")then
                  i4miB = nInactive( iIrrOrbKet )+1
                  i4maB = i4miB + nActive( iIrrOrbKet )
               else
                  i4miB = 1
                  i4maB = nc
               endif
               do iIrr4 = 1, size(vIrreps)
                  ivecmi(4) = nInactive( iIrr4 ) + 1
                  ivecma(4) = ivecmi(4) - 1 + nActive(   iIrr4 ) 
                  do iIrr2 = 1, size(vIrreps)
                     Irr1  = vIrreps(iIrr2) * vIrreps(iIrr4) * IrrOrbKet
                     iIrr1 = GlobalGroup%GetIrrepIndex(Irr1)
                     ivecmi(1) = nInactive( iIrr1 ) + 1
                     ivecma(1) = ivecmi(1) - 1 + nActive(   iIrr1 ) 
                     ivecmi(2) = nInactive( iIrr2 ) + 1
                     ivecma(2) = ivecmi(2) - 1 + nActive(   iIrr2 ) 
                     Q    => TDM_Manager%GetQ( IonBraLabel, IonKetLabel, iIrr1, iIrr2, iIrrOrbBra, S2 )
                     Biel => GlobalIntegral%GetBiel( IntLabel, iIrr1, iIrr2, iIrr4, iIrrOrbKet )
                     C    => Q%Contract_124_123(Biel,ivecmi,ivecma,i4miB,i4maB)
                     call self%Block%Add(C)
                  enddo
               enddo

             end BLOCK aic_vic_H

          endif

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(1,2)=TimeMatrix(1,2)+etime

       end select

    type is ( ClassVirtualInternalChannel )

       select type( pK => self%IonChKet )
       type is ( ClassActiveInternalChannel )

          call ElapsedTime()

          ! 1.  Write code ...........................................V
          ! 2.  Compile in debug .....................................V
          ! 3.  Runs w/o crashing ...?
          ! 4.  Unit testing
          !     4.a Dimensions are right
          !     4.b Right elements are zero/non-zero. 
          !     4.c Repres. matrix elem. are correct
          !     4.d Results match STEX benchmark
          ! 5.  Integration with other blocks
          !     5.a 4.a-c consistent upon integration
          !     5.b Results match STEX benchmark
          !
          !.. VIC-AIC  (CASE 2')
          call self%Transpose()
          call self%Build(bstat, switch)
          call self%Transpose()
          if(self%OpLabel(1:1)=="D") call self%Block%Multiply(-1.d0)

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(2,1)=TimeMatrix(2,1)+etime

       type is ( ClassVirtualInternalChannel )

          call ElapsedTime()

          !.. VIC-VIC  (CASE 3)

          if( IsOverlap .and. switch (2,2,1) ) then

             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing ....................................V
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   OVERLAP (0B)
             !   \langle A,i; S\Sigma | B,j; S\Sigma\rangle = \delta_{A B}\delta_{ij}
             vic_vic_S : BLOCK
               integer :: nc, nr
               call self%GetSize(nr,nc)
               call self%block%InitFull(nr,nc)
               self%block=0.d0
               if(IonBraLabel .is. IonKetLabel)call self%Block%AddDiagonal(1.d0)
             end BLOCK vic_vic_S
          endif


          if( IsOneBody .and. switch (2,2,2) )then
             call vxc_vxc_M()
          endif

          if( IsTwoBody .and. switch (2,2,3) )then

             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing ...?
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   HAMILTONIAN (2B)
             !   \langle A,i; S\Sigma | \hat{H} |B,j; S\Sigma\rangle =
             !   {\sum_{rs}}' \left( \delta_{S_AS_B} [ij|rs] A^{BA}_{sr} +
             !                                       [is|rj] B^{BA}_{sr} \right)
             vic_vic_H : BLOCK
               call xic_vic_B_Case12()
             end BLOCK vic_vic_H

          end if

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(2,2)=TimeMatrix(2,2)+etime

       type is ( ClassVirtualExternalChannel )

          call ElapsedTime()

          !.. VIC-VEC  (CASE 4)
          if( IsOverlap .and. switch (3,2,1) )then
             !
             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing ....................................V
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   OVERLAP (0B)
             !   \langle A,e; S\Sigma | B,e'; S\Sigma\rangle = \delta_{A B}S_{ee'}.
             vic_vec_S : BLOCK
               character(len=:), allocatable :: IntLabel
               integer                       :: nc, nr
               IntLabel=self%GetIntegralLabel1B()
               if( ( IonBraLabel .is. IonKetLabel ) .and. ( IntLabel .is. "HY_BS" ) )then
                  !***
                  !*** THE SECOND INDEX IS SUPPOSED TO BE THE LM PAIR INDEX OF THE KET
                  !***
                  !    CALL GLOBALINTEGRAL%lmPairIndex(..)
                  !    EVEN BETTER, GLOBALINTEGRAL MAY FIGURE OUT THE PAIR INDEX BY ITSELF
                  !    FROM THE VALUE OF L AND M
                  self%block = GlobalIntegral%Get1B( IntOpLabel,IntLabel, iIrrOrbBra, iIrrOrbKet )
               else
                  call self%GetSize(nr,nc)
                  call self%block%InitFull(nr,nc, 0,0)
               endif
             end BLOCK vic_vec_S
          endif

          if( IsOneBody .and. switch (3,2,2) )then
             call vxc_vxc_M()
          endif

          if( IsTwoBody .and. switch (3,2,3) )then
             call vxc_vec_2B()
          endif

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(2,3)=TimeMatrix(2,3)+etime

       end select

    type is ( ClassVirtualExternalChannel )
       select type( pK => self%IonChKet )
       type is ( ClassVirtualInternalChannel )

          !.. VEC-VIC   (CASE 4')
          call ElapsedTime()

          ! 1.  Write code ...........................................V
          ! 2.  Compile in debug .....................................V
          ! 3.  Runs w/o crashing ...?
          ! 4.  Unit testing
          !     4.a Dimensions are right
          !     4.b Right elements are zero/non-zero. 
          !     4.c Repres. matrix elem. are correct
          !     4.d Results match STEX benchmark
          ! 5.  Integration with other blocks
          !     5.a 4.a-c consistent upon integration
          !     5.b Results match STEX benchmark
          !
          call self%Transpose()
          call self%Build(bstat, switch)
          call self%Transpose()
          if(self%OpLabel(1:1)=="D") call self%Block%Multiply(-1.d0)

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(3,2)=TimeMatrix(3,2)+etime

       type is ( ClassVirtualExternalChannel )

          !.. VEC-VEC   (CASE 5)
          call ElapsedTime()

          if( IsOverlap .and. switch (3,3,1) )then
             ! 1.  Write code ...........................................V
             ! 2.  Compile in debug .....................................V
             ! 3.  Runs w/o crashing ....................................V
             ! 4.  Unit testing
             !     4.a Dimensions are right
             !     4.b Right elements are zero/non-zero. 
             !     4.c Repres. matrix elem. are correct
             !     4.d Results match STEX benchmark
             ! 5.  Integration with other blocks
             !     5.a 4.a-c consistent upon integration
             !     5.b Results match STEX benchmark
             !
             !   OVERLAP (0B)
             !   \langle A,e; S\Sigma | B,e'; S\Sigma\rangle = \delta_{A B}S_{ee'}.
             vec_vec_S : BLOCK
               character(len=:), allocatable :: IntLabel
               integer                       :: nc, nr
               IntLabel=self%GetIntegralLabel1B()
               if( IonBraLabel .is. IonKetLabel )then
                  if( IntLabel .is. "BS_BS" ) &
                       self%block = GlobalIntegral%Get1B( IntOpLabel,IntLabel, iIrrOrbBra, iIrrOrbKet )
               else
                  call self%GetSize(nr,nc)
                  call self%block%InitFull(nr,nc, 0,0)
               endif
             end BLOCK vec_vec_S
          endif

          if(IsOneBody .and. switch (3,3,2) )then
             call vxc_vxc_M()
          endif

          if(IsTwoBody .and. switch (3,3,3) )then
             call vxc_vec_2B()
          endif

          call ElapsedTime(etime)
          if(present(TimeMatrix)) TimeMatrix(3,3)=TimeMatrix(3,3)+etime

       end select
    end select
    bstat = 0
    !
  contains

    
    subroutine xic_vic_B_Case12()
      ! 1.  Write code ...........................................V
      ! 2.  Compile in debug .....................................V
      ! 3.  Runs w/o crashing ...?
      ! 4.  Unit testing
      !     4.a Dimensions are right
      !     4.b Right elements are zero/non-zero. 
      !     4.c Repres. matrix elem. are correct
      !     4.d Results match STEX benchmark
      ! 5.  Integration with other blocks
      !     5.a 4.a-c consistent upon integration
      !     5.b Results match STEX benchmark
      !
      !   HAMILTONIAN (2B)
      !   \langle A,i; S\Sigma | \hat{H} |B,j; S\Sigma\rangle_PART1 =
      !   {\sum_{rs}}' \delta_{S_AS_B} [ij|rs] A^{BA}_{sr}
      !
      type(ClassMatrix)  , pointer  :: A, B
      type(ClassMatrix)  , pointer  :: C
      type(ClassMatrix4D), pointer  :: Biel
      type(ClassIrrep)              :: Irr1
      character(len=:), allocatable :: IntLabel
      integer                       :: ivecmi(4), ivecma(4), nr, nc
      integer                       :: iIrr1, iIrr2, iIrr3, iIrr4

      call self%GetSize(nr,nc)

      !PART 1: A^BA_sr [sr|pi]
      if(SA2 == SB2)then
         IntLabel = self%GetIntegralLabel2B("12")
         if(IntLabel(4:4).is."L")then
            ivecmi(3) = nInactive( iIrrOrbBra )+1
            ivecma(3) = ivecmi(3) - 1 + nActive( iIrrOrbBra )
         else
            ivecmi(3) = 1
            ivecma(3) = nr
         endif
         if(IntLabel(5:5).is."L")then
            ivecmi(4) = nInactive( iIrrOrbKet )+1
            ivecma(4) = ivecmi(4) - 1 + nActive( iIrrOrbKet )
         else
            ivecmi(4) = 1
            ivecma(4) = nc
         endif
         do iIrr2 = 1, size(vIrreps)
            Irr1  = IrrAB*vIrreps(iIrr2)
            iIrr1 = GlobalGroup%GetIrrepIndex(Irr1)
            ivecmi(1) = nInactive( iIrr1 ) + 1
            ivecma(1) = ivecmi(1) - 1 + nActive( iIrr1 ) 
            ivecmi(2) = nInactive( iIrr2 ) + 1
            ivecma(2) = ivecmi(2) - 1 + nActive( iIrr2 )
            A    => TDM_Manager%GetA( IonBraLabel, IonKetLabel, Irr1 )
            Biel => GlobalIntegral%GetBiel( IntLabel, iIrr1, iIrr2, iIrrOrbBra, iIrrOrbKet  )
            C    => A%Contract4D12(Biel,ivecmi,ivecma)
            call self%Block%Add(C)
         enddo
      end if

      !.. PART B:  B^{BA}_sr [sp|ri]
      IntLabel = self%GetIntegralLabel2B("13")
      if(IntLabel(2:2).is."L")then
         ivecmi(2) = nInactive( iIrrOrbBra ) + 1
         ivecma(2) = ivecmi(2) - 1 + nActive( iIrrOrbBra )
      else
         ivecmi(2) = 1
         ivecma(2) = nr
      endif
      if(IntLabel(5:5).is."L")then
         ivecmi(4) = nInactive( iIrrOrbKet )+1
         ivecma(4) = ivecmi(4) - 1 + nActive( iIrrOrbKet )
      else
         ivecmi(4) = 1
         ivecma(4) = nc
      endif
      do iIrr3 = 1, size(vIrreps)
         Irr1  = IrrAB*vIrreps(iIrr3)
         iIrr1 = GlobalGroup%GetIrrepIndex(Irr1)
         ivecmi(1) = nInactive( iIrr1 ) + 1
         ivecma(1) = ivecmi(1) - 1 + nActive(   iIrr1 ) 
         ivecmi(3) = nInactive( iIrr3 ) + 1
         ivecma(3) = ivecmi(3) - 1 + nActive(   iIrr3 )
         ! write(*,*)
         ! write(*,*) "nInactive",IntLabel
         ! write(*,*) nInactive
         ! write(*,*) "--------nActive"
         ! write(*,*) nActive
         
         B    => TDM_Manager%GetB( IonBraLabel, IonKetLabel, Irr1 )
         ! write(*,*) "yyyy",B%Nrows(),B%Ncolumns() 
         Biel => GlobalIntegral%GetBiel( IntLabel, iIrr1, iIrrOrbBra, iIrr3, iIrrOrbKet  )
         C    => B%Contract4D13(Biel,ivecmi,ivecma)
         call self%Block%Add(C)
      enddo

    end subroutine xic_vic_B_Case12
    
    subroutine vxc_vxc_M()
      ! 1.  Write code ...........................................V
      ! 2.  Compile in debug .....................................V
      ! 3.  Runs w/o crashing ...?
      ! 4.  Unit testing
      !     4.a Dimensions are right
      !     4.b Right elements are zero/non-zero. 
      !     4.c Repres. matrix elem. are correct
      !     4.d Results match STEX benchmark
      ! 5.  Integration with other blocks
      !     5.a 4.a-c consistent upon integration
      !     5.b Results match STEX benchmark
      !
      !   1-BODY
      !   \langle A,i; S\Sigma | \hat{O} |B,j; S\Sigma\rangle = 
      !   \delta_{S_A S_B}\delta_{ij}O_{AB} +\delta_{AB} o_{ij}
      !
      !.. The One-body hamiltonian has the additional peculiarity
      !   that O_{AB} is diagonal, and o_{ij} is the h-tilde:
      !   \delta_{A B} ( \delta_{ij} E_A + \tilde{h}_{ij} )
      type(ClassMatrix)             :: Op
      real(kind(1d0))               :: O_AB
      character(len=:), allocatable :: IntLabel
      integer                       :: ndropBra, ndropKet!, iIrrBra, iIrrKet

      IntLabel=self%GetIntegralLabel1B()
      ndropBra = nInactive( iIrrOrbBra ) + nActive( iIrrOrbBra )
      ndropKet = nInactive( iIrrOrbKet ) + nActive( iIrrOrbKet )
      if(IonBraLabel .is. IonKetLabel)then
         Op= GlobalIntegral%Get1B( IntOpLabel,IntLabel, iIrrOrbBra, iIrrOrbKet )
         if(IntLabel(1:2).is."MO") call Op%Drop("rows"   ,1,ndropBra) 
         if(IntLabel(4:5).is."MO") call Op%Drop("columns",1,ndropKet) 
         self%block=Op
      endif

      O_AB=TDM_Manager%GetMatEl(IntOpLabel, IonBraLabel, IonKetLabel)
      if( abs(O_AB) > 0.d0 )then

         Op= GlobalIntegral%Get1B( OverlapLabel,IntLabel,iIrrOrbBra, iIrrOrbKet )
         if(IntLabel(1:2).is."MO") call Op%Drop("rows"   ,1,ndropBra) 
         if(IntLabel(4:5).is."MO") call Op%Drop("columns",1,ndropKet) 
         call Op%Multiply(O_AB)
         call self%Block%Add(Op)
      endif
    end subroutine vxc_vxc_M

    subroutine set_block_to_zero()
      integer :: nr, nc
      call self%GetSize(nr,nc)
      call self%block%InitFull(nr,nc,0,0)
    end subroutine set_block_to_zero

    subroutine vxc_vec_2B()
      !
      ! 1.  Write code ...........................................V
      ! 2.  Compile in debug .....................................V
      ! 3.  Runs w/o crashing ...?
      ! 4.  Unit testing
      !     4.a Dimensions are right
      !     4.b Right elements are zero/non-zero. 
      !     4.c Repres. matrix elem. are correct
      !     4.d Results match STEX benchmark
      ! 5.  Integration with other blocks
      !     5.a 4.a-c consistent upon integration
      !     5.b Results match STEX benchmark
      !
      !   HAMILTONIAN (2B)
      !   \langle A,p; S\Sigma | \hat{H} |B,q; S\Sigma\rangle =
      !   \sum_{\ell m}^{\ell>0} o^{(\ell m)}_{pq} M^{(\ell m)}_{AB}
      type(ClassMatrix)             :: Op
      real(kind(1d0))               :: M_AB
      character(len=:), allocatable :: IntLabel
      integer :: lket, mket, lmin, lmax, mmax, lbra, mbra, l, m, ilm
      integer :: indexBra, indexKet
      IntLabel=self%GetIntegralLabel1B()
      if(IntLabel(1:2).is."MO")then
         call set_block_to_zero()
         return
      endif
      select type ( ptr => self%IonChKet )
      type is ( ClassVirtualExternalChannel )
         lket     = ptr%GetL()
         mket     = ptr%GetM()
         indexKet = GlobalIntegral%lmPairIndex(lket,mket)
      class default
         return
      end select
      select type ( ptr => self%IonChBra )
      type is ( ClassVirtualInternalChannel )
         lmin   = 0
         lmax     = ptr%GetLmax() + lket
         mmax     = ptr%GetLmax() + abs(mket)
         indexBra = iIrrOrbBra
      type is ( ClassVirtualExternalChannel )
         lbra     = ptr%GetL()
         mbra     = ptr%GetM()
         indexBra = GlobalIntegral%lmPairIndex(lbra,mbra)
         lmin     = abs(lbra-lket)
         lmax     = lbra+lket
         mmax     = abs(mket)+abs(mbra)
      end select
      do l = lmin, lmax 
         do m = -mmax, mmax
            ilm = l**2 + l + m + 1
            M_AB=TDM_Manager%GetMatEl(IntOpLabel, IonBraLabel, IonKetLabel, ilm)
            if( abs(M_AB) > 0.d0 )then
               Op= GlobalIntegral%Get1B( OverlapLabel,IntLabel, indexBra, indexKet, ilm )
               call Op%Multiply(M_AB)
               call self%Block%Add(Op)
            endif
         enddo
      enddo
    end subroutine vxc_vec_2B


  end subroutine ClassCCBlockBuild


  function ClassCCBlockGetBlockName( self ) result(BlockName)
    class(ClassCCBlock)  , intent(in) :: self
    character(len=:), allocatable          :: BlockName
    allocate( BlockName, source = &
         AddSlash(self%GetStorageDir())//GetccExt(self%IonChBra)//"_"//GetccExt(self%IonChKet))
  end function ClassCCBlockGetBlockName

  function GetccExt(ic) result(Ext)
    class(ClassCloseCouplingChannel), intent(in) :: ic
    character(len=:), allocatable :: Ext
    select type ( ptr => ic )
    type is (ClassActiveInternalChannel)
       Ext=ptr%GetLabelExt()
    type is (ClassVirtualInternalChannel)
       Ext=ptr%GetLabelExt()
    type is (ClassVirtualExternalChannel)
       Ext=ptr%GetLabelExt()
    end select
  end function GetccExt

  function ClassCCBlockGetIntegralLabel1B(self) result( label )
    class(ClassCCBlock), intent(in) :: self
    character(len=5) :: label
    character(len=2) :: sIntBra, sIntKet
    character(len=:), allocatable :: labelBra, labelKet
    labelBra= GetccExt(self%IonChBra)
    if( labelBra .is. viMOc_ID )then
       sIntBra = "MO"
    elseif( labelBra .is. hiGOc_ID )then
       sIntBra = "HY"
       !*** Juan made a change here
    elseif( index(labelBra,".") .gt. 0 )then
    !elseif( labelBra .is. veSOc_ID )then
       sIntBra = "BS"
    endif
    labelKet= GetccExt(self%IonChKet)
    if( labelKet .is. viMOc_ID )then
       sIntKet = "MO"
    elseif( labelKet .is. hiGOc_ID )then
       sIntKet = "HY"
       !*** Juan made a change here
    elseif( index(labelKet,".") .gt. 0 )then
    !elseif( labelKet .is. veSOc_ID )then
       sIntKet = "BS"
    endif
    label = sIntBra//"_"//sIntKet
  end function ClassCCBlockGetIntegralLabel1B

  function ClassCCBlockGetIntegralLabel2B(self,case) result( label )
    class(ClassCCBlock), intent(in) :: self
    character(len=*)   , intent(in) :: case
    character(len=5) :: label
    character        :: sIntBra,sIntKet
    character(len=:), allocatable :: labelBra,labelKet
    labelBra= GetccExt(self%IonChBra)
    if( labelBra .is. viMOc_ID )then
       sIntBra = "L"
    elseif( labelBra .is. hiGOc_ID )then
       sIntBra = "S"
    endif
    labelKet= GetccExt(self%IonChKet)
    if( labelKet .is. viMOc_ID )then
       sIntKet = "L"
    elseif( labelKet .is. hiGOc_ID )then
       sIntKet = "S"
    endif
    if(case.is."12")label = "LL_"//sIntBra//sIntKet
    if(case.is."13")label = "L"//sIntBra//"_L"//sIntKet
  end function ClassCCBlockGetIntegralLabel2B

  
  integer function ClassCCBlockGetNRows( self ) result( N )
    class(ClassCCBlock), intent(in) :: self
    N = self%Block%NRows()
  end function ClassCCBlockGetNRows

  integer function ClassCCBlockGetNColumns( self ) result( N )
    class(ClassCCBlock), intent(in) :: self
    N = self%Block%NColumns()
  end function ClassCCBlockGetNColumns

  subroutine ClassCCBlockUnsetFormattedWrite( self )
    class(ClassCCBlock), intent(inout) :: self
    self%FormattedWrite = .false.
  end subroutine ClassCCBlockUnsetFormattedWrite

  subroutine ClassCCBlockSetFormattedWrite( self )
    class(ClassCCBlock), intent(inout) :: self
    self%FormattedWrite = .true.
  end subroutine ClassCCBlockSetFormattedWrite

  function ClassCCBlockFetchBlock( self ) result(NewBlock)
    class(ClassCCBlock), intent(in) :: self
    type(ClassMatrix) :: NewBlock
    NewBlock = self%Block
  end function ClassCCBlockFetchBlock

  subroutine ClassCCBlockFetchMatrix( self, Mat )
    class(ClassCCBlock), intent(in) :: self
    real(kind(1d0)) :: Mat(:,:)
    Mat = self%Block%A
  end subroutine ClassCCBlockFetchMatrix

  subroutine ClassCCBlockFetchMatrixCM( self, Mat )
    class(ClassCCBlock), intent(in)  :: self
    type(ClassMatrix)       , intent(out) :: Mat
    Mat = self%Block
  end subroutine ClassCCBlockFetchMatrixCM


  !-----------------------------------------------------------------------
  ! Symmetric Ionic Channel
  subroutine ChanContainerFinal( self )
    type(ChanContainer) :: self
    if(associated(self%Ch))then
       call self%Ch%Free()
    endif
  end subroutine ChanContainerFinal

  subroutine ClassSymIonCh_ParseConfLine( self, Strn )
    !
    use moduleParameterList
    use moduleString
    class(ClassSymIonCh), intent(inout) :: self

    !> String containing the specification of the parent ion and 
    !! of the Xlm waves that must be associated with it
    !! e.g.: 1Ag.1 ( aiM viM hiG beS: 0  0, 2  0, 2  2 )
    character(len=*)                    , intent(in)    :: Strn
    !
    integer                       :: ichar, dichar, ichar1, ichar2, iXlm
    character(len=:), allocatable :: PIonLabel, LMvecLabel, chLabel
    type(ClassGroup), pointer     :: Group
    integer                       :: NumLM, NAIMOS, NVIRTINTCH
    integer, allocatable          :: LMvector(:,:), ivBuf(:)
    type(ClassXlm)                :: Xlm
    character(len=8)              :: PionIrrep
    character(len=32)             :: lString, mString, LmaxString
    integer                       :: PionMult, PionN, iCh
    type(ClassIrrep), pointer     :: OrbitalIrrep
    type(ClassXlmSymmetricSet)    :: XlmSymSet
    logical                       :: UseAllXlm = .false.
    logical                       :: use_viMOc, use_hiGOc
    character(len=:), allocatable :: vic_type

    !.. 1. Extract parent-ion parameters
    ichar1 = index(Strn,'(')
    if ( ichar1 == 0 ) call Assert( "Invalid channel string." )
    allocate( PionLabel, source = trim(adjustl(Strn(:ichar1-1))) )

    !.. 1.a Check if the spin of the parent ion is compatible with the
    !       total spin of the space
    call ParsePionLabel( PionLabel, PionMult, PionIrrep, PionN )
    call CheckSpin( Self%Multiplicity, PionMult )

    !.. 2. Initialize the parent ion
    Group => Self%Irrep%GetGroup()
    call Self%ParentIon%Init( Group, PionLabel, Self%Charge + 1)

    OrbitalIrrep => Self%Irrep * Self%ParentIon%GetIrrep()
    call XlmSymSet%Init( Self%Lmax, OrbitalIrrep )

    !.. 1.b Initialize the name of the storage directory
    if(allocated(self%StorageDir))deallocate(self%StorageDir)
    allocate(self%StorageDir,source=self%Root//"/"//&
         self%ParentIon%GetLabel()//"_x_"//&
         OrbitalIrrep%GetName())

    !.. 3. Extract the list of channels
    ichar2 = index(Strn,')')
    if ( ichar2 == 0 ) call Assert( "Invalind channel string." )
    allocate( LMvecLabel, source = trim(adjustl(Strn(ichar1+1:ichar2-1))) )
    !.. 4. Search for the specification of the active internal MO channel
    NAIMOS=0
    chLabel = Capitalize(aiMOc_ID)
    ichar1  = index(LMvecLabel,chLabel)
    if(ichar1>0)then
       NAIMOS=NAIMOS+1
       LMvecLabel(ichar1:ichar1+len(chLabel)-1)=" "
    endif
    self%NAIMOS=NAIMOS

    !.. 5. Search for the specification of the virtual internal MO channel
    NVIRTINTCH=0
    chLabel = Capitalize(viMOc_ID)
    ichar1  = index(LMvecLabel,chLabel)
    use_viMOc = (ichar1>0)
    if(use_viMOc)then
       NVIRTINTCH=NVIRTINTCH+1
       LMvecLabel(ichar1:ichar1+len(chLabel)-1)=" "
    endif

    !.. Hybrid short-range channel
    chLabel = Capitalize(hiGOc_ID)
    ichar1  = index(LMvecLabel,chLabel)
    use_hiGOc = (ichar1>0)
    if(use_hiGOc)then
       NVIRTINTCH=NVIRTINTCH+1
       LMvecLabel(ichar1:ichar1+len(chLabel)-1)=" "
    endif
    self%NVIRTINTCH=NVIRTINTCH

    LMvecLabel=adjustl(LMvecLabel)

    !.. 6. Search for the external virtual functions
    chLabel = Capitalize( veSOc_ID )
    ichar1=index(LMvecLabel,chLabel)
    NumLM=0
    UseAllXlm=.FALSE.
    if( ichar1 > 0 )then

       LMvecLabel=adjustl(LMvecLabel(ichar1+len(chLabel):))
       ichar1=index(LMvecLabel,":")
       if(ichar1<0)then
          call ErrorMessage("Syntax Error in cc cfg file at "//trim(LMvecLabel))
          stop
       endif
       LMvecLabel=adjustl(LMvecLabel(ichar1+1:))


       UseAllXlm = ( index(LMvecLabel,USE_ALL_XLMSYM) > 0 )
       !.. Check if the specified Xlm are compatible with the current symmetry
       if ( .not. UseAllXlm ) then

          !.. Count number of specified Xlm channels 
          NumLM = 1
          ichar = 0
          XlmCycle : do

             dichar = index(LMvecLabel(ichar+1:),',')
             if(dichar<1)exit XlmCycle
             ichar = ichar + dichar
             LMvecLabel(ichar:ichar) = " "
             NumLM = NumLM + 1

          end do XlmCycle
          self%NVIRTEXTCH=NumLM

          !.. Parse list of specified Xlm channels
          allocate(ivBuf(2*NumLM))
          read(LMvecLabel,*) (ivBuf(iXlm),iXlm=1,2*NumLM) 
          allocate( LMVector( 2, NumLM ) )
          do iXlm=1,NumLM
             LMVector(1,iXlm)=ivBuf(2*iXlm-1)
             LMVector(2,iXlm)=ivBuf(2*iXlm)
          enddo
          deallocate(ivBuf)

       endif

    endif

    if( NumLM + NAIMOS + NVirtIntCh == 0 )then
       call ErrorMessage( "No channels found; nothing to do." )
       stop
    endif

    if( UseAllXlm )then
       NumLM = XlmSymSet%GetNXlm()
       Self%NVIRTEXTCH = NumLM
    endif
    allocate( Self%Chanv(NumLM + NAIMOS + NVIRTINTCH) )
    Self%NChans = size(self%chanv)

    do iXlm = 1, NumLM

       if( UseAllXlm )then
          Xlm = XlmSymSet%GetXlm(iXlm)
       else
          call Xlm%Init( LMvector(1,iXlm), LMvector(2,iXlm) )
          if ( .not. XlmSymSet%ValidXlm( Xlm ) ) then
             write(lString,*) Xlm%GetL()
             write(mString,*) Xlm%GetM()
             write(LmaxString,*) Self%Lmax
             call ErrorMessage( "The symmetry adapted spherical harmonic with l: "//&
                  trim(adjustl(lString))//" and m: "//trim(adjustl(mString))//&
                  " do not belong to the irreducible representation: "//&
                  OrbitalIrrep%GetName()//" with maximum angular momentum: "//&
                  trim(adjustl(LmaxString)) )
             stop
          end if
       endif

       Self%Chanv(NAIMOS+NVIRTINTCH+iXlm)%Ch => &
            virtual_external_channel_farm( &
            Self%ParentIon        , &
            Self%Multiplicity, &
            Self%Irrep       , &
            Self%StorageDir          , &
            Xlm                   )

    enddo

    do iCh = 1, NAIMOS
       Self%Chanv(iCh)%Ch => &
            Active_Internal_Channel_Farm(&
            Self%ParentIon        , &
            Self%Multiplicity, &
            Self%Irrep       , &
            Self%StorageDir          )
    enddo

    do iCh = 1, NVIRTINTCH
       if(use_viMOc.and.iCh==1)then
          vic_type = viMOc_ID
       else
          vic_type = hiGOc_ID
       endif
       Self%Chanv(NAIMOS+iCh)%Ch => &
            Virtual_Internal_Channel_Farm( &
            Self%ParentIon        , &
            vic_type              , &
            Self%Lmax             , &
            Self%Multiplicity, &
            Self%Irrep       , &
            Self%StorageDir          )
    enddo

  contains
    !
    subroutine CheckSpin( TotalMult, PionMult )
      integer, intent(in) :: TotalMult, PionMult
      if ( abs( TotalMult - PionMult ) /= 1 )then 
         write(ERROR_UNIT,"(a,i0)")"Total multiplicity =", TotalMult
         write(ERROR_UNIT,"(a,i0)")"P-ion multiplicity =", PionMult
         call Assert("Total spin and parent-ion spin are incompatible")
      endif
    end subroutine CheckSpin
    !
  end subroutine ClassSymIonCh_ParseConfLine

  subroutine ClassSymIonCh_Show( self, unit )
    class(ClassSymIonCh), intent(in) :: self
    integer, optional             ,         intent(in) :: unit
    integer :: outunit, iChan
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    call self%ClassSymESpace%show(unit)
    call self%ParentIon%Show()
    do iChan = 1, self%NChans
       call self%Chanv(iChan)%Ch%Show(unit)
    enddo
  end subroutine ClassSymIonCh_Show

  subroutine ClassSymIonCh_Free( self )
    class(ClassSymIonCh), intent(inout) :: self
    integer :: iCh
    call self%ParentIon%Free()
    self%NAIMOS=0
    self%NVIRTINTCH=0
    self%NVIRTEXTCH=0
    if ( allocated(self%Chanv) )then
       do iCh=1,self%NChans
          if(associated(self%Chanv(iCh)%Ch))then
             call self%Chanv(iCh)%Ch%Free()
          endif
       enddo
       deallocate( self%Chanv )
    endif
    self%NChans = 0
    !
    call self%ClassSymESpace%Free()
  end subroutine ClassSymIonCh_Free

  subroutine ClassSymIonCh_Final( self )
    type(ClassSymIonCh) :: self
    call self%free()
  end subroutine ClassSymIonCh_Final


  subroutine ClassSymIonCh_Init( self, &
       Lmax, Multiplicity, Irrep, Charge, Root )
    class(ClassSymIonCh)     , intent(inout) :: self
    integer                  , intent(in)    :: Lmax
    integer                  , intent(in)    :: Multiplicity
    class(ClassIrrep), target, intent(in)    :: Irrep
    integer                  , intent(in)    :: Charge
    character(len=*)         , intent(in)    :: Root
    self%Lmax         =  Lmax
    self%Multiplicity =  Multiplicity
    self%Irrep        => Irrep
    self%Charge       =  Charge
    self%Root         =  Root 
  end subroutine ClassSymIonCh_Init


  character(len=128) function ClassSymIonCh_GetPILabelFun( self ) result( Label )
    class(ClassSymIonCh), intent(in)  :: self
    Label = self%ParentIon%GetLabel() 
  end function ClassSymIonCh_GetPILabelFun

  !.. Return the size associated to a specified kind of channel
  integer function ClassSymIonCh_GetPartialSize( self, RefChan ) result(isize)
    class(ClassSymIonCh), intent(in) :: self
    Class(ClassCloseCouplingChannel)            , intent(in) :: RefChan
    integer :: iCh
    isize=0
    do iCh=1,self%NChans
       if( same_type_as(self%Chanv(iCh)%Ch,RefChan) ) isize = isize + self%Chanv(iCh)%Ch%GetSize()
    enddo
  end function ClassSymIonCh_GetPartialSize

  subroutine ClassSymIonCh_GetNChans( self, NAIMOS, NVIRTINTCH, NVIRTEXTCH, NTOT ) 
    class(ClassSymIonCh), intent(in)  :: self
    integer                               , intent(out) :: NAIMOS, NVIRTINTCH, NVIRTEXTCH, NTOT
    NAIMOS = self%NAIMOS
    NVIRTINTCH = self%NVIRTINTCH
    NVIRTEXTCH = self%NVIRTEXTCH
    NTOT = self%NAIMOS + self%NVIRTINTCH + self%NVIRTEXTCH
  end subroutine ClassSymIonCh_GetNChans

  integer function ClassSymIonCh_GetNumChannels( self ) result(nchans)
    class(ClassSymIonCh), intent(in)  :: self
    nchans = self%NAIMOS + self%NVIRTINTCH + self%NVIRTEXTCH
  end function ClassSymIonCh_GetNumChannels
  
  function ClassSymIonCh_GetChannel( Self, iChan ) result(pChan)
    class(ClassSymIonCh),  intent(in) :: Self
    integer             ,  intent(in) :: iChan
    class(ClassCloseCouplingChannel), pointer :: pChan
    pChan => NULL()
    if(allocated(self%Chanv))then
       if(self%nchans>=iChan) pChan => Self%Chanv(iChan)%Ch
    endif
  end function ClassSymIonCh_GetChannel

  subroutine ClassSymIonCh_Condition( self ) 
    !
    ! 1.  Write code ...........................................V
    ! 2.  Compile in debug .....................................V
    ! 3.  Runs w/o crashing ...?
    ! 4.  Unit testing
    !     4.a Dimensions are right
    !     4.b Right elements are zero/non-zero. 
    !     4.c Repres. matrix elem. are correct
    !     4.d Results make sense (S pos. def., fedvr localized)
    ! 5.  Integration with other blocks
    !     5.a 4.a-c consistent upon integration
    !     5.b Results match STEX benchmark
    !
    class(ClassSymIonCh), intent(inout) :: self
    type(ClassSESSESBlock)              :: SBlk, HBlk, RBlk
    logical                             :: lTasks(N_SESSES_ID)
    character(len=:), allocatable       :: IonLabel

    call SBlk%Free()
    call HBlk%Free()
    call RBlk%Free()

    call SBlk%Init(self,self,"S",self%Root)
    call HBlk%Init(self,self,"H",self%Root) 
    call RBlk%Init(self,self,"R",self%Root) 

    IonLabel = self%ParentIon%GetLabel()
    lTasks=.FALSE.; lTasks(SESSES_LOAD)=.TRUE.


    BLOCK
      !.. Diagonalize the Hamiltonian in the viM (+) hiG spaces.
      !   In absence of vec channels, this diagonalization can
      !   be used to construct H0 in the solution of the interacting-channel
      !   case, as well as for the TDSE propagation. If everything is coded
      !   right, this same diagonalization would be performed anyway later on,
      !   so this step is only for consistency checks.
      !..
      type(ClassVirtualInternalChannel) :: vic
      real(kind(1d0)), allocatable      :: heval(:)
      type(ClassMatrix)                 :: hmat, hevec
      call Execute_Command_Line("mkdir -p "//self%StorageDir//"/vic")
      call HBlk%Driver(lTasks,vic,vic,IonLabel,IonLabel)
      call HBlk%Assemble(HMat,vic,vic,IonLabel,IonLabel)
      call hmat%Write(self%StorageDir//"/viG/H_viG_viG","formatted")
      call hmat%Diagonalize(heval,hevec)
      call hmat%Free()
      call hevec%Write(self%StorageDir//"/viG/H_viG_viG_evec",  "formatted")
      call SaveVector( self%StorageDir//"/viG/H_viG_eval",heval,"formatted")
      call hevec%Free()
      deallocate(heval)
    END BLOCK


    !************ THE PART BELOW BELONGS TO THE INTEGRALS!
!!$    dvr_bsplines: BLOCK
!!$      !.. Determine the FEDVR version of the outer B-splines
!!$      real(kind(1d0)), parameter        :: INT_DVR_NRM1_THR=1.d-12
!!$      type(ClassCCBlock), pointer       :: scc, rcc, hcc
!!$      type(ClassVirtualExternalChannel) :: vec
!!$      real(kind(1d0)), allocatable      :: reval(:), mrevec(:,:), rgrid(:), heval(:)
!!$      type(ClassMatrix)                 :: smat, rmat, hmat, revec, intdvr, extdvr, hevec
!!$      integer                           :: ich, idvr, bw, n, ndvr_min
!!$      character(len=:), allocatable     :: dir
!!$      logical :: vec_is_present = .TRUE.
!!$
!!$      do ich = 1, self%GetNumChannels()
!!$         scc => SBlk%GetBlock(ich,ich)
!!$         vec_is_present = scc%IsType( vec, vec )
!!$         if( vec_is_present )exit
!!$      end do
!!$      if(.not. vec_is_present)exit dvr_bsplines
!!$      rcc => RBlk%GetBlock(ich,ich)
!!$
!!$      dir = self%StorageDir//"/vec"
!!$      call Execute_Command_Line("mkdir -p "//dir)
!!$
!!$      smat = scc%FetchBlock()
!!$      rmat = rcc%FetchBlock()
!!$      n    = smat%NRows()
!!$      bw   = smat%UpperBandwidth()
!!$
!!$      !.. Drops the last bspline
!!$      call smat%drop("columns",n,n)
!!$      call smat%drop("rows"   ,n,n)
!!$      call rmat%drop("columns",n,n)
!!$      call rmat%drop("rows"   ,n,n)
!!$
!!$      call smat%Write(dir//"/S","formatted")
!!$      call rmat%Write(dir//"/R","formatted")
!!$      call rmat%Diagonalize(smat,reval,revec)
!!$      call rmat%Free()
!!$
!!$      call SaveVector(dir//"/R_eval",reval,"formatted")
!!$      call revec%Write(dir//"R_evec",      "formatted")
!!$
!!$      !.. Determines the set of dvr-like splines that
!!$      !   do not penetrate the internal region
!!$      call revec%FetchMatrix(mrevec)
!!$      ndvr_min = size(reval) + 1
!!$      do 
!!$         if(ndvr_min==1)exit
!!$         if(sum(abs(mrevec(1:bw,ndvr_min-1)))>INT_DVR_NRM1_THR)exit
!!$         ndvr_min=ndvr_min-1
!!$      enddo
!!$      intdvr = mrevec(:,1:ndvr_min-1)
!!$      extdvr = mrevec(:,ndvr_min:)
!!$      allocate(rgrid,source=reval(ndvr_min:))
!!$      call revec%Free()
!!$      deallocate(reval)
!!$
!!$      call SaveVector(dir//"/deSp_R_grid",rgrid,"formatted")
!!$      call extdvr%Write(dir//"/T_beSp_deSp",    "formatted")
!!$      call intdvr%Write(dir//"/T_beSp_diSp",    "formatted")
!!$
!!$      !.. Diagonalize the Hamiltonian in the outer region
!!$      do ich = 1, self%GetNumChannels()
!!$         hcc => HBlk%GetBlock(ich,ich)
!!$         if(.not. hcc%IsType( vec, vec ))cycle
!!$         dir = self%StorageDir//"/vec/"//hcc%GetLabel("bra")
!!$         call Execute_Command_Line("mkdir -p "//dir)
!!$         hmat = hcc%FetchBlock()
!!$         call hmat%drop("columns",n,n)
!!$         call hmat%drop("rows"   ,n,n)
!!$         call hmat%Write(dir//"/H_beSp_beSp","formatted")
!!$         call hmat%Multiply(extdvr,"Right","N")
!!$         call hmat%Multiply(extdvr,"Left" ,"T")
!!$         call hmat%Write(dir//"/H_deSp_deSp","formatted")
!!$         call hmat%Diagonalize(heval,hevec)
!!$         call hevec%Write(dir//"/H_deSp_eeSp_evec","formatted")
!!$         call SaveVector( dir//"/H_eeSp_eval",heval,"formatted")
!!$         call hevec%Multiply(extdvr,"Left","N")
!!$         call hevec%Write(dir//"/H_beSp_eeSp_evec","formatted")
!!$      end do
!!$
!!$      !.. Orthonormalize the penetrating group of fedvr
!!$      !   to the internal region, and 
!!$
!!$    END BLOCK dvr_bsplines


!!$      dir = self%StorageDir//"/vec/"//scc%GetLabel("bra")
!!$         call smat%Diagonalize(seval,tmat)
!!$         call smat%free()
!!$         call SaveVector(dir//"/S_eval",seval,"formatted")
!!$         call tmat%Write(dir//"S_evec","formatted")
!!$         seval = 1.d0 / sqrt(seval)
!!$         call tmat%Multiply(seval,"right")
!!$         call tmat%Write(dir//"T","formatted")
!!$         call hmat%Multiply(tmat,"Right","N")
!!$         call hmat%Multiply(tmat,"Left","T")
!!$         call hmat%Diagonalize(heval,hevec)


  end subroutine ClassSymIonCh_Condition


  !-------------------------------------------------------------------------------
  ! SES - SES Block
  subroutine ClassSESSESBlock_Free( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    Self%SpaceBra => NULL()
    Self%SpaceKet => NULL()
    if ( allocated(Self%mBlock  ) ) deallocate( Self%mBlock  )
    if ( allocated(Self%OpLabel ) ) deallocate( Self%OpLabel )
    if ( allocated(Self%IDLabel ) ) deallocate( Self%IDLabel )
    if ( allocated(Self%StorageDir ) ) deallocate( Self%StorageDir )
    Self%AvailableBlocks    = .false.
    Self%INITIALIZED        = .false.
  end subroutine ClassSESSESBlock_Free
  subroutine ClassSESSESBlock_Final( Self )
    type(ClassSESSESBlock) :: Self
    call Self%Free()
  end subroutine ClassSESSESBlock_Final
  logical function ClassSESSESBlock_IsInitialized( Self ) result( IsInitialized )
    class(ClassSESSESBlock), intent(in) :: Self
    IsInitialized = Self%initialized
  end function ClassSESSESBlock_IsInitialized


  !> Initializes the ClassSESSESBlock class.
  subroutine ClassSESSESBlock_Init( Self, SpaceBra, SpaceKet, OpLabel, Root, IDLabel )
    class(ClassSESSESBlock)        , intent(inout) :: Self
    class(ClassSymESpace), target  , intent(in)    :: SpaceBra
    class(ClassSymESpace), target  , intent(in)    :: SpaceKet
    character(len=*)               , intent(in)    :: OpLabel
    character(len=*)               , intent(in)    :: Root
    character(len=*)     , optional, intent(in)    :: IDLabel
    !
    integer :: NChBra, NChKet
    integer :: i, j
    class(ClassCloseCouplingChannel), pointer :: pChBra, pChKet
    !
    call self%Free()
    !
    Self%SpaceBra => SpaceBra
    Self%SpaceKet => SpaceKet
    allocate(Self%OpLabel,source=trim(adjustl(OpLabel)))
    if(present(IDLabel)) allocate(Self%IDLabel,source=trim(adjustl(IDLabel)))
    allocate(self%StorageDir,source=AddSlash(Root)//&
         SpaceBra%GetLabel()//"_"//SpaceKet%GetLabel()//"/"//self%OpLabel)
    !
    NChBra = Self%SpaceBra%GetNumChannels()
    NChKet = Self%SpaceKet%GetNumChannels()
    allocate( Self%mBlock( NChBra, NChKet ) )
    !
    do j = 1, NChKet
       pChKet => Self%SpaceKet%GetChannel(j)
       do i = 1, NChBra
          pChBra => Self%SpaceBra%GetChannel(i)
          call Self%mBlock(i,j)%Init( pChBra, pChKet, OpLabel, self%StorageDir, IDLabel ) 
       end do
    end do
    !
    Self%INITIALIZED = .true.
    !
  end subroutine ClassSESSESBlock_Init

  !> Initializes the ClassSESSESBlock class.
  function ClassSESSESBlock_GetBlock( Self, i, j ) result( pblk )
    class(ClassSESSESBlock), target, intent(inout) :: Self
    integer                        , intent(in)    :: i,j
    type(ClassCCBlock), pointer :: pblk
    pblk => self%mBlock( i, j )
  end function ClassSESSESBlock_GetBlock

  !> Initializes the ClassSESSESBlock class.
  subroutine ClassSESSESBlock_Driver( Self, TaskList, bChan, kChan, bIonLabel, kIonLabel )
    class(ClassSESSESBlock)                   , intent(inout) :: Self
    logical                                   , intent(in)    :: TaskList(:)
    class(ClassCloseCouplingChannel), optional, intent(in)    :: bChan, kChan
    character(len=*)                , optional, intent(in)    :: bIonLabel,  kIonLabel
    type(ClassCCBlock)      :: blk
    integer                 :: i, j, bstat
    logical                 :: skip_LD
    logical                 :: totally_symmetric_operator
    !
    if(.not.Self%Initialized) call Assert("SESSESBlock not initialized at Driver")
    totally_symmetric_operator = self%SpaceBra%GetLabel() .is. self%SpaceKet%GetLabel()
    do j = 1, Self%SpaceKet%GetNumChannels()

       if(present(kChan))then
          if(.not.self%mBlock( 1, j )%IsType(bChan,kChan))cycle
       endif
       if(present(kIonLabel))then
          if( self%mBlock( 1, j )%GetIonKetLabel() .isnt. kIonLabel )cycle
       endif

       do i = 1, Self%SpaceBra%GetNumChannels()
          if(present(bChan))then
             if(.not.self%mBlock( i, j )%IsType(bChan))cycle
          endif
          if(present(bIonLabel))then
             if( self%mBlock( i, j )%GetIonBraLabel() .isnt. bIonLabel )cycle
          endif
          bstat = 1
          if(TaskList(SESSES_LOAD )              ) call self%mblock(i,j)%Load ( bstat )
          if(TaskList(SESSES_BUILD).and.bstat/=0 ) call self%mblock(i,j)%Build( bstat, SESSES_build_switch )
          !.. Avoid computing the block below diagonal, for totally symmetric operators
          skip_LD = .false.
          if( totally_symmetric_operator )then
             if( i < j )then
                blk = self%mblock(i,j)
                call blk%Transpose()
                if(self%OpLabel(1:1)=="D") call blk%Scale(-1.d0)
                if(TaskList(SESSES_SAVE ).and.bstat==0 ) call blk%Save ( )
                if(TaskList(SESSES_FREE )              ) call blk%Free ( )
             endif
             skip_LD = i==j
          endif
          if(TaskList(SESSES_SAVE ).and.bstat==0 ) call self%mblock(i,j)%Save ( )
          if(TaskList(SESSES_FREE )              ) call self%mblock(i,j)%Free ( )
          if(skip_LD)exit
          !
       end do
    end do
  end subroutine ClassSESSESBlock_Driver

  !> Assemble in one matrix all the blocks in ClassSESSESBlock class.
  !! *** MUST ADD THE OPTION TO CUT OFF THE LAST BSPLINE, IF REQUESTED
  !! This is most easily done by altering the number of columns / rows in the object
  !! matrix, without any need of reallocating the object.
  !!
  subroutine ClassSESSESBlock_Assemble( Self, OpMat, bChan, kChan, bIonLabel, kIonLabel )
    !> Class to be assemble.
    class(ClassSESSESBlock)                   , intent(inout) :: Self
    type(ClassMatrix)                         , intent(out)   :: OpMat
    class(ClassCloseCouplingChannel), optional, intent(in)    :: bChan, kChan
    character(len=*)                , optional, intent(in)    :: bIonLabel,  kIonLabel
    !
    type(ClassMatrix) :: dmat
    integer :: i, j, i0, j0, ni, nj, nBra, nKet
    !
    call self%GetSize( nBra, nKet, bChan, kChan, bIonLabel, kIonLabel )
    call OpMat%InitFull(nBra,nKet)
    !
    i0 = 0
    do i = 1, self%SpaceBra%GetNumChannels()

       if(present(bChan))then
          if(.not.self%mBlock( i, 1 )%IsType(bChan))cycle
       endif
       if(present(bIonLabel))then
          if( self%mBlock( i, 1 )%GetIonBraLabel() .isnt. bIonLabel )cycle
       endif

       j0 = 0 
       do j = 1, self%SpaceBra%GetNumChannels()

          if(present(kChan))then
             if(.not.self%mBlock( i, j )%IsType(bChan,kChan))cycle
          endif
          if(present(kIonLabel))then
             if( self%mBlock( i, j )%GetIonKetLabel() .isnt. kIonLabel )cycle
          endif

          ni = self%mBlock( i, j )%GetNrows()
          nj = self%mBlock( i, j )%GetNColumns()

          call self%mBlock( i, j )%FetchMatrix( dMat )
          call OpMat%Set( dmat, i0+1, j0+1 )

          j0 = j0 + nj
       enddo

       i0 = i0 + ni
    enddo

  end subroutine ClassSESSESBlock_Assemble


  subroutine ClassSESSESBlock_GetSize( Self, bSize, kSize, bChan, kChan, bIonLabel, kIonLabel )
    !> Class to be assemble.
    class(ClassSESSESBlock)                   , intent(inout) :: Self
    integer                                   , intent(out)   :: bSize, kSize
    class(ClassCloseCouplingChannel), optional, intent(in)    :: bChan, kChan
    character(len=*)                , optional, intent(in)    :: bIonLabel, kIonLabel
    integer :: i, j
    bSize = 0
    do i = 1, self%SpaceBra%GetNumChannels()
       if(present(bChan))then
          if(.not.self%mBlock( i, 1 )%IsType(bChan=bChan))cycle
       endif
       if(present(bIonLabel))then
          if( self%mBlock( i, 1 )%GetIonBraLabel() .isnt. bIonLabel )cycle
       endif
       bSize=bSize+self%mBlock( i, 1 )%GetNRows()
    enddo
    kSize = 0
    do j = 1, self%SpaceKet%GetNumChannels()
       if(present(kChan))then
          if(.not.self%mBlock( 1, j )%IsType(kChan=kChan))cycle
       endif
       if(present(kIonLabel))then
          if( self%mBlock( 1, j )%GetIonKetLabel() .isnt. kIonLabel )cycle
       endif
       kSize=kSize+self%mBlock( 1, j )%GetNColumns()
    enddo
  end subroutine ClassSESSESBlock_GetSize


  !> Initializes the ClassSESSESBlock class.
  subroutine ClassSESSESBlock_Save( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    logical :: TaskList(N_SESSES_ID)
    TaskList=.FALSE.
    TaskList(SESSES_SAVE)=.TRUE.
    call self%driver(TaskList)
  end subroutine ClassSESSESBlock_Save


  !> Initializes the ClassSESSESBlock class.
  subroutine ClassSESSESBlock_Load( Self )
    class(ClassSESSESBlock), intent(inout) :: Self
    logical :: TaskList(N_SESSES_ID)
    TaskList=.FALSE.
    TaskList(SESSES_LOAD)=.TRUE.
    call self%driver(TaskList)
  end subroutine ClassSESSESBlock_Load


  logical function ValidSymmetries( OpLabel, BraSymSpace, KetSymSpace ) result(ValidSym)
    character(len=*)             , intent(in) :: OpLabel
    class(ClassSymESpace), target, intent(in) :: BraSymSpace
    class(ClassSymESpace), target, intent(in) :: KetSymSpace    
    type(ClassIrrep), pointer :: BraIrrep, KetIrrep, OpIrrep
    type(ClassGroup), pointer :: Group
    !
    ValidSym = .false.
    !
    BraIrrep => BraSymSpace%GetIrrep()
    KetIrrep => KetSymSpace%GetIrrep()
    Group    => BraIrrep%GetGroup()
    OpIrrep  =  GetOpIrrep( Group, OpLabel )
    !
    if ( ValidIrreps(BraIrrep,OpIrrep,KetIrrep) ) ValidSym = .true.
    !
  end function ValidSymmetries

end module ModuleSymESpace
