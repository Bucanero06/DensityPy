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
module ModuleCloseCouplingChannels

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleString
  use ModuleErrorHandling
  use ModuleGroups
  use ModuleXlm
  use ModuleParentIons

  implicit none
  private

  character(len=*), public, parameter :: aiMOc_ID = "aiM"
  character(len=*), public, parameter :: viMOc_ID = "viM"
  character(len=*), public, parameter :: hiGOc_ID = "hiG"
  character(len=*), public, parameter :: veSOc_ID = "beS"

  character(len=*), private, parameter :: BSPLINE_ID_PREFIX= " " !"be"


  type, public :: ClassCloseCouplingChannel
     type(ClassParentIon), pointer :: ParentIon
     integer                       :: TotalMultiplicity
     type(ClassIrrep),     pointer :: TotalIrrep
     type(ClassIrrep),     pointer :: OrbitalIrrep
     character(len=:), allocatable :: StorageDir
   contains
     procedure, public :: init0              => ClassCloseCouplingChannelInit0
     procedure, public :: show               => ClassCloseCouplingChannelShow
     procedure, public :: free               => ClassCloseCouplingChannelfree
     procedure, public :: GetSize            => ClassCloseCouplingChannelGetSize
     procedure, public :: GetStorageDir      => ClassCloseCouplingChannelGetStorageDir
     procedure, public :: GetLabel           => ClassCloseCouplingChannelGetLabel
     procedure, public :: GetPILabel         => ClassCloseCouplingChannelGetPILabel
     procedure, public :: GetPI              => ClassCloseCouplingChannelGetPI
     procedure, public :: GetPICharge        => ClassCloseCouplingChannelGetPICharge
     procedure, public :: GetPIMultiplicity  => ClassCloseCouplingChannelGetPIMultiplicity
     procedure, public :: GetPIEnergy        => ClassCloseCouplingChannelGetPIEnergy
     procedure, public :: GetTotIrrepName    => ClassCloseCouplingChannelGetTotIrrepName
     procedure, public :: GetTotIrrep        => ClassCloseCouplingChannelGetTotIrrep
     procedure, public :: GetOrbIrrep        => ClassCloseCouplingChannelGetOrbIrrep
     procedure, public :: GetIonIrrep        => ClassCloseCouplingChannelGetIonIrrep
     procedure, public :: GetTotMultiplicity => ClassCloseCouplingChannelGetTotMultiplicity
     procedure, public :: GetSymLabel        => ClassCloseCouplingChannelGetSymLabel
  end type ClassCloseCouplingChannel


  type, public, extends(ClassCloseCouplingChannel) :: ClassActiveInternalChannel
     private
     integer :: norb
   contains
     procedure, public :: init        => ClassActiveInternalChannel_Init
     procedure, public :: show        => ClassActiveInternalChannel_Show
     procedure, public :: load        => ClassActiveInternalChannel_Load
     procedure, public :: save        => ClassActiveInternalChannel_save
     procedure, public :: free        => ClassActiveInternalChannel_Free
     procedure, public :: GetLabelExt => ClassActiveInternalChannel_GetLabelExt
  end type ClassActiveInternalChannel

  type, public, extends( ClassCloseCouplingChannel ) :: ClassVirtualInternalChannel
     private
     character(len=:), allocatable :: Label ! m -> virtual MOs, h -> hybrid orbitals
     integer   :: Lmax
     !*** We must provide a pointer to the expansion of the short-range orbitals
     !    in terms of the B-splines at the end. It is not clear whether we need
     !    to store here the maximum L in the expansion.
   contains
     procedure, public :: init        => ClassVirtualInternalChannelInit
     procedure, public :: show        => ClassVirtualInternalChannelShow
     procedure, public :: free        => ClassVirtualInternalChannelFree
     procedure, public :: GetLabelExt => ClassVirtualInternalChannelGetLabelExt
     procedure, public :: isMolecular => ClassVirtualInternalChannelIsMolecular
     procedure, public :: GetLmax     => ClassVirtualInternalChannelGetLmax

     final              :: ClassVirtualInternalChannelFinal
  end type ClassVirtualInternalChannel

  type, public, extends( ClassCloseCouplingChannel ) :: ClassVirtualExternalChannel
     private
     type(ClassXlm)                :: Xlm
     character(len=:), allocatable :: IDPrefix
   contains
     procedure, public  :: init          => ClassVirtualExternalChannelInit
     procedure, public  :: show          => ClassVirtualExternalChannelShow
     procedure, public  :: Save          => ClassVirtualExternalChannelSave
     procedure, public  :: Load          => ClassVirtualExternalChannelLoad
     procedure, public  :: free          => ClassVirtualExternalChannelFree
     procedure, public  :: GetL          => ClassVirtualExternalChannelGetL
     procedure, public  :: GetM          => ClassVirtualExternalChannelGetM
     procedure, public  :: GetLabelExt   => ClassVirtualExternalChannelGetLabelExt
     generic,   public  :: assignment(=) => ClassVirtualExternalChannelCopyVEC
     procedure, private ::                  ClassVirtualExternalChannelCopyVEC
  end type ClassVirtualExternalChannel

  
  public :: Active_Internal_Channel_Farm
  public :: Virtual_Internal_Channel_Farm
  public :: virtual_external_channel_farm

  
contains

  subroutine ClassCloseCouplingChannelInit0(self,&
       ParentIon        , &
       TotalMultiplicity, &
       TotalIrrep       , &
       StorageDir       )
    class(ClassCloseCouplingChannel)   , intent(inout) :: self
    class(ClassParentIon), target, intent(in)    :: ParentIon
    integer                      , intent(in)    :: TotalMultiplicity
    class(ClassIrrep),     target, intent(in)    :: TotalIrrep
    character(len=*)             , intent(in)    :: StorageDir
    self%ParentIon          =>  ParentIon
    self%TotalMultiplicity  =   TotalMultiplicity
    self%TotalIrrep         =>  TotalIrrep
    self%OrbitalIrrep       =>  TotalIrrep * ParentIon%GetIrrep()
    if(allocated(self%StorageDir))deallocate(self%StorageDir)
    allocate(self%StorageDir,source=StorageDir)
  end subroutine ClassCloseCouplingChannelInit0

  subroutine ClassCloseCouplingChannelFree( self )
    class(ClassCloseCouplingChannel), intent(inout) :: self
  end subroutine ClassCloseCouplingChannelFree

  subroutine ClassCloseCouplingChannelShow( Self, unit )
    class(ClassCloseCouplingChannel), intent(in) :: Self
    integer, optional         , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    call self%ParentIon%show()
    write(outunit,"(a,x,i0)")"2S+1=",self%TotalMultiplicity
    call self%TotalIrrep%Show()
    call self%OrbitalIrrep%Show()
    write(outunit,"(a)") "Storage = "//self%StorageDir
  end subroutine ClassCloseCouplingChannelShow

  function ClassCloseCouplingChannelGetStorageDir( Self )result(Dir)
    class(ClassCloseCouplingChannel), intent(in) :: Self
    character(len=:), allocatable :: Dir
    if ( .not.allocated(Self%StorageDir) ) call Assert( &
         'The storage directory is not allocated, impossible to fetch.' )
    allocate( Dir, source = Self%StorageDir )
  end function ClassCloseCouplingChannelGetStorageDir

  integer function ClassCloseCouplingChannelGetSize( Self )result(isize)
    class(ClassCloseCouplingChannel), intent(in) :: Self
    iSize=0
  end function ClassCloseCouplingChannelGetSize

  function ClassCloseCouplingChannelGetLabel( self ) result( Label )
    class(ClassCloseCouplingChannel), intent(in) :: SELF
    character(len=:), allocatable :: Label
    allocate( Label, source = self%GetPILabel()//"x"//self%OrbitalIrrep%GetName() )
  end function ClassCloseCouplingChannelGetLabel

  function ClassCloseCouplingChannelGetPILabel( SELF ) result( Label )
    class(ClassCloseCouplingChannel), intent(in) :: SELF
    character(len=:), allocatable :: Label
    allocate( Label, source = SELF%Parention%GetLabel() )
  end function ClassCloseCouplingChannelGetPILabel

  function ClassCloseCouplingChannelGetPI( SELF ) result( PI )
    class(ClassCloseCouplingChannel), target, intent(in) :: SELF
    type(ClassparentIon), pointer :: PI
    allocate( PI, source = SELF%Parention )
  end function ClassCloseCouplingChannelGetPI

  integer function ClassCloseCouplingChannelGetPICharge( SELF ) result( PICharge )
    class(ClassCloseCouplingChannel), intent(in) :: SELF
    PICharge = SELF%ParentIon%GetCharge()
  end function ClassCloseCouplingChannelGetPICharge

  integer function ClassCloseCouplingChannelGetPIMultiplicity( SELF ) result( PIMult )
    class(ClassCloseCouplingChannel), intent(in) :: SELF
    PIMult = SELF%ParentIon%GetMultiplicity()
  end function ClassCloseCouplingChannelGetPIMultiplicity

  real(kind(1d0)) function ClassCloseCouplingChannelGetPIEnergy( SELF ) result( PIE )
    class(ClassCloseCouplingChannel), intent(inout) :: SELF
    PIE = SELF%ParentIon%GetEnergy()
  end function ClassCloseCouplingChannelGetPIEnergy

  function ClassCloseCouplingChannelGetTotIrrepName( Self ) result(IrrepName)
    class(ClassCloseCouplingChannel), intent(in) :: Self
    character(len=:), allocatable :: IrrepName
    allocate( IrrepName, source = Self%TotalIrrep%GetName() )
  end function ClassCloseCouplingChannelGetTotIrrepName

  function ClassCloseCouplingChannelGetTotIrrep( Self ) result(Irrep)
    class(ClassCloseCouplingChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self%TotalIrrep
  end function ClassCloseCouplingChannelGetTotIrrep

  function ClassCloseCouplingChannelGetOrbIrrep( Self ) result(Irrep)
    class(ClassCloseCouplingChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self%OrbitalIrrep
  end function ClassCloseCouplingChannelGetOrbIrrep

  function ClassCloseCouplingChannelGetIonIrrep( Self ) result(Irrep)
    class(ClassCloseCouplingChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self%ParentIon%GetIrrep()
  end function ClassCloseCouplingChannelGetIonIrrep

  function ClassCloseCouplingChannelGetTotMultiplicity( Self ) result(Mult)
    class(ClassCloseCouplingChannel), intent(in) :: Self
    integer :: Mult
    Mult = Self%TotalMultiplicity
  end function ClassCloseCouplingChannelGetTotMultiplicity

  function ClassCloseCouplingChannelGetSymLabel( Self ) result(Label)
    class(ClassCloseCouplingChannel), intent(in) :: Self
    character(len=:), allocatable :: Label
    allocate( Label, source = AlphabeticNumber(Self%GetTotMultiplicity())//Self%GetTotIrrepName() )
  end function ClassCloseCouplingChannelGetSymLabel

  !------------------------------------------------
  ! AIC
  function Active_Internal_Channel_Farm(&
       ParentIon        , &
       TotalMultiplicity, &
       TotalIrrep       , &
       StorageDir       ) &
       result ( aic )
    class(ClassParentIon), target  , intent(in) :: ParentIon
    integer                        , intent(in) :: TotalMultiplicity
    class(ClassIrrep),     target  , intent(in) :: TotalIrrep
    character(len=*)               , intent(in) :: StorageDir
    type(ClassActiveInternalChannel)    , pointer    :: aic
    allocate(aic)
    call aic%init(          &
         ParentIon        , &
         TotalMultiplicity, &
         TotalIrrep       , &
         StorageDir       )
  end function Active_Internal_Channel_Farm

  subroutine ClassActiveInternalChannel_Init(self,&
       ParentIon        , &
       TotalMultiplicity, &
       TotalIrrep       , &
       StorageDir       )
    class(ClassActiveInternalChannel), intent(inout) :: self
    class(ClassParentIon), target  , intent(in) :: ParentIon
    integer                        , intent(in) :: TotalMultiplicity
    class(ClassIrrep),     target  , intent(in) :: TotalIrrep
    character(len=*)               , intent(in) :: StorageDir
    call self%Init0(ParentIon,TotalMultiplicity,TotalIrrep,StorageDir)
    self%norb=0 !*** < must pass these somehow
  end subroutine ClassActiveInternalChannel_Init

  subroutine ClassActiveInternalChannel_Show(self,unit)
    class(ClassActiveInternalChannel), intent(in) :: self
    integer, optional             , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    call self%ClassCloseCouplingChannel%show(unit)
    if(present(unit)) outunit = unit
    write(outunit,"(i0)") self%norb
  end subroutine ClassActiveInternalChannel_Show

  subroutine ClassActiveInternalChannel_Load(self)
    class(ClassActiveInternalChannel), intent(inout) :: self
  end subroutine ClassActiveInternalChannel_Load

  subroutine ClassActiveInternalChannel_Save(self)
    class(ClassActiveInternalChannel), intent(in) :: self
  end subroutine ClassActiveInternalChannel_Save

  subroutine ClassActiveInternalChannel_Free(self)
    class(ClassActiveInternalChannel), intent(inout) :: self
  end subroutine ClassActiveInternalChannel_Free

  function ClassActiveInternalChannel_GetLabelExt(self) result( label )
    class(ClassActiveInternalChannel), intent(in) :: self
    character(len=:), allocatable :: label
    allocate(label,source=aiMOc_ID)
  end function ClassActiveInternalChannel_GetLabelExt

  !------------------------------------------------
  ! VIC
  function Virtual_Internal_Channel_Farm( &
       ParentIon, Label, Lmax, TotalMultiplicity, TotalIrrep, StorageDir ) &
       result ( Chan )
    class(ClassParentIon), target  , intent(in) :: ParentIon
    character(len=*)               , intent(in) :: Label
    integer                        , intent(in) :: Lmax
    integer                        , intent(in) :: TotalMultiplicity
    class(ClassIrrep),     target  , intent(in) :: TotalIrrep
    character(len=*)               , intent(in) :: StorageDir
    type(ClassVirtualInternalChannel)   , pointer    :: Chan
    allocate(Chan)
    call Chan%init( ParentIon, Label, Lmax, TotalMultiplicity, TotalIrrep, StorageDir )
  end function Virtual_Internal_Channel_Farm

  subroutine ClassVirtualInternalChannelInit( self, ParentIon, Label, Lmax, TotalMultiplicity, &
       TotalIrrep, StorageDir )
    class(ClassVirtualInternalChannel) , intent(inout) :: self
    class(ClassParentIon), target , intent(in)    :: ParentIon
    character(len=*)              , intent(in)    :: Label
    integer                       , intent(in)    :: Lmax
    integer                       , intent(in)    :: TotalMultiplicity
    class(ClassIrrep),     target , intent(in)    :: TotalIrrep
    character(len=*)              , intent(in)    :: StorageDir
    call self%Init0( ParentIon, TotalMultiplicity, TotalIrrep, StorageDir )
    self%Lmax  = Lmax
    self%Label = Label
  end subroutine ClassVirtualInternalChannelInit

  subroutine ClassVirtualInternalChannelShow( self, unit )
    class(ClassVirtualInternalChannel) , intent(in) :: self
    integer, optional             , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    call self%ClassCloseCouplingChannel%show(unit)
    if(present(unit)) outunit = unit
    write(outunit,"(a)") "Short Range Channel Info : "
    call self%ParentIon%Show(unit)
    write(outunit,"(a,i4)"          ) "  Maximum Angular Momentum :",self%Lmax
    write(outunit,"(a,i4)"          ) "  Channel Multiplicity :",self%TotalMultiplicity
    write(outunit,"(a)",advance="no") "  Channel Symmetry     :"
    call self%TotalIrrep%show(unit)
  end subroutine ClassVirtualInternalChannelShow

  logical function ClassVirtualInternalChannelIsMolecular( self ) result( isMolecular )
    class(ClassVirtualInternalChannel), intent(in) :: self
    isMolecular = self%Label .is. viMOc_ID
  end function ClassVirtualInternalChannelIsMolecular

  function ClassVirtualInternalChannelGetLabelExt( self ) result( label )
    class(ClassVirtualInternalChannel) , intent(in) :: self
    character(len=:), allocatable :: label
    allocate(label,source=self%label)
  end function ClassVirtualInternalChannelGetLabelExt
     
  integer function ClassVirtualInternalChannelGetLmax( SRC ) result( LMax )
    class(ClassVirtualInternalChannel), intent(in) :: SRC
    LMax = SRC%Lmax
  end function ClassVirtualInternalChannelGetLmax

  subroutine ClassVirtualInternalChannelFree( self )
    class(ClassVirtualInternalChannel), intent(inout) :: self
    self%ParentIon  => NULL()
    self%TotalIrrep => NULL()
    self%Lmax = -1
    self%TotalMultiplicity = 0
  end subroutine ClassVirtualInternalChannelFree
  
  subroutine ClassVirtualInternalChannelFinal( Channel )
    type(ClassVirtualInternalChannel) :: Channel
    call Channel%free()
  end subroutine ClassVirtualInternalChannelFinal

  !------------------------------------------------
  ! VEC
  function virtual_external_channel_farm( &
       ParentIon        , &
       TotalMultiplicity, &
       TotalIrrep       , &
       StorageDir       , &
       Xlm              ) &
       result ( vec )
    class(ClassParentIon), target  , intent(in) :: ParentIon
    integer                        , intent(in) :: TotalMultiplicity
    class(ClassIrrep),     target  , intent(in) :: TotalIrrep
    character(len=*)               , intent(in) :: StorageDir
    class(ClassXlm)                , intent(in) :: Xlm
    type(ClassVirtualExternalChannel), pointer  :: vec
    allocate(vec)
    call vec%Init(&
         ParentIon        , &
         TotalMultiplicity, &
         TotalIrrep       , &
         StorageDir       , &
         Xlm              )
  end function Virtual_External_Channel_Farm
  
  subroutine ClassVirtualExternalChannelInit( &
       self          , &
       ParentIon        , &
       TotalMultiplicity, &
       TotalIrrep       , &
       StorageDir       , &
       Xlm              )
    class(ClassVirtualExternalChannel) , intent(inout) :: self
    class(ClassParentIon), target  , intent(in)    :: ParentIon
    integer                        , intent(in)    :: TotalMultiplicity
    class(ClassIrrep),     target  , intent(in)    :: TotalIrrep
    character(len=*)               , intent(in)    :: StorageDir
    class(ClassXlm)                , intent(in)    :: Xlm
    self%IDPrefix = trim(adjustl(BSPLINE_ID_PREFIX))
    call self%Init0(ParentIon,TotalMultiplicity, &
       TotalIrrep,StorageDir)
    self%Xlm = Xlm
  end subroutine ClassVirtualExternalChannelInit

  subroutine ClassVirtualExternalChannelShow( self, unit )
    class(ClassVirtualExternalChannel), intent(in) :: self
    integer, optional             , intent(in) :: unit
    integer :: outunit
    call self%ClassCloseCouplingChannel%show(unit)
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    write(outunit,"(a)") "Partial Wave Channels Info : "
    call self%ParentIon%Show(unit)
    write(outunit,"(a,i4)"          ) "  Total Multiplicity :",self%TotalMultiplicity
    write(outunit,"(a)",advance="no") "  Total Symmetry     :"
    call self%TotalIrrep%show(unit)
    write(outunit,"(a)",advance="no") "  Orbital Symmetry     :"
    call self%OrbitalIrrep%show(unit)
    call self%Xlm%show(unit)
  end subroutine ClassVirtualExternalChannelShow

  subroutine ClassVirtualExternalChannelSave( self, unit )
    class(ClassVirtualExternalChannel), intent(in) :: self
    integer                       , intent(in) :: unit
    call self%ParentIon%SaveData(unit)
    write(unit,*) self%TotalMultiplicity
    call self%Xlm%Save(unit)
  end subroutine ClassVirtualExternalChannelSave

  subroutine ClassVirtualExternalChannelLoad( self, unit )
    class(ClassVirtualExternalChannel), intent(inout) :: self
    integer                       , intent(in)    :: unit
    call self%ParentIon%LoadData(unit)
    read(unit,*) self%TotalMultiplicity
    call self%Xlm%Load(unit)
  end subroutine ClassVirtualExternalChannelLoad

  subroutine ClassVirtualExternalChannelFree( self )
    class(ClassVirtualExternalChannel), intent(inout) :: self
    self%ParentIon    => NULL()
    self%TotalIrrep   => NULL()
    self%OrbitalIrrep => NULL()
    self%TotalMultiplicity = 0
    call self%Xlm%init(0,0)
  end subroutine ClassVirtualExternalChannelFree

  integer function ClassVirtualExternalChannelGetL( self ) result( L )
    class(ClassVirtualExternalChannel), intent(in) :: self
    L = self%Xlm%GetL()
  end function ClassVirtualExternalChannelGetL

  integer function ClassVirtualExternalChannelGetM( self ) result( M )
    class(ClassVirtualExternalChannel), intent(in) :: self
    M = self%Xlm%GetM()
  end function ClassVirtualExternalChannelGetM

  function ClassVirtualExternalChannelGetLabelExt( self ) result( Label )
    class(ClassVirtualExternalChannel), intent(in) :: self
    character(len=:), allocatable :: Label
    allocate( Label, source = self%IDPrefix//self%Xlm%GetLabel() )
  end function ClassVirtualExternalChannelGetLabelExt

  !> Copies the information of VEC class to another.
  subroutine ClassVirtualExternalChannelCopyVEC( VECout, VECin )
    class(ClassVirtualExternalChannel)        , intent(inout) :: VECout
    class(ClassVirtualExternalChannel), target, intent(in)    :: VECin
    call VECout%Free()
    VECout%ParentIon         => VECin%ParentIon
    VECout%TotalMultiplicity =  VECin%TotalMultiplicity
    VECout%TotalIrrep        => VECin%TotalIrrep
    VECout%OrbitalIrrep      => VECin%OrbitalIrrep
    VECout%Xlm               =  VECin%Xlm
    allocate( VECout%StorageDir, source = VECin%StorageDir )
  end subroutine ClassVirtualExternalChannelCopyVEC

  subroutine ClassVirtualExternalChannelFinal( self )
    type(ClassVirtualExternalChannel) :: self
    call self%free()
  end subroutine ClassVirtualExternalChannelFinal
  
end module ModuleCloseCouplingChannels
