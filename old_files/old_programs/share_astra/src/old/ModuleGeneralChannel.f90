module ModuleGeneralChannel
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleParentIons
  use ModuleGroups

  implicit none
  private

  type, public :: ClassGeneralChannel
     type(ClassParentIon), pointer :: ParentIon
     integer                       :: TotalMultiplicity
     type(ClassIrrep),     pointer :: TotalIrrep
     type(ClassIrrep),     pointer :: OrbitalIrrep
     character(len=:), allocatable :: StorageDir
   contains
     procedure, public :: init0              => ClassGeneralChannelInit0
     procedure, public :: show               => ClassGeneralChannelShow
     procedure, public :: free               => ClassGeneralChannelfree
     procedure, public :: GetSize            => ClassGeneralChannelGetSize
     procedure, public :: GetStorageDir      => ClassGeneralChannelGetStorageDir
     procedure, public :: GetLabel           => ClassGeneralChannelGetLabel
!!$     procedure, public :: GetLabelExt        => ClassGeneralChannelGetLabelExt
     procedure, public :: GetPILabel         => ClassGeneralChannelGetPILabel
     procedure, public :: GetPI              => ClassGeneralChannelGetPI
     procedure, public :: GetPICharge        => ClassGeneralChannelGetPICharge
     procedure, public :: GetPIMultiplicity  => ClassGeneralChannelGetPIMultiplicity
     procedure, public :: GetPIEnergy        => ClassGeneralChannelGetPIEnergy
     procedure, public :: GetTotIrrepName    => ClassGeneralChannelGetTotIrrepName
     procedure, public :: GetTotIrrep        => ClassGeneralChannelGetTotIrrep
     procedure, public :: GetOrbIrrep        => ClassGeneralChannelGetOrbIrrep
     procedure, public :: GetTotMultiplicity => ClassGeneralChannelGetTotMultiplicity
     procedure, public :: GetSymLabel        => ClassGeneralChannelGetSymLabel
  end type ClassGeneralChannel

contains

  subroutine ClassGeneralChannelInit0(self,&
       ParentIon        , &
       TotalMultiplicity, &
       TotalIrrep       , &
       StorageDir       )
    class(ClassGeneralChannel)   , intent(inout) :: self
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
  end subroutine ClassGeneralChannelInit0

  subroutine ClassGeneralChannelFree( self )
    class(ClassGeneralChannel), intent(inout) :: self
  end subroutine ClassGeneralChannelFree

  subroutine ClassGeneralChannelShow( Self, unit )
    class(ClassGeneralChannel), intent(in) :: Self
    integer, optional         , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    call self%ParentIon%show()
    write(outunit,"(a,x,i0)")"2S+1=",self%TotalMultiplicity
    call self%TotalIrrep%Show()
    call self%OrbitalIrrep%Show()
    write(outunit,"(a)") "Storage = "//self%StorageDir
  end subroutine ClassGeneralChannelShow

  function ClassGeneralChannelGetStorageDir( Self )result(Dir)
    class(ClassGeneralChannel), intent(in) :: Self
    character(len=:), allocatable :: Dir
    if ( .not.allocated(Self%StorageDir) ) call Assert( &
         'The storage directory is not allocated, impossible to fetch.' )
    allocate( Dir, source = Self%StorageDir )
  end function ClassGeneralChannelGetStorageDir

  integer function ClassGeneralChannelGetSize( Self )result(isize)
    class(ClassGeneralChannel), intent(in) :: Self
    iSize=0
  end function ClassGeneralChannelGetSize

  function ClassGeneralChannelGetLabel( self ) result( Label )
    class(ClassGeneralChannel), intent(in) :: SELF
    character(len=:), allocatable :: Label
    allocate( Label, source = self%GetPILabel()//"x"//self%OrbitalIrrep%GetName() )
  end function ClassGeneralChannelGetLabel

!!$  function ClassGeneralChannelGetLabelExt( self ) result( Label )
!!$    class(ClassGeneralChannel), intent(in) :: SELF
!!$    character(len=:), allocatable :: Label
!!$    allocate( Label, source = "GEN" )
!!$  end function ClassGeneralChannelGetLabelExt
!!$
  function ClassGeneralChannelGetPILabel( SELF ) result( Label )
    class(ClassGeneralChannel), intent(in) :: SELF
    character(len=:), allocatable :: Label
    allocate( Label, source = SELF%Parention%GetLabel() )
  end function ClassGeneralChannelGetPILabel

  function ClassGeneralChannelGetPI( SELF ) result( PI )
    class(ClassGeneralChannel), target, intent(in) :: SELF
    type(ClassparentIon), pointer :: PI
    allocate( PI, source = SELF%Parention )
  end function ClassGeneralChannelGetPI

  integer function ClassGeneralChannelGetPICharge( SELF ) result( PICharge )
    class(ClassGeneralChannel), intent(in) :: SELF
    PICharge = SELF%ParentIon%GetCharge()
  end function ClassGeneralChannelGetPICharge

  integer function ClassGeneralChannelGetPIMultiplicity( SELF ) result( PIMult )
    class(ClassGeneralChannel), intent(in) :: SELF
    PIMult = SELF%ParentIon%GetMultiplicity()
  end function ClassGeneralChannelGetPIMultiplicity

  real(kind(1d0)) function ClassGeneralChannelGetPIEnergy( SELF ) result( PIE )
    class(ClassGeneralChannel), intent(inout) :: SELF
    PIE = SELF%ParentIon%GetEnergy()
  end function ClassGeneralChannelGetPIEnergy

  function ClassGeneralChannelGetTotIrrepName( Self ) result(IrrepName)
    class(ClassGeneralChannel), intent(in) :: Self
    character(len=:), allocatable :: IrrepName
    allocate( IrrepName, source = Self%TotalIrrep%GetName() )
  end function ClassGeneralChannelGetTotIrrepName

  function ClassGeneralChannelGetTotIrrep( Self ) result(Irrep)
    class(ClassGeneralChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self%TotalIrrep
  end function ClassGeneralChannelGetTotIrrep

  function ClassGeneralChannelGetOrbIrrep( Self ) result(Irrep)
    class(ClassGeneralChannel), target, intent(in) :: Self
    type(ClassIrrep), pointer :: Irrep
    Irrep => Self%OrbitalIrrep
  end function ClassGeneralChannelGetOrbIrrep

  function ClassGeneralChannelGetTotMultiplicity( Self ) result(Mult)
    class(ClassGeneralChannel), intent(in) :: Self
    integer :: Mult
    Mult = Self%TotalMultiplicity
  end function ClassGeneralChannelGetTotMultiplicity

  function ClassGeneralChannelGetSymLabel( Self ) result(Label)
    class(ClassGeneralChannel), intent(in) :: Self
    character(len=:), allocatable :: Label
    allocate( Label, source = AlphabeticNumber(Self%GetTotMultiplicity())//Self%GetTotIrrepName() )
  end function ClassGeneralChannelGetSymLabel
  
end module ModuleGeneralChannel
