module ModuleSymmetricElectronicSpace

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleGroups
  use ModuleParentIons
  use ModuleSymmetricElectronicChannel
  use ModuleSymmetricLocalizedStates
  use ModuleMatrix
  use ModuleConstants
  use ModulePartialWaveChannel
  use ModuleShortRangeChannel
  use ModuleIO

  implicit none

  private

  type, public :: ClassSymmetricElectronicSpace

     private
     !> Maximum angular momentum the space can represent.
     integer                                            :: Lmax
     !> Total number of electrons.
     integer                                            :: Ne
     !> Spin multiplicity
     !! If spin is not conserved / well-defined, then Multiplicity = 0.
     integer                                            :: Multiplicity

     !> Point Group
     type(ClassGroup), pointer                          :: Group

     !> Spatial symmetry 
     type(ClassIrrep), pointer                          :: Irrep

     !> Close Coupling Sector of the space
     integer                                            :: NChannels

     type(ClassSymmetricElectronicChannel), allocatable :: Chanv(:)

     !> Symmetric Localized States (these are the bound states that are
     !! not associated to any specific ion).
     type(ClassSymmetricLocalizedStates)                :: LS

     !> Directory in which the results will be stored
     character(len=:), allocatable                      :: StorageDir

     !> Label of Parametric configuration (e.g., nuclear geometry)
     character(len=:), allocatable                      :: NuclearLabel

     logical                                            :: BoxStatesOnly = .FALSE.

   contains
     generic, public :: parseConfigFile     =>  ClassSymmetricElectronicSpaceParseConfigFile
     generic, public :: SetBoxOnly          => ClassSymmetricElectronicSpaceSetBoxOnly
     generic, public :: IsBoxOnly           => ClassSymmetricElectronicSpaceIsBoxOnly
     generic, public :: free                => ClassSymmetricElectronicSpaceFree
     generic, public :: show                => ClassSymmetricElectronicSpaceShow
     generic, public :: SetGroup            =>  ClassSymmetricElectronicSpaceSetGroup
     generic, public :: SetIrrep            =>  ClassSymmetricElectronicSpaceSetIrrep
     generic, public :: SetLSIrrep          =>  ClassSymmetricElectronicSpaceSetLSIrrep
     generic, public :: SetMultiplicity     =>  ClassSymmetricElectronicSpaceSetMultiplicity
     generic, public :: SetLmax             =>  ClassSymmetricElectronicSpaceSetLmax
     generic, public :: GetLmax             => ClassSymmetricElectronicSpaceGetLmax
     !> Gets the maximum L defined in the Close-Coupling anzat
     !! for the current ClassSymmetricElectronicSpace%
     generic, public :: GetMaxLinCC        => ClassSymmetricElectronicSpaceGetMaxLinCC
     generic, public :: GetTotNumElectrons  => ClassSymmetricElectronicSpaceGetTotNumElectrons
     generic, public :: GetNumQCOrbitals    => ClassSymmetricElectronicSpaceGetNumQCOrbitals
     generic, public :: GetNumChannels      => ClassSymmetricElectronicSpaceGetNumChannels
     generic, public :: GetStorageDir       => ClassSymmetricElectronicSpaceGetStorageDir
     generic, public :: GetSymStorageDir    => ClassSymmetricElectronicSpaceGetSymStorageDir
     generic, public :: GetNumPWC           => ClassSymmetricElectronicSpaceGetNumPWC
     generic, public :: PWCChannelIndex     => ClassSymmetricElectronicSpacePWCChannelIndex
     generic, public :: GetLabel            => ClassSymmetricElectronicSpaceGetLabel
     generic, public :: GetIrrepLabel       => ClassSymmetricElectronicSpaceGetIrrepLabel
     generic, public :: GetIrrep            => ClassSymmetricElectronicSpaceGetIrrep
     generic, public :: GetPILabel          => ClassSymmetricElectronicSpaceGetPILabel
     generic, public :: GetPICharge         => ClassSymmetricElectronicSpaceGetPICharge
     generic, public :: GetEnergyPI         => ClassSymmetricElectronicSpaceGetEnergyPI
     generic, public :: GetChannelPIEnergy  => ClassSymmetricElectronicSpaceGetChannelPIEnergy
     generic, public :: GetChannelL         => ClassSymmetricElectronicSpaceGetChannelL
     generic, public :: GetL                => ClassSymmetricElectronicSpaceGetL
     generic, public :: GetM                => ClassSymmetricElectronicSpaceGetM
     generic, public :: GetChannelSRC       => ClassSymmetricElectronicSpaceGetChannelSRC
     generic, public :: GetChannelPWC       => ClassSymmetricElectronicSpaceGetChannelPWC, &
                        ClassSymmetricElectronicSpaceGetChannelPWCFromAbsoluteIndex
     generic, public :: GetChannelLS        => ClassSymmetricElectronicSpaceGetChannelLS
     generic, public :: UnsetSpin           =>  ClassSymmetricElectronicSpaceUnsetSpin
     generic, public :: SetStorage          =>  ClassSymmetricElectronicSpaceSetStorage
     procedure, public :: SetNuclearLabel          => ClassSymmetricElectronicSpaceSetNuclearLabel

     procedure, private :: ClassSymmetricElectronicSpaceParseConfigFile
     procedure, private :: ClassSymmetricElectronicSpaceShow
     procedure, private :: ClassSymmetricElectronicSpaceFree
     procedure, private :: ClassSymmetricElectronicSpaceSetBoxOnly
     procedure, private :: ClassSymmetricElectronicSpaceIsBoxOnly
     procedure, private :: ClassSymmetricElectronicSpaceSetGroup
     procedure, private :: ClassSymmetricElectronicSpaceSetIrrep
     procedure, private :: ClassSymmetricElectronicSpaceSetLSIrrep
     procedure, private :: ClassSymmetricElectronicSpaceSetMultiplicity
     procedure, private :: ClassSymmetricElectronicSpaceSetLmax
     procedure, private :: ClassSymmetricElectronicSpaceGetTotNumElectrons
     procedure, private :: ClassSymmetricElectronicSpaceGetLmax
     procedure, private :: ClassSymmetricElectronicSpaceGetMaxLinCC
     procedure, private :: ClassSymmetricElectronicSpaceGetStorageDir
     procedure, private :: ClassSymmetricElectronicSpaceGetSymStorageDir
     procedure, private :: ClassSymmetricElectronicSpaceGetNumQCOrbitals
     procedure, private :: ClassSymmetricElectronicSpaceGetNumChannels
     procedure, private :: ClassSymmetricElectronicSpaceGetNumPWC
     procedure, private :: ClassSymmetricElectronicSpacePWCChannelIndex
     procedure, private :: ClassSymmetricElectronicSpaceGetLabel
     procedure, private :: ClassSymmetricElectronicSpaceGetIrrepLabel
     procedure, private :: ClassSymmetricElectronicSpaceGetIrrep
     procedure, private :: ClassSymmetricElectronicSpaceGetPILabel
     procedure, private :: ClassSymmetricElectronicSpaceGetPICharge
     procedure, private :: ClassSymmetricElectronicSpaceGetEnergyPI
     procedure, private :: ClassSymmetricElectronicSpaceGetChannelPIEnergy
     procedure, private :: ClassSymmetricElectronicSpaceGetChannelL
     procedure, private :: ClassSymmetricElectronicSpaceGetL
     procedure, private :: ClassSymmetricElectronicSpaceGetM
     procedure, private :: ClassSymmetricElectronicSpaceGetChannelSRC
     procedure, private :: ClassSymmetricElectronicSpaceGetChannelPWC
     procedure, private :: ClassSymmetricElectronicSpaceGetChannelPWCFromAbsoluteIndex
     procedure, private :: ClassSymmetricElectronicSpaceGetChannelLS
     procedure, private :: ClassSymmetricElectronicSpaceUnsetSpin
     procedure, private :: ClassSymmetricElectronicSpaceSetStorage
     final              :: ClassSymmetricElectronicSpaceFinal
  end type ClassSymmetricElectronicSpace

contains

  subroutine ClassSymmetricElectronicSpaceParseConfigFile( Space, &
       FileName, iflag, OnlyPrepareData )
    use moduleParameterList
    use moduleString
    use moduleIO
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    character(len=*)                    , intent(in)    :: FileName
    integer                             , intent(out)   :: iflag ! 0 => OK, n => some error
    logical, optional                   , intent(in)    :: OnlyPrepareData   
    !
    integer :: j, ichar, ichar2, dichar
    character(len=1024)      :: PWCStrn
    character(len=:), allocatable :: FullText
    character(len=:), allocatable :: GroupLabel, GeneralLabel, SymStrn
    integer :: NumChannels
    integer :: ParentIonCharge
    !
    !.. Parse the configuration file assuming the following format
    !
    !   Group = D2h
    !   [2Ag]{ 
    !          1Ag.1  ( 0  0, 2 0, 2 2 )
    !          3B2u.1 ( 1 -1, 3 -1, 3 -3 )
    !   }
    !
    !   [3B2u]{ 
    !          1Ag.1  ( ... )
    !          3B2u.1 ( ... )
    !          ...
    !   }
    !..


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

    call FetchGlobalVariable( FullText, "PARENT_ION_CHARGE"    , GeneralLabel )
    if(.not.allocated(GeneralLabel)) call Assert("Charge label missing in "//trim(FileName))
    read(GeneralLabel,*) ParentIonCharge

    call FetchGlobalVariable( FullText, "LMAX"      , GeneralLabel )
    if(.not.allocated(GeneralLabel)) call Assert("lmax label missing in "//trim(FileName))
    read(GeneralLabel,*) Space%Lmax

    call FetchGlobalVariable( FullText, "NELECTRONS", GeneralLabel )
    if(.not.allocated(GeneralLabel)) call Assert("Nelectrons label missing in "//trim(FileName))
    read(GeneralLabel,*) Space%Ne

    !.. 4. Isolate the irrep in the file 
    !..
    SymStrn = Space%GetLabel()

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
    pwcCycle : do 
       dichar = index(FullText(ichar+1:),")")
       if(dichar<1)exit pwcCycle
       ichar = ichar + dichar
       NumChannels = NumChannels+1
    enddo pwcCycle
    Space%NChannels = NumChannels
    allocate( Space%Chanv(NumChannels) )

    write(OUTPUT_UNIT,"(a)")
    write(OUTPUT_UNIT,"(a,i0,a)") " "//Space%GetLabel()//" ",Space%NChannels," Ionic Channels "
    
    call Space%LS%SetStorageDir( Space%StorageDir )
    call Space%LS%SetIrrep( Space%Irrep )
    call Space%LS%SetMultiplicity( Space%Multiplicity )

    !.. Parse the description of each PWC
    !..
    ichar = 0
    do j = 1, Space%NChannels

       ichar2 = ichar
       ichar = ichar+index(FullText(ichar+1:),")")
       PWCStrn=adjustl(FullText(ichar2+1:ichar))

       call Space%Chanv(j)%SetLmax( Space%Lmax )
       call Space%Chanv(j)%SetTotalMultiplicity( Space%Multiplicity )
       call Space%Chanv(j)%SetNelectrons( Space%Ne )
       call Space%Chanv(j)%SetTotalIrrep( Space%Irrep )
       call Space%Chanv(j)%SetParentIonCharge( ParentIonCharge )
       call Space%Chanv(j)%SetStorageDir( Space%StorageDir )
       write(OUTPUT_UNIT,"(4x,a)") trim(PWCStrn)
       if(present(OnlyPrepareData))then
          call Space%Chanv(j)%ParseConfLine( PWCStrn, OnlyPrepareData )
       else
          call Space%Chanv(j)%ParseConfLine( PWCStrn )
       endif

    enddo

  end subroutine ClassSymmetricElectronicSpaceParseConfigFile

  subroutine ClassSymmetricElectronicSpaceShow( Space, unit )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer, optional                   , intent(in) :: unit
    integer :: outunit
    !
    integer :: iChan
    !
    outunit = OUTPUT_UNIT
    !
    if(present(unit)) outunit = unit
    !
    write(outunit,"(a)") "Symmetric Electronic Space Info : "
    !
    write(outunit,"(a,i4)"          ) "  Maximum Angular Momentum :",Space%Lmax
    write(outunit,"(a,i4)"          ) "  Total Number of Electrons :",Space%Ne
    write(outunit,"(a,i4)"          ) "  Space Multiplicity :",Space%Multiplicity
    !
    write(outunit,"(a)",advance="no") "  Space Symmetry     :"
    !
    call Space%Irrep%show( outunit )
    !
    do iChan = 1, Space%NChannels
       call Space%Chanv(iChan)%show( outunit )
    enddo
    !
    call Space%LS%show( outunit )
    !
    write(outunit,"(a,a)"          ) "  Space Storage Dir :",Space%StorageDir
    write(outunit,"(a,a)"          ) "  Nuclear Geom Label :",Space%NuclearLabel
    !
  end subroutine ClassSymmetricElectronicSpaceShow

  subroutine ClassSymmetricElectronicSpaceFree( Space )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    Space%Lmax = -1
    Space%Ne = -1
    Space%Multiplicity = 0
    Space%Group => NULL()
    Space%Irrep => NULL()
!!$    Space%Basis => NULL()
    Space%NChannels = 0
    if ( allocated(Space%Chanv) ) deallocate( Space%Chanv )
    call Space%LS%Free()
    if ( allocated(Space%StorageDir) ) deallocate( Space%StorageDir )
    if ( allocated(Space%NuclearLabel) ) deallocate( Space%NuclearLabel )
  end subroutine ClassSymmetricElectronicSpaceFree

  subroutine ClassSymmetricElectronicSpaceSetBoxOnly( Space )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    Space%BoxStatesOnly = .TRUE.
  end subroutine ClassSymmetricElectronicSpaceSetBoxOnly

  subroutine ClassSymmetricElectronicSpaceUnSetBoxOnly( Space )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    Space%BoxStatesOnly = .FALSE.
  end subroutine ClassSymmetricElectronicSpaceUnSetBoxOnly

  logical function ClassSymmetricElectronicSpaceIsBoxOnly( Space ) result( IsBoxOnly )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    IsBoxOnly = Space%BoxStatesOnly
  end function ClassSymmetricElectronicSpaceIsBoxOnly

  subroutine ClassSymmetricElectronicSpaceFinal( Space )
    type(ClassSymmetricElectronicSpace) :: Space
    call Space%free()
  end subroutine ClassSymmetricElectronicSpaceFinal

  subroutine ClassSymmetricElectronicSpaceSetGroup( Space, Group )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    type(ClassGroup), target,             intent(in)    :: Group
    if ( associated(Space%Group) ) Space%Group => NULL()
    Space%Group => Group
  end subroutine ClassSymmetricElectronicSpaceSetGroup

  subroutine ClassSymmetricElectronicSpaceSetIrrep( Space, Irrep )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    type(ClassIrrep), target,             intent(in)    :: Irrep
    if ( associated(Space%Irrep) ) Space%Irrep => NULL()
    Space%Irrep => Irrep
  end subroutine ClassSymmetricElectronicSpaceSetIrrep

  subroutine ClassSymmetricElectronicSpaceSetLSIrrep( Space, Irrep )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    type(ClassIrrep), target,             intent(in)    :: Irrep
    call Space%LS%SetIrrep( Irrep )
  end subroutine ClassSymmetricElectronicSpaceSetLSIrrep

  subroutine ClassSymmetricElectronicSpaceSetMultiplicity( Space, MUltiplicity )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    integer,                              intent(in)    :: Multiplicity
    Space%Multiplicity = Multiplicity
  end subroutine ClassSymmetricElectronicSpaceSetMultiplicity

  subroutine ClassSymmetricElectronicSpaceSetLmax( Space, Lmax )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    integer,                              intent(in)    :: Lmax
    Space%Lmax = Lmax
  end subroutine ClassSymmetricElectronicSpaceSetLmax

  integer function ClassSymmetricElectronicSpaceGetLmax( Space ) result( Lmax )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    Lmax = Space%Lmax 
  end function ClassSymmetricElectronicSpaceGetLmax

  !> Gets the maximum L defined in the Close-Coupling anzat
  !! for the current ClassSymmetricElectronicSpace%
  integer function ClassSymmetricElectronicSpaceGetMaxLinCC( Space ) result( L )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer :: i
    if ( .not.allocated(Space%Chanv) ) call Assert( &
         'ClassSymmetricElectronicSpace has not been '//&
         'initialized, impossible to get the maximum L in '//&
         'the Close-Coupling expansion, aborting ...' )
    L = 0
    do i = 1, Space%NChannels
       L = max(L,Space%Chanv(i)%GetMaxLinCC())
    end do
  end function ClassSymmetricElectronicSpaceGetMaxLinCC

  integer function ClassSymmetricElectronicSpaceGetTotNumElectrons( Space ) result( Ne )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    Ne = Space%Ne
  end function ClassSymmetricElectronicSpaceGetTotNumElectrons

  function ClassSymmetricElectronicSpaceGetStorageDir( Space ) result( Dir )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    character(len=:), allocatable :: Dir
    if ( .not.allocated(Space%StorageDir) ) call Assert( &
         'Impossible to get the storage dir from the symmetric '//&
         'electronic space because is not allocated.' )
    allocate( Dir, source = Space%StorageDir  )
  end function ClassSymmetricElectronicSpaceGetStorageDir

  function ClassSymmetricElectronicSpaceGetSymStorageDir( Space ) result( Dir )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    character(len=:), allocatable :: Dir
    allocate( Dir, source = Space%LS%GetStorageDir()  )
  end function ClassSymmetricElectronicSpaceGetSymStorageDir

  integer function ClassSymmetricElectronicSpaceGetNumQCOrbitals( Space, Channel ) result( NQCO )
    !
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer, optional,                    intent(in) :: Channel
    !
    integer :: i
    !
    if ( present(Channel) ) then
       !
       if ( (Channel<1) .or. (Channel>Space%NChannels) ) then
          call Assert( "The index of the channel is invalid, impossible to get the number of QC orbitals." )
       end if
       !
       NQCO = Space%Chanv(Channel)%GetNumQCOrbitals()
       !
    else
       !
       NQCO = 0
       do i = 1, Space%NChannels
          NQCO = NQCO + Space%Chanv(i)%GetNumQCOrbitals()
       end do
       !
    end if
    !
  end function ClassSymmetricElectronicSpaceGetNumQCOrbitals

  integer function ClassSymmetricElectronicSpaceGetNumChannels( Space ) result(NC)
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    NC = Space%NChannels
  end function ClassSymmetricElectronicSpaceGetNumChannels

  function ClassSymmetricElectronicSpaceGetLabel( Space ) result( Label )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    character(len=:), allocatable :: Label
    allocate( Label, source = AlphabeticNumber(Space%Multiplicity)//Space%Irrep%GetName() )
  end function ClassSymmetricElectronicSpaceGetLabel

  function ClassSymmetricElectronicSpaceGetIrrepLabel( Space ) result( Label )
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    character(len=:), allocatable :: Label
    allocate( Label, source = Space%Irrep%GetName() )
  end function ClassSymmetricElectronicSpaceGetIrrepLabel

  function ClassSymmetricElectronicSpaceGetIrrep( Space ) result( SymIrrep )
    class(ClassSymmetricElectronicSpace), target, intent(in) :: Space
    type(ClassIrrep), pointer :: SymIrrep
    allocate( SymIrrep, source = Space%Irrep )
  end function ClassSymmetricElectronicSpaceGetIrrep

  character(len=128) function ClassSymmetricElectronicSpaceGetPILabel( Space, Channel ) result( Label )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    integer,                              intent(in)    :: Channel
    Label = Space%Chanv(Channel)%GetPILabelFun( ) 
  end function ClassSymmetricElectronicSpaceGetPILabel

  real(kind(1d0)) function ClassSymmetricElectronicSpaceGetPICharge( Space, Channel ) result(Charge)
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer,                              intent(in) :: Channel
    Charge = Space%Chanv(Channel)%GetParentIonCharge()
  end function ClassSymmetricElectronicSpaceGetPICharge

  integer function ClassSymmetricElectronicSpaceGetL( Space, Channel, PWC ) result(L)
    class(ClassSymmetricElectronicSpace), intent(in)    :: Space
    integer,                              intent(in)    :: Channel
    integer,                              intent(in)    :: PWC
    L = Space%Chanv(Channel)%GetL(PWC)
  end function ClassSymmetricElectronicSpaceGetL

  integer function ClassSymmetricElectronicSpaceGetM( Space, Channel, PWC ) result(M)
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    integer,                              intent(in)    :: Channel
    integer,                              intent(in)    :: PWC
    M = Space%Chanv(Channel)%GetM(PWC)
  end function ClassSymmetricElectronicSpaceGetM

  function ClassSymmetricElectronicSpaceGetChannelPWC( Space, iChannel, iPWC ) result(PWC)
    class(ClassSymmetricElectronicSpace), target, intent(in) :: Space
    integer,                                      intent(in)    :: iChannel
    integer,                                      intent(in)    :: iPWC
    type(ClassPartialWaveChannel), pointer :: PWC
    PWC => Space%Chanv(iChannel)%GetPWC(iPWC)
  end function ClassSymmetricElectronicSpaceGetChannelPWC

  function ClassSymmetricElectronicSpaceGetChannelPWCFromAbsoluteIndex( Space, iChannel ) result(PWC)
    class(ClassSymmetricElectronicSpace), target, intent(in) :: Space
    integer,                                      intent(in)    :: iChannel
    type(ClassPartialWaveChannel), pointer :: PWC
    integer :: Index1, Index2, m
    Index1 = 0
    Index2 = 0
    do m = 1, Space%GetNumChannels()
       Index2 = Index2 + Space%GetNumPWC(m)
       if ( (iChannel>=Index1) .and. (iChannel<=Index2) ) then
          PWC => Space%GetChannelPWC(m, iChannel-Index1)
          exit
       end if
       Index1 = Index2
    end do
  end function ClassSymmetricElectronicSpaceGetChannelPWCFromAbsoluteIndex

  function ClassSymmetricElectronicSpaceGetChannelSRC( Space, iChannel ) result(SRC)
    class(ClassSymmetricElectronicSpace), target, intent(in) :: Space
    integer,                                      intent(in)    :: iChannel
    type(ClassShortRangeChannel), pointer :: SRC
    SRC => Space%Chanv(iChannel)%GetSRC()
  end function ClassSymmetricElectronicSpaceGetChannelSRC

  function ClassSymmetricElectronicSpaceGetChannelLS( Space ) result(LS)
    class(ClassSymmetricElectronicSpace), target, intent(in) :: Space
    type(ClassSymmetricLocalizedStates), pointer :: LS
    LS => Space%LS
  end function ClassSymmetricElectronicSpaceGetChannelLS

  integer function ClassSymmetricElectronicSpaceGetNumPWC( Space, Channel ) result(NPWC)
    !
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer, optional,                    intent(in) :: Channel
    !
    integer :: i
    !
    if ( present(Channel) ) then
       !
       if ( (Channel<1) .or. (Channel>Space%NChannels) ) then
          call Assert( "The index of the channel is invalid, impossible to get the number of PWC orbitals." )
       end if
       !
       NPWC = Space%Chanv(Channel)%GetNumPWC()
       !
    else
       !
       NPWC = 0
       do i = 1, Space%NChannels
          NPWC = NPWC + Space%Chanv(i)%GetNumPWC()
       end do
       !
    end if
    !
  end function ClassSymmetricElectronicSpaceGetNumPWC

  integer function ClassSymmetricElectronicSpacePWCChannelIndex( Space, IndexPI, IndexXlm ) result(IndexPWC)
    !
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer,                              intent(in) :: IndexPI
    integer,                              intent(in) :: IndexXlm
    integer :: i
    !
    if ( (IndexPI<1) .or. (IndexPI>Space%NChannels) ) then
       call Assert( "The index of the channel is invalid, impossible to get the PWC index." )
    end if
    !
    IndexPWC = IndexXlm
    do i = 1, IndexPI-1
       IndexPWC = IndexPWC + Space%GetNumPWC(i)
    end do
    !
  end function ClassSymmetricElectronicSpacePWCChannelIndex

  real(kind(1d0)) function ClassSymmetricElectronicSpaceGetEnergyPI( Space, PIindex ) result(Energy)
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    !> Parent-ion index.
    integer,                              intent(in) :: PIindex
    Energy = Space%Chanv(PIindex)%GetEnergyPI()
  end function ClassSymmetricElectronicSpaceGetEnergyPI

  real(kind(1d0)) function ClassSymmetricElectronicSpaceGetChannelPIEnergy( Space, Channel ) result(Energy)
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer,                              intent(in) :: Channel
    integer :: Index1, Index2, m
    Index1 = 0
    Index2 = 0
    do m = 1, Space%GetNumChannels()
       Index2 = Index2 + Space%GetNumPWC(m)
       if ( (Channel>=Index1) .and. (Channel<=Index2) ) then
          Energy = Space%GetEnergyPI(m)
          exit
       end if
       Index1 = Index2
    end do
  end function ClassSymmetricElectronicSpaceGetChannelPIEnergy

  integer function ClassSymmetricElectronicSpaceGetChannelL( Space, Channel ) result(L)
    class(ClassSymmetricElectronicSpace), intent(in) :: Space
    integer,                              intent(in) :: Channel
    integer :: Index1, Index2, m
    Index1 = 0
    Index2 = 0
    do m = 1, Space%GetNumChannels()
       Index2 = Index2 + Space%GetNumPWC(m)
       if ( (Channel>=Index1) .and. (Channel<=Index2) ) then
          L = Space%GetL(m,Channel-Index1)
          exit
       end if
       Index1 = Index2
    end do
  end function ClassSymmetricElectronicSpaceGetChannelL
  
  subroutine ClassSymmetricElectronicSpaceUnsetSpin( Space )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    Space%Multiplicity = 0
  end subroutine ClassSymmetricElectronicSpaceUnsetSpin

  subroutine ClassSymmetricElectronicSpaceSetStorage( Space, StorageDir )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    character(len=*),                     intent(in)    :: StorageDir
    if ( allocated(Space%StorageDir) ) deallocate(Space%StorageDir)
    allocate( Space%StorageDir, source = StorageDir )
  end subroutine ClassSymmetricElectronicSpaceSetStorage

  subroutine ClassSymmetricElectronicSpaceSetNuclearLabel( Space, Label )
    class(ClassSymmetricElectronicSpace), intent(inout) :: Space
    character(len=*),                     intent(in)    :: Label
    if ( allocated(Space%NuclearLabel) ) deallocate( Space%NuclearLabel )
    allocate( Space%NuclearLabel, source = Label )
  end subroutine ClassSymmetricElectronicSpaceSetNuclearLabel

end module ModuleSymmetricElectronicSpace
