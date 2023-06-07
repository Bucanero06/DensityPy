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
!> Define the classes for the individual parent ions, for the set of
!! parent ions and for the multipoles between them.
!! Furthermore, it defines a global variable for the complete set of
!! parent ions as well as
!! - a configuration file for the parent ions
!! - the files that contain information on individual parent ions (from interface)
!! - the files with the mutipoles (from interface)
module ModuleParentIons

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleIO
  use ModuleErrorHandling
  use ModuleString
  use ModuleGroups
  use ModuleXlm

  implicit none

  private

  character(len=*), parameter :: PARENTION_DIR_NAME  = "ParentIons"
  character(len=*), parameter :: PIEnergyLabel       = "En_"

  !> The parent ions are identified by 
  !! - point group (which must be set once and for all at the very beginning)
  !! - multiplicity (if the spin is well defined, multiplicity >=1, &
  !!                 otherwise multiplicity =0)
  !! - irreducible representation, in the order provided by the list of the 
  !!   group class
  !! - number, which is assigned externally by the quantum chemistry program
  !!   and which does not necessarily reflect the energy order of the states
  !!   (due to the fact that the group used here may be just a sub-group of
  !!   the real molecular group, states with the same symmetry could actually
  !!   cross, thus altering the order of their energies).
  type, public :: ClassParentIon
     private

     real(kind(1d0))           :: Energy = 0.d0
     integer                   :: Charge = 0
     integer                   :: Multiplicity = 1
     type(ClassIrrep), pointer :: Irrep
     integer                   :: N = 0

   contains

     procedure, public :: init            => ClassParentIon_Init
     procedure, public :: free            => ClassParentIon_Free
     procedure, public :: show            => ClassParentIon_Show
     procedure, public :: SetEnergy       => ClassParentIon_SetEnergy
     procedure, public :: GetEnergy       => ClassParentIon_GetEnergy
     procedure, public :: GetCharge       => ClassParentIon_GetCharge
     procedure, public :: GetMultiplicity => ClassParentIon_GetMultiplicity
     procedure, public :: GetIrrep        => ClassParentIon_GetIrrep
     procedure, public :: GetNumber       => ClassParentIon_GetNumber
     procedure, public :: GetLabel        => ClassParentIon_GetLabel
     procedure, public :: ReadEnergy      => ClassParentIon_ReadEnergy
     generic  , public :: LoadData        => ClassParentIon_LoadData, &
          ClassParentIon_LoadAllDataFromUnit
     generic  , public :: SaveData        => ClassParentIon_SaveData, &
          ClassParentIon_SaveAllDataToUnit
     procedure, private :: ClassParentIon_LoadData
     procedure, private :: ClassParentIon_LoadAllDataFromUnit
     procedure, private :: ClassParentIon_SaveData
     procedure, private :: ClassParentIon_SaveAllDataToUnit
     final :: ClassParentIon_Final
  end type ClassParentIon

  public :: ParsePionLabel

contains
  
  subroutine ClassParentIon_Init( Self, Group, Label, iCharge, Energy )
    class(ClassParentIon)    , intent(inout) :: Self
    type(ClassGroup)         , intent(in)    :: Group
    character(len=*)         , intent(in)    :: Label
    integer                  , intent(in)    :: iCharge
    real(kind(1d0)), optional, intent(in)    :: Energy
    !
    integer :: Multiplicity
    character(len=8) :: IrrepLabel
    integer :: PionN
    if(present(Energy)) Self%Energy=Energy
    Self%Charge=iCharge
    call ParsePionLabel(Label,Multiplicity,IrrepLabel,PionN)
    if( Group%definesIrrep( trim(IrrepLabel) ) )then
       Self%Multiplicity =  Multiplicity
       Self%Irrep        => Group%GetIrrep( trim(IrrepLabel) )
       Self%N            =  PionN
    else
       call ErrorMessage("Invalid irrep "//trim(IrrepLabel)//" for group "//Group%GetName())
       stop
    endif
  end subroutine ClassParentIon_Init

  subroutine ClassParentIon_Show( Self, unit )
    class(ClassParentIon), intent(in) :: Self
    integer, optional    , intent(in) :: unit
    integer :: outunit
    character(len=8) :: strn
    outunit=OUTPUT_UNIT
    if(present(unit))outunit=unit
    write(outunit,"(a)") "Parent Ion info: "
    strn=Self%GetLabel()
    write(outunit,"(a)"      ,advance="no")"  "//strn
    write(outunit,"(a,1x,i0)",advance="no")"  Mult =",Self%GetMultiplicity() 
    write(outunit,"(a,1x,i0)",advance="no")"  Q =",Self%GetCharge() 
    write(outunit,"(a,1x,d14.6)")          "  E =",Self%GetEnergy()
  end subroutine ClassParentIon_Show

  function ClassParentIon_GetLabel( Self ) result( label )
    class(ClassParentIon), intent(in) :: Self
    character(len=:)    , allocatable :: label
    character(len=16) :: mstrn,nstrn,istrn
    write(mstrn,*)Self%GetMultiplicity()
    mstrn=adjustl(mstrn)
    istrn=adjustl(Self%Irrep%GetName())
    write(nstrn,*)Self%GetNumber()
    nstrn=adjustl(nstrn)
    allocate(label,source=trim(mstrn)//trim(istrn)//"."//trim(nstrn))
  end function ClassParentIon_GetLabel
    
  function ClassParentIon_GetMultiplicity( Self ) result( multiplicity )
    class(ClassParentIon), intent(in) :: Self
    integer                           :: multiplicity
    multiplicity = Self%multiplicity
  end function ClassParentIon_GetMultiplicity

  function ClassParentIon_GetIrrep( Self ) result( irrep )
    class(ClassParentIon), intent(in) :: Self
    type(ClassIrrep), pointer         :: irrep
    irrep => Self%irrep
  end function ClassParentIon_GetIrrep

  function ClassParentIon_GetNumber( Self ) result( number )
    class(ClassParentIon), intent(in) :: Self
    integer                           :: number
    number = Self%N
  end function ClassParentIon_GetNumber

  function ClassParentIon_GetEnergy( Self ) result( energy )
    class(ClassParentIon), intent(in) :: Self
    real(kind(1d0))                   :: energy
    energy = Self%Energy
  end function ClassParentIon_GetEnergy

  subroutine ClassParentIon_SetEnergy( Self, Energy )
    class(ClassParentIon), intent(inout) :: Self
    real(kind(1d0))      , intent(in)    :: Energy
    Self%Energy = Energy
  end subroutine ClassParentIon_SetEnergy

  function ClassParentIon_ReadEnergy( Self, Dir ) result( energy )
    class(ClassParentIon), intent(in) :: Self
    character(len=*)     , intent(in) :: Dir
    real(kind(1d0))                   :: energy
    character(len=:), allocatable :: FileName, StorageDir
    integer :: uid
    allocate( StorageDir, source=AddSlash(dir)//AddSlash(PARENTION_DIR_NAME) )
    allocate( FileName, source = StorageDir//PIEnergyLabel//Self%GetLabel() )
    call OpenFile( FileName, uid, 'read', 'formatted' )
    read(uid,*) energy
    close( uid )
  end function ClassParentIon_ReadEnergy

  function ClassParentIon_GetCharge( Self ) result( charge )
    class(ClassParentIon), intent(in) :: Self
    integer                           :: charge
    charge = Self%Charge
  end function ClassParentIon_GetCharge

  subroutine ClassParentIon_Free( Self )
    class(ClassParentIon) :: Self
    self%Energy = 0.d0
    self%Charge = 0
    self%Multiplicity = 0
    self%Irrep  => NULL()
    self%N      = 0
  end subroutine ClassParentIon_Free

  subroutine ClassParentIon_SaveData( Self, dir ) 
    class(ClassParentIon), intent(inout) :: Self
    character(len=*)     , intent(in)    :: dir
    character(len=:)     , allocatable   :: FileName, StorageDir
    integer                              :: uid, iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    allocate( StorageDir, source=AddSlash(dir)//AddSlash(PARENTION_DIR_NAME) )
    call Execute_Command_Line("mkdir -p "//StorageDir)
    allocate( FileName, source = StorageDir//PIEnergyLabel//Self%GetLabel() )
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="unknown",&
         action  ="write",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)call Assert(iomsg)
    write(uid,*) Self%Energy
    close(uid)
  end subroutine ClassParentIon_SaveData

  subroutine ClassParentIon_SaveAllDataToUnit( Self, uid ) 
    class(ClassParentIon), intent(inout) :: Self
    integer              , intent(in)    :: uid
    !
    integer                       :: iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    write(uid,*) Self%Energy
    write(uid,*) Self%Charge
    write(uid,*) Self%Multiplicity
    write(uid,*) Self%Irrep%GetGroupName()
    write(uid,*) Self%Irrep%GetName()
    write(uid,*) Self%N
    !
  end subroutine ClassParentIon_SaveAllDataToUnit

  subroutine ClassParentIon_LoadAllDataFromUnit( Self, uid ) 
    class(ClassParentIon), intent(inout) :: Self
    integer              , intent(in)    :: uid
    !
    integer                       :: iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    character(len=256) :: IrrepName, GroupName
    type(ClassGroup) :: Group
    !
    read(uid,*) Self%Energy
    read(uid,*) Self%Charge
    read(uid,*) Self%Multiplicity
    read(uid,*) GroupName
    read(uid,*) IrrepName
    read(uid,*) Self%N
    !
    call Group%Init( trim(GroupName) )
    if( Group%DefinesIrrep( trim(IrrepName) ) )then
       Self%Irrep => Group%GetIrrep( trim(IrrepName) )
    else
       call Assert("Invalid irrep "//trim(IrrepName)//" for group "//Group%GetName())
    end if
    !
  end subroutine ClassParentIon_LoadAllDataFromUnit

  subroutine ClassParentIon_LoadData( Self, dir ) 
    class(ClassParentIon), intent(inout) :: Self
    character(len=*)     , intent(in)    :: dir
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid, iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    allocate( FileName, source = AddSlash(dir)//&
         AddSlash(PARENTION_DIR_NAME)//PIEnergyLabel//Self%GetLabel())
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="old",&
         action  ="read",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Cannot Find File "//FileName)
       call Assert(iomsg)
    endif
    !
    read(uid,fmt=*,iostat=iostat) Self%Energy
    close(uid)
    !
  end subroutine ClassParentIon_LoadData

  subroutine ClassParentIon_Final( Self )
    type(ClassParentIon) :: Self
    call Self%free()
  end subroutine ClassParentIon_Final

  subroutine ParsePionLabel( Label, Multiplicity, Irrep, N )
    character(len=*) , intent(in)  :: Label
    integer          , intent(out) :: Multiplicity
    character(len=*) , intent(out) :: Irrep
    integer          , intent(out) :: N
    !
    integer :: iChar, iPos
    character(len=16) :: strn
    character(len=:),allocatable :: tmpLabel
    !
    allocate(tmpLabel,source=trim(adjustl(Label)))
    strn=" "
    iChar=1
    do 
       if(iChar>len_trim(tmpLabel))exit
       iPos=index("1234567890",tmpLabel(iChar:iChar))
       if(iPos<=0)exit
       iChar=iChar+1
    enddo
    iChar=iChar-1
    if(iChar<=0)call Assert("Missing multiplicity of Parent ion")

    strn=tmpLabel(1:iChar)
    read(strn,*)Multiplicity
    tmpLabel=tmpLabel(iChar+1:)
    !
    iPos=index(tmpLabel,".")
    Irrep=tmpLabel(1:iPos-1)
    strn=tmpLabel(iPos+1:)
    read(strn,*)N
    !
  end subroutine ParsePionLabel

end module ModuleParentIons
