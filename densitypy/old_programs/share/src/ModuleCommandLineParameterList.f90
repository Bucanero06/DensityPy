!! CONFIDENTIAL
!! ModuleCommandLineParameterList, Copyright (C) 2020 by Luca Argenti, PhD - Some Rights Reserved
!! ModuleCommandLineParameterList is licensed under a
!! Creative Commons Attribution-ShareAlike 4.0 International License.
!! A copy of the license is available at <http://creativecommons.org/licenses/by-nd/4.0/>.
!!
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!! email: luca.argenti@ucf.edu

! {{{ Detailed description

!> \file
!!
!! Defines the package of classes for reading
!! run time parameters from the command line.
!!
! }}}
module ModuleCommandLineParameterList

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none

  private

  !> Creates a list with all the parameters to be read from the command line.
  type, public ::  ClassCommandLineParameterList
     !
     private
     !
     !> Stores detailed description about the parameters should be written in the command line.
     character(len=:)    , allocatable :: Description
     !> Points to the first command line parameter available in the analized interval.
     type(ClassParameter), pointer     :: First => NULL()
     !> Points to the last command line parameter available in the analized interval.
     type(ClassParameter), pointer     :: Last  => NULL()
     !
   contains
     !
     !> Sets the ClassCommandLineParameterList's Description variable, from an external string.
     procedure, public  :: SetDescription
     !> Prints the program usage.
     procedure, public  :: PrintUsage
     !> Adds any parameter to the list of run time parameters to be read.
     generic  , public  :: Add     => AddSwitchParameter, AddValuedParameter    
     !> Parses the command line.
     procedure, public  :: Parse   => ParseCommandLine
     !> Retrieves whether a parameter is present in the command line or not.
     procedure, public  :: Present => ParameterIsPresent
     !> Gets the parameter value.
     procedure, public  :: Get     => GetParameterFromList
     !> Frees ClassCommandLineParameterList.
     procedure, public  :: Free    => FreeParameterList 
     !> Prints on a unit all the parameters.
     procedure, public  :: PrintAll => ParameterListPrint
     !
     procedure, private :: AddValuedParameter
     procedure, private :: AddSwitchParameter
     procedure, private :: WhichParameter
     procedure, private :: Check  
     procedure, private :: Insert
     !
  end type ClassCommandLineParameterList
 

  !> Stores the properties of any run time parameter.
  type, private :: ClassParameter
     !
     !> Points to the previous run time parameter.
     type(ClassParameter), pointer     :: Prev => NULL()
     !> Points to the next run time parameter.
     type(ClassParameter), pointer     :: Next => NULL()
     !
     !> Run time parameter's name.
     character(len=:)    , allocatable :: Name
     !> Run time parameter's purpose.
     character(len=:)    , allocatable :: Purpose
     !> Run time parameter's value.
     class(*)            , allocatable :: Value
     !
     !> Indicates whether the parameter is required for the program execution or not.
     logical                           :: IsRequired = .FALSE.
     !> Indicates whether the parameter is present in the command line or not.
     logical                           :: IsPresent  = .FALSE.
     !> Indicates whether the parameter has a value or not.
     logical                           :: IsValued   = .FALSE.
     !
   contains
     !
     !> Gets the value of the parameter.
     procedure :: Get   => GetFromParameter
     !> Sets the ClassParameter objects corresponding to the requested parameter.
     procedure :: Set   => ParameterSet 
     !> Print the usage of a run time parameter on the screen, e.g.:
     !! 
     !! [--help (print usage)]
     !!  -ne <number of energies>
     !!
     procedure :: Print => ParameterPrint
     !> Frees the ClassParameter objects.
     procedure :: Free  => FreeParameter 
     !> Converts the parameter value from string format to the appropiated one (integer, double precision, logical  or character).
     procedure :: StringToValue          
     !> Gets the kind of the parameter: "Integer", "Double", "Logical", "String" or " " for other kind. 
     procedure :: Kind  => ParameterKind 
     !> Gets the parameter default value in string format.
     procedure :: DefaultValueString  => ParameterDefaultValueString
     !
  end type ClassParameter


contains


  !> Sets the ClassCommandLineParameterList's Description variable, from an external string.
  subroutine SetDescription( List, Description )
    class(ClassCommandLineParameterList), intent(out) :: List
    character(len=*)                    , intent(in)  :: Description
    if(allocated(List%Description))deallocate(List%Description)
    allocate(List%Description,source=trim(adjustl(Description)))
  end subroutine SetDescription


  subroutine AddSwitchParameter( List, Name, Purpose )
    class(ClassCommandLineParameterList), intent(inout) :: List
    character(len=*)                    , intent(in)    :: Name
    character(len=*)                    , intent(in)    :: Purpose
    !
    logical, parameter :: IS_VALUED = .FALSE.
    class(ClassParameter), pointer :: Parameter
    allocate(Parameter)
    call Parameter%Set(Name,"optional",Purpose,IS_VALUED)
    call List%Insert(Parameter)
  end subroutine AddSwitchParameter


  subroutine AddValuedParameter( List, Name, Purpose, Value, Condition )
    class(ClassCommandLineParameterList), intent(inout) :: List
    character(len=*)                    , intent(in)    :: Name
    character(len=*)                    , intent(in)    :: Purpose
    Class(*)                            , intent(in)    :: Value
    character(len=*)                    , intent(in)    :: Condition
    !
    logical, parameter :: IS_VALUED = .TRUE.
    class(ClassParameter), pointer :: Parameter
    Parameter => CreateParameter( Value )
    call Parameter%Set(Name,Condition,Purpose,IS_VALUED)
    call List%Insert(Parameter)
  end subroutine AddValuedParameter


  function CreateParameter( Value ) result(Parameter)
    class(*), intent(in) :: Value
    class(ClassParameter), pointer :: Parameter
    allocate(Parameter)
    allocate(Parameter%value,source=Value)
  end function CreateParameter


  !> Sets the ClassParameter objects corresponding to the requested parameter.
  subroutine ParameterSet( Parameter, Name, Condition, Purpose, IsValued )
    class(classParameter), intent(inout) :: Parameter
    character(len=*)     , intent(in)    :: Name
    character(len=*)     , intent(in)    :: Condition
    character(len=*)     , intent(in)    :: Purpose
    logical              , intent(in)    :: IsValued
    Parameter%Prev => Null()
    Parameter%Next => Null()
    allocate(Parameter%Name,source=trim(adjustl(Name)))
    Parameter%IsPresent=.FALSE.
    select case( Condition )
    case("optional")
       Parameter%IsRequired=.FALSE.
    case("required")
       Parameter%IsRequired=.TRUE.
    case DEFAULT
       call Assert("Invalid parameter condition")
    end select
    allocate( Parameter%Purpose, source = trim( adjustl( Purpose ) ) )
    Parameter%IsValued = IsValued
  end subroutine ParameterSet


  subroutine Insert( List, Parameter )
    class(ClassCommandLineParameterList), intent(inout) :: List
    class(ClassParameter), pointer      , intent(in)    :: Parameter
    if( .not.associated(List%First) )then
       List%First => Parameter
       List%Last  => Parameter
    else
       Parameter%Prev => List%Last
       List%Last%Next => Parameter
       List%Last      => Parameter
    endif
  end subroutine Insert


  !> Print the usage of a run time parameter on the screen, e.g.:
  !! 
  !! [--help (print usage)]
  !!  -ne <number of energies>
  !!
  subroutine ParameterPrint( Parameter, uid )
    class(ClassParameter), intent(in) :: Parameter
    integer              , intent(in) :: uid
    character(len=512) :: Synopsis
    Synopsis=" "
    Synopsis=trim(Synopsis)//trim(Parameter%Name)
    if( Parameter%IsValued )then
       Synopsis=trim(Synopsis)//" <"//Parameter%Purpose//"> "
       if(.not.Parameter%IsRequired)then
          Synopsis=trim(Synopsis)//", default="//Parameter%DefaultValueString()
       endif
       Synopsis=trim(Synopsis)//" ("//Parameter%Kind()//") "
    else
       Synopsis=trim(Synopsis)//" ("//Parameter%Purpose//") "
    endif
    if( Parameter%IsRequired )then
       Synopsis=" "//trim(Synopsis)
    else
       Synopsis="["//trim(Synopsis)//"]"
    endif
    write(uid,"(a)") trim(Synopsis)
  end subroutine ParameterPrint


  !> Gets the kind of the parameter: "Integer", "Double", "Logical", "String" or " " for other kind. 
  function ParameterKind( Parameter ) result( KindStrn )
    Class(ClassParameter), intent(in) :: Parameter
    character(len=:), allocatable     :: KindStrn
    select type(ptr=>Parameter%Value)
    type is(integer)
       allocate(KindStrn,source="Integer")
    type is(real(kind(1d0)))
       allocate(KindStrn,source="Double")
    type is(logical)
       allocate(KindStrn,source="Logical")
    type is(character(len=*))
       allocate(KindStrn,source="String")
    class DEFAULT
       allocate(KindStrn,source=" ")
    end select
  end function ParameterKind


  !> Gets the parameter default value in string format.
  function ParameterDefaultValueString( Parameter ) result( DefaultValueStrn )
    Class(ClassParameter), intent(in) :: Parameter
    character(len=:), allocatable     :: DefaultValueStrn
    character(len=512) :: strn
    if(.not.allocated(Parameter%Value))then
       strn = "[Not Valued]"
    else
       select type(ptr=>Parameter%Value)
       type is(integer)
          write(strn,"(i0)") ptr
       type is(real(kind(1d0)))
          write(strn,"(d11.3)") ptr
       type is(logical)
          write(strn,*) ptr
       type is(character(len=*))
          strn=trim(adjustl(ptr))
       class DEFAULT
          strn = "[Unrecognized Value Kind]"
       end select
    endif
    allocate(DefaultValueStrn,source=trim(adjustl(strn)))
  end function ParameterDefaultValueString


  !> Prints the program usage.
  subroutine PrintUsage( List, OutputUnit )
    !
    class(ClassCommandLineParameterList), intent(in) :: List
    integer, optional                   , intent(in) :: OutputUnit
    character(len=64) :: ProgramName
    integer :: stat, LastSlash, uid
    !
    !> Determines the current name of the executable
    call Get_Command_Argument(0,ProgramName,status=stat)
    if(stat/=0)ProgramName="<Program Name>"
    LastSlash=index(ProgramName,"/",back=.true.)
    ProgramName=adjustl(ProgramName(LastSlash+1:))
    !
    uid = OUTPUT_UNIT
    if(present(OutputUnit))uid=OutputUnit
    !> Print the underlined name of the program followed 
    !! by its description 
    write(uid,"(a)") 
    write(uid,"(a)") trim(ProgramName)
    write(uid,"(a)") repeat("=",Len_Trim(ProgramName))
    write(uid,"(a)") trim(List%Description)
    write(uid,"(a)") 
    call List%PrintAll( uid )
    write(uid,"(a)") 
    !
  end subroutine PrintUsage


  subroutine ParameterListPrint( List, uid )
    class(ClassCommandLineParameterList), intent(in) :: List
    integer                             , intent(in) :: uid
    class(ClassParameter), pointer :: Parameter
    Parameter=>List%First
    do while(associated(Parameter))
       call Parameter%Print(uid)
       Parameter=>Parameter%Next
    enddo
  end subroutine ParameterListPrint


  !> Gets the parameter value.
  subroutine GetParameterFromList( List, Name, Value )
    Class(ClassCommandLineParameterList), intent(in) :: List
    character(len=*)                    , intent(in) :: Name
    class(*)                                         :: Value
    Class(ClassParameter), pointer :: Parameter
    Parameter => List%WhichParameter( trim( Name ) )
    if(.not.Associated(Parameter)) call Assert("Unrecognized Parameter "//trim(Name))
    call Parameter%Get(Value)
  end subroutine GetParameterFromList


  function WhichParameter( List, Name ) result( Parameter )
    Class(ClassCommandLineParameterList), intent(in) :: List
    character(len=*)                    , intent(in) :: Name
    Class(ClassParameter)               , pointer    :: Parameter
    Parameter => List%First
    do while( associated( Parameter ) )
       if(trim(Parameter%Name)==trim(Name))return
       Parameter => Parameter%Next
    enddo
  end function WhichParameter

  !> Gets the value of the parameter.
  subroutine GetFromParameter( Parameter, Value )
    Class(ClassParameter), intent(in) :: Parameter
    class(*)                          :: Value
    !.. I would have liked sooo much that the following
    !   worked. But unfortunately it doesn't.
    !if(SAME_TYPE_AS(Parameter%value,value))then
    !   value=Parameter%value
    !end if
    select type(ptr=>Parameter%value)
    type is (Integer)
       select type(value)
       type is (Integer)
          value=ptr
       class DEFAULT
          call Assert("Invalid type request")
       end select
    type is (real(kind(1d0)))
       select type(value)
       type is (real(kind(1d0)))
          value=ptr
       class DEFAULT
          call Assert("Invalid type request")
       end select
    type is (logical)
       select type(value)
       type is (logical)
          value=ptr
       class DEFAULT
          call Assert("Invalid type request")
       end select
    type is (character(len=*))
       select type(value)
       type is (character(len=*))
          value=trim(ptr)
       class DEFAULT
          call Assert("Invalid type request")
       end select
    class DEFAULT
       call Assert("Unknown type")
    end select
  end subroutine GetFromParameter


  !> Parses the command line.
  subroutine ParseCommandLine( List )
    !
    class(ClassCommandLineParameterList), intent(inout) :: List
    !
    integer, parameter :: MAX_COMMAND_LINE_LENGTH = 2000
    character(len=*), parameter :: HERE = ":ParseCommandLine:"
    character(len=MAX_COMMAND_LINE_LENGTH) :: CommandLine
    character(len=MAX_COMMAND_LINE_LENGTH) :: ParameterLine
    integer :: status
    integer :: EndOfExecutableName
    type(ClassParameter), pointer :: Parameter
    integer :: ParameterPosition
    integer :: EndOfParameterName
    integer :: EndOfParameterSpec
    integer :: iostat
    logical :: SUCCESS

    SUCCESS = .TRUE.

    !.. Read the command line
    call Get_Command( CommandLine, status = status )
    if(status/=0)call ASSERT(HERE//" Internal Error")
    
    !.. Purge the name of the command from the command line
    CommandLine=adjustl(CommandLine)
    EndOfExecutableName=index(CommandLine," ")
    CommandLine=CommandLine(EndOfExecutableName:)

    !.. Cycle over Formal Run Time Parameters
    Parameter => List%First
    list_scan : do while( associated( Parameter ) )
       
       !.. Search for the parameter in the command line.
       !   The spaces around the name are necessary to distinguish 
       !   " -n 100" from " -next 2" and " -strn fool-name "
       ParameterPosition = index( CommandLine, " "//trim(Parameter%Name)//" " )
       Parameter%IsPresent = ParameterPosition > 0
       if( Parameter%IsPresent )then
          ParameterPosition = ParameterPosition + 1
       else
          !.. If the parameter is absent but required, issue a warning and stop,
          !   otherwise update the pointer and move on
          if( Parameter%IsRequired )then
             !call ErrorMessage("Parameter "//trim(Parameter%Name)//" is required")
             SUCCESS = .FALSE.
          endif
          Parameter => Parameter%Next
          cycle list_scan
       endif
       
       !.. If the parameter is valued, read the value.
       if( Parameter%IsValued )then
          !
          !.. Extract Parameter Specification
          ParameterLine = adjustl(CommandLine(ParameterPosition:))
          !
          !.. Purge Parameter Name 
          EndOfParameterName = Index(ParameterLine," ")
          ParameterLine = adjustl(ParameterLine(EndOfParameterName:))
          !
          !.. Extract the string that supposedly contains the parameter value
          EndOfParameterSpec = Index(ParameterLine," ")
          if( EndOfParameterSpec <= 0 )then
             EndOfParameterSpec=len_trim(ParameterLine)
          else
             EndOfParameterSpec=EndOfParameterSpec-1
          endif
          ParameterLine(EndOfParameterSpec+1:)=" "
          !
          !.. Assign the value 
          call Parameter%StringToValue( trim(ParameterLine), iostat )
          SUCCESS = SUCCESS .and. ( iostat == 0 )
          !
       endif
       
       Parameter => Parameter%Next
       cycle list_scan

    enddo list_scan

    call List%Check(iostat)
    SUCCESS = SUCCESS .and. ( iostat == 0 )

    if( .not. SUCCESS )then
       call List%PrintUsage( OUTPUT_UNIT )
       STOP
    endif

  end subroutine ParseCommandLine


  !> Converts the parameter value from string format to the appropiated one (integer, double precision, logical  or character).
  subroutine StringToValue( Parameter, ValueStrn, iostat )
    class(ClassParameter), intent(inout) :: Parameter
    character(len=*)     , intent(in)    :: ValueStrn
    integer              , intent(out)   :: iostat
    character(len=512) :: iomsg
    iostat=0
    iomsg=" "
    select type (ptr=>Parameter%value)
    type is(integer)
       read(ValueStrn,*,iostat=iostat,iomsg=iomsg) ptr
    type is(real(kind(1d0)))
       read(ValueStrn,*,iostat=iostat,iomsg=iomsg) ptr
    type is(logical)
       read(ValueStrn,*,iostat=iostat,iomsg=iomsg) ptr
    type is(character(len=*))
       deallocate(Parameter%value)
       allocate(Parameter%value,source=ValueStrn)
    class DEFAULT
       call ErrorMessage("Non-standard type")
    end select
    if(iostat/=0)then
       call ErrorMessage("Invalid parameter "//trim(Parameter%Name))
    endif
    Parameter%IsPresent=.TRUE.
  end subroutine StringToValue


  !> Retrieves whether a parameter is present in the command line or not.
  logical function ParameterIsPresent( List, Name ) result( Present )
    class(ClassCommandLineParameterList), intent(in) :: List
    character(len=*)                    , intent(in) :: Name
    class(ClassParameter), pointer :: Parameter
    Parameter=> List%WhichParameter(Name)
    Present=.FALSE.
    if(Associated(Parameter)) Present=Parameter%IsPresent
  end function ParameterIsPresent


  !> Frees the ClassParameter objects.
  subroutine FreeParameter( Parameter )
    class( ClassParameter ), intent(inout) :: Parameter
    if(allocated(Parameter%Name   ))deallocate(Parameter%Name)
    if(allocated(Parameter%Purpose))deallocate(Parameter%Purpose)
    if(allocated(Parameter%Value  ))deallocate(Parameter%Value)
    Parameter%IsRequired=.FALSE.
    Parameter%IsPresent =.FALSE.
    Parameter%IsValued  =.FALSE.
    Parameter%Prev=>NULL()
    Parameter%Next=>NULL()
  end subroutine FreeParameter
  

  !> Frees ClassCommandLineParameterList.
  subroutine FreeParameterList( List )
    class( ClassCommandLineParameterList ), intent(inout) :: List
    type( ClassParameter ), pointer :: Parameter, Next
    Parameter => List%First
    do while(associated(Parameter))
       Next => Parameter%Next
       call Parameter%Free()
       deallocate(Parameter)
       Parameter => Next
    enddo
    List%First => NULL()
    List%Last  => NULL()
  end subroutine FreeParameterList


  subroutine Check( List, iostat )
    Class(ClassCommandLineParameterList), intent(in) :: List
    integer                             , intent(out):: iostat
    Class(ClassParameter), pointer :: Parameter
    iostat=0
    Parameter => List%First
    do while( associated( Parameter ) )
       if( Parameter%IsRequired .and. .not. Parameter%IsPresent )then
          call ErrorMessage("Required Parameter "//trim(Parameter%Name)//" is missing.")
          iostat=-1
       endif
       Parameter => Parameter%Next
    enddo
  end subroutine Check


  subroutine ASSERT(Message)
    use ISO_FORTRAN_ENV
    character(len=*), intent(in) :: Message
    write(OUTPUT_UNIT,"(a)")trim(Message)
    stop
  end subroutine ASSERT


  subroutine ErrorMessage(Message)
    use ISO_FORTRAN_ENV
    character(len=*), intent(in) :: Message
    write(OUTPUT_UNIT,"(a)")trim(Message)
  end subroutine ErrorMessage


end module ModuleCommandLineParameterList
