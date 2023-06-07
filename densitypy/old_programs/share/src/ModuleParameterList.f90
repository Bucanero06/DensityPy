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
! {{{ Detailed information
!> \file
!!
!! Module containing the neccesary classes that manage the configuration files reading.
! }}}

module ModuleParameterList


  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling

  
  private


  type, public :: ClassParameterList
     !
     !> Points to the first parameter in a list.
     Class(ClassParameter), pointer, private :: First => NULL()
     !> Points to the last parameter in a list.
     Class(ClassParameter), pointer, private :: Last  => NULL()
     !
   contains
     !
     !> Adds a parameter to the parameter list.
     procedure, public  :: Add     => AddParameter
     !> Parses a file either passing its name or its unit number.
     generic  , public  :: Parse     => ParseFile, ParseUnit
     generic  , public  :: ParseStrn => ParseString
     !> Gets the parameter value from a list of parameters.
     procedure, public  :: Get     => GetParameterFromList
     !> Retrieves whether a requested parameter is present or not in a list.
     procedure, public  :: Present => ParameterIsPresent
     !> Prints the parameter list on a unit.
     procedure, public  :: Print   => ParameterListPrint
     !> Prints the parameter list on a unit.
     procedure, public  :: WriteDefault => ParameterListWriteDefault
     !> Deallocates the parameter list.
     procedure, public  :: Free    => FreeParameterList
     !
     procedure, private :: ParseFile
     procedure, private :: ParseUnit
     procedure, private :: ParseString
     procedure, private :: WhichParameter
     procedure, private :: Check 
     procedure, private :: Insert
     !
     !> Deallocates the ClassParameterList.
     final :: FinalizeParameterList
     !
  end type ClassParameterList


  type, private :: ClassParameter
     !
     class(ClassParameter), pointer     :: Prev => NULL()
     class(ClassParameter), pointer     :: Next => NULL()
     character(len=:)     , allocatable :: Name
     class(*)             , allocatable :: value
     logical                            :: IsRequired = .FALSE.
     logical                            :: IsAbsent   = .TRUE.
     !
   contains
     !
     procedure :: Get   => GetFromParameter
     procedure :: Set   => ParameterSet
     procedure :: Print => ParameterPrint
     procedure :: PrintDefault => ParameterPrintDefault
     procedure :: Free  => FreeParameter
     procedure :: StringToValue
     !
  end type ClassParameter
  

contains


  function CreateParameter( Value ) result(Parameter)
    class(*), intent(in) :: Value
    class(ClassParameter), pointer :: Parameter
    allocate(Parameter)
    allocate(Parameter%value,source=Value)
  end function CreateParameter


  subroutine AddParameter( List, Name, Value, Condition )
    class(ClassParameterList), intent(inout) :: List
    character(len=*)         , intent(in)    :: Name
    character(len=*)         , intent(in)    :: Condition
    Class(*)                 , intent(in)    :: Value
    class(ClassParameter), pointer :: Parameter
    Parameter => CreateParameter( Value )
    call Parameter%Set(Name,Condition)
    call List%Insert(Parameter)
  end subroutine AddParameter


  subroutine Insert( List, Parameter )
    class(ClassParameterList)     , intent(inout) :: List
    class(ClassParameter), pointer, intent(in)    :: Parameter
    if( .not.associated(List%First) )then
       List%First => Parameter
       List%Last  => Parameter
    else
       Parameter%Prev => List%Last
       List%Last%Next => Parameter
       List%Last      => Parameter
    endif
  end subroutine Insert
  

  subroutine ParameterSet( Parameter, Name, Condition )
    class(classParameter), intent(inout) :: Parameter
    character(len=*)     , intent(in)    :: Name
    character(len=*)     , intent(in)    :: Condition
    Parameter%Prev => Null()
    Parameter%Next => Null()
    allocate(Parameter%Name,source=trim(adjustl(Name)))
    Parameter%IsAbsent=.TRUE.
    select case( Condition )
    case("optional")
       Parameter%IsRequired=.FALSE.
    case("required")
       Parameter%IsRequired=.TRUE.
    case DEFAULT
       call Assert("Invalid parameter condition")
    end select
  end subroutine ParameterSet


  subroutine ParameterPrint( Parameter, uid )
    class(ClassParameter), intent(in) :: Parameter
    integer              , intent(in) :: uid
    write(uid,*)trim(Parameter%Name)
    write(uid,*)Parameter%IsRequired
    write(uid,*)Parameter%IsAbsent
    select type(ptr=>Parameter%Value)
    type is(integer)
       write(uid,*)ptr
    type is(real(kind(1d0)))
       write(uid,*)ptr
    type is(logical)
       write(uid,*)ptr
    type is(character(len=*))
       write(uid,*,DELIM="QUOTE")ptr
    class DEFAULT
       call ErrorMessage("Non-standard type")
    end select
  end subroutine ParameterPrint

  subroutine ParameterPrintDefault( Parameter, uid )
    class(ClassParameter), intent(in) :: Parameter
    integer              , intent(in) :: uid
    write(uid,"(a)",advance="no")trim(Parameter%Name)//"="
    select type(ptr=>Parameter%Value)
    type is(integer)
       write(uid,*)ptr
    type is(real(kind(1d0)))
       write(uid,*)ptr
    type is(logical)
       write(uid,*)ptr
    type is(character(len=*))
       write(uid,"(a)")ptr
    class DEFAULT
    end select
  end subroutine ParameterPrintDefault


  subroutine ParameterListWriteDefault( List, File )
    class(ClassParameterList), intent(in) :: List
    character(len=*)         , intent(in) :: File
    integer :: uid, iostat
    character(len=100) :: iomsg
    class(ClassParameter), pointer :: Parameter
    open(newunit= uid       , &
         file   = File      , &
         form   ="formatted", &
         status ="unknown"  , &
         action ="write"    , &
         iostat = iostat    , &
         iomsg  = iomsg     )
    if(iostat/=0)then
       call ErrorMessage(trim(iomsg))
       stop
    endif
    Parameter=>List%First
    do while(associated(Parameter))
       call Parameter%PrintDefault(uid)
       Parameter=>Parameter%Next
    enddo
    close(uid)
  end subroutine ParameterListWriteDefault

  
  subroutine ParameterListPrint( List, uid )
    class(ClassParameterList), intent(in) :: List
    integer                  , intent(in) :: uid
    class(ClassParameter), pointer :: Parameter
    Parameter=>List%First
    do while(associated(Parameter))
       call Parameter%Print(uid)
       Parameter=>Parameter%Next
    enddo
  end subroutine ParameterListPrint


  subroutine GetParameterFromList( List, Name, Value )
    Class(ClassParameterList), intent(in) :: List
    character(len=*)         , intent(in) :: Name
    class(*)                              :: Value
    Class(ClassParameter), pointer :: Parameter
    Parameter => List%WhichParameter( trim( Name ) )
    if(.not.Associated(Parameter)) call Assert("Unrecognized Parameter "//trim(Name))
    call GetFromParameter(Parameter,Value)
  end subroutine GetParameterFromList


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


  subroutine ParseUnit( List, uid )
    !
    class(ClassParameterList), intent(inout) :: List
    integer                  , intent(in)    :: uid
    !
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Readable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    character(len=500):: FileName
    !
    character(len=2000) :: Line
    integer             :: lineNumber, ichar
    character(len=16)   :: LineNumberStrn
    character(len=:), allocatable  :: NameStrn
    character(len=:), allocatable  :: ValueStrn
    class(ClassParameter), pointer :: Parameter
    !
    INQUIRE(&
         UNIT  = uid     , &
         OPENED= Opened  , &
         READ  = Readable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg   , &
         NAME  = FileName  &
         )
    if(iostat/=0)then
       call ErrorMessage("Unreadable Configuration File "//trim(FileName))
       call Assert(iomsg)
    endif
    if( .not. Opened            ) call Assert("Input Unit is not open")
    if( trim(Readable) /= "YES" ) call Assert("Input Unit can't be read")
    if( trim(Form)/="FORMATTED" ) call Assert("Unformatted Configuration File "//trim(FileName))

    lineNumber=0
    scanFile: do

       read(uid,"(a)",iostat=iostat,iomsg=iomsg)Line
       if(iostat/=0)exit scanFile
       lineNumber=lineNumber+1
       write(LineNumberStrn,*)lineNumber

       !.. Cycle if the comment-less line is empty
       ichar=index(Line,"#")
       if(ichar>0)Line(ichar:)=" "
       Line=adjustl(Line)
       if(len_trim(Line)==0)cycle scanFile

       ichar=index(Line,"=")
       if(ichar<=0)then
          call ErrorMessage("Sintax Error on line "//&
               trim(LineNumberStrn)//"in "//trim(FileName))
          call ErrorMessage("'=' sign required for assignment")
          cycle scanFile
       endif

       if(allocated(NameStrn))deallocate(NameStrn)
       allocate(NameStrn,source=adjustl(Line(:ichar-1)))
       if(allocated(ValueStrn))deallocate(ValueStrn)
       allocate(ValueStrn,source=adjustl(Line(ichar+1:)))

       Parameter => List%WhichParameter( trim( NameStrn ) )
       if(.not.Associated(Parameter))then
          call ErrorMessage("Unrecognized Parameter "//trim(NameStrn))
          cycle scanFile
       endif

       call Parameter%StringToValue( trim(ValueStrn) )

    enddo scanFile

    call List%Check()
    if(allocated(NameStrn))deallocate(NameStrn)
    if(allocated(ValueStrn))deallocate(ValueStrn)

  end subroutine ParseUnit


  subroutine ParseString( List, Strn )
    !
    class(ClassParameterList), intent(inout) :: List
    character(len=*)         , intent(in)    :: Strn
    !
    integer           :: iostat
    !
    character(len=2000) :: Line
    integer             :: ichar, ich2
    character(len=16)   :: LineNumberStrn
    character(len=:), allocatable  :: NameStrn
    character(len=:), allocatable  :: ValueStrn
    class(ClassParameter), pointer :: Parameter
    !
    Line=adjustl(Strn)
    ichar=index(Line,"#")
    if(ichar>0)Line(ichar:)=" "
    Line=adjustl(Line)
    !
    scanStrn: do
       !
       !.. Exit if the line is empty 
       if(len_trim(Line)==0)exit scanStrn

       !.. If there are no valid assignments, issue a warning and exit
       ichar=index(Line,"=")
       if(ichar<=0)then
          call ErrorMessage("Sintax Error in String parsing: '=' sign required")
          exit scanStrn
       elseif(ichar==1)then
          call ErrorMessage("Sintax Error in String parsing: Parameter name required")
          exit scanStrn
       endif

       if(allocated(NameStrn))deallocate(NameStrn)
       allocate(NameStrn,source=trim(adjustl(Line(:ichar-1))))
       Line=adjustl(Line(ichar+1:))

       ichar=index(Line," ") 
       if(ichar==1)then
          call ErrorMessage("Sintax Error in String parsing: a value is required for "//NameStrn)
          exit scanStrn
       endif
      
       if(allocated(ValueStrn))deallocate(ValueStrn)
       allocate(ValueStrn,source=adjustl(Line(1:ichar-1)))
       Line=adjustl(Line(ichar:))

       Parameter => List%WhichParameter( trim( NameStrn ) )
       if(.not.Associated(Parameter))then
          call ErrorMessage("Unrecognized Parameter "//trim(NameStrn))
          cycle scanStrn
       endif

       call Parameter%StringToValue( trim(ValueStrn) )

    enddo scanStrn

    call List%Check()
    if(allocated(NameStrn))deallocate(NameStrn)
    if(allocated(ValueStrn))deallocate(ValueStrn)

  end subroutine ParseString


  subroutine ParseFile( List, FileName )
    class(ClassParameterList), intent(inout) :: List
    character(len=*)         , intent(in)    :: FileName
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    open(newunit = uid        , &
         file    = FileName   , &
         action  = "Read"     , &
         form    = "Formatted", &
         status  = "Old"      , &
         iostat  = iostat     , &
         iomsg   = iomsg      )
    if(iostat/=0)then
       call ErrorMessage("File "//trim(FileName)//" is not accessible")
       call ErrorMessage("iomsg ="//trim(iomsg))
       stop
    endif
    call ParseUnit( List, uid )
  end subroutine ParseFile


  function WhichParameter( List, Name ) result( Parameter )
    Class(ClassParameterList), intent(in) :: List
    character(len=*)         , intent(in) :: Name
    Class(ClassParameter)    , pointer    :: Parameter
    Parameter => List%First
    do while( associated( Parameter ) )
       if(trim(Parameter%Name)==trim(Name))return
       Parameter => Parameter%Next
    enddo
  end function WhichParameter


  subroutine StringToValue( Parameter, ValueStrn )
    class(ClassParameter), intent(inout) :: Parameter
    character(len=*)                     :: ValueStrn
    integer :: iostat
    character(len=IOMSG_LENGTH) :: iomsg
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
       STOP
    endif
    Parameter%IsAbsent=.FALSE.
  end subroutine StringToValue


  subroutine Check( List )
    Class(ClassParameterList), intent(in) :: List
    Class(ClassParameter), pointer :: Parameter
    Parameter => List%First
    do while( associated( Parameter ) )
       if( Parameter%IsRequired .and. Parameter%IsAbsent )then
          call ErrorMessage("Required Parameter "//trim(Parameter%Name)//" is missing.")
          STOP
       endif
       Parameter => Parameter%Next
    enddo
  end subroutine Check


  logical function ParameterIsPresent( List, Name ) result( Present )
    class(ClassParameterList), intent(in) :: List
    character(len=*)         , intent(in) :: Name
    class(ClassParameter), pointer :: Parameter
    Parameter=> List%WhichParameter(Name)
    Present=.FALSE.
    if(Associated(Parameter)) Present=.not.Parameter%IsAbsent
  end function ParameterIsPresent


  subroutine FreeParameter( Parameter )
    class( ClassParameter ), intent(inout) :: Parameter
    if(allocated(Parameter%Name))deallocate(Parameter%Name)
    Parameter%IsRequired=.FALSE.
    Parameter%IsAbsent=.TRUE.
    Parameter%Prev=>NULL()
    Parameter%Next=>NULL()
    if(allocated(Parameter%Value))deallocate(Parameter%Value)
  end subroutine FreeParameter
  

  subroutine FreeParameterList( List )
    class( ClassParameterList ), intent(inout) :: List
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


  subroutine FinalizeParameterList( List )
    type( ClassParameterList ), intent(inout) :: List
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
  end subroutine FinalizeParameterList


end module ModuleParameterList
