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
!> \file
!!
!! Provides assertion and error handling functions.
!!
module ModuleErrorHandling


    !The ModuleErrorHandling in Fortran provides a set of functionalities to handle errors in your program. It defines several assertion levels and provides the ability to change these levels at runtime.
    !
    !Here's an overview of the functionalities in the module:
    !    Definition of Assertion Levels: Different levels of assertions are defined with corresponding names. They are represented by instances of the ClassAssertionLevel type. They include LEVEL_WARNING, LEVEL_SEVERE, LEVEL_FATAL, and LEVEL_DEFAULT, with LEVEL_DEFAULT initially set to LEVEL_FATAL.
    !    Assertion Handling: The assert subroutine is provided for issuing assertions at different levels and to different output units. Depending on the level of the assertion, it can either continue execution with a warning, stop execution with an error message, or forcefully crash the program.
    !    Output Units: You can set and reset the output unit where the assertions are printed with SetDefaultUnit and ReSetDefaultUnit subroutines. The default output unit is the standard error.
    !    Default Assertion Levels: You can change the default assertion level with SetDefaultLevelToWarning, SetDefaultLevelToSevere, SetDefaultLevelToFatal, and ResetDefaultLevel subroutines.
    !    Unit Status Information: The module provides functionality to check if an output unit is valid for formatted output and to print the status of a unit with Unit_is_Valid and Print_Unit_Status subroutines respectively.
    !    Error Messages: A simpler subroutine, ErrorMessage, is also provided for writing messages to the standard error.
    !    Crashing the Program: In situations where you want to crash the program to aid in debugging, you can use the crash subroutine. It intentionally provokes a segmentation fault which should produce a core dump that you can inspect with a debugger.

    use, intrinsic :: ISO_FORTRAN_ENV
    !use, intrinsic :: IEEE_ARITHMETIC
    !use, intrinsic :: IEEE_EXCEPTIONS
    !use, intrinsic :: IEEE_FEATURES

    implicit none

    !> The name of the module.
    character(len = *), parameter :: MODULE_NAME = "ModuleErrorHandling"
    !> The maximum length of the message to be printed.
    integer, public, parameter :: IOMSG_LENGTH = 512

    private

    !> Sets the level of the asserion level.
    type ClassAssertionLevel
        !> Assertion level.
        integer :: level
        !> Assertion name.
        character(len = 16) :: name
    end type ClassAssertionLevel
    !
    !.. Private parameters
    !..
    type(ClassAssertionLevel), parameter :: LEVEL_WARNING = ClassAssertionLevel (1, "WARNING")
    type(ClassAssertionLevel), parameter :: LEVEL_SEVERE = ClassAssertionLevel (2, "SEVERE")
    type(ClassAssertionLevel), parameter :: LEVEL_FATAL = ClassAssertionLevel (3, "FATAL")
    type(ClassAssertionLevel), parameter :: LEVEL_DEFAULT = LEVEL_FATAL
    !
    !
    !> SINGLETON Collection of enumerated types for different error levels.
    type class_assertion
        !
        type(ClassAssertionLevel) :: LEVEL_WARNING = LEVEL_WARNING
        type(ClassAssertionLevel) :: LEVEL_SEVERE = LEVEL_SEVERE
        type(ClassAssertionLevel) :: LEVEL_FATAL = LEVEL_FATAL
        type(ClassAssertionLevel) :: LEVEL_DEFAULT = LEVEL_DEFAULT
        integer :: UNIT_DEFAULT = ERROR_UNIT
        !
    contains
        !
        !> Sets the default level to Warning.
        procedure, nopass :: SetDefaultLevelToWarning
        !> Sets the default level to Severe.
        procedure, nopass :: SetDefaultLevelToSevere
        !> Sets the default level to Fatal.
        procedure, nopass :: SetDefaultLevelToFatal
        !> Resets the level to the default value.
        procedure, nopass :: ResetDefaultLevel
        !> Sets the default output assertion unit to unit.
        !! If unit is open to formatted output, return true, otherwise print an error message to standard error and return false.
        procedure, nopass :: SetDefaultUnit
        !> Resets the default unit.
        procedure, nopass :: ReSetDefaultUnit
        !
    end type class_assertion
    !
    type(class_assertion), protected :: assertion
    public :: Assertion

    interface operator (==)
        module procedure CompareAssertLevels
    end interface

    interface ConvertToString
        module procedure int2strn
        module procedure bool2strn
    end interface ConvertToString


    !.. Public routines
    public :: Assert
    public :: ErrorMessage
    public :: StopExecution


contains

    !> Stop the program
    !!
    subroutine StopExecution()
        STOP
    end subroutine StopExecution


    !> Issue an assertion and alter execution, if required.
    !!
    !! Usage:
    !! Example 1) issue an error message to standar error and crash
    !!            call assert( "myfun:init:invalid argument" )
    !! Example 2) issue a warning to standard error but does not stop execution
    !!            call assert( "myfun:iter:low accuracy", assertion%LEVEL_WARNING )
    !! Example 3) issue a warning to a user defined output unit
    !!            call assert( "myfun:eval:factorial too large", assertion%LEVEL_SEVERE, fpt_err )
    !<
    subroutine assert(Message, Level, Unit)
        !
        use, intrinsic :: ISO_FORTRAN_ENV
        !
        implicit none
        !
        !.. Dummy arguments
        !..
        !> Message to be printed
        character(len = *), intent(in) :: Message
        !> level of assertion, e.g. assertion%LEVEL_WARNING
        type(ClassAssertionLevel), optional, intent(in) :: Level
        !> output unit. Default = ISO_FORTRAN_ERR:ERROR_UNIT
        integer, optional, intent(in) :: Unit
        !
        !.. Local variables
        !..
        character(len = *), parameter :: LINE_HEADER = "*** "
        !
        integer :: current_output_unit
        type(ClassAssertionLevel) :: current_level


        !.. Set current unit to default, if valid,
        !   and to standard error otherwise.
        !..
        if(.not. Unit_is_Valid(assertion%UNIT_DEFAULT))then
            !
            write(ERROR_UNIT, "(a)") "Assert Internal Error: default unit became invalid"
            call Print_Unit_Status(assertion%UNIT_DEFAULT)
            write(ERROR_UNIT, "(a)") "Reset default unit to ERROR_UNIT"
            assertion%UNIT_DEFAULT = ERROR_UNIT
            !
        endif
        current_output_unit = assertion%UNIT_DEFAULT
        !
        !.. Set a non-default Output Unit, if required and valid
        !..
        if(present(Unit))then
            if(Unit_is_Valid(Unit))then
                current_output_unit = unit
            else
                write(current_output_unit, "(a)") "Required Assertion Unit is invalid"
                call Print_Unit_Status(Unit)
            endif
        endif


        !.. Set a non-default assertion level, if required
        !..
        current_level = assertion%LEVEL_DEFAULT
        if(present(Level)) Current_Level = Level


        !.. Write error message
        !..
        write(CURRENT_OUTPUT_UNIT, "(a)") &
                LINE_HEADER // &
                        trim(Current_Level%Name) // ": " // &
                        trim(Message)


        !.. Force crash if the error is fatal
        !..
        if(Current_level == assertion%LEVEL_FATAL)then
            !
            write(CURRENT_OUTPUT_UNIT, "(a)") LINE_HEADER // "CRASH FORCED"
            call crash()
            !
        endif

        return
        !
    end subroutine assert


    !> Sets the default output assertion unit to unit.
    !! If unit is open to formatted output, return true, otherwise print an error message to standard error and return false.
    !
    subroutine SetDefaultUnit(unit, IOSTAT)
        integer, intent(in) :: unit
        integer, optional, intent(out) :: iostat
        if(Unit_is_Valid(Unit))then
            assertion%UNIT_DEFAULT = unit
            if(present(iostat))iostat = 0
        else
            call ASSERT("Invalid unit")
            call Print_Unit_Status(unit)
            if(present(iostat))iostat = -1
        endif
        return
    end subroutine SetDefaultUnit
    !
    !> Resets the default unit.
    subroutine ReSetDefaultUnit()
        assertion%UNIT_DEFAULT = ERROR_UNIT
    end subroutine ReSetDefaultUnit


    !.. Default-Level Initialization utils
    !..
    !> Sets the default level to Warning.
    subroutine SetDefaultLevelToWarning
        assertion%LEVEL_DEFAULT = LEVEL_WARNING
    end subroutine SetDefaultLevelToWarning
    !
    !> Sets the default level to Severe.
    subroutine SetDefaultLevelToSevere
        assertion%LEVEL_DEFAULT = LEVEL_SEVERE
    end subroutine SetDefaultLevelToSevere
    !
    !> Sets the default level to Severe.
    subroutine SetDefaultLevelToFatal
        assertion%LEVEL_DEFAULT = LEVEL_FATAL
    end subroutine SetDefaultLevelToFatal
    !
    !> Resets the level to the default value.
    subroutine ResetDefaultLevel
        assertion%LEVEL_DEFAULT = LEVEL_DEFAULT
    end subroutine ResetDefaultLevel


    !> Ascertain if Unit is open to formatted writing, in
    !!  which case it returns true, otherwise issue a warning
    !!  to standard error and return false
    !<
    logical function Unit_is_Valid(Unit) result(Valid)
        !
        implicit none
        !
        integer, intent(in) :: Unit
        !
        character(len = 16) :: Unit_is_Formatted = " "
        character(len = 16) :: Unit_is_Writable = " "
        character(len = 16) :: Unit_IO_Message = " "
        logical :: Unit_is_Open = .TRUE.
        integer :: Unit_IO_Status = 0
        !
        INQUIRE(&
                UNIT = Unit, &
                IOSTAT = Unit_IO_Status, &
                OPENED = Unit_is_Open, &
                FORMATTED = Unit_is_Formatted, &
                WRITE = Unit_is_Writable, &
                IOMSG = Unit_IO_Message)
        !
        Valid = Unit_is_Open .and. &
                trim(Unit_is_Formatted) == "YES" !.and. &
        !trim( Unit_is_Writable )  == "YES"
        !
        return
        !
    end function Unit_is_Valid


    !> Print the status of a given unit.
    !
    !.. Todo : add optional output unit argument
    !..
    subroutine Print_Unit_Status(Unit)
        !
        integer, intent(in) :: Unit
        !
        character(len = 256) :: Unit_File_Name = " "
        logical :: Unit_is_Open = .TRUE.
        character(len = 16) :: Unit_Position = " "
        character(len = 16) :: Unit_Action = " "
        character(len = 16) :: Unit_Form = " "
        character(len = 16) :: Unit_Access = " "
        integer :: Unit_IO_Status = 0
        character(len = 256) :: Unit_IO_Message = " "
        !
        INQUIRE(&
                UNIT = Unit, &
                NAME = Unit_File_Name, &
                OPENED = Unit_is_Open, &
                POSITION = Unit_Position, &
                ACTION = Unit_Action, &
                FORM = Unit_Form, &
                ACCESS = Unit_Access, &
                IOSTAT = Unit_IO_Status, &
                IOMSG = Unit_IO_Message)
        !
        write(ERROR_UNIT, "(a)") "----------------------------------"
        write(ERROR_UNIT, "(a)") "UNIT STATUS REPORT:"
        write(ERROR_UNIT, "(a)") "UNIT......." // ConvertToString(unit)
        write(ERROR_UNIT, "(a)") "NAME.......'" // trim(Unit_File_Name) // "'"
        write(ERROR_UNIT, "(a)") "OPENED....." // ConvertToString(Unit_is_Open)
        write(ERROR_UNIT, "(a)") "POSITION..." // Unit_Position
        write(ERROR_UNIT, "(a)") "ACTION....." // Unit_Action
        write(ERROR_UNIT, "(a)") "FORM......." // Unit_Form
        write(ERROR_UNIT, "(a)") "ACCESS....." // Unit_Access
        write(ERROR_UNIT, "(a)") "IOSTAT....." // ConvertToString(Unit_IO_Status)
        write(ERROR_UNIT, "(a)") "IOMSG......" // trim(Unit_IO_Message)
        write(ERROR_UNIT, "(a)") "----------------------------------"
        !
        return
        !
    end subroutine Print_Unit_Status


    !> Convert an integer to a left-aligned string
    !! of minimal length
    !<
    character(len = 16) function int2strn(int_num) result(strn)
        integer, intent(in) :: int_num
        write(strn, *)int_num
        strn = trim(adjustl(strn))
    end function int2strn
    !
    !> Convert a boolean to a left-aligned string
    !! of minimal length
    !<
    character(len = 16) function bool2strn(bool) result(strn)
        logical, intent(in) :: bool
        write(strn, "(L)")bool
        strn = trim(adjustl(strn))
    end function bool2strn


    !> Determines if two assertion levels coincide
    !<
    logical function CompareAssertLevels(level1, level2) result(equal)
        implicit none
        type(ClassAssertionLevel), intent(in) :: level1, level2
        equal = .FALSE.
        if(level1%level /= level2%level)return
        if(trim(level1%name) /= trim(level2%name))return
        equal = .TRUE.
        return
    end function CompareAssertLevels


    !> Write a message to standard error
    subroutine ErrorMessage(Message)
        implicit none
        character(len = *), intent(in) :: Message
        call Assert(Message, assertion%LEVEL_WARNING)
    end subroutine ErrorMessage


    !> Cause a crash of the program by  provoking a segmentation fault.
    !! A call to crash helps trace where it was invoked.
    !<
    subroutine crash()
        integer, allocatable :: A(:)
        deallocate(A)
        return
    end subroutine crash


end module ModuleErrorHandling
