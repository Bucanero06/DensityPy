! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
Module ModuleMainInterface
  
  private

  public :: GetRunTimeParameters
  
contains
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters( ProgramInputFile )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    !.. Level 0.1
    use ModuleCommandLineParameterList
    !.. Level 2
    use ModuleAstraCredit
    !
    implicit none
    !
    character(len=:), allocatable, intent(out)   :: ProgramInputFile

    type( ClassCommandLineParameterList ) :: List
    character(len=100)            :: StrnBuf

    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Converts the integrals from UKRmol+ to the Astra format "

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help"  , "Print Command Usage" )
    call List%Add( "-gif"    , "Astra Main Config File", "ASTRA.INP", "optional" )
    call List%Add( "-dummy"  , "dummy" )

    call List%Parse()

    call PrintAstraName()
    call PrintAstraCredit()
    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-gif"  , StrnBuf ); ProgramInputFile = trim( StrnBuf )
    
    call List%Free()

    write(OUTPUT_UNIT,"(a)" )
    write(OUTPUT_UNIT,"(a)" )   "Read run time parameters :"
    write(OUTPUT_UNIT,"(a)" )
    write(OUTPUT_UNIT,"(a)" )   "Program Input File  : "//ProgramInputFile

  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
