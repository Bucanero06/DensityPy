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

  use, intrinsic :: ISO_FORTRAN_ENV

  private
  
  public :: GetRunTimeParameters

contains

    !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters( &
       AstraCfgFile, &
       TestCase    , &
       RunDalton   , &
       RunLuciaPIE , &
       RunLuciaHAIC, &
       RunLuciaTDMs, &
       RunScatci   , &
       RunConvTDMs , &
       RunConvInts )
    !
    use ModuleCommandLineParameterList
    use ModuleAstraConfigFile
    use ModuleAstraCredit
    use ModuleAstraEnv
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: AstraCfgFile
    character(len=:), allocatable, intent(out) :: TestCase
    logical                      , intent(out) :: RunDalton
    logical                      , intent(out) :: RunLuciaPIE
    logical                      , intent(out) :: RunLuciaHAIC
    logical                      , intent(out) :: RunLuciaTDMs
    logical                      , intent(out) :: RunScatci
    logical                      , intent(out) :: RunConvTDMs
    logical                      , intent(out) :: RunConvInts

    type( ClassCommandLineParameterList ) :: List
    character(len=100)            :: StrnBuf
    character(len=:), allocatable :: DefaultTestCase

    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Set up some input files, if requested, and carries out "//&
         "the calculations that are preliminary to running astra"

    DefaultTestCase = AstraEnv%GetDefaultTest()
    
    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help", "Print Command Usage" )
    call List%Add( "-gif"  , "Astra Config File", ASTRA_CFG_FILE_DEF, "optional" )
    call List%Add( "-test" , "Test Case"        , DefaultTestCase, "optional" )
    call List%Add( "-run"  , "which tasks run"  , "dehtiTI"   , "optional" )

    call List%Parse()

    call PrintAstraName()
    call PrintAstraCredit()
    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-gif" , StrnBuf ); AstraCfgFile = trim( StrnBuf )
    call List%Get( "-test", StrnBuf ); TestCase = trim( StrnBuf )
    call List%Get( "-run" , StrnBuf )
    RunDalton    = ( index(StrnBuf,"d") > 0 )
    RunLuciaPIE  = ( index(StrnBuf,"e") > 0 )
    RunLuciaHAIC = ( index(StrnBuf,"h") > 0 )
    RunLuciaTDMs = ( index(StrnBuf,"t") > 0 )
    RunScatci    = ( index(StrnBuf,"i") > 0 )
    RunConvTDMs  = ( index(StrnBuf,"T") > 0 )
    RunConvInts  = ( index(StrnBuf,"I") > 0 )
    call List%Free()

    write(OUTPUT_UNIT,"(a)"    ) " Astra Main Config File : "//AstraCfgFile
    write(OUTPUT_UNIT,"(a)"    ) " Default Test Case      : "//TestCase
    write(OUTPUT_UNIT,"(a,l1)" ) " Run Dalton             : ",RunDalton
    write(OUTPUT_UNIT,"(a,l1)" ) " Run Lucia PI Energies  : ",RunLuciaPIE
    write(OUTPUT_UNIT,"(a,l1)" ) " Run Lucia H Act.In.Ch. : ",RunLuciaHAIC
    write(OUTPUT_UNIT,"(a,l1)" ) " Run Lucia TDM 1B & 2B  : ",RunLuciaTDMs
    write(OUTPUT_UNIT,"(a,l1)" ) " Run ScarCI Integrals   : ",RunScatci
    write(OUTPUT_UNIT,"(a,l1)" ) " Run Convert TDMs       : ",RunConvTDMs
    write(OUTPUT_UNIT,"(a,l1)" ) " Run Convert Integrals  : ",RunConvInts

  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
