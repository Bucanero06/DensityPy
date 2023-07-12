Module ModuleMainInterface
  
  private

  public :: GetRunTimeParameters
  

contains
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters( ProgramInputFile, CheckConversion, STEXapproximation )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleCommandLineParameterList
    use ModuleAstraCredit
    !
    implicit none
    !
    character(len=:), allocatable, intent(out)   :: ProgramInputFile
    logical                      , intent(out)   :: CheckConversion
    logical                      , intent(out)   :: STEXApproximation
    
    type( ClassCommandLineParameterList ) :: List
    character(len=100)            :: StrnBuf

    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Converts the integrals from UKRmol+ to the Astra format "

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help"  , "Print Command Usage" )
    call List%Add( "-gif"    , "Astra Main Config File", "ASTRA.INP", "optional" )
    call List%Add( "-check"  , "Check Conversion"    )
    call List%Add( "-stex"   , "Use Static-Exchange Approximation")
    call List%Add( "-dummy"  , "dummy" )

    call List%Parse()

    call PrintAstraName()
    call PrintAstraCredit()
    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-gif"  , StrnBuf ); ProgramInputFile = trim( StrnBuf )
    CheckConversion   = List%Present("-check")
    STEXapproximation = List%Present("-stex")
    
    call List%Free()

    write(OUTPUT_UNIT,"(a)" )
    write(OUTPUT_UNIT,"(a)" )   "Read run time parameters :"
    write(OUTPUT_UNIT,"(a)" )
    write(OUTPUT_UNIT,"(a)" )   "Program Input File  : "//ProgramInputFile

  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
