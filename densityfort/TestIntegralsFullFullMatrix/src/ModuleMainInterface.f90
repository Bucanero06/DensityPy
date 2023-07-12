Module ModuleMainInterface
  
  private

  public :: GetRunTimeParameters
  public :: ParseProgramInputFile
  

contains
  
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters( ProgramInputFile ,AssignedSymmetry )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleCommandLineParameterList
    !
    implicit none
    !
    character(len=:), allocatable, intent(out)   :: ProgramInputFile
    character(len=:), allocatable, intent(out)   :: AssignedSymmetry

    type( ClassCommandLineParameterList ) :: List
    character(len=100)            :: StrnBuf

    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Converts the integrals from UKRmol+ to the Astra format "

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help"  , "Print Command Usage" )
    call List.Add( "-gif"    , "General Input File"   , "TESTINTEGRALS.INP", "optional" )
    call List.Add( "-dummy"  , "dummy" )
    call List.Add( "-sym"    , "Assigned Symmetry", "Ag", "optional")

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if

    call List.Get( "-gif"  , StrnBuf ); ProgramInputFile = trim( StrnBuf )
    call List.Get( "-sym"  , StrnBuf ); AssignedSymmetry = trim( StrnBuf )
    
    call List.Free()

    write(OUTPUT_UNIT,"(a)" )
    write(OUTPUT_UNIT,"(a)" )   "Read run time parameters :"
    write(OUTPUT_UNIT,"(a)" )
    write(OUTPUT_UNIT,"(a)" )   "Program Input File  : "//ProgramInputFile
    write(OUTPUT_UNIT,"(a)" )   "Assigned Symmetry   : "//AssignedSymmetry

  end subroutine GetRunTimeParameters


  !> Parses the program input file.
  subroutine ParseProgramInputFile( &
       ProgramInputFile, &
       StorageDir      , &
       QCDir           , &
       ccConfigFile    )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleParameterList
    !
    implicit none
    !
    character(len=*)             , intent(in)  :: ProgramInputFile
    character(len=:), allocatable, intent(out) :: StorageDir
    character(len=:), allocatable, intent(out) :: QCDir
    character(len=:), allocatable, intent(out) :: ccConfigFile

    character(len=100)       :: strnBuf
    type(ClassParameterList) :: List

    call List.Add( "StorageDir"      , "store/", "optional" )
    call List.Add( "QCDir"           , "Lucia/", "optional" )
    call List.Add( "ccConfigFile"    , "CLSCPLNG.INP", "optional" )

    call List.Parse( ProgramInputFile )

    call List.Get( "StorageDir", strnBuf )
    allocate( StorageDir, source = trim(adjustl(strnBuf)) )

    call List.Get( "QCDir", strnBuf )
    allocate( QCDir     , source = trim(adjustl(strnBuf)) )

    call List.Get( "ccConfigFile", strnBuf )
    allocate( ccConfigFile, source = trim(adjustl(strnBuf)) )


    write(OUTPUT_UNIT,"(a)") "Storage Dir.................."//StorageDir
    write(OUTPUT_UNIT,"(a)") "QC Dir......................."//QCDir
    write(OUTPUT_UNIT,"(a)") "CC Config File..............."//ccConfigFile

  end subroutine ParseProgramInputFile


end Module ModuleMainInterface
