Module ModuleMainInterface
  
  private

  public :: GetRunTimeParameters
  public :: ParseProgramInputFile
  

contains
  
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters( ProgramInputFile ,AssignedSymmetry,lSWITCH )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleCommandLineParameterList
    !
    implicit none
    !
    character(len=:), allocatable, intent(out)   :: ProgramInputFile
    character(len=:), allocatable, intent(out)   :: AssignedSymmetry
    logical                      , intent(out)   :: lSWITCH(3,3,3)
    

    type( ClassCommandLineParameterList ) :: List
    character(len=100)            :: StrnBuf
    character(len=9)              :: StrnVec(3)
    integer                       :: i, j, k, iOp

    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Run the test between the atra STEX matrix elements blocks and the ones from the benchmark program "

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help"  , "Print Command Usage" )
    call List.Add( "-gif"    , "General Input File"   , "TESTINTEGRALS.INP", "optional" )
    call List.Add( "-sym"    , "Assigned Symmetry", "Ag", "optional")
    call List%Add( "-switchOver" , "overlap cases to compute: AA AI AE IA ...", "000000000", "optional" )
    call List%Add( "-switchMono" , "overlap cases to compute: AA AI AE IA ...", "000000000", "optional" )
    call List%Add( "-switchBiel" , "overlap cases to compute: AA AI AE IA ...", "000000000", "optional" )
    
    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if

  !.. Test Overlap (one by one, and then add the rest progressively
  !   AA AI AE IA II IE EA EI EE
  !    *
  !    * * 
  !    .... 
  !.. Test Monoelectronic (one by one, and then add the rest progressively
  !   AA AI AE IA II IE EA EI EE
  !    *
  !    * * 
  !    .... 
  !.. Test Bielectronic (one by one, and then add the rest progressively
  !   AA AI AE IA II IE EA EI EE
  !    *
  !    * * 
  !    ....
    call List%Get( "-switchOver",StrnVec(1) ); 
    call List%Get( "-switchMono",StrnVec(2) ); 
    call List%Get( "-switchBiel",StrnVec(3) ); 
    k = 0
    do iOp = 1, 3
       do j = 1, 3
          do i = 1, 3
             k = k + 1
             lSwitch(i,j,iOp) = StrnVec(iOp)(k:k) == "1"
          enddo
       enddo
       k = 0
    enddo

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
