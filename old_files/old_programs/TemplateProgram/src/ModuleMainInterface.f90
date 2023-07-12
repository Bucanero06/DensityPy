!! CONFIDENTIAL
Module ModuleMainInterface

  private
  
  public :: GetRunTimeParameters

contains

  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( FileName, nSize )
    !
    use, intrinsic :: ISO_FORTRAN_ENV

    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: FileName
    integer                      , intent(out) :: nSize
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Template Programs which diagonalizes a random matrix"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help" , "Print Command Usage" )
    call List%Add( "-o"     , "Output File" ,"eval", "optional" )
    call List%Add( "-n"     , "Matrix dim"  , 100  , "required" )

    call List%Parse()

    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-o",  strnBuf  )
    allocate(FileName,source=trim(adjustl(strnBuf)))

    call List%Get( "-n", nSize )

    call List%Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
    write(OUTPUT_UNIT,"(a)"   ) "Output File : "//FileName
    write(OUTPUT_UNIT,"(a,i0)") "Matrix size : ",nSize
    !
  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
