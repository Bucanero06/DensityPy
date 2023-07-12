! {{{ Detailed description

!> \mainpage Program <ProgramName> <Insert here what the program does>
!! 
!! Synopsis:
!! ---------
!!
!!     <Program Name> <mandatory run-time parameters (RTP)> [<optional RTP>]
!!
!! ___
!! Description:
!! ------------
!!
!! Input parameters:      {#Input_Parameters}
!! =================
!! [...Input](@ref ...Input) as specified in the command line.
!!
!> \file
!!
!!
! }}}
program Programtestlibcint

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO
  use ModuleConstants

  implicit none

  !.. Run-time parameters
  !..
  integer                       :: nSize
  character(len=:), allocatable :: FileName
  
  !.. Local parameters
  real(kind(1d0)), allocatable :: dMat(:,:)
  real(kind(1d0)), allocatable :: dEval(:)
  integer :: i,j

  call GetRunTimeParameters( FileName, nSize )

  allocate(dMat(nSize,nSize))
  call random_number(dMat)
  dMat = dMat + transpose(dMat)

  allocate(dEval(nSize))
  dEval = dMat(:,1)

  call SaveVector( FileName, dEval, "formatted" )

  stop

contains


  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( FileName, nSize )
    !
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

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help" , "Print Command Usage" )
    call List.Add( "-o"     , "Output File" ,"eval", "optional" )
    call List.Add( "-n"     , "Matrix dim"  , 100  , "required" )

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if

    call List.Get( "-o",  strnBuf  )
    allocate(FileName,source=trim(adjustl(strnBuf)))

    call List.Get( "-n", nSize )

    call List.Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
    write(OUTPUT_UNIT,"(a)"   ) "Output File : "//FileName
    write(OUTPUT_UNIT,"(a,i0)") "Matrix size : ",nSize
    !
  end subroutine GetRunTimeParameters

end program Programtestlibcint

