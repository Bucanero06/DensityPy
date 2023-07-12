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
program ProgramTemplate

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO
  use ModuleConstants
  use ModuleDiagonalize

  implicit none

  !.. Run-time parameters
  !..
  integer                       :: nSize
  character(len=:), allocatable :: FileName
  
  !.. Local parameters
  real(kind(1d0)), allocatable :: dMat(:,:)
  real(kind(1d0)), allocatable :: dEval(:)
  integer :: i,j

  integer        , parameter :: BS_ORDER = 6
  real(kind(1d0)), parameter :: XMIN     = 0.d0
  real(kind(1d0)), parameter :: XMAX     = 1.d0
  
  integer        , parameter :: NNODES     = 101
  integer        , parameter :: NNODEB     =  51
  integer        , parameter :: NNODES_OUT = NNODES - NNODEB + 1 + BS_ORDER - 1
  
  real(kind(1d0)) :: Grid1(NNODES), GridOut(NNODES_OUT)
  real(kind(1d0)) :: dx

!!$  call GetRunTimeParameters( FileName, nSize )

  dx = (XMAX-XMIN)/dble(NNODES-1)
  
  !.. Define uniform Bspline grid
  !   Here we will eliminate the first and last Bspline
  Grid1( 1 ) = XMIN 
  do i = 2, NNODES
     Grid1( i ) = Grid1( i - 1 ) + dx 
     write(*,*) "grid 1", i, Grid1(i)
  enddo
  
  !.. Define Right Bspline Gridu
  !   For the second basis, we will use the Bsplines from the second to 
  !   NNODEB+BS_ORDER-2 Bsplines in the first grid
  !   and from the BS_ORDER-th to the one before last Bsplines in the second grid
  GridOut(1) = Grid1(NNODEB)
  do i = 2, NNODES_OUT
     if( i < 2 * BS_ORDER )then
        GridOut(i) = GridOut(i-1) + dx / 2.d0
     else
        GridOut(i) = GridOut(i-1) + dx
     endif
     write(*,*) "grid out", i, GridOut(i)
  enddo

!!$  call Short_Diag( nSize, dMat, dEval )
!!$
!!$  call SaveVector( FileName, dEval, "formatted" )

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
!!$    call List.Add( "-o"     , "Output File" ,"eval", "optional" )
!!$    call List.Add( "-n"     , "Matrix dim"  , 100  , "required" )

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if
    
    allocate(FileName,source="tmp")
    nSize=0
    
!!$    call List.Get( "-o",  strnBuf  )
!!$    allocate(FileName,source=trim(adjustl(strnBuf)))
!!$
!!$    call List.Get( "-n", nSize )

    call List.Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
!!$    write(OUTPUT_UNIT,"(a)"   ) "Output File : "//FileName
!!$    write(OUTPUT_UNIT,"(a,i0)") "Matrix size : ",nSize
    !
  end subroutine GetRunTimeParameters

end program ProgramTemplate

