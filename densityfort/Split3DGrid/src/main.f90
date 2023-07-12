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
  use ModuleBSpline
  use ModuleIO
  use ModuleConstants

  implicit none

  !.. Run-time parameters
  !..
  character(len=:), allocatable :: InpFile
  character(len=:), allocatable :: OutFile
  integer                       :: nval, iCoord
  real(kind(1d0))               :: CoordVal, CoordSpan

  integer                       :: i, iv, iSlow
  real(kind(1d0))               :: xv(3), PrevSlowCoord

  integer                       :: uidInp, uidOut, iostat
  character(len=100)            :: iomsg

  !!!!!!!!!!!!!!11
  integer                       :: nPoints,j,jumpvariable,NewNPoints, k
  real(kind(1d0)) , allocatable :: coordinates(:,:), dvec(:) , avecs(:)
  real(kind(1d0)), allocatable  :: editedcoordinates(:,:),editedd(:)
  character(len=10000)          :: line
  real(kind(1d0))               :: leftover


  call GetRunTimeParameters( InpFile, OutFile, nval, iCoord, CoordVal, CoordSpan )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!>.. Reads number of points
  !.. Open Input File
  open(newunit= uidInp    , &
       file   = InpFile   , &
       form   ="formatted", &
       status ="old"      , &
       action ="read"     , &
       iostat = iostat    , &
       iomsg  = iomsg     )
  if(iostat /= 0)then
     call ErrorMessage("Unable to open "//InpFile//" iomsg: "//trim(iomsg))
     stop
  endif

  !.. Number of filled lines = nPoints of grid
  nPoints=0
  do
     read(uidInp,"(a)",iostat=iostat) line
     if(iostat/= 0) exit
     if(len_trim(line)==0)cycle
     nPoints = nPoints+1
  enddo

  !.. Close input file
  close(uidInp)


!>.. Start
  !.. Open Input File
  open(newunit= uidInp    , &
       file   = InpFile   , &
       form   ="formatted", &
       status ="old"      , &
       action ="read"     , &
       iostat = iostat    , &
       iomsg  = iomsg     )
  if(iostat /= 0)then
     call ErrorMessage("Unable to open "//InpFile//" iomsg: "//trim(iomsg))
     stop
  endif

  !..Open Output File
  open(newunit= uidOut    , &
       file   = OutFile   , &
       form   ="formatted", &
       status ="unknown"  , &
       action ="write"    , &
       iostat = iostat    , &
       iomsg  = iomsg     )
  if(iostat /= 0)then
     call ErrorMessage("Unable to open "//OutFile//" iomsg: "//trim(iomsg))
     stop
  endif


  jumpvariable= 25
  NewNPoints=nPoints/jumpvariable
  allocate(dvec(nPoints))
  allocate(editedd(NewNPoints))
  allocate(coordinates(3,nPoints))
  allocate(editedcoordinates(3,NewNPoints))


!>.. Sums over z space
  !..  Read Values of columns
  do i=1,nPoints
    read(uidInp,*,iostat=iostat) coordinates(1,i),coordinates(2,i),coordinates(3,i),dvec(i)
    if(iostat/=0)exit
  end do

  !..  Summ over every jumpvariable th's element
  j=1
  do i=1, NewNPoints
    editedd(i)= sum(dvec(j: j+(jumpvariable-1)))
    j=(i*jumpvariable)+1
    if (j >=nPoints) exit
  enddo

  !..  Write to outfile
  j=1
  write(uidOut,*)
  do i =1,NewNPoints
    write(uidOut,"(*(x,e24.16))") coordinates(1,j),&
                                  coordinates(2,j),&
                                  coordinates(3,j),&
                                  editedd(i)
    j=(i*jumpvariable) + 1
    if (j >= nPoints) then
      exit
    else
      if (coordinates(1,j) .ne. coordinates(1,(j-24))) then
        write(uidOut,*)
      end if
    end if
  end do


  close( uidInp )
  close( uidOut )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  allocate(avecs(nval))
!  !3D plot
!  do i=1,nPoints
!    read(uidInp,*,iostat=iostat) coordinates(1,i),coordinates(2,i),coordinates(3,i), dvec(i)
!  end do
!
!  do i=1,nPoints
!    if (coordinates(1,i) .ne. coordinates(1,i-1)) then
!      write(uidOut,*)
!    end if
!    write(uidOut,"(*(x,e24.16))") coordinates(1,i),coordinates(2,i),coordinates(3,i), dvec(i)
!
!  end do
!
!  close( uidInp )
!  close( uidOut )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1





!!!!!!!!!!!!!!!!!!!!!!!!!Done With Argenti (Original)
!Cuts a slice of grid
!
!  allocate(val(nval))
!  iSlow = 1
!  if(iCoord==1)iSlow=2
!
!  PrevSlowCoord=-1.d50
!  do
!
!     read(uidInp,*,iostat=iostat) (xv(i),i=1,3), (val(iv),iv=1,nval)
!     if(iostat/=0)exit
!     if(abs(xv(iCoord)-CoordVal)>CoordSpan)cycle
!
!     if(xv(iSlow)/=PrevSlowCoord)then
!        write(uidOut,*)
!        PrevSlowCoord = xv(iSlow)
!     endif
!     write(uidOut,"(*(x,e24.16))") (xv(i),i=1,3), (val(iv),iv=1,nval)
!
!  enddo
!
!  close( uidInp )
!  close( uidOut )
!!!!!!!!!!!!!!!!!!!!!!!!!!
  stop

contains


  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( InpFile, OutFile, nval, iCoord, CoordVal, CoordSpan )
    !
    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: InpFile
    character(len=:), allocatable, intent(out) :: OutFile
    integer                      , intent(out) :: nval, iCoord
    real(kind(1d0))              , intent(out) :: CoordVal, CoordSpan
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Split a 3D grid along specific xy, xz, or yz planes"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help" , "Print Command Usage" )
    call List%Add( "-i"     , "Input  File" ,"eval", "required" )
    call List%Add( "-o"     , "Output File" ,"eval", "required" )
    call List%Add( "-n"     , "n values"    , 1    , "optional" )
    call List%Add( "-c"     , "coord index (1:x, 2:y, 3:z)" ,3, "optional" )
    call List%Add( "-cval"  , "coord value" ,0.d0, "optional" )
    call List%Add( "-cspan" , "coord span" , 1.d-4, "optional" )
    call List%Add( "-xxx"   , "placeholder" )

    call List%Parse()

    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-i",  strnBuf  )
    allocate(InpFile,source=trim(adjustl(strnBuf)))
    call List%Get( "-o",  strnBuf  )
    allocate(OutFile,source=trim(adjustl(strnBuf)))
    call List%Get( "-n",  nval  )
    call List%Get( "-c",  iCoord  )
    call List%Get( "-cval",  CoordVal  )
    call List%Get( "-cspan",  CoordSpan  )

    call List%Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
    write(OUTPUT_UNIT,"(a)"   ) "Input  File : "//InpFile
    write(OUTPUT_UNIT,"(a)"   ) "Output File : "//OutFile
    !
  end subroutine GetRunTimeParameters

end program ProgramTemplate

