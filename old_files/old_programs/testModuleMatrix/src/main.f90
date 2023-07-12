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

  use ModuleMainInterface

  implicit none

  !.. Run-time parameters
  !..
  integer                       :: nSize
  character(len=:), allocatable :: FileName

  type(ClassMatrix) :: bMat
  real(kind(1d0)), allocatable :: dMat1(:,:), dMat2(:,:)
  
  call GetRunTimeParameters( FileName, nSize )

  call bMat%Init(...)
  do i
     do j
     enddo
  enddo
  dMat1 = bMat
  call bMat%Write(FileName)
  call bMat%Free()
  call bMat%Read(FileName)
  
  dMat2=bMat
  write(*,*) sum(abs(dMat1-dMat2))

  !.. Test which elements are non-zero and if the dimensions match
  !***

  !.. Test access 
  do i = 1, n
     do j = max(1,i-k+1), min(n,i+k-1)
        dBuf = bMat%Element(i,j)
        if(abs(dBuf-dMat1(i,j))>THRESHOLD)then
           write(*,*) ...
        endif
     enddo
  enddo

  stop
  
end program ProgramTemplate

