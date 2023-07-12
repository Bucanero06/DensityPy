!! CONFIDENTIAL
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

  use ModuleMainInterface

  use ModuleBSpline
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

  call GetRunTimeParameters( FileName, nSize )

  allocate(dMat(nSize,nSize))
  call random_number(dMat)
  dMat = dMat + transpose(dMat)

  allocate(dEval(nSize))
  call Short_Diag( nSize, dMat, dEval )

  call SaveVector( FileName, dEval, "formatted" )

end program ProgramTemplate

