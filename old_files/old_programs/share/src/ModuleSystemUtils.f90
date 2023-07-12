!! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
!> Declares some utility routines based on system calls.
module ModuleSystemUtils

  implicit none

  private

  public :: GetTime
  public :: ElapsedTime

  interface 
     Pure DoublePrecision function D2DFun( x, parvec ) 
       DoublePrecision          , intent(in) :: x
       DoublePrecision, optional, intent(in) :: parvec(*)
     end function D2DFun
  end interface
  public D2DFun

  interface realloc
     module procedure realloc_integer_vec
     module procedure realloc_dble_vec
     module procedure realloc_complex_vec
     module procedure realloc_integer_mat2
     module procedure realloc_dble_mat2
     module procedure realloc_complex_mat2
     module procedure realloc_dble_mat4
  end interface realloc
  public realloc

contains

  subroutine realloc_integer_vec(vec,n)
    integer, allocatable, intent(inout) :: vec(:)
    integer             , intent(in) :: n
    if(allocated(vec))then
       if(size(vec,1)/=n)then
          deallocate(vec)
          allocate(vec(n))
       endif
    else
       allocate(vec(n))
    endif
    vec=0
  end subroutine realloc_integer_vec

  subroutine realloc_dble_vec(vec,n)
    real(kind(1d0)), allocatable, intent(inout) :: vec(:)
    integer             , intent(in) :: n
    if(allocated(vec))then
       if(size(vec,1)/=n)then
          deallocate(vec)
          allocate(vec(n))
       endif
    else
       allocate(vec(n))
    endif
    vec=0.d0
  end subroutine realloc_dble_vec
  
  subroutine realloc_complex_vec(vec,n)
    complex(kind(1d0)), allocatable, intent(inout) :: vec(:)
    integer             , intent(in) :: n
    if(allocated(vec))then
       if(size(vec,1)/=n)then
          deallocate(vec)
          allocate(vec(n))
       endif
    else
       allocate(vec(n))
    endif
    vec=(0.d0,0.d0)
  end subroutine realloc_complex_vec

  subroutine realloc_integer_mat2(mat,n1,n2)
    integer, allocatable, intent(inout) :: mat(:,:)
    integer             , intent(in) :: n1,n2
    if(allocated(mat))then
       if(size(mat,1)/=n1.or.size(mat,2)/=n2)then
          deallocate(mat)
          allocate(mat(n1,n2))
       endif
    else
       allocate(mat(n1,n2))
    endif
    mat=0
  end subroutine realloc_integer_mat2

  subroutine realloc_dble_mat2(mat,n1,n2,EXPAND)
    real(kind(1d0)), allocatable, intent(inout) :: mat(:,:)
    integer             , intent(in) :: n1,n2
    logical, optional   , intent(in) :: EXPAND
    logical :: EXPAND_, SET_EXACT_SIZE
    EXPAND_=.FALSE.
    if(present(EXPAND))EXPAND_=EXPAND
    SET_EXACT_SIZE=.not.EXPAND_
    if(allocated(mat))then
       if(EXPAND_ .AND. ( size(mat,1)<n1 .OR. size(mat,2)<n2 ) )then
          deallocate(mat)
          allocate(mat(n1,n2))
       elseif(SET_EXACT_SIZE .AND. ( size(mat,1)/=n1 .OR. size(mat,2)/=n2 ) )then
          deallocate(mat)
          allocate(mat(n1,n2))
       endif
    else
       allocate(mat(n1,n2))
    endif
    mat=0.d0
  end subroutine realloc_dble_mat2

  subroutine realloc_dble_mat4(mat,n1,n2,n3,n4,EXPAND)
    real(kind(1d0)), allocatable, intent(inout) :: mat(:,:,:,:)
    integer             , intent(in) :: n1,n2,n3,n4
    logical, optional   , intent(in) :: EXPAND
    logical :: EXPAND_, SET_EXACT_SIZE
    EXPAND_=.FALSE.
    if(present(EXPAND))EXPAND_=EXPAND
    SET_EXACT_SIZE=.not.EXPAND_
    if(allocated(mat))then
       if(EXPAND_ .AND. ( size(mat,1)<n1 .OR. size(mat,2)<n2 .OR. size(mat,3)<n3 .OR. size(mat,4)<n4 ) )then
          deallocate(mat)
          allocate(mat(n1,n2,n3,n4))
       elseif(SET_EXACT_SIZE .AND. ( size(mat,1)/=n1 .OR. size(mat,2)/=n2 .OR. size(mat,3)/=n3 .OR. size(mat,4)/=n4 ) )then
          deallocate(mat)
          allocate(mat(n1,n2,n3,n4))
       endif
    else
       allocate(mat(n1,n2,n3,n4))
    endif
    mat=0.d0
  end subroutine realloc_dble_mat4

  subroutine realloc_complex_mat2(mat,n1,n2)
    complex(kind(1d0)), allocatable, intent(inout) :: mat(:,:)
    integer             , intent(in) :: n1,n2
    if(allocated(mat))then
       if(size(mat,1)/=n1.or.size(mat,2)/=n2)then
          deallocate(mat)
          allocate(mat(n1,n2))
       endif
    else
       allocate(mat(n1,n2))
    endif
    mat=(0.d0,0.d0)
  end subroutine realloc_complex_mat2

  !> Returns the absolute CPU clock time in 
  !> seconds with a precision of (allegedly)
  !> one microsecond.
  DoublePrecision function GetTime()
    integer(8) :: count, count_rate
    call system_clock(count,count_rate)
    GetTime=dble(count)/dble(count_rate)
  end function GetTime


  !> Returns the seconds passed since 
  !! the previous call to the function
  subroutine ElapsedTime(etime)
    real(kind(1d0)), optional, intent(out) :: etime
    logical        , save :: FIRST_CALL=.TRUE.
    real(kind(1d0)), save :: tic_time
    integer(8)     , save :: count0
    integer(8)            :: count, count_rate
    if(FIRST_CALL)then
       call system_clock(count,count_rate)
       FIRST_CALL=.FALSE.
       tic_time=1.d0/dble(count_rate)
       count0=count
    else
       call system_clock(count)
    endif
    if(present(etime)) etime = dble(max(count-count0,0)) * tic_time
    count0=count
  end subroutine ElapsedTime


end module ModuleSystemUtils
