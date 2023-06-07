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
Module ModulePOSIX
  !   From www.tutorialspoint.com/inter_process_communication
  !.. Portable Operative System Interface (POSIx) defines the
  !   Application Programming Interfaces (API) needed to perform
  !   various system calls
  !
  !   FIFO − Communication between two unrelated processes. FIFO is a full duplex,
  !   meaning the first process can communicate with the second process and vice
  !   versa at the same time.
  !
  !.. Define some standard-posix system calls via c
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  private

  type, private :: ClassPosix
     private
   contains
     procedure, public :: chdir  => ClassPosix_chdir
     procedure, public :: getenv => ClassPosix_getenv
     procedure, public :: getpid => ClassPosix_getpid
  end type ClassPosix
  !.. Singleton
  type(ClassPosix), public :: Posix

  interface
     function getpid() bind(C, name='getpid')
       use, intrinsic :: iso_c_binding
       integer :: getpid
     end function getpid
     integer function chdir(path) bind(C,name="chdir")
       use, intrinsic :: iso_c_binding
       character(kind=c_char) :: path(*)
     end function chdir
     function getenv(var_name) result(path) bind(C,name="getenv")
       use, intrinsic :: iso_c_binding
       import
       character(kind=c_char) :: var_name(*)
       type(c_ptr) :: path
     end function getenv
     function strlen(str) result(isize) bind(C,name="strlen")
       use, intrinsic :: iso_c_binding
       import
       type   (c_ptr), value :: str
       integer(c_int)        :: isize
     end function strlen
  end interface
  
contains

  function ClassPosix_getpid( self ) result(pid)
    class(ClassPosix), intent(in) :: self
    integer :: pid
    pid = getpid()
  end function ClassPosix_getpid

  subroutine ClassPosix_chdir( self, path, istat )
    class(ClassPosix), intent(in) :: self
    character(len=*) , intent(in) :: path
    integer, optional, intent(out):: istat
    integer :: ierr
    ierr = chdir(path//c_null_char)
    if(present(istat))istat=ierr
  end subroutine ClassPosix_chdir

  subroutine ClassPosix_getenv( self, var_name, svar )
    class(ClassPosix)            , intent(in)  :: self
    character(len=*)             , intent(in)  :: var_name
    character(len=:), allocatable, intent(out) :: svar
    type   (c_ptr) :: cstr
    integer(c_int) :: n
    cstr=getenv(var_name//c_null_char)
    if ( c_associated(cstr) )then
       n =strlen(cstr)
       block
         character(kind=c_char,len=n+1), pointer :: s
         call c_f_pointer(cptr=cstr,fptr=s)
         allocate( svar, source = s(1:n) )
         nullify(s)
       end block
    else
       allocate( svar, source="")
    endif
  end subroutine ClassPosix_getenv

end Module ModulePOSIX
