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
!.. Contains constants and methods related to astra in relation
!   the its environment. This includes, for example, environmental
!   variables, data-tree structures, and specific system calls.
module ModuleAstraEnv

  use ModulePosix
  
  implicit none
  private

  character(len=*), private, parameter :: INSTALL_DIR = "ASTRA_DIR"
  character(len=*), private, parameter :: TEST_SUBDIR = "tests"
  character(len=*), private, parameter :: LOGO_SUBDIR = "documents/logos"
  character(len=*), private, parameter :: TEST_DEFAULT= "N2/default"
  character(len=*), private, parameter :: LOGO_DEFAULT= "Astra_StarWars.txt"

  !.. SINGLETON
  type, private :: ClassAstraEnv
   contains
     procedure, public :: GetInstallDir => ClassAstraEnv_GetInstallDir
     procedure, public :: GetTestDir    => ClassAstraEnv_GetTestDir
     procedure, public :: GetDefaultTest=> ClassAstraEnv_GetDefaultTest
     procedure, public :: GetLogoFile   => ClassAstraEnv_GetLogoFile
  end type ClassAstraEnv
  type(ClassAstraEnv), public :: AstraEnv 

contains

  function ClassAstraEnv_GetInstallDir( self ) result( sres )
    class(ClassAstraEnv), intent(in) :: self
    character(len=:),    allocatable :: sres
    call Posix%getenv( INSTALL_DIR, sres )
  end function ClassAstraEnv_GetInstallDir

  function ClassAstraEnv_GetTestDir( self ) result( sres )
    class(ClassAstraEnv), intent(in) :: self
    character(len=:),    allocatable :: sres
    character(len=:),    allocatable :: Dir
    Dir=self%GetInstallDir()
    allocate(sres,source=Dir//"/"//TEST_SUBDIR)
    deallocate(Dir)
  end function ClassAstraEnv_GetTestDir
  
  function ClassAstraEnv_GetDefaultTest( self ) result( sres )
    class(ClassAstraEnv), intent(in) :: self
    character(len=:),    allocatable :: sres
    allocate(sres,source=TEST_DEFAULT)
  end function ClassAstraEnv_GetDefaultTest
  
  function ClassAstraEnv_GetLogoFile( self ) result( sres )
    class(ClassAstraEnv), intent(in) :: self
    character(len=:),    allocatable :: sres
    character(len=:),    allocatable :: Dir
    Dir=self%GetInstallDir()
    allocate(sres,source=Dir//"/"//LOGO_SUBDIR//"/"//LOGO_DEFAULT)
    deallocate(Dir)
  end function ClassAstraEnv_GetLogoFile
  
end module ModuleAstraEnv
