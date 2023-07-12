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
Module ModuleCfgTemplates

  use, intrinsic :: ISO_FORTRAN_ENV
  use ModuleAstraConfigFile
  implicit none
  private
  
  public :: CheckOrTemplateCfgFiles

contains

  subroutine CheckOrTemplateCfgFiles( &
       TestCase         , &
       ccCfgFile        , &
       LuciaCfgFile     , &
       ScatciCfgFile    , &
       DaltonCfgFile    , &
       MoleculeCfgFile  )
    use ModuleAstraEnv
    character(len=*), intent(in)  :: TestCase
    character(len=*), intent(in)  :: ccCfgFile
    character(len=*), intent(in)  :: LuciaCfgFile
    character(len=*), intent(in)  :: ScatciCfgFile
    character(len=*), intent(in)  :: DaltonCfgFile
    character(len=*), intent(in)  :: MoleculeCfgFile

    character(len=:), allocatable :: TestDir
    logical                       :: Success = .TRUE.
    logical                       :: exist
    
    TestDir = AstraEnv%GetTestDir()
    
    call testCfgFile(       ccCfgFile, TestDir, TestCase,      ASTRA_CCFILE_DEF, Success )
    call testCfgFile(    LuciaCfgFile, TestDir, TestCase,    LUCIA_CFG_FILE_DEF, Success )
    call testCfgFile(   ScatciCfgFile, TestDir, TestCase,   SCATCI_CFG_FILE_DEF, Success )
    call testCfgFile(   DaltonCfgFile, TestDir, TestCase,   DALTON_CFG_FILE_DEF, Success )
    call testCfgFile( MoleculeCfgFile, TestDir, TestCase, MOLECULE_CFG_FILE_DEF, Success )

    if( .not. Success )then
       call WriteFramedMessage("Please check the config files and run astraSetup again")
       if(len_trim(TestDir)>0)then
          write(ERROR_UNIT,"(a)") " You can modify the current configuration files, "
          write(ERROR_UNIT,"(a)") " or you can delete them and fetch any of the following"
          write(ERROR_UNIT,"(a)") " other default cases from the runtime options"
          call Execute_Command_Line(&
               "find "//TestDir//"/* -mindepth 1 -maxdepth 1 -type d | " // &
               " awk -F/ '{print "//'"    "'//"$(NF-1)"//'"/"'//"$(NF)}'")
       endif
       stop
    endif

  end subroutine CheckOrTemplateCfgFiles
  
  subroutine testCfgFile( CfgFile, TestDir, TestCase, FileDefault, Success )
    character(len=*), intent(in)    :: CfgFile, TestDir, TestCase, FileDefault
    logical         , intent(inout) :: Success
    logical                         :: exist
    inquire( file = CfgFile, exist = exist )
    Success = Success .and. exist
    if(.not.exist)then
       call WriteFramedMessage("File "//CfgFile//" not found.")
       if(len_trim(TestDir)>0)then
          write(ERROR_UNIT,"(a)")" Copy template "
          call Execute_Command_Line( "cp "//TestDir//"/"//TestCase//"/"//&
               FileDefault//"  "//CfgFile)
       endif
    endif
  end subroutine testCfgFile
  
  subroutine WriteFramedMessage(strn)
    character(len=*), intent(in) :: strn
     write(ERROR_UNIT,"(a)") " "
     write(ERROR_UNIT,"(a)") " **********************************************************"
     write(ERROR_UNIT,"(a)") " * "//trim(adjustl(strn))
     write(ERROR_UNIT,"(a)") " **********************************************************"
  end subroutine WriteFramedMessage

end Module ModuleCfgTemplates
