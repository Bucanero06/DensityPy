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
module ModuleASTRAConfigFile
  
  use, intrinsic :: ISO_FORTRAN_ENV
  
  implicit none
  private

  character(len=*), public, parameter :: ASTRA_CFG_FILE_DEF    = "ASTRA.INP"
  character(len=*), public, parameter :: ASTRA_CCFILE_DEF      = "CLSCPLNG.INP"
  character(len=*), public, parameter :: SCATCI_CFG_FILE_DEF   = "SCATCI.INP"
  character(len=*), public, parameter :: LUCIA_CFG_FILE_DEF    = "LUCIA.INP"
  character(len=*), public, parameter :: DALTON_CFG_FILE_DEF   = "DALTON.INP"
  character(len=*), public, parameter :: MOLECULE_CFG_FILE_DEF = "MOLECULE.INP"
  character(len=*), public, parameter :: ASTRA_STORAGE_DEF     = "store"
  character(len=*), public, parameter :: ASTRA_QCDIR_DEF       = "QC"
  character(len=*), public, parameter :: MOLDEN_FILE_DEF       = "molden.inp"
  character(len=*), public, parameter :: SCATCI_INT_FILE_DEF   = "moints"
  
  character(len=*)        , parameter :: ASTRA_CCFILE          = "cc_cfg_file"
  character(len=*)        , parameter :: SCATCI_CFG_FILE       = "scatci_cfg_file"
  character(len=*)        , parameter :: LUCIA_CFG_FILE        = "lucia_cfg_file"
  character(len=*)        , parameter :: DALTON_CFG_FILE       = "dalton_cfg_file"
  character(len=*)        , parameter :: MOLECULE_CFG_FILE     = "molecule_cfg_file"
  character(len=*)        , parameter :: ASTRA_STORAGE         = "storage_dir"
  character(len=*)        , parameter :: ASTRA_QCDIR           = "qc_dir"
  character(len=*)        , parameter :: MOLDEN_FILE           = "molden_file"
  character(len=*)        , parameter :: SCATCI_INT_FILE       = "scatci_int_file"
  
  public :: ParseAstraConfigFile
  
contains

  !> Parses Astra General Config file.
  subroutine ParseAstraConfigFile( &
       AstraCfgFile    , &
       ccCfgFile       , &
       StorageDir      , &
       QCDir           , &
       MoldenFile      , &
       ScatciIntFile   , &
       LuciaCfgFile    , &
       ScatciCfgFile   , &
       DaltonCfgFile   , &
       MoleculeCfgFile )
    !
    use ModuleParameterList
    !
    implicit none
    character(len=*)                       , intent(in)  :: AstraCfgFile
    character(len=:)          , allocatable, intent(out) :: ccCfgFile
    character(len=:)          , allocatable, intent(out) :: StorageDir
    character(len=:), optional, allocatable, intent(out) :: QCDir
    character(len=:), optional, allocatable, intent(out) :: MoldenFile
    character(len=:), optional, allocatable, intent(out) :: ScatciIntFile
    character(len=:), optional, allocatable, intent(out) :: LuciaCfgFile
    character(len=:), optional, allocatable, intent(out) :: ScatciCfgFile
    character(len=:), optional, allocatable, intent(out) :: DaltonCfgFile
    character(len=:), optional, allocatable, intent(out) :: MoleculeCfgFile

    character(len=100)       :: strnBuf
    type(ClassParameterList) :: List
    logical                  :: exist

    call List%Add( ASTRA_CCFILE     , ASTRA_CCFILE_DEF     , "optional" )
    call List%Add( ASTRA_STORAGE    , ASTRA_STORAGE_DEF    , "optional" )
    call List%Add( ASTRA_QCDIR      , ASTRA_QCDIR_DEF      , "optional" )
    call List%Add( MOLDEN_FILE      , MOLDEN_FILE_DEF      , "optional" )
    call List%Add( SCATCI_INT_FILE  , SCATCI_INT_FILE_DEF  , "optional" )
    call List%Add( LUCIA_CFG_FILE   , LUCIA_CFG_FILE_DEF   , "optional" )
    call List%Add( SCATCI_CFG_FILE  , SCATCI_CFG_FILE_DEF  , "optional" )
    call List%Add( DALTON_CFG_FILE  , DALTON_CFG_FILE_DEF  , "optional" )
    call List%Add( MOLECULE_CFG_FILE, MOLECULE_CFG_FILE_DEF, "optional" )
    
    inquire( file = AstraCfgFile, exist = exist )
    if(.not.exist)then
       write(ERROR_UNIT,"(a)") " " 
       write(ERROR_UNIT,"(a)") " ********************************************************"
       write(ERROR_UNIT,"(a)") " * File "//AstraCfgFile//" not found. Generate template  "
       write(ERROR_UNIT,"(a)") " ********************************************************"
       write(ERROR_UNIT,"(a)") " " 
       call list%WriteDefault( AstraCfgFile )
    else
       call List%Parse( AstraCfgFile )
    endif
    
    call List%Get( ASTRA_CCFILE, strnBuf )
    allocate( ccCfgFile, source = trim(strnBuf) )
    write(OUTPUT_UNIT,"(a)")    " CC Config File      : "//ccCfgFile
    call List%Get( ASTRA_STORAGE, strnBuf )
    allocate( StorageDir, source = trim(strnBuf) )
    write(OUTPUT_UNIT,"(a)")    " Astra Storage Dir   : "//StorageDir
    if(present(QCDir))then
       call List%Get( ASTRA_QCDIR, strnBuf )
       allocate( QCDir, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " QC Storage Dir      : "//QCDir
    endif
       if(present(MoldenFile))then
       call List%Get( MOLDEN_FILE, strnBuf )
       allocate( MoldenFile, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " Molden File         : "//MoldenFile
    endif
    if(present(LuciaCfgFile))then
       call List%Get( LUCIA_CFG_FILE  , strnBuf )
       allocate( LuciaCfgFile, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " Lucia  Config File  : "//LuciaCfgFile
    endif
    if(present(ScatciIntFile))then
       call List%Get( SCATCI_INT_FILE  , strnBuf )
       allocate( ScatciIntFile, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " Scatci Int File     : "//ScatciIntFile
    endif
    if(present(ScatciCfgFile))then
       call List%Get( SCATCI_CFG_FILE  , strnBuf )
       allocate( ScatciCfgFile, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " Scatci Config File  : "//ScatciCfgFile
    endif
    if(present(DaltonCfgFile))then
       call List%Get( DALTON_CFG_FILE  , strnBuf )
       allocate( DaltonCfgFile, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " Dalton Config File  : "//DaltonCfgFile
    endif
    if(present(MoleculeCfgFile))then
       call List%Get( MOLECULE_CFG_FILE  , strnBuf )
       allocate( MoleculeCfgFile, source = trim(strnBuf) )
       write(OUTPUT_UNIT,"(a)") " Molecule Config File: "//MoleculeCfgFile
    endif

  end subroutine ParseAstraConfigFile
  
end module ModuleASTRAConfigFile
