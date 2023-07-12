! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
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
program ProgramAstraSetup

  use, intrinsic :: ISO_FORTRAN_ENV
  use ModulePOSIX
  use ModuleSystemUtils
  use ModuleErrorHandling
  use ModuleIO
  use ModuleString
  use ModuleConstants
  use ModuleAstraConfigFile
  use ModuleMainInterface
  use ModuleCfgTemplates

  implicit none

  !.. Run-time parameters
  character(len=:), allocatable :: AstraCfgFile
  character(len=:), allocatable :: TestCase
  logical                       :: RunDalton
  logical                       :: RunLuciaPIE
  logical                       :: RunLuciaHAIC
  logical                       :: RunLuciaTDMs
  logical                       :: RunScatci
  logical                       :: RunConvertTDMs
  logical                       :: RunConvertIntegrals
  !.. Config file parameters
  character(len=:), allocatable :: StorageDir
  character(len=:), allocatable :: QCDir
  character(len=:), allocatable :: ccCfgFile
  character(len=:), allocatable :: MoldenFile
  character(len=:), allocatable :: ScatciIntFile
  character(len=:), allocatable :: MoleculeCfgFile
  character(len=:), allocatable :: ScatciCfgFile
  character(len=:), allocatable :: DaltonCfgFile
  character(len=:), allocatable :: LuciaCfgFile
  !.. Local variables
  character(len=:), allocatable :: CommandDir, dch
  
  call GetRunTimeParameters( &
       AstraCfgFile, &
       TestCase    , &
       RunDalton   , &
       RunLuciaPIE , &
       RunLuciaHAIC, &
       RunLuciaTDMs, &
       RunScatci   , &
       RunConvertTDMs,&
       RunConvertIntegrals)
  
  !.. Parse the main config file, if present,
  !   or creates it with default entries
  call ParseAstraConfigFile( &
       AstraCfgFile     , &
       ccCfgFile        , &
       StorageDir       , &
       QCDir            , &
       MoldenFile       , &
       ScatciIntFile    , &
       LuciaCfgFile     , &
       ScatciCfgFile    , &
       DaltonCfgFile    , &
       MoleculeCfgFile  )

  !.. Check config files existence and compliance
  !   If files are broken or missing, it populates
  !   them with the default or the required options
  !   and it stops
  call CheckOrTemplateCfgFiles( &
       TestCase         , &
       ccCfgFile        , &
       LuciaCfgFile     , &
       ScatciCfgFile    , &
       DaltonCfgFile    , &
       MoleculeCfgFile  )

  !.. Move the QC config files to the QC directory
  write(OUTPUT_UNIT,"(a)") " Enter "//QCDir
  call Posix%getenv("PWD",CommandDir)
  call Execute_Command_line("mkdir -p "//QCDir)
  call Execute_Command_line("cp "// LuciaCfgFile   //" "//QCDir//"/LUCIA.INP")
  call Execute_Command_line("cp "// ScatciCfgFile  //" "//QCDir//"/SCATCI.INP")
  call Execute_Command_line("cp "// DaltonCfgFile  //" "//QCDir//"/DALTON.INP")
  call Execute_Command_line("cp "// MoleculeCfgFile//" "//QCDir//"/MOLECULE.INP")
  
  !.. Move into the QC directory and executes the required programs
  call Posix%chdir(QCDir)

  if(RunDalton)then
     write(OUTPUT_UNIT,"(a)") " Execute DALTON"
     call Execute_Command_line("dalton.x")
  endif
  
  if(RunLuciaPIE)then
     write(OUTPUT_UNIT,"(a)") " Run LUCIA to compute Parent-ion energies"
     call SetLuciaOption("LUCIA.INP","LUCIA_PIE.INP","ICDENS",", 0, 1")
     call Execute_Command_line("lucia.x < LUCIA_PIE.INP > Lucia_Parent_En.out 2>&1")
  endif
     
  if(RunLuciaHAIC)then
     write(OUTPUT_UNIT,"(a)") " Run LUCIA to compute < L | a H a+ | R >"
     call SetLuciaOption("LUCIA.INP","LUCIA_HAIC.INP","ICDENS",", -1, 1")
     call Execute_Command_line("lucia.x < LUCIA_HAIC.INP > Lucia_Loc_H.out  2>&1")
  endif

  if(RunLuciaTDMs)then
     write(OUTPUT_UNIT,"(a)") " Run LUCIA to compute TDMs < L | a+ a | R > and < L | a+ a+ a a | R >"
     call SetLuciaOption("LUCIA.INP","LUCIA_TDMs.INP","ICDENS",", 2")
     call Execute_Command_line("lucia.x < LUCIA_TDMs.INP > Lucia_TDM1-2B.out 2>&1")
  endif

  if(RunScatci)then
     write(OUTPUT_UNIT,"(a)") " Run SCATCI to compute hybrid integrals"
     call Execute_Command_line("rm -f "//ScatciIntFile)
     call Execute_Command_line("scatci_integrals SCATCI.INP > scatci_integrals.out 2>&1")
  endif

  !.. Return to the root directory
  call Posix%chdir(CommandDir)
  call execute_command_line("mkdir -p "//StorageDir//"/log")
  
  !.. Determines whether the preset command is run in debug mode
  call getDebugFlag( dch )

  if(RunConvertTDMs)then
     write(OUTPUT_UNIT,"(a)") " Convert the density matrices to reduced form in Astra"
     !call Execute_Command_line(dch//"astraConvertDensityMatrices"//&
     !     " -gif "//AstraCfgFile//" >  "//StorageDir//"/log/ConvertDensityMatrices.out 2>&1")
     write(*,*) dch//"astraConvertDensityMatrices"//&
          " -gif "//AstraCfgFile//" >  "//StorageDir//"/log/ConvertDensityMatrices.out 2>&1"
     call Execute_Command_line(dch//"astraConvertDensityMatrices"//&
          " -gif "//AstraCfgFile//" >  "//StorageDir//"/log/ConvertDensityMatrices.out 2>&1")
  endif
     
  if(RunConvertIntegrals)then
     write(OUTPUT_UNIT,"(a)") " Store the hybrid integrals in Astra"
     write(*,*) dch//"astraConvertIntegralsUKRmol"//&
          " -gif "//AstraCfgFile//" >  "//StorageDir//"/log/ConvertIntegrals.out 2>&1"
     call Execute_Command_line(dch//"astraConvertIntegralsUKRmol"//&
          " -gif "//AstraCfgFile//" >  "//StorageDir//"/log/ConvertIntegrals.out 2>&1")
  endif
     
contains

  subroutine GetDebugFlag(dch)
    character(len=:), allocatable, intent(inout) :: dch
    character(len=50) :: scom
    integer :: i
    call get_command_argument( 0, scom )
    scom=adjustl(scom)
    i=index(scom,"astraSetup")
    if(i<=0)then
       write(*,*) "Unrecognized program name: should be 'astraSetup'"
       stop
    endif
    if(i==1)then
       allocate(dch,source="")
       return
    endif
    if(scom(i-1:i-1)=="d")then
       allocate(dch,source="d")
    else
       allocate(dch,source="")
    end if
  end subroutine GetDebugFlag
  
  subroutine SetLuciaOption(File1,File2,Opt,val)
    character(len=*), intent(in) :: File1, File2, Opt, Val
    character(len=1000) :: line
    integer :: uid1, uid2,iostat
    open(newunit=uid1,file=trim(file1))
    open(newunit=uid2,file=trim(file2))
    do
       read(uid1,"(a)",iostat=iostat)line
       if(iostat/=0)exit
       line=adjustl(line)
       if(line(1:1)=="*")cycle
       if(index(line,Opt)>0)line=Opt//val
       write(uid2,"(a)")trim(line)
    enddo
    close(uid1)
    close(uid2)

  end subroutine SetLuciaOption
  
end program ProgramAstraSetup

