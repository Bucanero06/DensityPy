module ModuleRTP_p

  implicit none

  private
  public :: GetRunTimeParameters

contains

  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( InpDir, OutDir, FileGeometry, &
       StepTime, StepWidth, FieldFile, Verbous, &
       nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax , i_excitation, i_epsilon)
    !
    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: InpDir
    character(len=:), allocatable, intent(out) :: OutDir
    character(len=:), allocatable, intent(out) :: FileGeometry
    real(kind(1d0))              , intent(out) :: StepTime
    real(kind(1d0))              , intent(out) :: StepWidth
    character(len=:), allocatable, intent(out) :: FieldFile
    logical                      , intent(out) :: Verbous
    integer                      , intent(out) :: nOmegas
    real(kind(1d0))              , intent(out) :: OmegaMin
    real(kind(1d0))              , intent(out) :: OmegaMax
    integer                      , intent(out) :: nTauOmegas
    real(kind(1d0))              , intent(out) :: TauOmegaMin
    real(kind(1d0))              , intent(out) :: TauOmegaMax
    !
    integer                      , intent(out) :: i_excitation
    integer                      , intent(out) :: i_epsilon
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Computes the FT of the dipole and charge response"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf
    !
    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help" , "Print Command Usage" )
    call List%Add( "-i"     , "Input  Dir" ,"CD_inp", "optional" )
    call List%Add( "-o"     , "Output Dir" ,"CD_out", "optional" )
    call List%Add( "-xyz"   , "Mol Geom File in Inp Dir" ,"geom.xyz", "optional" )
    call List%Add( "-stept" , "step time"     , 100.d0, "optional" )
    call List%Add( "-stepw" , "step width"    ,  10.d0, "optional" )
    call List%Add( "-field" , "Field File" ,"Field" , "optional" )
    call List%Add( "-v"     , "verbous"   )
    call List%Add( "-nw"    , "n  omegas"     , 101   , "optional" )
    call List%Add( "-wmin"  , "min omega"     ,   0.d0, "optional" )
    call List%Add( "-wmax"  , "max omega"     ,   1.d0, "optional" )
    call List%Add( "-ntw"    , "n tau omegas" , 101   , "optional" )
    call List%Add( "-twmin"  , "min tau omega",   0.d0, "optional" )
    call List%Add( "-twmax"  , "max tau omega",   1.d0, "optional" )
    call List%Add( "-s"    , "State to Excite Iteration", 0, "optional" )
    call List%Add( "-e"    , "1: Real, 2: Imaginary", 0, "optional" )
    call List%Add( "-xxx"   , "workaround")
    !
    call List%Parse()
    !
    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if
    !
    call List%Get( "-o",  strnBuf  )
    allocate(OutDir,source=trim(adjustl(strnBuf)))
    call List%Get( "-i",  strnBuf  )
    allocate(InpDir,source=trim(adjustl(strnBuf)))
    call List%Get( "-xyz",  strnBuf  )
    allocate(FileGeometry,source=InpDir//"/"//trim(adjustl(strnBuf)))
    call List%Get( "-stept", StepTime  )
    call List%Get( "-stepw", StepWidth )
    !
    call List%Get( "-s", i_excitation  )
    call List%Get( "-e", i_epsilon  )
    Verbous = List%Present("-v")
    call List%Get( "-field",  strnBuf  )
    allocate(FieldFile,source=trim(adjustl(strnBuf)))
    !
    call List%Get( "-nw"  , nOmegas  )
    call List%Get( "-wmin", OmegaMin )
    call List%Get( "-wmax", OmegaMax )
    call List%Get( "-ntw"  , nTauOmegas  )
    call List%Get( "-twmin", TauOmegaMin )
    call List%Get( "-twmax", TauOmegaMax )
    !
    call List%Free()
    !
  end subroutine GetRunTimeParameters

end module ModuleRTP_p















