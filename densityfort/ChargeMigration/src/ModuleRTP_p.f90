module ModuleRTP_p

  implicit none

  private
  public :: GetRunTimeParameters


contains

  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( InpDir, OutDir, FileGeometry, &
                nTimes, Tmin, Tmax, FieldFile, Verbous, Weight_File, Weight_File_flag, &
                SaveDensity, ivorb, OrbNumber, Volume, dephasing_factor, RELAXATION_FACTOR, &
                BATH_TEMPERATURE, i_excitation, i_epsilon)
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
    integer                      , intent(out) :: nTimes
    real(kind(1d0))              , intent(out) :: Tmin
    real(kind(1d0))              , intent(out) :: Tmax
    character(len=:), allocatable, intent(out) :: FieldFile
    character(len=:), allocatable, intent(out) :: Weight_File
    logical                      , intent(out) :: Weight_File_flag
    logical                      , intent(out) :: Verbous
    logical                      , intent(out) :: SaveDensity
    integer         , allocatable, intent(out) :: ivorb(:)
    real(kind(1d0))              , intent(out) :: Volume
    real(kind(1d0))              , intent(out) :: dephasing_factor
    real(kind(1d0))              , intent(out) :: RELAXATION_FACTOR
    real(kind(1d0))              , intent(out) :: BATH_TEMPERATURE
    !
    integer                      , intent(out) :: i_excitation
    integer                      , intent(out) :: i_epsilon
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Computes the Charge Density on a spatial grid from "//&
         "tabulated orbitals, The Density Matrices and Transition Density Matrices"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf

    CHARACTER(10000)  :: orbital_string
    integer, intent(out) :: OrbNumber
    integer ::pos1, pos2 , i, nOrb

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help" , "Print Command Usage" )
    call List%Add( "-i"     , "Input  Dir" ,"CD_inp", "optional" )
    call List%Add( "-o"     , "Output Dir" ,"CD_out", "optional" )
    call List%Add( "-xyz"   , "Mol Geom File in Inp Dir" ,"geom.xyz", "optional" )
    call List%Add( "-nt"    , "n times"    ,    101 , "optional" )
    call List%Add( "-tmin"  , "min time"   ,    0.d0, "optional" )
    call List%Add( "-tmax"  , "max time"   ,  200.d0, "optional" )
    call List%Add( "-field" , "Field File" ,"Field" , "optional" )
    call List%Add( "-w"     , "Becke's Weights File" ,"Weights_File" , "optional" )
    call List%Add( "-v"     , "verbous"   )
    call List%Add( "-sden"  , "Save Charge Density to file (time consuming!)")
    call List%Add( "-iorb"  , "absolute index active orbitals ", "7,8,9,10,11,12,13", "optional" )
    call List%Add( "-vol"   , "Volume of unit for grid", 1.d-3, "optional" )
    call List%Add( "-bath"  , "BATH TEMPERATURE", 327175.d-2, "optional" )
    call List%Add( "-rf"    , "RELAXATION FACTOR", 1.d-3, "optional" )
    call List%Add( "-df"    , "DEPHASING FACTOR", 1.d-3, "optional" )
    call List%Add( "-s"    , "State to Excite Iteration", 0, "optional" )
    call List%Add( "-e"    , "1: Real, 2: Imaginary", 0, "optional" )

    call List%Add( "-xxx"   , "workaround")

    call List%Parse()

    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-o",  strnBuf  )
    allocate(OutDir,source=trim(adjustl(strnBuf)))
    call List%Get( "-i",  strnBuf  )
    allocate(InpDir,source=trim(adjustl(strnBuf)))
    call List%Get( "-xyz",  strnBuf  )
    allocate(FileGeometry,source=InpDir//"/"//trim(adjustl(strnBuf)))
    call List%Get( "-nt",  nTimes  )
    call List%Get( "-tmin", Tmin  )
    call List%Get( "-tmax", Tmax  )
    call List%Get( "-vol", Volume  )
    call List%Get( "-rf", RELAXATION_FACTOR  )
    call List%Get( "-df", dephasing_factor  )
    call List%Get( "-bath", BATH_TEMPERATURE  )
    !
    call List%Get( "-s", i_excitation  )
    call List%Get( "-e", i_epsilon  )
    Verbous = List%Present("-v")

    SaveDensity = List%Present("-sden")

    Weight_File_flag = List%Present("-w")
    call List%Get( "-w",  strnBuf  )

    allocate(Weight_File,source=trim(adjustl(strnBuf)))

    call List%Get( "-field",  strnBuf  )

    allocate(FieldFile,source=trim(adjustl(strnBuf)))

    !..Gets String with Orbitals
    pos1=1
    nOrb=0
    call List%Get( "-iorb", orbital_string)
    OrbNumber=1
    do
       pos1=index(orbital_string,",")
       if(pos1<=0)exit
       OrbNumber=OrbNumber+1
       orbital_string(pos1:pos1)=" "
    enddo
    allocate(ivorb(OrbNumber))
    read(orbital_string,*) (ivorb(i),i=1,OrbNumber)

    call List%Free()
    !
  end subroutine GetRunTimeParameters



end module ModuleRTP_p















