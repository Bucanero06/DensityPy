module ModuleRTP_p

  implicit none

  private
  public :: GetRunTimeParameters


contains

  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( input_directory, output_directory, molecular_geometry_file, &
                n_times, t_min, t_max, FieldFile, Verbous, Weight_File, read_precomputed_weights_flag, &
                save_charge_migration_flag, ivorb, counted_number_of_orbitals, volume, dephasing_factor, relaxation_factor, &
                bath_temperature, i_excitation, i_epsilon)
    !
    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: input_directory
    character(len=:), allocatable, intent(out) :: output_directory
    character(len=:), allocatable, intent(out) :: molecular_geometry_file
    integer                      , intent(out) :: n_times
    real(kind(1d0))              , intent(out) :: t_min
    real(kind(1d0))              , intent(out) :: t_max
    character(len=:), allocatable, intent(out) :: FieldFile
    character(len=:), allocatable, intent(out) :: Weight_File
    logical                      , intent(out) :: read_precomputed_weights_flag
    logical                      , intent(out) :: Verbous
    logical                      , intent(out) :: save_charge_migration_flag
    integer         , allocatable, intent(out) :: ivorb(:)
    real(kind(1d0))              , intent(out) :: volume
    real(kind(1d0))              , intent(out) :: dephasing_factor
    real(kind(1d0))              , intent(out) :: relaxation_factor
    real(kind(1d0))              , intent(out) :: bath_temperature
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
    integer, intent(out) :: counted_number_of_orbitals
    integer ::pos1, pos2 , i, number_of_orbitals

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help" , "Print Command Usage" )
    call List%Add( "-i"     , "Input  Dir" ,"CD_inp", "optional" )
    call List%Add( "-o"     , "Output Dir" ,"CD_out", "optional" )
    call List%Add( "-xyz"   , "Mol Geom File in Inp Dir" ,"geom.xyz", "optional" )
    call List%Add( "-nt"    , "n times"    ,    101 , "optional" )
    call List%Add( "-t_min"  , "min time"   ,    0.d0, "optional" )
    call List%Add( "-t_max"  , "max time"   ,  200.d0, "optional" )
    call List%Add( "-field" , "Field File" ,"Field" , "optional" )
    call List%Add( "-w"     , "Becke's Weights File" ,"Weights_File" , "optional" )
    call List%Add( "-v"     , "verbous"   )
    call List%Add( "-sden"  , "Save Charge Density to file (time consuming!)")
    call List%Add( "-iorb"  , "absolute index active orbitals ", "7,8,9,10,11,12,13", "optional" )
    call List%Add( "-vol"   , "volume of unit for grid", 1.d-3, "optional" )
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
    allocate(output_directory,source=trim(adjustl(strnBuf)))
    call List%Get( "-i",  strnBuf  )
    allocate(input_directory,source=trim(adjustl(strnBuf)))
    call List%Get( "-xyz",  strnBuf  )
    allocate(molecular_geometry_file,source=input_directory//"/"//trim(adjustl(strnBuf)))
    call List%Get( "-nt",  n_times  )
    call List%Get( "-t_min", t_min  )
    call List%Get( "-t_max", t_max  )
    call List%Get( "-vol", volume  )
    call List%Get( "-rf", relaxation_factor  )
    call List%Get( "-df", dephasing_factor  )
    call List%Get( "-bath", bath_temperature  )
    !
    call List%Get( "-s", i_excitation  )
    call List%Get( "-e", i_epsilon  )
    Verbous = List%Present("-v")

    save_charge_migration_flag = List%Present("-sden")

    read_precomputed_weights_flag = List%Present("-w")
    call List%Get( "-w",  strnBuf  )

    allocate(Weight_File,source=trim(adjustl(strnBuf)))

    call List%Get( "-field",  strnBuf  )

    allocate(FieldFile,source=trim(adjustl(strnBuf)))

    !..Gets String with Orbitals
    pos1=1
    number_of_orbitals=0
    call List%Get( "-iorb", orbital_string)
    counted_number_of_orbitals=1
    do
       pos1=index(orbital_string,",")
       if(pos1<=0)exit
       counted_number_of_orbitals=counted_number_of_orbitals+1
       orbital_string(pos1:pos1)=" "
    enddo
    allocate(ivorb(counted_number_of_orbitals))
    read(orbital_string,*) (ivorb(i),i=1,counted_number_of_orbitals)

    call List%Free()
    !
  end subroutine GetRunTimeParameters



end module ModuleRTP_p















