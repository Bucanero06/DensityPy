Module ModuleMainInterface
  
  private

  public :: GetRunTimeParameters
  

contains
  
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters(      &
       lmax         , delete    )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleCommandLineParameterList
    !
    implicit none
    !
    integer                      , intent(out)   :: lmax
    integer, optional            , intent(out)   :: delete
    
    type( ClassCommandLineParameterList ) :: List
    character(len=100)            :: StrnBuf

    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Test some subroutines in ModuleAngularMomentum "

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help"  , "Print Command Usage" )
    call List.Add( "-lmax" , "Max Ang Mom", 4, "optional" )
    call List.Add( "-delete" , "nonsense variable", 4, "optional" )

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if
    

    call List.Get( "-lmax"  , lmax )

!    call List.Get( "-delete"  , delete )
    
    call List.Free()

  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
