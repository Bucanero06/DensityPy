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
Module ModuleMainInterface

  use, intrinsic :: ISO_FORTRAN_ENV
  
  private

  public :: GetRunTimeParameters

contains
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters( &
       ProgramInputFile          , &
       Multiplicity              , &
       vSymLabel                 , &
       threshold                 )
    !
    !.. Level 0.0
    use ModuleString
    !.. Level 0.1
    use ModuleCommandLineParameterList
    use ModuleGroups
    !.. Level 2
    use ModuleAstraCredit
    !
    implicit none
    character(len=:), allocatable, intent(out)   :: ProgramInputFile
    integer                      , intent(out)   :: Multiplicity
    character(len=5), allocatable, intent(inout) :: vSymLabel(:)
    real(kind(1d0))              , intent(out)   :: threshold
    
    type( ClassCommandLineParameterList ) :: List
    character(len=:), allocatable :: tBuf
    character(len=100)            :: StrnBuf
    integer                       :: iT
    
    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Condition the overlap and Hamiltonian for the required symmetries"

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help", "Print Command Usage" )
    call List%Add( "-gif"  , "Astra Main Config File", "ASTRA.INP", "optional" )
    call List%Add( "-sym"  , "Irrep List (one mult), e.g.: 1Ag[,1B2u[,...]] or 2ALL", "1Ag", "required" )
    call List%Add( "-thr"  , "S diagonalization threshold", 1.d-8, "optional" )
 
    call List%Parse()

    call PrintAstraName()
    call PrintAstraCredit()
    if(List%Present("--help"))then
       call List%PrintUsage()
       stop
    end if

    call List%Get( "-gif"  , StrnBuf ); ProgramInputFile = trim( StrnBuf )

    call List%Get( "-sym",StrnBuf ); 
    read(StrnBuf(1:1),*) Multiplicity
    allocate(vSymLabel(nTokens(StrnBuf,",")))
    do iT=1,size(vSymLabel)
       call GetToken(StrnBuf,iT,tBuf,",")
       tBuf=adjustl(tBuf)
       vSymLabel(iT)=tBuf(2:)
    enddo
    call List%Get( "-thr",threshold ); 

    call List%Free()

    write(OUTPUT_UNIT,"(a)" )        " Program Input File  : "//ProgramInputFile
    write(OUTPUT_UNIT,"(a,i0)")      " Multiplicity        : ", Multiplicity
    write(OUTPUT_UNIT,"(*(a,x))")    " Spatial Symmetry    :",(trim(vSymLabel(iT)),iT=1,size(vSymLabel))
    write(OUTPUT_UNIT,"(a,x,e14.6)") " Threshold           :",threshold

  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
