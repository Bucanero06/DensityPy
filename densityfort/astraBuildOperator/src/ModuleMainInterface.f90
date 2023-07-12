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
Module ModuleMainInterface

  use, intrinsic :: ISO_FORTRAN_ENV
  
  private

  public :: GetRunTimeParameters

contains
  
  !> Reads the run time parameters from the command line.
  subroutine GetRunTimeParameters(      &
       ProgramInputFile             ,   &
       vOpLabel                     ,   &
       Multiplicity                 ,   &
       vKetSymLabel                 ,   &
       Verbous                      ,   &
       lSwitch                      )
    !
    use ModuleString
    use ModuleCommandLineParameterList
    use ModuleGroups
    use ModuleAstraCredit
    !
    implicit none
    character(len=:), allocatable, intent(out)   :: ProgramInputFile
    character(len=5), allocatable, intent(out)   :: vOpLabel(:)
    integer                      , intent(out)   :: Multiplicity
    character(len=5), allocatable, intent(inout) :: vKetSymLabel(:)
    logical                      , intent(out)   :: Verbous
    logical                      , intent(out)   :: lSWITCH(3,3,3)

    type( ClassCommandLineParameterList ) :: List
    character(len=:), allocatable :: tBuf
    character(len=100)            :: StrnBuf
    character(len=9)              :: StrnVec(3)
    integer                       :: LogoFormat, iTk, i, j, k, iOp
    
    character(len=*), parameter   :: PROGRAM_DESCRIPTION =&
         "Build the selected operator matrix, in the "//&
         "Close-Coupling basis."//&
         "Example:"//&
         "$> astraBuildOperator -op x -ketsym 3Ag"

    call List%SetDescription(PROGRAM_DESCRIPTION)
    call List%Add( "--help"  , "Print Command Usage" )
    call List%Add( "-gif"    , "Astra Main Config File", "ASTRA.INP", "optional" )
    call List%Add( "-op"     , "Operator List (e.g., 'H', 'x', 'H,S,y,z', ...)", "HS", "required" )
    call List%Add( "-ketsym" , "Ket Irrep List (one mult), e.g.: 1Ag[,1B2u[,...]] or 2ALL", "1Ag", "required" )
    call List%Add( "-v"      , "Verbous" )
    call List%Add( "-switchOver" , "overlap cases to compute: AA AI AE IA ...", "000000000", "optional" )
    call List%Add( "-switchMono" , "overlap cases to compute: AA AI AE IA ...", "000000000", "optional" )
    call List%Add( "-switchBiel" , "overlap cases to compute: AA AI AE IA ...", "000000000", "optional" )

    call List%Parse()

    call PrintAstraName()
    call PrintAstraCredit()
    if(List%Present("--help"))then
       call List%PrintUsage()
       
       stop
    end if

    call List%Get( "-gif"  , StrnBuf ); ProgramInputFile = trim( StrnBuf )
    call List%Get( "-op"   , StrnBuf )

    allocate(vOpLabel(nTokens(StrnBuf,",")))
    do iTk=1,size(vOpLabel)
       call GetToken(StrnBuf,iTk,tBuf,",")
       vOpLabel(iTk)=adjustl(tBuf)
    enddo
    
    call List%Get( "-ketsym",StrnBuf ); 
    read(StrnBuf(1:1),*) Multiplicity
       
       
    allocate(vKetSymLabel(nTokens(StrnBuf,",")))
    do iTk=1,size(vKetSymLabel)
       call GetToken(StrnBuf,iTk,tBuf,",")
       tBuf=adjustl(tBuf)
       vKetSymLabel(iTk)=tBuf(2:)
    enddo

    Verbous = List%Present("-v")
    

  !.. Test Overlap (one by one, and then add the rest progressively
  !   AA AI AE IA II IE EA EI EE
  !    *
  !    * * 
  !    .... 
  !.. Test Monoelectronic (one by one, and then add the rest progressively
  !   AA AI AE IA II IE EA EI EE
  !    *
  !    * * 
  !    .... 
  !.. Test Bielectronic (one by one, and then add the rest progressively
  !   AA AI AE IA II IE EA EI EE
  !    *
  !    * * 
  !    ....
    call List%Get( "-switchOver",StrnVec(1) ); 
    call List%Get( "-switchMono",StrnVec(2) ); 
    call List%Get( "-switchBiel",StrnVec(3) ); 
    
    do iOp = 1, 3
       k = 0
       do j = 1, 3
          do i = 1, 3
             k = k + 1
             lSwitch(i,j,iOp) = StrnVec(iOp)(k:k) == "1"
          enddo
       enddo
    enddo
    call List%Free()

    write(OUTPUT_UNIT,"(a)" )     " Program Input File  : "//ProgramInputFile
    write(OUTPUT_UNIT,"(*(a,x))") " Operator Name       :",(trim(vOpLabel(iTk)),iTk=1,size(vOpLabel))
    write(OUTPUT_UNIT,"(a,i0)")   " Multiplicity        : ", Multiplicity
    write(OUTPUT_UNIT,"(*(a,x))") " Ket Spatial Symmetry:",(trim(vKetSymLabel(iTk)),iTk=1,size(vKetSymLabel))
    write(OUTPUT_UNIT,"(a,l1)")   " Verbous             : ", Verbous

  end subroutine GetRunTimeParameters

end Module ModuleMainInterface
