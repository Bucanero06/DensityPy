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
!>  Insert Doxygen comments here
!!  
!!  This program builds an "arbitray" operator between
!!  symmetric close-coupling spaces of a molecule, within
!!  D2h and its subgroups.
!! 
program BuildSymmetricElectronicSpace

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleString
  use ModuleGroups
  use ModuleAstraConfigFile
  use ModuleOrbitalBasis
  use ModuleXlm
  use ModuleMolecularGeometry
  use ModuleIntegrals
  use ModuleDensityMatrices
  use ModuleSymESpace
  use ModuleESpace
  
  use ModuleMainInterface

  implicit none

  !.. Run time parameters
  character(len=:), allocatable :: AstraConfigFile
  character(len=5), allocatable :: vOpLabel(:)
  integer                       :: Multiplicity
  character(len=5), allocatable :: vKetSymLabel(:)
  logical                       :: Verbous
  logical                       :: lSWITCH(3,3,3)
  
  !.. Config file parameters
  character(len=:), allocatable :: StorageDir
  character(len=:), allocatable :: QCDir
  character(len=:), allocatable :: ccConfigFile
  character(len=:), allocatable :: MoleculeConfigFile

  type(ClassESpace)             :: Space
  character(len=:), allocatable :: BraSymLabel
  type(ClasssymCCSpace), pointer :: BraSymSpace, KetSymSpace
  type(ClassSESSESBlock)        :: OperatorMat

  type(ClassIrrep), pointer     :: KetIrrep, BraIrrep, OpIrrep
  character(len=*), parameter   :: IDLabel = " "
  logical                       :: TaskList(N_SESSES_ID)
  character(len=:), allocatable :: SpaceDir
  integer                       :: iOp, iSym
  integer         , allocatable :: vnActive(:), vnInActive(:)
  real(kind(1d0))               :: TotalNuclearCharge, TotalNElectrons

  ! !.. DEBUG PARAMETERS
  ! integer ::  SWITCH(3,3,3)
  ! SWITCH = reshape( (/ &
  !      !0,0,0, 0,0,0, 0,0,0, &
  !      1,1,1, 1,1,1, 1,1,1, &
  !      1,1,1, 1,1,1, 1,1,1, &
  !      0,0,0, 0,0,0, 0,0,0 /), (/3,3,3/) )
  ! lSWITCH = SWITCH > 0

  call GetRunTimeParameters( &
       AstraConfigFile     , &
       vOpLabel            , &
       Multiplicity        , &
       vKetSymLabel        , &
       Verbous             , &
       lSwitch             )

  call Set_SESSES_switch( "build", lSWITCH )
  
  call ParseAstraConfigFile( &
       AstraConfigFile     , &
       ccConfigFile        , &
       StorageDir          , &
       MoleculeCfgFile = MoleculeConfigFile  )

  !.. Setup the Electronic Space
  
  call Space%SetMultiplicity( Multiplicity )
  call Space%SetRootDir     ( StorageDir )
  call Space%ParseConfigFile( ccConfigFile )
  SpaceDir=Space%GetStorageDir()
  if(Verbous) call Space%Show()
  call MolecularGeometry%ParseFile(MoleculeConfigFile)
  TotalNuclearCharge = MolecularGeometry%getCharge()
  TotalNElectrons    = TotalNuclearCharge - space%GetPionCharge()

  !.. Setup the global group
  call GlobalGroup%Init(Space%GetGroupName())

  !.. Setup the global Xlm set
  call GlobalXlmSet%Init(GlobalGroup,Space%GetLmax())

  !.. Load the transition density matrices
  write(OUTPUT_UNIT,"(a)") " Loading Transition Density Matrices"
  call TDM_Manager%SetStorageDir( StorageDir )
  call TDM_Manager%SetParentIons( Space%GetPionList() )
  call TDM_Manager%Init()
  call TDM_Manager%IO( "1", "Read" )
  call TDM_Manager%IO( "2", "Read" )
  call TDM_Manager%IO( "H", "Read" )
  call TDM_Manager%IO( "E", "Read" )
  call TDM_Manager%GetNActiveEls(vnInactive,vnActive)

  !.. Load the Integrals
  write(OUTPUT_UNIT,"(a)") " Loading Mono and Bielectronic Integrals"
  call GlobalIntegral%SetStorage( AddSlash(StorageDir) )
  call GlobalIntegral%Setlmax( Space%GetLmax())
  call GlobalIntegral%SetnInactive( vnInactive )
  call GlobalIntegral%ReadFromFile()
  call GlobalIntegral%ConvertHinHtilde(TotalNElectrons, vnInactive) 


  !.. Load the parameters of the orbital basis
  !   and communicate to the Orbital Basis the
  !   vector of active orbitals (which is an information specific
  !   to the CI calculations, rather than inherent to the basis or
  !   to the integrals
  call OrbitalBasis%SetStorageDir(StorageDir)
  call OrbitalBasis%Load()
  call OrbitalBasis%SetNactive(vnInactive,vnActive)
 
  do iOp = 1, size(vOpLabel)

     call SetStringToUppercase( vOpLabel(iOp) )

 
     
     if(index(vKetSymLabel(1),"ALL")>0)then
        deallocate(vKetSymLabel)
        allocate(vKetSymLabel(GlobalGroup%GetNIrreps()))
        do iSym=1,GlobalGroup%GetNIrreps()
           vKetSymLabel(iSym)=GlobalGroup%GetIrrName(iSym)
        enddo
     endif

     do iSym = 1, size(vKetSymLabel)
        !.. Determines the bra and operator irreps
        call Space%CheckSymmetry( trim(vKetSymLabel(iSym)) )
        KetIrrep    => GlobalGroup%GetIrrep( trim(vKetSymLabel(iSym)) )
        OpIrrep     => GetOpIrrep( GlobalGroup, trim(vOpLabel(iOp)) )
        BraIrrep    => OpIrrep * KetIrrep
        BraSymLabel =  BraIrrep%GetName()
        write(OUTPUT_UNIT,"(a)",advance="no") " < "//BraSymLabel//" | "//trim(vOpLabel(iOp))//&
             " | "//trim(vKetSymLabel(iSym))//" > ..."
        !.. Set up the symmetric Bra and Ket electronic spaces
        BraSymSpace => Space%GetSymElectSpace(       BraSymLabel        ) 
        KetSymSpace => Space%GetSymElectSpace( trim(vKetSymLabel(iSym)) )
        call OperatorMat%Free()
        call OperatorMat%Init( BraSymSpace, KetSymSpace, trim(vOpLabel(iOp)), SpaceDir, IDLabel )
        TaskList( SESSES_LOAD  ) = .FALSE.
        TaskList( SESSES_BUILD ) = .TRUE.
        TaskList( SESSES_SAVE  ) = .TRUE.
        TaskList( SESSES_FREE  ) = .TRUE.

        call OperatorMat%Driver( TaskList )
        
        write(OUTPUT_UNIT,"(a)") "...V"
        
     enddo
     
  enddo

  stop

end program BuildSymmetricElectronicSpace
