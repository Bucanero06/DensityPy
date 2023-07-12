! CONFIDENTIAL
!!>  Insert Doxygen comments here
!!  
!!  This program read, write and / or generates the density matrices and does
!!  the several necessary transformations (spin coupling etc).
!! 
program ConvertDensityMatricesSTEX

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleMainInterface
  use ModuleOrbitalBasis
  use ModuleAstraConfigFile
  use ModuleString
  use ModuleIO
  use ModuleGroups
  use ModuleESpace
  use ModuleDensityMatrices
  use ModuleParentIons
  use ModuleIntegrals
  use ModuleMatrix



  implicit none

  !.. Run time parameters
  character(len=:), allocatable :: AstraConfigFile
  logical                       :: CheckConversion
  logical                       :: UseSTEX

  !.. Config file parameters
  character(len=:), allocatable :: StorageDir
  character(len=:), allocatable :: ccConfigFile
  character(len=:), allocatable :: QCDir

  type(ClassESpace)    :: Space
  type(ClassGroup)    , pointer :: Group
  type(ClassParentIon), pointer :: ParentIonList(:)

  !> index of the irrep in ASTRA order corresponding to
  !! a given index of the same irrep in DALTON order
  integer, allocatable :: ASTRAIrr_from_DALTONIrr(:)
  !> Number of inactive molecular orbitals per irreducible representation (valid for alpha and beta) 
  integer, allocatable :: ninactive(:)
  !> Number of active molecular orbitals per irreducible representation (valid for alpha and beta) 
  integer, allocatable :: nactel(:)
  !> Number of pairs of orbitals per irreducible representation (first index) with alpha-beta (second index equal to 1)
  !  and alpha-alpha or beta-beta spin projections (second index equal to 2) 
  integer, allocatable :: nactpairs(:,:)
  !> Number of different spin projections which are availables for each ion (size(vPInMs2) = nIons)
  integer, allocatable :: vPInMs2(:)
  !> As a function of the index wich runs in the number of available spin projections (first index) and the ion number (second index),
  !  mPIlistMs2 gives the spin projection value
  integer, allocatable :: mPIlistMs2(:,:)
  !> As a function of the index wich runs in the number of available spin projections (first index) and the ion number (second index),
  !> mLUCIAIndexPIMs2 gives the corresponding absolute state LUCIA index
  integer, allocatable :: mLUCIAIndexPIMs2(:,:)
  !> number of occupied orbitals per irrep in the reference configuration
  !! if it is defined. This quantity is only relevant in case the STEX
  !! run-time option is specified.
  integer, allocatable :: nRefOccupied(:)

  call GetRunTimeParameters( AstraConfigFile, CheckConversion, UseSTEX )

  call ParseAstraConfigFile( &
       AstraConfigFile     , &
       ccConfigFile        , &
       StorageDir          , &
       QCDir               )

  call Space%SetRootDir     ( StorageDir   )
  call Space%ParseConfigFile( ccConfigFile )

  Group         => Space%GetGroup()
  ParentIonList => Space%GetPionList()

  call GlobalGroup%Init(Group%GetName())

  
  !.. Reads the order of the irreps in Dalton
  call ReadDALTONIrreps( ASTRAIrr_from_DALTONIrr, QCDir )

  !.. Reads the List of LUCIA states for each parent ion
  !   as well as the number of active electrons per irrep
  !..
  call ReadLUCIAIonsAndNactel( ASTRAIrr_from_DALTONIrr, ParentIonList, &
       nactel, nactpairs, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir )

  !.. Pass the list of parent ions to the Transition-Density-Matrix Manager
  call TDM_Manager%SetParentIons( ParentIonList )
  call TDM_Manager%SetNActiveEls( ninactive, nactel )
  call TDM_Manager%SetStorageDir( StorageDir )
  call TDM_Manager%Init()

  nRefOccupied = Space%GetnRefOccupied()

  call GlobalIntegral%SetStorage( AddSlash(StorageDir) )
!  call GlobalIntegral%Setlmax( lmax )
  call GlobalIntegral%ReadFromFile()
  
  call SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, "1", QCDir, nRefOccupied )
  call SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, "H", QCDir, nRefOccupied )
  call SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, "E", QCDir, nRefOccupied )
  call SetupTDM2( ParentIonList, nactel, nactpairs, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir, nRefOccupied )

  !.. Save the TDM in the format of astra
  call TDM_Manager%IO(  "1", "Write" )
  call TDM_Manager%IO(  "2", "Write" )
  call TDM_Manager%IO(  "H", "Write" )
  call TDM_Manager%IO(  "E", "Write" )
  call TDM_Manager%Free()

  stop


contains


  subroutine SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, &
       ID_CHAR, QCDir, nRefOccupied )

    type(ClassParentIon), intent(in) :: ParentIonList(:)
    integer             , intent(in) :: ASTRAIrr_from_DALTONIrr(:)
    integer             , intent(in) :: nactel(:), vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    character           , intent(in) :: ID_CHAR
    character(len=*)    , intent(in) :: QCDir
    integer, optional   , intent(in) :: nRefOccupied(:)

    integer                          :: iIonA, iIonB, MultA, MultB, iStateA, iStateB, Sigma2, nIons
    integer                          :: numA, numB, iSymA, iSymB
    integer                          :: nIrreps, nOrbMax, iIrr1, nr, nc, i, j
    integer        , allocatable     :: vnr(:), vnc(:)
    logical        , allocatable     :: Irr1Loaded(:) 
    real(kind(1d0)), allocatable     :: Rho(:,:,:,:)

    nIons=size(ParentIonList)
    nIrreps=GlobalGroup%GetNirreps()
    allocate(Irr1Loaded(nIrreps))
    Irr1Loaded=.FALSE.

    nOrbMax=maxval(nactel)

    allocate(Rho(nOrbMax,nOrbMax,nIrreps,2))
    Rho=0.d0
    allocate(vnr(nIrreps),vnc(nIrreps))

    do iIonA = 1, nIons
       MultA = ParentIonList(iIonA)%GetMultiplicity() 
       iSymA = GlobalGroup%GetIrrepIndex( ParentIonList(iIonA)%GetIrrep() )
       NumA  = ParentIonList(iIonA)%GetNumber()

       do iIonB = 1, nIons
          MultB = ParentIonList(iIonB)%GetMultiplicity()
          iSymB = GlobalGroup%GetIrrepIndex( ParentIonList(iIonB)%GetIrrep() )
          NumB  = ParentIonList(iIonB)%GetNumber()
          if(abs(MultB-MultA)>2)cycle

          call GetIndexesWithMaxCommonSigma(iIonA,iIonB,vPInMs2,mPIlistMs2,mLUCIAIndexPIMs2,&
               iStateA, iStateB, Sigma2)

          !.. Note that, in the case of static exchange, all the ions obtained above have Sigma2 = 1

          if(    ID_CHAR.is."1")then
             call SetupSTEXRho(iSymA,NumA,iSymB,NumB,nactel,nRefOccupied,Irr1Loaded,vnr,vnc,Rho)
          elseif(ID_CHAR.is."H")then
             call SetupSTEXHab(iSymA,NumA,iSymB,NumB,nactel,nRefOccupied,Irr1Loaded,vnr,vnc,Rho)
          endif
          !call LoadRho(iStateA,iStateB,ASTRAIrr_from_DALTONIrr,Irr1Loaded,vnr,vnc,Rho,ID_CHAR,QCDir)


          do iIrr1=1,nIrreps
             if(.not.Irr1Loaded(iIrr1))cycle
             nr=vnr(iIrr1)
             nc=vnc(iIrr1)
             if(    ID_CHAR.is."1")then
                call TDM_Manager%SetTDM1( iIonA, iIonB, iIrr1, Sigma2, nr, nc, &
                     Rho(1,1,iIrr1,1), nOrbMax, Rho(1,1,iIrr1,2), nOrbMax )
             elseif(ID_CHAR.is."H")then
                call TDM_Manager%SetHAIC( iIonA, iIonB, iIrr1, Sigma2, nr, nc, &
                     Rho(1,1,iIrr1,1), nOrbMax, Rho(1,1,iIrr1,2), nOrbMax )
             elseif( (ID_CHAR.is."E") .and. iIonA == iIonB )then
                call TDM_Manager%SetPIEN( iIonA, Rho(1,1,iIrr1,1) )
             endif
                
          enddo

       enddo
    enddo
  end subroutine SetupTDM1


  !   
  !   2. Store the TDM in the ASTRA format, via the methods of the TDMManager methods.
  !..
  subroutine SetupTDM2( ParentIonList, nactel, nactpairs, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir, nRefOccupied )
    use ModuleMatrix

    type(ClassParentIon), intent(in) :: ParentIonList(:)
    integer             , intent(in) :: ASTRAIrr_from_DALTONIrr(:)
    integer             , intent(in) :: nactel(:), nactpairs(:,:)
    integer             , intent(in) :: vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    character(len=*)    , intent(in) :: QCDir
    integer, optional   , intent(in) :: nRefOccupied(:)

    integer, parameter :: ALPHA_BETA = 1
    integer, parameter :: ALPHA_ALPHA= 2
    integer, parameter :: BETA_BETA  = 3
    !
    integer :: iIonA, iIonB, MultA, MultB, Sigma2_case1, Sigma2_case2, nIons
    integer :: numA, numB, iSymA, iSymB
    integer :: nIrreps, iIrr1, iIrr2, iIrr3, iIrr4, n1, n2, n3, n4, iType, iMs
    integer :: iStateA, iStateB, Sigma2
    !
    !   PI => one of the possible spin projections for the ion
    !.. Index 1 - 3 = symmetry first - third orbitals
    !   Index 4 = type  ( 1, 2, 3 for alpha-beta, alpha-alpha, beta-beta )
    type(ClassMatrix4D), pointer :: Pi(:,:,:,:)
    type(ClassIrrep)   , pointer :: IrrepAB

    nIrreps=GlobalGroup%GetNirreps()

    nIons=size(ParentIonList)
    do iIonA = 1, nIons
       MultA = ParentIonList(iIonA)%GetMultiplicity()
       iSymA = GlobalGroup%GetIrrepIndex( ParentIonList(iIonA)%GetIrrep() )
       NumA  = ParentIonList(iIonA)%GetNumber()
       do iIonB = 1, nIons
          MultB = ParentIonList(iIonB)%GetMultiplicity()
          iSymB = GlobalGroup%GetIrrepIndex( ParentIonList(iIonB)%GetIrrep() )
          NumB  = ParentIonList(iIonB)%GetNumber()
          if(abs(MultB-MultA)>2)cycle
          IrrepAB => ParentIonList(iIonA)%GetIrrep() * ParentIonList(iIonB)%GetIrrep()

          allocate(Pi(nIrreps,nIrreps,nIrreps,3))
          do iIrr1=1,nIrreps
             n1=nactel(iIrr1)
             do iIrr2=1,nIrreps
                n2=nactel(iIrr2)
                do iIrr3=1,nIrreps
                   n3=nactel(iIrr3)
                   iIrr4 = GlobalGroup%GetIrrepIndex( IrrepAB * &
                        GlobalGroup%GetIrrep(iIrr1) * &
                        GlobalGroup%GetIrrep(iIrr2) * &
                        GlobalGroup%GetIrrep(iIrr3) )
                   n4=nactel(iIrr4)
                   do iType = 1, 3
                      call PI(iIrr1,iIrr2,iIrr3,iType)%Init(n1,n2,n3,n4)
                   enddo
                enddo
             enddo
          enddo
          call GetIndexesWithMaxCommonSigma( iIonA, iIonB, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, &
               iStateA, iStateB, Sigma2)
          !.. Note that, in the case of static exchange, all the ions obtained above have Sigma2 = 1
          call SetupSTEXPI(iSymA,NumA,iSymB,NumB,nactel,nRefOccupied,Pi)
          do iIrr1=1,nIrreps
             do iIrr2=1,nIrreps
                do iIrr3=1,nIrreps
                   if(Pi( iIrr1, iIrr2, iIrr3, ALPHA_BETA  )%TotalSize() == 0)cycle
                   call TDM_Manager%SetTDM2( iIonA, iIonB, iIrr1, iIrr2, iIrr3, Sigma2, &
                        Pi( iIrr1, iIrr2, iIrr3, ALPHA_BETA  ) , &
                        Pi( iIrr1, iIrr2, iIrr3, ALPHA_ALPHA ) , &
                        Pi( iIrr1, iIrr2, iIrr3, BETA_BETA   ) )
                   do iType=1,3
                      call Pi( iIrr1, iIrr2, iIrr3, iType  )%Free()
                   enddo
                enddo
             enddo
          enddo
          deallocate(Pi)
       enddo
    enddo
    !
  end subroutine SetupTDM2



  subroutine SetupSTEXRho( iSymA, NumA, iSymB, NumB, nactel, nRefOcc, Irr1Loaded, vnr, vnc, Rho )
    integer        , intent(in) :: iSymA, NumA, iSymB, NumB, nactel(:), nRefOcc(:)
    logical        , intent(out):: Irr1Loaded(:)
    integer        , intent(out):: vnr(:), vnc(:)
    real(kind(1d0)), intent(out):: Rho(:,:,:,:)

    integer         , parameter :: ALPHA = 1
    integer         , parameter :: BETA  = 2
    integer :: iIrrP,  iIrrQ
    integer :: iOrbP,  iOrbQ, iOrbA, iOrbB
    integer :: iA, iB, iP, iQ
    integer :: S2A, S2B, S2P, S2Q

    type(ClassIrrep), pointer    :: IrrepList(:)
    real(kind(1d0))              :: dval
    integer                      :: iType, nIrreps

    Rho = 0.d0
    vnr = nactel
    vnc = nactel
    IrrepList => GlobalGroup%GetIrrepList()
    nIrreps   = size(IrrepList)

    iOrbA = nRefOcc(iSymA) + 1 - NumA
    iOrbB = nRefOcc(iSymB) + 1 - NumB
    S2A   = -1
    S2B   = -1

    iA = AbsIndex( S2A, iSymA, iOrbA, nIrreps )
    iB = AbsIndex( S2B, iSymB, iOrbB, nIrreps )

    !.. Set \Pi^{AB}=0 
    do iIrrP = 1, nIrreps

       iIrrQ = GlobalGroup%GetIrrepIndex( &
            IrrepList(iSymA) * &
            IrrepList(iSymB) * &
            IrrepList(iIrrP) )

       do iType = 1, 2

          select case (iType)
          case(ALPHA)
             S2P =  1
             S2Q =  1
          case(BETA)
             S2P = -1
             S2Q = -1
          end select


          do iOrbP = 1, nactel(iIrrP)
             do iOrbQ = 1, nactel(iIrrQ)

                iP = AbsIndex( S2P, iIrrP, iOrbP, nIrreps )
                iQ = AbsIndex( S2Q, iIrrQ, iOrbQ, nIrreps )

                dval = TDM1STEX(iA,iB,iP,iQ)
                Rho(iOrbP, iOrbQ, iIrrP, iType) = TDM1STEX(iA,iB,iP,iQ)


             enddo
          enddo

       enddo
       Irr1Loaded(iIrrP) = .TRUE.
    enddo


  end subroutine SetupSTEXRho


    subroutine SetupSTEXHab( iSymA, NumA, iSymB, NumB, nactel, nRefOcc, Irr1Loaded, vnr, vnc, Hab )
    integer        , intent(in) :: iSymA, NumA, iSymB, NumB, nactel(:), nRefOcc(:)
    logical        , intent(out):: Irr1Loaded(:)
    integer        , intent(out):: vnr(:), vnc(:)
    real(kind(1d0)), intent(out):: Hab(:,:,:,:)

    integer         , parameter :: ALPHA = 1
    integer         , parameter :: BETA  = 2
    integer :: iIrrP, iIrrQ, iIrrR, iIrrS, iIrrT, iIrrU
    integer :: iOrbA, iOrbB, iOrbP, iOrbQ, iOrbR, iOrbS, iOrbT, iOrbU
    integer :: iA, iB, iP, iQ, iR, iS, iT, iU
    integer :: S2A, S2B, S2P, S2Q, S2S, S2R, S2T, S2U

    type(ClassIrrep), pointer    :: IrrepList(:)
    real(kind(1d0))              :: dval, dval0, dval1
    integer                      :: iType, nIrreps
    type(ClassMatrix), pointer   :: mat1, mat2, mat3
    character(len=16)            :: HAMILTO_Label


    HAMILTO_Label = I1B_ID_LIST(I1B_HAMILTO)

    Hab = 0.d0
    vnr = nactel
    vnc = nactel
    IrrepList => GlobalGroup%GetIrrepList()
    nIrreps   = size(IrrepList)

    iOrbA = nRefOcc(iSymA) + 1 - NumA
    iOrbB = nRefOcc(iSymB) + 1 - NumB
    S2A   = -1
    S2B   = -1

    iA = AbsIndex( S2A, iSymA, iOrbA, nIrreps )
    iB = AbsIndex( S2B, iSymB, iOrbB, nIrreps )


    !*** we need a funtion to order the indexes which werem't stored to avoid redundancies
    !*** we need to run iP and iQ outside the reference configuration except iP .eq. iA and iQ .eq. iB
    !*** we are assuming that the nRefOcc(i) orbitals are located at the begining of the nactel(i) list of orbitals.

    
    !.. Set \Pi^{AB}=0 
    do iIrrP = 1, nIrreps
       iIrrQ = GlobalGroup%GetIrrepIndex( &
            IrrepList(iSymA) * &
            IrrepList(iSymB) * &
            IrrepList(iIrrP) )
       call mat1%Free()
       call mat2%Free()
       call mat3%Free()
       mat1 = GlobalIntegral%Get1B(HAMILTO_Label,"MO_MO",iIrrP,iIrrQ)
       mat2 = GlobalIntegral%Get1B(HAMILTO_Label,"MO_MO",iIrrR,iIrrR)
       mat3 = GlobalIntegral%Get1B(HAMILTO_Label,"MO_MO",iSymA,iSymB)
       do iType = 1, 2

          select case (iType)
          case(ALPHA)
             S2P =  1
             S2Q =  1
          case(BETA)
             S2P = -1
             S2Q = -1
          end select

          do iOrbP = 1, nactel(iIrrP)
             do iOrbQ = 1, nactel(iIrrQ)
                iP = AbsIndex( S2P, iIrrP, iOrbP, nIrreps )
                iQ = AbsIndex( S2Q, iIrrQ, iOrbQ, nIrreps )
                dval = 0.d0
                !.. One body contribution
                if(iA.eq.iB)then
                   !if(.not.GlobalIntegral%MOMOIsInitialized(I1B_HAMILTO,iIrrP,iIrrQ,iOrbP,iOrbQ))cycle
                   if(.not.GlobalIntegral%OneBIsInitialized(HAMILTO_Label,"MO_MO",iIrrP,iIrrQ))cycle
                   !dval = GlobalIntegral%GetMOMO(I1B_HAMILTO,iIrrP,iIrrQ,iOrbP,iOrbQ)
                   dval = mat1%Element(iOrbP,iOrbQ)
                   if(iP.eq.iQ)then
                      dval0 = 0.d0
                      do S2R = 1, -1, -2
                         do iIrrR =1, nIrreps
                            do iOrbR =1, nRefOcc(iIrrR)
                               !if(.not.GlobalIntegral%MOMOIsInitialized(I1B_HAMILTO,iIrrR,iIrrR,iOrbR,iOrbR))cycle
                               !dval = dval + GlobalIntegral%GetMOMO(I1B_HAMILTO,iIrrR,iIrrR,iOrbR,iOrbR)
                               if(.not.GlobalIntegral%OneBIsInitialized(HAMILTO_Label,"MO_MO",iIrrR,iIrrR))cycle
                               dval = dval + mat2%Element(iOrbR,iOrbR)
                            enddo
                         enddo
                      enddo
                   endif
                end if
                if(iP.eq.iQ)then
                   !if(.not.GlobalIntegral%MOMOIsInitialized(I1B_HAMILTO,iSymA,iSymB,iOrbA,iOrbB))cycle
                   !dval = dval - GlobalIntegral%GetMOMO(I1B_HAMILTO,iSymA,iSymB,iOrbA,iOrbB)
                   if(.not.GlobalIntegral%OneBIsInitialized(HAMILTO_Label,"MO_MO",iSymA,iSymB))cycle
                   dval = dval - mat3%Element(iOrbA,iOrbB)
                endif
                
                do S2R = 1, -1, -2
                   do iIrrR =1, nIrreps
                      do iOrbR =1, nRefOcc(iIrrR)
                         do S2S = 1, -1, -2
                            do iIrrS =1, nIrreps
                               do iOrbS =1, nRefOcc(iIrrR)
                                  if(iP.eq.iQ)then
                                     do S2T = 1, -1, -2
                                        do iIrrT =1, nIrreps
                                           do iOrbT =1, nRefOcc(iIrrR)
                                              do S2U = 1, -1, -2
                                                 do iIrrU =1, nIrreps
                                                    do iOrbU =1, nRefOcc(iIrrR)
                                                       if(.not.GlobalIntegral%BielIsInitialized("LL_LL",iIrrR,iIrrS,iIrrT,iIrrU))cycle
                                                       dval0 = TDM2STEX( iA, iB, iR, iT, iS, iU )
                                                       !dval1 = GlobalIntegral%GetLLLL(iIrrR,iIrrS,iIrrT,iIrrU,iOrbR,iOrbS,iOrbT,iOrbU)
                                                       dval1 = GlobalIntegral%GetBiel("LL_LL",iIrrR,iIrrS,iIrrT,iIrrU,iOrbR,iOrbS,iOrbT,iOrbU)
                                                       dval = dval + dval0 * dval1
                                                    enddo
                                                 enddo
                                              enddo
                                           enddo
                                        enddo
                                     enddo
                                  endif
                                  if(.not.GlobalIntegral%BielIsInitialized("LL_LL",iIrrR,iIrrS,iIrrP,iIrrQ))cycle
                                  if(.not.GlobalIntegral%BielIsInitialized("LL_LL",iIrrR,iIrrQ,iIrrP,iIrrS))cycle
                                  dval0 = TDM1STEX(iA,iB,iR,iS)
                                  !dval1 = GlobalIntegral%GetLLLL(iIrrR,iIrrS,iIrrP,iIrrQ,iOrbR,iOrbS,iOrbP,iOrbQ)
                                  dval1 = GlobalIntegral%GetBiel("LL_LL",iIrrR,iIrrS,iIrrP,iIrrQ,iOrbR,iOrbS,iOrbP,iOrbQ)
                                  dval = dval + dval0 * dval1
                                  dval1 = GlobalIntegral%GetBiel("LL_LL",iIrrR,iIrrQ,iIrrP,iIrrS,iOrbR,iOrbQ,iOrbP,iOrbS)
                                  !dval1 = GlobalIntegral%GetLLLL(iIrrR,iIrrQ,iIrrP,iIrrS,iOrbR,iOrbQ,iOrbP,iOrbS)
                                  dval = dval - dval0 * dval1
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
                Hab(iOrbP, iOrbQ, iIrrP, iType) = dval
             enddo
          enddo
       enddo
       Irr1Loaded(iIrrP) = .TRUE.
    enddo

  end subroutine SetupSTEXHab

 




  subroutine SetupSTEXPI( iSymA, NumA, iSymB, NumB, nactel, nRefOcc, Pi )
    use ModuleMatrix
    integer            , intent(in) :: iSymA, NumA, iSymB, NumB
    integer            , intent(in) :: nactel(:)
    integer            , intent(in) :: nRefOcc(:)
    type(ClassMatrix4D), pointer, intent(inout):: Pi(:,:,:,:)

    character(len=*), parameter :: LUCIA_TDM2B="TDM2B_"
    integer         , parameter :: ALPHA_BETA = 1
    integer         , parameter :: ALPHA_ALPHA= 2
    integer         , parameter :: BETA_BETA  = 3

    integer :: iIrrP,  iIrrQ,  iIrrS, iIrrR
    integer :: iOrbP,  iOrbQ,  iOrbS,  iOrbR, iOrbA, iOrbB
    integer :: iA, iB, iP, iQ, iS, iR
    integer :: S2A, S2B, S2P, S2Q, S2S, S2R

    type(ClassIrrep), pointer    :: IrrepList(:)
    real(kind(1d0))              :: dval
    integer                      :: iType, nIrreps

    IrrepList => GlobalGroup%GetIrrepList()
    nIrreps = size(IrrepList)

    iOrbA = nRefOcc(iSymA) + 1 - NumA
    iOrbB = nRefOcc(iSymB) + 1 - NumB
    S2A=-1
    S2B=-1

    iA = AbsIndex( S2A, iSymA, iOrbA, nIrreps )
    iB = AbsIndex( S2B, iSymB, iOrbB, nIrreps )
    !.. Set \Pi^{AB}=0 
    do iIrrP = 1, size(Pi,1)
       do iIrrQ = 1, size(Pi,2)
          do iIrrS = 1, size(Pi,3)

             iIrrR = GlobalGroup%GetIrrepIndex( &
                  IrrepList(iSymA) * &
                  IrrepList(iSymB) * &
                  IrrepList(iIrrP) * &
                  IrrepList(iIrrQ) * &
                  IrrepList(iIrrS) )

             do iType = 1, size(Pi,4)

                select case (iType)
                case(ALPHA_BETA)
                   S2P =  1
                   S2Q = -1
                case(ALPHA_ALPHA)
                   S2P =  1
                   S2Q =  1
                case(BETA_BETA)
                   S2P = -1
                   S2Q = -1
                end select
                S2S = S2P
                S2R = S2R
                ! if(.not.Pi(iIrrP,iIrrQ,iIrrS,iType)%IsInitialized())then
                !    write(*,*)"No es"
                !    stop
                ! endif
                call Pi(iIrrP,iIrrQ,iIrrS,iType)%set(0.d0)

                do iOrbP = 1, nactel(iIrrP)
                   do iOrbQ = 1, nactel(iIrrQ)
                      do iOrbS = 1, nactel(iIrrS)
                         do iOrbR = 1, nactel(iIrrR)
                            iP = AbsIndex( S2P, iIrrP, iOrbP, nIrreps )
                            iQ = AbsIndex( S2Q, iIrrQ, iOrbQ, nIrreps )
                            iS = AbsIndex( S2S, iIrrS, iOrbS, nIrreps )
                            iR = AbsIndex( S2R, iIrrR, iOrbR, nIrreps )

                            dval = TDM2STEX(iA,iB,iP,iQ,iS,iR)
                            call Pi(iIrrP,iIrrQ,iIrrS,iType)%Set(iOrbP,iOrbQ,iOrbS,iOrbR,dval)

                         enddo
                      enddo

                   enddo
                enddo

             enddo
          enddo
       enddo
    enddo

  end subroutine SetupSTEXPI
  !
  integer function AbsIndex( S2, iSym, iOrb, nsym ) result( ind )
    integer,   intent(in) :: S2, iSym, iOrb, nsym
    ind = S2 * ( nsym * ( iOrb - 1 ) + iSym )
  end function AbsIndex
  !
  !>  Computes the 1B transition density matrix between static-exchange ionic states (holes on HF closed-shell det)
  !!   \rho^{BA}_{RS,PQ} =
  !!   \langle \mathbf{K}_{A}| a^{\dagger}_{P}a_{Q} |\mathbf{K}_{B}\rangle=
  !!    (1-\delta_{AP})(1-\delta_{BQ})\times(\delta_{AB}\delta_{PQ}-\delta_{AQ}\delta_{PB})
  !<
  real(kind(1d0)) function TDM1STEX( iA, iB, iP, iQ ) result( dv )
    integer, intent(in) :: iA, iB, iP, iQ
    dv = (1.d0-delta(iA,iP)) * (1.d0-delta(iB,iQ))    *               &
         ( delta(iA,iB) * delta(iP,iQ) - delta(iA,iQ) * delta(iP,iB) )
  end function TDM1STEX
  !
  !>  Computes the 2B transition density matrix between static-exchange ionic states (holes on HF closed-shell det)
  !!   \pi^{BA}_{RS,PQ} =
  !!   \langle \mathbf{K}_{A}| a^{\dagger}_{P}a^{\dagger}_{Q}a_{S}a_{R} |\mathbf{K}_{B}\rangle=
  !!   ( 1 - \delta_{AP} ) ( 1 - \delta_{AQ} ) ( 1 - \delta_{PQ} ) \times
  !!   ( 1 - \delta_{BR} ) ( 1 - \delta_{BS} ) ( 1 - \delta_{RS} ) \times
  !!   \left( \delta_{AB} \delta_{PR} \delta_{QS} - \delta_{AB} \delta_{QR} \delta_{PS} +
  !!          \delta_{QB} \delta_{AR} \delta_{PS} - \delta_{PB} \delta_{AR} \delta_{QS} +
  !!          \delta_{PB} \delta_{QR} \delta_{AS} - \delta_{QB} \delta_{PR} \delta_{AS} \right)
  !<
  real(kind(1d0)) function TDM2STEX( iA, iB, iP, iQ, iS, iR ) result( dv )
    integer, intent(in) :: iA, iB, iP, iQ, iS, iR
    dv = (1.d0-delta(iA,iP)) * (1.d0-delta(iA,iQ)) * (1.d0-delta(iP,iQ)) *   &
         (1.d0-delta(iB,iR)) * (1.d0-delta(iB,iS)) * (1.d0-delta(iR,iS)) * ( &
         delta(iA,iB) * ( delta(iP,iR) * delta(iQ,iS) - delta(iQ,iR) * delta(iP,iS) ) + &
         delta(iA,iR) * ( delta(iQ,iB) * delta(iP,iS) - delta(iP,iQ) * delta(iQ,iS) ) + &
         delta(iA,iS) * ( delta(iP,iB) * delta(iQ,iR) - delta(iQ,iB) * delta(iP,iR) ) )
  end function TDM2STEX
  !
    !
  !>  Computes the 3B transition density matrix between static-exchange ionic states (holes on HF closed-shell det)
  !!   \pi^{BA}_{RS,PQ} =
  !!   \langle \mathbf{K}_{A}| a^{\dagger}_{P}a^{\dagger}_{Q}a_{S}a_{R} |\mathbf{K}_{B}\rangle=
  !!   ( 1 - \delta_{AP} ) ( 1 - \delta_{AQ} ) ( 1 - \delta_{PQ} ) \times
  !!   ( 1 - \delta_{BR} ) ( 1 - \delta_{BS} ) ( 1 - \delta_{RS} ) \times
  !!   \left( \delta_{AB} \delta_{PR} \delta_{QS} - \delta_{AB} \delta_{QR} \delta_{PS} +
  !!          \delta_{QB} \delta_{AR} \delta_{PS} - \delta_{PB} \delta_{AR} \delta_{QS} +
  !!          \delta_{PB} \delta_{QR} \delta_{AS} - \delta_{QB} \delta_{PR} \delta_{AS} \right)
  !<
  real(kind(1d0)) function TDM3STEX( iA, iB, iO, iP, iQ, iT, iS, iR ) result( dv )
    integer, intent(in) :: iA, iB, iO, iP, iQ, iT, iS, iR
    dv = 0.d0
    ! (1.d0-delta(iA,iP)) * (1.d0-delta(iA,iQ)) * (1.d0-delta(iP,iQ)) *   &
    !      (1.d0-delta(iB,iR)) * (1.d0-delta(iB,iS)) * (1.d0-delta(iR,iS)) * ( &
    !      delta(iA,iB) * ( delta(iP,iR) * delta(iQ,iS) - delta(iQ,iR) * delta(iP,iS) ) + &
    !      delta(iA,iR) * ( delta(iQ,iB) * delta(iP,iS) - delta(iP,iQ) * delta(iQ,iS) ) + &
    !      delta(iA,iS) * ( delta(iP,iB) * delta(iQ,iR) - delta(iQ,iB) * delta(iP,iR) ) )
  end function TDM3STEX
  real(kind(1d0)) function delta(i,j) result(res)
    integer, intent(in) :: i, j
    res = 0.d0
    if(i == j) res = 1.d0
  end function delta


  subroutine GetIndexesWithMaxCommonSigma( iIonA, iIonB, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, &
       iStateA, iStateB, Sigma2)
    integer, intent(in)  :: iIonA, iIonB, vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    integer, intent(out) :: iStateA, iStateB, Sigma2

    integer :: iMs2A, iMs2B, Ms2A, Ms2B
    logical :: FOUND

    FOUND=.FALSE.
    Sigma2=-2
    do iMs2A = vPInMs2(iIonA),1,-1
       Ms2A = mPIlistMs2(iMs2A,iIonA)
       do iMs2B = vPInMs2(iIonB),1,-1
          Ms2B = mPIlistMs2(iMs2B,iIonB)
          if(Ms2B==Ms2A)then
             FOUND=.TRUE.
             if(Ms2B>Sigma2)then
                Sigma2=Ms2B
                iStateA=mLUCIAIndexPIMs2(iMs2A,iIonA)
                iStateB=mLUCIAIndexPIMs2(iMs2B,iIonB)
             endif
          endif
       enddo
    enddo
    if(.not.FOUND)then
       write(*,*) "Not found suitable common Sigma to compute Rho1 between ion states",iIonA, iIonB
       stop
    endif
  end subroutine GetIndexesWithMaxCommonSigma


!   subroutine ReadDALTONIrreps( ASTRAIrr_from_DALTONIrr, QCDir )

!     integer, allocatable, intent(out) :: ASTRAIrr_from_DALTONIrr(:)
!     character(len=*)    , intent(in)  :: QCDir
    
!     character(len=*), parameter :: DALTON_OUTPUT = "DALTON.OUT"
!     character(len=*), parameter :: IRREP_LINE_ID = "The irrep name for each symmetry:"

!     integer :: uid, iostat, i, iIrrASTRA, iIrrDALTON, nIrr
!     character(len=1000) :: line, iomsg
!     character(len=3) :: IrrName
!     type(ClassIrrep), pointer :: IrrepPtr

!     open(newunit = uid, &
!          file    = QCDir//"/"//DALTON_OUTPUT, &
!          status  ="old", &
!          form    ="formatted",&
!          action  ="read", &
!          iostat  = iostat, &
!          iomsg   = iomsg )
!     if(iostat/=0)then
!        call ErrorMessage("Error opening "//DALTON_OUTPUT//" "//trim(iomsg))
!        stop
!     endif
!     do
!        read(uid,"(a)",iostat=iostat)line
!        if(iostat/=0)then
!           call ErrorMessage("irrep names not found in "//DALTON_OUTPUT)
!           stop
!        endif
!        i=index(line,IRREP_LINE_ID)
!        if(i>0)exit
!     enddo
!     line=adjustl(line(i+len(IRREP_LINE_ID):))
! !!$    write(*,"(a)")trim(line)

!     nIrr = GlobalGroup%GetNIrreps()

!     if(allocated(ASTRAIrr_from_DALTONIrr))deallocate(ASTRAIrr_from_DALTONIrr)
!     allocate(ASTRAIrr_from_DALTONIrr(nIrr))
!     ASTRAIrr_from_DALTONIrr=0
!     write(*,"(a)") "n Irr DALTON,  n Irr ASTRA,   Irr Name"
!     do iIrrDALTON = 1, nIrr
!        i=index(line,":")
!        line=adjustl(line(i+1:))
!        IrrName = line(1:3)
!        IrrepPtr => GlobalGroup%GetIrrep(trim(line(1:3)))
!        iIrrASTRA= GlobalGroup%GetIrrepIndex(IrrepPtr)
!        write(*,"(2(x,i2),2x,a3)") iIrrDALTON, iIrrASTRA, IrrName
!        ASTRAIrr_from_DALTONIrr(iIrrDALTON)=iIrrASTRA
!     enddo
!     close(uid)


!   end subroutine ReadDALTONIrreps


  subroutine ReadLUCIAIonsAndNactel(ASTRAIrr_from_DALTONIrr, ParentIonList, &
       nactel, nactpairs, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir )

    integer                          , intent(in)  :: ASTRAIrr_from_DALTONIrr(*)
    integer             , allocatable, intent(out) :: nactel(:)
    integer             , allocatable, intent(out) :: nactpairs(:,:)
    type(ClassParentIon), pointer    , intent(in)  :: ParentIonList(:)
    integer             , allocatable, intent(out) :: vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    character(len=*)                 , intent(in)  :: QCDir
    
    integer         , parameter :: MAX_N_MS=3
    character(len=*), parameter :: LUCIA_INFO = "TDM_STATES.info"

    integer :: uid, iostat, iBuf, iIrrASTRA, iIrrDALTON, nIrr, iState, nStates, iState2, iMul
    character(len=1000) :: iomsg

    type(ClassIrrep), pointer :: IrrepPtr
    integer, allocatable :: vMul(:), vIrrDALTON(:), vIon(:), vMs2(:)

    integer :: nPions, iIon, iNum, nMs2, iIrr, iMs2
    logical :: IONS_NOT_FOUND = .FALSE.

    open(newunit = uid, &
         file    = QCDir//"/"//LUCIA_INFO, &
         status  ="old", &
         form    ="formatted",&
         action  ="read", &
         iostat  = iostat, &
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Error opening "//QCDir//"/"//LUCIA_INFO//" "//trim(iomsg))
       stop
    endif

    !.. Read the LUCIA states
    read(uid,*) nStates
    allocate( vMul(nStates), vIrrDALTON(nStates), vMs2(nStates) )
    do iState=1,nStates
       read(uid,*) iBuf, vMul(iState), vIrrDALTON(iState), vMs2(iState)
    enddo

    !.. Determines, for each LUCIA state, the order index of the ion within the symmetry
    !   and multiplicity of that state
    !..
    allocate( vIon(nStates) )
    vIon=0
    do iState=1,nStates
       do iState2=iState-1,1,-1
          if(  vMul      (iState2) == vMul      (iState) .and. &
               vIrrDALTON(iState2) == vIrrDALTON(iState) .and. &
               vMs2      (iState2) == vMs2      (iState) )then
             vIon(iState)=vIon(iState2)
             exit
          endif
       enddo
       vIon(iState)=vIon(iState)+1
       write(*,*) iState, vMul(iState), vIrrDALTON(iState), vIon(iState), vMs2(iState)
    enddo

    !.. Determines, for each ASTRA ion, how many Ms projections are available, from LUCIA,
    !   as well as the list of these Ms2, and the corresponding absolute state LUCIA index
    !..
    nPions = size(ParentIonList)
    if(allocated(vPInMs2)) deallocate(vPInMs2)
    allocate(vPInMs2(nPions))
    if(allocated(mPIlistMs2)) deallocate(mPIlistMs2)
    allocate(mPIlistMs2(MAX_N_MS,nPions))
    if(allocated(mLUCIAIndexPIMs2)) deallocate(mLUCIAIndexPIMs2)
    allocate(mLUCIAIndexPIMs2(MAX_N_MS,nPions))
    vPInMs2=0
    mPIlistMs2=0
    mLUCIAIndexPIMs2=0
    do iIon = 1, nPions

       iMul     =  ParentIonList(iIon)%GetMultiplicity()
       IrrepPtr => ParentIonList(iIon)%GetIrrep()
       iIrr     =  GlobalGroup%GetIrrepIndex(IrrepPtr)
       iNum     =  ParentIonList(iIon)%GetNumber()

       nMs2 = 0
       do iState = 1, nStates
          iIrrASTRA = ASTRAIrr_from_DALTONIrr( vIrrDALTON(iState) )
          if(  iMul == vMul(iState) .and. &
               iIrr == iIrrASTRA    .and. &
               iNum == vIon(iState) )then
             nMs2 = nMs2 + 1
             vPInMs2(iIon)=nMs2
             mPIlistMs2(nMs2,iIon)= vMs2(iState)
             mLUCIAIndexPIMs2(nMs2,iIon) = iState
          endif
       enddo
       if(nMs2==0)then
          IONS_NOT_FOUND = .TRUE.
          write(*,*) "Could not find in LUCIA the following ASTRA ion "
          call ParentIonList(iIon)%Show()
       endif

!!$       write(*,*) "###", iIon, iMul, iIrr, iNum, vPInMs2(iIon)
!!$       do iMs2 = 1, vPInMs2(iIon)
!!$          write(*,*) "     ",mPIlistMs2(iMs2, iIon), mLUCIAIndexPIMs2(iMs2, iIon)
!!$       enddo

    enddo
    if(IONS_NOT_FOUND) stop


!!$       

    !.. Read the number of string of operators of a given kind for each symmetry
    !   1  a_{\alpha\Gamma n}    !<--- destructor operator spin alpha, sym gamma, orbital n
    !   2  a_{\beta \Gamma n}    
    !   3  a_{\alpha \Gamma' n} a_{\beta \Gamma" m}  ! where sym gamma = gamma' x gamma"
    !   4  a_{\alpha \Gamma' n} a_{\alpha\Gamma" m}  
    !   5  a_{\beta  \Gamma' n} a_{\beta \Gamma" m}  
    !..
    nIrr = GlobalGroup%GetNIrreps()
    if(allocated(nactel))deallocate(nactel)
    allocate(nactel(nIrr))
    if(allocated(nactpairs))deallocate(nactpairs)
    allocate(nactpairs(nIrr,2))
    read(uid,*)
    do iIrrDALTON = 1, nIrr
       iIrrASTRA = ASTRAIrr_from_DALTONIrr( iIrrDALTON )
       read(uid,*) iBuf, nactel(iIrrASTRA), iBuf, &
            nactpairs(iIrrASTRA,1), nactpairs(iIrrASTRA,2), iBuf
    enddo
    close(uid)

  end subroutine ReadLUCIAIonsAndNactel




end program ConvertDensityMatricesSTEX
