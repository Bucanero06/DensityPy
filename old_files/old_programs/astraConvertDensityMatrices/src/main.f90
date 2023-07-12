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
!>  Insert Doxygen comments here
!!  
!!  This program read, write and / or generates the density matrices and does
!!  the several necessary transformations (spin coupling etc).
!! 
program ParseDensityMatrices

  use, intrinsic :: ISO_FORTRAN_ENV

  !.. Level 0.0
  use ModuleErrorHandling
  use ModuleString
  use ModuleIO
  !.. Level 0.1
  use ModuleGroups
  !.. Level 1
  use ModuleAstraConfigFile
  use ModuleOrbitalBasis
  !.. Level 2
  use ModuleParentIons
  !.. Level 3
  use ModuleDensityMatrices
  !.. Level 5
  use ModuleESpace

  use ModuleMainInterface

  implicit none

  !.. Run time parameters
  character(len=:), allocatable :: AstraConfigFile
  logical                       :: CheckConversion

  !.. Config file parameters
  character(len=:), allocatable :: StorageDir
  character(len=:), allocatable :: ccConfigFile
  character(len=:), allocatable :: QCDir

  type(ClassESpace)             :: Space
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

  call GetRunTimeParameters( AstraConfigFile, CheckConversion )

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
       ninactive, nactel, nactpairs, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir )

  !.. Pass the list of parent ions to the Transition-Density-Matrix Manager
  call TDM_Manager%SetParentIons( ParentIonList )
  call TDM_Manager%SetNActiveEls( ninactive, nactel )
  call TDM_Manager%SetStorageDir( StorageDir )
  call TDM_Manager%Init()

  call SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, "1", QCDir )
  call SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, "H", QCDir )
  call SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, "E", QCDir )
  call SetupTDM2( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir )

  !.. Save the TDM in the format of astra
  call TDM_Manager%IO( "1", "Write" )
  call TDM_Manager%IO( "2", "Write" )
  call TDM_Manager%IO( "H", "Write" )
  call TDM_Manager%IO( "E", "Write" )
  call TDM_Manager%Free()

  stop

contains


  subroutine SetupTDM1( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, &
       ID_CHAR, QCDir )

    type(ClassParentIon), intent(in) :: ParentIonList(:)
    integer             , intent(in) :: ASTRAIrr_from_DALTONIrr(:)
    integer             , intent(in) :: nactel(:), vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    character           , intent(in) :: ID_CHAR
    character(len=*)    , intent(in) :: QCDir

    integer                          :: iIonA, iIonB, MultA, MultB, iStateA, iStateB, Sigma2, nIons
    integer                          :: numA, numB, iSymA, iSymB
    integer                          :: nIrreps, nOrbMax, iIrr1, nr, nc
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

          call LoadRho(iStateA,iStateB,ASTRAIrr_from_DALTONIrr,Irr1Loaded,vnr,vnc,Rho,ID_CHAR,QCDir)

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


  !.. This subroutine performs two tasks:
  !   1. Loads the two-body transition density matrix between two ionic states in the 
  !   format provided by LUCIA, namely
  !   \f[
  !      \pi_(ij,kl;type)^{AMs,BMs} = < A, Ms | C_{ij;type} A_{kl;type} | B, Ms >
  !   \f]
  !   where C_{ij,type} and A_{kl;type} are pairs of constructor and destructor operators,
  !   respectively, 
  !   \f[
  !      A_{kl;type} = a_{k,Ms_{1,type}}a_{l,Ms_{2,type}} \\
  !      B_{kl;type} = a^\dagger_{k,Ms_{1,type}}a^\dagger_{l,Ms_{2,type}}.
  !   \f]
  !   In the case of two-body TDM, in Lucia terminology, there are three possible types:
  !   type = 3    =>  Ms_1, Ms_2   = \alpha, \beta
  !   type = 4    =>  Ms_1, Ms_2   = \alpha, \alpha
  !   type = 5    =>  Ms_1, Ms_2   = \beta , \beta
  !   which include all possible cases, since the 2BTDM is antisymmetric with respect
  !   to exchanges of either the two creators or the two constructors (type 1 and 2 concern
  !   one-body TDMs). Since the A and B ionic states have always the same spin projection, 
  !   the A and C operators are both of the same type.
  !   Thanks to permutation symmetry, for operators of type=4 and 5, LUCIA only provides
  !   orbital symmetries in incresing order within the pair in the A or C string. 
  !   In case both orbitals in an A or C string have the same symmetry, the orbital indexes 
  !   within that symmetry are also in ascending order.
  !   Finally, LUCIA lists the matrix elements of the TDM2B for each kind (ab,aa,bb) in 
  !   rectangular blocks identified by the symmetry of the creator string. The symmetry of
  !   the destructor string is uniquely identified from that of the creator and of the two
  !   ions'. The indexes of such a block run over all the possible symmetry pairs of the two
  !   creator operators. For type 3 (alpha,beta), for each symmetry of the alpha operator, cycles
  !   over the symmetry of the beta operator. For type 4 and 5 ***
  !
  !   ======== EXPLANATION FROM LUCIA note_order_genop_strings ===============================
  !   The individual string-part (ca, cb, aa, ab) for a given CAAB-type
  !   are organized so these are grouped according to symmetry.
  !   Furthermore, only the unique strings are included. This is ensured by requiring that the
  !   orbitals in a string are in strictly ascending order, when reading from left to right.
  !   Thus for the the strings
  !   a_{k\alpha} a_{l\alpha} the ordering k<l is used.
  !
  !   The symmetry ordering of a CAAB group strings in terms of its four components
  !   goes as - for strings with total symmetry ISM
  !
  !   Loop over symmetry over creator strings-parts (ISM_C)
  !   Symmetry ISM_A of annihilation strings is then obtained from the total symmetry(ISM) and ISM_C
  !     Loop over symmetry of creation alpha-stringpart, ISM_CA
  !     Symmetry of creation beta-stringparts can then be deduced from ISM_CA and ISM_C
  !      Loop over symmetry of annihilation alpha-stringparts, ISM_AA
  !      Symmetry of annihilation beta-stringparts can then be deduced from ISM_AA and ISM_A
  !    . We have now symmetry for the four stringparts of the total string defined, so
  !    . now we start to loop over the string-parts with given symmetries
  !    . The ordering corresponds to a matrix T(I_CA, I_CB, IAA, I_AB) in Fortran, so the
  !    . loop structure is
  !        Loop over ab stringparts with given symmetry
  !         Note: The ordering within a stringpart is described below
  !         loop over aa stringparts with given symmetry
  !          loop over cb stringparts of given symmetry
  !           loop over ca stringparts of given symmetry
  !
  !           end of loop over ca stringparts of given symmetry
  !          end of loop over cb stringparts of given symmetry
  !         end of loop over aa stringparts with given symmetry
  !        end of loop over ab stringparts with given symmetry
  !    
  !      end of Loop over symmetry of annihilation alpha-stringparts
  !     end of Loop over symmetry of creation alpha-strinparts
  !   end of Loop over symmetry over creator strinparts (ISM_C)
  !
  !   We now have to define the ordering of strings of a given components, 
  !   say the creation of alpha-electrons.
  !
  !   The general ordering is a bit hairy, so let us confine us for this time 
  !   around to string-parts containing 1 or 2 electrons.  Furthermore, we will 
  !   restrict us to discuss the case with a active orbital space, i.e. a CAS calculation.
  !
  !   For 1 electron, the scheme is simple, just loop over the orbitals
  !   of the given symmetry in ascending order.
  !
  !   For 2 electrons with the same spin, the string-part has the form
  !   a^\dagger_I a^\dagger_J, where I and J are spin-orbitals. If the total symmetry
  !   of the 2-electron string should be SYM_IJ, the loops goes as
  !
  !   Loop over symmetry SYM_I of I
  !     SYM_J = SYM_I x SYM_IJ.
  !     if( SYM_J < SYM_I ) cycle  (The loop is restricted so SYM_J geq SYM_I)
  !     !
  !     Loop over orbitals J of symmetry SYM_J
  !       Loop over orbitals I of symmetry SYM_I
  !       The loop is restricted so I<J if SYM_I = SYM_J
  !         !
  !         Next string
  !         !
  !       End of loop over orbitals I of symmetry SYM_I
  !     End of loop over orbitals J of symmetry SYM_J
  !     !
  !   End of loop over symmetry SYM_I
  !
  !   ==================================================================================
  !
  !   
  !   2. Store the TDM in the ASTRA format, via the methods of the TDMManager methods.
  !..
  subroutine SetupTDM2( ParentIonList, nactel, &
       ASTRAIrr_from_DALTONIrr, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir )
    use ModuleMatrix

    type(ClassParentIon), intent(in) :: ParentIonList(:)
    integer             , intent(in) :: ASTRAIrr_from_DALTONIrr(:)
    integer             , intent(in) :: nactel(:) !, nactpairs(:,:)
    integer             , intent(in) :: vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    character(len=*)    , intent(in) :: QCDir

    integer, parameter :: ALPHA_BETA = 1
    integer, parameter :: ALPHA_ALPHA= 2
    integer, parameter :: BETA_BETA  = 3
    !
    integer :: iIonA, iIonB, MultA, MultB, nIons
    integer :: numA, numB, iSymA, iSymB
    integer :: nIrreps, iIrr1, iIrr2, iIrr3, iIrr4, n1, n2, n3, n4, iType
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

          call LoadPI( iStateA, iStateB, nactel, ASTRAIrr_from_DALTONIrr, Pi, QCDir )

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


  subroutine LoadRho(iStateA,iStateB,ASTRAIrr_from_DALTONIrr,Irr1Loaded,vnr,vnc,Rho,ID_CHAR,QCDir)
    integer         , intent(in) :: iStateA, iStateB, ASTRAIrr_from_DALTONIrr(:)
    logical         , intent(out):: Irr1Loaded(:)
    integer         , intent(out):: vnr(:), vnc(:)
    real(kind(1d0)) , intent(out):: Rho(:,:,:,:)
    character       , intent(in) :: ID_CHAR
    character(len=*), intent(in) :: QCDir

    character(len=*), parameter   :: LUCIA_TDM1B="LUCIA_BLK+1-_"
    character(len=*), parameter   :: LUCIA_HAIC ="LUCIA_BLK-H+_"
    character(len=*), parameter   :: LUCIA_PIE  ="LUCIA_BLKH_"
    character(len=:), allocatable :: sLUCIA
    character(len=20) :: Suffix
    character(len=500):: iomsg
    integer :: uid, iostat, iIrrASTRA, iIrrDALTON, iType, iBuf, nr, nc, i, j

    if(ID_CHAR.is."1")then
       sLUCIA=LUCIA_TDM1B
    elseif(ID_CHAR.is."H")then
       sLUCIA=LUCIA_HAIC
    elseif(ID_CHAR.is."E")then
       sLUCIA=LUCIA_PIE
       if(iStateA/=iStateB)return
    endif

    Rho        = 0.d0
    vnr        = 0
    vnc        = 0
    Irr1Loaded = .FALSE.

    Suffix = adjustl(AlphabeticNumber(iStateA)//"."//AlphabeticNumber(iStateB))
    open(newunit = uid, &
         file    = QCDir//"/"//sLUCIA//trim(Suffix), &
         status  ="old", &
         form    ="formatted",&
         action  ="read", &
         iostat  = iostat, &
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Error opening "//QCDir//"/"//sLUCIA//trim(Suffix)//" "//trim(iomsg))
       stop
    endif
    do
       read(uid,*,iostat=iostat) iType, iIrrDALTON, iBuf, nr, nc
       if(iostat/=0)exit
       if(iType/=1.and.iType/=2)then
          write(*,*) "Inconsistent type in "//sLUCIA//trim(Suffix)
          stop
       endif
       iIrrASTRA=ASTRAIrr_from_DALTONIrr(iIrrDALTON)
       Irr1Loaded(iIrrASTRA)=.TRUE.
       vnr(iIrrASTRA)=nr
       vnc(iIrrASTRA)=nc
       do i=1,nr
          read(uid,*) ( Rho(i,j,iIrrASTRa,iType), j = 1, nc )
       enddo
    enddo
    close(uid)

  end subroutine LoadRho

  subroutine LoadPi( iStateA, iStateB, nactel, ASTRAIrr_from_DALTONIrr, Pi, QCDir )
    use ModuleMatrix
    integer            , intent(in) :: iStateA, iStateB
    integer            , intent(in) :: nactel(:)
    integer            , intent(in) :: ASTRAIrr_from_DALTONIrr(:)
    type(ClassMatrix4D), intent(inout):: Pi(:,:,:,:)
    character(len=*)   , intent(in) :: QCDir

    character(len=*), parameter :: LUCIA_TDM2B="LUCIA_BLK++1--_"
    integer         , parameter :: ALPHA_BETA = 1
    integer         , parameter :: ALPHA_ALPHA= 2
    integer         , parameter :: BETA_BETA  = 3

    character(len=20) :: Suffix
    character(len=500):: iomsg
    integer :: iIrr,   iIrr1,  iIrr2,  iIrr3
    integer :: iIrrAa, iIrrAb, iIrrA1, iIrrA2
    integer :: iIrrCa, iIrrCb, iIrrC1, iIrrC2
    integer :: iOrb1,  iOrb2,  iOrb3,  iOrb4
    integer :: iIrrDALTON_A, iIrrDALTON_Aa, iIrrDALTON_A1, iIrrDALTON_A2
    integer :: iIrrDALTON_C, iIrrDALTON_Ca, iIrrDALTON_C1, iIrrDALTON_C2
    integer :: iIrrASTRA_A,  iIrrASTRA_C
    integer, allocatable :: DALTONIrr_from_ASTRAIrr(:)

    type(ClassIrrep), pointer    :: IrrepList(:)
    type(ClassIrrep), pointer    :: Irr_A, Irr_Aa, Irr_Ab, Irr_A1, Irr_A2
    type(ClassIrrep), pointer    :: Irr_C, Irr_Ca, Irr_Cb, Irr_C1, Irr_C2

    real(kind(1d0)), allocatable :: dMat(:,:)
    integer                      :: i, j, nr, nc
    integer                      :: uid, iostat, iType, iTypeLUCIA, nIrreps

    IrrepList => GlobalGroup%GetIrrepList()
    nIrreps = size(IrrepList)

    allocate(DALTONIrr_from_AstraIrr(nIrreps))
    do iIrr = 1, nIrreps
       DALTONIrr_from_AstraIrr( ASTRAIrr_from_DALTONIrr( iIrr ) ) = iIrr
    enddo

    !.. Set \Pi^{AB}=0 
    do iIrr1 = 1, size(Pi,1)
       do iIrr2 = 1, size(Pi,2)
          do iIrr3 = 1, size(Pi,3)
             do iType = 1, size(Pi,4)
                call Pi(iIrr1,iIrr2,iIrr3,iType)%set(0.d0)
             enddo
          enddo
       enddo
    enddo

    Suffix = adjustl(AlphabeticNumber(iStateA)//"."//AlphabeticNumber(iStateB))
    open(newunit = uid, &
         file    = QCDir//"/"//LUCIA_TDM2B//trim(Suffix), &
         status  ="old", &
         form    ="formatted",&
         action  ="read", &
         iostat  = iostat, &
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Error opening "//QCDir//"/"//LUCIA_TDM2B//trim(Suffix)//" "//trim(iomsg))
       stop
    endif
    do

       read(uid,*,iostat=iostat) iTypeLUCIA, iIrrDALTON_C, iIrrDALTON_A, nr, nc
       if(iostat/=0)exit
       if(iTypeLUCIA<3)then
          write(*,*) "Inconsistent type in "//QCDir//"/"//LUCIA_TDM2B//trim(Suffix)
          stop
       endif
       if(iTypeLUCIA>5)exit
       if(allocated(dMat))deallocate(dMat)
       allocate(dMat(nr,nc))
       do i=1,nr
          read(uid,*) ( dMat(i,j), j = 1, nc )
       enddo

       !.. Determine the index and pointer to the irrep of the creation 
       !   operator string in ASTRA ordering
       iIrrASTRA_C =  ASTRAIrr_from_DALTONIrr(iIrrDALTON_C)
       Irr_C       => IrrepList( iIrrASTRA_C )

       !.. Determine the index and pointer to the irrep of the annihilation 
       !   operator string in ASTRA ordering
       iIrrASTRA_A =  ASTRAIrr_from_DALTONIrr(iIrrDALTON_A)
       Irr_A       => IrrepList( iIrrASTRA_A )

       select case(iTypeLUCIA)
       case(3)
          iType = ALPHA_BETA
       case(4)
          iType = ALPHA_ALPHA
       case(5)
          iType = BETA_BETA
       end select

       if(iType == ALPHA_BETA )then 

          i=0
          do iIrrDALTON_Ca = 1, nIrreps

             iIrrCa =  ASTRAIrr_from_DALTONIrr( iIrrDALTON_Ca )
             Irr_Ca => IrrepList( iIrrCa )
             Irr_Cb => Irr_C * Irr_Ca
             iIrrCb =  GlobalGroup%GetIrrepIndex( Irr_Cb )

             if(nactel(iIrrCa)*nactel(iIrrCb)==0)cycle

             do iOrb2 = 1, nactel(iIrrCb)
                do iOrb1 = 1, nactel(iIrrCa)
                   i=i+1

                   j = 0
                   do iIrrDALTON_Aa = 1, nIrreps

                      iIrrAa =  ASTRAIrr_from_DALTONIrr(iIrrDALTON_Aa)
                      Irr_Aa => IrrepList( iIrrAa )
                      Irr_Ab => Irr_A * Irr_Aa
                      iIrrAb =  GlobalGroup%GetIrrepIndex( Irr_Ab )

                      if(nactel(iIrrAa)*nactel(iIrrAb)==0)cycle

                      do iOrb4 = 1, nactel(iIrrAb)
                         do iOrb3 = 1, nactel(iIrrAa)
                            j = j+1

                            call Pi(iIrrCa,iIrrCb,iIrrAa,iType)%Set(iOrb1,iOrb2,iOrb3,iOrb4,dMat(i,j))

                         enddo
                      enddo

                   enddo

                enddo
             enddo
          enddo

       else !.. alpha-alpha or beta-beta

          i=0
          do iIrrDALTON_C2 = 1, nIrreps

             iIrrC2 =  ASTRAIrr_from_DALTONIrr(iIrrDALTON_C2)
             Irr_C2 => IrrepList( iIrrC2 )
             Irr_C1 => Irr_C * Irr_C2
             iIrrC1 =  GlobalGroup%GetIrrepIndex( Irr_C1 )
             iIrrDALTON_C1 = DALTONIrr_from_ASTRAIrr( iIrrC1 )

             if( iIrrDALTON_C1 < iIrrDALTON_C2    ) cycle
             if(nactel(iIrrC1)*nactel(iIrrC2)==0) cycle

             do iOrb2 = 1, nactel(iIrrC2)
                do iOrb1 = 1, nactel(iIrrC1)
                   if( iIrrC1 == iIrrC2 .and. iOrb1 <= iOrb2 ) cycle
                   i=i+1

                   j = 0
                   do iIrrDALTON_A2 = 1, nIrreps

                      iIrrA2 =  ASTRAIrr_from_DALTONIrr(iIrrDALTON_A2)
                      Irr_A2 => IrrepList( iIrrA2 )
                      Irr_A1 => Irr_A * Irr_A2
                      iIrrA1 =  GlobalGroup%GetIrrepIndex( Irr_A1 )
                      iIrrDALTON_A1 = DALTONIrr_from_ASTRAIrr( iIrrA1 )

                      if( iIrrDALTON_A1 < iIrrDALTON_A2    ) cycle
                      if(nactel(iIrrA1)*nactel(iIrrA2)==0) cycle

                      do iOrb4 = 1, nactel(iIrrA2)
                         do iOrb3 = 1, nactel(iIrrA1)
                            if( iIrrA1 == iIrrA2 .and. iOrb3 <= iOrb4 ) cycle
                            j = j+1

                            call Pi(iIrrC1,iIrrC2,iIrrA1,iType)%Set(iOrb1,iOrb2,iOrb3,iOrb4,  dMat(i,j) )
                            call Pi(iIrrC1,iIrrC2,iIrrA2,iType)%Set(iOrb1,iOrb2,iOrb4,iOrb3, -dMat(i,j) )
                            call Pi(iIrrC2,iIrrC1,iIrrA1,iType)%Set(iOrb2,iOrb1,iOrb3,iOrb4, -dMat(i,j) )
                            call Pi(iIrrC2,iIrrC1,iIrrA2,iType)%Set(iOrb2,iOrb1,iOrb4,iOrb3,  dMat(i,j) )

                         enddo
                      enddo

                   enddo

                enddo
             enddo
          enddo

       endif

    enddo
    close(uid)

  end subroutine LoadPi


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




  subroutine ReadLUCIAIonsAndNactel(ASTRAIrr_from_DALTONIrr, ParentIonList, &
       ninactive, nactel, nactpairs, vPInMs2, mPIlistMs2, mLUCIAIndexPIMs2, QCDir )

    integer                          , intent(in)  :: ASTRAIrr_from_DALTONIrr(*)
    integer             , allocatable, intent(out) :: ninactive(:)
    integer             , allocatable, intent(out) :: nactel(:)
    integer             , allocatable, intent(out) :: nactpairs(:,:)
    type(ClassParentIon), pointer    , intent(in)  :: ParentIonList(:)
    integer             , allocatable, intent(out) :: vPInMs2(:), mPIlistMs2(:,:), mLUCIAIndexPIMs2(:,:)
    character(len=*)                 , intent(in)  :: QCDir

    integer         , parameter :: MAX_N_MS=3
    character(len=*), parameter :: LUCIA_INFO = "LUCIA_STATES_2.info"
    character(len=*), parameter :: LUCIA_CFG  = "LUCIA.INP"

    integer :: uid, iostat, iBuf, iIrrASTRA, iIrrDALTON, nIrr, iState, nStates, iState2, iMul
    character(len=1000) :: iomsg

    type(ClassIrrep), pointer :: IrrepPtr
    integer, allocatable :: vMul(:), vIrrDALTON(:), vIon(:), vMs2(:)

    integer :: nPions, iIon, iNum, nMs2, iIrr
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

    enddo
    if(IONS_NOT_FOUND) stop

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

    BLOCK
      character(len=1000) :: line
      integer :: i
      !.. Read the inactive orbitals
      open(newunit = uid, &
           file    = QCDir//"/"//LUCIA_CFG, &
           status  ="old", &
           form    ="formatted",&
           action  ="read", &
           iostat  = iostat, &
           iomsg   = iomsg )
      if(iostat/=0)then
         call ErrorMessage("Error opening "//QCDir//"/"//LUCIA_CFG//" "//trim(iomsg))
         stop
      endif
      do
         call FetchLineStrn(uid,line,"*")
         if(len_trim(line)==0)exit
         i=index(line,"Inash")
         if(i>0)exit
      enddo
      allocate(ninactive(nIrr))
      ninactive=0
      if(i>0)then
         call FetchLineStrn(uid,line,"*")
         if(len_trim(line)==0)then
            write(ERROR_UNIT,"(a)") "Unexpected EOF after Inash in "//QCDir//"/"//LUCIA_CFG
            stop
         endif
      endif
      do
         i=index(line,",")
         if(i<=0)exit
         line(i:i)=" "
      enddo
      read(line,*) (ninactive(ASTRAIrr_from_DALTONIrr(iIrrDALTON)),iIrrDALTON=1,nirr)
      close(uid)
    END BLOCK

  end subroutine ReadLUCIAIonsAndNactel


end program ParseDensityMatrices