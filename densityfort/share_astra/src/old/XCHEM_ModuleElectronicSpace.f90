module ModuleElectronicSpace
  
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleIO
  use ModuleString
  use ModuleGroups
  use ModuleParentIons
  use ModuleSymmetricElectronicSpace
  use ModuleMatrix
  use ModuleShortRangeOrbitals
  
  implicit none

  private
  
  character(len=*), parameter :: CLOSE_COUPLING_DIR      = "CloseCoupling"
  character(len=*), parameter :: HamiltonianFileRootName = "H_"
  character(len=*), parameter :: KinEnergyFileRootName   = "K_"
  character(len=*), parameter :: OverlapFileRootName     = "S_"
  character(len=*), parameter :: DipoleLenFileRootName   = "DipoleLen_"
  character(len=*), parameter :: DipoleVelFileRootName   = "DipoleVel_"
  character(len=*), parameter :: QCFileExtension         = "QC"


  type, public :: ClassElectronicSpace
     
     private

     logical          :: MULTIPLICITY_IS_SET = .FALSE.

     integer          :: Multiplicity
     type(ClassGroup) :: Group
     integer          :: ParentIonCharge
     integer          :: Lmax
     integer          :: Ne

     integer          :: NIrreps

     integer                       :: NParentIons
     type(ClassParentIon), pointer :: ParentIonList(:)

     !.. Number of orbitals, per each irrep, that is doubly occupied
     !   in a closed-shell reference configuration. This parameter
     !   is useful only in the case of the static-exchange approximation
     !   where such reference configuration is needed to define the
     !   nature of all parent ions
     !..
     integer, allocatable :: nRefOccupied(:)

     type(ClassSymmetricElectronicSpace), allocatable :: IrrepSpaceVec(:)
     
     !.. This is where the close-coupling data for an assigned QCI are stored
     !   and it includes the nuclear configuration, if present
     !..
     character(len=:), allocatable :: StorageDir
     !> Label of Parametric configuration (e.g., nuclear geometry)
     character(len=:), allocatable :: NuclearLabel

   contains
     
     generic, public :: Initialized        => ClassElectronicSpaceInitialized
     generic, public :: ParseConfigFile    => ClassElectronicSpaceParseConfigFile
     generic, public :: SetMultiplicity    => ClassElectronicSpaceSetMultiplicity
     generic, public :: SetStorageDir      => ClassElectronicSpaceSetStorageDir
     generic, public :: SetNuclearLabel    => ClassElectronicSpaceSetNuclearLabel
     
     generic, public :: GetGroup           => ClassElectronicSpaceGetGroup
     generic, public :: GetGroupName       => ClassElectronicSpaceGetGroupName
     generic, public :: GetPionList        => ClassElectronicSpaceGetPionList
     !> Gets the maximum L defined in the Close-Coupling anzat
     !! for the current ClassElectronicSpace%
     generic, public :: GetMaxLinCC        => ClassElectronicSpaceGetMaxLinCC
     generic, public :: GetLmax            => ClassElectronicSpaceGetLmax
     generic, public :: GetStorageDir      => ClassElectronicSpaceGetStorageDir
     generic, public :: GetPionCharge      => ClassElectronicSpaceGetPionCharge
     generic, public :: GetNIrreps         => ClassElectronicSpaceGetNIrreps
     generic, public :: GetnRefOccupied    => ClassElectronicSpaceGetnRefOccupied
     generic, public :: CheckSymmetry      => ClassElectronicSpaceCheckSymmetry
     generic, public :: GetSymElectSpace   => ClassElectronicSpaceGetSymElectSpace, ClassElectronicSpaceGetSymElectSpaceByIrrepName

     generic, public :: AcquireOverlap     => ClassElectronicSpaceAcquireOverlap
     generic, public :: AcquireHamiltonian => ClassElectronicSpaceAcquireHamiltonian
     generic, public :: AcquireKinEnergy   => ClassElectronicSpaceAcquireKinEnergy
     generic, public :: AcquireDipoleLen   => ClassElectronicSpaceAcquireDipoleLen
     generic, public :: AcquireDipoleVel   => ClassElectronicSpaceAcquireDipoleVel
     generic, public :: Show                => ClassElectronicSpaceShow


     ! {{{ private procedures

     procedure, private :: ClassElectronicSpaceInitialized
     procedure, private :: ClassElectronicSpaceParseConfigFile
     procedure, private :: ClassElectronicSpaceSetMultiplicity
     procedure, private :: ClassElectronicSpaceSetStorageDir
     procedure, private :: ClassElectronicSpaceSetNuclearLabel

     procedure, private :: ClassElectronicSpaceGetPionList
     procedure, private :: ClassElectronicSpaceGetStorageDir
     procedure, private :: ClassElectronicSpaceGetGroup
     procedure, private :: ClassElectronicSpaceGetGroupName
     procedure, private :: ClassElectronicSpaceGetLmax
     procedure, private :: ClassElectronicSpaceGetMaxLinCC
     procedure, private :: ClassElectronicSpaceGetPionCharge
     procedure, private :: ClassElectronicSpaceGetNIrreps
     procedure, private :: ClassElectronicSpaceGetnRefOccupied
     procedure, private :: ClassElectronicSpaceCheckSymmetry
     procedure, private :: ClassElectronicSpaceGetSymElectSpace
     procedure, private :: ClassElectronicSpaceGetSymElectSpaceByIrrepName

     procedure, private :: ClassElectronicSpaceAcquireOverlap
     procedure, private :: ClassElectronicSpaceAcquireHamiltonian
     procedure, private :: ClassElectronicSpaceAcquireKinEnergy
     procedure, private :: ClassElectronicSpaceAcquireDipoleLen
     procedure, private :: ClassElectronicSpaceAcquireDipoleVel
     procedure, private :: ClassElectronicSpaceShow

     ! }}}
  end type ClassElectronicSpace


  public :: GetCloseCouplingDir
  public :: GetQCFileExtension

  

contains


  subroutine ClassElectronicSpaceAcquireOverlap( Space, SymLabel, PWCBraLabel, PWCKetLabel, inFileName )
    class(ClassElectronicSpace), intent(in) :: Space
    character(len=*)           , intent(in) :: SymLabel
    character(len=*)           , intent(in) :: PWCBraLabel
    character(len=*)           , intent(in) :: PWCKetLabel
    character(len=*)           , intent(in) :: inFileName
    character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
    logical :: exist
    type(ClassMatrix) :: MatA
    integer :: uid

    call ErrorMessage("Must check validity of "//trim(SymLabel)//", " //trim(PWCBraLabel)//", and " //trim(PWCKetLabel))

    allocate(outDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
         trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
         trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

    call Execute_Command_Line(" mkdir -p "//OutDir)

    allocate(outFileName,source=outDir//OverlapFileRootName//"0.0"//QCFileExtension)

    INQUIRE( file = inFileName, exist = exist )
    if ( exist ) then
       call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
       call CheckWrittenBlock( outFileName )
    else
       allocate(auxDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
            trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
            trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
       allocate(auxFileName,source=auxDir//OverlapFileRootName//"0.0"//QCFileExtension)
       call OpenFile( auxFileName, uid, 'read', 'formatted' )
       call MatA%Read( uid )
       close( uid )
       call CheckBlock( MatA )
       call MatA%Transpose()
       call OpenFile( outFileName, uid, 'write', 'formatted' )
       call MatA%Write( uid )
       close( uid )
    end if

  end subroutine ClassElectronicSpaceAcquireOverlap


  subroutine ClassElectronicSpaceAcquireHamiltonian( Space, SymLabel, PWCBraLabel, PWCKetLabel, inFileName )
    class(ClassElectronicSpace), intent(in) :: Space
    character(len=*)           , intent(in) :: SymLabel
    character(len=*)           , intent(in) :: PWCBraLabel
    character(len=*)           , intent(in) :: PWCKetLabel
    character(len=*)           , intent(in) :: inFileName
    character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
    logical :: exist
    type(ClassMatrix) :: MatA
    integer :: uid

    call ErrorMessage("Must check validity of "//trim(SymLabel)//", " //trim(PWCBraLabel)//", and " //trim(PWCKetLabel))

    allocate(outDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
         trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
         trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

    call Execute_Command_Line(" mkdir -p "//OutDir)

    allocate(outFileName,source=outDir//HamiltonianFileRootName//"0.0"//QCFileExtension)

    INQUIRE( file = inFileName, exist = exist )
    if ( exist ) then
       call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
       call CheckWrittenBlock( outFileName )
    else
       allocate(auxDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
            trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
            trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
       allocate(auxFileName,source=auxDir//HamiltonianFileRootName//"0.0"//QCFileExtension)
       call OpenFile( auxFileName, uid, 'read', 'formatted' )
       call MatA%Read( uid )
       close( uid )
       call CheckBlock( MatA )
       call MatA%Transpose()
       call OpenFile( outFileName, uid, 'write', 'formatted' )
       call MatA%Write( uid )
       close( uid )
    end if

  end subroutine ClassElectronicSpaceAcquireHamiltonian


  subroutine ClassElectronicSpaceAcquireKinEnergy( Space, SymLabel, PWCBraLabel, PWCKetLabel, inFileName )
    class(ClassElectronicSpace), intent(in) :: Space
    character(len=*)           , intent(in) :: SymLabel
    character(len=*)           , intent(in) :: PWCBraLabel
    character(len=*)           , intent(in) :: PWCKetLabel
    character(len=*)           , intent(in) :: inFileName
    character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
    logical :: exist
    type(ClassMatrix) :: MatA
    integer :: uid

    call ErrorMessage("Must check validity of "//trim(SymLabel)//", " //trim(PWCBraLabel)//", and " //trim(PWCKetLabel))

    allocate(outDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
         trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
         trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

    call Execute_Command_Line(" mkdir -p "//OutDir)

    allocate(outFileName,source=outDir//KinEnergyFileRootName//"0.0"//QCFileExtension)

    INQUIRE( file = inFileName, exist = exist )
    if ( exist ) then
       call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
       call CheckWrittenBlock( outFileName )
    else
       allocate(auxDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
            trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
            trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
       allocate(auxFileName,source=auxDir//KinEnergyFileRootName//"0.0"//QCFileExtension)
       INQUIRE( file = auxFileName, exist = exist )
       if ( .not.exist ) then
          call ErrorMessage( 'The file '//&
               auxFileName//' does not exist.')
          return
       end if
       call OpenFile( auxFileName, uid, 'read', 'formatted' )
       call MatA%Read( uid )
       close( uid )
       call CheckBlock( MatA )
       call MatA%Transpose()
       call OpenFile( outFileName, uid, 'write', 'formatted' )
       call MatA%Write( uid )
       close( uid )
    end if

  end subroutine ClassElectronicSpaceAcquireKinEnergy


  subroutine ClassElectronicSpaceAcquireDipoleLen( Space, SymBraLabel, PWCBraLabel, SymKetLabel, PWCKetLabel, inFileName, mDip )
    class(ClassElectronicSpace), intent(in) :: Space
    character(len=*)           , intent(in) :: SymBraLabel
    character(len=*)           , intent(in) :: SymKetLabel
    character(len=*)           , intent(in) :: PWCBraLabel
    character(len=*)           , intent(in) :: PWCKetLabel
    character(len=*)           , intent(in) :: inFileName
    integer                    , intent(in) :: mDip
    character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
    logical :: exist
    type(ClassMatrix) :: MatA
    integer :: uid

    call ErrorMessage("Must check validity of "//trim(SymBraLabel)//", " //trim(PWCBraLabel)//&
         ", "//trim(SymKetLabel)//", and " //trim(PWCKetLabel))

    allocate(outDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
         trim(adjustl(SymBraLabel))//"_"//AddSlash(trim(adjustl(SymKetLabel)))//&
         trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

    call Execute_Command_Line(" mkdir -p "//OutDir)

    allocate(outFileName,source=outDir//DipoleLenFileRootName//"1."//AlphabeticNumber(mDip)//QCFileExtension)

    INQUIRE( file = inFileName, exist = exist )
    if ( exist ) then
       call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
       call CheckWrittenBlock( outFileName )
    else
       allocate(auxDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
            trim(adjustl(SymBraLabel))//"_"//AddSlash(trim(adjustl(SymKetLabel)))//&
            trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
       allocate(auxFileName,source=auxDir//DipoleLenFileRootName//"1."//AlphabeticNumber(mDip)//QCFileExtension)
       call OpenFile( auxFileName, uid, 'read', 'formatted' )
       call MatA%Read( uid )
       close( uid )
       call CheckBlock( MatA )
       call MatA%Transpose()
       call OpenFile( outFileName, uid, 'write', 'formatted' )
       call MatA%Write( uid )
       close( uid )
    end if

  end subroutine ClassElectronicSpaceAcquireDipoleLen


  subroutine ClassElectronicSpaceAcquireDipoleVel( Space, SymBraLabel, PWCBraLabel, SymKetLabel, PWCKetLabel, inFileName, mDip )
    class(ClassElectronicSpace), intent(in) :: Space
    character(len=*)           , intent(in) :: SymBraLabel
    character(len=*)           , intent(in) :: SymKetLabel
    character(len=*)           , intent(in) :: PWCBraLabel
    character(len=*)           , intent(in) :: PWCKetLabel
    character(len=*)           , intent(in) :: inFileName
    integer                    , intent(in) :: mDip
    character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
    logical :: exist
    type(ClassMatrix) :: MatA
    integer :: uid

    call ErrorMessage("Must check validity of "//trim(SymBraLabel)//", " //trim(PWCBraLabel)//&
         ", "//trim(SymKetLabel)//", and " //trim(PWCKetLabel))

    allocate(outDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
         trim(adjustl(SymBraLabel))//"_"//AddSlash(trim(adjustl(SymKetLabel)))//&
         trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

    call Execute_Command_Line(" mkdir -p "//OutDir)

    allocate(outFileName,source=outDir//DipoleVelFileRootName//"1."//AlphabeticNumber(mDip)//QCFileExtension)

    INQUIRE( file = inFileName, exist = exist )
    if ( exist ) then
       call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
       call CheckWrittenBlock( outFileName )
    else
       allocate(auxDir,source=Space%GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
            trim(adjustl(SymBraLabel))//"_"//AddSlash(trim(adjustl(SymKetLabel)))//&
            trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
       allocate(auxFileName,source=auxDir//DipoleVelFileRootName//"1."//AlphabeticNumber(mDip)//QCFileExtension)
       call OpenFile( auxFileName, uid, 'read', 'formatted' )
       call MatA%Read( uid )
       close( uid )
       call CheckBlock( MatA )
       call MatA%Transpose()
       call MatA%Multiply(-1.d0)
       call OpenFile( outFileName, uid, 'write', 'formatted' )
       call MatA%Write( uid )
       close( uid )
    end if

  end subroutine ClassElectronicSpaceAcquireDipoleVel
  
  


  subroutine ClassElectronicSpaceSetMultiplicity( Space, Multiplicity )
    class(ClassElectronicSpace), intent(inout) :: Space
    integer                    , intent(in)    :: Multiplicity
    Space%Multiplicity = Multiplicity
    Space%MULTIPLICITY_IS_SET = .TRUE.
  end subroutine ClassElectronicSpaceSetMultiplicity


  subroutine ClassElectronicSpaceSetStorageDir( Space, StorageDir )
    class(ClassElectronicSpace), intent(inout) :: Space
    character(len=*)           , intent(in)    :: StorageDir
    allocate( Space%StorageDir, source = StorageDir )
  end subroutine ClassElectronicSpaceSetStorageDir



  subroutine ClassElectronicSpaceSetNuclearLabel( Space, Label )
    class(ClassElectronicSpace), intent(inout) :: Space
    character(len=*)           , intent(in)    :: Label
    if( allocated(Space%NuclearLabel) ) deallocate( Space%NuclearLabel )
    allocate( Space%NuclearLabel, source = Label )
  end subroutine ClassElectronicSpaceSetNuclearLabel
  

  function ClassElectronicSpaceGetGroup( Space ) result(Group)
    class(ClassElectronicSpace), intent(in) :: Space
    type(ClassGroup), pointer :: Group
    allocate(Group,source=Space%Group)
  end function ClassElectronicSpaceGetGroup


  function ClassElectronicSpaceGetPionList( Space ) result(IonList)
    class(ClassElectronicSpace), intent(in) :: Space
    type(ClassParentIon), pointer :: IonList(:)
    allocate(IonList,source=Space%ParentIonList)
  end function ClassElectronicSpaceGetPionList


  function ClassElectronicSpaceGetLmax( Space ) result(Lmax)
    class(ClassElectronicSpace), intent(in) :: Space
    integer :: Lmax
    Lmax=Space%LMax
  end function ClassElectronicSpaceGetLmax


  !> Gets the maximum L defined in the Close-Coupling anzat
  !! for the current ClassElectronicSpace%
  function ClassElectronicSpaceGetMaxLinCC( Space ) result(L)
    class(ClassElectronicSpace), intent(in) :: Space
    integer :: L
    integer :: i
    !
    if ( .not.allocated(Space%IrrepSpaceVec) ) call Assert( &
         'ClassElectronicSpace has not been '//&
         'initialized, impossible to get the maximum L in '//&
         'the Close-Coupling expansion, aborting ...' )
    !
    L = 0
    do i = 1, Space%NIrreps
       L = max(L,Space%IrrepSpaceVec(i)%GetMaxLinCC())
    end do
  end function ClassElectronicSpaceGetMaxLinCC


  function ClassElectronicSpaceGetPionCharge( Space ) result(PionCharge)
    class(ClassElectronicSpace), intent(in) :: Space
    integer :: PionCharge
    PionCharge=Space%ParentIonCharge
  end function ClassElectronicSpaceGetPionCharge
    
  function ClassElectronicSpaceGetStorageDir( Space ) result(StorageDir)
    class(ClassElectronicSpace), intent(in) :: Space
    character(len=:), allocatable :: StorageDir
    allocate(StorageDir,source=Space%StorageDir)
  end function ClassElectronicSpaceGetStorageDir


  subroutine ClassElectronicSpaceParseConfigFile( Space, &
       FileName, &
       OnlyPrepareData, VERBOUS_ )

    class(ClassElectronicSpace), intent(inout) :: Space
    character(len=*)           , intent(in)    :: FileName
    logical, optional          , intent(in)    :: OnlyPrepareData
    integer, optional          , intent(in)    :: VERBOUS_

    character(len=:), allocatable :: FullText, GroupLabel, GeneralLabel, irrepLabel, ParentIonList
    character(len=:), allocatable :: snRefOccupied
    type(ClassIrrep), pointer     :: irrepv(:)
    integer                       :: iIrrep, ichar, iflag, iSym, iIon, iIrr, iostat
    character(len=:), allocatable :: sIon, sOrb
    character(len=*), parameter   :: EOLN_STRN = "<\EOLN>"
    integer :: VERBOUS
    VERBOUS = 0
    if(present(VERBOUS_))VERBOUS=VERBOUS_
    

    !.. Load the Configuration file
    call GetFullText( FileName, FullText, EOLN_STRN )
    call SetStringToUppercase( FullText )

    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "GROUP", GroupLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GroupLabel)) call Assert("Group label missing in "//trim(FileName))
    write(*,"(a)") " Group               : "//GroupLabel
    call Space%Group%init( GroupLabel )

    call FetchGlobalVariable( FullText, "PARENT_ION_CHARGE"    , GeneralLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GeneralLabel)) call Assert("Charge label missing in "//trim(FileName))
    write(*,"(a)") " Parent-ion Charge   : "//GeneralLabel
    read(GeneralLabel,*) Space%ParentIonCharge

    call FetchGlobalVariable( FullText, "LMAX"      , GeneralLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GeneralLabel)) call Assert("lmax label missing in "//trim(FileName))
    write(*,"(a)") " Lmax                : "//GeneralLabel
    read(GeneralLabel,*) Space%Lmax

    call FetchGlobalVariable( FullText, "NELECTRONS", GeneralLabel, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(GeneralLabel)) call Assert("Nelectrons label missing in "//trim(FileName))
    write(*,"(a)") " Tot. N. Electrons   : "//GeneralLabel
    read(GeneralLabel,*) Space%Ne

    call FetchGlobalVariable( FullText, "PARENT_ION_LIST", ParentIonList, EOLN_STRN = EOLN_STRN )
    if(.not.allocated(ParentIonList)) call Assert("Parent-ion list missing in "//trim(FileName))
    write(*,"(a)") " List of parent ions : "//trim(ParentIonList)
    Space%NParentIons = nTokens(ParentIonList)
    if(associated(Space%ParentIonList))deallocate(Space%ParentIonList)
    allocate(Space%ParentIonList(Space%NParentIons))
    do iIon = 1, Space%NParentIons
       call GetToken( ParentIonList, iIon, sIon )
       call Space%ParentIonList(iIon)%Init(Space%Group,sIon,Space%ParentIonCharge)
       if(VERBOUS>0)call Space%ParentIonList(iIon)%Show()
    enddo

    call FetchGlobalVariable( FullText, "N_REF_OCCUPIED", snRefOccupied, EOLN_STRN = EOLN_STRN )
    allocate(Space%nRefOccupied(Space%Group%GetNIrreps()))
    if(.not.allocated(snRefOccupied))then
       Space%nRefOccupied=0
    else
       write(*,"(a)") " HF orbital occupancy: "//trim(snRefOccupied)
       if( Space%Group%GetNIrreps() /= nTokens(snRefOccupied) )then
          write(*,"(a)")"the number of N_REF_OCCUPIED orbitals does not match the number of irreps"
          stop
       endif
       do iIrr = 1, size(Space%nRefOccupied)
          call GetToken( snRefOccupied, iIrr, sOrb )
          read(sOrb,*,iostat=iostat)Space%nRefOccupied(iIrr)
          if(iostat/=0)then
             write(*,"(a)")"syntax error in N_REF_OCCUPIED"
             stop
          endif
       enddo
    endif

    if(Space%MULTIPLICITY_IS_SET)then


       !.. Determines the number of symmetric spaces
       irrepv => Space%Group%GetIrrepList()
       Space%NIrreps = 0
       write(*,"(a)",advance="no") " List of Symmetries  :"
       do iIrrep = 1, size( irrepv )
          irrepLabel = irrepv(iIrrep)%GetName()
          call SetStringToUppercase( irrepLabel )
          ichar = index( FullText, "["//AlphabeticNumber(Space%Multiplicity)//irrepLabel//"]" )
          if(ichar>0)then
             Space%NIrreps = Space%NIrreps + 1
             write(*,"(a)",advance="no") " "//AlphabeticNumber(Space%Multiplicity)//irrepLabel
          endif
       enddo
       write(*,*)
       
       !.. Parse each Symmetric Electronic Space
       allocate( Space%IrrepSpaceVec( Space%NIrreps ) )

       iSym = 0
       do iIrrep = 1, size( irrepv )

          irrepLabel = irrepv(iIrrep)%GetName()
          call SetStringToUppercase( irrepLabel )
          ichar = index( FullText, "["//AlphabeticNumber(Space%Multiplicity)//irrepLabel//"]" )
          if(ichar<1) cycle

          iSym = iSym + 1

          if ( allocated( Space%NuclearLabel ) ) &
               call Space%IrrepSpaceVec( iSym )%SetNuclearLabel( Space%NuclearLabel )
          call Space%IrrepSpaceVec( iSym )%SetGroup( Space%Group )
          call Space%IrrepSpaceVec( iSym )%SetIrrep( irrepv( iIrrep ) )
          call Space%IrrepSpaceVec( iSym )%SetLSIrrep( irrepv( iIrrep ) )
          call Space%IrrepSpaceVec( iSym )%SetMultiplicity( Space%Multiplicity )
          call Space%IrrepSpaceVec( iSym )%SetStorage( Space%StorageDir )
          !
          call Space%IrrepSpaceVec( iSym )%ParseConfigFile( FileName, iflag, OnlyPrepareData )
          if(iflag /= 0 )call Assert( "Invalid syntax for "//&
               Space%IrrepSpaceVec( iSym )%GetLabel()//" in "//trim(FileName) )

       enddo

    endif

  end subroutine ClassElectronicSpaceParseConfigFile


  !> Retrieves wether the electronic space has been initialized or not.
  logical function ClassElectronicSpaceInitialized( Self ) result(res)
    !> Class of the electronic space.
    class(ClassElectronicSpace), intent(in) :: Self
    !
    res = .false.
    !
    if ( allocated(Self%IrrepSpaceVec) ) res = .true.
    !
  end function ClassElectronicSpaceInitialized



  subroutine ClassElectronicSpaceShow( Self, unit )
    !> Class of the electronic space.
    class(ClassElectronicSpace), intent(in) :: Self
    integer, optional          , intent(in) :: unit
    integer :: i ,outunit
    !
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    !
    write(outunit,"(a,i4)"          ) "  Space Multiplicity :",Self%Multiplicity
    write(outunit,"(a,i4)"          ) "  Parent-ion Charge :",Self%ParentIonCharge
    write(outunit,"(a,i4)"          ) "  Maximum Angular Momentum :",Self%Lmax
    write(outunit,"(a,i4)"          ) "  Number of Electrons :",Self%Ne
    call Self%Group%Show()
    write(outunit,"(a,i4)"          ) "  Number of Irreps :",Self%NIrreps
    do i = 1, Self%NIrreps
       call Self%IrrepSpaceVec(i)%Show()
    end do
    write(outunit,"(a,a)"          ) "  Space Storage Dir :",Self%StorageDir
    write(outunit,"(a,a)"          ) "  Space Nuclear Label :",Self%NuclearLabel
    !
  end subroutine ClassElectronicSpaceShow


  !> Gets the number of irreducible representations of the point group used in the close-coupling expansion.
  integer function ClassElectronicSpaceGetNIrreps( Self ) result(N)
    !> Class of the electronic space.
    class(ClassElectronicSpace), intent(in) :: Self
    !
    if ( .not.Self%Initialized() ) call Assert( "Impossible to get the number of irreducible representations of the "//&
         "electronic space because hast not been initialized." )
    N = Self%NIrreps
    !
  end function ClassElectronicSpaceGetNIrreps


  !> Gets the number of doubly occupied orbitals per each irrep, in the reference configuration
  !! which is assumed to be closed shell
  function ClassElectronicSpaceGetnRefOccupied( Self ) result(nRefOccupied)
    !> Class of the electronic space.
    class(ClassElectronicSpace), intent(in) :: Self
    integer, allocatable :: nRefOccupied(:)
    !
    if ( .not.Self%Initialized() ) call Assert( "electronic space has not been initialized." )
    allocate(nRefOccupied, source = Self%nRefOccupied)
    !
  end function ClassElectronicSpaceGetnRefOccupied


  !> Gets the electronic space group name.
  function ClassElectronicSpaceGetGroupName( Self ) result(GName)
    !> Class of the electronic space.
    class(ClassElectronicSpace), intent(in) :: Self
    character(len=:), allocatable :: GName
    !
    if ( .not.Self%Initialized() ) call Assert( "Impossible to get the group name of the "//&
         "electronic space because hast not been initialized." )
    allocate( GName, source = Self%Group%GetName() )
    !
  end function ClassElectronicSpaceGetGroupName


  !> Get from the electronic space, the symmetric electronic space giving the corresponding index in the list.
  function ClassElectronicSpaceGetSymElectSpace( Self, iIrrep ) result( SymSpace )
    !> Class of the electronic space.
    class(ClassElectronicSpace), target, intent(in) :: Self
    !> Index correspondind to the symmetric electronic space in the list.
    integer                    , intent(in) :: iIrrep
    !> Class of the symmetric electronic space.
    type(ClassSymmetricElectronicSpace), pointer :: SymSpace
    !
    if ( .not.Self%Initialized() ) call Assert( "To get a symmetric electronic space, this must be initialized." )
    !
    allocate( SymSpace, source = Self%IrrepSpaceVec(iIrrep) )
    !
  end function ClassElectronicSpaceGetSymElectSpace



  !> Get from the electronic space, the symmetric electronic space giving the irreducible representation name.
  function ClassElectronicSpaceGetSymElectSpaceByIrrepName( Self, IrrepName ) result( SymSpace )
    !> Class of the electronic space.
    class(ClassElectronicSpace), target, intent(in) :: Self
    !> Index correspondind to the symmetric electronic space in the list.
    character(len=*)                   , intent(in) :: IrrepName
    !> Class of the symmetric electronic space.
    type(ClassSymmetricElectronicSpace), pointer :: SymSpace
    integer :: i
    !
    if ( .not.Self%Initialized() ) call Assert( "To get a symmetric electronic space, this must be initialized." )
    !
    do i = 1, Self%GetNIrreps()
       if ( Self%IrrepSpaceVec(i)%GetIrrepLabel() .is. IrrepName ) then
          allocate( SymSpace, source = Self%IrrepSpaceVec(i) )
          return
       end if
    end do
    !
  end function ClassElectronicSpaceGetSymElectSpaceByIrrepName


  subroutine ClassElectronicSpaceCheckSymmetry( Space, SymLabel )
    class(ClassElectronicSpace), target, intent(in) :: Space
    character(len=*)                   , intent(in) :: SymLabel
    !
    integer :: NumSymElectSpace, i
    type(ClassSymmetricElectronicSpace), pointer :: SymElectSpace
    logical :: PresentIrrep
    !
    NumSymElectSpace = Space%GetNIrreps()
    PresentIrrep = .false.
    do i = 1, NumSymElectSpace
       allocate( SymElectSpace, source = Space%GetSymElectSpace(i) )
       if ( SymLabel .is. SymElectSpace%GetIrrepLabel() ) then
          PresentIrrep = .true.
          return
       end if
       deallocate( SymElectSpace )
    end do
    if ( .not.PresentIrrep )then
       call ErrorMessage( &
         ' The irreducible representation '//SymLabel//&
         ' is not present in the close-coupling expansion for the '//Space%GetGroupName()//&
         ' group.' )
       stop
    endif
  end subroutine ClassElectronicSpaceCheckSymmetry

  

  function GetCloseCouplingDir() result(Dir)
    character(len=:), allocatable :: Dir
    allocate( Dir, source = CLOSE_COUPLING_DIR )
  end function GetCloseCouplingDir


  function GetQCFileExtension() result(Dir)
    character(len=:), allocatable :: Dir
    allocate( Dir, source = QCFileExtension )
  end function GetQCFileExtension


  !> By convention, if a block stored has a 0 dimension, this is
  !! changed by 1, anf fill with zeros.
  subroutine CheckReadBlock( Array, ArrayIsFine )
    real(kind(1d0)), allocatable, intent(inout) :: Array(:,:)
    logical, optional           , intent(out)   :: ArrayIsFine
    !
    integer :: NRows, NCols
    logical :: BlockIsZero
    integer, parameter :: NewDimension = 1
    !
    BlockIsZero = .false.
    NRows = size(Array,1)
    NCols = size(Array,2)
    !
    if ( NRows == 0 ) then
       NRows = NewDimension
       BlockIsZero = .true.
    end if
    !
    if ( NCols == 0 ) then
       NCols = NewDimension
       BlockIsZero = .true.
    end if
    !
    if ( BlockIsZero ) then
       deallocate( Array )
       allocate( Array(NRows,NCols) )
       Array = 0.d0
    end if
    !
    if ( present(ArrayIsFine) ) then
       if ( BlockIsZero ) then
          ArrayIsFine = .false.
       else
          ArrayIsFine = .true.
       end if
    end if
    !
  end subroutine CheckReadBlock

  
  !> By convention, if a block stored has a 0 dimension, this is
  !! changed by 1, anf fill with zeros.
  subroutine CheckBlock( Mat )
    class(ClassMatrix), intent(inout) :: Mat
    !
    integer :: NRows, NCols
    logical :: BlockIsZero
    integer, parameter :: NewDimension = 1
    !
    BlockIsZero = .false.
    NRows = Mat%NRows()
    NCols = Mat%NColumns()
    !
    if ( NRows == 0 ) then
       NRows = NewDimension
       BlockIsZero = .true.
    end if
    !
    if ( NCols == 0 ) then
       NCols = NewDimension
       BlockIsZero = .true.
    end if
    !
    if ( BlockIsZero ) then
       call Mat%Free()
       call Mat%InitFull( NRows, NCols )
    end if
    !
  end subroutine CheckBlock



  subroutine CheckWrittenBlock( FileName )
    character(len=*), intent(in) :: FileName
    real(kind(1d0)), allocatable :: Array(:,:)
    logical :: ArrayIsFine
    integer :: uid
    type(ClassMatrix) :: Mat
    call ReadMatrix( FileName , Array )
    call CheckReadBlock( Array )
    Mat = Array
    deallocate( Array )
    call OpenFile( FileName, uid, 'write', 'formatted' )
    call Mat%Write( uid )
    close( uid )
  end subroutine CheckWrittenBlock





end module ModuleElectronicSpace
