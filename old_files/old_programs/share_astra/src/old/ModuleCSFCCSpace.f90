module ModuleCSFCCSpace
  
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleIO
  use ModuleParentIons
  
  
  implicit none

  private

  character(len=*), parameter :: AUG_SUBDIR = "aug/"
  character(len=*), parameter :: ION_SUBDIR = "ion/"
  character(len=*), parameter :: ONE_SUBDIR = "one/"
  integer         , parameter :: MAX_NLW_PER_PION_TERM = 32
  
  type, public :: ClassCSFCCSpace
     
     !private
     
     !> Storage directory
     character(len=:), allocatable :: StoreDir

     !> Name of the system. E.g., Neon
     character(len=:), allocatable :: name

     !> multiplicity 2S+1
     integer  :: multiplicity
     
     !> Nuclear Charge 
     real(kind(1d0))  :: NuclearCharge

     !> maximum value for the principal quantum number in the active space
     !> it is assumed to be the same as the maximum n for the reference configurations
     integer  :: nmax

     !> maximum value for the orbital angular momentum in the active space
     integer  :: lmax

     !> number of shells nl, with n <= nmax, l <= min(lmax, n-1)
     integer  :: nshells

     !> maximum number of excitations from the references
     integer  :: nexcite

     !> minimum number of electrons per subshell: minnel(in), in = 1, 2, ..., nmax - 1
     !> sub-shell kind: closed ("c"), active ("c"), with at least n electrons ("<n>")
     !! sshKind( iShell ), iShell = 1, 2, ... , nshells  
     integer  , allocatable :: minnel(:)
     character, allocatable :: sshKind(:)


     !> number of reference configurations for the parent ions
     !> reference configurations for the parent ions
     !> number of parent-ion terms
     !> list of parent-ion terms
     integer              :: nIonRefs
     integer, allocatable :: ionRefs(:,:)
     integer              :: nIonTerms
     character(len=3), allocatable :: ionTerms(:)


     !> number of reference configurations for the neutrals
     !> reference configurations for the neutral terms
     !> number of neutral terms
     !> list of neutral terms 
     integer              :: nAugRefs
     integer, allocatable :: AugRefs(:,:)
     integer              :: nAugTerms
     character(len=3), allocatable :: AugTerms(:)

     !> maximum orbital angular momentum for the photoelectron
     integer  :: lwmax 

   contains
     
     generic, public :: ParseConfigFile    => ClassCSFCCSpaceParseConfigFile
     generic, public :: Show                => ClassCSFCCSpaceShow
     generic, public :: GenerateList        => ClassCSFCCSpaceGenerateList

     !.. Private procedures
     generic  , private :: GenerateListIons => ClassCSFCCSpaceGenerateListIons
     generic  , private :: GenerateListLCs => ClassCSFCCSpaceGenerateListLCs
     generic  , private :: GenerateListCSFCC => ClassCSFCCSpaceGenerateListCSFCC
     generic  , private :: AugmentShell => ClassCSFCCSpaceAugmentShell
     generic  , private :: SetupDirTree => ClassCSFCCSpaceSetupDirTree
     generic  , private :: IonTermDir   => ClassCSFCCSpaceIonTermDir
     generic  , private :: AugTermDir   => ClassCSFCCSpaceAugTermDir

     procedure, private :: ClassCSFCCSpaceParseConfigFile
     procedure, private :: ClassCSFCCSpaceShow
     procedure, private :: ClassCSFCCSpaceGenerateList

     procedure, private :: ClassCSFCCSpaceGenerateListIons
     procedure, private :: ClassCSFCCSpaceGenerateListLCs
     procedure, private :: ClassCSFCCSpaceGenerateListCSFCC
     procedure, private :: ClassCSFCCSpaceAugmentShell
     procedure, private :: ClassCSFCCSpaceSetupDirTree
     procedure, private :: ClassCSFCCSpaceIonTermDir
     procedure, private :: ClassCSFCCSpaceAugTermDir

  end type ClassCSFCCSpace


  
contains


  subroutine ClassCSFCCSpaceParseConfigFile( Space, &
       FileName )
    
    class(ClassCSFCCSpace), intent(inout) :: Space
    character(len=*)      , intent(in)    :: FileName
    
    character(len=:), allocatable :: FullText, strnBuf, subText
    integer :: in, il, iShell
    integer :: i1, i2
    integer :: iostat
    integer :: iRef, iTerm

    !.. Load the Configuration file
    call GetFullText( FileName, FullText )
    !call SetStringToUppercase( FullText )
    
    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "NAME", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("NAME missing in "//trim(FileName))
    allocate( Space.name, source=trim(adjustl(strnBuf)) )

    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "STORE_DIR", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("STORE_DIR missing in "//trim(FileName))
    i2=len_trim(strnBuf)
    if(strnBuf(i2:i2)=="/")then
       allocate( Space.StoreDir, source=trim(adjustl(strnBuf)) )
    else
       allocate( Space.StoreDir, source=trim(adjustl(strnBuf))//"/" )
    endif

    !.. Determines the general variables
    call FetchGlobalVariable( FullText, "NUCLEAR_CHARGE", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("NUCLEAR_CHARGE missing in "//trim(FileName))
    read(strnBuf,*) Space.NuclearCharge
    if(Space.NuclearCharge < 0) call Assert("Negative nuclear charge in "//trim(FileName))
       
    call FetchGlobalVariable( FullText, "NMAX", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("NMAX missing in "//trim(FileName))
    read(strnBuf,*) Space.nmax
    if(Space.nmax < 1) call Assert("nmax < 1 in "//trim(FileName))

    call FetchGlobalVariable( FullText, "LMAX", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("LMAX missing in "//trim(FileName))
    read(strnBuf,*) Space.lmax
    if(Space.lmax > Space.nmax-1) call Assert("lmax > nmax+1 in "//trim(FileName))
    iShell=0
    do in = 1, Space.nmax
       do il = 0, min( in - 1, Space.lmax )
          iShell=iShell+1
       enddo
    enddo
    Space.nshells=iShell

    call FetchGlobalVariable( FullText, "LWMAX", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("LWMAX missing in "//trim(FileName))
    read(strnBuf,*) Space.lwmax

    call FetchGlobalVariable( FullText, "NEXCITE", strnBuf )
    if(.not.allocated(strnBuf)) call Assert("NEXCITE missing in "//trim(FileName))
    read(strnBuf,*) Space.nexcite

    !.. Read MINNEL
    i1=index(FullText, "MINNEL")
    if(i1<=0) call Assert("MINNEL missing in "//trim(FileName))
    i2=index(FullText(i1:),"(")
    i1=i1+i2!first character after parenthesis
    i2=index(FullText(i1:),")")
    i2=i1+i2-2 !last character before parenthesis
    allocate(Space.minnel( Space.nmax - 1 ) )
    Space.minnel=0
    read(FullText(i1:i2),*,iostat=iostat) ( Space.minnel( in ), in = 1, Space.nmax - 1 )
    if(iostat/=0)call Assert("invalid format for MINNEL in "//trim(FileName))

    !.. Read SSHKIND
    i1=index(FullText, "SSHKIND")
    if(i1<=0) call Assert("SSHKIND missing in "//trim(FileName))
    i2=index(FullText(i1:),"(")
    i1=i1+i2!first character after parenthesis
    i2=index(FullText(i1:),")")
    i2=i1+i2-2 !last character before parenthesis
    allocate(Space.sshKind( Space.nshells ) )
    Space.sshKind=" "
    read(FullText(i1:i2),*,iostat=iostat) ( Space.sshKind( iShell ), iShell = 1, Space.nshells )
    if(iostat/=0)call Assert("invalid format for SSHKIND in "//trim(FileName))

    !.. Read IONREF
    i1=index(FullText, "IONREF")
    if(i1<=0) call Assert("IONREF missing in "//trim(FileName))
    i2=index(FullText(i1:),"{")
    i1=i1+i2!first character after parenthesis
    i2=index(FullText(i1:),"}")
    i2=i1+i2-2 !last character before parenthesis
    allocate(subText,source=FullText(i1:i2))
    
    !.. Count the number of references
    Space.nionrefs=0
    i1=1
    do 
       i2= index(subText(i1:),")")
       if(i2<=0)exit
       Space.nionrefs = Space.nionrefs+1
       i1=i1+i2
       if(i1>len_trim(subText))exit
    enddo

    !.. Load the ionic references
    allocate(Space.ionrefs( Space.nshells, Space.nionrefs ) )
    do iref=1, Space.nionrefs
       i2 = index(subText,")") - 1
       i1 = index(subText(:i2),"(") + 1
       read(subText(i1:i2),*,iostat=iostat) (Space.ionrefs( iShell, iRef ), iShell = 1, Space.nshells )
       subText(i1-1:i2+1)=" "
    enddo

    !.. Read IONTERMS
    i1=index(FullText, "IONTERMS")
    if(i1<=0) call Assert("IONTERMS missing in "//trim(FileName))
    i2=index(FullText(i1:),"(")
    i1=i1+i2!first character after parenthesis
    i2=index(FullText(i1:),")")
    i2=i1+i2-2 !last character before parenthesis
    !.. Count number of ion terms
    Space.nionterms=nTokens( FullText(i1:i2) )
    allocate( Space.ionTerms( Space.nionterms ) )
    do iTerm = 1, Space.nionterms
       call GetToken( FullText(i1:i2), iTerm, strnBuf )
       Space.ionTerms( iTerm ) = strnBuf
    enddo

    call Space.AugmentShell()
    
    !.. Read AUGTERMS
    i1=index(FullText, "AUGTERMS")
    if(i1<=0) call Assert("AUGTERMS missing in "//trim(FileName))
    i2=index(FullText(i1:),"(")
    i1=i1+i2!first character after parenthesis
    i2=index(FullText(i1:),")")
    i2=i1+i2-2 !last character before parenthesis
    !.. Count number of ion terms
    Space.naugterms=nTokens( FullText(i1:i2) )
    allocate( Space.augTerms( Space.naugterms ) )
    do iTerm = 1, Space.naugterms
       call GetToken( FullText(i1:i2), iTerm, strnBuf )
       Space.augTerms( iTerm ) = strnBuf
    enddo

  end subroutine ClassCSFCCSpaceParseConfigFile

  subroutine ClassCSFCCSpaceAugmentShell( Self )
    class(ClassCSFCCSpace), intent(inout) :: Self
    integer, allocatable :: tempRef(:,:), ivec(:)
    integer :: nRef, iShell, in, il, iAugRef, iRef,max_shell_el
    character(len=*), parameter :: digits="0123456789"
    allocate(tempRef( Self.nshells, Self.nionrefs * Self.nshells ) )
    allocate(ivec( Self.nshells ) )
    tempRef=0
    nRef=0
    do iRef=1,Self.nionrefs
       iShell=0
       do in=1,Self.nmax
          il_cycle : do il=0,min(in-1,Self.lmax)
             iShell=iShell+1
             if( Self.sshKind( iShell ) == "c" )cycle il_cycle
             !
             !*** Originally, the program did not augment the 
             !    ionic configuration in the inactive (computed by mchf) 
             !    orbitals. However, this is a mistake: the inactive 
             !    orbitals are still internal orbitals, and hence they
             !    are needed to define the penetration terms.
             !    On the other hand, the virtual orbitals, which for
             !    our purpose can be regarded as those orbitals that 
             !    have not been computed in mchf and that are therefore
             !    part of the external orbitals, should not be augmented
             !    because they are already represented by the [ion ;l] functions.
             !*** 
!!$          if( Self.sshKind( iShell ) == "i" )cycle il_cycle
             if( Self.sshKind( iShell ) == "I" )cycle il_cycle
             if( Self.sshKind( iShell ) == "v" )cycle il_cycle
             !
             max_shell_el=4*il+2
             !max_shell_el = index(digits,Self.sshKind( iShell ))-1
             !if(max_shell_el<0)max_shell_el=4*il+2
             if( Self.ionrefs( iShell, iRef )+1 > max_shell_el )cycle il_cycle
             ivec=Self.ionrefs(:,iRef)
             ivec(iShell)=ivec(iShell)+1
             do iAugRef =1, nRef
                if(sum(abs(tempRef(:,iAugRef)-ivec))==0)cycle il_cycle
             enddo
             nRef=nRef+1
             tempRef(:,nRef)=ivec
          enddo il_cycle
       enddo
    enddo
    Self.nAugRefs = nRef
    allocate(Self.AugRefs(self.nshells,self.nAugRefs))
    Self.AugRefs=tempRef(:,1:nRef)
    deallocate(ivec,tempRef)

  end subroutine ClassCSFCCSpaceAugmentShell
  

  subroutine ClassCSFCCSpaceShow( Self, unit )
    !> Class of the electronic space.
    class(ClassCSFCCSpace), intent(in) :: Self
    integer, optional     , intent(in) :: unit
    integer :: outunit, iRef, ish
    !
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    !
    write(outunit,"(a)"        ) "  System Name..............: "//Self.Name
    write(outunit,"(a)"        ) "  System Name..............: "//Self.StoreDir
    write(outunit,"(a,i0)"     ) "  CSF CC Multiplicity......: ", Self.Multiplicity
    write(outunit,"(a,f3.0)"   ) "  Nuclear Charge...........: ", Self.NuclearCharge
    write(outunit,"(a,i0)"     ) "  n max active space.......: ", Self.nmax
    write(outunit,"(a,i0)"     ) "  l max active space.......: ", Self.lmax
    write(outunit,"(a,i0)"     ) "  n subshells..............: ", Self.nshells
    write(outunit,"(a,*(i0,x))") "  min n electrons per shell: ", Self.minnel
    write(outunit,"(a,*(a,x))" ) "  subshell kind............: ", Self.sshKind
    write(outunit,"(a,i0)"     ) "  n ion reference configs..: ", Self.nionrefs
    do iRef = 1, Self.nionrefs
       write(outunit,"(a,i2,a,*(i0,x))") "  ion ref config [",iRef,"]......: ",&
            (Self.ionRefs(iSh,iRef),iSh=1,Self.nshells)
    enddo
    write(outunit,"(a,i0)"     ) "  n excitations............: ", Self.nexcite
    write(outunit,"(a,i0)"     ) "  n ion terms..............: ", Self.nIonTerms
    write(outunit,"(a,*(a,x))" ) "  Ion terms................: ", Self.ionTerms
    write(outunit,"(a,i0)"     ) "  n Augmented configs......: ", Self.nAugRefs
    write(outunit,"(a,i0)"     ) "  l max electron wave......: ", Self.lwmax
    do iRef = 1, Self.nAugRefs
       write(outunit,"(a,i2,a,*(i0,x))") "  Aug ref config [",iRef,"]......: ",&
            (Self.AugRefs(iSh,iRef),iSh=1,Self.nshells)
    enddo
    write(outunit,"(a,i0)"     ) "  n Aug terms..............: ", Self.nAugTerms
    write(outunit,"(a,*(a,x))" ) "  Aug terms................: ", Self.AugTerms
    !
  end subroutine ClassCSFCCSpaceShow


  function ClassCSFCCSpaceIonTermDir( Self, iTerm ) result(strnBuf)
    class(ClassCSFCCSpace), intent(in)  :: Self
    integer               , intent(in)  :: iTerm
    character(len=:)      , allocatable :: strnBuf
    allocate(strnBuf,source=Self.StoreDir//ION_SUBDIR//Self.ionTerms(iTerm)//"/")
  end function ClassCSFCCSpaceIonTermDir

  function ClassCSFCCSpaceAugTermDir( Self, iTerm ) result(strnBuf)
    class(ClassCSFCCSpace), intent(in)  :: Self
    integer               , intent(in)  :: iTerm
    character(len=:)      , allocatable :: strnBuf
    allocate(strnBuf,source=Self.StoreDir//AUG_SUBDIR//Self.AugTerms(iTerm)//"/")
  end function ClassCSFCCSpaceAugTermDir

  subroutine ClassCSFCCSpaceSetupDirTree( Self )
    class(ClassCSFCCSpace), intent(in) :: Self
    integer :: iTerm
    call system("mkdir -p "//Self.StoreDir)
    call system("mkdir -p "//Self.StoreDir//AUG_SUBDIR)
    call system("mkdir -p "//Self.StoreDir//ION_SUBDIR)
    call system("mkdir -p "//Self.StoreDir//ONE_SUBDIR)
    do iTerm = 1, Self.nIonTerms
       call system("mkdir -p "//Self.IonTermDir( iTerm ) )
    enddo
    do iTerm = 1, Self.nAugTerms
       call system("mkdir -p "//Self.AugTermDir( iTerm ) )
    enddo
  end subroutine ClassCSFCCSpaceSetupDirTree


  subroutine ClassCSFCCSpaceGenerateList( Self )
    class(ClassCSFCCSpace), intent(in) :: Self
    !
    call Self.SetupDirTree
    !
    call Self.GenerateListIons
    call Self.GenerateListLCs
    call Self.GenerateListCSFCC
    !
  end subroutine ClassCSFCCSpaceGenerateList


  subroutine ClassCSFCCSpaceGenerateListIons( Self )
    use ModuleCSFLists
    class(ClassCSFCCSpace), intent(in) :: Self
    
    character(len=3) :: act
    integer :: iTerm, iRef
    
    !.. Generate list for the parent ions
    do iTerm=1,Self.nIonTerms
       act="new"
       call system("rm "//Self.IonTermDir(iTerm)//Self.ionTerms( iTerm )//".cfg" )
       do iRef = 1, Self.nIonRefs
          call MakeCSFList( act, &
               Self.nmax, &
               Self.lmax, &
               Self.minnel, &
               Self.nmax, &
               Self.ionRefs(:,iRef), &
               Self.sshKind, &
               Self.ionTerms( iTerm ), &
               Self.nexcite, &
               Self.IonTermDir(iTerm)//Self.ionTerms( iTerm )//".cfg" )
          act="add"
       enddo

       !***
       !STOP
       !***

    enddo


  end subroutine ClassCSFCCSpaceGenerateListIons


  subroutine ClassCSFCCSpaceGenerateListLCs( Self )
    use ModuleCSFLists
    class(ClassCSFCCSpace), intent(in) :: Self
    
    character(len=3) :: act
    integer :: iTerm, iRef
    
    !.. Generate list for the localized channel
    do iTerm=1,Self.nAugTerms
       act="new"
       call system("rm "//Self.AugTermDir(iTerm)//Self.AugTerms( iTerm )//"_LC.cfg")
       do iRef = 1, Self.nAugRefs
          call MakeCSFList( act, &
               Self.nmax, &
               Self.lmax, &
               Self.minnel, &
               Self.nmax, &
               Self.AugRefs(:,iRef), &
               Self.sshKind, &
               Self.AugTerms( iTerm ), &
               Self.nexcite, &
               Self.AugTermDir(iTerm)//Self.AugTerms( iTerm )//"_LC.cfg" )
          act="add"
          !write(*,"(a,*(x,i0))") Self.AugTerms( iTerm ),Self.AugRefs( :, iRef )
       enddo
    enddo

  end subroutine ClassCSFCCSpaceGenerateListLCs


  subroutine ClassCSFCCSpaceGenerateListCSFCC( Self )
    use ModuleCSFLists
    class(ClassCSFCCSpace), intent(in) :: Self
    
    integer                       :: iAugTerm, iIonTerm
    character(len=:), allocatable :: FileIon, FileCC

    !.. Generate list for the localized channel
    do iAugTerm=1,Self.nAugTerms
       !
       do iIonTerm=1,Self.nIonTerms
          !
          allocate( FileIon, source = &
               Self.IonTermDir(iIonTerm)//&
               Self.ionTerms( iIonTerm )//".cfg")
          !
          allocate( FileCC,  source = &
               Self.AugTermDir(iAugTerm)//&
               Self.AugTerms( iAugTerm )//"_"//&
               Self.ionTerms( iIonTerm )//".cfg" )
          !
          call AugmentCSFList( &
               FileIon, &
               Self.IonTerms( iIonTerm ), &
               Self.lwmax, &
               Self.AugTerms( iAugTerm ), &
               FileCC )
          !
          deallocate(FileIon)
          deallocate(FileCC)
          !
       enddo
       !
    enddo

  end subroutine ClassCSFCCSpaceGenerateListCSFCC


  ! subroutine ClassCSFCCSpaceAcquireHamiltonian( Space, SymLabel, PWCBraLabel, PWCKetLabel, inFileName )
  !   class(ClassCSFCCSpace), intent(in) :: Space
  !   character(len=*)           , intent(in) :: SymLabel
  !   character(len=*)           , intent(in) :: PWCBraLabel
  !   character(len=*)           , intent(in) :: PWCKetLabel
  !   character(len=*)           , intent(in) :: inFileName
  !   character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
  !   logical :: exist
  !   type(ClassMatrix) :: MatA
  !   integer :: uid

  !   call ErrorMessage("Must check validity of "//trim(SymLabel)//", " //trim(PWCBraLabel)//", and " //trim(PWCKetLabel))

  !   allocate(outDir,source=Space.GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
  !        trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
  !        trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

  !   call Execute_Command_Line(" mkdir -p "//OutDir)

  !   allocate(outFileName,source=outDir//HamiltonianFileRootName//"0.0"//QCFileExtension)

  !   INQUIRE( file = inFileName, exist = exist )
  !   if ( exist ) then
  !      call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
  !      call CheckWrittenBlock( outFileName )
  !   else
  !      allocate(auxDir,source=Space.GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
  !           trim(adjustl(SymLabel))//"_"//AddSlash(trim(adjustl(SymLabel)))//&
  !           trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
  !      allocate(auxFileName,source=auxDir//HamiltonianFileRootName//"0.0"//QCFileExtension)
  !      call OpenFile( auxFileName, uid, 'read', 'formatted' )
  !      call MatA.Read( uid )
  !      close( uid )
  !      call CheckBlock( MatA )
  !      call MatA.Transpose()
  !      call OpenFile( outFileName, uid, 'write', 'formatted' )
  !      call MatA.Write( uid )
  !      close( uid )
  !   end if

  ! end subroutine ClassCSFCCSpaceAcquireHamiltonian


  ! subroutine ClassCSFCCSpaceAcquireDipoleLen( Space, SymBraLabel, PWCBraLabel, SymKetLabel, PWCKetLabel, inFileName, mDip )
  !   class(ClassCSFCCSpace), intent(in) :: Space
  !   character(len=*)           , intent(in) :: SymBraLabel
  !   character(len=*)           , intent(in) :: SymKetLabel
  !   character(len=*)           , intent(in) :: PWCBraLabel
  !   character(len=*)           , intent(in) :: PWCKetLabel
  !   character(len=*)           , intent(in) :: inFileName
  !   integer                    , intent(in) :: mDip
  !   character(len=:), allocatable :: outFileName, outDir, auxDir, auxFileName
  !   logical :: exist
  !   type(ClassMatrix) :: MatA
  !   integer :: uid

  !   call ErrorMessage("Must check validity of "//trim(SymBraLabel)//", " //trim(PWCBraLabel)//&
  !        ", "//trim(SymKetLabel)//", and " //trim(PWCKetLabel))

  !   allocate(outDir,source=Space.GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
  !        trim(adjustl(SymBraLabel))//"_"//AddSlash(trim(adjustl(SymKetLabel)))//&
  !        trim(adjustl(PWCBraLabel))//"_"//AddSlash(trim(adjustl(PWCKetLabel))))

  !   call Execute_Command_Line(" mkdir -p "//OutDir)

  !   allocate(outFileName,source=outDir//DipoleLenFileRootName//"1."//AlphabeticNumber(mDip)//QCFileExtension)

  !   INQUIRE( file = inFileName, exist = exist )
  !   if ( exist ) then
  !      call Execute_Command_Line(" cp "//trim(inFileName)//" "//outFileName)
  !      call CheckWrittenBlock( outFileName )
  !   else
  !      allocate(auxDir,source=Space.GetStorageDir()//AddSlash(CLOSE_COUPLING_DIR)//&
  !           trim(adjustl(SymBraLabel))//"_"//AddSlash(trim(adjustl(SymKetLabel)))//&
  !           trim(adjustl(PWCKetLabel))//"_"//AddSlash(trim(adjustl(PWCBraLabel))))
  !      allocate(auxFileName,source=auxDir//DipoleLenFileRootName//"1."//AlphabeticNumber(mDip)//QCFileExtension)
  !      call OpenFile( auxFileName, uid, 'read', 'formatted' )
  !      call MatA.Read( uid )
  !      close( uid )
  !      call CheckBlock( MatA )
  !      call MatA.Transpose()
  !      call OpenFile( outFileName, uid, 'write', 'formatted' )
  !      call MatA.Write( uid )
  !      close( uid )
  !   end if

  ! end subroutine ClassCSFCCSpaceAcquireDipoleLen


  ! subroutine ClassCSFCCSpaceSetMultiplicity( Space, Multiplicity )
  !   class(ClassCSFCCSpace), intent(inout) :: Space
  !   integer                    , intent(in)    :: Multiplicity
  !   Space.Multiplicity = Multiplicity
  ! end subroutine ClassCSFCCSpaceSetMultiplicity


  ! subroutine ClassCSFCCSpaceSetStorageDir( Space, StorageDir )
  !   class(ClassCSFCCSpace), intent(inout) :: Space
  !   character(len=*)           , intent(in)    :: StorageDir
  !   allocate( Space.StorageDir, source = StorageDir )
  ! end subroutine ClassCSFCCSpaceSetStorageDir



  ! subroutine ClassCSFCCSpaceSetNuclearLabel( Space, Label )
  !   class(ClassCSFCCSpace), intent(inout) :: Space
  !   character(len=*)           , intent(in)    :: Label
  !   if( allocated(Space.NuclearLabel) ) deallocate( Space.NuclearLabel )
  !   allocate( Space.NuclearLabel, source = Label )
  ! end subroutine ClassCSFCCSpaceSetNuclearLabel
  

  ! function ClassCSFCCSpaceGetGroup( Space ) result(Group)
  !   class(ClassCSFCCSpace), intent(in) :: Space
  !   type(ClassGroup), pointer :: Group
  !   allocate(Group,source=Space.Group)
  ! end function ClassCSFCCSpaceGetGroup


  ! function ClassCSFCCSpaceGetLmax( Space ) result(Lmax)
  !   class(ClassCSFCCSpace), intent(in) :: Space
  !   integer :: Lmax
  !   Lmax=Space.LMax
  ! end function ClassCSFCCSpaceGetLmax

  ! function ClassCSFCCSpaceGetPionCharge( Space ) result(PionCharge)
  !   class(ClassCSFCCSpace), intent(in) :: Space
  !   integer :: PionCharge
  !   PionCharge=Space.ParentIonCharge
  ! end function ClassCSFCCSpaceGetPionCharge
    
  ! function ClassCSFCCSpaceGetStorageDir( Space ) result(StorageDir)
  !   class(ClassCSFCCSpace), intent(in) :: Space
  !   character(len=:), allocatable :: StorageDir
  !   allocate(StorageDir,source=Space.StorageDir)
  ! end function ClassCSFCCSpaceGetStorageDir



  ! !> Retrieves wether the electronic space has been initialized or not.
  ! logical function ClassCSFCCSpaceInitialized( Self ) result(res)
  !   !> Class of the electronic space.
  !   class(ClassCSFCCSpace), intent(in) :: Self
  !   !
  !   res = .false.
  !   !
  !   if ( allocated(Self.IrrepSpaceVec) ) res = .true.
  !   !
  ! end function ClassCSFCCSpaceInitialized





  ! !> Gets the number of irreducible representations of the point group used in the close-coupling expansion.
  ! integer function ClassCSFCCSpaceGetNIrreps( Self ) result(N)
  !   !> Class of the electronic space.
  !   class(ClassCSFCCSpace), intent(in) :: Self
  !   !
  !   if ( .not.Self.Initialized() ) call Assert( "Impossible to get the number of irreducible representations of the "//&
  !        "electronic space because hast not been initialized." )
  !   N = Self.NIrreps
  !   !
  ! end function ClassCSFCCSpaceGetNIrreps


  ! !> Gets the electronic space group name.
  ! function ClassCSFCCSpaceGetGroupName( Self ) result(GName)
  !   !> Class of the electronic space.
  !   class(ClassCSFCCSpace), intent(in) :: Self
  !   character(len=:), allocatable :: GName
  !   !
  !   if ( .not.Self.Initialized() ) call Assert( "Impossible to get the group name of the "//&
  !        "electronic space because hast not been initialized." )
  !   allocate( GName, source = Self.Group.GetName() )
  !   !
  ! end function ClassCSFCCSpaceGetGroupName


  ! !> Get from the electronic space, the symmetric electronic space giving the corresponding index in the list.
  ! function ClassCSFCCSpaceGetSymElectSpace( Self, iIrrep ) result( SymSpace )
  !   !> Class of the electronic space.
  !   class(ClassCSFCCSpace), target, intent(in) :: Self
  !   !> Index correspondind to the symmetric electronic space in the list.
  !   integer                    , intent(in) :: iIrrep
  !   !> Class of the symmetric electronic space.
  !   type(ClassSymmetricElectronicSpace), pointer :: SymSpace
  !   !
  !   if ( .not.Self.Initialized() ) call Assert( "To get a symmetric electronic space, this must be initialized." )
  !   !
  !   allocate( SymSpace, source = Self.IrrepSpaceVec(iIrrep) )
  !   !
  ! end function ClassCSFCCSpaceGetSymElectSpace



  ! !> Get from the electronic space, the symmetric electronic space giving the irreducible representation name.
  ! function ClassCSFCCSpaceGetSymElectSpaceByIrrepName( Self, IrrepName ) result( SymSpace )
  !   !> Class of the electronic space.
  !   class(ClassCSFCCSpace), target, intent(in) :: Self
  !   !> Index correspondind to the symmetric electronic space in the list.
  !   character(len=*)                   , intent(in) :: IrrepName
  !   !> Class of the symmetric electronic space.
  !   type(ClassSymmetricElectronicSpace), pointer :: SymSpace
  !   integer :: i
  !   !
  !   if ( .not.Self.Initialized() ) call Assert( "To get a symmetric electronic space, this must be initialized." )
  !   !
  !   do i = 1, Self.GetNIrreps()
  !      if ( Self.IrrepSpaceVec(i).GetIrrepLabel() .is. IrrepName ) then
  !         allocate( SymSpace, source = Self.IrrepSpaceVec(i) )
  !         return
  !      end if
  !   end do
  !   !
  ! end function ClassCSFCCSpaceGetSymElectSpaceByIrrepName


  ! subroutine ClassCSFCCSpaceCheckSymmetry( Space, SymLabel )
  !   class(ClassCSFCCSpace), target, intent(in) :: Space
  !   character(len=*)                   , intent(in) :: SymLabel
  !   !
  !   integer :: NumSymElectSpace, i
  !   type(ClassSymmetricElectronicSpace), pointer :: SymElectSpace
  !   logical :: PresentIrrep
  !   !
  !   NumSymElectSpace = Space.GetNIrreps()
  !   PresentIrrep = .false.
  !   do i = 1, NumSymElectSpace
  !      allocate( SymElectSpace, source = Space.GetSymElectSpace(i) )
  !      if ( SymLabel .is. SymElectSpace.GetIrrepLabel() ) then
  !         PresentIrrep = .true.
  !         return
  !      end if
  !      deallocate( SymElectSpace )
  !   end do
  !   if ( .not.PresentIrrep ) call Assert( &
  !        'The irreducible representation '//SymLabel//&
  !        ' is not present in the close-coupling expansion for the '//Space.GetGroupName()//&
  !        ' group.' )
  ! end subroutine ClassCSFCCSpaceCheckSymmetry

  

  ! function GetCloseCouplingDir() result(Dir)
  !   character(len=:), allocatable :: Dir
  !   allocate( Dir, source = CLOSE_COUPLING_DIR )
  ! end function GetCloseCouplingDir


  ! function GetQCFileExtension() result(Dir)
  !   character(len=:), allocatable :: Dir
  !   allocate( Dir, source = QCFileExtension )
  ! end function GetQCFileExtension


  ! !> By convention, if a block stored has a 0 dimension, this is
  ! !! changed by 1, anf fill with zeros.
  ! subroutine CheckReadBlock( Array, ArrayIsFine )
  !   real(kind(1d0)), allocatable, intent(inout) :: Array(:,:)
  !   logical, optional           , intent(out)   :: ArrayIsFine
  !   !
  !   integer :: NRows, NCols
  !   logical :: BlockIsZero
  !   integer, parameter :: NewDimension = 1
  !   !
  !   BlockIsZero = .false.
  !   NRows = size(Array,1)
  !   NCols = size(Array,2)
  !   !
  !   if ( NRows == 0 ) then
  !      NRows = NewDimension
  !      BlockIsZero = .true.
  !   end if
  !   !
  !   if ( NCols == 0 ) then
  !      NCols = NewDimension
  !      BlockIsZero = .true.
  !   end if
  !   !
  !   if ( BlockIsZero ) then
  !      deallocate( Array )
  !      allocate( Array(NRows,NCols) )
  !      Array = 0.d0
  !   end if
  !   !
  !   if ( present(ArrayIsFine) ) then
  !      if ( BlockIsZero ) then
  !         ArrayIsFine = .false.
  !      else
  !         ArrayIsFine = .true.
  !      end if
  !   end if
  !   !
  ! end subroutine CheckReadBlock

  
  ! !> By convention, if a block stored has a 0 dimension, this is
  ! !! changed by 1, anf fill with zeros.
  ! subroutine CheckBlock( Mat )
  !   class(ClassMatrix), intent(inout) :: Mat
  !   !
  !   integer :: NRows, NCols
  !   logical :: BlockIsZero
  !   integer, parameter :: NewDimension = 1
  !   !
  !   BlockIsZero = .false.
  !   NRows = Mat.NRows()
  !   NCols = Mat.NColumns()
  !   !
  !   if ( NRows == 0 ) then
  !      NRows = NewDimension
  !      BlockIsZero = .true.
  !   end if
  !   !
  !   if ( NCols == 0 ) then
  !      NCols = NewDimension
  !      BlockIsZero = .true.
  !   end if
  !   !
  !   if ( BlockIsZero ) then
  !      call Mat.Free()
  !      call Mat.InitFull( NRows, NCols )
  !   end if
  !   !
  ! end subroutine CheckBlock



  ! subroutine CheckWrittenBlock( FileName )
  !   character(len=*), intent(in) :: FileName
  !   real(kind(1d0)), allocatable :: Array(:,:)
  !   logical :: ArrayIsFine
  !   integer :: uid
  !   type(ClassMatrix) :: Mat
  !   call ReadMatrix( FileName , Array )
  !   call CheckReadBlock( Array )
  !   Mat = Array
  !   deallocate( Array )
  !   call OpenFile( FileName, uid, 'write', 'formatted' )
  !   call Mat.Write( uid )
  !   close( uid )
  ! end subroutine CheckWrittenBlock





end module ModuleCSFCCSpace
