module ModuleCSFCCSymSpace
  
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleIO
  use ModuleParentIons
  
  
  implicit none

  private

  complex(kind(1d0)), parameter :: Z0 = (0.d0,0.d0)
  character(len=*)  , parameter :: AUG_SUBDIR = "aug/"
  character(len=*)  , parameter :: ION_SUBDIR = "ion/"
  character(len=*)  , parameter :: ONE_SUBDIR = "one/"
  integer           , parameter :: MAX_NLW_PER_PION_TERM = 32
  
  type, public :: ClassCSFCCSymSpace
     
     private

     character(len=:), allocatable :: CSF_Basis_File
     character(len=:), allocatable :: StoreDir
     character(len=3)              :: term
     integer                       :: multiplicity

     integer                       :: nCSFPWC !.. Total number of CSF PWC
     integer                       :: nCFG    !.. Number of CSF: ncfg = nref + nCSFPWC
     integer                       :: nREF    !.. Among the configs, how many are LC states
     integer         , allocatable :: nCSFPWC_per_CSF_pion(:)
     integer         , allocatable :: CSFPWC_ipion(:)
     integer         , allocatable :: CSFPWC_dim(:)
     integer         , allocatable :: CSFPWC_idim(:)
     integer         , allocatable :: CSFPWC_lw(:)
     integer                       :: CSFPWC_totdim

     integer                       :: tot_ncsf_pions !auxiliary
     integer         , allocatable :: nlw_per_pion_term(:)
     integer         , allocatable :: lwvec_per_pion_term(:,:)
     
   contains
     
       generic, public :: Load => ClassCSFCCSymSpaceLoadBasis
       generic, public :: Show => ClassCSFCCSymSpaceShow
       generic, public :: Get_lw_ppit   => ClassCSFCCSymSpaceGet_lw_ppit
       generic, public :: GetnRef => ClassCSFCCSymSpaceGetnRef
       generic, public :: Get_iTermOfPWCpion => ClassCSFCCSymSpaceGet_iTermOfPWCpion
       generic, public :: Get_lwOfPWC => ClassCSFCCSymSpaceGet_lwOfPWC
       generic, public :: Get_PWCdim => ClassCSFCCSymSpaceGet_PWCdim
       generic, public :: Get_PWCidim => ClassCSFCCSymSpaceGet_PWCidim
       generic, public :: Get_nPWC => ClassCSFCCSymSpaceGet_nPWC
       generic, public :: Get_TotalSize => ClassCSFCCSymSpaceGet_TotalSize
       generic, public :: LoadH0 => ClassCSFCCSymSpaceLoadH0
       generic, public :: LoadH0Row => ClassCSFCCSymSpaceLoadH0Row
       generic, public :: LoadBlockStripe => ClassCSFCCSymSpaceLoadBlockStripe
       generic, public :: SaveH0 => ClassCSFCCSymSpaceSaveH0
       generic, public :: GetH0Dir => ClassCSFCCSymSpaceGetH0Dir
       
       !.. Private procedures
       generic, private :: InitNames => ClassCSFCCSymSpaceInitNames
       generic, private :: GetAbsoluteIonIndexOfPWC => ClassCSFCCSymSpaceGetAbsoluteIonIndexOfPWC

       procedure, private :: ClassCSFCCSymSpaceGet_nPWC
       procedure, private :: ClassCSFCCSymSpaceGet_PWCdim
       procedure, private :: ClassCSFCCSymSpaceGet_PWCidim
       procedure, private :: ClassCSFCCSymSpaceGet_lwOfPWC
       procedure, private :: ClassCSFCCSymSpaceGet_iTermOfPWCpion
       procedure, private :: ClassCSFCCSymSpaceGetAbsoluteIonIndexOfPWC
       procedure, private :: ClassCSFCCSymSpaceLoadBasis
       procedure, private :: ClassCSFCCSymSpaceLoadH0
       procedure, private :: ClassCSFCCSymSpaceLoadH0Row
       procedure, private :: ClassCSFCCSymSpaceLoadBlockStripe
       procedure, private :: ClassCSFCCSymSpaceSaveH0
       procedure, private :: ClassCSFCCSymSpaceGetH0Dir
       procedure, private :: ClassCSFCCSymSpaceShow
       procedure, private :: ClassCSFCCSymSpaceInitNames
       procedure, private :: ClassCSFCCSymSpaceGet_lw_ppit
       procedure, private :: ClassCSFCCSymSpaceGetnRef
       procedure, private :: ClassCSFCCSymSpaceGet_TotalSize

  end type ClassCSFCCSymSpace

contains

  integer function ClassCSFCCSymSpaceGet_TotalSize( Self ) result( n )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    n = self.CSFPWC_totdim
  end function ClassCSFCCSymSpaceGet_TotalSize

  integer function ClassCSFCCSymSpaceGet_nPWC( Self ) result( n )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    n = self.nCSFPWC
  end function ClassCSFCCSymSpaceGet_nPWC

  integer function ClassCSFCCSymSpaceGet_PWCdim( Self, iPWC ) result( n )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    integer                  , intent(in) :: iPWC
    n = self.CSFPWC_dim( iPWC )
  end function ClassCSFCCSymSpaceGet_PWCdim

  integer function ClassCSFCCSymSpaceGet_PWCidim( Self, iPWC ) result( n )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    integer                  , intent(in) :: iPWC
    n = self.CSFPWC_idim( iPWC )
  end function ClassCSFCCSymSpaceGet_PWCidim

  integer function ClassCSFCCSymSpaceGet_lwOfPWC( Self, iPWC ) result( lw )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    integer                  , intent(in) :: iPWC
    lw = self.CSFPWC_lw( iPWC )
  end function ClassCSFCCSymSpaceGet_lwOfPWC

  integer function ClassCSFCCSymSpaceGet_iTermOfPWCpion( Self, iPWC ) result( iTerm )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    integer                  , intent(in) :: iPWC
    integer :: iPion
    iPion = self.CSFPWC_ipion( iPWC )
    iTerm = ParentIons.Get_iTerm( iPion )
  end function ClassCSFCCSymSpaceGet_iTermOfPWCpion
    
  integer function ClassCSFCCSymSpaceGetAbsoluteIonIndexOfPWC( Self, iPWC ) result( iPion )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    integer                  , intent(in) :: iPWC
    iPion = self.CSFPWC_ipion( iPWC )
  end function ClassCSFCCSymSpaceGetAbsoluteIonIndexOfPWC
    
  integer function ClassCSFCCSymSpaceGetnRef( Self ) result( n )
    class(ClassCSFCCSymSpace), intent(in)  :: Self
    n = self.nref
  end function ClassCSFCCSymSpaceGetnRef
    
  subroutine ClassCSFCCSymSpaceGet_lw_ppit( Self, vec, mat )
    class(ClassCSFCCSymSpace), intent(in)  :: Self
    integer, allocatable     , intent(out) :: vec(:), mat(:,:)
    allocate(vec,source=self.nlw_per_pion_term)
    allocate(mat,source=self.lwvec_per_pion_term)
  end subroutine ClassCSFCCSymSpaceGet_lw_ppit
    

  !.. Return the directory where the hamiltonian is stored
  subroutine ClassCSFCCSymSpaceGetH0Dir( Self, subdir, H0Dir )
    class(ClassCSFCCSymSpace)    , intent(inout) :: Self
    character(len=*)             , intent(in)    :: subdir
    character(len=:), allocatable, intent(out)   :: H0Dir
    allocate( H0Dir, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
         "H0_CSF_"//subdir//"/" )
  end subroutine ClassCSFCCSymSpaceGetH0Dir


  !.. Fetch the field-free hamiltonian from disk
  subroutine ClassCSFCCSymSpaceLoadH0Row( Self, subdir, zH0, iPWC )

    class(ClassCSFCCSymSpace)      , intent(inout) :: Self
    character(len=*)               , intent(in)    :: subdir
    complex(kind(1d0)), allocatable, intent(out)   :: zH0(:,:)
    integer                        , intent(in)    :: iPWC ! 0 => LC, 1,2,...,nCSFPWC => CSFPWCs
    
    character(len=:), allocatable   :: FileName, pwcLabel, ipwcLabel, jpwcLabel
    integer                         :: i, j, n, nlc
    integer                         :: di
    integer                         :: jpwc, dj, jOffset
    complex(kind(1d0)), allocatable :: zmat(:,:)
    

    if( iPWC < 0 .or. iPWC > Self.nCSFPWC )then
       call Assert("Invalid channel index in the call to LoadH0Row")
       stop
    endif

    n = Self.Get_TotalSize()
    nlc = self.nref
    if( iPWC == 0 )then
       
       if(allocated(zH0))deallocate(zH0)
       allocate(zH0(nlc,n))
       zH0=z0

       !.. Load LC-LC
       allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
            "H0_CSF_"//subdir//"/H0_LC_LC" )
       call LoadMatrix( FileName, zmat, "unformatted" )
       i = size( zmat, 1 )
       j = size( zmat, 2 )
       if( i /= nlc .or. j /= nlc )then
          write(ERROR_UNIT,*) "*** nlc = ",nlc,", but on file ("//FileName//") nref = ",i,j
          call Assert("Inconsistent LC size ")
       end if
       zH0(1:nlc,1:nlc)=zmat
       deallocate(FileName)

       !.. Get LC-PWC by loading the PWC-LC block and transposing those
       do jpwc = 1, self.nCSFPWC
          jOffset= self.CSFPWC_idim(jpwc-1)
          dj = self.CSFPWC_dim( jPWC )
          pwcLabel = AlphabeticNumber( jPWC, self.nCSFPWC, "0" )
          allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
               "H0_CSF_"//subdir//"/H0_"//pwcLabel//"_LC" )
          deallocate(pwcLabel)
          call LoadMatrix( FileName, zmat, "unformatted" )
          i = size( zmat, 1 )
          j = size( zmat, 2 )
          if( i /= dj .or. j /= nlc )then
             write(ERROR_UNIT,*) "*** dj, nlc = ",dj, nlc,", but on file ("//FileName//") dj,nlc = ",i,j
             call Assert("Inconsistent size ")
          end if
          zH0( 1:nlc, jOffset+1:jOffset+dj ) = transpose( zmat )
          deallocate(FileName)
       enddo
       
    else

       di = self.CSFPWC_dim( iPWC )
       ipwcLabel = AlphabeticNumber( iPWC, self.nCSFPWC, "0" )

       if(allocated(zH0))deallocate(zH0)
       allocate(zH0(di,n))
       zH0=z0
       
       !.. Load the PWC - LC block 
       allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
            "H0_CSF_"//subdir//"/H0_"//ipwcLabel//"_LC" )
       call LoadMatrix( FileName, zmat, "unformatted" )
       i = size( zmat, 1 )
       j = size( zmat, 2 )
       if( i /= di .or. j /= nlc )then
          write(ERROR_UNIT,*) "*** di, nlc = ",di, nlc,", but on file ("//FileName//") di,nlc = ",i,j
          call Assert("Inconsistent size ")
       end if
       zH0(1:di,1:nlc)=zmat
       deallocate(FileName)

       !.. Load the iPWC-jPWC blocks with jPWC <= iPWC 
       do jpwc = 1, ipwc
          !
          jOffset= self.CSFPWC_idim(jpwc-1)
          dj = self.CSFPWC_dim( jPWC )
          if(allocated(jpwcLabel)) deallocate(jpwcLabel)
          jpwcLabel = AlphabeticNumber( jPWC, self.nCSFPWC, "0" )
          !
          if(allocated(FileName)) deallocate(FileName)
          allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
               "H0_CSF_"//subdir//"/H0_"//ipwcLabel//"_"//jpwcLabel )
          !
          call LoadMatrix( FileName, zmat, "unformatted" )
          i = size( zmat, 1 )
          j = size( zmat, 2 )
          if( i /= di .or. j /= dj )then
             write(ERROR_UNIT,*) "*** di, dj = ",di, dj,", but on file ("//FileName//") di,dj = ",i,j
             call Assert("Inconsistent size ")
          end if
          zH0( 1:di, jOffset + 1 : jOffset + dj ) = zmat
       enddo

       !.. Get the iPWC-jPWC blocks with jPWC > iPWC by loading the jPWC-iPWC block and transposing it
       do jpwc = ipwc + 1, Self.nCSFPWC
          !
          jOffset= self.CSFPWC_idim(jpwc-1)
          dj = self.CSFPWC_dim( jPWC )
          if(allocated(jpwcLabel)) deallocate(jpwcLabel)
          jpwcLabel = AlphabeticNumber( jPWC, self.nCSFPWC, "0" )
          !
          if(allocated(FileName)) deallocate(FileName)
          allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
               "H0_CSF_"//subdir//"/H0_"//jpwcLabel//"_"//ipwcLabel )
          !
          call LoadMatrix( FileName, zmat, "unformatted" )
          i = size( zmat, 2 )
          j = size( zmat, 1 )
          if( i /= di .or. j /= dj )then
             write(ERROR_UNIT,*) "*** di, dj = ",di, dj,", but on file ("//FileName//") di,dj = ",i,j
             call Assert("Inconsistent size ")
          end if
          zH0( 1:di, jOffset + 1 : jOffset + dj ) = transpose(zmat)
          !
       enddo

    endif

    return
  end subroutine ClassCSFCCSymSpaceLoadH0Row



  !.. Fetch the field-free hamiltonian from disk
  subroutine ClassCSFCCSymSpaceLoadBlockStripe( Self, dim, dir, prefix, suffix, zA )

    class(ClassCSFCCSymSpace)      , intent(inout) :: Self
    integer                        , intent(in)    :: dim
    character(len=*)               , intent(in)    :: dir
    character(len=*)               , intent(in)    :: prefix
    character(len=*)               , intent(in)    :: suffix
    complex(kind(1d0)), allocatable, intent(out)   :: zA(:,:)
    
    integer                        , parameter     :: BRA_CSF_STRIPE = 1
    integer                        , parameter     :: KET_CSF_STRIPE = 2

    character(len=:), allocatable   :: FileName, pwcLabel
    integer                         :: ipwc, di, iOffset, n, n_var
    complex(kind(1d0)), allocatable :: zmat(:,:)
    
    n = Self.Get_TotalSize()
    if(dim /= 1 .and. dim /= 2 ) call Assert("Invalid value of dim in LoadBlockStripe")
    do ipwc = 0, self.nCSFPWC

       pwcLabel = AlphabeticNumber( iPWC, self.nCSFPWC, "0" )
       if(allocated(FileName))deallocate(FileName)
       allocate( FileName, source = dir // prefix // pwclabel // suffix )
       call LoadMatrix( FileName, zmat, "unformatted" )
       
       select case ( dim )
       case( BRA_CSF_STRIPE )

          if( ipwc == 0 )then
             n_var = size( zmat, 2 )
             allocate( zA( n, n_var ) )
             zA = Z0
          endif
          if( n_var /= size( zmat, 2 ) ) call Assert("Inconsistent size in "//FileName)

          iOffset = self.CSFPWC_idim(ipwc-1)
          di      = self.CSFPWC_dim( ipwc  )
          zA( iOffset + 1 : iOffset + di, : ) = zmat
          
       case( KET_CSF_STRIPE )

          if( ipwc == 0 )then
             n_var = size( zmat, 1 )
             allocate( zA( n_var, n ) )
             zA = Z0
          endif
          if( n_var /= size( zmat, 1 ) ) call Assert("Inconsistent size in "//FileName)

          iOffset = self.CSFPWC_idim(ipwc-1)
          di      = self.CSFPWC_dim( ipwc  )
          zA( :, iOffset + 1 : iOffset + di ) = zmat
          
       end select

    enddo

    return
  end subroutine ClassCSFCCSymSpaceLoadBlockStripe



  !.. Fetch the field-free hamiltonian from disk
  subroutine ClassCSFCCSymSpaceLoadH0( Self, subdir, zH0 )

    class(ClassCSFCCSymSpace)      , intent(inout) :: Self
    character(len=*)               , intent(in)    :: subdir
    complex(kind(1d0)), allocatable, intent(out)   :: zH0(:,:)
    
    character(len=:), allocatable   :: FileName, pwcLabel, ipwcLabel, jpwcLabel
    integer                         :: i, j, n, nlc
    integer                         :: ipwc, di, iOffset
    integer                         :: jpwc, dj, jOffset
    complex(kind(1d0)), allocatable :: zmat(:,:)
    

    if(allocated(zH0))deallocate(zH0)
    n = Self.Get_TotalSize()
    allocate(zH0(n,n))
    zH0=z0
    
    !.. Load LC-LC
    nlc = self.nref
    allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
         "H0_CSF_"//subdir//"/H0_LC_LC" )
    
    call LoadMatrix( FileName, zmat, "unformatted" )
    i = size( zmat, 1 )
    j = size( zmat, 2 )
    if( i /= nlc .or. j /= nlc )then
       write(ERROR_UNIT,*) "*** nlc = ",nlc,", but on file ("//FileName//") nref = ",i,j
       call Assert("Inconsistent LC size ")
    end if
    zH0(1:nlc,1:nlc)=zmat
    deallocate(FileName)


    !.. Load PWC-LC
    do ipwc = 1, self.nCSFPWC
       iOffset= self.CSFPWC_idim(ipwc-1)
       di = self.CSFPWC_dim( iPWC )
       pwcLabel = AlphabeticNumber( iPWC, self.nCSFPWC, "0" )
       allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
            "H0_CSF_"//subdir//"/H0_"//pwcLabel//"_LC" )
       deallocate(pwcLabel)
       call LoadMatrix( FileName, zmat, "unformatted" )
       i = size( zmat, 1 )
       j = size( zmat, 2 )
       if( i /= di .or. j /= nlc )then
          write(ERROR_UNIT,*) "*** di, nlc = ",di, nlc,", but on file ("//FileName//") di,nlc = ",i,j
          call Assert("Inconsistent size ")
       end if
       zH0(iOffset+1:iOffset+di,1:nlc)=zmat
       deallocate(FileName)
    enddo


    !.. Load PWC-PWC
    do ipwc = 1, self.nCSFPWC
       !
       iOffset= self.CSFPWC_idim(ipwc-1)
       di = self.CSFPWC_dim( iPWC )
       if(allocated(ipwcLabel)) deallocate(ipwcLabel)
       ipwcLabel = AlphabeticNumber( iPWC, self.nCSFPWC, "0" )
       !
       do jpwc = 1, ipwc
          !
          jOffset= self.CSFPWC_idim(jpwc-1)
          dj = self.CSFPWC_dim( jPWC )
          if(allocated(jpwcLabel)) deallocate(jpwcLabel)
          jpwcLabel = AlphabeticNumber( jPWC, self.nCSFPWC, "0" )
          !
          if(allocated(FileName)) deallocate(FileName)
          allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
               "H0_CSF_"//subdir//"/H0_"//ipwcLabel//"_"//jpwcLabel )
          !
          call LoadMatrix( FileName, zmat, "unformatted" )
          i = size( zmat, 1 )
          j = size( zmat, 2 )
          if( i /= di .or. j /= dj )then
             write(ERROR_UNIT,*) "*** di, dj = ",di, dj,", but on file ("//FileName//") di,dj = ",i,j
             call Assert("Inconsistent size ")
          end if
          zH0( iOffset + 1 : iOffset + di, jOffset + 1 : jOffset + dj ) = zmat
       enddo
       
    enddo

    !.. Symmetrizes the hamiltonian
    do i=1,n-1
       do j=i+1,n
          zH0(i,j) = zH0(j,i)
       enddo
    enddo

    return
  end subroutine ClassCSFCCSymSpaceLoadH0
    

  !.. Fetch the field-free hamiltonian from disk
  subroutine ClassCSFCCSymSpaceSaveH0( Self, subdir, zH0 )

    class(ClassCSFCCSymSpace)      , intent(inout) :: Self
    character(len=*)               , intent(in)    :: subdir
    complex(kind(1d0)), allocatable, intent(in)    :: zH0(:,:)
    
    character(len=:), allocatable   :: FileName, pwcLabel, ipwcLabel, jpwcLabel
    integer                         :: n, nlc, d_max
    integer                         :: ipwc, di, iOffset
    integer                         :: jpwc, dj, jOffset
    complex(kind(1d0)), allocatable :: zmat(:,:)
    

    n = Self.Get_TotalSize()
    
    !.. Save LC-LC
    nlc = self.nref
    allocate(zmat(nlc,nlc))
    zmat = zH0(1:nlc,1:nlc)
    allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
         "H0_CSF_"//subdir//"/H0_LC_LC" )
    call SaveMatrix( FileName, zmat, "unformatted" )
    deallocate(FileName)
    deallocate(zmat)

    !.. Find max channel size
    d_max=0
    do ipwc = 1, self.nCSFPWC
       d_max=max(d_max,self.CSFPWC_dim( iPWC ))
    enddo


    !.. Save PWC-LC
    allocate(zmat(d_max,nlc))
    zmat=Z0
    !
    do ipwc = 1, self.nCSFPWC
       !
       iOffset= self.CSFPWC_idim(ipwc-1)
       di = self.CSFPWC_dim( iPWC )
       !
       zmat(1:di,1:nlc) = zH0(iOffset+1:iOffset+di,1:nlc)
       !
       pwcLabel = AlphabeticNumber( iPWC, self.nCSFPWC, "0" )
       allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
            "H0_CSF_"//subdir//"/H0_"//pwcLabel//"_LC" )
       !
       call SaveMatrix( FileName, zmat, di, nlc, "unformatted" )
       !
       deallocate(pwcLabel)
       deallocate(FileName)
       !
    enddo
    deallocate(zmat)


    !.. Save PWC-PWC
    
    !.. Allocate utility block to the max channel size
    allocate(zmat(d_max,d_max))
    zmat=Z0
    !
    do ipwc = 1, self.nCSFPWC
       !
       iOffset= self.CSFPWC_idim(ipwc-1)
       di = self.CSFPWC_dim( iPWC )
       if(allocated(ipwcLabel)) deallocate(ipwcLabel)
       ipwcLabel = AlphabeticNumber( iPWC, self.nCSFPWC, "0" )
       !
       do jpwc = 1, ipwc
          !
          jOffset= self.CSFPWC_idim(jpwc-1)
          dj = self.CSFPWC_dim( jPWC )
          if(allocated(jpwcLabel)) deallocate(jpwcLabel)
          jpwcLabel = AlphabeticNumber( jPWC, self.nCSFPWC, "0" )
          !
          zmat(1:di,1:dj) = zH0( iOffset + 1 : iOffset + di, jOffset + 1 : jOffset + dj )
          !
          if(allocated(FileName)) deallocate(FileName)
          allocate( FileName, source = Self.StoreDir // AUG_SUBDIR // Self.term // "/" // &
               "H0_CSF_"//subdir//"/H0_"//ipwcLabel//"_"//jpwcLabel )
          !
          call SaveMatrix( FileName, zmat, di, dj, "unformatted" )
          !
          !
       enddo
       
    enddo
    deallocate(zmat)

    return
  end subroutine ClassCSFCCSymSpaceSaveH0
    

  subroutine ClassCSFCCSymSpaceLoadBasis( Self, store, term )

    class(ClassCSFCCSymSpace), intent(inout) :: Self
    character(len=*)         , intent(in) :: store
    character(len=*)         , intent(in) :: term

    integer              :: uid, ipion, iCSFPWC
    character(len=128)   :: strnBuf
    character(len=10000) :: longLine
    integer :: nBuf, iBuf, ilw
    integer :: nterms_pion !auxiliary
    integer :: iterm_pion

    call Self.InitNames( store, term )

    call OpenFile( Self.CSF_Basis_File, uid, "read", "formatted" )

    !.. Reads the total number of configurations, ncfg, 
    !   and the number of "reference" states, nref
    !   (i.e., the number of localized configs)
    !..
    read(uid,*) strnBuf, Self.ncfg, Self.nref
    read(uid,"(a)") longLine
    Self.tot_ncsf_pions = nTokens(trim(longLine))
    nBuf = ParentIons.GetTotNcsf()
    if( Self.tot_ncsf_pions /= nBuf )then
       call ErrorMessage("inconsistent total number of parent ions in "//Self.CSF_Basis_File)
       write(*,"(a,i0)") "ParentIons.GetTotNcsf ",nBuf
       write(*,"(a,i0)") "self.tot_ncsf_pions ",Self.tot_ncsf_pions
       write(*,"(a)")    "Line in basis file: "//trim(longLine)
       stop
    endif
    allocate( Self.nCSFPWC_per_CSF_pion( Self.tot_ncsf_pions ) )
    read(longLine,*) ( Self.nCSFPWC_per_CSF_pion( ipion ), ipion = 1, Self.tot_ncsf_pions )
    Self.nCSFPWC = sum( Self.nCSFPWC_per_CSF_pion )

    if( self.ncfg /= self.nref + self.nCSFPWC )then
       write(*,"(a,*(i0,a))") &
            "ncfg (",self.ncfg,")/= "//&
            "nref (",self.nref,") + "//&
            "nCSFPWC (",self.nCSFPWC,")"
       stop
    endif
    !
    !.. read CSFPWC_lw( iCSFPWC ): the angular momentum of the electron in channel iCSFPWC
    allocate( self.CSFPWC_lw  ( self.nCSFPWC ) )
    read(uid,*) ( self.CSFPWC_lw( iCSFPWC ), iCSFPWC = 1, self.nCSFPWC )
    !
    !.. read CSFPWC_dim( iCSFPWC ): the size of the iCSFPWC partial-wave channel
    allocate( self.CSFPWC_dim ( 0:self.nCSFPWC ) )
    self.CSFPWC_dim(0) = self.nref
    read(uid,*) ( self.CSFPWC_dim( iCSFPWC ), iCSFPWC = 1, self.nCSFPWC )
    close(uid)

    !.. Computes idim for CSF PWCs
    !..
    allocate( self.CSFPWC_idim( -1:self.nCSFPWC ) )
    self.CSFPWC_idim = 0
    do iCSFPWC = 0, self.nCSFPWC 
       self.CSFPWC_idim( iCSFPWC ) = self.CSFPWC_idim( iCSFPWC-1 ) + self.CSFPWC_dim( iCSFPWC ) 
    enddo
    self.CSFPWC_totdim = self.CSFPWC_idim( self.nCSFPWC )


    !.. Determines how many and which photoelectron angular momenta are coupled to each parent-ion term
    nterms_pion = ParentIons.GetNterms()
    allocate( self.nlw_per_pion_term( nterms_pion ) )
    allocate( self.lwvec_per_pion_term( MAX_NLW_PER_PION_TERM, nterms_pion ) )
    self.nlw_per_pion_term = 0
    self.lwvec_per_pion_term = 0
    do iterm_pion = 1, nterms_pion
       iPion = ParentIons.GetnCSFupToTerm( iterm_pion - 1 ) + 1 
       self.nlw_per_pion_term( iterm_pion ) = self.nCSFPWC_per_CSF_pion( iPion )
       iBuf = sum( self.nCSFPWC_per_CSF_pion( 1:iPion ) ) - self.nCSFPWC_per_CSF_pion( iPion )
       self.lwvec_per_pion_term( 1:self.nlw_per_pion_term( iterm_pion ), iterm_pion ) = &
            self.CSFPWC_lw( iBuf+1 : iBuf+self.nlw_per_pion_term( iterm_pion ) )
       write(*,"(2(x,a),(x,i0),3x,*(i0,x))") term, ParentIons.GetTermName( iterm_pion ), &
            self.nlw_per_pion_term( iterm_pion ), &
            self.lwvec_per_pion_term( 1:self.nlw_per_pion_term( iterm_pion ), iterm_pion )
    enddo
    
    !.. Establishes the correspondence between partial-wave channels and parent ions
    !   CSFPWC_ipion( iCSFPWC ): the absolute index of the parent ion associated to iCSFPWC
    allocate( self.CSFPWC_ipion( self.nCSFPWC ) )
    iCSFPWC=0
    !
    !.. Cycle over all the parent ions
    do iPion = 1, self.tot_ncsf_pions
       !
       !.. Cycle over all the channels generated
       !   by the current parent ion
       do ilw = 1, self.nCSFPWC_per_CSF_pion( iPion )
          !
          iCSFPWC = iCSFPWC + 1
          self.CSFPWC_ipion( iCSFPWC ) = iPion
          !
       enddo
       !
    enddo
    !
  end subroutine ClassCSFCCSymSpaceLoadBasis


  subroutine ClassCSFCCSymSpaceShow( Self )
    class(ClassCSFCCSymSpace), intent(in) :: Self
    integer :: iPion, iTerm_Pion, iCSFPWC
    write(*,*)
    write(*,*) "Properties of the CSF Partial Wave Channels in "//self.Term//" symmetry"
    write(*,*) "======================================================================="
    write(*,"(a,i0)")" N CSF                :  ", self.nCFG
    write(*,"(a,i0)")" N CSF REF            :  ", self.nref
    write(*,"(a,i0)")" N CSF ions           :  ", self.tot_ncsf_pions
    write(*,"(a,i0)")" N CSF PWC            :  ", self.nCSFPWC
    write(*,"(a,*(x,i0))")" N CSF PWC per CSF ion:  ", &
         ( self.nCSFPWC_per_CSF_pion( ipion ), ipion = 1, self.tot_ncsf_pions )
    write(*,"(a,i0)")" Number of CSF REF:  ", self.CSFPWC_dim(0)
    write(*,"(a)") " CSFPWC  ion  Sym_ion  l   dim   idim"
    do iCSFPWC = 1, self.nCSFPWC
       iPion      = self.CSFPWC_ipion( iCSFPWC )
       iTerm_Pion = ParentIons.Get_iTerm( iPion )
       write(*,"(x,i4,2x,i4,3x,a,i0,a,i2,x,i5,x,i6)") iCSFPWC, iPion, &
            ParentIons.GetTermName( iTerm_Pion )//" (",iTerm_Pion,") ", &
            self.CSFPWC_lw( iCSFPWC ), self.CSFPWC_dim( iCSFPWC ), self.CSFPWC_idim( iCSFPWC )
    enddo
  end subroutine ClassCSFCCSymSpaceShow


  subroutine ClassCSFCCSymSpaceInitNames( Self, store, term )
    class(ClassCSFCCSymSpace), intent(inout) :: Self
    character(len=*)         , intent(in)    :: store
    character(len=*)         , intent(in)    :: term
    character(len=1024) :: strnBuf_Store, strnBuf_term 
    !
    Self.StoreDir = FormatAsDir( store )
    Self.term     = trim(adjustl(term))
    strnBuf_Store = trim(adjustl(Self.StoreDir))
    strnBuf_term  = trim(adjustl(term))
    allocate( Self.CSF_Basis_File, &
         source = &
         trim(strnBuf_Store) // &
         AUG_SUBDIR    // &
         trim(strnBuf_term) // "/"   // &
         "CSF_Basis_"//trim(strnBuf_term)//".dat" )
    !allocate( Self.CSF_Basis_File, &
    !     source = &
    !     Self.StoreDir // &
    !     AUG_SUBDIR    // &
    !     Self.term // "/"   // &
    !     "CSF_Basis_"//Self.term//".dat" )
    !
  end subroutine ClassCSFCCSymSpaceInitNames

end module ModuleCSFCCSymSpace
