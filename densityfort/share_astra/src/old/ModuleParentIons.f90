module ModuleParentIons

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleIO

  implicit none

  private

  character(len=*), parameter :: ION_DIR  = "ion/"
  character(len=*), parameter :: ALL_IONS_FILE = "pion_evecs_all.dat"

  type, public :: ClassParentIons
     !
     private
     !
     character(len=:), allocatable :: StoreDir
     character(len=:), allocatable :: Pion_Evec_File
     integer                       :: nterms_pion
     character(len=3), allocatable :: term_name_pion(:)
     integer         , allocatable :: nlw_per_pion_term(:)
     integer         , allocatable :: lwvec_per_pion_term(:,:)
     !
     integer         , allocatable :: term_nCSF_pion(:)
     integer         , allocatable :: term_inCSF_pion(:) 
     integer                       :: max_ncsf_pion
     integer                       :: tot_ncsf_pions
     !
     integer         , allocatable :: term_nCI_pion(:)
     integer                       :: max_nCI_pion
     !
     integer         , allocatable :: ipion_to_iterm_pion(:)
     real(kind(1d0)) , allocatable :: CI_Ener_pion(:,:)
     real(kind(1d0)) , allocatable :: CI_Coef_pion(:,:,:)
     !
   contains
     !
     generic  , public  :: Load => ClassParentIonsLoad
     generic  , public  :: Show => ClassParentIonsShow
     generic  , public  :: GetNterms => ClassParentIonsGetNterms
     generic  , public  :: GetTotNcsf => ClassParentIonsGetTotNcsf
     generic  , public  :: Get_nCIperTerm => ClassParentIonsGetnCIperTerm
     generic  , public  :: GetnCSFupToTerm => ClassParentIonsGetnCSFupToTerm
     generic  , public  :: GetTermName => ClassParentIonsGetTermName
     generic  , public  :: Get_iTerm => ClassParentIonsGet_iTerm
     generic  , public  :: GetiTermFromTerm => ClassParentIonsGetiTermFromTerm
     generic  , public  :: GetMaxNCSF => ClassParentIonsGetMaxNCSF
     generic  , public  :: Get_CIEnergy => ClassParentIonsGet_CIEnergy, ClassParentIonsGet_CIEnergyFromTerm
     generic  , public  :: Get_CIVector => ClassParentIonsGet_CIVector, ClassParentIonsGet_CIVectorFromTerm
     generic  , public  :: Get_nCSFperTerm => ClassParentIonsGet_nCSFperTerm, &
          ClassParentIonsGet_nCSFperTermFromName
     !
     generic  , private :: InitNames => ClassParentIonsInitNames
     !
     procedure, private :: ClassParentIonsLoad
     procedure, private :: ClassParentIonsShow
     procedure, private :: ClassParentIonsGetNterms
     procedure, private :: ClassParentIonsGetTotNcsf
     procedure, private :: ClassParentIonsGetnCSFupToTerm
     procedure, private :: ClassParentIonsGetnCIperTerm
     procedure, private :: ClassParentIonsGetTermName
     procedure, private :: ClassParentIonsGet_iTerm
     procedure, private :: ClassParentIonsGetiTermFromTerm
     procedure, private :: ClassParentIonsGetMaxNCSF
     procedure, private :: ClassParentIonsGet_CIEnergy
     procedure, private :: ClassParentIonsGet_CIEnergyFromTerm
     procedure, private :: ClassParentIonsGet_CIVector
     procedure, private :: ClassParentIonsGet_CIVectorFromTerm
     procedure, private :: ClassParentIonsGet_nCSFperTerm
     procedure, private :: ClassParentIonsGet_nCSFperTermFromName
     !
     procedure, private :: ClassParentIonsInitNames
     !
     ! final :: ClassParentIonsFinal
     !
  end type ClassParentIons


  !.. SINGLETON
  type(ClassParentIons), public :: ParentIons


contains

  integer function ClassParentIonsGetiTermFromTerm( self, Term ) result( iTerm )
    class(ClassParentIons), intent(in) :: self
    character(len=*)      , intent(in) :: Term
    integer :: it
    iTerm = 0
    do it = 1, self.GetNTerms()
       if( self.GetTermName( it ) == Term ) iTerm = it
    enddo
  end function ClassParentIonsGetiTermFromTerm

  integer function ClassParentIonsGet_nCSFperTermFromName( self, Term ) result( n )
    class(ClassParentIons), intent(in) :: self
    character(len=*)      , intent(in) :: Term
    integer :: iTerm
    iTerm = self.GetiTermFromTerm( Term )
    n = self.Get_nCSFperTerm( iTerm )
  end function ClassParentIonsGet_nCSFperTermFromName

  integer function ClassParentIonsGet_nCSFperTerm( self, iTerm ) result( n )
    class(ClassParentIons), intent(in) :: self
    integer               , intent(in) :: iTerm
    n = self.term_nCSF_pion( iTerm )
  end function ClassParentIonsGet_nCSFperTerm

  real(kind(1d0)) function ClassParentIonsGet_CIEnergy( self, iPion, iTerm ) result( E )
    class(ClassParentIons), intent(in) :: self
    integer               , intent(in) :: iPion
    integer               , intent(in) :: iTerm
    E = self.CI_Ener_pion( iPion, iTerm )
  end function ClassParentIonsGet_CIEnergy

  real(kind(1d0)) function ClassParentIonsGet_CIEnergyFromTerm( self, iPion, Term ) result( E )
    class(ClassParentIons), intent(in) :: self
    character(len=*)      , intent(in) :: Term
    integer               , intent(in) :: iPion
    integer :: iTerm
    iTerm=self.GetiTermFromTerm(Term)
    E = self.Get_CIEnergy( iPion, iTerm )
  end function ClassParentIonsGet_CIEnergyFromTerm

  subroutine ClassParentIonsGet_CIVectorFromTerm( self, iPion, Term, n, civec )
    class(ClassParentIons)      , intent(in)  :: self
    integer                     , intent(in)  :: iPion
    character(len=*)            , intent(in)  :: Term
    integer                     , intent(out) :: n
    real(kind(1d0)), allocatable, intent(out) :: civec(:)
    integer :: iTerm
    iTerm=self.GetiTermFromTerm(Term)
    call self.Get_CIVector( iPion, iTerm, n, civec )
  end subroutine ClassParentIonsGet_CIVectorFromTerm

  subroutine ClassParentIonsGet_CIVector( self, iPion, iTerm, n, civec )
    class(ClassParentIons)      , intent(in)  :: self
    integer                     , intent(in)  :: iPion, iTerm
    integer                     , intent(out) :: n
    real(kind(1d0)), allocatable, intent(out) :: civec(:)
    n = self.Get_nCSFperTerm( iTerm )
    if(allocated(civec))then
       if(size(civec,1)<n)then
          deallocate(civec)
          allocate(civec(n))
       else
          civec=0.d0
       endif
    else
       allocate(civec(n))
    endif
    civec(1:n) = self.CI_Coef_pion( 1:n, iPion, iTerm )
  end subroutine ClassParentIonsGet_CIVector

  integer function ClassParentIonsGetNCSFUpToTerm( self, iTerm ) result( n )
    class(ClassParentIons), intent(in) :: self
    integer               , intent(in) :: iTerm
    n = self.term_inCSF_pion( iTerm )
  end function ClassParentIonsGetNCSFUpToTerm

  integer function ClassParentIonsGetMaxNCSF( self ) result( n )
    class(ClassParentIons), intent(in) :: self
    n = self.max_ncsf_pion
  end function ClassParentIonsGetMaxNCSF

  integer function ClassParentIonsGetnCIperTerm( self, iTerm ) result( n )
    class(ClassParentIons), intent(in) :: self
    integer               , intent(in) :: iTerm
    n = self.term_nCI_pion( iTerm )
  end function ClassParentIonsGetnCIperTerm

  function ClassParentIonsGetTermName( self, iTerm ) result( name )
    class(ClassParentIons), intent(in) :: self
    integer               , intent(in) :: iTerm
    character(len=:), allocatable :: name
    allocate(name,source=trim(adjustl(self.term_name_pion(iTerm))))
  end function ClassParentIonsGetTermName

  integer function ClassParentIonsGet_iTerm( self, iPion ) result( iTerm )
    class(ClassParentIons), intent(in) :: self
    integer               , intent(in) :: iPion
    iTerm = self.ipion_to_iterm_pion( iPion )
  end function ClassParentIonsGet_iTerm

  integer function ClassParentIonsGetNterms( self ) result( n )
    class(ClassParentIons), intent(in) :: self
    n = self.nterms_pion
  end function ClassParentIonsGetNterms

  integer function ClassParentIonsGetTotNcsf( self ) result( n )
    class(ClassParentIons), intent(in) :: self
    n = self.tot_ncsf_pions
  end function ClassParentIonsGetTotNcsf

  subroutine ClassParentIonsLoad( self, store )
    !
    class(ClassParentIons), intent(inout) :: self
    character(len=*)      , intent(in)    :: store
    !
    integer :: uid, iterm, imi, ima
    integer :: nAvailableCIPions
    integer :: iCIPion, jCIpion, iCSFpion
    !
    call Self.InitNames( store )

    !.. Load the spectral resolution of the parent ions
    call OpenFile( Self.Pion_Evec_File, uid, "read", "formatted" )
    read(uid,*) Self.nterms_pion
    allocate( self.term_name_pion( self.nterms_pion ) )
    allocate( self.term_nCSF_pion( self.nterms_pion ) )
    allocate( self.term_nCI_pion ( self.nterms_pion ) )
    read(uid,*) ( self.term_name_pion( iterm ), iterm = 1, self.nterms_pion )
    read(uid,*) ( self.term_nCSF_pion( iterm ), iterm = 1, self.nterms_pion )
    self.max_ncsf_pion = maxval( self.term_nCSF_pion )
    self.tot_ncsf_pions= sum   ( self.term_nCSF_pion )
    !
    !.. term_inCSF_pion(iterm-1)+1 and idim_pion(iterm) are the 
    !   absolute position of the first and last parent ion, respectively,
    !   for the parent-ion term iterm.
    !..
    allocate( self.term_inCSF_pion ( 0: self.nterms_pion ) )
    self.term_inCSF_pion(0)=0
    do iterm = 1, self.nterms_pion
       self.term_inCSF_pion( iterm ) = &
            self.term_inCSF_pion( iterm - 1 ) + &
            self.term_nCSF_pion( iterm )
    enddo
    !
    !.. ipion_to_pion_iterm( ipion )
    !   gives the iterm corresponding to the absolute index 
    !   of a given parent ion
    !..
    allocate( self.ipion_to_iterm_pion( self.tot_ncsf_pions ) )
    do iterm = 1, self.nterms_pion
       imi = self.term_inCSF_pion( iterm - 1 ) + 1
       ima = self.term_inCSF_pion( iterm )
       self.ipion_to_iterm_pion( imi:ima ) = iterm
    enddo
    !
    !.. Read the number of CI parent ions that we wish to use to
    !   transform the hamiltonian
    !..
    read(uid,*) ( self.term_nCI_pion (iterm), iterm = 1, self.nterms_pion )
    self.max_nCI_pion = maxval( self.term_nCI_pion )
    !
    !.. 
    allocate( self.CI_Ener_pion( self.max_nCI_pion, self.nterms_pion ) )
    allocate( self.CI_Coef_pion( self.max_ncsf_pion, self.max_nCI_pion, self.nterms_pion ) )
    self.CI_Ener_pion = 0.d0
    self.CI_Coef_pion = 0.d0
    do iterm = 1, self.nterms_pion
       !.. Read how many CI parent ions are actually listed
       read(uid,*) nAvailableCIPions
       if(nAvailableCIPions < self.term_nCI_pion( iterm ))then
          write(ERROR_UNIT,*) "Inconsistency in the number of requested and available parent ions"
          stop
       endif
       !.. read CI energies for current term
       read(uid,*) ( &
            self.CI_Ener_pion( iCIpion, iterm ), &
            iCIpion = 1, self.term_nCI_pion( iterm ) &
            )
       !.. read requested CI eigenvectors (they are saved as rows)
       do jCIpion = 1, self.term_nCI_pion( iterm )
          read(uid,*) ( &
               self.CI_Coef_pion( iCSFpion, jCIpion, iterm ), &
               iCSFpion = 1, self.term_nCSF_pion( iterm ) &
               )
       end do
       !.. skip the eigenvectors that are available but not needed
       do jCIpion = self.term_nCI_pion( iterm ) + 1, nAvailableCIPions
          read(uid,*)
       enddo
    end do
    close(uid)

  end subroutine ClassParentIonsLoad

  subroutine ClassParentIonsInitNames( Self, store )
    class(ClassParentIons), intent(inout) :: Self
    character(len=*)      , intent(in)    :: store
    character(len=1024) :: strnBuf
    Self.StoreDir = FormatAsDir( store )
    strnBuf=trim(adjustl(Self.StoreDir))
    allocate( Self.Pion_Evec_File, &
         source = &
         trim(strnBuf) // &
         ION_DIR       // &
         ALL_IONS_FILE  )
    !allocate( Self.Pion_Evec_File, &
    !     source = &
    !    Self.StoreDir // &
    !    ION_DIR       // &
    !    ALL_IONS_FILE  )
  end subroutine ClassParentIonsInitNames


  subroutine ClassParentIonsShow( self, uid )
    class(ClassParentIons), intent(in) :: self
    integer, optional     , intent(in) :: uid
    !
    integer          :: outunit
    integer          :: iterm, iCIpion, jCIpion, iCSFpion
    !
    outunit=OUTPUT_UNIT
    if(present(uid))outunit=uid
    write(outunit,*) 
    write(outunit,"(a)") " Parent Ion info: "
    write(outunit,"(a)") " ===============  " 
    write(outunit,"(a,x,i0)")   "Number of ion symmetries :", self.nterms_pion
    write(outunit,"(a,*(x,a3))")"List of ion symmetries   :", &
         ( self.term_name_pion(iterm), iterm = 1, self.nterms_pion )
    write(outunit,"(a,*(x,i0))")"N of CSF per symmetry    :", &
         ( self.term_nCSF_pion(iterm), iterm = 1, self.nterms_pion )
    write(outunit,"(a,x,i0)")   "Max N of CSF per symmetry:", self.max_ncsf_pion
    write(outunit,"(a,*(x,i0))")"N of CI ions per sym.    :", &
         ( self.term_nCI_pion (iterm), iterm = 1, self.nterms_pion )
    do iterm = 1, self.nterms_pion
       write(outunit,*)
       write(outunit,"(2x,a)") "CI Energies "//self.term_name_pion(iterm)
       write(outunit,"(*(x,f12.6))") &
            ( self.CI_Ener_pion( iCIpion, iterm ), &
            iCIpion = 1, self.term_nCI_pion( iterm ) )
       write(outunit,"(2x,a)") "CI Eigenvectors (one per row)"
       do jCIpion = 1, self.term_nCI_pion( iterm )
          write(outunit,"(*(x,f10.6))") &
               ( self.CI_Coef_pion( iCSFpion, jCIpion, iterm ), &
               iCSFpion = 1, self.term_nCSF_pion( iterm ) )
       end do
    end do
    !
  end subroutine ClassParentIonsShow


end module ModuleParentIons
