module ModuleQCBasis


  use, intrinsic :: ISO_FORTRAN_ENV
  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO

  implicit none

  private

  character(len=*), parameter, public :: QC_NEW_LINE="\NL\ "

  integer          , private, parameter :: iATOM_LB = -1
  integer          , private, parameter :: NATOMS = 18
  character(len=11), private, parameter :: LIST_ATOMIC_NAMES(iATOM_LB:NATOMS)=(/&
       "BARE"    , "VIRTUAL"  , &
       "HYDROGEN", "HELIUM"   , &
       "LITHIUM" , "BERYLLIUM", "BORON"   , "CARBON" , "NITROGEN"   , "OXYGEN", "FLUORINE", "NEON" , &
       "SODIUM"  , "MAGNESIUM", "ALUMINUM", "SILICON", "PHOSPHOROUS", "SULFUR", "CHLORINE", "ARGON" /) 
  integer          , private, parameter :: LIST_ATOMIC_CHARGES(iATOM_LB:NATOMS)=(/&
       0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 /)

  private


  !.. Set of atoms 
  type, public :: ClassQCAtomicSetBasis

     logical :: INITIALIZED = .FALSE.
     integer :: nat
     type(ClassQCAtomBasis), pointer :: vAtom(:)

   contains

     procedure, public  :: free      => ClassQCAtomicSetBasis_Free
     procedure, public  :: ParseFile => ClassQCAtomicSetBasis_ParseFile
     procedure, public  :: show      => ClassQCAtomicSetBasis_Show
     !
     !.. Accessors
     procedure, public  :: GetNsh    => ClassQCAtomicSetBasis_GetNsh
     procedure, public  :: GetLmax   => ClassQCAtomicSetBasis_GetLmax
     procedure, public  :: GetCont  => ClassQCAtomicSetBasis_GetCont
     !
     !
     procedure, private :: ClassQCAtomicSetBasis_Free
     procedure, private :: ClassQCAtomicSetBasis_ParseFile
     procedure, private :: ClassQCAtomicSetBasis_Show
     procedure, private :: ClassQCAtomicSetBasis_GetNsh
     procedure, private :: ClassQCAtomicSetBasis_GetLmax
     procedure, private :: ClassQCAtomicSetBasis_GetCont
     
  end type ClassQCAtomicSetBasis


  !.. Defines an atom
  !..
  type, public :: ClassQCAtomBasis

     logical :: INITIALIZED = .FALSE.
     
     character(len=:), allocatable :: name
     integer                       :: charge
     integer                       :: number
     integer                       :: lmax !.. Max Angular momentum
     type(ClassQCShell), pointer   :: vShell(:)    

   contains

     procedure, public  :: free  => ClassQCAtomBasis_Free
     procedure, public  :: show  => ClassQCAtomBasis_Show
     procedure, public  :: Parse => ClassQCAtomBasis_Parse
     !
     !.. Accessors
     procedure, public  :: GetNsh => ClassQCAtomBasis_GetNsh
     procedure, public  :: GetCont => ClassQCAtomBasis_GetCont
     !
     !
     procedure, private :: ClassQCAtomBasis_Free
     procedure, private :: ClassQCAtomBasis_Show
     procedure, private :: ClassQCAtomBasis_Parse
     procedure, private :: ClassQCAtomBasis_GetNsh
     procedure, private :: ClassQCAtomBasis_GetCont
     
  end type ClassQCAtomBasis


  !.. Contains the exponents and contraction
  !   coefficients for a shell with well defined
  !   angular momentum
  !..
  type, public :: ClassQCShell

     logical :: INITIALIZED = .FALSE.
     
     integer :: l    !.. Angular momentum
     integer :: nExp !.. Number of exponents
     integer :: nCon !.. Number of contractions

     real(kind(1d0)), allocatable :: vexp(:)   !.. Vector of Exponents
     real(kind(1d0)), allocatable :: mcoe(:,:) !.. Matrix of contr. coef.

   contains

     procedure, public  :: init  => ClassQCShell_Init
     procedure, public  :: free  => ClassQCShell_Free
     procedure, public  :: show  => ClassQCShell_Show
     procedure, public  :: Parse => ClassQCShell_Parse
     !
     !.. Accessors
     procedure, public  :: GetNcon => ClassQCShell_GetNcon
     procedure, public  :: GetCont => ClassQCShell_GetCont
     !
     !
     procedure, private :: ClassQCShell_Init
     procedure, private :: ClassQCShell_Free
     procedure, private :: ClassQCShell_Show
     procedure, private :: ClassQCShell_Parse
     procedure, private :: ClassQCShell_GetNcon
     procedure, private :: ClassQCShell_GetCont
     
  end type ClassQCShell
  
contains

  !.. Initialization of a shell
  subroutine ClassQCShell_Init( self, l, nexp, ncon, vexp, mcoe )

    class( ClassQCShell ), intent(inout) :: self
    integer            , intent(in)    :: l
    integer            , intent(in)    :: nexp
    integer            , intent(in)    :: ncon
    real(kind(1d0))    , intent(in)    :: vexp(:)
    real(kind(1d0))    , intent(in)    :: mcoe(:,:)

    if(l < 0)     call Assert("Invalid number of exponents")
    if(nexp<=0  ) call Assert("Invalid number of exponents")
    if(ncon<=0  ) call Assert("Invalid number of contraction coefficients")
    if(ncon>nexp) call Assert("Invalid number of contraction coefficients")
    if(size(vexp)  <nexp) call Assert("Inconsistent exponent-list size")
    if(size(mcoe,1)<nexp) call Assert("Inconsistent coeff-list size")
    if(size(mcoe,2)<ncon) call Assert("Inconsistent coeff-list size")

    call self.free()

    self.l    = l
    self.nexp = nexp
    self.ncon = ncon
    allocate(self.vexp(nexp))
    self.vexp=vexp(1:nexp)
    allocate(self.mcoe(nexp,ncon))
    self.mcoe=mcoe(1:nexp,1:ncon)
 
    self.INITIALIZED = .TRUE.

  end subroutine ClassQCShell_Init

  !.. Free a shell
  subroutine ClassQCShell_Free( self )
    class( ClassQCShell ), intent(inout) :: self
    self.l    = 0
    self.nexp = 0
    self.ncon = 0
    if(allocated(self.vexp)) deallocate(self.vexp)
    if(allocated(self.mcoe)) deallocate(self.mcoe)
    self.INITIALIZED = .FALSE.
  end subroutine ClassQCShell_Free

  !.. Write the content of the shell in formatted form
  !   to standard output or to an optional assigned unit
  !..
  integer function ClassQCShell_GetNcon( self ) result( res )
    class( ClassQCShell ), intent(inout) :: self
    res = ( 2 * self.l + 1 ) * self.ncon
  end function ClassQCShell_GetNcon

  subroutine ClassQCShell_GetCont( self, nExp, vExp, nCon, mCoe ) 
    !
    class( ClassQCShell )         , intent(in)    :: self
    integer                       , intent(out)   :: nExp
    real(kind(1d0)), allocatable  , intent(inout) :: vExp(:)
    integer                       , intent(out)   :: nCon
    real(kind(1d0)), allocatable  , intent(inout) :: mCoe(:,:)
    !
    if( .not. self.INITIALIZED )return
    nExp = self.nExp
    nCon = self.nCon
    if(allocated(vExp))deallocate(vExp)
    if(allocated(mCoe))deallocate(mCoe)
    allocate(vExp,source=self.vExp)
    allocate(mCoe,source=self.mCoe)
    !
  end subroutine ClassQCShell_GetCont


  !.. Write the content of the shell in formatted form
  !   to standard output or to an optional assigned unit
  !..
  subroutine ClassQCShell_Show( self, unit_ )
    class( ClassQCShell ), intent(inout) :: self
    integer, optional  , intent(in)    :: unit_
    integer :: unit, icon, iexp

    unit = OUTPUT_UNIT 
    if( present( unit_ ) ) unit = unit_

    write(unit, "(a)"    ) " = shell = "
    write(unit, "(a,i0)" ) "  l    = ",self.l
    write(unit, "(a,i0)" ) "  nexp = ",self.nexp
    write(unit, "(a,i0)" ) "  ncon = ",self.ncon
    do iexp = 1, self.nexp
       write(unit, "(a,*(d24.15,x))") "  ",&
            self.vexp(iexp), &
            (self.mcoe(iexp,icon),icon=1,self.ncon)
    enddo
    write(unit, "(a)"    ) " -------- "
  end subroutine ClassQCShell_Show


  !.. Parse a text with shell information with the Dalton convention
  subroutine ClassQCShell_Parse( self, text )

    class( ClassQCShell ), intent(inout) :: self
    character(len=*)     , intent(in)    :: text

    character(len=*), parameter :: LANGSTRN = "spdfghijklmno"

   integer                      :: l, nexp, ncon
    real(kind(1d0)), allocatable :: vexp(:)
    real(kind(1d0)), allocatable :: mcoe(:,:)
    integer             :: ich1, idch, ich, lch, icon, iexp
    character(len=10000) :: line

    ich1 = 1

    !.. Read the angular momentum
    idch = index(text(ich1:),QC_NEW_LINE) - 2
    line = text(ich1:ich1+idch)
    ich=index(line,"!")
    line=adjustl(line(ich+1:))
    l = index(LANGSTRN,line(1:1)) - 1
    if(l < 0) call Assert("Invalid number of exponents")
    ich1 = ich1 + idch + len(QC_NEW_LINE) + 1

    !.. Read the number of exponents and contracted functions
    idch = index(text(ich1:),QC_NEW_LINE) - 2
    line = text(ich1:ich1+idch)
    ich=index(line,"H")
    line=adjustl(line(ich+1:))
    read(line,*) nexp, ncon
    if(nexp<=0  ) call Assert("Invalid number of exponents")
    if(ncon<=0  ) call Assert("Invalid number of contraction coefficients")
    if(ncon>nexp) call Assert("Invalid number of contraction coefficients")
    ich1 = ich1 + idch + len(QC_NEW_LINE) + 1

    !.. Read the exponents and the contraction coefficents
    allocate(vexp(nexp),mcoe(nexp,ncon))
    do iexp = 1, nexp
       idch = index(text(ich1:),QC_NEW_LINE)
       if(idch<=0)idch=len_trim(text(ich1:))+1
       idch = idch - 2
       line=trim(adjustl(text(ich1:ich1+idch)))
       read(line,*) vexp(iexp), (mcoe(iexp,icon),icon=1,ncon)
       ich1 = ich1 + idch + len(QC_NEW_LINE) + 1
    enddo

    call self.init(l,nexp,ncon,vexp,mcoe)
    
  end subroutine ClassQCShell_Parse


  !.. Parse a text with the info for all the shells of an atom
  subroutine ClassQCAtomBasis_Parse( self, text, iChStart )

    class( ClassQCAtomBasis ), intent(inout) :: self
    character(len=*)         , intent(in)    :: text
    integer, optional        , intent(in)    :: iChStart

    integer :: ich1, ich2, idch, ich, lch, icon, iexp
    integer :: lmax, lsh, iAtom
    character(len=10000) :: line


    call self.free()

    ich1 = 1
    if( present( iChStart ) )then
       if( iChStart < 1 .or. iChStart > len_trim(text) )call Assert("Error in call to parse atom")
       ich1 = iChStart
    endif

    !.. Reads the header
    idch = index(text(ich1:),QC_NEW_LINE) - 2
    line = text(ich1:ich1+idch)
    ich=index(line,"!")
    line=adjustl(line(ich+1:))
    ich=index(line," ")
    allocate(self.name,source=line(1:ich-1))
    ich1 = ich1 + idch + len(QC_NEW_LINE) + 1
    
    !.. Determines the charge
    chargeCycle: do iAtom = iATOM_LB, NATOMS
       if(self.name  == LIST_ATOMIC_NAMES(   iAtom ) )then
          self.charge = LIST_ATOMIC_CHARGES( iAtom )
          exit chargeCycle
       endif
    enddo chargeCycle

    !.. Counts the number of shells
    lmax = -1
    ich=ich1
    do
       idch = index(text(ich:),"functions") - 2 
       if(idch <=0 ) exit
       lmax = lmax + 1
       ich = ich + idch + len("functions") + 1
    enddo
    if(lmax < 0) call Assert("Invalid number of shells")
    self.lmax = lmax

    !.. Cycle over the shells
    allocate(self.vshell(0:self.lmax))
    do lsh = 0, self.lmax
       

       !.. Identifies the beginning of the shell paragraph
       idch = index(text(ich1:),"functions") - 2
       ich  = index(text(1:ich1+idch),QC_NEW_LINE,BACK=.TRUE.) - 2
       ich1 = ich + len(QC_NEW_LINE) + 1

       !.. Identifies the end of the shell paragraph
       if(lsh < self.lmax)then
          ich  = index(text(ich1:),QC_NEW_LINE)+ich1-1
          idch = index(text(ich:),"functions") - 2
          ich2 = index(text(ich:ich+idch),QC_NEW_LINE,BACK=.TRUE.) - 2
          ich2 = ich + ich2 + len(QC_NEW_LINE)
       else
          ich2 = len_trim(text)
       endif

       line=text(ich1:ich2)
       call self.vshell(lsh).parse(line)
       
       ich1 = ich2 + 1
       
    enddo
    self.INITIALIZED = .TRUE.
    
  end subroutine ClassQCAtomBasis_Parse
  
  integer function ClassQCAtomBasis_GetNsh( self ) result( res )
    class( ClassQCAtomBasis ), intent(inout) :: self
    integer :: lsh
    res = 0
    if( .not. self.INITIALIZED ) return
    res = self.lmax + 1
  end function ClassQCAtomBasis_GetNsh

  subroutine ClassQCAtomBasis_GetCont( self, l, nExp, vExp, nCon, mCoe ) 
    !
    class( ClassQCAtomBasis )     , intent(in)    :: self
    integer                       , intent(in)    :: l
    integer                       , intent(out)   :: nExp
    real(kind(1d0)), allocatable  , intent(inout) :: vExp(:)
    integer                       , intent(out)   :: nCon
    real(kind(1d0)), allocatable  , intent(inout) :: mCoe(:,:)
    !
    if( .not. self.INITIALIZED )return
    call self.vShell( l ).GetCont( nExp, vExp, nCon, mCoe )
    !
  end subroutine ClassQCAtomBasis_GetCont


  subroutine ClassQCAtomBasis_Show( self, unit_ )
    class( ClassQCAtomBasis ), intent(inout) :: self
    integer, optional   , intent(in)    :: unit_
    integer :: unit, lsh

    if( .not. self.INITIALIZED ) return

    unit = OUTPUT_UNIT 
    if( present( unit_ ) ) unit = unit_
    
    write(unit, "(a)"    ) " = atom = "
    write(unit, "(a)"    ) "  name = "//self.name
    write(unit, "(a,i0)" ) "  Z    = ",self.charge
    write(unit, "(a,i0)" ) "  lmax = ",self.lmax
    do lsh = 0, self.lmax
       call self.vshell(lsh).show(unit)
    enddo
    write(unit, "(a)"    ) " -------- "
    
  end subroutine ClassQCAtomBasis_Show


  subroutine ClassQCAtomBasis_Free( self )
    class( ClassQCAtomBasis ), intent(inout) :: self
    integer :: lsh

    if(allocated(self.name))deallocate(self.name)
    self.number = 0
    if(associated(self.vshell))then
       do lsh = 0, self.lmax
          call self.vshell(lsh).free()
       enddo
       deallocate(self.vshell)
    endif
    self.lmax=-1
    self.INITIALIZED = .FALSE.
  end subroutine ClassQCAtomBasis_Free
  
  

!====================================================================
!    ATOMIC-SET BASIS
!====================================================================

  !.. Parse a text with the info for all the shells of an atom
  subroutine ClassQCAtomicSetBasis_ParseFile( self, FileName )

    class( ClassQCAtomicSetBasis ), intent(inout) :: self
    character(len=*)         , intent(in)    :: FileName

    character(len=:), allocatable :: text
    character(len=10000)          :: AtomicText
    integer                       :: iAtom, iAtom2, iChStart, DeltaCh, ich

    call self.free()
    call GetFullText( FileName, text, " "//QC_NEW_LINE )

    self.nat = NATOMS
    allocate(self.vAtom( self.nat ))
    do iAtom = iATOM_LB, self.nat
       
       !.. Determine the beginning of the atom specification
       iChStart = index( text, "! "//LIST_ATOMIC_NAMES( iAtom ) )
       if( iChStart <= 0 )cycle

       !.. Determines the end of the atom specification, by 
       !   cycling over all possible other atoms (ugh!)
       !..
       DeltaCh = len_trim(text) - iChStart
       do iAtom2 = iATOM_LB, self.nat
          if(iAtom2 == iAtom)cycle
          ich = index( text(iChStart:iChStart+DeltaCh), "! "//LIST_ATOMIC_NAMES( iAtom2 ) ) - 2
          if(ich < 0 )cycle
          DeltaCh = min( ich, DeltaCh ) 
       enddo
       AtomicText = text(iChStart:iChStart+DeltaCh)
       call self.vAtom( iAtom ).parse( AtomicText )
    enddo

    self.INITIALIZED = .TRUE.

  end subroutine ClassQCAtomicSetBasis_ParseFile
  

  integer function ClassQCAtomicSetBasis_GetNsh( self, AtomicNumber ) result( res )
    class( ClassQCAtomicSetBasis ), intent(in) :: self
    integer                       , intent(in) :: AtomicNumber
    res = 0
    if( .not. self.INITIALIZED )return
    res = self.vAtom( AtomicNumber ).GetNsh()
  end function ClassQCAtomicSetBasis_GetNsh

  integer function ClassQCAtomicSetBasis_GetLmax( self, AtomicNumber ) result( res )
    class( ClassQCAtomicSetBasis ), intent(in) :: self
    integer                       , intent(in) :: AtomicNumber
    res = 0
    if( .not. self.INITIALIZED )return
    res = self.vAtom( AtomicNumber ).lmax
  end function ClassQCAtomicSetBasis_GetLmax

  subroutine ClassQCAtomicSetBasis_GetCont( self, AtomicNumber, l, nExp, vExp, nCon, mCoe ) 
    !
    class( ClassQCAtomicSetBasis ), intent(in)    :: self
    integer                       , intent(in)    :: AtomicNumber
    integer                       , intent(in)    :: l
    integer                       , intent(out)   :: nExp
    real(kind(1d0)), allocatable  , intent(inout) :: vExp(:)
    integer                       , intent(out)   :: nCon
    real(kind(1d0)), allocatable  , intent(inout) :: mCoe(:,:)
    !
    if( .not. self.INITIALIZED )return
    call self.vAtom( AtomicNumber ).GetCont( l, nExp, vExp, nCon, mCoe )
    !
  end subroutine ClassQCAtomicSetBasis_GetCont


  subroutine ClassQCAtomicSetBasis_Show( self, unit_ )
    class( ClassQCAtomicSetBasis ), intent(in) :: self
    integer, optional             , intent(in) :: unit_

    integer :: unit, iAtom

    if( .not. self.INITIALIZED )return

    unit = OUTPUT_UNIT 
    if( present( unit_ ) ) unit = unit_

    do iAtom = 0, self.nat
       call self.vAtom( iAtom ).show( unit )
    enddo

  end subroutine ClassQCAtomicSetBasis_Show


  subroutine ClassQCAtomicSetBasis_Free( self )
    class( ClassQCAtomicSetBasis ), intent(inout) :: self

    integer :: iAtom

    if( .not. self.INITIALIZED ) return

    do iAtom = 0, self.nat
       call self.vAtom( iAtom ).free()
    enddo
    deallocate(self.vAtom)
    self.nat = -1
    self.INITIALIZED = .FALSE.

  end subroutine ClassQCAtomicSetBasis_Free
  
  


end module ModuleQCBasis
