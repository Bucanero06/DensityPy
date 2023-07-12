!! CONFIDENTIAL
! {{{ Detailed information
!> \file
!!
!! Module containing the neccesary classes that manage the scattering states.
!!
!! OUTPUT FILES:
!!
!! - See ClassScatteringStatePrintToUnit for the scattering observables (phaseshifts and eigenchannels)
!!
! }}}
module ModuleScatteringStates

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleIO
  use ModuleString
  use ModuleNameConventions
  use ModuleConstants
  use ModuleSystemUtils

  implicit none

  private

  integer, parameter :: MAX_CHANNEL_NAME_LENGHT = 32

  type, public :: ClassScatteringStateList
     !
     private
     Class(ClassScatteringState), pointer     :: First => NULL()
     Class(ClassScatteringState), pointer     :: Last  => NULL()
     Class(ClassScatteringState), pointer     :: Curr  => NULL()
     character(len=:)           , allocatable :: dir
     integer :: nItems
     !
   contains
     !
     procedure, public :: Init                  => ClassScatteringStateListInit
     procedure, public :: Add                   => ClassScatteringStateListAdd
     procedure, public :: EnforceContinuity     => ClassScatteringStateListEnforceContinuity
     procedure, public :: Save                  => ClassScatteringStateListSave
     procedure, public :: Load                  => ClassScatteringStateListLoad
     procedure, public :: Transform             => ClassScatteringStateListTransform
     procedure, public :: GetBestUniformFilling => ClassScatteringStateListGetBestUniformFilling
     procedure, public :: GetBestRydbergFilling => ClassScatteringStateListGetBestRydbergFilling
     procedure, public :: GetResolveFilling     => ClassScatteringStateListGetResolveFilling
     procedure, public :: GetRefineFilling      => ClassScatteringStateListGetRefineFilling
     procedure, public :: GetEnergyList         => ClassScatteringStateListGetEnergyList
     procedure, public :: GetnOpenList          => ClassScatteringStateListGetnOpenList
     procedure, public :: GetChParams           => ClassScatteringStateListGetChParams
     procedure, public :: Free                  => ClassScatteringStateListFree
     procedure, public :: ProjectOnBra          => ClassScatteringStateListProjectOnBra
     procedure, public :: GetPsiMinus           => ClassScatteringStateListGetPsiMinus
     procedure, public :: SetPsiMinus           => ClassScatteringStateListSetPsiMinus
     !
     final :: ClassScatteringStateListFinalize
     !
  end type ClassScatteringStateList
  

  type, public :: ClassScatteringState
     !
     private
     !
     class(ClassScatteringState), pointer :: Prev => NULL()
     class(ClassScatteringState), pointer :: Next => NULL()
     !
     integer                         :: nOpen
!!$     integer           , allocatable :: listOpen(:)
     logical                         :: solved
     real   (kind(1d0))              :: Energy
     integer                         :: nBoxStates
     !> RawVec are the scattering states that are
     !! fitted and behave asymptotically as 
     !! Psi = F A + G B
     integer           , allocatable :: lw(:)
     real   (kind(1d0)), allocatable :: IonEnergy(:)
     character(len=MAX_CHANNEL_NAME_LENGHT), allocatable :: Name(:)
     real   (kind(1d0)), allocatable :: A(:,:)
     real   (kind(1d0)), allocatable :: B(:,:)
     real   (kind(1d0)), allocatable :: Fitness(:,:)
     complex(kind(1d0)), allocatable :: S(:,:)
     real   (kind(1d0)), allocatable :: K(:,:)
     real   (kind(1d0)), allocatable :: PHI(:)
     real   (kind(1d0)), allocatable :: EIGCH(:,:)
     real   (kind(1d0)), allocatable :: RawVec(:,:)
     complex(kind(1d0)), allocatable :: PsiMinus(:,:)     
     !
   contains
     !
     procedure, public  :: Init     => ClassScatteringStateInit
     generic  , public  :: Save     => ClassScatteringStateSaveToUnit
     generic  , public  :: Print    => ClassScatteringStatePrintToUnit
     procedure, public  :: Solve    => ClassScatteringStateSolve
     procedure, public  :: Transform=> ClassScatteringStateTransform
     procedure, public  :: IsSolved => ClassScatteringStateIsSolved
     procedure, public  :: AlignToPrev => ClassScatteringStateAlignToPrev
     procedure, public  :: Setlw    => ClassScatteringStateSetlw
     procedure, public  :: SetAsymp => ClassScatteringStateSetAsymp
     procedure, public  :: SetThr   => ClassScatteringStateSetThr
     procedure, public  :: SetName => ClassScatteringStateSetName
     procedure, public  :: SetCoef  => ClassScatteringStateSetCoef
     procedure, public  :: GetTotPh => ClassScatteringStateGetTotPh
     procedure, public  :: GetEnergy => ClassScatteringStateGetEnergy
     procedure, public  :: ProjectOnBra => ClassScatteringStateProjectOnBra
     procedure, public  :: GetPsiMinus  => ClassScatteringStateGetPsiMinus
     procedure, public  :: SetPsiMinus  => ClassScatteringStateSetPsiMinus
     procedure, public  :: Free     => ClassScatteringStateFree
     procedure, public  :: PrintChannelInfo => ClassScatteringStatePrintChannelInfo
     !
     procedure, private :: ClassScatteringStateSaveToUnit
     procedure, private :: ClassScatteringStatePrintToUnit
     !
     final :: ClassScatteringStateFinalize
     !
  end type ClassScatteringState

  interface LoadScatState
     module procedure ClassScatteringStateLoadFromUnit
  end interface LoadScatState
  
  interface operator (.precedes.)
     module procedure ClassScatteringStatePrecedes
  end interface operator (.precedes.)

contains


  !===================================================================
  !===================================================================
  !!
  !!           CLASSSCATTERINGSTATELIST  PROCEDURES
  !!
  !===================================================================
  !===================================================================

  subroutine ClassScatteringStateListInit( self, dir )
    use ModuleString
    implicit none
    class(ClassScatteringStateList), intent(inout) :: self
    character(len=*)               , intent(in)    :: dir
    call self.free()
    self.dir=FormatAsDir(dir)
    self.nItems=0
  end subroutine ClassScatteringStateListInit
  
  subroutine ClassScatteringStateListEnforceContinuity( Self )
    implicit none
    class(ClassScatteringStateList)     , intent(inout) :: Self
    class(ClassScatteringState), pointer :: Ptr
    if( .not.associated(Self%First) )return
    Ptr=>Self%First
    do while(associated(Ptr%Next))
       Ptr=>Ptr%Next
       call Ptr.AlignToPrev()
    enddo
  end subroutine ClassScatteringStateListEnforceContinuity
    
  
  subroutine ClassScatteringStateListAdd( Self, Item )

    implicit none

    class(ClassScatteringStateList)    , intent(inout) :: Self
    type(ClassScatteringState), pointer, intent(in)    :: Item
    type(ClassScatteringState), pointer :: Ptr

    if( .not.associated(Self%Curr) )then

       Item%Prev  => Null()
       Item%Next  => Null()
       Self%First => Item
       Self%Last  => Item
       Self%Curr  => Item
       
    else

       Ptr => Self%Curr

       !.. If Item < Curr, move Curr left
       if( Item .precedes. Ptr )then
          
          do while( Item .precedes. Ptr )
             Ptr => Ptr.Prev
             if( .not. associated( Ptr ) ) exit
          enddo

          if(.not.associated(Ptr))then

             !.. Add Item at the beginning of the list
             Item%Prev => Null()
             Item%Next => Self%First
             Self%First%Prev => Item
             Self%First=> Item
             Self%Curr => Item

          else

             !.. Add Item right after Ptr, which now precedes Item
             Item%Prev => Ptr
             Item%Next => Ptr%Next
             Ptr%Next%Prev => Item
             Ptr%Next  => Item

          endif

       elseif( Ptr .precedes. Item )then

          do while( Ptr .precedes. Item )
             Ptr => Ptr.Next
             if( .not. associated( Ptr ) ) exit
          enddo

          if(.not.associated(Ptr))then

             !.. Add Item at the end of the list
             Item%Prev      => Self%Last
             Item%Next      => Null()
             Self%Last%Next => Item
             Self%Last      => Item
             Self%Curr      => Item

          else

             !.. Add Item right before Ptr, which now follows Item
             Item%Prev      => Ptr%Prev
             Item%Next      => Ptr
             Item%Prev%Next => Item
             Item%Next%Prev => Item
             Self%Curr      => Item

             do 
                Ptr => Ptr.Next
                if( .not. associated( Ptr ) ) exit
                !write(*,*) Ptr.GetEnergy()
             enddo
             
          endif

       endif

    endif

    Self.nItems=Self.nItems+1

  end subroutine ClassScatteringStateListAdd

  subroutine ClassScatteringStateListSave( self, label )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    character(len=*)                 , intent(in)    :: label
    !
    type( ClassScatteringState ), pointer :: Item
    character(len=:)        , allocatable :: fmt_File, bin_File, info_File
    integer                               :: fmt_uid,  bin_uid,  info_uid
    !
    call system("mkdir -p "//self.dir)
    allocate(  fmt_File, source = self.dir // SCAT_FILE_ROOT // label // FMT_SUFFIX )
    allocate(  bin_File, source = self.dir // SCAT_FILE_ROOT // label // BIN_SUFFIX )
    call OpenFile(  fmt_File,  fmt_uid, "write",   "formatted" )
    call OpenFile(  bin_File,  bin_uid, "write", "unformatted" )
    !
    Item => self%First
    do while( associated( Item ) )
       call Item%Print( fmt_uid )
       call Item%Save( bin_uid )
       Item => Item%Next
    enddo
    !
    allocate( info_File, source = self.dir // SCAT_FILE_ROOT // label // ".info" )
    call OpenFile( info_File, info_uid, "write",   "formatted" )
    Item => self%Last
    call Item%PrintChannelInfo( info_uid )
    !
  end subroutine ClassScatteringStateListSave


  subroutine ClassScatteringStateListTransform( self, dMat )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    real(kind(1d0))                  , intent(in)    :: dMat(:,:)
    !
    type( ClassScatteringState ), pointer :: Item
    !
    Item => self%First
    do while( associated( Item ) )
       call Item%Transform( dMat )
       Item => Item%Next
    enddo
    !
  end subroutine ClassScatteringStateListTransform


  !.. Evaluates \langle Bra|\Psi^-_{\alpha E}\rangle
  subroutine ClassScatteringStateListGetEnergyList( self, Evec )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    real   (kind(1d0)), allocatable  , intent(out)   :: Evec(:)  
    integer :: nEnergies, iItem
    type( ClassScatteringState ), pointer :: Item
    !
    nEnergies=self%nItems
    call realloc( Evec, nEnergies )
    Item => self%First
    iItem=0
    do while( associated( Item ) )
       iItem=iItem+1
       Evec( iItem ) = Item%GetEnergy()
       Item => Item%Next
    enddo
  end subroutine ClassScatteringStateListGetEnergyList

  !.. Evaluates \langle Bra|\Psi^-_{\alpha E}\rangle
  subroutine ClassScatteringStateListGetnOpenList( self, nOpen )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    integer           , allocatable  , intent(out)   :: nOpen(:)  
    integer :: nEnergies, iItem
    type( ClassScatteringState ), pointer :: Item
    !
    nEnergies=self%nItems
    call realloc( nOpen, nEnergies )
    Item => self%First
    iItem=0
    do while( associated( Item ) )
       iItem=iItem+1
       nOpen( iItem ) = Item%nOpen
       Item => Item%Next
    enddo
  end subroutine ClassScatteringStateListGetnOpenList

  !.. Returns the channel parameters 
  subroutine ClassScatteringStateListGetChParams( self, nChan, &
    Pion, Mion, Lion, Nion, Eion, lwav )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    integer                          , intent(out)   :: nChan
    integer           , allocatable  , intent(out)   :: Pion(:) 
    integer           , allocatable  , intent(out)   :: Mion(:) 
    integer           , allocatable  , intent(out)   :: Lion(:) 
    integer           , allocatable  , intent(out)   :: Nion(:) 
    real(kind(1d0))   , allocatable  , intent(out)   :: Eion(:) 
    integer           , allocatable  , intent(out)   :: lwav(:) 

    character(len=*), parameter :: Lstrn="SPDFGHIJKLMNO"
    character(len=*), parameter :: Pstrn="eo"
    integer :: iOpen, i 
    type( ClassScatteringState ), pointer :: Item
    character(len=10), allocatable :: chName
    !
    if(allocated(Pion))deallocate(Pion)
    if(allocated(Mion))deallocate(Mion)
    if(allocated(Lion))deallocate(Lion)
    if(allocated(Nion))deallocate(Nion)
    if(allocated(Eion))deallocate(Eion)
    if(allocated(lwav))deallocate(lwav)
    !
    Item => self%Last
    nChan = Item%nOpen
    allocate( Pion( nChan ) )
    allocate( Mion( nChan ) )
    allocate( Lion( nChan ) )
    allocate( Nion( nChan ) ) 
    allocate( Eion( nChan ) ) 
    allocate( lwav( nChan ) )
    !
    do iOpen = 1, nChan
       chName = Item.Name( iOpen )
       i=index(chName,".")
       read(chName(1:i-1),*) Nion(iOpen)
       chName=chName(i+1:)
       read(chName(1:1),*) Mion(iOpen)
       Lion(iOpen) = index(Lstrn,chName(2:2))-1
       Pion(iOpen) = index(Pstrn,chName(3:3))-1
       Eion(iOpen) = Item.IonEnergy( iOpen ) 
       lwav(iOpen) = Item.lw( iOpen )
    enddo
    !
  end subroutine ClassScatteringStateListGetChParams

  !.. Evaluates \langle Bra|\Psi^-_{\alpha E}\rangle
  subroutine ClassScatteringStateListProjectOnBra( self, zBra, Projection )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    complex(kind(1d0))               , intent(in)    :: zBra(:)
    complex(kind(1d0)), allocatable  , intent(out)   :: Projection(:,:)
    !
    integer :: nEnergies, nOpenMax, iItem
    type( ClassScatteringState ), pointer :: Item
    complex(kind(1d0)), allocatable :: zOpen(:)
    !
    nEnergies = self%nItems
    nOpenMax = self%last%nopen
    call realloc( Projection, nEnergies, nOpenMax )
    allocate( zOpen( nOpenMax ) )
    !
    Item => self%First
    iItem=0
    do while( associated( Item ) )
       iItem=iItem+1
       call Item%ProjectOnBra( zBra, zOpen )
       Projection( iItem, : ) = zOpen
       Item => Item%Next
    enddo
    deallocate( zOpen )
    !
  end subroutine ClassScatteringStateListProjectOnBra

  !.. Evaluates \langle Bra|\Psi^-_{\alpha E}\rangle
  subroutine ClassScatteringStateListGetPsiMinus( self, iEn, iCh, zKet )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    integer                          , intent(in)    :: iEn, iCh
    complex(kind(1d0))               , intent(out)   :: zKet(:)
    !
    integer :: iItem
    type( ClassScatteringState ), pointer :: Item
    !
    Item => self%First
    iItem=0
    do while( associated( Item ) )
       iItem=iItem+1
       if(iItem==iEn)then
          call Item%GetPsiMinus( iCh, zKet )
          exit
       endif
       Item => Item%Next
    enddo
    !
  end subroutine ClassScatteringStateListGetPsiMinus

  !.. Evaluates \langle Bra|\Psi^-_{\alpha E}\rangle
  subroutine ClassScatteringStateListSetPsiMinus( self, iEn, iCh, zKet )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    integer                          , intent(in)    :: iEn, iCh
    complex(kind(1d0))               , intent(in)    :: zKet(:)
    !
    integer :: iItem
    type( ClassScatteringState ), pointer :: Item
    !
    Item => self%First
    iItem=0
    do while( associated( Item ) )
       iItem=iItem+1
       if(iItem==iEn)then
          call Item%SetPsiMinus( iCh, zKet )
          exit
       endif
       Item => Item%Next
    enddo
    !
  end subroutine ClassScatteringStateListSetPsiMinus

  subroutine ClassScatteringStateListLoad( self, label )
    implicit none
    class( ClassScatteringStateList ), intent(inout) :: self
    character(len=*)                 , intent(in)    :: label
    !
    type( ClassScatteringState ), pointer :: Item
    character(len=:)        , allocatable :: bin_File
    integer                               :: bin_uid
    logical :: exist
    !
    allocate( bin_File, source = self.dir // SCAT_FILE_ROOT // label // BIN_SUFFIX )
    inquire( file=bin_file, exist=exist) 
    if(.not.exist)return
    call OpenFile( bin_File, bin_uid, "read", "unformatted" )
    !
    do 
       Item => LoadScatState( bin_uid )
       if( .not. associated( Item ) ) exit
       call self.Add( Item )
    enddo
    !
  end subroutine ClassScatteringStateListLoad

  subroutine ClassScatteringStateListFree( Self )
    class( ClassScatteringStateList ), intent(inout) :: Self
    type( ClassScatteringState ), pointer :: Item, Next
    Item => Self%First
    do while(associated(Item))
       Next => Item%Next
       call Item%Free()
       deallocate(Item)
       Item => Next
    enddo
    if(allocated(self.dir))deallocate(self.dir)
    Self%First => NULL()
    Self%Last  => NULL()
    Self%Curr  => NULL()
    Self%nItems=0
  end subroutine ClassScatteringStateListFree

  subroutine ClassScatteringStateListFinalize( self )
    type( ClassScatteringStateList ), intent(inout) :: self
    call self.free()
  end subroutine ClassScatteringStateListFinalize

  subroutine ClassScatteringStateListGetBestUniformFilling( Self, Emin, Emax, dEmax, nE, Evec )
    implicit none
    class(ClassScatteringStateList), intent(inout) :: Self
    real(kind(1d0))                , intent(in)    :: Emin
    real(kind(1d0))                , intent(in)    :: Emax
    real(kind(1d0))                , intent(in)    :: dEmax
    integer                        , intent(out)   :: nE
    real(kind(1d0)), allocatable   , intent(out)   :: Evec(:)
    class(ClassScatteringState), pointer :: Ptr1, Ptr2, Ptr
    integer :: iEn, iEn0, Delta_nE
    if( .not.associated( Self%First ) )then
       nE = max( nint( ( Emax - Emin ) / dEmax ) + 1, 2 )
       allocate( Evec, source = [ ( Emin + ( Emax - Emin ) / dble( nE - 1 ) * dble( iEn - 1 ), iEn = 1, nE ) ] )
       return
    endif
    if( Self%First%GetEnergy() > Emax .or. Self%Last%GetEnergy() < Emin )then
       nE = max( nint( ( Emax - Emin ) / dEmax ) + 1, 2 )
       allocate( Evec, source = [ ( Emin + ( Emax - Emin ) / dble( nE - 1 ) * dble( iEn - 1 ), iEn = 1, nE ) ] )
       return
    endif
    !.. Find smallest available energy above Emin
    Ptr1 => Self%First
    do while( Ptr1%GetEnergy() < Emin )
       Ptr1 => Ptr1%Next
    enddo
    !.. Find largest available energy below Emax
    Ptr2 => Self%Last
    do while( Ptr2%GetEnergy() > Emax )
       Ptr2 => Ptr2%Prev
    enddo
    !.. If there are no points between Emin and Emax, just takes a single uniform grid
    if( Ptr2%GetEnergy() < Ptr1%GetEnergy() )then
       nE = max( nint( ( Emax - Emin ) / dEmax ) + 1, 2 )
       allocate( Evec, source = [ ( Emin + ( Emax - Emin ) / dble( nE - 1 ) * dble( iEn - 1 ), iEn = 1, nE ) ] )
       return
    endif
    !
    !.. At this point, there is at least one available energy strictly between Emin and Emax
    !
    !.. Counts number of points
    !.. Left interval
    nE = int( ( Ptr1%GetEnergy() - Emin ) / dEmax ) + 1
    !.. Intermediate intervals
    Ptr => Ptr1
    do while( Ptr%GetEnergy() < Ptr2%GetEnergy() )
       nE= nE + int( ( Ptr%Next%GetEnergy() - Ptr%GetEnergy() ) / dEmax ) 
       Ptr => Ptr.Next
    enddo
    !.. Last interval
    nE = nE + int( ( Emax - Ptr2%GetEnergy() ) / dEmax ) + 1
    !
    !.. determines the vector of energies
    allocate( Evec( nE ) )
    !.. Left interval
    Delta_nE = int( ( Ptr1%GetEnergy() - Emin ) / dEmax ) + 1
    iEn0 = 0
    do iEn = 1, Delta_nE !.. Does not fill the right bound
       Evec( iEn0 + iEn ) = Emin + ( Ptr1%GetEnergy() - Emin ) / dble( Delta_nE ) * dble( iEn - 1 ) 
    enddo
    iEn0 = iEn0 + Delta_nE
    !.. Intermediate intervals
    Ptr => Ptr1
    do while( Ptr%GetEnergy() < Ptr2%GetEnergy() )
       Delta_nE = int( ( Ptr%Next%GetEnergy() - Ptr%GetEnergy() ) / dEmax ) 
       do iEn = 1, Delta_nE
          Evec( iEn0 + iEn ) = Ptr%GetEnergy() + ( Ptr%Next%GetEnergy() - Ptr%GetEnergy() ) / dble( Delta_nE + 1 ) * dble( iEn )
       enddo
       iEn0 = iEn0 + Delta_nE
       Ptr => Ptr.Next
    enddo
    !.. Last interval
    Delta_nE = int( ( Emax - Ptr2%GetEnergy() ) / dEmax ) + 1
    do iEn = 1, Delta_nE
       Evec( iEn0 + iEn ) = Ptr2%GetEnergy() + ( Emax - Ptr2%GetEnergy() ) / dble( Delta_nE ) * dble( iEn )
    enddo
  end subroutine ClassScatteringStateListGetBestUniformFilling


  subroutine ClassScatteringStateListGetBestRydbergFilling( Self, Emin, Emax, dpqnmax, IonCharge, Ethr, nE, Evec )
    implicit none
    class(ClassScatteringStateList), intent(inout) :: Self
    real(kind(1d0))                , intent(in)    :: Emin
    real(kind(1d0))                , intent(in)    :: Emax
    real(kind(1d0))                , intent(in)    :: dpqnmax
    real(kind(1d0))                , intent(in)    :: IonCharge
    real(kind(1d0))                , intent(in)    :: Ethr
    integer                        , intent(out)   :: nE
    real(kind(1d0)), allocatable   , intent(out)   :: Evec(:)
    real(kind(1d0)) :: pqnmin, pqnmax, pqn, pqn1, pqn2
    class(ClassScatteringState), pointer :: Ptr1, Ptr2, Ptr
    integer :: iEn, iEn0, Delta_nE
    pqnmin = E2pqn(Emin,IonCharge,Ethr)
    pqnmax = E2pqn(Emax,IonCharge,Ethr)
    if( .not.associated( Self%First ) )then
       nE = max( nint( ( pqnmax - pqnmin ) / dpqnmax ) + 1, 2 )
       allocate( Evec( nE ) )
       do iEn = 1, nE
          pqn = pqnmin + ( pqnmax - pqnmin ) / dble( nE - 1 ) * dble( iEn - 1 )
          Evec( iEn ) = pqn2E( pqn, IonCharge, Ethr )
       enddo
       return
    endif
    if( Self%First%GetEnergy() > Emax .or. Self%Last%GetEnergy() < Emin )then
       nE = max( nint( ( pqnmax - pqnmin ) / dpqnmax ) + 1, 2 )
       allocate( Evec( nE ) )
       do iEn = 1, nE
          pqn = pqnmin + ( pqnmax - pqnmin ) / dble( nE - 1 ) * dble( iEn - 1 )
          Evec( iEn ) = pqn2E( pqn, IonCharge, Ethr )
       enddo
       return
    endif
    !.. Find smallest available energy above Emin
    Ptr1 => Self%First
    do while( Ptr1%GetEnergy() < Emin )
       Ptr1 => Ptr1%Next
    enddo
    !.. Find largest available energy below Emax
    Ptr2 => Self%Last
    do while( Ptr2%GetEnergy() > Emax )
       Ptr2 => Ptr2%Prev
    enddo
    !.. If there are no points between Emin and Emax, just takes a single uniform grid
    if( Ptr2%GetEnergy() < Ptr1%GetEnergy() )then
       nE = max( nint( ( pqnmax - pqnmin ) / dpqnmax ) + 1, 2 )
       allocate( Evec( nE ) )
       do iEn = 1, nE
          pqn = pqnmin + ( pqnmax - pqnmin ) / dble( nE - 1 ) * dble( iEn - 1 )
          Evec( iEn ) = pqn2E( pqn, IonCharge, Ethr )
       enddo
       return
    endif
    !
    !.. At this point, there is at least one available energy strictly between Emin and Emax
    !
    !.. Counts number of points
    !.. Left interval
    pqn2=E2pqn( Ptr1%GetEnergy(), IonCharge, Ethr )
    nE = int( ( pqn2 - pqnmin ) / dpqnmax ) + 1
    !.. Intermediate intervals
    Ptr => Ptr1
    do while( Ptr%GetEnergy() < Ptr2%GetEnergy() )
       nE= nE + int( ( E2pqn( Ptr%Next%GetEnergy(), IonCharge, Ethr ) - E2pqn( Ptr%GetEnergy(), IonCharge, Ethr ) ) / dpqnmax ) 
       Ptr => Ptr.Next
    enddo
    !.. Last interval
    pqn1=E2pqn( Ptr2%GetEnergy(), IonCharge, Ethr )
    nE = nE + int( ( pqnmax - pqn1 ) / dpqnmax ) + 1
    !
    !.. determines the vector of energies
    allocate( Evec( nE ) )
    !.. Left interval
    pqn2=E2pqn( Ptr1%GetEnergy(), IonCharge, Ethr )
    Delta_nE = int( ( pqn2 - pqnmin ) / dpqnmax ) + 1
    iEn0 = 0
    do iEn = 1, Delta_nE !.. Does not fill the right bound
       Evec( iEn0 + iEn ) = pqn2E( pqnmin + ( pqn2 - pqnmin ) / dble( Delta_nE ) * dble( iEn - 1 ), IonCharge, Ethr )
    enddo
    iEn0 = iEn0 + Delta_nE
    !.. Intermediate intervals
    Ptr => Ptr1
    do while( Ptr%GetEnergy() < Ptr2%GetEnergy() )
       pqn1 = E2pqn( Ptr%GetEnergy()     , IonCharge, Ethr ) 
       pqn2 = E2pqn( Ptr%Next%GetEnergy(), IonCharge, Ethr ) 
       Delta_nE = int( ( pqn2 - pqn1 ) / dpqnmax )
       do iEn = 1, Delta_nE
          Evec( iEn0 + iEn ) = pqn2E( pqn1 + ( pqn2 - pqn1 ) / dble( Delta_nE + 1 ) * dble( iEn ), IonCharge, Ethr )
       enddo
       iEn0 = iEn0 + Delta_nE
       Ptr => Ptr.Next
    enddo
    !.. Last interval
    pqn1 = E2pqn( Ptr2%GetEnergy(), IonCharge, Ethr )
    Delta_nE = int( ( pqnmax - pqn1 ) / dpqnmax ) + 1
    do iEn = 1, Delta_nE
       Evec( iEn0 + iEn ) = pqn2E( pqn1 + ( pqnmax - pqn1 ) / dble( Delta_nE ) * dble( iEn ), IonCharge, Ethr )
    enddo
  contains
    real(kind(1d0)) function E2pqn(E,IonCharge,Eth) result(pqn)
      implicit none
      real(kind(1d0)), intent(in) :: E, IonCharge, Eth
      pqn = IonCharge / sqrt( 2.d0 * ( Eth - E ) )
    end function E2pqn
    real(kind(1d0)) function pqn2E(pqn,IonCharge,Eth) result(E)
      implicit none
      real(kind(1d0)), intent(in) :: pqn, IonCharge, Eth
      E = Eth - 0.5d0 * ( IonCharge / pqn )**2
    end function Pqn2E
  end subroutine ClassScatteringStateListGetBestRydbergFilling


  subroutine ClassScatteringStateListGetResolveFilling( Self, Emin, Emax, nE, Evec )
    implicit none
    class(ClassScatteringStateList), intent(inout) :: Self
    real(kind(1d0))                , intent(in)    :: Emin
    real(kind(1d0))                , intent(in)    :: Emax
    integer                        , intent(out)   :: nE
    real(kind(1d0)), allocatable   , intent(out)   :: Evec(:)
    class(ClassScatteringState), pointer :: Ptr1, Ptr
    real(kind(1d0)) :: tph1,tph2,tph3,tph4
    nE=0
    if( .not.associated( Self%First ) )return
    if( Self%First%GetEnergy() > Emax .or. Self%Last%GetEnergy() < Emin )return
    !.. Find smallest available energy above Emin
    Ptr1 => Self%First
    do while( Ptr1%GetEnergy() < Emin )
       Ptr1 => Ptr1%Next
    enddo
    !
    !.. At this point, there is at least one available energy strictly between Emin and Emax
    !
    !.. Counts number of points
    !.. Intermediate intervals
    Ptr => Ptr1
    if( .not. associated( Ptr%Next           ) ) return
    if( .not. associated( Ptr%Next%Next      ) ) return
    if( .not. associated( Ptr%Next%Next%Next ) ) return
    !
    tph1 = Ptr%GetTotPh()
    tph2 = Ptr%Next%GetTotPh()
    tph3 = Ptr%Next%Next%GetTotPh()
    do while( Ptr%Next%Next%GetEnergy() < Emax )
       tph4 = Ptr%Next%Next%Next%GetTotPh()
       !
       !.. Add a point if the total phase has an up-down-up behavior
       !   which may indicate the presence of a resonance.
       !..
       if(  tph2 > tph1 .and. &
            tph3 < tph2 .and. &
            tph4 > tph3 )then
          nE = nE + 1
       endif
       Ptr => Ptr.Next
       if( .not. associated( Ptr%Next%Next%Next ) ) exit
       tph1=tph2
       tph2=tph3
       tph3=tph4
    enddo
    if( nE == 0 ) return
    !
    !.. determines the vector of energies
    allocate( Evec( nE ) )
    !.. Intermediate intervals
    Ptr => Ptr1
    !
    tph1 = Ptr%GetTotPh()
    tph2 = Ptr%Next%GetTotPh()
    tph3 = Ptr%Next%Next%GetTotPh()
    nE=0
    do while( Ptr%Next%Next%GetEnergy() < Emax )
       tph4 = Ptr%Next%Next%Next%GetTotPh()
       if(  tph2 > tph1 .and. &
            tph3 < tph2 .and. &
            tph4 > tph3 )then
          nE = nE + 1
          Evec( nE ) = ( Ptr%Next%Next%GetEnergy() + Ptr%Next%GetEnergy() ) / 2.d0
       endif
       Ptr => Ptr.Next
       if( .not. associated( Ptr%Next%Next%Next ) ) exit
       tph1=tph2
       tph2=tph3
       tph3=tph4
    enddo
  end subroutine ClassScatteringStateListGetResolveFilling
    
  subroutine ClassScatteringStateListGetRefineFilling( Self, Emin, Emax, dPhi, nE, Evec )
    implicit none
    class(ClassScatteringStateList), intent(inout) :: Self
    real(kind(1d0))                , intent(in)    :: Emin
    real(kind(1d0))                , intent(in)    :: Emax
    real(kind(1d0))                , intent(in)    :: dPhi
    integer                        , intent(out)   :: nE
    real(kind(1d0)), allocatable   , intent(out)   :: Evec(:)
    class(ClassScatteringState), pointer :: Ptr1, Ptr
    real(kind(1d0)) :: tph1,tph2
    nE=0
    if( .not.associated( Self%First ) )return
    if( Self%First%GetEnergy() > Emax .or. Self%Last%GetEnergy() < Emin )return
    !.. Find smallest available energy above Emin
    Ptr1 => Self%First
    do while( Ptr1%GetEnergy() < Emin )
       Ptr1 => Ptr1%Next
    enddo
    !
    !.. At this point, there is at least one available energy strictly between Emin and Emax
    !
    !.. Counts number of points
    !.. Intermediate intervals
    Ptr => Ptr1
    if( .not. associated( Ptr%Next           ) ) return
    !
    tph1 = Ptr%GetTotPh()
    do while( Ptr%Next%GetEnergy() < Emax )
       tph2 = Ptr%Next%GetTotPh()
       if( abs( tph2 - tph1 ) > dPhi ) nE = nE + 1
       Ptr => Ptr.Next
       if( .not. associated( Ptr%Next%Next ) ) exit
       tph1=tph2
    enddo
    if( nE == 0 ) return
    !
    !.. determines the vector of energies
    allocate( Evec( nE ) )
    !.. Intermediate intervals
    Ptr => Ptr1
    !
    tph1 = Ptr%GetTotPh()
    nE=0
    do while( Ptr%Next%GetEnergy() < Emax )
       tph2 = Ptr%Next%GetTotPh()
       if( abs( tph2 - tph1 ) > dPhi )then
          nE = nE + 1
          Evec( nE ) = ( Ptr%Next%GetEnergy() + Ptr%GetEnergy() ) / 2.d0
       endif
       Ptr => Ptr.Next
       if( .not. associated( Ptr%Next%Next ) ) exit
       tph1=tph2
    enddo
  end subroutine ClassScatteringStateListGetRefineFilling
    
  !===================================================================
  !===================================================================
  !!
  !!              CLASSSCATTERINGSTATE PROCEDURES
  !!
  !===================================================================
  !===================================================================

  subroutine ClassScatteringStateInit(self,Energy,nOpen,n )
    class(ClassScatteringState), intent(inout) :: self
    real(kind(1d0))            , intent(in)    :: Energy
    integer                    , intent(in)    :: nOpen
    integer                    , intent(in)    :: n
    self.Energy=Energy
    self.nopen=nopen
    self.nBoxStates=n
    self.solved=.FALSE.
    allocate(self.IonEnergy(nOpen))
    allocate(self.Name(nOpen))
    allocate(self.lw(nOpen))
    allocate(self.RawVec(n,nOpen))
    allocate(self.PsiMinus(n,nOpen))
    allocate(self.A(nOpen,nOpen))
    allocate(self.B(nOpen,nOpen))
    allocate(self.Fitness(nOpen,nOpen))
    allocate(self.S(nOpen,nOpen))
    allocate(self.K(nOpen,nOpen))
    allocate(self.PHI(nOpen))
    allocate(self.EIGCH(nOpen,nOpen))
  end subroutine ClassScatteringStateInit

  subroutine ClassScatteringStateSaveToUnit(self,uid)
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: uid
    character(len=16) :: form
    integer :: iOpen, iBox, iOpen2
    inquire(uid,form=form)
    if(trim(adjustl(form)).is."unformatted")then
       write(uid)     self.Energy, self.nBoxStates, self.nOpen, self.solved
!!$       write(uid) (   self.listOpen ( iOpen )        , iOpen  = 1, self.nOpen )
       write(uid) (   self.Name     ( iOpen )        , iOpen  = 1, self.nOpen )
       write(uid) (   self.lw       ( iOpen )        , iOpen  = 1, self.nOpen )
       write(uid) (   self.IonEnergy( iOpen )        , iOpen  = 1, self.nOpen )
       write(uid) (   self.PHI      ( iOpen )        , iOpen  = 1, self.nOpen )
       write(uid) ( ( self.A        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid) ( ( self.B        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid) ( ( self.Fitness  ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid) ( ( self.S        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid) ( ( self.K        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid) ( ( self.EIGCH    ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       do iOpen = 1, self.nOpen
          write(uid) (        self.RawVec  ( iBox, iOpen )  , iBox = 1, self.nBoxStates )
          write(uid) (  dble( self.PsiMinus( iBox, iOpen ) ), iBox = 1, self.nBoxStates )
          write(uid) ( aimag( self.PsiMinus( iBox, iOpen ) ), iBox = 1, self.nBoxStates )
       enddo
    else
       write(uid,"("//EDBL_FMT//",2(x,i0),x,l1)") self.Energy, self.nBoxStates, self.nOpen, self.solved
!!$       write(uid) (   self.listOpen ( iOpen )        , iOpen  = 1, self.nOpen )
       write(uid,FULL_STRN_SEQ_FMT) (   self.Name     ( iOpen )         , iOpen  = 1, self.nOpen )
       write(uid,FULL_INT0_SEQ_FMT) (   self.lw       ( iOpen )         , iOpen  = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) (   self.IonEnergy( iOpen )         , iOpen  = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) (   self.PHI      ( iOpen )         , iOpen  = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( ( self.A        ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( ( self.B        ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( ( self.Fitness  ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( ( dble( self.S  ( iOpen2, iOpen )), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( (aimag( self.S  ( iOpen2, iOpen )), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( ( self.K        ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       write(uid,FULL_EDBL_SEQ_FMT) ( ( self.EIGCH    ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       do iOpen = 1, self.nOpen
          write(uid,FULL_EDBL_SEQ_FMT) (        self.RawVec   ( iBox, iOpen )  , iBox = 1, self.nBoxStates )
          write(uid,FULL_EDBL_SEQ_FMT) (  dble( self.PsiMinus ( iBox, iOpen ) ), iBox = 1, self.nBoxStates )
          write(uid,FULL_EDBL_SEQ_FMT) ( aimag( self.PsiMinus ( iBox, iOpen ) ), iBox = 1, self.nBoxStates )
       enddo
    endif
  end subroutine ClassScatteringStateSaveToUnit

  function ClassScatteringStateLoadFromUnit(uid) result(self)
    integer                    , intent(in) :: uid
    class(ClassScatteringState), pointer    :: self
    character(len=16) :: form
    real(kind(1d0))   :: Energy
    integer           :: nBoxStates, nOpen 
    logical           :: solved
    integer           :: iOpen, iBox, iOpen2, iostat
    real(kind(1d0)), allocatable :: dvec(:)
    real(kind(1d0)), allocatable :: dmat(:,:)
    inquire(uid,form=form)
    if(trim(adjustl(form)).is."unformatted")then
       read(uid,iostat=iostat)     Energy, nBoxStates, nOpen, solved
       if(iostat/=0)then
          self => NULL()
          return
       else
          allocate(self)
       endif
       call Self.Init(Energy,nOpen,nBoxStates)
       self.solved=solved
!!$       read(uid) (   self.listOpen ( iOpen )        , iOpen  = 1, self.nOpen )
       read(uid) (   self.Name     ( iOpen )        , iOpen  = 1, self.nOpen )
       read(uid) (   self.lw       ( iOpen )        , iOpen  = 1, self.nOpen )
       read(uid) (   self.IonEnergy( iOpen )        , iOpen  = 1, self.nOpen )
       read(uid) (   self.PHI      ( iOpen )        , iOpen  = 1, self.nOpen )
       read(uid) ( ( self.A        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid) ( ( self.B        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid) ( ( self.Fitness  ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid) ( ( self.S        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid) ( ( self.K        ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid) ( ( self.EIGCH    ( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       allocate(dvec(self.nBoxStates))
       do iOpen = 1, self.nOpen
          read(uid) ( self.RawVec( iBox, iOpen ), iBox = 1, self.nBoxStates )
          read(uid) ( dvec( iBox ), iBox = 1, self.nBoxStates )
          self.PsiMinus(:,iOpen) = Z1*dvec
          read(uid) ( dvec( iBox ), iBox = 1, self.nBoxStates )
          self.PsiMinus(:,iOpen) = self.PsiMinus(:,iOpen) + Zi * dvec
       enddo
       deallocate(dvec)
    else
       read(uid,*,iostat=iostat) Energy, nBoxStates, nOpen, solved
       if(iostat/=0)then
          self => NULL()
          return
       else
          allocate(self)
       endif
       call Self.Init(Energy,nOpen,nBoxStates)
       self.solved=solved
       allocate(dmat(nOpen,nOpen))
!!$    read(uid,*) (   self.listOpen ( iOpen )         , iOpen  = 1, self.nOpen )
       read(uid,*) (   self.Name     ( iOpen )         , iOpen  = 1, self.nOpen )
       read(uid,*) (   self.lw       ( iOpen )         , iOpen  = 1, self.nOpen )
       read(uid,*) (   self.IonEnergy( iOpen )         , iOpen  = 1, self.nOpen )
       read(uid,*) (   self.PHI      ( iOpen )         , iOpen  = 1, self.nOpen )
       read(uid,*) ( ( self.A        ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid,*) ( ( self.B        ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid,*) ( ( self.Fitness  ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid,*) ( ( dmat          ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       self.S = Z1 * dmat
       read(uid,*) ( ( dmat          ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       self.S = self.S + Zi * dmat
       read(uid,*) ( ( self.K        ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       read(uid,*) ( ( self.EIGCH    ( iOpen2, iOpen ) , iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
       deallocate(dmat)
       allocate(dvec(self.nBoxStates))
       do iOpen = 1, self.nOpen
          read(uid,*) ( self.RawVec( iBox, iOpen ), iBox = 1, self.nBoxStates )
          read(uid,*) ( dvec( iBox ), iBox = 1, self.nBoxStates )
          self.PsiMinus(:,iOpen) = Z1*dvec
          read(uid,*) ( dvec( iBox ), iBox = 1, self.nBoxStates )
          self.PsiMinus(:,iOpen) = self.PsiMinus(:,iOpen) + Zi * dvec
       enddo
       deallocate(dvec)
    endif
  end function ClassScatteringStateLoadFromUnit

  !> Save on disk, in formatted form, the information about the scattering matrix
  !!
  !! COL 1 : Energy E
  !! COL 2 : Number of open channels at E (NOPEN)
  !! COL 3 : total phaseshift $ \phi_{tot}(E) = \sum_{\alpha\in Open} \phi_\alpha(E)$
  !! COL 4 : total error 
  !! COL 5 to (NOPEN+4) : partial phaseshifts $\phi_\alpha(E)$
  !! COL NOPEN + 5 to NOPEN + 4 + NOPEN^2 : column real eigenchannels C_{\alpha\beta}(E)
  !!
  !! Explanation: the scattering matrix \f$ S(E) \f$ is
  !! \f[
  !! S_{\alpha\beta}(E) = \sum_\gamma C_{\alpha\gamma}(E) \exp[2 i\phi_{\gamma}(E) ] C_{\beta\gamma}^\dagger(E)
  !! \f]
  !! with \f$ C_{\alpha\beta}^* = C_{\alpha\beta} \f$.
  !<
  subroutine ClassScatteringStatePrintToUnit(self,uid)
    implicit none
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: uid
    integer :: iOpen, iOpen2
    write( uid,"("//EDBL_FMT//",x,i0)", advance = "no" )     self.Energy,   self.nOpen
    write( uid, FULL_EDBL_SEQ_FMT     , advance = "no" ) sum(self.PHI), sum(self.Fitness)
    write( uid, FULL_EDBL_SEQ_FMT     , advance = "no" ) (   self.PHI  ( iOpen )        , iOpen  = 1, self.nOpen )
    write( uid, FULL_EDBL_SEQ_FMT     , advance = "no" ) ( ( self.EIGCH( iOpen2, iOpen ), iOpen2 = 1, self.nOpen ), iOpen = 1, self.nOpen )
    write( uid, * )
  end subroutine ClassScatteringStatePrintToUnit

  !> Save info about the channels open at a given energy
  subroutine ClassScatteringStatePrintChannelInfo(self,uid)
    implicit none
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: uid
    integer :: iOpen
    write( uid,"(i0,x,"//EDBL_FMT//")")  self.nOpen, self.Energy
    do iOpen = 1, self.nOpen
       write(uid,"(a,x,i0,x,f24.16)") self.Name( iOpen ), self.lw( iOpen ), self.IonEnergy( iOpen )
    enddo
    write( uid, * )
  end subroutine ClassScatteringStatePrintChannelInfo
  
  logical function ClassScatteringStatePrecedes( ss1, ss2 ) result( ss1_precedes_ss2 )
    implicit none
    class(ClassScatteringState), intent(in) :: ss1,ss2
    ss1_precedes_ss2 = ss1.GetEnergy() < ss2.GetEnergy()
  end function ClassScatteringStatePrecedes
  
  subroutine ClassScatteringStateAlignToPrev(self)
    implicit none
    class(ClassScatteringState), intent(inout) :: self
    integer :: n1, n2, mv(1), iOpen1, iOpen2, ich
    real(kind(1d0)), allocatable :: S(:)
    real(kind(1d0)) :: phi
    if(.not.(self%IsSolved()))call self%Solve()
    if(.not.(self%Prev%IsSolved()))call self%Prev%Solve()

    n1=self.Prev.nopen
    n2=self.nopen
    allocate(S(n2))
    do iOpen1 = 1, n1
       S=0.d0
       do iOpen2 = iOpen1, n2
          S(iOpen2)=abs(dot_product(self.prev.eigch(1:n1,iOpen1),self.eigch(1:n1,iOpen2)))
       enddo
       mv=maxloc(S)
       ich=mv(1)
       if(ich/=iOpen1)then
          S=self.eigch(:,ich)
          self.eigch(:,ich)=self.eigch(:,iOpen1)
          self.eigch(:,iOpen1)=S
          phi=self.phi(ich)
          self.phi(ich)=self.phi(iOpen1)
          self.phi(iOpen1)=phi
       endif
       do 
          if(self.phi(iOpen1)>self.prev.phi(iOpen1)+1.58)then
             self.phi(iOpen1)=self.phi(iOpen1)-PI
             cycle
          endif
          exit
       enddo
       do 
          if(self.phi(iOpen1)<self.prev.phi(iOpen1)-1.58)then
             self.phi(iOpen1)=self.phi(iOpen1)+PI
             cycle
          endif
          exit
       enddo
    enddo
  end subroutine ClassScatteringStateAlignToPrev

  subroutine ClassScatteringStateSetlw(self,ich,lw)
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: ich
    integer                    , intent(in)    :: lw
    self.lw(ich)=lw
  end subroutine ClassScatteringStateSetlw

  subroutine ClassScatteringStateSetThr(self,ich,Ethreshold)
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: ich
    real(kind(1d0))            , intent(in)    :: Ethreshold
    self.IonEnergy(ich)=Ethreshold
  end subroutine ClassScatteringStateSetThr

  subroutine ClassScatteringStateSetName(self,ich,Name)
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: ich
    character(len=*)           , intent(in)    :: Name
    self.Name(ich)=trim(adjustl(Name))
  end subroutine ClassScatteringStateSetName

  subroutine ClassScatteringStateSetCoef(self,ich,vec)
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: ich
    real(kind(1d0))            , intent(in)    :: vec(:)
    self.RawVec(:,ich)=vec
  end subroutine ClassScatteringStateSetCoef

  subroutine ClassScatteringStateSetAsymp(self,ich,A,B,F)
    class(ClassScatteringState), intent(inout) :: self
    integer                    , intent(in)    :: ich
    real(kind(1d0))            , intent(in)    :: A(:),B(:),F(:)
    self.A(:,ich)=A
    self.B(:,ich)=B
    self.Fitness(:,ich)=F
    !
  end subroutine ClassScatteringStateSetAsymp


  !.. Analyze the asymptotics of the channel                              
  !   The solutions will have the asymptotic form                         
  !   $$                                                                  
  !      Psi = F A + G B = (G+iF)(B-iA)/2 + (G-iF)(B+iA)/2 =
  !                      = W+ (B-iA) + W- (B+iA)
  !   $$                                   
  !   where
  !           W+ = (G+iF)/2      W- = (G-iF)/2
  !
  !   so the solution satisfying incoming boundary conditions is          
  !   $$                                                                  
  !      Psi- = Psi (B-iA)^{-1} = W+ + W- (B+iA)/(B-iA)          
  !   $$                                                                  
  !   The on-shell scattering matrix is then given by                     
  !   $$                                                                  
  !       S = (B-iA)/(B+iA)                                               
  !   $$                                                                  
  !   And the reaction matrix,                                            
  !   $$                                                                  
  !      K = i/pi (S-1)/(S+1)                                             
  !      S = (1-i\pi K)/(1+i\pi K)                                        
  !   $$                                                                  
  !   is                                                                  
  !   $$                                                                  
  !      pi K = A / B                                                     
  !   $$                                                                  
  !   The condition for the results to be consistent, therefore, is       
  !   that the matrix A/B is Hermitean (and hence, the scattering matrix  
  !   is unitary).                                                        
  !..                                                                     
  subroutine ClassScatteringStateSolve(self)
    !
    use ModuleDiagonalize
    !
    implicit none
    !
    class(ClassScatteringState), intent(inout) :: self
    real   (kind(1d0)), allocatable :: B_1(:,:), Kmat(:,:), work(:)
    integer           , allocatable :: ipiv(:)
    complex(kind(1d0)), allocatable :: zmat(:,:), B_iA(:,:), zwork(:), zRaw(:,:)
    real   (kind(1d0)) ::  OPTIMAL_WORKSPACE_SIZE(1), normK, asymK
    complex(kind(1d0)) :: zOPTIMAL_WORKSPACE_SIZE(1)
    integer, parameter :: WORKSPACE_QUERY = -1
    integer :: i, nopen, lwork, info, nb, no
    !
    nopen=self.nopen
    allocate(B_1,source=self.B)
    allocate(ipiv(nopen))
    !
    !.. Compute the reaction matrix K = A/B /pi                                               
    !
    !.. 1. Compute LU factorization of B
    call dgetrf(nopen,nopen,B_1,nopen,ipiv,info)
    if(info/=0) write(ERROR_UNIT,"(a,i0)")"dgetrf fail, info = ",info
    !
    !.. 2. Compute B^-1
    call dgetri(nopen,B_1,nopen,ipiv,OPTIMAL_WORKSPACE_SIZE,WORKSPACE_QUERY,info)
    lwork=int(OPTIMAL_WORKSPACE_SIZE(1)+0.1)
    allocate(work(lwork))
    call dgetri(nopen,B_1,nopen,ipiv,work,lwork,info)
    deallocate(work)
    if(info/=0) write(ERROR_UNIT,"(a,i0)")"dgetri fail, info = ",info
    !
    !.. Compute the reaction matrix K = A / B / pi
    self.K=matmul(self.A,B_1)/PI
    deallocate(B_1)
    !
    normK=sum(abs(self.K))
    asymK=sum(abs(self.K-transpose(self.K)))
    !
    !.. symmetrizes the matrix, which should have been symmetric in the first place
    !self.K=0.5d0*(self.K+transpose(self.K))
    !
    write(*,*) "normK, asymK", normK, asymK
    !
    !.. Diagonalize K 
    allocate(Kmat,source=self.K)
    call Short_Diag(nopen,Kmat,self.PHI)
    self.EIGCH=Kmat
    deallocate(Kmat)
    self.PHI=-atan(PI*self.PHI)
    allocate(zmat,source=(1.d0,0.d0)*self.EIGCH)
    do i=1,size(zmat,1)
       zmat(:,i) = zmat(:,i) * exp( (0.d0,2.d0) * self.PHI(i) )
    enddo
    self.S=matmul(zmat,transpose(self.EIGCH))
    deallocate(zmat)
    !
    !.. Compute the scattering vectors Psi-
    !
    !.. Compute ( B - iA )^-1
    !
    !.. 1. Compute LU factorization of ( B - i A )
    allocate(B_iA(nopen,nopen))
    B_iA = Z1 * self.B - Zi * self.A
    call zgetrf(nopen,nopen,B_iA,nopen,ipiv,info)
    if(info/=0) write(ERROR_UNIT,"(a,i0)")"dgetrf fail, info = ",info
    !
    !.. 2. Compute 1 / ( B - iA )
    call zgetri(nopen,B_iA,nopen,ipiv,ZOPTIMAL_WORKSPACE_SIZE,WORKSPACE_QUERY,info)
    lwork=int(dble(ZOPTIMAL_WORKSPACE_SIZE(1))+0.1)
    allocate(zwork(lwork))
    call zgetri(nopen,B_iA,nopen,ipiv,zwork,lwork,info)
    if(info/=0) write(ERROR_UNIT,"(a,i0)")"dgetri fail, info = ",info
    deallocate(zwork,ipiv)
    !
    !.. 3. Compute Psi- = Psi / (B-iA)
    nb = self.nBoxStates
    no = self.nopen
    allocate(zRaw(nb,no))
    zRaw=Z1*self.RawVec
    call ZGEMM("N","N",nb,no,no,Z1,zRaw,nb,B_iA,no,Z0,self.PsiMinus,nb)
    deallocate(zRaw)
    self.solved=.TRUE.
  end subroutine ClassScatteringStateSolve

  subroutine ClassScatteringStateTransform( self, dMat )
    implicit none
    class( ClassScatteringState ), intent(inout) :: self
    real(kind(1d0))              , intent(in)    :: dMat(:,:)
    real(kind(1d0)), allocatable, save :: &
         dRaw   (:,:), &
         dPsiRe1(:,:), dPsiRe2(:,:), &
         dPsiIm1(:,:), dPsiIm2(:,:)
    integer :: nb, no
    nb = self.nBoxStates
    no = self.nopen
    call realloc( dRaw   , nb, no )
    call realloc( dPsiRe1, nb, no )
    call realloc( dPsiIm1, nb, no )
    call realloc( dPsiRe2, nb, no )
    call realloc( dPsiIm2, nb, no )
    dRaw    =        self.RawVec
    dPsiRe1 =  dble( self.PsiMinus )
    dPsiIm1 = aimag( self.PsiMinus )
    dRaw    = 0.d0
    dPsiRe2 = 0.d0
    dPsiIm2 = 0.d0
    call DGEMM( "N", "N", nb, no, nb, 1.d0, dMat, nb, dRaw,    nb, 0.d0, self.RawVec, nb )
    call DGEMM( "N", "N", nb, no, nb, 1.d0, dMat, nb, dPsiRe1, nb, 0.d0, dPsiRe2    , nb )
    call DGEMM( "N", "N", nb, no, nb, 1.d0, dMat, nb, dPsiIm1, nb, 0.d0, dPsiIm2    , nb )
    self.PsiMinus = Z1 * dPsiRe2 + Zi * dPsiIm2
  end subroutine ClassScatteringStateTransform

  subroutine ClassScatteringStateProjectOnBra( self, zBra, zOpen )
    implicit none
    class( ClassScatteringState ), intent(inout) :: self
    complex(kind(1d0))           , intent(in)    :: zBra(:)
    complex(kind(1d0))           , intent(out)   :: zOpen(:)
    integer :: iOpen, iBox
    zOpen=Z0
    do iOpen=1,self.nopen
       do iBox=1,self.nBoxStates
          !.. Assumes zBra does not need to be conjugated (i.e., it is already in c^\dagger form)
          zOpen( iOpen ) = zOpen( iOpen ) + zBra( iBox ) * self%PsiMinus( iBox, iOpen )
       enddo
    enddo
  end subroutine ClassScatteringStateProjectOnBra

  subroutine ClassScatteringStateGetPsiMinus( self, ich, zKet )
    implicit none
    class( ClassScatteringState ), intent(inout) :: self
    integer                      , intent(in)    :: ich
    complex(kind(1d0))           , intent(out)   :: zKet(:)
    integer :: iBox
    zKet=Z0
    if(ich<=0.or.ich>self.nopen)return
    zKet(:) = self%PsiMinus(:,ich)
  end subroutine ClassScatteringStateGetPsiMinus

  subroutine ClassScatteringStateSetPsiMinus( self, ich, zKet )
    implicit none
    class( ClassScatteringState ), intent(inout) :: self
    integer                      , intent(in)    :: ich
    complex(kind(1d0))           , intent(in)    :: zKet(:)
    integer :: iBox
    if(ich<=0.or.ich>self.nopen)return
    self%PsiMinus(:,ich) = zKet(:)
  end subroutine ClassScatteringStateSetPsiMinus

  logical function ClassScatteringStateIsSolved(self) result(solved)
    class(ClassScatteringState), intent(in) :: self
    solved=self.solved
  end function ClassScatteringStateIsSolved

  real(kind(1d0)) function ClassScatteringStateGetTotPh(self) result(totPh)
    class(ClassScatteringState), intent(in) :: self
    totPh=sum(self.PHI)
  end function ClassScatteringStateGetTotPh

  real(kind(1d0)) function ClassScatteringStateGetEnergy(self) result(E)
    class(ClassScatteringState), intent(in) :: self
    E=self.Energy
  end function ClassScatteringStateGetEnergy

  subroutine ClassScatteringStateFree( Self )
    class( ClassScatteringState ), intent(inout) :: Self
    !
    Self%Prev       => NULL()
    Self%Next       => NULL()
    self.nOpen      =  0
    self.solved     = .FALSE.
    self.Energy     =  0.d0
    self.nBoxStates =  0
    !
!!$    if(allocated(self.listOpen )) deallocate(self.listOpen)
    if(allocated(self.lw       )) deallocate(self.lw)
    if(allocated(self.IonEnergy)) deallocate(self.IonEnergy)
    if(allocated(self.RawVec   )) deallocate(self.RawVec)
    if(allocated(self.PsiMinus )) deallocate(self.PsiMinus)
    if(allocated(self.A        )) deallocate(self.A)
    if(allocated(self.B        )) deallocate(self.B)
    if(allocated(self.Fitness  )) deallocate(self.Fitness)
    if(allocated(self.S        )) deallocate(self.S)
    if(allocated(self.K        )) deallocate(self.K)
    if(allocated(self.PHI      )) deallocate(self.PHI)
    if(allocated(self.EIGCH    )) deallocate(self.EIGCH)
    !
  end subroutine ClassScatteringStateFree

  subroutine ClassScatteringStateFinalize( Self )
    type( ClassScatteringState ), intent(inout) :: Self
    call Self.free()
  end subroutine ClassScatteringStateFinalize

end module ModuleScatteringStates
