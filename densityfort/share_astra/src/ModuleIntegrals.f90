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
module ModuleIntegrals

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleSystemUtils
  use ModuleString
  use ModuleErrorHandling
  use ModuleMatrix
  use ModuleBspline
  use ModuleGroups
  use ModuleXlm
  use ModuleMolecularGeometry

  implicit none

  !..  The term MO here refers to all molecular orbitals: inactive,
  !    active, and virtual orbitals. In terms of the channel
  !    distinction in the close-coupling calculation, the most
  !    relevant aspect is whether an orbital is inactive (and
  !    assumed to be always doubly occupied in the ions), active,
  !    i.e., with variable occupancy in the ions, or virtual, i.e.,
  !    never occupied in the ion. In terms of many-body formulas,
  !    therefore, the virtual molecular orbitals are in the same
  !    category as the hybrid orbitals. On the other hand, the
  !    boundary between inactive, active, and virtual orbitals is
  !    movable without affecting the orbital identity or integrals.
  !    For this reason, it is best to keep all the MOs together in
  !    the integrals. But they should be referenced as such, MO,
  !    rather than with the name of a specific channel. This
  !    conflation of domains can only cause confusion.
  !    The labeling in the integrals should be MO, HY, BS.
  !..
  
  private


  integer, parameter, public :: I1B_N_TYPES =11 
  integer, parameter, public :: I1B_OVERLAP = 1
  integer, parameter, public :: I1B_KINETIC = 2
  integer, parameter, public :: I1B_NUCLATT = 3
  integer, parameter, public :: I1B_HAMILTO = 4
  integer, parameter, public :: I1B_COORD_X = 5
  integer, parameter, public :: I1B_COORD_Y = 6
  integer, parameter, public :: I1B_COORD_Z = 7
  integer, parameter, public :: I1B_NABLA_X = 8
  integer, parameter, public :: I1B_NABLA_Y = 9
  integer, parameter, public :: I1B_NABLA_Z =10
  integer, parameter, public :: I1B_MULTIPO =11

  character(len=*), parameter :: INTEGRAL_SUBDIR            = "integrals/"
  character(len=*), parameter :: BIELECTRONIC_INTEGRAL_FILE = "biel"
  character(len=*), parameter :: NUCLEAR_MOMENT_FILE        = "nuclear_moments"
  character(len=*), parameter :: LOC_ORBITAL_MOMENT_FILE    = "locorbs_moments"
  character(len=*), parameter :: BSPLINE_INTEGRAL_FILE      = "bspline_ints"

  character(len=*), parameter :: NUCLEAR_MOMENT_NAME = "NuclearMoments"
  integer, parameter :: I1B_ID_LEN = 16
  character(len=I1B_ID_LEN), public, parameter :: I2B_ID_NAME = "bielectronic"
  character(len=I1B_ID_LEN), public, parameter :: I1B_ID_LIST(I1B_N_TYPES) = [ &
       "Overlap", &
       "Kinetic", &
       "NuclAtt", &
       "Hamilto", &
       "Coord_X", &
       "Coord_Y", &
       "Coord_Z", &
       "Nabla_X", &
       "Nabla_Y", &
       "Nabla_Z", &
       "MP_    "]

  !!--------------------

!!$  integer :: MODULE_INTEGRALS_LMAX
!!$
  type ClassMoments
     private
     !.. Definition:
     !
     !   M(\alpha,\beta,lm)%A(i,j) = 
     !   \frac{4\pi}{2l+1} \int d^3r r^\ell X_{lm}(\hat{r}) 
     !                               \varphi^*_{\alpha i}(\vec{r}) \varphi_{\beta j}(\vec{r})
     !
     !   where lm is the index for the (l,m) pair
     logical :: INITIALIZED = .FALSE.
     !type(ClassIrrep), pointer     :: OpIrrep(:)
     integer                       :: lmax ! ??
     integer         , allocatable :: nMOv(:)

     !.. First  index: irrep of the Bra orbital
     !   Second index: irrep of the Ket orbital
     !   Third  index: pair (lm) of the multipole r^\ell X_{\ell m}(\hat{r})
     type(ClassMatrix), pointer    :: MO_MO(:,:,:)

   contains

     procedure, public :: init => ClassMomentsInit
     procedure, public :: free => ClassMomentsFree
     procedure, public :: set  => ClassMomentsSet
     procedure, public :: Write=> ClassMomentsWriteToFile
     procedure, public :: Read => ClassMomentsReadFromFile
!!$     procedure, public :: get  => ClassMoments_get !.. Returns pointer to a Moment with given l and m

  end type ClassMoments

  type ClassNuclearMoments
     !private
     !.. Definition:
     !
     !   M(l,m) = \frac{4\pi}{2l+1} \sum_{\alpha} Q_\alpha R_\alpha^l X_{lm}(\hat{R}_\alpha)
     !
     !   where \alpha, \vec{R}_\alpha, and Q_\alpha indicate a nucleus 
     !   index, position, and charge, respectively.
     !
     !   They are all totally symmetric, by definition of molecular symmetry
     !   which means that only the totally symmetric Xlm have non-zero moment.
     !..
     logical :: INITIALIZED = .FALSE.
     integer                       :: lmax
     character(len=:), allocatable :: Name
     real(kind(1d0)) , allocatable :: Moments(:,:)
   contains

     generic, public :: init       => ClassNuclearMomentsInit
     generic, public :: free       => ClassNuclearMomentsFree
     generic, public :: Set        => ClassNuclearMomentsSet
     generic, public :: Write      => ClassNuclearMomentsWriteToFile
     generic, public :: Read       => ClassNuclearMomentsReadFromFile
     generic, public :: GetMoments => ClassNuclearMomentsGetMoments
     
     procedure, private :: ClassNuclearMomentsInit
     procedure, private :: ClassNuclearMomentsFree
     procedure, private :: ClassNuclearMomentsSet
     procedure, private :: ClassNuclearMomentsReadFromFile
     procedure, private :: ClassNuclearMomentsWriteToFile
     procedure, private :: ClassNuclearMomentsGetMoments
     
  end type ClassNuclearMoments


  type ClassBsBlocks
     private
     !.. This module stores the ibtegrals between bsplines used 
     !   both in the VIC - VIC and in the VEC - VEC integrals.
     logical :: INITIALIZED = .FALSE.
     integer                       :: lmax
     integer                       :: NBsplines,BsOrder
     real(kind(1d0)) , allocatable :: Overlap(:,:)
     real(kind(1d0)) , allocatable :: SecondDeriv(:,:)
     real(kind(1d0)) , allocatable :: Centrifugal(:,:) !*** Could be replaced by Multipoles(:,:,1)
     character(len=:), allocatable :: Name
     !
     !.. Integral of r^{l+1}
     !   Multipoles(i,j,l) = N_i N_j \int dr B_i(r) B_j(r) / r^{l+1}
     !   Notice that, therefore, the third index starts at 0
     real(kind(1d0)) , allocatable :: Multipoles(:,:,:)

   contains

     generic, public :: Init  => ClassBsBlocksInit
     generic, public :: Free  => ClassBsBlocksFree
     generic, public :: Set   => ClassBsBlocksSet
     generic, public :: Write => ClassBsBlocksWriteToFile
     generic, public :: Read  => ClassBsBlocksReadFromFile
     generic, public :: GetOverlap     => ClassBsBlocksGetOverlap
     generic, public :: GetSecondDeriv    => ClassBsBlocksGetSecondDeriv
     generic, public :: GetCentrifugal => ClassBsBlocksGetCentrifugal
     generic, public :: GetMultipoles  => ClassBsBlocksGetMultipoles

     procedure, private :: ClassBsBlocksInit
     procedure, private :: ClassBsBlocksFree
     procedure, private :: ClassBsBlocksSet
     procedure, private :: ClassBsBlocksReadFromFile
     procedure, private :: ClassBsBlocksWriteToFile
     procedure, private :: ClassBsBlocksGetOverlap
     procedure, private :: ClassBsBlocksGetSecondDeriv
     procedure, private :: ClassBsBlocksGetCentrifugal
     procedure, private :: ClassBsBlocksGetMultipoles
  end type ClassBsBlocks


  type ClassIntegral1Body

     private
     logical :: INITIALIZED = .FALSE.
     character(len=:), allocatable :: Name
     type(ClassIrrep), pointer     :: OpIrrep
     integer                       :: lmax ! ??
     integer         , allocatable :: nMOv(:)
     integer         , allocatable :: nHYv(:)
     integer                       :: nBS !.. "Radial" index of the Virtual External Orbitals 

     !.. First  index: irrep of the Bra orbital
     !   Second index: irrep of the Ket orbital
     type(ClassMatrix), allocatable   :: MO_MO(:,:)
     type(ClassMatrix), allocatable   :: MO_HY(:,:)
     type(ClassMatrix), allocatable   :: HY_HY(:,:)

     !.. First  index: irrep of Bra orbital
     !   Second index: lm pair index of the Ket
     !                 note: whreas all "lm" pairs are present,
     !                       only the elements with lm compatible
     !                       with the Bra and Op irrep are allocated
     type(ClassMatrix), allocatable    :: HY_BS(:,:)

     !.. First  index: lm pair index of the Bra
     !.. Second index: lm pair index of the Ket
     type(ClassMatrix), allocatable    :: BS_BS(:,:)

   contains

     generic  , public  :: init => ClassIntegral1BodyInit
     generic  , public  :: free => ClassIntegral1BodyFree
     generic  , public  :: Set  => ClassIntegral1Body_Set_SubMatrix 
     generic  , public  :: Write  => ClassIntegral1Body_WriteToFile
     generic  , public  :: Read   => ClassIntegral1Body_ReadFromFile

     generic  , public  :: Set_MO_MO => &
          ClassIntegral1Body_Set_MO_MO_Element, &
          ClassIntegral1Body_Set_MO_MO_Matrix, &
          ClassIntegral1Body_Set_MO_MO_SubMatrix
     
     procedure, private :: ClassIntegral1BodyInit
     procedure, private :: ClassIntegral1BodyFree
     procedure, private :: ClassIntegral1Body_Set_SubMatrix
     procedure, private :: ClassIntegral1Body_WriteToFile
     procedure, private :: ClassIntegral1Body_ReadFromFile

     procedure, private :: ClassIntegral1Body_Set_MO_MO_Element
     procedure, private :: ClassIntegral1Body_Set_MO_MO_Matrix
     procedure, private :: ClassIntegral1Body_Set_MO_MO_SubMatrix


     
  end type ClassIntegral1Body

  type ClassIntegral2Body

     private
     logical :: INITIALIZED = .FALSE.
     type(ClassIrrep), pointer     :: OpIrrep
     character(len=:), allocatable :: Name
     integer                       :: lmax ! ??
     integer         , allocatable :: nMOv(:)
     integer         , allocatable :: nHYv(:)

     !.. Bielectronic integrals, indexed according to the density notation
     !   \[
     !   [ f_1 f_2 | f_3 f_4 ] = \int d^3 r_1 \int d^3 r_2 
     !                           \frac{\rho_{12}(\vec{r}_1 \rho_{34}(\vec{r}_2}{r_{12}}
     !   \]
     !   where \$ \rho_{ij}(\vec{r}) = f_i(\vec{r}) f_j(\vec{r}) \$
     ! 
     !.. First  index: irrep of the first  orbital of the left  density
     !   Second index: irrep of the second orbital of the left  density
     !   Third  index: irrep of the first  orbital of the right density
     !   Fourth index: irrep of the second orbital of the right density
     !
     !   Notice the eightfold symmetry Z2^3 = Z2 x Z2 x Z2 of the integrals:
     !
     !   [ 1 2 | 3 4 ] = [ 2 1 | 3 4 ]   Z2 
     !   [ 1 2 | 3 4 ] = [ 1 2 | 4 3 ]   Z2
     !   [ 1 2 | 3 4 ] = [ 4 3 | 1 2 ]   Z2
     !   
     !   In principle, therefore, for orbitals all in the same group (MO, HY, or BS)
     !   the irreps can be ordered as 
     !   i1 >= i2
     !   13 >= i4
     !   i1 >= i3 and if i1 = i3, then i2 >= i4
     !
     type(ClassMatrix4D), pointer    :: LL_LL(:,:,:,:)
     type(ClassMatrix4D), pointer    :: LL_LS(:,:,:,:)
     type(ClassMatrix4D), pointer    :: LL_SS(:,:,:,:)
     type(ClassMatrix4D), pointer    :: LS_LS(:,:,:,:)

     !.. Note: the integrals LL_SP and LL_PP are exactly decomposed in terms of multipoles,
     !         and hence there is no reason to store them with four indexes. Instead,
     !         it is much better to contrac the multipoles with the TDM and evaluate
     !         the bielectronic contribution by multiplying those transition multipoles
     !         with the 1B integrals of the corresponding 1/r^{l+1} between SP or PP.
     !         For single ionization, no other integrals are needed
     !..

   contains

     generic  , public  :: init => ClassIntegral2BodyInit
     generic  , public  :: free => ClassIntegral2BodyFree
     generic  , public  :: Set  => ClassIntegral2Body_Set4DMat 
     generic  , public  :: CheckIrrepOrder => ClassIntegral2BodyCheckIrrepOrder
     generic  , public  :: Write   => ClassIntegral2Body_WriteToFile
     generic  , public  :: Read    => ClassIntegral2Body_ReadFromFile


     procedure, private :: ClassIntegral2BodyInit
     procedure, private :: ClassIntegral2BodyFree
     procedure, private :: ClassIntegral2Body_Set4DMat
     procedure, public  :: ClassIntegral2BodyCheckIrrepOrder
     procedure, private :: ClassIntegral2Body_WriteToFile
     procedure, private :: ClassIntegral2Body_ReadFromFile

  end type ClassIntegral2Body

  !.. Singleton for the module
  type, private :: ClassIntegralSingleton

     private
     character(len=:), allocatable :: store
     integer                       :: lmax = 0
     integer         , allocatable :: nInactive(:)
     !*** IT WOULD BE BEST NOT TO MAKE INT_1B PUBLIC AND INSTEAD
     !    MAKE IT A VECTOR OF POINTER, WITH A FUNCTION THAT RETURNS
     !    THE POINTER TO THE CLASSINTEGRAL1BODY OBJECT THAT CORRESPONDS
     !    TO A CERTAIN INTEGRAL THAT IS REQUESTED.
     type(ClassIntegral1Body), public, allocatable :: Int_1B(:)
     type(ClassIntegral2Body), public              :: Int_2B
     !.. Evaluates the Nuclear moments
     type(ClassNuclearMoments), public             :: NucMoments
     type(ClassMoments)       , public             :: LocMoments
     !.. Stores the Bsplines blocks
     type(ClassBsBlocks)     , public              :: BsBlocks
     

   contains

     procedure, public  :: Init              => ClassIntegralSingleton_Init
     procedure, public  :: Get_Multipole_ID  => ClassIntegralSingleton_Get_Multipole_ID
     procedure, public  :: PairTolm          => ClassIntegralSingleton_PairTolm
     procedure, public  :: lmPairIndex       => ClassIntegralSingleton_lmPairIndex
     procedure, public  :: Free              => ClassIntegralSingleton_Free
     procedure, public  :: SetStorage        => ClassIntegralSingleton_SetStorage
     procedure, public  :: Setlmax           => ClassIntegralSingleton_Setlmax
     procedure, public  :: SetnInactive      => ClassIntegralSingleton_SetnInactive

     procedure, public  :: GetnMOv           => ClassIntegralSingleton_GetnMOv
     procedure, public  :: GetnHYv           => ClassIntegralSingleton_GetnHYv
     procedure, public  :: GetnBS            => ClassIntegralSingleton_GetnBS
     
     procedure, public  :: WriteToFile       => ClassIntegralSingleton_WriteToFile
     procedure, public  :: ReadFromFile      => ClassIntegralSingleton_ReadFromFile
     procedure, public  :: OneBIsInitialized => ClassIntegralSingleton_OneBIsInitialized
     procedure, public  :: Get1B             => ClassIntegralSingleton_Get1B
     procedure, public  :: GetLocMoments     => ClassIntegralSingleton_GetLocMoments
     procedure, public  :: LocMomentsIsInit  => ClassIntegralSingleton_LocMomentsIsInit
     procedure, public  :: BielIsInitialized => ClassIntegralSingleton_BielIsInitialized
     generic  , public  :: GetBiel           => ClassIntegralSingleton_GetBiel, &
          ClassIntegralSingleton_GetBiel_blk
     procedure, public  :: Condition_beS     => ClassIntegralSingleton_Condition_beS
     
     procedure, public  :: ConvertHinHtilde  => ClassIntegralSingleton_ConvertHinHtilde
     procedure, private :: ClassIntegralSingleton_ConvertHinHtilde
     procedure, private :: ClassIntegralSingleton_GetBiel
     procedure, private :: ClassIntegralSingleton_GetBiel_blk
     
  end type ClassIntegralSingleton
  type(ClassIntegralSingleton), public :: GlobalIntegral

  type, private :: DVecContainer
     real(kind(1d0)), allocatable :: v(:)
  end type DVecContainer
  
contains

  Pure DoublePrecision function Unity(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = 1.d0
  end function Unity

  Pure DoublePrecision function Power(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = x**parvec(1)
!!$    y = exp(log(x)*parvec(1))
  end function Power

  subroutine ClassMomentsFree( self )
    class(ClassMoments), intent(inout) :: self
    integer :: lmax
    integer :: iIrr1, iIrr2, l, m, ilm
    lmax = self%lmax
    if(allocated(self%nMOv))deallocate(self%nMOv)
    !.. Cycle over irreps Bra
    do l = 0, lmax
       do m = -l, l
          ilm = lmPairIndex(l,m)
          do iIrr1 = 1, GlobalGroup%GetnIrreps()
             !.. Cycle over irreps Ket
             do iIrr2 = 1, GlobalGroup%GetnIrreps()
                call self%MO_MO(iIrr1,iIrr2,ilm)%Free()
             enddo
          enddo
       end do
    end do
    deallocate(self%MO_MO)
    self%lmax=-1
    self%INITIALIZED=.FALSE.
  end subroutine ClassMomentsFree

  subroutine ClassMomentsInit( self, lmax, nMOv )
    use ModuleXlm
    class(ClassMoments)    , intent(inout) :: self
    integer                , intent(in)    :: lmax
    integer                , intent(in)    :: nMOv(:)
    !type(ClassIrrep), pointer, intent(in ) :: OpIrrep

    integer, allocatable :: mlist(:)
    type(ClassIrrep), pointer :: IrrPtr1, IrrPtr2, OpIrrep
    type(ClassIrrep), pointer :: IrrList(:)
    type(ClassXlmSymmetricSet), pointer :: XlmSymSet
    ! class(ClassXlm),  pointer :: Xlm
    !  class(ClassXlm) :: Xlm
    integer                   :: l, m, im, iIrr1, iIrr2, ilm, nlm, nIrreps

    nIrreps =  GlobalGroup%GetnIrreps()
    IrrList => GlobalGroup%GetIrrepList()
    allocate(self%nMOv(nIrreps))

    self%nMOv = nMOv
    self%lmax  = lmax
    nlm = lmPairIndex(2*lmax,2*lmax)

    allocate(self%MO_MO(nIrreps,nIrreps,nlm))

    !.. Cycle over the Bra irreps
    do iIrr1 = 1, nIrreps
       IrrPtr1 => IrrList( iIrr1 )
       !
       do iIrr2 = 1, nIrreps
          IrrPtr2 => IrrList( iIrr2 )
          !
          OpIrrep   => IrrPtr1 * IrrPtr2
          XlmSymSet => GlobalXlmSet%GetSymSet( OpIrrep )   
          do l = 0, lmax
             call XlmSymSet%GetMList( l, mlist )
             if(.not.allocated(mlist))cycle
             !
             do im = 1, size(mlist)
                m = mlist(im)
                ilm = lmPairIndex(l,m)
                call self%MO_MO(iIrr1,iIrr2,ilm)%InitFull( nMOv (iIrr1), nMOv (iIrr2) )
             end do
             !
          end do
       end do
    end do
    self%INITIALIZED=.TRUE.

  end subroutine ClassMomentsInit

  subroutine ClassMomentsSet( self, iIrr1, iIrr2, ilm, dMat )
    use ModuleXlm
    class(ClassMoments)    , intent(inout) :: self
    integer                , intent(in)    :: iIrr1, iIrr2, ilm
    real(kind(1d0))        , intent(in)    :: dMat(:,:)
    self%MO_MO(iIrr1,iIrr2,ilm) = dMat
  end subroutine ClassMomentsSet


  subroutine ClassMomentsWriteToFile( self, FileName )
    class(ClassMoments), intent(in) :: self
    character(len=*)   , intent(in) :: FileName

    integer :: uid, ilm, iostat, l, m
    integer :: nIrreps, iIrr, iIrr1, iIrr2
    character(len=500) :: iomsg

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="unknown", &
         action  ="write"  , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in Writing "//trim(iomsg)
       stop
    endif
    nIrreps = GlobalGroup%GetnIrreps()

    write(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif

    write(uid) self%lmax
    write(uid) (self%nMOv(iIrr),iIrr=1,nIrreps)
    
    do iIrr1 = 1, nIrreps
       do iIrr2 = 1, nIrreps
          do l = 0, 2*self%lmax
             do m = -l, l
                ilm = lmPairIndex(l,m)
                if(self%MO_MO(iIrr1,iIrr2,ilm)%IsInitialized())then
                   write(uid) iIrr1, iIrr2, ilm
                   call self%MO_MO(iIrr1,iIrr2,ilm)%Write(uid)
                   !write(*,*) self%MO_MO(iIrr1,iIrr2,ilm)%NRows(),self%MO_MO(iIrr1,iIrr2,ilm)%NColumns()
                   !write(*,*) iIrr1, iIrr2, ilm
                   !write(*,*) self%MO_MO(iIrr1,iIrr2,ilm)%A
                   
                endif
             enddo
          enddo
       enddo
    enddo
    write(uid) -1, -1, -1
    close(uid)
    !pause

  end subroutine ClassMomentsWriteToFile


  subroutine ClassMomentsReadFromFile( self, FileName )
    class(ClassMoments), intent(inout) :: self
    character(len=*)   , intent(in)    :: FileName

    integer :: uid, ilm, iostat
    integer :: nIrreps, iIrr, iIrr1, iIrr2
    character(len=500) :: iomsg

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="old"    , &
         action  ="read"   , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in reading "//trim(iomsg)
       stop
    endif
    read(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif
    
    nIrreps = GlobalGroup%GetnIrreps()
    allocate(self%nMOv(nIrreps))
    read(uid) self%lmax
    read(uid) (self%nMOv(iIrr),iIrr=1,nIrreps)
    allocate(self%MO_MO(nIrreps,nIrreps,lmPairIndex(2*self%lmax,2*self%lmax)))
    do 
       read(uid) iIrr1, iIrr2, ilm
       !pause
       if(iIrr1<0)exit
       call self%MO_MO(iIrr1, iIrr2,ilm)%read(uid)
       !write(*,*) self%MO_MO(iIrr1, iIrr2,ilm)%NRows(),self%MO_MO(iIrr1, iIrr2,ilm)%NColumns()
       !write(*,*) iIrr1, iIrr2, ilm
       !write(*,*) self%MO_MO(iIrr1, iIrr2,ilm)%isinitialized()
       !write(*,*) self%MO_MO(iIrr1, iIrr2,ilm)%A
    enddo

    close(uid)
    !stop
  end subroutine ClassMomentsReadFromFile


  subroutine ClassBsBlocksFree( self )
    class(ClassBsBlocks)  , intent(inout) :: self
    if(allocated(self%Overlap))deallocate(self%Overlap)
    if(allocated(self%SecondDeriv))deallocate(self%SecondDeriv)
    if(allocated(self%Centrifugal))deallocate(self%Centrifugal)
    if(allocated(self%Multipoles))deallocate(self%Multipoles)
    self%lmax=-1
    self%NBsplines=-1
    self%INITIALIZED=.FALSE.
  end subroutine ClassBsBlocksFree

  subroutine ClassBsBlocksInit( self, lmax, NBsplines, BsOrder )
    class(ClassBsBlocks), intent(inout) :: self
    integer             , intent(in)    :: lmax
    integer             , intent(in)    :: NBsplines
    integer             , intent(in)    :: BsOrder
    self%BsOrder   = BsOrder
    self%lmax      = lmax
    self%NBsplines = NBsplines
    allocate(self%Overlap(NBsplines,NBsplines))
    allocate(self%SecondDeriv(NBsplines,NBsplines))
    allocate(self%Centrifugal(NBsplines,NBsplines))
    allocate(self%Multipoles(NBsplines,NBsplines,0: 2 * lmax))

    self%Overlap    =0.d0
    self%SecondDeriv=0.d0
    self%Centrifugal=0.d0
    self%Multipoles =0.d0

    self%INITIALIZED=.TRUE.
  end subroutine ClassBsBlocksInit

  subroutine ClassBsBlocksSet( self , BsSet )
    class(ClassBsBlocks), intent(inout) :: self
    class(ClassBSpline)  , intent(in)    :: BsSet

    procedure(D2DFun), pointer :: fPtr
    real(kind(1d0))            :: parvec(2) = 0.d0
    integer                    :: l, lmax, BsOrder, NBsplines
    integer                    :: in1, in2
    real(kind(1d0))            :: norm

    if(.not.self%INITIALIZED)then
       write(*,*) "BsBlocksSet: BsBlocks is not initialized"
       stop
    endif
    
   BsOrder   = self%BsOrder
   lmax      = self%lmax
   NBsplines = self%NBsplines

    fPtr => Power
    do in1 = 1, NBsplines
       do in2 = 1, NBsplines
          if(abs(in1-in2)>=BsOrder)cycle
          norm = BsSet%GetNormFactor(in1) * BsSet%GetNormFactor(in2)

          ! parvec(1)=0.d0
          ! GlobalIntegral%BsBlocks%Overlap    (in1, in2) = norm * BsSet%Integral(fPtr,in1,in2,0,0,parvec=parvec)
          ! GlobalIntegral%BsBlocks%SecondDeriv(in1, in2) = norm * BsSet%Integral(fPtr,in1,in2,0,2,parvec=parvec)

          ! parvec(1)=-2.d0
          ! GlobalIntegral%BsBlocks%Centrifugal(in1, in2) = norm * BsSet%Integral(fPtr,in1,in2,0,0,parvec=parvec)

          ! do l = 0, lmax
          !    parvec(1)=-1.d0*( l + 1 )
          !    GlobalIntegral%BsBlocks%Multipoles(in1, in2, l) = norm * BsSet%Integral(fPtr,in1,in2,0,0,parvec=parvec)
          ! enddo
          
          parvec(1)=0.d0
          self%Overlap    (in1, in2) = norm * BsSet%Integral(fPtr,in1,in2,0,0,parvec=parvec)
          self%SecondDeriv(in1, in2) = norm * BsSet%Integral(fPtr,in1,in2,0,2,parvec=parvec)

          parvec(1)=-2.d0
          self%Centrifugal(in1, in2) = norm * BsSet%Integral(fPtr,in1,in2,0,0,parvec=parvec)

          do l = 0, 2*lmax
             parvec(1)=-1.d0*( l + 1 )
             self%Multipoles(in1, in2, l) = norm * BsSet%Integral(fPtr,in1,in2,0,0,parvec=parvec)
          enddo
       enddo
    enddo

  end subroutine ClassBsBlocksSet

  subroutine ClassBsBlocksWriteToFile( self, FileName )
    class(ClassBsBlocks), intent(in) :: self
    character(len=*)   , intent(in) :: FileName

    integer :: uid, iostat, n, k, lm, i, j, l
    character(len=500) :: iomsg

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="unknown", &
         action  ="write"  , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in Writing "//trim(iomsg)
       stop
    endif

    write(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif

    write(uid) self%BsOrder
    write(uid) self%lmax
    write(uid) self%NBsplines

    n  = self%NBsplines
    k  = self%BsOrder
    lm = self%lmax
    !
    write(uid) ((self%Overlap    (i,j)  ,i=max(1,j-k+1),min(n,j+k-1)),j=1,n)
    write(uid) ((self%SecondDeriv(i,j)  ,i=max(1,j-k+1),min(n,j+k-1)),j=1,n)
    write(uid) ((self%Centrifugal(i,j)  ,i=max(1,j-k+1),min(n,j+k-1)),j=1,n)
    write(uid)(((self%Multipoles (i,j,l),i=max(1,j-k+1),min(n,j+k-1)),j=1,n),l=0,2*lm)
    
    close(uid)

  end subroutine ClassBsBlocksWriteToFile


  subroutine ClassBsBlocksReadFromFile( self, FileName )
    class(ClassBsBlocks), intent(inout) :: self
    character(len=*)   , intent(in)    :: FileName

    integer :: NBsplines, lmax, n, k, i, j, l
    integer :: uid, iostat
    character(len=500) :: iomsg

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="old"    , &
         action  ="read"   , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in reading "//trim(iomsg)
       stop
    endif
    read(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif

    read(uid) self%BsOrder
    read(uid) self%lmax
    read(uid) self%NBsplines
    NBsplines = self%NBsplines
    lmax      = self%lmax
    allocate(self%Overlap(    NBsplines, NBsplines))
    allocate(self%SecondDeriv(NBsplines, NBsplines))
    allocate(self%Centrifugal(NBsplines, NBsplines))
    allocate(self%Multipoles( NBsplines, NBsplines,0: 2 * lmax))
    
    n  = self%NBsplines
    k  = self%BsOrder
    !
    read(uid) ((self%Overlap    (i,j)  ,i=max(1,j-k+1),min(n,j+k-1)),j=1,n)
    read(uid) ((self%SecondDeriv(i,j)  ,i=max(1,j-k+1),min(n,j+k-1)),j=1,n)
    read(uid) ((self%Centrifugal(i,j)  ,i=max(1,j-k+1),min(n,j+k-1)),j=1,n)
    read(uid)(((self%Multipoles (i,j,l),i=max(1,j-k+1),min(n,j+k-1)),j=1,n),l=0,2*lmax)

    close(uid)
  end subroutine ClassBsBlocksReadFromFile

  real(kind(1d0)) function ClassBsBlocksGetOverlap(self, i, j) result(res)
    class(ClassBsBlocks), intent(in)  :: self
    integer             , intent(in)  :: i, j
    res=self%Overlap(i,j)
  end function  ClassBsBlocksGetOverlap

  real(kind(1d0)) function ClassBsBlocksGetSecondDeriv(self, i, j) result(res)
    class(ClassBsBlocks), intent(in)  :: self
    integer             , intent(in)  :: i, j
    res=self%SecondDeriv(i,j)
  end function  ClassBsBlocksGetSecondDeriv

   real(kind(1d0)) function ClassBsBlocksGetCentrifugal(self, i, j) result(res)
    class(ClassBsBlocks), intent(in)  :: self
    integer             , intent(in)  :: i, j
    res=self%Centrifugal(i,j)
  end function  ClassBsBlocksGetCentrifugal

   real(kind(1d0)) function ClassBsBlocksGetMultipoles(self, i, j, l) result(res)
    class(ClassBsBlocks), intent(in)  :: self
    integer             , intent(in)  :: i, j, l
    res=self%Multipoles(i, j, l)
  end function  ClassBsBlocksGetMultipoles
  
  subroutine ClassNuclearMomentsFree( self )
    class(ClassNuclearMoments)  , intent(inout) :: self
    if(allocated(self%Moments))deallocate(self%Moments)
    self%INITIALIZED=.FALSE.
    self%lmax=0
    if(allocated(self%name))deallocate(self%name)
  end subroutine ClassNuclearMomentsFree

  subroutine ClassNuclearMomentsInit( self, lmax )
    use ModuleConstants
    use ModuleAngularMomentum
    class(ClassNuclearMoments)  , intent(inout) :: self
    integer                     , intent(in)    :: lmax
    call self%Free()
    self%lmax = lmax 
    allocate(self%Name,source=trim(adjustl(NUCLEAR_MOMENT_NAME)))
    self%Name = NUCLEAR_MOMENT_NAME
    allocate(self%Moments(0:2*lmax,-2*lmax:2*lmax))
    self%Moments = 0.d0
    self%INITIALIZED=.TRUE.
  end subroutine ClassNuclearMomentsInit



    subroutine ClassNuclearMomentsSet( self, MolGeom )
    use ModuleConstants
    use ModuleAngularMomentum
    class(ClassNuclearMoments)  , intent(inout) :: self
    type(ClassMolecularGeometry), intent(in)    :: MolGeom

    integer :: iAtom, l, im, m, nm, lmax
    integer, allocatable :: mlist(:)
    real(kind(1d0)) :: Q, R, theta, phi, Rvec(3), factor
    type(ClassXlmSymmetricSet) :: XlmSymSet
    type( ClassIrrep ), pointer :: TotSymIrrep

    if(.not.self%INITIALIZED)then
       write(*,*) "ClassNuclearMomentsSet: ClassNuclearMoments is not initialized"
       stop
    endif
    lmax = self%lmax 

    !.. Determines the totally symmetric m for the various m
    TotSymIrrep => GlobalGroup%GetTotSymIrrep()
    call XlmSymSet%Init(2*lmax,TotSymIrrep)

    do iAtom = 1, MolGeom%GetNat()
       Q = - MolGeom%GetCharge(iAtom)
       call MolGeom%GetCoords(iAtom,Rvec)
       call CartesianToSpherical(Rvec,R,theta,phi)
       do l = 0, 2*lmax
          nm = XlmSymSet%GetNm( l )
          if(nm==0)cycle
          factor = 4.d0 * PI / dble( 2 * l + 1 ) * Q * R**l
          !.. Only totally symmetric components contribute
          call XlmSymSet%GetMList( l, mlist )
          do im=1,nm
             m=mlist(im)
             self%Moments(l,m) = self%Moments(l,m) + factor * EvalXlm(l,m,theta,phi)

          enddo
       enddo
    enddo
  end subroutine ClassNuclearMomentsSet



  subroutine ClassNuclearMomentsWriteToFile( self, FileName )
    use ModuleConstants
    use ModuleAngularMomentum
    class(ClassNuclearMoments) , intent(in) :: self
    character(len=*)           , intent(in) :: FileName
    
    integer                    :: l, m
    integer                    :: uid, iostat
    character(len=500)         :: iomsg, sBuf


    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="unknown", &
         action  ="write"  , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in Writing "//trim(iomsg)
       stop
    endif
    
    write(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif

    write(uid) self%lmax
    sBuf = self%Name
    write(uid) sBuf
    
    write(uid) ((self%Moments(l,m),m=-l,l),l=0,2*self%lmax)
    
    close(uid)
  end subroutine ClassNuclearMomentsWriteToFile

  subroutine ClassNuclearMomentsReadFromFile( self, FileName )
    use ModuleConstants
    use ModuleAngularMomentum
    class(ClassNuclearMoments)  , intent(inout) :: self
    character(len=*)            , intent(in)    :: FileName
    
    integer                    :: l, m
    integer                    :: uid, iostat
    character(len=500)         :: iomsg, sBuf

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="old"    , &
         action  ="read"   , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in reading "//trim(iomsg)
       stop
    endif
  
    read(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif
    
    read(uid) self%lmax
    read(uid) sBuf
    allocate(self%Name,source=trim(adjustl(sBuf)))

    allocate(  self%Moments(0:2*self%lmax,-2*self%lmax:2*self%lmax))
    self%Moments=0.d0
    read(uid)((self%Moments(l,m),m=-l,l),l=0,2*self%lmax)
    close(uid)
  end subroutine ClassNuclearMomentsReadFromFile

  real(kind(1d0)) function ClassNuclearMomentsGetMoments(self, l, m) result(res)
    class(ClassNuclearMoments), intent(in)  :: self
    integer                   , intent(in)  :: l, m
    res=self%Moments( l, m )
  end function  ClassNuclearMomentsGetMoments

  integer function Get_Moment_ID(l,m) result(id)
    integer, intent(in) :: l, m
    id = lmPairIndex(l,m)
  end function Get_Moment_ID

  !Wrappers
  integer function ClassIntegralSingleton_Get_Multipole_ID(self,l,m) result(id)
    class(ClassIntegralSingleton), intent(in) :: self
    integer, intent(in) :: l, m
    id = Get_Multipole_ID(l,m)
  end function ClassIntegralSingleton_Get_Multipole_ID

  !.. From l in input, l=0,1,2,... returns the id of 
  !   the integral of 1/r^{l+1}
  integer function Get_Multipole_ID(l,m) result(id)
    integer, intent(in) :: l, m
    id = I1B_MULTIPO + lmPairIndex(l,m) - 1
  end function Get_Multipole_ID

  subroutine ClassIntegralSingleton_PairTolm( self, ind, l, m )
    class(ClassIntegralSingleton), intent(in) :: self
    integer, intent(in) :: ind
    integer, intent(out):: l,m
    call PairTolm( ind, l, m )
  end subroutine ClassIntegralSingleton_PairTolm


  integer function ClassIntegralSingleton_GetnMOv(self,iIrr) result(ires)
    class(ClassIntegralSingleton), intent(in) :: self
    integer                      , intent(in) :: iIrr
    ires = self%Int_1B(I1B_OVERLAP)%nMOv(iIrr)
  end function ClassIntegralSingleton_GetnMOv
  
  integer function ClassIntegralSingleton_GetnHYv(self,iIrr) result(ires)
    class(ClassIntegralSingleton), intent(in) :: self
    integer                      , intent(in) :: iIrr
    ires = self%Int_1B(I1B_OVERLAP)%nHYv(iIrr)
  end function ClassIntegralSingleton_GetnHYv
  
  integer function ClassIntegralSingleton_GetnBS(self) result(ires)
    class(ClassIntegralSingleton), intent(in) :: self
    ires = self%Int_1B(I1B_OVERLAP)%nBS
  end function ClassIntegralSingleton_GetnBS
  
  subroutine PairTolm( ind, l, m )
    integer, intent(in) :: ind
    integer, intent(out):: l,m
    l = int(sqrt(dble(ind-1)))
    m = ind - l**2 -l -1
  end subroutine PairTolm

  integer function ClassIntegralSingleton_lmPairIndex(self,l,m) result( lmPair )
    class(ClassIntegralSingleton), intent(in) :: self
    integer, intent(in) :: l, m
    lmPair = lmPairIndex(l,m)
  end function ClassIntegralSingleton_lmPairIndex

  integer function lmPairIndex(l,m) result( lmPair )
    integer, intent(in) :: l, m
    lmPair = 0
    if(abs(m)>l)return
    lmPair = l**2 + l + m + 1
  end function lmPairIndex

  subroutine ClassIntegralSingleton_SetStorage(self, StorageDir ) 
    class(ClassIntegralSingleton), intent(inout) :: self
    character(len=*)             , intent(in)    :: StorageDir
    allocate(self%store,source=trim(adjustl(StorageDir)))
  end subroutine ClassIntegralSingleton_SetStorage

  subroutine ClassIntegralSingleton_Setlmax(self, lmax ) 
    class(ClassIntegralSingleton), intent(inout) :: self
    integer                      , intent(in)    :: lmax
    self%lmax = lmax
  end subroutine ClassIntegralSingleton_Setlmax

  subroutine ClassIntegralSingleton_SetnInactive(self, nInactive ) 
    class(ClassIntegralSingleton), intent(inout) :: self
    integer                      , intent(in)    :: nInactive(:)
    allocate(self%nInactive,source=nInactive)
  end subroutine ClassIntegralSingleton_SetnInactive

  subroutine ClassIntegralSingleton_Free(self) 
    use, intrinsic :: ISO_FORTRAN_ENV
    class(ClassIntegralSingleton), intent(inout) :: self
    integer :: iType, nMultipoles, iMult

    if(allocated(self%store))deallocate(self%store)
    self%lmax=0

    if(allocated(self%Int_1B))then
       do iType=1,I1B_N_TYPES-1
          call self%Int_1B(iType)%free()
       enddo
       nMultipoles = size(self%Int_1B,1) - I1B_N_TYPES + 1
       do iMult = 1, nMultipoles
          call self%Int_1B( I1B_N_TYPES + iMult - 1 )%free()
       enddo
       deallocate(self%Int_1B)
    endif

    !.. Save to disk bielectronic integrals
    call self%Int_2B%free()

    !.. Save to disk nuclear moments
    call self%NucMoments%free()

    !.. Save to disk localized-orbital moments
    call self%LocMoments%free()
    
    !.. Save to disk B-spline integral blocks
    call self%BsBlocks%free()
    
  end subroutine ClassIntegralSingleton_Free



  subroutine ClassIntegralSingleton_Init(self, &
       lmax, nMOv, nHybv, NBsplines, BsOrder,  nInactive ) 
    use, intrinsic :: ISO_FORTRAN_ENV
    class(ClassIntegralSingleton), intent(inout) :: self
    integer                      , intent(in)    :: lmax
    integer                      , intent(in)    :: nMOv(:)  !vec num MO      per irrep
    integer                      , intent(in)    :: nHybv(:) !vec num Hyb Orb per irrep
    integer                      , intent(in)    :: NBsplines     !num ext rad Orb (e.g., ext BSplines)
    integer                      , intent(in)    :: BsOrder !for Bsplines, it is 2*kmax - 1
    integer, optional            , intent(in)    :: nInactive(:)!vec num inactive MO per irrep
!    type(ClassMolecularGeometry) , intent(in)    :: MolGeom

    integer             :: nExt, BandWidth
    character(len=10)   :: sBuf
    character(len=1000) :: FileName
    integer             :: iType, nMultipoles, iMult, l, m
    character(len=:), allocatable :: name !
    type(ClassIrrep), pointer     :: OpIrrep
    type(ClassXlm)                :: Xlm

    !BsOrder   = BsSet%GetOrder()
    !NBsplines = BsSet%GetNBsplines()
    nExt      = NBsplines - 2 * BsOrder + 1 
    !*** This change is also in ModuleUKRMolInterface, to avoid the last bsplines and perform the comparison with Juan's code
    nExt      = NBsplines - 2 * BsOrder + 1 - (BsOrder - 2)
    BandWidth =             2 * BsOrder - 1
    
    self%lmax = lmax
    allocate(self%Int_1B( I1B_N_TYPES - 1 + (2*self%lmax+1)**2 ))
    !.. Init monoelectronic integrals
    !    ( self, name, OpIrrep, lmax, nMOv, nHybv, nExt, BandWidth )
    call Xlm%Init(0,0)
    OpIrrep => Xlm%GetIrrep( GlobalGroup )
    !write(*,*) OpIrrep%getName()
    !pause
    do iType=1,I1B_N_TYPES-1
       OpIrrep => Xlm%GetIrrep( GlobalGroup )
       allocate(name,source=trim(adjustl(I1B_ID_LIST( iType ))))
       !*** we want to assign the opIrrep for x,y,z and nabla components:
       ! if     ((name.eq."Coord_X").or.(name.eq."Nabla_X"))then
       !    OpIrrep => XYZIrrep%GetxIrrep()
       ! elseif ((name.eq."Coord_Y").or.(name.eq."Nabla_Y"))then
       !    OpIrrep => XYZIrrep%GetyIrrep( )
       ! elseif ((name.eq."Coord_Z").or.(name.eq."Nabla_Z"))then
       !    OpIrrep => XYZIrrep%GetzIrrep()
       ! endif
       call self%Int_1B(iType)%Init(name, OpIrrep, lmax, nMOv, nHybv, nExt, BandWidth)
       deallocate(name)
    enddo
    nMultipoles = size(self%Int_1B,1) - I1B_N_TYPES + 1
    do iMult = 1, nMultipoles
       call PairTolm( iMult, l, m )
       call Xlm%Init(l,m)
       OpIrrep => Xlm%GetIrrep( GlobalGroup )
       write(sBuf,"(i5)") iMult
       FileName=trim( I1B_ID_LIST( I1B_N_TYPES ) ) // trim( adjustl( sBuf ) )
       allocate(name,source=trim(adjustl(FileName)))
       call self%Int_1B( I1B_N_TYPES + iMult - 1 )%Init(name, OpIrrep, lmax, nMOv, nHybv, nExt, BandWidth)
       deallocate(name)
    enddo
    !.. Init bielectronic integrals
    !( self, name, OpIrrep, nMOv, nHYv )
    call Xlm%Init(0,0)
    OpIrrep => Xlm%GetIrrep( GlobalGroup )
     allocate(name,source=trim(adjustl(I2B_ID_NAME)))
    call self%Int_2B%Init( name, OpIrrep, nMOv, nHybv)
    deallocate(name)
    !.. Init nuclear moments
    call self%NucMoments%Init( lmax )
    !.. Init localized-orbital moments
    call self%LocMoments%Init(lmax , nMOv )
    !.. Init BsBlocks
    call self%BsBlocks%Init( lmax, NBsplines, BsOrder )

    !.. Build h-tilde for all the block types.
    !***

    
    
  end subroutine ClassIntegralSingleton_Init

  
  subroutine ClassIntegralSingleton_WriteToFile(self) 
    use, intrinsic :: ISO_FORTRAN_ENV
    class(ClassIntegralSingleton), intent(in) :: self
    character(len=10)   :: sBuf
    character(len=1000) :: FileName
    integer :: iType, nMultipoles, iMult
    character(len=:), allocatable :: intdir

    if(.not.allocated(self%store))then
       write(ERROR_UNIT,"(A)") "INTEGRAL STORAGE DIRECTORY NOT DEFINED - HALT FORCED"
       stop
    endif

    call execute_command_line("mkdir -p "//self%store)
    allocate(intdir,source=self%store//INTEGRAL_SUBDIR)
    call execute_command_line("mkdir -p "//intdir)

    !.. Save to disk monoelectronic integrals
    do iType=1,I1B_N_TYPES-1
       call self%Int_1B(iType)%write( intdir // trim( I1B_ID_LIST( iType ) ) )
    enddo
    nMultipoles = size(self%Int_1B,1) - I1B_N_TYPES + 1
    !write(*,*) size(self%Int_1B,1),I1B_N_TYPES
    !pause
    do iMult = 1, nMultipoles
       !write(*,*) "nMultipolesnMultipoles",iMult,nMultipoles
       write(sBuf,"(i5)") iMult
       FileName=intdir // trim( I1B_ID_LIST( I1B_N_TYPES ) ) // trim( adjustl( sBuf ) )
       call self%Int_1B( I1B_N_TYPES + iMult - 1 )%write( trim( FileName ) )
    enddo

    !.. Save to disk bielectronic integrals
    call self%Int_2B%Write( intdir // BIELECTRONIC_INTEGRAL_FILE )

    !.. Save to disk nuclear moments
    call self%NucMoments%Write( intdir // NUCLEAR_MOMENT_FILE )

    !.. Save to disk localized-orbital moments
    call self%LocMoments%Write( intdir // LOC_ORBITAL_MOMENT_FILE )
    
    !.. Save to disk B-spline integral blocks
    call self%BsBlocks%Write( intdir // BSPLINE_INTEGRAL_FILE )

    !*** ADD SAVING OF THE CONDITIONING MATRICES AND VECTORS
!!$    call SaveVector(dir//"/R_eval",reval,"formatted")
!!$    call revec%Write(dir//"R_evec",      "formatted")
    !********************************
    !********************************
    !********************************
    !********************************
    !********************************
    !********************************
    !********************************
    
  end subroutine ClassIntegralSingleton_WriteToFile

  subroutine ClassIntegralSingleton_ReadFromFile(self) 
    use, intrinsic :: ISO_FORTRAN_ENV
    class(ClassIntegralSingleton), intent(inout) :: self
    character(len=10)   :: sBuf
    character(len=1000) :: FileName
    integer :: iType, nMultipoles, iMult
    character(len=:), allocatable :: intdir

    if(.not.allocated(self%store))then
       write(ERROR_UNIT,"(A)") "INTEGRAL STORAGE DIRECTORY NOT DEFINED - HALT FORCED"
       stop
    endif

    allocate(intdir,source=self%store//INTEGRAL_SUBDIR)

    allocate(self%Int_1B( I1B_N_TYPES - 1 + (2*self%lmax+1)**2 ))

    !.. Read from disk monoelectronic integrals
    do iType=1,I1B_N_TYPES-1
       call self%Int_1B(iType)%Read( intdir // trim( I1B_ID_LIST( iType ) ) )
    enddo

    nMultipoles = size(self%Int_1B,1) - I1B_N_TYPES + 1

    do iMult = 1, nMultipoles
       !write(*,*) "nMultipolesnMultipoles",iMult,nMultipoles
       write(sBuf,"(i5)") iMult
       FileName=intdir // trim( I1B_ID_LIST( I1B_N_TYPES ) ) // trim( adjustl( sBuf ) )
       call self%Int_1B( I1B_N_TYPES + iMult - 1 )%Read( trim( FileName ) )
    enddo
    
    !.. Read from disk bielectronic integrals
    call self%Int_2B%Read( intdir // BIELECTRONIC_INTEGRAL_FILE )

    !.. Read from disk nuclear moments
    call self%NucMoments%Read( intdir // NUCLEAR_MOMENT_FILE )

    !.. Read from disk localized-orbital moments
    call self%LocMoments%Read( intdir // LOC_ORBITAL_MOMENT_FILE )
    
    !.. Read from disk B-spline integral blocks
    call self%BsBlocks%Read( intdir // BSPLINE_INTEGRAL_FILE )
    
  end subroutine ClassIntegralSingleton_ReadFromFile

  !.. Orthonormalizes the external B-spline basis to the internal hybrid basis.
  !   MUST ASSUME AN ARBITRARY VALUE FOR THE MULTIPLICITY AT THE BOUNDARY
  !..
  subroutine ClassIntegralSingleton_Condition_beS(self) 
    class(ClassIntegralSingleton), target, intent(in) :: self
    real(kind(1d0)), parameter        :: INT_DVR_NRM1_THR=1.d-12
    integer                           :: i, l, m, n, lmpair, iIrr, nIrreps, im, n_Overlap_B_ext
    integer             , allocatable :: mlist(:)
    type(ClassMatrix)                 :: S_beSp_beSp, R_beSp_beSp, S_deSp_B, S_hig_B, revec, mat
    real(kind(1d0))     , allocatable :: reval(:), vint(:), vaux(:)
    type(DVecContainer) , pointer     :: Bnorm(:)
    type(ClassXlmSymmetricSet), pointer :: XlmSymSet
    type(ClassXlm)                    :: Xlm
    type(ClassIrrep)    , pointer     :: IrrList(:), Irr
    real(kind(1d0))                   :: dS_BB, normf, norm_ext
    
    lmpair      = self%lmPairIndex(0,0)
    S_beSp_beSp = self%Get1B("S","BS_BS",lmpair,lmpair)
    S_deSp_B    = S_beSp_beSp
    dS_BB       = S_beSp_beSp%Element(1,1)
    n           = S_beSp_beSp%NRows()

    call S_deSp_B%   drop("columns",2,n)
    call S_deSp_B%   drop("rows"   ,1,1)
    call S_deSp_B%   drop("rows"   ,n,n)
    
    call S_beSp_beSp%drop("columns",1,1)
    call S_beSp_beSp%drop("rows"   ,1,1)
    call S_beSp_beSp%drop("columns",n,n)
    call S_beSp_beSp%drop("rows"   ,n,n)

    R_beSp_beSp = self%Get1B("R","BS_BS",lmpair,lmpair)

    call R_beSp_beSp%drop("columns",1,1)
    call R_beSp_beSp%drop("rows"   ,1,1)
    call R_beSp_beSp%drop("columns",n,n)
    call R_beSp_beSp%drop("rows"   ,n,n)

    call R_beSp_beSp%Diagonalize(S_beSp_beSp,reval,revec)

    !.. Overlap between boundary and FEDVR functions
    vaux=S_deSp_B%FetchColumn(1)
    call S_deSp_B%Multiply(revec,"Left","T")
    norm_ext = S_deSp_B%Norm2()

    !.. Determines how many external FEDVR functions
    !   overlap with the boundary function
    n_Overlap_B_ext = S_deSp_B%NRows()
    do i = S_deSp_B%NRows(), 1, -1
       if(abs(S_deSp_B%Element(i,1))/sqrt(dS_BB) < INT_DVR_NRM1_THR)exit
       if(n_Overlap_B_ext<1)exit
       n_Overlap_B_ext = n_Overlap_B_ext - 1
    enddo
    
    allocate(Bnorm(self%lmPairIndex(self%lmax,self%lmax)))
    
    nIrreps =  GlobalGroup%GetnIrreps()
    IrrList => GlobalGroup%GetIrrepList()
    do iIrr = 1, nIrreps
       
       Irr       =  IrrList( iIrr )
       XlmSymSet => GlobalXlmSet%GetSymSet( Irr )
       
       S_hiG_B   =  self%Get1B("S","BS_HY",iIrr,1)
       
       do l=0,self%lmax
          
          call XlmSymSet%GetMList( l, mlist )
          do im = 1, size(mlist)
             m = mlist( im )
             
             lmpair = self%lmPairIndex(l,m)
             S_hiG_B   =  self%Get1B("S","BS_HY",iIrr,lmpair)
             call S_hig_B%drop("columns",2,n)
             vint = S_hig_B%FetchColumn(1)

             allocate( Bnorm( lmpair )%v( S_hiG_B%NRows() + 1 + S_deSp_B%NRows() ) )
             normf = 1.d0 / sqrt( dS_BB - sum(vint*vint) - norm_ext**2 )
             Bnorm(lmpair)%v = 0.d0
             Bnorm(lmpair)%v(1:size(vint)   ) = - vint * normf
             Bnorm(lmpair)%v(  size(vint)+1 ) =  dS_BB * normf
             Bnorm(lmpair)%v(  size(vint)+2:) = - vaux * normf

          enddo
       enddo
       
       !**** MUST KEEP A COPY OF BNORM IN THE SINGLETON,
       !     ALONGSIDE 
       deallocate(Bnorm)
       
    enddo
    
  end subroutine ClassIntegralSingleton_Condition_beS

  function ClassIntegralSingleton_Get1B(self, OpLabel, otypes, iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet, ilm) result(Mat)
    class(ClassIntegralSingleton), target, intent(in) :: self
    character(len=*)                     , intent(in) :: OpLabel, otypes
    integer                              , intent(in) :: iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet
    integer, optional                    , intent(in) :: ilm
    integer                    :: i, iKind, imp
    type(ClassMatrix), pointer :: mat

    !iIrrOrbBra = GlobalGroup%GetIrrepIndex(IrrOrbBra)
    !iIrrOrbKet = GlobalGroup%GetIrrepIndex(IrrOrbKet)

    iKind=-1
    do i = 1, I1B_N_TYPES
       if(OpLabel(1:3) .is. I1B_ID_LIST(i)(1:3))then
          iKind=i
          exit
       endif
    enddo
    if(iKind<0)then
       call ErrorMessage("Could not find the operator "//OpLabel//" in Get1B")
       ERROR STOP 1
    endif
    
    if(iKind.eq.I1B_MULTIPO)then
       if(.not.present(ilm))then
          write(*,*) "ClassIntegralSingleton_Get1B: ilm must be present in the input for the multipole integrals"
          stop
       endif
       imp = ilm - 1
    else
       if(present(ilm))then
          write(*,*) "Inconsistent apeareance of ilm optional variable for other than the multipole operator"
          stop
       endif
       imp = 0
    endif
    
    !.. If the bra or ket orbitals are of type  ClassVirtualExternalChannel,
    !   then iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet correspond to the integer
    !   label of the pairs lbra mbra and lket mket respectively.
    if(otypes.is."MO_MO")then
       mat => self%Int_1B(ikind+imp)%MO_MO(iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet)
    elseif(otypes.is."MO_HY")then
       mat => self%Int_1B(ikind+imp)%MO_HY(iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet)
    elseif(otypes.is."HY_HY")then
       mat => self%Int_1B(ikind+imp)%HY_HY(iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet)
    elseif(otypes.is."HY_BS")then
       mat => self%Int_1B(ikind+imp)%HY_BS(iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet)
    elseif(otypes.is."BS_BS")then
       mat => self%Int_1B(ikind+imp)%BS_BS(iSpatialSymmetryOrbBra, iSpatialSymmetryOrbKet) 
    endif

  end function  ClassIntegralSingleton_Get1B

  function ClassIntegralSingleton_GetLocMoments(self, iIrrOrbBra, iIrrOrbKet, ilm) result(Mat)
    class(ClassIntegralSingleton), target, intent(in) :: self
    integer                              , intent(in) :: iIrrOrbBra, iIrrOrbKet
    integer                              , intent(in)  :: ilm
    type(ClassMatrix), pointer :: mat
    mat => self%LocMoments%MO_MO(iIrrOrbBra, iIrrOrbKet, ilm)
  end function  ClassIntegralSingleton_GetLocMoments

  logical function ClassIntegralSingleton_LocMomentsIsInit(self, iIrrOrbBra, iIrrOrbKet, ilm) result(res)
    class(ClassIntegralSingleton), intent(in)  :: self
    integer                      , intent(in)  :: iIrrOrbBra, iIrrOrbKet
    integer                      , intent(in)  :: ilm
    type(ClassMatrix), pointer :: mat
    res = self%LocMoments%MO_MO(iIrrOrbBra, iIrrOrbKet, ilm)%IsInitialized()
  end function  ClassIntegralSingleton_LocMomentsIsInit
  
  subroutine ClassIntegralSingleton_ConvertHinHtilde(self, TotNElectrons, nInactive) 
    use, intrinsic :: ISO_FORTRAN_ENV
    class(ClassIntegralSingleton), intent(inout) :: self
    real(kind(1d0))              , intent(in)    :: TotNElectrons
    integer                      , intent(in)    :: nInactive(:)
    integer                       :: iType, hType, tType, oType, nr, nc, ix
    integer                       :: iI , iIrrOrbBra, iIrrOrbKet, nIrreps
    type(ClassMatrix4D), pointer  :: J4D, K4D
    type(ClassMatrix)             :: J2mK
    character(len=:), allocatable :: intdir
    
    if(.not.allocated(self%store))then
       write(ERROR_UNIT,"(A)") "INTEGRAL STORAGE DIRECTORY NOT DEFINED - HALT FORCED"
       stop
    endif
    allocate(intdir,source=self%store//INTEGRAL_SUBDIR)
    hType = I1B_HAMILTO
    tType = I1B_KINETIC
    oType = I1B_MULTIPO

    nIrreps = size(self%Int_1B(hType)%MO_MO,1)
    do iIrrOrbKet = 1, nIrreps
       do iIrrOrbBra = 1, nIrreps
          !.. Could it be the case that MO_MO is initialized but MO_HY not?
          !if(.not.self%Int_1B(hType)%MO_MO(iIrrOrbBra, iIrrOrbKet)%IsInitialized()) cycle
          if(self%Int_1B(hType)%MO_MO(iIrrOrbBra, iIrrOrbKet)%IsInitialized())then
             nr = self%Int_1B(hType)%MO_MO(iIrrOrbBra, iIrrOrbKet)%NRows()
             nc = self%Int_1B(hType)%MO_MO(iIrrOrbBra, iIrrOrbKet)%NColumns()
             call J2mK%InitFull(nr,nc)
             J2mK = 0.d0
             do iI = 1, nIrreps
                if(.not.self%BielIsInitialized("LL_LL", iI        , iI, iIrrOrbBra, iIrrOrbKet ))cycle
                if(.not.self%BielIsInitialized("LL_LL", iIrrOrbBra, iI, iI        , iIrrOrbKet ))cycle
                J4D => self%GetBiel( "LL_LL", iI        , iI, iIrrOrbBra, iIrrOrbKet )
                K4D => self%GetBiel( "LL_LL", iIrrOrbBra, iI, iI        , iIrrOrbKet )
                do ix =1, nInactive(iI)
                   J2mK%A(:,:) = J2mK%A(:,:) + 2 * J4D%A(ix,ix,:,:) - K4D%A(:,ix,ix,:)
                end do
             end do
             call self%Int_1B(hType)%MO_MO(iIrrOrbBra, iIrrOrbKet)%Add(J2mK)
             call J2mK%Free()
          endif
          
          if(self%Int_1B(hType)%MO_HY(iIrrOrbBra, iIrrOrbKet)%IsInitialized())then
             nr = self%Int_1B(hType)%MO_HY(iIrrOrbBra, iIrrOrbKet)%NRows()
             nc = self%Int_1B(hType)%MO_HY(iIrrOrbBra, iIrrOrbKet)%NColumns()
             call J2mK%InitFull(nr,nc)
             J2mK = 0.d0
             do iI = 1, nIrreps
                if(.not.self%BielIsInitialized("LL_LS", iI        , iI, iIrrOrbBra, iIrrOrbKet ))cycle
                if(.not.self%BielIsInitialized("LL_LS", iIrrOrbBra, iI, iI        , iIrrOrbKet ))cycle
                J4D => self%GetBiel( "LL_LS", iI        , iI, iIrrOrbBra, iIrrOrbKet )
                K4D => self%GetBiel( "LL_LS", iIrrOrbBra, iI, iI        , iIrrOrbKet )
                do ix =1, nInactive(iI)
                   J2mK%A(:,:) = J2mK%A(:,:) + 2 * J4D%A(ix,ix,:,:) - K4D%A(:,ix,ix,:)
                end do
             end do
             call self%Int_1B(hType)%MO_HY(iIrrOrbBra, iIrrOrbKet)%Add(J2mK)
             call J2mK%Free()
          endif
          
          if(self%Int_1B(hType)%HY_BS(iIrrOrbBra, iIrrOrbKet)%IsInitialized())then
             nr = self%Int_1B(hType)%HY_BS(iIrrOrbBra, iIrrOrbKet)%NRows()
             nc = self%Int_1B(hType)%HY_BS(iIrrOrbBra, iIrrOrbKet)%NColumns()
             call J2mK%InitFull(nr,nc)
             J2mK = self%Int_1B(oType)%HY_BS(iIrrOrbBra, iIrrOrbKet)
             call J2mK%Multiply(TotNElectrons)
             call self%Int_1B(hType)%HY_BS(iIrrOrbBra, iIrrOrbKet)%Add(J2mK)
             call J2mK%Free()
          endif
          
          if(self%Int_1B(hType)%BS_BS(iIrrOrbBra, iIrrOrbKet)%IsInitialized())then
             nr = self%Int_1B(hType)%BS_BS(iIrrOrbBra, iIrrOrbKet)%NRows()
             nc = self%Int_1B(hType)%BS_BS(iIrrOrbBra, iIrrOrbKet)%NColumns()
             call J2mK%InitFull(nr,nc)
             J2mK = self%Int_1B(oType)%BS_BS(iIrrOrbBra, iIrrOrbKet)
             call J2mK%Multiply(TotNElectrons)
             call self%Int_1B(hType)%BS_BS(iIrrOrbBra, iIrrOrbKet)%Add(J2mK)
             call J2mK%Free()
          endif

       end do
    end do
    
    !.. Save to disk (it overrides the original h)
    !   call self%Int_1B(hType)%write( intdir // trim( I1B_ID_LIST( hType ) ) )
    
  end subroutine ClassIntegralSingleton_ConvertHinHtilde

  logical function ClassIntegralSingleton_OneBIsInitialized(self, OpLabel, otypes, iIrrOrbBra, iIrrOrbKet, ilm) result(res)
    class(ClassIntegralSingleton), intent(in)  :: self
    character(len=*)             , intent(in)  :: OpLabel, otypes
    integer                      , intent(in)  :: iIrrOrbBra, iIrrOrbKet
    integer, optional            , intent(in)  :: ilm
    integer  :: i, iKind, imp
    res = .FALSE.
    
    iKind=-1
    do i = 1, I1B_N_TYPES
       if(OpLabel(1:3) .is. I1B_ID_LIST(i)(1:3))then
          iKind=i
          exit
       endif
    enddo
    if(iKind<0)then
       call ErrorMessage("Could not find the operator "//OpLabel//" in Get1B")
       ERROR STOP 1
    endif

    if(iKind.eq.I1B_MULTIPO)then
       if(.not.present(ilm))then
          write(*,*) "ClassIntegralSingleton_OneBIsInitialized: &
               ilm must be present in the input for the multipole integrals"
          stop
       endif
       imp = ilm - 1
    else
       if(present(ilm))then
          write(*,*) "ClassIntegralSingleton_OneBIsInitialized: &
               Inconsistent apeareance of ilm optional variable for other than the multipole operator"
          stop
       endif
       imp = 0
    endif
    if(otypes.is."MO_MO")then
       res = self%Int_1B(ikind+imp)%MO_MO(iIrrOrbBra, iIrrOrbKet)%IsInitialized()
    elseif(otypes.is."MO_HY")then
       res = self%Int_1B(ikind+imp)%MO_HY(iIrrOrbBra, iIrrOrbKet)%IsInitialized()
    elseif(otypes.is."HY_HY")then
       res = self%Int_1B(ikind+imp)%HY_HY(iIrrOrbBra, iIrrOrbKet)%IsInitialized()
    elseif(otypes.is."HY_BS")then
       res = self%Int_1B(ikind+imp)%HY_BS(iIrrOrbBra, iIrrOrbKet)%IsInitialized()
    elseif(otypes.is."BS_BS")then
       res = self%Int_1B(ikind+imp)%BS_BS(iIrrOrbBra, iIrrOrbKet)%IsInitialized()
    endif
  end function  ClassIntegralSingleton_OneBIsInitialized


  function ClassIntegralSingleton_GetBiel_blk(self, label, iI, jI, kI, lI ) result(p4Dmat)
    class(ClassIntegralSingleton), target, intent(in) :: self
    character(len=5)                     , intent(in) :: label
    integer                              , intent(in) :: iI, jI, kI, lI
    type(ClassMatrix4D), pointer :: p4Dmat
    if(    label.eq."LL_LL")then
       p4Dmat => self%Int_2B%LL_LL(iI, jI, kI, lI)
    elseif(label.eq."LL_LS")then
       p4Dmat => self%Int_2B%LL_LS(iI, jI, kI, lI)
    elseif(label.eq."LL_SS")then
       p4Dmat => self%Int_2B%LL_SS(iI, jI, kI, lI)
    elseif(label.eq."LS_LS")then
       p4Dmat => self%Int_2B%LS_LS(iI, jI, kI, lI)
    endif
  end function  ClassIntegralSingleton_GetBiel_Blk

  real(kind(1d0)) function ClassIntegralSingleton_GetBiel(self, label, iI, jI, kI, lI, i, j, k, l) result(res)
    class(ClassIntegralSingleton), intent(in) :: self
    character(len=5)             , intent(in) :: label
    integer                      , intent(in) :: iI, jI, kI, lI, i, j, k, l
    if(    label.eq."LL_LL")then
       res = self%Int_2B%LL_LL(iI, jI, kI, lI)%A( i, j, k, l)
    elseif(label.eq."LL_LS")then
       res = self%Int_2B%LL_LS(iI, jI, kI, lI)%A( i, j, k, l)
    elseif(label.eq."LL_SS")then
       res = self%Int_2B%LL_SS(iI, jI, kI, lI)%A( i, j, k, l)
    elseif(label.eq."LS_LS")then
       res = self%Int_2B%LS_LS(iI, jI, kI, lI)%A( i, j, k, l)
    endif
  end function  ClassIntegralSingleton_GetBiel

  logical function ClassIntegralSingleton_BielIsInitialized(self, label, iI, jI, kI, lI) result(res)
    class(ClassIntegralSingleton), intent(in)  :: self
    character(len=5)                           :: label
    integer                      , intent(in)  :: iI, jI, kI, lI
    if(    label.eq."LL_LL")then
       res = self%Int_2B%LL_LL(iI, jI, kI, lI)%IsInitialized()
    elseif(label.eq."LL_LS")then
       res = self%Int_2B%LL_LS(iI, jI, kI, lI)%IsInitialized()
    elseif(label.eq."LL_SS")then
       res = self%Int_2B%LL_SS(iI, jI, kI, lI)%IsInitialized()
    elseif(label.eq."LS_LS")then
       res = self%Int_2B%LS_LS(iI, jI, kI, lI)%IsInitialized()
    endif
  end function  ClassIntegralSingleton_BielIsInitialized

  
  subroutine ClassIntegral1BodyFree( self )
    class(ClassIntegral1Body), intent(inout) :: self

    integer :: iIrr1, iIrr2, ilm1, ilm2

    if(.not.self%INITIALIZED) return

    if(allocated(self%name)) deallocate(self%name)
    self%OpIrrep => NULL()

    !.. Cycle over irreps Bra
    do iIrr1 = 1, GlobalGroup%GetnIrreps()
       !.. Cycle over irreps Ket
       do iIrr2 = 1, GlobalGroup%GetnIrreps()
          call self%MO_MO(iIrr1,iIrr2)%Free()
          call self%MO_HY(iIrr1,iIrr2)%Free()
          call self%HY_HY(iIrr1,iIrr2)%Free()
       enddo
       !.. Cycle over Xlm Ket
       do ilm2 = 1, lmPairIndex(self%lmax,self%lmax)
          call self%HY_BS(iIrr1,ilm2)%Free()
       enddo
    enddo
    if(allocated(self%MO_MO)) deallocate(self%MO_MO)
    if(allocated(self%MO_HY)) deallocate(self%MO_HY)
    if(allocated(self%HY_HY)) deallocate(self%HY_HY)
    if(allocated(self%HY_BS)) deallocate(self%HY_BS)

    !.. Cycle over Xlm Bra
    do ilm1 = 1, lmPairIndex(self%lmax,self%lmax)
       !.. Cycle over Xlm Ket
       do ilm2 = 1, lmPairIndex(self%lmax,self%lmax)
          call self%BS_BS(ilm1,ilm2)%Free()
       enddo
    enddo
    if(allocated(self%BS_BS)) deallocate(self%BS_BS)

    self%lmax = -1
    if(allocated(self%nMOv)) deallocate(self%nMOv)
    if(allocated(self%nHYv)) deallocate(self%nHYv)
    self%nBS = 0
    self%INITIALIZED = .FALSE.

  end subroutine ClassIntegral1BodyFree


  subroutine ClassIntegral1BodyInit( self, name, OpIrrep, lmax, nMOv, nHybv, nExt, BandWidth )
    class(ClassIntegral1Body), intent(inout) :: self
    character(len=*)         , intent(in ) :: name
    type(ClassIrrep), pointer, intent(in ) :: OpIrrep
    integer                  , intent(in ) :: lmax
    integer                  , intent(in ) :: nMOv(:)  !vec num MO      per irrep
    integer                  , intent(in ) :: nHybv(:) !vec num Hyb Orb per irrep
    integer                  , intent(in ) :: nExt     !num ext rad Orb (e.g., ext BSplines)
    integer                  , intent(in ) :: BandWidth!for Bsplines, it is 2*kmax - 1

    integer :: nIrreps, nlm, im2, im1, m1, m2, l1, l2
    integer :: iIrr1, iIrr2, ilm1, ilm2
    integer :: nAboveDiag 
    integer, allocatable :: mlist1(:), mlist2(:)
    type(ClassIrrep), pointer :: IrrPtr1, IrrPtr2
    type(ClassIrrep), pointer :: IrrList(:)
    type(ClassXlmSymmetricSet), pointer :: XlmSymSet1, XlmSymSet2
    if(lmax<0)                             call Assert("Wrong lmax      in Int1BInit")
    if(BandWidth<0.or.mod(BandWidth,2)/=1) call Assert("Wrong BandWidth in Int1BInit")
    !*** etc. defensive programming
    
    call self%free()
    nAboveDiag = (BandWidth-1)/2
    allocate(self%name,source=trim(adjustl(name)))
    self%OpIrrep => OpIrrep
    self%lmax    =  lmax
    allocate(self%nMOv,source=nMOv)
    allocate(self%nHYv,source=nHybv)
    self%nBS=nExt
    nIrreps =  GlobalGroup%GetnIrreps()
    IrrList => GlobalGroup%GetIrrepList()
    nlm = lmPairIndex(lmax,lmax)
     allocate(self%MO_MO(nIrreps,nIrreps))
    allocate(self%MO_HY(nIrreps,nIrreps))
    allocate(self%HY_HY(nIrreps,nIrreps))
    allocate(self%HY_BS(nIrreps,nlm))
    allocate(self%BS_BS(nlm,nlm))
    !.. Cycle over irreps Bra
    do iIrr1 = 1, GlobalGroup%GetnIrreps()
       !.. Determines the only relevant Ket Irrep for the present operator 
       iIrr2 = GlobalGroup%GetIrrepIndex( IrrList(iIrr1) * OpIrrep )
       call self%MO_MO(iIrr1,iIrr2)%InitFull( nMOv (iIrr1), nMOv (iIrr2) )
       call self%MO_HY(iIrr1,iIrr2)%InitFull( nMOv (iIrr1), nHybv(iIrr2) )
       call self%HY_HY(iIrr1,iIrr2)%InitFull( nHybv(iIrr1), nHybv(iIrr2) )
       IrrPtr1 => IrrList( iIrr1 )
       IrrPtr2 => IrrList(iIrr1) * OpIrrep
       XlmSymSet1 => GlobalXlmSet%GetSymSet( IrrPtr1 )
       XlmSymSet2 => GlobalXlmSet%GetSymSet( IrrPtr2 )
       do l2 = 0, lmax
          call XlmSymSet2%GetMList( l2, mlist2 )
          if(.not.allocated(mlist2))cycle
          do im2 = 1, size(mlist2)
             m2 = mlist2(im2)
             ilm2 = lmPairIndex(l2,m2)
             call self%HY_BS(iIrr1,ilm2)%InitFull(nHybv(iIrr1),nAboveDiag)
          enddo
       enddo
       do l1 = 0, lmax
          call XlmSymSet1%GetMList( l1, mlist1 )
          if(.not.allocated(mlist1))cycle
          do im1 = 1, size(mlist1)
             m1 = mlist1(im1)
             ilm1 = lmPairIndex(l1,m1)

             do l2 = 0, lmax
                call XlmSymSet2%GetMList( l2, mlist2 )
                if(.not.allocated(mlist2))cycle
                do im2 = 1, size(mlist2)
                   m2 = mlist2(im2)
                   ilm2 = lmPairIndex(l2,m2)

                   !call self%BS_BS(ilm1,ilm2)%InitDBanded(nExt,nExt,nAboveDiag,nAboveDiag)
                   call self%BS_BS(ilm1,ilm2)%InitFull(nExt,nExt)
                enddo
             enddo

          enddo
       enddo

    enddo

    self%INITIALIZED = .TRUE.
    

  end subroutine ClassIntegral1BodyInit



  subroutine ClassIntegral1Body_Set_MO_MO_Element( self, i, j, iIrrBra, iIrrKet, val )
    class(ClassIntegral1Body), intent(inout) :: self
    integer                  , intent(in)    :: i, j, iIrrBra, iIrrKet
    real(kind(1d0))          , intent(in)    :: val
    call self%MO_MO(iIrrBra,iIrrKet)%Set(i,j,val)
  end subroutine ClassIntegral1Body_Set_MO_MO_Element

  subroutine ClassIntegral1Body_Set_MO_MO_Matrix( self, iIrrBra, iIrrKet, dmat )
    class(ClassIntegral1Body), intent(inout) :: self
    integer                  , intent(in)    :: iIrrBra, iIrrKet
    real(kind(1d0))          , intent(in)    :: dmat(:,:)
    call self%MO_MO(iIrrBra,iIrrKet)%Set(dmat)
  end subroutine ClassIntegral1Body_Set_MO_MO_Matrix

  subroutine ClassIntegral1Body_Set_MO_MO_subMatrix( self, iIrrBra, iIrrKet, dmat, nr, nc )
    class(ClassIntegral1Body), intent(inout) :: self
    integer                  , intent(in)    :: iIrrBra, iIrrKet
    real(kind(1d0))          , intent(in)    :: dmat(:,:)
    integer                  , intent(in)    :: nr, nc
    call self%MO_MO(iIrrBra,iIrrKet)%Set(dmat,nr,nc)
  end subroutine ClassIntegral1Body_Set_MO_MO_subMatrix

  real(kind(1d0)) function ClassIntegral1Body_GetMOMO( self, iIrrBra, iIrrKet, i, j) result(res)
     class(ClassIntegral1Body), intent(in) :: self
     integer                  , intent(in)    :: iIrrBra, iIrrKet, i, j
     res = self%MO_MO(iIrrBra,iIrrKet)%Element( i, j)
   end function ClassIntegral1Body_GetMOMO

   real(kind(1d0)) function ClassIntegral1Body_GetMOHY( self, iIrrBra, iIrrKet, i, j) result(res)
     class(ClassIntegral1Body), intent(in) :: self
     integer                  , intent(in)    :: iIrrBra, iIrrKet, i, j
     res = self%MO_HY(iIrrBra,iIrrKet)%Element( i, j)
   end function ClassIntegral1Body_GetMOHY

   real(kind(1d0)) function ClassIntegral1Body_GetHYHY( self, iIrrBra, iIrrKet, i, j) result(res)
     class(ClassIntegral1Body), intent(in) :: self
     integer                  , intent(in)    :: iIrrBra, iIrrKet, i, j
     res = self%HY_HY(iIrrBra,iIrrKet)%Element( i, j)
   end function ClassIntegral1Body_GetHYHY

   real(kind(1d0)) function ClassIntegral1Body_GetHYBS( self, iIrrBra, iIrrKet, i, j) result(res)
     class(ClassIntegral1Body), intent(in) :: self
     integer                  , intent(in)    :: iIrrBra, iIrrKet, i, j
     res = self%HY_BS(iIrrBra,iIrrKet)%Element( i, j)
   end function ClassIntegral1Body_GetHYBS

   real(kind(1d0)) function ClassIntegral1Body_GetBSBS( self, iIrrBra, iIrrKet, i, j) result(res)
     class(ClassIntegral1Body), intent(in) :: self
     integer                  , intent(in)    :: iIrrBra, iIrrKet, i, j
     res = self%BS_BS(iIrrBra,iIrrKet)%Element( i, j)
   end function ClassIntegral1Body_GetBSBS
   

  
  subroutine ClassIntegral1Body_Set_SubMatrix( self, BlockName, iIrrBra, iIrrKet, dmat, nr, nc, nAboveDiag )
    class(ClassIntegral1Body), intent(inout) :: self
    character(len=*)         , intent(in)    :: BlockName
    integer                  , intent(in)    :: iIrrBra, iIrrKet
    real(kind(1d0))          , intent(in)    :: dmat(:,:)
    integer                  , intent(in)    :: nr, nc
    integer, optional        , intent(in)    :: nAboveDiag
    integer                                  :: i, j
    
    select case( BlockName )
    case("MO_MO")
       !call self%MO_MO(iIrrBra,iIrrKet)%Set(dmat,nr,nc)
       self%MO_MO(iIrrBra,iIrrKet) = dmat(1:nr,1:nc)
    case("MO_HY")
       !call self%MO_HY(iIrrBra,iIrrKet)%Set(dmat,nr,nc)
       self%MO_HY(iIrrBra,iIrrKet) = dmat(1:nr,1:nc)
    case("HY_HY")
       !call self%HY_HY(iIrrBra,iIrrKet)%Set(dmat,nr,nc)
       self%HY_HY(iIrrBra,iIrrKet) = dmat(1:nr,1:nc)
    case("HY_BS")
       !call self%HY_BS(iIrrBra,iIrrKet)%Set(dmat,nr,nc)
       self%HY_BS(iIrrBra,iIrrKet) = dmat(1:nr,1:nc)
    case("BS_BS")
       !*** Build the routine to assign the matrix to the DBanded, with dmat with the appropiate structure
       !call self%BS_BS(iIrrBra,iIrrKet)%Set(dmat,nr,nc)
        ! write(*,*) "ClassIntegral1Body_Set_SubMatrix:", BlockName, iIrrBra, iIrrKet
        ! write(*,*) self%BS_BS(iIrrBra,iIrrKet)%NRows(),self%BS_BS(iIrrBra,iIrrKet)%NColumns()
        ! write(*,*) nr,nc,self%BS_BS(iIrrBra,iIrrKet)%IsDBanded(),self%BS_BS(iIrrBra,iIrrKet)%IsFull()
       !write(*,*) dmat(1,1:nc)
      ! write(*,*) nr,nc
       ! do i = 1, nr
       !    !write(*,*) i,nAboveDiag,max(1,i-nAboveDiag),min(nc,i+nAboveDiag)
       !     do j = max(1,i-nAboveDiag),min(nc,i+nAboveDiag)
       !       ! write(*,*) i,j,dmat(i,j)
       !        call self%BS_BS(iIrrBra,iIrrKet)%SetElement(i,j,dmat(i,j))
       !     enddo
       !  enddo
               !pause
        self%BS_BS(iIrrBra,iIrrKet) = dmat(1:nr,1:nc)
        !write(*,*) self%BS_BS(iIrrBra,iIrrKet)%IsDBanded(),self%BS_BS(iIrrBra,iIrrKet)%IsFull()
       !!--------------------
    end select
  end subroutine ClassIntegral1Body_Set_SubMatrix

  subroutine ClassIntegral1Body_WriteToFile( self, FileName )
    class(ClassIntegral1Body), intent(in) :: self
    character(len=*)         , intent(in) :: FileName
    integer :: uid, nIrreps, ilm1, ilm2, nlms, iIrr, iIrr1, iIrr2, iostat, i, j
    character(len=500) :: iomsg, sBuf

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="unknown", &
         action  ="write"  , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in Writing "//trim(iomsg)
       stop
    endif
    nIrreps = GlobalGroup%GetnIrreps()

    write(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif

    sBuf = self%Name
    write(uid) sBuf
    call self%OpIrrep%Write(uid)
    write(uid) self%lmax
    write(uid) (self%nMOv(iIrr),iIrr=1,nIrreps)
    write(uid) (self%nHYv(iIrr),iIrr=1,nIrreps)
    write(uid) self%nBS
    
    do iIrr1 = 1, nIrreps
       do iIrr2 = 1, nIrreps
          if(self%MO_MO(iIrr1,iIrr2)%IsInitialized())then
             write(uid) iIrr1, iIrr2
             call self%MO_MO(iIrr1,iIrr2)%Write(uid)
          endif
       enddo
    enddo

    write(uid) -1, -1

    do iIrr1 = 1, nIrreps
       do iIrr2 = 1, nIrreps
          if(self%MO_HY(iIrr1,iIrr2)%IsInitialized())then
             write(uid) iIrr1, iIrr2
             call self%MO_HY(iIrr1,iIrr2)%Write(uid)
          endif
       enddo
    enddo
    write(uid) -1, -1

    do iIrr1 = 1, nIrreps
       do iIrr2 = 1, nIrreps
          if(self%HY_HY(iIrr1,iIrr2)%IsInitialized())then
             write(uid) iIrr1, iIrr2
             call self%HY_HY(iIrr1,iIrr2)%Write(uid)
          endif
       enddo
    enddo
    write(uid) -1, -1

    nlms = size(self%HY_BS,2)
    write(uid) nlms
    do iIrr1 = 1, nIrreps
       do ilm2 = 1, nlms
          if(self%HY_BS(iIrr1,ilm2)%IsInitialized())then
             write(uid) iIrr1, ilm2
             call self%HY_BS(iIrr1,ilm2)%Write(uid)
          endif
       enddo
    enddo
    write(uid) -1, -1

    do ilm1 = 1, nlms
       do ilm2 = 1, nlms
          if(self%BS_BS(ilm1,ilm2)%IsInitialized())then
             write(uid) ilm1, ilm2
             call self%BS_BS(ilm1,ilm2)%Write(uid)
             ! write(*,*) "**********************",FileName
             ! write(*,*) ilm1,ilm2,self%BS_BS(ilm1,ilm2)%IsDBanded(),self%BS_BS(ilm1,ilm2)%IsFull()
             ! Do i=1,self%BS_BS(ilm1,ilm2)%NRows()
             !    Do j=1,self%BS_BS(ilm1,ilm2)%NColumns()
             !       write(*,*) i,j,self%BS_BS(ilm1,ilm2)%Element(i,j)
             !    enddo
             ! enddo
             ! write(*,*) self%BS_BS(ilm1,ilm2)%A
             ! pause
          endif
       enddo
    enddo
    write(uid) -1, -1

    close(uid)

    !pause

  end subroutine ClassIntegral1Body_WriteToFile


  subroutine ClassIntegral1Body_ReadFromFile( self, FileName )
    class(ClassIntegral1Body), intent(inout) :: self
    character(len=*)         , intent(in)    :: FileName
    integer :: uid, nIrreps, ilm1, ilm2, nlms, iIrr, iIrr1, iIrr2, iostat, i, j
    character(len=500) :: iomsg, sBuf

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="old"    , &
         action  ="read"   , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in reading "//trim(iomsg)
       stop
    endif
!write(*,*) Filename
    read(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       write(*,*) "Integral1Body is not Initializad"
       close(uid)
       return
    endif
    
    read(uid) sBuf
    allocate(self%Name,source=trim(adjustl(sBuf)))
    allocate(self%OpIrrep)
    call self%OpIrrep%read(uid)
    read(uid) self%lmax

    nIrreps = GlobalGroup%GetnIrreps()
    allocate(self%nMOv(nIrreps),self%nHYv(nIrreps))
    read(uid) (self%nMOv(iIrr),iIrr=1,nIrreps)
    read(uid) (self%nHYv(iIrr),iIrr=1,nIrreps)
    read(uid) self%nBS

    allocate(self%MO_MO(nirreps,nirreps))
    do
       read(uid) iIrr1, iIrr2
       if(iIrr1<0)exit
       call self%MO_MO(iIrr1,iIrr2)%read(uid)
    enddo

    allocate(self%MO_HY(nirreps,nirreps))
    do
       read(uid) iIrr1, iIrr2
       if(iIrr1<0)exit
       call self%MO_HY(iIrr1,iIrr2)%read(uid)
    enddo

    allocate(self%HY_HY(nirreps,nirreps))
    do 
       read(uid) iIrr1, iIrr2
       if(iIrr1<0)exit
       call self%HY_HY(iIrr1,iIrr2)%Read(uid)
    enddo

    read(uid) nlms

    allocate(self%HY_BS(nirreps,nlms))
    do
       read(uid) iIrr1, ilm2
       if(iIrr1<0)exit
       call self%HY_BS(iIrr1,ilm2)%Read(uid)
    enddo

    allocate(self%BS_BS(nlms,nlms))
    do
       read(uid) ilm1, ilm2
       if(ilm1<0)exit
       call self%BS_BS(ilm1,ilm2)%Read(uid)
       ! write(*,*) "*********************",FileName
       ! write(*,*) ilm1,ilm2
       ! write(*,*) self%BS_BS(ilm1,ilm2)%Element(4,1)
       ! Do i=1,self%BS_BS(ilm1,ilm2)%NRows()
       !    Do j=1,self%BS_BS(ilm1,ilm2)%NColumns()
       !       write(*,*) i,j,self%BS_BS(ilm1,ilm2)%Element(i,j)
       !    enddo
       ! enddo
       ! write(*,*) self%BS_BS(ilm1,ilm2)%A
       ! pause
    enddo

    close(uid)
!pause
  end subroutine ClassIntegral1Body_ReadFromFile


  subroutine ClassIntegral2BodyInit( self, name, OpIrrep, nMOv, nHYv )
    class(ClassIntegral2Body), intent(inout) :: self
    character(len=*)         , intent(in ) :: name
    type(ClassIrrep), pointer, intent(in ) :: OpIrrep
    integer                  , intent(in ) :: nMOv(:)  !vec num MO      per irrep
    integer                  , intent(in ) :: nHYv(:) !vec num Hyb Orb per irrep

    integer :: nIrreps
    integer :: iIrr1, iIrr2, iIrr3, iIrr4
    type(ClassIrrep), pointer :: IrrList(:)
    type(ClassIrrep), pointer :: IrrPtr


    call self%free()
    allocate(self%nMOv,source=nMOv)
    allocate(self%nHYv,source=nHYv)

    self%OpIrrep => OpIrrep
    nIrreps =  GlobalGroup%GetnIrreps()
    IrrList => GlobalGroup%GetIrrepList()
    allocate(self%name,source=trim(adjustl(name)))
    allocate(self%LL_LL(nIrreps,nIrreps,nIrreps,nIrreps))
    allocate(self%LL_LS(nIrreps,nIrreps,nIrreps,nIrreps))
    allocate(self%LL_SS(nIrreps,nIrreps,nIrreps,nIrreps))
    allocate(self%LS_LS(nIrreps,nIrreps,nIrreps,nIrreps))

    !   For orbitals all in the same group (MOs, Hybrid, or Bsplines)
    !   the irreps can be ordered as 
    do iIrr1 = 1, nIrreps
       !
       do iIrr2 = 1, nIrreps
          !
          do iIrr3 = 1, nIrreps
             !
             !.. Chose iIrr4 compatible with the symmetry.
             IrrPtr =>  IrrList(iIrr1) * IrrList(iIrr2) * IrrList(iIrr3)
             iIrr4 = GlobalGroup%GetIrrepIndex( IrrPtr )

             if( CheckIrrepOrder("LL_LL", iIrr1, iIrr2, iIrr3, iIrr4 ) ) &
                  call self%LL_LL(        iIrr1, iIrr2, iIrr3, iIrr4 )%Init( &
                  nMOv(iIrr1), nMOv(iIrr2), nMOv(iIrr3), nMOv(iIrr4) )

             if( CheckIrrepOrder("LL_LS", iIrr1, iIrr2, iIrr3, iIrr4 ) ) &
                  call self%LL_LS(        iIrr1, iIrr2, iIrr3, iIrr4 )%Init( &
                  nMOv(iIrr1), nMOv(iIrr2), nMOv(iIrr3), nHYv(iIrr4) )

             if( CheckIrrepOrder("LL_SS", iIrr1, iIrr2, iIrr3, iIrr4 ) ) &
                  call self%LL_SS(        iIrr1, iIrr2, iIrr3, iIrr4 )%Init( &
                  nMOv(iIrr1), nMOv(iIrr2), nHYv(iIrr3), nHYv(iIrr4) )

             if( CheckIrrepOrder("LS_LS", iIrr1, iIrr2, iIrr3, iIrr4 ) ) &
                  call self%LS_LS(        iIrr1, iIrr2, iIrr3, iIrr4 )%Init( &
                  nMOv(iIrr1), nHYv(iIrr2), nMOv(iIrr3), nHYv(iIrr4) )
          enddo
       enddo
    enddo

    self%INITIALIZED=.TRUE.
  end subroutine ClassIntegral2BodyInit

  subroutine ClassIntegral2BodyFree( self )
    class(ClassIntegral2Body), intent(inout) :: self

    integer :: iIrr1, iIrr2, iIrr3, iIrr4

    if(.not.self%INITIALIZED) return

    if(allocated(self%name))deallocate(self%name)

    if(associated(self%LL_LL))then
       !.. Cycle over irreps Bra1
       do iIrr1 = 1, GlobalGroup%GetnIrreps()
          !.. Cycle over irreps Bra2
          do iIrr2 = 1, GlobalGroup%GetnIrreps()
             !.. Cycle over irreps Ket1
             do iIrr3 = 1, GlobalGroup%GetnIrreps()
                !.. Cycle over irreps Ket2
                do iIrr4 = 1, GlobalGroup%GetnIrreps()
                   call self%LL_LL(iIrr1,iIrr2,iIrr3,iIrr4)%Free()
                   call self%LL_LS(iIrr1,iIrr2,iIrr3,iIrr4)%Free()
                   call self%LL_SS(iIrr1,iIrr2,iIrr3,iIrr4)%Free()
                   call self%LS_LS(iIrr1,iIrr2,iIrr3,iIrr4)%Free()
                enddo
             enddo
          enddo
       enddo

       deallocate(self%LL_LL)
       deallocate(self%LL_LS)
       deallocate(self%LL_SS)
       deallocate(self%LS_LS)
    endif

    self%lmax = -1
    if(allocated(self%nMOv)) deallocate(self%nMOv)
    if(allocated(self%nHYv)) deallocate(self%nHYv)

    self%INITIALIZED = .FALSE.

  end subroutine ClassIntegral2BodyFree



  subroutine ClassIntegral2Body_WriteToFile( self, FileName )
    class(ClassIntegral2Body), intent(in) :: self
    character(len=*)         , intent(in) :: FileName
    integer :: uid, iostat
    integer :: nIrreps, iIrr, iIrr1, iIrr2, iIrr3, iIrr4
    character(len=500) :: iomsg, sBuf
    
    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="unknown", &
         action  ="write"  , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in Writing "//trim(iomsg)
       stop
    endif
    nIrreps = GlobalGroup%GetnIrreps()
    
    write(uid) self%INITIALIZED
    
    if(.not.self%INITIALIZED)then
       close(uid)
       return
    endif

    sBuf = self%Name
    write(uid) sBuf
    call self%OpIrrep%Write(uid)
    write(uid) self%lmax
    write(uid) (self%nMOv(iIrr),iIrr=1,nIrreps)
    write(uid) (self%nHYv(iIrr),iIrr=1,nIrreps)

    do iIrr1 = 1, nIrreps
       do iIrr2 = 1, nIrreps
          do iIrr3 = 1, nIrreps
             do iIrr4 = 1, nIrreps
                if(self%LL_LL(iIrr1,iIrr2,iIrr3,iIrr4)%IsInitialized())then
                   write(uid) iIrr1, iIrr2, iIrr3, iIrr4
                   call self%LL_LL(iIrr1,iIrr2,iIrr3,iIrr4)%Write(uid)
                endif
             enddo
          enddo
       enddo
    enddo
    write(uid) -1, -1, -1, -1
    
    do iIrr1 = 1, nIrreps
       do iIrr2 = 1, nIrreps
          do iIrr3 = 1, nIrreps
             do iIrr4 = 1, nIrreps
                if(self%LL_LS(iIrr1,iIrr2,iIrr3,iIrr4)%IsInitialized())then
                   write(uid) iIrr1, iIrr2, iIrr3, iIrr4
                   call self%LL_LS(iIrr1,iIrr2,iIrr3,iIrr4)%Write(uid)
                endif
             enddo
          enddo
       enddo
    enddo
     write(uid) -1, -1, -1, -1
     
     do iIrr1 = 1, nIrreps
        do iIrr2 = 1, nIrreps
           do iIrr3 = 1, nIrreps
              do iIrr4 = 1, nIrreps
                 if(self%LL_SS(iIrr1,iIrr2,iIrr3,iIrr4)%IsInitialized())then
                    write(uid) iIrr1, iIrr2, iIrr3, iIrr4
                    call self%LL_SS(iIrr1,iIrr2,iIrr3,iIrr4)%Write(uid)
                 endif
              enddo
           enddo
        enddo
     enddo
     write(uid) -1, -1, -1, -1
     
     do iIrr1 = 1, nIrreps
        do iIrr2 = 1, nIrreps
           do iIrr3 = 1, nIrreps
              do iIrr4 = 1, nIrreps
                 if(self%LS_LS(iIrr1,iIrr2,iIrr3,iIrr4)%IsInitialized())then
                    write(uid) iIrr1, iIrr2, iIrr3, iIrr4
                    call self%LS_LS(iIrr1,iIrr2,iIrr3,iIrr4)%Write(uid)
                 endif
              enddo
           enddo
        enddo
     enddo
     write(uid) -1, -1, -1, -1

    close(uid)

  end subroutine ClassIntegral2Body_WriteToFile


  subroutine ClassIntegral2Body_ReadFromFile( self, FileName )
    class(ClassIntegral2Body), intent(inout) :: self
    character(len=*)         , intent(in)    :: FileName
    integer :: uid, iostat
    integer :: nIrreps, iIrr, iIrr1, iIrr2, iIrr3, iIrr4
    character(len=500) :: iomsg, sBuf
    logical :: ShouldBeInitialized

    open(newunit = uid     , &
         file    = FileName, &
         form    ="unformatted", &
         status  ="old"    , &
         action  ="read"   , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       write(ERROR_UNIT,"(a)") "Unable to open file "//FileName//" in reading "//trim(iomsg)
       stop
    endif

    read(uid) ShouldBeInitialized

    if(.not.ShouldBeInitialized)then
       close(uid)
       return
    endif
    
    read(uid) sBuf
    allocate(self%Name,source=trim(adjustl(sBuf)))
    allocate(self%OpIrrep)
    call self%OpIrrep%read(uid)
    read(uid) self%lmax

    nIrreps = GlobalGroup%GetnIrreps()
    allocate(self%nMOv(nIrreps),self%nHYv(nIrreps))
    read(uid) (self%nMOv(iIrr),iIrr=1,nIrreps)
    read(uid) (self%nHYv(iIrr),iIrr=1,nIrreps)

    allocate(self%LL_LL(nirreps,nirreps,nirreps,nirreps))
    do
       read(uid) iIrr1, iIrr2, iIrr3, iIrr4
       if(iIrr1<0)exit
       call self%LL_LL(iIrr1,iIrr2,iIrr3,iIrr4)%read(uid)
    enddo

    allocate(self%LL_LS(nirreps,nirreps,nirreps,nirreps))
    do
       read(uid) iIrr1, iIrr2, iIrr3, iIrr4
       if(iIrr1<0)exit
       call self%LL_LS(iIrr1,iIrr2,iIrr3,iIrr4)%read(uid)
    enddo

    allocate(self%LL_SS(nirreps,nirreps,nirreps,nirreps))
    do
       read(uid) iIrr1, iIrr2, iIrr3, iIrr4
       if(iIrr1<0)exit
       call self%LL_SS(iIrr1,iIrr2,iIrr3,iIrr4)%read(uid)
    enddo

    allocate(self%LS_LS(nirreps,nirreps,nirreps,nirreps))
    do
       read(uid) iIrr1, iIrr2, iIrr3, iIrr4
       if(iIrr1<0)exit
       call self%LS_LS(iIrr1,iIrr2,iIrr3,iIrr4)%read(uid)
    enddo

    close(uid)
    self%INITIALIZED = .TRUE.

  end subroutine ClassIntegral2Body_ReadFromFile

  

  logical function CheckIrrepOrder(BlockLabel, i1,i2,i3,i4) result(lRes)
    character(len=*), intent(in) :: BlockLabel
    integer         , intent(in) :: i1, i2, i3, i4
    select case( BlockLabel )
    case( "LL_LL" )
       lRes = CheckIrrepOrder_LL_LL(i1,i2,i3,i4)
    case( "LL_LS" )
       lRes = CheckIrrepOrder_LL_LS(i1,i2,i3,i4)
    case( "LL_SS" )
       lRes = CheckIrrepOrder_LL_SS(i1,i2,i3,i4)
    case( "LS_LS" )
       lRes = CheckIrrepOrder_LS_LS(i1,i2,i3,i4)
    end select
  end function CheckIrrepOrder

  logical function CheckIrrepOrder_LL_LL(i1,i2,i3,i4) result(lRes)
    integer, intent(in) :: i1, i2, i3, i4
    lRes = .TRUE.
    lRes = lRes .and. ( max(i1,i2) >= max(i3,i4) )
    lRes = lRes .and. ( i1 >= i2 )
    lRes = lRes .and. ( i3 >= i4 )
    lRes = lRes .and. ( i1 > i3 .or. i2 >= i4 )
  end function CheckIrrepOrder_LL_LL

  logical function CheckIrrepOrder_LL_LS(i1,i2,i3,i4) result(lRes)
    integer, intent(in) :: i1, i2, i3, i4
    lRes = ( i1 >= i2 )
  end function CheckIrrepOrder_LL_LS

  logical function CheckIrrepOrder_LL_SS(i1,i2,i3,i4) result(lRes)
    integer, intent(in) :: i1, i2, i3, i4
    lRes = .TRUE.
    lRes = lRes .and. ( i1 >= i2 )
    lRes = lRes .and. ( i3 >= i4 )
  end function CheckIrrepOrder_LL_SS

  logical function CheckIrrepOrder_LS_LS(i1,i2,i3,i4) result(lRes)
    integer, intent(in) :: i1, i2, i3, i4
    lRes = .TRUE.
    lRes = lRes .and. ( i1 >= i3 )
    lRes = lRes .and. ( i1 > i3 .or. i2 >= i4 )
  end function CheckIrrepOrder_LS_LS


  logical function ClassIntegral2BodyCheckIrrepOrder(self, BlockLabel, i1,i2,i3,i4) result(lRes)
    class(ClassIntegral2Body), intent(in) :: self
    character(len=*)         , intent(in) :: BlockLabel
    integer                  , intent(in) :: i1, i2, i3, i4
    lRes = CheckIrrepOrder(BlockLabel, i1,i2,i3,i4)
  end function ClassIntegral2BodyCheckIrrepOrder

  real(kind(1d0)) function ClassIntegral2Body_GetLLLL(self,iI1,iI2,iI3,iI4,i1,i2,i3,i4) result(res)
    class(ClassIntegral2Body), intent(in) :: self
    integer                  , intent(in) :: iI1,iI2,iI3,iI4,i1,i2,i3,i4
    res = self%LL_LL(iI1,iI2,iI3,iI4)%A(i1,i2,i3,i4)
  end function ClassIntegral2Body_GetLLLL

  real(kind(1d0)) function ClassIntegral2Body_GetLLLS(self,iI1,iI2,iI3,iI4,i1,i2,i3,i4) result(res)
    class(ClassIntegral2Body), intent(in) :: self
    integer                  , intent(in) :: iI1,iI2,iI3,iI4,i1,i2,i3,i4
    res = self%LL_LS(iI1,iI2,iI3,iI4)%A(i1,i2,i3,i4)
  end function ClassIntegral2Body_GetLLLS

  real(kind(1d0)) function ClassIntegral2Body_GetLLSS(self,iI1,iI2,iI3,iI4,i1,i2,i3,i4) result(res)
    class(ClassIntegral2Body), intent(in) :: self
    integer                  , intent(in) :: iI1,iI2,iI3,iI4,i1,i2,i3,i4
    res = self%LL_SS(iI1,iI2,iI3,iI4)%A(i1,i2,i3,i4)
  end function ClassIntegral2Body_GetLLSS

  real(kind(1d0)) function ClassIntegral2Body_GetLSLS(self,iI1,iI2,iI3,iI4,i1,i2,i3,i4) result(res)
    class(ClassIntegral2Body), intent(in) :: self
    integer                  , intent(in) :: iI1,iI2,iI3,iI4,i1,i2,i3,i4
    res = self%LS_LS(iI1,iI2,iI3,iI4)%A(i1,i2,i3,i4)
  end function ClassIntegral2Body_GetLSLS

  subroutine ClassIntegral2Body_Set4DMat( self, BlockLabel, iIrr1, iIrr2, iIrr3, iIrr4, d4mat )
    class(ClassIntegral2Body), intent(inout) :: self
    character(len=*)         , intent(in)    :: BlockLabel
    integer                  , intent(in)    :: iIrr1, iIrr2, iIrr3, iIrr4
    real(kind(1d0))          , intent(in)    :: d4mat(:,:,:,:)
    select case( BlockLabel )
    case("LL_LL")
       call self%LL_LL(iIrr1, iIrr2, iIrr3, iIrr4)%Set( d4mat )
    case("LL_LS")
       call self%LL_LS(iIrr1, iIrr2, iIrr3, iIrr4)%Set( d4mat )
    case("LS_LS")
       call self%LS_LS(iIrr1, iIrr2, iIrr3, iIrr4)%Set( d4mat )
    case("LL_SS")
       call self%LL_SS(iIrr1, iIrr2, iIrr3, iIrr4)%Set( d4mat )
    end select
  end subroutine ClassIntegral2Body_Set4DMat


end module ModuleIntegrals


