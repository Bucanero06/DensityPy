Module ModuleBasisJUAN
  use precisn
  use ukrmol_interface
  use mpi_mod 
  implicit none

  !.. expansion coefficients of the ukrmol elements in terms of the primitive basis
  real(kind(1d0)), allocatable     :: cf(:,:),cfi(:,:)      
  integer              :: SymmetryIdentifier
  integer              :: nIrreps
  integer, allocatable :: irrO(:) ! This vector is for changing the order of the irresps, in order to adapt the code to the ASTRA or DALTON ordering.
  integer, allocatable :: jrrO(:) ! The inverse of irrO
  integer, allocatable :: nLOC(:) ! Number of PI     or target    orbitals per irrep
  integer, allocatable :: nSRC(:) ! Number of hybrid or continuum orbitals per irrep
  integer, allocatable :: nPWC(:) ! Number of partial wave channels        per irrep
  integer, allocatable :: nTOT(:) ! Total Number of               orbitals per irrep
  integer, allocatable :: sTOT(:) ! Total Number of orbitals up to a given     irrep
  integer, allocatable :: iLOC(:) ! Counter 
  integer, allocatable :: iSRC(:) ! Counter

  integer                        :: n_bs,l_bs
  real(kind(1d0))                :: a,bspline_grid_start
  integer                        :: no_bsplines,bspline_order
  integer                        :: min_bspline_l,max_bspline_l,no_nodes,bspline_indices(1:2,0:20)
  DoublePrecision, allocatable   :: bsgrid(:)
  real(kind(1d0))                :: da  !ukrmol bsplines grid parameter
  real(kind(1d0)), allocatable   :: NB(:),NBin(:)  !NB is the factor to normalize Bsplines square integral to 1
  real(kind(1d0))                :: Rmin,Rmax
  integer                        :: Bs1,Bs2,dBra,dKet 
  integer                        :: NNodes,BsplineOrder,NBspl,Lmax
  integer                        :: NlastBspl,NfirstBspl,NExtBspl

  !.. this indexes will have the dimension of the matrix (minus 1, ground state), identifying the stex element
  integer, allocatable             :: iha(:),iep(:)
  !.. ukrmol orthonormalized basis and primitive sizes respectively
  integer                          :: Nbasis,Ndim,NprimUkrmol,Nukrmol,NExtOrb
  !.. location of the l and m of the multipolar expansion with the iflag of ukrmol
  integer, allocatable             :: l_mp(:),m_mp(:),flag_ml(:,:)
  !.. dimension of the above vectors
  integer                          :: Nmultipoles,Lmax_mp     

  !.. Extended orthonormalized bsplines coefficients and overlap marices
  real(kind(1d0)), allocatable   :: Xi(:,:),Xim1(:,:),Spqm1(:,:),Spq(:,:),Spq_aux(:,:),Spq2(:,:),Ono(:,:),KEno(:,:),H0no(:,:)
  real(kind(1d0)), allocatable   :: gijkl(:,:,:,:),fijkl(:,:,:,:)

  integer, allocatable :: MO(:),MeO(:),hA(:),xeP(:)
  integer, allocatable :: hB(:),xeQ(:),shA(:),sxeP(:),shB(:),sxeQ(:)
  integer, allocatable :: nirr_lm(:),irr_l(:,:),irr_m(:,:),irr_nr(:,:),irr_vs_lm(:,:)
  integer*8, allocatable    :: irr_series(:),irr_counter(:)
  integer, allocatable :: RefConf(:)
  !.. ukrmol scatci_integrals inp   definitions
  integer, allocatable :: nMO(:),iMO(:),irr_total(:),MOrel(:)
  integer              :: imin_irr(8),imax_irr(8)

  !.. For the basis indexing without spin (argument of GET_INTEGRALS) it returns the relative index within the symmetry and irr respectively
  integer, allocatable :: absTOrel(:),absTOirr(:),absTOl(:),absTOm(:),absTOnr(:)
  integer, allocatable :: absTOukr(:)
  character*3, allocatable :: absTOchar(:)

  
  !.. For the basis orthonormalization and hamiltonian
  real(kind(1d0)), allocatable   :: Hirr(:,:),Over(:,:),En(:),TMa(:,:)
  integer, allocatable           :: brr(:)
  real(kind(1d0)), allocatable   :: PrjMat(:,:),Xip(:,:),OB(:,:),iOB(:,:),OukEB(:,:,:,:),OXiEB(:,:)
  DoublePrecision, allocatable   :: arr(:)

  !.. For the lapack subroutines, and the unit of files
  real(kind(1d0)), allocatable   :: WORK(:)
  integer                        :: LWORK,INFO,uid
  integer, allocatable           :: iWORK(:)


    !.. Local parameters for READ_UKRMOLP_INTS
  logical, parameter      :: master_writes_to_stdout = .false.  ! all processe
  logical, parameter      :: allow_shared_memory     = .false.  ! not yet prop
  integer*8 :: nfte  = 10
  integer*8 :: nfti  =  1
  integer*8 :: lembf =  0
  integer*8 :: nocsf =  0
  integer*8 :: nfta  =  6
  integer*8 :: isymtp=  0
  integer*8 :: nsym1
  integer*8 :: iposit=  0
  logical*8 :: qmoln    = .false.
  real(kind=wp)           :: scalem
  character(len=line_len) :: ukrmolp_header
  character(len=120)      :: name
  integer*8                 :: nalm
  integer*8        :: nint1e,nint2e



  
  type, public :: BasisElementInfo
     !.. l,m,irr.  nr is the number of radial bspline or the number of molecular orbital counting from the first orbital
     integer      :: l,m,nr,irr,n_lm
     integer      :: irs     ! this is the index between the subset Loc, stc, or pwc
     integer      :: i_sf  !spin free element counter LOC-SRC-PWC iIrrep1, 2, etc.
     integer      :: cr !this is the index of the external Bspline to which it is identical
     integer      :: RCukLoc=0
     integer      :: RCukLocRel=0
     integer      :: spin               !spin  +1 or -1
     character*2  :: OType
     logical      :: updownS=.false.            !two different values of spin
     logical      :: MOtype=.false.             !Molecular orbita
     logical      :: RCtype=.false.             !Molecular orbital
     logical      :: OBtype=.false.            !orthonormalized Bspline
     logical      :: EBtype=.false.            !External Bspline
  End type BasisElementInfo

  type, public :: STEXElementInfo
     type(BasisElementInfo) :: ihole,iexci
     integer      :: l,m,nr,irr         
     integer      :: spin               !spin  +1 or -1
     logical      :: updownS=.false.            !two different values of spin
     logical      :: MOtype=.false.            !Molecular orbital
     logical      :: OBtype=.false.            !orthonormalized Bspline
     logical      :: EBtype=.false.            !External Bspline
  End type STEXElementInfo


end module ModuleBasisJUAN
