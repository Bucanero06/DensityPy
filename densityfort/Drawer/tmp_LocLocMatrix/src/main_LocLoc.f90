program read_moints
  use mpi_mod !UKRMOL
  !  use ukrmol_interface, ONLY : read_ukrmolp_ints,&
  !    data_file_obj_id, line_len,GET_INTEGRAL,READ_UKRMOLP_BASIS 
  use ukrmol_interface 
  use precisn !UKRMOL
  implicit none
  integer :: indx
  integer :: iexerror
  integer :: irrdim
  integer :: nfte=10
  integer :: nfti=1
  integer :: lembf=0
  integer :: nocsf=0
  integer :: nfta=6
  integer :: isymtp=0
  integer :: iposit=0
  integer :: nint1e,nint2e,ia,ib,ic,id,i,j,k,l,i1,i2,i3,i4
  integer :: ind(1), two_ind(1:2,1)
  real(kind=8) ivalue,ivalue2
  character(len=132)  name2
  integer  n_AO,LWORK,INFO
  !  real(kind=cfp), allocatable :: cf(:,:),Hirr(:,:),En(:)
  real(kind=8), allocatable :: cf(:,:),Hirr(:,:),En(:)
  integer, allocatable    :: irr_series(:),irr_counter(:)
  real*8, allocatable     :: WORK(:)
  logical, parameter      :: master_writes_to_stdout = .false.  ! all processe
  logical, parameter      :: allow_shared_memory     = .false.  ! not yet prop
  real(kind=wp)            :: scalem
  character(len=line_len) :: ukrmolp_header
  character(len=120)      :: name
  logical                 :: qmoln = .false., tmpexist = .false.
  integer                 :: nalm,Nelec,Nbasis,Norb,Ndim
  integer, allocatable    :: MO(:),MeO(:),hA(:),xeP(:),hB(:),xeQ(:),shA(:),sxeP(:),shB(:),sxeQ(:)
  integer, external       :: kron 
  real*8                  :: NucRep,H0d










end program read_moints



