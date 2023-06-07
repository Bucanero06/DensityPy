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
Module ModuleUKRmolInterface

  use, intrinsic :: ISO_FORTRAN_ENV

  !.. UKRMOL module
  use precisn, only : wp, cfp  
  use basis_data_generic_mod 

  !.. Level 0.0
  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleConstants
  !.. Level 0.1
  use ModuleBspline
  use ModuleGroups
  !.. Level 1
  use ModuleXlm
  use ModuleMolecularGeometry
  !.. Level 2
  use ModuleIntegrals

  implicit none

  private

  !*** This value is now read from ukrmol data and stored in Class_UKRmol_orbitals%Lmax_mp
  !*** fix dependence of Class_UKRmol_Moments on this parameter
  integer         , public , parameter :: UKRMOL_LMAX_MULTIPOLE = 20
  character(len=*), public , parameter :: UKRMOL_LOG_FILE       = "log_file_UK"

  !****
  integer, parameter :: N_EXT_BSPLINE_INTERVALS = 8
  !****


  !.. Local parameters for READ_UKRMOLP_INTS
  integer*8                     , parameter :: UKRMOL_NFTI  =  1
  integer                       , parameter :: UKRMOL_LINE_LEN = 207
  character(len=UKRMOL_LINE_LEN), parameter :: UKRMOL_data_file_obj_id = "DATA FILE version 1.0"
  logical                                   :: UKRMOL_BASIS_INITIALIZED = .FALSE.

  !.. This class is used to parse src bspline component to the SRC_PWC soubrutine
  type, private :: ClassSRCParam
     !.. Expansion coefficients of the Ukrmol orbitals. The second index correspond to the orbital.
     !   The first one runs over the atomic and bspline basis elements
     real(kind(1d0)), allocatable :: cf(:,:)
     !.. Position of the first (ni_bto) and last (nf_bto) bto in the basis. nlr_bs is the last radial bspline
     !   includded in UKRmol.
     integer              :: ni_bto, nf_bto, nlr_bs     
     !.. For each positon ni_bto<=i<=nf_bto, l_bto(i), m_bto(i) and nr_bto(i)<=nlr_bs give
     !   the angular and radial numbers of the bsplines.
     integer, allocatable :: l_bto(:), m_bto(:), nr_bto(:)
     !.. Order of the bsplines and first index of the external Bsplines.
     integer              :: BsOrder, ifEBs    
  end type ClassSRCParam



  
  type, private :: Class_UKRmol_orbitals
     !
     private
     integer              :: SymmetryIdentifier
     integer              :: nIrreps
     !.. Lmax_mp is associated to the property integrals max value set in Ukrmol
     integer              :: lmax,Lmax_mp,nlmq  !Lmax_mpx
     !.. Dimension of the primitive basis of ukrmol (atomic orbitals plus non orthogonalized bsplines)
     integer              :: NprimUkrmol
     !
     !.. BSPLINE DATA
     !.. Common parameters
     integer              :: order
     !
     !.. UKRMOL+ BSPLINE parameters
     integer              :: nBs_UKRmol
     integer              :: nnodes_UKRmol
     real(kind(1d0))      :: BsA_UKRmol
     real(kind(1d0))      :: BsB_UKRmol
     !
     !.. BSPLINE parameters of the atomic UKRMOL basis
     !.. The first and last coefficient indexes corresponding to the bto orbitals (the last is the basis dimension)
     integer              :: ni_bto, nf_bto, nlr_bs
     !.. Values of l, m and radial index nr associated to a given coefficient in the expansion of SRC
     integer, allocatable :: l_bto(:), m_bto(:), nr_bto(:) 
     !
     !.. ASTRA BSPLINE parameters
     integer :: nNodesAstra

     integer :: nBsplines
     integer :: iBs_LastExt !.. Index of the last external Bspline, i.e., the (nBsplines-1) in the new set.

     real(kind(1d0)) :: Bs_Rmin_Astra, Bs_Rmax_Astra
     integer         :: nExtBs_Astra
     integer         :: iBs_FirstExt !.. Index of the first external Bspline, i.e., the (2*order-1) in the new set.
     !
     integer, allocatable :: nLOC(:) ! Number of PI     or target    orbitals per irrep
     integer, allocatable :: nSRC(:) ! Number of hybrid or continuum orbitals per irrep
     integer, allocatable :: nTOT(:) ! Total Number of               orbitals per irrep
     integer, allocatable :: iLOC(:) ! Number of PI     or target    orbitals per irrep
     integer, allocatable :: iSRC(:) ! Number of hybrid or continuum orbitals per irrep
     integer, allocatable :: sTOT(:) ! Total Number of orbitals up to a given     irrep
     integer, allocatable :: iseries(:)
     integer, allocatable :: indexASTRA(:)

     !
     real(kind(1d0)), allocatable :: cf(:,:)
     !
   contains
     !
     generic, public  :: Init        => Class_UKRmol_orbitals_InitFromData 
     generic, public  :: Free        => Class_UKRmol_orbitals_Free
     generic, public  :: Get_Sym     => Class_UKRmol_orbitals_Get_Sym
     generic, public  :: Get_lmax    => Class_UKRmol_orbitals_Get_lmax
     generic, public  :: Get_Lmax_mp => Class_UKRmol_orbitals_Get_Lmax_mp
     generic, public  :: Get_order   => Class_UKRmol_orbitals_Get_order
     generic, public  :: Get_nIrr    => Class_UKRmol_orbitals_Get_nIrr
     generic, public  :: Get_nLOC    => Class_UKRmol_orbitals_Get_nLOC
     generic, public  :: Get_nSRC    => Class_UKRmol_orbitals_Get_nSRC
     generic, public  :: Get_nTOT    => Class_UKRmol_orbitals_Get_nTOT
     generic, public  :: Get_iLOC    => Class_UKRmol_orbitals_Get_iLOC
     generic, public  :: Get_iSRC    => Class_UKRmol_orbitals_Get_iSRC
     generic, public  :: Get_iseries => Class_UKRmol_orbitals_Get_iseries
     generic, public  :: Get_sTOT    => Class_UKRmol_orbitals_Get_sTOT
     generic, public  :: Get_nEBs    => Class_UKRmol_orbitals_Get_nEBs

     generic, public  :: Get_l_bto       => Class_UKRmol_orbitals_Get_l_bto
     generic, public  :: Get_nr_bto      => Class_UKRmol_orbitals_Get_nr_bto
     generic, public  :: Get_ni_bto      => Class_UKRmol_orbitals_Get_ni_bto
     generic, public  :: Get_nf_bto      => Class_UKRmol_orbitals_Get_nf_bto
     generic, public  :: Get_nlr_bs      => Class_UKRmol_orbitals_Get_nlr_bs
     generic, public  :: Get_cf          => Class_UKRmol_orbitals_Get_cf
     generic, public  :: Get_NprimUKRmol => Class_UKRmol_orbitals_Get_NprimUKRmol
     
     generic, public  :: Get_nExBsAstra  => Class_UKRmol_orbitals_Get_nExBsAstra
     generic, public  :: Get_BsFirst     => Class_UKRmol_orbitals_Get_BsFirst
     generic, public  :: Get_BsLast      => Class_UKRmol_orbitals_Get_BsLast
     generic, public  :: Get_nNodesAstra => Class_UKRmol_orbitals_Get_nNodesAstra
     generic, public  :: Get_BsRminAstra => Class_UKRmol_orbitals_Get_BsRminAstra
     generic, public  :: Get_BsRmaxAstra => Class_UKRmol_orbitals_Get_BsRmaxAstra
     generic, public  :: Get_nBs         => Class_UKRmol_orbitals_Get_nBs
     generic, public  :: Get_Bsorder     => Class_UKRmol_orbitals_Get_Bsorder
     generic, public  :: Get_BsA_UKRmol  => Class_UKRmol_orbitals_Get_BsA_UKRmol
     generic, public  :: Get_BsB_UKRmol  => Class_UKRmol_orbitals_Get_BsB_UKRmol
     generic, public  :: Get_nBs_UKRmol  => Class_UKRmol_orbitals_Get_nBs_UKRmol
     generic, public  :: Get_IntIndex    => Class_UKRmol_orbitals_Get_IntIndex
     !
     procedure, private :: Class_UKRmol_orbitals_InitFromData
     procedure, private :: Class_UKRmol_orbitals_Free
     procedure, private :: Class_UKRmol_orbitals_Get_Sym
     procedure, private :: Class_UKRmol_orbitals_Get_lmax
     procedure, private :: Class_UKRmol_orbitals_Get_Lmax_mp
     procedure, private :: Class_UKRmol_orbitals_Get_order
     procedure, private :: Class_UKRmol_orbitals_Get_nIrr
     procedure, private :: Class_UKRmol_orbitals_Get_nLOC
     procedure, private :: Class_UKRmol_orbitals_Get_nSRC
     procedure, private :: Class_UKRmol_orbitals_Get_nTOT
     procedure, private :: Class_UKRmol_orbitals_Get_iLOC
     procedure, private :: Class_UKRmol_orbitals_Get_iSRC
     procedure, private :: Class_UKRmol_orbitals_Get_iseries
     procedure, private :: Class_UKRmol_orbitals_Get_sTOT
     procedure, private :: Class_UKRmol_orbitals_Get_nEBs

     procedure, private :: Class_UKRmol_orbitals_Get_l_bto
     procedure, private :: Class_UKRmol_orbitals_Get_nr_bto
     procedure, private :: Class_UKRmol_orbitals_Get_ni_bto
     procedure, private :: Class_UKRmol_orbitals_Get_nf_bto
     procedure, private :: Class_UKRmol_orbitals_Get_nlr_bs
     procedure, private :: Class_UKRmol_orbitals_Get_cf
     procedure, private :: Class_UKRmol_orbitals_Get_NprimUKRmol
     
     procedure, private :: Class_UKRmol_orbitals_Get_nExBsAstra
     procedure, private :: Class_UKRmol_orbitals_Get_BsFirst
     procedure, private :: Class_UKRmol_orbitals_Get_BsLast
     procedure, private :: Class_UKRmol_orbitals_Get_nNodesAstra
     procedure, private :: Class_UKRmol_orbitals_Get_BsRminAstra
     procedure, private :: Class_UKRmol_orbitals_Get_BsRmaxAstra
     procedure, private :: Class_UKRmol_orbitals_Get_nBs
     procedure, private :: Class_UKRmol_orbitals_Get_Bsorder
     procedure, private :: Class_UKRmol_orbitals_Get_BsA_UKRmol
     procedure, private :: Class_UKRmol_orbitals_Get_BsB_UKRmol
     procedure, private :: Class_UKRmol_orbitals_Get_nBs_UKRmol
     procedure, private :: Class_UKRmol_orbitals_Get_IntIndex

     !
  end type Class_UKRmol_orbitals
  !.. SINGLETON
  type( Class_UKRmol_orbitals), public :: UKRmol_orbitals
  type( ClassSRCParam )       , public :: srcBspar     


  type, public :: Class_UKRmol_Moments
     !
     private
     integer  :: Nmoments 
     integer              :: i
     real(kind(1d0)), allocatable :: A(:) 
     !
   contains
     !
     generic  , public  :: Init   => Class_UKRmol_Moments_Init
     generic  , public  :: lm2ind => Class_UKRmol_Moments_lm2ind
     generic  , public  :: ind2l  => Class_UKRmol_Moments_ind2l
     generic  , public  :: ind2m  => Class_UKRmol_Moments_ind2m
     !
     procedure, private :: Class_UKRmol_Moments_Init
     procedure, private :: Class_UKRmol_Moments_lm2ind
     procedure, private :: Class_UKRmol_Moments_ind2l
     procedure, private :: Class_UKRmol_Moments_ind2m
     !
  end type Class_UKRmol_Moments


  !.. Identifying labels of the one-body integrals in UKRmol interface
  character(len=*), parameter :: UKRMOL_OVERLAP     = 'Overlap integrals'
  character(len=*), parameter :: UKRMOL_KINETIC     = 'Kinetic energy integrals'
  character(len=*), parameter :: UKRMOL_MOMENTS  = 'Property integrals'
  character(len=*), parameter :: UKRMOL_NUCLEAR_ATT = 'Nuclear attraction integrals'
  character(len=*), parameter :: UKRMOL_HAMILTONIAN = 'One electron Hamiltonian integrals'

  !.. Column index, which coincides with the Integral ID, for the UKRmol 
  !   One-body integrals as returned by UKRMol_Get_Integral_Index
  integer*8 :: UKRMOL_1BINT_ID_OVERLAP
  integer*8 :: UKRMOL_1BINT_ID_KINETIC
  integer*8 :: UKRMOL_1BINT_ID_MOMENTS
  integer*8 :: UKRMOL_1BINT_ID_NUCLEAR_ATT
  integer*8 :: UKRMOL_1BINT_ID_HAMILTONIAN


  !*** MUST READ INFO OF BSPLINES

  public :: UKRMol_Init_Integrals

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
  end function Power


  subroutine Init_UKRMOLP_BASIS( IntegralFile )
    
    use mpi_mod
    use ukrmol_interface
    character(len=*), optional, intent(in) :: IntegralFile
    
    logical*8          , parameter :: MASTER_WRITES_TO_STDOUT = .False.
    character(len=UKRMOL_LINE_LEN) :: ukrmolp_header
    integer :: uid, iostat
    character(len=500) :: iomsg
    
    if( UKRMOL_BASIS_INITIALIZED ) return
    call mpi_mod_start( master_writes_to_stdout )
    if(present(IntegralFile))then
       open(newunit= uid                        , &
            file   = trim(adjustl(IntegralFile)), &
            form   ='unformatted'               , &
            status ='old'                       , &
            action ='read'                      , &
            access ='stream'                    , &
            iostat = iostat                     , &
            iomsg  = iomsg                      )
       if(iostat/=0)then
          call ErrorMessage("Error opening : "//trim(IntegralFile))
          call ErrorMessage("iomsg         : "//trim(iomsg))
          stop
       endif
       read (uid) ukrmolp_header
       close(uid)
    else
       open( UKRMOL_NFTI , form='unformatted', access='stream')
       read( UKRMOL_NFTI) ukrmolp_header
       close(UKRMOL_NFTI)
    endif
    !.. Check header conformity
    if ( ukrmolp_header .ne. UKRMOL_data_file_obj_id ) then
       WRITE(*,'(a,";",a)') ukrmolp_header, UKRMOL_data_file_obj_id
       STOP "An unknown header on the integrals file."
    end if
    call READ_UKRMOLP_BASIS(UKRMOL_NFTI,IntegralFile)
    UKRMOL_BASIS_INITIALIZED = .TRUE.

  end subroutine Init_UKRMOLP_BASIS


  subroutine UKRMol_Init_Integrals( QCDir, MoldenFile, IntegralFile )

    use ModuleIntegrals
    use ModuleOrbitalBasis
    use ukrmol_interface 
    use precisn

    character(len=*)          , intent(in) :: QCDir
    character(len=*)          , intent(in) :: MoldenFile
    character(len=*), optional, intent(in) :: IntegralFile

    integer*8, parameter :: NFTE  = 10
    integer*8, parameter :: LEMBF =  0
    integer*8, parameter :: NOCSF =  0
    integer*8, parameter :: NFTA  =  6
    integer*8, parameter :: ISYMTP=  0
    integer*8, parameter :: IPOSIT=  0
    character(len=120), parameter :: NAME = "MYNICENAME"
    logical*8, parameter :: QMOLN = .false.
    integer*8, parameter :: IWRITE = 3

    integer*8          :: nint1e, nint2e
    integer*8          :: irrdim
    integer*8          :: irr_series(8)
    real(kind=wp)      :: scalem = 1.0_cfp
    integer*8          :: nalm
    integer            :: iIrr

    type(ClassXlm)     :: Xlm
    type(ClassIrrep), pointer :: OpIrrep
    type(ClassIrrep), pointer :: IrrList(:)
    integer                   :: nIrreps

    integer            :: lmax, order, BandWidth, I1B_ID
    integer, allocatable, target :: nLOCv(:), nSRCv(:), iLOCv(:), iSRCv(:)


    !*** 
    ! Number of external radial numerical functions (e.g., Bsplines)
    !    integer, parameter :: NextF = 50
    integer                      :: nExtF
    procedure(D2DFun) , pointer  :: fptr,gptr
    integer                      :: NBsplines,ifBs,ilBs
    real(kind(1d0)), parameter   :: RNUC = 1.0369872111d0
    real(kind(1d0)), parameter   :: ZNit = 7.d0
    character(len=7)             :: Int1BLabel, Int2BLabel

    integer                      :: maxsize, l, m
    integer*8                    :: UKR_ID
    type(ClassIrrep), pointer    :: IrrBra, IrrKet
    real(kind(1d0)), allocatable :: dmat(:,:)
    character(len=5)             :: BlockLabel

    integer, pointer             :: ivec1(:), ivec2(:), nvec1(:), nvec2(:)

    type(ClassBSpline) :: BsSet
    integer            :: Nnodes, BsOrder, iNode
    integer, allocatable :: ASTRAIrr_from_DALTONIrr(:)
    real(kind(1d0))    :: BsRmin, BsRmax
    real(kind(1d0)), allocatable :: Bsgrid(:)

    !.. This object parses the molecular geometry from the molden file.
    type(ClassMolecularGeometry) :: MolGeom



    call Init_UKRMOLP_BASIS(IntegralFile)

    irrdim = UKRmol_orbitals%Get_nIrr()
    lmax   = UKRmol_orbitals%Get_lmax()
    order  = UKRmol_orbitals%Get_order()
    BandWidth = 2 * order - 1
    allocate(nLOCv(irrdim), nSRCv(irrdim), iLOCv(irrdim), iSRCv(irrdim) )

    iLOCv(1)=0
    do iIrr = 1, irrdim
       irr_series( iIrr ) = UKRmol_orbitals%Get_iseries( iIrr )
       nLOCv     ( iIrr ) = UKRmol_orbitals%Get_nLOC( iIrr )
       nSRCv     ( iIrr ) = UKRmol_orbitals%Get_nSRC( iIrr )
       iLOCv     ( iIrr ) = UKRmol_orbitals%Get_iLOC( iIrr )
       iSRCv     ( iIrr ) = UKRmol_orbitals%Get_iSRC( iIrr )
       ! if(iIrr>1) iLOCv(iIrr) = iLOCv(iIrr-1) + nLOCv(iIrr-1) + nSRCv(iIrr-1)

       ! iSRCv(iIrr) = iLOCv(iIrr)   + nLOCv(iIrr)
    enddo

    call READ_UKRMOLP_INTS( NFTE, UKRMOL_NFTI, LEMBF, nint1e, nint2e, NOCSF,&
         NFTA, ISYMTP, irrdim, irr_series, IPOSIT, scalem, NAME, nalm, QMOLN, &
         IntegralFile )
    call READ_UKRMOLP_PROPERTY_INTS( UKRMOL_NFTI, IWRITE, IntegralFile )

    !.. Determines the identifying index of the one-body integrals
    !..
    UKRMOL_1BINT_ID_OVERLAP     = UKRMol_Get_Integral_Index( UKRMOL_OVERLAP     )    
    UKRMOL_1BINT_ID_KINETIC     = UKRMol_Get_Integral_Index( UKRMOL_KINETIC     )    
    UKRMOL_1BINT_ID_MOMENTS     = UKRMol_Get_Integral_Index( UKRMOL_MOMENTS     )    
    UKRMOL_1BINT_ID_NUCLEAR_ATT = UKRMol_Get_Integral_Index( UKRMOL_NUCLEAR_ATT )    
    UKRMOL_1BINT_ID_HAMILTONIAN = UKRMol_Get_Integral_Index( UKRMOL_HAMILTONIAN )    

    write(*,*) "UKRMOL_1BINT_ID_MOMENTS", UKRMOL_1BINT_ID_MOMENTS

    !*** this has been shortened to keep only regular bsplines and facilitate comparison with Juan's code
    nExtF   = UKRmol_orbitals%nExtBs_Astra !***
    !write(*,*) nExtF
    !nExtF   = UKRmol_orbitals%nExtBs_Astra- (UKRmol_orbitals%order - 1)
    write(*,*) "nExtF",nExtF
    !pause
    ifBs    = UKRmol_orbitals%iBs_FirstExt
    ilBs    = UKRmol_orbitals%iBs_LastExt
    maxsize = max(maxval(nLOCv), maxval(nSRCv), nExtF)
    allocate(dMat(maxsize,maxsize))


    
    Nnodes = UKRmol_orbitals%nNodesAstra  !***
    BsOrder= UKRmol_orbitals%order        !***
    BsRmin = UKRmol_orbitals%Bs_Rmin_Astra!***
    BsRmax = UKRmol_orbitals%Bs_Rmax_Astra!***

    NBsplines = UKRmol_orbitals%nBsplines         !***


    ! nExtF   = UKRmol_orbitals%Get_nExBsAstra() !***
    ! ifBs    = UKRmol_orbitals%Get_BsFirst()
    ! ilBs    = UKRmol_orbitals%Get_BsLast()
    ! maxsize = max(maxval(nLOCv), maxval(nSRCv), nExtF)
    ! allocate(dMat(maxsize,maxsize))

    ! Nnodes = UKRmol_orbitals%Get_nNodesAstra()
    ! BsOrder= UKRmol_orbitals%Get_order()
    ! BsRmin = UKRmol_orbitals%Get_BsRminAstra()
    ! BsRmax = UKRmol_orbitals%Get_BsRmaxAstra()

    ! NBsplines = UKRmol_orbitals%Get_nBs()
    
    allocate(Bsgrid(Nnodes))
    do iNode = 1, Nnodes
       Bsgrid(iNode) = BsRmin + (BsRmax - BsRmin) / dble( Nnodes - 1 ) * dble( iNode - 1 )
    enddo

    call BsSet%Init(Nnodes,BsOrder,Bsgrid)
    !*** MUST BE GENERALIZED TO GENERAL PAIRS OF BRA AND KET ANGULAR MOMENTA
    !    allocate(BsMat(2*BsOrder-2,2*BsOrder-2))

    !*** USE NORMALIZATION CONSTANT THAT IS ALREADY GIVEN IN THE BSPLINE CLASS
    !.. Here we allocate the object src_Bspar which contains the information related to the bspline part
    !   of the src orbital. Te object is used to parse it to ComputeSRCPWCBlocks routine.
    srcBspar%ni_bto  = UKRmol_orbitals%ni_bto  !.. Starting position of the bsplines in the coefficient list
    srcBspar%nf_bto  = UKRmol_orbitals%nf_bto  !.. Final     position of the bsplines in the coefficient list
    srcBspar%nlr_bs  = UKRmol_orbitals%nlr_bs  !.. Last radial bspline includded in UKRmol.
    allocate(srcBspar%l_bto(srcBspar%ni_bto:srcBspar%nf_bto))   !.. Angular quantum number l
    allocate(srcBspar%m_bto(srcBspar%ni_bto:srcBspar%nf_bto))   !.. Angular quantum number m
    allocate(srcBspar%nr_bto(srcBspar%ni_bto:srcBspar%nf_bto))  !.. Radial  Bspline index  nr
    srcBspar%l_bto   = UKRmol_orbitals%l_bto
    srcBspar%m_bto   = UKRmol_orbitals%m_bto
    srcBspar%nr_bto  = UKRmol_orbitals%nr_bto
    allocate(srcBspar%cf(srcBspar%nf_bto,srcBspar%nf_bto))      !.. Coefficients.
    !*** THIS MUST BE WRONG, THEIR DIMENSIONS DO NOT MATCH
    srcBspar%cf      = UKRmol_orbitals%cf
    srcBspar%BsOrder = BsOrder
    srcBspar%ifEBs   = ifBs

    ! allocate(GlobalIntegral%Int_1B( I1B_N_TYPES - 1 + (2*lmax+1)**2 ))
    
    !    call GlobalIntegral%NucMoments%Init( lmax )                       !.. Allocate
    !    call GlobalIntegral%LocMoments%Init( lmax, nLOCv )                !.. Only allocate
    !    call GlobalIntegral%BsBlocks%Init( lmax, NBsplines, BsOrder )     !.. Allocate
    call GlobalIntegral%Init( lmax, nLOCv, nSRCv , NBsplines, BsOrder )

    call ReadDALTONIrreps( ASTRAIrr_from_DALTONIrr, QCDir )
    
    call OrbitalBasis%Init( lmax, nLOCv, nSRCv, BsSet )
    
    !*** TEST IF IT WORKS!
    call MolGeom%ParseMolden( MoldenFile )

    call GlobalIntegral%NucMoments%Set( MolGeom )              !.. Compute

    call GlobalIntegral%BsBlocks%Set( BsSet )                  !.. Compute


        
    !.. Converts the one-body integrals to the local format
    nIrreps =  GlobalGroup%GetnIrreps()
    IrrList => GlobalGroup%GetIrrepList()
    call Xlm%Init(0,0)  
    OpIrrep => Xlm%GetIrrep( GlobalGroup )
 

    !.. Overlap
    I1B_ID     = I1B_OVERLAP
    UKR_ID     = UKRMOL_1BINT_ID_OVERLAP
    Int1BLabel = I1B_ID_LIST( I1B_ID )
    call TransferOperator()


    !.. Kinetic Energy
    I1B_ID = I1B_KINETIC
    UKR_ID = UKRMOL_1BINT_ID_KINETIC
    Int1BLabel = I1B_ID_LIST( I1B_ID )
    call TransferOperator()


    !.. Nuclear Attraction
    I1B_ID = I1B_NUCLATT
    UKR_ID = UKRMOL_1BINT_ID_NUCLEAR_ATT
    Int1BLabel = I1B_ID_LIST( I1B_ID )
    call TransferOperator()

    !.. Hamiltonian
    I1B_ID = I1B_HAMILTO
    UKR_ID = UKRMOL_1BINT_ID_HAMILTONIAN
    Int1BLabel = I1B_ID_LIST( I1B_ID )
    call TransferOperator()

    
    !.. Multipoles, i.e., X_{lm}(\hat{r}) / r^{l+1}
    do l = 0, 2*lmax 
       do m = -l, l
          call Xlm%Init(l,m)
          OpIrrep => Xlm%GetIrrep( GlobalGroup )
          I1B_ID = GlobalIntegral%Get_Multipole_ID( l, m )
          Int1BLabel=I1B_ID_LIST(I1B_MULTIPO)
          write(Int1BLabel(4:),'(I0.4)') I1B_ID - I1B_MULTIPO + 1
          call TransferOperator()
       enddo
    end do

    !.. LOC-LOC Moments, i.e., the ME of \frac{4\pi}{2\ell+1} r^l X_{lm}(\hat{r}) 
    !   on the basis of the molecular orbitals
    !   LOC (MO) - LOC (MO) 
    !..
    call TransferMoments( lmax, nLOCv, iLOCv, IrrList )


    deallocate(dMat)
    
    !.. Bielectronic Integrals
    call Xlm%Init(0,0)
    OpIrrep => Xlm%GetIrrep( GlobalGroup )
    Int2BLabel = "bielectronic"
    call Transfer2BOperator()

    !.. Empties the UKRMOL+ variables
    call DESTROY_UKRMOLP()

    call MolGeom%Free()

    !*** This is commented since it has to be done after writing the matrix, otherwise is set as not initialized when readed
    !call GlobalIntegral%BsBlocks%Free( )


  contains


    subroutine TransferOperator()

      integer :: in1, in2, iBs1, iBs2


      !call GlobalIntegral%Int_1B(I1B_ID)%Init( Int1BLabel, &
      !     OpIrrep, lmax, nLOCv, nSRCv, nExtF, BandWidth )

      if(Int1BLabel(1:3)/="MP_")then

         !.. LOC (MO) - LOC (MO) 
         !BlockLabel = "LOC_LOC"
         BlockLabel = "MO_MO"
         nvec1 => nLOCv
         nvec2 => nLOCv
         ivec1 => iLOCv
         ivec2 => iLOCv
         dMat=0.d0
         !
         call TransferBlock()

         !.. LOC (MO) - SRC (HYB)
!         BlockLabel = "LOC_SRC"
         BlockLabel = "MO_HY"
         nvec1 => nLOCv
         nvec2 => nSRCv
         ivec1 => iLOCv
         ivec2 => iSRCv
         dMat=0.d0
         !
         call TransferBlock()

         !.. SRC (HYB)- SRC (HYB) 
         !BlockLabel = "SRC_SRC"
         BlockLabel = "HY_HY"
         nvec1 => nSRCv
         nvec2 => nSRCv
         ivec1 => iSRCv
         ivec2 => iSRCv
         dMat=0.d0
         !
         call TransferBlock()

      endif

      !.. SRC (HYB)- PWC (NUM) 
      !.. Fill the Bspline - Bspline matrix
      BlockLabel = "HY_BS"
      nvec1 => nSRCv
      ivec1 => iSRCv
      dMat=0.d0
      !
      call ComputeSRCPWCBlock( lmax, BsOrder, I1B_ID, nIrreps, IrrList, OpIrrep, nvec1, ivec1 )


      !.. PWC (NUM)- PWC (NUM) 
!!$      BlockLabel = "PWC_PWC"
      BlockLabel = "BS_BS"
      dMat=0.d0
      !
      call ComputePWCPWCBlock( nExtf, lmax, BsOrder, I1B_ID, nIrreps, IrrList, OpIrrep, ifBs )



    end subroutine TransferOperator


    subroutine Transfer2BOperator()
      implicit none

      integer   :: iIrr1, iIrr2, iIrr3, iIrr4
      integer   :: n1, n2, n3, n4, i, n
      integer*8 :: i1, i2, i3, i4
      real(kind(1d0)), allocatable :: D4mat(:,:,:,:)
      type(ClassIrrep), pointer :: IrrPtr
      
      integer            :: iBlock
      integer, parameter :: NBlocks = 4
      character(len=5)   :: BlockVec(1:4) = [ "LL_LL", "LL_LS", "LL_SS", "LS_LS" ], BlockLabel
      integer            :: nmat(nIrreps,4), imat(nIrreps,4)
      
      call GlobalIntegral%Int_2B%init( Int2BLabel, OpIrrep, nLOCv, nSRCv )
      
      BlockLoop : do iBlock = 1, 4
         
         BlockLabel = BlockVec(iBlock)
         do i = 1, 4
            n=i
            if(i>2)n=i+1
            if(BlockLabel(n:n)=="L")then
               nmat(:,i)=nLOCv
               imat(:,i)=iLOCv
            else
               nmat(:,i)=nSRCv
               imat(:,i)=iSRCv
            endif
         enddo
         
         do iIrr1 = 1, nIrreps
            !
            do iIrr2 = 1, nIrreps
               !
               do iIrr3 = 1, nIrreps
                  !
                  !.. Chose iIrr4 compatible with the symmetry.
                  IrrPtr =>  IrrList(iIrr1) * IrrList(iIrr2) * IrrList(iIrr3)
                  iIrr4 = GlobalGroup%GetIrrepIndex( IrrPtr )
                  
                  if( GlobalIntegral%Int_2B%CheckIrrepOrder(BlockLabel, iIrr1, iIrr2, iIrr3, iIrr4 ) )then
                     
                     allocate( D4mat( nmat(iIrr1,1), nmat(iIrr2,2), nmat(iIrr3,3), nmat(iIrr4,4) ) )
                     do n1 = 1, nmat( iIrr1, 1 )
                        i1 = imat( iIrr1, 1 ) + n1
                        do n2 = 1, nmat( iIrr2, 2 )
                           i2 = imat( iIrr2, 2 ) + n2
                           do n3 = 1, nmat( iIrr3, 3 )
                              i3 = imat( iIrr3, 3 ) + n3
                              do n4 = 1, nmat( iIrr4, 4 )
                                 i4 = imat( iIrr4, 4 ) + n4
                                 d4Mat(n1,n2,n3,n4) = GET_2B_INTEGRAL_ASTRA(i1, i2, i3, i4)
                              enddo
                           enddo
                        enddo
                     enddo
                     
                     call GlobalIntegral%Int_2B%Set( BlockLabel, iIrr1, iIrr2, iIrr3, iIrr4, d4Mat )
                     
                     deallocate( D4mat )
                     
                  endif
                  
               enddo
            enddo
         enddo
      enddo BlockLoop


    end subroutine Transfer2BOperator


    subroutine TransferBlock()
      implicit none
      integer :: i, j, iIrrBra, iIrrKet, iInt, jInt
      integer*8 :: diInt, djInt, l_mp, m_mp
      ! write(*,*)
      ! write(*,*) "====="//BlockLabel//"=============================",UKR_ID
      do iIrrBra = 1, nIrreps
         IrrBra => IrrList( iIrrBra ) 
         IrrKet => IrrBra * OpIrrep
         iIrrKet = GlobalGroup%GetIrrepIndex( IrrKet )
         dMat = 0.d0
         !***  MUST MAKE SURE THE IRREP ORDER IN UKRMOL+ IS THE SAME AS IN MODULEGROUPS
         do i = 1, nvec1( iIrrBra )
            iInt = ivec1( iIrrBra ) + i
            !write(*,"(i5)",advance="no")iInt
            do j = 1, nvec2( iIrrKet )
               jInt = ivec2( iIrrKet ) + j
               diInt = iInt
               djInt = jInt
               dMat(i,j) = GET_1B_INTEGRAL_ASTRA( diInt, djInt, UKR_ID)
               !write(*,"(x,f8.2)",advance="no") dMat(i,j) 
            enddo
         enddo
         call GlobalIntegral%Int_1B(I1B_ID)%Set( BlockLabel, iIrrBra, iIrrKet, dMat, nvec1( iIrrBra ), nvec2( iIrrKet ) )
         !write(*,*) iIrrBra, sum(abs(dMat(1:nvec1(iIrrBra),1:nvec2(iIrrKet))))
         !pause
         !write(*,*) iIrrBra,iIrrKet,dMat(1,1)
      enddo
    end subroutine TransferBlock

  end subroutine UKRMol_Init_Integrals


  subroutine TransferMoments( lmax, nLOCv, iLOCv, IrrList )
    use ukrmol_interface 

    implicit none
    integer                  , intent(in) :: lmax, nLOCv(:), iLOCv(:)
    type(ClassIrrep), pointer, intent(in) :: IrrList(:)

    integer                               :: iIrrBra, iIrrKet, l, m, im, ilm
    integer                               :: i, j, iInt, jInt, nIrreps
    integer*8                             :: diInt, djInt, l_mp, m_mp
    type(ClassIrrep)          , pointer   :: IrrBra, IrrKet, ProdIrrep
    type(ClassXlmSymmetricSet), pointer   :: XlmSymSet
    integer        , allocatable          :: mlist(:)
    real(kind(1d0)), allocatable          :: dMat(:,:)
    real(kind(1d0))                       :: factor

    i=maxval(nLOCv)
    
    nIrreps = GlobalGroup%GetnIrreps()

    do iIrrBra = 1, nIrreps
       IrrBra => IrrList( iIrrBra ) 

       do iIrrKet = 1, nIrreps
          IrrKet => IrrList( iIrrKet ) 
          ProdIrrep => IrrBra * IrrKet 
          XlmSymSet => GlobalXlmSet%GetSymSet( ProdIrrep )   

          allocate(dMat(nLOCv(iIrrBra),nLOCv(iIrrKet)))
          dMat=0.d0

          do l = 0, 2*lmax
             call XlmSymSet%GetMList( l, mlist )
             if(.not.allocated(mlist))cycle
             !
             do im = 1, size(mlist)
                m = mlist(im)
                ilm = GlobalIntegral%lmPairIndex(l,m)

                l_mp = l
                m_mp = m
                factor = sqrt(4.d0 * PI / dble(2 * l + 1 ))

                dMat = 0.d0
                do i = 1, nLOCv( iIrrBra )
                   iInt = iLOCv( iIrrBra ) + i
                   do j = 1, nLOCv( iIrrKet )
                      jInt = iLOCv( iIrrKet ) + j
                      diInt = iInt
                      djInt = jInt
                      dMat(i,j) = factor * GET_1B_INTEGRAL_ASTRA( diInt, djInt, UKRMOL_1BINT_ID_MOMENTS, l_mp, m_mp )
                      
                   enddo
                enddo
                call GlobalIntegral%LocMoments%Set( iIrrBra, iIrrKet, ilm, dMat )

             enddo
          enddo

          deallocate(dMat)

       enddo
    enddo
  end subroutine TransferMoments


  subroutine ComputeSRCPWCBlock( lmax, BsOrder, I1B_ID, nIrreps, IrrList, OpIrrep, nvec1,ivec1)
    implicit none
    integer                  , intent(in) :: lmax, BsOrder, I1B_ID, nIrreps
    integer                  , intent(in) :: nvec1(:), ivec1(:)
    type(ClassIrrep), pointer, intent(in) :: IrrList(:), OpIrrep
    type(ClassXlmSymmetricSet), pointer   :: XlmSymSet

    type(ClassIrrep), pointer             :: IrrBra, IrrKet
    character(len=*), parameter           :: BlockLabel = "HY_BS"

    integer              :: i, j
    integer              :: iIrrBra, iIrrKet, iInt
    integer              :: l, m, lp, mp, lk, mk, ilm, im
    integer              :: in1, in2, jn1, iBs1, iBs2
    integer, allocatable :: mlist(:)
    real(kind(1d0))      :: rvalue, rvalue2, factor1
    real(kind(1d0)), allocatable :: dMat(:,:)

    !write(*,*) "====="//BlockLabel//"============================="

    if( I1B_ID >= I1B_MULTIPO ) call GlobalIntegral%PairTolm( I1B_ID - I1B_MULTIPO + 1, lk, mk )
    allocate(dMat(maxval(nvec1),BsOrder - 1))

    do iIrrBra = 1, nIrreps
       IrrBra => IrrList( iIrrBra ) 
       IrrKet => IrrBra * OpIrrep
       iIrrKet = GlobalGroup%GetIrrepIndex( IrrKet )
       XlmSymSet => GlobalXlmSet%GetSymSet( IrrKet )   

       do l = 0, lmax
          call XlmSymSet%GetMList( l, mlist )
          if(.not.allocated(mlist))cycle
          !
          do im = 1, size(mlist)
             m = mlist(im)
             ilm = GlobalIntegral%lmPairIndex(l,m)
             !
             dMat = 0.d0
             do i = 1, nvec1( iIrrBra )
                iInt = ivec1( iIrrBra ) + i
                !
                do in2  = BsOrder ,2*BsOrder - 2
                   iBs2 = in2 + BsOrder - 1
                   j    = in2 - BsOrder + 1
                   rvalue = 0.d0
                   !
                   SRCBSExpansionLoop: do jn1  = srcBspar%ni_bto, srcBspar%nf_bto
                      !.. Transform the bs index to the external indexing
                      iBs1 = srcBspar%nr_bto(jn1) - srcBspar%nlr_bs + 2 * BsOrder - 2
                      in1  = iBs1 - BsOrder + 1
                      lp   = srcBspar%l_bto(jn1) 
                      mp   = srcBspar%m_bto(jn1)
                      if(abs(in1-in2)>=BsOrder)cycle
                      if(abs(srcBspar%cf(jn1,iInt))<= 1.0E-15 )cycle
                      select case( I1B_ID )
                      case( I1B_OVERLAP )
                         if((l.ne.lp).or.(m.ne.mp))cycle
                         rvalue = rvalue + srcBspar%cf(jn1,iInt) * GlobalIntegral%BsBlocks%GetOverlap(iBs1, iBs2)
                      case( I1B_KINETIC )
                         if((l.ne.lp).or.(m.ne.mp))cycle
                         rvalue2 =  GlobalIntegral%BsBlocks%GetSecondDeriv(iBs1, iBs2)
                         rvalue = rvalue +  srcBspar%cf(jn1,iInt) * &
                              (-0.5d0)*(rvalue2 - l*(l+1)*GlobalIntegral%BsBlocks%GetCentrifugal(iBs1, iBs2))
                      case( I1B_NUCLATT )
                         do lk = abs (l - lp ), l + lp
                            factor1 = 0.d0
                            do mk = -lk, lk
                               factor1 = factor1 +  &
                                    ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                                    GlobalIntegral%NucMoments%GetMoments( lk, mk )
                            end do
                            rvalue = rvalue + srcBspar%cf(jn1,iInt) * factor1 * &
                                 GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, lk)       
                         end do
                      case( I1B_HAMILTO )
                         if((l.eq.lp).and.(m.eq.mp))then
                            rvalue2 = GlobalIntegral%BsBlocks%GetSecondDeriv(iBs1, iBs2)
                            rvalue2 = (-0.5d0)*(rvalue2 - &
                                 l*(l+1)*GlobalIntegral%BsBlocks%GetCentrifugal(iBs1, iBs2))
                         else
                            rvalue2 = 0.0d0
                         endif

                         !.. Since the multipolar interaction will be added separately
                         !   later on, there is no need to have it here, and hence
                         !   only the monopolar nuclear attraction will be added.
                         !   The electron repulsion, which typically offsets the nuclear
                         !   attraction, will be added later, when computing the h-tilde
                         !   effective hamiltonian.
                         !..
                         do lk = abs (l - lp ), l + lp
                            factor1 = 0.d0
                            do mk = -lk, lk
                               factor1 = factor1 +  &
                                    ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                                    GlobalIntegral%NucMoments%GetMoments( lk, mk )
                            end do
                            rvalue2 = rvalue2 +  factor1 * &
                                 GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, lk)       
                         end do

                         !.. Add monopolar nuclear attraction term
                         ! if(lp==l.and.mp==m) &
                         !      rvalue2 = rvalue2 + &
                         !      ThreeXlmIntegral( l, m, 0, 0, l, m ) * &
                         !      GlobalIntegral%NucMoments%GetMoments( 0, 0 ) * &
                         !      GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, 0)  
                         

                         rvalue = rvalue + srcBspar%cf(jn1,iInt) * rvalue2
                         !pause
                      case DEFAULT
                         if( I1B_ID >= I1B_MULTIPO )then
                            rvalue = rvalue + srcBspar%cf(jn1,iInt) * &
                                 ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                                 GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, lk)             
                         endif
                          
                      end select
                   end do SRCBSExpansionLoop
                   dMat( i , j ) = rvalue
                end do
             enddo
             call GlobalIntegral%Int_1B(I1B_ID)%Set( &
                  BlockLabel, iIrrBra, ilm, dMat, nvec1( iIrrBra ), BsOrder-1 )
          enddo
       enddo
    enddo
    deallocate(dMat)
    
  end subroutine ComputeSRCPWCBlock

  subroutine ComputePWCPWCBlock( nExtF, lmax, BsOrder, I1B_ID, nIrreps, IrrList, OpIrrep, ifBs)
    implicit none
    integer                  , intent(in) :: nExtf, lmax, ifBs, BsOrder, I1B_ID, nIrreps
    type(ClassIrrep), pointer, intent(in) :: IrrList(:), OpIrrep
    type(ClassXlmSymmetricSet), pointer   :: XlmSymSet,XlmSymSetp
    type(ClassIrrep), pointer             :: IrrBra, IrrKet
    character(len=*), parameter           :: BlockLabel = "BS_BS"

    integer                      :: iIrrBra, iIrrKet
    integer                      :: l, m, lp, mp, lk, mk, ilm, ilmp, im, imp
    integer                      :: in1, in2, iBs1, iBs2
    real(kind(1d0))              :: rvalue, rvalue2, factor1, func1, func2
    integer        , allocatable :: mlist(:),mplist(:)
    real(kind(1d0)), allocatable :: dMat(:,:)

    if( I1B_ID >= I1B_MULTIPO ) call GlobalIntegral%PairTolm( I1B_ID - I1B_MULTIPO + 1, lk, mk )
    allocate(dMat(nExtF,nExtF))

    !write(*,*) "====="//BlockLabel//"=============================",I1B_ID_LIST( I1B_ID )
    do iIrrBra = 1, nIrreps
       IrrBra => IrrList( iIrrBra ) 
       IrrKet => IrrBra * OpIrrep
       iIrrKet = GlobalGroup%GetIrrepIndex( IrrKet )
       XlmSymSet  => GlobalXlmSet%GetSymSet( IrrBra )   
       XlmSymSetp => GlobalXlmSet%GetSymSet( IrrKet )
       !
       do l = 0, lmax
          call XlmSymSet%GetMList( l, mlist )
          if(.not.allocated(mlist))cycle
          !
          do im = 1, size(mlist)
             m = mlist(im)
             ilm = GlobalIntegral%lmPairIndex(l,m)
             !
             do lp = 0, lmax
                call XlmSymSetp%GetMList( lp, mplist )
                if(.not.allocated(mplist))cycle
                !
                do imp = 1, size(mplist)
                   mp = mplist(imp)
                   ilmp = GlobalIntegral%lmPairIndex(lp,mp)
                   !
                   dMat = 0.d0
                   do in1 = 1, nExtF
                      iBs1 = in1 + ifBs - 1
                      do in2 = 1, nExtF
                         if(abs(in1-in2)>=BsOrder)cycle
                         iBs2 = in2 + ifBs - 1
                         rvalue = 0.d0
                         select case( I1B_ID )
                         case( I1B_OVERLAP )
                            if((l.ne.lp).or.(m.ne.mp))cycle
                            rvalue =  GlobalIntegral%BsBlocks%GetOverlap(iBs1, iBs2)
                         case( I1B_KINETIC )
                            if((l.ne.lp).or.(m.ne.mp))cycle
                            rvalue =  GlobalIntegral%BsBlocks%GetSecondDeriv(iBs1, iBs2)
                            rvalue =  (-0.5d0)*(rvalue - &
                                 l*(l+1)*GlobalIntegral%BsBlocks%GetCentrifugal(iBs1, iBs2))
                         case( I1B_NUCLATT )
                            do lk = abs (l - lp ), l + lp
                               factor1 = 0.d0
                               do mk = -lk, lk
                                  factor1 = factor1 +  &
                                       ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                                       GlobalIntegral%NucMoments%GetMoments( lk, mk )
                               end do
                               rvalue = rvalue + factor1 * &
                                    GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, lk)       
                            end do
                         case( I1B_HAMILTO )
                            if((l.eq.lp).and.(m.eq.mp))then
                               rvalue2 = GlobalIntegral%BsBlocks%GetSecondDeriv(iBs1, iBs2)
                               rvalue2 = (-0.5d0)*(rvalue2 - &
                                    l * (l+1) * GlobalIntegral%BsBlocks%GetCentrifugal(iBs1, iBs2))
                            else
                               rvalue2 = 0.0d0
                            endif


                            !.. Since the multipolar interaction will be added separately
                            !   later on, there is no need to have it here, and hence
                            !   only the monopolar nuclear attraction will be added.
                            !   The electron repulsion, which typically offsets the nuclear
                            !   attraction, will be added later, when computing the h-tilde
                            !   effective hamiltonian.
                            !..
                            ! do lk = abs (l - lp ), l + lp
                            !    factor1 = 0.d0
                            !    do mk = -lk, lk
                            !       factor1 = factor1 +  &
                            !            ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                            !            GlobalIntegral%NucMoments%GetMoments( lk, mk )
                            !    end do
                            !    rvalue2 = rvalue2 +  factor1 * &
                            !         GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, lk)       
                            ! end do

                            !.. Add monopolar nuclear attraction term
                            if(lp==l.and.mp==m) &
                                 rvalue2 = rvalue2 + &
                                 ThreeXlmIntegral( l, m, 0, 0, l, m ) * &
                                 GlobalIntegral%NucMoments%GetMoments( 0, 0 ) * &
                                 GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, 0)
                         
                             rvalue = rvalue + rvalue2
                         case DEFAULT
                            if( I1B_ID >= I1B_MULTIPO )then
                               rvalue = ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                                    GlobalIntegral%BsBlocks%GetMultipoles(iBs1, iBs2, lk)  
                            endif
                         end select
                         dMat( in1 , in2 ) = rvalue
                      end do
                   end do
                   !call GlobalIntegral%Int_1B(I1B_ID)%Set( BlockLabel, ilm, ilmp, dMat, nExtF, nExtF, BsOrder - 1 )
                   call GlobalIntegral%Int_1B(I1B_ID)%Set( BlockLabel, ilm, ilmp, dMat, nExtF, nExtF )
                   ! write(*,*) dMat
                   ! pause
                end do
             enddo
          end do
       enddo
    enddo

  end subroutine ComputePWCPWCBlock

  real(kind(1d0)) function GET_1B_INTEGRAL_ASTRA( a, b, Integral_ID, l, m ) result( qres )
    use ukrmol_interface 
    implicit none
    integer*8,         intent(in) :: a, b
    integer*8,           intent(in) :: Integral_ID
    integer*8, optional, intent(in) :: l, m
    integer*8 :: c, d
    c = UKRmol_orbitals%indexASTRA(a)
    d = UKRmol_orbitals%indexASTRA(b)
    if(present(l).and.present(m))then
       qres = GET_1B_INTEGRAL( c, d, Integral_ID, l, m )
    else
       qres = GET_1B_INTEGRAL( c, d, Integral_ID)
    endif
  end function GET_1B_INTEGRAL_ASTRA

  real(kind(1d0)) function GET_2B_INTEGRAL_ASTRA(a,b,c,d) result( qres )
    use ukrmol_interface 
    implicit none
    integer*8, intent(in) :: a,b,c,d
    integer*8 :: i1, i2, i3, i4
    i1 = UKRmol_orbitals%indexASTRA(a)
    i2 = UKRmol_orbitals%indexASTRA(b)
    i3 = UKRmol_orbitals%indexASTRA(c)
    i4 = UKRmol_orbitals%indexASTRA(d)
    qres = GET_2B_INTEGRAL(i1, i2, i3, i4)
  end function GET_2B_INTEGRAL_ASTRA
     
  subroutine Class_UKRmol_Moments_Init( self )
    class(Class_UKRmol_Moments), intent(out) :: self
    self%Nmoments = ( UKRMOL_LMAX_MULTIPOLE + 1 )**2
  end subroutine Class_UKRmol_Moments_Init

  integer function Class_UKRmol_Moments_lm2ind( self, l, m ) result( ind )
    class(Class_UKRmol_Moments), intent(in) :: self
    integer                       , intent(in) :: l, m
    ind = l**2+l+m+1
  end function Class_UKRmol_Moments_lm2ind

  integer function Class_UKRmol_Moments_ind2l( self, ind ) result( l )
    class(Class_UKRmol_Moments), intent(in) :: self
    integer                       , intent(in) :: ind
    l = int(sqrt(dble(ind)-1.d0))
  end function Class_UKRmol_Moments_ind2l

  integer function Class_UKRmol_Moments_ind2m( self, ind ) result( m )
    class(Class_UKRmol_Moments), intent(in) :: self
    integer                       , intent(in) :: ind
    integer :: l
    l = self%ind2l(ind)
    m = ind - l**2 - l - 1
  end function Class_UKRmol_Moments_ind2m

  subroutine Class_UKRmol_orbitals_Free( self )
    class( Class_UKRmol_orbitals ), intent(inout) :: self

    if(allocated(self%nLOC)) deallocate(self%nLOC)
    if(allocated(self%nSRC)) deallocate(self%nSRC)
    self%SymmetryIdentifier = -1
    self%nIrreps            =  0
    self%lmax               = -1
    self%order              =  0
  end subroutine Class_UKRmol_orbitals_Free


!!$  !*** To initialize the info about the number of orbitals in each symmetry
!!$  !    it is worth checking the function
!!$  !    READ_UKRMOP_BASIS in ukrmol_interface.f90
!!$  !    and the methods of molecular_orbital_basis_obj in molecular_basis_mod.f90
!!$  subroutine Class_UKRmol_orbitals_InitFromFile( self, LogFile )
!!$    class( Class_UKRmol_orbitals ), intent(out) :: self
!!$    character(len=*)              , intent(in)  :: LogFile
!!$
!!$    integer             :: uid, iostat, i, iIrrep
!!$    character(len=500)  :: iomsg
!!$    character(len=1000) :: line
!!$    call self%free()
!!$
!!$    open(newunit = uid       , &
!!$         file    = LogFile   , &
!!$         status  ='old'      , &
!!$         form    ='formatted', &
!!$         action  ='read'     , &
!!$         iostat  = iostat    , &
!!$         iomsg   = iomsg     )
!!$    if(iostat/=0)then
!!$       call ErrorMessage("File "//UKRMOL_LOG_FILE//" not found. Stop forced.")
!!$       stop
!!$    endif
!!$    do
!!$       read(uid,"(a)",iostat=iostat) line
!!$       if(iostat/=0)then
!!$          call ErrorMessage("orbital and sym parameters not found in "//UKRMOL_LOG_FILE)
!!$          stop
!!$       endif
!!$       line=adjustl(line)
!!$       i=index(line,"Point-group symmetry identifier: ")
!!$       if(i>0)exit
!!$    enddo
!!$    i=len_trim("Point-group symmetry identifier: ")
!!$    read(line(i+1:),*) self%SymmetryIdentifier
!!$    read(uid,"(a)") line
!!$    i=len_trim("Number of irreducible representations: ")
!!$    read(line(i+1:),*) self%nIrreps
!!$    read(uid,*)
!!$    read(uid,*)
!!$    read(uid,*)
!!$    allocate(self%nLOC(self%nIrreps))
!!$    read(uid,*) (self%nLOC(iIrrep),iIrrep=1,self%nIrreps)
!!$    read(uid,*)
!!$    allocate(self%nSRC(self%nIrreps))
!!$    read(uid,*) (self%nSRC(iIrrep),iIrrep=1,self%nIrreps)
!!$    close(uid)
!!$    allocate(self%nTOT(self%nIrreps))
!!$    self%nTOT = self%nLOC + self%nSRC
!!$    allocate(self%sTOT(0:self%nIrreps))
!!$    self%sTOT=0
!!$    do iIrrep = 1,self%nIrreps
!!$       self%sTOT(iIrrep) = self%sTOT(iIrrep-1) + self%nTOT(iIrrep)
!!$    enddo
!!$  end subroutine Class_UKRmol_orbitals_InitFromFile


  subroutine Class_UKRmol_orbitals_InitFromData( self, IntegralFile )

    use ukrmol_interface
    use const, only:  D2h_names, &
         D2_names , &
         C2h_names, &
         C2_names , &
         C2v_names, &
         Cs_names , &
         Ci_names , &
         C1_names

    implicit none

    class( Class_UKRmol_orbitals ), intent(out) :: self
    character(len=*), optional    , intent(in)  :: IntegralFile
    
    integer, parameter :: M_MULTIPLICITY_COLUMN = 5
    integer, parameter :: L_COLUMN              = 6

    integer                :: i, j, k, l, m, iIrrep, nIrreps, NprimUkrmol, iIrr, jIrr
    integer*8              :: in8, number_of_bto_shells, number_of_shells, number_of_cgto_shells
    logical*8, allocatable :: list(:)
    type( BTO_shell_data_obj ), allocatable :: dummy_bto(:)
    real(kind(1d0))               :: DeltaR
    real(kind(1d0)), allocatable  :: cf(:,:) 
    character(len=132)            :: name
    integer*8                     :: no_sym, nob(8), nlmq
    character(len=3), allocatable :: GChar(:)
    character(len=12)             :: IrrName,GrpName
    !.. Given iIrr =1, nIrreps ordering of irreps in ASTRA, irrO(iIrr)
    !.. returns the place of iIrr considered in the UKRmol modules,
    !.. while jrrO is the inverse of irrO.
    !.. arrO and brrO have the dimension of the total number of states
    !.. in Ukrmol and are used to point the ASTRA absolute indexes to the UKRmol ones.
    integer, allocatable :: irrO(:) , jrrO(:) , arrO(:) , brrO(:)
    integer, allocatable :: nVEC1(:), nVEC2(:), nVEC3(:), iLOC0(:), iLOC(:), iSRC(:), nTOT0(:)
    
    call self%free()

    call Init_UKRMOLP_BASIS( IntegralFile )

    self%SymmetryIdentifier = molecular_orbital_basis%pg
    self%nIrreps = molecular_orbital_basis%no_irr
    nIrreps = self%nIrreps
    allocate(self%nTOT(self%nIrreps))
    allocate(self%nLOC(self%nIrreps))
    allocate(self%nSRC(self%nIrreps))
    allocate(self%iLOC(self%nIrreps))
    allocate(self%iSRC(self%nIrreps))
    allocate(self%sTOT(0:self%nIrreps))
    self%nLOC=0
    self%nSRC=0
    self%sTOT=0

    do in8 = 1, self%nIrreps
       self%nTOT( int( in8 ) ) = molecular_orbital_basis%get_number_of_orbitals(in8)
       allocate( list( self%nTOT( int( in8 ) ) ) )
       call molecular_orbital_basis%get_continuum_flags(in8,list)
       self%nSRC( int( in8 ) ) = count( list )
       deallocate( list )
       self%nLOC( int( in8 ) ) = self%nTOT( int( in8 )    ) - self%nSRC( int( in8 ) ) 
       self%sTOT( int( in8 ) ) = self%sTOT( int( in8 ) -1 ) + self%nTOT( int( in8 ) )
    end do

    write(*,"(a,i0)"     ) "Symmetry Identifier:",  self%Get_Sym()
    write(*,"(a,*(x,i3))") "nIrrep             :",  self%nIrreps
    write(*,"(a,*(x,i3))") "nLOC               :", (self%nLOC(iIrrep),iIrrep=1,self%nIrreps)
    write(*,"(a,*(x,i3))") "nSRC               :", (self%nSRC(iIrrep),iIrrep=1,self%nIrreps)
    write(*,"(a,*(x,i3))") "nTOT               :", (self%nTOT(iIrrep),iIrrep=1,self%nIrreps)
    write(*,"(a,*(x,i3))") "sTOT               :", (self%sTOT(iIrrep),iIrrep=1,self%nIrreps)

    write(*,*) "Number of functions          ",molecular_orbital_basis%ao_basis%number_of_functions
    write(*,*) "Number of continuum functions",molecular_orbital_basis%ao_basis%n_cont_fns
    write(*,*) "Number of molecular functions",molecular_orbital_basis%ao_basis%n_cgto_functions
    write(*,*) "Numberof shells              ",molecular_orbital_basis%ao_basis%number_of_shells

    allocate( iLOC0(nIrreps) , nTOT0(nIrreps) )
    nTOT0        = self%nTOT
    self%iseries = self%nTOT
    iLOC0(1)=0
    do iIrrep = 1, nIrreps
       if(iIrrep>1) iLOC0(iIrrep) = &
            iLOC0(iIrrep - 1) + self%nLOC(iIrrep - 1) + self%nSRC(iIrrep - 1)
    enddo
    !write(*,*) "iLOC0",iLOC0
    
    ! do i=1,molecular_orbital_basis%ao_basis%number_of_shells
    !    write(*,'(7(3x,I3))') i,molecular_orbital_basis%ao_basis%shell_descriptor(1,i),&
    !         & molecular_orbital_basis%ao_basis%shell_descriptor(2,i),&
    !         & molecular_orbital_basis%ao_basis%shell_descriptor(3,i),&
    !         & molecular_orbital_basis%ao_basis%shell_descriptor(4,i),&
    !         & molecular_orbital_basis%ao_basis%shell_descriptor(5,i),&
    !         & molecular_orbital_basis%ao_basis%shell_descriptor(6,i)
    ! end do
    ! stop

    !> For each shell (2nd dimension) the first dimension contains the following info: shell type (CGTO: 1, BTO: 2)
    !row 1, index of the shell within its own type, row 2,
    !> is a continuum shell (1) or not (0), row 3, starting index for the functions in the shell, row 4,
    !the number of functions in the shell, row 5 and the angular momentum of the shell, row 6.
    !> The order in which the mapping is stored corresponds to the order in which the shells were added to the basis.
    !       integer, allocatable :: shell_descriptor(:,:)


    call molecular_orbital_basis%ao_basis%get_all_BTO_shells(dummy_bto,number_of_bto_shells)
    write(*,*) size(dummy_bto)
    !write(*,*) (dummy_bto(i)%bspline_index,i=1,number_of_bto_shells)
    write(*,"(a,f14.6,a,f14.6,a)") &
         "Radial Interval = [ ",&
         dummy_bto(1)%bspline_grid%A,", ", &
         dummy_bto(1)%bspline_grid%B,"]" 
    write(*,*) dummy_bto(1)%bspline_grid%order
    !
    !************ GET THE BSPLINE PARAMETERS *****************
    write(*,*) dummy_bto(1)%bspline_grid%n     !.. Number of Bsplines
!!$    write(*,*) dummy_bto(1)%bspline_grid%knots 


    number_of_shells      = molecular_orbital_basis%ao_basis%number_of_shells
    number_of_cgto_shells = number_of_shells - number_of_bto_shells

    !!vvvvvvv Juan vvvvvvv
    !NprimUkrmol = sum( molecular_orbital_basis%ao_basis%shell_descriptor( M_MULTIPLICITY_COLUMN, 1:number_of_shells ) )
    NprimUkrmol = molecular_orbital_basis%ao_basis%number_of_functions
    self%nf_bto = NprimUkrmol
    self%NprimUkrmol = NprimUkrmol
    self%ni_bto = molecular_orbital_basis%ao_basis%n_cgto_functions + 1
    allocate(self%l_bto(self%ni_bto:self%nf_bto))
    allocate(self%m_bto(self%ni_bto:self%nf_bto))
    allocate(self%nr_bto(self%ni_bto:self%nf_bto))
    do i = 1, number_of_shells
       if(molecular_orbital_basis%ao_basis%shell_descriptor(3,i).eq.1)then
          j = molecular_orbital_basis%ao_basis%shell_descriptor(4,i)
          l = molecular_orbital_basis%ao_basis%shell_descriptor(6,i)
          k = i - number_of_cgto_shells
          do m = -l, l
             self%l_bto( j )  = l
             self%m_bto( j )  = m
             self%nr_bto( j ) = dummy_bto(k)%bspline_index
             j = j + 1
          end do
       end if
    end do
    self%nlr_bs = self%nr_bto(self%nf_bto)  !Index of the last radial bspline
    !!--------------------

    self%lmax = molecular_orbital_basis%ao_basis%shell_descriptor( L_COLUMN, number_of_shells )

    !.. UKRMOL+ Bspline parameters
    self%order         = dummy_bto(1)%bspline_grid%order
    self%BsA_UKRmol    = dummy_bto(1)%bspline_grid%A
    self%BsB_UKRmol    = dummy_bto(1)%bspline_grid%B
    self%nBs_UKRmol    = dummy_bto(1)%bspline_grid%n
    self%nnodes_UKRmol = self%nBs_UKRmol - self%order + 2

    !.. ASTRA Bspline parameters
    DeltaR = ( self%BsB_UKRmol-self%BsA_UKRmol ) / dble(self%nnodes_UKRmol - 1) 
    self%Bs_Rmin_Astra = self%BsB_UKRmol - DeltaR * 2.d0 * ( self%order - 1)
    self%Bs_Rmax_Astra = self%BsB_UKRmol + DeltaR * N_EXT_BSPLINE_INTERVALS
    self%nNodesAstra   = 2 * self%order - 1 + N_EXT_BSPLINE_INTERVALS
    self%iBs_FirstExt  = 2 * self%order - 1


    !!vvvvvvv Juan vvvvvvv
    self%nBsplines    = self%nNodesAstra + self%order - 2
    self%iBs_LastExt  = self%nBsplines - ( self%order - 1)
    self%nExtBs_Astra = self%iBs_LastExt - self%iBs_FirstExt + 1 ! Minus the last one
    write(*,*) self%nBsplines
    write(*,*) self%iBs_LastExt
    write(*,*) self%iBs_FirstExt
    write(*,*) self%nExtBs_Astra
    write(*,*) "astra bspline parameters"

! write(*,*) "Bs_Rmin_Astra",self%Bs_Rmin_Astra
! write(*,*) "Bs_Rmax_Astra",self%Bs_Rmax_Astra
! write(*,*) "nNodesAstra",self%nNodesAstra
! write(*,*) "Order:",self%order
! write(*,*) "self%nBsplines",self%nBsplines
! write(*,*) "iBs_LastExt",self%iBs_LastExt
! write(*,*) "iBs_FirstExt",self%iBs_FirstExt
! write(*,*) "self%nExtBs_Astra",self%nExtBs_Astra
! write(*,*) "DeltaR",DeltaR
! write(*,*) "N_EXT_BSPLINE_INTERVALS",N_EXT_BSPLINE_INTERVALS
! pause

    !!-------------------
    allocate(self%cf(1:NprimUkrmol,1:sum(self%nTOT)))
    self%cf=0.d0

   
    call molecular_orbital_basis%get_orbital_coefficient_matrix(self%cf)
    call READ_UKRMOLP_PROPERTY_INTS(UKRMOL_NFTI,3_8,IntegralFile)

    ! this gives the info of the total number "nlmq" of property integrals in moints ukrmol file.
    call GET_NAME_SYM(name,no_sym,nob,nlmq)
    self%nlmq    = nlmq
    self%Lmax_mp = int(sqrt(dble(nlmq))) - 1


    !.. From now on we set everithing to change the Irrep order used in Ukrmol to the one set through the ASTRA modules.
    !   The irreps ordering in Ukrmol can be accesed through the Cons

    
    allocate( GChar( nIrreps ) )  
    GrpNAme = GlobalGroup%GetName()
    if(     GrpName .eq. "D2h" )then
       GChar = D2h_names
    elseif( GrpName .eq. "D2"  )then
       GChar = D2_names
    elseif( GrpName .eq. "C2h" )then
       GChar = C2h_names
    elseif( GrpName .eq. "C2v" )then
       GChar = C2v_names
    elseif( GrpName .eq. "C2"  )then
       GChar = C2_names
    elseif( GrpName .eq. "Cs"  )then
       GChar = Cs_names
    elseif( GrpName .eq. "Ci"  )then
       GChar = Ci_names
    elseif( GrpName .eq. "C1"  )then
       GChar = C1_names
    endif

    !.. irr0(iIrr) is the ukrmol irrep index jIrr as a function of the astra  irrep index iIrr
    !.. jrr0(jIrr) is the astra  irrep index iIrr as a function of the ukrmol irrep index jIrr
    allocate( irrO(nIrreps) , jrrO(nIrreps) )
    do iIrr = 1, nIrreps
       IrrName = GlobalGroup%getIrrName(iIrr)
       do jIrr = 1, nIrreps
          if( IrrName .eq. GChar(jIrr) )then
             irrO(iIrr) = jIrr
             jrrO(jIrr) = iIrr
             cycle
          end if
       enddo
    enddo
    ! do i = 1, nIrreps
    !    write(*,*) i,irro(i),jrro(i)
    ! end do
    ! pause
    
    allocate( arrO(sum(self%nTOT)) , brrO(sum(self%nTOT)) )
    do i = 1, sum(self%nTOT)
       arrO(i) = i
    end do

    allocate(nVEC1(nIrreps),nVEC2(nIrreps),nVEC3(nIrreps))
    nVEC1 = self%nTOT
    nVEC2 = self%nLOC
    nVEC3 = self%nSRC
    do iIrrep = 1, nIrreps
       self%nTOT(iIrrep) = nVEC1(irrO(iIrrep))
       self%nLOC(iIrrep) = nVEC2(irrO(iIrrep))
       self%nSRC(iIrrep) = nVEC3(irrO(iIrrep))
    end do
    deallocate(nVEC1,nVEC2,nVEC3)
    self%sTOT = 0
    do iIrrep = 1, nIrreps
       self%sTOT( iIrrep ) = self%sTOT( iIrrep -1 ) + self%nTOT( iIrrep )
    end do

    allocate( iLOC(nIrreps), iSRC(nIrreps) )

    self%iLOC(1) = 0
    do iIrrep = 1, nIrreps
       if(iIrrep>1) self%iLOC(iIrrep) = &
            self%iLOC(iIrrep - 1) + &
            self%nLOC(iIrrep - 1) + &
            self%nSRC(iIrrep - 1)
       self%iSRC(iIrrep) = self%iLOC(iIrrep) + self%nLOC(iIrrep    )
    enddo
    !write(*,*) "iLOC",iLOC
    do iIrrep = 1, nIrreps
       brrO(self%iLOC ( iIrrep  ) + 1 : &
            self%iLOC ( iIrrep  ) + self%nTOT( iIrrep ) ) = &
            arrO( &
            iLOC0( irrO(iIrrep) ) + 1 : &
            iLOC0( irrO(iIrrep) ) + nTOT0( irrO(iIrrep) ) )
    enddo
    deallocate(iLOC0,nTOT0)
    allocate(self%indexASTRA(sum(self%nTOT)))
    self%indexASTRA = brrO


    allocate(cf(1:NprimUkrmol,1:sum(self%nTOT)))
    do i = 1, sum(self%nTOT)
       cf(:,i) = self%cf(:,brrO(i))
    end do
    self%cf = cf
    deallocate(cf)
    deallocate(brrO)
  end subroutine Class_UKRmol_orbitals_InitFromData


  integer function Class_UKRmol_orbitals_Get_BsFirst( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%iBs_FirstExt
  end function Class_UKRmol_orbitals_Get_BsFirst

  integer function Class_UKRmol_orbitals_Get_BsLast( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%iBs_LastExt
  end function Class_UKRmol_orbitals_Get_BsLast
  
  integer function Class_UKRmol_orbitals_Get_nExBsAstra( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nExtBs_Astra
  end function Class_UKRmol_orbitals_Get_nExBsAstra
  
  integer function Class_UKRmol_orbitals_Get_nBs( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nBsplines 
  end function Class_UKRmol_orbitals_Get_nBs
  
  integer function Class_UKRmol_orbitals_Get_Bsorder( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%order
  end function Class_UKRmol_orbitals_Get_Bsorder

  integer function Class_UKRmol_orbitals_Get_BsA_UKRmol( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%BsA_UKRmol
  end function Class_UKRmol_orbitals_Get_BsA_UKRmol

  integer function Class_UKRmol_orbitals_Get_IntIndex( self, i ) result( j )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer, intent(in)  :: i
    j = self%indexASTRA(i)
  end function Class_UKRmol_orbitals_Get_IntIndex
  
  integer function Class_UKRmol_orbitals_Get_BsB_UKRmol( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%BsB_UKRmol
  end function Class_UKRmol_orbitals_Get_BsB_UKRmol

  integer function Class_UKRmol_orbitals_Get_nBs_UKRmol( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nBs_UKRmol
  end function Class_UKRmol_orbitals_Get_nBs_UKRmol
  
  integer function Class_UKRmol_orbitals_Get_BsRmaxAstra( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%Bs_Rmax_Astra
  end function Class_UKRmol_orbitals_Get_BsRmaxAstra
    
  integer function Class_UKRmol_orbitals_Get_BsRminAstra( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%Bs_Rmin_Astra
  end function Class_UKRmol_orbitals_Get_BsRminAstra

  integer function Class_UKRmol_orbitals_Get_nNodesAstra( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nNodesAstra
  end function Class_UKRmol_orbitals_Get_nNodesAstra
  
  integer function Class_UKRmol_orbitals_Get_order( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%order
  end function Class_UKRmol_orbitals_Get_order

  integer function Class_UKRmol_orbitals_Get_nEBs( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nExtBs_Astra
  end function Class_UKRmol_orbitals_Get_nEBs

  integer function Class_UKRmol_orbitals_Get_lmax( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%lmax
  end function Class_UKRmol_orbitals_Get_lmax

  !.. This gives the Lmas_mp corresponding to the maximum number of property integrals evaluated by Ukrmol.
  integer function Class_UKRmol_orbitals_Get_Lmax_mp( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%Lmax_mp
  end function Class_UKRmol_orbitals_Get_Lmax_mp

  integer function Class_UKRmol_orbitals_Get_l_bto( self, i ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: i
    res = self%l_bto( i )
  end function Class_UKRmol_orbitals_Get_l_bto

  integer function Class_UKRmol_orbitals_Get_nr_bto( self, i ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: i
    res = self%nr_bto( i )
  end function Class_UKRmol_orbitals_Get_nr_bto

  integer function Class_UKRmol_orbitals_Get_ni_bto( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%ni_bto
  end function Class_UKRmol_orbitals_Get_ni_bto
  
  integer function Class_UKRmol_orbitals_Get_nf_bto( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nf_bto
  end function Class_UKRmol_orbitals_Get_nf_bto
  
  integer function Class_UKRmol_orbitals_Get_nlr_bs( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nlr_bs
  end function Class_UKRmol_orbitals_Get_nlr_bs

  real(kind(1d0)) function Class_UKRmol_orbitals_Get_cf( self , ir, ic ) result( cf)
    class( Class_UKRmol_orbitals ), intent(in)  :: self
    integer                       , intent(in)  :: ir, ic
    cf = self%cf(ir,ic)
  end function Class_UKRmol_orbitals_Get_cf

  !*** DOES NOT WORK, CANT BE SEEN FROM MAIN
  ! subroutine Class_UKRmol_orbitals_Get_cf( self , cf, nr, nc )
  !   class( Class_UKRmol_orbitals ), intent(in)  :: self
  !   integer                       , intent(in)  :: nr, nc
  !   real(kind(1d0))               , intent(out) :: cf( nr, nc )
  !   cf(:,:) = self%cf(:,:)
  ! end subroutine Class_UKRmol_orbitals_Get_cf

  integer function Class_UKRmol_orbitals_Get_NprimUKRmol( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%NprimUKRmol
  end function Class_UKRmol_orbitals_Get_NprimUKRmol

  integer function Class_UKRmol_orbitals_Get_Sym( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%SymmetryIdentifier
  end function Class_UKRmol_orbitals_Get_Sym

  integer function Class_UKRmol_orbitals_Get_nIrr( self ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    res = self%nIrreps
  end function Class_UKRmol_orbitals_Get_nIrr

  integer function Class_UKRmol_orbitals_Get_nLOC( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%nLOC(iIrr)
  end function Class_UKRmol_orbitals_Get_nLOC

  integer function Class_UKRmol_orbitals_Get_nSRC( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%nSRC(iIrr)
  end function Class_UKRmol_orbitals_Get_nSRC

  integer function Class_UKRmol_orbitals_Get_iLOC( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%iLOC(iIrr)
  end function Class_UKRmol_orbitals_Get_iLOC

  integer function Class_UKRmol_orbitals_Get_iSRC( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%iSRC(iIrr)
  end function Class_UKRmol_orbitals_Get_iSRC

  integer function Class_UKRmol_orbitals_Get_nTOT( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%nTOT(iIrr)
  end function Class_UKRmol_orbitals_Get_nTOT

  !.. This is the same as nTOT but with the original ordering, needed to initialize the ukrmol integrals.
  integer function Class_UKRmol_orbitals_Get_iseries( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%iseries(iIrr)
  end function Class_UKRmol_orbitals_Get_iseries
  
  integer function Class_UKRmol_orbitals_Get_sTOT( self, iIrr ) result( res )
    class( Class_UKRmol_orbitals ), intent(in) :: self
    integer                       , intent(in) :: iIrr
    res = self%sTOT(iIrr)
  end function Class_UKRmol_orbitals_Get_sTOT

  !> Irr_lm returns the irrep corresponding to a given pair of lm,
  !! according to the correspondence in ukrmol (same order than Dalton):
  !!
  !! *** MUST EXTEND TO OTHER GROUPS
  !!
  !!      Ag    Xlm       l even m even  N =  l/2+1      1  ->  1
  !!      B3u   Xlm       l odd  m odd   N = (l+1)/2     2  ->  8
  !!      B2u   Xl-m      l odd  m odd   N = (l+1)/2     3  ->  7
  !!      B1g   Xl-m      l even m even  N =  l/2        4  ->  2 
  !!      B1u   Xlm       l odd  m even  N = (l+1)/2     5  ->  6
  !!      B2g   Xlm       l even m odd   N =  l/2        6  ->  3
  !!      B3g   Xl-m      l even m odd   N =  l/2        7  ->  4
  !!      Au    Xl-m      l odd  m even  N = (l-1)/2     8  ->  5
  !! 
  !! which differs from the one used in CloseCouplingBasis (XCHEM)
  !!      Ag    Xlm       l even m even  N =  l/2+1
  !!      B1g   Xl-m      l even m even  N =  l/2
  !!      B2g   Xlm       l even m odd   N =  l/2
  !!      B3g   Xl-m      l even m odd   N =  l/2
  !!      Au    Xl-m      l odd  m even  N = (l-1)/2
  !!      B1u   Xlm       l odd  m even  N = (l+1)/2
  !!      B2u   Xl-m      l odd  m odd   N = (l+1)/2
  !!      B3u   Xlm       l odd  m odd   N = (l+1)/2
  !<
  function Irr_lm( l, m) result( i )
    implicit none
    integer, intent(in) :: l, m 
    integer             :: i    
    i = 0
    if ( ( mod( l, 2 ) .eq. 0 ) .and. ( mod( m, 2 ) .eq. 0 ) ) then
       if ( m .ge. 0 ) then
          i = 1 
       else
          i = 4 
       end if
    elseif( ( mod( l, 2 ) .eq. 0 ) .and. ( abs( mod( m, 2 ) ) .eq. 1 ) ) then
       if ( m .ge. 0 ) then
          i = 6 
       else
          i = 7 
       end if
    elseif ( ( mod( l, 2 ) .eq. 1 ) .and. ( mod( m, 2) .eq. 0 ) ) then
       if ( m .ge. 0 ) then
          i = 5 
       else
          i = 8  
       end if
    elseif ( ( mod( l, 2 ) .eq. 1 ) .and. ( abs( mod( m, 2 ) ) .eq. 1 ) ) then
       if ( m .ge. 0 ) then
          i = 2 
       else
          i = 3 
       end if
    end if
  end function Irr_lm


  ! !This is the integral of three spherical harmonics (none of them conjugate)
  function ThreeYlmIntegral( l1, m1, l2, m2, l3, m3 ) result(inte)

    use ModuleAngularMomentum
    implicit none

    integer, intent(in)        :: l1, m1, l2, m2, l3, m3 ! input
    real(kind(1d0))            :: inte ! output
    real(kind(1d0)), parameter :: PI = 3.14159265358979323844d0

    inte = 0.d0
    if ( ( m1 + m2 ) .eq. ( -m3 ) ) then
       inte = ( (-1.d0)**( m3 ) ) * sqrt( ( (2*l1+1)*(2*l2+1) ) / ( (4*PI)*(2*l3+1) ) )
       inte = inte * &
            ClebschGordanCoefficient(l1, l2, l3, m1, m2) * &
            ClebschGordanCoefficient(l1, l2, l3, 0, 0)
    end iF
  end function ThreeYlmIntegral

  !This is the integral of Xl1m1  Xl2m2 Yl3m3 
  ! function TwoXlmOneYlmIntegral( l1, m1, l2, m2, l3, m3 ) result( inte )

  !   implicit none

  !   integer, intent(in)         :: l1,m1,l2,m2,l3,m3 ! input
  !   real(kind(1d0))             :: inte,amp1,fas1,amp2,fas2 ! output
  !   real(kind(1d0)), parameter  :: PI=3.14159265358979323844d0

  !   inte = 0.d0
  !   if( m1 .eq. 0 )then
  !      amp1 = 1.d0
  !      fas1 = 0.d0
  !   elseif( m1 .lt. 0 )then
  !      amp1 = 1.d0 / sqrt(2.d0)
  !      fas1 = amp1 * ((-1.d0)**abs(m1))
  !   else
  !      amp1 =  1.d0/sqrt(2.d0)
  !      fas1 = -amp1*((-1.d0)**abs(m1))
  !   end if

  !   if( m2 .eq. 0 )then
  !      amp2 = 1.d0
  !      fas2 = 0.d0
  !   elseif( m2 .lt. 0 )then
  !      amp2 = 1.d0/sqrt(2.d0)
  !      fas2 = amp2*((-1.d0)**abs(m2))
  !   else
  !      amp2 =  1.d0/sqrt(2.d0)
  !      fas2 = -amp2*((-1.d0)**abs(m2))
  !   end if
  !   inte = 0.d0
  !   inte = amp1 * amp2 *        ThreeYlmIntegral( l1,  abs(m1), l2,  abs(m2), l3, m3)
  !   inte = inte + amp1 * fas2 * ThreeYlmIntegral( l1,  abs(m1), l2, -abs(m2), l3, m3)
  !   inte = inte + fas1 * amp2 * ThreeYlmIntegral( l1, -abs(m1), l2,  abs(m2), l3, m3)
  !   inte = inte + fas1 * fas2 * ThreeYlmIntegral( l1, -abs(m1), l2, -abs(m2), l3, m3)

  ! end function TwoXlmOneYlmIntegral

  !This is the integral of Xl1m1  Xl2m2 Yl3m3 
    function ThreeXlmIntegral( l1, m1, l2, m2, l3, m3) result( inteR )

      implicit none

      integer, intent(in)    :: l1, m1, l2, m2, l3, m3
      real(kind(1d0))        :: inteR
      complex(kind(1d0))     :: inte, amp1, fas1, amp2, fas2, amp3, fas3
      complex(kind(1d0)), parameter  :: ci = (0.d0,1.d0), c1 = (1.d0,0.d0), c0 = (0.d0,0.d0)
      real   (kind(1d0)), parameter  :: PI = 3.14159265358979323844d0

      inte=c0
      if ( m1 .eq. 0 )then
         amp1 = c1*1.d0
         fas1 = c1*0.d0
      elseif( m1 .lt. 0 )then
         amp1 = -ci*1.d0/sqrt(2.d0)
         fas1 = -amp1*((-1.d0)**abs(m1))
      else
         amp1 = c1*1.d0/sqrt(2.d0)
         fas1 = amp1*((-1.d0)**abs(m1))
      end if
      if( m2 .eq. 0 )then
         amp2 = c1*1.d0
         fas2 = c1*0.d0
      elseif( m2 .lt. 0 )then
         amp2 = -ci*1.d0/sqrt(2.d0)
         fas2 = -amp2*((-1.d0)**abs(m2))
      else
         amp2 = c1*1.d0/sqrt(2.d0)
         fas2 = amp2*((-1.d0)**abs(m2))
      end if
      if( m3 .eq. 0 )then
         amp3 = c1*1.d0
         fas3 = c1*0.d0
      elseif( m3 .lt. 0 )then
         amp3 = -ci*1.d0/sqrt(2.d0)
         fas3 = -amp3*((-1.d0)**abs(m3))
      else
         amp3 = c1*1.d0/sqrt(2.d0)
         fas3 = amp3*((-1.d0)**abs(m3))
      end if
      inte = c0
     ! m1=m2=m3=0
      inte =        amp1 * amp2 * amp3 * ThreeYlmIntegral( l1,  abs(m1), l2,  abs(m2), l3,  abs(m3) )
      inte = inte + fas1 * fas2 * fas3 * ThreeYlmIntegral( l1, -abs(m1), l2, -abs(m2), l3, -abs(m3) ) !this is exactly 0 because fas1,2,3=0

  !    If m1=0 this is multiplied by i^{-1} because either m2 or m3 are negative
      inte = inte + amp1 * amp2 * fas3 * ThreeYlmIntegral( l1,  abs(m1), l2,  abs(m2), l3, -abs(m3) ) 
      inte = inte + amp1 * fas2 * amp3 * ThreeYlmIntegral( l1,  abs(m1), l2, -abs(m2), l3,  abs(m3) )

      inte = inte + amp1 * fas2 * fas3 * ThreeYlmIntegral( l1,  abs(m1), l2, -abs(m2), l3, -abs(m3) )

      inte = inte + fas1 * amp2 * amp3 * ThreeYlmIntegral( l1, -abs(m1), l2,  abs(m2), l3,  abs(m3) ) 
      inte = inte + fas1 * amp2 * fas3 * ThreeYlmIntegral( l1, -abs(m1), l2,  abs(m2), l3, -abs(m3) )

      inte = inte + fas1 * fas2 * amp3 * ThreeYlmIntegral( l1, -abs(m1), l2, -abs(m2), l3,  abs(m3) )

      inteR = real(inte)
      if( abs( aimag(inte) ) .gt. ( 10.d0**(-10)) )then
         write(*,*) "matrix element is complex, something is wrong",inte
         write(*,*) l1,m1,l2,m2,l3,m3
         stop
      end if

    end function ThreeXlmIntegral

  ! Evaluates  the real spherical harmonics
  ! function Xlm( theta, phi, l, m ) result( inteR )

  !   use ModuleAngularMomentum
  !   implicit none

  !   real(kind(1d0)), intent(in) :: theta, phi
  !   integer, intent(in)         :: l, m 
  !   real(kind(1d0))             :: inteR
  !   complex(kind(1d0))          :: inte, amp1, fas1, amp2, fas2, amp3, fas3 
  !   complex(kind(1d0)), parameter  :: ci = (0.d0,1.d0), c1 = (1.d0,0.d0), c0 = (0.d0,0.d0)
  !   real   (kind(1d0)), parameter  :: PI = 3.14159265358979323844d0

  !   inte = c0
  !   if( m .eq. 0 )then
  !      amp1 = c1*1.d0
  !      fas1 = c1*0.d0
  !   elseif(m .lt. 0 )then
  !      amp1 = -ci*1.d0/sqrt(2.d0)
  !      fas1 =  amp1*((-1.d0)**abs(m)) 
  !   else
  !      amp1 =  c1*1.d0/sqrt(2.d0)
  !      fas1 = -amp1*((-1.d0)**abs(m))
  !   end if

  !   inte = c0
  !  ! m1=m2=m3=0
  !   inte = amp1 * Ylm( theta, phi, l, abs(m) ) + fas1 * Ylm( theta, phi, l, -abs(m) )

  !   inteR = real( inte )
  !   if( abs( aimag( inte ) ) .gt. 0.d0 )then
  !      write(*,*) "Xlm is complex, something is wrong"
  !      write(*,*) theta,phi,l,m
  !      stop
  !   end if

  ! end function Xlm

end Module ModuleUKRmolInterface
