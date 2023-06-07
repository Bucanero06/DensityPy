
! {{{ Detailed description

!> \file
!!
!! Contains all the needed attributes and procedures to manage an arbitrary radial basis.

! }}}
module ModuleBasis

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use ModuleSystemUtils
  use ModuleErrorHandling
  use ModuleParameterList
  use ModuleBSpline
  use ModuleString
  use ModuleGaussian
  use ModuleGaussianBSpline
  use ModuleMatrix

  implicit none
  private

  !> Enumerates the radial basis' possible kinds.
  enum, bind(c)
     enumerator :: BASIS_KIND_BSPLINE
     enumerator :: BASIS_KIND_GAUSSIANBSPLINE
     enumerator :: BASIS_KIND_GAUSSIAN
     enumerator :: BASIS_KIND_UNDEFINED
  end enum


  !> Specifier of the subset of basis on which a 
  !!   vector/matrix dimension is expressed.
  !..
  enum, bind(c)
     enumerator :: REP_LEVEL_TOTAL_INIT
     enumerator :: REP_LEVEL_REGULAR_INIT
     enumerator :: REP_LEVEL_PRECONDITIONED_INIT
     enumerator :: REP_LEVEL_PRECONDITIONED_CONTINUOUS_INIT
     enumerator :: REP_LEVEL_ORTHONORMAL_INIT
  end enum

  !> TOTAL:
  !! Complete set of radial functions
  integer, public, parameter :: &
       REP_LEVEL_TOTAL = &
       REP_LEVEL_TOTAL_INIT

  !> REGULAR:
  !! Subset of TOTAL compatible with a given pair of angular momenta.
  !! (i.e.: \f$ u_l(r) ~ r^{l+1} \f$ close to the origin)
  integer, public, parameter :: &
       REP_LEVEL_REGULAR = &
       REP_LEVEL_REGULAR_INIT

  !> PRECONDITIONED:
  !! Linearly independent system of functions build from the REGULAR
  !! set according to criteria specified by the Basis module.
  !! In the B-spline case, it is required that the last B-spline is 
  !! assumed to be linearly independent of all other functions and 
  !! not mixed with them.
  integer, public, parameter :: &
       REP_LEVEL_PRECONDITIONED = &
       REP_LEVEL_PRECONDITIONED_INIT

  !> PRECONDITIONED_CONTINUOUS:
  !! Continuous subset of PRECONDITIONED. In practice, if the basis 
  !! comprises a set of B-splines, PRECONDITIONED_CONTINUOUS is 
  !! obtaiend from PRECONDITIONED by excluding the last B-spline.
  integer, public, parameter :: &
       REP_LEVEL_PRECONDITIONED_CONTINUOUS = &
       REP_LEVEL_PRECONDITIONED_CONTINUOUS_INIT

  !> ORTHONORMAL:
  !! Orthonormal set build from the set of PRECONDITIONED_CONTINUOUS functions
  integer, public, parameter :: &
       REP_LEVEL_ORTHONORMAL = &
       REP_LEVEL_ORTHONORMAL_INIT
  

  character(len=*), public, parameter :: BANDED_PATTERN_IDENTIFIER = "BANDED"
  character(len=*), public, parameter :: FULL_PATTERN_IDENTIFIER   = "FULL"

  character(len=*), parameter :: BSPLINE_STRN_IDENTIFIER = "BSpline"
  character(len=*), parameter :: BSPLINE_BEST_PATTERN    =  BANDED_PATTERN_IDENTIFIER

  character(len=*), parameter :: GAUSSIAN_STRN_IDENTIFIER= "Gaussian"
  character(len=*), parameter :: GAUSSIAN_BEST_PATTERN   =  FULL_PATTERN_IDENTIFIER

  character(len=*), parameter :: GAUSSIANBSPLINE_STRN_IDENTIFIER     = "GaussianBSpline"
  character(len=*), parameter :: GAUSSIANBSPLINE_BEST_PATTERN        =  FULL_PATTERN_IDENTIFIER

  integer,          parameter :: DEFAULT_OUTPUT_UNIT     =  OUTPUT_UNIT




  type, public :: ClassBasis

     private

     !> Directory in which all the radial basis information will be stored.
     character(len=:), allocatable :: Dir

     !> Pointer that makes the reference to the proper radial basis kind 
     !! as specified by the configuration files.
     class(*)        , pointer     :: RadialBasis => NULL()

     !> The number of points that will be used to plot the radial basis functions.
     integer                       :: NPlot = 0

     !> Indicates whether the radial basis functions will be plotted or not.
     logical                       :: PlotRequired = .FALSE.

     !> Class with the aim to manage the configuration 
     !! file reading ([ModuleParameterList](@ref ModuleParameterList.f90)).
     type(ClassParameterList)      :: List
     !> Distance at which the plot will start.
     real(kind(1d0)) :: RPlotMin
     !> Distance at which the plot will finish.
     real(kind(1d0)) :: RPlotMax

   contains

     !> Frees the ClassBasis attributes.
     procedure :: Free => BasisFree

     !> Initializes ClassBasis reading the information from a configuration file.
     procedure :: Init => BasisInitReadConfigurationFile

     !> Checks whether ClassBasis has been initialized or not.
     procedure :: Initialized => BasisIsInitialized

     !> Gets the number of radial basis functions used in the computation of the operators matrix elements.
     procedure :: GetNFun

     !> Gets the number of B-splines to be skipped.
     procedure :: GetNBSplineSkip

     !> Returns the basis function normalization factor.
     procedure :: GetNormFactor

     !> Retrieves the directory in wich the radial basis information will be stored.
     procedure :: GetDir

     !> Retrieves the number of points that will be used to plot the radial basis functions.
     procedure :: GetNPlot

     !> Returns a linear grid on the radial domain appropriate for the current radial basis.
     procedure :: GetPlotGrid

     !> Gets the maximum angular momentum in the Gaussian expansion.
     procedure :: GetGaussianLMax

     !> Gets the number of exponents in the Gaussian expansion.
     procedure :: GetGaussianNumExponents

     !> Gets the number of exponents in the Gaussian expansion.
     procedure :: RemoveGaussianRowNormalization

     !> Gets a vector with all the monomials exponents of the Gaussian functions.
     procedure, public :: GetAllMonomials

     !> Gets the minimum and maximum index corresponding to a given angular momentum
     procedure, public :: GetLRowInterval

     !> Gets a vector with all the exponents in the exponential of the Gaussian functions.
     procedure, public :: GetAllExponents

     !> Get maximum radius for the current basis.
     !! Coincides with the maximum plot radius.
     procedure :: GetRMax
     
     ! Get a given node position for B-splines.
     procedure :: GetNodePosition

     !> Gets the radial distance at which the plot will start.
     procedure :: GetRPlotMin

     !> Gets the radial distance at which the plot will finish.
     procedure :: GetRPlotMax

     !> Gets the information of the linear independent gaussian function (preconditioned): monomial exponent, gaussian exponent and normalization factor.
     procedure :: GetLinIndGaussian

     !> Retrieves whether the radial basis functions will be plotted or not.
     procedure :: ToBePlotted

     !> Performs the plot of the radial basis.
     procedure :: Plot => PlotBasis

     !> Computes the  k-th derivatives either of a choosen radial basis function, 
     !! or of a function expressed in terms of the radial basis, both evaluated 
     !! in certain position. 
     generic   :: Eval => BasisFunctionEval, BasisExpansionEval

     !> Computes the  k-th derivatives either of a choosen radial basis function, 
     !! or of a function expressed in terms of the radial basis, both evaluated 
     !! in a position array.
     generic   :: Tabulate =>  BasisFunctionTabulate, BasisExpansionTabulate

     !> Computes the integral
     !! \f[
     !! \int_{a}^{b}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr
     !!\f]
     !! Where \f$f(r)\f$ is a local operator and \f$G_{i,j}\f$ radial basis functions. 
     !! If a break point \f$BP\f$ is introduced, then the integral is splitted in two parts
     !! \f[
     !!   \int_{a}^{BP}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr+
     !!   \int_{BP}^{b}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr
     !! \f]
     procedure :: Integral 
     !
     !> Computes the integrals
     !! \f[
     !! \int_{a}^{b}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr
     !!\f]
     !! For \f$j\f$ from some value \f$j_{min}\f$ to \f$j_{max}\f$, and stores the results in a vector with \f$j_{max}-j_{min}+1\f$ components. Where \f$f(r)\f$ is a local operator and \f$G_{i,j}\f$ radial basis functions. 
     !! If a break point \f$BP\f$ is introduced, then the integral is splitted in two parts
     !! \f[
     !!   \int_{a}^{BP}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr+
     !!   \int_{BP}^{b}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr
     !! \f]
     procedure :: IntegralRow
     !
     !> Computes the integrals
     !! \f[
     !! \int_{a}^{b}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr
     !!\f]
     !! For \f$i\f$ from some value \f$i_{min}\f$ to \f$i_{max}\f$, and stores the results in a vector with \f$i_{max}-i_{min}+1\f$ components. Where \f$f(r)\f$ is a local operator and \f$G_{i,j}\f$ radial basis functions. 
     !! If a break point \f$BP\f$ is introduced, then the integral is splitted in two parts
     !! \f[
     !!   \int_{a}^{BP}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr+
     !!   \int_{BP}^{b}\frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)\frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr
     !! \f]
     procedure :: IntegralCol

     !> Computes the integral
     !! \f[
     !!    \int_a^b dr 
     !!    \frac{d^{n_1}G(r)}{dr^{n_1}} f(r)
     !! \f]
     !! where \f$f(r)\f$ is a *smooth* function in the radial basis support.
     !! If the first (last) boundary of the integration interval is not specified, 
     !! then the value \f$a=-\infty\f$ and \f$b=\infty\f$ is assumed.
     procedure :: MonoIntegral

     !> Computes the principal-part radial basis functions integral
     !! \f[
     !!  \int dr G_n(r)\,\frac{\mathcal{P}}{r-r_0}\,G_m(r)
     !! \f]
     procedure :: CauchyIntegral

     !> Computes the monoelectronic or bielectronic multipole.
     generic   :: Multipole  =>  MonoelectronicMultipole, BielectronicMultipole

     !> Computes the "Spherical Bessel Transform" of a radial basis function 
     !! \f$G(r)\f$, defined as:
     !! \f[
     !!   F_{l}(p)=\int\frac{G(r)}{r}j_{l}(p\cdot r)r^{2}dr,
     !! \f]
     !! for all the angular momenta \f$\ell\f$ up to a specified lmax as 
     !! \f$ res(\ell+1)=F_\ell(p)\f$.
     procedure :: SphericalBesselTransform

     !> Checks whether tha radial basis characteristics (fingerprint) saved previously 
     !! in a unit, coincides or not with the current radial basis in used.
     procedure :: MatchFingerPrintOnUnit

     !> Writes in a unit the radial basis kind (B-splines, Gaussian, etc.).
     procedure :: WriteKind => BasisWriteKind

     !> Writes in a unit the adequate radial basis information (fingerprint) 
     !! that allows to identify it from other basis that could possible be used.
     procedure :: WriteFingerprint 

     !> Depending on the radial basis kind, selects the best matrix pattern to 
     !! store the local operators matrix elements: Full or Banded.
     procedure :: BestPattern

     !> Determines the number of subdiagonals of the operator matrices computed 
     !! using the radial basis.
     procedure :: LowerBandwidth

     !> Determines the number of superdiagonals of the operator matrices computed using the radial basis.
     procedure :: UpperBandwidth
     !> Creates in a directory all the configuration files that define the current basis.
     procedure :: SaveConfig

     !> Gets the total number of radial basis functions.
     procedure :: GetTotalNumBasisFunctions


     !> Retrieves the set of regular indexes of the radial basis functions
     !! for a given angular momentum
     procedure :: GetRegularIndexes
     !> Retrieves the number of regular radial basis functions.
     procedure :: NRegular
     !> Retrieves the number of preconditioned radial basis functions.
     procedure :: NPreconditioned
     !> Retrieves the number of continuous preconditioned radial basis functions.
     procedure :: NPreconditionedContinuous
     !> Builds the directory in which the results will be stored.
     procedure :: SetRoot

     !> Retrieves whether the preconditioner is available for certain angular momentum or not.
     procedure :: PreconditionerIsAvailable
     !> Loads the preconditioner.
     procedure :: LoadPreconditioner
     !> Computes the preconditioner.
     procedure :: ComputePreconditioner
     !> Saves the preconditioner.
     procedure :: SavePreconditioner

     !> Casts a matrix which is in a certain representation to obtain the matrix belonging to other requested one.
     procedure :: CastMatrix

     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to other bras' representation.
     generic   :: CastMatrixBra => &
          CastMatrixBra_FromClassMatrix,&
          CastMatrixBra_FromMatrix, &
          CastMatrixBra_FromVector
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to other bras' representation. The original matrix is in [ClassMatrix](@ref ModuleMatrix.f90) format.
     procedure :: CastMatrixBra_FromClassMatrix
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to the total basis bras' representation. The original matrix is in [ClassMatrix](@ref ModuleMatrix.f90) format.
     procedure :: CastMatrixBraToTotal
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to the regular basis bras' representation. The original matrix is in [ClassMatrix](@ref ModuleMatrix.f90) format.
     procedure :: CastMatrixBraToRegular
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to the preconditioned basis bras' representation. The original matrix is in [ClassMatrix](@ref ModuleMatrix.f90) format.
     procedure :: CastMatrixBraToPreconditioned
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to the continuous preconditioned basis bras' representation. The original matrix is in [ClassMatrix](@ref ModuleMatrix.f90) format.
     procedure :: CastMatrixBraToPreconditionedContinuous
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M' x N) matrix belonging to other bras' representation. The original matrix is stored as a 2-dimensional array.
     procedure :: CastMatrixBra_FromMatrix
     !
     !> Casts an original (M x 1) vector in a certain representation to obtain the (M' x 1) vector belonging to other bras' representation.
     procedure :: CastMatrixBra_FromVector

     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to other kets' representation.
     generic   :: CastMatrixKet => &
          CastMatrixKet_FromClassMatrix, &
          CastMatrixKet_FromMatrix, &
          CastMatrixKet_FromVector
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to other kets' representation. The original matrix is in [ClassMatrix](@ref ModuleMatrix.f90) format.
     procedure :: CastMatrixKet_FromClassMatrix
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to the total basis kets' representation.
     procedure :: CastMatrixKetToTotal
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to the regular basis kets' representation.
     procedure :: CastMatrixKetToRegular
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to the preconditioned basis kets' representation.
     procedure :: CastMatrixKetToPreconditioned
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to the continuous preconditioned basis kets' representation.
     procedure :: CastMatrixKetToPreconditionedContinuous
     !
     !> Casts an original (M x N) matrix in a certain representation to obtain the (M x N') matrix belonging to other kets' representation. The original matrix is stored as a 2-dimensional array.
     procedure :: CastMatrixKet_FromMatrix
     !
     !> Casts an original (1 x N) vector in a certain representation to obtain the (1 x N') vector belonging to other kets' representation.
     procedure :: CastMatrixKet_FromVector

     procedure, private :: BasisFunctionEval
     procedure, private :: BasisExpansionEval
     procedure, private :: BasisFunctionTabulate
     procedure, private :: BasisExpansionTabulate
     procedure, private :: MonoelectronicMultipole
     procedure, private :: BielectronicMultipole
     procedure, private :: MatchKindOnUnit

  end type ClassBasis


contains

  !> Retrieves the number of points that will be used 
  !! to plot the radial basis functions.
  integer function GetNPlot(Basis)
    class(ClassBasis), intent(in) :: Basis
    GetNPlot=Basis.NPlot
  end function GetNPlot


  !> Performs the plot of the radial basis.
  subroutine PlotBasis(Basis, Dir)
    !
    class(ClassBasis), intent(inout) :: Basis
    character(len=*) , intent(in)    :: Dir
    !
    DoublePrecision, allocatable :: XGrid(:)
    DoublePrecision, allocatable :: YGrid(:)
    integer             :: ix, npts, iFun
    integer             :: uid, uidAll
    character(len=1000) :: FileName
    character(len=IOMSG_LENGTH) :: iomsg
    integer                     :: iostat

    !.. Set plot grid
    !..
    call Basis.GetPlotGrid( XGRID )

    !.. Tabulate functions and write them on disk
    !..
    FileName=trim(Dir)//"All"
    open(Newunit =  uidAll    , &
         File    =  FileName  , &
         Status  = "unknown"  , &
         Action  = "write"    , &
         Form    = "formatted", &
         iostat  =  iostat    , &
         iomsg   =  iomsg     )
    if(iostat/=0)then
       call ErrorMessage(iomsg)
       return
    endif
    !
    npts = Basis.GetNPlot()
    allocate( YGrid, source = XGrid )
    YGrid = 0.d0
    !
    do iFun = 1, Basis.GetNFun()
       !
       FileName=trim(Dir)//AlphabeticNumber(iFun,Basis.GetNFun(),"0")
       open(Newunit =  uid       , &
            File    =  FileName  , &
            Status  = "unknown"  , &
            Action  = "write"    , &
            Form    = "formatted", &
            iostat  =  iostat    , &
            iomsg   =  iomsg     )
       if(iostat/=0)then
          call ErrorMessage(iomsg)
          exit
       endif
       !
       call Basis.Tabulate( npts, XGrid, YGrid, iFun )
       !
       do ix=1,npts
          write(uid   ,"(2(x,e24.16E3))") XGrid(ix), YGrid(ix)
          write(uidAll,"(2(x,e24.16E3))") XGrid(ix), YGrid(ix)
       enddo
       !
       write(uidAll,*)
       !
       close(uid)
       !
    enddo
    !
    close(uidAll)
    !
    deallocate( XGrid, YGrid )
    !
  end subroutine PlotBasis


  !> Retrieves whether the radial basis functions will be plotted or not.
  logical function ToBePlotted(Basis)
    class(ClassBasis), intent(in) :: Basis
    ToBePlotted=Basis.PlotRequired
  end function ToBePlotted

  !> Builds the directory in which the results will be stored.
  subroutine SetRoot( self, RootDir )
    class(ClassBasis), intent(inout) :: self
    character(len=*),  intent(in) :: RootDir
    character(len=:), allocatable :: strn
    allocate(strn,source=AddSlash(RootDir)//self.Dir)
    self.Dir=strn
  end subroutine SetRoot

  !> Retrieves the directory in wich the radial basis 
  !! information will be stored.
  function GetDir(Basis) result(Dir)
    class(ClassBasis), intent(in) :: Basis
    character(len=:), allocatable :: Dir
    integer :: len
    len=len_trim(adjustl(Basis.Dir))
    if(Basis.Dir(len:len)/="/")then
       allocate(Dir,source=trim(adjustl(Basis.Dir))//"/")
    else
       allocate(Dir,source=trim(adjustl(Basis.Dir)))
    endif
  end function GetDir


  !> Determines the number of subdiagonals of the operator 
  !! matrices computed using the radial basis.
  integer function LowerBandwidth(Basis) result( Bandwidth )
    class(ClassBasis), intent(in) :: Basis
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       Bandwidth=ptr.GetOrder()-1
    class is (ClassGaussian)
       Bandwidth=ptr.GetNFun()-1
    class is (ClassGaussianBSpline)
       Bandwidth=ptr.GetNFun()-1
    class DEFAULT
    end select
  end function LowerBandwidth


  !> Gets the number of B-splines to be skipped.
  integer function GetNBSplineSkip(Basis) result( res )
    class(ClassBasis), intent(in) :: Basis
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       res = 1
    class is (ClassGaussian)
       res = 0
    class is (ClassGaussianBSpline)
       res = ptr.GetNBSplineSkip()
    class DEFAULT
    end select
  end function GetNBSplineSkip


  !> Determines the number of superdiagonals of the operator
  !! matrices computed using the radial basis.
  integer function UpperBandwidth(Basis) result( Bandwidth )
    class(ClassBasis), intent(in) :: Basis
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       Bandwidth=ptr.GetOrder()-1
    class is (ClassGaussian)
       Bandwidth=ptr.GetNFun()-1
    class is (ClassGaussianBSpline)
       Bandwidth=ptr.GetNFun()-1
    class DEFAULT
    end select
  end function UpperBandwidth


  !> Depending on the radial basis kind, selects the best 
  !! matrix pattern to store the local operators matrix 
  !! elements: Full or Banded.
  function BestPattern(Basis) result(Pattern)
    class(ClassBasis), intent(in) :: Basis
    character(len=:), allocatable :: Pattern
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       allocate(Pattern,source=BSPLINE_BEST_PATTERN)
    class is (ClassGaussian)
       allocate(Pattern,source=GAUSSIAN_BEST_PATTERN)
    class is (ClassGaussianBSpline)
       allocate(Pattern,source=GAUSSIANBSPLINE_BEST_PATTERN)
    class DEFAULT
    end select
  end function BestPattern


  !> Checks whether ClassBasis has been initialized or not.
  logical function BasisIsInitialized(Basis) &
       result(Initialized)
    class(ClassBasis), intent(in) :: Basis
    Initialized=associated(Basis.RadialBasis)
  end function BasisIsInitialized


  !> Computes the "Spherical Bessel Transform" of a radial basis function \f$G(r)\f$, defined as:
  !! \f[
  !!   F_{l}(p)=\int\frac{G(r)}{r}j_{l}(p\cdot r)r^{2}dr,
  !! \f]
  !! for all the angular momenta \f$\ell\f$ up to a specified lmax as \f$ res(\ell+1)=F_\ell(p)\f$.
  subroutine SphericalBesselTransform(Basis,p,Element,lmax,res)
    Class(ClassBasis), intent(in) :: Basis
    DoublePrecision  , intent(in) :: p
    integer          , intent(in) :: Element,lmax
    DoublePrecision  , intent(out):: res(*)
    integer :: Bs
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs=Element+1
       call ptr.SphericalBesselTransform(p,Bs,lmax,res)
    class is(ClassGaussian)
       res(1:LMAX+1)=0.d0
!!$       call ptr.SphericalBesselTransform(p,Element,lmax,res)
    class is(ClassGaussianBSpline)
       res(1:LMAX+1)=0.d0
!!$       call ptr.SphericalBesselTransform(p,Element,lmax,res)
    class DEFAULT
    end select
  end subroutine SphericalBesselTransform


  !> Computes the monoelectronic multipole.
  !! \f[
  !!    \int_0^\infty dr G_n(R)\frac{r_<^\ell}{r_>^{\ell+1}}G_m(r)
  !! \f]
  !! where \f$r_<=\min(r,R)\f$ and \f$r_>=\max(r,R)\f$. 
  DoublePrecision function MonoElectronicMultipole(Basis,iBra,iKet,R,l)&
       result(Multipole)
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: iBra,iKet
    DoublePrecision  , intent(in) :: R
    integer          , intent(in) :: l
    integer :: Bs1, Bs2
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs1 = iBra + 1
       Bs2 = iKet + 1
       Multipole = ptr.Multipole(Bs1,Bs2,R,l)
    class is(ClassGaussian)
       Multipole = ptr.Multipole(iBra,iKet,R,l)
    class is(ClassGaussianBSpline)
       Multipole = ptr.Multipole(iBra,iKet,R,l)
    class DEFAULT
    end select
  end function MonoElectronicMultipole


  !> Computes the bielectronic multipole
  !! \f[
  !!     \int_0^\infty dr\int_0^{\infty} dr' 
  !!      G_{n_1}(r) G_{n_2}(r') 
  !!      \frac{r_<^\ell}{r_>^{\ell+1}}
  !!      G_{n_3}(r) G_{n_4}(r') 
  !! \f]
  !! between non-normalized \f$G(r)f$ radial basis functions.
  DoublePrecision function BiElectronicMultipole(Basis,n1,n2,n3,n4,l)&
       result(Multipole)
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: n1,n2,n3,n4
    integer          , intent(in) :: l
    integer :: Bs1, Bs2, Bs3, Bs4
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs1=n1+1
       Bs2=n2+1
       Bs3=n3+1
       Bs4=n4+1
       Multipole=ptr.Multipole(Bs1,Bs2,Bs3,Bs4,l)
    class is(ClassGaussian)
       Multipole=0.d0
    class is(ClassGaussianBSpline)
       Multipole=0.d0
    class DEFAULT
    end select
  end function BiElectronicMultipole


  !> Return the basis function normalization factor.
  DoublePrecision function GetNormFactor(Basis,Element) &
       result( NormFactor )
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: Element
    integer :: Bs
    NormFactor=0.d0
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs=Element+1
       NormFactor=ptr.GetNormFactor(Bs)
    class is(ClassGaussian)
       NormFactor=ptr.GetNormFactor(Element)
    class is(ClassGaussianBSpline)
       NormFactor=ptr.GetNormFactor(Element)
    class DEFAULT
    end select
  end function GetNormFactor


  !> Computes the principal-part radial basis functions integral
  !! \f[
  !!  \int dr G_n(r)\,\frac{\mathcal{P}}{r-r_0}\,G_m(r)
  !! \f]
  DoublePrecision function CauchyIntegral(Basis,iBra,iKet,SingularPoint) &
       result( Integral )
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: iBra, iKet
    DoublePrecision  , intent(in) :: SingularPoint
    !
    integer :: Bs1,Bs2
    !
    Integral=0.d0
    !
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs1=iBra+1
       Bs2=iKet+1
       Integral=ptr.CauchyIntegral(Bs1,Bs2,SingularPoint)
    class is(ClassGaussian)
       Integral=0.d0
    class is(ClassGaussianBSpline)
       Integral=0.d0
    class DEFAULT
    end select
    !
  end function CauchyIntegral



  !> Computes the integral
  !! \f[
  !!    \int_a^b dr 
  !!    \frac{d^{n_1}G(r)}{dr^{n_1}} f(r)
  !! \f]
  !! where \f$f(r)\f$ is a *smooth* function in the radial basis support.
  !! If the first (last) boundary of the integration interval is not specified, 
  !! then the value \f$a=-\infty\f$ and \f$b=\infty\f$ is assumed.
  DoublePrecision function MonoIntegral(&
       Basis          , &
       FunPtr         , &
       Element        , &
       DerivativeOrder, &
       LowerBound     , &
       UpperBound     ) &
       result( Integral )
    !
    Class(ClassBasis)        , intent(in) :: Basis
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: Element
    integer        , optional, intent(in) :: DerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    !
    DoublePrecision :: a,b
    integer :: Bs,n
    !
    Integral=0.d0
    !
    n = 0; if(present(DerivativeOrder)) n = DerivativeOrder
    a = -huge(1.d0); if(present(LowerBound)) a = LowerBound
    b =  huge(1.d0); if(present(UpperBound)) b = UpperBound
    !
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs=Element+1
       Integral=ptr.MonoIntegral(FunPtr,Bs,n,a,b)
    class is(ClassGaussian)
       !
       ! ***
       !
       Integral=0.d0
    class is(ClassGaussianBSpline)
       !
       ! ***
       !
       Integral=0.d0
    class DEFAULT
    end select
    !
  end function MonoIntegral


  !> Computes the  k-th derivatives of a choosen radial basis function evaluated in a position array. 
  subroutine BasisFunctionTabulate(Basis,ndata,xvec,yvec,Element,nDer) 
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: ndata
    DoublePrecision  , intent(in) :: xvec(1:ndata)
    DoublePrecision  , intent(out):: yvec(1:ndata)
    integer          , intent(in) :: Element
    integer, optional, intent(in) :: nDer
    !
    integer :: n,Bs
    n=0;if(present(nDer))n=nDer
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       Bs=Element+1
       call ptr.Tabulate(ndata,xvec,yvec,Bs,n)
    class is (ClassGaussian)
       call ptr.Tabulate(ndata,xvec,yvec,Element,n)
    class is (ClassGaussianBSpline)
       call ptr.Tabulate(ndata,xvec,yvec,Element,n)
    class DEFAULT
       yvec=0.d0
    end select
  end subroutine BasisFunctionTabulate


  !> Computes the  k-th derivative of a function expressed in terms of 
  !! the radial basis and evaluated in a position array.
  subroutine BasisExpansionTabulate(Basis,ndata,xvec,yvec,FunVec,nDer) 
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: ndata
    DoublePrecision  , intent(in) :: xvec(1:ndata)
    DoublePrecision  , intent(out):: yvec(1:ndata)
    DoublePrecision  , intent(in) :: FunVec(:)
    integer, optional, intent(in) :: nDer
    !
    integer :: n
    n = 0; if( present( nDer ) ) n = nDer
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       call ptr.Tabulate(ndata,xvec,yvec,FunVec,n,SKIP_FIRST=.TRUE.)
    class is (ClassGaussian)
       call ptr.Tabulate(ndata,xvec,yvec,FunVec,n)
    class is (ClassGaussianBSpline)
       call ptr.Tabulate(ndata,xvec,yvec,FunVec,n)
    class DEFAULT
       yvec=0.d0
    end select
  end subroutine BasisExpansionTabulate

  !> Gets the number of radial basis used in the computation of the operators matrix elements.
  integer function GetNFun(Basis, Identifier)
    class(ClassBasis), intent(in) :: Basis
    character(len=*), optional, intent(in):: Identifier
    GetNFun=0
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       if ( (present(Identifier)) .and. (Identifier .is. "B") ) then
          !
          GetNFun = ptr.GetNBSplines() - 1
          !
       elseif ( (present(Identifier)) .and. (Identifier .is. "G") ) then
          !
          GetNFun = 0
          !
       else
          !
          GetNFun=ptr.GetNBSplines()-1
          !
       end if
    class is (ClassGaussian)
       if ( (present(Identifier)) .and. (Identifier .is. "B") ) then
          !
          GetNFun = 0
          !
       elseif ( (present(Identifier)) .and. (Identifier .is. "G") ) then
          !
          GetNFun=ptr.GetNFun()
          !
       else
          !
          GetNFun=ptr.GetNFun()
          !
       end if
    class is (ClassGaussianBSpline)
       if ( (present(Identifier)) .and. (Identifier .is. "B") ) then
          !
          GetNFun = ptr.GetBsplinesNFun()
          !
       elseif ( (present(Identifier)) .and. (Identifier .is. "G") ) then
          !
          GetNFun=ptr.GetGaussianNFun()
          !
       else
          !
          GetNFun=ptr.GetNFun()
          !
       end if
    class DEFAULT
    end select
  end function GetNFun




  !> Gets the information of the linear independent gaussian function (preconditioned): monomial exponent, gaussian exponent and normalization factor.
  subroutine GetLinIndGaussian( Basis, MonExp, GaussExp, NormFactor )
    !> Radial basis class.
    class(ClassBasis), intent(in) :: Basis
    !> Monomial exponents.
    integer, allocatable, intent(out) :: MonExp(:)
    !> Gaussian exponents.
    real(kind(1d0)), allocatable, intent(out) :: GaussExp(:)
    !> Normalization factors.
    real(kind(1d0)), allocatable, intent(out) :: NormFactor(:)
    !
    select type(ptr=>Basis.RadialBasis)
    class is (ClassBSpline)
       call Assert("The linear independent Gaussian functions cannot be fetched from a B-splines radial basis.")
    class is (ClassGaussian)
       call ptr.GetLIGaussian( MonExp, GaussExp, NormFactor )
    class is (ClassGaussianBSpline)
       call ptr.GetLinIndGaussian( MonExp, GaussExp, NormFactor )
    class DEFAULT
    end select
  end subroutine GetLinIndGaussian


  !> Computes the  k-th derivatives of a choosen radial basis function evaluated in certain position. 
  DoublePrecision function BasisFunctionEval(Basis,x,Element,nDer) 
    Class(ClassBasis), intent(in) :: Basis
    DoublePrecision  , intent(in) :: x
    integer          , intent(in) :: Element
    integer, optional, intent(in) :: nDer
    integer :: Bs,n
    BasisFunctionEval=0.d0
    n=0;if(present(nDer))n=nDer
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs = Element + 1
       BasisFunctionEval=ptr.Eval(x,Bs,n)
    class is(ClassGaussian)
       BasisFunctionEval=ptr.Eval(x,Element,n)
    class is(ClassGaussianBSpline)
       BasisFunctionEval=ptr.Eval(x,Element,n)
    class DEFAULT
    end select
  end function BasisFunctionEval


  !> Computes the  k-th derivatives of a function expressed in terms 
  !! of the radial basis evaluated in certain position.
  DoublePrecision function BasisExpansionEval(Basis,x,fc,nDer, Identifier)
    class(ClassBasis), intent(in) :: Basis
    DoublePrecision  , intent(in) :: x
    DoublePrecision  , intent(in) :: fc(:)
    integer, optional, intent(in) :: nDer
    character(len=*), optional, intent(in) :: Identifier
    integer :: n
    BasisExpansionEval=0.d0
    n=0;if(present(nDer))n=nDer
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       BasisExpansionEval=ptr.Eval(x,fc,n,SKIP_FIRST = .TRUE. )
    class is(ClassGaussian)
       BasisExpansionEval=ptr.Eval(x,fc,n)
    class is(ClassGaussianBSpline)
       BasisExpansionEval=ptr.Eval(x,fc,n,Identifier)
    class DEFAULT
    end select
  end function BasisExpansionEval


  !> Frees the ClassBasis attributes.
  subroutine BasisFree(Basis)
    Class(ClassBasis), intent(inout) :: Basis
    if(allocated(Basis.Dir))deallocate(Basis.Dir)
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       call ptr.Free()
    class is(ClassGaussian)
       call ptr.Free()
    class is(ClassGaussianBSpline)
       call ptr.Free()
    class DEFAULT
    end select
    if(associated(Basis.RadialBasis))then
       deallocate(Basis.RadialBasis)
       Basis.RadialBasis => NULL()
    endif
    Basis.NPlot=0
    Basis.PlotRequired=.FALSE.
  end subroutine BasisFree



  !> Computes the integral
  !! \f[
  !!   \int_{a}^{b}
  !!     \frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)
  !!     \frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr,
  !!\f]
  !! where \f$f(r)\f$ is a local operator and \f$G_{i,j}\f$ 
  !! radial basis functions. If a break point \f$BP\f$ is 
  !! introduced, then the integral is splitted in two parts
  !! \f[
  !!   \int_{a}^{BP}
  !!     \frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)
  !!     \frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr + 
  !!   \int_{BP}^{b}
  !!     \frac{d^{n_{1}}G_{i}(r)}{dr^{n_{1}}}f(r)
  !!     \frac{d^{n_{2}}G_{j}(r)}{dr^{n_{2}}}r^{2}dr.
  !! \f]
  DoublePrecision function Integral( &
       Basis             , &
       FunPtr            , & 
       iBra              , &
       iKet              , &
       BraDerivativeOrder, &
       KetDerivativeOrder, &
       LowerBound        , &
       UpperBound        , &
       BreakPoint        )
    !
    Class(ClassBasis)        , intent(in) :: Basis
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: iBra
    integer                  , intent(in) :: iKet
    integer        , optional, intent(in) :: BraDerivativeOrder
    integer        , optional, intent(in) :: KetDerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    DoublePrecision, optional, intent(in) :: BreakPoint
    !
    integer :: Bs1,Bs2,n1,n2
    DoublePrecision :: a,b,c
    !
    Integral=0.d0
    !
    n1 = 0; if(present(BraDerivativeOrder)) n1 = BraDerivativeOrder
    n2 = 0; if(present(KetDerivativeOrder)) n2 = KetDerivativeOrder
    a = -huge(1.d0); if(present(LowerBound)) a = LowerBound
    b =  huge(1.d0); if(present(UpperBound)) b = UpperBound
    c = -huge(1.d0); if(present(BreakPoint)) c = BreakPoint
    !
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs1=iBra+1
       Bs2=iKet+1
       Integral=ptr.Integral(FunPtr,Bs1,Bs2,n1,n2,a,b,c)
    class is(ClassGaussian)
       Integral=ptr.Integral(FunPtr,iBra,iKet,n1,n2,a,b,c)
    class is(ClassGaussianBSpline)
       Integral=ptr.Integral(FunPtr,iBra,iKet,n1,n2,a,b,c)
    class DEFAULT
    end select
    !
  end function Integral


  subroutine IntegralRow( &
       Basis             , &
       FunPtr            , & 
       iBra              , &
       iKetMin           , &
       iKetMax           , &
       IntegralVec       , &
       BraDerivativeOrder, &
       KetDerivativeOrder, &
       LowerBound        , &
       UpperBound        , &
       BreakPoint        )
    !
    Class(ClassBasis)        , intent(in) :: Basis
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: iBra
    integer                  , intent(in) :: iKetMin
    integer                  , intent(in) :: iKetMax
    DoublePrecision          , intent(out):: IntegralVec(:)
    integer        , optional, intent(in) :: BraDerivativeOrder
    integer        , optional, intent(in) :: KetDerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    DoublePrecision, optional, intent(in) :: BreakPoint
    !
    integer :: Bs1,Bs2,n1,n2,iKet
    DoublePrecision :: a,b,c, Integral
    !
    IntegralVec(1:iKetMax-iKetMin+1)=0.d0
    !
    n1 = 0; if(present(BraDerivativeOrder)) n1 = BraDerivativeOrder
    n2 = 0; if(present(KetDerivativeOrder)) n2 = KetDerivativeOrder
    a = -huge(1.d0); if(present(LowerBound)) a = LowerBound
    b =  huge(1.d0); if(present(UpperBound)) b = UpperBound
    c = -huge(1.d0); if(present(BreakPoint)) c = BreakPoint
    !
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs1=iBra+1
       do iKet=iKetMin,iKetMax
          Bs2=iKet+1
          IntegralVec(iKet-iKetMin+1)=ptr.Integral(FunPtr,Bs1,Bs2,n1,n2,a,b,c)
       enddo
    class is(ClassGaussian)
       do iKet=iKetMin,iKetMax
          IntegralVec(iKet-iKetMin+1)=ptr.Integral(FunPtr,iBra,iKet,n1,n2,a,b,c)
       enddo
    class is(ClassGaussianBSpline)
       do iKet=iKetMin,iKetMax
          Integral=ptr.Integral(FunPtr,iBra,iKet,n1,n2,a,b,c)
          IntegralVec(iKet-iKetMin+1)=Integral
          if( ptr.RemainingKetsGiveZero( iBra, iKet, Integral ) )exit
       enddo
    class DEFAULT
    end select
    !
  end subroutine IntegralRow



  subroutine IntegralCol( &
       Basis             , &
       FunPtr            , & 
       iBraMin           , &
       iBraMax           , &
       iKet              , &
       IntegralVec       , &
       BraDerivativeOrder, &
       KetDerivativeOrder, &
       LowerBound        , &
       UpperBound        , &
       BreakPoint        )
    !
    Class(ClassBasis)        , intent(in) :: Basis
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: iBraMin
    integer                  , intent(in) :: iBraMax
    integer                  , intent(in) :: iKet
    DoublePrecision          , intent(out):: IntegralVec(:)
    integer        , optional, intent(in) :: BraDerivativeOrder
    integer        , optional, intent(in) :: KetDerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    DoublePrecision, optional, intent(in) :: BreakPoint
    !
    integer :: Bs1,Bs2,n1,n2,iBra
    DoublePrecision :: a,b,c, Integral
    !
    IntegralVec(1:iBraMax-iBraMin+1)=0.d0
    !
    n1 = 0; if(present(BraDerivativeOrder)) n1 = BraDerivativeOrder
    n2 = 0; if(present(KetDerivativeOrder)) n2 = KetDerivativeOrder
    a = -huge(1.d0); if(present(LowerBound)) a = LowerBound
    b =  huge(1.d0); if(present(UpperBound)) b = UpperBound
    c = -huge(1.d0); if(present(BreakPoint)) c = BreakPoint
    !
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Bs2=iKet+1
       do iBra=iBraMin,iBraMax
          Bs1=iBra+1
          IntegralVec(iBra-iBraMin+1)=ptr.Integral(FunPtr,Bs1,Bs2,n1,n2,a,b,c)
       enddo
    class is(ClassGaussian)
       do iBra=iBraMin,iBraMax
          IntegralVec(iBra-iBraMin+1)=ptr.Integral(FunPtr,iBra,iKet,n1,n2,a,b,c)
       enddo
    class is(ClassGaussianBSpline)
       do iBra=iBraMin,iBraMax
          Integral=ptr.Integral(FunPtr,iBra,iKet,n1,n2,a,b,c)
          IntegralVec(iBra-iBraMin+1)=Integral
          if( ptr.RemainingBrasGiveZero( iBra, iKet, Integral ) )exit
       enddo
    class DEFAULT
    end select
    !
  end subroutine IntegralCol



  !> Creates in a directory all the configuration
  !!    files that define the current basis.
  subroutine SaveConfig(Basis, ConfigurationFile, Dir)
    !
    Class(ClassBasis), intent(inout) :: Basis
    character(len=*) , intent(in)    :: ConfigurationFile
    character(len=*) , intent(in)    :: Dir
    !
    character(len=1000) :: RadialBasisFile
    !
    call Basis.List.Add( "RadialBasisFile", "A"   , "required" )
    call Basis.List.Parse( ConfigurationFile )
    call Basis.List.Get( "RadialBasisFile", RadialBasisFile    ) 
    !
    RadialBasisFile = adjustl(RadialBasisFile)
    !
    select type( ptr => Basis.RadialBasis )
    class is( ClassBSpline )
       call ptr.SaveConfig(trim(RadialBasisFile), trim(Dir))
    class is( ClassGaussian )
       call ptr.SaveConfig(trim(RadialBasisFile), trim(Dir))
    class is( ClassGaussianBSpline )
       call ptr.SaveConfig(trim(RadialBasisFile), trim(Dir))
    end select
    !
  end subroutine SaveConfig


  !> Initializes ClassBasis reading the information from a configuration file.
  subroutine BasisInitReadConfigurationFile(Basis,ConfigurationFile,IOSTAT)
    !
    Class(ClassBasis), intent(inout)         :: Basis
    character(len=*) , intent(in)            :: ConfigurationFile
    integer          , intent(out), optional :: IOSTAT
    !
    character(len=1000) :: BasisDir
    character(len=100)  :: RadialBasisKind
    character(len=1000) :: RadialBasisFile
    integer :: status
    !
    if(present(IOSTAT))IOSTAT=0
    !
    !
    !.. Parse the configuration file
    !..
    call Basis.List.Add( "BasisKind"      , "A"   , "required" )
    call Basis.List.Add( "RadialBasisFile", "A"   , "required" )
    call Basis.List.Add( "BasisDir"       , "A"   , "optional" )
    call Basis.List.Add( "PlotBasis"      , .TRUE., "optional" )
    call Basis.List.Add( "Nplot"          , 500   , "optional" )
    call Basis.List.Add( "RPlotMin"       , 0.d0  , "optional" )
    call Basis.List.Add( "RPlotMax"       , 1.d0  , "optional" )
    !
    call Basis.List.Parse( ConfigurationFile )
    !
    call Basis.List.Get( "BasisKind"      , RadialBasisKind    )
    call Basis.List.Get( "RadialBasisFile", RadialBasisFile    ) 
    call Basis.List.Get( "BasisDir"       , BasisDir           )
    call Basis.List.Get( "PlotBasis"      , Basis.PlotRequired )
    call Basis.List.Get( "Nplot"          , Basis.Nplot        )
    allocate(Basis.Dir,source=trim(adjustl(BasisDir)))
    !
    !.. Determines which kind the basis is
    !..
    RadialBasisKind=adjustl(RadialBasisKind)
    select case( trim(adjustl(RadialBasisKind)) )
    case( BSPLINE_STRN_IDENTIFIER )
       allocate(ClassBSpline::Basis.RadialBasis)
    case( GAUSSIAN_STRN_IDENTIFIER )
       allocate(ClassGaussian::Basis.RadialBasis)
    case( GAUSSIANBSPLINE_STRN_IDENTIFIER )
       allocate(ClassGaussianBSpline::Basis.RadialBasis)
    case DEFAULT
       call ErrorMessage("Unrecognized Basis Kind")
       STOP
    end select
    !
    select type( ptr => Basis.RadialBasis )
    class is( ClassBSpline )
       call ptr.Init(RadialBasisFile,status)
    class is( ClassGaussian )
       call ptr.Init(RadialBasisFile,status)
    class is( ClassGaussianBSpline )
       call ptr.Init(RadialBasisFile,status)
    end select

    Basis.RPlotMin = 0.d0
    if( Basis.List.present("RPlotMin"))then
       call Basis.List.Get("RPlotMin", Basis.RPlotMin)
    endif
    Basis.RPlotMax = Basis.GetRMax()
    if( Basis.List.present("RPlotMax"))then
       call Basis.List.Get("RPlotMax", Basis.RPlotMax)
    endif

    if(status/=0)then
       if(present(IOSTAT))then
          IOSTAT=status
       else
          call ErrorMessage("Basis Initialization failed")
       end if
    endif
    !
  end subroutine BasisInitReadConfigurationFile


  !> Writes in a unit the adequate radial basis information 
  !! (fingerprint) that allows to identify it from other 
  !! basis that could possible be used.
  subroutine WriteFingerPrint(Basis,uid)
    class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: uid
    call Basis.WriteKind(uid)
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       call ptr.WriteFingerprint(uid)
    class is(ClassGaussian)
       call ptr.WriteFingerprint(uid)
    class is(ClassGaussianBSpline)
       call ptr.WriteFingerprint(uid)
    class DEFAULT
    end select
  end subroutine WriteFingerPrint


  !> Checks whether tha radial basis characteristics (fingerprint) 
  !! saved previously in a unit, coincides or not with the current 
  !! radial basis in used.
  logical function MatchFingerPrintOnUnit(Basis,uid) result(Match)
    Class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: uid
    Match=Basis.MatchKindOnUnit(uid)
    if(.not.Match)return
    select type(ptr=>Basis.RadialBasis)
    class is(ClassBSpline)
       Match=ptr.MatchFingerPrintOnUnit(uid)
    class is(ClassGaussian)
       Match=ptr.MatchFingerPrintOnUnit(uid)
    class is(ClassGaussianBSpline)
       Match=ptr.MatchFingerPrintOnUnit(uid)
    class DEFAULT
       Match=.FALSE.
    end select
  end function MatchFingerPrintOnUnit


  !> Writes in a unit the radial basis kind (B-splines, Gaussian, etc.).
  subroutine BasisWriteKind(Basis,OutUnit)
    class(ClassBasis), intent(in) :: Basis
    integer, optional, intent(in) :: OutUnit
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Writable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    INQUIRE(&
         UNIT  = OutUnit , &
         OPENED= Opened  , &
         WRITE = Writable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Output Unit is closed")
    if( trim(Writable) /= "YES" ) call Assert("Output Unit can't be written")
    select case( Form )
    case( "FORMATTED" )
       select type(ptr => Basis.RadialBasis)
       class is (ClassBSpline)
          write(OutUnit,"(i)") BASIS_KIND_BSPLINE
       class is(ClassGaussian)
          write(OutUnit,"(i)") BASIS_KIND_GAUSSIAN
       class is(ClassGaussianBSpline)
          write(OutUnit,"(i)") BASIS_KIND_GAUSSIANBSPLINE
       class DEFAULT
       end select
    case( "UNFORMATTED" )
       select type(ptr => Basis.RadialBasis)
       class is (ClassBSpline)
          write(OutUnit) BASIS_KIND_BSPLINE
       class is (ClassGaussian)
          write(OutUnit) BASIS_KIND_GAUSSIAN
       class is (ClassGaussianBSpline)
          write(OutUnit) BASIS_KIND_GAUSSIANBSPLINE
       class DEFAULT
       end select
    case DEFAULT
       call Assert("Invalid unit format")
    end select
  end subroutine BasisWriteKind


  !> Checks whether the radial basis kind saved previously 
  !! in a unit, coincides or not with the kind currently used.
  logical function MatchKindOnUnit(Basis,InUnit) &
       result(Match)
    class(ClassBasis), intent(in) :: Basis
    integer, optional, intent(in) :: InUnit
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Readable
    character(len=16) :: Form
    integer :: iostat
    logical :: Opened
    integer :: KindOnUnit
    Match=.FALSE.
    INQUIRE(&
         UNIT  = InUnit  , &
         OPENED= Opened  , &
         READ  = Readable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Input Unit is closed")
    if( trim(Readable) /= "YES" ) call Assert("Input Unit can't be read")
    select case( Form )
    case( "FORMATTED" )
       read(InUnit,"(i)",iostat=iostat,iomsg=iomsg) KindOnUnit
       if(iostat/=0) call Assert(iomsg)
       select type(ptr => Basis.RadialBasis)
       class is (ClassBSpline)
          Match = KindOnUnit == BASIS_KIND_BSPLINE
       class is (ClassGaussian)
          Match = KindOnUnit == BASIS_KIND_GAUSSIAN
       class is (ClassGaussianBSpline)
          Match = KindOnUnit == BASIS_KIND_GAUSSIANBSPLINE
       class DEFAULT
       end select
    case( "UNFORMATTED" )
       read(InUnit,iostat=iostat,iomsg=iomsg) KindOnUnit
       if(iostat/=0) call Assert(iomsg)
       select type(ptr => Basis.RadialBasis)
       class is (ClassBSpline)
          Match = KindOnUnit == BASIS_KIND_BSPLINE
       class is (ClassGaussian)
          Match = KindOnUnit == BASIS_KIND_GAUSSIAN
       class is (ClassGaussianBSpline)
          Match = KindOnUnit == BASIS_KIND_GAUSSIANBSPLINE
       class DEFAULT
       end select
    case DEFAULT
       call Assert("Invalid unit format")
    end select
  end function MatchKindOnUnit


  !> Gets the total number of radial basis functions.
  integer function GetTotalNumBasisFunctions( Basis ) result(TotalNumBasisFunctions)
    class(ClassBasis), intent(in) :: Basis
    select type ( ptr => Basis.RadialBasis )
    class is ( ClassBSpline )
       TotalNumBasisFunctions = ptr.GetNBSplines()
    class is ( ClassGaussian )
       TotalNumBasisFunctions = ptr.GetNFun()
    class is ( ClassGaussianBSpline )
       TotalNumBasisFunctions = ptr.GetNFun()
    class DEFAULT
       ! Fill with the gaussian part
    end select
  end function GetTotalNumBasisFunctions


  !> Retrieves the minimum regular index for BSplines starting
  !! at the origin. Apply to the pure B-spline basis only.
  integer function MinimumRegularBSplineIndex( BSplineBasis, L)
    type(ClassBSpline), intent(in) :: BSplineBasis
    integer,            intent(in) :: L
    MinimumRegularBSplineIndex = min( L+1, BSplineBasis.GetOrder() )
  end function MinimumRegularBSplineIndex
  !> Retrieves the maximum regular index for BSplines 
  integer function MaximumRegularBSplineIndex( BSplineBasis )
    type(ClassBSpline), intent(in) :: BSplineBasis
    MaximumRegularBSplineIndex = BSplineBasis.GetNBSplines()-1 
  end function MaximumRegularBSplineIndex




  !> Retrieves the vector of regular indexes for a given angular momentum
  subroutine GetRegularIndexes( Basis, RegularIndexes, NRegularIndexes, L )
    class(ClassBasis)   , intent(in)  :: Basis
    integer, allocatable, intent(out) :: RegularIndexes(:)
    integer             , intent(out) :: NRegularIndexes
    integer             , intent(in)  :: L
    !
    integer :: i,imi,ima
    !
    select type ( ptr => Basis.RadialBasis )

    class is ( ClassBSpline )

       imi = MinimumRegularBSplineIndex(ptr,L)
       ima = MaximumRegularBSplineIndex(ptr)
       NRegularIndexes = ima - imi + 1
       
       if(allocated(RegularIndexes))deallocate(RegularIndexes)
       allocate(RegularIndexes(NRegularIndexes))
       RegularIndexes=[(i,i=imi,ima)]
       
    class is ( ClassGaussian )
       
       imi = ptr.GetMinimumRegularIndex(L)
       ima = ptr.GetMaximumRegularIndex(L)
       NRegularIndexes = ima - imi + 1
       
       if(allocated(RegularIndexes))deallocate(RegularIndexes)
       allocate(RegularIndexes(NRegularIndexes))
       RegularIndexes=[(i,i=imi,ima)]       
       
    class is ( ClassGaussianBSpline )
       
       call ptr.GetRegularIndexes( RegularIndexes, NRegularIndexes, L )

    class DEFAULT
    end select

  end subroutine GetRegularIndexes


  !> Retrieves the vector of regular indexes for a given angular momentum
  integer function NRegular( Basis, L ) 
    class(ClassBasis), intent(in)  :: Basis
    integer          , intent(in)  :: L
    !
    integer :: imi, ima
    select type ( ptr => Basis.RadialBasis )
       
    class is ( ClassBSpline )
       
       !.. For Bsplines, NRegular is indeed 1 unit lower than the number 
       !   of columns of the corresponding operator's representation matrix.
       !..
       imi = MinimumRegularBSplineIndex(ptr,L)
       ima = MaximumRegularBSplineIndex(ptr)
       NRegular = ima - imi + 1
       
    class is ( ClassGaussian )
       
       imi = ptr.GetMinimumRegularIndex(L)
       ima = ptr.GetMaximumRegularIndex(L)
       NRegular = ima - imi + 1
       
    class is ( ClassGaussianBSpline )
       
       NRegular = ptr.NRegular(L)
       
    class DEFAULT
    end select

  end function NRegular


  !> Retrieves the vector of regular indexes for a given angular momentum
  integer function NPreconditioned( self, L, Identifier )
    class(ClassBasis), intent(in)  :: self
    integer          , intent(in)  :: L
    character(len=*), optional, intent(in) :: Identifier 
    !
    select type ( ptr => Self.RadialBasis )
       class is ( ClassBSpline )
          NPreconditioned = self.NRegular(L)
       class is ( ClassGaussian )
          NPreconditioned = ptr.NPreconditioned( L )
       class is ( ClassGaussianBSpline )
             NPreconditioned = ptr.NPreconditioned( L, Identifier )
       class DEFAULT
    end select

  end function NPreconditioned


  !> Retrieves the vector of regular indexes for a given angular momentum
  integer function NPreconditionedContinuous( self, L, Identifier )
    class(ClassBasis), intent(in)  :: self
    integer          , intent(in)  :: L
    character(len=*), optional, intent(in) :: Identifier 
    !
    select type ( ptr => Self.RadialBasis )
       class is ( ClassBSpline )
          !.. The last B-spline was already eliminated when NRegular(L) 
          !   was called, but here NPreconditionedContinuous is lowered down one unit again. 
          NPreconditionedContinuous = self.NPreconditioned(L) - 1
       class is ( ClassGaussian )
          NPreconditionedContinuous = self.NPreconditioned( L )
       class is ( ClassGaussianBSpline )
          if ( (present(Identifier)) .and. (Identifier .is. "G") ) then
             !
             NPreconditionedContinuous = self.NPreconditioned( L,Identifier )
             !
          else
             !
             NPreconditionedContinuous = self.NPreconditioned( L ) - 1
             !
          end if
       class DEFAULT
    end select

  end function NPreconditionedContinuous


  !> Returns a linear grid on the radial domain appropriate for the current radial basis.
  subroutine GetPlotGrid( Basis, Xgrid )
    !
    class(ClassBasis),            intent(in)  :: Basis
    DoublePrecision, allocatable, intent(out) :: Xgrid(:)
    !
    DoublePrecision :: Xmin, Xmax
    integer         :: ix, npts
    !
    Xmin = Basis.RPlotMin
    Xmax = Basis.RPlotMax
    !
    npts = Basis.NPlot
    !
    if(allocated( XGRID )) deallocate( XGRID )
    allocate( XGRID( npts ) )
    do ix = 1, npts
       XGRID(ix) = Xmin+(Xmax-Xmin)*dble(ix-1)/dble(npts-1)
    end do
    !
  end subroutine GetPlotGrid


  !> Gets the maximum angular momentum in the Gaussian expansion.
  Integer function GetGaussianLMax( Basis ) result (Lmax)
    class(ClassBasis), intent(in)  :: Basis
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBSpline)
       call ErrorMessage( &
            "There are not Gaussian functions in this radial "//&
            "basis, only B-splines, a dummy value of 1 will be "//&
            "passed as Lmax.")
       Lmax = 1
    class is (ClassGaussianBSpline)
       Lmax = ptr.GetLMax()
    class is (ClassGaussian)
       Lmax = ptr.GetLMax()
    class DEFAULT
       call Assert("Invalid radial basis type")
    end select
  end function GetGaussianLMax


  !> Gets the maximum angular momentum in the Gaussian expansion.
  Integer function GetGaussianNumExponents( Basis ) result (NumExp)
    class(ClassBasis), intent(in)  :: Basis
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBSpline)
       call Assert("There are not Gaussian functions in this radial basis, only B-splines")
    class is (ClassGaussianBSpline)
       NumExp = ptr.GetNumExponents()
    class is (ClassGaussian)
       NumExp = ptr.GetNAlpha()
    class DEFAULT
       call Assert("Invalid radial basis type")
    end select
  end function GetGaussianNumExponents



  !> Get maximum radius for the current basis.
  !! Coincides with the maximum plot radius.
  DoublePrecision function GetRMax( Basis ) result (R)
    class(ClassBasis), intent(in)  :: Basis
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBSpline)
       R = ptr.GetNodePosition( ptr.GetNNodes() )
    class is (ClassGaussianBSpline)
       R = ptr.GetMaxPlotRadius()
    class is (ClassGaussian)
       R = ptr.GetMaxPlotRadius()
    class DEFAULT
       call ErrorMessage("Node position is not available for Gaussian Basis.")
       call ErrorMessage("You are possibly using a basis unsuitable to the current calculation.")
       call Assert("This error should not take place at this level, but at the input parsing level!")
    end select
  end function GetRMax



  !> Get a given node position for B-splines.
  DoublePrecision function GetNodePosition( Basis, Node ) result (R)
    class(ClassBasis), intent(in)  :: Basis
    integer,           intent(in)  :: Node
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBSpline)
       R = ptr.GetNodePosition( Node )
    class is (ClassGaussianBSpline)
       R = ptr.GetNodePosition( Node )
    class is (ClassGaussian)
       call Assert("Node position is not available for Gaussian Basis.")
    class DEFAULT
       call ErrorMessage("Node position is not available for Gaussian Basis.")
       call ErrorMessage("You are possibly using a basis unsuitable to the current calculation.")
       call Assert("This error should not take place at this level, but at the input parsing level!")
    end select
  end function GetNodePosition



  logical function PreconditionerIsAvailable( Basis, L ) result( Available )
    class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: L
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBspline)
       Available = ptr.PreconditionerIsAvailable( L, Basis.GetDir() )
    class is (ClassGaussian)
       Available = ptr.PreconditionerIsAvailable( L, Basis.GetDir() )
    class is (ClassGaussianBSpline)
       Available = ptr.PreconditionerIsAvailable( L, Basis.GetDir() )
    class DEFAULT
    end select
  end function PreconditionerIsAvailable


  subroutine ComputePreconditioner( Basis, L )
    class(ClassBasis),              intent(inout) :: Basis
    integer          ,              intent(in)    :: L
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBspline)
       call ptr.ComputePreconditioner( L )
    class is (ClassGaussian)
       call ptr.ComputePreconditioner( L )
    class is (ClassGaussianBSpline)
       call ptr.ComputePreconditioner( L )
    class DEFAULT
    end select
  end subroutine ComputePreconditioner


  subroutine SavePreconditioner( Basis, L )
    class(ClassBasis), intent(in) :: Basis
    integer          , intent(in) :: L
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBspline)
       call ptr.SavePreconditioner( L, Basis.GetDir() )
    class is (ClassGaussian)
       call ptr.SavePreconditioner( L, Basis.GetDir() )
    class is (ClassGaussianBSpline)
       call ptr.SavePreconditioner( L, Basis.GetDir() )
    class DEFAULT
    end select
  end subroutine SavePreconditioner




  subroutine LoadPreconditioner( Basis, L )
    class(ClassBasis), intent(inout) :: Basis
    integer          , intent(in)    :: L
    !
    select type( ptr => Basis.RadialBasis )
    class is (ClassBspline)
       call ptr.LoadPreconditioner( L, Basis.GetDir() )
    class is (ClassGaussian)
       call ptr.LoadPreconditioner( L, Basis.GetDir() )
    class is (ClassGaussianBSpline)
       call ptr.LoadPreconditioner( L, Basis.GetDir() )
    class DEFAULT
    end select
  end subroutine LoadPreconditioner


  subroutine CastMatrixBra_FromMatrix( self, MatIn, MatOut, InLevel, OutLevel, L )
    class(ClassBasis), intent(in)    :: self
    real(kind(1d0))  , intent(inout) :: MatIn(:,:)
    real(kind(1d0)), allocatable , intent(out) :: MatOut(:,:)
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    integer          , intent(in)    :: L
    !
    type(ClassMatrix) :: CMat
    CMat=MatIn
    call self.CastMatrixBra( CMat, InLevel, OutLevel, L )
    allocate(MatOut,source=CMat.A)
  end subroutine CastMatrixBra_FromMatrix

  subroutine CastMatrixBra_FromVector( self, VecIn, VecOut, InLevel, OutLevel, L )
    class(ClassBasis), intent(in)    :: self
    real(kind(1d0))  , intent(inout) :: VecIn(:)
    real(kind(1d0)), allocatable, intent(out) :: VecOut(:)
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    integer          , intent(in)    :: L
    !
    type(ClassMatrix) :: CMat
    CMat=VecIn
    call self.CastMatrixBra( CMat, InLevel, OutLevel, L )
    allocate(VecOut,source=CMat.A(:,1))
  end subroutine CastMatrixBra_FromVector

  subroutine CastMatrix( self, Mat, LBra, LKet, InLevel, OutLevel )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: LBra
    integer          , intent(in)    :: LKet
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    call self.CastMatrixBra( Mat, InLevel, OutLevel, LBra )
    call self.CastMatrixKet( Mat, InLevel, OutLevel, LKet )
  end subroutine CastMatrix


  subroutine CastMatrixBra_FromClassMatrix( self, Mat, InLevel, OutLevel, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    integer          , intent(in)    :: L
    select case( OutLevel )
    case( REP_LEVEL_TOTAL )
       call self.CastMatrixBraToTotal( Mat, InLevel, L )
    case( REP_LEVEL_REGULAR )
       call self.CastMatrixBraToRegular( Mat, InLevel, L )
    case( REP_LEVEL_PRECONDITIONED )
       call self.CastMatrixBraToPreconditioned( Mat, InLevel, L )
    case( REP_LEVEL_PRECONDITIONED_CONTINUOUS )
       call self.CastMatrixBraToPreconditionedContinuous( Mat, InLevel, L )
    end select
  end subroutine CastMatrixBra_FromClassMatrix


  subroutine CastMatrixKet_FromClassMatrix( self, Mat, InLevel, OutLevel, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    integer          , intent(in)    :: L
    select case( OutLevel )
    case( REP_LEVEL_TOTAL )
       call self.CastMatrixKetToTotal( Mat, InLevel, L )
    case( REP_LEVEL_REGULAR )
       call self.CastMatrixKetToRegular( Mat, InLevel, L )
    case( REP_LEVEL_PRECONDITIONED )
       call self.CastMatrixKetToPreconditioned( Mat, InLevel, L )
    case( REP_LEVEL_PRECONDITIONED_CONTINUOUS )
       call self.CastMatrixKetToPreconditionedContinuous( Mat, InLevel, L )
    end select
  end subroutine CastMatrixKet_FromClassMatrix

  subroutine CastMatrixKet_FromMatrix( self, MatIn, MatOut, InLevel, OutLevel, L )
    class(ClassBasis), intent(in)    :: self
    real(kind(1d0))  , intent(inout) :: MatIn(:,:)
    real(kind(1d0)), allocatable, intent(out) :: MatOut(:,:)
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    integer          , intent(in)    :: L
    !
    type(ClassMatrix) :: CMat
    CMat=MatIn
    call self.CastMatrixKet( CMat, InLevel, OutLevel, L )
    allocate(MatOut,source=CMat.A)
  end subroutine CastMatrixKet_FromMatrix

  subroutine CastMatrixKet_FromVector( self, VecIn, VecOut, InLevel, OutLevel, L )
    class(ClassBasis), intent(in)    :: self
    real(kind(1d0))  , intent(inout) :: VecIn(:)
    real(kind(1d0)), allocatable, intent(out) :: VecOut(:)
    integer          , intent(in)    :: InLevel
    integer          , intent(in)    :: OutLevel
    integer          , intent(in)    :: L
    !
    type(ClassMatrix) :: CMat
    call CMat.InitFull(1,UBOUND(VecIn,1)-LBOUND(VecIn,1)+1)
    CMat.A(1,:)=VecIn
    call self.CastMatrixKet( CMat, InLevel, OutLevel, L )
    allocate(VecOut,source=CMat.A(1,:))
  end subroutine CastMatrixKet_FromVector


  subroutine CastMatrixBraToTotal( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L
    
    character(len=*), parameter :: HERE="ClassBasis::CastMatrixBraToTotal : "
    integer :: NTotal, MinReg, MaxReg, NAdd
    character :: NewLine
    character(len=1000) :: ErrMsg
    NewLine = New_Line("a")

    if( Level == REP_LEVEL_TOTAL   ) return
    if( Level  > REP_LEVEL_REGULAR ) call self.CastMatrixBraToRegular( Mat, Level, L )

    Ntotal = self.GetNFun()
    if( Mat.NRows() /= self.NRegular(L) )then
       write(ErrMsg,*)&
            "Inconsistent matrix size "  , NewLine, &
            "L      = ", L               , NewLine, &
            "NReg   = ", self.NRegular(L), NewLine, &
            "Ntotal = ", Ntotal          , NewLine, &
            "NRows  = ", Mat.NRows()
       call Assert(HERE//trim(ErrMsg))
    endif

    !.. In the case of either Gaussian or B-spline, the total
    !   representation is obtained by just adding a few empty
    !   rows at the beginning and at the end of the matrix
    !.. 
    select type( ptr => self.RadialBasis )
    class is (ClassBSpline)

       MinReg = MinimumRegularBSplineIndex( ptr, L )
       MaxReg = MaximumRegularBSplineIndex( ptr )
       call Mat.AddRows(MinReg-1,"START")

    class is(ClassGaussian)

       MinReg = ptr.GetMinimumRegularIndex( L )
       MaxReg = ptr.GetMaximumRegularIndex( L )
       call Mat.AddRows(MinReg-1,"START")
       call Mat.AddRows(Ntotal-MaxReg,"END")

    class is(ClassGaussianBSpline)

       MinReg = ptr.GetMinimumGaussianRegularIndex( L )
       call Mat.AddRows(MinReg-1,"START")

       MaxReg = ptr.GetMaximumGaussianRegularIndex( L )
       NAdd = ptr.GetGaussianNFun() - MaxReg + ptr.GetNBSplineSkip()
       call Mat.AddRows( Nadd, After = MaxReg )

    class DEFAULT
    end select
  end subroutine CastMatrixBraToTotal


  subroutine CastMatrixBraToRegular( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L

    real(kind(1d0)), allocatable :: Precon(:,:)
    integer :: MinReg, MaxReg, NTotal, NRem
    integer :: NRegular, NPrecon
    
    Ntotal = self.GetNFun()

    if( Level < REP_LEVEL_REGULAR )then

       !.. Downstram conversion
       !..
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)

          MinReg = MinimumRegularBSplineIndex( ptr, L )
          MaxReg = MaximumRegularBSplineIndex( ptr )
          call Mat.RemoveRows(MinReg-1,"START")

       class is(ClassGaussian)

          MinReg = ptr.GetMinimumRegularIndex( L )
          MaxReg = ptr.GetMaximumRegularIndex( L )
          call Mat.RemoveRows(MinReg-1,"START")
          call Mat.RemoveRows(Ntotal-MaxReg,"END")

       class is(ClassGaussianBSpline)

          MinReg = ptr.GetMinimumGaussianRegularIndex( L )
          MaxReg = ptr.GetMaximumGaussianRegularIndex( L )
          NRem   = ptr.GetGaussianNFun() - MaxReg + ptr.GetNBSplineSkip()
          !.. Beware that the two removals don't commute
          call Mat.RemoveRows(NRem,After=MaxReg)
          call Mat.RemoveRows(MinReg-1,"START")

       class DEFAULT
       end select

    elseif( Level > REP_LEVEL_REGULAR )then

       !.. Upstream conversion
       !..
       call self.CastMatrixBraToPreconditioned( Mat, Level, L )
       !.. Convert up to regular
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)
          !
          !.. Do nothing: for B-spline, preconditioned = regular already
          !
       class is(ClassGaussian)
          !
          NRegular = ptr.NRegular( L )
          NPrecon  = ptr.NPreconditioned( L )
          allocate(Precon(NRegular,NPrecon))
          Precon=0.d0
          call ptr.FetchPreconditioner( Precon, L )
          !
          call Mat.Multiply(Precon,"Left","N")
          !
       class is(ClassGaussianBSpline)
          !
          !.. At the moment, returns a preconditioning
          !   matrix which is block-diagonal and 
          !   which coincides with the identity for 
          !   the non-overlapping BSplines.
          !   
          call ptr.FetchPreconditioner( Precon, L )
          call Mat.Multiply(Precon,"Left","N")
          !
       class DEFAULT
       end select

    end if

  end subroutine CastMatrixBraToRegular


  subroutine CastMatrixBraToPreconditioned( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L

    integer                      :: NRegular, NPrecon
    real(kind(1d0)), allocatable :: Precon(:,:)

    if( Level < REP_LEVEL_PRECONDITIONED )then

       !.. Down-stream conversion
       !..
       call self.castMatrixBraToRegular( Mat, Level, L )
       !
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)
          !
          !.. Do nothing: for B-splines, regular = preconditioned
          !
       class is(ClassGaussian)

          NRegular = ptr.NRegular( L )
          NPrecon  = ptr.NPreconditioned( L )
          allocate(Precon(NRegular,NPrecon))
          Precon=0.d0
          call ptr.FetchPreconditioner( Precon, L )
          call Mat.Multiply(Precon,"Left","T")

       class is(ClassGaussianBSpline)
          !
          !.. At the moment, returns a preconditioning
          !   matrix which is block-diagonal and 
          !   which coincides with the identity for 
          !   the non-overlapping BSplines.
          !   
          call ptr.FetchPreconditioner( Precon, L )
          call Mat.Multiply(Precon,"Left","T")
          !
       class DEFAULT
       end select

    elseif( Level >  REP_LEVEL_PRECONDITIONED )then
       
       !.. Up-stream conversion
       !..
       call self.castMatrixBraToPreconditionedContinuous( Mat, Level, L )
       !
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)

          call Mat.AddRows(1,"END")

       class is(ClassGaussian)
          !
          !.. Do nothing: all Gaussian functions are already continuous
          !
       class is(ClassGaussianBSpline)
          
          call Mat.AddRows(1,"END")
          
       class DEFAULT
       end select

    endif

  end subroutine CastMatrixBraToPreconditioned


  subroutine CastMatrixBraToPreconditionedContinuous( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L

    character(len=*), parameter :: HERE="ClassBasis::CastMatrixBraToPreconditionedContinuous"

    if( Level < REP_LEVEL_PRECONDITIONED_CONTINUOUS )then

       !.. Down-stream conversion
       !..
       call self.castMatrixBraToPreconditioned( Mat, Level, L )
       !
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)
          
          call Mat.RemoveRows(1,"END")
          
       class is(ClassGaussian)
          !
          !.. Do nothing: all Gaussian are already Continuous
          !
       class is(ClassGaussianBSpline)
          
          call Mat.RemoveRows(1,"END")
          
       class DEFAULT
       end select
       
    elseif( Level >  REP_LEVEL_PRECONDITIONED_CONTINUOUS )then
       !
       !.. There currently are no implemented levels
       !   downstram to preconditioned continuous
       !
       call Assert(HERE//"Upconversion *to* PRECONDITIONED_CONTINUOUS level not implemented")
       !
    endif
    
  end subroutine CastMatrixBraToPreconditionedContinuous


  subroutine CastMatrixKetToTotal( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L

    character(len=*), parameter :: HERE="ClassBasis::CastMatrixKetToTotal : "
    integer :: NTotal, MinReg, MaxReg, NAdd
    !
    if( Level == REP_LEVEL_TOTAL   ) return
    if( Level  > REP_LEVEL_REGULAR ) call self.CastMatrixKetToRegular( Mat, Level, L )

    Ntotal = self.GetNFun()


    !.. In the case of either Gaussian or B-spline, the total
    !   representation is obtained by just adding a few empty
    !   rows at the beginning and at the end of the matrix
    !.. 
    select type( ptr => self.RadialBasis )
    class is (ClassBSpline)

       MinReg = MinimumRegularBSplineIndex( ptr, L )
       MaxReg = MaximumRegularBSplineIndex( ptr )
       call Mat.AddColumns(MinReg-1,"START")

    class is(ClassGaussian)

       MinReg = ptr.GetMinimumRegularIndex( L )
       MaxReg = ptr.GetMaximumRegularIndex( L )
       call Mat.AddColumns(MinReg-1,"START")
       call Mat.AddColumns(Ntotal-MaxReg,"END")

    class is(ClassGaussianBSpline)

       MinReg = ptr.GetMinimumGaussianRegularIndex( L )
       call Mat.AddColumns(MinReg-1,"START")

       MaxReg = ptr.GetMaximumGaussianRegularIndex( L )
       NAdd = ptr.GetGaussianNFun() - MaxReg + ptr.GetNBSplineSkip()
       call Mat.AddColumns( Nadd, After = MaxReg )

    class DEFAULT
    end select

  end subroutine CastMatrixKetToTotal


  subroutine CastMatrixKetToRegular( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L

    real(kind(1d0)), allocatable :: Precon(:,:)
    integer :: MinReg, MaxReg, NTotal, NRem
    integer :: NRegular, NPrecon
    
    Ntotal = self.GetNFun()

    if( Level < REP_LEVEL_REGULAR )then
       

       !.. Convert downstream
       !..
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)

          MinReg = MinimumRegularBSplineIndex( ptr, L )
          MaxReg = MaximumRegularBSplineIndex( ptr )
          call Mat.RemoveColumns(MinReg-1,"START")

       class is(ClassGaussian)

          MinReg = ptr.GetMinimumRegularIndex( L )
          MaxReg = ptr.GetMaximumRegularIndex( L )
          call Mat.RemoveColumns(MinReg-1,"START")
          call Mat.RemoveColumns(Ntotal-MaxReg,"END")

       class is(ClassGaussianBSpline)

          MinReg = ptr.GetMinimumGaussianRegularIndex( L )
          MaxReg = ptr.GetMaximumGaussianRegularIndex( L )
          NRem   = ptr.GetGaussianNFun() - MaxReg + ptr.GetNBSplineSkip()
          !.. Beware that the two removals don't commute
          call Mat.RemoveColumns(NRem,After=MaxReg)
          call Mat.RemoveColumns(MinReg-1,"START")

       class DEFAULT
       end select

    elseif( Level > REP_LEVEL_REGULAR )then

       !.. Convert upstream
       !..
       call self.CastMatrixKetToPreconditioned( Mat, Level, L )
       !.. Convert up to regular
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)
          !
          !.. Do nothing: for B-spline, preconditioned = regular already
          !
       class is(ClassGaussian)
          !
          NRegular = ptr.NRegular( L )
          NPrecon  = ptr.NPreconditioned( L )
          allocate(Precon(NRegular,NPrecon))
          Precon=0.d0
          call ptr.FetchPreconditioner( Precon, L )
          !
          call Mat.Multiply(Precon,"Right","T")
          !
       class is(ClassGaussianBSpline)
          !
          call ptr.FetchPreconditioner( Precon, L )
          call Mat.Multiply(Precon,"Right","T")
          !
       class DEFAULT
       end select

    end if

  end subroutine CastMatrixKetToRegular


  subroutine CastMatrixKetToPreconditioned( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L
    
    integer                      :: NRegular, NPrecon
    real(kind(1d0)), allocatable :: Precon(:,:)

    if( Level < REP_LEVEL_PRECONDITIONED )then

       !.. Down-stream conversion
       !..
       call self.castMatrixKetToRegular( Mat, Level, L )
       !
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)
          !
          !.. Do nothing: for B-splines, regular = preconditioned
          !
       class is(ClassGaussian)

          NRegular = ptr.NRegular( L )
          NPrecon  = ptr.NPreconditioned( L )
          allocate(Precon(NRegular,NPrecon))
          Precon=0.d0
          call ptr.FetchPreconditioner( Precon, L )
          call Mat.Multiply(Precon,"Right","N")
          
       class is(ClassGaussianBSpline)

          !.. At the moment, returns a preconditioning
          !   matrix which is block-diagonal and 
          !   which coincides with the identity for 
          !   the non-overlapping BSplines.
          !   
          call ptr.FetchPreconditioner( Precon, L )
          call Mat.Multiply(Precon,"Right","N")

       class DEFAULT
       end select

    elseif( Level >  REP_LEVEL_PRECONDITIONED )then
       
       !.. Up-stream conversion
       !..
       call self.castMatrixKetToPreconditionedContinuous( Mat, Level, L )
       !
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)

          call Mat.AddColumns(1,"END")

       class is(ClassGaussian)
          !
          !.. Do nothing: all Gaussian functions are already continuous
          !
       class is(ClassGaussianBSpline)
          
          call Mat.AddColumns(1,"END")
          
       class DEFAULT
       end select

    endif

    
  end subroutine CastMatrixKetToPreconditioned


  subroutine CastMatrixKetToPreconditionedContinuous( self, Mat, Level, L )
    class(ClassBasis), intent(in)    :: self
    type(ClassMatrix), intent(inout) :: Mat
    integer          , intent(in)    :: Level
    integer          , intent(in)    :: L
    
    character(len=*), parameter :: HERE="ClassBasis::CastMatrixKetToPreconditionedContinuous"

    if( Level < REP_LEVEL_PRECONDITIONED_CONTINUOUS )then

       !.. Down-stream conversion
       !..
       call self.castMatrixKetToPreconditioned( Mat, Level, L )
       !
       select type( ptr => self.RadialBasis )
       class is (ClassBSpline)
          
          call Mat.RemoveColumns(1,"END")
          
       class is(ClassGaussian)
          !
          !.. Do nothing: all Gaussian are already Continuous
          !
       class is(ClassGaussianBSpline)
          call Mat.RemoveColumns(1,"END")
       class DEFAULT
       end select
       
    elseif( Level >  REP_LEVEL_PRECONDITIONED_CONTINUOUS )then
       !
       !.. There currently are no implemented levels
       !   downstram to preconditioned continuous
       !
       call Assert(HERE//"Upconversion *to* PRECONDITIONED_CONTINUOUS level not implemented")
       !
    endif

  end subroutine CastMatrixKetToPreconditionedContinuous


  !> Gets the radial distance at which the plot will start.
  real(kind(1d0)) function GetRPlotMin( Basis ) result (Rmin)
    !
    class(ClassBasis), intent(in)  :: Basis
    !
    Rmin = Basis.RPlotMin
    !
  end function GetRPlotMin


  !> Gets the radial distance at which the plot will finish.
  real(kind(1d0)) function GetRPlotMax( Basis ) result (Rmax)
    !
    class(ClassBasis), intent(in)  :: Basis
    !
    Rmax = Basis.RPlotMax
    !
  end function GetRPlotMax



  !> Gets a vector with all the monomials exponents of the Gaussian functions.
  subroutine GetAllMonomials( self, Vec )
    !
    class(ClassBasis),    intent(in)  :: self
    integer, allocatable, intent(out) :: Vec(:)
    !
    select type( ptr => self.RadialBasis )
    class is (ClassBSpline)
       call Assert( "The stand-alone B-splines radial basis does not have monomials, Gaussian functions do have." )
    class is(ClassGaussian)
       call ptr.GetAllMonomials( Vec )
    class is(ClassGaussianBSpline)
       call ptr.GetAllMonomials( Vec )
    class DEFAULT
       call Assert( "Invalid radial basis type." )
    end select
    !
  end subroutine GetAllMonomials


  subroutine RemoveGaussianRowNormalization( Basis, GABSMat )
    !
    class(ClassBasis) , intent(in)    :: Basis
    class(ClassMatrix), intent(inout) :: GABSMat
    !
    integer :: i, j, NumGaussianGABS, NumFunGABS
    real(kind(1d0)) :: GaussNorm, NewElemVal
    !
    NumFunGABS      = Basis.GetNFun()
    if ( GABSMat.NRows() /= NumFunGABS ) call Assert( "To remove the Gaussian normalization the full GABS matrix must be passed." )
    if ( GABSMat.NColumns() /= NumFunGABS ) call Assert( "To remove the Gaussian normalization the full GABS matrix must be passed." )

    NumGaussianGABS = Basis.GetNFun("G")
    !
    do i = 1, NumGaussianGABS
       !
       GaussNorm = Basis.GetNormFactor( i )
       !
       do j = 1, NumFunGABS
          NewElemVal = GABSMat.Element(i,j) / GaussNorm
          call GABSMat.SetElement(i,j,NewElemVal)
          !
       end do
    end do
    !
  end subroutine RemoveGaussianRowNormalization




  !> Retrieves the minimum and maximum index corresponding to a given angular momentum.
  subroutine GetLRowInterval( self, L, Nmin, Nmax )
    !
    class(ClassBasis),    intent(in)  :: self
    integer, intent(in)  :: L
    integer, intent(out) :: Nmin
    integer, intent(out) :: Nmax
    !
    integer :: i
    logical :: SetNmin
    !
    integer, allocatable  :: Vec(:)

    call self.GetAllMonomials(Vec)
    SetNmin = .false.
    do i = 1, size(Vec)
       if ( (Vec(i) >= L) .and. (mod(Vec(i)-L,2) == 0) ) then
          if ( .not. SetNmin ) then
             Nmin = i
             Nmax = Nmin
             SetNmin = .true.
          else
             Nmax = Nmax + 1
          end if
       end if
    end do
    !
  end subroutine GetLRowInterval




     !> Gets a vector with all the monomials exponents of the Gaussian functions.
  subroutine GetAllExponents( self, Vec )
    !
    class(ClassBasis),    intent(in)  :: self
    real(kind(1d0)), allocatable, intent(out) :: Vec(:)
    !
    select type( ptr => self.RadialBasis )
    class is (ClassBSpline)
       call Assert( "The stand-alone B-splines radial basis does not have exponents, Gaussian functions do have." )
    class is(ClassGaussian)
       call ptr.GetAllExponents( Vec )
    class is(ClassGaussianBSpline)
       call ptr.GetAllExponents( Vec )
    class DEFAULT
       call Assert( "Invalid radial basis type." )
    end select
    !
  end subroutine GetAllExponents



end module ModuleBasis
