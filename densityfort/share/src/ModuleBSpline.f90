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

! {{{ Detailed description

!> \file
!!
!! Defines the variables and methods for B-splines basis.

! }}}
module ModuleBSpline

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleParameterList
  use ModuleString
  use ModuleSystemUtils

  implicit none

  private


  ! {{{ Private Attributes

  !
  logical                    :: INITIALIZED = .FALSE.
  !
  !.. Local Gaussian weights
  DoublePrecision, parameter :: PI = 3.14159265358979323846d0
  DoublePrecision, parameter :: GAUSS_POINTS_CONVERGENCE_THRESHOLD = 2.d-15
  integer        , parameter :: NGAUSS = 64
  DoublePrecision            :: Gauss_points(NGAUSS)
  DoublePrecision            :: Gauss_weight(NGAUSS)
  !
  !.. Local Factorial Parameters and variables
  integer, parameter :: DIMFACT = 127
  DoublePrecision    :: fact(0:DIMFACT)
  DoublePrecision    :: tar(0:((DIMFACT+1)*(DIMFACT+2)/2-1))
  DoublePrecision    :: parfactMat(0:DIMFACT,0:DIMFACT)
  !
  !.. B-spline parameters
  integer        , parameter :: MAX_NUMBER_NODES = 10000
  DoublePrecision, parameter :: NOD_THRESHOLD    = 1.d-15
  integer        , parameter :: MAXNR            = 20
  !

  ! }}}

  !> Set of B-splines (spline basis functions)
  !!
  ! {{{ Detailed Description

  !> The spline are defined as piecewise polynomials of order \f$k\f$ 
  !! (maximum degree \f$k-1\f$), \f$\mathcal{C}^\infty\f$ everywhere, except at a fixed set of knots \f$\{t_i\}\f$ 
  !! where they are at least \f$\mathcal{C}^{k_i-2}\f$, with \f$k \geq k_i\geq 1\f$. 
  !! \f$\nu_i=k-k_i+1\f$ is called the knot multiplicity because the same spline basis can be obtained 
  !! from a basis with \f$\sum_i \nu_i\f$ maximally regular knots in the limit where sets of \f$\nu_i\f$ 
  !! contiguous knots coalesce to \f$i\f$-th knots. In other terms, the spline space is specified by an order 
  !! \f$k\f$ and a non-decreasing set of knots. The total dimension of spline space is \f$n+k\f$, but if 
  !! we specialize to the subset which is zero below the lowest and above the highest knots we are left 
  !! with a space of dimension \f$n-k\f$. As a consequence every spline extends at least over \f$k\f$ adjacent 
  !! intervals (\f$k+1\f$ consecutive knots, when counted with their multiplicity). The \f$n-k\f$ independent 
  !! splines restricted to \f$k\f$ intervals, clearly a basis of the spline space, are commonly called B-splines.
  !! Following deBoor\cite{deBoor}, we define the \f$i\f$-th B-spline \f$B_i^k(x)\f$ of order \f$k\f$ and which 
  !! extends from the knot \f$t_i\f$ to the knot \f$t_{i+k}\f$, as follows:
  !! \f{eqnarray}
  !! B_i^1(x)&=&\theta(x-t_i)\,\cdot\,\theta(t_{i+1}-x)\\
  !! B_i^k(x)&=&\frac{x-t_i}{t_{i+k-1}-t_i}\,\,
  !!            B_i^{k-1}(x)+\frac{t_{i+k}-x}{t_{i+k}-t_{i+1}}\,\,
  !!            B_{i+1}^{k-1}(x).
  !! \f}
  !! In the following, unless otherwise stated, we shall refer to standard set of knots where the 
  !! first and last breakpoints are \f$k\f$ times degenerate, while the other nodes are non degenerate:
  !! \f[
  !! t_1=t_2=\ldots=t_k\leq t_{k+1}\leq\ldots\leq t_{n-k}\leq t_{n-k+1}=t_{n-k+2}=\ldots=t_{n}
  !! \f]
  !! The use of B-splines in atomic and molecular physics calculations is reviewed in 
  !! [\cite{rpp.64.1815}](http://iopscience.iop.org/0034-4885/64/12/205/). 
  !! 
  !! B-splines are invariant upon affine transformations: that is, if
  !! \f{equation}
  !! \mathrm{if}\quad x\to x'=a\,x +b,\quad t_i\to t_i'=a\,t_i+b\quad\mathrm{then}
  !! \quad{B_i^{k}}'(x')=B_i^k(x).
  !! \f}
  !! 
  !! It is useful to define the \f$L^2\f$-normalized B-splines as
  !! \f{equation}
  !! \bar{B}_i(x)\equiv B_i(x)/\|B_i\|_{L^2}
  !! \f}
  !! 
  !! Finally, the fact that each couple of different B-splines has a 
  !! significant portion of their domain where they do not superpose
  !! gives rise to a strong linear independence. 

  ! }}} 
  type, public :: ClassBSpline
     !
     private
     !> Number of non-degenerate nodes.
     integer :: NNodes
     !
     !> Spline order.
     !! A B-spline of order \f$ k\f$ is a piecewise polynomials
     !! of order \f$k-1\f$ which extends over \f$k\f$ intervals (including
     !! degenerate intervals, i.e., those with zero length).
     !! The quality of the Bspline evaluation is severely compromised
     !! for high order. At order around 35, in the case of just two nodes
     !! (which is the most demanding condition), the numerical error starts
     !! being of the same order of magnitude of the function itself.
     !! For order \f$\sim\f$20, there are occasional points of irregularity
     !! for a reasonable number of intervals (e.g.,\f$\geq 10\f$). Anyhow,
     !! these caveats arise only in cases which already lie very much outside
     !! the range of applicability of the BSplines as a variational basis, 
     !! which is above order 12 or so.
     integer :: Order
     !
     !> Total number of B-spline. 
     !! \f[
     !! NBsplines = NNodes + Order - 2
     !!\f]
     integer :: NBSplines
     !
     !> grid which contains the node positions. 
     !! - \f$g_{-k+1}\,=\,g_{-k+2}\,=\,\cdots\,=\,g_{-1}\,=\,g_0\f$ 
     !!   is the first \f$ k\f$ times degenerate node.
     !! - \f$g_{n-1}\,=\,g_{n}\,=\,\cdots\,=\,g_{n+k-2}\f$ is the 
     !!   last \f$ k\f$ times degenerate node.
     !! The set of distinct nodes is mapped to the positions [0:n-1] of grid;  
     !! and hence, the upper bound of the i-th BSpline support is grid(i)
     real(kind(1d0)),  allocatable :: grid(:)
     !
     !> Normalization constants.
     !! \f[
     !!    \mathrm{f}_i=\|B_i\|_{L^2}^{-1}=\frac{1}{ \sqrt{ \langle B_i | B_i \rangle } }
     !! \f]
     real(kind(1d0)), private, allocatable :: f(:)
     !
     !
     !> Matrix of polynomial coefficients.
     !! The support $\mathcal{D}_i$ of the \f$i\f$-th B-spline comprises \f$k\f$ consecutive intervals:
     !! \f[
     !!      \mathcal{D}_i = [g_{i-k},g_i] = [g_{i-k+1},g_{i-k+1}] \cup \cdots \cup [g_{i-1}:g_i].
     !! \f]
     !! If we enumerate the individual intervals from \f$0\f$ to \f$k-1\f$, the explicit expression of 
     !! the \f$i\f$-th B-spline in interval $j$ is
     !! \f[
     !!     B_i(x) = \sum_{m=0}^{k-1}
     !!        \left(
     !!                \frac{x-g_{i-k+j}}{g_{i-k+j+1}-g_{i-k+j}}
     !!        \right)^m \,\,  c(m,j,i)
     !! \f]
     real(kind(1d0)), private, allocatable :: c(:,:,:)
     !
!!$     !> Node Characteristic functions
!!$     logical        , private              :: CARFUNINIT=.FALSE.
!!$     DoublePrecision, private, allocatable :: carfun(:,:)
!!$     DoublePrecision, private, allocatable :: carfunt(:,:)
!!$     DoublePrecision, private, allocatable :: vec(:)
     !
     !> Indicates whether ClassBSpline has been initialized or not.
     logical, private :: INITIALIZED=.FALSE.
     !
   contains

     !.. Basics
     !> Initializes ClassBSpline.
     generic   :: Init =>  &
          BSplineInit, &
          BSplineInitReadGridFromFile, &
          BSplineInitReadConfigurationFile
     !> Deallocates the ClassBSpline atributes.
     procedure :: Free =>  BSplineFree
     procedure :: Save =>  ClassBSpline_Save
     procedure :: Load =>  ClassBSpline_Load

!!$     generic   :: assignment(=) => ClassBspline_Copy
!!$     procedure :: ClassBspline_Copy 
     
     !.. Accessors
     !> Gets the number of nodes.
     procedure :: GetNNodes         =>  BSplineGetNNodes
     !> Gets the B-splines order.
     procedure :: GetOrder          =>  BSplineGetOrder
     !> Gets the number of B-splines.
     procedure :: GetNBSplines      =>  BSplineGetNBSplines
     !> Gets the position of a requested node.
     procedure :: GetNodePosition   =>  BSplineNodePosition
     !> Gets the normalization factor of a requested B-spline.
     procedure :: GetNormFactor     =>  BSplineNormalizationFactor

     !> Evaluates either a choosen B-spline function, 
     !! or a linear combination of B-splines, both evaluated in certain position. 
     generic   :: Eval              =>  BSplineEval, BSplineFunctionEval

     !> Tabulates in a slected domain either  a choosen B-spline function, 
     !! or a linear combination of B-splines.
     generic   :: Tabulate          =>  BSplineTabulate, BSplineFunctionTabulate

     !> Computes the integral
     !! \f[
     !!    \int_{a}^{b} r^{2}dr
     !!       \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}}
     !!\f]
     !! Where \f$f(r)\f$ is a local operator, and if a break point \f$BP\f$ is introduced, 
     !! then the integral is splitted in two parts
     !! \f[
     !!    \int_{a}^{BP} r^{2}dr
     !!       \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}} + 
     !!    \int_{BP}^{b} r^{2}dr
     !!       \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}}
     !! \f]
     generic   :: Integral =>  BSplineIntegral
     procedure :: BSplineIntegral

     !> Computes the integral
     !! \f[
     !!    \int_a^b dr 
     !!    \frac{d^{n_1}B_1(r)}{dr^{n_1}} f(r)
     !! \f]
     !! where \f$f(r)\f$ is a *smooth* function in the B-splines support.
     !! If the first (last) boundary of the integration interval is not 
     !! specified, then the value \f$a=-\infty\f$ and \f$b=\infty\f$ is assumed.
     procedure :: MonoIntegral      =>  BSplineMonoIntegral

     !> Computes the principal-part B-spline integral
     !! \f[
     !!  \int dr Bs_1(r)\,\frac{\mathcal{P}}{r-r_0}\,Bs_2(r)
     !! \f]
!!$     procedure :: CauchyIntegral    =>  BSplineCauchyIntegral

     !> Computes the monoelectronic or bielectronic multipole.
     generic   :: Multipole         =>  MonoelectronicMultipole, BielectronicMultipole

     !> Computes the Spherical Bessel Transform of a B-spline, defined as:
     !! \f[
     !!   F_{l}(p)=\int\frac{Bs_{n}(r)}{r}j_{l}(p\cdot r)r^{2}dr,
     !! \f]
     !! for all the angular momenta \f$\ell\f$ up to a specified lmax as 
     !! \f$ res(\ell+1)=F_\ell(p)\f$.
     procedure :: SphericalBesselTransform

     !> Logical function that retrieves whether the preconditioner is available or not.
     procedure :: PreconditionerIsAvailable
     !> Saves the preconditioner.
     procedure :: SavePreconditioner
     !> Loads the preconditioner.
     procedure :: LoadPreconditioner
     !> Computes the preconditioner.
     procedure :: ComputePreconditioner
     !> Creates in a directory all the configuration
     !!    files that define the current basis.
     procedure :: SaveConfig
     

     !.. Assignment
     !> Copies the B-spline set from one ClassBSpline to another. 
     generic, public :: ASSIGNMENT(=) => CopyBSplineSet

     !.. Private Interface
     procedure, private :: MonoElectronicMultipole
     procedure, private :: BiElectronicMultipole
     procedure, private :: BSplineInit
     procedure, private :: BsplineInitReadGridFromFile
     procedure, private :: BsplineInitReadConfigurationFile
     procedure, private :: CopyBSplineSet
     procedure, private :: BSplineEval
     procedure, private :: BSplineFunctionEval
     procedure, private :: BSplineTabulate
     procedure, private :: BSplineFunctionTabulate

  end type ClassBSpline



  public :: MixedBSplineIntegral

contains


  subroutine ClassBspline_Copy( self, bSet )
    class(ClassBspline), intent(out) :: self
    type (ClassBspline), intent(in)  :: bSet
    call self%Free()
    self%Nnodes    = bSet%Nnodes
    self%Order     = bSet%Order
    self%NBsplines = bSet%NBsplines
    allocate(self%grid, source=bSet%grid )
    allocate(self%f   , source=bSet%f )
    allocate(self%c   , source=bSet%c )
    self%INITIALIZED = .TRUE.
  end subroutine ClassBspline_Copy
   
  subroutine ClassBspline_Save( self, uid )
    class(ClassBspline), intent(in) :: self
    integer            , intent(in)  :: uid
    integer :: i
    write(uid,"(*(x,i0))") &
         self%Nnodes, self%Order, self%NBsplines, &
         size(self%grid), size(self%f), &
         (size(self%c,i),i=1,3) 
    write(uid,"(*(x,e24.16))") self%grid
    write(uid,"(*(x,e24.16))") self%f
    write(uid,"(*(x,e24.16))") self%c
  end subroutine ClassBspline_Save
   
  subroutine ClassBspline_Load( self, uid )
    class(ClassBspline), intent(out) :: self
    integer            , intent(in)  :: uid
    integer :: ng, nf, nc1, nc2, nc3
    call self%free()
    read(uid,*) self%Nnodes, self%Order, self%NBsplines, &
         ng, nf, nc1, nc2, nc3
    allocate(self%grid(ng))
    allocate(self%f   (nf))
    allocate(self%c(nc1,nc2,nc3))
    read(uid,*) self%grid
    read(uid,*) self%f
    read(uid,*) self%c
  end subroutine ClassBspline_Load


  !> Initialize module basic constants
  subroutine Init_Spline_Module()
    call InitGaussPoints
    call InitFactorials
    INITIALIZED=.TRUE.
  end subroutine Init_Spline_Module


  subroutine CheckInitialization(s)
    Class(ClassBSpline), intent(in) :: s
    if(.not.s%INITIALIZED) call Assert("Spline set is not initialized")
  end subroutine CheckInitialization


  !> Initialize the points and weights for 
  !! Gauss integration
  subroutine InitGaussPoints()
    integer         :: i,j
    DoublePrecision :: p1, p2, p3, pp, z, z1
    do i=1,(NGAUSS+1)/2
       z=cos(PI*(i-.25d0)/(NGAUSS+.5d0))
       inna : do
          p1=1.d0
          p2=0.d0
          do j=1,NGAUSS
             p3=p2
             p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
          enddo
          pp=NGAUSS*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
          if(abs(z-z1)<=GAUSS_POINTS_CONVERGENCE_THRESHOLD)exit inna
       enddo inna
       Gauss_Points(i)=(1.d0-z)/2.d0
       Gauss_Points(NGAUSS+1-i)=(1.d0+z)/2.d0
       Gauss_Weight(i)=1.d0/((1.d0-z*z)*pp*pp)
       Gauss_Weight(NGAUSS+1-i)=1.d0/((1.d0-z*z)*pp*pp)
    enddo
  end subroutine InitGaussPoints


  subroutine InitFactorials()
    integer :: i,j,k,n
    !
    !Initialize Factorials
    fact(0)=1.d0
    do i = 1, DIMFACT
       fact(i)=fact(i-1)*dble(i)
    enddo
    !
    !.. Initialize Binomials
    k=0
    tar=0.d0
    do i = 0,DIMFACT
       do j = 0,i
          tar(k)=fact(i)/fact(i-j)/fact(j)
          k=k+1
       enddo
    enddo
    !
    parfactmat=1.d0
    do n=1,DIMFACT
       do k=1,n
          do i=n-k+1,n
             parfactmat(n,k)=parfactmat(n,k)*dble(i)
          enddo
       enddo
    enddo
    !
  end subroutine InitFactorials

  ! }}}


  ! {{{ BSpline Accessor Methods

  integer function BSplineGetOrder(s) result(Order)
    class(ClassBSpline), intent(in) :: s
    Order=s%Order
  end function BSplineGetOrder

  integer function BSplineGetNNodes(s) result(NNodes)
    class(ClassBSpline), intent(in) :: s
    NNodes=s%NNodes
  end function BSplineGetNNodes

  integer function BSplineGetNBSplines(s) result(NBSplines)
    class(ClassBSpline), intent(in) :: s
    NBSplines=s%NBSplines
  end function BSplineGetNBSplines

  !> Return the position of a node
  DoublePrecision function BSplineNodePosition(SplineSet,Node) result(Position)
    Class(ClassBSpline), intent(in) :: SplineSet
    integer           , intent(in) :: Node
    if( Node <= 1 )then
       Position = SplineSet%Grid(0)
    elseif( Node >= SplineSet%NNodes )then
       Position = SplineSet%Grid( SplineSet%NNodes-1 )
    else
       Position = SplineSet%Grid( Node-1 )
    endif
    return
  end function BSplineNodePosition

  !> Return the position of a node
  DoublePrecision function BSplineNormalizationFactor(SplineSet,Bs) &
       result(Normalization)
    Class(ClassBSpline), intent(in) :: SplineSet
    integer            , intent(in) :: Bs
    Normalization = 0.d0
    if(Bs<1.or.Bs>SplineSet%NBSplines)return
    Normalization=SplineSet%f(Bs)
  end function BSplineNormalizationFactor

  ! }}}


  !> Copies the B-spline set sOrigin to sDestination 
  subroutine CopyBSplineSet(sDestination,sOrigin)
    Class(ClassBSpline), intent(inout):: sDestination
    Class(ClassBSpline), intent(in)   :: sOrigin
    !
    sDestination%NNodes=sOrigin%NNodes
    sDestination%Order=sOrigin%Order
    sDestination%NBSplines=sOrigin%NBSplines
    !
    if(allocated(sDestination%Grid))deallocate(sDestination%Grid)
    allocate(sDestination%Grid(1-sDestination%Order:sDestination%NBSplines))
    sDestination%Grid=sOrigin%Grid
    !
    if(allocated(sDestination%c))deallocate(sDestination%c)
    allocate(sDestination%c(&
         0:sDestination%Order-1 ,&
         0:sDestination%Order-1 ,&
         1:sDestination%NBSplines+sDestination%Order-2))
    sDestination%c=sOrigin%c
    !
    if(allocated(sDestination%f))deallocate(sDestination%f)
    allocate(sDestination%f(sDestination%NBSplines))
    sDestination%f=sOrigin%f
    !
    sDestination%INITIALIZED=.TRUE.
  end subroutine CopyBSplineSet


  !> Initialize the set of B-spline after reading the parameteres from a file
  subroutine BSplineInitReadConfigurationFile( self, ConfigurationFile, IOSTAT )
    !
    Class(ClassBSpline), intent(inout)         :: self
    character(len=*)   , intent(in)            :: ConfigurationFile
    integer            , intent(out), optional :: IOSTAT
    !
    type(ClassParameterList) :: List
    integer :: NNodes, Order, status, inod
    real(kind(1d0)), allocatable :: LocalGrid(:)
    real(kind(1d0)) :: Rmin, Rmax 
    character(len=512) :: NodeFile, Method
    !
    !.. Define parameters and parse configuration file
    call List%Add( "NumberNodes",            101, "required" )
    call List%Add( "SplineOrder",              7, "required" )
    call List%Add( "Method"     ,         "File", "required" )
    call List%Add( "NodeFile"   , "BSplineNodes", "optional" )
    call List%Add( "Rmin"       ,         0.d0  , "optional" )
    call List%Add( "Rmax"       ,       100.d0  , "optional" )
    call List%Parse( ConfigurationFile )
    call List%Get( "NumberNodes", Nnodes )
    call List%Get( "SplineOrder", Order )
    call List%Get( "NodeFile"   , NodeFile )
    call List%Get( "Method"     , Method )
    call List%Get( "Rmin"       , Rmin )
    call List%Get( "Rmax"       , Rmax )
    
    if( Nnodes < 2 )then
       call ErrorMessage(" NumberNodes parameter in "//trim(ConfigurationFile)//" must be greater than 2.")
       STOP
    endif
    if( Order < 4 )then
       call ErrorMessage(" SplineOrder parameter in "//trim(ConfigurationFile)//" is anomalously small.")
    endif
    if( trim(Method) .is. "File")then

       if( List%Present("NodeFile") )then 
          call self%Init( Nnodes, Order, NodeFile, iostat = status )
       else
          call ErrorMessage("NodeFile missing in "//trim(ConfigurationFile)//" while Method = File")
          STOP
       endif

    elseif( trim(Method) .is. "UniformGrid" )then

       if(Rmin<0.d0)call ErrorMessage("Negative Lower radial bound")
       if(Rmax<Rmin)then
          call ErrorMessage("Rmax parameter in "//trim(ConfigurationFile)//" must be larger than Rmin")
          STOP
       endif
       allocate(LocalGrid(NNodes))
       LocalGrid=[(Rmin+(Rmax-Rmin)*dble(inod-1)/dble(NNodes-1),inod=1,NNodes)]
       call self%Init(NNodes,Order,LocalGrid,status)  

    else

       call ErrorMessage(" '"//trim(Method)//"' is entry for the parameter Method in "//trim(ConfigurationFile))
       call ErrorMessage(" Valid Method Options: File, UniformGrid")
       STOP

    endif
    call List%Free( )

    if(present(IOSTAT))then
       IOSTAT=status
    elseif(status/=0)then
       call ErrorMessage("Spline Initialization failed")
    end if

  end subroutine BSplineInitReadConfigurationFile



  subroutine BSplineInitReadGridFromFile(s,NumberOfNodes,order,GridFile,IOSTAT)
    !
    Class(ClassBSpline), intent(inout)         :: s
    integer            , intent(in)            :: NumberOfNodes
    integer            , intent(in)            :: order
    character(len=*)   , intent(in)            :: GridFile
    integer            , intent(out), optional :: IOSTAT

    integer :: i, uid, Status, FileStatus, ReadStatus
    character(len=IOMSG_LENGTH)  :: iomsg
    DoublePrecision, allocatable :: LocalGrid(:)

    if(present(IOSTAT))IOSTAT=0

    if(.not.INITIALIZED) call Init_Spline_Module()

    !.. Check input
    call CheckNodesAndOrder(NumberOfNodes,Order,Status)
    if(Status/=0)then
       if(present(IOSTAT))then
          IOSTAT=Status
          return
       endif
       call Assert("Invalid Number of Nodes / Order")
    endif

    !.. Read the node grid from file
    open(NEWUNIT= uid,&
         FILE   = GridFile   ,&
         STATUS = "Old"      ,&
         FORM   = "Formatted",&
         IOSTAT = FileStatus ,&
         IOMSG  = iomsg)
    if(FileStatus/=0)then
       if(present(IOSTAT))then
          IOSTAT = -5
          return
       endif
       call Assert(iomsg)   
    endif
    !
    allocate(LocalGrid(NumberOfNodes))
    LocalGrid=0.d0
    !
    do i = 1, NumberOfNodes
       read(uid,*,&
            IOSTAT=ReadStatus,&
            IOMSG=iomsg) LocalGrid(i)
       if(ReadStatus/=0)exit
    enddo
    close(uid,IOSTAT=FileStatus)
    if(ReadStatus/=0)then
       deallocate(LocalGrid)
       if(present(IOSTAT))then
          IOSTAT = -5
          return
       else
          call Assert(iomsg)
       endif
    else
       if(FileStatus/=0) call Assert("Close error")
    endif


    call BSplineInit(s,NumberOfNodes,Order,LocalGrid,status)
    if(allocated(LocalGrid))deallocate(LocalGrid) 
    if(status/=0)then
       if(present(IOSTAT))then
          IOSTAT=status
          return
       else
          call Assert("Call to BSplineInit failed")
       end if
    endif

  end subroutine BSplineInitReadGridFromFile

  
  subroutine CheckNodesAndOrder(NumberOfNodes,Order,STAT)
    integer, intent(in)  :: NumberOfNodes, Order
    integer, intent(out) :: STAT
    STAT = -1; if ( NumberOfNodes < 2 ) return
    STAT = -2; if (    order      < 1 ) return
    STAT = 0
  end subroutine CheckNodesAndOrder
  

  !> Initialize a BSpline set
  Subroutine BSplineInit(s,NumberOfNodes,order,grid,IOSTAT)
    !
    !> Bspline basis set to be initialized
    Class(ClassBSpline), intent(inout)           :: s
    !> Number of Nodes
    integer            , intent(in)              :: NumberOfNodes
    !> Spline order
    integer            , intent(in)              :: order
    !> Grid of nodes (without repetition of degenerate nodes)
    DoublePrecision    , intent(in)              :: grid(1:NumberOfNodes)
    !> IOSTAT = 0 on successful exit, IOSTAT/=0 otherwise
    !! If IOSTAT is not specified, the programs terminates with
    !! an Assertion on error.
    integer            , intent(out),   optional :: IOSTAT

    !Local variables
    integer         :: i, Bs, Status
    procedure(D2DFUN), pointer :: FunPtr

    if(present(IOSTAT)) IOSTAT=0
    if(.not.INITIALIZED)call Init_Spline_Module()


    !.. Check input
    call CheckNodesAndOrder(NumberOfNodes,Order,Status)
    if(Status/=0)then
       if(present(IOSTAT))then
          IOSTAT=Status
          return
       endif
       call Assert("Invalid Number of Nodes / Order")
    endif


!!$    s%CARFUNINIT= .FALSE.
    s%NNodes    = NumberOfNodes
    s%Order     = order
    s%NBSplines = NumberOfNodes + order - 2


    !.. Check Grid format
    if(  LBOUND(Grid,1) > 1             .or. &
         UBOUND(Grid,1) < NumberOfNodes ) then
       if(present(IOSTAT))then
          IOSTAT=-4
          return
       endif
       call Assert("Invalid node list bounds")
    endif
    do i=2,NumberOfNodes
       if( Grid(i) < Grid(i-1) )then
          if(present(IOSTAT))then
             IOSTAT=1
             return
          else
             call Assert("Decreasing list of nodes!")
          endif
       endif
    enddo


    !.. Maps the Grid to [0:n-1]; hence, the upper bound
    !   of the i-th BSpline support is g(i)
    if( allocated( s%Grid ) ) deallocate( s%Grid )
    allocate( s%Grid( -(order-1) : NumberOfNodes-1 + order-1 ) )
    s%Grid( 1-order : -1 ) = Grid(1)
    s%Grid( 0:NumberOfNodes-1 ) = Grid(1:NumberOfNodes)
    s%Grid( NumberOfNodes:NumberOfNodes+order-2 ) = Grid(NumberOfNodes)


    !.. Initializes the BSpline coefficients normalized
    !   so that the B-spline set is a partition of unity
    !..
    if( allocated( s%c ) ) deallocate( s%c )
    allocate( s%c( 0:order-1, 0:order-1, 1:s%NBSplines ) )
    s%c = 0.d0
    do Bs=1,s%NBSplines
       call ComputeCoefficientsSingleBSpline(s%Order,s%Grid(Bs-s%Order),s%c(0,0,Bs))
    enddo
    s%INITIALIZED=.TRUE.


    !.. Initialize the normalization factors
    !   ||Bs_i(x)*s%f(i)||_L2=1.d0
    !..
    if( allocated( s%f ) ) deallocate( s%f )
    allocate( s%f( 1 : s%NBSplines ) )
    s%f = 0.d0
    funPtr=>Unity
    do i = 1, s%NBSplines
       s%f(i)=1.d0/sqrt(s%Integral(FunPtr,i,i))
    enddo


  end Subroutine BSplineInit


  Pure DoublePrecision function Unity(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y=1.d0
  end function Unity


  !> Given the order and a set of order+1 consecutive nodes,
  !> builds the corresponding (uniquely defined) B-spline.
  subroutine ComputeCoefficientsSingleBSpline(order,vec,Coefficients)
    !
    !> Order of the B-spline = number of (possibly degenerate) 
    !> intervals that form the B-spline support = Polynomials degree + 1
    !> (e.g., order=4 means piecewise cubic polynomials, which require
    !> four parameters: \f$ a_0 + a_1 x + a_2 * x^2 + a_3 * x^3 \f$.
    integer        , intent(in)  :: order
    !
    !> vector of order+1 nodes. For Order = 5, for example, vec contains 
    !> the six bounds of the five consecutive intervals that form the 
    !> support of the B-spline, indexed from ONE to SIX.
    !
    !  vec(1)  vec(2)    vec(3)  vec(4) vec(5)    vec(6)
    !    |-------|---------|-------|------|---------|
    DoublePrecision, intent(in)  :: vec(1:order+1)
    !
    !> Matrix of polynomial coefficients in each interval.
    DoublePrecision, intent(out) :: Coefficients(0:order-1,0:order-1)
    !
    !Local variables
    integer         :: k, inter, tempBs
    DoublePrecision :: SpanBSpline
    DoublePrecision :: SpanCurrentInterval
    DoublePrecision :: SpanPriorToInterval
    DoublePrecision :: SpanToTheEnd
    DoublePrecision :: FactorConstant
    DoublePrecision :: FactorMonomial
    DoublePrecision, allocatable :: coe(:,:,:,:)
    !
    !.. Computes the B-spline in each available
    !   sequence of order consecutive intervals
    !
    !
    allocate( coe(0:order-1,0:order-1,order,order) )
    coe=0.d0
    coe(0,0,1,:)=1.d0
    !
    !.. Cycle over the spline order
    !
    do k = 2, order
       !
       !.. Cycle over the (temporary) Bsplines available 
       !   for any given order, which must all be computed
       !
       do tempBs = 1, order - k + 1
          !
          !.. Each B-spline of order k is obtained as the sum of 
          !   the contribution from the two B-splines of order k-1
          !   whose support is comprised in that of the final 
          !   B-spline.
          !
          !.. Contribution of the first B-spline of order k-1
          !   counted only if the B-spline is not degenerate
          !..
          SpanBSpline = vec(tempBs+k-1) - vec(tempBs)
          if( SpanBSpline > NOD_THRESHOLD )then
             do inter=0,k-2
                !
                !.. Constant Component
                SpanPriorToInterval = vec(tempBs+inter) - vec(tempBs)
                FactorConstant = SpanPriorToInterval / SpanBSpline 
                coe(0:k-2,inter,k,tempBs) = coe(0:k-2,inter,k,tempBs) + &
                     coe(0:k-2,inter,k-1,tempBs) * FactorConstant
                !
                !.. Monomial Component
                SpanCurrentInterval = vec(tempBs+inter+1) - vec(tempBs+inter)
                FactorMonomial = SpanCurrentInterval / SpanBSpline
                coe(1:k-1,inter,k,tempBs) = coe(1:k-1,inter,k,tempBs) + &
                     coe(0:k-2,inter,k-1,tempBs) * FactorMonomial
                !
             enddo
          endif
          !
          !.. Contribution of the second B-spline of order k-1
          !   counted only if the B-spline is not degenerate
          !..
          SpanBSpline = vec(tempBs+k) - vec(tempBs+1)
          if( SpanBSpline > NOD_THRESHOLD )then
             do inter=1,k-1
                !
                !.. Constant Component
                SpanToTheEnd = vec(tempBs+k)-vec(tempBs+inter)
                FactorConstant = SpanToTheEnd / SpanBSpline
                coe(0:k-2,inter,k,tempBs) = coe(0:k-2,inter,k,tempBs) + &
                     coe(0:k-2,inter-1,k-1,tempBs+1) * FactorConstant
                !
                !.. Monomial Component
                SpanCurrentInterval = vec(tempBs+inter+1) - vec(tempBs+inter)
                FactorMonomial = SpanCurrentInterval / SpanBSpline
                coe(1:k-1,inter,k,tempBs) = coe(1:k-1,inter,k,tempBs) - &
                     coe(0:k-2,inter-1,k-1,tempBs+1) * FactorMonomial
                !
             enddo
          endif
          !
       enddo
    enddo
    !
    Coefficients=coe(:,:,order,1)
    deallocate(coe)
    !
  end subroutine ComputeCoefficientsSingleBSpline


  !> Deallocates a ClassBSpline variable
  !*** Should be changed to the FinalizeClassBSpline
  subroutine BSplineFree(s)
    Class(ClassBSpline), intent(inout) :: s
    if(allocated(s%Grid))deallocate(s%Grid)
    if(allocated(s%f))deallocate(s%f)
    if(allocated(s%c))deallocate(s%c)
!!$    if(allocated(s%carfun))deallocate(s%carfun)
!!$    if(allocated(s%carfunt))deallocate(s%carfunt)
!!$    if(allocated(s%vec))deallocate(s%vec)
    s%NNodes=0
    s%Order=0
    s%NBSplines=0
    s%INITIALIZED=.FALSE.
!!$    s%CARFUNINIT=.FALSE.
    return
  end subroutine BSplineFree


  !> Returns n if x is in \f$(g_n,g_{n+1}]\f$.
  integer function which_interval(x,s)
    !
    DoublePrecision   , intent(in) :: x
    type(ClassBSpline), intent(in) :: s
    !
    integer, save :: i1 = 0
    integer       :: i2
    !
    if(x>s%Grid(i1).and.x<=s%Grid(i1+1))then
       which_interval=i1
       return
    endif
    if(x>s%Grid(i1+1).and.x<=s%Grid(min(i1+2,s%NNodes+s%Order-1)))then
       which_interval=i1+1
       return
    endif
    i1=0
    i2=0
    !
    if(x<s%Grid(0))then
       which_interval=-1
       return
    elseif(x>s%Grid(s%NNodes-1))then
       which_interval=s%NNodes!-1
       return
    endif
    if(x<s%Grid(i1))i1=0
    if(x>s%Grid(i2))i2=s%NNodes-1
    do while(i2-i1>1)
       if(x>s%Grid((i2+i1)/2))then
          i1=(i1+i2)/2
       else
          i2=(i1+i2)/2
       endif
    enddo
    which_interval=i1
    return
  end function which_interval


  !> Computes the n-th derivative of B-spline Bs, \f$ n=0,1,2,\ldots\f$
  !!
  ! {{{ Detailed Description:

  !! ---------------------
  !! B-splines are positive: expansion coefficients of represented functions are 
  !! akin to the functions themselves. Without rapid oscillations of coefficients, 
  !! the cancellation errors are kept at minimum. 
  !! In particular \cite{deBoor},
  !! if \f$f=\sum_i B_ic_i\f$
  !! \f{equation}
  !! \left|c_i-\frac{M+m}{2}\right|\leq D_{k}\frac{M-m}{2},\quad m
  !!  =\min_{x\in[a,b]} f(x),\,\,M=\max_{x\in[a,b]}f(x)
  !! \f}
  !! where \f$D_k\f$ is a constant that depends only on \f$k\f$ and not on the
  !! particular partition of the \f$[a,b]\f$ interval.
  !! At any point, the evaluation of a linear combination of B-splines 
  !! requires the evaluation of \f$k\f$ basis functions only. 

  ! }}}
  DoublePrecision function BSplineEval(s,x,Bs,n_)
    Class(ClassBSpline), intent(in) :: s
    DoublePrecision    , intent(in) :: x
    integer            , intent(in) :: Bs
    integer, optional  , intent(in) :: n_
    !
    integer         :: i,j,k,n
    DoublePrecision :: r,a,w,OneOverInterval,OneOverIntervalToTheN
    !
    call CheckInitialization(s)
    !
    BSplineEval=0.d0
    if(Bs<1.or.Bs>s%NBSplines)return
    !
    n=0; if(present(n_)) n=n_
    if(n<0) call Assert("Invalid derivative order in BSpline Evaluation")
    if(n>=s%Order)return
    !
    i=which_interval(x,s)
    if(  i<max(0,Bs-s%Order).or.&
         i>min(s%NNodes-1,Bs-1))return
    !
    j=i-Bs+s%Order
    OneOverInterval=1.d0/(s%Grid(i+1)-s%Grid(i))
    r=(x-s%Grid(i))*OneOverInterval
    !
    a=1.d0
    w=0.d0
    do k=n,s%Order-1
       w=w+a*parfactmat(k,n)*s%c(k,j,Bs)
       a=a*r
    enddo
    !
    OneOverIntervalToTheN=1.d0
    do k=1,n
       OneOverIntervalToTheN=OneOverIntervalToTheN*OneOverInterval
    enddo
    BSplineEval=w*OneOverIntervalToTheN
    !
    return
  end function BSplineEval


  !> Computes the k-th derivative of a function expressed in terms of B-splines.
  DoublePrecision function BSplineFunctionEval(s,x,fc,k,SKIP_FIRST,Bsmin,Bsmax)
    class(ClassBSpline), intent(in) :: s
    DoublePrecision    , intent(in) :: x
    DoublePrecision    , intent(in) :: fc(:)
    integer, optional  , intent(in) :: k
    logical, optional  , intent(in) :: SKIP_FIRST
    integer, optional  , intent(in) :: Bsmin,Bsmax
    integer :: Bsmi, Bsma
    integer :: Bs,i,nskip
    call CheckInitialization(s)
    BSplineFunctionEval=0.d0
    i=which_interval(x,s)
    if(i<0.or.i>s%NNodes-2)return
    Bsmi=1
    nskip=0
    if(present(SKIP_FIRST))then
       if(SKIP_FIRST)then
          Bsmi=2
          nskip=1
       endif
    endif
    if(present(Bsmin))then
       Bsmi = max(Bsmin,Bsmi)
       nskip= max(Bsmi-1,nskip)
    endif
    Bsmi=max(Bsmi,i+1)
    Bsma = min(s%NBSplines,i+s%Order)
    if(present(Bsmax)) Bsma = min(Bsmax,Bsma)
    do Bs=Bsmi,Bsma
       BSplineFunctionEval=BSplineFunctionEval+s%Eval(x,Bs,k)*fc(Bs-nskip)
    enddo
    return
  end function BSplineFunctionEval
  

  !> Computes the  k-th derivatives of a choosen B-spline evaluated in a position array. 
  subroutine BSplineTabulate(s,ndata,xvec,yvec,Bs,n_) 
    Class(ClassBSpline), intent(in) :: s
    integer            , intent(in) :: ndata
    DoublePrecision    , intent(in) :: xvec(1:ndata)
    DoublePrecision    , intent(out):: yvec(1:ndata)
    integer            , intent(in) :: Bs
    integer, optional  , intent(in) :: n_
    !
    integer :: i,n
    n=0;if(present(n_))n=n_
    do i=1,ndata
       yvec(i)=s%Eval(xvec(i),Bs,n)
    enddo
  end subroutine BSplineTabulate


  !> Computes the  k-th derivatives of a function expressed in terms of B-splines and evaluated in an array of positions. 
  subroutine BSplineFunctionTabulate(s,ndata,xvec,yvec,FunVec,n_,SKIP_FIRST) 
    Class(ClassBSpline), intent(in) :: s
    integer            , intent(in) :: ndata
    DoublePrecision    , intent(in) :: xvec(1:ndata)
    DoublePrecision    , intent(out):: yvec(1:ndata)
    DoublePrecision    , intent(in) :: FunVec(:)
    integer, optional  , intent(in) :: n_
    logical, optional  , intent(in) :: SKIP_FIRST
    logical :: SKIP_FIRST_LOC
    !
    integer :: i,n
    SKIP_FIRST_LOC=.FALSE.
    if(present(SKIP_FIRST))then
       SKIP_FIRST_LOC=SKIP_FIRST
    endif
    n=0;if(present(n_))n=n_
    do i=1,ndata
       yvec(i)=s%Eval(xvec(i),FunVec,n,SKIP_FIRST=SKIP_FIRST_LOC)
    enddo
  end subroutine BSplineFunctionTabulate


  !> Computes the integral
  !! \f[
  !!    \int_{a}^{b} \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
  !!                 \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}} r^{2}dr
  !! \f]
  !! Where \f$f(r)\f$ is a local operator, and if a break point \f$BP\f$ 
  !! is introduced, then the integral is splitted in two parts
  !! \f[
  !!    \int_{a}^{BP} \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
  !!                  \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}} r^{2}dr + 
  !!    \int_{BP}^{b} \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
  !!                  \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}}r^{2}dr
  !! \f]
  DoublePrecision function BSplineIntegral( &
       s                 , &
       FunPtr            , & 
       Bs1               , &
       Bs2               , &
       BraDerivativeOrder, &
       KetDerivativeOrder, &
       LowerBound        , &
       UpperBound        , &
       BreakPoint        , &
       parvec            ) &
       result( Integral )
    !
    Class(ClassBSpline)      , intent(in) :: s
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: Bs1
    integer                  , intent(in) :: Bs2
    integer        , optional, intent(in) :: BraDerivativeOrder
    integer        , optional, intent(in) :: KetDerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    DoublePrecision, optional, intent(in) :: BreakPoint
    DoublePrecision, optional, intent(in) :: Parvec(*)
    !
    integer         :: n1, n2, iostat
    DoublePrecision :: a, b
    real(kind(1d0)) :: NewBreakPoint
    !
    Integral=0.d0
    call CheckInitialization(s)
    call CheckParameters( iostat )
    if( iostat/=0 )return
    if( present(BreakPoint) )then
       !***
       if ( BreakPoint < a ) then
          NewBreakPoint = a + epsilon(1d0)
       elseif (  BreakPoint > b ) then
          NewBreakPoint = b - epsilon(1d0)
       end if
       !***
       Integral = &
            BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2, a, BreakPoint,    Parvec ) + &
            BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2, NewBreakPoint, b, Parvec )
    else
       Integral = BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2,a,b, Parvec )
    endif
    !
    return
    !
  contains
    !
    subroutine CheckParameters( IOStat )
      integer, intent(out) :: IOStat
      call CheckBSplineIndexes( IOStat ); if( IOStat /= 0 ) return
      call CheckIntegralBounds( IOStat ); if( IOStat /= 0 ) return
      call CheckDerivativeOrder
    end subroutine CheckParameters
    !
    subroutine CheckBSplineIndexes( IOStat )
      integer, intent(out) :: IOStat
      IOStat=1
      if( min(Bs1,Bs2) <  1           ) return      
      if( max(Bs1,Bs2) >  s%NBSplines ) return
      if( abs(Bs1-Bs2) >= s%Order     ) return
      IOStat=0
    end subroutine CheckBSplineIndexes
    !
    subroutine CheckIntegralBounds( IOStat )
      integer, intent(out) :: IOStat
      DoublePrecision :: aPlus, bMinus
      IOStat=1
      a = s%Grid( max(Bs1,Bs2) - s%Order )
      b = s%Grid( min(Bs1,Bs2) )
      if(present(LowerBound))then
         if( LowerBound >= b )return
         a=max(a,LowerBound)
      endif
      if(present(UpperBound))then
         if( UpperBound <= a )return
         b=min(b,UpperBound)
      endif
      IOStat=0
    end subroutine CheckIntegralBounds
    !
    subroutine CheckDerivativeOrder
      n1=0; if(present(BraDerivativeOrder)) n1=BraDerivativeOrder
      n2=0; if(present(KetDerivativeOrder)) n2=KetDerivativeOrder
    end subroutine CheckDerivativeOrder
    !
  end function BSplineIntegral



  ! {{{ Detailed Description

  !> Compute the integral
  !! \f[
  !!    \int_a^b dr 
  !!    \frac{d^{n_1}B_1(r)}{dr^{n_1}} f(r)
  !!    \frac{d^{n_2}B_2(r)}{dr^{n_2}},
  !! \f]
  !! where \f$f(r)\f$ is a *smooth* function in the
  !! intersection of the supports of the two B-splines.
  !! If the first (last) boundary of the integration 
  !! interval is not specified, then the value a=-oo
  !! (b=+oo) is assumed.
  !!
  !!
  !! Properties of B-spline integrals:
  !! ---------------------------------
  !! Since the support \f$D_i\f$ of the \f$i\f$-th B-spline is formed by
  !! \f$k\f$ consecutive intervals, the integrals between two B-splines and a 
  !! local operator are zero unless their indices differ less than the order \f$k\f$:
  !! \f{equation}
  !! \langle B_i|O|B_j\rangle=\int_{D_i \cap D_j} B_i(x) o(x) B_j(x) dx.
  !! \f}
  !! As a consequence, matrices are sparse and often band-diagonal. 
  !! On uniform grids, matrices of translationally invariant operators, that 
  !! is with kernel \f$o(x,y)=o(x-y)\f$, are Toeplitz:
  !! \f{eqnarray}
  !! \langle B_i|O|B_j\rangle=\langle B_{i+n}|O|B_{j+n}\rangle.
  !! \f}
  !! If an hermitian operator \f$O\f$ is both local and translationally invariant, 
  !! like the identity and the kinetic energy, its matrix on a uniform grid
  !! is Toeplitz and banded, in other terms it is defined by just \f$k\f$ numbers. 
  !! 
  !! With Gauss-type integration technique, the matrix elements of operators 
  !! with polynomial kernels are exact. The error \f$\epsilon\f$ in the approximation
  !! of a \f$\mathcal{C}^k\f$ function \f$f\f$ with B-splines, is bounded by
  !! \f[
  !!     \epsilon\leq \mathrm{const}_k|\mathbf{t}|^k\|D^kf\|,\quad |\mathbf{t}|=
  !!         \max_i(t_{i+1}-t_i),\quad \|f\|=\max_{x\in[a,b]}|f(x)|
  !! \f]
  !!

  ! }}}
  DoublePrecision function BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2,a,b,parvec) &
       result( Integral )
    !
    !.. Assumes that the input data have been sanitized
    !
    Class(ClassBSpline), intent(in) :: s
    procedure(D2DFun)  , pointer    :: FunPtr
    integer            , intent(in) :: Bs1, Bs2
    integer            , intent(in) :: n1, n2
    DoublePrecision    , intent(in) :: a, b
    DoublePrecision, optional, intent(in) :: parvec(*)
    !
    integer         :: Interval, IntervalMin, IntervalMax
    DoublePrecision :: LowerBound, UpperBound
    !
    DoublePrecision :: PartialIntegral
    DoublePrecision :: IntervalWidth
    DoublePrecision :: Radius
    integer         :: iGauss
    !
    DoublePrecision :: aPlus, bMinus
    !
    Integral=0.d0
    !
    aPlus  = UpperLimitTo( a )
    bMinus = LowerLimitTo( b )
    IntervalMin = which_interval( aPlus,  s )
    IntervalMax = which_interval( bMinus, s )
    !
    do Interval = IntervalMin, IntervalMax
       LowerBound = max( a, s%Grid( Interval   ) )
       UpperBound = min( b, s%Grid( Interval+1 ) )
       !
       IntervalWidth = UpperBound - LowerBound
       if( IntervalWidth < NOD_THRESHOLD ) cycle
       !
       PartialIntegral = 0.d0
       do iGauss=1,NGAUSS
          !
          Radius = LowerBound + IntervalWidth * Gauss_Points( iGauss )
          PartialIntegral = PartialIntegral + &
               FunPtr(Radius, Parvec) * &
               s%Eval(Radius,Bs1,n1)  * &
               s%Eval(Radius,Bs2,n2)  * &
               Gauss_weight( iGauss )
          !
       enddo
       PartialIntegral = PartialIntegral * IntervalWidth
       !
       Integral = Integral + PartialIntegral
       !
    enddo
    !
  end function BSplineDriverIntegral



  !> Computes the integral
  !! \f[
  !!    \int_a^b dr 
  !!    \frac{d^{n_1}B_1(r)}{dr^{n_1}} f(r)
  !! \f]
  !! where \f$f(r)\f$ is a *smooth* function in the B-splines support.
  !! If the first (last) boundary of the integration interval is not specified,
  !! then the value \f$a=-\infty\f$ and \f$b=\infty\f$ is assumed.
  DoublePrecision function BSplineMonoIntegral( &
       s              , &
       FunPtr         , &
       Bs1            , &
       DerivativeOrder, &
       LowerBound     , &
       UpperBound     , &
       BreakPoint     , &
       Parvec         ) &
       result( Integral )
    !
    Class(ClassBSpline)      , intent(in) :: s
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: Bs1
    integer        , optional, intent(in) :: DerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    DoublePrecision, optional, intent(in) :: BreakPoint
    DoublePrecision, optional, intent(in) :: Parvec(*)
    !
    integer         :: n1
    DoublePrecision :: a, b
    !
    integer         :: IOStat
    real(kind(1d0)) :: NewBreakPoint
    !
    Integral=0.d0
    call CheckInitialization(s)
    call CheckParameters( IOStat )
    if( IOStat/=0 )return
    if( present(BreakPoint) )then
    !***
    if ( BreakPoint < a ) then
       NewBreakPoint = a + epsilon(1d0)
    elseif (  BreakPoint > b ) then
       NewBreakPoint = b - epsilon(1d0)
    end if
    !***
       Integral = &
            BSplineMonoIntegralDriver(s,FunPtr,Bs1,n1,a,NewBreakPoint,Parvec) + &
            BSplineMonoIntegralDriver(s,FunPtr,Bs1,n1,NewBreakPoint,b,Parvec)
    else
       Integral = &
            BSplineMonoIntegralDriver(s,FunPtr,Bs1,n1,a,b,Parvec)
    endif
    !
  contains
    !
    subroutine CheckParameters( IOStat )
      integer, intent(out) :: IOStat
      call CheckBSplineIndexes( IOStat ); if( IOStat /= 0 ) return
      call CheckIntegralBounds( IOStat ); if( IOStat /= 0 ) return
      call CheckDerivativeOrder
    end subroutine CheckParameters
    !
    subroutine CheckBSplineIndexes( IOStat )
      integer, intent(out) :: IOStat
      IOStat=1
      if( Bs1 <  1           ) return      
      if( Bs1 >  s%NBSplines ) return
      IOStat=0
    end subroutine CheckBSplineIndexes
    !
    subroutine CheckIntegralBounds( IOStat )
      integer, intent(out) :: IOStat
      DoublePrecision :: aPlus, bMinus
      IOStat=1
      a = s%Grid( Bs1 - s%Order )
      b = s%Grid( Bs1 )
      if(present(LowerBound))then
         if( LowerBound >= b )return
         a=max(a,LowerBound)
      endif
      if(present(UpperBound))then
         if( UpperBound <= a )return
         b=min(b,UpperBound)
      endif
      IOStat=0
    end subroutine CheckIntegralBounds
    !
    subroutine CheckDerivativeOrder
      n1=0; if(present(DerivativeOrder)) n1=DerivativeOrder
    end subroutine CheckDerivativeOrder
    !
  end function BSplineMonoIntegral




  !> Computes the integral
  !! \f[
  !!    \int_a^b dr 
  !!    \frac{d^{n_1}B_1(r)}{dr^{n_1}} f(r)
  !! \f]
  !! where \f$f(r)\f$ is a *smooth* function in the B-splines support.
  !! If the first (last) boundary of the integration interval is not 
  !! specified, then the value \f$a=-\infty\f$ and \f$b=\infty\f$ is assumed.
  DoublePrecision function BSplineMonoIntegralDriver(s,FunPtr,Bs1,n1,a,b,Parvec) &
       result( Integral )
    !
    Class(ClassBSpline)      , intent(in) :: s
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: Bs1
    integer        , optional, intent(in) :: n1
    DoublePrecision, optional, intent(in) :: a,b
    DoublePrecision, optional, intent(in) :: parvec(*)
    !
    integer         :: Interval, IntervalMin, IntervalMax
    DoublePrecision :: LocalLowerBound, LocalUpperBound
    !
    DoublePrecision :: PartialIntegral
    DoublePrecision :: IntervalWidth
    DoublePrecision :: Radius
    integer         :: iGauss
    !
    DoublePrecision :: aPlus, bMinus
    !
    Integral=0.d0
    !
    aPlus  = UpperLimitTo( a )
    bMinus = LowerLimitTo( b )
    IntervalMin = which_interval( aPlus,  s )
    IntervalMax = which_interval( bMinus, s )
    !
    do Interval = IntervalMin, IntervalMax
       !
       LocalLowerBound = max( a, s%Grid( Interval   ) )
       LocalUpperBound = min( b, s%Grid( Interval+1 ) )
       !
       IntervalWidth = LocalUpperBound - LocalLowerBound
       if( IntervalWidth < NOD_THRESHOLD ) cycle
       !
       PartialIntegral = 0.d0
       do iGauss=1,NGAUSS
          !
          Radius = LocalLowerBound + IntervalWidth * Gauss_Points( iGauss )
          PartialIntegral = PartialIntegral + &
               FunPtr(Radius, parvec) * &
               s%Eval(Radius,Bs1,n1)  * &
               Gauss_weight( iGauss )
          !
       enddo
       PartialIntegral = PartialIntegral * IntervalWidth
       !
       Integral = Integral + PartialIntegral
       !
    enddo
    !
  end function BSplineMonoIntegralDriver


!!$  !> Computes the principal-part B-spline integral
!!$  !! \f[
!!$  !!  \int dr Bs_1(r)\,\frac{\mathcal{P}}{r-r_0}\,Bs_2(r)
!!$  !! \f]
!!$  DoublePrecision function BSplineCauchyIntegral(s,Bs1,Bs2,SingularPoint) &
!!$       result( Integral )
!!$    Class(ClassBSpline), intent(in) :: s
!!$    integer            , intent(in) :: Bs1, Bs2
!!$    DoublePrecision    , intent(in) :: SingularPoint
!!$    !
!!$    integer :: is,id,Bsma,indis
!!$    DoublePrecision :: ErrorEstimate
!!$    DoublePrecision :: LowerBound
!!$    DoublePrecision :: UpperBound
!!$    DoublePrecision :: SingularPoint_
!!$    !
!!$    Integral=0.d0
!!$    call CheckInitialization(s)
!!$    call CheckParameters()
!!$    !
!!$    LowerBound=s%Grid(Bsma-s%Order)
!!$    UpperBound=s%Grid(Bsma-s%Order+1+min(s%Order-1-indis,s%NBSplines-Bsma))
!!$    !
!!$    SingularPoint_=SingularPoint
!!$    if( abs( SingularPoint - LowerBound ) < UpperLimitTo( LowerBound ) - LowerBound )&
!!$         SingularPoint_ = LowerLimitTo( LowerBound )
!!$    if( abs( SingularPoint - UpperBound ) < UpperLimitTo( UpperBound ) - UpperBound )&
!!$         SingularPoint_ = UpperLimitTo( UpperBound )
!!$    !
!!$    call LocalIntegrationRoutine(   &
!!$         LowerBound    , &
!!$         UpperBound    , &
!!$         SingularPoint_, &
!!$         Integral      , &
!!$         ErrorEstimate )
!!$    !
!!$  contains
!!$    !
!!$    subroutine CheckParameters()
!!$      if(min(Bs1,Bs2)<1)return
!!$      Bsma=max(Bs1,Bs2)
!!$      if(Bsma>s%NBSplines)return
!!$      if(abs(Bs1-Bs2)>=s%Order)return
!!$      is=max(0,Bs2-Bs1)
!!$      id=max(0,Bs1-Bs2)
!!$      indis=abs(Bs1-Bs2)
!!$    end subroutine CheckParameters
!!$    !
!!$    Pure DoublePrecision function LocalFunction(x) result(y)
!!$      DoublePrecision, intent(in) :: x
!!$      y = s%Eval(x,Bs1) * s%Eval(x,Bs2)
!!$    end function LocalFunction
!!$    !
!!$    !
!!$    subroutine LocalIntegrationRoutine(&
!!$         LowerBound    , &
!!$         UpperBound    , &
!!$         SingularPoint , &
!!$         Integral      , &
!!$         ErrorEstimate )
!!$      !
!!$      implicit none
!!$      !
!!$      DoublePrecision           , intent(in)  :: LowerBound
!!$      DoublePrecision           , intent(in)  :: UpperBound
!!$      DoublePrecision           , intent(in)  :: SingularPoint
!!$      DoublePrecision           , intent(out) :: Integral
!!$      DoublePrecision           , intent(out) :: ErrorEstimate
!!$      !
!!$      !.. Error Codes
!!$      integer, parameter :: NORMAL_AND_RELIABLE_TERMINATION_ACHIEVED = 0
!!$      integer, parameter :: MAXIMUM_NUMBER_OF_SUBDIVISIONS_ACHIEVED  = 1
!!$      integer, parameter :: OCCURRENCE_OF_ROUNDOFF_ERROR_DETECTED    = 2
!!$      integer, parameter :: EXTREMELY_BAD_INTEGRAND_BEHAVIOR_OCCURS  = 3
!!$      integer, parameter :: THE_INPUT_IS_INVALID                     = 6
!!$      !
!!$      !.. Accuracy Parameters
!!$      DoublePrecision, parameter :: ABSOLUTE_ACCURACY_REQUESTED     = 1.d-8 
!!$      DoublePrecision, parameter :: RELATIVE_ACCURACY_REQUESTED     = 1.d-8 
!!$      integer        , parameter :: INITIAL_NUMBER_OF_SUBINTERVALS  = 7
!!$      integer        , parameter :: INCREASE_NUMBER_OF_SUBINTERVALS = 4
!!$      integer        , parameter :: MAXIMUM_NUMBER_OF_SUBINTERVALS  = 31
!!$      integer        , parameter :: LENGTH_WORK = 4 * MAXIMUM_NUMBER_OF_SUBINTERVALS
!!$      !
!!$      integer :: ierror
!!$      integer :: ActualNumberOfSubintervals
!!$      integer :: NumberOfIntegralEvaluations
!!$      integer :: CurrentMaxNumSubintervals
!!$      !
!!$      integer        , save :: iWork( MAXIMUM_NUMBER_OF_SUBINTERVALS ) = 0.d0
!!$      DoublePrecision, save :: Work( LENGTH_WORK ) = 0
!!$      !
!!$      interface
!!$         subroutine dqawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,ier,&
!!$              limit,lenw,last,iwork,work)
!!$           double precision a,abserr,b,c,epsabs,epsrel,f,result,work
!!$           integer ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
!!$           dimension iwork(limit),work(lenw)
!!$           external f
!!$         end subroutine dqawc
!!$      end interface
!!$
!!$      CurrentMaxNumSubintervals = INITIAL_NUMBER_OF_SUBINTERVALS
!!$      do
!!$         ierror=0
!!$         call dqawc(                       & 
!!$              LocalFunction              , &
!!$              LowerBound                 , &
!!$              UpperBound                 , &
!!$              SingularPoint              , &
!!$              ABSOLUTE_ACCURACY_REQUESTED, &
!!$              RELATIVE_ACCURACY_REQUESTED, &
!!$              Integral                   , &
!!$              ErrorEstimate              , &
!!$              NumberOfIntegralEvaluations, &
!!$              ierror                     , &
!!$              CurrentMaxNumSubintervals  , &
!!$              LENGTH_WORK                , &
!!$              ActualNumberOfSubintervals , &
!!$              iWork                      , &
!!$              Work                       )
!!$         select case( ierror )
!!$         case( NORMAL_AND_RELIABLE_TERMINATION_ACHIEVED )
!!$            return
!!$         case( MAXIMUM_NUMBER_OF_SUBDIVISIONS_ACHIEVED )
!!$            CurrentMaxNumSubintervals = &
!!$                 CurrentMaxNumSubintervals + &
!!$                 INCREASE_NUMBER_OF_SUBINTERVALS
!!$         case( OCCURRENCE_OF_ROUNDOFF_ERROR_DETECTED )
!!$            CurrentMaxNumSubintervals = &
!!$                 CurrentMaxNumSubintervals + &
!!$                 INCREASE_NUMBER_OF_SUBINTERVALS
!!$         case( EXTREMELY_BAD_INTEGRAND_BEHAVIOR_OCCURS )
!!$            exit
!!$         case( THE_INPUT_IS_INVALID )
!!$            exit
!!$         end select
!!$         if( CurrentMaxNumSubintervals > MAXIMUM_NUMBER_OF_SUBINTERVALS ) exit
!!$      enddo
!!$      return
!!$    end subroutine LocalIntegrationRoutine
!!$    !
!!$  end function BSplineCauchyIntegral



  DoublePrecision function MixedBSplineIntegral( &
       s1, s2, FunPtr, Bs1, Bs2, BraDerivativeOrder, KetDerivativeOrder, LowerBound, UpperBound, Parvec ) &
       result( Integral )
    !
    Class(ClassBSpline)      , intent(in) :: s1, s2
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: Bs1, Bs2
    integer        , optional, intent(in) :: BraDerivativeOrder, KetDerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound, UpperBound
    DoublePrecision, optional, intent(in) :: Parvec(*)
    !
    integer         :: n1, n2
    DoublePrecision :: a, b
    !
    integer         :: Interval
    integer         :: IntervalMin1, IntervalMax1
    integer         :: IntervalMin2, IntervalMax2
    DoublePrecision :: LocalLowerBound, LocalUpperBound
    !
    DoublePrecision :: PartialIntegral
    DoublePrecision :: IntervalWidth
    DoublePrecision :: Radius
    integer         :: iGauss
    integer         :: IOStat
    integer                      :: TotalIntervals
    DoublePrecision, allocatable :: MergedGrid(:)
    !
    Integral=0.d0
    call CheckInitialization(s1)
    call CheckInitialization(s2)
    call CheckParameters( IOStat )
    call MergeGrid
    if( IOStat/=0 )return
    !
    do Interval = 1, TotalIntervals
       LocalLowerBound = max( a, MergedGrid( Interval   ) )
       LocalUpperBound = min( b, MergedGrid( Interval+1 ) )
       !
       IntervalWidth = LocalUpperBound - LocalLowerBound
       if( IntervalWidth < NOD_THRESHOLD ) cycle
       !
       PartialIntegral = 0.d0
       do iGauss=1,NGAUSS
          !
          Radius = LocalLowerBound + IntervalWidth * Gauss_Points( iGauss )
          PartialIntegral = PartialIntegral + &
               FunPtr(Radius, Parvec) * &
               s1%Eval(Radius,Bs1,n1) * &
               s2%Eval(Radius,Bs2,n2) * &
               Gauss_weight( iGauss )
          !
       enddo
       PartialIntegral = PartialIntegral * IntervalWidth
       !
       Integral = Integral + PartialIntegral
       !
    enddo
    !
    return
    !
  contains
    !
    subroutine CheckParameters( IOStat )
      integer, intent(out) :: IOStat
      call CheckBSplineIndexes( IOStat ); if( IOStat /= 0 ) return
      call CheckIntegralBounds( IOStat ); if( IOStat /= 0 ) return
      call CheckDerivativeOrder
    end subroutine CheckParameters
    !
    subroutine CheckBSplineIndexes( IOStat )
      !
      integer, intent(out) :: IOStat
      !
      IOStat=1
      !
      if( Bs1 <  1            ) return      
      if( Bs1 >  s1%NBSplines ) return
      !
      if( Bs2 <  1            ) return      
      if( Bs2 >  s2%NBSplines ) return
      !
      IOStat=0
      !
    end subroutine CheckBSplineIndexes
    !
    subroutine CheckIntegralBounds( IOStat )
      !
      integer, intent(out) :: IOStat
      !
      DoublePrecision :: aPlus, bMinus
      DoublePrecision :: a1,b1,a2,b2
      !
      IOStat=1
      !
      a1 = s1%Grid( Bs1 - s1%Order )
      b1 = s1%Grid( Bs1 )
      !
      a2 = s2%Grid( Bs2 - s2%Order )
      b2 = s2%Grid( Bs2 )
      !
      a=max(a1,a2)
      b=min(b1,b2)
      !
      if(present(LowerBound))then
         if( LowerBound >= b )return
         a=max(a,LowerBound)
      endif
      if(present(UpperBound))then
         if( UpperBound <= a )return
         b=min(b,UpperBound)
      endif
      IOStat=0
      aPlus =UpperLimitTo(a)
      bMinus=LowerLimitTo(b)
      !
      IntervalMin1=which_interval(aPlus ,s1)
      IntervalMax1=which_interval(bMinus,s1)
      !
      IntervalMin2=which_interval(aPlus ,s2)
      IntervalMax2=which_interval(bMinus,s2)
      !
    end subroutine CheckIntegralBounds
    !
    subroutine CheckDerivativeOrder
      n1=0; if(present(BraDerivativeOrder)) n1=BraDerivativeOrder
      n2=0; if(present(KetDerivativeOrder)) n2=KetDerivativeOrder
    end subroutine CheckDerivativeOrder
    !
    subroutine MergeGrid
      !
      integer :: i1,i2,i
      !
      TotalIntervals=&
           IntervalMax1-IntervalMin1+1+&
           IntervalMax2-IntervalMin2+1
      allocate(MergedGrid(TotalIntervals+1))
      i1=IntervalMin1
      i2=IntervalMin2
      do i=1,TotalIntervals
         if(s1%Grid(i1)<=s2%Grid(i2).and.i1<=IntervalMax1)then
            MergedGrid(i)=s1%Grid(i1)
            i1=i1+1
         else
            MergedGrid(i)=s2%Grid(i2)
            i2=i2+1
         endif
      enddo
      MergedGrid(TotalIntervals+1)=min(s1%Grid(IntervalMax1+1),s2%Grid(IntervalMax2+1))
      !
    end subroutine MergeGrid
    !
  end function MixedBSplineIntegral


  Pure DoublePrecision function UpperLimitTo(x) result(xPlus)
    DoublePrecision, intent(in) :: x
    DoublePrecision :: tin,eps
    eps=epsilon(1.d0)
    tin=tiny(x)*(1.d0+eps)
    xPlus=(x+tin)*(1.d0+eps)
  end function UpperLimitTo


  Pure DoublePrecision function LowerLimitTo(x) result(xMinus)
    DoublePrecision, intent(in) :: x
    DoublePrecision :: tin,eps
    eps=epsilon(1.d0)
    tin=tiny(x)*(1.d0+eps)
    xMinus=(x-tin)*(1.d0-eps)
  end function LowerLimitTo



  !> Computes the "Spherical Bessel Transform" of a B-spline, defined as:
  !! \f[
  !!   F_{l}(p)=\int\frac{Bs_{n}(r)}{r}j_{l}(p\cdot r)r^{2}dr,
  !! \f]
  !! for all the angular momenta \f$\ell\f$ up to a specified lmax as \f$ res(\ell+1)=F_\ell(p)\f$.
  subroutine SphericalBesselTransform(s,p,Bs,lmax,res)
    Class(ClassBSpline), intent(in) :: s
    DoublePrecision    , intent(in) :: p
    integer            , intent(in) :: Bs
    integer            , intent(in) :: lmax
    DoublePrecision    , intent(out):: res(*)
    !Local variables
    DoublePrecision, parameter :: MIN_P=1.d-10
    integer         :: interval, pow, m
    DoublePrecision :: vjl(lmax+1),del,spline,a,x,pr,p_

    call CheckInitialization(s)

    if(lmax<0)return
    p_=max(p,MIN_P)
    res(1:lmax+1)=0.d0
    if(Bs>s%NBSplines.or.Bs<1)return
    !Ciclo sugli intervalli
    do interval=0,min(s%Order-1,s%NBSplines-Bs)
       a=s%Grid(Bs-s%Order+interval)
       del=s%Grid(Bs-s%Order+interval+1)-a
       if(del<NOD_THRESHOLD)cycle
       !Ciclo sui punti di Gauss
       do m=1,NGAUSS
          x=Gauss_Points(m)
          pr=p_*(a+del*x)
          call Spherical_Besselxj(pr,lmax,vjl)
          spline=0.d0
          do pow=0,s%Order-1
             spline=spline+s%c(pow,interval,Bs)*(x**(pow))
          enddo
          res(1:lmax+1)=res(1:lmax+1)+&
               Gauss_weight(m)*del*spline*vjl(1:lmax+1)
       enddo
    enddo
    res(1:lmax+1)=res(1:lmax+1)/p_

  contains

    subroutine Spherical_Besselxj(x,lmax,res)
      !Returns the x*j_l(x) for all l=0,1,...,lmax
      !as res(l+1)=x*j_l(x)
      DoublePrecision, intent(in) :: x
      integer        , intent(in) :: lmax
      DoublePrecision, intent(out):: res(lmax+1)
      DoublePrecision :: y
      integer         :: l
      if(lmax<0)return
      res(1:lmax+1)=0.d0
      if(x==0.d0)return
      y=1.d0/x
      res(1)=sin(x)
      if(lmax>=1)res(2)=y*res(1)-cos(x)
      do l=2,lmax
         res(l+1)=dble(2*l-1)*y*res(l)-res(l-1)
      enddo
      return
    end subroutine Spherical_Besselxj

  end subroutine SphericalBesselTransform


  ! {{{ Multipoles routines

  !> Computes the monoelectronic multipole
  !! \f[
  !!    \int_0^\infty dr B_1(R)\frac{r_<^\ell}{r_>^{\ell+1}}B_2(r)
  !! \f]
  !! where \f$r_<=\min(r,R)\f$ and \f$r_>=\max(r,R)\f$. 
  !!
  DoublePrecision function MonoElectronicMultipole(s,Bs1,Bs2,R,l) result(Multipole)
    Class(ClassBSpline), intent(in) :: s
    integer            , intent(in) :: Bs1, Bs2
    DoublePrecision    , intent(in) :: R
    integer            , intent(in) :: l
    !
    procedure(D2DFun), pointer :: FunPtr
    FunPtr => MultipolarPotential
    Multipole = s%Integral( FunPtr, Bs1, Bs2, BreakPoint = R )
  contains
    Pure DoublePrecision function MultipolarPotential(x,parvec) result(y)
      DoublePrecision, intent(in) :: x
      DoublePrecision, optional, intent(in) :: parvec(*)
      y = Pow(min(x,R),l)*Pow(max(x,R),-l-1)
    end function MultipolarPotential
  end function MonoElectronicMultipole


  Pure DoublePrecision function Pow(x,n)
    DoublePrecision, intent(in) :: x
    integer        , intent(in) :: n
    integer :: i
    Pow=1.d0
    if(n==0)return
    Pow=0.d0
    if(abs(x)<tiny(x))return
    Pow=1.d0
    do i=1,n
       Pow=Pow*x
    enddo
    do i=n,-1
       Pow=Pow/x
    enddo
  end function Pow


  !> Computes the bielectronic multipole
  !! \f[
  !!     \int_0^\infty dr\int_0^{\infty} dr' 
  !!      B_{n_1}(r) B_{n_2}(r') 
  !!      \frac{r_<^\ell}{r_>^{\ell+1}}
  !!      B_{n_3}(r) B_{n_4}(r') 
  !! \f]
  !! between non-normalized B-splines
  DoublePrecision function BielectronicMultipole(s,n1,n3,n2,n4,ll)&
       result( Multipole )
    !
    class(ClassBSpline), intent(in) :: s
    integer            , intent(in) :: n1,n2,n3,n4,ll

    DoublePrecision :: z12mi(0:2*s%Order),z12ma(0:2*s%Order) 
    DoublePrecision :: z34mi(0:2*s%Order),z34ma(0:2*s%Order)
    DoublePrecision :: zz12 (0:2*s%Order),zz34 (0:2*s%Order)
    DoublePrecision :: v1(0:2*(s%Order-1)),v2(0:2*(s%Order-1))
    DoublePrecision :: ain,w1,w2,x,d
    character       :: cas
    integer         :: ka,ntot,jm,nua,nub,nu,ma,mb,mc,md,mu,ima,jmi

    call CheckInitialization(s)

    Multipole = 0.d0

    zz12 =0.d0
    zz34 =0.d0
    z12mi=0.d0
    z12ma=0.d0
    z34mi=0.d0
    z34ma=0.d0

    ka=s%Order
    ntot=s%NBSplines

    jm=2*ka-2

    nua=n1-n2
    nub=n3-n4
    nu =n1-n3

    ma=n3
    mb=min(n4+ka,n1)
    if(n4+ka<=n1)then 
       cas="A"
       mc=max(n4+ka,n1)
       md=n2+ka
    elseif(n4<=n2)then
       cas="B"	  
       mc=max(n4+ka,n1)
       md=n2+ka
    else
       cas="C"
       mc=n2+ka
       md=n4+ka
    endif

    select case (cas)

    case("A")

       w1=0.d0
       do mu=0,mb-ma-1
          x=s%Grid(n3+mu-ka)
          d=s%Grid(n3+mu+1-ka)-x
          if(d<NOD_THRESHOLD) cycle
          if((n3+mu<ka).or.(n3+mu>ntot))cycle
          call coff(s%c(0,mu,n3),s%c(0,mu+nub,n4),ka-1,v1)
          w1=w1+relle(v1,x,d,jm,ll)
       enddo

       w2=0.d0
       do mu=0,md-mc-1
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v1)
          z12ma(mu)=relmu(v1,ll,x,d,jm)
          w2=w2+z12ma(mu)
       enddo

       ain=w1*w2

    case("B")

       do mu=0,mc-mb-1
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v1)
          call coff(s%c(0,mu+nu,n3),s%c(0,mu+nu+nub,n4),ka-1,v2)
          zz12(mu+nu)=rmarmi(v1,v2,ll,x,d,jm)
          zz34(mu+nu)=rmarmi(v2,v1,ll,x,d,jm)
       enddo

       do mu=0,mc-mb-2 
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v1)
          z12mi(mu+nu)=relle(v1,x,d,jm,ll)
       enddo

       do mu=0,md-mb-1
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v1)
          z12ma(nu+mu)=relmu(v1,ll,x,d,jm)
       enddo

       do mu=0,mc-ma-1 
          x=s%Grid(n3+mu-ka)
          d=s%Grid(n3+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n3+mu<ka).or.(n3+mu>ntot))cycle
          if((mc==md).and.(n3+mu+1.eq.md))cycle
          call coff(s%c(0,mu,n3),s%c(0,mu+nub,n4),ka-1,v1)
          z34mi(mu)=relle(v1,x,d,jm,ll)
       enddo

       do mu=nu+1,mc-ma-1
          x=s%Grid(n3+mu-ka)
          d=s%Grid(n3+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n3+mu<ka).or.(n3+mu>ntot))cycle
          call coff(s%c(0,mu,n3),s%c(0,mu+nub,n4),ka-1,v1)
          z34ma(mu)=relmu(v1,ll,x,d,jm)
       enddo

       ain=0.d0
       do ima=nu,nu+md-mb-1
          do jmi=0,min(ima-1,mc-ma-1)
             ain=ain+z34mi(jmi)*z12ma(ima)
          enddo
       enddo
       do ima=nu+1,mc-ma-1
          do jmi=nu,ima-1
             ain=ain+z34ma(ima)*z12mi(jmi)
          enddo
       enddo
       do ima=nu,mc-ma-1
          ain=ain+zz12(ima)+zz34(ima)
       enddo

    case("C")

       do mu=0,mc-mb-1
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v1)
          z12mi(nu+mu)=relle(v1,x,d,jm,ll)
       enddo

       do mu=0,mc-mb-1
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v1)
          z12ma(nu+mu)=relmu(v1,ll,x,d,jm)
       enddo

       do mu=nu+1,md-ma-1
          x=s%Grid(n3+mu-ka)
          d=s%Grid(n3+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n3+mu<ka).or.(n3+mu>ntot))cycle
          call coff(s%c(0,mu,n3),s%c(0,mu+nub,n4),ka-1,v1)
          z34ma(mu)=relmu(v1,ll,x,d,jm)
       enddo

       do mu=0,mc-ma-2
          x=s%Grid(n3+mu-ka)
          d=s%Grid(n3+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n3+mu<ka).or.(n3+mu>ntot))cycle
          call coff(s%c(0,mu,n3),s%c(0,mu+nub,n4),ka-1,v1)
          z34mi(mu)=relle(v1,x,d,jm,ll)
       enddo

       do mu=0,mc-mb-1
          x=s%Grid(n1+mu-ka)
          d=s%Grid(n1+mu+1-ka)-x
          if(d<NOD_THRESHOLD)cycle
          if((n1+mu<ka).or.(n1+mu>ntot))cycle
          call coff(s%c(0,mu+nu,n3),s%c(0,mu+nu+nub,n4),ka-1,v1)
          call coff(s%c(0,mu,n1),s%c(0,mu+nua,n2),ka-1,v2)
          zz12(nu+mu)=rmarmi(v1,v2,ll,x,d,jm)
          zz34(nu+mu)=rmarmi(v2,v1,ll,x,d,jm)
       enddo

       ain=0.d0
       do ima=nu,nu+mc-mb-1
          do jmi=0,ima-1
             ain=ain+z12ma(ima)*z34mi(jmi)
          enddo
       enddo
       do ima=nu+1,md-ma-1
          do jmi=nu,min(ima-1,mc-ma-1)
             ain=ain+z34ma(ima)*z12mi(jmi)
          enddo
       enddo
       ain=ain+sum(zz12(nu:nu+mc-mb-1)+zz34(nu:nu+mc-mb-1))

    end select

    Multipole = ain

    return
  end function BielectronicMultipole


  subroutine coff(a,b,k,c)
    integer        , intent(in) :: k
    DoublePrecision, intent(in) :: a(0:*), b(0:*)
    DoublePrecision, intent(out):: c(0:*)
    integer :: i,j
    c(0:2*k)=0.d0
    do i=0,k
       do j=0,k
          c(i+j)=c(i+j)+a(i)*b(j)
       enddo
    enddo
    return
  end subroutine coff


  !> Computes
  !! \f[ 
  !!    d^{l+1}\, \sum_{j=0}^l\, \sum_{i=0}^{j_m}\,
  !!             \frac{l}{j}\,
  !!             \left(\frac{x}{d}\right)^j\,
  !!             \frac{a_i}{i+l-j+1}
  !! \f]
  DoublePrecision function relle(a,x,d,jm,l)
    !
    DoublePrecision, intent(in) :: a(0:*)
    DoublePrecision, intent(in) :: x,d
    integer        , intent(in) :: jm,l
    !
    DoublePrecision :: y,yj,aa
    integer         :: i,j
    !
    relle=0.d0
    y=x/d
    do j=0,l
       if(j==0) then
          yj=1.d0
       else
          yj=(y**j)
       endif
       aa=0.d0
       do i=0,jm
          aa=aa+a(i)/(i+l+1-j)
       enddo
       relle=relle+aa*tar(l*(l+1)/2+j)*yj
    enddo
    relle=relle*(d**(l+1))
    return
  end function relle


  !> Computes
  !! \f[
  !! \delta \int \frac{dz}{(x+z \delta)^{l+1}\,\sum_k\left(a_k z^k\right)}
  !! \f]
  DoublePrecision function relmu(a,l,x,d,jm)
    DoublePrecision, intent(in) :: a(0:*)
    DoublePrecision, intent(in) :: x,d
    integer        , intent(in) :: l,jm
    DoublePrecision :: y,yy,pr
    integer         :: i,k
    yy=x/d
    relmu=0.d0
    do i=1,NGAUSS      
       y=Gauss_Points(i)
       pr=a(jm)
       do k=1,jm
          pr=pr*y+a(jm-k)
       enddo
       relmu=relmu+Gauss_Weight(i)*pr/((y+yy)**(l+1))
    enddo
    relmu=relmu/(d**l)
    return
  end function relmu


  !> Computes
  !! \f[
  !! d \int_0^1 dz \sum_i a_i z^i\,
  !!   \left(\frac{d}{x+z d}\right)^{l+1}
  !! \int_0^z dy \sum_j b_j y^j (x/d+y)**l
  !! \f]
  DoublePrecision function rmarmi(a,b,l,x,d,jm)
    DoublePrecision, intent(in) :: a(0:*), b(0:*)
    DoublePrecision, intent(in) :: x,d
    integer        , intent(in) :: l,jm
    DoublePrecision :: xd,z,w,dxz,pra,aa,seco,xdlj,bb,prb
    integer         :: i,j,k,jmi,ind,nfx
    xd=x/d
    rmarmi=0.d0
    do i=1,NGAUSS
       z=Gauss_Points(i)
       w=Gauss_Weight(i)*d
       dxz=d/(x+z*d)
       pra=a(jm)
       do k=1,jm
          pra=pra*z+a(jm-k)
       enddo
       pra=pra*z
       aa=pra*(dxz**(l+1))  
       ind=l*(l+1)/2
       seco=0.d0
       jmi=0
       if(xd<1.d-9)jmi=l
       do j=jmi,l
          if(j==l)then
             xdlj=1.d0
          else
             xdlj=xd**(l-j)
          endif
          nfx=jm+j+1   
          bb=tar(ind+j)*xdlj*(z**j)
          prb=b(jm)/nfx
          do k=1,jm
             prb=prb*z+b(jm-k)/(nfx-k)
          enddo
          seco=seco+bb*prb
       enddo
       rmarmi=rmarmi+aa*seco*w
    enddo
    return
  end function rmarmi


  !.. Set of do-nothing subroutines for the 
  !   preconditioning of the B-spline basis
  !..
  logical function PreconditionerIsAvailable(BSSet, L, dir ) result(Available)
    class(ClassBSpline), intent(in) :: BSSet
    integer            , intent(in) :: L
    character(len=*)   , intent(in) :: dir
    Available = L >= 0
  end function PreconditionerIsAvailable
  !
  subroutine ComputePreconditioner(BSSet, L ) 
    class(ClassBSpline), intent(inout) :: BSSet
    integer            , intent(in)    :: L
  end subroutine ComputePreconditioner
  !
  subroutine SavePreconditioner(BSSet, L, dir )
    class(ClassBSpline), intent(in) :: BSSet
    integer            , intent(in) :: L
    character(len=*)   , intent(in) :: dir
  end subroutine SavePreconditioner
  !
  subroutine LoadPreconditioner(BSSet, L, dir )
    class(ClassBSpline), intent(inout) :: BSSet
    integer            , intent(in)    :: L
    character(len=*)   , intent(in)    :: dir
  end subroutine LoadPreconditioner


!!$  !> (was cosmu)
!!$  !! Computes all non-vanishing multipoles l
!!$  !! between **NORMALIZED** B-splines,
!!$  !! \f[
!!$  !!    \bar{B}_i(r) \equiv B_i(r) / \|B_i\|_{L^2}.
!!$  !! \f]
!!$  subroutine ComputeSplineMultipole(s,ssss,ll)
!!$    type(ClassBSpline), intent(in)  :: s
!!$    DoublePrecision   , intent(out) :: ssss(0:)
!!$    integer           , intent(in)  :: ll
!!$    integer :: n1,n2,n3,n4,ind,n4ma
!!$    call CheckInitialization(s)
!!$    ssss=0.d0
!!$    ind=0
!!$    do n1=1,s%NBSplines
!!$       do n3=max(1,n1-s%Order+1),n1
!!$          do n2=1,n1
!!$             n4ma=n2
!!$             if(n2==n1)n4ma=n3
!!$             do n4=max(1,n2-s%Order+1),n4ma
!!$                ind=ind+1
!!$                ssss(ind)=s%Multipole(n1,n2,n3,n4,ll)*s%f(n1)*s%f(n2)*s%f(n3)*s%f(n4)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$    return
!!$  end subroutine ComputeSplineMultipole

  ! }}}


  ! {{{ Interpolation Routines (currently inactive)

!!$  !*** Routines of dubious current value follow which should be
!!$  !*** seriously considered for complete deletion. If not, they
!!$  !*** should probably be rewritten.
!!$
!!$  !> Computes the matrix of charactersitic functions \f$\phi_i(x)\f$ 
!!$  !> of the grid points \f$\{g_j\}_{j=0,1,\ldots,Nnodes-1}\f$
!!$  !> \f[
!!$  !>  \phi_i(g_j)=\delta_{ij}
!!$  !> \f]
!!$  !> for a set of CUBIC splines (order=4), under the constraint
!!$  !> that the third derivative at the two boundaries is unitary.
!!$  integer function Init_CarFun(s)
!!$    type(ClassBSpline), intent(inout) :: s
!!$    DoublePrecision, allocatable :: A(:,:)
!!$    DoublePrecision :: x
!!$    integer :: i,j
!!$    Init_CarFun=0
!!$    call CheckInitialization(s)
!!$    if(s%Order/=4)then
!!$       Init_CarFun=-1
!!$       return
!!$    endif
!!$    allocate(A(s%NBSplines,s%NBSplines))
!!$    A=0.d0
!!$    do i=1,s%NNodes
!!$       x=s%Grid(i-1)
!!$       do j=1,s%NBSplines
!!$          A(i,j)=s%Eval(x,j)
!!$       enddo
!!$    enddo
!!$    x=s%Grid(0)
!!$    do j=1,s%NBSplines
!!$       A(s%NNodes+1,j)=s%Eval(x,j,3)
!!$    enddo
!!$    x=s%Grid(s%NNodes-1)
!!$    do j=1,s%NBSplines
!!$       A(s%NNodes+2,j)=s%Eval(x,j,3)
!!$    enddo
!!$    call DLINRGspline(s%NBSplines,A,s%NBSplines,A,s%NBSplines)
!!$    allocate(s%carfun(s%NBSplines,s%NNodes))
!!$    allocate(s%carfunt(s%NNodes,s%NBSplines))
!!$    s%carfun(1:s%NBSplines,1:s%NNodes)=A(1:s%NBSplines,1:s%NNodes)
!!$    s%carfunt(1:s%NNodes,1:s%NBSplines)=transpose(A(1:s%NBSplines,1:s%NNodes))
!!$    deallocate(A)
!!$    s%CARFUNINIT=.TRUE.
!!$    return
!!$  end function Init_CarFun
!!$
!!$
!!$  !.. As Init_CarFun, with the difference that for the first node it does not 
!!$  !   compute the characteristic function, but one with unitary derivative
!!$  !..
!!$  integer function Init_CarFun2(s)
!!$    type(ClassBSpline), intent(inout) :: s
!!$    DoublePrecision, allocatable :: A(:,:)
!!$    DoublePrecision :: x
!!$    integer :: i,j
!!$    Init_CarFun2=0
!!$    call CheckInitialization(s)
!!$    if(s%Order/=4)then
!!$       Init_CarFun2=-1
!!$       return
!!$    endif
!!$    allocate(A(s%NBSplines,s%NBSplines))
!!$    A=0.d0
!!$    do i=1,s%NNodes-1
!!$       do j=1,s%NBSplines
!!$          x=s%Grid(i)
!!$          A(i,j)=s%Eval(x,j)
!!$       enddo
!!$    enddo
!!$    x=s%Grid(0)
!!$    do j=1,s%NBSplines
!!$       A(s%NNodes,j)=s%Eval(x,j,3)
!!$    enddo
!!$    do j=1,s%NBSplines
!!$       A(s%NNodes+1,j)=s%Eval(x,j,2)
!!$    enddo
!!$    x=s%Grid(s%NNodes-1)
!!$    do j=1,s%NBSplines
!!$       A(s%NNodes+2,j)=s%Eval(x,j,3)
!!$    enddo
!!$    call DLINRGspline(s%NBSplines,A,s%NBSplines,A,s%NBSplines)
!!$    allocate(s%carfun(s%NBSplines,s%NNodes-1))
!!$    allocate(s%carfunt(s%NNodes-1,s%NBSplines))
!!$    s%carfun(1:s%NBSplines,1:s%NNodes-1)=A(1:s%NBSplines,1:s%NNodes-1)
!!$    s%carfunt(1:s%NNodes-1,1:s%NBSplines)=transpose(A(1:s%NBSplines,1:s%NNodes-1))
!!$    deallocate(A)
!!$    s%CARFUNINIT=.TRUE.
!!$    return
!!$  end function Init_CarFun2
!!$
!!$
!!$  !> Interpolate a function on the basis of its value on the
!!$  !> grid of nodes of a pre-defined cubic BSpline set.
!!$  DoublePrecision function Evaluate_Function_ongrid(x,fc,s)
!!$    implicit none
!!$    DoublePrecision  , intent(in) :: x
!!$    DoublePrecision  , intent(in) :: fc(:)
!!$    type(ClassBSpline), intent(inout) :: s
!!$    call CheckInitialization(s)
!!$    if(s%Order/=4)call Assert("Evaluate_Function_ongrid available for cubic splines only")
!!$    if(.not.s%CARFUNINIT)then
!!$       if(Init_CarFun(s)/=0)call Assert("Initialization of characteristic functions failed")
!!$    endif
!!$    if(.not.allocated(s%vec))allocate(s%vec(s%NBSplines))
!!$    s%vec=matmul(s%carfun,fc(1:s%NNodes))
!!$    Evaluate_Function_ongrid=s%FEval(x,s%vec(1:s%NBSplines))
!!$    return
!!$  end function Evaluate_Function_ongrid
!!$
!!$
!!$  DoublePrecision function Evaluate_Function_ongrid2(x,fc,s)
!!$    !Valida esclusivamente per spline cubiche (s%Order=4)
!!$    !Prende in ingresso il campionamento di una funzione
!!$    !in corrispondenza dei nodi della griglia:
!!$    ! fc(i)=f(xi)
!!$    implicit none
!!$    DoublePrecision  , intent(in) :: x
!!$    DoublePrecision  , intent(in) :: fc(:)
!!$    type(ClassBSpline), intent(inout) :: s
!!$    call CheckInitialization(s)
!!$    if(s%Order/=4)call Assert("Evaluate_Function_ongrid available for cubic splines only")
!!$    if(.not.s%CARFUNINIT)then
!!$       if(Init_CarFun2(s)/=0)call Assert("Initialization of characteristic functions failed")
!!$    endif
!!$    if(.not.allocated(s%vec))allocate(s%vec(s%NBSplines))
!!$    s%vec=matmul(s%carfun,fc(1:s%NNodes-1))
!!$    Evaluate_Function_ongrid2=s%FEval(x,s%vec(1:s%NBSplines))
!!$    return
!!$  end function Evaluate_Function_ongrid2
!!$  !
!!$
!!$  subroutine Evaluate_G0matrix(G0,E0,s)
!!$    implicit none
!!$    DoublePrecision   , intent(inout) :: G0(:,:)
!!$    DoublePrecision   , intent(in)    :: E0
!!$    type(ClassBSpline), intent(inout) :: s
!!$    !
!!$    DoublePrecision, allocatable :: A(:,:)
!!$    integer :: i,j
!!$
!!$    call CheckInitialization(s)
!!$    if(s%Order/=4)call Assert("Evaluate_G0matrix available for cubic splines only")
!!$    if(.not.s%CARFUNINIT)then
!!$       if(Init_CarFun(s)/=0)call Assert("Initialization of Characteristic functions failed")
!!$    endif
!!$
!!$    allocate(A(s%NBSplines,s%NBSplines))
!!$    A=0.d0
!!$    do i = 1, s%NBSplines
!!$       do j = max(1,i-s%Order+1),i
!!$          A(i,j)=s%CauchyIntegral(i,j,E0)
!!$          A(j,i)=A(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    G0=0.d0
!!$    !.. This is certainly a cause for inefficiency
!!$    G0=matmul(s%carfunt,matmul(A,s%carfun))
!!$    deallocate(A)
!!$
!!$    return
!!$  end subroutine Evaluate_G0matrix
!!$
!!$
!!$  subroutine dlinrgspline(n, a, lda, ainv, ldainv)
!!$    implicit none
!!$    integer, intent(in)::n,lda,ldainv
!!$    DoublePrecision:: a(lda,*),ainv(ldainv,*)
!!$    DoublePrecision,allocatable::wk(:)
!!$    integer ::ipiv(n)
!!$    integer:: info,lw
!!$    !
!!$    ainv(1:n,1:n)=a(1:n,1:n)
!!$    call dgetrf(ldainv,n,ainv,ldainv,ipiv,info)
!!$    !
!!$    allocate(wk(n))
!!$    lw=-1
!!$    call dgetri(n,ainv,ldainv,ipiv,wk,lw,info)
!!$    lw=idint(wk(1))+1
!!$    deallocate(wk)
!!$    !
!!$    allocate(wk(lw))
!!$    call dgetri(n,ainv,ldainv,ipiv,wk,lw,info)
!!$    deallocate(wk)
!!$    !
!!$    return
!!$  end subroutine dlinrgspline

! }}}


  !> Creates in a directory all the configuration
  !!    files that define the current basis.
  subroutine SaveConfig(self, ConfigurationFile, Dir)
    !
    Class(ClassBSpline), intent(in)  :: self
    character(len=*)   , intent(in)  :: ConfigurationFile
    character(len=*)   , intent(in)  :: Dir
    !
    call SYSTEM("mkdir -p "//trim(Dir))
    call SYSTEM("cp "//trim(ConfigurationFile)//" "//trim(Dir))
    !
  end subroutine SaveConfig



end module ModuleBSpline

! }}}
