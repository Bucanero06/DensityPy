 
! {{{ Detailed description

!> \file
!!
!! Defines a group of representation for some basic monoelectronic operators and the methods required to compute, load and store them.

! }}}
Module ModuleBasicMatrixElements

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleRepresentation
  use ModuleBasis
  use ModuleMatrix

  implicit none
  
  private



  !> Defines a group of representations of basic operators:
  !! S : overlap     
  !!     \f$ \mathbf{S}_{ij}=\langle\chi_i|\chi_j\rangle \f$
  !! D : derivative  
  !!     \f$ \mathbf{D}_{ij}=\langle\chi_i|\frac{d}{dr}|\chi_j\rangle \f$
  !! K : Radial kinetic Energy
  !!     \f$ \mathbf{K}_{ij}=-\frac{1}{2}\langle\chi_i|\frac{d^2}{dr^2}|\chi_j\rangle \f$
  !! L : Centrifugal potential
  !!     \f$ \mathbf{L}_{ij}=\langle\chi_i|r^{-2}|\chi_j\rangle \f$
  !! R : Radius 
  !!     \f$ \mathbf{R}_{ij}=\langle\chi_i|r|\chi_j\rangle \f$
  !! RD: Radius Times Derivative
  !!     \f$ \mathbf{RD}_{ij}=\langle\chi_i|r\frac{d}{dr}|\chi_j\rangle \f$
  !! C : Coulomb Potential
  !!     \f$\mathbf{C}_{ij}=\langle\chi_i|r^{-1}|\chi_j\rangle \f$
  !! EffPot : Effective Potential
  !!     \f$\mathbf{V_{eff}}_{ij}=\langle\chi_i|V_{eff}|\chi_j\rangle \f$
  !!
  type, public :: ClassBasicMatrixElements
     !> Class of representations associated to the  \f$ \mathbf{S}_{ij}\f$ matrix.
     type(ClassRepresentation) :: S
     !> Class of representations associated to the  \f$ \mathbf{D}_{ij}\f$ matrix.
     type(ClassRepresentation) :: D
     !> Class of representations associated to the  \f$ \mathbf{K}_{ij}\f$ matrix.
     type(ClassRepresentation) :: K
     !> Class of representations associated to the  \f$ \mathbf{L}_{ij}\f$ matrix.
     type(ClassRepresentation) :: L
     !> Class of representations associated to the  \f$ \mathbf{R}_{ij}\f$ matrix.
     type(ClassRepresentation) :: R
     !> Class of representations associated to the  \f$ \mathbf{RD}_{ij}\f$ matrix.
     type(ClassRepresentation) :: RD
     !> Class of representations associated to the  \f$ \mathbf{C}_{ij}\f$ matrix.
     type(ClassRepresentation) :: C
     !> Class of representations associated to the  \f$ \mathbf{V_{eff}}_{ij}\f$ matrix.
     type(ClassRepresentation) :: Veff
     !> Class of representations associated to the  \f$ \mathbf{Multipole_{eff}}_{ij}\f$ matrix.
     type(ClassRepresentation), allocatable :: Multipole(:)

   contains
     
     !> Initializes the ClassBasicMatrixElements.
     procedure :: Init

     !> Computes the radial Matrix Elements of several local operators 
     !! using the choosed basis kind. More precisely computes 
     !! - the Overlap Matrix (S), 
     !! - the matrices asociated with: 
     !! -- the Nabla Operator(D), 
     !! -- the Kinetic Energy Operator(K), 
     !! -- the \f$ 1/r^2\f$ operator (L), 
     !! -- the r operator (R), 
     !! -- the r*nabla operator (RD), and 
     !! -- the 1/r operator (C).
     !! -- the Veff operator (Veff).
     procedure :: Compute

     !> Computes the matrix of the operator \f$r^(-L)\f$
!!$     procedure :: ComputeRToMinusLMat => BasicMatElemComputeRToMinusLMat

     !> Saves the computed local operator's matrix elements.
     procedure :: Save => SaveBasicMatrixElements

     !> Reads from a file the local operator's matrix elements.
     procedure :: Read => ReadBasicMatrixElements

     !> Fetches the Metric matrix 
     !! \f[
     !!    \mathbf{S}_{ij}
     !! \f] 
     !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
     procedure :: FetchMetricMatrix

     !> Fetches the Multipole matrix 
     !! \f[
     !!    \mathbf{\frac{1}{r^{l+1}}}_{ij}
     !! \f] 
     !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
     procedure :: FetchMultipoleMatrix

     !> Fetches the metric matrix defined on the conditioned
     !! basis for a specified angular momentum
     procedure :: FetchMetric
     !> Fetches the Preconditioned Metric matrix.   
     procedure :: FetchPreconditionedMetric
     !> Fetches the Preconditioned Continuous Metric matrix.
     procedure :: FetchPreconditionedContinuousMetric

     !> Fetches the Length matrix 
     !! \f[
     !!    \mathbf{R}_{ij}
     !! \f] 
     !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
     procedure :: FetchLengthMatrix

     !> Fetches the Velocity matrix 
     !! \f[
     !!    \mathbf{V}_{ij}=-\sqrt{l+1}(\mathbf{D}_{ij}+(l+1)\mathbf{C}_{ij})
     !! \f] 
     !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
     procedure :: FetchVelocityMatrix

     procedure :: FetchGeneralVelocityMatrix
     procedure :: FetchGeneralLengthMatrix

     !> Fetches the Acceleration matrix 
     !! \f[
     !!    \mathbf{L}_{ij}
     !! \f] 
     !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
     procedure :: FetchAccelerationMatrix

     !> Fetches the Hamiltonian matrix 
     !! \f[
     !!    \mathbf{H}_{ij}=\mathbf{K}_{ij}-Z\mathbf{C}_{ij}+\frac{1}{2}l(l+1)\mathbf{L}_{ij}
     !! \f] 
     !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
     procedure :: FetchHamiltonianMatrix

     procedure :: FetchKineticMatrix

     procedure :: FetchSpecialHamiltonianMatrix

     !> Fetches the metric matrix defined on the pre-conditioned
     !! continuous basis for a specified angular momentum
     procedure :: FetchPreconditionedHamiltonian
     procedure :: FetchPreconditionedContinuousHamiltonian
     !> Fetches the Hamiltonian matrix in real or complex format and in some specified representation level.
     generic :: FetchHamiltonian => &
          FetchHamiltonianDouble,&
          FetchHamiltonianComplex

     !> Return .TRUE. if a file supposed to contained previously computed 
     !! monolelectronic matrix elements is readable, and its content is compatible 
     !! with the current [radial basis](@ref ModuleBasis.f90), return .FALSE. otherwise.
     procedure :: Saved

     procedure, private :: FetchHamiltonianDouble
     procedure, private :: FetchHamiltonianComplex
!!$     procedure, private :: BasicMatElemComputeRToMinusLMat

  end type ClassBasicMatrixElements

  character(len=*), parameter, private :: BASIC_MATRIX_ELEMENTS_FILE = "matint"
  character(len=*), parameter, public :: DEFAULT_SAVE_FORMAT = "UNFORMATTED"

  public :: BasicMatElemComputeAndSaveRToMinusNMat
  public :: BasicMatElemComputeRToMinusNMat


contains


  !> Return .TRUE. if FileName is readable, and its content is compatible 
  !! with the current Basis associated to MESet, return .FALSE. otherwise.
  logical function Saved( MESet, DirName )
    !> Class of basic matrix elements.
    class(ClassBasicMatrixElements), intent(in) :: MESet
    !> Name of the file.
    character(len=*)               , intent(in) :: DirName
    !
    character(len=*), parameter :: Form = DEFAULT_SAVE_FORMAT
    character(len=:), allocatable :: FileName
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    Saved=.FALSE.
    !
    allocate( FileName, source = trim(adjustl(DirName))//BASIC_MATRIX_ELEMENTS_FILE )
    open(&
         Newunit =  uid      , &
         File    =  FileName , &
         Status  = "old"     , &
         Action  = "read"    , &
         Form    =  Form     , &
         iostat  =  iostat   , &
         iomsg   =  iomsg      &
         )
    if(iostat/=0)return
    !
    Saved = MESet.S.Basis.MatchFingerprintOnUnit( uid )
    !
    close(uid)
    !
  end function Saved


  !> Initializes the ClassBasicMatrixElements.
  subroutine Init( MESet, Basis )
    !
    class(ClassBasicMatrixElements), intent(inout) :: MESet
    class(ClassBasis)              , intent(in)    :: Basis
    !
    integer :: i, NumMultipoles
    ! Maximum angular momentum that can be represented by the Gaussian functions subset.
    integer ::Lmax
    !
    Lmax = Basis.GetGaussianLMax( )
    NumMultipoles = 2*Lmax
    !
    call MESet%S%Init( Basis )
    call MESet%D%Init( Basis )
    call MESet%K%Init( Basis )
    call MESet%L%Init( Basis )
    call MESet%R%Init( Basis )
    call MESet%C%Init( Basis )
    call MESet%RD%Init( Basis )
    call MESet%Veff%Init( Basis )
    !
    allocate( MESet%Multipole(NumMultipoles) )
    do i = 1, NumMultipoles
       call MESet%Multipole(i)%Init( Basis )
    end do
    !
  end subroutine Init


  !> Computes the radial Matrix Elements of several local operators using the 
  !! choosed basis kind. More precisely computes the Overlap Matrix (S), the 
  !! matrices asociated with: the Nabla Operator(D), the Kinetic Energy 
  !! Operator(K), the \f$ 1/r^2\f$ operator (L), the r operator (R), the 
  !! r*nabla operator (RD), and the 1/r operator (C).
  !!
  subroutine Compute(MESet, FileVeff) 
    !
    class(ClassBasicMatrixElements), intent(inout) :: MESet
    character(len=*), optional,      intent(in)    :: FileVeff
    !
    procedure(D2DFun), pointer :: fptr
    ! Nuclear charge.
    integer :: Z
    ! Number of electrons.
    integer :: N
    ! power of 'r' in Veff.
    integer, allocatable :: b(:)
    ! coefficients in the exponents in Veff.
    real(kind(1d0)), allocatable :: d(:)
    ! coefficients in the summation in Veff.
    real(kind(1d0)), allocatable :: c(:)
    integer :: i, NumMultipoles, MultipoleOrder
    !
    fptr => CentPot; call MESet%L%Integral( fptr )
    fptr => CoulPot; call MESet%C%Integral( fptr )
    fptr => Unity  ; call MESet%S%Integral( fptr ) 
    fptr => Radius ; call MESet%R%Integral( fptr )
    fptr => Unity  ; call MESet%D%Integral( fptr,  KetDerivativeOrder = 1 )
    fptr => Radius ; call MESet%RD%Integral( fptr, KetDerivativeOrder = 1 )
    fptr => Unity  ; call MESet%K%Integral( fptr,  KetDerivativeOrder = 2 )
    call MESet%K%Multiply(-0.5d0) 
    !
    if ( present(FileVeff) ) then
       call LoadEffectivePotential( &
            FileVeff, &
            Z, N, &
            b, d, c )
       fptr => Veffective; call MESet%Veff%Integral( fptr )
    else
       ! The effective portential will be equal to Coulomb potential.
       fptr => CoulPot; call MESet%Veff%Integral( fptr )
    end if
    !
    !
    NumMultipoles = size(MESet%Multipole)
    do i = 1, NumMultipoles
       MultipoleOrder = i + 1
       fptr => MultipolePot; call MESet%Multipole(i)%Integral( fptr )
    end do
   !
  contains
    !
    Pure DoublePrecision function CentPot( r ) result( V )
      DoublePrecision, intent(in) :: r
      DoublePrecision, parameter :: MINR = sqrt(tiny(1.d0))
      V = 0.d0; if(abs(r)>MINR) V = 1.d0 / ( r * r )
    end function CentPot
    Pure DoublePrecision function CoulPot( r ) result( V )
      DoublePrecision, intent(in) :: r
      V = 0.d0; if(abs(r)>tiny(1.d0)) V = 1.d0 / r
    end function CoulPot
    Pure DoublePrecision function Unity( r )
      DoublePrecision, intent(in) :: r
      Unity = 1.d0
    end function Unity
    Pure DoublePrecision function Radius( r )
      DoublePrecision, intent(in) :: r
      Radius = r
    end function Radius
    Pure DoublePrecision function Veffective( r )
      DoublePrecision, intent(in) :: r
      integer :: i
      Veffective = 0.d0
      do i = 1, size(b)
         Veffective = Veffective + c(i)*r**b(i)*exp(-d(i)*r)
      end do
      Veffective = -1.d0/r*( dble(Z)-dble(N)+1.d0+(dble(N)-1.d0)*Veffective )
    end function Veffective
    Pure DoublePrecision function MultipolePot( r ) result( V )
      DoublePrecision, intent(in) :: r
      DoublePrecision, parameter :: MINR = sqrt(tiny(1.d0))
      V = 0.d0; if(abs(r)>MINR) V = 1.d0 / ( r ** MultipoleOrder )
    end function MultipolePot
    !
  end subroutine Compute




  subroutine LoadEffectivePotential( &
       FileVeff, &
       Z, N, &
       b, d, c )
    !
    use ModuleParameterList
    !
    character(len=*),             intent(in)  :: FileVeff
    integer,                      intent(out) :: Z, N
    integer, allocatable,         intent(out) :: b(:)
    real(kind(1d0)), allocatable, intent(out) :: d(:), c(:)
    !
    integer :: uid, iostat, NumSum, i
    character(len=IOMSG_LENGTH) :: iomsg, File
    type(ClassParameterList) :: List
    !
    call List%Add( "Z", 1 , "required" )
    call List%Add( "N", 1 , "required" )
    call List%Add( "ParametersFile", "File" , "required" )
    call List%Parse( FileVeff )
    call List%Get( "Z", Z  )
    call List%Get( "N", N  )
    call List%Get( "ParametersFile", File  )
    !
    open(Newunit =  uid     , &
         File    =  trim(adjustl(File)), &
         Status  = "old"    , &
         Action  = "read"   , &
         Form    =  "formatted"    , &
         iostat  =  iostat  , &
         iomsg   =  iomsg   )
    if(iostat/=0)call Assert(iomsg)
    !
    read(uid,*) NumSum
    allocate( b(NumSum) )
    allocate( d(NumSum) )
    allocate( c(NumSum) )
    !
    do i = 1, NumSum
       read(uid,*) b(i), d(i), c(i)
    end do
    !
    close( uid )
    !
  end subroutine LoadEffectivePotential


  !> Computes the matrix of the operator \f$r^(-N)\f$
  subroutine BasicMatElemComputeAndSaveRToMinusNMat( Basis, N, File, BraDer, KetDer, LowerLim, UpperLim, BreakPoint )
    !
!!$    class(ClassBasicMatrixElements), intent(in) :: MESet
    class(ClassBasis)              , intent(in)    :: Basis
    integer,                         intent(in)    :: N
    character(len=*),                intent(in)    :: File
    integer,         optional,       intent(in)    :: BraDer
    integer,         optional,       intent(in)    :: KetDer
    real(kind(1d0)), optional,       intent(in)    :: LowerLim
    real(kind(1d0)), optional,       intent(in)    :: UpperLim
    real(kind(1d0)), optional,       intent(in)    :: BreakPoint
    !
    type(ClassMatrix) :: Mat
    type(ClassRepresentation) :: RMinusL
    procedure(D2DFun), pointer :: fptr
    real(kind(1d0)) :: a, b, bp
    integer :: BD, KD
    !
    if ( present(BraDer) ) then
       BD = BraDer
    else
       BD = 0
    end if
    !
    if ( present(KetDer) ) then
       KD = KetDer
    else
       KD = 0
    end if
    !
    if ( present(LowerLim) ) then
       a = LowerLim
    else
       a = -huge(1.d0)
    end if
    !
    if ( present(UpperLim) ) then
       b = UpperLim
    else
       b = huge(1.d0)
    end if
    !
    if ( present(BreakPoint) ) then
       bp = BreakPoint
    else
       bp = a + epsilon(1.d0)
    end if
    !
    !
    call RMinusL.Init( Basis )
    !
    fptr => Pot
!!$    call RMinusL.Integral( fptr )
    call RMinusL.Integral( fptr, BD, KD, a, b, bp )
    !
    Mat = RMinusL.Matrix
    !
    call Mat.Write( File )
    !
  contains
    !
    Pure DoublePrecision function Pot( r ) result( V )
      DoublePrecision, intent(in) :: r
      DoublePrecision, parameter :: MINR = sqrt(tiny(1.d0))
      V = 0.d0; if(abs(r)>MINR) V = 1.d0 / ( r ** N )
    end function Pot 
    !
  end subroutine BasicMatElemComputeAndSaveRToMinusNMat



  subroutine BasicMatElemComputeRToMinusNMat( Basis, N, Array, BraDer, KetDer, LowerLim, UpperLim, BreakPoint )
    !
!!$    class(ClassBasicMatrixElements), intent(in) :: MESet
    class(ClassBasis)              , intent(in)    :: Basis
    integer,                         intent(in)    :: N
    class(ClassMatrix),              intent(out)   :: Array
    integer,         optional,       intent(in)    :: BraDer
    integer,         optional,       intent(in)    :: KetDer
    real(kind(1d0)), optional,       intent(in)    :: LowerLim
    real(kind(1d0)), optional,       intent(in)    :: UpperLim
    real(kind(1d0)), optional,       intent(in)    :: BreakPoint
    !
    type(ClassMatrix) :: Mat
    type(ClassRepresentation) :: RMinusL
    procedure(D2DFun), pointer :: fptr
    real(kind(1d0)) :: a, b, bp
    integer :: BD, KD
    !
    if ( present(BraDer) ) then
       BD = BraDer
    else
       BD = 0
    end if
    !
    if ( present(KetDer) ) then
       KD = KetDer
    else
       KD = 0
    end if
    !
    if ( present(LowerLim) ) then
       a = LowerLim
    else
       a = -huge(1.d0)
    end if
    !
    if ( present(UpperLim) ) then
       b = UpperLim
    else
       b = huge(1.d0)
    end if
    !
    if ( present(BreakPoint) ) then
       bp = BreakPoint
    else
       bp = a + epsilon(1.d0)
    end if
    !
    !
    call RMinusL.Init( Basis )
    !
    fptr => Pot
!!$    call RMinusL.Integral( fptr )
    call RMinusL.Integral( fptr, BD, KD, a, b, bp )
    !
    Mat = RMinusL.Matrix
    !
    Array = Mat
    !
  contains
    !
    Pure DoublePrecision function Pot( r ) result( V )
      DoublePrecision, intent(in) :: r
      DoublePrecision, parameter :: MINR = sqrt(tiny(1.d0))
      integer :: i
      V = 0.d0
      if(abs(r)>MINR) then
         V = 1.d0 / ( r ** N )
!!$         V = 1.d0
!!$         do i = 1, N
!!$            V = V/r  
!!$         end do
      end if
    end function Pot
    !
  end subroutine BasicMatElemComputeRToMinusNMat





  !> Saves the computed local operator's matrix elements.
  subroutine SaveBasicMatrixElements ( MESet, DirName ) 
    !
    class(ClassBasicMatrixElements), intent(in) :: MESet
    character(len=*),                intent(in) :: DirName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=*), parameter :: Form = DEFAULT_SAVE_FORMAT
    integer :: NumMultipoles, i
    character(len=:), allocatable :: FileName
    !
    allocate( FileName, source = trim(adjustl(DirName))//BASIC_MATRIX_ELEMENTS_FILE )
    NumMultipoles = size(MESet%Multipole)
    !
    open(Newunit =  uid     , &
         File    =  FileName, &
         Status  = "unknown", &
         Action  = "write"  , &
         Form    =  Form    , &
         iostat  =  iostat  , &
         iomsg   =  iomsg   )
    if(iostat/=0)call Assert(iomsg)
    !
    call MESet%S%Write( uid )
    call MESet%L%Write( uid )
    call MESet%C%Write( uid )
    call MESet%R%Write( uid )
    call MESet%D%Write( uid )
    call MESet%K%Write( uid )
    call MESet%RD%Write( uid )
    call MESet%Veff%Write( uid )
    do i = 1, NumMultipoles
       call MESet%Multipole(i)%Write( uid )
    end do
    !
    close(uid)
    !
  end subroutine SaveBasicMatrixElements


  !> Reads from a file the local operator's matrix elements.
  subroutine ReadBasicMatrixElements ( MESet, DirName ) 
    !
    class(ClassBasicMatrixElements), intent(inout) :: MESet
    character(len=*),                intent(in)    :: DirName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=*), parameter :: Form = DEFAULT_SAVE_FORMAT
    integer :: i, NumMultipoles
    character(len=:), allocatable :: FileName
    !
    allocate( FileName, source = trim(adjustl(DirName))//BASIC_MATRIX_ELEMENTS_FILE )
    NumMultipoles = size(MESet%Multipole)
    !
    open(Newunit =  uid     , &
         File    =  FileName, &
         Status  = "old"    , &
         Action  = "read"   , &
         Form    =  Form    , &
         iostat  =  iostat  , &
         iomsg   =  iomsg   )
    if(iostat/=0)call Assert(iomsg)
    !
    call MESet.S.Read( uid )
    call MESet.L.Read( uid )
    call MESet.C.Read( uid )
    call MESet.R.Read( uid )
    call MESet.D.Read( uid )
    call MESet.K.Read( uid )
    call MESet.RD.Read( uid )
    call MESet.Veff.Read( uid )
    do i = 1, NumMultipoles
       call MESet.Multipole(i).Read( uid )
    end do
    !
    close(uid)
    !
  end subroutine ReadBasicMatrixElements


  !> Fetches the Full Metric matrix in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchMetricMatrix( MESet, Metric )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Metric
    Metric = MESet.S.Matrix
  end subroutine FetchMetricMatrix


  !> Fetches the Full Multipole matrix /f$ \frac{1}{r^{l+1}}/f$ in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchMultipoleMatrix( MESet, AngularMomentum, Multipole )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    integer,                         intent(in)  :: AngularMomentum
    type(ClassMatrix)              , intent(out) :: Multipole
    Multipole = MESet.Multipole(AngularMomentum).Matrix
  end subroutine FetchMultipoleMatrix


  !> Fetches the Preconditioned Metric matrix in form of 
  !! [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchPreconditionedMetric( MESet, Metric, L )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Metric
    integer                        , intent(in)  :: L
    call MESet.S.GetPreconditionedMatrix( Metric, L, L )
  end subroutine FetchPreconditionedMetric

  !> Fetches the Preconditioned Continuous Metric matrix 
  !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchPreconditionedContinuousMetric( self, Metric, L )
    class(ClassBasicMatrixElements), intent(in)  :: self
    type(ClassMatrix)              , intent(out) :: Metric
    integer                        , intent(in)  :: L
    call self.S.GetPreconditionedContinuousMatrix( Metric, L, L )
  end subroutine FetchPreconditionedContinuousMetric

  
  subroutine FetchMetric( self, Metric, L, RepLevel )
    class(ClassBasicMatrixElements), intent(in)  :: self
    type(ClassMatrix)              , intent(out) :: Metric
    integer                        , intent(in)  :: L
    integer                        , intent(in)  :: RepLevel
    call self.S.GetMatrix( Metric, L, L, RepLevel, RepLevel )
  end subroutine FetchMetric


  !> Fetches the Length matrix \f$\mathbf{R}_{ij}\f$ in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchLengthMatrix( MESet, Length, Param1 )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Length
    DoublePrecision,                 intent(in)  :: Param1
    !
    Length = MESet.R.Matrix
    call Length.Multiply( Param1 )
  end subroutine FetchLengthMatrix



  !> Fetches the Velocity matrix 
  !! \f[
  !!    \mathbf{V}_{ij} = -\sqrt{l+1}\left[ \mathbf{D}_{ij} + (l+1) \mathbf{C}_{ij} \right]
  !! \f]
  !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchVelocityMatrix( MESet, Velocity, Param1, Param2 )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Velocity
    DoublePrecision,                 intent(in)  :: Param1, Param2
    !
    type(ClassMatrix) :: R_Derivative
    type(ClassMatrix) :: CoulombPotential
    !
    R_Derivative     = MESet.D.Matrix
    CoulombPotential = MESet.C.Matrix
    !
    call R_Derivative.Multiply( Param1 )
    call CoulombPotential.Multiply( Param2 )
    !
    Velocity = R_Derivative
    call Velocity.Add( CoulombPotential )
  end subroutine FetchVelocityMatrix


  subroutine FetchGeneralVelocityMatrix( MESet, Velocity, LBra, LKet )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Velocity
    integer,                         intent(in)  :: LBra
    integer,                         intent(in)  :: LKet
    !
    type(ClassMatrix) :: R_Derivative
    type(ClassMatrix) :: CoulombPotential
    real(kind(1d0)) :: D, C
    !
    R_Derivative     = MESet.D.Matrix
    !
    if ( abs(LKet-LBra) /= 1 .or. LKet < 0 .or. LBra < 0 ) then
       Velocity = R_Derivative
       Velocity = 0.d0
       return
    end if
    !
    CoulombPotential = MESet.C.Matrix
    !
    !
    if ( LBra > LKet ) then
       D = sqrt(dble(LBra))
       C = D*(-dble(LBra))
    else
       D = -sqrt(dble(LKet))
       C = D*(dble(LKet))
    end if
    !
    call R_Derivative.Multiply( D )
    call CoulombPotential.Multiply( C )
    Velocity = R_Derivative
    call Velocity.Add( CoulombPotential )
    !
  end subroutine FetchGeneralVelocityMatrix


  subroutine FetchGeneralLengthMatrix( MESet, Length, LBra, LKet )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Length
    integer,                         intent(in)  :: LBra
    integer,                         intent(in)  :: LKet
    !
    real(kind(1d0)) :: C
    !
    Length = MESet.R.Matrix
    !
    if ( abs(LKet-LBra) /= 1 .or. LKet < 0 .or. LBra < 0 ) then
       Length = 0.d0
       return
    end if
    !
    if ( LBra > LKet ) then
       C = sqrt(dble(LBra))
    else
       C = -sqrt(dble(LKet))
    end if
    !
    call Length.Multiply( C )
    !
  end subroutine FetchGeneralLengthMatrix


  real(kind(1d0)) function delta( x, y ) result(d)
    !
    integer, intent(in) :: x, y
    !
    if ( x == y ) then
       d = 1.d0
    else
       d = 0.d0
    end if
    !
  end function delta
  


  !> Fetches the Acceleration matrix 
  !! \f[
  !!     \mathbf{L}_{ij}
  !! \f] 
  !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchAccelerationMatrix( MESet, Acceleration, Param1 )
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix)              , intent(out) :: Acceleration
    DoublePrecision,                 intent(in)  :: Param1
    !
    Acceleration = MESet.L.Matrix
    call Acceleration.Multiply( Param1 )
  end subroutine FetchAccelerationMatrix



  !> Fetches the Hamiltonian matrix 
  !! \f[ 
  !!     \mathbf{H}_{ij} = \mathbf{K}_{ij} - 
  !!                       \mathbf{C}_{ij} Z + 
  !!                       \mathbf{L}_{ij} \frac{\ell(\ell+1)}{2}
  !! \f] 
  !! in form of [ClassMatrix](@ref ModuleMatrix.f90).
  subroutine FetchHamiltonianMatrix( MESet, Hamiltonian, L, Z, UseEffectivePot )
    !
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix),               intent(out) :: Hamiltonian
    integer,                         intent(in)  :: L, Z
    logical, optional,               intent(in)  :: UseEffectivePot
    !
    type(ClassMatrix) :: RadialKineticEnergy
    type(ClassMatrix) :: AngularKineticEnergy
    type(ClassMatrix) :: CoulombPotential
    type(ClassMatrix) :: EffectivePotential
    !
    DoublePrecision   :: Par_L, Par_Z
    !
    RadialKineticEnergy = MESet.K.Matrix
    !
    AngularKineticEnergy = MESet.L.Matrix
    Par_L = 0.5d0*dble(L)*(dble(L)+1.d0)
    call AngularKineticEnergy.Multiply( Par_L )
    !
    if ( Present(UseEffectivePot) .and. UseEffectivePot ) then
       CoulombPotential = MESet.Veff.Matrix ! rename to coulomb
    else
       CoulombPotential = MESet.C.Matrix
       Par_Z = -dble(Z)
       call CoulombPotential.Multiply( Par_Z )
    end if
    !
    Hamiltonian = RadialKineticEnergy
    call Hamiltonian.Add( AngularKineticEnergy )
    call Hamiltonian.Add( CoulombPotential )
    !
  end subroutine FetchHamiltonianMatrix


  subroutine FetchKineticMatrix( MESet, Kinetic, L )
    !
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix),               intent(out) :: Kinetic
    integer,                         intent(in)  :: L
    !
    type(ClassMatrix) :: RadialKineticEnergy
    type(ClassMatrix) :: AngularKineticEnergy
    !
    DoublePrecision   :: Par_L
    !
    RadialKineticEnergy = MESet.K.Matrix
    !
    AngularKineticEnergy = MESet.L.Matrix
    Par_L = 0.5d0*dble(L)*(dble(L)+1.d0)
    call AngularKineticEnergy.Multiply( Par_L )
    !
    Kinetic = RadialKineticEnergy
    call Kinetic.Add( AngularKineticEnergy )
    !
  end subroutine FetchKineticMatrix



  subroutine FetchSpecialHamiltonianMatrix( MESet, Hamiltonian, L, Z, Basis )
    !
    class(ClassBasicMatrixElements), intent(in)  :: MESet
    type(ClassMatrix),               intent(out) :: Hamiltonian
    integer,                         intent(in)  :: L, Z
    class(ClassBasis),               intent(in)  :: Basis
    !
    type(ClassMatrix) :: RadialKineticEnergy
    type(ClassMatrix) :: AngularKineticEnergy
    type(ClassMatrix) :: CoulombPotential
    type(ClassMatrix) :: RepulsivePotential
    
    !
    DoublePrecision   :: Par_L, Par_Z
    !
    call BasicMatElemComputeRToMinusNMat( Basis, 1, RepulsivePotential, 0, 0, Basis.GetNodePosition(2) )
    !
    RadialKineticEnergy = MESet.K.Matrix
    !
    AngularKineticEnergy = MESet.L.Matrix
    Par_L = 0.5d0*dble(L)*(dble(L)+1.d0)
    call AngularKineticEnergy.Multiply( Par_L )
    !
    CoulombPotential = MESet.C.Matrix
    Par_Z = -dble(Z+1)
    call CoulombPotential.Multiply( Par_Z )
    !
    Hamiltonian = RadialKineticEnergy
    call Hamiltonian.Add( AngularKineticEnergy )
    call Hamiltonian.Add( CoulombPotential )
    call Hamiltonian.Add( RepulsivePotential )
    !
  end subroutine FetchSpecialHamiltonianMatrix


  !> Fetches the Hamiltonian matrix 
  !! \f[ 
  !!     \mathbf{H}_{ij} = \mathbf{K}_{ij} - 
  !!                       \mathbf{C}_{ij} Z + 
  !!                       \mathbf{L}_{ij} \frac{\ell(\ell+1)}{2}
  !! \f] 
  !! in form of [ClassMatrix](@ref ModuleMatrix.f90), on the
  !! Preconditioned Continuous basis for a given angular momentum
  subroutine FetchPreconditionedContinuousHamiltonian( self, Hamiltonian, L, Z, UseEffectivePot )
    !
    class(ClassBasicMatrixElements), intent(in)  :: self
    type(ClassMatrix),               intent(out) :: Hamiltonian
    integer,                         intent(in)  :: L, Z
    logical, optional,               intent(in)  :: UseEffectivePot
    !
    if ( present(UseEffectivePot) .and. UseEffectivePot ) then
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z, UseEffectivePot )
    else
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z )
    end if
    call self.S.Basis.CastMatrixBra( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         REP_LEVEL_PRECONDITIONED_CONTINUOUS, &
         L )
    call self.S.Basis.CastMatrixKet( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         REP_LEVEL_PRECONDITIONED_CONTINUOUS, &
         L )
    !
  end subroutine FetchPreconditionedContinuousHamiltonian


  !> Fetches the Hamiltonian matrix in real format and in some specified representation level.
  subroutine FetchHamiltonianDouble( self, Hamiltonian, L, Z, RepLevel, UseEffectivePot )
    !
    class(ClassBasicMatrixElements), intent(in)  :: self
    type(ClassMatrix),               intent(out) :: Hamiltonian
    integer,                         intent(in)  :: L, Z
    integer,                         intent(in)  :: RepLevel
    logical, optional,               intent(in)  :: UseEffectivePot
    !
    if ( present(UseEffectivePot) .and. UseEffectivePot ) then
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z, UseEffectivePot )
    else
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z )
    end if
       !
    call self.S.Basis.CastMatrixBra( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         RepLevel, &
         L )
    call self.S.Basis.CastMatrixKet( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         RepLevel, &
         L )
    !
  end subroutine FetchHamiltonianDouble


  !> Fetches the Hamiltonian matrix in complex format and in some specified representation level.
  subroutine FetchHamiltonianComplex( self, ZHamiltonian, L, Z, RepLevel, UseEffectivePot )
    !
    class(ClassBasicMatrixElements), intent(in)  :: self
    type(ClassComplexMatrix),        intent(out) :: ZHamiltonian
    integer,                         intent(in)  :: L, Z
    integer,                         intent(in)  :: RepLevel
    logical, optional,               intent(in)  :: UseEffectivePot
    !
    type(ClassMatrix) :: Hamiltonian
    !
    if ( present(UseEffectivePot) .and. UseEffectivePot ) then
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z, UseEffectivePot )
    else
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z )
    end if
    !
    call self.S.Basis.CastMatrixBra( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         RepLevel, &
         L )
    call self.S.Basis.CastMatrixKet( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         RepLevel, &
         L )
    !
    ZHamiltonian=Hamiltonian
    !
  end subroutine FetchHamiltonianComplex


  !> Fetches the Hamiltonian matrix 
  !! \f[ 
  !!     \mathbf{H}_{ij} = \mathbf{K}_{ij} - 
  !!                       \mathbf{C}_{ij} Z + 
  !!                       \mathbf{L}_{ij} \frac{\ell(\ell+1)}{2}
  !! \f] 
  !! in form of [ClassMatrix](@ref ModuleMatrix.f90), on the
  !! Preconditioned basis for a given angular momentum
  subroutine FetchPreconditionedHamiltonian( self, Hamiltonian, L, Z, UseEffectivePot )
    !
    class(ClassBasicMatrixElements), intent(in)  :: self
    type(ClassMatrix),               intent(out) :: Hamiltonian
    integer,                         intent(in)  :: L, Z
    logical, optional,               intent(in)  :: UseEffectivePot
    !
    if ( present(UseEffectivePot) .and. UseEffectivePot ) then
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z, UseEffectivePot )
    else
       call self.FetchHamiltonianMatrix( Hamiltonian, L, Z )
    end if
    call self.S.Basis.CastMatrixBra( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         REP_LEVEL_PRECONDITIONED, &
         L )
    call self.S.Basis.CastMatrixKet( &
         Hamiltonian, &
         REP_LEVEL_TOTAL, &
         REP_LEVEL_PRECONDITIONED, &
         L )
    !
  end subroutine FetchPreconditionedHamiltonian



end Module ModuleBasicMatrixElements
