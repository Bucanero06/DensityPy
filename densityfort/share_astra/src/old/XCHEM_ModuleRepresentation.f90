
! {{{ Detailed description

!> \file
!!
!! This module contains all the attributes and procedures neccesary to deal with the most appropiated matrix representation for a general monoelectronic local operator on a given basis.

! }}}

Module ModuleRepresentation

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleMatrix
  use ModuleBasis

  implicit none

  private
  
  
  
  !> Defines the matrix representation
  !> of an unspecified operator on a 
  !> given basis
  type, public :: ClassRepresentation
     !
     !private
     !
     !.. LBRA (LKET) are referenced only if the basis level is not TOTAL
     !> Operator matrix bras components representation. Initializes as REP_LEVEL_TOTAL, could be REP_LEVEL_REGULAR, REP_LEVEL_PRECONDITIONED or REP_LEVEL_PRECONDITIONED_CONTINUOUS.
     integer :: BRA_LEVEL = REP_LEVEL_TOTAL
     !> Operator matrix kets components representation. Initializes as REP_LEVEL_TOTAL, could be REP_LEVEL_REGULAR, REP_LEVEL_PRECONDITIONED or REP_LEVEL_PRECONDITIONED_CONTINUOUS.
     integer :: KET_LEVEL = REP_LEVEL_TOTAL
     !> Angular momentum associated to the bras operator matrix components.
     integer :: LBRA = 0
     !> Angular momentum associated to the kets operator matrix components.
     integer :: LKET = 0
     !
     !> Points to the [radial basis class](@ref Basis).
     type(ClassBasis), pointer :: Basis => NULL() !<-- Aggregation
     !> Points to the [matrix class](@ref ModuleMatrix.f90).
     type(ClassMatrix)         :: Matrix          !<-- Composition
     !
   contains

     !> Initializes ClassRepresentation.
     procedure :: Init     =>  RepresentationInit
     !> Deallocates ClassRepresentation.
     procedure :: Free     =>  RepresentationFree
     !> Retrieves the requested operator matrix element.
     procedure :: Element  =>  RepresentationElement
     !> Writes the representation relevant information (radial basis fingerprint and operator matrix) to a file either through its name or unit number.
     generic   :: Write    =>  WriteOnFile, WriteOnUnit
     !> Reads the operator matrix from a file if the current radial basis fingerprint matches with the stored in the same file of the matrix.
     procedure :: Read     =>  RepresentationRead
     !> Computes the \f$f(r)\f$ local operator matrix elements in the general form:
     !! \f[
     !!    \int_{a}^{BP} r^{2}dr
     !!       \frac{d^{n_{1}}\psi_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}\psi_{j}(r)}{dr^{n_{2}}} + 
     !!    \int_{BP}^{b} r^{2}dr
     !!       \frac{d^{n_{1}}\psi_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}\psi_{j}(r)}{dr^{n_{2}}}
     !! \f]
     !! Where \f$BP\f$ is an optional break point and \f$\left\{\psi_{i}(r)\right\}\f$ the radial basis functions. 
     procedure :: Integral =>  RepresentationIntegral
     !> Multiplies the operator representation matrix by a real number.
     procedure :: Multiply =>  Double_x_Matrix
     !> Returns .TRUE. if the requested file is readable, and its content 
     !! is compatible with the current [Basis](@ref Basis) associated to the operator representation;
     !! return .FALSE. otherwise.
     procedure :: Saved
     !> Converts the operator matrix stored in TOTAL representation to one of the others: REGULAR, PRECONDITIONED, PRECONDITIONED_CONTINUOUS.
     procedure :: GetMatrix => RepresentationGetMatrixDriver
     !> Gets the operator preconditioned matrix.
     procedure :: GetPreconditionedMatrix => RepresentationGetPreconditionedMatrix
     !> Gets the operator preconditioned continuous matrix.
     procedure :: GetPreconditionedContinuousMatrix => RepresentationGetPreconditionedContinuousMatrix

     !.. Conversion methods
     !> Retrieves whether the bras representation is TOTAL or not.
     procedure :: BraIsTotal
     !> Retrieves whether the kets representation is TOTAL or not.
     procedure :: KetIsTotal

     !> Casts the operator matrix to TOTAL representation.
     procedure :: CastToTotal
     !> Casts the operator matrix to REGULAR representation.
     procedure :: CastToRegular
     !> Casts the operator matrix to PRECONDITIONED representation.
     procedure :: CastToPreConditioned
     !> Casts the operator matrix to PRECONDITIONED_CONTINUOUS representation.
     procedure :: CastToPreConditionedContinuous

     !> Casts the operator matrix bras components to TOTAL representation.
     procedure :: CastBraToTotal
     !> Casts the operator matrix bras components to REGULAR representation.
     procedure :: CastBraToRegular
     !> Casts the operator matrix bras components to PRECONDITIONED representation.
     procedure :: CastBraToPreConditioned
     !> Casts the operator matrix bras components to PRECONDITIONED_CONTINUOUS representation.
     procedure :: CastBraToPreConditionedContinuous

     !> Casts the operator matrix kets components to TOTAL representation.
     procedure :: CastKetToTotal
     !> Casts the operator matrix kets components to REGULAR representation.
     procedure :: CastKetToRegular
     !> Casts the operator matrix kets components to PRECONDITIONED representation.
     procedure :: CastKetToPreConditioned
     !> Casts the operator matrix kets components to PRECONDITIONED_CONTINUOUS representation.
     procedure :: CastKetToPreConditionedContinuous

!!$     procedure :: MatrixMult
     !
!!$     procedure :: FirstFunctionInvolved
     procedure :: AllocateClassMatrix
     !
     procedure, private :: WriteOnUnit
     procedure, private :: WriteOnFile
     !
  end type ClassRepresentation



contains


  !> Return .TRUE. if FileName is readable, and its content 
  !! is compatible with the current Basis associated to Rep;
  !! return .FALSE. otherwise.
  logical function Saved( Rep, FileName )
    !
    class(ClassRepresentation), intent(in) :: Rep
    character(len=*)          , intent(in) :: FileName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Readable
    character(len=16) :: Form
    !
    Saved=.FALSE.
    !
    !*** Must eliminate enquire, because it does not 
    !    work the way it is expected to if the file is
    !    not already open!
    INQUIRE(&
         File  = FileName, &
         Read  = Readable, &
         Form  = Form    , &
         iostat= iostat    &
         )
    if( iostat /= 0 )return
    if( trim(Readable) /= "YES" )return
    !
    open(&
         Newunit =  uid       , &
         File    =  FileName  , &
         Status  = "old"      , &
         Action  = "read"     , &
         Form    =  trim(Form), &
         iostat  =  iostat    , &
         iomsg   =  iomsg       &
         )
    if(iostat/=0)call Assert(iomsg)
    !
    Saved = Rep.Basis.MatchFingerprintOnUnit( uid )
    !
    close(uid)
    !
  end function Saved


  subroutine Double_x_Matrix( Rep, Factor )
    class(ClassRepresentation), intent(inout) :: Rep
    DoublePrecision           , intent(in)    :: Factor
    call Rep.Matrix.Multiply( Factor )
  end subroutine Double_x_Matrix


  subroutine RepresentationIntegral( &
       Representation    , &
       FunPtr            , & 
       BraDerivativeOrder, &
       KetDerivativeOrder, &
       LowerBound        , &
       UpperBound        , &
       BreakPoint        )
    !
    Class(ClassRepresentation), intent(inout) :: Representation
    procedure(D2DFun)         , pointer       :: FunPtr
    integer        , optional , intent(in)    :: BraDerivativeOrder
    integer        , optional , intent(in)    :: KetDerivativeOrder
    DoublePrecision, optional , intent(in)    :: LowerBound
    DoublePrecision, optional , intent(in)    :: UpperBound
    DoublePrecision, optional , intent(in)    :: BreakPoint
    !
    integer         :: n1, n2, iRow, iCol, NFun, LBW, UBW
    integer         :: iRowMin, iRowMax, iColMin, iColMax
    DoublePrecision, allocatable :: IntegralVec(:)
    DoublePrecision :: a, b, c
    !
    n1 = 0; if(present(BraDerivativeOrder)) n1 = BraDerivativeOrder
    n2 = 0; if(present(KetDerivativeOrder)) n2 = KetDerivativeOrder
    a = -huge(1.d0); if(present(LowerBound)) a = LowerBound
    b =  huge(1.d0); if(present(UpperBound)) b = UpperBound
    c = -huge(1.d0); if(present(BreakPoint)) c = BreakPoint
    !
    NFun = Representation.Basis.GetNFun()
    LBW  = Representation.Basis.LowerBandWidth()
    UBW  = Representation.Basis.UpperBandWidth()
    !
    Representation.Matrix = 0.d0
    !
    !  **** TEST CAREFULLY *****
    !
    !
!!$    !  **** OLD SLOW VERSION *****
!!$    !
!!$    do iCol = 1, NFun
!!$       do iRow = max( 1, iCol - UBW ), min( NFun, iCol + LBW )
!!$          !
!!$          Integral = Representation.Basis.Integral( FunPtr, iRow, iCol, n1, n2, a, b, c )
!!$          call Representation.Matrix.SetElement( iRow, iCol, Integral )
!!$          !
!!$       enddo
!!$    enddo
    !
    ! **** NEW VERSION ****
    !
    allocate(IntegralVec(NFun))
    IntegralVec=0.d0
    do iCol = 1, NFun
       iRowMin = iCol
       iRowMax = min( NFun, iCol + LBW )
       call Representation.Basis.IntegralCol( FunPtr, iRowMin, iRowMax, iCol, IntegralVec, n1, n2, a, b, c )
       do iRow = iRowMin, iRowMax
          call Representation.Matrix.SetElement( iRow, iCol, IntegralVec(iRow-iRowMin+1) )
       enddo
    enddo
    do iRow = 1, NFun-1
       iColMin = iRow+1
       iColMax = min( NFun, iRow + UBW )
       call Representation.Basis.IntegralRow( FunPtr, iRow, iColMin, iColMax, IntegralVec, n1, n2, a, b, c )
       do iCol = iColMin, iColMax
          call Representation.Matrix.SetElement( iRow, iCol, IntegralVec(iCol-iColMin+1) )
       enddo
    enddo
    deallocate(IntegralVec)

    !  **** TEST CAREFULLY *****
    !
    !
  end subroutine RepresentationIntegral


  subroutine RepresentationFree( Representation )
    Class(ClassRepresentation), intent(inout) :: Representation
    if(associated(Representation.Basis))then
       call Representation.Basis.Free()
       Representation.Basis => NULL()
    endif
    call Representation.Matrix.Free()
  end subroutine RepresentationFree


  subroutine RepresentationInit( Representation, Basis )
    Class( ClassRepresentation ), intent(inout) :: Representation
    Class( ClassBasis ), target , intent(in)    :: Basis
    integer :: NFun,Lbw,Ubw
    call Representation.Free()
    if( .not. Basis.Initialized() ) call Assert("Basis not initialized")
    Representation.Basis => Basis
    NFun = Basis.GetNFun()
    if( NFun <= 0 ) call Assert("Invalid Representation Size")
    select case(Basis.BestPattern())
    case( BANDED_PATTERN_IDENTIFIER )
       Lbw  = Basis.LowerBandwidth()
       Ubw  = Basis.UpperBandwidth()
       call Representation.Matrix.InitBanded(NFun,NFun,Lbw,Ubw)
    case DEFAULT
       call Representation.Matrix.InitFull(NFun,NFun)
    end select
  end subroutine RepresentationInit


  subroutine WriteOnFile( Rep, FileName, RequiredForm )
    Class( ClassRepresentation ), intent(in) :: Rep
    character(len=*)            , intent(in) :: FileName
    character(len=*), optional  , intent(in) :: RequiredForm
    !
    character(len=*), parameter :: DEFAULT_FORM = "unformatted"
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=20) :: Form
    integer           :: iostat, uid
    !
    Form = DEFAULT_FORM
    if(present(RequiredForm)) Form = RequiredForm
    !
    open(NewUnit = uid     , &
         File    = FileName, &
         Form    = Form    , &
         Action  = "write" , &
         iostat  = iostat  , &
         iomsg   = iomsg   )
    if(iostat/=0)then
       call ErrorMessage("File "//trim(FileName)//" can't be created")
       call ErrorMessage(iomsg)
       STOP
    endif
    !
    call Rep.Write( uid )
    !
    close( uid )
    !
  end subroutine WriteOnFile


  subroutine WriteOnUnit( Rep, OutputUnit )
    Class( ClassRepresentation ), intent(in) :: Rep
    integer                     , intent(in) :: OutputUnit
    call Rep.Basis.WriteFingerprint( OutputUnit )
    call Rep.Matrix.Write( OutputUnit )
  end subroutine WriteOnUnit


  subroutine RepresentationRead( Representation, InputUnit )
    Class( ClassRepresentation ), intent(inout):: Representation
    integer                     , intent(in)   :: InputUnit
    if(.not.Representation.Basis.MatchFingerprintOnUnit( InputUnit ))call Assert("Invalid Fingerprint")
    call Representation.Matrix.Read( InputUnit )
  end subroutine RepresentationRead


  DoublePrecision function RepresentationElement( Representation, iRow, iCol ) result( Element )
    Class( ClassRepresentation ), intent(in) :: Representation
    integer                     , intent(in) :: iRow, iCol
    Element=Representation.Matrix.Element(iRow,iCol)
  end function RepresentationElement


  subroutine CastToTotal(self)
    class(ClassRepresentation), intent(inout) :: self
    call self.CastBraToTotal()
    call self.CastKetToTotal()
  end subroutine CastToTotal

  subroutine CastToRegular(self,LBra,LKet)
    class(ClassRepresentation), intent(inout) :: self
    integer, intent(in) :: LBra, LKet
    call self.CastBraToRegular(LBra)
    call self.CastKetToRegular(LKet)
  end subroutine CastToRegular

  subroutine CastToPreConditioned(self,LBra,LKet)
    class(ClassRepresentation), intent(inout) :: self
    integer, intent(in) :: LBra, LKet
    call self.CastBraToPreConditioned(LBra)
    call self.CastKetToPreConditioned(LKet)
  end subroutine CastToPreConditioned

  subroutine CastToPreConditionedContinuous(self,LBra,LKet)
    class(ClassRepresentation), intent(inout) :: self
    integer, intent(in) :: LBra, LKet
    call self.CastBraToPreConditionedContinuous(LBra)
    call self.CastKetToPreConditionedContinuous(LKet)
  end subroutine CastToPreConditionedContinuous


  logical function BraIsTotal(self)
    class(ClassRepresentation), intent(in) :: self
    BraIsTotal = self.BRA_LEVEL == REP_LEVEL_TOTAL
  end function BraIsTotal
  logical function BraIsRegular(self)
    class(ClassRepresentation), intent(in) :: self
    BraIsRegular = self.BRA_LEVEL == REP_LEVEL_REGULAR
  end function BraIsRegular
  logical function BraIsPreConditioned(self)
    class(ClassRepresentation), intent(in) :: self
    BraIsPreconditioned = self.BRA_LEVEL == REP_LEVEL_PRECONDITIONED
  end function BraIsPreConditioned
  logical function BraIsPreConditionedContinuous(self)
    class(ClassRepresentation), intent(in) :: self
    BraIsPreconditionedContinuous = self.BRA_LEVEL == REP_LEVEL_PRECONDITIONED_CONTINUOUS
  end function BraIsPreConditionedContinuous
  logical function KetIsTotal(self)
    class(ClassRepresentation), intent(in) :: self
    KetIsTotal = self.BRA_LEVEL == REP_LEVEL_TOTAL
  end function KetIsTotal
  logical function KetIsRegular(self)
    class(ClassRepresentation), intent(in) :: self
    KetIsRegular = self.BRA_LEVEL == REP_LEVEL_REGULAR
  end function KetIsRegular
  logical function KetIsPreconditioned(self)
    class(ClassRepresentation), intent(in) :: self
    KetIsPreconditioned = self.BRA_LEVEL == REP_LEVEL_PRECONDITIONED
  end function KetIsPreconditioned
  logical function KetIsPreConditionedContinuous(self)
    class(ClassRepresentation), intent(in) :: self
    KetIsPreconditionedContinuous = self.BRA_LEVEL == REP_LEVEL_PRECONDITIONED_CONTINUOUS
  end function KetIsPreConditionedContinuous
  
  subroutine CastBraToTotal(self)
    class(ClassRepresentation), intent(inout) :: self
    if(self.BraIsTotal())return
    call self.Basis.CastMatrixBraToTotal( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_TOTAL
  end subroutine CastBraToTotal

  subroutine CastBraToRegular(self,L)
    class(ClassRepresentation), intent(inout) :: self
    integer, optional         , intent(in)    :: L
    !
    if(self.BraIsTotal().and..not.present(L))&
         call Assert("L required in cast of Total Rep to Regular")
    if(present(L).and.self.LBRA/=L)&
         call Assert("L mismatch in Rep cast")
    call self.Basis.CastMatrixBraToRegular( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_REGULAR
  end subroutine CastBraToRegular

  subroutine CastBraToPreconditioned(self,L)
    class(ClassRepresentation), intent(inout) :: self
    integer, optional         , intent(in)    :: L
    !
    if(self.BraIsTotal().and..not.present(L))&
         call Assert("L required in cast of Total Rep to Preconditioned")
    if(present(L).and.self.LBRA/=L)&
         call Assert("L mismatch in Rep cast")
    call self.Basis.CastMatrixBraToPreconditioned( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_PRECONDITIONED
  end subroutine CastBraToPreconditioned

  subroutine CastBraToPreconditionedContinuous(self,L)
    class(ClassRepresentation), intent(inout) :: self
    integer, optional         , intent(in)    :: L
    !
    if(self.BraIsTotal().and..not.present(L))&
         call Assert("L required in cast of Total Rep to PreconditionedContinuous")
    if(present(L).and.self.LBRA/=L)&
         call Assert("L mismatch in Rep cast")
    call self.Basis.CastMatrixBraToPreconditionedContinuous( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_PRECONDITIONED_CONTINUOUS
  end subroutine CastBraToPreconditionedContinuous

  subroutine CastKetToTotal(self)
    class(ClassRepresentation), intent(inout) :: self
    if(self.KetIsTotal())return
    call self.Basis.CastMatrixKetToTotal( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_TOTAL
  end subroutine CastKetToTotal

  subroutine CastKetToRegular(self,L)
    class(ClassRepresentation), intent(inout) :: self
    integer, optional         , intent(in)    :: L
    !
    if(self.KetIsTotal().and..not.present(L))&
         call Assert("L required in cast of Total Rep to Regular")
    if(present(L).and.self.LBRA/=L)&
         call Assert("L mismatch in Rep cast")
    call self.Basis.CastMatrixKetToRegular( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_REGULAR
  end subroutine CastKetToRegular

  subroutine CastKetToPreconditioned(self,L)
    class(ClassRepresentation), intent(inout) :: self
    integer, optional         , intent(in)    :: L
    !
    if(self.KetIsTotal().and..not.present(L))&
         call Assert("L required in cast of Total Rep to Preconditioned")
    if(present(L).and.self.LBRA/=L)&
         call Assert("L mismatch in Rep cast")
    call self.Basis.CastMatrixKetToPreconditioned( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_PRECONDITIONED
  end subroutine CastKetToPreconditioned

  subroutine CastKetToPreconditionedContinuous(self,L)
    class(ClassRepresentation), intent(inout) :: self
    integer, optional         , intent(in)    :: L
    !
    if(self.KetIsTotal().and..not.present(L))&
         call Assert("L required in cast of Total Rep to PreconditionedContinuous")
    if(present(L).and.self.LBRA/=L)&
         call Assert("L mismatch in Rep cast")
    call self.Basis.CastMatrixKetToPreconditionedContinuous( self.Matrix, self.BRA_LEVEL, self.LBRA )
    self.BRA_LEVEL = REP_LEVEL_PRECONDITIONED_CONTINUOUS
  end subroutine CastKetToPreconditionedContinuous


  subroutine RepresentationGetRegularMatrix( Representation, Matrix, LBra, LKet )
    implicit none
    class( ClassRepresentation ), intent(in) :: Representation
    type(ClassMatrix)           , intent(out):: Matrix
    integer                     , intent(in) :: LBra, LKet
    !
    integer, allocatable :: RegularIndexesBra(:)
    integer              :: NRegularIndexesBra
    integer, allocatable :: RegularIndexesKet(:)
    integer              :: NRegularIndexesKet
    !
    call Representation.Basis.GetRegularIndexes( RegularIndexesBra, NRegularIndexesBra, LBra )
    call Representation.Basis.GetRegularIndexes( RegularIndexesKet, NRegularIndexesKet, LKet )
    if(NRegularIndexesBra<=0.or.NRegularIndexesKet<=0)then
       call ErrorMessage("ModuleRepresentation::GetRegularMatrix : Invalid angular momentum")
       call StopExecution()
    endif
    !*** I had to call GetAsymmetricallyIndexedSubmatrix explicitly because it seems that
    !*** the generic GetSubmatrix is unable to find a specific matching subroutine based
    !*** on the arguments (even if it should). Similar problems have been reported for 
    !*** ComposerXE2013.
    call Representation.Matrix.GetAsymmetricallyIndexedSubmatrix( &
         Matrix, RegularIndexesBra, NRegularIndexesBra, RegularIndexesKet, NRegularIndexesKet )
  end subroutine RepresentationGetRegularMatrix
  
  
  subroutine RepresentationGetPreconditionedMatrix( self, Matrix, Lbra, Lket )
    class( ClassRepresentation ), intent(in) :: Self
    type(ClassMatrix)           , intent(out):: Matrix
    integer                     , intent(in) :: LBra
    integer                     , intent(in) :: LKet
    call self.GetMatrix(&
         Matrix, LBra, LKet, REP_LEVEL_PRECONDITIONED, REP_LEVEL_PRECONDITIONED )
  end subroutine RepresentationGetPreconditionedMatrix


  subroutine RepresentationGetPreconditionedContinuousMatrix( self, Matrix, Lbra, Lket )
    class( ClassRepresentation ), intent(in) :: Self
    type(ClassMatrix)           , intent(out):: Matrix
    integer                     , intent(in) :: LBra
    integer                     , intent(in) :: LKet
    call self.GetMatrix(&
         Matrix, LBra, LKet, &
         REP_LEVEL_PRECONDITIONED_CONTINUOUS, &
         REP_LEVEL_PRECONDITIONED_CONTINUOUS )
  end subroutine RepresentationGetPreconditionedContinuousMatrix


  subroutine RepresentationGetMatrixDriver( &
       Self          , &
       Matrix        , &
       LBra          , &
       LKet          , &
       TargetLevelBra, &
       TargetLevelKet  &
       )
    class( ClassRepresentation ), intent(in) :: Self
    type(ClassMatrix)           , intent(out):: Matrix
    integer                     , intent(in) :: LBra
    integer                     , intent(in) :: LKet
    integer                     , intent(in) :: TargetLevelBra
    integer                     , intent(in) :: TargetLevelKet
    !
    character(len=*), parameter :: HERE="ClassRepresentation::GetMatrixDriver : "
    integer :: CurrentBraLevel
    integer :: CurrentKetLevel
    !
    CurrentBraLevel = Self.BRA_LEVEL
    CurrentKetLevel = Self.KET_LEVEL
    !
    !.. Conversions can only be one directional. Therefore, if we are on an L-branch, it is only
    !   legitimate to ask to move either to REP_LEVEL_TOTAL, where the angular momentum is irrelevant,
    !   or to a point within the same branch. Request to convert a representation from one angular
    !   momentum to another is most likely an error, and is thus treated.
    !..
    if( CurrentBraLevel /= REP_LEVEL_TOTAL .and. TargetLevelBra /= REP_LEVEL_TOTAL )then
       if( Self.LBRA /= LBRA ) call Assert(HERE//" Invalid request of inter-L-branch bra conversion.")
    endif
    if( CurrentKetLevel /= REP_LEVEL_TOTAL .and. TargetLevelKet /= REP_LEVEL_TOTAL )then
       if( Self.LKET /= LKET ) call Assert(HERE//" Invalid request of inter-L-branch ket conversion.")
    endif
    !
    Matrix = Self.Matrix
    call Self.Basis.CastMatrixBra( Matrix, Self.BRA_LEVEL, TargetLevelBra, LBra )

    call Self.Basis.CastMatrixKet( Matrix, Self.KET_LEVEL, TargetLevelKet, LKet )
    !
  end subroutine RepresentationGetMatrixDriver
  

!!$  !> ...
!!$  subroutine RepresentationBoxEigInit( &
!!$       Representation, &
!!$       Basis, &
!!$       L, &
!!$       TotalNumBasisFunctions, &
!!$       MinimumRegularIndex, &
!!$       MaximumRegularIndex, &
!!$       NumberRegularFunctions, &
!!$       Lbw, &
!!$       Ubw, &
!!$       Mat, &
!!$       Vec )
!!$    Class( ClassRepresentation ), intent(inout) :: Representation
!!$    Class( ClassBasis ), target , intent(inout) :: Basis
!!$    integer,                      intent(in)    :: L
!!$    integer,                      intent(out)   :: TotalNumBasisFunctions
!!$    integer,                      intent(out)   :: MinimumRegularIndex
!!$    integer,                      intent(out)   :: MaximumRegularIndex
!!$    integer,                      intent(out)   :: NumberRegularFunctions
!!$    integer,                      intent(out)   :: Lbw, Ubw
!!$    DoublePrecision, allocatable, intent(out)   :: Mat(:,:)
!!$    DoublePrecision, allocatable, intent(out)   :: Vec(:)
!!$    !
!!$    Representation.Basis => Basis
!!$    !
!!$    TotalNumBasisFunctions = Basis.GetTotalNumBasisFunctions( )
!!$    MinimumRegularIndex    = Basis.MinimumRegularIndex( L )
!!$    MaximumRegularIndex    = Basis.MaximumRegularIndex( L )
!!$    NumberRegularFunctions = MaximumRegularIndex - MinimumRegularIndex + 1
!!$    !
!!$    select case(Basis.BestPattern())
!!$    case( BANDED_PATTERN_IDENTIFIER )
!!$       Lbw  = Basis.LowerBandwidth()
!!$       Ubw  = Basis.UpperBandwidth()
!!$       call Representation.Matrix.InitBanded( &
!!$            NumberRegularFunctions,&
!!$            NumberRegularFunctions,&
!!$            Lbw,&
!!$            Ubw )
!!$    case DEFAULT
!!$       call Representation.Matrix.InitFull( NumberRegularFunctions,&
!!$            NumberRegularFunctions )
!!$    end select
!!$    call Representation.Matrix.InitBoxEigMatrix( Mat, Vec )  
!!$    !
!!$  end subroutine RepresentationBoxEigInit


!!$  subroutine GetSubMatrix( Representation, MinimumRegularIndex, MaximumRegularIndex, EigenVec )
!!$    Class( ClassRepresentation ), intent(inout) :: Representation
!!$    integer,                      intent(in)    :: MinimumRegularIndex
!!$    integer,                      intent(in)    :: MaximumRegularIndex
!!$    DoublePrecision,              intent(in)    :: EigenVec(:,:)
!!$    call Representation.Matrix.GetSubMatrix( MinimumRegularIndex, MaximumRegularIndex, EigenVec )
!!$  end subroutine GetSubMatrix


!!$  subroutine FetchMatrix( Representation, Mat )
!!$    Class( ClassRepresentation ), intent(inout) :: Representation
!!$    DoublePrecision, allocatable, intent(out)    :: Mat(:,:)
!!$    call Representation.Matrix.InitBoxEigMatrix( Mat )
!!$  end subroutine FetchMatrix


!!$  subroutine FirstFunctionInvolved( Representation, ipot, r, MaximumRegularIndex, i1, i2 )
!!$    class(ClassRepresentation), intent(inout) :: Representation
!!$    integer,                    intent(in)    :: ipot
!!$    DoublePrecision,            intent(in)    :: r(:)
!!$    integer,                    intent(in)    :: MaximumRegularIndex
!!$    integer,                    intent(out)   :: i1, i2
!!$    !
!!$    call Representation.Basis.FirstFunctionInvolved( ipot, r, MaximumRegularIndex, i1, i2 )
!!$  end subroutine FirstFunctionInvolved


!!$  subroutine MatrixMult( Rep, TransA, TRansB, A, B, C )
!!$    class(ClassRepresentation), intent(inout) :: Rep
!!$    character(len=1),          intent(in) :: TransA, TransB
!!$    DoublePrecision,           intent(in)    :: A(:,:), B(:,:)
!!$    DoublePrecision,           intent(out)   :: C(:,:)
!!$    call Rep.Matrix.MatrixMult( TransA, TRansB, A, B, C )
!!$  end subroutine MatrixMult


  subroutine AllocateClassMatrix( Rep, NewClassMatrix )
    class(ClassRepresentation),     intent(in)  :: Rep
    type(ClassMatrix), allocatable, intent(out) :: NewClassMatrix 
    allocate( NewClassMatrix, source = Rep.Matrix )
  end subroutine AllocateClassMatrix


end Module ModuleRepresentation
