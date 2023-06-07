    
! {{{ Detailed description

!> \file
!!
!! Defines a general matrix kind and basic operations on it. Also defines a spectral resolution class devoted to the eigenvalues and eigenvectors management.

! }}}
module ModuleMatrix
        
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV
       
  use ModuleErrorHandling
  use ModuleString
    
  implicit none
  private

  integer         , parameter :: DEFAULT_OUTPUT_UNIT = OUTPUT_UNIT
  integer         , parameter :: DEFAULT_INPUT_UNIT  = INPUT_UNIT

  integer         , parameter :: FILE_NAME_LENGTH    = 512

  real(kind(1d0)) , parameter :: COMPUTATION_THRESHOLD    = 1.d-12
 
  enum, bind(c) 
     enumerator :: MATRIX_PATTERN_FULL
     enumerator :: MATRIX_PATTERN_BANDED
     enumerator :: MATRIX_PATTERN_MIXED
     enumerator :: MATRIX_PATTERN_UNDEFINED
     enumerator :: MATRIX_PATTERN_DIAGONAL
  end enum

  character(len=*), parameter :: DOUBLE_PRINT_FORMAT  = "(d24.16)"!"(x,e14.6)"! "(d24.16)"
  character(len=*), parameter :: COMPLEX_PRINT_FORMAT = "(2(x,d24.16))"

  !> Defines a general DoublePrecision matrix
  type, public :: ClassMatrix
     !
     private
     !
     !.. Pattern : Filled pattern
     !.. NR      : Number of Rows
     !.. NC      : Number of Columns
     !.. NU      : Number of Upper diagonals in banded matrix.
     !.. NL      : Number of Lower diagonals in banded matrix
     !
     !> Pattern which will be used to stored the general matrix, could be 
     !! - Full storage (same array as the original matrix) or 
     !! - Banded storage (takes advantage of banded matrices saving memory 
     !!   space, that would not be possible using the Full storage 
     !!   representation of the original matrix), or 
     !! - Diagonal storage (takes advantage of dagonal matrices).  
     integer :: Pattern = MATRIX_PATTERN_UNDEFINED
     !> Number of matrix rows.
     integer :: NR = 0
     !> Number of matrix columns.
     integer :: NC = 0
     !> Number of matrix superdiagonals.
     integer :: NU = 0
     !> Number of matrix subdiagonals.
     integer :: NL = 0

     !> Minimum matrix row index
     integer :: NRMin = 0

     !> Maximum matrix row index
     integer :: NRMax = 0

     !> Minimum matrix column index
     integer :: NCMin = 0

     !> Maximum matrix column index
     integer :: NCMax = 0

     !> The matrix in the requested representation( Full, Banded or Diagonal).
     Real(kind(1d0)), public, allocatable :: A(:,:) 

   contains
     
     !> Checks if thw ClassMatrix has been initialized
     procedure :: IsInitialized  => ClassMatrixIsInitialized
     !> Initializes the ClassMatrix with the matrix 
     !! stored in the Full representation.
     procedure :: InitFull       => ClassMatrixInitFull

     !> Initializes the ClassMatrix with the matrix 
     !! stored in the Banded representation.
     procedure :: InitBanded     => ClassMatrixInitBanded

     !> Frees all the attributes of the ClassMatrix.
     procedure :: Free           => ClassMatrixFree

     !> Returns True if the matrix pattern is Full.
     procedure :: IsFull         => ClassMatrixIsFull

     !> Returns True if the matrix pattern is Banded.
     procedure :: IsBanded       => ClassMatrixIsBanded

     !> Returns True if the matrix pattern is Diagonal.
     procedure :: IsDiagonal     => ClassMatrixIsDiagonal

     !> Retrieves the number of  subdiagonals in the matrix. 
     procedure :: LowerBandWidth => ClassMatrixLowerBandwidth

     !> Retrieves the number of superdiagonals in the matrix. 
     procedure :: UpperBandWidth => ClassMatrixUpperBandwidth

     !> Performs the summation of the matrix elements absolute value.
     procedure :: Norm1          => ClassMatrixNorm1

     !> Performs the squared root , of the matrix elements squared value summation.
     procedure :: Norm2          => ClassMatrixNorm2

     !> Retrieves the matrix number of rows.
     procedure :: NRows          => ClassMatrixNRows

     !> Retrieves the matrix number of columns.
     procedure :: NColumns       => ClassMatrixNColumns

     !> Fetches a requested matrix elements.
     procedure :: Element        => ClassMatrixElement

     !> Sets a requested matrix elements.
     procedure :: SetElement     => ClassMatrixSetElement

     !> Sets a complete matrix to the matrix class.
     procedure :: SetMatrix      => ClassMatrixSetMatrix

     !> Fetches a complete matrix representation.
     procedure :: FetchMatrix

     !> Fetches from a complete matrix representation, the
     !! submatrix whose rows or columns corresponds to the
     !! indices specified.
     procedure :: FetchMatrixByIndices  =>  FetchClassMatrixByIndices

     !> Fetches the vector in a ClassMatrix when it consists in a one dimensional array.
     procedure :: FetchVector

     !>
     procedure :: Symmetrize

     procedure :: AntiSymmetrize

     procedure :: FindZeroLines => ClassMatrixFindZeroLines

     !> Extracts from a matrix representation corresponding to an 
     !! original matrix, the representation corresponding to an original submatrix.
     !! *** Warning *** I had to make GetAsymmetricallyIndexedSubmatrix public because
     !! GetSubMatrix refuses to resolve the name based on the arguments, even if in 
     !! principle it should.
     generic :: GetSubMatrix   => & 
          ClassMatrixGetSubmatrix_D, &
          ClassMatrixGetSubmatrix_C, &
          ClassMatrixGetSubmatrix_Rect, &
          ClassMatrixGetIndexedSubmatrix,&
          GetAsymmetricallyIndexedSubmatrix

     !> Retrieves whether two matrix classes share the same size or not.
     generic :: operator(.SameSize.) => ClassMatrixSameSize

     !> Retrieves if two matrix classes share the same size or not.
     procedure, private :: ClassMatrixSameSize

     !> Takes the Banded representation pertaining to an original 
     !! squared matrix (with the same number of lower and superdiagonals) 
     !! and turns it in the Full representation.
     procedure :: ConvertToSquared

     !> Writes the relevant information of the matrix class to a unit.
     generic :: Write => &
          WriteClassMatrixToUnit, &
          WriteClassMatrixToFile

     !> Reads the relevant information of the matrix class from a unit.
     generic :: Read  => &
          ReadClassMatrixFromUnit, &
          ReadClassMatrixFromFile

     !> Performs the matrix multiplication, of the array contained in the 
     !! matrix class, either by a number or by another matrix pertaining 
     !! to other matrix class.
     generic   :: Multiply     => &
          Double_x_ClassMatrix, &
          ClassMatrix_x_ClassMatrix, &
          ClassMatrix_x_Matrix

     !> Adds up two matrices pertaining to two different matrix classes.
     procedure :: Add            => AddMatrixToMatrix

     !> Starting from two matrices (\f$A\f$, \f$B\f$) and a number (\f$\alpha\f$), 
     !! constructs a new one (\f$G\f$) more appropiated to solve the system 
     !! of linear equations.
     !!
     !!\f[
     !! G=A-\alpha\cdot B
     !!\f]
     procedure :: Compose      => ClassMatrixCompose

     procedure :: ComposeFromBlocks      => ClassMatrixComposeFromBlocks
     
     procedure :: BuildUpMatrix      => ClassMatrixBuildUpMatrix

     !> Performs the LU factorization of a matrix.
     procedure :: Factorize      => ClassMatrixFactorize

     !> Solves a system of linear equations once the matrix has been previously LU factorized.
     generic :: LinEqSolver      => ClassMatrixLinEqSolver, ClassMatrixLinEqSolverGeneral

     !> Selects which kind of linear equation solver will be used related with 
     !! the matrix pattern, when the method used is Inverse Iteration.
     procedure :: ClassMatrixSelectLinEqSolver  

     !> Checks if the Matrix has at least one zero row.
     procedure :: CheckZeroRowOrColumn

     !> Computes the inverse of the matrix.
     procedure, public :: Inverse => ClassMatrixInverse

     procedure, public :: IsIdentity => ClassMatrixIsIdentity

     procedure, public :: MatrixIsDiagonal => ClassMatrixMatrixIsDiagonal

     procedure, public :: IsZero => ClassMatrixIsZero

     procedure, public :: IsSymmetric => ClassMatrixIsSymmetric

     procedure, public :: IsAntiSymmetric => ClassMatrixIsAntiSymmetric

     generic, public :: Transpose => ClassMatrixTranspose, ClassMatrixTransposeOld

     !> Allows to a matrix class, to inherit its attributes from another, or to 
     !! set its matrix equal to a number (integer or real).
     generic :: assignment(=) => &
          AssignDoubleToClassMatrix, &
          AssignIntegerToClassMatrix, &
          CopyClassMatrixToClassMatrix, &
          CopyMatrixToClassMatrix, &
          CopyVectorToClassMatrix

     !> Solves the eigenvalues and eigenvectors problem for a matrix.
     generic   :: Diagonalize => GeneralEigenValueSolver, EigenvalueSolver, GeneralEigenValueSolverConditionNumber
     !> Compute the transformation matrix that condition an original matrix containing linear dependencies.
     procedure :: ConditionMatrix => ClassMatrixConditionMatrix
     !> Gets the minimum eigenvalue of a matrix.
     procedure :: GetMinEigVal    => ClassMatrixGetMinEigVal
     !> Gets the maxiimum eigenvalue of a matrix.
     procedure :: GetMaxEigVal    => ClassMatrixGetMaxEigVal
     !> Adds new rows full with zeros to an existent matrix, in the position specified.
     procedure :: AddRows => ClassMatrixAddRows
     !> Adds new columns full with zeros to an existent matrix, in the position specified.
     procedure :: AddColumns => ClassMatrixAddColumns
     !> Removes the requested amount of rows from a matrix, in the position specified.
     procedure :: RemoveRows => ClassMatrixRemoveRows
     !> Removes the requested amount of columns from a matrix, in the position specified.
     procedure :: RemoveColumns => ClassMatrixRemoveColumns

     procedure, private :: ClassMatrixIsInitialized
     procedure, private :: ReadClassMatrixFromUnit
     procedure, private :: ReadClassMatrixFromFile
     procedure, private :: WriteClassMatrixToUnit
     procedure, private :: WriteClassMatrixToFile
     procedure, private :: ClassMatrixGetSubmatrix_D
     procedure, private :: ClassMatrixGetSubmatrix_C
     procedure, private :: ClassMatrixGetSubmatrix_Rect
     procedure, private :: ClassMatrixGetIndexedSubmatrix
     !> Gets a rectangular submatrix from the ClassMatrix's matrix 
     !! defined by two independent vectors of row and column indexes
     procedure, public  :: GetAsymmetricallyIndexedSubmatrix
     procedure, private :: AssignDoubleToClassMatrix
     procedure, private :: AssignIntegerToClassMatrix
     procedure, private :: CopyClassMatrixToClassMatrix
     procedure, private :: CopyMatrixToClassMatrix
     procedure, private :: CopyVectorToClassMatrix
     procedure, private :: GeneralEigenValueSolver
     procedure, private :: GeneralEigenValueSolverConditionNumber
     procedure, private :: EigenvalueSolver
     procedure, private :: ClassMatrixConditionMatrix
     procedure, private :: HasSameShapeAs
     procedure, private :: Double_x_ClassMatrix
     procedure, private :: ClassMatrix_x_ClassMatrix
     procedure, private :: ClassMatrix_x_Matrix
     procedure, private :: ClassMatrixLinEqSolver
     procedure, private :: ClassMatrixLinEqSolverGeneral
     procedure, private :: ClassMatrixTranspose
     procedure, private :: ClassMatrixTransposeOld

!!$     generic   :: operator(+)   =>   ClassMatrixPlusClassMatrix
!!$          ClassMatrixPlusInteger    , &
!!$          ClassMatrixPlusDouble     , &
!!$          ClassMatrixPlusMatrix     , &
!!$          IntegerPlusClassMatrix    , &
!!$          DoublePlusClassMatrix     , &
!!$          MatrixPlusClassMatrix     
!!$     !
!!$     procedure :: ClassMatrixPlusClassMatrix
!!$     generic   :: operator(*)   =>    &
!!$          ClassMatrixTimesDouble, &
!!$          DoubleTimesClassMatrix
!!$          ClassMatrixTimesClassMatrix, &
!!$          ClassMatrixTimesInteger    , &
!!$          ClassMatrixTimesMatrix     , &
!!$          IntegerTimesClassMatrix    , &
!!$          DoubleTimesClassMatrix     , &
!!$          MatrixTimesClassMatrix  

     !> Deallocates all ClassMatrix atributes.
     final :: ClassMatrixFinalize

  end type ClassMatrix


  !> Adequate class for the eigenvalues and eigenvectors management derived from real matrices.
  type, public :: ClassSpectralResolution
     !
     private
     !> Number of eigenvalues.
     integer                      :: NEigenvalues
     !> Number of Eigenvectors matrix rows. 
     !! Dim <= NEigenvalues
     integer                      :: Dim
     !> A vector to store the eigenvalues.
     Real(kind(1d0)), allocatable :: EigenValues(:)
     !> A matrix to store the eigenvectors.
     Real(kind(1d0)), public,allocatable :: EigenVectors(:,:)
     !
   contains
     !
     !> Retrieves the number of eigenvalues in a spectral resolution class.
     procedure :: NEval =>  NevalSpectralResolution
     !> Retrieves the number of Eigenvectors matrix rows (variable Dim) 
     !! in a spectral resolution class.
     procedure :: Size  =>  SizeSpectralResolution
     !> Initializes the spectral resolution class.
     generic   :: Init  =>  InitSpectralResolutionFull, InitSpectralResolutionReduced
     !> Retrieves whether a previously stored spectral resolution class 
     !! information in a unit, is consistent or not with a new one available.
     procedure :: IsConsistent =>  SpectralResolutionIsConsistent
     !> Reads from a file or a unit, the previously stored spectral resolution 
     !! class information.
     generic   :: Read =>  ReadSpectralResolutionFromUnit, ReadSpectralResolutionFromFile
     !> Writes in a file or a unit, the spectral resolution class information.
     generic   :: Write =>  WriteSpectralResolutionToUnit, WriteSpectralResolutionToFile
     !> Writes in a file or a unit, the spectral resolution class eigenvalues.
     generic   :: WriteEigenvalues =>  WriteEigenvaluesToUnit, WriteEigenvaluesToFile
     !> Can fetch from a spectral resolution class either the eigenvalues, 
     !! the eigenvectors, a single eigenvector or generate a matrix class 
     !! with the eigenvector matrix.
     generic   :: Fetch =>  &
          FetchEigenValues, &
          FetchEigenVectors, &
          FetchSingleEigenvector, &
          FetchSingleEigenvectorMat, &
          FetchClassMatrixEigenvectors
     !> Given two spectral resolutions with eigenvector matrices \f$C\f$ and 
     !! \f$U\f$, transforms the former representation via the latter to obtain 
     !! the new one \f$C'\f$:
     !!
     !!\f[
     !! C'=U^{T}\cdot C,
     !!\f]
     !!
     !! if the metric is not taken into account; or
     !!
     !!\f[
     !! C'=U^{T}\cdot S\cdot C,
     !!\f]
     !!
     !! if the metric is present.
     generic   :: Transform        =>  TransformEigenvectors, TransformEigenvectorsWithMetric
     !> Set to + the sign of the first entry in each eigenvector
     procedure :: SyncFirstSign
     !> Eliminates the null eigenspace
     procedure :: PurgeNull        =>  SpectralResolutionPurgeNull
     !> Sort eigenvalues and and its corresponding eigenvectors in ascending order of the former.
     procedure :: Sort             =>  ClassSpectralResolutionSort
     !> Intended to remove the last N basis functions.
     procedure :: ReduceDimension  =>  ClassSpectralResolutionReduceDimension
     !> Set a choosen spectral resolution eigenvalue equal to some external one, or the entire vector of eigenvalues.
     generic :: SetEigenValues     =>  SetOneEigenValue, SetAllEigenValues
     generic :: SetEigenVectors    =>  ClassSpectralResolutionSetEigenVectors, ClassSpectralResolutionSetOneEigenVector, ClassSpectralResolutionSetEigenVectorsMat
     !> Frees the spectral resolution class attributes.
     procedure :: Free             => FreeSpectralResolution
     !
     procedure, private :: FetchEigenValues
     procedure, private :: FetchEigenVectors
     procedure, private :: FetchSingleEigenvector
     procedure, private :: FetchSingleEigenvectorMat
     procedure, private :: FetchClassMatrixEigenvectors
     procedure, private :: InitSpectralResolutionFull
     procedure, private :: InitSpectralResolutionReduced
     procedure, private :: WriteSpectralResolutionToUnit
     procedure, private :: WriteSpectralResolutionToFile
     procedure, private :: ReadSpectralResolutionFromUnit
     procedure, private :: ReadSpectralResolutionFromFile
     procedure, private :: WriteEigenvaluesToUnit
     procedure, private :: WriteEigenvaluesToFile
     procedure, private :: TransformEigenvectors
     procedure, private :: TransformEigenvectorsWithMetric
     procedure, private :: SetOneEigenValue
     procedure, private :: SetAllEigenValues
     procedure, private :: ClassSpectralResolutionSetEigenVectors
     procedure, private :: ClassSpectralResolutionSetEigenVectorsMat
     procedure, private :: ClassSpectralResolutionSetOneEigenVector
     procedure, private :: ClassSpectralResolutionReduceDimension
     !
  end type ClassSpectralResolution


  !> Defines a general Complex(kind(1d0)) matrix class.
  type, public :: ClassComplexMatrix
     !
     private
     !
     !.. Pattern : Filled pattern
     !.. NR      : Number of Rows
     !.. NC      : Number of Columns
     !.. NU      : Number of Upper diagonals in banded matrix.
     !.. NL      : Number of Lower diagonals in banded matrix.
     !
     !> Pattern which will be used to stored the general matrix, could be 
     !! - Full storage (same array as the original matrix) or 
     !! - Banded storage (takes advantage of banded matrices saving memory 
     !!   space, that would not be possible using the Full storage 
     !!   representation of the original matrix), or 
     !! - Diagonal storage (takes advantage of dagonal matrices).  
     integer :: Pattern = MATRIX_PATTERN_UNDEFINED
     !
     !> Number of matrix rows.
     integer :: NR = 0
     !> Number of matrix columns.
     integer :: NC = 0
     !> Number of matrix superdiagonals.
     integer :: NU = 0
     !> Number of matrix subdiagonals.
     integer :: NL = 0
     !
     !> Minimum matrix row index
     integer :: NRMin = 0
     !> Maximum matrix row index
     integer :: NRMax = 0
     !> Minimum matrix column index
     integer :: NCMin = 0
     !> Maximum matrix column index
     integer :: NCMax = 0
     !
     !> The matrix in the requested representation( Full, Banded or Diagonal).
     Complex(kind(1d0)), public, allocatable :: A(:,:)
     !
   contains
     !
     !> Checks if thw ClassMatrix has been initialized
     procedure :: IsInitialized  => ClassComplexMatrixIsInitialized

     !> Initializes the ClassComplexMatrix with 
     !! the matrix stored in the Full representation.
     procedure :: InitFull       => ClassComplexMatrixInitFull

     !> Initializes the ClassComplexMatrix with 
     !! the matrix stored in the Banded representation.
     procedure :: InitBanded     => ClassComplexMatrixInitBanded

     !> Frees all the attributes of the ClassComplexMatrix.
     procedure :: Free           => ClassComplexMatrixFree

     !> Returns True if the matrix pattern is Full.
     procedure :: IsFull         => ClassComplexMatrixIsFull

     !> Returns True if the matrix pattern is Banded.
     procedure :: IsBanded       => ClassComplexMatrixIsBanded

     !> Returns True if the matrix pattern is Diagonal.
     procedure :: IsDiagonal     => ClassComplexMatrixIsDiagonal

     !> Retrieves the number of  subdiagonals in the matrix. 
     procedure :: LowerBandWidth => ClassComplexMatrixLowerBandwidth

     !> Retrieves the number of superdiagonals in the matrix. 
     procedure :: UpperBandWidth => ClassComplexMatrixUpperBandwidth
 
     procedure :: ConvertToFull  => ClassComplexMatrixConvertToFull

     !> Performs the summation of the matrix elements absolute value.
     procedure :: Norm1          => ClassComplexMatrixNorm1

     !> Performs the squared root , of the matrix elements absolute squared value summation.
     procedure :: Norm2          => ClassComplexMatrixNorm2

     !> Retrieves the matrix number of rows.
     procedure :: NRows          => ClassComplexMatrixNRows

     !> Retrieves the matrix number of columns.
     procedure :: NColumns       => ClassComplexMatrixNColumns

     !> Fetches a requested matrix elements.
     procedure :: Element        => ClassComplexMatrixElement

     !> Fetches a complete matrix representation.
     procedure :: FetchMatrix    => FetchComplexMatrix

     !> Fetches from a complete complex matrix representation, the
     !! submatrix whose rows or columns corresponds to the
     !! indices specified.
     procedure :: FetchMatrixByIndices => FetchClassComplexMatrixByIndices

     !> Adds new rows full with complex zeros to an existent complex matrix, in the position specified.
     procedure :: AddRows => ClassComplexMatrixAddRows
     !> Adds new columns full with complex zeros to an existent complex matrix, in the position specified.
     procedure :: AddColumns => ClassComplexMatrixAddColumns

     !> Removes the requested amount of rows from a complex matrix, in the position specified.
     generic :: RemoveRows => ClassComplexMatrixRemoveRowsWhere, ClassComplexMatrixRemoveRowsAfter
     !> Removes the requested amount of columns from a complex matrix, in the position specified.
     generic :: RemoveColumns => ClassComplexMatrixRemoveColumnsWhere, ClassComplexMatrixRemoveColumnsAfter

     !> Removes all the zero rows of a complex matrix
     procedure :: RemoveAllZeroRows => ClassComplexMatrixRemoveAllZeroRows
     !> Removes all the zero columns of a complex matrix
     procedure :: RemoveAllZeroColumns => ClassComplexMatrixRemoveAllZeroColumns

     !> Providing a complex matrix class, it is placed in the correct position inside a higher complex matrix class.
     procedure :: ComposeFromBlocks => ClassComplexMatrixComposeFromBlocks

     procedure :: BuildUpMatrix => ClassComplexMatrixBuildUpMatrix

     !> Sets a requested matrix elements.
     procedure :: SetElement     => ClassComplexMatrixSetElement

     procedure :: Compose        => ClassComplexMatrixCompose

     procedure :: Factorize      => ClassComplexMatrixFactorize

     !> Solves a system of linear equations once the complex matrix has been previously LU factorized.
     generic :: LinEqSolver      => ClassComplexMatrixLinEqSolver, ClassComplexMatrixLinEqSolverGeneral

     procedure :: ClassComplexMatrixSelectLinEqSolver  

     !> Extracts from a matrix representation corresponding to 
     !! an original matrix, the representation corresponding 
     !! to an original submatrix.
     generic :: GetSubMatrix =>  ClassComplexMatrixGetSubmatrix, ClassComplexMatrixGetSubmatrix_Rect

     !> Writes the relevant information of the complex matrix class to a unit.
     procedure :: Write          => ClassComplexMatrixWriteToUnit

     !> Reads the relevant information of the complex matrix class from a unit.
     procedure :: Read           => ClassComplexMatrixReadFromUnit

     !> Performs the matrix multiplication, of the array contained in the complex matrix class by a complex number.
     generic   :: Multiply       => ClassComplexMatrixMultiplyByComplex, ClassComplexMatrixMultiplyByClassMatrix, ClassComplexMatrixMultiplyByClassComplexMatrix

     !> Adds up two matrices pertaining to two different complex matrix classes.
     generic :: Add            => ClassComplexMatrixAddClassComplexMatrix, ClassComplexMatrixAddClassMatrix
     
     !> Returns the real part of a complex matrix through a subroutine
     generic   :: RePart             => ClassComplexMatrixRealPartSub

     !> Returns the real part of a complex matrix through a function
     generic   :: Re             => ClassComplexMatrixRealPartFun

     !> Returns the imaginary part of a complex matrix through a subroutine
     generic   :: ImPart             => ClassComplexMatrixImaginaryPartSub

     !> Returns the imaginary part of a complex matrix through a function
     generic   :: Im             =>  ClassComplexMatrixImaginaryPartFun

     !> Allows to a complex matrix class, either to inherit its attributes 
     !! from another complex matrix class, or simply from a matrix class. 
     !! Also permits to set its matrix equal to a number (integer, real or complex).
     generic   :: assignment(=) => &
          ClassComplexMatrixAssignComplex, &
          ClassComplexMatrixAssignDouble, &
          ClassComplexMatrixAssignInteger, &
          ClassComplexMatrixAssignComplexMatrixFull, &
          ClassComplexMatrixCopyToClassComplexMatrix, &
          ClassComplexMatrixCopyClassMatrix, &
          ClassComplexMatrixAssignRealMatrix

     !> Solves the eigenvalues and eigenvectors problem for a complex matrix.
     generic   :: Diagonalize => &
          GeneralComplexEigenvalueSolver,&
          GeneralComplexEigenvalueSolverRealMetric,&
          ComplexEigenvalueSolver, &
          GeneralComplexEigenvalueSolverWithRegularization

     generic :: Regularize => GeneralComplexMatrixRegularization

     generic, public :: TransposeConjugate => ClassComplexMatrixTransposeConjugate, ClassComplexMatrixTransposeConjugateOld

     procedure :: Symmetrize => ClassComplexMatrixSymmetrize

     procedure :: AntiSymmetrize => ClassComplexMatrixAntiSymmetrize
     
     procedure :: RemoveEpsilon

     procedure, public :: IsHermitian => ClassComplexMatrixIsHermitian
     procedure, public :: IsSymmetric => ClassComplexMatrixIsSymmetric
     procedure, public :: IsAntiSymmetric => ClassComplexMatrixIsAntiSymmetric
     procedure, public :: IsReal => ClassComplexMatrixIsReal

     !> Computes the inverse of the complex matrix.
     procedure, public :: Inverse => ClassComplexMatrixInverse

     procedure, public :: IsIdentity => ClassComplexMatrixIsIdentity

     procedure, public :: MatrixIsDiagonal => ClassComplexMatrixMatrixIsDiagonal

     procedure, public :: IsUnitary => ClassComplexMatrixIsUnitary

     procedure, public :: IsZero => ClassComplexMatrixIsZero

     procedure, public :: IsComplex => ClassComplexMatrixIsComplex

     procedure, public :: Transpose => ClassComplexMatrixTranspose

     procedure, public :: SolveHomogScattExpert => ClassComplexMatrixSolveHomogScattExpert

!!$     procedure :: ClassComplexMatrixAllZeros



!!$     procedure, private :: ClassComplexMatrixAllZeros
     procedure, private :: ClassComplexMatrixRealPartSub
     procedure, private :: ClassComplexMatrixRealPartFun
     procedure, private :: ClassComplexMatrixImaginaryPartSub     
     procedure, private :: ClassComplexMatrixImaginaryPartFun     
     procedure, private :: ClassComplexMatrixGetSubmatrix
     procedure, private :: ClassComplexMatrixGetSubmatrix_Rect

     procedure, private :: ClassComplexMatrixAssignComplex
     procedure, private :: ClassComplexMatrixAssignDouble
     procedure, private :: ClassComplexMatrixAssignInteger
     procedure, private :: ClassComplexMatrixAssignComplexMatrixFull

     procedure, private :: ClassComplexMatrixCopyToClassComplexMatrix
     procedure, private :: ClassComplexMatrixCopyClassMatrix
     procedure, private :: ClassComplexMatrixAssignRealMatrix
     procedure, private :: GeneralComplexEigenvalueSolverRealMetric
     procedure, private :: GeneralComplexEigenvalueSolver
     procedure, private :: ComplexEigenvalueSolver
     procedure, private :: GeneralComplexEigenvalueSolverWithRegularization
     procedure, private :: GeneralComplexMatrixRegularization
     procedure, private :: ClassComplexMatrixAddRows
     procedure, private :: ClassComplexMatrixAddColumns
     procedure, private :: HasSameShapeAs => ClassComplexMatrixHasSameShapeAs
     procedure, private :: ClassComplexMatrixMultiplyByComplex
     procedure, private :: ClassComplexMatrixMultiplyByClassMatrix
     procedure, private :: ClassComplexMatrixMultiplyByClassComplexMatrix
     procedure, private :: ClassComplexMatrixRemoveRowsWhere
     procedure, private :: ClassComplexMatrixRemoveRowsAfter
     procedure, private :: ClassComplexMatrixRemoveColumnsWhere
     procedure, private :: ClassComplexMatrixRemoveColumnsAfter
     procedure, private :: ClassComplexMatrixLinEqSolver
     procedure, private :: ClassComplexMatrixLinEqSolverGeneral
     procedure, private :: ClassComplexMatrixTranspose
     procedure, private :: ClassComplexMatrixTransposeConjugate
     procedure, private :: ClassComplexMatrixTransposeConjugateOld
     procedure, private :: ClassComplexMatrixAddClassComplexMatrix
     procedure, private :: ClassComplexMatrixAddClassMatrix

     !> Deallocates all ClassComplexMatrix atributes.
     final :: ClassComplexMatrixFinalize

  end type ClassComplexMatrix

  !> Adequate class for the eigenvalues and eigenvectors management derived from complex matrices.
  type, public :: ClassComplexSpectralResolution
     !
     private
     !> The number of eigenvalues.
     integer                         :: NEigenvalues
     !> The number of rows that the eigenvector matrices will have. 
     !! Usually will coincide with the number of eigenvalues (number 
     !! of eigenvector matrices columns) but not always.
     integer                         :: Dim
     !> A vector to store the eigenvalues.
     Complex(kind(1d0)), allocatable :: EigenValues(:)
     !
     !*** We should seriously consider whether to change this to
     !    a ClassComplexMatrix 
     !> A matrix to store the left eigenvectors.
     Complex(kind(1d0)), allocatable :: LeftEigenVectors(:,:)
     !> A matrix to store the right eigenvectors.
     Complex(kind(1d0)), allocatable :: RightEigenVectors(:,:)
     !
   contains
     !
     !> Retrieves the number of eigenvalues in a complex spectral resolution class.
     procedure :: NEval            =>  NevalComplexSpectralResolution
     !> Retrieves the number of eigenvectors matrices rows (variable Dim) in a complex spectral resolution class.
     procedure :: Size             =>  SizeComplexSpectralResolution
     !> Initializes the complex spectral resolution class.
     generic   :: Init             =>  InitComplexSpectralResolutionFull, InitComplexSpectralResolutionReduced
     !> Frees the complex spectral resolution class attributes.
     procedure :: Free             =>  ClassComplexSpectralResolutionFree
     !> Retrieves whether a previously stored complex spectral resolution class information in a unit, is consistent or not with a new one available.
     procedure :: IsConsistent     =>  ComplexSpectralResolutionIsConsistent
     !> Reads from a file or a unit, the previously stored complex spectral resolution class information.
     generic   :: Read             =>  ReadComplexSpectralResolutionFromUnit, ReadComplexSpectralResolutionFromFile
     !> Writes in a file or a unit, the complex spectral resolution class information.
     generic   :: Write            =>  WriteComplexSpectralResolutionToUnit, WriteComplexSpectralResolutionToFile
     !> Writes in a file or a unit, the complex spectral resolution class eigenvalues.
     generic   :: WriteEigenvalues =>  WriteComplexEigenvaluesToUnit, WriteComplexEigenvaluesToFile
     !> Can fetch from a complex spectral resolution class either the eigenvalues, the eigenvectors or a single eigenvector.
     generic   :: Fetch            =>  FetchComplexEigenvalues, FetchComplexEigenvectors, FetchComplexEigenvectorsMat, FetchSingleComplexEigenvector, FetchSingleComplexEigenvectorMat
     !> Sets the eigenvalues of the complex spectral resolution.
     generic :: SetEigenValues   =>  SetComplexSpectralResolutionEigenValues
     !> Sets the eigenvectors of the complex spectral resolution.
     generic :: SetEigenVectors   =>  SetComplexSpectralResolutionEigenVectors, SetComplexSpectralResolutionEigenVectorsMat
     !> Given two spectral resolutions (complex or not)  with eigenvector matrices \f$C\f$ and \f$U\f$ (in the case of complex spectral resolutions the metrices are the corresponding to the right eigenvectors), transforms the former representation via the latter to obtain the new one \f$C'\f$:
     !!
     !!\f[
     !! C'=U^{T}\cdot C,
     !!\f]
     !!
     !! if the metric is not taken into account; or
     !!
     !!\f[
     !! C'=U^{T}\cdot S\cdot C,
     !!\f]
     !!
     !! if the metric is present, which could be in ClassMatrix or ClassComplexMatrix.
     generic   :: Transform        =>  &
          TransformComplexEigenvectors_C, &
          TransformComplexEigenvectors_D, &
          TransformComplexEigenvectorsWithMetric_CC, &
          TransformComplexEigenvectorsWithMetric_CD, &
          TransformComplexEigenvectorsWithMetric_DC, &
          TransformComplexEigenvectorsWithMetric_DD, &
          TransformComplexEigenvectors_ClassMatrix, &
          TransformComplexEigenvectors_ClassComplexMatrix
     !> Set to + the sign of the first entry in each eigenvector
     procedure :: SyncFirstSign    => ComplexSpectralResolutionSyncFirstSign
     !> Eliminates the null eigenspace
     procedure :: PurgeNull        => ComplexSpectralResolutionPurgeNull
     !> Sort eigenvalues and and its corresponding eigenvectors in ascending order of the former's real part.
     procedure :: Sort             => ClassComplexSpectralResolutionSort
     !> Intended to remove the last N basis functions.
     procedure :: ReduceDimension  =>  ClassComplexSpectralResolutionReduceDimension
     !
     !> The complex spectral resolution inherits its attributes from a non-complex spectral resolution.
     generic   :: assignment(=) => DoubleSpectralResolutionToComplexSpectralResolution
     !
     procedure, private :: FetchComplexEigenvalues
     procedure, private :: FetchComplexEigenvectors
     procedure, private :: FetchComplexEigenvectorsMat
     procedure, private :: FetchSingleComplexEigenvector
     procedure, private :: FetchSingleComplexEigenvectorMat
     procedure, private :: InitComplexSpectralResolutionFull
     procedure, private :: InitComplexSpectralResolutionReduced
     procedure, private :: WriteComplexSpectralResolutionToUnit
     procedure, private :: WriteComplexSpectralResolutionToFile
     procedure, private :: ReadComplexSpectralResolutionFromUnit
     procedure, private :: ReadComplexSpectralResolutionFromFile
     procedure, private :: WriteComplexEigenvaluesToUnit
     procedure, private :: WriteComplexEigenvaluesToFile
     procedure, private :: TransformComplexEigenvectors_C
     procedure, private :: TransformComplexEigenvectors_D
     procedure, private :: TransformComplexEigenvectorsWithMetric_CC
     procedure, private :: TransformComplexEigenvectorsWithMetric_CD
     procedure, private :: TransformComplexEigenvectorsWithMetric_DC
     procedure, private :: TransformComplexEigenvectorsWithMetric_DD
     procedure, private :: TransformComplexEigenvectors_ClassMatrix
     procedure, private :: TransformComplexEigenvectors_ClassComplexMatrix
     procedure, private :: DoubleSpectralResolutionToComplexSpectralResolution
     procedure, private :: SetComplexSpectralResolutionEigenValues
     procedure, private :: SetComplexSpectralResolutionEigenVectors
     procedure, private :: SetComplexSpectralResolutionEigenVectorsMat
     procedure, private :: ClassComplexSpectralResolutionReduceDimension
     !
     !> Deallocates all ClassComplexSpectralResoluion atributes.
     final :: ClassComplexSpectralResolutionFinal
     !
  end type ClassComplexSpectralResolution


  interface AllocateMatrix
     module procedure AllocateDoubleMatrix, AllocateComplexMatrix
  end interface AllocateMatrix


  public :: ReScaleRightEigVectors
  public :: NormalizedAsHermitian
  


contains


  !> Returns True if the matrix pattern is Full.
  logical function ClassMatrixIsFull(Matrix) result(IsFull)
    Class(ClassMatrix), intent(in) :: Matrix
    IsFull = ( Matrix.Pattern == MATRIX_PATTERN_FULL )
  end function ClassMatrixIsFull
  !
  !> Returns True if the matrix pattern is Banded.
  logical function ClassMatrixIsBanded(Matrix) result(IsBanded)
    Class(ClassMatrix), intent(in) :: Matrix
    IsBanded = ( Matrix.Pattern == MATRIX_PATTERN_BANDED )
  end function ClassMatrixIsBanded
  !
  !> Returns True if the matrix pattern is Diagonal.
  logical function ClassMatrixIsDiagonal(Matrix) result(IsDiagonal)
    Class(ClassMatrix), intent(in) :: Matrix
    IsDiagonal = ( Matrix.Pattern == MATRIX_PATTERN_DIAGONAL )
  end function ClassMatrixIsDiagonal

  !> Retrieves the number of  subdiagonals in the matrix. 
  integer function ClassMatrixLowerBandwidth(Matrix) result(LowerBandwidth)
    Class(ClassMatrix), intent(in) :: Matrix
    LowerBandwidth = Matrix.NL
  end function ClassMatrixLowerBandwidth

  !> Retrieves the number of superdiagonals in the matrix. 
  integer function ClassMatrixUpperBandwidth(Matrix) result(UpperBandwidth)
    Class(ClassMatrix), intent(in) :: Matrix
    UpperBandwidth = Matrix.NU
  end function ClassMatrixUpperBandwidth


  !> Write matrix to file
  subroutine WriteClassMatrixToFile( Matrix, FileName, Form )
    class(ClassMatrix)        , intent(in) :: Matrix
    character(len=*)          , intent(in) :: FileName
    character(len=*), optional, intent(in) :: Form

    character(len=*), parameter   :: HERE = "ClassMatrix::WriteClassMatrixToFile : "
    character(len=*), parameter   :: DEFAULT_FORM = "unformatted"
    integer                       :: iostat, uid, islash
    character(len=:), allocatable :: FileForm
    character(len=IOMSG_LENGTH)   :: iomsg

    if(present(Form))then
       allocate(FileForm,source=Form)
    else
       allocate(FileForm,source=DEFAULT_FORM)
    endif

    !***
    !***
    islash = index( trim(FileName), "/", BACK = .TRUE. )
    call system("mkdir -p "//FileName(:islash-1))
    !***
    !***

    open(NewUnit =  uid     , &
         File    =  FileName, &
         form    =  FileForm, &
         status  = "unknown", &
         action  = "write"  , &
         iostat  = iostat   , &
         iomsg   = iomsg    )
    if(iostat/=0)call Assert(HERE//"Open of '"//trim(FileName)//"' failed : "//trim(iomsg))
    call Matrix.Write(uid)
    close(uid)

  end subroutine WriteClassMatrixToFile


  !> Writes the relevant information of the matrix class to a unit.
  subroutine WriteClassMatrixToUnit(Matrix,OutputUnit)
    Class(ClassMatrix), intent(in) :: Matrix
    integer, optional , intent(in) :: OutputUnit
    !
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Writable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    integer           :: OutUnit 
    if( present(OutputUnit) )then
       OutUnit=OutputUnit
    else
       OutUnit = DEFAULT_OUTPUT_UNIT
    endif

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
    !
    select case (trim(FORM))
    case("FORMATTED")
       call WriteClassMatrixToFormattedUnit(Matrix,OutUnit)
    case("UNFORMATTED")
       call WriteClassMatrixToUnformattedUnit(Matrix,OutUnit)
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
    !
  end subroutine WriteClassMatrixToUnit


  subroutine WriteClassMatrixToFormattedUnit(Matrix,OutUnit)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: OutUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=*), parameter :: FORMAT_INTS="(*(x,i))"
    character(len=*), parameter :: FORMAT_MAT="(*"//DOUBLE_PRINT_FORMAT//")"
    !
    !.. Write the Matrix attributes and dimensions
    write(OutUnit      , &
         FORMAT_INTS   , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix.Pattern, &
         Matrix.NR     , &
         Matrix.NC     , &
         Matrix.NU     , &
         Matrix.NL     , &
         Matrix.NRmin  , &
         Matrix.NRmax  , &
         Matrix.NCmin  , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    !.. Write Matrix
    do j=Matrix.NCmin,Matrix.NCmax
       write(OutUnit      , &
            FORMAT_MAT    , &
            IOSTAT=iostat , &
            IOMSG =iomsg  ) &
            (Matrix.A(i,j),&
            i=Matrix.NRmin,Matrix.NRmax)
       if(iostat/=0)call Assert(iomsg)
    enddo
    !
  end subroutine WriteClassMatrixToFormattedUnit


  subroutine WriteClassMatrixToUnformattedUnit(Matrix,OutUnit)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: OutUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    !
    !.. Write the Matrix attributes and dimensions
    write(OutUnit      , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix.Pattern, &
         Matrix.NR     , &
         Matrix.NC     , &
         Matrix.NU     , &
         Matrix.NL     , &
         Matrix.NRmin  , &
         Matrix.NRmax  , &
         Matrix.NCmin  , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    !.. Write Matrix
    write(OutUnit      , &
         IOSTAT=iostat , &
         IOMSG =iomsg  ) &
         ((Matrix.A(i,j),&
         i=Matrix.NRmin,Matrix.NRmax),&
         j=Matrix.NCmin,Matrix.NCmax)
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine WriteClassMatrixToUnformattedUnit


  subroutine ReadClassMatrixFromFile( Matrix, FileName, Form, iostat )
    class(ClassMatrix)        , intent(out):: Matrix
    character(len=*)          , intent(in) :: FileName
    character(len=*), optional, intent(in) :: Form
    integer         , optional, intent(out):: iostat

    character(len=*), parameter   :: HERE = "ClassMatrix::ReadClassMatrixFromFile : "
    character(len=*), parameter   :: DEFAULT_FORM = "unformatted"
    integer                       :: uid
    character(len=:), allocatable :: FileForm
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    if(present(Form))then
       allocate(FileForm,source=Form)
    else
       allocate(FileForm,source=DEFAULT_FORM)
    endif
    !
    call Matrix.Free()
    open(NewUnit =  uid     , &
         File    =  FileName, &
         form    =  FileForm, &
         status  = "unknown", &
         action  = "read"   , &
         iostat  = iostat   , &
         iomsg   = iomsg    )
    if(iostat/=0) then
       call Assert(iomsg)
       return
    end if
    call Matrix.Read(uid)
    close(uid)

  end subroutine ReadClassMatrixFromFile


  !> Reads the relevant information of the matrix class from a unit.
  subroutine ReadClassMatrixFromUnit(Matrix,InputUnit)
    Class(ClassMatrix), intent(inout):: Matrix
    integer, optional , intent(in)   :: InputUnit
    !
    integer :: InUnit
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=FILE_NAME_LENGTH)     :: FileName
    character(len=16) :: Readable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    if(present(InputUnit))then
       InUnit = InputUnit
    else
       InUnit= DEFAULT_INPUT_UNIT
    end if
    call Matrix.Free()
    INQUIRE(&
         UNIT  = InUnit  , &
         NAME  = FileName, &
         OPENED= Opened  , &
         READ  = Readable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Input Unit is closed")
    if( trim(Readable) /= "YES" ) call Assert("Input Unit can't be read")
    !
    select case (trim(FORM))
    case("FORMATTED")
       call ReadClassMatrixFromFormattedUnit(Matrix,InUnit)
    case("UNFORMATTED")
       call ReadClassMatrixFromUnformattedUnit(Matrix,InUnit)
    case DEFAULT
       call Assert("Invalid Input Unit Format")
    end select
    !
  end subroutine ReadClassMatrixFromUnit
  !
  subroutine ReadClassMatrixFromFormattedUnit(Matrix,InUnit)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: InUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg=" "
    character(len=*), parameter :: FORMAT_INTS="(*(x,i))"
    character(len=*), parameter :: FORMAT_MAT="(*"//DOUBLE_PRINT_FORMAT//")"
    !
    !.. Read the Matrix attributes and dimensions
    read(InUnit,FMT=FORMAT_INTS, &
         IOSTAT=iostat   , &
         IOMSG=iomsg     ) &
         Matrix.Pattern  , &
         Matrix.NR       , &
         Matrix.NC       , &
         Matrix.NU       , &
         Matrix.NL       , &
         Matrix.NRmin    , &
         Matrix.NRmax    , &
         Matrix.NCmin    , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    call AllocateMatrix( Matrix.A  , &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    !
    !.. Read Matrix
    do j=Matrix.NCmin,Matrix.NCmax
       read(InUnit,FMT=FORMAT_MAT  , &
            IOSTAT=iostat , &
            IOMSG =iomsg  ) &
            (Matrix.A(i,j), &
            i=Matrix.NRmin,Matrix.NRmax)
       if(iostat/=0)call Assert(iomsg)
    enddo
    !
  end subroutine ReadClassMatrixFromFormattedUnit
  !
  subroutine ReadClassMatrixFromUnformattedUnit(Matrix,InUnit)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: InUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    ! 
    !.. Write the Matrix attributes and dimensions
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix.Pattern, &
         Matrix.NR     , &
         Matrix.NC     , &
         Matrix.NU     , &
         Matrix.NL     , &
         Matrix.NRmin  , &
         Matrix.NRmax  , &
         Matrix.NCmin  , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    call AllocateMatrix( Matrix.A  , &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    !
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG =iomsg  ) &
         ((Matrix.A(i,j),&
         i=Matrix.NRmin,Matrix.NRmax),&
         j=Matrix.NCmin,Matrix.NCmax)
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine ReadClassMatrixFromUnformattedUnit


  !> Set all the elements of the matrix equal to real number.
  subroutine AddDoubleToClassMatrix(Matrix,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: x
    if(.not.allocated(Matrix.A))then
       call Assert("Matrix not initialized")
    else
       Matrix.A=x
    endif
  end subroutine AddDoubleToClassMatrix
  !
  !> Set all the elements of the matrix equal to an integer number.
  subroutine AssignIntegerToClassMatrix(Matrix,i)
    Class(ClassMatrix), intent(inout) :: Matrix
    Integer           , intent(in)    :: i
    call AssignDoubleToClassMatrix(Matrix,dble(i))
  end subroutine AssignIntegerToClassMatrix

  !> Multiplies the matrix by a real number.
  subroutine Double_x_ClassMatrix(Matrix,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: x
    if(.not.allocated(Matrix.A))then
       call Assert("Matrix not initialized")
    else
       Matrix.A=Matrix.A*x
    endif
  end subroutine Double_x_ClassMatrix




  !> Multiplies on the left a ClassMatrix object by a
  !! double matrix
  subroutine ClassMatrix_x_Matrix(self, MatFactor, Side, MatFactorTRANS )
    Class(ClassMatrix), intent(inout) :: self
    DoublePrecision   , intent(in)    :: MatFactor(:,:)
    character(len=*),   intent(in)    :: Side
    character(len=*),   intent(in)    :: MatFactorTRANS
    type(ClassMatrix) :: B
    character(len=*), parameter :: HERE = "ClassMatrix::ClassMatrix_x_Matrix : "
    if(.not.allocated(self.A))then
       call Assert(HERE//"Matrix not initialized")
    else
       B=MatFactor
       call self.Multiply(B,Side,MatFactorTRANS)
    endif
  end subroutine ClassMatrix_x_Matrix


  !> Multiplies the matrices belonging to two different matrix classes
  subroutine ClassMatrix_x_ClassMatrix( MatrixA, MatrixB, Side, MatrixBType )
    Class(ClassMatrix), target, intent(inout) :: MatrixA
    Class(ClassMatrix), target, intent(in)    :: MatrixB
    character(len=*),           intent(in)    :: Side
    character(len=*),           intent(in)    :: MatrixBType
    !
    Class(ClassMatrix), pointer :: PtrA, PtrB
    character         :: TypeA, TypeB
    Type(ClassMatrix) :: MatrixC
    integer :: m, n, k, lda, ldb, ldc
    !
    if( ( MatrixBType .isnt. "N" ) .and. ( MatrixBType .isnt. "T" ) ) &
         call Assert( 'Improper matrix type for multiplication, it must be "N" ("n") or "T" ("t")' )
    !
    if( Side .is. "Right" )then
       PtrA  => MatrixA
       TypeA = "N"
       PtrB  => MatrixB
       TypeB =  MatrixBType
    elseif( Side .is. "Left" )then
       PtrA  => MatrixB
       TypeA =  MatrixBType
       PtrB  => MatrixA
       TypeB = "N"
    else
       call Assert( ' improper multiplication side, must be "right" or "left"' )
    endif

    lda = size( PtrA.A, 1 )
    if( TypeA .is. "N" )then
       m = size( PtrA.A, 1 )
       k = size( PtrA.A, 2 )
    else
       m = size( PtrA.A, 2 )
       k = size( PtrA.A, 1 )
    endif

    ldb = size( PtrB.A, 1 )
    if(TypeB .is. "N")then
       n = size( PtrB.A, 2 )
    else
       n = size( PtrB.A, 1 )
    endif

    ldc = m
    call MatrixC.InitFull( m, n )
    !
    call DGEMM( &
         TypeA, TypeB, m, n, k, 1.d0, &
         PtrA.A,    lda, &
         PtrB.A,    ldb, 0.d0, &
         MatrixC.A, ldc )    
    !
    MatrixA = MatrixC
    !
    call MatrixC.Free()
    !
  end subroutine ClassMatrix_x_ClassMatrix


  !> Assign a real number to every element of the ClassMatrix's matrix.
  subroutine AssignDoubleToClassMatrix(Matrix,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: x
    if(.not.allocated(Matrix.A))then
       call Assert("Matrix not initialized")
    else
       Matrix.A=x
    endif
  end subroutine AssignDoubleToClassMatrix


  !> Copies the information of a matrix class to another.
  subroutine CopyClassMatrixToClassMatrix(MatrixOut,MatrixInp)
    Class(ClassMatrix), intent(inout) :: MatrixOut
    type (ClassMatrix), intent(in)    :: MatrixInp
    integer :: LBR, UBR, LBC, UBC, iRow,iCol
    !
    call MatrixOut.Free()
    !
    MatrixOut.Pattern = MatrixInp.Pattern
    !
    MatrixOut.NR = MatrixInp.NR
    MatrixOut.NC = MatrixInp.NC
    MatrixOut.NL = MatrixInp.NL
    MatrixOut.NU = MatrixInp.NU
    !
    MatrixOut.NRMin = MatrixInp.NRMin
    MatrixOut.NRMax = MatrixInp.NRMax
    MatrixOut.NCMin = MatrixInp.NCMin
    MatrixOut.NCMax = MatrixInp.NCMax
    !
    LBR=LBOUND(MatrixInp.A,1)
    UBR=UBOUND(MatrixInp.A,1)
    LBC=LBOUND(MatrixInp.A,2)
    UBC=UBOUND(MatrixInp.A,2)
!!$    allocate(MatrixOut.A,source=MatrixInp.A)
    allocate(MatrixOut.A(LBR:UBR,LBC:UBC))
    MatrixOut.A=0.d0
    do iCol=LBC, UBC
       do iRow=LBR,UBR
          MatrixOut.A(iRow,iCol)=MatrixInp.A(iRow,iCol)
       enddo
    enddo
    !
  end subroutine CopyClassMatrixToClassMatrix


  !> Assign to a ClassMatrix with Full pattern,
  !! the content of a bidimensional array
  subroutine CopyMatrixToClassMatrix(MatrixOut,MatrixInp)
    Class(ClassMatrix), intent(inout) :: MatrixOut
    DoublePrecision   , intent(in)    :: MatrixInp(:,:)
    integer :: NR, NC
    call MatrixOut.Free()
    NR=UBOUND(MatrixInp,1)-LBOUND(MatrixInp,1)+1
    NC=UBOUND(MatrixInp,2)-LBOUND(MatrixInp,2)+1
    call MatrixOut.initFull(NR,NC)
    MatrixOut.A=MAtrixInp
  end subroutine CopyMatrixToClassMatrix


  !> Assign to a ClassMatrix with Full pattern,
  !! the content of a bidimensional array
  subroutine CopyVectorToClassMatrix(MatrixOut,VectorInp)
    Class(ClassMatrix), intent(inout) :: MatrixOut
    DoublePrecision   , intent(in)    :: VectorInp(:)
    integer :: NR, NC
    call MatrixOut.Free()
    NR=UBOUND(VectorInp,1)-LBOUND(VectorInp,1)+1
    NC=1
    call MatrixOut.initFull(NR,NC)
    MatrixOut.A(:,1)=VectorInp
  end subroutine CopyVectorToClassMatrix


  !> Create and allocate a ClassMatrix object with Full 
  !! pattern, if it is not defined, and free it before, if it is.
  subroutine ClassMatrixInitFull( Matrix, NR, NC )
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR
    integer           , intent(in)    :: NC
    !
    call Matrix.Free()
    call SetMatrixNominalSize( Matrix, NR, NC )
    Matrix.NL=NR-1
    Matrix.NU=NC-1
    call SetMatrixPhysicalSize( Matrix, 1, NR, 1, NC )
    call AllocateMatrix( Matrix.A,   &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    Matrix.Pattern = MATRIX_PATTERN_FULL
    Matrix.A = 0.d0
    !
  end subroutine ClassMatrixInitFull


  !> Create and allocate a ClassMatrix object with Banded  pattern, if it is not defined, and free it before, if it is not.
  subroutine ClassMatrixInitBanded( Matrix, NR, NC, NL, NU )
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR
    integer           , intent(in)    :: NC
    integer           , intent(in)    :: NL
    integer           , intent(in)    :: NU
    !
    call Matrix.Free()
    call SetMatrixNominalSize( Matrix, NR, NC )
    Matrix.NL=max(0,min(NL,NR-1))
    Matrix.NU=max(0,min(NU,NC-1))
    call SetMatrixPhysicalSize( Matrix, -Matrix.NU, Matrix.NL, 1, Matrix.NC )
    call AllocateMatrix( Matrix.A,   &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    Matrix.Pattern = MATRIX_PATTERN_BANDED
    Matrix.A = 0.d0
    !
  end subroutine ClassMatrixInitBanded


  !> Set the rows and comlumns minimum and maximum indexes.
  subroutine SetMatrixPhysicalSize( Matrix, NRMin, NRMax, NCMin, NCMax )
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NRMin, NRMax, NCMin, NCMax
    Matrix.NRmin = NRMin
    Matrix.NRmax = NRMax
    Matrix.NCmin = NCMin
    Matrix.NCmax = NCMax
  end subroutine SetMatrixPhysicalSize


  !> Fetches the requested matrix element.
  DoublePrecision function ClassMatrixElement(Matrix,i,j) result(Element) 
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    Element=0.d0
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       Element = ClassMatrixFullElement(Matrix,i,j)
    case( MATRIX_PATTERN_BANDED )
       Element = ClassMatrixBandedElement(Matrix,i,j)
!!$    case( MATRIX_PATTERN_DIAGONAL )
!!$       Element = ClassMatrixDiagonalElement(Matrix,i,j)
    case DEFAULT
    end select
  end function ClassMatrixElement

  !> Sets the requested matrix element.
  subroutine ClassMatrixSetElement(Matrix,i,j,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    DoublePrecision   , intent(in)    :: x
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call ClassMatrixFullSetElement(Matrix,i,j,x)
    case( MATRIX_PATTERN_BANDED )
       call ClassMatrixBandedSetElement(Matrix,i,j,x)
    case DEFAULT
    end select
  end subroutine ClassMatrixSetElement

  !> Sets the requested matrix element when the array is in Full representation.
  subroutine ClassMatrixFullSetElement(Matrix,i,j,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    DoublePrecision   , intent(in)    :: x
    call CheckIndexBounds(Matrix,i,j)
    Matrix.A(i,j)=x
  end subroutine ClassMatrixFullSetElement

  !> Sets the requested matrix element when the array is in Banded representation.
  subroutine ClassMatrixBandedSetElement(Matrix,i,j,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    DoublePrecision   , intent(in)    :: x
    integer :: iPhys
    call CheckIndexBounds(Matrix,i,j)
    iPhys=i-j
    if( iPhys < -Matrix.NU .or. iPhys > Matrix.NL )&
         call Assert("Invalid banded matrix index")
    Matrix.A(iPhys,j)=x
  end subroutine ClassMatrixBandedSetElement


  !> Checks that the matrix element indexed have valid values.
  subroutine CheckIndexBounds(Matrix,i,j)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    if( i < 1 .or. i > Matrix.NR )call Assert("Row index is off bound")
    if( j < 1 .or. j > Matrix.NC )call Assert("Column index is off bound")
  end subroutine CheckIndexBounds


  !> Fetches the requested matrix elements when the array is stored in Full representation.
  DoublePrecision function ClassMatrixFullElement(Matrix,i,j) result(Element)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    integer  :: iPhys,jPhys
    call CheckIndexBounds(Matrix,i,j)
    iPhys=Matrix.NRMin-1+i
    jPhys=Matrix.NCMin-1+j
    Element=Matrix.A(iPhys,jPhys)
  end function ClassMatrixFullElement


  !> Fetches the requested matrix elements when the array is stored in Banded representation.
  DoublePrecision function ClassMatrixBandedElement(Matrix,i,j) result(Element)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    integer :: iPhys
    call CheckIndexBounds(Matrix,i,j)
    iPhys=i-j
    if( iPhys >= -Matrix.NU .and. iPhys <= Matrix.NL )then
       Element=Matrix.A(iPhys,j)
    else
       Element=0.d0
    end if
  end function ClassMatrixBandedElement


  !> Sets the number of rows and columns of the matrix.
  subroutine SetMatrixNominalSize(Matrix,NR,NC)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR,NC
    if(NR<=0)call Assert("Invalid Number of Rows")
    if(NC<=0)call Assert("Invalid Number of Columns")
    Matrix.NR=NR
    Matrix.NC=NC
  end subroutine SetMatrixNominalSize


  !> Allocates a real matrix and initializes it to zero.
  subroutine AllocateDoubleMatrix(A,NRMin,NRMax,NCMin,NCMax)
    DoublePrecision, allocatable, intent(inout) :: A(:,:)
    integer                     , intent(in)    :: NRMin, NRmax
    integer                     , intent(in)    :: NCMin, NCmax
    integer :: Stat=0
    character(len=IOMSG_LENGTH) :: iomsg=" "
    if( NRMax < NRmin .or. NCMax < NCMin )return
    if(allocated(A))deallocate(A)
    allocate(A(NRMin:NRMax,NCMin:NCMax),STAT=Stat,ERRMSG=iomsg)
    if(Stat/=0)call Assert(iomsg)
    A=0.d0
  end subroutine AllocateDoubleMatrix

  !> Allocates a complex matrix and initializes it to zero.
  subroutine AllocateComplexMatrix(A,NRMin,NRMax,NCMin,NCMax)
    Complex(kind(1d0)), allocatable, intent(inout) :: A(:,:)
    integer                     , intent(in)    :: NRMin, NRmax
    integer                     , intent(in)    :: NCMin, NCmax
    integer :: Stat=0
    character(len=IOMSG_LENGTH) :: iomsg=" "
    if( NRMax < NRmin .or. NCMax < NCMin )return
    if(allocated(A))deallocate(A)
    allocate(A(NRMin:NRMax,NCMin:NCMax),STAT=Stat,ERRMSG=iomsg)
    if(Stat/=0)call Assert(iomsg)
    A=0.d0
  end subroutine AllocateComplexMatrix


  !> Free allocatable space and sets to zero all the members of a ClassMatrix object.
  subroutine ClassMatrixFree( Matrix )
    Class(ClassMatrix), intent(inout) :: Matrix
    if(allocated(Matrix.A)) deallocate(Matrix.A)
    Matrix.NR=0
    Matrix.NC=0
    Matrix.Pattern=MATRIX_PATTERN_UNDEFINED
    Matrix.NU=0
    Matrix.NL=0
  end subroutine ClassMatrixFree


  !> Free allocatable space and sets to zero all the members of a ClassMatrix object.
  subroutine ClassMatrixFinalize( Matrix )
    type(ClassMatrix), intent(inout) :: Matrix
    call Matrix.Free()
  end subroutine ClassMatrixFinalize


  !> Retrieves the number of rows of the matrix in ClassMatrix.
  integer function ClassMatrixNRows( Matrix ) result( NRows )
    Class(ClassMatrix), intent(in) :: Matrix
    NRows = Matrix.NR
  end function ClassMatrixNRows
  !
  !> Retrieves the number of columns of the matrix in ClassMatrix.
  integer function ClassMatrixNColumns( Matrix ) result( NColumns )
    Class(ClassMatrix), intent(in) :: Matrix
    NColumns = Matrix.NC
  end function ClassMatrixNColumns


  !> Sort eigenvalues and and its corresponding eigenvectors in ascending order of the former.
  subroutine ClassSpectralResolutionSort( SpecRes )
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    !
    integer, parameter :: INCREASING_ORDER_SORT_VECTOR      =  2
    integer, parameter :: INCREASING_ORDER_DONT_SORT_VECTOR =  1
    integer, parameter :: DECREASING_ORDER_DONT_SORT_VECTOR = -1
    integer, parameter :: DECREASING_ORDER_SORT_VECTOR      = -2
    integer, parameter :: SORT_CRITERION = INCREASING_ORDER_DONT_SORT_VECTOR
    !
    DoublePrecision, allocatable :: EvalRe(:)
    DoublePrecision, allocatable :: TmpEigenvalues(:)
    DoublePrecision, allocatable :: TmpEigenvectors(:,:)
    integer        , allocatable :: IPerm(:)
    integer :: Neval,dim,iEval,ier
    !
    Neval=SpecRes.Neval()
    dim=SpecRes.Size()
    if(Neval<=0)return
    allocate(EvalRe(Neval),IPerm(Neval))
    IPerm=[(iEval,iEval=1,Neval)]
    EvalRe=dble(SpecRes.Eigenvalues(1:Neval))
    !
    call DPSORT(EvalRe,Neval,IPerm,SORT_CRITERION,ier)
    if(ier/=0)call Assert("DPSORT: invalid input")
    deallocate(EvalRe)
    !
    allocate(TmpEigenvalues(Neval),TmpEigenvectors(dim,Neval))
    do iEval = 1, Neval
       TmpEigenvalues(   iEval)=SpecRes.Eigenvalues(   IPerm(iEval))
       TmpEigenvectors(:,iEval)=SpecRes.Eigenvectors(:,IPerm(iEval))
    enddo
    deallocate(IPerm)
    SpecRes.Eigenvalues(   1:Neval)=TmpEigenvalues
    deallocate(TmpEigenvalues)
    SpecRes.Eigenvectors(:,1:Neval)=TmpEigenvectors
    deallocate(TmpEigenvectors)
    !
  end subroutine ClassSpectralResolutionSort


  subroutine ClassMatrixAddRows( self, NNewRows, where_, After )
    class(ClassMatrix), intent(inout) :: self
    integer           , intent(in)    :: NNewRows
    !> where is to be either "START" or "END"
    character(len=*), optional , intent(in) :: where_
    !> After can be any row in the input matrix.
    !! If specified, the function adds NNewRows between
    !! After and After+1.
    !! One can specify either where_ or After, but not both.
    integer         , optional , intent(in) :: After
    type(ClassMatrix) :: mat
    if(NNewRows<=0)return
    if(present(where_))then

       if(present(After))call Assert("In AddRows, you can't specify both 'where_' and 'After'")

       if(where_.is."START")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows()+NNewRows, self.NColumns() )
             mat.A(NNewRows+1:,:)=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             !.. The absolute indexing of the matrix changes
             !   therefore it is better to make a full copy
             !   even if the allocated space wouldn't change
             call mat.initBanded(&
                  self.NRows()+NNewRows,&
                  self.NColumns(),&
                  self.LowerBandwidth()+NNewRows,&
                  self.UpperBandwidth()-NNewRows)
             mat.A=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select
       elseif(where_.is."END")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows()+NNewRows, self.NColumns() )
             mat.A(:self.NRows(),:)=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             !.. In this case both the absolute indexing and
             !   the storage size do not change
             !..
             self.NR=self.NR+NNewRows
          case DEFAULT 
             call Assert("Unrecognized matrix pattern")
          end select
       endif

    elseif( present(After) )then

       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows()+NNewRows, self.NColumns() )
          mat.A(1:After,:)=self.A(1:After,:)
          mat.A(After+1+NNewRows:,:)=self.A(After+1:,:)
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call Assert("'After' option not implemented in AddRows for banded matrices")
       case DEFAULT 
          call Assert("Unrecognized matrix pattern")
       end select

    else

       call Assert("In AddRows, either 'where_' or 'After' must be specified")

    endif

  end subroutine ClassMatrixAddRows



  subroutine ClassMatrixAddColumns( self, NNewCols, where_, After )
    class(ClassMatrix), intent(inout) :: self
    integer           , intent(in) :: NNewCols
    character(len=*), optional, intent(in) :: where_
    integer, optional , intent(in) :: After
    type(ClassMatrix) :: mat
    if(NNewCols<=0)return
    if(present(where_))then

       if(present(After))call Assert("In AddColumns, you can't specify both 'where_' and 'After'")

       if(where_.is."START")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows(), self.NColumns()+Nnewcols )
             mat.A(:,Nnewcols+1:)=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             call mat.initBanded(&
                  self.NRows(),&
                  self.NColumns()+Nnewcols,&
                  self.LowerBandwidth()-NNewCols,&
                  self.UpperBandwidth()+NNewCols)
             mat.A(:,Nnewcols+1:)=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select
       elseif(where_.is."END")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows(), self.NColumns()+Nnewcols )
             mat.A(:,:self.NColumns())=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             call mat.initBanded(&
                  self.NRows(),&
                  self.NColumns()+Nnewcols,&
                  self.LowerBandwidth(),&
                  self.UpperBandwidth())
             mat.A(:,:self.NColumns())=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select
       endif

    elseif( present(After) )then

       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows(), self.NColumns()+NNewCols )
          mat.A(:,1:After)=self.A(:,1:After)
          mat.A(:,After+1+NNewCols:)=self.A(:,After+1:)
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call Assert("'After' option not implemented in AddColumns for banded matrices")
       case DEFAULT 
          call Assert("Unrecognized matrix pattern")
       end select

    else

       call Assert("In AddColumns, either 'where_' or 'After' must be specified")

    endif

  end subroutine ClassMatrixAddColumns


  subroutine ClassMatrixRemoveRows( self, NKillRows, where_, After )
    class(ClassMatrix)        , intent(inout) :: self
    integer                   , intent(in)    :: NKillRows
    character(len=*), optional, intent(in)    :: where_
    integer         , optional, intent(in)    :: After

    character(len=*), parameter :: HERE="ClassMatrix::RemoveRows : "
    type(ClassMatrix) :: mat
    integer :: iRow, iColumn

    if(NKillRows<=0) return
    if(NKillRows>self.NRows()) call Assert(HERE//"too many rows to be eliminated")

    if(present(where_))then

       if(where_.is."START")then

          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows()-NKillRows, self.NColumns() )
             do iColumn=1,self.NColumns()
                do iRow = 1, self.NRows() - NkillRows
                   mat.A(iRow,iColumn)=self.A(NKillRows+iRow,iColumn)
                enddo
             enddo
             self=mat
          case( MATRIX_PATTERN_BANDED )
             !set to zero the elements in the rows eliminated
             do iRow = 1, NKillRows
                do iColumn=&
                     max(1,iRow-self.LowerBandwidth()),&
                     min(self.NColumns(),iRow+self.UpperBandwidth())
                   call self.SetElement(iRow,iColumn,0.d0)
                enddo
             enddo
             !.. The absolute indexing of the matrix changes
             !   therefore it is better to make a full copy
             !   even if the allocated space wouldn't change
             call mat.initBanded(&
                  self.NRows()-NKillRows,&
                  self.NColumns(),&
                  self.LowerBandwidth()-NKillRows,&
                  self.UpperBandwidth()+NKillRows)
             mat.A=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select

       elseif(where_.is."END")then

          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows()-NKillRows, self.NColumns() )
             do iColumn =1, self.NColumns()
                do iRow = 1, self.NRows()-NKillRows
                   mat.A(iRow,iColumn)=self.A(iRow,iColumn)
                enddo
             enddo
             self=mat
          case( MATRIX_PATTERN_BANDED )
             !.. In this case both the absolute indexing and
             !   the storage size do not change
             !..
             !set to zero the elements in the rows eliminated
             do iRow = self.NRows()-NKillRows+1,self.NRows()
                do iColumn=&
                     max(1,iRow-self.LowerBandwidth()),&
                     min(self.NColumns(),iRow+self.UpperBandwidth())
                   call self.SetElement(iRow,iColumn,0.d0)
                enddo
             enddo
             self.NR=self.NR-NKillRows
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select

       endif

    elseif(present(After))then

       if(NKillRows>self.NRows()-After) call Assert(HERE//"too many rows to be eliminated")

       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows()-NKillRows, self.NColumns() )
          do iColumn=1,self.NColumns()
             do iRow=1,After
                mat.A(iRow,iColumn)=self.A(iRow,iColumn)
             enddo
             do iRow=After+1,self.NRows()-NKillRows
                mat.A(iRow,iColumn)=self.A(iRow+NKillRows,iColumn)
             enddo
          enddo
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call Assert("'After' option not implemented in RemoveRows for banded matrices")
       case DEFAULT 
          call Assert("Unrecognized matrix pattern")
       end select

    else

       call Assert("In RemoveRows, either 'where_' or 'After' must be specified")

    endif
  end subroutine ClassMatrixRemoveRows


  subroutine ClassMatrixRemoveColumns( self, Nkillcols, where_, After )
    class(ClassMatrix)        , intent(inout) :: self
    integer                   , intent(in)    :: Nkillcols
    character(len=*), optional, intent(in)    :: where_
    integer         , optional, intent(in)    :: After

    character(len=*), parameter :: HERE="ClassMatrix::RemoveColumns : "
    type(ClassMatrix) :: mat
    integer :: i, j

    if( Nkillcols <= 0 )return
    if( NKillCols > self.NColumns() ) call Assert(HERE//"too many Columns to be eliminated")

    if( present(where_) )then 

       if(where_.is."START")then

          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows(), self.NColumns()-Nkillcols )
! This line introduces core dumped segmetation error.
!!$             mat.A=self.A(:,Nkillcols+1:)
do j = 1, mat.NColumns()
   !
   do i = 1, mat.NRows()
      !
      mat.A(i,j) = self.A(i,Nkillcols+j)
      !
   end do
   !
end do
             self=mat
          case( MATRIX_PATTERN_BANDED )
             call mat.initBanded(&
                  self.NRows(),&
                  self.NColumns()-Nkillcols,&
                  self.LowerBandwidth()+NKillCols,&
                  self.UpperBandwidth()-NKillCols)
             mat.A=self.A(:,Nkillcols+1:)
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select

       elseif(where_.is."END")then

          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows(), self.NColumns()-Nkillcols )
             ! This lines causes segmentation fault when many basis are present.
!!$             mat.A=self.A(:,:mat.NColumns())
             do j = 1, mat.NColumns()
                !
                do i = 1, mat.NRows()
                   !
                   mat.A(i,j) = self.A(i,j)
                   !
                end do
                !
             end do
             self=mat
          case( MATRIX_PATTERN_BANDED )
             call mat.initBanded(&
                  self.NRows(),&
                  self.NColumns()-Nkillcols,&
                  self.LowerBandwidth(),&
                  self.UpperBandwidth())
             mat.A=self.A(:,:mat.NColumns())
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select

       endif

    elseif( present( After ) )then

       if( NKillCols > self.NColumns() - After ) call Assert(HERE//"too many Columns to be eliminated")

       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows(), self.NColumns()-NKillCols )
          mat.A(:,1:After)=self.A(:,1:After)
          ! This line causes errors when many basis functions are present.
!!$          mat.A(:,After+1:)=self.A(:,After+1+NKillCols:)
          do j = After+1, mat.NColumns()
             !
             do i = 1, mat.NRows()
                !
                mat.A(i,j) = self.A(i,j+NKillCols)
                !
             end do
             !
          end do
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call Assert("'After' option not implemented in RemoveCols for banded matrices")
       case DEFAULT 
          call Assert(HERE//"Unrecognized matrix pattern")
       end select

    else

       call Assert("In RemoveRows, either 'where_' or 'After' must be specified")

    endif

  end subroutine ClassMatrixRemoveColumns


  !> Sort eigenvalues and and its corresponding eigenvectors in ascending order of the former's real part.
  subroutine ClassComplexSpectralResolutionSort( SpecRes )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    !
    integer, parameter :: INCREASING_ORDER_SORT_VECTOR      =  2
    integer, parameter :: INCREASING_ORDER_DONT_SORT_VECTOR =  1
    integer, parameter :: DECREASING_ORDER_DONT_SORT_VECTOR = -1
    integer, parameter :: DECREASING_ORDER_SORT_VECTOR      = -2
    integer, parameter :: SORT_CRITERION = INCREASING_ORDER_DONT_SORT_VECTOR
    !
    DoublePrecision   , allocatable :: EvalRe(:)
    Complex(kind(1d0)), allocatable :: TmpEigenvalues(:)
    Complex(kind(1d0)), allocatable :: TmpEigenvectors(:,:)
    integer           , allocatable :: IPerm(:)
    integer :: Neval,dim,iEval,ier
    !
    Neval=SpecRes.Neval()
    dim=SpecRes.Size()
    if(Neval<=0)return
    allocate(EvalRe(Neval),IPerm(Neval))
    IPerm=[(iEval,iEval=1,Neval)]
    EvalRe=dble(SpecRes.Eigenvalues(1:Neval))
    !
    call DPSORT(EvalRe,Neval,IPerm,SORT_CRITERION,ier)
    if(ier/=0)call Assert("DPSORT: invalid input")
    deallocate(EvalRe)
    !
    allocate(TmpEigenvalues(Neval),TmpEigenvectors(dim,Neval))
    do iEval = 1, Neval
       TmpEigenvalues(   iEval)=SpecRes.Eigenvalues(   IPerm(iEval))
       TmpEigenvectors(:,iEval)=SpecRes.RightEigenvectors(:,IPerm(iEval))
    enddo
    deallocate(IPerm)
    SpecRes.Eigenvalues(   1:Neval)=TmpEigenvalues
    deallocate(TmpEigenvalues)
    SpecRes.RightEigenvectors(:,1:Neval)=TmpEigenvectors
    deallocate(TmpEigenvectors)
    !
  end subroutine ClassComplexSpectralResolutionSort


  !> Solves the general eigenvalue problem:
  !!\f[
  !! \left(\mathbb{A}-\alpha\cdot\mathbb{S}\right)\cdot\mathbb{c}=0
  !!\f]
  !! Assumes Metric (\f$\mathbb{S}\f$) is definite positive.
  subroutine GeneralEigenValueSolver( Matrix, Metric, SpecRes )
    Class(ClassMatrix),            intent(in)  :: Matrix
    type(ClassMatrix),             intent(in)  :: Metric
    type(ClassSpectralResolution), intent(out) :: SpecRes
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call FullGeneralEigenvalueSolver( Matrix, Metric, SpecRes )
    case( MATRIX_PATTERN_BANDED )
       call BandGeneralEigenvalueSolver( Matrix, Metric, SpecRes )
    case DEFAULT 
       call Assert('Error: unrecognized matrix pattern')
    end select
    !
  end subroutine GeneralEigenValueSolver


  subroutine GeneralEigenValueSolverConditionNumber( Matrix, Metric, SpecRes, CN, MaxEigVal )
    Class(ClassMatrix),            intent(inout)  :: Matrix
    type(ClassMatrix),             intent(in)     :: Metric
    type(ClassSpectralResolution), intent(out)    :: SpecRes
    real(kind(1d0)),               intent(in)     :: CN
    real(kind(1d0)),  optional   , intent(in)     :: MaxEigVal
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call FullGeneralEigenvalueSolverConditionNumber( Matrix, Metric, SpecRes, CN, MaxEigVal )
    case( MATRIX_PATTERN_BANDED )
       call Assert( "The general real matrix diagonalization using a condition number is only implemented for full pattern matrices." )
    case DEFAULT 
       call Assert('Error: unrecognized matrix pattern')
    end select
    !
  end subroutine GeneralEigenValueSolverConditionNumber



  subroutine FullGeneralEigenvalueSolverConditionNumber( Matrix, Metric, SpecRes, CN, MaxEigVal )
    !
    Class(ClassMatrix),            intent(inout)  :: Matrix
    type(ClassMatrix),             intent(in)     :: Metric
    type(ClassSpectralResolution), intent(out)    :: SpecRes 
    real(kind(1d0)),               intent(in)     :: CN
    real(kind(1d0)),  optional   , intent(in)     :: MaxEigVal
    !
    type(ClassSpectralResolution) :: MetricSpecRes
    real(kind(1d0)), allocatable :: MetricEigVal(:), MetricEigVec(:,:), ArrayMetricEigVec(:,:), FirstMatTransf(:,:), PristineMat(:,:), SecondMatTransf(:,:), EigVec(:,:), NewEigVec(:,:)
    integer :: n, i, NumZeroEigVal, m, IndexBigEigValCut
    type(ClassMatrix) :: NewMatrix
    real(kind(1d0)) :: HighestEigVal, HighestAllowedEigVal
    real(kind(1d0)), parameter :: InnerThr = 1.d-6
    !
    if ( present(MaxEigVal) ) then
       HighestEigVal = MaxEigVal
    else
       HighestEigVal = huge(1d0)
    end if
    !
    call Metric.Diagonalize( MetricSpecRes )
    !
    call MetricSpecRes.Fetch( MetricEigVal )
    call MetricSpecRes.Fetch( MetricEigVec )
    !
    call MetricSpecRes.Free()
    n = size(MetricEigVal)
    !
    NumZeroEigVal = 0
    IndexBigEigValCut = 0
    !
    HighestAllowedEigVal = min( MetricEigVal(n),HighestEigVal)
    ! Supposes eigenvalues in ascending order.
    do i = 1, n
       !
       if ( MetricEigVal(i)/HighestAllowedEigVal < CN ) then
          NumZeroEigVal = NumZeroEigVal + 1
       end if
       if ( (MetricEigVal(i) > (HighestAllowedEigVal + InnerThr)) .and. &
            (IndexBigEigValCut == 0) ) then
          IndexBigEigValCut = i-1
       end if
       !
    end do
    write(output_unit,*) "NumZeroEigVal", NumZeroEigVal
    write(output_unit,*) "HighestAllowedEigVal", HighestAllowedEigVal
    !
    if ( IndexBigEigValCut == 0 ) IndexBigEigValCut = n
    !
    m = n - NumZeroEigVal - (n-IndexBigEigValCut)
    !
    MetricEigVec(:,1:NumZeroEigVal) = 0.d0
    if ( IndexBigEigValCut < n ) then
       MetricEigVec(:,IndexBigEigValCut+1:) = 0.d0
    end if
    !
    do i = NumZeroEigVal+1, IndexBigEigValCut
       MetricEigVec(:,i) = MetricEigVec(:,i)/sqrt(MetricEigVal(i))
    enddo
    !
    allocate( ArrayMetricEigVec(n,m) )
    ArrayMetricEigVec = MetricEigVec(:,NumZeroEigVal+1:IndexBigEigValCut)
    !
    allocate(FirstMatTransf(n,m))
    !
    FirstMatTransf = 0.d0
    !
    call Matrix.FetchMatrix( PristineMat )
    !
    call dgemm( "N", "N", n, m, n, 1.d0, PristineMat, n, &
         ArrayMetricEigVec, n, 0.d0, FirstMatTransf, n )
    !
    deallocate(PristineMat)
    !
    allocate(SecondMatTransf(m,m))
    !
    call dgemm( "T", "N", m, m, n, 1.d0, ArrayMetricEigVec, n, &
         FirstMatTransf, n, 0.d0, SecondMatTransf, m )
    !
    deallocate(FirstMatTransf)
    !
    NewMatrix = SecondMatTransf
    call NewMatrix.Diagonalize( SpecRes )
    !
    deallocate(SecondMatTransf)
    !
    call SpecRes.Fetch( EigVec )
    !
    allocate( NewEigVec(n,m) )
    NewEigVec = 0.d0
    !
    call dgemm( "N", "N", n, m, m, 1.d0, ArrayMetricEigVec, n, &
         EigVec, m, 0.d0, NewEigVec, n )
    !
    deallocate( SpecRes.Eigenvectors )
    allocate( SpecRes.Eigenvectors, source = NewEigVec )
    SpecRes.Dim = size(NewEigVec,1)
    SpecRes.NEigenvalues = size(NewEigVec,2)
    deallocate( MetricEigVec, ArrayMetricEigVec )
    deallocate( EigVec )
    !
  end subroutine FullGeneralEigenvalueSolverConditionNumber





  subroutine ClassMatrixConditionMatrix( Matrix, CN, TransfMat )
    !
    Class(ClassMatrix),            intent(inout)  :: Matrix
    real(kind(1d0)),               intent(in)     :: CN
    type(ClassMatrix),             intent(out)    :: TransfMat
    !
    type(ClassSpectralResolution) :: SpecRes
    real(kind(1d0)), allocatable :: EigVal(:), EigVec(:,:)
    integer :: n, i, NumZeroEigVal
    !
    call Matrix.Diagonalize( SpecRes )
    !
    call SpecRes.Fetch( EigVal )
    call SpecRes.Fetch( EigVec )
    !
    call SpecRes.Free()
    n = size(EigVal)
    !
    NumZeroEigVal = 0
    !
    ! Supposes eigenvalues in ascending order.
    do i = 1, n
       if ( EigVal(i)/EigVal(n) < CN ) then
          NumZeroEigVal = NumZeroEigVal + 1
       end if
    end do
    write(*,*) "NumZeroEigVal", NumZeroEigVal
    !
    EigVec(:,1:NumZeroEigVal) = 0.d0
    !
    do i = NumZeroEigVal+1, n
       EigVec(:,i) = EigVec(:,i)/sqrt(EigVal(i))
    enddo
    !
    TransfMat = EigVec(:,NumZeroEigVal+1:)
    !
  end subroutine ClassMatrixConditionMatrix
       





  !> Solves the general eigenvalue problem when the matrices are in 
  !! Full representation. On output, the Metric S is overwritten by the
  !! triangular factor U or L from the Cholesky factorization 
  !! S=U**T*U or S=L*L**T
  subroutine FullGeneralEigenvalueSolver ( Matrix, Metric, SpecRes )
    Class(ClassMatrix),            intent(in)  :: Matrix
    type(ClassMatrix),             intent(in)  :: Metric
    type(ClassSpectralResolution), intent(out) :: SpecRes 
    !
    integer, parameter :: Ax_EQ_Lambda_Bx = 1
    integer, parameter :: PROBLEM_TYPE = Ax_EQ_Lambda_Bx
    !
    integer :: n, lda, ldb, INFO
    integer(kind=8) :: LWORK
    DoublePrecision, allocatable :: WORK(:)
    !
    n = Matrix.NRows()
    call SpecRes.Init( n )
    !
    SpecRes.EigenVectors = Matrix.A
    !
    lda = n
    ldb = n
    !
    !.. Optimal work-size query
    LWORK=-1
    allocate( WORK( 1 ) )
    call DSYGV( PROBLEM_TYPE, 'V', 'U', n, &
         SpecRes.EigenVectors, lda, Metric.A, ldb, &
         SpecRes.EigenValues, WORK, LWORK, INFO )
    LWORK = max(3*n,int(WORK(1)))

    deallocate( WORK )
    allocate( WORK( LWORK ) )
    !
    call DSYGV( PROBLEM_TYPE, 'V', 'U', n, &
         SpecRes.EigenVectors, lda, &
         Metric.A, ldb, &
         SpecRes.EigenValues, &
         WORK, LWORK, INFO )
    deallocate( WORK )
    !
    if(info<0)call ErrorMessage(&
         "DSYGV: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0.and.info<=n)call ErrorMessage(&
         "DSYGV: DSYEV failed to converge, and "//AlphabeticNumber(info)//" off-diagonal"//&
         "elements of an intermediate tridiagonal did not converge to zero.")
    if(info>n.and.info<=2*n)call ErrorMessage(&
         "DSYGV: the leading minor of order "//AlphabeticNumber(info-n)//&
         " of the metric is not positive-definite."//&
         "The factorization of the metric could not be completed and no"//&
         "eigenvalues or eigenvectors were computed.")
    !
    if(info==0)then
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
    endif
    !
  end subroutine FullGeneralEigenvalueSolver

  !> Solves the general eigenvalue problem when the matrices are in Banded representation.
  subroutine BandGeneralEigenvalueSolver ( Matrix, Metric, SpecRes )
    Class(ClassMatrix),            intent(in)  :: Matrix
    type(ClassMatrix),             intent(in)  :: Metric
    type(ClassSpectralResolution), intent(out) :: SpecRes
    !
    integer :: n, ka, kb, ldab, ldbb, ldz, INFO, liwork
    integer(kind=8) :: lwork
    DoublePrecision, allocatable :: work(:)
    DoublePrecision, allocatable :: ABand(:,:), BBand(:,:)
    !
    n  = size(Matrix.A,2)
    call SpecRes.Init(n)
    !
    ka = Matrix.NU
    kb = Matrix.NU
    lwork = 2*n**2 + 5*n + 1
    if ( (abs(lwork) >= huge(n)) .or. (lwork < 0) ) call Assert( &
         'The lwork variable exceeds the integer kind=4 limits for DSYEVD.' )
    liwork = 5*n+3
    ldab = ka+1
    ldbb = ldab
    ldz  = n
    !
    allocate( WORK( 3*n ) )
    allocate( ABand, source = Matrix.A(-ka:0,:) )
    allocate( BBand, source = Metric.A(-kb:0,:) )
    call DSBGV( 'V', 'U', n, ka, kb, &
         ABand, ldab, &
         BBand, ldbb, &
         SpecRes.EigenValues, SpecRes.EigenVectors, ldz, &
         WORK, INFO )  
    !
    if(info<0)call ErrorMessage(&
         "DSBGV: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0.and.info<=n)call ErrorMessage(&
         "DSBGV: DSYEV failed to converge, and "//AlphabeticNumber(info)//" off-diagonal"//&
         "elements of an intermediate tridiagonal did not converge to zero.")
    if(info>n.and.info<=2*n)call ErrorMessage(&
         "DSBGV: DPBSTF returned info="//AlphabeticNumber(info-n)//" which means that the metric is not positive-definite."//&
         "The factorization of the metric could not be completed and no eigenvalues or eigenvectors were computed.")
    !
    if(info==0)then
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
    endif
    !
    deallocate(work,ABand,BBand)
    !
  end subroutine BandGeneralEigenvalueSolver


  !> Solves the eigenvalues eigenvectors problem for a ClassMatrix's matrix and stores the results in ClassSpectralResolution.
  subroutine EigenvalueSolver( Matrix, SpecRes )
    Class(ClassMatrix),            intent(in)  :: Matrix
    type(ClassSpectralResolution), intent(out) :: SpecRes 
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call FullEigenvalueSolver( Matrix, SpecRes )
    case( MATRIX_PATTERN_BANDED )
       call BandEigenvalueSolver( Matrix, SpecRes )
    case DEFAULT 
       call Assert('Error: non-proper matrix pattern')
    end select
    !
  end subroutine EigenvalueSolver

  !> Solves the eigenvalues eigenvectors problem for a ClassMatrix's matrix in Full representation and stores the results in ClassSpectralResolution.
  subroutine FullEigenvalueSolver( Matrix, SpecRes ) 
    Class(ClassMatrix),            intent(in)    :: Matrix
    type(ClassSpectralResolution), intent(out) :: SpecRes
    !
    integer :: lda, n
    integer(kind=8) :: lwork, liwork, INFO
    DoublePrecision, allocatable :: work(:)
    integer(kind=8) , allocatable :: iwork(:)
    !
    n = Matrix.NRows()
    call SpecRes.Init( n )
    SpecRes.EigenVectors = Matrix.A
    lda = n
allocate( work(1), iwork(1) )
    call DSYEVD( 'V', 'U', n, &
         SpecRes.EigenVectors, lda, &
         SpecRes.EigenValues, &
         work, -1, iwork, -1, INFO )
    lwork = max(int(work(1)),2*n**2+6*n+1)
    if ( (abs(lwork) >= huge(n)) .or. (lwork < 0) ) call Assert( &
         'The lwork variable exceeds the integer kind=4 limits for DSYEVD.' )
    liwork = max(int(iwork(1)),5*n+3)
    deallocate( work, iwork )
    allocate( work(lwork), iwork(liwork) )
    !
    call DSYEVD( 'V', 'U', n, &
         SpecRes.EigenVectors, lda, &
         SpecRes.EigenValues, &
         work, lwork, iwork, liwork, INFO )
    !
    if(info<0)call ErrorMessage(&
         "DSYEVD: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0.and.info<=n)call ErrorMessage(&
         "DSYEVD: failed to converge since "//AlphabeticNumber(info)//" off-diagonal"//&
         "elements of an intermediate tridiagonal did not converge to zero.")
    if(info>n)call ErrorMessage("DSYEVD: failed to compute some eigenvalues")
    !
    if(info==0)then
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
    endif
    !
    deallocate(work,iwork)
    !
  end subroutine FullEigenvalueSolver

  !> Solves the eigenvalues eigenvectors problem for a ClassMatrix's matrix in Banded representation and stores the results in ClassSpectralResolution.
  subroutine BandEigenvalueSolver( Matrix, SpecRes )
    Class(ClassMatrix),            intent(in)    :: Matrix
    type(ClassSpectralResolution), intent(out) :: SpecRes
    !
    integer :: n, lda, kd, ldz, INFO
    DoublePrecision, allocatable :: work(:)
    !
    n = Matrix.NRows()
    call SpecRes.Init( n )
    SpecRes.EigenVectors = Matrix.A
    kd = Matrix.NU
    lda = size(Matrix.A,1)
    ldz = n
    !
    allocate( work( 3*n - 2 ) )
    !
    call DSBEV( 'V', 'U', n, kd,&
         Matrix.A, lda,&
         SpecRes.EigenValues,&
         SpecRes.EigenVectors, ldz,&
         work, INFO )
    if(info<0)call ErrorMessage(&
         "DSBEV: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0.and.info<=n)call ErrorMessage(&
         "DSBEV: failed to converge since "//AlphabeticNumber(info)//" off-diagonal"//&
         "elements of an intermediate tridiagonal did not converge to zero.")
    !
    if(info==0)then
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
    endif
    !
    deallocate( work )
    !
  end subroutine BandEigenvalueSolver


  !> Initializes the spectral resolution class for the case when the matrix of eigenvectors is squared and the number of rows is equal to the number of eigenvalues.
  subroutine InitSpectralResolutionFull( SpecRes, n )
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    integer,                        intent(in)    :: n
    !
    SpecRes.NEigenvalues = n
    SpecRes.Dim = n
    !
    if(allocated(SpecRes.EigenValues))then
       if(size(SpecRes.Eigenvalues,1)/=n)then
          deallocate( SpecRes.EigenValues )
          allocate( SpecRes.EigenValues( n ))
       endif
    else
       allocate( SpecRes.EigenValues( n ))
    endif
    SpecRes.EigenValues = 0.d0
    !
    if(allocated(SpecRes.EigenVectors))then
       if(  size(SpecRes.Eigenvectors,1) /= n  .or. &
            size(SpecRes.Eigenvectors,2) /= n  )then
          deallocate(SpecRes.EigenVectors)
          allocate(SpecRes.EigenVectors( n, n ) )
       endif
    else
       allocate(SpecRes.EigenVectors( n, n ) )
    endif
    SpecRes.EigenVectors = 0.d0
    !
  end subroutine InitSpectralResolutionFull



  !> Initializes the spectral resolution class for the case when the matrix 
  !! of eigenvectors is rectangular and the number of rows larger than the 
  !! number of eigenvalues.
  subroutine InitSpectralResolutionReduced( SpecRes, Dim, NEigenvalues )
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    integer,                        intent(in)    :: Dim, NEigenvalues
    !
    SpecRes.Dim = Dim
    SpecRes.NEigenvalues = NEigenvalues
    !
    if(allocated(SpecRes.EigenValues))then
       if(size(SpecRes.Eigenvalues,1)/=NEigenvalues)then
          deallocate( SpecRes.EigenValues )
          allocate( SpecRes.EigenValues( NEigenvalues ))
       endif
    else
       allocate( SpecRes.EigenValues( NEigenvalues ))
    endif
    SpecRes.EigenValues = 0.d0
    !
    if(allocated(SpecRes.EigenVectors))then
       if(  size(SpecRes.Eigenvectors,1) /= Dim  .or. &
            size(SpecRes.Eigenvectors,2) /= NEigenvalues  )then
          deallocate(SpecRes.EigenVectors)
          allocate(SpecRes.EigenVectors( Dim, NEigenvalues ) )
       endif
    else
       allocate(SpecRes.EigenVectors( Dim, NEigenvalues ) )
    endif
    SpecRes.EigenVectors = 0.d0
    !
  end subroutine InitSpectralResolutionReduced

  !> Retrieves the number of eigenvalues.
  integer function NevalSpectralResolution( SpecRes ) &
       result( Neval )
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    Neval = SpecRes.NEigenvalues
  end function NevalSpectralResolution

  !> Set to + the sign of the first entry in each eigenvector
  subroutine SyncFirstSign( SpecRes ) 
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    integer :: iEval
    if(.not.allocated(SpecRes.Eigenvectors))return
    do iEval = 1, SpecRes.Neval()
       if(SpecRes.Eigenvectors(1,iEval)<0.d0)then
          SpecRes.Eigenvectors(:,iEval)=-SpecRes.Eigenvectors(:,iEval)
       endif
    enddo
  end subroutine SyncFirstSign

  !> Eliminates the null eigenspace
  subroutine SpectralResolutionPurgeNull( SpecRes, Threshold )
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    DoublePrecision, optional     , intent(in)    :: Threshold
    DoublePrecision, parameter :: DEFAULT_THRESHOLD = 1.d-16
    integer, allocatable :: ivec(:)
    integer :: iEval
    integer :: NValidEval
    DoublePrecision :: ActualThreshold
    ActualThreshold=DEFAULT_THRESHOLD
    if(present(Threshold))ActualThreshold=Threshold
    allocate(ivec(SpecRes.NEigenvalues))
    ivec=0
    NValidEval=0
    do iEval = 1, SpecRes.Neval()
       if( abs(SpecRes.EigenValues(iEval))<=Threshold )cycle
       NValidEval=NValidEval+1
       ivec(NValidEval)=iEval
    enddo
    SpecRes.NEigenvalues=NValidEval
    do iEval = 1, NValidEval
       SpecRes.EigenValues(   iEval) = SpecRes.EigenValues(   ivec(iEval))
       SpecRes.EigenVectors(:,iEval) = SpecRes.EigenVectors(:,ivec(iEval))
    enddo
    deallocate(ivec)
  end subroutine SpectralResolutionPurgeNull



  !> Retrieves whether a previously stored spectral resolution class information 
  !! in a unit, is consistent or not with a new one available.
  logical function SpectralResolutionIsConsistent( SpecRes, FileName )&
       result( IsConsistent )
    !
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    character(len=*)              , intent(in)    :: FileName
    !
    integer :: uid, iostat, Dim, Neval
    !
    IsConsistent = .FALSE.
    OPEN(NewUnit =  uid         , &
         File    =  FileName    , &
         Form    = "unformatted", &
         Status  = "old"        , &
         Action  = "read"       , &
         iostat  = iostat       )
    if( iostat /= 0 )return
    !
    read(uid,iostat=iostat) Dim, NEval
    if(iostat==0)then
       IsConsistent = ( Dim == SpecRes.Dim .and. NEval <= Dim )
    endif
    close(uid)
    !
  end function SpectralResolutionIsConsistent


  !> Writes in a unit the spectral resolution class information.
  subroutine WriteSpectralResolutionToUnit( SpecRes, uid )
    !
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    integer,              optional, intent(in) :: uid
    !
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Writable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    integer           :: OutUnit 
    if( present(uid) )then
       OutUnit=uid
    else
       OutUnit=DEFAULT_OUTPUT_UNIT
    endif

    INQUIRE(&
         UNIT  = OutUnit , &
         OPENED= Opened  , &
         WRITE = Writable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Output Unit is not open")
    if( trim(Writable) /= "YES" ) call Assert("Output Unit can't be written")
    !
    select case (trim(FORM))
    case("FORMATTED")
       call WriteSpectralResolutionToFormattedUnit(SpecRes,OutUnit)
    case("UNFORMATTED")
       call WriteSpectralResolutionToUnformattedUnit(SpecRes,OutUnit)
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
    !
  end subroutine WriteSpectralResolutionToUnit


  !> Writes in a unit the spectral resolution class information.
  subroutine WriteSpectralResolutionToUnformattedUnit( SpecRes, uid )
    !
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    integer,                        intent(in) :: uid
    integer :: iostat1, iostat2, iostat3
    integer :: i, j
    !
    write( unit=uid, iostat=iostat1 )     SpecRes.Dim, SpecRes.NEigenvalues
    write( unit=uid, iostat=iostat2 ) (   SpecRes.EigenValues (j)  , j=1, SpecRes.NEigenvalues )
    write( unit=uid, iostat=iostat3 ) ( ( SpecRes.EigenVectors(i,j), i=1, SpecRes.Dim ), j=1, SpecRes.NEigenvalues )
    !
    if ( iostat1 /=0 ) call ErrorMessage( 'Error trying to write the new upper index and the number of regular functions' )
    if ( iostat2 /=0 ) call ErrorMessage( 'Error trying to write the hamiltonian eigenvalues on file'  )
    if ( iostat3 /=0 ) call ErrorMessage( 'Error trying to write the hamiltonian eigenvectors on file' )
    !
  end subroutine WriteSpectralResolutionToUnformattedUnit


  !> Writes in a unit the spectral resolution class information.
  subroutine WriteSpectralResolutionToFormattedUnit( SpecRes, uid )
    !
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    integer,                        intent(in) :: uid
    integer :: iostat1, iostat2, iostat3
    integer :: i, j
    !
    write( uid, "(2(x,i8))", iostat=iostat1 )     SpecRes.Dim, SpecRes.NEigenvalues
    write( uid, "(*(x,e24.16))", iostat=iostat2 ) (   SpecRes.EigenValues (j)  , j=1, SpecRes.NEigenvalues )
    write( uid, "(a)" ) " "
    do j=1, SpecRes.NEigenvalues
       write( uid, "(*(x,e24.16))", iostat=iostat3 ) ( SpecRes.EigenVectors(i,j), i=1, SpecRes.Dim )
       if(iostat3/=0)exit
    enddo
    !
    if ( iostat1 /=0 ) call ErrorMessage( 'Error trying to write the new upper index and the number of regular functions' )
    if ( iostat2 /=0 ) call ErrorMessage( 'Error trying to write the hamiltonian eigenvalues on unit'  )
    if ( iostat3 /=0 ) call ErrorMessage( 'Error trying to write the hamiltonian eigenvectors on unit' )
    !
  end subroutine WriteSpectralResolutionToFormattedUnit


  !> Writes in a file the spectral resolution class information.
  subroutine WriteSpectralResolutionToFile( SpecRes, FileName )
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    character(len=*)              , intent(in) :: FileName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    open(newunit =  uid         , &
         file    =  FileName    , &
         form    = "unformatted", &
         status  = "unknown"    , &
         action  = "write"      , &
         iostat  =  iostat      , &
         iomsg   =  iomsg       )
    if(iostat/=0)call Assert(iomsg)
    !
    call SpecRes.Write( uid )
    !
    close( uid )
    !
  end subroutine WriteSpectralResolutionToFile


  !> Reads from a file the previously stored spectral resolution class information.
  subroutine ReadSpectralResolutionFromFile( SpecRes, FileName )
    Class(ClassSpectralResolution), intent(out) :: SpecRes
    character(len=*)              , intent(in)  :: FileName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    open(Newunit =  uid         , &
         File    =  FileName    , &
         Status  = "old"        , &
         Action  = "read"       , &
         Form    = "unformatted", &
         iostat  =  iostat      , &
         iomsg   =  iomsg       )
    if ( iostat /= 0 ) then
       call ErrorMessage(iomsg)
       return
    endif
    !
    call SpecRes.Read( uid )
    !
    close( uid )
    !
  end subroutine ReadSpectralResolutionFromFile


  !> Reads from a unit the previously stored spectral resolution class information.
  subroutine ReadSpectralResolutionFromUnit( SpecRes, uid )
    Class(ClassSpectralResolution), intent(out) :: SpecRes
    integer                       , intent(in)  :: uid
    !
    integer :: i, j, iostat, Dim, Neval
    character(len=IOMSG_LENGTH) :: iomsg
    !
    read( uid, iostat=iostat, iomsg=iomsg ) Dim, Neval
    if (iostat/=0) call ErrorMessage(iomsg )
    call SpecRes.Init( Dim, Neval )
    read( uid, iostat=iostat, iomsg=iomsg )  ( SpecRes.Eigenvalues(j), j=1, Neval )             
    if (iostat/=0) call ErrorMessage(iomsg )
    read( uid, iostat=iostat, iomsg=iomsg ) (( SpecRes.Eigenvectors(i,j), i=1, Dim ), j=1, Neval )
    if (iostat/=0) call ErrorMessage(iomsg )
    !
  end subroutine ReadSpectralResolutionFromUnit


  !> Writes in a unit the spectral resolution class eigenvalues.
  subroutine WriteEigenvaluesToUnit( SpecRes, uid, NColumns, NumberFormat )!@
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    integer                       , intent(in) :: uid
    integer, optional             , intent(in) :: NColumns
    character(len=*), optional    , intent(in) :: NumberFormat
    !
    character(len=*), parameter :: DEFAULT_FORMAT = "f34.17"
    integer, parameter :: DEFAULT_NCOLUMNS = 5
    integer :: NCol, iEval
    integer :: iostat
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=32) :: form
    character(len=:), allocatable :: ActualFormat
    !
    if(present(NumberFormat))then
       allocate(ActualFormat,source=trim(NumberFormat))
    else
       allocate(ActualFormat,source=DEFAULT_FORMAT)
    endif
    NCol=DEFAULT_NCOLUMNS
    if(present(NColumns))then
       if(NColumns<1)then
          call ErrorMessage("WriteEigenvalues: Wrong number of columns")
          return
       endif
       NCol=NColumns
    endif
    !
    INQUIRE( uid, form=form, iostat=iostat, iomsg=iomsg )
    if(iostat/=0)then
       call ErrorMessage("WriteEigenvalues: "//trim(iomsg))
       return
    endif
    if(trim(form)/="FORMATTED")then
       call ErrorMessage("WriteEigenvalues: wrong uid format")
       return
    endif
    !
    write(uid,*) SpecRes.Dim,SpecRes.NEigenvalues
    do iEval = 1, SpecRes.NEigenvalues
       !
       if(mod(iEval-1,NCol)==0)then
          if(iEval/=1) write(uid,*)
          write(uid,"(i5)",advance="no") iEval
       endif
       write(uid,"(x,"//ActualFormat//")",advance="no") SpecRes.EigenValues(iEval) 
       !
    enddo
    write(uid,*)
    !
    deallocate(ActualFormat)
  end subroutine WriteEigenvaluesToUnit


  !> Writes in a file the spectral resolution class eigenvalues.
  subroutine WriteEigenvaluesToFile( SpecRes, FileName, NColumns, NumberFormat )
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    character(len=*)              , intent(in) :: FileName
    integer, optional             , intent(in) :: NColumns
    character(len=*), optional    , intent(in) :: NumberFormat
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    open(newunit =  uid        , &
         file    =  FileName   , &
         form    = "formatted" , &
         status  = "unknown"   , &
         action  = "write"     , &
         iostat  =  iostat     , &
         iomsg   =  iomsg      )
    if(iostat/=0)call assert(iomsg)
    if(present(NColumns))then
       if(present(NumberFormat))then
          call SpecRes.WriteEigenvalues( uid, NColumns = NColumns, NumberFormat = NumberFormat )
       else
          call SpecRes.WriteEigenvalues( uid, NColumns = NColumns )
       endif
    else
       if(present(NumberFormat))then
          call SpecRes.WriteEigenvalues( uid, NumberFormat = NumberFormat )
       else
          call SpecRes.WriteEigenvalues( uid )
       endif
    endif
    close( uid, iostat=iostat, iomsg=iomsg )
    if(iostat/=0)call assert(iomsg)
    !
  end subroutine WriteEigenvaluesToFile


  !> Gets a square submatrix from the ClassMatrix's matrix defined by a 
  !! vector of indexes
  subroutine ClassMatrixGetIndexedSubMatrix( Matrix, SubMatrix, Indexes, NIndexes )!@
    !
    Class(ClassMatrix), intent(in)  :: Matrix  
    type(ClassMatrix) , intent(out) :: SubMatrix
    integer           , intent(in)  :: Indexes(:)
    integer           , intent(in)  :: NIndexes
    !
    integer         :: NL, NU
    integer         :: iRow, iCol, iSubRow, iSubCol
    integer         :: MaxAbsIndex
    DoublePrecision :: Element
    character(len=*), parameter :: HERE = "ClassMatrix::GetIndexedSubmatrix : "
    !
    if( ubound( Indexes, 1 ) < NIndexes )call Assert(HERE//"Invalid Index vector")
    MaxAbsIndex = min(Matrix.NR, Matrix.NC)
    if( Indexes(1) < 1 )call Assert(HERE//"index below lower bound")
    if( Indexes(1) > MaxAbsIndex )call Assert(HERE//"index above upper bound")
    do iRow = 2, NIndexes
       if( Indexes(iRow) <= Indexes(iRow-1) )call Assert(HERE//"Non-monotonous index vector")
       if( Indexes(iRow) > MaxAbsIndex) call Assert(HERE//"index above upper bound")
    enddo
    !
    NL = Matrix.LowerBandwidth()
    NU = Matrix.UpperBandwidth()
    !
    if( Matrix.IsFull() )then
       call SubMatrix.InitFull( NIndexes, NIndexes )
    elseif( Matrix.IsBanded() )then
       call SubMatrix.InitBanded( NIndexes, NIndexes, NL, NU )
    else
       call Assert("Unrecognized Matrix type")
    endif
    !
    do iSubCol = 1, NIndexes
       iCol = Indexes( iSubCol )
       do iSubRow = 1, NIndexes
          iRow = Indexes( iSubRow )
          Element = Matrix.Element( iRow, iCol )
          call SubMatrix.SetElement( iSubRow, iSubCol, Element )
       enddo
    enddo
    !
  end subroutine ClassMatrixGetIndexedSubMatrix


  !> Gets a rectangular submatrix from the ClassMatrix's matrix 
  !! defined by two independent vectors of row and column indexes
  subroutine GetAsymmetricallyIndexedSubMatrix( &
       Matrix, SubMatrix, IndexesBra , NIndexesBra, IndexesKet, NIndexesKet )!@
    Class(ClassMatrix), intent(in)  :: Matrix  
    type(ClassMatrix) , intent(out) :: SubMatrix
    integer           , intent(in)  :: IndexesBra(:)
    integer           , intent(in)  :: NIndexesBra
    integer           , intent(in)  :: IndexesKet(:)
    integer           , intent(in)  :: NIndexesKet
    !
    !.. Indexes are required to be in strictly increasing order and
    !   within the range of indexes of the matrix to be transformed. 
    !   Furthermore, in order to preserve matrix structure, or at 
    !   least keep track of it in a reasonably limited number of ways, 
    !   additional constraints may be required. 
    !..
    !
    integer         :: NL, NU, VerticalShift
    integer         :: iRow, iCol, iSubRow, iSubCol
    DoublePrecision :: Element
    character(len=*), parameter :: HERE = "ClassMatrix::GetAsymmetricallyIndexedSubmatrix : "
    !
    call CheckIndexes(NIndexesBra,IndexesBra,Matrix.NR)
    call CheckIndexes(NIndexesKet,IndexesKet,Matrix.NC)
    !
    if( Matrix.IsFull() )then
       call SubMatrix.InitFull( NIndexesBra, NIndexesKet )
    elseif( Matrix.IsBanded() )then
       call AdditionalCheckIndexesBanded(NIndexesBra,IndexesBra)
       call AdditionalCheckIndexesBanded(NIndexesKet,IndexesKet)
       VerticalShift = IndexesBra(1) - IndexesKet(1)
       NL = Matrix.LowerBandwidth() - VerticalShift
       NU = Matrix.UpperBandwidth() + VerticalShift
       call SubMatrix.InitBanded( NIndexesBra, NIndexesKet, NL, NU )
    else
       call Assert("Unrecognized Matrix type")
    endif
    !
    do iSubCol = 1, NIndexesKet
       iCol = IndexesKet( iSubCol )
       do iSubRow = 1, NIndexesBra
          iRow = IndexesBra( iSubRow )
          Element = Matrix.Element( iRow, iCol )
          call SubMatrix.SetElement( iSubRow, iSubCol, Element )
       enddo
    enddo
    !
  contains
    subroutine CheckIndexes(NIndexes,Indexes,MaxAbsIndex)
      integer, intent(in) :: NIndexes
      integer, intent(in) :: Indexes(:)
      integer, intent(in) :: MaxAbsIndex
      integer :: iRow
      if( ubound( Indexes, 1 ) < NIndexes )call Assert(HERE//"Invalid Index vector")
      if( Indexes(1) < 1 )call Assert(HERE//"index below lower bound")
      if( Indexes(1) > MaxAbsIndex )call Assert(HERE//"index above upper bound")
      do iRow = 2, NIndexes
         if( Indexes(iRow) <= Indexes(iRow-1) )call Assert(HERE//"Non-monotonous index vector")
         if( Indexes(iRow) > MaxAbsIndex) call Assert(HERE//"index above upper bound")
      enddo
    end subroutine CheckIndexes
    !.. For banded matrices, it is required that only extremal values are
    !   dropped, i.e., that the set of indexes is contiguous. In this way,
    !   the resulting matrix is still banded.
    subroutine AdditionalCheckIndexesBanded(NIndexes,Indexes)
      integer, intent(in) :: NIndexes
      integer, intent(in) :: Indexes(:)
      integer :: i
      do i=1,NIndexes-1
         if(Indexes(i+1)/=Indexes(i)+1)&
              call Assert(HERE//"Only contiguous set of subindexes are allowed for banded matrices")
      enddo
    end subroutine AdditionalCheckIndexesBanded
  end subroutine GetAsymmetricallyIndexedSubMatrix


  !> Gets a square submatrix from the ClassMatrix's matrix defined by a minimum
  !! and a maximum index which are assumed to be equal for both rows and columns.
  subroutine ClassMatrixGetSubMatrix_D( Matrix, MinSubIndex, MaxSubIndex, SubMatrix )!@
    !
    Class(ClassMatrix), intent(in) :: Matrix  
    integer,            intent(in)    :: MinSubIndex
    integer,            intent(in)    :: MaxSubIndex
    type(ClassMatrix),  intent(out)   :: SubMatrix
    !
    integer         :: SubDim, NL, NU
    integer         :: iRow, iCol, iSubRow, iSubCol, iSubRowMin, iSubRowMax
    DoublePrecision :: Element
    !
    if( MinSubIndex <= 0 )call Assert("Invalid MinSubIndex")
    if( MaxSubIndex > min( Matrix.NR, Matrix.NC ) ) &
         call Assert("Invalid MaxSubIndex")
    !
    SubDim = MaxSubIndex - MinSubIndex + 1
    NL = Matrix.LowerBandwidth()
    NU = Matrix.UpperBandwidth()
    !
    if( Matrix.IsFull() )then
       call SubMatrix.InitFull( SubDim, SubDim )
    elseif( Matrix.IsBanded() )then
       call SubMatrix.InitBanded( SubDim, SubDim, NL, NU )
    else
       call Assert("Unrecognized Matrix type")
    endif
    !
    do iSubCol = 1, SubDim
       iCol = ( MinSubIndex - 1 ) + iSubCol
       iSubRowMin = max(1,iSubCol-NU)
       iSubRowMax = min(SubDim,iSubCol+NL)
       do iSubRow = iSubRowMin, iSubRowMax 
          iRow = ( MinSubIndex - 1 ) + iSubRow
          Element = Matrix.Element( iRow, iCol )
          call SubMatrix.SetElement( iSubRow, iSubCol, Element )
       enddo
    enddo
    !
  end subroutine ClassMatrixGetSubMatrix_D

  !> Gets a squared submatrix from the ClassMatrix's matrix defined by a minimum 
  !! and a maximum index which are assumed to be equal for both rows and columns.
  subroutine ClassMatrixGetSubMatrix_C( Matrix, MinSubIndex, MaxSubIndex, SubMatrix )!@
    !
    Class(ClassMatrix), intent(in) :: Matrix  
    integer,            intent(in)    :: MinSubIndex
    integer,            intent(in)    :: MaxSubIndex
    type(ClassComplexMatrix),  intent(out)   :: SubMatrix
    !
    type(ClassMatrix) :: TempMatrix
    !
    call Matrix.GetSubmatrix(MinSubIndex,MaxSubIndex,TempMatrix)
    SubMatrix=TempMatrix
    call TempMatrix.Free()
    !
  end subroutine ClassMatrixGetSubMatrix_C


  !> Gets a rectangular matrix which is extracted from an original one, depending on the four indices specified, two that set the row interval and the other two that set the column interval.
  subroutine ClassMatrixGetSubMatrix_Rect( &
       Matrix     , &
       MinRowIndex, &
       MaxRowIndex, &
       MinColIndex, &
       MaxColIndex, &
       SubMatrix    )
    !
    Class(ClassMatrix), intent(in)    :: Matrix  
    integer,            intent(in)    :: MinRowIndex
    integer,            intent(in)    :: MaxRowIndex
    integer,            intent(in)    :: MinColIndex
    integer,            intent(in)    :: MaxColIndex
    type(ClassMatrix),  intent(out)   :: SubMatrix
    !
    integer :: i, j
    !
    call SubMatrix.InitFull( MaxRowIndex-MinRowIndex+1,MaxColIndex-MinColIndex+1 )
    do j = MinColIndex, MaxColIndex
       do i = MinRowIndex, MaxRowIndex
          SubMatrix.A(i-MinRowIndex+1,j-MinColIndex+1) = Matrix.A(i,j)
       end do
    end do
    !
  end subroutine ClassMatrixGetSubMatrix_Rect


  !> Muliplies the matrix contained in a matrix class by a real number.
  subroutine ClassMatrixTimesDouble( Matrix, Number )
    Class(ClassMatrix), intent(inout) :: Matrix 
    DoublePrecision,    intent(in)    :: Number
    Matrix.A = Number * Matrix.A
  end subroutine ClassMatrixTimesDouble


  !> Adds up two matrices belonging to two different matrix classes.
  subroutine AddMatrixToMatrix( Matrix, DeltaMatrix, Strength ) 
    Class(ClassMatrix), intent(inout) :: Matrix
    Class(ClassMatrix), intent(in)    :: DeltaMatrix
    logical, optional,  intent(in)    :: Strength
    !
    integer :: LBR, UBR, LBC, UBC, iRow,iCol
    !
    LBR=LBOUND(Matrix.A,1)
    UBR=UBOUND(Matrix.A,1)
    LBC=LBOUND(Matrix.A,2)
    UBC=UBOUND(Matrix.A,2)
    !
    if( DeltaMatrix.HasSameShapeAs(Matrix,Strength) )then
!!$       Matrix.A = Matrix.A + DeltaMatrix.A
       do iCol=LBC, UBC
          do iRow=LBR,UBR
             Matrix.A(iRow,iCol)=Matrix.A(iRow,iCol)+DeltaMatrix.A(iRow,iCol)
          enddo
       enddo
    else
       call Assert("Incompatible ClassMatrix shapes")
    endif
  end subroutine AddMatrixToMatrix


  !> Returns True if two matrices have the same number of rows, columns, lower subdiagonal and superdiagonals.
  logical function HasSameShapeAs( Mat1, Mat2, Strength ) result(Same)
    Class(ClassMatrix), intent(in) :: Mat1, Mat2
    logical, optional,  intent(in) :: Strength
    Same=.FALSE.
    if ( present(Strength) .and. Strength==.false. ) then
       if( Mat1.Pattern /= Mat2.Pattern ) return
       if(  Mat1.NR /= Mat2.NR .or. &
            Mat1.NC /= Mat2.NC ) return
       Same=.TRUE.
    else
       if( Mat1.Pattern /= Mat2.Pattern ) return
       if(  Mat1.NR /= Mat2.NR .or. &
            Mat1.NC /= Mat2.NC .or. &
            Mat1.NL /= Mat2.NL .or. &
            Mat1.NU /= Mat2.NU )return
       Same=.TRUE.
    end if
    return
  end function HasSameShapeAs


  !> Fetches the eigenvalues from ClassSpectralResolution.
  subroutine FetchEigenValues( SpecRes, Vector )
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    DoublePrecision, allocatable,   intent(out) :: Vector(:)
    allocate( Vector, source = SpecRes.EigenValues )
  end subroutine FetchEigenValues


  !> Fetches the eigenvectors from ClassSpectralResolution.
  subroutine FetchEigenVectors( SpecRes, Mat )
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    DoublePrecision, allocatable,   intent(out) :: Mat(:,:)
    allocate( Mat, source = SpecRes.EigenVectors )
  end subroutine FetchEigenVectors


  !> Fetches the requested single eigenvector from ClassSpectralResolution.
  subroutine FetchSingleEigenVector( SpecRes, Vec, n )
    Class(ClassSpectralResolution), intent(in)  :: SpecRes
    DoublePrecision, allocatable  , intent(out) :: Vec(:)
    integer                       , intent(in)  :: n
    integer :: Dim
    if(n>SpecRes.Neval())call ErrorMessage("Requested eigenvector doesn't exist")
    Dim = size(SpecRes.Eigenvectors,1)
    if(allocated(Vec))then
       if(size(Vec,1)/=Dim)then
          deallocate(Vec)
          allocate(Vec(Dim))
       endif
    else
       allocate(Vec(Dim))
    endif
    Vec = SpecRes.Eigenvectors(:,n)
  end subroutine FetchSingleEigenVector




  !> Fetches the requested single eigenvector from ClassSpectralResolution in ClassMatrix form.
  subroutine FetchSingleEigenVectorMat( SpecRes, Mat, n )
    Class(ClassSpectralResolution), intent(in)  :: SpecRes
    class(ClassMatrix)            , intent(out) :: Mat
    integer                       , intent(in)  :: n
    integer :: Dim
    real(kind(1d0)), allocatable :: Array(:,:)
    if(n>SpecRes.Neval())call ErrorMessage("Requested eigenvector doesn't exist")
    allocate( Array(SpecRes.Dim,1) )
    Array(:,1) = SpecRes.Eigenvectors(:,n)
    Mat = Array
    deallocate( Array )
  end subroutine FetchSingleEigenVectorMat


  !> Transforms the eigenvectors C of the spectral resolution
  !! on the basis defined by the columns U of NewBasis:
  !!\f[
  !! C'=U^{T}\cdot C,
  !!\f]
  subroutine TransformEigenvectors( SpecRes, NewBasis )
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    Class(ClassSpectralResolution), intent(in)    :: NewBasis
    DoublePrecision, allocatable :: Matrix(:,:)
    integer :: NR,NC,NK
    NR=NewBasis.NEigenvalues
    NC=SpecRes.NEigenvalues
    NK=min(SpecRes.Dim,NewBasis.Dim)
    allocate(Matrix(NR,NC))
    Matrix=0.d0
    call DGEMM( "T", "N", NR, NC, NK, 1.d0,&
         NewBasis.Eigenvectors, NewBasis.Dim, &
         SpecRes.Eigenvectors,  SpecRes.Dim, &
         0.d0, Matrix, NR )
    deallocate(SpecRes.Eigenvectors)
    allocate(SpecRes.Eigenvectors,source=Matrix)
    SpecRes.Dim=NR
    deallocate(Matrix)
  end subroutine TransformEigenvectors


  !> transform the eigenvectors C of the spectral resolution
  !! on the basis defined by the columns U of NewBasis, taking
  !! into account the metric S of the basis:
  !!\f[
  !! C'=U^{T}\cdot S\cdot C,
  !!\f]
  subroutine TransformEigenvectorsWithMetric( SpecRes, NewBasis, Metric )
    Class(ClassSpectralResolution), intent(inout) :: SpecRes
    Class(ClassSpectralResolution), intent(in)    :: NewBasis
    Class(ClassMatrix)            , intent(in)    :: Metric
    !
    DoublePrecision, allocatable :: Matrix(:,:)
    DoublePrecision, allocatable :: vec(:)
    integer :: iEval
    !
    !.. Pre-multiplies the eigenvectors by the metric, using
    !   the fact that the metric is symmetric
    select case( Metric.Pattern )
    case( MATRIX_PATTERN_FULL )
       allocate(Matrix(SpecRes.Dim,SpecRes.NEigenvalues))
       Matrix=0.d0
       call DSYMM("L","U",SpecRes.Dim,SpecRes.NEigenvalues,1.d0,&
            Metric.A,Metric.NR,&
            SpecRes.Eigenvectors,SpecRes.Dim,&
            0.d0,Matrix,SpecRes.Dim)
       deallocate(SpecRes.Eigenvectors)
       allocate(SpecRes.Eigenvectors,source=Matrix)
       deallocate(Matrix)
    case( MATRIX_PATTERN_BANDED )
       allocate(Vec(SpecRes.Dim))
       Vec=0.d0
       do iEval=1,SpecRes.NEigenvalues
          call DSBMV("U",SpecRes.Dim,Metric.NU,1.d0,&
               Metric.A,Metric.NL+Metric.NU+1,&
               SpecRes.Eigenvectors(1,iEval),1,&
               0.d0,Vec,1)
          SpecRes.Eigenvectors(:,iEval)=Vec
       enddo
       deallocate(Vec)
    case DEFAULT
       call Assert('Error: unrecognized matrix pattern')
    end select
    call SpecRes.Transform( NewBasis )
  end subroutine TransformEigenvectorsWithMetric


  !> Returns True if the matrix pattern is Full.
  logical function ClassComplexMatrixIsFull(Matrix) result(IsFull)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    IsFull = ( Matrix.Pattern == MATRIX_PATTERN_FULL )
  end function ClassComplexMatrixIsFull
  !
  !> Returns True if the matrix pattern is Banded.
  logical function ClassComplexMatrixIsBanded(Matrix) result(IsBanded)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    IsBanded = ( Matrix.Pattern == MATRIX_PATTERN_BANDED )
  end function ClassComplexMatrixIsBanded
  !
  !> Returns True if the matrix pattern is Diagonal.
  logical function ClassComplexMatrixIsDiagonal(Matrix) result(IsDiagonal)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    IsDiagonal = ( Matrix.Pattern == MATRIX_PATTERN_DIAGONAL )
  end function ClassComplexMatrixIsDiagonal

  !> Retrieves the number of  subdiagonals in the matrix. 
  integer function ClassComplexMatrixLowerBandwidth(Matrix) result(LowerBandwidth)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    LowerBandwidth = Matrix.NL
  end function ClassComplexMatrixLowerBandwidth

  !> Retrieves the number of superdiagonals in the matrix. 
  integer function ClassComplexMatrixUpperBandwidth(Matrix) result(UpperBandwidth)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    UpperBandwidth = Matrix.NU
  end function ClassComplexMatrixUpperBandwidth


  !> Writes the relevant information of the complex matrix class to a unit.
  subroutine ClassComplexMatrixWriteToUnit(Matrix,OutputUnit)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer, optional , intent(in) :: OutputUnit
    !
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=16) :: Writable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    integer           :: OutUnit 
    if(present(OutputUnit))then
       OutUnit=OutputUnit
    else
       OutUnit = DEFAULT_OUTPUT_UNIT
    endif
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
    !
    select case (trim(FORM))
    case("FORMATTED")
       call WriteClassComplexMatrixToFormattedUnit(Matrix,OutUnit)
    case("UNFORMATTED")
       call WriteClassComplexMatrixToUnformattedUnit(Matrix,OutUnit)
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
    !
  end subroutine ClassComplexMatrixWriteToUnit


  subroutine WriteClassComplexMatrixToFormattedUnit(Matrix,OutUnit)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer           , intent(in) :: OutUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=*), parameter :: FORMAT_ZMAT="(*"//COMPLEX_PRINT_FORMAT//")"
    character(len=*), parameter :: FORMAT_INTS="(*(x,i))"
    !
    !.. Write the Matrix attributes and dimensions
    write(OutUnit      , &
         FORMAT_INTS   , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix.Pattern, &
         Matrix.NR     , &
         Matrix.NC     , &
         Matrix.NU     , &
         Matrix.NL     , &
         Matrix.NRmin  , &
         Matrix.NRmax  , &
         Matrix.NCmin  , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    !.. Write Matrix
    do j=Matrix.NCmin,Matrix.NCmax
       write(OutUnit      , &
            FORMAT_ZMAT, &
            IOSTAT=iostat , &
            IOMSG =iomsg  ) &
            ( Matrix.A(i,j), i = Matrix.NRmin, Matrix.NRmax )
       if(iostat/=0)call Assert(iomsg)
    enddo
    !
  end subroutine WriteClassComplexMatrixToFormattedUnit
  !
  subroutine WriteClassComplexMatrixToUnformattedUnit(Matrix,OutUnit)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer           , intent(in) :: OutUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    !
    !.. Write the Matrix attributes and dimensions
    write(OutUnit      , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix.Pattern, &
         Matrix.NR     , &
         Matrix.NC     , &
         Matrix.NU     , &
         Matrix.NL     , &
         Matrix.NRmin  , &
         Matrix.NRmax  , &
         Matrix.NCmin  , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    !.. Write Matrix
    write(OutUnit      , &
         IOSTAT=iostat , &
         IOMSG =iomsg  ) &
         ((Matrix.A(i,j),&
         i=Matrix.NRmin,Matrix.NRmax),&
         j=Matrix.NCmin,Matrix.NCmax)
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine WriteClassComplexMatrixToUnformattedUnit


  !> Reads the relevant information of the complex matrix class from a unit.
  subroutine ClassComplexMatrixReadFromUnit(Matrix,InputUnit)
    Class(ClassComplexMatrix), intent(inout):: Matrix
    integer, optional , intent(in)   :: InputUnit
    !
    integer :: InUnit
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=FILE_NAME_LENGTH)     :: FileName
    character(len=16) :: Readable
    character(len=16) :: Form
    integer           :: iostat
    logical           :: Opened
    if(present(InputUnit))then
       InUnit = InputUnit
    else
       InUnit= DEFAULT_INPUT_UNIT
    end if
    INQUIRE(&
         UNIT  = InUnit  , &
         NAME  = FileName, &
         OPENED= Opened  , &
         READ  = Readable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Input Unit is closed")
    if( trim(Readable) /= "YES" ) call Assert("Input Unit can't be read")
    !
    select case (trim(FORM))
    case("FORMATTED")
       call ReadClassComplexMatrixFromFormattedUnit(Matrix,InUnit)
    case("UNFORMATTED")
       call ReadClassComplexMatrixFromUnformattedUnit(Matrix,InUnit)
    case DEFAULT
       call Assert("Invalid Input Unit Format")
    end select
    !
  end subroutine ClassComplexMatrixReadFromUnit
  !
  subroutine ReadClassComplexMatrixFromFormattedUnit(Matrix,InUnit)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: InUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg=" "
    character(len=*), parameter :: FORMAT_ZMAT="(*"//COMPLEX_PRINT_FORMAT//")"
    character(len=*), parameter :: FORMAT_INTS="(*(x,i))"
    !
    !.. Read the Matrix attributes and dimensions
    read(InUnit,FMT=FORMAT_INTS    , &
         IOSTAT=iostat   , &
         IOMSG=iomsg     ) &
         Matrix.Pattern  , &
         Matrix.NR       , &
         Matrix.NC       , &
         Matrix.NU       , &
         Matrix.NL       , &
         Matrix.NRmin    , &
         Matrix.NRmax    , &
         Matrix.NCmin    , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    call AllocateMatrix( Matrix.A  , &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    !
    !.. Read Matrix
    do j=Matrix.NCmin,Matrix.NCmax
       read(InUnit,FMT=FORMAT_ZMAT  , &
            IOSTAT=iostat , &
            IOMSG =iomsg  ) &
            (Matrix.A(i,j), &
            i=Matrix.NRmin,Matrix.NRmax)
       if(iostat/=0)call Assert(iomsg)
    enddo
    !
  end subroutine ReadClassComplexMatrixFromFormattedUnit
  !
  subroutine ReadClassComplexMatrixFromUnformattedUnit(Matrix,InUnit)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: InUnit
    integer :: iostat, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    !
    !.. Write the Matrix attributes and dimensions
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix.Pattern, &
         Matrix.NR     , &
         Matrix.NC     , &
         Matrix.NU     , &
         Matrix.NL     , &
         Matrix.NRmin  , &
         Matrix.NRmax  , &
         Matrix.NCmin  , &
         Matrix.NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    call AllocateMatrix( Matrix.A  , &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    !
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG =iomsg  ) &
         ((Matrix.A(i,j),&
         i=Matrix.NRmin,Matrix.NRmax),&
         j=Matrix.NCmin,Matrix.NCmax)
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine ReadClassComplexMatrixFromUnformattedUnit


  !> Multiplies the matrix by a complex number.
  subroutine ClassComplexMatrixMultiplyByComplex(Matrix,x)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    Complex(kind(1d0)), intent(in)    :: x
    if(.not.allocated(Matrix.A))then
       call Assert("Matrix not initialized")
    else
       Matrix.A=Matrix.A*x
    endif
  end subroutine ClassComplexMatrixMultiplyByComplex


  !> Multiplies the complex matrix by a real matrix
  subroutine ClassComplexMatrixMultiplyByClassMatrix( MatrixA, MatrixB, Side, MatrixBType )
    Class(ClassComplexMatrix), intent(inout) :: MatrixA
    Class(ClassMatrix),        intent(in)    :: MatrixB
    character(len=*),          intent(in)    :: Side
    character(len=*),          intent(in)    :: MatrixBType
    !
    type(ClassComplexMatrix) :: NewMatrixB
    !
    NewMatrixB = MatrixB
    !
    call ClassComplexMatrixMultiplyByClassComplexMatrix( MatrixA, NewMatrixB, Side, MatrixBType )
    !
  end subroutine ClassComplexMatrixMultiplyByClassMatrix


  !> Multiplies the complex matrix by a complex matrix
  subroutine ClassComplexMatrixMultiplyByClassComplexMatrix( MatrixA, MatrixB, Side, MatrixBType )
    Class(ClassComplexMatrix), target, intent(inout) :: MatrixA
    Class(ClassComplexMatrix), target, intent(in)    :: MatrixB
    character(len=*),                  intent(in)    :: Side
    character(len=*),                  intent(in)    :: MatrixBType
    !
    Class(ClassComplexMatrix), pointer :: PtrA, PtrB
    character         :: TypeA, TypeB
    Type(ClassComplexMatrix) :: MatrixC
    integer :: m, n, k, lda, ldb, ldc
    !
    if( ( MatrixBType .isnt. "N" ) .and. ( MatrixBType .isnt. "T" ) ) &
         call Assert( 'Improper matrix type for multiplication, it must be "N" ("n") or "T" ("t")' )
    !
    if( Side .is. "Right" )then
       PtrA  => MatrixA
       TypeA = "N"
       PtrB  => MatrixB
       TypeB =  MatrixBType
    elseif( Side .is. "Left" )then
       PtrA  => MatrixB
       TypeA =  MatrixBType
       PtrB  => MatrixA
       TypeB = "N"
    else
       call Assert( ' improper multiplication side, must be "right" or "left"' )
    endif

    lda = size( PtrA.A, 1 )
    if( TypeA .is. "N" )then
       m = size( PtrA.A, 1 )
       k = size( PtrA.A, 2 )
    else
       m = size( PtrA.A, 2 )
       k = size( PtrA.A, 1 )
    endif

    ldb = size( PtrB.A, 1 )
    if(TypeB .is. "N")then
       n = size( PtrB.A, 2 )
    else
       n = size( PtrB.A, 1 )
    endif

    ldc = m
    call MatrixC.InitFull( m, n )
    !
    call ZGEMM( &
         TypeA, TypeB, m, n, k, (1.d0,0.d0), &
         PtrA.A,    lda, &
         PtrB.A,    ldb, (0.d0,0.d0), &
         MatrixC.A, ldc )    
    !
    MatrixA = MatrixC
    !
    call MatrixC.Free()
    !
  end subroutine ClassComplexMatrixMultiplyByClassComplexMatrix


  !> Set all the elements of the matrix equal to an integer number.
  subroutine ClassComplexMatrixAssignInteger(Matrix,i)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    Integer           , intent(in)    :: i
    call ClassComplexMatrixAssignDouble(Matrix,dble(i))
  end subroutine ClassComplexMatrixAssignInteger

  !> Set all the elements of the matrix equal to real number.
  subroutine ClassComplexMatrixAssignDouble(Matrix,x)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    real(kind(1d0))   , intent(in)    :: x
    if(.not.allocated(Matrix.A))then
       call Assert("Matrix not initialized")
    else
       Matrix.A=x
    endif
  end subroutine ClassComplexMatrixAssignDouble


  !> Set all the elements of the matrix equal to a complex number.
  subroutine ClassComplexMatrixAssignComplex(Matrix,x)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    Complex(kind(1d0))   , intent(in)    :: x
    if(.not.allocated(Matrix.A))then
       call Assert("Matrix not initialized")
    else
       Matrix.A=x
    endif
  end subroutine ClassComplexMatrixAssignComplex


  !> Sets de complex matrix of the ClassComplexMatrix equal to an external one which is in full storage.
  subroutine ClassComplexMatrixAssignComplexMatrixFull(Matrix,FullMatrix)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    Complex(kind(1d0))   ,     intent(in)    :: FullMatrix(:,:)
    if(.not.allocated(Matrix.A))then
       call Matrix.InitFull( size(FullMatrix,1), size(FullMatrix,2) )
    endif
    Matrix.A = FullMatrix
  end subroutine 


  !> Copies the attributes from a ClassMatrix to a ClassComplexMatrix.
  subroutine ClassComplexMatrixCopyClassMatrix( MatrixOut, MatrixInp )
    Class(ClassComplexMatrix), intent(out) :: MatrixOut
    Class(ClassMatrix)       , intent(in)  :: MatrixInp
    integer :: LBR, UBR, LBC, UBC
    !
    call MatrixOut.Free()
    !
    MatrixOut.Pattern = MatrixInp.Pattern
    !
    MatrixOut.NR = MatrixInp.NR
    MatrixOut.NC = MatrixInp.NC
    MatrixOut.NL = MatrixInp.NL
    MatrixOut.NU = MatrixInp.NU
    !
    MatrixOut.NRMin = MatrixInp.NRMin
    MatrixOut.NRMax = MatrixInp.NRMax
    MatrixOut.NCMin = MatrixInp.NCMin
    MatrixOut.NCMax = MatrixInp.NCMax
    !
    LBR=LBOUND(MatrixInp.A,1)
    UBR=UBOUND(MatrixInp.A,1)
    LBC=LBOUND(MatrixInp.A,2)
    UBC=UBOUND(MatrixInp.A,2)
    allocate(MatrixOut.A(LBR:UBR,LBC:UBC))
    !
    MatrixOut.A = (1.d0,0.d0) * MatrixInp.A
    !
  end subroutine ClassComplexMatrixCopyClassMatrix


  subroutine ClassComplexMatrixAssignRealMatrix( MatrixOut, ArrayIn )
    !
    Class(ClassComplexMatrix), intent(out) :: MatrixOut
    real(kind(1d0))          , intent(in)  :: ArrayIn(:,:)
    !
    type(ClassMatrix) :: Matrix
    !
    Matrix = ArrayIn
    MatrixOut = Matrix
    call Matrix.Free()
    !
  end subroutine ClassComplexMatrixAssignRealMatrix



  !> Returns the real part of a complex matrix
  subroutine ClassComplexMatrixRealPartSub( ComplexMatrix, DoubleMatrix )
    Class(ClassComplexMatrix), intent(in)  :: ComplexMatrix
    Class(ClassMatrix)       , intent(out) :: DoubleMatrix
    integer :: LBR, UBR, LBC, UBC
    !
    call DoubleMatrix.Free()
    !
    DoubleMatrix.Pattern = ComplexMatrix.Pattern
    !
    DoubleMatrix.NR = ComplexMatrix.NR
    DoubleMatrix.NC = ComplexMatrix.NC
    DoubleMatrix.NL = ComplexMatrix.NL
    DoubleMatrix.NU = ComplexMatrix.NU
    !
    DoubleMatrix.NRMin = ComplexMatrix.NRMin
    DoubleMatrix.NRMax = ComplexMatrix.NRMax
    DoubleMatrix.NCMin = ComplexMatrix.NCMin
    DoubleMatrix.NCMax = ComplexMatrix.NCMax
    !
    LBR=LBOUND(ComplexMatrix.A,1)
    UBR=UBOUND(ComplexMatrix.A,1)
    LBC=LBOUND(ComplexMatrix.A,2)
    UBC=UBOUND(ComplexMatrix.A,2)
    allocate(DoubleMatrix.A(LBR:UBR,LBC:UBC))
    !
    DoubleMatrix.A = dble(ComplexMatrix.A)
    !
  end subroutine ClassComplexMatrixRealPartSub
  !
  function ClassComplexMatrixRealPartFun( ComplexMatrix ) result( DoubleMatrix )
    Class(ClassComplexMatrix), intent(in)  :: ComplexMatrix
    type(ClassMatrix)  :: DoubleMatrix
    call ComplexMatrix.RePart(DoubleMatrix)
  end function ClassComplexMatrixRealPartFun


  !> Returns the imaginary part of a complex matrix
  subroutine ClassComplexMatrixImaginaryPartSub( ComplexMatrix, DoubleMatrix )
    Class(ClassComplexMatrix), intent(in)  :: ComplexMatrix
    Class(ClassMatrix) :: DoubleMatrix
    integer :: LBR, UBR, LBC, UBC
    !
    call DoubleMatrix.Free()
    !
    DoubleMatrix.Pattern = ComplexMatrix.Pattern
    !
    DoubleMatrix.NR = ComplexMatrix.NR
    DoubleMatrix.NC = ComplexMatrix.NC
    DoubleMatrix.NL = ComplexMatrix.NL
    DoubleMatrix.NU = ComplexMatrix.NU
    !
    DoubleMatrix.NRMin = ComplexMatrix.NRMin
    DoubleMatrix.NRMax = ComplexMatrix.NRMax
    DoubleMatrix.NCMin = ComplexMatrix.NCMin
    DoubleMatrix.NCMax = ComplexMatrix.NCMax
    !
    LBR=LBOUND(ComplexMatrix.A,1)
    UBR=UBOUND(ComplexMatrix.A,1)
    LBC=LBOUND(ComplexMatrix.A,2)
    UBC=UBOUND(ComplexMatrix.A,2)
    allocate(DoubleMatrix.A(LBR:UBR,LBC:UBC))
    !
    DoubleMatrix.A = aimag(ComplexMatrix.A)
    !
  end subroutine ClassComplexMatrixImaginaryPartSub
  !
  function ClassComplexMatrixImaginaryPartFun( ComplexMatrix ) result( DoubleMatrix )
    Class(ClassComplexMatrix), intent(in)  :: ComplexMatrix
    type(ClassMatrix) :: DoubleMatrix
    call ComplexMatrix.ImPart(DoubleMatrix)
  end function ClassComplexMatrixImaginaryPartFun


  !> Copies the information from a ClassComplexMatrix to another ClassComplexMatrix.
  subroutine ClassComplexMatrixCopyToClassComplexMatrix(MatrixOut,MatrixInp)
    Class(ClassComplexMatrix), intent(inout) :: MatrixOut
    Class(ClassComplexMatrix), intent(in)    :: MatrixInp
    integer :: LBR, UBR, LBC, UBC, iRow, iCol
    !
    call MatrixOut.Free()
    !
    MatrixOut.Pattern = MatrixInp.Pattern
    !
    MatrixOut.NR = MatrixInp.NR
    MatrixOut.NC = MatrixInp.NC
    MatrixOut.NL = MatrixInp.NL
    MatrixOut.NU = MatrixInp.NU
    !
    MatrixOut.NRMin = MatrixInp.NRMin
    MatrixOut.NRMax = MatrixInp.NRMax
    MatrixOut.NCMin = MatrixInp.NCMin
    MatrixOut.NCMax = MatrixInp.NCMax
    !
    LBR=LBOUND(MatrixInp.A,1)
    UBR=UBOUND(MatrixInp.A,1)
    LBC=LBOUND(MatrixInp.A,2)
    UBC=UBOUND(MatrixInp.A,2)
    allocate(MatrixOut.A(LBR:UBR,LBC:UBC))
    do iCol=LBC, UBC
       do iRow=LBR,UBR
          MatrixOut.A(iRow,iCol)=MatrixInp.A(iRow,iCol)
       enddo
    enddo
!!$    MatrixOut.A=MatrixInp.A
    !
  end subroutine ClassComplexMatrixCopyToClassComplexMatrix


  !> Creates and allocates a ClassComplexMatrix object with Full pattern, if it is not defined, and free it before, if it is.
  subroutine ClassComplexMatrixInitFull( Matrix, NR, NC )
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR
    integer           , intent(in)    :: NC
    !
    call Matrix.Free()
    call SetComplexMatrixNominalSize( Matrix, NR, NC )
    call SetComplexMatrixPhysicalSize( Matrix, 1, NR, 1, NC )
    call AllocateMatrix( Matrix.A,   &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    Matrix.Pattern = MATRIX_PATTERN_FULL
    Matrix.A = 0.d0
    !
  end subroutine ClassComplexMatrixInitFull


  !> Create and allocate a ClassComplexMatrix object with Banded pattern, if it is not defined, and free it before, if it is not.
  subroutine ClassComplexMatrixInitBanded( Matrix, NR, NC, NL, NU )
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer                  , intent(in)    :: NR
    integer                  , intent(in)    :: NC
    integer                  , intent(in)    :: NL
    integer                  , intent(in)    :: NU
    !
    call Matrix.Free()
    call SetComplexMatrixNominalSize( Matrix, NR, NC )
    Matrix.NL=max(0,min(NL,NR-1))
    Matrix.NU=max(0,min(NU,NC-1))
    call SetComplexMatrixPhysicalSize( Matrix, -Matrix.NU, Matrix.NL, 1, Matrix.NC )
    call AllocateMatrix( Matrix.A,   &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    Matrix.Pattern = MATRIX_PATTERN_BANDED
    Matrix.A = 0.d0
    !
  end subroutine ClassComplexMatrixInitBanded


  !> Sets the rows' and columns' minimum and maximum indexes.
  subroutine SetComplexMatrixPhysicalSize( Matrix, NRMin, NRMax, NCMin, NCMax )
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NRMin, NRMax, NCMin, NCMax
    Matrix.NRmin = NRMin
    Matrix.NRmax = NRMax
    Matrix.NCmin = NCMin
    Matrix.NCmax = NCMax
  end subroutine SetComplexMatrixPhysicalSize


  !> Fetches the requested matrix element.
  Complex(kind(1d0)) function ClassComplexMatrixElement(Matrix,i,j) result(Element) 
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    Element=0.d0
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       Element = ClassComplexMatrixFullElement(Matrix,i,j)
    case( MATRIX_PATTERN_BANDED )
       Element = ClassComplexMatrixBandedElement(Matrix,i,j)
    case DEFAULT
    end select
  end function ClassComplexMatrixElement


  !> Sets the requested element.
  subroutine ClassComplexMatrixSetElement(Matrix,i,j,x)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    Complex(kind(1d0))   , intent(in)    :: x
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call ClassComplexMatrixFullSetElement(Matrix,i,j,x)
    case( MATRIX_PATTERN_BANDED )
       call ClassComplexMatrixBandedSetElement(Matrix,i,j,x)
    case DEFAULT
    end select
  end subroutine ClassComplexMatrixSetElement

  !> Sets the requested element when the matrix has Full pattern.
  subroutine ClassComplexMatrixFullSetElement(Matrix,i,j,x)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    Complex(kind(1d0))   , intent(in)    :: x
    call CheckComplexMatrixIndexBounds(Matrix,i,j)
    Matrix.A(i,j)=x
  end subroutine ClassComplexMatrixFullSetElement

  !> Sets the requested element when the matrix has Banded pattern.
  subroutine ClassComplexMatrixBandedSetElement(Matrix,i,j,x)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    Complex(kind(1d0))   , intent(in)    :: x
    integer :: iPhys
    call CheckComplexMatrixIndexBounds(Matrix,i,j)
    iPhys=i-j
    if( iPhys < -Matrix.NU .or. iPhys > Matrix.NL )&
         call Assert("Invalid banded matrix index")
    Matrix.A(iPhys,j)=x
  end subroutine ClassComplexMatrixBandedSetElement


  !> Checks if requested matrix element indexes are valid.
  subroutine CheckComplexMatrixIndexBounds(Matrix,i,j)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    if( i < 1 .or. i > Matrix.NR )call Assert("Row index is off bound")
    if( j < 1 .or. j > Matrix.NC )call Assert("Column index is off bound")
  end subroutine CheckComplexMatrixIndexBounds


  !> Fetches the requested matrix element when the matrix has Full pattern.
  Complex(kind(1d0)) function ClassComplexMatrixFullElement(Matrix,i,j) result(Element)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    integer  :: iPhys,jPhys
    call CheckComplexMatrixIndexBounds(Matrix,i,j)
    iPhys=Matrix.NRMin-1+i
    jPhys=Matrix.NCMin-1+j
    Element=Matrix.A(iPhys,jPhys)
  end function ClassComplexMatrixFullElement


  !> Fetches the requested matrix element when the matrix has Banded pattern.
  Complex(kind(1d0)) function ClassComplexMatrixBandedElement(Matrix,i,j) result(Element)
    Class(ClassComplexMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    integer :: iPhys
    call CheckComplexMatrixIndexBounds(Matrix,i,j)
    iPhys=i-j
    if( iPhys >= -Matrix.NU .and. iPhys <= Matrix.NL )then
       Element=Matrix.A(iPhys,j)
    else
       Element=0.d0
    end if
  end function ClassComplexMatrixBandedElement
 

  !> Sets the number of rows and columns of the complex matrix.
  subroutine SetComplexMatrixNominalSize(Matrix,NR,NC)
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR,NC
    if(NR<=0)call Assert("Invalid Number of Rows")
    if(NC<=0)call Assert("Invalid Number of Columns")
    Matrix.NR=NR
    Matrix.NC=NC
  end subroutine SetComplexMatrixNominalSize


  !> Frees allocatable space and sets to zero all the members of a ClassComplexMatrix object.
  subroutine ClassComplexMatrixFree( Matrix )
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    if(allocated(Matrix.A))deallocate(Matrix.A)
    Matrix.NR=0
    Matrix.NC=0
    Matrix.Pattern=MATRIX_PATTERN_UNDEFINED
    Matrix.NU=0
    Matrix.NL=0
  end subroutine ClassComplexMatrixFree


  !> Frees allocatable space and sets to zero all the members of a ClassComplexMatrix object.
  subroutine ClassComplexMatrixFinalize( Matrix )
    type(ClassComplexMatrix), intent(inout) :: Matrix
    call Matrix.Free()
  end subroutine ClassComplexMatrixFinalize


  !> Gets the number of rows of the complex matrix.
  integer function ClassComplexMatrixNRows( Matrix ) result( NRows )
    Class(ClassComplexMatrix), intent(in) :: Matrix
    NRows = Matrix.NR
  end function ClassComplexMatrixNRows
  !
  !> Gets the number of columns of the complex matrix.
  integer function ClassComplexMatrixNColumns( Matrix ) result( NColumns )
    Class(ClassComplexMatrix), intent(in) :: Matrix
    NColumns = Matrix.NC
  end function ClassComplexMatrixNColumns

  !> Performs the summation of the matrix elements absolute value.
  DoublePrecision function ClassComplexMatrixNorm1( Matrix ) result( Norm1 )
    Class(ClassComplexMatrix), intent(in) :: Matrix
    Norm1 = sum(abs(Matrix.A))
  end function ClassComplexMatrixNorm1

  !> Performs the squared root , of the matrix elements absolute squared value summation.
  DoublePrecision function ClassComplexMatrixNorm2( Matrix ) result( Norm2 )
    Class(ClassComplexMatrix), intent(in) :: Matrix
    Norm2 = sqrt(sum(abs(Matrix.A)**2))
  end function ClassComplexMatrixNorm2

  !> Performs the summation of the matrix elements absolute value.
  DoublePrecision function ClassMatrixNorm1( Matrix ) result( Norm1 )
    Class(ClassMatrix), intent(in) :: Matrix
    Norm1 = sum(abs(Matrix.A))
  end function ClassMatrixNorm1

  !> Performs the squared root , of the matrix elements absolute squared value summation.
  DoublePrecision function ClassMatrixNorm2( Matrix ) result( Norm2 ) 
    Class(ClassMatrix), intent(in) :: Matrix
    Norm2 = sqrt(sum(abs(Matrix.A)**2))
  end function ClassMatrixNorm2


  !> Solves the general eigenvalue problem:
  !!
  !!\f[
  !! \left(\mathbb{A}-\alpha\cdot\mathbb{S}\right)\cdot\mathbb{c}=0
  !!\f] 
  !! Assumes Metric (\f$\mathbb{S}\f$) is definite positive.
  !! But also eliminates the first N eigenvectors that do not fulfill the condition number specified.
  subroutine GeneralComplexEigenvalueSolverWithRegularization(  Matrix, Metric, SpecRes, ConditionNumber)
    !
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassComplexMatrix),             intent(in)  :: Metric
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    real(kind(1d0)),                      intent(in)  :: ConditionNumber
    !
    type(ClassComplexMatrix) :: NewMatrix
    type(ClassComplexSpectralResolution) :: MetricSpecRes
    integer :: NumZeroEigVal, n, m, i
    complex(kind(1d0)), allocatable :: MetricEigVal(:), MetricEigVec(:,:), FirstMatTransf(:,:), PristineMat(:,:), SecondMatTransf(:,:), EigVec(:,:), NewEigVec(:,:), ArrayMetricEigVec(:,:)
    complex(kind(1d0)), parameter :: Z0 = (0.d0,0.d0)
    complex(kind(1d0)), parameter :: Z1 = (1.d0,0.d0)
    ! 
    call Metric.Diagonalize( MetricSpecRes )
    !
    call MetricSpecRes.Fetch( MetricEigVal )
    call MetricSpecRes.Fetch( MetricEigVec )
    !
    call MetricSpecRes.Free()
    !
    n = size(MetricEigVal)
    !
    NumZeroEigVal = 0
    !
    ! Supposes eigenvalues in ascending order.
    do i = 1, n
       !
       if ( dble(MetricEigVal(i))/dble(MetricEigVal(n)) < ConditionNumber ) then
          NumZeroEigVal = NumZeroEigVal + 1
       end if
       !
    end do
    write(*,*) "NumZeroEigVal", NumZeroEigVal
    !
    m = n - NumZeroEigVal
    !
    MetricEigVec(:,1:NumZeroEigVal) = Z0
    !
    do i = NumZeroEigVal+1, n
       MetricEigVec(:,i) = MetricEigVec(:,i)/sqrt(dble(MetricEigVal(i)))
    enddo
    !
    allocate( ArrayMetricEigVec(n,m) )
    ArrayMetricEigVec = MetricEigVec(:,NumZeroEigVal+1:)
    !
    !
    allocate(FirstMatTransf(n,m))
    !
    FirstMatTransf = Z0
    !
    call Matrix.FetchMatrix( PristineMat )
    !
!!$    call zgemm( "N", "N", n, m, n, Z1, PristineMat, n, &
!!$         MetricEigVec(1,n-m+1), n, Z0, FirstMatTransf, n )
    call zgemm( "N", "N", n, m, n, Z1, PristineMat, n, &
         ArrayMetricEigVec, n, Z0, FirstMatTransf, n )
    !
    deallocate(PristineMat)
    !
    allocate(SecondMatTransf(m,m))
    !
!!$    call zgemm( "C", "N", m, m, n, Z1, MetricEigVec(1,n-m+1), n, &
!!$         FirstMatTransf, n, Z0, SecondMatTransf, m )
    call zgemm( "C", "N", m, m, n, Z1, ArrayMetricEigVec, n, &
         FirstMatTransf, n, Z0, SecondMatTransf, m )
    !
    deallocate(FirstMatTransf)
    !
    NewMatrix = SecondMatTransf
    call NewMatrix.Diagonalize( SpecRes )
    !
    deallocate(SecondMatTransf)
    !
    call SpecRes.Fetch( EigVec )
    !
    allocate( NewEigVec(n,m) )
    NewEigVec = Z0
    !
!!$    call zgemm( "N", "N", n, m, m, Z1, MetricEigVec(1,n-m+1), n, &
!!$         EigVec, m, Z0, NewEigVec, n )
    call zgemm( "N", "N", n, m, m, Z1, ArrayMetricEigVec, n, &
         EigVec, m, Z0, NewEigVec, n )
    !
    deallocate( SpecRes.RightEigenvectors )
    allocate( SpecRes.RightEigenvectors, source = NewEigVec )
    SpecRes.Dim = size(NewEigVec,1)
    SpecRes.NEigenvalues = size(NewEigVec,2)
    deallocate( MetricEigVec, ArrayMetricEigVec )
    deallocate( EigVec )
    !
  end subroutine GeneralComplexEigenvalueSolverWithRegularization



  !! Preconditions the metric and the hamiltonian without diagonalizing the Hamiltonian. Also returns the transformation T matrix.
  !! Assumes Metric (\f$\mathbb{S}\f$) is definite positive.
  !! But also eliminates the first N eigenvectors that do not fulfill the condition number specified.
  subroutine GeneralComplexMatrixRegularization(  Matrix, Metric, ConditionNumber, TMat )
    !
    Class(ClassComplexMatrix),            intent(inout)  :: Matrix
    type(ClassComplexMatrix),             intent(inout)  :: Metric
    real(kind(1d0)),                      intent(in)  :: ConditionNumber
    class(ClassComplexMatrix),            intent(out) :: TMat
    !
    type(ClassComplexMatrix) :: NewMatrix, THMat
    type(ClassComplexSpectralResolution) :: MetricSpecRes
    integer :: NumZeroEigVal, n, m, i
    complex(kind(1d0)), allocatable :: MetricEigVal(:), MetricEigVec(:,:), FirstMatTransf(:,:), PristineMat(:,:), SecondMatTransf(:,:), EigVec(:,:), NewEigVec(:,:), ArrayMetricEigVec(:,:)
    complex(kind(1d0)), parameter :: Z0 = (0.d0,0.d0)
    complex(kind(1d0)), parameter :: Z1 = (1.d0,0.d0)
    ! 
    call Metric.Diagonalize( MetricSpecRes )
    !
    call MetricSpecRes.Fetch( MetricEigVal )
    call MetricSpecRes.Fetch( MetricEigVec )
    !
    call MetricSpecRes.Free()
    !
    n = size(MetricEigVal)
    !
    NumZeroEigVal = 0
    !
    ! Supposes eigenvalues in ascending order.
    do i = 1, n
       !
       if ( dble(MetricEigVal(i))/dble(MetricEigVal(n)) < ConditionNumber ) then
          NumZeroEigVal = NumZeroEigVal + 1
       end if
       !
    end do
    write(*,*) "NumZeroEigVal", NumZeroEigVal
    !
    m = n - NumZeroEigVal
    !
    MetricEigVec(:,1:NumZeroEigVal) = Z0
    !
    do i = NumZeroEigVal+1, n
       MetricEigVec(:,i) = MetricEigVec(:,i)/sqrt(dble(MetricEigVal(i))) 
    enddo
    !
    allocate( ArrayMetricEigVec(n,m) )
    ArrayMetricEigVec = MetricEigVec(:,NumZeroEigVal+1:)
    !
    TMat = ArrayMetricEigVec
    !
    call TMat.TransposeConjugate( THMat )
    !
    call Metric.Multiply( TMat, "Right", "N" )
    call Metric.Multiply( THMat, "Left", "N" )
    !
    !
    allocate(FirstMatTransf(n,m))
    !
    FirstMatTransf = Z0
    !
    call Matrix.FetchMatrix( PristineMat )
    !
!!$    call zgemm( "N", "N", n, m, n, Z1, PristineMat, n, &
!!$         MetricEigVec(1,n-m+1), n, Z0, FirstMatTransf, n )
    call zgemm( "N", "N", n, m, n, Z1, PristineMat, n, &
         ArrayMetricEigVec, n, Z0, FirstMatTransf, n )
    !
    deallocate(PristineMat)
    !
    allocate(SecondMatTransf(m,m))
    !
!!$    call zgemm( "C", "N", m, m, n, Z1, MetricEigVec(1,n-m+1), n, &
!!$         FirstMatTransf, n, Z0, SecondMatTransf, m )
    call zgemm( "C", "N", m, m, n, Z1, ArrayMetricEigVec, n, &
         FirstMatTransf, n, Z0, SecondMatTransf, m )
    !
    deallocate(FirstMatTransf)
    !
    NewMatrix = SecondMatTransf
    Matrix = NewMatrix
    call NewMatrix.Free()
    !
    deallocate(SecondMatTransf)
    !
    deallocate( MetricEigVec, ArrayMetricEigVec )
    !
  end subroutine GeneralComplexMatrixRegularization



  subroutine ClassComplexMatrixTransposeConjugateOld( Mat, OutMat )
    !
    class(ClassComplexMatrix), intent(in)    :: Mat
    class(ClassComplexMatrix), intent(inout) :: OutMat
    !
    complex(kind(1d0)), allocatable :: Array(:,:), OutArray(:,:)
    integer :: i, j
    !
    call Mat.FetchMatrix( Array )
    allocate ( OutArray(size(Array,2),size(Array,1)) )
    !
    do j = 1, size(OutArray,2)
       do i = 1, size(OutArray,1)
          OutArray(i,j) = conjg(Array(j,i))
       end do
    end do
    !
    OutMat = OutArray
    deallocate( Array )
    deallocate( OutArray )
    !
  end subroutine ClassComplexMatrixTransposeConjugateOld



  subroutine ClassComplexMatrixTransposeConjugate( Mat )
    !
    class(ClassComplexMatrix), intent(inout) :: Mat
    !
    complex(kind(1d0)), allocatable :: Array(:,:), OutArray(:,:)
    integer :: i, j
    !
    call Mat.FetchMatrix( Array )
    allocate ( OutArray(size(Array,2),size(Array,1)) )
    !
    do j = 1, size(OutArray,2)
       do i = 1, size(OutArray,1)
          OutArray(i,j) = conjg(Array(j,i))
       end do
    end do
    !
    call Mat.Free()
    Mat = OutArray
    deallocate( Array )
    deallocate( OutArray )
    !
  end subroutine ClassComplexMatrixTransposeConjugate



  subroutine ClassMatrixTransposeOld( Mat, OutMat )
    !
    class(ClassMatrix), intent(inout)    :: Mat
    class(ClassMatrix), intent(inout) :: OutMat
    !
    real(kind(1d0)), allocatable :: Array(:,:), OutArray(:,:)
    integer :: i, j
    !
    call Mat.FetchMatrix( Array )
    allocate ( OutArray(size(Array,2),size(Array,1)) )
    !
    do j = 1, size(OutArray,2)
       do i = 1, size(OutArray,1)
          OutArray(i,j) = Array(j,i)
       end do
    end do
    !
    OutMat = OutArray
    deallocate( Array )
    deallocate( OutArray )
    !
  end subroutine ClassMatrixTransposeOld



  subroutine ClassMatrixTranspose( Mat )
    !
    class(ClassMatrix), intent(inout)    :: Mat
    !
    real(kind(1d0)), allocatable :: Array(:,:), OutArray(:,:)
    integer :: i, j
    !
    call Mat.FetchMatrix( Array )
    call Mat.Free()
    allocate ( OutArray(size(Array,2),size(Array,1)) )
    !
    do j = 1, size(OutArray,2)
       do i = 1, size(OutArray,1)
          OutArray(i,j) = Array(j,i)
       end do
    end do
    !
    Mat = OutArray
    deallocate( Array )
    deallocate( OutArray )
    !
  end subroutine ClassMatrixTranspose



  subroutine ClassComplexMatrixTranspose( Mat, OutMat )
    !
    class(ClassComplexMatrix), intent(inout) :: Mat
    class(ClassComplexMatrix), intent(inout) :: OutMat
    !
    complex(kind(1d0)), allocatable :: Array(:,:), OutArray(:,:)
    integer :: i, j
    !
    call Mat.FetchMatrix( Array )
    allocate ( OutArray(size(Array,2),size(Array,1)) )
    !
    do j = 1, size(OutArray,2)
       do i = 1, size(OutArray,1)
          OutArray(i,j) = Array(j,i)
       end do
    end do
    !
    OutMat = OutArray
    deallocate( Array )
    deallocate( OutArray )
    !
  end subroutine ClassComplexMatrixTranspose



  !> Solves the general eigenvalue problem:
  !!
  !!\f[
  !! \left(\mathbb{A}-\alpha\cdot\mathbb{S}\right)\cdot\mathbb{c}=0
  !!\f]
  !! Assumes Metric (\f$\mathbb{S}\f$) is definite positive.
  subroutine GeneralComplexEigenvalueSolver( Matrix, Metric, SpecRes )
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassComplexMatrix),             intent(in)  :: Metric
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call FullGeneralComplexEigenvalueSolver( Matrix, Metric, SpecRes )
    case( MATRIX_PATTERN_BANDED )
       call BandGeneralComplexEigenvalueSolver( Matrix, Metric, SpecRes )
    case DEFAULT 
       call Assert('Error: unrecognized matrix pattern')
    end select
    !
  end subroutine GeneralComplexEigenvalueSolver
  !
  !.. Assumes Metric is definite positive
  subroutine GeneralComplexEigenvalueSolverRealMetric( Matrix, Metric, SpecRes )
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassMatrix),                    intent(in)  :: Metric
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    type(ClassComplexMatrix) :: ComplexMetric
    ComplexMetric = Metric
    call Matrix.Diagonalize(ComplexMetric,SpecRes)
  end subroutine GeneralComplexEigenvalueSolverRealMetric


  !.. Assumes the matrix is symmetric, 
  !   so that one needs only right eigenvectors
  subroutine FullGeneralComplexEigenvalueSolver ( Matrix, Metric, SpecRes )
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassComplexMatrix),             intent(in)  :: Metric
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes 
    !
    integer :: n, info, iEval, i, j
    integer(kind=8) :: lwork
    Complex(kind(1d0)), allocatable :: work(:)
    DoublePrecision   , allocatable :: rwork(:)
    Complex(kind(1d0)), allocatable :: A(:,:), B(:,:), ScalingDenominators(:)
    !
    !.. FROM LAPACK MANUAL:
    !   ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices 
    !   (A,B), the generalized eigenvalues, and optionally, the left and/or 
    !   right generalized eigenvectors. A generalized eigenvalue for a pair 
    !   of matrices (A,B) is a scalar lambda or a ratio alpha/beta = lambda, 
    !   such that A - lambda*B is singular. It is usually represented as the 
    !   pair (alpha,beta), as there is a reasonable interpretation for beta=0, 
    !   and even for both being zero. The right generalized eigenvector v(j) 
    !   corresponding to the generalized eigenvalue lambda(j) of (A,B) satisfies
    !   A * v(j) = lambda(j) * B * v(j).
    !   The left generalized eigenvector u(j) corresponding to the generalized 
    !   eigenvalues lambda(j) of (A,B) satisfies 
    !               u(j)**H * A = lambda(j) * u(j)**H * B
    !   where u(j)**H is the conjugate-transpose of u(j).
    !
    n = Matrix.NRows()
    call SpecRes.Init( n )
    !
    !.. Size Query
    lwork=-1 
    allocate(work(1))
    allocate(rwork(8*n))
    allocate(A,source=Matrix.A)
    allocate(B,source=Metric.A)
    allocate(ScalingDenominators(n))
    call ZGGEV("N","V",n,A,n,B,n,&
         SpecRes.Eigenvalues,ScalingDenominators,&
         SpecRes.LeftEigenvectors ,1,&
         SpecRes.RightEigenvectors,n,&
         work,lwork,rwork,info)
    if(info/=0)then
       call ErrorMessage("ZGGEV: Dimension query failed")
       lwork=2*n
    endif
    lwork=max(int(dble(work(1)))+1,2*n+1)
    deallocate(work)
    !
    !.. Prepare work space
    allocate(work(lwork))
    work=(0.d0,0.d0)
    rwork=0.d0
    ScalingDenominators=(0.d0,0.d0)
    !
    !.. Diagonalization
    !
    call ZGGEV("N","V",n,A,n,B,n,&
         SpecRes.Eigenvalues,ScalingDenominators,&
         SpecRes.LeftEigenvectors ,1,&
         SpecRes.RightEigenvectors,n,&
         work,lwork,rwork,info)
    !
    if(info<0)call ErrorMessage(&
         "ZGGEV: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0.and.info<=n)call ErrorMessage(&
         "ZGGEV: The QR algorithm failed to compute all the eigenvalues,"//&
         "and no eigenvectors have been computed; the eigenvalues "//&
         AlphabeticNumber(1)//":"//AlphabeticNumber(info)//" did not converge.")
    if(info==n+1)call ErrorMessage("ZGGEV: Other than QZ iteration failed in zhgeqz")
    if(info==n+2)call ErrorMessage("ZGGEV: Error return from ztgevc")
    if(info> n+2)call ErrorMessage("ZGGEV: unspecified LAPACK failure")
    !
    if(info>=0)then
       !
       do iEval=info+1,n
!!$          write(*,*) 'iEval',iEval
!!$          write(*,*) '   ', SpecRes.Eigenvalues(iEval)
!!$          write(*,*) '   ', ScalingDenominators(iEval)
          SpecRes.Eigenvalues(iEval-info)=SpecRes.Eigenvalues(iEval)/ScalingDenominators(iEval)
          SpecRes.RightEigenvectors(:,iEval-info)=SpecRes.RightEigenvectors(:,iEval)
       enddo
       SpecRes.NEigenvalues=n-info
       !
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
       !.. Sort eigenvalues and eigenvectors in ascending order of the real part
    endif
    !
    deallocate(work,rwork,A,B,ScalingDenominators)
    !
    !
  end subroutine FullGeneralComplexEigenvalueSolver
  !
  !.. There's no banded solver for general complex matrices.
  subroutine BandGeneralComplexEigenvalueSolver ( Matrix, Metric, SpecRes )
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassComplexMatrix),             intent(in)  :: Metric
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    type(ClassComplexMatrix) :: FullMatrix, FullMetric
    !
    call Matrix.ConvertToFull(FullMatrix)
    call Metric.ConvertToFull(FullMetric)
    call FullGeneralComplexEigenvalueSolver( FullMatrix, FullMetric, SpecRes )
    !
  end subroutine BandGeneralComplexEigenvalueSolver


  !> Converts the complex matrix in Banded representation to Full representation.
  subroutine ClassComplexMatrixConvertToFull( MatrixIn, MatrixOut )
    Class(ClassComplexMatrix),            intent(in)  :: MatrixIn
    type(ClassComplexMatrix),             intent(out) :: MatrixOut
    integer :: iRow, iCol, UBW, LBW, n
    complex(kind(1d0)) :: Element
    call MatrixOut.InitFull(MatrixIn.NRows(),MatrixIn.NColumns())
    N=MatrixIn.NColumns()
    UBW=MatrixIn.UpperBandWidth()
    LBW=MatrixIn.LowerBandWidth()
    do iCol=1,N
       do iRow=max(1,iCol-UBW),min(N,iCol+UBW)
          Element=MatrixIn.Element(iRow,iCol)
          call MatrixOut.SetElement(iRow,iCol,Element)
       enddo
    enddo
  end subroutine ClassComplexMatrixConvertToFull


  !> Solves the eigenvalues eigenvectors problem for a ClassComplexMatrix's matrix and stores the results in ClassComplexSpectralResolution.
  !! Assumes the matrix is symmetric? (maybe it was intended to write hermitic), so that one needs only right eigenvectors.
  subroutine ComplexEigenvalueSolver( Matrix, SpecRes )
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes 
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       if ( Matrix.IsHermitian(COMPUTATION_THRESHOLD) ) then
          call FullHermitianEigenvalueSolver( Matrix, SpecRes )
       else
          call FullComplexEigenvalueSolver( Matrix, SpecRes )
       end if
!!$          call FullComplexEigenvalueSolver( Matrix, SpecRes )
    case( MATRIX_PATTERN_BANDED )
       call BandComplexEigenvalueSolver( Matrix, SpecRes )
    case DEFAULT 
       call Assert('Error: non-proper matrix pattern')
    end select
    !
  end subroutine ComplexEigenvalueSolver


  !
  !> Solves the eigenvalues eigenvectors problem for a ClassComplexMatrix's matrix in Full representation and stores the results in ClassComplexSpectralResolution.
  subroutine FullComplexEigenvalueSolver( Matrix, SpecRes ) 
    Class(ClassComplexMatrix),            intent(in)  :: Matrix
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    !
    integer :: n, info, i, j
    integer(kind=8) :: lwork
    Complex(kind(1d0)), allocatable :: work(:)
    DoublePrecision   , allocatable :: rwork(:)
    Complex(kind(1d0)), allocatable :: A(:,:)
    complex(kind(1d0)), parameter :: Z1 = dcmplx(1.d0,0.d0)
    complex(kind(1d0)), parameter :: Zi = dcmplx(0.d0,1.d0)
    complex(kind(1d0)), parameter :: Z0 = dcmplx(0.d0,0.d0)
    complex(kind(1d0)) :: zw1
    complex(kind(1d0)), allocatable :: CR(:,:), CL(:,:), Prod(:,:)
    type(ClassComplexMatrix) :: Aux,Aux2
    !
    n = Matrix.NRows()
    call SpecRes.Init( n )
    !
    !.. Size Query
    lwork=-1 
    allocate(work(1))
    allocate(rwork(2*n))
    allocate(A,source=Matrix.A)
    call ZGEEV("N","V",n,A,n,&
         SpecRes.Eigenvalues,&
         SpecRes.LeftEigenvectors ,1,&
         SpecRes.RightEigenvectors,n,&
         work,lwork,rwork,info)
    if(info/=0)then
       call ErrorMessage("ZGEEV: Dimension query failed")
       lwork=2*n
    endif
    lwork=max(int(dble(work(1)))+1,2*n+1)
    deallocate(work)
    !
    !.. Prepare work space
    allocate(work(lwork))
    work=(0.d0,0.d0)
    rwork=0.d0
    !
    !.. Diagonalization
    call ZGEEV("N","V",n,A,n,&
         SpecRes.Eigenvalues,&
         SpecRes.LeftEigenvectors ,1,&
         SpecRes.RightEigenvectors,n,&
         work,lwork,rwork,info)
    !
    if(info<0)call ErrorMessage(&
         "ZGEEV: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0)call ErrorMessage(&
         "ZGEEV: The QR algorithm failed to compute all the eigenvalues,"//&
         "and no eigenvectors have been computed; the eigenvalues "//&
         AlphabeticNumber(1)//":"//AlphabeticNumber(info)//" did not converge.")
    !
    if(info==0)then
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
    endif
    !
    deallocate(work,rwork,A)
    !
    if ( Matrix.IsSymmetric(COMPUTATION_THRESHOLD) ) then
       !..Rescale the right eigenvectors CR to have the proper normalization:
       ! CL^{\dagga} * CR = 1. 
       allocate( CR, source = SpecRes.RightEigenvectors )
       call ReScaleRightEigVectors( CR )
       !
       deallocate( SpecRes.RightEigenvectors )
       allocate( SpecRes.RightEigenvectors, source = CR )
       deallocate( CR )
    end if
       !
  end subroutine FullComplexEigenvalueSolver
  


  !
  !> Solves the eigenvalues eigenvectors problem for a ClassComplexMatrix's matrix in Banded representation and stores the results in ClassComplexSpectralResolution.
  subroutine BandComplexEigenvalueSolver( Matrix, SpecRes )
    Class(ClassComplexMatrix),            intent(in)    :: Matrix
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    type(ClassComplexMatrix) :: FullMatrix
    !
    call Matrix.ConvertToFull(FullMatrix)
    if ( Matrix.IsHermitian(COMPUTATION_THRESHOLD) ) then
       call FullHermitianEigenvalueSolver( FullMatrix, SpecRes )
    else
       call FullComplexEigenvalueSolver( FullMatrix, SpecRes )
    end if
    !
  end subroutine BandComplexEigenvalueSolver



  ! Diagonalizes a hermitian complex matrix in full format.
  subroutine FullHermitianEigenvalueSolver( Matrix, SpecRes )
    !
    Class(ClassComplexMatrix),            intent(in)    :: Matrix
    type(ClassComplexSpectralResolution), intent(out) :: SpecRes
    !
    integer :: n, info
    integer(kind=8) :: lwork
    complex(kind(1d0)), allocatable :: work(:),wcplx(:), Array(:,:)
    real(kind(1d0)), allocatable :: rwork(:), w(:)
    complex(kind(1d0)), parameter :: Z1 = (1.d0,0.d0)
    !
!!$    write(*,*) "Diagonalizing hermitian matrix..."
    call Matrix.FetchMatrix( Array )
    ! 
    allocate( work(1) )
    n = size(Array,1)
    allocate( rwork(3*n-2), w(n) )
    call ZHEEV( 'V', 'U', &
         n, Array, n, &
         w, work, -1, rwork, info ) 
    lwork = max(int(work(1)),2*n-1)
    deallocate( work, rwork, w )
    allocate( work(lwork), rwork(3*n-2), w(n) )
    !
    call ZHEEV( 'V', 'U', &
         n, Array, n, &
         w, work, lwork, rwork, info ) 
    !
    allocate( wcplx(n) )
    wcplx = w * Z1
    !
    call SpecRes.SetEigenValues( wcplx )
    !
    if ( allocated(SpecRes.RightEigenvectors) ) then
       deallocate( SpecRes.RightEigenvectors )
    end if
    allocate( SpecRes.RightEigenvectors, source = Array )
    !
    if(info==0)then
       !
       call SpecRes.SyncFirstSign()
       call SpecRes.Sort()
       !
    else
       !
       call Assert( "Some error has ocurred diagonalizing the hermitian matrix, info = "//AlphabeticNumber(info) )
       !
    endif
    !
  end subroutine FullHermitianEigenvalueSolver


  !> The complex spectral resolution inherits its attributes from a non-complex spectral resolution.
  subroutine DoubleSpectralResolutionToComplexSpectralResolution( SpecRes_C, SpecRes_D )
    Class(ClassComplexSpectralResolution), intent(out) :: SpecRes_C
    Class(ClassSpectralResolution)       , intent(in)  :: SpecRes_D
    integer :: Dim, Neval
    call SpecRes_C.Free()
    Dim  =SpecRes_D.Size()
    Neval=SpecRes_D.Neval()
    call SpecRes_C.Init(Dim,Neval)
    SpecRes_C.Eigenvalues(1:Neval)=(1.d0,0.d0)*SpecRes_D.Eigenvalues(1:Neval)
    SpecRes_C.RightEigenvectors(1:Dim,1:Neval)=SpecRes_D.Eigenvectors(1:Dim,1:Neval)
  end subroutine DoubleSpectralResolutionToComplexSpectralResolution


  !> Initializes the complex spectral resolution class for the case when the matrix of eigenvectors is squared and the number of rows is equal to the number of eigenvalues.
  subroutine InitComplexSpectralResolutionFull( SpecRes, n )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    integer,                        intent(in)    :: n
    !
    SpecRes.NEigenvalues = n
    SpecRes.Dim = n
    !
    if(allocated(SpecRes.EigenValues))then
       if(size(SpecRes.Eigenvalues,1)/=n)then
          deallocate( SpecRes.EigenValues )
          allocate( SpecRes.EigenValues( n ))
       endif
    else
       allocate( SpecRes.EigenValues( n ))
    endif
    SpecRes.EigenValues = 0.d0
    !
    if(allocated(SpecRes.RightEigenvectors))then
       if(  size(SpecRes.RightEigenvectors,1) /= n  .or. &
            size(SpecRes.RightEigenvectors,2) /= n  )then
          deallocate(SpecRes.RightEigenvectors)
          allocate(SpecRes.RightEigenvectors( n, n ) )
       endif
    else
       allocate(SpecRes.RightEigenvectors( n, n ) )
    endif
    SpecRes.RightEigenvectors = 0.d0
    !
    if(allocated(SpecRes.LeftEigenvectors))deallocate(SpecRes.LeftEigenvectors)
    allocate(SpecRes.LeftEigenvectors(1,n))
    SpecRes.LeftEigenvectors=(0.d0,0.d0)
    !
  end subroutine InitComplexSpectralResolutionFull


  !> Initializes the complex spectral resolution class for the case when the matrix of eigenvectors is rectangular and the number of rows lower than the number of eigenvalues.
  subroutine InitComplexSpectralResolutionReduced( SpecRes, Dim, NEigenvalues )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    integer,                        intent(in)    :: Dim, NEigenvalues
    !
    SpecRes.Dim = Dim
    SpecRes.NEigenvalues = NEigenvalues
    !
    if(allocated(SpecRes.EigenValues))then
       if(size(SpecRes.Eigenvalues,1)/=NEigenvalues)then
          deallocate( SpecRes.EigenValues )
          allocate( SpecRes.EigenValues( NEigenvalues ))
       endif
    else
       allocate( SpecRes.EigenValues( NEigenvalues ))
    endif
    SpecRes.EigenValues = 0.d0
    !
    if(allocated(SpecRes.RightEigenvectors))then
       if(  size(SpecRes.RightEigenvectors,1) /= Dim  .or. &
            size(SpecRes.RightEigenvectors,2) /= NEigenvalues  )then
          deallocate(SpecRes.RightEigenvectors)
          allocate(SpecRes.RightEigenvectors( Dim, NEigenvalues ) )
       endif
    else
       allocate(SpecRes.RightEigenvectors( Dim, NEigenvalues ) )
    endif
    SpecRes.RightEigenvectors = 0.d0
    !
    if(allocated(SpecRes.LeftEigenvectors))deallocate(SpecRes.LeftEigenvectors)
    allocate(SpecRes.LeftEigenvectors(1,1))
    SpecRes.LeftEigenvectors=(0.d0,0.d0)
    !
  end subroutine InitComplexSpectralResolutionReduced

  !> Frees the ClassComplexSpectralResolution attributes.
  subroutine ClassComplexSpectralResolutionFree( SpecRes )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    SpecRes.NEigenvalues=0
    SpecRes.Dim=0
    if(allocated(SpecRes.Eigenvalues))deallocate(SpecRes.Eigenvalues)
    if(allocated(SpecRes.LeftEigenvectors))deallocate(SpecRes.LeftEigenvectors)
    if(allocated(SpecRes.RightEigenvectors))deallocate(SpecRes.RightEigenvectors)
  end subroutine ClassComplexSpectralResolutionFree

  !> Calls the Free subroutine.
  subroutine ClassComplexSpectralResolutionFinal( SpecRes )
    Type(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    call SpecRes.Free()
  end subroutine ClassComplexSpectralResolutionFinal


  !> Retrieves the number of eigenvalues.
  integer function NevalComplexSpectralResolution( SpecRes ) &
       result( Neval )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    Neval = SpecRes.NEigenvalues
  end function NevalComplexSpectralResolution

  !> Retrieves the number of rows of the eigenvectors matrices.
  integer function SizeComplexSpectralResolution( SpecRes ) &
       result( Size )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    Size = SpecRes.Dim
  end function SizeComplexSpectralResolution

  !> Retrieves the number of rows of the eigenvectors matrix.
  integer function SizeSpectralResolution( SpecRes ) result( Size )
    Class(ClassSpectralResolution), intent(in) :: SpecRes
    Size = SpecRes.Dim
  end function SizeSpectralResolution

  !> Set to + the sign of the first entry in each eigenvector
  subroutine ComplexSpectralResolutionSyncFirstSign( SpecRes ) 
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    integer :: iEval
    complex(kind(1d0)) :: UnitaryFactor
    if(.not.allocated(SpecRes.RightEigenvectors))return
    do iEval = 1, SpecRes.Neval()
       if(abs(SpecRes.RightEigenvectors(1,iEval))>0.d0)then
          UnitaryFactor=conjg(SpecRes.RightEigenvectors(1,iEval))/abs(SpecRes.RightEigenvectors(1,iEval))
          SpecRes.RightEigenvectors(:,iEval)=UnitaryFactor*SpecRes.RightEigenvectors(:,iEval)
       endif
    enddo
  end subroutine ComplexSpectralResolutionSyncFirstSign

  !> Eliminates the null eigenspace
  subroutine ComplexSpectralResolutionPurgeNull( SpecRes, Threshold )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Real(kind(1d0)), optional            , intent(in)    :: Threshold
    Complex(kind(1d0)), parameter :: DEFAULT_THRESHOLD = 1.d-16
    integer, allocatable :: ivec(:)
    integer :: iEval
    integer :: NValidEval
    Complex(kind(1d0)) :: ActualThreshold
    ActualThreshold=DEFAULT_THRESHOLD
    if(present(Threshold))ActualThreshold=Threshold
    allocate(ivec(SpecRes.NEigenvalues))
    ivec=0
    NValidEval=0
    do iEval = 1, SpecRes.Neval()
       if( abs(SpecRes.EigenValues(iEval))<=Threshold )cycle
       NValidEval=NValidEval+1
       ivec(NValidEval)=iEval
    enddo
    SpecRes.NEigenvalues=NValidEval
    do iEval = 1, NValidEval
       SpecRes.EigenValues(iEval) = SpecRes.EigenValues(ivec(iEval))
       SpecRes.RightEigenvectors(:,iEval) = SpecRes.RightEigenvectors(:,ivec(iEval))
    enddo
    deallocate(ivec)
  end subroutine ComplexSpectralResolutionPurgeNull


  !> Retrieves whether a previously stored complex spectral resolution class information in a unit, is consistent or not with a new one available.
  logical function ComplexSpectralResolutionIsConsistent( SpecRes, FileName )&
       result( IsConsistent )
    !
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    character(len=*)              , intent(in)    :: FileName
    !
    integer :: uid, iostat, Dim, Neval
    !
    IsConsistent = .FALSE.
    OPEN(NewUnit =  uid         , &
         File    =  FileName    , &
         Form    = "unformatted", &
         Status  = "old"        , &
         Action  = "read"       , &
         iostat  = iostat       )
    if( iostat /= 0 )return
    !
    read(uid,iostat=iostat) Dim, NEval
    if(iostat==0)then
       IsConsistent = ( Dim == SpecRes.Dim .and. NEval <= Dim )
    endif
    close(uid)
    !
  end function ComplexSpectralResolutionIsConsistent


  !> Writes in a unit, the complex spectral resolution class information.
  subroutine WriteComplexSpectralResolutionToUnit( SpecRes, uid )
    !
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    integer,                        intent(in) :: uid
    integer :: iostat1, iostat2, iostat3
    integer :: i, j
    !
    write( unit=uid, iostat=iostat1 )     SpecRes.Dim, SpecRes.NEigenvalues
    write( unit=uid, iostat=iostat2 ) (   SpecRes.EigenValues (j)  , j=1, SpecRes.NEigenvalues )
    write( unit=uid, iostat=iostat3 ) ( ( SpecRes.RightEigenvectors(i,j), i=1, SpecRes.Dim ), j=1, SpecRes.NEigenvalues )
    !
    if ( iostat1 /=0 ) call ErrorMessage( 'Error trying to write the new upper index and the number of regular functions' )
    if ( iostat2 /=0 ) call ErrorMessage( 'Error trying to write the hamiltonian eigenvalues on file'  )
    if ( iostat3 /=0 ) call ErrorMessage( 'Error trying to write the hamiltonian eigenvectors on file' )
    !
  end subroutine WriteComplexSpectralResolutionToUnit


  !> Writes in a file, the complex spectral resolution class information.
  subroutine WriteComplexSpectralResolutionToFile( SpecRes, FileName )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    character(len=*)              , intent(in) :: FileName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    open(newunit =  uid         , &
         file    =  FileName    , &
         form    = "unformatted", &
         status  = "unknown"    , &
         action  = "write"      , &
         iostat  =  iostat      , &
         iomsg   =  iomsg       )
    if(iostat/=0)call Assert(iomsg)
    !
    call SpecRes.Write( uid )
    !
    close( uid )
    !
  end subroutine WriteComplexSpectralResolutionToFile


  !> Reads from a file, the previously stored complex spectral resolution class information.
  subroutine ReadComplexSpectralResolutionFromFile( SpecRes, FileName )
    Class(ClassComplexSpectralResolution), intent(out) :: SpecRes
    character(len=*)              , intent(in)  :: FileName
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    open(Newunit =  uid         , &
         File    =  FileName    , &
         Status  = "old"        , &
         Action  = "read"       , &
         Form    = "unformatted", &
         iostat  =  iostat      , &
         iomsg   =  iomsg       )
    if ( iostat /= 0 ) then
       call ErrorMessage(iomsg)
       return
    endif
    !
    call SpecRes.Read( uid )
    !
    close( uid )
    !
  end subroutine ReadComplexSpectralResolutionFromFile

  !> Reads from a unit, the previously stored complex spectral resolution class information.
  subroutine ReadComplexSpectralResolutionFromUnit( SpecRes, uid )
    Class(ClassComplexSpectralResolution), intent(out) :: SpecRes
    integer                       , intent(in)  :: uid
    !
    integer :: i, j, iostat, Dim, Neval
    character(len=IOMSG_LENGTH) :: iomsg
    !
    read( uid, iostat=iostat, iomsg=iomsg ) Dim, Neval
    if (iostat/=0) call ErrorMessage(iomsg )
    call SpecRes.Init( Dim, Neval )
    read( uid, iostat=iostat, iomsg=iomsg )  ( SpecRes.Eigenvalues(j), j=1, Neval )             
    if (iostat/=0) call ErrorMessage(iomsg )
    read( uid, iostat=iostat, iomsg=iomsg ) (( SpecRes.RightEigenvectors(i,j), i=1, Dim ), j=1, Neval )
    if (iostat/=0) call ErrorMessage(iomsg )
    !
  end subroutine ReadComplexSpectralResolutionFromUnit


  !> Writes in a unit, the complex spectral resolution class eigenvalues.
  subroutine WriteComplexEigenvaluesToUnit( SpecRes, uid, NColumns, NumberFormat )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    integer                       , intent(in) :: uid
    integer, optional             , intent(in) :: NColumns
    character(len=*), optional    , intent(in) :: NumberFormat
    !
    character(len=*), parameter :: DEFAULT_FORMAT = "(2(x,f24.16))"
    integer, parameter :: DEFAULT_NCOLUMNS = 5
    integer :: NCol, iEval
    integer :: iostat
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=32) :: form
    character(len=:), allocatable :: ActualFormat
    !
    if(present(NumberFormat))then
       allocate(ActualFormat,source=trim(NumberFormat))
    else
       allocate(ActualFormat,source=DEFAULT_FORMAT)
    endif
    NCol=DEFAULT_NCOLUMNS
    if(present(NColumns))then
       if(NColumns<1)then
          call ErrorMessage("WriteEigenvalues: Wrong number of columns")
          return
       endif
       NCol=NColumns
    endif
    !
    INQUIRE( uid, form=form, iostat=iostat, iomsg=iomsg )
    if(iostat/=0)then
       call ErrorMessage("WriteEigenvalues: "//trim(iomsg))
       return
    endif
    if(trim(form)/="FORMATTED")then
       call ErrorMessage("WriteEigenvalues: wrong uid format")
       return
    endif
    !
    write(uid,*) SpecRes.Dim,SpecRes.NEigenvalues
    do iEval = 1, SpecRes.NEigenvalues
       !
       if(mod(iEval-1,NCol)==0)then
          if(iEval/=1) write(uid,*)
          write(uid,"(i5)",advance="no") iEval
       endif
       write(uid,"(x,"//ActualFormat//")",advance="no") SpecRes.EigenValues(iEval) 
       !
    enddo
    write(uid,*)
    !
    deallocate(ActualFormat)
  end subroutine WriteComplexEigenvaluesToUnit

  !> Writes in a file the complex spectral resolution class eigenvalues.
  subroutine WriteComplexEigenvaluesToFile( SpecRes, FileName, NColumns, NumberFormat )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    character(len=*)              , intent(in) :: FileName
    integer, optional             , intent(in) :: NColumns
    character(len=*), optional    , intent(in) :: NumberFormat
    !
    integer :: uid, iostat
    character(len=IOMSG_LENGTH) :: iomsg
    !
    open(newunit =  uid        , &
         file    =  FileName   , &
         form    = "formatted" , &
         status  = "unknown"   , &
         action  = "write"     , &
         iostat  =  iostat     , &
         iomsg   =  iomsg      )
    if(iostat/=0)call assert(iomsg)
    if(present(NColumns))then
       if(present(NumberFormat))then
          call SpecRes.WriteEigenvalues( uid, NColumns = NColumns, NumberFormat = NumberFormat )
       else
          call SpecRes.WriteEigenvalues( uid, NColumns = NColumns )
       endif
    else
       if(present(NumberFormat))then
          call SpecRes.WriteEigenvalues( uid, NumberFormat = NumberFormat )
       else
          call SpecRes.WriteEigenvalues( uid )
       endif
    endif
    close( uid, iostat=iostat, iomsg=iomsg )
    if(iostat/=0)call assert(iomsg)
    !
  end subroutine WriteComplexEigenvaluesToFile


  !> Gets a square submatrix from ClassComplexMatrix's matrix defined by a minimum and a maximum index which are assumed to be equal for both rows and columns.
  subroutine ClassComplexMatrixGetSubMatrix( Matrix, MinSubIndex, MaxSubIndex, SubMatrix )!@
    !
    Class(ClassComplexMatrix), intent(in) :: Matrix  
    integer,            intent(in)    :: MinSubIndex
    integer,            intent(in)    :: MaxSubIndex
    type(ClassComplexMatrix),  intent(out)   :: SubMatrix
    !
    integer         :: SubDim, NL, NU
    integer         :: iRow, iCol, iSubRow, iSubCol, iSubRowMin, iSubRowMax
    Complex(kind(1d0)) :: Element
    !
    if( MinSubIndex <= 0 )call Assert("Invalid MinSubIndex")
    if( MaxSubIndex > min( Matrix.NR, Matrix.NC ) ) &
         call Assert("Invalid MaxSubIndex")
    !
    SubDim = MaxSubIndex - MinSubIndex + 1
    NL = Matrix.LowerBandwidth()
    NU = Matrix.UpperBandwidth()
    !
    if( Matrix.IsFull() )then
       call SubMatrix.InitFull( SubDim, SubDim )
    elseif( Matrix.IsBanded() )then
       call SubMatrix.InitBanded( SubDim, SubDim, NL, NU )
    else
       call Assert("Unrecognized Matrix type")
    endif
    !
    do iSubCol = 1, SubDim
       iCol = ( MinSubIndex - 1 ) + iSubCol
       iSubRowMin = max(1,iSubCol-NU)
       iSubRowMax = min(SubDim,iSubCol+NL)
       do iSubRow = iSubRowMin, iSubRowMax 
          iRow = ( MinSubIndex - 1 ) + iSubRow
          Element = Matrix.Element( iRow, iCol )
          call SubMatrix.SetElement( iSubRow, iSubCol, Element )
       enddo
    enddo
    !
  end subroutine ClassComplexMatrixGetSubMatrix



  !> Gets a rectangular matrix which is extracted from an original complex one, depending on the four indices specified, two that set the row interval and the other two that set the column interval.
  subroutine ClassComplexMatrixGetSubMatrix_Rect( &
       Matrix     , &
       MinRowIndex, &
       MaxRowIndex, &
       MinColIndex, &
       MaxColIndex, &
       SubMatrix    )
    !
    Class(ClassComplexMatrix), intent(in) :: Matrix  
    integer,            intent(in)    :: MinRowIndex
    integer,            intent(in)    :: MaxRowIndex
    integer,            intent(in)    :: MinColIndex
    integer,            intent(in)    :: MaxColIndex
    type(ClassComplexMatrix),  intent(out)   :: SubMatrix
    !
    integer :: i, j
    !
    call SubMatrix.InitFull( MaxRowIndex-MinRowIndex+1,MaxColIndex-MinColIndex+1 )
    do j = MinColIndex, MaxColIndex
       do i = MinRowIndex, MaxRowIndex
          SubMatrix.A(i-MinRowIndex+1,j-MinColIndex+1) = Matrix.A(i,j)
       end do
    end do
    !
  end subroutine ClassComplexMatrixGetSubMatrix_Rect



  !> Multiplies the complex matrix by complex number.
  subroutine ClassComplexMatrixTimesDouble( Matrix, Number )
    Class(ClassComplexMatrix), intent(inout) :: Matrix 
    Complex(kind(1d0)),    intent(in)    :: Number
    Matrix.A = Number * Matrix.A
  end subroutine ClassComplexMatrixTimesDouble


  !> Adds up two complex matrices belonging to two different complex matrix classes.
  subroutine ClassComplexMatrixAddClassComplexMatrix( Matrix, DeltaMatrix ) 
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    Class(ClassComplexMatrix), intent(in)    :: DeltaMatrix
    !
    integer :: LBR, UBR, LBC, UBC, iRow,iCol
    !
    LBR=LBOUND(Matrix.A,1)
    UBR=UBOUND(Matrix.A,1)
    LBC=LBOUND(Matrix.A,2)
    UBC=UBOUND(Matrix.A,2)
    !
    if( DeltaMatrix.HasSameShapeAs(Matrix) )then
!!$       Matrix.A = Matrix.A + DeltaMatrix.A
       do iCol=LBC, UBC
          do iRow=LBR,UBR
             Matrix.A(iRow,iCol)=Matrix.A(iRow,iCol)+DeltaMatrix.A(iRow,iCol)
          enddo
       enddo
    else
       call Assert("Incompatible ClassComplexMatrix shapes")
    endif
  end subroutine ClassComplexMatrixAddClassComplexMatrix



  !> Adds up two complex matrices belonging to two different matrix classes.
  subroutine ClassComplexMatrixAddClassMatrix( Matrix, DeltaMatrix ) 
    Class(ClassComplexMatrix), intent(inout) :: Matrix
    Class(ClassMatrix)       , intent(in)    :: DeltaMatrix
    !
    type(ClassComplexMatrix) :: AuxMat
    AuxMat = DeltaMatrix
    call Matrix.Add( AuxMat )
    call AuxMat.Free()
    !
  end subroutine ClassComplexMatrixAddClassMatrix


  !> Returns True if two complex matrices have the same number of rows, columns, lower subdiagonal and superdiagonals.
  logical function ClassComplexMatrixHasSameShapeAs( Mat1, Mat2 ) result(Same)
    Class(ClassComplexMatrix), intent(in) :: Mat1, Mat2
    Same=.FALSE.
    if( Mat1.Pattern /= Mat2.Pattern ) return
    if(  Mat1.NR /= Mat2.NR .or. &
         Mat1.NC /= Mat2.NC .or. &
         Mat1.NL /= Mat2.NL .or. &
         Mat1.NU /= Mat2.NU )return
    Same=.TRUE.
    return
  end function ClassComplexMatrixHasSameShapeAs


  !> Fetches the eigenvalues from ClassComplexSpectralResolution.
  subroutine FetchComplexEigenvalues( SpecRes, Vector )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    Complex(kind(1d0)), allocatable,   intent(out) :: Vector(:)
    allocate( Vector, source = SpecRes.EigenValues )
  end subroutine FetchComplexEigenvalues


  !> Fetches the right eigenvectors from ClassComplexSpectralResolution.
  subroutine FetchComplexEigenvectors( SpecRes, Mat )
    Class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    Complex(kind(1d0)), allocatable,   intent(out) :: Mat(:,:)
    allocate( Mat, source = SpecRes.RightEigenvectors )
  end subroutine FetchComplexEigenvectors


  !> Fetches the right eigenvectors from ClassComplexSpectralResolution in a ClassComplexMatrix form.
  subroutine FetchComplexEigenvectorsMat( SpecRes, Mat )
    class(ClassComplexSpectralResolution), intent(in) :: SpecRes
    class(ClassComplexMatrix)            , intent(out) :: Mat
    Mat = SpecRes.RightEigenvectors
  end subroutine FetchComplexEigenvectorsMat

  
  !> Fetches the requested single eigenvector from ClassComplexSpectralResolution.
  subroutine FetchSingleComplexEigenvector( SpecRes, Vec, n )
    Class(ClassComplexSpectralResolution), intent(in)  :: SpecRes
    Complex(kind(1d0)), allocatable      , intent(out) :: Vec(:)
    integer                              , intent(in)  :: n
    integer :: Dim
    if(n>SpecRes.Neval())call ErrorMessage("Requested eigenvector doesn't exist")
    Dim = size(SpecRes.RightEigenvectors,1)
    if(allocated(Vec))then
       if(size(Vec,1)/=Dim)then
          deallocate(Vec)
          allocate(Vec(Dim))
       endif
    else
       allocate(Vec(Dim))
    endif
    Vec = SpecRes.RightEigenvectors(:,n)
  end subroutine FetchSingleComplexEigenvector


  !> Fetches the requested single eigenvector from ClassComplexSpectralResolution as a ClassComplexMatrix.
  subroutine FetchSingleComplexEigenvectorMat( SpecRes, Mat, n )
    Class(ClassComplexSpectralResolution), intent(in)  :: SpecRes
    Class(ClassComplexMatrix)            , intent(out) :: Mat
    integer                              , intent(in)  :: n
    integer :: Dim
    if(n>SpecRes.Neval())call ErrorMessage("Requested eigenvector doesn't exist")
    Dim = size(SpecRes.RightEigenvectors,1)
    call Mat.InitFull( Dim, 1 )
    Mat.A(:,1) = SpecRes.RightEigenvectors(:,n)
  end subroutine FetchSingleComplexEigenvectorMat


  !> Transforms the eigenvectors C of the spectral resolution on the basis defined by the columns U of NewBasis:
  !!\f[
  !! C'=U^{T}\cdot C,
  !!\f]
  subroutine TransformComplexEigenvectors_C( SpecRes, NewBasis )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassComplexSpectralResolution), intent(in)    :: NewBasis
    Complex(kind(1d0)), allocatable :: Matrix(:,:)
    integer :: NR,NC,NK
    NR=NewBasis.NEigenvalues
    NC=SpecRes.NEigenvalues
    NK=min(SpecRes.Dim,NewBasis.Dim)
    allocate(Matrix(NR,NC))
    Matrix=0.d0
    call ZGEMM( "C", "N", NR, NC, NK, (1.d0,0.d0),&
         NewBasis.RightEigenvectors, NewBasis.Dim, &
         SpecRes.RightEigenvectors,  SpecRes.Dim, &
         (0.d0,0.d0), Matrix, NR )
    deallocate(SpecRes.RightEigenvectors)
    allocate(SpecRes.RightEigenvectors,source=Matrix)
    SpecRes.Dim=NR
    deallocate(Matrix)
  end subroutine TransformComplexEigenvectors_C
  !
  subroutine TransformComplexEigenvectors_D( SpecRes, NewBasis_D )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassSpectralResolution), intent(in) :: NewBasis_D
    Complex(kind(1d0)), allocatable :: Matrix(:,:)
    integer :: NR,NC,NK
    NR=NewBasis_D.NEigenvalues
    NC=SpecRes.NEigenvalues
    NK=min(SpecRes.Dim,NewBasis_D.Dim)
    allocate(Matrix(NR,NC))
    Matrix=0.d0
    call ZGEMM( "T", "N", NR, NC, NK, (1.d0,0.d0),&
         NewBasis_D.Eigenvectors, NewBasis_D.Dim, &
         SpecRes.RightEigenvectors,  SpecRes.Dim, &
         (0.d0,0.d0), Matrix, NR )
    deallocate(SpecRes.RightEigenvectors)
    allocate(SpecRes.RightEigenvectors,source=Matrix)
    SpecRes.Dim=NR
    deallocate(Matrix)
  end subroutine TransformComplexEigenvectors_D


  !> transform the eigenvectors C of the spectral resolution
  !! on the basis defined by the columns U of NewBasis, taking
  !! into account the metric S of the basis:
  !!\f[
  !! C'=U^{T}\cdot S\cdot C,
  !!\f]
  subroutine TransformComplexEigenvectorsWithMetric_CC( SpecRes, NewBasis_C, Metric_C )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassComplexSpectralResolution), intent(in)    :: NewBasis_C
    Class(ClassComplexMatrix)            , intent(in)    :: Metric_C
    call SpecRes.Transform(Metric_C)
    call SpecRes.Transform(NewBasis_C)
  end subroutine TransformComplexEigenvectorsWithMetric_CC
  !
  subroutine TransformComplexEigenvectorsWithMetric_CD( SpecRes, NewBasis_C, Metric_D )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassComplexSpectralResolution), intent(in)    :: NewBasis_C
    Class(ClassMatrix)                   , intent(in)    :: Metric_D
    call SpecRes.Transform(Metric_D)
    call SpecRes.Transform(NewBasis_C)
  end subroutine TransformComplexEigenvectorsWithMetric_CD
  !
  subroutine TransformComplexEigenvectorsWithMetric_DC( SpecRes, NewBasis_D, Metric_C )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassSpectralResolution)       , intent(in)    :: NewBasis_D
    Class(ClassComplexMatrix)            , intent(in)    :: Metric_C
    call SpecRes.Transform(Metric_C)
    call SpecRes.Transform(NewBasis_D)
  end subroutine TransformComplexEigenvectorsWithMetric_DC
  !
  subroutine TransformComplexEigenvectorsWithMetric_DD( SpecRes, NewBasis_D, Metric_D )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassSpectralResolution)       , intent(in)    :: NewBasis_D
    Class(ClassMatrix)                   , intent(in)    :: Metric_D
    call SpecRes.Transform(Metric_D)
    call SpecRes.Transform(NewBasis_D)
  end subroutine TransformComplexEigenvectorsWithMetric_DD



  !> Multiply from the left the right eigenvectors of spectral resolution by  a matrix in ClassComplexMatrix form.
  subroutine TransformComplexEigenvectors_ClassComplexMatrix( SpecRes, LeftFactor )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassComplexMatrix)            , intent(in)    :: LeftFactor
    !
    Complex(kind(1d0)), allocatable :: Matrix(:,:)
    Complex(kind(1d0)), allocatable :: vec(:)
    integer :: iEval
    !
    if( LeftFactor.NColumns() /= SpecRes.Size() )call Assert("Unmatching dimensions in ClassMatrix x SpectralResolution")
    !
    select case( LeftFactor.Pattern )
    case( MATRIX_PATTERN_FULL )
       allocate(Matrix(SpecRes.Dim,SpecRes.NEigenvalues))
       Matrix=0.d0
       call ZGEMM("N","N",SpecRes.Dim,SpecRes.NEigenvalues,SpecRes.Dim,(1.d0,0.d0),&
            LeftFactor.A,LeftFactor.NR,&
            SpecRes.RightEigenvectors,SpecRes.Dim,&
            (0.d0,0.d0),Matrix,SpecRes.Dim)
       deallocate(SpecRes.RightEigenvectors)
       allocate(SpecRes.RightEigenvectors,source=Matrix)
       deallocate(Matrix)
    case( MATRIX_PATTERN_BANDED )
       allocate(Vec(SpecRes.Dim))
       Vec=0.d0
       do iEval=1,SpecRes.NEigenvalues
          call ZGBMV("N",SpecRes.Dim,SpecRes.Dim,LeftFactor.NL,LeftFactor.NU,1.d0,&
               LeftFactor.A,LeftFactor.NL+LeftFactor.NU+1,&
               SpecRes.RightEigenvectors(1,iEval),1,&
               0.d0,Vec,1)
          SpecRes.RightEigenvectors(:,iEval)=Vec
       enddo
       deallocate(Vec)
    case DEFAULT
       call Assert('Error: unrecognized matrix pattern')
    end select
  end subroutine TransformComplexEigenvectors_ClassComplexMatrix


  !> Multiply from the left the right eigenvectors of spectral resolution by  a matrix in ClassMatrix form.
  subroutine TransformComplexEigenvectors_ClassMatrix( SpecRes, LeftFactor )
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    Class(ClassMatrix)                   , intent(in)    :: LeftFactor
    !
    Complex(kind(1d0)), allocatable :: Matrix(:,:)
    Type(ClassComplexMatrix) :: LeftFactor_C, CompMat
    !
    if( LeftFactor.NColumns() /= SpecRes.Size() )call Assert("Unmatching dimensions in ClassMatrix x SpectralResolution")
    !
    select case( LeftFactor.Pattern )
    case( MATRIX_PATTERN_FULL )
       CompMat = SpecRes.RightEigenvectors
       call CompMat.Multiply( LeftFactor, 'Left', 'N' )
       call CompMat.FetchMatrix( Matrix )
       call CompMat.Free()
       !.. For some reason these commented lines don't work when it should.
       ! It has been replaced by the lines above
!!$       allocate(Matrix(SpecRes.Dim,SpecRes.NEigenvalues))
!!$       Matrix=0.d0
!!$       call ZGEMM("N","N",SpecRes.Dim,SpecRes.NEigenvalues,SpecRes.Dim,(1.d0,0.d0),&
!!$            LeftFactor.A,LeftFactor.NR,&
!!$            SpecRes.RightEigenvectors,SpecRes.Dim,&
!!$            (0.d0,0.d0),Matrix,SpecRes.Dim)
       deallocate(SpecRes.RightEigenvectors)
       allocate(SpecRes.RightEigenvectors,source=Matrix)
       deallocate(Matrix)
    case( MATRIX_PATTERN_BANDED )
       LeftFactor_C = LeftFactor
       call SpecRes.Transform(LeftFactor_C)
       call LeftFactor_C.Free()
    case DEFAULT
       call Assert('Error: unrecognized matrix pattern')
    end select
  end subroutine TransformComplexEigenvectors_ClassMatrix


  !> Returns True if both matrices have the same number of rows and columns.
  logical function ClassMatrixSameSize( Matrix1, Matrix2 )
    Class(ClassMatrix), intent(in)  :: Matrix1
    Class(ClassMatrix), intent(in)  :: Matrix2
    !
    ClassMatrixSameSize = .FALSE.
    if ( (Matrix1.NR == Matrix2.NR) .and. (Matrix1.NC == Matrix2.NC) ) then
       ClassMatrixSameSize = .TRUE. 
    end if
  end function ClassMatrixSameSize


  !> Saves the eigenvectors in ClassSpectralResolution as a matrix in ClassMatrix.  Only for Full pattern.
  subroutine FetchClassMatrixEigenvectors( SpecRes, Matrix )
    Class(ClassSpectralResolution), intent(in)  :: SpecRes
    Type(ClassMatrix),              intent(out) :: Matrix
    !
    integer :: Dim, Neval
    Dim   = SpecRes.Size()
    Neval = SpecRes.Neval()
    call Matrix.InitFull( Dim, Neval )
    Matrix.A = SpecRes.EigenVectors(1:Dim,1:Neval)
  end subroutine FetchClassMatrixEigenvectors



  !> Convert a banded stored rectangular array to the original squared  matrix.
  !! It is considered that the number of subdiagonals is equal to the number of superdiagonals.
  subroutine ConvertToSquared( Matrix, ResMat )
    Class(ClassMatrix), intent(inout) :: Matrix
    Class(ClassMatrix), optional, intent(out) :: ResMat
    !
    type(ClassMatrix) :: B
    integer :: i, j
    logical :: OnlySymmetric
    !
    OnlySymmetric = .false.
    ! 
    call B.InitFull( Matrix.NC, Matrix.NC )
    if ( Matrix.Pattern == MATRIX_PATTERN_BANDED ) then
       !
       if ( OnlySymmetric ) then
          ! only for symmetric original squared matrices
          do j = 1, Matrix.NC
             do i = j, j + min( Matrix.NL, Matrix.NC-j )
                !
                B.A(i,j) = Matrix.A(i-j,j)
                B.A(j,i) = B.A(i,j)
                !
             end do
          end do
       else
          do j = 1, Matrix.NC
             do i = max(1,j-Matrix.NL),  min(Matrix.NC,j+Matrix.NL)
                !
                B.A(i,j) = Matrix.A(i-j,j)
                !
             end do
          end do
          !
       end if
    elseif ( Matrix.Pattern == MATRIX_PATTERN_FULL ) then
!!$       B.A = Matrix.A
       B = Matrix
    else
       call Assert( 'Non-identified matrix type' )
    end if
    !
    if ( present(ResMat) ) then
       ResMat = B
    end if
    !
    if ( .not.present(ResMat) ) then
       call Matrix.Free
       Matrix = B
    end if
    !
    call B.Free(  )
  end subroutine ConvertToSquared




  !> Sets a choosen spectral resolution eigenvalue equal to some external one.
  subroutine SetOneEigenValue( SpecRes, N, Value )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    integer,                        intent(in)    :: N
    DoublePrecision,                intent(in)    :: Value
    !
    SpecRes.EigenValues(N) = Value
  end subroutine SetOneEigenValue


  !> Sets the complete vector of eigenvalue equal to some external one.
  subroutine SetAllEigenValues( SpecRes, Vector )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    real(kind(1d0)),                intent(in)    :: Vector(:)
    !
    if ( allocated(SpecRes.EigenValues) ) deallocate( SpecRes.EigenValues )
    allocate( SpecRes.EigenValues, source = Vector )
  end subroutine SetAllEigenValues



  subroutine ClassSpectralResolutionSetEigenVectors( SpecRes, EigVec )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    real(kind(1d0)),                intent(in)    :: EigVec(:,:)
    if ( allocated(SpecRes.EigenVectors) ) deallocate( SpecRes.EigenVectors )
    allocate( SpecRes.EigenVectors, source = EigVec )
  end subroutine ClassSpectralResolutionSetEigenVectors


  subroutine ClassSpectralResolutionSetEigenVectorsMat( SpecRes, EigVecMat )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    class(ClassMatrix),           intent(in)      :: EigVecMat
    if ( allocated(SpecRes.EigenVectors) ) deallocate( SpecRes.EigenVectors )
    call EigVecMat.FetchMatrix( SpecRes.EigenVectors )
  end subroutine ClassSpectralResolutionSetEigenVectorsMat



  !> Builds a Matrix to be factorized that further will be used to find its null space.
  subroutine ClassMatrixCompose( Matrix, MatA, Num, MatB )
    class(ClassMatrix), intent(inout) :: Matrix
    type(ClassMatrix),  intent(in)    :: MatA
    DoublePrecision,    intent(in)    :: Num
    type(ClassMatrix),  intent(in)    :: MatB
    !
    Matrix = MatB
    call Matrix.Multiply( -Num )
    call Matrix.Add( MatA )
  end subroutine ClassMatrixCompose



  !> Builds a Complex Matrix to be factorized that further will be used to find its null space.
  subroutine ClassComplexMatrixCompose( Matrix, MatA, Num, MatB )
    class(ClassComplexMatrix), intent(inout) :: Matrix
    type(ClassComplexMatrix),  intent(in)    :: MatA
    DoublePrecision,    intent(in)    :: Num
    type(ClassComplexMatrix),  intent(in)    :: MatB
    !
    complex(kind(1d0)) :: ComplexNum
    !
    Matrix = MatB
    ComplexNum = dcmplx(-Num,0.d0)
    call Matrix.Multiply( ComplexNum )
    call Matrix.Add( MatA )
  end subroutine ClassComplexMatrixCompose



  !> Performs the LU factorization of the matrix. If the variable Homogeneous = True, as well as in general the matrix's rows and/or columns might be linear independent, the last row is set equal to zero to garantee linear dependence and solve the system directly. If not then the Inverse Iteration Method will be used later to solve the linear equations system.
  subroutine ClassMatrixFactorize( Matrix, ipiv, Homogeneous, equed, Rvec, Cvec, AFB )
    class(ClassMatrix),                      intent(inout) :: Matrix
    integer,                    allocatable, intent(out)   :: ipiv(:)
    logical,          optional,              intent(inout) :: Homogeneous
    character(len=1), optional,              intent(out)   :: equed
    DoublePrecision,  optional, allocatable, intent(out)   :: Rvec(:)
    DoublePrecision,  optional, allocatable, intent(out)   :: Cvec(:)
    DoublePrecision,  optional, allocatable, intent(out)   :: AFB(:,:)
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call LUFullClassMatrixFactorize( Matrix, ipiv, Homogeneous )
    case( MATRIX_PATTERN_BANDED )
       if ( .not.present(equed) ) then
          call LUBandClassMatrixFactorize( Matrix, ipiv, Homogeneous )
       else
          call OldLUBandFactorize( Matrix, ipiv, equed, Rvec, Cvec, AFB )
       end if
    case DEFAULT 
       call Assert('Error: unrecognized matrix pattern')
    end select
  end subroutine ClassMatrixFactorize


  !> Performs the LU factorization of the complex matrix. If the variable Homogeneous = True, as well as in general the matrix's rows and/or columns might be linear independent, the last row is set equal to zero to garantee linear dependence and solve the system directly. If not then the Inverse Iteration Method will be used later to solve the linear equations system.
  subroutine ClassComplexMatrixFactorize( Matrix, ipiv, Homogeneous )
    class(ClassComplexMatrix),               intent(inout) :: Matrix
    integer,                    allocatable, intent(out)   :: ipiv(:)
    logical,          optional,              intent(inout) :: Homogeneous
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call LUFullClassComplexMatrixFactorize( Matrix, ipiv, Homogeneous )
    case( MATRIX_PATTERN_BANDED )
       call LUBandClassComplexMatrixFactorize( Matrix, ipiv, Homogeneous )
    case DEFAULT 
       call Assert('Error: unrecognized matrix pattern')
    end select
    !
  end subroutine ClassComplexMatrixFactorize


  !> Performs the LU factorization of the matrix when it is in Full representation.
  subroutine LUFullClassMatrixFactorize( Matrix, ipiv, Homogeneous )
    class(ClassMatrix),         intent(inout) :: Matrix
    integer, allocatable,       intent(out)   :: ipiv(:)
    logical,          optional, intent(inout) :: Homogeneous
    !
    integer :: m, n, lda, info
    character(len=100) :: strn
    logical :: IsHomogeneous
    !
    if (.not. present(Homogeneous) ) then
       IsHomogeneous = .false.
    else
       IsHomogeneous = Homogeneous
    end if
    !
    m = size( Matrix.A, 1 )
    n = size( Matrix.A, 2 )
    lda = m
    allocate( ipiv(max(1,min(m,n))) )
    !
    if ( IsHomogeneous ) then
       ! Make the last row equal to zero. (Assure the solution of the homogeneous problem.
       Matrix.A(n,:) = 0.d0
    end if
    !
    call dgetrf( m, &
         n, &
         Matrix.A, &
         lda, &
         ipiv, &
         info  )
    !
    if ( info < 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    elseif ( info > 0 ) then
!!$      write(strn,*) info
!!$      call ErrorMessage(" The diagonal upper triangular matrix element "//trim(adjustl(strn))//" is 0. A system of linear equations cannot be solved using the 'U' factor. Changing the '0' for an infinitesimal ...  ")         
!!$      write(OUTPUT_UNIT,"(a)")
!!$      write(OUTPUT_UNIT,"(a)") ' Approximating factorization to avoid singular matrices...'
       Matrix.A(info, info) = 1.d-50
    end if
  end subroutine LUFullClassMatrixFactorize


  !> Performs the LU factorization of the complex matrix when it is in Full representation.
  subroutine LUFullClassComplexMatrixFactorize( Matrix, ipiv, Homogeneous )
    class(ClassComplexMatrix),  intent(inout) :: Matrix
    integer, allocatable,       intent(out)   :: ipiv(:)
    logical,          optional, intent(inout) :: Homogeneous
    !
    integer :: m, n, lda, info
    character(len=100) :: strn
    logical :: IsHomogeneous
    !
    if (.not. present(Homogeneous) ) then
       IsHomogeneous = .false.
    else
       IsHomogeneous = Homogeneous
    end if
    !
    m = size( Matrix.A, 1 )
    n = size( Matrix.A, 2 )
    lda = m
    allocate( ipiv(max(1,min(m,n))) )
    !
    if ( IsHomogeneous ) then
       ! Make the last row equal to zero. (Assure the solution of the homogeneous problem.
       Matrix.A(n,:) = dcmplx(0.d0,0.d0)
    end if
    !
    call zgetrf( m, &
         n, &
         Matrix.A, &
         lda, &
         ipiv, &
         info  )
    !
    if ( info < 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    elseif ( info > 0 ) then
!!$      write(strn,*) info
!!$      call ErrorMessage(" The diagonal upper triangular matrix element "//trim(adjustl(strn))//" is 0. A system of linear equations cannot be solved using the 'U' factor. Changing the '0' for an infinitesimal ...  ")         
!!$      write(OUTPUT_UNIT,"(a)")
!!$      write(OUTPUT_UNIT,"(a)") ' Approximating factorization to avoid singular matrices...'
       Matrix.A(info, info) = dcmplx(1.d-50,0.d0)
    end if
    !
  end subroutine LUFullClassComplexMatrixFactorize



  !> Performs the LU factorization of the matrix when it is in Banded representation.
  subroutine LUBandClassMatrixFactorize( Matrix, ipiv, Homogeneous )
    class(ClassMatrix),         intent(inout) :: Matrix
    integer, allocatable,       intent(out)   :: ipiv(:)
    logical, optional,          intent(inout) :: Homogeneous
    !
    integer :: m, n, kl, ku, lda, info, i
    character(len=100) :: strn
    DoublePrecision, allocatable :: B(:,:)
    logical :: IsHomogeneous
    !
    if (.not. present(Homogeneous) ) then
       IsHomogeneous = .false.
    else
       IsHomogeneous = Homogeneous
    end if
    !
    n = size( Matrix.A, 2 )
    kl = Matrix.NL
    ku = Matrix.NU
    m = n
    lda = 2*kl+ku+1
    !
    if ( IsHomogeneous ) then
       ! Make zero the corresponding elements in band storage in a way that the last row of the original squared matrix is zero.(Assure the solution of the homogeneous problem)
       do i = 0, kl
          Matrix.A(i,n-i) = 0.d0 
       end do
    end if
    !
    ! Rearrange the matrix to be passed to the LAPACK subroutine
    allocate( B(2*kl+ku+1, n) )
    B(kl+1:2*kl+ku+1,:) = Matrix.A(:,:)
    deallocate( Matrix.A )
    Matrix.NRmin = 1
    Matrix.NRmax = 2*kl+ku+1
    call AllocateMatrix( Matrix.A,   &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    Matrix.A(kl+1:2*kl+ku+1,:) = B(kl+1:2*kl+ku+1,:)
    deallocate( B )
    allocate( ipiv(max(1,min(m,n))) )
    !
    call dgbtrf( m, &
         n, &
         kl, &
         ku, &
         Matrix.A, &
         lda, &
         ipiv, &
         info  )
    !
    if ( info < 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    elseif ( info > 0 ) then
!!$      write(strn,*) info
!!$      call ErrorMessage(" The diagonal upper triangular matrix element "//trim(adjustl(strn))//" is 0. A system of linear equations cannot be solved using the 'U' factor. Changing the '0' for an infinitesimal ... ")
!!$      write(OUTPUT_UNIT,"(a)")
!!$      write(OUTPUT_UNIT,"(a)") ' Approximating factorization to avoid singular matrices...'
       Matrix.A(kl+ku+1, info) = 1.d-50
    end if
  end subroutine LUBandClassMatrixFactorize



  !> Performs the LU factorization of the complex matrix when it is in Banded representation.
  subroutine LUBandClassComplexMatrixFactorize( Matrix, ipiv, Homogeneous )
    class(ClassComplexMatrix),  intent(inout) :: Matrix
    integer, allocatable,       intent(out)   :: ipiv(:)
    logical, optional,          intent(inout) :: Homogeneous
    !
    integer :: m, n, kl, ku, lda, info, i
    character(len=100) :: strn
    complex(kind(1d0)), allocatable :: B(:,:)
    logical :: IsHomogeneous
    !
    if (.not. present(Homogeneous) ) then
       IsHomogeneous = .false.
    else
       IsHomogeneous = Homogeneous
    end if
    !
    n = size( Matrix.A, 2 )
    kl = Matrix.NL
    ku = Matrix.NU
    m = n
    lda = 2*kl+ku+1
    !
    if ( IsHomogeneous ) then
       ! Make zero the corresponding elements in band storage in a way that the last row of the original squared matrix is zero.(Assure the solution of the homogeneous problem)
       do i = 0, kl
          Matrix.A(i,n-i) = dcmplx(0.d0,0.d0) 
       end do
    end if
    !
    ! Rearrange the matrix to be passed to the LAPACK subroutine
    allocate( B(2*kl+ku+1, n) )
    B(kl+1:2*kl+ku+1,:) = Matrix.A(:,:)
    deallocate( Matrix.A )
    Matrix.NRmin = 1
    Matrix.NRmax = 2*kl+ku+1
    call AllocateMatrix( Matrix.A,   &
         Matrix.NRmin, Matrix.NRmax, &
         Matrix.NCmin, Matrix.NCmax  )
    Matrix.A(kl+1:2*kl+ku+1,:) = B(kl+1:2*kl+ku+1,:)
    deallocate( B )
    allocate( ipiv(max(1,min(m,n))) )
    !
    call zgbtrf( m, &
         n, &
         kl, &
         ku, &
         Matrix.A, &
         lda, &
         ipiv, &
         info  )
    !
    if ( info < 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    elseif ( info > 0 ) then
!!$      write(strn,*) info
!!$      call ErrorMessage(" The diagonal upper triangular matrix element "//trim(adjustl(strn))//" is 0. A system of linear equations cannot be solved using the 'U' factor. Changing the '0' for an infinitesimal ... ")
!!$      write(OUTPUT_UNIT,"(a)")
!!$      write(OUTPUT_UNIT,"(a)") ' Approximating factorization to avoid singular matrices...'
       Matrix.A(kl+ku+1, info) = dcmplx(1.d-50,0.d0)
    end if
  end subroutine LUBandClassComplexMatrixFactorize



  !> Performs the LU factorization of the matrix and at the same time solves the linear equations system. Routine used in the old version of the program.
  subroutine OldLUBandFactorize( Matrix, ipiv, equed, Rvec, Cvec, AFB )
    class(ClassMatrix),           intent(inout) :: Matrix
    integer, allocatable,         intent(out)   :: ipiv(:)
    character(len=1),             intent(out)   :: equed
    DoublePrecision, allocatable, intent(out)   :: Rvec(:)
    DoublePrecision, allocatable, intent(out)   :: Cvec(:)
    DoublePrecision, allocatable, intent(out)   :: AFB(:,:)
    !
    integer :: n, kl, ku, lda, info, nrhs, ldafb, ldb, ldx
    character(len=100) :: strn
    DoublePrecision, allocatable :: Bmat(:,:), Xmat(:,:), work(:)
    DoublePrecision :: Rcond, Ferr, Berr
    integer, allocatable :: iwork(:)
    !
    n = size( Matrix.A, 2 ) 
    kl = Matrix.NL
    ku = Matrix.NU
    nrhs = 0
    lda = kl + ku +1
    ldafb = 2*kl + ku + 1
    ldb = n
    ldx = n
    allocate( AFB(ldafb,n), Rvec(n), Cvec(n), Bmat(n,1), Xmat(n,1) )
    allocate( ipiv(n), work(3*n), iwork(n) )
    AFB  = 0.d0
    Rvec = 0.d0
    Cvec = 0.d0
    Bmat = 0.d0
    Xmat = 0.d0
    !
    call DGBSVX( 'E',       &
         'N',       &
         n,        &
         kl,       &
         ku,       &
         nrhs,     &
         Matrix.A, &
         lda,      &
         AFB,      &
         ldafb,    &
         ipiv,     &
         equed,    &
         Rvec,     &
         Cvec,     &
         Bmat,     &
         ldb,      &
         Xmat,     &
         ldx,      &
         Rcond,    &
         Ferr,     &
         Berr,     &
         work,     &
         iwork,    &
         INFO      )
    !
    !
    if ( info < 0 ) then
       !
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
       !
    elseif( info == size(AFB,2) + 1 ) then
       !
       write(*,"(a,d12.3,a)") " RCOND = ",RCOND," is smaller than machine precision"
       !
    elseif ( info > 0 ) then
       !
       AFB(kl+ku+1, info) = 1.d-50
       !
    end if
    !
    deallocate( Bmat, Xmat, work, iwork )
    !
  end subroutine OldLUBandFactorize



  !> Solves the linear equations system. 
  !! Supposes that the matrix was previously LU factorized. 
  !! If method is 
  !!  - "II", uses Inverse Iteration
  !!  - "LU", directly solves the homogeneous system.
  subroutine ClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, Method, EigenVector, equed, Rvec, Cvec, AFB )
    class(ClassMatrix),            intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    DoublePrecision,               intent(in)    :: VectorRHS(:,:)
    DoublePrecision,  allocatable, intent(out)   :: EigenVector(:)
    character(len=*),              intent(in)    :: Method
    character(len=1), optional,    intent(in)    :: equed
    DoublePrecision,  optional,    intent(in)    :: Rvec(:)
    DoublePrecision,  optional,    intent(in)    :: Cvec(:)
    DoublePrecision,  optional,    intent(in)    :: AFB(:,:)
    !
    character(len=*), parameter :: HERE="ClassMatrix::LinEqSolver : "
    !
    if ( Method .is. "LU" ) then
       select case( Matrix.Pattern )
       case( MATRIX_PATTERN_FULL )
          call FullClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
       case( MATRIX_PATTERN_BANDED )
          if ( .not.present(equed) ) then
             call BandClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
          else
             call OldBandLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector, equed, Rvec, Cvec, AFB )
          end if
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
    elseif ( Method .is. "II" ) then
       call ClassMatrixInverIterSolver( Matrix, ipiv, EigenVector, equed, Rvec, Cvec, AFB )
    else
       call Assert(HERE//" Invalid method '"//Method//"'")
    end if
    !
  end subroutine ClassMatrixLinEqSolver



  !> Solves the linear equations system. 
  !! Supposes that the matrix was previously LU factorized. 
  !! If method is 
  !!  - "II", uses Inverse Iteration
  !!  - "LU", directly solves the homogeneous system.
  subroutine ClassMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, Method, EigenVectors, equed, Rvec, Cvec, AFB )
    class(ClassMatrix),            intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    DoublePrecision,               intent(in)    :: VectorRHS(:,:)
    DoublePrecision,  allocatable, intent(out)   :: EigenVectors(:,:)
    character(len=*),              intent(in)    :: Method
    character(len=1), optional,    intent(in)    :: equed
    DoublePrecision,  optional,    intent(in)    :: Rvec(:)
    DoublePrecision,  optional,    intent(in)    :: Cvec(:)
    DoublePrecision,  optional,    intent(in)    :: AFB(:,:)
    !
    character(len=*), parameter :: HERE="ClassMatrix::LinEqSolver : "
    !
    if ( Method .is. "LU" ) then
       select case( Matrix.Pattern )
       case( MATRIX_PATTERN_FULL )
          call FullClassMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, EigenVectors )
       case( MATRIX_PATTERN_BANDED )
          if ( .not.present(equed) ) then
             call BandClassMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, EigenVectors )
          else
!!$             call OldBandLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector, equed, Rvec, Cvec, AFB )
          end if
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
    elseif ( Method .is. "II" ) then
!!$       call ClassMatrixInverIterSolver( Matrix, ipiv, EigenVector, equed, Rvec, Cvec, AFB )
    else
       call Assert(HERE//" Invalid method '"//Method//"'")
    end if
    !
  end subroutine ClassMatrixLinEqSolverGeneral



  !> Solves the linear equations system. 
  !! Supposes that the complex matrix was previously LU factorized. 
  !! If method is 
  !!  - "II", uses Inverse Iteration
  !!  - "LU", directly solves the homogeneous system.
  subroutine ClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS,Method, EigenVector )
    class(ClassComplexMatrix),     intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    complex(kind(1d0)),            intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),allocatable,intent(out)   :: EigenVector(:)
    character(len=*),              intent(in)    :: Method
    !
    character(len=*), parameter :: HERE="ClassMatrix::LinEqSolver : "
    !
    if ( Method .is. "LU" ) then
       !
       select case( Matrix.Pattern )
       case( MATRIX_PATTERN_FULL )
          call FullClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
       case( MATRIX_PATTERN_BANDED )
          call BandClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
       !
    elseif ( Method .is. "II" ) then
       !
       call ClassComplexMatrixInverIterSolver( Matrix, ipiv, EigenVector )
       !
    else
       !
       call Assert(HERE//" Invalid method '"//Method//"'")
       !
    end if
    !
  end subroutine ClassComplexMatrixLinEqSolver



  !> Solves the linear equations system. 
  !! Supposes that the complex matrix was previously LU factorized. 
  !! If method is 
  !!  - "II", uses Inverse Iteration
  !!  - "LU", directly solves the homogeneous system.
  subroutine ClassComplexMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS,Method, EigenVectors )
    class(ClassComplexMatrix),     intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    complex(kind(1d0)),            intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),allocatable,intent(out)   :: EigenVectors(:,:)
    character(len=*),              intent(in)    :: Method
    !
    character(len=*), parameter :: HERE="ClassMatrix::LinEqSolver : "
    !
    if ( Method .is. "LU" ) then
       !
       select case( Matrix.Pattern )
       case( MATRIX_PATTERN_FULL )
          call FullClassComplexMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, EigenVectors )
       case( MATRIX_PATTERN_BANDED )
          call BandClassComplexMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, EigenVectors )
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
       !
    elseif ( Method .is. "II" ) then
       !
!!$       call ClassComplexMatrixInverIterSolver( Matrix, ipiv, EigenVectors )
       !
    else
       !
       call Assert(HERE//" Invalid method '"//Method//"'")
       !
    end if
    !
  end subroutine ClassComplexMatrixLinEqSolverGeneral



  !> Solves the linear equations system. Supposes that the matrix was previously 
  !! LU factorized, that it is in Full representation and that it is squared.
  subroutine FullClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector ) 
    class(ClassMatrix),           intent(inout) :: Matrix
    integer,                      intent(in)    :: ipiv(:)
    DoublePrecision,              intent(in)    :: VectorRHS(:,:)
    DoublePrecision, allocatable, intent(out)   :: EigenVector(:)
    !   
    integer :: n, nrhs, lda, ldb, info
    character(len=100) :: strn
    DoublePrecision, allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 1 )
    nrhs = size(VectorRHS,2)
    lda = n
    ldb = lda
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call dgetrs( 'N', &
         n, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVector, source = VectorRHS2(:,1) )
    !
  end subroutine FullClassMatrixLinEqSolver



  !> Solves the linear equations system. Supposes that the matrix was previously 
  !! LU factorized, that it is in Full representation and that it is squared.
  subroutine FullClassMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, EigenVectors ) 
    class(ClassMatrix),           intent(inout) :: Matrix
    integer,                      intent(in)    :: ipiv(:)
    DoublePrecision,              intent(in)    :: VectorRHS(:,:)
    DoublePrecision, allocatable, intent(out)   :: EigenVectors(:,:)
    !   
    integer :: n, nrhs, lda, ldb, info
    character(len=100) :: strn
    DoublePrecision, allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 1 )
    nrhs = size(VectorRHS,2)
    lda = n
    ldb = lda
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call dgetrs( 'N', &
         n, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVectors, source = VectorRHS2 )
    !
  end subroutine FullClassMatrixLinEqSolverGeneral



  !> Solves the linear equations system. Supposes that the complex matrix was previously 
  !! LU factorized, that it is in Full representation and that it is squared.
  subroutine FullClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector ) 
    class(ClassComplexMatrix),     intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    complex(kind(1d0)),            intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),allocatable,intent(out)   :: EigenVector(:)
    !   
    integer :: n, nrhs, lda, ldb, info, i
    character(len=100) :: strn
    complex(kind(1d0)), allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 1 )
    nrhs = size(VectorRHS,2)
    lda = n
    ldb = lda
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call zgetrs( 'N', &
         n, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVector, source = VectorRHS2(:,1) )
    !
  end subroutine FullClassComplexMatrixLinEqSolver



  !> Solves the linear equations system. Supposes that the complex matrix was previously 
  !! LU factorized, that it is in Full representation and that it is squared.
  subroutine FullClassComplexMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, EigenVectors ) 
    class(ClassComplexMatrix),     intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    complex(kind(1d0)),            intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),allocatable,intent(out)   :: EigenVectors(:,:)
    !   
    integer :: n, nrhs, lda, ldb, info, i
    character(len=100) :: strn
    complex(kind(1d0)), allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 1 )
    nrhs = size(VectorRHS,2)
    lda = n
    ldb = lda
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call zgetrs( 'N', &
         n, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVectors, source = VectorRHS2 )
    !
  end subroutine FullClassComplexMatrixLinEqSolverGeneral



  ! Solves the linear equations system. Supposes that the matrix was previously LU factorized and that has Banded pattern.
  subroutine BandClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, Eigenvector ) 
    class(ClassMatrix),           intent(inout) :: Matrix
    integer,                      intent(in)    :: ipiv(:)
    DoublePrecision,              intent(in)    :: VectorRHS(:,:)
    DoublePrecision, allocatable, intent(out)   :: EigenVector(:)
    !   
    integer :: n, kl, ku, nrhs, lda, ldb, info
    character(len=100) :: strn
    DoublePrecision, allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 2 )
    kl = Matrix.NL
    ku = Matrix.NU
    nrhs = size(VectorRHS,2)
    lda = 2*kl+ku+1
    ldb = n
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call dgbtrs( 'N', &
         n, &
         kl, &
         ku, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVector, source = VectorRHS2(:,1) )
    !
  end subroutine BandClassMatrixLinEqSolver



  ! Solves the linear equations system. Supposes that the matrix was previously LU factorized and that has Banded pattern.
  subroutine BandClassMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, Eigenvectors ) 
    class(ClassMatrix),           intent(inout) :: Matrix
    integer,                      intent(in)    :: ipiv(:)
    DoublePrecision,              intent(in)    :: VectorRHS(:,:)
    DoublePrecision, allocatable, intent(out)   :: EigenVectors(:,:)
    !   
    integer :: n, kl, ku, nrhs, lda, ldb, info
    character(len=100) :: strn
    DoublePrecision, allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 2 )
    kl = Matrix.NL
    ku = Matrix.NU
    nrhs = size(VectorRHS,2)
    lda = 2*kl+ku+1
    ldb = n
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call dgbtrs( 'N', &
         n, &
         kl, &
         ku, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVectors, source = VectorRHS2 )
    !
  end subroutine BandClassMatrixLinEqSolverGeneral



  ! Solves the linear equations system. Supposes that the complex matrix was previously LU factorized and that has Banded pattern.
  subroutine BandClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS, Eigenvector ) 
    class(ClassComplexMatrix),     intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    complex(kind(1d0)),            intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),allocatable,intent(out)   :: EigenVector(:)
    !   
    integer :: n, kl, ku, nrhs, lda, ldb, info, i
    character(len=100) :: strn
    complex(kind(1d0)), allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 2 )
    kl = Matrix.NL
    ku = Matrix.NU
    nrhs = size(VectorRHS,2)
    lda = 2*kl+ku+1
    ldb = n
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call zgbtrs( 'N', &
         n, &
         kl, &
         ku, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVector, source = VectorRHS2(:,1) )
    !
  end subroutine BandClassComplexMatrixLinEqSolver



  ! Solves the linear equations system. Supposes that the complex matrix was previously LU factorized and that has Banded pattern.
  subroutine BandClassComplexMatrixLinEqSolverGeneral( Matrix, ipiv, VectorRHS, Eigenvectors ) 
    class(ClassComplexMatrix),     intent(inout) :: Matrix
    integer,                       intent(in)    :: ipiv(:)
    complex(kind(1d0)),            intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),allocatable,intent(out)   :: EigenVectors(:,:)
    !   
    integer :: n, kl, ku, nrhs, lda, ldb, info, i
    character(len=100) :: strn
    complex(kind(1d0)), allocatable :: VectorRHS2(:,:)
    !
    n = size( Matrix.A, 2 )
    kl = Matrix.NL
    ku = Matrix.NU
    nrhs = size(VectorRHS,2)
    lda = 2*kl+ku+1
    ldb = n
    !
    allocate(VectorRHS2(n,nrhs)) 
    VectorRHS2 = VectorRHS
    call zgbtrs( 'N', &
         n, &
         kl, &
         ku, &
         nrhs, &
         Matrix.A, &
         lda, &
         ipiv, &
         VectorRHS2, &
         ldb, &
         info   )
    !
    if (info /= 0 ) then
       write(strn,*) -info
       call Assert(" The "//trim(adjustl(strn))//" parameter had an illegal value.")
    end if
    !
    allocate( EigenVectors, source = VectorRHS2 )
    !
  end subroutine BandClassComplexMatrixLinEqSolverGeneral



  ! Solves the linear equations system. Supposes that the matrix was previously LU factorized. Routine used in the old program version.
  subroutine OldBandLinEqSolver( Matrix, ipiv, VectorRHS, Eigenvector, equed, Rvec, Cvec, AFB ) 
    class(ClassMatrix),           intent(inout) :: Matrix
    integer,                      intent(in)    :: ipiv(:)
    DoublePrecision,              intent(in)    :: VectorRHS(:,:)
    DoublePrecision, allocatable, intent(out)   :: EigenVector(:)
    character(len=1),             intent(in)    :: equed
    DoublePrecision,              intent(in)    :: Rvec(:)
    DoublePrecision,              intent(in)    :: Cvec(:)
    DoublePrecision,              intent(in)    :: AFB(:,:)
    !
    integer :: n, kl, ku, lda, info, nrhs, ldafb, ldb, ldx
    DoublePrecision, allocatable :: work(:)
    DoublePrecision :: Rcond, Ferr, Berr
    integer, allocatable :: iwork(:)
    !
    n = size( Matrix.A, 2 ) 
    kl = Matrix.NL
    ku = Matrix.NU
    nrhs = 1
    lda = kl + ku +1
    ldafb = 2*kl + ku + 1
    ldb = n
    ldx = n
    allocate( work(3*n), iwork(n), EigenVector(n) )
    !
    call DGBSVX( 'E',       &
         'N',       &
         n,        &
         kl,       &
         ku,       &
         nrhs,     &
         Matrix.A, &
         lda,      &
         AFB,      &
         ldafb,    &
         ipiv,     &
         equed,    &
         Rvec,     &
         Cvec,     &
         VectorRHS,     &
         ldb,      &
         EigenVector,     &
         ldx,      &
         Rcond,    &
         Ferr,     &
         Berr,     &
         work,     &
         iwork,    &
         INFO      )


  end subroutine OldBandLinEqSolver


  !> Set a choosen spectral resolution eigenvector equal to an external value.
  subroutine ClassSpectralResolutionSetOneEigenVector( SpecRes, N, Vector )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    integer,                        intent(in)    :: N
    DoublePrecision,                intent(in)    :: Vector(:)
    !
    SpecRes.EigenVectors(:,N) = 0.d0 
    SpecRes.EigenVectors(1:size(Vector),N) = Vector 
  end subroutine ClassSpectralResolutionSetOneEigenVector


  !> Solves the linear equations system using Inverse Iteration Method for Class Matrix.
  subroutine ClassMatrixInverIterSolver( Matrix, ipiv, EigenVector, equed, Rvec, Cvec, AFB )
    class(ClassMatrix),                      intent(inout) :: Matrix
    integer,                                 intent(in)    :: ipiv(:)
    DoublePrecision,            allocatable, intent(out)   :: EigenVector(:)
    character(len=1), optional,              intent(in)    :: equed
    DoublePrecision,  optional,              intent(in)    :: Rvec(:)
    DoublePrecision,  optional,              intent(in)    :: Cvec(:)
    DoublePrecision,  optional,              intent(in)    :: AFB(:,:)

    integer        , parameter :: MIN_NITER_TO_STATIONARY_POINT = 200
    real(kind(1d0)), parameter :: NORM_TOLERANCE = 1.d-20

    DoublePrecision, allocatable :: RandVec(:), IterVec(:,:)
    DoublePrecision              :: NormRandVec, NormEigenVector, w1, w2, w3
    Integer                      :: n, i
    !
    n = size( Matrix.A, 2 )
    allocate( RandVec(n), IterVec(n,1) )
    IterVec = 0.d0
    call random_number( RandVec )
    !
    NormRandVec = sqrt( sum( RandVec*RandVec ) )
    RandVec = RandVec/NormRandVec
    IterVec(:,1) = RandVec
    !
    i = 0
    w3 = 1.d0
    !
    do
       i = i+1
       !.. Back substitution
       !.. 
       call ClassMatrixSelectLinEqSolver( Matrix, ipiv, IterVec, EigenVector, equed, Rvec, Cvec, AFB )
       !.. Normalization
       !..
       NormEigenVector = sqrt(dot_product(EigenVector,EigenVector))
       EigenVector = EigenVector/NormEigenVector
       w1 = dot_product(IterVec(:,1),EigenVector)
       IterVec(:,1) = IterVec(:,1) - w1 * EigenVector
       w2 = sqrt( dot_product( IterVec(:,1), IterVec(:,1) ) )
       if( w2 < NORM_TOLERANCE )exit
       if( (i >= MIN_NITER_TO_STATIONARY_POINT ) .and. (w2 >= w3) ) then
          !
          call ErrorMessage("Inverse Iteration has not converged:")
          write(OUTPUT_UNIT,*) "i =", i, "norm =", w2
          return
          !
       end if
       w3 = w2
       IterVec(:,1) = EigenVector
       deallocate( EigenVector )
       !
    end do

    deallocate(RandVec,IterVec)

  end subroutine ClassMatrixInverIterSolver




  !> Solves the linear equations system using Inverse Iteration Method for Class Complex Matrix. Only for one right side vector at the time.
  subroutine ClassComplexMatrixInverIterSolver( Matrix, ipiv, EigenVector )
    class(ClassComplexMatrix),               intent(inout) :: Matrix
    integer,                                 intent(in)    :: ipiv(:)
    complex(kind(1d0)),         allocatable, intent(out)   :: EigenVector(:)

    integer        , parameter :: MIN_NITER_TO_STATIONARY_POINT = 200
    real(kind(1d0)), parameter :: NORM_TOLERANCE = 1.d-20

    complex(kind(1d0)), allocatable :: ComplexRandVec(:), IterVec(:,:)
    real(kind(1d0)),    allocatable :: RandVec1(:), RandVec2(:)
    DoublePrecision              :: NormComplexRandVec, NormEigenVector, w1, w2, w3
    Integer                      :: n, i
    !
    n = size( Matrix.A, 2 )
    allocate( RandVec1(n), RandVec2(n),ComplexRandVec(n),IterVec(n,1) )
    IterVec = dcmplx(0.d0,0.d0)
    !
    call random_number( RandVec1 )
    call random_number( RandVec2 )
    !
    do i =1, n
       ComplexRandVec(i) = dcmplx( RandVec1(i), RandVec2(i) )
    end do
    !
    NormComplexRandVec = sqrt( sum( conjg(ComplexRandVec)*ComplexRandVec ) )
    ComplexRandVec = ComplexRandVec/NormComplexRandVec
    !
    IterVec(:,1) = ComplexRandVec
    !
    i = 0
    w3 = 1.d0
    !
    do
       i = i+1
       !.. Back substitution
       !.. 
       call ClassComplexMatrixSelectLinEqSolver( Matrix, ipiv, IterVec, EigenVector )
       !.. Normalization
       !..
       NormEigenVector = sqrt( sum( conjg(EigenVector)*EigenVector ) )
       EigenVector = EigenVector/NormEigenVector
       w1 = sum( conjg(IterVec(:,1))*EigenVector(:) )
       IterVec(:,1) = IterVec(:,1) - w1 * EigenVector(:)
       w2 = sqrt( sum( conjg(IterVec(:,1)) * IterVec(:,1) ) )
       if( w2 < NORM_TOLERANCE )exit
       if( (i >= MIN_NITER_TO_STATIONARY_POINT ) .and. (w2 >= w3) ) then
          !
          call ErrorMessage("Inverse Iteration has not converged:")
          write(OUTPUT_UNIT,*) "i =", i, "norm =", w2
          return
          !
       end if
       w3 = w2
       IterVec(:,1) = EigenVector(:)
       deallocate( EigenVector )
       !
    end do
    !
    deallocate(RandVec1, RandVec2, ComplexRandVec, IterVec)
    !
  end subroutine ClassComplexMatrixInverIterSolver




  !> Selects which kind of linear equation solver will be used 
  !! related with the matrix pattern, when the method used is Inverse Iteration.
  subroutine ClassMatrixSelectLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector, equed, Rvec, Cvec, AFB )
    class(ClassMatrix),                      intent(inout) :: Matrix
    integer,                                 intent(in)    :: ipiv(:)
    DoublePrecision,                         intent(in)    :: VectorRHS(:,:)
    DoublePrecision,            allocatable, intent(out)   :: EigenVector(:)
    character(len=1), optional,              intent(in)    :: equed
    DoublePrecision,  optional,              intent(in)    :: Rvec(:)
    DoublePrecision,  optional,              intent(in)    :: Cvec(:)
    DoublePrecision,  optional,              intent(in)    :: AFB(:,:)   !   

    character(len=*), parameter :: HERE = "ClassMatrix::SelectLinEqSolver : "

    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call FullClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
    case( MATRIX_PATTERN_BANDED )
       if ( .not.present(equed) ) then
          call BandClassMatrixLinEqSolver( Matrix, ipiv, VectorRHS, Eigenvector )
       else
          call OldBandLinEqSolver( Matrix, ipiv, VectorRHS, Eigenvector, equed, Rvec, Cvec, AFB )
       end if
    case DEFAULT
       call Assert(HERE//"only full and banded pattern are permitted")
    end select

  end subroutine ClassMatrixSelectLinEqSolver



  !> Selects which kind of linear equation solver will be used 
  !! related with the complex matrix pattern, when the method used is Inverse Iteration.
  subroutine ClassComplexMatrixSelectLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
    class(ClassComplexMatrix),               intent(inout) :: Matrix
    integer,                                 intent(in)    :: ipiv(:)
    complex(kind(1d0)),                      intent(in)    :: VectorRHS(:,:)
    complex(kind(1d0)),         allocatable, intent(out)   :: EigenVector(:)
    !
    character(len=*), parameter :: HERE = "ClassMatrix::SelectLinEqSolver : "
    !
    select case( Matrix.Pattern )
    case( MATRIX_PATTERN_FULL )
       call FullClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS, EigenVector )
    case( MATRIX_PATTERN_BANDED )
       call BandClassComplexMatrixLinEqSolver( Matrix, ipiv, VectorRHS, Eigenvector )
    case DEFAULT
       call Assert(HERE//"only full and banded pattern are permitted")
    end select
    !
  end subroutine ClassComplexMatrixSelectLinEqSolver



  !> initialize a ClassMatrix with Full pattern and get as matrix an external one.
  subroutine ClassMatrixSetMatrix( Matrix, VectorArray, MatrixArray )
    Class(ClassMatrix),        intent(inout) :: Matrix
    DoublePrecision, optional, intent(in)    :: VectorArray(:)
    DoublePrecision, optional, intent(in)    :: MatrixArray(:,:)
    !
    if ( present( VectorArray ) .and. .not.present( MatrixArray ) ) then
       call Matrix.InitFull( size(VectorArray,1), 1 )
       Matrix.A(:,1) = VectorArray
    elseif ( .not.present( VectorArray ) .and. present( MatrixArray ) ) then
       call Matrix.InitFull( size(MatrixArray,1), size(MatrixArray,2) )
       Matrix.A = MatrixArray
    else
       call Assert("There must be passsed one and only one array")
    end if
  end subroutine ClassMatrixSetMatrix



  !> Fetches the vector in a ClassMatrix when it consists in a one dimensional array.
  subroutine FetchVector( Matrix, Array )
    Class(ClassMatrix),                     intent(inout) :: Matrix
    DoublePrecision, allocatable,           intent(out)   :: Array(:)
    !
    allocate( Array, source =  Matrix.A(:,1) )
  end subroutine FetchVector



  !> Fetches the matrix in a ClassMatrix.
  subroutine FetchMatrix( Matrix, Array )
    Class(ClassMatrix),                     intent(in) :: Matrix
    real(kind(1d0)), allocatable,           intent(out):: Array(:,:)
    !
    allocate( Array, source = Matrix.A )
  end subroutine FetchMatrix


  !> Fetches from a complete matrix representation, the
  !! submatrix whose rows or columns corresponds to the
  !! indices specified.
  subroutine FetchClassMatrixByIndices( Matrix, Indices, Label, OutMat )
    !> Original ClassMatrix.
    class(ClassMatrix), intent(in)  :: Matrix
    !> Vector of indices to extract the rows or
    !! columns and build the new ClassMatrix.
    integer           , intent(in)  :: Indices(:)
    !> It has to be either 'rows' or 'columns', to
    !! extract those given in the Indices vector.
    character(len=*)  , intent(in)  :: Label
    !> New ClassMatrix.
    class(ClassMatrix), intent(out) :: OutMat
    !
    real(kind(1d0)), allocatable :: OldArray(:,:), NewArray(:,:)
    integer :: NewDim, i
    !
    NewDim = size(Indices)
    call Matrix.FetchMatrix( OldArray )
    !
    if ( Label .is. 'ROWS' ) then
       !
       allocate( NewArray(NewDim,Matrix.NColumns()) )
       do i = 1, NewDim
          NewArray(i,:) = OldArray(Indices(i),:)
       end do
       !
    elseif ( Label .is. 'COLUMNS' ) then
       !
       allocate( NewArray(Matrix.NRows(),NewDim) )
       do i = 1, NewDim
          NewArray(:,i) = OldArray(:,Indices(i))
       end do
       !
    else
       call Assert( 'Invalid label, it has to be either '//&
            '"rows" or "columns".')
    end if
    !
    deallocate( OldArray )
    OutMat = NewArray
    deallocate( NewArray )
    !
  end subroutine FetchClassMatrixByIndices


  !> Frees the ClassSpectralResolution attributes.
  subroutine FreeSpectralResolution( SpecRes )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    !
    if ( allocated(SpecRes.EigenValues) ) deallocate(SpecRes.EigenValues)  
    if ( allocated(SpecRes.EigenVectors) ) deallocate(SpecRes.EigenVectors)
    SpecRes.Dim = 0
    SpecRes.NEigenvalues = 0
  end subroutine FreeSpectralResolution



  !> Fetches the complex matrix representation
  subroutine FetchComplexMatrix( Matrix, Array )
    Class(ClassComplexMatrix),      intent(in)  :: Matrix
    complex(kind(1d0)), allocatable, intent(out) :: Array(:,:)
    !
    allocate( Array, source = Matrix.A )
  end subroutine FetchComplexMatrix


  !> Fetches from a complete complex matrix representation, the
  !! submatrix whose rows or columns corresponds to the
  !! indices specified.
  subroutine FetchClassComplexMatrixByIndices( Matrix, Indices, Label, OutMat )
    !> Original ClassComplexMatrix.
    class(ClassComplexMatrix), intent(in)  :: Matrix
    !> Vector of indices to extract the rows or
    !! columns and build the new ClassMatrix.
    integer                  , intent(in)  :: Indices(:)
    !> It has to be either 'rows' or 'columns', to
    !! extract those given in the Indices vector.
    character(len=*)         , intent(in)  :: Label
    !> New ClassComplexMatrix.
    class(ClassComplexMatrix), intent(out) :: OutMat
    !
    complex(kind(1d0)), allocatable :: OldArray(:,:), NewArray(:,:)
    integer :: NewDim, i
    !
    NewDim = size(Indices)
    call Matrix.FetchMatrix( OldArray )
    !
    if ( Label .is. 'ROWS' ) then
       !
       allocate( NewArray(NewDim,Matrix.NColumns()) )
       do i = 1, NewDim
          NewArray(i,:) = OldArray(Indices(i),:)
       end do
       !
    elseif ( Label .is. 'COLUMNS' ) then
       !
       allocate( NewArray(Matrix.NRows(),NewDim) )
       do i = 1, NewDim
          NewArray(:,i) = OldArray(:,Indices(i))
       end do
       !
    else
       call Assert( 'Invalid label, it has to be either '//&
            '"rows" or "columns".')
    end if
    !
    deallocate( OldArray )
    OutMat = NewArray
    deallocate( NewArray )
    !
  end subroutine FetchClassComplexMatrixByIndices



  !> Checks if the Matrix has at least one zero row.
  logical function CheckZeroRowOrColumn(Matrix) result(res)
    !
    Class(ClassMatrix), intent(in)  :: Matrix
    !
    type(ClassMatrix) :: CopyMatrix
    integer :: i, j
    integer :: NumberZeros
    real(kind(1d0)), parameter :: Param = epsilon(1d0)
    !
    CopyMatrix = Matrix
    !
    call CopyMatrix.ConvertToSquared()
    !
    res = .false.
    !
    do j = 1, CopyMatrix.NColumns() 
       !
       NumberZeros = 0
       !
       do i = 1, CopyMatrix.NRows()
          !
          if ( abs(CopyMatrix.A(i,j)) < Param ) then
             !
             NumberZeros = NumberZeros + 1
             !
          end if
          !
       end do
       !
       if ( NumberZeros == CopyMatrix.NRows() ) then
          !
          res = .true.
          write(*,*) ' the column ', j, 'is zero'
          return
          !
       end if
    end do
    !
    !
    !
    do i = 1, CopyMatrix.NRows()
       !
       NumberZeros = 0
       !
       do j = 1, CopyMatrix.NColumns() 
          !
          if ( abs(CopyMatrix.A(i,j)) < Param ) then
             !
             NumberZeros = NumberZeros + 1
             !
          end if
          !
       end do
       !
       if ( NumberZeros == CopyMatrix.NColumns() ) then
          !
          res = .true.
          write(*,*) ' the row ', i, 'is zero'
          return
          !
       end if
    end do
    !
    !
  end function CheckZeroRowOrColumn



  subroutine Symmetrize( Matrix )
    !
    class(ClassMatrix), intent(inout) :: Matrix
    !
    integer :: i,j
    real(kind(1d0)) :: Average
    !
    if ( Matrix.NRows() /= Matrix.NColumns() ) then
       call Assert("The number of rows and the number of columns differ, in symmetrization.")
    end if
    !
    do j = 1, Matrix.NColumns()
       do i = j+1, Matrix.NRows()
          Average = 0.5d0 * ( &
               Matrix.Element(i,j) + &
               Matrix.Element(j,i) )
          call Matrix.SetElement(i,j,Average)
          call Matrix.SetElement(j,i,Average)
       end do
    end do
    !
  end subroutine Symmetrize



  subroutine ClassComplexMatrixSymmetrize( Matrix )
    !
    class(ClassComplexMatrix), intent(inout) :: Matrix
    !
    integer :: i,j
    complex(kind(1d0)) :: Average
    !
    if ( Matrix.NRows() /= Matrix.NColumns() ) then
       call Assert("The number of rows and the number of columns differ, in symmetrization.")
    end if
    !
    do j = 1, Matrix.NColumns()
       do i = j+1, Matrix.NRows()
          Average = 0.5d0 * ( &
               Matrix.Element(i,j) + &
               Matrix.Element(j,i) )
          call Matrix.SetElement(i,j,Average)
          call Matrix.SetElement(j,i,Average)
       end do
    end do
    !
  end subroutine ClassComplexMatrixSymmetrize



  subroutine ClassComplexMatrixAntiSymmetrize( Matrix )
    !
    class(ClassComplexMatrix), intent(inout) :: Matrix
    !
    integer :: i,j
    complex(kind(1d0)) :: Average
    !
    if ( Matrix.NRows() /= Matrix.NColumns() ) then
       call Assert("The number of rows and the number of columns differ, in symmetrization.")
    end if
    !
    do j = 1, Matrix.NColumns()
       do i = j+1, Matrix.NRows()
          Average = 0.5d0 * ( &
               Matrix.Element(i,j) - &
               Matrix.Element(j,i) )
          call Matrix.SetElement(i,j,Average)
          call Matrix.SetElement(j,i,-Average)
       end do
    end do
    !
  end subroutine ClassComplexMatrixAntiSymmetrize




  subroutine AntiSymmetrize( Matrix )
    !
    class(ClassMatrix), intent(inout) :: Matrix
    !
    integer :: i,j
    real(kind(1d0)) :: Average
    !
    if ( Matrix.NRows() /= Matrix.NColumns() ) then
       call Assert("The number of rows and the number of columns differ, in symmetrization.")
    end if
    !
    do j = 1, Matrix.NColumns()
       do i = 1, j-1
          Average = Matrix.Element(j,i)
          call Matrix.SetElement(i,j,-Average)
       end do
    end do
    !
  end subroutine AntiSymmetrize




  subroutine RemoveEpsilon( Matrix, Param )
    !
    class(ClassComplexMatrix),  intent(inout) :: Matrix
    real(kind(1d0)), optional, intent(in)    :: Param
    !
    integer :: i, j
    real(kind(1d0)), parameter :: DefaultParam = 1.d-14
    real(kind(1d0)) :: ActualPar
    complex(kind(1d0)), parameter :: z0 = (0.d0,0.d0)
    !
    if ( present(Param) ) then
       ActualPar = Param
    else
       ActualPar = DefaultParam
    end if
    !
    do j = 1, Matrix.NColumns()
       do i = 1, Matrix.NRows()
          !
          if ( abs(Matrix.Element(i,j)) < ActualPar ) then
             call Matrix.SetElement( i,j,z0 )
          end if
          !
       end do
    end do
    !
  end subroutine RemoveEpsilon


  subroutine SetComplexSpectralResolutionEigenValues( SpecRes, EigValVector )
    !
    Class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    complex(kind(1d0)),                    intent(in)    :: EigValVector(:)
    if ( allocated(SpecRes.EigenValues) ) deallocate( SpecRes.EigenValues )
    allocate( SpecRes.EigenValues, source = EigValVector )
  end subroutine SetComplexSpectralResolutionEigenValues




  subroutine SetComplexSpectralResolutionEigenVectors( SpecRes, EigVecMat, Identifier )
    class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    complex(kind(1d0)),                    intent(in)    :: EigVecMat(:,:)
    logical,   optional,                   intent(in)    :: Identifier
    if ( allocated(SpecRes.RightEigenVectors) ) deallocate( SpecRes.RightEigenVectors )
    allocate( SpecRes.RightEigenVectors, source = EigVecMat )
  end subroutine SetComplexSpectralResolutionEigenVectors



  subroutine SetComplexSpectralResolutionEigenVectorsMat( SpecRes, EigVecMat, Identifier )
    class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    class(ClassComplexMatrix)            , intent(in)    :: EigVecMat
    logical,   optional,                   intent(in)    :: Identifier
    if ( allocated(SpecRes.RightEigenVectors) ) deallocate( SpecRes.RightEigenVectors )
    call EigVecMat.FetchMatrix( SpecRes.RightEigenVectors )
  end subroutine SetComplexSpectralResolutionEigenVectorsMat


  !> Intended to remove the last N basis functions.
  subroutine ClassComplexSpectralResolutionReduceDimension( SpecRes, NumFunc )
    class(ClassComplexSpectralResolution), intent(inout) :: SpecRes
    integer                              , intent(in)    :: NumFunc
    complex(kind(1d0)), allocatable :: EigVec(:,:)
    SpecRes.Dim = SpecRes.Dim - NumFunc
    call SpecRes.Fetch( EigVec )
    deallocate( SpecRes.RightEigenVectors )
    allocate( SpecRes.RightEigenVectors, source = EigVec(1:SpecRes.Dim,:) )
    deallocate( EigVec )
  end subroutine ClassComplexSpectralResolutionReduceDimension


  !> Intended to remove the last N basis functions.
  subroutine ClassSpectralResolutionReduceDimension( SpecRes, NumFunc )
    class(ClassSpectralResolution), intent(inout) :: SpecRes
    integer                       , intent(in)    :: NumFunc
    real(kind(1d0)), allocatable :: EigVec(:,:)
    SpecRes.Dim = SpecRes.Dim - NumFunc
    call SpecRes.Fetch( EigVec )
    deallocate( SpecRes.EigenVectors )
    allocate( SpecRes.EigenVectors, source = EigVec(1:SpecRes.Dim,:) )
    deallocate( EigVec )
  end subroutine ClassSpectralResolutionReduceDimension




  subroutine ClassComplexMatrixAddRows( self, NNewRows, where_, After )
    class(ClassComplexMatrix), intent(inout) :: self
    integer           , intent(in)    :: NNewRows
    !> where is to be either "START" or "END"
    character(len=*), optional , intent(in) :: where_
    !> After can be any row in the input matrix.
    !! If specified, the function adds NNewRows between
    !! After and After+1.
    !! One can specify either where_ or After, but not both.
    integer         , optional , intent(in) :: After
    type(ClassComplexMatrix) :: mat
    if(NNewRows<=0)return
    if(present(where_))then

       if(present(After))call Assert("In AddRows, you can't specify both 'where_' and 'After'")

       if(where_.is."START")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows()+NNewRows, self.NColumns() )
             mat.A(NNewRows+1:,:)=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             !.. The absolute indexing of the matrix changes
             !   therefore it is better to make a full copy
             !   even if the allocated space wouldn't change
             call mat.initBanded(&
                  self.NRows()+NNewRows,&
                  self.NColumns(),&
                  self.LowerBandwidth()+NNewRows,&
                  self.UpperBandwidth()-NNewRows)
             mat.A=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select
       elseif(where_.is."END")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows()+NNewRows, self.NColumns() )
             mat.A(:self.NRows(),:)=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             !.. In this case both the absolute indexing and
             !   the storage size do not change
             !..
             self.NR=self.NR+NNewRows
          case DEFAULT 
             call Assert("Unrecognized matrix pattern")
          end select
       endif

    elseif( present(After) )then

       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows()+NNewRows, self.NColumns() )
          mat.A(1:After,:)=self.A(1:After,:)
          mat.A(After+1+NNewRows:,:)=self.A(After+1:,:)
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call Assert("'After' option not implemented in AddRows for banded matrices")
       case DEFAULT 
          call Assert("Unrecognized matrix pattern")
       end select

    else

       call Assert("In AddRows, either 'where_' or 'After' must be specified")

    endif

  end subroutine ClassComplexMatrixAddRows




  subroutine ClassComplexMatrixAddColumns( self, NNewCols, where_, After )
    class(ClassComplexMatrix), intent(inout) :: self
    integer           , intent(in) :: NNewCols
    character(len=*), optional, intent(in) :: where_
    integer, optional , intent(in) :: After
    type(ClassComplexMatrix) :: mat
    if(NNewCols<=0)return
    if(present(where_))then

       if(present(After))call Assert("In AddColumns, you can't specify both 'where_' and 'After'")

       if(where_.is."START")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows(), self.NColumns()+Nnewcols )
             mat.A(:,Nnewcols+1:)=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             call mat.initBanded(&
                  self.NRows(),&
                  self.NColumns()+Nnewcols,&
                  self.LowerBandwidth()-NNewCols,&
                  self.UpperBandwidth()+NNewCols)
             mat.A(:,Nnewcols+1:)=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select
       elseif(where_.is."END")then
          select case( self.Pattern )
          case( MATRIX_PATTERN_FULL )
             call mat.initFull( self.NRows(), self.NColumns()+Nnewcols )
             mat.A(:,:self.NColumns())=self.A
             self=mat
          case( MATRIX_PATTERN_BANDED )
             call mat.initBanded(&
                  self.NRows(),&
                  self.NColumns()+Nnewcols,&
                  self.LowerBandwidth(),&
                  self.UpperBandwidth())
             mat.A(:,:self.NColumns())=self.A
             self=mat 
          case DEFAULT 
             call Assert('Error: unrecognized matrix pattern')
          end select
       endif

    elseif( present(After) )then

       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows(), self.NColumns()+NNewCols )
          mat.A(:,1:After)=self.A(:,1:After)
          mat.A(:,After+1+NNewCols:)=self.A(:,After+1:)
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call Assert("'After' option not implemented in AddColumns for banded matrices")
       case DEFAULT 
          call Assert("Unrecognized matrix pattern")
       end select

    else

       call Assert("In AddColumns, either 'where_' or 'After' must be specified")

    endif

  end subroutine ClassComplexMatrixAddColumns





  subroutine ClassComplexMatrixRemoveRowsWhere( self, NKillRows, where_ )
    class(ClassComplexMatrix) , intent(inout) :: self
    integer                   , intent(in)    :: NKillRows
    character(len=*)          , intent(in)    :: where_

    character(len=*), parameter :: HERE="ClassMatrix::RemoveRows : "
    type(ClassComplexMatrix) :: mat
    integer :: iRow, iColumn
    complex(kind(1d0)), parameter :: Z0 = dcmplx(0.d0,0.d0)
    !
    if(NKillRows<=0) return
    if(NKillRows>self.NRows()) call Assert(HERE//"too many rows to be eliminated")
    !
    if(where_.is."START")then
       !
       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows()-NKillRows, self.NColumns() )
          do iColumn=1,self.NColumns()
             do iRow = 1, self.NRows() - NkillRows
                mat.A(iRow,iColumn)=self.A(NKillRows+iRow,iColumn)
             enddo
          enddo
          self=mat
       case( MATRIX_PATTERN_BANDED )
          !set to zero the elements in the rows eliminated
          do iRow = 1, NKillRows
             do iColumn=&
                  max(1,iRow-self.LowerBandwidth()),&
                  min(self.NColumns(),iRow+self.UpperBandwidth())
                call self.SetElement(iRow,iColumn,Z0)
             enddo
          enddo
          !.. The absolute indexing of the matrix changes
          !   therefore it is better to make a full copy
          !   even if the allocated space wouldn't change
          call mat.initBanded(&
               self.NRows()-NKillRows,&
               self.NColumns(),&
               self.LowerBandwidth()-NKillRows,&
               self.UpperBandwidth()+NKillRows)
          mat.A=self.A
          self=mat 
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
       !
    elseif(where_.is."END")then
       !
       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows()-NKillRows, self.NColumns() )
          do iColumn =1, self.NColumns()
             do iRow = 1, self.NRows()-NKillRows
                mat.A(iRow,iColumn)=self.A(iRow,iColumn)
             enddo
          enddo
          self=mat
       case( MATRIX_PATTERN_BANDED )
          !.. In this case both the absolute indexing and
          !   the storage size do not change
          !..
          !set to zero the elements in the rows eliminated
          do iRow = self.NRows()-NKillRows+1,self.NRows()
             do iColumn=&
                  max(1,iRow-self.LowerBandwidth()),&
                  min(self.NColumns(),iRow+self.UpperBandwidth())
                call self.SetElement(iRow,iColumn,Z0)
             enddo
          enddo
          self.NR=self.NR-NKillRows
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
       !
    endif
    !
  end subroutine ClassComplexMatrixRemoveRowsWhere




  subroutine ClassComplexMatrixRemoveRowsAfter( self, NKillRows, After )
    class(ClassComplexMatrix) , intent(inout) :: self
    integer                   , intent(in)    :: NKillRows
    integer                   , intent(in)    :: After

    character(len=*), parameter :: HERE="ClassMatrix::RemoveRows : "
    type(ClassComplexMatrix) :: mat
    integer :: iRow, iColumn
    complex(kind(1d0)), parameter :: Z0 = dcmplx(0.d0,0.d0)
    !
    if(NKillRows<=0) return
    if(NKillRows>self.NRows()) call Assert(HERE//"too many rows to be eliminated")
    !
    if(NKillRows>self.NRows()-After) call Assert(HERE//"too many rows to be eliminated")
    !
    select case( self.Pattern )
    case( MATRIX_PATTERN_FULL )
       call mat.initFull( self.NRows()-NKillRows, self.NColumns() )
       do iColumn=1,self.NColumns()
          do iRow=1,After
             mat.A(iRow,iColumn)=self.A(iRow,iColumn)
          enddo
          do iRow=After+1,self.NRows()-NKillRows
             mat.A(iRow,iColumn)=self.A(iRow+NKillRows,iColumn)
          enddo
       enddo
       self=mat
    case( MATRIX_PATTERN_BANDED )
       call Assert("'After' option not implemented in RemoveRows for banded matrices")
    case DEFAULT 
       call Assert("Unrecognized matrix pattern")
    end select
    !
  end subroutine ClassComplexMatrixRemoveRowsAfter






  subroutine ClassComplexMatrixRemoveColumnsWhere( self, Nkillcols, where_ )
    class(ClassComplexMatrix) , intent(inout) :: self
    integer                   , intent(in)    :: Nkillcols
    character(len=*)          , intent(in)    :: where_

    character(len=*), parameter :: HERE="ClassMatrix::RemoveColumns : "
    type(ClassComplexMatrix) :: mat
    integer :: i, j
    !
    if( Nkillcols <= 0 )return
    if( NKillCols > self.NColumns() ) call Assert(HERE//"too many Columns to be eliminated")
    !
    if(where_.is."START")then
       !
       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows(), self.NColumns()-Nkillcols )
          ! This line introduces core dumped segmetation error.
!!$             mat.A=self.A(:,Nkillcols+1:)
          do j = 1, mat.NColumns()
             !
             do i = 1, mat.NRows()
                !
                mat.A(i,j) = self.A(i,Nkillcols+j)
                !
             end do
             !
          end do
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call mat.initBanded(&
               self.NRows(),&
               self.NColumns()-Nkillcols,&
               self.LowerBandwidth()+NKillCols,&
               self.UpperBandwidth()-NKillCols)
          mat.A=self.A(:,Nkillcols+1:)
          self=mat 
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
       !
    elseif(where_.is."END")then
       
       select case( self.Pattern )
       case( MATRIX_PATTERN_FULL )
          call mat.initFull( self.NRows(), self.NColumns()-Nkillcols )
          ! This lines causes segmentation fault when many basis are present.
!!$             mat.A=self.A(:,:mat.NColumns())
          do j = 1, mat.NColumns()
             !
             do i = 1, mat.NRows()
                !
                mat.A(i,j) = self.A(i,j)
                !
             end do
             !
          end do
          self=mat
       case( MATRIX_PATTERN_BANDED )
          call mat.initBanded(&
               self.NRows(),&
               self.NColumns()-Nkillcols,&
               self.LowerBandwidth(),&
               self.UpperBandwidth())
          mat.A=self.A(:,:mat.NColumns())
          self=mat 
       case DEFAULT 
          call Assert('Error: unrecognized matrix pattern')
       end select
       !
    endif
    !
  end subroutine ClassComplexMatrixRemoveColumnsWhere





  subroutine ClassComplexMatrixRemoveColumnsAfter( self, Nkillcols, After )
    class(ClassComplexMatrix) , intent(inout) :: self
    integer                   , intent(in)    :: Nkillcols
    integer                   , intent(in)    :: After

    character(len=*), parameter :: HERE="ClassMatrix::RemoveColumns : "
    type(ClassComplexMatrix) :: mat
    integer :: i, j
    !
    if( Nkillcols <= 0 )return
    if( NKillCols > self.NColumns() ) call Assert(HERE//"too many Columns to be eliminated")
    !
    if( NKillCols > self.NColumns() - After ) call Assert(HERE//"too many Columns to be eliminated")
    !
    select case( self.Pattern )
    case( MATRIX_PATTERN_FULL )
       call mat.initFull( self.NRows(), self.NColumns()-NKillCols )
       mat.A(:,1:After)=self.A(:,1:After)
       ! This line causes errors when many basis functions are present.
!!$          mat.A(:,After+1:)=self.A(:,After+1+NKillCols:)
       do j = After+1, mat.NColumns()
          !
          do i = 1, mat.NRows()
             !
             mat.A(i,j) = self.A(i,j+NKillCols)
             !
          end do
          !
       end do
       self=mat
    case( MATRIX_PATTERN_BANDED )
       call Assert("'After' option not implemented in RemoveCols for banded matrices")
    case DEFAULT 
       call Assert(HERE//"Unrecognized matrix pattern")
    end select
    !
  end subroutine ClassComplexMatrixRemoveColumnsAfter



  subroutine ClassComplexMatrixRemoveAllZeroRows( self, Tolerance )
    !
    class(ClassComplexMatrix), intent(inout) :: self
    real(kind(1d0)), optional, intent(in)    :: Tolerance
    !
    integer :: i, j, Counter, NumBlocks
    real(kind(1d0)), parameter :: TolDefault = 1.d-10
    real(kind(1d0)) :: Tol
    logical :: RowIsZero
    integer, allocatable :: ZeroRowIndex(:), NewZeroRowIndex(:), NumRowBlock(:)
    !
    if ( present(Tolerance) ) then
       Tol = Tolerance
    else
       Tol = TolDefault
    end if
    !
    allocate( ZeroRowIndex(self.NRows()) )
    ZeroRowIndex = 0
    !
    Counter = 0
    !
    select case( self.Pattern )
    case( MATRIX_PATTERN_FULL )
       !
       do i = 1, self.NRows()
          !
          RowIsZero = .true.
          !
          do j = 1, self.NColumns()
             !
             if( abs(self.Element(i,j)) > Tol ) then
                RowIsZero = .false.
                exit
             end if
             !
          end do
          !
          if ( RowIsZero ) then
             Counter = Counter + 1
             ZeroRowIndex(Counter) = i 
          end if
          !
       end do
       !
       allocate( NewZeroRowIndex(Counter) )
       NewZeroRowIndex = 0
       allocate( NumRowBlock(Counter) )
       NumRowBlock = 1
       NumBlocks = 1
       NewZeroRowIndex(NumBlocks) = ZeroRowIndex(1)
       ! Determines the number of contiguous rows.
       do i = 1, Counter - 1
          !
          if ( (ZeroRowIndex(i+1) - ZeroRowIndex(i)) == 1 ) then
             NumRowBlock(NumBlocks) = NumRowBlock(NumBlocks) + 1
          else
             NumBlocks = NumBlocks + 1
             NewZeroRowIndex(NumBlocks) = ZeroRowIndex(i+1)
          end if
          !
       end do
       !
       do i = 1, NumBlocks
          write(*,*) "remove after row", NewZeroRowIndex(i)-1, "num of rows ", NumRowBlock(i), "tot rows to remove ", Counter
          call self.RemoveRows( NumRowBlock(i), NewZeroRowIndex(i)-1 )
          NewZeroRowIndex = NewZeroRowIndex - NumRowBlock(i)
       end do
       !
    case( MATRIX_PATTERN_BANDED )
       call Assert("Remove all zero rows not implemented for banded matrices")
    case DEFAULT 
       call Assert("Unrecognized matrix pattern")
    end select
    !
  end subroutine ClassComplexMatrixRemoveAllZeroRows



  subroutine ClassComplexMatrixRemoveAllZeroColumns( self, Tolerance )
    !
    class(ClassComplexMatrix), intent(inout) :: self
    real(kind(1d0)), optional, intent(in)    :: Tolerance
    !
    integer :: i, j, Counter, NumBlocks
    real(kind(1d0)), parameter :: TolDefault = 1.d-10
    real(kind(1d0)) :: Tol
    logical :: ColIsZero
    integer, allocatable :: ZeroColIndex(:), NewZeroColIndex(:), NumColBlock(:)
    !
    if ( present(Tolerance) ) then
       Tol = Tolerance
    else
       Tol = TolDefault
    end if
    !
    allocate( ZeroColIndex(self.NColumns()) )
    ZeroColIndex = 0
    !
    Counter = 0
    !
    select case( self.Pattern )
    case( MATRIX_PATTERN_FULL )
       !
       do j = 1, self.NColumns()
          !
          ColIsZero = .true.
          !
          do i = 1, self.NRows()
             !
             if( abs(self.Element(i,j)) > Tol ) then
                ColIsZero = .false.
                exit
             end if
             !
          end do
          !
          if ( ColIsZero ) then
             Counter = Counter + 1
             ZeroColIndex(Counter) = j
          end if
          !
       end do
       !
       allocate( NewZeroColIndex(Counter) )
       NewZeroColIndex = 0
       allocate( NumColBlock(Counter) )
       NumColBlock = 1
       NumBlocks = 1
       NewZeroColIndex(NumBlocks) = ZeroColIndex(1)
       ! Determines the number of contiguous rows.
       do i = 1, Counter - 1
          !
          if ( (ZeroColIndex(i+1) - ZeroColIndex(i)) == 1 ) then
             NumColBlock(NumBlocks) = NumColBlock(NumBlocks) + 1
          else
             NumBlocks = NumBlocks + 1
             NewZeroColIndex(NumBlocks) = ZeroColIndex(i+1)
          end if
          !
       end do
       !
       do j = 1, NumBlocks
write(*,*) "remove after col ", NewZeroColIndex(j)-1, "num of cols ", NumColBlock(j), "tot cols to remove ", Counter
          call self.RemoveColumns( NumColBlock(j), NewZeroColIndex(j)-1 )
          NewZeroColIndex = NewZeroColIndex - NumColBlock(j)
       end do
       !
    case( MATRIX_PATTERN_BANDED )
       call Assert("Remove all zero columns not implemented for banded matrices")
    case DEFAULT 
       call Assert("Unrecognized matrix pattern")
    end select
    !
  end subroutine ClassComplexMatrixRemoveAllZeroColumns




  !> Providing a complex matrix class, it is placed in the correct position inside a higher complex matrix class.
  subroutine ClassComplexMatrixComposeFromBlocks( self, RowBeg, RowEnd, ColBeg, ColEnd, BlockMat )
    !> Higher matrix.
    class(ClassComplexMatrix), intent(inout) :: self
    !> Row of the higher matrix at which the smaller one will begin.
    integer,                   intent(in)    :: RowBeg
    !> Row of the higher matrix at which the smaller one will end.
    integer,                   intent(in)    :: RowEnd
    !> Column of the higher matrix at which the smaller one will begin.
    integer,                   intent(in)    :: ColBeg
    !> Column of the higher matrix at which the smaller one will end.
    integer,                   intent(in)    :: ColEnd
    !> Smaller matrix.
    type(ClassComplexMatrix),  intent(in)    :: BlockMat
    !
    integer :: i, j
    !
    if ( .not. allocated(self.A) ) then
       call Assert( "The matrix has not been initialize, imposible to compose from blocks." )
    end if 
    !
    if ( BlockMat.NRows() /= (RowEnd-RowBeg+1) ) then
       call Assert( "The number of rows of the block does not match with the rows limits specified." )
    end if
    !
    if ( BlockMat.NColumns() /= (ColEnd-ColBeg+1) ) then
       call Assert( "The number of columns of the block does not match with the comlumns limits specified." )
    end if 
    !
    do j = ColBeg, ColEnd
       do i = RowBeg, RowEnd
          call self.SetElement( i, j, BlockMat.Element(i-RowBeg+1,j-ColBeg+1) )
       end do
    end do
    !
  end subroutine ClassComplexMatrixComposeFromBlocks



  subroutine ClassComplexMatrixBuildUpMatrix( Self, Blocks )
    class(ClassComplexMatrix), intent(inout) :: Self
    class(ClassComplexMatrix), intent(in)    :: Blocks(:,:)
    !
    integer :: RowBeg, RowEnd, ColBeg, ColEnd, TotNumColumns, TotNumRows
    integer :: i, j, NumBra, NumKet
    !
    call Self.Free()
    !
    NumBra = size(Blocks,1)
    NumKet = size(Blocks,2)
    !
    TotNumColumns = 0
    do j = 1, NumKet
       TotNumColumns = TotNumColumns + Blocks(1,j).NColumns()
    end do
    !
    TotNumRows = 0
    do i = 1, NumBra
       TotNumRows = TotNumRows + Blocks(i,1).NRows()
    end do
    ! 
    call Self.InitFull( TotNumRows, TotNumColumns )
    !
    RowBeg = 1
    ColBeg = 1
    !
    do j = 1, NumKet
       do  i = 1, NumBra 
          !
          RowEnd = RowBeg + Blocks(i,j).NRows() - 1
          ColEnd = ColBeg + Blocks(i,j).NColumns() - 1
          !
          call Self.ComposeFromBlocks( &
               RowBeg, &
               RowEnd, &
               ColBeg, &
               ColEnd, &
               Blocks(i,j) )
          !
          RowBeg = RowEnd + 1
          if ( i == NumBra ) then
             ColBeg = ColEnd + 1
             RowBeg = 1
          end if
          !
       end do
    end do
    !
  end subroutine ClassComplexMatrixBuildUpMatrix



  subroutine ClassMatrixBuildUpMatrix( Self, Blocks )
    class(ClassMatrix), intent(inout) :: Self
    class(ClassMatrix), intent(in)    :: Blocks(:,:)
    !
    integer :: RowBeg, RowEnd, ColBeg, ColEnd, TotNumColumns, TotNumRows
    integer :: i, j, NumBra, NumKet
    !
    call Self.Free()
    !
    NumBra = size(Blocks,1)
    NumKet = size(Blocks,2)
    !
    TotNumColumns = 0
    do j = 1, NumKet
       TotNumColumns = TotNumColumns + Blocks(1,j).NColumns()
    end do
    !
    TotNumRows = 0
    do i = 1, NumBra
       TotNumRows = TotNumRows + Blocks(i,1).NRows()
    end do
    ! 
    call Self.InitFull( TotNumRows, TotNumColumns )
    !
    RowBeg = 1
    ColBeg = 1
    !
    do j = 1, NumKet
       do  i = 1, NumBra 
          !
          RowEnd = RowBeg + Blocks(i,j).NRows() - 1
          ColEnd = ColBeg + Blocks(i,j).NColumns() - 1
          !
          call Self.ComposeFromBlocks( &
               RowBeg, &
               RowEnd, &
               ColBeg, &
               ColEnd, &
               Blocks(i,j) )
          !
          RowBeg = RowEnd + 1
          if ( i == NumBra ) then
             ColBeg = ColEnd + 1
             RowBeg = 1
          end if
          !
       end do
    end do
    !
  end subroutine ClassMatrixBuildUpMatrix



  !> Providing a matrix class, it is placed in the correct position inside a higher  matrix class.
  subroutine ClassMatrixComposeFromBlocks( self, RowBeg, RowEnd, ColBeg, ColEnd, BlockMat )
    !> Higher matrix.
    class(ClassMatrix), intent(inout) :: self
    !> Row of the higher matrix at which the smaller one will begin.
    integer,                   intent(in)    :: RowBeg
    !> Row of the higher matrix at which the smaller one will end.
    integer,                   intent(in)    :: RowEnd
    !> Column of the higher matrix at which the smaller one will begin.
    integer,                   intent(in)    :: ColBeg
    !> Column of the higher matrix at which the smaller one will end.
    integer,                   intent(in)    :: ColEnd
    !> Smaller matrix.
    type(ClassMatrix),  intent(in)    :: BlockMat
    !
    integer :: i, j
    !
    if ( .not. allocated(self.A) ) then
       call Assert( "The matrix has not been initialize, imposible to compose from blocks." )
    end if 
    !
    if ( BlockMat.NRows() /= (RowEnd-RowBeg+1) ) then
       call Assert( "The number of rows of the block does not match with the rows limits specified." )
    end if
    !
    if ( BlockMat.NColumns() /= (ColEnd-ColBeg+1) ) then
       call Assert( "The number of columns of the block does not match with the comlumns limits specified." )
    end if 
    !
    do j = ColBeg, ColEnd
       do i = RowBeg, RowEnd
          call self.SetElement( i, j, BlockMat.Element(i-RowBeg+1,j-ColBeg+1) )
       end do
    end do
    !
  end subroutine ClassMatrixComposeFromBlocks



  subroutine ClassComplexMatrixAllZeros( Mat, IsZero, Threshold )
    !
    class(ClassComplexMatrix), intent(in)  :: Mat
    logical,                   intent(out) :: IsZero
    real(kind(1d0)), optional, intent(in)  :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    IsZero = .false.
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    do j = 1, Mat.NColumns()
       do i = 1, Mat.NRows()
          if ( abs(Mat.Element(i,j)) > Thr ) then
             return
          end if
       end do
    end do
    !
    IsZero = .true.
    !
  end subroutine ClassComplexMatrixAllZeros




  logical function ClassComplexMatrixIsHermitian( Mat,Threshold ) result( IsHerm )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    integer :: i, j
    real(kind(1d0)) :: Val, Thr
    !
    if ( Mat.NRows() /= Mat.NColumns() ) call Assert( &
         'To check the hermiticity of a matrix this has to be squared.' )
    !
    IsHerm = .false.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    do j = 1, Mat.NColumns()
       do i = 1, j
          Val = abs( Mat.Element(i,j)-conjg(Mat.Element(j,i)) )
          if ( Val > Thr ) then
!!$             write(OUTPUT_UNIT,*) "Row index ", i, "Columns index ", j
!!$             write(OUTPUT_UNIT,*) "Elem i,j ", Mat.Element(i,j), "Elem j, i ", Mat.Element(j,i), "Abs Difference", Val, "Tolerance ", Thr
             return
          end if
       end do
    end do
    !
    IsHerm = .true.
    !
  end function ClassComplexMatrixIsHermitian



  logical function ClassComplexMatrixIsSymmetric( Mat, Threshold ) result( IsSym )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    integer :: i, j
    real(kind(1d0)) :: Val, Thr
    !
    IsSym = .false.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    do j = 1, Mat.NColumns()
       do i = 1, j-1
          Val = abs( Mat.Element(i,j)-Mat.Element(j,i) )
          if ( Val > Thr ) then
             write(OUTPUT_UNIT,*) "Row index ", i, "Columns index ", j
             write(OUTPUT_UNIT,*) "Elem i,j ", Mat.Element(i,j), "Elem j, i ", Mat.Element(j,i), "Abs Difference", Val, "Tolerance ", Thr
             return
          end if
       end do
    end do
    !
    IsSym = .true.
    !
  end function ClassComplexMatrixIsSymmetric



  logical function ClassMatrixIsSymmetric( Mat, Threshold ) result( IsSym )
    !
    class(ClassMatrix),        intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    integer :: i, j
    real(kind(1d0)) :: Val, Thr
    !
    IsSym = .false.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    do j = 1, Mat.NColumns()
       do i = 1, j-1
          Val = abs( Mat.Element(i,j)-Mat.Element(j,i) )
          if ( Val > Thr ) then
             write(OUTPUT_UNIT,*) "Row index ", i, "Columns index ", j
             write(OUTPUT_UNIT,*) "Elem i,j ", Mat.Element(i,j), "Elem j, i ", Mat.Element(j,i), "Abs Difference", Val, "Tolerance ", Thr
             return
          end if
       end do
    end do
    !
    IsSym = .true.
    !
  end function ClassMatrixIsSymmetric



  logical function ClassComplexMatrixIsAntiSymmetric( Mat, Threshold ) result( IsASym )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    integer :: i, j
    real(kind(1d0)) :: Val, Thr
    !
    IsASym = .false.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    do j = 1, Mat.NColumns()
       do i = 1, j-1
          Val = abs( Mat.Element(i,j)+Mat.Element(j,i) )
          if ( Val > Thr ) then
             write(OUTPUT_UNIT,*) "Row index ", i, "Columns index ", j
             write(OUTPUT_UNIT,*) "Elem i,j ", Mat.Element(i,j), "Elem j, i ", Mat.Element(j,i), "Abs Difference", Val, "Tolerance ", Thr
             return
          end if
       end do
    end do
    !
    IsASym = .true.
    !
  end function ClassComplexMatrixIsAntiSymmetric



  logical function ClassMatrixIsAntiSymmetric( Mat, Threshold ) result( IsASym )
    !
    class(ClassMatrix),        intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    integer :: i, j
    real(kind(1d0)) :: Val, Thr
    !
    IsASym = .false.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    do j = 1, Mat.NColumns()
       do i = 1, j-1
          Val = abs( Mat.Element(i,j)+Mat.Element(j,i) )
          if ( Val > Thr ) then
             write(OUTPUT_UNIT,*) "Row index ", i, "Columns index ", j
             write(OUTPUT_UNIT,*) "Elem i,j ", Mat.Element(i,j), "Elem j, i ", Mat.Element(j,i), "Abs Difference", Val, "Tolerance ", Thr
             return
          end if
       end do
    end do
    !
    IsASym = .true.
    !
  end function ClassMatrixIsAntiSymmetric




  subroutine ClassMatrixInverse( Mat, ResMat )
    !
    class(ClassMatrix), intent(inout)  :: Mat
    class(ClassMatrix), intent(out) :: ResMat
    !
    type(ClassMatrix) :: AuxMat
    real(kind(1d0)), allocatable :: AuxArray(:,:), work(:)
    integer :: n, info
    integer(kind=8) :: lwork
    integer, allocatable :: ipiv(:)
    !
    if ( Mat.IsBanded() ) then
       call Mat.ConvertToSquared( AuxMat )
    elseif( Mat.IsFull() ) then
       AuxMat = Mat
    end if
    !
    if ( AuxMat.NRows() /= AuxMat.NColumns() ) then
       call Assert( "Impossible to compute the matrix inverse because is rectangular." )
    end if
    !
    call AuxMat.FetchMatrix( AuxArray )
    n = AuxMat.NRows()
    allocate( ipiv(n) )
    call AuxMat.Free()
    ! 
    call dgetrf( n, n, AuxArray, n, ipiv, info )
    !
    if ( info /= 0 ) then
       write(output_unit,*) "The info parameter in LU factorization is no zero, info =", info
       call Assert( "Program aborted due to deficient LU factorization." )
    end if
    !
    allocate( work(1) )
    call dgetri( n, AuxArray, n, ipiv, work, -1, info )
    lwork = max(int(work(1)),n)
    deallocate( work )
    allocate( work(lwork) )
    !
    call dgetri( n, AuxArray, n, ipiv, work, lwork, info )
    !
    ResMat = AuxArray
    !
    deallocate( AuxArray )
    !
  end subroutine ClassMatrixInverse




  subroutine ClassComplexMatrixInverse( Mat, ResMat )
    !
    class(ClassComplexMatrix), intent(in)  :: Mat
    class(ClassComplexMatrix), intent(out) :: ResMat
    !
    type(ClassComplexMatrix) :: AuxMat
    complex(kind(1d0)), allocatable :: AuxArray(:,:), work(:)
    integer :: n, info
    integer(kind=8) :: lwork
    integer, allocatable :: ipiv(:)
    !
    if ( Mat.IsBanded() ) then
       call Mat.ConvertToFull( AuxMat )
    elseif( Mat.IsFull() ) then
       AuxMat = Mat
    end if
    !
    if ( AuxMat.NRows() /= AuxMat.NColumns() ) then
       call Assert( "Impossible to compute the matrix inverse because is rectangular." )
    end if
    ! 
    call AuxMat.FetchMatrix( AuxArray )
    n = AuxMat.NRows()
    allocate( ipiv(n) )
    call AuxMat.Free()
    !  
    call zgetrf( n, n, AuxArray, n, ipiv, info )
    !
    if ( info /= 0 ) then
       write(output_unit,*) "The info parameter in LU factorization is not zero, info =", info
       call Assert( "Program aborted due to deficient LU factorization." )
    end if
    !
    lwork = -1
    allocate(work(1))
    !
    call zgetri( n, AuxArray, n, ipiv, work, lwork, info )
    !
    lwork = max(int(work(1)),n)
    deallocate( work )
    allocate( work(lwork) )
    !
    call zgetri( n, AuxArray, n, ipiv, work, lwork, info )
    !
    if ( info /= 0 ) then
       write(output_unit,*) "The info parameter in the invertion is not zero, info =", info
       call Assert( "Program aborted due to deficient invertion." )
    end if
    ResMat = AuxArray
    !
    deallocate( AuxArray )
    !
  end subroutine ClassComplexMatrixInverse




  logical function ClassComplexMatrixIsIdentity( Mat, Threshold ) result( Id )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    Id = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is the Identity Matrix the format must be Full." )
    end if
    !
    if ( Mat.NRows() /= Mat.NColumns() ) then
       call Assert( "To determine if the matrix is the Identity it must be squared." )
    end if
    !
    !
    !..Check the diagonal
    do i = 1, Mat.NRows()
       if ( abs(dble(Mat.Element(i,i)) - 1.d0) > Thr ) then
          Id = .false.
          write(output_unit,*) "Re(Element)", i, i, dble(Mat.Element(i,i))
          return
       elseif ( abs(aimag(Mat.Element(i,i)) - 0.d0) > Thr ) then
          Id = .false.
          write(output_unit,*) "Imag(Element)", i, i, aimag(Mat.Element(i,i))
          return
       end if
    end do
    !
    !.. Check off-diafonal
    do j = 1, Mat.NRows()
       do i = 1, Mat.NRows()
          if ( i /= j ) then
             if ( abs(dble(Mat.Element(i,j)) - 0.d0) > Thr ) then
                Id = .false.
                write(output_unit,*) "Re(Element)", i, j, dble(Mat.Element(i,j))
                return
             elseif ( abs(aimag(Mat.Element(i,j)) - 0.d0) > Thr ) then
                Id = .false.
                write(output_unit,*) "Imag(Element)", i, j, aimag(Mat.Element(i,j))
                return
             end if
          end if
       end do
    end do
    !
  end function ClassComplexMatrixIsIdentity



  logical function ClassComplexMatrixMatrixIsDiagonal( Mat, Threshold ) result( Diag )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    Diag = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is a diagonal Matrix the format must be Full." )
    end if
    !
    if ( Mat.NRows() /= Mat.NColumns() ) then
       call Assert( "To determine if the matrix is a diagonal Identity it must be squared." )
    end if
    !
    !.. Check off-diafonal
    do j = 1, Mat.NRows()
       do i = 1, Mat.NRows()
          if ( i /= j ) then
             if ( abs(dble(Mat.Element(i,j)) - 0.d0) > Thr ) then
                Diag = .false.
                write(output_unit,*) "Re(Element)", i, j, dble(Mat.Element(i,j))
                return
             elseif ( abs(aimag(Mat.Element(i,j)) - 0.d0) > Thr ) then
                Diag = .false.
                write(output_unit,*) "Imag(Element)", i, j, aimag(Mat.Element(i,j))
                return
             end if
          end if
       end do
    end do
    !
  end function ClassComplexMatrixMatrixIsDiagonal




  logical function ClassComplexMatrixIsUnitary( Mat, Threshold ) result( Unitary )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    type(ClassComplexMatrix) :: AuxMat
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is unitary the Matrix format must be Full." )
    end if
    !
    if ( Mat.NRows() /= Mat.NColumns() ) then
       call Assert( "To determine if the matrix is unitary it must be squared." )
    end if
    !
    call Mat.TransposeConjugate( AuxMat )
    call AuxMat.Multiply( Mat, 'Right', 'N' )
    Unitary = AuxMat.IsIdentity( Thr )
    !
  end function ClassComplexMatrixIsUnitary




  logical function ClassComplexMatrixIsReal( Mat, Threshold ) result( Id )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    Id = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is real the format must be Full." )
    end if
    !
    !
    do j = 1, Mat.NColumns()
       do i = 1, Mat.NRows()
          if ( abs(aimag(Mat.Element(i,j))) > Thr ) then
             Id = .false.
             write(output_unit,*) "Imag(Element)", i, j, aimag(Mat.Element(i,j))
             return
          end if
       end do
    end do
    !
  end function ClassComplexMatrixIsReal




  logical function ClassMatrixIsIdentity( Mat, Threshold ) result( Id )
    !
    class(ClassMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    Id = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is the Identity Matrix the format must be Full." )
    end if
    !
    if ( Mat.NRows() /= Mat.NColumns() ) then
       call Assert( "To determine if the matrix is the Identity it must be squared." )
    end if
    !
    !
    !..Check the diagonal
    do i = 1, Mat.NRows()
       if ( abs(Mat.Element(i,i) - 1.d0) > Thr ) then
          Id = .false.
          write(OUTPUT_UNIT,*) "Element", i, i, Mat.Element(i,i)
          return
       end if
    end do
    !
    !.. Check off-diafonal
    do j = 1, Mat.NRows()
       do i = 1, Mat.NRows()
          if ( i /= j ) then
             if ( abs(Mat.Element(i,j) - 0.d0) > Thr ) then
                Id = .false.
                write(OUTPUT_UNIT,*) "Element", i, j, Mat.Element(i,j)
                return
             end if
          end if
       end do
    end do
    !
  end function ClassMatrixIsIdentity



  logical function ClassMatrixMatrixIsDiagonal( Mat, Threshold ) result( Diag )
    !
    class(ClassMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    Diag = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is a diagonal Matrix the format must be Full." )
    end if
    !
    if ( Mat.NRows() /= Mat.NColumns() ) then
       call Assert( "To determine if the matrix is diagonal it must be squared." )
    end if
    !
    !.. Check off-diafonal
    do j = 1, Mat.NRows()
       do i = 1, Mat.NRows()
          if ( i /= j ) then
             if ( abs(Mat.Element(i,j) - 0.d0) > Thr ) then
                Diag = .false.
                write(OUTPUT_UNIT,*) "Element", i, j, Mat.Element(i,j)
                return
             end if
          end if
       end do
    end do
    !
  end function ClassMatrixMatrixIsDiagonal




  logical function ClassMatrixIsZero( Mat, Threshold ) result( Zero )
    !
    class(ClassMatrix),        intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr
    !
    Zero = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is the Identity Matrix the format must be Full." )
    end if
    !
    !
    do j = 1, Mat.NColumns()
       do i = 1, Mat.NRows()
          if ( abs(Mat.Element(i,j) - 0.d0) > Thr ) then
             Zero = .false.
             write(output_unit,*) "Element", i, j, Mat.Element(i,j)
             return
          end if
       end do
    end do
    !
  end function ClassMatrixIsZero




  logical function ClassComplexMatrixIsZero( Mat, Threshold ) result( Zero )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr, Val
    !
    Zero = .true.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is the Identity Matrix the format must be Full." )
    end if
    !
    !
    do j = 1, Mat.NColumns()
       do i = 1, Mat.NRows()
          if ( abs(dble(Mat.Element(i,j)) - 0.d0) > Thr ) then
             Zero = .false.
             Val = dble(Mat.Element(i,j))
             write(*,*) "Re(Element)", i, j, Val
             return
          elseif ( abs(aimag(Mat.Element(i,j)) - 0.d0) > Thr ) then
             Zero = .false.
             Val = aimag(Mat.Element(i,j))
             write(*,*) "Imag(Element)", i, j, Val
             return
          end if
       end do
    end do
    !
  end function ClassComplexMatrixIsZero



  logical function ClassComplexMatrixIsComplex( Mat, Threshold ) result( Comp )
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    real(kind(1d0)), optional, intent(in) :: Threshold
    !
    integer :: i, j
    real(kind(1d0)) :: Thr, Val
    !
    Comp = .false.
    !
    Thr = tiny(1d0)
    if ( present(Threshold) ) Thr = Threshold
    !
    if ( .not.Mat.IsFull() ) then
       call Assert( "To test if the matrix is the Identity Matrix the format must be Full." )
    end if
    !
    !
    do j = 1, Mat.NColumns()
       do i = 1, Mat.NRows()
          if ( abs(aimag(Mat.Element(i,j))) > Thr ) then
             Comp = .true.
             Val = aimag(Mat.Element(i,j))
             write(*,*) "Imag(Element)", i, j, Val
             return
          end if
       end do
    end do
    !
  end function ClassComplexMatrixIsComplex




  subroutine ClassMatrixFindZeroLines( Mat, Indexes, Identifier, Tolerance )
    !
    class(ClassMatrix), intent(in) :: Mat
    integer, allocatable, intent(out) :: Indexes(:)
    character(len=*),     intent(in) :: Identifier
    real(kind(1d0)), optional, intent(in) :: Tolerance
    !
    real(kind(1d0)) :: Tol
    integer :: i, j, Counter
    integer, allocatable :: AuxInd(:)
    logical :: LineIsZero
    !
    Tol = tiny(1d0)
    if ( present(Tolerance) ) Tol = Tolerance
    !
    if ( Identifier .is. "row" ) then
       !
       allocate( AuxInd(Mat.NRows()) )
       AuxInd = 0
       !
       Counter = 0
       do i = 1, Mat.NRows()
          !
          LineIsZero = .true.
          !
          do j = 1, Mat.NColumns()
             if ( abs(Mat.Element(i,j)) > Tol ) then
                LineIsZero = .false.
                exit
             end if
          end do
          !
          if ( LineIsZero ) then
             Counter = Counter + 1
             AuxInd(Counter) = i
          end if
          !
       end do
       !
    elseif ( Identifier .is. "column" ) then
       !
       allocate( AuxInd(Mat.NColumns()) )
       AuxInd = 0
       !
       Counter = 0
       do j = 1, Mat.NColumns()
          !
          LineIsZero = .true.
          !
          do i = 1, Mat.NRows()
             if ( abs(Mat.Element(i,j)) > Tol ) then
                LineIsZero = .false.
                exit
             end if
          end do
          !
          if ( LineIsZero ) then
             Counter = Counter + 1
             AuxInd(Counter) = j
          end if
          !
       end do
       !
    else
       !
       call Assert( "The identifier for finding the lines that are zero must be 'row' or 'column', aborting program." )
       !
    end if
    !
    !
    if ( Counter > 0 ) then
       allocate( Indexes, source = AuxInd(1:Counter) )
    else
       allocate( Indexes(1) )
       Indexes(1) = 0
!!$       write(output_unit,*) "There are no zero "//trim(adjustl(Identifier))
    end if
    !
    deallocate( AuxInd )
    !
  end subroutine ClassMatrixFindZeroLines





  subroutine ClassComplexMatrixSolveHomogScattExpert( Mat, ArrayLU )
    !
    !> Rectangular scattering matrix H-E*S, NC = NR + NumSol.
    class(ClassComplexMatrix),       intent(in)  :: Mat
    complex(kind(1d0)), allocatable, intent(out) :: ArrayLU(:,:)
    !
    !> Squared scattering matrix
    type(ClassComplexMatrix) :: SqrMat
    complex(kind(1d0)), allocatable :: a(:,:), af(:,:), b(:,:), work(:), x(:,:)
    real(kind(1d0)), allocatable :: r(:), c(:), rwork(:), ferr(:), berr(:)
    real(kind(1d0)) :: rcond
    integer :: NumSol, info, n, nrhs, lda, ldaf, ldb, ldx
    integer, allocatable :: ipiv(:)
    character(len=:), allocatable :: fact, trans, equed
    !
    complex(kind(1d0)), parameter :: Z0 = dcmplx(0.d0,0.d0)
    complex(kind(1d0)), parameter :: Z1 = dcmplx(1.d0,0.d0)
    complex(kind(1d0)), parameter :: Zi = dcmplx(0.d0,1.d0)
    !
    NumSol = Mat.NColumns()-Mat.NRows()
    allocate( b(Mat.NColumns(),Mat.NColumns()) )
!!$    b = Z0
    b = Z1*1.d-50
    !
    SqrMat = Mat
    call SqrMat.AddRows( NumSol, "END" )
    ! 
    !.. Routine arguments
    allocate( fact,  source = "N" )
    allocate( trans, source = "N" )
    allocate( equed, source = "I" ) ! only to reserve memory
    n = Mat.NColumns()
    nrhs = NumSol
    call SqrMat.FetchMatrix( a )
    allocate( work(2*n) )
    lda  = n
    ldaf = n
    ldb  = n
    ldx  = n
    allocate( af(ldaf,ldaf) )
    allocate( rwork(2*n) ) 
    allocate( ipiv(n) )
    allocate( r(n), c(n) )
    allocate( x(ldx,nrhs) )
    allocate( ferr(nrhs) )
    allocate( berr(nrhs) )
    !
    !
    call zgesvx( fact, trans, n, nrhs, a, lda, &
         af, ldaf, ipiv, equed, r, c, b, ldb, &
         x, ldx, rcond, ferr, berr, work, rwork, info )  
    !
    write(output_unit,*) "nrhs = ", nrhs
    write(output_unit,*) "info = ", info
    write(output_unit,*) "rcond = ", rcond
    write(output_unit,*) "ferr: ", ferr
    write(output_unit,*) "berr: ", berr
    write(output_unit,*) "equed = ", equed
    write(output_unit,*) "rwork(1) = ", rwork(1)
    !
    allocate( ArrayLU, source = a )
    !
  end subroutine ClassComplexMatrixSolveHomogScattExpert



  logical function ClassMatrixIsInitialized( Mat ) result(IsInit)
    !
    class(ClassMatrix), intent(in) :: Mat
    !
    if ( allocated(Mat.A) ) then
       IsInit = .true.
    else
       IsInit = .false.
    end if
    !
  end function ClassMatrixIsInitialized



  logical function ClassComplexMatrixIsInitialized( Mat ) result(IsInit)
    !
    class(ClassComplexMatrix), intent(in) :: Mat
    !
    if ( allocated(Mat.A) ) then
       IsInit = .true.
    else
       IsInit = .false.
    end if
    !
  end function ClassComplexMatrixIsInitialized


  
  !..Rescale the right eigenvectors CR to have the proper normalization:
  ! CL^{\dagga} * CR = 1. 
  subroutine ReScaleRightEigVectors( CR )
    !
    complex(kind(1d0)), intent(inout) :: CR(:,:)
    !
    integer :: j
    complex(kind(1d0)) :: zw1
    type(ClassComplexMatrix) :: Aux, Aux2
    complex(kind(1d0)), parameter :: Z1=dcmplx(1.d0,0.d0)
    !
    do j = 1, size(CR,2)
       zw1 = Z1/sqrt(sum(CR(:,j)**2))
       CR(:,j) = zw1 * CR(:,j)
    end do
    !
    !
    Aux = CR
    call Aux.Transpose( Aux2 )
    call Aux2.Multiply( Aux, "Right", "N" )
    if ( .not.Aux2.IsIdentity( COMPUTATION_THRESHOLD ) ) call ErrorMessage( &
         "Right eigenvectors are not well normalized" )
    !
  end subroutine ReScaleRightEigVectors



  !..Rescale the eigenvectors CR to have the normalization:
  ! CR^{\dagga} * CR = 1. 
  subroutine NormalizedAsHermitian( CR )
    !
    complex(kind(1d0)), intent(inout) :: CR(:,:)
    !
    integer :: j
    complex(kind(1d0)) :: zw1
    type(ClassComplexMatrix) :: Aux, Aux2
    complex(kind(1d0)), parameter :: Z1=dcmplx(1.d0,0.d0)
    !
    do j = 1, size(CR,2)
       zw1 = Z1/sqrt(sum(CR(:,j)*conjg(CR(:,j))))
       CR(:,j) = zw1 * CR(:,j)
    end do
    !
    !
    Aux = CR
    call Aux.TransposeConjugate( Aux2 )
    call Aux2.Multiply( Aux, "Right", "N" )
    if ( .not.Aux2.IsIdentity( COMPUTATION_THRESHOLD ) ) call ErrorMessage( &
         "The hermitian operator's eigenvectors are not well normalized" )
    !
  end subroutine NormalizedAsHermitian



  !> Gets the maximum eigenvalue of a matrix.
  real(kind(1d0)) function ClassMatrixGetMaxEigVal( Mat ) result(EigVal)
    !
    class(ClassMatrix), intent(in) :: Mat
    type(ClassSpectralResolution) :: SpecRes
    real(kind(1d0)), allocatable :: EigValVec(:)
    !
    call Mat.Diagonalize( SpecRes )
    call SpecRes.Fetch( EigValVec )
    call SpecRes.Free()
    EigVal = maxval(EigValVec) 
    !
  end function ClassMatrixGetMaxEigVal


  !> Gets the minimum eigenvalue of a matrix.
  real(kind(1d0)) function ClassMatrixGetMinEigVal( Mat ) result(EigVal)
    !
    class(ClassMatrix), intent(in) :: Mat
    type(ClassSpectralResolution) :: SpecRes
    real(kind(1d0)), allocatable :: EigValVec(:)
    !
    call Mat.Diagonalize( SpecRes )
    call SpecRes.Fetch( EigValVec )
    call SpecRes.Free()
    EigVal = minval(EigValVec) 
    !
  end function ClassMatrixGetMinEigVal



end module ModuleMatrix
