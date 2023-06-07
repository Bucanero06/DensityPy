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
module ModuleMatrix

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString

  implicit none
  private

  real(kind(1d0)) , parameter :: THRESHOLD_MATRIX_WRITE = 1.d-20
  integer         , parameter :: DEFAULT_OUTPUT_UNIT    = OUTPUT_UNIT
  integer         , parameter :: DEFAULT_INPUT_UNIT     = INPUT_UNIT
  integer         , parameter :: FILE_NAME_LENGTH       = 512
  real(kind(1d0)) , parameter :: COMPUTATION_THRESHOLD  = 1.d-12

  enum, bind(c) 
     enumerator :: MATRIX_PATTERN_FULL
     enumerator :: MATRIX_PATTERN_DBANDED
     enumerator :: MATRIX_PATTERN_SUBMATRIX
     enumerator :: MATRIX_PATTERN_MIXED
     enumerator :: MATRIX_PATTERN_UNDEFINED
     enumerator :: MATRIX_PATTERN_DIAGONAL
  end enum

  character(len=*), parameter :: DOUBLE_PRINT_FORMAT  = "(d24.16)"

  !> Defines a general DoublePrecision matrix
  type, public :: ClassMatrix
     private
     !.. Pattern : Filled pattern
     !.. NR      : Number of Rows
     !.. NC      : Number of Columns
     !.. NU      : Number of Upper diagonals in dbanded matrix.
     !.. NL      : Number of Lower diagonals in dbanded matrix
     !> Pattern which will be used to stored the general matrix, could be 
     !! - Full storage (same array as the original matrix) or 
     !! - Dbanded storage (takes advantage of dbanded matrices saving memory 
     !!   space, that would not be possible using the Full storage 
     !!   representation of the original matrix), or 
     !! - Diagonal storage (takes advantage of dagonal matrices)%  
     integer :: Pattern = MATRIX_PATTERN_UNDEFINED
     integer :: NR = 0
     integer :: NC = 0
     integer :: NU = 0
     integer :: NL = 0
     !*** MUST INTEGRATE THE NON-ZERO DIMENSIONS IN THE REST OF THE CODE
     integer :: NNZC = 0
     integer :: NNZR = 0
     !.. For SUBMATRIX, [R|C]Mask(i)=.TRUE. if the row|columned is stored
     logical, allocatable :: RMAsk(:)
     logical, allocatable :: CMAsk(:)
     !> Minimum matrix row index
     integer :: NRMin = 0
     !> Maximum matrix row index
     integer :: NRMax = 0
     !> Minimum matrix column index
     integer :: NCMin = 0
     !> Maximum matrix column index
     integer :: NCMax = 0
     !> Allocated size (potentially non-zero)
     integer :: Size  = 0
     !> The matrix in the requested representation( Full, Dbanded or Diagonal)
     Real(kind(1d0)), public, allocatable :: A(:,:) 
   contains
     procedure, public :: IsInitialized  => Classmatrix_IsInitialized
     procedure, public :: InitFull       => Classmatrix_InitFull
     procedure, public :: InitDbanded    => Classmatrix_InitDbanded
     procedure, public :: InitSubmatrix  => Classmatrix_InitSubmatrix
     procedure, public :: Free           => Classmatrix_Free
     procedure, public :: IsFull         => Classmatrix_IsFull
     procedure, public :: IsDbanded      => Classmatrix_IsDbanded
     procedure, public :: IsDiagonal     => Classmatrix_IsDiagonal
     procedure, public :: LowerBandWidth => Classmatrix_LowerBandwidth
     procedure, public :: UpperBandWidth => Classmatrix_UpperBandwidth
     procedure, public :: Norm1          => Classmatrix_Norm1
     procedure, public :: Norm2          => Classmatrix_Norm2
     procedure, public :: NRows          => Classmatrix_NRows
     procedure, public :: NColumns       => Classmatrix_NColumns
     procedure, public :: Element        => Classmatrix_Element
     procedure, public :: SetElement     => Classmatrix_SetElement
     procedure, public :: FetchMatrix    => Classmatrix_FetchMatrix
     procedure, public :: FetchColumn    => Classmatrix_FetchColumn
     procedure, public :: Drop           => ClassMatrix_Drop
     procedure, public :: GetSubMatrix   => ClassMatrix_GetSubmatrix
     procedure, public :: Add            => ClassMatrix_AddClassMatrix
     procedure, public :: SetDiagonal    => ClassMatrix_SetDiagonal
     procedure, public :: AddDiagonal    => ClassMatrix_AddDiagonal
     procedure, public :: Plug           => Classmatrix_Plug
     procedure, public :: BuildUpMatrix  => Classmatrix_BuildUpMatrix
     procedure, public :: Transpose      => Classmatrix_Transpose
     generic  , public :: Write          => WriteClassMatrixToUnit, WriteClassMatrixToFile
     generic  , public :: Read           => ReadClassMatrixFromUnit, ReadClassMatrixFromFile
     procedure, public :: Contract4D12   => ClassMatrix_x_ClassMatrix4D12 
     procedure, public :: Contract4D13   => ClassMatrix_x_ClassMatrix4D13 
     generic  , public :: Diagonalize    => &
          ClassMatrix_Diagonalize, &
          ClassMatrix_Diagonalize_BANDED
     generic  , public :: Set      => &
          Classmatrix_SetMatrix_    , &
          Classmatrix_SetSubMatrix  , &
          Classmatrix_SetSubMatrixCM, &
          Classmatrix_SetVector     , &
          Classmatrix_SetElement 
     generic  , public :: Multiply => &
          Double_x_ClassMatrix      , &
          ClassMatrix_x_ClassMatrix , &
          ClassMatrix_x_Matrix      , &
          ClassMatrix_x_vec
     generic :: assignment(=)      => &
          AssignDoubleToClassMatrix , &
          AssignIntegerToClassMatrix, &
          CopyClassMatrixToClassMatrix, &
          CopyMatrixToClassMatrix   , &
          CopyVectorToClassMatrix

     procedure, private :: Classmatrix_SetElement
     procedure, private :: Classmatrix_SetMatrix_
     procedure, private :: Classmatrix_SetSubMatrix
     procedure, private :: Classmatrix_SetSubMatrixCM
     procedure, private :: Classmatrix_SetVector
     procedure, private :: ReadClassMatrixFromUnit
     procedure, private :: ReadClassMatrixFromFile
     procedure, private :: WriteClassMatrixToUnit
     procedure, private :: WriteClassMatrixToFile
     procedure, private :: AssignDoubleToClassMatrix
     procedure, private :: AssignIntegerToClassMatrix
     procedure, private :: CopyClassMatrixToClassMatrix
     procedure, private :: CopyMatrixToClassMatrix
     procedure, private :: CopyVectorToClassMatrix
     procedure, private :: Double_x_ClassMatrix
     procedure, private :: ClassMatrix_x_ClassMatrix
     procedure, private :: ClassMatrix_x_Vec
     procedure, private :: ClassMatrix_x_Matrix
     procedure, private :: ClassMatrix_Diagonalize
     procedure, private :: ClassMatrix_Diagonalize_Banded

     !> Deallocates all ClassMatrix atributes.
     final :: Classmatrix_Finalize

  end type ClassMatrix
  
  public :: ClassMatrixFullContraction


  !> Defines a general DoublePrecision matrix
  type, public :: ClassMatrix4D
     private
     Real(kind(1d0)), public, allocatable :: A(:,:,:,:) 
   contains
     private
     procedure, public :: IsInitialized  => ClassMatrix4DIsInitialized
     procedure, public :: Init           => ClassMatrix4DInit
     procedure, public :: Free           => ClassMatrix4DFree
     procedure, public :: Size           => ClassMatrix4DSize
     procedure, public :: TotalSize      => ClassMatrix4DTotalSize
     procedure, public :: Get            => ClassMatrix4DGetValue
     procedure, public :: Write          => WriteClassMatrix4DToUnit
     procedure, public :: Read           => ReadClassMatrix4DFromUnit
     procedure, public :: Multiply       =>  Double_x_ClassMatrix4D
     procedure, public :: Contract_124_123 => ClassMatrix4D_x_ClassMatrix4D_124_123
     generic  , public :: Add            => ClassMatrix4DAddValue, ClassMatrix4DAddMatrix
     generic  , public :: Set            => ClassMatrix4DSetValue, &
          ClassMatrix4DSet4DMat, ClassMatrix4DSetAllMatrix
     procedure, private :: ClassMatrix4DAddValue
     procedure, private :: ClassMatrix4DAddMatrix
     procedure, private :: ClassMatrix4DSetAllMatrix
     procedure, private :: ClassMatrix4DSetValue
     procedure, private :: ClassMatrix4DSet4DMat

     final       :: ClassMatrix4DFinalize

  end type ClassMatrix4D

  interface AllocateMatrix
     module procedure AllocateDoubleMatrix!, AllocateComplexMatrix
     module procedure AllocateDoubleMatrix4D
  end interface AllocateMatrix

contains

  !> Returns True if the matrix pattern is Full.
  logical function Classmatrix_IsFull(Matrix) result(IsFull)
    Class(ClassMatrix), intent(in) :: Matrix
    IsFull = ( Matrix%Pattern == MATRIX_PATTERN_FULL )
  end function Classmatrix_IsFull
  !
  !> Returns True if the matrix pattern is Dbanded.
  logical function Classmatrix_IsDbanded(Matrix) result(IsDbanded)
    Class(ClassMatrix), intent(in) :: Matrix
    IsDbanded = ( Matrix%Pattern == MATRIX_PATTERN_DBANDED )
  end function Classmatrix_IsDbanded
  !
  !> Returns True if the matrix pattern is Diagonal.
  logical function Classmatrix_IsDiagonal(Matrix) result(IsDiagonal)
    Class(ClassMatrix), intent(in) :: Matrix
    IsDiagonal = ( Matrix%Pattern == MATRIX_PATTERN_DIAGONAL )
  end function Classmatrix_IsDiagonal

  !> Retrieves the number of  subdiagonals in the matrix. 
  integer function Classmatrix_LowerBandwidth(Matrix) result(LowerBandwidth)
    Class(ClassMatrix), intent(in) :: Matrix
    LowerBandwidth = Matrix%NL
  end function Classmatrix_LowerBandwidth

  !> Retrieves the number of superdiagonals in the matrix. 
  integer function Classmatrix_UpperBandwidth(Matrix) result(UpperBandwidth)
    Class(ClassMatrix), intent(in) :: Matrix
    UpperBandwidth = Matrix%NU
  end function Classmatrix_UpperBandwidth

  !> Retrieves the number of superdiagonals in the matrix. 
  real(kind(1d0)) function ClassmatrixFullContraction(A,B) result(res)
    type(ClassMatrix), intent(in) :: A,B
    integer :: imi, ima, jmi, jma
    !************* BUGGED - PROBABLY WORKS ONLY FOR FULL MATRICES *********
    res=0.d0
    if(.not. ( allocated(A%A) .and. allocated(B%A) ) ) return
    if(min(A%NNZR,A%NNZC,B%NNZR,B%NNZC)==0)return
    if((A%NR.ne.B%NR).or.(A%NC.ne.B%NC))then
       write(*,*) A%NR,B%NR,A%NC,B%NC
       write(*,*) "ClassmatrixFullContraction: Matrices are not of the same size"
       stop
    endif
    ! jmi = max(LBOUND(A%A,2),LBOUND(B%A,2))
    ! jma = min(UBOUND(A%A,2),UBOUND(B%A,2))
    ! imi = max(LBOUND(A%A,1),LBOUND(B%A,1))
    ! ima = min(UBOUND(A%A,1),UBOUND(B%A,1))
    imi = 1
    ima = A%NRmax
    jmi = 1
    jma = A%NCmax
    res = sum( A%A(imi:ima,jmi:jma) * B%A(imi:ima,jmi:jma) )
  end function ClassmatrixFullContraction

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

    islash = index( trim(FileName), "/", BACK = .TRUE. )
    call Execute_Command_Line("mkdir -p "//FileName(:islash-1))

    open(NewUnit =  uid     , &
         File    =  FileName, &
         form    =  FileForm, &
         status  = "unknown", &
         action  = "write"  , &
         iostat  = iostat   , &
         iomsg   = iomsg    )
    if(iostat/=0)call Assert(HERE//"Open of '"//trim(FileName)//"' failed : "//trim(iomsg))
    call Matrix%Write(uid)
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
    character(len=*), parameter :: FORMAT_INTS="(*(x,i0))"
    character(len=*), parameter :: FORMAT_MAT="(*"//DOUBLE_PRINT_FORMAT//")"
    !
    !.. Write the Matrix attributes and dimensions
    write(OutUnit      , &
         FORMAT_INTS   , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix%Pattern, &
         Matrix%NR     , &
         Matrix%NC     , &
         Matrix%NU     , &
         Matrix%NL     , &
         Matrix%NNZR   , &
         Matrix%NNZC   , &
         Matrix%NRmin  , &
         Matrix%NRmax  , &
         Matrix%NCmin  , &
         Matrix%NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    if(Matrix%NR==0.or.Matrix%NC==0)return
    if(sum(abs(Matrix%A))<THRESHOLD_MATRIX_WRITE)then
       write(OutUnit,"(i0)") 0
       return
    endif
    write(OutUnit,"(i0)") 1
    !.. Write Matrix
    do j=Matrix%NCmin,min(Matrix%NCmax,Matrix%NCmin+Matrix%NNZC-1)
       write(OutUnit      , &
            FORMAT_MAT    , &
            IOSTAT=iostat , &
            IOMSG =iomsg  ) &
            (Matrix%A(i,j),&
            i=Matrix%NRmin,min(Matrix%NRmax,Matrix%NRmin+Matrix%NNZR-1))
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
         Matrix%Pattern, &
         Matrix%NR     , &
         Matrix%NC     , &
         Matrix%NU     , &
         Matrix%NL     , &
         Matrix%NNZR   , &
         Matrix%NNZC   , &
         Matrix%NRmin  , &
         Matrix%NRmax  , &
         Matrix%NCmin  , &
         Matrix%NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    if(Matrix%NR==0.or.Matrix%NC==0)return
    if(sum(abs(Matrix%A))<THRESHOLD_MATRIX_WRITE)then
       write(OutUnit) 0
       return
    endif
    write(OutUnit) 1
    !.. Write Matrix
    write(OutUnit      , &
         IOSTAT=iostat , &
         IOMSG =iomsg  ) &
         ((Matrix%A(i,j),&
         i=Matrix%NRmin,min(Matrix%NRmax,Matrix%NRmin+Matrix%NNZR-1)),&
         j=Matrix%NCmin,min(Matrix%NCmax,Matrix%NCmin+Matrix%NNZC-1))
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
    call Matrix%Free()
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
    call Matrix%Read(uid)
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
    call Matrix%Free()
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
    integer :: iostat, i, j, izero
    character(len=IOMSG_LENGTH) :: iomsg=" "
    character(len=*), parameter :: FORMAT_INTS="(*(x,i0))"
    character(len=*), parameter :: FORMAT_MAT="(*"//DOUBLE_PRINT_FORMAT//")"
    !
    !.. Read the Matrix attributes and dimensions
    read(InUnit,*, &
         IOSTAT=iostat   , &
         IOMSG=iomsg     ) &
         Matrix%Pattern  , &
         Matrix%NR       , &
         Matrix%NC       , &
         Matrix%NU       , &
         Matrix%NL       , &
         Matrix%NNZR     , &
         Matrix%NNZC     , &
         Matrix%NRmin    , &
         Matrix%NRmax    , &
         Matrix%NCmin    , &
         Matrix%NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    if(Matrix%NR==0.or.Matrix%NC==0)return
    call AllocateMatrix( Matrix%A  , &
         Matrix%NRmin, Matrix%NRmax, &
         Matrix%NCmin, Matrix%NCmax  )
    read(InUnit,*) izero
    if(izero==0)then
       matrix%A=0.d0
       return
    endif
    !
    !.. Read Matrix
    do j=Matrix%NCmin,min(Matrix%NCmax,Matrix%NCmin+Matrix%NNZC-1)
       read(InUnit,FMT=FORMAT_MAT  , &
            IOSTAT=iostat , &
            IOMSG =iomsg  ) &
            (Matrix%A(i,j), &
            i=Matrix%NRmin,min(Matrix%NRmax,Matrix%NRmin+Matrix%NNZR-1))
       if(iostat/=0)call Assert(iomsg)
    enddo
    !
  end subroutine ReadClassMatrixFromFormattedUnit
  !
  subroutine ReadClassMatrixFromUnformattedUnit(Matrix,InUnit)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: InUnit
    integer :: iostat, i, j, izero
    character(len=IOMSG_LENGTH) :: iomsg
    ! 
    !.. Write the Matrix attributes and dimensions
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         Matrix%Pattern, &
         Matrix%NR     , &
         Matrix%NC     , &
         Matrix%NU     , &
         Matrix%NL     , &
         Matrix%NNZR   , &
         Matrix%NNZC   , &
         Matrix%NRmin  , &
         Matrix%NRmax  , &
         Matrix%NCmin  , &
         Matrix%NCmax 
    if(iostat/=0)call Assert(iomsg)
    !
    if(Matrix%NR==0.or.Matrix%NC==0)return
    call AllocateMatrix( Matrix%A  , &
         Matrix%NRmin, Matrix%NRmax, &
         Matrix%NCmin, Matrix%NCmax  )
    read(InUnit) izero
    if(izero==0)then
       matrix%A=0.d0
       return
    endif
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG =iomsg  ) &
         ((Matrix%A(i,j),&
         i=Matrix%NRmin,min(Matrix%NRmax,Matrix%NRmin+Matrix%NNZR-1)),&
         j=Matrix%NCmin,min(Matrix%NCmax,Matrix%NCmin+Matrix%NNZC-1))
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine ReadClassMatrixFromUnformattedUnit


  !> Set all the elements of the matrix equal to real number.
  subroutine AddDoubleToClassMatrix(Matrix,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: x
    if(.not.allocated(Matrix%A))then
       call Assert("Matrix not initialized")
    else
       Matrix%A=x
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
    if(.not.allocated(Matrix%A))then
       call Assert("Matrix not initialized")
    else
       Matrix%A=Matrix%A*x
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
    if(.not.allocated(self%A))then
       call Assert(HERE//"Matrix not initialized")
    else
       B=MatFactor
       call self%Multiply(B,Side,MatFactorTRANS)
    endif
  end subroutine ClassMatrix_x_Matrix

  !> Multiplies ClassMatrix by a diagonal matrix, expressed as an array of doubles,
  !! from either the left or the right. It currently works only with full matrices
  subroutine ClassMatrix_x_vec(self, vec, Side )
    Class(ClassMatrix), intent(inout) :: self
    DoublePrecision   , intent(in)    :: Vec(:)
    character(len=*),   intent(in)    :: Side
    integer :: i
    if(side.is."left")then
       do i=1,self%NR
          self%A(i,:)=self%A(i,:)*vec(i)
       enddo
    else
       do i=1,self%NC
          self%A(:,i)=self%A(:,i)*vec(i)
       enddo
    endif
  end subroutine ClassMatrix_x_vec


  !> Multiplies the matrices belonging to two different matrix classes
  !*** Does not treat consistently the multiplication with dbanded matrices
  subroutine ClassMatrix_x_ClassMatrix( MatrixA, MatrixB, Side, MatrixBType )
    !
    !********** MUST DEBUG !!!! *******************
    !
    Class(ClassMatrix), target, intent(inout) :: MatrixA
    Class(ClassMatrix), target, intent(in)    :: MatrixB
    character(len=*),           intent(in)    :: Side
    character(len=*),           intent(in)    :: MatrixBType
    !
    Class(ClassMatrix), pointer :: PtrA, PtrB
    character         :: TypeA, TypeB
    Type(ClassMatrix) :: MatrixC
    integer :: m, n, k, lda, ldb, ldc, nr, nc
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

    lda = size( PtrA%A, 1 )
    if( TypeA .is. "N" )then
       !m = min( size( PtrA%A, 1 ), PtrA%NNZR )
       !write(*,*) m,size( PtrA%A, 1 ), PtrA%NNZR
       m = size( PtrA%A, 1 )
       nr= PtrA%NR
       if(TypeB .is. "N")then
          !k = min( size( PtrA%A, 2 ), size( PtrB%A, 1 ), PtrA%NNZC, PtrB%NNZR )
          k = min( size( PtrA%A, 2 ), size( PtrB%A, 1 ) )
       else
          !k = min( size( PtrA%A, 2 ), size( PtrB%A, 2 ), PtrA%NNZC, PtrB%NNZC )
          k = min( size( PtrA%A, 2 ), size( PtrB%A, 2 ))
       endif
    else
       !m = min( size( PtrA%A, 2 ), PtrA%NNZC )
       m = size( PtrA%A, 2 )
       nr= PtrA%NC
       if(TypeB .is. "N")then
          !k = min( size( PtrA%A, 1 ), size( PtrB%A, 1 ), PtrA%NNZR, PtrB%NNZR )
          k = min( size( PtrA%A, 1 ), size( PtrB%A, 1 ))
       else
          !k = min( size( PtrA%A, 1 ), size( PtrB%A, 2 ), PtrA%NNZR, PtrB%NNZC )
          k = min( size( PtrA%A, 1 ), size( PtrB%A, 2 ))
       endif
    endif

    ldb = size( PtrB%A, 1 )
    if(TypeB .is. "N")then
       !n = min( size( PtrB%A, 2 ), PtrB%NNZC )
       n =  size( PtrB%A, 2 )
       nc= PtrB%NC
    else
       !n = min( size( PtrB%A, 1 ), PtrB%NNZR )
       n = size( PtrB%A, 1 )
       nc= PtrB%NR
    endif


    ldc = m
    call MatrixC%InitFull( nr, nc, m, n )
    !
    ! write(*,*) "lda",TypeA, TypeB, m, n, k, 1.d0, &
    !      PtrA%A,    lda, &
    !      PtrB%A,    ldb, 0.d0, &
    !      MatrixC%A, ldc
    call DGEMM( &
         TypeA, TypeB, m, n, k, 1.d0, &
         PtrA%A,    lda, &
         PtrB%A,    ldb, 0.d0, &
         MatrixC%A, ldc )
    !write(*,*) "ldc",ldc
    !pause
    !
    MatrixA = MatrixC
    !
    call MatrixC%Free()
    !
  end subroutine ClassMatrix_x_ClassMatrix


    !> Multiplies the matrices belonging to two different matrix classes
  !*** Does not treat consistently the multiplication with dbanded matrices
  function ClassMatrix_x_ClassMatrix4D12( MatrixA, MatrixB, vimi, vima ) result( MatrixC )
    !
    Class(ClassMatrix)  , intent(in) :: MatrixA
    type (ClassMatrix4D), intent(in) :: MatrixB
    type (ClassMatrix)  , pointer    :: MatrixC
    integer             , intent(in) :: vimi(4), vima(4)
    !
    integer         :: i1,i2,i3,i4,j1,j2,j3,j4
    real(kind(1d0)) :: dBuf
    !
    allocate(MatrixC)
    call MatrixC%InitFull( vima(3)-vimi(3)+1, vima(4)-vimi(4)+1 )
    MatrixC = 0.d0
    !write(*,*)
    do i4 = vimi(4), vima(4)
       j4 = i4 - vimi(4) + 1
       do i3 = vimi(3), vima(3)
          j3 = i3 - vimi(3) + 1
          dBuf = 0.d0
          do i2 = vimi(2), vima(2)
             j2 = i2 - vimi(2) + 1
             do i1 = vimi(1), vima(1)
                j1 = i1 - vimi(1) + 1
                dBuf = dBuf + MatrixA%A(j1,j2) * MatrixB%A(j1,j2,j3,j4)
             enddo
          enddo
          MatrixC%A(j3,j4)=dBuf
         ! MatrixC%A(i3-vimi(3)+1,i4-vimi(4)+1)=dBuf
       enddo
    enddo
    !
  end function ClassMatrix_x_ClassMatrix4D12

    !> Multiplies the matrices belonging to two different matrix classes
  !*** Does not treat consistently the multiplication with dbanded matrices
  function ClassMatrix_x_ClassMatrix4D13( MatrixA, MatrixB, vimi, vima ) result( MatrixC )
    !
    Class(ClassMatrix)  , intent(in) :: MatrixA
    type (ClassMatrix4D), intent(in) :: MatrixB
    type (ClassMatrix)  , pointer    :: MatrixC
    integer             , intent(in) :: vimi(4), vima(4)
    !
    integer         :: i1,i2,i3,i4,j1,j2,j3,j4
    real(kind(1d0)) :: dBuf
    !
    allocate(MatrixC)
    call MatrixC%InitFull( vima(2)-vimi(2)+1, vima(4)-vimi(4)+1 )
    MatrixC = 0.d0
    do i4 = vimi(4), vima(4)
       j4 = i4 - vimi(4) + 1
       do i2 = vimi(2), vima(2)
          j2 = i2 - vimi(2) + 1
          dBuf = 0.d0
          do i3 = vimi(3), vima(3)
             j3 = i3 - vimi(3) + 1
             do i1 = vimi(1), vima(1)
                j1 = i1 - vimi(1) + 1
                dBuf = dBuf + MatrixA%A(j1,j3) * MatrixB%A(j1,j2,j3,j4)
             enddo
          enddo
          !MatrixC%A(i2-vimi(2)+1,i4-vimi(4)+1)=dBuf
          MatrixC%A(j2,j4)=dBuf
       enddo
    enddo
    !
  end function ClassMatrix_x_ClassMatrix4D13


  !> Assign a real number to every element of the ClassMatrix's matrix.
  subroutine AssignDoubleToClassMatrix(Matrix,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: x
    if(.not.allocated(Matrix%A))then
       call Assert("Matrix not initialized")
    else
       Matrix%A=x
    endif
  end subroutine AssignDoubleToClassMatrix


  !> Copies the information of a matrix class to another.
  subroutine CopyClassMatrixToClassMatrix(MatrixOut,MatrixInp)
    Class(ClassMatrix), intent(inout) :: MatrixOut
    type (ClassMatrix), intent(in)    :: MatrixInp
    integer :: LBR, UBR, LBC, UBC, iRow,iCol
    !
    call MatrixOut%Free()
    !
    MatrixOut%Pattern = MatrixInp%Pattern
    !
    MatrixOut%NR = MatrixInp%NR
    MatrixOut%NC = MatrixInp%NC
    MatrixOut%NL = MatrixInp%NL
    MatrixOut%NU = MatrixInp%NU
    !
    MatrixOut%NRMin = MatrixInp%NRMin
    MatrixOut%NRMax = MatrixInp%NRMax
    MatrixOut%NCMin = MatrixInp%NCMin
    MatrixOut%NCMax = MatrixInp%NCMax
    !
    LBR=LBOUND(MatrixInp%A,1)
    UBR=UBOUND(MatrixInp%A,1)
    LBC=LBOUND(MatrixInp%A,2)
    UBC=UBOUND(MatrixInp%A,2)
!!$    allocate(MatrixOut%A,source=MatrixInp%A)
    allocate(MatrixOut%A(LBR:UBR,LBC:UBC))
    MatrixOut%A=0.d0
    do iCol=LBC, UBC
       do iRow=LBR,UBR
          MatrixOut%A(iRow,iCol)=MatrixInp%A(iRow,iCol)
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
    call MatrixOut%Free()
    NR=UBOUND(MatrixInp,1)-LBOUND(MatrixInp,1)+1
    NC=UBOUND(MatrixInp,2)-LBOUND(MatrixInp,2)+1
    call MatrixOut%initFull(NR,NC)
    MatrixOut%A=MAtrixInp
  end subroutine CopyMatrixToClassMatrix


  !> Assign to a ClassMatrix with Full pattern,
  !! the content of a bidimensional array
  subroutine CopyVectorToClassMatrix(MatrixOut,VectorInp)
    Class(ClassMatrix), intent(inout) :: MatrixOut
    DoublePrecision   , intent(in)    :: VectorInp(:)
    integer :: NR, NC
    call MatrixOut%Free()
    NR=UBOUND(VectorInp,1)-LBOUND(VectorInp,1)+1
    NC=1
    call MatrixOut%initFull(NR,NC)
    MatrixOut%A(:,1)=VectorInp
  end subroutine CopyVectorToClassMatrix


  !> Create and allocate a ClassMatrix object with Full 
  !! pattern, if it is not defined, and free it before, if it is.
  subroutine Classmatrix_InitFull( Matrix, NR, NC, NNZR, NNZC )
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR, NC
    integer, optional , intent(in)    :: NNZR, NNZC
    !
    call Matrix%Free()
    call SetMatrixNominalSize( Matrix, NR, NC )
    Matrix%NL=NR-1
    Matrix%NU=NC-1
    if(present(NNZR))then
       Matrix%NNZR=NNZR
    else
       Matrix%NNZR=NR
    endif
    if(present(NNZC))then
       Matrix%NNZC=NNZC
    else
       Matrix%NNZC=NC
    endif
    call SetMatrixPhysicalSize( Matrix, 1, Matrix%NNZR, 1, Matrix%NNZC )
    call AllocateMatrix( Matrix%A,   &
         Matrix%NRmin, Matrix%NRmax, &
         Matrix%NCmin, Matrix%NCmax  )
    Matrix%Pattern = MATRIX_PATTERN_FULL
    if(allocated(Matrix%A)) Matrix%A = 0.d0
    !
  end subroutine Classmatrix_InitFull


  !> Create and allocate a ClassMatrix object with Dbanded  pattern, if it is not defined, and free it before, if it is not.
  subroutine Classmatrix_InitDbanded( Matrix, NR, NC, NL, NU )
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR
    integer           , intent(in)    :: NC
    integer           , intent(in)    :: NL
    integer           , intent(in)    :: NU
    !
    call Matrix%Free()
    call SetMatrixNominalSize( Matrix, NR, NC )
    Matrix%NL=max(0,min(NL,NR-1))
    Matrix%NU=max(0,min(NU,NC-1))
    !*** CHECK
    Matrix%NNZR=NR
    Matrix%NNZC=NC
    call SetMatrixPhysicalSize( Matrix, -Matrix%NU, Matrix%NL, 1, Matrix%NC )
    call AllocateMatrix( Matrix%A,   &
         Matrix%NRmin, Matrix%NRmax, &
         Matrix%NCmin, Matrix%NCmax  )
    Matrix%Pattern = MATRIX_PATTERN_DBANDED
    Matrix%A = 0.d0
    !
  end subroutine Classmatrix_InitDbanded


  !> Create and allocate a ClassMatrix object with Dbanded  pattern, if it is not defined, and free it before, if it is not.
  subroutine Classmatrix_InitSubmatrix( Matrix, RMask, CMask )
    Class(ClassMatrix), intent(inout) :: Matrix
    logical           , intent(in)    :: RMask(:), CMask(:)
    !
    call Matrix%Free()
    Matrix%NR    = size(RMask)
    Matrix%NC    = size(CMask)
    !*** CHECK
    Matrix%NNZR=Matrix%NR
    Matrix%NNZC=Matrix%NC
    Matrix%NL    = Matrix%NR-1
    Matrix%NU    = Matrix%NC-1
    Matrix%NRMin = 1
    Matrix%NCMin = 1
    Matrix%NRMax = count(RMask)
    Matrix%NCMax = count(CMask)
    allocate(Matrix%RMask,source=RMask)
    allocate(Matrix%CMask,source=CMask)
    call AllocateMatrix( Matrix%A,   &
         Matrix%NRmin, Matrix%NRmax, &
         Matrix%NCmin, Matrix%NCmax  )
    Matrix%Pattern = MATRIX_PATTERN_SUBMATRIX
    Matrix%size=size(Matrix%A,1)*size(Matrix%A,1)
    if(Matrix%size>0)Matrix%A = 0.d0
    !
  end subroutine Classmatrix_InitSubmatrix


  !> Set the rows and comlumns minimum and maximum indexes.
  subroutine SetMatrixPhysicalSize( Matrix, NRMin, NRMax, NCMin, NCMax )
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NRMin, NRMax, NCMin, NCMax
    Matrix%NRmin = NRMin
    Matrix%NRmax = NRMax
    Matrix%NCmin = NCMin
    Matrix%NCmax = NCMax
  end subroutine SetMatrixPhysicalSize


  !> Fetches the requested matrix element.
  DoublePrecision function Classmatrix_Element(Matrix,i,j) result(Element) 
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    Element=0.d0
    select case( Matrix%Pattern )
    case( MATRIX_PATTERN_FULL )
       Element = ClassMatrix_FullElement(Matrix,i,j)
    case( MATRIX_PATTERN_DBANDED )
       Element = ClassMatrix_DbandedElement(Matrix,i,j)
    case DEFAULT
    end select
  end function Classmatrix_Element

  !> Sets the requested matrix element.
  subroutine Classmatrix_SetElement(Matrix,i,j,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    DoublePrecision   , intent(in)    :: x
    select case( Matrix%Pattern )
    case( MATRIX_PATTERN_FULL )
       call Classmatrix_FullSetElement(Matrix,i,j,x)
    case( MATRIX_PATTERN_DBANDED )
       call Classmatrix_DbandedSetElement(Matrix,i,j,x)
    case DEFAULT
    end select
  end subroutine Classmatrix_SetElement

  !> Sets the diagonal of the matrix to a given scalar
  subroutine Classmatrix_SetDiagonal(self,x)
    Class(ClassMatrix), intent(inout) :: self
    DoublePrecision   , intent(in)    :: x
    integer :: i
    do i=1,min(self%NRows(),self%NColumns())
       call self%Set(i,i,x)
    enddo
  end subroutine Classmatrix_SetDiagonal

  !> Add to the diagonal of the matrix to a given scalar
  subroutine Classmatrix_AddDiagonal(self,x)
    Class(ClassMatrix), intent(inout) :: self
    DoublePrecision   , intent(in)    :: x
    real(kind(1d0)) :: val
    integer :: i
    do i=1,min(self%NRows(),self%NColumns())
       val = self%Element(i,i) + x
       call self%Set(i,i,val)
    enddo
  end subroutine Classmatrix_AddDiagonal

  !> Sets the requested matrix element when the array is in Full representation.
  subroutine Classmatrix_FullSetElement(Matrix,i,j,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    DoublePrecision   , intent(in)    :: x
    call CheckIndexBounds(Matrix,i,j)
    Matrix%A(i,j)=x
  end subroutine Classmatrix_FullSetElement

  !> Sets the requested matrix element when the array is in Dbanded representation.
  subroutine Classmatrix_DbandedSetElement(Matrix,i,j,x)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: i,j
    DoublePrecision   , intent(in)    :: x
    integer :: iPhys
    call CheckIndexBounds(Matrix,i,j)
    iPhys=i-j
    if( iPhys < -Matrix%NU .or. iPhys > Matrix%NL )&
         call Assert("Invalid dbanded matrix index")
    Matrix%A(iPhys,j)=x
  end subroutine Classmatrix_DbandedSetElement


  !> Checks that the matrix element indexed have valid values.
  subroutine CheckIndexBounds(Matrix,i,j)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    if( i < 1 .or. i > Matrix%NR )call Assert("Row index is off bound")
    if( j < 1 .or. j > Matrix%NC )call Assert("Column index is off bound")
  end subroutine CheckIndexBounds


  !> Fetches the requested matrix elements when the array is stored in Full representation.
  DoublePrecision function ClassMatrix_FullElement(Matrix,i,j) result(Element)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    integer  :: iPhys,jPhys
    call CheckIndexBounds(Matrix,i,j)
    iPhys=Matrix%NRMin-1+i
    jPhys=Matrix%NCMin-1+j
    Element=Matrix%A(iPhys,jPhys)
  end function ClassMatrix_FullElement


  !> Fetches the requested matrix elements when the array is stored in Dbanded representation.
  DoublePrecision function ClassMatrix_DbandedElement(Matrix,i,j) result(Element)
    Class(ClassMatrix), intent(in) :: Matrix
    integer           , intent(in) :: i,j
    integer :: iPhys
    call CheckIndexBounds(Matrix,i,j)
    iPhys=i-j
    if( iPhys >= -Matrix%NU .and. iPhys <= Matrix%NL )then
       Element=Matrix%A(iPhys,j)
    else
       Element=0.d0
    end if
  end function ClassMatrix_DbandedElement


  !> Sets the number of rows and columns of the matrix.
  subroutine SetMatrixNominalSize(Matrix,NR,NC)
    Class(ClassMatrix), intent(inout) :: Matrix
    integer           , intent(in)    :: NR,NC
    if(NR<0)call Assert("Invalid Number of Rows")
    if(NC<0)call Assert("Invalid Number of Columns")
    Matrix%NR=NR
    Matrix%NC=NC
  end subroutine SetMatrixNominalSize


  !> Allocates a real matrix and initializes it to zero.
  subroutine AllocateDoubleMatrix(A,NRMin,NRMax,NCMin,NCMax)
    DoublePrecision, allocatable, intent(inout) :: A(:,:)
    integer                     , intent(in)    :: NRMin, NRmax
    integer                     , intent(in)    :: NCMin, NCmax
    integer :: Stat=0
    character(len=IOMSG_LENGTH) :: iomsg=" "
    if( NRMax  < NRMin-1 .or. NCMax  < NCMin-1 )return
    if( NRMax == NRMin-1 .or. NCMax == NCMin-1 )then
       allocate(A(NRMax-NRmin+1,NCMax-NCmin+1))
       return
    endif
    if(allocated(A))deallocate(A)
    allocate(A(NRMin:NRMax,NCMin:NCMax),STAT=Stat,ERRMSG=iomsg)
    if(Stat/=0)call Assert(iomsg)
    if(size(A,1)*size(A,2)>0)A=0.d0
  end subroutine AllocateDoubleMatrix

  !> Allocates a 4D real matrix and initializes it to zero.
  subroutine AllocateDoubleMatrix4D(A,NMin1,NMax1,NMin2,NMax2,NMin3,NMax3,NMin4,NMax4)
    DoublePrecision, allocatable, intent(inout) :: A(:,:,:,:)
    integer                     , intent(in)    :: NMin1,NMax1
    integer                     , intent(in)    :: NMin2,NMax2
    integer                     , intent(in)    :: NMin3,NMax3
    integer                     , intent(in)    :: NMin4,NMax4
    integer :: Stat=0
    character(len=IOMSG_LENGTH) :: iomsg=" "
    if( NMax1 < NMin1 .or. NMax2 < NMin2 )return
    if( NMax3 < NMin3 .or. NMax4 < NMin4 )return
    if(allocated(A))deallocate(A)
    allocate(A(NMin1:NMax1,NMin2:NMax2,NMin3:NMax3,NMin4:NMax4),STAT=Stat,ERRMSG=iomsg)
    if(Stat/=0)call Assert(iomsg)
    A=0.d0
  end subroutine AllocateDoubleMatrix4D

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
  subroutine Classmatrix_Free( Matrix )
    Class(ClassMatrix), intent(inout) :: Matrix
    if(allocated(Matrix%A)) deallocate(Matrix%A)
    Matrix%NR=0
    Matrix%NC=0
    Matrix%Pattern=MATRIX_PATTERN_UNDEFINED
    Matrix%NU=0
    Matrix%NL=0
  end subroutine Classmatrix_Free


  !> Free allocatable space and sets to zero all the members of a ClassMatrix object.
  subroutine Classmatrix_Finalize( Matrix )
    type(ClassMatrix), intent(inout) :: Matrix
    call Matrix%Free()
  end subroutine Classmatrix_Finalize


  !> Retrieves the number of rows of the matrix in ClassMatrix%
  integer function Classmatrix_NRows( Matrix ) result( NRows )
    Class(ClassMatrix), intent(in) :: Matrix
    NRows = Matrix%NR
  end function Classmatrix_NRows
  !
  !> Retrieves the number of columns of the matrix in ClassMatrix%
  integer function Classmatrix_NColumns( Matrix ) result( NColumns )
    Class(ClassMatrix), intent(in) :: Matrix
    NColumns = Matrix%NC
  end function Classmatrix_NColumns


  !> Gets a rectangular submatrix from the ClassMatrix's matrix 
  !! defined by two independent vectors of row and column indexes
  subroutine ClassMatrix_GetSubmatrix( &
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
    character(len=*), parameter :: HERE = "ClassMatrix::ClassMatrix_GetSubmatrix : "
    !
    call CheckIndexes(NIndexesBra,IndexesBra,Matrix%NR)
    call CheckIndexes(NIndexesKet,IndexesKet,Matrix%NC)
    !
    if( Matrix%IsFull() )then
       call SubMatrix%InitFull( NIndexesBra, NIndexesKet )
    elseif( Matrix%IsDbanded() )then
       call AdditionalCheckIndexesDbanded(NIndexesBra,IndexesBra)
       call AdditionalCheckIndexesDbanded(NIndexesKet,IndexesKet)
       VerticalShift = IndexesBra(1) - IndexesKet(1)
       NL = Matrix%LowerBandwidth() - VerticalShift
       NU = Matrix%UpperBandwidth() + VerticalShift
       call SubMatrix%InitDbanded( NIndexesBra, NIndexesKet, NL, NU )
    else
       call Assert("Unrecognized Matrix type")
    endif
    !
    do iSubCol = 1, NIndexesKet
       iCol = IndexesKet( iSubCol )
       do iSubRow = 1, NIndexesBra
          iRow = IndexesBra( iSubRow )
          Element = Matrix%Element( iRow, iCol )
          call SubMatrix%SetElement( iSubRow, iSubCol, Element )
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
    !.. For dbanded matrices, it is required that only extremal values are
    !   dropped, i.e., that the set of indexes is contiguous. In this way,
    !   the resulting matrix is still dbanded.
    subroutine AdditionalCheckIndexesDbanded(NIndexes,Indexes)
      integer, intent(in) :: NIndexes
      integer, intent(in) :: Indexes(:)
      integer :: i
      do i=1,NIndexes-1
         if(Indexes(i+1)/=Indexes(i)+1)&
              call Assert(HERE//"Only contiguous set of subindexes are allowed for dbanded matrices")
      enddo
    end subroutine AdditionalCheckIndexesDbanded
  end subroutine ClassMatrix_GetSubmatrix


  !> Muliplies the matrix contained in a matrix class by a real number.
  subroutine Classmatrix_TimesDouble( Matrix, Number )
    Class(ClassMatrix), intent(inout) :: Matrix 
    DoublePrecision,    intent(in)    :: Number
    Matrix%A = Number * Matrix%A
  end subroutine Classmatrix_TimesDouble


  !> Adds up two matrices belonging to two different matrix classes.
  subroutine ClassMatrix_AddClassMatrix( Matrix, DeltaMatrix ) 
    Class(ClassMatrix), intent(inout) :: Matrix
    Class(ClassMatrix), intent(in)    :: DeltaMatrix
    !
    integer :: LBR, UBR, LBC, UBC, iRow,iCol
    !
    LBR=LBOUND(Matrix%A,1)
    UBR=UBOUND(Matrix%A,1)
    LBC=LBOUND(Matrix%A,2)
    UBC=UBOUND(Matrix%A,2)
    !
    do iCol=LBC, UBC
       do iRow=LBR,UBR
          Matrix%A(iRow,iCol)=Matrix%A(iRow,iCol)+DeltaMatrix%A(iRow,iCol)
       enddo
    enddo
  end subroutine ClassMatrix_AddClassMatrix


  !> Performs the summation of the matrix elements absolute value.
  DoublePrecision function Classmatrix_Norm1( Matrix ) result( Norm1 )
    Class(ClassMatrix), intent(in) :: Matrix
    Norm1 = sum(abs(Matrix%A))
  end function Classmatrix_Norm1

  !> Performs the squared root , of the matrix elements absolute squared value summation.
  DoublePrecision function Classmatrix_Norm2( Matrix ) result( Norm2 ) 
    Class(ClassMatrix), intent(in) :: Matrix
    Norm2 = sqrt(sum(abs(Matrix%A)**2))
  end function Classmatrix_Norm2

  subroutine Classmatrix_Transpose( self )
    class(ClassMatrix), intent(inout)    :: self
    integer :: i,j,nr, nc, nu, nl, nrmi, nrma, ncmi, ncma
    real(kind(1d0)), allocatable :: A(:,:)
    nr         = self%nc
    nc         = self%nr
    nu         = self%nl
    nl         = self%nu
    nrmi       = self%nrmax
    nrma       = self%nrmin
    ncmi       = self%ncmax
    ncma       = self%ncmin
    self%nc    = nc
    self%nr    = nr
    self%nl    = nl
    self%nu    = nu
    self%nrmax = nrma
    self%nrmin = nrmi
    self%ncmax = ncma
    self%ncmin = ncmi
    if(allocated(self%A))then
       allocate(A(nr,nc))
       do i=1,nr
          do j=1,nc
             A(i,j)=self%A(j,i)
          enddo
       enddo
       deallocate(self%A)
       allocate(self%A,source=A)
       deallocate(A)
    endif
  end subroutine Classmatrix_Transpose


  !> initialize a ClassMatrix with Full pattern and get as matrix an external one.
  subroutine Classmatrix_SetMatrix_( Matrix, MatrixArray )
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: MatrixArray(:,:)
    call Matrix%InitFull( size(MatrixArray,1), size(MatrixArray,2) )
    Matrix%A = MatrixArray
  end subroutine Classmatrix_SetMatrix_

  !> Sets a subblock of a matrix
  !! The receiving matrix is assumed to be initialized already. 
  subroutine Classmatrix_SetSubMatrix( self, Mat, rmin, cmin )
    Class(ClassMatrix), intent(inout) :: self
    DoublePrecision   , intent(in)    :: Mat(:,:)
    integer           , intent(in)    :: rmin, cmin
    integer :: nr,nc
    nr=size(Mat,1)
    nc=size(Mat,2)
    self%A(:,:) = Mat(1:rmin,1:cmin)
    self%A(rmin:rmin+nr-1,cmin:cmin+nc-1) = Mat
  end subroutine Classmatrix_SetSubMatrix

  !> Sets a subblock of a matrix. Works only with full matrices, so far
  !! The receiving matrix is assumed to be initialized already. 
  subroutine Classmatrix_SetSubMatrixCM( self, Mat, rmin, cmin )
    Class(ClassMatrix), intent(inout) :: self
    type (ClassMatrix), intent(in)    :: Mat
    integer           , intent(in)    :: rmin, cmin
    integer :: nr,nc
    nr=size(Mat%A,1)
    nc=size(Mat%A,2)
    self%A(rmin:rmin+nr-1,cmin:cmin+nc-1) = Mat%A
  end subroutine Classmatrix_SetSubMatrixCM

  !> initialize a ClassMatrix with Full pattern and get as matrix an external one.
  subroutine Classmatrix_SetVector( Matrix, VectorArray )
    Class(ClassMatrix), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: VectorArray(:)
    call Matrix%InitFull( size(VectorArray,1), 1 )
    Matrix%A(:,1) = VectorArray
  end subroutine Classmatrix_SetVector

  subroutine Classmatrix_FetchMatrix( Matrix, Array )
    Class(ClassMatrix),           intent(in) :: Matrix
    real(kind(1d0)), allocatable, intent(out):: Array(:,:)
    allocate( Array, source = Matrix%A )
  end subroutine Classmatrix_FetchMatrix

  function Classmatrix_FetchColumn( self, iCol ) result(vec)
    Class(ClassMatrix),           intent(in) :: self
    integer           ,           intent(in) :: iCol
    real(kind(1d0)), allocatable             :: vec(:)
    if(allocated(vec))then
       if(size(vec)/=size(self%A,1))deallocate(vec)
    endif
    if(.not.allocated(vec)) allocate(vec(size(self%A,1)))
    vec = self%A(:,iCol)
  end function Classmatrix_FetchColumn

  subroutine Classmatrix_BuildUpMatrix( Self, Blocks )
    class(ClassMatrix), intent(inout) :: Self
    class(ClassMatrix), intent(in)    :: Blocks(:,:)
    integer :: RowBeg, RowEnd, ColBeg, ColEnd, TotNumRows
    integer :: i, j, NumBra, NumKet, TotNumColumns
    call Self%Free()
    NumBra = size(Blocks,1)
    NumKet = size(Blocks,2)
    TotNumColumns = 0
    do j = 1, NumKet
       TotNumColumns = TotNumColumns + Blocks(1,j)%NColumns()
    end do
    TotNumRows = 0
    do i = 1, NumBra
       TotNumRows = TotNumRows + Blocks(i,1)%NRows()
    end do
    call Self%InitFull( TotNumRows, TotNumColumns )
    RowBeg = 1
    ColBeg = 1
    do j = 1, NumKet
       do  i = 1, NumBra 
          RowEnd = RowBeg + Blocks(i,j)%NRows() - 1
          ColEnd = ColBeg + Blocks(i,j)%NColumns() - 1
          call Self%Plug( Blocks(i,j), RowBeg, ColBeg )
          RowBeg = RowEnd + 1
          if ( i == NumBra ) then
             ColBeg = ColEnd + 1
             RowBeg = 1
          end if
       end do
    end do
  end subroutine Classmatrix_BuildUpMatrix

  
  subroutine Classmatrix_Plug( self, BlockMat, RowBeg, ColBeg  )
    class(ClassMatrix), intent(inout) :: self
    type(ClassMatrix) , intent(in)    :: BlockMat
    integer           , intent(in)    :: RowBeg, ColBeg
    integer :: i, j
    do j = ColBeg, ColBeg + BlockMat%NColumns() - 1
       do i = RowBeg, RowBeg + BlockMat%NRows() - 1
          call self%SetElement( i, j, BlockMat%Element(i-RowBeg+1,j-ColBeg+1) )
       end do
    end do
  end subroutine Classmatrix_Plug


  logical function Classmatrix_IsInitialized( Mat ) result(IsInit)
    !
    class(ClassMatrix), intent(in) :: Mat
    !
    if ( allocated(Mat%A) ) then
       IsInit = .true.
    else
       IsInit = .false.
    end if
    !
  end function Classmatrix_IsInitialized


  !*** MUST BE ADAPTED TO BANDED MATRICES!!!!
  subroutine Classmatrix_Drop( self, sdim, mi,ma ) 
    Class(ClassMatrix), intent(inout) :: self
    character(len=*)  , intent(in)    :: sdim
    integer           , intent(in)    :: mi, ma
    integer                           :: i, ndrop
    ndrop = ma-mi+1
    if(sdim.is."columns")then
       do i=mi,self%NC-ndrop
          self%A(:,i)=self%A(:,i+ndrop)
       enddo
       self%NCMax = self%NCMax - ndrop
    else
       do i=mi,self%NR-ndrop
          self%A(i,:)=self%A(i+ndrop,:)
       enddo
       self%NRMax = self%NRMax - ndrop
    endif
  end subroutine ClassMatrix_Drop


  !============= 4D Matrices =======================================

  logical function ClassMatrix4DIsInitialized( Mat ) result(IsInit)
    !
    class(ClassMatrix4D), intent(in) :: Mat
    !
    if ( allocated(Mat%A) ) then
       IsInit = .true.
    else
       IsInit = .false.
    end if
    !
  end function ClassMatrix4DIsInitialized

  subroutine ClassMatrix4DInit( self, n1, n2, n3, n4 )
    class(ClassMatrix4D), intent(inout) :: self
    integer             , intent(in)    :: n1, n2, n3, n4
    call self%Free()
    allocate(self%A(n1,n2,n3,n4))
    self%A=0.d0
  end subroutine ClassMatrix4DInit

  integer function ClassMatrix4DSize( self, iDim ) result( iSize )
    class(ClassMatrix4D), intent(in) :: self
    integer             , intent(in) :: iDim
    iSize = size( self%A, iDim )
  end function ClassMatrix4DSize

  integer function ClassMatrix4DTotalSize( self ) result( iSize )
    class(ClassMatrix4D), intent(in) :: self
    integer :: iDim
    iSize = 1
    do iDim=1,4
       iSize = iSize * size( self%A, iDim )
    enddo
  end function ClassMatrix4DTotalSize

  subroutine ClassMatrix4DFree( self )
    class(ClassMatrix4D), intent(inout) :: self
    if(allocated(self%A)) deallocate(self%A)
  end subroutine ClassMatrix4DFree

  subroutine ClassMatrix4DFinalize( self )
    type(ClassMatrix4D), intent(inout) :: self
    call self%Free()
  end subroutine ClassMatrix4DFinalize

  subroutine ClassMatrix4DSetValue( self, n1, n2, n3, n4, dval )
    class(ClassMatrix4D), intent(inout) :: self
    integer             , intent(in)    :: n1, n2, n3, n4
    real(kind(1d0))     , intent(in)    :: dval
    self%A(n1,n2,n3,n4) = dval
  end subroutine ClassMatrix4DSetValue

  subroutine ClassMatrix4DAddValue( self, n1, n2, n3, n4, dval )
    class(ClassMatrix4D), intent(inout) :: self
    integer             , intent(in)    :: n1, n2, n3, n4
    real(kind(1d0))     , intent(in)    :: dval
    self%A(n1,n2,n3,n4) = self%A(n1,n2,n3,n4) + dval
  end subroutine ClassMatrix4DAddValue

  subroutine ClassMatrix4DAddMatrix( self, Alpha, dMat )
    class(ClassMatrix4D), intent(inout) :: self
    real(kind(1d0))     , intent(in)    :: dMat(:,:,:,:), Alpha
    logical :: SAME
    SAME=size(dMat,1) == size(self%A,1) .and. &
         size(dMat,2) == size(self%A,2) .and. &
         size(dMat,3) == size(self%A,3) .and. &
         size(dMat,4) == size(self%A,4)
    if( .not. SAME )then
       write(*,*) "size 1", size(dMat,1), size(self%A,1) 
       write(*,*) "size 2", size(dMat,2), size(self%A,2) 
       write(*,*) "size 3", size(dMat,3), size(self%A,3) 
       write(*,*) "size 4", size(dMat,4), size(self%A,4) 
       call Assert("In 4DAddMatrix, sizes do not match")
       stop
    endif
    self%A = self%A + dMat * Alpha
  end subroutine ClassMatrix4DAddMatrix

  subroutine ClassMatrix4DSetAllMatrix( self, dval )
    class(ClassMatrix4D), intent(inout) :: self
    real(kind(1d0))     , intent(in)    :: dval
    if(self%TotalSize()<=0)return
    self%A = dval
  end subroutine ClassMatrix4DSetAllMatrix

  subroutine ClassMatrix4DSet4DMat( self, D4Mat )
    class(ClassMatrix4D), intent(inout) :: self
    real(kind(1d0))     , intent(in)    :: D4Mat(:,:,:,:)
    integer :: iDim
    do iDim = 1, 4
       if( size(D4Mat,iDim) /= self%Size( iDim ))then
          call assert( "Incommensurate size in D4Mat assignment" )
          stop
       endif
    enddo
    self%A = D4Mat
  end subroutine ClassMatrix4DSet4DMat

  subroutine ClassMatrix4DGetValue( self, n1, n2, n3, n4, dval )
    class(ClassMatrix4D), intent(inout) :: self
    integer             , intent(in)    :: n1, n2, n3, n4
    real(kind(1d0))     , intent(out)    :: dval
    dval = self%A(n1,n2,n3,n4)
  end subroutine ClassMatrix4DGetValue

  !> Multiplies the matrix by a real number.
  subroutine Double_x_ClassMatrix4D(Matrix,x)
    Class(ClassMatrix4D), intent(inout) :: Matrix
    DoublePrecision   , intent(in)    :: x
    if(.not.allocated(Matrix%A))then
       call Assert("Matrix not initialized")
    else
       Matrix%A=Matrix%A*x
    endif
  end subroutine Double_x_ClassMatrix4D

  !> Writes the relevant information of the matrix class to a unit.
  subroutine WriteClassMatrix4DToUnit(Matrix,OutputUnit)
    Class(ClassMatrix4D), intent(in) :: Matrix
    integer, optional   , intent(in) :: OutputUnit
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
       call WriteClassMatrix4DToFormattedUnit(Matrix,OutUnit)
    case("UNFORMATTED")
       call WriteClassMatrix4DToUnformattedUnit(Matrix,OutUnit)
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
    !
  end subroutine WriteClassMatrix4DToUnit


  subroutine WriteClassMatrix4DToFormattedUnit(self,OutUnit)
    Class(ClassMatrix4D), intent(in) :: self
    integer             , intent(in) :: OutUnit
    integer :: iostat, i, j, k, l, n1, n2, n3, n4
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=*), parameter :: FORMAT_INTS="(*(x,i0))"
    character(len=*), parameter :: FORMAT_MAT="(*"//DOUBLE_PRINT_FORMAT//")"
    !
    !.. Write the Matrix attributes and dimensions
    N1=self%Size(1)
    N2=self%Size(2)
    N3=self%Size(3)
    N4=self%Size(4)
    write(OutUnit      , &
         FORMAT_INTS   , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         N1, N2, N3, N4     
    if(iostat/=0)call Assert(iomsg)
    !
    !.. Write Matrix
    do l=1,N4
       do k=1,N3
          do j=1,N2
             write(OutUnit                    , &
                  FORMAT_MAT                  , &
                  IOSTAT=iostat               , &
                  IOMSG =iomsg                ) &
                  (self%A(i,j,k,l),i=1,N1 )
             if(iostat/=0)call Assert(iomsg)
          enddo
       enddo
    enddo
    !
  end subroutine WriteClassMatrix4DToFormattedUnit


  subroutine WriteClassMatrix4DToUnformattedUnit(self,OutUnit)
    Class(ClassMatrix4D), intent(in) :: self
    integer             , intent(in) :: OutUnit
    integer :: iostat, i, j, k, l, n1, n2, n3, n4
    character(len=IOMSG_LENGTH) :: iomsg
    !
    !.. Write the Matrix attributes and dimensions
    N1=self%Size(1)
    N2=self%Size(2)
    N3=self%Size(3)
    N4=self%Size(4)
    write(OutUnit      , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         N1, N2, N3, N4
    if(iostat/=0)call Assert(iomsg)
    !
    !.. Write Matrix
    write(OutUnit       , &
         IOSTAT=iostat  , &
         IOMSG =iomsg  )  &
         ( ( ( (          &
         self%A(i,j,k,l), &
         i = 1, N1 )    , &
         j = 1, N2 )    , &
         k = 1, N3 )    , &
         l = 1, N4 )
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine WriteClassMatrix4DToUnformattedUnit

  !> Reads the relevant information of the matrix class from a unit.
  subroutine ReadClassMatrix4DFromUnit(self,InputUnit)
    Class(ClassMatrix4D), intent(inout):: self
    integer, optional   , intent(in)   :: InputUnit
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
    call self%Free()
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
       call ReadClassMatrix4DFromFormattedUnit(self,InUnit)
    case("UNFORMATTED")
       call ReadClassMatrix4DFromUnformattedUnit(self,InUnit)
    case DEFAULT
       call Assert("Invalid Input Unit Format")
    end select
    !
  end subroutine ReadClassMatrix4DFromUnit

  subroutine ReadClassMatrix4DFromFormattedUnit(self,InUnit)
    Class(ClassMatrix4D), intent(inout) :: self
    integer             , intent(in)    :: InUnit
    integer :: iostat, i, j, k, l, n1, n2, n3, n4
    character(len=IOMSG_LENGTH) :: iomsg=" "
    character(len=*), parameter :: FORMAT_INTS="(*(x,i))"
    character(len=*), parameter :: FORMAT_MAT="(*"//DOUBLE_PRINT_FORMAT//")"
    !
    !.. Read the Matrix attributes and dimensions
    read(InUnit,*, &
         IOSTAT=iostat   , &
         IOMSG=iomsg     ) &
         N1, N2, N3, N4     
    if(iostat/=0)call Assert(iomsg)
    !
    allocate( self%A( N1, N2, N3, N4 ) )
    self%A=0.d0
    !
    !.. Read Matrix
    do l=1,N4
       do k=1,N3
          do j=1,N2
             read(InUnit,FMT=FORMAT_MAT       , &
                  IOSTAT=iostat               , &
                  IOMSG =iomsg  )               &
                  (self%A(i,j,k,l),i=1,N1)
             if(iostat/=0)call Assert(iomsg)
          enddo
       enddo
    enddo
    !
  end subroutine ReadClassMatrix4DFromFormattedUnit
  !
  subroutine ReadClassMatrix4DFromUnformattedUnit(self,InUnit)
    Class(ClassMatrix4D), intent(inout) :: self
    integer             , intent(in)    :: InUnit
    integer :: iostat, i, j, k, l, N1, N2, N3, N4
    character(len=IOMSG_LENGTH) :: iomsg
    ! 
    !.. Write the Matrix attributes and dimensions
    read(InUnit        , &
         IOSTAT=iostat , &
         IOMSG=iomsg   ) &
         N1, N2, N3, N4     
    if(iostat/=0)call Assert(iomsg)
    !
    allocate( self%A( N1, N2, N3, N4 ) )
    self%A=0.d0
    !
    read(InUnit             , &
         IOSTAT=iostat      , &
         IOMSG =iomsg  )      &
         ((((self%A(i,j,k,l), &
         i = 1, N1 )        , &
         j = 1, N2 )        , &
         k = 1, N3 )        , &
         l = 1, N4 )
    if(iostat/=0)call Assert(iomsg)
    !
  end subroutine ReadClassMatrix4DFromUnformattedUnit

  function ClassMatrix4D_x_ClassMatrix4D_124_123( self, B, vimiA, vimaA, i4miB, i4maB ) result( C )
    !
    Class(ClassMatrix4D), intent(in) :: self
    type (ClassMatrix4D), intent(in) :: B
    type (ClassMatrix)  , pointer    :: C
    integer             , intent(in) :: vimiA(4), vimaA(4), i4miB, i4maB
    !
    integer         :: i1,i2,i3,i4,i4_
    real(kind(1d0)) :: dBuf
    !
    allocate(C)
    call C%InitFull( vimaA(3)-vimiA(3)+1, i4maB-i4miB + 1 )
    C = 0.d0
    !C_{_3_,_4_} = A_{1'2',_3_ 4'} B_{1'2',4' _4_}
    do i4_ = i4miB, i4maB
       do i3 = vimiA(3), vimaA(3)
          dBuf = 0.d0
          do i4 = vimiA(4), vimaA(4)
             do i2 = vimiA(2), vimaA(2)
                do i1 = vimiA(1), vimaA(1)
                   dBuf = dBuf + self%A(i1,i2,i3,i4) * B%A(i1,i2,i4,i4_)
                enddo
             enddo
          enddo
          C%A(i3-vimiA(3)+1,i4_-i4miB+1)=dBuf
       enddo
    enddo
    !
  end function ClassMatrix4D_x_ClassMatrix4D_124_123


  
  !----------------------------------------------------------------
  !.. Spectral analysis

  !> Diagonalization of Class Matrix with D&C, 
  !! assuming the matrix is symmetric
  !! Returns a diagonal matrix with the eigenvalues, and
  !! a matrix with the eigenvectors
  subroutine ClassMatrix_Diagonalize( Matrix, Eval, Evec ) 
    Class(ClassMatrix)          , intent(in)  :: Matrix
    real(kind(1d0)), allocatable, intent(out) :: Eval(:)
    type(ClassMatrix)           , intent(out) :: Evec
    !
    integer                      :: lda, n, i
    integer                      :: lwork, liwork, INFO
    real(kind(1d0)), allocatable :: work(:)
    integer(kind=8), allocatable :: iwork(:)

    n    = Matrix%NRows()
    lda  = Matrix%NRMax - Matrix%NRMin + 1

    Evec = Matrix
    allocate(Eval(n))
    Eval = 0.d0
    allocate( work(1), iwork(1) )
    call DSYEVD( 'V', 'U', n, Evec%A, lda, Eval, work, -1, iwork, -1, INFO )
    lwork = max(int(work(1)),2*n**2+6*n+1)
    if ( (abs(lwork) >= huge(n)) .or. (lwork < 0) ) call Assert( &
         'The lwork variable exceeds the integer kind=4 limits for DSYEVD.' )
    liwork = max(int(iwork(1)),5*n+3)
    deallocate( work, iwork )
    allocate( work(lwork), iwork(liwork) )
    !
    call DSYEVD( 'V', 'U', n, Evec%A, lda, Eval, work, lwork, iwork, liwork, INFO )
    !
    if(info<0)call ErrorMessage(&
         "DSYEVD: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
    if(info>0.and.info<=n)call ErrorMessage(&
         "DSYEVD: failed to converge since "//AlphabeticNumber(info)//" off-diagonal"//&
         "elements of an intermediate tridiagonal did not converge to zero.")
    if(info>n)call ErrorMessage("DSYEVD: failed to compute some eigenvalues")
    !.. 
    do i=1,n
       if(Evec%A(1,i)<0.d0) Evec%A(:,i)=-Evec%A(:,i)
    enddo
    !
    deallocate(work,iwork)
    !
  end subroutine ClassMatrix_Diagonalize

  !> Diagonalization of Class Matrix with D&C, 
  !! assuming the matrix is symmetric, returns a diagonal matrix
  !! with the eigenvalues, and a matrix with the eigenvectors of
  !! a generalized banded matrix secular problem
  subroutine ClassMatrix_Diagonalize_Banded( Matrix, SMatrix, Eval, Evec ) 
    Class(ClassMatrix)          , intent(in)  :: Matrix
    Class(ClassMatrix)          , intent(in)  :: SMatrix
    real(kind(1d0)), allocatable, intent(out) :: Eval(:)
    type(ClassMatrix)           , intent(out) :: Evec
    !
    integer                      :: lda, n, i
    integer                      :: lwork, liwork, INFO
    real(kind(1d0)), allocatable :: work(:)
    integer(kind=8), allocatable :: iwork(:)

!!$******** MUST IMPLEMENT EX NOVO WITH A CALL TO DSBGVD
!!$    
!!$    n    = Matrix%NRows()
!!$    lda  = Matrix%NRMax - Matrix%NRMin + 1
!!$
!!$    Evec = Matrix
!!$    allocate(Eval(n))
!!$    Eval = 0.d0
!!$    allocate( work(1), iwork(1) )
!!$    call DSYEVD( 'V', 'U', n, Evec%A, lda, Eval, work, -1, iwork, -1, INFO )
!!$    lwork = max(int(work(1)),2*n**2+6*n+1)
!!$    if ( (abs(lwork) >= huge(n)) .or. (lwork < 0) ) call Assert( &
!!$         'The lwork variable exceeds the integer kind=4 limits for DSYEVD.' )
!!$    liwork = max(int(iwork(1)),5*n+3)
!!$    deallocate( work, iwork )
!!$    allocate( work(lwork), iwork(liwork) )
!!$    !
!!$    call DSYEVD( 'V', 'U', n, Evec%A, lda, Eval, work, lwork, iwork, liwork, INFO )
!!$    !
!!$    if(info<0)call ErrorMessage(&
!!$         "DSYEVD: The "//AlphabeticNumber(-info)//"-th parameter has an illegal value")
!!$    if(info>0.and.info<=n)call ErrorMessage(&
!!$         "DSYEVD: failed to converge since "//AlphabeticNumber(info)//" off-diagonal"//&
!!$         "elements of an intermediate tridiagonal did not converge to zero.")
!!$    if(info>n)call ErrorMessage("DSYEVD: failed to compute some eigenvalues")
!!$    !.. 
!!$    do i=1,n
!!$       if(Evec%A(1,i)<0.d0) Evec%A(:,i)=-Evec%A(:,i)
!!$    enddo
!!$    !
!!$    deallocate(work,iwork)
    !
  end subroutine ClassMatrix_Diagonalize_Banded

end module ModuleMatrix
