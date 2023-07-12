module ModuleSymmetricLocalizedStates
  
  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleGroups
  use ModuleSymmetricElectronicChannel
  use ModuleGeneralChannel
  
  implicit none

  private


  !> States which are not associated to any given
  !! parent ion, and which are expressed entirely
  !! in terms of molecular orbitals. Therefore, they 
  !! do not give rise to any matrix element with 
  !! partial wave channels. Examples: core-hole states
  !! and multiply-excited configurations to describe 
  !! short-range correlation.
  type, public :: ClassSymmetricLocalizedStates
     ! {{{ private attributes
     private
     integer                       :: Multiplicity
     type(ClassIrrep), pointer     :: Irrep
     integer                       :: NStates
     character(len=:), allocatable :: StorageDir
     ! }}}
   contains
     generic, public :: SetIrrep            => ClassSymmetricLocalizedStatesSetIrrep
     generic, public :: SetMultiplicity     => ClassSymmetricLocalizedStatesSetMultiplicity
     generic, public :: SetNStates          => ClassSymmetricLocalizedStatesSetNStates
     generic, public :: SetStorageDir       => ClassSymmetricLocalizedStatesSetStorageDir
     generic, public :: GetStorageDir       => ClassSymmetricLocalizedStatesGetStorageDir
     generic, public :: GetSymLabel         => ClassSymmetricLocalizedStatesGetSymLabel
     generic, public :: GetPILabel          => ClassSymmetricLocalizedStatesGetPILabel
     generic, public :: GetDiffuseIrrepName => ClassSymmetricLocalizedStatesGetDiffuseIrrepName
     generic, public :: show => ClassSymmetricLocalizedStatesShow
     generic, public :: free => ClassSymmetricLocalizedStatesFree
     ! {{{ private procedures

     !generic :: assignment(=) => &
     !        CopySLStoSLS

     procedure, private :: ClassSymmetricLocalizedStatesSetIrrep
     procedure, private :: ClassSymmetricLocalizedStatesSetMultiplicity
     procedure, private :: ClassSymmetricLocalizedStatesSetNStates
     procedure, private :: ClassSymmetricLocalizedStatesSetStorageDir
     procedure, private :: ClassSymmetricLocalizedStatesGetStorageDir
     procedure, private :: ClassSymmetricLocalizedStatesGetSymLabel
     procedure, private :: ClassSymmetricLocalizedStatesGetPILabel
     procedure, private :: ClassSymmetricLocalizedStatesGetDiffuseIrrepName
     procedure, private :: ClassSymmetricLocalizedStatesShow
     procedure, private :: ClassSymmetricLocalizedStatesFree
     !procedure, private :: CopySLStoSLS 
     final              :: ClassSymmetricLocalizedStatesFinal

     ! }}} 
  end type ClassSymmetricLocalizedStates


contains  
  

  subroutine ClassSymmetricLocalizedStatesShow( LS, unit )
    class(ClassSymmetricLocalizedStates), intent(in) :: LS
    integer, optional                   , intent(in) :: unit
    integer :: outunit
    outunit = OUTPUT_UNIT
    if(present(unit)) outunit = unit
    write(outunit,"(a)",advance="no") "  Localized States Multiplicity :"
    write(outunit,"(i0)") LS%Multiplicity
    if ( associated(LS%irrep) ) then
       write(outunit,"(a)",advance="no") "  Localized States Symmetry :"
       call LS%irrep%show(unit)
    end if
    write(outunit,"(a,i4)"          ) "  Number Localized States  :",LS%NStates
  end subroutine ClassSymmetricLocalizedStatesShow


  subroutine ClassSymmetricLocalizedStatesFree( LS )
    class(ClassSymmetricLocalizedStates), intent(inout) :: LS
    LS%Multiplicity = 0
    LS%Irrep   => NULL()
    LS%NStates =  0
    if ( allocated(LS%StorageDir) ) deallocate( LS%StorageDir )
  end subroutine ClassSymmetricLocalizedStatesFree

  
  subroutine ClassSymmetricLocalizedStatesFinal( LS )
    type(ClassSymmetricLocalizedStates) :: LS
    call LS%free()
  end subroutine ClassSymmetricLocalizedStatesFinal


  subroutine  ClassSymmetricLocalizedStatesSetIrrep( LS, Irrep )
    class(ClassSymmetricLocalizedStates), intent(inout) :: LS
    class(ClassIrrep), target,            intent(in)    :: Irrep
    if ( associated(LS%Irrep) ) LS%Irrep => NULL()
    LS%Irrep => Irrep
  end subroutine ClassSymmetricLocalizedStatesSetIrrep



  subroutine  ClassSymmetricLocalizedStatesSetNStates( LS, N )
    class(ClassSymmetricLocalizedStates), intent(inout) :: LS
    integer,                              intent(in)    :: N
    LS%NStates = N
  end subroutine ClassSymmetricLocalizedStatesSetNStates



  subroutine ClassSymmetricLocalizedStatesSetStorageDir( LS, InDir )
    class(ClassSymmetricLocalizedStates), intent(inout) :: LS
    character(len=*)                    , intent(in)    :: InDir
    if ( allocated(LS%StorageDir) ) deallocate( LS%StorageDir )
    allocate( LS%StorageDir, source = InDir )
  end subroutine ClassSymmetricLocalizedStatesSetStorageDir



  subroutine ClassSymmetricLocalizedStatesSetMultiplicity( LS, Multiplicity )
    class(ClassSymmetricLocalizedStates), intent(inout) :: LS
    integer                             , intent(in)    :: Multiplicity
    LS%Multiplicity = Multiplicity
  end subroutine ClassSymmetricLocalizedStatesSetMultiplicity



  function ClassSymmetricLocalizedStatesGetStorageDir( LS ) result(Dir)
    class(ClassSymmetricLocalizedStates), intent(in) :: LS
    character(len=:), allocatable :: Dir
    if ( .not.allocated(LS%StorageDir) ) call Assert( &
         'The storage dir is not allocated for Loc States, impossible to fetch.' )
    allocate( Dir, source = LS%StorageDir )
  end function ClassSymmetricLocalizedStatesGetStorageDir



  function ClassSymmetricLocalizedStatesGetSymLabel( LS ) result(SymLabel)
    class(ClassSymmetricLocalizedStates), intent(in) :: LS
    character(len=:), allocatable :: SymLabel
    allocate( SymLabel, source = AlphabeticNumber(LS%Multiplicity)//LS%Irrep%GetName() )
  end function ClassSymmetricLocalizedStatesGetSymLabel


  !> The purpose is to be symmetric with respect to the other calls for
  !! SRC and PWC channels.
  function ClassSymmetricLocalizedStatesGetPILabel( LS ) result(PILabel)
    class(ClassSymmetricLocalizedStates), intent(in) :: LS
    character(len=:), allocatable :: PILabel
    allocate( PILabel, source = 'Loc' )
  end function ClassSymmetricLocalizedStatesGetPILabel


  !> The purpose is to be symmetric with respect to the other calls for
  !! SRC and PWC channels.
  function ClassSymmetricLocalizedStatesGetDiffuseIrrepName( LS ) result(DLabel)
    class(ClassSymmetricLocalizedStates), intent(in) :: LS
    character(len=:), allocatable :: DLabel
    allocate( DLabel, source = 'DummyString' )
  end function ClassSymmetricLocalizedStatesGetDiffuseIrrepName

end module ModuleSymmetricLocalizedStates
