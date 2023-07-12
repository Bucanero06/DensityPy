! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
!>  Devises a series of transformations of the close-coupling
!!  basis that are necessary for the spectral analysis and
!!  time propagation of the system:
!!  1. elimination of linear dependencies in the cc active sector
!!  2. orthonormalization of the outer functions
!!  3. diagonalization of the internal and external-spherical
!!     diagonal blocks
!!  Point 1 is specific to a given selection of ions in the
!!  close coupling space. It gives rise to a basis that mixes
!!  the ions together, and which may be smaller than the original
!!
program ConditionSymESpace

  use ModuleString
  use ModuleMatrix
  use ModuleGroups
  use ModuleAstraConfigFile
  use ModuleXlm
  use ModuleSymESpace
  use ModuleESpace
  use ModuleMainInterface

  implicit none

  character(len=:), allocatable :: AstraConfigFile
  integer                       :: Multiplicity
  character(len=5), allocatable :: vSymLabel(:)

  character(len=:), allocatable :: ccConfigFile, RootDir
  type(ClassESpace)             :: Space
  type(ClasssymCCSpace), pointer :: SymSpace
  integer                       :: iSym

  call GetRunTimeParameters( AstraConfigFile, Multiplicity, vSymLabel )
  call ParseAstraConfigFile( AstraConfigFile, ccConfigFile, RootDir   )

  call Space%SetMultiplicity( Multiplicity )
  call Space%SetRootDir     ( RootDir      )
  call Space%ParseConfigFile( ccConfigFile )

  call GlobalGroup%Init(Space%GetGroupName())
  call GlobalXlmSet%Init(GlobalGroup,Space%GetLmax())

  call Execute_Command_Line("mkdir -p test/")
  do iSym = 1, size( vSymLabel )

     call Space%CheckSymmetry( trim(vSymLabel(iSym)) )
     SymSpace => Space%GetSymElectSpace( trim(vSymLabel(iSym)) ) 

     call Execute_Command_Line("Run STEX program with symmetry Sym and spit result in test/Hmat and test/Smat")
     call TestMatrix(SymSpace,"H","test/Hmat")
     call TestMatrix(SymSpace,"S","test/Smat")

  enddo

contains

  subroutine TestMatrix(SymSpace,OpLabel,STEXFileName)
    type(ClasssymCCSpace), intent(inout) :: SymSpace
    character(len=*)     , intent(in)    :: OpLabel, STEXFileName
    type(ClassMatrix) :: Mat, Mat_stex
    real(kind(1d0)), parameter :: THRESHOLD = 1.d-10
    real(kind(1d0)), parameter :: ELEMENT_THRESHOLD = 1.d-14
    real(kind(1d0)) :: Astra_Value, STEX_Value
    integer         :: i, j
    call AssembleMatrix(SymSpace,OpLabel,Mat)
    call Mat_stex%Read(STEXFileName)
    call Mat_stex%Multiply(-1.d0)
    call Mat%Add(Mat_stex)
    write(*,"(a)",advance="no") "["//OpLabel//"] "
    if(Mat%Norm1() > THRESHOLD)then
       write(*,*) "FAIL"
       !..Check matrix elements individually and print out
       !  data of those that differ
    else
       write(*,*) "pass"
    endif
    do i = 1, Mat%NRows()
       do j = 1, Mat%NColumns()
          Astra_Value = Mat%Element(i,j)
          STEX_Value = -Mat_stex%Element(i,j) 
          if( Astra_Value - STEX_Value > ELEMENT_THRESHOLD )then
             write(*,*) " Element i,j=",i,j," differ: astra Hij=",Astra_Value, " STEX Hij =",STEX_Value
          endif
       enddo
    enddo
  end subroutine TestMatrix
  
  subroutine AssembleMatrix(SymSpace,OpLabel,Mat)
    type(ClasssymCCSpace), intent(inout) :: SymSpace
    character(len=*)     , intent(in)    :: OpLabel
    type(ClassMatrix)    , intent(out)   :: Mat
    type(ClassSESSESBlock) :: Blk
    logical                :: lTasks(N_SESSES_ID) = .FALSE.
    lTasks(SESSES_LOAD)=.TRUE.
    call Blk%Free()
    call Blk%Init(SymSpace,SymSpace,OpLabel,SymSpace%GetRoot())
    call Blk%Driver(lTasks)
    call Blk%Assemble(Mat)
  end subroutine AssembleMatrix
  
  
end program ConditionSymESpace
