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
  real(kind(1d0))               :: threshold

  character(len=:), allocatable :: ccConfigFile, RootDir
  type(ClassESpace)             :: Space
  type(ClasssymCCSpace), pointer :: SymSpace
  integer                       :: iSym

  call GetRunTimeParameters( AstraConfigFile, Multiplicity, vSymLabel, threshold )
  call ParseAstraConfigFile( AstraConfigFile, ccConfigFile, RootDir   )

  call Space%SetMultiplicity( Multiplicity )
  call Space%SetRootDir     ( RootDir      )
  call Space%ParseConfigFile( ccConfigFile )

  call GlobalGroup%Init(Space%GetGroupName())
  call GlobalXlmSet%Init(GlobalGroup,Space%GetLmax())

  do iSym = 1, size( vSymLabel )

     call Space%CheckSymmetry( trim(vSymLabel(iSym)) )
     SymSpace => Space%GetSymElectSpace( trim(vSymLabel(iSym)) ) 
     call SymSpace%Condition(threshold)
     
  enddo

end program ConditionSymESpace
