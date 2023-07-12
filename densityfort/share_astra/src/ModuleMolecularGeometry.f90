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
module ModuleMolecularGeometry

  use, intrinsic :: ISO_FORTRAN_ENV

  !L0.0
  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO

  implicit none

  private

  integer, parameter :: MAX_NATOMS = 100

  type, public :: ClassMolecularGeometry
     !
     private
     !
     logical         :: INITIALIZED = .FALSE.
     integer         :: nat = 0
     integer         :: AtNumb(    1:MAX_NATOMS) = 0
     character(len=4):: Symbol(    1:MAX_NATOMS) =" "
     real(kind(1d0)) :: Charge(    1:MAX_NATOMS) = 0.d0
     real(kind(1d0)) :: Coords(1:3,1:MAX_NATOMS) = 0.d0
     !
   contains
     !
     procedure, public  :: free        => ClassMolecularGeometry_Free
     procedure, public  :: ParseFile   => ClassMolecularGeometry_ParseFile
     procedure, public  :: ParseMolden => ClassMolecularGeometry_ParseMoldenFile
     procedure, public  :: show        => ClassMolecularGeometry_Show
     !
     !.. Accessors routines
     procedure, public  :: getNat    => ClassMolecularGeometry_GetNat
     procedure, public  :: getAtNumb => ClassMolecularGeometry_GetAtNumb  
     generic  , public  :: getCharge => ClassMolecularGeometry_GetCharge, ClassMolecularGeometry_GetTotalCharge
     procedure, public  :: getCoords => ClassMolecularGeometry_GetCoords
     !
     !
     procedure, private :: ClassMolecularGeometry_Free
     procedure, private :: ClassMolecularGeometry_ParseFile
     procedure, private :: ClassMolecularGeometry_ParseMoldenFile
     procedure, private :: ClassMolecularGeometry_Show
     procedure, private :: ClassMolecularGeometry_GetNat
     procedure, private :: ClassMolecularGeometry_GetAtNumb
     procedure, private :: ClassMolecularGeometry_GetCharge
     procedure, private :: ClassMolecularGeometry_GetTotalCharge
     procedure, private :: ClassMolecularGeometry_GetCoords
     !
  end type ClassMolecularGeometry
  type(ClassMolecularGeometry), public :: MolecularGeometry
  
!!$  I must extend the parser of parameters to the case of inline list, with and without values
!!$  type(ClassParameterList) :: BasisParameters


contains

  !> Return the total number of atoms 
  integer function ClassMolecularGeometry_GetNat( self ) result( nat )
    class( ClassMolecularGeometry ), intent(in) :: self
    nat = 0
    if(self%INITIALIZED) nat = self%nat
  end function ClassMolecularGeometry_GetNat

  !> Return the atomic number of any given atom
  integer function ClassMolecularGeometry_GetAtNumb( self, iAtom ) result( res )
    class( ClassMolecularGeometry ), intent(in) :: self
    integer                        , intent(in) :: iAtom
    res=0
    if(iAtom<1.or.iAtom>self%nat)return
    if(self%INITIALIZED) res = self%AtNumb(iAtom)
  end function ClassMolecularGeometry_GetAtNumb

  !> Return the atomic number of any given atom
  real(kind(1d0)) function ClassMolecularGeometry_GetCharge( self, iAtom ) result( res )
    class( ClassMolecularGeometry ), intent(in) :: self
    integer                        , intent(in) :: iAtom
    res=0.d0
    if(iAtom<1.or.iAtom>self%nat)return
    if(self%INITIALIZED) res = self%Charge(iAtom)
  end function ClassMolecularGeometry_GetCharge

  !> Return the atomic number of any given atom
  real(kind(1d0)) function ClassMolecularGeometry_GetTotalCharge( self ) result( res )
    class( ClassMolecularGeometry ), intent(in) :: self
    integer :: iAtom
    res=0.d0
    if(.not.self%INITIALIZED)return
    do iAtom=1,self%nat
       res = res + self%Charge(iAtom)
    end do
  end function ClassMolecularGeometry_GetTotalCharge

  !> Return the atomic number of any given atom
  subroutine ClassMolecularGeometry_GetCoords( self, iAtom, xv )
    class( ClassMolecularGeometry ), intent(in)  :: self
    integer                        , intent(in)  :: iAtom
    real(kind(1d0))                , intent(out) :: xv(3)
    xv=0.d0
    if(iAtom<1.or.iAtom>self%nat)return
    if(self%INITIALIZED) xv(1:3) = self%Coords(1:3,iAtom)
  end subroutine ClassMolecularGeometry_GetCoords


  !.. Free the molecular geometry 
  subroutine ClassMolecularGeometry_Free( self )
    class( ClassMolecularGeometry ), intent(inout) :: self
    self%nat    = 0
    self%AtNumb = 0
    self%Symbol = " "
    self%Charge = 0.d0
    self%Coords = 0.d0
    self%INITIALIZED = .FALSE.
  end subroutine ClassMolecularGeometry_Free


  !.. Print the molecular geometry on standard output or 
  !   on an alternative optional formatted unit
  subroutine ClassMolecularGeometry_Show( self, unit_ )
    class( ClassMolecularGeometry ), intent(inout) :: self
    integer, optional              , intent(in)    :: unit_
    integer :: unit, iAtom
    if( .not. self%INITIALIZED )return
    unit = OUTPUT_UNIT 
    if( present( unit_ ) ) unit = unit_
    write(unit,"(a)") "Molecular Geometry"
    do iAtom = 1, self%nat
       write(unit,"(i3,x,a,x,i3,x,f4.0,3(x,f12.6))") &
            iAtom , &
            self%Symbol(  iAtom), &
            self%AtNumb(  iAtom), &
            self%Charge(  iAtom), &
            self%Coords(1,iAtom), &
            self%Coords(2,iAtom), &
            self%Coords(3,iAtom)
    enddo
  end subroutine ClassMolecularGeometry_Show


  !.. Parse a file with the description of the molecular geometry
  subroutine ClassMolecularGeometry_ParseFile( self, FileName )

    class( ClassMolecularGeometry ), intent(inout) :: self
    character(len=*)               , intent(in)    :: FileName

    integer             :: uid, iostat, ich
    character(len=50)   :: skey
    character(len=1000) :: line, line2
    integer :: nAtomTypes, nAtoms, iAtomType, nSameAtoms, iAtom, iCoord
    real(kind(1d0)) :: Charge

    call self%free()
    call OpenFile( FileName, uid, "read", "formatted" )

    !.. Skip the header of MOLECULE.INP
    do 
       read(uid,"(a)",iostat=iostat) line
       if( iostat /= 0 )then
          call Assert("Error Reading "//FileName//". Stop Forced")
          stop
       endif
       ich=index(trim(line),"Atomtypes")
       if(ich>0)exit
    enddo
    
    !.. Parse the parameter line
    
    !.. Atomtypes
    skey="Atomtypes="
    ich=index(trim(line),trim(skey))
    if(ich>0)then
       line2=adjustl(line(ich+len_trim(skey):))
       ich=index(line2," ")
       read(line2(1:ich),*) nAtomTypes
    endif

    write(*,*) "nAtomTypes =",nAtomTypes

    !.. Generators
    ich=index(trim(line),"Generators")
    if(ich>0)then
       !***
       !    DO NOTHING AT THE MOMENT, AS WE HAVEN'T INCLUDED SYMMETRY YET
       !***
    endif

    !.. Nosymmetry
    ich=index(trim(line),"Nosymmetry")
    if(ich>0)then
       !***
       !    DO NOTHING AT THE MOMENT, AS WE HAVEN'T INCLUDED SYMMETRY YET
       !***
    endif
    
    nAtoms=0
    !
    do iAtomType = 1, nAtomTypes

       call FetchLineStrn( uid, line )

       ich=index(trim(line),"Charge=")
       if(ich<0) call Assert("Atomic Charge missing in Molecular-Geometry File")
       read(line(ich+len("Charge="):),*)Charge

       ich=index(trim(line),"Atoms=")
       if(ich<0) call Assert("Number of atoms missing in Molecular-Geometry File")
       read(line(ich+len("Atoms="):),*)nSameAtoms
       
       do iAtom = 1, nSameAtoms
          
          call FetchLineStrn( uid, line )
          line=adjustl(line)
          ich=index(line," ")
          if(ich<0) call Assert("Atomic parameters missing in Molecular-Geometry File")
          
          nAtoms = nAtoms + 1
          
          self%nat  =  nAtoms
          self%AtNumb( nAtoms ) = int( Charge + 0.1d0 )
          self%Symbol( nAtoms ) = line(:ich-1)
          self%Charge( nAtoms ) = Charge
          read(line(ich+1:),*) (self%Coords( iCoord, nAtoms ), iCoord = 1, 3 )

       enddo
       
    end do
    
    self%INITIALIZED = .TRUE.
    
  end subroutine ClassMolecularGeometry_ParseFile


  !.. Same as ClassMolecularGeometry_ParseFile but for the molden format file.
  subroutine ClassMolecularGeometry_ParseMoldenFile( self, FileName )

    class( ClassMolecularGeometry ), intent(inout) :: self
    character(len=*)               , intent(in)    :: FileName

    integer             :: uid, iostat, ich
    character(len=50)   :: skey
    character(len=1000) :: line, line2
    integer :: nAtomTypes, nAtoms, iAtomType, nSameAtoms, iAtom, iCoord, iCharge
    real(kind(1d0)) :: Charge, xCoord, yCoord, zCoord
    CHARACTER*6     :: name

    call self%free()
    call OpenFile( FileName, uid, "read", "formatted" )

    !.. Skip the header of MOLECULE.INP
    do 
       read(uid,"(a)",iostat=iostat) line
       if( iostat /= 0 )then
          call Assert("Error Reading "//FileName//". [Atoms] not found. Stop Forced")
          stop
       endif
       ich=index(trim(line),"[Atoms]")
       if(ich>0)exit
    enddo
    
    nAtoms=0
    do 
       read(uid,"(a)",iostat=iostat) line
       if( iostat /= 0 )then
          call Assert("Error Reading "//FileName//". [Atoms] not found. Stop Forced")
          stop
       endif

       ich=index(trim(line),"[")
       if(ich>0)exit
       !.. This format matchs the one used in DALTON. It is abamolden.F file which writes the molden file.
       read ( line, '(A,1X,I5,1X,I5,3(1X,F20.10))') &
       &  name, iAtom, icharge, xCoord, yCoord, zCoord
       
       nAtoms = nAtoms + 1

       self%nat  =     nAtoms
       self%AtNumb(    nAtoms ) = iCharge
       self%Symbol(    nAtoms ) = name(:4)
       self%Charge(    nAtoms ) = dble(iCharge)
       self%Coords( 1, nAtoms ) = xCoord
       self%Coords( 2, nAtoms ) = yCoord
       self%Coords( 3, nAtoms ) = zCoord

       ! write(*,*) self%AtNumb(    nAtoms ) 
       ! write(*,*) self%Symbol(    nAtoms )
       ! write(*,*) self%Charge(    nAtoms )
       ! write(*,*) self%Coords( 1, nAtoms )
       ! write(*,*) self%Coords( 2, nAtoms )
       ! write(*,*) self%Coords( 3, nAtoms )
        
       
    enddo
    
    self%INITIALIZED = .TRUE.
    
  end subroutine ClassMolecularGeometry_ParseMoldenFile

end module ModuleMolecularGeometry
