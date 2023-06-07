!> Define the classes for the individual parent ions, for the set of
!! parent ions and for the multipoles between them.
!! Furthermore, it defines a global variable for the complete set of
!! parent ions as well as
!! - a configuration file for the parent ions
!! - the files that contain information on individual parent ions (from interface)
!! - the files with the mutipoles (from interface)
module ModuleParentIons

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString
  use ModuleGroups
  use ModuleSymmetryAdaptedSphericalHarmonics
  use ModuleIO

  implicit none

  private

  character(len=*), parameter :: PARENTION_DIR_NAME  = "ParentIons"
  character(len=*), parameter :: PIEnergyLabel       = "En_"

  character(len=*), parameter :: MULTIPOLE_ROOT_NAME = "CartesianMultipoles_"
  character(len=*), parameter :: MultipolesLabel = "multipoles.dat"
  character(len=*), parameter :: MultLabelStore  = "mult"
  character(len=*), parameter :: KinEnergyLabelStore  = "kin"
  character(len=*), parameter :: MULTIPOLE_IO_FORMAT = "(x,d24.16)"


  !> The parent ions are identified by 
  !! - point group (which must be set once and for all at the very beginning)
  !! - multiplicity (if the spin is well defined, multiplicity >=1, &
  !!                 otherwise multiplicity =0)
  !! - irreducible representation, in the order provided by the list of the 
  !!   group class
  !! - number, which is assigned externally by the quantum chemistry program
  !!   and which does not necessarily reflect the energy order of the states
  !!   (due to the fact that the group used here may be just a sub-group of
  !!   the real molecular group, states with the same symmetry could actually
  !!   cross, thus altering the order of their energies).
  type, public :: ClassParentIon
     ! {{{ private attributes

     private
     real(kind(1d0))           :: Energy
     integer                   :: Charge
     integer                   :: Multiplicity
     type(ClassIrrep), pointer :: Irrep
     integer                   :: N
     !> Number of electrons in the parent ion.
     integer, public           :: NePI

     ! }}}
   contains
     generic, public :: init            => ClassParentIonInit
     generic, public :: free            => ClassParentIonFree
     generic, public :: show            => ClassParentIonShow
     generic, public :: GetEnergy       => ClassParentIonGetEnergy
     generic, public :: SetEnergy       => ClassParentIonSetEnergy
     generic, public :: ReadEnergy      => ClassParentIonReadEnergy
     generic, public :: GetCharge       => ClassParentIonGetCharge
     generic, public :: GetMultiplicity => ClassParentIonGetMultiplicity
     generic, public :: GetIrrep        => ClassParentIonGetIrrep
     generic, public :: GetNumber       => ClassParentIonGetNumber
     generic, public :: GetNelect       => ClassParentIonGetNelect
     generic, public :: SetNelect       => ClassParentIonSetNelect
     generic, public :: GetLabel        => ClassParentIonGetLabel
     generic, public :: LoadData        => ClassParentIonLoadData, ClassParentIonLoadAllDataFromUnit
     generic, public :: SaveData        => ClassParentIonSaveData, ClassParentIonSaveAllDataToUnit
     ! {{{ private procedures

     procedure, private :: ClassParentIonInit
     procedure, private :: ClassParentIonFree
     procedure, private :: ClassParentIonShow
     procedure, private :: ClassParentIonGetEnergy
     procedure, private :: ClassParentIonSetEnergy
     procedure, private :: ClassParentIonReadEnergy
     procedure, private :: ClassParentIonGetCharge
     procedure, private :: ClassParentIonGetMultiplicity
     procedure, private :: ClassParentIonGetIrrep
     procedure, private :: ClassParentIonGetNumber
     procedure, private :: ClassParentIonGetNelect
     procedure, private :: ClassParentIonSetNelect
     procedure, private :: ClassParentIonGetLabel
     procedure, private :: ClassParentIonLoadData
     procedure, private :: ClassParentIonLoadAllDataFromUnit
     procedure, private :: ClassParentIonSaveData
     procedure, private :: ClassParentIonSaveAllDataToUnit
     final :: ClassParentIonFinal

     ! }}}
  end type ClassParentIon


  type, private :: ClassComplexVec
     complex(kind(1d0)), allocatable  :: v(:)
  end type ClassComplexVec

  !> Value of the multipoles with a given angular momentum
  !! and symmetry, between two fixed parent ions.
  type, public :: ClassMultipole
     ! {{{ private attributes

     private
     integer                            :: Lmax
     type(ClassParentIon), pointer      :: PionBra
     type(ClassParentIon), pointer      :: PionKet
     type(ClassIrrep), pointer          :: Irrep
     character(len=:), allocatable      :: Dir
     character(len=:), allocatable      :: File
     character(len=:), allocatable      :: CartFile
     integer                            :: Nmultipoles
     type(ClassCartesianSymmetricSet)   :: CartesianSet
     real(kind(1d0)), allocatable       :: CartesianMultipole(:)
     type(ClassXlmSymmetricSet)         :: XlmSymSet
     !> XlmMultipole(0:LMAX).v(1:...)
     type(ClassComplexVec), allocatable :: XlmMultipole(:)
     !
     !.. Comment: the use of objects for seemingly simple variables
     !   like the Xlm multipoles or the cartesian multipoles is 
     !   justified on the basis that the conversion between them
     !   ought to be well defined and managed by a separate module.
     !
     ! }}}
   contains
     generic, public :: setLMax                       => ClassMultipoleSetLMax
     generic, public :: setIrrep                      => ClassMultipoleSetIrrep
     generic, public :: setDir                        => ClassMultipoleSetDir
     generic, public :: setFile                       => ClassMultipoleSetFile
     generic, public :: setCartFile                   => ClassMultipoleSetCartFile
     generic, public :: setCartEl                     => ClassMultipoleSetCartEl
     ! Return " x  xxx  xyy  xzz " etc.
     generic, public :: FetchCartesianStrn            => ClassMultipoleFetchCartesianStrn
     generic, public :: setEmpty                      => ClassMultipoleSetEmpty
     generic, public :: Init                          => ClassMultipoleInit
     generic, public :: LoadCartesianData             => ClassMultipoleLoadCartesianData
     generic, public :: GetN                          => ClassMultipoleGetN
     generic, public :: GetMindex                     => ClassMultipoleGetMindex
     generic, public :: GetMultipole                  =>  ClassMultipoleGetXlmMultipole 
     generic, public :: SetCartesianValue             => ClassMultipoleSetCartesianValue
     generic, public :: GetCartesianValue             => ClassMultipoleGetCartesianValue
     generic, public :: FetchCartesianMultipole       => ClassMultipoleFetchCartesianMultipole
     generic, public :: show                          => ClassMultipoleShow
     generic, public :: ConvertCartesianToSpherical   => ClassMultipoleConvertCartesianToSpherical
     generic, public :: Save                          => ClassMultipoleSave
     generic, public :: Read                          => ClassMultipoleRead
     generic, private :: FileName                     => ClassMultipoleFileName
     generic, public :: free                          => ClassMultipoleFree
     ! {{{ private procedures

     procedure, private :: ClassMultipoleFileName
     procedure, private :: ClassMultipoleSetLMax
     procedure, private :: ClassMultipoleSetIrrep
     procedure, private :: ClassMultipoleSetDir
     procedure, private :: ClassMultipoleSetFile
     procedure, private :: ClassMultipoleSetCartFile
     procedure, private :: ClassMultipoleSetEmpty
     procedure, private :: ClassMultipoleInit
     procedure, private :: ClassMultipoleLoadCartesianData
     procedure, private :: ClassMultipoleConvertCartesianToSpherical
     procedure, private :: ClassMultipoleGetN
     procedure, private :: ClassMultipoleGetMindex
     procedure, private :: ClassMultipoleGetXlmMultipole
     procedure, private :: ClassMultipoleSetCartesianValue
     procedure, private :: ClassMultipoleGetCartesianValue
     procedure, private :: ClassMultipoleSetCartEl
     procedure, private :: ClassMultipoleFetchCartesianStrn
     procedure, private :: ClassMultipoleFetchCartesianMultipole
     procedure, private :: ClassMultipoleSave
     procedure, private :: ClassMultipoleRead
     procedure, private :: ClassMultipoleShow
     procedure, private :: ClassMultipoleFree
     final :: ClassMultipoleFinal

     ! }}}
  end type ClassMultipole


  !> Value of the multipoles with a given angular momentum
  !! and symmetry, between two fixed parent ions.
  type, public :: ClassKineticEnergy
     ! {{{ private attributes

     private
     type(ClassParentIon), pointer      :: PionBra
     type(ClassParentIon), pointer      :: PionKet
     character(len=:), allocatable      :: Dir
     character(len=:), allocatable      :: File
     real(kind(1d0))                    :: KinEnergy = 0.d0
     !
     ! }}}
   contains
     generic, public :: setDir         => ClassKineticEnergySetDir
     generic, public :: setFile        => ClassKineticEnergySetFile
     generic, public :: Init           => ClassKineticEnergyInit
     generic, public :: GetKinEnergy   => ClassKineticEnergyGetKinEnergy
     generic, public :: FetchKinEnergy => ClassKineticEnergyFetchKinEnergy
     generic, public :: show           => ClassKineticEnergyShow
     generic, public :: Save           => ClassKineticEnergySave
     generic, public :: Read           => ClassKineticEnergyRead
     generic, private :: FileName      => ClassKineticEnergyFileName
     generic, public :: free           => ClassKineticEnergyFree
     ! {{{ private procedures

     procedure, private :: ClassKineticEnergySetDir
     procedure, private :: ClassKineticEnergySetFile
     procedure, private :: ClassKineticEnergyInit
     procedure, private :: ClassKineticEnergyGetKinEnergy
     procedure, private :: ClassKineticEnergyFetchKinEnergy
     procedure, private :: ClassKineticEnergyShow
     procedure, private :: ClassKineticEnergySave
     procedure, private :: ClassKineticEnergyRead
     procedure, private :: ClassKineticEnergyFileName
     procedure, private :: ClassKineticEnergyFree
     final :: ClassKineticEnergyFinal

     ! }}}
  end type ClassKineticEnergy

  public :: ParsePionLabel


contains


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !         Procedures   ClassParentIon
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine ClassParentIonInit( ParentIon, Group, Label, iCharge, Energy )
    class(ClassParentIon), intent(inout) :: ParentIon
    type(ClassGroup)     , intent(in)    :: Group
    character(len=*)     , intent(in)    :: Label
    integer              , intent(in)    :: iCharge
    real(kind(1d0)), optional, intent(in) :: Energy
    !
    integer :: Multiplicity
    character(len=8) :: IrrepLabel
    integer :: PionN
    if(present(Energy)) ParentIon.Energy=Energy
    ParentIon.Charge=iCharge
    call ParsePionLabel(Label,Multiplicity,IrrepLabel,PionN)
    if( Group.definesIrrep( trim(IrrepLabel) ) )then
       ParentIon.Multiplicity =  Multiplicity
       ParentIon.Irrep        => Group.GetIrrep( trim(IrrepLabel) )
       ParentIon.N            =  PionN
    else
       call ErrorMessage("Invalid irrep "//trim(IrrepLabel)//" for group "//Group.GetName())
       stop
    endif
  end subroutine ClassParentIonInit



  subroutine ParsePionLabel( Label, Multiplicity, Irrep, N )
    character(len=*) , intent(in)  :: Label
    integer          , intent(out) :: Multiplicity
    character(len=*) , intent(out) :: Irrep
    integer          , intent(out) :: N
    !
    integer :: iChar, iPos
    character(len=16) :: strn
    character(len=:),allocatable :: tmpLabel
    !
    allocate(tmpLabel,source=trim(adjustl(Label)))
    strn=" "
    iChar=1
    do 
       if(iChar>len_trim(tmpLabel))exit
       iPos=index("1234567890",tmpLabel(iChar:iChar))
       if(iPos<=0)exit
       iChar=iChar+1
    enddo
    iChar=iChar-1
    if(iChar<=0)call Assert("Missing multiplicity of Parent ion")

    strn=tmpLabel(1:iChar)
    read(strn,*)Multiplicity
    tmpLabel=tmpLabel(iChar+1:)
    !
    iPos=index(tmpLabel,".")
    Irrep=tmpLabel(1:iPos-1)
    strn=tmpLabel(iPos+1:)
    read(strn,*)N
    !
  end subroutine ParsePionLabel



  subroutine ClassParentIonShow( ParentIon, unit )
    class(ClassParentIon), intent(in) :: ParentIon
    integer, optional    , intent(in) :: unit
    integer :: outunit
    character(len=8) :: strn
    outunit=OUTPUT_UNIT
    if(present(unit))outunit=unit
    write(outunit,"(a)") "Parent Ion info: "
    strn=ParentIon.GetLabel()
    write(outunit,"(a)"      ,advance="no")"  "//strn
    write(outunit,"(a,1x,i2)",advance="no")"  Ne Parent Ion =",ParentIon.NePI 
    write(outunit,"(a,1x,i2)",advance="no")"  Mult =",ParentIon.GetMultiplicity() 
    write(outunit,"(a,1x,i2)",advance="no")"  Q =",ParentIon.GetCharge() 
    write(outunit,"(a,1x,d14.6)")          "  E =",ParentIon.GetEnergy()
  end subroutine ClassParentIonShow


  function ClassParentIonGetLabel( ParentIon ) result( label )
    class(ClassParentIon), intent(in) :: ParentIon
    character(len=:)    , allocatable :: label
    character(len=16) :: mstrn,nstrn,istrn
    write(mstrn,*)ParentIon.GetMultiplicity()
    mstrn=adjustl(mstrn)
    istrn=adjustl(ParentIon.irrep.GetName())
    write(nstrn,*)ParentIon.GetNumber()
    nstrn=adjustl(nstrn)
    allocate(label,source=trim(mstrn)//trim(istrn)//"."//trim(nstrn))
  end function ClassParentIonGetLabel

    
  function ClassParentIonGetMultiplicity( ParentIon ) result( multiplicity )
    class(ClassParentIon), intent(in) :: ParentIon
    integer                           :: multiplicity
    multiplicity = ParentIon.multiplicity
  end function ClassParentIonGetMultiplicity


  function ClassParentIonGetIrrep( ParentIon ) result( irrep )
    class(ClassParentIon), intent(in) :: ParentIon
    type(ClassIrrep), pointer         :: irrep
    irrep => ParentIon.irrep
  end function ClassParentIonGetIrrep


  function ClassParentIonGetNumber( ParentIon ) result( number )
    class(ClassParentIon), intent(in) :: ParentIon
    integer                           :: number
    number = ParentIon.N
  end function ClassParentIonGetNumber


  function ClassParentIonGetNelect( ParentIon ) result( number )
    class(ClassParentIon), intent(in) :: ParentIon
    integer                           :: number
    number = ParentIon.NePI
  end function ClassParentIonGetNelect


  subroutine ClassParentIonSetNelect( ParentIon, Ne )
    class(ClassParentIon), intent(inout) :: ParentIon
    integer              , intent(in)    :: Ne
    ParentIon.NePI = Ne
  end subroutine ClassParentIonSetNelect


  function ClassParentIonGetEnergy( ParentIon ) result( energy )
    class(ClassParentIon), intent(in) :: ParentIon
    real(kind(1d0))                   :: energy
    energy = ParentIon.Energy
  end function ClassParentIonGetEnergy


  subroutine ClassParentIonSetEnergy( ParentIon, Energy )
    class(ClassParentIon), intent(inout) :: ParentIon
    real(kind(1d0))      , intent(in)    :: Energy
    ParentIon.Energy = Energy
  end subroutine ClassParentIonSetEnergy


  function ClassParentIonReadEnergy( PIon, Dir ) result( energy )
    class(ClassParentIon), intent(in) :: PIon
    character(len=*)     , intent(in) :: Dir
    real(kind(1d0))                   :: energy
    character(len=:), allocatable :: FileName, StorageDir
    integer :: uid
    !
    allocate( StorageDir, source=AddSlash(dir)//AddSlash(PARENTION_DIR_NAME) )
    allocate( FileName, source = StorageDir//PIEnergyLabel//PIon.GetLabel() )
    call OpenFile( FileName, uid, 'read', 'formatted' )
    read(uid,*) energy
    close( uid )
    !
  end function ClassParentIonReadEnergy


  function ClassParentIonGetCharge( ParentIon ) result( charge )
    class(ClassParentIon), intent(in) :: ParentIon
    integer                           :: charge
    charge = ParentIon.Charge
  end function ClassParentIonGetCharge


  subroutine ClassParentIonFree( ParentIon )
    class(ClassParentIon) :: ParentIon
    parentIon.Energy = 0.d0
    parentIon.Charge = 0
    parentIon.Multiplicity = 0
    parentIon.Irrep  => NULL()
    parentIon.N      = 0
    parentIon.NePI   = -1
  end subroutine ClassParentIonFree



  subroutine ClassParentIonSaveData( Pion, dir ) 
    class(ClassParentIon), intent(inout) :: Pion
    character(len=*)     , intent(in)    :: dir
    !
    character(len=:), allocatable :: FileName, StorageDir
    integer                       :: uid, iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    allocate( StorageDir, source=AddSlash(dir)//AddSlash(PARENTION_DIR_NAME) )
    call Execute_Command_Line("mkdir -p "//StorageDir)
    allocate( FileName, source = StorageDir//PIEnergyLabel//Pion.GetLabel() )
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="unknown",&
         action  ="write",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)call Assert(iomsg)
    write(uid,*) Pion.Energy
    close(uid)
    !
  end subroutine ClassParentIonSaveData


  subroutine ClassParentIonSaveAllDataToUnit( Pion, uid ) 
    class(ClassParentIon), intent(inout) :: Pion
    integer              , intent(in)    :: uid
    !
    integer                       :: iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    write(uid,*) Pion.Energy
    write(uid,*) Pion.Charge
    write(uid,*) Pion.Multiplicity
    write(uid,*) Pion.Irrep.GetGroupName()
    write(uid,*) Pion.Irrep.GetName()
    write(uid,*) Pion.N
    write(uid,*) Pion.NePI
    !
  end subroutine ClassParentIonSaveAllDataToUnit


  subroutine ClassParentIonLoadAllDataFromUnit( Pion, uid ) 
    class(ClassParentIon), intent(inout) :: Pion
    integer              , intent(in)    :: uid
    !
    integer                       :: iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    character(len=256) :: IrrepName, GroupName
    type(ClassGroup) :: Group
    !
    read(uid,*) Pion.Energy
    read(uid,*) Pion.Charge
    read(uid,*) Pion.Multiplicity
    read(uid,*) GroupName
    read(uid,*) IrrepName
    read(uid,*) Pion.N
    read(uid,*) Pion.NePI
    !
    call Group.Init( trim(GroupName) )
    if( Group.DefinesIrrep( trim(IrrepName) ) )then
       Pion.Irrep => Group.GetIrrep( trim(IrrepName) )
    else
       call Assert("Invalid irrep "//trim(IrrepName)//" for group "//Group.GetName())
    end if
    !
  end subroutine ClassParentIonLoadAllDataFromUnit


  subroutine ClassParentIonLoadData( Pion, dir ) 
    class(ClassParentIon), intent(inout) :: Pion
    character(len=*)     , intent(in)    :: dir
    !
    character(len=:), allocatable :: FileName
    integer                       :: uid, iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    allocate( FileName, source = AddSlash(dir)//AddSlash(PARENTION_DIR_NAME)//PIEnergyLabel//Pion.GetLabel())
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="old",&
         action  ="read",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Cannot Find File "//FileName)
       call Assert(iomsg)
    endif
    !
    read(uid,fmt=*,iostat=iostat) Pion.Energy
    close(uid)
    !
  end subroutine ClassParentIonLoadData


  subroutine ClassParentIonFinal( ParentIon )
    type(ClassParentIon) :: ParentIon
    call ParentIon.free()
  end subroutine ClassParentIonFinal


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !         Procedures   ClassMultipoles
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  
  subroutine ClassMultipoleSetLMax( Multipole, LMax ) 
    class(ClassMultipole), intent(inout) :: Multipole
    integer              , intent(in)    :: LMax
    Multipole.LMax = LMax
  end subroutine ClassMultipoleSetLMax


  subroutine ClassMultipoleSetIrrep( Multipole, IrrepPtr ) 
    class(ClassMultipole)    , intent(inout) :: Multipole
    class(ClassIrrep), target, intent(in)    :: IrrepPtr
    Multipole.Irrep => IrrepPtr
  end subroutine ClassMultipoleSetIrrep


  subroutine ClassMultipoleSetDir( Multipole, DirName ) 
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: DirName
    if(allocated(Multipole.Dir))deallocate(Multipole.Dir)
    allocate(Multipole.Dir,source=DirName)
  end subroutine ClassMultipoleSetDir


  subroutine ClassMultipoleSetFile( Multipole, FileName ) 
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: FileName
    if(allocated(Multipole.File))deallocate(Multipole.File)
    allocate(Multipole.File,source=FileName)
  end subroutine ClassMultipoleSetFile


  subroutine ClassMultipoleSetCartFile( Multipole, CartFileName ) 
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: CartFileName
    if(allocated(Multipole.CartFile))deallocate(Multipole.CartFile)
    allocate(Multipole.CartFile,source=CartFileName)
  end subroutine ClassMultipoleSetCartFile


  subroutine ClassMultipoleSetEmpty( Multipole ) 
    class(ClassMultipole), intent(inout) :: Multipole
    Multipole.NMultipoles = 0
  end subroutine ClassMultipoleSetEmpty


  subroutine ClassMultipoleInit( Multipole, LMax, PionBra, PionKet ) 
    class(ClassMultipole), intent(inout) :: Multipole
    integer              , intent(in)    :: LMax
    class(ClassParentIon), target, intent(in) :: PionBra, PionKet
    !
    call Multipole.SetLMax( LMax )
    Multipole.PionBra => PionBra
    Multipole.PionKet => PionKet
    call Multipole.SetIrrep( PionBra.GetIrrep() * PionKet.GetIrrep() )
    !
    call Multipole.CartesianSet.Init( LMax, Multipole.Irrep )
    call Multipole.CartesianSet.ComputeTransfMat()
    !
    call Multipole.XlmSymSet.Init( LMax, Multipole.Irrep )
    !
    Multipole.NMultipoles = Multipole.CartesianSet.GetN()
    if(allocated(Multipole.CartesianMultipole))deallocate(Multipole.CartesianMultipole)
    allocate(Multipole.CartesianMultipole(Multipole.NMultipoles))
    Multipole.CartesianMultipole = 0.d0
  end subroutine ClassMultipoleInit



  subroutine ClassMultipoleFree( Multipole ) 
    class(ClassMultipole), intent(inout) :: Multipole
    !
    Multipole.LMax = -1
    Multipole.PionBra  => NULL()
    Multipole.PionKet  => NULL()
    Multipole.Irrep    => NULL()
    !
    if ( allocated(Multipole.Dir) ) deallocate( Multipole.Dir )
    if ( allocated(Multipole.File) ) deallocate( Multipole.File )
    if ( allocated(Multipole.CartFile) ) deallocate( Multipole.CartFile )
    !
    Multipole.Nmultipoles = -1
    call Multipole.CartesianSet.Free()
    call Multipole.XlmSymSet.Free( )
    if ( allocated(Multipole.XlmMultipole) ) deallocate( Multipole.XlmMultipole )
  end subroutine ClassMultipoleFree


  
  subroutine ClassMultipoleFinal( Multipole )
    type(ClassMultipole) :: Multipole
    call Multipole.Free()
  end subroutine ClassMultipoleFinal



  !> Get the Xlm multipole between parent ions.
  real(kind(1d0)) function ClassMultipoleGetXlmMultipole( Multipole, l, m ) result(XlmMult)
    class(ClassMultipole), intent(inout) :: Multipole
    integer,               intent(in) :: l
    integer,               intent(in) :: m
    integer :: Mindex
    !
    XlmMult = 0.d0
    if ( .not.allocated(Multipole.XlmMultipole) ) call Assert( 'Impossible to get parent ion Xlm multipole, multipoles have not been allocated' )
    if ( .not.allocated(Multipole.XlmMultipole(l).v) ) then
       call ErrorMessage( 'The parent ion multipole with l and m: '//AlphabeticNumber(l)//' and '//AlphabeticNumber(m)//' is not available, put to zero.' )
       return
    end if
    !
    Mindex = Multipole.GetMindex( l, m )
    if ( Mindex < 0 ) return
    XlmMult = dble(Multipole.XlmMultipole(l).v(Mindex))
    !
  end function ClassMultipoleGetXlmMultipole



  !> Gets for a given angular momentum l, its projection m corresponding index in the vector of Xlm multipoles.
  integer function ClassMultipoleGetMindex( Multipole, l, m ) result(Mindex)
    class(ClassMultipole), intent(inout) :: Multipole
    integer,               intent(in) :: l
    integer,               intent(in) :: m
    !
    integer, allocatable :: VectorOfM(:)
    complex(kind(1d0)), allocatable :: CartExpansionMat(:,:)
    integer :: i
    !
    Mindex = -1
    call Multipole.XlmSymSet.GetCartesianExpan(  Multipole.CartesianSet,  l,  VectorOfM,  CartExpansionMat  )
    !
    do i = 1, size( VectorOfM )
       if ( m == VectorOfM(i) ) then
          Mindex = i
          return
       end if
    end do
    !
!!$    if ( Mindex < 0 ) call Assert( 'Error finding the m index of Xlm parent ion multipole.' )
    !
  end function ClassMultipoleGetMindex




  subroutine ClassMultipoleLoadCartesianData( Multipole, dir ) 
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: dir
    character(len=:), allocatable :: CartFileName
    character(len=IOMSG_LENGTH)::iomsg
    character(len=512) :: line
    integer :: iostat, iChar
    integer :: uid, ix,iy,iz
    logical :: DataReadCorrectly
    real(kind(1d0)) :: dMult
    !
    call SetCartFileName()
    !
    open(newunit = uid, &
         file    = CartFileName, &
         form    ="formatted", &
         status  ="old",&
         action  ="read",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)then
       call ErrorMessage("Cannot Find File "//CartFileName)
       call Assert(iomsg)
    endif
    !
    DataReadCorrectly = .FALSE.
    !
    do
       read(uid,"(a)",iostat=iostat)line
       if(iostat/=0)exit
       iChar=index(line,"#")
       if(iChar>=1)line=line(:iChar-1)
       if(len_trim(line)==0)cycle
       read(line,*,iostat=iostat)ix,iy,iz,dMult
       if(iostat>0)call Assert("Error in multipole value in "//CartFileName)
       if ( iostat < 0 ) exit
       call Multipole.SetCartesianValue( ix, iy, iz, dMult )
       DataReadCorrectly = .TRUE.
    enddo
    !
    close(uid)
    !
    if( .not. DataReadCorrectly ) call Assert("Couldn't read multipoles in "//CartFileName)
    !
  contains
    !
    subroutine SetCartFileName()
      integer :: iChar
      character(len=1024) :: tmpDir
      tmpDir=adjustl(dir)
      call FormatDirectoryName( tmpDir )
      allocate( CartFileName, source = trim(tmpDir)//Multipole.CartFile )
    end subroutine SetCartFileName
    !
  end subroutine ClassMultipoleLoadCartesianData




  subroutine ClassMultipoleConvertCartesianToSpherical( Multipole )
    !
    class(ClassMultipole), intent(inout) :: Multipole
    !
    integer :: l,m
    integer, allocatable :: VectorOfM(:)
    complex(kind(1d0)), allocatable :: CartExpansionMat(:,:)
    !
    if( allocated( Multipole.XlmMultipole ) ) deallocate( Multipole.XlmMultipole )
    allocate( Multipole.XlmMultipole( 0 : Multipole.Lmax ) )
    !
    do l = 0, Multipole.Lmax
       call Multipole.XlmSymSet.GetCartesianExpan(  Multipole.CartesianSet,  l,  VectorOfM,  CartExpansionMat  )
       if( .not.allocated(VectorOfM) ) cycle
       allocate(  Multipole.XlmMultipole( l ).v( size( VectorOfM ) ) )
       !
       do m = 1, size( VectorOfM )
          Multipole.XlmMultipole( l ).v( m ) = sum( CartExpansionMat( m,:) * Multipole.CartesianMultipole(:) )
       end do
       !
    enddo
    !
  end subroutine ClassMultipoleConvertCartesianToSpherical



  function ClassMultipoleFileName( Multipole ) result( FileName )
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=:), allocatable :: FileName
    allocate(FileName,source=MultLabelStore//Multipole.PionBra.GetLabel()//"-"//Multipole.PionKet.GetLabel())
  end function ClassMultipoleFileName



  subroutine ClassMultipoleSave( Multipole, Storage )
    !
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: Storage
    !
    character(len=:), allocatable :: FileName, StorageDir
    integer                       :: uid, iostat, i, l
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    allocate( StorageDir, source=AddSlash(Storage)//AddSlash(PARENTION_DIR_NAME) )
    call Execute_Command_Line("mkdir -p "//StorageDir)
    allocate( FileName, source = StorageDir//Multipole.FileName() )
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="unknown",&
         action  ="write",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)call Assert(iomsg)
    !
    write(uid,*) Multipole.lmax
    write(uid,*) Multipole.Nmultipoles
    write(uid,"(*"//MULTIPOLE_IO_FORMAT//")") (Multipole.CartesianMultipole(i),i=1,Multipole.Nmultipoles)
    write(uid,"(L)") allocated(Multipole.XlmMultipole)
    if(allocated(Multipole.XlmMultipole))then
       do l=0,Multipole.Lmax
          write(uid,"(L)") allocated(Multipole.XlmMultipole(l).v)
          if( allocated(Multipole.XlmMultipole(l).v) )then
             write(uid,*) size(Multipole.XlmMultipole(l).v)
             write(uid,"(*"//MULTIPOLE_IO_FORMAT//")") (Multipole.XlmMultipole(l).v(i),i=1,size(Multipole.XlmMultipole(l).v))
          endif
       enddo
    endif
    !
    close(uid)
    !
  end subroutine ClassMultipoleSave



  subroutine ClassMultipoleRead( Multipole, Storage )
    !
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: Storage
    !
    character(len=:), allocatable :: FileName, StorageDir
    integer                       :: uid, iostat, i, l
    character(len=IOMSG_LENGTH)   :: iomsg
    logical                       :: AvailXlmMult, AvailVec
    integer                       :: DimVec, Lmax
    logical                       :: exist
    !
    allocate( StorageDir, source=AddSlash(Storage)//AddSlash(PARENTION_DIR_NAME) )
    allocate( FileName, source = StorageDir//Multipole.FileName() )
    INQUIRE( file = FileName, exist = exist )
    if ( .not.exist ) then
       call ErrorMessage( 'The multipole: '//FileName//' has not been saved previously, impossible to read.' )
       return
    end if
    !
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="old",&
         action  ="read",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)call Assert(iomsg)
    !
    read(uid,*) Lmax
    if(Lmax /= Multipole.Lmax) call Assert( "incompatible lmax in "//FileName )
    Multipole.Lmax = Lmax
    read(uid,*) Multipole.Nmultipoles
    if ( allocated(Multipole.CartesianMultipole) ) deallocate( Multipole.CartesianMultipole )
    allocate( Multipole.CartesianMultipole(Multipole.Nmultipoles) )
    read(uid,"(*"//MULTIPOLE_IO_FORMAT//")") (Multipole.CartesianMultipole(i),i=1,Multipole.Nmultipoles)
    read(uid,"(L)") AvailXlmMult
    if ( AvailXlmMult ) then
       if ( allocated(Multipole.XlmMultipole) ) deallocate( Multipole.XlmMultipole )
       allocate( Multipole.XlmMultipole(0:Multipole.Lmax) )
       do l=0,Multipole.Lmax
          read(uid,"(L)") AvailVec
          if( AvailVec )then
             read(uid,*) DimVec
             if ( allocated(Multipole.XlmMultipole(l).v) ) deallocate( Multipole.XlmMultipole(l).v )
             allocate( Multipole.XlmMultipole(l).v(DimVec) )
             read(uid,"(*"//MULTIPOLE_IO_FORMAT//")") (Multipole.XlmMultipole(l).v(i),i=1,DimVec)
          endif
       enddo
    endif
    !
    close(uid)
    !
  end subroutine ClassMultipoleRead



  !> Returns the number of non-vanishing symmetric multipoles
  function ClassMultipoleGetN( Multipole ) result( N )
    class(ClassMultipole), intent(in) :: Multipole
    integer :: N
    N = Multipole.NMultipoles
  end function ClassMultipoleGetN



  !> Set the value of the cartesian multipole introducing the x, y and z exponents.
  subroutine ClassMultipoleSetCartesianValue( Multipole, ix, iy, iz, dValue ) 
    class(ClassMultipole), intent(inout) :: Multipole
    integer              , intent(in)    :: ix, iy, iz
    real(kind(1d0))      , intent(in)    :: dValue
    integer :: iMult
    iMult = Multipole.CartesianSet.GetIndex(ix,iy,iz)
    if(iMult==0)call ErrorMessage("unmatching triplet")
    Multipole.CartesianMultipole(iMult)=dValue
  end subroutine ClassMultipoleSetCartesianValue


  !> Get the value of the cartesian multipole introducing the x, y and z exponents.
  function ClassMultipoleGetCartesianValue( Multipole, ix, iy, iz ) result(dValue) 
    class(ClassMultipole), intent(inout)  :: Multipole
    integer              , intent(in)  :: ix, iy, iz
    real(kind(1d0))                    :: dValue
    integer :: iMult
    dValue = 0.d0
    iMult = Multipole.CartesianSet.GetIndex(ix,iy,iz)
    if(iMult==0)call ErrorMessage("unmatching triplet")
    if ( iMult > 0 ) dValue = Multipole.CartesianMultipole(iMult)
  end function ClassMultipoleGetCartesianValue



  !> Set the value of the cartesian multipole introducing the label, i.e: xx, xyzz, etc.
  subroutine ClassMultipoleSetCartEl( Multipole, Label, dValue ) 
    class(ClassMultipole), intent(inout) :: Multipole
    character(len=*)     , intent(in)    :: Label
    real(kind(1d0))      , intent(in)    :: dValue
    integer :: ix, iy, iz
    !
    !..Determine the number of x, y and z
    call FindRepeatingNumberOfPatternInString( Label, 'x', ix )
    call FindRepeatingNumberOfPatternInString( Label, 'y', iy )
    call FindRepeatingNumberOfPatternInString( Label, 'z', iz )
    call Multipole.SetCartesianValue( ix, iy, iz, dValue )
  end subroutine ClassMultipoleSetCartEl


  !> Get the value of the cartesian multipole introducing the label, i.e: xx, xyzz, etc.
  function ClassMultipoleFetchCartesianMultipole( Multipole, Label ) result(dValue) 
    class(ClassMultipole), intent(inout)  :: Multipole
    character(len=*)     , intent(in)  :: Label
    real(kind(1d0))                    :: dValue
    integer :: ix, iy, iz
    !
    !..Determine the number of x, y and z
    call FindRepeatingNumberOfPatternInString( Label, 'x', ix )
    call FindRepeatingNumberOfPatternInString( Label, 'y', iy )
    call FindRepeatingNumberOfPatternInString( Label, 'z', iz )
    dValue = Multipole.GetCartesianValue( ix, iy, iz )
  end function ClassMultipoleFetchCartesianMultipole


  !> Returns " x  xxx  xyy  xzz " etc.
  subroutine ClassMultipoleFetchCartesianStrn( Multipole, Strn ) 
    class(ClassMultipole),         intent(inout) :: Multipole
    character(len=:), allocatable, intent(out)   :: Strn
    call Multipole.CartesianSet.FetchMonomialStrn( Strn )
  end subroutine ClassMultipoleFetchCartesianStrn


  !> Returns the number of non-vanishing symmetric multipoles
  subroutine ClassMultipoleShow( Multipole, unit ) 
    class(ClassMultipole), intent(in) :: Multipole
    integer, optional    , intent(in) :: unit
    integer :: outunit, iMult, l, m
    integer, allocatable :: mvec(:)
    outunit=OUTPUT_UNIT
    if(present(unit))outunit=unit
    write(outunit,"(a,i2)") "  Lmax =",Multipole.LMax
    write(outunit,"(a,i3)") "  NMult=",Multipole.NMultipoles
    if(Multipole.NMultipoles<=0)return
    if(allocated(Multipole.File    )) write(outunit,"(a)")    "  File ="//Multipole.File
    if(allocated(Multipole.CartFile)) write(outunit,"(a)")    "  CartFile ="//Multipole.CartFile
    call Multipole.Irrep.Show(unit)
    call Multipole.CartesianSet.show(unit)
    call Multipole.XlmSymSet.Show()
    !
    if(allocated(Multipole.CartesianMultipole))then
       do iMult=1,Multipole.NMultipoles
          write(outunit,"(a,i3,a,d14.6)")&
               " Cart Mult(",iMult,") =",Multipole.CartesianMultipole(iMult)
       enddo
    endif
    !
    if(allocated( Multipole.XlmMultipole ))then
       !
       do l = 0, Multipole.Lmax
          if(.not.allocated(Multipole.XlmMultipole(l).v))cycle
          call Multipole.XlmSymSet.GetMlist(l,mvec)
          !
          do m = 1, size(Multipole.XlmMultipole(l).v)
             write(outunit,"(a,i2,a,i2,a,*(x,d14.6))")&
                  " Sph  Mult(",l,",",mvec(m),") = ",Multipole.XlmMultipole(l).v(m)
          enddo
          !
       enddo
    endif
    !
  end subroutine ClassMultipoleShow



  subroutine ParseParentIonsString( NPions, vPions, PIonsStrn )
    integer                      , intent(out) :: NPions
    character(len=*), allocatable, intent(out) :: VPions(:)
    character(len=*)             , intent(in)  :: PionsStrn
    character(len=:), allocatable :: Strn
    character(len=8) :: PionLabel
    integer :: i
    !
    allocate(strn,source=trim(adjustl(PionsStrn))//" ")
    !
    !.. Determines the number of parent ions
    i=1
    NPions=0
    do
       if(len_trim(strn)<=0)exit
       read(strn,*) PionLabel
       NPions=NPions+1
       i=index(strn," ")
       strn=adjustl(strn(i:))
    enddo
    !
    allocate(vPions(NPions))
    !.. Read the parent ions
    strn=trim(adjustl(PionsStrn))//" "
    i=1
    NPions=0
    do
       if(len_trim(strn)<=0)exit
       read(strn,*) PionLabel
       NPions=NPions+1
       VPions(NPions)=trim(adjustl(PionLabel))
       i=index(strn," ")
       strn=adjustl(strn(i:))
    enddo
    !
  end subroutine ParseParentIonsString


!--------------------------------------------
!--- Procedures ClassKineticEnergy
!--------------------------------------------

  subroutine ClassKineticEnergySetDir( Self, DirName ) 
    class(ClassKineticEnergy), intent(inout) :: Self
    character(len=*)         , intent(in)    :: DirName
    if(allocated(Self.Dir))deallocate(Self.Dir)
    allocate(Self.Dir,source=DirName)
  end subroutine ClassKineticEnergySetDir


  subroutine ClassKineticEnergySetFile( Self, FileName ) 
    class(ClassKineticEnergy), intent(inout) :: Self
    character(len=*)         , intent(in)    :: FileName
    if(allocated(Self.File))deallocate(Self.File)
    allocate(Self.File,source=FileName)
  end subroutine ClassKineticEnergySetFile


  subroutine ClassKineticEnergyInit( Self, PionBra, PionKet ) 
    class(ClassKineticEnergy)    , intent(inout) :: Self
    class(ClassParentIon), target, intent(in)    :: PionBra, PionKet
    Self.PionBra => PionBra
    Self.PionKet => PionKet
  end subroutine ClassKineticEnergyInit
  
  
  real(kind(1d0)) function ClassKineticEnergyGetKinEnergy( Self ) result(KE)
    class(ClassKineticEnergy), intent(in) :: Self
    KE = Self.KinEnergy
  end function ClassKineticEnergyGetKinEnergy
  
  subroutine ClassKineticEnergyShow( Self, unit )
    class(ClassKineticEnergy), intent(in) :: Self
    integer, optional    , intent(in) :: unit
    !
    integer :: outunit
    !
    outunit=OUTPUT_UNIT
    if(present(unit))outunit=unit
    !
    if(allocated(Self.File    )) write(outunit,"(a)")    "  File ="//Self.File
    write(outunit,"(a,f24.16)") "  Kin Energy =",Self.KinEnergy
    !
  end subroutine ClassKineticEnergyShow

  

  subroutine ClassKineticEnergySave( Self, Storage )
    !
    class(ClassKineticEnergy), intent(inout) :: Self
    character(len=*)         , intent(in)    :: Storage
    !
    character(len=:), allocatable :: FileName, StorageDir
    integer                       :: uid, iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    allocate( StorageDir, source=AddSlash(Storage)//AddSlash(PARENTION_DIR_NAME) )
    call Execute_Command_Line("mkdir -p "//StorageDir)
    allocate( FileName, source = StorageDir//Self.FileName() )
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="unknown",&
         action  ="write",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)call Assert(iomsg)
    !
    write(uid,*) Self.KinEnergy
    !
    close(uid)
    !
  end subroutine ClassKineticEnergySave



  subroutine ClassKineticEnergyRead( Self, Storage )
    !
    class(ClassKineticEnergy), intent(inout) :: Self
    character(len=*)         , intent(in)    :: Storage
    !
    character(len=:), allocatable :: FileName, StorageDir
    integer                       :: uid, iostat
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    allocate( StorageDir, source=AddSlash(Storage)//AddSlash(PARENTION_DIR_NAME) )
    allocate( FileName, source = StorageDir//Self.FileName() )
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="old",&
         action  ="read",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0)call Assert(iomsg)
    !
    read(uid,*) Self.KinEnergy
    !
    close(uid)
    !
  end subroutine ClassKineticEnergyRead


  function ClassKineticEnergyFileName( Self ) result( FileName )
    class(ClassKineticEnergy), intent(inout) :: Self
    character(len=:), allocatable :: FileName
    allocate(FileName,source=KinEnergyLabelStore//Self.PionBra.GetLabel()//"-"//Self.PionKet.GetLabel())
  end function ClassKineticEnergyFileName


  subroutine ClassKineticEnergyFree( Self ) 
    class(ClassKineticEnergy), intent(inout) :: Self
    !
    Self.PionBra  => NULL()
    Self.PionKet  => NULL()
    if ( allocated(Self.Dir) ) deallocate( Self.Dir )
    if ( allocated(Self.File) ) deallocate( Self.File )
    Self.KinEnergy = 0.d0
  end subroutine ClassKineticEnergyFree


  subroutine ClassKineticEnergyFinal( Self )
    type(ClassKineticEnergy) :: Self
    call Self.Free()
  end subroutine ClassKineticEnergyFinal


  subroutine ClassKineticEnergyFetchKinEnergy( Self, FileName, IBra, IKet )
    !
    class(ClassKineticEnergy), intent(inout) :: Self
    character(len=*)         , intent(in)    :: Filename
    integer                  , intent(in)    :: IBra, IKet
    !
    real(kind(1d0)) :: KE
    integer :: uid, iostat, i, j
    character(len=IOMSG_LENGTH)   :: iomsg
    !
    open(newunit = uid, &
         file    = FileName, &
         form    ="formatted", &
         status  ="old",&
         action  ="read",&
         iostat  = iostat,&
         iomsg   = iomsg )
    if(iostat/=0) then
       call ErrorMessage(iomsg)
       return
    end if
    !
    ! Skip first lines
    do i = 1, IBra
       read(uid,*)
    end do
    read(uid,*) ( KE, j = 1,IKet )
    !
    close( uid )
    !
    Self.KinEnergy = KE
    !
  end subroutine ClassKineticEnergyFetchKinEnergy


end module ModuleParentIons
