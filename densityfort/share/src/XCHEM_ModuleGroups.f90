!! CONFIDENTIAL
!> Defines the composition rules and irreducible representation
!! of the abelian group D2h and its subgroups 
!!
!! C1, Cs, C2, Ci, C2v, C2h, and D2.
!! 
!! Since these groups are commutative, the number of classes and, 
!! thus, of irreducible representations, coincides with the order 
!! of the group.
!!
module ModuleGroups
  ! {{{ Header 

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString

  implicit none

  private

  logical :: MODULE_IS_INITIALIZED = .FALSE.

  ! }}}


  !> Symmetry operations. E.g.: rotations, reflection planes, inversion, etc.
  type, public :: SymmetryOperation
     ! {{{ private attributes
     private
     character(len=:), allocatable :: Name
     real(kind(1d0)) :: R(3,3)
     ! }}}
   contains
     generic  , public  :: init           => SymmetryOperationInit
     generic  , public  :: free           => SymmetryOperationFree
     generic  , public  :: show           => SymmetryOperationShow
     generic  , public  :: write          => SymmetryOperationWriteUnit
     generic  , public  :: read           => SymmetryOperationReadUnit
     generic  , public  :: assignment(=)  => SymmetryOperationAssign
     ! {{{ private procedures
     procedure, private :: SymmetryOperationAssign
     procedure, private :: SymmetryOperationInit
     procedure, private :: SymmetryOperationFree
     procedure, private :: SymmetryOperationShow
     procedure, private :: SymmetryOperationWriteUnit
     procedure, private :: SymmetryOperationReadUnit
     final :: SymmetryOperationFinal
     ! }}}
  end type SymmetryOperation
  type(SymmetryOperation), public, protected :: R_E, R_I, R_C2x, R_C2y, R_C2z, R_Qxy, R_Qxz, R_Qyz


  !> Irreducible representations. A.k.a., symmetries.
  type, public :: ClassIrrep
     ! {{{ private attributes

     private
     type(ClassGroup), pointer     :: Group => NULL()
     character(len=:), allocatable :: Name
     integer                       :: NClasses
     real(kind(1d0)) , allocatable :: Characters(:)

     ! }}} 
   contains
     generic  , public  :: init           => ClassIrrepInit
     generic  , public  :: free           => ClassIrrepFree
     generic  , public  :: show           => ClassIrrepShow
     generic  , public  :: write          => ClassIrrepWriteUnit
     generic  , public  :: read           => ClassIrrepReadUnit
     generic  , public  :: assignment(=)  => ClassIrrepAssign
     generic  , public  :: operator(*)    => ClassIrrepOTimes
     generic  , public  :: setGroup       => ClassIrrepSetGroup
     generic  , public  :: getName        => ClassIrrepGetName
     generic  , public  :: getGroupName   => ClassIrrepGetGroupName
     generic  , public  :: getGroup       => ClassIrrepGetGroup
     generic  , public  :: getIdStrn      => ClassIrrepGetIdStrn
     generic  , public  :: nameIs         => ClassIrrepNameIs
     ! {{{ private procedures

     procedure, private :: ClassIrrepAssign
     procedure, private :: ClassIrrepOTimes
     procedure, private :: ClassIrrepInit
     procedure, private :: ClassIrrepFree
     procedure, private :: ClassIrrepShow
     procedure, private :: ClassIrrepWriteUnit
     procedure, private :: ClassIrrepReadUnit
     procedure, private :: ClassIrrepNameIs
     procedure, private :: ClassIrrepSetGroup
     procedure, private :: ClassIrrepGetName
     procedure, private :: ClassIrrepGetGroupName
     procedure, private :: ClassIrrepGetGroup
     procedure, private :: ClassIrrepGetIdStrn
     final :: ClassIrrepFinal

     ! }}}
  end type ClassIrrep


  !> (Abelian) groups
  type, public :: ClassGroup
     ! {{{ private attributes
     private
     logical                              :: Initialized = .FALSE.
     character(len=:)       , allocatable :: Name
     integer                              :: Order
     type(SymmetryOperation), allocatable :: Elements(:)
     type(ClassIrrep)       , allocatable :: Irreps(:)
     ! }}} 
   contains
     generic, public :: init          => ClassGroupInitDefault, ClassGroupInitNew
     generic, public :: IsInitialized => ClassGroupIsInitialized
     generic, public :: free          => ClassGroupFree
     generic, public :: show          => ClassGroupShow
     generic, public :: getName       => ClassGroupGetName
     generic, public :: definesIrrep  => ClassGroupDefinesIrrepName
     generic, public :: getnIrreps    => ClassGroupGetnIrreps
     generic, public :: getIrrep      => ClassGroupGetIrrepPtrFromName, ClassGroupGetIrrepPtrFromIndex 
     generic, public :: getTotSymIrrep=> ClassGroupGetTotSymIrrepPtr
     generic, public :: getIrrepIndex => ClassGroupGetIrrepIndex
     generic, public :: getIrrepList  => ClassGroupGetIrrepList
     generic, public :: GetMonomialIrrep=> ClassGroupGetMonomialIrrep
     generic, public :: GetTotallySymmetricIrrep=> ClassGroupGetTotallySymmetricIrrep
     generic, public :: GetxIrrep     => ClassGroupGetxIrrep
     generic, public :: GetyIrrep     => ClassGroupGetyIrrep
     generic, public :: GetzIrrep     => ClassGroupGetzIrrep
     generic, public :: getIrrName    => ClassGroupGetIrrepName
     ! {{{ private procedures

     procedure, private :: ClassGroupInitDefault
     procedure, private :: ClassGroupInitNew
     procedure, private :: ClassGroupIsInitialized
     procedure, private :: ClassGroupFree
     procedure, private :: ClassGroupShow
     procedure, private :: ClassGroupGetName
     procedure, private :: ClassGroupDefinesIrrepName
     procedure, private :: ClassGroupGetIrrepPtrFromName
     procedure, private :: ClassGroupGetIrrepPtrFromIndex
     procedure, private :: ClassGroupGetTotSymIrrepPtr
     procedure, private :: ClassGroupGetnIrreps
     procedure, private :: ClassGroupGetIrrepList
     procedure, private :: ClassGroupGetMonomialIrrep
     procedure, private :: ClassGroupGetTotallySymmetricIrrep
     procedure, private :: ClassGroupGetxIrrep
     procedure, private :: ClassGroupGetyIrrep
     procedure, private :: ClassGroupGetzIrrep
     procedure, private :: ClassGroupGetIrrepIndex
     procedure, private :: ClassGroupGetIrrepName
     final :: ClassGroupFinal

     ! }}}
  end type ClassGroup
  type(ClassGroup), public, target :: GlobalGroup


  public :: ValidIrreps


contains 



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Module procedures
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  

  subroutine InitModuleGroups
    if( MODULE_IS_INITIALIZED )return
    call InitElements
    MODULE_IS_INITIALIZED = .TRUE.
  end subroutine InitModuleGroups


  subroutine InitElements
    R_E   = SymmetryOperation( "E"  , RESHAPE( [ 1, 0, 0,  0, 1, 0,  0, 0, 1 ],  [3,3] ) )
    R_I   = SymmetryOperation( "I"  , RESHAPE( [-1, 0, 0,  0,-1, 0,  0, 0,-1 ],  [3,3] ) )
    R_C2x = SymmetryOperation( "C2x", RESHAPE( [ 1, 0, 0,  0,-1, 0,  0, 0,-1 ],  [3,3] ) )
    R_C2y = SymmetryOperation( "C2y", RESHAPE( [-1, 0, 0,  0, 1, 0,  0, 0,-1 ],  [3,3] ) )
    R_C2z = SymmetryOperation( "C2z", RESHAPE( [-1, 0, 0,  0,-1, 0,  0, 0, 1 ],  [3,3] ) )
    R_Qxy = SymmetryOperation( "Qxy", RESHAPE( [ 1, 0, 0,  0, 1, 0,  0, 0,-1 ],  [3,3] ) )
    R_Qxz = SymmetryOperation( "Qxz", RESHAPE( [ 1, 0, 0,  0,-1, 0,  0, 0, 1 ],  [3,3] ) )
    R_Qyz = SymmetryOperation( "Qyz", RESHAPE( [-1, 0, 0,  0, 1, 0,  0, 0, 1 ],  [3,3] ) )
  end subroutine InitElements


  ! {{{ C1 description
  !> \f$ C_1 \f$ group.
  ! C1     E
  ! 
  ! E      E
  !
  !
  ! C1     E 
  !
  ! A      1
  !
  !
  ! C1     A
  ! 
  ! A      A
  ! }}}
  subroutine InitGroupC1(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0 !.. u stands for unity
    call Group%free
    call Group%init( "C1", [ R_E ], [ &
         ClassIrrepFarm( "A", [ u ] ) ] )
  end subroutine InitGroupC1


  ! {{{ C2 description
  ! C2     E  C2z
  ! 
  ! E      E  C2z
  ! C2z   C2z  E
  !
  !
  ! C2     E  C2z
  !
  ! A      1   1
  ! B      1  -1
  !
  !
  ! C2     A   B
  !
  ! A      A   B
  ! B      B   A
  ! }}}
  subroutine InitGroupC2(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "C2", [ R_E, R_C2z ], [ &
         ClassIrrepFarm( "A", [ u, u ] ), &
         ClassIrrepFarm( "B", [ u,-u ] ) ] )
  end subroutine InitGroupC2


  ! {{{ Cs description 

  ! Cs     E  Qxy
  ! 
  ! E      E  Qxy
  ! Qxy   Qxy  E
  !
  !
  ! Cs     E  Qxy
  ! 
  ! A'     1   1
  ! A"     1  -1
  !
  !
  ! Cs     A'  A"
  ! 
  ! A'     A'  A"
  ! A"     A"  A'

  ! }}}
  subroutine InitGroupCs(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "Cs", [ R_E, R_Qxy ], [ &
         ClassIrrepFarm( "Ap" , [ u, u ] ), &
         ClassIrrepFarm( "App", [ u,-u ] ) ] )
  end subroutine InitGroupCs


  ! {{{ Ci description
  ! Ci     E   i
  !
  ! Ag     1   1
  ! Au     1  -1
  !
  !
  ! Ci     E   i
  !
  ! E      E   i
  ! i      i   E
  !
  !
  ! Ci     Ag  Au
  !
  ! Ag     Ag  Au
  ! Au     Au  Ag
  ! }}}
  subroutine InitGroupCi(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "Ci", [ R_E, R_I ], [ &
         ClassIrrepFarm( "Ag", [ u, u ] ), &
         ClassIrrepFarm( "Au", [ u,-u ] ) ] )
  end subroutine InitGroupCi


  ! {{{ C2v description 
  ! C2v    E  C2z Qxz Qyz
  !
  ! E      E  C2z Qxz Qyz
  ! C2z   C2z  E  Qxz Qyz 
  ! Qxz   Qxz Qxz  E  C2z
  ! Qyz   Qyz Qyz C2z  E
  !
  !
  ! C2v    E  C2z Qxz Qyz
  !
  ! A1     1   1   1   1
  ! A2     1   1  -1  -1
  ! B1     1  -1   1  -1
  ! B2     1  -1  -1   1
  !
  !
  ! C2v   A1  A2  B1  B2
  !
  ! A1    A1  A2  B1  B2
  ! A2    A2  A1  B2  B1
  ! B1    B1  B2  A1  A2
  ! B2    B2  B1  A2  A1
  ! }}}
  subroutine InitGroupC2v(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "C2v", [ R_E, R_C2z, R_Qxz, R_Qyz ], [ &
         ClassIrrepFarm( "A1" , [ u, u, u, u ] ), &
         ClassIrrepFarm( "A2" , [ u, u,-u,-u ] ), &
         ClassIrrepFarm( "B1" , [ u,-u, u,-u ] ), &
         ClassIrrepFarm( "B2" , [ u,-u,-u, u ] ) ] )
  end subroutine InitGroupC2v


  ! {{{ C2h description 
  ! C2h    E  C2z  I  Qxy
  !
  ! E      E  C2z  I  Qxy
  ! C2z   C2z  E  Qxy  I
  ! I      I  Qxy  E  C2z
  ! Qxy   Qxy  I  C2z  E
  !
  !
  ! C2h    E  C2z  I  Qxy
  !
  ! Ag     1   1   1   1
  ! Bg     1  -1   1  -1
  ! Au     1   1  -1  -1
  ! Bu     1  -1  -1   1
  !
  !
  ! C2h   Ag  Bg  Au  Bu
  !
  ! Ag    Ag  Bg  Au  Bu
  ! Bg    Bg  Ag  Bu  Au
  ! Au    Au  Bu  Ag  Bg
  ! Bu    Bu  Au  Bg  Ag
  ! }}}
  subroutine InitGroupC2h(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "C2h", [ R_E, R_C2z, R_I, R_Qxy ], [ &
         ClassIrrepFarm( "Ag" , [ u, u, u, u ] ), &
         ClassIrrepFarm( "Bg" , [ u,-u, u,-u ] ), &
         ClassIrrepFarm( "Au" , [ u, u,-u,-u ] ), &
         ClassIrrepFarm( "Bu" , [ u,-u,-u, u ] ) ] )
  end subroutine InitGroupC2h


  ! {{{ D2 description 
  ! D2     E  C2z C2y C2x  
  !
  ! E      E  C2z C2y C2x
  ! C2z   C2z  E  C2x C2y
  ! C2y   C2y C2x  E  C2z
  ! C2x   C2x C2y C2z  E
  !
  !
  ! D2     E  C2z C2y C2x  
  !
  ! A      1   1   1   1   
  ! B1     1   1  -1  -1   
  ! B2     1  -1   1  -1   
  ! B3     1  -1  -1   1   
  !
  !
  ! D2    A   B1  B2  B3
  !
  ! A     A   B1  B2  B3   
  ! B1    B1  A   B3  B2  
  ! B2    B2  B3  A   B1 
  ! B3    B3  B2  B1  A  
  ! }}}
  subroutine InitGroupD2(Group)
    type( ClassGroup ), intent(out) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "D2", [ R_E, R_C2z, R_C2y, R_C2x ], [ &
         ClassIrrepFarm( "A" , [ u, u, u, u ] ), &
         ClassIrrepFarm( "B1", [ u, u,-u,-u ] ), &
         ClassIrrepFarm( "B2", [ u,-u, u,-u ] ), &
         ClassIrrepFarm( "B3", [ u,-u,-u, u ] )  ] )
  end subroutine InitGroupD2
 
  ! {{{ D2h description
  ! D2h    E  C2z C2y C2x    I  Qxy Qxz Qyz
  !
  ! E      E  C2z C2y C2x    I  Qxy Qxz Qyz
  ! C2z   C2z  E  C2x C2y   Qxy  I  Qyz Qxz
  ! C2y   C2y C2x  E  C2z   Qxz Qyz  I  Qxy
  ! C2x   C2x C2y C2z  E    Qyz Qxz Qxy  I
  ! 
  ! I      I  Qxy Qxz Qyz    E  C2z C2y C2x
  ! Qxy   Qxy  I  Qyz Qxz   C2z  E  C2x C2y
  ! Qxz   Qxz Qyz  I  Qxy   C2y C2x  E  C2z
  ! Qyz   Qyz Qxz Qxy  I    C2x C2y C2z  E
  !
  !
  ! D2h    E  C2z C2y C2x    I  Qxy Qxz Qyz
  !
  ! Ag     1   1   1   1     1   1   1   1 
  ! B1g    1   1  -1  -1     1   1  -1  -1 
  ! B2g    1  -1   1  -1     1  -1   1  -1 
  ! B3g    1  -1  -1   1     1  -1  -1   1 
  !
  ! Au     1   1   1   1    -1  -1  -1  -1 
  ! B1u    1   1  -1  -1    -1  -1   1   1 
  ! B2u    1  -1   1  -1    -1   1  -1   1 
  ! B3u    1  -1  -1   1    -1   1   1  -1 
  !
  !
  ! D2h   Ag  B1g B2g B3g   Au  B1u B2u B3u
  !
  ! Ag    Ag  B1g B2g B3g   Au  B1u B2u B3u
  ! B1g   B1g Ag  B3g B2g   B1u Au  B3u B2u
  ! B2g   B2g B3g Ag  B1g   B2u B3u Au  B1u
  ! B3g   B3g B2g B1g Ag    B3u B2u B1u Au
  !
  ! Au    Au  B1u B2u B3u   Ag  B1g B2g B3g
  ! B1u   B1u Au  B3u B2u   B1g Ag  B3g B2g
  ! B2u   B2u B3u Au  B1u   B2g B3g Ag  B1g
  ! B3u   B3u B2u B1u Au    B3g B2g B1g Ag
  !
  !
  !  COMPARISON WITH SYMMETRIES IN DALTON
  !  ====================================
  !  
  !  Ref: Dalton 2018 manual, Chapter 24
  !  In Dalton, the following convention for the order of the symmetry operations
  !  in D2h is used:
  !
  !   1    2    3    4    5    6    7    8  
  !   E   G1   G2   G1G2 G3  G1G3 G2G3 G1G2G3
  !
  !  So, in the case of the standard choice for the generators X Y Z, the operations are
  !
  !   1    2    3    4    5    6    7    8  
  !   E   Qyz  Qxz  C2z  Qxy  C2y  C2x   I
  !
  !  Similarly, here's the irreps' parity properties alongside their standard name when
  !  the generators are X Y Z:
  !
  !   _____________  G1 G2 G3    
  !   1  |   E    |  +  +  +    Ag 
  !   2  |   G1   |  -  +  +    B3u
  !   3  |   G2   |  +  -  +    B2u
  !   4  |  G1G2  |  -  -  +    B1g
  !   5  |   G3   |  +  +  -    B1u
  !   4  |  G1G3  |  -  +  -    B2g
  !   4  |  G2G3  |  +  -  -    B3g
  !   8  | G1G2G3 |  -  -  -    Au
  !
  !  which is definitely weird ...
  !   
  !  The XCHEM ordering of the operations would have been obtained with the generators
  !  
  !  XY XZ XYZ
  !  
  !  in which case the correspondence between the irreps 1..8 in Dalton and the irreps in XCHEM would be
  !  Dalton:  1   2   3   4   5   6   7   8    
  !  XCHEM :  Ag B2g B1g B3g  Au B2u B1u B3u
  !
  ! }}}
  subroutine InitGroupD2h( Group )
    type( ClassGroup ), intent(inout) :: Group
    real(kind(1d0))   , parameter   :: u = 1.d0
    call Group%free
    call Group%init( "D2h", [ R_E, R_C2z, R_C2y, R_C2x, R_I, R_Qxy, R_Qxz, R_Qyz ], [ &
         ClassIrrepFarm( "Ag" , [ u, u, u, u,   u, u, u, u ] ), &
         ClassIrrepFarm( "B1g", [ u, u,-u,-u,   u, u,-u,-u ] ), &
         ClassIrrepFarm( "B2g", [ u,-u, u,-u,   u,-u, u,-u ] ), &
         ClassIrrepFarm( "B3g", [ u,-u,-u, u,   u,-u,-u, u ] ), &
         ClassIrrepFarm( "Au" , [ u, u, u, u,  -u,-u,-u,-u ] ), &
         ClassIrrepFarm( "B1u", [ u, u,-u,-u,  -u,-u, u, u ] ), &
         ClassIrrepFarm( "B2u", [ u,-u, u,-u,  -u, u,-u, u ] ), &
         ClassIrrepFarm( "B3u", [ u,-u,-u, u,  -u, u, u,-u ] ) ] )
  end subroutine InitGroupD2h


  function ClassIrrepFarm( Name, Characters, Group ) result( irrep )
    character(len=*), intent(in) :: Name
    real(kind(1d0)) , intent(in) :: Characters(:)
    type(ClassGroup), optional, target, intent(in) :: Group
    type(ClassIrrep), pointer :: irrep
    allocate(irrep)
    call irrep%init(Name,Characters)
    if(present(Group))then
       irrep%Group => Group
    endif
  end function ClassIrrepFarm


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Class SymmetryOperation Methods
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  subroutine SymmetryOperationAssign( OpOut, OpIn )
    class( SymmetryOperation ), intent(inout) :: OpOut
    class( SymmetryOperation ), intent(in)    :: OpIn
    call OpOut%init(OpIn%Name,OpIn%R)
  end subroutine SymmetryOperationAssign
  
  subroutine SymmetryOperationInit( Op, Name, R )
    class( SymmetryOperation ), intent(inout) :: Op
    character(len=*)          , intent(in)    :: Name
    real(kind(1d0))           , intent(in)    :: R(3,3)
    call Op%free
    allocate(Op%Name,source=Name)
    Op%R=R
  end subroutine SymmetryOperationInit

  subroutine SymmetryOperationFree( Op )
    class( SymmetryOperation ), intent(inout) :: Op
    if(allocated(Op%Name)) deallocate(Op%Name)
    Op%R=0.d0
  end subroutine SymmetryOperationFree

  subroutine SymmetryOperationShow( Op, unit )
    class( SymmetryOperation ), intent(in) :: Op
    integer, optional         , intent(in) :: unit
    integer :: outunit, iRow, iCol
    outunit=OUTPUT_UNIT
    if(present(unit))outunit=unit
    write(outunit,"(a)") "        "//Op%Name
    do iRow=1,3
       write(outunit,"(*(1x,f4.0))")(Op%R(iRow,iCol),iCol=1,3)
    enddo
  end subroutine SymmetryOperationShow

  subroutine SymmetryOperationWriteUnit( Op, unit )
    class( SymmetryOperation ), intent(in) :: Op
    integer                   , intent(in) :: unit
    character(len=64) :: iomsg, Writable, Form
    integer           :: iostat
    logical           :: Opened
    INQUIRE(&
         UNIT  = unit    , &
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
       call SymmetryOperationIOFmt(Op,unit,"WRITE")
    case("UNFORMATTED")
       call SymmetryOperationIOUnf(Op,unit,"WRITE")
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
    !
  end subroutine SymmetryOperationWriteUnit

  subroutine SymmetryOperationReadUnit( Op, unit )
    class( SymmetryOperation ), intent(out) :: Op
    integer                   , intent(in)    :: unit
    character(len=64) :: iomsg, Readable, Form
    integer           :: iostat
    logical           :: Opened
    INQUIRE(&
         UNIT  = unit    , &
         OPENED= Opened  , &
         READ  = Readable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Output Unit is not open")
    if( trim(Readable) /= "YES" ) call Assert("Output Unit can't be read")
    !
    select case (trim(FORM))
    case("FORMATTED")
       call SymmetryOperationIOFmt(Op,unit,"READ")
    case("UNFORMATTED")
       call SymmetryOperationIOUnf(Op,unit,"READ")
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
    !
  end subroutine SymmetryOperationReadUnit

  subroutine SymmetryOperationIOFmt(Op,unit,action)
    class( SymmetryOperation )   :: Op
    integer         , intent(in) :: unit
    character(len=*), intent(in) :: action
    !
    character(len=64) :: iomsg
    integer           :: iostat
    integer :: iRow, iCol
    character(len=64) :: Name
    real(kind(1d0)) :: R(3,3)
    if(action.is."WRITE")then
       write(unit,"(a,*(1x,f4.0))",iostat=iostat,iomsg=iomsg) Op%Name,((Op%R(iRow,iCol),iRow=1,3),iCol=1,3)
       if(iostat/=0)call Assert(iomsg)
    elseif(action.is."READ")then
       read (unit,"(a,*(1x,f4.0))",iostat=iostat,iomsg=iomsg) Name,((R(iRow,iCol),iRow=1,3),iCol=1,3)
       if(iostat/=0)call Assert(iomsg)
       call Op%init(trim(adjustl(Name)),R)
    else
       call Assert("Invalid action")
    endif
  end subroutine SymmetryOperationIOFmt
  
  subroutine SymmetryOperationIOUnf(Op,unit,action)
    class( SymmetryOperation )   :: Op
    integer         , intent(in) :: unit
    character(len=*), intent(in) :: action
    !
    character(len=64) :: iomsg
    integer           :: iostat
    integer :: iRow, iCol, NameLength
    character(len=:), allocatable :: Name
    real(kind(1d0)) :: R(3,3)
    if(action=="WRITE")then
       write(unit,iostat=iostat,iomsg=iomsg) len(Op%Name)
       if(iostat/=0)call Assert(iomsg)
       write(unit,iostat=iostat,iomsg=iomsg) Op%Name, ((Op%R(iRow,iCol),iRow=1,3),iCol=1,3)
       if(iostat/=0)call Assert(iomsg)
    elseif(action=="READ" )then
       read (unit,iostat=iostat,iomsg=iomsg) NameLength
       if(iostat/=0)call Assert(iomsg)
       allocate(Name,source=REPEAT("#",NameLength))
       read (unit,iostat=iostat,iomsg=iomsg) Name, ((R(iRow,iCol),iRow=1,3),iCol=1,3)
       if(iostat/=0)call Assert(iomsg)
       call Op%init(Name,R)
    else
       call Assert("Invalid action")
    endif
  end subroutine SymmetryOperationIOUnf
  
  subroutine SymmetryOperationFinal( Op )
    type( SymmetryOperation ), intent(inout) :: Op
    call Op%Free
  end subroutine SymmetryOperationFinal



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Class Irrep Methods
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  subroutine ClassIrrepAssign( IrrepOut, IrrepIn )
    class( ClassIrrep ), intent(inout) :: IrrepOut
    class( ClassIrrep ), intent(in)    :: IrrepIn
    call IrrepOut%Init( IrrepIn%Name, IrrepIn%Characters )
    IrrepOut%Group => IrrepIn%Group
  end subroutine ClassIrrepAssign

  subroutine ClassIrrepInit( Irrep, Name, Characters, Group )
    class( ClassIrrep )          , intent(inout) :: Irrep
    character(len=*)             , intent(in)    :: Name
    real(kind(1d0))              , intent(in)    :: Characters(:)
    class( ClassGroup ), optional, target, intent(in) :: Group
    call Irrep%free
    allocate(Irrep%Name,source=Name)
    Irrep%NClasses=size(Characters)
    allocate(Irrep%Characters,source=Characters)
    Irrep%Group=>NULL()
    if(present(Group)) Irrep%Group => Group
  end subroutine ClassIrrepInit

  subroutine ClassIrrepFree( irrep )
    class( ClassIrrep ), intent(inout) :: irrep
    if(allocated(irrep%Name      )) deallocate(irrep%Name)
    if(allocated(irrep%Characters)) deallocate(irrep%Characters)
    irrep%NClasses=0
  end subroutine ClassIrrepFree

  subroutine ClassIrrepFinal( irrep )
    type( ClassIrrep ), intent(inout) :: irrep
    call irrep%free
  end subroutine ClassIrrepFinal

  subroutine ClassIrrepShow( irrep, unit )
    class( ClassIrrep ), intent(in) :: irrep
    integer,  optional , intent(in) :: unit
    integer :: outunit, iClass
    outunit=OUTPUT_UNIT
    if(present(unit)) outunit=unit
    write(outunit,"(a3,*(x,f4.0))") irrep%Name, (irrep%Characters(iClass),iClass=1,irrep%NClasses)
  end subroutine ClassIrrepShow

  subroutine ClassIrrepWriteUnit( irrep, unit )
    class( ClassIrrep ), intent(in) :: irrep
    integer            , intent(in) :: unit
    character(len=64) :: iomsg, Writable, Form
    integer           :: iostat
    logical           :: Opened
    INQUIRE(&
         UNIT  = unit    , &
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
       call ClassIrrepIOFmt(irrep,unit,"WRITE")
    case("UNFORMATTED")
       call ClassIrrepIOUnf(irrep,unit,"WRITE")
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
  end subroutine ClassIrrepWriteUnit

  subroutine ClassIrrepReadUnit( irrep, unit )
    class( ClassIrrep ), intent(out) :: irrep
    integer            , intent(in)  :: unit
    character(len=64) :: iomsg, Readable, Form
    integer           :: iostat
    logical           :: Opened
    INQUIRE(&
         UNIT  = unit    , &
         OPENED= Opened  , &
         READ  = Readable, &
         FORM  = Form    , &
         IOSTAT= iostat  , &
         IOMSG = iomsg     )
    if(iostat/=0) call Assert(iomsg)
    if( .not. Opened            ) call Assert("Output Unit is not open")
    if( trim(Readable) /= "YES" ) call Assert("Output Unit can't be Read")
    !
    select case (trim(FORM))
    case("FORMATTED")
       call ClassIrrepIOFmt(irrep,unit,"READ")
    case("UNFORMATTED")
       call ClassIrrepIOUnf(irrep,unit,"READ")
    case DEFAULT
       call Assert("Invalid Output Unit Format")
    end select
  end subroutine ClassIrrepReadUnit

  subroutine ClassIrrepIOFmt(irrep,unit,action)
    class( ClassIrrep )          :: irrep
    integer         , intent(in) :: unit
    character(len=*), intent(in) :: action
    !
    character(len=64) :: iomsg
    integer           :: iostat, iClass, NClasses
    character(len=64) :: Name
    real(kind(1d0)) , allocatable :: Characters(:)
    if(action.is."WRITE")then
       write(unit,"(a,1x,i4)",iostat=iostat,iomsg=iomsg) irrep%Name, irrep%NClasses
       if(iostat/=0)call Assert(iomsg)
       write(unit,"(*(1x,f4.0))",iostat=iostat,iomsg=iomsg)(irrep%Characters(iClass),iClass=1,irrep%NClasses)
       if(iostat/=0)call Assert(iomsg)
    elseif(action.is."READ")then
       read(unit,*,iostat=iostat,iomsg=iomsg) Name, NClasses
       if(iostat/=0)call Assert(iomsg)
       allocate(Characters(NClasses))
       read(unit,*,iostat=iostat,iomsg=iomsg)(Characters(iClass),iClass=1,NClasses)
       if(iostat/=0)call Assert(iomsg)
       call irrep%init(trim(adjustl(Name)),Characters)
    else
       call Assert("Invalid action")
    endif
  end subroutine ClassIrrepIOFmt

  subroutine ClassIrrepIOUnf(irrep,unit,action)
    class( ClassIrrep )          :: irrep
    integer         , intent(in) :: unit
    character(len=*), intent(in) :: action
    !
    character(len=64) :: iomsg
    integer           :: iostat, iClass, NClasses, NameLength
    character(len=:), allocatable :: Name
    real(kind(1d0)) , allocatable :: Characters(:)
    if(action.is."WRITE")then
       write(unit,iostat=iostat,iomsg=iomsg) len(irrep%Name), irrep%NClasses
       if(iostat/=0)call Assert(iomsg)
       write(unit,iostat=iostat,iomsg=iomsg) irrep%Name, (irrep%Characters(iClass),iClass=1,irrep%NClasses)
       if(iostat/=0)call Assert(iomsg)
    elseif(action.is."READ")then
       read(unit,iostat=iostat,iomsg=iomsg) NameLength, NClasses
       if(iostat/=0)call Assert(iomsg)
       allocate(Name,source=REPEAT("#",NameLength))
       allocate(Characters(NClasses))
       read(unit,iostat=iostat,iomsg=iomsg) Name, (Characters(iClass),iClass=1,NClasses)
       if(iostat/=0)call Assert(iomsg)
       call irrep%init(Name,Characters)
    else
       call Assert("Invalid action")
    endif
  end subroutine ClassIrrepIOUnf

  logical function ClassIrrepNameIs(irrep,Name) result(IsName)
    class( ClassIrrep ), intent(in) :: irrep
    character(len=*)   , intent(in) :: Name
    IsName = Name .is. irrep%Name
  end function ClassIrrepNameIs

  subroutine ClassIrrepSetGroup(irrep,Group) 
    class( ClassIrrep )        , intent(inout) :: irrep
    class( ClassGroup ), target, intent(in)    :: Group
    irrep%Group => Group
  end subroutine ClassIrrepSetGroup

  function ClassIrrepGetName(irrep) result(Name)
    class( ClassIrrep ), intent(in) :: irrep
    character(:), allocatable   :: Name
    allocate(Name,source=irrep%Name)
  end function ClassIrrepGetName

  function ClassIrrepGetIdStrn(irrep) result(IdStrn)
    class( ClassIrrep ), intent(in) :: irrep
    character(len=:), allocatable   :: IdStrn
    allocate(IdStrn,source=irrep%Group%GetName()//"."//irrep%Name)
  end function ClassIrrepGetIdStrn

  function ClassIrrepGetGroupName(irrep) result(GroupName)
    class( ClassIrrep ), intent(in) :: irrep
    character(len=:), allocatable   :: GroupName
    allocate(GroupName,source=irrep%Group%GetName())
  end function ClassIrrepGetGroupName

  function ClassIrrepGetGroup(irrep) result(Group)
    class( ClassIrrep ), intent(in) :: irrep
    type( ClassGroup ), pointer    :: Group
    Group => irrep%Group
  end function ClassIrrepGetGroup

  function ClassIrrepOTimes(irrep1,irrep2) result(irrepProd)
    class( ClassIrrep ), intent(in) :: irrep1, irrep2
    class( ClassIrrep ), pointer    :: irrepProd
    real(kind(1d0)), parameter :: THRESHOLD = 1.d-5
    real(kind(1d0)), allocatable :: vec(:)
    class( ClassGroup ), pointer :: Group
    integer :: iClass
    real(kind(1d0)) :: ScalarProduct
    !if(.not. (irrep1%Group%GetName .is. irrep2%Group%GetName ) ) call Assert("Incompatible irrep")
    allocate(vec(irrep1%NClasses))
    vec = irrep1%Characters * irrep2%Characters
    Group => irrep1%Group
    irrepProd => NULL()
    do iClass = 1, irrep1%NClasses
       ScalarProduct = sum( vec * Group%Irreps(iClass)%Characters ) / Group%Order
       if( abs(ScalarProduct - 1.d0) < THRESHOLD )then
          irrepProd => Group%Irreps(iClass)
          exit
       endif
    enddo
  end function ClassIrrepOTimes


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Class Group Methods
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  subroutine ClassGroupInitDefault( Group, Name )
    use ModuleString
    class( ClassGroup ), intent(out) :: Group
    character(len=*)   , intent(in)  :: Name
    !
    character(len=:), allocatable :: UCName
    call Group%free
    allocate(UCName,source=Name)
    call SetStringToUppercase( UCName )
    call InitModuleGroups
    select case ( UCName )
    case( "C1" )
       call InitGroupC1 ( Group )
    case( "C2" )
       call InitGroupC2 ( Group )
    case( "CS" )
       call InitGroupCs ( Group )
    case( "CI" )
       call InitGroupCi ( Group )
    case( "C2V" )
       call InitGroupC2v( Group )
    case( "C2H" )
       call InitGroupC2h( Group )
    case( "D2" )
       call InitGroupD2 ( Group )
    case( "D2H" )
       call InitGroupD2h( Group )
    case default
       call Assert("Unrecognized group "//trim(Name))
    end select
    Group%Initialized=.TRUE.
  end subroutine ClassGroupInitDefault

  subroutine ClassGroupInitNew( Group, Name, Elements, Irreps )
    use ModuleString
    class( ClassGroup )       , intent(out) :: Group
    character(len=*)          , intent(in)  :: Name
    class( SymmetryOperation ), intent(in)  :: Elements(:)
    class( ClassIrrep )       , intent(in)  :: Irreps(:)
    !
    integer :: iOrd
    call InitModuleGroups
    call Group%free()
    allocate(Group%Name,source=Name)
    Group%Order=size(Elements)
    allocate(Group%Elements(Group%Order))
    allocate(Group%Irreps(Group%Order))
    do iOrd = 1, Group%Order
       Group%Elements(iOrd)=Elements(iOrd)
       Group%Irreps(iOrd)=Irreps(iOrd)
       call Group%Irreps(iOrd)%SetGroup(Group)
    enddo
    Group%Initialized=.TRUE.
  end subroutine ClassGroupInitNew
  
  subroutine ClassGroupFree( Group )
    class( ClassGroup ), intent(inout) :: Group
    if(allocated(Group%Name    )) deallocate(Group%Name)
    if(allocated(Group%Elements)) deallocate(Group%Elements)
    Group%Order=0
    Group%Initialized=.FALSE.
  end subroutine ClassGroupFree
  
  subroutine ClassGroupFinal( Group )
    type( ClassGroup ), intent(inout) :: Group
    call Group%Free
  end subroutine ClassGroupFinal
  
  subroutine ClassGroupShow( Group, unit )
    class( ClassGroup ), intent(in) :: Group
    integer, optional  , intent(in) :: unit
    !
    integer :: iOrd, outunit
    outunit=OUTPUT_UNIT
    if(present(unit))outunit=unit
    if(.not.Group%IsInitialized())call Assert("Trying to show uninitialized group")
    write(outunit,"(a)"        ) " Name     : "//Group%Name
    write(outunit,"(a,i2)"     ) " Order    :", Group%Order
    write(outunit,"(a,*(1x,a))") " Elements :",(Group%Elements( iOrd )%Name, iOrd = 1, Group%Order )
    do iOrd = 1, Group%Order
       write(outunit,*)
       call Group%Elements(iOrd)%show(unit)
    enddo
    write(outunit,*)
    do iOrd = 1, Group%Order
       call Group%Irreps(iOrd)%show(unit)
    enddo
    write(outunit,*)
  end subroutine ClassGroupShow

  logical function ClassGroupDefinesIrrepName( Group, IrrepName ) result( IrrepIsLegit )
    class( ClassGroup ), intent(in) :: Group
    character(len=*)   , intent(in) :: IrrepName
    !
    integer :: iClass
    IrrepIsLegit = .False.
    do iClass=1, Group%Order
       IrrepIsLegit = Group%Irreps(iClass)%NameIs(IrrepName)
       if(IrrepIsLegit)exit
    enddo
  end function ClassGroupDefinesIrrepName

  function ClassGroupGetIrrepPtrFromName( Group, irrepName ) result( irrepPtr )
    class( ClassGroup ), target , intent(in)  :: Group
    character(len=*)            , intent(in)  :: irrepName
    type( ClassIrrep ), pointer :: irrepPtr
    !
    integer :: iClass
    irrepPtr => NULL()
    do iClass=1, Group%Order
       if( Group%Irreps( iClass )%NameIs( IrrepName ) ) IrrepPtr => Group%Irreps( iClass )
    enddo
  end function ClassGroupGetIrrepPtrFromName

  function ClassGroupGetIrrepPtrFromIndex( Group, iIrrep ) result( irrepPtr )
    class( ClassGroup ), target , intent(in)  :: Group
    integer                     , intent(in)  :: iIrrep
    type( ClassIrrep ), pointer :: irrepPtr
    !
    integer :: iClass
    IrrepPtr => Group%Irreps( iIrrep )
  end function ClassGroupGetIrrepPtrFromIndex


  function ClassGroupGetIrrepName( Group , iIrrep ) result( irrepName )
    class( ClassGroup ), target , intent(in)  :: Group
    integer                     , intent(in)  :: iIrrep
    character(len=:)           , allocatable  :: irrepName
    allocate(irrepName,source=Group%Irreps(iIrrep)%getName())
  end function ClassGroupGetIrrepName


  function ClassGroupGetTotSymIrrepPtr( Group ) result( irrepPtr )
    class( ClassGroup ), target , intent(in)  :: Group
    type( ClassIrrep ), pointer :: irrepPtr
    IrrepPtr => Group%Irreps( 1 )
  end function ClassGroupGetTotSymIrrepPtr

  integer function ClassGroupGetIrrepIndex( Group, Irrep ) result( Index )
    !
    class(ClassGroup), intent(in) :: Group
    class(ClassIrrep), intent(in) :: Irrep
    !
    integer :: i
    !
    if ( .not. Group%IsInitialized() ) then
       call Assert( "The point group has not been initialized, imposible to get the irreducible representation index." )
    end if
    !
    do i = 1, Group%Order
       if ( Irrep%GetName() == Group%Irreps(i)%GetName() ) then
          Index = i
          exit
       end if
    end do
    !
  end function ClassGroupGetIrrepIndex


  function ClassGroupIsInitialized( Group ) result( Initialized )
    class( ClassGroup ), intent(in) :: Group
    logical                         :: Initialized
    Initialized = Group%Initialized
  end function ClassGroupIsInitialized

  function ClassGroupGetName( Group ) result( Name )
    class( ClassGroup ), target , intent(in)  :: Group
    character(len=:), allocatable             :: Name
    if(Group%IsInitialized())then
       allocate(Name,source=Group%Name)
    else
       allocate(Name,source="")
    endif
  end function ClassGroupGetName  

  integer function ClassGroupGetnIrreps( Group ) result( nIrreps )
    class( ClassGroup ), target , intent(in)  :: Group
    nIrreps = 0
    if(Group%IsInitialized()) nIrreps = size( Group%Irreps )
  end function ClassGroupGetnIrreps

  function ClassGroupGetIrrepList( Group ) result( irrepv )
    class( ClassGroup ), target , intent(in)  :: Group
    type( ClassIrrep ), pointer, dimension(:) :: irrepv
    if(Group%IsInitialized())then
       Irrepv => Group%Irreps
    else
       Irrepv => NULL()
    endif
  end function ClassGroupGetIrrepList

  function ClassGroupGetTotallySymmetricIrrep( Group ) result( irrep )
    class( ClassGroup ), target , intent(in)  :: Group
    type( ClassIrrep ), pointer :: irrep
    if(Group%IsInitialized())then
       Irrep => Group%Irreps(1)
    else
       Irrep => NULL()
    endif
  end function ClassGroupGetTotallySymmetricIrrep

  !> Returns a pointer to the irreducible representation of x
  function ClassGroupGetxIrrep( Group ) result( irrep )
    class( ClassGroup ), target , intent(in)  :: Group
    type( ClassIrrep ), pointer :: irrep
    if(Group%IsInitialized())then
       if( Group%GetName() .is. "C1" ) irrep => Group%GetIrrep("A"  )
       if( Group%GetName() .is. "Cs" ) irrep => Group%GetIrrep("Ap" )
       if( Group%GetName() .is. "C2" ) irrep => Group%GetIrrep("B"  )
       if( Group%GetName() .is. "Ci" ) irrep => Group%GetIrrep("Au" )
       if( Group%GetName() .is. "C2v") irrep => Group%GetIrrep("B1" )
       if( Group%GetName() .is. "C2h") irrep => Group%GetIrrep("Bu" )
       if( Group%GetName() .is. "D2" ) irrep => Group%GetIrrep("B3" )
       if( Group%GetName() .is. "D2h") irrep => Group%GetIrrep("B3u")
    else
       Irrep => NULL()
    endif
  end function ClassGroupGetxIrrep

  !> Returns a pointer to the irreducible representation of y
  function ClassGroupGetyIrrep( Group ) result( irrep )
    class( ClassGroup ), target , intent(in)  :: Group
    type( ClassIrrep ), pointer :: irrep
    if(Group%IsInitialized())then
       if( Group%GetName() .is. "C1" ) irrep => Group%GetIrrep("A"  )
       if( Group%GetName() .is. "Cs" ) irrep => Group%GetIrrep("Ap" )
       if( Group%GetName() .is. "C2" ) irrep => Group%GetIrrep("B"  )
       if( Group%GetName() .is. "Ci" ) irrep => Group%GetIrrep("Au" )
       if( Group%GetName() .is. "C2v") irrep => Group%GetIrrep("B2" )
       if( Group%GetName() .is. "C2h") irrep => Group%GetIrrep("Bu" )
       if( Group%GetName() .is. "D2" ) irrep => Group%GetIrrep("B2" )
       if( Group%GetName() .is. "D2h") irrep => Group%GetIrrep("B2u")
    else
       Irrep => NULL()
    endif
  end function ClassGroupGetyIrrep

  !> Returns a pointer to the irreducible representation of z
  function ClassGroupGetzIrrep( Group ) result( irrep )
    class( ClassGroup ), target , intent(in)  :: Group
    type( ClassIrrep ), pointer :: irrep
    if(Group%IsInitialized())then
       if( Group%GetName() .is. "C1" ) irrep => Group%GetIrrep("A"  )
       if( Group%GetName() .is. "Cs" ) irrep => Group%GetIrrep("App")
       if( Group%GetName() .is. "C2" ) irrep => Group%GetIrrep("A"  )
       if( Group%GetName() .is. "Ci" ) irrep => Group%GetIrrep("Au" )
       if( Group%GetName() .is. "C2v") irrep => Group%GetIrrep("A1" )
       if( Group%GetName() .is. "C2h") irrep => Group%GetIrrep("Au" )
       if( Group%GetName() .is. "D2" ) irrep => Group%GetIrrep("B1" )
       if( Group%GetName() .is. "D2h") irrep => Group%GetIrrep("B1u")
    else
       Irrep => NULL()
    endif
  end function ClassGroupGetzIrrep

  function ClassGroupGetMonomialIrrep( Group, ix, iy, iz ) result( irrep )
    class( ClassGroup ), target , intent(in) :: Group
    integer                     , intent(in) :: ix, iy, iz
    type( ClassIrrep ), pointer              :: irrep, irrep2
    irrep => Group%GetTotallySymmetricIrrep()
    if(mod(ix,2)==1)then
       irrep2 => irrep * Group%GetxIrrep()
       irrep  => irrep2
    endif
    if(mod(iy,2)==1)then
       irrep2 => irrep * Group%GetyIrrep()
       irrep  => irrep2
    endif
    if(mod(iz,2)==1)then
       irrep2 => irrep * Group%GetzIrrep()
       irrep  => irrep2
    endif
  end function ClassGroupGetMonomialIrrep

  
  logical function ValidIrreps( BraIrrep, OpIrrep, KetIrrep ) result(ValidIrr)
    class(ClassIrrep), target, intent(in) :: BraIrrep
    class(ClassIrrep), target, intent(in) :: OpIrrep
    class(ClassIrrep), target, intent(in) :: KetIrrep
    type(ClassIrrep), pointer :: AuxIrrep
    ValidIrr = .false.
    AuxIrrep => OpIrrep * KetIrrep
    if ( BraIrrep%NameIs( AuxIrrep%GetName() ) ) ValidIrr = .true.
  end function ValidIrreps
  

end module ModuleGroups
