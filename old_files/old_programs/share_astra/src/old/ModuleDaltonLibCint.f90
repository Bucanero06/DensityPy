module ModuleDaltonLibCint
  
  use, intrinsic :: ISO_FORTRAN_ENV
  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO
  use ModuleConstants

  use ModuleMolecularGeometry
  use ModuleQCBasis
  use ModuleLibCint

  implicit none

  private
  
  public :: DaltonToCintData
  public :: plotOrbitals

  integer, public, protected :: TEST_CHARGE

  integer        , parameter :: DELTA_ORB_FLAG = 1
  real(kind(1d0)), parameter :: DELTA_ALPHA    = 1.d8
  integer                    :: DELTA_ORB
  integer                    :: DELTA_SHELL
  

contains
  

  subroutine DaltonToCintData( MolGeom, atomicSet, CintData )
    type(ClassMolecularGeometry), intent(in)  :: MolGeom
    type(ClassQCAtomicSetBasis) , intent(in)  :: atomicSet
    type(ClassCintData)         , intent(out) :: CintData

    integer :: iAtom, nat, totNsh, icharge, AtomicNumber, l, ish
    integer :: nExp, nCon
    real(kind(1d0)), allocatable :: vexp(:), mcoe(:,:)
    real(kind(1d0)) :: xv(3)

    call CintData.free()

    !.. Count the total number of shells
    nat = MolGeom.GetNat()
    totNsh = 0
    do iAtom = 1, nat

       !.. Identify the atom 
       AtomicNumber = MolGeom.GetAtNumb( iAtom )
       
       !.. Accumulates the number of shell for the current atom
       totNsh = totNsh + atomicSet.GetNsh( AtomicNumber )

    enddo

    !.. Initialize Cint
    call CintData.init( nat + 1 + DELTA_ORB_FLAG, totNsh + DELTA_ORB_FLAG )

    !.. Fills Cint with the contraction parameters
    ish = 0
    do iAtom = 1, nat

       !.. Set the charge
       !icharge = int( MolGeom.GetCharge( iAtom ) + 0.1d0 )
       icharge = 0 !.. To compute the potential, we need all nuclei to be neutral
       call CintData.setCharge( iAtom, icharge ) 

       !.. Set the coordinates
       call MolGeom.GetCoords( iAtom, xv )
       call CintData.setCoords( iAtom, xv(1), xv(2), xv(3) )

       !.. Set the contraction coefficients
       AtomicNumber = MolGeom.GetAtNumb( iAtom )

       do l = 0, atomicSet.GetLmax( AtomicNumber )
          ish = ish + 1
          call atomicSet.GetCont( AtomicNumber, l, nExp, vExp, nCon, mCoe )
          call CintData.setShell(   iAtom, ish, l, nExp, vExp, nCon, mCoe )
       enddo

    enddo

    !.. Set the test charge
    TEST_CHARGE = nat+1
    call CintData.setCharge( TEST_CHARGE, 0 )
    call CintData.setCoords( TEST_CHARGE, 0.d0, 0.d0, 0.d0 )
    !.. Set the test charge

    !.. Set the delta-like orbital
    if(DELTA_ORB_FLAG==1)then
       DELTA_ORB   = nat+2
       DELTA_SHELL = totNsh+1 
       deallocate(vExp);allocate(vExp(1))
       vExp=DELTA_ALPHA
       deallocate(mCoe);allocate(mCoe(1,1))
       !
       !.. Normalize the function in such a way that it supposedly 
       !   becomes a representation of the delta function
       !.. 
       mCoe=exp(0.75d0*log(DELTA_ALPHA)-0.25d0*log(2.d0*PI))
       call CintData.setCharge( DELTA_ORB, 0 )
       call CintData.setCoords( DELTA_ORB, 0.d0, 0.d0, 0.d0 )
       call CintData.setShell ( DELTA_ORB, totNsh+1, 0, 1, vExp, 1, mCoe )
    endif

  end subroutine DaltonToCintData


  subroutine plotOrbitals(CintData,ish,xv,dbuf)

    type(ClassCintData), intent(inout) :: CintData
    integer            , intent(in)    :: ish
    real(kind(1d0))    , intent(in)    :: xv(3)
    real(kind(1d0))    , intent(out)   :: dbuf(:)
    
    logical        , save              :: FIRST_CALL = .TRUE.
    real(kind(1d0)), save, allocatable :: orbv(:,:,:)
    integer :: ds, iShell, info
    
    if(FIRST_CALL)then
       ds = CintData.GetMaxShellSize()
       allocate(orbv(ds,1,1))
       orbv=0.d0
       FIRST_CALL = .FALSE.
    endif

    call CintData.SetCoords( DELTA_ORB, xv(1), xv(2), xv(3) )
    call CintData.int1e( "ovlp", ish, DELTA_SHELL, orbv )

    ds = CintData.GetShellSize( ish )
    
    dbuf(1:ds) = orbv(1:ds,1,1)
    
  end subroutine plotOrbitals


end module ModuleDaltonLibCint
