! {{{ Detailed description

!> \mainpage Program TabulateOrbitals tabulates the orbitals on a grid that can 
!! be loaded from paraview
!! 
!! Synopsis:
!! ---------
!!
!!     <Program Name> <mandatory run-time parameters (RTP)> [<optional RTP>]
!!
!! ___
!! Description:
!! ------------
!!
!! Input parameters:      {#Input_Parameters}
!! =================
!! [...Input](@ref ...Input) as specified in the command line.
!!
!> \file
!!
!!
! }}}
program ProgramTabulateOrbitals

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO

  !.. Application Modules
  use ModuleLibCint
  use ModuleQCBasis
  use ModuleMolecularGeometry
  use ModuleDaltonLibCint

  implicit none

  character(len=*), parameter :: GRID_FILE = "Grid"
  
!!$  !.. Local parameters
!!$  integer                       :: uid, ish, ish1, ish2, i, j, di, dj, size
!!$  character(len=:), allocatable :: sish
!!$  real(kind(1d0)), allocatable  :: buf1e_tot(:,:), orb(:)
!!$  real(kind(1d0)), allocatable  :: buf1e(:,:,:)
!!$
  type(ClassMolecularGeometry):: MolGeom
  type(ClassQCAtomicSetBasis) :: atomicSet
  type(ClassCintData)         :: CintData
  
  !.. Grid parameters
  integer                      :: ipt, npt
  real(kind(1d0))              :: x, y, z, xv(3)
  real(kind(1d0)), allocatable :: grid(:,:)

  call MolGeom.parseFile("MOLECULE.INP")
  call atomicSet.parseFile( "cc-PVTZ-Dalton.bas" )
  call DaltonToCintData( MolGeom, atomicSet, CintData )

  call tmpGenerateGrid( GRID_FILE )
  call LoadGrid( GRID_FILE, npt, grid )

  !.. Load the orbitals
  

  
  Size = CintData.GetMaxShellSize()
  allocate(orb(Size))
  !
  do ish = 1, CintData.GetNshells()
     
     di = CintData.GetShellSize(ish)
     sish=AlphabeticNumber(ish,100,"0")
     call OpenFile( "Orb"//sish, uid, "write", "formatted" )
     do ipt = 1, npt
        !
        xv=grid(:,ipt)
        call plotOrbitals(CintData,ish,xv,orb)
        !
        if( ipt > 1 )then
           if( xv(1) /= grid(1,ipt-1) ) write(uid,*) 
        endif
        write(uid,"(*(x,e14.6))") xv(1), xv(2), xv(3), ( orb(i), i = 1, di )
        !
     enddo
     close(uid)
  enddo


  stop


  Size = CintData.GetBasisSize()
  allocate(buf1e_tot(Size,Size))
  !
  call CintData.SetCharge( TEST_CHARGE, 1 )
  !
  call OpenFile( "Pot", uid, "write", "formatted" )
  do ipt = 1, npt
     !
     x = grid(1,ipt)
     y = grid(2,ipt)
     z = grid(3,ipt)
     !
     call CintData.SetCoords( TEST_CHARGE, x, y, z )
     call CintData.int1e_nuc( buf1e_tot )
     buf1e_tot = - buf1e_tot !.. Must revert the sign because the nucleus is positive
     !
     if( ipt > 1 )then
        if( x /= grid(1,ipt-1) ) write(uid,*) 
     endif
     write(uid,"(*(x,e14.6))") x, y, z, ( ( buf1e_tot(i,j), i = 1, Size ), j = 1, Size )
     !
  enddo
  close(uid)

  


  stop
  
  !.. Evaluate the potential of the product of contracted
  !   functions in shells 1 and 4.
  !..
  ish1 = 2
  ish2 = 5
  di = CintData.GetShellSize( ish1 )
  dj = CintData.GetShellSize( ish2 )
  allocate(buf1e(di,dj,1))
  !
  call CintData.SetCharge( TEST_CHARGE, 1 )
  !
  call OpenFile( "Pot", uid, "write", "formatted" )
  do ipt = 1, npt
     !
     x = grid(1,ipt)
     y = grid(2,ipt)
     z = grid(3,ipt)
     !
     call CintData.SetCoords( TEST_CHARGE, x, y, z )
     call CintData.int1e_nuc( ish1, ish2, buf1e )
     buf1e = -buf1e
     !
     if( ipt > 1 .and. x /= grid(1,ipt-1) ) write(uid,*) 
     write(uid,"(*(x,e14.6))") x, y, z, ( ( buf1e(i,j,1), i = 1, di ), j = 1, dj )
     !
  enddo
  close(uid)

  deallocate(buf1e)

contains

  
  subroutine LoadGrid(FileName,npts,grid)
    character(len=*)            , intent(in)  :: FileName
    integer                     , intent(out) :: npts
    real(kind(1d0)), allocatable, intent(out) :: grid(:,:)
    integer :: uid,iostat,ipt, ix
    real(kind(1d0)) :: w
    call OpenFile(FileName,uid,"read","formatted")
    npts=0
    do
       read(uid,*,iostat=iostat) w,w,w
       if(iostat/=0)exit
       npts=npts+1
    enddo
    allocate(grid(3,npts))
    rewind(uid)
    do ipt=1,npts
       read(uid,*) (grid(ix,ipt),ix=1,3)
    enddo
    close(uid)
  end subroutine LoadGrid
  

  subroutine tmpGenerateGrid( FileName )
    character(len=*), intent(in) :: FileName

    !---------- Grid data ----------------------
    integer                     :: ix
    integer        , parameter  :: nxpt = 101
    real(kind(1d0)), parameter  :: xmin = -4.d0
    real(kind(1d0)), parameter  :: xmax =  4.d0
    !
    integer                     :: iz
    integer        , parameter  :: nzpt = 101
    real(kind(1d0)), parameter  :: zmin = -4.d0
    real(kind(1d0)), parameter  :: zmax =  4.d0
    !
    real(kind(1d0))             :: x, y, z
    !-------------------------------------------
    integer :: uid

    call OpenFile( FileName, uid, "write", "formatted" )
    y = 0.d0
    do ix = 1, nxpt
       x = xmin + (xmax - xmin) * dble(ix-1)/dble(nxpt-1)
       do iz = 1, nzpt
          z = zmin + (zmax - zmin) * dble(iz-1)/dble(nzpt-1)
          ipt = ipt + 1
          write(uid,"(*(x,d24.16))") x,y,z
       enddo
    enddo
    close(uid)
end subroutine TmpGenerateGrid


end program ProgramTabulateOrbitals

