program TabulatePotential

  !.. System Modules
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

  !.. Run-time parameters
  character(len=:), allocatable :: MoleculeInputFile
  character(len=:), allocatable :: GTOBasisInputFile
  character(len=:), allocatable :: LUCIA_CMOAO_NOSYM
  character(len=:), allocatable :: MOAO_InfoFile

  !.. Local parameters
  integer                       :: uid, iostat, ish, ish1, ish2, i, j, di, dj, size
  character(len=:), allocatable :: sish
  character(len=500)            :: iomsg
  real(kind(1d0)) , allocatable :: buf1e_tot(:,:), orb(:)
  real(kind(1d0)) , allocatable :: buf1e(:,:,:)

  type(ClassMolecularGeometry)  :: MolGeom
  type(ClassQCAtomicSetBasis)   :: atomicSet
  type(ClassCintData)           :: CintData
  
  !.. Grid parameters
  integer                       :: ipt, npt
  real(kind(1d0))               :: x, y, z, xv(3)
  real(kind(1d0)), allocatable  :: grid(:,:)
  
  !.. Molecular Orbitals
  integer :: nsym, isym
  integer :: nsymAO(8)
  integer :: nsymMO(8)
  integer :: nAO, nMO, iAO, iMO

  real(kind(1d0)), allocatable :: dmat(:,:)

  type ClassDMat2
     integer :: n1,n2
     real(kind(1d0)), allocatable :: A(:,:)
  end type ClassDMat2
  type(ClassDMat2), allocatable :: CAOMO(:)

  call GetRunTimeParameters( &
       MoleculeInputFile   , &
       GTOBasisInputFile   , &
       MOAO_InfoFile       , &
       LUCIA_CMOAO_NOSYM   )

  call MolGeom.parseFile( MoleculeInputFile )
  call atomicSet.parseFile( GTOBasisInputFile )
  call DaltonToCintData( MolGeom, atomicSet, CintData )

  call MolGeom.show()
  call atomicSet.show()
  call CintData.Print( OUTPUT_UNIT )

  !.. Load information on the size of the molecular orbitals
  nsymAO=0
  nsymMO=0
  open(newunit = uid          , &
       file    = MOAO_InfoFile, &
       form    ="formatted"   , &
       status  ="old"         , &
       action  ="read"        , &
       iostat  = iostat       , &
       iomsg   = iomsg        )
  if( iostat /= 0 )then
     write(ERROR_UNIT,"(a)") trim(iomsg)
     stop
  end if
  read(uid,*) nsym
  read(uid,*) (nsymAO(isym),isym=1,nsym)
  read(uid,*) (nsymMO(isym),isym=1,nsym)
  nAO = sum(nsymAO)
  nMO = sum(nsymMO)
  if( nMO /= nAO )then
     write(ERROR_UNIT,"(a,i0)") "#MO = ",nMO
     write(ERROR_UNIT,"(a,i0)") "#AO = ",nAO
     write(ERROR_UNIT,"(a)"   ) "#MO /= #AO : STOP FORCED"
     STOP
  endif
  !
  !*** MUST HACK FROM HERE THE ORDERING OF THE ORBITALS FROM 
  !    FORT.91 WHERE THERE IS THE LIST OF THE ORBITALS
  !***

  close(uid)


  !.. Load the raw molecular orbitals
  open(newunit = uid              , &
       file    = LUCIA_CMOAO_NOSYM, &
       form    ="unformatted"     , &
       status  ="old"             , &
       action  ="read"            , &
       iostat  = iostat           , &
       iomsg   = iomsg            )
  if( iostat /= 0 )then
     write(ERROR_UNIT,"(a)") trim(iomsg)
     stop
  end if
  allocate(dmat(nAO,nMO))
  dmat=0.d0
  read(uid) ( ( dmat(iAO,iMO), iAO = 1, nAO ), iMO = 1, nMO )
  close(uid)

  
  
  !.. Convert the MO to a set of matrices separated by symmetry
  allocate( CAOMO( nsym ) )
  iMO = 0
  do isym = 1, nsym
     if( nsymMO( isym ) <= 0 )cycle
     CAOMO(isym)%n1 = nAO
     CAOMO(isym)%n2 = nsymMO(isym)
     allocate( CAOMO(isym)%A( nAO, nsymMO( isym ) ) )
     CAOMO(isym)%A=dmat(:,iMO+1:iMO+nsymMO(isym))
     iMO=iMO+nsymMO(isym)
  enddo
  
  !*** We must compare the order of the shells in libcint, compared
  !   with the orbitals in Dalton
  !
  !  To a first sight, it seems we should be able to extract this information
  !  from the file fort.91 and subsequently determine the permutation needed
  !  to load the orbitals. In fact, we could apply the permutation already at
  !  the level of loading or re-assigning from dmat to CAOMO the orbitals above.
  !***
  
  
  stop

  call tmpGenerateGrid( GRID_FILE )
  call LoadGrid( GRID_FILE, npt, grid )

  !*** Add the contraction of the potential on the molecular orbitals
  !*** Include the treatment of symmetry

  !.. Plot the normalization of the orbitals
  !   It is apparent that the primitive functions DO NOT get normalized
  !..
  do ish = 1, CintData.GetNshells()
     
     di = CintData.GetShellSize(ish)
     allocate(buf1e(di,di,1))

     call CintData.int1e("ovlp",ish,ish,buf1e)
     do i = 1, di
        write(*,"(2(i,x),d24.16)") ish, i, buf1e(i,i,1)
     enddo

     deallocate(buf1e)

  enddo
  

  !.. Plot the orbitals on an (x,z) grid
  !..
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

  
  !.. Load the coefficients of the molecular orbitals
  !..


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

  
  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( &
       MoleculeInputFile, &
       GTOBasisInputFile, &
       MOAO_InfoFile    , &
       LUCIA_CMOAO_NOSYM   )
    !
    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: MoleculeInputFile
    character(len=:), allocatable, intent(out) :: GTOBasisInputFile
    character(len=:), allocatable, intent(out) :: MOAO_InfoFile
    character(len=:), allocatable, intent(out) :: LUCIA_CMOAO_NOSYM
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Compute the electrostatic potential due to the product of arbitrary molecular orbitals"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help" , "Print Command Usage" )
    call List.Add( "-minp"  , "molecular input file" ,"MOLECULE.INP"      , "optional" )
    call List.Add( "-gtob"  , "GTO Basis File"       ,"cc-PVDZ-Dalton.bas", "optional" )
    call List.Add( "-moaoi" , "MO AO sym info file"  ,"fort.91"           , "optional" )
    call List.Add( "-cmoao" , "LUCIA CMOAO Nosym"    ,"LUCIA_CMOAO_NOSYM" , "optional" )

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if

    call List.Get( "-minp",  strnBuf  )
    allocate( MoleculeInputFile, source=trim(adjustl(strnBuf)))

    call List.Get( "-gtob",  strnBuf  )
    allocate( GTOBasisInputFile, source=trim(adjustl(strnBuf)))

    call List.Get( "-moaoi",  strnBuf  )
    allocate( MOAO_InfoFile, source=trim(adjustl(strnBuf)))

    call List.Get( "-cmoao",  strnBuf  )
    allocate( LUCIA_CMOAO_NOSYM, source=trim(adjustl(strnBuf)))


    call List.Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
    write(OUTPUT_UNIT,"(a)"   ) "Molecular File : "//MoleculeInputFile
    write(OUTPUT_UNIT,"(a)"   ) "GTO Basis File : "//GTOBasisInputFile
    write(OUTPUT_UNIT,"(a)"   ) "MOAO Info File : "//MOAO_InfoFile
    write(OUTPUT_UNIT,"(a)"   ) "LUCIA MO  File : "//LUCIA_CMOAO_NOSYM
    !
  end subroutine GetRunTimeParameters

  
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


end program TabulatePotential
