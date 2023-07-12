
program TabulatePotential

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleLibCint
  use ModuleQCBasis
  
  implicit none
  
  !--------- Molecular data ------------------
  integer         :: NATOMS  = 3
  integer         :: NSHELLS = 4

  real(kind(1d0)) :: vExpS(1:3)     = [ 6.d0, 2.d0, 0.8d0 ]
  real(kind(1d0)) :: vCoeS(1:3,1:2) = RESHAPE( [ 0.7d0, 0.6d0, 0.5d0, 0.4d0, 0.3d0, 0.2d0 ], [3,2] ) 
  real(kind(1d0)) :: vExpP(1:1)     = [ 0.9d0 ]
  real(kind(1d0)) :: vCoeP(1:1,1:1) = RESHAPE( [ 1.d0 ], [1,1] ) 
  !-------------------------------------------
  
  
  !---------- Grid data ----------------------
  integer                     :: iy
  integer        , parameter  :: nypt = 101
  real(kind(1d0)), parameter  :: ymin = -2.d0
  real(kind(1d0)), parameter  :: ymax =  2.d0
  !
  integer                     :: iz
  integer        , parameter  :: nzpt = 101
  real(kind(1d0)), parameter  :: zmin = -2.d0
  real(kind(1d0)), parameter  :: zmax =  2.d0
  !
  real(kind(1d0))             :: x, y, z
  !-------------------------------------------
  
  !.. Local parameters
  integer                      :: uid, ish1, ish2, i, j, di, dj
  real(kind(1d0)), allocatable :: buf1e(:,:,:)

  !.. All-important molecular data type
  type(ClassCintData) :: CintData

  !.. test class shell
  type(ClassQCShell)     :: shell
  type(ClassQCAtomBasis) :: atom
  type(ClassQCAtomicSetBasis) :: atomicSet
  
  !*** MUST ADD THE PARSING OF THE WHOLE DALTON BASIS
  !    AND SUBSEQUENTLY THE DEFINITION OF THE MOLECULAR
  !    GEOMETRY, OF THE MOLECULAR ORBITALS, AND FINALLY
  !    THE CALCULATION OF THE POTENTIAL.
  !***


  character(len=*), parameter :: ATOM_HEADER = &
       "! HYDROGEN       (5s,2p,1d) -> [3s,2p,1d]"
  character(len=*), parameter :: S_SHELL_TEXT = &
       "! s functions "//QC_NEW_LINE// &
       "H    5    3   "//QC_NEW_LINE// &
       "33.8700000  0.0060680  0.0000000  0.0000000 "//QC_NEW_LINE// &
       " 5.0950000  0.0453080  0.0000000  0.0000000 "//QC_NEW_LINE// &
       " 1.1590000  0.2028220  0.0000000  0.0000000 "//QC_NEW_LINE// &
       " 0.3258000  0.0000000  1.0000000  0.0000000 "//QC_NEW_LINE// &
       " 0.1027000  0.0000000  0.0000000  1.0000000 "
  character(len=*), parameter :: P_SHELL_TEXT = &
       "! p functions "//QC_NEW_LINE// &
       "H    2    2   "//QC_NEW_LINE// &
       " 1.4070000  1.0000000  0.0000000 "//QC_NEW_LINE// &
       " 0.3880000  0.0000000  1.0000000 "
  character(len=*), parameter :: D_SHELL_TEXT = &
       "! d functions "//QC_NEW_LINE// &
       "H    1    1   "//QC_NEW_LINE// &
       " 1.0570000  1.0000000 "

  character(len=*), parameter :: ATOM_TEXT = &
       ATOM_HEADER  // QC_NEW_LINE // &
       ATOM_HEADER  // QC_NEW_LINE // &
       S_SHELL_TEXT // QC_NEW_LINE // &
       P_SHELL_TEXT // QC_NEW_LINE // &
       D_SHELL_TEXT 

  character(len=:), allocatable :: text
  
  allocate( text, source = P_SHELL_TEXT )
  !call shell.init(0,3,2,vExpS,vCoeS)

  call shell.parse(text)
  call shell.show()

  deallocate(text)
  allocate( text, source = ATOM_TEXT ) 
  call atom.parse(text)
  call atom.show()

  write(*,"(a)") ATOM_TEXT
  call atomicSet.parseFile( "cc-PVTZ-Dalton.bas" )
  call atomicSet.show()

  stop

  
  !.. Initializes the molecular basis
  !..
  call CintData.init( NATOMS, NSHELLS )
  !
  call CintData.setCharge( 1, 0 )                        ! ( Atom, Charge )
  call CintData.setCharge( 2, 0 )
  call CintData.setCharge( 3, 1 )                        ! <- atom #3 is a test charge
  !
  call CintData.setCoords( 1, 0.d0, 0.d0, -0.8d0 )       ! ( Atom, x, y, z )
  call CintData.setCoords( 2, 0.d0, 0.d0,  0.8d0 )
  call CintData.setCoords( 3, 0.d0, 0.d0,  0.0d0 )
  !
  call CintData.setShell ( 1, 1, 0, 3, vExpS, 2, vCoeS ) ! ( Atom, sh, l, nPrim, vec Exp, nContr, mat Coef )
  call CintData.setShell ( 1, 2, 1, 1, vExpP, 1, vCoeP )
  call CintData.setShell ( 2, 3, 0, 3, vExpS, 2, vCoeS )
  call CintData.setShell ( 2, 4, 1, 1, vExpP, 1, vCoeP )
  !
  call CintData.Print( OUTPUT_UNIT )

  
  !.. Evaluate the potential of the product of contracted
  !   functions in shells 1 and 4.
  !..
  ish1 = 1
  ish2 = 4
  di = CintData.GetShellSize( ish1 )
  dj = CintData.GetShellSize( ish2 )
  allocate(buf1e(di,dj,1))
  !
  open(newunit = uid       , &
       file    ="Pot"      , &
       form    ="formatted", &
       status  ="unknown"  )
  !
  x = 0.3d0
  do iy = 1, nypt
     y = ymin + (ymax - ymin) * dble(iy-1)/dble(nypt-1)
     do iz = 1, nzpt
        z = zmin + (zmax - zmin) * dble(iz-1)/dble(nzpt-1)
        !
        call CintData.SetCoords( 3, x, y, z )
        call CintData.int1e_nuc( ish1, ish2, buf1e )
        !
        write(uid,"(*(x,e14.6))") y, z, ( ( buf1e(i,j,1), i = 1, di ), j = 1, dj )
     enddo
     write(uid,*) 
  enddo
  deallocate(buf1e)
  close(uid)

  
end program TabulatePotential
