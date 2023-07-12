module ModuleLibCint

  
  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none

  private
  
  !.. External LibCint functions
  !..
  double precision, external :: CINTgto_norm
  integer         , external :: CINTcgto_spheric
  integer         , external :: cint1e_nuc_sph

  !.. Static parameters
  !..
  integer, parameter :: PTR_ENV_START    = 20
  integer, parameter :: MAX_N_SHELLS     = 100
  integer, parameter :: MAX_N_PRIMITIVES = 50
  !
  integer, parameter :: ATM_SLOTS  = 6
  integer, parameter :: CHARGE_OF  = 1
  integer, parameter :: PTR_COORD  = 2
  !
  integer, parameter :: BAS_SLOTS  = 8
  integer, parameter :: ATOM_OF    = 1
  integer, parameter :: ANG_OF     = 2
  integer, parameter :: NPRIM_OF   = 3
  integer, parameter :: NCTR_OF    = 4
  integer, parameter :: KAPPA_OF   = 5
  integer, parameter :: PTR_EXP    = 6
  integer, parameter :: PTR_COEFF  = 7


  !.. Wrapper for the data used by LibCint
  !.. 
  type, public :: ClassCintData
     !
     private
     !
     logical :: INITIALIZED = .FALSE.
     !
     !.. natm: Number of atoms
     !   nbas: total number of shells
     !   atm : some integer atomic parameters
     !   bas : some integer basis parameters
     !   env : some real basis parameters
     !   shls: temporary pointer to shell groups
     !..
     integer                      :: natm
     integer                      :: nbas
     integer                      :: nenv
     integer        , allocatable :: atm(:,:)
     integer        , allocatable :: bas(:,:)
     real(kind(1d0)), allocatable :: env(:)
     integer                      :: shls(4)
     !
   contains
     !
     procedure, public  :: int1e_nuc        => ClassCintData_int1e_nuc
     !
     procedure, public  :: SetCharge        => ClassCintData_SetCharge
     procedure, public  :: SetCoords        => ClassCintData_SetCoords
     procedure, public  :: SetShell         => ClassCintData_SetShell
     !
     procedure, public  :: GetNAtoms        => ClassCintData_GetNAtoms
     procedure, public  :: GetNShells       => ClassCintData_GetNShells
     procedure, public  :: GetShellSize     => ClassCintData_GetShellSize
     procedure, public  :: GetCharge        => ClassCintData_GetCharge
     procedure, public  :: GetCoords        => ClassCintData_GetCoords
     !
     procedure, public  :: Init             => ClassCintData_Init
     procedure, public  :: Free             => ClassCintData_Free
     procedure, public  :: IsInitialized    => ClassCintData_IsInitialized
     procedure, public  :: IsNotInitialized => ClassCintData_IsNotInitialized
     procedure, public  :: Print            => ClassCintData_Print
     !
     procedure, private :: ClassCintData_int1e_nuc
     procedure, private :: ClassCintData_SetCharge
     procedure, private :: ClassCintData_SetCoords
     procedure, private :: ClassCintData_SetShell
     procedure, private :: ClassCintData_GetNAtoms
     procedure, private :: ClassCintData_GetNShells
     procedure, private :: ClassCintData_GetShellSize
     procedure, private :: ClassCintData_GetCharge
     procedure, private :: ClassCintData_GetCoords
     procedure, private :: ClassCintData_Init 
     procedure, private :: ClassCintData_Free
     procedure, private :: ClassCintData_IsInitialized
     procedure, private :: ClassCintData_IsNotInitialized
     procedure, private :: ClassCintData_Print
     !
  end type ClassCintData

  
contains

  
  !=========== INTERFACES TO LIBCINT ROUTINES ===========!
  
  !.. Interface to int1e_nuc_sph
  subroutine ClassCintData_int1e_nuc( self, iShell1, iShell2, Vnuc, info )
    
    class(ClassCintData), intent(in)    :: self
    integer             , intent(in)    :: iShell1
    integer             , intent(in)    :: iShell2
    real(kind(1d0))     , intent(out)   :: Vnuc(:,:,:)
    integer, optional   , intent(out)   :: info
    
    integer :: shls(4), info_

    if ( self.IsNotInitialized()     ) return
    if ( iShell1 < 0                 ) return
    if ( iShell1 > self.GetNShells() ) return
    if ( iShell2 < 0                 ) return
    if ( iShell2 > self.GetNShells() ) return

    shls(1) = ishell1 - 1
    shls(2) = ishell2 - 1
    shls(3) = 0
    shls(4) = 0
    
    info_ = cint1e_nuc_sph( Vnuc, shls, self.atm, self.natm, self.bas, self.nbas, self.env )
    if(present(info)) info = info_
    
  end subroutine ClassCintData_int1e_nuc
  

  !============ SETUP ROUTINES ============!
  
  subroutine ClassCintData_SetShell( self, iAtom, iShell, lorb, nExp, vExp, nCtr, mCoe )

    class(ClassCintData), intent(inout) :: self
    integer             , intent(in)    :: iAtom    !.. atom on which the shell is centered
    integer             , intent(in)    :: iShell   !.. shell number
    integer             , intent(in)    :: lorb     !.. shell orbital angular momentum
    integer             , intent(in)    :: nExp     !.. Number of Exponents
    real(kind(1d0))     , intent(in)    :: vExp(:)  !.. Vector of Exponents
    integer             , intent(in)    :: nCtr     !.. Number of contractions
    real(kind(1d0))     , intent(in)    :: mCoe(:,:)!.. Matrix of Contraction coefficients

    integer :: offset
    integer :: iExp, iCtr, iBuf
    real(kind(1d0)) :: normCoef( nExp )
    
    !.. Check input consistency
    !*** Should issue a warning ***
    !..
    if ( self.IsNotInitialized()    ) return
    if ( iAtom  < 0                 ) return
    if ( iAtom  > self.GetNAtoms()  ) return
    if ( iShell < 0                 ) return
    if ( iShell > self.GetNShells() ) return
    if ( lorb < 0                   ) return
    if ( nExp < 0                   ) return
    if ( nExp > MAX_N_PRIMITIVES    ) return
    if ( nCtr > nExp                ) return
    if ( size( vExp, 1 ) < nExp     ) return
    if ( size( mCoe, 1 ) < nExp     ) return
    if ( size( mCoe, 2 ) < nCtr     ) return
    
    offset = PTR_ENV_START + 3 * self.GetNAtoms() + &
         ( iAtom - 1 ) * MAX_N_SHELLS * MAX_N_PRIMITIVES * ( MAX_N_PRIMITIVES + 1 ) + &
         ( iShell- 1 ) * MAX_N_PRIMITIVES * ( MAX_N_PRIMITIVES + 1 )

    self.bas( ATOM_OF  , iShell ) = iAtom - 1 !.. index is zero-based in this case
    self.bas( ANG_OF   , iShell ) = lorb
    self.bas( NPRIM_OF , iShell ) = nExp
    self.bas( NCTR_OF  , iShell ) = nCtr
    self.bas( PTR_EXP  , iShell ) = offset
    self.bas( PTR_COEFF, iShell ) = offset + nExp

    !.. Initialize the exponents
    offset = self.bas( PTR_EXP, iShell )
    do iExp = 1, nExp 
       self.env( offset + iExp ) = vExp( iExp )
    enddo

    !.. Initialize the normalization coefficients
    do iExp = 1, nExp
       normCoef( iExp ) = CINTgto_norm( lorb, vExp( iExp ) )
    enddo
    
    !.. Initialize the coefficients
    offset = self.bas( PTR_COEFF, iShell )
    iBuf = 0
    do iCtr = 1, nCtr
       do iExp = 1, nExp
          iBuf = iBuf + 1
          self.env( offset + iBuf ) = mCoe( iExp, iCtr ) * normCoef( iExp )
       enddo
    enddo
    
    return
    !
  end subroutine ClassCintData_SetShell


  subroutine ClassCintData_SetCoords( self, iAtom, x, y, z )

    class(ClassCintData), intent(inout) :: self
    integer             , intent(in)    :: iAtom
    real(kind(1d0))     , intent(in)    :: x,y,z

    integer :: offset
    
    if ( self.IsNotInitialized()   ) return
    if ( iAtom  < 0                ) return
    if ( iAtom  > self.GetNAtoms() ) return

    offset = PTR_ENV_START + ( iAtom - 1 ) * 3
    self.atm( PTR_COORD, iAtom ) = offset 

    self.env( offset + 1 ) = x
    self.env( offset + 2 ) = y
    self.env( offset + 3 ) = z
    
    return
    !
  end subroutine ClassCintData_SetCoords


  subroutine ClassCintData_SetCharge( self, iAtom, iCharge )

    class(ClassCintData), intent(inout) :: self
    integer             , intent(in)    :: iAtom
    integer             , intent(in)    :: iCharge
    
    if ( self.IsNotInitialized()   ) return
    if ( iAtom  < 0                ) return
    if ( iAtom  > self.GetNAtoms() ) return

    self.atm( CHARGE_OF, iAtom ) = iCharge

    return
    !
  end subroutine ClassCintData_SetCharge

  
  !=========== ACCESSOR ROUTINES ============!

  integer function ClassCintData_GetNAtoms( self ) result( nAtoms )
    !
    class(ClassCintData), intent(in) :: self
    !
    nAtoms = 0
    if( self.IsNotInitialized() )return
    nAtoms = self.natm
    !
  end function ClassCintData_GetNAtoms


  integer function ClassCintData_GetNShells( self ) result( nShells )
    !
    class(ClassCintData), intent(in) :: self
    !
    nShells = 0
    if( self.IsNotInitialized() )return
    nShells = self.nbas
    !
  end function ClassCintData_GetNShells
  

  integer function ClassCintData_GetShellSize( self, iShell ) result( iSize )
    !
    class(ClassCintData), intent(in) :: self
    integer             , intent(in) :: iShell
    !
    iSize = 0
    if( self.IsNotInitialized()    ) return
    if( iShell < 0                 ) return
    if( iShell > self.GetNShells() ) return
    !
    isize = CINTcgto_spheric( iShell - 1, self.bas )
    !
  end function ClassCintData_GetShellSize
  

  integer function ClassCintData_GetCharge( self, iAtom ) result( iCharge )
    !
    class(ClassCintData), intent(in) :: self
    integer             , intent(in) :: iAtom    
    !
    iCharge = 0
    !
    if( self.IsNotInitialized() )return
    if ( iAtom  < 0                ) return
    if ( iAtom  > self.GetNAtoms() ) return
    !
    iCharge = self.atm( CHARGE_OF, iAtom )
    !
  end function ClassCintData_GetCharge
  

  subroutine ClassCintData_GetCoords( self, iAtom, x, y, z ) 
    !
    class(ClassCintData), intent(in) :: self
    integer             , intent(in) :: iAtom
    real(kind(1d0))     , intent(out):: x,y,z
    !
    integer :: offset
    !
    if ( self.IsNotInitialized()   ) return
    if ( iAtom  < 0                ) return
    if ( iAtom  > self.GetNAtoms() ) return
    !
    offset = self.atm( PTR_COORD, iAtom )
    !
    x = self.env( offset + 1 )
    y = self.env( offset + 2 )
    z = self.env( offset + 3 )
    !
  end subroutine ClassCintData_GetCoords
  

  subroutine ClassCintData_Print( self, uid )

    class(ClassCintData), intent(in) :: self
    integer, optional   , intent(in) :: uid

    integer :: iAtom, iShell, iExp, iCtr, iCharge
    real(kind(1d0)) :: x, y, z

    integer :: offset
    
    if ( self.IsNotInitialized()   ) return
    !
    !.. Must inquiry the unit format
    write(uid,*) " N ATOMS  = ",self.GetNAtoms() 
    write(uid,*) " N SHELLS = ",self.GetNShells()
    !
    do iAtom = 1, self.GetNAtoms()

       write(uid,"(a)") "---------------------"
       write(uid,"(a,x,i0)") " ATOM #",iAtom
       write(uid,"(a,x,i0)") " CHARGE =",self.GetCharge( iAtom )
       call self.GetCoords( iAtom, x, y, z )
       write(uid,"(a,*(x,e14.6))") " COORDS =",x,y,z
       
    enddo
    !
    do iShell = 1, self.GetNShells()
       write(uid,"(a)") "---------------------"
       write(uid,"(a,x,i0)") " SHELL #", iShell
       write(uid,"(a,x,i0)") " SIZE  =", self.GetShellSize(iShell)
       write(uid,"(a,x,i0)") " ATOM  =", self.bas( ATOM_OF  , iShell )
       write(uid,"(a,x,i0)") " L     =", self.bas( ANG_OF   , iShell )
       write(uid,"(a,x,i0)") " NEXP  =", self.bas( NPRIM_OF , iShell )
       write(uid,"(a,x,i0)") " NCTR  =", self.bas( NCTR_OF  , iShell )
       write(uid,"(a,x,i0)") " EOFF  =", self.bas( PTR_EXP  , iShell )
       write(uid,"(a,x,i0)") " COFF  =", self.bas( PTR_COEFF, iShell )

       offset = self.bas( PTR_EXP, iShell )
       write(uid,"(a,*(x,e14.6))") " EXP  =", self.env( offset + 1 : offset +  self.bas( NPRIM_OF , iShell ) )

       offset = self.bas( PTR_COEFF, iShell )
       do iCtr = 1,  self.bas( NCTR_OF  , iShell )
          write(uid,"(a,i0,*(x,e14.6))") " COEF ", iCtr, &
               self.env( offset + 1 : offset + self.bas( NPRIM_OF , iShell ) )
          offset = offset + self.bas( NPRIM_OF , iShell )
       enddo
                
    enddo
    !
  end subroutine ClassCintData_Print


  logical function ClassCintData_IsNotInitialized( self ) result( IsNotInitialized )
    !
    class(ClassCintData), intent(in) :: self
    !
    IsNotInitialized = .not. ( self.IsInitialized() )
    !
  end function ClassCintData_IsNotInitialized


  logical function ClassCintData_IsInitialized( self ) result( IsInitialized )
    !
    class(ClassCintData), intent(in) :: self
    !
    IsInitialized = self.INITIALIZED
    !
  end function ClassCintData_IsInitialized

  
  subroutine ClassCintData_Free( self )
    !
    class(ClassCintData), intent(inout) :: self
    !
    self.INITIALIZED = .FALSE.
    !
    if(allocated(self.atm)) deallocate(self.atm)
    if(allocated(self.bas)) deallocate(self.bas)
    if(allocated(self.env)) deallocate(self.env)
    !
    self.natm=0
    self.nbas=0
    self.nenv=0
    self.shls=0
    !
  end subroutine ClassCintData_Free

  
  subroutine ClassCintData_Init( self, natoms, nshells )
      
    class(ClassCintData), intent(inout) :: self
    integer             , intent(in)    :: natoms
    integer             , intent(in)    :: nshells
    
    call self.free()
    !
    self.natm = natoms
    self.nbas = nshells
    self.nenv = &
         PTR_ENV_START + &
         NATOMS * 3    + &  ! atomic coordinates
         NATOMS * MAX_N_SHELLS * MAX_N_PRIMITIVES + & ! exponents                    
         NATOMS * MAX_N_SHELLS * MAX_N_PRIMITIVES**2  ! contraction coefficients
    !
    allocate( self.atm( ATM_SLOTS, natoms  ) ); self.atm = 0
    allocate( self.bas( BAS_SLOTS, nshells ) ); self.bas = 0
    allocate( self.env( self.nenv          ) ); self.env = 0.d0
    !
    self.INITIALIZED = .TRUE.
    !
    return
    !
  end subroutine ClassCintData_Init


end module ModuleLibCint
