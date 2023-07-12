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
Module ModuleDensityMatrices

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleSystemUtils
  use ModuleAngularMomentum
  use ModuleString
  use ModuleGroups
  use ModuleMatrix
  use ModuleParentIons
  use ModuleIntegrals

  implicit none

  private

  integer         , parameter :: N_REDUCED_TDM1B_TYPES = 2
  integer         , parameter :: N_REDUCED_TDM2B_TYPES = 6
  character(len=1), parameter :: REDUCED_TDM1B_TYPES(*)=["0","1"]
  character(len=3), parameter :: REDUCED_TDM2B_TYPES(*)=["000","011","101","110","111","112"]
  character(len=*), parameter :: TDM_GLOBAL_DIR        ="TDM/"

  character(len=*), parameter :: TDM1_ROOT = "TDM1_"
  character(len=*), parameter :: TDM2_ROOT = "TDM2_"
  character(len=*), parameter :: HAIC_ROOT = "HAIC_"
  character(len=*), parameter :: PIEN_ROOT = "PIEN_"

  type, public :: ClassTDMGeneral
     logical :: Initialized = .FALSE.
   contains
     procedure, public  :: IsInitialized => ClassTDMGeneralIsInitialized
     procedure, public  :: Write => ClassTDMGeneralDoNothing1arg
     procedure, public  :: Read  => ClassTDMGeneralDoNothing1arg
     procedure, public  :: Free  => ClassTDMGeneralDoNothing0arg
  end type ClassTDMGeneral

  type, private, extends(ClassTDMGeneral) :: ClassScalar
     private
     real(kind(1d0)) :: Energy
   contains
     procedure, public  :: Free      => ClassScalarFree
     procedure, public  :: Init      => ClassScalarInit
     procedure, public  :: Set       => ClassScalarSet
     procedure, public  :: Write     => ClassScalarWriteToUnit
     procedure, public  :: Read      => ClassScalarReadFromUnit
  end type ClassScalar
  
  !> One-body transition density matrix
  !! The reduced 1BTDM are defined as
  !! \[
  !!    \mathsf{R}^{BA}_{[q,p]_T}=\langle A \| [b^\dagger_p\otimes b_q]_T\|B\rangle
  !! \]
  !! \[
  !!    \mathsf{R}^{BA}_{[q,p]_{T}}C_{S_B \Sigma_B,T\tau}^{S_A\Sigma_A}=
  !!      \Pi_{S_A}\sum_{\pi\theta}C_{\frac{1}{2}\pi,\frac{1}{2}-\theta}^{T\tau} (-1)^{-\frac{1}{2}+\theta}\rho^{BA}_{QP}
  !! Let's consider $T=0,1$ separately. For $T=0$, the matrix $\mathsf{R}^{BA}_{[q,p]_{0}}$ is non zero only if
  !! the two ions have the same multiplicity. For any value of their projection $\Sigma_A=\Sigma_B$, we can compute
  !! the reduced matrix element as
  !! \[
  !!    \mathsf{R}^{BA}_{[q,p]_{0}}&=\sqrt{S_A+1/2}\,\sum_{\sigma}\rho^{BA}_{q_\sigma p_\sigma}.
  !! \]
  !! For $T=1$, the reduced matrix element is zero if both $A$ and $B$ are singlet states. For parent ions with an
  !! odd number of electrons, it is possible to compute $\mathsf{R}^{BA}_{[q,p]_{1}}$ from the uncoupled density
  !! matrix between the states with $\Sigma_A=\Sigma_B=\Sigma=1/2$,
  !! \[
  !!    \mathsf{R}^{BA}_{[q,p]_{1}} =   \frac{\sqrt{S_A+1/2}}{C_{S_B 1/2,10}^{S_A1/2}} ( \rho^{BA}_{q_\alpha p_\alpha}
  !!                                  - \rho^{BA}_{q_\beta p_\beta}).
  !! \]
  !! For parent ions with an even number of electrons, the matrix element between ions with different spin
  !! can be computed choosing $\Sigma_A=\Sigma_B=\Sigma=0$ (which is required if one of the two spin is zero).
  !! For parent ions with the same, non-zero spin, however, it is necessary to use any other spin projection
  !! (e.g., $\Sigma_A=\Sigma_B=\Sigma=1$)
  !! \[
  !!    \mathsf{R}^{BA}_{[q,p]_{1}} &= \frac{\sqrt{S_A+1/2}}{C_{S_B \Sigma,10}^{S_A\Sigma}}
  !!    ( \rho^{BA}_{q_\alpha p_\alpha} - \rho^{BA}_{q_\beta p_\beta}).
  !! \]
  !!
  !! \[
  !!    A^{BA}_{sr} = \sqrt{\frac{2}{2S_A+1}}\,\mathsf{R}^{BA}_{[s,r]_0}, \qquad (S_B=S_A),
  !! \]
  !! \[
  !!    B^{BA}_{pr}=-(-1)^{S_B+1/2+S}\sum_T \Pi_{T}\sjs{T}{1/2}{1/2}{S}{S_B}{S_A}\mathsf{R}^{BA}_{[p,r]_T}.
  !! \]
  !<
  type, public, extends(ClassTDMGeneral) :: ClassTDM1

     private

     !.. Irreducible representation of the AB product
     type(ClassIrrep), pointer :: IrrepAB
     !.. SA2 = 2 * S_A, SB2 = 2 * S_B
     integer :: SA2, SB2
     
     !.. Number of active orbitals per irreducible representation
     integer          , allocatable :: nactive(:)

     !.. Matrices R_[pq]0, R_[pq]1
     !   the first index is for the irrep of the first orbital. The irrep of the second
     !   is given by IrrepAB * irrep. 
     !   The second index is for the spin-coupling type, i.e., T=0 or T=1 
     !   (indexes 1 and 2, respectively)
     !   The R matrices do not depend on the total spin
     type(ClassMatrix), pointer     :: R(:,:)

     !.. Matrix A_pq (proportional to singlet spin coupling)
     !   the index is for the irrep of the first orbital. The irrep of the second
     !   is given by IrrepAB * irrep
     logical          , allocatable :: A_IS_COMPUTED(:)
     type(ClassMatrix), pointer     :: A(:)

     !.. Matrix B_pq 
     !   the index is for the irrep of the first orbital. The irrep of the second
     !   is given by IrrepAB * irrep
     !   The B matrix *does* depend on the total spin coupling,
     !   it is computed only upon request, and it is not recomputed if the
     !   requested spin coupling is already available 
     !.. S2 = 2 * S, where S is the total spin of the augmented state "ion+electron"
     integer                        :: S2
     logical                        :: B_IS_COMPUTED
     type(ClassMatrix), pointer     :: B(:)

   contains

     procedure, public  :: Free        => ClassTDM1Free
     procedure, public  :: Init        => ClassTDM1Init
     procedure, public  :: SetBlock    => ClassTDM1SetBlock
     procedure, public  :: Write       => ClassTDM1WriteToUnit
     procedure, public  :: Read        => ClassTDM1ReadFromUnit
     procedure, public  :: GetABlock   => ClassTDM1GetABlock
     procedure, public  :: GetBBlock   => ClassTDM1GetBBlock
     procedure, public  :: GetABlockf  => ClassTDM1GetABlockFun
     procedure, public  :: GetBBlockf  => ClassTDM1GetBBlockFun
     procedure, public  :: GetRTBlock  => ClassTDM1GetRTBlock
     procedure, public  :: GetRTBlockf => ClassTDM1GetRTBlockFun
     procedure, private :: ComputeA    => ClassTDM1ComputeA
     procedure, private :: ComputeB    => ClassTDM1ComputeB

     final :: ClassTDM1Final

  end type ClassTDM1

  !> Two-body transition density matrix
  !! The reduced 2BTDM are defined as
  !! \[
  !!    \mathsf{\Pi}^{BA}_{[[sr]_J,[pq]_T]_K} =\frac{-(-1)^{\sigma+\rho} \Pi_{S_A}}{C_{S_B\Sigma_B,K\kappa}^{S_A\Sigma_A}}
  !!                                           \sum_{\tau,\mu}C_{T\tau,J-\mu}^{K\kappa}
  !!                                           \sum_{\sigma,\rho,\pi,\theta}\pi^{BA}_{RS,PQ}  \times
  !!                                           C_{\frac{1}{2}-\sigma,\frac{1}{2}-\rho}^{J-\mu}C_{\frac{1}{2}\pi,\frac{1}{2}\theta}^{T\tau}
  !! \]
  !!     P^{BA}_{ps,qr} = -(-1)^{S+S_B+1/2}\sum_{KTJ}(-1)^{J+K}\Pi_{KTJ}\sjs{S_A}{K}{S_B}{1/2}{S}{1/2} \sjs{J}{T}{K}{1/2}{1/2}{1/2} \mathsf{\Pi}^{BA}_{[[ps]_J,[qr]_T]_K}
  !! \[
  !!
  !! \]
  !!     Q^{BA}_{ps,qr} = -(-1)^{S+S_B+1/2}\sum_{KTJ}(-1)^{T+K}\Pi_{KTJ}\sjs{S_A}{K}{S_B}{1/2}{S}{1/2} \sjs{J}{T}{K}{1/2}{1/2}{1/2} \mathsf{\Pi}^{BA}_{[[ps]_J,[qr]_T]_K}
  !! \[
  !!
  !<
  type, public, extends(ClassTDMGeneral) :: ClassTDM2

     private

     !.. Irreducible representation of the AB product
     type(ClassIrrep) :: IrrepAB
     !.. SA2 = 2 * S_A, SB2 = 2 * S_B
     integer :: SA2, SB2

     !.. Matrices Pi_{[[sr]_J,[pq]_T]_K} 
     !   where JTK = 000, 110, 101, 011, 111, 112. 
     !   The first three indexes are for the irrep of 
     !   the first three orbitals. The irrep of the fourth one
     !   is given by Irrep_q = Irrep_s * Irrep_r * Irrep_p. 
     !   The last is for the spin-coupling type, i.e., 
     !   for the index pointing to the values of the JTK triad 
     !   as described by the vector REDUCED_TDM2B_TYPES(*)
     !   The Pi matrices do not depend on the total spin.
     !..
     type(ClassMatrix4D), pointer :: Pi(:,:,:,:)

     !.. Number of active orbitals per irreducible representation
     integer          , allocatable :: nactive(:)

     !.. Matrices P_pq and  Q_pq
     !   The irrep of the fourth orbital 
     !   is given by Irrep_q = IrrepAB * Irrep_s * Irrep_r * Irrep_p
     integer                        :: S2
     !if SA == SB  =/= 0 
     !   S = SA-1/2   S_<
     !   S = SA+1/2   S_>
     !if SA == SB == 0
     !   S = 1/2      
     !if SA == SB +/- 1
     !   S = SA +/- 1/2
     !
     logical          , allocatable :: Q_IS_COMPUTED(:,:,:,:)
     type(ClassMatrix4D), pointer   :: Q(:,:,:,:)!<=== must have index 2 and 3 flipped 

   contains

     procedure, public  :: Free      => ClassTDM2Free
     generic,   public  :: Init      => ClassTDM2Init
     generic,   public  :: SetBlock  => ClassTDM2SetBlock
     procedure, public  :: Write     => ClassTDM2WriteToUnit
     procedure, public  :: Read      => ClassTDM2ReadFromUnit
     generic,   public  :: GetPiBlock  => ClassTDM2GetPiBlock

!!     procedure, public  :: GetQBlock   => ClassTDM2GetQBlock
     procedure, public  :: GetQBlockf  => ClassTDM2GetQBlockFun
     procedure, private :: ComputeQ    => ClassTDM2ComputeQ
     

     procedure, private :: ClassTDM2Free
     procedure, private :: ClassTDM2Init 
     procedure, private :: ClassTDM2SetBlock
     procedure, private :: ClassTDM2WriteToUnit
     procedure, private :: ClassTDM2ReadFromUnit
     procedure, private :: ClassTDM2GetPiBlock
!!     procedure, private :: ClassTDM2GetQBlock
     procedure, private :: ClassTDM2GetQBlockFun
     procedure, private :: ClassTDM2ComputeQ
     

     final :: ClassTDM2Final

  end type ClassTDM2

  type, private :: ClassTDMManager
     private
     integer                           :: Nions = 0
     integer             , allocatable :: ninactive(:)
     integer             , allocatable :: nactive(:)
     character(len=:)    , allocatable :: StorageDir 
     type(ClassParentIon), pointer     :: Ionv(:)
     type(ClassScalar)   , pointer     :: PIEN(:)
     !.. The first (second) index correspond to the location in the Ionv of the ion in the bra (ket)
     type(ClassTDM1)     , pointer     :: HAIC(:,:)
     type(ClassTDM1)     , pointer     :: TDM1(:,:)
     type(ClassTDM2)     , pointer     :: TDM2(:,:)
!!$     type(ClassTDM3)     , pointer :: TDM3(:,:)

     type(ClassMatrix), allocatable :: Multipoles(:,:)

   contains

     generic,   public  :: Init    => ClassTDMManager_Init
     generic,   public  :: SetParentIons => ClassTDMManager_SetParentIons
     procedure, public  :: SetStorageDir => ClassTDMManager_SetStorageDir
     procedure, public  :: SetNActiveEls => ClassTDMManager_Setnactive
     procedure, public  :: GetNActiveEls => ClassTDMManager_Getnactive
     generic,   public  :: Free    => ClassTDMManager_Free
     procedure, public  :: SetTDM1 => ClassTDMManager_SetTDM1
     procedure, public  :: SetTDM2 => ClassTDMManager_SetTDM2
     procedure, public  :: SetHAIC => ClassTDMManager_SetHAIC
     procedure, public  :: SetPIEN => ClassTDMManager_SetPIEN
     procedure, public  :: IO      => ClassTDMManager_IO
     procedure, public  :: GetHAIC => ClassTDMManager_GetHAIC

     procedure, private :: GetIonIndex => ClassTDMManager_GetIonIndex
     procedure, private :: ClassTDMManager_SetParentIons
     procedure, private :: ClassTDMManager_Init
     procedure, private :: ClassTDMManager_Free

!!$     procedure, public  :: ClassTDMManager_ComputeMoments
     procedure, public  :: GetMatEl => ClassTDMManager_GetMatEl  !<--- <A|X|B> for a required lm, and A&B labels
!!$     procedure, public  :: ClassTDMManager_GetDipole !<--- <A|X|B> for a required lm, and A&B labels
!!$     procedure, public  :: ClassTDMManager_GetMoment !<--- <A|r^lXlm|B> for a required lm, and A&B labels
     procedure, public :: GetB => ClassTDMManager_GetB !<-function that returns pointer to matrix B between two ions
     procedure, public :: GetA => ClassTDMManager_GetA !<-function that returns pointer to matrix A between two ions
     procedure, public :: GetQ => ClassTDMManager_GetQ
     
     final :: ClassTDMManager_Final

  end type ClassTDMManager

  type(ClassTDMManager), public :: TDM_Manager

contains

  function ClassTDMManager_GetQ( self, IonBraLabel, IonKetLabel, iIrrP, iIrrQ, iIrrS, S2 ) result(Q)
    class(ClassTDMManager), intent(inout) :: self
    character(len=*)      , intent(in)    :: IonBraLabel, IonKetLabel
    integer               , intent(in)    :: iIrrP, iIrrQ, iIrrS
    integer               , intent(in)    :: S2
    type(ClassMatrix4D), pointer :: Q
    integer                    :: iIonBra, iIonKet
    iIonBra    = self%GetIonIndex(IonBraLabel)
    iIonKet    = self%GetIonIndex(IonKetLabel)
    Q => self%TDM2(iIonBra,iIonKet)%GetQBlockF(iIrrP, iIrrQ, iIrrS, S2)
  end function ClassTDMManager_GetQ

  !.. <A|O|B> = \Pi_{S_A}^{-1}   <A||O||B>
  real(kind(1d0)) function ClassTDMManager_GetMatEl( self, OpLabel, IonBraLabel, IonKetLabel, ilm ) result(O_AB)
    class(ClassTDMManager), intent(inout) :: self
    character(len=*)      , intent(in)    :: OpLabel, IonBraLabel, IonKetLabel
    integer, optional     , intent(in)    :: ilm
    integer                    :: iIonBra, iIonKet, iIrrOrbBra, iIrrOrbKet, nIrreps, SA2, SB2
    type(ClassMatrix)          :: OneB
    type(ClassMatrix), pointer :: A
    type(ClassIrrep) , pointer :: IrrIonBra, IrrIonKet, IrrOrbBra
    type(ClassIrrep) , pointer :: irreplist(:)

    irreplist => GlobalGroup%GetIrrepList()
    nIrreps   =  GlobalGroup%GetNIrreps()
    iIonBra   = self%GetIonIndex(IonBraLabel)
    iIonKet   = self%GetIonIndex(IonKetLabel)
    IrrIonBra => self%Ionv(iIonBra)%GetIrrep()
    IrrIonKet => self%Ionv(iIonKet)%GetIrrep()
    O_AB = 0.d0
    !..  SB must be equal to SA otherwhise A is 0
    SA2=self%Ionv(iIonBra)%GetMultiplicity()-1
    SB2=self%Ionv(iIonKet)%GetMultiplicity()-1
    if(SA2.ne.SB2)return
    do iIrrOrbBra = 1, nIrreps
       iIrrOrbKet = GlobalGroup%GetIrrepIndex((irreplist(iIrrOrbBra) * IrrIonBra) * IrrIonKet )
       if(present(ilm))then
          OneB = GlobalIntegral%Get1B(OpLabel, "MO_MO", iIrrOrbBra, iIrrOrbKet, ilm)
       else
          OneB = GlobalIntegral%Get1B(OpLabel, "MO_MO", iIrrOrbBra, iIrrOrbKet)
       endif
       if(OneB%NRows()*OneB%NColumns().eq.0)cycle
       call OneB%Drop("rows"   ,1,self%ninactive(iIrrOrbBra)) 
       call OneB%Drop("columns",1,self%ninactive(iIrrOrbKet))
       call OneB%Transpose()
       !**** POSSIBLY, GETA DOES NOT RETURN THE A MATRIX
       !  write(*,*) "LLLLEEEEGGGOOOO",IonBraLabel, IonKetLabel, iIrrOrbBra!,self%TDM1(iIonBra,iIonKet)%A_IS_COMPUTED(iIrrOrbBra)
       A => self%GetA( IonBraLabel, IonKetLabel, irreplist(iIrrOrbBra) )
       if(.not.self%TDM1(iIonBra,iIonKet)%A_IS_COMPUTED(iIrrOrbBra))cycle
       !write(*,*) "PASOOOOOOOOO",IonBraLabel, IonKetLabel, iIrrOrbBra,allocated(self%TDM1(iIonBra,iIonKet)%A_IS_COMPUTED)
       O_AB = O_AB + ClassMatrixFullContraction(A,OneB)
       call OneB%Free()
       nullify(A)
    enddo
  end function ClassTDMManager_GetMatEl

  function ClassTDMManager_GetA( self, IonBraLabel, IonKetLabel, IrrOrbBra ) result(A)
    class(ClassTDMManager), intent(inout) :: self
    character(len=*)      , intent(in)    :: IonBraLabel, IonKetLabel
    type(ClassIrrep)      , intent(in)    :: IrrOrbBra
    type(ClassMatrix), pointer :: A
    integer                    :: iIonBra, iIonKet, iIrrOrbBra, SA2, SB2
    iIonBra    = self%GetIonIndex(IonBraLabel)
    iIonKet    = self%GetIonIndex(IonKetLabel)
    SA2=self%Ionv(iIonBra)%GetMultiplicity()-1
    SB2=self%Ionv(iIonKet)%GetMultiplicity()-1
    if(SA2.ne.SB2)return
    iIrrOrbBra = GlobalGroup%GetIrrepIndex(IrrOrbBra)
    A => self%TDM1(iIonBra,iIonKet)%GetABlockf(iIrrOrbBra)
  end function ClassTDMManager_GetA

  function ClassTDMManager_GetB( self, IonBraLabel, IonKetLabel, IrrOrbBra ) result(B)
    class(ClassTDMManager), intent(inout) :: self
    character(len=*)      , intent(in)    :: IonBraLabel, IonKetLabel
    type(ClassIrrep)      , intent(in)    :: IrrOrbBra
    type(ClassMatrix), pointer :: B
    integer                    :: iIonBra, iIonKet, iIrrOrbBra
    iIonBra    = self%GetIonIndex(IonBraLabel)
    iIonKet    = self%GetIonIndex(IonKetLabel)
    iIrrOrbBra = GlobalGroup%GetIrrepIndex(IrrOrbBra)
    B => self%TDM1(iIonBra,iIonKet)%GetBBlockf(iIrrOrbBra)
  end function ClassTDMManager_GetB

  subroutine ClassTDMManager_SetStorageDir( self, StorageDir )
    class(ClassTDMManager), intent(inout) :: self
    character(len=*)      , intent(in)    :: StorageDir
    allocate(self%StorageDir,source=StorageDir)
    call Execute_Command_Line("mkdir -p "//FormatAsDir(self%StorageDir) // TDM_GLOBAL_DIR)
  end subroutine ClassTDMManager_SetStorageDir

  subroutine ClassTDMManager_Setnactive( self, ninactive, nactel )
    class(ClassTDMManager), intent(inout) :: self
    integer               , intent(in)    :: ninactive(:)
    integer               , intent(in)    :: nactel(:)
    allocate(self%ninactive,source=ninactive)
    allocate(self%nactive  ,source=nactel)
  end subroutine ClassTDMManager_Setnactive

  subroutine ClassTDMManager_Getnactive( self, ninactive, nactel )
    class(ClassTDMManager), intent(in)  :: self
    integer, allocatable  , intent(out) :: ninactive(:)
    integer, allocatable  , intent(out) :: nactel(:)
    allocate(ninactive,source=self%ninactive)
    allocate(nactel   ,source=self%nactive)
  end subroutine ClassTDMManager_Getnactive

  subroutine ClassTDMManager_SetParentIons( self, PionList )
    class(ClassTDMManager),          intent(inout) :: self
    type (ClassParentIon ), pointer, intent(in)    :: PionList(:)
    call self%Free()
    self%Nions=size(PionList)
!!$    allocate(self%Ionv(self%Nions))
!!$    do iIon = 1, self%Nions
!!$       self%Ionv(iIon) => PionList(iIon)
!!$    enddo
    allocate(self%Ionv,source=PionList)
  end subroutine ClassTDMManager_SetParentIons

  subroutine ClassTDMManager_SetHAIC( self, iIonA, iIonB, iIrr1, Sigma2, nr, nc, RhoUp, LDU, RhoDown, LDD )
    class(ClassTDMManager), intent(inout) :: self
    integer               , intent(in)    :: iIrr1, Sigma2, iIonA, iIonB, nr, nc, LDU, LDD
    real(kind(1d0))       , intent(in)    :: RhoUp(LDU,*), RhoDown(LDD,*)
    call self%HAIC(iIonA,iIonB)%SetBlock(iIrr1, Sigma2, nr, nc, RhoUp, LDU, RhoDown, LDD)
  end subroutine ClassTDMManager_SetHAIC

  subroutine ClassTDMManager_GetHAIC( self, IonALabel, IonBLabel, Irr1, H0, H1 )
    class(ClassTDMManager), intent(inout) :: self
    character(len=*)      , intent(in)    :: IonALabel, IonBLabel
    type(ClassIrrep)      , intent(in)    :: Irr1
    type(ClassMatrix)     , intent(out)   :: H0, H1
    integer :: iIrr1, iIonA, iIonB
    iIonA = self%GetIonIndex( IonALabel )
    iIonB = self%GetIonIndex( IonBLabel )
    iIrr1 = GlobalGroup%GetIrrepIndex(Irr1)
    call self%HAIC(iIonA,iIonB)%GetRTBlock(iIrr1, "0", H0 )
    call self%HAIC(iIonA,iIonB)%GetRTBlock(iIrr1, "1", H1 )
  end subroutine ClassTDMManager_GetHAIC

  integer function ClassTDMManager_GetIonIndex( self, IonLabel ) result(iIon)
    class(ClassTDMManager)      , intent(inout) :: self
    character(len=*)            , intent(in)    :: IonLabel
    do iIon = 1, self%Nions
       if( IonLabel .is. ( self%Ionv(iIon)%GetLabel() ) ) exit
    enddo
  end function ClassTDMManager_GetIonIndex

  subroutine ClassTDMManager_SetPIEN( self, iIon, Energy )
    class(ClassTDMManager), intent(inout) :: self
    integer               , intent(in)    :: iIon
    real(kind(1d0))       , intent(in)    :: Energy
    call self%PIEN(iIon)%Set(Energy)
  end subroutine ClassTDMManager_SetPIEN

  subroutine ClassTDMManager_SetTDM1( self, iIonA, iIonB, iIrr1, Sigma2, nr, nc, RhoUp, LDU, RhoDown, LDD )
    class(ClassTDMManager), intent(inout) :: self
    integer               , intent(in)    :: iIrr1, Sigma2, iIonA, iIonB, nr, nc, LDU, LDD
    real(kind(1d0))       , intent(in)    :: RhoUp(LDU,*), RhoDown(LDD,*)
    call self%TDM1(iIonA,iIonB)%SetBlock(iIrr1, Sigma2, nr, nc, RhoUp, LDU, RhoDown, LDD)
  end subroutine ClassTDMManager_SetTDM1

  subroutine ClassTDMManager_SetTDM2( self, iIonA, iIonB, iIrr1, iIrr2, iIrr3, Sigma2, PIAB, PIAA, PIBB )
    class(ClassTDMManager), intent(inout) :: self
    integer               , intent(in)    :: iIrr1, iIrr2, iIrr3, Sigma2, iIonA, iIonB
    type(ClassMatrix4D)   , intent(in)    :: PIAB, PIAA, PIBB
    call self%TDM2(iIonA,iIonB)%SetBlock( iIrr1, iIrr2, iIrr3, Sigma2, PIAB, PIAA, PIBB )
  end subroutine ClassTDMManager_SetTDM2

  subroutine ClassTDMManager_IO( self, ID_CHAR, Action )
    class(ClassTDMManager), intent(inout) :: self
    character             , intent(in)    :: ID_CHAR
    character(len=*)      , intent(in)    :: Action
    !
    character(len=1000) :: FileName, iomsg
    integer :: uid, iIon1, iIon2, iostat
    class(ClassTDMGeneral), pointer :: pTDM
    character(:), allocatable :: sROOT

    call Execute_Command_Line("mkdir -p "//&
         FormatAsDir(trim(self%StorageDir))//TDM_GLOBAL_DIR )

    select case ( ID_CHAR )
    case( "1" )
       allocate(sROOT,source=TDM1_ROOT)
    case( "2" )
       allocate(sROOT,source=TDM2_ROOT)
    case( "H" )
       allocate(sROOT,source=HAIC_ROOT)
    case( "E" )
       allocate(sROOT,source=PIEN_ROOT)
    end select

    do iIon1 = 1, self%Nions
       do iIon2 = 1, self%Nions

          select case ( ID_CHAR )
          case( "1" )
             pTDM => self%TDM1(iIon1,iIon2)
          case( "2" )
             pTDM => self%TDM2(iIon1,iIon2)
          case( "H" )
             pTDM => self%HAIC(iIon1,iIon2)
          case( "E" )
             if(iIon1/=iIon2)cycle
             pTDM => self%PIEN(iIon1)
          end select

          if( ( Action .is. "Write" ) .and. .not. pTDM%IsInitialized() ) cycle

          FileName = FormatAsDir(Self%StorageDir) // TDM_GLOBAL_DIR // sROOT // &
               self%Ionv(iIon1)%GetLabel() // "_" // &
               self%Ionv(iIon2)%GetLabel()
          !
          open(newunit = uid, &
               file    = trim(FileName), &
               status  ="unknown", &
               form    ="formatted",&
               action  = Action, &
               iostat  = iostat, &
               iomsg   = iomsg )
          if(iostat/=0)then
             write(*,"(a)") "Error opening "//Filename//" : "//trim(iomsg)
             stop
          endif

          if( Action .is. "Write" ) call pTDM%Write(uid)
          if( Action .is. "Read"  ) call pTDM%Read(uid)

          close(uid)

       enddo
    enddo
    
  end subroutine ClassTDMManager_IO

!!$  subroutine ClassTDMManager_SetNActiveEls( self, ninactive, nactel )
!!$    class(ClassTDMManager), intent(inout) :: self
!!$    integer               , intent(in)    :: ninactive(:) 
!!$    integer               , intent(in)    :: nactel(:) 
!!$    allocate(self%ninactive,source=ninactive)
!!$    allocate(self%nactive,source=nactel)
!!$  end subroutine ClassTDMManager_SetNActiveEls
!!$
  subroutine ClassTDMManager_Init( self )
    class(ClassTDMManager),           intent(inout) :: self
    !
    type(ClassIrrep) :: IrrepAB
    integer :: iIon1, iIon2, nIons, SA2, SB2, nirrep, uid

    if(.not.allocated(self%nactive))then
       open( newunit=uid, file=FormatAsDir(self%StorageDir) // TDM_GLOBAL_DIR // "MOset", &
            form="formatted", status="unknown", action="read")
       read(uid,*) nirrep
       allocate(self%ninactive(nirrep))
       allocate(self%nactive(nirrep))
       read(uid,*) self%ninactive
       read(uid,*) self%nactive
       close(uid)
    else
       open( newunit=uid, file=FormatAsDir(self%StorageDir) // TDM_GLOBAL_DIR // "MOset", &
            form="formatted", status="unknown", action="write")
       write(uid,"(i0)") size(self%nactive)
       write(uid,"(*(x,i0))") self%ninactive
       write(uid,"(*(x,i0))") self%nactive
       close(uid)
    endif

    nIons = size(self%ionv)

    !.. Initialize size of TDMn
    allocate(self%TDM1(nIons,nIons))
    allocate(self%TDM2(nIons,nIons))
    allocate(self%HAIC(nIons,nIons))
    allocate(self%PIEN(nIons))

    do iIon1 = 1, nIons
       do iIon2 = 1, nIons
          IrrepAB = self%ionv(iIon1)%GetIrrep() *  self%ionv(iIon2)%GetIrrep()
          SA2 = self%ionv(iIon1)%GetMultiplicity() - 1
          SB2 = self%ionv(iIon2)%GetMultiplicity() - 1
          call self%TDM1(iIon1,iIon2)%Init( IrrepAB, SA2, SB2, self%nactive )
          call self%HAIC(iIon1,iIon2)%Init( IrrepAB, SA2, SB2, self%nactive )
          call self%TDM2(iIon1,iIon2)%Init( IrrepAB, SA2, SB2, self%nactive )
       enddo
       call self%PIEN(iIon1)%Init()
    enddo

  end subroutine ClassTDMManager_Init

  subroutine ClassTDMManager_Free( self )
    class(ClassTDMManager),  intent(inout) :: self
    integer :: i, j
    !
    if(.not.associated(self%Ionv))return
    do i=1,size(self%Ionv)
       call self%Ionv(i)%Free()
    enddo
    nullify(self%Ionv)
    !
    do i=1,self%nions
       do j=1,self%nions
          call self%TDM1(i,j)%Free()
          call self%TDM2(i,j)%Free()
          call self%HAIC(i,j)%Free()
       enddo
       call self%PIEN(i)%Free()
    enddo
    deallocate(self%TDM1)
    deallocate(self%TDM2)
    deallocate(self%HAIC)
    deallocate(self%PIEN)
    self%Nions=0
    !
  end subroutine ClassTDMManager_Free

  subroutine ClassTDMManager_Final( self )
    type(ClassTDMManager),  intent(inout) :: self
    !
    call self%Free()
    !
  end subroutine ClassTDMManager_Final

  subroutine ClassScalarFree( self )
    class(ClassScalar), intent(inout) :: self
    self%Energy=0.d0
  end subroutine ClassScalarFree

  subroutine ClassScalarInit( self )
    class(ClassScalar), intent(inout) :: self
    self%Energy=0.d0
    self%Initialized=.TRUE.
  end subroutine ClassScalarInit

  subroutine ClassScalarSet( self, Energy )
    class(ClassScalar), intent(inout) :: self
    real(kind(1d0))  , intent(in)    :: Energy
    self%Energy=Energy
  end subroutine ClassScalarSet

  subroutine ClassScalarWriteToUnit( self, uid )
    class(ClassScalar), intent(inout) :: self
    integer          , intent(in)    :: uid
    character(len=30) :: Form
    logical           :: Formatted
    INQUIRE( UNIT = uid, FORM  = Form )
    Formatted = trim(Form) .is. "FORMATTED" 
    if(Formatted)then
       write(uid,"(x,e24.16)") self%Energy
    else
       write(uid) self%Energy
    end if
  end subroutine ClassScalarWriteToUnit

  subroutine ClassScalarReadFromUnit( self, uid )
    class(ClassScalar), intent(inout) :: self
    integer          , intent(in)    :: uid
    character(len=30) :: Form
    logical           :: Formatted
    INQUIRE( UNIT = uid, FORM  = Form )
    Formatted = trim(Form) .is. "FORMATTED" 
    if(Formatted)then
       read(uid,*) self%Energy
    else
       read(uid) self%Energy
    end if
  end subroutine ClassScalarReadFromUnit

  subroutine ClassTDMGeneralDoNothing0arg( self )
    Class(ClassTDMGeneral), intent(inout) :: self
  end subroutine ClassTDMGeneralDoNothing0arg

  subroutine ClassTDMGeneralDoNothing1arg( self, uid )
    Class(ClassTDMGeneral), intent(inout) :: self
    integer               , intent(in)    :: uid
  end subroutine ClassTDMGeneralDoNothing1arg

  subroutine ClassTDM1Final( self )
    type(ClassTDM1), intent(inout) :: self
    call self%Free()
  end subroutine ClassTDM1Final

  logical function ClassTDMGeneralIsInitialized( self ) result( res )
    class(ClassTDMGeneral), intent(in) :: self
    res = self%Initialized
  end function ClassTDMGeneralIsInitialized

  subroutine ClassTDM1Free( self )
    class(ClassTDM1), intent(inout) :: self
    integer :: i, iType
    if(associated(self%IrrepAB)) self%IrrepAB => NULL()
    if( allocated(self%nactive)) deallocate(self%nactive)
    if(associated(self%R))then
       do iType = 1, size(self%R,2)
          do i = 1, size(self%R,1)
             call self%R(i,iType)%Free()
          enddo
       enddo
       deallocate(self%R)
    endif
    if(associated(self%A))then
       do i = 1, size(self%A,1)
          call self%A(i)%Free()
       enddo
       deallocate(self%A)
    endif
    if(associated(self%B))then
       do i = 1, size(self%B,1)
          call self%B(i)%Free()
       enddo
       deallocate(self%B)
    endif
    if(allocated(self%A_IS_COMPUTED))deallocate(self%A_IS_COMPUTED)
    self%B_IS_COMPUTED=.FALSE.
    self%SA2=-1
    self%SB2=-1
    self%S2=-1
    self%Initialized = .FALSE.
  end subroutine ClassTDM1Free

  integer function GET_TDM1B_TYPE(strn) result(iType)
    character(len=*), intent(in) :: strn
    do iType=1,N_REDUCED_TDM1B_TYPES
       if(trim(adjustl(strn))==REDUCED_TDM1B_TYPES(iType))exit
    enddo
  end function GET_TDM1B_TYPE

  integer function GET_TDM2B_TYPE(strn) result(iType)
    character(len=*), intent(in) :: strn
    do iType=1,N_REDUCED_TDM2B_TYPES
       if(trim(adjustl(strn))==REDUCED_TDM2B_TYPES(iType))exit
    enddo
  end function GET_TDM2B_TYPE

  subroutine ClassTDM1Init( self, IrrepAB, SA2, SB2, nactive )

    class(ClassTDM1), intent(inout) :: self
    type(ClassIrrep), target, intent(in) :: IrrepAB
    integer         , intent(in)    :: SA2, SB2
    integer         , intent(in)    :: nactive(*)

    integer :: iIrrBra, iIrrKet, nIrr, nbra, nket
    type(ClassIrrep), pointer :: irreplist(:)
    type(ClassIrrep), pointer :: irrepKet

    if(SA2<0.or.SB2<0)then
       write(ERROR_UNIT,"(a)") "Wrong multiplicity passed to TDM1 Init"
       stop
    endif

    call self%free()
    self%Initialized=.FALSE.

    allocate(self%IrrepAB,source=IrrepAB)
    self%SA2 = SA2
    self%SB2 = SB2
    self%S2  = -1

    if( abs( SA2 - SB2 ) > 2 ) return

    nIrr = GlobalGroup%GetNIrreps()
    allocate(self%nactive(nIrr))
    self%nactive = nactive(1:nIrr)
    irreplist => GlobalGroup%GetIrrepList()
    if(SA2==SB2)then
       allocate(self%A(nIrr))
       allocate(self%A_IS_COMPUTED(nIrr))
       self%A_IS_COMPUTED=.FALSE.
    endif
    allocate(self%B(nIrr))
    self%B_IS_COMPUTED=.FALSE.
    allocate(self%R(nIrr,N_REDUCED_TDM1B_TYPES))

    do iIrrBra = 1, nIrr
       irrepKet => IrrepAB * irreplist( iIrrBra )
       iIrrKet = GlobalGroup%GetIrrepIndex( irrepKet )
       nbra = nactive(iIrrBra)
       nket = nactive(iIrrKet)
       call self%B(iIrrBra)%InitFull(nbra,nket)
       if(SA2==SB2) call self%A(iIrrBra)%InitFull(nbra,nket)
       if(SA2==SB2) call self%R(iIrrBra,GET_TDM1B_TYPE("0"))%InitFull(nbra,nket)
       if(SA2+SB2>0)call self%R(iIrrBra,GET_TDM1B_TYPE("1"))%InitFull(nbra,nket)
    enddo

  end subroutine ClassTDM1Init

  subroutine ClassTDM1SetBlock( self, iIrr1, Sigma2, nr, nc, RhoUp, LDU, RhoDown, LDD )

    class(ClassTDM1), intent(inout) :: self
    integer         , intent(in)    :: iIrr1, Sigma2, nr, nc, LDU, LDD
    real(kind(1d0)) , intent(in)    :: RhoUp(LDU,*), RhoDown(LDD,*)
    !
    real(kind(1d0)) :: f1, f2, fc, Sa, Sb, Sigma
    integer :: iType

    Sa = dble(self%SA2)/2.d0
    Sb = dble(self%SB2)/2.d0
    f1 = sqrt(Sa+0.5d0)

    if( self%SA2 == self%SB2 )then
       iType = GET_TDM1B_TYPE("0")
       self%R(iIrr1,iType)%A = RhoUp(1:nr,1:nc) + RhoDown(1:nr,1:nc)
       self%R(iIrr1,iType)%A = f1 * self%R(iIrr1,iType)%A
    endif

    Sigma = dble(Sigma2)/2.d0
    if( self%SA2 == self%SB2 .and. self%SA2 /= 0 .and. Sigma2 == 0 )then
       write(ERROR_UNIT,"(a)") "Invalid Spin projection for the ions, in the calculation of R1 TDM"
       stop
    endif

    if( self%SA2 + self%SB2 > 0 )then
       iType = GET_TDM1B_TYPE("1")
       call clebsch(Sb,1.d0,Sa,Sigma,0.d0,Sigma,fc)
       f2 = f1 / fc
       self%R(iIrr1,iType)%A = RhoUp(1:nr,1:nc) - RhoDown(1:nr,1:nc)
       self%R(iIrr1,iType)%A = f2 * self%R(iIrr1,iType)%A
    endif

    self%Initialized=.TRUE.

  end subroutine ClassTDM1SetBlock

  subroutine ClassTDM1GetRTBlock( self, iIrr1, strnT, Mat )
    class(ClassTDM1) , intent(in) :: self
    integer          , intent(in) :: iIrr1
    character(len=*) , intent(in) :: strnT
    type(ClassMatrix), intent(out):: Mat
    integer :: nr, nc, iType

    iType = GET_TDM1B_TYPE(strnT)
    nr=self%R(iIrr1,iType)%NRows()
    nc=self%R(iIrr1,iType)%NColumns()
    if(nr==0.or.nc==0)return    
    call Mat%InitFull(nr,nc)
    if(allocated(self%R(iIrr1,iType)%A))then
       Mat = self%R(iIrr1,iType)%A
    else
       Mat = 0.d0
    endif
  end subroutine ClassTDM1GetRTBlock

  function   ClassTDM1GetRTBlockFun( self, iIrr1, strnT) result(Mat)
    class(ClassTDM1) , intent(in) :: self
    integer          , intent(in) :: iIrr1
    character(len=*) , intent(in) :: strnT
    type(ClassMatrix), pointer :: Mat
    integer :: nr, nc, iType

    iType = GET_TDM1B_TYPE(strnT)
    nr=self%R(iIrr1,iType)%NRows()
    nc=self%R(iIrr1,iType)%NColumns()
    if(nr==0.or.nc==0)return    
    call Mat%InitFull(nr,nc)
    !if(allocated(self%R(iIrr1,iType)%A))then
       Mat => self%R(iIrr1,iType)
    !else
       !Mat => 0.d0
    !endif
  end function ClassTDM1GetRTBlockFun

  subroutine ClassTDM1ComputeA( self, iIrr )
    class(ClassTDM1), intent(inout) :: self
    integer         , intent(in)    :: iIrr
    self%A(iIrr)%A = self%R(iIrr,GET_TDM1B_TYPE("0"))%A / sqrt(dble(self%SA2)/2.d0+0.5d0)
    self%A_IS_COMPUTED(iIrr)=.TRUE.
  end subroutine ClassTDM1ComputeA

  subroutine ClassTDM1GetABlock( self, iIrr, Mat )
    class(ClassTDM1)             , intent(inout):: self
    integer                      , intent(in)   :: iIrr
    real(kind(1d0)) , allocatable, intent(out)  :: Mat(:,:)
    integer :: nr, nc
    nr=self%A(iIrr)%NRows()
    nc=self%A(iIrr)%NColumns()
    if(nr==0.or.nc==0)return    
    call realloc(Mat,nr,nc,EXPAND=.TRUE.)
    if(.not.self%A_IS_COMPUTED(iIrr))call self%ComputeA(iIrr)
    Mat(1:nr,1:nc) = self%A(iIrr)%A
  end subroutine ClassTDM1GetABlock

  function   ClassTDM1GetABlockFun( self, iIrr ) result(Mat)
    class(ClassTDM1)             , intent(inout):: self
    integer                      , intent(in)   :: iIrr
    type(ClassMatrix), pointer :: Mat
    integer :: nr, nc
    nr=self%A(iIrr)%NRows()
    nc=self%A(iIrr)%NColumns()
    if(nr==0.or.nc==0)return    
    if(.not.self%A_IS_COMPUTED(iIrr))call self%ComputeA(iIrr)
    Mat => self%A(iIrr)
  end function ClassTDM1GetABlockFun
  
  subroutine ClassTDM1ComputeB( self, S2 )

    class(ClassTDM1), intent(inout) :: self
    integer         , intent(in)    :: S2

    integer         :: iIrr, nIrr, iType0, iType1
    real(kind(1d0)) :: sgn, Sa, Sb, S
    self%S2 = S2

    Sa = dble(self%SA2) / 2.d0
    Sb = dble(self%SB2) / 2.d0
    S  = dble(self%S2 ) / 2.d0
    nIrr = GlobalGroup%GetNIrreps()
    iType0=GET_TDM1B_TYPE("0")
    iType1=GET_TDM1B_TYPE("1")
    do iIrr = 1, nIrr
       self%B(iIrr)%A = 0.d0
       if(allocated(self%R(iIrr, iType0)%A))  &
            self%B(iIrr)%A = self%B(iIrr)%A + &
            SixJSymbol_HalfHalf( 0.d0, Sb, Sa, S ) * self%R(iIrr, iType0)%A
       if(allocated(self%R(iIrr, iType1)%A))  &
            self%B(iIrr)%A = self%B(iIrr)%A + &
            sqrt(3.d0) * SixJSymbol_HalfHalf( 1.d0, Sb, Sa, S ) * self%R(iIrr, iType1)%A
       sgn = 1.d0 - 2.d0 * mod(nint(dble(self%SB2+self%S2+1)/2.d0)+1,2)
       self%B(iIrr)%A = sgn * self%B(iIrr)%A
    enddo
    self%B_IS_COMPUTED=.TRUE.
  end subroutine ClassTDM1ComputeB

  subroutine ClassTDM1GetBBlock( self, iIrr, S2, Mat )
    class(ClassTDM1)             , intent(inout):: self
    integer                      , intent(in)   :: iIrr
    integer                      , intent(in)   :: S2
    real(kind(1d0)) , allocatable, intent(out)  :: Mat(:,:)
    integer :: nr, nc
    nr=self%A(iIrr)%NRows()
    nc=self%A(iIrr)%NColumns()
    if(nr==0.or.nc==0)return    
    call realloc(Mat,nr,nc,EXPAND=.TRUE.)
    if( .not. ( S2 == self%S2 .AND. self%B_IS_COMPUTED ) ) call self%ComputeB(S2)
    Mat(1:nr,1:nc) = self%B(iIrr)%A
  end subroutine ClassTDM1GetBBlock

  function   ClassTDM1GetBBlockFun( self, iIrr ) result(Mat)
    class(ClassTDM1)             , intent(inout):: self
    integer                      , intent(in)   :: iIrr
    type(ClassMatrix), pointer :: Mat
    integer :: nr, nc
    nr=self%B(iIrr)%NRows()
    nc=self%B(iIrr)%NColumns()
    if(nr==0.or.nc==0)return    
    if(.not.self%B_IS_COMPUTED)call self%ComputeB(iIrr)
    Mat => self%B(iIrr)
  end function ClassTDM1GetBBlockFun

  subroutine ClassTDM1WriteToUnit( self, uid )

    class(ClassTDM1), intent(inout) :: self
    integer         , intent(in)    :: uid

    integer :: nIrr, iIrr, iType
    character(len=16) :: Form
    logical :: Formatted, Initialized

    nIrr=GlobalGroup%GetNIrreps()
    !
    INQUIRE( UNIT = uid, FORM  = Form )
    Formatted = trim(Form) .is. "FORMATTED" 
    if(Formatted)then
       write(uid,"(*(x,i4))") self%SA2, self%SB2, nIrr
       write(uid,"(*(x,i4))") self%nactive
    else
       write(uid) self%SA2, self%SB2, nIrr
       write(uid) self%nactive
    end if
    !
    call self%IrrepAB%Write(uid)
    do iIrr = 1, nIrr
       do iType = 1, N_REDUCED_TDM1B_TYPES
          Initialized = self%R(iIrr,iType)%IsInitialized()
          if(Formatted)then
             write(uid,*) Initialized
          else
             write(uid)   Initialized
          endif
          if( Initialized )then
             call self%R(iIrr,iType)%Write(uid)
          endif
       enddo
    enddo
  end subroutine ClassTDM1WriteToUnit

  subroutine ClassTDM1ReadFromUnit( self, uid )

    class(ClassTDM1), intent(inout) :: self
    integer         , intent(in)    :: uid

    integer :: nIrr, iIrr, iType
    character(len=16) :: Form
    logical :: Formatted, Initialized

    integer :: SA2, SB2
    integer, allocatable :: nactive(:)
    type(ClassIrrep) :: IrrepAB

    INQUIRE( UNIT = uid, FORM  = Form )
    Formatted = trim(Form) .is. "FORMATTED" 
    if(Formatted)then
       read(uid,*) SA2, SB2, nIrr
       allocate(nactive(nIrr))
       nactive=0
       read(uid,*) nactive
    else
       read(uid) SA2, SB2, nIrr
       allocate(nactive(nIrr))
       nactive=0
       read(uid) nactive
    end if
    !
    call IrrepAB%Read(uid)
    call IrrepAB%SetGroup(GlobalGroup)
    call self%init(IrrepAB,SA2,SB2,nactive)
    do iIrr = 1, nIrr
       do iType = 1, N_REDUCED_TDM1B_TYPES
          if(Formatted)then
             read(uid,*) Initialized
          else
             read(uid)   Initialized
          endif
          if(Initialized)then
             call self%R(iIrr,iType)%Read(uid)
          endif
       enddo
    enddo
    self%Initialized = .TRUE.
    deallocate(nactive)
  end subroutine ClassTDM1ReadFromUnit

  subroutine ClassTDM2Final( self )
    type(ClassTDM2), intent(inout) :: self
    call self%Free()
  end subroutine ClassTDM2Final

  subroutine ClassTDM2Free( self )
    class(ClassTDM2), intent(inout) :: self
    integer :: i, j, k, s, jtk
    call self%IrrepAB%Free()
    if( allocated(self%nactive)) deallocate(self%nactive)
    if(associated(self%Pi))then
       do jtk = 1, N_REDUCED_TDM2B_TYPES
          do i = 1, size(self%Pi,1)
             do j = 1, size(self%Pi,2)
                do k = 1, size(self%Pi,3)
                   call self%Pi(i,j,k,jtk)%Free()
                enddo
             enddo
          enddo
       enddo
       deallocate(self%Pi)
    endif
    if(associated(self%Q))then
       do i = 1, size(self%Q,1)
          do j = 1, size(self%Q,2)
             do k = 1, size(self%Q,3)
                do s = 1, size(self%Q,4)
                   call self%Q(i,j,k,s)%Free()
                enddo
             enddo
          enddo
       enddo
       deallocate(self%Q)
    endif
    self%SA2=0
    self%SB2=0
    self%S2 =0
    if(allocated(self%Q_IS_COMPUTED))deallocate(self%Q_IS_COMPUTED)
    self%Initialized  =.FALSE.
  end subroutine ClassTDM2Free

  subroutine ClassTDM2Init( self, IrrepAB, SA2, SB2, nactive )

    class(ClassTDM2), intent(inout) :: self
    type(ClassIrrep), target, intent(in)    :: IrrepAB
    integer         , intent(in)    :: SA2, SB2 
    integer         , intent(in)    :: nactive(*)

    integer :: iIrr_s, iIrr_r, iIrr_p, iIrr_q, nIrr, iStot, S2
    integer :: ndimr, ndims, ndimp, ndimq, jtk, J, T, K, nScouplings
    type(ClassIrrep), pointer :: irreplist(:)
    type(ClassIrrep), pointer :: irrepAB1, irrepAB12, irrepAB123

    if(SA2<0.or.SB2<0)then
       write(ERROR_UNIT,"(a)") "Wrong multiplicity passed to TDM2 Init"
       stop
    endif

    call self%free()
    self%Initialized = .FALSE.

    self%IrrepAB = IrrepAB
    self%SA2 = SA2
    self%SB2 = SB2
    self%S2  = -1

    nScouplings=1
    if(SA2==SB2.and. SA2/=0)nScouplings=2

    if( abs( SA2 - SB2 ) > 2 ) return

    nIrr = GlobalGroup%GetNIrreps()
    allocate(self%nactive(nIrr))
    self%nactive = nactive(1:nIrr)
    irreplist => GlobalGroup%GetIrrepList()

    allocate(self%Pi(nIrr,nIrr,nIrr,N_REDUCED_TDM2B_TYPES))
    allocate(self%Q( nIrr,nIrr,nIrr,nScouplings))
    allocate(self%Q_IS_COMPUTED(nIrr,nIrr,nIrr,nScouplings))
    do iIrr_p = 1, nIrr
       irrepAB1 => IrrepAB * irreplist( iIrr_p )
       ndimp = nactive(iIrr_p)
       if(ndimp==0)cycle
       do iIrr_r = 1, nIrr
          irrepAB12 => irrepAB1 * irreplist( iIrr_r )
          ndimr = nactive(iIrr_r)
          if(ndimr==0)cycle
          do iIrr_s = 1, nIrr
             irrepAB123 => irrepAB12 * irreplist( iIrr_s )
             ndims = nactive(iIrr_s)
             if(ndims==0)cycle
             iIrr_q = GlobalGroup%GetIrrepIndex( irrepAB123 )
             ndimq = nactive(iIrr_q)
             if(ndimq==0)cycle
             do jtk = 1, N_REDUCED_TDM2B_TYPES
                call GetJTKindexes( jtk, J, T, K )
                if( 2 * K < abs( SA2 - SB2 ) .or. 2 * K > SA2 + SB2 ) cycle 
                call self%Pi(iIrr_s,iIrr_r,iIrr_p,jtk  )%Init(ndims,ndimr,ndimp,ndimq)
                iStot = 0
                do S2 = max(abs(SA2-1),abs(SB2-1)), min(SA2+1,SB2+1), 2
                   iStot=iStot+1
                   call self%Q( iIrr_s,iIrr_r,iIrr_p,iStot)%Init(ndims,ndimr,ndimp,ndimq)
                enddo
             enddo
          enddo
       enddo
    enddo
    self%Q_IS_COMPUTED=.FALSE.
  end subroutine ClassTDM2Init

  subroutine GetJTKindexes(iJTK,J,T,K)
    integer, intent(in) :: iJTK
    integer, intent(out):: J, T, K
    character(len=*), parameter :: LIST_INTEGERS = "012"
    character(len=3):: sJTK
    sJTK = REDUCED_TDM2B_TYPES(iJTK)
    J = index(LIST_INTEGERS, sJTK(1:1)) - 1
    T = index(LIST_INTEGERS, sJTK(2:2)) - 1
    K = index(LIST_INTEGERS, sJTK(3:3)) - 1
  end subroutine GetJTKindexes

  !> Converts the uncoupled TDM2B to coupled reduced form
  !>  \f[
  !>      \mathsf{\Pi}^{BA}_{[[sr]_J,[pq]_{T}]_K}&=\frac{(-1)^{J} \Pi_{S_A}}{C_{S_B\Sigma,K0}^{S_A\Sigma}}\times\\
  !> \times\Big[&
  !>   - C_{T1,J-1}^{K0}C_{\frac{1}{2}\frac{1}{2},\frac{1}{2}\frac{1}{2}}^{J1}
  !>     C_{\frac{1}{2}\frac{1}{2},\frac{1}{2}\frac{1}{2}}^{T 1}\pi^{BA}_{r_\alpha s_\alpha,p_\alpha q_\alpha}\\
  !>  &- C_{T-1,J1}^{K0}C_{\frac{1}{2}-\frac{1}{2},\frac{1}{2}-\frac{1}{2}}^{J-1}
  !>     C_{\frac{1}{2}-\frac{1}{2},\frac{1}{2}-\frac{1}{2}}^{T -1}\pi^{BA}_{r_\beta s_\beta,p_\beta q_\beta}\\
  !>  & -C_{T 0,J0}^{K0}C_{\frac{1}{2}-\frac{1}{2},\frac{1}{2}\frac{1}{2}}^{J0}
  !>     C_{\frac{1}{2}-\frac{1}{2},\frac{1}{2}\frac{1}{2}}^{T 0}\pi^{BA}_{r_\alpha s_\beta, q_\alpha p_\beta}\\
  !>  &  C_{T 0,J0}^{K0}C_{\frac{1}{2}-\frac{1}{2},\frac{1}{2}\frac{1}{2}}^{J0}
  !>     C_{\frac{1}{2}\frac{1}{2},\frac{1}{2}-\frac{1}{2}}^{T 0}\pi^{BA}_{r_\alpha s_\beta,p_\alpha q_\beta}\\
  !>  &  C_{T 0,J0}^{K0}C_{\frac{1}{2}\frac{1}{2},\frac{1}{2}-\frac{1}{2}}^{J0}
  !>     C_{\frac{1}{2}-\frac{1}{2},\frac{1}{2}\frac{1}{2}}^{T 0}\pi^{BA}_{s_\alpha r_\beta,q_\alpha p_\beta}\\
  !>  & -C_{T 0,J0}^{K0}C_{\frac{1}{2}\frac{1}{2},\frac{1}{2}-\frac{1}{2}}^{J0}
  !>     C_{\frac{1}{2}\frac{1}{2},\frac{1}{2}-\frac{1}{2}}^{T 0}\pi^{BA}_{s_\alpha r_\beta,p_\alpha q_\beta}
  !>     \Big]
  !>
  subroutine ClassTDM2SetBlock( self, iIrr1, iIrr2, iIrr3, iSigma2, PIAB, PIAA, PIBB )

    class(ClassTDM2)   , intent(inout) :: self
    integer            , intent(in)    :: iIrr1, iIrr2, iIrr3, iSigma2
    type(ClassMatrix4D), intent(in)    :: PIAB, PIAA, PIBB

    real(kind(1d0)) :: sgn, GeneralFactor, c1, c2, c3, alpha
    real(kind(1d0)) :: Sa, Sb, SigmaB, SigmaA
    integer         :: n1, n2, n3, n4, i1, i2, i3, i4 
    integer         :: K, T, J, iJTK, iIrr4, S2
    real(kind(1d0)), allocatable :: dmat2134_AA(:,:,:,:)
    real(kind(1d0)), allocatable :: dmat2134_BB(:,:,:,:)
    real(kind(1d0)), allocatable :: dmat2134_AB(:,:,:,:)
    real(kind(1d0)), allocatable :: dmat2143_AB(:,:,:,:)
    real(kind(1d0)), allocatable :: dmat1234_AB(:,:,:,:)
    real(kind(1d0)), allocatable :: dmat1243_AB(:,:,:,:)

    iIrr4 = GlobalGroup%GetIrrepIndex( self%IrrepAB * &
         GlobalGroup%GetIrrep(iIrr1) * &
         GlobalGroup%GetIrrep(iIrr2) * &
         GlobalGroup%GetIrrep(iIrr3) )

    Sa   = dble(self%SA2)/2.d0
    Sb   = dble(self%SB2)/2.d0

    SigmaA = 0.5d0 * iSigma2
    SigmaB = 0.5d0 * iSigma2

    n1 = PIAB%size(1)
    n2 = PIAB%size(2)
    n3 = PIAB%size(3)
    n4 = PIAB%size(4)

    allocate(dmat2134_AA(n2,n1,n3,n4))
    allocate(dmat2134_BB(n2,n1,n3,n4))
    allocate(dmat2134_AB(n2,n1,n3,n4))
    allocate(dmat2143_AB(n2,n1,n4,n3))
    allocate(dmat1234_AB(n1,n2,n3,n4))
    allocate(dmat1243_AB(n1,n2,n4,n3))

    dmat2134_AA=0.d0
    dmat2134_BB=0.d0
    dmat2134_AB=0.d0
    dmat2143_AB=0.d0
    dmat1234_AB=0.d0
    dmat1243_AB=0.d0

    do i1 = 1, n1
       do i2 = 1, n2
          dmat2134_AA(i2,i1,:,:) = PIAA%A( i1, i2, :, : )
          dmat2134_AB(i2,i1,:,:) = PIAB%A( i1, i2, :, : )
          dmat2134_BB(i2,i1,:,:) = PIBB%A( i1, i2, :, : )
          do i3 = 1, n3
             dmat2143_AB(i2,i1,:,i3) = PIAB%A( i1, i2, i3, : )
          enddo
       enddo
    enddo
    dmat1234_AB=PIAB%A
    do i3 = 1, n3
       do i4 = 1, n4
          dmat1243_AB(:,:,i4,i3)=PIAB%A(:,:,i3,i4)
       enddo
    enddo

    do iJTK = 1, N_REDUCED_TDM2B_TYPES

       call getJTKindexes( iJTK, J, T, K )

       if( 2 * K < abs( self%SA2 - self%SB2 ) .or. 2 * K > self%SA2 + self%SB2 ) cycle 

       !.. Computes the general factor in the expression of the reduced coupled TDM2B
       !   in terms of the uncoupled TDM2B
       !..
       call clebsch(Sb, 1.d0*K, Sa,  0.5d0*iSigma2, 0.d0, 0.5d0*iSigma2, c1)
       if(abs(c1)<1.d-10)then
          write(*,*) "HORROR! VANISHING CLEBSCH GORDAN IN TDM2B CONVERSION!"
          stop
       endif
       sgn = (1.d0-2.d0*mod(J,2))
       GeneralFactor = - sgn * sqrt(2.d0*Sa+1.d0) / c1

       if( abs( self%SA2 - self%SB2 ) > 2*K ) cycle

       if( J == 1 .and. T==1 )then

          call clebsch(1.d0*T, 1.d0*J, 1.d0*K,  1.d0, -1.d0, 0.d0, c1)
          call clebsch( 0.5d0,  0.5d0, 1.d0*J, 0.5d0, 0.5d0, 1.d0, c2)
          call clebsch( 0.5d0,  0.5d0, 1.d0*T, 0.5d0, 0.5d0, 1.d0, c3)
          alpha = c1 * c2 * c3 * GeneralFactor

          !.. Formula 1
          call self%Pi( iIrr2, iIrr1, iIrr3, iJTK )%Add( alpha, dMat2134_AA )

          !.. Formula 2
          sgn = 1.d0 - 2.d0 * mod(K,2)
          call self%Pi( iIrr2, iIrr1, iIrr3, iJTK )%Add( sgn * alpha, dMat2134_BB )

       endif

       call clebsch(1.d0*T, 1.d0*J, 1.d0*K,  0.d0,  0.d0, 0.d0, c1)
       call clebsch( 0.5d0,  0.5d0, 1.d0*J,-0.5d0, 0.5d0, 0.d0, c2)
       call clebsch( 0.5d0,  0.5d0, 1.d0*T,-0.5d0, 0.5d0, 0.d0, c3)
       alpha = c1 * c2 * c3 * GeneralFactor 

       !.. Formula 3
       call self%Pi( iIrr2, iIrr1, iIrr4, iJTK )%Add( alpha, dMat2143_AB )

       !.. Formula 4
       sgn = 1.d0 - 2.d0 * mod(T,2)
       call self%Pi( iIrr2, iIrr1, iIrr3, iJTK )%Add( sgn * alpha, dMat2134_AB )

       !.. Formula 5
       sgn = 1.d0 - 2.d0 * mod(J,2)
       call self%Pi( iIrr1, iIrr2, iIrr4, iJTK )%Add( sgn * alpha, dMat1243_AB )

       !.. Formula 6
       sgn = 1.d0 - 2.d0 * mod(T+J,2)
       call self%Pi( iIrr1, iIrr2, iIrr3, iJTK )%Add( sgn * alpha, dMat1234_AB )

    enddo

    !.. We must be careful here since**** subroutine ClassTDM2ComputeQ( self, iIrr_p, iIrr_q, iIrr_s , S2 ) cannot be called here
    !   since self%Pi(iIrr1, iIrr2, iIrr3) is not completely filled. It should work in the following way.
    do S2 = max(self%SA2-1,self%SB2-1), min(self%SA2+1,self%SB2+1)
       if(S2.ge.0)then
          call self%ComputeQ(iIrr2, iIrr3, iIrr1, S2) !..  Formulas 1 and 2 partial contributions.
          call self%ComputeQ(iIrr2, iIrr4, iIrr1, S2) !..  Formula  3       partial contribution.
          call self%ComputeQ(iIrr2, iIrr3, iIrr1, S2) !..  Formula  4       partial contribution.
          call self%ComputeQ(iIrr1, iIrr4, iIrr2, S2) !..  Formula  5       partial contribution.
          call self%ComputeQ(iIrr1, iIrr3, iIrr2, S2) !..  Formula  6       partial contribution. 
       endif
    enddo
    
    deallocate(dMat2134_AA,dMat2134_BB,dMat2134_AB,dMat2143_AB,dMat1234_AB,dMat1243_AB)

    self%Initialized = .TRUE.

  end subroutine ClassTDM2SetBlock

  subroutine ClassTDM2GetPiBlock( self, iIrr_s, iIrr_r, iIrr_p, strnT, Mat )

    class(ClassTDM2)             , intent(in)  :: self
    integer                      , intent(in)  :: iIrr_s, iIrr_r, iIrr_p
    character(len=*)             , intent(in) :: strnT
    ! character(len=*)             , intent(in)  :: sBlock
    real(kind(1d0)) , allocatable, intent(out):: Mat(:,:,:,:)
    integer :: N1, N2, N3, N4, i, j, k, l, iType

    N1=self%Pi(iIrr_s, iIrr_r, iIrr_p, 1)%Size(1)
    N2=self%Pi(iIrr_s, iIrr_r, iIrr_p, 1)%Size(2)
    N3=self%Pi(iIrr_s, iIrr_r, iIrr_p, 1)%Size(3)
    N4=self%Pi(iIrr_s, iIrr_r, iIrr_p, 1)%Size(4)


    if(N1==0.or.N2==0.or.N3==0.or.N4==0)return
    call realloc(Mat,N1,N2,N3,N4,EXPAND=.TRUE.)
    iType = GET_TDM2B_TYPE(strnT)
    do l = 1, N4
       do k = 1, N3
          do j = 1, N2
             do i = 1, N1
                Mat(i,j,k,l) = self%Pi(iIrr_s, iIrr_r, iIrr_p, iType)%A(i,j,k,l)
             enddo
          enddo
       enddo
    enddo

  end subroutine ClassTDM2GetPiBlock

  subroutine ClassTDM2ComputeQ( self, iIrr_p, iIrr_q, iIrr_s , S2 )

    class(ClassTDM2)   , intent(inout) :: self
    integer            , intent(in)    :: iIrr_p,iIrr_q,iIrr_s
    integer            , intent(in)    :: S2
    real(kind(1d0)) :: factor, fase, piKTJ
    real(kind(1d0)) :: S, Sa, Sb
    integer         :: i4,i3,i2
    integer         :: K, T, J, jtk, iS
    character(len=3):: jtkCh

    Sa   = dble(self%SA2)/2.d0
    Sb   = dble(self%SB2)/2.d0
    S    = dble(S2)/2.d0
    iS=1
    if( self%SA2==self%SB2 .and. self%SA2/=0 .and. S2==self%SA2+1 ) iS=2

    fase = -1.d0 * (1.d0 - 2.d0 * mod( nint(S+Sb+0.5d0), 2 ) )

    do jtk = 1, N_REDUCED_TDM2B_TYPES
       jtkCh = REDUCED_TDM2B_TYPES(jtk)
       read (jtkCh(1:1),'(I10)') J
       read (jtkCh(1:1),'(I10)') T
       read (jtkCh(1:1),'(I10)') K
       !
       piKTJ = sqrt((2.d0*T+1.d0)*(2.d0*J+1.d0)*(2.d0*K+1.d0))
       !.. Properties of the the six-j symbol
       !  sjs_hh( c, d, e, f )
       !  | 1/2  1/2  c |       |  d    e   c  |       |  d    c   e  |
       ! <               >  =  <                >  =  <                >
       !  |  d   e    f |       | 1/2  1/2  f  |       | 1/2   f  1/2 |
       factor = piKTJ * fase * &
            SixJSymbol_HalfHalf( K*1.d0, Sa    , Sb    , S     ) * &
            SixJSymbol_HalfHalf( K*1.d0, J*1.d0, T*1.d0, 0.5d0 )
       if(mod(T+K,2)==1) factor = - factor
       do i4 = 1, size(self%Pi(iIrr_p,iIrr_s,iIrr_q,jtk)%A,4)
          do i3 = 1, size(self%Pi(iIrr_p,iIrr_s,iIrr_q,jtk)%A,3)
             do i2 = 1, size(self%Pi(iIrr_p,iIrr_s,iIrr_q,jtk)%A,2)
                self%Q(iIrr_p,iIrr_q,iIrr_s,iS)%A(:,i3,i2,i4) = &
                     self%Q(iIrr_p,iIrr_q,iIrr_s,iS)%A(:,i3,i2,i4) + &
                     factor * self%Pi(iIrr_p,iIrr_s,iIrr_q,jtk)%A(:,i2,i3,i4)
             enddo
          enddo
       enddo
    enddo
    self%Q_IS_COMPUTED(iIrr_p,iIrr_q,iIrr_s,iS)=.TRUE.
  end subroutine ClassTDM2ComputeQ

  function   ClassTDM2GetQBlockFun( self, iIrr_s,iIrr_r,iIrr_p,S2 ) result(Mat)
    class(ClassTDM2)             , intent(inout):: self
    integer                      , intent(in)   :: iIrr_s,iIrr_r,iIrr_p
    integer                      , intent(in)   :: S2
    type(ClassMatrix4D), pointer :: Mat
    integer :: iS
    iS=1
    if( self%SA2==self%SB2 .and. self%SA2/=0 .and. S2==self%SA2+1 ) iS=2
    if(.not.self%Q_IS_COMPUTED(iIrr_s,iIrr_r,iIrr_p,iS))call self%ComputeQ(iIrr_s,iIrr_r,iIrr_p,S2)
    Mat => self%Q(iIrr_s,iIrr_r,iIrr_p,iS)
  end function ClassTDM2GetQBlockFun

  subroutine ClassTDM2WriteToUnit( self, uid )
    class(ClassTDM2), intent(inout) :: self
    integer         , intent(in)    :: uid
    integer :: nIrr, iIrr_s, iIrr_r, iIrr_p, jtk
    character(len=16) :: Form
    logical :: Formatted
    nIrr=size(self%nactive)

    INQUIRE( UNIT = uid, FORM  = Form )
    Formatted = trim(Form) .is. "FORMATTED" 
    if(Formatted)then
       write(uid,"(*(x,i4))") self%SA2, self%SB2, nIrr
       write(uid,"(*(x,i4))") self%nactive
    else
       write(uid) self%SA2, self%SB2, nIrr
       write(uid) self%nactive
    end if

    call self%IrrepAB%Write(uid)
    do iIrr_s = 1, nIrr
       do iIrr_r = 1, nIrr
          do iIrr_p = 1, nIrr
             do jtk = 1, N_REDUCED_TDM2B_TYPES
                if(self%Pi(iIrr_s,iIrr_r,iIrr_p,jtk)%TotalSize()<=0)cycle
                if(Formatted)then
                   write(uid,"(*(x,i4))") iIrr_s, iIrr_r, iIrr_p, jtk
                else
                   write(uid) iIrr_s, iIrr_r, iIrr_p, jtk
                endif
                call self%Pi(iIrr_s,iIrr_r,iIrr_p,jtk)%Write(uid)
             enddo
          enddo
       enddo
    enddo
  end subroutine ClassTDM2WriteToUnit

  subroutine ClassTDM2ReadFromUnit( self, uid )
    class(ClassTDM2), intent(inout) :: self
    integer         , intent(in)    :: uid
    integer :: nIrr, iIrr_s, iIrr_r, iIrr_p, jtk
    integer :: SA2, SB2, iostat
    integer, allocatable :: nactive(:)
    type(ClassIrrep) :: IrrepAB
    character(len=16) :: Form
    logical :: Formatted

    INQUIRE( UNIT = uid, FORM  = Form )
    Formatted = trim(Form) .is. "FORMATTED" 
    if(Formatted)then
       read(uid,*) SA2, SB2, nIrr
       allocate(nactive(nIrr))
       nactive=0
       read(uid,*) nactive
    else
       read(uid) SA2, SB2, nIrr
       allocate(nactive(nIrr))
       nactive=0
       read(uid) nactive
    end if
    call IrrepAB%Read(uid)
    call IrrepAB%SetGroup(GlobalGroup)
    call self%init(IrrepAB,SA2,SB2,nactive)
    do
       if(Formatted)then
          read(uid,*,iostat=iostat) iIrr_s, iIrr_r, iIrr_p, jtk
       else
          write(uid,iostat=iostat) iIrr_s, iIrr_r, iIrr_p, jtk
       endif
       if(iostat/=0)exit
       call self%Pi(iIrr_s,iIrr_r,iIrr_p,jtk)%Read(uid)
    enddo
    self%initialized = .TRUE.
    deallocate(nactive)

  end subroutine ClassTDM2ReadFromUnit

end module ModuleDensityMatrices
