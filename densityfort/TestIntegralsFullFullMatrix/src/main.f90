! {{{ Detailed description

!> \mainpage Program <ProgramName> <Insert here what the program does>
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
program TestIntegrals

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleString
  use ModuleGroups
  use ModuleESpace
  use ModuleXlm

  use ModuleBspline
  use ModuleBasisJUAN
  use ModuleParametersJUAN
  use ModuleDensityMatricesJUAN
  !  use bspline_grid_mod
  use ModuleIntegrals
  use ModuleMainInterface
  !.. 
  ! use mpi_mod 
  ! use ukrmol_interface
  ! use symmetry
  use precisn
  ! use ModuleUKRmolInterface
  !..

  implicit none

  !.. Run-time parameters
  character(len=:), allocatable :: ProgramInputFile
  character(len=:), allocatable :: ccConfigFile
  character(len=:), allocatable :: AssignedSymmetry

  !.. Config file parameters
  character(len=:), allocatable :: StorageDir
  character(len=:), allocatable :: QCDir

  type(ClassESpace)             :: Space
  type(ClassGroup), pointer     :: Group
  type(ClassIrrep), pointer     :: TotIrrepPtr
  integer                       :: TotIrrepIndex


  !.. Local parameters
  type(ClassBspline)                  :: Bspline,Bsplinex
  integer                             :: ia,i,j,k,l,iflag
  real(kind(1d0))                     :: ivalue,ivalue2,ivalue3,r
  real(kind(1d0))                     :: H0d, time1, time2
  type(BasisElementInfo), allocatable :: ielement(:),ielement0(:),iprim(:)
  real(kind(1d0))                     :: OneBody, Biel, Overlap, fun0, fun1, fun2
  integer                             :: Ndim_sym
  integer, allocatable                :: SymIndexes(:)


  !  External Functions
  integer, external           :: Kron,irr_lm
  real(kind(1d0)), external   :: GET_INTEGRALS_E,GET_INTEGRALS_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX,TwoBody_STEX
  real(kind(1d0)), external   :: OneBody_STEX_ASTRA,TwoBody_STEX_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX_NO_ASTRA,TwoBody_STEX_NO_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX_TDMATRIX ,TwoBody_STEX_TDMATRIX  


  call GetRunTimeParameters( ProgramInputFile, AssignedSymmetry )

  call ParseProgramInputFile( &
       ProgramInputFile     , &
       StorageDir           , &
       QCDir                , &
       ccConfigFile         )

  call Space%SetRootDir  ( AddSlash(StorageDir) )

  !call Space%SetNuclearLabel( QCDir )
  call Space%ParseConfigFile( ccConfigFile , 1)

  Group => Space%GetGroup()
  call GlobalGroup%Init(Group%GetName())
  call GlobalXlmSet%Init(Group,lmax)

  TotIrrepPtr   => GlobalGroup%GetIrrep(AssignedSymmetry)
  TotIrrepIndex =  GlobalGroup%getIrrepIndex(TotIrrepPtr)
  write(*,*) "TotIrrepIndex",TotIrrepIndex

  !.. Initialize the parameters from UKRmol file
  call InitFromFiles
  write(*,*) "lmax",lmax

  !.. This routine looks set variables for the change of irreps order fro ukrmol to astra ones 
 
  !  call UKRmol_Orbitals%Init()
  !  call UKRMol_Init_Integrals()


  !Here we define bsplines which correspond to an extension of the ukrmol bspline basis.
  !If we want to set the same basis we must consider the following:
  !The number  "no_bsplines" in scatci_integrals input corresponds to NNodes+BsplineOrder-2
  ! "a", the R-matrix radius is our Rmax (bspline_grid_start = 0.0)
  !  bspline_order is BsplineOrder. The  bspline_indices(1,l)  and bspline_indices(2,0) are the first and last bspline elements considered. We must remove the first one for l=0 (bspline_indices(1,0)=2) and the first and second from l>0  (bspline_indices(1,l>0)=3). bspline_indices(2,l) is chossen in order to exclude the bsplines which reach the radious "a", corresponding to the last BsplineOrder-1, so it has to be set to bspline_indices(2,l)=no_bsplines-(bspline_order-1) if we want to remove them from the calculation.
  !Once we do the calculation with ukrmol we have to set our bsplines in order to exactly match the bsplines up to bspline_indices(2,l). This can be set by making a grid in a larger domain where it mathc the ukrmol internal grid. There is an overlap region defined by the last bspline_order-1 orbitals considered in ukrmol.
  !The bsplines which are not in the internal basis and do overlap with the last elements of the inner basis, the elment number no_bsplines-bspline_order-1, are the ones which go from no_bsplines-bspline_order to no_bsplines-2. Since this last element of the inner basis is connected to  all basis elements due to the orthonormalization processes, these lastmentioned elements are coupled with all the inner ukrmol elements.

  !Initialize external bsplines values
  NBspl      = 0
  NlastBspl  = -1
  NfirstBspl = 0
  NExtBspl   = 0
  Rmax = a  !our new redius. If we do not want the external bsplines we set Rmax=a, ELSE we define a larger radius:
  Rmax = 30.1333333333333d0
  !.. Initialize both internal (Bspline) and external (Bsplinex) Bsplines
  call InitBsplines(Bspline )
  call InitBsplinex(Bsplinex)
  
  write(*,*)  NBspl
  write(*,*)  NlastBspl 
  write(*,*)  NfirstBspl
  write(*,*)  NExtBspl

  write(*,*) "===================================="
  call GlobalIntegral%SetStorage( AddSlash(StorageDir) )
  call GlobalIntegral%Setlmax( lmax )
  ! call GlobalIntegral%Init( lmax, nLOC, nSRC , no_bsplines, BsplineOrder , Bsplinex)
  call GlobalIntegral%ReadFromFile()
  write(*,*) "===================================="


  !.. Set up basis parameters
  call InitBasisParameters

  !.. Initialize the object that identify the ukrmol basis elements (symmetry, l, m, nr, etc)
  allocate(iprim(1:NprimUkrmol))
  call InitIprim(iprim)
  !Here we organize the one particle basis elements
  Nbasis=sum(irr_total)*2
  write(*,*) "New dimension of the basis:",sum(irr_total)
  write(*,*) "New dimension of the basis including spin :",Nbasis

  !.. Initialization of the object that comprises the information of the channel (hole, excited state, etc)
  allocate(ielement(1:Nbasis),ielement0(sum(nTOT+nPWC)))
  call InitIelement(ielement,ielement0)

  !.. Overlap matrix needed for the extended basis
  allocate(Spq(1:Nbasis,1:Nbasis))
  Spq=0.d0
  Do j=1,Nbasis
     Do i=1,Nbasis
        Spq(j,i)=GET_INTEGRALS_E(j,i,0,0,1,Bsplinex,ielement,iprim)
     End Do
  End Do
  call STEX_BASIS(ielement)

  !.. Now we build the system of equations
  allocate(Hirr(1:Ndim,1:Ndim),Over(1:Ndim,1:Ndim),En(1:Ndim))  
  Over=0.d0
  Hirr=0.d0
  En=0.d0

  !.. Hartree Fock Energy
  NucRep=(7*7.d0/(2*0.54875d0))*(1.d0/1.8897259886)
  write(*,*) "Nuclear repulsion",NucRep
  !NucRep=23.6261352

  !.. HF energy
  iflag=4
  !ivalue2=OneBody_STEX_ASTRA(1,1,iflag,Bsplinex,ielement,iprim)
  !ivalue3=TwoBody_STEX_ASTRA(1,1,Bsplinex,ielement,iprim)
  
  ivalue2=OneBody_STEX_TDMATRIX(1,1,iflag,Bsplinex,ielement,iprim)
  ivalue3=TwoBody_STEX_TDMATRIX(1,1,Bsplinex,ielement,iprim)

  ivalue=ivalue2+ivalue3
  write(*,*) "Hartree Fock Energy",ivalue+NucRep,ivalue2,ivalue3,NucRep
  write(*,*) "--------------"

  Onebody = 0.d0
  Biel    = 0.d0
  Overlap = 0.d0

  
  iflag=4

  !.. File with the matrix elements comming from previous program, which read the ukrmol integrals and evaluates
  !   the Bs contribution on the fly
  open(unit=101,file="ELEMENTOS.dat", ACTION='READ')
  
  Do j=1,Ndim
     !
     call cpu_time(time1)
     write(*,*) "j",j,Ndim
     
     Do i=1,Ndim
        Over(i,j)=1.d0*Kron(i,j)

        !.. STEX MATRIX using ukrmol basis only
        !Overlap = (1.d0/dble(NRefConf))*OneBody_STEX(j,i,1,Bsplinex,ielement,iprim)   
        !OneBody = OneBody_STEX(j,i,iflag,Bsplinex,ielement,iprim)
        !Biel    = TwoBody_STEX(j,i,Bsplinex,ielement,iprim)

        !.. STEX MATRIX using ukrmol plus extended (optional)
        Overlap = (1.d0/dble(NRefConf))*OneBody_STEX_NO_ASTRA(j,i,1,Bsplinex,ielement,iprim)
        OneBody = OneBody_STEX_NO_ASTRA(j,i,iflag,Bsplinex,ielement,iprim)
        Biel    = TwoBody_STEX_NO_ASTRA(j,i,Bsplinex,ielement,iprim)
        
        !.. STEX TDM  using ukrmol basis plus extended (optional)
        !Overlap = (1.d0/dble(NRefConf))*OneBody_STEX_TDMATRIX(j,i,1,Bsplinex,ielement,iprim)
        !OneBody = OneBody_STEX_TDMATRIX(j,i,iflag,Bsplinex,ielement,iprim)
        !Biel    = TwoBody_STEX_TDMATRIX(j,i,Bsplinex,ielement,iprim)


        !================================================================
        !.. Comparison of the matrix elements between different programs
        read(101,*) k, l
        if((j.ne.k).or.(i.ne.l))then
           write(*,*) "dimension of the cmatrices are not the same",j,k,"----",i,l
           stop
        endif
        read(101,*) fun0,fun1,fun2
        if((abs(Onebody+Biel-fun0-fun1)+abs(Overlap-fun2)).gt.0.00000001d0)then
           write(*,*) "in pair j,i",j,i,(abs(Onebody+Biel-fun0-fun1)+abs(Overlap-fun2))
           write(*,'(a,3F24.13)') "Our values :",Onebody,Biel,Overlap
           write(*,'(a,3F24.13)') "in the file:",fun0,fun1,fun2
           pause
        endif
        !================================================================
        ! write(*,*) j,i
        ! write(*,*) Onebody,Biel,Overlap,"Ham:",Onebody+Biel
        ! pause
        
        Over(i,j) = Overlap
        Hirr(j,i) = OneBody + Biel + Overlap*NucRep
        if(isnan(Hirr(j,i)))then
           write(*,*)Over(i,j)
           write(*,*) j,i,iflag,(1.d0/dble(NRefConf))
           stop
        end if

     End Do
     call cpu_time(time2)
     write(*,*) "Time",(time2-time1),(time2-time1)*Ndim/3600.d0
  End Do
  
  Close(101)


  !check symmetry
  H0d=0.d0
  ivalue=0.d0
  Do j=1,Nbasis
     Do i=1,Nbasis
        H0d=H0d+abs(Hirr(j,i)-Hirr(i,j))
        ivalue=ivalue+abs(Over(j,i)-Over(i,j))
     End DO
  End Do
  write(*,*) "assymetry", H0d,ivalue*0.5/Nbasis,Nbasis



  !.. Evaluate the eigenvalues
  LWORK=Ndim*3-1
  allocate(WORK(1:LWORK))
  WORK=0.d0
  INFO=0
  call dsygv(1,'N','U',Ndim,Hirr,Ndim,Over,Ndim,En,WORK,LWORK,INFO)
  write(*,*) "energies INFOpp",INFO
  !pause
  Do i=1,10!Ndim
     Write(*,*) "eigen",i,En(i)
  End do
  !.. File for the comparision of the energies
  !open(unit=2, file="energy_ukrmol.dat")
  open(unit=2, file="energy_ukrmol_extended_NO_ASTRA.dat")
  !open(unit=2, file="energy_ukrmol_extended_TDMATRIX.dat")
  Do i=1,Ndim
     write(2,*) En(i)
  End Do
  close(2)
  deallocate(WORK)

  !.. Comparison of the energies
  open(unit=1, file="energy_ukrmol_reference.dat")

  !open(unit=2, file="energy_ukrmol.dat")
  open(unit=2, file="energy_ukrmol_extended_NO_ASTRA.dat")
  !open(unit=2, file="energy_ukrmol_extended_TDMATRIX.dat")
  H0d=0.d0
  r=0.d0
  Do i=1,Ndim
     read(2,*) ivalue2
     read(1,*) ivalue
     If(abs(ivalue-ivalue2).gt.r)then
        r=abs(ivalue-ivalue2)
        ia=i
     End IF
     H0d=H0d+abs(ivalue-ivalue2)
  End Do
  close(1)
  close(2)

  write(*,*) "difference between eigenvalues",H0d
  write(*,*) "worst term",ia,r




  call mpi_mod_finalize




end program TestIntegrals
