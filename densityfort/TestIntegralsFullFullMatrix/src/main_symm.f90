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

  !.. Config file parameters
  character(len=:), allocatable :: StorageDir
  character(len=:), allocatable :: QCDir

  type(ClassESpace)             :: Space
  type(ClassGroup), pointer     :: Group



  !.. Local parameters
  type(ClassBspline)                  :: Bspline,Bsplinex
  integer                             :: ia,i,j,k,l,iflag
  real(kind(1d0))                     :: ivalue,ivalue2,ivalue3,r
  real(kind(1d0))                     :: H0d, time1, time2
  type(BasisElementInfo), allocatable :: ielement(:),ielement0(:),iprim(:)
  real(kind(1d0))                     :: OneBody, Biel, Overlap, fun0, fun1, fun2


  !  External Functions
  integer, external           :: Kron,irr_lm
  real(kind(1d0)), external   :: GET_INTEGRALS_E,GET_INTEGRALS_Eno,GET_INTEGRALS_Xi,GET_INTEGRALS_ASTRA
  real(kind(1d0)), external   :: GET_1B_INTEGRAL_ASTRA,GET_1B_INTEGRAL_UKRMOL,GET_1B_INTEGRAL_EXTENDED
  real(kind(1d0)), external   :: GET_2B_INTEGRAL_ASTRA,GET_2B_INTEGRAL_EXTENDED
  real(kind(1d0)), external   :: OneBody_STEX,TwoBody_STEX
  real(kind(1d0)), external   :: OneBody_STEX_ASTRA,TwoBody_STEX_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX_NO_ASTRA,TwoBody_STEX_NO_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX_TDMATRIX ,TwoBody_STEX_TDMATRIX  

  call GetRunTimeParameters( ProgramInputFile )

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
  !.. 
  IrrepList = GlobalGroup%GetIrrepList()

  
  call InitFromFiles
  write(*,*) "lmax",lmax
  !pause

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
  call InitBsplines(Bspline )
  call InitBsplinex(Bsplinex)
  !make another subroutine for initializing the external bsplines
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


  !pause
  call InitBasisParameters

  allocate(iprim(1:NprimUkrmol))
  call InitIprim(iprim)
  !Here we organize the one particle basis elements
  Nbasis=sum(irr_total)*2 !irr_total(1)!sum(irr_total)*2!irr_total(1)!irr_total(1)*2      ! sum(irr_total)*2   !        !Spin
  write(*,*) "New dimension of the basis:",sum(irr_total)
  write(*,*) "New dimension of the basis including spin :",Nbasis




  !Here we organize the one particle basis elements
  Nbasis=sum(irr_total)*2 !irr_total(1)!sum(irr_total)*2!irr_total(1)!irr_total(1)*2      ! sum(irr_total)*2   !        !Spin
  write(*,*) "New dimension of the basis:",sum(irr_total)
  !write(*,*) "New dimension of the basis including spin :",Nbasis

  !NExtOrb=0
  allocate(ielement(1:Nbasis),ielement0(sum(nTOT+nPWC)))
  call InitIelement(ielement,ielement0)


  ! .. Basis orthonormalization: Xi
  ! write(*,*) "inn"
  ! call Get_Xi(Bsplinex,ielement,iprim)
  ! write(*,*) "out"

  ! !.. Transformation Matrix 
  ! allocate(TMa(1:Ndim,1:Ndim))
  ! allocate(arr(1:Ndim))
  ! arr=0.d0
  ! TMa=0.d0
  ! TMa(1,1)=1.d0
  ! Do j=2,Ndim  !UkXi
  !    Do i=2,Ndim  !UkBE
  !       If(iha(i).eq.iha(j))then
  !          If(ielement(iep(i))%spin.eq.ielement(iep(j))%spin)Then
  !             TMa(i,j)=Xi(iep(i),iep(j))
  !          End IF
  !       End IF
  !    End Do
  ! End Do

! iflag=4
! fun0=GET_INTEGRALS_ASTRA(155,155,0,0,iflag,Bsplinex,ielement,iprim)
! write(*,*) fun0
! write(*,*) "---------------"
! fun0=GET_INTEGRALS_E(155,155,0,0,iflag,Bsplinex,ielement,iprim)
! write(*,*) fun0

! stop


  allocate(Spq(1:Nbasis,1:Nbasis))
  Spq=0.d0
  Do j=1,Nbasis
     Do i=1,Nbasis
        Spq(j,i)=GET_INTEGRALS_E(j,i,0,0,1,Bsplinex,ielement,iprim)
     End Do
  End Do
  !pause
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

  ! iflag=UKRMol_Get_Integral_Index('Property integrals')
  ! write(*,*) "iflag",iflag
  ! ivalue = GET_1B_INTEGRAL_ASTRA(1,1,iflag,0,0)
  ! write(*,*) ivalue,sqrt(4*3.1415d0)
  ! ivalue = GET_1B_INTEGRAL_UKRMOL(1,1,iflag,0,0)
  ! write(*,*) ivalue,sqrt(4*3.1415d0)
  ! ivalue = GET_1B_INTEGRAL_EXTENDED(1,1,iflag,0,0, Bsplinex, iprim)
  ! write(*,*) ivalue




  iflag=4
  ivalue2=OneBody_STEX_ASTRA(1,1,iflag,Bsplinex,ielement,iprim)
  write(*,*) "ivalue2",ivalue2
  ivalue2=OneBody_STEX_TDMATRIX(1,1,iflag,Bsplinex,ielement,iprim)
  write(*,*) "ivalue2",ivalue2
  !pause
  !ivalue3=TwoBody_STEX_ASTRA(1,1,Bsplinex,ielement,iprim)
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
 

  !.. Determine dimension submatrix with assigned symmetry (which you will from command line)
  !..
  integer :: Ndim_sym
  integer, allocatable :: SymIndexes(:)
  Ndim_sym=0
  do j=1,Ndim
     if(Symmetry(j) == AssignedSymmetry ) Ndim_sym=Ndim_sym+1
  enddo
  allocate(SymIndexes(Ndim_sym))
  i=0
  do j=1,Ndim
     if(Symmetry(j) == AssignedSymmetry )then
        SymIndexes(i)=j
     endif
  enddo
     
  !.. Define here Hirr_sym, etc.

  Do j_sym=1,Ndim_sym
     call cpu_time(time1)
     j= SymIndexes(j_sym)
     write(*,*) "j_sym, j",j_sym,j,Ndim_sym
          
     Do i_sym=1,Ndim_sym
        Over(i,j)=1.d0*Kron(i,j)
        i=SymIndexes(i_sym)

        !.. STEX MATRIX using ukrmol basis only
        !Overlap = (1.d0/dble(NRefConf))*OneBody_STEX(j,i,1,Bsplinex,ielement,iprim)   
        !OneBody = OneBody_STEX(j,i,iflag,Bsplinex,ielement,iprim)
        !Biel    = TwoBody_STEX(j,i,Bsplinex,ielement,iprim)

        !.. STEX TDM  using ukrmol basis plus extended (optional)
        !Overlap = (1.d0/dble(NRefConf))*OneBody_STEX_TDMATRIX(j,i,1,Bsplinex,ielement,iprim)
        !OneBody = OneBody_STEX_TDMATRIX(j,i,iflag,Bsplinex,ielement,iprim)
        !Biel    = TwoBody_STEX_TDMATRIX(j,i,Bsplinex,ielement,iprim)

        !.. STEX MATRIX using ukrmol plus extended (optional)
        Overlap = (1.d0/dble(NRefConf))*OneBody_STEX_NO_ASTRA(j,i,1,Bsplinex,ielement,iprim)
        OneBody = OneBody_STEX_NO_ASTRA(j,i,iflag,Bsplinex,ielement,iprim)
        Biel    = TwoBody_STEX_NO_ASTRA(j,i,Bsplinex,ielement,iprim)


        
        

        ! ivalue2=TwoBody_STEX_TDMATRIX(j,i,Bsplinex,ielement,iprim)
        ! ivalue3=TwoBody_STEX_NO_ASTRA(j,i,Bsplinex,ielement,iprim)
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        read(101,*) k, l
        if((j.ne.k).or.(i.ne.l))then
           write(*,*) "dimension of the cmatrices are not the same",j,k,"----",i,l
           stop
        endif
        read(101,*) fun0,fun1,fun2
        if((abs(Onebody+Biel-fun0-fun1)+abs(Overlap-fun2)).gt.0.001d0)then
           write(*,*) "in pair j,i",j,i,(abs(Onebody+Biel-fun0-fun1)+abs(Overlap-fun2))
           
           write(*,*) "Our values :",Onebody,Biel,Overlap
           write(*,*) "in the file:",fun0,fun1,fun2
           pause
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        ! write(*,*) j,i
        ! write(*,*) Onebody,Biel,Overlap,"Ham:",Onebody+Biel
        ! pause
        Over_sym(i_sym,j_sym) = Overlap
        Hirr_sym(i_sym,j_sym) = OneBody + Biel + Overlap*NucRep
        if(isnan(Hirr(j_sym,i_sym)))then
           write(*,*)Over(i_sym,j_sym)
           write(*,*) j,i,iflag,(1.d0/dble(NRefConf))
           stop
        end if
        ! if(i.eq.j)then
        !    write(*,*) i,j,Hirr(j,i),Over(i,j)
        !    endif
     End Do
     call cpu_time(time2)
     write(*,*) "Time",(time2-time1),(time2-time1)*Ndim/3600.d0
     !pause
  End Do
  Close(101)


  !check hermiticity
  H_Av_Err=sum(abs(Hirr-transpose(Hirr))/dble(Ndim_sym)**2
  S_Av_Err=sum(abs(Over-transpose(Over))/dble(Ndim_sym)**2
  H_oo_Err=maxval(abs(Hirr-transpose(Hirr))
  S_oo_Err=maxval(abs(Over-transpose(Over))
  write(*,*) "H average and peak asymmetry", H_Av_Err, H_oo_Err
  write(*,*) "S average and peak asymmetry", S_Av_Err, S_oo_Err


  !.. Save basis and matrices 
  do i=1,Ndim_sym
     write(uid_basis,"(*(i0))") i,hole_sym,hole_n,Orb_kind,Orb_sym,Orb_lm,Orb_index
  enddo
  type(ClassMatrix)::Smat,Hmat
  Smat=Over
  Hmat=Hirr
  call Smat%write("test/Smat")
  call Hmat%write("test/Hmat")

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
  !.. File for the comparision of the energyes
  !open(unit=2, file="energy_ukrmol.dat")
  open(unit=2, file="energy_ukrmol_extended_NO_ASTRA.dat")
  Do i=1,Ndim
     write(2,*) En(i)
  End Do
  close(2)
  deallocate(WORK)

  !.. Comparison of the energies
  open(unit=1, file="energy_ukrmol_reference.dat")

  !open(unit=2, file="energy_ukrmol.dat")
  open(unit=2, file="energy_ukrmol_extended.dat")
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

