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

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO
  use ModuleGroups
  use ModuleSymESpace
  use ModuleESpace
  use ModuleXlm
  use ModuleConstants
  use ModuleDiagonalize
  use ModuleBspline
  use ModuleAngularMomentum
  use ModuleBasisJUAN
  use ModuleParametersJUAN
  use ModuleDensityMatricesJUAN
  use bspline_grid_mod
  use ModuleIntegrals
  use ModuleUKRmolInterface
  use ModuleMolecularGeometry
  use ModuleMainInterface
  !.. 
  use mpi_mod 
  use ukrmol_interface
  use symmetry
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

  type(ClassESpace)    :: Space
  type(ClassGroup), pointer     :: Group
  type(ClassMolecularGeometry)  :: MolGeom

  
  integer*8                 :: NFTINT
  integer*8                 :: irrdim
  type(ClassBspline)        :: Bspline,Bsplinex
  
  !.. Local parameters
  integer*8        :: in8

  !.. this refers to ukrmol interface variables
  type(bspline_grid_obj) :: bspline_grid
  integer*8 :: lunit, posit
  integer*8 :: lb, bspline_index, number_of_functions, pos_after_rw
  real(kind=cfp)  :: norm
  logical*8 :: non_zero_at_boundary
  
  integer         :: ia,ib,ic,id,i,j,k,l,m,lp,mp,i1,i2,i3,i4,iflag
  real(kind(1d0)) :: ivalue,ivalue2,ivalue3,dr,r,time1,time2
  
  integer         :: iIrr1, iIrr2
  
  integer                        :: IOSTAT
  real(kind(1d0))                :: Integral
  character(len=12)              :: GridFile
  character(len=12)              :: IrrName
  character(len=2)               :: PB
  integer                 :: indx
  real(kind(1d0))         :: H0d
  type(BasisElementInfo), allocatable :: ielement(:),ielement0(:),iprim(:)

  integer                        ::  LowInBspl
  integer                        ::  LasExBspl
  type(ClassDensityMatrices)             :: STEX_DM

!  External Functions
  integer, external :: Kron,irr_lm
  real(kind(1d0)), external   :: GET_INTEGRALS_E,GET_INTEGRALS_Eno,GET_INTEGRALS_Xi
  real(kind(1d0)), external   :: GET_1B_INTEGRAL_ASTRA,GET_1B_INTEGRAL_UKRMOL,GET_1B_INTEGRAL_EXTENDED
  real(kind(1d0)), external   :: GET_2B_INTEGRAL_ASTRA,GET_2B_INTEGRAL_EXTENDED
  real(kind(1d0)), external   :: OneBody_STEX,TwoBody_STEX
  real(kind(1d0)), external   :: OneBody_STEX_ASTRA,TwoBody_STEX_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX_NO_ASTRA,TwoBody_STEX_NO_ASTRA
  real(kind(1d0)), external   :: OneBody_STEX_TDMATRIX ,TwoBody_STEX_TDMATRIX  

! MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT

! BASIS
!   use precisn
!   use ukrmol_interface
!   use mpi_mod
  
!  Module ModuleDensityMatricesJUAN

!   use ModuleBasisJUAN

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

  call InitFromFiles
  write(*,*) "lmax",lmax
!pause

  
!  call UKRmol_Orbitals%Init()
!  call UKRMol_Init_Integrals()



! nIrreps = GlobalGroup%GetnIrreps()
  ! do iIrr1 = 1, nIrreps
  !    do iIrr2 = 1, nIrreps
  !       write(*,*) iIrr1,iIrr2,allocated(GlobalIntegral%Int_1B(1)%LOC_LOC(iIrr1,iIrr2)%A)
  !    end do
  ! enddo

 
!  call UKRmol_Orbitals%Init()
!  call UKRMol_Init_Integrals()


  ! irr_series  = nTOT
  ! irr_counter = sTOT
  ! irrdim=size(irr_series)
  ! Nukrmol=sum(irr_series)
  ! Do i=1,8
  !    write(*,*) irr_series(i)
  ! End Do   
  ! write(*,*) "Nukrmol",Nukrmol
  ! write(*,*) "ukrmol    ^ ^ ^ ^ "

!pause

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


! j=45
! i=273
! write(*,*) Nbasis,j,i
! ivalue=GET_INTEGRALS_E(j,i,0,0,1,Bsplinex,ielement,iprim)
! !pause


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
!ivalue2=OneBody_STEX_ASTRA(1,1,iflag,Bsplinex,ielement,iprim)
ivalue2=OneBody_STEX_TDMATRIX(1,1,iflag,Bsplinex,ielement,iprim)
!ivalue3=TwoBody_STEX_ASTRA(1,1,Bsplinex,ielement,iprim)
ivalue3=TwoBody_STEX_TDMATRIX(1,1,Bsplinex,ielement,iprim)

ivalue=ivalue2+ivalue3
write(*,*) "Hartree Fock Energy",ivalue+NucRep,ivalue2,ivalue3,NucRep
!


!.. STEX MATRIX using ukrmol basis only
iflag=4

Do j=6,Ndim
   !if(mod(j,100).eq.0)then
   !call cpu_time(time1)
      write(*,*) "j",j,Ndim
    !  endif
   Do i=48,Ndim
      Over(i,j)=1.d0*Kron(i,j)

       !Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX(j,i,1,Bsplinex,ielement,iprim)   
     ! ivalue=OneBody_STEX(j,i,iflag,Bsplinex,ielement,iprim)
!       ivalue=ivalue+TwoBody_STEX(j,i,Bsplinex,ielement,iprim)

      write(*,*) "j,i",j,i
       Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX_TDMATRIX(j,i,1,Bsplinex,ielement,iprim)   
       ivalue =OneBody_STEX_TDMATRIX(j,i,iflag,Bsplinex,ielement,iprim)
       ! ivalue2=OneBody_STEX_NO_ASTRA(  j,i,iflag,Bsplinex,ielement,iprim)
       ! if(abs(ivalue-ivalue2).gt.0.0001d0)then
       !    write(*,*) "ONEBODY"
       !    write(*,*) ivalue,ivalue2, j,i
       !    pause
       ! endif
       

       ivalue3=TwoBody_STEX_TDMATRIX(j,i,Bsplinex,ielement,iprim)
       ! ivalue2=TwoBody_STEX_NO_ASTRA(j,i,Bsplinex,ielement,iprim)
       !  if(abs(ivalue2-ivalue3).gt.0.0001d0)then
       !     write(*,*) "TWOBODY"
       !     write(*,*) ivalue2,ivalue3, j,i
       !     pause
       !  endif
       
        ivalue=ivalue+ivalue3
       
       !Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX_ASTRA(j,i,1,Bsplinex,ielement,iprim)   
       !ivalue=OneBody_STEX_ASTRA(j,i,iflag,Bsplinex,ielement,iprim)
       !ivalue=ivalue+TwoBody_STEX_ASTRA(j,i,Bsplinex,ielement,iprim)

      !.. The STEX_NO_ASTRA includes the STEX_ASTRA as a particular case. Perhaps it is slower.
      !Over(i,j)=(1.d0/dble(NRefConf))*OneBody_STEX_NO_ASTRA(j,i,1,Bsplinex,ielement,iprim)
      !ivalue=OneBody_STEX_NO_ASTRA(j,i,iflag,Bsplinex,ielement,iprim)
      !ivalue=ivalue+TwoBody_STEX_NO_ASTRA(j,i,Bsplinex,ielement,iprim)
      
       Hirr(j,i)=ivalue+Over(i,j)*NucRep
       if(isnan(Hirr(j,i)))then
          write(*,*)Over(i,j)
          write(*,*) j,i,iflag,(1.d0/dble(NRefConf))
          stop
       end if
       ! if(i.eq.j)then
       !    write(*,*) i,j,Hirr(j,i),Over(i,j)
       !    endif
    End Do
    !call cpu_time(time2)
    !write(*,*) "Time",time2-time1
   !pause
End Do

!pause
! write(*,*) Ndim
! open(unit=10,file="matrizH.dat")
! do j=1,Ndim
! do i=1,Ndim
!    write(10,*) Hirr(j,i)
! End Do
! End Do

  
!check symmetry
H0d=0.d0
ivalue=0.d0
Do j=1,Nbasis!NprimUkrmol
   Do i=1,Nbasis!NprimUkrmol
      H0d=H0d+abs(Hirr(j,i)-Hirr(i,j))
      !      write(*,*) j,i,abs(Hirr(j,i)-Hirr(i,j))
      ivalue=ivalue+abs(Over(j,i)-Over(i,j))
   End DO
   !      write(*,*) H0d
   !     pause
End Do
write(*,*) "assymetry", H0d,ivalue*0.5/Nbasis,Nbasis


! Hirr=matmul(Hirr,TMa)
! Hirr=matmul(transpose(TMa),Hirr)
! Over=matmul(Over,TMa)
! Over=matmul(transpose(TMa),Over)


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
!open(unit=2, file="energy_ukrmol.dat")
open(unit=2, file="energy_ukrmol_extended.dat")
Do i=1,Ndim
   write(2,*) En(i)
End Do
close(2)
deallocate(WORK)
!pause

!.. Comparison of the energies
open(unit=1, file="energy_ukrmol_reference.dat")


!open(unit=2, file="energy_ukrmol.dat")
open(unit=2, file="energy_ukrmol_extended.dat")
H0d=0.d0
r=0.d0
Do i=1,Ndim
   read(2,*) ivalue2
   read(1,*) ivalue
   !write(*,*) ivalue-ivalue2,iha(i),iep(i),i
   If(abs(ivalue-ivalue2).gt.r)then
      r=abs(ivalue-ivalue2)
      ia=i
   End IF
   H0d=H0d+abs(ivalue-ivalue2)
End Do
close(1)
close(2)

 write(*,*) "difference between eigenvalues",H0d
 write(*,*) ia,r


  

call mpi_mod_finalize


!!$  call GetRunTimeParameters( FileName, nSize )

stop

!contains






  
end program TestIntegrals


!!! EXTENDED 1B INTEGRALS EXTENDED 1B INTEGRALS EXTENDED 1B INTEGRALS EXTENDED 1B INTEGRALS  
!  iflag=UKRMol_Get_Integral_Index('Property integrals')
!  write(*,*) "iflag",iflag
!  iflag=3
!  i=1
!  j=8
!  l=0
!  m=1
! do i=1,sum(nTOT+nPWC*1)
!    do j=1,sum(nTOT+nPWC*1)
!       !if((absTOchar(i).ne."LOC").or.(absTOchar(j).ne."LOC"))cycle
!       do l=0,2
!          do m=-l,l!
!      ! if((absTOchar(i).eq."LOC").and.(absTOchar(j).eq."LOC"))then  !for the moments
!          ivalue  = GET_1B_INTEGRAL_ASTRA( i, j, iflag,l,m)
!          ivalue2 = GET_1B_INTEGRAL_EXTENDED( i, j, iflag,l,m, Bsplinex, iprim)
!          write(*,*) absTOchar(i),absTOchar(j),ivalue,ivalue2
!          if(abs(ivalue-ivalue2).gt.0.0001d0)then
!             write(*,*) l,m
!             write(*,*) i,j,ivalue,ivalue2
!             write(*,*) "<<<<<<<<"
!             write(*,*) absTOchar(i),absTOirr(i),absTOrel(i)
!             write(*,*) absTOchar(j),absTOirr(j),absTOrel(j),absTOl(j),absTOm(j),absTOnr(j)
!             write(*,*) ">>>>>>>>>"
!             stop
!          endif
!       !endif
!       end do
!     enddo
!  enddo
! write(*,*) "No problems 1B flag:",iflag
! pause

 ! ivalue = GET_1B_INTEGRAL_EXTENDED( 45,45, 5,0,0, Bsplinex, iprim )
 ! write(*,*) ivalue
 ! pause
! i=3
! j=136
! k=53
! l=155

! do i=187,sum(nTOT+nPWC)
!    do j=1,sum(nTOT+nPWC)
!       write(*,*) i,j
!       do k=1,sum(nTOT+nPWC)
!          do l=1,sum(nTOT+nPWC)
! ia=0
! if(absTOchar(i).eq."LOC")ia=ia+1
! if(absTOchar(j).eq."LOC")ia=ia+1
! if(absTOchar(k).eq."LOC")ia=ia+1
! if(absTOchar(l).eq."LOC")ia=ia+1
! if(ia.lt.2)cycle

! if(ia.lt.2)then
!    write(*,*) "pero ia lt 2"
! endif
            
!             ivalue = GET_2B_INTEGRAL_ASTRA( i, j, k, l)
!             ivalue2 = GET_2B_INTEGRAL_EXTENDED( i, j, k, l, Bsplinex, iprim)
!             if(abs(ivalue-ivalue2).gt.0.00001d0)then
!                write(*,*) "-----------------",ivalue,ivalue2
!                write(*,*) i,j,k,l
!                stop
!             endif
!             !write(*,*) ivalue,ivalue-ivalue2
!          end do
!       end do
!    end do
! end do
! write(*,*) "No problems 2B"
!stop
!!! EXTENDED 1B INTEGRALS EXTENDED 1B INTEGRALS EXTENDED 1B INTEGRALS EXTENDED 1B INTEGRALS  


!tires plus

! UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY
! iflag=UKRMol_Get_Integral_Index('Property integrals')
! write(*,*) "iflag",iflag
! do i=1,sum(nTOT+nPWC)
!    do j=1,sum(nTOT+nPWC)
!       if((absTOchar(i).ne."LOC").or.(absTOchar(j).ne."LOC"))cycle
!       do l=0,2
!          do m=-l,l
!             write(*,*) l,m,i,j

!             ivalue = GET_1B_INTEGRAL_ASTRA( i, j, iflag,l,m)
!             !ivalue2 = GET_1B_INTEGRAL_UKRMOL( i, j, iflag,l,m)
!             ivalue2 =GET_INTEGRALS_E(i,j,0,0,iflag,Bsplinex,ielement,iprim)
!             if(abs(ivalue-ivalue2).gt.0.0001d0)then
!                write(*,*) i,j,l,m,ivalue,ivalue2
               
!                write(*,*) "<<<<<<<<",i,j
!                write(*,*) absTOrel(i),nLOC(absTOirr(i)),absTOrel(j),nLOC(absTOirr(j))
!                write(*,*) ">>>>>>>>>"
!                stop
!             endif
!             !write(*,*) ivalue,ivalue2
!          end do
!       end do
!    enddo
! enddo
! write(*,*) "No problems 1B flag:",iflag
! pause
! do i=1,sum(nTOT)
!    do j=1,sum(nTOT)
!       write(*,*) i,j
!       do k=1,sum(nTOT)
!          do l=1,sum(nTOT)
            
!             ivalue = GET_2B_INTEGRAL_ASTRA( i, j, k, l)
!             ivalue2 = GET_2B_INTEGRAL_EXTENDED( i, j, k, l)
!             if(abs(ivalue-ivalue2).gt.0.d0)then
!                write(*,*) i,j,k,l,ivalue,ivalue2
!                stop
!             endif
!             !write(*,*) ivalue,ivalue-ivalue2
!          end do
!       end do
!    end do
! end do
! write(*,*) "No problems 2B"
! pause
! UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY UKRMOL ONLY
