!> Reads the run time parameters specified in the command line.
! Flags: flag=0  one body hamiltonian
!        flag=1  overlap
function OneBody_STEX(i,j,flag,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleBspline
  use ModuleSystemUtils
  use mpi_mod 
  use ukrmol_interface 
  use precisn
  use ModuleIntegrals
  implicit none
  integer, intent(in) :: i,j,flag
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,lk
  integer*8 :: k,i8flag
  real(kind(1d0)) :: ivalue,ivalue2,Integral
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL

  !write(*,*) "Sdas",i,j,flag
  ME=0.d0
  i8flag=int(flag,kind(i8flag))
  !id=0
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           k=ielement(ia)%RCukLoc
           
           Integral=GET_INTEGRALS(k,k,0,0,i8flag)
           !Integral=GET_1B_INTEGRAL_ASTRA(k,k,flag,0,1)
          ! write(*,*) k,Integral
           ivalue=ivalue+Integral
           !id=id+1
        End IF
     End Do
    ! write(*,*) "Sdas",ivalue
     ! pause
     ME=ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     lk=max(i,j)

       ME=GET_INTEGRALS_E(iep(lk),iha(lk),0,0,flag,Bsplinex,ielement,iprim)
  ElseIf((i.gt.1).and.(j.gt.1))then
     If(iha(i).eq.iha(j))then
      
       ME=GET_INTEGRALS_E(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)
     End IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     If(i.eq.j)then
        ivalue2=0.d0
        Do ia=1,Nbasis
           If(ielement(ia)%RCtype)then
              k=ielement(ia)%RCukLoc
              ivalue2=ivalue2+GET_INTEGRALS(k,k,0,0,i8flag)
              !ivalue2=ivalue2+GET_1B_INTEGRAL_ASTRA(k,k,flag,0,1)
              !ivalue=GET_1B_INTEGRAL_ASTRA(k,k,flag,0,1)
              !write(*,*) k,ivalue
           End IF
        End Do
        ME=ME+ivalue2
     ENd IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     If(iep(i).eq.iep(j))then
       ivalue=GET_INTEGRALS_E(iha(j),iha(i),0,0,flag,Bsplinex,ielement,iprim)
       ME=ME-ivalue
     End IF
  End If

End function OneBody_STEX


function OneBody_STEX_ASTRA(i,j,flag,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleBspline
  use ModuleSystemUtils
  use mpi_mod 
  use ukrmol_interface 
  use precisn
  use ModuleIntegrals
  implicit none
  integer, intent(in) :: i,j,flag
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,lk
  integer :: k,i8flag
  real(kind(1d0)) :: ivalue,ivalue2,Integral
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL
  

  ME=0.d0
  i8flag=int(flag,kind(i8flag))
  !id=0q
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           k=ielement(ia)%RCukLoc
           Integral=GET_INTEGRALS_ASTRA(ia,ia,0,0,flag,Bsplinex,ielement,iprim)
           ivalue=ivalue+Integral
           !id=id+1
        End IF
     End Do
     ME=ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     lk=max(i,j)

       ME=GET_INTEGRALS_ASTRA(iep(lk),iha(lk),0,0,flag,Bsplinex,ielement,iprim)
  ElseIf((i.gt.1).and.(j.gt.1))then
     If(iha(i).eq.iha(j))then
      
       ME=GET_INTEGRALS_ASTRA(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)
     End IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     If(i.eq.j)then
        ivalue2=0.d0
        Do ia=1,Nbasis
           If(ielement(ia)%RCtype)then
              k=ielement(ia)%RCukLoc
              ivalue2=ivalue2+GET_INTEGRALS_ASTRA(ia,ia,0,0,flag,Bsplinex,ielement,iprim)
           End IF
        End Do
        ME=ME+ivalue2
     ENd IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     If(iep(i).eq.iep(j))then
       ivalue=GET_INTEGRALS_ASTRA(iha(j),iha(i),0,0,flag,Bsplinex,ielement,iprim)
       ME=ME-ivalue
     End IF
  End If

End function OneBody_STEX_ASTRA




!these are the one body matix elements for the STEX hamiltonian where the orbitals are in general not orthonormal
function OneBody_STEX_NO(i,j,flag,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleBspline
  use ukrmol_interface
  implicit none
  integer, intent(in) :: i,j,flag
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,ib
  integer k
  integer*8 i8,i8flag
  real(kind(1d0)) :: ivalue
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL
  integer, external :: kron
  ME=0.d0
  i8flag=int(flag,kind(i8flag))
  !id=0
  !This term remains the same as the reference configuration orbitals are orthonormal to the rest of states
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           ib=ielement(ia)%RCukLoc
           i8=ib
           ivalue=ivalue+GET_INTEGRALS(i8,i8,0,0,i8flag)
           !ivalue=ivalue+GET_1B_INTEGRAL_UKRMOL(ib,ib,flag,0,1)
           !write(*,*) ib,GET_1B_INTEGRAL_ASTRA(ib,ib,flag,0,1)
           !id=id+1
        End IF
     End Do
     ME=ivalue
     !write(*,*) "Sdas",ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     k=max(i,j)
        ME=GET_INTEGRALS_E(iha(k),iep(k),0,0,flag,Bsplinex,ielement,iprim)

  ElseIf((i.gt.1).and.(j.gt.1))then
     If(iha(i).eq.iha(j))then
        ME=GET_INTEGRALS_E(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)

        ivalue=0.d0
        If(abs(Spq(iep(j),iep(i))).gt.0.d0)Then
           Do ia=1,Nbasis
              If(ielement(ia)%RCtype)then
                 ib=ielement(ia)%RCukLoc
                 i8=ib
                 ivalue=ivalue+GET_INTEGRALS(i8,i8,0,0,i8flag)
                 !ivalue=ivalue+GET_1B_INTEGRAL_UKRMOL(ib,ib,flag,0,1)
                 !write(*,*) ib,GET_1B_INTEGRAL_ASTRA(ib,ib,flag,0,1),Spq(iep(j),iep(i))
              End IF
           End Do
        End If
        ME=ME+ivalue*Spq(iep(j),iep(i))
     ENd IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ME=ME-Spq(iep(j),iep(i))*GET_INTEGRALS_E(iha(j),iha(i),0,0,flag,Bsplinex,ielement,iprim)
  End If


End function OneBody_STEX_NO

function OneBody_STEX_TDMATRIX(i,j,flag,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleDensityMatricesJUAN
  use ModuleBspline
  use ukrmol_interface
  use ModuleParametersJUAN
  implicit none
  integer                , intent(in) :: i,j,flag
  type(ClassBspline)     , intent(in)             :: Bsplinex
  type(BasisElementInfo) , intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  type(ClassDensityMatrices) :: DM
  integer :: ir, is, ia, ip, ib, iq
  integer :: i1, i2
  real(kind(1d0)) :: ivalue,ivalue2
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL
  integer, external :: kron

  ME=0.d0
  if((i.eq.1).and.(j.eq.1))then
     do ia=1,Nelec
        ir=RefConf(ia)
        ME=ME+GET_INTEGRALS_ASTRA(ir,ir,0,0,flag,Bsplinex,ielement,iprim)
     enddo
  endif
  
  if(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     ia = iha(max(i,j))
     ip = iep(max(i,j))
     ivalue = 0.d0
     do i1=1,Nelec
        ir=RefConf(i1)
        do i2=1,Nelec
           is=RefConf(i2)
           
           ivalue2=DM%rho_stex(is,ia,ip,ir)
           if(abs(ivalue2).lt.0.00000001d0)cycle
           ME=ME-ivalue2*GET_INTEGRALS_ASTRA(ir,is,0,0,flag,Bsplinex,ielement,iprim)

           ! ME=ME-DM%rho_stex(is,ia,ip,ir)*GET_INTEGRALS_ASTRA(ir,is,0,0,flag,Bsplinex,ielement,iprim)
        enddo
     enddo
     ME=ME+GET_INTEGRALS_ASTRA(ip,ia,0,0,flag,Bsplinex,ielement,iprim)
  endif
   
  if((i.gt.1).and.(j.gt.1))then
     ivalue = 0.d0
     ia = iha(j)
     ip = iep(j)
     ib = iha(i)
     iq = iep(i)
     do i1=1,Nelec
        ir=RefConf(i1)
        do i2=1,Nelec
           is=RefConf(i2)
           ivalue2=DM%rho_stex(ib,ia,is,ir)
           if(abs(ivalue2).lt.0.00000001d0)cycle
           ivalue=ivalue+ivalue2*GET_INTEGRALS_ASTRA(ir,is,0,0,flag,Bsplinex,ielement,iprim)
           !ivalue=ivalue+DM%rho_stex(ib,ia,is,ir)*GET_INTEGRALS_ASTRA(ir,is,0,0,flag,Bsplinex,ielement,iprim)
        enddo
     enddo
     ME=ME+ivalue*Spq(ip,iq)
     if(ia.eq.ib)then
        ME=ME+GET_INTEGRALS_ASTRA(ip,iq,0,0,flag,Bsplinex,ielement,iprim)
     endif
  endif

  
 !     do ia = 1, Nbasis
 !        if(ielement(ia)%RCtype)then
 !           do ib = 1, Nbasis
 !              if(ielement(ib)%RCtype)then
 !                 ivalue=DM%rho_stex(iha(j),iha(i),ia,ib)
 !                 ivalue2=GET_INTEGRALS_ASTRA(ia,ib,0,0,flag,Bsplinex,ielement,iprim)
 !                 ME=ME+ivalue2*ivalue
 !              endif
 !           enddo
 !        endif
 !     enddo
 !     ME = ME * Spq(iep(j),iep(i))
 !     if(iha(i).eq.iha(j))then
 !        ME = ME + GET_INTEGRALS_ASTRA(iep(i),iep(j),0,0,flag,Bsplinex,ielement,iprim)
 !     endif


 !  !..  The conditional on j and i is to avoid the calculation of a null contribution
 !  ! if(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
 !     k=max(i,j)
 !     l=min(i,j)
 !     ivalue = 0.d0
 !     do ia=1,Nbasis
 !        if(ielement(ia)%RCtype)then
 !           if(ielement(iep(l))%RCtype)then
 !              ivalue=ivalue-DM%rho_stex(iha(k),iha(l),iep(l),ia)*GET_INTEGRALS_ASTRA(ia,iep(k),0,0,flag,Bsplinex,ielement,iprim)
 !           endif
 !        endif
 !     enddo
 !     ME=ME+ivalue
 !  ! endif

 !  !..  The conditional on j and i is to avoid the calculation of a null contribution
 ! ! if((i.eq.1).and.(j.eq.1))then
 !     ivalue = 0.d0
 !     Do ia=1,Nbasis
 !        If(ielement(ia)%RCtype)then
 !           Do ib=1,Nbasis
 !              If(ielement(ib)%RCtype)then
 !                 if((ielement(iep(i))%RCtype).and.(ielement(iep(j))%RCtype))then
 !                    ivalue2 = DM%pi_stex(iha(i),iha(j),ia,iep(i),iep(j),ib)
 !                    ivalue  = ivalue+ivalue2*GET_INTEGRALS_ASTRA(ia,ib,0,0,flag,Bsplinex,ielement,iprim)
 !                 endif
 !              End IF
 !           End Do
 !           if((ielement(iep(i))%RCtype).and.(ielement(iep(j))%RCtype))then
 !              ivalue2 = DM%rho_stex(iha(i),iha(j),ia,iep(j))
 !              ivalue  = ivalue-ivalue2*GET_INTEGRALS_ASTRA(ia,iep(i),0,0,flag,Bsplinex,ielement,iprim)
 !           endif
 !        End IF
 !     End Do
 !     ME=ME+ivalue
 ! ! endif
  
End function OneBody_STEX_TDMATRIX


function TwoBody_STEX_TDMATRIX(i,j,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleDensityMatricesJUAN
  use ModuleBspline
  use ukrmol_interface
  use ModuleParametersJUAN
  implicit none
  integer, intent(in) :: i,j
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  type(ClassDensityMatrices) :: DM
  integer :: i1, i2, i3, i4
  integer :: ir, is, it, iu, ib, iq, ia, ip
  real(kind(1d0)) :: ivalue,ivalue2
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  integer, external :: kron

  ME=0.d0
  if((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     do i1 = 1, Nelec
        ir = RefConf(i1)
        do i2 = 1, Nelec
           is = RefConf(i2)
           do i3 = 1, Nelec
              it = RefConf(i3)
              do i4 = 1, Nelec
                 iu = RefConf(i4)
                 ivalue2 = DM%rho_stex(is,ir,iu,it)
                 if(abs(ivalue2).lt.0.00000001d0)cycle
                 ivalue=ivalue+0.5d0*ivalue2*GET_INTEGRALS_ASTRA(ir,is,it,iu,0,Bsplinex,ielement,iprim)
              End Do
           End Do
        End Do
     End Do
     ME=ivalue
  endif

  if(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     ib=iha(max(i,j))
     iq=iep(max(i,j))
     ivalue = 0.d0
     do i1 = 1, Nelec
        ir = RefConf(i1)
        do i2 = 1, Nelec
           is = RefConf(i2)
           do i3 = 1, Nelec
              it = RefConf(i3)
              ivalue2 = DM%rho_stex(ib,ir,is,it)
              if(abs(ivalue2).lt.0.00000001d0)cycle
              ivalue = ivalue + 0.5d0*ivalue2*GET_INTEGRALS_ASTRA(ir,iq,it,is,0,Bsplinex,ielement,iprim)
              ivalue = ivalue - 0.5d0*ivalue2*GET_INTEGRALS_ASTRA(ir,is,it,iq,0,Bsplinex,ielement,iprim)
           enddo
        enddo
     enddo
     ME = ivalue
  endif

  if((i.gt.1).and.(j.gt.1))then
     ME=0.d0
     ia=iha(i)
     ip=iep(i)
     ib=iha(j)
     iq=iep(j)

     if(abs(Spq(ip,iq)).gt.1E-8)Then
        ivalue = 0.d0
        do i1 = 1, Nelec
           ir = RefConf(i1)
           do i2 = 1, Nelec
              is = RefConf(i2)
              do i3 = 1, Nelec
                 it = RefConf(i3)
                 do i4 = 1, Nelec
                    iu = RefConf(i4)
                    ivalue2 = DM%pi_stex(ib,ia,is,iu,ir,it)
                    if(abs(ivalue2).lt.0.000000000001d0)cycle
                    ivalue = ivalue + 0.5d0*ivalue2*GET_INTEGRALS_ASTRA(ir,is,it,iu,0,Bsplinex,ielement,iprim)
                 enddo
              enddo
           enddo
        enddo
        ME=ME+ivalue*Spq(ip,iq)
     endif

     ivalue = 0.d0
     do i1 = 1, Nelec
        ir = RefConf(i1)
        do i2 = 1, Nelec
           is = RefConf(i2)
           ivalue2 = DM%rho_stex(ib,ia,is,ir)
           if(abs(ivalue2).lt.0.000000000001d0)cycle
           ivalue = ivalue + ivalue2*GET_INTEGRALS_ASTRA(ir,is,ip,iq,0,Bsplinex,ielement,iprim)
           ivalue = ivalue - ivalue2*GET_INTEGRALS_ASTRA(ir,iq,ip,is,0,Bsplinex,ielement,iprim)
        enddo
     enddo
     ME = ME + ivalue
  endif

 contains


   real(kind(1d0)) function kdelta(i,j) result(res)
     integer, intent(in) :: i, j
     res = 1.d0
     if(i.ne.j) res = 0.d0
   end function kdelta

    
End function TwoBody_STEX_TDMATRIX

!these are the one body matix elements for the STEX hamiltonian where the orbitals are in general not orthonormal
function OneBody_STEX_NO_ASTRA(i,j,flag,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleBspline
  use ukrmol_interface
  implicit none
  integer, intent(in) :: i,j,flag
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,ib
  integer k
  integer*8 i8,i8flag
  real(kind(1d0)) :: ivalue
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  integer, external :: kron
  ME=0.d0
  i8flag=int(flag,kind(i8flag))
  !id=0
  !This term remains the same as the reference configuration orbitals are orthonormal to the rest of states
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           ib=ielement(ia)%RCukLoc
           i8=ib
           ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,ia,0,0,flag,Bsplinex,ielement,iprim)
           !ivalue=ivalue+GET_1B_INTEGRAL_UKRMOL(ib,ib,flag,0,1)
           !write(*,*) ib,GET_1B_INTEGRAL_ASTRA(ib,ib,flag,0,1)
           !id=id+1
        End IF
     End Do
     ME=ivalue
     !write(*,*) "Sdas",ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     k=max(i,j)

     ME=GET_INTEGRALS_ASTRA(iha(k),iep(k),0,0,flag,Bsplinex,ielement,iprim)
     ! write(*,*) k,iha(k),iep(k),ME
     ! ME=GET_INTEGRALS_E(iha(k),iep(k),0,0,flag,Bsplinex,ielement,iprim)
     ! write(*,*) k,iha(k),iep(k),ME,"ME@"
        
  ElseIf((i.gt.1).and.(j.gt.1))then
     !write(*,*) iha(i),iha(j),Spq(iep(j),iep(i))
     If(iha(i).eq.iha(j))then
        ME=GET_INTEGRALS_ASTRA(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)

        !write(*,*) iep(j),iep(i),ME
      ! ME=GET_INTEGRALS_E(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)
      ! write(*,*) ME,"ME@"
        
        ivalue=0.d0 
        If(abs(Spq(iep(j),iep(i))).gt.0.d0)Then
           Do ia=1,Nbasis
              If(ielement(ia)%RCtype)then
                 ib=ielement(ia)%RCukLoc
                 i8=ib
                 ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,ia,0,0,flag,Bsplinex,ielement,iprim)
                 !ivalue=ivalue+GET_1B_INTEGRAL_UKRMOL(ib,ib,flag,0,1)
                 !write(*,*) ib,GET_1B_INTEGRAL_ASTRA(ib,ib,flag,0,1),Spq(iep(j),iep(i))
              End IF
           End Do
        End If
        ME=ME+ivalue*Spq(iep(j),iep(i))
     ENd IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ME=ME-Spq(iep(j),iep(i))*GET_INTEGRALS_ASTRA(iha(j),iha(i),0,0,flag,Bsplinex,ielement,iprim)
     !ivalue=Spq(iep(j),iep(i))*GET_INTEGRALS_ASTRA(iha(j),iha(i),0,0,flag,Bsplinex,ielement,iprim)
     !write(*,*) ivalue
  End If


End function OneBody_STEX_NO_ASTRA

!these are the one body matix elements for the STEX hamiltonian where the orbitals are in general not orthonormal
function TwoBody_STEX_NO_ASTRA(i,j,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleBspline
  use ukrmol_interface
  implicit none
  integer, intent(in) :: i,j
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,ib,ic,id
  integer*8 :: i8,j8
  real(kind(1d0)) :: ivalue,ivalue2
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  integer, external :: kron

  
  ME=0.d0
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           ic=ielement(ia)%RCukLoc
           i8=ic
           Do ib=1,Nbasis
              If(ielement(ib)%RCtype)then
                 id=ielement(ib)%RCukLoc
                 j8=ic
                 If(ia.ne.ib)then
                    ivalue=ivalue+0.5d0*GET_INTEGRALS_ASTRA(ia,ia,ib,ib,0,Bsplinex,ielement,iprim)
                    ivalue=ivalue-0.5d0*GET_INTEGRALS_ASTRA(ia,ib,ia,ib,0,Bsplinex,ielement,iprim)
                    
                    !ivalue=ivalue+0.5d0*GET_2B_INTEGRAL_ASTRA(ic,ic,id,id)
                    !ivalue=ivalue-0.5d0*GET_2B_INTEGRAL_ASTRA(ic,id,ic,id)*kron(ielement(ia)%spin,ielement(ib)%spin)
                 EndIf
              End IF
           End Do
        End IF
     End Do
     ME=ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     ib=max(i,j)
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           If(ia.ne.iha(ib))then
              ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,ia,iha(ib),iep(ib),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_ASTRA(ia,iha(ib),ia,iep(ib),0,Bsplinex,ielement,iprim)
           End IF
        End If
     End Do
     ME=ivalue
  ElseIf((i.gt.1).and.(j.gt.1))Then
     !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
     ivalue=0.d0
     !write(*,*) "ia ia",iha(i),iha(j)
     If(iha(i).eq.iha(j))then
        ivalue2=0.d0

        Do ia=1,Nbasis
           If(ielement(ia)%RCtype.and.(ia.ne.iha(i)))then
              !ME=GET_INTEGRALS_ASTRA(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              !write(*,*) "ME1",ME
              ivalue2=ivalue2+GET_INTEGRALS_ASTRA(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              ivalue2=ivalue2-GET_INTEGRALS_ASTRA(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
              !ME=GET_INTEGRALS_ASTRA(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
              !write(*,*) "ME2",ME
           End If
        End DO
         !write(*,*) ivalue2
         !pause
        ivalue=ivalue2
        ivalue2=0.d0
        If(abs(Spq(iep(j),iep(i))).gt.0.d0)Then
           Do ia=1,Nbasis
              If(ielement(ia)%RCtype)then
                 Do ib=1,Nbasis
                    If(ielement(ib)%RCtype)then
                       If((ia.ne.ib).and.(ia.ne.iha(i)).and.(ib.ne.iha(i)))Then
                          ivalue2=ivalue2+GET_INTEGRALS_ASTRA(ia,ia,ib,ib,0,Bsplinex,ielement,iprim)
                          ivalue2=ivalue2-GET_INTEGRALS_ASTRA(ia,ib,ia,ib,0,Bsplinex,ielement,iprim)
                       End If
                    End If
                 End Do
              End If
           End Do
           !write(*,*) ivalue2,Spq(iep(j),iep(i))
           ivalue=ivalue+Spq(iep(j),iep(i))*0.5d0*ivalue2
        End IF

        !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
     Else
        ivalue=ivalue+GET_INTEGRALS_ASTRA(iha(i),iep(i),iha(j),iep(j),0,Bsplinex,ielement,iprim)
        ivalue=ivalue-GET_INTEGRALS_ASTRA(iha(i),iha(j),iep(i),iep(j),0,Bsplinex,ielement,iprim)
        ivalue2=0.d0
        If(abs(Spq(iep(i),iep(j))).gt.0.d0)Then
           Do ia=1,Nbasis
              If((ielement(ia)%RCtype).and.(ia.ne.iha(i)).and.(ia.ne.iha(j)))then
                 ivalue2=ivalue2+GET_INTEGRALS_ASTRA(ia,iha(i),ia,iha(j),0,Bsplinex,ielement,iprim)
                 ivalue2=ivalue2-GET_INTEGRALS_ASTRA(ia,ia,iha(i),iha(j),0,Bsplinex,ielement,iprim)
              End IF
           End Do
           ivalue=ivalue+Spq(iep(i),iep(j))*ivalue2
        End If
     End IF
     ME=ivalue
  End If


End function TwoBody_STEX_NO_ASTRA

function TwoBody_STEX(i,j,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleSystemUtils
  use ModuleBspline
  use ukrmol_interface
  use ModuleUKRmolInterface
  implicit none
  integer, intent(in) :: i,j
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,ib,ic,id
  integer :: k,l,icount
  integer*8 :: i8,j8
  real(kind(1d0)) :: ivalue,ivalue2
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  integer, external :: kron

  ME=0.d0

  icount=0
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           i8=ielement(ia)%RCukLoc
           ic=ielement(ia)%RCukLoc
           Do ib=1,Nbasis
              If(ielement(ib)%RCtype)then
                 j8=ielement(ib)%RCukLoc
                 id=ielement(ib)%RCukLoc
                 If(ia.ne.ib)then
                    !ivalue=ivalue+0.5d0*GET_INTEGRALS(i8,i8,j8,j8,0)
                    !ivalue=ivalue-0.5d0*GET_INTEGRALS(i8,j8,i8,j8,0)*kron(ielement(ia)%spin,ielement(ib)%spin)
!                    ivalue2 =           0.5d0*GET_2B_INTEGRAL_ASTRA(ic,ic,id,id)
!                    ivalue2 = ivalue2 - 0.5d0*GET_2B_INTEGRAL_ASTRA(ic,id,ic,id)*kron(ielement(ia)%spin,ielement(ib)%spin)
                    ivalue2 =           0.5d0!*  GET_INTEGRALS_ASTRA(ia,ib,ic,id,flag,Bsplinex,ielement,iprim)
                    ivalue2 = ivalue2 - 0.5d0*GET_2B_INTEGRAL_ASTRA(ic,id,ic,id)*kron(ielement(ia)%spin,ielement(ib)%spin)
                    ivalue=ivalue+ivalue2
                    !write(*,*) k,ielement(ia)%spin,l,ielement(ib)%spin
                 EndIf
              End IF
           End Do
        End IF
     End Do
     ME=ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then

     ib=max(i,j)
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           k=ielement(ia)%RCukLoc
           If(ia.ne.iha(ib))then
              ivalue=ivalue+GET_INTEGRALS_E(ia,ia,iha(ib),iep(ib),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_E(ia,iha(ib),ia,iep(ib),0,Bsplinex,ielement,iprim)
              !  pause
           End IF
        End If
     End Do
     ME=ivalue

  ElseIf((i.gt.1).and.(j.gt.1))Then
     ivalue=0.d0
     !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
     If((iha(i).ne.iha(j)).and.(iep(i).ne.iep(j)))then
        ivalue=ivalue+GET_INTEGRALS_E(iha(j),iep(j),iha(i),iep(i),0,Bsplinex,ielement,iprim)
        ivalue=ivalue-GET_INTEGRALS_E(iha(j),iha(i),iep(j),iep(i),0,Bsplinex,ielement,iprim)
     End IF
     !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
     If((iha(i).eq.iha(j)).and.(iep(i).ne.iep(j)))then
        Do ia=1,Nbasis
           If((ielement(ia)%RCtype).and.(ia.ne.iha(i)))then
              k=ielement(ia)%RCukLoc
              !      write(*,*) ia,iep(j),iep(i)
              ivalue=ivalue+GET_INTEGRALS_E(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_E(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
           End IF
        End Do
     End IF
     !33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
     If((iha(i).ne.iha(j)).and.(iep(i).eq.iep(j)))then
        ivalue=ivalue+GET_INTEGRALS_E(iha(i),iep(j),iha(j),iep(i),0,Bsplinex,ielement,iprim)
        ivalue=ivalue-GET_INTEGRALS_E(iha(i),iha(j),iep(j),iep(i),0,Bsplinex,ielement,iprim)
        Do ia=1,Nbasis
           If((ielement(ia)%RCtype).and.(ia.ne.iha(i)).and.(ia.ne.iha(j)))then
              k=ielement(ia)%RCukLoc
              ivalue=ivalue+GET_INTEGRALS_E(ia,iha(j),ia,iha(i),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_E(iha(j),iha(i),ia,ia,0,Bsplinex,ielement,iprim)
           End IF
        End Do
     End IF
     !444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
     If((iha(i).eq.iha(j)).and.(iep(i).eq.iep(j)))then
        Do ia=1,Nbasis
           If((ielement(ia)%RCtype).and.(ia.ne.iha(i)))then
              k=ielement(ia)%RCukLoc
              ivalue=ivalue+GET_INTEGRALS_E(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_E(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
              Do ib=1,Nbasis
                 If((ielement(ib)%RCtype).and.(ib.ne.iha(i)).and.(ib.ne.ia))then
                    l=ielement(ib)%RCukLoc
                    ivalue=ivalue+0.5d0*GET_INTEGRALS_E(ia,ia,ib,ib,0,Bsplinex,ielement,iprim)
                    ivalue=ivalue-0.5d0*GET_INTEGRALS_E(ia,ib,ia,ib,0,Bsplinex,ielement,iprim)
                 End If
              End Do
           End IF
        End Do
     End IF
     ME=ivalue
  End If

End function TwoBody_STEX



!these are the one body matix elements for the STEX hamiltonian where the orbitals are in general not orthonormal
function TwoBody_STEX_NO(i,j,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleBspline
  use ukrmol_interface
  implicit none
  integer, intent(in) :: i,j
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,ib,ic,id
  integer*8 :: i8,j8
  real(kind(1d0)) :: ivalue,ivalue2
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  integer, external :: kron

  
  ME=0.d0
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           ic=ielement(ia)%RCukLoc
           i8=ic
           Do ib=1,Nbasis
              If(ielement(ib)%RCtype)then
                 id=ielement(ib)%RCukLoc
                 j8=ic
                 If(ia.ne.ib)then
                    ivalue=ivalue+0.5d0*GET_INTEGRALS(i8,i8,j8,j8,0)
                    ivalue=ivalue-0.5d0*GET_INTEGRALS(i8,j8,i8,j8,0)*kron(ielement(ia)%spin,ielement(ib)%spin)
                    
                    !ivalue=ivalue+0.5d0*GET_2B_INTEGRAL_ASTRA(ic,ic,id,id)
                    !ivalue=ivalue-0.5d0*GET_2B_INTEGRAL_ASTRA(ic,id,ic,id)*kron(ielement(ia)%spin,ielement(ib)%spin)
                 EndIf
              End IF
           End Do
        End IF
     End Do
     ME=ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
     ib=max(i,j)
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           If(ia.ne.iha(ib))then
              ivalue=ivalue+GET_INTEGRALS_E(ia,ia,iha(ib),iep(ib),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_E(ia,iha(ib),ia,iep(ib),0,Bsplinex,ielement,iprim)
           End IF
        End If
     End Do
     ME=ivalue
  ElseIf((i.gt.1).and.(j.gt.1))Then
     !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
     ivalue=0.d0
     If(iha(i).eq.iha(j))then
        ivalue2=0.d0

        Do ia=1,Nbasis
           If(ielement(ia)%RCtype.and.(ia.ne.iha(i)))then
              ivalue2=ivalue2+GET_INTEGRALS_E(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              ivalue2=ivalue2-GET_INTEGRALS_E(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
           End If
        End DO
        ivalue=ivalue2
        ivalue2=0.d0
        If(abs(Spq(iep(j),iep(i))).gt.0.d0)Then
           Do ia=1,Nbasis
              If(ielement(ia)%RCtype)then
                 Do ib=1,Nbasis
                    If(ielement(ib)%RCtype)then
                       If((ia.ne.ib).and.(ia.ne.iha(i)).and.(ib.ne.iha(i)))Then
                          ivalue2=ivalue2+GET_INTEGRALS_E(ia,ia,ib,ib,0,Bsplinex,ielement,iprim)
                          ivalue2=ivalue2-GET_INTEGRALS_E(ia,ib,ia,ib,0,Bsplinex,ielement,iprim)
                       End If
                    End If
                 End Do
              End If
           End Do
           ivalue=ivalue+Spq(iep(j),iep(i))*0.5d0*ivalue2
        End IF

        !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
     Else
        ivalue=ivalue+GET_INTEGRALS_E(iha(i),iep(i),iha(j),iep(j),0,Bsplinex,ielement,iprim)
        ivalue=ivalue-GET_INTEGRALS_E(iha(i),iha(j),iep(i),iep(j),0,Bsplinex,ielement,iprim)
        ivalue2=0.d0
        If(abs(Spq(iep(i),iep(j))).gt.0.d0)Then
           Do ia=1,Nbasis
              If((ielement(ia)%RCtype).and.(ia.ne.iha(i)).and.(ia.ne.iha(j)))then
                 ivalue2=ivalue2+GET_INTEGRALS_E(ia,iha(i),ia,iha(j),0,Bsplinex,ielement,iprim)
                 ivalue2=ivalue2-GET_INTEGRALS_E(ia,ia,iha(i),iha(j),0,Bsplinex,ielement,iprim)
              End IF
           End Do
           ivalue=ivalue+Spq(iep(i),iep(j))*ivalue2
        End If
     End IF
     ME=ivalue
  End If


End function TwoBody_STEX_NO





function TwoBody_STEX_ASTRA(i,j,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN
  use ModuleSystemUtils
  use ModuleBspline
  use ukrmol_interface
  use ModuleUKRmolInterface
  implicit none
  integer, intent(in) :: i,j
  type(ClassBspline)    , intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
  integer :: ia,ib,ic,id
  integer :: k,l,icount
  integer*8 :: i8,j8
  real(kind(1d0)) :: ivalue,ivalue2
  real(kind(1d0)) :: ME
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  real(kind(1d0)), external  :: GET_INTEGRALS_ASTRA
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  integer, external :: kron

   ME=0.d0

  icount=0
  If((i.eq.1).and.(j.eq.1))then
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           i8=ielement(ia)%RCukLoc
           ic=ielement(ia)%RCukLoc
           Do ib=1,Nbasis
              If(ielement(ib)%RCtype)then
                 j8=ielement(ib)%RCukLoc
                 id=ielement(ib)%RCukLoc
                 If(ia.ne.ib)then
                    !ivalue=ivalue+0.5d0*GET_INTEGRALS(i8,i8,j8,j8,0)
                    !ivalue=ivalue-0.5d0*GET_INTEGRALS(i8,j8,i8,j8,0)*kron(ielement(ia)%spin,ielement(ib)%spin)
!                    ivalue2 =           0.5d0*GET_2B_INTEGRAL_ASTRA(ic,ic,id,id)
!                    ivalue2 = ivalue2 - 0.5d0*GET_2B_INTEGRAL_ASTRA(ic,id,ic,id)*kron(ielement(ia)%spin,ielement(ib)%spin)
                    ivalue2 =           0.5d0*GET_INTEGRALS_ASTRA(ia,ia,ib,ib,0,Bsplinex,ielement,iprim)
                    ivalue2 = ivalue2 - 0.5d0*GET_INTEGRALS_ASTRA(ia,ib,ia,ib,0,Bsplinex,ielement,iprim)
                    ivalue=ivalue+ivalue2
                    !write(*,*) k,ielement(ia)%spin,l,ielement(ib)%spin
                 EndIf
              End IF
           End Do
        End IF
     End Do
     ME=ivalue
  ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then

     ib=max(i,j)
     ivalue=0.d0
     Do ia=1,Nbasis
        If(ielement(ia)%RCtype)then
           k=ielement(ia)%RCukLoc
           If(ia.ne.iha(ib))then
              ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,ia,iha(ib),iep(ib),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_ASTRA(ia,iha(ib),ia,iep(ib),0,Bsplinex,ielement,iprim)
              !  pause
           End IF
        End If
     End Do
     ME=ivalue

  ElseIf((i.gt.1).and.(j.gt.1))Then
     ivalue=0.d0
     !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
     If((iha(i).ne.iha(j)).and.(iep(i).ne.iep(j)))then
        ivalue=ivalue+GET_INTEGRALS_ASTRA(iha(j),iep(j),iha(i),iep(i),0,Bsplinex,ielement,iprim)
        ivalue=ivalue-GET_INTEGRALS_ASTRA(iha(j),iha(i),iep(j),iep(i),0,Bsplinex,ielement,iprim)
     End IF
     !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
     If((iha(i).eq.iha(j)).and.(iep(i).ne.iep(j)))then
        Do ia=1,Nbasis
           If((ielement(ia)%RCtype).and.(ia.ne.iha(i)))then
                    write(*,*) ia,iep(j),iep(i)
              ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_ASTRA(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
           End IF
        End Do
     End IF
     !33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
     If((iha(i).ne.iha(j)).and.(iep(i).eq.iep(j)))then
        ivalue=ivalue+GET_INTEGRALS_ASTRA(iha(i),iep(j),iha(j),iep(i),0,Bsplinex,ielement,iprim)
        ivalue=ivalue-GET_INTEGRALS_ASTRA(iha(i),iha(j),iep(j),iep(i),0,Bsplinex,ielement,iprim)
        Do ia=1,Nbasis
           If((ielement(ia)%RCtype).and.(ia.ne.iha(i)).and.(ia.ne.iha(j)))then
              k=ielement(ia)%RCukLoc
              ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,iha(j),ia,iha(i),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_ASTRA(iha(j),iha(i),ia,ia,0,Bsplinex,ielement,iprim)
           End IF
        End Do
     End IF
     !444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
     If((iha(i).eq.iha(j)).and.(iep(i).eq.iep(j)))then
        Do ia=1,Nbasis
           If((ielement(ia)%RCtype).and.(ia.ne.iha(i)))then
              k=ielement(ia)%RCukLoc
              ivalue=ivalue+GET_INTEGRALS_ASTRA(ia,ia,iep(j),iep(i),0,Bsplinex,ielement,iprim)
              ivalue=ivalue-GET_INTEGRALS_ASTRA(ia,iep(j),ia,iep(i),0,Bsplinex,ielement,iprim)
              Do ib=1,Nbasis
                 If((ielement(ib)%RCtype).and.(ib.ne.iha(i)).and.(ib.ne.ia))then
                    l=ielement(ib)%RCukLoc
                    ivalue=ivalue+0.5d0*GET_INTEGRALS_ASTRA(ia,ia,ib,ib,0,Bsplinex,ielement,iprim)
                    ivalue=ivalue-0.5d0*GET_INTEGRALS_ASTRA(ia,ib,ia,ib,0,Bsplinex,ielement,iprim)
                 End If
              End Do
           End IF
        End Do
     End IF
 
     ME=ivalue
  End If

End function TwoBody_STEX_ASTRA







