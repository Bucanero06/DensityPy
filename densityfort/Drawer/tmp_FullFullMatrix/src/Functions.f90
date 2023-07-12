  !> Reads the run time parameters specified in the command line.
  ! Flags: flag=0  one body hamiltonian
  !        flag=1  overlap
 ! function OneBody_STEX_DM(i,j,flag,Bsplinex,ielement,iprim,DM) result(ME)
 !   use ModuleBasisJUAN
 !   use ModuleBspline
 !   use ukrmol_interface 
 !   use ModuleDensityMatricesJUAN
 !   implicit none
 !   class(ClassDensityMatrices), intent(inout) :: DM
 !   integer,                intent(in) :: i,j,flag
 !   type(ClassBspline)    , intent(in)             :: Bsplinex
 !   type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
    
 !   integer :: ia,ib,ic,id,l,m,lp,mp,lk,mk,Bra,Ket
 !   integer :: k,symI,symJ
 !   integer*8 i8,j8,i8flag
 !   real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
 !   real(kind(1d0)) :: ME,nor
 !   real(kind(1d0)), external  :: GET_INTEGRALS_E
 !   integer, external :: kron
!

 !   ME=0.d0
 !   i8flag=int(flag,kind(i8flag))
 !   !id=0
    !This term remains the same as the reference configuration orbitals are orthonormal to the rest of states

  !   Do symI=1,nIrreps
  !      symJ=irrmult(DM.Orbitals(iha(i)).irr,DM.Orbitals(iha(j)).irr)
  !      allocate(DM.RhoDM(1:DM.irr_sizes(symI),1:DM.irr_sizes(symJ)))
  !     ! call ReadSTEXDensityMatrixRho(DM, "rho_dir", iha(i), iha(j), symI)
  !      symJ=irrmult(symI,symJ)
  !      ic=0
  !      Do ia=DM.imin_irr(symI),DM.imax_irr(symI)
  !         ic=ic+1
  !         id=0
  !         Do ib=DM.imin_irr(symJ),DM.imax_irr(symJ)
  !            id=id+1
  !            ME=ME+DM.RhoDM(ic,id)*GET_INTEGRALS_E(ic,id,0,0,flag,Bsplinex,ielement,iprim)
  !         End Do
  !      end do
  !      deallocate(DM.RhoDM)
  !   End Do


  !   ME=Spq(iep(j),iep(i))*ME
  !   If(iha(i).eq.iha(j))then
  !      ME=ME+GET_INTEGRALS_E(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)
  !   End IF
    

  ! End function OneBody_STEX_DM



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
    implicit none
    integer, intent(in) :: i,j,flag
    type(ClassBspline)    , intent(in)             :: Bsplinex
    type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
    integer :: ia,ib,ic,id,m,lp,mp,lk,mk,Bra,Ket
    integer*8 :: k,l,i8flag
    real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
    real(kind(1d0)) :: ME,nor
    real(kind(1d0)), external  :: GET_INTEGRALS_E
    ME=0.d0
    i8flag=int(flag,kind(i8flag))

    !id=0
    If((i.eq.1).and.(j.eq.1))then
       ivalue=0.d0
       Do ia=1,Nbasis
          If(ielement(ia)%RCtype)then
             k=int(ielement(ia)%RCukLoc,kind(k))
             ivalue=ivalue+GET_INTEGRALS(k,k,0,0,i8flag)
             id=id+1
          End IF
       End Do
      ! write(*,*) id,ivalue
      ! pause
       ME=ivalue
    ElseIf(((i.eq.1).and.(j.gt.1)).or.((j.eq.1).and.(i.gt.1)))then
       lk=max(i,j)
       
       ME=GET_INTEGRALS_E(iep(lk),iha(lk),0,0,i8flag,Bsplinex,ielement,iprim)
    ElseIf((i.gt.1).and.(j.gt.1))then
       If(iha(i).eq.iha(j))then
          ivalue=GET_INTEGRALS_E(iep(j),iep(i),0,0,flag,Bsplinex,ielement,iprim)
          ME=ivalue
       End IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       If(i.eq.j)then
          ivalue2=0.d0
          Do ia=1,Nbasis
             If(ielement(ia)%RCtype)then
                k=int(ielement(ia)%RCukLoc,kind(k))
                ivalue2=ivalue2+GET_INTEGRALS(k,k,0,0,i8flag)
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

   function TwoBody_STEX(i,j,flag,Bsplinex,ielement,iprim) result(ME)
     use ModuleBasisJUAN
     use ModuleSystemUtils
     use ModuleBspline
     use ukrmol_interface 
     implicit none
     integer, intent(in) :: i,j,flag
     type(ClassBspline)    , intent(in)             :: Bsplinex
    type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
    integer :: ia,ib,ic,id,m,lp,mp,lk,mk,Bra,Ket
    integer :: k,l,icount
    integer*8 :: i8,j8
    real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
    real(kind(1d0)) :: ME,nor
    real(kind(1d0)), external  :: GET_INTEGRALS_E
    integer, external :: kron

    icount=0
    If((i.eq.1).and.(j.eq.1))then
       ivalue=0.d0
       Do ia=1,Nbasis
          If(ielement(ia)%RCtype)then
             i8=ielement(ia)%RCukLoc
             Do ib=1,Nbasis
                If(ielement(ib)%RCtype)then
                   j8=ielement(ib)%RCukLoc
                   If(ia.ne.ib)then
                      !write(*,*) ielement(ia)%RCukLoc,ielement(ia)%RCukLoc
                     ! write(*,*) ielement(ia)%RCukLoc,ielement(ib)%RCukLoc
                      ivalue2=0.5d0*GET_INTEGRALS(i8,i8,j8,j8,0)
                      ivalue2=ivalue2-0.5d0*GET_INTEGRALS(i8,j8,i8,j8,0)*kron(ielement(ia)%spin,ielement(ib)%spin)
                      ivalue=ivalue+ivalue2
                      !ivalue=ivalue+0.5d0*GET_INTEGRALS(i8,i8,j8,j8,0)
                      !ivalue=ivalue-0.5d0*GET_INTEGRALS(i8,j8,i8,j8,0)*kron(ielement(ia)%spin,ielement(ib)%spin)
                     ! write(*,*) ivalue,ivalue2,icount,kron(ielement(ia)%spin,ielement(ib)%spin),ielement(ia)%spin,ielement(ib)%spin
                   ! icount=icount+1
                      !write(*,*) i8,j8
                      !write(*,*) GET_INTEGRALS(i8,i8,j8,j8,0),GET_INTEGRALS(i8,j8,i8,j8,0),kron(ielement(ia)%spin,ielement(ib)%spin)
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
              !  write(*,*) GET_INTEGRALS_E(ia,ia,iha(ib),iep(ib),0,Bsplinex,ielement,iprim)
              !   write(*,*)GET_INTEGRALS_E(ia,iha(ib),ia,iep(ib),0,Bsplinex,ielement,iprim)
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
    function OneBody_STEX_NO(i,j,flag,Bsplinex,ielement,iprim) result(ME)
    use ModuleBasisJUAN
    use ModuleBspline
        use ukrmol_interface 

    implicit none
    integer, intent(in) :: i,j,flag
    type(ClassBspline)    , intent(in)             :: Bsplinex
    type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
    integer :: ia,ib,ic,id,l,m,lp,mp,lk,mk,Bra,Ket
    integer k
    integer*8 i8,j8,i8flag
    real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
    real(kind(1d0)) :: ME,nor
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
             i8=ielement(ia)%RCukLoc
             ivalue=ivalue+GET_INTEGRALS(i8,i8,0,0,i8flag)
             id=id+1
          End IF
       End Do
       ME=ivalue
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
                   i8=ielement(ia)%RCukLoc
                   ivalue=ivalue+GET_INTEGRALS(i8,i8,0,0,i8flag)
                End IF
             End Do
          End If
          ME=ME+ivalue*Spq(iep(j),iep(i))
       ENd IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ME=ME-Spq(iep(j),iep(i))*GET_INTEGRALS_E(iha(j),iha(i),0,0,flag,Bsplinex,ielement,iprim)
       End If


     End function OneBody_STEX_NO


         !these are the one body matix elements for the STEX hamiltonian where the orbitals are in general not orthonormal
  function TwoBody_STEX_NO(i,j,flag,Bsplinex,ielement,iprim) result(ME)
    use ModuleBasisJUAN
    use ModuleBspline
        use ukrmol_interface 

    implicit none
    integer, intent(in) :: i,j,flag
    type(ClassBspline)    , intent(in)             :: Bsplinex
    type(BasisElementInfo), intent(in)             :: ielement(1:Nbasis),iprim(1:NprimUkrmol)
    integer :: ia,ib,ic,id,m,lp,mp,lk,mk,Bra,Ket
    integer :: k,l
    integer*8 :: i8,j8,k8,l8,i8flag
    real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
    real(kind(1d0)) :: ME,nor
    real(kind(1d0)), external  :: GET_INTEGRALS_E
    integer, external :: kron
 
    i8flag=int(flag,kind(i8flag))
    ME=0.d0
    If((i.eq.1).and.(j.eq.1))then
       ivalue=0.d0
       Do ia=1,Nbasis
          If(ielement(ia)%RCtype)then
             i8=ielement(ia)%RCukLoc
             Do ib=1,Nbasis
                If(ielement(ib)%RCtype)then
                   j8=ielement(ib)%RCukLoc
                   If(ia.ne.ib)then
                      ivalue=ivalue+0.5d0*GET_INTEGRALS(i8,i8,j8,j8,0)
                      ivalue=ivalue-0.5d0*GET_INTEGRALS(i8,j8,i8,j8,0)*kron(ielement(ia)%spin,ielement(ib)%spin)
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



   !position_flag=-4   ===>  positron Hamiltonian
   !position_flag=-2   ===>  Kinetic Energy
   !position_flag=-1   ===>  Nuclear attraction
   !position_flag=0    ===>  electron Hamiltonian
   !positron_flag>0    ===>  property integrals
   !last index Get_integrals   l  m
   !position_flag=1            0  0    Since we are dealing with solid spherica it is also equal to the overlap. 
   !position_flag=2            1 -1
   !position_flag=3            1  0
   !position_flag=4            1  1
   !position_flag=5            2 -2
   !position_flag=6            2 -1
   !position_flag=7            2  0
   !position_flag=8            2  1
   !position_flag=9            2  2    and so on up to l=20 which is set at line 265 of scatci_integrals.f90
  function GET_INTEGRALS_Xi(ia,ib,ic,id,flag,Bsplinex,ielement,iprim) result(ME)
    use ModuleBasisJUAN
    use ModuleBspline
    implicit none

    !    REAL(kind(1d0)), INTENT(in) :: r
    integer, intent(in) :: ia,ib,ic,id,flag
    type(ClassBspline)    , intent(in)             :: Bsplinex
    !integer  :: BSplineGetOrder
    type(BasisElementInfo), intent(in)             :: iprim(1:NprimUkrmol),ielement(1:Nbasis)   !
    integer :: ja,jb,jc,jd,i,j,nxi,nxj,ixi,jxi,k,l,m,lp,mp,lk,mk,Bra,Ket
    real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
    real(kind(1d0)) :: ME,nor
    logical         :: twoBcase
    real(kind(1d0)), external  :: GET_INTEGRALS_E


  !   n_bs=(1)
  ! write(*,*) "bspline"
  !                                 summ=Bsplinex.Integral(fptr1,11,6,0,0)
  !                                 write(*,*) summ
    !                                 stop
    If(ic*id.ne.0)then
       write(*,*) "GET_INTEGRALS_Xi:  yet not implemented"
       stop
    End IF
    
    ME=0.d0

    Do j=1,Nbasis
       Do i=1,Nbasis
          ME=ME+Xi(j,ia)*Xi(i,ib)*GET_INTEGRALS_E(j,i,ic,id,flag,Bsplinex,ielement,iprim)
       End DO
    End Do
    
       
  End function GET_INTEGRALS_Xi




  
   !position_flag=-4   ===>  positron Hamiltonian
   !position_flag=-2   ===>  Kinetic Energy
   !position_flag=-1   ===>  Nuclear attraction
   !position_flag=0    ===>  electron Hamiltonian
   !positron_flag>0    ===>  property integrals
   !last index Get_integrals   l  m
   !position_flag=1            0  0    Since we are dealing with solid spherica it is also equal to the overlap. 
   !position_flag=2            1 -1
   !position_flag=3            1  0
   !position_flag=4            1  1
   !position_flag=5            2 -2
   !position_flag=6            2 -1
   !position_flag=7            2  0
   !position_flag=8            2  1
   !position_flag=9            2  2    and so on up to l=20 which is set at line 265 of scatci_integrals.f90
  function GET_INTEGRALS_E(ia,ib,ic,id,flag,Bsplinex,ielement,iprim) result(ME)


    use ModuleBasisJUAN   
    use ModuleParametersJUAN
    use ModuleBspline
    use ModuleSystemUtils
    use ModuleConstants
    use mpi_mod 
    use ukrmol_interface 
    use precisn 
    implicit none

    !    REAL(kind(1d0)), INTENT(in) :: r
    integer, intent(in) :: ia,ib,ic,id
    integer, intent(in) :: flag
    procedure(D2DFun) , pointer    :: fptr1,fptr2
    type(ClassBspline)    , intent(in)             :: Bsplinex
    !integer  :: BSplineGetOrder
    type(BasisElementInfo), intent(in)             :: iprim(1:NprimUkrmol),ielement(1:Nbasis)   !
    integer :: ja,jb,jc,jd,i,j,k,l,m,lp,mp,lk,mk,Bra,Ket
    integer*8 :: i8,j8,k8,l8,flagd
    real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
    real(kind(1d0)) :: ME,nor
    real(kind(1d0)), external :: ThreeXlmIntegral,Xlm
    real(kind(1d0))                    :: parvec(2) = 0.d0
    logical         :: twoBcase

!write(*,*) cf(:,8)
    
  !   n_bs=(1)
  ! write(*,*) "bspline"
  !                                 summ=Bsplinex.Integral(fptr1,11,6,0,0)
  !                                 write(*,*) summ
  !                                 stop
    flagd=flag
    ME=0.d0
    fptr1 => power_r
    fptr2 => partial_r
    
    
    If((ic.eq.0).and.(id.eq.0))Then
       j=ia
       i=ib
       ! do m=1,Nbasis
       !    write(*,*) m,ielement(m)%EBtype
       !    end do
       ! pause
       !    write(*,*) i,j,ielement(j)%EBtype,ielement(i)%EBtype,flag
       If(ielement(j)%spin.eq.ielement(i)%spin)Then
          If((.not.ielement(j)%EBtype).and.(.not.ielement(i)%EBtype))Then
             i8=ielement(j)%nr
             j8=ielement(i)%nr
             ME=GET_INTEGRALS(i8,j8,0,0,flagd)
            ! write(*,*) "aca1"
          End IF
          If(((.not.ielement(j)%EBtype).and.ielement(i)%EBtype).or.(ielement(j)%EBtype.and.(.not.ielement(i)%EBtype)))Then
             ivalue=0.d0
         !    write(*,*) "aca2",i,j,ielement(j)%EBtype,ielement(i)%EBtype
             Do k=sum(nMO)+1,NprimUkrmol
                If(ielement(j)%EBtype.and.(.not.ielement(i)%EBtype))Then
                   jc=j
                   Bra=ielement(jc)%nr
                   lp=ielement(jc)%l
                   mp=ielement(jc)%m
                   jd=i
                   Ket=iprim(k)%cr
                   l=iprim(k)%l
                   m=iprim(k)%m
                ElseIf((.not.ielement(j)%EBtype).and.ielement(i)%EBtype)Then
                   jc=i
                   Ket=ielement(jc)%nr
                   l=ielement(jc)%l
                   m=ielement(jc)%m
                   jd=j
                   Bra=iprim(k)%cr
                   lp=iprim(k)%l
                   mp=iprim(k)%m
                End IF
!write(*,*) "l,lp",l,lp
                If(abs(cf(k,ielement(jd)%nr)).gt.0.d0)then
                   nor=NB(Bra)*NB(Ket)
                   If(iprim(k)%cr.ge.BsplineOrder)then
                      If((lp.eq.l).and.(mp.eq.m))Then
                         If((flag.eq.0).or.(flag.eq.-2))then
                            !RADIAL KINETIC ENERGY
                            parvec(1)=0.d0
                            ivalue2=-Bsplinex.Integral(fptr1,Bra,Ket,0,2,parvec=parvec)                         
                            parvec(1)=-2.d0
                            !CENTRIFUGAL BARRIER
                            ivalue2=ivalue2+dble(l*(l+1))*Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)                      
                            ivalue=ivalue+cf(k,ielement(jd)%nr)*ivalue2*nor*0.5d0
                         End IF
                         If(flag.eq.1)then
                            !OVERLAP
                            parvec(1)=0.d0
                            ivalue2=Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)*nor

                            ivalue=ivalue+cf(k,ielement(jd)%nr)*ivalue2
                            !  write(*,*) ivalue,cf(k,ielement(jd)%nr),ivalue2,ielement(jd)%nr,k
                         End If
                      End IF
                      If(flag.eq.0)then
                         !CENTRAL POTENTIAL
                         ivalue2=0.d0
                         Do lk=abs(l-lp),l+lp
                            parvec(1)=Rn2
                            parvec(2)=lk*1.d0
                            suml=-Zn2*((4*PI)/dble(2*lk+1))*Bsplinex.Integral(fptr2,Bra,Ket,0,0,parvec=parvec)*nor
                            summ=0.d0
                            Do mk=-lk,lk
                               summ=summ+ThreeXlmIntegral(lp,mp,lk,mk,l,m)*(Xlm(0.d0,0.d0,lk,mk)+Xlm(PI,0.d0,lk,mk))
                            End Do
                            ivalue2=ivalue2+suml*summ
                         End Do
                         ivalue=ivalue+cf(k,ielement(jd)%nr)*ivalue2
                      End If
                      !! MULTIPOLE INTEGRALS    !Need to be tested
                      If(flag.gt.1)then                                               
                         lk= l_mp(flag)
                         mk= m_mp(flag)
                         parvec(1)=1.d0*lk                     
                         ivalue2=Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)*sqrt(4.D0*PI/dble(2*lk+1))*nor                
                         ivalue=ivalue+cf(k,ielement(jd)%nr)*ivalue2*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
                      End if
                   End IF
                End IF
             End do
             ME=ivalue
          End If
          ivalue=0.d0


          If(ielement(j)%EBtype.and.ielement(i)%EBtype)Then
           !  write(*,*) "aca4"
             jc=i
             jd=j
             nor=NB(ielement(jd)%nr)*NB(ielement(jc)%nr)
             lp=ielement(j)%l
             mp=ielement(j)%m
             l=ielement(i)%l
             m=ielement(i)%m
             If((lp.eq.l).and.(mp.eq.m))Then
                If((flag.eq.0).or.(flag.eq.-2))then
                   !RADIAL KINETIC ENERGY
                   parvec(1)=0.d0
                   ivalue2=-Bsplinex.Integral(fptr1,ielement(jd)%nr,ielement(jc)%nr,0,2,parvec=parvec)                               
                   parvec(1)=-2.d0
                   !CENTRIFUGAL BARRIER
                   ivalue2=ivalue2+dble(l*(l+1))*Bsplinex.Integral(fptr1,ielement(jd)%nr,ielement(jc)%nr,0,0,parvec=parvec)
                   ivalue=ivalue+ivalue2*0.5d0*nor
                End IF
                If(flag.eq.1)then
                   !OVERLAP
                   parvec(1)=0.d0
                   ivalue2=Bsplinex.Integral(fptr1,ielement(jd)%nr,ielement(jc)%nr,0,0,parvec=parvec)*nor                        
                   ivalue=ivalue+ivalue2
                End IF
             End IF

             If(flag.eq.0)then
                !CENTRAL POTENTIAL
                ivalue2=0.d0
                Do lk=abs(l-lp),l+lp                                                                                             
                   parvec(1)=Rn2
                   parvec(2)=lk*1.d0
                   suml=-Zn2*((4*PI)/(dble(2*lk+1)))*Bsplinex.Integral(fptr2,ielement(jd)%nr,ielement(jc)%nr,0,0,parvec=parvec)*nor
                   summ=0.d0
                   Do mk=-lk,lk
                      summ=summ+ThreeXlmIntegral(lp,mp,l,m,lk,mk)*(Xlm(0.d0,0.d0,lk,mk)+Xlm(PI,0.d0,lk,mk))
                   End Do
                   ivalue2=ivalue2+suml*summ
                End Do
                ivalue=ivalue+ivalue2
             End If
            !! MULTIPOLE INTEGRALS     !Need to be tested
             If(flag.gt.1)then                                               
                lk= l_mp(flag)
                mk= m_mp(flag)
                parvec(1)=1.d0*lk
                ivalue2=Bsplinex.Integral(fptr1,ielement(jd)%nr,ielement(jc)%nr,0,0,parvec=parvec)*sqrt(4.D0*PI/dble(2*lk+1))*nor
                ivalue=ivalue+ivalue2*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
             End if
             ME=ivalue
          End IF

       End IF
    End IF

    !TWO BODY INTEGRALS

       
    If((ic.ne.0).and.(id.ne.0))then



       !       !Estos valores son identicos
       !       !ivalue=GET_INTEGRAL(i1,i2,i4,i3,0)
       !       !ivalue=GET_INTEGRAL(i2,i1,i3,i4,0)
       !       !ivalue=GET_INTEGRAL(i2,i1,i4,i3,0)
       !       !ivalue=GET_INTEGRAL(i3,i4,i1,i2,0)
       !       !ivalue=GET_INTEGRAL(i3,i4,i2,i1,0)
       !       !ivalue=GET_INTEGRAL(i4,i3,i2,i1,0)
       !       !ivalue=GET_INTEGRAL(i4,i3,i1,i2,0)


       ja=min(ia,ib)
       jb=max(ia,ib)
       jc=min(ic,id)
       jd=max(ic,id)
      ! write(*,*) ielement(ja)%OType,ielement(jb)%OType,ielement(jc)%OType,ielement(jd)%OType
       If((ielement(ja)%spin.eq.ielement(jb)%spin).and.(ielement(jc)%spin.eq.ielement(jd)%spin))Then
       twoBcase=.false.
          If((.not.ielement(ja)%EBtype).and.(.not.ielement(jb)%EBtype))then
             If((.not.ielement(jc)%EBtype).and.(.not.ielement(jd)%EBtype))then
                i8=ielement(ja)%nr
                j8=ielement(jb)%nr
                k8=ielement(jc)%nr
                l8=ielement(jd)%nr
                ME=GET_INTEGRALS(i8,j8,k8,l8,0)
                twoBcase=.true.
             End If
          End If
          !          !Terms involved in <K|g|K_{B}^{T}>: does not have contribution for T in the external Bsplines
          !Positivos:
          !EBtype EBtype MOtype MOtype
          If(ielement(ja)%EBtype.and.ielement(jb)%EBtype)then
             If(ielement(jc)%MOtype.and.ielement(jd)%MOtype)then
                nor=NB(ielement(ja)%nr)*NB(ielement(jb)%nr)
                !puedo chequar aqui si hay overlap entre lso bsplines
                lp=ielement(ja)%l
                mp=ielement(ja)%m
                l=ielement(jb)%l
                m=ielement(jb)%m
                Bra=ielement(ja)%nr
                Ket=ielement(jb)%nr
                ME=0.d0
                i8=ielement(jc)%nr
                j8=ielement(jd)%nr
                Do lk=abs(l-lp),l+lp
                   parvec(1)=-1.d0*(lk+1)
                   summ=0.d0
                   Do mk=-lk,lk
                      k8=flag_ml(lk,mk)
                      summ=summ+GET_INTEGRALS(i8,j8,0,0,k8)*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
                   End Do
                   ME=ME+summ*sqrt((4*PI)/dble(2*lk+1))*Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)*nor
                End Do
                twoBcase=.true.
             End If
          End If
          !MOtype MOtype EBtype EBtype         
          If(ielement(ja)%MOtype.and.ielement(jb)%MOtype)then
             If(ielement(jc)%EBtype.and.ielement(jd)%EBtype)then
                nor=NB(ielement(jc)%nr)*NB(ielement(jc)%nr)
                !puedo chequar aqui si hay overlap entre lso bsplines
                lp=ielement(jc)%l
                mp=ielement(jc)%m
                l=ielement(jd)%l
                m=ielement(jd)%m
                Bra=ielement(jc)%nr
                Ket=ielement(jd)%nr
                ME=0.d0
                i8=ielement(ja)%nr
                j8=ielement(jb)%nr
                Do lk=abs(l-lp),l+lp
                   parvec(1)=-1.d0*(lk+1)
                   summ=0.d0
                   Do mk=-lk,lk
                      k8=flag_ml(lk,mk)
                      summ=summ+GET_INTEGRALS(i8,j8,0,0,k8)*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
                   End Do
                   ME=ME+summ*sqrt((4*PI)/dble(2*lk+1))*Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)*nor
                End Do
               twoBcase=.true.
             End If
          End If
          !OBtype EBtype MOtype MOtype
          If(ielement(ja)%OBtype.and.ielement(jb)%EBtype)then
             If(ielement(jb)%nr.lt.(NfirstBspl+BsplineOrder-1))Then
                If(ielement(jc)%MOtype.and.ielement(jd)%MOtype)then
                   ME=0.d0
                   Do k=sum(nMO)+1,NprimUkrmol
                      If(abs(cf(k,ielement(ja)%nr)).gt.0.d0)then
                         If(iprim(k)%cr.ge.BsplineOrder)then
                            Bra=ielement(jb)%nr
                            lp=ielement(jb)%l
                            mp=ielement(jb)%m
                            Ket=iprim(k)%cr
                            l=iprim(k)%l
                            m=iprim(k)%m
                            nor=NB(iprim(k)%cr)*NB(ielement(jb)%nr)
                            i8=ielement(jc)%nr
                            j8=ielement(jd)%nr
                            Do lk=abs(l-lp),l+lp
                               parvec(1)=-1.d0*(lk+1)
                               summ=0.d0
                               Do mk=-lk,lk
                                  k8=flag_ml(lk,mk)
                                  summ=summ+GET_INTEGRALS(i8,j8,0,0,k8)*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
                               End Do
                               ME=ME+summ*sqrt((4*PI)/dble(2*lk+1))*cf(k,ielement(ja)%nr)*Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)*nor
                            End Do
                         End IF
                      End IF
                   End Do
                End If
                twoBcase=.true.
             End If
          End If
          
          !MOtype MOtype OBtype EBtype         
          If(ielement(ja)%MOtype.and.ielement(jb)%MOtype)then
             If(ielement(jc)%OBtype.and.ielement(jd)%EBtype)then
                If(ielement(jd)%nr.lt.(NfirstBspl+BsplineOrder-1))Then
                   ! write(*,*) ielement(ja)%OType,ielement(jb)%OType,ielement(jc)%OType,ielement(jd)%OType
                   ME=0.d0
                   Do k=sum(nMO)+1,NprimUkrmol
                      If(abs(cf(k,ielement(jc)%nr)).gt.0.d0)then
                         If(iprim(k)%cr.ge.BsplineOrder)then
                            Bra=ielement(jd)%nr
                            lp=ielement(jd)%l
                            mp=ielement(jd)%m
                            Ket=iprim(k)%cr
                            l=iprim(k)%l
                            m=iprim(k)%m
                            nor=NB(iprim(k)%cr)*NB(ielement(jd)%nr)
                            i8=ielement(ja)%nr
                            j8=ielement(jb)%nr
                            Do lk=abs(l-lp),l+lp
                               parvec(1)=-1.d0*(lk+1)
                               summ=0.d0
                               Do mk=-lk,lk
                                  k8=flag_ml(lk,mk)
                                  summ=summ+GET_INTEGRALS(i8,j8,0,0,k8)*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
                               End Do
                                !  write(*,*) "aca3",lk,ielement(jc)%nr,Bra,Ket
                                ME=ME+summ*sqrt((4*PI)/dble(2*lk+1))*cf(k,ielement(jc)%nr)*Bsplinex.Integral(fptr1,Bra,Ket,0,0,parvec=parvec)*nor
                                ! write(*,*) "aca3"
                            End Do
                         End IF
                      End IF
                   End Do
                End If
                twoBcase=.true.
             End If
          End If
          
                     !  If(ielement(jd)%nr.lt.(NfirstBspl+BsplineOrder-1))Then 

          ! If(.not.twoBcase)then
          !    write(*,*) "this case does not correspond to a physical situation in STEX Hamiltonian "
          !    stop
          ! End IF


       End IF
    End If


    
  contains

  Pure fUNCTION partial_r( r , parvec ) RESULT(y)
      implicit none
    REAL(kind(1d0)), INTENT(in) :: r
    DoublePrecision, optional, intent(in) :: parvec(*)
    REAL(kind(1d0)) :: y
    y = (r**(-0))*(min(r,parvec(1))**parvec(2))/(max(r,parvec(1))**(parvec(2)+1.d0))
  END FUNCTION partial_r




  Pure DoublePrecision function power_r(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = x**parvec(1)
  end function power_r
  

End function GET_INTEGRALS_E


  



  Function Irr_lm(l,m) result(i)
      implicit none
    integer, intent(in) :: l,m ! input
    integer             :: i ! output
    i=0
    If((mod(l,2).eq.0).and.(mod(m,2).eq.0))then
       If(m.ge.0)then
          i=1 
       Else
          i=4 
       End IF
    ElseIf((mod(l,2).eq.0).and.(abs(mod(m,2)).eq.1))then
       If(m.ge.0)then
          i=6 
       Else
          i=7 
       End IF
    ElseIf((mod(l,2).eq.1).and.(mod(m,2).eq.0))then
       If(m.ge.0)then
          i=5 
       Else
          i=8 
       End IF
    ElseIf((mod(l,2).eq.1).and.(abs(mod(m,2)).eq.1))then
       If(m.ge.0)then
          i=2 
       Else
          i=3 
       End IF
    End IF
  end function irr_lm

! !This is the order found in ukrmol (same order than dalton):
! !#      Ag    Xlm       l even m even  N =  l/2+1      1  ->  1
! !#      B3u   Xlm       l odd  m odd   N = (l+1)/2     2  ->  8
! !#      B2u   Xl-m      l odd  m odd   N = (l+1)/2     3  ->  7
! !#      B1g   Xl-m      l even m even  N =  l/2        4  ->  2 
! !#      B1u   Xlm       l odd  m even  N = (l+1)/2     5  ->  6
! !#      B2g   Xlm       l even m odd   N =  l/2        6  ->  3
! !#      B3g   Xl-m      l even m odd   N =  l/2        7  ->  4
! !#      Au    Xl-m      l odd  m even  N = (l-1)/2     8  ->  5
! !This is the order found in CloseCouplingBasis file
! !#      Ag    Xlm       l even m even  N =  l/2+1
! !#      B1g   Xl-m      l even m even  N =  l/2
! !#      B2g   Xlm       l even m odd   N =  l/2
! !#      B3g   Xl-m      l even m odd   N =  l/2
! !#      Au    Xl-m      l odd  m even  N = (l-1)/2
! !#      B1u   Xlm       l odd  m even  N = (l+1)/2
! !#      B2u   Xl-m      l odd  m odd   N = (l+1)/2
! !#      B3u   Xlm       l odd  m odd   N = (l+1)/2

  !This is the integral of three spherical harmonics (none of them conjugate)
function ThreeYlmIntegral(l1,m1,l2,m2,l3,m3) result(inte)
  use ModuleAngularMomentum
      implicit none
    integer, intent(in) :: l1,m1,l2,m2,l3,m3 ! input
    real(kind(1d0))     :: inte ! output
    real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
    inte=0.d0
    If((m1+m2).eq.(-m3))then
   inte=((-1.d0)**(m3))*sqrt(((2*l1+1)*(2*l2+1))/((4*PI)*(2*l3+1)))
   inte=inte*ClebschGordanCoefficient(l1,l2,l3,m1,m2)*ClebschGordanCoefficient(l1,l2,l3,0,0)
   End IF
  end function ThreeYlmIntegral  

!This is the integral of Xl1m1  Xl2m2 Yl3m3 
  function TwoXlmOneYlmIntegral(l1,m1,l2,m2,l3,m3) result(inte)
      implicit none
    integer, intent(in) :: l1,m1,l2,m2,l3,m3 ! input
    real(kind(1d0))     :: inte,amp1,fas1,amp2,fas2 ! output
    real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
    real   (kind(1d0)), external :: ThreeYlmIntegral
    write(*,*) "TwoXlmOneYlmIntegral must be complex"
    stop
    inte=0.d0
    If(m1.eq.0)then
       amp1=1.d0
       fas1=0.d0
    ElseIf(m1.lt.0)then
       amp1=1.d0/sqrt(2.d0)
       fas1=amp1*((-1.d0)**abs(m1))
    Else
       amp1=1.d0/sqrt(2.d0)
       fas1=-amp1*((-1.d0)**abs(m1))
    End IF
    If(m2.eq.0)then
       amp2=1.d0
       fas2=0.d0
    ElseIf(m2.lt.0)then
       amp2=1.d0/sqrt(2.d0)
       fas2=amp2*((-1.d0)**abs(m2))
    Else
       amp2=1.d0/sqrt(2.d0)
       fas2=-amp2*((-1.d0)**abs(m2))
    End IF
    inte=0.d0
    inte=     amp1*amp2*ThreeYlmIntegral(l1,abs(m1),l2, abs(m2),l3,m3)
    inte=inte+amp1*fas2*ThreeYlmIntegral(l1,abs(m1),l2,-abs(m2),l3,m3)
    inte=inte+fas1*amp2*ThreeYlmIntegral(l1,-abs(m1),l2, abs(m2),l3,m3)
    inte=inte+fas1*fas2*ThreeYlmIntegral(l1,-abs(m1),l2,-abs(m2),l3,m3)    
  end function TwoXlmOneYlmIntegral


!This is the integral of Xl1m1  Xl2m2 Yl3m3 
  function ThreeXlmIntegral(l1,m1,l2,m2,l3,m3) result(inteR)
    implicit none
    integer, intent(in) :: l1,m1,l2,m2,l3,m3 ! input
    ! real(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 ! output
    real(kind(1d0))     :: inteR
    complex(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 ! output
    complex(kind(1d0)), parameter  :: ci=(0.d0,1.d0),c1=(1.d0,0.d0),c0=(0.d0,0.d0)
    real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
    real   (kind(1d0)), external :: ThreeYlmIntegral
    inte=c0
    If(m1.eq.0)then
       amp1=c1*1.d0
       fas1=c1*0.d0
    ElseIf(m1.lt.0)then
       amp1=-ci*1.d0/sqrt(2.d0)
       fas1=-amp1*((-1.d0)**abs(m1))
    Else
       amp1=c1*1.d0/sqrt(2.d0)
       fas1=amp1*((-1.d0)**abs(m1))
    End IF
    If(m2.eq.0)then
       amp2=c1*1.d0
       fas2=c1*0.d0
    ElseIf(m2.lt.0)then
       amp2=-ci*1.d0/sqrt(2.d0)
       fas2=-amp2*((-1.d0)**abs(m2))
    Else
       amp2=c1*1.d0/sqrt(2.d0)
       fas2=amp2*((-1.d0)**abs(m2))
    End IF
    If(m3.eq.0)then
       amp3=c1*1.d0
       fas3=c1*0.d0
    ElseIf(m3.lt.0)then
       amp3=-ci*1.d0/sqrt(2.d0)
       fas3=-amp3*((-1.d0)**abs(m3))
    Else
       amp3=c1*1.d0/sqrt(2.d0)
       fas3=amp3*((-1.d0)**abs(m3))
    End IF
    inte=c0
   ! m1=m2=m3=0
    inte=     amp1*amp2*amp3*ThreeYlmIntegral(l1, abs(m1),l2, abs(m2),l3, abs(m3))
    inte=inte+fas1*fas2*fas3*ThreeYlmIntegral(l1,-abs(m1),l2,-abs(m2),l3,-abs(m3)) !this is exactly 0 because fas1,2,3=0

!    If m1=0 this is multiplied by i^{-1} because either m2 or m3 are negative
    inte=inte+amp1*amp2*fas3*ThreeYlmIntegral(l1, abs(m1),l2, abs(m2),l3,-abs(m3)) 
    inte=inte+amp1*fas2*amp3*ThreeYlmIntegral(l1, abs(m1),l2,-abs(m2),l3, abs(m3))
  
    inte=inte+amp1*fas2*fas3*ThreeYlmIntegral(l1, abs(m1),l2,-abs(m2),l3,-abs(m3))
 
    inte=inte+fas1*amp2*amp3*ThreeYlmIntegral(l1,-abs(m1),l2, abs(m2),l3, abs(m3)) 
    inte=inte+fas1*amp2*fas3*ThreeYlmIntegral(l1,-abs(m1),l2, abs(m2),l3,-abs(m3))
 
    inte=inte+fas1*fas2*amp3*ThreeYlmIntegral(l1,-abs(m1),l2,-abs(m2),l3, abs(m3))

    inteR=real(inte)
    If(abs(aimag(inte)).gt.(10.d0**(-10)))then

       write(*,*) "matrix element is complex, something is wrong",inte
       write(*,*) l1,m1,l2,m2,l3,m3
       stop
    End If

    
  end function ThreeXlmIntegral

!This is the integral of Xl1m1  Xl2m2 Yl3m3 
  function Xlm(theta,phi,l,m) result(inteR)
    use ModuleAngularMomentum
    implicit none
    real   (kind(1d0)), intent(in) :: theta,phi
    integer, intent(in) :: l,m ! input
    ! real(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 
    real(kind(1d0))     :: inteR! output
    complex(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 
    complex(kind(1d0)), parameter  :: ci=(0.d0,1.d0),c1=(1.d0,0.d0),c0=(0.d0,0.d0)
    real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
    inte=c0
    If(m.eq.0)then
       amp1=c1*1.d0
       fas1=c1*0.d0
    ElseIf(m.lt.0)then
       amp1=-ci*1.d0/sqrt(2.d0)
       fas1=amp1*((-1.d0)**abs(m))
    Else
       amp1=c1*1.d0/sqrt(2.d0)
       fas1=-amp1*((-1.d0)**abs(m))
    End IF
   
    inte=c0
   ! m1=m2=m3=0
    inte=amp1*Ylm(theta,phi,l,abs(m))+fas1*Ylm(theta,phi,l,-abs(m))

    inteR=real(inte)
    If(abs(aimag(inte)).gt.0.d0)then

       write(*,*) "Xlm is complex, something is wrong"
       write(*,*) theta,phi,l,m
       stop
    End If

    
  end function Xlm
  
  function kron(i,j) result(k)
      implicit none
    integer, intent(in) :: i,j ! input
    integer    :: k ! output
    k=0
    If(i.eq.j)then
       k=1
    End IF
  end function kron  


