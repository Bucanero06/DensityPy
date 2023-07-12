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
  Module ParametersBsplines
  implicit none

  
  integer n_bs
  
end module ParametersBsplines

program ProgramTemplate

  use, intrinsic :: ISO_FORTRAN_ENV

  use ParametersBsplines
  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleBSpline
  use ModuleIO
  use ModuleConstants
  use ModuleDiagonalize

  implicit none

  !.. Run-time parameters
  !..
  ! integer                        :: nSize !
  ! character(len=:), allocatable  :: FileName
  
  ! !.. Local parameters
  ! real(kind(1d0)), allocatable   :: dMat(:,:),dOve(:,:)
  ! real(kind(1d0)), allocatable   :: dEval(:),Xi(:,:),Xip(:),Bo(:,:)
  ! integer                        :: i,j,k,l,m,n,ia,ib
  ! real(kind(1d0)), allocatable   :: WORK(:)
  ! integer                        :: LWORK,INFO,uid,Ndim
  ! integer, allocatable           :: iWORK(:),brr(:)
  

  ! type(ClassBspline)             :: Bsp_in,Bsp_out
  ! integer, parameter             :: Bsorder = 3
  ! integer                        :: Ng_in,Ng_out,Nbs_in,Nbs_out,Nextra,Nsub,Nl_Bs_in,Nf_Bs_out,NLAST
  ! real(kind(1d0)), allocatable   :: Gr_in(:),Gr_out(:),MaKe(:,:),MaOe(:,:)
  ! real(kind(1d0))                :: Ra,dx,r,dr,ivalue,ivalue2
  ! character*12                   :: GridFile
  ! procedure(D2DFun) , pointer    :: fptr,sptr

  integer                        :: i,j
  real(kind(1d0))                :: Ra,ivalue
  integer                        :: Ng_in,Ng_out
  real(kind(1d0)), allocatable   :: Gr_in(:),Gr_out(:)
  procedure(D2DFun) , pointer    :: fptr
  type(ClassBspline)             :: Bsp_in,Bsp_out
  integer, parameter             :: Bsorder = 3

  
  

  Ra=1.d0                  !Radius
  
  Ng_in = 10+1 
  write(*,*) "Ng_in",Ng_in
  write(*,*) "inner grid"
  allocate(Gr_in(1:Ng_in))
  Do i=1,Ng_in
     Gr_in(i)=Ra*dble(i-1)/dble(Ng_in-1)
     write(*,*) i, Gr_in(i)
  End Do
  call Bsp_in.Init(Ng_in,Bsorder,Gr_in)

  Ng_out = 5+1
  write(*,*) "N grid out",Ng_out
  write(*,*) "outer grid"
  allocate(Gr_out(1:Ng_out))
  Do i=1, Ng_out
     Gr_out(i)=0.5d0*Ra+0.5d0*Ra*dble(i-1)/dble(Ng_out-1)
     write(*,*) Gr_out(i)
  End Do
  call Bsp_out.Init(Ng_out,Bsorder,Gr_out)

  fptr => powers_r

  write(*,*) "Here it is the problem"
  
  ivalue=Bsp_out.Integral(fptr,1,1,0,0)
  write(*,*) "VALOR",ivalue
  pause
  ivalue=Bsp_in%Eval(1.0d0,1)
  ivalue=Bsp_out.Integral(fptr,1,1,0,0)
  write(*,*) "VALOR",ivalue


  
  stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  
!  Ng_in = 10+1             !Number of regular bsplines grid betwen 0.d0 and 1.d0

!  Nl_Bs_in=7               !This number is the last inner bspline used in the mixed expansion,
                           !the last with nonzero probability in (0.0,0.5) (automate this).
                           !we use Nl_Bs_in+1 to refer to the corresponding matrix element (first element is removed)
  
!  Nsub=0                   !Number of times we divide the firsts (Bsorder-1) intervals in order to increase the density of bsplines
!  Nextra=(Bsorder-1)*Nsub  !this is a multiple of (Bsorder-1), the Nsub being the number of times we divide each single interval
!  Ng_out = 5+1   !+Nextra      !Total number of external bsplines

!  Nbs_in=Ng_in+Bsorder-2
!  Nbs_out=Ng_out+Bsorder-2 !Ideal case





  
  !  n_bs=0
 ! The external grid has a bunch of grid points to match the basis
  ! Do i=1, Ng_out-Nextra
  !    Gr_out(i)=0.5d0*Ra+0.5d0*Ra*dble(i-1)/dble(Ng_out-Nextra-1)
  !    write(*,*) Gr_out(i)
  ! End Do
  ! dx=Gr_in(2)/dble(Nsub+1)
  ! j=Ng_out-Nextra+1
  ! Do i=1,(Bsorder-1)
  !    Do k=1,Nsub
  !    Gr_out(j)=Gr_out(i)+k*dx
  !    j=j+1
  !    End Do
  ! End Do
  ! allocate(brr(1:Ng_out))
  ! Do i=1,Ng_out
  !    brr(i)=i
  ! End Do
  ! call sort_0(Ng_out,Gr_out,brr)
  ! Do i=1, Ng_out
  !    write(*,*) i,Gr_out(i)
  ! End Do

  !Initialize the outher bsplines

  !write the inner basis
  ! GridFile="Bs_in000.dat"
  ! Do i=2,Nbs_in-1  !Internal is from 2 to 30, we plot also the 31
  !    write(GridFile(6:8),'(I3.3)') i
  !    open(newunit = uid     , &
  !         file    = GridFile)
  !    dr=0.001d0
  !    Do r=0.d0,Ra,dr
  !       ivalue=Bsp_in%Eval(r,i)
  !       If(abs(ivalue).gt.0.d0)then
  !          write(uid,*) r,ivalue
  !       End If
  !    End do
  !    close(uid)
  ! End do
  
  !write the outher basis
  ! GridFile="Bsout000.dat"
  ! Do i=Bsorder,Nbs_out-1
  !    write(GridFile(6:8),'(I3.3)') i
  !    open(newunit = uid     , &
  !         file    = GridFile)
  !    dr=0.001d0
  !    Do r=0.d0,Ra,dr
  !       ivalue=Bsp_out%Eval(r,i)
  !       If(abs(ivalue).gt.0.d0)then
  !          write(uid,*) r,ivalue
  !       End If
  !    End do
  !    close(uid)
  ! End do


  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  orthonormalization of the internal basis
  ! Ndim=Nl_Bs_in-1
  ! allocate(Bo(1:Ndim,1:Ndim))
  ! allocate(dMat(1:Ndim,1:Ndim),dOve(1:Ndim,1:Ndim),dEval(1:Ndim))
  ! dMat=0.d0
  ! dOve=0.d0
  ! dEval=0.d0
  ! Bo=0.d0
  ! Do i=1,Ndim
  !    Do j=1,Ndim
  !       dMat(i,j)=Bsp_in.Integral(fptr,j+1,i+1,0,0)
  !    End Do
  ! End Do
  ! dOve=dMat
  ! LWORK=Ndim*3-1
  ! allocate(WORK(1:LWORK))
  ! WORK=0.d0
  ! INFO=0
  ! call dsygv(1,'V','U',Ndim,dMat,Ndim,dOve,Ndim,dEval,WORK,LWORK,INFO)
  ! write(*,*) "INFO Bo in",INFO
  ! deallocate(WORK)
  ! Bo=dMat
  ! deallocate(dMat,dOve,dEval)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !Reference values for the kinetic energy, with the uniform basis in (0.0,1.0)
!   Ndim=Nbs_in-2 ! we remove the first and last elements
!   write(*,*) "Ndim in",Ndim
!   allocate(dMat(1:Ndim,1:Ndim),dOve(1:Ndim,1:Ndim),dEval(1:Ndim))
!   dMat=0.d0
!   dOve=0.d0
!   dEval=0.d0
!   Do i=1,Ndim
!      Do j=1,Ndim
!         dMat(j,i)=-0.5d0*Bsp_in.Integral(fptr,j+1,i+1,0,2)
!         dOve(j,i)=Bsp_in.Integral(fptr,j+1,i+1,0,0)
!      End Do
!   End Do
  
!  LWORK=Ndim*3-1
!  allocate(WORK(1:LWORK))
!  WORK=0.d0
!  INFO=0
!  call dsygv(1,'V','U',Ndim,dMat,Ndim,dOve,Ndim,dEval,WORK,LWORK,INFO)
!  write(*,*) "INFO eigen in",INFO
!  write(*,*) "Energy values vs analitical result"
!  Do i=1,Ndim
!     write(*,*) i,dEval(i),0.5*((i*3.14159265359d0)**2.d0)
!  End Do
! deallocate(WORK)

!  !We plot the eigenfunctions
!   GridFile="KE_e1000.dat"
!   Do j=1,5!Ndim
!      write(GridFile(6:8),'(I3.3)') j
!      open(newunit = uid     , &
!           file    = GridFile)
!      dr=0.01d0
!      Do r=0.d0,Ra,dr
!         ivalue=0.d0
!          Do i=1,Nbs_in-2
!             ivalue=ivalue+dMat(i,j)*Bsp_in%Eval(r,i+1)
!          End Do
!         write(uid,*) r,ivalue
!      End do
!      close(uid)
!      !  Integral = Bspline.Integral(fptr,i,Bs2,dBra,dKet)
!      !  write(*,*) "Integral",i,Integral
!   End do
!   deallocate(dMat,dOve,dEval)
!   pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   !The eigenvalues of the projector
!   Nf_Bs_out=Bsorder  !the first external Bspline considered, Tipically Bsorder
!   Ndim=Nbs_out-(Nf_Bs_out-1)-1 !we remove the first k-1 and the last one
!   write(*,*) "Ndim out",Ndim

!   allocate(dMat(1:Ndim,1:Ndim),dOve(1:Ndim,1:Ndim),dEval(1:Ndim))
!   dMat=0.d0
!   dOve=0.d0
!   dEval=0.d0
!   Do i=1,Ndim
! !     l=i+1 !position of the bspline
!      l=i+(Nf_Bs_out-1) !position of the bspline
!      write(*,*) i,l
!      Do j=i,Ndim
! !        m=j+1
!         m=j+(Nf_Bs_out-1)
!         Do k=Nl_Bs_in-(Bsorder-1)+1,Nl_Bs_in,1 !Do over the projector
!          !  write(*,*)i,j, k
!            ivalue=       MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,k,l,0,0,0.5d0,1.d0)
!            ivalue=ivalue*MixedBSplineIntegral(Bsp_out,Bsp_in,fptr,m,k,0,0,0.5d0,1.d0)
!            dMat(j,i)=dMat(j,i)+ivalue
!         End Do
!         !dOve(j,i)=MixedBSplineIntegral(Bsp_out,Bsp_out,fptr,m,l, 0,0,0.5d0,1.d0)  !it give same results than the line below
!         dOve(j,i)=Bsp_out.Integral(fptr,l,m,0,0)  !overlap
!         dOve(i,j)=dOve(j,i)                       
!         dMat(i,j)=dMat(j,i)
!      End Do
!   End Do
! !stop
!   write(*,*) "for the non orthonormalized internal basis"
!   allocate(Xi(1:Ndim,1:Ndim))  !the vectors with the coefficients. The Bsplines inside are not orthonormalized (we could do it)
!   Xi=0.d0
!   LWORK=Ndim*3-1
!   allocate(WORK(1:LWORK))
!   WORK=0.d0
!   INFO=0
!   call dsygv(1,'V','U',Ndim,dMat,Ndim,dOve,Ndim,dEval,WORK,LWORK,INFO)
!   write(*,*) "INFO eigen proyector 1",INFO,WORK(1)
!   Do i=1,Ndim!-(Bsorder-1)
!      write(*,*) i,dEval(i)
!      Xi(:,i)=dMat(:,i)
!   End Do

!   write(*,*) "for the orthonormalized internal basis"

!   allocate(MaKe(Nl_Bs_in-(Bsorder-1)+1:Nl_Bs_in,Nf_Bs_out:Nbs_out-1))
!   MaKe=0.d0
!   Do i=1,Ndim
!      n=i+(Nf_Bs_out-1)
!      Do k=Nl_Bs_in-(Bsorder-1)+1,Nl_Bs_in,1
!     !    write(*,*) i,k
!       MaKe(k,n)=MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,k,n,0,0,0.5d0,1.d0)
!    End Do
!    End Do
!   dMat=0.d0
!   dOve=0.d0
!   dEval=0.d0
!   Do i=1,Ndim
!      n=i+(Nf_Bs_out-1) !position of the bspline
!     !  write(*,*) i,n,Ndim
!      Do j=i,Ndim
!         m=j+(Nf_Bs_out-1)
!         Do ia=1,Nl_Bs_in-1 !Do over the projector
!            ib=ia
!          ! Do ib=1,Nl_Bs_in-1
!               Do k=Nl_Bs_in-(Bsorder-1)+1,Nl_Bs_in,1
!                  Do l=Nl_Bs_in-(Bsorder-1)+1,Nl_Bs_in,1
!                     dMat(j,i)=dMat(j,i)+MaKe(k,n)*MaKe(l,m)*Bo(k-1,ia)*Bo(l-1,ib)
!                  End Do
!               End Do
!           ! End Do
!         End Do

!         !dOve(j,i)=MixedBSplineIntegral(Bsp_out,Bsp_out,fptr,m,l, 0,0,0.5d0,1.d0)  !it give same results than the line below
!         dOve(j,i)=Bsp_out.Integral(fptr,n,m,0,0)  !overlap
!         dOve(i,j)=dOve(j,i)                       
!         dMat(i,j)=dMat(j,i)
!      End Do
!   End Do
! deallocate(MaKe)


!   WORK=0.d0
!   INFO=0
!   call dsygv(1,'V','U',Ndim,dMat,Ndim,dOve,Ndim,dEval,WORK,LWORK,INFO)
!   write(*,*) "INFO eigen proyector 2",INFO,WORK(1)
!   Do i=1,Ndim!-(Bsorder-1)
!      write(*,*) i,dEval(i)
!     Xi(:,i)=dMat(:,i)
!   End Do
  
!   deallocate(WORK)
!   deallocate(dMat,dOve,dEval)

!  pause


! write(*,*) "!Check orthonormality of the obtained states"
! Do j=1,Ndim
!    Do i=1,Ndim
!       ivalue=0.d0
!       Do ia=1,Ndim
!          l=ia+(Bsorder-1)
!          Do ib=1,Ndim
!             m=ib+(Bsorder-1)
!             ivalue=ivalue+Xi(ia,i)*Xi(ib,j)*Bsp_out.Integral(fptr,l,m,0,0)
!          End Do
!       End Do
!       write(*,*) i,j,ivalue
!    End Do
! End Do
!   pause
! write(*,*) "!Check orthogonality to the inner bsplines"
! ivalue=0.d0
! Do k=Nl_Bs_in-(Bsorder-1)+1,Nl_Bs_in,1
!    Do j=1,Ndim
!      ivalue=0.d0
!       Do ia=1,Ndim
!           l=ia+(Bsorder-1)
!           ivalue=ivalue+MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,k,l, 0,0,0.d0,1.d0)*Xi(ia,j)
!        End Do
!     write(*,*) k,j,ivalue
!     End Do
!  End Do
!  pause
! write(*,*) " !Check orthogonality to the inner orthonormalized bsplines"
! Do j=1,Ndim
!    Do ib=1,Nl_Bs_in-1
!       ivalue=0.d0
!       !      Do k=Nl_Bs_in-(Bsorder-1)+1,Nl_Bs_in,1
!        Do k=2,Nl_Bs_in,1
!       Do ia=1,Ndim
!          l=ia+(Bsorder-1)
         
!          ivalue=ivalue+MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,k,l, 0,0,0.d0,1.d0)*Xi(ia,j)*Bo(k-1,ib)
!        !  write(*,*) ia,k
!        End Do
!        End Do
!        write(*,*) j,ib,ivalue
!       ! pause
!     End Do
!  End Do
!  pause
! !! Plot the orthonormalized external basis
!   GridFile="Xi_e1000.dat"
!   Do j=Ndim,1,-1
!      write(GridFile(6:8),'(I3.3)') j
!      open(newunit = uid     , &
!           file    = GridFile)
!      dr=0.001d0
!      Do r=0.d0,Ra,dr
!         ivalue=0.d0
!          Do i=1,Ndim
!             ivalue=ivalue+Xi(i,j)*Bsp_out%Eval(r,i+Nf_Bs_out-1)
!          End Do
!         write(uid,*) r,ivalue
!      End do
!      close(uid)
!   End do
!  pause


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Lets try to represent each of the original bsplines between 0.0 and 1.0
  !with the combination of the two basis elements
  
 !  NLAST=Nbs_out-1-(Nf_Bs_out-1)
 !  NLAST=NLAST-0  !how many do we remove
 !  Ndim=Nl_Bs_in-1+NLAST 
 !  write(*,*) "Ndim",Ndim,Nl_Bs_in,NLAST
 !  allocate(Xip(1:Ndim))
 !     ivalue=Bsp_in%Eval(1.0d0,2)
 !     ivalue=Bsp_out.Integral(fptr,1,1,0,0)
 !     write(*,*) "VALOR",ivalue
 !     pause

 !  GridFile="RB_in000.dat"
 !  Do i=2,Nbs_in-1  !Do in bsplines
 !      Xip=0.d0
 !     Do j=1,Nl_Bs_in-1
 !        Do ia=1,Nl_Bs_in-1
 !    !       Xip(j)=Xip(j)+Bo(ia,j)*Bsp_in.Integral(fptr,ia+1,i,0,0)
 !        End Do
 !     ENd Do
 !     Do ia=1,NLAST
 !        j=Nl_Bs_in-1+ia
 !        Do l=Nf_Bs_out,Nbs_out-1  !do the in external bspline index
 !        !   Xip(j)=Xip(j)+MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,i,l,0,0,0.d0,1.d0)*Xi(l-Nf_Bs_out+1,ia)
 !        End Do
 !     End Do
 !     !pause
     
 !     ! Do l=1,Nbs_out  !do in external bspline index
 !     !    Do m=1,Nbs_out !do in external bspline index
 !     !       ivalue=Bsp_out.Integral(fptr,l,m,0,0)
 !     !       write(*,*) ivalue,l,m
 !     !    End Do
 !     ! End Do
 !     write(*,*) "-----------------------------"
 !     ivalue=Bsp_out.Integral(fptr,1,1,0,0)
 !     write(*,*) i,Nbs_in,ivalue
 ! !pause
 !     write(GridFile(6:8),'(I3.3)') i
 !     open(newunit = uid     , &
 !          file    = GridFile)
 !     dr=0.001d0
 !     Do r=0.000d0,Ra*1.0001d0,Ra*1.0001d0!dr
 !        ivalue=0.d0
 !        Do j=1,Nl_Bs_in-1
 !           Do ia=1,Nl_Bs_in-1
 !                                 ivalue=Bsp_out.Integral(fptr,1,1,0,0)
 !     write(*,*) ivalue
 !              write(*,*) ia,j,size(Xip),ia+1,r
 !              ivalue=ivalue+Xip(j)*Bo(ia,j)*Bsp_in%Eval(r,ia+1)
 !                   ivalue=Bsp_out.Integral(fptr,1,1,0,0)
 !     write(*,*) ivalue
 !           End Do
 !        End Do
 !        write(*,*) "aaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
 !        pause
 !        Do ia=1,NLAST
 !           j=Nl_Bs_in-1+ia
 !           Do l=Nf_Bs_out,Nbs_out-1  !do the in external bspline index
 !              write(*,*) j,l-Nf_Bs_out+1,ia,l
 !              ivalue=ivalue+Xip(j)*Xi(l-Nf_Bs_out+1,ia)*Bsp_out%Eval(r,l)
 !                                 ivalue=Bsp_out.Integral(fptr,1,1,0,0)
 !     write(*,*) ivalue
 !           End Do
 !        ENd Do
 !        write(*,*) "bbbbbbbbbbbbbbbbb"
 !        pause
 !        If(abs(ivalue).gt.0.d0)then
 !           write(uid,*) r,ivalue
 !        End If
 !     End do
 !     close(uid)
 !  End Do

 !  !  Do l=1,Nbs_out  !do in external bspline index
 !  !    Do m=1,Nbs_out !do in external bspline index
 !  !      ivalue=Bsp_out.Integral(fptr,l,m,0,0)
 !  !       write(*,*) ivalue,l,m
 !  !    End Do
 !  ! End Do
 ! stop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Kinetic Energy with the extended basis

!   allocate(MaKe(Nf_Bs_out:Nbs_out-1,Nf_Bs_out:Nbs_out-1),MaOe(Nf_Bs_out:Nbs_out-1,Nf_Bs_out:Nbs_out-1))
!   MaKe=0.d0
!   MaOe=0.d0
!   Do l=Nf_Bs_out,Nbs_out-1  !do in external bspline index
!      Do m=Nf_Bs_out,Nbs_out-1  !do in external bspline index
!         MaKe(l,m)=Bsp_out.Integral(fptr,l,m,0,2)
!         MaOe(l,m)=Bsp_out.Integral(fptr,l,m,0,0)
!        ! write(*,*) MaOe(l,m),l,m
!      End Do
!   End Do

!   Do l=1,Nbs_out  !do in external bspline index
!      Do m=1,Nbs_out !do in external bspline index
!        ivalue=Bsp_out.Integral(fptr,l,m,0,0)
!         write(*,*) ivalue,l,m
!      End Do
!   End Do
  
!   NLAST=Nbs_out-1-(Nf_Bs_out-1)
!   NLAST=NLAST-0  !how many do we remove
!   Ndim=Nl_Bs_in-1+NLAST 

!   write(*,*) Nf_Bs_out,Nbs_out-1,Nlast,Ndim
!   pause

  
!   allocate(dMat(1:Ndim,1:Ndim),dOve(1:Ndim,1:Ndim),dEval(1:Ndim))
!   dMat=0.d0
!   dOve=0.d0
!   dEval=0.d0

  

!   Do j=1,Nl_Bs_in-1
!      Do i=1,Nl_Bs_in-1
        
!         dMat(j,i)=-0.5d0*Bsp_in.Integral(fptr,j+1,i+1,0,2)
!         dOve(j,i)=       Bsp_in.Integral(fptr,j+1,i+1,0,0)
!         write(*,*) j,i,"a",ivalue2
!      End Do
!      Do ia=1,NLAST
!         i=Nl_Bs_in-1+ia
                
!         ivalue=0.d0
!         ivalue2=0.d0
!         Do l=Nf_Bs_out,Nbs_out-1  !do the in external bspline index
!            ivalue  = ivalue  + MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,j+1,l,0,2,0.d0,1.d0)*Xi(l-Nf_Bs_out+1,ia)
!            ivalue2 = ivalue2 + MixedBSplineIntegral(Bsp_in,Bsp_out,fptr,j+1,l,0,0,0.d0,1.d0)*Xi(l-Nf_Bs_out+1,ia)
!         End Do
!         dMat(j,i)=-0.5d0*ivalue
!         dOve(j,i)=ivalue2
!         write(*,*) j,i,"b",ivalue2
!      End Do
!   End Do
!   Do ib=1,NLAST
!      j=Nl_Bs_in-1+ib
!      Do i=1,Nl_Bs_in-1
        
!         ivalue=0.d0
!         ivalue2=0.d0
!         Do m=Nf_Bs_out,Nbs_out-1  !do in external bspline index
!            ivalue  = ivalue  + MixedBSplineIntegral(Bsp_out,Bsp_in,fptr,m,i+1,0,2,0.d0,1.d0)*Xi(m-Nf_Bs_out+1,ib)
!            ivalue2 = ivalue2 + MixedBSplineIntegral(Bsp_out,Bsp_in,fptr,m,i+1,0,0,0.d0,1.d0)*Xi(m-Nf_Bs_out+1,ib)
!         End Do
!         dMat(j,i)=-0.5d0*ivalue
!         dOve(j,i)=ivalue2
!         write(*,*) j,i,"c",ivalue2
!      End Do
!      Do ia=1,NLAST
!         i=Nl_Bs_in-1+ia
        
!         ivalue=0.d0
!         ivalue2=0.d0
!         Do l=Nf_Bs_out,Nbs_out-1  !do in external bspline index
!            Do m=Nf_Bs_out,Nbs_out-1  !do in external bspline index
!            ivalue  = ivalue  + Xi(l-Nf_Bs_out+1,ia)*Xi(m-Nf_Bs_out+1,ib)*MaKe(m,l)
!            ivalue2 = ivalue2 + Xi(l-Nf_Bs_out+1,ia)*Xi(m-Nf_Bs_out+1,ib)*MaOe(m,l)
!            End Do
!         End Do
!         dMat(j,i)=-0.5d0*ivalue
!         dOve(j,i)=ivalue2
!         write(*,*) j,i,"d",ivalue2,ia,ib
!      End Do
!   End Do

!      write(*,*) MaOe
!   pause

!   ! write(*,*) dMat
!   ! pause
!   ! write(*,*) dOve
!   ! pause

  
!  LWORK=Ndim*3-1
!  allocate(WORK(1:LWORK))
!  WORK=0.d0
!  INFO=0
!  call dsygv(1,'V','U',Ndim,dMat,Ndim,dOve,Ndim,dEval,WORK,LWORK,INFO)
!  write(*,*) "INFO eigen proyector",INFO
!  Do i=1,10
!     write(*,*) i,dEval(i),0.5*(i*3.14159265359d0)**2
!  End Do
!  deallocate(WORK)


!  GridFile="Ke_e2000.dat"
!  Do j=1,5!Ndim
!     write(GridFile(6:8),'(I3.3)') j
!     open(newunit = uid     , &
!          file    = GridFile)
!     dr=0.01d0
!     Do r=0.d0,Ra,dr
!        ivalue=0.d0
!        Do i=1,Nl_Bs_in-1
!           ivalue=ivalue+dMat(i,j)*Bsp_in%Eval(r,i+1)
!        End Do
!        Do ia=1,NLAST
!           i=Nl_Bs_in-1+ia
!           Do m=Nf_Bs_out,Nbs_out-1  
!              ivalue=ivalue+dMat(i,j)*Xi(m-Bsorder+1,ia)*Bsp_out%Eval(r,m)
!           End Do
!        End Do
!        write(uid,*) r,ivalue
!     End do
!     close(uid)
!  End do


 
! deallocate(dMat,dOve,dEval)





      
!  call GetRunTimeParameters( FileName, nSize )

  ! allocate(dMat(nSize,nSize))
  ! call random_number(dMat)
  ! dMat = dMat + transpose(dMat)

  ! allocate(dEval(nSize))
  ! call Short_Diag( nSize, dMat, dEval )

  ! call SaveVector( FileName, dEval, "formatted" )

  stop

contains
  
  Pure fUNCTION powers_r(r) RESULT(y)
  !  use ParametersBsplines
      implicit none
    REAL(kind(1d0)), INTENT(in) :: r
    REAL(kind(1d0)) :: y
    y = 1.d0!r**(n_bs-0)
  END FUNCTION powers_r

    Pure fUNCTION sin_r(r) RESULT(y)
    use ParametersBsplines
      implicit none
    REAL(kind(1d0)), INTENT(in) :: r
    REAL(kind(1d0)) :: y
    y = 1.d0*sin(r*3.14159265359d0)
  END FUNCTION sin_r

  subroutine sort_0(n,arr,brr)
    implicit none
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(inout) :: brr(n)
    REAL(kind(1d0)), INTENT(inout) :: arr(n)
    REAL(kind(1d0)) :: a
    integer :: b,i,j
    !Sorts an array arr(1:n) into ascending numerical order, by straight insertion, while making
    !the corresponding rearrangement of the array brr(1:n).
    do j=2,n		!Pick out each element in turn.
       a=arr(j)
       b=brr(j)
       do i=j-1,1,-1	!Look for the place to insert it.
          if(arr(i).le.abs(a))goto 10
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
       enddo
       i=0
10     arr(i+1)=a	!Insert it.
       brr(i+1)=b
    enddo
    return
  END subroutine sort_0

  !> Reads the run time parameters specified in the command line.
  subroutine GetRunTimeParameters( FileName, nSize )
    !
    use ModuleErrorHandling
    use ModuleCommandLineParameterList
    use ModuleString
    !
    implicit none
    !
    character(len=:), allocatable, intent(out) :: FileName
    integer                      , intent(out) :: nSize
    !
    character(len=*), parameter :: PROGRAM_DESCRIPTION=&
         "Template Programs which diagonalizes a random matrix"
    type( ClassCommandLineParameterList ) :: List
    character(len=512) :: strnBuf

    call List.SetDescription(PROGRAM_DESCRIPTION)
    call List.Add( "--help" , "Print Command Usage" )
    call List.Add( "-o"     , "Output File" ,"eval", "optional" )
    call List.Add( "-n"     , "Matrix dim"  , 100  , "required" )

    call List.Parse()

    if(List.Present("--help"))then
       call List.PrintUsage()
       stop
    end if

    call List.Get( "-o",  strnBuf  )
    allocate(FileName,source=trim(adjustl(strnBuf)))

    call List.Get( "-n", nSize )

    call List.Free()

    write(OUTPUT_UNIT,"(a)"   ) "Run time parameters :"
    write(OUTPUT_UNIT,"(a)"   ) "Output File : "//FileName
    write(OUTPUT_UNIT,"(a,i0)") "Matrix size : ",nSize
    !
  end subroutine GetRunTimeParameters

end program ProgramTemplate

