! Module dealing with the EM-field - atom interaction radial integrals
!
! First  version (memint.f90) by Thomas Carette, March 2012
! Second version by Luca Argenti, July 2017.

!..
!  The mmltpol module deals with computing multipole related things.
!  After a LONG hunt in codes breaching into atsp2k, Luca and I managed
!  to identify ang(j1,j2)%c being set to the angular part of the reduced
!  matrix element.
!  
!  /Tor, Oct 2016
!..

module ModuleCSFDipole

  use ModuleIO
  use ModuleErrorHandling
  use ModuleString

  use mradint
  use mcspline
  use mmltpol

  implicit none

  !.. Orbital integrals
  type(radint)      , allocatable :: rint(:)
  type(radint)      , allocatable :: vint(:)

  !.. b-spline integrals
  complex(kind(1d0)), allocatable, private :: rr(:)
  complex(kind(1d0)), allocatable, private :: Br(:,:)
  complex(kind(1d0)), allocatable, private :: Bd(:,:)
  complex(kind(1d0)), allocatable, private :: Bc(:,:)

  private contl


contains


  subroutine ComputeAndSaveDipoleBlocks( Dir, SaveFormatted )
    !
    ! written by TC, March 2012
    ! updated by LA, July  2017
    !
    ! for the velocity gauge, only the Vinti integrals for el left
    ! larger than el right have been calculated. Since they are
    ! antisymmetric, there is an extra -1 phase if el1>el2.
    !
    !=== NAME CONVENTIONS  ===========================
    character(len=*), parameter :: BIN_SUFFIX = ".bin"
    character(len=*), parameter :: FMT_SUFFIX = ".txt"
    !=================================================
    
    character(len=*)  ,  intent(in) :: Dir
    logical           ,  intent(in) :: SaveFormatted

    !.. Local parameters
    complex(kind(1d0)), parameter   :: Z0  = (0.d0,0.d0)
    real   (kind(1d0)), parameter   :: EPS =  1.d-12

    integer                         :: i, el1, el2, j1, j2, dim1, dim2, m
    real   (kind(1d0))              :: sgn
    complex(kind(1d0)), allocatable :: blk_l(:,:), blk_v(:,:)
    character(len=:)  , allocatable :: FileName_l, FileName_v, pwclabel1, pwclabel2
    
    !.. The matrix element  between pwc states has the form
    !    < a l || Ov || b l' > = Cab < l | l' > + Dab < l || O || l' >
    !   Here we define the two coefficients Cab and Dab, which we will
    !   save for later use.
    !..
    complex(kind(1d0)) :: zCab, zDab
    integer            :: uid

    !.. Here    "1" stands for the inital state (ket)
    !   whereas "2" stnads for the final  state (bra)

    !.. LC - LC block
    !..
    allocate( blk_l( nref(2), nref(1) ) )
    allocate( blk_v( nref(2), nref(1) ) )
    !
    blk_l = Z0
    blk_v = Z0
    !
    do j1 = 1, nref( 1 )
       do j2 = 1, nref( 2 )

          if( abs( ang( J1, J2 )%c ) < EPS ) cycle

          el1 = ang ( J1, J2)%i( 1 )
          el2 = ang ( J1, J2)%i( 2 )
          i   = angi(el1,el2)

          sgn = 1.d0; if( el1 > el2 ) sgn = -1.d0

          blk_l( j2, j1 ) = ang( J1, J2 )%c * rint( i )%I(1,1)
          blk_v( j2, j1 ) = ang( J1, J2 )%c * vint( i )%I(1,1) * sgn

       enddo
    enddo
    !
    !.. Save LC - LC blocks
    !..
    pwclabel1 = AlphabeticNumber( 0, ncfg( 1 ) - nref( 1 ), "0" )
    pwclabel2 = AlphabeticNumber( 0, ncfg( 2 ) - nref( 2 ), "0" )
    allocate( FileName_l, source = Dir // "HL_" // pwclabel2 // "_" // pwclabel1 )
    allocate( FileName_v, source = Dir // "HV_" // pwclabel2 // "_" // pwclabel1 )
    !
    call SaveMatrix( FileName_l // BIN_SUFFIX, blk_l, "unformatted" )
    call SaveMatrix( FileName_v // BIN_SUFFIX, blk_v, "unformatted" )
    !
    if( SaveFormatted )then
       call SaveMatrix( FileName_l // FMT_SUFFIX, blk_l, "formatted" )
       call SaveMatrix( FileName_v // FMT_SUFFIX, blk_v, "formatted" )
    endif
    !
    deallocate( blk_l, blk_v )


    !.. PWC - LC block
    !..
    allocate( blk_l( ndim, nref(1) ) )
    allocate( blk_v( ndim, nref(1) ) )
    !
    do j2 = nref( 2 ) + 1, ncfg( 2 )

       dim2 = ndim - nspecl( contl( 2, j2 ) )
       !
       blk_l = Z0
       blk_v = Z0
       !
       do j1 = 1, nref( 1 ) 

          if( abs( ang( J1, J2 )%c ) < EPS ) cycle

          el1 = ang(   J1,  J2 )%i( 1 )
          el2 = ang(   J1,  J2 )%i( 2 )
          i   = angi( el1, el2 )

          sgn = 1.d0; if( el1 > el2 ) sgn = -1.d0

          blk_l( 1 : dim2, j1 ) = ang( J1, J2 )%c * rint(i)%I(1,1:dim2)
          blk_v( 1 : dim2, j1 ) = ang( J1, J2 )%c * vint(i)%I(1,1:dim2) * sgn

       enddo

       !.. Save PWC - LC block
       !..
       pwclabel1 = AlphabeticNumber( 0             , ncfg( 1 ) - nref( 1 ), "0" )
       pwclabel2 = AlphabeticNumber( j2 - nref( 2 ), ncfg( 2 ) - nref( 2 ), "0" )
       if( allocated( FileName_l ) ) deallocate( FileName_l )
       if( allocated( FileName_v ) ) deallocate( FileName_v )
       allocate( FileName_l, source = Dir // "HL_" // pwclabel2 // "_" // pwclabel1 )
       allocate( FileName_v, source = Dir // "HV_" // pwclabel2 // "_" // pwclabel1 )
       call SaveMatrix( FileName_l // BIN_SUFFIX, blk_l, dim2, nref(1), "unformatted" )
       call SaveMatrix( FileName_v // BIN_SUFFIX, blk_v, dim2, nref(1), "unformatted" )

       if( SaveFormatted )then
          call SaveMatrix( FileName_l // FMT_SUFFIX, blk_l, dim2, nref(1), "formatted" )
          call SaveMatrix( FileName_v // FMT_SUFFIX, blk_v, dim2, nref(1), "formatted" )
       endif

    enddo
    deallocate(blk_l,blk_v)


    !.. LC - PWC blocks
    !
    allocate( blk_l( nref(2), ndim ) )
    allocate( blk_v( nref(2), ndim ) )
    !
    do j1 = nref( 1 ) + 1, ncfg( 1 )

       dim1 = ndim - nspecl( contl( 1, j1 ) )
       !
       blk_l = Z0
       blk_v = Z0
       !
       do j2 = 1, nref(2)

          if( abs( ang( J1, J2 )%c ) < EPS ) cycle

          el1 = ang(  J1, J2)%i( 1 )
          el2 = ang(  J1, J2)%i( 2 )
          i   = angi(el1,el2)

          sgn = 1.d0; if( el1 > el2 ) sgn = -1.d0

          blk_l( j2, 1:dim1 ) = ang( J1, J2 )%c * rint(i)%I(1,1:dim1)
          blk_v( j2, 1:dim1 ) = ang( J1, J2 )%c * vint(i)%I(1,1:dim1) * sgn

       enddo

       !.. Save PWC - LC block
       !..
       pwclabel1 = AlphabeticNumber( j1 - nref( 1 ), ncfg( 1 ) - nref( 1 ), "0" )
       pwclabel2 = AlphabeticNumber( 0             , ncfg( 2 ) - nref( 2 ), "0" )
       if( allocated( FileName_l ) ) deallocate( FileName_l )
       if( allocated( FileName_v ) ) deallocate( FileName_v )
       allocate( FileName_l, source = Dir // "HL_" // pwclabel2 // "_" // pwclabel1 )
       allocate( FileName_v, source = Dir // "HV_" // pwclabel2 // "_" // pwclabel1 )
       call SaveMatrix( FileName_l // BIN_SUFFIX, blk_l, nref(2), dim1, "unformatted" )
       call SaveMatrix( FileName_v // BIN_SUFFIX, blk_v, nref(2), dim1, "unformatted" )

       if( SaveFormatted )then
          call SaveMatrix( FileName_l // FMT_SUFFIX, blk_l, nref(2), dim1, "formatted" )
          call SaveMatrix( FileName_v // FMT_SUFFIX, blk_v, nref(2), dim1, "formatted" )
       endif

    enddo
    deallocate(blk_l,blk_v)


    !.. PWC - PWC blocks
    !
    allocate( blk_l( ndim, ndim ) )
    allocate( blk_v( ndim, ndim ) )
    !
    do j1 = nref( 1 ) + 1, ncfg( 1 )
       !
       dim1 = ndim - nspecl( contl( 1, j1 ) )
       !
       do j2=nref(2)+1,ncfg(2)

          dim2 = ndim - nspecl( contl( 2, j2 ) )

          blk_l = Z0
          blk_v = Z0

          zCab = Z0
          zDab = Z0

          if( abs( ang( J1, J2 )%c ) > EPS ) then

             el1 = ang(J1,J2)%i(1)
             el2 = ang(J1,J2)%i(2)
             i   = angi( min( el1, el2 ), max( el1, el2 ) )

             sgn = 1.d0; if( el1 > el2 ) sgn = -1.d0

             if( rint(i)%t == 0 )then
                do m=1,dim2
                   blk_l(m,m) = ang(J1,J2)%c * rint(i)%I(1,1)
                   blk_v(m,m) = ang(J1,J2)%c * vint(i)%I(1,1) * sgn
                enddo
                zCab = ang(J1,J2)%c * sgn
             else
                if (el1<el2)then
                   blk_l(1:dim2,1:dim1) = ang(J1,J2)%c * transpose( rint(i)%I(1:dim1,1:dim2) )
                   blk_v(1:dim2,1:dim1) = ang(J1,J2)%c * transpose( vint(i)%I(1:dim1,1:dim2) )
                   zDab = ang(J1,J2)%c 
                else
                   blk_l(1:dim2,1:dim1) = ang(J1,J2)%c * rint(i)%I(1:dim2,1:dim1)
                   blk_v(1:dim2,1:dim1) = ang(J1,J2)%c * vint(i)%I(1:dim2,1:dim1) * sgn
                   zDab = ang(J1,J2)%c * sgn 
                endif
             endif

          endif

          !.. Save PWC - PWC block
          !..
          pwclabel1 = AlphabeticNumber( j1 - nref( 1 ), ncfg( 1 ) - nref( 1 ), "0" )
          pwclabel2 = AlphabeticNumber( j2 - nref( 2 ), ncfg( 2 ) - nref( 2 ), "0" )
          if( allocated( FileName_l ) ) deallocate( FileName_l )
          if( allocated( FileName_v ) ) deallocate( FileName_v )
          allocate( FileName_l, source = Dir // "HL_" // pwclabel2 // "_" // pwclabel1 )
          allocate( FileName_v, source = Dir // "HV_" // pwclabel2 // "_" // pwclabel1 )
          call SaveMatrix( FileName_l // BIN_SUFFIX, blk_l, dim2, dim1, "unformatted" )
          call SaveMatrix( FileName_v // BIN_SUFFIX, blk_v, dim2, dim1, "unformatted" )

          !.. Save the coefficient of the asymptotic one-body dipole and overlap 
          !   contributions to the many-body dipole transition
          !..
          open(newunit = uid,&
               file    = FileName_v//"_coef"//FMT_SUFFIX,&
               form    ="formatted", &
               status  ="unknown")
          write(uid,"(4(x,e24.16))") dble(zCab), aimag(zCab), dble(zDab), aimag(zDab)
          close(uid)
               

          if( SaveFormatted )then
             call SaveMatrix( FileName_l // FMT_SUFFIX, blk_l, dim2, dim1, "formatted" )
             call SaveMatrix( FileName_v // FMT_SUFFIX, blk_v, dim2, dim1, "formatted" )
          endif

       enddo
    enddo
    !
    deallocate(blk_l,blk_v)


  end subroutine ComputeAndSaveDipoleBlocks


  subroutine rade1(orbs) !,enorbs)
    IMPLICIT NONE
    complex(kind(1d0)), intent(in) :: orbs(ndim,ndim,0:lmax)
    !    complex(kind(1d0)),intent(in) :: enorbs(ndim,0:lmax)
    integer l

    ! calculate the radial integral in B-spline basis

    allocate( rr((nnod-1)*(kord+1)))
    call makerr(t,nknot,kord,rr)
    !
    ! Length
    !
    allocate(Br(ndim,ndim))
    call brint_ext(Br,nknot,kord,ilast,t,rr)
    allocate(rint(nrint))
    call compute_radint(orbs,Br,nrint,ossn,ossl,lmax,nwf,&
         nspec,nspecl,inptr(1:nrint,:),rint)
    !
    ! Velocity
    !
    allocate(Bd(ndim,ndim))
    allocate(Bc(ndim,ndim))
    call bderiv_ext(nknot,kord,ilast,t,rr,Bd,Bc)
    allocate(vint(nrint))
    call compute_vint(orbs) !,enorbs)

  end subroutine rade1

  subroutine compute_vint(orbs) !,enorbs)
    use mcspline
    IMPLICIT NONE

    complex(kind(1d0)), parameter :: Z1=(1.d0,0.d0)
    complex(kind(1d0)), parameter :: Z0=(0.d0,0.d0)
    complex(kind(1d0)), parameter :: Zi=(0.d0,1.d0)

    complex(kind(1d0)),intent(in) :: orbs(ndim,ndim,0:lmax)
    !    complex(kind(1d0)),intent(in) :: enorbs(ndim,0:lmax)

    complex(kind(1d0)) :: mat(2*kord-1,ndim,1:lmax)

    integer m,i,j,e1,e2,lm,l1,l2,p1,p2,p1_,p2_
    character*1 tr

    complex(kind(1d0)) c,ZDOTU
    complex(kind(1d0)) vec3(ndim)

    !
    !   1:    Vinti J(l-1,l  ) ~       <l-1 || p || l  > 
    !
    ! This routine is made so that it yields the same integrals as grad.f in
    ! ATSP2K-librad
    !

    do lm=1,lmax
       do j=1,ndim
          m=kord-j
          do i=max(1,j-kord+1),min(ndim,j+kord-1)
             mat(m+i,j,lm)= Bd(i,j)+dble(lm)*Bc(i,j)
          enddo
       enddo
    enddo

    ! integrals

    ! According to mmltpol, ossn(e1)<=ossn(e2)

    e1=0
    e2=0
    do i=1,nrint
       e1=inptr(i,1)
       e2=inptr(i,2)
       l1=ossl(e1)
       l2=ossl(e2)

       p1=ossn(e1)-l1
       if(p1<=nspecl(l1))then
          lm=l1
          tr='T'
          c=-Z1
          if(l2>l1)then
             lm=l2
             tr='N'
             c=Z1
          endif
          call ZGBMV(tr,ndim,ndim,kord-1,kord-1,c,&
               mat(:,:,lm), &
               2*kord-1,orbs(:,p1,l1),1,Z0,vec3,1)

          p2=ossn(e2)-l2
          if(p2<=nspecl(l2))then
             vint(i)%t=0

             allocate(vint(i)%I(1,1))
             vint(i)%I(1,1)=ZDOTU(ndim,vec3,1,orbs(:,p2,l2),1) 

          else
             vint(i)%t=1

             allocate(vint(i)%I(1,1:ndim-nspecl(l2)))

             ! TODO: BLAS that
             p2=nspecl(l2)
             do p2_=1,ndim-nspecl(l2)
                p2=p2+1
                vint(i)%I(1,p2_)=ZDOTU(ndim,vec3,1,orbs(:,p2,l2),1)

             enddo
          endif
       else
          vint(i)%t=2

          allocate(vint(i)%I(1:ndim-nspecl(l1),1:ndim-nspecl(l2)))

          ! The reordering is thrown to the garbage ! so be it
          ! TODO?: BLAS that (matrix band * matrix)

          lm=l1
          tr='T'
          c=-Z1
          if(l2>l1)then
             lm=l2
             tr='N'
             c=Z1
          endif

          p2=nspecl(l2)
          do p2_=1,ndim-nspecl(l2)
             p2=p2+1

             call ZGBMV(tr,ndim,ndim,kord-1,kord-1,c,&
                  mat(:,:,lm), &
                  2*kord-1,orbs(:,p2,l2),1,Z0,vec3,1)

             ! TODO: BLAS that
             p1=nspecl(l1)
             do p1_=1,ndim-nspecl(l1)
                p1=p1+1
                vint(i)%I(p1_,p2_)=ZDOTU(ndim,vec3,1,orbs(:,p1,l1),1) !/(enorbs(p1,l1)-enorbs(p2,l2) )
             enddo
          enddo
       endif

    enddo

  end subroutine compute_vint

  integer function contl(i,ii)
    IMPLICIT NONE
    integer,intent(in) :: i,ii
    integer,save :: n
    n=nwf
    do
       if(occ(i)%w(n,ii)>0) exit
       n=n-1
    enddo
    contl=ossl(n)
  end function contl

end module ModuleCSFDipole
