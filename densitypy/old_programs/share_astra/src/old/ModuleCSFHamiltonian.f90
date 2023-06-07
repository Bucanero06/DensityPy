!
!> This module deals with the definition, calculation, and I/O procedures
!! for the non-relativistic field-free Hamiltonian in a CSF CC basis
!!
!! The routines have been adapated from those written by Thomas Carette in 2012
!< 
module ModuleCSFHamiltonian

  use mdata
  use msint
  use m1eint
  use mcsfl

  IMPLICIT NONE

  character(len=*), private , parameter :: THIS_MODULE  = "ModuleCSFHamiltonian"
  real(kind(1d0)) , private , parameter :: H0_THRESHOLD = 1.d-20

contains


  subroutine ComputeAndStore_Hamiltonian(&
       Hdir                   , &
       Force_LC_LC_Calculation, &
       Force_PW_LC_Calculation, &
       Force_PW_PW_Calculation, &
       Eval_LC_LC             , &
       Eval_PW_LC             , &
       Eval_PW_PW             , &
       Save_Formatted           &
       )
    !
    use ModuleErrorHandling
    use ModuleString
    use ModuleIO
    !
    implicit none
    !
    !.. Local constant parameters
    character(len=*)  , parameter :: ROUTINE_NAME ="ComputeAndStore_Hamiltonian"
    character(len=*)  , parameter :: THIS_ROUTINE = THIS_MODULE // "_" // ROUTINE_NAME
    complex(kind(1d0)), parameter :: Z0 = (0.d0,0.d0)
    !
    !.. Assumes there is one and only one symmetry
    !   Here, BL (block) refers to the term (e.g., 1Se) 
    !   rather than to a block of matrix elements between
    !   two sets of configurations
    integer           , parameter :: BL =  1 
    logical           , parameter :: PRINT_UNSPARSE_INFO = .FALSE.
    logical           , parameter :: STORE_UNSPARSE_COEF = .TRUE.

    !.. subroutine parameters
    character(len=*), intent(in) :: Hdir
    logical         , intent(in) :: Force_LC_LC_Calculation
    logical         , intent(in) :: Force_PW_LC_Calculation
    logical         , intent(in) :: Force_PW_PW_Calculation
    logical         , intent(in) :: Eval_LC_LC
    logical         , intent(in) :: Eval_PW_LC
    logical         , intent(in) :: Eval_PW_PW
    logical         , intent(in) :: Save_Formatted
    !
    !.. Local variables
    integer, allocatable :: ndim_true_cfg(:)
    integer, allocatable :: lorb_cfg(:)
    !
    complex(kind(1d0)), allocatable :: H0_blk(:,:)
    !
    integer            :: i, j, n, ri, k, nri
    integer            :: iPWC, di, li, p1
    integer            :: jPWC, dj, lj, p2
    real   (kind(1d0)) :: ang
    complex(kind(1d0)) :: z
    !
    character(len=:), allocatable :: FileName, pwcName, ipwcName, jpwcName
    integer :: uid
    logical :: file_exists


    !.. Initialize the angular coefficients
    !..
    call unsparse( PRINT_UNSPARSE_INFO, STORE_UNSPARSE_COEF )


    !.. Save on file the number of states in H0 per configuration
    !..
    allocate( ndim_true_cfg( ncfg_bl( BL ) ) )
    allocate( lorb_cfg(      ncfg_bl( BL ) ) )
    lorb_cfg(      1 : nref( BL ) ) = -1
    ndim_true_cfg( 1 : nref( BL ) ) =  1
    do j = nref( BL ) + 1, ncfg_bl( BL )
       lorb_cfg(      j ) = contl( j, BL )
       ndim_true_cfg( j ) = ndim - nspecl( contl( j, BL ) )
    end do
    !
    fileName = "CSF_Basis_"//term_bl(BL)//".dat"
    call OpenFile( FileName, uid, "write", "formatted" )
    deallocate(fileName)
    write(uid,"(a,2(x,i0))") term_bl(BL), ncfg_bl(BL), nref(BL)
    write(uid,"(*(x,i0))"  ) ncsf(:,BL)
    write(uid,"(*(x,i0))"  ) lorb_cfg(nref(BL)+1:)
    write(uid,"(*(x,i0))"  ) ndim_true_cfg(nref(BL)+1:)
    close(uid)
    deallocate(lorb_cfg)
    deallocate(ndim_true_cfg)


    !..  LC - LC  block
    if(allocated(FileName)) deallocate(FileName)
    allocate( FileName, source = trim(Hdir) // "H0_LC_LC" )
    !
    INQUIRE( file = FileName, exist = file_exists )
    if( ( ( .not. file_exists ) .or. Force_LC_LC_Calculation ) .and. Eval_LC_LC )then

       allocate(H0_blk(nref(BL),nref(BL)))
       H0_blk=Z0
       !
       do i = 1, nref(BL) !.. Cycles over the lower-diagonal CFG pairs 
          do j = 1, i     !   
             do n = 1, tot_hang(bl)%h(i,j)%ni !.. Cycles over the GUGA contributions
                ang =  tot_hang(BL)%h(i,j)%int(n)%c
                ri  =  tot_hang(BL)%h(i,j)%int(n)%inptr
                if(ri<=intptr(3))then
                   nri=iord(1,ri)
                   k  =iord(2,ri)
                   H0_blk(i,j) = H0_blk(i,j) + ang * sintval(nri,k)%I(1,1,1,1)
                else
                   H0_blk(i,j) = H0_blk(i,j) + ang * Iint(ri)%I(1,1)
                endif
             enddo
             H0_blk(j,i)=H0_blk(i,j)
          enddo
       enddo

       !.. Save LC - LC block
       call SaveMatrix( FileName, H0_blk, "unformatted" )
       !
       if ( Save_Formatted ) then
          if(allocated(FileName)) deallocate(FileName)
          allocate(FileName,source=trim(Hdir)//"H0_LC_LC.txt")
          call SaveMatrix( FileName, H0_blk, "formatted" )
       endif
       !
       deallocate(H0_blk)

    endif



    !.. PWC - LC  blocks
    !
    if( Eval_PW_LC )then
       !
       allocate(H0_blk(ndim,nref(BL)))
       !
       do iPWC = nref(BL)+1, ncfg_bl(BL)  !.. Cycles over the CSF CC channels
          !
          li   = contl( iPWC, BL )        !.. Orbital angular momentum of added electron
          di   = ndim - nspecl( li )      !.. Total number of CSF CC states in the channel
          !
          pwcName = AlphabeticNumber( iPWC - nref(BL), ncfg_bl(BL) - nref(BL), "0" )
          if(allocated(FileName)) deallocate(FileName)
          allocate(FileName,source=trim(Hdir)//"H0_"//pwcName//"_LC")
          !
          INQUIRE( file = FileName, exist = file_exists )
          if( file_exists .and. .not. Force_PW_LC_Calculation ) cycle
          !
          H0_blk=Z0
          !
          do j = 1, nref( BL )            !.. Cycles over the CSF LC configurations
             do n = 1, tot_hang( BL )%h( iPWC, j )%ni !.. Cycles over GUGA contributions
                ang =  tot_hang( BL )%h( iPWC, j )%int(n)%c
                ri  =  tot_hang( BL )%h( iPWC, j )%int(n)%inptr
                if( ri <= intptr(3) )then !.. slater integral
                   nri = iord( 1, ri )
                   k   = iord( 2, ri )
                   if(sintval(nri,k)%t==0) call Assert("Internal Error in "//THIS_ROUTINE//" 001")
                   do p2 = 1, di          !.. Cycle over the channel radial index
                      H0_blk( p2, j ) = H0_blk( p2, j ) + ang * sintval(nri,k)%I(1,1,1,p2)
                   enddo
                else
                   if(Iint(ri)%t==0) call Assert("Internal Error in "//THIS_ROUTINE//" 002")
                   do p2 = 1, di          !.. Cycle over the channel radial index
                      H0_blk( p2, j ) = H0_blk( p2, j ) + ang * Iint(ri)%I(1,p2)
                   enddo
                endif
             enddo
          enddo

          !.. Save PWC - LC block
          call SaveMatrix( FileName, H0_blk, di, nref( BL ), "unformatted" )
          !
          if ( Save_Formatted ) then
             if(allocated(FileName)) deallocate(FileName)
             allocate(FileName,source=trim(Hdir)//"H0_"//pwcName//"_LC.txt")
             call SaveMatrix( FileName, H0_blk, di, nref( BL ), "formatted" )
          endif

       enddo
       !
       deallocate(H0_blk)
       !
    endif


    !.. PWC - PWC blocks
    !
    if( Eval_PW_PW )then
       !
       allocate( H0_blk( ndim, ndim ) )
       !
       do iPWC = nref(BL)+1, ncfg_bl(BL)  !.. Cycles over the CSF CC bra channels
          !
          li   = contl( iPWC, BL )        !.. Orbital angular momentum of added electron in bra PWC
          di   = ndim - nspecl( li )      !.. Total number of CSF CC states in the bra PWC
          ipwcName = AlphabeticNumber( iPWC - nref(BL), ncfg_bl(BL) - nref(BL), "0" )
          !
          do jPWC = nref(BL)+1, iPWC      !.. Cycles over the CSF CC ket channels (lower triangular pairs)
             !
             lj   = contl( jPWC, BL )     !.. Orbital angular momentum of added electron in ket PWC
             dj   = ndim - nspecl( lj )   !.. Total number of CSF CC states in the ket PWC
             jpwcName = AlphabeticNumber( jPWC - nref(BL), ncfg_bl(BL) - nref(BL), "0" )

             if(allocated(FileName)) deallocate(FileName)
             allocate(FileName,source=trim(Hdir)//"H0_"//ipwcName//"_"//jpwcName)
             INQUIRE( file = FileName, exist = file_exists )
             if( file_exists .and. .not. Force_PW_PW_Calculation ) cycle

             H0_blk = Z0 

             do n = 1, tot_hang( BL )%h( iPWC, jPWC )%ni !.. Cycles over GUGA contributions
                ang =  tot_hang( BL )%h( iPWC, jPWC )%int(n)%c
                ri  =  tot_hang( BL )%h( iPWC, jPWC )%int(n)%inptr

                if(ri<=intptr(3))then !.. slater integral
                   !
                   nri=iord(1,ri)
                   k  =iord(2,ri)
                   !
                   if( sintval(nri,k)%t == 0 )then !.. Diagonal part
                      !
                      if( li /= lj ) call Assert("Internal Error in "//THIS_ROUTINE//" 004")
                      z = ang * sintval( nri, k )%I(1,1,1,1)
                      do p2 = 1, di
                         H0_blk(p2,p2) = H0_blk(p2,p2) + z
                      enddo
                      !
                   elseif( sintval(nri,k)%t == 1 )then
                      !
                      if( li /= lj ) call Assert("Internal Error in "//THIS_ROUTINE//" 005")
                      if( iPWC == jPWC )then
                         do p2 = 1, di
                            H0_blk(p2,p2) = H0_blk(p2,p2) + ang * sintval(nri,k)%I(1,1,1,p2)
                         enddo
                      endif
                      !
                   elseif( sintval(nri,k)%t == 21 )then
                      !
                      if( lj == ossl(intord(4,nri,k)) .or. ( iPWC == jPWC ) )then
                         do p2 = 1, dj
                            do p1 = 1, di
                               H0_blk(p1,p2) = H0_blk(p1,p2) + ang * sintval(nri,k)%I(1,p1,1,p2)
                            enddo
                         enddo
                      else !.. Executed only when iPWC /= jPWC
                         do p2 = 1, dj
                            do p1 = 1, di
                               H0_blk(p1,p2) = H0_blk(p1,p2) + ang * sintval(nri,k)%I(1,p2,1,p1)
                            enddo
                         enddo
                      endif
                      !
                   else ! R(a,b,el,e'l')
                      !
                      if( lj == ossl(intord(4,nri,k)) .or. ( iPWC == jPWC ) )then
                         H0_blk(1:di,1:dj) = H0_blk(1:di,1:dj) + ang * sintval(nri,k)%I(1,1,1:di,1:dj)
                      else !.. Executed only when iPWC /= jPWC
                         do p2 = 1, dj
                            H0_blk(1:di,p2) = H0_blk(1:di,p2) + ang * sintval(nri,k)%I(1,1,p2,1:di)
                         enddo
                      endif
                      !
                   endif
                   !
                else !.. single electron (?)
                   !
                   if(Iint(ri)%t==1) call Assert("Internal Error in "//THIS_ROUTINE//" 006")
                   if(Iint(ri)%t==0)then !.. diagonal terms, li=lj
                      if( li /= lj ) call Assert("Internal Error in "//THIS_ROUTINE//" 007")
                      z = ang * Iint(ri)%I(1,1)
                      do p2 = 1, dj
                         H0_blk( p2, p2 ) = H0_blk( p2, p2 ) + z
                      enddo
                   else !.. off-diagonal terms
                      H0_blk(1:di,1:dj) = H0_blk(1:di,1:dj) + ang*Iint(ri)%I(1:di,1:dj)
                   endif
                   !
                endif
                !
             enddo

             !.. Save PWC - PWC block
             call SaveMatrix( FileName, H0_blk, di, dj, "unformatted" )
             !
             if ( Save_Formatted ) then
                if(allocated(FileName)) deallocate(FileName)
                allocate(FileName,source=trim(Hdir)//"H0_"//ipwcName//"_"//jpwcName//".txt")
                call SaveMatrix( FileName, H0_blk, di, dj, "formatted" )
             endif

          enddo
       enddo
       !
       deallocate( H0_blk )
       !
    endif

    return
  end subroutine ComputeAndStore_Hamiltonian
  
  
end module ModuleCSFHamiltonian
