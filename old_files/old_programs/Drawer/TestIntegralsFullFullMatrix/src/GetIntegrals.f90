!position_flag=1     ===>  Overlap
!position_flag=2     ===>  Kinetic Energy
!position_flag=3     ===>  Nuclear attraction
!position_flag=4     ===>  electron Hamiltonian
function GET_INTEGRALS_ASTRA(ia,ib,ic,id,flag,Bsplinex,ielement,iprim) result(ME)
  use ModuleBasisJUAN   
  use ModuleParametersJUAN
  use ModuleBspline
  use ModuleSystemUtils
  use ModuleConstants
  use mpi_mod 
  use ukrmol_interface 
  use precisn
  use ModuleIntegrals
  use ModuleUKRmolInterface
  implicit none

  !    REAL(kind(1d0)), INTENT(in) :: r
  integer, intent(in) :: ia,ib,ic,id
  integer, intent(in) :: flag
  procedure(D2DFun) , pointer    :: fptr1,fptr2
  type(ClassBspline)    , intent(in)             :: Bsplinex
  !integer  :: BSplineGetOrder
  type(BasisElementInfo), intent(in)             :: iprim(1:NprimUkrmol),ielement(1:Nbasis)   !
  integer :: ja,jb,jc,jd,i,j,k,l,m,lp,mp,lk,mk,Bra,Ket, nr, nrp
  integer :: iflag_mp, i_sf
  integer*8 :: i8,j8,k8,l8,flagd
  real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
  real(kind(1d0)) :: ME,nor,rvalue, factor
  real(kind(1d0)), external :: ThreeXlmIntegral,Xlm
  real(kind(1d0))                    :: parvec(2) = 0.d0
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_EXTENDED
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_EXTENDED
  logical         :: twoBcase
  integer   :: iIrra, iIrrb, aia, aib
  character*3 :: cha,chb
  !iflag_mp=UKRMol_Get_Integral_Index('Property integrals')

  flagd=flag
  ME=0.d0
  !ONE-BODY INTEGRALS
  If((ic.eq.0).and.(id.eq.0))Then
     i=ielement(ia)%i_sf
     j=ielement(ib)%i_sf
     iIrra = absTOirr(i)
     iIrrb = absTOirr(j)
     aia   = absTOrel(i)
     aib   = absTOrel(j)
     cha   = absTOchar(i)
     chb   = absTOchar(j)
     !write(*,*) i,j,ia,ib
     call Order_Indexes_1B(i,j)
     If(ielement(ia)%spin.eq.ielement(ib)%spin)Then
!        ME=GET_1B_INTEGRAL_EXTENDED( i , j, flag, 0, 1, Bsplinex, iprim )
        ME=GET_1B_INTEGRAL_ASTRA( i , j, flag, 0, 1 )


        if((cha.eq."SRC").and.(chb.eq."PWC"))then   
           l   = absTOl(i)
           m   = absTOm(i)
           nr  = absTOnr(i)-NfirstBspl+1
           i_sf = ielement(j)%i_sf
           rvalue =0.d0
           Do k=sum(nMO)+1,NprimUkrmol
              if(abs(cf(k,i_sf))<= 1.0E-15 )cycle
              nrp = iprim(k)%cr
              if(abs(nr-nrp)>=BsplineOrder)cycle
              !write(*,*) nr,nrp
              lp  = iprim(k)%l
              mp  = iprim(k)%m

              do lk = abs (l - lp ), l + lp
                 factor = 0.d0
                 do mk = -lk, lk
                    factor = factor +  &
                         ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                         GlobalIntegral%NucMoments%GetMoments( lk, mk )
                 end do
                 rvalue = rvalue +  cf(k,i_sf)*factor * &
                      GlobalIntegral%BsBlocks%GetMultipoles(nr, nrp, lk)       
              end do

              if(lp==l.and.mp==m)then
                 rvalue = rvalue - cf(k,i_sf)*ThreeXlmIntegral( lp, mp, 0, 0, l, m ) * &
                      GlobalIntegral%NucMoments%GetMoments( 0, 0 ) * &
                      GlobalIntegral%BsBlocks%GetMultipoles(nr, nrp, 0)
              endif
           enddo
           ME=ME+rvalue
        endif


        
        if((cha.eq."PWC").and.(chb.eq."PWC"))then
           l   = absTOl(i)
           m   = absTOm(i)
           nr  = absTOnr(i)!-NfirstBspl+1
           lp  = absTOl(j)
           mp  = absTOm(j)
           nrp = absTOnr(j)!-NfirstBspl+1
           !write(*,*) nr,l,m,nrp,lp,mp
           rvalue = 0.d0
           do lk = abs (l - lp ), l + lp
              factor = 0.d0
              do mk = -lk, lk
                 factor = factor +  &
                      ThreeXlmIntegral( lp, mp, lk, mk, l, m ) * &
                      GlobalIntegral%NucMoments%GetMoments( lk, mk )
              end do
              rvalue = rvalue +  factor * &
                   GlobalIntegral%BsBlocks%GetMultipoles(nr, nrp, lk)       
           end do

           if(lp==l.and.mp==m)then
              rvalue = rvalue - ThreeXlmIntegral( lp, mp, 0, 0, l, m ) * &
                   GlobalIntegral%NucMoments%GetMoments( 0, 0 ) * &
                   GlobalIntegral%BsBlocks%GetMultipoles(nr, nrp, 0)
           endif

           ME=ME+rvalue
        endif
        
     End IF
  End IF

  !TWO-BODY INTEGRALS


  If((ic.ne.0).and.(id.ne.0))then
     i=ielement(ia)%i_sf
     j=ielement(ib)%i_sf
     k=ielement(ic)%i_sf
     l=ielement(id)%i_sf
     If((ielement(ia)%spin.eq.ielement(ib)%spin).and.(ielement(ic)%spin.eq.ielement(id)%spin))Then
       !ME= GET_2B_INTEGRAL_EXTENDED(i,j,k,l, Bsplinex, iprim)
       ME= GET_2B_INTEGRAL_ASTRA(i,j,k,l)
     End IF
  End If
end function GET_INTEGRALS_ASTRA

real(kind(1d0)) function GET_1B_INTEGRAL_ASTRA( ia, ib, Integral_ID, li, mi ) result( qres )
  use ModuleUKRmolInterface
  use ukrmol_interface
  use ModuleIntegrals
  use ModuleBasisJUAN
  use ModuleMatrix
  implicit none
  integer, intent(in) :: ia, ib, Integral_ID
  integer, intent(in) :: li, mi
!  integer, optional, intent(in) :: lll, mmm
  character(len=:), allocatable :: OpLabel
  integer*8 :: da8, db8, dflag, dl, dm
  integer   :: ic, iIrra, iIrrb, aia, aib
  integer   :: il, ir, ilm, ilmp, l, m, lp, mp, nr, nrp
  character*3 :: cha,chb
  type(ClassMatrix), pointer :: Mat
  dflag = int(Integral_ID,kind(dflag))

  !set OpLabel vs Integral_ID
  OpLabel = I1B_ID_LIST(Integral_ID)
  
  !  da8    = UKRmol_orbitals%Get_IntIndex(a)
  !  db8    = UKRmol_orbitals%Get_IntIndex(b)
   il = ia
   ir = ib
  ! if(absTOchar(il).eq.absTOchar(ir))then
  !    if(absTOirr(il).gt.absTOirr(ir))then
  !       il = ib
  !       ir = ia
  !    endif
  ! endif
  ! if((absTOchar(ir).eq."LOC").and.(absTOchar(il).ne."LOC"))then
  !    il = ib
  !    ir = ia
  ! endif
  ! if((absTOchar(ib).eq."SRC").and.(absTOchar(ia).eq."PWC"))then
  !    il = ib
  !    ir = ia
  ! endif
  !call Order_Indexes_1B(il,ir)
  
  
  iIrra = absTOirr(il)
  iIrrb = absTOirr(ir)
  aia   = absTOrel(il)
  aib   = absTOrel(ir)
  cha   = absTOchar(il)
  chb   = absTOchar(ir)
  !write(*,*) cha,chb,iIrra,iIrrb,ia,ib,aia,aib!,li,mi
  qres = 0.d0
  
  if(abs(mi)<=li)then
     ilm = li**2 + li + mi + 1
     if((cha.eq."LOC").and.(chb.eq."LOC"))then
        if(.not.GlobalIntegral%LocMomentsIsInit(iIrra,iIrrb,ilm))return
        Mat => GlobalIntegral%GetLocMoments(iIrra,iIrrb,ilm)
        
     endif
  else     

     if((cha.eq."LOC").and.(chb.eq."LOC"))then
        if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"MO_MO",iIrra,iIrrb))return
        Mat => GlobalIntegral%Get1B(OpLabel,"MO_MO",iIrra,iIrrb)
     endif
     if((cha.eq."LOC").and.(chb.eq."SRC"))then


        if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"MO_HY",iIrra,iIrrb))return
        Mat => GlobalIntegral%Get1B(OpLabel,"MO_HY",iIrra,iIrrb)
     endif
     if((chb.eq."LOC").and.(cha.eq."SRC"))then
        if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"MO_HY",iIrrb,iIrra))return
        Mat => GlobalIntegral%Get1B(OpLabel,"MO_HY",iIrrb,iIrra)
     endif
     if((cha.eq."SRC").and.(chb.eq."SRC"))then
        if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"HY_HY",iIrrb,iIrra))return
        Mat => GlobalIntegral%Get1B(OpLabel,"HY_HY",iIrrb,iIrra)     
     endif
     if((cha.eq."SRC").and.(chb.eq."PWC"))then
        l   = absTOl(ir)
        m   = absTOm(ir)
        nr  = absTOnr(ir)-NfirstBspl+1
        ilm = l**2 + l + m + 1
        if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"HY_BS",iIrra,ilm))return
        Mat => GlobalIntegral%Get1B(OpLabel,"HY_BS",iIrra,ilm)

        !***I DON TLIKE THIS: THE MATRIX HAS NOT THE CORRECT PHYSICAL DIMENSIONS
        if(aib.gt.Mat%NColumns())return
        qres = Mat%Element(aia,aib)
        return
     endif
     if((cha.eq."PWC").and.(chb.eq."PWC"))then
        l   = absTOl(ir)
        m   = absTOm(ir)
        nr  = absTOnr(ir)-NfirstBspl+1
        lp  = absTOl(il)
        mp  = absTOm(il)
        nrp = absTOnr(il)-NfirstBspl+1
        ilm  = l  **2 + l  + m  + 1
        ilmp = lp **2 + lp + mp + 1
        if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"BS_BS",ilmp,ilm))return
        Mat => GlobalIntegral%Get1B(OpLabel,"BS_BS",ilmp,ilm)
     endif
     if((cha.eq."PWC").and.(chb.eq."LOC"))return
     if((cha.eq."LOC").and.(chb.eq."PWC"))return

  endif
  if(associated(Mat)) qres = Mat%Element(aia,aib)


end function GET_1B_INTEGRAL_ASTRA

subroutine Order_Indexes_1B(i,j)
  use ModuleBasisJUAN
  implicit none
  integer, intent(inout) :: i , j
  integer                :: il, ir
  il = i
  ir = j
  if(absTOchar(il).eq.absTOchar(ir))then
     if(absTOirr(il).gt.absTOirr(ir))then
        il = j
        ir = i
     endif
  endif
  if((absTOchar(ir).eq."LOC").and.(absTOchar(il).ne."LOC"))then
     il = j
     ir = i
  endif
  if((absTOchar(j).eq."SRC").and.(absTOchar(i).eq."PWC"))then
     il = j
     ir = i
  endif
  i = il
  j = ir
end subroutine Order_Indexes_1B

real(kind(1d0)) function GET_2B_INTEGRAL_ASTRA(ia,ib,ic,id) result( qres )
  use ModuleUKRmolInterface
  use ukrmol_interface
  use ModuleBasisJUAN
  use ModuleIntegrals
  Use ModuleMatrix
  use special_functions, only: ipair
  implicit none
  integer, intent(in) :: ia,ib,ic,id
  integer :: ind(4),iind,i,lk,mk,l,m,lp,mp
  integer :: iIrr1, iIrr2, iIrr3, iIrr4
  integer :: ai1, ai2, ai3, ai4, nL, nS
  integer :: int_ID, ilm, ilmp, ilmk
  character*3 :: chi1, chi2, chi3, chi4
  real(kind(1d0)) :: mult, moments, ivalue
  character(len=:), allocatable :: OpLabel
  type(ClassMatrix) :: Mat
  
  nL=0
  nS=0
  ind(1) = ia
  ind(2) = ib
  ind(3) = ic
  ind(4) = id
  call IndexPermutation(ind)
  
  qres = 0.d0
  
  iIrr1 = absTOirr(ind(1))
  iIrr2 = absTOirr(ind(2))
  iIrr3 = absTOirr(ind(3))
  iIrr4 = absTOirr(ind(4))
  ai1   = absTOrel(ind(1))
  ai2   = absTOrel(ind(2))
  ai3   = absTOrel(ind(3))
  ai4   = absTOrel(ind(4))
  chi1  = absTOchar(ind(1))
  chi2  = absTOchar(ind(2))
  chi3  = absTOchar(ind(3))
  chi4  = absTOchar(ind(4))
!  write(*,*) iIrr1,iIrr2,iIrr3,iIrr4
!  write(*,*) chi1,chi2,chi3,chi4,ai1,ai2,ai3,ai4

  if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi3.eq."LOC").and.(chi4.eq."LOC"))then
     
     if(.not.GlobalIntegral%BielIsInitialized("LL_LL",iIrr1,iIrr2,iIrr3,iIrr4))then
        !write(*,*) GlobalIntegral%BielIsInitialized("LL_LL",iIrr1,iIrr2,iIrr3,iIrr4),iIrr1,iIrr2,iIrr3,iIrr4
        return
     endif
     qres = GlobalIntegral%GetBiel("LL_LL",iIrr1,iIrr2,iIrr3,iIrr4,ai1, ai2, ai3, ai4)
  endif
  if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi3.eq."LOC").and.(chi4.eq."SRC"))then
     if(.not.GlobalIntegral%BielIsInitialized("LL_LS",iIrr1,iIrr2,iIrr3,iIrr4))then
        !write(*,*) GlobalIntegral%BielIsInitialized("LL_LS",iIrr1,iIrr2,iIrr3,iIrr4),iIrr1,iIrr2,iIrr3,iIrr4
        return
     endif
     qres = GlobalIntegral%GetBiel("LL_LS",iIrr1,iIrr2,iIrr3,iIrr4,ai1, ai2, ai3, ai4)
  endif
  if((chi1.eq."LOC").and.(chi2.eq."SRC").and.(chi3.eq."LOC").and.(chi4.eq."SRC"))then
     if(.not.GlobalIntegral%BielIsInitialized("LS_LS",iIrr1,iIrr2,iIrr3,iIrr4))then
        !write(*,*) GlobalIntegral%BielIsInitialized("LS_LS",iIrr1,iIrr2,iIrr3,iIrr4),iIrr1,iIrr2,iIrr3,iIrr4
        return
     endif
     qres = GlobalIntegral%GetBiel("LS_LS",iIrr1,iIrr2,iIrr3,iIrr4,ai1, ai2, ai3, ai4)
  endif
  if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi3.eq."SRC").and.(chi4.eq."SRC"))then
     if(.not.GlobalIntegral%BielIsInitialized("LL_SS",iIrr1,iIrr2,iIrr3,iIrr4))then
        !write(*,*) GlobalIntegral%BielIsInitialized("LL_SS",iIrr1,iIrr2,iIrr3,iIrr4),iIrr1,iIrr2,iIrr3,iIrr4
        return
     endif
     qres = GlobalIntegral%GetBiel("LL_SS",iIrr1,iIrr2,iIrr3,iIrr4,ai1, ai2, ai3, ai4)
  endif




  
     !   if(.not.GlobalIntegral%LocMomentsIsInit(iIrra,iIrrb,ilm))return
     !    Mat = GlobalIntegral%GetLocMoments(iIrra,iIrrb,ilm)
         
     ! if((cha.eq."LOC").and.(chb.eq."LOC"))then
     !    if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"MO_MO",iIrra,iIrrb))return
     !    Mat = GlobalIntegral%Get1B(OpLabel,"MO_MO",iIrra,iIrrb)


  
  if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi4.eq."PWC"))then
     lp = absTOl(ind(4))
     mp = absTOm(ind(4))
     ilmp=lp**2 + lp + mp + 1
     ivalue = 0.d0
     OpLabel="MP_    "
     if(chi3.eq."SRC")then
        do lk = 0, 2*lmax
           do mk = -lk, lk
              int_ID = GlobalIntegral%Get_Multipole_ID( lk, mk )
              write(OpLabel(4:),'(I0.4)') int_ID - I1B_MULTIPO + 1
              ilmk = lk**2 + lk + mk + 1
              if(.not.GlobalIntegral%LocMomentsIsInit(iIrr1,iIrr2,ilmk))then
                ! write(*,*) "LocMoments NO",GlobalIntegral%LocMomentsIsInit(iIrr1,iIrr2,ilmk)
                 cycle
              else
                 !write(*,*) "LocMoments SI",GlobalIntegral%LocMomentsIsInit(iIrr1,iIrr2,ilmk)
              endif
              if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"HY_BS",iIrr3,ilmp,ilmk))then
                ! write(*,*) iIrr3,ilmp,lk,OpLabel
                ! write(*,*) "OneB NO",GlobalIntegral%OneBIsInitialized(OpLabel,"HY_BS",iIrr3,ilmp,ilmk)
                 cycle
              else
                ! write(*,*) "OneB SI",GlobalIntegral%OneBIsInitialized(OpLabel,"HY_BS",iIrr3,ilmp,ilmk)
              endif
              !write(*,*) "PASOOOO"
              Mat = GlobalIntegral%GetLocMoments(iIrr1,iIrr2,ilmk)
              
              moments = Mat%Element(ai1, ai2)
              !write(*,*) "moments",moments
              call Mat%Free()
              Mat = GlobalIntegral%Get1B(OpLabel,"HY_BS",iIrr3,ilmp,ilmk)
              !***I DON TLIKE THIS: THE MATRIX HAS NOT THE CORRECT PHYSICAL DIMENSIONS
              if(ai4.gt.Mat%NColumns())cycle
              mult    = Mat%Element(ai3, ai4)
             ! write(*,*) "mult",mult
              call Mat%Free()
              ivalue  = ivalue + moments * mult
           enddo
        enddo
        !write(*,*) "ACA1"
        !pause
     elseif(chi3.eq."PWC")then

        l = absTOl(ind(3))
        m = absTOm(ind(3))
        ilm =l**2 + l + m + 1
        do lk = 0, 2*lmax
           do mk = -lk, lk
              int_ID = GlobalIntegral%Get_Multipole_ID( lk, mk )
              write(OpLabel(4:),'(I0.4)') int_ID - I1B_MULTIPO + 1
              ilmk = lk**2 + lk + mk + 1
!              write(*,*) GlobalIntegral%LocMoments%LOC_LOC(iIrr1,iIrr2,ilmk)%ISINITIALIZED(),GlobalIntegral%Int_1B(int_ID)%PWC_PWC(ilm,ilmp)%ISINITIALIZED(),iIrr1,iIrr2
              if(.not.GlobalIntegral%LocMomentsIsInit(iIrr1,iIrr2,ilmk))cycle
              if(.not.GlobalIntegral%OneBIsInitialized(OpLabel,"BS_BS",ilm,ilmp,ilmk))cycle
              !write(*,*) "PASOOOO",iIrr1,iIrr2,ilmk
              Mat = GlobalIntegral%GetLocMoments(iIrr1,iIrr2,ilmk)
              moments = Mat%Element(ai1, ai2)
              !write(*,*) "moments",moments,OpLabel
              call Mat%Free()
              Mat = GlobalIntegral%Get1B(OpLabel,"BS_BS",ilm,ilmp,ilmk)
              mult    = Mat%Element(ai3, ai4)
              !write(*,*) "mult",mult,Mat%NRows(),Mat%NColumns(),ai3,ai4
              !write(*,*) Mat%IsDBanded(),Mat%IsFull()
              call Mat%Free()
              ivalue  = ivalue + moments * mult
             
           enddo
        enddo
     endif
     qres = ivalue
  endif

  
  ! if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi3.eq."SRC").and.(chi4.eq."PWC"))then
  ! endif
  ! if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi3.eq."PWC").and.(chi4.eq."PWC"))then
  ! endif
  if((chi1.eq."LOC").and.(chi2.eq."SRC").and.(chi3.eq."LOC").and.(chi4.eq."PWC"))then
     qres = 0.d0
  endif
  if((chi1.eq."LOC").and.(chi2.eq."PWC").and.(chi3.eq."LOC").and.(chi4.eq."PWC"))then
     qres = 0.d0
  endif


  
end function GET_2B_INTEGRAL_ASTRA

!*** THIS MUST NOT BEUSED AS IT IS BUILD WITH A DIFFERENT IRREP ORDER
real(kind(1d0)) function GET_1B_INTEGRAL_UKRMOL( ia, ib, Integral_ID, l, m ) result( qres )
  use ModuleUKRmolInterface
  use ukrmol_interface
  use ModuleIntegrals
  use ModuleBasisJUAN
  use ModuleConstants
  implicit none
  integer, intent(in) :: ia, ib, Integral_ID
  integer, intent(in) :: l, m

  integer*8       :: da8, db8, dflag, dl, dm
  integer         :: iIrra, iIrrb, aia, aib
  character*3     :: cha,chb
  real(kind(1d0)) :: factor
  dflag = int(Integral_ID,kind(dflag))

  aia = absTOukr( ia )
  aib = absTOukr( ib )
  if((aia.eq.-1).or.(aib.eq.-1))then
      write(*,*) "nbad index in GET_1B_INTEGRAL_EXTENDED, must not belong to the pwc",aia,aib
      stop
   endif
     
  da8    = UKRmol_orbitals%Get_IntIndex(aia)
  db8    = UKRmol_orbitals%Get_IntIndex(aib)
!  if(present(l).and.present(m))then
  if(abs(m)<=l)then
     dl=int(l,kind(dl))
     dm=int(m,kind(dm))
     ! write(*,*) "HOLA2",da8, db8, dflag , dl, dm,Integral_ID
     factor = sqrt(4.d0 * PI / dble(2 * dl + 1 ))
     qres = factor * GET_1B_INTEGRAL( da8, db8, dflag , dl, dm )
    ! write(*,*) "HOLA3",qres
  else
    ! write(*,*) "Astra2,",da8,db8
     qres = GET_1B_INTEGRAL( da8, db8, dflag )
  endif
end function GET_1B_INTEGRAL_UKRMOL



real(kind(1d0)) function GET_1B_INTEGRAL_EXTENDED( ia, ib, Integral_ID, li, mi, Bsplinex, iprim ) result( qres )
  use ModuleUKRmolInterface
  use ModuleSystemUtils
  use ukrmol_interface
  use ModuleIntegrals
  use ModuleBasisJUAN
  use ModuleBspline
  use ModuleParametersJUAN
  use ModuleConstants
  implicit none
  integer, intent(in) :: ia, ib, Integral_ID
  integer, intent(in) :: li, mi
  type(ClassBspline), intent(in)             :: Bsplinex
  type(BasisElementInfo), intent(in)             :: iprim(1:NprimUkrmol)
  integer*8   :: da8, db8, dflag, dl, dm
  integer     :: fa , fb, flag, il, ir
  integer     :: iIrra, iIrrb, aia, aib
  character*3 :: cha, chb
  integer     :: l, m, lp, mp, nr, nrp, lk, mk, k
  real(kind(1d0)) :: ivalue, ivalue2, suml, summ, sumk,nmom
  real(kind(1d0)) :: nor,factor
  real(kind(1d0)), external :: ThreeXlmIntegral, Xlm
  real(kind(1d0))                    :: parvec(2) = 0.d0
  procedure(D2DFun) , pointer    :: fptr1, fptr2

  ! NOTE THAT MOMENTS IS DIFFEREMT FROM MULTIPOLES. THE PROGRAM IS SET TO GET THE MOMENTS 4*PI/(2*L+1) r^(l) Xlm(\hat{r}) FOR THE
  ! LOC-LOC, LOC-SRC AND SRC-SRC ORBITALS AND THE MULTIPOLES r^(-l-1) Xlm(\hat{r}) FOR THE
  ! SRC-PWC AND PWC-PWC ONES.

  
  flag = Integral_ID
  il = ia
  ir = ib
  dflag = int(Integral_ID,kind(dflag))
  qres  = 0.d0
  aia   = absTOukr( il )
  aib   = absTOukr( ir )
  if((aia.ne.-1).and.(aib.ne.-1))then
     da8    = UKRmol_orbitals%Get_IntIndex(aia)
     db8    = UKRmol_orbitals%Get_IntIndex(aib)
     if(abs(mi)<=li)then  !MOMENTS 
        dl=int(li,kind(dl))
        dm=int(mi,kind(dm))
       ! write(*,*)da8, db8, dflag , dl, dm
        factor = sqrt(4.d0 * PI / dble(2 * dl + 1 ))
        qres = factor * GET_1B_INTEGRAL( da8, db8, dflag , dl, dm )
     else
        qres = GET_1B_INTEGRAL( da8, db8, dflag )
     endif
  endif


  if((aia.eq.-1).and.(aib.ne.-1))then
     il = ib
     ir = ia
  endif
  aia   = absTOukr( il )
  aib   = absTOukr( ir )
  cha   = absTOchar( il )
  chb   = absTOchar( ir )
  
  if((aia.ne.-1).and.(aib.eq.-1))then
     if(cha.ne."LOC")then
        l   = absTOl(ir)
        m   = absTOm(ir)
        nr  = absTOnr(ir)
        fa  = aia!UKRmol_orbitals%Get_IntIndex(aia)
        ivalue = 0.d0
        fptr1 => power_r
        fptr2 => partial_r
        !write(*,*) cf(:,fa)
        do k=sum(nMO)+1,NprimUkrmol
           if(abs(cf(k,fa)).gt.0.d0)then
              lp  = iprim(k)%l
              mp  = iprim(k)%m
              nrp = iprim(k)%cr
              if(abs(nrp-nr).ge.BsplineOrder)cycle
                 nor=NB(nrp)*NB(nr)
                  parvec(1)=0.d0
                 if((lp.eq.l).and.(mp.eq.m))then
                    if((flag.eq.4).or.(flag.eq.2))then
                       !write(*,*) ivalue2*nor
                       !pause
                       !RADIAL KINETIC ENERGY
                       parvec(1)=0.d0
                       ivalue2=-Bsplinex.Integral(fptr1,nrp,nr,0,2,parvec=parvec)                         
                       parvec(1)=-2.d0
                       !CENTRIFUGAL BARRIER
                       ivalue2=ivalue2+dble(l*(l+1))*Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)                      
                       ivalue=ivalue+cf(k,fa)*ivalue2*nor*0.5d0
                    end IF
                    If(flag.eq.1)then
                       !OVERLAP
                       parvec(1)=0.d0
                       ivalue2=Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)*nor
                       ivalue=ivalue+cf(k,fa)*ivalue2
                    end if
                 endIf
                 if((flag.eq.4).or.(flag.eq.3))then
                    !CENTRAL POTENTIAL
                    ivalue2=0.d0
                    do lk=abs(l-lp),l+lp
                       parvec(1)=-(lk+1)*1.d0
                       suml= Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)*nor
                       nmom = -Zn2*(Rn2**lk)*((4*PI)/dble(2*lk+1))
                       summ=0.d0
                       do mk=-lk,lk
                          summ=summ+ThreeXlmIntegral(lp,mp,lk,mk,l,m)*(Xlm(0.d0,0.d0,lk,mk)+Xlm(PI,0.d0,lk,mk))*nmom
                       end do
                       ivalue2=ivalue2+suml*summ
                    end do
                    ivalue=ivalue+cf(k,fa)*ivalue2
                 endif
                 !! MULTIPOLE INTEGRALS    !Need to be tested
                 if(flag.eq.5)then                                               
                    lk= li
                    mk= mi
                    parvec(1)=-1.d0*(lk+1)                     
                    ivalue2=Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)*nor
                   ! write(*,*) k,fa,lp,mp,l,m,nrp,nr

                    ivalue=ivalue+cf(k,fa)*ivalue2*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
                   ! write(*,*) cf(k,fa),ThreeXlmIntegral(lp,mp,lk,mk,l,m),ivalue2,ivalue
                 end if
           endif
        enddo
        qres = ivalue
     end if
  endif

  if((aia.eq.-1).and.(aib.eq.-1))then
     l   = absTOl(ib)
     m   = absTOm(ib)
     nr  = absTOnr(ib)
     lp   = absTOl(ia)
     mp   = absTOm(ia)
     nrp  = absTOnr(ia)
     fptr1 => power_r
     fptr2 => partial_r
     ivalue = 0.d0
     nor=NB(nrp)*NB(nr)
     if((lp.eq.l).and.(mp.eq.m))then
        if((flag.eq.4).or.(flag.eq.2))then
           !RADIAL KINETIC ENERGY
           parvec(1)=0.d0
           ivalue2=-Bsplinex.Integral(fptr1,nrp,nr,0,2,parvec=parvec)                         
           parvec(1)=-2.d0
           !CENTRIFUGAL BARRIER
           ivalue2=ivalue2+dble(l*(l+1))*Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)                      
           ivalue=ivalue+ivalue2*nor*0.5d0
        end IF
        If(flag.eq.1)then
           !OVERLAP
           parvec(1)=0.d0
           ivalue2=Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)*nor
           ivalue=ivalue+ivalue2
        end if
     endIf
     if((flag.eq.4).or.(flag.eq.3))then
        !CENTRAL POTENTIAL
        ivalue2=0.d0
        do lk=abs(l-lp),l+lp
           parvec(1)=-(lk+1)*1.d0
           suml= Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)*nor
           nmom = -Zn2*(Rn2**lk)*((4*PI)/dble(2*lk+1))
           summ=0.d0
           do mk=-lk,lk
              summ=summ+ThreeXlmIntegral(lp,mp,lk,mk,l,m)*(Xlm(0.d0,0.d0,lk,mk)+Xlm(PI,0.d0,lk,mk))*nmom
           end do
           ivalue2=ivalue2+suml*summ
        end do
        ivalue=ivalue+ivalue2
     endif
     !! MULTIPOLE INTEGRALS    !Need to be tested
     if(flag.eq.5)then                                               
        lk= li
        mk= mi
        parvec(1)=-1.d0*(lk+1)                       
        ivalue2=Bsplinex.Integral(fptr1,nrp,nr,0,0,parvec=parvec)*nor                
        ivalue=ivalue+ivalue2*ThreeXlmIntegral(lp,mp,lk,mk,l,m)
     end if
     qres = ivalue
  endif

contains

  Pure fUNCTION partial_r( r , parvec ) RESULT(y)
    implicit none
    REAL(kind(1d0)), INTENT(in) :: r
    DoublePrecision, optional, intent(in) :: parvec(*)
    REAL(kind(1d0)) :: y
    y = (min(r,parvec(1))**parvec(2))/(max(r,parvec(1))**(parvec(2)+1.d0))
  END FUNCTION partial_r




  Pure DoublePrecision function power_r(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = x**parvec(1)
  end function power_r

  
end function GET_1B_INTEGRAL_EXTENDED




real(kind(1d0)) function GET_2B_INTEGRAL_EXTENDED(ia,ib,ic,id, Bsplinex, iprim) result( qres )
  use ModuleUKRmolInterface
  use ukrmol_interface
  use ModuleBasisJUAN
  use ModuleBspline
  use ukrmol_interface
  implicit none
  integer               , intent(in) :: ia,ib,ic,id
  type(ClassBspline)    , intent(in) :: Bsplinex
  type(BasisElementInfo), intent(in) :: iprim(1:NprimUkrmol)
  integer*8   :: i1, i2, i3, i4
  integer     :: ind(4),nL,nP,i,lk,mk,int_ID
  integer     :: iIrr1, iIrr2, iIrr3, iIrr4 
  character*4 :: chi1, chi2,chi3,chi4
  real(kind(1d0)) :: moments, mult, ivalue
  real(kind(1d0)), external :: GET_1B_INTEGRAL_EXTENDED


  ind(1) = ia
  ind(2) = ib
  ind(3) = ic
  ind(4) = id
  
  nL = 0
  nP = 0

  qres = 0.d0
  do i =1, 4
     if(absTOchar(ind(i)).eq."LOC")then
        nL = nL + 1
     endif
     if(absTOchar(ind(i)).eq."PWC")then
        nP = nP + 1
     endif
  enddo
  if(nL.lt.2)then
     write(*,*) "GET_2B_INTEGRAL_EXTENDED: the minimum number of LOC orbitals must be 2... it is",nL
     do i =1, 4
        write(*,*) absTOirr(ind(i))
     enddo
     pause
     return
  endif

  if(nP.eq.0)then
     i1 = UKRmol_orbitals%Get_IntIndex( absTOukr(ind(1 )))
     i2 = UKRmol_orbitals%Get_IntIndex( absTOukr(ind(2 )))
     i3 = UKRmol_orbitals%Get_IntIndex( absTOukr(ind(3 )))
     i4 = UKRmol_orbitals%Get_IntIndex( absTOukr(ind(4 )))
     qres = GET_2B_INTEGRAL(i1, i2, i3, i4)
     return
  endif

  call IndexPermutation(ind)
  
  iIrr1 = absTOirr(ind(1))
  iIrr2 = absTOirr(ind(2))
  iIrr3 = absTOirr(ind(3))
  iIrr4 = absTOirr(ind(4))
  i1   = absTOrel(ind(1))
  i2   = absTOrel(ind(2))
  i3   = absTOrel(ind(3))
  i4   = absTOrel(ind(4))
  chi1  = absTOchar(ind(1))
  chi2  = absTOchar(ind(2))
  chi3  = absTOchar(ind(3))
  chi4  = absTOchar(ind(4))

!  write(*,*) chi1,chi2,chi3,chi4,i1,i2,i3,i4
  if((chi1.eq."LOC").and.(chi2.eq."SRC").and.(chi3.eq."LOC").and.(chi4.eq."PWC"))return
  if((chi1.eq."LOC").and.(chi2.eq."PWC").and.(chi3.eq."LOC").and.(chi4.eq."PWC"))return

  if((chi1.eq."LOC").and.(chi2.eq."LOC").and.(chi4.eq."PWC"))then
     ivalue = 0.d0
     int_ID = UKRMol_Get_Integral_Index('Property integrals')
     do lk = 0, 2*lmax
        do mk = -lk, lk
           moments = GET_1B_INTEGRAL_EXTENDED( ind(1), ind(2), int_ID, lk, mk, Bsplinex, iprim )
           mult    = GET_1B_INTEGRAL_EXTENDED( ind(3), ind(4), int_ID, lk, mk, Bsplinex, iprim )
           
           !write(*,*) moments,mult,lk,mk, i3,i4,ind(3), ind(4)
           ivalue  = ivalue + moments * mult
        enddo
     enddo
     qres = ivalue
  endif
  
  
end function GET_2B_INTEGRAL_EXTENDED



Subroutine IndexPermutation(ind)
!  use ModuleUKRmolInterface
!  use ukrmol_interface
  use ModuleBasisJUAN
!  use ModuleIntegrals
!  use special_functions, only: ipair
  implicit none
  integer, intent(inout) :: ind(4)
  integer :: i, nL, nS

  nL=0
  nS=0
  
  do i =1, 4
     if(absTOchar(ind(i)).eq."LOC")nL = nL + 1
  enddo
  if(nL.lt.2)then
     !return
     write(*,*) "IndexPermutation from GET_2B_INTEGRAL_ASTRA: the minimum number of LOC orbitals must be 2",nL
     do i =1, 4
        write(*,*) absTOirr(ind(i))
     enddo
  endif

  !.. Ordering irreps for the case LL_LL
  if(nL.eq.4)then
     if(max(absTOirr(ind(3)),absTOirr(ind(4))) > max(absTOirr(ind(1)),absTOirr(ind(2))))then
        call permutation4(ind)
     endif
     do i=1,3,2
        if(absTOirr(ind(i)) < absTOirr(ind(i+1)))then
           call permutation2( ind , i)
        endif
     enddo
  endif
  !.. This work for LL_LX and LX_LX  (X = S or P) to put L on the left.
  do i = 1, 3, 2
     if(((absTOchar(ind(i)).eq."SRC").or.(absTOchar(ind(i)).eq."PWC")).and.(absTOchar(ind(i+1)).eq."LOC"))then
        call permutation2( ind , i)
     endif
  enddo
  if(nL.eq.3)then
     if((absTOchar(ind(2)).eq."SRC").or.(absTOchar(ind(2)).eq."PWC"))then
        call permutation4(ind)         
     end if
     if(absTOirr(ind(1)) < absTOirr(ind(2)))then
        call permutation2( ind , 1)
     endif
  endif
  if(nL.eq.2)then
     !LS_LS  LP_LP
     if((absTOchar(ind(1)).eq."LOC").and.(absTOchar(ind(3)).eq."LOC"))then
        if(absTOchar(ind(2)).eq.absTOchar(ind(4)))then
           if((absTOirr(ind(1)).lt.absTOirr(ind(3))).or.((absTOirr(ind(1)).eq.absTOirr(ind(3))) &
                .and.(absTOirr(ind(4)).gt.absTOirr(ind(2)))))then
              call permutation4(ind)
           endif
        else
           ! LS_LP  LP_LS
           if(absTOchar(ind(2)).eq."PWC")then
              call permutation4(ind)
           endif
        endif
     endif
     !  SS_LL  SP_LL   PS_LL PP_LL ---> LL_XW
    ! write(*,*) absTOchar(ind(1)),absTOchar(ind(2)),absTOchar(ind(3)),absTOchar(ind(4))
    ! write(*,*) absTOirr(ind(1)),absTOirr(ind(2)),absTOirr(ind(3)),absTOirr(ind(4))
     !write(*,*) "aca 1"
     if((absTOchar(ind(3)).eq."LOC").and.(absTOchar(ind(4)).eq."LOC"))then
        
       ! write(*,*) "aca 2"
        
        call permutation4(ind)
        
        !write(*,*) "aca 3",absTOirr(ind(1)),absTOirr(ind(2)),absTOirr(ind(3)),absTOirr(ind(4))
        
        if(absTOirr(ind(1)) < absTOirr(ind(2)))then
                   call permutation2( ind , 1)
        endif
     endif
     if((absTOchar(ind(1)).eq."LOC").and.(absTOchar(ind(2)).eq."LOC"))then
        if(absTOirr(ind(1)) < absTOirr(ind(2)))then
                   call permutation2( ind , 1)
        endif
     endif




     
     if(absTOchar(ind(3)).eq.absTOchar(ind(4)))then
        if(absTOirr(ind(3)) < absTOirr(ind(4)))then
           call permutation2( ind , 3)
        endif
     else
        if(absTOchar(ind(3)).eq."PWC")then
           call permutation2( ind , 3)
        endif
     endif
  endif
  
contains

  subroutine permutation4( vec )
    implicit none
    integer, INTENT(inout) :: vec(4)
    integer :: ij
    ij = vec(1)
    vec(1) = vec(3)
    vec(3) = ij
    ij = vec(2)
    vec(2) = vec(4)
    vec(4) = ij
  END subroutine permutation4

  subroutine permutation2( vec , i)
    implicit none
    integer, INTENT(inout) :: vec(4)
    integer, INTENT(in) :: i
    integer :: ij
    ij = vec(i)
    vec(i) = vec(i+1)
    vec(i+1) = ij
  END subroutine permutation2
  
end subroutine IndexPermutation



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









  
!position_flag=4     ===>  electron Hamiltonian
!position_flag=2     ===>  Kinetic Energy
!position_flag=3     ===>  Nuclear attraction
!positron_flag=11    ===>  property integrals

!***NOTWORKING
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
  use ModuleIntegrals
  use ModuleUKRmolInterface
  implicit none

  !    REAL(kind(1d0)), INTENT(in) :: r
  integer, intent(in) :: ia,ib,ic,id
  integer, intent(in) :: flag
  procedure(D2DFun) , pointer    :: fptr1,fptr2
  type(ClassBspline)    , intent(in)             :: Bsplinex
  !integer  :: BSplineGetOrder
  type(BasisElementInfo), intent(in)             :: iprim(1:NprimUkrmol),ielement(1:Nbasis)   !
  integer :: ja,jb,jc,jd,i,j,k,l,m,lp,mp,lk,mk,Bra,Ket
  integer :: iflag_mp
  integer*8 :: i8,j8,k8,l8,flagd
  real(kind(1d0)) :: ivalue,ivalue2,suml,summ,sumk
  real(kind(1d0)) :: ME,nor
  real(kind(1d0)), external :: ThreeXlmIntegral,Xlm
  real(kind(1d0))                    :: parvec(2) = 0.d0
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_ASTRA
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_UKRMOL
  real(kind(1d0)), external  :: GET_1B_INTEGRAL_EXTENDED
  real(kind(1d0)), external  :: GET_2B_INTEGRAL_ASTRA
  logical         :: twoBcase

  iflag_mp=UKRMol_Get_Integral_Index('Property integrals')

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
    ! write(*,*) ielement(j)%spin,ielement(i)%spin,.not.ielement(j)%EBtype,.not.ielement(i)%EBtype
     ! do m=1,Nbasis
     !    write(*,*) m,ielement(m)%EBtype
     !    end do
     ! pause
     !    write(*,*) i,j,ielement(j)%EBtype,ielement(i)%EBtype,flag
     If(ielement(j)%spin.eq.ielement(i)%spin)Then
        If((.not.ielement(j)%EBtype).and.(.not.ielement(i)%EBtype))Then
          ! write(*,*)"AAA"
           i8=ielement(j)%nr
           j8=ielement(i)%nr
           
          ! write(*,*) ielement(j)%nr,ielement(i)%nr,flagd,GET_INTEGRALS(i8,j8,0,0,flagd)
           ME=GET_INTEGRALS(i8,j8,0,0,flagd)
           !ME=GET_1B_INTEGRAL_EXTENDED(ielement(j)%nr,ielement(i)%nr,flag,0,1)
           ! write(*,*)  ME,i8,j8,0,0,flagd
           ! ME=GET_1B_INTEGRAL(i8,j8,flagd)
           ! write(*,*) "aca2",ME
        End IF
        If(((.not.ielement(j)%EBtype).and.ielement(i)%EBtype).or.(ielement(j)%EBtype.and.(.not.ielement(i)%EBtype)))Then
           ivalue=0.d0
           !write(*,*) "aca2",i,j,ielement(j)%EBtype,ielement(i)%EBtype
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
              !write(*,*) "l,lp",l,lp,ielement(jd)%nr,jd,NprimUkrmol
              If(abs(cf(k,ielement(jd)%nr)).gt.0.d0)then
                 If(iprim(k)%cr.ge.BsplineOrder)then
                    nor=NB(Bra)*NB(Ket)
                    If((lp.eq.l).and.(mp.eq.m))Then
                       If((flag.eq.4).or.(flag.eq.2))then
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

                    If((flag.eq.4).or.(flag.eq.3))then
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
                    If(flag.eq.11)then                                               
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
              If((flag.eq.4).or.(flag.eq.2))then

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

           If((flag.eq.4).or.(flag.eq.3))then
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
           If(flag.eq.11)then                                               
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
              !ME=GET_2B_INTEGRAL_ASTRA(ielement(ja)%nr,ielement(jb)%nr,ielement(jc)%nr,ielement(jd)%nr)
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
                    !summ=summ+GET_1B_INTEGRAL_ASTRA(ielement(jc)%nr,ielement(jd)%nr, iflag_mp, lk, mk )

                    
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
                    !summ=summ+GET_1B_INTEGRAL_ASTRA(ielement(ja)%nr,ielement(jb)%nr, iflag_mp, lk, mk )
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
                                !summ=summ+GET_1B_INTEGRAL_ASTRA(ielement(jc)%nr,ielement(jd)%nr, iflag_mp, lk, mk )
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
                                !summ=summ+GET_1B_INTEGRAL_ASTRA(ielement(ja)%nr,ielement(jb)%nr, iflag_mp, lk, mk )
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



