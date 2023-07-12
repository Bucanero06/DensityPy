subroutine InitFromFiles
  use ModuleErrorHandling
  use ModuleBasisJUAN
  use basis_data_generic_mod
  implicit none

  integer             :: iostat, i,j, iIrrep ,l,m,iexerror
  character(len=500)  :: iomsg 
  character(len=1000) :: line 
  character(len=*),  parameter :: UKRMOL_LOG_FILE = "log_file.1"
  character(len=*),  parameter :: UKRMOL_INP_FILE = "inp"
  integer*8  ::in8,jn8,number_of_bto_shells,number_of_cgto_shells
  type(BTO_shell_data_obj), allocatable :: dummy_bto(:)
  type(CGTO_shell_data_obj), allocatable :: dummy_cgto(:)
  logical*8, allocatable :: list(:)
  character(len=132)  :: names
  integer*8 :: no_sym, nob(8), nlmq
  ! MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT MODULO INPUT
  call execute_command_line("cp ./moints ./fort.1",exitstat=iexerror)
  if(iexerror/= 0)  stop 'error occured trying to cp moints file to fort.1!'
  call mpi_mod_start(master_writes_to_stdout)!, allow_shared_memory)
  scalem=  1.0_cfp
  print*, 'loaded ukrmol_interface successfully'
  open(nfti,form='unformatted',access='stream')
  read(nfti) ukrmolp_header
  close(nfti)
  print*, ukrmolp_header, data_file_obj_id
  IF (ukrmolp_header .eq. data_file_obj_id) THEN 
     print*, 'UKRMol+ format found on the input integrals file'!,irr_series
       
     call READ_UKRMOLP_BASIS(NFTI)
     nIrreps=molecular_orbital_basis%no_irr
     allocate(irr_series(nIrreps),irr_counter(nIrreps))
     nsym1=molecular_orbital_basis%no_irr
     Do in8=1,nIrreps
        i=in8
        irr_series(i)=molecular_orbital_basis%get_number_of_orbitals(in8)
     End Do
     call READ_UKRMOLP_INTS(nfte,nfti,lembf,nint1e,nint2e,&
          nocsf,nfta,isymtp,nsym1,irr_series,iposit,scalem,name,&
          nalm,qmoln)
     !  RETURN
  ELSE
     WRITE(*,'(a,";",a)')ukrmolp_header,data_file_obj_id
     STOP "An unknown header on the integrals file."
  ENDIF

  !  call  READ_UKRMOLP_BASIS(NFTINT)
  write(*,*) "read",NFTI
  allocate(cf(1:NprimUkrmol,1:NprimUkrmol))
  ! allocate(cf(1:Nukrmol,1:Nukrmol))
  cf=0.d0
  !cfi=0.d0
  call molecular_orbital_basis%get_orbital_coefficient_matrix(cf)
  call READ_UKRMOLP_PROPERTY_INTS(NFTI,3)

  SymmetryIdentifier=molecular_orbital_basis%pg
  allocate(nTOT(nIrreps))
  nTOT = irr_series
  allocate(sTOT(1:nIrreps))
  sTOT=0
  do iIrrep = 1,nIrreps
     sTOT(iIrrep) = sum(nTOT(1:iIrrep))
  enddo
  ! Do in8=1,sum(nTOT)
  !    write(*,*)in8,molecular_orbital_basis%get_orbital_symmetry(in8)
  ! End Do
  allocate(nLOC(nIrreps))
  allocate(nSRC(nIrreps))
  nLOC=0
  nSRC=0
  Do in8=1,nIrreps
     write(*,*) in8,molecular_orbital_basis%get_number_of_orbitals(in8)
     allocate(list(1:molecular_orbital_basis%get_number_of_orbitals(in8)))
     call molecular_orbital_basis%get_continuum_flags(in8,list)
     i=in8
     Do jn8=1,molecular_orbital_basis%get_number_of_orbitals(in8)
        If(list(jn8))Then
           nSRC(i)=nSRC(i)+1
        Else
           nLOC(i)=nLOC(i)+1
        End IF
     End Do
     deallocate(list)
  End Do
  allocate(nMO(1:nIrreps))
  nMO=nLOC
  write(*,*) SymmetryIdentifier
  write(*,*) nIrreps
  write(*,*) "-------"
  write(*,*) nLOC
  write(*,*) "-------"
  write(*,*) nSRC
  write(*,*) "-------"
  write(*,*) nTOT
  write(*,*) "-----------"
!pause
  write(*,*) "Number of functions",molecular_orbital_basis%ao_basis%number_of_functions
  write(*,*) "Number of continuum functions",molecular_orbital_basis%ao_basis%n_cont_fns
  write(*,*) "Number of molecular functions",molecular_orbital_basis%ao_basis%n_cgto_functions
  write(*,*) "Numberof shells",molecular_orbital_basis%ao_basis%number_of_shells


  call GET_NAME_SYM(names,no_sym,nob,nlmq)
  Lmax_mp=int(sqrt(dble(nlmq)))-1
  Nmultipoles=(Lmax_mp+1)**2     !dimension of the above vectors
  allocate(l_mp(1:Nmultipoles),m_mp(1:Nmultipoles))  !values of the l and m as a function of the flag value in ukrmol
  allocate(flag_ml(0:Lmax_mp,-Lmax_mp:Lmax_mp))      !the flag which is associated to a given pair {l,m}
  j=0
  Do l=0,Lmax_mp
     Do m=-l,l
        j=j+1
        l_mp(j)=l
        m_mp(j)=m
        flag_ml(l,m)=j
     End Do
  End Do

  !number_of_bto_shells=molecular_orbital_basis%ao_basis%number_of_shells   !this gives the total
  call molecular_orbital_basis%ao_basis%get_all_CGTO_shells(dummy_cgto,number_of_cgto_shells)
  call molecular_orbital_basis%ao_basis%get_all_BTO_shells(dummy_bto,number_of_bto_shells)
  !!!For our code it would be :
  
  Lmax=molecular_orbital_basis%ao_basis%shell_descriptor(6,molecular_orbital_basis%ao_basis%number_of_shells)
  min_bspline_l = 0
  max_bspline_l = Lmax
  bspline_indices=0
  bspline_indices(:,0:Lmax)=dummy_bto(number_of_bto_shells)%bspline_index  !!!! <<<===  THIS IS THE NUMBER
  Do in8=number_of_cgto_shells+1,number_of_cgto_shells+number_of_bto_shells
     l=molecular_orbital_basis%ao_basis%shell_descriptor(6,in8)
     jn8=molecular_orbital_basis%ao_basis%shell_descriptor(2,in8)
     bspline_indices(1,l)=min(bspline_indices(1,l),dummy_bto(jn8)%bspline_index)
  End Do

! write(*,*) bspline_indices(1,0:)
! write(*,*) bspline_indices(2,0:)
 ! pause
   bspline_grid_start=dummy_bto(1)%bspline_grid%A
   bspline_order=dummy_bto(1)%bspline_grid%order
   no_bsplines=dummy_bto(1)%bspline_grid%n
! !  write(*,*) bspline_grid_start
! !  write(*,*) bspline_order
! !  write(*,*) no_bsplines
   a=dummy_bto(1)%bspline_grid%B
   no_nodes=no_bsplines-bspline_order+2
   da=(a-bspline_grid_start)/(no_nodes-1)
!   ! write(*,*) a,da
!   ! write(*,*) bspline_grid_start
!   ! write(*,*) bspline_order
!   ! write(*,*) no_bsplines




  


  ! write(*,*) "------------",size(molecular_orbital_basis%ao_basis%indices_to_shells(1,:)) 
  ! write(*,*) molecular_orbital_basis%ao_basis%indices_to_shells(2,:)
  ! ! De aqui sale iprim, solo debemos definir el orden para el numero angular m
  ! Do in8=1,molecular_orbital_basis%ao_basis%number_of_shells
  !    write(*,'(I3.3,1x,I3.3,1x,I3.3,1x,I3.3,1x,I3.3,1x,I3.3,1x,I3.3)') in8, &
  !         molecular_orbital_basis%ao_basis%shell_descriptor(1,in8), &
  !         molecular_orbital_basis%ao_basis%shell_descriptor(2,in8), &
  !         molecular_orbital_basis%ao_basis%shell_descriptor(3,in8), &
  !         molecular_orbital_basis%ao_basis%shell_descriptor(4,in8), &
  !         molecular_orbital_basis%ao_basis%shell_descriptor(5,in8), &
  !         molecular_orbital_basis%ao_basis%shell_descriptor(6,in8)
  ! End Do

  ! number_of_bto_shells=molecular_orbital_basis%ao_basis%number_of_shells
  ! call molecular_orbital_basis%ao_basis%get_all_BTO_shells(dummy_bto,number_of_bto_shells)
  ! write(*,*) size(dummy_bto)
  ! write(*,*) dummy_bto(:)%bspline_index
  ! write(*,*) dummy_bto(1)%bspline_grid%A
  ! write(*,*) dummy_bto(1)%bspline_grid%order
  ! write(*,*) dummy_bto(1)%bspline_grid%n
  ! write(*,*) size(molecular_orbital_basis%so2mo_range(:,1)),size(molecular_orbital_basis%so2mo_range(1,:))
  ! !write(*,*) molecular_orbital_basis%so2mo_range(1,:)
  ! !write(*,*) "-----"
  ! write(*,*) molecular_orbital_basis%so2mo_range(2,:)
  ! pause

  ! !> The range of the radial basis of which this function is a member.
  !    real(kind=cfp) :: A, B = -1
  !    !> Number of break-points in the basis.
  !    integer :: no_bps = -1
  !    !> Order of the radial B-spline basis.
  !    integer :: order = -1
  !    !> Number of B-splines.
  !    integer :: n
  !    !> Index of the first radial B-spline whose first derivative at point r=A is 0.
  !    integer :: ind_0_der = 0
  !    !> Array of coefficients in the B-spline basis that can be used to evaluate B-spline basis functions.
  !    real(kind=cfp), allocatable :: bcoef(:)
  !    !> Number of knots.
  !    integer :: no_knots = 0
  !    !> Array of knots.
  !    real(kind=cfp), allocatable :: knots(:)
  !    !> Work array used for evaluation of this spline.
  !    real(kind=cfp), allocatable :: work(:)
  !    !> Helper variable used for evaluation of B-splines.
  !    integer :: inbv = 0
  !    !> Relative precision (tolerance) for numerical calculation of integrals involving this B-spline grid.
  !    real(kind=cfp) :: tol = 0.0_cfp





  ! !.. Here we read log_file.1 and search for the basis size, number of irr, states, property integrals, etc.
  ! open(newunit = uid               , & 
  !      file    = UKRMOL_LOG_FILE   , & 
  !      status  ='old'              , & 
  !      form    ='formatted'        , & 
  !      action  ='read'             , & 
  !      iostat  = iostat            , & 
  !      iomsg   = iomsg             ) 
  ! if(iostat/=0)then 
  !    call ErrorMessage("File "//UKRMOL_LOG_FILE//" not found. Stop forced.") 
  !    stop 
  ! endif

  ! do 
  !    read(uid,"(a)",iostat=iostat) line 
  !    if(iostat/=0)then 
  !       call ErrorMessage("orbital and sym parameters not found in "//UKRMOL_LOG_FILE) 
  !       stop 
  !    endif
  !    line=adjustl(line) 
  !    i=index(line,"--------->couplings_type: precalculating the Gaunt coefficients for L_max = ")
  !    if(i>0)exit 
  ! enddo
  ! i=len_trim("--------->couplings_type: precalculating the Gaunt coefficients for L_max =")
  ! read(line(i+1:),*) Lmax_mp
  !This is set in the ukrmol code where it is indicated the numbers of solid harmonics included as one body operators
  ! write(*,*) "Lmax_mp",Lmax_mp !l=20 is set at line 265 of scatci_integrals.f90 (the original value was l=2) 
  ! Nmultipoles=(Lmax_mp+1)**2     !dimension of the above vectors
  ! allocate(l_mp(1:Nmultipoles),m_mp(1:Nmultipoles))  !values of the l and m as a function of the flag value in ukrmol
  ! allocate(flag_ml(0:Lmax_mp,-Lmax_mp:Lmax_mp))      !the flag which is associated to a given pair {l,m}
  ! j=0
  ! Do l=0,20
  !    Do m=-l,l
  !       j=j+1
  !       l_mp(j)=l
  !       m_mp(j)=m
  !       flag_ml(l,m)=j
  !    End Do
  ! End Do
  ! do 
  !    read(uid,"(a)",iostat=iostat) line 
  !    if(iostat/=0)then 
  !       call ErrorMessage("orbital and sym parameters not found in "//UKRMOL_LOG_FILE) 
  !       stop 
  !    endif
  !    line=adjustl(line) 
  !    i=index(line,"Point-group symmetry identifier: ")
  !    if(i>0)exit 
  ! enddo
  ! ! i=len_trim("Point-group symmetry identifier: ") 
  ! ! read(line(i+1:),*) SymmetryIdentifier 
  ! read(uid,"(a)") line 
  ! i=len_trim("Number of irreducible representations: ")
  ! read(line(i+1:),*) nIrreps 
  ! read(uid,*) 
  ! read(uid,*) 
  ! read(uid,*) 
  ! allocate(nLOC(nIrreps)) 
  ! read(uid,*) (nLOC(iIrrep),iIrrep=1,nIrreps) 
  ! read(uid,*) 
  ! allocate(nSRC(nIrreps)) 
  ! read(uid,*) (nSRC(iIrrep),iIrrep=1,nIrreps) 
  ! close(uid) 
  ! ! allocate(nTOT(nIrreps)) 
  ! ! nTOT = nLOC + nSRC
  ! allocate(sTOT(1:nIrreps))
  ! sTOT=0
  ! do iIrrep = 1,nIrreps
  !    sTOT(iIrrep) = sum(nTOT(1:iIrrep))
  ! enddo
  ! !.. Molecular orbitals at each irr comming from Dalton. This is passsed to scatci_integrals through molden and inp file
  ! allocate(nMO(1:nIrreps))
  ! nMO=nLOC
  ! write(*,*) SymmetryIdentifier
  ! write(*,*) nIrreps
  ! write(*,*) "-------"
  ! write(*,*) nLOC
  ! write(*,*) "-------"
  ! write(*,*) nSRC
  ! write(*,*) "-------"
  ! write(*,*) nTOT
  ! write(*,*) "-----------"
  ! close(uid)

!   !.. Here we read ukrmol input to get Bsplines grid parameters
!   open(newunit = uid               , & 
!        file    = UKRMOL_INP_FILE   , & 
!        status  ='old'              , & 
!        form    ='formatted'        , & 
!        action  ='read'             , & 
!        iostat  = iostat            , & 
!        iomsg   = iomsg             ) 
!   if(iostat/=0)then 
!      call ErrorMessage("File "//UKRMOL_INP_FILE//" not found. Stop forced.") 
!      stop 
!   endif
!   do 
!      read(uid,"(a)",iostat=iostat) line 
!      if(iostat/=0)then 
!         call ErrorMessage("orbital and sym parameters not found in "//UKRMOL_INP_FILE) 
!         stop 
!      endif
!      line=adjustl(line) 
!      i=index(line,"&target_data") 
!      if(i>0)exit 
!   enddo
!   read(uid,"(a)") line
!   i=index(line,"=")
!   j=index(line,",")
! !  read(line(i+1:j-1),*) a
!   do 
!      read(uid,"(a)",iostat=iostat) line 
!      if(iostat/=0)then 
!         call ErrorMessage("orbital and sym parameters not found in "//UKRMOL_INP_FILE) 
!         stop 
!      endif
!      line=adjustl(line) 
!      i=index(line,"bspline_grid_start") 
!      if(i>0)exit 
!   enddo
!   i=index(line,"=")
! !  read(line(i+1:),*) bspline_grid_start
!   read(uid,"(a)") line
!   i=index(line,"=")
! !  read(line(i+1:),*) bspline_order
!   read(uid,"(a)") line
!   i=index(line,"=")
! !  read(line(i+1:),*) no_bsplines
!   read(uid,*)
!   read(uid,"(a)") line
!   i=index(line,"=")
!   j=index(line,",")
! !  read(line(i+1:j-1),*) min_bspline_l
!   line(i:i)=" "
!   line(j:j)=" "
!   i=index(line,"=")
!   j=index(line,",")
! !  read(line(i+1:j-1),*) max_bspline_l
!   write(*,*) bspline_grid_start
!   write(*,*) bspline_order
!   write(*,*) no_bsplines
!   write(*,*) min_bspline_l,max_bspline_l
! !  no_nodes=no_bsplines-bspline_order+2
! !  da=(a-bspline_grid_start)/(no_nodes-1)
!   !bspline_indices=0
!   do l=min_bspline_l,max_bspline_l
!      read(uid,"(a)") line
!      i=index(line,"=")
!      !write(*,*) l,i,bspline_indices(1,l)
!    !  read(line(i+1:),*) bspline_indices(1,l)
!      !write(*,*) l,i,bspline_indices(1,l)
!      read(uid,"(a)") line
!      i=index(line,"=")
!   !   read(line(i+1:),*) bspline_indices(2,l)
!      read(uid,*)
!     ! write(*,*) "l=",l,bspline_indices(1,l),bspline_indices(2,l)
!   End Do
!   close(uid)

!  ! write(*,*) bspline_indices(1,0:)
!  ! write(*,*) bspline_indices(2,0:)
!  ! pause

end subroutine InitFromFiles



subroutine InitBsplines(Bspline,Bsplinex)
  use ModuleErrorHandling
  use ModuleBasisJUAN
  use ModuleBspline
  use ModuleSystemUtils

  implicit none

  type(ClassBspline), intent(inout)  :: Bspline,Bsplinex
  integer                            :: iostat, i, iIrrep ,l
  character(len=500)                 :: iomsg 
  character(len=1000)                :: line
  character(len=12)                  :: GridFile
  real(kind(1d0))                    :: Integral,r,dr,ivalue
  procedure(D2DFun), pointer         :: fptr
  real(kind(1d0))                    :: parvec(2) = 0.d0
  !write(*,*) "InitBsplines"
  
  !..  THIS IS DEFINED IN THE INPUT OF SCATCI_INTEGRALS. DEFINITIONS HERE MUST MATCH WITH THAT INPUT  vvvvvvvvvvvvv 
  da=(a-bspline_grid_start)/(no_nodes-1)  !this is the spacement of the radial grid in  scatci_integrals
  no_nodes=no_bsplines-bspline_order+2
  allocate(bsgrid(1:no_nodes))
  Do i=1,no_nodes
     bsgrid(i)=bspline_grid_start+(i-1)*da
  End do

  !.. Initialize Bspline: these are the ukrmol Bsplines
  call Bspline.Init(no_nodes,bspline_order,bsgrid)
  fptr => Power
  parvec(1)=0.d0
  ! GridFile="Bukrm000.dat"
  ! Do i=1,no_bsplines
  !    write(GridFile(6:8),'(I3.3)') i
  !    open(newunit = uid     , &
  !         file    = GridFile)
  !    dr=0.1d0
  !    Do r=0.d0,a,dr
  !       ivalue=Bspline%EvalBS(r,i)
  !       write(uid,*) r,ivalue
  !    End do
  !    close(uid)
  !    !  Integral = Bspline.Integral(fptr,i,Bs2,dBra,dKet)
  !    !  write(*,*) "Integral",i,Integral
  ! End do

  
  deallocate(bsgrid)

  !..  Normalization of internal Bsplines
  allocate(NBin(1:no_bsplines))
  NBin=0.d0
  fptr => Power
  parvec(1)=0.d0
  Do i=1,no_bsplines
     Integral = Bspline.Integral(fptr,i,i,0,0,parvec=parvec)
     NBin(i)=1.d0/sqrt(Integral)
  End Do

  !.. HERE IT IS THE DEFINITION OF OUR EXTERNAL BSPLINES GRID
  !NOW: the last "free" bsplines (the last one that do not see the border "a") extend over bspline_order+1
  !knot points, that is a spatial region equal to
  !da*bspline_order. It overlap to the left with other bspline_order-1 bsplines, the first of them starting at the point a-(2*bspline_order-1)
  !which is the minimum radius we have to set in the new grid with a spacing da up to a desired Rmax >a.

  Rmin=a-(2*bspline_order-2)*da
!  Rmax=a  !our new redius. If we do not want the external bsplines we set Rmax=a, ELSE we define a larger radius:
!  Rmax=30.1333333333333d0 
  dr=da  !we use same step than in ukrmol, except for the last node where we match Rmax
  !.. we check if (Rmax-Rmin)/dr corresponds to an integer. This is to correctly define the number of nodes of the external bsplines
  ivalue=(Rmax-Rmin)/dr
  If((ivalue-1.d0*int(ivalue)).gt.0.d0)then
     NNodes=Int(ivalue)+2
     allocate(bsgrid(1:NNodes))
     i=1
     Do r=Rmin,Rmax,dr
        bsgrid(i)=r
        i=i+1
     End do
     bsgrid(i)=Rmax
  Else
     NNodes=Int(ivalue)+1
     allocate(bsgrid(1:NNodes))
     i=1
     Do r=Rmin,Rmax,dr
        bsgrid(i)=r
        i=i+1
     End do
  End If
  BsplineOrder=bspline_order
  NBspl=NNodes+BsplineOrder-2
  Lmax=max_bspline_l
  !.. the size of the external basis is the total number of external bsplines, minus the first 2*(BsplineOrder-1) and minus the last one
  NlastBspl=NBspl-(BsplineOrder-1) !this definition is for comparison with ukrmol
  NfirstBspl=(2*BsplineOrder-1)
  NExtBspl=NlastBspl-NfirstBspl+1
  write(*,*) "NNodes      ",NNodes
  write(*,*) "BsplineOrder",BsplineOrder
  write(*,*) "NBspl       ",NBspl
  write(*,*) "Lmax        ",Lmax
  write(*,*) "NlastBspl   ",NlastBspl
  write(*,*) "NfirstBspl  ",NfirstBspl
  write(*,*) "NExtBspl    ",NExtBspl

  !.. Initialize external Bsplines
  call Bsplinex.Init(NNodes,BsplineOrder,bsgrid)
  n_bs=0   !the exponent of r in function powers_r
  fptr => Power
  parvec(1)=0.d0
  ! GridFile="BExte000.dat"
  ! Do i=1,NBspl!no_bsplines
  !    write(GridFile(6:8),'(I3.3)') i
  !    open(newunit = uid     , &
  !         file    = GridFile)
  !    dr=0.1d0
  !    Do r=dr,Rmax,dr
  !       ivalue=Bsplinex%EvalBS(r,i)
  !       write(uid,*) r,ivalue*1.0
  !    End do
  !    close(uid)
  ! End do

  !..  Normalization of external Bsplines
  allocate(NB(1:NBspl))
  NB=0.d0
  Do i=1,NBspl
     Integral = Bsplinex.Integral(fptr,i,i,0,0,parvec=parvec)
     NB(i)=1.d0/sqrt(Integral)
  End Do

!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
contains

  Pure DoublePrecision function Unity(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = 1.d0
  end function Unity

  Pure DoublePrecision function Power(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = x**parvec(1)
  end function Power

end subroutine InitBsplines


subroutine InitBasisParameters
  use ModuleErrorHandling
  use ModuleBasisJUAN
  implicit none

  integer           :: l,m,i,j,k,iIrreps
  integer, external :: irr_lm
  !.. This is valid for D2h. We have to implement the general case.
  !.. Here we clasify the angular basis between the different irr. This will allow us to locate the bsplibe basis in the ukrmol matrix
  allocate(nirr_lm(1:nIrreps))  !the number of different angular pairs belonging to each irr
  nirr_lm=0
  Do l=0,Lmax
     If(mod(l,2).eq.0)then
        nirr_lm(1)=nirr_lm(1)+(l/2+1)   
        nirr_lm(4)=nirr_lm(4)+(l/2)
        nirr_lm(6)=nirr_lm(6)+(l/2)
        nirr_lm(7)=nirr_lm(7)+(l/2)
     Else
        nirr_lm(8)=nirr_lm(8)+((l-1)/2)
        nirr_lm(5)=nirr_lm(5)+((l+1)/2)
        nirr_lm(3)=nirr_lm(3)+((l+1)/2)
        nirr_lm(2)=nirr_lm(2)+((l+1)/2)
     End If
  End do

  Do iIrreps=1,nIrreps
     write(*,*) iIrreps,nirr_lm(iIrreps)
  End Do
  write(*,*) "^^^^  number of partial waves per irr for the selected lmax"


  !.. Counter for l,m and number of radial elements nr at each irr. The lenght of those vectors is equal or lower than maxval(nirr_lm)
  allocate(irr_l(1:nIrreps,1:maxval(nirr_lm)),irr_m(1:nIrreps,1:maxval(nirr_lm)),irr_nr(1:nIrreps,1:maxval(nirr_lm)))
  irr_l=0
  irr_m=0
  irr_nr=0
  nirr_lm=0
  Do l=0,Lmax
     Do m=-l,l
        iIrreps=irr_lm(l,m)                !this is a function which gives irr as a function of (l,m)
        nirr_lm(iIrreps)=nirr_lm(iIrreps)+1      !we define this again and use it also as a counter
        irr_l(iIrreps,nirr_lm(iIrreps))=l
        irr_m(iIrreps,nirr_lm(iIrreps))=m
        irr_nr(iIrreps,nirr_lm(iIrreps))=bspline_indices(2,l)-bspline_indices(1,l)+1
       
     End Do
  End Do

  ! !New dimensions of the basis, divided by irr
   allocate(irr_total(1:nIrreps))
  irr_total=nTOT
  Do i=1,nIrreps
     irr_total(i)=irr_total(i)+NExtBspl*nirr_lm(i)
  End Do
  Do i=1,nIrreps
     write(*,*)i, irr_total(i)
  End Do
  write(*,*) "^^^^  extended basis",NExtBspl

  !.. Size of the primitive ukrmol basis (without removing states after orthonormalization)
  NprimUkrmol=0
  Do i=1,nIrreps
     l=0
      write(*,*) "a",nirr_lm(i)
      Do j=1,nirr_lm(i)
         write(*,*) "b",irr_nr(i,j)
        l=irr_nr(i,j)+l
     End do
     NprimUkrmol=NprimUkrmol+nLOC(i)+l
  End Do
  write(*,*) "NprimUkrmol",NprimUkrmol
!pause
end subroutine InitBasisParameters


subroutine InitIprim(iprim)
  use ModuleErrorHandling
  use ModuleBasisJUAN
  implicit none

  integer           :: l,m,i,j,k,iIrreps
  integer, external :: irr_lm
  type(BasisElementInfo),  intent(inout)  :: iprim(1:NprimUkrmol)

  Do k=1,sum(nMO)
     iprim(k)%irr=1  !need to be found. We do not need it by the moment
     iprim(k)%nr=k
     iprim(k)%MOtype=.true.
     iprim(k)%OType="MO"
  End Do
  Do l=0,Lmax
     Do j=bspline_indices(1,l),bspline_indices(2,l)
        Do m=-l,l
           iprim(k)%irr=irr_lm(l,m)  
           iprim(k)%nr=j
           iprim(k)%cr=j-bspline_indices(2,l)+2*(BsplineOrder-1)
           iprim(k)%l=l
           iprim(k)%m=m
           iprim(k)%OBtype=.true.
           iprim(k)%OType="OB"
           k=k+1
        End Do
     End Do
  End Do
  write(*,*) "primitive ukrmol basis size:",NprimUkrmol,"(it should be equal to",k-1,")"
end subroutine InitIprim


subroutine InitIelement(ielement)
  use ModuleErrorHandling
  use ModuleBasisJUAN
  use ModuleParametersJUAN
  implicit none

  integer           :: l,m,i,j,k,ia,ib,ic,iIrreps
  integer, external :: irr_lm
  type(BasisElementInfo),  intent(inout)  :: ielement(1:Nbasis)


  
Nelec=14
allocate(MO(1:Nelec))

MO(1)=1
MO(2)=2
MO(3)=3
MO(4:6)=MO(1:3)
Mo(7)=sum(irr_series(1:1))+1
MO(8)=MO(7)
MO(9)=sum(irr_series(1:2))+1!137
MO(10)=MO(9)
MO(11)=sum(irr_series(1:4))+1!203
MO(12)=sum(irr_series(1:4))+2!204
MO(13:14)=MO(11:12)

allocate(RefConf(1:Nelec))
RefConf(1)=1
RefConf(3)=2
RefConf(5)=3
RefConf(7)=2*sum(irr_total(1:1))+1!102
RefConf(9)=2*sum(irr_total(1:2))+1!137
RefConf(11)=2*sum(irr_total(1:4))+1!203
RefConf(13)=2*sum(irr_total(1:4))+2!204

RefConf(2)=1+irr_total(1)
RefConf(4)=2+irr_total(1)
RefConf(6)=3+irr_total(1)
RefConf(8)=2*sum(irr_total(1:1))+1+irr_total(2)
RefConf(10)=2*sum(irr_total(1:2))+1+irr_total(3)
RefConf(12)=2*sum(irr_total(1:4))+1+irr_total(5)
RefConf(14)=2*sum(irr_total(1:4))+2+irr_total(5)
NExtOrb=0
!write(*,*) irr_total
  l=1
  k=1
  Do i=1,8
     imin_irr(i)=k
     Do ia=1,-1,-2
   !     Do j=irr_counter(i-1)+1,irr_counter(i-1)+nMO(i)
                   Do j=sTOT(i-1)+1,sTOT(i-1)+nMO(i)
           ielement(k)%irr=i
           ielement(k)%nr=j
           ielement(k)%spin=ia
           ielement(k)%MOtype=.true.
           ielement(k)%OType="MO"
           If(any(RefConf==k))then
              ielement(k)%RCtype=.true.
              ielement(k)%RCukLoc=MO(l)
              write(*,*) k,ielement(k)%RCukLoc,l,ia
              ielement(k)%OType="RC"
              l=l+1
           End IF
           k=k+1
        End DO
        Do j=sTOT(i-1)+nMO(i)+1,sTOT(i-1)+nTOT(i)
           ielement(k)%irr=i
           ielement(k)%nr=j
           ielement(k)%spin=ia
           ielement(k)%OBtype=.true.
           ielement(k)%OType="OB"
           k=k+1
        End Do
        Do ib=1,nirr_lm(i)
           Do ic=NfirstBspl,NlastBspl
              ielement(k)%nr=ic
              ielement(k)%n_lm=ib
              ielement(k)%l=irr_l(i,ib)
              !  write(*,*) ic,NExtBspl,NBspl,irr_l(i,ib),irr_m(i,ib)
              ielement(k)%m=irr_m(i,ib)
              ielement(k)%irr=i
              ielement(k)%spin=ia
              ielement(k)%EBtype=.true.
              ielement(k)%OType="EB"
              k=k+1
              NExtOrb=NExtOrb+1
           End Do
        End Do
     End Do
     imax_irr(i)=k-1
  ENd Do
  write(*,*) "Extended Bsplines basis size Nbasis:",Nbasis,"(it should be equal to",k-1,")","Reference Conf:",l-1
  write(*,*) "Number of added Orbitals:",NExtOrb
  Nbasis=k-1 !this redefinition is for considering a unique splin and irr, etc.
  NRefConf=l-1
  
end subroutine InitIelement

subroutine STEX_BASIS(ielement)
  use ModuleErrorHandling
  use ModuleBasisJUAN
  use ModuleParametersJUAN
  implicit none

  integer           :: l,m,i,j,k,ia,ib,ic,iIrreps
  integer, external :: irr_lm
  type(BasisElementInfo),  intent(in)  :: ielement(1:Nbasis)
!Here we define the single excitation states in terms of the possitions defined forin ielements
Norb=NRefConf!l-1!Nelec!
Ndim=1+Norb*(Nbasis-Norb)
allocate(iha(1:Ndim),iep(1:Ndim))
iha=0
iep=0

k=1
! iha(k)=RefConf(1)
! iep(k)=iha(k)

k=2
l=0
Do i=1,Nbasis
   If(ielement(i)%RCtype)Then
      l=l+1
      Do j=1,Nbasis
         If(.not.ielement(j)%RCtype)Then
            If(ielement(i)%spin.eq.ielement(j)%spin)Then
               ! If(.not.ielement(j)%EBtype)Then
               iha(k)=i
               iep(k)=j
               k=k+1
               ! End IF
            End If
         End IF
      End Do
   End if
End Do
write(*,*) "STEX basis size Ndim:",Ndim,"(it should be",k-1,")","Reference Conf:",l
Ndim=k-1

end subroutine STEX_BASIS

subroutine Get_Xi(Bsplinex,ielement,iprim)
  use ModuleErrorHandling
  use ModuleBasisJUAN
  use ModuleBspline
  use ModuleSystemUtils
  implicit none

  type(ClassBspline),      intent(in)  :: Bsplinex
  type(BasisElementInfo),  intent(in)  :: ielement(1:Nbasis)
  type(BasisElementInfo),  intent(in)  :: iprim(1:NprimUkrmol)
  integer           :: ia,ib,ic,id,i,j,k,l
  real(kind(1d0))   :: ivalue,ivalue2
  procedure(D2DFun), pointer :: fptr
  real(kind(1d0))            :: parvec(2) = 0.d0
  real(kind(1d0)), external  :: GET_INTEGRALS_E
  
  fptr => Power
  parvec(1)=1.d0
  allocate(PrjMat(1:NExtBspl,1:NExtBspl),Over(1:NExtBspl,1:NExtBspl),En(1:NExtBspl))
  PrjMat=0.d0
  Over=0.d0
  En=0.d0
  n_bs=0
  Do i=NfirstBspl,NlastBspl
     ia=i-NfirstBspl+1
     Do j=NfirstBspl,NlastBspl
        ib=j-NfirstBspl+1
        ivalue2=0.d0
        Do k=BsplineOrder,2*(BsplineOrder-1)
           ivalue=NB(k)*NB(k)*NB(i)*NB(j)
           ivalue2=ivalue2+ivalue*Bsplinex.Integral(fptr,j,k,0,0,parvec=parvec)*Bsplinex.Integral(fptr,k,i,0,0,parvec=parvec)
        End Do
        PrjMat(ia,ib)=ivalue2
        Over(ia,ib)=NB(i)*NB(j)*Bsplinex.Integral(fptr,j,i,0,0)
     End Do
  End Do
  LWORK=NExtBspl*3-1
  allocate(WORK(1:LWORK))
  WORK=0.d0
  INFO=0
  call dsygv(1,'V','U',NExtBspl,PrjMat,NExtBspl,Over,NExtBspl,En,WORK,LWORK,INFO)
  allocate(brr(1:NExtBspl),arr(1:NExtBspl))
  arr=abs(En)*1.d0
  Do i=1,NExtBspl
     brr(i)=i
  End Do
  write(*,*) "aca1"
  call sort_0(NExtBspl,arr,brr)

Over=PrjMat
Do i=1,NExtBspl
    write(*,*) i,arr(i),En(brr(NExtBspl-i+1))
    PrjMat(:,i)=Over(:,brr(NExtBspl-i+1))
End Do

write(*,*) "INFO",INFO,NExtBspl
deallocate(Over,En,WORK,arr,brr)


! !Write functions
!   GridFile="OBExt000.dat"
!   Do i=1,NExtBspl!no_bsplines
!      write(GridFile(6:8),'(I3.3)') i
!      open(newunit = uid     , &
!           file    = GridFile)
!      dr=0.1d0
!      Do r=dr,Rmax,dr
!         ivalue=0.d0
!         Do j=1,NExtBspl
!            ia=NfirstBspl+j-1
!            ivalue=ivalue+NB(ia)*Bsplinex%EvalBS(r,ia)*PrjMat(i,j)
!         End Do
!         write(uid,*) r,ivalue*1.0
!      End do
!      close(uid)
!   End do
  
!!Check overlap with the internal basis
! Do i=1,NExtBspl
!    ivalue=0.d0
!    Do k=BsplineOrder,2*(BsplineOrder-1)
!         ivalue2=0.d0
!         Do j=1,NExtBspl
!            ia=NfirstBspl+j-1
!            ivalue2=ivalue2+NB(k)*NB(ia)*Bsplinex.Integral(fptr,ia,k,0,0)*PrjMat(j,i)
!         End Do
!        ! write(*,*)i, ivalue2
!         ivalue=ivalue+abs(ivalue2)
!      End Do
!      write(*,*) "------"
!      write(*,*) "ivalue",i,k,ivalue
!     ! pause
!  End Do
 
! !Check orthonormality
! Do i=1,NExtBspl
!    Do j=1,NExtBspl
!       ivalue2=0.d0
!       Do k=1,NExtBspl
!          Do l=1,NExtBspl
!             ib=NfirstBspl+l-1
!             ia=NfirstBspl+k-1
!             ivalue2=ivalue2+NB(ib)*NB(ia)*Bsplinex.Integral(fptr,ia,ib,0,0)*PrjMat(k,i)*PrjMat(l,j)
!          End Do
!       End Do
!       write(*,*)i,j,ivalue2
!    End Do
! End Do


allocate(Xi(1:Nbasis,1:Nbasis),Xip(1:Nbasis,1:Nbasis))
!Put the coefficient in vectors of the size Nbasis, containing the information
  Xi=0.d0
  Xip=0.d0
  Do j=1,Nbasis  !Do in element Xi
     Do i=1,Nbasis ! Do in components
        If(ielement(j)%EBtype.and.ielement(i)%EBtype)Then
           If(ielement(j)%spin.eq.ielement(i)%spin)Then
              If((ielement(j)%l.eq.ielement(i)%l).and.(ielement(j)%m.eq.ielement(i)%m))Then
                 ic=ielement(i)%nr
                 id=ielement(j)%nr
                 ia=ic-NfirstBspl+1
                 ib=id-NfirstBspl+1
                 Xi(i,j)=PrjMat(ia,ib)
              End If
           End If
        End If
     End Do
  End Do
  Xip=Xi



    !Put the overlap matrix the Ukrmol+EB basis in order to make the projections as matrix multiplication
  write(*,*) "-----"
  n_bs=0
  allocate(Spq(1:Nbasis,1:Nbasis))
  Spq=0.d0
  Do j=1,Nbasis
     Do i=1,Nbasis
        Spq(j,i)=GET_INTEGRALS_E(j,i,0,0,1,Bsplinex,ielement,iprim)
     End Do
  End Do
  write(*,*) "-----"

  !Gram-Schmith between the states and the internal ukrmol:
  allocate(arr(1:Nbasis))
  arr=0.d0
  Do j=1,Nbasis
     If(ielement(j)%EBtype)Then
        If(ielement(j)%nr.lt.NfirstBspl+(BsplineOrder-1))then
           !project out the ukrmol basis
           Do i=1,Nbasis
              If(ielement(i)%OBtype)Then
                 If(ielement(j)%spin.eq.ielement(i)%spin)Then
                    If(ielement(j)%irr.eq.ielement(i)%irr)Then
                       arr(:)=matmul(Spq,Xi(:,j))
                       Xi(i,j)=Xi(i,j)-arr(i)
                    End If
                 End If
              End If
           End Do
           !project out the j-1 already orthonormalized vectors
           Do i=1,j-1
              If(ielement(i)%EBtype)Then
                 If(ielement(i)%nr.lt.NfirstBspl+(BsplineOrder-1))then
                    If(ielement(j)%spin.eq.ielement(i)%spin)Then
                       If((ielement(j)%l.eq.ielement(i)%l).and.(ielement(j)%m.eq.ielement(i)%m))Then
                          ivalue=sum(Xip(:,j)*matmul(Spq,Xi(:,i)))
                          Xi(:,j)=Xi(:,j)-ivalue*Xi(:,i)
                       End If
                    End If
                 End If
              End If
           End Do
           !Normalize
           ivalue=sum(Xip(:,j)*matmul(Spq,Xi(:,j)))
           Xi(:,j)=Xi(:,j)/sqrt(ivalue)
        End If
     End If
  End Do
  deallocate(arr)
  !Let's put the Ukrmol basis in the same matrix:
  Do i=1,Nbasis
     If(.not.ielement(i)%EBtype)Then
        Xi(i,i)=1.d0
     End If
  End Do
  allocate(Xim1(1:Nbasis,1:Nbasis))
  Xim1=Xi
  allocate(iWORK(1:Nbasis))
  LWORK=Nbasis
  allocate(WORK(1:LWORK))
  call dgetrf(Nbasis,Nbasis,Xim1,Nbasis,iWORK,INFO)
    write(*,*) "INFO",INFO
  call dgetri(Nbasis,Xim1,Nbasis,iWORK,WORK,LWORK,INFO)
  write(*,*) "INFO",INFO,WORK(1)
  deallocate(WORK,iWORK)

 !Check orthonormality 
  ivalue2=0.d0
  do j=1,Nbasis
     do i=1,Nbasis
        ivalue=sum(Xi(:,j)*matmul(Spq,Xi(:,i)))
        ivalue2=ivalue2+abs(ivalue)
        If((abs(ivalue).gt.0.01d0).and.(abs(ivalue).lt.0.99d0))then
           write(*,*) j,i,ivalue,ielement(i)%OType,ielement(j)%OType
          pause
        End If
           ! If(abs(ivalue).gt.0.99d0)then
           ! write(*,*) j,i,ivalue,ielement(i)%EBtype
           ! End If
     End Do
  End Do

  


contains



  Pure DoublePrecision function Power(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y = x**parvec(1)
  end function Power
  
end subroutine Get_Xi





!.. This is the sort routine from numerical recipes. We keep it as we found a problem with the one in the library 
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
10   arr(i+1)=a	!Insert it.
     brr(i+1)=b
  enddo
  return
END subroutine sort_0


!.. This is based on the first routine introduced by Heman. Now it is obsolete.
! subroutine DimensionUkrmol(irr_series,irr_counter)
!   !
!   use ModuleBasisJUAN

!   implicit none
!   !
!   logical                             :: tmpexist = .false.
!   integer                             :: irrdim,indx,iexerror,uid
!   integer*8, intent(INOUT) :: irr_series(8),irr_counter(8)

!   inquire(file='tmpfile', exist=tmpexist)
!   write(*,*) "aca1",tmpexist
!   if( tmpexist ) then
!      open(newunit = uid     , &
!           file    ='tmpfile', &
!           status  ='old'    , &
!           form    ='formatted')
!      read(uid,*) irrdim
!      write(*,*) "aca2",irrdim
!      ! allocate(irr_series(1:irrdim))
!      write(*,*) "aca3",irr_series
!      read(uid,*) irr_series(:)
!      close(uid)
!   else  
!      write(*,*) "aca4"
!      open(newunit = uid, &
!           file    ="command.bash", &
!           form    ="formatted"   , &
!           status  ="unknown"     )
!      write(uid,'(A3,A57,A11,A,A13)') "awk", &
!           " '$0 ~ /Number of irreducible representations/{print $5 >", &
!           ' "tmpfile"}',"'"," ./log_file.0"
!      write(uid,'(A93,A11,A,A13)') "awk '$0 ~ /Number of molecular orbitals in each irreducible representation:/{getline;print >>",&
!           ' "tmpfile"}',"'"," ./log_file.0"
!      close(uid)

!      call execute_command_line("bash ./command.bash",exitstat=iexerror)
!      if(iexerror /= 0) stop 'error occured trying to read pattern#1 in log_file!'
!      write(*,*) "aca5"
!      inquire(file='tmpfile', exist=tmpexist)
!      if( .not. tmpexist ) stop 'awk failed to read keywords from log_file!'

!      open(newunit = uid     , &
!           file    ='tmpfile', &
!           status  ='old'    )
!      read(uid,*) irrdim
!      ! allocate(irr_series(irrdim)) 
!      read(uid,*) irr_series(:)
!      close(uid)
!      call execute_command_line("rm command.bash")

!   endif
!   write(*,*) "aca6"

!   !  allocate(irr_counter(0:irrdim+1))
!   irr_counter=0
!   do indx=1,irrdim
!      irr_counter(indx)=sum(irr_series(1:indx))
!   End do

!   write(*,*) "aca7"
!   print*,'! Number!  of irreducible representations:', irrdim

!   do indx= 1, irrdim
!      print*, 'Number of molecular orbitals in irreducible representation #',indx,":", irr_series(indx),sum(irr_series(1:indx))
!   end do

!   call execute_command_line("cp ./moints ./fort.1",exitstat=iexerror)
!   if(iexerror/= 0)  stop 'error occured trying to cp moints file to fort.1!'
! end subroutine DimensionUkrmol
