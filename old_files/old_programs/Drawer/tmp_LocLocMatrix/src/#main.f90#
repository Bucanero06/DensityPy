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
program ProgramTemplate

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleSystemUtils
  use ModuleString
  use ModuleIO
  use ModuleConstants
  use ModuleDiagonalize

  !.. 
  use mpi_mod 
  use ukrmol_interface 
  use precisn 
  !..

  implicit none

  !.. Run-time parameters
  !..
  integer                       :: nSize
  character(len=:), allocatable :: FileName


  !.. Local parameters
  logical, parameter      :: master_writes_to_stdout = .false.  ! all processe
  logical, parameter      :: allow_shared_memory     = .false.  ! not yet prop
  integer :: nfte  = 10
  integer :: nfti  =  1
  integer :: lembf =  0
  integer :: nocsf =  0
  integer :: nfta  =  6
  integer :: isymtp=  0
  integer :: iposit=  0
  logical :: qmoln    = .false.
  logical :: tmpexist = .false.

  integer :: indx
  integer :: iexerror
  integer :: irrdim
  integer :: nint1e,nint2e,ia,ib,ic,id,i,j,k,l,i1,i2,i3,i4
  integer :: ind(1), two_ind(1:2,1)
  real(kind(1d0)) ivalue,ivalue2
  character(len=132)  name2
  integer  n_AO,LWORK,INFO
  real(kind(1d0)), allocatable :: cf(:,:),Hirr(:,:),En(:)
  integer, allocatable    :: irr_series(:),irr_counter(:)
  real(kind(1d0)), allocatable     :: WORK(:)
  real(kind=wp)           :: scalem
  character(len=line_len) :: ukrmolp_header
  character(len=120)      :: name
  integer                 :: nalm,Nelec,Nbasis,Norb,Ndim
  integer, allocatable    :: MO(:),MeO(:),hA(:),xeP(:),hB(:),xeQ(:),shA(:),sxeP(:),shB(:),sxeQ(:)
  real(kind(1d0))         :: NucRep,H0d
  integer :: uid,iflag

  inquire(file='tmpfile', exist=tmpexist)
  if( tmpexist ) then
     open(newunit = uid     , &
          file    ='tmpfile', &
          status  ='old'    , &
          form    ='formatted')
     read(uid,*) irrdim
     allocate(irr_series(irrdim)) 
     read(uid,*) irr_series(:)
     close(uid)
  else  

     open(newunit = uid, &
          file    ="command.bash", &
          form    ="formatted"   , &
          status  ="unknown"     )
     write(uid,'(A3,A57,A11,A,A13)') "awk", &
          " '$0 ~ /Number of irreducible representations/{print $5 >", &
          ' "tmpfile"}',"'"," ./log_file.0"
     write(uid,'(A93,A11,A,A13)') "awk '$0 ~ /Number of molecular orbitals in each irreducible representation:/{getline;print >>",&
          ' "tmpfile"}',"'"," ./log_file.0"
     close(uid)

     call execute_command_line("bash ./command.bash",exitstat=iexerror)
     if(iexerror /= 0) stop 'error occured trying to read pattern#1 in log_file!'

     inquire(file='tmpfile', exist=tmpexist)
     if( .not. tmpexist ) stop 'awk failed to read keywords from log_file!'

     open(newunit = uid     , &
          file    ='tmpfile', &
          status  ='old'    )
     read(uid,*) irrdim
     allocate(irr_series(irrdim)) 
     read(uid,*) irr_series(:)
     close(uid)
     call execute_command_line("rm command.bash")

  endif


  allocate(irr_counter(0:irrdim+1))
  irr_counter=0
  do indx=1,irrdim
     irr_counter(indx)=irr_series(indx)
  End do


  print*,'Number of irreducible representations:', irrdim

  do indx= 1, irrdim
     print*, 'Number of molecular orbitals in irreducible representation #',indx,":", irr_series(indx),sum(irr_series(1:indx))
  end do

  call execute_command_line("cp ./moints ./fort.1",exitstat=iexerror)
  if(iexerror /= 0)  stop 'error occured trying to cp moints file to fort.1!'


  call mpi_mod_start(master_writes_to_stdout)!, allow_shared_memory)

  scalem=  1.0_cfp
  print*, 'loaded ukrmol_interface successfully'
  open(nfti,form='unformatted',access='stream')
  read(nfti) ukrmolp_header
  close(nfti)
  print*, ukrmolp_header, data_file_obj_id

  IF (ukrmolp_header .eq. data_file_obj_id) THEN 
     print*, 'UKRMol+ format found on the input integrals file'
     call READ_UKRMOLP_INTS(nfte,nfti,lembf,nint1e,nint2e,&
          nocsf,nfta,isymtp,irrdim,irr_series,iposit,scalem,name,&
          nalm,qmoln)
     !  RETURN
  ELSE
     WRITE(*,'(a,";",a)')ukrmolp_header,data_file_obj_id
     STOP "An unknown header on the integrals file."
  ENDIF
  allocate(cf(1:sum(irr_series),1:sum(irr_series)))
! allocate(cf(1:Nukrmol,1:Nukrmol))
  cf=0.d0
  call molecular_orbital_basis%get_orbital_coefficient_matrix(cf)
  call READ_UKRMOLP_PROPERTY_INTS(NFTI,3)

   Nbasis=sum(irr_series)
!     allocate(MeO(1:26))

!     Do i=1,7
!        MeO(i)=i
!     End Do
!     Do i=1,3
!        MeO(7+i)=sum(irr_series(1:1))+i
!     End Do
!     Do i=1,3
!        MeO(10+i)=sum(irr_series(1:2))+i
!     End Do
!     Do i=1,7
!        MeO(13+i)=sum(irr_series(1:4))+i
!     End Do
!         Do i=1,3
!        MeO(20+i)=sum(irr_series(1:5))+i
!     End Do
!         Do i=1,3
!        MeO(23+i)=sum(irr_series(1:6))+i
!     End Do
    

Ndim=Nbasis
  allocate(Hirr(1:Ndim,1:Ndim),En(1:Ndim))  
  Hirr=0.d0
  En=0.d0

   Do j=1,Ndim!10!Nbasis
      Do i=1,Ndim!10!Nbasis
         ! If((i.ne.47).or.(j.ne.47))then
         !     If((i.ne.64).or.(j.ne.64))then
         Hirr(i,j)=GET_INTEGRALS(i,j,0,0,0)!,cf(i,45)
         !Hirr(i,j)=GET_INTEGRALS(MeO(i),MeO(j),0,0,0)
     ! End If
     ! End If
End Do
write(*,*) j,GET_INTEGRALS(j,j,0,0,0)
End Do
  

  LWORK=3*Ndim-1
  allocate(WORK(1:LWORK))
  WORK=0.d0
 ! INFO=0
  write(*,*) "LWORK",LWORK,Ndim,size(Hirr),size(En)
  call dsyev('V','U',Ndim,Hirr,Ndim,En,WORK,LWORK,INFO)       

  write(*,*) "energies",INFO
  Do i=1,10!Ndim
     Write(*,*) i,En(i)
  End do
  
  ! write(*,*) GET_INTEGRALS(2,2,0,0,1)
  

  stop

  
  Nbasis=sum(irr_series)
  Write(*,*) "Nbasis",Nbasis
  Norb=7
  Nelec=2*Norb
  allocate(MO(1:7),MeO(1:26))
  MO(1)=1
  MO(2)=2
  MO(3)=3
  MO(4)=sum(irr_series(1:1))+1!102
  MO(5)=sum(irr_series(1:2))+1!137
  MO(6)=sum(irr_series(1:4))+1!203
  MO(7)=sum(irr_series(1:4))+2!204
  MeO=0
  MeO(1)=MO(3)+1
  MeO(2)=MO(3)+2
  MeO(3)=MO(3)+3
  MeO(4)=MO(3)+4
  MeO(5)=MO(4)+1
  MeO(6)=MO(4)+2
  MeO(7)=MO(5)+1
  MeO(8)=MO(5)+2
  MeO(9)=MO(7)+1
  MeO(10)=MO(7)+2
  MeO(11)=MO(7)+3
  MeO(12)=MO(7)+4
  MeO(13)=MO(7)+5
  MeO(14)=sum(irr_series(1:5))+1
  MeO(15)=sum(irr_series(1:5))+2
  MeO(16)=sum(irr_series(1:5))+3
  MeO(17)=sum(irr_series(1:6))+1
  MeO(18)=sum(irr_series(1:6))+2
  MeO(19)=sum(irr_series(1:6))+3

  !MO(18)=sum(irr_series(1:6))+2
  !MO(19)=sum(irr_series(1:6))+3
  ! 3 1 1 0 2 0 0  =  7
  ! 7 3 3 0 7 3 3  = 26
  Ndim=1+2*Norb*(Nbasis-Norb)
  Do i=1,26
     write(*,*)i, MeO(i)
  End Do
  write(*,*) "Ndim",Ndim,MO
  allocate(hA(Ndim),xeP(Ndim),hB(Ndim),xeQ(Ndim))
  allocate(shA(Ndim),sxeP(Ndim),shB(Ndim),sxeQ(Ndim))
  hA=0
  xeP=0
  hB=0
  xeQ=0
  shA=0
  sxeP=0
  shB=0
  sxeQ=0
  l=1
  Do ia=-1,1,2
     Do j=1,Norb
        Do i=1,Nbasis
           If(any(MO==i).eq..false.)Then
               If(any(MeO==i).eq..true.)Then  !considers only single excitations  to the molecular orbitals the molecular orbitals
 !             If(i.lt.30)then
                 l=l+1
                 hA(l)=MO(j)
                 xeP(l)=i
               !  write(*,*) l,hA(l),xeP(l)
                 shA(l)=ia
              End If
           End If
        End Do
     End Do
  End Do
  write(*,*) "Redefino Ndim",Ndim,"==>",l
  Ndim=l







  
  pause
  hB=hA
  xeQ=xeP
  sxeP=shA
  sxeQ=sxeP
  shB=shA
  allocate(Hirr(1:Ndim,1:Ndim),En(1:Ndim))  
  Hirr=0.d0
  En=0.d0
  !pause
  !Elemento primero de la diagonal: Energ'ia Hartree Fock
  NucRep=(7*7.d0/(2*0.54875d0))*(1.d0/1.8897259886)
  write(*,*) "Nuclear repulsion",NucRep
  !NucRep=23.6261352

  iflag=0
  i=1
  j=i
  ivalue=0.d0
  H0d=0.d0
  Do ia=-1,1,2
     Do k=1,Norb
        H0d=H0d+GET_INTEGRALS(MO(k),MO(k),0,0,iflag)
        ivalue=ivalue+GET_INTEGRALS(MO(k),MO(k),0,0,iflag)
        Do ib=-1,1,2
           Do l=1,Norb
              If((k.ne.l).or.(ia.ne.ib))then
                 ivalue=ivalue+0.5d0*GET_INTEGRALS(MO(k),MO(k),MO(l),MO(l),0)
                 ivalue=ivalue-0.5d0*GET_INTEGRALS(MO(k),MO(l),MO(l),MO(k),0)*kron(ia,ib)
              EndIf
           End Do
        End Do
     End Do
  End Do

  Hirr(i,j)=ivalue+NucRep
  write(*,*) "Hartree Fock Energy", Hirr(i,j),H0d


  Do j=2,Ndim
     ivalue=0.d0
    ivalue=ivalue+GET_INTEGRALS(hA(j),xeP(j),0,0,iflag)

     !     write(*,*) "chiches",GET_INTEGRALS(hA(j),xeP(j),0,0,0)
     Do ia=-1,1,2
        Do k=1,Norb
           If((MO(k).ne.hA(j)).or.(ia.ne.shA(j)))then
               ivalue=ivalue+GET_INTEGRALS(MO(k),MO(k),hA(j),xeP(j),0)
               ivalue=ivalue-GET_INTEGRALS(MO(k),hA(j),MO(k),xeP(j),0)*kron(ia,shA(j))*kron(ia,sxeP(j))

           End IF
        End Do
     End Do
     Hirr(i,j)=ivalue
     Hirr(j,i)=ivalue
     !   write(*,*) Hirr(i,j)
  End DO
  
  ! Do i=1,Ndim
  !    write(*,*) i,Hirr(1,i),hA(i),xeP(i),shA(i),sxeP(i)
  !    End Do
! write(*,*) Hirr(1,:)
! pause

  !j=1
  !Do i1=1,7
  !Do i2=1,i1
  !write(*,*) j,i2+i1*(i1-1)/2
  !j=j+1
  !End Do
  !End Do
  !i1=1
  !i2=4
  !i3=3
  !i4=4
  !ivalue=GET_INTEGRAL(i1,i2,i3,i4,0)
  !ivalue=GET_INTEGRAL(i1,i2,i4,i3,0)
  !ivalue=GET_INTEGRAL(i2,i1,i3,i4,0)
  !ivalue=GET_INTEGRAL(i2,i1,i4,i3,0)
  !ivalue=GET_INTEGRAL(i3,i4,i1,i2,0)
  !ivalue=GET_INTEGRAL(i3,i4,i2,i1,0)
  !ivalue=GET_INTEGRAL(i4,i3,i2,i1,0)
  !ivalue=GET_INTEGRAL(i4,i3,i1,i2,0)
  !write(*,*) "----"
  !ivalue=GET_INTEGRAL(i1,i3,i2,i4,0)
  !ivalue=GET_INTEGRAL(i1,i3,i4,i2,0)
  !ivalue=GET_INTEGRAL(i3,i1,i2,i4,0)
  !ivalue=GET_INTEGRAL(i3,i1,i4,i2,0)
  !ivalue=GET_INTEGRAL(i4,i2,i3,i1,0)
  !ivalue=GET_INTEGRAL(i4,i2,i1,i3,0)
  !ivalue=GET_INTEGRAL(i2,i4,i3,i1,0)
  !ivalue=GET_INTEGRAL(i2,i4,i1,i3,0)
  !write(*,*) "sdfghj"
  !ivalue=GET_INTEGRAL(i1,i4,i3,i2,0)
  !ivalue=GET_INTEGRAL(i1,i4,i2,i3,0)
  !ivalue=GET_INTEGRAL(i2,i3,i1,i4,0)
  !ivalue=GET_INTEGRAL(i2,i3,i4,i1,0)
  !ivalue=GET_INTEGRAL(i3,i2,i1,i4,0)
  !ivalue=GET_INTEGRAL(i3,i2,i4,i1,0)
  !ivalue=GET_INTEGRAL(i4,i1,i3,i2,0)
  !ivalue=GET_INTEGRAL(i4,i1,i2,i3,0)
  !pause



  Do i=2,Ndim
     Do j=2,Ndim
        !write(*,*) i,j
        ivalue=0.d0
        If((xeP(j).eq.xeQ(i)).and.(sxeP(j).eq.sxeQ(i)))then     
           If((hA(j).eq.hB(i)).and.(shA(j).eq.shB(i)))then
             ivalue=ivalue+H0d+NucRep
           End If
           ivalue=ivalue-GET_INTEGRALS(hA(j),hB(i),0,0,iflag)*kron(shA(j),shB(i))
           
        End IF
        If((hA(j).eq.hB(i)).and.(shA(j).eq.shB(i)))then
         !  write( *,*) xeP(j),xeQ(i)
           ivalue=ivalue+GET_INTEGRALS(xeP(j),xeQ(i),0,0,iflag)*kron(sxeP(j),sxeQ(i))
        End If
        ! !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
        If((xeP(j).ne.xeQ(i)).or.(sxeP(j).ne.sxeQ(i)))then
           If((hA(j).ne.hB(i)).or.(shA(j).ne.shB(i)))then
               ivalue=ivalue+GET_INTEGRALS(hB(i),xeQ(i),hA(j),xeP(j),0)*kron(shA(j),sxeP(j))*kron(shB(i),sxeQ(i))
               ivalue=ivalue-GET_INTEGRALS(hB(i),hA(j),xeQ(i),xeP(j),0)*kron(shA(j),shB(i))*kron(sxeP(j),sxeQ(i))
           End IF
        End IF
        ! !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
        If((xeP(j).ne.xeQ(i)).or.(sxeP(j).ne.sxeQ(i)))then
           If((hA(j).eq.hB(i)).and.(shA(j).eq.shB(i)))then
              Do ia=-1,1,2
                 Do k=1,Norb
                    If((ia.ne.shA(j)).or.(MO(k).ne.hA(j)))then
                        ivalue=ivalue+GET_INTEGRALS(MO(k),MO(k),xeQ(i),xeP(j),0)*kron(sxeP(j),sxeQ(i))
                        ivalue=ivalue-GET_INTEGRALS(MO(k),xeQ(i),MO(k),xeP(j),0)*kron(sxeP(j),ia)*kron(ia,sxeQ(i))
                    End IF
                 End Do
              End Do
           End IF
        End IF
        ! !33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
        If((xeP(j).eq.xeQ(i)).and.(sxeP(j).eq.sxeQ(i)))then
           If((hA(j).ne.hB(i)).or.(shA(j).ne.shB(i)))then
              Do ia=-1,1,2
                 Do k=1,Norb
                    If((ia.ne.shA(j)).or.(MO(k).ne.hA(j)))then
                       If((ia.ne.shB(i)).or.(MO(k).ne.hB(i)))then
                           ivalue=ivalue+GET_INTEGRALS(MO(k),hB(i),MO(k),hA(j),0)*kron(shB(i),ia)*kron(shA(j),ia)
                           ivalue=ivalue-GET_INTEGRALS(MO(k),MO(k),hB(i),hA(j),0)*kron(shB(i),shA(j))
                       End If
                    End IF
                 End Do
              End Do
               ivalue=ivalue+GET_INTEGRALS(xeP(j),hB(i),xeP(j),hA(j),0)*kron(shB(i),sxeP(j))*kron(shA(j),sxeP(j))
               ivalue=ivalue-GET_INTEGRALS(xeP(j),xeP(j),hB(i),hA(j),0)*kron(shB(i),shA(j))
           End IF
        End IF
        ! !444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
         If((xeP(j).eq.xeQ(i)).and.(sxeP(j).eq.sxeQ(i)))then
            If((hA(j).eq.hB(i)).and.(shA(j).eq.shB(i)))then
              Do ia=-1,1,2
                 Do k=1,Norb
                    Do ib=-1,1,2
                       Do l=1,Norb
                          If((ia.ne.shA(j)).or.(MO(k).ne.hA(j)))then
                             If((ib.ne.shA(j)).or.(MO(l).ne.hA(j)))then
                                If((ia.ne.ib).or.(MO(k).ne.MO(l)))then
                                    ivalue=ivalue+0.5d0*GET_INTEGRALS(MO(k),MO(k),MO(l),MO(l),0)
                                    ivalue=ivalue-0.5d0*GET_INTEGRALS(MO(k),MO(l),MO(k),MO(l),0)*kron(ia,ib)
                                End If
                             End If
                          End IF
                       End Do
                    End Do
                    If((ia.ne.shA(j)).or.(MO(k).ne.hA(j)))then
                        ivalue=ivalue+GET_INTEGRALS(MO(k),MO(k),xeP(j),xeP(j),0)
                        ivalue=ivalue-GET_INTEGRALS(MO(k),xeP(j),MO(k),xeP(j),0)*kron(ia,sxeP(j))
                    End IF
                 End Do
              End Do
               !ivalue=ivalue+NucRep
            End IF
         End IF

         
        Hirr(i,j)=ivalue
        !      Hirr(j,i)=ivalue
     End DO
     !Hirr(i,i)=Hirr(i,i)+NucRep
     !pause
  End DO

  ivalue=0.d0
  Do i=1,Ndim
     Do j=1,Ndim
        ivalue=ivalue+abs(Hirr(i,j)-Hirr(j,i))*1.d0
        !write(*,*)ivalue,i,j, Hirr(i,j)-Hirr(j,i)
     End Do
!write(*,*) i,Hirr(i,i)
     !pause
  End Do

  write(*,*) "assimetry",ivalue
  !pause

  !call molecular_orbital_basis%get_orbital_coefficient_matrix(cf)

  !call molecular_orbital_basis%print_orbitals

  !Do i+1,32i
  !write(*,*) sum(cf(:,14)*cf(:,15))

  !Do ib=1,sum(irr_series)
  !Do ia=1,40
  !write(*,*)ib, ia,cf(ia,ib)
  !End Do
  !pause
  !End DO


  !write(*,*) "paso",Hirr

  LWORK=3*Ndim-1
  allocate(WORK(1:LWORK))
  WORK=0.d0
  INFO=0
  write(*,*) "LWORK",LWORK,Ndim,size(Hirr),size(En)
  call dsyev('V','U',Ndim,Hirr,Ndim,En,WORK,LWORK,INFO)       

  write(*,*) "energies",INFO
  Do i=1,10!Ndim
     Write(*,*) i,En(i)
  End do

  ! open(unit=2, file="energ_lucia.dat")
  ! Do i=1,Ndim
  !    read(2,*) ivalue
  !    write(*,*) ivalue-En(i)
  ! End Do
  ! close(2)

  call mpi_mod_finalize


!!$  call GetRunTimeParameters( FileName, nSize )

  stop

contains

  function kron(i,j) result(k)
    integer, intent(in) :: i,j ! input
    integer             :: k ! output
    k=0
    If(i.eq.j)then
       k=1
    End IF
  end function kron  

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

