module ModuleCSFLists

  private
    
  public :: MakeCSFList
  public :: AugmentCSFList

  character(len=*), parameter :: lstrnUC="SPDFGHIKLMNOQRTUVWXYZ"
  character(len=*), parameter :: lstrn="spdfghiklmnoqrtuvwxyz"
  character(len=*), parameter :: pstrn="eo"

contains

  integer function mchar2m( mchar ) result( m )
    character, intent(in) :: mchar
    read(mchar,*)m
  end function mchar2m

  integer function lchar2l( lchar ) result( l )
    use ModuleString
    character, intent(in) :: lchar
    character :: lchar_
    lchar_=lchar
    call SetStringToUpperCase(lchar_)
    l=index(lstrnUC,lchar_)-1
  end function lchar2l

  integer function pchar2p( pchar ) result( p )
    character, intent(in) :: pchar
    p=index(pstrn,pchar)-1
  end function pchar2p

  character function m2mchar( m ) result( mchar )
    integer, intent(in) :: m
    character(len=8) :: strnBuf
    write(strnBuf,"(i0)")m
    strnBuf=adjustl(strnBuf)
    mchar=strnBuf(1:1)
  end function m2mchar

  character function l2lchar( l ) result( lchar )
    integer, intent(in) :: l
    lchar=lstrn(l+1:l+1)
  end function l2lchar

  character function l2lcharUC( l ) result( lchar )
    integer, intent(in) :: l
    lchar=lstrnUC(l+1:l+1)
  end function l2lcharUC

  character function p2pchar( p ) result( pchar )
    integer, intent(in) :: p
    pchar=pstrn(p+1:p+1)
  end function p2pchar

  !.. returns the multiplicity, angular momentum, and parity
  !   of a given term
  subroutine ExtractTermQN( term, m, l, p )
    implicit none
    character(len=*), intent(in) :: term
    integer         , intent(out):: m, l, p
    m=mchar2m(term(1:1))
    l=lchar2l(term(2:2))
    p=pchar2p(term(3:3))
  end subroutine ExtractTermQN


  !.. Starting from the list of N-electron configurations for an ionic 
  !   term IonTerm, contained in the file FileIon, generate a list of 
  !   states with N+1 electrons and term AugTerm, by augmenting the ion 
  !   configuration with an electron outside the active space, with 
  !   maximum orbital angular momentum lwmax. The result is saved to FileCC.
  !..
  subroutine AugmentCSFList( FileIon, IonTerm, lwmax, AugTerm, FileCC )

    use ModuleErrorHandling
    use ModuleString
    use ModuleIO

    implicit none

    character(len=*), parameter :: lstrn="spdfghiklmnoqrtuvwxyz"
    character(len=*), parameter :: pstrn="eo"

    character(len=*), intent(in) :: FileIon
    character(len=*), intent(in) :: IonTerm
    integer         , intent(in) :: lwmax
    character(len=*), intent(in) :: AugTerm
    character(len=*), intent(in) :: FileCC

    integer :: ion_uid, aug_uid, ncsfion, iIon
    character(len=512) :: line1, line2
    character          :: firstChar
    integer :: mion, lion, pion
    integer :: maug, laug, paug
    integer :: lw, lwmi, lwma, nlw
    integer :: nsubShells, icharShellma

    !.. Check compatibility of the two terms
    !..
    call ExtractTermQN( IonTerm, mion, lion, pion )
    call ExtractTermQN( AugTerm, maug, laug, paug )
    !
    !
    !.. Determines the range of orbital angular momenta 
    !   that are compatible with the parent-ion and 
    !   augemented-state angular momenta and parity
    lwmi = abs( lion - laug )
    lwmi = lwmi + mod( pion + lwmi + paug, 2 )
    lwma = min( lwmax, lion + laug )
    lwma = lwma - mod( lwmi + lwma, 2 )
    nlw  = ( lwma - lwmi ) / 2 + 1

    if( lwma < lwmi ) nlw=0 

    !.. Check if their multiplicity is compatible
    if( abs( mion - maug ) /= 1 ) nlw=0 

    call openFile(FileIon,ion_uid,"read" ,"formatted")
    call openFile(FileCC ,aug_uid,"write","formatted")
    !
    !.. The first two lines are there by default
    !   It assumes that the first line contains the number
    !   of ionic configurations
    read (ion_uid,"(a)") line1
    read (ion_uid,"(a)") line2
    read (line1,*)ncsfion
    do iIon=1,ncsfion-1
       write(aug_uid,"(2x,i0,a)",advance="no") nlw,","
    enddo
    write(aug_uid,"(2x,i0,a)") nlw,","
    write(aug_uid,"(a)") trim(line2)
    if(nlw>0)then
    do
       !.. For example, line1 and line2 could read
       !  2s( 2)  2p( 3)  3s( 1)  3p( 1)
       ! 1S0 2P1 2S1 2P1 2P  3P  2P 
       !|_sub-sh coupl__|_sh coupl_|
       !..
       read(ion_uid,"(a)") line1
       firstChar=adjustl(line1)
       if(firstChar=="*")exit
       read(ion_uid,"(a)") line2
       !
       nsubShells = ( nTokens( line2 ) + 1 ) / 2
       icharShellma = 4 * nsubShells
       do lw = lwmi, lwma, 2
          write(aug_uid,"(a)")trim(line1)//"  ;"//l2lchar(lw)//"( 1)"
          if(len_trim(line2(icharShellma+1:))==0)then
             write(aug_uid,"(a)")&
                  trim(line2(1:icharShellma))//&
                  " 2"//l2lcharUC(lw)//"1"//&
                  " "//m2mchar(maug)//l2lcharUC(laug)
          else
             write(aug_uid,"(a)")&
                  trim(line2(1:icharShellma))//&
                  " 2"//l2lcharUC(lw)//"1"//&
                  trim(line2(icharShellma+1:))//&
                  "  "//m2mchar(maug)//l2lcharUC(laug)
          endif
       enddo
    enddo
    endif
    write(aug_uid,"(a)")"*"
    close( ion_uid )
    close( aug_uid )
    
  end subroutine AugmentCSFList


  subroutine MakeCSFList( &
       act, nmax, lmax, minnel, nmaxref, refconfig, sshkind, term, nex, io_file )
    use ModuleErrorHandling
    use ModuleIO

    implicit none

    character(len=*), parameter :: com_file="comscr"
    character(len=*), parameter :: lstrn="spdfghiklmnoqrtuvwxyz"
    character(len=*), parameter :: pstrn="eo"

    character(len=*), intent(in) :: act
    integer         , intent(in) :: nmax, lmax

    !> minnel(in) is the minimum number of electrons in
    !! subshells with principal quantum number n <= in.
    integer         , intent(in) :: minnel(:)
    integer         , intent(in) :: nmaxref
    integer         , intent(in) :: refconfig(:)
    !> For subshells in the REFERENCE: closed, inactive, active, or minimum ( c/i/a/0..<n>/ )
    !> For subshells in RAS3: inactive, active, doubled excited or virtual? (i/a/d/v)
    character       , intent(in) :: sshkind(:)
    character(len=*), intent(in) :: term
    !> number of excitations
    integer         , intent(in) :: nex
    character(len=*), intent(in) :: io_file
    integer, parameter :: NRANCHAR=16
    character(len=NRANCHAR) :: ranstrn
    character(len=32) :: nstrn

    integer   :: uid, in, il, iSubShell, ipar, iostat
    character :: slmax, RefCFGParity
    logical   :: Make_New_List, Add_to_List, io_file_exists
    logical   :: Virtual_present

    if( nmax < 1 )      call Assert("nmax < 1 in NewList")
    if( lmax < 0 )      call Assert("lmax < 0 in NewList")
    if( lmax > nmax-1 ) call Assert("lmax > nmax - 1 in NewList")
    slmax = lstrn(lmax+1:lmax+1)

    Virtual_present=.FALSE.
    do iSubShell=1,size(sshKind)
       in=index(sshKind(iSubShell),"v")
       if(in>0)then
          Virtual_present=.TRUE.
          exit
       endif
    enddo

    !.. Check that the parity of the reference configuration is right
    iPar=0 !iPar=0 => even, iPar=1 => odd
    iSubShell=0
    do in=1,nmaxref
       do il=0,in-1
          iSubShell=iSubShell+1
          iPar = mod( iPar + il * refconfig(iSubShell), 2 )
       enddo
    enddo
    RefCFGParity = pstrn(iPar+1:iPar+1)
    if( term(3:3) /= RefCFGParity )return

    Make_New_List =  trim(act) == "new"
    Add_to_List   =  trim(act) == "add"

    !.. Check if the input file exists. If not, set Make_New_List to true
    inquire( file=io_file, exist = io_file_exists )
    if( io_file_exists .and. Make_New_List )then
       call System( "rm    " // io_file )
       call System( "touch " // io_file )
    endif
    if( .not. io_file_exists )then
       Make_New_List = .TRUE.
       Add_to_List   = .FALSE.
       call System( "touch " // io_file )
    endif

    call openFile(com_file,uid,"write","formatted")
    if( Make_New_List )then
       write(uid,*)
    elseif( Add_to_List )then
       write(uid,"(a)") "a"
    endif
    write(uid,*) !MCHF (instead of Breit)
    write(uid,*) !Default Symmetry
    if( Add_to_List ) write(uid,*) !Use new reference set
    write(uid,"(i0)") nmax  !highest principal quantum number
    write(uid,"(a)") slmax !highest orbital angular momentum
    write(uid,"(a)") "n"
    write(uid,"(a)") "y" !limitation to the population of subshells
    do in=1,nmax-1
       write(uid,"(i0)") minnel(in) !minimum number of electrons for n <= in
    enddo
    write(uid,"(i0)") nmaxref !maximum n in reference config
    iSubShell=0
    do in=1,nmaxref
       do il=0,min(in-1,lmax)
          iSubShell=iSubShell+1
          if( Add_to_List .and. sshKind(iSubShell) == "c" )cycle 
          write(uid,"(i0)") refconfig(iSubShell)
          write(uid,"(a)" ) sshKind(iSubShell)
       enddo
    enddo
    do in=nmaxref+1,nmax
       do il=0,min(in-1,lmax)
          iSubShell=iSubShell+1
          write(uid,"(a)") sshKind(iSubShell)
       enddo
    enddo
    if(Make_New_List) write(uid,"(a)")trim(adjustl(term(1:2)))
    write(uid,"(i0)") nex
    if(virtual_present) write(uid,"(i0)") 0 ! no excitations to the virtual set
    write(uid,"(i0)") !Do not generate another list
    close(uid)
    
    call sublsgen( com_file, io_file, io_file )
    
!!$    !***
!!$    STOP
!!$    !***

    !.. add header to the list, with the total number of configurations
    !
    call OpenFile( io_file, uid, "read", "formatted" )
    if(uid/=0)then
       in=0
       do
          read(uid,*,iostat=iostat)
          if(iostat/=0)exit
          in=in+1
       enddo
       close(uid)
       in=(in-3)/2
       write(nstrn,"(i5)")in
       ranstrn = randomStrn( NRANCHAR )
       call system("echo '"//trim(nstrn)//"' > "//ranstrn)
       call system("tail -n +2 "//io_file//" >> "//ranstrn)
       call system("mv "//ranstrn//" "//io_file)
    endif

    call system("rm "//com_file)

    return
  end subroutine MakeCSFList
  

  subroutine sublsgen( command_stack_file, input_file, output_file )
    ! *
    ! *    Written by:
    ! *    Lennart Sturesson
    ! *    Department of Physics
    ! *    Lund Institute of Technology
    ! *    P.O. Box 118, S-221 00  Lund, Sweden
    ! *
    ! *    Modified in April 2017 by
    ! *    Luca Argenti (luca.argenti@ucf.edu)
    ! *    Deptartment of Physics & CREOL
    ! *    University of Central Florida
    ! *    Physics Science Building 304
    ! *    4111 Libra Drive
    ! *    32816 Orlando (FL)
    ! *    United States
    ! *

    use ModuleCom

    use ModuleIO
    use ModuleErrorHandling

    logical, parameter :: DYN = .TRUE.
    integer, parameter :: RAN_STRN_LEN=12

    character(len=*), intent(in) :: command_stack_file
    character(len=*), intent(in) :: input_file
    character(len=*), intent(in) :: output_file

    integer org(1:15,0:10),varmax,resS,resL,skal,anel,par,&
         low(1:15,0:10),posn(110),posl(110),nmax,cfmax,&
         J2min,J2max,minS,minL,virmax,lim(15),nold
    logical lock(1:15,0:10),closed(1:15,0:10),second,slut,breit,&
         virtu(1:15,0:10),vir,dubbel(1:15,0:10),extra,advexp
    character X*1,Y*1
    logical :: UseTempFile
    integer :: uid, iostat
    character(len=128) :: strnBuf

    !.. Set the communication
    ! file names in ModuleCom
    !======================

    if(allocated(CLIST_INP)) deallocate(CLIST_INP)
    if(trim(adjustl(input_file))==trim(adjustl(output_file)))then
       inquire( file=input_file, exist = UseTempFile )
       allocate(CLIST_INP,source="__tmp__"//randomStrn(RAN_STRN_LEN))
       call system("cp "//input_file//" "//CLIST_INP)
    else
       allocate(CLIST_INP,source=trim(adjustl(input_file)))
       UseTempFile=.FALSE.
    endif
    !
    if(allocated(CLIST_OUT)) deallocate(CLIST_OUT)
    allocate( CLIST_OUT, source = trim(adjustl(output_file)) )
    !
    if(allocated(CLIST_NEW)) deallocate(CLIST_NEW)
    allocate( CLIST_NEW, source = 'clist.new' )
    !
    if(allocated(CLIST_LOG)) deallocate(CLIST_LOG)
    allocate( CLIST_LOG, source = 'clist.log' )

    !======================

    call openFile( command_stack_file, inp_com_uid, "read" , "formatted" )

    if(LOG_INFO) call openFile(CLIST_LOG,log_com_uid,"write","formatted")
    if(TER_INFO) write(*,200) 'New list, add to existing list, expand existing ',&
         'list? (*/a/e/s/r)'
    read(inp_com_uid,100) X
    if(TER_INFO) write(*,*) X
    advexp = X.EQ.'e' .OR. X.EQ.'E'
    if (advexp) then
       if(TER_INFO) write(*,*)
       if(TER_INFO) write(*,*) 'This option is only running in the MCHF-mode!'
       if(TER_INFO) write(*,*)
       breit = .FALSE.
    else
       if(TER_INFO) write(*,200) 'Breit or MCHF? (B/*)'
       read(inp_com_uid,100) Y
       if(TER_INFO) write(*,*) Y
       breit = Y.EQ.'b' .OR. Y.EQ.'B'
    endif
    call Reffa(posn,posl,dyn)
    if (X.EQ.'a' .OR. X.EQ.'A' .OR. advexp) then
       call Adder(closed,slut,resS,resL,anel,par,dyn,advexp)
       if (slut) then
          if(TER_INFO) write(*,200)
          if(TER_INFO) write(*,200) 'The clist.inp-file is not readable! Run the "r"-option if the file is a copy of a cfg-file!' 
          stop
       endif
       if (.NOT. advexp) then
          if(TER_INFO) write(*,200)
          if(TER_INFO) write(*,200) 'Use new reference set or make simple expansion of the active sets in input-file? (*/e)'
          read(inp_com_uid,100) Y
          if(TER_INFO) write(*,*) Y
          extra = Y.EQ.'e' .OR. Y.EQ.'E'
          if (extra) then
             call Scan(nold)
             if (nold.LT.10) then
                if(TER_INFO) write(*,300) 'The highest identified n-number is',nold,'.'
             else
                if(TER_INFO) write(*,400) 'The highest identified n-number is',nold,'.' 
             endif
             if (nold.GT.14) then
                if(TER_INFO) write(*,200) 'The identified n-number is to high. - Highest possible n-number is 14!'
                stop
             endif
             call Expand(nold,dyn)
          endif
       else
          extra = .FALSE.
       endif
       second = .TRUE.
    else
       call Matain(org,lock,closed,resS,resL,varmax,virmax,skal,nmax,&
            anel,par,dyn,low,breit,J2min,J2max,minS,minL,virtu,&
            vir,lim,dubbel)
       call Twolines(org,closed,.TRUE.)
       call Blanda(org,varmax,virmax,lock,resS,resL,skal,nmax,dyn,&
            .TRUE.,low,posn,posl,breit,J2min,J2max,minS,minL,&
            virtu,vir,lim,dubbel)
       second = .FALSE.
       extra  = .FALSE.
    endif
    if (.NOT.extra) then
       if (advexp) then
          call Matcin(lock,closed,varmax,cfmax,nmax)
          call Blandc(varmax,cfmax,lock,resS,resL,nmax,dyn,posn,posl)
          second = .FALSE.
       else
          call Matbin(org,lock,closed,varmax,virmax,skal,second,anel,&
               par,dyn,low,breit,J2min,J2max,resS,resL,minS,minL,&
               virtu,vir,nmax,lim,dubbel)
       endif
    endif
    if (second) then
       if (.NOT.extra) then
          call Twolines(org,closed,.FALSE.)
          call Blanda(org,varmax,virmax,lock,resS,resL,skal,nmax,dyn,&
               .FALSE.,low,posn,posl,breit,J2min,J2max,minS,&
               minL,virtu,vir,lim,dubbel)
       endif
       call Merge(.FALSE.,dyn)
       if(TER_INFO) write(*,200) 'The merged file is called '//CLIST_OUT
    else
       call Merge(.TRUE.,dyn)
       if(TER_INFO) write(*,200) 'The generated file is called '//CLIST_OUT
    endif

    if(LOG_INFO)close(log_com_uid)
    call system("rm "//CLIST_NEW)
    if(UseTempFile) call system("rm "//CLIST_INP)
    
    !.. Check that the output file contains at least one valid configuration.
    !   If not, it cancels it
    !..
    call openFile(CLIST_OUT, uid, "read", "formatted")
    read(uid,*,iostat=iostat)
    read(uid,*,iostat=iostat)
    read(uid,*,iostat=iostat) strnBuf
    strnBuf=adjustl(strnBuf)
    if(strnBuf(1:1)=="*".or.iostat/=0) call system("rm "//CLIST_OUT)


100 format(A)
200 format(' ',10A)
300 format(' ',A,I2,A)
400 format(' ',A,I3,A)

  end subroutine sublsgen

end module ModuleCSFLists
