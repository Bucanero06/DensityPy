!! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
module ModuleIO

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleErrorHandling
  use ModuleString

  implicit none

  private

  character(len=*), public, parameter :: STRN_FMT="a"
  character(len=*), public, parameter :: STRN_SEQ_FMT="*(x,"//STRN_FMT//")"
  character(len=*), public, parameter :: FULL_STRN_SEQ_FMT="("//STRN_SEQ_FMT//")"

  character(len=*), public, parameter :: INT0_FMT="i0"
  character(len=*), public, parameter :: INT0_SEQ_FMT="*(x,"//INT0_FMT//")"
  character(len=*), public, parameter :: FULL_INT0_SEQ_FMT="("//INT0_SEQ_FMT//")"

  character(len=*), public, parameter :: EDBL_FMT="e24.16"
  character(len=*), public, parameter :: EDBL_SEQ_FMT="*(x,"//EDBL_FMT//")"
  character(len=*), public, parameter :: FULL_EDBL_SEQ_FMT="("//EDBL_SEQ_FMT//")"

  public :: OpenFile
  public :: GetFullText
  public :: WriteFullText
  public :: FetchGlobalVariable
  public :: FetchLine
  public :: FetchLineStrn
  public :: CheckFileMustBePresent
  public :: ReadMatrix
  public :: WriteMatrix

  interface SaveMatrix
     module procedure SaveComplexMatrix
     module procedure SaveComplexMatrix_fixedSize
     module procedure SaveComplexMatrix_nonzero
  end interface SaveMatrix
  public SaveMatrix
  interface LoadMatrix
     module procedure LoadComplexMatrix
  end interface LoadMatrix
  public LoadMatrix

  interface SaveVector
     module procedure SaveComplexVector
     module procedure SaveDbleVector
  end interface SaveVector
  public SaveVector
  interface LoadVector
     module procedure LoadComplexVector
     module procedure LoadDbleVector
  end interface LoadVector
  public LoadVector


  public :: randomStrn


contains
  
  function randomStrn( n ) result( strnBuf )
    character(len=*), parameter :: chars  ="qwertyuiopasdfghjklzxcvbnm1234567890"
    integer         , parameter :: nchars = 36
    character(len=*), parameter :: master ="0000 0000 0000 0000 0000 0000 0000 0000 "
    integer         , parameter :: nmax   = 40
    integer, intent(in) :: n
    character(len=:), allocatable :: strnBuf
    integer :: n_, ich, iran
    real(kind(1d0)) :: dran
    logical, save :: First_call = .TRUE.
    if(First_call)then
       call Random_Seed
       First_call=.FALSE.
    endif
    if(n<=0)return
    n_=min(n,nmax)
    allocate(strnBuf,source=master(1:n_))
    do ich=1,n_
       call random_number(dran)
       iran=int((nchars-1.d0-1.d-14)*dran)+1
       strnBuf(ich:ich)=chars(iran:iran)
    enddo
  end function randomStrn

  subroutine GetFullText( FileName, FullText_, Newline_Char_ )
    character(len=*),              intent(in)    :: FileName
    character(len=:), allocatable, intent(inout) :: FullText_
    character(len=*), optional   , intent(in)    :: Newline_Char_

    character, parameter          :: NEWLINE_CHAR_DEFAULT = " "
    character(len=1024*1024)      :: FullText
    character(len=1024)           :: LineStrn
    character(len=:), allocatable :: Newline_Char
    integer                       :: uid, iostat, ichar

    call OpenFile( FileName, uid, "read", "formatted" )

    if(present(Newline_Char_))then
       allocate(Newline_Char, source = Newline_Char_ )
    else
       allocate(Newline_Char, source = NEWLINE_CHAR_DEFAULT )
    endif

    FullText=" "
    lineCycle : do

       !.. Read one line from file. Exit if the file is over
       read(uid,"(a)",IOSTAT=iostat) LineStrn
       if( iostat /= 0 ) exit lineCycle

       !.. Skips comented (#) and blanck lines.
       ichar = index(LineStrn,"#")
       if ( ichar > 0 ) LineStrn(ichar:) = " "
       LineStrn = adjustl(LineStrn)
       if ( len_trim(LineStrn) == 0 ) cycle lineCycle

       !.. Accumulates the new line on the full text variable
       FullText=trim(FullText)//Newline_Char//trim(LineStrn)

    end do lineCycle
    close(uid)
    FullText=trim(FullText)//Newline_Char

    if(allocated(FullText_)) deallocate(FullText_)
    allocate(FullText_,source=trim(adjustl(FullText)))

  end subroutine GetFullText


  subroutine WriteFullText( FileName, FullText, NewlineChar )
    character(len=*), intent(in)  :: FileName
    character(len=*), intent(in)  :: FullText
    character(len=*), intent(in)  :: NewlineChar
    integer                       :: uid, ichar, ichar1
    call OpenFile( FileName, uid, "write", "formatted" )
    ichar1=1
    do
       ichar = index(FullText,NewLineChar)
       if(ichar<1)then
          write(uid,"(a)")FullText(ichar1:)
          exit
       endif
       write(uid,"(a)")FullText(ichar1:ichar-1)
       ichar1=ichar+len_trim(NewLineChar)
    end do
    close(uid)
  end subroutine WriteFullText


  subroutine FetchGlobalVariable( Text, VariableName, VariableLabel, EOLN_STRN )
    character(len=*),              intent(in)  :: Text
    character(len=*),              intent(in)  :: VariableName
    character(len=:), allocatable, intent(out) :: VariableLabel
    character(len=*), optional   , intent(in)  :: EOLN_STRN
    !
    character(len=512) :: strn
    integer :: ichar, icharEnd
    !
    ichar = index( Text, VariableName )
    if(ichar<1) return
    !
    ichar = ichar + index( Text(ichar+1:),"=" )
!!$    write(*,"(a,x,i0,x,a)") VariableName, ichar, trim(Text(ichar:))
    if(present(EOLN_STRN))then
       iCharEnd=ichar + index( Text(ichar+1:), EOLN_STRN ) - 1
       allocate( VariableLabel, source = trim(adjustl(Text(ichar+1:iCharEnd))))
    else
       allocate( VariableLabel, source = trim(adjustl(Text(ichar+1:))))
    endif
    !
  end subroutine FetchGlobalVariable


  subroutine FetchLineStrn( uid, LineStrn, ComCh_ )
    integer                   , intent(in)  :: uid
    character(len=*)          , intent(out) :: LineStrn
    character(len=*), optional, intent(in)  :: ComCh_
    !
    integer             :: ichar, iostat
    character(len=:), allocatable :: ComCh
    !
    if(present(ComCh_))then
       allocate(ComCh,source=ComCh_)
    else
       allocate(ComCh,source="#")
    endif
    lineCycle : do

       read(uid,"(a)",IOSTAT=iostat) LineStrn
       if( iostat /= 0 )then
          LineStrn=" "
          exit lineCycle
       endif

       !.. Skips comented (#) and blanck lines.
       ichar = index(LineStrn,ComCh)
       if ( ichar > 0 ) LineStrn(ichar:) = " "
       LineStrn = adjustl(LineStrn)
       if ( len_trim(LineStrn) == 0 )then
          cycle lineCycle
       else
          exit lineCycle
       endif
    end do lineCycle
    
  end subroutine FetchLineStrn


  subroutine FetchLine( uid, Line )
    integer                      , intent(in)  :: uid
    character(len=:), allocatable, intent(out) :: Line
    character(len=1024) :: LineStrn
    call FetchLineStrn(uid,LineStrn)
    if( len_trim(LineStrn) == 0 )then
       allocate(Line,source=" ")
    else
       allocate(Line,source=trim(adjustl(LineStrn)))
    endif
  end subroutine FetchLine


!!$  subroutine FetchLine( uid, Line )
!!$    integer                      , intent(in)  :: uid
!!$    character(len=:), allocatable, intent(out) :: Line
!!$    !
!!$    character(len=1024) :: LineStrn
!!$    integer             :: ichar, iostat
!!$    !
!!$    lineCycle : do
!!$
!!$       read(uid,"(a)",IOSTAT=iostat) LineStrn
!!$       if( iostat /= 0 )then
!!$          LineStrn=" "
!!$          exit lineCycle
!!$       endif
!!$
!!$       !.. Skips comented (#) and blanck lines.
!!$       ichar = index(LineStrn,"#")
!!$       if ( ichar > 0 ) LineStrn(ichar:) = " "
!!$       LineStrn = adjustl(LineStrn)
!!$       if ( len_trim(LineStrn) == 0 )then
!!$          cycle lineCycle
!!$       else
!!$          exit lineCycle
!!$       endif
!!$    end do lineCycle
!!$    
!!$    if( len_trim(LineStrn) == 0 )then
!!$       allocate(Line,source=" ")
!!$    else
!!$       allocate(Line,source=trim(adjustl(LineStrn)))
!!$    endif
!!$  end subroutine FetchLine


  subroutine OpenFile( FileName, uid, Purpose, Form )
    !> The name of the file to be opened.
    character(len=*),           intent(in)    :: FileName
    !> Unit number associated to the opened file.
    integer,                    intent(out)   :: uid
    !> Can be "write" or "read" depending on the action needed.
    character(len=*), optional, intent(in)    :: Purpose
    !> Can be "formatted" or "unformatted".
    character(len=*), optional, intent(in)    :: Form
    !
    logical :: opened
    integer :: iostat
    character(len=IOMSG_LENGTH) :: iomsg
    character(len=:), allocatable :: NewPurpose, NewForm
    !
    if ( .not.present(Purpose) ) then
       allocate( NewPurpose, source = "read" )
    else
       allocate( NewPurpose, source = Purpose )
    end if
    if ( .not.present(Form) ) then
       allocate( NewForm, source = "unformatted" )
    else
       allocate( NewForm, source = Form )
    end if
    !
    OPEN(NewUnit =  uid,          &
         File    =  FileName,     &
         Form    =  NewForm,      &
         Status  = "unknown",     &
         Action  = NewPurpose,    &
         iostat  = iostat ,       &
         iomsg   = iomsg    )
    if( iostat /= 0 )then
       !call ErrorMessage(iomsg)
       uid=0
       return
    endif

    INQUIRE( unit = uid, &
         iomsg = iomsg, &
         iostat = iostat, &
         opened = opened )
    if ( iostat /= 0 ) call Assert(iomsg)

  end subroutine OpenFile


  subroutine SaveComplexMatrix_nonzero( FileName, A, form, THRESHOLD )
    character(len=*)  , intent(in) :: FileName, form
    complex(kind(1d0)), intent(in) :: A(:,:)
    real(kind(1d0))   , intent(in) :: THRESHOLD
    integer :: uid, i, j, ni, nj
    call OpenFile( FileName, uid, "write", form )
    ni=size(A,1)
    nj=size(A,2)
    
    if( form .is. "formatted" )then
       write(uid,"(a)") "Z"
       write(uid,"(2(x,i0),x,d24.16)") ni, nj,THRESHOLD
       do i = 1, ni
          do j = 1, nj
             if( abs(A(i,j)) > THRESHOLD ) write( uid, "(2(x,i0),2(x,d24.16))" ) i,j,dble(A(i,j)),aimag(A(i,j))
          enddo
       enddo
    elseif( form .is. "unformatted" )then
       write(uid) "Z"
       write(uid) ni, nj, THRESHOLD
       do i=1,ni
          do j=1,nj
             write(uid) i,j,A(i,j)
          enddo
       enddo
    endif
    close(uid)
  end subroutine SaveComplexMatrix_nonzero

  subroutine SaveComplexMatrix( FileName, A, form )
    character(len=*)  , intent(in) :: FileName, form
    complex(kind(1d0)), intent(in) :: A(:,:)
    integer :: uid, i, j, ni, nj
    call OpenFile( FileName, uid, "write", form )
    ni=size(A,1)
    nj=size(A,2)
    if( form .is. "formatted" )then
       write(uid,"(a)") "C"
       write(uid,"(2(x,i0))") ni, nj
       do i = 1, ni
          write( uid, "(*(x,d24.16))" ) ( dble(A(i,j)), j = 1, nj )
          write( uid, "(*(x,d24.16))" ) (aimag(A(i,j)), j = 1, nj )
       enddo
    elseif( form .is. "unformatted" )then
       if( sum(abs(A(1:ni,1:nj))) == 0.d0 )then
          write(uid) "0"
          write(uid) ni, nj
       else
          write(uid) "C"
          write(uid) ni, nj
          write(uid) ( ( A(i,j), i = 1, ni ), j = 1, nj )
       endif
    endif
    close(uid)
  end subroutine SaveComplexMatrix

  subroutine SaveComplexMatrix_fixedsize( FileName, A, di, dj, form )
    character(len=*)  , intent(in) :: FileName, form
    complex(kind(1d0)), intent(in) :: A(:,:)
    integer           , intent(in) :: di, dj
    integer :: uid, i, j, ni, nj
    call OpenFile( FileName, uid, "write", form )
    ni=size(A,1)
    nj=size(A,2)
    if(ni<di.or.nj<dj)then
       write(ERROR_UNIT,"(a)") "Size inconsistency in SaveComplexMatrix_fixedsize"
       stop
    endif
    if( form .is. "formatted" )then
       write(uid,"(a)") "C"
       write(uid,"(2(x,i0))") di, dj
       do i = 1, di
          write( uid, "(*(x,d24.16))" ) ( dble(A(i,j)), j = 1, dj )
          write( uid, "(*(x,d24.16))" ) (aimag(A(i,j)), j = 1, dj )
       enddo
    elseif( form .is. "unformatted" )then
       if( sum(abs(A(1:di,1:dj))) == 0.d0 )then
          write(uid) "0"
          write(uid) di, dj
       else
          write(uid) "C"
          write(uid) di, dj
          write(uid) ( ( A(i,j), i = 1, di ), j = 1, dj )
       endif
    endif
    close(uid)
  end subroutine SaveComplexMatrix_fixedsize

  subroutine LoadComplexMatrix( FileName, A, form, maxncols )
    character(len=*)               , intent(in)   :: FileName, form
    complex(kind(1d0)), allocatable, intent(inout):: A(:,:)
    !.. The optional argument maxncols is ignored in formatted form
    integer           , optional   , intent(in)   :: maxncols
    integer   :: uid, i, j, ni, nj, iostat
    character :: ch
    logical   :: exist
    complex(kind(1d0)) :: zw
    real   (kind(1d0)) :: wRe, wIm
    real   (kind(1d0)), allocatable :: ArowRe(:), ArowIm(:)

    
    INQUIRE( file = FileName, exist = exist )
    if (.not.exist ) call Assert( 'File: '//FileName//' does not exist.' )
    call OpenFile( FileName, uid, "read", form )
    if( form .is. "formatted" )then
       read(uid,"(a)")ch
       if(ch.is."C")then
          read(uid,*) ni, nj
          if(allocated(A))deallocate(A)
          allocate(A(ni,nj))
          A=(0.d0,0.d0)
          allocate(ArowRe(nj))
          allocate(ArowIm(nj))
          do i = 1, ni
             read( uid, "(*(x,d24.16))" ) ( ArowRe(j), j = 1, nj )
             read( uid, "(*(x,d24.16))" ) ( ArowIm(j), j = 1, nj )
             A(i,:) = (1.d0,0.d0) * ArowRe + (0.d0,1.d0) * ArowIm
          enddo
          deallocate(ArowRe,ArowIm)
       elseif(ch.is."Z")then
          read(uid,*) ni, nj
          if(allocated(A))deallocate(A)
          allocate(A(ni,nj))
          A=(0.d0,0.d0)
          do 
             read( uid, *, iostat=iostat ) i,j,wRe,wIm
             if(iostat/=0)exit
             A(i,j) = (1.d0,0.d0) * wRe + (0.d0,1.d0) * wIm
          enddo
       else
          call ErrorMessage("Kind mismatch in "//Filename)
          stop
       endif
    elseif( form .is. "unformatted" )then

       read(uid)ch
       if(ch.is."0")then
          read(uid) ni, nj
          if(present(maxncols)) nj = max( min( maxncols, nj ), 0 )
          if(allocated(A))deallocate(A)
          allocate(A(ni,nj))
          A=(0.d0,0.d0)
       elseif(ch.is."C")then
          read(uid) ni, nj
          if(present(maxncols)) nj = max( min( maxncols, nj ), 0 )
          if(allocated(A))deallocate(A)
          allocate(A(ni,nj))
          A=(0.d0,0.d0)
          read(uid) ( ( A(i,j), i = 1, ni ), j = 1, nj )
       elseif(ch.is."Z")then
          read(uid) ni, nj
          if(present(maxncols)) nj = max( min( maxncols, nj ), 0 )
          if(allocated(A))deallocate(A)
          allocate(A(ni,nj))
          A=(0.d0,0.d0)
          do
             read(uid,iostat=iostat) i,j,zw
             if(iostat/=0)exit
             if(j>nj)cycle
             A(i,j)=zw
          enddo
       else
          call ErrorMessage("Kind mismatch in "//Filename)
          stop
       endif
    endif
    close(uid)
  end subroutine LoadComplexMatrix

  subroutine SaveDbleVector( FileName, A, form )
    character(len=*), intent(in) :: FileName, form
    real(kind(1d0)) , intent(in) :: A(:)
    integer :: uid, i, ni
    call OpenFile( FileName, uid, "write", form )
    ni=size(A)
    if( form .is. "formatted" )then
       write(uid,"(a)") "D"
       write(uid,"(x,i0)") ni
       do i = 1, ni
          write( uid, "(x,d24.16)" )  A(i)
       enddo
    elseif( form .is. "unformatted" )then
       write(uid) "D"
       write(uid) ni
       write(uid) ( A(i), i = 1, ni )
    endif
    close(uid)
  end subroutine SaveDbleVector

  subroutine LoadDbleVector( FileName, A, form )
    character(len=*)            , intent(in)   :: FileName, form
    real(kind(1d0)), allocatable, intent(inout):: A(:)
    integer :: uid, i, ni
    character :: ch
    call OpenFile( FileName, uid, "read", form )
    if( form .is. "formatted" )then
       read(uid,"(a)")ch
       if(ch.isnt."D")then
          call ErrorMessage("Kind mismatch in "//Filename)
          stop
       endif
       read(uid,*) ni
       if(allocated(A))deallocate(A)
       allocate(A(ni))
       A=0.d0
       do i = 1, ni
          read( uid, "(x,d24.16)" ) A(i)
       enddo
    elseif( form .is. "unformatted" )then
       read(uid)ch
       if(ch.isnt."D")then
          call ErrorMessage("Kind mismatch in "//Filename)
          stop
       endif
       read(uid) ni
       if(allocated(A))deallocate(A)
       allocate(A(ni))
       A=0.d0
       read(uid) ( A(i), i = 1, ni )
    endif
    close(uid)
  end subroutine LoadDbleVector

  subroutine SaveComplexVector( FileName, A, form )
    character(len=*)  , intent(in) :: FileName, form
    complex(kind(1d0)), intent(in) :: A(:)
    integer :: uid, i, ni
    call OpenFile( FileName, uid, "write", form )
    ni=size(A)
    if( form .is. "formatted" )then
       write(uid,"(a)") "C"
       write(uid,"(x,i0)") ni
       do i = 1, ni
          write( uid, "(2(x,d24.16))" )  dble(A(i)), aimag(A(i))
       enddo
    elseif( form .is. "unformatted" )then
       write(uid) "C"
       write(uid) ni
       write(uid) ( A(i), i = 1, ni )
    endif
    close(uid)
  end subroutine SaveComplexVector

  subroutine LoadComplexVector( FileName, A, form )
    character(len=*)               , intent(in)   :: FileName, form
    complex(kind(1d0)), allocatable, intent(inout):: A(:)
    integer :: uid, i, ni
    character :: ch
    real(kind(1d0)) :: ArowRe, ArowIm
    call OpenFile( FileName, uid, "read", form )
    if( form .is. "formatted" )then
       read(uid,"(a)")ch
       if(ch.isnt."C")then
          call ErrorMessage("Kind mismatch in "//Filename)
          stop
       endif
       read(uid,*) ni
       if(allocated(A))deallocate(A)
       allocate(A(ni))
       A=(0.d0,0.d0)
       do i = 1, ni
          read( uid, "(2(x,d24.16))" ) ArowRe, ArowIm 
          A(i) = (1.d0,0.d0) * ArowRe + (0.d0,1.d0) * ArowIm
       enddo
    elseif( form .is. "unformatted" )then
       read(uid)ch
       if(ch.isnt."C")then
          call ErrorMessage("Kind mismatch in "//Filename)
          stop
       endif
       read(uid) ni
       if(allocated(A))deallocate(A)
       allocate(A(ni))
       A=(0.d0,0.d0)
       read(uid) ( A(i), i = 1, ni )
    endif
    close(uid)
  end subroutine LoadComplexVector

  subroutine ReadMatrix( FileName, Array, OnlyDimension )
    !
    character(len=*)            , intent(in)  :: FileName
    real(kind(1d0)), allocatable, intent(out) :: Array(:,:)
    logical        , optional   , intent(in)  :: OnlyDimension
    !
    integer :: uid, iostat, NumRows, NumColumns, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    logical :: opened
    !
    NumRows = 0
    NumColumns = 0
    !
    call OpenFile( FileName, uid, "read", "formatted" )
    INQUIRE( unit = uid, &
         iomsg = iomsg, &
         iostat = iostat, &
         opened = opened )
    if ( iostat /= 0 ) call Assert(iomsg)
    !
    if ( opened ) then
       !
       read(uid,*,iostat=iostat) NumRows, NumColumns
       if(iostat/=0)then
          call Assert( "Invalid dimensions on the first line in "//FileName )
       endif
       !
       allocate( Array(NumRows,NumColumns) )
       Array = 0.d0
       !
       if ( present(OnlyDimension) .and. OnlyDimension ) return
       !
       !.. Read by rows.
       do i = 1, NumRows
          read(uid,*,iostat=iostat) ( Array(i,j), j=1,NumColumns )
          if(iostat/=0)then
             call ErrorMessage( "Read error on line "//AlphabeticNumber(i+2)//" in file "//FileName )
          endif
       end do
       !
       close(uid)
       !
    else
       !
       call Assert( "The file: "//FileName//" does not exist." )
       !
    end if
    !
  end subroutine ReadMatrix



  subroutine WriteMatrix( FileName, Array )
    !
    character(len=*), intent(in) :: FileName
    real(kind(1d0)) , intent(in) :: Array(:,:)
    !
    integer :: uid, iostat, NumRows, NumColumns, i, j
    character(len=IOMSG_LENGTH) :: iomsg
    !
    NumRows    = size(Array,1)
    NumColumns = size(Array,2)
    !
    call OpenFile( FileName, uid, "write", "formatted" )
    !
    write(uid,IOSTAT=iostat,fmt=*) NumRows, NumColumns
    if ( iostat /= 0 )call Assert( "Error writing  block dimensions in "//FileName//"." )
    !
    !.. Write by rows.
    do i = 1, NumRows
       write(uid,IOSTAT=iostat,fmt=*) ( Array(i,j), j=1,NumColumns )
       if(iostat/=0) call Assert( "Write error on line "//AlphabeticNumber(j+2)//" in file "//FileName )
    end do
    !
    close(uid)
    !
  end subroutine WriteMatrix


  subroutine CheckFileMustBePresent( FileName )
    character(len=*), intent(in) :: FileName
    logical :: exist
    INQUIRE( file = FileName, exist = exist )
    if (.not.exist ) call Assert( &
         'File: '//FileName//' is not present and it must be.' )
  end subroutine CheckFileMustBePresent


end module ModuleIO
