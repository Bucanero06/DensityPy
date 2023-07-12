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
!> \file
!!
!! Defines some subroutines appropiated to work with strings.


module ModuleString

  private

  character(len=*), parameter :: UPPERCASE_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=*), parameter :: LOWERCASE_ALPHABET = "abcdefghijklmnopqrstuvwxyz" 
  character(len=*), parameter :: DIGITS_STRING = "0123456789"
 
  character, parameter, public :: ASCII_BACKSPACE = ACHAR(  8)
  character, parameter, public :: ASCII_TAB       = ACHAR(  9)
  character, parameter, public :: ASCII_NEW_LINE  = ACHAR( 10)
  character, parameter, public :: ASCII_VTAB      = ACHAR( 11)

  interface operator (.is.)
     module procedure CaseInsensitiveStrnCmp
  end interface

  interface operator (.isnt.)
     module procedure CaseInsensitiveStrnCmpNot
  end interface

  
  public :: operator(.is.)
  public :: operator(.isnt.)
  public :: SetStringToUppercase
  public :: FormatDirectoryName
  public :: AlphabeticNumber
  public :: Capitalize
  public :: ConvertToStrn
  public :: AddSlash
  public :: FindRepeatingNumberOfPatternInString
  public :: StringHasNumbers
  public :: nTokens
  public :: GetToken
  public :: FormatAsDir

contains

  function Capitalize( strn ) result (cstrn)
    character(len=*), intent(in)  :: strn
    character(len=:), allocatable :: cstrn
    cstrn=strn
    call SetStringToUppercase( cstrn )
  end function Capitalize

  !> Counts the number of tokens (a token is defined here as a sequence of non-space characters)
  integer function nTokens( strn, separator_list_ ) result( n )
    implicit none
    character(len=*)          , intent(in) :: strn
    character(len=*), optional, intent(in) :: separator_list_
    !
    character       , parameter   :: SEPARATOR_LIST_DEFAULT = " "
    character(len=:), allocatable :: separator_list
    !
    integer :: i,j
    !
    n=0
    if(len_trim( strn ) == 0)return
    !
    if(present(separator_list_))then
       allocate(separator_list,source=separator_list_)
    else
       allocate(separator_list,source=SEPARATOR_LIST_DEFAULT)
    endif
    !
    i=1
    do
       !j=verify(strn(i:)," ")
       j=verify(strn(i:),separator_list)
       if(j<=0)exit
       n=n+1
       j=i-1+j
       !i=index(strn(j:)," ")
       i=scan(strn(j:),separator_list)
       if(i<=0)exit
       i=j-1+i
       !n=n+1
       if(i>len(strn))exit
    enddo
    !
    if(allocated(separator_list))deallocate(separator_list)
    !
  end function nTokens


  subroutine GetToken( strn, iToken, token, separator_list_ )
    implicit none
    character(len=*),              intent(in) :: strn
    integer         ,              intent(in) :: iToken
    character(len=:), allocatable, intent(out):: token
    character(len=*), optional   , intent(in) :: separator_list_
    !
    character       , parameter   :: SEPARATOR_LIST_DEFAULT = " "
    character(len=:), allocatable :: separator_list
    !
    integer :: i,j,n
    !
    if(present(separator_list_))then
       allocate(separator_list,source=separator_list_)
    else
       allocate(separator_list,source=SEPARATOR_LIST_DEFAULT)
    endif
    !
    if(iToken<1)return
    if(iToken>nTokens(strn,separator_list))return
    if(allocated(token))deallocate(token)
    i=1
    n=0
    do 
       !j=verify(strn(i:)," ")
       j=verify(strn(i:),separator_list)
       if(j<=0)exit
       n=n+1
       j=i-1+j
       !i=index(strn(j:)," ")
       i=scan(strn(j:),separator_list)
       if(i<=0)then
         i=len_trim(strn)+1
       else
         i=j-1+i
       endif
       if(n==iToken)then
          allocate(token,source=strn(j:i-1))
          exit
       endif
    enddo
  end subroutine GetToken


  !> Converts a double precision number to a string.
  function ConvertToStrn(InputObject) result(strn)
    implicit none
    !> Input number.
    class(*), intent(in) :: InputObject
    character(len=:), allocatable :: tmpStrn
    character(len=:), allocatable :: Strn
    integer         :: iBuffer
    DoublePrecision :: dBuffer
    allocate(tmpStrn,source="000000000000000000000000")
    select type( Ptr => InputObject)
    type is ( Integer )
       iBuffer=Ptr
       write(tmpStrn,"(i0)") iBuffer
    type is( Real(kind(1d0)) )
       dBuffer=Ptr
       write(tmpStrn,"(d24.16)") dBuffer
    class default
       tmpStrn=" "
    end select
    allocate(Strn,source=trim(adjustl(tmpStrn)))
  end function ConvertToStrn

  !> Given a positive integer number n and a maximum number Nmax, 
  !> returns a string with as many non-blank characters as digits 
  !> in the maximum number: the last are the digits of n, the first
  !> are a repetition of a filling character. E.g.,
  !>  n=42 Nmax=10000 FillingCharacter="*" => Returns "***42".
  function AlphabeticNumber(n,Nmax,FillingCharacter) result(strn)
    implicit none
    !> Original positive integer number.
    integer  , intent(in) :: n
    !> Number of non-blanck characters in the new string containing the original number and the filling characters.
    integer  , optional, intent(in) :: Nmax
    !> The filling character.
    character, optional, intent(in) :: FillingCharacter
    character(len=:), allocatable :: Strn
    character, parameter :: DEFAULT_FILLING_CHARACTER="0"
    character :: FillingCharacter_
    !
    integer :: i,NDigitsN, NDigitsNmax
    allocate(Strn,source="000000000000000000000000")
    write(Strn,*)n
    Strn=trim(adjustl(Strn))
    NDigitsN=NDigitsInteger(n)
    if(present(Nmax))then
       NDigitsNmax=NDigitsInteger(Nmax)
    else
       NDigitsNmax=NDigitsN
    endif
    FillingCharacter_=DEFAULT_FILLING_CHARACTER
    if(present(FillingCharacter))FillingCharacter_=FillingCharacter
    do i=NDigitsN+1,NDigitsNmax
       Strn=FillingCharacter_//trim(Strn)
    enddo
  end function AlphabeticNumber


  !> Returns the number of digits of a positive integer number.
  !>
  integer function NDigitsInteger( n ) result( NDigit)
    implicit none
    integer, intent(in) :: n
    NDigit=int(log(dble(abs(n))+0.5d0)/log(10.d0))
    return
  end function NDigitsInteger


  !> Complete Dir with a "/" character at the end, if absent
  function AddSlash( Dir ) result( DirSlash )
    !> Directory name.
    character(len=*), intent(in) :: Dir
    character(len=:), allocatable :: DirSlash
    if(Dir(len_trim(Dir):len_trim(Dir))/="/")then
       allocate(DirSlash,source=trim(adjustl(Dir))//"/")
    else
       allocate(DirSlash,source=trim(adjustl(Dir)))
    endif
  end function AddSlash


  !> Complete Dir with a "/" character at the end, if absent
  subroutine  FormatDirectoryName( DirectoryPath )
    implicit none
    !> Directory path.
    character(len=*), intent(inout) :: DirectoryPath
    integer :: i
    DirectoryPath=adjustl(DirectoryPath)
    i=len_trim(DirectoryPath)
    if(DirectoryPath(i:i)/="/")DirectoryPath(i+1:i+1)="/"
  end subroutine FormatDirectoryName


  !> Complete Dir with a "/" character at the end, if absent
  function  FormatAsDir( inpDir ) result( outDir ) 
    implicit none
    character(len=*), intent(in)  :: inpDir
    character(len=:), allocatable :: outDir
    integer :: i
    i = len_trim( inpDir )
    if( i == 0 )then
       allocate( outDir, source = trim(adjustl( inpDir )) )
       return
    endif
    if( inpDir( i : i ) == "/" )then
       allocate( outDir, source = trim(adjustl( inpDir )) )
    else
       allocate( outDir, source = trim(adjustl( inpDir ))//"/" )
    endif
  end function FormatAsDir


  logical function CaseInsensitiveStrnCmp( StringA, StringB ) result( Equivalent )
    !
    character(len=*), intent(in) :: StringA, StringB
    !
    integer   :: CharacterPosition
    integer   :: LengthOfStringA
    integer   :: LengthOfStringB
    character :: charA
    character :: charB
    !
    Equivalent = .TRUE.
    !
    LengthOfStringA = len_trim( StringA )
    LengthOfStringB = len_trim( StringB )
    Equivalent = Equivalent .and. ( LengthOfStringA == LengthOfStringB )
    !
    if(.not.Equivalent) return
    !
    if( LengthOfStringA == 0 )return
    !
    do CharacterPosition = 1, LengthOfStringA
       !
       charA = UppercaseCharacter( StringA( CharacterPosition : CharacterPosition ) )
       charB = UppercaseCharacter( StringB( CharacterPosition : CharacterPosition ) )
       !
       Equivalent = Equivalent .and. ( charA == charB )
       if( .not.Equivalent ) exit
       !
    enddo
    !
  end function CaseInsensitiveStrnCmp


  logical function CaseInsensitiveStrnCmpNot( StringA, StringB ) result( NotEquivalent )
    character(len=*), intent(in) :: StringA, StringB
    NotEquivalent = .not. ( StringA .is. StringB )
  end function CaseInsensitiveStrnCmpNot


  character function UppercaseCharacter( Char ) result( UCChar )
    implicit none
    character(len=*), intent(in) :: Char
    integer :: LetterPosition
    LetterPosition = index( LOWERCASE_ALPHABET, Char )
    if( LetterPosition > 0 )then
       UCChar = UPPERCASE_ALPHABET( LetterPosition : LetterPosition )
    else
       UCChar = Char
    endif
  end function UppercaseCharacter
  
  !> Converts all the letters in a string to upper case.
  subroutine SetStringToUppercase( String ) 
    implicit none
    !> In the input is the string to be transformed, in the ourput the upper case string.
    character(len=*), intent(inout) :: String
    integer :: ichar
    character :: char
    do ichar=1,len_trim( String )
       char=String(ichar:ichar)
       String(ichar:ichar)=UppercaseCharacter(char)
    enddo
  end subroutine SetStringToUppercase
  

  !> Find how many times a substring is repeated in a string
  subroutine FindRepeatingNumberOfPatternInString( String, Pattern, N )
    implicit none
    !> Full string.
    character(len=*), intent(in)  :: String
    !> Sub-string to be found.
    character(len=*), intent(in)  :: Pattern
    !> Times the sub.string is found in the string
    integer,          intent(out) :: N
    !
    integer :: ichar, ichar2
    character(len=:), allocatable :: CopyString, CopyPattern
    !
    N = 0
    !
    allocate( CopyString, source = String )
    allocate( CopyPattern, source = Pattern )
    call SetStringToUppercase(CopyString)
    call SetStringToUppercase(CopyPattern)
    !
    ichar2 = 0
    do
       ichar = index( CopyString(ichar2+1:), CopyPattern )
       if ( ichar < 1 ) exit
       N = N + 1
       ichar2 = ichar2 + ichar
    end do
  end subroutine FindRepeatingNumberOfPatternInString
    


  logical function StringHasNumbers( Strn ) result(HasNum)
    character(len=*), intent(in) :: Strn
    integer :: i, ichar
    HasNum = .false.
    do i = 1, LEN(DIGITS_STRING)
       ichar = INDEX(Strn,DIGITS_STRING(i:i))
       if ( ichar > 0 ) then
          HasNum = .true.
          return
       end if
    end do
  end function StringHasNumbers



end module ModuleString
