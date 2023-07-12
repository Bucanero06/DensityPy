! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
 module utils
 use precisn

 private

 !> Records the number of messages of each type written out to stdout. For list
 !> of types of messages see XERMSG.
 integer :: n_messages(-1:2) = 0

 public xermsg, report_statistics_error_messages, search_string, delete_file

contains

!> \warning This routine is not thread-safe
!> \verbatim
!>***PURPOSE  Process error messages for SLATEC and other libraries.
!>***LIBRARY   SLATEC (XERROR)
!>***CATEGORY  R3C
!>***TYPE      ALL (XERMSG-A)
!>***KEYWORDS  ERROR MESSAGE, XERROR
!>***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!>***DESCRIPTION
!>
!>   XERMSG processes a diagnostic message in a manner determined by the
!>   value of LEVEL and the current value of the library error control
!>   flag, KONTRL.  See subroutine XSETF for details.
!>
!>    LIBRAR   A character constant (or character variable) with the name
!>             of the library.  This will be 'SLATEC' for the SLATEC
!>             Common Math Library.  The error handling package is
!>             general enough to be used by many libraries
!>             simultaneously, so it is desirable for the routine that
!>             detects and reports an error to identify the library name
!>             as well as the routine name.
!>
!>    SUBROU   A character constant (or character variable) with the name
!>             of the routine that detected the error.  Usually it is the
!>             name of the routine that is calling XERMSG.  There are
!>             some instances where a user callable library routine calls
!>             lower level subsidiary routines where the error is
!>             detected.  In such cases it may be more informative to
!>             supply the name of the routine the user called rather than
!>             the name of the subsidiary routine that detected the
!>             error.
!>
!>    MESSG    A character constant (or character variable) with the text
!>             of the error or warning message.  In the example below,
!>             the message is a character constant that contains a
!>             generic message.
!>
!>                   CALL XERMSG ('SLATEC', 'MMPY',
!>                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!>                  *3, 1)
!>
!>             It is possible (and is sometimes desirable) to generate a
!>             specific message--e.g., one that contains actual numeric
!>             values.  Specific numeric values can be converted into
!>             character strings using formatted WRITE statements into
!>             character variables.  This is called standard Fortran
!>             internal file I/O and is exemplified in the first three
!>             lines of the following example.  You can also catenate
!>             substrings of characters to construct the error message.
!>             Here is an example showing the use of both writing to
!>             an internal file and catenating character strings.
!>
!>                   CHARACTER*5 CHARN, CHARL
!>                   WRITE (CHARN,10) N
!>                   WRITE (CHARL,10) LDA
!>                10 FORMAT(I5)
!>                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!>                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!>                  *   CHARL, 3, 1)
!>
!>             There are two subtleties worth mentioning.  One is that
!>             the // for character catenation is used to construct the
!>             error message so that no single character constant is
!>             continued to the next line.  This avoids confusion as to
!>             whether there are trailing blanks at the end of the line.
!>             The second is that by catenating the parts of the message
!>             as an actual argument rather than encoding the entire
!>             message into one large character variable, we avoid
!>             having to know how long the message will be in order to
!>             declare an adequate length for that large character
!>             variable.  XERMSG calls XERPRN to print the message using
!>             multiple lines if necessary.  If the message is very long,
!>             XERPRN will break it into pieces of 72 characters (as
!>             requested by XERMSG) for printing on multiple lines.
!>             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!>             so that the total line length could be 76 characters.
!>             Note also that XERPRN scans the error message backwards
!>             to ignore trailing blanks.  Another feature is that
!>             the substring '$$' is treated as a new line sentinel
!>             by XERPRN.  If you want to construct a multiline
!>             message without having to count out multiples of 72
!>             characters, just use '$$' as a separator.  '$$'
!>             obviously must occur within 72 characters of the
!>             start of each line to have its intended effect since
!>             XERPRN is asked to wrap around at 72 characters in
!>             addition to looking for '$$'.
!>
!>    NERR     An integer value that is chosen by the library routine's
!>             author.  It must be in the range -99 to 999 (three
!>             printable digits).  Each distinct error should have its
!>             own error number.  These error numbers should be described
!>             in the machine readable documentation for the routine.
!>             The error numbers need be unique only within each routine,
!>             so it is reasonable for each routine to start enumerating
!>             errors from 1 and proceeding to the next integer.
!>
!>    LEVEL    An integer value in the range 0 to 2 that indicates the
!>             level (severity) of the error.  Their meanings are
!>
!>            -1  An error message. An error mesage is reported but no
!>                abort is performed. This is used in mpi_mod before calling
!>                MPI_ABORT.
!>
!>             0  A warning message.  This is used if it is not clear
!>                that there really is an error, but the user's attention
!>                may be needed.
!>
!>             1  A recoverable error.  This is used even if the error is
!>                so serious that the routine cannot return any useful
!>                answer.  If the user has told the error package to
!>                return after recoverable errors, then XERMSG will
!>                return to the Library routine which can then return to
!>                the user's routine.  The user may also permit the error
!>                package to terminate the program upon encountering a
!>                recoverable error.
!>
!>             2  A fatal error.  XERMSG will not return to its caller
!>                after it receives a fatal error.  This level should
!>                hardly ever be used; it is much better to allow the
!>                user a chance to recover.  An example of one of the few
!>                cases in which it is permissible to declare a level 2
!>                error is a reverse communication Library routine that
!>                is likely to be called repeatedly until it integrates
!>                across some interval.  If there is a serious error in
!>                the input such that another step cannot be taken and
!>                the Library routine is called again without the input
!>                error having been corrected by the caller, the Library
!>                routine will probably be called forever with improper
!>                input.  In this case, it is reasonable to declare the
!>                error to be fatal.
!>
!>    Each of the arguments to XERMSG is input; none will be modified by
!>    XERMSG.  A routine may make multiple calls to XERMSG with warning
!>    level messages; however, after a call to XERMSG with a recoverable
!>    error, the routine should return to the user.  Do not try to call
!>    XERMSG with a second recoverable error after the first recoverable
!>    error because the error package saves the error number.  The user
!>    can retrieve this error number by calling another entry point in
!>    the error handling package and then clear the error number when
!>    recovering from the error.  Calling XERMSG in succession causes the
!>    old error number to be overwritten by the latest error number.
!>    This is considered harmless for error numbers associated with
!>    warning messages but must not be done for error numbers of serious
!>    errors.  After a call to XERMSG with a recoverable error, the user
!>    must be given a chance to call NUMXER or XERCLR to retrieve or
!>    clear the error number.
!>***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!>                 Error-handling Package, SAND82-0800, Sandia
!>                 Laboratories, 1982.
!> \endverbatim
!
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!
         use const, only: stdout, line_len, stderr
         implicit none
         character(len=*), intent(in) :: librar, subrou, messg
         integer, intent(in) :: nerr, level

         integer :: i, length

            if (LEVEL < -1 .or. LEVEL > 2) then
               write(stdout,'(1X,"***MESSAGE FROM ROUTINE",1X,a,1X,"IN LIBRARY",1X,a)') "XERMSG", "utils"
               write(stdout,'(1X,"***FATAL ERROR, PROG ABORTED, TRACEBACK REQUESTED")')
               write(stdout,'(1X,"*  ",a)') "Invalid error level number on input."
               stop 1
            endif

            write(stdout,'(1X,"***MESSAGE FROM ROUTINE",1X,a,1X,"IN LIBRARY",1X,a)') trim(SUBROU), trim(LIBRAR)

            if (LEVEL .eq. 0) write(stdout,'(1X,"***INFORMATIVE MESSAGE, PROG CONTINUES, TRACEBACK REQUESTED")')
            if (LEVEL .eq. 2 .or. LEVEL .eq. -1) write(stdout,'(1X,"***FATAL ERROR, PROG ABORTED, TRACEBACK REQUESTED")')
            if (LEVEL .eq. 1) write(stdout,'(1X,"***POTENTIALLY RECOVERABLE ERROR, PROG CONTINUES, TRACEBACK REQUESTED")')

            !Print the message in chunks of line_len-4 length
            length = len(trim(MESSG))
            do i=1,length,line_len-4
               write(stdout,'(1X,"*  ",a)') MESSG(i:min(i-1 + line_len-4,length))
            enddo
            write(stdout,'(1X,"*")')

            if (LEVEL .eq. 0) then
               write(stdout,'(1X,"*  WARNING NUMBER = ",i0)') NERR
            else
               write(stdout,'(1X,"*  ERROR NUMBER = ",i0)') NERR
            endif

            write(stdout,'(1X,"***END OF MESSAGE")')

            n_messages(LEVEL) = n_messages(LEVEL) + 1

            !todo in this case the program should continue according to the XERMGS specifications but I am using it to terminate the program completely.
            if (LEVEL .eq. 1) then
               write(stderr,'(1X,"***DUE TO AN ERROR THE PROGRAM HAS BEEN TERMINATED***")')
               STOP 1
            endif

            if (LEVEL .eq. 2) then
               write(stderr,'(1X,"***DUE TO AN ERROR THE PROGRAM HAS BEEN TERMINATED***")')
               STOP 1
            endif

            if (LEVEL .eq. -1) then
               write(stderr,'(1X,"***DUE TO AN ERROR THE PROGRAM HAS BEEN TERMINATED***")')
            endif

            ! Now, flush the I/O buffers so that we can safely call MPI_Abort right after XERMSG ends,
            ! without worries that the message just printed disappears.

            flush(stdout)
            flush(stderr)

      END SUBROUTINE XERMSG

   !> Writes to stderr unit the number of messages of each type that were written out to stdout.
   subroutine report_statistics_error_messages(myrank)
      use const, only: stderr
      implicit none
      integer, intent(in) :: myrank

         write(stderr,'(1X,"***ERROR AND WARNING MESSAGE COUNTS FOR RANK: ",i5)') myrank
         write(stderr,'(1X,"NUMBER OF WARNING MESSAGES: ",i0)') n_messages(-1)+n_messages(0)
         write(stderr,'(1X,"NUMBER OF ERROR MESSAGES: ",i0)') n_messages(1)+n_messages(2)
         write(stderr,'(1X,"THE MESSAGES CAN BE FOUND INSPECTING THE STDOUT UNIT FOR THIS RANK.")')

   end subroutine report_statistics_error_messages

   !> This function searches an input unit for the given string and returns the logical value .true. if the string has been found.
   !> In that case the input unit is positioned one line after the sought string.
   !> If the string is not found then the returned value is .false. and the input unit is positioned at its end.
   !> \param[in] unit_no unit_no Integer value corresponding to the unit to be searched.
   !> \param[in] string A string which is to be searched for on the input unit.
   !> \param[in] rew Logical value. If set to .true. then the input unit is rewound before the search is commenced.
   !> \param[in] fmted Optional logical input variable. If present and set to .false. then the reading will be done using
   !>                  unformatted read. The default is that the input unit is open for formatted reading.
   !> \todo Add a test if the input unit is open for formatted reading.
   !> \todo Add a test for unformatted unit which should search for a string not present on the file. I encountered a case
   !>       in which err > 0 was returned rather than err < 0 (end of file reached).
   function search_string(unit_no,string,rew,fmted)
      use const, only: line_len
      logical :: search_string
      integer, intent(in) :: unit_no
      character(*), intent(in) :: string
      logical, intent(in) :: rew
      logical, optional, intent(in) :: fmted
   
      integer :: err
      character(len=line_len) :: line
   
         search_string = .false.
   
         if (rew) rewind unit_no !rewind the unit if required
   
         do
   
            if (present(fmted) .and. .not.(fmted)) then !unformatted read
               read(unit_no,iostat=err,end=10,err=20) line !read-in one line
            else !formatted read
               read(unit_no,'(a)',iostat=err,end=10,err=20) line !read-in one line
            endif
   
            if (index(line,string) > 0) then !finish if the string appears on the line
               search_string = .true.
               return
            endif
   
         enddo

   10 return !end-of-file reached, i.e. the string has not been found

   20 print *,'An error occured during read.',err,line
      return
   
   end function search_string


   !> \brief   Delete a named file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Removes a named disk file by opening and closing with status = 'delete'.
   !> Does nothing if the file does not exist (or cannot be opened).
   !>
   subroutine delete_file (filename)

      character(len=*), intent(in) :: filename

      integer :: u, ierr

      open (newunit = u, iostat = ierr, file = filename, status = 'old')
      if (ierr == 0) close (u, status = 'delete')

   end subroutine delete_file

 end module utils
