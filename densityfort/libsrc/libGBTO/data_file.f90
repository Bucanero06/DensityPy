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
module data_file
use precisn
use const, only: line_len, stdout, no_header, data_file_obj_id, mib
use utils, only: xermsg

   !> \class <data_header_obj>
   !> This structure contains the data defining the header for one data set stored on the data file.
   !> \warning Only the master process actually performs the read/write operations. In case of reading the master reads the data and then broadcasts to every process. However, in case of writing every other
   !> routine that is accessing the same file as that which was written into by this object must take this into account: the current position in the file will be different for other process than the
   !> master (the other ones will be typically "behind").
   type :: data_header_obj
      !> Description (name) for the stored data.
      character(line_len) :: name = no_header
      !> Positions, within the file storing the data, of the first and the last record for this data set. The positions are given in bytes.
      integer :: first_record = -1, last_record = -1
      !> Position, within the file storing the data, of the header for the header for the next data set stored on the file.
      integer :: next_data = -1
      !> Position, within the file storing the data, of the first line of the header information.
      integer :: header_pos = -1
   contains
      !> Reads-in the header data given the starting position of the header in the file and the unit number.
      procedure :: read => read_header_obj
      !> Writes the header data into the unit number at a position given by the member variable header_pos.
      procedure :: write => write_header_obj
      !> Prints the header data into stdout.
      procedure :: print => print_header_obj
   end type data_header_obj

   private read_header_obj, write_header_obj, print_header_obj

   !> \class <data_file_obj>
   !> This object provides access to the data stored on the data file. The first record on the file is the identifier given by the constant data_file_obj_id. The next record is the number of data sets 
   !> (headers) stored on the file. After that the header of the first (if any) data set follows.
   !> \warning All read/write operations on the associated file are performed only by the master process. What the master has read-in is broadcast to the other processes. Any other routines added to this
   !> object in the future must work in the same way.
   type data_file_obj
      !> Headers and locations of the data sets stored on the file.
      type(data_header_obj), allocatable, private :: data_header(:)
      !> Number of headers on the file, i.e. the size of the array data_header.
      integer, private :: no_headers = 0
      !> Position, within the file, of the record containing the value no_headers.
      integer, private :: pos_no_headers = 0
      !> Path to the file used to store the data. Set following the call to the method open.
      character(line_len), private :: path = ''
      !> Unit number connected to the file. Set by open.
      integer, private :: unit = 0
      !> Set to .true. once the file has been opened using the method open.
      logical, private :: opened = .false.
      !> This is set to true once all records have been completed, i.e. the values of first_record, last_record have been written to the file. It is also set by load_headers.
      logical, private :: complete_records = .false.
   contains
      !> Opens the file for stream access given its path. If the file does not exist then a new one is created, otherwise the file is opened in the append mode and the existing headers for the data sets 
      !> are read-in and printed out on screen.
      !> \todo The unit number associated to the file should be returned by this function - this will allow me to get rid of get_unit_no below.
      procedure :: open => open_data_file
      !> Returns the unit number associated with the file.
      procedure :: get_unit_no => return_unit_no
      !> Returns the path associated with the file.
      procedure :: get_file_name => return_path
      !> This routine goes through the opened file and scans it for header information which is loaded into the data_header structure.
      procedure, private :: load_headers => load_data_headers
      !> Writes a new header on the file and returns the value first_record which belongs to the new header. This value can be used by other routines to write the actual data to the file starting
      !> on the position first_record. Each record must be terminated by a call to close_record with the values corresponding to the starting and final positions of the actual data record.
      procedure :: start_record => start_header_record
      !> Completes the given record given its header and the positions of the start and end of the data record.
      procedure :: close_record => close_header_record
      !> Removes the last header from the file: this can be useful when working with files that contain a broken last record.
      procedure :: remove_last_record
      !> Prints on the screen the header data.
      procedure :: print_headers => print_header_records
      !> Attempts to find the specified header and returns the positions of start and end of the corresponding data record.
      procedure :: find_header => find_data_start
      !> Returns the header containing the specified strings. Returns 0 if a unique matching header has been found and returns a non-zero number (error number) if not.
      procedure :: get_header_containing => get_data_header_containing_strings
      !> Returns the header with the given sequence number.
      procedure :: get_header => get_data_header
      !> Returns the number of headers on the file.
      procedure :: get_no_headers => get_headers
      !> Closes the data file.
      procedure :: close => close_data_file
   end type data_file_obj

   private open_data_file, load_data_headers, start_header_record, close_header_record, close_data_file, print_header_records, &
            find_data_start, return_unit_no, get_data_header, get_headers, return_path
   private get_data_header_containing_strings

contains

   function read_header_obj(this,lunit,before_header)
      use mpi_mod
      implicit none
      class(data_header_obj) :: this
      integer, intent(in) :: lunit, before_header
      integer :: read_header_obj
 
         read_header_obj = 0

         if (before_header .le. 0) read_header_obj = 1 !call xermsg('data_file','read_header_obj','The given header position is .le. 0.',1,1)

         !only master reads followed by broadcast of the data to everyone
         this%header_pos = before_header 
         if (myrank .eq. master) then
            read(lunit,pos=before_header,err=10) this%name
            read(lunit,err=10) this%first_record, this%last_record
            read(lunit,err=10) this%next_data
         endif

         call mpi_mod_bcast(this%name,master)

         call mpi_mod_bcast(this%first_record,master)

         call mpi_mod_bcast(this%last_record,master)

         call mpi_mod_bcast(this%next_data,master)

         return
 
 10      read_header_obj = 2

   end function read_header_obj

   subroutine write_header_obj(this,lunit)
      use mpi_mod
      implicit none
      class(data_header_obj) :: this
      integer, intent(in) :: lunit

         if (this%header_pos .le. 0) call xermsg('data_file','read_header_obj','The header position is set to .le. 0.',1,1)

         !Only master writes: it is very important to keep on mind that this influences the current position in the file that would be obtained using inquire by each process.
         !Therefore every other routine interacting with the file must take this into account.
         if (myrank .eq. master) then
            write(lunit,pos=this%header_pos) this%name
            write(lunit) this%first_record, this%last_record
            write(lunit) this%next_data
         endif

   end subroutine write_header_obj

   subroutine print_header_obj(this)
      implicit none
      class(data_header_obj) :: this

         write(stdout,'(10x,"Contents of Data Header:")')

         write(stdout,'("Header position: ",i0)') this%header_pos
         write(stdout,'("Header: ",a)') this%name
         write(stdout,'("Data set record start, end: ",i0," ",i0)') this%first_record, this%last_record
         write(stdout,'("Data set size [Mib]: ",f25.15)') (this%last_record-this%first_record)/real(mib,kind=cfp)
         write(stdout,'("Position of the next header: ",i0)') this%next_data

   end subroutine print_header_obj

   subroutine open_data_file(this,path)
      use mpi_mod
      implicit none
      class(data_file_obj) :: this
      character(len=*), intent(in) :: path

      integer :: err, i
      logical :: ex, op
      character(len=len('UNKNOWN')) :: rw, frm
      character(line_len) :: identifier
      integer(kind=selected_int_kind(9)) :: flt_precision, int_bit_size

         if (this%opened) call xermsg('data_file','open_data_file','The object is already associated with a file.',1,1)

         ex = .false.
         if (myrank .eq. master) then
            inquire(file=path,iostat=err,readwrite=rw,exist=ex,opened=op,unformatted=frm)
            if (err .ne. 0) then
               call xermsg('data_file','open_data_file','Error executing the inquire statement.',2,1)
            endif

            if (ex) then
               if (op) then
                  call xermsg('data_file','open_data_file','The input file is already opened.',3,1)
               end if
               if (frm == 'NO') then
                  call xermsg('data_file','open_data_file','The input file cannot be opened for unformatted access.',4,1)
               end if
               if (rw == 'NO') then
                  call xermsg('data_file','open_data_file','The input file cannot be used for read/write.',5,1)
               end if
            end if
            
            open(file=path,newunit=this%unit,access='stream',form='unformatted',position='append',iostat=err)
            if (err .ne. 0) then
               call xermsg('data_file','open_data_file','Error opening the data file.',6,1)
            endif
         endif

         !Broadcast the values that we need to all other processes; similarly below
         call mpi_mod_bcast(ex,master)

         call mpi_mod_bcast(this%unit,master)

         !if the file exists then read-in the first line and compare it against the expected identifier to see if the file can be read by this object.
         if (ex) then

            if (myrank .eq. master) then
               read(this%unit,pos=1,iostat=err) identifier
               read(this%unit) flt_precision, int_bit_size

               !Read the precision parameters in the form of 32 bit integers and compare them with the current precision to see if they are compatible
               if (flt_precision .ne. precision(cfp_dummy) .or. int_bit_size .ne. bit_size(cfp)) then
                  call xermsg ('data_file', 'open_data_file', &
                               'Attempt to read a data file that was created with a code compiled &
                               &for different precision parameters.', 7, 1)
               endif
   
               !error reading the file
               if (err .ne. 0) then
                  close(this%unit)
                  call xermsg('data_file','open_data_file','Error reading the existing data file.',7,1)
               endif
   
               if (identifier .ne. data_file_obj_id) then
                  close(this%unit)
                  call xermsg ('data_file', 'open_data_file', &
                               'The file specified does not have the correct identifier. &
                               &The input file is corrupt or incompatible.', 8, 1)
               endif
   
               inquire(this%unit,pos=i)
               this%pos_no_headers = i
            endif

            call mpi_mod_bcast(this%pos_no_headers,master)

         else !if a new file has been opened then write into it the file identifier as the first record and the number of data sets stored on the file
        
            this%no_headers = 0
            if (myrank .eq. master) then
               write(this%unit,pos=1) data_file_obj_id

               !Write the precision parameters in the form of 32 bit integers
               flt_precision = precision(cfp_dummy)
               int_bit_size = bit_size(cfp)
               write(this%unit) flt_precision, int_bit_size

               inquire(this%unit,pos=i)
               this%pos_no_headers = i

               write(this%unit) this%no_headers
            endif

            call mpi_mod_bcast(this%pos_no_headers,master)

         endif

         !reset the values of all member variables
         this%opened = .true.
         this%path = path
         if (allocated(this%data_header)) deallocate(this%data_header)
         this%complete_records = .false.

         !finally, read-in all existing (if any) headers from the file and set the value of complete_records
         call this%load_headers

         write(stdout,'(/,"DATA FILE OPENED:",a)') path
         write(stdout,'("File identifier: ",a)') data_file_obj_id

         if (this%no_headers > 0) then
            call this%print_headers
         endif

   end subroutine open_data_file

   function return_unit_no(this)
      implicit none
      class(data_file_obj) :: this
      integer :: return_unit_no

         if (.not.(this%opened)) then
            call xermsg('data_file','return_unit_no','The file has not been opened.',1,1)
         endif

         return_unit_no = this%unit

   end function return_unit_no

   function return_path(this)
      implicit none
      class(data_file_obj) :: this
      character(len=line_len) :: return_path

         if (.not.(this%opened)) then
            call xermsg('data_file','return_path','The file has not been opened.',1,1)
         endif

         return_path = this%path

   end function return_path

   subroutine load_data_headers(this)
      use mpi_mod
      implicit none
      class(data_file_obj) :: this
      character(line_len) :: identifier

      integer :: err, i, current_pos
      integer(kind=selected_int_kind(9)) :: flt_precision, int_bit_size

         if (.not.(this%opened)) then
            call xermsg('data_file','load_data_headers','The file has not been opened.',1,1)
         endif

         !read the first two lines identifying the file and giving the number of data sets stored on the file
         if (myrank .eq. master) then
            read(this%unit,pos=1) identifier
            read(this%unit) flt_precision, int_bit_size

            !Read the precision parameters in the form of 32 bit integers and compare them with the current precision to see if they are compatible
            if (flt_precision .ne. precision(cfp_dummy) .or. int_bit_size .ne. bit_size(cfp)) then
               print *,flt_precision,int_bit_size,precision(cfp_dummy),bit_size(cfp)
               call xermsg ('data_file', 'load_data_headers', &
                            'Attempt to read a data file that was created with a code compiled &
                            &for different precision parameters.', 2, 1)
            endif

            read(this%unit) this%no_headers
         endif

         !Broadcast the values read-in to the other processes; similarly below
         call mpi_mod_bcast(this%no_headers,master)

         !read-in all existing header information from the input file
         if (this%no_headers > 0) then

            allocate(this%data_header(this%no_headers),stat=err)
            if (err .ne. 0) call xermsg('data_file','load_data_headers','Memory allocation has failed.',err,1)

            if (myrank .eq. master) then
               inquire(this%unit,pos=current_pos)
            endif

            call mpi_mod_bcast(current_pos,master)

            !read-in all headers; note that the header i-1 points to the position of the header i.
            this%complete_records = .true.
            do i=1,this%no_headers
               if (i .eq. 1) then 
                  err = this%data_header(1)%read(this%unit,current_pos)
               else
                  err = this%data_header(i)%read(this%unit,this%data_header(i-1)%next_data)
               endif

               if (err .ne. 0) call xermsg('data_file','load_data_headers','header_data_obj returned an error code.',err,1)

               if (this%data_header(i)%first_record .le. 0 .or. this%data_header(i)%last_record .le. 0) then 
                  this%complete_records = .false.
               endif
            enddo

            if (.not. this % complete_records) then
                call xermsg ('data_file', 'load_data_headers', &
                             'Some header data stored on the input file are incomplete or corrupt.', 3, 1)
            end if

            !The last header should terminate the string of headers pointing to each other so the value next_header must not be a positive number.
            if (this%data_header(this%no_headers)%next_data > 0) then
               call xermsg ('data_file', 'load_data_headers', &
                            'The last data header is pointing to a location within the file but this is not allowed.', 4, 1)
            endif

         endif

   end subroutine load_data_headers

   function start_header_record(this,header)
      use mpi_mod
      implicit none
      class(data_file_obj) :: this
      character(len=*), intent(in) :: header
      integer :: start_header_record
 
      integer :: i, err
      type(data_header_obj), allocatable :: tmp_headers(:)
      character(line_len) :: head

         if (.not.(this%opened)) then
            call xermsg('data_file','start_header_record','The file has not been opened.',1,1)
         endif

         head = header
         do i=1,this%no_headers
            if (this%data_header(i)%name .eq. head) then
               print *,i,head
               call xermsg('data_file','start_header_record','The file already contains a header of the same name.',2,1)
            endif
         enddo

         !determine the position on the file where the new record will start
         if (this%no_headers > 0) then

            !we allow creation of a new header only if all previous ones have been completed.
            if (.not. this % complete_records) then
                call xermsg ('data_file', 'start_header_record', &
                             'Some records are not complete. A new record cannot be created.', 3, 1)
            end if

            start_header_record = this%data_header(this%no_headers)%last_record
            if (start_header_record <= 0) then
                call xermsg ('data_file', 'start_header_record', &
                             'The position of the new record is .le. 0. Error in data_header.', 4, 1)
            end if

            !update the last record: make it point to the location of the new header.
            this%data_header(this%no_headers)%next_data = start_header_record
            call this%data_header(this%no_headers)%write(this%unit)

            !copy the old data into a temporary structure
            call move_alloc(this%data_header,tmp_headers)

            !resize the data_header array and copy the old data back
            allocate(this%data_header(this%no_headers+1),stat=err)
            if (err .ne. 0) call xermsg('data_file','start_header_record','Memory allocation 1 failed.',err,1)

            this%data_header(1:this%no_headers) = tmp_headers(1:this%no_headers)

            this%no_headers = this%no_headers+1

         else !this will be the first header on the file

            allocate(this%data_header(1),stat=err)
            if (err .ne. 0) call xermsg('data_file','start_header_record','Memory allocation 2 failed.',err,1)

            this%no_headers = 1

            !the first record is located right after the line containing the total number of headers on the file.
            if (myrank .eq. master) then
               read(this%unit,pos=this%pos_no_headers) i
               inquire(this%unit,pos=start_header_record)
            endif

            !Broadcast this value to the other processes; similarly below
            call mpi_mod_bcast(start_header_record,master)

         endif
 
         !set-up the new header and write it to the file
         this%data_header(this%no_headers)%header_pos = start_header_record
         this%data_header(this%no_headers)%name = header
         this%data_header(this%no_headers)%first_record = -1
         this%data_header(this%no_headers)%last_record = -1
         this%data_header(this%no_headers)%next_data = -1

         call this%data_header(this%no_headers)%write(this%unit)

         !The first byte following the header data is where the data record starts and we determine it below in start_header_record.
         if (myrank .eq. master) then
            inquire(this%unit,pos=start_header_record)

            !finally, increment the total number of data sets stored on the file
            write(this%unit,pos=this%pos_no_headers) this%no_headers
         endif

         !All processes update the value of start_header_record
         call mpi_mod_bcast(start_header_record,master)

         this%data_header(this%no_headers)%first_record = start_header_record

         this%complete_records = .false.

   end function start_header_record

   subroutine close_header_record(this,header,first_record,last_record)
      use mpi_mod
      implicit none
      class(data_file_obj) :: this
      character(len=*), intent(in) :: header
      integer, intent(in) :: first_record, last_record

      character(line_len) :: head
      integer :: j, i, cnt

         if (last_record <= first_record) then
            call xermsg ('data_file', 'close_header_record', 'On input last_record .le. first_record which is invalid.', 1, 1)
         end if

         if (.not.(this%opened)) then
            call xermsg('data_file','close_header_record','The file has not been opened.',2,1)
         endif

         if (this%no_headers .le. 0) call xermsg('data_file','close_header_record','There are no records on the file.',3,1)

         !find the header that matches the one we want to close
         head = header
         cnt = 0
         do j=1,this%no_headers
            if (head .eq. this%data_header(j)%name) then 
               i = j
               cnt = cnt + 1
            endif
         enddo

         if (cnt > 1) then
            call xermsg ('data_file', 'close_header_record', 'The file contains more than one record with the same header.', 4, 1)
         end if
         if (cnt < 1) then
            call xermsg ('data_file', 'close_header_record', 'The file contains no header matching the one given.', 5, 1)
         end if
         
         if (first_record /= this % data_header(i) % first_record) then
            call xermsg ('data_file', 'close_header_record', 'The value of first_record is incompatible with the one saved.', 6, 1)
         end if

         !update the record data saved on the file
         this%data_header(i)%last_record = last_record
         call this%data_header(i)%write(this%unit)

         !check if all other headers have been closed
         this%complete_records = .true.
         do i=1,this%no_headers
            if (this%data_header(i)%first_record .le. 0 .or. this%data_header(i)%last_record .le. 0) then
               this%complete_records = .false.
            endif
         enddo

   end subroutine close_header_record

   subroutine remove_last_record(this)
      use mpi_mod
      implicit none
      class(data_file_obj) :: this

      character(line_len) :: identifier
      integer(kind=selected_int_kind(9)) :: flt_precision, int_bit_size

         if (.not.(this%opened)) then
            call xermsg('data_file','remove_last_record','The file has not been opened.',1,1)
         endif

         if (this%no_headers <= 1) then
            call xermsg ('data_file', 'remove_last_record', &
                         'There are no records on the file or only one is present: &
                         &this routine works with at least two records.', 2, 1)
         end if

         this%no_headers = this%no_headers-1

         !Remove the pointer to the last header from the one before last header.
         this%data_header(this%no_headers)%next_data = -1

         !update the data file's information on the number of headers
         if (myrank .eq. master) then
            read(this%unit,pos=1) identifier
            read(this%unit) flt_precision, int_bit_size

            write(this%unit) this%no_headers
         endif

         !update the record data saved on the file
         call this%data_header(this%no_headers)%write(this%unit)

   end subroutine remove_last_record

   subroutine print_header_records(this)
      implicit none
      class(data_file_obj) :: this
      integer :: i

         if (.not.(this%opened)) then
            call xermsg('data_file','print_header_records','The file has not been opened.',1,1)
         endif

         if (this%no_headers .eq. 0) write(stdout,'("The file is empty.")')

         write(stdout,'("Number of data sets stored on the file: ",i0)') this%no_headers
         do i=1,this%no_headers
            call this%data_header(i)%print
         enddo

   end subroutine print_header_records

   function find_data_start(this,header,first_record,last_record)
      implicit none
      class(data_file_obj) :: this
      character(len=*), intent(in) :: header
      integer :: find_data_start
      integer, intent(out) :: first_record,last_record

      character(line_len) :: head
      integer :: i, cnt, j

         if (.not.(this%opened)) then
            find_data_start = 1
            return
         endif

         find_data_start = 0
         first_record = -1
         last_record = -1
    
         head = header
         cnt = 0
         do j=1,this%no_headers
            if (this%data_header(j)%name .eq. head) then
               cnt = cnt + 1
               i = j
            endif
         enddo

         if (cnt > 1) then 
            find_data_start = 2 !call xermsg('data_file','find_data_start','The file contains more than one record with the same header.',2,1)
            return
         endif
         if (cnt < 1) then
            find_data_start = 3 !call xermsg('data_file','find_data_start','The file contains no header matching the one given.',3,1)
            return
         endif

         first_record = this%data_header(i)%first_record
         last_record = this%data_header(i)%last_record

   end function find_data_start

   subroutine close_data_file(this)
      use mpi_mod
      implicit none
      class(data_file_obj) :: this

         if (.not.(this%opened)) then
            call xermsg('data_file','close_data_file','The file has not been opened.',1,1)
         endif 

         if (myrank .eq. master) close(this%unit)
         this%opened = .false.
         this%path = ''
         this%no_headers = 0
         this%pos_no_headers = 0
         this%unit = 0
         this%opened = .false.
         this%complete_records = .false.
         if (allocated(this%data_header)) deallocate(this%data_header)

         write(stdout,'(/,"DATA FILE CLOSED")')

   end subroutine close_data_file

   subroutine get_data_header(this,i,data_header)
      implicit none
      class(data_file_obj) :: this
      integer, intent(in) :: i
      type(data_header_obj), intent(out) :: data_header

         if (.not.(this%opened)) then
            call xermsg('data_file','get_data_header','The file has not been opened.',1,1)
         endif

         if (i <= 0 .or. i > this % no_headers) then
            call xermsg ('data_file', 'get_data_header', 'The sequence number of the requested header is out of range.', 2, 1)
         end if

         data_header = this%data_header(i)

   end subroutine get_data_header

   function get_data_header_containing_strings(this,data_header,str1,str2)
      implicit none
      class(data_file_obj) :: this
      integer :: get_data_header_containing_strings
      type(data_header_obj), intent(out) :: data_header
      character(len=*), intent(in) :: str1
      character(len=*), optional, intent(in) :: str2

      integer :: i, matching, no_matched
      logical :: found

         get_data_header_containing_strings = 0

         if (.not.(this%opened)) then
            call xermsg('data_file','get_data_header_containing_strings','The file has not been opened.',1,1)
         endif

         !Find all headers that contain str1 (and str2)
         found = .false.
         matching = 0
         no_matched = 0
         do i=1,this%no_headers
            !match two strings
            if (present(str2)) then
               if (index(this % data_header(i) % name, adjustl(trim(str1))) > 0 .and. &
                   index(this % data_header(i) % name, adjustl(trim(str2))) > 0) then 
                  found = .true.
                  matching = i
                  no_matched = no_matched + 1
               endif
            else !match only one string
               if (index(this%data_header(i)%name,adjustl(trim(str1))) > 0) then
                  found = .true.
                  matching = i
                  no_matched = no_matched + 1
               endif
            endif
         enddo !i

         if (found) then
            if (no_matched > 1) then !error condition number 1: more than one matching header
               get_data_header_containing_strings = 1
               return
            else
               data_header = this%data_header(matching)
            endif
         else !error condition number 2: no matching header
            get_data_header_containing_strings = 2
            return
         endif

   end function get_data_header_containing_strings
  
   function get_headers(this)
      implicit none
      class(data_file_obj) :: this
      integer :: get_headers

         if (.not.(this%opened)) then
            call xermsg('data_file','get_headers','The file has not been opened.',1,1)
         endif

         get_headers = this%no_headers

   end function get_headers

end module data_file
