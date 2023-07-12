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
!> This module implements the p2d_array_obj which is used throughout the code to store and manipulate the calculated integrals.
module parallel_arrays
   use precisn
   use mpi_mod
   use mpi_memory_mod
   use const, only: stdout, Mib, line_len

   private

   public p2d_array_obj

   !> This is an array object which is used to store the atomic and molecular integrals in memory.
   !> The primary data structure is the 2D array this%a(1:this%d1,1:this%d2). The object may use the auxiliary integer array this%block_offset(1:this%no_blocks)
   !> as in the case of computation of the atomic integrals in the MPI regime. The use of the block_offset array is controled by the initialization routine init.
   !> Each column of this%a has a name which is contained in this%column_descriptor(1:this%d2). The intention is that the different types of integrals (e.g. overlap, kinetic energy, etc.) are kept in
   !> different columns of the array this%a.
   type :: p2d_array_obj
      !> Set to .true. following a call to init.
      logical, private :: initialized = .false.
      !> Set to .true. by init if the object has been initialized with no_blocks > 0. This implies the use of the block_offset array.
      logical, private :: indexed = .false.
      !> Number of blocks.
      integer, private :: no_blocks = -1
      !> Local portion of the global array 'a': a_loc(my_start:my_end,1:d2). In the local mode my_start = 1 on all processes.
      !> The array can be shared on a single node in order to reduce memory usage in MPI-Scatci (Ahmed's implementation).
      !> \warning Never change the dimensions of this array by hand: if you need to do it then implement a method for this object otherwise you're making the d1,d2,my_start,etc. values inconsistent.
      !> \warning The type of this array variable cannot be ``real(kind=cfp), allocatable, target" since that is not allowed by Fortran 2003 standard. 
      !> \warning Therefore if MPI 3.0 is to be used this object must be rewritten to work with array a(:,:) as an external input declared with the ``target" attribute or using module-local variables.
      !> \warning I.e. the type-bound array a(:,:) should be replaced by ``real(kind=cfp), pointer :: a(:,:)" and the pointer associated when ``init" or ``read" is called.
      real(kind=cfp), pointer :: a(:,:) => null()
      !> Used in case of shared memory allocation of this%a. It stores the MPI window that this%a is stored in. If it is -1 then we are in the non-shared mode for this%a.
      integer, private :: shared_window = -1
      !> Dimensions of the array 'a'.
      integer, private :: d1 = -1, d2 = -1
      !> Pointers to start of blocks of data within the array 'a'. block_offset has dimension no_blocks.
      integer, allocatable :: block_offset(:)
      !> This array is used to name the data stored in the columns of 'a'.
      character(len=line_len), allocatable, private :: column_descriptor(:)
   contains
      !> Initializes the object. 
      !> Arguments:
      !> \param[in] n
      !> n is integer: number of rows of this%a
      !> \param[in] d2_end
      !> d2_end is integer: the number of columns of this%a
      !> \param[in] no_blocks:
      !> no_blocks is integer: number of elements of this%block_offset
      !> \param[in] column_descriptor
      !> column_descriptor is a character array of dimension d2_end: contains descriptions of the values stored in the columns of this%a.
      procedure :: init
      !> Deallocates all arrays and restores the object to an uninitialized state. It is a blocking routine.
      procedure :: final
      !> If the array is initialized with no_blocks > 0 this routine can be used to copy the contents of an external integer array to the array this%block_offset which must have the same size.
      procedure :: set_block_offset
      !> Returns the values d1,d2,no_blocks.
      procedure :: get_array_dim
      !> Sorts the array 'a' using keys given as an array of integers of the same dimension as 'a'.
      !> It is a non-blocking routine.
      procedure :: sort_array_with_keys
      !> Prints out the non-zero contents of the array this%a (and this%block_offset, if this%no_blocks > 0). It is a non-blocking routine.
      procedure :: print
      !> Writes to the disk the arrays this%a,this%block_offset from the process specified. It is a blocking routine.
      procedure :: write
      !> The arrays this%a,this%block_offset are read-in from the disk by the master task and redistributed to other processes. 
      !> If shared-memory MPI is used then each NODE keeps only one copy of the arrays this%a,this%block_offset otherwise the array is kept by every MPI task. 
      !> It is a blocking routine.
      procedure :: read
      !> In-place reduction of the local arrays this%a on all processes using MPI_ALLREDUCE. It is a blocking routine.
      procedure :: reduce_a_local
      !> Reduction of the local arrays this%a on all processes to the array kept by the master. It is a blocking routine.
      procedure :: reduce_a_to_master
      !> Logical function: does the object contain the array this%block_offset? It is a non-blocking routine.
      procedure :: have_offsets
      !> Returns the column index in the array this%a corresponding to the name of the column specified on input.
      procedure :: find_column_matching_name
      !> Returns the character array column_descriptor.
      procedure :: get_column_descriptor
      !> Combines two columns of this%a into a third one using a given operation: '+','-','*','/'.
      procedure :: combine_columns
      !> Multiplies a given column of this%a by a real number.
      procedure :: multiply_column
   end type p2d_array_obj

contains

   function have_offsets(this)
      implicit none
      class(p2d_array_obj) :: this
      logical :: have_offsets

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'have_offsets', 'The object has not been initialized.', 1, 1)
         end if

         have_offsets = this%indexed

   end function have_offsets

   subroutine sort_array_with_keys(this,n,d2,keys,sort_method)
      use sort, only: heap_sort_int_float, sort_int_float
      use omp_lib
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: n, d2, sort_method
      integer, intent(inout) :: keys(n,d2)

      real(kind=wp) :: start_t, end_t

         start_t = omp_get_wtime()

         write(stdout,'("--------->","p2d_array_obj:sort_array_with_keys")')

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'sort_array_with_keys', 'The object has not been initialized.', 1, 1)
         end if

         if (n > this % d1 .or. n <= 0 .or. d2 <= 0 .or. d2 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'sort_array_with_keys', &
                             'On input n and/or d2 were out of allowed limits for the array this%a.', 2, 1)
         end if

         select case (sort_method)
            case (1)
               call heap_sort_int_float(n,d2,keys,this%a)
            case (2)
               call sort_int_float(n,d2,keys,this%a)
            case default
               call mpi_xermsg ('parallel_arrays', 'sort_array_with_keys', &
                                'On input sort_method was out of the allowed range [1,2].', 3, 1)
         end select

         write(stdout,'("<---------","done:p2d_array_obj:sort_array_with_keys")')

         end_t = omp_get_wtime()

         write(stdout,'("sort_array_with_keys took [s]",f25.15)') end_t-start_t

   end subroutine sort_array_with_keys

   function init(this,n,d2_end,no_blocks,column_descriptor)
      implicit none
      integer :: init
      class(p2d_array_obj) :: this
      integer, intent(in) :: n, d2_end,no_blocks
      character(len=line_len), intent(in) :: column_descriptor(:)

      integer :: err, i, j, memory
      integer, allocatable :: nelem(:)

         init = 0

         call check_mpi_running

         call mpi_mod_barrier(err)

         if (this%initialized) then
            err = this%final()
            if (err .ne. 0) call mpi_xermsg('parallel_arrays','init','Error finalizing the object.',1,1)
         endif

         if (d2_end .le. 0 .or. n .le. 0) then
            print *,d2_end,n
            call mpi_xermsg('parallel_arrays','init','d2_end .le. 0 .or. n .le. 0 but both must be .ge. 1.',2,1)
         endif

         if (size(column_descriptor) .ne. d2_end) then
            print *,size(column_descriptor),d2_end
            call mpi_xermsg ('parallel_arrays', 'init', &
                             'size of the input array column_descriptor .ne. d2_end but they must be equal.', 3, 1)
         endif

         if (no_blocks > 0) then
            write(stdout,'("--------->","p2d_array_obj:init; indexed array")')
            this%indexed = .true.
         else
            write(stdout,'("--------->","p2d_array_obj:init; unindexed array")')
            this%indexed = .false.
         endif

         !Transfer the column labels
         allocate(this%column_descriptor,source=column_descriptor,stat=err)
         if (err .ne. 0) call mpi_xermsg('parallel_arrays','init','Memory allocation 0 failed.',err,1)

         allocate(nelem(1:nprocs),stat=err)
         if (err .ne. 0) call mpi_xermsg('parallel_arrays','init','Memory allocation 1 failed.',err,1)

         !Check that all processes request the same type of array (shared or local)
         nelem(:) = 0
         call mpi_mod_allgather(no_blocks,nelem)

         do i=1,nprocs
            do j=1,i
               if (nelem(i) > 0 .and. nelem(j) <= 0) then
                    call mpi_xermsg ('parallel_arrays', 'init', &
                                     'At least one process requests shared array while another one requests local array.', 4, 1)
               end if
            enddo
         enddo

         !Check that the d2_end dimension is the same on all processes: choose arbitrarilly the root's d2_start value as the reference.
         nelem(:) = 0
         call mpi_mod_allgather(d2_end,nelem)

         if (nelem(myrank+1) /= nelem(1)) then
            call mpi_xermsg ('parallel_arrays', 'init', 'd2_start value must be the same for all processes.', 5, 1)
         end if

         !Send the size of the local array 'a' to everyone so that each process knows exactly how many elements every other process keeps.
         nelem(:) = 0
         call mpi_mod_allgather(n,nelem)

         call mpi_mod_barrier(err) !putting the barrier here ensures that we catch any errors in the dimensioning before the allocations take place.

         this%d1 = n
         this%d2 = d2_end
         this%no_blocks = no_blocks

         write(stdout,'(10X,"rank ",i0,", array dimensions: ",i10,1X,i10)') myrank, n, d2_end
         if (this%indexed) then
            write(stdout,'(10X,"rank ",i0,", number of blocks: ",i10)') myrank, no_blocks
         endif

         !count the number of bytes
         memory = n*d2_end*cfp_bytes
         write(stdout,'("Memory to be allocated for the array [Mib]: ",f12.3)') memory/(Mib*1.0_cfp)

         if (associated(this%a)) deallocate(this%a)
         if (allocated(this%block_offset)) deallocate(this%block_offset)

         if (this%indexed) then
            memory = no_blocks*bit_size(memory)/8
            write(stdout,'("Memory to be allocated for the block offset [Mib]: ",f12.3)') memory/(Mib*1.0_cfp)
            allocate(this%block_offset(this%no_blocks),stat=err)
            this % shared_window = mpi_memory_allocate_real_2dim(this % a, n, d2_end)
         else
            this % shared_window = mpi_memory_allocate_real_2dim(this % a, n, d2_end)
         endif

         this%a = 0.0_cfp
         if (this%indexed) this%block_offset = -1

         this%initialized = .true.

         call mpi_mod_barrier(err)

         write(stdout,'("<---------","done:p2d_array_obj:init")')

   end function init

   function final(this)
      implicit none
      integer :: final
      class(p2d_array_obj) :: this

      integer :: ierr

         final = 0

         call check_mpi_running

         call mpi_mod_barrier(ierr)

         write(stdout,'("--------->","p2d_array_obj:final")')

         ierr = 0
         if (allocated(this%block_offset)) deallocate(this%block_offset,stat=ierr)
         if (ierr .ne. 0) call mpi_xermsg('parallel_arrays','final','Memory deallocation 3 failed.',ierr,1)
         ierr = 0
         if (allocated(this%column_descriptor)) deallocate(this%column_descriptor,stat=ierr)
         if (ierr .ne. 0) call mpi_xermsg('parallel_arrays','final','Memory deallocation 4 failed.',ierr,1)
         
         if (this%shared_window /= -1) then
            call mpi_memory_deallocate_real_2dim(this%a,size(this%a),this%shared_window)
            this%shared_window = -1
         else
            ierr = 0
            if (associated(this%a)) deallocate(this%a,stat=ierr)
            if (ierr .ne. 0) call mpi_xermsg('parallel_arrays','final','Memory deallocation 1 failed.',ierr,1)
         endif

         nullify(this % a)

         this%initialized = .false.
         this%indexed = .false.
         this%d1 = -1
         this%d2 = -1
         this%no_blocks = -1

         call mpi_mod_barrier(ierr)

         write(stdout,'("<---------","done:p2d_array_obj:final")')

   end function final

   subroutine set_block_offset(this,block_offset)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: block_offset(:)

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'set_block_offset', 'The object has not been initialized.', 1, 1)
         end if

         if (this % no_blocks /= size(block_offset)) then
            call mpi_xermsg ('parallel_arrays', 'set_block_offset', &
                             'On input size of block_offset does not match the size for which the object was initialized.', 2, 1)
         end if

         this%block_offset(1:this%no_blocks) = block_offset(1:this%no_blocks)

   end subroutine set_block_offset

   subroutine get_array_dim(this,d1,d2,no_blocks)
      use omp_lib
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(out) :: d1,d2,no_blocks

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','get_array_dim','The object has not been initialized.',1,1)

         d1 = this%d1
         d2 = this%d2
         no_blocks = this%no_blocks

   end subroutine get_array_dim

   subroutine print(this,only_non_zero)
      implicit none
      class(p2d_array_obj) :: this
      logical, intent(in) :: only_non_zero

      integer :: i, j

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','print','The object has not been initialized.',1,1)

         write(stdout,'("--------->","p2d_array_obj:print")')

         if (this%indexed) then
            write(stdout,'("Number of block offsets: ",i10)') j, this%no_blocks
            if (only_non_zero) then
               do i=1,this%no_blocks
                  if (this%block_offset(i) .ne. 0) write(stdout,'(i10,1X,i20)') i,this%block_offset(i)
               enddo
            else
               do i=1,this%no_blocks
                  write(stdout,'(i10,1X,i20)') i,this%block_offset(i)
               enddo
            endif
         endif

         do j=1,this%d2
            write(stdout,'("Column ",i5," descriptor: ",a)') j, this%column_descriptor(j)
            if (only_non_zero) then
               do i=1,this%d1
                  if (this%a(i,j) .ne. 0.0_cfp) write(stdout,'(i0,1x,e25.15)') i,this%a(i,j)
               enddo !i
            else
               do i=1,this%d1
                  write(stdout,'(i0,1x,e25.15)') i,this%a(i,j)
               enddo !i
            endif
         enddo !j

         write(stdout,'("<---------","done:p2d_array_obj:print")')

   end subroutine print

   function find_column_matching_name(this,column_name)
      implicit none
      class(p2d_array_obj) :: this
      character(len=*), intent(in) :: column_name
      integer :: find_column_matching_name

      integer :: i
      integer :: a_size
      logical :: found

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'find_column_matching_name', 'The object has not been initialized.', 1, 1)
         end if
         a_size = size(this%a,2)
         found = .false.
         do i=1,a_size
            !write(stdout,"(a,' ',a)") trim(this%column_descriptor(i)) ,' ',trim(column_name)
            if (trim(this%column_descriptor(i)) .eq. trim(column_name)) then
               find_column_matching_name = i
               found = .true.
               exit
            endif
         enddo !i

         if (.not.(found)) then
            print *,'required column:',column_name
            print *,'available columns:',this%column_descriptor(:)
            call mpi_xermsg ('parallel_arrays', 'find_column_matching_name', &
                             'A column matching the input name has not been found in the array.', 2, 1)
         endif

   end function find_column_matching_name

   subroutine get_column_descriptor(this,column_descriptor)
      implicit none
      class(p2d_array_obj) :: this
      character(len=line_len), allocatable :: column_descriptor(:)

      integer :: err

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'get_column_descriptor', 'The object has not been initialized.', 1, 1)
         end if

         if (allocated(column_descriptor)) deallocate(column_descriptor)

         allocate(column_descriptor,source=this%column_descriptor,stat=err)
         if (err .ne. 0) call mpi_xermsg('parallel_arrays','get_column_descriptor','Memory allocation failed.',err,1)

   end subroutine get_column_descriptor

   !> \todo the 'rank' parameter should be allowed to be < 0 in which case all processes write equal portion of the integrals array thereby speeding up the write on Lustre system.
   subroutine write(this,lunit,record_start,position_after_write,rank)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: lunit, rank
      integer, intent(in) :: record_start
      integer, intent(out) :: position_after_write

      integer :: err

      real :: start_t, end_t

         call cpu_time(start_t)

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','write','The object has not been initialized.',1,1)

         write(stdout,'("--------->","p2d_array_obj:write")')

         call mpi_mod_barrier(err)

         if (rank < 0 .or. rank > nprocs-1) call mpi_xermsg('parallel_arrays','write','On input rank is out of range.',3,1)

         if (myrank .eq. rank) then
            write(lunit,pos=record_start,err=10) size(this%a,1),size(this%a,2)
            write(lunit) this%column_descriptor
            write(lunit,err=10) this%a

            write(lunit) this%no_blocks
            if (this%no_blocks > 0) then
               write(lunit) this%block_offset
            endif
            inquire(lunit,pos=position_after_write)
         endif

         !get the position_after_write on all processes
         call mpi_mod_bcast(position_after_write,int(rank,kind=mpiint))

         write(stdout,'("<---------","done:p2d_array_obj:write")')

         call cpu_time(end_t)
         
         write(stdout,'("write took [s] ",f25.15)') end_t-start_t

         return

 10      call mpi_xermsg ('parallel_arrays', 'write', &
                          'Error executing the write command while writing the array data to the disk.', 4, 1)

   end subroutine write

   subroutine read(this,file_name,lunit,record_start,position_after_read)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: lunit
      integer, intent(in) :: record_start
      character(len=*), intent(in) :: file_name
      integer, intent(out) :: position_after_read

      integer :: err, buflen, fh, disp, d1, d2, no_blocks,i,j
      integer, allocatable :: nelem(:)
      integer(kind=mpiint) :: test
      real(kind=cfp), allocatable :: buffer(:)
      character(len=line_len), allocatable :: column_descriptor(:)
      integer(mpiint)  :: local_master_array(nprocs)
      real :: start_t, end_t

        
         call cpu_time(start_t)
         write(stdout,'("--------->","p2d_array_obj:read")')

         call mpi_mod_barrier(err)

         if (this%initialized) then
            err = this%final()
         endif

         if (myrank .eq. master) then
            read(lunit,pos=record_start,err=50) this%d1,this%d2
         endif
         !Everyone needs this before hand
         call mpi_mod_bcast(this%d1,master)
         call mpi_mod_bcast(this%d2,master)
         
         !Allocate the shared memory   
         this%shared_window = mpi_memory_allocate_real_2dim(this%a,this%d1,this%d2)
         !now we get process one in the local node to do this
         !Fence off the memory for it

         call mpi_memory_synchronize(this%shared_window)

         if (myrank .eq. master) then    
            allocate(this%column_descriptor(this%d2),stat=err)
            if (err .ne. 0) then
               print *,this%d1,this%d2
               call mpi_xermsg('parallel_arrays','read','Master memory allocation failed.',err,1)
            endif

            read(lunit) this%column_descriptor
            inquire(unit=lunit,pos=position_after_read)
            read(lunit) this%a
            read(lunit) this%no_blocks

            if (this%no_blocks > 0) then
               print *,'blocks',this%no_blocks
               allocate(this%block_offset(this%no_blocks),stat=err)
               if (err .ne. 0) call mpi_xermsg('parallel_arrays','read','Memory allocation 1 failed.',err,1)
               read(lunit) this%block_offset
            endif

            inquire(unit=lunit,pos=position_after_read)
         endif

         !Update every other process with the memory
         call mpi_memory_synchronize(this%shared_window)

         !Master must ensure every process knows what the array's dimensions are.
         call mpi_mod_bcast(this%no_blocks,master)

         !Master must ensure every process knows where in the file the data record ends.
         call mpi_mod_bcast(position_after_read,master)
        
         !Allocate everthing else except for a
         if (myrank .ne. master) then
            allocate(this%column_descriptor(this%d2),stat=err)
            if (err .ne. 0) call mpi_xermsg('parallel_arrays','read','Memory allocation 2 failed.',err,1)
            if (this%no_blocks > 0) then
               allocate(this%block_offset(this%no_blocks),stat=err)
               if (err .ne. 0) call mpi_xermsg('parallel_arrays','read','Memory allocation 3 failed.',err,1)
            endif
         endif

         call mpi_mod_bcast(this%column_descriptor,master)

         if (shared_enabled) then !Every node has only one copy of the array this%a

            local_master_array = 1
            call mpi_mod_allgather(local_rank,local_master_array)
            write(stdout,"('Gethered all proc ids')")
            call mpi_memory_synchronize(this%shared_window)
            write(stdout,"('Fencing')")
            if (myrank .eq. master) then
               do i=1,nprocs
                  if (i-1 .eq. master) cycle
                  write(stdout,"(i4)") i
                  if (local_master_array(i) == local_master) then 
                     write(stdout,"('Sending array to ',i4)") (i-1)
                     do j=1,this%d2
                        call mpi_mod_send(int(i-1,mpiint),this%a(1:this%d1,j),1,this%d1)
                        if (this % no_blocks > 0) then
                            call mpi_mod_send (int(i - 1, mpiint), this % block_offset(1 : this % no_blocks), 1, this % no_blocks)
                        end if
                     enddo   
                  endif
               enddo
            endif
            
            if (local_rank == local_master .and. myrank /= master) then 
               do j=1,this%d2
                  write(stdout,"('Gathering info my array is size ',i12)") size(this%a,1)
                  call mpi_mod_recv(master,1,this%a(1:this%d1,j),this%d1)
                  if (this%no_blocks > 0) call mpi_mod_recv(master,1,this%block_offset(1:this%no_blocks),this%no_blocks)
               enddo
            endif
    
            call mpi_memory_synchronize(this%shared_window)

         else !every MPI task keeps a copy of the arrays this%a,this%block_offset

            if (this%no_blocks > 0) call mpi_mod_bcast(this%block_offset,master)

            !Broadcast the array one dimension at a time - this can be improved
            !to read the full 2d array at once but for most applications what
            !follows should be enough.
            do d2=1,this%d2
               call mpi_mod_bcast(this%a(1:this%d1,d2),master)
            enddo

         endif

         this%initialized = .true.

         call mpi_mod_barrier(err)

         write(stdout,'("<---------","done:p2d_array_obj:read")')

         call cpu_time(end_t)

         write(stdout,'("read took [s] ",f25.15)') end_t-start_t

         return
   
 50      call mpi_xermsg ('parallel_arrays', 'read', &
                          'Error executing the read command while reading the array data from the disk.', 3, 1)
   
   end subroutine read
   
   subroutine reduce_a_local(this)
      implicit none
      class(p2d_array_obj) :: this

      integer :: err, i, nelem
      real :: start_t, end_t

         call cpu_time(start_t)

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','reduce_a_local','The object has not been initialized.',1,1)

         write(stdout,'("--------->","p2d_array_obj:reduce_a_local")')

         call mpi_mod_barrier(err)

         !For the reduction to work all processes must keep the same number of elements of the array to be reduced. The index range for 'a' is guaranteed to be the same if the number of elements 
         !agree since the local array 'a' on each process starts with index 1. To be absolutely sure that the dimension is the same on all processes we arbitrarily select the master as the reference
         !process and use its value of size(this%a,1) to check that all processes's local arrays 'a' have the same dimension.
         nelem = size(this%a,1)
         call mpi_mod_bcast(nelem,master)
         
         if (nelem .ne. size(this%a,1)) then
            print *,myrank,nelem,size(this%a,1)
            call mpi_xermsg ('parallel_arrays', 'reduce_a_local', &
                             'The local dimensions of the array a to be reduced must be the same on all processes but isnt.', 3, 1)
         endif

         write(stdout,'(10X,"Reducing (summing) ",i0," rows of the local array a accross all processes...")') nelem
       
         !This can be rewritten so that the reduction is on the whole 2D array but for most applications this should be good enough. 
         do i=1,size(this%a,2)
            call mpi_reduceall_inplace_sum_cfp(this%a(1:size(this%a,1),i), nelem)
         enddo

         write(stdout,'(10X,"The local array on all processes has been replaced by its reduced version &
                            &obtained by summing contributions from all processes.")')

         write(stdout,'("--------->","done:p2d_array_obj:reduce_a_local")')

         call cpu_time(end_t)

         write(stdout,'("reduce_a_local took [s] ",f25.15)') end_t-start_t

   end subroutine reduce_a_local

   subroutine reduce_a_to_master(this)
      implicit none
      class(p2d_array_obj) :: this

      integer :: err, i, nelem
      real :: start_t, end_t

      integer, parameter :: tag = 1

         call cpu_time(start_t)

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'reduce_a_to_master', 'The object has not been initialized.', 1, 1)
         end if

         write(stdout,'("--------->","p2d_array_obj:reduce_a_to_master")')

         call mpi_mod_barrier(err)

         !For the reduction to work all processes must keep the same number of elements of the array to be reduced. The index range for 'a' is guaranteed to be the same if the number of elements 
         !agree since the local array 'a' on each process starts with index 1. To be absolutely sure that the dimension is the same on all processes we arbitrarily select the master as the reference
         !process and use its value of size(this%a,1) to check that all processes's local arrays 'a' have the same dimension.
         nelem = size(this%a,1)
         call mpi_mod_bcast(nelem,master)
         
         if (nelem .ne. size(this%a,1)) then
            print *,myrank,nelem,size(this%a,1)
            call mpi_xermsg ('parallel_arrays', 'reduce_a_to_master', &
                             'The local dimensions of the array a to be reduced must be the same on all processes but isnt.', 3, 1)
         endif

         write(stdout,'(10X,"Master is reducing (summing) ",i0," rows of the local array a using &
                            &contributions from all processes...")') nelem

         !This can be rewritten so that the reduction is on the whole 2D array but for most applications this should be good enough. 
         do i=1,size(this%a,2)
            call mpi_reduce_inplace_sum_cfp(this%a(1:size(this%a,1),i), nelem)
         enddo

         write(stdout,'(10X,"The local array on the Master process has been replaced by its reduced version obtained &
                            &by summing contributions from all processes.")')

         write(stdout,'("--------->","done:p2d_array_obj:reduce_a_to_master")')

         call cpu_time(end_t)

         write(stdout,'("reduce_a_to_master took [s] ",f25.15)') end_t-start_t

   end subroutine reduce_a_to_master

   subroutine combine_columns(this,col1,op,col2,col3)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: col1,col2,col3
      character(len=1), intent(in) :: op

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'The object has not been initialized.', 1, 1)
         end if

         if (col1 <= 0 .or. col1 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'On input col1 was out of range.', 2, 1)
         end if
         if (col2 <= 0 .or. col2 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'On input col2 was out of range.', 3, 1)
         end if
         if (col3 <= 0 .or. col3 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'On input col3 was out of range.', 4, 1)
         end if

         if (shared_enabled) then

            call mpi_memory_synchronize(this%shared_window)
            if (local_rank == local_master) then
               select case (op)
                  case ('+')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) + this%a(1:this%d1,col2)
                  case ('-')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) - this%a(1:this%d1,col2)
                  case ('*')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) * this%a(1:this%d1,col2)
                  case ('/')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) / this%a(1:this%d1,col2)
                  case default
                     call mpi_xermsg('parallel_arrays','combine_columns','On input op was not one of: +,-,*,/.',5,1)
               end select
            endif
            call mpi_memory_synchronize(this%shared_window)

          else !each MPI task has its own copy of this%a:

            select case (op)
               case ('+')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) + this%a(1:this%d1,col2)
               case ('-')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) - this%a(1:this%d1,col2)
               case ('*')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) * this%a(1:this%d1,col2)
               case ('/')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) / this%a(1:this%d1,col2)
               case default
                  call mpi_xermsg('parallel_arrays','combine_columns','On input op was not one of: +,-,*,/.',6,1)
            end select

         endif
   end subroutine combine_columns

   subroutine multiply_column(this,col,fac)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: col
      real(kind=cfp), intent(in) :: fac

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'multiply_column', 'The object has not been initialized.', 1, 1)
         end if
         if (col <= 0 .or. col > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'multiply_column', 'On input col was out of range.', 2, 1)
         end if

         if (shared_enabled) then

            call mpi_memory_synchronize(this%shared_window)

            if (local_rank == local_master) then
               this%a(1:this%d1,col) = this%a(1:this%d1,col)*fac
            endif

            call mpi_memory_synchronize(this%shared_window)
         else
            this%a(1:this%d1,col) = this%a(1:this%d1,col)*fac
         endif

   end subroutine multiply_column

end module parallel_arrays
