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

!> \brief   MPI library interfaces
!> \authors Z Masin, J Benda
!> \date    2014 - 2019
!>
!> Convenience generic interfaces for MPI calls. Specific implementations available for 4-byte and 8-byte integers
!> and 8-byte and 16-byte real data types.
!>
!> \warning This module is also used by SCATCI but always in DOUBLE PRECISION. However, SCATCI must always be linked
!> to the integral library of precision used to calculate the integrals file.
!> Therefore, in order to ensure correct functionality of SCATCI we need to ensure this module supports calls
!> to both DOUBLE and QUAD precision MPI routines within the same program run.
!> This is achieved by interfacing all routines manipulating floating point variables for both wp and ep precisions.
!>
module mpi_mod

   use iso_fortran_env, only: int32, int64, real64, real128
   use precisn,         only: shortint, longint, wp, ep, ep1, cfp, storage_unit_int, wp_bytes, ep1_bytes, cfp_bytes, &
                              print_precision_params
   use utils,           only: xermsg, report_statistics_error_messages
   use const,           only: stdout, stdout_unit_base, redirect_stdout_to_file

#ifdef usempi
   use mpi
#endif

   implicit none

   !integer types that will be used in calls to the MPI_ routines:
#ifdef usempi
   integer, parameter, public :: mpiint  = kind(MPI_COMM_WORLD)
   integer, parameter, public :: mpiaddr = MPI_ADDRESS_KIND
   integer, parameter, public :: mpicnt  = MPI_COUNT_KIND
   integer, parameter, public :: mpiofs  = MPI_OFFSET_KIND
#else
   integer, parameter, public :: mpiint  = kind(0)
   integer, parameter, public :: mpiaddr = longint
   integer, parameter, public :: mpicnt  = longint
   integer, parameter, public :: mpiofs  = longint
#endif

   !> This can be set only once at the beginning of the program run by a call to MPI_MOD_START. The value of this variable
   !> is accessed by MPI_MOD_CHECK which aborts the program if MPI is not running. This variable stays .false. if the
   !> library was compiled without MPI support.
   logical, public, protected :: mpi_running = .false.

   !> This is similar to \ref mpi_running, but will be always set to true once MPI_MOD_START is called, whether
   !> or not the library was built with MPI support. So, while \ref mpi_running can be used to determine whether
   !> the program is running under MPI, the variable \ref mpi_started indicates whether the call to MPI_MOD_START
   !> had been already done (even in serial build).
   logical, public, protected :: mpi_started = .false.

   !> Total number of processes. Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: nprocs = -1

   !> The local process number (rank). Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: myrank = -1

   !> Total number of processes on the node this task is bound to. Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: local_nprocs = 1

   !> The local process number (rank) on the node this task is bound to. Set by mpi_mod_start.
   integer(kind=mpiint), public, protected :: local_rank = 0

   !> The routine mpi_start creates the shared_communicator which groups all tasks sitting on the same node. If we're using
   !> MPI 3.0 standard the shared_communicator enables creation of shared memory areas. In this case shared_enabled is set
   !> to .true. If we're not using MPI 3.0 then shared_enabled is set to .false.
   logical, public, protected :: shared_enabled = .false.

   !> ID of the master process.
   !> \warning This must be set to 0 so don't change it.
   integer(kind=mpiint), parameter, public :: master = 0

   !> This variable is used only in max_data_count below but it MUST be 32 bit integer.
   integer(kind=shortint), private, parameter :: dummy_32bit_integer = 1
   !> This variable is used only in max_data_count below but it MUST be 64 bit integer.
   integer(kind=longint), private, parameter :: dummy_64bit_integer = 1
#ifdef usempi
   !> Type value corresponding to the default integer type. Set by mpi_mod_start.
   integer(kind=mpiint), private :: mpi_mod_int = -1

   !> Type value corresponding to the MPI integer type. Set by mpi_mod_start.
   integer(kind=mpiint), private :: mpi_mod_mpiint = -1

   !> MPI type value corresponding to the default (cfp) float. Set by mpi_mod_start.
   integer(kind=mpiint), private :: mpi_mod_cfp = -1

   !> MPI type value corresponding to the double precision float. Set by mpi_mod_start.
   integer(kind=mpiint), private :: mpi_mod_wp = -1

   !> MPI type value corresponding to the quad precision float. Set by mpi_mod_start.
   integer(kind=mpiint), private :: mpi_mod_ep = -1

   !> Name of the processor on which the current process is running.
   character(len=MPI_MAX_PROCESSOR_NAME), public, protected :: procname

   !> Intra-node communicator created by mpi_mod_start.
   integer(kind=mpiint), public, protected :: shared_communicator

   !> Integer identifying the node the MPI task belongs to.
   integer(kind=mpiint), private :: node_colour
#else
   character(len=*), parameter, public :: procname = "N/A"
   integer(kind=mpiint), public, protected :: shared_communicator = -1
#endif
#if defined(usempi)
   !> The largest allowed data type count (i.e. the number of array elements) that can be passed to the MPI routines for all available data types. 
   !> This limitation exists in Intel MPI but perhaps also in other MPI libraries? The maximum size of the message is given in bytes by the value of huge(dummy_32bit_integer).
   !> These values are used in the routines mpi_mod_bcast_array_default_integer, mpi_mod_bcast_wp_array, mpi_mod_bcast_ep_array to split large messages into smaller chunks.
   !> In principle this should be implemented for the other MPI routines which work with arrays. This will be done only when needed.
   !> todo the same principle is implemented separately in the routines naive_mpi_reduce_inplace_sum_* so these routines should be changed to use the parameters defined here.
   integer, parameter, private :: max_data_count_longint = huge(dummy_32bit_integer)/(bit_size(dummy_64bit_integer)/8)
   integer, parameter, private :: max_data_count_wp = huge(dummy_32bit_integer)/wp_bytes
   integer, parameter, private :: max_data_count_ep = huge(dummy_32bit_integer)/ep1_bytes
#endif

   interface mpi_mod_bcast
      module procedure mpi_mod_bcast_logical
      module procedure mpi_mod_bcast_character, mpi_mod_bcast_character_array
      module procedure mpi_mod_bcast_int32,     mpi_mod_bcast_int32_array
      module procedure mpi_mod_bcast_int64,     mpi_mod_bcast_int64_array
      module procedure mpi_mod_bcast_wp,        mpi_mod_bcast_wp_array
      module procedure mpi_mod_bcast_ep,        mpi_mod_bcast_ep_array
   end interface mpi_mod_bcast

   interface mpi_mod_send
      module procedure mpi_mod_send_int32_array
      module procedure mpi_mod_send_int64_array
      module procedure mpi_mod_send_wp_array
      module procedure mpi_mod_send_ep_array
   end interface mpi_mod_send

   interface mpi_mod_send_dynamic
      module procedure mpi_mod_send_wp_array_dynamic
      module procedure mpi_mod_send_ep_array_dynamic
   end interface mpi_mod_send_dynamic

   interface mpi_mod_isend
      module procedure mpi_mod_isend_int32_array
      module procedure mpi_mod_isend_int64_array
   end interface mpi_mod_isend

   interface mpi_mod_recv
      module procedure mpi_mod_recv_int32_array
      module procedure mpi_mod_recv_int64_array
      module procedure mpi_mod_recv_wp_array
      module procedure mpi_mod_recv_ep_array
   end interface mpi_mod_recv

   interface mpi_mod_recv_dynamic
      module procedure mpi_mod_recv_wp_array_dynamic
      module procedure mpi_mod_recv_ep_array_dynamic
   end interface mpi_mod_recv_dynamic

   interface mpi_mod_file_read_cfp
      module procedure mpi_mod_file_read_wp
      module procedure mpi_mod_file_read_ep
   end interface mpi_mod_file_read_cfp

   interface mpi_mod_file_write
      module procedure mpi_mod_file_write_int32, mpi_mod_file_write_array1d_int32, mpi_mod_file_write_array2d_int32
      module procedure mpi_mod_file_write_real64, mpi_mod_file_write_array1d_real64, mpi_mod_file_write_array2d_real64, &
                       mpi_mod_file_write_array3d_real64, mpi_mod_file_write_darray2d_real64
   end interface mpi_mod_file_write

   interface mpi_mod_allgather
      module procedure mpi_mod_allgather_character
      module procedure mpi_mod_allgather_int32
      module procedure mpi_mod_allgather_int64
   end interface mpi_mod_allgather

   interface mpi_mod_file_set_view
      module procedure mpi_mod_file_set_view_wp
      module procedure mpi_mod_file_set_view_ep
   end interface mpi_mod_file_set_view

   interface mpi_reduce_inplace_sum_cfp
      module procedure mpi_reduce_inplace_sum_wp
      module procedure mpi_reduce_inplace_sum_ep
   end interface mpi_reduce_inplace_sum_cfp

   interface mpi_reduce_sum
      module procedure mpi_reduce_sum_wp
      module procedure mpi_reduce_sum_ep
   end interface mpi_reduce_sum

   interface mpi_reduceall_sum_cfp
      module procedure mpi_reduceall_sum_wp
      module procedure mpi_reduceall_sum_ep
   end interface mpi_reduceall_sum_cfp

   interface naive_mpi_reduce_inplace_sum
      module procedure naive_mpi_reduce_inplace_sum_wp
      module procedure naive_mpi_reduce_inplace_sum_ep
   end interface naive_mpi_reduce_inplace_sum

   interface mpi_reduceall_inplace_sum_cfp
      module procedure mpi_reduceall_inplace_sum_wp
      module procedure mpi_reduceall_inplace_sum_ep
   end interface mpi_reduceall_inplace_sum_cfp

   interface mpi_mod_rotate_arrays_around_ring
      module procedure mpi_mod_rotate_arrays_around_ring_wp
      module procedure mpi_mod_rotate_arrays_around_ring_ep
   end interface mpi_mod_rotate_arrays_around_ring

   interface mpi_mod_rotate_cfp_arrays_around_ring
      module procedure mpi_mod_rotate_wp_arrays_around_ring
      module procedure mpi_mod_rotate_ep_arrays_around_ring
   end interface mpi_mod_rotate_cfp_arrays_around_ring

   !> \brief   Convenience MPI type wrapper with garbage collection
   !> \authors J Benda
   !> \date    2019
   !>
   !> Derived type that wraps construction and destruction of a MPI type describing a blocked array.
   !>
   type :: MPIBlockArray
      integer(mpiint) :: mpitype = 0
   contains
      procedure :: setup => mpi_mod_create_block_array
      final     :: mpi_mod_free_block_array
   end type MPIBlockArray

contains

   !> \brief   Comm or default comm
   !> \authors J Benda
   !> \date    2019
   !>
   !> Auxiliary function used in other subroutines for brevity. If an argument is given, it will return its value. If no
   !> argument is given, the function returns MPI_COMM_WORLD.
   !>
   integer(mpiint) function mpi_mod_comm (comm_opt) result (comm)

      integer(mpiint), intent(in), optional :: comm_opt
#ifdef usempi
      if (present(comm_opt)) then
         comm = comm_opt
      else
         comm = MPI_COMM_WORLD
      end if
#else
      comm = 0
#endif
   end function mpi_mod_comm


   !> \brief   Wrapper around MPI_Comm_rank
   !> \authors J Benda
   !> \date    2019
   !>
   !> Wrapper around a standard MPI rank query which returns zero when compiled
   !> without MPI.
   !>
   subroutine mpi_mod_rank (rank, comm)

      integer(mpiint),           intent(out) :: rank
      integer(mpiint), optional, intent(in)  :: comm
#ifdef usempi
      integer(mpiint) :: ierr

      call MPI_Comm_rank(mpi_mod_comm(comm), rank, ierr)
#else
      rank = master
#endif
   end subroutine mpi_mod_rank


   !> Interface to the routine MPI_BARRIER.
   subroutine mpi_mod_barrier (error, comm)

      integer,         intent(out)          :: error
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr

      call MPI_Barrier(mpi_mod_comm(comm), ierr)
      error = ierr
#else
      error = 0
#endif
   end subroutine mpi_mod_barrier


   !> Analogue of the xermsg routine. This routine uses xermsg to first print out the error message and then either aborts 
   !> the program or exits the routine depending on the level of the message (warning/error) which is defined in the same way as in xermsg.
   subroutine mpi_xermsg(mod_name,routine_name,err_msg,err,level)
      implicit none
      character(len=*), intent(in) :: mod_name,routine_name,err_msg
      integer, intent(in) :: err, level
#ifdef usempi
      integer(kind=mpiint) :: ierr, error

         if (level > 0) then
            !This is an error that should cause abort but we want to abort MPI
            !in a clean way: by -1 we force xermsg to print an error message but
            !don't abort using the stop command. We abort below using MPI_ABORT.
            call xermsg(mod_name,routine_name,err_msg,err,-1)
         else
            call xermsg(mod_name,routine_name,err_msg,err,level)
         endif
         error = err
         if (level .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
#else
         call xermsg(mod_name,routine_name,err_msg,err,level)
#endif

   end subroutine mpi_xermsg

   !> This is a lightweight routine which aborts the program if MPI is not found running.
   subroutine check_mpi_running
      implicit none

#ifdef usempi
         if (.not.(mpi_running)) call xermsg('mpi_mod','mpi_mod_start','MPI has not been initialized. The program will abort.',1,1)
#endif

   end subroutine check_mpi_running


   !> \brief   Display information about current MPI setup
   !> \authors Z Masin, J Benda
   !> \date    2014 - 2019
   !>
   !> Originally part of \ref mpi_mod_start, but separated for better application-level control of output from the library.
   !>
   !> \param[in] u  Unit for the text output (mostly stdout).
   !>
   subroutine mpi_mod_print_info (u)

      integer, intent(in) :: u

      write (u, *)

      if (.not. mpi_running) then
         write (u, '(1X,"MPI not running.")')
      else
         write (u, '(1X,"MPI running with ",I0," processes.")') nprocs
         write (u, '(1X,"My rank is ",I0)') myrank
         write (u, '(1X,"I am running on the processor with name: ",A)') procname

         call print_precision_params(stdout)

         if (shared_enabled) then
            write (u, *)
            write (u, '(1X,"Using MPI-3 shared memory")')
            write (u, '(1X,"Number of MPI tasks on this node ",I0)') local_nprocs
            write (u, '(1X,"My rank on this node is ",I0)') local_rank
         end if
      end if

      write (u, *)

   end subroutine mpi_mod_print_info


   !> Initializes MPI, assigns the rank for each process and finds out how many processors we're using. It also maps the current
   !> floating point precision (kind=cfp) to the corresponding MPI numeric type. This routine also sets the unit for standard
   !> output that each rank will use (the value of stdout). In case of serial run the value stdout is not set here and is left
   !> to the default value input_unit as specified in the module const.
   !> \warning This must the first statement in every level3 program. 
   subroutine mpi_mod_start (do_stdout, allow_shared_memory)
      implicit none
      logical, optional, intent(in) :: do_stdout, allow_shared_memory
      logical :: do_stdout_
#ifdef usempi
      integer(kind=mpiint) :: ierr=0, error=0, isize=0, zero=0
      integer :: myint, imyrank
      character(len=MPI_MAX_ERROR_STRING) :: estring
      integer(kind=mpiint) :: error_class,error2,error_length

        do_stdout_ = .false.
        if(present(do_stdout)) do_stdout_ = do_stdout


         if (mpi_running) then
   
            call mpi_xermsg('mpi_mod', 'mpi_mod_start', &
                            'Attempt to start MPI while it has been initialized before. The program will abort.', 1, 1)
   
         else

            !Initialize MPI without threading
            call MPI_INIT(error)

            if (error .ne. MPI_SUCCESS) then
               call xermsg('mpi_mod','mpi_mod_start','MPI initialization has failed. The program will abort.',2,0)
               call MPI_Error_class(error, error_class,error2);
               call MPI_Error_string(error, estring,error_length,error2);
               write(*,*) estring(1:error_length)
               write(*,*) 'Error class',error_class
               call MPI_ABORT(MPI_COMM_WORLD,ierr,error)
            else

               !Determine basic properties of the parallelism: number of processes and rank of my MPI process.
               call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,error)
               if (error .ne. MPI_SUCCESS) then
                  call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'MPI_COMM_SIZE has failed.', 3, 1)
               end if

               call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,error)
               if (error .ne. MPI_SUCCESS) then
                  call mpi_xermsg('mpi_mod', 'mpi_mod_start','MPI_COMM_RANK has failed.', 5, 1)
               end if

               !This is where we redirect all standard output to a file. By standard output we mean all output using the write statement where the unit number used is 'stdout'.
               !Following the call to redirect_stdout_to_file the value of the variable const/stdout will be associated with a file intended for
               !standard output for the process with rank=myrank. Finally, remember that redirect_stdout_to_file can be called only once!
               imyrank = myrank
               call redirect_stdout_to_file(imyrank, do_stdout_)

               call MPI_GET_PROCESSOR_NAME(procname,isize,error)
               if (error .ne. MPI_SUCCESS) then
                  call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'MPI_GET_PROCESSOR_NAME has failed.', 6, 1)
               end if

               mpi_mod_wp = MPI_DOUBLE_PRECISION
               mpi_mod_ep = MPI_REAL16
   
               !Determine the MPI type corresponding to the default real floating point numbers.
               if (cfp .eq. wp) then
                  mpi_mod_cfp = MPI_DOUBLE_PRECISION
               elseif (cfp .eq. ep1) then
                  mpi_mod_cfp = MPI_REAL16
               else
                  call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'Using an unsupported floating point type.', 7, 1)
               end if
   
               !Type of the default integer.
               isize = bit_size(myint)/8
               call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER,isize,mpi_mod_int,ierr)

               !Type of the MPI integer.
               mpi_mod_mpiint = MPI_INTEGER
   
               mpi_started = .true.
               mpi_running = .true.
#   ifdef mpithree
               shared_enabled = .false.

               if (present(allow_shared_memory)) then
                  shared_enabled = allow_shared_memory
               end if

               if (shared_enabled) then
                  call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, zero, MPI_INFO_NULL, shared_communicator, error)
                  call MPI_Comm_size(shared_communicator, local_nprocs, error)
                  call MPI_Comm_rank(shared_communicator, local_rank, error)
               else
#   endif
                  shared_enabled = .false.
                  shared_communicator = MPI_COMM_WORLD
                  local_nprocs = nprocs
                  local_rank = myrank
#   ifdef mpithree
               end if
#   endif

            endif
   
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#else
        do_stdout_ = .false.
        if(present(do_stdout)) do_stdout_ = do_stdout

         call print_precision_params(stdout)

         nprocs = 1
         myrank = 0

         !This is where we redirect all standard output to a file. Remember that redirect_stdout_to_file can be called only once!
         if (.not. (mpi_started)) call redirect_stdout_to_file(myrank,do_stdout_)

         mpi_started = .true.
#endif

   end subroutine mpi_mod_start

   !> It is not used anywhere in the code yet but it might become useful later on. It creates intra-node communicators without the memory sharing capability.
   subroutine setup_intranode_local_mem_communicator
      implicit none
#ifdef usempi
      integer, parameter :: nmaxp = MPI_MAX_PROCESSOR_NAME
      character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: pnames(:)
      integer(kind=mpiint), allocatable :: colours(:)
      logical :: match      
      integer :: n,i,j
      integer(kind=mpiint) :: error,ierr,zero=0

         !Generate a new intra-node communicator
         allocate(pnames(nprocs),colours(nprocs))
         call mpi_mod_allgather(procname,pnames)

         !Find the tasks that sit on the node with the same name
         n = 1
         colours = -1
         do i=1,nprocs
            match = .false.
            do j=1,nprocs
               if (strings_are_same(pnames(i),pnames(j),nmaxp) .and. colours(j) .eq. -1) then
                  colours(j) = n
                  match = .true.
               endif
            enddo
            if (match) n = n + 1
         enddo

         do i=1,nprocs
            if (strings_are_same(pnames(i),procname,nmaxp)) node_colour = colours(i)
         enddo

         call MPI_COMM_SPLIT(MPI_COMM_WORLD,node_colour,zero,shared_communicator,ierr)
         if (ierr .ne. MPI_SUCCESS) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'MPI_COMM_SPLIT has failed.', 8, 1)
         end if

         !Determine basic properties of the parallelism for the intra-node communicator: number of processes and rank of the MPI process.
         call MPI_COMM_SIZE(shared_communicator,local_nprocs,error)
         if (error .ne. MPI_SUCCESS) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'MPI_COMM_SIZE has failed.', 9, 1)
         end if

         call MPI_COMM_RANK(shared_communicator,local_rank,error)
         if (error .ne. MPI_SUCCESS) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'MPI_COMM_RANK has failed.', 10, 1)
         end if

         write(stdout,'(/,10X,"Total number of different nodes ",i0)') maxval(colours)
         write(stdout,'(10X,"Node colour ",i0)') node_colour
         write(stdout,'(10X,"Number of MPI tasks on this node ",i0)') local_nprocs
         write(stdout,'(10X,"My rank on this node is ",i0)') local_rank
#endif
   end subroutine setup_intranode_local_mem_communicator

   function strings_are_same(str1,str2,length)
      implicit none
      integer, intent(in) :: length
      character(len=length), intent(in) :: str1, str2
      logical :: strings_are_same

      integer :: i

         strings_are_same = .true.

         do i=1,length
            if (str1(i:i) .ne. str2(i:i)) strings_are_same = .false.
         enddo

   end function strings_are_same

   !> Terminates the MPI session and stops the program. It is a blocking routine.
   subroutine mpi_mod_finalize
      implicit none
#ifdef usempi
      integer(kind=mpiint) :: ierr, error
      integer :: imyrank

         call check_mpi_running

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         imyrank = myrank
         call report_statistics_error_messages(imyrank)

         call MPI_FINALIZE(error)

         if (error .ne. MPI_SUCCESS) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_start', 'MPI finalization has failed.', 1, 1)
         end if
#endif
         write(stdout,'("<-------->","done:mpi_mod:mpi_mod_finalize")')

         myrank = -1
         nprocs = -1

         mpi_running = .false.

         stop

   end subroutine mpi_mod_finalize


   !> \brief   Set up a new MPI type describing a long array composed of blocks
   !> \authors J Benda
   !> \date    2019
   !>
   !> Creates a MPI type handle describing a structured data type that can be used to transfer very long arrays of data.
   !> The array is represented as a set of blocks, whose element count does not exceed the maximal value of the MPI
   !> integer type. To avoid resource leaks, the returned handle should be released by a call to `MPI_Type_free`.
   !>
   !> \param[in] this      Block array object to initialize.
   !> \param[in] n         Number of elements in the array.
   !> \param[in] elemtype  MPI datatype handle of elements of the array.
   !>
   subroutine mpi_mod_create_block_array (this, n, elemtype)

      use iso_c_binding, only: c_int

      class(MPIBlockArray), intent(inout) :: this

      integer,         intent(in) :: n
      integer(mpiint), intent(in) :: elemtype
#ifdef usempi
      integer, parameter :: blocksize = huge(0_c_int)  ! largest array addressable by C integer (standard MPI interface)

      integer(mpiint),  allocatable :: blocklengths(:), types(:)
      integer(mpiaddr), allocatable :: offsets(:)
      integer(mpiint)               :: nblocks, lastblock, ierr, i
      integer(mpicnt)               :: elembytes

      nblocks = n / blocksize
      lastblock = mod(n, blocksize)
      if (lastblock /= 0) nblocks = nblocks + 1
      allocate (blocklengths(nblocks), offsets(nblocks), types(nblocks))

      ! determine byte size of a single element of the array
      call MPI_Type_size_x(elemtype, elembytes, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_create_block_array', 'Failed determine byte size of MPI type.', int(ierr), 1)
      end if

      ! elements of all blocks have the same type
      types(:) = elemtype

      ! all blocks have the same size, except for the last one (unless the array length is a multiple of 'blocksize')
      blocklengths(:) = blocksize
      if (lastblock /= 0) then
        blocklengths(nblocks) = lastblock
      end if

      ! set up block offsets in bytes
      offsets(:) = 0
      do i = 1, nblocks - 1
         offsets(i + 1) = offsets(i) + blocklengths(i) * int(elembytes, mpiaddr)
      end do

      ! create the block array MPI derived type
      call MPI_Type_create_struct(nblocks, blocklengths, offsets, types, this % mpitype, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_create_block_array', 'Failed to create blocked array MPI type.', int(ierr), 1)
      end if

      ! finalize the derived type for use
      call MPI_Type_commit(this % mpitype, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_create_block_array', 'Failed to commit blocked array MPI type.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_create_block_array


   !> \brief   Release MPI handle of the block array type
   !> \authors J Benda
   !> \date    2019
   !>
   !> This is a destructor for MPIBlockArray. It frees the handle for the MPI derived type (if any).
   !>
   subroutine mpi_mod_free_block_array (this)

      type(MPIBlockArray), intent(inout) :: this
#ifdef usempi
      integer(mpiint) :: ierr

      if (this % mpitype /= 0) then
         call MPI_Type_free(this % mpitype, ierr)
      end if
#endif
   end subroutine mpi_mod_free_block_array


   subroutine mpi_mod_bcast_logical (val, from, comm)

      logical,         intent(inout)        :: val
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr, one = 1

      if (nprocs == 1) return

      call MPI_Bcast(val, one, MPI_LOGICAL, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_logical', 'MPI_BCAST has failed.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_logical


   subroutine mpi_mod_bcast_int32 (val, from, comm)

      integer(int32),  intent(inout)        :: val
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr, one = 1

      if (nprocs == 1) return

      call MPI_Bcast(val, one, MPI_INTEGER4, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_int32', 'MPI_BCAST has failed.', int(ierr), 1)
      end if
#endif

   end subroutine mpi_mod_bcast_int32


   subroutine mpi_mod_bcast_int64 (val, from, comm)

      integer(int64),  intent(inout)        :: val
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr, one = 1

      if (nprocs == 1) return

      call MPI_Bcast(val, one, MPI_INTEGER8, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_int64', 'MPI_BCAST has failed.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_int64


   subroutine mpi_mod_bcast_int32_array (val, from, comm)

      integer(int32),  intent(inout)        :: val(:)
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      if (nprocs == 1 .or. size(val) == 0) return

      call blockarray % setup(size(val), MPI_INTEGER4)
      call MPI_Bcast(val, one, blockarray % mpitype, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_int32_array', 'Failed to broadcast int32 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_int32_array


   subroutine mpi_mod_bcast_int64_array (val, from, comm)

      integer(int64),  intent(inout)        :: val(:)
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      if (nprocs == 1 .or. size(val) == 0) return

      call blockarray % setup(size(val), MPI_INTEGER8)
      call MPI_Bcast(val, one, blockarray % mpitype, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_int64_array', 'Failed to broadcast int64 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_int64_array


   subroutine mpi_mod_bcast_wp (val, from, comm)

      real(wp),        intent(inout)        :: val
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(kind=mpiint) :: ierr, one = 1

      if (nprocs == 1) return

      call MPI_Bcast(val, one, mpi_mod_wp, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_wp', 'MPI_BCAST has failed.', int(ierr), 0)
      end if
#endif
   end subroutine mpi_mod_bcast_wp


   subroutine mpi_mod_bcast_ep (val, from, comm)

      real(ep),        intent(inout)        :: val
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(kind=mpiint) :: ierr, one = 1

      if (nprocs == 1) return

      call MPI_Bcast(val, one, mpi_mod_ep, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_ep', 'MPI_BCAST has failed.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_ep


   subroutine mpi_mod_bcast_wp_array (val, from, comm)

      real(wp),        intent(inout)        :: val(:)
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      if (nprocs == 1 .or. size(val) == 0) return

      call blockarray % setup(size(val), mpi_mod_wp)
      call MPI_Bcast(val, one, blockarray % mpitype, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_wp_array', 'Failed to broadcast wp array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_wp_array


   subroutine mpi_mod_bcast_ep_array (val, from, comm)

      real(ep),        intent(inout)        :: val(:)
      integer(mpiint), intent(in)           :: from
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      if (nprocs == 1 .or. size(val) == 0) return

      call blockarray % setup(size(val), mpi_mod_ep)
      call MPI_Bcast(val, one, blockarray % mpitype, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_ep_array', 'Failed to broadcast ep array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_ep_array


   subroutine mpi_mod_bcast_character (val, from, comm)

      character(len=*), intent(inout)        :: val
      integer(mpiint),  intent(in)           :: from
      integer(mpiint),  intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr, l

      if (nprocs == 1) return

      l = len(val)
      call MPI_Bcast(val, l, MPI_CHARACTER, from, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_bcast_character', 'MPI_BCAST has failed.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_character


   subroutine mpi_mod_bcast_character_array (val, from, comm)

      character(len=*), intent(inout)        :: val(:)
      integer(mpiint),  intent(in)           :: from
      integer(mpiint),  intent(in), optional :: comm
#ifdef usempi
      integer             :: n
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      if (nprocs == 1 .or. size(val) == 0) return

      n = len(val(1)) * size(val)

      call blockarray % setup(n, MPI_CHARACTER)
      call MPI_Bcast(val, one, blockarray % mpitype, from, mpi_mod_comm(comm), ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_send_character_array', 'Failed to send character array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_bcast_character_array


   subroutine mpi_mod_send_wp_array (to, buffer, tag, n, comm)

      integer(mpiint), intent(in) :: to
      real(wp),        intent(in) :: buffer(:)
      integer,         intent(in) :: tag, n
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, mpi_mod_wp)
      call MPI_Send(buffer, one, blockarray % mpitype, to, int(tag, mpiint), mpi_mod_comm(comm), ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_send_wp_array', 'Failed to send wp array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_send_wp_array


   subroutine mpi_mod_send_ep_array (to, buffer, tag, n, comm)

      integer(mpiint), intent(in) :: to
      real(ep),        intent(in) :: buffer(:)
      integer,         intent(in) :: tag, n
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, mpi_mod_ep)
      call MPI_Send(buffer, one, blockarray % mpitype, to, int(tag, mpiint), mpi_mod_comm(comm), ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_send_ep_array', 'Failed to send ep array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_send_ep_array


   subroutine mpi_mod_send_int32_array (to, buffer, tag, n, comm)

      integer(mpiint), intent(in) :: to
      integer(int32),  intent(in) :: buffer(:)
      integer,         intent(in) :: tag, n
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, MPI_INTEGER4)
      call MPI_Send(buffer, one, blockarray % mpitype, to, int(tag, mpiint), mpi_mod_comm(comm), ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_send_int32_array', 'Failed to send int32 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_send_int32_array


   subroutine mpi_mod_send_int64_array (to, buffer, tag, n, comm)

      integer(mpiint), intent(in) :: to
      integer(int64),  intent(in) :: buffer(:)
      integer,         intent(in) :: tag, n
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, MPI_INTEGER8)
      call MPI_Send(buffer, one, blockarray % mpitype, to, int(tag, mpiint), mpi_mod_comm(comm), ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_send_int64_array', 'Failed to send int64 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_send_int64_array


   subroutine mpi_mod_isend_int32_array (to, buffer, tag, n, comm)

      integer(mpiint), intent(in) :: to
      integer(int32),  intent(in) :: buffer(:)
      integer,         intent(in) :: tag, n
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1, req
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, MPI_INTEGER4)
      call MPI_Isend(buffer, one, blockarray % mpitype, to, int(tag, mpiint), mpi_mod_comm(comm), req, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_isend_int32_array', 'Failed to isend int32 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_isend_int32_array


   subroutine mpi_mod_isend_int64_array (to, buffer, tag, n, comm)

      integer(mpiint), intent(in) :: to
      integer(int64),  intent(in) :: buffer(:)
      integer,         intent(in) :: tag, n
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1, req
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, MPI_INTEGER8)
      call MPI_Isend(buffer, one, blockarray % mpitype, to, int(tag, mpiint), mpi_mod_comm(comm), req, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_isend_int64_array', 'Failed to isend int64 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_isend_int64_array


   subroutine mpi_mod_recv_wp_array (from, tag, buffer, n, comm)

      integer(mpiint), intent(in)  :: from
      integer,         intent(in)  :: tag
      integer,         intent(in)  :: n
      real(wp),        intent(out) :: buffer(:)
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, mpi_mod_wp)
      call MPI_Recv(buffer, one, blockarray % mpitype, from, int(tag, mpiint), mpi_mod_comm(comm), MPI_STATUS_IGNORE, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_recv_wp_array', 'Failed to recv wp array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_recv_wp_array


   subroutine mpi_mod_recv_ep_array (from, tag, buffer, n, comm)

      integer(mpiint), intent(in)  :: from
      integer,         intent(in)  :: tag
      integer,         intent(in)  :: n
      real(ep),        intent(out) :: buffer(:)
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, mpi_mod_ep)
      call MPI_Recv(buffer, one, blockarray % mpitype, from, int(tag, mpiint), mpi_mod_comm(comm), MPI_STATUS_IGNORE, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_recv_ep_array', 'Failed to recv ep array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_recv_ep_array


   subroutine mpi_mod_recv_int32_array (from, tag, buffer, n, comm)

      integer(mpiint), intent(in)  :: from
      integer,         intent(in)  :: tag
      integer,         intent(in)  :: n
      integer(int32),  intent(out) :: buffer(:)
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, MPI_INTEGER4)
      call MPI_Recv(buffer, one, blockarray % mpitype, from, int(tag, mpiint), mpi_mod_comm(comm), MPI_STATUS_IGNORE, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_recv_int32_array', 'Failed to recv int32 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_recv_int32_array


   subroutine mpi_mod_recv_int64_array (from, tag, buffer, n, comm)

      integer(mpiint), intent(in)  :: from
      integer,         intent(in)  :: tag
      integer,         intent(in)  :: n
      integer(int64),  intent(out) :: buffer(:)
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint)     :: ierr, one = 1
      type(MPIBlockArray) :: blockarray

      call blockarray % setup(n, MPI_INTEGER8)
      call MPI_Recv(buffer, one, blockarray % mpitype, from, int(tag, mpiint), mpi_mod_comm(comm), MPI_STATUS_IGNORE, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_recv_int64_array', 'Failed to recv int64 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_recv_int64_array


   subroutine mpi_mod_send_wp_array_dynamic (to, buffer, tag, comm)

      integer(mpiint), intent(in) :: to
      real(wp),        intent(in) :: buffer(:)
      integer,         intent(in) :: tag
      integer(mpiint), intent(in), optional :: comm

      call mpi_mod_send(to, buffer, tag, size(buffer), comm)

   end subroutine mpi_mod_send_wp_array_dynamic


   subroutine mpi_mod_send_ep_array_dynamic (to, buffer, tag, comm)

      integer(mpiint), intent(in) :: to
      real(ep),        intent(in) :: buffer(:)
      integer,         intent(in) :: tag
      integer(mpiint), intent(in), optional :: comm

      call mpi_mod_send(to, buffer, tag, size(buffer), comm)

   end subroutine mpi_mod_send_ep_array_dynamic


   subroutine mpi_mod_recv_wp_array_dynamic (from, tag, buffer, n, comm)

      integer(mpiint), intent(in)  :: from
      integer,         intent(in)  :: tag
      integer,         intent(out) :: n
      integer(mpiint), intent(in),  optional    :: comm
      real(wp),        intent(out), allocatable :: buffer(:)
#ifdef usempi
      integer             :: error
      integer(mpiint)     :: ierr, one = 1, stat(MPI_STATUS_SIZE)
      integer(mpicnt)     :: sz
      type(MPIBlockArray) :: blockarray

      call MPI_Probe(from, int(tag, mpiint), mpi_mod_comm(comm), stat, ierr)
      call MPI_Get_elements_x(stat, mpi_mod_wp, sz, ierr)

      n = sz
      if (allocated(buffer)) deallocate (buffer)
      allocate (buffer(1:n), stat = error)
      if (error /= 0) call mpi_xermsg('mpi_mod', 'mpi_mod_recv_wp_array_dynamic', 'Memory allocation failed.', error, 1)

      call blockarray % setup(n, mpi_mod_wp)
      call MPI_Recv(buffer, one, blockarray % mpitype, from, int(tag, mpiint), mpi_mod_comm(comm), stat, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_recv_wp_array_dynamic', 'Failed to recv int64 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_recv_wp_array_dynamic


   subroutine mpi_mod_recv_ep_array_dynamic (from, tag, buffer, n, comm)

      integer(mpiint), intent(in)  :: from
      integer,         intent(in)  :: tag
      integer,         intent(out) :: n
      integer(mpiint), intent(in),  optional    :: comm
      real(ep),        intent(out), allocatable :: buffer(:)
#ifdef usempi
      integer             :: error
      integer(mpiint)     :: ierr, one = 1, stat(MPI_STATUS_SIZE)
      integer(mpicnt)     :: sz
      type(MPIBlockArray) :: blockarray

      call MPI_Probe(from, int(tag, mpiint), mpi_mod_comm(comm), stat, ierr)
      call MPI_Get_elements_x(stat, mpi_mod_ep, sz, ierr)

      n = sz
      if (allocated(buffer)) deallocate (buffer)
      allocate (buffer(1:n), stat = error)
      if (error /= 0) call mpi_xermsg('mpi_mod', 'mpi_mod_recv_ep_array_dynamic', 'Memory allocation failed.', error, 1)

      call blockarray % setup(n, mpi_mod_ep)
      call MPI_Recv(buffer, one, blockarray % mpitype, from, int(tag, mpiint), mpi_mod_comm(comm), stat, ierr)

      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_recv_ep_array_dynamic', 'Failed to recv int64 array.', int(ierr), 1)
      end if
#endif
   end subroutine mpi_mod_recv_ep_array_dynamic


   subroutine mpi_mod_file_open_read (filename, fh, ierr, comm)

      character(len=*), intent(in)           :: filename
      integer(mpiint),  intent(out)          :: fh
      integer,          intent(out)          :: ierr
      integer(mpiint),  intent(in), optional :: comm
#ifdef usempi      
      integer(mpiint) :: error

      call MPI_File_open(mpi_mod_comm(comm), filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, error)
      ierr = error
#else
      open (newunit = fh, file = filename, status = 'old', action = 'read', &
            form = 'unformatted', access = 'stream', position = 'rewind', iostat = ierr)
#endif
   end subroutine mpi_mod_file_open_read


   !> Sets view for a file containing floating point numbers with kind=wp.
   subroutine mpi_mod_file_set_view_wp(fh, disp, ierr, wp_dummy)
      implicit none
      integer, intent(in) :: fh, disp
      integer, intent(out) :: ierr
      real(kind=wp), intent(in) :: wp_dummy
#ifdef usempi
      integer(kind=MPI_OFFSET_KIND) :: disp_conv
      integer(kind=mpiint) :: fh_conv, ierr_conv, etype, filetype, error
      character(len=*), parameter :: datarep = 'native'

         fh_conv = fh
         disp_conv = disp
         if (fh_conv .ne. fh .or. disp_conv .ne. disp) then
            print *,fh_conv,fh,disp_conv,disp
            call mpi_xermsg('mpi_mod', 'mpi_mod_file_set_view_wp', &
                            'At least one of fh,disp arguments are too large for the mpi integer kind.', 1, 1)
         endif
         etype = mpi_mod_wp
         filetype = mpi_mod_wp
         call MPI_FILE_SET_VIEW(fh_conv, disp_conv, etype, filetype, datarep, MPI_INFO_NULL, ierr_conv)
         ierr = ierr_conv
         if (ierr /= MPI_SUCCESS) then
            call mpi_xermsg ('mpi_mod', 'mpi_mod_file_set_view_wp', 'MPI_FILE_SET_VIEW failed with an error message.', ierr, 1)
         end if
#else
         call xermsg('mpi_mod','mpi_mod_file_set_view_wp','Not yet implemented for serial compilation.',1,1)
#endif

   end subroutine mpi_mod_file_set_view_wp


   !> Sets view for a file containing floating point numbers with kind=ep.
   subroutine mpi_mod_file_set_view_ep(fh, disp, ierr, ep_dummy)
      implicit none
      integer, intent(in) :: fh, disp
      integer, intent(out) :: ierr
      real(kind=ep), intent(in) :: ep_dummy
#ifdef usempi
      integer(kind=MPI_OFFSET_KIND) :: disp_conv
      integer(kind=mpiint) :: fh_conv, ierr_conv, etype, filetype, error
      character(len=*), parameter :: datarep = 'native'

         fh_conv = fh
         disp_conv = disp
         if (fh_conv .ne. fh .or. disp_conv .ne. disp) then
            print *,fh_conv,fh,disp_conv,disp
            call mpi_xermsg('mpi_mod', 'mpi_mod_file_set_view_ep', &
                            'At least one of fh,disp arguments are too large for the mpi integer kind.', 1, 1)
         endif
         etype = mpi_mod_ep
         filetype = mpi_mod_ep
         call MPI_FILE_SET_VIEW(fh_conv, disp_conv, etype, filetype, datarep, MPI_INFO_NULL, ierr_conv)
         ierr = ierr_conv
         if (ierr /= MPI_SUCCESS) then
            call mpi_xermsg ('mpi_mod', 'mpi_mod_file_set_view_ep', 'MPI_FILE_SET_VIEW failed with an error message.', ierr, 1)
         end if
#else
         call xermsg('mpi_mod','mpi_mod_file_set_view_ep','Not yet implemented for serial compilation.',1,1)
#endif

   end subroutine mpi_mod_file_set_view_ep

   subroutine mpi_mod_file_read_wp(fh,buffer,buflen,ierr)
      implicit none
      integer, intent(in) :: fh, buflen
      integer, intent(out) :: ierr
      real(kind=wp), intent(out) :: buffer(:)
#ifdef usempi
      integer(kind=mpiint) :: fh_conv, ierr_conv, buflen_conv, etype, error, test

         etype = mpi_mod_wp
         fh_conv = fh
         buflen_conv = buflen
         if (fh_conv .ne. fh .or. buflen_conv .ne. buflen) then
            print *,fh_conv,fh,buflen_conv,buflen
            call mpi_xermsg('mpi_mod', 'mpi_mod_file_read_wp', &
                            'At least one of fh,buflen arguments are too large for the mpi integer kind.', 1, 1)
         endif

         test = buflen_conv*wp_bytes
         if (test .ne. buflen_conv*wp_bytes) then
            call mpi_xermsg ('mpi_mod', 'mpi_mod_file_read_wp', &
                             'The array to be read-in is too large to be read-in using MPI_FILE_READ.', 2, 1)
         endif
         
         call MPI_FILE_READ(fh_conv, buffer, buflen_conv, etype, MPI_STATUS_IGNORE, ierr_conv)
         ierr = ierr_conv
         if (ierr /= MPI_SUCCESS) then
            call mpi_xermsg ('mpi_mod', 'mpi_mod_file_read_wp', 'MPI_FILE_READ failed with an error message.', ierr, 1)
         end if
#else
         call xermsg('mpi_mod','mpi_mod_file_read_wp','Not yet implemented for serial compilation.',1,1)
#endif

   end subroutine mpi_mod_file_read_wp

   subroutine mpi_mod_file_read_ep(fh,buffer,buflen,ierr)
      implicit none
      integer, intent(in) :: fh, buflen
      integer, intent(out) :: ierr
      real(kind=ep1), intent(out) :: buffer(:)
#ifdef usempi
      integer(kind=mpiint) :: fh_conv, ierr_conv, buflen_conv, etype, error, test

         etype = mpi_mod_ep
         fh_conv = fh
         buflen_conv = buflen
         if (fh_conv .ne. fh .or. buflen_conv .ne. buflen) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_file_read_ep', &
                            'At least one of fh,buflen arguments are too large for the mpi integer kind.', 1, 1)
         end if

         test = buflen_conv*ep1_bytes
         if (test .ne. buflen_conv*ep1_bytes) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_file_read_ep', &
                            'The array to be read-in is too large to be read-in using MPI_FILE_READ.', 2, 1)
         end if

         call MPI_FILE_READ(fh_conv, buffer, buflen_conv, etype, MPI_STATUS_IGNORE, ierr_conv)
         ierr = ierr_conv
         if (ierr /= MPI_SUCCESS) then
            call mpi_xermsg('mpi_mod', 'mpi_mod_file_read_ep', 'MPI_FILE_READ failed with an error message.', ierr, 1)
         end if
#else
         call xermsg('mpi_mod','mpi_mod_file_read_ep','Not yet implemented for serial compilation.',1,1)
#endif
   end subroutine mpi_mod_file_read_ep


   subroutine mpi_mod_file_close (fh, ierr)

      integer(mpiint), intent(inout) :: fh
      integer,         intent(out)   :: ierr
#ifdef usempi
      integer(mpiint) :: error

      call MPI_File_close(fh, error)
      ierr = error
#else
      close (fh)
#endif
   end subroutine mpi_mod_file_close


   subroutine mpi_reduce_inplace_sum_wp(buffer,nelem)
      use omp_lib
      implicit none
      real(kind=wp), intent(inout) :: buffer(:)
      integer, intent(in) :: nelem
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error, block_size
      integer :: n_times_reduce, i, first, last
      real(kind=wp) :: start_t, end_t, total_start, total_end
         if(nprocs <= 1) return
         write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem) then
            print *,nelem_conv,nelem
            call mpi_xermsg('mpi_mod', 'mpi_mod_reduce_inplace_sum_wp', &
                            'The nelem argument is too large for the mpi integer kind. The program will abort.', 1, 1)
         end if

         !MPI_REDUCE can handle at most max_data_count array elements to reduce
         !so we need to determine how many times we need to call it to reduce
         !the whole array which can have more than max_data_count array
         !elements. 
         n_times_reduce = ceiling(nelem_conv/(max_data_count_wp*1.0_wp))

         block_size = nelem_conv/n_times_reduce 

         !In each iteration we reduce the elements buffer(first:last).
         first = 0 
         last = 0
         do i=1,n_times_reduce
            if (i .eq. n_times_reduce) block_size = nelem_conv - last !reduce what's left 
            first = last + 1
            last = last + block_size
            if (n_times_reduce > 1) then
               write(stdout,'("Reducing elements: ",i0," to ",i0)') first, last
            endif
            start_t = omp_get_wtime()
            write(stdout,"('Start reduce')")
            if (myrank .eq. master) then
                write(stdout,"('Start reduce master')")
               call MPI_REDUCE(MPI_IN_PLACE, buffer(first:last), block_size, mpi_mod_wp, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            else
                write(stdout,"('Start reduce slave')")
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_REDUCE(buffer(first:last), block_size,   block_size, mpi_mod_wp, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            endif
            end_t = omp_get_wtime()
            write(stdout,'("Reduction took: ",f8.3," [s]")') end_t-start_t
         enddo !i

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && !defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         call naive_mpi_reduce_inplace_sum_wp(buffer,nelem,.false.)
#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         !Additionally we split the reduction into separate reductions each not larger than ~2GiB.
         call naive_mpi_reduce_inplace_sum_wp(buffer,nelem,.true.)
#endif

   end subroutine mpi_reduce_inplace_sum_wp

   subroutine mpi_reduce_inplace_sum_ep(buffer,nelem)
      use omp_lib
      implicit none
      real(kind=ep), intent(inout) :: buffer(:)
      integer, intent(in) :: nelem
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error, block_size
      integer :: n_times_reduce, i, first, last
      real(kind=ep) :: start_t, end_t, total_start, total_end
         if(nprocs <= 1) return
         write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem) then
            print *,nelem_conv,nelem
            call mpi_xermsg('mpi_mod', 'mpi_mod_reduce_inplace_sum_ep', &
                            'The nelem argument is too large for the mpi integer kind. The program will abort.', 1, 1)
         end if

         !MPI_REDUCE can handle at most max_data_count array elements to reduce
         !so we need to determine how many times we need to call it to reduce
         !the whole array which can have more than max_data_count array
         !elements. 
         n_times_reduce = ceiling(nelem_conv/(max_data_count_ep*1.0_ep))

         block_size = nelem_conv/n_times_reduce 

         !In each iteration we reduce the elements buffer(first:last).
         first = 0 
         last = 0
         do i=1,n_times_reduce
            if (i .eq. n_times_reduce) block_size = nelem_conv - last !reduce what's left 
            first = last + 1
            last = last + block_size
            if (n_times_reduce > 1) then
               write(stdout,'("Reducing elements: ",i0," to ",i0)') first, last
            endif
            start_t = omp_get_wtime()
            write(stdout,"('Start reduce')")
            if (myrank .eq. master) then
                write(stdout,"('Start reduce master')")
               call MPI_REDUCE(MPI_IN_PLACE, buffer(first:last), block_size, mpi_mod_ep, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            else
                write(stdout,"('Start reduce slave')")
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_REDUCE(buffer(first:last), block_size,   block_size, mpi_mod_ep, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            endif
            end_t = omp_get_wtime()
            write(stdout,'("Reduction took: ",f8.3," [s]")') end_t-start_t
         enddo !i

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && !defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         call naive_mpi_reduce_inplace_sum_ep(buffer,nelem,.false.)
#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         !Additionally we split the reduction into separate reductions each not larger than ~2GiB.
         call naive_mpi_reduce_inplace_sum_ep(buffer,nelem,.true.)
#endif

   end subroutine mpi_reduce_inplace_sum_ep

   subroutine mpi_reduce_sum_wp(src,dest,nelem)
      use omp_lib
      implicit none
      real(kind=wp), intent(in) :: src(:)
      real(kind=wp), intent(in) :: dest(:)
      integer, intent(in) :: nelem
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error, block_size
      integer :: n_times_reduce, i, first, last
      real(kind=wp) :: start_t, end_t, total_start, total_end

         write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem) then
            print *,nelem_conv,nelem
            call mpi_xermsg('mpi_mod', 'mpi_mod_reduce_inplace_sum_wp', &
                            'The nelem argument is too large for the mpi integer kind. The program will abort.', 1, 1)
         end if

         !MPI_REDUCE can handle at most max_data_count array elements to reduce
         !so we need to determine how many times we need to call it to reduce
         !the whole array which can have more than max_data_count array
         !elements. 
         n_times_reduce = ceiling(nelem_conv/(max_data_count_wp*1.0_wp))

         block_size = nelem_conv/n_times_reduce 

         !In each iteration we reduce the elements buffer(first:last).
         first = 0 
         last = 0
         do i=1,n_times_reduce
            if (i .eq. n_times_reduce) block_size = nelem_conv - last !reduce what's left 
            first = last + 1
            last = last + block_size
            if (n_times_reduce > 1) then
               write(stdout,'("Reducing elements: ",i0," to ",i0)') first, last
            endif
            start_t = omp_get_wtime()
            if (myrank .eq. master) then
               call MPI_REDUCE(src(first:last), dest(first:last), block_size, mpi_mod_wp, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            else
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_REDUCE(src(first:last), dest(first:last),   block_size, mpi_mod_wp, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            endif
            end_t = omp_get_wtime()
            write(stdout,'("Reduction took: ",f8.3," [s]")') end_t-start_t
         enddo !i

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && !defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
        ! call naive_mpi_reduce_inplace_sum_wp(buffer,nelem,.false.)
#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         !Additionally we split the reduction into separate reductions each not larger than ~2GiB.
         !call naive_mpi_reduce_inplace_sum_wp(buffer,nelem,.true.)
#endif

   end subroutine mpi_reduce_sum_wp

   subroutine mpi_reduce_sum_ep(src,dest,nelem)
      use omp_lib
      implicit none
      real(kind=ep), intent(in) :: src(:)
      real(kind=ep), intent(in) :: dest(:)
      integer, intent(in) :: nelem
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error, block_size
      integer :: n_times_reduce, i, first, last
      real(kind=ep) :: start_t, end_t, total_start, total_end

         write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem) then
            print *,nelem_conv,nelem
            call mpi_xermsg('mpi_mod', 'mpi_mod_reduce_inplace_sum_ep', &
                            'The nelem argument is too large for the mpi integer kind. The program will abort.', 1, 1)
         end if

         !MPI_REDUCE can handle at most max_data_count array elements to reduce
         !so we need to determine how many times we need to call it to reduce
         !the whole array which can have more than max_data_count array
         !elements. 
         n_times_reduce = ceiling(nelem_conv/(max_data_count_ep*1.0_ep))

         block_size = nelem_conv/n_times_reduce 

         !In each iteration we reduce the elements buffer(first:last).
         first = 0 
         last = 0
         do i=1,n_times_reduce
            if (i .eq. n_times_reduce) block_size = nelem_conv - last !reduce what's left 
            first = last + 1
            last = last + block_size
            if (n_times_reduce > 1) then
               write(stdout,'("Reducing elements: ",i0," to ",i0)') first, last
            endif
            start_t = omp_get_wtime()
            if (myrank .eq. master) then
               call MPI_REDUCE(src(first:last), dest(first:last), block_size, mpi_mod_ep, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            else
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_REDUCE(src(first:last), dest(first:last),   block_size, mpi_mod_ep, MPI_SUM, master, MPI_COMM_WORLD, ierr)
            endif
            end_t = omp_get_wtime()
            write(stdout,'("Reduction took: ",f8.3," [s]")') end_t-start_t
         enddo !i

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && !defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
        ! call naive_mpi_reduce_inplace_sum_ep(buffer,nelem,.false.)
#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         !Additionally we split the reduction into separate reductions each not larger than ~2GiB.
         !call naive_mpi_reduce_inplace_sum_ep(buffer,nelem,.true.)
#endif

   end subroutine mpi_reduce_sum_ep

        
   subroutine mpi_reduceall_sum_wp (src, dest, nelem, comm)
      use omp_lib
      real(kind=wp), intent(in) :: src(:)
      real(kind=wp), intent(in) :: dest(:)
      integer, intent(in) :: nelem
      integer(mpiint), intent(in), optional :: comm
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error, block_size
      integer :: n_times_reduce, i, first, last
      real(kind=wp) :: start_t, end_t, total_start, total_end

        if(nprocs <= 1) return
         !write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem) then
            print *,nelem_conv,nelem
            call mpi_xermsg('mpi_mod', 'mpi_mod_reduce_inplace_sum_wp', &
                            'The nelem argument is too large for the mpi integer kind. The program will abort.', 1, 1)
         end if

         !MPI_REDUCE can handle at most max_data_count array elements to reduce
         !so we need to determine how many times we need to call it to reduce
         !the whole array which can have more than max_data_count array
         !elements. 
         n_times_reduce = ceiling(nelem_conv/(max_data_count_wp*1.0_wp))

         block_size = nelem_conv/n_times_reduce 

         !In each iteration we reduce the elements buffer(first:last).
         first = 0 
         last = 0
         do i=1,n_times_reduce
            if (i .eq. n_times_reduce) block_size = nelem_conv - last !reduce what's left 
            first = last + 1
            last = last + block_size
            if (n_times_reduce > 1) then
               write(stdout,'("Reducing elements: ",i0," to ",i0)') first, last
            endif
            start_t = omp_get_wtime()
            if (myrank .eq. master) then
               call MPI_Allreduce(src(first:last), dest(first:last), block_size, mpi_mod_wp, MPI_SUM, mpi_mod_comm(comm), ierr)
            else
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_Allreduce(src(first:last), dest(first:last),   block_size, mpi_mod_wp, MPI_SUM, mpi_mod_comm(comm), ierr)
            endif
            end_t = omp_get_wtime()
           ! write(stdout,'("Reduction took: ",f8.3," [s]")') end_t-start_t
         enddo !i

         total_end = omp_get_wtime()
       !  write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

        ! write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && !defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
        ! call naive_mpi_reduce_inplace_sum_wp(buffer,nelem,.false.)
#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         !Additionally we split the reduction into separate reductions each not larger than ~2GiB.
         !call naive_mpi_reduce_inplace_sum_wp(buffer,nelem,.true.)
#endif

   end subroutine mpi_reduceall_sum_wp

   subroutine mpi_reduceall_sum_ep (src, dest, nelem, comm)
      use omp_lib
      real(kind=ep), intent(in) :: src(:)
      real(kind=ep), intent(in) :: dest(:)
      integer, intent(in) :: nelem
      integer(mpiint), intent(in), optional :: comm
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error, block_size
      integer :: n_times_reduce, i, first, last
      real(kind=ep) :: start_t, end_t, total_start, total_end

        if(nprocs <= 1) return
         !write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem) then
            print *,nelem_conv,nelem
            call mpi_xermsg('mpi_mod', 'mpi_mod_reduce_inplace_sum_ep', &
                            'The nelem argument is too large for the mpi integer kind. The program will abort.', 1, 1)
         end if

         !MPI_REDUCE can handle at most max_data_count array elements to reduce
         !so we need to determine how many times we need to call it to reduce
         !the whole array which can have more than max_data_count array
         !elements. 
         n_times_reduce = ceiling(nelem_conv/(max_data_count_ep*1.0_ep))

         block_size = nelem_conv/n_times_reduce 

         !In each iteration we reduce the elements buffer(first:last).
         first = 0 
         last = 0
         do i=1,n_times_reduce
            if (i .eq. n_times_reduce) block_size = nelem_conv - last !reduce what's left 
            first = last + 1
            last = last + block_size
            if (n_times_reduce > 1) then
               write(stdout,'("Reducing elements: ",i0," to ",i0)') first, last
            endif
            start_t = omp_get_wtime()
            if (myrank .eq. master) then
               call MPI_Allreduce(src(first:last), dest(first:last), block_size, mpi_mod_ep, MPI_SUM, mpi_mod_comm(comm), ierr)
            else
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_Allreduce(src(first:last), dest(first:last),   block_size, mpi_mod_ep, MPI_SUM, mpi_mod_comm(comm), ierr)
            endif
            end_t = omp_get_wtime()
           ! write(stdout,'("Reduction took: ",f8.3," [s]")') end_t-start_t
         enddo !i

         total_end = omp_get_wtime()
       !  write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

        ! write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && !defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
        ! call naive_mpi_reduce_inplace_sum_ep(buffer,nelem,.false.)
#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks) && defined(splitreduce)
         !We use the naive implementation of MPI_REDUCE in case of quad precision calculation and when the MPI library doesn't provide working MPI_REDUCE.
         !Additionally we split the reduction into separate reductions each not larger than ~2GiB.
         !call naive_mpi_reduce_inplace_sum_ep(buffer,nelem,.true.)
#endif
   end subroutine mpi_reduceall_sum_ep


   subroutine mpi_reduceall_max_default_integer(src,dest,communicator)
      use omp_lib
      implicit none
      integer, intent(in) :: src
      integer, intent(inout) :: dest
      integer(mpiint), intent(in), optional :: communicator
#ifdef usempi
!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error
      integer :: n_times_reduce, i, first, last, block_size
      real(kind=wp) :: start_t, end_t, total_start, total_end
      
      integer(kind=mpiint)      ::      mpi_comm
      integer(kind=mpiint)      ::      mpi_rank,mpi_nprocs
      
      if(present(communicator)) then
        mpi_comm = communicator
                 call MPI_COMM_SIZE(mpi_comm,mpi_nprocs,ierr)
                 call MPI_COMM_RANK( mpi_comm,mpi_rank,ierr)
      else
                mpi_comm = MPI_COMM_WORLD
                mpi_nprocs = nprocs
                mpi_rank = myrank
      
      endif

        if(kind(src) .eq. longint) then
                 call MPI_ALLREDUCE(src, dest,1_mpiint, MPI_INTEGER8, MPI_MAX,  mpi_comm, ierr)
        else
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_ALLREDUCE(src, dest,1_mpiint, MPI_INTEGER4, MPI_MAX,  mpi_comm, ierr)
        endif

        return
#else

        dest = src         
         
#endif
   end subroutine mpi_reduceall_max_default_integer


   subroutine mpi_reduceall_min_default_integer(src,dest,communicator)
      use omp_lib
      implicit none
      integer, intent(in) :: src
      integer, intent(inout) :: dest
      integer(mpiint), intent(in), optional :: communicator
#ifdef usempi

!We use MPI_REDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_REDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv, ierr, error
      integer :: n_times_reduce, i, first, last, block_size
      real(kind=wp) :: start_t, end_t, total_start, total_end
      integer(kind=mpiint)      ::      mpi_comm
      integer(kind=mpiint)      ::      mpi_rank,mpi_nprocs
      
      if(present(communicator)) then
        mpi_comm = communicator
                 call MPI_COMM_SIZE(mpi_comm,mpi_nprocs,ierr)
                 call MPI_COMM_RANK( mpi_comm,mpi_rank,ierr)
      else
                mpi_comm = MPI_COMM_WORLD
                mpi_nprocs = nprocs
                mpi_rank = myrank
      
      endif

        if(kind(src) .eq. longint) then
                 call MPI_ALLREDUCE(src, dest,1_mpiint, MPI_INTEGER8, MPI_MIN, mpi_comm, ierr)
        else if(kind(src) .eq. shortint) then
               !on other processes the reference to 'buffer' cannot be repeated (where the receive buffer is expected) when MPI_IN_PLACE is used so we put there a random variable instead (nelem_conv).
               call MPI_ALLREDUCE(src, dest,1_mpiint, MPI_INTEGER4, MPI_MIN,  mpi_comm, ierr)
        endif

        return
#else

        dest = src         
         
#endif

   end subroutine mpi_reduceall_min_default_integer

   subroutine naive_mpi_reduce_inplace_sum_wp(buffer,nelem,split)
      use omp_lib
      implicit none
      real(kind=wp), intent(inout) :: buffer(:)
      integer, intent(in) :: nelem
      logical, intent(in) :: split

      real(kind=wp), allocatable :: buffer_temp(:)
      integer, parameter :: tag = 1
      integer, parameter :: max_bytes = huge(dummy_32bit_integer)-2
      integer(kind=mpiint) :: conv
      integer :: j, k, err, chunk, chunk_size, no_chunks, k_min, k_max, current_nelem
      real(kind=wp) :: total_start, total_end

         write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

         call mpi_mod_barrier(err)

         total_start = omp_get_wtime()

         if (split) then
            if (nelem*wp_bytes > max_bytes) then
               write(stdout,'("The reduction of the array will be split since the array has size (bytes): ",&
                    &i0," which exceeds the limit (bytes): ",i0)') nelem*wp_bytes, max_bytes
               chunk_size = max_bytes/wp_bytes !number of elements that can be reduced in one go
               no_chunks = ceiling(1.0_wp*nelem/chunk_size)
               write(stdout,'("Number of chunks and the number of elements in one chunk: ",i0," ",i0)') no_chunks, chunk_size
            else
               chunk_size = nelem
               no_chunks = 1
            endif
         else
            chunk_size = nelem
            no_chunks = 1
         endif

         if (myrank .eq. master) then
            allocate(buffer_temp(chunk_size),stat=err)
            if (err .ne. 0) call mpi_xermsg('mpi_mod','mpi_mod_reduce_inplace_sum_wp','Memory allocation on master failed.',err,1)
         endif
 
         k_max = 0
         do chunk=1,no_chunks

            k_min = k_max+1
            if (chunk .eq. no_chunks) then
               k_max = nelem
            else
               k_max = k_max + chunk_size
            endif

            current_nelem = k_max-k_min+1

            if (split) write(stdout,'("Reducing chunk number: ",i5,";array element range: ",i0," ",i0)') chunk, k_min,k_max

            !This can be rewritten so that the reduction is on the whole 2D array but for most applications this should be good enough. 
            if (myrank .eq. master) then
               !receive from all and sum
               do j=1,nprocs-1
                  conv = j
                  call mpi_mod_recv(conv,tag,buffer_temp,current_nelem)
                  !Note that there will be some false sharing in the loop below but it should be insignificant in cases nelem >> nprocs.
                  !$OMP PARALLEL DEFAULT(NONE) SHARED(buffer_temp,nelem,buffer,k_min,k_max) PRIVATE(k)
                  !$OMP DO
                  do k=k_min,k_max
                     buffer(k) = buffer(k) + buffer_temp(k-k_min+1)
                  enddo !k
                  !$OMP END DO
                  !$OMP END PARALLEL
               enddo !j
            else
               !send my array to the master
               call mpi_mod_send(master,buffer(k_min:k_max),tag,current_nelem)
            endif
   
            call mpi_mod_barrier(err)

         enddo !chunk

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_wp")')

   end subroutine naive_mpi_reduce_inplace_sum_wp

   subroutine naive_mpi_reduce_inplace_sum_ep(buffer,nelem,split)
      use omp_lib
      implicit none
      real(kind=ep), intent(inout) :: buffer(:)
      integer, intent(in) :: nelem
      logical, intent(in) :: split

      real(kind=ep), allocatable :: buffer_temp(:)
      integer, parameter :: tag = 1
      integer, parameter :: max_bytes = huge(dummy_32bit_integer)-2
      integer(kind=mpiint) :: conv
      integer :: j, k, err, chunk, chunk_size, no_chunks, k_min, k_max, current_nelem
      real(kind=ep) :: total_start, total_end

         write(stdout,'("--------->","mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

         call mpi_mod_barrier(err)

         total_start = omp_get_wtime()

         if (split) then
            if (nelem*ep1_bytes > max_bytes) then
               write(stdout,'("The reduction of the array will be split since the array has size (bytes): ",i0,&
                    &" which exceeds the limit (bytes): ",i0)') nelem*ep1_bytes, max_bytes
               chunk_size = max_bytes/ep1_bytes !number of elements that can be reduced in one go
               no_chunks = ceiling(1.0_ep*nelem/chunk_size)
               write(stdout,'("Number of chunks and the number of elements in one chunk: ",i0," ",i0)') no_chunks, chunk_size
            else
               chunk_size = nelem
               no_chunks = 1
            endif
         else
            chunk_size = nelem
            no_chunks = 1
         endif

         if (myrank .eq. master) then
            allocate(buffer_temp(chunk_size),stat=err)
            if (err .ne. 0) call mpi_xermsg('mpi_mod','mpi_mod_reduce_inplace_sum_ep','Memory allocation on master failed.',err,1)
         endif
 
         k_max = 0
         do chunk=1,no_chunks

            k_min = k_max+1
            if (chunk .eq. no_chunks) then
               k_max = nelem
            else
               k_max = k_max + chunk_size
            endif

            current_nelem = k_max-k_min+1

            if (split) write(stdout,'("Reducing chunk number: ",i5,";array element range: ",i0," ",i0)') chunk, k_min,k_max

            !This can be rewritten so that the reduction is on the whole 2D array but for most applications this should be good enough. 
            if (myrank .eq. master) then
               !receive from all and sum
               do j=1,nprocs-1
                  conv = j
                  call mpi_mod_recv(conv,tag,buffer_temp,current_nelem)
                  !Note that there will be some false sharing in the loop below but it should be insignificant in cases nelem >> nprocs.
                  !$OMP PARALLEL DEFAULT(NONE) SHARED(buffer_temp,nelem,buffer,k_min,k_max) PRIVATE(k)
                  !$OMP DO
                  do k=k_min,k_max
                     buffer(k) = buffer(k) + buffer_temp(k-k_min+1)
                  enddo !k
                  !$OMP END DO
                  !$OMP END PARALLEL
               enddo !j
            else
               !send my array to the master
               call mpi_mod_send(master,buffer(k_min:k_max),tag,current_nelem)
            endif
   
            call mpi_mod_barrier(err)

         enddo !chunk

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_mod_reduce_inplace_sum_ep")')

   end subroutine naive_mpi_reduce_inplace_sum_ep

   subroutine mpi_reduceall_inplace_sum_wp(buffer, nelem, comm)
      use omp_lib
      real(kind=wp), intent(inout) :: buffer(:)
      integer, intent(in) :: nelem
      integer(mpiint), intent(in), optional :: comm
      integer(kind=mpiint) :: ierr, error
#ifdef usempi
      integer(kind=mpiint) :: nelem_conv
      logical :: split = .false.
      integer :: block, no_blocks, first, last, block_size
      real(kind=wp) :: total_start, total_end

         write(stdout,'("--------->","mpi_mod:mpi_reduceall_inplace_sum_wp")')
         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem .or. nelem_conv > max_data_count_wp) then
            print *,nelem_conv,nelem,max_data_count_wp
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_wp', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
         endif

         if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_wp
            if (no_blocks*max_data_count_wp < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_wp-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_Allreduce(MPI_IN_PLACE, buffer(first:last), nelem_conv, mpi_mod_wp, MPI_SUM, mpi_mod_comm(comm), ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            call MPI_Allreduce(MPI_IN_PLACE, buffer, nelem_conv, mpi_mod_wp, MPI_SUM, mpi_mod_comm(comm), ierr)
         endif

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_reduceall_inplace_sum_wp")')
#endif
   end subroutine mpi_reduceall_inplace_sum_wp

   subroutine mpi_reduceall_inplace_sum_ep (buffer, nelem, comm)
      use omp_lib
      real(kind=ep), intent(inout) :: buffer(:)
      integer, intent(in) :: nelem
      integer(mpiint), intent(in), optional :: comm
      integer(kind=mpiint) :: ierr, error
#if (defined(usempi) && !defined(usequadprec)) || (defined(usempi) && defined(usequadprec) && defined(quadreduceworks))
!We use MPI_ALLREDUCE in these cases:
!A: double precision
!B: quad precision but only if MPI_ALLREDUCE is guaranteed to work.
      integer(kind=mpiint) :: nelem_conv
      logical :: split = .false.
      integer :: block, no_blocks, first, last, block_size
      real(kind=ep) :: total_start, total_end

         write(stdout,'("--------->","mpi_mod:mpi_reduceall_inplace_sum_ep")')
         total_start = omp_get_wtime()

         nelem_conv = nelem
         if (nelem_conv .ne. nelem .or. nelem_conv > max_data_count_ep) then
            print *,nelem_conv,nelem,max_data_count_ep
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_ep', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
         endif

         if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_ep
            if (no_blocks*max_data_count_ep < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_ep-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_Allreduce(MPI_IN_PLACE, buffer(first:last), nelem_conv, mpi_mod_ep, MPI_SUM, mpi_mod_comm(comm), ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            call MPI_Allreduce(MPI_IN_PLACE, buffer, nelem_conv, mpi_mod_ep, MPI_SUM, mpi_mod_comm(comm), ierr)
         endif

         total_end = omp_get_wtime()
         write(stdout,'("Complete reduction took: ",f8.3," [s]")') total_end-total_start

         write(stdout,'("<---------","done:mpi_mod:mpi_reduceall_inplace_sum_ep")')

#elif defined(usempi) && defined(usequadprec) && !defined(quadreduceworks)
!We don't have any working implementation of MPI_ALLREDUCE
         call mpi_xermsg('mpi_mod', 'mpi_mod_reduceall_inplace_sum_ep', &
                         'Quad precision version of MPI_ALLREDUCE would not work (FPP directive usequadprec not defined). ', 2, 1)
#endif

   end subroutine mpi_reduceall_inplace_sum_ep


   !>This will rotate a combination of the number of elements, an integer array and float array once around in a ring formation
   subroutine mpi_mod_rotate_arrays_around_ring_wp (elem_count, int_array, wp_array, max_num_elements, comm_opt)

      integer(longint), intent(in)           :: max_num_elements
      integer(longint), intent(inout)        :: elem_count, int_array(:)
      real(wp),         intent(inout)        :: wp_array(:)
      integer(mpiint),  intent(in), optional :: comm_opt
#ifdef usempi
      integer(kind=mpiint) :: comm, comm_size, comm_rank, proc_left, proc_right, last_proc, ierr
      integer(kind=mpiint) :: stat(MPI_STATUS_SIZE)
      integer(kind=mpiint) :: nelem_conv

      logical :: split = .false.
      integer :: block, no_blocks, first, last, block_size,nelem

      if (present(comm_opt)) then
         comm = comm_opt
         call MPI_Comm_size(comm, comm_size, ierr)
         call MPI_Comm_rank(comm, comm_rank, ierr)
      else
         comm = MPI_COMM_WORLD
         comm_size = nprocs
         comm_rank = myrank
      end if

      if (comm_size == 1) return

      proc_right = mod((comm_rank + 1_mpiint), comm_size)
      proc_left  = mod((comm_rank - 1_mpiint + comm_size), comm_size)

      call MPI_Sendrecv_replace(elem_count, 1_mpiint, MPI_INTEGER8, proc_right, 1_mpiint, proc_left, 1_mpiint, comm, stat, ierr)

     nelem = max_num_elements
     !Do wp array first
     nelem_conv = nelem
     if (nelem_conv .ne. max_num_elements .or. nelem_conv > max_data_count_wp) then
            print *,nelem_conv,nelem,max_data_count_wp
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_wp', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif
     
     
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_wp
            if (no_blocks*max_data_count_wp < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_wp-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (wp_array(first:last), nelem_conv, mpi_mod_wp, proc_right, &
                                          3_mpiint, proc_left, 3_mpiint, comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = max_num_elements
            call MPI_sendrecv_replace (wp_array, nelem_conv, mpi_mod_wp, proc_right, 3_mpiint, proc_left, &
                                       3_mpiint, comm, stat, ierr)
     endif
     
     split = .false.
     
     nelem = max_num_elements*2
     !Now int array
       !Do wp array first
     nelem_conv = nelem
     if (nelem_conv .ne.  nelem .or. nelem_conv > max_data_count_wp) then
            print *,nelem_conv,nelem,max_data_count_wp
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_wp', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif   
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_wp
            if (no_blocks*max_data_count_wp < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_wp-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (int_array(first:last), nelem_conv, MPI_INTEGER8, proc_right, 2_mpiint, proc_left, &
                                          2_mpiint, comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = max_num_elements * 2
            call MPI_sendrecv_replace (int_array, nelem_conv, MPI_INTEGER8, proc_right, 2_mpiint, proc_left, 2_mpiint, &
                                       comm, stat, ierr)
     endif     
                
                
      !Everyone starts reciving
      !if(myrank /= master) then
        !       call MPI_Recv(&token, 1, MPI_INT, world_rank - 1, 0,
         !                      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
          
         
         
#endif

   end subroutine mpi_mod_rotate_arrays_around_ring_wp

   !>This will rotate a combination of the number of elements, an integer array and float array once around in a ring formation
   subroutine mpi_mod_rotate_arrays_around_ring_ep (elem_count, int_array, ep_array, max_num_elements, comm_opt)

      integer(longint), intent(in)           :: max_num_elements
      integer(longint), intent(inout)        :: elem_count,int_array(:)
      real(ep),         intent(inout)        :: ep_array(:)
      integer(mpiint),  intent(in), optional :: comm_opt
#ifdef usempi
      integer(kind=mpiint) :: comm, comm_rank, comm_size, proc_left, proc_right, last_proc, ierr
      integer(kind=mpiint) :: stat(MPI_STATUS_SIZE)
      integer(kind=mpiint) :: nelem_conv

      logical :: split = .false.
      integer :: block, no_blocks, first, last, block_size,nelem

      if (present(comm_opt)) then
         comm = comm_opt
         call MPI_Comm_size(comm, comm_size, ierr)
         call MPI_Comm_rank(comm, comm_rank, ierr)
      else
         comm = MPI_COMM_WORLD
         comm_size = nprocs
         comm_rank = myrank
      end if

      if (comm_size == 1) return

      proc_right = mod((comm_rank + 1_mpiint), comm_size)
      proc_left  = mod((comm_rank - 1_mpiint + comm_size), comm_size)

      call MPI_Sendrecv_replace(elem_count, 1_mpiint, MPI_INTEGER8, proc_right, 1_mpiint, proc_left, 1_mpiint, comm, stat, ierr)

     nelem = max_num_elements
     !Do ep array first
     nelem_conv = nelem
     if (nelem_conv .ne. max_num_elements .or. nelem_conv > max_data_count_ep) then
            print *,nelem_conv,nelem,max_data_count_ep
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_ep', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif
     
     
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_ep
            if (no_blocks*max_data_count_ep < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_ep-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (ep_array(first:last), nelem_conv, mpi_mod_ep, proc_right, 3_mpiint, &
                                          proc_left, 3_mpiint, comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = max_num_elements
            call MPI_sendrecv_replace (ep_array, nelem_conv, mpi_mod_ep, proc_right, 3_mpiint, proc_left, &
                                       3_mpiint, comm, stat, ierr)
     endif
     
     split = .false.
     
     nelem = max_num_elements*2
     !Now int array
       !Do ep array first
     nelem_conv = nelem
     if (nelem_conv .ne.  nelem .or. nelem_conv > max_data_count_ep) then
            print *,nelem_conv,nelem,max_data_count_ep
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_ep', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif   
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_ep
            if (no_blocks*max_data_count_ep < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_ep-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (int_array(first:last), nelem_conv, MPI_INTEGER8, proc_right, 2_mpiint, proc_left,&
                                          2_mpiint, comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = max_num_elements * 2
            call MPI_sendrecv_replace (int_array, nelem_conv, MPI_INTEGER8, proc_right, 2_mpiint, proc_left, &
                                       2_mpiint, comm, stat, ierr)
     endif     
                
                
      !Everyone starts reciving
      !if(myrank /= master) then
        !       call MPI_Recv(&token, 1, MPI_INT, world_rank - 1, 0,
         !                      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
          
         
         
#endif

   end subroutine mpi_mod_rotate_arrays_around_ring_ep

   !>This will rotate a combination of the number of elements, an integer array and float array once around in a ring formation
   subroutine mpi_mod_rotate_wp_arrays_around_ring(elem_count,wp_array,max_num_elements,communicator)
      implicit none
      integer(longint),intent(in)        :: max_num_elements
      integer(longint), intent(inout) :: elem_count
      real(wp),intent(inout) :: wp_array(:)
      integer(mpiint), optional, intent(in) :: communicator
#ifdef usempi
      integer(kind=mpiint)      ::      proc_left,proc_right,last_proc,ierr
      integer(kind=mpiint)     ::      stat(MPI_STATUS_SIZE)
      integer(kind=mpiint) :: nelem_conv
      logical :: split = .false.
      integer(kind=mpiint)      ::      comm_rank,comm_nprocs,mpi_comm
      integer :: block, no_blocks, first, last, block_size,nelem
      
      
      if (nprocs .eq. 1) return
    
      if(present(communicator)) then
                 mpi_comm = communicator
                 call MPI_COMM_SIZE(mpi_comm,comm_nprocs,ierr)
                 call MPI_COMM_RANK( mpi_comm,comm_rank,ierr)
     else
        mpi_comm = MPI_COMM_WORLD
        comm_rank = myrank
        comm_nprocs = nprocs
     endif
      
      proc_right = mod((comm_rank + 1_mpiint), comm_nprocs)
      proc_left  = mod((comm_rank - 1_mpiint + comm_nprocs), comm_nprocs)


     call MPI_sendrecv_replace(elem_count,1_mpiint,MPI_INTEGER8,proc_right,1_mpiint,proc_left,1_mpiint,mpi_comm ,stat,ierr)
     
     nelem = max_num_elements
     !Do wp array first
     nelem_conv = nelem
     if (nelem_conv .ne. max_num_elements .or. nelem_conv > max_data_count_wp) then
            print *,nelem_conv,nelem,max_data_count_wp
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_wp', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif
     
     
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_wp
            if (no_blocks*max_data_count_wp < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_wp-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (wp_array(first:last), nelem_conv, mpi_mod_wp, proc_right, 3_mpiint, &
                                          proc_left, 3_mpiint, mpi_comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = nelem
            call MPI_sendrecv_replace(wp_array, nelem_conv,mpi_mod_wp,proc_right,3_mpiint,proc_left,3_mpiint,mpi_comm ,stat,ierr)
     endif
     
                
#endif

   end subroutine mpi_mod_rotate_wp_arrays_around_ring

   !>This will rotate a combination of the number of elements, an integer array and float array once around in a ring formation
   subroutine mpi_mod_rotate_ep_arrays_around_ring(elem_count,ep_array,max_num_elements,communicator)
      implicit none
      integer(longint),intent(in)        :: max_num_elements
      integer(longint), intent(inout) :: elem_count
      real(ep),intent(inout) :: ep_array(:)
      integer(mpiint), optional, intent(in) :: communicator
#ifdef usempi
      integer(kind=mpiint)      ::      proc_left,proc_right,last_proc,ierr
      integer(kind=mpiint)     ::      stat(MPI_STATUS_SIZE)
      integer(kind=mpiint) :: nelem_conv
      logical :: split = .false.
      integer(kind=mpiint)      ::      comm_rank,comm_nprocs,mpi_comm
      integer :: block, no_blocks, first, last, block_size,nelem
      
      
      if (nprocs .eq. 1) return
    
      if(present(communicator)) then
                 mpi_comm = communicator
                 call MPI_COMM_SIZE(mpi_comm,comm_nprocs,ierr)
                 call MPI_COMM_RANK( mpi_comm,comm_rank,ierr)
     else
        mpi_comm = MPI_COMM_WORLD
        comm_rank = myrank
        comm_nprocs = nprocs
     endif
      
      proc_right = mod((comm_rank + 1_mpiint),comm_nprocs)
      proc_left  = mod((comm_rank - 1_mpiint + comm_nprocs),comm_nprocs)


     call MPI_sendrecv_replace(elem_count,1_mpiint,MPI_INTEGER8,proc_right,1_mpiint,proc_left,1_mpiint,mpi_comm ,stat,ierr)
     
     nelem = max_num_elements
     !Do ep array first
     nelem_conv = nelem
     if (nelem_conv .ne. max_num_elements .or. nelem_conv > max_data_count_ep) then
            print *,nelem_conv,nelem,max_data_count_ep
            call xermsg ('mpi_mod', 'mpi_mod_reduceall_inplace_sum_ep', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif
     
     
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_ep
            if (no_blocks*max_data_count_ep < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_ep-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (ep_array(first:last), nelem_conv, mpi_mod_ep, proc_right, 3_mpiint, &
                                          proc_left, 3_mpiint, mpi_comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = nelem
            call MPI_sendrecv_replace(ep_array, nelem_conv,mpi_mod_ep,proc_right,3_mpiint,proc_left,3_mpiint,mpi_comm ,stat,ierr)
     endif
     
                
#endif

   end subroutine mpi_mod_rotate_ep_arrays_around_ring



   !>This will rotate a combination of the number of elements, an integer array and float array once around in a ring formation
   subroutine mpi_mod_rotate_int_arrays_around_ring(elem_count,int_array,max_num_elements,communicator)
      implicit none
      integer(longint),intent(in)        :: max_num_elements
      integer(longint), intent(inout) :: elem_count
      integer(longint),intent(inout) :: int_array(:)
      integer(mpiint), optional, intent(in) :: communicator
#ifdef usempi
      integer(kind=mpiint)      ::      proc_left,proc_right,last_proc,ierr
      integer(kind=mpiint)     ::      stat(MPI_STATUS_SIZE)
      integer(kind=mpiint) :: nelem_conv
      logical :: split = .false.
      integer(kind=mpiint)      ::      comm_rank,comm_nprocs,mpi_comm
      integer :: block, no_blocks, first, last, block_size,nelem
      
      
      if (nprocs .eq. 1) return
    
      if(present(communicator)) then
                 mpi_comm = communicator
                 call MPI_COMM_SIZE(mpi_comm,comm_nprocs,ierr)
                 call MPI_COMM_RANK( mpi_comm,comm_rank,ierr)
     else
        mpi_comm = MPI_COMM_WORLD
        comm_rank = myrank
        comm_nprocs = nprocs
     endif
      
      proc_right = mod((comm_rank + 1_mpiint),comm_nprocs)
      proc_left  = mod((comm_rank - 1_mpiint + comm_nprocs),comm_nprocs)


     call MPI_sendrecv_replace(elem_count,1_mpiint,MPI_INTEGER8,proc_right,1_mpiint,proc_left,1_mpiint,mpi_comm ,stat,ierr)
     
     nelem = max_num_elements
     !Do int array first
     nelem_conv = nelem
     if (nelem_conv .ne. max_num_elements .or. nelem_conv > max_data_count_longint) then
            print *,nelem_conv,nelem,max_data_count_longint
            call xermsg ('mpi_mod', 'mpi_mod_rotate_int_arrays_around_ring', &
                         'The nelem argument is too large for the mpi integer kind. &
                         &The reduction will be split into several bits.', 1, 0)
            split = .true.
     endif
     
     
     if (split) then !split the reduction into bits of length at most max_data_count
            no_blocks = nelem/max_data_count_longint
            if (no_blocks*max_data_count_longint < nelem) no_blocks = no_blocks+1
            first=1
            do block=1,no_blocks
               last = min(first+max_data_count_longint-1,nelem)
               nelem_conv=last-first+1
               write(stdout,'("Reducing block no., index range, number of elements:",4(i0,1x))') block,first,last,nelem_conv
               call MPI_sendrecv_replace (int_array(first:last), nelem_conv, MPI_INTEGER8, proc_right, 3_mpiint, &
                                          proc_left, 3_mpiint, mpi_comm, stat, ierr)
               first=last+1
            enddo !block             
         else !reduce all elements at once
            nelem_conv = nelem
            call MPI_sendrecv_replace (int_array, nelem_conv, MPI_INTEGER8, proc_right, 3_mpiint, &
                                       proc_left, 3_mpiint, mpi_comm, stat, ierr)
     endif
     
                
#endif

   end subroutine mpi_mod_rotate_int_arrays_around_ring
   

   !> All gather for one integer send from every process.
   subroutine mpi_mod_allgather_int32 (send, receive, comm)

      integer(int32),  intent(in)           :: send
      integer(int32),  intent(out)          :: receive(:)
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr, one = 1

      call MPI_Allgather(send, one, MPI_INTEGER4, receive, one, MPI_INTEGER4, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_allgather_int32', 'All-gather failed.', int(ierr), 1)
      end if
#else
      receive(1) = send
#endif
   end subroutine mpi_mod_allgather_int32


   subroutine mpi_mod_allgather_int64 (send, receive, comm)

      integer(int64),  intent(in)           :: send
      integer(int64),  intent(out)          :: receive(:)
      integer(mpiint), intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: ierr, one = 1

      call MPI_Allgather(send, one, MPI_INTEGER8, receive, one, MPI_INTEGER8, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_allgather_int64', 'All-gather failed.', int(ierr), 1)
      end if
#else
      receive(1) = send
#endif
   end subroutine mpi_mod_allgather_int64


   !> All gather for one string of length MPI_MAX_PROCESSOR_NAME send from every process.
   subroutine mpi_mod_allgather_character (send, receive, comm)
#ifdef usempi
      character(len=MPI_MAX_PROCESSOR_NAME), intent(in) :: send
      character(len=MPI_MAX_PROCESSOR_NAME), intent(out) :: receive(:)
      integer(mpiint), intent(in), optional :: comm
      integer(mpiint) :: ierr, n

      n = MPI_MAX_PROCESSOR_NAME
      call MPI_Allgather(send, n, MPI_CHARACTER, receive, n, MPI_CHARACTER, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_allgather_int64', 'All-gather failed.', int(ierr), 1)
      end if
#else
      character(len=*), intent(in) :: send
      character(len=*), intent(out) :: receive(:)
      integer(mpiint),  intent(in), optional :: comm
      receive(1) = send
#endif
   end subroutine mpi_mod_allgather_character


   !> Interface to the routine MPI_BARRIER.
  function mpi_mod_wtime() result(t)
      implicit none
      real(wp) :: t

      integer         :: count, count_rate, count_max
      real(wp), save :: overflow   =  0
      integer, save   :: last_count = -1
#ifdef usempi
        t = MPI_WTIME()
        
#else
        call system_clock(count,count_rate,count_max)
        ! Try to detect a rollover
        if (count < last_count) then
           overflow = overflow + count_max
        end if
        last_count = count
        ! Convert to seconds
        t = (overflow+count)/count_rate

#endif

   end function mpi_mod_wtime


   subroutine mpi_mod_file_open_write (filename, fh, ierr, comm)

      character(len=*), intent(in)           :: filename
      integer(mpiint),  intent(out)          :: fh
      integer,          intent(out)          :: ierr
      integer(mpiint),  intent(in), optional :: comm
#ifdef usempi
      integer(mpiint) :: mpierr

      call MPI_File_open(mpi_mod_comm(comm), filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpierr)
      ierr = mpierr
#else
      open (newunit = fh, file = filename, iostat = ierr, access = 'stream')
#endif
   end subroutine mpi_mod_file_open_write


   !> \brief   Let master process write one 4-byte integer to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Convenience wrapper around \ref mpi_mod_file_write_array1d_int32, which does the real work.
   !>
   !> \param[in] fh  File handle as returned by mpi_mod_file_open_*.
   !> \param[in] n   Integer to write.
   !>
   subroutine mpi_mod_file_write_int32 (fh, n)

      integer(mpiint), intent(in) :: fh
      integer(int32),  intent(in) :: n

      integer(int32) :: buf(1)

      ! convert to rank-1 array and delegate work
      buf(1) = n
      call mpi_mod_file_write_array1d_int32(fh, buf, 1)

   end subroutine mpi_mod_file_write_int32


   !> \brief   Let master process write array of 4-byte integers to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Uses blocking collective MPI output routine to write array of integers. Only master process
   !> will do the writing, other processes behave as if the provided length parameter was zero.
   !> The integers are written at the master's file position, which gets updated. Positions of other
   !> processes in the file remain as they were.
   !>
   !> \param[in] fh      File handle as returned by mpi_mod_file_open_*.
   !> \param[in] array   Integer array to write.
   !> \param[in] length  Number of elements in the array.
   !>
   subroutine mpi_mod_file_write_array1d_int32 (fh, array, length)

      integer(mpiint), intent(in) :: fh
      integer,         intent(in) :: length
      integer(int32),  intent(in) :: array(length)
#ifdef usempi
      integer(mpiint) :: mpilen, mpierr, stat(MPI_STATUS_SIZE)
      integer(mpiofs) :: pos, offs

      mpilen = merge(length, 0, myrank == master)  ! only master will write something

      call MPI_File_get_position(fh, pos, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_int32', &
                         'Failed to obtain file position.', int(mpierr), 1)
      end if

      call MPI_File_get_byte_offset(fh, pos, offs, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_int32', &
                         'Failed to convert file position to byte offset.', int(mpierr), 1)
      end if

      call MPI_File_set_view(fh, offs, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_int32', &
                         'Failed to define file view.', int(mpierr), 1)
      end if

      call MPI_File_write_all(fh, array, mpilen, MPI_INTEGER4, stat, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_int32', &
                         'Failed to write file.', int(mpierr), 1)
      end if
#else
      write (fh) array(1:length)
#endif
   end subroutine mpi_mod_file_write_array1d_int32


   !> \brief   Let master process write 2d array of 4-byte integers to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Convenience wrapper around \ref mpi_mod_file_write_array1d_int32, which does the real work.
   !>
   !> \param[in] fh       File handle as returned by mpi_mod_file_open_*.
   !> \param[in] array    Two-dimensional integer array to write.
   !> \param[in] length1  Number of rows (= leading dimension).
   !> \param[in] length2  Number of columns.
   !>
   subroutine mpi_mod_file_write_array2d_int32 (fh, array, length1, length2)

      integer(mpiint), intent(in) :: fh
      integer,         intent(in) :: length1, length2
      integer(int32),  intent(in) :: array(length1, length2)

      ! convert to rank-1 array and delegate work
      call mpi_mod_file_write_array1d_int32(fh, array, length1 * length2)

   end subroutine mpi_mod_file_write_array2d_int32


   !> \brief   Let master process write an 8-byte real to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Convenience wrapper around \ref mpi_mod_file_write_array1d_real64, which does the real work.
   !>
   !> \param[in] fh  File handle as returned by mpi_mod_file_open_*.
   !> \param[in] x   Real number to write.
   !>
   subroutine mpi_mod_file_write_real64 (fh, x)

      integer(mpiint), intent(in) :: fh
      real(real64),    intent(in) :: x

      real(real64) :: buf(1)

      ! convert to rank-1 array and delegate work
      buf(1) = x
      call mpi_mod_file_write_array1d_real64(fh, buf, 1)

   end subroutine mpi_mod_file_write_real64


   !> \brief   Let master process write array of 8-byte reals to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Uses blocking collective MPI output routine to write array of reals. Only master process
   !> will do the writing, other processes behave as if the provided length parameter was zero.
   !> The reals are written at the master's file position, which gets updated. Positions of other
   !> processes in the file remain as they were.
   !>
   !> \param[in] fh      File handle as returned by mpi_mod_file_open_*.
   !> \param[in] array   Real array to write.
   !> \param[in] length  Number of elements in the array.
   !>
   subroutine mpi_mod_file_write_array1d_real64 (fh, array, length)

      integer(mpiint), intent(in) :: fh
      integer,         intent(in) :: length
      real(real64),    intent(in) :: array(length)
#ifdef usempi
      integer(mpiint) :: mpifh, mpilen, mpierr, stat(MPI_STATUS_SIZE)
      integer(mpiofs) :: pos, offs

      mpifh = fh
      mpilen = merge(length, 0, myrank == master)  ! only master will write something

      call MPI_File_get_position(fh, pos, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_real64', &
                         'Failed to obtain file position.', int(mpierr), 1)
      end if

      call MPI_File_get_byte_offset(fh, pos, offs, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_real64', &
                         'Failed to convert file position to byte offset.', int(mpierr), 1)
      end if

      call MPI_File_set_view(mpifh, offs, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_real64', &
                         'Failed to define file view.', int(mpierr), 1)
      end if

      call MPI_File_write_all(mpifh, array, mpilen, MPI_REAL8, stat, mpierr)
      if (mpierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_array1d_real64', &
                         'Failed to write file.', int(mpierr), 1)
      end if
#else
      write (fh) array(1:length)
#endif
   end subroutine mpi_mod_file_write_array1d_real64


   !> \brief   Let master process write 2d array of 8-byte reals to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Convenience wrapper around \ref mpi_mod_file_write_array1d_real64, which does the real work.
   !>
   !> \param[in] fh       File handle as returned by mpi_mod_file_open_*.
   !> \param[in] array    Two-dimensional real array to write.
   !> \param[in] length1  Number of rows (= leading dimension).
   !> \param[in] length2  Number of columns.
   !>
   subroutine mpi_mod_file_write_array2d_real64 (fh, array, length1, length2)

      integer(mpiint), intent(in) :: fh
      integer,         intent(in) :: length1, length2
      real(real64),    intent(in) :: array(length1, length2)

      ! convert to rank-1 array and delegate work
      call mpi_mod_file_write_array1d_real64(fh, array, length1 * length2)

   end subroutine mpi_mod_file_write_array2d_real64


   !> \brief   Let master process write 3d array of 8-byte reals to stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Convenience wrapper around \ref mpi_mod_file_write_array1d_real64, which does the real work.
   !>
   !> \param[in] fh       File handle as returned by mpi_mod_file_open_*.
   !> \param[in] array    Three-dimensional real array to write.
   !> \param[in] length1  Leading dimension.
   !> \param[in] length2  Next dimension.
   !> \param[in] length3  Next dimension.
   !>
   subroutine mpi_mod_file_write_array3d_real64 (fh, array, length1, length2, length3)

      integer(mpiint), intent(in) :: fh
      integer,         intent(in) :: length1, length2, length3
      real(real64),    intent(in) :: array(length1, length2, length3)

      ! convert to rank-1 array and delegate work
      call mpi_mod_file_write_array1d_real64(fh, array, length1 * length2 * length3)

   end subroutine mpi_mod_file_write_array3d_real64


   !> \brief   Collectively resize file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Change size of a file; mostly used to truncate the file to zero size. Only available when compiled with MPI.
   !>
   !> \param[in] fh  File handle as returned by mpi_mod_file_open_*.
   !> \param[in] sz  New size of the file.
   !>
   subroutine mpi_mod_file_set_size (fh, sz)

      integer(mpiint), intent(in) :: fh
      integer,         intent(in) :: sz
#ifdef usempi
      integer(mpiint) :: ierr

      call MPI_File_set_size(fh, int(sz, mpiofs), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_set_size', 'Failed to set file size.', int(ierr), 1)
      end if
#else
      ! not implemented
#endif
   end subroutine mpi_mod_file_set_size


   !> \brief   Write block-cyclically distributed matrix to a stream file
   !> \authors J Benda
   !> \date    2019
   !>
   !> Assumes ScaLAPACK-compatible row-major (!) block-cyclic distribution. Writes the data at the master process
   !> file position, which it first broadcasts to other processes.
   !>
   !> \warning At the moment, there is a limitation on the element count of the local portion of the distributed matrix, which
   !>          mustn't exceed 2**31 elements (i.e. 16 GiB).
   !>
   !> \param[in] fh        File handle as returned by mpi_mod_file_open_*.
   !> \param[in] m         Number of rows in the distributed matrix.
   !> \param[in] n         Number of columns in the distributed matrix.
   !> \param[in] nprow     Number of rows in the process grid.
   !> \param[in] npcol     Number of columns in the process grid.
   !> \param[in] rb        Block row count.
   !> \param[in] cb        Block column column.
   !> \param[in] locA      Local portion of the distributed matrix.
   !> \param[in] myrows    Number of rows in locA.
   !> \param[in] mycols    Number of columns in locA.
   !> \param[in] comm_opt  Optional communicator on which the file is shared (MPI world used if not provided).
   !>
   subroutine mpi_mod_file_write_darray2d_real64 (fh, m, n, nprow, npcol, rb, cb, locA, myrows, mycols, comm)

      integer(mpiint), intent(in)           :: fh
      integer,         intent(in)           :: m, n, nprow, npcol, rb, cb, myrows, mycols
      integer(mpiint), intent(in), optional :: comm
      real(real64),    intent(in)           :: locA(myrows, mycols)
#ifdef usempi
      integer(mpiint) :: datype, ierr, comm_size, comm_rank, stat(MPI_STATUS_SIZE), ndims, &
                         gsizes(2), distrb(2), dargs(2), psizes(2), one = 1, locsiz
      integer(mpiofs) :: pos, offs

      offs   = 0
      ndims  = 2
      gsizes = (/ m, n /)
      distrb = (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC /)
      dargs  = (/ rb, cb /)
      psizes = (/ nprow, npcol /)
      locsiz = myrows * mycols

      ! MPI environment info
      call MPI_Comm_size(mpi_mod_comm(comm), comm_size, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to get comunicator size.', int(ierr), 1)
      end if
      call MPI_Comm_rank(mpi_mod_comm(comm), comm_rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to get communicator rank.', int(ierr), 1)
      end if

      ! impose master's file position on other processes
      call MPI_File_sync(fh, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to synchronize file.', int(ierr), 1)
      end if
      call MPI_File_get_position(fh, pos, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to get file position.', int(ierr), 1)
      end if
      call MPI_File_get_byte_offset(fh, pos, offs, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to convert file position to byte offset.', int(ierr), 1)
      end if
      call MPI_Bcast(offs, one, MPI_OFFSET, master, mpi_mod_comm(comm), ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to broadcast file offset from comm master.', int(ierr), 1)
      end if

      ! create MPI datatype that mirrors layout of your ScaLAPACK block-cyclic decomposition
      call MPI_Type_create_darray(comm_size, comm_rank, ndims, gsizes, distrb, dargs, &
                                  psizes, MPI_ORDER_FORTRAN, MPI_REAL8, datype, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to set up distributed array type.', int(ierr), 1)
      end if
      call MPI_Type_commit(datype, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to commit distributed array type.', int(ierr), 1)
      end if

      ! define data view compatible with the custom datatype
      call MPI_File_set_view(fh, offs, MPI_REAL8, datype, 'native', MPI_INFO_NULL, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to define file view. Crappy MPI implementation?', int(ierr), 1)
      end if

      ! collective write
      call MPI_File_write_all(fh, locA, locsiz, MPI_REAL8, stat, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to write file.', int(ierr), 1)
      end if
      call MPI_File_sync(fh, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to synchronize file after write.', int(ierr), 1)
      end if

      ! cleanup
      call MPI_Type_free(datype, ierr)
      if (ierr /= MPI_SUCCESS) then
         call mpi_xermsg('mpi_mod', 'mpi_mod_file_write_darray2d_real64', &
                         'Failed to release derived MPI type.', int(ierr), 1)
      end if
#else
      write (fh) locA(1:myrows,1:mycols)
#endif
   end subroutine mpi_mod_file_write_darray2d_real64


end module mpi_mod
