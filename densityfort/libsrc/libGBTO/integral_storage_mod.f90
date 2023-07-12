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
!> The object integral_calculation controls the integral evaluation. The user first declares object (type) integral_storage_obj type
!> in his/her program, fills it in with the required parameters and then passes this object as an argument to the add_int method.
!> This adds the specified integral calculation to the list of the integrals to calculate. The actual calculation is then performed
!> using the calculate method.
module integral_storage_mod
use precisn
use const, only: line_len, stdout, int_rel_prec, int_del_thr, no_header
!use common_obj, only: darray_2d
use parallel_arrays, only: p2d_array_obj
use data_file
use utils, only: xermsg

   !> \class <integral_data>
   !> Contains options that control numerical details of the integral calculation. This type can be extended for a particular type
   !> of integral calculations to include the options relevant only for those.
   !> \todo the first line of read/write should include version number of the object to ensure I can tell
   !>       if I am reading things correctly.
   type integral_options_obj
      !> Rel. precision for the resulting integrals. Error will be triggered if there is a suspicion that the integral fails to meet
      !> the precision specified. This is implemented for the GTO/BTO integrals.
      real(kind=cfp) :: prec = int_rel_prec
      !> Threshold magnitude for the integrals: integrals with absolute values smaller than this will be zeroed out.
      real(kind=cfp) :: tol = int_del_thr
      !> R-matrix radius. This value must be always consistent with the extent of the BTO basis. The last endpoint of the radial
      !> B-spline grid must be the value of the R-matrix radius. If no BTOs are present in the continuum basis then this value can
      !> be set to an arbitrary value.
      real(kind=cfp) :: a = 0.0_cfp
      !> Maximum allowed size (in Mib) of the ijrs as allocated by molecular_orbital_basis_obj%two_electron_integrals.
      real(kind=cfp) :: max_ijrs_size = -1.0_cfp
      !> If set to .true. then 2p integrals of the type [continuum,continuum|continuum,continuum] will be calculated.
      !> If set to .false. (default) then 2p integrals with only one particle in the continuum will be calculated.
      logical :: two_p_continuum = .false.
      !> Set to .true. to use the CGTO-only integral algorithm based on the spherical CGTOs.
      !> This is numerically stable even for high L.
      logical :: use_spherical_cgto_alg = .false.
      !> Method to calculate the mixed BTO/CGTO integrals: 1 = Legendre expansion, 2 = Lebedev quadrature.
      !> todo broadcast it too in the parallel calculation!
      integer :: mixed_ints_method = -1
      !> If this variable is set then all intermediate integrals Y_lm needed for calculation of the BTO/GTO integrals
      !> will be saved to disk in the directory specified.
      character(len=line_len) :: scratch_directory = ''
      !> Maximum value of L to be used in the Legendre expansion when calculating the 1-electron, resp. 2-electron mixed integrals.
      !> todo broadcast it too in the parallel calculation!
      integer :: max_l_legendre_1el = -1, max_l_legendre_2el = -1
      !> Grid step for the r1 coordinate used in the calculation of the mixed BTO/GTO integrals.
      real(kind=cfp) :: delta_r1 = 0.25_cfp
      !> The maximum value of the angular momentum of the property operator.
      !> The default is 2, i.e. property integrals up to quadrupoles will be calculated.
      integer :: max_property_l = -1
      !> Set to .true. if calculation of overlap and kinetic energy integrals is required.
      logical :: calculate_overlap_ints = .false.
      !> Set to .true. if calculation of kinetic energy integrals is required.
      logical :: calculate_kinetic_energy_ints = .false.
      !> Set to .true. if calculation of nuclear attraction integrals is required.
      logical :: calculate_nuclear_attraction_ints = .false.
      !> Set to .true. if calculation of one electron Hamiltonian integrals are required
      !> (sum of kinetic energy and nuclear attraction integrals).
      logical :: calculate_one_el_hamiltonian_ints = .false.
      !> Set to .true. if calculation of property integrals (with max_property_l) is required.
      logical :: calculate_property_ints = .false.
      !> Set to .true. if calculation of two electron integrals is required.
      logical :: calculate_two_el_ints = .false.
      !> Request to print the evaluated integrals at the end of the calculation. This parameter is not saved to/read from the disk.
      logical :: print_integrals = .false.
   contains
      !> Checks that the precision values are within a reasonable range.
      procedure :: check => check_integral_options_obj
      !> Reads the integral_options data structure from the disk given the starting position of the record.
      procedure :: read => read_integral_options_obj
      !> Writes the integral_options data structure to the disk given the starting position of the record.
      procedure :: write => write_integral_options_obj
      !> Prints to screen the contents of this object.
      procedure :: print => print_integral_options_obj
   end type integral_options_obj 

   !private routines
   private check_integral_options_obj, read_integral_options_obj, write_integral_options_obj, print_integral_options_obj

   !> \class <integral_storage_obj>
   !> This object provides access to the results of an integral calculation. The object must be associated with the data which can
   !> be either in memory or on the disk. The association with the target is performed calling 'init' with the appropriate arguments.
   !> If the results are stored in memory then the array pointed to by the variable 'integral' must be associated with a target.
   !> If the integrals are stored on the disk then the object integral_file is used to provide access to the data.
   !> Note that neither integral_file nor integral are private variables. This is neccessary so that the integral routines can have
   !> direct (i.e. fast access) to the data.
   !> \warning The type-bound methods which require on input/output an integral file (or integrals in memory) should always pass the
   !> reference to this file (memory location) via the integral_storage_obj associated with this file (memory). The file name or
   !> integral array for input/output should never be passed directly to the methods. The reason is that the input/output may be
   !> associated with another integral_storage_obj. Modifying the file (memory) without involving the associated integral_storage_obj
   !> in the process can result in the associated integral_storage_obj getting confused at some point or even corrupting the data
   !> on the file! In any case it is strongly reccommended that all interaction with the stored data is performed implementing or
   !> extending type-bound methods for this object rather than accessing the disk file directly.
   type integral_storage_obj
      !> Object used to manage access to the data stored on the disk.
      type(data_file_obj) :: integral_file
      !> Pointer to the array where the integrals are stored.
      !> The second dimension is used to distinguish integrals of different type. E.g. we usually calculate the overlap and
      !> the kinetic energy integrals at once, so the first column may be used for the overlap integrals, while the second one
      !> for the kinetic energy ones.
      class(p2d_array_obj), pointer :: integral => null()
      !> Header for the data stored in memory.
      type(data_header_obj) :: data_header
      !> This variable is set following the call to 'init'.
      logical, private :: initialized = .false.
      !> These variables are set following the call to 'init'. They tell me whether the data are stored in memory or on the disk.
      logical, private :: memory = .false., disk = .false.
   contains
      !> This routines initializes the object. Depending on whether the pointer to the output array or a file name are given
      !> on input the variable 'integral' or 'integral_file' is set. The routine returns 0 on successful initialization.
      !> Non-zero result means an error.
      procedure :: init => init_integral_storage_obj
      !> This routine restores the state of the object to the default state, i.e. if the storage is in memory it nullifies
      !> the association of the 'integral' pointer. If the storage is on disk then the integral file is closed.
      procedure :: final => final_integral_storage_obj
      !> This function checks if the data have meaningful values.
      procedure :: check => check_integral_storage_obj
      !> This routine reads into memory the integrals with the given header and the integral options from the integral_storage_obj
      !> given. There is a good reason why the integral file is not given directly but via the integral_storage_obj - see the object
      !> description. Note that by default the integral arrays are read-in using the shared-memory feature so in the MPI mode there
      !> will be only one copy of the integral array PER NODE. If the shared-memory capability is not available then the code will
      !> resort to the local mode where each MPI task keeps its own copy of the integral array.
      procedure, non_overridable :: read => read_integrals
      !> This routine writes the integrals and the integral options to the disk. There is a good reason why the path to the integral
      !> file is not given directly as an argument but indirectly via the integral_storage_obj - see the object description.
      procedure, non_overridable :: write => write_integrals
      !> Returns the value of 'memory'.
      procedure, non_overridable :: in_memory => get_memory
      !> Returns the value of 'disk'.
      procedure, non_overridable :: on_disk => get_disk
      !> Given the string defining a basis set and a string defining the integral type it returns a single string which merges these
      !> two together without additional spaces in between.
      procedure, non_overridable :: contruct_header_string
   end type integral_storage_obj

   private check_integral_storage_obj, init_integral_storage_obj, read_integrals, write_integrals, final_integral_storage_obj
   private get_memory, get_disk

contains

   !todo put the 10d-10 value to the const module
   function check_integral_options_obj(this)
      implicit none
      class(integral_options_obj) :: this
      integer :: check_integral_options_obj

         check_integral_options_obj = 0
     
         ! relative precision of the integrals must be a positive number not greater than 10d-10:
         ! we don't allow precision worse than that
         if (this%prec .le. 0.0_cfp .or. this%prec > 10e-10_cfp) then
            check_integral_options_obj = 1
            return
         endif

         ! absolute value of the smallest integral to retain must be positive or zero (retain all integrals).
         if (this%tol < 0.0_cfp) then
            check_integral_options_obj = 2
            return
         endif

         if (this%max_property_l < 0 .and. this%calculate_property_ints) then
            check_integral_options_obj = 3
            return
         elseif (this%max_property_l .ge. 0 .and. .not.(this%calculate_property_ints)) then
            check_integral_options_obj = 4
            return
         endif

   end function check_integral_options_obj

   subroutine read_integral_options_obj(this,unit,record_start,position_after_read)
      use mpi_mod
      implicit none
      class(integral_options_obj) :: this
      integer, intent(in) :: unit, record_start
      integer, intent(out) :: position_after_read
      integer :: dummy

         if (record_start <= 0) then
            call xermsg ('integral_options_obj', 'read_integral_options_obj', 'On input the start of the record was .le. 0', 1, 1)
         end if

         if (myrank .eq. master) then
            read(unit,pos=record_start,err=10) this%prec, this%tol
            read(unit) this%a
            read(unit) this%max_ijrs_size
            read(unit) dummy; this%two_p_continuum = dummy/=0
            read(unit) dummy; this%use_spherical_cgto_alg = dummy/=0
            read(unit) this%max_property_l
            read(unit) dummy; this%calculate_overlap_ints = dummy/=0
            read(unit) dummy; this%calculate_kinetic_energy_ints = dummy/=0
            read(unit) dummy; this%calculate_nuclear_attraction_ints = dummy/=0
            read(unit) dummy; this%calculate_one_el_hamiltonian_ints = dummy/=0
            read(unit) dummy; this%calculate_property_ints = dummy/=0
            read(unit) dummy; this%calculate_two_el_ints = dummy/=0

            inquire(unit,pos=position_after_read)
         endif

         call mpi_mod_bcast(this%prec,master)
         call mpi_mod_bcast(this%tol,master)
         call mpi_mod_bcast(this%a,master)
         call mpi_mod_bcast(this%two_p_continuum,master)
         call mpi_mod_bcast(this%use_spherical_cgto_alg,master)
         call mpi_mod_bcast(this%max_property_l,master)
         call mpi_mod_bcast(this%calculate_overlap_ints,master)
         call mpi_mod_bcast(this%calculate_kinetic_energy_ints,master)
         call mpi_mod_bcast(this%calculate_nuclear_attraction_ints,master)
         call mpi_mod_bcast(this%calculate_one_el_hamiltonian_ints,master)
         call mpi_mod_bcast(this%calculate_property_ints,master)
         call mpi_mod_bcast(this%calculate_two_el_ints,master)
         call mpi_mod_bcast(position_after_read,master)
         call mpi_mod_bcast(this%max_ijrs_size,master)

         return

10       call xermsg('integral_options_obj','read_integral_options_obj','Error while executing the read command.',2,1)

   end subroutine read_integral_options_obj

   subroutine write_integral_options_obj(this,unit,record_start,position_after_write)
      use mpi_mod
      implicit none
      class(integral_options_obj) :: this
      integer, intent(in) :: unit, record_start
      integer, intent(out) :: position_after_write

         if (record_start <= 0) then
            call xermsg ('integral_options_obj', 'write_integral_options_obj', 'On input the start of the record was .le. 0', 1, 1)
         end if

         if (myrank .eq. master) then
            write(unit,pos=record_start,err=10) this%prec, this%tol
            write(unit) this%a
            write(unit) this%max_ijrs_size
            write(unit) merge(1, 0, this%two_p_continuum)
            write(unit) merge(1, 0, this%use_spherical_cgto_alg)
            write(unit) this%max_property_l
            write(unit) merge(1, 0, this%calculate_overlap_ints)
            write(unit) merge(1, 0, this%calculate_kinetic_energy_ints)
            write(unit) merge(1, 0, this%calculate_nuclear_attraction_ints)
            write(unit) merge(1, 0, this%calculate_one_el_hamiltonian_ints)
            write(unit) merge(1, 0, this%calculate_property_ints)
            write(unit) merge(1, 0, this%calculate_two_el_ints)

            inquire(unit,pos=position_after_write)
         endif

         !send the current position in the file (after write) to all processes
         call mpi_mod_bcast(position_after_write,master)

         return

10       call xermsg('integral_options_obj','write_integral_options_obj','Error while executing the write command.',2,1)

   end subroutine write_integral_options_obj

   subroutine print_integral_options_obj(this)
      implicit none
      class(integral_options_obj) :: this
         
         write(stdout,'("--------->","integral_options_obj:print")') 
         write(stdout,'("Requested minimal relative precision: ",e25.15)') this%prec
         write(stdout,'("Requested threshold for the magnitude of the retained integrals: ",e25.15)') this%tol
         write(stdout,'("Requested value of the R-matrix radius: ",e25.15)') this%a
         write(stdout,'("Maximum allowed size of the temporary ijrs array (Mib): ",f10.3)') this%max_ijrs_size
         if (this%two_p_continuum) then
            write(stdout,'("Integrals including TWO PARTICLES in the continuum selected.")')
         else
            write(stdout,'("Integrals including ONE PARTICLE in the continuum selected.")')
         endif
         write(stdout,'("<---------","done:integral_options_obj:print")')         

   end subroutine print_integral_options_obj

   function init_integral_storage_obj(this,memory,disk)
      implicit none
      class(integral_storage_obj) :: this
      class(p2d_array_obj), target, optional, intent(in) :: memory
      character(len=*), optional, intent(in) :: disk
      integer :: init_integral_storage_obj

         this % initialized = .false.
 
         init_integral_storage_obj = 0

         !Check that the required storage location is given consistently
         if (present(memory) .and. present(disk) .or. ( .not.(present(memory)) .and. .not.(present(disk)) )) then
            init_integral_storage_obj = 1
            return
         endif

         if (present(memory)) then
            this%integral => memory
            this%memory = .true.
            this%disk = .false.
         endif

         if (present(disk)) then

            !open the file for stream access and read-in the existing (if any) headers on the file.
            call this%integral_file%open(disk)

            this%memory = .false.
            this%disk = .true.
         endif

         this%initialized = .true.

   end function init_integral_storage_obj

   subroutine final_integral_storage_obj(this)
      implicit none
      class(integral_storage_obj) :: this
      type(data_header_obj) :: default_header

         if (.not. this % initialized) then
            call xermsg ('integral_storage_obj', 'final_integral_storage_obj', 'The object has not been initialized.', 1, 1)
         end if

         if (this%disk) call this%integral_file%close
         if (this%memory) this%integral => null()

         this%disk = .false.
         this%memory = .false.
         this%data_header = default_header
         this%initialized = .false.

   end subroutine final_integral_storage_obj

   function check_integral_storage_obj(this)
      implicit none
      class(integral_storage_obj) :: this
      integer :: check_integral_storage_obj

         check_integral_storage_obj = 0

         !we can have data either on the disk or in the memory but not both
         if (this%memory .and. this%disk) then
            check_integral_storage_obj = 1
            return
         endif

         if (this%memory) then
            !The integral array must be associated with a target
            if (.not.associated(this%integral)) then
               check_integral_storage_obj = 2
               return
            endif

            !The header must be specified
            if (this%data_header%name .eq. no_header) then
               check_integral_storage_obj = 3
               return
            endif
         endif

   end function check_integral_storage_obj

   subroutine read_integrals(this,src,header,tgt_int_opt,tgt_is_local)
      use mpi_mod
      implicit none
      class(integral_storage_obj) :: this
      character(len=*), intent(in) :: header
      class(integral_options_obj), intent(out) :: tgt_int_opt
      class(integral_storage_obj), intent(in) :: src
      logical, intent(in) :: tgt_is_local

      integer :: err, first_record, last_record, lunit, pos_in_file
      character(line_len) :: head, file_name

         if (.not. this % initialized) then
            call xermsg ('integral_storage_obj', 'read_integrals', 'The object has not been initialized.', 1, 1)
         end if

         if (this % disk) then
            call xermsg ('integral_storage_obj', 'read_integrals', &
                         'The storage has been associated with disk but reading into memory has been requested.', 2, 1)
         end if

         err = src%check()
         if (err /= 0) then
            call xermsg ('integral_storage_obj', 'read_integrals', 'Check of the source integral storage has failed.', 3, 1)
         end if

         head = header

         err = src%integral_file%find_header(head,first_record,last_record)
         if (err /= 0) then
            call xermsg ('integral_storage_obj', 'read_integrals', &
                'Searching for the requested header has been terminated with an error code. See data_file_obj for details.', 4, 1)
         end if

         !The following may occur if writing of the data record has not been finished or the record is corrupt.
         if (first_record <= 0 .or. last_record <= 0) then
            call xermsg ('integral_storage_obj', 'read_integrals', 'The requested data record is missing or not complete.', 5, 1)
         end if

         lunit = src%integral_file%get_unit_no()
         file_name = src%integral_file%get_file_name()

         !HERE WE READ-IN THE ACTUAL INTEGRAL DATA INTO THE integral_storage_obj GIVEN BY src

         !the first record are the integral_options
         call tgt_int_opt%read(lunit,first_record,pos_in_file)

         err = tgt_int_opt%check()
         if (err /= 0) then
            call xermsg ('integral_storage_obj', 'read_integrals', &
                         'Check of the read integral options structure has failed.', err, 1)
         end if

         call tgt_int_opt%print

         !the second record are the actual integrals: these are either copied to each process or scattered if sharing is required.
         first_record = pos_in_file
         if (tgt_is_local) then
            call this%integral%read(file_name,lunit,first_record,pos_in_file)
         else
            call xermsg('integral_storage_obj','read_integrals','Use of tgt_is_local is obsolete and not implemented anymore.',7,1)
         endif

         if (pos_in_file .ne. last_record) then
            print *,pos_in_file,last_record
            call xermsg ('integral_storage_obj', 'read_integrals', &
                         'The header data is corrupt: the position of the last record given by data_header is not correct.', 8, 1)
         endif

         this%data_header%name = head

         return

 10      call xermsg ('integral_storage_obj', 'read_integrals', &
                      'Error executing the read command while reading the integral data.', 9, 1)

   end subroutine read_integrals

   subroutine write_integrals(this,src,header,src_int_opt)
      use mpi_mod
      implicit none
      class(integral_storage_obj) :: this
      character(len=*), intent(in) :: header
      class(integral_storage_obj), intent(in) :: src
      class(integral_options_obj), intent(in) :: src_int_opt

      integer :: first_record, last_record, lunit, i, err, pos_in_file
      character(line_len) :: head

         if (.not. this % initialized) then
            call xermsg ('integral_storage_obj', 'write_integrals', &
                         'The object has not been initialized.', 1, 1)
         end if

         if (this % memory) then
            call xermsg ('integral_storage_obj', 'write_integrals', &
                         'The storage has been associated with memory but writing into disk has been requested.', 2, 1)
         end if

         err = src%check()
         if (err /= 0) then
            call xermsg ('integral_storage_obj', 'write_integrals', &
                         'Check of the source integral storage has failed.', 3, 1)
         end if

         err = src_int_opt%check()
         if (err /= 0) then
            call xermsg ('integral_storage_obj', 'write_integrals', &
                         'Check of the input integral options structure has failed.', 4, 1)
         end if
 
         head = header

         if (head /= src % data_header % name) then
            call xermsg ('integral_storage_obj', 'write_integrals', &
                         'The requested header has not been found on the source integral storage.', 5, 1)
         end if

         lunit = this%integral_file%get_unit_no()

         !HERE WE WRITE THE ACTUAL INTEGRAL DATA INTO THE integral_storage_obj (associated with a disk file) GIVEN BY tgt
         first_record = this%integral_file%start_record(head)

         !The first record are the integral_options
         call src_int_opt%write(lunit,first_record,pos_in_file)

         ! The second record are the actual integrals. These are written in the mode corresponding to the mode
         ! of the src%integral object that holds the integrals. Only master writes its data to the disk.
         i = master
         call src%integral%write(lunit,pos_in_file,last_record,i)

         call this%integral_file%close_record(head,first_record,last_record)

         return

 10      call xermsg ('integral_storage_obj', 'write_integrals', &
                      'Error executing the write command while writing the integral data.', 8, 1)

   end subroutine write_integrals

   function get_memory(this)
      implicit none
      class(integral_storage_obj) :: this
      logical :: get_memory

         if (.not.(this%initialized)) call xermsg('integral_storage_obj','get_memory','The object has not been initialized.',1,1)

         get_memory = this%memory

   end function get_memory

   function get_disk(this)
      implicit none
      class(integral_storage_obj) :: this
      logical :: get_disk

         if (.not.(this%initialized)) call xermsg('integral_storage_obj','get_disk','The object has not been initialized.',1,1)

         get_disk = this%disk

   end function get_disk

   function contruct_header_string(this,basis_name,integral_name)
      use const, only: line_len
      implicit none
      class(integral_storage_obj) :: this
      character(len=*), intent(in) :: basis_name,integral_name
      character(len=line_len) :: contruct_header_string

         contruct_header_string = trim(basis_name)//':'//trim(integral_name) 

   end function contruct_header_string

end module integral_storage_mod
