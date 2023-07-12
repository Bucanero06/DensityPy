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
module atomic_basis_mod
   use precisn
   use const, only: line_len, no_header, stdout, sym_op_nam_len
   use symmetry
   use integral_storage_mod
   use mpi_mod
   use basis_data_generic_mod

   implicit none

   private

   !Only the objects themselves are visible from outside of this module.
   public atomic_orbital_basis_obj, CGTO_shell_data_obj, BTO_shell_data_obj

   !> This object is used to include all kinds of auxiliary data needed during the integral calculation.
   type :: integral_data_obj
       !> This can be arbitrary but in practice centers other than CMS are rarely required.
       real(kind=cfp) :: property_center(1:3) = 0.0_cfp
       !> Number of different types of integrals to calculate.
       integer :: number_of_types = 0
       !> This array is used to identify the types of integrals to calculate. It has size number_of_types.
       character(len=line_len), allocatable :: column_descriptor(:)
       !> The largest L-value in the BTO basis and the smallest index of the radial B-spline included.
       integer :: max_bspline_l = -1, first_bspline_index = -1
       !> Column indices corresponding to the different types of integrals to be calculated.
       integer :: olap_column = 0, kei_column = 0, prop_column = 0, nai_column = 0, one_elham_column = 0, two_el_column = 0
   end type integral_data_obj

   type, extends(basis_data_generic_obj) :: atomic_orbital_basis_obj
       !> n(1) gives the number of CGTO shells in the basis and n(2) gives the number of BTOs in the basis.
       integer, private :: n(2) = 0
       !> Shell data for CGTOs. This array must have size n(1).
       type(CGTO_shell_data_obj), allocatable, private :: CGTO_shell_data(:)
       !> Shell data for BTOs. This array must have size n(2).
       type(BTO_shell_data_obj), allocatable, private :: BTO_shell_data(:)
       !> Array of dimensions (2,number_of_basis_functions) containing mapping of basis function indices to shell indices (1) and indices within the shells (2).
       !> This table is useful for looking up which shells need to be involved in calculation of a given integral and for picking out the requested integral
       !> from the integral batch generated over shells of functions.
       integer, allocatable :: indices_to_shells(:,:)
       !> For each shell (2nd dimension) the first dimension contains the following info: shell type (CGTO: 1, BTO: 2), row 1, index of the shell within its own type, row 2,
       !> is a continuum shell (1) or not (0), row 3, starting index for the functions in the shell, row 4, the number of functions in the shell, row 5 and the angular momentum of the shell, row 6.
       !> The order in which the mapping is stored corresponds to the order in which the shells were added to the basis.
       integer, allocatable :: shell_descriptor(:,:)
       !> ordered_pairs: Indices mapping pairs of functions to a single index. This mapping is required to index properly 2p integrals with only 1p in the continuum.
       !> ordered_pairs_to_indices: mapping of the ordered pair indices to the indices of the two corresponding basis functions.
       integer, allocatable :: ordered_pairs(:,:), ordered_pairs_to_indices(:,:)
       !> Number of target functions, index of the last CT pair of functions and other values needed to index optimally 2-electron 1p continuum integrals.
       integer :: n_target_fns = -1, n_cont_fns = -1, last_CT_fn = -1, n_prec_ints = -1, n_TT_pairs = -1
       !> Indices mapping pairs of shells to a single index. This mapping is required to index properly 2p integrals with only 1p in the continuum in case of multiple MPI tasks being used.
       integer, allocatable :: ordered_shell_pairs(:,:)
       !> Number of target shells, index of the last CT pair of shells and other values needed to index optimally 2-electron 1p continuum integrals in case of multiple MPI tasks are used.
       integer :: n_target_sh = -1, n_cont_sh = -1, last_CT_sh = -1, n_prec_sh = -1, n_TT_sh_pairs = -1
       !> shell_pair_indices: Mapping between indices of pairs of shells with the indices of the shells within their own type of shells of basis functions.
       !> shell_pair_type: Mapping between pairs of shells and their type: GG (1), BG (2), GB (3), BB (4).
       integer, allocatable :: shell_pair_indices(:,:), shell_pair_type(:)
       !> Total number of CGTO functions in the basis.
       integer :: n_cgto_functions = -1
       !> This data structure is used to include all kinds of auxiliary data needed during the integral calculation.
       type(integral_data_obj), private :: integral_data
       !> Auxiliary arrays used by eval_orbital.
       real(kind=cfp), allocatable :: ao_basis_at_r(:,:), r(:,:)
       !> Used only during initialization while the shells are being added.
       integer, private :: space_allocated = 0
       !> Set to .true. following input of all basis functions for which space has been allocated.
       logical, private :: initialized = .false.
       !> Set to .true. following a call to init.
       logical, private :: init_called = .false.
       !> Set to .true. in add_shell each time a continuum function is added. This serves the purpose of checking that all target functions preceed the continuum ones.
       logical, private :: continuum_added = .false.
   contains
       !> Allocates space for a given number and type of shells of basis functions.
       procedure :: init
       !> Finalizes the basis set.
       procedure :: final
       !> Adds any function which is a primitive gto (single exponent) located at the center of mass, uses add_shell.
       procedure :: add_cms_gtos
       !> Adds data for one shell into the basis set. All target shells must be added before the continuum ones. This is important for the integral generation algorithms where this situation simplifies a lot
       !> what would otherwise be a complicated book-keeping. Internally, the norms of the contracted GTOs are multiplied in with the primitive norms. This strategy is used in the whole basis set to simplify
       !> things.
       procedure :: add_shell
       !> Routine that calculates all 1-electron integrals over a given combination of shells of functions. Note that this routine can only be called once the this%integral_data has been prepared.
       procedure, private :: shell_pair_one_electron_integrals
       !> Routine that calculates 2-electron integrals over a given combination of shells of functions. Note that this routine can only be called once the this%integral_data has been prepared.
       procedure, private :: shell_quartet_two_electron_integrals
       !> Generates a starting index for each shell in the basis.
       procedure, private :: generate_shell_indices
       !> Generates the indices n_target_fns, n_cont_fns, etc. and prints the summary table of basis function indices and their types.
       procedure, private :: generate_fn_indices
       !> Calculates and stores 1-electron integrals for all pairs of shells in the basis.
       procedure :: one_electron_integrals
       !> Calculates and stores 2-electron integrals for all pairs of shells in the basis.
       procedure :: two_electron_integrals
       !> Calculates index for the given type of integral. The integral type is given by one of the strings definied in the const module and by the number of basis functions on input.
       procedure :: integral_index
       !> Takes the p2d_array of integrals on input containing the overlap integrals and returns the 2D array corresponding to the AO overlap matrix.
       procedure :: assemble_overlap_matrix
       !> Returns the name of the basis set.
       procedure :: get_basis_name
       !> Returns the name of the i-th shell in the basis set.
       procedure :: get_shell_name
       !> Returns the shell data for the i-th shell in the basis set.
       procedure :: get_shell_data
       !> Returns an array containing data for all shells of CGTOs in the basis.
       procedure :: get_all_CGTO_shells
       !> Returns an array containing data for all shells of BTOs in the basis.
       procedure :: get_all_BTO_shells
       !> Returns the object bspline_grid_obj if the basis contains shells of BTO functions.
       procedure :: get_bspline_grid
       !> Returns the number of functions in a shell with the given index.
       procedure :: get_number_of_functions_in_shell
       !> Returns the value of initialized.
       procedure :: is_initialized
       !> For a given IRR it returns a list of logical values of size equal to the number of functions in the basis.
       !> i-th element of the output array is set to .true. if the i-th function in the basis has symmetry irr andcoms from a shell between indices from and to.
       procedure :: get_symmetry_flags
       !> For a given IRR it returns a list of logical values of size equal to the number of functions in the basis.
       !> i-th element of the output array is set to .true. if the i-th function in the basis is a continuum function with the given symmetry.
       procedure :: get_continuum_flags
       !> Returns the values of the smallest and largest L of the continuum functions in the basis. If the basis doesn't contain continuum then it terminates with an error.
       procedure :: get_continuum_l_range
       !> Prints the basis set data to stdout.
       procedure :: print => print_atomic_orbital_basis_obj
       !> Calculates amplitudes for all basis functions in the basis and builds the channel numbers. The continuum is first normalized to the requested radius.
       !> The additional normalization factor for CGTOs is multiplied in with the primitive norms.
       procedure :: calculate_amplitudes
       !> Returns .true. if the basis contains the BTOs.
       procedure :: contains_btos
       !> Normalizes the atomic continuum functions to the R-matrix sphere with a given radius. The normalization is only done if a > 0.0_cfp.
       procedure, private :: normalize_continuum
       !> Evaluates the AO basis for the given set of points in space.
       procedure :: eval_basis
   end type atomic_orbital_basis_obj

 contains

   function init(this,n,geometry)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: n
      class(geometry_obj), intent(in) :: geometry
      integer :: init

         write(stdout,'("--------->","atomic_orbital_basis_obj:init")')

         init = 0

         if (this%initialized) then
            init = this%final()
            if (init .ne. 0) then
               call xermsg('atomic_orbital_basis_obj', 'init', &
                           'Finalization has failed. See atomic_orbital_basis_obj%final for details.', init, 0)
               return
            endif
         endif

         if (n < 0) then
            init = 1
            call xermsg('atomic_orbital_basis_obj','init','On input the value of n was out of range.',init,0)
            return
         endif

         init = this%symmetry_data%init(geometry)
         if (init .ne. 0) then
            call xermsg('atomic_orbital_basis_obj', 'init', &
                        'Symmetry initialization failed. See symmetry_obj%init for details.', init, 0)
            return
         endif

         this%space_allocated = n
         allocate(this%shell_descriptor(6,this%space_allocated),stat=init)
         if (init .ne. 0) then
            call xermsg('atomic_orbital_basis_obj','init','Memory allocation of this%shell_descriptor has failed.',init,0)
            return
         endif
         this%shell_descriptor = 0

         this%n = 0
         this%number_of_shells = 0
         this%init_called = .true.
         this%continuum_added = .false.
         this%number_of_functions = 0
         this%n_cgto_functions = 0

         write(stdout,'("<---------","atomic_orbital_basis_obj:init")')

   end function init

   function final(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer :: final

         write(stdout,'("--------->","atomic_orbital_basis_obj:final")')

         final = 0

         if (allocated(this%CGTO_shell_data)) deallocate(this%CGTO_shell_data)
         if (allocated(this%BTO_shell_data)) deallocate(this%BTO_shell_data)
         this%number_of_shells = 0
         this%number_of_functions = 0
         if (allocated(this%shell_descriptor)) deallocate(this%shell_descriptor)
         if (allocated(this%indices_to_shells)) deallocate(this%indices_to_shells)
         if (allocated(this%ordered_pairs)) deallocate(this%ordered_pairs)
         if (allocated(this%ordered_pairs_to_indices)) deallocate(this%ordered_pairs_to_indices)
         if (allocated(this%ordered_shell_pairs)) deallocate(this%ordered_shell_pairs)
         if (allocated(this%ao_basis_at_r)) deallocate(this%ao_basis_at_r)
         if (allocated(this%r)) deallocate(this%r)
         this%n = 0
         this%space_allocated = 0
         this%initialized = .false.
         this%continuum_added = .false.
         this%init_called = .false.

         write(stdout,'("<---------","atomic_orbital_basis_obj:final")')

   end function final

   subroutine add_cms_gtos(this,min_l,max_l,max_num,no_exps,exponents,n_gto,non_zero_at_boundary)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      type(CGTO_shell_data_obj)  :: shell
      real(kind=cfp), intent(in) :: exponents(1:max_num,min_l:max_l)
      integer, intent(in)    :: min_l, max_l, no_exps(min_l:max_l),max_num
      integer, intent(out) :: n_gto
      integer :: i, j
      logical, intent(in) :: non_zero_at_boundary

         write(stdout,'("--------->","atomic_orbital_basis_obj:add_cms_gtos")')

         ! Add the gto functions as single primitive CGTOs to this object
         n_gto = 0
         do i=min_l,max_l
            if (i < 0) exit

            shell%l = i
            shell%number_of_functions = 2*shell%l+1
            shell%center(1:3) = (/0.0_cfp,0.0_cfp,0.0_cfp/)
            shell%number_of_primitives = 1
            call shell%make_space(shell%number_of_primitives)

            do j=1,no_exps(i)
               if (exponents(j,i) < 0.0_cfp ) cycle
               shell%exponents(1:shell%number_of_primitives) = exponents(j,i)
               shell%contractions(1:shell%number_of_primitives) = 1.0_cfp
               shell%non_zero_at_boundary = non_zero_at_boundary
               call this%add_shell(shell)
               n_gto = n_gto + shell%number_of_functions
            enddo !j

         enddo !i

         write(stdout,'("<---------","atomic_orbital_basis_obj:add_cms_gtos")')

   end subroutine add_cms_gtos

   subroutine add_shell(this,shell_data)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      class(shell_data_obj), intent(inout) :: shell_data

      integer :: i, err
      type(CGTO_shell_data_obj), allocatable :: temp_CGTO_shell_data(:)
      type(BTO_shell_data_obj), allocatable :: temp_BTO_shell_data(:)

         write(stdout,'("--------->","atomic_orbital_basis_obj:add_shell")')

         if (.not. this % init_called) then
            call xermsg('atomic_orbital_basis_obj', 'add_shell', 'Attempt to call add_shell before calling init.', 2, 1)
         end if

         if (this % initialized) then
            call xermsg('atomic_orbital_basis_obj', 'add_shell', &
                        'All shells of functions for which space has been allocated have already been supplied.', 1, 1)
         end if

         call shell_data%normalize
         call shell_data%print

         select type (shell => shell_data)
            type is (CGTO_shell_data_obj)

               !Resize the basis of CGTOs:
               call move_alloc(this%CGTO_shell_data,temp_CGTO_shell_data)
               this%n(1) = this%n(1)+1
               allocate(this%CGTO_shell_data(this%n(1)),stat=err)
               if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'add_shell', 'Memory allocation 1 failed.',err, 1)

               !Copy the data for the previous CGTOs:
               do i=1,this%n(1)-1
                  this%CGTO_shell_data(i) = temp_CGTO_shell_data(i)
               enddo !i
               if (this%n(1) > 1) deallocate(temp_CGTO_shell_data)

               if (shell%non_zero_at_boundary) then
                  write(stdout,'("Atomic continuum shell of type CGTO_shell_data_obj has been added &
                                 &to the atomic basis as shell number: ",i4)') this % n(1)
                  this%continuum_added = .true.
               else
                  if (this % continuum_added) then
                     call xermsg('atomic_orbital_basis_obj', 'add_shell', &
                                 'Attempt to add a CGTO target shell beyond a continuum shell.', 2, 1)
                  end if
                  write(stdout,'("Atomic target shell of type CGTO_shell_data_obj has been added &
                                 &to the atomic basis as shell number: ",i4)') this % n(1)
               endif

               !Add the new shell data to the list:
               this%CGTO_shell_data(this%n(1)) = shell

               !Multiply-in the contraction norm with the primitive norms:
               this%CGTO_shell_data(this % n(1)) % norms(:) = this % CGTO_shell_data(this % n(1)) % norm &
                                                            * this % CGTO_shell_data(this % n(1)) % norms(:)
               this%CGTO_shell_data(this%n(1))%norm = 1.0_cfp

               !Save the relative index of this shell and mark it as a CGTO shell
               this%number_of_shells = this%number_of_shells+1

               this%shell_descriptor(1,this%number_of_shells) = 1
               this%shell_descriptor(2,this%number_of_shells) = this%n(1)
               if (this%CGTO_shell_data(this%n(1))%is_continuum()) this%shell_descriptor(3,this%number_of_shells) = 1
               this%shell_descriptor(6,this%number_of_shells) = this%CGTO_shell_data(this%n(1))%l

               this%n_cgto_functions = this%n_cgto_functions + shell%number_of_functions
            type is (BTO_shell_data_obj)

               !Resize the basis of BTOs:
               call move_alloc(this%BTO_shell_data,temp_BTO_shell_data)
               this%n(2) = this%n(2)+1
               allocate(this%BTO_shell_data(this%n(2)),stat=err)
               if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'add_shell', 'Memory allocation 2 failed.',err, 1)

               !Copy the data for the previous BTOs:
               do i=1,this%n(2)-1
                  this%BTO_shell_data(i) = temp_BTO_shell_data(i)
               enddo !i
               if (this%n(2) > 1) deallocate(temp_BTO_shell_data)

               if (shell%non_zero_at_boundary) then
                  write(stdout,'("Atomic continuum shell of type BTO_shell_data_obj has been added to the atomic basis.")')
                  this%continuum_added = .true.
               else
                  if (this%continuum_added) then
                     call xermsg('atomic_orbital_basis_obj', 'add_shell', &
                                 'Attempt to add a BTO target shell beyond a continuum shell.', 3, 1)
                  end if
                  write(stdout,'("Atomic target shell of type BTO_shell_data_obj has been added to the atomic basis.")')
               endif

               !Add the shell data to the list
               !todo we should check here that the B-spline grid is always the same and equal to the grid set for the first BTO shell.
               !todo check for duplicities in the bspline_index
               this%BTO_shell_data(this%n(2)) = shell
               !Save the relative index of this shell and mark it as a BTO shell
               this%number_of_shells = this%number_of_shells+1

               this%shell_descriptor(1,this%number_of_shells) = 2
               this%shell_descriptor(2,this%number_of_shells) = this%n(2)
               if (this%BTO_shell_data(this%n(2))%is_continuum()) this%shell_descriptor(3,this%number_of_shells) = 1
               this%shell_descriptor(6,this%number_of_shells) = this%BTO_shell_data(this%n(2))%l
            class default
               call xermsg ('atomic_orbital_basis_obj', 'add_shell', &
                  'The shell type must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 4, 1)
         end select

         this%shell_descriptor(5,this%number_of_shells) = shell_data%number_of_functions
         this%number_of_functions = this%number_of_functions + shell_data%number_of_functions

         if (this%number_of_shells .eq. this%space_allocated) then
            this%initialized = .true.

            call this%generate_shell_indices
            call this%generate_fn_indices

         endif

         write(stdout,'("<---------","atomic_orbital_basis_obj:add_shell")')

   end subroutine add_shell

   subroutine generate_fn_indices(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this

      integer :: i,j,err,ij,ind
      logical :: found
      logical, allocatable :: is_continuum(:), is_continuum_irr(:)

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'generate_fn_indices', &
                         'The object has not been initialized or not all shells have been added.', 1, 1)
         end if

         write(stdout,'(/,"All shells have been added. &
                          &Mapping of basis function indices to shells and relative indices within shells follows:")')

         !Determine which functions are the continuum ones:
         allocate(is_continuum(this%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'generate_fn_indices', 'Memory allocation 3 failed.',err, 1)

         is_continuum = .false.
         do i=1,this%symmetry_data%get_no_irrep(this%symmetry_data%get_pg())
            call this%get_continuum_flags(i,is_continuum_irr)
            do j=1,this%number_of_functions
               if (is_continuum_irr(j)) is_continuum(j) = .true.
            enddo
         enddo !i

         !Generate mapping of basis function indices to shell indices (1) and indices within the shells (2):
         if (allocated(this%indices_to_shells)) deallocate(this%indices_to_shells)
         allocate(this%indices_to_shells(2,this%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'generate_fn_indices', 'Memory allocation 3 failed.',err, 1)

         do i=1,this%number_of_functions
            found = .false.
            do j=1,this%number_of_shells-1
               !Find the shell of which this function is a member:
               if (i .ge. this%shell_descriptor(4,j) .and. i < this%shell_descriptor(4,j+1)) then
                  this%indices_to_shells(1,i) = j
                  this%indices_to_shells(2,i) = i-this%shell_descriptor(4,j)+1
                  found = .true.
                  exit
               endif
            enddo !j
            if (.not.(found)) then
               !Is this function member of the last shell?
               if (i .ge. this%shell_descriptor(4,this%number_of_shells-1) .and. i .le. this%number_of_functions) then
                  this%indices_to_shells(1,i) = this%number_of_shells
                  this%indices_to_shells(2,i) = i-this%shell_descriptor(4,this%number_of_shells)+1
               else
                  call xermsg ('atomic_orbital_basis_obj', 'generate_fn_indices', &
                               'Error generating mapping from basis function indices to shells.', 3, 1)
               endif
            endif
            write(stdout,'(i5,1X,i5,1X,i5,1X,l)') i,this%indices_to_shells(1:2,i),is_continuum(i)
         enddo !i

         !Calculate indices for pairs of functions such that the CC pairs
         !have the largest indices: this ensures optimal indexing in case
         !2p integrals with only 1p in the continuum integrals are requested.
         if (allocated(this%ordered_pairs)) deallocate(this%ordered_pairs)
         if (allocated(this%ordered_pairs_to_indices)) deallocate(this%ordered_pairs_to_indices)
         i = this%number_of_functions*(this%number_of_functions+1)/2
         allocate(this%ordered_pairs(2,i),this%ordered_pairs_to_indices(2,i),stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'generate_fn_indices', 'Memory allocation 4 failed.',err, 1)

         this%n_cont_fns = count(is_continuum)
         this%n_target_fns = this%number_of_functions - this%n_cont_fns
         this%n_TT_pairs = this%n_target_fns*(this%n_target_fns+1)/2

         this%ordered_pairs = 0
         this%ordered_pairs_to_indices = 0
         this%last_CT_fn = 0
         ind = 0
         ij = 0
         do i=1,this%number_of_functions
            do j=1,i
               ij = ij + 1

               if (is_continuum(i) .and. is_continuum(j)) cycle !Skip the CC functions for now

               ind = ind + 1
               this%ordered_pairs(1,ij) = ind !index of this pair of functions
               this%ordered_pairs_to_indices(1:2,ind) = (/i,j/)

               !type of this pair of functions: TT = 1, CT = 2, CC = 3
               if (.not.(is_continuum(i)) .and. .not.(is_continuum(j))) then
                  this%ordered_pairs(2,ij) = 1
               elseif (is_continuum(i) .or. is_continuum(j)) then
                  this%ordered_pairs(2,ij) = 2
               endif

               !Find the index of the last CT-type function
               if (is_continuum(i) .or. is_continuum(j)) then
                  this%last_CT_fn = max(this%last_CT_fn,ind)
               endif

            enddo !q
         enddo !p

         this%n_prec_ints = this%last_CT_fn*(this%last_CT_fn+1)/2

         ij = 0
         do i=1,this%number_of_functions
            do j=1,i
               ij = ij + 1

               if (is_continuum(i) .and. is_continuum(j)) then !index only the CC pairs of functions
                  ind = ind + 1
                  this%ordered_pairs(1,ij) = ind !index of this pair of functions
                  this%ordered_pairs(2,ij) = 3
                  this%ordered_pairs_to_indices(1:2,ind) = (/i,j/)
               endif

            enddo !j
         enddo !i

   end subroutine generate_fn_indices

   subroutine print_atomic_orbital_basis_obj(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'print_atomic_orbital_basis_obj', &
                         'The object has not been initialized or not all shells have been added.', 1, 1)
         end if

!         write(stdout,'(/,"--------->","atomic_orbital_basis_obj:print")')
!         write(stdout,'("GTO basis shell data:")')
!         write(stdout,'("shell; index of the basis function starting the shell; nucleus index; GTO center; GTO L; number of primitives; norm of the contracted GTO")')
!         write(stdout,'("exponents of the primitive GTOs; contraction coefficients for the primitive GTOs; norms of the primitive GTOs")')
!
!         do s=1,this%shell
!            p = this%n_prim(s)
!            write(stdout,'(i4,1x,i4,1x,i2,1x,3e25.15,1x,i1,1x,i2,1x,e25.15)') s, this%cgto_bas_ind(s), this%nuc(s), this%RA(1:3,s), this%cgto_l(s), this%n_prim(s), this%cgto_norm(s)
!            do i=1,p
!               write(stdout,'(5X,e25.15,1X,e25.15,1X,e25.15)') this%alp(i,s), this%ccf(i,s), this%prim_norm(i,s)
!            enddo
!         enddo
!         write(stdout,'("<---------","done:atomic_orbital_basis_obj:print")')

   end subroutine print_atomic_orbital_basis_obj

   subroutine shell_pair_one_electron_integrals(this,i,j,integral_options,int_index,integrals,number_of_shell_integrals)
      use const, only: number_of_types_el_ints
      use cgto_integrals, only: GG_shell_integrals
      use bto_integrals_mod, only: BB_shell_integrals
      use bto_gto_integrals_mod, only: BG_shell_integrals
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: i, j
      class(integral_options_obj), intent(in) :: integral_options
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)
      integer, intent(out) :: number_of_shell_integrals

      integer :: ind_Gi, ind_Gj, ind_Bi, ind_Bj

         if (this%shell_descriptor(1,i) .eq. 1) then    !the i-shell is a CGTO shell
            ind_Gi = this%shell_descriptor(2,i)
            if (this%shell_descriptor(1,j) .eq. 1) then     !the j-shell is a CGTO shell
               !print *,'1'
               ind_Gj = this%shell_descriptor(2,j)
               number_of_shell_integrals = this % CGTO_shell_data(ind_Gi) % number_of_functions &
                                         * this % CGTO_shell_data(ind_Gj) % number_of_functions
               call GG_shell_integrals(this % CGTO_shell_data(ind_Gi), this % CGTO_shell_data(ind_Gj), ind_Gi, ind_Gj, &
                                       this % shell_descriptor(4,i), this % shell_descriptor(4,j), &
                                       integral_options % use_spherical_cgto_alg, integral_options % max_property_l, &
                                       this % integral_data % property_center, this % symmetry_data, &
                                       this % integral_data % olap_column, this % integral_data % kei_column, &
                                       this % integral_data % prop_column, this % integral_data % nai_column, &
                                       this % integral_data % one_elham_column, int_index, integrals)
            elseif (this%shell_descriptor(1,j) .eq. 2) then !the j-shell is a BTO shell
               !print *,'2'
               ind_Bj = this%shell_descriptor(2,j)
               number_of_shell_integrals = this % CGTO_shell_data(ind_Gi) % number_of_functions &
                                         * this % BTO_shell_data(ind_Bj) % number_of_functions
               call BG_shell_integrals(this%CGTO_shell_data(ind_Gi),this%BTO_shell_data(ind_Bj),&
                                      &this%shell_descriptor(4,i),this%shell_descriptor(4,j),ind_Gi,&
                                      &this%integral_data%olap_column,this%integral_data%kei_column,this%integral_data%prop_column,&
                                      &this%integral_data%nai_column,this%integral_data%one_elham_column,int_index,integrals)
            else !error
               call xermsg ('atomic_orbital_basis_obj', 'shell_pair_one_electron_integrals', &
                            'The shell type B must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 2, 1)
            endif
         elseif (this%shell_descriptor(1,i) .eq. 2) then !the i-shell is a BTO shell
            ind_Bi = this%shell_descriptor(2,i)
            if (this%shell_descriptor(1,j) .eq. 1) then     !the j-shell is a CGTO shell
               !print *,'3'
               ind_Gj = this%shell_descriptor(2,j)
               number_of_shell_integrals = this % CGTO_shell_data(ind_Gj) % number_of_functions &
                                         * this % BTO_shell_data(ind_Bi) % number_of_functions
               call BG_shell_integrals(this%CGTO_shell_data(ind_Gj),this%BTO_shell_data(ind_Bi),&
                                      &this%shell_descriptor(4,j),this%shell_descriptor(4,i),ind_Gj,&
                                      &this%integral_data%olap_column,this%integral_data%kei_column,this%integral_data%prop_column,&
                                      &this%integral_data%nai_column,this%integral_data%one_elham_column,int_index,integrals)
            elseif (this%shell_descriptor(1,j) .eq. 2) then !the j-shell is a BTO shell
               !print *,'4'
               ind_Bj = this%shell_descriptor(2,j)
               number_of_shell_integrals = this % BTO_shell_data(ind_Bi) % number_of_functions &
                                         * this % BTO_shell_data(ind_Bj) % number_of_functions
               call BB_shell_integrals(this % BTO_shell_data(ind_Bi), this % BTO_shell_data(ind_Bj), &
                                       this % shell_descriptor(4,i), this % shell_descriptor(4,j), &
                                       integral_options % a, integral_options % max_property_l, &
                                       this % integral_data % property_center, this % symmetry_data, &
                                       this % integral_data % olap_column, this % integral_data % kei_column, &
                                       this % integral_data % prop_column, this % integral_data % nai_column, &
                                       this % integral_data % one_elham_column, int_index, integrals)
            else !error
               call xermsg ('atomic_orbital_basis_obj', 'shell_pair_one_electron_integrals', &
                            'The shell type B must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 3, 1)
            endif
         else !error
            call xermsg ('atomic_orbital_basis_obj', 'shell_pair_one_electron_integrals', &
                         'The shell type A must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 1, 1)
         endif

   end subroutine shell_pair_one_electron_integrals

   subroutine generate_shell_indices(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this

      integer :: err, i, j, ij, ind, number_of_pairs

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'generate_shell_indices', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         !The shells are indexed in the order in which they were added to the basis: this order is stored in this%shell_descriptor.
         this%shell_descriptor(4,1) = 1
         do i=2,this%number_of_shells
            if (this%shell_descriptor(1,i-1) .eq. 1) then     !the i-shell is a CGTO shell
               this % shell_descriptor(4,i) = this % shell_descriptor(4,i-1) &
                                            + this % CGTO_shell_data(this % shell_descriptor(2,i-1)) % number_of_functions
            elseif (this%shell_descriptor(1,i-1) .eq. 2) then !the i-shell is a BTO shell
               this % shell_descriptor(4,i) = this % shell_descriptor(4,i-1) &
                                            + this % BTO_shell_data(this % shell_descriptor(2,i-1)) % number_of_functions
            else !error
               call xermsg ('atomic_orbital_basis_obj', 'generate_shell_indices', &
                            'The shell type B must be one of:CGTO_shell_data_obj, BTO_shell_data_obj.', 2, 1)
            endif
         enddo !i

         this%n_target_sh = 0
         this%n_cont_sh = 0

         this%last_CT_sh = 0
         this%n_prec_sh = 0
         this%n_TT_sh_pairs = 0
         !Count the number of target and continuum shells
         do i=1,this%number_of_shells
            if (this%shell_descriptor(3,i) .eq. 1) then
               this%n_cont_sh = this%n_cont_sh + 1
            else
               this%n_target_sh = this%n_target_sh + 1
            endif
         enddo

         this%n_TT_sh_pairs = this%n_target_sh*(this%n_target_sh+1)/2

         !Calculate indices for pairs of shells such that the CC pairs
         !have the largest indices: this ensures optimal indexing of shell
         !quartets in case 2p integrals with only 1p in the continuum integrals are requested and more than one MPI task is used.
         number_of_pairs = this%number_of_shells*(this%number_of_shells+1)/2
         if (allocated(this%ordered_shell_pairs)) deallocate(this%ordered_shell_pairs)
         allocate(this%ordered_shell_pairs(2,number_of_pairs),stat=err)
         if (err /= 0) then
            call xermsg ('atomic_orbital_basis_obj', 'generate_shell_indices', &
                         'Memory allocation 4 failed.',err, 1)
         end if

         this%last_CT_sh = 0
         this%ordered_shell_pairs = 0
         ind = 0
         ij = 0
         do i=1,this%number_of_shells
            do j=1,i
               ij = ij + 1

               if (this%shell_descriptor(3,i) .eq. 1 .and. this%shell_descriptor(3,j) .eq. 1) cycle !Skip the CC shells for now

               ind = ind + 1
               this%ordered_shell_pairs(1,ij) = ind !index of this pair of shells

               !type of this pair of functions: TT = 1, CT = 2, CC = 3
               if (this%shell_descriptor(3,i) .eq. 0 .and. this%shell_descriptor(3,j) .eq. 0) then
                  this%ordered_shell_pairs(2,ij) = 1
               elseif (this%shell_descriptor(3,i) .eq. 1.or. this%shell_descriptor(3,j) .eq. 1) then
                  this%ordered_shell_pairs(2,ij) = 2
               endif

               !Find the index of the last CT-type shell
               if (this%shell_descriptor(3,i) .eq. 1 .or. this%shell_descriptor(3,j) .eq. 1) then
                  this%last_CT_sh = max(this%last_CT_sh,ind)
               endif

            enddo !q
         enddo !p

         this%n_prec_sh = this%last_CT_sh*(this%last_CT_sh+1)/2

         ij = 0
         do i=1,this%number_of_shells
            do j=1,i
               ij = ij + 1

               if (this%shell_descriptor(3,i) .eq. 1 .and. this%shell_descriptor(3,j) .eq. 1) then !index only the CC pairs of shells
                  ind = ind + 1
                  this%ordered_shell_pairs(1,ij) = ind !index of this pair of shells
                  this%ordered_shell_pairs(2,ij) = 3
               endif

            enddo !j
         enddo !i

         if (allocated(this%shell_pair_indices)) deallocate(this%shell_pair_indices)
         if (allocated(this%shell_pair_type)) deallocate(this%shell_pair_type)
         allocate(this%shell_pair_indices(4,number_of_pairs),this%shell_pair_type(number_of_pairs),stat=err)
         if (err /= 0) then
            call mpi_xermsg ('atomic_orbital_basis_obj', 'generate_shell_indices', &
                             'Memory allocation 2 failed.', err, 1)
         end if

         !Note that the same loop structure for pairs must be maintained in the routine GGGG_initialize and wherever else loop over pairs of shells is encountered.
         !This is especially important for the calculation of the tail integrals where the target shells have the same angular momentum.
         ij = 0
         do i=1,this%number_of_shells
            do j=1,i
               ij = ij+1

               !The type of the pair of functions is computed as an index of the
               !combination of pairs in a square matrix with two rows (1=G,2=B).
               this%shell_pair_type(ij) = this%shell_descriptor(1,i) + 2*(this%shell_descriptor(1,j)-1)

               this%shell_pair_indices(1,ij) = this%shell_descriptor(2,i)
               this%shell_pair_indices(2,ij) = this%shell_descriptor(2,j)
               this%shell_pair_indices(3,ij) = i
               this%shell_pair_indices(4,ij) = j

            enddo !j
         enddo !i

   end subroutine generate_shell_indices

   function get_number_of_functions_in_shell(this,i)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: i
      integer :: get_number_of_functions_in_shell

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_number_of_functions_in_shell', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this % number_of_shells) then
            call xermsg ('atomic_orbital_basis_obj', 'get_number_of_functions_in_shell', &
                         'On input i was out of range.', 2, 1)
         end if

         if (this%shell_descriptor(1,i) .eq. 1) then     !the i-shell is a CGTO shell
            get_number_of_functions_in_shell = this%CGTO_shell_data(this%shell_descriptor(2,i))%number_of_functions
         elseif (this%shell_descriptor(1,i) .eq. 2) then !the i-shell is a BTO shell
            get_number_of_functions_in_shell = this%BTO_shell_data(this%shell_descriptor(2,i))%number_of_functions
         else !error
            call xermsg ('atomic_orbital_basis_obj', 'get_number_of_functions_in_shell', &
                         'The shell type B must be one of:CGTO_shell_data_obj, BTO_shell_data_obj. &
                         &Error in this%shell_descriptor.', 3, 1)
         endif

   end function get_number_of_functions_in_shell

   subroutine one_electron_integrals(this,integral_storage,integral_options)
      use const
      use bto_integrals_mod, only: BB_initialize
      use bto_gto_integrals_mod, only: BG_initialize, max_l_BG
      use cgto_integrals, only: GG_initialize
      use omp_lib
      implicit none
      class(atomic_orbital_basis_obj) :: this
      class(integral_options_obj), intent(in) :: integral_options
      class(integral_storage_obj), intent(inout) :: integral_storage

      !Input/output of the integral routines:
      integer :: i, j, k, ind, integral_type, number_of_shell_integrals,number_of_types, err
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: shell_integrals(:,:)

      !Input/output of the calculated integrals:
      integer :: number_of_integrals, lunit, first_record, number_of_zero_ints, current_pos, a, last_record
      type(p2d_array_obj), target :: temp_p_integral_array
      type(p2d_array_obj), pointer :: p_integral_array
      integer, parameter :: number_of_blocks = 0
      character(len=line_len) :: column_descriptor(number_of_types_el_ints), header

      !Auxiliary:
      real(kind=wp) :: start_t, end_t
      integer :: max_l_cgto

         start_t = omp_get_wtime()

         call mpi_mod_barrier(err)

         write(stdout,'("--------->","atomic_orbital_basis_obj:one_electron_integrals")')

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         call integral_options%print

         err = integral_options%check()
         if (err /= 0) then
            call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                         'integral_options%check failed with an error message; see integral_options_obj%check.', err, 1)
         end if

         if (integral_options % use_spherical_cgto_alg) then
            call xermsg('atomic_orbital_basis_obj', 'one_electron_integrals', &
                        'Not enabled yet for nuclear attraction and property integrals.', -2, 0)
         end if

         !Determine the total number of pairs of functions for which the integrals will be calculated (this is equal to the number of rows in the final array of integrals).
         number_of_integrals = this%number_of_functions*(this%number_of_functions+1)/2

         if (number_of_integrals == 0) then
            call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                         'The estimated number of integrals to be generated is zero.', 3, 1)
         end if

         !Determine the types of integrals to calculate (the number of types of integrals is equal to the number of columns of the final array of integrals).
         number_of_types = 0
         if (integral_options%calculate_overlap_ints) then
            number_of_types = number_of_types + 1
            column_descriptor(number_of_types) = overlap_ints
         endif
         if (integral_options%calculate_kinetic_energy_ints) then
            number_of_types = number_of_types + 1
            column_descriptor(number_of_types) = kinetic_ints
         endif

         if (integral_options % calculate_overlap_ints .neqv. integral_options % calculate_kinetic_energy_ints) then
            call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                         'If one of overlap integrals/kinetic energy integrals are selected, &
                         &the other one must be selected too.', 4, 1)
         end if

         if (integral_options%calculate_nuclear_attraction_ints) then
            number_of_types = number_of_types + 1
            column_descriptor(number_of_types) = nuc_rep_att_ints
         endif
         if (integral_options%calculate_one_el_hamiltonian_ints) then
            number_of_types = number_of_types + 1
            column_descriptor(number_of_types) = one_elham
            if (integral_options % calculate_kinetic_energy_ints .neqv. integral_options % calculate_nuclear_attraction_ints) then
               call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                            'If the 1-electron Hamiltonian integrals are selected the kinetic energy and nuclear att. &
                            &integrals must be selected too.', 5, 1)
            end if
         endif
         !Transfer all labels to the this%integral_data%column_descriptor data structure. The property integrals are always stored in the last columns of the final integral array.
         if (integral_options%calculate_property_ints) then
            i = (integral_options%max_property_l+1)**2 !number of extra columns required to store all property integrals

            allocate(this%integral_data%column_descriptor(number_of_types+i),stat=err)
            if (err /= 0) then
               call mpi_xermsg('atomic_orbital_basis_obj', 'one_electron_integrals', &
                               'Memory allocation 1 failed.', err, 1)
            end if

            this%integral_data%column_descriptor(1:number_of_types) = column_descriptor(1:number_of_types)
            do j=number_of_types+1,number_of_types+i
               this%integral_data%column_descriptor(j) = property_ints
            enddo !j
            number_of_types = number_of_types+i
         else
            allocate(this%integral_data%column_descriptor(number_of_types),stat=err)
            if (err /= 0) then
               call mpi_xermsg('atomic_orbital_basis_obj', 'one_electron_integrals', &
                               'Memory allocation 2 failed.', err, 1)
            end if

            this%integral_data%column_descriptor(1:number_of_types) = column_descriptor(1:number_of_types)
         endif

         this%integral_data%number_of_types = number_of_types

         write(stdout,'("List of 1-electron integrals to calculate and their column indices in the final array of integrals: ")')
         do j=1,number_of_types
            write(stdout,'(i0,": ",a207)') j, adjustl(this%integral_data%column_descriptor(j))
         enddo !i

         if (number_of_types == 0) then
            call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                         'No integrals to calculate have been selected.', 6, 1)
         end if

         header = integral_storage%contruct_header_string(this%get_basis_name(),one_electron_ints)

         !allocate the output arrays if we request the output to be stored in memory
         if (integral_storage%in_memory()) then
            integral_storage%data_header%name = header

            p_integral_array => integral_storage%integral
            !we allocate space for a non-indexed (that is purely local) array with number_of_types columns and number_of_integrals rows
            err = p_integral_array%init(number_of_integrals,number_of_types,number_of_blocks,this%integral_data%column_descriptor)
            if (err /= 0) then
               call mpi_xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                                'Array initialization 1 failed; see p2d_array_obj%init.', err, 1)
            end if
         endif

         !if we request the output to be stored on disk then start a new record on the data file that will contain the integrals
         if (integral_storage%on_disk()) then

            !temporary storage for the integrals
            !we allocate space for a non-indexed (that is purely local) array with number_of_types columns and number_of_integrals rows
            err = temp_p_integral_array % init(number_of_integrals, &
                                               number_of_types,     &
                                               number_of_blocks,    &
                                               this % integral_data % column_descriptor)
            if (err /= 0) then
               call mpi_xermsg('atomic_orbital_basis_obj', 'one_electron_integrals', &
                               'Array initialization 2 failed; see p2d_array_obj%init.', err, 1)
            end if
            p_integral_array => temp_p_integral_array

            lunit = integral_storage%integral_file%get_unit_no() !unit that is associated to the file opened
            first_record = integral_storage%integral_file%start_record(header) !position, within the data file, of the first record available for the integral data
         endif
!
!------- PREPARE THE DATA NEEDED FOR THE CALCULATION:
!
         ! Determine the column indices of the integrals to calculate: note that the routines from shell_pair_one_electron_integrals
         ! must place the respective integrals in the same columns!!!
         this % integral_data % olap_column = 0
         this % integral_data % kei_column  = 0
         this % integral_data % prop_column = 0
         this % integral_data % nai_column  = 0
         this % integral_data % one_elham_column = 0

         do i=1,this%integral_data%number_of_types
            select case (this%integral_data%column_descriptor(i))
            case (overlap_ints) !Overlap integrals:
               this%integral_data%olap_column = i
            case (kinetic_ints) !Kinetic energy integrals:
               this%integral_data%kei_column = i
            case (nuc_rep_att_ints) !Nuclear attraction integrals:
               this%integral_data%nai_column = i
            case (one_elham) !One electron Hamiltonian integrals:
               this%integral_data%one_elham_column = i
            case (property_ints)
               this%integral_data%prop_column = i
               exit !we stop as soon as we encounter the first property integrals header since all types of properties must be calculated at once.
            end select
         enddo !i

         !todo check that the last BTO doesn't leak outside of the sphere r=integral_options%a.
         !Generate the value of the starting index for each shell in the basis: result in this%shell_descriptor(4,:)
         call this%generate_shell_indices

         !Initialize the integrals module for CGTO/CGTO integrals.
         if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
            call GG_initialize(this%CGTO_shell_data,integral_options%a)
         endif

         if (integral_options%a > 0.0_cfp) then
            call this%normalize_continuum(integral_options%a)
         endif

         if (allocated(this % BTO_shell_data) .and. size(this % BTO_shell_data) > 0) then
            this%integral_data%max_bspline_l = maxval(this%BTO_shell_data(:)%l)
            this%integral_data%first_bspline_index = minval(this%BTO_shell_data(:)%bspline_index)

            !Initialize quadratures needed to compute the BTO/BTO integrals:
            call BB_initialize(this % BTO_shell_data(1) % bspline_grid, &
                               this % integral_data % max_bspline_l, &
                               integral_options % max_property_l, &
                               this % symmetry_data)

            !Initialize quadratures needed to compute the BTO/CGTO integrals:
            if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
               max_l_cgto = maxval(this%CGTO_shell_data(:)%l)
               call BG_initialize(integral_options % max_l_legendre_1el, &
                                  this % BTO_shell_data(1) % bspline_grid, &
                                  this % integral_data % first_bspline_index, &
                                  this % integral_data % max_bspline_l, &
                                  integral_options % max_property_l, &
                                  max_l_cgto, &
                                  this % symmetry_data % nucleus, &
                                  integral_options % delta_r1, &
                                  this % integral_data % nai_column, &
                                  integral_options % mixed_ints_method)
            endif

         endif

         !todo it should be made possible to change this via integral_options
         this%integral_data%property_center(1:3) = 0.0_cfp
!
!-------
!
         !Allocate enough space for int_index,shell_integrals:
         number_of_shell_integrals = 0
         do i=1,this%number_of_shells
            do j=1,i
               number_of_shell_integrals = max(number_of_shell_integrals,this%shell_descriptor(5,i)*this%shell_descriptor(5,j))
            enddo !j
         enddo !i

         allocate(int_index(2, number_of_shell_integrals), &
                  shell_integrals(number_of_shell_integrals, this % integral_data % number_of_types), &
                  stat = err)

         if (err .ne. 0) call xermsg('atomic_orbital_basis_obj', 'one_electron_integrals', 'Memory allocation 3 failed.',err,1)
         int_index = 0; shell_integrals = 0.0_cfp

         !Calculate and store the integrals for all pairs of shells in the basis:
         number_of_zero_ints = 0

         !We loop over the upper triangular part of the shell pair matrix since this is the optimal order for the algorithm for the mixed BG integrals.
         do i=1,this%number_of_shells
            do j=i,this%number_of_shells

               !Calculate all types of 1-electron integrals for this shell-pair:
               call this%shell_pair_one_electron_integrals(i,j,integral_options,int_index,shell_integrals,number_of_shell_integrals)

               !Index and store the integrals:
               do integral_type=1,number_of_types
                  do k=1,number_of_shell_integrals

                     !ignore integrals smaller than tol
                     if (abs(shell_integrals(k,integral_type)) < integral_options%tol) then
                        number_of_zero_ints = number_of_zero_ints + 1
                        cycle
                     endif

                     !Note that this must be consistent with the indexing scheme in the routine integral_index.
                     ind = int_index(1,k)*(int_index(1,k)-1)/2 + int_index(2,k)
                     p_integral_array%a(ind,integral_type) = shell_integrals(k,integral_type)
                  enddo !k
               enddo !integral_type

            enddo !j
         enddo !i

         if (allocated(this % BTO_shell_data) .and. size(this % BTO_shell_data) > 0) then
            write(stdout,'("Maximum L needed to converge the mixed nuclear attraction integrals: ",i4)') max_l_BG
         end if

         !If requested print the non-zero integrals
         if (integral_options%print_integrals) then
            call p_integral_array%print(.true.)
         endif

         !Dump all integrals to disk and close the record.
         !Only master writes to the file but every process gets value of: current_pos (integral_options%write) and
         !last_record (p_integral_array%write), i.e. the positions in the file where
         !the master ended up following the calls to these two routines.
         if (integral_storage%on_disk()) then

            write(stdout,'("Saving integrals to disk...")')

            !The first record are the integral options.
            call integral_options%write(lunit,first_record,current_pos)

            !The second record are the ordered integrals.
            a = master
            call p_integral_array%write(lunit,current_pos,last_record,a)

            !Every process closes the record so that they all keep identical header information.
            call integral_storage%integral_file%close_record(header,first_record,last_record)

            err = p_integral_array%final()
            if (err /= 0) then
               call xermsg ('atomic_orbital_basis_obj', 'one_electron_integrals', &
                            'Deallocation of the temporary integral array failed.', 7, 1)
            end if

            nullify(p_integral_array)
            !FROM NOW ON p_integral_array is nullified

            write(stdout,'("...done")')

         endif

         write(stdout,'("<---------","atomic_orbital_basis_obj:one_electron_integrals")')

         end_t = omp_get_wtime()
         write(stdout,'("One_electron_integrals took [s]: ",f25.15)') end_t-start_t

         call mpi_mod_barrier(err)

   end subroutine one_electron_integrals

   subroutine two_electron_integrals(this,integral_storage,integral_options)
      use const
      use special_functions, only: ipair
      use cgto_integrals, only: GGGG_initialize, GGGG_final
      use bto_gto_integrals_mod, only: BBGG_shell_integrals, BG_mixed_2el_initialize, BGGG_shell_integrals, &
                                       BGBG_shell_integrals, max_l_BGGG, max_l_BGBG, BG_final
      use gto_routines, only: index_2el_drv, index_1p_continuum
      use sort, only: sort_int_float
      use omp_lib
      implicit none
      class(atomic_orbital_basis_obj) :: this
      class(integral_options_obj), intent(in) :: integral_options
      class(integral_storage_obj), intent(inout) :: integral_storage

      !Input/output of the integral routines:
      integer :: i, j, k, l, ind, integral_type, number_of_shell_integrals,number_of_types, err
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: shell_integrals(:,:)

      !Input/output of the calculated integrals:
      integer :: number_of_integrals, lunit, first_record, current_pos, a, last_record, max_functions, max_functions_continuum
      type(p2d_array_obj), target :: temp_p_integral_array
      type(p2d_array_obj), pointer :: p_integral_array
      integer :: number_of_blocks = 0
      character(len=line_len) :: column_descriptor(number_of_types_el_ints), header

      !Book-keeping:
      logical, allocatable :: is_CC(:), is_TT(:), is_CT(:)
      integer :: ij, kl, number_of_pairs, number_of_cgto_pairs, dim_shell_integrals, thread_id, n_threads
      logical :: keep_ab_cd_order, storage_in_memory, is_CCTT

      !Auxiliary:
      real(kind=wp) :: start_t, end_t, t1, t2
      integer, allocatable :: quartet_offset(:)
      integer :: pq, rs, rank, no_quartets, offset, quartet_ind, i_abs, j_abs, k_abs, l_abs, integral_class
      integer :: kl_min, kl_max, ij_min, ij_max, number_of_classes, indexing_method

         start_t = omp_get_wtime()

         call mpi_mod_barrier(err)

         write(stdout,'("--------->","atomic_orbital_basis_obj:two_electron_integrals")')

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         call integral_options%print

         err = integral_options%check()
         if (err /= 0) then
            call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                         'integral_options%check failed with an error message; see integral_options_obj%check.', err, 1)
         end if

         if (.not. integral_options % calculate_two_el_ints) then
            call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                         'Inconsistency with integral_options: calculate_two_el_ints = .false.', 2, 1)
         end if

         !Determine the types of integrals to calculate (the number of types of integrals is equal to the number of columns of the final array of integrals).
         number_of_types = 1 !we implement only one type of 2-electron integral
         column_descriptor(1) = two_el_ints

         !Transfer all labels to the this%integral_data%column_descriptor data structure.
         if (allocated(this%integral_data%column_descriptor)) deallocate(this%integral_data%column_descriptor)
         allocate(this%integral_data%column_descriptor(number_of_types),stat=err)
         if (err /= 0) then
            call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                             'Memory allocation 1 failed.', err, 1)
         end if

         this%integral_data%column_descriptor(1:number_of_types) = column_descriptor(1:number_of_types)
         this%integral_data%number_of_types = number_of_types

         write(stdout,'("List of 2-electron integrals to calculate and their column indices in the final array of integrals: ")')
         do j=1,number_of_types
            write(stdout,'(i0,": ",a207)') j, adjustl(this%integral_data%column_descriptor(j))
         enddo !i
!
!------- PREPARE THE DATA NEEDED FOR THE CALCULATION:
!
         !Generate the value of the starting index for each shell in the basis: result in this%shell_descriptor(4,:)
         call this%generate_shell_indices

         !Generate the logical arrays telling me what type each shell-pair is:
         number_of_pairs = this%number_of_shells*(this%number_of_shells+1)/2

         allocate(is_CC(number_of_pairs),is_CT(number_of_pairs),is_TT(number_of_pairs),stat=err)
         if (err /= 0) then
            call mpi_xermsg('atomic_orbital_basis_obj', 'two_electron_integrals', &
                            'Memory allocation 2 failed.', err, 1)
         end if
         is_CC = .false.; is_CT = .false.; is_TT = .false.

         !Note that the same loop structure for pairs must be maintained in the routine GGGG_initialize and wherever else loop over pairs of shells is encountered.
         !This is especially important for the calculation of the tail integrals where the target shells have the same angular momentum.
         ij = 0
         do i=1,this%number_of_shells
            do j=1,i
               ij = ij+1

               if (this%shell_descriptor(3,i) .eq. 1 .and. this%shell_descriptor(3,j) .eq. 1) is_CC(ij) = .true.
               if (this%shell_descriptor(3,i) .eq. 0 .and. this%shell_descriptor(3,j) .eq. 1) is_CT(ij) = .true.
               if (this%shell_descriptor(3,i) .eq. 1 .and. this%shell_descriptor(3,j) .eq. 0) is_CT(ij) = .true.
               if (this%shell_descriptor(3,i) .eq. 0 .and. this%shell_descriptor(3,j) .eq. 0) is_TT(ij) = .true.

            enddo !j
         enddo !i

         !Decide how the 2-electron integrals within each quartet of shells will be ordered
         if (nprocs > 1) then
            !Integrals for each quartet of shells are ordered so that the 1st column corresponds to functions in the shell with the largest starting index,
            !the 2nd column corresponds to functions in the shell with the 2nd largest starting index, etc.
            keep_ab_cd_order = .true.
            indexing_method = 2
         else
            !Natural ordering as output by the integral algorithms
            keep_ab_cd_order = .false.
            indexing_method = 1
         endif

         if (integral_options%a > 0.0_cfp) then
            call this%normalize_continuum(integral_options%a)
         endif

         number_of_cgto_pairs = size(this%CGTO_shell_data)
         number_of_cgto_pairs = number_of_cgto_pairs*(number_of_cgto_pairs+1)/2
         if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
            !Initialize all data needed to compute the 2-electron integrals (and tails) over GTOs-only
            call GGGG_initialize (this % CGTO_shell_data, this % shell_descriptor(4,:), integral_options % tol, &
                                  integral_options % a, keep_ab_cd_order, integral_options % two_p_continuum, indexing_method)
         endif

         if (allocated(this % BTO_shell_data) .and. size(this % BTO_shell_data) > 0) then
            this%integral_data%max_bspline_l = maxval(this%BTO_shell_data(:)%l)
            this%integral_data%first_bspline_index = minval(this%BTO_shell_data(:)%bspline_index)

            if (integral_options%two_p_continuum) then
               call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                                '2-el ints for 2p in the continuum involving B-splines have not been implemented.', 2, 1)
            end if

            !Initialize quadratures needed to compute the BTO/CGTO integrals:
            if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
               if (nprocs > 1) then
                  call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                                   'Parallel mixed 2-el ints not implemented: indexing needs to be amended.', 3, 1)
               end if
               if (integral_options%mixed_ints_method .le. 0 .or. integral_options%mixed_ints_method > 3) then
                  print *,integral_options%mixed_ints_method
                  call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                                   'Select a valid method for the calculation of the mixed BTO/CGTO integrals.', 4, 1)
               endif
               call BG_mixed_2el_initialize(integral_options % max_l_legendre_2el, &
                                            this % BTO_shell_data(1) % bspline_grid, &
                                            this % integral_data % first_bspline_index, &
                                            this % integral_data % max_bspline_l, &
                                            integral_options % delta_r1, &
                                            this % CGTO_shell_data, &
                                            keep_ab_cd_order, &
                                            integral_options % mixed_ints_method, &
                                            integral_options % scratch_directory, &
                                            indexing_method)
            endif
         endif

         !Set the column index of the 2-el integrals to calculate: note that the routines from shell_quartet_two_electron_integrals
         !must place the respective integrals in the same column!!!
         this%integral_data%two_el_column = 1
!
!------- Allocate enough space for int_index,shell_integrals:
         if (integral_options%two_p_continuum) then
            max_functions = 0
            do i=1,this%number_of_shells
               max_functions = max(max_functions,this%shell_descriptor(5,i))
            enddo !i
            number_of_shell_integrals = max_functions**4
         else !discard the (CC|CC) and (CC|CT) combinations of integrals from the counting

            !What is the largest number of functions in the target and continuum shells?
            max_functions = 0 !target
            max_functions_continuum = 0 !continuum

            do i=1,this%number_of_shells
               if (this%shell_descriptor(3,i) .eq. 1) then
                  max_functions_continuum = max(max_functions_continuum,this%shell_descriptor(5,i))
               else
                  max_functions = max(max_functions,this%shell_descriptor(5,i))
               endif
            enddo !i
            !Find out what is the largest number of integrals in all types of quartets of shells for which the integrals will be calculated.
            number_of_shell_integrals = max_functions**4 !(TT|TT)
            number_of_shell_integrals = max(number_of_shell_integrals,max_functions**3*max_functions_continuum) !(TT|TC)
            number_of_shell_integrals = max(number_of_shell_integrals,max_functions**2*max_functions_continuum**2) !(TC|TC) and (TT|CC)
         endif

         !Save the number of shell integrals so that each thread can later allocate its own copy of int_index and shell_integrals
         dim_shell_integrals = number_of_shell_integrals
!
!------- WORK DISTRIBUTION: perform the initial sweep of all unique quartets of shells to determine which shell-pairs this process must calculate integrals for.
         !For each ij the integrals are cyclically distributed among the processes.
         if (nprocs > 1) then

            if (.not.(integral_options%two_p_continuum)) then
               if (this%n_cont_sh .eq. 0) then !the basis contains only target shells
                  no_quartets = this%n_TT_sh_pairs*(this%n_TT_sh_pairs+1)/2
               else !the basis contains continuum shells
                  no_quartets = this%n_prec_sh !all TTTT, CTTT and CTCT quartets
                  no_quartets = no_quartets + this%n_cont_sh*(this%n_cont_sh+1)/2*this%n_TT_sh_pairs !add the CCTT class
               endif
            else
               no_quartets = number_of_pairs*(number_of_pairs+1)/2 !number of unique pairs of pairs of shells including all unique quartets of continuum shells
            endif

            number_of_blocks = no_quartets
            allocate(quartet_offset(no_quartets),stat=err)
            if (err /= 0) then
               call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                            'Memory allocation 4 failed.', err, 1)
            end if

            t1 = omp_get_wtime()
            no_quartets = 0
            number_of_integrals = 0
            !Note that the order of the ij,kl loops must be the same as in the actual calculation below!!!
            !WARNING: only the kl=1,ij order gives the correct results. For some reason the order kl=ij,number_of_pairs gives incorrect results.
            !In serial (nprocs .eq. 1) everything works fine no matter what order is used.
            offset = 0
            quartet_offset = -1
            do ij=1,number_of_pairs
               i = this%shell_pair_indices(3,ij)
               j = this%shell_pair_indices(4,ij)
               do kl=1,ij !ij,number_of_pairs
                  !quartet_ind = quartet_ind + 1
                  k = this%shell_pair_indices(3,kl)
                  l = this%shell_pair_indices(4,kl)
                  !In case we have only 1p in the continuum we don't need 2p integrals of the type [CC|CC] and [CC|CT].
                  is_CCTT = .false.
                  if (.not.(integral_options%two_p_continuum)) then
                     if (is_CC(ij) .and. .not.(is_TT(kl))) cycle
                     if (is_CC(kl) .and. .not.(is_TT(ij))) cycle
                     if (is_CC(ij) .and. is_TT(kl)) is_CCTT = .true.
                     if (is_TT(ij) .and. is_CC(kl)) is_CCTT = .true.
                  endif
                  i_abs = this%shell_pair_indices(3,ij)
                  j_abs = this%shell_pair_indices(4,ij)
                  k_abs = this%shell_pair_indices(3,kl)
                  l_abs = this%shell_pair_indices(4,kl)
                  quartet_ind = index_1p_continuum(this % ordered_shell_pairs, &
                                                   i_abs, j_abs, k_abs, l_abs, is_CCTT, &
                                                   this % last_CT_sh, &
                                                   this % n_prec_sh, &
                                                   this % n_TT_sh_pairs)

                  rank = mod(quartet_ind,int(nprocs)) != cyclic redistribution of work
                  if (rank .eq. myrank) then !this combination of pairs of shells is the one for me
                     no_quartets = no_quartets + 1
                     !total number of integrals within this shell (ab|cd)
                     number_of_shell_integrals = this % shell_descriptor(5,i) &
                                               * this % shell_descriptor(5,j) &
                                               * this % shell_descriptor(5,k) &
                                               * this % shell_descriptor(5,l)
                     quartet_offset(quartet_ind) = offset !offset+1 is the position in p_integral_array%(:,two_el_column) where integrals for the shell quartet [ij|kl] start.
                     offset = offset + number_of_shell_integrals
                     number_of_integrals = number_of_integrals + number_of_shell_integrals
                  endif
               enddo !kl
            enddo !ij

            write(stdout,'("Rank ",i0," is processing ",i0," quartets of shells corresponding to ",i0," integrals.")') &
                myrank, no_quartets, number_of_integrals

            t2 = omp_get_wtime()
            write(stdout,'("Work redistribution took: ",f8.3," [s]")') t2-t1

         else !SERIAL

            if (.not.(integral_options%two_p_continuum)) then
               if (this%n_cont_fns .eq. 0) then !the basis contains only target functions
                  number_of_integrals = this%n_TT_pairs*(this%n_TT_pairs+1)/2
               else !the basis contains continuum functions
                  number_of_integrals = this%n_prec_ints !all TTTT, CTTT and CTCT pairs
                  number_of_integrals = number_of_integrals + this%n_cont_fns*(this%n_cont_fns+1)/2*this%n_TT_pairs !add the CCTT class
               endif
            else
               number_of_integrals = this%number_of_functions*(this%number_of_functions+1)/2 !number of unique pairs of functions
               number_of_integrals = number_of_integrals*(number_of_integrals+1)/2 !number of unique pairs of pairs of functions including all unique quartets of continuum functions
            endif

            !todo put this into parallel_array?
            !For serial run there is no need to maintain the block indexing structure so we must set number_of_blocks = 0.
            number_of_blocks = 0

         endif !nprocs > 1
         !END: WORK DISTRIBUTION
!
!------- Allocate space for all integrals to be calculated by this process.
         !Name of the integrals to be stored in memory/disk.
         header = integral_storage%contruct_header_string(this%get_basis_name(),two_el_ints)

         !allocate the output arrays if we request the output to be stored in memory
         if (integral_storage%in_memory()) then
            storage_in_memory = .true.
            integral_storage%data_header%name = header

            p_integral_array => integral_storage%integral
            !we allocate space for a non-indexed (that is purely local) array with number_of_types columns and number_of_integrals rows
            err = p_integral_array%init(number_of_integrals,number_of_types,number_of_blocks,this%integral_data%column_descriptor)
            if (err /= 0) then
               call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                                'Array initialization 1 failed; see p2d_array_obj%init.', err, 1)
            end if
         endif

         !if we request the output to be stored on disk then start a new record on the data file that will contain the integrals
         if (integral_storage%on_disk()) then
            storage_in_memory = .false.
            !temporary storage for the integrals
            !we allocate space for a non-indexed (that is purely local) array with number_of_types columns and number_of_integrals rows
            err = temp_p_integral_array % init(number_of_integrals, number_of_types, number_of_blocks, &
                                               this % integral_data % column_descriptor)
            if (err /= 0) then
               call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                                'Array initialization 2 failed; see p2d_array_obj%init.', err, 1)
            end if
            p_integral_array => temp_p_integral_array

            lunit = integral_storage%integral_file%get_unit_no() !unit that is associated to the file opened
            first_record = integral_storage%integral_file%start_record(header) !position, within the data file, of the first record available for the integral data
         endif

         if (number_of_blocks > 0) then
            call p_integral_array%set_block_offset(quartet_offset)
            deallocate(quartet_offset)
         endif
!
!------- Calculate and store the CGTO/CGTO and CGTO/BTO integrals for the (possibly subset) of the unique quartets of shells in the basis:
!
         if (allocated(this % BTO_shell_data) .and. size(this % BTO_shell_data) > 0) then
            number_of_classes = 4
            !we split the integrals into several classes:
            !GGGG: type = 1
            !BGGG: type = 2
            !BGBG: type = 3
            !BBGG: type = 4
         else
            number_of_classes = 1
            !we have only one class of integrals to evaluate: GGGG
         endif

         do integral_class=1,number_of_classes

            t1 = omp_get_wtime()

            select case(integral_class)
               case (1) !GGGG
                  write(stdout,'(/,10X,"Calculating GGGG integrals...")')
               case (2) !BGGG
                  write(stdout,'(/,10X,"Calculating BGGG integrals...")')
               case (3) !BGBG
                  write(stdout,'(/,10X,"Calculating BGBG integrals...")')
               case (4) !BBGG
                  write(stdout,'(/,10X,"Calculating BBGG integrals...")')
               case default
                  call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', 'This class should not arise.',5,1)
            end select

            ij_min = 1; ij_max = number_of_pairs

            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP & PRIVATE(ij,kl,kl_min,kl_max,quartet_ind,is_CCTT,integral_type,k,pq,rs,i,ind,shell_integrals,int_index, &
            !$OMP &         number_of_shell_integrals,thread_id,n_threads,i_abs,j_abs,k_abs,l_abs)

            thread_id = omp_get_thread_num()
            n_threads = omp_get_num_threads()

            if (omp_in_parallel()) then
               !$OMP SINGLE
               write(stdout,'("Threading over quartets of shells using: ",i2," threads.")') n_threads
               !$OMP END SINGLE
            endif

            !Each thread allocates its own copy of int_index, shell_integrals:
            if (.not.(allocated(int_index))) then
               allocate(int_index(4,dim_shell_integrals),stat=err)
               if (err /= 0) then
                  call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                               'Memory allocation 3a failed.', err,1)
               end if
               int_index = 0
            endif

            if (.not.(allocated(shell_integrals))) then
               allocate(shell_integrals(dim_shell_integrals,this%integral_data%number_of_types),stat=err)
               if (err /= 0) then
                  call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                               'Memory allocation 3b failed.', err, 1)
               end if
               shell_integrals = 0.0_cfp
            endif

            !$OMP DO SCHEDULE(DYNAMIC)
            do ij=ij_min,ij_max

               select case(integral_class)
                  case (1) !GGGG
                     if (this%shell_pair_type(ij) .ne. 1) cycle
                     kl_min = 1
                     kl_max = ij
                  case (2) !BGGG
                     if (this%shell_pair_type(ij) .ne. 1) cycle
                     kl_min = number_of_cgto_pairs+1
                     kl_max = number_of_pairs
                  case (3) !BGBG
                     if (this%shell_pair_type(ij) .ne. 2) cycle
                     kl_min = number_of_cgto_pairs+1
                     kl_max = number_of_pairs
                  case (4) !BBGG
                     if (this%shell_pair_type(ij) .ne. 1) cycle
                     kl_min = number_of_cgto_pairs+1
                     kl_max = number_of_pairs
                  case default
                     call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', 'This class should not arise.',6,1)
               end select

!               kl_min = 1; kl_max = ij

               !WARNING: only the kl=1,ij order gives the correct results. For some reason the order kl=ij,number_of_pairs gives incorrect results.
               !In serial (nprocs .eq. 1) everything works fine no matter what order is used.
               do kl=kl_min,kl_max

                  select case(integral_class)
                     case (1) !GGGG
                        if (this%shell_pair_type(kl) .ne. 1) cycle
                     case (2) !BGGG
                        if (this%shell_pair_type(kl) .ne. 2) cycle
                     case (3) !BGBG
                        if (this%shell_pair_type(kl) .ne. 2) cycle
                     case (4) !BBGG
                        if (this%shell_pair_type(kl) .ne. 4) cycle
                     case default
                        call mpi_xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', 'This class should not arise.',6,1)
                  end select

                  !In case we have only 1p in the continuum we don't need 2p integrals of the type [CC|CC] and [CC|CT].
                  is_CCTT = .false.
                  if (.not.(integral_options%two_p_continuum)) then
                     if (is_CC(ij) .and. .not.(is_TT(kl))) cycle
                     if (is_CC(kl) .and. .not.(is_TT(ij))) cycle
                     if (is_CC(ij) .and. is_TT(kl)) is_CCTT = .true.
                     if (is_TT(ij) .and. is_CC(kl)) is_CCTT = .true.
                  endif

                  if (nprocs > 1) then
                     !MPI parallel calculation: calculate integrals only for those quartets of shells which have been marked for me.
                     !Note that the if statements resulting in 'cycle' must be placed after the line above which increments quartet_ind.

                     i_abs = this%shell_pair_indices(3,ij)
                     j_abs = this%shell_pair_indices(4,ij)
                     k_abs = this%shell_pair_indices(3,kl)
                     l_abs = this%shell_pair_indices(4,kl)
                     !It is assumed in index_1p_continuum that i_abs.ge.j_abs, k_abs.ge.l_abs.
                     quartet_ind = index_1p_continuum(this % ordered_shell_pairs, &
                                                      i_abs, j_abs, k_abs, l_abs, is_CCTT, &
                                                      this % last_CT_sh, &
                                                      this % n_prec_sh, &
                                                      this % n_TT_sh_pairs)
                     if (p_integral_array%block_offset(quartet_ind) .eq. -1) cycle
                  endif

                  !Calculate all types of 2-electron integrals for this shell-quartet (this includes tail subtraction where needed):
                  call this % shell_quartet_two_electron_integrals(ij, kl, integral_options, int_index, &
                                                                   shell_integrals, number_of_shell_integrals)

                  !Order and store the calculated integrals:
                  !A) Serial mode:
                  !integrals for each quartet of shells are ordered in the way
                  !which is most advantageous for the integral algorithms, i.e.
                  !no post-ordering is applied. The indexing below assumes that
                  !the indices 1,2,3,4 corresponding to the functions in the
                  !shells a,b,c,d are ordered so that the triangular index of
                  !the AB pair is .ge. index of the CD pair. This implies that
                  !the indices ind_a,ind_b,ind_c,ind_d of the functions in
                  !(ab|cd) must be ordered: ind_a.ge.ind_b, ind_c.ge.ind_d, ind_a.ge.ind_c.
                  !B) MPI parallel:
                  !for each quartet of shells the ``columns" 1,2,3,4 of the linear array shell_integrals(:,integral_type) are permuted according to the starting index of the functions in each shell so that
                  !the starting indices a,b,c,d corresponding to the columns 1,2,3,4 are a.ge.b,c.ge.d,a.ge.c. In other words for each quartet of shells there is a unique order in which the integrals
                  !are stored. This method allows me to pick out the integrals during the transformation without having to store any additional information except the list of quartets of shells that each
                  !MPI task has calculated the integrals for. The index of the quartet of shells is calculated using indexing with or without 2p in the continuum just like index of the actual
                  !integrals in the serial case.
                  !todo skip integrals smaller than integral_options%tol; similarly compute the index only if the integral is significant
                  if (nprocs > 1) then !PARALLEL RUN, RESULTS IN MEMORY OR TO DISK
                     i = p_integral_array%block_offset(quartet_ind)
                     do integral_type=1,number_of_types
                        do k=1,number_of_shell_integrals
                           p_integral_array%a(i+k,integral_type) = shell_integrals(k,integral_type)
                        enddo
                     enddo
                  else !SERIAL RUN, RESULTS IN MEMORY OR TO DISK
                     !This is essentially an inlined version of the function index_1p_continuum:
                     if (is_CCTT) then
                        do integral_type=1,number_of_types
                           do k=1,number_of_shell_integrals
                              !pq = ipair(int_index(1,k)) + int_index(2,k)
                              !ind = ipair(pq) + ipair(int_index(3,k)) + int_index(4,k)
                              pq = this%ordered_pairs(1,ipair(int_index(1,k)) + int_index(2,k))
                              rs = this%ordered_pairs(1,ipair(int_index(3,k)) + int_index(4,k))
                              if (pq > this%last_CT_fn) then !pq is CC pair
                                 ind = this%n_prec_ints + rs + this%n_TT_pairs*(pq-this%last_CT_fn-1)
                              else !rs is CC pair
                                 ind = this%n_prec_ints + pq + this%n_TT_pairs*(rs-this%last_CT_fn-1)
                              endif
                              if (ind > number_of_integrals) stop "indexing error CCTT"
                              p_integral_array%a(ind,integral_type) = shell_integrals(k,integral_type)
                           enddo
                        enddo
                     else !TTTT, CTTT, CTCT group or 2p in the continuum
                        !In this case we simply calculate the global index for each (ab|cd) integral and put it into p_integral_array%a.
                        do integral_type=1,number_of_types
                           do k=1,number_of_shell_integrals
                              !pq = ipair(int_index(1,k)) + int_index(2,k)
                              !ind = ipair(pq) + ipair(int_index(3,k)) + int_index(4,k)
                              pq = this%ordered_pairs(1,ipair(int_index(1,k)) + int_index(2,k))
                              rs = this%ordered_pairs(1,ipair(int_index(3,k)) + int_index(4,k))
                              ind = ipair(max(pq,rs)) + min(pq,rs)
                              if (ind > number_of_integrals) stop "indexing error"
                              p_integral_array%a(ind,integral_type) = shell_integrals(k,integral_type)
                           enddo
                        enddo
                        !Note that using the index function here is considerably slower. That's why the indices are generated directly above.
                        !ind(1:number_of_shell_integrals) = this%integral_index(int_index)
                     endif
                  endif

               enddo !kl
            enddo !ij

            if (integral_class .eq. 1) then
               call GGGG_final !deallocate all intermediate arrays needed for calculation of the 2-electron integrals over GTOs
               !$OMP BARRIER
               !$OMP SINGLE
               write(stdout,'("GGGG interm. arrays have been deallocated")')
               !$OMP END SINGLE
            endif

            !$OMP END PARALLEL

            t2 = omp_get_wtime()
            write(stdout,'("Integral calculation took: ",f8.3," [s]")') t2-t1

            if (integral_class .eq. 4) then
               write(stdout,'("Maximum L needed to converge BGGG class: ",i4)') max_l_BGGG
               write(stdout,'("Maximum L needed to converge BGBG class: ",i4)') max_l_BGBG

               call BG_final
               write(stdout,'(10X,"bto_gto_integrals_mod finalized")')
            endif

         enddo !integral_class

         !If requested print the non-zero integrals
         if (integral_options%print_integrals) then
            call p_integral_array%print(.true.)
         endif

         !Dump all integrals to disk and close the record.
         !Only master writes to the file but every process gets value of: current_pos (integral_options%write) and last_record (p_integral_array%write), i.e. the positions in the file where
         !the master ended up following the calls to these two routines.
         if (integral_storage%on_disk()) then

            write(stdout,'("Saving integrals to disk...")')

            !The first record are the integral options.
            call integral_options%write(lunit,first_record,current_pos)

            !The second record are the ordered integrals: we write them using a method whose choice depends on whether each process calculated all AO integrals or not.
            if (nprocs > 1) then !PARALLEL RUN
               !todo I should gather first the missing parts of the this%a array and write the whole thing including the recalculated offsets to disk.
               call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                            'Writing to disk not implemented for more than 1 MPI task.', 5, 1)
            else !SERIAL RUN
               a = master
               call p_integral_array%write(lunit,current_pos,last_record,a)
            endif

            !Every process closes the record so that they all keep identical header information.
            call integral_storage%integral_file%close_record(header,first_record,last_record)

            err = p_integral_array%final()
            if (err == 0) then
               call xermsg ('atomic_orbital_basis_obj', 'two_electron_integrals', &
                            'Deallocation of the temporary integral array failed.', 5, 1)
            end if

            nullify(p_integral_array)
            !FROM NOW ON p_integral_array is nullified

            write(stdout,'("...done")')

         endif

         write(stdout,'("<---------","atomic_orbital_basis_obj:two_electron_integrals")')

         end_t = omp_get_wtime()
         write(stdout,'("Two_electron_integrals took [s]: ",f25.15)') end_t-start_t

         call mpi_mod_barrier(err)

   end subroutine two_electron_integrals

   subroutine shell_quartet_two_electron_integrals(this,ij,kl,integral_options,int_index,integrals,number_of_shell_integrals)
      use const, only: number_of_types_el_ints
      use cgto_integrals, only: GGGG_shell_integrals
      use bto_gto_integrals_mod, only: BBGG_shell_integrals, BGGG_shell_integrals, &
                                       BGBG_shell_integrals, lebedev_BGGG_shell_integrals
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: ij, kl
      class(integral_options_obj), intent(in) :: integral_options
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)
      integer, intent(out) :: number_of_shell_integrals

      integer :: i, j, k, l, shell_quartet_type, ab, cd, i_abs, j_abs, k_abs, l_abs

        !Relative indices (i.e. within their own shell type) of the shells: i,j,k,l.
        !Absolute indices of the shells: i_abs,j_abs,k_abs,l_abs
        if (this%shell_pair_type(ij) .ge. this%shell_pair_type(kl)) then
           ab = this%shell_pair_type(ij)
           cd = this%shell_pair_type(kl)
           i = this%shell_pair_indices(1,ij)
           j = this%shell_pair_indices(2,ij)
           k = this%shell_pair_indices(1,kl)
           l = this%shell_pair_indices(2,kl)
           i_abs = this%shell_pair_indices(3,ij)
           j_abs = this%shell_pair_indices(4,ij)
           k_abs = this%shell_pair_indices(3,kl)
           l_abs = this%shell_pair_indices(4,kl)
        else
           ab = this%shell_pair_type(kl)
           cd = this%shell_pair_type(ij)
           i = this%shell_pair_indices(1,kl)
           j = this%shell_pair_indices(2,kl)
           k = this%shell_pair_indices(1,ij)
           l = this%shell_pair_indices(2,ij)
           i_abs = this%shell_pair_indices(3,kl)
           j_abs = this%shell_pair_indices(4,kl)
           k_abs = this%shell_pair_indices(3,ij)
           l_abs = this%shell_pair_indices(4,ij)
        endif

        !Calculate the index corresponding to the type of the quartet. shell_pair_type = 1 (GG), 2 (BG), 3 (GB), 4 (BB)
        !The shell pairs have been ordered above so that the pair-type index AB is never smaller than the pair-type index CD.
        shell_quartet_type = ab*(ab-1)/2+cd !e.g. 1 = [GG|GG], 3,5,6 = [BG|BG], etc.

!        write(*,'("t",2i3,i3)') ab,cd,shell_quartet_type
        !write(*,'("pairs",2i)') ij,kl

        !Calculte the number of integrals in this quartet of shells by taking the product of the number of functions in each shell.
        number_of_shell_integrals = this % shell_descriptor(5, this % shell_pair_indices(3, ij)) &
                                  * this % shell_descriptor(5, this % shell_pair_indices(4, ij)) &
                                  * this % shell_descriptor(5, this % shell_pair_indices(3, kl)) &
                                  * this % shell_descriptor(5, this % shell_pair_indices(4, kl))

        select case (shell_quartet_type)
           case (1) ![GG|GG]
              !Note that the value of the R-matrix radius was set using GGGG_initialize
              call GGGG_shell_integrals (this % CGTO_shell_data(i), &
                                         this % CGTO_shell_data(j), &
                                         this % CGTO_shell_data(k), &
                                         this % CGTO_shell_data(l), &
                                         i, j, k, l, &
                                         this % shell_descriptor(4,i_abs), &
                                         this % shell_descriptor(4,j_abs), &
                                         this % shell_descriptor(4,k_abs), &
                                         this % shell_descriptor(4,l_abs), &
                                         integral_options % use_spherical_cgto_alg, &
                                         this % integral_data % two_el_column, int_index, integrals)
           case (2) ![BG|GG]
              !print *,'BGGG 2',ij,kl
              if (integral_options%mixed_ints_method .eq. 2) then
                 call lebedev_BGGG_shell_integrals (this % BTO_shell_data(i),  &
                                                    this % CGTO_shell_data(j), &
                                                    this % CGTO_shell_data(k), &
                                                    this % CGTO_shell_data(l), &
                                                    i, j, k, l, &
                                                    this % shell_descriptor(4,i_abs), &
                                                    this % shell_descriptor(4,j_abs), &
                                                    this % shell_descriptor(4,k_abs), &
                                                    this % shell_descriptor(4,l_abs), &
                                                    this % integral_data % two_el_column, int_index, integrals)
              elseif (integral_options%mixed_ints_method .eq. 1 .or. &
                      integral_options%mixed_ints_method .eq. 3) then
                 call BGGG_shell_integrals (this % BTO_shell_data(i),  &
                                            this % CGTO_shell_data(j), &
                                            this % CGTO_shell_data(k), &
                                            this % CGTO_shell_data(l), &
                                            i, j, k, l, &
                                            this % shell_descriptor(4,i_abs), &
                                            this % shell_descriptor(4,j_abs), &
                                            this % shell_descriptor(4,k_abs), &
                                            this % shell_descriptor(4,l_abs), &
                                            this % integral_data % two_el_column, int_index, integrals)
              endif
              !stop "test"
           case (3) ![BG|BG]
              !print *,'BGBG 3',ij,kl
              call BGBG_shell_integrals (this % BTO_shell_data(i),  &
                                         this % CGTO_shell_data(j), &
                                         this % BTO_shell_data(k),  &
                                         this % CGTO_shell_data(l), &
                                         i, j, k, l, &
                                         this % shell_descriptor(4,i_abs), &
                                         this % shell_descriptor(4,j_abs), &
                                         this % shell_descriptor(4,k_abs), &
                                         this % shell_descriptor(4,l_abs), &
                                         this % integral_data % two_el_column, int_index, integrals)
           case (4) ![GB|GG]
              !print *,'GBGG 4',ij,kl
              if (integral_options%mixed_ints_method .eq. 2) then
                 call lebedev_BGGG_shell_integrals (this % BTO_shell_data(j),  &
                                                    this % CGTO_shell_data(i), &
                                                    this % CGTO_shell_data(k), &
                                                    this % CGTO_shell_data(l), &
                                                    j, i, k, l, &
                                                    this % shell_descriptor(4,j_abs), &
                                                    this % shell_descriptor(4,i_abs), &
                                                    this % shell_descriptor(4,k_abs), &
                                                    this % shell_descriptor(4,l_abs), &
                                                    this % integral_data % two_el_column, int_index, integrals)
              elseif (integral_options%mixed_ints_method .eq. 1 .or. &
                      integral_options%mixed_ints_method .eq. 3) then
                 call BGGG_shell_integrals (this % BTO_shell_data(j),  &
                                            this % CGTO_shell_data(i), &
                                            this % CGTO_shell_data(k), &
                                            this % CGTO_shell_data(l), &
                                            j, i, k, l, &
                                            this % shell_descriptor(4,j_abs), &
                                            this % shell_descriptor(4,i_abs), &
                                            this % shell_descriptor(4,k_abs), &
                                            this % shell_descriptor(4,l_abs), &
                                            this % integral_data % two_el_column, int_index, integrals)
              endif
              !stop "test"
           case (5) ![GB|BG]
              !print *,'GBBG 5',ij,kl
              call BGBG_shell_integrals (this % BTO_shell_data(j),  &
                                         this % CGTO_shell_data(i), &
                                         this % BTO_shell_data(k),  &
                                         this % CGTO_shell_data(l), &
                                         j, i, k, l, &
                                         this % shell_descriptor(4,j_abs), &
                                         this % shell_descriptor(4,i_abs), &
                                         this % shell_descriptor(4,k_abs), &
                                         this % shell_descriptor(4,l_abs), &
                                         this % integral_data % two_el_column, int_index, integrals)
           case (6) ![GB|GB]
              !print *,'GBGB 6',ij,kl
              call BGBG_shell_integrals (this % BTO_shell_data(j),  &
                                         this % CGTO_shell_data(i), &
                                         this % BTO_shell_data(l),  &
                                         this % CGTO_shell_data(k), &
                                         j, i, l, k, &
                                         this % shell_descriptor(4,j_abs), &
                                         this % shell_descriptor(4,i_abs), &
                                         this % shell_descriptor(4,l_abs), &
                                         this % shell_descriptor(4,k_abs), &
                                         this % integral_data % two_el_column, int_index, integrals)
           case (7) ![BB|GG]
              !print *,'BBGG 7',ij,kl
              call BBGG_shell_integrals (this % BTO_shell_data(i),  &
                                         this % BTO_shell_data(j),  &
                                         this % CGTO_shell_data(k), &
                                         this % CGTO_shell_data(l), &
                                         i, j, k, l, &
                                         this % shell_descriptor(4,i_abs), &
                                         this % shell_descriptor(4,j_abs), &
                                         this % shell_descriptor(4,k_abs), &
                                         this % shell_descriptor(4,l_abs), &
                                         this % integral_data % two_el_column, int_index, integrals)
           case (8) ![BB|BG]
              call xermsg ('atomic_orbital_basis_obj', 'shell_quartet_two_electron_integrals', &
                           '[BB|BG] type not implemented: use two_p_continuum = .false.', 1, 1)
           case (9) ![BB|GB]
              call xermsg ('atomic_orbital_basis_obj', 'shell_quartet_two_electron_integrals', &
                           '[BB|GB] type not implemented: use two_p_continuum = .false.', 2, 1)
           case(10) ![BB|BB]
              call xermsg ('atomic_orbital_basis_obj', 'shell_quartet_two_electron_integrals', &
                           '[BB|BB] type not implemented: use two_p_continuum = .false.', 3, 1)
           case default
              call xermsg ('atomic_orbital_basis_obj', 'shell_quartet_two_electron_integrals', &
                           'Error in shell_quartet_type: the allowed values are 1-10.', 4, 1)
        end select

   end subroutine shell_quartet_two_electron_integrals

   subroutine normalize_continuum(this,a)
      use gto_routines, only: cms_gto_norm
      implicit none
      class(atomic_orbital_basis_obj) :: this
      real(kind=cfp), intent(in) :: a

      integer :: i
      real(kind=cfp) :: norm, r1, r2

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj','normalize_continuum', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (a > 0.0_cfp) then
            !Calculate norms of all CGTO continuum shells to the R-matrix sphere (if there are any). This is used when calculating the tails.
            if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
               do i=1,size(this%CGTO_shell_data)
                  !Normalize all continuum CGTO shells to the R-matrix sphere:
                  if (this%CGTO_shell_data(i)%is_continuum()) then
                     norm = cms_gto_norm(a, &
                                         this % CGTO_shell_data(i) % l, &
                                         this % CGTO_shell_data(i) % number_of_primitives, &
                                         this % CGTO_shell_data(i) % exponents, &
                                         this % CGTO_shell_data(i) % contractions, &
                                         this % CGTO_shell_data(i) % norm, &
                                         this % CGTO_shell_data(i) % norms)
                     this%CGTO_shell_data(i)%norms(:) = norm*this%CGTO_shell_data(i)%norm*this%CGTO_shell_data(i)%norms(:)
                     this%CGTO_shell_data(i)%norm = 1.0_cfp
                  endif
               enddo !i
            endif
            !Check that the value of the R-matrix radius is compatible with the radial extent of the BTOs (if there are any).
            if (allocated(this % BTO_shell_data) .and. size(this % BTO_shell_data) > 0) then
               do i=1,size(this%BTO_shell_data)
                  call this%BTO_shell_data(i)%bspline_grid%bspline_range(this%BTO_shell_data(i)%bspline_index,r1,r2)
                  if (r2 > a) then
                     call xermsg ('atomic_orbital_basis_obj', 'normalize_continuum', &
                                  'The BTOs included in the basis must not leak outside of the R-matrix sphere.', 2, 1)
                  end if
               enddo !i
            endif
         endif

   end subroutine normalize_continuum

   function integral_index(this,integral_type,bf_indices)
      use const
      use special_functions, only: ipair
      implicit none
      class(atomic_orbital_basis_obj) :: this
      character(len=*), intent(in) :: integral_type
      integer, intent(in) :: bf_indices(:,:)
      integer :: integral_index(size(bf_indices,2))

      integer :: ind,i,j,k, iAB, iCD

         !todo chcek for initialization

         if (size(bf_indices,1) .eq. 2) then !1-electron integral index

            !todo this can be improved by swapping the do-loop with select case for
            !each integral type. Best to pack the indexing computation into an
            !elementary routine.
            do k=1,size(bf_indices,2)
               i = maxval(bf_indices(1:2,k))
               j = minval(bf_indices(1:2,k))

               ind = ipair(i)+j

               !At the moment all 1-electron atomic integrals are index in the same way.
               select case (integral_type)

                  case (overlap_ints)
                     integral_index(k) = ind
                  case (kinetic_ints)
                     integral_index(k) = ind
                  case (property_ints)
                     integral_index(k) = ind
                  case (nuc_rep_att_ints)
                     integral_index(k) = ind
                  case (one_elham)
                     integral_index(k) = ind
                  case default
                  call xermsg ('atomic_orbital_basis_obj', 'integral_index', &
                               'Unrecognized one electron atomic integral type on input.', 1, 1)

               end select
            enddo !k

         elseif (size(bf_indices,1) .eq. 4) then !2-electron integral index
            do i=1,size(bf_indices,2)

               !todo the basis_set_data (i.e. cgto_basis_data_obj) should contain information on whether the index is for 2p in the cont. or 1p in the continuum

               !If needed, permute the indices into the standard order and calculate indices of the AB, CD pairs:
               if (bf_indices(1,i) .ge. bf_indices(2,i)) then !a .ge. b
                  iAB = ipair(bf_indices(1,i))+bf_indices(2,i) !a*(a-1)/2+b
               else
                  iAB = ipair(bf_indices(2,i))+bf_indices(1,i) !b*(b-1)/2+a
               endif
               if (bf_indices(3,i) .ge. bf_indices(4,i)) then !c .ge. d
                  iCD = ipair(bf_indices(3,i))+bf_indices(4,i) !c*(c-1)/2+d
               else
                  iCD = ipair(bf_indices(4,i))+bf_indices(3,i) !d*(d-1)/2+c
               endif

               iAB = this%ordered_pairs(1,iAB)
               iCD = this%ordered_pairs(1,iCD)

               !Standard indexing for TTTT CTCT CTTT classes and for 2p in the continuum.
               !Swap the (ab) and (cd) pairs of functions if neccessary.
               if (iAB < iCD) then
                  integral_index(i) = ipair(iCD)+iAB
               else
                  integral_index(i) = ipair(iAB)+iCD
               endif

               !todo check for 2p in the continuum or not!
               !Special indexing for CCTT class in case 1p in the continuum
               if (iAB .le. this%n_TT_pairs .or. iCD .le. this%n_TT_pairs) then !AB and/or CD are TT
                  if (iAB > this%last_CT_fn) then !iAB is CC pair
                     integral_index(i) = this%n_prec_ints + iCD + this%n_TT_pairs*(iAB-this%last_CT_fn-1)
                  elseif(iCD > this%last_CT_fn) then !iCD is CC pair
                     integral_index(i) = this%n_prec_ints + iAB + this%n_TT_pairs*(iCD-this%last_CT_fn-1)
                  endif
               endif

            enddo
         else
            call xermsg ('atomic_orbital_basis_obj', 'integral_index', &
                         'On input the number of basis function indices per integral must be either 2 or 4.', 3, 1)
         endif

   end function integral_index

   subroutine assemble_overlap_matrix(this,ao_integrals,overlap_matrix)
      use const, only: overlap_ints
      implicit none
      class(atomic_orbital_basis_obj) :: this
      type(p2d_array_obj) :: ao_integrals
      real(kind=cfp), allocatable :: overlap_matrix(:,:)

      integer :: ao_indices(2,1), ind(1), err, i, j, overlaps_column

         write(stdout,'("--------->","atomic_orbital_basis_obj:assemble_overlap_matrix")')

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'assemble_overlap_matrix', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (.not.(associated(ao_integrals%a))) then
            call xermsg ('atomic_orbital_basis_objj', 'assemble_overlap_matrix', &
                         'The input array of overlap integrals has not been allocated.', 11, 1)
         else

            !Transfer the overlap_integrals into overlap_matrix
            overlaps_column = ao_integrals%find_column_matching_name(overlap_ints)
            write(stdout,'("Atomic overlap integrals found in column number: ",i4)') overlaps_column

            if (allocated(overlap_matrix)) deallocate(overlap_matrix)
            allocate(overlap_matrix(this%number_of_functions,this%number_of_functions),stat=err)
            if (err /= 0) then
               call xermsg ('atomic_orbital_basis_obj', 'assemble_overlap_matrix', &
                            'Memory allocation 0 failed.', err, 1)
            end if

            write(stdout,'("Transforming the AO integrals into the atomic overlap matrix...")')

            do i=1,this%number_of_functions
               do j=1,i
                  ao_indices(1:2,1) = (/i,j/)
                  ind(1:1) = this%integral_index(overlap_ints,ao_indices)
                  overlap_matrix(j,i) = ao_integrals%a(ind(1),overlaps_column)
                  overlap_matrix(i,j) = overlap_matrix(j,i)
               enddo !j
            enddo !i

             write(stdout,'("...done")')

         endif

         write(stdout,'("<---------","atomic_orbital_basis_obj:assemble_overlap_matrix")')

   end subroutine assemble_overlap_matrix

   function get_basis_name(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      character(len=line_len) :: get_basis_name

         get_basis_name = "atomic_orbital_basis_obj"

   end function get_basis_name

   function get_shell_name(this,i)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      character(len=line_len) :: get_shell_name
      integer, intent(in) :: i

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_shell_name', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this % number_of_shells) then
            call xermsg ('atomic_orbital_basis_obj', 'get_shell_name', &
                         'On input the value of i was out of range.', 2, 1)
         end if

         if (this%shell_descriptor(1,i) .eq. 1) then
            get_shell_name = this%CGTO_shell_data(this%shell_descriptor(2,i))%name()
         elseif (this%shell_descriptor(1,i) .eq. 2) then
            get_shell_name = this%BTO_shell_data(this%shell_descriptor(2,i))%name()
         else
            call xermsg ('atomic_orbital_basis_obj', 'get_shell_name', 'Error in this%shell_descriptor.',3,1)
         endif

   end function get_shell_name

   subroutine get_shell_data(this,i,shell_data)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: i
      class(shell_data_obj), intent(out) :: shell_data

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_shell_data', 'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this%number_of_shells) then
            call xermsg ('atomic_orbital_basis_obj', 'get_shell_data', 'On input the value of i was out of range.', 2, 1)
         end if

         select type (shell => shell_data)
            type is (CGTO_shell_data_obj)
               if (this % shell_descriptor(1,i) == 2) then
                  call xermsg ('atomic_orbital_basis_obj', 'get_shell_data', &
                               'Requested a BTO shell but on input shell_data was CGTO_shell_data_obj.', 3, 1)
               end if
               shell = this%CGTO_shell_data(this%shell_descriptor(2,i))
            type is (BTO_shell_data_obj)
               if (this%shell_descriptor(1,i) == 1) then
                  call xermsg ('atomic_orbital_basis_obj', 'get_shell_data', &
                               'Requested a CGTO shell but on input shell_data was BTO_shell_data_obj.', 4, 1)
               end if
               shell = this%BTO_shell_data(this%shell_descriptor(2,i))
            class default
               call xermsg ('atomic_orbital_basis_obj', 'get_shell_data', 'Not implemented for this shell type.',5,1)
         end select

   end subroutine get_shell_data

   subroutine get_all_CGTO_shells(this,CGTO_shells,number_of_cgto_shells)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(out) :: number_of_cgto_shells
      type(CGTO_shell_data_obj), allocatable :: CGTO_shells(:)

      integer :: err

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_all_CGTO_shells', 'The basis set has not been initialized.', 1, 1)
         end if

         if (this % n(1) == 0) then
            call xermsg ('atomic_orbital_basis_obj', 'get_all_CGTO_shells', &
                         'This basis set does not contain any CGTO shells.', 2, 1)
         end if

         if (allocated(CGTO_shells)) deallocate(CGTO_shells)
         allocate(CGTO_shells,source=this%CGTO_shell_data,stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'get_all_CGTO_shells', 'Memory allocation failed.',err,1)

         number_of_cgto_shells = size(this%CGTO_shell_data)

   end subroutine get_all_CGTO_shells

   subroutine get_all_BTO_shells(this,BTO_shells,number_of_bto_shells)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(out) :: number_of_bto_shells
      type(BTO_shell_data_obj), allocatable :: BTO_shells(:)

      integer :: err

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_all_BTO_shells', 'The basis set has not been initialized.', 1, 1)
         end if

         if (this % n(2) == 0) then
            call xermsg ('atomic_orbital_basis_obj', 'get_all_BTO_shells', 'This basis set does not contain any BTO shells.', 2, 1)
         end if

         if (allocated(BTO_shells)) deallocate(BTO_shells)
         allocate(BTO_shells,source=this%BTO_shell_data,stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'get_all_BTO_shells', 'Memory allocation failed.',err,1)

         number_of_bto_shells = size(this%BTO_shell_data)

   end subroutine get_all_BTO_shells

   subroutine get_bspline_grid(this,bspline_grid)
      use bspline_grid_mod
      implicit none
      class(atomic_orbital_basis_obj) :: this
      type(bspline_grid_obj) :: bspline_grid

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_bspline_grid', 'The basis set has not been initialized.', 1, 1)
         end if

         if (this % n(2) == 0) then
            call xermsg ('atomic_orbital_basis_obj', 'get_bspline_grid', 'This basis set does not contain any BTO shells.', 2, 1)
         end if

         bspline_grid = this%BTO_shell_data(1)%bspline_grid

   end subroutine get_bspline_grid

   function is_initialized(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      logical :: is_initialized

         is_initialized = this%initialized

   end function is_initialized

   subroutine get_symmetry_flags(this,irr,list,from,to)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: irr, from, to ! indexes from and to refer to the range of atomic shells between which to check for symmetry conformity
      logical, allocatable :: list(:)

      integer :: err, i, ind, n, l, m, max_l
      integer, allocatable :: lm_irr(:)
      character(len=sym_op_nam_len) :: nam

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_symmetry_flags', &
                         'The object has not been initialized or not all shells have been read-in.', 1, 1)
         end if

         if (irr .le. 0 .or. irr > this%symmetry_data%get_no_irrep(this%symmetry_data%get_pg())) &
            call xermsg ('atomic_orbital_basis_obj', 'get_symmetry_flags', 'On input the value of irr was out of range.', 2, 1)

         max_l = 0
         if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
            max_l = max(max_l,maxval(this % CGTO_shell_data(:) % l))
         end if
         if (allocated(this % BTO_shell_data)  .and. size(this % BTO_shell_data)  > 0) then
            max_l = max(max_l,maxval(this % BTO_shell_data(:)  % l))
         end if

         !remember that we're ignoring the irr parameter since we don't implement AO symmetry yet!
         if (allocated(list)) deallocate(list)
         allocate(list(this%number_of_functions),lm_irr((max_l+1)**2),stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'get_symmetry_flags', 'Memory allocation error.',err, 1)
         list = .false.

         if ( from.gt.to ) return

         do l=0,max_l
            do m=-l,l
               i = l*l+l+m+1
               lm_irr(i) = this%symmetry_data%sph_harm_pg_sym(l,m,nam)
            enddo
         enddo !l

         ind = 0
         do i=1,this%number_of_shells,1

            if (this%shell_descriptor(1,i) .eq. 1) then    !the i-shell is a CGTO shell
               n = this%CGTO_shell_data(this%shell_descriptor(2,i))%number_of_functions
               l = this%CGTO_shell_data(this%shell_descriptor(2,i))%l
            elseif (this%shell_descriptor(1,i) .eq. 2) then !the i-shell is a BTO shell
               n = this%BTO_shell_data(this%shell_descriptor(2,i))%number_of_functions
               l = this%BTO_shell_data(this%shell_descriptor(2,i))%l
            else !error
               call xermsg ('atomic_orbital_basis_obj', 'get_symmetry_flags', &
                            'The shell type A must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 3, 1)
            endif

            if ( i < from .or. i > to ) then
               ind = ind + n
            else
               do m=-l,l
                  ind = ind + 1
                  if (lm_irr(l*l+l+m+1) .eq. irr) list(ind) = .true.
               enddo !m
            endif

         enddo !i

         if (ind /= this % number_of_functions) then
            call xermsg ('atomic_orbital_basis_obj', 'get_symmetry_flags', &
                         'Inconsistency in internal data: programming error or data corruption.', 4, 1)
         end if

   end subroutine get_symmetry_flags

   subroutine get_continuum_flags(this,irr,list)
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: irr
      logical, allocatable :: list(:)

      integer :: err, i, ind, n, l, m, max_l
      integer, allocatable :: lm_irr(:)
      character(len=sym_op_nam_len) :: nam

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_continuum_flags', &
                         'The object has not been initialized or not all shells have been read-in.', 1, 1)
         end if

         if (irr .le. 0 .or. irr > this%symmetry_data%get_no_irrep(this%symmetry_data%get_pg())) &
         &call xermsg ('atomic_orbital_basis_obj', 'get_continuum_flags', 'On input the value of irr was out of range.', 2, 1)

         max_l = 0
         if (allocated(this % CGTO_shell_data) .and. size(this % CGTO_shell_data) > 0) then
            max_l = max(max_l,maxval(this % CGTO_shell_data(:) % l))
         end if
         if (allocated(this % BTO_shell_data)  .and. size(this % BTO_shell_data)  > 0) then
            max_l = max(max_l,maxval(this % BTO_shell_data(:)  % l))
         end if

         !remember that we're ignoring the irr parameter since we don't implement AO symmetry yet!
         if (allocated(list)) deallocate(list)
         allocate(list(this%number_of_functions),lm_irr((max_l+1)**2),stat=err)
         if (err .ne. 0) call xermsg ('atomic_orbital_basis_obj', 'get_continuum_flags', 'Memory allocation error.',err, 1)
         list = .false.

         do l=0,max_l
            do m=-l,l
               i = l*l+l+m+1
               lm_irr(i) = this%symmetry_data%sph_harm_pg_sym(l,m,nam)
            enddo
         enddo !l

         ind = 0
         do i=1,this%number_of_shells

            if (this%shell_descriptor(1,i) .eq. 1) then    !the i-shell is a CGTO shell
               n = this%CGTO_shell_data(this%shell_descriptor(2,i))%number_of_functions
               if (this%shell_descriptor(3,i) .eq. 1) then
                  l = this%CGTO_shell_data(this%shell_descriptor(2,i))%l
                  do m=-l,l
                     ind = ind + 1
                     if (lm_irr(l*l+l+m+1) .eq. irr) list(ind) = .true.
                  enddo !m
               else
                  ind = ind + n
               endif
            elseif (this%shell_descriptor(1,i) .eq. 2) then !the i-shell is a BTO shell
               n = this%BTO_shell_data(this%shell_descriptor(2,i))%number_of_functions
               if (this%BTO_shell_data(this%shell_descriptor(2,i))%non_zero_at_boundary) then
                  l = this%BTO_shell_data(this%shell_descriptor(2,i))%l
                  do m=-l,l
                     ind = ind + 1
                     if (lm_irr(l*l+l+m+1) .eq. irr) list(ind) = .true.
                  enddo !m
               else
                  ind = ind + n
               endif
            else !error
               call xermsg ('atomic_orbital_basis_obj', 'get_continuum_flags', &
                            'The shell type A must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 3, 1)
            endif

         enddo !i

         if (ind /= this % number_of_functions) then
            call xermsg ('atomic_orbital_basis_obj', 'get_continuum_flags', &
                         'Inconsistency in internal data: programming error or data corruption.', 4, 1)
         end if

   end subroutine get_continuum_flags

   subroutine get_continuum_l_range(this,min_l,max_l)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(out) :: min_l, max_l

      integer :: ind, i

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'get_continuum_l_range', &
                         'The object has not been initialized or not all shells have been read-in.', 1, 1)
         end if

         max_l = -1
         min_l = huge(min_l)

         do i=1,this%number_of_shells

            if (this%shell_descriptor(1,i) .eq. 1) then    !the i-shell is a CGTO shell
               ind = this%shell_descriptor(2,i)
               if (this%CGTO_shell_data(ind)%is_continuum()) then
                  max_l = max(max_l,this%CGTO_shell_data(ind)%l)
                  min_l = min(min_l,this%CGTO_shell_data(ind)%l)
               endif
            elseif (this%shell_descriptor(1,i) .eq. 2) then !the i-shell is a BTO shell
               ind = this%shell_descriptor(2,i)
               if (this%BTO_shell_data(ind)%non_zero_at_boundary) then
                  max_l = max(max_l,this%BTO_shell_data(ind)%l)
                  min_l = min(min_l,this%BTO_shell_data(ind)%l)
               endif
            else !error
               call xermsg ('atomic_orbital_basis_obj', 'get_continuum_l_range', &
                            'The shell type A must be one of: CGTO_shell_data_obj, BTO_shell_data_obj.', 2, 1)
            endif

         enddo !i

         if (min_l > max_l) then
            call xermsg ('atomic_orbital_basis_obj', 'get_continuum_l_range', &
                         'The atomic basis does not include any continuum functions.', 3, 1)
         end if

   end subroutine get_continuum_l_range

   subroutine calculate_amplitudes(this,a,normalize_to_a,amplitudes,continuum_channels)
      use gto_routines, only: CGTO_amplitude
      implicit none
      class(atomic_orbital_basis_obj) :: this
      real(kind=cfp), intent(in) :: a
      logical, intent(in) :: normalize_to_a
      integer, allocatable :: continuum_channels(:,:)
      real(kind=cfp), allocatable :: amplitudes(:,:)

      integer :: i, ind, cnt, l, m, m_cgto, m_bto, err, l_max, l_min, n_channels, channel_index
      real(kind=cfp) :: cgto_red_b_amp, bto_red_b_amp
      character(len=sym_op_nam_len) :: nam
      integer, parameter :: der = 0, der1 = 1

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'calculate_amplitudes', &
                         'The object has not been initialized or not all shells have been read-in.', 1, 1)
         end if

         write(stdout,'("--------->","atomic_orbital_basis_obj:calculate_amplitudes")')

         if (normalize_to_a) then
            write(stdout,'("Continuum functions will be normalized to the requested radius: ",e25.15)') a
            call this%normalize_continuum(a)
         else
            write(stdout,'("Continuum functions will NOT be normalized to the requested radius: ",e25.15)') a
         endif

         call this%get_continuum_l_range(l_min,l_max)

         n_channels = (l_max+1)**2-(l_min)**2
         if (allocated(continuum_channels)) deallocate(continuum_channels)
         allocate(continuum_channels(3,n_channels),stat=err)
         if (err .ne. 0) call xermsg('atomic_orbital_basis_obj', 'continuum_channels', 'Memory allocation failed.',err, 1)

         i = 0
         do l=l_min,l_max
            do m=-l,l
               i = i + 1
               continuum_channels(1:2,i) = (/m,l/)
               continuum_channels(3,i) = this%symmetry_data%sph_harm_pg_sym(l,m,nam)
            enddo !m
         enddo !l

         if (allocated(amplitudes)) deallocate(amplitudes)
         allocate(amplitudes(n_channels,this%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('atomic_orbital_basis_obj', 'calculate_amplitudes', 'Memory allocation failed.',err, 1)

         amplitudes = 0.0_cfp

         cnt = 0
         do i=1,this%number_of_shells
            if (this%shell_descriptor(1,i) .eq. 1) then !CGTO shell
               ind = this%shell_descriptor(2,i)

               !Skip all target shells
               if (this%CGTO_shell_data(ind)%is_continuum()) then
                  l = this%CGTO_shell_data(ind)%l
                  channel_index = l**2 - l_min**2
                  cgto_red_b_amp = CGTO_amplitude(a, l, this % CGTO_shell_data(ind) % number_of_primitives, &
                                                        this % CGTO_shell_data(ind) % norm, &
                                                        this % CGTO_shell_data(ind) % norms, &
                                                        this % CGTO_shell_data(ind) % contractions, &
                                                        this % CGTO_shell_data(ind) % exponents)
                  do m_cgto=-l,l
                     m = m_cgto ! = delta_{m,mp}
                     amplitudes(channel_index+l+m+1,cnt+m_cgto+l+1) = cgto_red_b_amp
                     if (cgto_red_b_amp /= 0.0_cfp) then
                        write(stdout,'(i0,1x,i0,e25.15)') &
                            channel_index + l + m + 1, cnt + m_cgto + l + 1, amplitudes(channel_index+l+m+1,cnt+m_cgto+l+1)
                     end if
                  enddo !m_cgto
               endif

               cnt = cnt+this%CGTO_shell_data(ind)%number_of_functions

            elseif (this%shell_descriptor(1,i) .eq. 2) then !BTO shell
               ind = this%shell_descriptor(2,i)

               if (this%BTO_shell_data(ind)%non_zero_at_boundary) then
                  l = this%BTO_shell_data(ind)%l
                  channel_index = l**2 - l_min**2
                  bto_red_b_amp = this % BTO_shell_data(ind) % bspline_grid % bspline_amplitude(&
                        a, &
                        this % BTO_shell_data(ind) % norm, &
                        this % BTO_shell_data(ind) % bspline_index, &
                        der)
                  do m_bto=-l,l
                     m = m_bto ! = delta_{m,mp}
                     amplitudes(channel_index+l+m+1,cnt+m_bto+l+1) = bto_red_b_amp
                     if (bto_red_b_amp .ne. 0.0_cfp) then
                        write(stdout,'(i0,1x,i0,e25.15)') &
                            channel_index + l + m + 1, cnt + m_bto + l + 1, amplitudes(channel_index+l+m+1,cnt+m_bto+l+1)
                     end if
                  enddo !m_bto
               endif

               cnt = cnt+this%BTO_shell_data(ind)%number_of_functions

            else
               call xermsg ('atomic_orbital_basis_obj', 'calculate_amplitudes', 'Error in this%shell_descriptor.',2,1)
            endif
         enddo !i

         write(stdout,'("<---------","atomic_orbital_basis_obj:calculate_amplitudes")')

   end subroutine calculate_amplitudes

   function contains_btos(this)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      logical :: contains_btos

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'contains_btos', &
                         'The object has not been initialized or not all shells have been read-in.', 1, 1)
         end if

         contains_btos = .false.
         if (allocated(this % BTO_shell_data) .and. size(this % BTO_shell_data) > 0) contains_btos = .true.

   end function contains_btos

   subroutine eval_basis(this,r,n_points,ao_basis_at_r)
      implicit none
      class(atomic_orbital_basis_obj) :: this
      integer, intent(in) :: n_points
      real(kind=cfp), intent(in) :: r(3,n_points)
      real(kind=cfp), allocatable :: ao_basis_at_r(:,:)

      real(kind=cfp), allocatable :: shell_val(:,:)
      integer :: i, j, n_max, err, n, starting_index
      logical :: different_r

         if (.not. this % initialized) then
            call xermsg ('atomic_orbital_basis_obj', 'eval_basis', &
                         'The object has not been initialized or not all shells have been read-in.', 1, 1)
         end if

         write(stdout,'("--------->","atomic_orbital_basis_obj:eval_basis")')

         !First check if the points on input are different to the ones used in
         !the previous call (if any).
         different_r = .false.
         if (size(r,2) .ne. size(this%r,2)) different_r = .true.

         if (.not.different_r) then
            do i=1,size(r,2)
               do j=1,3
                  if (r(j,i) .ne. this%r(j,i)) then
                     different_r = .true.
                     exit
                  endif
               enddo
               if (different_r) exit
            enddo !i
         endif

         if (different_r) then
            if (allocated(this%r)) deallocate(this%r)
            if (allocated(ao_basis_at_r)) deallocate(ao_basis_at_r)
            if (allocated(this%ao_basis_at_r)) deallocate(this%ao_basis_at_r)

            allocate(this%r,source=r,stat=err)
            if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj','eval_orbital', 'Memory allocation 1 failed.', err, 1)
            allocate(ao_basis_at_r(this % number_of_functions, n_points), &
                     this % ao_basis_at_r(this % number_of_functions, n_points), stat = err)
            if (err .ne. 0) call xermsg('atomic_orbital_basis_obj','eval_basis','Memory allocation 2 error.',err,1)

            n_max = maxval(this%shell_descriptor(5,:)) !maximum number of functions in a single shell
            this%ao_basis_at_r = 0.0_cfp

            ! Note: This parallel section needs to use DEFAULT(SHARED) to allow work with polymorphic objects with gfortran.
            !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,n,j,starting_index,shell_val,err) SHARED(this,n_max,n_points,r)
            allocate(shell_val(n_max,n_points),stat=err)
            if (err .ne. 0) call xermsg('atomic_orbital_basis_obj','eval_basis','Memory allocation 3 error.',err,1)
            shell_val = 0.0_cfp
            !$OMP DO
            do i=1,this%number_of_shells
               n = this%shell_descriptor(5,i) !number of functions in this shell
               starting_index = this%shell_descriptor(4,i) !starting index for the functions in this shell

               j = this%shell_descriptor(2,i) !index of the shell within its own type
               if (this%shell_descriptor(1,i) .eq. 1) then !CGTO shell
                  shell_val(1:n,1:n_points) = this%CGTO_shell_data(j)%eval(r,n_points)
               elseif (this%shell_descriptor(1,i) .eq. 2) then !BTO shell
                  shell_val(1:n,1:n_points) = this%BTO_shell_data(j)%eval(r,n_points)
               else
                  call xermsg('atomic_orbital_basis_obj','eval_basis','Error or unimplemented for this shell type.',2,1)
               endif

               this%ao_basis_at_r(starting_index:starting_index+n-1,1:n_points) = shell_val(1:n,1:n_points)
            enddo !i
            !$OMP END DO
            !$OMP END PARALLEL

         elseif (.not.(allocated(ao_basis_at_r))) then
            allocate(ao_basis_at_r(this%number_of_functions,n_points),stat=err)
            if (err .ne. 0) call xermsg('atomic_orbital_basis_obj','eval_basis','Memory allocation 4 error.',err,1)
         endif

         ao_basis_at_r = this%ao_basis_at_r

         write(stdout,'("<---------","atomic_orbital_basis_obj:eval_basis")')

   end subroutine eval_basis

end module atomic_basis_mod
