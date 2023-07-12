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
module molecular_basis_mod
   use precisn
   use const, only: line_len, no_header, stdout
   use symmetry
   use integral_storage_mod
   use mpi_mod
   use atomic_basis_mod
   use basis_data_generic_mod
   use mpi_memory_mod
 
   implicit none
 
   private
 
   !Only the objects themselves are visible from outside of this module.
   public molecular_orbital_basis_obj, orbital_data_obj, sym_ortho_io
 
   type, extends(basis_data_generic_obj) :: molecular_orbital_basis_obj
       !> Pointer to the AO basis object in terms of which this MO set is specified. This pointer must be associated before calling init.
       class(atomic_orbital_basis_obj), pointer :: ao_basis => null()
       !> This pointer must be associated by the user before calling one_electron_integrals or two_electron_integrals to the set of atomic integrals to be transformed.
       class(integral_storage_obj), pointer :: ao_integral_storage => null()
       !> Orbital data sets (typically one for each IRR). This array will have size number_of_shells but only during initialization, i.e. while the orbital sets are being added.
       !> This array is discarded following initialization.
       type(orbital_data_obj), allocatable, private :: orbital_data(:)
       !> Array which maps the absolute index of the orbital in the basis to the index of the orbital within the orbital_data array (absolute_to_relative).
       !> relative_to_absolute is the inverse of absolute_to_relative.
       !> absolute_to_relative(2,absolute_index) = IRR of the orbital with absolute index absolute_index.
       !> absolute_to_relative(1,absolute_index) = index of the orbital within the set of orbitals with the same IRR.
       integer, allocatable, private :: absolute_to_relative(:,:), relative_to_absolute(:,:)
       !> Point group symmetry of the orbitals in this set.
       integer :: pg = -1
       !> Number of irreducible representations.
       integer :: no_irr = 0
       !> so2mo_range: indices of the first and the last MO to which a given AO contributes. mo2so_range is the opposite for the MOs.
       !> These arrays have a meaning identical to the one explained in:
       !> S.Yamamoto, U. Nagashima, CPC 166 (2005) 58-65.
       integer, allocatable :: so2mo_range(:,:), mo2so_range(:,:)
       !> Arrays used to compute the index for the 2-particle symmetric (AB|CD) integrals. For details see add_function.
       integer, allocatable :: block_offset(:), sym_offset(:)
       !> For each orbital in the basis this list contains the value .true./.false. depending on whether the orbital is continuum or not.
       logical, private, allocatable :: is_continuum(:)
       !> Indices of the 2-electron integrals for each type of integral. Determined by two_electron_integrals.
       !> The pointer attribute is necessary for shared memory usage.
       integer, pointer :: ijkl_indices(:,:) => null()
       !> Used in case of shared memory allocation of ijkl_indices. It stores the MPI window that ijkl_indices is stored in.
       !> If it is -1 then we are in the non-shared mode for ijkl_indices.
       integer, private :: shared_window_ijkl_indices  = -1
       !> Index of the last integral stored in ijkl_indices.
       integer :: ind_ijkl_integral = 0
       !> Used only during initialization while the orbital sets are being added.
       logical, private, allocatable :: sets_added(:)
       !> Set to .true. following input of all basis functions for which space has been allocated.
       logical, private :: initialized = .false.
       !> Set to .true. following a call to init.
       logical, private :: init_called = .false.
   contains
       !> Allocates space for a given number and type of shells of basis functions.
       procedure :: init
       !> Finalizes the basis set.
       procedure :: final
       !> Adds data for one shell into the basis set.
       procedure :: add_shell
       !> Prints the orbital basis set data to the stdout unit.
       procedure :: print => print_molecular_orbital_basis_obj
       !> Prints the orbital coefficients to the stdout unit.
       procedure :: print_orbitals
       !> Prints an orbital table showing the orbitals sorted in energy.
       procedure :: print_energy_sorted_orbital_table
       !> Calculates the values in arrays so2mo, mo2so, absolute_to_relative, relative_to_absolute, block_offset, sym_offset, is_continuum.
       procedure, private :: determine_auxiliary_indices
       !> Performs Gramm-Schmidt or Symmetric orthogonalization. The type of orthogonalization is selected using the logical
       !> parameters 'gramm_schmidt' or 'symmetric'. By default for G-S orthogonalization all orbitals are orthogonalized
       !> starting from the orbital with index 1 in the basis set. By default for symmetric orthogonalization all orbitals
       !> are orthogonalized. The range of 'active' and 'passive' orbitals can be selected by specifying the optional integers
       !> active_start, active_end and passive_start, passive_end which specify the range of indices for the orbitals to orthogonalize
       !> and not to orthogonalize. The format for the indices is (/num,sym/) where num is the (external) number of the orbital
       !> and sym is its symmetry. For symmetric orthogonalization the set of 'passive' orbitals is used only to check at the end
       !> that all 'active' orbitals are orthogonal to the 'passive' orbitals. The orthogonalization requires the AO overlap integrals array. 
       !> If symmetric orthogonalization is required then the data structure of type sym_ortho_io must be also present on input.
       !> It specifies the deletion thresholds for each symmetry. On output it contains the list of orbitals to delete
       !> (i.e. those that didn't pass the deletion threshold criterion) in the form of a logical array which marks the orbitals
       !> for deletion.
       procedure :: orthogonalize
       !> Deletes specified orbitals of a given symmetry. The orbitals to delete are marked in an input logical array which must
       !> have size equal to the number of orbitals in the given symmetry.
       procedure :: delete_orbitals
       !> Calculates and stores 1-electron integrals for all pairs of shells in the basis. The atomic integrals to be transformed
       !> are input via the type-bound pointer ao_integral_storage.
       procedure :: one_electron_integrals
       !> Calculates and stores 2-electron integrals for all pairs of shells in the basis. The atomic integrals to be transformed
       !> are input via the type-bound pointer ao_integral_storage.
       procedure :: two_electron_integrals
       !> Variant of "two_electron_integrals" useful for B-spline-only continuum basis. The subroutine makes use of sparsity of
       !> the two-electron integral matrix.
       procedure :: two_electron_integrals_sparse
       !> Transforms two of four atomic orbital indices to molecular ones.
       procedure :: transform_two_indices
       !> Retrieves a block of atomic integrals from the integral array.
       procedure :: fetch_atomic_integrals_block
       !> Moves temporary integral arrays to this object's data arrays and stores the integrals to disk.
       procedure :: finalize_two_electron_integrals
       !> Calculates indices for 1- or 2-electron integrals given their type and the number of basis function pairs/quartets.
       procedure :: integral_index
       !> Returns the name of the basis set.
       procedure :: get_basis_name
       !> Returns the name of the i-th shell in the basis set.
       procedure :: get_shell_name
       !> Returns the shell data for the i-th shell in the basis set.
       procedure :: get_shell_data
       !> Returns an array containing data for all shells of CGTOs in the basis.
       procedure :: get_all_orbital_sets
       !> Returns the matrix containing coefficients for all orbitals in the basis.
       procedure :: get_orbital_coefficient_matrix
       !> Returns the value of initialized.
       procedure :: is_initialized
       !> Same as for atomic_orbital_basis_obj but here we use the symmetry information for the orbitals.
       procedure :: get_continuum_flags
       !> Returns the number of orbitals in a given symmetry.
       procedure :: get_number_of_orbitals
       !> Returns the index within its own symmetry of a given orbital.
       procedure :: get_index_within_symmetry
       !> Returns the symmetry of a given orbital.
       procedure :: get_orbital_symmetry
       !> Given the pair of numbers num,sym it returns the absolute index of the orbital within the whole orbital set.
       procedure :: get_absolute_index
       !> Calculates orbital amplitudes for all continuum channels.
       procedure :: calculate_amplitudes
       !> Calculates radial charge densities of all orbitals in the basis. If the input value of rmat_radius is > 0.0_cfp then
       !> the continuum functions will be normalized to the R-matrix radius rmat_radius. If rmat_radius .le. 0 then normalization
       !> of the continuum functions will not be done. This is useful in case the orbital set corresponds to the Dyson orbitals obtained from CDENPROP: in this
       !> case the orbital coefficients already include the continuum normalization factors and therefore rmat_radius must be set to < 0.0_cfp.
       procedure :: radial_charge_density => orbital_radial_charge_density
       !> Deletes orbital coefficients with magnitude smaller than thrs_orb_cf.
       procedure :: delete_small_coefficients
       !> Writes to disk the array ijkl_indices and the value ind_ijkl_integral. As usual, only master writes its own array.
       procedure :: write_ijkl_indices
       !> Reads from the disk the array ijkl_indices and the value ind_ijkl_integral. The reading is done by the master task which also perform redistribution to other processes.
       !> If shared-memory MPI is used then each NODE keeps only one copy of the array this%ijkl_indices, otherwise the array is kept by every MPI task.
       procedure :: read_ijkl_indices
       !> Evaluates a given orbital (specified by its absolute index) at a set of points in space. The sign of the orbital
       !> at the corresponding points is output too in a separate array. No normalization of the continuum functions is performed
       !> so make sure you're either using atomic basis set whose functions have been normalized to the required R-matrix radius or 
       !> use orbital coefficients which include the continuum normalization factors (as is the case for the Dyson orbitals produced by CDENPROP).
       procedure :: eval_orbital
   end type molecular_orbital_basis_obj

   !> \class <sym_ortho_io>
   !> This data structure is used for input/output by the method orthogonalize of the molecular_orbital_basis_obj.
   type sym_ortho_io
      !> On input: deletion threshold for linearly dependent orbitals.
      real(kind=cfp) :: del_thrs(1:8) = (/-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp/)
      !> On output: number of functions 'deleted' in each symmetry.
      logical, allocatable :: to_delete(:)
   contains
      !> Checks that the input data is OK. We check: del_thrs.
      procedure :: check => check_sym_ortho_io
   end type sym_ortho_io

 contains
 
   function init(this,n,geometry)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: n
      class(geometry_obj), intent(in) :: geometry
      integer :: init
 
         write(stdout,'("--------->","molecular_orbital_basis_obj:init")')
 
         init = 0
 
         if (this%initialized) then
            init = this%final()
            if (init .ne. 0) then
               call xermsg ('molecular_orbital_basis_obj', 'init', &
                            'Finalization has failed. See molecular_orbital_basis_obj%final for details.', init, 0)
               return
            endif
         endif
 
         if (n < 0) then
            init = 1
            call xermsg('molecular_orbital_basis_obj','init','On input the value of n was out of range.',init,0)
            return
         endif
 
         init = this%symmetry_data%init(geometry)
         if (init .ne. 0) then
            call xermsg ('molecular_orbital_basis_obj', 'init', &
                         'Symmetry initialization failed. See symmetry_obj%init for details.', init, 0)
            return
         endif

         this%pg = this%symmetry_data%get_pg()
         this%no_irr = this%symmetry_data%get_no_irrep(this%pg)

         !Note that in this case we don't try to exit the routine with an error
         !code. On run-time unassociated pointer can cause all kinds of mess and
         !this may make it hard to track the source of the problem if the error
         !codes from this routine are not properly processed.
         if (.not. associated(this % ao_basis)) then
            call xermsg ('molecular_orbital_basis_obj', 'init', &
                         'The ao_basis pointer to the AO basis set has not been associated. Fatal error.', 1, 1)
         end if

         !todo test that the atomic and the molecular symmetry data are the same!!!

         if (n /= this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'init', &
                         'The number n on input must be equal to the number of symmetries.', 2, 1)
         end if

         allocate(this%orbital_data(this%no_irr),this%sets_added(this%no_irr),stat=init)
         if (init .ne. 0) then
            call xermsg('molecular_orbital_basis_obj','init','Memory allocation of this%orbital_data has failed.',init,0)
            return
         endif

         this%sets_added = .false.
         this%init_called = .true.
         this%number_of_functions = 0
         this%number_of_shells = 0

         write(stdout,'("<---------","molecular_orbital_basis_obj:init")')
 
   end function init
 
   function final(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer :: final

      integer :: err
 
         write(stdout,'("--------->","molecular_orbital_basis_obj:final")')

         final = 0
 
         this%number_of_shells = 0
         this%number_of_functions = 0
         if (allocated(this%orbital_data)) deallocate(this%orbital_data)
         if (allocated(this%absolute_to_relative)) deallocate(this%absolute_to_relative)
         if (allocated(this%relative_to_absolute)) deallocate(this%relative_to_absolute)
         if (allocated(this%sets_added)) deallocate(this%sets_added)
         if (this%shared_window_ijkl_indices /= -1) then
            call mpi_memory_deallocate_integer_2dim(this%ijkl_indices,size(this%ijkl_indices),this%shared_window_ijkl_indices)
            this%shared_window_ijkl_indices = -1
         else
            if (associated(this%ijkl_indices)) deallocate(this%ijkl_indices)
         endif
         !this%ijkl_indices => null()
         
         
         this%init_called = .false.

         this%pg = -1
         this%no_irr = 0

         deallocate(this%so2mo_range,this%mo2so_range,this%block_offset,this%sym_offset,this%is_continuum,stat=err)
         if (err .ne. 0) final = 1
        
         this%initialized = .false.
         !WE MUST NOT NULLIFY THE POINTER HERE: IF WE DO IT HERE THEN READING OF
         !THE BASIS WILL FAIL SINCE THE POINTER WILL NOT BE ASSOCIATED WHEN
         !ADDING SETS OF ORBITALS. (FINALIZATION OF THE BASIS IS ALWAYS DONE
         !BEFORE READING ANY BASIS).
         !nullify(this%ao_basis)
 
         write(stdout,'("<---------","molecular_orbital_basis_obj:final")')
 
   end function final
 
   subroutine add_shell(this,shell_data)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(shell_data_obj), intent(inout) :: shell_data

      type(orbital_data_obj), allocatable :: temp_orbital_data(:)
      integer :: err, i

         write(stdout,'("--------->","molecular_orbital_basis_obj:add_shell")')

         if (.not. this % init_called) then
            call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                         'Attempt to call add_shell before calling init.', 1, 1)
         end if
 
         if (this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                         'All orbital sets for which space has been allocated have already been supplied.', 2, 1)
         end if

         call shell_data%normalize
         call shell_data%print

         select type (orbital_set => shell_data)
            type is (orbital_data_obj)

               !No duplicities allowed and symmetries of all orbital sets must be consistent:
               if (this % pg /= orbital_set % point_group) then
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'All sets of orbitals in the basis must belong to the same point group symmetry.', 3, 1)
               end if

               if (orbital_set%irr <= 0 .or. orbital_set%irr > this % no_irr) then
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The IRR of the orbital set on input is out of range.', 4, 1)
               end if

               if (this % sets_added(orbital_set % irr)) then
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The basis set already contains orbitals for this symmetry.', 5, 1)
               end if

               if (orbital_set%number_of_coefficients .ne. this%ao_basis%number_of_functions) then
                  print *,orbital_set%number_of_coefficients,this%ao_basis%number_of_functions
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The number of orbital coefficients is not compatible with the AO basis set.', 6, 1)
               endif

               !Resize the basis of orbitals:
               call move_alloc(this%orbital_data,temp_orbital_data)
               this%number_of_shells = this%number_of_shells + 1
               allocate(this%orbital_data(this%number_of_shells),stat=err)
               if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'add_shell', 'Memory allocation 1 failed.',err, 1)

               !Copy the data for the previous sets of orbitals:
               do i=1,this%number_of_shells-1
                  this%orbital_data(i) = temp_orbital_data(i)
               enddo !i
               if (this%number_of_shells > 1) deallocate(temp_orbital_data)

               !Add the new orbital_set data to the list:
               this%orbital_data(this%number_of_shells) = orbital_set
               write(stdout,'("Orbital set of type orbital_data_obj has been added to the molecular basis.")')
               this%sets_added(orbital_set%irr) = .true.
            class default
               call xermsg ('molecular_orbital_basis_obj', 'add_shell', 'The shell type must be orbital_data_obj.',1, 1)
         end select

         this%number_of_functions = this%number_of_functions + shell_data%number_of_functions
 
         if (this%number_of_shells .eq. this%no_irr) then

            write(stdout,'(/,"Orbitals for all symmetries have been supplied. &
                             &Generating indices and analyzing orbital coefficients...")')

            call this%determine_auxiliary_indices
   
            this%initialized = .true.

         endif
 
         write(stdout,'("<---------","molecular_orbital_basis_obj:add_shell")')
 
   end subroutine add_shell

   subroutine print_molecular_orbital_basis_obj(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: i, n_tgt(this%no_irr), n_cnt(this%no_irr)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'print_molecular_orbital_basis_obj', &
                         'The object has not been initialized or not all orbitals have been read-in.', 1, 1)
         end if

         write(stdout,'(/,"--------->","molecular_orbital_basis_obj:print",/)')

         write(stdout,'("Point-group symmetry identifier: ",i0)') this%pg
         write(stdout,'("Number of irreducible representations: ",i0)') this%no_irr
         write(stdout,'("Number of molecular orbitals in each irreducible representation: ")') 
         write(stdout,'(8(i0,1x))') this%orbital_data(1:this%no_irr)%number_of_functions

         n_tgt = 0
         n_cnt = 0
         do i=1,this%number_of_functions
            if (this%is_continuum(i)) then
               n_cnt(this%absolute_to_relative(2,i)) = n_cnt(this%absolute_to_relative(2,i)) + 1
            else
               n_tgt(this%absolute_to_relative(2,i)) = n_tgt(this%absolute_to_relative(2,i)) + 1
            endif
         enddo !i

         write(stdout,'("Number of target orbitals: ")')
         write(stdout,'(8(i0,1x))') n_tgt(1:this%no_irr)

         write(stdout,'("Number of continuum orbitals: ")')
         write(stdout,'(8(i0,1x))') n_cnt(1:this%no_irr)

         write(stdout,'("Name of the associated AO basis: ",a)') trim(this%ao_basis%get_basis_name())
         write(stdout,'("Number of AO basis functions for each irreducible representation: ")') 
         write(stdout,'(8(i0,1x))') this%orbital_data(1:this%no_irr)%number_of_coefficients

         write(stdout,'(/,"Symmetries and indices of the orbitals:")')

         write(stdout,'("Index within symmetry, Orbital symmetry, Overall index, Is continuum")')
         do i=1,this%number_of_functions
            write(stdout,'(3(i0,1x),1X,l)') this%absolute_to_relative(1:2,i), i, this%is_continuum(i)
         enddo

         write(stdout,'("<---------","done:molecular_orbital_basis_obj:print")')

   end subroutine print_molecular_orbital_basis_obj

   subroutine print_orbitals(this)
      use const, only: orbs_line
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: i, j, k, m, mx, symmetry
      real(kind=cfp) :: cf_tmp(orbs_line)
!      real(kind=cfp), allocatable :: cf(:,:)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'print_orbitals', &
                         'The object has not been initialized or not all orbitals have been read-in.', 1, 1)
         end if

         write(stdout,'(/,"--------->","molecular_orbital_basis_obj:print_orbitals",/)')

         write(stdout,'(/,"so2mo_range:")')
         mx = this%ao_basis%number_of_functions
         do i=1,mx
            write(stdout,'("AO ",i0," MO start, end: ",i0,1x,i0)') i, this%so2mo_range(1:2,i)
         enddo

         write(stdout,'(/,"mo2so_range:")')
         do i=1,this%number_of_functions
            write(stdout,'("MO ",i0," AO start, end: ",i0,1x,i0)') i, this%mo2so_range(1:2,i)
         enddo

         write(stdout,'(/,10X,"Orbital coefficients follow")')

!         call this%get_orbital_coefficient_matrix(cf)
!         write(stdout,'(/,10X,"Orbital coefficients follow")')
!         k = 0
!         do i=1,this%number_of_functions/orbs_line
!            write(stdout,'(/,10X,50(i,2X))') (j,j=k+1,k+orbs_line)
!            do j=1,mx
!               write(stdout,'(i,50e25.15)') j, cf(j,k+1:k+orbs_line)
!            enddo
!            k = k + orbs_line
!         enddo
!
!         m = mod(this%number_of_functions,orbs_line)
!         if (m > 0) then
!            write(stdout,'(/,10X,50(i,2X))') (j,j=k+1,k+m)
!            do j=1,mx
!               write(stdout,'(i,50e25.15)') j, cf(j,k+1:k+m)
!            enddo
!         endif

         do symmetry=1,this%no_irr
            write(stdout,'(/,10X,"Symmetry: ",i4)') symmetry
            k = 0
            do i=1,this%orbital_data(symmetry)%number_of_functions/orbs_line
               write(stdout,'(/,10X,50(i0,2X))') (this%relative_to_absolute(j,symmetry),j=k+1,k+orbs_line)
               do j=1,mx
                  cf_tmp(1:orbs_line) = this%orbital_data(symmetry)%coefficients(j,k+1:k+orbs_line)
                  write(stdout,'(i0,50e25.15)') j, cf_tmp(1:orbs_line) !this%orbital_data(symmetry)%coefficients(j,k+1:k+orbs_line)
               enddo
               k = k + orbs_line
            enddo
   
            m = mod(this%orbital_data(symmetry)%number_of_functions,orbs_line)
            if (m > 0) then
               write(stdout,'(/,10X,50(i0,2X))') (this%relative_to_absolute(j,symmetry),j=k+1,k+m)
               do j=1,mx
                  cf_tmp(1:m) = this%orbital_data(symmetry)%coefficients(j,k+1:k+m)
                  write(stdout,'(i0,50e25.15)') j, cf_tmp(1:m) !this%orbital_data(symmetry)%coefficients(j,k+1:k+m)
               enddo
            endif
         enddo !symmetry

         write(stdout,'("<---------","done:molecular_orbital_basis_obj:print_orbitals")')

   end subroutine print_orbitals

   subroutine print_energy_sorted_orbital_table(this)
      use common_obj, only: print_orbital_table
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: i, j, n, err, n_tgt, num
      integer, allocatable :: num_sym(:,:)
      real(kind=cfp), allocatable :: energies(:)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'print_energy_sorted_orbital_table', &
                         'The object has not been initialized or not all orbitals have been read-in.', 1, 1)
         end if

         write(stdout,'(/,"--------->","molecular_orbital_basis_obj:print_energy_sorted_orbital_table",/)')

         write(stdout,'(/,10X,"Continuum orbitals will be ignored.")')

         n_tgt = this%number_of_functions - count(this%is_continuum)

         allocate(energies(n_tgt),num_sym(2,n_tgt),stat=err)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'print_energy_sorted_orbital_table', 'Memory allocation failed.', err, 1)
         end if

         n = 0
         do i=1,this%number_of_shells !over all symmetries
            do j=1,this%orbital_data(i)%number_of_functions
               num = this%relative_to_absolute(j,i)
               if (this%is_continuum(num)) cycle
               n = n + 1
               energies(n) = this%orbital_data(i)%energy(j)
               num_sym(1,n) = j
               num_sym(2,n) = i
            enddo !i
         enddo

         call print_orbital_table(energies,num_sym,n_tgt,this%number_of_shells,.true.)

         write(stdout,'("<---------","done:molecular_orbital_basis_obj:print_energy_sorted_orbital_table")')

   end subroutine print_energy_sorted_orbital_table
 
   subroutine one_electron_integrals(this,integral_storage,integral_options)
      use const
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(integral_options_obj), intent(in) :: integral_options
      class(integral_storage_obj), intent(inout) :: integral_storage
 
      !Input/output of the calculated integrals:
      type(integral_storage_obj) :: ao_integrals_disk
      type(integral_options_obj) :: ao_int_opt
      integer :: number_of_integrals, lunit, first_record, current_pos, last_record, d1, d2, no_blocks
      integer :: p, q, i, j, ij, err, no_ao, no_mo, int_type
      integer, allocatable :: int_index(:,:), ind(:)
      real(kind=cfp), allocatable :: cf(:,:), cf_t(:,:), iq(:,:), iq_t(:,:), ao_int(:)
      real(kind=cfp) :: mo_int
      type(p2d_array_obj), target :: integral_src, integral_tgt !we really need two of these in case disk-to-disk AO-MO run is required
      type(p2d_array_obj), pointer :: ao_integrals, mo_integrals
      logical, parameter :: ao_is_local = .true.
      integer, parameter :: number_of_blocks = 0
      character(len=line_len), allocatable :: column_descriptor(:)
      character(len=line_len) :: ao_header, mo_header
 
         call mpi_mod_barrier(err)
 
         write(stdout,'("--------->","molecular_orbital_basis_obj:one_electron_integrals")')
 
         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', 'The basis set has not been initialized.', 1, 1)
         end if

         if (this % number_of_functions == 0) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                         'Number of molecular orbitals on input is zero.', 2, 1)
         end if

         if (.not.associated(this%ao_integral_storage) .or. .not.associated(this%ao_basis)) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                         'On input at least one of this%ao_integral_storage, this%ao_basis have not been associated.', 4, 1)
         endif

         !Header for the AO integrals that we're looking for
         ao_header = this%ao_integral_storage%contruct_header_string(this%ao_basis%get_basis_name(),one_electron_ints)

         !Header for the MO integrals that will be calculated
         mo_header = integral_storage%contruct_header_string(this%get_basis_name(),one_p_sym_ints)

         !In this section we associate ao_integrals which contains the input AO integrals with the appropriate source.
         !In case the AO integrals are stored in memory then we point directly to the array
         !holding them. If the AO integrals are on disk then we load them into the local array 'integral' and set the pointer ao_integrals to that.
         !At the moment 1p integral transform using SHARED input AO integrals is not supported.
         if (this%ao_integral_storage%in_memory()) then
            ao_integrals => this%ao_integral_storage%integral
            if (this%ao_integral_storage%data_header%name .ne. ao_header) then
               call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                            'The AO integrals on input are not compatible with the AO basis set for the MO basis set.', 5, 1)
            endif
            write(stdout,'("AO integrals with header: ",a)') ao_header
            if (ao_integrals % have_offsets()) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'The AO integrals on input must be LOCAL, i.e. not shared accross processes.', 6, 1)
            end if
         endif

         !load the AO integrals into memory as LOCAL arrays (if they are stored on disk)
         if (this%ao_integral_storage%on_disk()) then
            write(stdout,'("Loading AO integrals from the disk...")')

            err = ao_integrals_disk%init(memory=integral_src)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', 'Memory allocation 3 failed.', err, 1)
            end if

            call ao_integrals_disk%read(this%ao_integral_storage,ao_header,ao_int_opt,ao_is_local)
            ao_integrals => ao_integrals_disk%integral !this points to the local array 'integral_src'
            write(stdout,'("AO integrals with header: ",a)') ao_header
            write(stdout,'("...done")')
         endif

         !BEYOND THIS POINT ao_integrals POINTS TO AN ARRAY CONTAINING THE AO INTEGRALS TO BE TRANSFORMED

         call ao_integrals%get_array_dim(d1,d2,no_blocks) !This gives the dimensions of the AO integrals array
         call ao_integrals%get_column_descriptor(column_descriptor)

         write(stdout,'("On input there is ",i0," types of AO integrals")') d2
         write(stdout,'("Number of AO integrals of each type: ",i0)') d1

         !note that instead of loading the AO integrals we can calculate them now since we have the pointer to the AO basis set (this%ao_basis) and the AO integral routine...

         number_of_integrals = this%number_of_functions*(this%number_of_functions+1)/2 !total number of integrals to calculate; note that this does not include symmetry reduction that we can acheive if we know the symmetry of the 1p operator.

         !Allocate the output arrays if we request the output to be stored in memory.
         if (integral_storage%in_memory()) then
            integral_storage%data_header%name = mo_header

            mo_integrals => integral_storage%integral
            !We allocate space for a non-indexed (that is purely local) array with d2 columns and number_of_integrals rows.
            !The columns in the mo_integrals array correspond to the types of AO integrals we have on input.
            err = mo_integrals%init(number_of_integrals,d2,number_of_blocks,column_descriptor)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'Array initialization failed; see p2d_array_obj%init.', err, 1)
            end if
         endif

         !If we request the output to be stored on disk then start a new record on the data file that will contain the integrals.
         !We also allocate temporary storage for the transformed integrals.
         if (integral_storage%on_disk()) then

            !temporary storage for the integrals
            mo_integrals => integral_tgt
            !We allocate space for a non-indexed (that is purely local) array with d2 columns and number_of_integrals rows.
            !The columns in the mo_integrals array correspond to the types of AO integrals we have on input.
            err = mo_integrals%init(number_of_integrals,d2,number_of_blocks,column_descriptor)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'Array initialization 2 failed; see p2d_array_obj%init.', err, 1)
            end if

            lunit = integral_storage%integral_file%get_unit_no()             !unit that is associated to the file opened
            first_record = integral_storage%integral_file%start_record(mo_header) !position, within the data file, of the first record available for the integral data
         endif

         !BEYOND THIS POINT mo_integrals POINTS TO AN ARRAY CONTAINING THE TRANSFORMED MO INTEGRALS

         !1-PARTICLE INTEGRAL TRANSFORM STARTS HERE: the AO integrals are accessed through the pointer ao_integrals; the MO integrals are accessed through the pointer mo_integrals

         no_ao = this%ao_basis%number_of_functions !total number of AOs
         no_mo = this%number_of_functions !total number of MOs

         allocate(cf_t(no_mo,no_ao),int_index(1:2,no_ao),ind(no_ao),iq(this%number_of_functions,no_ao), &
                  iq_t(no_ao,this%number_of_functions),ao_int(no_ao),stat=err)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', 'Memory allocation 4 failed.', err, 1)
         end if

         !Copy the orbital coefficients to one array: this relies on the fact that the molecular orbitals are indexed symmetry by symmetry.
         call this%get_orbital_coefficient_matrix(cf)

         !transpose the MO coefficient matrix: we use it in the first step where we iterate over the MOs.
         cf_t = transpose(cf)

         !iterate over all types of AO integrals
         do int_type=1,d2

            write(stdout,'("Transforming AO integral type: ",a," ...")') adjustl(column_descriptor(int_type))

            !iterate over all AO integrals (p|O|q) and accumulate their contributions to (i|O|q) where i is MO.
            iq(:,:) = 0.0_cfp
            do p=1,no_ao
               !construct the list of AO indices corresponding to all unique pairs (pq) of the AOs.
               do q=1,p
                  int_index(1,q) = p
                  int_index(2,q) = q
               enddo !q
   
               !calculate indices of the AO integrals corresponding to all unique pairs (pq) of the AOs.
               ind(1:p) = this%ao_basis%integral_index(column_descriptor(int_type),int_index(1:2,1:p))

               !load the AO integrals (p|O|q) with GLOBAL indices 'ind' into ao_int.
               do j=1,p
                  ao_int(j) = ao_integrals%a(ind(j),int_type)
               enddo

               !transform the first index: p->i 
               do j=1,p
                  q = int_index(2,j)
                  if (ao_int(j) .ne. 0.0_cfp) then
                     if (p .eq. q) then
                        do i = this%so2mo_range(1,p), this%so2mo_range(2,p)
                           iq(i,q) = iq(i,q) + ao_int(j)*cf_t(i,p)
                        enddo !i
                     else !we assume the integral is symmetric: (p|O|q) = (q|O|p) so we calculate at once contributions of (p|O|q) to iq(i,q) and iq(i,p)
                        do i = this%so2mo_range(1,p), this%so2mo_range(2,p)
                           iq(i,q) = iq(i,q) + ao_int(j)*cf_t(i,p)
                        enddo !i
                        do i = this%so2mo_range(1,q), this%so2mo_range(2,q)
                           iq(i,p) = iq(i,p) + ao_int(j)*cf_t(i,q)
                        enddo !i
                     endif
                  endif
               enddo !j
            enddo !p

            iq_t = transpose(iq)

            !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,ij,mo_int,q) SHARED(this,iq_t,cf,integral_options,mo_integrals,int_type)
            !$OMP DO SCHEDULE(DYNAMIC)
            do i=1,this%number_of_functions
               do j=1,i
                  !the ij indexing is equivalent to the index function
                  ij = i*(i-1)/2+j

                  mo_int = sum(iq_t(this%mo2so_range(1,j):this%mo2so_range(2,j),i) &
                               * cf(this%mo2so_range(1,j):this%mo2so_range(2,j),j))

                  !The line above is the equivalent of:
                  !mo_int = 0.0_cfp
                  !do q=this%mo2so_range(1,j), this%mo2so_range(2,j)
                  !   mo_int = mo_int + iq_t(q,i)*cf(q,j)
                  !enddo !q

                  if (abs(mo_int) < integral_options%tol) then
                     mo_int = 0.0_cfp
                  else
                     mo_integrals%a(ij,int_type) = mo_int
                  endif

               enddo !j
            enddo !i
            !$OMP END DO
            !$OMP END PARALLEL

            write(stdout,'("...done")')

         enddo !int_type

         deallocate(cf_t,int_index,iq,ind)

         nullify(ao_integrals)
         !FROM NOW ON ao_integrals is nullified

         !If requested print the non-zero integrals
         if (integral_options%print_integrals) then
            call mo_integrals%print(.true.)
         endif

         !dump all integrals to disk and close the record
         if (integral_storage%on_disk()) then

            write(stdout,'("Saving integrals to disk...")')

            !The first record are the integral options.
            call integral_options%write(lunit,first_record,current_pos)

            ! The second record are the ordered integrals: we write them using a method whose choice
            ! depends on whether each process keeps the copy of the full integral array or not.
            i = master
            call mo_integrals%write(lunit,current_pos,last_record,i)
         
            !Every process closes the record so that they all keep identical header information.
            call integral_storage%integral_file%close_record(mo_header,first_record,last_record)
               
            err = mo_integrals%final()
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'Deallocation of the temporary integral array failed.', 5, 1)
            end if

            nullify(mo_integrals)
            !FROM NOW ON mo_integrals is nullified

            write(stdout,'("...done")')

         else
 
            nullify(mo_integrals)
            !FROM NOW ON mo_integrals is nullified

         endif
 
         write(stdout,'("<---------","molecular_orbital_basis_obj:one_electron_integrals")')
 
         call mpi_mod_barrier(err)
 
   end subroutine one_electron_integrals


    !> \brief   Transform atomic 2-electron integrals to molecular ones
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Calculates the molecular-orbital 2-electron integrals (ij|kl) from the known atomic-orbital 2-electron integrals [pq|rs].
    !> The algorithm proceeds in several steps:
    !>  1. Expand all atomic integrals to a working array, including all those redundant in the sense of index symmetries.
    !>  2. Sort the integral array so that all integrals with common p,q are clustered together, resulting in a sequence of
    !>     sparse matrices whose elements are indexed by atomic indices r,s.
    !>  3. Transform both indices of these sparse matrices to molecular ones by application of the expansion coefficient matrices,
    !>     yielding object [pq|kl).
    !>  4. Sort the integral array again, but now group together elements with the same k,l, resulting in a sequence of sparse
    !>     matrices whose elements are indexed by atomic indices p,q. Symbolically: (kl|pq].
    !>  5. Transform both indices p,q by application of the coefficient matrices, yielding (kl|ij).
    !>  6. Finally, sort the integrals to satisfy expectations of the library.
    !>
    !> In the present implementation, steps 1, 2 and 3 are fused to avoid large memory use; only a few blocks are expanded
    !> at a time.
    !>
    !> Some work can be saved with the knowledge of index and orbital symmetries:
    !>  - In the step 1 discard CCCC and CCCT types if not needed.
    !>  - In the step 3 it is possible to calculate the numbers [pq|kl) just for k >= l, due to symmetry k <-> l. In addition
    !>    to this, the CCCC and CCCT combinations can be ignored if two electrons in continuum are not required.
    !>  - In the step 5 it is possible to used the same, ignoring anything else than i >= j, and one can skip also all i,j
    !>    pairs, whose combined orbital symmetry is different than the combined orbital symmetry of the pair k,l. This is
    !>    due to the fact that the product of the two pairs i,j and k,l must be totally symmetric. Furthermore, one can
    !>    skip all i,j pairs that violate [ij] >= [kl], due to the symmetry i,j <-> k,l. Here [ij] is the standard triangular
    !>    multi-index function. Again, as before, CCCC and CCCT combination can be skipped.
    !>
    !> Note that this subroutine uses rectangular indexing function and transforms the indices to the triangular ones only
    !> at the very end.
    !>
    subroutine two_electron_integrals_sparse(this, integral_storage, integral_options)

        use omp_lib, only: omp_get_wtime
        use sort,    only: heap_sort_int_float

        class(molecular_orbital_basis_obj)         :: this
        class(integral_options_obj), intent(in)    :: integral_options
        class(integral_storage_obj), intent(inout) :: integral_storage

        type(p2d_array_obj),     pointer     :: ao_integrals
        character(len=line_len), allocatable :: column_descriptor(:)

        real(kind=cfp), allocatable :: cf(:,:), cft(:,:), Rv(:,:), Cv(:)
        integer,        allocatable :: Cp(:), Cj(:), Ri(:,:)

        real(kind=cfp) :: tolerance
        real(kind=wp)  :: start_t, end_t
        integer        :: d1, d2, no_blocks_ao, no_ao, no_mo, Rn, ierr, i, j
        logical        :: skip2ec

        write(stdout,'("--------->","molecular_orbital_basis_obj:two_electron_integrals_sparse")')

        if (nprocs > 1) then
           call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals_sparse', &
                        'Parallelization using MPI not implemented yet.', 1, 1)
        endif

        start_t = omp_get_wtime()

        ao_integrals => this % ao_integral_storage % integral

        call ao_integrals % get_column_descriptor(column_descriptor)
        call ao_integrals % get_array_dim(d1, d2, no_blocks_ao)

        no_ao = this % ao_basis % number_of_functions  ! total number of AOs
        no_mo = this % number_of_functions             ! total number of MOs

        write(stdout, '("On input there is ",I0," types of AO integrals")') d2
        write(stdout, '("Total number of AO integrals of each type: ",I0)') d1
        write(stdout, '("Every process keeps the full set of AO integrals.")')
        write(stdout, '("Number of atomic orbitals: ",I0)') no_ao
        write(stdout, '("Number of molecular orbitals: ",I0)') no_mo

        ! construct coefficient matrix C(no_ao,no_mo)
        call this % get_orbital_coefficient_matrix(cf)
        call sparse_from_dense(Cp, Cj, Cv, cf)
        write(stdout, '("Coefficient matrix density: ",I0," %")') (Cp(no_mo + 1) - 1) * 100 / (no_ao * no_mo)

        ! transform the two summation indices r,s to k,l
        tolerance = 0
        skip2ec = (this % ao_basis % n_cont_fns > 0) .and. (.not. integral_options % two_p_continuum)

        call this % transform_two_indices(no_ao, no_ao, no_mo, Rn, Ri, Rv, Cp, Cj, Cv, &
                                          1, tolerance, skip2ec)

        ! invert indexing (ie. make p,q the summation indices)
        call invert_indexing(Rn, Ri, no_ao, no_mo)

        ! sort the integral array (cluster by k,l)
        call sort_intermediate_integrals(Rn, Ri, Rv, no_ao * no_ao)

        ! transform the two summation indices p,q to i,j
        tolerance = integral_options % tol
        call this % transform_two_indices(no_mo, no_ao, no_mo, Rn, Ri, Rv, Cp, Cj, Cv, &
                                          2, tolerance, skip2ec)

        ! reindex the array to the expected triangular format
        call rect_index_to_tri_index(Rn, Ri, no_mo, 1)

        ! sort the integral array in accord with the triangular indexing
        call heap_sort_int_float(Rn, 1, Ri, Rv)

        ! store integrals to disk and finalize the calculation process
        call this % finalize_two_electron_integrals(Rn, Ri, Rv, integral_storage, integral_options, column_descriptor)

        end_t = omp_get_wtime()

        write(stdout,'("two_electron_integrals_sparse took [s]: ",F25.15)') end_t - start_t
        write(stdout,'("<---------","molecular_orbital_basis_obj:two_electron_integrals_sparse")')

    end subroutine two_electron_integrals_sparse


    !> \brief   Copy atomic 2-electron integrals to a sparse array
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Retrieve all [pq|rs] atomic integrals for given block index (ie. for given p,q). Write them into
    !> the supplied (pre-allocated) arrays Ri(:), Rv(:) as elements of sparse 4-index tensor; Ri(:) will
    !> contains rectangular zero-based multi-indices of the non-zero elements, Rv(:) will contain the
    !> values of the integrals.
    !>
    !> \todo Do not repeatedly recalculate `no_T_ao`.
    !>
    !> \param this  Reference to the parent type.
    !> \param iblk  Zero-based rectangular multi-index of the block to retrive.
    !> \param rbeg  Starting index for arrays Ri, Rv.
    !> \param Rn    On return, number of elements written into the arrays Rv and Ri.
    !> \param Ri    Zero-based multi-index for each element in Rv.
    !> \param Rv    Array of non-zero atomic two-electron integrals.
    !> \param skip2ec  Whether to skip [CC|CT] and [CC|CC] integrals (mostly .TRUE.).
    !>
    subroutine fetch_atomic_integrals_block(this, iblk, rbeg, Rn, Ri, Rv, skip2ec)

        class(molecular_orbital_basis_obj) :: this

        integer,                    intent(in)    :: iblk, rbeg
        integer,                    intent(out)   :: Rn
        real(kind=cfp), allocatable               :: Rv(:,:) !intent(inout)
        integer,        allocatable               :: Ri(:,:) !intent(inout)
        logical,                    intent(in)    :: skip2ec

        type(p2d_array_obj), pointer :: ao_integrals

        integer :: no_ao, no_T_ao, last_CT_fn, n_prec_ints, n_TT_pairs, p, q, r, s, u, v, pq, pqm, rs, rsm, I
        integer :: ptype, qtype, rtype, stype
        logical :: two_p_continuum

        ao_integrals => this % ao_integral_storage % integral

        no_ao       = this % ao_basis % number_of_functions
        n_prec_ints = this % ao_basis % n_prec_ints
        n_TT_pairs  = this % ao_basis % n_TT_pairs
        last_CT_fn  = this % ao_basis % last_CT_fn

        ! find the last atomic T orbital
        no_T_ao = 0
        do u = 1, this % ao_basis % number_of_functions
            v = this % ao_basis % indices_to_shells(1, u)
            if (this % ao_basis % shell_descriptor(3, v) == 0) then
                no_T_ao = u
            end if
        end do

        Rn = 0
        p = 1 + iblk / no_ao
        q = 1 + mod(iblk, no_ao)

        ptype = merge(0, 1, p <= no_T_ao)  ! T (= 0) or C (= 1)
        qtype = merge(0, 1, q <= no_T_ao)  ! T (= 0) or C (= 1)

        pq = max(p,q) * (max(p,q) - 1) / 2 + min(p,q)

        do r = 1, no_ao

            rtype = merge(0, 1, r <= no_T_ao)  ! T (= 0) or C (= 1)

            ! skip CCCT and CCCC combinations
            if (skip2ec .and. ptype + qtype + rtype == 3) then
                exit
            end if

            do s = 1, no_ao

                stype = merge(0, 1, s <= no_T_ao)  ! T (= 0) or C (= 1)
                rs = max(r,s) * (max(r,s) - 1) / 2 + min(r,s)

                ! skip TCCC, CTCC and CCTC combinations
                if (skip2ec .and. ptype + qtype + rtype + stype == 3) then
                    exit
                end if

                ! get mapped indices
                rsm = this % ao_basis % ordered_pairs(1, rs)
                pqm = this % ao_basis % ordered_pairs(1, pq)

                ! compactify the overall non-redundant multi-index
                u = max(pqm, rsm) ; v = min(pqm, rsm) ; I = u * (u - 1) / 2 + v

                ! special indexing for CCTT class when 2p continuum is disabled
                if (skip2ec) then
                    if (v <= n_TT_pairs .and. u > last_CT_fn) I = n_prec_ints + v + n_TT_pairs * (u - last_CT_fn - 1)
                    if (u <= n_TT_pairs .and. v > last_CT_fn) I = n_prec_ints + u + n_TT_pairs * (v - last_CT_fn - 1)
                end if

                ! sanity check - the index must not exceed the size of the atomic integral storage
                if (I > size(ao_integrals % a, 1)) then
                    write (stdout, '("Inconsistency: Index overflow!")')
                    write (stdout, '("- Number of AOs: ",I0)') no_ao
                    write (stdout, '("- n_prec_ints: ",I0)') n_prec_ints
                    write (stdout, '("- n_TT_pairs: ",I0)') n_TT_pairs
                    write (stdout, '("- last_CT_fn: ",I0)') last_CT_fn
                    write (stdout, '("- p, q, r, s: ",I0,1x,I0,1x,I0,1x,I0)') p, q, r, s
                    write (stdout, '("- u, v, I: ",I0,1x,I0,1x,I0)') u, v, I
                    flush (stdout)
                    stop 1
                end if

                ! store this integral into the sparse array
                if (ao_integrals % a(I, 1) /= 0) then
                    Rn = Rn + 1
                    Ri(Rn + rbeg - 1, 1) = (((p-1) * no_ao + (q-1)) * no_ao + (r-1)) * no_ao + (s-1)  ! zero-based multi-index
                    Rv(Rn + rbeg - 1, 1) = ao_integrals % a(I, 1)
                end if

            end do ! s

        end do  ! r

    end subroutine fetch_atomic_integrals_block


    !> \brief   Extract non-zero elements from a dense matrix
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Given a dense matrix A, copy its non-zero elements into the array Sv, with corresponding CSC indices in Sp and Sj.
    !>
    !> \param Sp  Index of the first element in Sj/Sv in the given column. If two consecutive elements in Sp have the same
    !>            value, it means that the column has zero length.
    !> \param Sj  Row index of the corresponding element in Sv.
    !> \param Sv  Non-zero elements of A.
    !> \param A   Dense matrix (column major storage).
    !>
    subroutine sparse_from_dense(Sp, Sj, Sv, A)

        real(kind=cfp), allocatable, intent(in)  :: A(:,:)
        real(kind=cfp), allocatable, intent(out) :: Sv(:)
        integer,        allocatable, intent(out) :: Sp(:), Sj(:)

        integer :: row, col, ierr, m, n, pos, nnz

        m   = size(A, 1)     ! leading dimension of the matrix
        n   = size(A, 2)     ! other dimension
        nnz = count(A /= 0)  ! number of non-zero elements

        allocate(Sp(n + 1), Sj(nnz), Sv(nnz), stat = ierr)
        if (ierr /= 0) call xermsg ('molecular_orbital_basis_obj', 'sparse_from_dense', 'Memory allocation failure.', 1, 1)

        pos = 0
        Sp(1) = 1

        do col = 1, n
            do row = 1, m
                if (A(row,col) /= 0) then
                    pos = pos + 1
                    Sj(pos) = row
                    Sv(pos) = A(row,col)
                end if
            end do
            Sp(col + 1) = pos + 1
        end do

    end subroutine sparse_from_dense


    !> \brief   Transform both indices of each block of a large sparse matrix by the (sparse) coefficient matrix
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> In the first step of the transformation (ie. [pq|rs] -> [pq|kl)), the input arrays Ri, Rv are ignored and the sparse
    !> blocks Rb of the 4-index tensor R are retrieved directly from the atomic integrals storage. Ri, Rv are then allocated
    !> and filled with the transformed data.
    !>
    !> In the second step, the sparse blocks Rb are truly read from the input arrays Ri, Rv. These are then destroyed and
    !> reallocated to the proper size before returning the transformed data.
    !>
    !> Transformation of every block amonts to two sparse matrix multiplications (with some restrictions on non-zero
    !> elements). The result of each such transformation is a subset of a dense block. Each transformed block is stored in
    !> a separate array in a linked list data structure to avoid frequent reallocations of a growing array of elements.
    !> These are merged to Ri, Rv at the very end of the transformation step.
    !>
    !> \todo The merging could be parallelized, too, if each block also stored its offset.
    !>
    !> \param Nbk Number of blocks along the diagonal (= nAO in the first step, nMO in the second step).
    !> \param nAO Number of atomic orbitals (dimension of the index to transform).
    !> \param nMO Number of molecular orbitals (dimension of the index after transformation).
    !> \param Rn  Currently used portion of the arrays Ri, Rv.
    !> \param Ri  Zero-based multiindex for each element in Rv.
    !> \param Rv  Structurally non-zero elements of the tensor R (integral storage).
    !> \param Cp  Column pointers in Ci (ie. where the coefficients corresponding to i-th MO start in Cv).
    !> \param Cj  Row indices corresponding to elements in Cv.
    !> \param Cv  Structurally non-zero elements of C (the MO <- AO coefficient matrix).
    !> \param step       First or second step of the transformation.
    !> \param tolerance  Minimal absolute value of integral to keep it.
    !> \param skip2ec    Whether to skip CCCT and CCCC combinations.
    !>
    subroutine transform_two_indices(this, Nbk, nAO, nMO, Rn, Ri, Rv, Cp, Cj, Cv, step, tolerance, skip2ec)

        use const,   only: abel_prod_tab
        use omp_lib, only: omp_get_num_threads, omp_get_thread_num, omp_get_wtime

        ! helper structure to hold transformed blocks without the need of a frequent re-allocation
        type :: LinkedList
            integer,        allocatable :: Wi(:,:)
            real(kind=cfp), allocatable :: Wv(:,:)
            type(LinkedList),   pointer :: next => null()
            type(LinkedList),   pointer :: back => null()
        end type LinkedList

        class(molecular_orbital_basis_obj) :: this

        real(kind=cfp), allocatable, intent(inout) :: Rv(:,:)
        real(kind=cfp), allocatable, intent(in)    :: Cv(:)
        integer,        allocatable, intent(inout) :: Ri(:,:)
        integer,        allocatable, intent(in)    :: Cp(:), Cj(:)
        integer,                     intent(in)    :: Nbk, nAO, nMO, step
        integer,                     intent(inout) :: Rn
        logical,                     intent(in)    :: skip2ec
        real(kind=cfp),              intent(in)    :: tolerance

        real(wp) :: start_t, end_t
        real(kind=cfp) :: x
        type(LinkedList),   pointer :: buffptr => null(), myptr => null()
        real(kind=cfp), allocatable :: Rv_out(:), Vv(:,:), Wv(:,:)
        integer,        allocatable :: Ri_out(:), Vp(:), Vj(:,:), Wq(:), Wj(:,:), Rp(:), Rj(:,:), blocks(:)
        integer :: a, b, c, d, ab, cd, i, j, iblk, nblk, mx, nnz, nonzeros, discarded, ithread, nthreads, ierr
        integer :: rbeg, rcur, rend, blocksym, blockcnt, TA_orb(9), CA_orb(9), TM_orb(9), CM_orb(9), iT, iC, iRp, orb

        start_t = omp_get_wtime()

        blockcnt = 0    ! combined CC/CT/TT type of the orbitals with the block indices
        blocksym = 0    ! combined symmetry of the orbitals with the block indices (only used in the second stage)

        ! Get starting (and first continuum) molecular orbitals of individual symmetries
        !
        ! - For example, when the number of target+continuum orbitals per symmetry is 5+4,0+3,1+3,0+0,0+3,0+0,0+0,0+0,
        !   so in total 19 orbitals, then the two helper offset arrays will look like this:
        !
        !      TM_orb = 1  10  13  17  17  20  20  20  20
        !      CM_orb = 6  10  14  17  17  20  20  20  20

        iT = 1 ; TM_orb = this % number_of_functions + 1
        iC = 1 ; CM_orb = this % number_of_functions + 1
        sym_loop: do i = 1, this % no_irr
            orb_loop: do j = 1, this % orbital_data(i) % number_of_functions
                orb = this % relative_to_absolute(j, i)
                ! first orbital of each symmetry (if there is any) will be stored in TM_orb for that symmetry and
                ! also all previous empty symmetries; it will be also stored in CM_orb for all previous empty symmetries
                if (j == 1) then
                    do while (iT <= i) ; TM_orb(iT) = orb ; iT = iT + 1 ; end do
                    do while (iC <  i) ; CM_orb(iC) = orb ; iC = iC + 1 ; end do
                end if
                ! first continuum orbital of each symmetry (if there is any) will be stored in CM_orb
                if (this % is_continuum(orb)) then ; CM_orb(iC) = orb ; iC = iC + 1 ; exit orb_loop ; end if
            end do orb_loop
        end do sym_loop

        ! Now do the same for atomic orbitals.
        !
        ! This is simpler, because the atomic orbitals are not sorted by symmetry; instead
        ! they are stored in one consecutive array with target orbitals at the beginning.

        TA_orb = this % ao_basis % number_of_functions + 1
        CA_orb = this % ao_basis % number_of_functions + 1

        TA_orb(1)  = 1
        do j = 1, this % ao_basis % number_of_functions
            i = this % ao_basis % indices_to_shells(1, j)
            if (this % ao_basis % shell_descriptor(3, i) == 1) CA_orb(1) = min(j, CA_orb(1))
        end do

        write (stdout, '("==================================================================")')
        if (step == 1) then
            write (stdout, '("Transforming the first pair of atomic indices, [AA|AA] -> [AA|MM)")')
        else
            write (stdout, '("Transforming the second pair of atomic indices, (MM|AA] -> (MM|MM)")')
        end if
        write (stdout, '("Atomic orbital offsets: ")')
        write (stdout, '(5x,"T",I6,"   (",I0,")")') TA_orb(1), TA_orb(2)
        write (stdout, '(5x,"C",I6)') CA_orb(1)
        write (stdout, '("Molecular orbital offsets per symmetry: ")')
        write (stdout, '(5x,"T",8I6,"   (",I0,")")') TM_orb(1:8), TM_orb(9)
        write (stdout, '(5x,"C",8I6)') CM_orb(1:8)

        ! First, find out the number of blocks so that we can process them concurrently.
        !
        ! - In the first step (ie. [pq|rs] -> [pq|kl)), there will be at most "nthreads" blocks in memory at a single time; we
        !   allocate space for them here. The total number of blocks is known, too: one for each pair of atomic orbitals.
        ! - In the second step (ie. (kl|pq] -> (kl|ij)), the number of blocks is retrieved from the R-tensor itself, and blocks
        !   offsets in the linear storage arrays are stored for later convenience.

        !$omp parallel shared(nthreads)
        !$omp master
        nthreads = omp_get_num_threads()
        !$omp end master
        !$omp end parallel

        if (step == 1) then
            nblk = nAO * nAO
            allocate (Ri(nthreads * nAO * nAO, 1),stat=ierr)  ! This is enough to hold even dense blocks, ...
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 1 failure.', ierr, 1)
            endif
            allocate (Rv(nthreads * nAO * nAO, 1),stat=ierr)  ! ... so any sparse block should fit in without problems.
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 2 failure.', ierr, 1)
            endif
        else
            allocate (blocks(Nbk * Nbk + 1),stat=ierr)
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 3 failure.', ierr, 1)
            endif
            nblk = 0 ; mx = -1
            do rcur = 1, Rn
                ab = Ri(rcur, 1) / (nAO * nAO)
                if (mx /= ab) then
                    nblk = nblk + 1
                    blocks(nblk) = rcur
                end if
                mx = ab
            end do
            blocks(nblk + 1) = Rn + 1  ! points after the end of Ri
        end if

        ! And here starts the parallel transformation itself. The threads sequentially pick one block at a time and process
        ! it, transforming it from atomic to molecular indices.
        !
        ! - In the first step (ie. [pq|rs] -> [pq|kl)), the blocks to transform are constructed on the fly by the very thread
        !   that is about to do the transformation.
        ! - In the second step (ie. (kl|pq] -> (kl|ij)), the blocks are retrieved from the R-tensor as obtained from the subroutine
        !   arguments.

        discarded = 0
        nonzeros = 0

        !$omp parallel &
        !$omp& shared  (this, Ri, Ri_out, Rv, Rv_out, Nbk, nblk, blocks, step, nAO, nMO, TA_orb, CA_orb, &
        !$omp&          TM_orb, CM_orb, skip2ec, Cp, Cj, Cv, buffptr) &
        !$omp& private (iblk, rbeg, rcur, rend, ab, cd, a, b, c, d, i, j, blocksym, blockcnt, orb, &
        !$omp&          iRp, Rp, Rj, Vp, Vv, Vj, Wq, Wv, Wj, nnz, ithread, myptr, ierr) &
        !$omp& reduction (+ : discarded, nonzeros)

        ithread = omp_get_thread_num()

        ! allocate workspaces (again, use enough space for dense blocks)
        allocate(Rp(nAO + 1), Rj(nAO * nAO, 1), stat=ierr)
        if (ierr /= 0) then
           call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 4 failure.', ierr, 1)
        endif
        allocate(Vp(nMO + 1), Vv(nAO * nMO, 1), Vj(nAO * nMO, 1), stat=ierr)
        if (ierr /= 0) then
           call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 5 failure.', ierr, 1)
        endif
        allocate(Wq(nMO + 1), Wv(nMO * nMO, 1), Wj(nMO * nMO, 1), stat=ierr)
        if (ierr /= 0) then
           call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 6 failure.', ierr, 1)
        endif

        !$omp do schedule (dynamic, 1)

        ! loop over all blocks of the R tensor
        R_loop: do iblk = 1, nblk

            ! get section of data corresponding to this block
            if (step == 1) then
                rbeg = nAO * nAO * ithread + 1
                rend = nAO * nAO * (ithread + 1)
                call this % fetch_atomic_integrals_block(iblk - 1, rbeg, nnz, Ri, Rv, skip2ec)
                rend = rbeg + nnz - 1
            else
                rbeg = blocks(iblk)
                rend = blocks(iblk + 1) - 1
            end if

            ! decode the multi-index
            ab = Ri(rbeg, 1) / (nAO * nAO)
            a = 1 + ab / Nbk
            b = 1 + mod(ab, Nbk)

            ! get CC/CT/TT (= 2/1/0) type of the two (atomic or molecular) orbitals
            blockcnt = 0
            if (step == 1) then
                ! type of atomic orbital pair
                orb = max(a,b) * (max(a,b) - 1) / 2 + min(a,b)
                blockcnt = this % ao_basis % ordered_pairs(2, orb) - 1
            else
                ! type of molecular orbital pair
                if (this % is_continuum(a)) blockcnt = blockcnt + 1
                if (this % is_continuum(b)) blockcnt = blockcnt + 1
            end if

            ! get combined symmetry of the two (molecular) orbitals (only used in the last step)
            if (step == 2) then
                blocksym = abel_prod_tab(count(TM_orb <= a), count(TM_orb <= b))
            end if

            ! unpack multi-index Ri to two CSC split indices Rp, Rj
            iRp = 1
            do rcur = rbeg, rend
                ! decode the multi-index
                cd = mod(Ri(rcur, 1), nAO * nAO)
                c = 1 + cd / nAO
                d = 1 + mod(cd, nAO)

                ! fill all column pointers up to the current column and set the row index for this element
                do while (iRp <= c) ; Rp(iRp) = rcur - rbeg + 1 ; iRp = iRp + 1 ; end do
                Rj(rcur - rbeg + 1, 1) = d
            end do

            ! finalize the column pointer array (point beyond the storage)
            Rp(iRp:nAO + 1) = (rend + 1) - rbeg + 1

            ! multiply C . R -> V, use neither triangular reduction, nor pyramidal reduction
            call transform_one_index(nAO, nMO, Cp, Cj, Cv, &
                                     nAO, nAO, Rp, Rj, rbeg, Rv, &
                                               Vp, Vj, Vv, &
                                    .false., -1, blocksym, blockcnt, &
                                    TM_orb, CM_orb, TA_orb, CA_orb, skip2ec)

            ! multiply C . V -> W, use triangular reduction, and (in second step) also pyramidal reduction
            call transform_one_index(nAO, nMO, Cp, Cj, Cv, &
                                     nAO, nMO, Vp, Vj, 1, Vv, &
                                               Wq, Wj, Wv, &
                                    .true.,  merge(-1, ab, step == 1), blocksym, blockcnt, &
                                    TM_orb, CM_orb, TM_orb, CM_orb, skip2ec)

            ! number of non-zero elements (integrals) in the transformed block
            nnz = Wq(nMO + 1) - 1

            ! reconstruct rectangular multi-indices from CSC indices and trim integrals by the given threshold
            nnz = 0
            do i = 1, nMO
                do j = Wq(i), Wq(i + 1) - 1
                    if (abs(Wv(j, 1)) >= tolerance) then
                        nnz = nnz + 1
                        Wj(nnz, 1) = (ab * nMO + (i - 1)) * nMO + (Wj(j, 1) - 1)
                        Wv(nnz, 1) = Wv(j, 1)
                    end if
                end do
            end do
            discarded = discarded + Wq(nMO + 1) - 1 - nnz
            nonzeros = nonzeros + nnz

            ! add a new element to the auxiliary buffer list in a thread-safe way
            !$omp critical
            if (.not. associated(buffptr)) then
                allocate (buffptr, stat=ierr)
                if (ierr /= 0) then
                   call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 7 failure.', ierr, 1)
                endif
            else
                allocate (buffptr % next, stat=ierr)
                if (ierr /= 0) then
                   call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 8 failure.', ierr, 1)
                endif
                buffptr % next % back => buffptr
                buffptr => buffptr % next
            end if
            myptr => buffptr
            !$omp end critical

            ! copy this thread's data to the auxiliary buffer list
            allocate (myptr % Wi(1:nnz, 1:1), source = Wj(1:nnz, 1:1), stat=ierr)
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 9 failure.', ierr, 1)
            endif

            allocate (myptr % Wv(1:nnz, 1:1), source = Wv(1:nnz, 1:1), stat=ierr)
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 10 failure.', ierr, 1)
            endif

        end do R_loop

        !$omp end parallel

        write (stdout, '("Elements discarded due to threshold: ",I0)') discarded

        ! drop atomic integrals, we will no longer need them (and every bit of free memory is needed below)
        if (step == 1) then
            if (this % ao_integral_storage % integral % final() /= 0) then
                write (stdout, '("WARNING: Finalization of the atomic integral storage failed!")')
            end if
        end if

        ! re-allocate work arrays to the new size
        if (allocated(Ri)) deallocate (Ri) ; allocate (Ri(nonzeros, 1), stat=ierr)
        if (ierr /= 0) then
           call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 11 failure.', ierr, 1)
        endif

        if (allocated(Rv)) deallocate (Rv) ; allocate (Rv(nonzeros, 1), stat=ierr)
        if (ierr /= 0) then
           call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 12 failure.', ierr, 1)
        endif

        ! copy data from the auxiliary buffer list into the work arrays; release buffers once processed
        Rn = 0
        do while (associated(buffptr))
            nnz = size(buffptr % Wv, 1)
            Ri(Rn + 1 : Rn + nnz, 1) = buffptr % Wi(1:nnz, 1)
            Rv(Rn + 1 : Rn + nnz, 1) = buffptr % Wv(1:nnz, 1)
            Rn = Rn + nnz
            if (associated(buffptr % back)) then
                buffptr => buffptr % back
                deallocate (buffptr % next)
            else
                deallocate (buffptr)
                buffptr => null()
            end if
        end do

        write (stdout, '("Non-zero elements on exit: ",I0)') Rn

        end_t = omp_get_wtime()

        write (stdout, '("Time spent in transform_two_indices [s]: ",F25.15,/)') end_t - start_t

    end subroutine transform_two_indices


    !> \brief   Multiplication of two sparse matrices (somewhat tweaked)
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Multiply two sparse matrices, producing a sparse matrix as a result.
    !>
    !> The subroutine allows restricting the operation to some elements to allow making use of the index symmetries
    !> of two-particle integrals:
    !>  - If "triangle" is true, then only a triangular part of the resulting matrix will be computed.
    !>  - If "pyramid" is non-negative, then only multi-indices smaller than the given value will be computed.
    !>
    !> If both the above are true, which indicates the second stage of transformation, ie. (kl|pq] -> (kl|ij), then also:
    !>  - Skip combinations of orbitals with incompatible symmetries. The integral can be only nonzero when the overall symmetry
    !>    of the four orbitals is totally symmetric.
    !>  - If "skip2ec" is true, do not calculate molecular integrals of CCCC and CCCT types.
    !>
    !> \param nAr     Number of rows in A (leading dimension).
    !> \param nAc     Number of cols in A (the other dimension).
    !> \param Ap      Positions in Av of the first element of each column of A.
    !> \param Aj      Row indices corresponding to elements in Av.
    !> \param Av      Non-zero elements of the matrix A.
    !> \param nBr     Number of rows in B (leading dimension).
    !> \param nBc     Number of cols in B (the other dimension).
    !> \param Bp      Positions in Bv of the first element of each column of B.
    !> \param Bj      Row indices corresponding to elements in Bv.
    !> \param Bv_beg  Starting index in array Bv.
    !> \param Bv      Non-zero elements of the matrix B.
    !> \param Cp      On return, Positions in Cv of the first element of each column of C.
    !> \param Cj      Row indices corresponding to elements in Cv.
    !> \param Cv      Non-zero elements of the sparse matrix product.
    !> \param Ci      Zero-based multi-indices corresponding to the elements of Cv.
    !> \param triangle  Calculate only a triangular subset of C.
    !> \param pyramid   Further restriction on calculated elements of C (limit on multi-index).
    !> \param blocksym  Combined symmetry of the other pair of orbitals.
    !> \param blockcnt  CC/CT/TT (= 2/1/0) type of the other pair of orbitals.
    !> \param TA_orb    Helper array with the index of the first orbital per symmetry. Corresponds to column index of A.
    !> \param CA_orb    Helper array with the index of the first continuum orbital per symmetry. Corresponds to column index of A.
    !> \param TB_orb    Helper array with the index of the first orbital per symmetry. Corresponds to column index of B.
    !> \param CB_orb    Helper array with the index of the first continuum orbital per symmetry. Corresponds to column index of B.
    !> \param skip2ec   Whether to skip CCCC and CCCT combinations of molecular orbitals.
    !>
    subroutine transform_one_index (nAr, nAc, Ap, Aj, Av, nBr, nBc, Bp, Bj, Bv_beg, Bv, Cp, Cj, Cv, &
                                    triangle, pyramid, blocksym, blockcnt, TA_orb, CA_orb, TB_orb, CB_orb, skip2ec)

        use const, only: abel_prod_tab

        real(kind=cfp), allocatable :: Av(:), Bv(:,:) !intent(in)
        real(kind=cfp), allocatable :: Cv(:,:) !intent(inout)
        integer, allocatable :: Ap(:), Aj(:), Bp(:), Bj(:,:) !intent(in)
        integer, allocatable :: Cp(:), Cj(:,:)               !intent(inout)
        integer, intent(in)  :: nAr, nAc, nBr, nBc, Bv_beg, pyramid, blocksym, blockcnt
        integer, intent(in)  :: TA_orb(9), CA_orb(9), TB_orb(9), CB_orb(9)
        logical, intent(in)  :: triangle, skip2ec

        real(kind=cfp) :: x
        integer :: apos, acur, aend, arow, acol, bpos, bcur, brow, bend, bcol, cpos
        integer :: asym,  bsym    ! symmetry m-value of the current A/B orbital (column index)
        integer :: atype, btype   ! T/C (= 0/1) type of the current A/B orbital (column index)
        integer :: sym            ! combined orbital symmetry of the block and of the A orbital
        integer :: A_T(9), B_T(9) ! zero-based indices of the first A/B orbital in given symmetry
        integer :: A_C(9), B_C(9) ! zero-based indices of the first A/B continuum orbital in given symmetry
        integer :: iA_T, iB_T     ! current positions in arrays A_T, B_T
        integer :: iA_C, iB_C     ! current positions in arrays A_C, B_C
        logical :: last_step      ! transformation of the last index

        ! sanity check - the multiplication goes in sync along A and B storage (ie. A has to be transposed)
        if (nAr /= nBr) then
            call xermsg ('molecular_orbital_basis_obj', 'sparse_mmul', 'Non-conformant matrices passed as arguments.', 1, 1)
        end if

        ! determine if this is the last step (and column indices of both A and B correspond to molecular orbitals)
        last_step = triangle .and. pyramid >= 0

        ! scan the matrices to obtain important column pointers to speed up jumping through the columns
        call extract_orb_data(TA_orb, CA_orb, nAr, nAc, Ap, A_T, A_C)
        call extract_orb_data(TB_orb, CB_orb, nBr, nBc, Bp, B_T, B_C)

        ! loop through rows of A
        iA_T = 1; iA_C = 1; asym = 0; atype = 0; cpos = 0
        A_loop: do acol = 1, nAc

            ! starting position of this column in C
            Cp(acol) = cpos + 1

            ! first and last element of the column
            apos = Ap(acol)
            aend = Ap(acol + 1) - 1

            ! exit if at end of matrix
            if (apos == Ap(nAc + 1)) exit A_loop

            ! update current orbital symmetry / type
            do while (A_T(iA_T) <= apos .or. A_C(iA_C) <= apos)
                if (iA_C < iA_T) then
                    atype = 1 ; iA_C = iA_C + 1
                else
                    atype = 0 ; iA_T = iA_T + 1 ; asym = asym + 1
                end if
            end do

            ! skip CCCT and CCCC combinations
            if (skip2ec .and. blockcnt + atype >= 3) cycle A_loop

            ! in the last step calculate the combined [kli] orbital symmetry
            if (last_step) sym = abel_prod_tab(blocksym, asym)

            ! loop through elements of B
            iB_T = 1; iB_C = 1; bsym = 0; btype = 0
            B_loop: do bcol = 1, nBc

                ! first and last element of the column
                bpos = Bp(bcol)
                bend = Bp(bcol + 1) - 1

                ! exit if at end of matrix
                if (bpos == Bp(nBc + 1)) exit B_loop

                ! take into account (i.e. skip) the i <-> j, k <-> l and [ij] <-> [kl] symmetry
                if ((triangle .and. bcol > acol) .or. &
                    (pyramid >= 0 .and. (acol - 1) * nBc + (bcol - 1) > pyramid)) then
                    exit B_loop
                end if

                ! update current orbital symmetry and type
                do while (B_T(iB_T) <= bpos .or. B_C(iB_C) <= bpos)
                    if (iB_C < iB_T) then
                        btype = 1 ; iB_C = iB_C + 1
                    else
                        btype = 0 ; iB_T = iB_T + 1 ; bsym = bsym + 1
                    end if
                end do

                ! skip combinations CTCC
                if (skip2ec .and. blockcnt + atype + btype >= 3) cycle B_loop

                ! during the last step skip (non-contributing) orbitals of symmetries other than sym_i * sym_j * sym_k
                if (last_step .and. sym /= bsym) cycle B_loop

                ! calculate dot product of the two sparse columns
                x = 0
                acur = apos
                bcur = bpos
                do while (acur <= aend .and. bcur <= bend)
                    arow = Aj(acur)
                    brow = Bj(bcur, 1)
                    if      (arow < brow) then ; acur = acur + 1
                    else if (arow > brow) then ; bcur = bcur + 1
                    else
                        x = x + Av(acur) * Bv(Bv_beg-1 + bcur, 1)
                        acur = acur + 1
                        bcur = bcur + 1
                    end if
                end do

                ! store the result in the next free cell of C
                if (x /= 0) then
                    cpos = cpos + 1
                    Cj(cpos, 1) = bcol              ! set row index
                    Cv(cpos, 1) = x                 ! set value
                end if

            end do B_loop

        end do A_loop

        ! finalize the C column pointer array
        Cp(nAc + 1) = cpos + 1

    end subroutine transform_one_index


    !> \brief  Get special sparse matrix column pointers
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Constructs two sets of column pointers: One pointing to the first column corresponding to the first molecular
    !> orbital in every symmetry, the other pointing to the first continuum molecular orbital in every symmetry.
    !>
    !> \param T_orb  Absolute index of the first molecular orbital per symmetry.
    !> \param C_orb  Absolute index of the first continuum molecular orbital per symmetry.
    !> \param nAr    Number of rows in Ai.
    !> \param nAc    Number of columns in Ai (must be equal to the total number of molecular orbitals).
    !> \param An     Number of elements in Ai.
    !> \param A_T    On output, positions in Ai of the column corresponding to the first molecular orbital per symmetry.
    !> \param A_C    On output, positions in Ai of the column corresponding to the first continuum molecular orbital per symmetry.
    !>
    subroutine extract_orb_data(T_orb, C_orb, nAr, nAc, Ap, A_T, A_C)

        integer, allocatable :: Ap(:) !intent(in)
        integer, intent(in)  :: T_orb(9), C_orb(9), nAr, nAc
        integer, intent(out) :: A_T(9), A_C(9)

        integer :: orb, pos, iT, iC

        A_T = Ap(nAc + 1)  ! by default point beyond the storage
        A_C = Ap(nAc + 1)  ! by default point beyond the storage

        iT = 1
        iC = 1

        do orb = 1, nAc
            pos = Ap(orb)
            do while (orb >= T_orb(iT)) ; A_T(iT) = pos ; iT = iT + 1 ; end do
            do while (orb >= C_orb(iC)) ; A_C(iC) = pos ; iC = iC + 1 ; end do
        end do

    end subroutine extract_orb_data


    !> \brief   Twice cycle 4-tensor indices
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Transforms packed multi-indices of a 4-tensor, rotating them twice.
    !> This way, the same element will change its multi-index from [abcd] to [cdab].
    !>
    !> \param Rn   Number of values in Ri.
    !> \param Ri   Zero-based multiindices to transform; originally (((a-1)*nAO + (b-1))*nMO + (c-1))*nMO + (d-1).
    !> \param nAO  Maximal value of the two outer indices (a,b).
    !> \param nMO  Maximal value of the two inner indices (c,d).
    !>
    subroutine invert_indexing(Rn, Ri, nAO, nMO)

        integer, intent(inout), allocatable :: Ri(:,:)
        integer, intent(in) :: Rn, nAO, nMO

        integer :: I, J, p, q, k, l

        write (stdout, '("==================================================================")')
        write (stdout, '("Transposing [AA|MM) -> (MM|AA]",/)')

        !$omp parallel do default(none) private(I, J, k, l, p, q) shared(Rn, Ri) firstprivate(nAO, nMO)
        do I = 1, Rn

            ! decode rectangular index
            J = Ri(I,1) ; l = mod(J, nMO)
            J = J / nMO ; k = mod(J, nMO)
            J = J / nMO ; q = mod(J, nAO)
            J = J / nAO ; p = J

            ! construct zero-based multi-index of a transposed element
            Ri(I,1) = ((k * nMO + l) * nAO + p) * nAO + q

        end do

    end subroutine invert_indexing


    !> \brief   Sort the re-indexed integrals between transformation steps
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Attempts to allocate work arrays that have the same size as input (ie. BIG ones).
    !> If successful, will sort the integrals using a buffered parallel merge sort.
    !> Otherwise, serial in-place heap sort is used as a fallback option.
    !>
    !> The advantage of the merge sort is that the input array already consists of concatenated
    !> sorted sub-arrays by construction; there are nAO * nAO sorted segments, roughly of comparable size.
    !>
    !> The parallel merge sort has two stages: In the first stage, adjacent sorted sub-arrays are,
    !> in parallel, periodically collapsed (merged) together to form larger sorted segments, until
    !> the remaining number of (by now fairly large) segments becomes comparable to the number of threads.
    !> In the second stage, only one pair of segments is merged at a time (the one that will produce the smallest
    !> combined segment), but now using a parallel merge. This way, good parallel scaling is achieved
    !> in both stages, at least as far as memory bandwidth permits.
    !>
    !> \param Rn  Size of Ri, Rv.
    !> \param Ri  Integer array to sort.
    !> \param Rv  Real array to sort.
    !> \param n   Number of sorted sequences in Ri.
    !>
    subroutine sort_intermediate_integrals (Rn, Ri, Rv, n)

        use omp_lib, only: omp_get_num_threads, omp_get_wtime
        use sort,    only: heap_sort_int_float

        integer,                intent(in)    :: Rn, n
        integer,   allocatable, intent(inout) :: Ri(:,:)
        real(cfp), allocatable, intent(inout) :: Rv(:,:)

        integer, parameter :: k = 10  ! more or less arbitrary small number

        integer :: ierr0, ierr1, ierr2, nthreads, i, j, nseg, a, b, c, m
        integer,   allocatable :: Wi(:,:), offsets(:)
        real(cfp), allocatable :: Wv(:,:)
        real(wp) :: start_t, end_t

        start_t = omp_get_wtime()

        allocate (offsets(n + 1), stat = ierr0)
        allocate(Wi(size(Ri, 1), size(Ri, 2)), stat = ierr1) !F2008 version: allocate (Wi, mold = Ri, stat = ierr1)
        allocate(Wv(size(Rv, 1), size(Rv, 2)), stat = ierr2) !F2008 version: allocate (Wv, mold = Rv, stat = ierr2)

        ! fall back to the serial in-place algorithm in case of a lack of memory
        if (ierr0 /= 0 .or. ierr1 /= 0 .or. ierr2 /= 0) then
            write (stdout, '("Using serial heap sort")')
            call heap_sort_int_float(Rn, 1, Ri, Rv)
            end_t = omp_get_wtime()
            write (stdout, '("Time spent in sort_intermediate_integrals [s]: ",F25.15)') end_t - start_t
            return
        end if

        !$omp parallel
        !$omp master
        nthreads = omp_get_num_threads()
        !$omp end master
        !$omp end parallel

        write (stdout, '("Using parallel merge sort with ",I0," threads")') nthreads

        ! get segment offsets
        nseg = 0
        do i = 1, Rn
            if (i == 1) then
                nseg = nseg + 1
                offsets(nseg) = i
            else if (Ri(i - 1, 1) > Ri(i, 1)) then
                nseg = nseg + 1
                offsets(nseg) = i
            end if
        end do
        offsets(nseg + 1) = Rn + 1

        ! reduce the number of segments to a more manageable size
        !$omp parallel default(none) firstprivate(Rn, nthreads) shared(nseg, offsets, Ri, Rv, Wi, Wv) private(i, a, b, c)
        do while (nseg > k * nthreads)
            !$omp do schedule(dynamic, 1)
            do i = 1, nseg / 2
                a = offsets(2*i - 1)
                b = offsets(2*i)
                c = offsets(2*i + 1)
                call parallel_merge_sorted_int_float(1, Ri   , Ri   , Wi   , &
                                                        Rv   , Rv   , Wv   , &
                                                        a,b-1, b,c-1, a,c-1, &
                                                     1)
                Ri(a:c-1, 1) = Wi(a:c-1, 1)
                Rv(a:c-1, 1) = Wv(a:c-1, 1)
            end do
            !$omp single
            if (mod(nseg, 2) == 1) then
                Wi(offsets(nseg):Rn, 1) = Ri(offsets(nseg):Rn, 1)
                Wv(offsets(nseg):Rn, 1) = Rv(offsets(nseg):Rn, 1)
            end if
            do i = 1, (nseg + 1) / 2
                offsets(i) = offsets(2*i - 1)
            end do
            offsets((nseg + 1) / 2 + 1) = offsets(nseg + 1)
            nseg = (nseg + 1) / 2
            !$omp end single
        end do
        !$omp end parallel

        ! merge segments further until there is just one of them (= the sorted array)
        do while (nseg > 1)
            ! find the smallest double segment
            i = 1
            m = offsets(3) - offsets(1)
            do j = 2, nseg - 1
                if (m > offsets(j + 2) - offsets(j)) then
                    m = offsets(j + 2) - offsets(j)
                    i = j
                end if
            end do
            a = offsets(i)
            b = offsets(i + 1)
            c = offsets(i + 2)
            ! perform the multi-threaded merge
            call parallel_merge_sorted_int_float(nthreads, Ri   , Ri   , Wi   , &
                                                           Rv   , Rv   , Wv   , &
                                                           a,b-1, b,c-1, a,c-1, &
                                                         1)
            ! copy merged segment back to source
            Ri(a:c-1, 1) = Wi(a:c-1, 1)
            Rv(a:c-1, 1) = Wv(a:c-1, 1)
            ! collapse offsets
            offsets(i + 1 : nseg) = offsets(i + 2 : nseg + 1)
            nseg = nseg - 1
        end do

        ! check that the sort went well
        do i = 2, Rn
            if (Ri(i - 1, 1) > Ri(i, 1)) then
                write(stdout, '(I0,",",I0,": ",I0," > ",I0)') i - 1, i, Ri(i - 1, 1), Ri(i, 1)
                call xermsg ('molecular_orbital_basis_obj', 'sort_intermediate_integrals', 'Sort failed.', 1, 1)
            end if
        end do

        end_t = omp_get_wtime()
        write (stdout, '("Time spent in sort_intermediate_integrals [s]: ",F25.15)') end_t - start_t

    end subroutine sort_intermediate_integrals


    !> \brief  Merge sorted indexed arrays
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Merge two sorted integer arrays, and also corresponding real arrays. Uses the given number of threads.
    !>
    !> \param nthreads  Number of OpenMP threads to use.
    !> \param Ai  First integer array to merge.
    !> \param Bi  Second integer array to merge.
    !> \param Ci  Destination integer array (combined length of Ai, Bi).
    !> \param Af  First real array to merge.
    !> \param Bf  Second real array to merge.
    !> \param Cf  Destination real array (combined length of Af, Bf).
    !> \param As, Ae  Indices defining sections of the arrays Ai,Af to use: Ai(As:Ae), Af(As:Ae).
    !> \param Bs, Be  Indices defining sections of the arrays Bi,Bf to use: Bi(Bs:Be), Bf(Bs:Be).
    !> \param Cs, Ce  Indices defining sections of the arrays Ci,Cf to use: Ci(Cs:Ce), Cf(Cs:Ce).
    !>
    recursive subroutine parallel_merge_sorted_int_float (nthreads, Ai    , Bi    , Ci    , &
                                                                    Af    , Bf    , Cf    , &
                                                                    As, Ae, Bs, Be, Cs, Ce, &
                                                                    icol )

        integer,               intent(in)  :: nthreads, As, Ae, Bs, Be, Cs, Ce, icol
        integer,   allocatable :: Ai(:,:), Bi(:,:)  !intent(in)  
        integer,   allocatable :: Ci(:,:)           !intent(out) 
        real(cfp), allocatable :: Af(:,:), Bf(:,:)  !intent(in)  
        real(cfp), allocatable :: Cf(:,:)           !intent(out) 

        integer :: i, An, Bn, Cn, A_offsets(nthreads + 1), B_offsets(nthreads + 1), C_offsets(nthreads + 1)

        An = Ae-As+1 ! = size(Ai)
        Bn = Be-Bs+1 ! = size(Bi)
        Cn = Ce-Cs+1 ! = size(Ci)

        ! let A be larger
        if (An < Bn) then
            call parallel_merge_sorted_int_float(nthreads, Bi    , Ai    , Ci    , &
                                                           Bf    , Af    , Cf    , &
                                                           Bs, Be, As, Ae, Cs, Ce, &
                                                           icol )
            return
        end if

        ! partition A into "nthreads" segments
        do i = 1, nthreads + 1
            A_offsets(i) = (i - 1) * An / nthreads + 1
        end do

        ! bisect-find corresponding segments in B
        B_offsets(1) = Bs
        do i = 2, nthreads
            B_offsets(i) = lower_bound(Bi, icol, B_offsets(i - 1), Be, Ai(As - 1 + A_offsets(i), icol))
        end do
        B_offsets(nthreads + 1) = Be + 1

        ! here B_offsets must be relative wrt Bs
        B_offsets = B_offsets - Bs + 1

        ! determine segment sizes in C
        do i = 1, nthreads + 1
            C_offsets(i) = A_offsets(i) + B_offsets(i) - 1
        end do

        ! transform relative offsets to absolute offsets
        A_offsets = A_offsets + As - 1
        B_offsets = B_offsets + Bs - 1
        C_offsets = C_offsets + Cs - 1

        ! merge the segments independently
        !$omp parallel do num_threads(nthreads)
        do i = 1, nthreads
            call serial_merge_sorted_int_float (  &
                Ai                           , Bi                           , Ci                           ,  &
                Af                           , Bf                           , Cf                           ,  &
                A_offsets(i),A_offsets(i+1)-1, B_offsets(i),B_offsets(i+1)-1, C_offsets(i),C_offsets(i+1)-1,  &
                icol )
        end do

    end subroutine parallel_merge_sorted_int_float


    !> \brief  Merge sorted indexed arrays
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Merge two sorted integer arrays, and also corresponding real arrays.
    !>
    !> \param Ai  First integer array to merge.
    !> \param Bi  Second integer array to merge.
    !> \param Ci  Destination integer array (combined length of Ai, Bi).
    !> \param Af  First real array to merge.
    !> \param Bf  Second real array to merge.
    !> \param Cf  Destination real array (combined length of Af, Bf).
    !> \param As, Ae  Indices defining sections of the arrays Ai,Af to use: Ai(As:Ae), Af(As:Ae).
    !> \param Bs, Be  Indices defining sections of the arrays Bi,Bf to use: Bi(Bs:Be), Bf(Bs:Be).
    !> \param Cs, Ce  Indices defining sections of the arrays Ci,Cf to use: Ci(Cs:Ce), Cf(Cs:Ce).
    !>
    subroutine serial_merge_sorted_int_float (Ai    , Bi    , Ci    , &
                                              Af    , Bf    , Cf    , &
                                              As, Ae, Bs, Be, Cs, Ce, &
                                              icol)

        integer,   allocatable :: Ai(:,:), Bi(:,:)  ! intent(in)   
        integer,   allocatable :: Ci(:,:)           ! intent(inout)
        real(cfp), allocatable :: Af(:,:), Bf(:,:)  ! intent(in)   
        real(cfp), allocatable :: Cf(:,:)           ! intent(inout)
        integer,    intent(in) :: As, Ae, Bs, Be, Cs, Ce, icol

        integer :: i, j, k, An, Bn, Cn

        An = Ae-As+1 !size(Ai)
        Bn = Be-Bs+1 !size(Bi)
        Cn = Ce-Cs+1 !size(Ci)  ! == An + Bn

        i = 1
        j = 1

        do k = 1, Cn
            if (i <= An .and. (j > Bn .or. Ai(As-1 + min(i,An), icol) <= Bi(Bs-1 + min(j,Bn), icol))) then
                Ci(Cs-1 + k, icol) = Ai(As-1 + i, icol)
                Cf(Cs-1 + k, icol) = Af(As-1 + i, icol)
                i = i + 1
            else
                Ci(Cs-1 + k, icol) = Bi(Bs-1 + j, icol)
                Cf(Cs-1 + k, icol) = Bf(Bs-1 + j, icol)
                j = j + 1
            end if
        end do

    end subroutine serial_merge_sorted_int_float


    !> \brief  Index of first element that is not less
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Funny name of this function comes from the name of the standard C++ function "std::lower_bound".
    !>
    !> Returns index of the first element in a column of sorted array "A" (from "i" to "j") that is not less than
    !> the given value "v", or "j+1" if there is no such element. Takes advantage of the fact that "A"
    !> is sorted and uses interval halving.
    !>
    !> \param A    Sorted (ascending) integer array to bisect-search.
    !> \param icol Column of the array A
    !> \param i    Start of interval to search.
    !> \param j    End of interval to search.
    !> \param v    Value to compare with.
    !>
    integer function lower_bound (A, icol, i, j, v) result(right)

        integer              :: icol, i, j, v
        integer, allocatable :: A(:,:) !intent(in)

        integer :: left, mid

        left  = i
        right = j + 1

        do while (left /= right)
            mid = (left + right) / 2
            if (A(mid,icol) < v) then
                left = mid + 1
            else
                right = mid
            end if
        end do

    end function lower_bound


    !> \brief   Fix indexing of the integral storage
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Convert the greedy (rectangular) working indexing scheme used in two_electron_integrals_sparse to the compact (triangular) one.
    !>
    !> \param Rn  Size of the intex array.
    !> \param Ri  Index array.
    !> \param nMO Maximal value of individual sub-indices (rectangular size).
    !>
    subroutine rect_index_to_tri_index(Rn, Ri, nMO, icol)

        integer, allocatable :: Ri(:,:) !intent(inout)
        integer, intent(in)  :: Rn, nMO, icol

        integer :: i, n, a, b, c, d, u, v

        !$omp parallel do default(none) private(i, n, a, b, c, d, u, v) shared(Rn, Ri, icol) firstprivate(nMO)
        do i = 1, Rn

            ! decode rectangular index
            n = Ri(i, icol)   ; d = 1 + mod(n, nMO)
            n = n / nMO ; c = 1 + mod(n, nMO)
            n = n / nMO ; b = 1 + mod(n, nMO)
            n = n / nMO ; a = 1 + n

            ! construct triangular index
            u     = max(a,b) * (max(a,b) - 1) / 2 + min(a,b)
            v     = max(c,d) * (max(c,d) - 1) / 2 + min(c,d)
            Ri(i, icol) = max(u,v) * (max(u,v) - 1) / 2 + min(u,v)

        end do

    end subroutine rect_index_to_tri_index


    !> \brief   Final stage of calculation of the 2e integrals
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Writes integrals and indices to disk and releases allocated memory.
    !>
    subroutine finalize_two_electron_integrals(this, Rn, Ri, Rv, integral_storage, integral_options, column_descriptor)

        use const, only: two_p_sym_ints, ijkl_indices_header

        class(molecular_orbital_basis_obj)         :: this
        class(integral_options_obj), intent(in)    :: integral_options
        class(integral_storage_obj), intent(inout) :: integral_storage

        character(len=line_len), allocatable, intent(in) :: column_descriptor(:)

        real(kind=cfp), allocatable, intent(inout), target :: Rv(:,:)
        integer,        allocatable, intent(inout), target :: Ri(:,:)
        integer,                     intent(in)            :: Rn

        type(p2d_array_obj), pointer :: mo_integrals
        type(p2d_array_obj), target  :: tmp_mo_integrals
        real(kind=cfp),      pointer :: backup_a(:,:)

        character(len=line_len) :: mo_header, ind_header

        integer :: I, J, u, v, a, b, c, d, current_pos, first_record, last_record, lunit, ierr, imaster

        mo_integrals => tmp_mo_integrals

        ! initialize integral storage
        if (mo_integrals % init(Rn, 1, 0, column_descriptor) /= 0) then
                call xermsg ('molecular_basis_mod', 'finalize_two_electron_integrals', &
                             'Molecular integral storage initialization failed.', 1, 1)
        end if

!         open (11, form = 'formatted')
!         do I = 1, Rn
!             J = Ri(I,1)
!             u = int(0.5 + sqrt(2.0 * J)) ; v = J - u * (u - 1) / 2
!             a = int(0.5 + sqrt(2.0 * u)) ; b = u - a * (a - 1) / 2
!             c = int(0.5 + sqrt(2.0 * v)) ; d = v - c * (c - 1) / 2
!             write (11, '(I15,1x,I5,1x,I5,1x,I5,1x,I5,1x,E25.15)') J, a, b, c, d, Rv(i,1)
!         end do
!         close (11)

        ! push the sorted integrals and indices to the storage
        backup_a => mo_integrals % a; mo_integrals % a => Rv
        this % ijkl_indices => Ri
        this % ind_ijkl_integral = Rn

        write (stdout, '("Saving integrals to disk...")')

        imaster    = master
        lunit      = integral_storage % integral_file % get_unit_no()
        mo_header  = integral_storage % contruct_header_string(this % get_basis_name(), two_p_sym_ints)
        ind_header = integral_storage % contruct_header_string(this % get_basis_name(), ijkl_indices_header)

        integral_storage % data_header % name = mo_header

        first_record = integral_storage % integral_file % start_record(mo_header)
        call integral_options % write(lunit, first_record, current_pos)
        call mo_integrals % write(lunit, current_pos, last_record, imaster)
        call integral_storage % integral_file % close_record(mo_header, first_record, last_record)

        first_record = integral_storage % integral_file % start_record(ind_header)
        call this % write_ijkl_indices(lunit, first_record, last_record)
        call integral_storage % integral_file % close_record(ind_header, first_record, last_record)
        this % ijkl_indices => null()

        write (stdout, '("...done")')

        ! release molecular integrals
        mo_integrals % a => backup_a
        if (mo_integrals % final() /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals', &
                         'Deallocation of the temporary integral array failed.', 5, 1)
        end if

    end subroutine finalize_two_electron_integrals


   subroutine two_electron_integrals(this,integral_storage,integral_options)
      use const, only: abel_prod_tab, Mib, cache_line_size, two_el_ints, two_p_sym_ints, ijkl_indices_header
      use mpi_mod 
      use omp_lib
      use parallel_arrays
      use special_functions, only: ipair, unpack_pq
      use sort, only: sort_int_float, heap_sort_int_float
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(integral_options_obj), intent(in) :: integral_options
      class(integral_storage_obj), intent(inout) :: integral_storage

      type(integral_storage_obj) :: ao_integrals_disk
      type(integral_options_obj) :: ao_int_opt
      integer :: p, q, pq, r, s, rs, i, j, ij, lunit, no_int, err, first_record, last_record, no_ao, int_type
      integer :: no_pairs, no_tot, my_start, my_end, block, n_threads, thread_id
      integer :: sym_i, sym_j, orb_i, orb_j, orb_j_it, d1, d2, current_pos
      integer :: rs_start, rs_end, rs_block, no_rs_blocks, rs_block_size, last_tgt_ao, n_ints, no_blocks_ao
      integer, allocatable :: rs_ind(:,:), ij_orbital_range(:,:), last_tgt_mo(:), ind_ijkl_integral(:), ijkl_indices(:,:), &
                              ij_offset(:), n_integrals(:)
      integer(kind=1), allocatable :: ij_type(:)
      type(p2d_array_obj), target :: integral_src, integral_tgt !we really need two of these in case disk-to-disk AO-MO run is required
      type(p2d_array_obj), pointer :: ao_integrals, mo_integrals
      real(kind=cfp), allocatable :: cf(:,:), cf_t(:,:), iqrs(:,:), ijrs(:,:), ijks_t(:,:), hlp(:), ijkl_integrals(:,:,:), &
                                     cf_t_non_zero(:,:)
      integer, allocatable :: mo_indices(:,:), n_cf_t_non_zero(:)
      real(kind=wp) :: trmo_start, trmo_end, start_t, end_t, total_trmo, total_ijrs, total_ijkl, val
      logical :: i_is_continuum, j_is_continuum
      logical :: ao_is_local = .true.
      character(len=line_len), allocatable :: column_descriptor(:)
      character(len=line_len) :: ao_header, mo_header, ind_header

      integer, parameter :: no_blocks = 0, padding = 2*cache_line_size/cfp_bytes
      logical, parameter :: mo_is_local = .true. !same meaning as no_blocks = 0

         start_t = omp_get_wtime()
 
         call mpi_mod_barrier(err)
 
         write(stdout,'("--------->","molecular_orbital_basis_obj:two_electron_integrals")')
 
         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', 'The basis set has not been initialized.',1, 1)
         end if

         if (this % number_of_functions == 0) then
            call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                         'Number of molecular orbitals on input is zero.', 2, 1)
         end if

         if (.not.associated(this%ao_integral_storage) .or. .not.associated(this%ao_basis)) then
            call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                         'On input at least one of this%ao_integral_storage, this%ao_basis have not been associated.', 4, 1)
         endif

         !Header for the AO integrals that we're looking for
         ao_header = this%ao_integral_storage%contruct_header_string(this%ao_basis%get_basis_name(),two_el_ints)

         !Header for the MO integrals that will be calculated
         mo_header = integral_storage%contruct_header_string(this%get_basis_name(),two_p_sym_ints)

         !Header for the MO integrals that will be calculated
         ind_header = integral_storage%contruct_header_string(this%get_basis_name(),ijkl_indices_header)

         ! In this section we associate ao_integrals which contains the input AO integrals with the appropriate source.
         ! In case the AO integrals are stored in memory then we point directly to the array holding them.
         ! If the AO integrals are on disk then we load them into the local array 'integral' and set the pointer ao_integrals to that.
         ! At the moment 1p integral transform using SHARED input AO integrals is not supported.
         if (this%ao_integral_storage%in_memory()) then
            ao_integrals => this%ao_integral_storage%integral
            if (this%ao_integral_storage%data_header%name .ne. ao_header) then
               call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                            'The AO integrals on input are not compatible with the AO basis set for the MO basis set.', 5, 1)
            endif
            write(stdout,'("AO integrals with header: ",a)') ao_header
            ao_is_local = .not.(ao_integrals%have_offsets()) !find out if the AO integrals in the memory are shared or local
         endif

         !load the AO integrals into memory as LOCAL arrays (if they are stored on disk)
         if (this%ao_integral_storage%on_disk()) then
            write(stdout,'("Loading AO integrals from the disk...")')

            err = ao_integrals_disk%init(memory=integral_src)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', 'Memory allocation 3 failed.', err, 1)
            end if

            call ao_integrals_disk%read(this%ao_integral_storage,ao_header,ao_int_opt,ao_is_local)
            ao_integrals => ao_integrals_disk%integral !this points to the local array 'integral_src'
            write(stdout,'("AO integrals with header: ",a)') ao_header
            write(stdout,'("...done")')
         endif

         !BEYOND THIS POINT ao_integrals POINTS TO AN ARRAY CONTAINING THE AO INTEGRALS TO BE TRANSFORMED

         call ao_integrals%get_array_dim(d1,d2,no_blocks_ao) !This gives the dimensions of ao_integrals%a(1:d1,1:d2)
         call ao_integrals%get_column_descriptor(column_descriptor)

         write(stdout,'("On input there is ",i0," types of AO integrals")') d2
         write(stdout,'("Total number of AO integrals of each type: ",i0)') d1
         if (.not.(ao_is_local)) then
            write(stdout,'("AO integrals are scattered among all processes.")')
         else
            write(stdout,'("Every process keeps the full set of AO integrals.")')
         endif

         ! Note that instead of loading the AO integrals we can calculate them now since we have the pointer to the AO basis set
         ! (this%ao_basis) and the AO integral routine (this%ao_integrals)...

         !Assign a pointer to the requested array given by the user.
         if (integral_storage%in_memory()) then
            integral_storage%data_header%name = mo_header

            mo_integrals => integral_storage%integral
         endif

         !If we request the output to be stored on disk then start a new record on the data file that will contain the integrals.
         !Assign a pointer to a temporary array.
         if (integral_storage%on_disk()) then
            !temporary storage for the integrals
            mo_integrals => integral_tgt

            lunit = integral_storage%integral_file%get_unit_no()             !unit that is associated to the file opened
            first_record = integral_storage%integral_file%start_record(mo_header) !position, within the data file, of the first record available for the integral data
         endif

         !BEYOND THIS POINT mo_integrals POINTS TO AN ARRAY CONTAINING THE TRANSFORMED MO INTEGRALS

         ! 2-PARTICLE INTEGRAL TRANSFORM STARTS HERE: the AO integrals are accessed through the pointer ao_integrals;
         ! the MO integrals are accessed through the pointer mo_integrals

         no_ao = this%ao_basis%number_of_functions !total number of AOs

         !todo The ijrs array can be made much smaller once symmetry for the AOs has been implemented
         rs = no_ao*(no_ao+1)/2
         no_pairs = rs

         ij = this%number_of_functions !total number of MOs
         ij = ij*(ij+1)/2

         allocate(cf_t(this%number_of_functions,no_ao),rs_ind(1:2,rs),ij_type(ij),last_tgt_mo(this%no_irr),stat=err)
         if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 4 failed.', err, 1)

         !Copy the orbital coefficients to one array: this relies on the fact that the molecular orbitals are indexed symmetry by symmetry.
         call this%get_orbital_coefficient_matrix(cf)

         !Transposing the MO coefficient matrix and the other matrices speeds up the computation significantly.
         cf_t = transpose(cf)

         !Get index of the last AO representing the target electrons: this procedure assumes all target functions preceed the continuum ones.
         !todo once symmetry for the AOs has been implemented I should get this for each symmetry, i.e. last_tgt_ao will be an array
         last_tgt_ao = 0
         do i=1,no_ao
            j = this%ao_basis%indices_to_shells(1,i) !shell index containing the i-th function
            if (this%ao_basis%shell_descriptor(3,j) .eq. 1) cycle !skip continuum functions
            last_tgt_ao = max(last_tgt_ao,i)
         enddo

         !Prepare the set of (rs)-indices that will be split among all threads
         rs = 0
         do r=1,no_ao
            do s=1,r
               rs = rs + 1
               rs_ind(1:2,rs) = (/r,s/)
            enddo !s
         enddo !r

         !Determine the type of the orbital pairs
         do i=1,this%no_irr
            !Find the last target orbital within each symmetry: this procedure assumes all target orbitals preceed the continuum ones.
            last_tgt_mo(i) = 0
            do j=1,this%orbital_data(i)%number_of_functions
               if (this%is_continuum(this%relative_to_absolute(j,i))) cycle
               last_tgt_mo(i) = max(last_tgt_mo(i),this%relative_to_absolute(j,i))
            enddo !j
         enddo

         ij = 0
         do i=1,this%number_of_functions
            p = this%absolute_to_relative(2,i)
            i_is_continuum = .false.
            if (i > last_tgt_mo(p) .and. i <= this % relative_to_absolute(this % orbital_data(p) % number_of_functions, p)) then
                i_is_continuum = .true.
            end if
            do j=1,i
               q = this%absolute_to_relative(2,j)
               j_is_continuum = .false.
               if (j > last_tgt_mo(q) .and. j <= this % relative_to_absolute(this % orbital_data(q) % number_of_functions, q)) then
                  j_is_continuum = .true.
               end if
               ij = ij + 1
               ij_type(ij) = 1 !TT
               if (i_is_continuum .neqv. j_is_continuum) then
                  ij_type(ij) = 2 !TC or CT
               elseif (i_is_continuum .and. j_is_continuum) then
                  ij_type(ij) = 3 !CC
               endif
            enddo !j
         enddo !i

         !Loop over all types of AO integrals to transform.
         do int_type=1,d2

            trmo_start = omp_get_wtime()

            write(stdout,'("Transforming AO integral type ",i0," ...")') int_type

            !We have to do this since ijrs and ij_orbital_range are being allocated below (i.e. inside the int_type loop).
            if (allocated(ijrs)) deallocate(ijrs)
            if (allocated(ij_orbital_range)) deallocate(ij_orbital_range)
            no_int = 0

            ! Note: This parallel section needs to use DEFAULT(SHARED) to allow work with polymorphic objects with gfortran.
            !$OMP PARALLEL &
            !$OMP & DEFAULT(SHARED) &
            !$OMP & PRIVATE(rs,r,s,p,q,ij,i,j,pq,block,n_threads,thread_id,err,rs_block_size,rs_start,rs_end,rs_block, &
            !$OMP &         no_rs_blocks,ij_orbital_range,sym_i,sym_j,orb_i,orb_j,orb_j_it,my_start,my_end,iqrs,ijks_t,hlp, &
            !$OMP &         ij_offset,val,n_ints,no_tot) &
            !$OMP & SHARED(no_pairs,ao_is_local,no_ao,this,int_type,cf,cf_t,ijrs,ao_integrals,rs_ind,integral_options, &
            !$OMP &        mo_integrals,total_trmo,trmo_start,trmo_end,total_ijrs,total_ijkl,stdout,last_tgt_ao,ij_type, &
            !$OMP &        ijkl_integrals,ind_ijkl_integral,ijkl_indices,n_integrals,nprocs,cf_t_non_zero, &
            !$OMP &        mo_indices,n_cf_t_non_zero,myrank,column_descriptor,d2) &
            !$OMP & REDUCTION(+:no_int)
            n_threads = omp_get_num_threads()
            thread_id = omp_get_thread_num()

            p = this%no_irr*(this%no_irr+1)/2  !total number of symmetry pairs

            !The ij_orbital_range array is used in the second part of the integral transform to split the (i,j) indices in each symmetry combination among the threads.
            !In the buffers allocated below (iqrs,ijrs,ijks_t) the letters q,r,s stand for atomic orbitals and the letters i,j,k stand for the molecular orbitals.
            !It proved necessary to perform the allocations of the arrays iqrs,ijks_t,hlp explicitly here rather than allocating them before the parallel region and making
            !them PRIVATE. For large calculations (large AO basis) these arrays were causing SEGFAULTs even with unlimited stack set on the compute nodes (the stack can never be
            !unlimited in practice so there is always a maximum that can be breached).
            ij = this%number_of_functions*(this%number_of_functions+1)/2
            allocate(ij_orbital_range(4,p),iqrs(this%number_of_functions,no_ao),&
                        ijks_t(no_ao,this%number_of_functions),hlp(no_ao),ij_offset(ij),stat=err)
            if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 5 failed.', err, 1)
            ij_orbital_range = 0

            !split the unique (i,j) orbital indices in each combination of symmetries among the threads
            do sym_i=1,this%no_irr
               do sym_j=1,sym_i

                  p = sym_i*(sym_i-1)/2 + sym_j !the index of the symmetry block for the pair of symmetries corresponding to the pair (sym_i,sym_j), sym_i .g.e sym_j

                  !calculate the total number of orbital pairs for this combination of symmetries
                  if (sym_i .eq. sym_j) then
                     ij = this%orbital_data(sym_i)%number_of_functions*(this%orbital_data(sym_i)%number_of_functions+1)/2
                  else
                     ij = this%orbital_data(sym_i)%number_of_functions*this%orbital_data(sym_j)%number_of_functions
                  endif

                  ij_orbital_range(1,p) = 1; ij_orbital_range(3,p) = -1
                  if (ij .eq. 0) then 
                     cycle
                  endif

                  if (ij .le. n_threads) then
                     my_start = thread_id+1
                     my_end = thread_id+1
                  else
                     block = ceiling(ij/(n_threads*1.0))
                     my_start = block*thread_id+1 !ij-index for my first orbital pair
                     my_end = min(ij,block*(thread_id+1)) !ij-index for my last orbital pair
                  endif

                  !since my_end is always at most equal to ij we must make sure that if this thread cannot encounter the 'if (ij .eq. my_end)' below 
                  !in case it is not supposed to process any orbitals for this pair of symmetries.
                  if (my_start > my_end) cycle

                  ij = 0
                  do orb_i=1,this%orbital_data(sym_i)%number_of_functions   !over all MOs in symmetry sym_i
                     if (sym_i .eq. sym_j) then
                        orb_j_it = orb_i          !both orbitals are from the same symmetry and therefore the loop over the second orbital must be only from orb_j_it = 1 to orb_j_it = orb_i.
                     else
                        orb_j_it = this%orbital_data(sym_j)%number_of_functions !both orbitals come from different symmetries so the loop over the second loop must be over all orbitals in that symmetry.
                     endif

                     do orb_j=1,orb_j_it
                        ij = ij + 1

                        !Determine the first and the last pair of orbital indices that I'll be transforming, for each pair of symmetries.
                        if (ij .eq. my_start) then
                           ij_orbital_range(1,p) = orb_i
                           ij_orbital_range(2,p) = orb_j
                        endif
                        if (ij .eq. my_end) then
                           ij_orbital_range(3,p) = orb_i
                           ij_orbital_range(4,p) = orb_j
                        endif
                        if (ij > my_end) then
                           exit
                        endif
                     enddo !orb_j
                  enddo !orb_i
               enddo !sym_j
            enddo !sym_i

            ij = this%number_of_functions*(this%number_of_functions+1)/2 !number of (i,j) combinations of orbitals

            if (integral_options%max_ijrs_size > 0.0_cfp) then
               !Limit the size of the ijrs array to the maximum size (in Mib) supplied by the user.
               i = integral_options%max_ijrs_size*Mib/(cfp_bytes*1.0_cfp) !total number of elements allowed in the ijrs array
               rs_block_size = min(i/(ij+padding),no_pairs)
               !$OMP SINGLE
               write(stdout,'("Size of the ijrs array limited to (Mib): ",f8.3)') integral_options%max_ijrs_size
               !$OMP END SINGLE
            else
               no_tot = this%block_offset(size(this%block_offset))  !the sum of integrals in each of the unique symmetry blocks
               !Set-up the size of the (rs)-block (i.e. the number of (r,s)-indices in one block) relatively to the size of the array of the transformed integrals.
               rs_block_size = min(no_tot/ij,no_pairs)
            endif
            
            no_rs_blocks = no_pairs/rs_block_size
            if (no_rs_blocks*rs_block_size < no_pairs) no_rs_blocks = no_rs_blocks + 1

            call generate_ij_offset(this,ij_type,ij_orbital_range,integral_options%two_p_continuum,ij_offset,n_ints)

            !$OMP SINGLE
            allocate(n_integrals(n_threads),stat=err)
            if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 6a failed.', err, 1)
            n_integrals = 0
            !$OMP END SINGLE

            !Total number of ijkl integrals generated by this thread:
            n_integrals(thread_id+1) = n_ints
            !$OMP BARRIER

            !$OMP SINGLE
            write(stdout,'(/,"Total number of (rs)-indices: ",i0)') no_pairs
            write(stdout,'("Number of (rs)-index blocks: ",i0,", number of (rs)-indices in block: ",i0)') &
                no_rs_blocks, rs_block_size
            write(stdout,'(/,"Size of the ijrs buffer to be allocated [Mib]: ",f25.15)') &
                (rs_block_size*(ij+padding)*cfp_bytes*1.0_cfp)/Mib

            !Allocate the output arrays if we request the output to be stored in memory or on the disk.
            !Note that we enforce the mo_integrals array to be LOCAL on each process (no_blocks=0), i.e. every process gets
            !the full integrals array. The present transformation algorithm doesn't 
            !implement splitting of the transformed (ij|kl) integrals among the processes. We allocate space for a non-indexed
            !(that is purely local) array with d2 columns and no_tot rows.
            !The columns in the mo_integrals array correspond to the types of AO integrals we have on input.
            if (int_type .eq. 1) then
               no_tot = sum(n_integrals)
               write(stdout,'("Total number of MO integrals of each type that will be obtained: ",i0)') no_tot
               err = mo_integrals%init(no_tot,d2,no_blocks,column_descriptor)
               if (err /= 0) then
                  call xermsg ('molecular_basis_mod', 'two_electron_integrals', &
                               'Array initialization 2 failed; see p2d_array_obj%init.', err, 1)
               end if
               !Allocate the array of indices for the transformed integrals.
               if (associated(this%ijkl_indices)) deallocate(this%ijkl_indices)
               allocate(this%ijkl_indices(no_tot,d2),stat=err)
               if (err /= 0) then
                  call xermsg ('molecular_basis_mod', 'two_electron_integrals', &
                               'Array initialization 3 failed; see p2d_array_obj%init.', err, 1)
               end if
               this%ijkl_indices = 0
               this%ind_ijkl_integral = 0 !index of the last integral stored in mo_integrals%a
            endif

            i = maxval(n_integrals)
            !todo The ijrs array can be made much smaller once symmetry for the AOs has been implemented
            allocate(ijrs(ij+padding,rs_block_size),ijkl_integrals(i,n_threads,int_type:int_type), &
                     ind_ijkl_integral(n_threads),ijkl_indices(i,n_threads),stat=err)
            if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 6 failed.', err, 1)
            ijrs(:,:) = 0.0_cfp
            total_ijrs = 0.0_cfp
            total_ijkl = 0.0_cfp
            ijkl_integrals = 0.0_cfp
            ind_ijkl_integral = 0
            ijkl_indices = 0
            val = 2*i*n_threads*cfp_bytes/(Mib*1.0_cfp)
            write(stdout,'("Memory storage for temporary arrays ijkl_integrals,ijkl_indices (MiB): ",f8.3)') val
            write(stdout,'("Memory has been successfuly allocated.",/)')

            call extract_non_zero_cf(cf_t,cf_t_non_zero,mo_indices,n_cf_t_non_zero)
            !$OMP END SINGLE

            !Total number of ijkl integrals generated by this thread:
            ind_ijkl_integral(thread_id+1) = n_integrals(thread_id+1)

            !We divide the whole set of (r,s)-indices into blocks and for each block we transform the first two indices (pq|rs)->(ij|rs) and then calculate their contribution 
            !to the final set of transformed integrals: (ij|rs) -> (ij|kl).
            do rs_block=1,no_rs_blocks

               !Determine the range of rs-indices in the current block
               rs_start = rs_block_size*(rs_block-1)+1
               rs_end = min(rs_block_size*rs_block,no_pairs)

               !$OMP SINGLE
               trmo_start = omp_get_wtime()
               !$OMP END SINGLE

               !I. THE FIRST PART OF THE INTEGRAL TRANSFORM:
               !Transform the first two indices: (pq|rs) -> (ij|rs) for all (i,j)-indices (pairs of orbitals) and the current
               !block of (rs)-indices (rs_start:rs_end). Here we loop over the (rs) indices of the AO integrals (pq|rs) and
               !obtain the partially transformed integrals (ij|rs), where i,j are the molecular orbitals. It is only this part
               !of the integral transformation algorithm where need to accommodate the fact that the AO integrals may be incomplete
               !if they were split among the processes (ao_is_local).
               if (ao_is_local) then !every process keeps all AO integrals
                  call omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local ( &
                            ao_integrals,int_type,rs_start,rs_end,ij_type,iqrs,hlp,this,cf,cf_t,ijrs,no_ao,last_tgt_ao, &
                            no_int,no_pairs,rs_ind,integral_options%two_p_continuum &
                  )
               else !in case the AO integrals are scattered randomly among the processes
                  call omp_two_p_transform_pqrs_block_to_ijrs_AO_is_not_local ( &
                            ao_integrals,int_type,rs_start,rs_end,ij_type,iqrs,hlp,this,cf,cf_t,ijrs,no_ao,last_tgt_ao, &
                            no_int,no_pairs,rs_ind,integral_options%two_p_continuum &
                  )
               endif
               !Remember that we don't need BARRIER here only as long as the loop over rs-indices inside the preceeding routines
               !is split using OMP DO since that implies synchronization.

               !Putting the barrier here ensures accurate timing is reported below: we want the timing for the slowest thread to be reported.
               !$OMP BARRIER

               !$OMP SINGLE
               trmo_end = omp_get_wtime()
               write(stdout,'("(rs)-index block: ",i0)') rs_block
               write(stdout,'("(pq|rs) -> (ij|rs) for rs-block: ",i0,"->",i0," done in ",f8.3," [s]")') &
                    rs_start,rs_end, trmo_end-trmo_start
               total_ijrs = total_ijrs + trmo_end-trmo_start
               trmo_start = omp_get_wtime()
               !$OMP END SINGLE

               !II. THE SECOND PART OF THE INTEGRAL TRANSFORM:
               !Transform the last two indices: (ij|rs) -> (ij|kl)
               !Note that we use as the ijks array the array iqrs since they have the same size.
               call omp_two_p_transform_ijrs_block_to_ijkl ( &
                        ijrs,iqrs,ijks_t,cf,cf_t_non_zero,mo_indices,n_cf_t_non_zero,no_ao,this,rs_ind,ij_type,ao_is_local, &
                        integral_options%tol,ijkl_integrals,ijkl_indices,int_type,rs_start,rs_end,no_pairs,ij_orbital_range, &
                        integral_options%two_p_continuum,thread_id,ij_offset &
               )

               !Putting the barrier here ensures accurate timing is reported below: we want the timing for the slowest thread to be reported.
               !$OMP BARRIER

               !$OMP SINGLE
               trmo_end = omp_get_wtime()
               write(stdout,'("(ij|rs) -> (ij|kl) contribution done in ",f8.3," [s]")') trmo_end-trmo_start
               total_ijkl = total_ijkl + trmo_end-trmo_start
               !$OMP END SINGLE

            enddo !rs_block

            !$OMP SINGLE
            deallocate(ijrs)
            !$OMP END SINGLE

            deallocate(ijks_t,iqrs,hlp)

            !$OMP SINGLE
            !todo put the integrals into blocks using ij_offset but generated
            !for the global array not for individual threads.
            this%ind_ijkl_integral = 0
            do i=1,n_threads
               do j=1,ind_ijkl_integral(i)
                  !At this stage we cannot skip integrals when more processes are used since the integral array must be reduced first.
                  if (nprocs .eq. 1) then
                     if (ijkl_integrals(j,i,int_type) .eq. 0.0_cfp) cycle
                  endif
                  !WARNING: if saving only the non-zero integrals any indexing
                  !making use of ij_offset becomes wrong!
                  this%ind_ijkl_integral = this%ind_ijkl_integral + 1
                  if (this%ind_ijkl_integral > size(mo_integrals%a,1)) then
                     print *,this%ind_ijkl_integral,size(mo_integrals%a,1)
                     stop "error 2:buffer too small"
                  endif
                  mo_integrals%a(this%ind_ijkl_integral,int_type) = ijkl_integrals(j,i,int_type)
                  !WARNING: note that the present implementation may not work correctly if int_type > 1 since the this%ind_ijkl_integral should depend on int_type too!
                  this%ijkl_indices(this%ind_ijkl_integral,int_type) = ijkl_indices(j,i)
!                  print *,'ind',this%ind_ijkl_integral,this%ijkl_indices(this%ind_ijkl_integral,int_type)
               enddo !j
            enddo !i

            deallocate(ind_ijkl_integral,ijkl_integrals,ijkl_indices,ij_offset)
            !$OMP END SINGLE

            !$OMP END PARALLEL

            write(stdout,'(/,"Sorting the final array of integrals...")')
            !save the unsorted data so we can recover it later if needed
            !open(newunit=i,file='unsorted_integrals',status='replace',form='unformatted')
            !write(i) this%ind_ijkl_integral
            !write(i) this%ijkl_indices(1:this%ind_ijkl_integral,int_type)
            !write(i) this%ind_ijkl_integral
            !write(i) mo_integrals%a(1:this%ind_ijkl_integral,int_type)
            !close(i)
            trmo_start = omp_get_wtime()
            !Heap sort is generally faster than quicksort for the arrays that
            !we're dealing with. Actually, this may be caused by the current quicksort implementation not being very good.
            call heap_sort_int_float(this%ind_ijkl_integral,int_type,this%ijkl_indices,mo_integrals%a)
            !call sort_int_float(this%ind_ijkl_integral,int_type,this%ijkl_indices,mo_integrals%a) !Quick sort
            trmo_end = omp_get_wtime()
            write(stdout,'("...done and took [s]: ",f8.3)') trmo_end-trmo_start

!            do i=1,this%ind_ijkl_integral
!               !todo temp debug
!               write(50,'(2i15,e25.15)') i,this%ijkl_indices(i,int_type),mo_integrals%a(i,int_type)
!            enddo

            total_trmo = total_ijrs + total_ijkl
            write(stdout,'(/,"Number of AO integrals used: ",i0)') no_int
            write(stdout,'("Total time for (pq|rs) -> (ij|rs) ",f8.3," [s]")') total_ijrs
            write(stdout,'("Total time for (ij|rs) -> (ij|kl) ",f8.3," [s]")') total_ijkl
            write(stdout,'("(pq|rs) -> (ij|kl) done in ",f8.3," [s]")') total_trmo
            write(stdout,'("Number of ijkl integrals: ",i15,1X,f8.3)') &
                this % ind_ijkl_integral, this % ind_ijkl_integral * cfp_bytes / (Mib * 1.0_cfp)

         enddo !int_type

         deallocate(cf_t)

         !Get rid of the AO integrals now: in case of very large MO integrals array the MPI reduction below may need a lot of memory.
         err = ao_integrals%final()
         if (err /= 0) then
            call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Finalizing the ao_integrals array failed.', err, 1)
         end if

         !If the AO integrals were stored only in memory then we have to inform the user that we've destroyed them now.
         if (this % ao_integral_storage%in_memory()) then
            write(stdout,'("WARNING: The AO integrals array supplied on input has been destroyed.")')
         end if

         !THIS IS AN IMPORTANT STEP: in case the AO integral array was shared (and MO integral array is local) then
         !mo_integrals%a contains only the CONTRIBUTION of each process to the final transformed 
         !integrals. Therefore we need to REUDCE (sum) the contributions to the MO integrals from each process.
         !The results are communicated to all processes so at the end each process keeps its own copy
         !of the full MO integrals array (i.e. this is compatible with the LOCAL attribute of the mo_integrals array).
         if (.not.(ao_is_local) .and. mo_is_local) then
            if (cfp .eq. ep1) then
               call mpi_xermsg ('molecular_basis_mod', 'two_electron_integrals', &
                                'In quad precision alltoall reduction not implemented. &
                                &Only master will have the final MO integrals.', 1, 0)
               call mo_integrals%reduce_a_to_master
            else
               call mo_integrals%reduce_a_local
            endif

            !Get rid of the transformed integrals smaller than the threshold. Note that in the serial run this has been done above.
            !todo this should only be done on the master process if cfp .eq. ep1
            write(stdout,'("Deleting integrals smaller than the threshold...")')
            do int_type=1,d2
               !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i) SHARED(mo_integrals,int_type,integral_options)
               !$OMP DO
               do i=1,size(mo_integrals%a,1)
                  if (abs(mo_integrals%a(i,int_type)) < integral_options%tol) mo_integrals%a(i,int_type) = 0.0_cfp
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
            enddo
            write(stdout,'("...done")')
         endif

         !If requested print the non-zero integrals
         if (integral_options%print_integrals) then
            if (cfp == ep1) then
               !Only the master process has the final integrals
               if (myrank == master) then
                  call mo_integrals%print(.true.)
               else
                  write(stdout,'("The molecular integrals are printed only by the master process, see log_file.0")')
               endif
            else
               call mo_integrals%print(.true.)
            endif
         endif

         !Dump all integrals to disk and close the record.
         if (integral_storage%on_disk()) then

            write(stdout,'("Saving integrals to disk...")')

            !The first record are the integral options.
            call integral_options%write(lunit,first_record,current_pos)

            !The second record are the ordered integrals: only master writes them to disk.
            i = master
            call mo_integrals%write(lunit,current_pos,last_record,i)

            !Every process closes the record so that they all keep identical header information.
            call integral_storage%integral_file%close_record(mo_header,first_record,last_record)

            err = mo_integrals%final()

            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                             'Deallocation of the temporary integral array failed.', 5, 1)
            end if

            !Write out the ijkl indices
            first_record = integral_storage%integral_file%start_record(ind_header) !position, within the data file, of the first record available for the integral data

            call this%write_ijkl_indices(lunit,first_record,last_record)

            call integral_storage%integral_file%close_record(ind_header,first_record,last_record)

            write(stdout,'("...done")')

         endif
 
         write(stdout,'("<---------","molecular_orbital_basis_obj:two_electron_integrals")')
 
         call mpi_mod_barrier(err)

         end_t = omp_get_wtime()

         write(stdout,'("Two_electron_integrals took [s]: ",f25.15)') end_t-start_t
 
   end subroutine two_electron_integrals

   function integral_index(this,integral_type,bf_indices)
      use const
      use special_functions, only: ipair
      implicit none
      class(molecular_orbital_basis_obj) :: this
      character(len=*), intent(in) :: integral_type
      integer, intent(in) :: bf_indices(:,:)
      integer :: integral_index(size(bf_indices,2))

      integer :: ind,i,j,k,iAB, iCD, a,b,c,d
      logical :: found

         if (size(bf_indices,1) .eq. 2) then !1-electron integral index

            !todo this can be improved by swapping the do-loop with select case for
            !each integral type. Best to pack the indexing computation into an elementary routine.
            do k=1,size(bf_indices,2)
               i = maxval(bf_indices(1:2,k))
               j = minval(bf_indices(1:2,k))

               ind = i*(i-1)/2+j

               !At the moment all 1-electron atomic integrals are index in the
               !same way.
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
                  call xermsg ('molecular_orbital_basis_obj', 'integral_index', &
                               'Unrecognized one electron molecular integral type on input.', 1, 1)

               end select
            enddo !k

         elseif (size(bf_indices,1) .eq. 4) then !2-electron integral index
            do i=1,size(bf_indices,2)
               a = maxval(bf_indices(1:2,i))

               b = minval(bf_indices(1:2,i))
               c = maxval(bf_indices(3:4,i))
               d = minval(bf_indices(3:4,i))
!write(*,'(4i3)') a,b,c,d
               if (a .ge. c) then
                  iAB = ipair(a)+b
                  iCD = ipair(c)+d
               else
                  iAB = ipair(c)+d
                  iCD = ipair(a)+b
               endif
               ind = ipair(iAB)+iCD
               ind = ipair(max(iAB,iCD))+min(iAB,iCD)
!write(*,'(4i3)') ipair(a),ipair(c),iAB,iCD
!write(*,'(4i3)') ind,ipair(iAB),iCD
               !Indexing relying on the sorted array of molecular integrals
               call bisect_index(ind,this%ijkl_indices,1,this%ind_ijkl_integral,j,found)
               
               if (found) then
                  integral_index(i) = j
               else
                  integral_index(i) = -1 !this integral is zero
               endif
            enddo !i
         else
            call xermsg ('molecular_orbital_basis_obj', 'integral_index', &
                         'The number of orbital indices on input must be either 2 or 4.', 3, 1)
         endif

   end function integral_index

   function get_basis_name(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      character(len=line_len) :: get_basis_name

         get_basis_name = "molecular_orbital_basis_obj"
 
   end function get_basis_name

   function get_shell_name(this,i)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      character(len=line_len) :: get_shell_name
      integer, intent(in) :: i

      type(orbital_data_obj) :: dummy

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_name', 'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this % number_of_shells) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_name', 'On input the value of i was out of range.', 2, 1)
         end if

         get_shell_name = dummy%name()

   end function get_shell_name

   subroutine get_shell_data(this,i,shell_data)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(shell_data_obj), intent(out) :: shell_data
      integer, intent(in) :: i

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_data', 'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this % number_of_shells) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_data', 'On input the value of i was out of range.', 2, 1)
         end if
         
         select type (shell => shell_data)
            type is (orbital_data_obj)
               shell = this%orbital_data(i)
            class default
               call xermsg ('molecular_orbital_basis_obj', 'get_shell_data', &
                            'On input shell_data must be of orbital_data_obj type.', 3, 1)
         end select

   end subroutine get_shell_data

   subroutine get_all_orbital_sets(this,orbital_sets,number_of_orbital_sets)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(out) :: number_of_orbital_sets
      type(orbital_data_obj), allocatable :: orbital_sets(:)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_all_orbital_sets', 'The basis set has not been initialized.', 1, 1)
         end if

         number_of_orbital_sets = 0

         !todo this must involve transfering all orbital data from the type-bound arrays into the orbital_sets structure.
         call xermsg ('molecular_orbital_basis_obj', 'get_all_orbital_sets', 'Not implemented yet.',1,1)
 
   end subroutine get_all_orbital_sets

   subroutine get_orbital_coefficient_matrix(this,cf)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      real(kind=cfp), allocatable :: cf(:,:)

      integer :: err, cnt, i

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_coefficient_matrix', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (allocated(cf)) deallocate(cf)
         allocate(cf(this%ao_basis%number_of_functions,this%number_of_functions),stat=err)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_coefficient_matrix', 'Memory allocation failed.', err, 1)
         end if

         !Place the orbital coefficients into columns of cf symmetry by symmetry.
         cnt = 0
         do i=1,this%no_irr
            cf(1:this%ao_basis%number_of_functions,cnt+1:cnt+this%orbital_data(i)%number_of_functions) = &
            &this%orbital_data(i)%coefficients(1:this%ao_basis%number_of_functions,1:this%orbital_data(i)%number_of_functions)
            cnt = cnt + this%orbital_data(i)%number_of_functions
         enddo !i

   end subroutine get_orbital_coefficient_matrix

   function is_initialized(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      logical :: is_initialized

         is_initialized = this%initialized

   end function is_initialized

   !todo change this work just like the atomic one!!!
   subroutine get_continuum_flags(this,irr,list)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: irr
      logical, allocatable :: list(:)

      integer :: i, err

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_continuum_flag', 'The basis set has not been initialized.', 1, 1)
         end if

         if (irr <= 0 .or. irr > this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'get_continuum_flag', 'On input the value of irr was out of range.', 2, 1)
         end if

         if (allocated(list)) deallocate(list)
         allocate(list(this%orbital_data(irr)%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj', 'get_continuum_flag', 'Memory allocation failed.',err,1)

         list = .false.
         do i=1,this%orbital_data(irr)%number_of_functions
            list(i) = this%is_continuum(this%relative_to_absolute(i,irr))
         enddo !i

   end subroutine get_continuum_flags

   function get_number_of_orbitals(this,irr)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: irr
      integer :: get_number_of_orbitals

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_number_of_orbitals', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (irr <= 0 .or. irr > this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'get_number_of_orbitals', &
                         'On input the value of irr was out of range.', 2, 1)
         end if

         get_number_of_orbitals = this%orbital_data(irr)%number_of_functions

   end function get_number_of_orbitals

   function get_index_within_symmetry(this,absolute_index)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: absolute_index
      integer :: get_index_within_symmetry

         if (.not. this % initialized) then
            call xermsg('molecular_orbital_basis_obj', 'get_index_within_symmetry', &
                        'The basis set has not been initialized.', 1, 1)
         end if

         if (absolute_index <= 0 .or. absolute_index > this % number_of_functions) then
            call xermsg ('molecular_orbital_basis_obj', 'get_index_within_symmetry', &
                         'On input absolute_index was out of range.', 2, 1)
         end if

         get_index_within_symmetry = this%absolute_to_relative(1,absolute_index)

   end function get_index_within_symmetry

   function get_orbital_symmetry(this,absolute_index)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: absolute_index
      integer :: get_orbital_symmetry

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_symmetry', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (absolute_index <= 0 .or. absolute_index > this % number_of_functions) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_symmetry', &
                         'On input absolute_index was out of range.', 2, 1)
         end if

         get_orbital_symmetry = this%absolute_to_relative(2,absolute_index)

   end function get_orbital_symmetry

   function get_absolute_index(this,num,sym)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: num,sym
      integer :: get_absolute_index

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_absolute_index', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (sym > this % no_irr .or. sym <= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'get_absolute_index', &
                         'On input orbital symmetry was out of range.', 2, 1)
         end if
         if (num > this % orbital_data(sym) % number_of_functions .or. num <= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'get_absolute_index', &
                         'On input orbital index within the given symmetry was out of range.', 3, 1)
         end if

         get_absolute_index = this%relative_to_absolute(num,sym) !the index of the MO within the whole orbital set
  
   end function get_absolute_index 

   subroutine calculate_amplitudes(this,a,normalize_to_a,amplitudes,continuum_channels)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      real(kind=cfp), intent(in) :: a
      logical, intent(in) :: normalize_to_a
      integer, allocatable :: continuum_channels(:,:)
      real(kind=cfp), allocatable :: amplitudes(:,:)

      real(kind=cfp), allocatable :: ao_amplitudes(:,:)
      integer :: n_channels, n_orbitals, err, cnt, i, j, k

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'amplitudes', 'The basis set has not been initialized.', 1, 1)
         end if

         write(stdout,'("--------->","molecular_orbital_basis_obj:calculate_amplitudes")')

         if (a .le. 0.0_cfp) call xermsg('molecular_orbital_basis_obj', 'amplitudes', 'On input a .le. 0.0_cfp: invalid input.',2,1)

         !Calculate amplitudes of the atomic functions:
         call this%ao_basis%calculate_amplitudes(a,normalize_to_a,ao_amplitudes,continuum_channels)
         n_channels = size(ao_amplitudes,1)

         !Contract with the molecular orbital coefficients:
         if (allocated(amplitudes)) deallocate(amplitudes)
         allocate(amplitudes(n_channels,this%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj', 'amplitudes', 'Memory allocation failed.',err, 1)

         cnt = 0
         do i=1,this%no_irr

            if (this % orbital_data(i) % number_of_coefficients /= this % ao_basis % number_of_functions) then
                call xermsg ('molecular_orbital_basis_obj', 'amplitudes', 'The AO and MO basis sets are incompatible.', 3, 1)
            end if

            n_orbitals = this%orbital_data(i)%number_of_functions

            amplitudes(1:n_channels,cnt+1:cnt+n_orbitals) = matmul(ao_amplitudes,this%orbital_data(i)%coefficients)

            do j=1,n_orbitals
               do k=1,n_channels
                  if (amplitudes(k,cnt+j) .ne. 0.0_cfp) write(stdout,'(3(i0,1x),e25.15)') i,j,k, amplitudes(k,cnt+j)
               enddo !k
            enddo !j

            cnt = cnt + n_orbitals

         enddo !i

         write(stdout,'("<---------","molecular_orbital_basis_obj:calculate_amplitudes")')

   end subroutine calculate_amplitudes

   subroutine orthogonalize(this,overlap_matrix,symmetry,gramm_schmidt,gramm_schmidt_one_by_one,symmetric,sym_ortho_data,&
                            active_start,active_end,passive_start,passive_end,check_overlaps)
      use orthogonalization, only: GS_ortho_routine, SYM_ortho_routine
      use const, only: thrs_gs_ortho, thrs_symm_ortho, overlap_ints, int_del_thr, thrs_lin_dep_gs_ortho, thrs_cf_sym_ortho_trans
      use parallel_arrays, only: p2d_array_obj
      implicit none
      class(molecular_orbital_basis_obj) :: this
      real(kind=cfp), allocatable :: overlap_matrix(:,:)
      type(sym_ortho_io), intent(inout), optional :: sym_ortho_data
      integer, intent(in) :: symmetry
      logical, intent(in), optional :: gramm_schmidt, gramm_schmidt_one_by_one, symmetric, check_overlaps
      integer, intent(in), optional :: active_start, active_end, passive_start, passive_end

      integer :: i, j, k, no_cf, err, a_start, a_end, p_start, p_end, p_start_abs, p_end_abs
      real(kind=cfp), allocatable :: olap_orbs(:,:), thrs_ao(:)
      logical :: do_gs_ortho, do_sym_ortho, all_active, check
      real(kind=cfp) :: thrs, bto_thrs, cgto_thrs

         write(stdout,'("--------->","molecular_orbital_basis_obj:orthogonalize")')

         write(stdout,'(/,"Orbital orthogonalization")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         do_gs_ortho = .false.
         do_sym_ortho = .false.

         if (present(gramm_schmidt) .and. present(symmetric)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Specify only one of: gramm_schmidt, symmetric but not both.', 2, 1)
         end if
         if (present(gramm_schmidt_one_by_one) .and. (present(gramm_schmidt) .or. present(symmetric))) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'If only gramm_schmidt_one_by_one is given then none of gramm_schmidt and symmetric can be given.', 2, 1)
         end if

         if (present(gramm_schmidt)) do_gs_ortho = gramm_schmidt
         if (present(gramm_schmidt_one_by_one)) do_gs_ortho = gramm_schmidt_one_by_one
         if (present(symmetric)) do_sym_ortho = symmetric

         if (do_gs_ortho .and. present(sym_ortho_data)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'On input sym_ortho_data given but Gramm-Schmidt ortho. requested.', 3, 1)
         end if
         if (do_sym_ortho .and. .not. present(sym_ortho_data)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Symmetric ortho. requested but input sym_ortho_data not given.', 4, 1)
         end if

         if (present(active_start) .neqv. present(active_end)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Only one of active_start,active_end has been specified &
                         &but both are required if one of them is given.', 5, 1)
         endif

         if (present(passive_start) .neqv. present(passive_end)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Only one of passive_start,passive_end has been specified &
                         &but both are required if one of them is given.', 6, 1)
         endif

         if (symmetry <= 0 .or. symmetry > this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'On input, the value of symmetry was out of range.', 7, 1)
         end if
!
!--------Determine and check absolute indices of the active and passive orbitals. Transfer the values of active_* and passive_* to the variables a_* and p_*.
!
         if (present(active_start)) then

            if (active_end < active_start .or. active_start <= 0 .or. &
                active_start > this % orbital_data(symmetry) % number_of_functions) then
               print *,active_start,active_end,this%orbital_data(symmetry)%number_of_functions
               call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                            'On input, active_start, active_end was incorrect', 8, 1)
            endif

            write(stdout,'("Symmetry: ",i1)') symmetry
            write(stdout,'("Range of active orbitals within the symmetry: ",2i5)') active_start,active_end

            a_start = active_start
            a_end = active_end

         else
            a_start = 1
            a_end = this%orbital_data(symmetry)%number_of_functions
         endif

         if (present(passive_start)) then

            if (passive_end < passive_start .or. passive_start <= 0 .or. &
                passive_start > this % orbital_data(symmetry) % number_of_functions) then
               print *,passive_start,passive_end
               call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                            'On input, passive_start, passive_end was incorrect', 9, 1)
            endif

            write(stdout,'("Symmetry: ",i1)') symmetry
            write(stdout,'("Range of passive orbitals within the symmetry: ",2i5)') passive_start,passive_end

            p_start = passive_start
            p_end = passive_end

         else
            p_start = 0
            p_end = 0
         endif

         if (a_start .eq. 1 .and. a_end .eq. this%orbital_data(symmetry)%number_of_functions) then
            all_active = .true.
         else
            all_active = .false.
         endif

         if ((p_end .ge. a_start .and. p_start .le. a_start) .or. (p_end .ge. a_end .and. p_start .le. a_end)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'The ranges of active and passive orbitals must be disjunct.', 10, 1)
         endif

         !by default we request checking of the overlaps but this can be overrided by the user
         check = .true.
         if (present(check_overlaps)) check = check_overlaps

         no_cf = this%ao_basis%number_of_functions       !the number of AOs in this symmetry
!
!------- Determine the deletion thresholds for the orbital coefficients.
!
         allocate(thrs_ao(no_cf),stat=err)
         if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 3 failed.', err, 1)

         if (do_gs_ortho) then
            bto_thrs = F1MACH(4,cfp_dummy)
            cgto_thrs = thrs_lin_dep_gs_ortho
         elseif (do_sym_ortho) then
            bto_thrs = F1MACH(4,cfp_dummy) !thrs_cf_sym_ortho_trans
            cgto_thrs = thrs_cf_sym_ortho_trans
         endif

         !If the basis contains BTOs then typically there will be many small coefficients 
         !which must not be neglected otherwise the resulting orbitals will not be orthogonal with a sufficient precision.
         if (this%ao_basis%contains_btos()) then
            cgto_thrs = bto_thrs
            write(stdout,'(10x,"BTO coefficients with magnitude smaller than: ",e25.15," will be deleted.")') bto_thrs
            write(stdout,'(10x,"CGTO coefficients with magnitude smaller than: ",e25.15," will be deleted.")') cgto_thrs
            j = 0
            thrs_ao = cgto_thrs
            do i=1,this%ao_basis%number_of_shells
               k = this%ao_basis%shell_descriptor(5,i) !number of functions in the shell

               !this%ao_basis%shell_descriptor(1,i) == 2 for BTO shells
               if (this%ao_basis%shell_descriptor(1,i) .eq. 2) thrs_ao(j+1:j+k) = bto_thrs

               j = j + k
            enddo !i
         else
            thrs_ao = cgto_thrs
            write(stdout,'(10x,"Coefficients with magnitude smaller than: ",e25.15," will be deleted.")') cgto_thrs
         endif
!
!--------Gramm-Schmidt orthogonalization
!
         if (do_gs_ortho) then

            write(stdout,'(/,10x,"Gramm-Schmidt orthogonalization requested")')

            !Column indices in the mo2so_range are absolute
            if (p_end .ge. p_start .and. p_end > 0) then
               p_start_abs = this%relative_to_absolute(p_start,symmetry)
               p_end_abs = this%relative_to_absolute(p_end,symmetry)
            else
               p_start_abs = 1
               p_end_abs = 0
               p_start = 1 !setting p_e smaller than p_s guarantees that no orbitals will be treated as fixed in GS_ortho_routine
               p_end = 0
            endif

            if (present(gramm_schmidt_one_by_one)) then
               write(stdout,'(10x,"The active orbitals will be orthogonalized one-by-one way.")')
               !Orthogonalize a set of orbitals in the one-by-one fashion (this is most likely to be used for the continuum orthogonalization)
               !$OMP PARALLEL DEFAULT(NONE) &
               !$OMP & PRIVATE(i) &
               !$OMP & SHARED(a_start,a_end,p_start_abs,p_end_abs,no_cf,p_start,p_end,this,overlap_matrix,symmetry,thrs_ao)
               !$OMP DO
               do i=a_start,a_end
                  call GS_ortho_routine (no_cf, p_start, p_end, i, i, this % orbital_data(symmetry) % coefficients, &
                                         overlap_matrix, symmetry, this % mo2so_range(1:2, p_start_abs:max(p_end_abs, 1)), thrs_ao)
               enddo !i
               !$OMP END DO
               !$OMP END PARALLEL
            else
               call GS_ortho_routine (no_cf, p_start, p_end, a_start, a_end, this % orbital_data(symmetry) % coefficients, &
                                      overlap_matrix, symmetry, this % mo2so_range(1:2, p_start_abs:max(p_end_abs,1)), thrs_ao)
            endif

            write(stdout,'("...finished.")')

         endif
!
!--------Symmetric orthogonalization
!
         if (do_sym_ortho) then

            err = sym_ortho_data%check(this%no_irr)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                             'sym_ortho_data%check failed. See the routine for details.', err, 1)
            end if

            if (allocated(sym_ortho_data%to_delete)) deallocate(sym_ortho_data%to_delete)
            allocate(sym_ortho_data%to_delete(1:this%orbital_data(symmetry)%number_of_functions),stat=err)
            if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 1 failed.', err, 1)

            !This array will contain .true. values for those absolute orbital indices which should be deleted from the orbital basis.
            sym_ortho_data%to_delete(:) = .false.

            !Symmetrically orthogonalize a subset of orbitals in the set. Typically, this would be used following G-S orthogonalization of the continuum against the target orbitals 
            !to orthogonalize the subset of continuum orbitals among themselves.
            write(stdout,'(/,10x,"Symmetric orthogonalization requested")')

            write(stdout,'("Deletion threshold: ",e25.15)') sym_ortho_data%del_thrs(symmetry)

            allocate(olap_orbs(a_start:a_end,a_start:a_end),stat=err)
            if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 2 failed.', err, 1)

            !Form the overlap matrix in the basis of the orbitals to orthogonalize
            !olaps = c**T * S * c
            olap_orbs(a_start:a_end,a_start:a_end) = matmul( &
                matmul( &
                    transpose(this % orbital_data(symmetry) % coefficients(1:no_cf,a_start:a_end)), &
                    overlap_matrix(1:no_cf,1:no_cf) &
                ), &
                this % orbital_data(symmetry) % coefficients(1:no_cf, a_start:a_end) &
            )

            call SYM_ortho_routine (no_cf, a_start, a_end, this % orbital_data(symmetry) % coefficients, &
                                    olap_orbs, sym_ortho_data % del_thrs(symmetry), thrs_ao, sym_ortho_data % to_delete)

            deallocate(olap_orbs)

            k = count(sym_ortho_data%to_delete)
            if (k > 0) then
               write(stdout,'(/,"Number of orbitals marked for deletion: ",i0)') k

               write(stdout,'("Orbital indices: number.symmetry, absolute index: ")')
               do i=a_start,a_end
                  if (sym_ortho_data % to_delete(i)) then
                     write(stdout,'(2i5,i0)') this % absolute_to_relative(1,i), this % absolute_to_relative(2,i), i
                  end if
               enddo
            else
               write(stdout,'(/,"There are no orbitals to delete in this symmetry.")')
            endif
         endif
!
!--------The orbitals have changed so recalculate the so2mo,mo2so indices.
!
         call this%determine_auxiliary_indices
!
!--------Optional checking of orbital orthogonality
!
         if (check) then

            if (do_gs_ortho) thrs = thrs_gs_ortho
            if (do_sym_ortho) thrs = thrs_symm_ortho

            if (do_sym_ortho) then !if we performed the symmetric orthogonalization then we must not check orthogonality between the orbitals that are to be deleted
               do i=a_start,a_end
                  if (sym_ortho_data%to_delete(i)) a_end = a_end - 1 !the orbitals to delete are always the last ones in each symmetry so we just lower the index of the last active orbital
               enddo
            endif

            no_cf = this%ao_basis%number_of_functions !the number of AOs in this symmetry

            if (a_end .ge. a_start .and. a_end > 0) then

               write(stdout,'(/,10x,"Checking orthogonality of the active orbitals in symmetry: ",i0,/)') symmetry
   
               allocate(olap_orbs(a_start:a_end,a_start:a_end),stat=err)
               if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 3 failed.', err, 1)
   
               !form the overlap matrix in the basis of the orthogonalized orbitals
               !olaps = c**T * S * c
               associate(orbital => this%orbital_data(symmetry))
                  olap_orbs(a_start:a_end,a_start:a_end) = matmul(matmul(transpose(orbital%coefficients(1:no_cf,a_start:a_end)),&
                                                                              &overlap_matrix(1:no_cf,1:no_cf)),&
                                                                              &orbital%coefficients(1:no_cf,a_start:a_end))
               end associate
      
               do i=a_start,a_end
                  write(stdout,'(i5,1x,e25.15)') i,olap_orbs(i,i)
                  if (abs(olap_orbs(i,i)-1.0_cfp) .ge. thrs) then
                     write(stdout,'(2e25.15)') thrs,olap_orbs(i,i)
                     call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                                  'Norm of an orbital is not 1 within accuracy given by thrs.', 14, 0)
                  endif
                  !If we orthogonalize in the one-by-one way then we cannot expect that the active orbitals will be orthogonal with each other so this check is skipped.
                  if (.not.(present(gramm_schmidt_one_by_one))) then
                     do j=i+1,a_end
                        if (olap_orbs(j,i) .ge. thrs) then
                           write(stdout,'(i0,1x,i0,2e25.15)') j,i,olap_orbs(j,i),thrs
                           call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                                        'Overlap between two different orbitals is .ge. thrs.', 15, 0)
                        endif
                     enddo
                  endif
               enddo
   
               deallocate(olap_orbs)

               if (p_end .ge. p_start .and. p_end > 0) then
                  write(stdout,'(/,10x,"Checking orthogonality of the active orbitals &
                                        &with respect to the fixed orbitals in symmetry: ",i0,/)') symmetry
   
                  allocate(olap_orbs(p_start:p_end,a_start:a_end),stat=err)
                  if (err /= 0) then
                     call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 4 failed.', err, 1)
                  end if
      
                  !form the overlap matrix in the basis of the orthogonalized orbitals
                  !olaps = c**T * S * c
                  associate(orbital => this%orbital_data(symmetry))
                     olap_orbs(p_start:p_end,a_start:a_end) = matmul(matmul(transpose(orbital%coefficients(1:no_cf,p_start:p_end)),&
                                                                                   &overlap_matrix(1:no_cf,1:no_cf)),&
                                                                                   &orbital%coefficients(1:no_cf,a_start:a_end))
                  end associate

                  do i=a_start,a_end
                     do j=p_start,p_end
                        if (olap_orbs(j,i) .ge. thrs) then
                           write(stdout,'(i0,1x,i0,2e25.15)') j,i,olap_orbs(j,i),thrs
                           call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                                        'Overlap between an active and a passive orbital is .ge. thrs.', 17, 0)
                        endif
                     enddo
                  enddo
      
                  deallocate(olap_orbs)
               endif

            endif

            write(stdout,'(/,"Orthogonalization complete.")')

         else

            write(stdout,'(/,"Checking of orbital overlaps not requested. Orthogonalization complete.")')

         endif

         write(stdout,'("<---------","done:molecular_orbital_basis_obj:orthogonalize")')
 
   end subroutine orthogonalize

   subroutine delete_small_coefficients(this)
      use const, only: thrs_orb_cf
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: symmetry, j, k

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'delete_small_coefficients', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         write(stdout,'("--------->","molecular_orbital_basis_obj:delete_small_coefficients")')

         write(stdout,'("Removing orbital coefficients with magnitude smaller than: ",e25.15)') thrs_orb_cf

         do symmetry=1,size(this%orbital_data)
            do k=1,size(this%orbital_data(symmetry)%coefficients,2)
               do j=1,this%orbital_data(symmetry)%number_of_coefficients
                  if (abs(this % orbital_data(symmetry) % coefficients(j,k)) < thrs_orb_cf) then
                     this % orbital_data(symmetry) % coefficients(j,k) = 0.0_cfp
                  end if
               enddo !j
            enddo !k
         enddo !i

         write(stdout,'("<---------","done:molecular_orbital_basis_obj:delete_small_coefficients")')

   end subroutine delete_small_coefficients

   subroutine delete_orbitals(this,symmetry,to_delete)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: symmetry
      logical, intent(in) :: to_delete(:)

      type(orbital_data_obj) :: temp_orbitals
      integer :: n_orbitals_to_keep, err, i, cnt

         write(stdout,'("--------->","molecular_orbital_basis_obj:delete_orbitals")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'delete_orbitals', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         if (size(to_delete) /= this % orbital_data(symmetry) % number_of_functions) then
            call xermsg ('molecular_orbital_basis_obj', 'delete_orbitals', &
                         'The size of the array on input is not equal to the number of orbitals in the given symmetry.', 2, 1)
         end if

         temp_orbitals = this%orbital_data(symmetry)

         n_orbitals_to_keep = this%orbital_data(symmetry)%number_of_functions - count(to_delete)

         deallocate(this % orbital_data(symmetry) % energy, &
                    this % orbital_data(symmetry) % occup, &
                    this % orbital_data(symmetry) % spin, &
                    this % orbital_data(symmetry) % coefficients)

         allocate(this % orbital_data(symmetry) % energy(n_orbitals_to_keep), &
                  this % orbital_data(symmetry) % occup(n_orbitals_to_keep), &
                  this % orbital_data(symmetry) % spin(n_orbitals_to_keep), &
                  this % orbital_data(symmetry) % coefficients(this % orbital_data(symmetry) % number_of_coefficients, &
                    n_orbitals_to_keep), stat = err)
         if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'delete_orbitals', 'Memory allocation failed.', err, 1)
 
         this%orbital_data(symmetry)%number_of_functions = n_orbitals_to_keep

         !Keep only the orbitals that are not marked for deletion
         cnt = 0
         do i=1,temp_orbitals%number_of_functions
            if (.not.(to_delete(i))) then
               cnt = cnt + 1
               this%orbital_data(symmetry)%energy(cnt) = temp_orbitals%energy(i)
               this%orbital_data(symmetry)%occup(cnt) = temp_orbitals%occup(i)
               this%orbital_data(symmetry)%spin(cnt) = temp_orbitals%spin(i)
               this%orbital_data(symmetry)%coefficients(1:temp_orbitals%number_of_coefficients,cnt) &
                    = temp_orbitals%coefficients(1:temp_orbitals%number_of_coefficients,i)
            else
               write(stdout,'("Orbital ",i0," from symmetry ",i1," has been deleted.")') i, symmetry
            endif
         enddo !i

         this%number_of_functions = sum(this%orbital_data(:)%number_of_functions)

         call this%determine_auxiliary_indices

         write(stdout,'("<---------","molecular_orbital_basis_obj:delete_orbitals")')

   end subroutine delete_orbitals

   function check_sym_ortho_io(this,last_irr)
      implicit none
      class(sym_ortho_io) :: this
      integer :: check_sym_ortho_io
      integer, intent(in) :: last_irr

      integer :: i

         check_sym_ortho_io = 0

         do i=1,last_irr
            if (this%del_thrs(i) .le. 0.0_cfp) check_sym_ortho_io = i
         enddo
      
   end function check_sym_ortho_io

   subroutine determine_auxiliary_indices(this)
      use const, only: abel_prod_tab
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: no_tot, no_ao, i, j, k, l, col, err, ij_irr, p, q, kl_irr, kl, ij, no_sym_pairs, n_p, no_pairs
      integer, allocatable :: sym_seq_no(:)
      logical, allocatable :: ao_is_continuum(:)
      logical :: all_zero

         no_tot = this%number_of_functions
         no_ao = this%ao_basis%number_of_functions

         no_pairs = this%number_of_functions*(this%number_of_functions+1)/2
         no_sym_pairs = this%no_irr*(this%no_irr+1)/2   !number of unique pairs of symmetries
         no_sym_pairs = no_sym_pairs*(no_sym_pairs+1)/2 !number of unique pairs of pairs of symmetries 

         n_p = maxval(this%orbital_data(:)%number_of_functions)
         n_p = n_p*(n_p+1)/2

         if (allocated(this%block_offset)) deallocate(this%block_offset)
         if (allocated(this%sym_offset)) deallocate(this%sym_offset)
         if (allocated(this%is_continuum)) deallocate(this%is_continuum)
         if (allocated(this%absolute_to_relative)) deallocate(this%absolute_to_relative)
         if (allocated(this%relative_to_absolute)) deallocate(this%relative_to_absolute)
         if (allocated(this%so2mo_range)) deallocate(this%so2mo_range)
         if (allocated(this%mo2so_range)) deallocate(this%mo2so_range)

         allocate(this % block_offset(no_sym_pairs + 1), &
                  this % sym_offset(n_p), &
                  this % is_continuum(no_tot), &
                  this % absolute_to_relative(2,no_tot), &
                  this % relative_to_absolute(no_ao, this % no_irr), &
                  stat = err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj','add_shell','Memory allocation 2 failed.',err,1)

         this%is_continuum = .false.
         this%block_offset = 0
         this%absolute_to_relative = 0

         allocate(this%so2mo_range(1:2,no_ao),this%mo2so_range(1:2,no_tot),stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj','add_shell','Memory allocation failed.',err,1)

         !determine the so2mo_range, mo2so_range scanning the whole orbital coefficient data and looking for orbitals which have non-zero AO coefficients.
         this%so2mo_range(1,:) = no_tot+1
         this%so2mo_range(2,:) = 0
         this%mo2so_range(1,:) = no_ao+1
         this%mo2so_range(2,:) = 0
         col = 0
         !loop over all symmetries:
         do i=1,this%number_of_shells
            associate(orbital => this%orbital_data(i))
            !Get continuum flags for each AO for this symmetry.
            call this%ao_basis%get_continuum_flags(this%orbital_data(i)%irr,ao_is_continuum)
            !all orbitals within a symmetry:
            do j=1,orbital%number_of_functions
               col = sum(this%orbital_data(1:i-1)%number_of_functions) + j !the MOs are stored sequentially in columns symmetry by symmetry.
               !write(stdout,'("Orbital with num.sym ",i,".",i," has index",i)') orbital%sym_ind, orbital%symmetry, col
               this%absolute_to_relative(1:2,col) = (/j,orbital%irr/) !index within its own symmetry (j) and symmetry of the orbital (orbital%irr)
               this%relative_to_absolute(j,orbital%irr) = col !the index of the MO within the whole orbital set
               all_zero = .true.
               do k=1,no_ao
                  if (orbital%coefficients(k,j) .ne. 0.0_cfp) then
                     all_zero = .false.
                     if (col < this%so2mo_range(1,k)) this%so2mo_range(1,k) = col
                     if (col > this%so2mo_range(2,k)) this%so2mo_range(2,k) = col
                     if (k < this%mo2so_range(1,col)) this%mo2so_range(1,col) = k
                     if (k > this%mo2so_range(2,col)) this%mo2so_range(2,col) = k
                     if (ao_is_continuum(k)) this%is_continuum(col) = .true. !this orbital must represent the continuum since it contains non-zero coefficients for at least one continuum AO.
                  endif
               enddo !k
               ! If the basis includes trivial orbitals then the identification of target and continuum orbitals becomes ambiguous.
               ! This is not necessarily a problem but there is no physical reason why such orbitals should be included in the basis.
               if (all_zero) then
                  print *,i,j
                  print *,orbital%coefficients(:,j)
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The orbital basis includes at least one orbital with only zero-containing coefficient array. &
                               &This is not allowed.', 7, 1)
               endif
            enddo !j
            end associate
         enddo !i

         !Determine the number of transformed MO 2-particle symmetric integrals for each unique quadruplet of symmetries (ij|kl)
         allocate(sym_seq_no,source=this%block_offset,stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj','add_shell','Memory allocation 2 failed.',err,1)

         this%block_offset = 0
         sym_seq_no = 0
         p = 0
         do i=1,this%no_irr
            do j=1,i
               !How many unique pairs of MOs there are for the symmetry pair (i,j)
               if (i .eq. j) then
                  ij = this%orbital_data(i)%number_of_functions*(this%orbital_data(i)%number_of_functions+1)/2
               else
                  ij = this%orbital_data(i)%number_of_functions*this%orbital_data(j)%number_of_functions
               endif 
               ij_irr = abel_prod_tab(i,j) !IRR of the product of symmetries (i,j)
               p = i*(i-1)/2+j             !index of the unique pair of symmetries (i,j)
               q = 0
               do k=1,i
                  do l=1,k
                     q = k*(k-1)/2+l
                     kl_irr = abel_prod_tab(k,l)
                     if (p .ge. q .and. ij_irr .eq. kl_irr) then !only unique quartets of symmetries; additonally we have the symmetry restriction IRR_ij .eq. IRR_kl
                        if (k .eq. l) then
                           kl = this%orbital_data(k)%number_of_functions*(this%orbital_data(k)%number_of_functions+1)/2
                        else
                           kl = this%orbital_data(k)%number_of_functions*this%orbital_data(l)%number_of_functions
                        endif
                        !below we compute the number of integrals for the quartet of symmetries (p,q), p .ge. q.
                        if (p .eq. q) then  !(i,j) pair = (k,l) pair
                           sym_seq_no(p*(p-1)/2+q) = ij*(ij+1)/2 !only unique quartets of MOs
                        else
                           sym_seq_no(p*(p-1)/2+q) = ij*kl
                        endif
                        !write(stdout,'("block",4i,2i)') max(i,j),min(i,j),max(k,l),min(k,l), p*(p-1)/2+q, sym_seq_no(p*(p-1)/2+q)
                     endif
                  enddo !l
               enddo !k
            enddo !j
         enddo !i

         !Clearly, from the construction below, block_offset(p) contains the total number of integrals that preceed the given symmetry block 'p=s*(s-1)/2+t'. The values s,t are defined as:
         !s=i*(i-1)/2+j, t=k*(k-1)/2+l where i,j,k,l are the IRRs of the MOs in (ab|cd) and we assume i .ge. j, k .ge. l.
         !The 'integral block' is defined here as the set of integrals (ab|cd) where the MOs a,b,c,d come from the symmetry block 'p'.
         !Within each integral block the index of the integral is given depending on the symmetries of the MOs: (II|II), (IJ|IJ), (II|JJ), (IJ|KL).
         !The full index therefore is: block_offset(p) + block_index; the expressions for the block_index for the four cases above can be found in the relevent 2-particle indexing routines.

         j = size(this%block_offset)-1 !total number of unique quartets of symmetries; note that the last value in block_offset is the total number of 2-particle symmetric integrals.
         this%block_offset(j+1) = sum(sym_seq_no)

         do i=2,j
            this%block_offset(i) = sum(sym_seq_no(1:i-1)) !compute the number of integrals preceeding the block i.
         enddo !i
         this%block_offset(1) = 0

         deallocate(sym_seq_no)

    end subroutine determine_auxiliary_indices

!---- ROUTINES FOR 2P TRANSFORMATION:
    !> This routine transposes the 2D matrix on input while determining whether all elements of the array are zero.
    subroutine transpose_2d(batch_in,batch_t,nrow,ncol)
       use const, only: tile
       implicit none
       integer, intent(in) :: nrow,ncol
       real(kind=cfp), allocatable :: batch_in(:,:) !(nrow,ncol)
       real(kind=cfp), allocatable :: batch_t(:,:) !(ncol,nrow)
       logical :: all_zero
 
       integer :: i, j, ii, jj, iend, jend
 
          all_zero = .true.
 
          do jj = 1,nrow,tile
             jend = min (nrow,jj+tile-1)
             do ii = 1,ncol,tile
                iend = min (ncol,ii+tile-1)
 
                if (.not.(all_zero)) then
                   do j = jj,jend
                      do i = ii,iend
                         batch_t(i,j) = batch_in(j,i)
                      enddo
                   enddo
                else
                   do j = jj,jend
                      do i = ii,iend
                         batch_t(i,j) = batch_in(j,i)
                         if (batch_t(i,j) .ne. 0.0_cfp) all_zero = .false.
                      enddo
                   enddo
                endif
 
             enddo
          enddo
       
    end subroutine transpose_2d
 
    !> \warning Note that the use of the 'allocatable' attribute for the argument arrays is key for performance since this
    !> attribute allows the compiler to assume that the arrays are contiguous in memory.
    !> Alternatively the 'contiguous' attribute can be used but it is a F2008 feature so we omit it here.
    !> It is assumed that the threads have been launched outside of this routine.
    subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local(ao_integrals,int_type,rs_start,rs_end,ij_type,iqrs,&
                                        hlp,mobas,cf,cf_t,ijrs,no_ao,last_tgt_ao,no_int,no_pairs,rs_ind,two_p_continuum)
       use omp_lib
       implicit none
       class(molecular_orbital_basis_obj), intent(in) :: mobas
       type(p2d_array_obj), intent(inout) :: ao_integrals
       integer, intent(in) :: int_type,rs_start,rs_end,no_ao,last_tgt_ao,no_pairs,rs_ind(2,no_pairs)
       integer(kind=1), allocatable :: ij_type(:)
       integer, intent(inout) :: no_int
       logical, intent(in) :: two_p_continuum
 
       !we use the allocatable attribute so that the compiler can assume that the arrays are unit-stride:
       !this has ~20% impact on the performance of this routine.
       real(kind=cfp), allocatable :: iqrs(:,:), hlp(:), ijrs(:,:), cf(:,:), cf_t(:,:)
 
       integer :: rs, pq, p, q, ij, i, j, max_p, max_q, ijrs_type, pq_mapped, rs_mapped, rs_tmp, iam, r, s
       real(kind=cfp) :: mo_int
       logical :: ao_contains_continuum, dont_two_p_continuum

          if (mobas%ao_basis%n_cont_fns > 0) then
             ao_contains_continuum = .true.
          else
             ao_contains_continuum = .false.
          endif

          dont_two_p_continuum = ao_contains_continuum .and. .not.(two_p_continuum)

          !Split the rs-indices in the block among threads
          !$OMP DO SCHEDULE(DYNAMIC)
          do rs=rs_start,rs_end !rs is the index of the (rs)-pair corresponding to the AO integrals (pq|rs) that we want to process here.
             iam = omp_get_thread_num()
       
             ! Pre-load all (available) AO integrals (pq|O|rs) (== (pq|rs)) with the given (rs) index
             ! into hlp array - this speeds up the 1st step a little bit compared with accessing ao_integrals%a
             ! in the j-loop below. We preload the integrals always only one row at a time so that we don't need
             ! an extra array of size ~no_ao**2.
    
             !1st step - transform the first index: p->i
             ! We loop over all AO integrals (pq|rs) with the given rs index: if integrals for only 1p in the continuum
             ! are required then we skip the AO integrals with 2p in the continuum.
 
             max_p = no_ao
             if (dont_two_p_continuum) then
                !rs pair is CC => we need only (TT|CC) AO integrals
                if (mobas%ao_basis%ordered_pairs(2,rs) .eq. 3) max_p = last_tgt_ao
             endif

             r = rs_ind(1,rs)
             s = rs_ind(2,rs)

             !Determine the ordered index of this AO pair
             rs_mapped = mobas%ao_basis%ordered_pairs(1,rs)
 
             iqrs = 0.0_cfp
             do p=1,max_p
                max_q = p
                if (dont_two_p_continuum) then
                   !rs pair is CT and p is C => we must load only (TT|CT), (CT|CT) AO integrals
                   if (mobas%ao_basis%ordered_pairs(2,rs) .eq. 2 .and. p > last_tgt_ao) max_q = last_tgt_ao
                endif
                ij = p*(p-1)/2
                !print *,max_p,max_q,ij
                !Preload the integrals (p,q) for the current p into an intermediate buffer hlp.
                do q=1,max_q
                   pq = ij + q

                   !Determine the ordered index of this AO pair
                   pq_mapped = mobas%ao_basis%ordered_pairs(1,pq)

                   ! WARNING: note that we have inlined here the indexing function for AOs (i.e. computation of 'i')!
                   ! We are not using the generic index function for performance reasons.
                   pq = max(pq_mapped,rs_mapped)
                   rs_tmp = min(pq_mapped,rs_mapped)

                   i = pq*(pq-1)/2+rs_tmp !standard indexing for TTTT CTCT CTTT classes and for 2p in the continuum

                   if (dont_two_p_continuum) then !Special indexing for CCTT class in case 1p in the continuum
                      if (pq .le. mobas%ao_basis%n_TT_pairs .or. rs_tmp .le. mobas%ao_basis%n_TT_pairs) then !is pq or rs TT pair?
                         if (pq > mobas%ao_basis%last_CT_fn) then !pq=CC, rs=TT
                            i = mobas%ao_basis%n_prec_ints + rs_tmp + mobas%ao_basis%n_TT_pairs*(pq-mobas%ao_basis%last_CT_fn-1)
                         elseif(rs_tmp > mobas%ao_basis%last_CT_fn) then !pq=TT, rs=CC pair
                            i = mobas%ao_basis%n_prec_ints + pq + mobas%ao_basis%n_TT_pairs*(rs_tmp-mobas%ao_basis%last_CT_fn-1)
                         endif
                      endif
                   endif
                   !if (pq .ge. rs) then
                   !   i = pq*(pq-1)/2+rs
                   !else
                   !   i = rs*(rs-1)/2+pq
                   !endif
                   if (i > size(ao_integrals%a,1)) then
                      print *,p,q,pq,rs_tmp,i
                   endif
                   hlp(q) = ao_integrals%a(i,int_type)
                   !write(100+myrank,'(4i4,e25.15)') p,q,r,s,hlp(q)
                enddo !q
 
                !Loop over the preloaded integrals and for each of them calculate its contribution to iqrs.
                do q=1,min(p-1,max_q)
                   if (hlp(q) .ne. 0.0_cfp) then
                      !We assume that the AO hlp(q) is symmetric: (pq|O|rs) = (qp|O|rs) so we calculate at once contributions of (pq|O|rs) to iqrs(i,q) and iqrs(i,p)
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the p-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                      do i = mobas%so2mo_range(1,q), mobas%so2mo_range(2,q) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,p) = iqrs(i,p) + hlp(q)*cf_t(i,q)
                      enddo !i
                   endif
                enddo !q
 
                if (max_q .eq. p) then
                   q = p
                   if (hlp(q) .ne. 0.0_cfp) then
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                   endif
                endif
             enddo !p
       
             no_int = no_int + no_pairs
       
             !2nd step - transform the second index: q->j
             !This step takes a little over half of the compute time for the whole first step (for no symmetry case)
             ij = 0
             p = rs-rs_start+1 !relative rs-index
             do i=1,mobas%number_of_functions
                !Here transpose one row of iqrs: this allows vectorization in the q-loop below and ensures cache locality.
                do q=1,no_ao
                   hlp(q) = iqrs(i,q)
                enddo
                do j=1,i
                   ij = ij + 1
                   !todo implement the same thing as in the not_local equivalent
                   if (dont_two_p_continuum) then
                      !Each pair, (i,j) or (r,s), is assigned a numeric value defining its type: TT=1, TC=2, CC=3.
                      !Within this scheme the pairs of types of [MO,AO] that can occur are: [1,1], [2,1], [1,2], [3,1], [1,3], [2,2].
                      !If we want integrals only for 1p in the continuum then the types that we want to skip have indices: [2,3], [3,2], [3,3].
                      !These pairs can be identified using the sum of their types: 5 and 6. This is what we test below.
                      ijrs_type = ij_type(ij)+mobas%ao_basis%ordered_pairs(2,rs)
                      if (ijrs_type .eq. 5 .or. ijrs_type .eq. 6) cycle
                   endif
                   !todo once AO symmetry has been implemented then I can check here that sym(i*j) = sym(r*s). If not then skip this j.
                   mo_int = 0.0_cfp
                   do q=mobas%mo2so_range(1,j), mobas%mo2so_range(2,j) !over all atomic functions which contribute to the j-th molecular orbital.
                      mo_int = mo_int + hlp(q)*cf(q,j)
                   enddo !q
                   ijrs(ij,p) = mo_int !(ij|O|rs), where i,j are MOs and r,s are AOs.
                enddo !j
             enddo !i
          enddo !rs
          !$OMP END DO
 
    end subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local

    !> \warning Note that the use of the 'allocatable' attribute for the argument arrays is key for performance since this
    !>          attribute allows the compiler to assume that the arrays are contiguous in memory.
    !> Alternatively the 'contiguous' attribute can be used but it is a F2008 feature so we omit it here.
    !> It is assumed that the threads have been launched outside of this routine.
    subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_not_local(ao_integrals,int_type,rs_start,rs_end,ij_type,&
                                iqrs,hlp,mobas,cf,cf_t,ijrs,no_ao,last_tgt_ao,no_int,no_pairs,rs_ind,two_p_continuum)
       use omp_lib
       use gto_routines, only: find_mapping, index_1p_continuum
       implicit none
       class(molecular_orbital_basis_obj), intent(in) :: mobas
       type(p2d_array_obj), intent(inout) :: ao_integrals
       integer, intent(in) :: int_type,rs_start,rs_end,no_ao,last_tgt_ao,no_pairs
       integer, allocatable :: rs_ind(:,:)
       integer(kind=1), allocatable :: ij_type(:)
       integer, intent(inout) :: no_int
       logical, intent(in) :: two_p_continuum
 
       !we use the allocatable attribute so that the compiler can assume that the arrays are unit-stride: this has ~20% impact on the performance of this routine.
       real(kind=cfp), allocatable :: iqrs(:,:), hlp(:), ijrs(:,:), cf(:,:), cf_t(:,:)
 
       integer :: rs, pq, p, q, ij, i, j, k, l, max_p, max_q, ijrs_type, iam, n_shell_r,n_shell_s,n_shell_p,n_shell_q
       integer :: p_rel,q_rel,r_rel,s_rel,ind_start_r,ind_start_s,ind_start_p,ind_start_q,i1,i2,i3,i4
       integer :: ind_orig(4), n(4), map(4), n_map(3), r_shell,s_shell,p_shell,q_shell,r,s
       real(kind=cfp) :: mo_int
       logical :: ao_contains_continuum, dont_two_p_continuum, rs_is_CC, rs_is_TT, pq_is_CC, pq_is_TT, is_CCTT

          if (mobas%ao_basis%n_cont_fns > 0) then
             ao_contains_continuum = .true.
          else
             ao_contains_continuum = .false.
          endif

          dont_two_p_continuum = ao_contains_continuum .and. .not.(two_p_continuum)

          !Split the rs-indices in the block among threads
          !todo this could be simplified by looping instead over all my quartets
          !of shells in ao_integrals%block_offset(:) and picking out only those
          ![pq|rs] where the pq or rs part falls within rs_start,rs_end
          !$OMP DO SCHEDULE(DYNAMIC)
          do rs=rs_start,rs_end !rs is the index of the (rs)-pair corresponding to the AO integrals (pq|rs) that we want to process here.
             iam = omp_get_thread_num()
       
             ! Pre-load all (available) AO integrals (pq|O|rs) (== (pq|rs)) with the given (rs) index into hlp array - this
             ! speeds up the 1st step a little bit compared with accessing ao_integrals%a in the j-loop below.
             ! We preload the integrals always only one row at a time so that we don't need an extra array of size ~no_ao**2.
    
             !1st step - transform the first index: p->i
             ! We loop over all AO integrals (pq|rs) with the given rs index: if integrals for only 1p in the continuum
             ! are required then we skip the AO integrals with 2p in the continuum.

             rs_is_CC = .false.
             rs_is_TT = .false.
             max_p = no_ao
             if (dont_two_p_continuum) then
                !rs pair is CC => we need only (TT|CC) AO integrals
                if (mobas%ao_basis%ordered_pairs(2,rs) .eq. 3) then
                   max_p = last_tgt_ao
                   rs_is_CC = .true.
                elseif (mobas%ao_basis%ordered_pairs(2,rs) .eq. 1) then
                   rs_is_TT = .true.
                endif
             endif

             r = rs_ind(1,rs)
             s = rs_ind(2,rs)

             r_shell = mobas%ao_basis%indices_to_shells(1,r) !index of the shell the r-th function is part of
             s_shell = mobas%ao_basis%indices_to_shells(1,s) !index of the shell the s-th function is part of
             k = max(r_shell,s_shell)
             l = min(r_shell,s_shell)

             r_rel = mobas%ao_basis%indices_to_shells(2,r) !index of the r-th function within the shell to which it belongs
             s_rel = mobas%ao_basis%indices_to_shells(2,s) !index of the s-th function within the shell to which it belongs
             ind_start_r = mobas%ao_basis%shell_descriptor(4,r_shell) !starting index for the shell the r-function is part of
             ind_start_s = mobas%ao_basis%shell_descriptor(4,s_shell) !starting index for the shell the s-function is part of
             n_shell_r = mobas%ao_basis%shell_descriptor(5,r_shell) !number of functions in the shell of which the r-function is part of 
             n_shell_s = mobas%ao_basis%shell_descriptor(5,s_shell) !number of functions in the shell of which the s-function is part of

             iqrs = 0.0_cfp
             do p=1,max_p
                max_q = p

                if (dont_two_p_continuum) then
                   !rs pair is CT and p is C => we must load only (TT|CT), (CT|CT) AO integrals
                   if (mobas%ao_basis%ordered_pairs(2,rs) .eq. 2 .and. p > last_tgt_ao) max_q = last_tgt_ao
                endif

                pq = p*(p-1)/2

                p_shell = mobas%ao_basis%indices_to_shells(1,p) !index of the shell the p-th function is part of
                ind_start_p = mobas%ao_basis%shell_descriptor(4,p_shell) !starting index for the shell the p-function is part of
                n_shell_p = mobas%ao_basis%shell_descriptor(5,p_shell) !number of functions in the shell of which the p-function is part of 
                p_rel = mobas%ao_basis%indices_to_shells(2,p) !index of the p-th function within the shell to which it belongs

                !Preload the integrals [pq|rs] for fixed p,r,s into the intermediate buffer hlp.
                do q=1,max_q
                   pq = pq + 1

                   q_shell = mobas%ao_basis%indices_to_shells(1,q) !index of the shell the q-th function is part of
                   i = max(p_shell,q_shell)
                   j = min(p_shell,q_shell)

                   pq_is_CC = .false.
                   pq_is_TT = .false.
                   is_CCTT = .false.
                   if (dont_two_p_continuum) then
                      if (mobas%ao_basis%ordered_pairs(2,pq) .eq. 1) then
                         pq_is_TT = .true.
                      elseif (mobas%ao_basis%ordered_pairs(2,pq) .eq. 3) then
                         pq_is_CC = .true.
                      endif
                      if (pq_is_TT .and. rs_is_CC) then
                         is_CCTT = .true.
                      elseif (pq_is_CC .and. rs_is_TT) then
                         is_CCTT = .true.
                      endif
                   endif

                   !i = index of the quartet of shells to which the [pq|rs] integral belongs
                   i = index_1p_continuum(mobas % ao_basis % ordered_shell_pairs, &
                                          i, j, k, l, is_CCTT, &
                                          mobas % ao_basis % last_CT_sh, &
                                          mobas % ao_basis % n_prec_sh, &
                                          mobas % ao_basis % n_TT_sh_pairs)

                   if (ao_integrals%block_offset(i) .eq. -1) then !The requested integral is not kept by this task
                      hlp(q) = 0.0_cfp
                      cycle
                   endif

                   !Each function p,q,r,s corresponds to a given shell. In the
                   !section below we permute the p,q,r,s indices to the standard order (determined by the starting indices of the functions in each shell)
                   !which was used to save the integrals (see atomic_orbital_basis_obj%two_electron_integrals).
                   ind_start_q = mobas%ao_basis%shell_descriptor(4,q_shell) !starting index for the shell the q-function is part of
                   n_shell_q = mobas%ao_basis%shell_descriptor(5,q_shell) !number of functions in the shell of which the q-function is part of
                   q_rel = mobas%ao_basis%indices_to_shells(2,q) !index of the q-th function within the shell to which it belongs

                   ind_orig(1:4) = (/ind_start_p,ind_start_q,ind_start_r,ind_start_s/)
                   n(1:4) = (/n_shell_p,n_shell_q,n_shell_r,n_shell_s/)

                   !Map the order of the corresponding p,q,r,s shells to the order in which the integrals within the shell are saved.
                   call find_mapping(ind_orig,n,n_map,map)

                   !Permute the relative indices into the order in which the integrals were saved
                   ind_orig(1:4) = (/p_rel,q_rel,r_rel,s_rel/)
                   i1 = ind_orig(map(1))
                   i2 = ind_orig(map(2))
                   i3 = ind_orig(map(3))
                   i4 = ind_orig(map(4))

                   !Compute the index of the integral [pq|rs] within its own
                   !quartet of shells and add to it the offset for the corresponding quartet of shells.
                   i = ao_integrals%block_offset(i) + i1 + n_map(1)*(i2-1) + n_map(2)*(i3-1) + n_map(3)*(i4-1)

                   if (i > size(ao_integrals%a,1)) then
                      print *,'indexing error:',p,q,r,s,i
                   endif
                   hlp(q) = ao_integrals%a(i,int_type)
                   !if (hlp(q) .ne. 0.0_cfp) write(100+myrank,'(4i4,e25.15)') p,q,r,s,hlp(q)
                enddo !q

                !Loop over the preloaded integrals and for each of them calculate its contribution to iqrs.
                do q=1,min(p-1,max_q)
                   if (hlp(q) .ne. 0.0_cfp) then
                      !We assume that the AO hlp(q) is symmetric: (pq|O|rs) = (qp|O|rs) so we calculate at once contributions of (pq|O|rs) to iqrs(i,q) and iqrs(i,p)
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the p-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                      do i = mobas%so2mo_range(1,q), mobas%so2mo_range(2,q) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,p) = iqrs(i,p) + hlp(q)*cf_t(i,q)
                      enddo !i
                   endif
                enddo !q
 
                if (max_q .eq. p) then
                   q = p
                   if (hlp(q) .ne. 0.0_cfp) then
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                   endif
                endif
             enddo !p
       
             no_int = no_int + no_pairs

             !2nd step - transform the second index: q->j
             !This step takes a little over half of the compute time for the whole first step (for no symmetry case)
             ij = 0
             p = rs-rs_start+1 !relative rs-index
             do i=1,mobas%number_of_functions
                !Here transpose one row of iqrs: this allows vectorization in the q-loop below and ensures cache locality.
                do q=1,no_ao
                   hlp(q) = iqrs(i,q)
                enddo
                do j=1,i
                   ij = ij + 1
                   !todo implement the same thing as in the not_local equivalent
                   if (dont_two_p_continuum) then
                      !Each pair, (i,j) or (r,s), is assigned a numeric value defining its type: TT=1, TC=2, CC=3.
                      !Within this scheme the pairs of types of [MO,AO] that can occur are: [1,1], [2,1], [1,2], [3,1], [1,3], [2,2].
                      !If we want integrals only for 1p in the continuum then the types that we want to skip have indices: [2,3], [3,2], [3,3].
                      !These pairs can be identified using the sum of their types: 5 and 6. This is what we test below.
                      ijrs_type = ij_type(ij)+mobas%ao_basis%ordered_pairs(2,rs)
                      if (ijrs_type .eq. 5 .or. ijrs_type .eq. 6) cycle
                   endif
                   !todo once AO symmetry has been implemented then I can check here that sym(i*j) = sym(r*s). If not then skip this j.
                   mo_int = 0.0_cfp
                   do q=mobas%mo2so_range(1,j), mobas%mo2so_range(2,j) !over all atomic functions which contribute to the j-th molecular orbital.
                      mo_int = mo_int + hlp(q)*cf(q,j)
                   enddo !q
                   ijrs(ij,p) = mo_int !(ij|O|rs), where i,j are MOs and r,s are AOs.
                enddo !j
             enddo !i

          enddo !rs
          !$OMP END DO
 
    end subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_not_local
 
    subroutine generate_ij_offset(mobas,ij_type,ij_orbital_range,two_p_continuum,ij_offset,n_integrals)
       use const, only: abel_prod_tab
       use mpi_mod
       use omp_lib
       use special_functions, only: ipair
       use sort, only: sort_int_float
       implicit none
       integer, allocatable :: ij_orbital_range(:,:), ij_offset(:)
       integer, intent(out) :: n_integrals
       integer(kind=1) :: ij_type(:)
       class(molecular_orbital_basis_obj) :: mobas
       logical, intent(in) :: two_p_continuum
 
       integer :: ij, sym_block, sym_i, sym_j, ij_irr, p, orb_i, orb_j, i, j, k, sym_k, sym_l, q, kl_irr
       integer :: orb_k_it, orb_k, orb_l, orb_kl, l, kl, orb_l_it, ijkl_type, offset, orb_ij, orb_j_start, orb_j_end, err
       logical :: case_iiii, case_iijj, case_ijij, case_ijkl 
       integer, allocatable :: tmp(:)

       ij_offset = 0
       ij = 0
       sym_block = 0
       do sym_i=1,mobas%no_irr
          do sym_j=1,sym_i
 
             ij_irr = abel_prod_tab(sym_i,sym_j) !IRR of the product of symmetries of orbitals i,j
             p = sym_i*(sym_i-1)/2 + sym_j       !the index of the symmetry block for the pair of symmetries corresponding to the pair (sym_i,sym_j), sym_i .g.e sym_j
 
             !We loop over those (i,j) orbital indices that have been assigned to this thread for the current combination of symmetries (sym_i,sym_j).
             do orb_i=ij_orbital_range(1,p),ij_orbital_range(3,p)
 
                orb_j_start = 1
                if (sym_i .eq. sym_j) then 
                   orb_j_end = orb_i          !both orbitals are from the same symmetry and therefore the loop over the second orbital must be only from 1 to orb_i.
                else
                   orb_j_end = mobas%orbital_data(sym_j)%number_of_functions !both orbitals come from different symmetries so the loop over the second loop must be over all orbitals in that symmetry.
                endif
 
                if (orb_i .eq. ij_orbital_range(1,p)) orb_j_start = ij_orbital_range(2,p) !my first 'j' orbital index
                if (orb_i .eq. ij_orbital_range(3,p)) orb_j_end = ij_orbital_range(4,p)   !my last 'j' orbital index
                
                do orb_j=orb_j_start,orb_j_end
 
                   i = max(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_i
                   j = min(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_j
                   ij = ipair(i) + j

                   do sym_k=1,sym_i
                      do sym_l=1,sym_k
                         q = sym_k*(sym_k-1)/2 + sym_l
       
                         if (q > p) cycle !symmetry of the integral (ij|O|kl) => we want only unique pairs of pairs of symmetries (p,q) with p=(sym_i,sym_j),q=(sym_k,sym_l)
                         kl_irr = abel_prod_tab(sym_k,sym_l) !IRR of the product of symmetries of orbitals k,l
       
                         !symmetry restriction assuming that the operator O in (ij|O|kl) is totally symmetric; hence the symmetry restriction is given only by the IRRs of the MOs.
                         if (kl_irr .ne. ij_irr) cycle
       
                         sym_block = p*(p-1)/2 + q !index of the unique quartet of symmetries (sym_i sym_j|sym_k sym_l) which defines the symmetry block
                         offset = mobas%block_offset(sym_block) !the total number of integrals preceeding the current symmetry block; this is the base for the index of the quartet of MOs
       
                         case_iiii = .false. 
                         case_iijj = .false.
                         case_ijij = .false.
                         case_ijkl = .false.
       
                         !determine the combination of symmetries we are dealing with; this is used below to ensure only unique combinations of orbitals are produced
                         if (sym_i .eq. sym_j .and. sym_k .eq. sym_l .and. sym_i .eq. sym_k) then
                            case_iiii = .true. !(II|II)
                         else if (sym_i .eq. sym_j .and. sym_k .eq. sym_l) then
                            case_iijj = .true. !(II|JJ)
                         else if (sym_i .eq. sym_k .and. sym_j .eq. sym_l) then
                            case_ijij = .true. !(IJ|IJ)
                         else
                            case_ijkl = .true. !(IJ|KL)
                         endif
       
                         !Compute indices of the (IJ| or (II| bra pairs so that we can compare them against the indices of the ket pairs and make sure orb_ij .ge. orb_kl.
                         !Case (II|: the index of a unique pair of orbitals is a simple triangularization
                         !Case (IJ|: the index of a pair of orbitals (each from a different symmetry) is the same as the sequence number of the element (orb_j,orb_i) in a linearized 2D array.
                         if (case_iiii .or. case_iijj) orb_ij = ipair(orb_i) + orb_j
                         if (case_ijij) orb_ij = orb_j + mobas%orbital_data(sym_j)%number_of_functions*(orb_i-1)
 
                         if (case_iiii .or. case_ijij) orb_k_it = orb_i           !we need (ij|kl) with i .ge. k
                         if (case_iijj .or. case_ijkl) orb_k_it = mobas%orbital_data(sym_k)%number_of_functions !k comes from a symmetry different to i so we have to loop over all orbitals in that symmetry.
 
                         do orb_k=1,orb_k_it
 
                            if (case_iiii .or. case_iijj) orb_l_it = orb_k          !we need (ij|kl) with k .ge. l
                            if (case_ijij .or. case_ijkl) orb_l_it = mobas%orbital_data(sym_l)%number_of_functions !l comes from a symmetry different to k so we have to loop over all orbitals in that symmetry.
 
                            do orb_l=1,orb_l_it
 
                               k = mobas%relative_to_absolute(orb_k,sym_k) !overall index of the orb_k orbital
                               l = mobas%relative_to_absolute(orb_l,sym_l) !overall index of the orb_l orbital
                               kl = ipair(max(k,l)) + min(k,l)
 
                               if (.not.(two_p_continuum)) then
                                  !If integrals for only 1p in the continuum are required then we skip integrals of the type (TC|CC), (CT|CC) and (CC|CC).
                                  ijkl_type = ij_type(ij)+ij_type(kl)
                                  if (ijkl_type .eq. 5 .or. ijkl_type .eq. 6) cycle
                               endif
 
                               !For (II|II) and (IJ|IJ) cases we need to make sure that ij .ge. kl: compute indices of the |IJ) or |II) ket pairs and compare them to the bra indices.
                               if (case_iiii) then
                                  orb_kl = ipair(orb_k) + orb_l
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               if (case_ijij) then 
                                  orb_kl = orb_l + mobas%orbital_data(sym_l)%number_of_functions*(orb_k-1)
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               !Order in which I will loop over the ijkl integrals: total number of integrals for each ij pair.
                               ij_offset(ij) = ij_offset(ij) + 1
 
                            enddo !orb_l
                         enddo !orb_k
 
                      enddo !sym_l
                   enddo !sym_k
                enddo !orb_j
             enddo !orb_i
          enddo !sym_j
       enddo !sym_i

       call move_alloc(ij_offset,tmp)
       allocate(ij_offset(size(tmp)),stat=err)
       if (err .ne. 0) call xermsg ('molecular_basis_mod','generate_ij_offset','Memory allocation failed.',err, 1)

       ij_offset = 0
       do ij=2,size(tmp)
          !Total number of integrals preceeding the current ij pair:
          ij_offset(ij) = ij_offset(ij-1) + tmp(ij-1)
       enddo

       !Total number of integrals generated by this thread
       n_integrals = ij_offset(size(tmp)) + tmp(size(tmp))

    end subroutine generate_ij_offset

    subroutine extract_non_zero_cf(cf_t,cf_t_non_zero,mo_indices,n_non_zero)
       implicit none
       real(kind=cfp), allocatable :: cf_t(:,:), cf_t_non_zero(:,:)
       integer, allocatable :: mo_indices(:,:), n_non_zero(:)

       integer :: i, j, p, n_mo, n_ao, err

          n_mo = size(cf_t,1)
          n_ao = size(cf_t,2)

          if (allocated(cf_t_non_zero)) deallocate(cf_t_non_zero)
          if (allocated(n_non_zero)) deallocate(n_non_zero)
          if (allocated(mo_indices)) deallocate(mo_indices)

          allocate(cf_t_non_zero(n_mo,n_ao),mo_indices(n_mo,n_ao),n_non_zero(n_ao),stat=err)
          if (err .ne. 0) call xermsg ('molecular_basis_mod','extract_non_zero_cf','Memory allocation failed.',err, 1)

          n_non_zero = 0
          cf_t_non_zero = 0.0_cfp
          mo_indices = 0
          do p=1,n_ao
             j = 0
             do i=1,n_mo
                if (cf_t(i,p) .ne. 0.0_cfp) then
                   j = j + 1
                   n_non_zero(p) = n_non_zero(p) + 1
                   cf_t_non_zero(j,p) = cf_t(i,p)
                   mo_indices(j,p) = i
                endif
             enddo !p
          enddo !i

    end subroutine extract_non_zero_cf
 
    !> The routine assumes that the array for output (ijkl_integrals) has been zeroed-out before the algorithm starts and that the threads have been launched outside of this routine.
    !> \warning Note the use of the 'allocatable' attributes for some argument arrays - this strongly affects performance as explained in the comment for omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local.
    subroutine omp_two_p_transform_ijrs_block_to_ijkl(ijrs,ijks,ijks_t,cf,cf_t_non_zero,mo_indices, &
                n_cf_t_non_zero,no_ao,mobas,rs_ind,ij_type,ao_is_local,tol,&
              &ijkl_integrals,ijkl_indices,int_type,rs_start,rs_end,no_pairs,ij_orbital_range,two_p_continuum,thread_id,ij_offset)
       use const, only: abel_prod_tab
       use mpi_mod
       use omp_lib
       use special_functions, only: ipair
       use sort, only: sort_int_float
       implicit none
       real(kind=cfp), allocatable :: ijrs(:,:), cf(:,:), cf_t_non_zero(:,:), ijks(:,:), ijks_t(:,:), ijkl_integrals(:,:,:)
       real(kind=cfp) :: tol
       integer, allocatable :: rs_ind(:,:), ij_orbital_range(:,:), ijkl_indices(:,:), ij_offset(:), &
                                mo_indices(:,:), n_cf_t_non_zero(:)
       integer, intent(in) :: int_type, rs_start, rs_end, no_ao, no_pairs, thread_id
       integer(kind=1) :: ij_type(:)
       class(molecular_orbital_basis_obj) :: mobas
       logical, intent(in) :: ao_is_local, two_p_continuum
 
       integer :: ij, sym_block, sym_i, sym_j, ij_irr, p, orb_i, orb_j, i, j, rs, r, s, k, sym_k, sym_l, q, &
                    kl_irr, orb_ij, orb_j_start, orb_j_end, a,b,t
       integer :: cnt, orb_k_it, orb_k, orb_l, orb_kl, l, kl, orb_l_it, rs_relative, ijrs_type, ijkl_type
       real(kind=cfp) :: mo_int, threshold
       logical :: case_iiii, case_iijj, case_ijij, case_ijkl, check_small

       !If we will be calculating contributions to the transformed integrals of the last set of rs-indices then we can delete the final integrals smaller than a given threshold.
       !However, this can only be done if all AO integrals are kept by each process - otherwise the accumulated transformed integrals come only from contributions of those AO integrals
       !that are kept by this process. This prevents me from deleting the MO integrals here - this can only be done once all MO integrals have been accumulated from all processes. 
       if (rs_end .eq. no_pairs .and. ao_is_local) then
          check_small = .true.
       else
          check_small = .false.
       endif

       !The value 'threshold' is used to neglect small ijrs integrals. The use of this threshold approximately halves the compute time for the
       !second step (at least for the case: pyrazine 6-311+G_dp/all HF orbitals). Note that we can use the 'tol' value only if this process keeps all
       !AO integrals (ao_is_local .eq. true.)! If the AO integrals are scattered among all processes then ijrs(:,:) contains only PARTIAL contributions to the full
       !integral so we cannot neglect these partial contributions on the 'tol' level here.
       if (ao_is_local) then
          threshold = tol
       else
          threshold = 0.0_cfp
       endif
 
       !II. THE SECOND PART OF THE INTEGRAL TRANSFORM:
       !We loop over all unique symmetry blocks and in each symmetry block we loop over all unique combinations of orbitals. We precompute, for each symmetry block, the ij part of the full integral 
       !index for (ij|kl). Only the kl part is computed in the inner-most loop. Note that the loops over the orbital indices orb_k and orb_l start with 1. These orb_ indices correspond to the internal
       !orbital indices as defined in molecular_orbital_basis_obj. These indices always start with 1 no matter what the actual, i.e. external, indices of the orbitals are. In other words the algorithm works
       !even for cases in which the external indices of the orbitals don't start with 1.
       ij = 0
       sym_block = 0
       do sym_i=1,mobas%no_irr
          do sym_j=1,sym_i
 
             ij_irr = abel_prod_tab(sym_i,sym_j) !IRR of the product of symmetries of orbitals i,j
             p = sym_i*(sym_i-1)/2 + sym_j       !the index of the symmetry block for the pair of symmetries corresponding to the pair (sym_i,sym_j), sym_i .g.e sym_j
 
             !We loop over those (i,j) orbital indices that have been assigned to this thread for the current combination of symmetries (sym_i,sym_j).
             do orb_i=ij_orbital_range(1,p),ij_orbital_range(3,p)
 
                orb_j_start = 1
                if (sym_i .eq. sym_j) then 
                   orb_j_end = orb_i          !both orbitals are from the same symmetry and therefore the loop over the second orbital must be only from 1 to orb_i.
                else
                   orb_j_end = mobas%orbital_data(sym_j)%number_of_functions !both orbitals come from different symmetries so the loop over the second loop must be over all orbitals in that symmetry.
                endif
 
                if (orb_i .eq. ij_orbital_range(1,p)) orb_j_start = ij_orbital_range(2,p) !my first 'j' orbital index
                if (orb_i .eq. ij_orbital_range(3,p)) orb_j_end = ij_orbital_range(4,p)   !my last 'j' orbital index
                
                do orb_j=orb_j_start,orb_j_end
 
                   i = max(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_i
                   j = min(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_j
                   ij = ipair(i) + j
 
                   !3rd step - transform the third index: r->k
                   ijks(:,:) = 0.0_cfp
                   do rs=rs_start,rs_end
 
                      if (.not.(two_p_continuum)) then
                         !If integrals for only 1p in the continuum are required then we skip half-transformed integrals of the type (TC|CC], (CT|CC] and (CC|CC].
                         !For explanation of ijrs_type see omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local
                         ijrs_type = ij_type(ij) + mobas%ao_basis%ordered_pairs(2,rs)
                         if (ijrs_type .eq. 5 .or. ijrs_type .eq. 6) cycle
                      endif
 
                      rs_relative = rs-rs_start +1 !the indexing in ijrs is relative to rs_start
                      r = rs_ind(1,rs)
                      s = rs_ind(2,rs)
                      mo_int = ijrs(ij,rs_relative) != (ij|O|rs]
                      !write(100+myrank,'(2i10,e25.15)') ij,rs,mo_int
 
                      !The encompassing if statement makes efficient use of the branch predictor so it should use 
                      !only a negligible amount of compute time.
                      if (ao_is_local) then
                         if (abs(mo_int) .le. threshold) cycle
                      endif
 
                      !todo once AO symmetry has been implemented I can check first if a given (ij|rs) integral, ijrs(rs,ij), contributes to the (ij|ks) integrals.
                      if (s .eq. r) then
!                         do k=mobas%so2mo_range(1,r), mobas%so2mo_range(2,r)
!                            ijks(k,s) = ijks(k,s) + mo_int*cf_t(k,r)
!                         enddo !k
                          do k=1,n_cf_t_non_zero(r)
                             a = mo_indices(k,r)
                             ijks(a,s) = ijks(a,s) + mo_int*cf_t_non_zero(k,r)
                          enddo !k
                      else
!                         do k=mobas%so2mo_range(1,r), mobas%so2mo_range(2,r)
!                            ijks(k,s) = ijks(k,s) + mo_int*cf_t(k,r)
!                         enddo !k
                         do k=1,n_cf_t_non_zero(r)
                             a = mo_indices(k,r)
                             ijks(a,s) = ijks(a,s) + mo_int*cf_t_non_zero(k,r)
                         enddo !k
!                         do k=mobas%so2mo_range(1,s), mobas%so2mo_range(2,s)
!                            ijks(k,r) = ijks(k,r) + mo_int*cf_t(k,s)
!                         enddo !k
                         do k=1,n_cf_t_non_zero(s)
                            a = mo_indices(k,s)
                            ijks(a,r) = ijks(a,r) + mo_int*cf_t_non_zero(k,s)
                         enddo !k
                      endif
                   enddo !rs
 
                   !4th step - transform the fourth index: s->l
                   call transpose_2d(ijks,ijks_t,mobas%number_of_functions,no_ao) !ijks_t = ijks**T: this pays off in cache locality (speed up) below
 
                   !Total number of integrals preceeding the current ij pair:
                   cnt = ij_offset(ij)
 
                   do sym_k=1,sym_i
                      do sym_l=1,sym_k
                         q = sym_k*(sym_k-1)/2 + sym_l
       
                         if (q > p) cycle !symmetry of the integral (ij|O|kl) => we want only unique pairs of pairs of symmetries (p,q) with p=(sym_i,sym_j),q=(sym_k,sym_l)
                         kl_irr = abel_prod_tab(sym_k,sym_l) !IRR of the product of symmetries of orbitals k,l
       
                         !symmetry restriction assuming that the operator O in (ij|O|kl) is totally symmetric; hence the symmetry restriction is given only by the IRRs of the MOs.
                         if (kl_irr .ne. ij_irr) cycle
       
                         sym_block = p*(p-1)/2 + q !index of the unique quartet of symmetries (sym_i sym_j|sym_k sym_l) which defines the symmetry block
       
                         case_iiii = .false. 
                         case_iijj = .false.
                         case_ijij = .false.
                         case_ijkl = .false.
       
                         !determine the combination of symmetries we are dealing with; this is used below to ensure only unique combinations of orbitals are produced
                         if (sym_i .eq. sym_j .and. sym_k .eq. sym_l .and. sym_i .eq. sym_k) then
                            case_iiii = .true. !(II|II)
                            orb_ij = ipair(orb_i) + orb_j
                            orb_k_it = orb_i           !we need (ij|kl) with i .ge. k
                         else if (sym_i .eq. sym_j .and. sym_k .eq. sym_l) then
                            case_iijj = .true. !(II|JJ)
                            orb_k_it = mobas%orbital_data(sym_k)%number_of_functions !k comes from a symmetry different to i so we have to loop over all orbitals in that symmetry.
                         else if (sym_i .eq. sym_k .and. sym_j .eq. sym_l) then
                            case_ijij = .true. !(IJ|IJ)
                            orb_ij = orb_j + mobas%orbital_data(sym_j)%number_of_functions*(orb_i-1)
                            orb_k_it = orb_i           !we need (ij|kl) with i .ge. k
                         else
                            case_ijkl = .true. !(IJ|KL)
                            orb_k_it = mobas%orbital_data(sym_k)%number_of_functions !k comes from a symmetry different to i so we have to loop over all orbitals in that symmetry.
                         endif
 
                         do orb_k=1,orb_k_it
 
                            if (case_iiii .or. case_iijj) orb_l_it = orb_k          !we need (ij|kl) with k .ge. l
                            if (case_ijij .or. case_ijkl) orb_l_it = mobas%orbital_data(sym_l)%number_of_functions !l comes from a symmetry different to k so we have to loop over all orbitals in that symmetry.
 
                            do orb_l=1,orb_l_it
 
                               k = mobas%relative_to_absolute(orb_k,sym_k) !overall index of the orb_k orbital
                               l = mobas%relative_to_absolute(orb_l,sym_l) !overall index of the orb_l orbital
                               kl = ipair(max(k,l)) + min(k,l)
 
                               if (.not.(two_p_continuum)) then
                                  !If integrals for only 1p in the continuum are required then we skip integrals of the type (TC|CC), (CT|CC) and (CC|CC).
                                  ijkl_type = ij_type(ij)+ij_type(kl)
                                  if (ijkl_type .eq. 5 .or. ijkl_type .eq. 6) cycle
                               endif
 
                               !For (II|II) and (IJ|IJ) cases we need to make sure that ij .ge. kl: compute indices of the |IJ) or |II) ket pairs and compare them to the bra indices.
                               if (case_iiii) then
                                  orb_kl = ipair(orb_k) + orb_l
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               if (case_ijij) then 
                                  orb_kl = orb_l + mobas%orbital_data(sym_l)%number_of_functions*(orb_k-1)
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               !In the other cases we just compute indices of the ket pairs.
                               if (case_iijj) orb_kl = ipair(orb_k) + orb_l
                               if (case_ijkl) orb_kl = orb_l+(orb_k-1)*mobas%orbital_data(sym_l)%number_of_functions
 
                               !transform the fourth index: s->l
                               !mo_int = (ij|kl), where i,j,k,l are the molecular orbitals.
                               a = mobas%mo2so_range(1,l); b = mobas%mo2so_range(2,l)
                               mo_int = sum(ijks_t(a:b,k)*cf(a:b,l))
                               !do s=mobas%mo2so_range(1,l), mobas%mo2so_range(2,l)
                               !   mo_int = mo_int + ijks_t(s,k)*cf(s,l)
                               !enddo !s
 
                               !INTEGRAL INDEXING FOLLOWS: this must be consistent with the index function
 
                               a = max(ij,kl)
                               b = min(ij,kl)
                               t = ipair(a) + b

                               cnt = cnt + 1
                               ijkl_indices(cnt,thread_id+1) = t
                               ijkl_integrals(cnt,thread_id+1,int_type) = ijkl_integrals(cnt,thread_id+1,int_type) + mo_int

                               !Delete final MO integrals smaller than threshold. Split the if statement into two to make efficient
                               !use of the branch predictor.
                               if (check_small) then
                                  if (abs(ijkl_integrals(cnt,thread_id+1,int_type)) < tol) then
                                     ijkl_integrals(cnt,thread_id+1,int_type) = 0.0_cfp
                                  end if
                               endif
          
                            enddo !orb_l
                         enddo !orb_k
 
                      enddo !sym_l
                   enddo !sym_k
                enddo !orb_j
             enddo !orb_i
          enddo !sym_j
       enddo !sym_i

    end subroutine omp_two_p_transform_ijrs_block_to_ijkl

    subroutine bisect_index(val,ijkl_indices,col,last_index,ind,found)
       implicit none
       integer, intent(in) :: val, last_index, col
       integer, pointer     :: ijkl_indices(:,:)
       integer, intent(out) :: ind
       logical, intent(out) :: found

       integer :: hlf, A, B, i, j

          found = .false.
          if (last_index .eq. 0) return

          if (last_index .le. 10) then
             do i=1,last_index
                if (ijkl_indices(i,col) .eq. val) then
                   found = .true.
                   ind = i
                   exit
                endif
             enddo !i
             return
          endif
          if (val > ijkl_indices(last_index,col) .or. val < ijkl_indices(1,col)) return

!          write(100+col,*) 'in bisect',val,last_index
!          write(100+col,*) 'max',ijkl_indices(1,col),ijkl_indices(last_index,col) 
!          do i=1,last_index
!             write(100+col,*) 'c',i,ijkl_indices(i,col)
!          enddo
          A = 1
          B = last_index
          i = 0
          do
             i = i + 1
             if (B .le. A+10) then
                do j=A,B
                   if (ijkl_indices(j,col) .eq. val) then
                      found = .true.
                      ind = j
                      exit
                   endif
                enddo !j
                exit
             endif
             hlf = A+(B-A)/2
!             write(100+col,'(3i5,3i10)') A,hlf,B,ijkl_indices(A,col),ijkl_indices(hlf,col),ijkl_indices(B,col)
             if (ijkl_indices(hlf,col) .eq. val) then
                found = .true.
                ind = hlf
                exit
             endif
             if (val < ijkl_indices(hlf,col)) then
                B = hlf-1
             elseif (val > ijkl_indices(hlf,col)) then
                A = hlf+1
             endif
!             if (ijkl_indices(A,col) .eq. val) then
!                found = .true.
!                ind = A
!                exit
!             endif
!             if (ijkl_indices(B,col) .eq. val) then
!                found = .true.
!                ind = B
!                exit
!             endif
!             if (B .eq. A+1) exit !further division is not possible: the index is not there
          enddo

!          write(100+col,*) 'end',found,val,i
                 
    end subroutine bisect_index

    subroutine write_ijkl_indices(this,lunit,record_start,position_after_write)
       use mpi_mod
       implicit none
       class(molecular_orbital_basis_obj) :: this
       integer, intent(in) :: lunit, record_start
       integer, intent(out) :: position_after_write

       integer :: err
       
         write(stdout,'("--------->","molecular_orbital_basis_obj:write_ijkl_indices")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'write_ijkl_indices', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         call mpi_mod_barrier(err)

         if (.not. associated(this % ijkl_indices)) then
            call xermsg ('molecular_orbital_basis_obj', 'write_ijkl_indices', &
                         'The ijkl_indices array has not been allocated.', 2, 1)
         end if

         if (myrank .eq. master) then
            write(lunit,pos=record_start,err=10) this%ind_ijkl_integral, size(this%ijkl_indices,2)
            write(lunit,err=10) this%ijkl_indices
            inquire(lunit,pos=position_after_write)
         endif

         !get the position_after_write on all processes
         call mpi_mod_bcast(position_after_write,master)

         write(stdout,'("<---------","molecular_orbital_basis_obj:write_ijkl_indices")')

         return

 10      call mpi_xermsg ('molecular_orbital_basis_obj', 'write_ijkl_indices', &
                          'Error executing the write command while writing the array data to the disk.', 2, 1)

    end subroutine write_ijkl_indices

    subroutine read_ijkl_indices(this,lunit,record_start,position_after_read)
       implicit none
       class(molecular_orbital_basis_obj) :: this
       integer, intent(in) :: lunit, record_start
       integer, intent(out) :: position_after_read
       integer :: err, d2, i,j
       integer(mpiint)  :: local_master_array(nprocs)

         write(stdout,'("--------->","molecular_orbital_basis_obj:read_ijkl_indices")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'read_ijkl_indices', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if
       
         call mpi_mod_barrier(err)

         if (this%shared_window_ijkl_indices /= -1) then
            d2 = size(this%ijkl_indices,2)
            call mpi_memory_deallocate_integer_2dim(this%ijkl_indices,this%ind_ijkl_integral*d2,this%shared_window_ijkl_indices)
            this%shared_window_ijkl_indices = -1
            !this%ijkl_indices => null()
         endif

         write(stdout,*) record_start,local_rank,lunit

         if (myrank .eq. master) then
            read(lunit,pos=record_start,err=50) this%ind_ijkl_integral, d2
         endif
         call mpi_mod_bcast(d2,master)
         call mpi_mod_bcast(this%ind_ijkl_integral,master)

         write(stdout,*) this%ind_ijkl_integral,d2 

         !All local masters allocate memory for the index array: in non-shared mode every process is the local master.
         this%shared_window_ijkl_indices = mpi_memory_allocate_integer_2dim(this%ijkl_indices,this%ind_ijkl_integral,d2)
         
         call mpi_memory_synchronize(this%shared_window_ijkl_indices)

         if (myrank .eq. master) then
            read(lunit,err=50) this%ijkl_indices
            inquire(lunit,pos=position_after_read)
         endif

         call mpi_memory_synchronize(this%shared_window_ijkl_indices)

         !get the position_after_read on all processes
         call mpi_mod_bcast(position_after_read,master)
         call mpi_mod_bcast(d2,master)

         if (shared_enabled) then !Every node has only one copy of the indexing array

            local_master_array = 1
            
            call mpi_mod_allgather(local_rank,local_master_array)
            
            call mpi_memory_synchronize(this%shared_window_ijkl_indices)

            if (myrank .eq. master) then
               do i=1,nprocs
                  if (i-1 .eq. master) cycle
                  if (local_master_array(i) == local_master) then 
                     write(stdout,"('Sending array to ',i4)") (i-1)
                     do j=1,d2
                        call mpi_mod_send(int(i-1,mpiint),this%ijkl_indices(:,j),1,this%ind_ijkl_integral)
                     enddo   
                  endif
               enddo
            endif
            
            if (local_rank == local_master .and. myrank /= master) then 
               do j=1,d2
                  write(stdout,"('Gathering info my array is size ',i12)") size(this%ijkl_indices,1)
                  call mpi_mod_recv(master,1,this%ijkl_indices(:,j),this%ind_ijkl_integral)
               enddo
            endif
    
            call mpi_memory_synchronize(this%shared_window_ijkl_indices)

         else !every MPI task keeps a copy of the indexing array

            call mpi_mod_barrier(err)

            do i=1,d2
                call mpi_mod_bcast(this%ijkl_indices(1:this%ind_ijkl_integral,d2),master)
            enddo !i

         endif

         call mpi_mod_barrier(err)
 
         write(stdout,'("<---------","molecular_orbital_basis_obj:read_ijkl_indices")')
 
         return
 
 50      call mpi_xermsg ('molecular_orbital_basis_obj', 'read_ijkl_indices', &
                          'Error executing the read command while reading the array data to the disk.', 2, 1)

    end subroutine read_ijkl_indices

    subroutine orbital_radial_charge_density(this,rmat_radius,A,B,delta_r,save_to_disk,charge_densities)
       use cgto_pw_expansions_mod
       use phys_const, only: fourpi
       use gto_routines, only: cms_gto_norm
       use const, only: line_len, fmat
       implicit none
       class(molecular_orbital_basis_obj) :: this
       real(kind=cfp), intent(in) :: rmat_radius,A,B,delta_r
       logical, intent(in) :: save_to_disk
       real(kind=cfp), allocatable :: charge_densities(:,:)
 
       integer :: i, j, k, err, n, s1, CGTO1_M, s2, CGTO2_M, ao1, ao2, BA_ind, n_shell_pairs, no_orbitals,&
                    lu, number_of_cgto_shells, max_bspline_l, max_prop_l, p
       integer :: bto1_index, bto2_index, BTO_M1, BTO_M2, number_of_bto_shells, lb, lbmb, cnt
       type(CGTO_shell_data_obj), allocatable :: dummy_cgto(:)
       type(BTO_shell_data_obj), allocatable :: dummy_bto(:)
       type(CGTO_shell_pair_pw_expansion_obj), allocatable :: CGTO_shell_pair_pw_expansion(:)
       type(CGTO_shell_pw_expansion_obj), allocatable :: CGTO_shell_pw_expansion(:)
       type(pw_expansion_obj) :: grid
       real(kind=cfp), allocatable :: orb_cf(:,:), bto_amplitude(:,:), jacobian_r(:)
       real(kind=cfp) :: fac, cf, norm, r
       character(len=line_len) :: file_name
       integer, parameter :: der = 0
 
          write(stdout,'("--------->","molecular_orbital_basis_obj:radial_charge_density")')
 
          if (.not. this % initialized) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', &
                          'The object has not been initialized or not all orbitals have been added.', 1, 1)
          end if
 
          if (A < 0.0_cfp) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input A < 0.', 2, 1)
          end if
          if (B <= 0.0_cfp) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input B .le. 0.', 3, 1)
          end if
          if (delta_r <= 0.0_cfp) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input delta_r .le. 0.', 4, 1)
          end if
          if (B < A) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input B < A.', 5, 1)
          end if

          !Print out distances of atoms from CMS: useful for analyzing the
          !radial charge density plots.
          write(stdout,'(/,"Atom distances from CMS")')
          do i=1,this%symmetry_data%no_nuc
             r = sqrt(dot_product(this%symmetry_data%nucleus(i)%center,this%symmetry_data%nucleus(i)%center))
             call this%symmetry_data%nucleus(i)%print
             write(stdout,'(5X,"Distance from CMS = ",e25.15," Bohr")') r
          enddo

          !Construct grid for [A,B],delta_r
          call grid%eval_regular_grid(A,B,delta_r)

          n = grid%n_total_points
          allocate(jacobian_r(n),stat=err)
          if (err /= 0) then
             call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 1 error.', err, 1)
          end if

          !write(stdout,'(/,"Radial grid:")')
          do i=1,n
             !write(stdout,'(e25.15)') grid%r_points(i)
             jacobian_r(i) = grid%r_points(i)**2
          enddo !i

          if (this%ao_basis%n_cgto_functions < this%ao_basis%number_of_functions) then
             call this%ao_basis%get_all_BTO_shells(dummy_bto,number_of_bto_shells)
             max_prop_l = 0
             max_bspline_l = maxval(dummy_bto(:)%l)
          else
             number_of_bto_shells = 0
             max_prop_l = 0
             max_bspline_l = -1
          endif

          if (this%ao_basis%n_cgto_functions > 0) then
 
             call this%ao_basis%get_all_CGTO_shells(dummy_cgto,number_of_cgto_shells)
   
             !Calculate radial charge_densities for all pairs of the CGTOs and multiply them in
             !with the orbital coefficients to obtain radial charge_densities of the orbitals.
             n_shell_pairs = number_of_cgto_shells*(number_of_cgto_shells+1)/2
             allocate(CGTO_shell_pair_pw_expansion(n_shell_pairs),CGTO_shell_pw_expansion(number_of_cgto_shells),stat=err)
             if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 3 error.', err, 1)
             end if
    
             do i=1,number_of_cgto_shells
                if (rmat_radius > 0.0_cfp .and. this%ao_basis%shell_descriptor(3,i) .eq. 1) then !Normalize the continuum to the R-matrix sphere
                   fac = cms_gto_norm(rmat_radius, dummy_cgto(i) % l, dummy_cgto(i) % number_of_primitives, &
                                      dummy_cgto(i) % exponents, dummy_cgto(i) % contractions, dummy_cgto(i) % norm, &
                                      dummy_cgto(i) % norms)
                   write(stdout,'("Continuum normalization factor",i0,e25.15)') i,fac
                   dummy_cgto(i)%norm = dummy_cgto(i)%norm*fac
                elseif (this%ao_basis%shell_descriptor(3,i) .eq. 1 .and. rmat_radius .le. 0.0_cfp) then
                   write(stdout,'("Continuum functions in this shell will not be normalized &
                                    &to the R-matrix sphere since rmat_radius .le. 0.0_cfp.")')
                endif
             enddo !i
    
             write(stdout,'("Precalculating radial amplitudes for all pairs of shells of CGTOs...")')
   
             call init_CGTO_pw_expansions_mod(0,maxval(dummy_cgto(:)%l))
   
             cnt = 0
             do i=1,number_of_cgto_shells
                do j=1,i
                   k = i*(i-1)/2+j
                   call CGTO_shell_pair_pw_expansion(k)%init_CGTO_shell_pair_pw_expansion(dummy_cgto(i),i,dummy_cgto(j),j)
                   call CGTO_shell_pair_pw_expansion(k)%assign_grid(grid%r_points,grid%weights)
                   call CGTO_shell_pair_pw_expansion(k)%eval_CGTO_shell_pair_pw_expansion
                   fac = k/real(n_shell_pairs)*100
                   if (mod(nint(fac),5) .eq. 0 .and. nint(fac) .ne. cnt) then
                      write(stdout,'(f8.3,"% done")') real(fac)
                      cnt = nint(fac)
                   endif
                enddo !j
             enddo !i
   
             if (number_of_cgto_shells < this%ao_basis%number_of_shells) then
                call init_CGTO_pw_expansions_mod(max_bspline_l,maxval(dummy_cgto(:)%l))
                do i=1,number_of_cgto_shells
                   call CGTO_shell_pw_expansion(i)%init_CGTO_shell_pw_expansion(dummy_cgto(i),i)
                   call CGTO_shell_pw_expansion(i)%assign_grid(grid%r_points,grid%weights)
                   call CGTO_shell_pw_expansion(i) % eval_CGTO_shell_pw_expansion (dummy_bto(1) % bspline_grid % knots, &
                                                                                   max_bspline_l, max_prop_l, 0)
                enddo !i
             endif
   
             write(stdout,'("...done")')
          else
             number_of_cgto_shells = 0
          endif

          if (number_of_cgto_shells < this%ao_basis%number_of_shells) then
             !Calculate amplitudes of the radial parts of the BTOs: B(r)/r
             allocate(bto_amplitude(n,dummy_bto(1)%bspline_grid%n),stat=err)
             if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 5 error.', err, 1)
             end if
             bto_amplitude = 0.0_cfp
             norm = 1.0_cfp
             do bto1_index=1,dummy_bto(1)%bspline_grid%n
                do j=1,n
                   bto_amplitude(j,bto1_index) = dummy_bto(1) % bspline_grid % bspline_amplitude(grid % r_points(j), &
                                                                                    norm, bto1_index, der) / grid % r_points(j)
                enddo !j
             enddo !bto1_index
          else
             number_of_bto_shells = 0
          endif

          !Copy the orbital coefficients to one array: this relies on the fact that the molecular orbitals are indexed symmetry by symmetry.
          call this%get_orbital_coefficient_matrix(orb_cf)
    
          !For each orbital calculate the radial charge_densities from the radial amplitudes of the GTOs.
          no_orbitals = size(orb_cf,2)
          if (allocated(charge_densities)) deallocate(charge_densities)
          allocate(charge_densities(n,no_orbitals),stat=err)
          if (err /= 0) then
             call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 5 error.', err, 1)
          end if
    
          charge_densities = 0.0_cfp
          do i=1,no_orbitals
             !Loop over all unique pairs of the CGTO shells
             ao1 = 0
             do s1=1,number_of_cgto_shells
                do CGTO1_M=1,dummy_cgto(s1)%number_of_functions
                   ao1 = ao1 + 1
                   ao2 = 0
                   !CGTO/CGTO contribution
                   do s2=1,s1
                      k = s1*(s1-1)/2+s2
                      do CGTO2_M=1,dummy_cgto(s2)%number_of_functions
                         ao2 = ao2 + 1

                         if (ao2 > ao1) cycle

                         BA_ind = CGTO2_M+dummy_cgto(s2)%number_of_functions*(CGTO1_M-1)
                         if (CGTO_shell_pair_pw_expansion(k)%neglect_m_lm(BA_ind,1)) cycle

                         fac = 2.0_cfp !off diagonal terms contribute twice since we loop only over the unique pairs of CGTOs
                         if (ao2 .eq. ao1) fac = 1.0_cfp !diagonal terms contribute only once
                         !the amplitudes of the CGTO pairs correspond to
                         !projections on X_{0,0} = 1/sqrt(4*pi) so we have to
                         !get rid of the pi-related factor.
                         cf = fac*orb_cf(ao1,i)*orb_cf(ao2,i)*sqrt(fourpi)
                         if (cf .eq. 0.0_cfp) cycle

                         charge_densities(1:n,i) = charge_densities(1:n,i) &
                                                + cf * CGTO_shell_pair_pw_expansion(k) % angular_integrals(1:n,BA_ind,1)

                      enddo !CGTO2_M
                   enddo !s2

                   !CGTO1/BTO contribution
                   fac = 2.0_cfp !off diagonal terms contribute twice since we loop only over the unique pairs of CGTO/BTO
                   ao2 = this%ao_basis%n_cgto_functions
                   do s2=1,number_of_bto_shells
                      lb = dummy_bto(s2)%l
                      bto1_index = dummy_bto(s2)%bspline_index
                      do BTO_M1 = -lb,lb
                         lbmb = lb*lb+lb+BTO_M1+1
                         ao2 = ao2 + 1

                         !ind = CGTO_shell_pw_expansion(s1)%non_neg_indices(CGTO1_M,1,lbmb)
                         !if (ind .eq. 0) cycle

                         !the amplitudes of the CGTOs correspond to projections
                         !on X_{l,m} of the BTO so there is no pi-related factor
                         !to get rid of but we have to multiply-in the BTO norm.
                         cf = fac*orb_cf(ao1,i)*orb_cf(ao2,i)*dummy_bto(s2)%norm
                         if (cf .eq. 0.0_cfp) cycle
                     
                         p = CGTO_shell_pw_expansion(s1)%non_neg_indices_l(CGTO1_M,lbmb)
                         if (p .eq. 0) cycle

                         charge_densities(1:n,i) = charge_densities(1:n,i) &
                                    + cf*CGTO_shell_pw_expansion(s1)%angular_integrals(1:n,p)*bto_amplitude(1:n,bto1_index)

                      enddo !BTO_M1
                   enddo !s2

                enddo !CGTO1_M
             enddo !s1

             !BTO/BTO contribution
             ao1 = this%ao_basis%n_cgto_functions
             do s1=1,number_of_bto_shells
                bto1_index = dummy_bto(s1)%bspline_index
                do BTO_M1=-dummy_bto(s1)%l,dummy_bto(s1)%l
                   ao1 = ao1 + 1
                   ao2 = this%ao_basis%n_cgto_functions
                   do s2=1,s1
                      bto2_index = dummy_bto(s2)%bspline_index
                      do BTO_M2=-dummy_bto(s2)%l,dummy_bto(s2)%l
                         ao2 = ao2 + 1
                         if (ao2 > ao1) cycle

                         !The only angular integrals which are non-zero are
                         !those where the angular parts of the BTOs are the same.
                         if ((dummy_bto(s1)%l .ne. dummy_bto(s2)%l) .or. (BTO_M1 .ne. BTO_M2)) cycle

                         fac = 2.0_cfp !off diagonal terms contribute twice since we loop only over the unique pairs of CGTOs
                         if (ao2 .eq. ao1) fac = 1.0_cfp !diagonal terms contribute only once
                         cf = fac*orb_cf(ao1,i)*orb_cf(ao2,i)*dummy_bto(s1)%norm*dummy_bto(s2)%norm
                         if (cf .eq. 0.0_cfp) cycle

                         charge_densities(1:n,i) = charge_densities(1:n,i) &
                                    + cf*bto_amplitude(1:n,bto1_index)*bto_amplitude(1:n,bto2_index)

                      enddo !BTO_M2
                   enddo !s2
                enddo !BTO_M1
             enddo !s1
 
             if (save_to_disk) then
                !Save the orbital density into a file with name: orb_rad_den_num.sym
                write(file_name,'(i6,".",i1)') this%absolute_to_relative(1,i),this%absolute_to_relative(2,i)
                if (this%is_continuum(i)) then
                   file_name = trim("orb_rad_den_continuum_"//adjustl(file_name))
                else
                   file_name = trim("orb_rad_den_target_"//adjustl(file_name))
                endif
                open(file=file_name,newunit=lu,status='replace',form=fmat,iostat=err)
                if (err /= 0) then
                    call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', &
                                 'Error opening the file for orbital charge density.', 6, 1)
                end if
 
                write(lu,'("#Radial charge density for orbital: ",3(i0,1x))') &
                            i, this%absolute_to_relative(1,i),this%absolute_to_relative(2,i)
                charge_densities(1:n,i) = charge_densities(1:n,i)*jacobian_r(1:n)
                fac = sum(charge_densities(1:n,i)*grid%weights(1:n))
                write(lu,'("#Integral over the charge density: ",e25.15)') fac
                do j=1,n
                   write(lu,'(2e25.15)') grid%r_points(j),charge_densities(j,i)
                enddo !j
                close(lu)
                write(stdout,'("Radial charge density for orbital ",i0," has been written to file: ",a)') &
                            i, trim(adjustl(file_name))
                write(stdout,'("Integral over the charge density: ",e25.15)') fac
             endif
          enddo !i
 
    end subroutine orbital_radial_charge_density

    subroutine eval_orbital(this,orb_i,r,n_points,orbital_at_r,sign_at_r)
       use cgto_pw_expansions_mod
       use gto_routines, only: cms_gto_norm
       use const, only: line_len, fmat
       implicit none
       class(molecular_orbital_basis_obj) :: this
       integer, intent(in) :: n_points, orb_i
       real(kind=cfp), intent(in) :: r(3,n_points)
       real(kind=cfp), allocatable :: orbital_at_r(:)
       integer, allocatable :: sign_at_r(:)

       integer :: i, j, n, irr, err
       real(kind=cfp), allocatable :: ao_basis_at_r(:,:)

          write(stdout,'("--------->","molecular_orbital_basis_obj:eval_orbital")')
 
          if (.not. this % initialized) then
             call xermsg ('molecular_orbital_basis_obj', 'eval_orbital', &
                          'The object has not been initialized or not all orbitals have been added.', 1, 1)
          end if

          call this%ao_basis%eval_basis(r,n_points,ao_basis_at_r)

          if (allocated(orbital_at_r)) deallocate(orbital_at_r)
          if (allocated(sign_at_r)) deallocate(sign_at_r)
          allocate(orbital_at_r(n_points),sign_at_r(n_points),stat=err)
          if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj','eval_orbital', 'Memory allocation 2 failed.', err, 1)

          !Contract the orbital coefficients with the values of the AO functions at the specified points.
          j = this%absolute_to_relative(1,orb_i)
          irr = this%absolute_to_relative(2,orb_i)
          n = this%ao_basis%number_of_functions
          do i=1,n_points
             orbital_at_r(i) = sum(this%orbital_data(irr)%coefficients(1:n,j)*ao_basis_at_r(1:n,i))
             sign_at_r(i) = nint(sign(1.0_cfp,orbital_at_r(i)))
          enddo !i

          write(stdout,'("<---------","molecular_orbital_basis_obj:eval_orbital")')

    end subroutine eval_orbital

end module molecular_basis_mod
