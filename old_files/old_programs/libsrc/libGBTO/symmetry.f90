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
module symmetry
use precisn
use utils, only: xermsg
use const, only: stdout, sym_op_nam_len
use common_obj, only: nucleus_type

   private

   public geometry_obj, symmetry_obj, determine_pg

   !> \class <geometry_obj>
   !> Input data for the symmetry_obj%init routine.
   type geometry_obj
      !> Number of symmetry operations. The maximum number of symmetry operations is 3.
      integer :: no_sym_op = 0
      !> Symmetry operations. These can be one of: X, Y, Z, XY, YZ, XZ, XYZ and determine which axes change sign under the symmetry operation.
      character(len=sym_op_nam_len) :: sym_op(1:3) = (/'  ','  ','  '/)
      !> Total number of nuclei.
      integer :: no_nuc = 0
      !> List of all nuclei. This array must be of size no_nuc.
      type(nucleus_type), allocatable :: nucleus(:)
      !> If set to true (default) then the nuclear symmetry information will be used to construct the list of symmetrically equivalent nuclei, etc. If set to false the molecular symmetry is still going to
      !> be determined and used, but the list of symmetrically equivalent nuclei, etc. will be constructed as if the molecule had no symmetry.
      !> \todo This switch will become obsolete once the list of symmetrically equivalent nuclei is always constructed (even if no_sym_op == 0). 
      logical :: use_symmetry = .true.
   contains
      !> Checks that the symmetry operations are consistent with the given list of nuclei.
      procedure :: check => check_geometry_obj
      !> Reads-in the symmetry input data from the file and position given.
      procedure :: read => read_geometry_obj
      !> Writes the symmetry input data into the file and position given.
      procedure :: write => write_geometry_obj
      !> This routine includes in the geometry data the continuum scattering centre.
      procedure :: add_scattering_centre => add_continuum
   end type geometry_obj

   !> \class <symmetry_obj>
   !> This structure holds all necessary symmetry information on the symmetry of the molecule that can be used in the integral calculation.
   !> \todo implement a member function calculating the nuclear repulsion energy
   !> \todo test symmetry determination properly.
   type symmetry_obj
      !> Integer identifier of the point group symmetry of the molecule. These are defined in the module const. This number is set by set_pg in case no_sym_op > 0.
      integer, private :: pg = 0
      !> Number of IRRs in the present point group. Set upon initialization.
      integer, private :: no_irrep = 0
      !> Number of symmetry operations. The maximum number of symmetry operations is 3.
      integer :: no_sym_op = 0
      !> Symmetry operations. These can be one of: X, Y, Z, XY, YZ, XZ, XYZ and determine which axes change sign under the symmetry operation.
      character(len=sym_op_nam_len) :: sym_op(1:3) = (/'  ','  ','  '/)
      !> If set to true (default) then the nuclear symmetry information will be used to construct the list of symmetrically equivalent nuclei, etc. If set to false the molecular symmetry is still going to
      !> be determined and used, but the list of symmetrically equivalent nuclei, etc. will be constructed as if the molecule had no symmetry.
      logical :: use_symmetry = .true.
      !> Number of the symmetrically non-redundant nuclei.
      integer :: no_sym_nuc = 0
      !> Total number of nuclei.
      integer :: no_nuc = 0
      !> List of all nuclei. This array must be of size no_nuc.
      type(nucleus_type), allocatable :: nucleus(:)
      !> List of the symmetrically non-redundant nuclei. This array must be of size no_sym_nuc
      type(nucleus_type), allocatable :: non_red_nuc(:)
      !> How many symmetrically equivalent partners each symmetrically non-redundant nucleus has. Must be of size no_sym_nuc.
      integer, allocatable :: no_eqv_nuc(:)
      !> List of the nuclei that are symmetrically equivalent to the nuclei from the array 'nucleus'. This array must be of the size (no_sym_nuc,no_nuc).
      type(nucleus_type), allocatable :: eqv_nuc(:,:)
      !> Names of the IRRs for the present orientation of the molecule.
      character(len=sym_op_nam_len), allocatable, private :: irr_names(:)
      !> Set to .true. after initialization
      logical, private :: initialized = .false. 
   contains
      !> Takes as an argument variable of the type geometry_obj, checks it and transfers the data here. Returns 0 if the init was successfull.
      procedure :: init => init_symmetry_obj
      !> This function checks that the values above are consistent with the sizes of the corresponding arrays, etc.
      !> \todo add checking for nuclei that are very close to the origin (dist < d1mach(4))
      procedure :: check => check_symmetry_obj
      !> Returns the irreducible representation of a given spherical harmonic in the point group symmetry of the molecule. It uses the point group symmetry identifier of the molecule and the list of 
      !> symmetry operations.
      !> \warning This note is aimed at developers. Before the routine set_irr_numbers is run sph_harm_pg_sym gives 0. 
      !> However, this should never be a problem since set_irr_numbers is always run upon calling init and sph_harm_pg_sym is not accessible before init is ran.
      procedure :: sph_harm_pg_sym => get_sph_harm_pg_sym
      !> This function returns the integer identifier of the point group symmetry corresponding to the list of the symmetry operations sym_op. It also sets the value pg.
      procedure :: get_pg => determine_pg_symmetry
      !> Given an integer number identifying the point group of the molecule, this function returns the number of irreducible representations of that point group.
      procedure :: get_no_irrep => return_no_irrep
      !> Sets the symmetry data: non_red_nuc, no_eqv_nuc, eqv_nuc
      procedure :: set_symmetry => set_symmetry_data
      !> Returns the symmetry data in the form of the geometry_obj data structure.
      procedure :: get_geometry => get_geometry_obj
      !> Given an integer number identifying the IRR this function returns the appropriate name of the IRR.
      procedure :: get_irr_name => get_name_of_irr
      !> This is run upon initialization to determine the sequence numbers of IRRs corresponding to the orientation of the molecule.
      procedure, private :: set_irr_numbers
   end type symmetry_obj

contains

   function check_geometry_obj(this)
      implicit none
      class(geometry_obj) :: this
      integer :: check_geometry_obj

      integer :: i, j, k
      character(len=sym_op_nam_len) :: tmp
      logical :: known_sym_op
      real(kind=cfp) :: test(1:3), test_sym_op(1:3), r(1:3), dist
      real(kind=cfp) :: cfp_dummy

         write(stdout,'("--------->","geometry_obj:check")')

         check_geometry_obj = 0

         if (this % no_nuc <= 0) then
            call xermsg ('symmetry', 'check_geometry_obj', 'Number of nuclei (no_nuc) is .le. 0.', 1, 1)
         end if

         if (.not. allocated(this % nucleus)) then
            call xermsg ('symmetry', 'check_geometry_obj', 'The nucleus array has not been allocated.', 2, 1)
         end if

         if (this % no_nuc /= size(this % nucleus)) then
            call xermsg ('symmetry', 'check_geometry_obj', 'Number of nuclei (no_nuc) .ne. size of the nucleus array.', 3, 1)
         end if

         if (this % no_sym_op < 0) then
            call xermsg ('symmetry', 'check_geometry_obj', 'Number of symmetry operations (no_sym_op) is < 0.', 4, 1)
         end if

         !check that we recognize all symmetry operations
         do i=1,this%no_sym_op
            known_sym_op = .false.
            tmp = trim(this%sym_op(i))
            if (tmp .eq. 'X') known_sym_op = .true.
            if (tmp .eq. 'Y') known_sym_op = .true.
            if (tmp .eq. 'Z') known_sym_op = .true.
            if (tmp .eq. 'XY') known_sym_op = .true.
            if (tmp .eq. 'YZ') known_sym_op = .true.
            if (tmp .eq. 'XZ') known_sym_op = .true.
            if (tmp .eq. 'XYZ') known_sym_op = .true.

            if (.not.(known_sym_op)) then
               print *,tmp
               call xermsg ('symmetry', 'check_geometry_obj', &
                            'An unknown symmetry element is present in the list of the symmetry operations.', 5, 1)
            endif

            !make sure that the symmetry element is not repeated
            do j=1,this%no_sym_op
               if (j .ne. i) then
                  if (tmp .eq. trim(this%sym_op(j))) then
                     call xermsg ('symmetry', 'check_geometry_obj', &
                                  'Each symmetry element can be present only once in the list of the symmetry elements.', 6, 1)
                  endif
               endif
            enddo
         enddo

         !check that the symmetry elements are consistent with the molecular geometry:
         !The symmetry operations are consistent with the geometry only if for each symmetrically non-redundant nucleus there exists a combination of symmetry elements which allows me to obtain 
         !from the symmetrically non-redundant ones all symmetrically equivalent nuclei.
         do j=1,this%no_nuc

            !coordinates of the nucleus with the symmetry operation this%sym_op(i) applied on them
            test(1:3) = this%nucleus(j)%center(1:3)

            do i=1,this%no_sym_op

               !test_sym_op is used to test whether each symmetry element generates at least one equivalent nucleus; if not then it is not an allowed symmetry operation for this molecule
               test_sym_op(1:3) = this%nucleus(j)%center(1:3)
   
               !apply the next symmetry element on test(1:3)
               call apply_sym_op(test,this%sym_op(i))

               !apply the symmetry operation on test_sym_op(1:3)
               call apply_sym_op(test_sym_op,this%sym_op(i))

               !test whether 'this%sym_op(i)' is an allowed symmetry element: the position given by test_sym_op must correspond to a nucleus of the molecule
               known_sym_op = .false.
               do k=1,this%no_nuc
                  r = this%nucleus(k)%center - test_sym_op
                  dist = sqrt(dot_product(r,r))
                  if (dist .le. f1mach(3,cfp_dummy) .and. this%nucleus(k)%name .eq. this%nucleus(j)%name) then
                     known_sym_op = .true.
                  endif
               enddo
               
               if (.not.(known_sym_op)) then
                  print *,this%sym_op(i)
                  call this%nucleus(j)%print
                  call xermsg('symmetry','check_geometry_obj','A symmetry operation is inconsistent with molecular geometry.',7,1)
               endif

            enddo !i 

         enddo !j

         write(stdout,'("<---------","done:geometry_obj:check")')

   end function check_geometry_obj

   subroutine read_geometry_obj(this,lunit,posit,pos_after_read)
      use mpi_mod
      implicit none
      class(geometry_obj) :: this
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: pos_after_read

      integer :: err, i

         write(stdout,'("--------->","geometry_obj:read")')

         pos_after_read = 0
         if (myrank .eq. master) then
            read(lunit,pos=posit,err=10) this%no_sym_op
            read(lunit,err=10) this%sym_op(1:3)
            read(lunit,err=10) this%no_nuc
   
            allocate(this%nucleus(1:this%no_nuc),stat=err)
            if (err .ne. 0) call xermsg('symmetry','read_geometry_obj', 'Memory allocation failed.',err,1)
   
            do i=1,this%no_nuc
               read(lunit,err=10) this%nucleus(i)%center(1:3), this%nucleus(i)%charge, this%nucleus(i)%nuc, this%nucleus(i)%name
            enddo
   
            read(lunit,err=10) this%use_symmetry
            inquire(lunit,pos=pos_after_read)
         endif

         !master broadcasts all its data to the other processes
         call mpi_mod_bcast(pos_after_read,master)
         call mpi_mod_bcast(this%no_sym_op,master)
         call mpi_mod_bcast(this%sym_op,master)
         call mpi_mod_bcast(this%no_nuc,master)

         if (myrank .ne. master) then
            allocate(this%nucleus(1:this%no_nuc),stat=err)
            if (err .ne. 0) call xermsg('symmetry','read_geometry_obj', 'Memory allocation 2 failed.',err,1)
         endif

         do i=1,this%no_nuc
            call mpi_mod_bcast(this%nucleus(i)%center(1:3),master)
            call mpi_mod_bcast(this%nucleus(i)%charge,master)
            call mpi_mod_bcast(this%nucleus(i)%nuc,master)
            call mpi_mod_bcast(this%nucleus(i)%name,master)
         enddo

         call mpi_mod_bcast(this%use_symmetry,master)

         err = this%check()
         if (err /= 0) then
            call xermsg ('symmetry', 'read_geometry_obj', 'geometry_obj data read-in but geometry_obj%check failed.', err, 1)
         end if

         write(stdout,'("<---------","done:geometry_obj:read")')

         return

 10      call xermsg('symmetry','read_geometry_obj','Error reading the geometry_obj data from the file and position given.',1,1)

   end subroutine read_geometry_obj

   subroutine write_geometry_obj(this,lunit,posit,pos_after_write)
      use mpi_mod
      implicit none
      class(geometry_obj) :: this
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: pos_after_write

      integer :: err, i

         write(stdout,'("--------->","geometry_obj:write")')

         err = this%check()
         if (err /= 0) then
            call xermsg ('symmetry', 'read_geometry_obj', 'geometry_obj%check failed. Erroneous data will not be written.', 1, 1)
         end if

         pos_after_write = 0
         if (myrank .eq. master) then
            write(lunit,pos=posit,err=10) this%no_sym_op
            write(lunit,err=10) this%sym_op(1:3)
            write(lunit,err=10) this%no_nuc

            do i=1,this%no_nuc
               write(lunit,err=10) this%nucleus(i)%center(1:3), this%nucleus(i)%charge, this%nucleus(i)%nuc, this%nucleus(i)%name
            enddo

            write(lunit,err=10) this%use_symmetry
            inquire(lunit,pos=pos_after_write)
         endif

         !master ensures all processes know where the record ends
         call mpi_mod_bcast(pos_after_write,master)

         write(stdout,'("<---------","done:geometry_obj:write")')

         return

 10      call xermsg('symmetry','write_geometry_obj','Error writing the geometry_obj data into the file and position given.',2,1)

   end subroutine write_geometry_obj

   function init_symmetry_obj(this,symmetry_data)
      implicit none
      class(symmetry_obj) :: this
      type(geometry_obj) :: symmetry_data
      integer :: init_symmetry_obj

      integer :: err, i

         write(stdout,'("--------->","symmetry_obj:init")')

         init_symmetry_obj = 0

         !check the input data
         err = symmetry_data%check()
         if (err /= 0) then
            call xermsg ('symmetry', 'init_symmetry_obj', &
                         'Check of input symmetry and nuclear data. failed; see geometry_obj%check for details.', err, 1)
         end if

         !transfer the nuclear data and symmetry elements
         if (allocated(this%nucleus)) deallocate(this%nucleus)
         allocate(this%nucleus(1:symmetry_data%no_nuc),stat=err)
         if (err .ne. 0) call xermsg('symmetry','init_symmetry_obj','Memory allocation failed.',err,1)

         do i=1,symmetry_data%no_nuc
            this%nucleus(i) = symmetry_data%nucleus(i)
         enddo

         this%no_nuc = symmetry_data%no_nuc
         this%no_sym_op = symmetry_data%no_sym_op
         this%sym_op = symmetry_data%sym_op
         this%use_symmetry = symmetry_data%use_symmetry

         this%initialized = .true.

         !Set the point group identifier and the number of IRRs
         i = this%get_pg()
         this%no_irrep = this%get_no_irrep(i)

         !Set the names of the IRRs (in Molpro convention
         call this%set_irr_numbers

         !Determine the symmetry information about the nuclei
         call this%set_symmetry(use_symmetry=this%use_symmetry)

         init_symmetry_obj = this%check()

         write(stdout,'("<---------","done:symmetry_obj:init")')
         
   end function init_symmetry_obj

   !> \warning Note that before this routine is called sph_harm_pg_sym gives always 0.
   subroutine set_irr_numbers(this)
      use const
      implicit none
      class(symmetry_obj) :: this

      integer :: irr_number, i, err, j, m
      character(len=sym_op_nam_len) :: name_of_irr

         if (allocated(this%irr_names)) deallocate(this%irr_names)
         allocate(this%irr_names(this%no_irrep),stat=err)
         if (err .ne. 0) call xermsg('symmetry','init_symmetry_obj','Memory allocation failed.',err,1)
         this%irr_names = ''

         this%irr_names = '?'

         !In Molpro the order of the IRRs is fixed regardless of the choice of orientation of the molecule wrt coordinate axes.
         select case (this%pg)
         case (C2v_id)
            this%irr_names = C2v_names
         case (D2h_id)
            this%irr_names = D2h_names
         case (C2h_id)
            this%irr_names = C2h_names
         case (D2_id)
            this%irr_names = D2_names
         case (C2_id)
            this%irr_names = C2_names
         case (Cs_id)
            this%irr_names = Cs_names
         case (Ci_id)
            this%irr_names = Ci_names
         case(C1_id)
            this%irr_names = C1_names
         case default
            call xermsg('symmetry','set_irr_numbers','Unknown point group identifier.',2,1)
         end select

         write(stdout,'("Order of IRRs is: ",8(i1,1X,a3,1X))') (i,this%irr_names(i),i=1,this%no_irrep)

   end subroutine set_irr_numbers

   function check_symmetry_obj(this)
      implicit none
      class(symmetry_obj) :: this
      integer :: check_symmetry_obj

      integer :: i, j, k, err, cnt1, cnt2

         write(stdout,'("--------->","symmetry_obj:check")')

         if (.not. this % initialized) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'The object has not been initialized.', 1, 1)
         end if

         check_symmetry_obj = 0

         if (this % no_sym_nuc <= 0 .or. this % no_nuc <= 0) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'Number of symmetrically non-redundant nuclei on input is .le. 0.', 1, 1)
         end if

         if (this % no_sym_nuc > this % no_nuc) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'Number of symmetrically non-redundant nuclei on input is greater than the total number of nuclei.', 0, 1)
         end if

         if (.not. allocated(this % non_red_nuc)) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'The nuclei have not been given any parameters.', 2, 1)
         end if

         if (size(this % non_red_nuc) /= this % no_sym_nuc) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'Inconsistence between number of symmetrically non-redundant nuclei and size of the corresp. array.', 3, 1)
         end if

         if (.not. allocated(this % no_eqv_nuc)) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'Array no_eqv_nuc has not been allocated.', 4, 1)
         end if

         if (size(this % no_eqv_nuc) /= this % no_sym_nuc) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'Incorrect size of the array no_eqv_nuc.', 5, 1)
         end if

         !check that the number of symmetrically equivalent nuclei and their parents give the total number of nuclei.
         if (sum(this % no_eqv_nuc) + this % no_sym_nuc /= this % no_nuc) then
            call xermsg ('symmetry', 'check_symmetry_obj', &
                         'The number of nuclei given in no_eqv_nuc is inconsistent with no_nuc.', 6, 1)
         end if

         !if no_sym_nuc .ne. no_nuc then the molecule has a symmetry and we expect a list of symmetrically equivalent nuclei on input.
         if (this%no_sym_nuc .ne. this%no_nuc) then
            if (.not. allocated(this % eqv_nuc)) then
                call xermsg ('symmetry', 'check_symmetry_obj', &
                             'The symmetrically equivalent nuclei have not been given any parameters.', 4, 1)
            end if
            if (size(this % eqv_nuc(:,1)) /= this % no_sym_nuc) then
                call xermsg ('symmetry', 'check_symmetry_obj', &
                             'Inconsistence between number of symmetrically equivalent nuclei &
                             &and size of the corresp. array.', 5, 1)
            end if
            if (size(this % eqv_nuc(1,:)) /= this % no_nuc) then
                call xermsg ('symmetry', 'check_symmetry_obj', &
                             'Inconsistence (2) between number of symmetrically equivalent nuclei *&
                             &and size of the corresp. array.', 6, 1)
            end if
         endif

         cnt1 = 0 !counts the number of nuclei centered on the CMS with Z=0, i.e. the nuclei for the basis of the continuum. Only one such center is allowed.
         cnt2 = 0 !counts the number of nuclei centered on the CMS with Z .ne. 0, i.e. genuine nuclei. Only one such nucleus is allowed.

         do i=1,this%no_sym_nuc
            !make sure that the nucleus itself is OK
            err = this%non_red_nuc(i)%check()
            if (err /= 0) then
                call xermsg ('symmetry', 'check_symmetry_obj', &
                             'Nucleus did not pass check; see check_symmetry_obj for error description.', err, 1)
            end if

!            if (this%non_red_nuc(i)%nuc == 0 .and. this%non_red_nuc(i)%charge == 0.0_cfp) cnt1 = cnt1 + 1
!            if (this%non_red_nuc(i)%nuc == 0 .and. this%non_red_nuc(i)%charge .ne. 0.0_cfp) cnt2 = cnt2 + 1
!            if (cnt1 > 1 .or. cnt2 > 1) call xermsg('symmetry','check_symmetry_obj','Too many nuclei centered on the CMS.',13,1)

            !now check that this nucleus has not been repeated in the list of the symmetrically non-redundant nuclei.
            do j=1,this%no_sym_nuc
               if (i .ne. j) then
                  !it is not allowed to repeat a nucleus if it is not centered on the CMS.
                  if (this % non_red_nuc(i) % nuc /= 0 .and. this % non_red_nuc(i) % nuc == this % non_red_nuc(j) % nuc) then
                      call xermsg ('symmetry', 'check_symmetry_obj', &
                                   'Repetition of nucleus in the list of symmetrically non-redundant nuclei.', 7, 1)
                  end if
               endif
            enddo

            !check the array of the equivalent nuclei for repetitions and errors
            write (stdout, '("no_eqv_nuc ALLOCATED:",L)') allocated(this%no_eqv_nuc)
            do j=1,this%no_eqv_nuc(i)
               !make sure that the nucleus itself is OK
               err = this%eqv_nuc(i,j)%check()
               if (err /= 0) then
                  call xermsg ('symmetry', 'check_symmetry_obj', &
                               'eqv_nuc did not pass check; see check_symmetry_obj for error description.', err, 1)
               end if

               if (this % eqv_nuc(i,j) % nuc == this % non_red_nuc(i) % nuc .and. this % non_red_nuc(i) % nuc /= 0) then
                  call xermsg ('symmetry', 'check_symmetry_obj', &
                               'A nucleus not sitting on the CMS cannot be symmetrically equivalent with itself.', 8, 1)
               end if
 
               if (this % eqv_nuc(i,j) % charge /= this % non_red_nuc(i) % charge) then
                  call xermsg ('symmetry', 'check_symmetry_obj', &
                               'Charge on a supposedly symmetrically equivalent nucleus .ne. to that of its &
                               &symmetrically non-redundant partner.', 10, 1)
               end if

               if (dot_product(this % eqv_nuc(i,j) % center, this % eqv_nuc(i,j) % center) /= &
                   dot_product(this % non_red_nuc(i) % center, this % non_red_nuc(i) % center)) then
                  call xermsg ('symmetry', 'check_symmetry_obj', &
                               'Distance from the origin of a supposedly symmetrically equivalent nucleus .ne. &
                               &to that of its symmetrically non-redundant partner.', 11, 1)
               end if

               !check whether the symetrically equivalent nucleus that is not centered on CMS has a position different from all symmetrically non-redundant nuclei.
               do k=1,this%no_sym_nuc
                  if (this % eqv_nuc(i,j) % nuc /= 0 .and. &
                      this % eqv_nuc(i,j) % center(1) == this % non_red_nuc(k) % center(1) .and. &
                      this % eqv_nuc(i,j) % center(2) == this % non_red_nuc(k) % center(2) .and. &
                      this % eqv_nuc(i,j) % center(3) == this % non_red_nuc(k) % center(3)) then
                      call xermsg ('symmetry', 'check_symmetry_obj', &
                                   'A symmetrically equivalent nucleus has the same position as its &
                                   &symmetrically non-redundant partner.', 12, 1)
                  end if

               enddo

            enddo
         enddo

         if (.not. allocated(this % nucleus)) then
            call xermsg ('symmetry', 'check_symmetry_obj', 'The number nucleus array has not been allocated.', 18, 1)
         end if

         if (size(this % nucleus) /= this % no_nuc) then
            call xermsg ('symmetry', 'check_symmetry_obj', 'The size of the nucleus array .ne. number of nuclei.', 19, 1)
         end if

         if (this % no_sym_op < 0 .or. this % no_sym_op > 3) then
            call xermsg ('symmetry', 'check_symmetry_obj', 'The number of symmetry elements is out of range.', 13, 1)
         end if

         if (this % no_sym_op > 0 .and. (this % pg <= 0 .or. this % pg > 8)) then
            call xermsg ('symmetry', 'check_symmetry_obj', 'The point group symmetry has not been set or is out of range.', 14, 1)
         end if

         if (this % no_irrep <= 0) then
            call xermsg ('symmetry', 'check_symmetry_obj', 'The value of no_irrep is either wrong or it has not been set.', 15, 1)
         end if

         if (size(this % irr_names) /= this % no_irrep) then
            call xermsg ('symmetry', 'check_symmetry_obj', 'Size of the array irr_names is wrong or no_irrep is wrong.', 16, 1)
         end if

         write(stdout,'("<---------","done:symmetry_obj:check")')

   end function check_symmetry_obj

   !> This routine determines the non-redundant nuclei. It should be run only after the symmetry has been thoroughly checked!
   subroutine set_symmetry_data(this,use_symmetry)
      implicit none
      class(symmetry_obj), intent(inout) :: this
      logical, optional, intent(in) :: use_symmetry

      integer :: i, j, k, err, non_red, no_nuc, ind
      real(kind=cfp) :: test(1:3), test_sym_op(1:3), r(1:3), dist
      logical :: known, symmetry
      logical, allocatable :: are_sym_eqv(:,:)
      integer, allocatable :: sym_eqv(:)
      real(kind=cfp) :: cfp_dummy

         write(stdout,'("--------->","symmetry_obj:set_symmetry")')

         if (.not.(this%initialized)) call xermsg('symmetry','set_symmetry_data','The object has not been initialized.',1,1)

         if (this%no_nuc .le. 0 .or. .not.allocated(this%nucleus)) then
            call xermsg('symmetry','set_symmetry_data','Nuclear data for symmetry processing are not ready.',2,1)
         endif

         symmetry = .true. !by default we will use the symmetry data to construct the list of symmetrically equivalent nuclei, etc.

         if (present(use_symmetry)) then !however if the user doesn't want it we will ignore it
            if (.not.(use_symmetry)) symmetry = .false.
         endif

         !are_sym_eqv(i,j) will be set to true below if the nuclei i,j are symmetrically equivalent.
         !sym_eqv(:) holds sequence numbers of the symmetrically non-redundant nuclei
         allocate(are_sym_eqv(1:this%no_nuc,1:this%no_nuc),sym_eqv(1:this%no_nuc),stat=err)
         if (err .ne. 0) call xermsg('symmetry','set_symmetry_data','Memory allocation 1 failed.',err,1)
         are_sym_eqv(:,:) = .false.
         sym_eqv(:) = -1

         !determine the symmetrically equivalent nuclei
         do j=1,this%no_nuc

            !coordinates of the nucleus with the symmetry operation this%sym_op(i) applied on them
            test(1:3) = this%nucleus(j)%center(1:3)

            !the continuum scattering center never has a symmetrical partner and that is OK
            if (this%nucleus(j)%is_continuum()) then
               write(stdout,'("Nucleus ",i2," is the scattering center.")') this%nucleus(j)%nuc
            endif

            do i=1,this%no_sym_op

               !test_sym_op is used to test whether each symmetry element generates at least one equivalent nucleus; if not then it is not an allowed symmetry operation for this molecule
               test_sym_op(1:3) = this%nucleus(j)%center(1:3)
   
               !apply the next symmetry element on test(1:3)
               call apply_sym_op(test,this%sym_op(i))

               !apply the symmetry operation on test_sym_op(1:3)
               call apply_sym_op(test_sym_op,this%sym_op(i))

               !find out whether we find a symmetrical partner for nucleus j that is sitting on the position given by test(1:3)
               do k=j+1,this%no_nuc
                  !first test application of the symmetry operation this%sym_op(i) alone using test_sym_op(1:3)
                  r = this%nucleus(k)%center - test_sym_op
                  dist = sqrt(dot_product(r,r))
                  if (dist .le. f1mach(3,cfp_dummy) .and. this%nucleus(k)%name .eq. this%nucleus(j)%name) then
                     write(stdout,'("Nuclei ",i2," and ",i2," are symmetrically equivalent.")') &
                        this%nucleus(j)%nuc, this%nucleus(k)%nuc
                     !set this pair of nuclei as symmetrically equivalent:
                     are_sym_eqv(j,k) = .true.
                     are_sym_eqv(k,j) = .true.
                  else !test the current combination of the symmetry elements in 'test'
                     r = this%nucleus(k)%center - test
                     dist = sqrt(dot_product(r,r))
                     if (dist .le. f1mach(3,cfp_dummy) .and. this%nucleus(k)%name .eq. this%nucleus(j)%name) then
                        write(stdout,'("Nuclei ",i2," and ",i2," are symmetrically equivalent.")') &
                            this%nucleus(j)%nuc, this%nucleus(k)%nuc
                        !set this pair of nuclei as symmetrically equivalent:
                        are_sym_eqv(j,k) = .true.
                        are_sym_eqv(k,j) = .true.
                     endif
                  endif
               enddo !k
               
            enddo !i 

         enddo !j

         !determine the symmetrically non-redundant nuclei
         non_red = 0 !number of symmetrically non-redundant nuclei
         do j=1,this%no_nuc

            known = .false. !is this nucleus on the list of symmetrically non-redundant nuclei?
            do i=1,non_red
               if (sym_eqv(i) .eq. j) known = .true.
            enddo

            if (.not.(known)) then !check whether this nucleus is not symmetrically equivalent with one of the nuclei already on the list

               known = .false.
               do i=1,non_red
                  if (are_sym_eqv(sym_eqv(i),j)) known = .true. !yes the nucleus j is symmetrically equivalent to one of the symmetrically non-redundant nuclei on the list
               enddo

               if (.not.(known)) then !add this nucleus on the list of symmetrically non-redundant nuclei
                  non_red = non_red + 1
                  sym_eqv(non_red) = j
                  write(stdout,'("Nucleus ",i2,1x,a2," is symmetrically non-redundant.")') this%nucleus(j)%nuc, this%nucleus(j)%name
               endif
            endif
         enddo

         !make sure that the symmetrically non-redundant nuclei are indeed linked to all atoms in the molecule
         no_nuc = 0 !total number of nuclei
         do i=1,non_red
            ind = sym_eqv(i)
            no_nuc = no_nuc + count(are_sym_eqv(ind,ind+1:this%no_nuc)) !sum of .true. values in the upper triangle of are_sym_eqv
         enddo
         no_nuc = no_nuc + non_red !add the symmetrically non-redundant nuclei and we should get the number of atoms in the molecule

         if (no_nuc .ne. this%no_nuc) then
            print *,no_nuc,this%no_nuc
            call xermsg ('symmetry', 'set_symmetry_data', &
                         'The symmetry information is incomplete: some nuclei cannot be reached using &
                         &the given symmetry operations.', 3, 1)
         endif

         !Build the nuclear symmetry information
         if (symmetry) then !we want to use the symmetry data
            this%no_sym_nuc = non_red !number of symmetrically non-redundant centers.
   
            !get rid of possilbe old symmetry data and allocate space for the nuclei:
            if (allocated(this%non_red_nuc)) deallocate(this%non_red_nuc)
            if (allocated(this%eqv_nuc)) deallocate(this%eqv_nuc)
            if (allocated(this%no_eqv_nuc)) deallocate(this%no_eqv_nuc)
            allocate(this % non_red_nuc(1 : this % no_sym_nuc), &
                     this % no_eqv_nuc(1 : this % no_sym_nuc), &
                     this % eqv_nuc(1 : this % no_sym_nuc, 1 : this % no_nuc), stat = err)
            if (err .ne. 0) call xermsg('symmetry','set_symmetry_data','Memory allocation 1 failed.',err,1)
   
            !construct the data on the symmetrically equivalent nuclei
            this%no_eqv_nuc(:) = 0
            do i=1,non_red
   
               !the list of the non_redundant nuclei
               ind = sym_eqv(i)
               this%non_red_nuc(i) = this%nucleus(ind)
               this%no_eqv_nuc(i) = count(are_sym_eqv(ind,ind+1:this%no_nuc)) !how many symmetrically equivalent nuclei each non-redundant has
               write(stdout,'("Non-redundant nucleus ",i2,1x,a," has ",i2," symmetrically equivalent nuclei.")') &
                    this % nucleus(sym_eqv(i)) % nuc, this % nucleus(sym_eqv(i)) % name, this % no_eqv_nuc(i)
   
               !transfer the nuclei this%nucleus(:) that are equivalent to this%nucleus(sym_eqv(i)) to the this%eqv_nuc(i,:) array
               k = 0
               do j=ind+1,this%no_nuc
                  if (are_sym_eqv(ind,j)) then
                     k = k + 1
                     this%eqv_nuc(i,k) = this%nucleus(j)
                     write(stdout,'(5x,"Equivalent nucleus: ",i2,1x,a)') this%nucleus(j)%nuc, this%nucleus(j)%name
                  endif
               enddo
   
               if (k /= this % no_eqv_nuc(i)) then
                  call xermsg ('symmetry', 'set_symmetry_data', 'Error in the symmetry information transfer.', 4, 1)
               end if
   
            enddo

         else !the user wants to ingore the symmetry data

            write(stdout,'("The target symmetry will be ignored for construction of the list of symmetrically equivalent nuclei.")')

            !Transfer the atom data into the 'this' structure. Note that we transfer the atom data into the symmetry_data as if no molecular symmetry was present.
            this%no_sym_nuc = this%no_nuc !number of symmetrically non-redundant centers. We always choose this%no_atoms, i.e. no symmetry at all.
   
            !get rid of possilbe old symmetry data and allocate space for the nuclei:
            if (allocated(this%non_red_nuc)) deallocate(this%non_red_nuc)
            if (allocated(this%eqv_nuc)) deallocate(this%eqv_nuc)
            if (allocated(this%no_eqv_nuc)) deallocate(this%no_eqv_nuc)
            allocate(this % non_red_nuc(1 : this % no_sym_nuc), &
                     this % no_eqv_nuc(1 : this % no_sym_nuc), &
                     this % eqv_nuc(1 : this % no_sym_nuc, 1 : this % no_nuc), stat = err)
            if (err .ne. 0) call xermsg('symmetry','set_symmetry_data','Memory allocation 1 failed.',err,1)
   
            this%no_eqv_nuc(:) = 0          !we don't specify the molecular symmetry here, i.e. all nuclei have 0 symmetrically equivalent partners
   
            !finally, copy the atom data to the non_red_nuc structure
            do i=1,this%no_nuc
               this%non_red_nuc(i) = this%nucleus(i)
            enddo

          endif

         write(stdout,'("<---------","done:symmetry_obj:set_symmetry")')
      
   end subroutine set_symmetry_data

   !> This routine takes the symmetry data which is assumed to describe only symmetry of the target molecule and includes in it symmetry data for the continuum scattering center.
   !> If the continuum center is already in the list then the output structure is identical to symmetry_data on input.
   !> \todo Test that the molecular symmetry data is preserved correctly.
   subroutine add_continuum(this)
      use const, only: nam_scattering_centre, id_scattering_centre
      implicit none
      class(geometry_obj) :: this

      integer :: i, err
      type(nucleus_type), allocatable :: nucleus(:)

         write(stdout,'("--------->","geometry_obj:add_scattering_centre")')

         !Search for the continuum center. If it is already there then quit.
         do i=1,this%no_nuc
            if (trim(this%nucleus(i)%name) .eq. nam_scattering_centre .and. this%nucleus(i)%nuc .eq. 0) return
         enddo

         if (this%no_nuc > 0) then !save the existing atom data
            allocate(nucleus(this%no_nuc),stat=err)
            if (err .ne. 0) call xermsg('geometry_obj','add_continuum','Memory allocation 1 failed.',err,1)
            nucleus(1:this%no_nuc) = this%nucleus(1:this%no_nuc)
         endif

         !make space for the continuum center
         if (allocated(this%nucleus)) deallocate(this%nucleus)
         this%no_nuc = this%no_nuc + 1
         allocate(this%nucleus(1:this%no_nuc),stat=err)
         if (err .ne. 0) call xermsg('geometry_obj','add_continuum','Memory allocation 2 failed.',err,1)

         if (this%no_nuc > 0) then !transfer the existing atom data
            this%nucleus(1:this%no_nuc-1) = nucleus(1:this%no_nuc-1)
         endif

         !add data for the scattering center
         this%nucleus(this%no_nuc)%center = (/0.0_cfp,0.0_cfp,0.0_cfp/)
         this%nucleus(this%no_nuc)%charge = 0.0_cfp
         this%nucleus(this%no_nuc)%nuc = id_scattering_centre
         this%nucleus(this%no_nuc)%name = nam_scattering_centre

         write(stdout,'("<---------","done:geometry_obj:add_scattering_centre")')

   end subroutine add_continuum

   !> This routine applies a symmetry operation given by the string sym_op on the vector r.
   subroutine apply_sym_op(r,sym_op)
      implicit none
      real(kind=cfp), intent(inout) :: r(1:3)
      character(len=sym_op_nam_len), intent(in) :: sym_op

      character(len=sym_op_nam_len) :: tmp

         tmp = trim(sym_op)
         if (tmp .eq. 'X') then
            r(1) = -r(1)
            return
         endif
         if (tmp .eq. 'Y') then
            r(2) = -r(2)
            return
         endif
         if (tmp .eq. 'Z') then
            r(3) = -r(3)
            return
         endif
         if (tmp .eq. 'XY') then
            r(1) = -r(1)
            r(2) = -r(2)
            return
         endif
         if (tmp .eq. 'YZ') then
            r(2) = -r(2)
            r(3) = -r(3)
            return
         endif
         if (tmp .eq. 'XZ') then
            r(1) = -r(1)
            r(3) = -r(3)
            return
         endif
         if (tmp .eq. 'XYZ') then
            r(1) = -r(1)
            r(2) = -r(2)
            r(3) = -r(3)
            return
         endif

         print *,tmp
         call xermsg('symmetry','apply_sym_op','An unknown symmetry element on input.',1,1)
      
   end subroutine apply_sym_op

   !> \todo add checking of 'this'.
   function determine_pg_symmetry(this,print_to_stdout)
   use const
      implicit none
      class(symmetry_obj) :: this
      integer :: determine_pg_symmetry
      logical, optional :: print_to_stdout

         if (.not.(this%initialized)) call xermsg('symmetry','determine_pg_symmetry','The object has not been initialized.',1,1)

         this % pg = determine_pg(this % no_sym_op, this % sym_op, print_to_stdout)
         determine_pg_symmetry = this % pg

   end function determine_pg_symmetry


   !> \brief   Find point group by symmetry operations
   !> \authors Zdenek Masin
   !> \date    2016
   !>
   !> Returns the group identifier as defined in "const.f90", based on the set of symmetry operations
   !> passed as arguments (each one of 'X', 'Y', 'Z', 'XY', 'YZ', 'XZ', 'XYZ', indicating simple and
   !> combined plane reflections that produce no extra sign).
   !>
   !> \param no_sym_op        Number of symmetry operations set in \c sym_op.
   !> \param sym_op           A triplet of strings, of which only \c no_sym_op need to be set.
   !> \param print_to_stdout  Print the group name to stdout.
   !>
   function determine_pg (no_sym_op, sym_op, print_to_stdout)
      use const
      integer :: no_sym_op, determine_pg
      logical, optional :: print_to_stdout
      character(len=sym_op_nam_len) :: sym_op(3), tmp, tmp_1, tmp_2, tmp_3
      logical :: t1_is_x_y_z, t1_is_xy_yz_xz, t2_is_x_y_z, t2_is_xy_yz_xz, t3_is_x_y_z, print_to_screen

         print_to_screen = .false.
         if (present(print_to_stdout)) print_to_screen = print_to_stdout

         determine_pg = 0

         !determine the symmetry based on the number of symmetry operations

         if (no_sym_op == 0) then
            if (print_to_screen) write(stdout,'(/,"C1 symmetry")')
            determine_pg = C1_id
            return
         endif

         if (no_sym_op == 1) then
            tmp = trim(sym_op(1))
            if (tmp .eq. 'X' .or. tmp .eq. 'Y' .or. tmp .eq. 'Z') then
               if (print_to_screen) write(stdout,'(/,"Cs symmetry")')
               determine_pg = Cs_id
               return
            endif
            if (tmp .eq. 'XY' .or. tmp .eq. 'YZ' .or. tmp .eq. 'XZ') then
               if (print_to_screen) write(stdout,'(/,"C2 symmetry")')
               determine_pg = C2_id
               return
            endif
            if (tmp .eq. 'XYZ') then
               if (print_to_screen) write(stdout,'(/,"Ci symmetry")')
               determine_pg = Ci_id
               return
            endif
         endif

         if (no_sym_op == 2) then
            tmp_1 = trim(sym_op(1))
            tmp_2 = trim(sym_op(2))

            t1_is_x_y_z = .false.
            t2_is_x_y_z = .false.
            if (tmp_1 .eq. 'X' .or. tmp_1 .eq. 'Y' .or. tmp_1 .eq. 'Z') t1_is_x_y_z = .true.
            if (tmp_2 .eq. 'X' .or. tmp_2 .eq. 'Y' .or. tmp_2 .eq. 'Z') t2_is_x_y_z = .true.
            t1_is_xy_yz_xz = .false.
            t2_is_xy_yz_xz = .false.
            if (tmp_1 .eq. 'XY' .or. tmp_1 .eq. 'YZ' .or. tmp_1 .eq. 'XZ') t1_is_xy_yz_xz = .true.
            if (tmp_2 .eq. 'XY' .or. tmp_2 .eq. 'YZ' .or. tmp_2 .eq. 'XZ') t2_is_xy_yz_xz = .true.

            !test for C2v
            if (t1_is_x_y_z .and. t2_is_x_y_z) then
               if (print_to_screen) write(stdout,'(/,"C2v symmetry")')
               determine_pg = C2v_id
               return
            endif

            !test for C2h
            if ((t1_is_x_y_z .and. t2_is_xy_yz_xz) .or. (t1_is_xy_yz_xz .and. t2_is_x_y_z)) then
               if (print_to_screen) write(stdout,'(/,"C2h symmetry")')
               determine_pg = C2h_id
               return
            endif

            !test for D2
            if (t1_is_xy_yz_xz .and. t2_is_xy_yz_xz) then
               if (print_to_screen) write(stdout,'(/,"D2 symmetry")')
               determine_pg = D2_id
               return
            endif
         endif

         if (no_sym_op == 3) then
            tmp_1 = trim(sym_op(1))
            tmp_2 = trim(sym_op(2))
            tmp_3 = trim(sym_op(3))

            t1_is_x_y_z = .false.
            t2_is_x_y_z = .false.
            t3_is_x_y_z = .false.
            if (tmp_1 .eq. 'X' .or. tmp_1 .eq. 'Y' .or. tmp_1 .eq. 'Z') t1_is_x_y_z = .true.
            if (tmp_2 .eq. 'X' .or. tmp_2 .eq. 'Y' .or. tmp_2 .eq. 'Z') t2_is_x_y_z = .true.
            if (tmp_3 .eq. 'X' .or. tmp_3 .eq. 'Y' .or. tmp_3 .eq. 'Z') t3_is_x_y_z = .true.

            !test for D2h
            if (t3_is_x_y_z .and. t1_is_x_y_z .and. t2_is_x_y_z) then
               if (print_to_screen) write(stdout,'(/,"D2h symmetry")')
               determine_pg = D2h_id
               return
            endif
         endif

         !we'll reach this bit only if the point group has not been identified
         call xermsg ('symmetry', 'determine_pg', &
                      'The input symmetry operations are wrong or do not correspond to any of the supported point groups.', 2, 1)

   end function determine_pg


   function return_no_irrep(this,pg)
   use const
      implicit none
      class(symmetry_obj) :: this
      integer, intent(in) :: pg
      integer :: return_no_irrep

         return_no_irrep = 0

         if (pg .le. 0 .or. pg > 8) then
            call xermsg('symmetry','return_no_irrep','The input point group identifier is out of range.',2,1)
         endif

         if (pg .eq. C1_id) return_no_irrep = 1
         if (pg .eq. Cs_id) return_no_irrep = 2
         if (pg .eq. C2_id) return_no_irrep = 2
         if (pg .eq. Ci_id) return_no_irrep = 2
         if (pg .eq. C2v_id) return_no_irrep = 4
         if (pg .eq. C2h_id) return_no_irrep = 4
         if (pg .eq. D2_id) return_no_irrep = 4
         if (pg .eq. D2h_id) return_no_irrep = 8

   end function return_no_irrep

   subroutine get_geometry_obj(this,geometry)
      implicit none
      class(symmetry_obj) :: this
      type(geometry_obj) :: geometry

      integer :: err

         if (.not.(this%initialized)) call xermsg('symmetry','get_geometry_obj','The object has not been initialized.',1,1)

         geometry%no_sym_op = this%no_sym_op
         geometry%sym_op = this%sym_op
         geometry%no_nuc = this%no_nuc
         geometry%use_symmetry = this%use_symmetry
         allocate(geometry%nucleus,source=this%nucleus,stat=err)
         if (err .ne. 0) call xermsg('symmetry','get_geometry_obj','Memory allocation error.',err,1)

   end subroutine get_geometry_obj

   !> \param[in] l L angular number of the real spherical harmonic
   !> \param[in] m M angular number of the real spherical harmonic
   !> \param[in] pg Integer identifier of the point group for which we want to obtain the point-group symmetry of the real spherical harmonic \f$X_{L,M}\f$.
   !> The module const contains the (Abelian) point-group identifiers and the table of characters for each point group as taken from Atkins, Tables for Group Theory, 1990.
   !> We assume that the spherical harmonics are defined using the usual choice of the spherical polar coordinates: i.e. the positive unit vector Z in the z-direction is obtained as a vector product of 
   !> the unit vectors X and Y pointing along the positive x and y axes: Z = X x Y; Z=cos(theta), X=sin(theta)cos(phi), Y=sin(theta)sin(phi).
   !> \param[out] sph_harm_pg_sym Integer corresponding to the irreducible representation of the point group into which the given spherical harmonic \f$X_{L,M}\f$ belongs. The identification is based on 
   !> determining the transformation properties of the spherical harmonic with respect to the symmetry elements of the point group and matching these with a row of the character table.
   !> \param[out] nam Character of length 2 labeling the irreducible representation.
   !> \todo Make sure that the coefficients c(:) in the cartesian->spherical transform are used correctly: accumulate first the contributions to each combination of exponents x,y,z and then test if the
   !> coefficient is nonzero. (For larger l the coefficients sometimes compensate to produce a zero resulting coefficient)!
   function get_sph_harm_pg_sym(this,l,m,nam)
   use const
   use special_functions, only: cfp_sph_to_cart_mapping
      implicit none
      class(symmetry_obj) :: this
      integer :: get_sph_harm_pg_sym
      integer, intent(in) :: l, m
      character(len=sym_op_nam_len), intent(out) :: nam

      integer :: s_x, s_y, s_z !signs of the spherical harmonic L,M after application of the symmetry operations corresponding to reflection in the planes x, y, z.
      integer :: c2_x, c2_y, c2_z !signs of the spherical harmonic L,M after application of the symmetry operations corresponding to rotations by pi around the axes x, y, z
      integer :: c_i !sign after application of the inversion operation
      integer :: s_h !sign after application of the reflection in the plane of the molecule
      integer :: c2, s_v, s_vp, s_3, s_2, s_1, c2_2, c2_1

      integer :: i, no_terms, chars(8)
      integer, allocatable :: x_exp(:),y_exp(:),z_exp(:)
      real(kind=cfp), allocatable :: c(:)
      logical :: is_x, is_y, is_z

         if (.not.(this%initialized)) call xermsg('symmetry','get_sph_harm_pg_sym','The object has not been initialized.',1,1)

         if (abs(m) > l) then
            print *,l,m
            call xermsg('symmetry','get_sph_harm_pg_sym','abs(m) > l',2,1)
         endif

         get_sph_harm_pg_sym = 0

         if (l > 0) then
            !set the transformation properties of the spherical harmonics under the symmetry operations
            call cfp_sph_to_cart_mapping(l,m,c,x_exp,y_exp,z_exp)
         endif

         s_x = 1
         s_y = 1
         s_z = 1

         !set the transformation properties of the spherical harmonics under the symmetry operations
         no_terms = merge(size(c), 0, allocated(c))
         do i=1,no_terms
            if (c(i) .ne. 0.0_cfp) then !if an odd exponent is present for a given coordinate then we know that this spherical harmonic will change sign under reflection
               if (mod(x_exp(i),2) .ne. 0) s_x = -1
               if (mod(y_exp(i),2) .ne. 0) s_y = -1
               if (mod(z_exp(i),2) .ne. 0) s_z = -1
            endif
         enddo

         !rotations can be expressed using reflections
         c2_z = s_x*s_y
         c2_y = s_x*s_z
         c2_x = s_z*s_y

         !inversion is also simple
         c_i = s_x*s_y*s_z

         !compare the transformation properties with definitions of the irreducible representations of the given point group symmetry.
         select case (this%pg)
         case (C1_id)
            chars(1) = 1
            get_sph_harm_pg_sym = match_irr(chars(1:1),C1_char_tab)
            nam = C1_names(get_sph_harm_pg_sym)

         case (Ci_id)
            chars(1) = c_i
            get_sph_harm_pg_sym = match_irr(chars(1:1),Ci_char_tab)
            nam = Ci_names(get_sph_harm_pg_sym)

         case (Cs_id)

            !define the sigma_h operation; the canonical one is Z
            if (trim(this%sym_op(1)) .eq. 'X') s_h = s_x
            if (trim(this%sym_op(1)) .eq. 'Y') s_h = s_y
            if (trim(this%sym_op(1)) .eq. 'Z') s_h = s_z

            chars(1) = s_h
            get_sph_harm_pg_sym = match_irr(chars(1:1),Cs_char_tab)
            nam = Cs_names(get_sph_harm_pg_sym)

         case (C2h_id)

            !the canonical type of the symmetry elements is: XY, Z

            !define the sigma_h operation
            if (trim(this%sym_op(1)) .eq. 'X' .or. trim(this%sym_op(2)) .eq. 'X') s_h = s_x
            if (trim(this%sym_op(1)) .eq. 'Y' .or. trim(this%sym_op(2)) .eq. 'Y') s_h = s_y
            if (trim(this%sym_op(1)) .eq. 'Z' .or. trim(this%sym_op(2)) .eq. 'Z') s_h = s_z

            !define the C2 operation
            if (trim(this%sym_op(1)) .eq. 'XY' .or. trim(this%sym_op(2)) .eq. 'XY') c2 = c2_z
            if (trim(this%sym_op(1)) .eq. 'YZ' .or. trim(this%sym_op(2)) .eq. 'YZ') c2 = c2_x
            if (trim(this%sym_op(1)) .eq. 'XZ' .or. trim(this%sym_op(2)) .eq. 'XZ') c2 = c2_y

            !use the inversion, C2 and sigma_h operations to determine the irreducible representation
            chars(1:3) = (/c2,c_i,s_h/)
            get_sph_harm_pg_sym = match_irr(chars(1:3),C2h_char_tab)
            nam = C2h_names(get_sph_harm_pg_sym)

         case (C2v_id)
            !the order and type of the symmetry elements matters; the canonical order and type is: X, Y
            is_x = .false.
            is_y = .false.
            is_z = .false.
           
            if (trim(this%sym_op(1)) .eq. 'X' .or. trim(this%sym_op(2)) .eq. 'X') is_x = .true.
            if (trim(this%sym_op(1)) .eq. 'Y' .or. trim(this%sym_op(2)) .eq. 'Y') is_y = .true.
            if (trim(this%sym_op(1)) .eq. 'Z' .or. trim(this%sym_op(2)) .eq. 'Z') is_z = .true.

            !determine the C2 axis
            if (is_x .and. is_y) c2 = c2_z
            if (is_y .and. is_z) c2 = c2_x
            if (is_x .and. is_z) c2 = c2_y

            !the first symmetry element defines the sigma_v' operation, e.g. for H2O this would be the plane of the molecule
            if (trim(this%sym_op(1)) .eq. 'X') s_vp = s_x
            if (trim(this%sym_op(1)) .eq. 'Y') s_vp = s_y
            if (trim(this%sym_op(1)) .eq. 'Z') s_vp = s_z

            !the second symmetry element defines the sigma_v operation
            if (trim(this%sym_op(2)) .eq. 'X') s_v = s_x
            if (trim(this%sym_op(2)) .eq. 'Y') s_v = s_y
            if (trim(this%sym_op(2)) .eq. 'Z') s_v = s_z

            chars(1:2) = (/s_v,s_vp/)
            get_sph_harm_pg_sym = match_irr(chars(1:2),C2v_char_tab)
            nam = C2v_names(get_sph_harm_pg_sym)

         case (D2h_id)
            !the order of the X,Y,Z symmetry elements defines the orientation of the molecule with respect to the X,Y,Z planes of reflection; the canonical order is: X,Y,Z

            if (trim(this%sym_op(1)) .eq. 'X') s_1 = s_x !for pyrazine this would define the plane of the molecule
            if (trim(this%sym_op(1)) .eq. 'Y') s_1 = s_y
            if (trim(this%sym_op(1)) .eq. 'Z') s_1 = s_z

            if (trim(this%sym_op(2)) .eq. 'X') s_2 = s_x
            if (trim(this%sym_op(2)) .eq. 'Y') s_2 = s_y
            if (trim(this%sym_op(2)) .eq. 'Z') s_2 = s_z

            if (trim(this%sym_op(3)) .eq. 'X') s_3 = s_x
            if (trim(this%sym_op(3)) .eq. 'Y') s_3 = s_y
            if (trim(this%sym_op(3)) .eq. 'Z') s_3 = s_z

            chars(1:3) = (/s_3,s_2,s_1/)
            get_sph_harm_pg_sym = match_irr(chars(1:3),D2h_char_tab)
            nam = D2h_names(get_sph_harm_pg_sym)

         case (D2_id)
            !the order and type of the rotational symmetry elements matters; the canonical order and type is: XZ, YZ

            if (trim(this%sym_op(1)) .eq. 'XZ') c2_2 = c2_y
            if (trim(this%sym_op(1)) .eq. 'XY') c2_2 = c2_z
            if (trim(this%sym_op(1)) .eq. 'YZ') c2_2 = c2_x

            if (trim(this%sym_op(2)) .eq. 'XZ') c2_1 = c2_y
            if (trim(this%sym_op(2)) .eq. 'XY') c2_1 = c2_z
            if (trim(this%sym_op(2)) .eq. 'YZ') c2_1 = c2_x

            chars(1:2) = (/c2_2,c2_1/)
            get_sph_harm_pg_sym = match_irr(chars(1:2),D2_char_tab)
            nam = D2_names(get_sph_harm_pg_sym)

         case (C2_id)
            is_x = .false.
            is_y = .false.
            is_z = .false.

            if (trim(this%sym_op(1)) .eq. 'XY') is_z = .true.
            if (trim(this%sym_op(1)) .eq. 'YZ') is_x = .true.
            if (trim(this%sym_op(1)) .eq. 'XZ') is_y = .true.

            !determine the C2 axis
            if (is_x) c2 = c2_x
            if (is_y) c2 = c2_y
            if (is_z) c2 = c2_z

            chars(1:1) = (/c2/)
            get_sph_harm_pg_sym = match_irr(chars(1:1),C2_char_tab)
            nam = C2_names(get_sph_harm_pg_sym)

         case default
            print *,this%pg,this%sym_op(:)
            call xermsg('symmetry','get_sph_harm_pg_sym','Unknown pg symmetry identifier.',3,1)
         end select

   end function get_sph_harm_pg_sym

   function match_irr(chars,char_tab)
      implicit none
      integer, intent(in) :: chars(:), char_tab(:,:)
      integer :: match_irr
      
      integer :: irr, sym_el, n_irr, n_sym_els
      logical :: match

         n_irr = size(char_tab,1)
         n_sym_els = size(char_tab,2)

         if (n_sym_els .ne. size(chars)) then
            call xermsg('symmetry', 'match_irr', 'Mismatch between size of character table and vector of characters.', 1, 1)
         end if
         match_irr = 0

         do irr = 1,n_irr
            match = .true.
            do sym_el = 1,n_sym_els
               if (chars(sym_el) .ne. char_tab(irr,sym_el)) match = .false.
            enddo
            if (match) then
               match_irr = irr
               return
            endif
         enddo

         call xermsg('symmetry','match_irr','There is no IRR that matches the input vector of characters',2,1)

   end function match_irr

   function get_name_of_irr(this,irr)
      use const
      implicit none
      class(symmetry_obj) :: this
      character(len=sym_op_nam_len) :: get_name_of_irr
      integer, intent(in) :: irr

         if (.not. this % initialized) then
            call xermsg ('symmetry', 'get_name_of_irr', 'The object has not been initialized.', 1, 1)
         end if
         if (irr <= 0 .or. irr > this % get_no_irrep(this % pg)) then
            call xermsg ('symmetry', 'get_name_of_irr', 'On input irr is out of range.', 2, 1)
         end if

         get_name_of_irr = this%irr_names(irr)

   end function get_name_of_irr

end module symmetry
