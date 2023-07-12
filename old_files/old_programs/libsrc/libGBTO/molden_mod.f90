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
!> \ingroup Molden
!> This module contains the object molden_input_obj which controls operations on a Molden input file. It reads, on demand, the basis set, nuclear data and molecular orbitals into the appropriate data
!> structures. The same data structures can be also written out to a file. The Molden file can hold basis set data for GTOs and STOs. Currently, reading of the GTO basis and the molecular orbitals 
!> has been implemented. The Molden format allows for use of flags which specify that the molecular orbitals are given in terms of coefficients for the spherical GTOs. 
!> These flags are: [5D],[5D10F],[7F],[5D7F],[5D],[9G] and currently ARE NOT supported by the molden_input_obj object. If any of these flags are found on the Molden file an error message is issued.
!> However, implementing the use of these flags is trivial and in fact simpler than having to work with orbital coefficients corresponding to the contracted cartesian GTOs which is the default.
!> When reading the MO coefficients we read-in the cartesian coefficients first and then convert them to the coefficients for the corresponding spherical GTOs. (We do the reverse when writing out orbital
!> coefficients for the spherical GTOs). This is done, for each shell, calculating the overlap integrals between the normalized contracted cartesian GTO basis functions and the normalized contracted 
!> spherical GTO corresponding to the same shell. The Molden file must be (according to its format specification) a formatted file. The exact format of all numbers and data is taken from Molpro.
!> All format specifications are contained in the module molden_const. There should be no constants and format definitions in this module - everything must be defined in the module molden_const. This 
!> ensures that if the input/output format needs to be changed that this can be done consistently for input/output only in one place.
module molden_mod
   use precisn
   use const, only: stdout, line_len, fmat
   use common_obj
   use molden_const
   use basis_data_generic_mod

   private

   public molden_input_obj

   !> \class <molden_input_obj>
   !> This object contains the data and methods necessary to control input and output of basis set and MO data in the Molden format.
   !> The user starts interaction with this object by calling the init method. After that the type-bound methods (e.g. read, write, final) become available.
   !> \todo Only master process should read/write everything. At the moment it must be ensured externally that each process writes into a file with a different name otherwise there will be conflicts.
   type molden_input_obj
      !> Path to the Molden file.
      character(len=line_len), private :: input_file
      !> Number of lines on the Molden file.
      integer, private :: n_lines = 0
      !> Line numbers for the array molden_file where the individual sections of the Molden file start.
      integer, private :: atoms_line = 0, gto_line = 0, mo_line = 0
      !> All lines of the Molden file: used only if the file has been initialized for input.
      character(len=line_len), allocatable, private :: molden_file(:)
      !> CGTO shells read-in from the Molden file.
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
      !> Molecular orbitals read-in from the Molden file.
      type(orbital_data_obj), allocatable :: orbital_data(:)
      !> List of atoms of size no_atoms. Set by the read_atoms procedure.
      type(nucleus_type), allocatable, private :: atom(:)
      !> This variable is set to 1 when the object has been initialized, i.e. the Molden file has been opened and checked.
      logical, private :: initialized = .false.
      !> If the Molden file specified in input_file is to be used for input then this variable will be set to 1. If the file is used for output then this variable will be set to 2.
      integer, private :: io = 0
      !> Unit number for the Molden file. It is set upon opening the input/output Molden file by the procedure 'init'.
      integer, private :: mldu = 0
      !> This variable is set to .true. if the Molden file contains molecular geometry data.
      logical, private :: contains_geometry = .false.
      !> This variable is set to .true. if the geometry data is in atomic units. If it is in Angstrom units then it remains 0.
      logical, private :: geom_units_au = .false.
      !> This variable is set to .true. if the Molden file contains GTO basis set data.
      logical, private :: contains_basis = .false.
      !> Point group ID of the molecule (if known; 0 otherwise).
      integer, private :: pg = 0
      !> Process or ignore alpha and beta orbitals: 0 = use all orbitals, 1 = alpha only (default), 2 = beta only.
      integer, private :: alpha_or_beta = 1
      !> If the Molden file contains GTO basis set data, then this variable is set (by get_molden_info) to the total number of shells of contracted spherical GTOs in the basis set.
      integer, private :: no_cgto_shells = 0
      !> If the Molden file contains GTO basis set data, then this variable is set (by get_molden_info) to the total number of contracted spherical GTOs in the basis set.
      integer, private :: no_cgto = 0
      !> If the Molden file contains GTO basis set data, then this variable is set (by get_molden_info) to the total number of contracted cartesian GTOs in the basis set.
      integer, private :: no_cgto_cart = 0
      !> If the Molden file contains GTO basis set data, then this variable is set (by get_molden_info) to the total number of primitive spherical GTOs in the basis set.
      integer, private :: no_pgto = 0
      !> This variable is set to .true. if the Molden file contains molecular orbitals.
      logical, private :: contains_mo = .false.
      !> This is set to true if the molecular orbital data contains the Sym flag specifying the orbital number and symmetry.
      logical, private :: contains_sym = .false.
      !> If the Molden file contains molecular orbitals then this number contains the value of the largest irreducible representation of an orbital on the file, i.e. the Molden file does not neccessarily
      !> include all orbitals for the system so the value max_irr cannot be generally used to determine the actual number of IRRs corresponding to the point group.
      integer, private :: max_irr = 0
      !> If the Molden file contains molecular orbitals then this number contains the number of orbitals in each of the irreducible representations.
      integer, private :: no_orbs(1:8) = 0
      !> Number of atoms. Set by the read_atoms procedure.
      integer, private :: no_atoms = 0
      !> How many molecular orbitals are on the Molden file.
      integer, private :: tot_orbs = 0
      !> Stores the number of molecular orbitals that have been read-in so far.
      integer, private :: read_orbs = 0
      !> If the Molden file is open for output then this variable is set to true once the geometry data have been written to the file using the method write_geometry.
      logical, private :: written_geom = .false.
      !> If the Molden file is open for output then this variable is set to true once the basis set data have been written to the file using the method write_basis.
      logical, private :: written_basis = .false.
      !> If the Molden file is open for output then this variable is set to true once the first molecular orbital has been written to the file using the method write_orbitals. This flag allows me to decide
      !> whether it is neccessary to write the header for the MO section to the Molden file or whether this has been already done.
      logical, private :: written_first_mo = .false.
      !> spherical_coeffs(l) == .true. if the MO coefficients corresponding to CGTOs with angular momentum l are for spherical CGTOs as opposed to cartesian CGTOs (which is the default). 
      !> The spherical flags are: [5D], [5D10F], [5D7F], [7F], [9G].
      logical, private :: spherical_coeffs(0:max_molden_l) = .false.
   contains
      !> Opens the input/output Molden file given in the first argument. This method requires as an argument an integer value specifying whether the Molden file will be used for input or output of data.
      !> For input set io=1, for output set io=2.
      procedure :: init => init_molden_file
      !> Transfers all data read-in from the Molden file from the type-bound structures atom,CGTO_shell_data,orbital_data into the equivalent allocatable structures supplied by the user.
      procedure :: read
      !> Writes all input data geometry,CGTO basis and molecular orbitals to the Molden file.
      procedure :: write
      !> Prints out the list of energy-sorted orbitals as read-in from the Molden file.
      procedure :: print_energy_sorted_orbital_table
      !> Returns in one integer array the number of molecular orbitals of each IRR saved on the file.
      procedure :: get_mo_num => get_mo_number_sym
      !> Finalization procedure that correctly closes the input/output Molden file and resets the object's variables.
      procedure :: final => close_molden_file
      !> Returns the number of CGTOs on the Molden file. Available following a successful call to init and in the reading mode.
      procedure :: get_no_cgto => get_number_of_cgtos
      !> Returns the number of CGTOs on the Molden file. Available following a successful call to init and in the reading mode.
      procedure :: get_no_cgto_shells => get_number_of_cgto_shells
      !> Private procedure used to determine what type of data a Molden file contains. It determines most of the logical variables of the molden_input_obj. Currently we check for presence of only some of
      !> the data Molden file can hold. Refer to the Molden website for full list of data the Molden format can hold: http://www.cmbi.kun.nl/molden/molden_format.html
      procedure, private :: get_molden_info
      !> Private procedure that reads in the geometry information into the private data variables: no_atoms, atom(:).
      procedure, private :: read_atoms
      !> Reads geometry data from the Molden file into a supplied array of type nucleus. It also gives the number of atoms of the molecule.
      procedure, private :: read_geometry => read_molden_geometry
      !> This procedure reads the GTO basis set from the Molden file into the user-supplied basis data structure of the type basis_data. Note that the supplied basis set can be also the mixed one.
      !> Keeping the input structure general we can initialize the cgto_basis_data part of the mixed basis set too.
      procedure, private :: read_basis => read_gto_basis
      !> Writes geometry data into the Molden file from a supplied array of type nucleus. We assume that the coordinates of the atoms on input are in atomic units. These are converted to Angstroms when
      !> writing the Molden file. If the scattering centre is present (i.e. atom with nuc == 0) then this atom must be the last one in the list of atoms. This is due to the need to
      !> output into Molden file only positive indices of the nuclei. On output to the Molden file the scattering centre is assigned index of the last atom +1: that's why it has to be
      !> the last one.
      procedure, private :: write_geometry => write_molden_geometry
      !> Write GTO basis set data in the Molden format. If the scattering centre is present (i.e. atom with nuc == 0) then this atom must be the last one in the list of atoms. 
      !> This is due to the need to output into Molden file only positive indices of the nuclei. On output to the Molden file the scattering centre is assigned index of the 
      !> last atom +1: that's why it has to be the last one.
      procedure, private :: write_basis => write_gto_basis
      !> Reads-in the selected molecular orbitals from a Molden file into the orbital_data_obj provided by the user. On input a contracted GTO basis set is required corresponding to the orbital being read-in. 
      !> The assumption is that the order of the contracted GTOs in the basis set corresponds to the order of the corresponding MO coefficients on the Molden file. 
      !> This assumption is always satisfied as long as the GTO basis set has been obtained from the same Molden file using the routine read_basis. 
      !> Otherwise you have to know what you're doing. However, the orbital_obj allows for calculation of total charge in the orbital using the density matrix allowing the user to check whether the orbitals
      !> have been read-in as expected. A further requirement on the input basis set is that the contracted spherical GTOs are sorted in such a way that for each shell with angular momentum L their order 
      !> is: M= -L, -L+1,...,L-1,L.
      !> If the user attempts to read more molecular orbitals then there are on the file an error message is issued and the program is terminated.
      procedure, private :: read_orbitals => read_selected_mo
      !> Writes into the Molden file data for a set of molecular orbitals. The input is the orbital data (orbital_data_obj) and the corresponding atomic basis set of contracted spherical GTOs.
      !> The orbital coefficients are converted to the Molden coefficients for contracted cartesian GTO basis and then written out.
      procedure, private :: write_orbitals => write_all_orbitals
   end type molden_input_obj

contains

   subroutine init_molden_file (this, inp, io, alpha_or_beta, pg)
      use utils, only: xermsg
      implicit none
      class(molden_input_obj) :: this
      integer, intent(in) :: io
      integer, intent(in), optional :: alpha_or_beta, pg
      character(len=line_len), intent(in) :: inp

      integer :: err
      logical :: ex, op
      character(len=len('UNKNOWN')) :: rw, frm

         write(stdout,'("--------->","molden_input_obj:init")')

         if (io > 2 .or. io <= 0) then
            call xermsg ('molden_mod', 'init_molden_file', &
                         'The input argument io has an incorrect value. The allowed values are: 1, 2.', 0, 1)
         end if

         !if the object has been initialized then it means that the corresponding Molden file is open -> close it now and reset the default values.
         if (this%initialized) then
            close(unit=this%mldu)
            this % initialized       = .false.
            this % mldu              = 0
            this % io                = 0
            this % n_lines           = 0
            this % atoms_line        = 0
            this % gto_line          = 0
            this % mo_line           = 0
            this % contains_geometry = .false.
            this % contains_basis    = .false.
            this % contains_mo       = .false.
            this % contains_sym      = .false.
            this % geom_units_au     = .false.
            this % no_atoms          = 0
            this % read_orbs         = 0
            this % tot_orbs          = 0
            this % max_irr           = 0
            this % no_orbs           = 0
            this % no_cgto           = 0
            this % no_cgto_cart      = 0
            this % no_cgto_shells    = 0
            this % no_pgto           = 0
            this % written_geom      = .false.
            this % written_basis     = .false.
            this % written_first_mo  = .false.
            if (allocated(this%molden_file)) deallocate(this%molden_file)
            if (allocated(this%CGTO_shell_data)) deallocate(this%CGTO_shell_data)
            if (allocated(this%orbital_data)) deallocate(this%orbital_data)
            if (allocated(this%atom)) deallocate(this%atom)
        endif

        ! store/reset the point group ID
        if (present(pg)) then
           this % pg = pg
        else
           this % pg = 0  ! unknown
        end if

        ! store/reset the alpha-beta orbital selection flag
        if (present(alpha_or_beta)) then
           this % alpha_or_beta = alpha_or_beta
        else
           this % alpha_or_beta = 1  ! alpha only
        end if

        inquire(file=inp,iostat=err,readwrite=rw,exist=ex,opened=op,formatted=frm)

        if (err /= 0) call xermsg('molden_mod','init_molden_file','Error executing the inquire statement.',1,1)
        if (op) call xermsg('molden_mod','init_molden_file','The input file is already opened by another application.',2,1)
        if (frm .eq. 'NO') call xermsg('molden_mod','init_molden_file','The input file cannot be opened for formatted access.',3,1)
        if (ex .and. rw .eq. 'NO') call xermsg('molden_mod','init_molden_file','The input file cannot be used for read/write.',4,1)
        if (io == 2 .and. ex) call xermsg('molden_mod','init_molden_file','The output file already exists.',5,1)
        if (io == 1 .and. .not.(ex)) call xermsg('molden_mod','init_molden_file','The input file does not exist.',6,1)

        this%io = io        

        if (io == 1) then !open the Molden file for formatted read
           open(file=inp,newunit=this%mldu,status='old',form=fmat,iostat=err)
           if (err /= 0) call xermsg('molden_mod','init_molden_file','Error opening the input file.',err,1)
        endif

        if (io == 2) then !open the Molden file for formatted write
           open(file=inp,newunit=this%mldu,status='new',form=fmat,iostat=err)
           if (err /= 0) call xermsg('molden_mod','init_molden_file','Error opening the output file.',err,1)
        endif

        this%input_file = inp

        this%initialized = .true.

        write(stdout,'("molden_input_obj has been initialized with the file: ",a)') trim(adjustl(this%input_file))
   
        if (io .eq. 1) then
           !now load all lines from the file into this%molden_file array and check what data does the Molden file contain - geometry, basis set, molecular orbitals
           call this%get_molden_info

           !read-in all available information
           if (this%contains_geometry) call this%read_atoms
           if (this%contains_basis) call this%read_basis(this%CGTO_shell_data)
           if (this%contains_mo) then
              call this%read_orbitals(this%CGTO_shell_data,this%no_orbs,this%orbital_data)
              !call this%print_energy_sorted_orbital_table
           endif
        endif

        write(stdout,'("<---------","done:molden_input_obj:init")')

   end subroutine init_molden_file

   subroutine close_molden_file(this)
      use utils, only: xermsg
      implicit none
      class(molden_input_obj) :: this

      integer :: err
      logical :: op

         write(stdout,'("--------->","molden_input_obj:final")')

         if (.not.(this%initialized)) then
            call xermsg('molden_mod','close_molden_file','The object has not been initialized.',1,1)
         else !close the Molden file and set the default values of the variables for this object.
            inquire(unit=this%mldu,opened=op)
            if (op) then
               close(unit=this%mldu,iostat=err)
               if (err /= 0) call xermsg('molden_mod','close_molden_file','Error occured closing the Molden file.',2,1)
            endif
            this % initialized       = .false.
            this % mldu              = 0
            this % io                = 0
            this % n_lines           = 0
            this % atoms_line        = 0
            this % gto_line          = 0
            this % mo_line           = 0
            this % contains_geometry = .false.
            this % contains_basis    = .false.
            this % contains_mo       = .false.
            this % contains_sym      = .false.
            this % geom_units_au     = .false.
            this % pg                = 0
            this % no_atoms          = 0
            this % read_orbs         = 0
            this % tot_orbs          = 0
            this % max_irr           = 0
            this % no_orbs           = 0
            this % no_cgto           = 0
            this % no_cgto_cart      = 0
            this % no_cgto_shells    = 0
            this % no_pgto           = 0
            this % no_atoms          = 0
            this % written_geom      = .false.
            this % written_basis     = .false.
            this % written_first_mo  = .false.
            this % spherical_coeffs(:) = .false.
            if (allocated(this%molden_file)) deallocate(this%molden_file)
            if (allocated(this%CGTO_shell_data)) deallocate(this%CGTO_shell_data)
            if (allocated(this%orbital_data)) deallocate(this%orbital_data)
            if (allocated(this%atom)) deallocate(this%atom)
            write(stdout,'(/,"The Molden file ",a," has been closed.")') trim(adjustl(this%input_file))
            this%input_file = ''
         endif

         write(stdout,'("<---------","done:molden_input_obj:final")')

   end subroutine close_molden_file

   subroutine get_molden_info(this)
      use const, only: pg_irr_names
      use utils, only: xermsg, search_string
      implicit none
      class(molden_input_obj) :: this

      integer :: err, pos_angs, pos_au, l, no_contr, num, sym, i, dot
      character(len=line_len) :: one_line, sym_flag, sym_val
      character(len=1) :: ang_type
      logical :: found
      real :: one

         write(stdout,'("--------->","molden_input_obj:get_molden_info")')

         !Count the number of lines of the molden file:
         this%n_lines = 0
         do
            read(this%mldu,'(a)',end=100) one_line
            this%n_lines = this%n_lines + 1
         enddo
     100 rewind this%mldu

         write(stdout,'("Number of lines of the Molden file: ",i15)') this%n_lines
 
         !Slurp the whole file line-by-line into the array this%molden_file:
         allocate(this%molden_file(this%n_lines),stat=err)
         if (err /= 0) call xermsg('molden_mod','get_molden_info','Memory allocation failed.',err,1)
 
         found = .false.
         this%atoms_line = 0
         this%gto_line = 0
         this%mo_line = 0
         do i=1,this%n_lines
            read(this%mldu,'(a)') one_line
            this%molden_file(i) = one_line
            if (index(one_line,header_molden) > 0) found = .true.
            if (index(one_line,header_atoms)  > 0) this % atoms_line = i
            if (index(one_line,header_gto)    > 0) this % gto_line   = i
            if (index(one_line,header_mo)     > 0) this % mo_line    = i
         enddo

         write(stdout,'("The Molden file has been read into the memory.")')

         if (.not. found) then
            call xermsg ('molden_mod', 'get_molden_info', &
                         'The input file does not contain the Molden file header.', 1, 1)
         end if
         if (this % atoms_line >= this % gto_line) then
            call xermsg ('molden_mod', 'get_molden_info', &
                         'Input file error: the header [ATOMS] must preceed the header [GTO].', 2, 1)
         end if
         if (this % gto_line >= this % mo_line) then
            call xermsg ('molden_mod', 'get_molden_info', &
                         'Input file error: the header [GTO] must preceed the header [MO].', 3, 1)
         end if

         this%contains_geometry = .false.
         this%contains_basis = .false.
         this%contains_mo = .false.

         if (this%atoms_line > 0) this%contains_geometry = .true.
         if (this%gto_line > 0) this%contains_basis = .true.
         if (this%mo_line > 0) this%contains_mo = .true.

         if (this%contains_geometry) then
            write(stdout,'("The Molden file contains the header for geometry data.")')

            !Examine the header line to see whether it contains the Angs or AU string so that we know which units the geometry is in.
            one_line = this%molden_file(this%atoms_line)
            pos_angs = index(one_line,str_angs) !pos_angs >0 if the 'Angs' substring appears on the line
            pos_au = index(one_line,str_au) !pos_au >0 if the 'AU' substring appears on the line
            if (pos_au > 0 .and. pos_angs > 0) then
                call xermsg ('molden_mod', 'get_molden_info', &
                             'Both length units appear on the line with the header [Atoms]. I am confused.', 4, 1)
            end if
            if (pos_au == 0 .and. pos_angs == 0) then
                call xermsg ('molden_mod', 'get_molden_info', &
                             'Units of length are not specified on the line with the header [Atoms].', 5, 1)
            end if
            if (pos_au > 0 .and. pos_angs .eq. 0) then !the units of length are atomic units
               this%geom_units_au = .true.
               write(stdout,'("The units of length are atomic units.")')
            endif
            if (pos_au .eq. 0 .and. pos_angs > 0) then !the units of length are Angstroms
               this%geom_units_au = .false.
               write(stdout,'("The units of length are Angstroms.")')
            endif

         end if

         if (this%contains_basis) then
            write(stdout,'("The Molden file contains the header for GTO basis data.")')

            if (.not. this % contains_geometry) then
                call xermsg ('molden_mod', 'get_molden_info', &
                             'Input file error: the header [ATOMS] is missing but is needed to read-in the basis set data.', 6, 1)
            end if

            !check for flags whose use we have not implemented yet.
            do i=this%gto_line,this%mo_line-1
               one_line = this%molden_file(i)
               if (index(one_line, str_sp) > 0) then
                  call xermsg ('molden_mod', 'get_molden_info', &
                               'The use of the sp flag has not been implemented yet.', 7, 1) !sp = contraction of s and p primitive GTO functions into one
               end if
            enddo !i

            !Determine the number of contracted spherical GTOs in the basis set
            this%no_cgto_shells = 0
            this%no_cgto = 0 !total number of contracted spherical GTOs in the basis set
            this%no_cgto_cart = 0 !total number of contracted cartesian GTOs in the basis set
            this%no_pgto = 0 !total number of primitive spherical GTOs in the basis set
            do i=this%gto_line+1,this%mo_line-1
               one_line = this%molden_file(i)

               !find out if the line contains GTO info
               do l=0,max_molden_l !over all supported GTO angular types
                  if (index(one_line,gto_typ(l)) .ne. 0) then !the line contains one of the letters 's,p,d,f,g' and hence info on the next shell of GTOs
                     this%no_cgto_shells = this%no_cgto_shells + 1
                     this%no_cgto = this%no_cgto + 2*l+1 !increment by the number of contracted spherical GTOs corresponding to this shell
                     this%no_cgto_cart = this%no_cgto_cart + (l+1)*(l+2)/2 !increment by the number of contracted cartesian GTOs corresponding to this shell

                     read(one_line,*,iostat=err) ang_type, no_contr, one !contracted GTO angular type, no. of primitives, 1.00
                     if (err /= 0) then
                        call xermsg ('molden_mod', 'get_molden_info', &
                                     'Problem parsing the string containing the contracted GTO info.', 7, 1)
                     end if
                     this%no_pgto = this%no_pgto + no_contr*(2*l+1) !increment by the number of primitive spherical GTOs corresponding to this shell
                     exit
                  endif
               enddo !l
            enddo !i

            write(stdout,'("The GTO basis set contains ",i5," contracted and ",i5," corresponding primitive spherical GTOs.")') &
                this%no_cgto, this%no_pgto
            write(stdout,'("The corresponding number of contracted cartesian GTOs: ",i5)') this%no_cgto_cart
            write(stdout,'("Number of shells of contracted GTOs in the basis: ",i5)') this%no_cgto_shells
         else
            this%no_cgto = 0
            this%no_cgto_cart = 0
            this%no_cgto_shells = 0
            this%no_pgto = 0
         end if

         if (this%contains_mo) then
            write(stdout,'("The Molden file contains the header for Molecular orbitals data.")')

            this%tot_orbs = 0

            !count the number of orbitals on the file and check for unsupported flags:
            this%no_orbs(:) = 0
            this%max_irr = 0
            do i=this%gto_line+1,this%n_lines
               one_line = this%molden_file(i)

               if (index(one_line,str_ene) > 0) this%tot_orbs = this%tot_orbs + 1

               !Extract the symmetry information (if available)
               if (index(one_line,str_sym) > 0) then
                  this%contains_sym = .true.
                  !line containing "Sym="
                  read(one_line,*) sym_flag, sym_val
                  sym = 0
                  dot = index(sym_val, '.')
                  if (dot /= 0) then
                     read(sym_val(1:dot-1), *) num  ! number before dot
                     read(sym_val(dot+1:),  *) sym  ! number after dot
                  else if (this % pg >= 1) then
                    !sym = findloc(pg_irr_names(:, this % pg), sym_val, 1)                ! clean way (Fortran 2008)
                     sym = maxloc(merge(1, 0, pg_irr_names(:, this % pg) == sym_val), 1)  ! hack for compatibility
                     num = this % no_orbs(sym) + 1
                  end if

                  if (sym == 0) then
                     call xermsg('molden_mod', 'get_molden_info', 'Incomprehensible symmetry notation "' // sym_val // '".', 1, 1)
                  end if

                  if (sym > this%max_irr) this%max_irr = sym
                  this%no_orbs(sym) = this%no_orbs(sym) + 1
               endif

               !Look for the flags which denote spherical (as opposed to cartesian) MO coefficients for particular angular momentum are saved on the file 
               if (index(one_line,str_5d10f) > 0) this%spherical_coeffs(2) = .true.
               if (index(one_line,str_5d7f) > 0) then
                  this%spherical_coeffs(0) = .true.
                  this%spherical_coeffs(1) = .true.
                  this%spherical_coeffs(2) = .true.
               end if
               if (index(one_line,str_7f) > 0) this%spherical_coeffs(3) = .true.
               if (index(one_line,str_5d) > 0 .or. index(one_line,str_5d7f) > 0) then
                  this%spherical_coeffs(2) = .true.
                  this%spherical_coeffs(3) = .true.
               endif
               if (index(one_line,str_9g) > 0) then
                  this%spherical_coeffs(4) = .true.
               endif

            enddo !i
            write(stdout,'("Number of molecular orbitals on the file: ",i5)') this%tot_orbs
            if (this%tot_orbs > 0) then
               if (this%contains_sym) then
                  write(stdout,'("The symmetry information for the orbitals is available.")')
               else
                  write(stdout,'("The symmetry information for the orbitals is NOT available.")')
               endif
            endif

         end if

         write(stdout,'("<---------","done:molden_input_obj:get_molden_info")')

   end subroutine get_molden_info

   subroutine read(this,nucleus,CGTO_shell_data,orbital_data)
      implicit none
      class(molden_input_obj) :: this
      type(nucleus_type), allocatable :: nucleus(:)
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
      type(orbital_data_obj), allocatable :: orbital_data(:)

      integer :: err

         write(stdout,'("--------->","molden_input_obj:read")')

         if (.not.(this%initialized)) call xermsg('molden_mod','read','The Molden object has not been initialized.',1,1)

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'read', &
                         'The Molden object has not been initialized for input (io .eq. 1).', 2, 1)
         end if

         if (this%contains_geometry) then
            if (allocated(nucleus)) deallocate(nucleus)
            allocate(nucleus,source=this%atom,stat=err)
            if (err /= 0) call xermsg('molden_mod','read','Transfer of geometry data has failed.',3,1)
            write(stdout,'("Molecular geometry has been transfered.")')
         endif

         if (this%contains_basis) then
            if (allocated(CGTO_shell_data)) deallocate(CGTO_shell_data)
            allocate(CGTO_shell_data,source=this%CGTO_shell_data,stat=err)
            if (err /= 0) call xermsg('molden_mod','read','Transfer of basis data has failed.',4,1)
            write(stdout,'("CGTO basis data have been transfered.")')
         endif

         if (this%contains_mo) then
            if (allocated(orbital_data)) deallocate(orbital_data)
            allocate(orbital_data,source=this%orbital_data,stat=err)
            if (err /= 0) call xermsg('molden_mod','read','Transfer of molecular orbitals has failed.',5,1)
            write(stdout,'("Molecular orbital data have been transfered.")')
         endif

         write(stdout,'("<---------","molden_input_obj:read")')

   end subroutine read

   subroutine write(this,nucleus,CGTO_shell_data,orbital_data)
      implicit none
      class(molden_input_obj) :: this
      type(nucleus_type), allocatable :: nucleus(:)
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
      type(orbital_data_obj), allocatable :: orbital_data(:)

         write(stdout,'("--------->","molden_input_obj:write")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'write', 'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 2) then
            call xermsg ('molden_mod', 'write', 'The Molden object has not been initialized for output (io .eq. 2).', 2, 1)
         end if

         if (allocated(nucleus)) then
            call this%write_geometry(nucleus)
         endif

         if (allocated(CGTO_shell_data)) then
            if (.not. allocated(nucleus)) then
                call xermsg ('molden_mod', 'write', &
                             'The CGTO basis data cannot be written since the nucleus array on input has not been allocated.', 3, 1)
            end if
            call this%write_basis(CGTO_shell_data)
         endif

         if (allocated(orbital_data)) then
            if (.not. allocated(CGTO_shell_data)) then
                call xermsg ('molden_mod', 'write', &
                             'The orbital basis data cannot be written since the nucleus array on input &
                             &has not been allocated.', 4, 1)
            end if
            call this%write_orbitals(CGTO_shell_data,orbital_data)
         endif

         write(stdout,'("<---------","molden_input_obj:write")')

   end subroutine write

   subroutine read_atoms(this)
      use utils, only: xermsg, search_string
      use molden_const, only: cms_thrs
      implicit none
      class(molden_input_obj) :: this

      integer :: err, i, k, Z
      character(len=line_len) :: line
      character(len=2) :: element_name
      real(kind=cfp) :: cms(3), total_mass, d

         write(stdout,'("--------->","molden_input_obj:read_atoms")')

         if (.not.(this%initialized)) call xermsg('molden_mod','read_atoms','The Molden object has not been initialized.',1,1)

         if (.not. this % contains_geometry) then
            call xermsg ('molden_mod', 'read_atoms', 'The Molden file does not contain molecular geometry info.', 2, 1)
         end if

         !find out how many atoms there are in total
         this%no_atoms = this%gto_line-1 - this%atoms_line
         write(stdout,'(/,"Number of atoms: ",i5)') this%no_atoms

         !allocate memory for the atom data
         if (allocated(this%atom)) deallocate(this%atom)
         allocate(this%atom(1:this%no_atoms),stat=err)
         if (err /= 0) call xermsg('molden_mod','read_atoms','Memory allocation failed.',4,1)

         !now read-in the atoms one-by-one
         do i=1,this%no_atoms
            line = this%molden_file(this%atoms_line+i)
            read(line,*) element_name, this%atom(i)%nuc, Z, (this%atom(i)%center(k),k=1,3) !transfer the data to the the atom(:) array.

            if (.not.(this%geom_units_au)) then !convert the distances from Angstroms to a.u. if needed
               this%atom(i)%center(:) = this%atom(i)%center(:)*angstrom_to_au
            endif

            this%atom(i)%name = element_name
            !According to molden format what we use here as Z is the element's atomic number
            this%atom(i)%charge = Z  !we have to do it like this since this%atom(i)%charge is of type real.
            err = this%atom(i)%check() !check that the atom data are plausible.
            !call this%atom(i)%print !print the nuclear data
         enddo

         call calculate_cms(this%no_atoms,this%atom,cms,total_mass)

         d = sqrt(dot_product(cms,cms))

         if (d > cms_thrs) then

            call xermsg ('molden_mod', 'read_molden_geometry', &
                         'Molecule is not centered on the center of mass or atomic masses are not as in the program (Molpro?) &
                         &used to generate the molden file.', 3, 0)

            write(stdout,'("Original center of mass coordinates are: ",3e25.15)') cms(1:3)
            write(stdout,'("Total mass of the molecule is (a.u.): ",e25.15)') total_mass

            !The loop below performs recentering to CMS but we probably don't want to enforce it so I've commented it out.
            !do i=1,this%no_atoms
            !   this%atom(i)%center(1:3) = this%atom(i)%center(1:3) - cms(1:3)
            !enddo

         endif

         write(stdout,'("<---------","molden_input_obj:read_atoms")')

   end subroutine read_atoms

   subroutine read_molden_geometry(this,nucleus,n)
      use utils, only: xermsg
      implicit none
      class(molden_input_obj) :: this
      type(nucleus_type), allocatable, intent(out) :: nucleus(:)
      integer, intent(out) :: n

      integer :: i, err

         write(stdout,'("--------->","molden_input_obj:read_geometry")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'read_molden_geometry', &
                         'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'read_molden_geometry', &
                         'Attempt to read from a Molden file which has been associated for output only.', 2, 1)
         end if

         if (allocated(nucleus)) deallocate(nucleus)
         allocate(nucleus(1:this%no_atoms),stat=err)
         if (err /= 0) call xermsg('molden_mod','read_molden_geometry','Memory allocation failed.',3,1)

         !transfer the atom data
         do i=1,this%no_atoms
            nucleus(i) = this%atom(i)
         enddo

         n = this%no_atoms

         write(stdout,'("Molecular geometry transferred.")')

         write(stdout,'("<---------","done:molden_input_obj:read_geometry")')

   end subroutine read_molden_geometry

   subroutine write_molden_geometry(this,nucleus)
      use utils, only: xermsg
      implicit none
      class(molden_input_obj) :: this
      class(nucleus_type), intent(in) :: nucleus(:)

      integer :: i, j, k, err, no_atoms, Z, nuc

         write(stdout,'("--------->","molden_input_obj:write_geometry")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'write_molden_geometry', 'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 2) then
            call xermsg ('molden_mod', 'write_molden_geometry', &
                         'Attempt to write into a Molden file which has been associated for input only.', 2, 1)
         end if
 
         if (this % written_geom) then
            call xermsg ('molden_mod', 'write_molden_geometry', &
                         'The geometry data have been written into the Molden file already.', 3, 1)
         end if

         no_atoms = size(nucleus)

         !Check for duplicities of the %nuc attribute of the nuclei, i.e. make sure each nucleus has a distinct index
         do i=1,no_atoms
            do j=1,no_atoms
               if (j .ne. i) then
                  if (nucleus(j)%nuc .eq. nucleus(i)%nuc) then
                     call xermsg ('molden_mod', 'write_molden_geometry', &
                                  'Duplicity in the %nuc attribute for at least one pair of nuclei.', 4, 1)
                  endif
               endif
            enddo !j
         enddo !i

         !allocate space for the local atoms(:) structure that will hold the nuclear data as they were exactly written to the file.
         if (allocated(this%atom)) deallocate(this%atom)

         allocate(this%atom(no_atoms),stat=err)
         if (err /= 0) call xermsg('molden_mod','write_molden_geometry','Memory allocation error.',5,1)

         this%no_atoms = no_atoms

         !Write the Molden format header
         write(this%mldu,'(a)') header_molden

         !Write the geometry header
         write(this%mldu,'(/,a,1x,a)') header_atoms, str_angs

         do j=1,no_atoms
            !write the nuclear data one nucleus at a time and save the data into the local atoms(:) structure for reference.

            !convert the nuclear charge value from float to integer
            Z = nucleus(j)%charge
            !if the nuclear charge is fractional then inform the user that it has been converted to integer due to Molden format requirements
            if (Z .ne. nucleus(j)%charge) then
               call xermsg ('molden_mod', 'write_molden_geometry', &
                            'The nuclear charge is fractional but had to be converted to an integer value &
                            &due to Molden format requirements.', 6, 2)
            endif

            !write the line containing information on one atom
            nuc = nucleus(j)%nuc
            if (nuc .eq. 0) then
               if (j .ne. no_atoms) then
                  call xermsg ('molden_mod', 'write_molden_geometry', &
                               'If the scattering centre is present then it must be the last atom in the list.', 7, 1)
               else
                  !The nuclear index on the Molden file must be larger than 0. (the same is done in write_gto_basis)
                  if (j .eq. 1) then
                     nuc = 1
                  else
                     nuc = max(1,this%atom(j-1)%nuc+1)
                  endif
               endif
            endif

            write(this%mldu, form_atom) nucleus(j)%name, nuc, Z, (nucleus(j)%center(k)/angstrom_to_au, k=1,3)

            !save the nuclear data as they were exactly written to the file.
            this%atom(j) = nucleus(j)
            this%atom(j)%charge = Z

         enddo !j

         !set the flag that is telling me whether the geometry has been written to the output or not.
         this%written_geom = .true.

         write(stdout,'("Molecular geometry written to the Molden file.")')

         write(stdout,'("<---------","done:molden_input_obj:read_geometry")')

   end subroutine write_molden_geometry

   function get_number_of_cgtos(this)
      implicit none
      class(molden_input_obj) :: this
      integer :: get_number_of_cgtos

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'get_number_of_cgtos', 'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'get_number_of_cgtos', &
                         'The Molden object has not been initialized for input: no_cgto value is not available.', 2, 1)
         end if

         if (.not. this % contains_basis) then
            call xermsg ('molden_mod', 'get_number_of_cgtos', &
                         'The Molden does not contain basis set data: no_cgto value is not available.', 3, 1)
         end if

         get_number_of_cgtos = this%no_cgto

   end function get_number_of_cgtos

   function get_number_of_cgto_shells(this)
      implicit none
      class(molden_input_obj) :: this
      integer :: get_number_of_cgto_shells

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'get_number_of_cgto_shells', 'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'get_number_of_cgto_shells', &
                         'The Molden object has not been initialized for input: no_cgto_shells value is not available.', 2, 1)
         end if

         if (.not. this % contains_basis) then
            call xermsg ('molden_mod', 'get_number_of_cgto_shells', &
                         'The Molden does not contain basis set data: no_cgto_shells value is not available.', 3, 1)
         end if

         get_number_of_cgto_shells = this%no_cgto_shells

   end function get_number_of_cgto_shells

   subroutine read_gto_basis(this,CGTO_shell_data)
      use utils, only: xermsg, search_string
      use basis_data_generic_mod
      implicit none
      class(molden_input_obj) :: this
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)

      real(kind=cfp) :: one
      integer :: err, empty_line_cnt, no_contr, i, k, l, at_i, shell_index, n
      character(len=line_len) :: line
      character(len=1) :: ang_type

         write(stdout,'("--------->","molden_input_obj:read_basis")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'read_gto_basis', 'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'read_gto_basis', &
                         'Attempt to read from a Molden file which has been associated for output only.', -1, 1)
         end if
   
         if (.not.(this%contains_basis)) then
            call xermsg('molden_mod','read_gto_basis','The Molden file does not contain GTO basis info.',2,1)
         endif
   
         if (.not.(this%contains_geometry)) then
            call xermsg ('molden_mod', 'read_gto_basis', &
                         'The Molden file does not contain molecular geometry info which is &
                         &needed in order to read the GTO basis.', 3, 1)
         endif
         
         !print the nuclear data
         write(stdout,'(/"Nuclear data from the Molden file ",a," (coordinates in a.u.): ")') trim(adjustl(this%input_file))
         write(stdout,'("Index of the nucleus; Z; Element name; Coordinates of the nucleus")')
         do i=1,this%no_atoms
            write(stdout,'(5X,i2,1X,f5.1,1X,a2,1X,3f25.15)') &
                this % atom(i) % nuc, &
                this % atom(i) % charge, &
                this % atom(i) % name, &
                (this % atom(i) % center(k), k = 1,3)
         enddo

         write(stdout,'(/,"Reading GTO basis set from the Molden file...")')

         !Allocate space for the CGTO shells that will be read-in
         if (allocated(CGTO_shell_data)) deallocate(CGTO_shell_data)
         allocate(CGTO_shell_data(this%no_cgto_shells),stat=err)
         if (err /= 0) call xermsg('molden_mod','read_gto_basis','Memory allocation failed.',err,1)

         !reading of the GTO basis set starts
         empty_line_cnt = 0 !this counts how many empty lines we have encountered since the end of the data for the last GTO.
         shell_index = 0 !this counts the shells of contracted GTOs
         n = this%gto_line
         do !this is a do loop over the atoms
            n = n + 1
            line = this%molden_file(n) !the whole line containing the atom sequence number and 0

            !now check whether this line contains GTO data
            if (index(line,'[') .ne. 0 .or. index(line,']') .ne. 0) exit   !finish: we have reached the end of the GTO info since another card begins

            if (len_trim(line) .eq. 0) then !the line is empty
               empty_line_cnt = empty_line_cnt + 1 !if the line is empty then increment the counter
               if (empty_line_cnt .eq. 2) exit !finish if we have encountered a second empty line

            else !the line is not empty, i.e. a GTO data for another atom follows. In this section we read-in data for all GTOs on this atom

               read(line,*) at_i, i !get the atom sequence number at_i

               empty_line_cnt = 0
               do !continue reading until we reach the end of the GTO data for this atom, i.e. until we reach an empty line
                  n = n +1
                  line = this%molden_file(n)

                  if (len_trim(line) .eq. 0) then
                     empty_line_cnt = empty_line_cnt + 1
                     exit !terminate the do loop if the line is empty
                  else !the line is not empty and contains the contraction info: in this section we read-in all primitive GTOs corresponding to one contracted GTO
                     read(line,*,iostat=err) ang_type, no_contr, one !contracted GTO angular type, no. of primitives, 1.00
                     if (err /= 0) then
                        call xermsg ('molden_mod', 'read_gto_basis', &
                                     'Problem parsing the string containing the contracted GTO info.', 5, 1)
                     end if

                     !determine which angular type we have:
                     k = max_molden_l+1
                     do !loop over all Molden-supported GTO angular types
                        k = k - 1
                        if (k .eq. -1) call xermsg('molden_mod','read_gto_basis','An unsupported GTO angular type detected.',6,1) !i.e. we have underrun the lower bound of gto_typ
                        err = index(line,gto_typ(k)) !find out whether the angular type string gto_typ(k) is present on the line
                        if (err > 0) exit !we've found which angular type we are reading-in now
                     enddo
                     l = k !L of the GTO
      
                     !Allocate space for no_contr primitives.
                     shell_index = shell_index + 1
                     call CGTO_shell_data(shell_index)%make_space(no_contr)
      
                     !now read-in the primitives: their exponents and contractions
                     do i=1,no_contr
                        n = n + 1
                        line = this%molden_file(n)

                        read(line,*) CGTO_shell_data(shell_index)%exponents(i), CGTO_shell_data(shell_index)%contractions(i) !exponent, contraction coeff
                        if (CGTO_shell_data(shell_index)%exponents(i) .le. 0.0_cfp) then
                           call xermsg('molden_mod','read_gto_basis','Error in Molden file: exponent .le. 0.',7,1)
                        endif
                        if (CGTO_shell_data(shell_index)%contractions(i) .eq. 0.0_cfp) then
                           call xermsg('molden_mod','read_gto_basis','Error in Molden file: contraction coefficient is zero.',8,0)
                        endif
                     enddo

                     !now we have all data required to construct the shell of contracted GTOs and add it to the atomic basis:
                     CGTO_shell_data(shell_index)%l = l !GTO L
                     CGTO_shell_data(shell_index)%number_of_functions = 2*CGTO_shell_data(shell_index)%l+1
                     CGTO_shell_data(shell_index)%center(1:3) = this%atom(at_i)%center(1:3) !atom coordinates
                     CGTO_shell_data(shell_index)%non_zero_at_boundary = .false.
                     call CGTO_shell_data(shell_index)%normalize
                     !CGTO_shell_data(shell_index)%nuc = at_i !number of the atom on which this GTO is sitting

                  endif !if the line is empty

               enddo !over all GTOs on this atoms

            endif !if the line is empty

         enddo !over all atoms

         write(stdout,'("finished")')

         write(stdout,'("<---------","molden_input_obj:read_basis")')

   end subroutine read_gto_basis

   subroutine write_gto_basis(this,CGTO_shell_data)
      use utils, only: xermsg
      implicit none
      class(molden_input_obj) :: this
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)

      integer :: i, j, k, shell, a, prims, i_zero, err, nuc, no_cgto_shells, max_l
      logical, allocatable :: shell_written(:)
      character(len=1) :: ang_type
      real :: one
      logical :: match

         write(stdout,'("--------->","molden_input_obj:write_basis")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'write_gto_basis', &
                         'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 2) then
            call xermsg ('molden_mod', 'write_gto_basis', &
                         'Attempt to write into a Molden file which has been associated for input only.', 2, 1)
         end if

         if (.not. this % written_geom) then
            call xermsg ('molden_mod', 'write_gto_basis', &
                         'Attempt to write the basis set data before the geometry data.', 3, 1)
         end if

         if (this % written_basis) then
            call xermsg ('molden_mod', 'write_gto_basis', &
                         'The basis set data have been written into the file already.', 4, 1)
         end if

         !Check that the centers of the GTOs are compatible with the centers of the nuclei as they were written to the file.
         max_l = 0
         no_cgto_shells = size(CGTO_shell_data)
         do i=1,no_cgto_shells
  
            !Find the largest angular momentum in the basis along the way
            max_l = max(max_l,CGTO_shell_data(i)%l)

            if (CGTO_shell_data(i)%l > max_molden_l) then
               call xermsg ('molden_mod', 'write_gto_basis', &
                            'The GTO basis set contains functions with L larger than what the MOLDEN format can handle. &
                            &GTO basis functions with L > max_molden_l will be ignored.', 6, 0)
            endif

            !Determine on which nucleus the current shell of GTOs is sitting
            k = 0
            j = size(this%atom)
            do 
               k = k + 1
               if (k > j) then
                  call xermsg ('molden_mod', 'write_gto_basis', &
                               'Nucleus matching the center of the GTO shell could not be found.', 7, 1)
               end if
               match = centers_match(CGTO_shell_data(i)%center,this%atom(k)%center)
               if (match) exit !we've found the matching nucleus
            enddo !k

         enddo !i

         !Write the GTO basis header
         write(this%mldu,'(a)') header_gto

         !Write the GTO basis functions in the order given by the MOLDEN format which is given by: 1) order of the atoms, 2) shell angular momentum.
         allocate(shell_written(no_cgto_shells),stat=err)
         if (err /= 0) call xermsg('molden_mod','write_gto_basis','Memory allocation error.',err,1)

         !As we write the GTO data we mark the shells which have been written out.
         shell_written = .false.

         !1) loop over all atoms sequentially
         do k=1,this%no_atoms

            nuc = this%atom(k)%nuc
            if (nuc .eq. 0) then
               if (k .ne. this%no_atoms) then
                  call xermsg ('molden_mod', 'write_gto_basis', &
                               'If the scattering centre is present then it must be the last atom in the list.', 9, 1)
               else
                  !The nuclear index on the Molden file must be larger than 0. (the same is done in write_geometry)
                  if (k .eq. 1) then
                     nuc = 1
                  else
                     nuc = max(1,this%atom(k-1)%nuc+1)
                  endif
               endif
            endif

            i_zero = 0
            write(this%mldu, form_atom_id) nuc, i_zero !write out the atom sequence number
      
            !2) loop over all shells of functions that we have in the basis: s,p,d,...,max_l
            do shell=0,max_l

               if (shell > max_molden_l) then
                  write(stdout,'("Ignoring GTO basis functions with L > ",i0)') shell
                  exit
               endif
      
               !character string corresponding to the current angular momentum: s,p,d,...,
               ang_type = gto_typ(shell)
      
               !Find out how many contracted GTOs with L=shell sit on the atom k.
               j = 0
               do i=1,no_cgto_shells
                  !Make sure that we do not confuse the continuum shells with
                  !shells corresponding to an atom sitting on CMS
                  if (CGTO_shell_data(i)%is_continuum() .and. k < this%no_atoms) cycle

                  match = centers_match(CGTO_shell_data(i)%center,this%atom(k)%center)
                  if (CGTO_shell_data(i)%l .eq. shell .and. match) then
                     if (CGTO_shell_data(i)%l .eq. shell) j = j +1
                  endif
               enddo !i
      
               !Loop over all contracted GTOs with the angular momentum L=shell
               do a=1,j
      
                  !Find the GTO shells that sit on the nucleus k and have the desired angular momentum and have not been written out yet.
                  do i=1,no_cgto_shells
                     !Make sure that we do not confuse the continuum shells with
                     !shells corresponding to an atom sitting on CMS
                     if (CGTO_shell_data(i)%is_continuum() .and. k < this%no_atoms) cycle

                     match = centers_match(CGTO_shell_data(i)%center,this%atom(k)%center)
                     if (match .and. CGTO_shell_data(i)%l .eq. shell .and. .not.(shell_written(i))) then
      
                        !Write out the header for this shell
                        one = 1.000000
                        write(this%mldu, form_gto_head) ang_type, CGTO_shell_data(i)%number_of_primitives, one !contracted GTO angular type, no. of primitives, 1.00
      
                        !Write out the primitive GTO data
                        do prims=1,CGTO_shell_data(i)%number_of_primitives
      
                           !Multiply in the overall GTO norms with the contraction coefficients to get the contractions in the MOLDEN convention.
                           write(this%mldu, form_gto_prim) CGTO_shell_data(i) % exponents(prims), &
                                                           CGTO_shell_data(i) % norm * CGTO_shell_data(i) % contractions(prims) !exponent, contraction coeff
      
                        enddo !prims
      
                        !mark this shell as written out so that we don't write it out again.
                        shell_written(i) = .true.
                        
                     endif
                  enddo !i
       
               enddo !a
            
            enddo !shell
      
            write(this%mldu,'("")') !empty line after each atom GTO data

         enddo !k

         !Finally, make sure that we have written out all shells of GTO function. If this test fails then it is either caused by some incompleteness of the nuclear data or a programming error in the above. 
         !Note that we perform this check only in case we should be writing out all GTO basis functions (that is in the case the basis doesn't contain functions with L>max_molden_l).
         if (max_l .le. max_molden_l .and. count(shell_written) .ne. no_cgto_shells) then
            print *,count(shell_written),no_cgto_shells
            call xermsg ('molden_mod', 'write_gto_basis', &
                         'The number GTO shells written out does not match the number of GTO shells in the basis.', 8, 1)
         endif

         this%written_basis = .true.

         write(stdout,'("GTO basis set written to the Molden file.")')

         write(stdout,'("<---------","molden_input_obj:write_basis")')

   end subroutine write_gto_basis

   subroutine read_selected_mo(this,CGTO_shell_data,nob,orbital_data)
      use utils, only: xermsg, search_string
      use basis_data_generic_mod
      use gto_routines, only: cart_cf_sph_cf,sph_cf_cart_cf 
      use const, only: thrs_lin_dep_gs_ortho, pg_irr_names
      implicit none
      class(molden_input_obj) :: this
      type(CGTO_shell_data_obj), intent(in) :: CGTO_shell_data(:)
      type(orbital_data_obj), allocatable :: orbital_data(:)
      integer, intent(in) :: nob(8)

      integer :: i, j, err, no_bf, no_cart_gtos, cnt, last_ind, shell, l, m, n_preceeding, sym_ind, irr, number_of_shells, n, dot
      integer :: no_orbs(8),no_sph_gtos
      character(len=line_len) :: line, energy_flag, occ_flag, spin_flag, sym_flag, sym_val
      character(len=max(len(str_alpha),len(str_beta))) :: spin
      integer, allocatable :: mo_index(:), i_exp(:), j_exp(:), k_exp(:)
      real(kind=cfp), allocatable :: cart_cf(:), sph_cf(:), temp_coefficients(:,:)
      logical :: found

         write(stdout,'("--------->","molden_input_obj:read_mo")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'read_selected_mo', &
                         'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'read_selected_mo', &
                         'Attempt to read from a Molden file which has been associated for output only.', -1, 1)
         end if

         if (.not. this % contains_mo) then
            call xermsg ('molden_mod', 'read_selected_mo', &
                         'The Molden input file does not contain molecular orbitals.', 2, 1)
         end if
 
         do i=1,8
            if (nob(i) > this%no_orbs(i)) then
               print *,i,nob(i),this%no_orbs(i)
               call xermsg ('molden_mod', 'read_selected_mo', &
                            'The number of requested orbitals to read is larger than the number of orbitals available.', 4, 1)
            endif
         enddo !i

         !Get information on the basis set:
         number_of_shells = size(CGTO_shell_data)
         j = 0
         do i=1,number_of_shells
            j = j + CGTO_shell_data(i)%number_of_functions
         enddo
         !if the basis set was not read-in using the molden_mod then this%no_cgto would not have been set!
         if (this%contains_basis) then !get the number of basis functions from the value that we've read in from the Molden file
            no_bf = this%no_cgto !number of contracted spherical GTOs. This is also the maximum number of orbitals on the file.
            if (no_bf .ne. j) then 
               print *,no_bf,j
               call xermsg ('molden_mod', 'read_selected_mo', &
                            'The number of basis function in the input atomic orbital basis does not equal the number &
                            &of basis functions from the basis set on the Molden file.', 5, 1)
            endif
         else
            no_bf = j
         endif

         ! there may be at most twice as many orbitals as basis functions (for UHF)
         if (2*no_bf < this%tot_orbs) then

            if (2*this%no_cgto_cart < this%tot_orbs) then
               call xermsg ('molden_mod', 'read_selected_mo', &
                            'The number of basis functions is smaller than the number of orbitals stored on the Molden file. &
                            &Basis not compatible?', 6, 1)
            else
               write(stdout,'("I am suspecting cartesian GTO basis was used to generate the molecular orbitals.")')
            endif
         end if

         !Allocate space for orbitals in all symmetries:
         if (allocated(orbital_data)) deallocate(orbital_data)
         allocate(orbital_data(8),stat=err) !we allocate the structure to the maximum number of symmetries possible so that outside this routine orbital_data for other symmetries can be added easilly.
         if (err /= 0) call xermsg('molden_mod','read_selected_mo','Memory allocation 1 failed.',err,1)

         !Set elementary orbital data:
         do i=1,8
            orbital_data(i)%irr = i
            orbital_data(i)%point_group = -1  !The point group ID must be set externally
            orbital_data(i)%number_of_coefficients = no_bf
            orbital_data(i)%number_of_functions = nob(i)
            if (nob(i) > 0) then
               allocate(orbital_data(i) % coefficients(no_bf,nob(i)), &
                        orbital_data(i) % energy(nob(i)), &
                        orbital_data(i) % occup(nob(i)), &
                        orbital_data(i) % spin(nob(i)), stat = err)
               if (err /= 0) call xermsg('molden_mod','read_selected_mo','Memory allocation 2 failed.',err,1)
               orbital_data(i) % coefficients = 0.0_cfp
               orbital_data(i) % energy = 0.0_cfp
               orbital_data(i) % occup = 0.0_cfp
               orbital_data(i) % spin = 0
            endif
         enddo !i

         this%read_orbs = 0 !number of orbitals read-in

         !Go through the [MO] section of the file and pick out the orbitals we want to read-in.
         n = this%mo_line
         no_orbs(:) = 0  ! tracks number of orbitals per IRR read from the file
         do
            if (this%read_orbs .eq. sum(nob)) exit !All orbitals have been read-in
            
            n = n + 1
            line = this%molden_file(n)

            found = .false.
            if (index(line,str_sym) > 0) found = .true. !we don't accept orbitals with no symmetry information:
            if (.not. found) then
                call xermsg ('molden_mod', 'read_selected_mo', &
                             'Problem in the Molden file: an input line containing Sym= is not present.', 8, 1)
            end if
   
            !Extract the symmetry information:
            if (index(line,str_sym) > 0) then
                !line containing "Sym="
                read(line, *) sym_flag, sym_val
                irr = 0
                dot = index(sym_val, '.')
                if (dot /= 0) then
                   read(sym_val(1:dot-1), *) sym_ind  ! number before dot
                   read(sym_val(dot+1:),  *) irr      ! number after dot
                else if (this % pg >= 1) then
                  !irr     = findloc(pg_irr_names(:, this % pg), sym_val, 1)                ! clean way (Fortran 2008)
                   irr     = maxloc(merge(1, 0, pg_irr_names(:, this % pg) == sym_val), 1)  ! hack for compatibility
                   sym_ind = no_orbs(irr) + 1
                   no_orbs(irr) = sym_ind
                end if
                if (irr == 0) then
                   call xermsg('molden_mod', 'get_molden_info', 'Incomprehensible symmetry notation "' // sym_val // '".', 1, 1)
                end if
            else !the Molden file does not contain symmetry information
                print *,line,str_sym
                call xermsg('molden_mod','read_selected_mo','Error parsing the line with orbital symmetry information.',9,1)
            endif

            !Skip this orbital if it is not one of the ones we want to read-in:
            if (sym_ind > nob(irr)) cycle
            write(stdout,'("Accepted orbital number.symmetry: ",i4,".",i1)') sym_ind,irr

            !Next line is energy information:
            n = n + 1
            line = this%molden_file(n)

            found = .false.
            if (index(line,str_ene) > 0) found = .true. !we don't accept orbitals with no energy information:
            if (.not. found) then
                call xermsg ('molden_mod', 'read_selected_mo', &
                             'Problem in the Molden file: an input line containing Ene= is not present.', 10, 1)
            end if

            !extract the energy information; line containing "Ene="
            read(line, *) energy_flag, orbital_data(irr)%energy(sym_ind)
            !extract spin information; line containing "Spin="
            n = n + 1
            line = this%molden_file(n)
            read(line, *) spin_flag, spin
            if (trim(adjustl(spin)) .eq. str_alpha) orbital_data(irr)%spin(sym_ind) = 1
            if (trim(adjustl(spin)) .eq. str_beta) orbital_data(irr)%spin(sym_ind) = 2
            !extract the occupation number; line containing "Occup="
            n = n + 1
            line = this%molden_file(n)
            read(line, *) occ_flag, orbital_data(irr)%occup(sym_ind)

            !the following lines on the Molden file contain the orbital coefficients:
            !loop over all spherical contracted basis set functions
              !for each shell: read at once coefficients for the contracted cartesian functions for the whole shell
                !transform the cartesian coefficients into coefficients for the contracted spherical GTOs.
   
            write(stdout,'("Index of the contracted cartesian/spherical GTO; cartesian/spherical MO coefficient; &
                           &Index of the contracted spherical GTO; spherical MO coefficient")')
   
            !loop over all spherical contracted basis set functions
            cnt = 0 !total number of contracted cartesian GTOs in the basis; we use this to check whether we've read-in all orbital coefficients
            shell = 1 !index of the shell of the GTO functions
            m = -CGTO_shell_data(shell)%l !the M value of the GTO that we should be processing next.
            do i=1,no_bf
               !see if we have exhausted all functions from the current shell; if yes then set the M value of the first GTO in the shell; prepare the GTO index to be incremented
               if (m > CGTO_shell_data(shell)%l) then
                  shell = shell + 1
                  if (shell > number_of_shells) then
                    call xermsg ('molden_mod', 'read_selected_mo', &
                                 'Internal inconsistency in the GTO basis set data.', -8, 1)
                  end if
                  m = -CGTO_shell_data(shell)%l
               endif
   
               !now read and convert the orbital coefficients from the Molden file
               l = CGTO_shell_data(shell)%l !current CGTO shell angular momentum
               if (l > max_molden_l) then
                  call xermsg ('molden_mod', 'read_selected_mo', &
                               'A GTO is present in the basis whose L exceeds maximum L allowed by the Molden format.', 8, 1)
               end if
               if (m .eq. -l) then !read-in the cartesian coefficients if and only if we encounter a new shell. Here we use the assumption about the order of the GTOs with a different M within one shell.
                  !If the coefficients for this shell are for spherical CGTOs then we don't need to perform the conversion from cartesian->spherical coefficient
                  if (this%spherical_coeffs(l)) then
            !         allocate(mo_index(1:2*l+1),sph_cf(1:2*l+1),cart_cf(1:2*l+1),stat=err)
no_cart_gtos=(l+1)*(l+2)/2
no_sph_gtos=2*l+1
             allocate(mo_index(1:no_sph_gtos),stat=err)
             allocate(cart_cf(1:no_sph_gtos),stat=err)
             allocate(sph_cf(1:no_sph_gtos),stat=err)
             allocate(i_exp(1:no_cart_gtos),stat=err)
             allocate(j_exp(1:no_cart_gtos),stat=err)
             allocate(k_exp(1:no_cart_gtos),stat=err)
                     if (err /= 0) call xermsg('molden_mod','read_selected_mo','Memory allocation 3 failed.',err,1)
                     do j=1,2*l+1
                        n = n + 1
                        if (n > this % n_lines) then
                            call xermsg ('molden_mod', 'read_selected_mo', &
                                         'End of file reached too soon. Basis set does not correspond to this MO?', 9, 1)
                        end if
                        line = this%molden_file(n)

                        if (err < 0) then
                            call xermsg ('molden_mod', 'read_selected_mo', &
                                         'End of file reached too soon. Basis set does not correspond to this MO?', 9, 1)
                        end if
                        read(line,*,iostat=err) mo_index(j), cart_cf(j) !index of the contracted GTO and the corresponding spherical MO coefficient
!                        read(line,*,iostat=err) mo_index(j), sph_cf(j) !index of the contracted GTO and the corresponding spherical MO coefficient
                        if (err /= 0) then
                            call xermsg ('molden_mod', 'read_selected_mo', &
                                         'Error parsing the line containing coefficients of the molecular orbital. &
                                         &Reading past the last coefficient of the orbital?', 10, 1)
                        end if
                        if (cart_cf(j) .ne. 0.0_cfp .and. abs(cart_cf(j)) < 10e-10_cfp) then
                           write(stdout,'("Coefficient ",e25.15," has been neglected.")') cart_cf(j)
                           cart_cf(j) = 0.0_cfp
                        endif
                        if (sph_cf(j) .ne. 0.0_cfp .and. abs(sph_cf(j)) < 10e-10_cfp) then
                           write(stdout,'("Coefficient ",e25.15," has been neglected.")') sph_cf(j)
                           sph_cf(j) = 0.0_cfp
                        endif
                        write(stdout,'("spherical",e25.15)') cart_cf(j)
                     enddo
                     !The spherical coefficients are ordered in the following way: M = 0,+1,-1,+2,-2,...,+L,-L. This order has been obtained from Molden/gaussian.f source file.
                     !Hay que chequear este orden mencionado en el formato
                     !molden: dice que para D, F y G es: 0, +1 -1,+2,-2, etc. Asumiendo ese orden, aca el los ordena de -l,...,l
                     !en dalton para P: pz esta en el tercer lugar en los
                     !coeficientes del output (no genera molden en cartesianas),
                     !y el coeficiente coincide con el tercer lugar del
                     !esferico (m=0), de lo que se deduce que para este caso p el
                     !orden puede ser distinto. Esto se discute en https://github.com/MolSSI/QCSchema/issues/45
                     !donde se menciona que para el orbital p la convencion se
                     !ajusta a xyz para igualar a las cartesianas, en este caso de dalton creemos que es: +1,-1,0

                     If(l.eq.1)then
                     sph_cf(1)=cart_cf(2)
                     sph_cf(2)=cart_cf(3)
                     sph_cf(3)=cart_cf(1)
                     else
                     sph_cf(l+0+1) = cart_cf(1) !m=0
                     do j=1,l
                        n_preceeding = (j-1)*2+1
                        sph_cf(l+j+1) = cart_cf(n_preceeding+1) !+m
                        sph_cf(l-j+1) = cart_cf(n_preceeding+2) !-m
                     enddo
                     end if

Do j=1,cart_shell(l)
call process_ang_fact(ang_fact(l,j),i_exp(j),j_exp(j),k_exp(j))
End Do


                     call sph_cf_cart_cf (l, i_exp, j_exp, k_exp, &
                                          CGTO_shell_data(shell) % exponents, &
                                          CGTO_shell_data(shell) % contractions, &
                                          CGTO_shell_data(shell) % norm, &
                                          sph_cf, cart_cf)
                     last_ind = mo_index(no_cart_gtos)
                     cnt = cnt + 2*l+1 !increment the total number of coefficients we've read-in
                     last_ind = mo_index(2*l+1)
                  else !the MO coefficients are for cartesian GTOs
                     no_cart_gtos = cart_shell(l) !(l+1)*(l+2)/2 !how many cartesians correspond to the shell with angular momentum L.
                     !read-in coefficients of the cartesian GTOs and use them to construct the spherical GTO coefficients
                     allocate(mo_index(1:no_cart_gtos),cart_cf(1:no_cart_gtos),sph_cf(1:2*l+1),i_exp(1:no_cart_gtos),&
                              j_exp(1:no_cart_gtos),k_exp(1:no_cart_gtos),stat=err)
                     if (err /= 0) call xermsg('molden_mod','read_selected_mo','Memory allocation 4 failed.',err,1)
                     do j=1,no_cart_gtos !loop over the contracted cartesian GTOs
                        n = n + 1
                        if (n > this % n_lines) then
                            call xermsg ('molden_mod', 'read_selected_mo', &
                                         'End of file reached too soon. Basis set does not correspond to this MO?', 11, 1)
                        end if
                        line = this%molden_file(n)
                        read(line,*,iostat=err) mo_index(j), cart_cf(j) !index of the contracted cartesian GTO and the corresponding MO coefficient
                        !read(line,'(i4,1x,f15.11)',iostat=err) mo_index(j), cart_cf(j) !NOTE: '(i4,1x,f15.11)' has been replaced in Molpro by a flexible format, so we put * above.
                        if (err /= 0) then
                            call xermsg ('molden_mod', 'read_selected_mo', &
                                         'Error parsing the line containing coefficients of the molecular orbital. &
                                         &Reading past the last coefficient of the orbital?', 12, 1)
                        end if
                        call process_ang_fact(ang_fact(l,j),i_exp(j),j_exp(j),k_exp(j)) !store the exponents of the contracted cartesian GTOs corresponding to this shell
                        if (cart_cf(j) .ne. 0.0_cfp .and. abs(cart_cf(j)) < 10e-10_cfp) then
                           write(stdout,'("Coefficient ",e25.15," has been neglected.")') cart_cf(j)
                           cart_cf(j) = 0.0_cfp
                        endif
                     enddo
                     cnt = cnt + no_cart_gtos !increment the total number of cartesian coefficients we've read-in
                     !calculate the contracted spherical GTO coefficients (sph_cf) for all M using the (contracted cartesian GTO)x(contracted spherical GTO) overlaps
                     call cart_cf_sph_cf (l, i_exp, j_exp, k_exp, &
                                          CGTO_shell_data(shell) % exponents, &
                                          CGTO_shell_data(shell) % contractions, &
                                          CGTO_shell_data(shell) % norm, &
                                          cart_cf, sph_cf)
                     last_ind = mo_index(no_cart_gtos)
                  endif
                  !transfer the spherical MO coefficients to the orbital_data structure
!write(*,*) sph_cf
!write(*,*) "-----"
!write(*,*) cart_cf
!pause
                  do j=-l,l
                     !again, we assume here that the basis functions are ordered such that M=-L,-L+1,...,L-1,L (see the do loop over j above)                       
                     orbital_data(irr)%coefficients(i-1 + j+l+1,sym_ind) = sph_cf(j+l+1)
                  enddo !j
                  !for l .le. 1 the cartesian GTOs and the spherical GTOs are the same, so we check that we have indeed obtained almost identical (up to a rounding error) coefficients:
                  if (l .eq. 0) then
                     if (abs(sph_cf(l+1) - cart_cf(1)) > s_p_conv_tol * abs(cart_cf(1))) then
                        call xermsg ('molden_mod', 'read_selected_mo', 'Error in spherical -> cartesian conversion.', 13, 1)
                     end if
                  endif
                  if (l .eq. 1) then
                     if (abs(sph_cf(1+l+1) - cart_cf(1)) > s_p_conv_tol * abs(cart_cf(1))) then
                        call xermsg ('molden_mod', 'read_selected_mo', 'Error in spherical -> cartesian conversion.', 14, 1) !x
                     end if
                     if (abs(sph_cf(-1+l+1) - cart_cf(2)) > s_p_conv_tol * abs(cart_cf(2))) then
                        call xermsg ('molden_mod', 'read_selected_mo', 'Error in spherical -> cartesian conversion.', 15, 1) !y
                     end if
                     if (abs(sph_cf(0+l+1) - cart_cf(3)) > s_p_conv_tol * abs(cart_cf(3))) then
                        call xermsg ('molden_mod', 'read_selected_mo', 'Error in spherical -> cartesian conversion.', 16, 1) !z
                     end if
                  endif

                  if (this%spherical_coeffs(l)) then
                     do j=1,2*l+1
                        write(stdout,'(i5,e25.15,i5,e25.15)') mo_index(j), sph_cf(j), i-1+j, sph_cf(j)
                     enddo
                     deallocate(mo_index,sph_cf,cart_cf,i_exp,j_exp,k_exp)
                  else
                     do j=1,no_cart_gtos
                        if (j .le. 2*l+1) then
                           write(stdout,'(i5,e25.15,i5,e25.15)') mo_index(j), cart_cf(j), i-1+j, sph_cf(j)
                        else
                           write(stdout,'(i5,e25.15)') mo_index(j), cart_cf(j)
                        endif
                     enddo
                     deallocate(mo_index,cart_cf,sph_cf,i_exp,j_exp,k_exp)
                  endif 
   
               endif
   
               !write(stdout,'(i5,e)') i, orbital_data%cf(i)
   
               m = m + 1 !set the M value of the next GTO from the basis that should be processed; we use it to determine if the next shell starts.
   
            enddo !i
   
            !the index of the last cartesian coefficient must match the number of cartesian coefficients actually read-in:
            if (cnt /= last_ind) then
                call xermsg ('molden_mod', 'read_selected_mo', &
                             'The expected number of MO coefficients does not match the number of coefficients read-in.', 17, 1)
            end if
   
            !finally, increment the counter for the number of orbitals that we've read-in
            this%read_orbs = this%read_orbs + 1

         enddo !all orbitals on the Molden file

         ! filter read orbitals by the requested spin
         if (this % alpha_or_beta == 1 .or. this % alpha_or_beta == 2) then
            do i = 1, 8
               ! skip IRRs with no orbitals
               if (orbital_data(i) % number_of_functions == 0) cycle

               ! compress arrays by discarding unwanted spin-orbitals
               n = 0
               do j = 1, orbital_data(i) % number_of_functions
                  if (orbital_data(i) % spin(j) == this % alpha_or_beta) then
                     n = n + 1
                     orbital_data(i) % spin(n)   = orbital_data(i) % spin(j)
                     orbital_data(i) % energy(n) = orbital_data(i) % energy(j)
                     orbital_data(i) % occup(n)  = orbital_data(i) % occup(j)
                     orbital_data(i) % coefficients(:,n) = orbital_data(i) % coefficients(:,j)
                  end if
               end do
               orbital_data(i) % number_of_functions = n

               ! quench the array of orbital coefficients
               call move_alloc(orbital_data(i) % coefficients, temp_coefficients)

               allocate(orbital_data(i) % coefficients&
                        (orbital_data(i)%number_of_coefficients, orbital_data(i) % number_of_functions),stat=err)
               if (err /= 0) call xermsg('molden_mod','read_selected_mo','Memory allocation 5 failed.',err,1)

               orbital_data(i) % coefficients(:, 1:orbital_data(i) % number_of_functions) = &
                            temp_coefficients(:, 1:orbital_data(i) % number_of_functions)

               deallocate(temp_coefficients)
            end do
         end if

         write(stdout,'("<---------","done:molden_input_obj:read_mo")')

   end subroutine read_selected_mo

   subroutine get_mo_number_sym(this,n)
      use utils, only: xermsg
      implicit none
      class(molden_input_obj) :: this
      integer, intent(out) :: n(1:8)

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'get_mo_number_sym', 'The Molden object has not been initialized.', 1, 1)
         end if

         if (.not. this % contains_sym) then
            call xermsg ('molden_mod', 'get_mo_number_sym', 'The Molden file does not contain orbital symmetry data.', 2, 1)
         end if

         n(1:8) = this%no_orbs(1:8)

   end subroutine get_mo_number_sym

   subroutine write_all_orbitals(this,CGTO_shell_data,orbital_data)
      use utils, only: xermsg
      use atomic_basis_mod
      use gto_routines, only: sph_cf_cart_cf
      use basis_data_generic_mod
      implicit none
      class(molden_input_obj) :: this
      type(orbital_data_obj) :: orbital_data(:)
      type(CGTO_shell_data_obj) :: CGTO_shell_data(:)

      character(len=max(len(str_alpha),len(str_beta))) :: spin
      integer :: i, j, k, l, err, sph_ind, ind, no_cgto_shells, max_l, n, irr, sym_ind, number_of_functions
      integer, allocatable :: ix(:),iy(:),iz(:)
      real(kind=cfp), allocatable :: molden_cf(:)
      real(kind=cfp) :: sph_cgto_norm

         write(stdout,'("--------->","molden_input_obj:write_orbitals")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'write_all_orbitals', &
                         'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 2) then
            call xermsg ('molden_mod', 'write_all_orbitals', &
                         'Attempt to write into a Molden file which has been associated for input only.', 2, 1)
         end if

         if (.not. this % written_geom) then
            call xermsg ('molden_mod', 'write_all_orbitals', &
                         'Attempt to write the basis set data before the geometry data.', 3, 1)
         end if

         if (.not. this % written_basis) then
            call xermsg ('molden_mod', 'write_all_orbitals', &
                         'Attempt to write the orbital data before the basis set data.', 4, 1)
         end if

         no_cgto_shells = size(CGTO_shell_data)
 
         number_of_functions = 0
         do i=1,no_cgto_shells
            number_of_functions = number_of_functions + CGTO_shell_data(i)%number_of_functions
         enddo !i

         do i=1,size(orbital_data)
            err = orbital_data(i)%check()
            if (err /= 0) then
               print *,i
               call xermsg ('molden_mod', 'write_all_orbitals', &
                            'Checking of the input orbital_data failed. See orbital_data_obj%check for details.', err, 1)
            endif
            if (orbital_data(i)%number_of_coefficients .ne. number_of_functions) then
               call xermsg ('molden_mod', 'write_all_orbitals', &
                            'The CGTO_shell_data and orbital_data correspond to basis sets with a different &
                            &number of basis functions.', 6, 1)
            endif
         enddo !i

         !Find the largest angular momentum in the CGTO basis:
         max_l = 0
         do i=1,no_cgto_shells
            max_l = max(max_l,CGTO_shell_data(i)%l)
         enddo
         if (max_l > max_molden_l) then
            call xermsg ('molden_mod', 'write_all_orbitals', &
                         'The GTO basis set contains functions with L larger than what the MOLDEN format can handle.&
                         & MO coefficients corresponding to GTO basis functions with L > max_molden_l will be ignored.', 8, 0)
         end if

         !write (if neccessary) the header for the Molecular orbital section
         if (.not.(this%written_first_mo)) then
            write(this%mldu,'(/,a)') header_mo
         endif

         !Loop over the orbitals on input and save them one-by-one
         do irr=1,size(orbital_data)
            do sym_ind = 1,orbital_data(irr)%number_of_functions

               !write the symmetry information; line containing "Sym="
               if (orbital_data(irr)%point_group > 0) then
                  write(this%mldu, form_sym) sym_ind, irr
                  write(stdout,'("Orbital number.symmetry: ",i4,".",i1)') sym_ind, irr
               else
                  write(stdout,'("Orbital symmetry information is not available!")')
               endif
      
               !write the energy information; line containing "Ene="
               write(this%mldu, form_ene) orbital_data(irr)%energy(sym_ind)
               
               !write the spin information; line containing "Spin="
               spin = 'N/A'
               if (orbital_data(irr)%spin(sym_ind) .eq. 1) spin = str_alpha
               if (orbital_data(irr)%spin(sym_ind) .eq. 2) spin = str_beta
               write(this%mldu, form_spin) spin
      
               !write the occupation number; line containing "Occup="
               write(this%mldu, form_occ) orbital_data(irr)%occup(sym_ind)

               write(stdout,'("Index of the contracted cartesian GTO; cartesian MO coefficient; &
                              &Index of the contracted spherical GTO; spherical MO coefficient")')
               !
               !------- Contracted spherical -> contracted cartesian GTO transformation
               !        Space for the cartesian exponents and for the spherical coefficients for one shell
               !
               l = min(max_l,max_molden_l) !min is here since we must make sure that we don't overrun the bounds in cart_shell below which ends at max_molden_l.
               j = cart_shell(l)
               allocate(molden_cf(j),ix(j),iy(j),iz(j),stat=err)
               if (err /= 0) call xermsg('molden_mod','write_all_orbitals','Memory allocation error.',err,1)

               j = 0
               ind = 0
               sph_ind = 0  !starting index for the spherical orbital coefficients in each shell
               do i=1,no_cgto_shells !over all shells of spherical GTOs
      
                  sph_cgto_norm = CGTO_shell_data(i)%norm
      
                  !get the cartesian exponents for this shell in the Molden order
                  l = CGTO_shell_data(i)%l
                  if (l > max_molden_l) then
                     write(stdout,'("Skipping MO coefficients for GTO functions with L > ",i0)') max_molden_l
                     cycle
                  endif
                  do j=1,cart_shell(l)
                     call process_ang_fact(ang_fact(l,j),ix(j),iy(j),iz(j))
                  enddo !j
      
                  k = 2*l+1 !how many spherical angular behaviours there are in this shell
                 
                  !Convert the cartesian orbital coefficients for this shell into coefficients for the spherical GTOs. The result is in molden_cf
                  n = CGTO_shell_data(i)%number_of_primitives
                  call sph_cf_cart_cf (l, ix, iy, iz, &
                                       CGTO_shell_data(i) % exponents(1:n), &
                                       CGTO_shell_data(i) % contractions(1:n), &
                                       sph_cgto_norm, &
                                       orbital_data(irr) % coefficients(sph_ind+1:sph_ind+k,sym_ind), &
                                       molden_cf)
      
                  !print out the coefficients for this shell and increment the index of the last spherical coefficient treated.
                  do j=1,cart_shell(l)
                     ind = ind + 1
                     write(this%mldu,'(i5,e25.15)') ind, molden_cf(j) !index of the contracted cartesian GTO and the corresponding MO coefficient
      
                     if (j .le. 2*l+1) then
                        sph_ind = sph_ind + 1
                        write(stdout,'(i5,e25.15,i5,e25.15)') &
                            ind, molden_cf(j), sph_ind, orbital_data(irr)%coefficients(sph_ind,sym_ind)
                     else
                        write(stdout,'(i5,e25.15)') ind, molden_cf(j)
                     endif
      
                     !todo for l .le. 1 the cartesian GTOs and the spherical GTOs are the same, so we check that we have indeed obtained almost identical (up to a rounding error) coefficients
                     !(see read_selected_mo)
      
                  enddo !j
      
               enddo !i
      
               !update the flag if this was the first molecular orbital written out
               if (.not.(this%written_first_mo)) this%written_first_mo = .true.

               deallocate(molden_cf,ix,iy,iz)
 
            enddo !sym_ind

         enddo !irr

         write(stdout,'("<---------","molden_input_obj:write_orbitals")')

   end subroutine write_all_orbitals

   !> Returns .true./.false. depending on whether the coordinates of the two centers match or not.
   function centers_match(center_A,center_B)
      implicit none
      real(kind=cfp), intent(in) :: center_A(3), center_B(3)
      logical :: centers_match
     
         centers_match = .true.
         if (center_A(1) .ne. center_B(1)) centers_match = .false.
         if (center_A(2) .ne. center_B(2)) centers_match = .false.
         if (center_A(3) .ne. center_B(3)) centers_match = .false.
    
   end function centers_match

   subroutine calculate_cms(no_atoms,atom,cms,total_mass)
      use const, only: amass
      implicit none
      integer, intent(in) :: no_atoms
      type(nucleus_type) :: atom(no_atoms)
      real(kind=cfp), intent(out) :: cms(3), total_mass

      integer :: i
      real(kind=cfp) :: atom_mass

         cms = 0.0_cfp
         total_mass = 0.0_cfp
         do i=1,no_atoms
            atom_mass = amass(nint(atom(i)%charge))
            total_mass = total_mass + atom_mass
            cms(1:3) = cms(1:3) + atom(i)%center(1:3)*atom_mass
         enddo

         cms = cms/total_mass

   end subroutine calculate_cms

   subroutine print_energy_sorted_orbital_table(this)
      use common_obj, only: print_orbital_table
      implicit none
      class(molden_input_obj) :: this

      integer :: i, j, n, err, n_tgt
      integer, allocatable :: num_sym(:,:)
      real(kind=cfp), allocatable :: energies(:)

         write(stdout,'("--------->","molden_input_obj:print_energy_sorted_orbital_table")')

         if (.not. this % initialized) then
            call xermsg ('molden_mod', 'print_energy_sorted_orbital_table', &
                         'The Molden object has not been initialized.', 1, 1)
         end if

         if (this % io /= 1) then
            call xermsg ('molden_mod', 'print_energy_sorted_orbital_table', &
                         'The Molden file has been associated for output only so no orbital data is available.', 2, 1)
         end if

         if (.not. this % contains_mo) then
            call xermsg ('molden_mod', 'print_energy_sorted_orbital_table', &
                         'The Molden file does not contain molecular orbitals.', 3, 1)
         end if

         n_tgt = sum(this%no_orbs)

         allocate(energies(n_tgt),num_sym(2,n_tgt),stat=err)
         if (err /= 0) call xermsg ('molden_mod', 'print_energy_sorted_orbital_table', 'Memory allocation failed.',err,1)

         n = 0
         do i=1,this%max_irr !over all symmetries
            do j=1,this%no_orbs(i)
               n = n + 1
               energies(n) = this%orbital_data(i)%energy(j)
               num_sym(1,n) = j
               num_sym(2,n) = i
            enddo !i
         enddo

         call print_orbital_table(energies,num_sym,n_tgt,this%max_irr,.true.)

         write(stdout,'("--------->","done:molden_input_obj:print_energy_sorted_orbital_table")')

   end subroutine print_energy_sorted_orbital_table

end module molden_mod
