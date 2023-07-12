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
module ukrmol_interface
   use const,   only: line_len, data_file_obj_id, sym_op_nam_len, one_electron_ints, overlap_ints, kinetic_ints, &
                      property_ints, nuc_rep_att_ints, one_elham, two_el_ints, one_p_sym_ints, two_p_sym_ints, ijkl_indices_header
   use mpi_mod, only: master, mpi_mod_bcast, mpi_mod_start, mpi_mod_finalize
   use precisn, only: wp
   use utils,   only: xermsg

   use atomic_basis_mod
   use molecular_basis_mod
   use parallel_arrays
   use integral_storage_mod
   use mpi_memory_mod

   private

   !> These two variables are used in SCATCI_ROUTINES to identify the format of the integral file on input.
   public data_file_obj_id, line_len

   !> The following routines are used by UKRmol-in an UKRmol-out programs (SCATCI, DENPROP, CDENPROP) and routines (SWINTERF).
   public READ_UKRMOLP_BASIS, GET_GEOM, GET_NAME_SYM, READ_UKRMOLP_INTS, READ_UKRMOLP_PROPERTY_INTS, TMG_UKPLUS, EVAL_AMPLITUDES
   public UKP_PREAMP, UKP_READAMP, GET_INTEGRAL, GET_INTEGRALS
   public CONSTRUCT_PINTEGRALS, WRITE_MOLDEN_DYSON_ORBITALS,DESTROY_UKRMOLP, START_MPI, FINALIZE_MPI
   public UKRMol_Get_Integral_Index
   public Get_1B_Integral, GET_2B_INTEGRAL

   !> These arrays are visible but not changeable from outside this module.
   type(p2d_array_obj), target :: one_electron_integrals, one_positron_integrals, two_electron_integrals

   !> Integral options as read-in by READ_UKRMOLP_PROPERTY_INTS or READ_UKRMOLP_INTS.
   type(integral_options_obj) :: options

   !> Numbers mapping the type of integral with its column number in one_electron_integrals%a(:,:) and one_positron_integrals%a(:,:)
   integer :: one_elham_column = -1, kei_column        = -1, nai_column     = -1
   integer :: property_column  = -1, one_posham_column = -1, overlap_column = -1
   
   !> Number of different (l,m) properties stored in one_electron_integrals.
   integer :: number_of_properties = -1

   !> Basis set atomic orbitals. Read-in by READ_UKRMOLP_BASIS.
   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis

   !> Basis set of molecular orbitals. Read-in by READ_UKRMOLP_BASIS.
   type(molecular_orbital_basis_obj), public :: molecular_orbital_basis

   !> Number of molecular orbitals in each symmetry. Determined by READ_UKRMOLP_BASIS.
   integer, allocatable :: num_orbs_sym(:)

   !> Array of dimensions (1:3,1:lmq) where lmq is the total number of (l,m) properties for which the integrals have been calculated. It is determined by READ_UKRMOLP_PROPERTY_INTS.
   integer, allocatable :: msym(:,:)

   !> Channel information (l,m,IRR) and the number of orbitals in each IRR. Determined by EVAL_AMPLITUDES.
   integer, allocatable :: continuum_channels_m_l_irr(:,:), ncont(:)

   !> Amplitudes at r=a for each orbital. Determined by EVAL_AMPLITUDES.
   real(kind=cfp), allocatable :: amplitudes(:,:)

   character(len=100), allocatable :: Table_1BINT_Labels(:)


contains

   !> This is used to start MPI in SWINTERF and HAMDIAG. In the serial integral code this is just a dummy routine that does nothing.
   subroutine START_MPI
      use mpi_mod
      implicit none

         if (.not. mpi_started) call mpi_mod_start

   end subroutine START_MPI

   !> This is used to start MPI in SWINTERF and HAMDIAG. In the serial integral code this is just a dummy routine that does nothing.
   subroutine FINALIZE_MPI
      use mpi_mod
      implicit none

         call mpi_mod_finalize

   end subroutine FINALIZE_MPI

   !> Converts the integer number XXX corresponding to the fort.XXX unit into actual string 'fort.XXX'.
   subroutine unit_number_2_file_name(lu,file_name)
      implicit none
      integer, intent(in) :: lu
      character(len=line_len), intent(out) :: file_name

      character(len=line_len) :: cnum

         file_name = 'fort.'
         write(cnum,'(i0)') lu 
         file_name = trim(file_name)//trim(adjustl(cnum))

   end subroutine unit_number_2_file_name

   subroutine DESTROY_UKRMOLP()
     implicit none
     integer :: dum

        dum = molecular_orbital_basis%final()
        dum = two_electron_integrals%final()
        dum = one_electron_integrals%final()

   end subroutine DESTROY_UKRMOLP

   !> Reads-in the AO and MO basis sets from the fort file with a given number.
   subroutine READ_UKRMOLP_BASIS(NFTINT,InputFileName)
      implicit none
      integer                   , intent(in) :: NFTINT
      character(len=*), optional, intent(in) :: InputFileName
      
      type(integral_storage_obj)        :: mo_ijkl_storage
      type(data_header_obj)             :: header
      character(len=line_len)           :: file_name
      integer :: irr, err, pos, lunit, uid
      logical :: op

      if(present(InputFileName))then
         file_name=trim(adjustl(InputFileName))
      else
         !we assume that the integrals are saved on the unit number 'nfti'
         !so we need to construct the appropriate file name.
         call unit_number_2_file_name(NFTINT,file_name)
      endif
         
      !make sure the integral file is closed before we attempt to open it: this can be
      !an issue when calling swinterf repeatedly in one program.
      inquire(file=file_name,opened=op,number=uid)
      if (op) close(uid)
      
      !read-in all basis sets and orbital data
      call atomic_orbital_basis%read(file_name)
      molecular_orbital_basis%ao_basis => atomic_orbital_basis
      call molecular_orbital_basis%read(file_name)
      
      err = mo_ijkl_storage%init(disk=file_name)
      if (err /= 0) then
         call xermsg ('ukrmol_interface', 'READ_UKRMOLP_BASIS', 'Error initializing the integrals storage object.', err, 1)
      end if
      
      !we need to get the full header corresponding to the ijkl_indices.
      err = mo_ijkl_storage % integral_file % get_header_containing (header, &
           molecular_orbital_basis % get_basis_name(), ijkl_indices_header)
      if (err /= 0) then
         call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', &
              'Error locating ijkl_indices; '//&
              'see data_header_obj%get_header_containing for details.', err, 1)
      end if
      
      call mpi_mod_bcast(header%first_record,master)
      lunit = mo_ijkl_storage%integral_file%get_unit_no()
      !write(stdout,*) lunit
      call mpi_mod_bcast(lunit,master)
      !write(stdout,"('LETS DO THIS BAYBE Lunit = ',2i16)") lunit,header%first_record
      call molecular_orbital_basis%read_ijkl_indices(lunit,header%first_record,pos)
      
      call mo_ijkl_storage%final

      if (allocated(num_orbs_sym)) deallocate(num_orbs_sym)
      allocate(num_orbs_sym(molecular_orbital_basis%no_irr),stat=err)
      if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_BASIS','Memory allocation failed.',err,1)
      
      do irr=1,molecular_orbital_basis%no_irr
         num_orbs_sym(irr) = molecular_orbital_basis%get_number_of_orbitals(irr)
      enddo !i
      
    end subroutine READ_UKRMOLP_BASIS

   !> Returns the nuclear data for the problem as contained in the molecular_orbital_basis object.
   !> \warning This routine ALWAYS returns double precision values due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   subroutine GET_GEOM(nnuc,cname,xnuc,ynuc,znuc,charge)
      implicit none
      integer, intent(out) :: nnuc
      real(kind=wp), allocatable, dimension(:), intent(out) :: charge, xnuc, ynuc, znuc
      character(len=8), allocatable, dimension(:), intent(out) :: cname

      integer :: i, j, k, err

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'GET_GEOM', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         nnuc = molecular_orbital_basis%symmetry_data%no_nuc
         allocate(charge(nnuc), xnuc(nnuc), ynuc(nnuc), znuc(nnuc), cname(nnuc),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','GET_GEOM','Memory allocation failed.',err,1)

         associate(nucleus => molecular_orbital_basis%symmetry_data%nucleus, nuclear_data => molecular_orbital_basis%symmetry_data)
            k = 0
            do i=1,nuclear_data%no_sym_nuc       !over all symmetrically non-redundant nuclei
               k = k + 1
               !data for the non-redundant nucleus itself
               write(cname(k),'(a2,i2)') nuclear_data%non_red_nuc(i)%name, i
               xnuc(k) = nuclear_data%non_red_nuc(i)%center(1)
               ynuc(k) = nuclear_data%non_red_nuc(i)%center(2)
               znuc(k) = nuclear_data%non_red_nuc(i)%center(3)
               charge(k) = nuclear_data%non_red_nuc(i)%charge
               do j=1,nuclear_data%no_eqv_nuc(i) !over all symmetrical partners
                  k = k + 1
                  !data for the symmetrical partner
                  write(cname(k),'(a2,i2)') nuclear_data%eqv_nuc(i,j)%name, i
                  xnuc(k) = nuclear_data%eqv_nuc(i,j)%center(1)
                  ynuc(k) = nuclear_data%eqv_nuc(i,j)%center(2)
                  znuc(k) = nuclear_data%eqv_nuc(i,j)%center(3)
                  charge(k) = nuclear_data%eqv_nuc(i,j)%charge
               enddo !j
            enddo !i
         end associate

   end subroutine GET_GEOM

   !> Returns the number of symmetries and some other basic parameters as obtained from the initialized object molecular_orbital_basis and property_integrals.
   subroutine GET_NAME_SYM(name,no_sym,nob,nlmq)
      implicit none
      character(len=132), intent(out) :: name
      integer, intent(out) :: no_sym, nob(:), nlmq

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'GET_NO_SYM', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         no_sym = molecular_orbital_basis%no_irr
         name = property_ints

         if (size(nob) < no_sym) call xermsg('ukrmol_interface','GET_NO_SYM','The array NOB on input is too small.',2,1)
         nob(1:no_sym) = num_orbs_sym(1:no_sym)

         if (number_of_properties > 0) then
            nlmq = number_of_properties !how many (l,m) property integrals we have (with l=0,1,...).
         else
            call xermsg('ukrmol_interface','GET_NO_SYM','Cannot supply nlmq: property have not been read-in.',3,1)
         endif

   end subroutine GET_NAME_SYM

   !> This routine reads-in all one electron and two electron integrals needed for SCATCI calculations. It is an analogoue of the INTSIN routine from SCATCI.
   !> \warning This routine ALWAYS requires double precision value of scalem on input due to the need for compatibility with UKRmol-in/out which is always compiled using double precision. Note that 
   !> this results in reducing the values in ke_integrals%a into double precision.
   !> NEW FEATURE: WORKS WITH SHARED MEMORY OF THE MPI-3.0 STANDARD. This means that all integrals are read-in the shared mode so in the MPI regime there is only one copy of the integral arrays PER NODE.
   !> If we do not have the shared memory capability then the code resorts to the local mode, i.e. each MPI task keeps its own copy of all integral arrays.
   subroutine READ_UKRMOLP_INTS( &
        nfte, nfti, lembf, nint1e, nint2e, nocsf, nfta, &
        isymtp, nsym1, nob1, iposit, scalem, name, nalm, qmoln, &
        InputFileName )
     use symmetry
     implicit none
     integer                   , intent(in)    :: nfte,nfti,lembf,nocsf,nfta
     integer                   , intent(in)    :: isymtp,nsym1,nob1(nsym1),iposit
     integer                   , intent(inout) :: nint1e, nint2e
     real(kind=wp)             , intent(in)    :: scalem
     character(len=120)        , intent(in)    :: name
     logical                   , intent(in)    :: qmoln
     integer                   , intent(inout) :: nalm
     character(len=*), optional, intent(in)    :: InputFileName
     
     !local variables
     character(len=line_len) :: file_name
     type(integral_storage_obj) :: mo_integral_storage
     type(integral_storage_obj), target :: tgt_storage
     type(data_header_obj) :: header
     integer :: i, j, err, nnuc, min_ind, max_ind
     real(kind=cfp) :: inv_scale
     real(kind=wp) :: nuc_repulsion, distance
     real(kind=wp), dimension(41) :: dtnuc
     integer, dimension(20) :: nhe
     integer, parameter :: number_of_blocks = 0
     logical, parameter :: tgt_is_local = .true.
     character(len=line_len), allocatable :: column_descriptor(:)

     !error reporting flag used in UKRmol but not here: we report errors immediately
     nalm = 0

     if (qmoln) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS',&
          'Use of the QMOLN flag has not been implemented.',1,1)

     !we assume that the integrals are saved on the unit number 'nfti'
     call READ_UKRMOLP_BASIS(NFTI,InputFileName)

     !check for data input errors
     !nsym1 is the number of the last IRR where we have some orbitals
     if (molecular_orbital_basis%no_irr < nsym1) then 
        print *,molecular_orbital_basis%no_irr, nsym1
        call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', &
             'The number of IRRs obtained from the orbital basis data is not compatible &
             &with the one given by SCATCI.', 2, 1)
     endif

     do i=1,nsym1
        if (num_orbs_sym(i) .ne. nob1(i)) then
           write(nfta,451)
           write(nfta,452) (num_orbs_sym(j),j=1,nsym1)
           write(nfta,452) (nob1(j),j=1,nsym1)
451        format(/' MISMATCH IN NOB BETWEEN INTEGRALS AND FORMULAE')
452        format(' NOB=',20i5)
           call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', &
                'The number of orbitals contained in the orbital basis data is not equal &
                &to the one expected by SCATCI.', 3, 1)
        endif
     enddo

     !open the data file and analyse the headers manually
     if(present(InputFileName))then
        file_name=trim(adjustl(InputFileName))
     else
        call unit_number_2_file_name(NFTI,file_name)
     endif
     err = mo_integral_storage%init(disk=file_name)
     if (err /= 0) then
        call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', 'Error initializing the integrals storage object.', err, 1)
     end if

     write(nfta,1600) nint1e, nint2e
1600 format(/,' On input SCATCI expects:',/,' 1-electron integrals, NINT1e =',i11,/,' 2-electron integrals, NINT2e =',i11)

     !read-in the 1-electron integrals into one_electron_integrals
     err = tgt_storage%init(memory=one_electron_integrals)
     if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error initializing tgt_storage.',err,1)

     !we need to get the full header corresponding to the transformed overlap and kinetic energy integrals.
     err = mo_integral_storage % integral_file % get_header_containing(header, &
          molecular_orbital_basis % get_basis_name(), &
          one_p_sym_ints)
     if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error locating the transformed 1-electron integrals; &
          &see data_header_obj%get_header_containing for details.',err,1)

     call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
     call tgt_storage%final

     write(nfta,'(/,"One electron integrals read-in.")')

     !Column number in one_electron_integrals%a corresponding to the 1-electron Hamiltonian integrals: these integrals must ALWAYS be calculated.
     one_elham_column = one_electron_integrals%find_column_matching_name(one_elham)
     kei_column = one_electron_integrals%find_column_matching_name(kinetic_ints)
     nai_column = one_electron_integrals%find_column_matching_name(nuc_rep_att_ints)

     property_column =one_electron_integrals%find_column_matching_name(property_ints)
     overlap_column =one_electron_integrals%find_column_matching_name(overlap_ints)
     
     !.. Initialize the table with the correspondence between integrals IDs and labels
     min_ind = min(one_elham_column,kei_column,nai_column,property_column,overlap_column)
     max_ind = max(one_elham_column,kei_column,nai_column,property_column,overlap_column)
     if(.not.allocated(Table_1BINT_Labels))allocate(Table_1BINT_Labels(min_ind:max_ind))
     Table_1BINT_Labels = " "
     Table_1BINT_Labels( one_elham_column) = one_elham
     Table_1BINT_Labels( kei_column      ) = kinetic_ints
     Table_1BINT_Labels( nai_column      ) = nuc_rep_att_ints
     Table_1BINT_Labels( property_column ) = property_ints
     Table_1BINT_Labels( overlap_column  ) = overlap_ints

     if (scalem .ne. 1.0_cfp) then
        write(nfta,'("Kinetic energy integrals will be rescaled (divided) by: ",e25.15)') scalem
        !Column number in one_electron_integrals%a corresponding to the kinetic energy integrals.
        kei_column = one_electron_integrals%find_column_matching_name(kinetic_ints)

        !Column number in one_electron_integrals%a corresponding to the nuclear attraction integrals.
        nai_column = one_electron_integrals%find_column_matching_name(nuc_rep_att_ints)

        !Rescale the kinetic energy integrals and recompute the 1-electron Hamiltonian integrals for the electron case.
        inv_scale = 1.0_cfp/scalem

        call one_electron_integrals%multiply_column(kei_column,inv_scale)
        !Direct summation below implies that the indexing of the kinetic and nuclear integrals is the same
        call one_electron_integrals%combine_columns(kei_column,'+',nai_column,one_elham_column)

     endif

     !Construct the one-positron integrals: K - V
     if (iposit .ne. 0) then
        write(nfta,'("Recomputing the 1-particle Hamiltonian integrals for the positron case...")')

        !Column number in one_electron_integrals%a corresponding to the kinetic energy integrals.
        kei_column = one_electron_integrals%find_column_matching_name(kinetic_ints)

        !Column number in one_electron_integrals%a corresponding to the nuclear attraction integrals.
        nai_column = one_electron_integrals%find_column_matching_name(nuc_rep_att_ints)

        call one_electron_integrals%get_column_descriptor(column_descriptor)

        !Allocate space for the 1p Hamiltonian integrals
        err = one_positron_integrals % init(size(one_electron_integrals % a, 1), &
             size(one_electron_integrals % a, 2), &
             number_of_blocks, &
             column_descriptor) 
        if (err /= 0) then
           call xermsg ('molecular_orbital_basis_obj', 'READ_UKRMOLP_INTS', &
                'Initialization of one_positron_integrals failed; see p2d_array%init.', err, 1)
        end if

        one_posham_column = one_electron_integrals%find_column_matching_name(one_elham)

        one_positron_integrals%a = one_electron_integrals%a

        !Direct subtraction below implies that the indexing of the kinetic and nuclear integrals is the same
        call one_positron_integrals%combine_columns(kei_column,'-',nai_column,one_posham_column)

        write(nfta,'("...done")')

     endif

     !read-in the 2p integrals
     err = tgt_storage%init(memory=two_electron_integrals)
     if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error initializing tgt_storage.',err,1)

     !we need to get the full header corresponding to the transformed 2-electron integrals.
     err = mo_integral_storage % integral_file % get_header_containing(header, &
          molecular_orbital_basis % get_basis_name(), &
          two_p_sym_ints)
     if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error locating the transformed 2-electron integrals; &
          &see data_header_obj%get_header_containing for details.',err,1)

     call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
     call tgt_storage%final

     write(nfta,'(/,"Two electron integrals read-in.")')

     !Now update the integral counters by the actual number of integrals read-in
     nint1e = size(one_electron_integrals%a,1)
     nint2e = size(two_electron_integrals%a,1)
     write(nfta,1700) nint1e, nint2e
1700 format(/,' UKRMol+ integrals read successfully:',&
          /,' 1-electron integrals, NINT1e =',i10,&
          /,' 2-electron integrals, NINT2e =',i10)

     !Compute the nuclear repulsion energy
     associate(nucleus => molecular_orbital_basis%symmetry_data%nucleus, nuclear_data => molecular_orbital_basis%symmetry_data)
     nuc_repulsion = 0.0_cfp
     nnuc = nuclear_data%no_nuc
     do j=2,nnuc
        do i=1, j-1
           if (nucleus(j)%charge*nucleus(i)%charge .eq. 0.0_cfp) cycle
           distance = sqrt(dot_product(nucleus(j)%center-nucleus(i)%center,nucleus(j)%center-nucleus(i)%center))
           if (distance .eq. 0.0_cfp) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Nuclei too close.',10,1)
           nuc_repulsion = nuc_repulsion + nucleus(j)%charge*nucleus(i)%charge/distance
        end do
     end do
     !The assignments to dtnuc effectively replicate the equivalence statements in INTSIN in scatci_build. 
     !However we drop the assignment to dtnuc(22) which is set only when Alchemy integrals are used and seems completely redundant.
     dtnuc(:) = 0.0_cfp
     dtnuc(1) = nuc_repulsion
     dtnuc(2) = nucleus(1)%charge !is this redundant?
   end associate

   write(nfta,1510) nuc_repulsion
1510 format(/,10x,'Nuclear repulsion energy = ',f15.7)

   !The assignment to nhe effectively replicates the equivalence statement in INTSIN in scatci_build
   nhe(:) = 0
   nhe(1) = num_orbs_sym(1) !is this redundant?

   !Write header on energy matrix file  
   open(unit=nfte,form='unformatted')
   write(nfte) nocsf, lembf, 0, nocsf, 0, nsym1, 0, 0, 0, 0, nnuc, 0, name, nhe, dtnuc

   call mo_integral_storage%final

   write(nfta,'(/,10x,"Finished reading UKRMol+ integrals.")')

 end subroutine READ_UKRMOLP_INTS


   !> LA: Added August 3, 2020
   integer function UKRMol_Get_Integral_Index( Integral_Label ) result( Integral_Index )
     implicit none
     character(len=*), intent(in) :: Integral_Label
     !Column number in one_electron_integrals%a corresponding to the given integral label
     Integral_Index = one_electron_integrals%find_column_matching_name( Integral_Label )
   end function UKRMol_Get_Integral_Index
   !====


   !> Reads-in the basis sets and the property integrals.
   subroutine READ_UKRMOLP_PROPERTY_INTS(NFTINT,IWRITE,InputFileName)
      use symmetry
      implicit none
      integer                   , intent(in) :: NFTINT, IWRITE
      character(len=*), optional, intent(in) :: InputFileName

      !local variables
      character(len=line_len)              :: file_name
      type(integral_storage_obj)           :: mo_integral_storage
      type(integral_storage_obj), target   :: tgt_storage
      type(data_header_obj)                :: header
      integer                              :: i, err, l, m
      character(len=sym_op_nam_len)        :: nam
      integer, parameter                   :: no_blocks = 0
      logical, parameter                   :: tgt_is_local = .true.
      character(len=line_len), allocatable :: column_descriptor(:)

      !we assume that the integrals are saved on the unit number 'nfti'
      call READ_UKRMOLP_BASIS(NFTINT,InputFileName)
      
      !open the data file and analyse the headers manually
      if(present(InputFileName))then
         file_name=trim(adjustl(InputFileName))
      else
         call unit_number_2_file_name(NFTINT,file_name)
      endif
      err = mo_integral_storage%init(disk=file_name)
      if (err /= 0) then
         call xermsg ('ukrmol_interface', 'READ_UKRMOLP_PROPERTY_INTS', &
              'Error initializing the integrals storage object.', err, 1)
      end if
      
      !we need to get the full header corresponding to the transformed property integrals.
      err = mo_integral_storage % integral_file % get_header_containing(header, &
           molecular_orbital_basis % get_basis_name(), &
           one_p_sym_ints)
      if (err /= 0) then
         call xermsg ('ukrmol_interface', 'READ_UKRMOLP_PROPERTY_INTS', &
              'Error locating the transformed 1-electron integrals; &
              &see data_header_obj%get_header_containing for details.', err, 1)
      end if
      
      !read-in all 1-electron integrals
      err = tgt_storage%init(memory=one_electron_integrals)
      if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_PROPERTY_INTS',&
           'Error initializing tgt_storage.',err,1)
      
      call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
      call tgt_storage%final
      
      write(IWRITE,'(/,"1-electron integrals read-in.")')
      
      property_column = one_electron_integrals%find_column_matching_name(property_ints)
      
      write(*,*) property_column,"property_column",property_ints
      
      call one_electron_integrals%get_column_descriptor(column_descriptor)
      
      number_of_properties = 0
      do i=1,size(column_descriptor)
         if (column_descriptor(i) .eq. property_ints) number_of_properties = number_of_properties + 1
      enddo
      
      write(IWRITE,'(/,"Number of types of property integrals: ",i4)') number_of_properties
      write(*,'(/,"Number of types of property integrals: ",i4)') number_of_properties
      
      if (allocated(msym)) deallocate(msym)
      allocate(msym(1:3,number_of_properties),stat=err)
      if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_PROPERTY_INTS',&
           'Memory allocation 2 error',err,1)
      
      !build the (l,m) property numbers up to the last property.
      i = 1
      l = 0
      m = 0 
      do
         msym(1,i) = m
         msym(2,i) = l
         
         !Find out the symmetry of the property operator (l,m) which is the same
         !as the IRR of the spherical harmonic in the P-G of the molecule.
         msym(3,i) = molecular_orbital_basis%symmetry_data%sph_harm_pg_sym(l,m,nam)
         
         write(stdout,'(/,10x,"Property (l,m): (",i4,",",i4,")")') l,m
         write(stdout,'(10x,"Property symmetry: ",i2)') msym(3,i)
         
         if (i .eq. number_of_properties) exit
         
         if (m .eq. l) then
            l = l + 1
            m = -l
         else
            m = m + 1
         endif
         
         i = i + 1
      enddo
      
      !Make sure we always close the integrals file to ensure we don't get errors in CDENPROP when the properties are being read-in multiple times.
      call mo_integral_storage%final

   end subroutine READ_UKRMOLP_PROPERTY_INTS

   !> Returns a one electron/positron integral (K+-V) or two-electron integral given the basis function indices and the positron flag. If c,d == 0 then the one electron/positron integral is returned.
   !> If positron_flag > 0 then the positronic one-particle integral (K-V) is returned.
   !> \warning This routine ALWAYS returns double precision value due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   !> NEW FEATURE: WORKS WITH SHARED MEMORY OF THE MPI-3.0 STANDARD
   function GET_INTEGRAL(a,b,c,d,positron_flag)
      implicit none
      integer, intent(in) :: a,b,c,d, positron_flag
      real(kind=wp) :: GET_INTEGRAL

      integer :: ind(1), two_ind(1:2,1), four_ind(1:4,1)

         GET_INTEGRAL = 0.0_wp
         if (c .le. 0 .or. d .le. 0) then !one particle case
            two_ind(1,1) = a
            two_ind(2,1) = b
            ind(1:1) = molecular_orbital_basis%integral_index(one_elham,two_ind)
write(*,*) "hola",size(one_positron_integrals%a(1,:)),size(one_positron_integrals%a(:,1))
            if (positron_flag .ne. 0) then !positron case: K-V, where V is calculated for electrons
               GET_INTEGRAL = one_positron_integrals%a(ind(1),one_posham_column)
            else !electron case: K+V, where V is calculated for electrons
               GET_INTEGRAL = one_electron_integrals%a(ind(1),one_elham_column)
            endif
            !write(*,'("1el",2i3,i,e25.15)') a,b,ind(1),GET_INTEGRAL
            !write(100,'("1el",2i3,e25.15)') a,b,GET_INTEGRAL
         else !two electron case
            four_ind(1,1) = a
            four_ind(2,1) = b
            four_ind(3,1) = c
            four_ind(4,1) = d
            ind(1:1) = molecular_orbital_basis%integral_index(two_el_ints,four_ind)
            if (ind(1) > 0) GET_INTEGRAL = two_electron_integrals%a(ind(1),1)
            write(*,'("2el",4i3,i0,e25.15)') a,b,c,d,ind(1),GET_INTEGRAL
           ! write(*,*) two_el_ints,four_ind
            !write(101,'("2el",4i3,e25.15)') a,b,c,d,GET_INTEGRAL
         endif

   end function GET_INTEGRAL

   !This routine is an extension of the GET_INTEGRAL. it allows to acces also to
   !the multipole elements. Note that positron_flag has changed to
   !position_flag:
   !
   !position_flag=-4   ===>  positron Hamiltonian
   !position_flag=-2   ===>  Kinetic Energy
   !position_flag=-1   ===>  Nuclear attraction
   !position_flag=0    ===>  electron Hamiltonian
   !positron_flag>0    ===>  property integrals
   !last index Get_integrals   l  m
   !position_flag=1            0  0    Since we are dealing with solid spherica it is also equal to the overlap. 
   !position_flag=2            1 -1
   !position_flag=3            1  0
   !position_flag=4            1  1
   !position_flag=5            2 -2
   !position_flag=6            2 -1
   !position_flag=7            2  0
   !position_flag=8            2  1
   !position_flag=9            2  2    and so on up to l=20 which is set at line 265 of scatci_integrals.f90 

   function GET_INTEGRALS(a,b,c,d,position_flag)
      implicit none
      integer, intent(in) :: a,b,c,d,position_flag
      real(kind=wp) :: GET_INTEGRALS

      integer :: ind(1), two_ind(1:2,1), four_ind(1:4,1)
!         write(*,*) "one_elham_column",one_elham_column,property_column,one_elham,kei_column,kinetic_ints
         GET_INTEGRALS = 0.0_wp
         if (c .le. 0 .or. d .le. 0) then !one particle case
            two_ind(1,1) = a
            two_ind(2,1) = b


!one_electron_ints, overlap_ints, kinetic_ints, &
!                      property_ints, nuc_rep_att_ints, one_elham

            if (position_flag .eq. -4) then !positron case: K-V, where V is calculated for electrons
            ind(1:1) = molecular_orbital_basis%integral_index(one_elham,two_ind)
               GET_INTEGRALS = one_positron_integrals%a(ind(1),one_posham_column)
            endif
            if (position_flag .eq. -2) then  !Kinetic energy integrals
            ind(1:1) = molecular_orbital_basis%integral_index(kinetic_ints,two_ind)
               GET_INTEGRALS = one_electron_integrals%a(ind(1),kei_column)
            endif
            if (position_flag .eq. -1) then  !Nuclear attraction integral
            ind(1:1) = molecular_orbital_basis%integral_index(nuc_rep_att_ints,two_ind)
               GET_INTEGRALS = one_electron_integrals%a(ind(1),nai_column)
            endif
            if (position_flag .eq. 0) then  !electron case: K+V, where V is calculated for electrons
               ind(1:1) = molecular_orbital_basis%integral_index(one_elham,two_ind)
               GET_INTEGRALS = one_electron_integrals%a(ind(1),one_elham_column) 
            endif
            if (position_flag .gt. 0) then  !electron case: property integrals with number qlm=position_flag
            ind(1:1) = molecular_orbital_basis%integral_index(property_ints,two_ind)
               GET_INTEGRALS = one_electron_integrals%a(ind(1),property_column-1+position_flag)
            endif
            !write(*,'("1el",2i3,i,e25.15)') a,b,ind(1),GET_INTEGRALS,property_column,position_flag
            !write(*,*) property_column,position_flag,one_posham_column,kei_column,nai_column,one_elham_column
                !write(100,'("1el",2i3,e25.15)') a,b,GET_INTEGRAL
         else !two electron case
            four_ind(1,1) = a
            four_ind(2,1) = b
            four_ind(3,1) = c
            four_ind(4,1) = d
            ind(1:1) = molecular_orbital_basis%integral_index(two_el_ints,four_ind)
            if (ind(1) > 0) GET_INTEGRALS = two_electron_integrals%a(ind(1),1)
            !write(*,'("2el",4i3,i,e25.15)') a,b,c,d,ind(1),GET_INTEGRALS
            !write(101,'("2el",4i3,e25.15)') a,b,c,d,GET_INTEGRAL
         endif

   end function GET_INTEGRALS


   !.. Returns a general one-body integrals. In the case of Multiples, one must specify
   !   the optional parameters l, m
   !..
   real(kind=wp) function GET_1B_INTEGRAL( a, b, Integral_ID, l, m ) result( qres )
     implicit none
     integer,           intent(in) :: a, b, Integral_ID
     integer, optional, intent(in) :: l, m
     integer :: ind(1), two_ind(1:2,1), lmIndex
     qres = 0.0_wp
     two_ind(1,1) = a
     two_ind(2,1) = b
     lmIndex  = 0
     if(present(l).and.present(m)) lmIndex  = l**2 + l + m 
     ind(1:1) = molecular_orbital_basis%integral_index(trim(Table_1BINT_Labels(Integral_ID)),two_ind)
     qres     = one_electron_integrals%a( ind(1), Integral_ID + lmIndex )
   end function GET_1B_INTEGRAL


   real(kind=wp) function GET_2B_INTEGRAL(a,b,c,d) result( qres )
     implicit none
     integer, intent(in) :: a,b,c,d
     integer :: ind(1), four_ind(1:4,1)
     qres = 0.0_wp
     four_ind(1,1) = a
     four_ind(2,1) = b
     four_ind(3,1) = c
     four_ind(4,1) = d
     ind(1:1) = molecular_orbital_basis%integral_index(two_el_ints,four_ind)
     if (ind(1) > 0) qres = two_electron_integrals%a(ind(1),1)
   end function GET_2B_INTEGRAL




   !> Multiplies the density matrix elements with the set of property integrals with (l,m) given by the linear index lmq = l*l+l+m+1. This routine may only be called following a call to 
   !> READ_UKRMOLP_PROPERTY_INTS. We assume that this routine is called with the first call having lmq=1. This is needed for the index function initialization. The input values are: the wavefunction symmetry
   !> difference (mdel), the number of elements in the density matrix (lden), the property l,m,q sequence number (lmq) and the block data constructed in DENPROP/TMGP. The output are: the l,m,q values 
   !> corresponding to lmq and the value for the electronic contribution to the property (prlmq).
   !> \warning This routine ALWAYS uses double precision values of the property integrals due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   subroutine TMG_UKPLUS(mdel,block_data,no_blocks,nob,den,lden,prlmq,lmq,l,m,q)
      implicit none
      integer, intent(in) :: mdel, lden, lmq, block_data(:,:), nob(:), no_blocks
      integer, intent(out) :: l,m,q
      real(kind=wp), dimension(lden) :: den
      real(kind=wp) :: prlmq
 
      integer :: isym, jsym, nsym, i, orb_i, orb_j, last, o_sym_i_start, o_sym_j_start, two_ind(1:2,1:1), ind(1:1), block

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'TMG_UKPLUS', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         l = msym(2,lmq)
         m = abs(msym(1,lmq))
         if (m .eq. 0) then 
            q = 0
         else
            q = sign(1,msym(1,lmq))
         endif

         nsym = molecular_orbital_basis%no_irr

         prlmq = 0.0_wp

         !Here we loop over the pairs of symmetries of orbitals that contribute to the density matrix. The order of the blocks of symmetries is important since this order corresponds to the order of the
         !blocks (and elements) of the density matrix (indexed by 'i').
         i = 0
         do block = 1,no_blocks
            if (block_data(1,block) .eq. mdel) then

               isym = block_data(2,block)
               jsym = block_data(3,block)

               !offsets for the orbital indices from the symmetries isym,jsym
               o_sym_i_start = sum(num_orbs_sym(1:isym-1))
               o_sym_j_start = sum(num_orbs_sym(1:jsym-1))
   
               !Note that the loops over orbitals are loops only over the orbitals that are present in the CONGEN input, i.e. we don't loop over all orbitals in the basis set: that number might be different
               !to the one used to construct the wavefunctions. The order of the loops orb_j, orb_i is also important:
               !we loop over the property integrals in the order in which the density matrix elements DEN(1:lden) were constructed (see DENPROP routines SETMBA, RWDIJ, INDEX1).
               last = nob(isym)
               do orb_j = 1,nob(jsym)
                  if (isym .eq. jsym) last = orb_j !only unique combinations of orbitals are used in the density matrix expressions
                  do orb_i = 1,last
                     i = i + 1
                     if (i .le. lden .and. DEN(i) .ne. 0.0_wp) then !there may be a non-zero contribution to the property
                        two_ind(1,1) = o_sym_i_start + orb_i !the index of the orbital in the set of orbitals for all symmetries, sim. for the second orbital.
                        two_ind(2,1) = o_sym_j_start + orb_j
                        ind(1:1) = molecular_orbital_basis%integral_index(property_ints,two_ind) !index of the property integral <orb_i|(l,m)|orb_j> for the orbitals from the symmetries isym, jsym.
                        prlmq = prlmq + DEN(i)*one_electron_integrals%a(ind(1),property_column-1+lmq)
                     endif
                  enddo !orb_j
               enddo !orb_i

            endif
         enddo

   end subroutine TMG_UKPLUS

   !> This routine interfaces with CDENPROP. It transfers the necessary orbital data into the argument variables. Most importantly it transfers the property integrals and their
   !> indices. We use READ_UKRMOLP_PROPERTY_INTS to read the integrals from disk. CDENPROP uses GAUSPROP indexing of the property integrals. Therefore we need to make sure that the
   !> indices of the property integrals obtained using this code match the GAUSPROP indexing scheme. The principle of computing the indices here is to loop over all pairs of orbitals
   !> in the basis and for each pair calculate the pair index in the GAUSPROP style (index_ij). This index is then used to store in inverted_indexv(index_ij,property_lm) the index of
   !> the corresponding property integral which is stored in xintegrals.
   !> \warning This routine ALWAYS returns double precision values due to the need for compatibility with CDENPROP which is always compiled using double precision.
   subroutine construct_pintegrals (nftprop,iwrite,no_of_integrals,no_of_properties,nilmq,lp,mp,qp,property_name,&
                                    nob,mob,mpob,inverted_indexv,xintegrals)
      use special_functions, only: ipair
      implicit none
      integer, intent(in) :: nftprop, iwrite
      integer, intent(out) :: no_of_integrals,no_of_properties
      integer, allocatable :: nilmq(:), lp(:), mp(:), qp(:), nob(:), mob(:), mpob(:), inverted_indexv(:,:)
      character(len=8), allocatable :: property_name(:)
      real(kind=wp), allocatable, intent(out) :: xintegrals(:,:)

      integer :: nobt, lmq, i, j, two_ind(1:2,1:1), ind(1:1), last, k, m, integral_index, sym_i, sym_j, iorb, jorb, &
                 relative_index, err
      integer, allocatable :: istart(:,:)

         !Read all basis sets and integrals into memory
         call READ_UKRMOLP_PROPERTY_INTS(nftprop,iwrite)

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'construct_pintegrals', &
                         'The orbital basis set data has not been initialized.', 1, 1)
         end if

         nobt = sum(num_orbs_sym(:))
         allocate(nob(molecular_orbital_basis%no_irr),mob(nobt),mpob(nobt),&
                  istart(molecular_orbital_basis%no_irr,molecular_orbital_basis%no_irr),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','construct_pintegrals','Memory allocation 1 failed.',err,1)
         istart(:,:) = 0

         !Compute pointers to start of each symmetry block - GAUSPROP style. Adopted from CDENPROP routine: symmetry_block_to_integral_table
         last=0
         k=1
         do m=1,molecular_orbital_basis%no_irr
            istart(m,m)=last+1
            last=last+num_orbs_sym(m)*(num_orbs_sym(m)+1)/2  
         end do

         do k=1,molecular_orbital_basis%no_irr-1
            do m=1,molecular_orbital_basis%no_irr-k
               istart(m,m+k)=last+1
               last=last+num_orbs_sym(m)*num_orbs_sym(m+k)
            end do
         end do

         !get dimensions of the integral array:
         no_of_properties = number_of_properties
         no_of_integrals = size(one_electron_integrals%a,1)

         allocate(lp(no_of_properties), mp(no_of_properties), nilmq(no_of_properties), &
                  qp(no_of_properties),property_name(no_of_properties),&
                  inverted_indexv(no_of_integrals,no_of_properties),xintegrals(no_of_integrals,no_of_properties),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','construct_pintegrals','Memory allocation 2 failed.',err,1)
         inverted_indexv(:,:) = 0; xintegrals(:,:) = 0.0_wp
         nilmq(:) = no_of_integrals

         nob(1:molecular_orbital_basis%no_irr) = num_orbs_sym(1:molecular_orbital_basis%no_irr)
         do i=1,nobt
            mpob(i) = molecular_orbital_basis%get_index_within_symmetry(i)  !orbital index within its symmetry
            mob(i) = molecular_orbital_basis%get_orbital_symmetry(i)-1 != symmetry of the orbital - 1
         enddo

         !Construct indices of the integrals in the GAUSPROP style and put them into inverted_indexv.
         do i=1,nobt
            two_ind(1,1) = i
            do j=1,i
               two_ind(2,1) = j

               !Compute the index for the corresponding property integrals: one_electron_integrals%a(ind(1),:)
               ind(1:1) = molecular_orbital_basis%integral_index(property_ints,two_ind)
               if (ind(1) > no_of_integrals) print *,'error',i,j,ind(1),no_of_integrals

               !Use GAUSPROP ordering of symmetry blocks.
               if (mob(i) > mob(j)) then
                  iorb = j
                  jorb = i
               else
                  iorb = i
                  jorb = j
               endif

               sym_i = mob(iorb)
               sym_j = mob(jorb)
               if (sym_i .eq. sym_j) then
                  relative_index = ipair( max(mpob(iorb),mpob(jorb)) ) + min(mpob(iorb),mpob(jorb))
                  integral_index = istart(sym_i+1,sym_j+1) + relative_index-1
                  if (integral_index > size(inverted_indexv,1)) print *,'inderror',integral_index,size(inverted_indexv,1)
                  do lmq=1,no_of_properties
                     if (one_electron_integrals % a(ind(1), property_column-1+lmq) /= 0.0_cfp) then
                        inverted_indexv(integral_index,lmq) = ind(1)
                     end if
                  enddo
               elseif (sym_i < sym_j) then
                  relative_index = mpob(iorb) + (mpob(jorb)-1)*nob(sym_i+1)
                  integral_index = istart(sym_i+1,sym_j+1) + relative_index-1
                  if (integral_index > size(inverted_indexv,1)) then
                    print *,'inderror',integral_index,size(inverted_indexv,1),relative_index,istart(sym_i+1,sym_j+1)
                  end if
                  do lmq=1,no_of_properties
                     if (one_electron_integrals % a(ind(1), property_column-1+lmq) /= 0.0_cfp) then
                        inverted_indexv(integral_index,lmq) = ind(1)
                     end if
                  enddo
               endif
               
            enddo !j  
         enddo !i


         do lmq=1,no_of_properties

            lp(lmq) = msym(2,lmq) 
            mp(lmq) = abs(msym(1,lmq))
            if (mp(lmq) .eq. 0) then
               qp(lmq) = 0
            else 
               qp(lmq) = sign(1,msym(1,lmq))
            endif
            if (lp(lmq) .eq. 0) property_name(lmq) = 'overlap'
            if (lp(lmq) .eq. 1) property_name(lmq) = 'dipole'
            if (lp(lmq) .eq. 2) property_name(lmq) = 'quadpole'
            if (lp(lmq) > 2) property_name(lmq) = '>2 mpole'

            !We just copy the property integrals: no reordering is needed since the correct pointers to the integrals (xintegrals) are stored in inverted_indexv computed above.
            xintegrals(1:no_of_integrals,lmq) = one_electron_integrals%a(1:no_of_integrals,property_column-1+lmq)
         enddo !lmq

   end subroutine construct_pintegrals

   !> This routine takes the Dyson orbitals produced by CDENPROP and writes them out in the Molden file format along with the GTO basis set. The CDENPROP Dyson's are expressed as
   !> a linear combination of the target and continuum MOs. Additionally to producing the Molden file this routine always saves the complete set of Dyson orbitals (along with the AO
   !> basis set) into an object of orbital_data_obj type and then writes this to disk. This ensures that we always retain the full information about the Dyson orbitals (in case
   !> the AO basis contains high L that Molden cannot handle).
   !> The value of the R-matrix radius stored in 'options' is used to normalize the continuum functions to the R-matrix radius. These normalization coefficients are multiplied in with the orbital coefficients
   !> to produce the Molden file. The orbital coefficients saved in the UKRmol+ format are modified, too. Therefore when calculating the radial densities of the Dysons from the UKRmol+ file using the
   !> routine orbital_basis_data_obj%radial_charge_density the input value of rmat_radius for this routine must be set to a value .le. 0.0_cfp. See also orbital_basis_data_obj%radial_charge_density.
   !> \warning This routine multiplies in the spin-related CG coefficient (see commit 879 to cdenprop_io) but it is hard coded to the value 1/sqrt(2). This must be changed for your particular case.
   !> \warning This routine ALWAYS REQUIRES double precision values of the coefficients on input due to the need for compatibility with CDENPROP which is always compiled 
   !> using double precision.
   subroutine write_molden_dyson_orbitals(dyson_orbitals,idyson_orbital_irrep,itarget_irrep,itarget_spin,dyson_orbital_norms)
      use molden_mod
      use symmetry
      use mpi_mod, only: myrank
      use gto_routines, only: cms_gto_norm
      use utils, only: delete_file
      implicit none
      real(kind=wp), intent(in) :: dyson_orbitals(:,:,:)
      integer, intent(in) :: idyson_orbital_irrep(:),itarget_irrep(:),itarget_spin(:)
      real(kind=wp), intent(in) :: dyson_orbital_norms(:,:)

      integer :: no_orbitals, no_target_states, no_bound_states
      type(geometry_obj) :: geometry
      type(molecular_orbital_basis_obj) :: Dyson_orbital_basis
      type(orbital_data_obj), allocatable :: orbital_data(:)
      type(molden_input_obj) :: molden_file
      integer :: i, j, k, n_cf, ao, err, irr, number_of_cgto_shells
      character(len=line_len) :: molden_file_name, cnum, dyson_basis_set_file
      integer, allocatable :: sym_ind(:), n_dyson_orbitals(:)
      type(CGTO_shell_data_obj), allocatable :: CGTO_shells(:)
      real(kind=cfp), allocatable :: cf(:,:), normalization_factor(:)
      real(kind=cfp) :: rmat_radius, spin_cg, fac

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'The orbital basis set data has not been initialized.', 1, 1)
         end if

         no_orbitals = size(dyson_orbitals,1)
         no_target_states = size(dyson_orbitals,2)
         no_bound_states = size(dyson_orbitals,3)

         if (size(dyson_orbital_norms,1) /= no_target_states .or. &
             size(dyson_orbital_norms,2) /= no_bound_states) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'On input the array dyson_orbital_norms does not have the expected dimensions.', 2, 1)
         end if

         if (sum(num_orbs_sym(:)) /= no_orbitals) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'Number of orbitals in CDENPROP is different to the number of orbitals in the orbital basis set.', 3, 1)
         end if

         allocate(orbital_data(1:molecular_orbital_basis%no_irr),&
                  sym_ind(1:molecular_orbital_basis%no_irr),&
                  n_dyson_orbitals(1:molecular_orbital_basis%no_irr),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation 1 failed',err,1)

         !Count the number of Dyson orbitals in each symmetry (this replicates the double loop from below)
         n_dyson_orbitals = 0
         do k=1,no_bound_states
            do j=1,no_target_states
               !the Dyson orbital symmetries given by CDENPROP start with 0
               n_dyson_orbitals(idyson_orbital_irrep(j)+1) = n_dyson_orbitals(idyson_orbital_irrep(j)+1) + 1
            enddo !j
         enddo !k

         !Note that this must be amended once symmetry for the AOs is used!
         n_cf = molecular_orbital_basis%ao_basis%number_of_functions
         do i=1,molecular_orbital_basis%no_irr

            allocate(orbital_data(i) % coefficients(n_cf,n_dyson_orbitals(i)), &
                     orbital_data(i) % energy(n_dyson_orbitals(i)), &
                     orbital_data(i) % spin(n_dyson_orbitals(i)), &
                     orbital_data(i) % occup(n_dyson_orbitals(i)), stat = err)
            if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation 2 failed',err,1)

            orbital_data(i)%number_of_coefficients = n_cf
            orbital_data(i)%irr = i
            orbital_data(i)%number_of_functions = n_dyson_orbitals(i)
            orbital_data(i)%point_group = molecular_orbital_basis%pg
            orbital_data(i)%occup = 1.0_cfp
            orbital_data(i)%spin = 1
            orbital_data(i)%energy = 0.0_cfp
         enddo

         call molecular_orbital_basis%ao_basis%get_all_CGTO_shells(CGTO_shells,number_of_cgto_shells)

         !Calculate normalization factors for the continuum GTOs and multiply them in with the orbital coefficients.
         allocate(normalization_factor(molecular_orbital_basis%ao_basis%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation -1 failed.',err,1)
         normalization_factor = 1.0_cfp

         rmat_radius = options%a
         ao = 0
         do i=1,number_of_cgto_shells
            if (rmat_radius > 0.0_cfp .and. molecular_orbital_basis%ao_basis%shell_descriptor(3,i) .eq. 1) then !Normalize the continuum to the R-matrix sphere
               fac = cms_gto_norm (rmat_radius, &
                                   CGTO_shells(i) % l, &
                                   CGTO_shells(i) % number_of_primitives, &
                                   CGTO_shells(i) % exponents, &
                                   CGTO_shells(i) % contractions, &
                                   CGTO_shells(i) % norm, &
                                   CGTO_shells(i) % norms)
               write(stdout,'("Continuum normalization factor",i0,e25.15)') i,fac
               normalization_factor(ao+1:ao+2*CGTO_shells(i)%l+1) = fac
            elseif (molecular_orbital_basis%ao_basis%shell_descriptor(3,i) .eq. 1 .and. rmat_radius .le. 0.0_cfp) then
               write(stdout,'("Continuum functions in this shell will not be normalized to &
                                &the R-matrix sphere since rmat_radius .le. 0.0_cfp.")')
            endif
            ao = ao + 2*CGTO_shells(i)%l+1
         enddo !i

         call molecular_orbital_basis%get_orbital_coefficient_matrix(cf)
         !
         !TRANSFORM THE DYSON ORBITALS FROM THE MOLECULAR ORBITAL BASIS INTO THE ATOMIC ORBITAL BASIS:
         !
         spin_cg = 1.0_cfp !1.0_cfp/sqrt(2.0_cfp) !Spin related Clebsch-Gordan coefficient (see commit 879 to trunk:cdenprop_io).
         write(stdout,'("The value of the spin related Clebsch-Gordan coefficient is: ",e25.15)') spin_cg

         write(stdout,'("Dyson orbital index to target state mapping: bound state,index within sym.,symmetry,target &
            &state sym.,target state spin")')
         sym_ind(:) = 0
         do k=1,no_bound_states
            do j=1,no_target_states
               irr = idyson_orbital_irrep(j)+1 !the Dyson orbital and target state symmetries given by CDENPROP start with 0
               sym_ind(irr) = sym_ind(irr) + 1
               write(stdout,'("Mapping: ",i3,1X,i5,1X,i1,1X,i1,1X,i1)') k,sym_ind(irr),irr, itarget_irrep(j)+1, itarget_spin(j)
               orbital_data(irr)%coefficients(:,sym_ind(irr)) = 0.0_wp
               orbital_data(irr)%occup(sym_ind(irr)) = dyson_orbital_norms(j,k)

               !For this Dyson orbital accumulate the contributions to the Dyson orbital's AO coefficients from all MO orbitals.
               do i=1,no_orbitals
                  if (dyson_orbitals(i,j,k) .eq. 0.0_cfp) cycle
                  do ao=1,n_cf
                     orbital_data(irr)%coefficients(ao,sym_ind(irr)) = orbital_data(irr)%coefficients(ao,sym_ind(irr)) &
                                                                        + cf(ao,i)*dyson_orbitals(i,j,k)
                  enddo !ao
               enddo !i

               !Include normalization of the continuum functions (needed for Molden) and the spin-related CG coefficient.
               do ao=1,n_cf
                  orbital_data(irr)%coefficients(ao,sym_ind(irr)) = orbital_data(irr)%coefficients(ao,sym_ind(irr)) &
                                                                    *normalization_factor(ao)*spin_cg
               enddo !ao

            enddo !j
         enddo !k
         !
         !OUTPUT THE DYSONS TO DISK:
         !
         dyson_basis_set_file = './dyson_orbitals.ukrmolp'
         molden_file_name = './dyson_orbitals.molden'
         write(cnum,'(i0)') myrank
         molden_file_name = trim(molden_file_name)//trim(adjustl(cnum))
         dyson_basis_set_file = trim(dyson_basis_set_file)//trim(adjustl(cnum))

         !Get the geometry and symmetry information from the GTO basis set
         call molecular_orbital_basis%ao_basis%symmetry_data%get_geometry(geometry)

         if (molecular_orbital_basis%ao_basis%contains_btos()) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'The Molden file cannot be written since the basis contains BTOs which are not supported by Molden.', 3, 0)
         else
            !Delete any already existing file with the target file name
            call delete_file(molden_file_name)

            !Write the CGTO basis and all Dyson orbitals into the Molden file
            call molden_file%init(molden_file_name,2)
            call molden_file%write(geometry%nucleus,CGTO_shells,orbital_data)
            call molden_file%final()
         endif

         !Initialize the basis set of Dyson orbitals - this will store all Dyson orbital coefficients and the full GTO basis. We are doing this so that
         !we always keep all information about the Dyson orbitals (the Molden file may not be able to store information for all L in the basis).
         Dyson_orbital_basis%ao_basis => molecular_orbital_basis%ao_basis
         err = Dyson_orbital_basis%init(molecular_orbital_basis%no_irr,geometry)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'Error initializing the Dyson orbital basis data structure.', 3, 1)
         end if

         do irr=1,molecular_orbital_basis%no_irr
            !Remove the normalization factor for the continuum functions from the orbital coefficients: that's only needed for Molden since it does not normalize to the R-matrix radius.
            do ao=1,n_cf
               orbital_data(irr)%coefficients(ao,:) = orbital_data(irr)%coefficients(ao,:)/normalization_factor(ao)
            enddo !ao
            !Add the Dyson orbitals for each symmetry to the full set of Dyson orbitals
            call Dyson_orbital_basis%add_shell(orbital_data(irr))
         enddo

         !Delete any already existing file with the same name as dyson_basis_set_file
         call delete_file(dyson_basis_set_file)

         !Finally, save the full set of Dysons and the AO basis set to a separate file that can be used for later processing
         call molecular_orbital_basis%ao_basis%write(dyson_basis_set_file)
         call Dyson_orbital_basis%write(dyson_basis_set_file)

         err = Dyson_orbital_basis%final()
         if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Error finalizing Dyson_orbital_basis.',4,1)

   end subroutine write_molden_dyson_orbitals

!OUTER REGION (SWINTERF) INTERFACING ROUTINES:

   !> Calculates the orbital amplitudes and determines the channel information. The result is stored in the module private variables amplitudes,continuum_channels_m_l_irr for later use by UKP_PREAMP, UKP_READAMP.
   !> Call to this routine must be preceeded by calling e.g. READ_UKRMOLP_BASIS in order to initialize the molecular_orbital_basis object which performs the amplitude evaluation.
   !> \warning This routine ALWAYS requires double precision value of the R-matrix radius on input due to the need for compatibility with UKRmol-out which is always compiled using double precision.
   subroutine EVAL_AMPLITUDES(a,normalize_to_a)
      implicit none
      real(kind=wp), intent(in) :: a
      logical, intent(in) :: normalize_to_a
      real(kind=cfp) :: a_conv

      logical, allocatable :: is_continuum(:)
      integer :: err, irrcont

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'eval_amplitudes', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

!         if (a .eq. real(options%a,kind=wp)) then
!            normalize_to_a = .true.
!         else
!            normalize_to_a = .false.
!         endif

         a_conv = a
         call molecular_orbital_basis%calculate_amplitudes(a_conv,normalize_to_a,amplitudes,continuum_channels_m_l_irr)

         if (.not.allocated(ncont)) then
            allocate(ncont(molecular_orbital_basis%no_irr),stat=err)
            if (err .ne. 0) call xermsg('ukrmol_interface','EVAL_AMPLITUDES','Memory allocation failed.',err,1)
            ncont = 0
         endif

         do irrcont=1,molecular_orbital_basis%no_irr
            call molecular_orbital_basis%get_continuum_flags(irrcont,is_continuum)
            ncont(irrcont) = count(is_continuum)

            write(stdout,'(/,10X,"IRR = ",i2," number of continuum orbitals = ",i0)') irrcont, ncont(irrcont)
         enddo !irrcont

   end subroutine EVAL_AMPLITUDES

   !> Transfers the channel angular and IRR numbers to the channel variables from SWINTERF. It is the equivalent of the PREAMP routine from SWINTERF. The call to this routine must be preceeded by call to
   !> EVAL_AMPLITUDES.
   subroutine UKP_PREAMP(IRRCONT,IX,LCHL,MCHL,QCHL,NCH,MAXNMO)
      implicit none
      integer, intent(in) :: IRRCONT,IX
      integer, intent(out) :: LCHL(*),MCHL(*),QCHL(*),NCH
      integer, intent(inout) :: MAXNMO

      integer :: lm, nch_lm, irr

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'UKP_PREAMP', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         nch_lm = size(continuum_channels_m_l_irr,2)

         nch = 0
         do lm=1,nch_lm
            irr = continuum_channels_m_l_irr(3,lm)
            if (irr .eq. irrcont) then
               nch = nch + 1
               LCHL(IX+nch-1) = continuum_channels_m_l_irr(2,lm)
               MCHL(IX+nch-1) = abs(continuum_channels_m_l_irr(1,lm))
               QCHL(IX+nch-1) = sign(1,continuum_channels_m_l_irr(1,lm))
               if (MCHL(IX+nch-1) .eq. 0) QCHL(IX+nch-1) = 0
            endif
         enddo

         maxnmo = max(maxnmo,ncont(irrcont))

   end subroutine UKP_PREAMP

   !> This routine transfers into WAMPS the R-matrix amplitudes for the orbitals. It is equivalent of the READAMP routine. Note that the orbital amplitudes are divided by sqrt(2). 
   !> See the routine outerio/READRM routine and the equation for the R-matrix there. It implies that the 1/2 factor from the R-matrix formula is absorbed into the product of the two amplitudes 
   !> in the summation over the R-matrix poles. Hence the factor 1/sqrt(2) by which the amplitudes must be multiplied.
   !> The call to this routine must be preceeded by call to UKP_PREAMP.
   !> \warning This routine ALWAYS outputs double precision values of the amplitudes in WAMPS due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   subroutine UKP_READAMP(WAMPS,NCHAN,IRRCHL,LCHL,MCHL,QCHL,NCONTCSF,MCONT,IPRINT)
      implicit none
      integer :: IRRCHL(*), LCHL(*), MCHL(*), QCHL(*), NCONTCSF(8), MCONT(:), NCHAN, IPRINT
      real(kind=wp) :: WAMPS(NCHAN,*)

      integer :: i, l, lm, nch_lm, orbs_start, orbs_end, nfound, nunique, ntgtsym
      logical :: match

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'UKP_READAMP', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         if (iprint > 0) WRITE(stdout,1000)
         if (iprint > 1) WRITE(stdout,1010) NCHAN,(I,IRRCHL(I),I=1,NCHAN)

         nch_lm = size(continuum_channels_m_l_irr,2) !number of different (l,m) channel numbers
         ntgtsym = size(mcont)    !number of symmetries of the target states

         if (iprint > 0) then
            do i=1,molecular_orbital_basis%no_irr
               nunique = 0
               do lm=1,nch_lm
                  if (continuum_channels_m_l_irr(3,lm) .eq. i) nunique = nunique + 1
               enddo
               WRITE(stdout,2010) i,molecular_orbital_basis%ao_basis%number_of_functions,num_orbs_sym(i),nunique
            enddo
         endif

         !loop over all channels in the problem
         nfound = 0
         do L=1,NCHAN

            if (IRRCHL(L) > molecular_orbital_basis%no_irr .or. IRRCHL(L) <= 0) then
                call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                             'IRRCHL value is out of range and is incompatible with the orbital basis data.', 2, 1)
            end if

            !make sure that for each target state symmetry there really is a continuum channel with symmetry that matches the contiuum symmetry given by MCONT (from SCATCI).
            match = .false.
            do i=1,ntgtsym
               if (IRRCHL(L) .eq. mcont(i)+1) match = .true.
            enddo
            if (.not. match) then
                call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                             'Mismatch between orbital basis and CI vector data: &
                             &symmetry of at least one continuum channel lies outside of the set of continuum symmetries &
                             &given by the MCONT values obtained from CI (SCATCI) data.', 3, 1)
            end if

            !loop over all channels from the boundary amplitude data and look for the matching channel data
            do lm=1,nch_lm
               if (IRRCHL(L) == continuum_channels_m_l_irr(3,lm) .and. &
                   MCHL(L)*QCHL(L) == continuum_channels_m_l_irr(1,lm) .and. &
                   LCHL(L) == continuum_channels_m_l_irr(2,lm)) then
                  nfound = nfound + 1
                  orbs_end = sum(num_orbs_sym(1:IRRCHL(L)))
                  orbs_start = orbs_end-ncont(IRRCHL(L))+1 !index of the first continuum orbital in symmetry IRRCHL(L).
                  i = sum( num_orbs_sym( 1:(IRRCHL(L)-1) ) )+1 !index of the first orbital in symmetry IRRCHL(L).
                  if (orbs_start < i) then
                     print *,orbs_start,i
                     call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                                  'Wrong number of continuum orbitals or an error in orbital data.', 3, 1)
                  endif
                  WAMPS(L,1:ncont(IRRCHL(L))) = amplitudes(lm,orbs_start:orbs_end)/sqrt(2.0_cfp)
                  if (iprint > 0) then
                      write(stdout,2810) L, continuum_channels_m_l_irr(3,lm), continuum_channels_m_l_irr(2,lm), &
                                            continuum_channels_m_l_irr(1,lm), (WAMPS(L,I),I=1,ncont(IRRCHL(L)))
                  end if
                  exit
               endif
            enddo !lm

         enddo !L

         ncontcsf = 0
         ncontcsf(1:molecular_orbital_basis%no_irr) = ncont(1:molecular_orbital_basis%no_irr)

         if (nfound /= nchan) then
            call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                         'The number of channels for which boundary amplitudes have been transfered &
                         &does not match the actual number of channels for which the amplitudes are required.', 4, 1)
         end if

 1000 FORMAT(/,10X,'====> UKP_READAMP - READ BOUNDARY AMPS <====',/)
 1010 FORMAT(/,10X,'Searching for data on ',I3,' channels ',//,&
     &         10X,'Seq No.  IRRCHL',/,&
     &         10X,'-------  ------',//,&
     &        (10X,I6,2X,I4))
 2010 FORMAT(/,10X,'Symmetry number = ',I3,//,10X,'No. of basis functions = ',I3,/,&
               10X,'No. of orbitals        = ',I3,/,10X,'No. of angular parts   = ',I3,/)
 2810 FORMAT(/1X,'Chan No. ',I2,' Symmetry = ',I2,' Angular behaviour (l,m)',' = ',2I5,/,(1X,5(F12.8,1X)))

   end subroutine UKP_READAMP

end module ukrmol_interface
