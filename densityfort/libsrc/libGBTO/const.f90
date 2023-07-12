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
module const
   use precisn
   use iso_fortran_env, only: error_unit, input_unit, output_unit, character_storage_size

   !> unit for standard output.
   !> This value is set by the routine redirect_stdout_to_file but it can be done only once during the course of the program.
   !> Typically, the value of stdout will be chosen different for each rank (in case of MPI runs) so that each rank will use
   !> its own unit for std. output. This also ensures that no modifications in the code are needed concerning standard output.
   integer, protected :: stdout = output_unit !, parameter :: stdout = output_unit
   !> unit for standard input
   integer, parameter :: stdin  = input_unit
   !> unit for standard error
   integer, parameter :: stderr = error_unit
   !> The size in bits of a character storage unit.
   integer, parameter :: char_len = character_storage_size

   !> Leading part of the file-name that will be used to collect standard output for each process. Each process gets its own file
   !> for standard output, e.g. "log_file.2" will contain standard output from MPI process with rank = 2.
   character(len=*), parameter :: stdout_file_name = "log_file."

   !> The base number that will stand for the unit number for each processes's standard output (see mpi_mod). This number must be
   !> chosen large enough so it cannot collide with any other file units that may be open during the program run.
   integer, parameter :: stdout_unit_base = 10000

   !> Parameter used in redirect_stdout_to_file to compute the unit number for each processes's standard output:
   !>     stdout = stdout_unit_base+myrank*stdout_unit_step
   !> This parameter allows to choose e.g. master's stdout to be unit 6 while the other processes will not interfere with any other
   !> fort units, e.g.: stdout_unit_base=6, stdout_unit_step=1000.
   integer, parameter :: stdout_unit_step = 1

   !> length of one line
   integer, parameter :: line_len = 207 !132

   integer, parameter :: len_ufmat = 11
   character(len_ufmat-2), parameter :: fmat = 'formatted'
   character(len_ufmat), parameter :: ufmat = 'unformatted'

   !> Default character string used for the int_input_output%name variable.
   character(len=19), parameter :: no_header = 'No header specified'

   !> How many different types of 1-electron integrals have been implemented so far: Overlap (1), Kinetic energy (2),
   !> Kinetic energy+Nuclear attraction (3), Nuclear attraction (4), property integrals (5).
   integer, parameter :: number_of_types_el_ints = 5

   !> Character parameters used by the method l_mol_basis%molecular_integrals to identify the type of molecular integral requested.
   !> These headers are used by the user to request the specific integral.
   character(len=*), parameter :: one_electron_ints = 'One electron integrals'
   character(len=*), parameter :: overlap_ints = 'Overlap integrals'
   character(len=*), parameter :: kinetic_ints = 'Kinetic energy integrals'
   character(len=*), parameter :: property_ints = 'Property integrals'
   character(len=*), parameter :: nuc_rep_att_ints = 'Nuclear attraction integrals'
   character(len=*), parameter :: one_elham = 'One electron Hamiltonian integrals'
   character(len=*), parameter :: two_el_ints = 'Two-electron integrals'
   character(len=*), parameter :: one_p_sym_ints = 'One-particle AO->MO transformed integrals for symmetric operators'
   character(len=*), parameter :: two_p_sym_ints = 'Two-particle AO->MO transformed integrals for symmetric operators'
   character(len=*), parameter :: ijkl_indices_header = 'ijkl indices'

   !> Character string present on the first line of a file which identifies the file as readable by the data_file_obj object.
   !> The length of this string must be line_len since that length is used in the data_file module when reading
   !> the header from the file.
   character(line_len), parameter :: data_file_obj_id = "DATA FILE version 1.0"

   !> Number of molecular orbitals for which orbital coefficients will be printed on one line by the routine
   !> molecular_orbital_basis_obj%print_orbitals.
   integer, parameter :: orbs_line = 6

   !> Maximum allowed number of contraction coefficients defining the contracted GTO function. This parameter can be
   !> increased/decreased as needed. It is used throught the program to define dimensions of some arrays.
   integer, parameter :: max_contr_len = 20

   !> The coefficients below are coefficients in the fit: mmax(T) = -18.206_cfp + 0.42734_cfp*T + 0.00075795_cfp*T**2
   !> This fit was obtained using Mathematica and estimates, for each T, the maximum value of m for which Fm(T) can be calculated
   !> using the ASYMPTOTIC formula to full double precision accuracy rel. error .le. 10^(-15)). Fm(T) is the Boys function.
   !> The values are stored in single precision since that is enough to calculate mmax(T) with high enough accuracy sufficient
   !> for the estimate.
   integer, parameter :: fit_order = 2
   real(kind=cfp), parameter :: wp_fit_terms(1:fit_order+1) = (/-18.206_cfp,0.42734_cfp,0.00075795_cfp/)

   !> The same as above (ASYMPTOTIC formula limit) but this time for full quadruple precision of Fm(T).
   real(kind=cfp), parameter :: ep_fit_terms(1:fit_order+1) = (/-21.408_cfp,0.19460_cfp,0.00093852_cfp/)

   !>The same as above but this time for the UPWARD recursion to full double precision accuracy.
   real(kind=cfp), parameter :: wp_fit_terms_up(1:fit_order+1) = (/-18.561_cfp,0.3748_cfp,0.00088246_cfp/)

   !> The same as above (UPWARD recursion formula limit) but this time for full quadruple precision of Fm(T).
   real(kind=cfp), parameter :: ep_fit_terms_up(1:fit_order+1) = (/-19.331_cfp,0.14211_cfp,0.00106170_cfp/)

   !> Convergence parameter for the calculation of the Boys function using boys_function.
   !> The convergence is chosen as the machine epsilon for real(kind=cfp) reals.
   real(kind=cfp), parameter :: boys_tol = epsilon(1.0_cfp)

   !> Convergence parameter for the calculation of the Boys function using boys_function.
   !> We need to define the quad precision separately due to the routine boys_function_quad which is explicitly written
   !> in quad precision.
   real(kind=ep1), parameter :: boys_tol_qprec = epsilon(1.0_ep1)

   !> Default value of Tmax which is used in init_boys to precalculate Fm(T) on the grid of T and m values. T=60 corresponds to
   !> the value for which, in double precision, the Boys function is calculated using the series expansion to full double precision
   !> for mmax=24. The corresponding value of T for quad precision is T=140. This parameter can be varied at will - it has only
   !> a minor effect on efficiency of the Taylor-based evaluation of Fm(T).
   real(kind=cfp), parameter :: boys_f_dprec_asym_thr = 60.0_cfp

   !> Default value of mmax for which is used in cgto_hgp to precalculate Fm(T) on the grid of T and m values. The minimum value
   !> of mmax is 4*L, where L is the maximum GTO L used in the calculations. The rest of the Boys function parameters have been
   !> set-up to give fully accurate results for m up to mmax=100. The value 24 below is just a sensible value that will be enough
   !> for most calculations. If the integral algorithm fails complaining that mmax is too small then set it here to the value at
   !> least 4*L where is the maximum L in your GTO atomic basis.
   integer, parameter :: mmax = 24

   !> Step in T in the grid of values of Fm(T) calculated for Taylor-expansion based evaluation of the Boys function
   !> for double precision calculations.
   real(kind=cfp), parameter :: boys_f_grid_step_wp = 0.1_cfp

   !> Step in T in the grid of values of Fm(T) calculated for Taylor-expansion based evaluation of the Boys function
   !> for quad precision calculations.
   real(kind=cfp), parameter :: boys_f_grid_step_ep = 0.01_cfp

   !> Empirically found value for the largest power of T required to calculate the Boys function \f$F_{m}(T)\f$ for
   !> \f$T=0,\dots,60\f$ and \f$m=0,\dots,24\f$ with the accuracy epsilon(1.0_wp).
   integer, parameter :: imax_wp = 140

   !> The same as imax_wp but for quad precision result.
   integer, parameter :: imax_ep = 340

   !> Order of the Taylor expansion for Taylor-expansion-based evaluation of the Boys function Fm(T) for double precision calculations.
   integer, parameter :: taylor_k_wp = 8

   !> Order of the Taylor expansion for Taylor-expansion-based evaluation of the Boys function Fm(T) for quad precision calculations.
   integer, parameter :: taylor_k_ep = 10

   !> Molecular integrals with relative precision last than or equal to this value will trigger an error. This parameter is used
   !> only as a default in integral_options%prec, so this parameter can be effectively adjusted on run-time. 
   !> \warning Not all integral calculation routines are necessarily using this parameter.
   real(kind=cfp), parameter :: int_rel_prec = 10**(-(precision(cfp_dummy)-5.0_cfp)) !10E-10_cfp for double precision

   !> Molecular integrals (contracted) smaller than this value will be neglected. Similarly to int_rel_prec this value is used
   !> only as a sensible default in integral_options%tol.
   real(kind=cfp), parameter :: int_del_thr = 10**(-(precision(cfp_dummy)-4.0_cfp)) !10E-11_cfp for double precision

   !> Maximum allowed value for a cross overlap of two orbitals after the symmetric orthogonalization has been performed.
   !> \warning This value should be actually equal to the threshold value for the integrals.
   real(kind=cfp), parameter :: thrs_symm_ortho = 10**(-(precision(cfp_dummy)-9.0_cfp)) !10E-6_cfp for double precision

   !> Maximum allowed value for a cross overlap of two orbitals after the Gramm-Schmidt orthogonalization has been performed.
   real(kind=cfp), parameter :: thrs_gs_ortho = 10**(-(precision(cfp_dummy)-5.0_cfp)) !10E-10_cfp for double precision

   !> Coefficients in the transformation matrix for the symmetric orthogonalization smaller than this value will be neglected.
   real(kind=cfp), parameter :: thrs_cf_sym_ortho_trans = 10**(-(precision(cfp_dummy)-5.0_cfp)) !10E-10_cfp for double precision

   !> Threshold value for self-overlap of an orbital (before normalization) for the Gramm-Schmidt orthogonalization. 
   !> If an orbital has a self-overlap (before normalization) smaller than this value then we assume that linear dependency in
   !> the orbital basis is present. In other words we decide that self-overlaps
   !> smaller than this value would lead to numerical problems.
   real(kind=cfp), parameter :: thrs_lin_dep_gs_ortho = 10**(-(precision(cfp_dummy)-8.0_cfp)) !1.0E-07_cfp for double precision

   !> Threshold for magnitude of final orbital coefficients as used by molecular_orbital_basis_obj%delete_small_coefficients.
   real(kind=cfp), parameter :: thrs_orb_cf = 10**(-(precision(cfp_dummy)-3.0_cfp))

   !> Absolute precision for the numerical quadrature routine dqags.
   real(kind=cfp), parameter :: EPSABS = 10**(-(precision(cfp_dummy)-5.0_cfp)) !10E-10_cfp for double precision
   !> Relative precision for the numerical quadrature routine dqags.
   real(kind=cfp), parameter :: EPSREL = 10**(-(precision(cfp_dummy)-5.0_cfp)) !10E-10_cfp for double precision
   !> Determines the maximum number of subintervals in the partition of the given integration interval in dqags.
   integer, parameter :: LIMIT = 1000
   !> Dimensioning parameter for dqags.
   integer, parameter :: LENW = 4*LIMIT

   !> Dimensioning parameter in DQELG which controls the maximum length of the vector that can be extrapolated.
   !> The original SLATEC routine used 52.
   integer, parameter :: max_epstab = 152

   !> Multiplication table for Abelian point groups.
   integer, parameter :: abel_prod_tab(8,8) &
     & =RESHAPE( (/        &
     & 1,2,3,4,5,6,7,8,    &
     & 2,1,4,3,6,5,8,7,    &
     & 3,4,1,2,7,8,5,6,    &
     & 4,3,2,1,8,7,6,5,    &
     & 5,6,7,8,1,2,3,4,    &
     & 6,5,8,7,2,1,4,3,    &
     & 7,8,5,6,3,4,1,2,    &
     & 8,7,6,5,4,3,2,1 /), (/8,8/) )

   !> Numerical identifier of the C1 point group-symmetry
   integer, parameter :: C1_id = 1
   !> Numerical identifier of the Cs point group-symmetry
   integer, parameter :: Cs_id = 2
   !> Numerical identifier of the C2 point group-symmetry
   integer, parameter :: C2_id = 3
   !> Numerical identifier of the Ci point group-symmetry
   integer, parameter :: Ci_id = 4
   !> Numerical identifier of the C2v point group-symmetry
   integer, parameter :: C2v_id = 5
   !> Numerical identifier of the C2h point group-symmetry
   integer, parameter :: C2h_id = 6
   !> Numerical identifier of the D2 point group-symmetry
   integer, parameter :: D2_id = 7
   !> Numerical identifier of the D2h point group-symmetry
   integer, parameter :: D2h_id = 8

   !> Length of the character variable 'name' in the nucleus_type object.
   integer, parameter :: nuc_nam_len = 2

   !> Name of the scattering centre.
   character(len=nuc_nam_len) :: nam_scattering_centre = 'sc'

   !> ID assigned to the scattering centre. Do not change this value.
   integer, parameter :: id_scattering_centre = 0

   !> Length of the character variable sym_op specifying the symmetry operation in the geometry_obj object.
   integer, parameter :: sym_op_nam_len = 3

   !> Names and order of IRR for C2v point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: C2v_names(4) = (/'A1','B1','B2','A2'/)
   !> s_v, s_vp symmetry elements. Atkins, Tables for Group theory.
   integer, parameter :: C2v_char_tab(4,2) &
     & =RESHAPE( (/        &
     & 1, 1,-1,-1,&
     & 1,-1, 1,-1 /), (/4,2/) )

   !> Names and order of IRR for D2h point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: D2h_names(8) = (/'Ag ','B3u','B2u','B1g','B1u','B2g','B3g','Au '/)
   !> s_3, s_2, s_1 symmetry elements. Atkins, Tables for Group theory.
   integer, parameter :: D2h_char_tab(8,3) &
     & =RESHAPE( (/        &
     & 1, 1, 1, 1,-1,-1,-1,-1, &
     & 1, 1,-1,-1, 1, 1,-1,-1, &
     & 1,-1, 1,-1, 1,-1, 1,-1 /), (/8,3/) )

   !> Names and order of IRR for C2h point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: C2h_names(4) = (/'Ag','Au','Bu','Bg'/)
   !> c2, c_i, s_h symmetry elements. Atkins, Tables for Group theory.
   integer, parameter :: C2h_char_tab(4,3) &
     & =RESHAPE( (/        &
     & 1, 1,-1,-1, &
     & 1,-1,-1, 1, &
     & 1,-1, 1,-1 /), (/4,3/) )

   !> Names and order of IRR for D2 point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: D2_names(4) = (/'A ','B3','B2','B1'/)
   !> c2_2, c2_1 symmetry elements. Atkins, Tables for Group theory.
   integer, parameter :: D2_char_tab(4,2) &
     & =RESHAPE( (/        &
     & 1,-1, 1,-1, &
     & 1, 1,-1,-1 /), (/4,2/) )

   !> Names and order of IRR for C2 point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: C2_names(2) = (/'A','B'/)
   !> c2 symmetry element. Atkins, Tables for Group theory.
   integer, parameter :: C2_char_tab(2,1) &
     & =RESHAPE( (/        &
     & 1,-1 /), (/2,1/) )

   !> Names and order of IRR for Cs point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: Cs_names(2) = (/'Ap ','App'/)
   !> s_h symmetry element. Atkins, Tables for Group theory.
   integer, parameter :: Cs_char_tab(2,1) &
     & =RESHAPE( (/        &
     & 1,-1 /), (/2,1/) )

   !> Names and order of IRR for Ci point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: Ci_names(2) = (/'Ag','Au'/)
   !> c_i symmetry element. Atkins, Tables for Group theory.
   integer, parameter :: Ci_char_tab(2,1) &
     & =RESHAPE( (/        &
     & 1,-1 /), (/2,1/) )

   !> IRR for Ci point group. Molpro order.
   character(len=sym_op_nam_len), parameter :: C1_names(1) = (/'A'/)
   !> Idetity symmetry element. Atkins, Tables for Group theory.
   integer, parameter :: C1_char_tab(1,1) &
     & =RESHAPE( (/        &
     & 1 /), (/1,1/) )

   !> Combined irreducible representation labels used in Psi4 Molden files.
   character(len=sym_op_nam_len), parameter :: pg_irr_names(8,8) = RESHAPE &
     ((/'A  ','   ','   ','   ','   ','   ','   ','   ',                   &  ! (:,1) ~ C1
        "A' ",'A" ','   ','   ','   ','   ','   ','   ',                   &  ! (:,2) ~ Cs
        'A  ','B  ','   ','   ','   ','   ','   ','   ',                   &  ! (:,3) ~ C2
        'Ag ','Au ','   ','   ','   ','   ','   ','   ',                   &  ! (:,4) ~ Ci
        'A1 ','B1 ','B2 ','A2 ','   ','   ','   ','   ',                   &  ! (:,5) ~ C2v
        'Ag ','Au ','Bu ','Bg ','   ','   ','   ','   ',                   &  ! (:,6) ~ C2h
        'A  ','B3 ','B2 ','B1 ','   ','   ','   ','   ',                   &  ! (:,7) ~ D2
        'Ag ','B3u','B2u','B1g','B1u','B2g','B3g','Au '/), (/8,8/))           ! (:,8) ~ D2h

   !> If it is requested that the Y_lm functions (evaluated by eval_BTO_CGTO_Y_lm) are saved to disk then all Y_lm with size (in MiB) greater than this value will be saved to disk.
   !> If you want to save all Y_lm regardless of their size then set this parameter to 0.0_cfp.
   real(kind=cfp), parameter :: Y_lm_size_threshold = 100.0_cfp

   !> Number of bytes in one Mib.
   integer, parameter :: Mib = 1048576

   !> Cache line size (bytes).
   integer, parameter :: cache_line_size = 64

   !> Tile size used for cache blocking. This is derived from cacheline size.
   !> \todo define l1_cache_size and use it to manage sizes of the buffers: check if the source-target buffers fit inside the cache.
   !>       If not then reallocate them (if possible) to their minimal sizes.
   integer, parameter :: tile = 64

   !> This value is controlled by the routine redirect_stdout_to_file.
   logical, private :: stdout_set = .false.

   !> Masses of the elements in the periodic table. Taken from DENPROP.
   !> \warning The masses of some elements (e.g. Boron) are way off from the NIST values!
   real(kind=cfp), parameter ::   amass(103) = (/ &
          1.0078246_cfp,    4.002601_cfp,     7.01600_cfp,      9.01218_cfp,     11.009307_cfp,     &
         12.000000_cfp,    14.0030738_cfp,   15.9949141_cfp,   18.9984022_cfp,   19.992441_cfp,     &
         22.9898_cfp,      23.98504_cfp,     26.98153_cfp,     27.976929_cfp,    30.973764_cfp,     &
         31.9720727_cfp,   34.9688531_cfp,   39.962386_cfp,    38.96371_cfp,     39.96259_cfp,      &
         44.95592_cfp,     48._cfp,          50.9440_cfp,      51.9405_cfp,      54.9380_cfp,       &
         55.9349_cfp,      58.9332_cfp,      57.9353_cfp,      62.9296_cfp,      63.9291_cfp,       &
         68.9257_cfp,      73.9219_cfp,      74.9216_cfp,      79.9165_cfp,      78.91839_cfp,      &
         83.91151_cfp,     84.9117_cfp,      87.9056_cfp,      88.9059_cfp,      89.9043_cfp,       &
         92.9060_cfp,      97.9055_cfp,      98._cfp,         101.9037_cfp,     102.9048_cfp,       &
        107.90389_cfp,    106.90509_cfp,    113.9036_cfp,     114.9041_cfp,     120._cfp,           &
        120.9038_cfp,     129.9067_cfp,     126.90466_cfp,    131.90416_cfp,    132.9051_cfp,       &
        137.9050_cfp,     138.9061_cfp,     139.9053_cfp,     140.9074_cfp,     141.9075_cfp,       &
        145._cfp,         151.9195_cfp,     152.9209_cfp,     157.9241_cfp,     159.9250_cfp,       &
        163.9288_cfp,     164.9303_cfp,     165.9304_cfp,     168.9344_cfp,     173.9390_cfp,       &
        174.9409_cfp,     179.9468_cfp,     180.9480_cfp,     183.9510_cfp,     186.9560_cfp,       &
        192._cfp,         192.9633_cfp,     194.9648_cfp,     196.9666_cfp,     201.970625_cfp,     &
        204.9745_cfp,     207.9766_cfp,     208.9804_cfp,     209._cfp,         210._cfp,           &
        222._cfp,         223._cfp,         226._cfp,         227._cfp,         232._cfp,           &
        231._cfp,         238._cfp,         237._cfp,         244._cfp,         243._cfp,           &
        247._cfp,         247._cfp,         251._cfp,         252._cfp,         257._cfp,           &
        258._cfp,         259._cfp,         260._cfp                                                &
   /)

contains

  !> This routine can be called only once. It redirects standard output into the file with name "stdout_file_name"//myrank
  !> where // stands for string concatenation and myrank is integer. The idea is that the mpi_mod module routine mpi_start
  !> redirects the standard output for each process into its own file by calling this routine with the rank of the process
  !> as argument. See also description of the parameters stdout_unit_base, stdout_unit_step.
  subroutine redirect_stdout_to_file(myrank,master_stdout)
     implicit none
     integer, intent(in) :: myrank
     character(len=line_len) :: file_name, num
     logical,optional,intent(in)        ::      master_stdout
     
     logical                            ::      master_stdout_
    
     master_stdout_ = .false.
     if(present(master_stdout)) master_stdout_ = master_stdout

        if (stdout_set) then
           stop "const/redirect_stdout_to_file: stdout has been set before: stdout can be set only once."
        else
           if (stdout < output_unit) then
              stop "const/redirect_stdout_to_file: &
                   &on input val < output_unit, where output_unit is the default unit for std. output."
           end if

           !We open a file for each process separately where standard output for each process will go.
           write(num,'(i0)') myrank
           file_name = trim(stdout_file_name)//trim(adjustl(num))
           
           if(master_stdout_ .and. myrank == 0) then
                stdout = 6
           else
                stdout = stdout_unit_base+myrank*stdout_unit_step
                open(unit=stdout, file=file_name, status="replace")
           endif
           !The following commented line can be used to ignore the standard output completely (on unix systems).
           !This may be useful for calculations which use a large number of MPI processes. In that case
           !I may want to keep output of the master process only and igonore stdout of the rest.
           !open(unit=stdout, file="/dev/null", status="old")

           stdout_set = .true.
        endif
     
  end subroutine redirect_stdout_to_file

end module
