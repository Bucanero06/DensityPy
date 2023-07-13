
ifeq ($(FC),ifort)
   FC_OPTS_GBTO=-O2 -c -i8
   #LINK_OPTS_GBTO=-mkl -i8 -Dusempi -Dsplitreduce -Dmpi64bitinteger
   LINK_OPTS_GBTO=-qmkl
   MOD_OPTS_GBTO=-module 
else ifeq($(FC),gfortran) 
   FC_OPTS_GBTO=-O2 -c -fdefault-integer-8 -ffree-line-length-132 
   LINK_OPTS_GBTO=
   MOD_OPTS_GBTO=-J
else
   echo "Unrecognized compiler"
endif

MOD_DIR_GBTO=./
OBJ_DIR_GBTO=./
LIB_DIR_GBTO=./

clean :
	rm *.o *.mod 

GBTO_lib_f90=atomic_basis_mod  common_obj  integral_storage_mod  pco_mod  basis_data_generic_mod  const lebedev  phys_const  bspline_base  coupling_obj  quadrature_module  bspline_grid_mod  data_file  molden_const  sort  bto_gto_integrals_mod   eri_sph_coord  molden_mod  special_functions  bto_integrals_mod  free_scattering_mod  molecular_basis_mod  symmetry  cgto_hgp  function_integration  orbital_routines  ukrmol_interface  cgto_integrals  general_quadrature  orthogonalization  utils  cgto_pw_expansions_mod  gto_routines  parallel_arrays  wigner_cf 

GBTO_lib_F90=blas_lapack  mpi_memory_mod  mpi_mod  precisn

GBTO_lib_f90o=$(foreach a,$(GBTO_lib_f90),$(a).o)
GBTO_lib_F90o=$(foreach a,$(GBTO_lib_F90),$(a).o)

$(GBTO_lib_f90o) : %.o : %.f90
	$(FC) $(FC_OPTS_GBTO) $(LINK_OPTS_GBTO) $< -o $@  $(MOD_OPTS_GBTO) $(MOD_DIR_GBTO) -I $(MOD_DIR_GBTO) 

$(GBTO_lib_F90o) : %.o : %.F90
	$(FC) $(FC_OPTS_GBTO) $(LINK_OPTS_GBTO) $< -o $@  $(MOD_OPTS_GBTO) $(MOD_DIR_GBTO) -I $(MOD_DIR_GBTO) 

install : $(GBTO_lib_f90o) $(GBTO_lib_F90o)
	ar -rcs libGBTO.a $^
	ranlib libGBTO.a 

const.o : precisn.o

ukrmol_interface.o : const.o mpi_mod.o precisn.o utils.o atomic_basis_mod.o molecular_basis_mod.o parallel_arrays.o integral_storage_mod.o mpi_memory_mod.o molden_mod.o

molden_mod.o : molden_const.o

atomic_basis_mod.o : precisn.o symmetry.o integral_storage_mod.o mpi_mod.o basis_data_generic_mod.o cgto_integrals.o bto_integrals_mod.o bto_gto_integrals_mod.o

bto_gto_integrals_mod.o : general_quadrature.o basis_data_generic_mod.o bspline_grid_mod.o common_obj.o coupling_obj.o cgto_pw_expansions_mod.o gto_routines.o const.o precisn.o utils.o

cgto_pw_expansions_mod.o : lebedev.o

cgto_integrals.o : cgto_hgp.o eri_sph_coord.o

free_scattering_mod.o : molecular_basis_mod.o blas_lapack.o

molecular_basis_mod.o : orthogonalization.o

orthogonalization.o : blas_lapack.o

symmetry.o : utils.o const.o common_obj.o mpi_mod.o special_functions.o

mpi_mod.o : precisn.o utils.o const.o

utils.o : const.o

common_obj.o : sort.o phys_const.o

integral_storage_mod.o : data_file.o parallel_arrays.o

parallel_arrays.o : mpi_memory_mod.o

basis_data_generic_mod.o : bspline_grid_mod.o gto_routines.o orbital_routines.o

bspline_grid_mod.o : bspline_base.o function_integration.o

function_integration.o : quadrature_module.o

quadrature_module.o : general_quadrature.o

gto_routines.o : coupling_obj.o

coupling_obj.o : wigner_cf.o
