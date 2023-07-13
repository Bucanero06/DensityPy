# Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
# email: luca.argenti@gmail.com
# email: luca.argenti@ucf.edu
# Luca Argenti is Associate Professor of Physics, Optics and Photonics
# at the Department of Physics and the College of Optics
# of the University of Central Florida
# 4111 Libra Drive
# Orlando, Florida, USA
#
CC=gcc
FC=ifort
#FC=gfortran
#..
DEB_FLAG=
HOST=$(shell hostname)
MAKE=/usr/bin/make
HOMEDIR=$(shell pwd)
COM_SHARE=${HOMEDIR}/etc/CommonShare.mk
COM_ASTRA=${HOMEDIR}/etc/CommonAstra.mk
COM=${HOMEDIR}/etc/CommonRules.mk
SHARE=${HOMEDIR}/share
SHARE_ASTRA=${HOMEDIR}/share_astra
LIBDIR=${HOMEDIR}/${DEB_FLAG}lib_${FC}
BIN=${HOMEDIR}/../bin_${FC}
OBJ=$(DEB_FLAG)obj_$(FC)
MOD=$(DEB_FLAG)mod_$(FC)
INCLUDE=${HOMEDIR}/inc
DEFAULT_INT=
#DEFAULT_INT=-i8

export HOST MAKE HOMEDIR LIBDIR INCLUDE COM COM_SHARE COM_ASTRA SHARE SHARE_ASTRA CC FC DEB BIN OBJ MOD DEB_FLAG 

ifeq ($(FC),ifort)

	#  MKL_PATH=/opt/intel/mkl/lib/intel64
	#  MKL_INCLUDE=/opt/intel/mkl/include

	# Via Intel's oneAPI Math Kernel Library (oneMKL)
	MKL_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64
	MKL_INCLUDE=/opt/intel/oneapi/mkl/latest/include


  LINK_OPTS=-I$(MKL_PATH) -qmkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -I$(MKL_INCLUDE)
  MOD_OPTION=-I$(MKL_INCLUDE) -module
  ifeq ($(DEB_FLAG),d)
    F77_OPTS=-c -g -traceback -safe-cray-ptr ${DEFAULT_INT} 
    FC_OPTS=-c -p  ${DEFAULT_INT}  -warn all -check all -debug all -g -traceback -heap-arrays -safe-cray-ptr -I$(MKL_INCLUDE)
    CC_OPTS = -c 
    FC_OPTS_ATSP=-c -p  ${DEFAULT_INT}  -g -traceback -safe-cray-ptr
    #F77_OPTS=-c -g -traceback -safe-cray-ptr ${DEFAULT_INT} 
    #FC_OPTS=-c -p  ${DEFAULT_INT}  -warn all -check all -debug all -g -traceback -safe-cray-ptr -I$(MKL_INCLUDE)
    #CC_OPTS = -c 
    #FC_OPTS_ATSP=-c -p  ${DEFAULT_INT}  -g -traceback -safe-cray-ptr
  else
    F77_OPTS=-c -O2  ${DEFAULT_INT}  
    FC_OPTS=-c -O2  ${DEFAULT_INT} -heap-arrays -I$(MKL_INCLUDE)
    FC_OPTS_ATSP=-c -O2  ${DEFAULT_INT} 
    #F77_OPTS=-c -O2  ${DEFAULT_INT}  -heap-arrays  
    #FC_OPTS=-c -O2  ${DEFAULT_INT}   -heap-arrays  -I$(MKL_INCLUDE)
    #FC_OPTS_ATSP=-c -O2  ${DEFAULT_INT} -heap-arrays 
    CC_OPTS = -c -O3 
  endif
  FC_MALLOC=LINUX_64
  #FC_MALLOC=LINUX

else ifeq ($(FC),gfortran)

  MKL_PATH=/opt/intel/mkl/lib/intel64
  #/intel64
  LINK_OPTS=-fopenmp -L$(MKL_PATH) -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
  MOD_OPTION=-J
  ifeq ($(DEB_FLAG),d)
    F77_OPTS=-c -Wall -fcray-pointer
    FC_OPTS=-c -Wall -fcheck=all -fcray-pointer -fbacktrace -g  
    CC_OPTS = -c
    FC_OPTS_ATSP=-c -Wall -fcray-pointer
  else
    #F77_OPTS=-c -O2 -heap-arrays
    F77_OPTS=-c -O2 -fcray-pointer
    FC_OPTS=-c -O2 -fcray-pointer
    FC_OPTS_ATSP=-c -O2 -fcray-pointer
    CC_OPTS = -c -O2
  endif
  FC_MALLOC=LINUX_64

else

  echo "unrecognized compiler"

endif



export FC_OPTS F77_OPTS CC_OPTS LINK_OPTS MOD_OPTION MKL_PATH FC_OPTS_ATSP FC_MALLOC

LIBSRCDIR=${HOMEDIR}/libsrc
libinstall:
	mkdir -p $(BIN)

	for i in $(LIBSRCDIR) ; \
	do \
		cd $$i ; \
		echo "Entering: $$i" ; \
		$(MAKE) install ; \
		cd .. ; \
		echo "Leaving: $$i" ; \
		echo; \
		echo;\
	done
	# for i in $(shell echo $(LIBDIR)/*.a); \
	# do \
	# 	ar -x $$i ;\
	# done


#.. Program update rules
#..
#.PHONY : $(ALL_PRGS) TDSE_PETSc
.PHONY : $(ALL_PRGS) 
$(ALL_PRGS) : export PRG_NAME=$(@F)
$(ALL_PRGS) :
	$(MAKE) $(BIN)/$(DEB_FLAG)$(@) -C $(PRG_NAME)
	rm -f ../bin
	ln -s $(BIN) ../bin



#all : $(ALL_PRGS) TDSE_PETSc
all : libinstall $(ALL_PRGS)

TDSE_PETSc : export PRG_NAME=$(@F)
TDSE_PETSc :
	$(MAKE) $(@F) -C $(@F)

#.. Clean rules
#..
CLEAN=$(foreach a,$(ALL_PRGS),clean_$(a))
#
$(CLEAN):
	$(MAKE) clean -C $(subst clean_,,$@)
	rm -f ../bin_$(FC)/$(DEB_FLAG)$(subst clean_,,$@)
#
clean : $(CLEAN)
	$(MAKE) clean -C $(notdir $(SHARE))
	$(MAKE) clean -C $(notdir $(SHARE_ASTRA))
	$(MAKE) clean -C $(notdir $(LIBSRCDIR))
	rm -f $(LIBDIR)/*.a
	rm -f inc/*.mod

