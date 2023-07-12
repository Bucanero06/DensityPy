# Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
# email: luca.argenti@gmail.com
# email: luca.argenti@ucf.edu
# Luca Argenti is Associate Professor of Physics, Optics and Photonics
# at the Department of Physics and the College of Optics
# of the University of Central Florida
# 4111 Libra Drive
# Orlando, Florida, USA
#
com_gen_ff77 = ipsort dpsort dquadpack xerror dgamma 

com_gen_ff90_level0 = \
ModuleErrorHandling   \
ModuleSystemUtils     \
ModulePOSIX           \
ModuleString          \
ModuleCommandLineParameterList \
ModuleAngularMomentum \
ModuleConstants       

com_gen_ff90_level1 = \
ModuleBSpline         \
ModuleMatrix          \
ModuleParameterList   \
ModuleIO              \
XCHEM_ModuleGroups

com_gen_ff90 = $(com_gen_ff90_level0) $(com_gen_ff90_level1)


$(SHARE)/$(OBJ)/ModuleBSpline.o : \
	                     $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
	                     $(SHARE)/$(OBJ)/ModuleParameterList.o \
	                     $(SHARE)/$(OBJ)/ModuleSystemUtils.o   \
	                     $(SHARE)/$(OBJ)/ModuleString.o 
$(SHARE)/$(OBJ)/ModuleParameterList.o : \
	                     $(SHARE)/$(OBJ)/ModuleErrorHandling.o
$(SHARE)/$(OBJ)/specfun.o : \
	                     $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
	                     $(SHARE)/$(OBJ)/cnormena1.o
$(SHARE)/$(OBJ)/ModuleDiagonalize.o  : \
	                     $(SHARE)/$(OBJ)/ModuleErrorHandling.o 
$(SHARE)/$(OBJ)/ModuleIO.o  : \
	                     $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                             $(SHARE)/$(OBJ)/ModuleString.o
$(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o : \
                             $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                             $(SHARE)/$(OBJ)/ModuleString.o
$(SHARE)/$(OBJ)/ModuleMatrix.o : \
                             $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                             $(SHARE)/$(OBJ)/ModuleString.o

