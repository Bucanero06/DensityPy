# Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
# email: luca.argenti@gmail.com
# email: luca.argenti@ucf.edu
# Luca Argenti is Associate Professor of Physics, Optics and Photonics
# at the Department of Physics and the College of Optics
# of the University of Central Florida
# 4111 Libra Drive
# Orlando, Florida, USA
#
#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=GBTO
#EXTRA_FC_OPTS=-i8

#.. General Shared Files
#.. 
include $(COM_SHARE)
gen_ff77 = $(com_gen_ff77)
gen_ff90 = $(com_gen_ff90) 

#.. Astra Shared Files
include $(COM_ASTRA)
gen_astra_ff90 = $(com_gen_astra_ff90)\
ModuleUKRmolInterface

#.. Local files
#..
loc_ff90 =  main ModuleMainInterface
loc_fc   = 

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(OBJ)/ModuleMainInterface.o                     \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o            \
                    $(SHARE)/$(OBJ)/ModuleString.o                   \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o             \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleESpace.o             \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleSymESpace.o          \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraConfigFile.o    \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleUKRmolInterface.o    \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o

$(OBJ)/ModuleMainInterface.o: \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraCredit.o       \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o\
                    $(SHARE)/$(OBJ)/ModuleParameterList.o           \
                    $(SHARE)/$(OBJ)/ModuleString.o 

#.. Commands with common syntax
#..
include $(COM) 
