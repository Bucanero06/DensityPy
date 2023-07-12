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
EXTLIBS=

#.. General files
#.. 
gen_ff77 = ipsort dpsort dquadpack xermsg 

gen_ff90 = \
ModuleCommandLineParameterList \
ModuleConstants \
ModuleErrorHandling \
ModuleIO \
ModuleParameterList \
ModulePOSIX \
ModuleString \
ModuleSystemUtils

gen_astra_ff90 = \
ModuleAstraConfigFile \
ModuleAstraCredit \
ModuleAstraEnv \
ModuleCfgTemplates 


#.. Local files
#..
loc_ff90=main ModuleMainInterface
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(OBJ)/ModuleMainInterface.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraConfigFile.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleCfgTemplates.o \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModulePOSIX.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o 
 
$(OBJ)/ModuleMainInterface.o: \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraConfigFile.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraEnv.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraCredit.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o 

#.. Commands with common syntax
#..
include $(COM)

