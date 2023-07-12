#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=GBTO
#EXTRA_FC_OPTS=-i8

#.. General files
#.. 
include $(COM_SHARE)
gen_ff77 = $(com_gen_ff77)
gen_ff90 = $(com_gen_ff90)

#.. Astra shared files
#..
include $(COM_ASTRA)
gen_astra_ff90 = $(com_gen_astra_ff90)

#.. Local files
#..
loc_ff90 =  main ModuleMainInterface
loc_fc   =

#.. Local files
#..
loc_ff90 =  main ModuleMainInterface
loc_fc   = 

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(OBJ)/ModuleMainInterface.o \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleMatrix.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleParentIons.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraConfigFile.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleDensityMatrices.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleESpace.o \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleSymESpace.o \

$(OBJ)/ModuleMainInterface.o: \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraCredit.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleString.o



#.. Commands with common syntax
#..
include $(COM) 
