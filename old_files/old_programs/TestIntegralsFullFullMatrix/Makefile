#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=GBTO
EXTRA_FC_OPTS=-i4

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
loc_ff90= main ModuleMainInterface ModuleParametersJUAN InitRoutines ModuleBasisJUAN GetIntegrals STEXMatrixElements ModuleDensityMatricesJUAN AngularFunctions 
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(SHARE)/$(OBJ)/ModuleString.o                   \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o             \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o                  \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o          \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleESpace.o             \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o                \
                    $(OBJ)/ModuleBasisJUAN.o                         \
                    $(OBJ)/ModuleParametersJUAN.o                    \
                    $(OBJ)/ModuleDensityMatricesJUAN.o               \
                    $(OBJ)/ModuleMainInterface.o                     \
                    $(OBJ)/InitRoutines.o                            \
                    $(OBJ)/STEXMatrixElements.o                      \
                    $(OBJ)/GetIntegrals.o                            \
                    $(OBJ)/AngularFunctions.o

$(OBJ)/ModuleMainInterface.o: \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleString.o                   \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o

$(OBJ)/ModuleBasisJUAN.o: \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o

$(OBJ)/InitRoutines.o: \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o            \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o             \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o                  \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o              \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o             \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o                \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleUKRmolInterface.o    \
                    $(OBJ)/ModuleBasisJUAN.o                         \
                    $(OBJ)/ModuleParametersJUAN.o

$(OBJ)/GetIntegrals.o: \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleUKRmolInterface.o    \
                    $(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o          \
                    $(SHARE)/$(OBJ)/ModuleMatrix.o                   \
                    $(SHARE)/$(OBJ)/ModuleConstants.o                \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o              \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o                  \
                    $(OBJ)/ModuleParametersJUAN.o                    \
                    $(OBJ)/ModuleBasisJUAN.o 

$(OBJ)/ModuleDensityMatricesJUAN.o: \
                    $(OBJ)/ModuleBasisJUAN.o

$(OBJ)/AngularFunctions.o: \
                    $(SHARE)/$(OBJ)/ModuleAngularMomentum.o

$(OBJ)/STEXMatrixElements.o: \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o                  \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o              \
                    $(OBJ)/ModuleBasisJUAN.o                         \
                    $(OBJ)/ModuleDensityMatricesJUAN.o               \
                    $(OBJ)/ModuleParametersJUAN.o 
 
#.. Commands with common syntax
#..
include $(COM)

