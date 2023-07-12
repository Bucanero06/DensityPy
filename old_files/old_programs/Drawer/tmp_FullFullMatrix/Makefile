#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=GBTO
EXTRA_FC_OPTS=-i4

#.. General files
#.. 
gen_ff77 = ipsort dpsort dquadpack xermsg 

gen_ff90 = \
ModuleBSpline \
ModuleErrorHandling \
ModuleSystemUtils \
ModuleString \
ModuleIO \
ModuleCommandLineParameterList \
ModuleParameterList \
ModuleConstants \
ModuleIntegrals \
ModuleMolecularGeometry \
ModuleAngularMomentum \
ModuleDiagonalize \
ModuleUKRmolInterface \
ModuleMatrix \
XCHEM_ModuleGroups \
XCHEM_ModuleSymmetryAdaptedSphericalHarmonics \
symba


#ModuleDiagonalize  
#.. Local files
#..
loc_ff90= main Subroutines Functions ModuleParametersJUAN ModuleBasisJUAN
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(SHARE)/$(OBJ)/ModuleIntegrals.o \
                    $(SHARE)/$(OBJ)/ModuleAngularMomentum.o \
                    $(SHARE)/$(OBJ)/ModuleDiagonalize.o \
                    $(OBJ)/ModuleDensityMatricesJUAN.o  \
                    $(SHARE)/$(OBJ)/ModuleUKRmolInterface.o  \
                    $(SHARE)/$(OBJ)/ModuleMatrix.o \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o \
                    $(SHARE)/$(OBJ)/XCHEM_ModuleSymmetryAdaptedSphericalHarmonics.o \
                    $(SHARE)/$(OBJ)/symba.o \
                    $(SHARE)/$(OBJ)/dpsort.o  \
                    $(OBJ)/ModuleParametersJUAN.o \
                    $(OBJ)/ModuleBasisJUAN.o \
                    $(OBJ)/Subroutines.o \
                    $(OBJ)/Functions.o


$(OBJ)/Subroutines.o: \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o \
                    $(SHARE)/$(OBJ)/ModuleAngularMomentum.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(OBJ)/ModuleParametersJUAN.o \
                    $(OBJ)/ModuleBasisJUAN.o \
                    $(SHARE)/$(OBJ)/ModuleUKRmolInterface.o  \

$(OBJ)/Functions.o: \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o \
                    $(SHARE)/$(OBJ)/ModuleAngularMomentum.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(OBJ)/ModuleParametersJUAN.o \
                    $(OBJ)/ModuleBasisJUAN.o \
                    $(SHARE)/$(OBJ)/ModuleUKRmolInterface.o  \


$(OBJ)/ModuleParametersJUAN.o: 

$(OBJ)/ModuleBasisJUAN.o: \
                    $(SHARE)/$(OBJ)/ModuleUKRmolInterface.o

$(OBJ)/ModuleDensityMatricesJUAN.o: \
                    $(OBJ)/ModuleBasisJUAN.o


#                    $(SHARE)/$(OBJ)/ModuleDiagonalize.o  
#.. Commands with common syntax
#..
include $(COM)

