#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=cint

#.. General files
#.. 
gen_ff77 = ipsort dpsort dquadpack xermsg 

gen_ff90 = \
ModuleErrorHandling \
ModuleSystemUtils \
ModuleString \
ModuleIO \
ModuleCommandLineParameterList \
ModuleConstants \
ModuleQCBasis \
ModuleMolecularGeometry \
ModuleDaltonLibCint \
ModuleLibCint

#.. Local files
#..
loc_ff90=TabulatePotential 
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/TabulatePotential.o: \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(SHARE)/$(OBJ)/ModuleQCBasis.o \
                    $(SHARE)/$(OBJ)/ModuleMolecularGeometry.o \
                    $(SHARE)/$(OBJ)/ModuleLibCint.o \
                    $(SHARE)/$(OBJ)/ModuleDaltonLibCint.o 
 

#.. Commands with common syntax
#..
include $(COM)

