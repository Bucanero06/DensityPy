#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=

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
ModuleDiagonalize 

#.. Local files
#..
loc_ff90=TabulateOrbitals
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/TabulateOrbitals.o: \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(SHARE)/$(OBJ)/ModuleDiagonalize.o 
 

#.. Commands with common syntax
#..
include $(COM)

