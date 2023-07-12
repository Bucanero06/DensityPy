#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=
#EXTRA_FC_OPTS=-i8

#.. General files
#.. 
#   Must replace symba with ModuleAngularMomentum
#..
gen_ff77 = dpsort dquadpack xerror dgamma 
gen_ff90 = \
ModuleSystemUtils \
ModuleErrorHandling \
ModuleString \
ModuleCommandLineParameterList \
ModuleConstants \
ModuleAngularMomentum 

#.. Local files
#..
loc_ff90 = main ModuleMainInterface 
loc_fc   = 

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleAngularMomentum.o \
                    $(OBJ)/ModuleMainInterface.o 

$(OBJ)/ModuleMainInterface.o: \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleString.o 

#.. Commands with common syntax
#..
include $(COM) 
