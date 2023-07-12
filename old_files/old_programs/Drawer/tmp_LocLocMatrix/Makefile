#.. List of External libraries that are needed.
#   if a library depends on others, it must come first
EXTLIBS=GBTO
EXTRA_FC_OPTS=-i8

#.. General files
#.. 
gen_ff77 = ipsort dpsort dquadpack xermsg 

gen_ff90 = \
ModuleBSpline2 \
ModuleErrorHandling \
ModuleSystemUtils \
ModuleString \
ModuleIO \
ModuleCommandLineParameterList \
ModuleParameterList \
ModuleConstants \
ModuleAngularMomentum \
ModuleDiagonalize 



#.. Local files
#..
loc_ff90=main
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(SHARE)/$(OBJ)/ModuleBSpline2.o \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(SHARE)/$(OBJ)/ModuleAngularMomentum.o \
                    $(SHARE)/$(OBJ)/ModuleDiagonalize.o 
 

#.. Commands with common syntax
#..
include $(COM)

