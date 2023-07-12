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
ModuleMatrix

#.. Local files
#..
loc_ff90=main ModuleMainInterface
loc_fc  =

#.. Local dependencies
#.. 
$(OBJ)/main.o: \
                    $(OBJ)/ModuleMainInterface.o \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleMatrix.o \
 
$(OBJ)/ModuleMainInterface.o: \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleCommandLineParameterList.o 

#.. Commands with common syntax
#..
include $(COM)

