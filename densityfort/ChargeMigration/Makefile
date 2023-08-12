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
ModuleBSpline \
ModuleParameterList \
ModulePulses_3D \
ncdphi \
ModuleDiagonalize 

#.. Local files
#..
#loc_ff90=main_old_trying_to_reconcile ModuleRTP Module_CD_IO Module_Becke #Module_Liouvillian_Propagator
loc_ff90=main ModuleRTP Module_CD_IO Module_Becke #Module_Liouvillian_Propagator
loc_fc  =

#.. Local dependencies
#..
#$(OBJ)/main_old_trying_to_reconcile.o: \

$(OBJ)/main.o: \
                    $(SHARE)/$(OBJ)/ModuleErrorHandling.o \
                    $(SHARE)/$(OBJ)/ModuleSystemUtils.o \
                    $(SHARE)/$(OBJ)/ModuleString.o \
                    $(SHARE)/$(OBJ)/ModuleIO.o \
                    $(SHARE)/$(OBJ)/ModuleBSpline.o \
                    $(SHARE)/$(OBJ)/ModuleParameterList.o \
                    $(SHARE)/$(OBJ)/ModuleConstants.o \
                    $(SHARE)/$(OBJ)/ModulePulses_3D.o \
                    $(SHARE)/$(OBJ)/ModuleDiagonalize.o \
                    $(OBJ)/ModuleRTP.o \
                    $(OBJ)/Module_CD_IO.o \
                    $(OBJ)/Module_Becke.o
                    #$(OBJ)/Module_Liouvillian_Propagator.o

$(OBJ)/ModuleRTP.o:

$(OBJ)/Module_CD_IO.o:

$(OBJ)/Module_Becke.o:

#$(OBJ)/Module_Liouvillian_Propagator.o:


#.. Commands with common syntax
#..
include $(COM)

