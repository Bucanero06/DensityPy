# Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
# email: luca.argenti@gmail.com
# email: luca.argenti@ucf.edu
# Luca Argenti is Associate Professor of Physics, Optics and Photonics
# at the Department of Physics and the College of Optics
# of the University of Central Florida
# 4111 Libra Drive
# Orlando, Florida, USA
#
#.. Compile shared files
#..
gen_of77=$(foreach a,$(gen_ff77),$(SHARE)/$(OBJ)/$(a).o)
gen_of90=$(foreach a,$(gen_ff90),$(SHARE)/$(OBJ)/$(a).o)
gen_astra_of90=$(foreach a,$(gen_astra_ff90),$(SHARE_ASTRA)/$(OBJ)/$(a).o)
$(gen_of77) : $(SHARE)/$(OBJ)/%.o : $(SHARE)/src/%.f
	$(FC) $(F77_OPTS) -extend-source $(MOD_OPTION) $(SHARE)/$(MOD) $< -o $@
$(gen_of90) : $(SHARE)/$(OBJ)/%.o : $(SHARE)/src/%.f90
	$(FC) $(FC_OPTS) $(MOD_OPTION) $(SHARE)/$(MOD) -I$(INCLUDE) $(SHARE)/$(MOD) $< -o $@
$(gen_astra_of90) : $(SHARE_ASTRA)/$(OBJ)/%.o : $(SHARE_ASTRA)/src/%.f90
	$(FC) $(FC_OPTS) $(MOD_OPTION) $(SHARE_ASTRA)/$(MOD) -I$(INCLUDE) -I$(SHARE)/$(MOD) $< -o $@


#.. Compile local files
#..
loc_of77=$(foreach a,$(loc_ff77),$(OBJ)/$(a).o)
loc_of90=$(foreach a,$(loc_ff90),$(OBJ)/$(a).o)
loc_of90=$(foreach a,$(loc_ff90),$(OBJ)/$(a).o)
$(loc_of77) : $(OBJ)/%.o : src/%.f
	$(FC) $(F77_OPTS) -extend-source  $(MOD_OPTION) $(MOD) -I$(SHARE)/$(MOD) -I$(SHARE_ASTRA)/$(MOD) $< -o $@
$(loc_of90) : $(OBJ)/%.o : src/%.f90
	$(FC) $(FC_OPTS) $(EXTRA_FC_OPTS) $(MOD_OPTION) $(MOD) -I$(SHARE)/$(MOD) -I$(SHARE_ASTRA)/$(MOD) -I$(INCLUDE) $< -o $@

loc_oc  =$(foreach a,$(loc_fc),$(OBJ)/$(a).o)
$(loc_oc) : $(OBJ)/%.o : src/%.c
	$(CC) $(CC_OPTS) $< -o $@

#.. Link
#..
LIBLIST := $(foreach lib, $(EXTLIBS),-l$(lib))
$(BIN)/$(DEB_FLAG)$(PRG_NAME) : $(loc_of90) $(loc_of77) $(gen_of77) $(gen_of90) $(gen_astra_of90) $(loc_oc)
	$(FC) $(LINK_OPTS) -o $(BIN)/$(DEB_FLAG)$(PRG_NAME) $(loc_of90) $(loc_of77) $(gen_of77) $(gen_of90) $(gen_astra_of90) $(loc_oc) -L$(LIBDIR) $(LIBLIST) 

#.. Clean
#..
clean :
	rm -f $(OBJ)/*.o $(MOD)/*.mod



