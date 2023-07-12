# Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
# email: luca.argenti@gmail.com
# email: luca.argenti@ucf.edu
# Luca Argenti is Associate Professor of Physics, Optics and Photonics
# at the Department of Physics and the College of Optics
# of the University of Central Florida
# 4111 Libra Drive
# Orlando, Florida, USA
#
com_gen_astra_ff90 =        \
ModuleESpace                \
ModuleSymESpace             \
ModuleCloseCouplingChannels \
ModuleDensityMatrices       \
ModuleIntegrals             \
ModuleParentIons            \
ModuleAstraCredit           \
ModuleAstraConfigFile       \
ModuleAstraEnv              \
ModuleMolecularGeometry     \
ModuleOrbitalBasis          \
ModuleXlm 


#-----------------------------------------------------------------------
#.. LEVEL 5
$(SHARE_ASTRA)/$(OBJ)/ModuleESpace.o :                               \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleSymESpace.o  \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleParentIons.o \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleMatrix.o           \
                            $(SHARE)/$(OBJ)/ModuleString.o

#-----------------------------------------------------------------------
#.. LEVEL 4
$(SHARE_ASTRA)/$(OBJ)/ModuleSymESpace.o :                            \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleCloseCouplingChannels.o \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleDensityMatrices.o \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleParentIons.o \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o  \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleOrbitalBasis.o  \
                            $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleAngularMomentum.o  \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleIO.o               \
                            $(SHARE)/$(OBJ)/ModuleString.o           \
                            $(SHARE)/$(OBJ)/ModuleParameterList.o    \
                            $(SHARE)/$(OBJ)/ModuleMatrix.o           \
                            $(SHARE)/$(OBJ)/ModuleConstants.o        

#-----------------------------------------------------------------------
#.. LEVEL 3
$(SHARE_ASTRA)/$(OBJ)/ModuleCloseCouplingChannels.o :                \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o        \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleParentIons.o \
                            $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleString.o           \
                            $(SHARE)/$(OBJ)/ModuleConstants.o        \
                            $(SHARE)/$(OBJ)/ModuleMatrix.o           \
                            $(SHARE)/$(OBJ)/ModuleIO.o
$(SHARE_ASTRA)/$(OBJ)/ModuleDensityMatrices.o  :                     \
                             $(SHARE_ASTRA)/$(OBJ)/ModuleParentIons.o\
                             $(SHARE_ASTRA)/$(OBJ)/ModuleMolecularGeometry.o \
                             $(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o \
	                     $(SHARE)/$(OBJ)/ModuleMatrix.o          \
                             $(SHARE)/$(OBJ)/ModuleAngularMomentum.o \
		             $(SHARE)/$(OBJ)/ModuleString.o          \
                             $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o 
$(SHARE_ASTRA)/$(OBJ)/ModuleUKRmolInterface.o :                      \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o  \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleOrbitalBasis.o \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleMolecularGeometry.o \
                            $(SHARE)/$(OBJ)/ModuleBSpline.o          \
                            $(SHARE)/$(OBJ)/ModuleConstants.o        \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleSystemUtils.o      \
                            $(SHARE)/$(OBJ)/ModuleAngularMomentum.o  

#-----------------------------------------------------------------------
#.. LEVEL 2
$(SHARE_ASTRA)/$(OBJ)/ModuleParentIons.o :                           \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o        \
                            $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleString.o           \
                            $(SHARE)/$(OBJ)/ModuleIO.o               \
                            $(SHARE)/$(OBJ)/ModuleParameterList.o
$(SHARE_ASTRA)/$(OBJ)/ModuleAstraCredit.o  :                         \
	                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraEnv.o
$(SHARE_ASTRA)/$(OBJ)/ModuleCfgTemplates.o  :                        \
	                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraConfigFile.o \
	                    $(SHARE_ASTRA)/$(OBJ)/ModuleAstraEnv.o 
$(SHARE_ASTRA)/$(OBJ)/ModuleIntegrals.o : \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleOrbitalBasis.o \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleMolecularGeometry.o \
                            $(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o        \
                            $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleString.o           \
                            $(SHARE)/$(OBJ)/ModuleMatrix.o           \
                            $(SHARE)/$(OBJ)/ModuleBSpline.o          \
                            $(SHARE)/$(OBJ)/ModuleSystemUtils.o      \
                            $(SHARE)/$(OBJ)/ModuleAngularMomentum.o  \
                            $(SHARE)/$(OBJ)/ModuleIO.o 

#-----------------------------------------------------------------------
#.. LEVEL 1 - Files depending only on general shared modules (LEVEL 0)
$(SHARE_ASTRA)/$(OBJ)/ModuleXlm.o  :                                 \
                            $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleString.o           \
                            $(SHARE)/$(OBJ)/ModuleConstants.o        \
                            $(SHARE)/$(OBJ)/ModuleAngularMomentum.o
$(SHARE_ASTRA)/$(OBJ)/ModuleMolecularGeometry.o :                    \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleSystemUtils.o      \
                            $(SHARE)/$(OBJ)/ModuleString.o           \
                            $(SHARE)/$(OBJ)/ModuleIO.o
$(SHARE_ASTRA)/$(OBJ)/ModuleAstraConfigFile.o  :                     \
	                    $(SHARE)/$(OBJ)/ModuleParameterList.o
$(SHARE_ASTRA)/$(OBJ)/ModuleOrbitalBasis.o  :                        \
	                    $(SHARE)/$(OBJ)/XCHEM_ModuleGroups.o     \
                            $(SHARE)/$(OBJ)/ModuleErrorHandling.o    \
                            $(SHARE)/$(OBJ)/ModuleBSpline.o          \
			    $(SHARE)/$(OBJ)/ModuleString.o
$(SHARE_ASTRA)/$(OBJ)/ModuleAstraEnv.o  :                            \
	                    $(SHARE)/$(OBJ)/ModulePOSIX.o

#======================================================================
