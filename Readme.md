This document describes the steps to prepare the input files for a complete calculation of the charge migration in the CO2 molecule, subject to external pulses. This case study can be used as a reference for other systems as well, since each Project follows the same steps. 
The process is highly automatized, with python scripts generating most of the necessary input files for the molcas and chargemigration programs. However, some of the parameters must be modified by the user both at the beginning and during the process.

*** The commands in the present manual refer to the distribution in fermi. Other
    workstation may require further adjustments.
    
To run the scripts successfully, python3.6 (or later) must be installed, alongside recents version of pip3.6 (or later), seaborn, matplotlib, and any other module that may be necessary

  sudo pip3.6 install --upgrade pip
  sudo python3.6 -m pip install seaborn matplotlib natsort
 
  or 

  python3.6 -m pip install --user seaborn matplotlib natsort

you must also install OpenMolcas, with the support libraries hdf5 and hdf5-devel (in CentOS: dnf install hdf5 hdf5-devel), and optionally its external graphical interface Luscus 


1. The first step in the process is the definition of the molecular geometry and electronic space and the use of OpenMolcas to compute the necessary structural observables for the system. In the present case study, we will focus on the input files necessary for the CO2 molecule. However, OpenMolcas offers a wide range of methods to compute the same quantity with various degree of accuracy, which may be appropriate for other systems. The interested reader will find plenty of information in OpenMolcas documentation (https://molcas.gitlab.io/OpenMolcas/sphinx/index.html, https://molcas.gitlab.io/OpenMolcas/Manual.pdf for the pdf form).

1a. Molecular Geometry File
Every project begins with a file containing the cartesian coordinates for the molecule in study, which must be in .xyz format. This file can be made either using a text editor or one of the many molecular graphical interface programs such as Luscus or Gabedit. The geometry file must be placed in the working directory. It is a good practice to use for the root name of the geometry file the name or formula of the molecule under study (e.g. CO2.xyz).

The xyz file specifies the number of atoms in the molecules, the molecule name, and the list of atoms with their symbol and cartesian coordinates in atomic units. The example below for the CO2 molecule, is taken from https://pubchem.ncbi.nlm.nih.gov/compound/Carbon-dioxide, where more examples of molecules at their equilibrium position can be found.


    -------------------- CO2.xyz --------------------
    3									     #Number of Atoms
    Carbon Dioxide							     #Name of Molecule
    O         -1.19700        0.00000        0.00000	     #Atoms
    O          1.19700        0.00000        0.00000
    C          0.00000        0.00000        0.00000
    -------------------------------------------------


1b) To compute the energies and properties of the electronic states of the molecule, we need to specify the type of ab initio calculation Molcas should perform (e.g., HF, CASSCF, etc.), the orbital basis, etc. Finally, we need to run Molcas. The calculations are launched with a python script, which is based on the configuration file "chargemigration.ini".

If you have not already prepared the configuration file "chargemigration.ini", you can generate a template with default values by running the script main.py

$./main.py

or

$main.py

if main.py is in path, and modify it afterwards. 

Once you have created the template configuration file, you should edit the entries under the sections "Project settings" and "Grid settings" , both will be used during the run of Molcas. The configuration file contains also several other parameters, which will be needed only later to direct the Charge Migration code (e.g. Charge Migration settings and below).

You can get an idea of the SCF orbitals beforehand using MolCalc, a handy online application. Otherwise it might be usefull to run a SCF calculation and take advantage of Molcas GUI Luscus to select the active space to be used during the calculations. Of course this is not required, and active space can easily be specified in the input file needed for Molcas; to help the user a template can be prepared by running the script main.py with arguement -i?

$main.py -i?


    ----------------------inputhelp.input---------------------------
    &GATEWAY 
     Title = _name_of_project_ 
     coord = _geometry_file_ 
     basis = _name_of_basis_ 
     group = c1 //this prevents use of symmetry
    &SEWARD
    &SCF
    &RASSCF
     Ras1   = _number_of_RAS1_orbitals_
     Ras2   = _number_of_RAS2_orbitals_
     RAS3   = _number_of_RAS3_orbitals_
     Inactive = _number_of_Inactive_orbitals_
     nactel   = _total_ _RAS1_ _RAS3_ //# of active electrons, allowed holes, allowed excitations
     TDM       //writes DM and TDM to RASSCF.h5
     CIRoot
     _highest_root_needed_ _dimention_of_CI_matrix_ _weight_ //e.g. 10 10 1 will do a state average calculation of the ground state and the lower 4 excited states
    &RASSI
    Nr of JobIph
     1 ALL     //e.g. # of job files, # of states in that file ==>> 1 10 or 1 ALL 
    MEES //writes one-electron properties
    -----------------------------------------------------------------


This contains no real input values for Molcas, only keywords needed and comments, therefore the user can edit it as needed. However, keywords used to print values under RASSCF and RASSI sections shall not be altered since these values will be used by the Charge-Migration program. 

The input file must  also be changed and renamed. # MUST NOT HAVE THE NAME <PROJECT>.input


1c) Selecting Active Space using Luscus (optional, if active space is known then skip)

If the active space is to be determined using luscus, then the input file must be updated in the Project-settings section

    ___INSERT PICTURE OF FILE____
    chargemigration.ini
    [Project settings]
    project_name = CO2 
    xyz molecule geometry = CO2.xyz
    basis = ANO-L-MB

and the option -lus should be specified (e.g. main.py -lus). The program will run an SCF calculation and launch the LUSCUS GUI application which will directly allow the user to select the active orbitals and RAS spaces for the subsequent RASSCF calculation.

Click on the "orbitals" tab and order them according to their index. At this stage, all orbitals should be labelled as either inactive or secondary. To define the active space for a RAS2 calculation, you can change the labels of the orbitals you want to active by clicking on their label and selecting RAS2 from the small scroll menu that will appear.

After choosing your active space, go to the "File" tab, click "Save As", and choose the name to be the project name with file format .GvOrb, for example "CO2.GvOrb", then choose to save as "Molcas orbital data" on the bottom right corner. Make sure the saved file is present in the current directory, if so, then LUSCUS can be closed. If you find the need to make any last minute changes to the active space after this point, instead of of repeating the previous steps, simply call $luscus using the terminal, open your $Project.lus file (e.g. CO2.lus) and proceed to make your changes.  

Before continuing, open your configuration file and fill the necessary settings. Keep in mind that the number of states asked for will also represent the highest state computed. For example, for our CO2 calculations we ask for "5" states with 1 being the ground and the rest excited states. If you need to change the size of the grid or the number of points, please do so as well. Lastly to force Molcas to use the active space found inside <ProjectName>.GvOrb then under the RASSCF section of the input file add "fileorb= <ProjectGvOrb_File>.GvOrb" and "typeindex".

    e.g. &RASSCF
        fileorb= CO2.GvOrb
        typeindex
        nactel   = 10 0 1 //# of active electrons, allowed holes, allowed excitations
        TDM       //writes DM and TDM to RASSCF.h5
        CIRoot
        10 10 1 
 
As seen above, the sections for RAS1, RAS2, RAS3 are not needed when using typeindex.


1d) Now that the configuration file and Molcas's input have both been properly filled, we can finally run Molcas. 


$main.py -i <input_file>    # THE INPUT FILE MUST NOT HAVE THE NAME <PROJECT>.input 


This run can take a while, depending on the size of the active space and the number of active electrons. The script and molcas will notify possible input or execution errors. If no errors are found, once completed, the main.py script will generate the necessary data for the subsequent run of the charge-migration program. These data are stored in the directory <Project>_output.


    ----------------------CO2_output Directory----------------------
         
         From Molcas output logfile
            X_DIPOLE
            Y_DIPOLE
            Z_DIPOLE
            
         From CO2.grid inside molcasworkdir directory
            grid7
            grid8
            grid9
            grid10
            grid11
            grid12
            grid13
            
         From CO2.rasscf.h5 inside molcasworkdir directory
            DENSITY_MATRIX
            TRANSITION_DENSITY_MATRIX
            ROOT_ENERGIES
            AO_MLTPL_X
            AO_MLTPL_Y
            AO_MLTPL_Z
            MO_ENERGIES
            MO_VECTORS
            
    ----------------------------------------------------------------

!Running Charge Migraiton

After confirming that all the expected files are found inside the <ProjectName>_output director, fill in the Charge Migration and Pulses settings in the configuration file with the parameters needed (or provide a Field Pulse file using "-field"), the proceed to call main.py -cm [-sden, -field].



	 
	 
	 


#added this comment just to test if i (ruben) am using the correct branch