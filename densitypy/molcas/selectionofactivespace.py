#! /usr/bin/env python
# >import selectionofactivespace as sas
from configparser import ConfigParser
from os import path
import sys

from densitypy.project_utils.def_functions import ExecutePymolcasWithErrorPrint, ExecuteNoWrite, Find


def SelectionOfActiveSpace(ini_file):
    # Read Project settings from the configuration file
    parser = ConfigParser()

    if ini_file:
        parser.read(ini_file)
    else:
        parser.read('chargemigration.ini')

    Project_Name = parser.get('Project settings', 'Project_Name')
    XYZ_Geometry = parser.get('Project settings', 'XYZ Molecule Geometry')
    basis_definition = parser.get('Project settings', 'Basis')

    # Write/Prepare input file used to select Active Space using luscus
    with open(Project_Name + '.input', 'w') as pyinput:
        pyinput.write("&GATEWAY"
                      " \n Title = " + Project_Name +
                      " \n coord = " + XYZ_Geometry +
                      "; basis = " + basis_definition +
                      "; group = c1"
                      "\n&SEWARD"
                      "\n&SCF"
                      "\n&GRID_IT\nALL")

    # Run Pymolcas using the input file previously created
    print("Wait while Molcas runs a SCF calculations and creates a Luscus file you can use to analyze the Orbitals")
    ExecutePymolcasWithErrorPrint("pymolcas " + Project_Name + ".input -f", Project_Name)

    # Look for <ProjectName>.lus file to be used by LUSCUS GUI
    pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")
    # molcas_workdir = GetValueOfAsString(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
    molcas_workdir = "MOLCAS_WORKDIR"
    luscusfiledirectory = Find(Project_Name + ".lus", ".", molcas_workdir)

    # Call luscus to allow selection of Active Space
    print("Opening Molcas GUI Luscus")

    ExecuteNoWrite("luscus " + luscusfiledirectory + "/" + Project_Name + ".lus")

    # Remove unnessesary created files just created
    ExecuteNoWrite("rm "
                   + Project_Name + ".input ")

    # Simply a reminder to the user to check input file
    # this of course could be done automatically,
    # however was not done here and the process is non time consuming
    print("Active space selected, please update the input file accordingly before proceeding."
          "\n To use the " + Project_Name + ".GvOrb file containning the active space that was "
          "just created using Luscus, add fileorb= " + Project_Name + ".GvOrb to the input file "
          "under Module RASSCF as well as include KeyWord 'typeindex'\n RASSCF will begin from "
          "these Orbital file, thus keywords 'RAS(1-3Zafr) =' should be removed to avoid "
          "confilics while declaring the active space")
    sys.exit()
