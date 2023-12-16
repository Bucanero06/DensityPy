#! /usr/bin/env python
## THIS FILE SHOULD BE CONSIDERED DEPRECATED
# >import selectionofactivespace as sas
from configparser import ConfigParser
from os import path
import sys

from densitypy.project_utils.command_execution import execute_command, print_molcas_log_errors
from densitypy.project_utils.file_directory_ops import find


def SelectionOfActiveSpace(ini_file):
    """
    (DEPRECATED) This function is used to select the active space using Luscus GUI (DEPRECATED)

    :param ini_file:
    :return:
    """
    # Read Project settings from the configuration file
    parser = ConfigParser()

    if ini_file:
        parser.read(ini_file)
    else:
        parser.read('chargemigration.ini')

    project_name = parser.get('Project settings', 'project_name')
    XYZ_Geometry = parser.get('Project settings', 'XYZ Molecule Geometry')
    basis_definition = parser.get('Project settings', 'Basis')

    # Write/Prepare input file used to select Active Space using luscus
    with open(project_name + '.input', 'w') as pyinput:
        pyinput.write("&GATEWAY"
                      " \n Title = " + project_name +
                      " \n coord = " + XYZ_Geometry +
                      "; basis = " + basis_definition +
                      "; group = c1"
                      "\n&SEWARD"
                      "\n&SCF"
                      "\n&GRID_IT\nALL")

    # Run Pymolcas using the input file previously created
    print("Wait while Molcas runs a SCF calculations and creates a Luscus file you can use to analyze the Orbitals")
    # execute_pymolcas_with_error_print("pymolcas " + project_name + ".input -f", project_name)
    try:
        execute_command("pymolcas " + project_name + ".input -f", write_stream=True)
    except Exception as e:
        print(e)
        print_molcas_log_errors(project_name + ".log", "Timing")
    # Look for <ProjectName>.lus file to be used by LUSCUS GUI
    pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")
    # molcas_workdir = GetValueOfAsString(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
    molcas_workdir = "MOLCAS_WORKDIR"
    luscusfiledirectory = find(project_name + ".lus", ".", molcas_workdir)

    # Call luscus to allow selection of Active Space
    print("Opening Molcas GUI Luscus")

    execute_command(f"luscus {luscusfiledirectory}/{project_name}.lus", write_stream=True)

    # Remove unnessesary created files just created
    execute_command("rm " + project_name + ".input ", write_stream=True)

    # Simply a reminder to the user to check input file
    # this of course could be done automatically,
    # however was not done here and the process is non time consuming
    print("Active space selected, please update the input file accordingly before proceeding."
          "\n To use the " + project_name + ".GvOrb file containning the active space that was "
          "just created using Luscus, add fileorb= " + project_name + ".GvOrb to the input file "
          "under Module RASSCF as well as include KeyWord 'typeindex'\n RASSCF will begin from "
          "these Orbital file, thus keywords 'RAS(1-3Zafr) =' should be removed to avoid "
          "confilics while declaring the active space")
    sys.exit()
