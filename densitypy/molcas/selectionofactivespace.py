#! /usr/bin/env python
## THIS FILE SHOULD BE CONSIDERED DEPRECATED
# >import selectionofactivespace as sas
from configparser import ConfigParser
from os import path
import sys

from densitypy.molcas.molcasscripts import copy_and_prepare_molcas_input_file_for_run, call_open_molcas
from densitypy.project_utils.command_execution import execute_command, print_molcas_log_errors
from densitypy.project_utils.file_directory_ops import find, make_directory, copy_to


def SelectionOfActiveSpace(json_config, **kwargs):
    """
    (DEPRECATED) This function is used to select the active space using Luscus GUI (DEPRECATED)

    :param ini_file:
    :return:
    """
    # Read Project settings from the configuration file
    parser = ConfigParser()

    # Split JSON config into sections
    project_settings = json_config['projectsettings']

    # Project Settings
    project_name = project_settings['projectname']
    xyz_geometry_path = project_settings['xyzmoleculegeometry']
    list_of_orbitals = project_settings['listofactiveorbitals']
    molcas_output_directory = project_settings['molcasoutputdirectory']

    # If basis set is not defined, use 6-31G* or if is passed in kwargs, use that
    if not kwargs.get('basis', False):
        basis = '6-31G*'

    # Prepare Input
    make_directory(molcas_output_directory)
    copy_to(xyz_geometry_path, molcas_output_directory)
    temp_project_name = f'selecting_space_{project_name}'
    molcas_input_path = path.join(molcas_output_directory, f'{temp_project_name}.input')


    # Write/Prepare input file used to select Active Space using luscus
    with open(molcas_input_path, "w") as pyinput:
        pyinput.write("&GATEWAY"
                      " \n Title = " + project_name +
                      " \n coord = " + xyz_geometry_path +
                      "; basis = " + basis +
                      "; group = c1"
                      "\n&SEWARD"
                      "\n&SCF"
                      # "\n&GRID_IT\nALL"
                      )

    # Run Pymolcas using the input file previously created
    print("Wait while Molcas runs a SCF calculations and creates the orbital files which you can use to analyze the Orbitals")
    try:
        call_open_molcas(temp_project_name, molcas_output_directory)
    except Exception as e:
        print(e)
        log_file_path = path.join(molcas_output_directory, f'selecting_space_{project_name}.log')
        if path.exists(log_file_path):
            try:
                print_molcas_log_errors(log_file_path, "Timing")
            except Exception as e:
                print(e)
        exit()

    scf_h5_path = f"{find(f'selecting_space_{project_name}.scf.h5',molcas_output_directory)}/selecting_space_{project_name}.scf.h5"

    if not scf_h5_path:
        print("SCF failed, check log file for errors")
        exit()

    # todo provide helpful insight into what orbitals are recommended prior to even loading pegamoid

    # Find what is the absolute path of pegamoid.py so we can call it from the terminal
    pegamoid_path = path.abspath(__file__).split('selectionofactivespace.py')[0] + 'pegamoid.py'
    print(f'pegamoid_path: {pegamoid_path}')
    if not pegamoid_path:
        print(f"pegamoid.py not found in {pegamoid_path}, please check if it is in the correct location")
        exit()
    # Call luscus to allow selection of Active Space
    print("Opening Molcas GUI Pegamoid")
    execute_command(f"python {pegamoid_path} {scf_h5_path}")

    sys.exit()
