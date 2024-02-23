#! /usr/bin/env python
## THIS FILE SHOULD BE CONSIDERED DEPRECATED
# >import selectionofactivespace as sas
from configparser import ConfigParser
from os import path
import sys

from densitypy.molcas.molcasscripts import copy_and_prepare_molcas_input_file_for_run, call_open_molcas
from densitypy.molcas.pegamoid import run_pegamoid
from densitypy.project_utils.command_execution import execute_command, print_molcas_log_errors
from densitypy.project_utils.file_directory_ops import find, make_directory, copy_to

def SelectionOfActiveSpace(xyz_file: str):
    """Do ground state cas"""

    xyz_file=path.abspath(xyz_file)
    # create a molecule
    from densitypy.autocas.scine_autocas.autocas_utils.molecule import Molecule
    molecule = Molecule(xyz_file)

    # initialize autoCAS and Molcas interface
    from densitypy.autocas.scine_autocas import Autocas
    autocas = Autocas(molecule)
    from densitypy.autocas.scine_autocas.interfaces.molcas import Molcas
    molcas = Molcas([molecule])
    # make initial active space and evaluate initial DMRG calculation
    occ_initial, index_initial = autocas.make_initial_active_space()
    print(f"occ_initial: {occ_initial}")
    print(f"index_initial: {index_initial}")

    # setup interface
    molcas.project_name = "example"
    # molcas.settings.work_dir =
    # molcas.environment.molcas_scratch_dir =
    molcas.settings.xyz_file = xyz_file



    # cas and hyphen do not matter for method names
    molcas.settings.method = "DMRGCI"

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 250
    molcas.settings.dmrg_sweeps = 5

    # make initial active space and evaluate initial DMRG calculation
    occ_initial, index_initial = autocas.make_initial_active_space()

    # no input means HF calculation
    molcas.calculate()

    # do cas calculation
    cas_results = molcas.calculate(occ_initial, index_initial)

    # energy = cas_results[0]
    s1_entropy = cas_results[1]
    # s2_entropy = cas_results[2]
    mut_inf = cas_results[3]

    # plot entanglement diagram
    from densitypy.autocas.scine_autocas.plots.entanglement_plot import EntanglementPlot
    plot = EntanglementPlot()
    plt = plot.plot(s1_entropy, mut_inf)  # type: ignore
    plt.savefig(molcas.settings.work_dir + "/entang.pdf")  # type: ignore

    # make active space based on single orbital entropies
    cas_occ, cas_index = autocas.get_active_space(
        occ_initial, s1_entropy   # type: ignore
    )

    # cas and hyphen do not matter for method names
    molcas.settings.method = "dmrg-scf"

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 2000
    molcas.settings.dmrg_sweeps = 20

    # Do a calculation with this CAS
    final_energy, final_s1, final_s2, final_mut_inf = molcas.calculate(cas_occ, cas_index)

    # use results
    n_electrons = sum(cas_occ)
    n_orbitals = len(cas_occ)
    print(f"final energy:      {final_energy}")
    print(f"final CAS(e, o):  ({n_electrons}, {n_orbitals})")
    print(f"final cas indices: {cas_index}")
    print(f"final occupation:  {cas_occ}")
    print(f"final s1:          {final_s1}")
    print(f"final s2: \n{final_s2}")
    print(f"final mut_inf: \n{final_mut_inf}")

    return cas_occ, cas_index, final_energy

def deprecated_SelectionOfActiveSpace(json_config, **kwargs):
    """
    (DEPRECATED) This function is used to select the active space using Pegamoid GUI (DEPRECATED)

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
    # try:
    #     call_open_molcas(temp_project_name, molcas_output_directory)
    # except Exception as e:
    #     print(e)
    #     log_file_path = path.join(molcas_output_directory, f'selecting_space_{project_name}.log')
    #     if path.exists(log_file_path):
    #         try:
    #             print_molcas_log_errors(log_file_path, "Timing")
    #         except Exception as e:
    #             print(e)
    #     exit()

    scf_h5_path = f"{find(f'selecting_space_{project_name}.scf.h5',molcas_output_directory)}/selecting_space_{project_name}.scf.h5"

    if not scf_h5_path:
        print("SCF failed, check log file for errors")
        exit()
    scf_h5_abs_path = path.abspath(scf_h5_path)
    # todo provide helpful insight into what orbitals are recommended prior to even loading pegamoid

    run_pegamoid(f1=scf_h5_abs_path)

    sys.exit()
