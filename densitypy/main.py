#!/usr/bin/env python3.10
# Made by Ruben
import os
import time
from multiprocessing import Pool, cpu_count
from os import path, system

from densitypy.charge_migration.chargemigratonscripts import Write_FieldHelp, Write_Pulses, Call_Charge_Migration, \
    Save_Previous_FT, Call_Charge_MigrationFT, Call_Spectrum_Reconstruction_n_Difference, Save_Spectrum_Difference
from densitypy.molcas.molcas_log_handler import DipolesProject
from densitypy.molcas.molcasscripts import CreateHelpInputFile, CopyInputFileToEdit, Make_Better_Grid, \
    Make_Grid_Coordinates, AddGridItToManualInputFile, Call_OpenMolcas, GetDipolesFromLogFile, Make_MU_HeatMap, \
    ExtractGridDensity, LoadFromh5File
from densitypy.molcas.selectionofactivespace import SelectionOfActiveSpace
from densitypy.project_utils.configuration_parser import parse_configuration_file
from densitypy.project_utils.def_functions import float_range, make_directory, file_lenth, copy_file_to, \
    find, check_compatibility_of_arguements, get_value_of_as_string
from densitypy.project_utils.logger import logger


# TODO FIXME this is not used and .Molcas is not created by default meaning that this will fail unless user correctly sets up OpenMolcas
def get_molcas_molcasrc_variables():
    """Get the environment variables for Molcas"""
    pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")

    if path.exists(pathtomolcasrc):
        molcas_workdir = get_value_of_as_string(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
        print("Molcas Work Directory = " + molcas_workdir)
    else:
        print("Molcas configuration file ~/.Molcas/molcasrc not found. Use pymolcas --setup to create one")
        exit()
    exit()


def run(json_config_path, pymolcas_input, study_directory,
        run_Molcas=False,
        run_ChargeMigration=True, run_ChargeMigrationFT=False,
        run_SpectrumReconstruction=False, lus=False, debug_mode=False,
        pymolcas_input_help=False, justh5=False, justgetdipoles=False,
        justgetdensity=False, writeCM=None, givenfieldfile=None,
        weights_file=None, fieldfilehelp=False,
        gridflag: bool or str = True,

        old_main=False, parallel=False,
        save_previous=False):
    # Change directory to study_directory
    orig_directory = os.getcwd()
    os.chdir(study_directory)

    if pymolcas_input_help:
        CreateHelpInputFile()

    # Load Values # TODO too explicit, does not allow for new args to be added easily
    COMMAND_ARGS = [pymolcas_input, lus, run_ChargeMigration, run_ChargeMigrationFT,
                    pymolcas_input_help, justh5, justgetdensity, justgetdipoles, fieldfilehelp,
                    run_SpectrumReconstruction]
    json_config = parse_configuration_file(json_config_path)  # If not passed or does not exist will write example file
    runnable_command_available_bool = True if any(COMMAND_ARGS) else False

    if json_config and not runnable_command_available_bool:
        logger.warning(
            f'Found {json_config_path} but nothing was ran. No command line action arguments used. Use -h for help.')

    if lus:
        # Selection of Active Space using lus argument. Requires Luscus
        pymolcas_input = check_compatibility_of_arguements("pymolcas_input", "lus")
        SelectionOfActiveSpace(json_config)

    # Name of Project
    project_settings = json_config['project_settings']
    project_name = project_settings['project_name']
    xyz_geometry = project_settings['xyz_molecule_geometry']
    number_of_states = project_settings['number_of_states']
    list_of_orbitals = project_settings['list_of_active_orbitals']
    molcas_output_directory = project_settings['molcas_output_directory']
    bin_directory = project_settings['bin_directory']

    # Grid Settings
    grid_settings = json_config['grid_settings']
    nx = grid_settings['number_of_points']['x_axis']
    ny = grid_settings['number_of_points']['y_axis']
    nz = grid_settings['number_of_points']['z_axis']
    xmin = grid_settings['x_range']['min']
    xmax = grid_settings['x_range']['max']
    ymin = grid_settings['y_range']['min']
    ymax = grid_settings['y_range']['max']
    zmin = grid_settings['z_range']['min']
    zmax = grid_settings['z_range']['max']
    step_size = grid_settings['step_size']
    boundary = grid_settings['boundary']

    # Charge Migration Parameters
    charge_migration_settings = json_config['charge_migration_settings']
    output_directory = charge_migration_settings['output_directory']
    config_field_file = charge_migration_settings['field_file']
    number_of_times = charge_migration_settings['number_of_times']
    min_time = charge_migration_settings['time_range']['min']
    max_time = charge_migration_settings['time_range']['max']
    bath_temp = charge_migration_settings['bath_temperature']
    dephasing_factor = charge_migration_settings['dephasing_factor']
    relaxing_factor = charge_migration_settings['relaxation_factor']

    # Pump Settings
    pump_settings = json_config['pulses_settings']['pump']
    type_of_pulse_pump = pump_settings['type_of_pulse']
    start_time = pump_settings['start_time']
    pump_central_frequency = pump_settings['central_frequency']
    pump_periods = pump_settings['periods']
    pump_phase = pump_settings['phase']
    pump_intensity = pump_settings['intensity']
    pump_polarization = pump_settings['polarization']

    # Probe Settings
    probe_settings = json_config['pulses_settings']['probe']
    type_of_pulse_probe = probe_settings['type_of_pulse']
    time_delay_start = probe_settings['time_delay_range']['start']
    time_delay_stop = probe_settings['time_delay_range']['stop']
    number_of_pp = probe_settings['number_of_pp']
    probe_central_frequency = probe_settings['central_frequency']
    probe_periods = probe_settings['periods']
    probe_phase = probe_settings['phase']
    probe_intensity = probe_settings['intensity']
    probe_polarization = probe_settings['polarization']

    # Charge Migration FT Parameters
    charge_migration_ft_settings = json_config['charge_migration_ft_settings']
    number_of_omegas = charge_migration_ft_settings['number_of_omegas']
    Min_Omegas = charge_migration_ft_settings['omega_range']['min']
    Max_Omegas = charge_migration_ft_settings['omega_range']['max']
    Number_of_TauOmegas = charge_migration_ft_settings['number_of_tauomegas']
    Min_TauOmega = charge_migration_ft_settings['tauomega_range']['min']
    Max_TauOmega = charge_migration_ft_settings['tauomega_range']['max']
    TimeStep_FT = charge_migration_ft_settings['timestep_ft']
    WidthStep_FT = charge_migration_ft_settings['widthstep_ft']
    Volume = step_size * step_size * step_size

    logger.info("Unit of Volume = " + str(Volume))

    # Not sure how to handle the 'writeCM' variable, assuming it is a flag in JSON config

    if writeCM:
        time_delay_range = writeCM
    else:
        change_in_delay = (time_delay_stop - time_delay_start) / (number_of_pp - 1)
        time_delay_range = list(float_range(time_delay_start,
                                            time_delay_stop + change_in_delay, change_in_delay))

    # >OpenMolcas
    if run_Molcas:
        N_points = 0
        # Prepare Input
        make_directory(molcas_output_directory)
        # MakeDirectoryNoDelete(Molcas_Output_Directory)
        #
        # Copies input file and returns keywords in input file
        true_values = CopyInputFileToEdit(pymolcas_input, project_name, molcas_output_directory)

        # gridflag can be
        #   True = use gridit on molcas with with uniform grid points
        #   False = Do not use or create any grid or add gridit to input file
        #   'smart' =  Use our smart grid algorithm to create a grid based on the molecule and step size # todo
        #   'limited' = Use our truncated grid algorithm to create a grid based on the molecule and config settings
        #   'manual' = Use a grid file to create a grid, this is the only option that requires a grid file to be specified in the config file
        if gridflag in ['smart', 'limited']:
            limitedgrid = True if gridflag == 'limited' else False
            N_points = Make_Better_Grid(molcas_output_directory, xyz_geometry, step_size, boundary, limitedgrid)
        elif gridflag == 'manual':
            gridfile = json_config['grid_settings']['default_grid_file']
            system(f"cp {gridfile} {molcas_output_directory}/gridcoord")
            N_points = file_lenth(f'{molcas_output_directory}/gridcoord')
            logger.info(f"Number of Points = {N_points}")
        elif gridflag:
            Make_Grid_Coordinates(molcas_output_directory, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
            N_points = nx * ny * nz
        else:
            assert gridflag is False, f"gridflag is {gridflag} but should be False, True, 'smart', 'limited', or 'manual'"

        if gridflag:
            AddGridItToManualInputFile(pymolcas_input, project_name, molcas_output_directory, list_of_orbitals,
                                       N_points)
        # >Run Open Molcas
        Call_OpenMolcas(project_name, molcas_output_directory)
        # Extract Data Required for Charge Migration
        copy_file_to(pymolcas_input, molcas_output_directory + "/")
        logfilepath = find(project_name + ".log", ".", molcas_output_directory)
        copy_file_to(logfilepath + "/" + project_name + ".log ", molcas_output_directory + "/")
        copy_file_to(xyz_geometry, molcas_output_directory + "/")

        if "RASSI" in true_values:
            GetDipolesFromLogFile(project_name, molcas_output_directory, number_of_states, molcas_output_directory)
            Make_MU_HeatMap(molcas_output_directory)
            project = DipolesProject(log_path=f'{molcas_output_directory}/{project_name}.log')
            project.write_csvs(directory=molcas_output_directory)
            fig = project.get_mu_heatmaps()
            fig.savefig(f'{molcas_output_directory}/dipole-heatmap.png')

        if gridflag:
            ExtractGridDensity(list_of_orbitals, project_name, molcas_output_directory, molcas_output_directory)

        if "RASSCF" in true_values:
            LoadFromh5File(project_name, molcas_output_directory, molcas_output_directory, true_values,
                           "DENSITY_MATRIX",
                           "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y",
                           "AO_MLTPL_Z", "MO_ENERGIES", "MO_VECTORS", justh5)

    # Useful Flags for e.g. debugging
    if justh5:
        true_values = []
        LoadFromh5File(project_name, molcas_output_directory, molcas_output_directory, true_values, "DENSITY_MATRIX",
                       "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y", "AO_MLTPL_Z",
                       "MO_ENERGIES", "MO_VECTORS", justh5)
    if justgetdensity:
        ExtractGridDensity(list_of_orbitals, project_name, molcas_output_directory, molcas_output_directory)

    if justgetdipoles:
        GetDipolesFromLogFile(project_name, molcas_output_directory, number_of_states, molcas_output_directory)
        Make_MU_HeatMap(molcas_output_directory)

    # >ChargeMigration
    if fieldfilehelp:
        Write_FieldHelp()
    if run_ChargeMigration:
        cm_time = time.time()
        if writeCM:
            pass
        elif (weights_file == None):
            if path.exists("Weights_File"):
                logger.info("Found Weight_File in curent directory, did you forget to use \"-w\"")
                userinput_weights = input("\"Y\" to use Weights_File:").upper().replace(
                    " ", "")
                yes_input = ["yes", "y", "Yes", "Y", "YES"]
                if (userinput_weights in yes_input):
                    weights_file = "Weights_File"
                    logger.info("Using \"-w\" flag")
        if givenfieldfile:
            Field_File = givenfieldfile
        else:
            Field_File = f"{molcas_output_directory}/Field_pulses"
            # Make Field File
            Write_Pulses(Field_File, type_of_pulse_pump, start_time, pump_central_frequency,
                         pump_periods, pump_phase, pump_intensity, pump_polarization,
                         type_of_pulse_probe, time_delay_range, probe_central_frequency,
                         probe_periods, probe_phase, probe_intensity, probe_polarization, writeCM)

        # Run Charge Migration
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        iExcitation = -1  # Just Ground and Prepare
        iEPSILON = -1  # Just Ground and Prepare
        Call_Charge_Migration(bin_directory, molcas_output_directory, output_directory, number_of_times, min_time,
                              max_time,
                              Field_File, TimeStep_FT, WidthStep_FT, xyz_geometry, list_of_orbitals, writeCM, Volume,
                              debug_mode, weights_file, dephasing_factor, relaxing_factor, bath_temp, iExcitation,
                              iEPSILON)

        if old_main:
            pass
        #
        elif parallel:
            #
            INDEX_MAT = []
            #
            for iExcitation in range(1, number_of_states + 1):  # begin iteration
                for iEPSILON in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((iExcitation, iEPSILON))
                    #
            Number_CPU = min(((number_of_states) * 2), os.cpu_count())
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {os.cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_Migration,
                               [(molcas_output_directory, output_directory, number_of_times,
                                 min_time, max_time, Field_File, TimeStep_FT, WidthStep_FT, xyz_geometry,
                                 list_of_orbitals, writeCM, Volume, debug_mode, weights_file, dephasing_factor,
                                 relaxing_factor, bath_temp, iExcitation, iEPSILON) for iExcitation, iEPSILON in
                                INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            Call_Charge_Migration(bin_directory, molcas_output_directory, output_directory, number_of_times, min_time,
                                  max_time, Field_File, TimeStep_FT, WidthStep_FT, xyz_geometry, list_of_orbitals,
                                  writeCM, Volume, debug_mode, weights_file, dephasing_factor, relaxing_factor,
                                  bath_temp, iExcitation, iEPSILON)
        copy_file_to(json_config, f"{output_directory}/")
    # >ChargeMigrationFT
    if run_ChargeMigrationFT:
        ##
        if save_previous:
            Save_Previous_FT(output_directory, dephasing_factor, relaxing_factor, time_delay_start,
                             time_delay_stop,
                             Min_Omegas, Max_Omegas, pump_periods, probe_periods, pump_intensity,
                             probe_intensity, number_of_pp, pump_phase, probe_phase, TimeStep_FT, WidthStep_FT,
                             pump_polarization, probe_polarization, number_of_omegas, Min_TauOmega,
                             Max_TauOmega)
        if givenfieldfile:
            Field_File = givenfieldfile
        else:
            Field_File = f"{molcas_output_directory}/Field_pulses"
            # Make Field File
            Write_Pulses(Field_File, type_of_pulse_pump, start_time, pump_central_frequency,
                         pump_periods, pump_phase, pump_intensity, pump_polarization,
                         type_of_pulse_probe, time_delay_range, probe_central_frequency,
                         probe_periods, probe_phase, probe_intensity, probe_polarization, writeCM)
        # Run Charge Migration FT
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if parallel:
            #
            INDEX_MAT = []
            #
            for iExcitation in range(1, number_of_states + 1):  # begin iteration
                for iEPSILON in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((iExcitation, iEPSILON))
                    #
            # system(f'rm OutDir/Dipole/Dipole_State_iExcitation')
            Number_CPU = min(((number_of_states) * 2), (cpu_count()) / 2)
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_MigrationFT,
                               [(molcas_output_directory, output_directory, xyz_geometry, number_of_omegas,
                                 Min_Omegas, Max_Omegas, Number_of_TauOmegas, Min_TauOmega, Max_TauOmega,
                                 TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON) for
                                iExcitation, iEPSILON in INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            Call_Charge_MigrationFT(bin_directory, molcas_output_directory, output_directory, xyz_geometry,
                                    number_of_omegas,
                                    Min_Omegas, Max_Omegas, Number_of_TauOmegas, Min_TauOmega, Max_TauOmega,
                                    TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON)
    # >SpectrumReconstruction
    if run_SpectrumReconstruction:
        if save_previous:
            Save_Spectrum_Difference(bin_directory, output_directory, 'difference_' + output_directory,
                                     dephasing_factor, relaxing_factor, time_delay_start, time_delay_stop, Min_Omegas,
                                     Max_Omegas, pump_periods, probe_periods, pump_intensity, probe_intensity,
                                     number_of_pp, pump_phase, probe_phase, TimeStep_FT, WidthStep_FT,
                                     pump_polarization, probe_polarization, number_of_omegas, Min_TauOmega,
                                     Max_TauOmega)
        #
        Call_Spectrum_Reconstruction_n_Difference(bin_directory, molcas_output_directory, output_directory,
                                                  xyz_geometry,
                                                  number_of_omegas, Min_Omegas, Max_Omegas, Number_of_TauOmegas,
                                                  Min_TauOmega, Max_TauOmega, debug_mode)

    # Change back to original directory
    os.chdir(orig_directory)


if __name__ == "__main__":
    run(json_config_path="chargemigration.json",
        study_directory="/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/",

        # Molcas Group
        run_Molcas=False,
        pymolcas_input='inputhelp.input',  # todo add the same mechanism that json_config has for this
        pymolcas_input_help=False,  # todo ^^^^ but do allow for explicit flag to be set
        gridflag=True,  # todo test this for the other options
        lus=False,  # todo test this option... think remove entirely for pymol or molcas-lus

        # Charge Migration Group
        run_ChargeMigration=True,
        run_ChargeMigrationFT=False,
        writeCM=None,  # todo test this option
        weights_file=None,  # todo test this option
        parallel=False,  # todo test this option
        old_main=True,
        givenfieldfile=None,
        fieldfilehelp=False,
        save_previous=False,  # todo test this option

        # Analysis Group
        run_SpectrumReconstruction=False,
        justh5=False, justgetdipoles=False, justgetdensity=False,

        # Development Group
        debug_mode=True,

        )
# todo
#   -update documentation
#   -analytics module (feature-detection, spectrum reconstruction, etc...)
#       - completely redo the parallelization of the code, this goes hand in hadn with the C_i coefficients reconstruction module (see documentation)
#   -refactor from execute_command to f2py calls for fotran code
#       -does not compile all the dependencies, need to write a script to do so
#   -automate f2py compilation
#   -name args when passing to functions
#   -modularize fortran code
#       -clean output files and directories (molcas,chargemigration (ft) done, rest ... to-do)
#   -change functions to lower case, change all variables to lower case, classes to upper case, etc... I followed
#   -update cli_v1.py with new changes to main.py
#   -make/generate tests
#   -main.run(...)-->GUI (Streamlit/Wave)
#
# todo
#   reduce code to main__
