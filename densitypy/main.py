#!/usr/bin/env python3.10
# Made by Ruben
import os
from multiprocessing import Pool, cpu_count
from os import path, system

from densitypy.Default_Settings.default_config import DEFAULT_BIN_FILE_PATH
from densitypy.charge_migration.chargemigratonscripts import Write_FieldHelp, Write_Pulses, Call_Charge_Migration, \
    Save_Previous_FT, Call_Charge_MigrationFT, Call_Spectrum_Reconstruction_n_Difference, Save_Spectrum_Difference, \
    Dipole_Charge_Comparison
from densitypy.molcas.molcas_log_handler import DipolesProject
from densitypy.molcas.molcasscripts import create_help_input_file, copy_input_file_to_edit, make_better_grid, \
    make_grid_coordinates, add_grid_it_to_manual_input_file, call_open_molcas, get_dipoles_from_log_file, \
    make_mu_heat_map, \
    extract_grid_density, load_fromh5_file
from densitypy.molcas.selectionofactivespace import SelectionOfActiveSpace
from densitypy.project_utils.configuration_parser import parse_configuration_file
from densitypy.project_utils.def_functions import float_range, make_directory, file_lenth, copy_file_to, find
from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])

def run(json_config_path, study_directory, pymolcas_input,
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

    # Load Values # TODO too explicit, does not allow for new args to be added easily
    COMMAND_ARGS = [pymolcas_input, lus, run_ChargeMigration, run_ChargeMigrationFT,
                    pymolcas_input_help, justh5, justgetdensity, justgetdipoles, fieldfilehelp,
                    run_SpectrumReconstruction]
    json_config = parse_configuration_file(json_config_path)  # If not passed or does not exist will write example file
    runnable_command_available_bool = True if any(COMMAND_ARGS) else False

    if json_config and not runnable_command_available_bool:
        logger.warning(
            f'Found {json_config_path} but nothing was ran. No command line action arguments used. Use -h for help.')

    if pymolcas_input_help or not path.exists(pymolcas_input):
        create_help_input_file()

    if lus:
        # Selection of Active Space using lus argument. Requires Luscus
        SelectionOfActiveSpace(json_config)  # todo need to update and most likely will switch programs

    # Split JSON config into sections
    project_settings = json_config['projectsettings']
    grid_settings = json_config['gridsettings']
    charge_migration_settings = json_config['chargemigrationsettings']
    pump_settings = json_config['pumppulsessettings']
    probe_settings = json_config['probepulsessettings']
    charge_migration_ft_settings = json_config['chargemigrationftsettings']

    # Project Settings
    project_name = project_settings['projectname']
    xyz_geometry = project_settings['xyzmoleculegeometry']
    number_of_states = project_settings['numberofstates']
    list_of_orbitals = project_settings['listofactiveorbitals']
    molcas_output_directory = project_settings['molcasoutputdirectory']

    # Grid Settings
    nx = grid_settings['numberofpointsxaxis']
    ny = grid_settings['numberofpointsyaxis']
    nz = grid_settings['numberofpointszaxis']
    xmin = grid_settings['xmin']
    xmax = grid_settings['xmax']
    ymin = grid_settings['ymin']
    ymax = grid_settings['ymax']
    zmin = grid_settings['zmin']
    zmax = grid_settings['zmax']
    step_size = grid_settings['stepsize']
    boundary = grid_settings['boundary']

    # Charge Migration Parameters
    output_directory = charge_migration_settings['outputdirectory']
    field_file = charge_migration_settings['fieldfile']
    number_of_times = charge_migration_settings['numberoftimes']
    min_time = charge_migration_settings['mintime']
    max_time = charge_migration_settings['maxtime']
    bath_temperature = charge_migration_settings['bathtemperature']
    dephasing_factor = charge_migration_settings['dephasingfactor']
    relaxation_factor = charge_migration_settings['relaxationfactor']

    # Pump Settings
    type_of_pulse_pump = pump_settings['typeofpulse']
    start_time = pump_settings['starttime']
    pump_central_frequency = pump_settings['pumpcentralfrequency']
    pump_periods = pump_settings['pumpperiods']
    pump_phase = pump_settings['pumpphase']
    pump_intensity = pump_settings['pumpintensity']
    pump_polarization = pump_settings['pumppolarization']

    # Probe Settings
    type_of_pulse_probe = probe_settings['typeofpulse']
    time_delay_start = probe_settings['timedelaystart']
    time_delay_stop = probe_settings['timedelaystop']
    number_of_pp = probe_settings['numberofpp']
    probe_central_frequency = probe_settings['probecentralfrequency']
    probe_periods = probe_settings['probeperiods']
    probe_phase = probe_settings['probephase']
    probe_intensity = probe_settings['probeintensity']
    probe_polarization = probe_settings['probepolarization']

    # Charge Migration FT Parameters
    number_of_omegas = charge_migration_ft_settings['numberofomegas']
    min_omegas = charge_migration_ft_settings['minomega']
    max_omegas = charge_migration_ft_settings['maxomega']
    number_of_tau_omegas = charge_migration_ft_settings['numberoftauomegas']
    min_tau_omega = charge_migration_ft_settings['mintauomega']
    max_tau_omega = charge_migration_ft_settings['maxtauomega']
    ft_time_step = charge_migration_ft_settings['fttimestep']
    ft_width_step = charge_migration_ft_settings['ftwidthstep']
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
        #
        # Copies input file and returns keywords in input file #todo should actually parse input , use pymolcas
        true_values = copy_input_file_to_edit(pymolcas_input, project_name, molcas_output_directory)

        copy_file_to(xyz_geometry, molcas_output_directory)

        # gridflag can be
        #   True = use gridit on molcas with with uniform grid points
        #   False = Do not use or create any grid or add gridit to input file
        #   'smart' =  Use our smart grid algorithm to create a grid based on the molecule and step size # todo
        #   'limited' = Use our truncated grid algorithm to create a grid based on the molecule and config settings
        #   'manual' = Use a grid file to create a grid, this is the only option that requires a grid file to be specified in the config file
        if gridflag in ['smart', 'limited']:
            limitedgrid = True if gridflag == 'limited' else False
            N_points = make_better_grid(molcas_output_directory, xyz_geometry, step_size, boundary, limitedgrid)
        elif gridflag == 'manual':
            gridfile = json_config['grid_settings']['default_grid_file']
            system(f"cp {gridfile} {molcas_output_directory}/gridcoord")
            N_points = file_lenth(f'{molcas_output_directory}/gridcoord')
            logger.info(f"Number of Points = {N_points}")
        elif gridflag:
            make_grid_coordinates(molcas_output_directory, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
            N_points = nx * ny * nz
        else:
            assert gridflag is False, f"gridflag is {gridflag} but should be False, True, 'smart', 'limited', or 'manual'"

        if gridflag:
            add_grid_it_to_manual_input_file(pymolcas_input, project_name, molcas_output_directory, list_of_orbitals,
                                             N_points)
        # >Run Open Molcas
        call_open_molcas(project_name, molcas_output_directory)
        # Extract Data Required for Charge Migration
        copy_file_to(pymolcas_input, molcas_output_directory + "/")
        logfilepath = find(project_name + ".log", ".", molcas_output_directory)
        copy_file_to(logfilepath + "/" + project_name + ".log ", molcas_output_directory + "/")
        copy_file_to(xyz_geometry, molcas_output_directory + "/")

        if "RASSI" in true_values:
            get_dipoles_from_log_file(project_name, molcas_output_directory, number_of_states, molcas_output_directory)
            make_mu_heat_map(molcas_output_directory)
            project = DipolesProject(log_path=f'{molcas_output_directory}/{project_name}.log')
            project.write_csvs(directory=molcas_output_directory)
            fig = project.get_mu_heatmaps()
            fig.savefig(f'{molcas_output_directory}/dipole-heatmap.png')

        if gridflag:
            extract_grid_density(list_of_orbitals, project_name, molcas_output_directory, molcas_output_directory)

        if "RASSCF" in true_values:
            load_fromh5_file(project_name, molcas_output_directory, molcas_output_directory, true_values,
                             "DENSITY_MATRIX",
                             "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y",
                             "AO_MLTPL_Z", "MO_ENERGIES", "MO_VECTORS", justh5)

    # Useful Flags for e.g. debugging
    if justh5:
        true_values = []
        load_fromh5_file(project_name, molcas_output_directory, molcas_output_directory, true_values, "DENSITY_MATRIX",
                         "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y", "AO_MLTPL_Z",
                         "MO_ENERGIES", "MO_VECTORS", justh5)
    if justgetdensity:
        extract_grid_density(list_of_orbitals, project_name, molcas_output_directory, molcas_output_directory)

    if justgetdipoles:
        get_dipoles_from_log_file(project_name, molcas_output_directory, number_of_states, molcas_output_directory)
        make_mu_heat_map(molcas_output_directory)

    # >ChargeMigration
    if fieldfilehelp:
        Write_FieldHelp()
    if run_ChargeMigration:
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
        Call_Charge_Migration(DEFAULT_BIN_FILE_PATH, molcas_output_directory, output_directory, number_of_times,
                              min_time,
                              max_time,
                              Field_File, ft_time_step, ft_width_step, xyz_geometry, list_of_orbitals, writeCM, Volume,
                              debug_mode, weights_file, dephasing_factor, relaxation_factor, bath_temperature,
                              iExcitation, iEPSILON)

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
                                 min_time, max_time, Field_File, ft_time_step, ft_width_step, xyz_geometry,
                                 list_of_orbitals, writeCM, Volume, debug_mode, weights_file, dephasing_factor,
                                 relaxation_factor, bath_temperature, iExcitation, iEPSILON) for iExcitation, iEPSILON
                                in
                                INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            Call_Charge_Migration(DEFAULT_BIN_FILE_PATH, molcas_output_directory, output_directory, number_of_times,
                                  min_time,
                                  max_time, Field_File, ft_time_step, ft_width_step, xyz_geometry, list_of_orbitals,
                                  writeCM, Volume, debug_mode, weights_file, dephasing_factor, relaxation_factor,
                                  bath_temperature, iExcitation, iEPSILON)
        copy_file_to(json_config, f"{output_directory}/")
    # >ChargeMigrationFT
    if run_ChargeMigrationFT:
        ##
        if save_previous:
            Save_Previous_FT(output_directory, dephasing_factor, relaxation_factor, time_delay_start,
                             time_delay_stop,
                             min_omegas, max_omegas, pump_periods, probe_periods, pump_intensity,
                             probe_intensity, number_of_pp, pump_phase, probe_phase, ft_time_step, ft_width_step,
                             pump_polarization, probe_polarization, number_of_omegas, min_tau_omega,
                             max_tau_omega)
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
            Number_CPU = min(((number_of_states) * 2), (cpu_count()) / 2)
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_MigrationFT,
                               [(molcas_output_directory, output_directory, xyz_geometry, number_of_omegas,
                                 min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                 ft_time_step, ft_width_step, Field_File, debug_mode, iExcitation, iEPSILON) for
                                iExcitation, iEPSILON in INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            Call_Charge_MigrationFT(DEFAULT_BIN_FILE_PATH, molcas_output_directory, output_directory, xyz_geometry,
                                    number_of_omegas,
                                    min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                    ft_time_step, ft_width_step, Field_File, debug_mode, iExcitation, iEPSILON)
    # >SpectrumReconstruction
    if run_SpectrumReconstruction:
        if save_previous:
            Save_Spectrum_Difference(DEFAULT_BIN_FILE_PATH, output_directory, 'difference_' + output_directory,
                                     dephasing_factor, relaxation_factor, time_delay_start, time_delay_stop, min_omegas,
                                     max_omegas, pump_periods, probe_periods, pump_intensity, probe_intensity,
                                     number_of_pp, pump_phase, probe_phase, ft_time_step, ft_width_step,
                                     pump_polarization, probe_polarization, number_of_omegas, min_tau_omega,
                                     max_tau_omega)
        #
        Call_Spectrum_Reconstruction_n_Difference(DEFAULT_BIN_FILE_PATH, molcas_output_directory, output_directory,
                                                  xyz_geometry,
                                                  number_of_omegas, min_omegas, max_omegas, number_of_tau_omegas,
                                                  min_tau_omega, max_tau_omega, debug_mode)

        # logger.info("Creating difference_" + SimOut + " to compare the Dipole and Charge Spectra")
        # Dipole_Charge_Comparison(SimOut + '/Dipole/DipoleFT_ww', SimOut +
        #                          '/Dipole/DipoleFT_ww_reconstructed', 'difference_' + SimOut)

    # Change back to original directory
    os.chdir(orig_directory)


if __name__ == "__main__":
    run(
        study_directory="/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/",
        json_config_path="configuration_help.json",
        # Molcas Group
        run_Molcas=True,
        pymolcas_input='molcas_input_help.input',  # todo add the same mechanism that json_config has for this
        pymolcas_input_help=False,  # todo ^^^^ but do allow for explicit flag to be set
        gridflag=True,  # todo test this for the other options
        lus=False,  # todo test this option... think remove entirely for pymol or molcas-lus

        # Charge Migration Group
        run_ChargeMigration=True,
        run_ChargeMigrationFT=True,
        writeCM=None,  # todo test this option
        weights_file=None,  # todo test this option
        parallel=False,  # todo test this option
        old_main=True,
        givenfieldfile=None,
        fieldfilehelp=False,
        save_previous=False,  # todo test this option

        # Analysis Group
        run_SpectrumReconstruction=True,
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
#   -implement more cross validation checks of settings between runs, between modules, and between configurations of different programs like molcas etc...
#   -make sure inputs are valid and available (check for existence of files, etc...) before running anything if their appropriate flags are set
#       - very important specifically is dealing with orbitals, stats, CI roots etc... e.g. only 4 active orbitals with 6 electrons, 3 filled and one excited should limit the max CIroots depending on number of configurations available
#   reduce code to main__
