#!/usr/bin/env python3.10
# Made by Ruben
import os
from multiprocessing import Pool, cpu_count
from os import path, system

from densitypy.Default_Settings.default_config import DEFAULT_BIN_FILE_PATH
from densitypy.charge_migration.chargemigratonscripts import Write_FieldHelp, Write_Pulses, Call_Charge_Migration, \
    Save_Previous_FT, Call_Charge_MigrationFT, Call_Spectrum_Reconstruction_n_Difference, Save_Spectrum_Difference
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


def run_densitypy(json_config_path, study_directory, molcas_input,
                  run_molcas=False, run_charge_migration=False, run_charge_migration_ft=False,
                  run_spectrum_reconstruction=False,
                  field_file_help=False, molcas_input_help=False,
                  lus=False, gridflag: bool or str = True, write_charge_migration=None,
                  debug_mode=False, justh5=False, justgetdipoles=False, justgetdensity=False,
                  weights_file=None, givenfieldfile=None,

                  old_main=False, parallel=False,
                  save_previous=False, make_fortran=False, make_fortran_config=None):
    if make_fortran:
        from densitypy.compiler.fortran_compilation_handler import compile_fortran_code
        compile_fortran_code(**make_fortran_config)

    # TODO too explicit, does not allow for new args to be added easily
    COMMAND_ARGS = [molcas_input, lus, run_charge_migration, run_charge_migration_ft,
                    molcas_input_help, justh5, justgetdensity, justgetdipoles, field_file_help,
                    run_spectrum_reconstruction]

    # Change directory to study_directory
    orig_directory = os.getcwd()
    os.chdir(study_directory)

    # Load Values # todo fixup the json_config and runnable command behavior, is ugly
    json_config = parse_configuration_file(json_config_path)  # If not passed or does not exist will write example file
    runnable_command_available_bool = True if any(COMMAND_ARGS) else False

    if json_config and not runnable_command_available_bool:
        logger.warning(
            f'Found {json_config_path} but nothing was ran. No command line action arguments used. Use -h for help.')

    if molcas_input_help or (molcas_input and not path.exists(molcas_input)):
        create_help_input_file()
        exit()

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
    experiment_directory = charge_migration_settings['experimentdirectory']
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

    # Not sure how to handle the 'write_charge_migration' variable, assuming it is a flag in JSON config

    if write_charge_migration:
        time_delay_range = write_charge_migration
    else:
        change_in_delay = (time_delay_stop - time_delay_start) / (number_of_pp - 1)
        time_delay_range = list(float_range(time_delay_start,
                                            time_delay_stop + change_in_delay, change_in_delay))

    # >OpenMolcas
    if run_molcas:

        n_points = 0
        # Prepare Input
        make_directory(molcas_output_directory)
        #
        # Copies input file and returns keywords in input file #todo should actually parse input , use pymolcas
        true_values = copy_input_file_to_edit(pymolcas_input=molcas_input, project_name=project_name,
                                              molcas_directory=molcas_output_directory)

        copy_file_to(xyz_geometry, molcas_output_directory)

        # gridflag can be
        #   True = use gridit on molcas with with uniform grid points
        #   False = Do not use or create any grid or add gridit to input file
        #   'smart' =  Use our smart grid algorithm to create a grid based on the molecule and step size # todo
        #   'limited' = Use our truncated grid algorithm to create a grid based on the molecule and config settings
        #   'manual' = Use a grid file to create a grid, this is the only option that requires a grid file to be specified in the config file
        if gridflag in ['smart', 'limited']:
            limitedgrid = True if gridflag == 'limited' else False
            n_points = make_better_grid(molcas_output_directory, xyz_geometry, step_size, boundary, limitedgrid)
        elif gridflag == 'manual':
            gridfile = json_config['grid_settings']['default_grid_file']
            system(f"cp {gridfile} {molcas_output_directory}/gridcoord")
            n_points = file_lenth(f'{molcas_output_directory}/gridcoord')
            logger.info(f"Number of Points = {n_points}")
        elif gridflag:
            make_grid_coordinates(molcas_output_directory, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
            n_points = nx * ny * nz
        else:
            assert gridflag is False, f"gridflag is {gridflag} but should be False, True, 'smart', 'limited', or 'manual'"

        if gridflag:
            add_grid_it_to_manual_input_file(molcas_input, project_name, molcas_output_directory, list_of_orbitals,
                                             n_points)
        # >Run Open Molcas
        call_open_molcas(project_name, molcas_output_directory)
        # Extract Data Required for Charge Migration
        copy_file_to(molcas_input, molcas_output_directory + "/")
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

    if field_file_help:
        Write_FieldHelp()

    # >ChargeMigration
    if run_charge_migration:
        if not path.exists(molcas_output_directory):
            logger.error(f"Could not find {molcas_output_directory}")
            exit()
        make_directory(experiment_directory)

        if write_charge_migration:
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
            field_file = givenfieldfile
        else:
            Write_Pulses(f"{experiment_directory}/{field_file}", type_of_pulse_pump, start_time, pump_central_frequency,
                         pump_periods, pump_phase, pump_intensity, pump_polarization,
                         type_of_pulse_probe, time_delay_range, probe_central_frequency,
                         probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration)

        # Run Charge Migration
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        i_excitation = -1  # Just Ground and Prepare
        i_epsilon = -1  # Just Ground and Prepare
        Call_Charge_Migration(DEFAULT_BIN_FILE_PATH, molcas_output_directory, experiment_directory, number_of_times,
                              min_time,
                              max_time,
                              f"{experiment_directory}/{field_file}", ft_time_step, ft_width_step,
                              f"{molcas_output_directory}/{xyz_geometry}", list_of_orbitals, write_charge_migration,
                              Volume,
                              debug_mode, weights_file, dephasing_factor, relaxation_factor, bath_temperature,
                              i_excitation, i_epsilon)

        if old_main:
            pass
        #
        elif parallel:
            #
            INDEX_MAT = []
            #
            for i_excitation in range(1, number_of_states + 1):  # begin iteration
                for i_epsilon in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((i_excitation, i_epsilon))
                    #
            Number_CPU = min(((number_of_states) * 2), os.cpu_count())
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {os.cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_Migration,
                               [(molcas_output_directory, experiment_directory, number_of_times,
                                 min_time, max_time, f"{experiment_directory}/{field_file}", ft_time_step,
                                 ft_width_step, f"{molcas_output_directory}/{xyz_geometry}",
                                 list_of_orbitals, write_charge_migration, Volume, debug_mode, weights_file,
                                 dephasing_factor,
                                 relaxation_factor, bath_temperature, i_excitation, i_epsilon) for
                                i_excitation, i_epsilon
                                in
                                INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            i_excitation = 0
            i_epsilon = 0
            Call_Charge_Migration(DEFAULT_BIN_FILE_PATH, molcas_output_directory, experiment_directory, number_of_times,
                                  min_time,
                                  max_time, f"{experiment_directory}/{field_file}", ft_time_step, ft_width_step,
                                  f"{molcas_output_directory}/{xyz_geometry}", list_of_orbitals,
                                  write_charge_migration, Volume, debug_mode, weights_file, dephasing_factor,
                                  relaxation_factor,
                                  bath_temperature, i_excitation, i_epsilon)
        copy_file_to(json_config_path, f"{experiment_directory}/")
    # >ChargeMigrationFT
    if run_charge_migration_ft:
        ##
        if save_previous:
            Save_Previous_FT(experiment_directory, dephasing_factor, relaxation_factor, time_delay_start,
                             time_delay_stop,
                             min_omegas, max_omegas, pump_periods, probe_periods, pump_intensity,
                             probe_intensity, number_of_pp, pump_phase, probe_phase, ft_time_step, ft_width_step,
                             pump_polarization, probe_polarization, number_of_omegas, min_tau_omega,
                             max_tau_omega)
        if givenfieldfile:
            field_file = givenfieldfile
        else:
            # Make Field File
            Write_Pulses(f"{experiment_directory}/{field_file}", type_of_pulse_pump, start_time, pump_central_frequency,
                         pump_periods, pump_phase, pump_intensity, pump_polarization,
                         type_of_pulse_probe, time_delay_range, probe_central_frequency,
                         probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration)
        # Run Charge Migration FT
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if parallel:
            #
            INDEX_MAT = []
            #
            for i_excitation in range(1, number_of_states + 1):  # begin iteration
                for i_epsilon in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((i_excitation, i_epsilon))
                    #
            Number_CPU = min(((number_of_states) * 2), (cpu_count()) / 2)
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_MigrationFT,
                               [(molcas_output_directory, experiment_directory,
                                 f"{molcas_output_directory}/{xyz_geometry}", number_of_omegas,
                                 min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                 ft_time_step, ft_width_step, f"{experiment_directory}/{field_file}", debug_mode,
                                 i_excitation, i_epsilon) for
                                i_excitation, i_epsilon in INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            i_excitation = 0
            i_epsilon = 0
            Call_Charge_MigrationFT(DEFAULT_BIN_FILE_PATH, molcas_output_directory, experiment_directory,
                                    f"{molcas_output_directory}/{xyz_geometry}",
                                    number_of_omegas,
                                    min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                    ft_time_step, ft_width_step, f"{experiment_directory}/{field_file}", debug_mode,
                                    i_excitation, i_epsilon)
    # >SpectrumReconstruction
    if run_spectrum_reconstruction:
        if save_previous:
            Save_Spectrum_Difference(DEFAULT_BIN_FILE_PATH, experiment_directory, 'difference_' + experiment_directory,
                                     dephasing_factor, relaxation_factor, time_delay_start, time_delay_stop, min_omegas,
                                     max_omegas, pump_periods, probe_periods, pump_intensity, probe_intensity,
                                     number_of_pp, pump_phase, probe_phase, ft_time_step, ft_width_step,
                                     pump_polarization, probe_polarization, number_of_omegas, min_tau_omega,
                                     max_tau_omega)
        #
        Call_Spectrum_Reconstruction_n_Difference(DEFAULT_BIN_FILE_PATH, molcas_output_directory, experiment_directory,
                                                  f"{molcas_output_directory}/{xyz_geometry}", number_of_omegas,
                                                  min_omegas, max_omegas,
                                                  number_of_tau_omegas, min_tau_omega, max_tau_omega, debug_mode)

        # logger.info("Creating difference_" + experiment_directory + " to compare the Dipole and Charge Spectra")
        # Dipole_Charge_Comparison(experiment_directory + '/Dipole/DipoleFT_ww', experiment_directory +
        #                          '/Dipole/DipoleFT_ww_reconstructed', 'difference_' + experiment_directory)

    # Change back to original directory
    os.chdir(orig_directory)


if __name__ == "__main__":
    make_fortran_config = dict(directory='/home/ruben/PycharmProjects/DensityPy/densityfort',
                               make_flags='all DEB_FLAG=d',
                               # ... pther settings but this is mvp
                               )

    run_densitypy(json_config_path="configuration_help.json",
                  study_directory="/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/",
                  molcas_input='molcas_input_help.input', run_molcas=False, run_charge_migration=True,
                  run_charge_migration_ft=False, run_spectrum_reconstruction=False, field_file_help=False,
                  molcas_input_help=False, lus=False, gridflag=True, write_charge_migration=None, debug_mode=True,
                  justh5=False, justgetdipoles=False, justgetdensity=False, weights_file=None, givenfieldfile=None,
                  old_main=True, parallel=False, save_previous=False,
                  make_fortran=True, make_fortran_config=make_fortran_config
                  )
# Todo:
#     refactor the code to avoid using os.chdir entirely, since it changes the state of the Python process and can potentially lead to confusing bugs. This can be achieved by using absolute paths in the shell commands instead of changing the working directory. However, the feasibility of this refactoring depends on the specific codebase and build system.
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
# Document my subroutine, keep in mind that the current comments might not actually represent the subroutine, so use them as context but the code is the final word when it comes to documentation:
