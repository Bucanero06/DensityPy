#!/usr/bin/env python3.10
# Made by Ruben
import os
from multiprocessing import Pool, cpu_count
from os import path, system

from densitypy.Default_Settings.default_config import DEFAULT_BIN_FILE_PATH
from densitypy.charge_migration.chargemigratonscripts import Write_FieldHelp, Write_Pulses, Call_Charge_Migration, \
    Call_Charge_MigrationFT, Call_Spectrum_Reconstruction_n_Difference,  generate_time_delays
from densitypy.molcas.DipolesLogParser import DipolesLogParser
from densitypy.molcas.molcasscripts import create_help_input_file, copy_and_parse_molcas_input_file_to_edit, \
    make_better_grid, \
    make_grid_coordinates, add_grid_it_to_manual_input_file, call_open_molcas, parse_project_grid_file, \
    load_project_rasscf_h5_file, write_grid_density_file
from densitypy.molcas.selectionofactivespace import SelectionOfActiveSpace
from densitypy.post_processing.plotting_module import plot_dipoles_v_time, plot_ft_pulses, plot_pulses, \
    plot_2d_spectrum, \
    plot_ft_all_dipoles_v_time, plot_2d_spectrum_peak_analysis, plot_2d_spectrum_interactive, \
    plot_atomic_dipoles_v_time, \
    difference_between_dipole_and_atomic_charges_v_time, plot_ft_all_atomic_dipoles_v_time
from densitypy.project_utils.configuration_parser import parse_configuration_file
from densitypy.project_utils.file_directory_ops import change_directory_manager, make_directory, copy_file_to, \
    file_lenth, find
from densitypy.project_utils.logger import setup_logger, DEFAULT_LOGGIN_VALUES

logger = setup_logger(__name__.split('.')[-1])


def run_densitypy(json_config_path, study_directory, molcas_input,
                  run_charge_migration=False, run_charge_migration_ft=False,
                  run_spectrum_reconstruction=False,
                  field_file_help=False, molcas_input_help=False,
                  lus=False, gridit: bool or str = True, write_charge_migration=None,
                  debug_mode=False, justh5=False, justgetdipoles=False, justgetdensity=False,
                  weights_file=None, givenfieldfile=None,
                   make_fortran=False, make_fortran_config=None,

                  plot=True  # todo usage to be changed, here for testing

                  ):
    if make_fortran:
        from densitypy.project_utils.fortran_compilation_handler import compile_ifort_fortran_code
        return_code, compiler_output_df = compile_ifort_fortran_code(**make_fortran_config)

        # todo either way should postprocess output for better user experience
        if return_code != 0:

            if compiler_output_df is not None:
                logger.error(f"Compilation failed with return code {return_code}.")
                # compiler_output_df = compiler_output_df.sort_values(by='Type',
                #                                                     key=lambda x: x.map(DEFAULT_LOGGIN_VALUES))
                # logger.(f'\n{compiler_output_df}')
            exit(return_code)

    # TODO too explicit, does not allow for new args to be added easily
    COMMAND_ARGS = [molcas_input, lus, run_charge_migration, run_charge_migration_ft,
                    molcas_input_help, justh5, justgetdensity, justgetdipoles, field_file_help,
                    run_spectrum_reconstruction, plot]

    # Change directory to study_directory and back to original directory when done or if crashes
    with change_directory_manager(study_directory):
        # Load Values # todo fixup the json_config and runnable command behavior, is ugly
        json_config = parse_configuration_file(
            json_config_path)  # If not passed or does not exist will write example file
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
        xyz_geometry_path = project_settings['xyzmoleculegeometry']
        number_of_states = project_settings['numberofstates']
        list_of_orbitals = project_settings['listofactiveorbitals']
        molcas_output_directory = project_settings['molcasoutputdirectory']
        experiment_directory = project_settings['experimentdirectory']

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
        time_Delay_weight_factor = probe_settings['timedelayweightfactor']
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

        time_delay_range = generate_time_delays(number_of_pp, time_delay_start, 0.0, time_delay_stop,
                                                time_Delay_weight_factor)

        # Useful Flags for e.g. debugging
        if justh5:
            load_project_rasscf_h5_file(project_name, molcas_output_directory, molcas_output_directory)
            exit()
        if justgetdensity:
            grid_file_data_dict = parse_project_grid_file(project_name, molcas_output_directory)
            write_grid_density_file(grid_file_data_dict, molcas_output_directory)
            exit()

        if justgetdipoles:
            project = DipolesLogParser(log_path=f'{molcas_output_directory}/{project_name}.log')
            project.write_csvs(directory=molcas_output_directory)
            fig = project.get_mu_heatmaps()
            fig.savefig(f'{molcas_output_directory}/dipole-heatmap.png')
            exit()

        if field_file_help:
            Write_FieldHelp()

        # >OpenMolcas
        if molcas_input:

            n_points = 0
            # Prepare Input
            make_directory(molcas_output_directory)
            #
            # Copies input file and returns keywords in input file #todo should actually parse input , use pymolcas
            keywords_needed_found = copy_and_parse_molcas_input_file_to_edit(pymolcas_input=molcas_input,
                                                                             project_name=project_name,
                                                                             molcas_directory=molcas_output_directory)

            copy_file_to(xyz_geometry_path, molcas_output_directory)

            # gridit can be fixme deprecated until updated to use new grid algorithm
            #   True = use gridit on molcas with with uniform grid points
            #   False = Do not use or create any grid or add gridit to input file
            #   'smart' =  Use our smart grid algorithm to create a grid based on the molecule and step size
            #   'limited' = Use our truncated grid algorithm to create a grid based on the molecule and config settings
            #   'manual' = Use a grid file to create a grid, this is the only option that requires a grid file to be specified in the config file
            if gridit in ['smart', 'limited']:
                limitedgrid = True if gridit == 'limited' else False
                n_points = make_better_grid(molcas_output_directory, xyz_geometry_path, step_size, boundary,
                                            limitedgrid)
            elif gridit == 'manual':
                gridfile = json_config['grid_settings']['default_grid_file']
                system(f"cp {gridfile} {molcas_output_directory}/gridcoord")
                n_points = file_lenth(f'{molcas_output_directory}/gridcoord')
                logger.info(f"Number of Points = {n_points}")
            elif gridit:
                n_points = make_grid_coordinates(molcas_output_directory, nx, ny, nz, xmin, xmax, ymin, ymax, zmin,
                                                 zmax)
            else:
                assert gridit is False, f"gridit is {gridit} but should be False, True, 'smart', 'limited', or 'manual'"

            if gridit:
                add_grid_it_to_manual_input_file(molcas_input, project_name, molcas_output_directory, list_of_orbitals,
                                                 n_points)
            # >Run Open Molcas
            call_open_molcas(project_name, molcas_output_directory)
            # Extract Data Required for Charge Migration
            copy_file_to(molcas_input, molcas_output_directory + "/")
            logfilepath = find(project_name + ".log", ".", molcas_output_directory)
            copy_file_to(logfilepath + "/" + project_name + ".log ", molcas_output_directory + "/")
            copy_file_to(xyz_geometry_path, molcas_output_directory + "/")

            if "RASSI" in keywords_needed_found:  # todo add to documentation that dipoles are only parsed for RASSI
                project = DipolesLogParser(log_path=f'{molcas_output_directory}/{project_name}.log')
                project.write_csvs(directory=molcas_output_directory)
                fig = project.get_mu_heatmaps()
                fig.savefig(f'{molcas_output_directory}/dipole-heatmap.png')

            if gridit:
                grid_file_data_dict = parse_project_grid_file(project_name, molcas_output_directory)
                write_grid_density_file(grid_file_data_dict, molcas_output_directory)

            if "RASSCF" in keywords_needed_found:
                load_project_rasscf_h5_file(project_name, molcas_output_directory, molcas_output_directory)

        # >ChargeMigration TODO(Subject to Change based on C_i reconstruction implementation)
        if run_charge_migration:
            if not path.exists(molcas_output_directory):
                logger.error(f"Could not find {molcas_output_directory}")
                exit()
            make_directory(experiment_directory)

            if givenfieldfile:
                field_file = givenfieldfile
            else:
                Write_Pulses(f"{experiment_directory}/{field_file}", type_of_pulse_pump, start_time,
                             pump_central_frequency,
                             pump_periods, pump_phase, pump_intensity, pump_polarization,
                             type_of_pulse_probe, time_delay_range, probe_central_frequency,
                             probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration)

            # Run Charge Migration TODO(Subject to Change based on C_i reconstruction implementation)
            Call_Charge_Migration(DEFAULT_BIN_FILE_PATH, molcas_output_directory, experiment_directory, number_of_times,
                                  min_time, max_time, f"{experiment_directory}/{field_file}", ft_time_step,
                                  ft_width_step, f"{molcas_output_directory}/{xyz_geometry_path}", list_of_orbitals,
                                  write_charge_migration, Volume, debug_mode, weights_file, dephasing_factor,
                                  relaxation_factor, bath_temperature)

            copy_file_to(json_config_path, f"{experiment_directory}")

        # >ChargeMigrationFT
        if run_charge_migration_ft:
            if givenfieldfile:
                field_file = givenfieldfile
            else:
                # Make Field File
                Write_Pulses(f"{experiment_directory}/{field_file}", type_of_pulse_pump, start_time,
                             pump_central_frequency,
                             pump_periods, pump_phase, pump_intensity, pump_polarization,
                             type_of_pulse_probe, time_delay_range, probe_central_frequency,
                             probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration)
            # Run Charge Migration FT
            Call_Charge_MigrationFT(DEFAULT_BIN_FILE_PATH, molcas_output_directory, experiment_directory,
                                    f"{molcas_output_directory}/{xyz_geometry_path}",
                                    number_of_omegas,
                                    min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                    ft_time_step, ft_width_step, f"{experiment_directory}/{field_file}",
                                    debug_mode)
        # >SpectrumReconstruction
        if run_spectrum_reconstruction:
            Call_Spectrum_Reconstruction_n_Difference(DEFAULT_BIN_FILE_PATH, molcas_output_directory,
                                                      experiment_directory,
                                                      f"{molcas_output_directory}/{xyz_geometry_path}",
                                                      number_of_omegas,
                                                      min_omegas, max_omegas,
                                                      number_of_tau_omegas, min_tau_omega, max_tau_omega, debug_mode)

        if plot:
            # Becke Weights
            # plot_becke_weights(study_directory, experiment_directory, xyz_geometry_path, weights_file)

            # Lets plot the Pulses
            plot_pulses(study_directory, experiment_directory, time_delay_range, min_time, max_time, plot_all=False)
            plot_ft_pulses(study_directory, experiment_directory, time_delay_range, plot_all=False)

            # Lets plot the Dipolar Reponse vs Time (t)
            plot_dipoles_v_time(study_directory, experiment_directory, time_delay_range, min_time, max_time,
                                plot_all=False)
            plot_atomic_dipoles_v_time(study_directory, experiment_directory, time_delay_range, min_time, max_time,
                                       xyz_geometry_path, plot_all=False)
            difference_between_dipole_and_atomic_charges_v_time(study_directory, experiment_directory, time_delay_range,
                                                                min_time, max_time, xyz_geometry_path)

            # Lets plot the Dipolar Reponse vs Time (t) in the Frequency Domain (w)
            plot_ft_all_dipoles_v_time(study_directory, experiment_directory, dephasing_factor, relaxation_factor,
                                       pump_settings, probe_settings, charge_migration_ft_settings,
                                       try_color_maps=False)
            plot_ft_all_atomic_dipoles_v_time(study_directory, experiment_directory, dephasing_factor,
                                              relaxation_factor,
                                              pump_settings, probe_settings, charge_migration_ft_settings,
                                              xyz_geometry_path, try_color_maps=False)

            # Lets plot the 2D Spectrum
            plot_2d_spectrum(study_directory, experiment_directory, dephasing_factor, relaxation_factor,
                             pump_settings, probe_settings, charge_migration_ft_settings)
            plot_2d_spectrum_peak_analysis(study_directory, experiment_directory)
            plot_2d_spectrum_interactive(study_directory, experiment_directory) # not working, shifts the axis


if __name__ == "__main__":
    make_fortran_config = dict(directory='/home/ruben/PycharmProjects/DensityPy/densityfort',
                               make_flags='all DEB_FLAG=d',
                               # ... other settings but this is mvp
                               )

    run_densitypy(json_config_path="xy_polarized_config.json",
                  study_directory="/home/ruben/PycharmProjects/DensityPy/Studies/ExampleStudy",
                  molcas_input=False,  # 'molcas_input_help.input', # False just means not running molcas
                  run_charge_migration=False,
                  run_charge_migration_ft=False,
                  run_spectrum_reconstruction=False,
                  plot=True,
                  #
                  field_file_help=False, molcas_input_help=False,
                  lus=False, gridit=True, write_charge_migration=None, debug_mode=False,
                  justh5=False, justgetdipoles=False, justgetdensity=False, weights_file=None, givenfieldfile=None,
                  make_fortran=False, make_fortran_config=make_fortran_config
                  )

# counting them down
# fixme most of the data is gotten from RASSCF H5 file only, but I believe we are looking for the integrals, RASSI?
#   -fix manipulation of Atomic Charges as they require a flipped sign or rotation of matrix e.g. [::-1] to correctly
#       match the dipole
#   -implement more  validation checks of settings between runs, between modules, and between configurations of
#       different programs like molcas etc...
#   -update documentation
#   -add more tests
#   -update gridit to use the new grid algorithm
