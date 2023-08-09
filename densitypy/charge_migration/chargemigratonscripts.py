#! /usr/bin/env python3.6
# >import chargemigratonscripts as rcms
from os import path, system

from densitypy.project_utils.def_functions import uniquify, execute_command, delete_files_or_directories

from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])



def Write_FieldHelp():
    # >Creates helpful field file guide for Pulse3d module
    with open('Field_pulses_template', 'w') as fout:
        fout.write("#"
                   "\n# A Pulse sequence is initiated with its name, which is a string within square parentheses."
                   "\n# e.g., [XUV]. The Pulse sequence is determined by a list of either names of other pulses "
                   "sequences, and/or any number of individual pulses, separated by semicolon \";\", and the last "
                   "element in the sequence must also end with a semicolon. The individual pulse is idenfied "
                   "by round parentheses, which contain:"
                   "\n# Type of pulse: G for Gaussian, C for cosine square"
                   "\n# central time, in atomic units# central frequency, in atomic units"
                   "\n# number of full periods in the time intervals where the pulse envelope is above half of its "
                   "maxiumum (full-width at half maximum)"
                   "\n# carrier-envelope phase# average intensity of the field if it were with constant envelope equal"
                   " to the maximum, measured in peta-Watt per square cm  (peta = 10^15)"
                   "\n# Optionally, one can indicate the polarization with its spherical angles"
                   "\n# theta (in degrees)"
                   "\n# phi (in degrees)"
                   "\n# e.g.  z => 0 0, x => 90 0, y => 90 90"
                   "\n\n[XUV]{( G  0.000 0.30 5 180  0.001 90 0 );}"
                   "\n[PP010]{ XUV; ( G 25.00 0.18 5 0 0.001 90 0 ); }"
                   "\nEXECUTE{PP010;}")
    logger.info("Field_pulses_template created")
    exit()


def generate_time_delays(number_of_pp, time_delay_start, time_delay_stop):
    if number_of_pp == 0:
        raise ValueError("Number of points must be greater than 0.")

    if number_of_pp == 1:
        time_delay_range = [0.0]
    elif number_of_pp % 2 == 0:  # even case
        change_in_delay = (time_delay_stop - time_delay_start) / (number_of_pp // 2)
        half_point = number_of_pp // 2
        time_delay_range = list(range(-half_point, half_point))
        time_delay_range = [x * change_in_delay for x in time_delay_range]
    else:  # odd case
        change_in_delay = (time_delay_stop - time_delay_start) / ((number_of_pp - 1) // 2)
        half_point = (number_of_pp - 1) // 2
        time_delay_range = list(range(-half_point, half_point + 1))
        time_delay_range = [x * change_in_delay for x in time_delay_range]

    return time_delay_range
# >Writes pulses to Field File ,Pump Probe Experiment (0)
def Write_Pulses(field_file, type_of_pulse_pump, start_time, pump_central_frequency,
                 pump_periods, pump_phase, pump_intensity, pump_polarization,
                 type_of_pulse_probe, time_delay_range, probe_central_frequency,
                 probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration):

    print(f'{field_file = }')

    with open(field_file, 'w') as fout:
        pump_polarization = " ".join(map(str, pump_polarization))
        probe_polarization = " ".join(map(str, probe_polarization))

        # >Writes Pump Pulse
        fout.write("[XUV]{( " + str(type_of_pulse_pump) + " " +
                   str(start_time) + " " +
                   str(pump_central_frequency) + " " +
                   str(pump_periods) + " " +
                   str(pump_phase) + "  " +
                   str(pump_intensity) + " " +
                   str(pump_polarization) + " );}")
        for Time_Delay in time_delay_range:
            # split probe_polarization from list to string of the form "theta phi"

            # >Writes Probe Pulse
            # fout.write("\n[PP" + str(Time_Delay) + "]{ XUV; ( " + str(type_of_pulse_probe) + " " +
            #            str(Time_Delay) + " " +
            #            str(probe_central_frequency) + " " +
            #            str(probe_periods) + " " + str(probe_phase) + " " +
            #            str(probe_intensity) + " " + " " + probe_polarization + " );}")
            fout.write(f"\n[PP{Time_Delay}]{{ XUV; ( {type_of_pulse_probe} {Time_Delay} {probe_central_frequency} "
                       f"{probe_periods} {probe_phase} {probe_intensity} {probe_polarization} );}}")

        fout.write("\nEXECUTE{")
        for Time_Delay in time_delay_range:
            fout.write("PP" + str(Time_Delay) + "; ")
        if write_charge_migration:
            fout.write("}")
        else:
            fout.write("XUV;}")


# # >THIS IS AN EDIT IN PROGRESS DO NOT USE< USE THE ABOVE
# def Write_Pulses(field_file, type_of_pulse_pump, start_time, pump_central_frequency,
#                  pump_periods, pump_phase, pump_intensity, pump_polarization,
#                  type_of_pulse_probe, time_delay_range, probe_central_frequency,
#                  probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration):
#     with open(field_file, 'w') as fout:
#         for Time_Delay in time_delay_range:
#             # >Writes Probe Pulse
#             fout.write("\n[PP" + str(Time_Delay) + "]{( " + str(type_of_pulse_probe) + " " +
#                        str(Time_Delay) + " " +
#                        str(probe_central_frequency) + " " +
#                        str(probe_periods) + " " + str(probe_phase) + " " +
#                        str(probe_intensity) + " " + " " + probe_polarization + " );}")
#         fout.write("\nEXECUTE{")
#         for Time_Delay in time_delay_range:
#             fout.write("PP" + str(Time_Delay) + "; ")
#         else:
#             fout.write("}")


# >Calls Charge Migration Code and gives command line arguements (1)
def Call_Charge_Migration(Bin_Directory, input_directory, experiment_directory, number_of_times,
                          min_time, max_time, field_file, stept, stepw,
                          geometry, orbital_list, write_charge_migrationflag, Volume, debug_mode, weights_file, dephasing_factor,
                          relaxation_factor, bath_temp, i_excitation, i_epsilon):

    weights_file_decoy = ""
    if weights_file:
        weights_file_decoy = "-w " + weights_file

    write_charge_migrationflag_decoy = ""
    if write_charge_migrationflag:
        write_charge_migrationflag_decoy = "-sden"

    screen_print = "Calling ChargeMigration"
    debug_mode_decoy = ""
    if debug_mode:
        debug_mode_decoy = "d"
        screen_print = "Calling ChargeMigration in Debug Mode (DEB_FLAG=d)"

    logger.info(screen_print)
    # >Runs Charge Migration (Fortran Code)
    execute_command(
        # FIXME : This is a temporary fix for the issue with the MKL library and my pycharm env
        f"export LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2023.2.0/lib/intel64:/opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH' &&"
        f"{Bin_Directory}/"
        f"{debug_mode_decoy}ChargeMigration -i {input_directory} -o {experiment_directory} -nt {str(number_of_times)} "
        f"-tmin {str(min_time)} -tmax {str(max_time)} -field {field_file} -vol {str(Volume)} -stept {str(stept)} "
        f"-stepw {str(stepw)} -xyz {str(geometry)} {weights_file_decoy} {write_charge_migrationflag_decoy} -iorb "
        f"{','.join(map(str, orbital_list))} -rf {str(relaxation_factor)} -bath {str(bath_temp)} -df "
        f"{str(dephasing_factor)} -s {i_excitation} -e {i_epsilon}"
    , _logger=logger)

    logger.info("Finished Executing ChargeMigration")


def Call_Charge_MigrationFT(Bin_Directory, input_directory, experiment_directory, geometry, Number_of_Omegas,
                            min_omegas, max_omegas, number_of_tau_omega, min_tau_omega, max_tau_omega,
                            ft_time_step, ft_width_step, field_file, debug_mode, i_excitation, i_epsilon):
    screen_print = "Running ChargeMigrationFT"
    debug_mode_decoy = ""
    if debug_mode:
        debug_mode_decoy = "d"
        screen_print = "Running ChargeMigrationFT in Debug Mode (DEB_FLAG=d)"

    logger.info(screen_print)
    # >Runs Charge Migration FT (Fortran Code)
    execute_command(
        # FIXME : This is a temporary fix for the issue with the MKL library and pycharm
        f"export LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2023.2.0/lib/intel64:/opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH' &&"
        f"{Bin_Directory}/"
        f'{debug_mode_decoy}ChargeMigrationFT -i {str(input_directory)} -o {str(experiment_directory)} -xyz '
        f'{str(geometry)} -stept {str(ft_time_step)} -stepw {str(ft_width_step)} -field {str(field_file)} -nw '
        f'{str(Number_of_Omegas)} -wmax {str(min_omegas)} -wmin {str(max_omegas)} -ntw {str(number_of_tau_omega)} '
        f'-twmax {str(min_tau_omega)} -twmin {str(max_tau_omega)} -s {i_excitation} -e {i_epsilon}'
    , _logger=logger)

    logger.info("Finished Executing ChargeMigrationFT")


def Save_Spectrum_Difference(Bin_Directory, experiment_directory, difference_file, dephasing_factor, relaxation_factor,
                             time_delay_start,
                             time_delay_stop, min_omegas,
                             max_omegas, pump_periods, probe_periods, pump_intensity,
                             probe_intensity, Number_Of_PP, pump_phase, probe_phase, ft_time_step, ft_width_step,
                             pump_polarization, probe_polarization, Number_of_Omegas, min_tau_omega,
                             max_tau_omega):
    FT_WW_TAIL = f'DephasingFactor_{dephasing_factor}_' \
                 f'RelaxingFactor_{relaxation_factor}_Tau_{time_delay_start}_{time_delay_stop}_' \
                 f'W_{min_omegas}_{max_omegas}_PumpPeriods_{pump_periods}_ProbePeriods_{probe_periods}_' \
                 f'PumpIntensity_{pump_intensity}_ProbeIntensity_{probe_intensity}_' \
                 f'PumpPhase_{pump_phase}_ProbePhase_{probe_phase}_ft_time_step_{ft_time_step}_ft_width_step_{ft_width_step}'

    oldpath = f'{experiment_directory}/DipoleFT_ww_reconstructed'
    if path.exists(oldpath):
        #
        newfilepath = f'{experiment_directory}/Dipole/DipoleFT_ww_reconstructed_{FT_WW_TAIL}'
        newfilepath = uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')

    oldpath = f'{difference_file}'
    if path.exists(oldpath):
        #
        newfilepath = f'{difference_file}_{FT_WW_TAIL}'
        newfilepath = uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')


def Dipole_Charge_Comparison(dipole_file, charge_file, output_file):
    delete_files_or_directories()
    system('cat ' + str(
        dipole_file) + ' | awk \'{print $1\" \"$2\" \"$3**2+$4**2\" \"$5**2+$6**2\" \"$7**2+$8**2}\' > temp_dipole')
    system('cat ' + str(charge_file) + ' | awk \'{print $3**2+$4**2\" \"$5**2+$6**2\" \"$7**2+$8**2}\' > temp_charge')
    system(
        'paste temp_dipole temp_charge| awk \'{print $1\" \"$2\" \"$3\" \"$4\" \"$5\" \"$6\" \"$7\" \"$8\" \"$3-$6\" \"$4-$7\" \"$5-$8}\' > temp_difference')

    system("perl -pe \"s/0 0 0 0 0 0   0 0 0/ /g\" temp_difference>" + str(output_file))
    system('rm temp_charge temp_dipole temp_difference')


def Call_Spectrum_Reconstruction_n_Difference(Bin_Directory, molcas_output_directory, experiment_directory, xyz_file_path, Number_of_Omegas,
                                              min_omegas, max_omegas, number_of_tau_omega, min_tau_omega, max_tau_omega,
                                              debug_mode):
    screen_print = "Running Spectrum Reconstruction"
    debug_mode_decoy = ""
    if debug_mode:
        debug_mode_decoy = "d"
        screen_print = "Running Spectrum Reconstruction in Debug Mode (DEB_FLAG=d)"
    # Momentary , needs work

    logger.info("Running SpectrumReconstruction ...")
    execute_command(
        # FIXME : This is a temporary fix for the issue with the MKL library and pycharm
        f"export LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2023.2.0/lib/intel64:/opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH' &&"
        f"{Bin_Directory}/"
        f'{debug_mode_decoy}SpectrumReconstruction -i {str(molcas_output_directory)} -o {str(experiment_directory)} -xyz '
        f'{str(xyz_file_path)} -nw {str(Number_of_Omegas)} -wmax {str(max_omegas)} -wmin {str(min_omegas)} -ntw {str(number_of_tau_omega)} -twmax {str(max_tau_omega)} -twmin {str(min_tau_omega)}'
        , _logger=logger
    )

    logger.info("Finished Executing SpectrumReconstruction")


def Save_Previous_FT(experiment_directory, dephasing_factor, relaxation_factor, time_delay_start, time_delay_stop, min_omegas,
                     max_omegas, pump_periods, probe_periods, pump_intensity,
                     probe_intensity, Number_Of_PP, pump_phase, probe_phase, ft_time_step, ft_width_step,
                     pump_polarization, probe_polarization, Number_of_Omegas, min_tau_omega,
                     max_tau_omega):
    #
    FT_WW_TAIL = f'DephasingFactor_{dephasing_factor}_' \
                 f'RelaxingFactor_{relaxation_factor}_Tau_{time_delay_start}_{time_delay_stop}_' \
                 f'W_{min_omegas}_{max_omegas}_PumpPeriods_{pump_periods}_ProbePeriods_{probe_periods}_' \
                 f'PumpIntensity_{pump_intensity}_ProbeIntensity_{probe_intensity}_' \
                 f'PumpPhase_{pump_phase}_ProbePhase_{probe_phase}_ft_time_step_{ft_time_step}_ft_width_step_{ft_width_step}'
    #
    FT_ALL_TAIL = f'DephasingFactor_{dephasing_factor}_' \
                  f'RelaxingFactor_{relaxation_factor}_tau_{time_delay_start}_{time_delay_stop}_' \
                  f'w_{min_omegas}_{max_omegas}_PumpPeriods_{pump_periods}_ProbePeriods_{probe_periods}_' \
                  f'PumpIntensity_{pump_intensity}_ProbeIntensity_{probe_intensity}_' \
                  f'PumpPhase_{pump_phase}_ProbePhase_{probe_phase}_ft_time_step_{ft_time_step}_ft_width_step_{ft_width_step}'
    #
    oldpath = f'{experiment_directory}/AtomicCharge/AtomicChargeFT_ww'
    if path.exists(oldpath):
        #
        newfilepath = f'{experiment_directory}/AtomicCharge/AtomicChargeFT_ww_{FT_WW_TAIL}'
        newfilepath = uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
    #
    oldpath = f'{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL'
    if path.exists(oldpath):
        #
        newfilepath = f'{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL_{FT_ALL_TAIL}'
        newfilepath = uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
    #
    oldpath = f'{experiment_directory}/Dipole/DipoleFT_ww'
    if path.exists(oldpath):
        #
        newfilepath = f'{experiment_directory}/Dipole/DipoleFT_ww_{FT_WW_TAIL}'
        newfilepath = uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
    #
    oldpath = f'{experiment_directory}/Dipole/DipoleFT_ALL'
    if path.exists(oldpath):
        #
        newfilepath = f'{experiment_directory}/Dipole/DipoleFT_ALL_{FT_ALL_TAIL}'
        newfilepath = uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
