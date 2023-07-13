#! /usr/bin/env python3.6
# >import chargemigratonscripts as rcms
import def_functions as rfunc
import sys
from os import path, system


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
    print("Field_pulses_template created")
    sys.exit()


# >Writes pulses to Field File ,Pump Probe Experiment (0)
def Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
                 Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
                 Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
                 Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization, writeCM):
    with open(Field_File, 'w') as fout:
        # >Writes Pump Pulse
        fout.write("[XUV]{( " + str(Type_of_Pulse_Pump) + " " +
                   str(Start_Time) + " " +
                   str(Pump_Central_Frequency) + " " +
                   str(Pump_Periods) + " " +
                   str(Pump_Phase) + "  " +
                   str(Pump_Intensity) + " " +
                   str(Pump_Polarization) + " );}")
        for Time_Delay in Time_Delay_Range:
            # >Writes Probe Pulse
            fout.write("\n[PP" + str(Time_Delay) + "]{ XUV; ( " + str(Type_of_Pulse_Probe) + " " +
                       str(Time_Delay) + " " +
                       str(Probe_Central_Frequency) + " " +
                       str(Probe_Periods) + " " + str(Probe_Phase) + " " +
                       str(Probe_Intensity) + " " + " " + Probe_Polarization + " );}")
        fout.write("\nEXECUTE{")
        for Time_Delay in Time_Delay_Range:
            fout.write("PP" + str(Time_Delay) + "; ")
        if writeCM:
            fout.write("}")
        else:
            fout.write("XUV;}")


# # >THIS IS AN EDIT IN PROGRESS DO NOT USE< USE THE ABOVE
# def Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
#                  Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
#                  Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
#                  Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization, writeCM):
#     with open(Field_File, 'w') as fout:
#         for Time_Delay in Time_Delay_Range:
#             # >Writes Probe Pulse
#             fout.write("\n[PP" + str(Time_Delay) + "]{( " + str(Type_of_Pulse_Probe) + " " +
#                        str(Time_Delay) + " " +
#                        str(Probe_Central_Frequency) + " " +
#                        str(Probe_Periods) + " " + str(Probe_Phase) + " " +
#                        str(Probe_Intensity) + " " + " " + Probe_Polarization + " );}")
#         fout.write("\nEXECUTE{")
#         for Time_Delay in Time_Delay_Range:
#             fout.write("PP" + str(Time_Delay) + "; ")
#         else:
#             fout.write("}")


# >Calls Charge Migration Code and gives command line arguements (1)
def Call_Charge_Migration(Input_Directory, Output_Directory, Number_Of_Times,
                          Min_Time, Max_Time, Field_File, stept, stepw,
                          geometry, orbital_list, writeCMflag, Volume, debug_mode, weights_file, dephasing_factor,
                          relaxing_factor, bath_temp, iExcitation, iEPSILON):
    weights_file_decoy = ""
    if weights_file:
        weights_file_decoy = "-w " + weights_file

    writeCMflag_decoy = ""
    if writeCMflag:
        writeCMflag_decoy = "-sden"

    screen_print = "Running ChargeMigration"
    debug_mode_decoy = ""
    if debug_mode:
        debug_mode_decoy = "d"
        screen_print = "Running ChargeMigration in Debug Mode (DEB_FLAG=d)"

    print(screen_print)
    # >Runs Charge Migration (Fortran Code)
    rfunc.Execute(
        f"{debug_mode_decoy}ChargeMigration -i {Input_Directory} -o {Output_Directory} -nt {str(Number_Of_Times)} "
        f"-tmin {str(Min_Time)} -tmax {str(Max_Time)} -field {Field_File} -vol {str(Volume)} -stept {str(stept)} "
        f"-stepw {str(stepw)} -xyz {str(geometry)} {weights_file_decoy} {writeCMflag_decoy} -iorb "
        f"{','.join(map(repr, orbital_list))} -rf {str(relaxing_factor)} -bath {str(bath_temp)} -df "
        f"{str(dephasing_factor)} -s {iExcitation} -e {iEPSILON}")

    print(f"\n{debug_mode_decoy}ChargeMigration -i {Input_Directory} -o {Output_Directory} -nt {str(Number_Of_Times)} "
          f"-tmin {str(Min_Time)} -tmax {str(Max_Time)} -field {Field_File} -vol {str(Volume)} -stept {str(stept)} "
          f"-stepw {str(stepw)} -xyz {str(geometry)} {weights_file_decoy} {writeCMflag_decoy} -iorb "
          f"{','.join(map(repr, orbital_list))} -rf {str(relaxing_factor)} -bath {str(bath_temp)} -df "
          f"{str(dephasing_factor)} -s {iExcitation} -e {iEPSILON}")


def Call_Charge_MigrationFT(Input_Directory, Output_Directory, geometry, Number_of_Omegas,
                            Min_Omegas, Max_Omegas, Number_of_TauOmega, Min_TauOmega, Max_TauOmega,
                            TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON):
    screen_print = "Running ChargeMigrationFT"
    debug_mode_decoy = ""
    if debug_mode:
        debug_mode_decoy = "d"
        screen_print = "Running ChargeMigrationFT in Debug Mode (DEB_FLAG=d)"

    print(screen_print)
    # >Runs Charge Migration FT (Fortran Code)
    rfunc.Execute(
        f'{debug_mode_decoy}ChargeMigrationFT -i {str(Input_Directory)} -o {str(Output_Directory)} -xyz '
        f'{str(geometry)} -stept {str(TimeStep_FT)} -stepw {str(WidthStep_FT)} -field {str(Field_File)} -nw '
        f'{str(Number_of_Omegas)} -wmax {str(Min_Omegas)} -wmin {str(Max_Omegas)} -ntw {str(Number_of_TauOmega)} '
        f'-twmax {str(Min_TauOmega)} -twmin {str(Max_TauOmega)} -s {iExcitation} -e {iEPSILON}')

    print(f'\n{debug_mode_decoy}ChargeMigrationFT -i {str(Input_Directory)} -o {str(Output_Directory)} -xyz '
          f'{str(geometry)} -stept {str(TimeStep_FT)} -stepw {str(WidthStep_FT)} -field {str(Field_File)} -nw '
          f'{str(Number_of_Omegas)} -wmax {str(Min_Omegas)} -wmin {str(Max_Omegas)} -ntw {str(Number_of_TauOmega)} '
          f'-twmax {str(Min_TauOmega)} -twmin {str(Max_TauOmega)} -s {iExcitation} -e {iEPSILON}')


def Save_Spectrum_Difference(Output_Directory, difference_file, Dephasing_Factor, Relaxing_Factor, Time_Delay_Start,
                             Time_Delay_Stop, Min_Omegas,
                             Max_Omegas, Pump_Periods, Probe_Periods, Pump_Intensity,
                             Probe_Intensity, Number_Of_PP, Pump_Phase, Probe_Phase, TimeStep_FT, WidthStep_FT,
                             Pump_Polarization, Probe_Polarization, Number_of_Omegas, Min_TauOmega,
                             Max_TauOmega):
    FT_WW_TAIL = f'DephasingFactor_{Dephasing_Factor}_' \
                 f'RelaxingFactor_{Relaxing_Factor}_Tau_{Time_Delay_Start}_{Time_Delay_Stop}_' \
                 f'W_{Min_Omegas}_{Max_Omegas}_PumpPeriods_{Pump_Periods}_ProbePeriods_{Probe_Periods}_' \
                 f'PumpIntensity_{Pump_Intensity}_ProbeIntensity_{Probe_Intensity}_' \
                 f'PumpPhase_{Pump_Phase}_ProbePhase_{Probe_Phase}_TimeStep_FT_{TimeStep_FT}_WidthStep_FT_{WidthStep_FT}'

    oldpath = f'{Output_Directory}/DipoleFT_ww_reconstructed'
    if path.exists(oldpath):
        #
        newfilepath = f'{Output_Directory}/Dipole/DipoleFT_ww_reconstructed_{FT_WW_TAIL}'
        newfilepath = rfunc.uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')

    oldpath = f'{difference_file}'
    if path.exists(oldpath):
        #
        newfilepath = f'{difference_file}_{FT_WW_TAIL}'
        newfilepath = rfunc.uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')


def Dipole_Charge_Comparison(dipole_file, charge_file, output_file):
    system('rm ' + str(output_file))
    system('cat ' + str(
        dipole_file) + ' | awk \'{print $1\" \"$2\" \"$3**2+$4**2\" \"$5**2+$6**2\" \"$7**2+$8**2}\' > temp_dipole')
    system('cat ' + str(charge_file) + ' | awk \'{print $3**2+$4**2\" \"$5**2+$6**2\" \"$7**2+$8**2}\' > temp_charge')
    system(
        'paste temp_dipole temp_charge| awk \'{print $1\" \"$2\" \"$3\" \"$4\" \"$5\" \"$6\" \"$7\" \"$8\" \"$3-$6\" \"$4-$7\" \"$5-$8}\' > temp_difference')

    system("perl -pe \"s/0 0 0 0 0 0   0 0 0/ /g\" temp_difference>" + str(output_file))
    system('rm temp_charge temp_dipole temp_difference')


def Call_Spectrum_Reconstruction_n_Difference(MolcasOut, SimOut, Geometry, Number_of_Omegas,
                                              Min_Omegas, Max_Omegas, Number_of_TauOmega, Min_TauOmega, Max_TauOmega,
                                              debug_mode):
    screen_print = "Running Spectrum Reconstruction"
    debug_mode_decoy = ""
    if debug_mode:
        debug_mode_decoy = "d"
        screen_print = "Running Spectrum Reconstruction in Debug Mode (DEB_FLAG=d)"
    # Momentary , needs work

    print("Running SpectrumReconstruction ...")
    system(
        debug_mode_decoy + 'SpectrumReconstruction -i ' + MolcasOut + ' -o ' + SimOut + ' -xyz ' + Geometry +
        ' -nw ' + str(Number_of_Omegas) + ' -wmax ' + str(Max_Omegas) + ' -wmin ' + str(Min_Omegas) + ' -ntw ' + str(
            Number_of_TauOmega) + ' -twmax ' + str(Max_TauOmega) + ' -twmin ' + str(Min_TauOmega))

    print("Creating difference_" + SimOut + " to compare the Dipole and Charge Spectra")
    Dipole_Charge_Comparison(SimOut + '/Dipole/DipoleFT_ww', SimOut +
                             '/Dipole/DipoleFT_ww_reconstructed', 'difference_' + SimOut)
    print("Done")


def Save_Previous_FT(Output_Directory, Dephasing_Factor, Relaxing_Factor, Time_Delay_Start, Time_Delay_Stop, Min_Omegas,
                     Max_Omegas, Pump_Periods, Probe_Periods, Pump_Intensity,
                     Probe_Intensity, Number_Of_PP, Pump_Phase, Probe_Phase, TimeStep_FT, WidthStep_FT,
                     Pump_Polarization, Probe_Polarization, Number_of_Omegas, Min_TauOmega,
                     Max_TauOmega):
    #
    FT_WW_TAIL = f'DephasingFactor_{Dephasing_Factor}_' \
                 f'RelaxingFactor_{Relaxing_Factor}_Tau_{Time_Delay_Start}_{Time_Delay_Stop}_' \
                 f'W_{Min_Omegas}_{Max_Omegas}_PumpPeriods_{Pump_Periods}_ProbePeriods_{Probe_Periods}_' \
                 f'PumpIntensity_{Pump_Intensity}_ProbeIntensity_{Probe_Intensity}_' \
                 f'PumpPhase_{Pump_Phase}_ProbePhase_{Probe_Phase}_TimeStep_FT_{TimeStep_FT}_WidthStep_FT_{WidthStep_FT}'
    #
    FT_ALL_TAIL = f'DephasingFactor_{Dephasing_Factor}_' \
                  f'RelaxingFactor_{Relaxing_Factor}_tau_{Time_Delay_Start}_{Time_Delay_Stop}_' \
                  f'w_{Min_Omegas}_{Max_Omegas}_PumpPeriods_{Pump_Periods}_ProbePeriods_{Probe_Periods}_' \
                  f'PumpIntensity_{Pump_Intensity}_ProbeIntensity_{Probe_Intensity}_' \
                  f'PumpPhase_{Pump_Phase}_ProbePhase_{Probe_Phase}_TimeStep_FT_{TimeStep_FT}_WidthStep_FT_{WidthStep_FT}'
    #
    oldpath = f'{Output_Directory}/AtomicCharge/AtomicChargeFT_ww'
    if path.exists(oldpath):
        #
        newfilepath = f'{Output_Directory}/AtomicCharge/AtomicChargeFT_ww_{FT_WW_TAIL}'
        newfilepath = rfunc.uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
    #
    oldpath = f'{Output_Directory}/AtomicCharge/AtomicChargeFT_ALL'
    if path.exists(oldpath):
        #
        newfilepath = f'{Output_Directory}/AtomicCharge/AtomicChargeFT_ALL_{FT_ALL_TAIL}'
        newfilepath = rfunc.uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
    #
    oldpath = f'{Output_Directory}/Dipole/DipoleFT_ww'
    if path.exists(oldpath):
        #
        newfilepath = f'{Output_Directory}/Dipole/DipoleFT_ww_{FT_WW_TAIL}'
        newfilepath = rfunc.uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')
    #
    oldpath = f'{Output_Directory}/Dipole/DipoleFT_ALL'
    if path.exists(oldpath):
        #
        newfilepath = f'{Output_Directory}/Dipole/DipoleFT_ALL_{FT_ALL_TAIL}'
        newfilepath = rfunc.uniquify(newfilepath)
        system(f'mv {oldpath} {newfilepath}')

