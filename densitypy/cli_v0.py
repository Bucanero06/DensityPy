#! /usr/bin/env python3.6
# Made by Ruben
import argparse
import os
import sys
import time
from configparser import ConfigParser
from multiprocessing import Pool, cpu_count
from os import remove, rmdir, path, system

from chargemigratonscripts import Write_FieldHelp, Write_Pulses, Call_Charge_Migration, Save_Previous_FT, \
    Call_Charge_MigrationFT, Save_Spectrum_Difference, Call_Spectrum_Reconstruction_n_Difference
from configurationfileparser import ParseConfigurationFile
from def_functions import Float_Range, GetValueOfAsString, MakeDirectory, File_Lenth, SaveFile, Find
from errorcheck import CheckCompatibilityOfArguements
from molcasscripts import CreateHelpInputFile, CopyInputFileToEdit, Make_Better_Grid, Make_Grid_Coordinates, \
    AddGridItToManualInputFile, Call_OpenMolcas, GetDipolesFromLogFile, ExtractGridDensity, LoadFromh5File, \
    Make_MU_HeatMap
from selectionofactivespace import SelectionOfActiveSpace


def run(args):
    # Gets Run_Time_Parameters
    pymolcas_input = args.input
    run_ChargeMigration = args.chargemigrationtrue
    run_ChargeMigrationFT = args.chargemigrationFTtrue
    run_SpectrumReconstruction = args.SpectrumTrue
    lus = args.lustrue
    debug_mode = args.debug
    input_help = args.input_help
    justh5 = args.justh5
    justgetdipoles = args.justgetdipoles
    justgetdensity = args.justgetdensity
    write_charge_migration = args.write_charge_migration
    givenfieldfile = args.givenfieldfile
    weights_file = args.weights_file
    fieldfilehelp = args.fieldfilehelp
    limitedgrid = args.limitedgrid
    smartgrid = args.smartgrid
    gridflag = args.gridflag
    gridfile = args.gridfile
    ini_file = args.ini_file
    cleandirectory = args.cleandirectory
    just_run_input = args.just_run_input
    old_main = args.old_main
    parallel = args.parallel
    save_previous = args.save_previous
    if input_help:
        CreateHelpInputFile()
    ParseConfigurationFile(pymolcas_input, lus, run_ChargeMigration, run_ChargeMigrationFT, input_help, justh5,
                           justgetdensity, justgetdipoles, fieldfilehelp, run_SpectrumReconstruction,
                           ini_file)
    if lus:
        # Selection of Active Space using lus argument. Requires Luscus
        pymolcas_input = CheckCompatibilityOfArguements("pymolcas_input", "lus")
        SelectionOfActiveSpace(ini_file)

    # >Load Values
    parser = ConfigParser()
    if ini_file:
        parser.read(ini_file)
    else:
        parser.read('chargemigration.ini')

    # Name of Project
    project_name = parser.get('Project settings', 'project_name')
    XYZ_Geometry = parser.get('Project settings', 'XYZ Molecule Geometry')
    Number_of_States = parser.getint('Project settings', 'Number of States')
    List_of_Orbitals = list(map(int, parser.get('Project settings', 'List_of_Active_Orbitals').split()))
    molcas_directory = parser.get('Project settings', 'Molcas Output Directory')

    # Grid Settings
    Nx = parser.getint('GRID settings', 'Number of Points, X axis')
    Ny = parser.getint('GRID settings', 'Number of Points, Y axis')
    Nz = parser.getint('GRID settings', 'Number of Points, Z axis')
    xmin = parser.getfloat('GRID settings', 'X MIN')
    xmax = parser.getfloat('GRID settings', 'X MAX')
    ymin = parser.getfloat('GRID settings', 'Y MIN')
    ymax = parser.getfloat('GRID settings', 'Y MAX')
    zmin = parser.getfloat('GRID settings', 'Z MIN')
    zmax = parser.getfloat('GRID settings', 'Z MAX')
    step_size = parser.getfloat('GRID settings', 'Step Size')
    Boundary = parser.getfloat('GRID settings', 'Boundary')

    # Charge Migration Parameters
    output_directory = parser.get('Charge Migration settings', 'Output Directory')
    Config_field_file = parser.get('Charge Migration settings', 'Field File')
    number_of_times = parser.getfloat('Charge Migration settings', 'Number of Times')
    min_time = parser.getfloat('Charge Migration settings', 'Min Time')
    max_time = parser.getfloat('Charge Migration settings', 'Max Time')
    Bath_Temp = parser.getfloat('Charge Migration settings', 'Bath Temperature')
    dephasing_factor = parser.getfloat('Charge Migration settings', 'Dephasing Factor')
    relaxation_factor = parser.getfloat('Charge Migration settings', 'Relaxation Factor')

    # Pump Settings
    type_of_pulse_pump = parser.get('Pulses settings (Pump)', 'Type of Pulse')
    start_time = parser.getfloat('Pulses settings (Pump)', 'start_time')
    pump_central_frequency = parser.getfloat('Pulses settings (Pump)', 'Pump Central Frequency')
    pump_periods = parser.getfloat('Pulses settings (Pump)', 'Pump Periods')
    pump_phase = parser.getfloat('Pulses settings (Pump)', 'Pump Phase')
    pump_intensity = parser.getfloat('Pulses settings (Pump)', 'Pump Intensity')
    pump_polarization = parser.get('Pulses settings (Pump)', 'Pump Polarization')
    # Probe Settings
    type_of_pulse_probe = parser.get('Pulses settings (Probe)', 'Type of Pulse')
    time_delay_start = parser.getint('Pulses settings (Probe)', 'Time Delay Start')
    time_delay_stop = parser.getint('Pulses settings (Probe)', 'Time Delay Stop')
    Number_Of_PP = parser.getint('Pulses settings (Probe)', 'Number Of PP')
    probe_central_frequency = parser.getfloat('Pulses settings (Probe)', 'Probe Central Frequency')
    probe_periods = parser.getfloat('Pulses settings (Probe)', 'Probe Periods')
    probe_phase = parser.getfloat('Pulses settings (Probe)', 'Probe Phase')
    probe_intensity = parser.getfloat('Pulses settings (Probe)', 'Probe Intensity')
    probe_polarization = parser.get('Pulses settings (Probe)', 'Probe Polarization')

    if write_charge_migration:
        time_delay_range = write_charge_migration
    else:
        Change_In_Delay = (time_delay_stop - time_delay_start) / (Number_Of_PP - 1)
        time_delay_range = list(Float_Range(time_delay_start,
                                                  time_delay_stop + Change_In_Delay, Change_In_Delay))

    # Charge Migration FT Parameters
    Number_of_Omegas = parser.getfloat('Charge Migration FT settings', 'Number of Omegas')
    min_omegas = parser.getfloat('Charge Migration FT settings', 'Min Omega')
    max_omegas = parser.getfloat('Charge Migration FT settings', 'Max Omega')
    number_of_tau_omegas = parser.getfloat('Charge Migration FT settings', 'Number of TauOmegas')
    min_tau_omega = parser.getfloat('Charge Migration FT settings', 'Min TauOmega')
    max_tau_omega = parser.getfloat('Charge Migration FT settings', 'Max TauOmega')
    ft_time_step = parser.getfloat('Charge Migration FT settings', 'TimeStep (FT)')
    ft_width_step = parser.getfloat('Charge Migration FT settings', 'WidthStep (FT)')
    Volume = step_size * step_size * step_size
    # print("Unit of Volume = " + str(Volume))

    # # Sanity Check
    # pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")
    # if path.exists(pathtomolcasrc):
    #     molcas_workdir = GetValueOfAsString(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
    # else:
    #     print("Molcas configuration file ~/.Molcas/molcasrc not found. Use pymolcas --setup to create one")
    #     sys.exit()
    molcas_workdir = "MOLCAS_WORKDIR"

    if path.exists("Sim_Time_Log"):
        remove("Sim_Time_Log")
    if limitedgrid :  # If this flag is on then it must means smartgrid is automatically True
        smartgrid = True
        gridflag= True
    # >OpenMolcas
    if pymolcas_input:
        n_points = 0
        # Prepare Input
        MakeDirectory(molcas_directory)
        true_values = CopyInputFileToEdit(pymolcas_input,
                                                 project_name)  # Copies input file and returns keywords in input file

        # MakeDirectoryNoDelete(molcas_directory)
        grid_time = time.time()
        if smartgrid and gridflag:
            n_points = Make_Better_Grid(molcas_directory, XYZ_Geometry, step_size, Boundary, limitedgrid)
        elif gridfile:
            system(f"cp {gridfile} {molcas_directory}/gridcoord")
            n_points = File_Lenth(f'{molcas_directory}/gridcoord')
            print(f"Number of Points = {n_points}")
        elif gridflag:
            Make_Grid_Coordinates(molcas_directory, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax)
            n_points = Nx * Ny * Nz
        grid_time = (time.time() - grid_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write(f"Time to Make Grid {grid_time}                          time = {time.time()}")
        if gridflag or gridfile:
            AddGridItToManualInputFile(pymolcas_input, project_name, molcas_directory, List_of_Orbitals,
                                              n_points)
        # >Run Open Molcas
        pymol_time = time.time()
        Call_OpenMolcas(project_name)
        pymol_time = (time.time() - pymol_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write("\nTime to Run Pymolcas" + str(pymol_time) +
                       "                          time = " + str(time.time()))
        # Extract Data Required for Charge Migration
        extract_time = time.time()
        SaveFile(pymolcas_input, molcas_directory + "/")
        logfilepath = Find(project_name + ".log", ".", molcas_workdir)
        SaveFile(logfilepath + "/" + project_name + ".log ", molcas_directory + "/")
        SaveFile(XYZ_Geometry, molcas_directory + "/")

        if "RASSI" in true_values:
            GetDipolesFromLogFile(project_name, molcas_directory, Number_of_States, molcas_workdir)
        # Make_MU_HeatMap(molcas_directory)

        if gridflag:
            ExtractGridDensity(List_of_Orbitals, project_name, molcas_workdir, molcas_directory)

        if "RASSCF" in true_values:
            LoadFromh5File(project_name, molcas_directory, molcas_workdir, true_values, "DENSITY_MATRIX",
                                  "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y",
                                  "AO_MLTPL_Z", "MO_ENERGIES", "MO_VECTORS", justh5)
        extract_time = (time.time() - extract_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write("\nTime to Extract Information from Pymolcas" + str(extract_time) +
                       "                          time = " + str(time.time()))

    # Useful Flags for e.g. debugging
    if justh5:
        true_values = []
        LoadFromh5File(project_name, molcas_directory, molcas_workdir, true_values, "DENSITY_MATRIX",
                              "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y", "AO_MLTPL_Z",
                              "MO_ENERGIES", "MO_VECTORS", justh5)
    if justgetdensity:
        ExtractGridDensity(List_of_Orbitals, project_name, molcas_workdir, molcas_directory)

    if justgetdipoles:
        GetDipolesFromLogFile(project_name, molcas_directory, Number_of_States, molcas_workdir)
        Make_MU_HeatMap(molcas_directory)

    if cleandirectory:  ######################################################
        remove(project_name + '.status')
        remove(project_name + '.err')
        remove(project_name + '.log')
        remove(project_name + '.input')
        rmdir(f'{output_directory}/Pulses')
        rmdir(molcas_workdir)

    # >ChargeMigration
    if fieldfilehelp:
        Write_FieldHelp()
    if run_ChargeMigration:
        cm_time = time.time()
        # if write_charge_migration:
        #     pass
        # elif (weights_file == None):
        #     if path.exists("Weights_File"):
        #         print("Found Weight_File in curent directory, did you forget to use \"-w\"")
        #         userinput_weights = input("\"Y\" to use Weights_File:").upper().replace(
        #             " ", "")
        #         yes_input = ["yes", "y", "Yes", "Y", "YES"]
        #         if (userinput_weights in yes_input):
        #             weights_file = "Weights_File"
        #             print("Using \"-w\" flag")
        if givenfieldfile:
            field_file = givenfieldfile
        else:
            field_file = "Field_pulses"
            # Make Field File
            Write_Pulses(field_file, type_of_pulse_pump, start_time, pump_central_frequency,
                              pump_periods, pump_phase, pump_intensity, pump_polarization,
                              type_of_pulse_probe, time_delay_range, probe_central_frequency,
                              probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration)
        # Run Charge Migration
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        i_excitation = -1  # Just Ground and Prepare
        i_epsilon = -1  # Just Ground and Prepare
        Call_Charge_Migration(molcas_directory, output_directory, number_of_times,
                                   min_time, max_time, field_file, ft_time_step, ft_width_step, XYZ_Geometry,
                                   List_of_Orbitals, write_charge_migration, Volume, debug_mode, weights_file, dephasing_factor,
                                   relaxation_factor, Bath_Temp, i_excitation, i_epsilon)

        if old_main:
            pass
        #
        elif parallel:
            #
            INDEX_MAT = []
            #
            for i_excitation in range(1, Number_of_States + 1):  # begin iteration
                for i_epsilon in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((i_excitation, i_epsilon))
                    #
            Number_CPU = min(((Number_of_States) * 2), cpu_count())
            print(f'\nRunning in Parallel')
            print(f'Using {Number_CPU} CPU\'s out of {cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_Migration,
                               [(molcas_directory, output_directory, number_of_times,
                                 min_time, max_time, field_file, ft_time_step, ft_width_step, XYZ_Geometry,
                                 List_of_Orbitals, write_charge_migration, Volume, debug_mode, weights_file, dephasing_factor,
                                 relaxation_factor, Bath_Temp, i_excitation, i_epsilon) for i_excitation, i_epsilon in
                                INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            i_excitation = 0
            i_epsilon = 0
            Call_Charge_Migration(molcas_directory, output_directory, number_of_times,
                                       min_time, max_time, field_file, ft_time_step, ft_width_step, XYZ_Geometry,
                                       List_of_Orbitals, write_charge_migration, Volume, debug_mode, weights_file, dephasing_factor,
                                       relaxation_factor, Bath_Temp, i_excitation, i_epsilon)
        cm_time = (time.time() - cm_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write(f"\nTime to Run ChargeMigration {str(cm_time / 60)} minutes")
        SaveFile(ini_file, f"{output_directory}/")
    # >ChargeMigrationFT
    if run_ChargeMigrationFT:
        cmft_time = time.time()
        ##
        if save_previous:
            Save_Previous_FT(output_directory, dephasing_factor, relaxation_factor, time_delay_start,
                                  time_delay_stop,
                                  min_omegas, max_omegas, pump_periods, probe_periods, pump_intensity,
                                  probe_intensity, Number_Of_PP, pump_phase, probe_phase, ft_time_step, ft_width_step,
                                  pump_polarization, probe_polarization, Number_of_Omegas, min_tau_omega,
                                  max_tau_omega)
        if givenfieldfile:
            field_file = givenfieldfile
        else:
            field_file = Config_field_file
            # Make Field File
            Write_Pulses(field_file, type_of_pulse_pump, start_time, pump_central_frequency,
                              pump_periods, pump_phase, pump_intensity, pump_polarization,
                              type_of_pulse_probe, time_delay_range, probe_central_frequency,
                              probe_periods, probe_phase, probe_intensity, probe_polarization, write_charge_migration)
        # Run Charge Migration FT
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if parallel:
            #
            INDEX_MAT = []
            #
            for i_excitation in range(1, Number_of_States + 1):  # begin iteration
                for i_epsilon in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((i_excitation, i_epsilon))
                    #
            # system(f'rm OutDir/Dipole/Dipole_State_i_excitation')
            Number_CPU = min(((Number_of_States) * 2), (cpu_count()) / 2)
            print(f'\nRunning in Parallel')
            print(f'Using {Number_CPU} CPU\'s out of {cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_MigrationFT,
                               [(molcas_directory, output_directory, XYZ_Geometry, Number_of_Omegas,
                                 min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                 ft_time_step, ft_width_step, field_file, debug_mode, i_excitation, i_epsilon) for
                                i_excitation, i_epsilon in INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            i_excitation = 0
            i_epsilon = 0
            Call_Charge_MigrationFT(molcas_directory, output_directory, XYZ_Geometry, Number_of_Omegas,
                                         min_omegas, max_omegas, number_of_tau_omegas, min_tau_omega, max_tau_omega,
                                         ft_time_step, ft_width_step, field_file, debug_mode, i_excitation, i_epsilon)
        cmft_time = (time.time() - cmft_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write(f"\nTime to Run ChargeMigrationFT {str(cmft_time / 60)} minutes")
    # >SpectrumReconstruction
    if run_SpectrumReconstruction:
        if save_previous:
            Save_Spectrum_Difference(output_directory, 'difference_' + output_directory, dephasing_factor,
                                          relaxation_factor, time_delay_start, time_delay_stop, min_omegas, max_omegas,
                                          pump_periods, probe_periods, pump_intensity, probe_intensity, Number_Of_PP,
                                          pump_phase, probe_phase, ft_time_step, ft_width_step, pump_polarization,
                                          probe_polarization, Number_of_Omegas, min_tau_omega, max_tau_omega)
        #
        Call_Spectrum_Reconstruction_n_Difference(molcas_directory, output_directory, XYZ_Geometry,
                                                       Number_of_Omegas, min_omegas, max_omegas, number_of_tau_omegas,
                                                       min_tau_omega, max_tau_omega, debug_mode)


def cli_main():
    parser = argparse.ArgumentParser(description="Help main.py script used for Molcas and ChargeMigration")

    groupc = parser.add_argument_group(description="Configuration File")
    groupc.add_argument("-ini", help="Specify configuration file, otherwise chargemigration.ini will be used as "
                                     "default or created",
                        dest="ini_file", type=str, action='store')

    # ..Command line arguements for Molcas
    group1 = parser.add_argument_group(description="Command line arguements for Molcas")
    group1.add_argument("-i", help="input file for OpenMolcas (pymolcas command). "
                                   "Used to indicate input file for Molcas. Name should differ from "
                                   "<ProjectName>.input to avoid overwriting.",
                        dest="input", type=str, action='store')
    group1.add_argument("-i?", help="Helps user by creating a starting input file for OpenMolcas (pymolcas command),"
                                    " which can be edited. To use open molcas with this file use -i "
                                    "followed by filename",
                        dest="input_help", action='store_true')
    group1.add_argument("-lus", help="If true, Uses Luscus to select active space file "
                                     "($ProjectName.GvOrb)",
                        dest="lustrue", action='store_true')
    group1.add_argument("-h5", help="This will only extract data found inside <Project>.rasscf.h5 "
                                    "from a previous run of Molcas.(Looks first inside of <Molcas_WORKDIR>)",
                        dest="justh5", action='store_true')
    group1.add_argument("-getdensities", help="This will only extract data found inside "
                                              "<Project>.grid from a previous run of Molcas."
                                              "(Looks first in the current directory)",
                        dest="justgetdensity", action='store_true')
    group1.add_argument("-getdipoles", help="This will only extract the dipoles inside of "
                                            "<Project>.log from a previous run of Molcas."
                                            "(Looks first in the current directory)",
                        dest="justgetdipoles", action='store_true')
    group1.add_argument("-gridflag", help="Required if grid is needed for post-molcas orbital density calculations",
                        dest="gridflag", action='store_true')
    group1.add_argument("-smartgrid", help="Automatically makes a grid for an arbitrary molecule using the size steps "
                                           "and Boundary(distance from the nuclei to us as cutoff for grid) input. Can be used with optional flag \"-limitedgrid\"."
                                           "Typically results in a cubic grid with unequal side lengths",
                        dest="smartgrid", action='store_true')
    group1.add_argument("-limitedgrid",
                        help="When used with flag \"-smartgrid\", new grid will only be composed of points \"Boundary\" distance away from each nuclei (reduces # points)."
                             "Typically results in oval grid",
                        dest="limitedgrid", action='store_true')
    group1.add_argument("-g", help="input grid for pymolcas. (optional)",
                        dest="gridfile", type=str, action='store')
    group1.add_argument("-clean",
                        help="Clean files created by molcas that are not used anymore in the current directory",
                        dest="cleandirectory", action='store_true')
    group1.add_argument("-r",
                        help="Just run the pymolcas input file",
                        dest="just_run_input", action='store_true')
    group1.add_argument("-old",
                        help="use old main, debugging, rm",
                        dest="old_main", action='store_true')

    # ..Command line arguements for ChargeMigration
    group2 = parser.add_argument_group(description="Command line arguements for ChargeMigration")
    group2.add_argument("-cm", help="If True, calls Charge Migration "
                                    "using the parameters in chargemigration.ini ",
                        dest="chargemigrationtrue", action='store_true')
    group2.add_argument("-field", help="Specifies field file containing the sequences of pulses to be "
                                       "used by ChargeMigration; otherwise the default is to look for \'Field_pulses\' "
                                       "or make one using infromation from chargemigration.ini. Works for both "
                                       "ChargeMigration or ChargeMigrationFT",
                        dest="givenfieldfile", type=str, action='store')
    group2.add_argument("-w", help="Feed ChargeMigration code a premade File containing Becke's Weights for the Grid, "
                                   "saves time in subsequent calculations after already computed",
                        dest="weights_file", type=str, action='store')
    group2.add_argument("-field?", help="Creates exemple template file which can be used with \'-field\'to specify the "
                                        "sequences of pulses to be used by ChargeMigration.",
                        dest="fieldfilehelp", action='store_true')
    group2.add_argument("-sden", help="Save Charge Density to file when running '-cm'(time consuming !)",
                        dest="write_charge_migration", nargs="+", type=float, action='store')
    group2.add_argument("-p", help="run in parallel",
                        dest="parallel", action='store_true')
    group2.add_argument("-debug", help="Runs ChargeMigration or ChargeMigrationFT in debug mode",
                        dest="debug", action='store_true')

    # ..Command line arguements for ChargeMigrationFT
    group3 = parser.add_argument_group(description="Command line arguements for ChargeMigrationFT")
    group3.add_argument("-ft", help="If True, calls Charge Migration FT"
                                    "using the parameters in chargemigration.ini ",
                        dest="chargemigrationFTtrue", action='store_true')
    group3.add_argument("-s", help="Calls Reconstruction of Spectrum code and scripts",
                        dest="SpectrumTrue", action='store_true')
    group3.add_argument("-save",
                        help="Saves FT and difference in spectra files from previous runs and renames them with useful information if "
                             "they are found within the current {SIM} directory in use or the working director in the case of the difference",
                        dest="save_previous", action='store_true', default=True)

    # group3.add_argument("-field", help="Specifies field file containing the sequences of pulses to be "
    #                                    "used by ChargeMigration; otherwise the default is to use the parameters "
    #                                    "found in chargemigration.ini to make and name the file",
    #                     dest="givenfieldfile", type=str, action='store')
    # group3.add_argument("-field?", help="Creates exemple template file which can be used with \'-field\'to specify the "
    #                                     "sequences of pulses to be used by ChargeMigration.",
    #                     dest="fieldfilehelp", action='store_true')
    # group3.add_argument("-sden", help="Save Charge Density to file when running '-cm'(time consuming !)",
    #                     dest="write_charge_migration", action='store_true')

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    cli_main()



###
