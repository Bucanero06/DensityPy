#! /usr/bin/env python3.6
# Made by Ruben
import argparse
import def_functions as rfunc
# import configurationfileparser as rparse
from configurationfileparser import ParseConfigurationFile
# import errorcheck as rcheck
from errorcheck import CheckCompatibilityOfArguements
import molcasscripts as rmolcs
import chargemigratonscripts as rcms
from selectionofactivespace import SelectionOfActiveSpace
from configparser import ConfigParser
from os import path, system, remove, rmdir
import sys
import time
import multiprocessing as mp


# Main
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
    writeCM = args.writeCM
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
        rmolcs.CreateHelpInputFile()
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
    Project_Name = parser.get('Project settings', 'Project_Name')
    XYZ_Geometry = parser.get('Project settings', 'XYZ Molecule Geometry')
    Number_of_States = parser.getint('Project settings', 'Number of States')
    List_of_Orbitals = list(map(int, parser.get('Project settings', 'List_of_Active_Orbitals').split()))
    Molcas_Directory = parser.get('Project settings', 'Molcas Output Directory')

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
    Output_Directory = parser.get('Charge Migration settings', 'Output Directory')
    Config_Field_File = parser.get('Charge Migration settings', 'Field File')
    Number_of_Times = parser.getfloat('Charge Migration settings', 'Number of Times')
    Min_Time = parser.getfloat('Charge Migration settings', 'Min Time')
    Max_Time = parser.getfloat('Charge Migration settings', 'Max Time')
    Bath_Temp = parser.getfloat('Charge Migration settings', 'Bath Temperature')
    Dephasing_Factor = parser.getfloat('Charge Migration settings', 'Dephasing Factor')
    Relaxing_Factor = parser.getfloat('Charge Migration settings', 'Relaxation Factor')

    # Pump Settings
    Type_of_Pulse_Pump = parser.get('Pulses settings (Pump)', 'Type of Pulse')
    Start_Time = parser.getfloat('Pulses settings (Pump)', 'Start_Time')
    Pump_Central_Frequency = parser.getfloat('Pulses settings (Pump)', 'Pump Central Frequency')
    Pump_Periods = parser.getfloat('Pulses settings (Pump)', 'Pump Periods')
    Pump_Phase = parser.getfloat('Pulses settings (Pump)', 'Pump Phase')
    Pump_Intensity = parser.getfloat('Pulses settings (Pump)', 'Pump Intensity')
    Pump_Polarization = parser.get('Pulses settings (Pump)', 'Pump Polarization')
    # Probe Settings
    Type_of_Pulse_Probe = parser.get('Pulses settings (Probe)', 'Type of Pulse')
    Time_Delay_Start = parser.getint('Pulses settings (Probe)', 'Time Delay Start')
    Time_Delay_Stop = parser.getint('Pulses settings (Probe)', 'Time Delay Stop')
    Number_Of_PP = parser.getint('Pulses settings (Probe)', 'Number Of PP')
    Probe_Central_Frequency = parser.getfloat('Pulses settings (Probe)', 'Probe Central Frequency')
    Probe_Periods = parser.getfloat('Pulses settings (Probe)', 'Probe Periods')
    Probe_Phase = parser.getfloat('Pulses settings (Probe)', 'Probe Phase')
    Probe_Intensity = parser.getfloat('Pulses settings (Probe)', 'Probe Intensity')
    Probe_Polarization = parser.get('Pulses settings (Probe)', 'Probe Polarization')

    if writeCM:
        Time_Delay_Range = writeCM
    else:
        Change_In_Delay = (Time_Delay_Stop - Time_Delay_Start) / (Number_Of_PP - 1)
        Time_Delay_Range = list(rfunc.Float_Range(Time_Delay_Start,
                                                  Time_Delay_Stop + Change_In_Delay, Change_In_Delay))

    # Charge Migration FT Parameters
    Number_of_Omegas = parser.getfloat('Charge Migration FT settings', 'Number of Omegas')
    Min_Omegas = parser.getfloat('Charge Migration FT settings', 'Min Omega')
    Max_Omegas = parser.getfloat('Charge Migration FT settings', 'Max Omega')
    Number_of_TauOmegas = parser.getfloat('Charge Migration FT settings', 'Number of TauOmegas')
    Min_TauOmega = parser.getfloat('Charge Migration FT settings', 'Min TauOmega')
    Max_TauOmega = parser.getfloat('Charge Migration FT settings', 'Max TauOmega')
    TimeStep_FT = parser.getfloat('Charge Migration FT settings', 'TimeStep (FT)')
    WidthStep_FT = parser.getfloat('Charge Migration FT settings', 'WidthStep (FT)')
    Volume = step_size * step_size * step_size
    # print("Unit of Volume = " + str(Volume))

    # Sanity Check
    pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")
    if path.exists(pathtomolcasrc):
        molcas_workdir = rfunc.GetValueOfAsString(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
    else:
        print("Molcas configuration file ~/.Molcas/molcasrc not found. Use pymolcas --setup to create one")
        sys.exit()

    if path.exists("Sim_Time_Log"):
        remove("Sim_Time_Log")
    if limitedgrid :  # If this flag is on then it must means smartgrid is automatically True
        smartgrid = True
        gridflag= True
    # >OpenMolcas
    if pymolcas_input:
        N_points = 0
        # Prepare Input
        rfunc.MakeDirectory(Molcas_Directory)
        true_values = rmolcs.CopyInputFileToEdit(pymolcas_input,
                                                 Project_Name)  # Copies input file and returns keywords in input file

        # rfunc.MakeDirectoryNoDelete(Molcas_Directory)
        grid_time = time.time()
        if smartgrid and gridflag:
            N_points = rmolcs.Make_Better_Grid(Molcas_Directory, XYZ_Geometry, step_size, Boundary, limitedgrid)
        elif gridfile:
            system(f"cp {gridfile} {Molcas_Directory}/gridcoord")
            N_points = rfunc.File_Lenth(f'{Molcas_Directory}/gridcoord')
            print(f"Number of Points = {N_points}")
        elif gridflag:
            rmolcs.Make_Grid_Coordinates(Molcas_Directory, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax)
            N_points = Nx * Ny * Nz
        grid_time = (time.time() - grid_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write(f"Time to Make Grid {grid_time}                          time = {time.time()}")
        if gridflag or gridfile:
            rmolcs.AddGridItToManualInputFile(pymolcas_input, Project_Name, Molcas_Directory, List_of_Orbitals,
                                              N_points)
        # >Run Open Molcas
        pymol_time = time.time()
        rmolcs.Call_OpenMolcas(Project_Name)
        pymol_time = (time.time() - pymol_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write("\nTime to Run Pymolcas" + str(pymol_time) +
                       "                          time = " + str(time.time()))
        # Extract Data Required for Charge Migration
        extract_time = time.time()
        rfunc.SaveFile(pymolcas_input, Molcas_Directory + "/")
        logfilepath = rfunc.Find(Project_Name + ".log", ".", molcas_workdir)
        rfunc.SaveFile(logfilepath + "/" + Project_Name + ".log ", Molcas_Directory + "/")
        rfunc.SaveFile(XYZ_Geometry, Molcas_Directory + "/")

        if "RASSI" in true_values:
            rmolcs.GetDipolesFromLogFile(Project_Name, Molcas_Directory, Number_of_States, molcas_workdir)
        # rmolcs.Make_MU_HeatMap(Molcas_Directory)

        if gridflag:
            rmolcs.ExtractGridDensity(List_of_Orbitals, Project_Name, molcas_workdir, Molcas_Directory)

        if "RASSCF" in true_values:
            rmolcs.LoadFromh5File(Project_Name, Molcas_Directory, molcas_workdir, true_values, "DENSITY_MATRIX",
                                  "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y",
                                  "AO_MLTPL_Z", "MO_ENERGIES", "MO_VECTORS", justh5)
        extract_time = (time.time() - extract_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write("\nTime to Extract Information from Pymolcas" + str(extract_time) +
                       "                          time = " + str(time.time()))

    # Useful Flags for e.g. debugging
    if justh5:
        true_values = []
        rmolcs.LoadFromh5File(Project_Name, Molcas_Directory, molcas_workdir, true_values, "DENSITY_MATRIX",
                              "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y", "AO_MLTPL_Z",
                              "MO_ENERGIES", "MO_VECTORS", justh5)
    if justgetdensity:
        rmolcs.ExtractGridDensity(List_of_Orbitals, Project_Name, molcas_workdir, Molcas_Directory)

    if justgetdipoles:
        rmolcs.GetDipolesFromLogFile(Project_Name, Molcas_Directory, Number_of_States, molcas_workdir)
        rmolcs.Make_MU_HeatMap(Molcas_Directory)

    if cleandirectory:  ######################################################
        remove(Project_Name + '.status')
        remove(Project_Name + '.err')
        remove(Project_Name + '.log')
        remove(Project_Name + '.input')
        rmdir(f'{Output_Directory}/Pulses')
        rmdir(molcas_workdir)

    # >ChargeMigration
    if fieldfilehelp:
        rcms.Write_FieldHelp()
    if run_ChargeMigration:
        cm_time = time.time()
        # if writeCM:
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
            Field_File = givenfieldfile
        else:
            Field_File = "Field_pulses"
            # Make Field File
            rcms.Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
                              Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
                              Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
                              Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization, writeCM)
        # Run Charge Migration
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        iExcitation = -1  # Just Ground and Prepare
        iEPSILON = -1  # Just Ground and Prepare
        rcms.Call_Charge_Migration(Molcas_Directory, Output_Directory, Number_of_Times,
                                   Min_Time, Max_Time, Field_File, TimeStep_FT, WidthStep_FT, XYZ_Geometry,
                                   List_of_Orbitals, writeCM, Volume, debug_mode, weights_file, Dephasing_Factor,
                                   Relaxing_Factor, Bath_Temp, iExcitation, iEPSILON)

        if old_main:
            pass
        #
        elif parallel:
            #
            INDEX_MAT = []
            #
            for iExcitation in range(1, Number_of_States + 1):  # begin iteration
                for iEPSILON in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((iExcitation, iEPSILON))
                    #
            Number_CPU = min(((Number_of_States) * 2), mp.cpu_count())
            print(f'\nRunning in Parallel')
            print(f'Using {Number_CPU} CPU\'s out of {mp.cpu_count()} available\n')
            pool = mp.Pool(Number_CPU)
            pool.starmap_async(rcms.Call_Charge_Migration,
                               [(Molcas_Directory, Output_Directory, Number_of_Times,
                                 Min_Time, Max_Time, Field_File, TimeStep_FT, WidthStep_FT, XYZ_Geometry,
                                 List_of_Orbitals, writeCM, Volume, debug_mode, weights_file, Dephasing_Factor,
                                 Relaxing_Factor, Bath_Temp, iExcitation, iEPSILON) for iExcitation, iEPSILON in
                                INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            rcms.Call_Charge_Migration(Molcas_Directory, Output_Directory, Number_of_Times,
                                       Min_Time, Max_Time, Field_File, TimeStep_FT, WidthStep_FT, XYZ_Geometry,
                                       List_of_Orbitals, writeCM, Volume, debug_mode, weights_file, Dephasing_Factor,
                                       Relaxing_Factor, Bath_Temp, iExcitation, iEPSILON)
        cm_time = (time.time() - cm_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write(f"\nTime to Run ChargeMigration {str(cm_time / 60)} minutes")
        rfunc.SaveFile(ini_file, f"{Output_Directory}/")
    # >ChargeMigrationFT
    if run_ChargeMigrationFT:
        cmft_time = time.time()
        ##
        if save_previous:
            rcms.Save_Previous_FT(Output_Directory, Dephasing_Factor, Relaxing_Factor, Time_Delay_Start,
                                  Time_Delay_Stop,
                                  Min_Omegas, Max_Omegas, Pump_Periods, Probe_Periods, Pump_Intensity,
                                  Probe_Intensity, Number_Of_PP, Pump_Phase, Probe_Phase, TimeStep_FT, WidthStep_FT,
                                  Pump_Polarization, Probe_Polarization, Number_of_Omegas, Min_TauOmega,
                                  Max_TauOmega)
        if givenfieldfile:
            Field_File = givenfieldfile
        else:
            Field_File = Config_Field_File
            # Make Field File
            rcms.Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
                              Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
                              Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
                              Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization, writeCM)
        # Run Charge Migration FT
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if parallel:
            #
            INDEX_MAT = []
            #
            for iExcitation in range(1, Number_of_States + 1):  # begin iteration
                for iEPSILON in range(1, 3):  # begin iteration
                    #
                    INDEX_MAT.append((iExcitation, iEPSILON))
                    #
            # system(f'rm OutDir/Dipole/Dipole_State_iExcitation')
            Number_CPU = min(((Number_of_States) * 2), (mp.cpu_count()) / 2)
            print(f'\nRunning in Parallel')
            print(f'Using {Number_CPU} CPU\'s out of {mp.cpu_count()} available\n')
            pool = mp.Pool(Number_CPU)
            pool.starmap_async(rcms.Call_Charge_MigrationFT,
                               [(Molcas_Directory, Output_Directory, XYZ_Geometry, Number_of_Omegas,
                                 Min_Omegas, Max_Omegas, Number_of_TauOmegas, Min_TauOmega, Max_TauOmega,
                                 TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON) for
                                iExcitation, iEPSILON in INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            rcms.Call_Charge_MigrationFT(Molcas_Directory, Output_Directory, XYZ_Geometry, Number_of_Omegas,
                                         Min_Omegas, Max_Omegas, Number_of_TauOmegas, Min_TauOmega, Max_TauOmega,
                                         TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON)
        cmft_time = (time.time() - cmft_time)
        with open("Sim_Time_Log", 'a') as fout:
            fout.write(f"\nTime to Run ChargeMigrationFT {str(cmft_time / 60)} minutes")
    # >SpectrumReconstruction
    if run_SpectrumReconstruction:
        if save_previous:
            rcms.Save_Spectrum_Difference(Output_Directory, 'difference_' + Output_Directory, Dephasing_Factor,
                                          Relaxing_Factor, Time_Delay_Start, Time_Delay_Stop, Min_Omegas, Max_Omegas,
                                          Pump_Periods, Probe_Periods, Pump_Intensity, Probe_Intensity, Number_Of_PP,
                                          Pump_Phase, Probe_Phase, TimeStep_FT, WidthStep_FT, Pump_Polarization,
                                          Probe_Polarization, Number_of_Omegas, Min_TauOmega, Max_TauOmega)
        #
        rcms.Call_Spectrum_Reconstruction_n_Difference(Molcas_Directory, Output_Directory, XYZ_Geometry,
                                                       Number_of_Omegas, Min_Omegas, Max_Omegas, Number_of_TauOmegas,
                                                       Min_TauOmega, Max_TauOmega, debug_mode)


def main():
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
                        dest="writeCM", nargs="+", type=float, action='store')
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
    #                     dest="writeCM", action='store_true')

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
###
