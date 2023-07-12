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
from densitypy.project_utils.configurationfileparser import ParseConfigurationFile
from densitypy.project_utils.def_functions import load_json_file, Float_Range, MakeDirectory, File_Lenth, SaveFile, \
    Find, CheckCompatibilityOfArguements, GetValueOfAsString
from densitypy.project_utils.logger import logger


def get_molcas_molcasrc_variables():
    """Get the environment variables for Molcas"""
    pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")

    if path.exists(pathtomolcasrc):
        molcas_workdir = GetValueOfAsString(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
        print("Molcas Work Directory = " + molcas_workdir)
    else:
        print("Molcas configuration file ~/.Molcas/molcasrc not found. Use pymolcas --setup to create one")
        exit()
    exit()


def run(json_config, pymolcas_input, study_directory,
        run_Molcas=False,
        run_ChargeMigration=True, run_ChargeMigrationFT=False,
        run_SpectrumReconstruction=False, lus=False, debug_mode=False,
        pymolcas_input_help=False, justh5=False, justgetdipoles=False,
        justgetdensity=False, writeCM=None, givenfieldfile=None,
        weights_file=None, fieldfilehelp=False,
        gridflag: bool or str = True,

        old_main=False, parallel=False,
        save_previous=False, run_Analysis=False, run_Plot=False,
        run_Kinetics=False):
    start_time = time.time()

    # Change directory to study_directory
    orig_directory = os.getcwd()  # fixme needs to be used lol
    os.chdir(study_directory)

    if pymolcas_input_help:
        CreateHelpInputFile()
    ParseConfigurationFile(pymolcas_input, lus, run_ChargeMigration, run_ChargeMigrationFT, pymolcas_input_help, justh5,
                           justgetdensity, justgetdipoles, fieldfilehelp, run_SpectrumReconstruction,
                           json_config)
    if lus:
        # Selection of Active Space using lus argument. Requires Luscus
        pymolcas_input = CheckCompatibilityOfArguements("pymolcas_input", "lus")
        SelectionOfActiveSpace(json_config)

    # Load Values
    config = load_json_file(json_config)

    # Name of Project
    project_settings = config['project_settings']
    Project_Name = project_settings['project_name']
    XYZ_Geometry = project_settings['xyz_molecule_geometry']
    Number_of_States = project_settings['number_of_states']
    List_of_Orbitals = project_settings['list_of_active_orbitals']
    Molcas_Directory = project_settings['molcas_output_directory']

    # Grid Settings
    grid_settings = config['grid_settings']
    Nx = grid_settings['number_of_points']['x_axis']
    Ny = grid_settings['number_of_points']['y_axis']
    Nz = grid_settings['number_of_points']['z_axis']
    xmin = grid_settings['x_range']['min']
    xmax = grid_settings['x_range']['max']
    ymin = grid_settings['y_range']['min']
    ymax = grid_settings['y_range']['max']
    zmin = grid_settings['z_range']['min']
    zmax = grid_settings['z_range']['max']
    step_size = grid_settings['step_size']
    Boundary = grid_settings['boundary']

    # Charge Migration Parameters
    charge_migration_settings = config['charge_migration_settings']
    Output_Directory = charge_migration_settings['output_directory']
    Config_Field_File = charge_migration_settings['field_file']
    Number_of_Times = charge_migration_settings['number_of_times']
    Min_Time = charge_migration_settings['time_range']['min']
    Max_Time = charge_migration_settings['time_range']['max']
    Bath_Temp = charge_migration_settings['bath_temperature']
    Dephasing_Factor = charge_migration_settings['dephasing_factor']
    Relaxing_Factor = charge_migration_settings['relaxation_factor']

    # Pump Settings
    pump_settings = config['pulses_settings']['pump']
    Type_of_Pulse_Pump = pump_settings['type_of_pulse']
    Start_Time = pump_settings['start_time']
    Pump_Central_Frequency = pump_settings['central_frequency']
    Pump_Periods = pump_settings['periods']
    Pump_Phase = pump_settings['phase']
    Pump_Intensity = pump_settings['intensity']
    Pump_Polarization = pump_settings['polarization']

    # Probe Settings
    probe_settings = config['pulses_settings']['probe']
    Type_of_Pulse_Probe = probe_settings['type_of_pulse']
    Time_Delay_Start = probe_settings['time_delay_range']['start']
    Time_Delay_Stop = probe_settings['time_delay_range']['stop']
    Number_Of_PP = probe_settings['number_of_pp']
    Probe_Central_Frequency = probe_settings['central_frequency']
    Probe_Periods = probe_settings['periods']
    Probe_Phase = probe_settings['phase']
    Probe_Intensity = probe_settings['intensity']
    Probe_Polarization = probe_settings['polarization']

    # Not sure how to handle the 'writeCM' variable, assuming it is a flag in JSON config
    writeCM = config.get('writeCM', False)  # fixme
    if writeCM:
        Time_Delay_Range = writeCM
    else:
        Change_In_Delay = (Time_Delay_Stop - Time_Delay_Start) / (Number_Of_PP - 1)
        Time_Delay_Range = list(Float_Range(Time_Delay_Start,
                                            Time_Delay_Stop + Change_In_Delay, Change_In_Delay))

    # Charge Migration FT Parameters
    charge_migration_ft_settings = config['charge_migration_ft_settings']
    Number_of_Omegas = charge_migration_ft_settings['number_of_omegas']
    Min_Omegas = charge_migration_ft_settings['omega_range']['min']
    Max_Omegas = charge_migration_ft_settings['omega_range']['max']
    Number_of_TauOmegas = charge_migration_ft_settings['number_of_tauomegas']
    Min_TauOmega = charge_migration_ft_settings['tauomega_range']['min']
    Max_TauOmega = charge_migration_ft_settings['tauomega_range']['max']
    TimeStep_FT = charge_migration_ft_settings['timestep_ft']
    WidthStep_FT = charge_migration_ft_settings['widthstep_ft']
    Volume = step_size * step_size * step_size

    logger.info("Unit of Volume = " + str(Volume))

    # >OpenMolcas
    if run_Molcas:
        N_points = 0
        # Prepare Input
        MakeDirectory(Molcas_Directory)
        # Copies input file and returns keywords in input file
        true_values = CopyInputFileToEdit(pymolcas_input, Project_Name, Molcas_Directory)

        # MakeDirectoryNoDelete(Molcas_Directory)

        # gridflag can be
        #   True = use gridit on molcas with with uniform grid points
        #   False = Do not use or create any grid or add gridit to input file
        #   'smart' =  Use our smart grid algorithm to create a grid based on the molecule and step size # todo
        #   'limited' = Use our truncated grid algorithm to create a grid based on the molecule and config settings
        #   'manual' = Use a grid file to create a grid, this is the only option that requires a grid file to be specified in the config file
        if gridflag in ['smart', 'limited']:
            limitedgrid = True if gridflag == 'limited' else False
            N_points = Make_Better_Grid(Molcas_Directory, XYZ_Geometry, step_size, Boundary, limitedgrid)
        elif gridflag == 'manual':
            gridfile = config['grid_settings']['default_grid_file']
            system(f"cp {gridfile} {Molcas_Directory}/gridcoord")
            N_points = File_Lenth(f'{Molcas_Directory}/gridcoord')
            logger.info(f"Number of Points = {N_points}")
        elif gridflag:
            Make_Grid_Coordinates(Molcas_Directory, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax)
            N_points = Nx * Ny * Nz
        else:
            assert gridflag is False, f"gridflag is {gridflag} but should be False, True, 'smart', 'limited', or 'manual'"

        if gridflag:
            AddGridItToManualInputFile(pymolcas_input, Project_Name, Molcas_Directory, List_of_Orbitals, N_points)
        # >Run Open Molcas
        Call_OpenMolcas(Project_Name, Molcas_Directory)
        # Extract Data Required for Charge Migration
        SaveFile(pymolcas_input, Molcas_Directory + "/")
        logfilepath = Find(Project_Name + ".log", ".", Molcas_Directory)
        SaveFile(logfilepath + "/" + Project_Name + ".log ", Molcas_Directory + "/")
        SaveFile(XYZ_Geometry, Molcas_Directory + "/")

        if "RASSI" in true_values:
            GetDipolesFromLogFile(Project_Name, Molcas_Directory, Number_of_States, Molcas_Directory)
            Make_MU_HeatMap(Molcas_Directory)
            project = DipolesProject(log_path=f'{Molcas_Directory}/{Project_Name}.log')
            project.write_csvs(directory=Molcas_Directory)
            fig = project.get_mu_heatmaps()
            fig.savefig(f'{Molcas_Directory}/dipole-heatmap.png')

        if gridflag:
            ExtractGridDensity(List_of_Orbitals, Project_Name, Molcas_Directory, Molcas_Directory)

        if "RASSCF" in true_values:
            LoadFromh5File(Project_Name, Molcas_Directory, Molcas_Directory, true_values, "DENSITY_MATRIX",
                           "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y",
                           "AO_MLTPL_Z", "MO_ENERGIES", "MO_VECTORS", justh5)

    # Useful Flags for e.g. debugging
    if justh5:
        true_values = []
        LoadFromh5File(Project_Name, Molcas_Directory, Molcas_Directory, true_values, "DENSITY_MATRIX",
                       "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y", "AO_MLTPL_Z",
                       "MO_ENERGIES", "MO_VECTORS", justh5)
    if justgetdensity:
        ExtractGridDensity(List_of_Orbitals, Project_Name, Molcas_Directory, Molcas_Directory)

    if justgetdipoles:
        GetDipolesFromLogFile(Project_Name, Molcas_Directory, Number_of_States, Molcas_Directory)
        Make_MU_HeatMap(Molcas_Directory)

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
            Field_File = "Field_pulses"
            # Make Field File
            Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
                         Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
                         Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
                         Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization, writeCM)

        # Run Charge Migration
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        iExcitation = -1  # Just Ground and Prepare
        iEPSILON = -1  # Just Ground and Prepare
        Call_Charge_Migration(Molcas_Directory, Output_Directory, Number_of_Times,
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
            Number_CPU = min(((Number_of_States) * 2), os.cpu_count())
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {os.cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_Migration,
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
            Call_Charge_Migration(Molcas_Directory, Output_Directory, Number_of_Times,
                                  Min_Time, Max_Time, Field_File, TimeStep_FT, WidthStep_FT, XYZ_Geometry,
                                  List_of_Orbitals, writeCM, Volume, debug_mode, weights_file, Dephasing_Factor,
                                  Relaxing_Factor, Bath_Temp, iExcitation, iEPSILON)
        SaveFile(json_config, f"{Output_Directory}/")
    # >ChargeMigrationFT
    if run_ChargeMigrationFT:
        cmft_time = time.time()
        ##
        if save_previous:
            Save_Previous_FT(Output_Directory, Dephasing_Factor, Relaxing_Factor, Time_Delay_Start,
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
            Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
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
            Number_CPU = min(((Number_of_States) * 2), (cpu_count()) / 2)
            logger.info(f'\nRunning in Parallel')
            logger.info(f'Using {Number_CPU} CPU\'s out of {cpu_count()} available\n')
            pool = Pool(Number_CPU)
            pool.starmap_async(Call_Charge_MigrationFT,
                               [(Molcas_Directory, Output_Directory, XYZ_Geometry, Number_of_Omegas,
                                 Min_Omegas, Max_Omegas, Number_of_TauOmegas, Min_TauOmega, Max_TauOmega,
                                 TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON) for
                                iExcitation, iEPSILON in INDEX_MAT]).get()
            pool.close()
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            iExcitation = 0
            iEPSILON = 0
            Call_Charge_MigrationFT(Molcas_Directory, Output_Directory, XYZ_Geometry, Number_of_Omegas,
                                    Min_Omegas, Max_Omegas, Number_of_TauOmegas, Min_TauOmega, Max_TauOmega,
                                    TimeStep_FT, WidthStep_FT, Field_File, debug_mode, iExcitation, iEPSILON)
    # >SpectrumReconstruction
    if run_SpectrumReconstruction:
        if save_previous:
            Save_Spectrum_Difference(Output_Directory, 'difference_' + Output_Directory, Dephasing_Factor,
                                     Relaxing_Factor, Time_Delay_Start, Time_Delay_Stop, Min_Omegas, Max_Omegas,
                                     Pump_Periods, Probe_Periods, Pump_Intensity, Probe_Intensity, Number_Of_PP,
                                     Pump_Phase, Probe_Phase, TimeStep_FT, WidthStep_FT, Pump_Polarization,
                                     Probe_Polarization, Number_of_Omegas, Min_TauOmega, Max_TauOmega)
        #
        Call_Spectrum_Reconstruction_n_Difference(Molcas_Directory, Output_Directory, XYZ_Geometry,
                                                  Number_of_Omegas, Min_Omegas, Max_Omegas, Number_of_TauOmegas,
                                                  Min_TauOmega, Max_TauOmega, debug_mode)

    # Change back to original directory
    os.chdir(orig_directory)


if __name__ == "__main__":
    run(json_config="chargemigration.json",
        study_directory="/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/",

        # Molcas Group
        run_Molcas=False,
        pymolcas_input='inputhelp.input',
        pymolcas_input_help=False,
        gridflag=True,
        lus=False,

        # Charge Migration Group
        run_ChargeMigration=True,
        run_ChargeMigrationFT=False,
        writeCM=None,
        weights_file=None,
        parallel=False,
        old_main=True,
        givenfieldfile=None,
        fieldfilehelp=False,
        save_previous=False,

        # Analysis Group
        run_SpectrumReconstruction=False,
        justh5=False, justgetdipoles=False, justgetdensity=False,
        run_Plot=False,
        run_Kinetics=False,

        # Development Group
        debug_mode=False,

        )
