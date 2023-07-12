#! /usr/bin/env python3.6
# Made by Ruben
# >import configurationfileparser as rparse

import configparser
import sys
from configparser import ConfigParser
from os import path


def ParseConfigurationFile(option1, option2, option3, option4, option5, option6, option7, option8, option9, option10,
                           configuration_file):
    data = [option1, option2, option3, option4, option5, option6, option7, option8, option9, option10,
            configuration_file]
    if configuration_file:
        if path.exists(configuration_file):
            if not any(data):
                print("Found " + str(configuration_file) + " but nothing was ran: "
                                                           "\n   No command line  action arguments used.")
                sys.exit()
            else:
                print('Found chargemigration.ini')
                print(str(configuration_file))


    elif path.exists("chargemigration.ini"):
        if not any(data):
            print("Found chargemigration.ini but nothing was ran: "
                  "\n   No command line  action arguments used.")
            sys.exit()
        else:
            print('Found chargemigration.ini')
    else:
        print("No configuration file found"
              "\nMaking chargemigration.ini using defaults. "
              "\nCheck/Edit configuration file before running"
              "\nExiting")
        config = ConfigParser()

        # Project Settings
        config['Project settings'] = {
            'Project_Name': 'NMA',
            'XYZ Molecule Geometry': 'NMA.xyz',
            'Basis': 'ANO-L-MB',  # only needed for use with -sas arguement. LUSUCS GUI
            'Number of States': '5',
            'List_of_Active_Orbitals': '17 18 19 20 21',  # Select which Orbitals to work with
            'Molcas Output Directory': 'NMA_output'
            # "Output Directory for Molcas, will be used also by other programs"
        }

        # >Settings for Grid (Number of points per axis and size of Cube)
        config['GRID settings'] = {
            'Number of Points, X axis': '15',
            'Number of Points, Y axis': '15',
            'Number of Points, Z axis': '15',
            'X MIN': '-5',
            'X MAX': '6',
            'Y MIN': '-4',
            'Y MAX': '5',
            'Z MIN': '-4',
            'Z MAX': '4',
            'Step Size': '0.1',
            'Boundary': '2'}

        # >Settings for Charge Migration
        config['Charge Migration settings'] = {
            'Output Directory': 'sim',
            'Field File': 'Field_File',
            'Number of Times': '1001',
            'Min Time': '-4000',
            'Max Time': '5000',
            'Bath Temperature': '3273.75',
            'Dephasing Factor': '0.001',
            'Relaxation Factor': '0.001',
        }

        # >Pump Settings
        config['Pulses settings (Pump)'] = {
            'Type of Pulse': 'G',
            'Start_Time': '0.000',
            'Pump Central Frequency': '0.4537',
            'Pump Periods': '5',
            'Pump Phase': '0',
            'Pump Intensity': '0.00001',
            'Pump Polarization': '90 90'}

        # >Probe Settings
        config['Pulses settings (Probe)'] = {
            'Type of Pulse': 'G',
            'Time Delay Start': '-2000',
            'Time Delay Stop': '3000',
            'Number Of PP': '101',
            'Probe Central Frequency': '0.07',
            'Probe Periods': '5',
            'Probe Phase': '0',
            'Probe Intensity': '0.00001',
            'Probe Polarization': '90 90'}

        # >Settings for Charge Migration FT
        config['Charge Migration FT settings'] = {
            'Number of Omegas': '101',
            'Min Omega': '0.2',
            'Max Omega': '0.65',
            'Number of TauOmegas': '101',
            'Min TauOmega': '-0.20',
            'Max TauOmega': '0.20',
            'TimeStep (FT)': '3200',
            'WidthStep (FT)': '400'}

        # Write Configuration FIle with Defaults
        with open('chargemigration.ini', 'w') as configfile:
            config.write(configfile)

        # Exit to allow editing of Configuration File before running
        sys.exit()
