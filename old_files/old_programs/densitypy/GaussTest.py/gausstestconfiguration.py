#! /usr/bin/env python3.6
#Made by Ruben
#>import configurationfileparser as rparse

import configparser
import sys
from configparser import ConfigParser
from os import path

def ParseConfigurationFile():
        config = ConfigParser()

        # Project Settings
        config['Project settings'] = {
            'Project_Name': 'CO2',
            'XYZ Molecule Geometry': 'CO2.xyz',
            'Basis':'ANO-L-MB', #only needs to be included here for use with -sas arguement. LUSUCS GUI
            # 'Number of Active Electrons':'10',
            'Number of States':'5',
            'List_of_Active_Orbitals': '7 8 9 10 11 12 13'}# Select which Orbitals to work with

        # >Settings for Grid (Number of points per axis and size of Cube)
        config['GRID settings'] = {
            'Number of Points, X axis': '15',
            'Number of Points, Y axis': '15',
            'Number of Points, Z axis': '15',
            'X MIN': '-4',
            'Y MIN': '-2',
            'Z MIN': '-2',
            'X MAX': '4',
            'Y MAX': '2',
            'Z MAX': '2'}

        # >Settings for Charge Migration
        config['Charge Migration settings'] = {
            'OutDirectory': 'sim',
            'Number of Times': '1001',
            'Min Time': '-500',
            'Max Time': '500',
            # 'Field File': 'Field',
            'TimeStep (FT)': '800',
            'WidthStep (FT)': '50'}

        # >Pump Settings
        config['Pulses settings (Pump)'] = {
            'Type of Pulse': 'G',
            'Start_Time': '0.000',
            'Pump Central Frequency': '0.30',
            'Pump Periods': '5',
            'Pump Phase': '180',
            'Pump Intensity': '0.001',
            'Pump Polarization': '90 0'}
        
        # >Probe Settings
        config['Pulses settings (Probe)'] = {
            'Type of Pulse': 'G',
            'Time Delay Start': '-50',
            'Time Delay Stop': '100',
            'Number Of PP': '151',
            'Probe Central Frequency': '0.30',
            'Probe Periods': '5',
            'Probe Phase': '0',
            'Probe Intensity': '0.001',
            'Probe Polarization': '90 0'}

        # Write Configuration FIle with Defaults
        with open('chargemigration.ini', 'w') as configfile:
            config.write(configfile)

        # Exit to allow editing of Configuration File before running
        sys.exit()