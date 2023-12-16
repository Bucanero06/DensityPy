#! /usr/bin/env python3.6

# def GetValueAfterString(filename,string,delimiter, afterline):
#     with open(filename, 'r') as fin:
#         for num, line in enumerate(fin):
#             if num > afterline:
#                 if string in line:
#                     option_value = (line.partition(delimiter)[2]).strip()
#         return option_value
#
#
#
# def FindLineForTitle0(filename,string,delimiter):
#     with open(filename, 'r') as fin:
#         reset_lines = []
#         reset_lines.append(0)
#         for num, line in enumerate(fin, 1):
#             if string in line:
#                 value = (line.partition(delimiter)[2]).strip()
#                 if value == '0':
#                     reset_lines.append(num)
#         return reset_lines
#                 # elif value == orbital:
#                 #     print("same")
#                 # option_value.append(value)
#         # return option_value
# /home/ruben/Documents/fermi/home/git/mojo/programs/densitypy/configurationfileparser.py
# def FormatDipoleFiles(numberofstates):
#     range_states = range(0, numberofstates, 1)
#     with open('../CO2_output/X_DIPOLE', 'r') as fin:
#         with open('../CO2_output/editedX_dipole', 'w') as fout:
#             for states in range_states:
#                 states += 1
#                 string = " " + str(states) + " "
#                 # stringstop = ["Title=", str(states + 1)]
#                 value = GetValueAfterString("../CO2_output/workingCO2Dipole", string, "      ")
#                 fout.write( ' '.join([str(f) for f in value]) + "\n")
# def GetValueAfterString(filename,string,delimiter):
#     with open(filename, 'r') as fin:
#         value=[]
#         for line in fin:
#             if string in line:
#                 option_value = (line.partition(delimiter)[2]).strip()
#                 value.append(option_value)
#         return value

# def Find(filename, *args):
#     directories=[*args]
#     foundfile = False
#     for searchdirectory in directories:
#         if path.exists(searchdirectory + "/" + filename):
#             foundfile = True
#             if searchdirectory == ".":
#                 print("Found " + str(filename) + " inside the current directory")
#             else:
#                 print("Found " + str(filename) + " inside " + str(searchdirectory) + " directory")
#             return searchdirectory
#             exit()
#         else:
#             pass
#     # if not exited by now it means that the file was not found in any of the given directories thus rise error
#     if foundfile != True:
#         print(str(filename) + " not found inside " + str(directories) + "\n exiting...")
#         sys.exit()
#
#
# Project_Name = "CO2"
#
# #Sanity Check
# pathtomolcasrc = path.expanduser("~/.Molcas/molcasrc")
# if path.exists(pathtomolcasrc):
#     molcas_workdir = rfunc.GetValueOfAsString(pathtomolcasrc, "MOLCAS_WORKDIR", "=")
# else:
#     print("Molcas configuration file ~/.Molcas/molcasrc not found. Use pymolcas --setup to create one")
#     sys.exit()
#
#
#
# logfilepath = rfunc.Find(Project_Name + ".log", ".", molcas_workdir)
# import time
# import progressbar


# f=[1,2,3,4,5,6,7,8,9,11,14]
#
# with progressbar.ProgressBar(max_value=18) as bar:
#     for i in f:
#         print(f[i])


# for i in progressbar.progressbar(f, redirect_stdout=True):
#     print(i)
#     # time.sleep(0.1)

# >Writes pulses to Field File ,Pump Probe Experiment (0)
# def Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
#                  Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
#                  Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
#                  Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization):
#     with open(Field_File, 'w') as fout:
#         # >Writes Pump Pulse
#         fout.write("[XUV]{( " + str(Type_of_Pulse_Pump) + " " +
#                    str(Start_Time) + " " +
#                    str(Pump_Central_Frequency) + " " +
#                    str(Pump_Periods) + " " +
#                    str(Pump_Phase) + "  " +
#                    str(Pump_Intensity) + " " +
#                    str(Pump_Polarization) + " );}")
#         for Time_Delay in Time_Delay_Range:
#             # >Writes Probe Pulse
#             fout.write("\n[PP" + str(Time_Delay) + "]{ XUV; ( " + str(Type_of_Pulse_Probe) + " " +
#                        str(Time_Delay) + " " +
#                        str(Probe_Central_Frequency) + " " +
#                        str(Probe_Periods) + " " + str(Probe_Phase) + " " +
#                        str(Probe_Intensity) + " " + " " + Probe_Polarization + " );}")
#         fout.write("\nEXECUTE{")
#         for Time_Delay in Time_Delay_Range:
#             fout.write("PP" + str(Time_Delay) + "; ")
#         fout.write("XUV;}")  # i ruben added "XUV;"
#
#
#     # Pump Settings
#     Type_of_Pulse_Pump = parser.get('Pulses settings (Pump)', 'Type of Pulse')
#     Start_Time = parser.getfloat('Pulses settings (Pump)', 'Start_Time')
#     Pump_Central_Frequency = parser.getfloat('Pulses settings (Pump)', 'Pump Central Frequency')
#     Pump_Periods = parser.getfloat('Pulses settings (Pump)', 'Pump Periods')
#     Pump_Phase = parser.getfloat('Pulses settings (Pump)', 'Pump Phase')
#     Pump_Intensity = parser.getfloat('Pulses settings (Pump)', 'Pump Intensity')
#     Pump_Polarization = parser.get('Pulses settings (Pump)', 'Pump Polarization')
#     # Probe Settings
#     Type_of_Pulse_Probe = parser.get('Pulses settings (Probe)', 'Type of Pulse')
#     Time_Delay_Start = parser.getint('Pulses settings (Probe)', 'Time Delay Start')
#     Time_Delay_Stop = parser.getint('Pulses settings (Probe)', 'Time Delay Stop')
#     Number_Of_PP = parser.getint('Pulses settings (Probe)', 'Number Of PP')
#     Probe_Central_Frequency = parser.getfloat('Pulses settings (Probe)', 'Probe Central Frequency')
#     Probe_Periods = parser.getfloat('Pulses settings (Probe)', 'Probe Periods')
#     Probe_Phase = parser.getfloat('Pulses settings (Probe)', 'Probe Phase')
#     Probe_Intensity = parser.getfloat('Pulses settings (Probe)', 'Probe Intensity')
#     Probe_Polarization = parser.get('Pulses settings (Probe)', 'Probe Polarization')
#     Change_In_Delay = (Time_Delay_Stop - Time_Delay_Start) / (Number_Of_PP - 1)
#     Time_Delay_Range = list(rfunc.Float_Range(Time_Delay_Start,
#                                               Time_Delay_Stop + Change_In_Delay, Change_In_Delay))
#
#
# Write_Pulses(Field_File, Type_of_Pulse_Pump, Start_Time, Pump_Central_Frequency,
#              Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
#              Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
#              Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization)

import pandas as pd
import sys
from os import remove, path
import numpy as np


def Read_xyz(xyz_file):
    atoms = []
    x = []
    y = []
    z = []
    atomic_coordinates = [] * 3
    with open(xyz_file, 'r') as fin:
        n_atoms = int(fin.readline())
        title = fin.readline()
        stop = 0
        for line in fin:
            stop = stop + 1
            atom_, x_, y_, z_ = line.split()
            atoms.append(atom_)
            x.append([float(x_)])
            y.append([float(y_)])
            z.append([float(z_)])
            if (stop == n_atoms):
                break

    for i in range(0, len(atoms)):
        atomic_coordinates.append((x[i], y[i], z[i]))

    print("filename:         %s" % xyz_file)
    print("title:            %s" % title)
    print("number of atoms:  %d" % n_atoms)

    return atomic_coordinates

#
# xyz_file = "MgP_.xyz"
# coords, atoms = Read_xyz(xyz_file)
# for i in coords:
#     print(i)
