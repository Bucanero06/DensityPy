#! /usr/bin/env python3.6
import numpy as np
import matplotlib.pyplot as plt
import decimal
import time
import fortranformat as ff
import math
import os
import sys

import gausstestfunc as gfunc
from configparser import ConfigParser

# lineformat = ff.FortranRecordWriter('(1E26.16)')
#
# # >Load Values
# parser = ConfigParser()
# parser.read('chargemigration.ini')
#
# # Charge Migration Parameters
# Number_of_Times = parser.getint('Charge Migration settings', 'Number of Times')
# Min_Time = parser.getint('Charge Migration settings', 'Min Time')
# Max_Time = parser.getint('Charge Migration settings', 'Max Time')
# Change_In_Time = (Max_Time - Min_Time) / (Number_of_Times - 1)
# Time_Range = list(gfunc.Float_Range(Min_Time, Max_Time + Change_In_Time, Change_In_Time))
#
# # Pump Settings
# Type_of_Pulse_Pump = parser.get('Pulses settings (Pump)', 'Type of Pulse')
# Start_Pump_Time = parser.getfloat('Pulses settings (Pump)', 'Start_Time')
# Pump_Central_Frequency = parser.getfloat('Pulses settings (Pump)', 'Pump Central Frequency')
# Pump_Periods = parser.getfloat('Pulses settings (Pump)', 'Pump Periods')
# Pump_Phase = parser.getfloat('Pulses settings (Pump)', 'Pump Phase')
# Pump_Intensity = parser.getfloat('Pulses settings (Pump)', 'Pump Intensity')
# Pump_Polarization = parser.get('Pulses settings (Pump)', 'Pump Polarization')
# # Probe Settings
# Type_of_Pulse_Probe = parser.get('Pulses settings (Probe)', 'Type of Pulse')
# Time_Delay_Start = parser.getint('Pulses settings (Probe)', 'Time Delay Start')
# Time_Delay_Stop = parser.getint('Pulses settings (Probe)', 'Time Delay Stop')
# Number_Of_PP = parser.getint('Pulses settings (Probe)', 'Number Of PP')
# Probe_Central_Frequency = parser.getfloat('Pulses settings (Probe)', 'Probe Central Frequency')
# Probe_Periods = parser.getfloat('Pulses settings (Probe)', 'Probe Periods')
# Probe_Phase = parser.getfloat('Pulses settings (Probe)', 'Probe Phase')
# Probe_Intensity = parser.getfloat('Pulses settings (Probe)', 'Probe Intensity')
# Probe_Polarization = parser.get('Pulses settings (Probe)', 'Probe Polarization')
# Change_In_Delay = (Time_Delay_Stop - Time_Delay_Start) / (Number_Of_PP - 1)
# Time_Delay_Range = list(gfunc.Float_Range(Time_Delay_Start, Time_Delay_Stop + Change_In_Delay, Change_In_Delay))
#
# # ..Parameters of Gaussian
# timemean = 0
# sigma_time = 40
# frequency = 1.2
# frequencywithdelay = 1.2
# Time_Axis = np.linspace(Min_Time, Max_Time, Number_of_Times)
# oscillation = np.cos(frequency * Time_Axis)
# Decay_Factor = 0.05
#
# # * (1 / (2 * np.pi * sigma_time)))
#
# # >Calls
# if os.path.exists("GaussTestDipole"):
#     gfunc.Execute("rm -r GaussTestDipole")
# os.mkdir("GaussTestDipole")
# time.sleep(1)
#
# for iPulse in Time_Delay_Range:
#     print("writing DipolePP" + iPulse)
#     oscillation_delay = np.cos(frequencywithdelay * float(iPulse))
#     # decay = gfunc.ExpDecay(1, float(iPulse), Decay_Factor) * oscillation_delay
#     smoothdecay = np.exp(-0.0015 * (float(iPulse))**2) * oscillation_delay
#
#     z = (np.exp(-((Time_Axis - float(iPulse)) ** 2 / (2 * sigma_time ** 2))) * oscillation)
#     Z = (z/10000) * smoothdecay #decay
#     with open('GaussTestDipole/DipolePP' + iPulse, 'w') as fout:
#         for iTime in range(0, len(Time_Axis)):
#             fout.write(
#                 "{} {} {} {} {} {} {} {}\n".format(str(iTime), float(Time_Range[iTime]), (float(Z[iTime])), 0.0000,
#                                                    0.0000, 0.0000, 0.0000, 0.0000))
# print('Done Writing Dipoles')
#
# # >Write Pulse File
# gfunc.Write_Pulses("Field_pulses_Test", Type_of_Pulse_Pump, Start_Pump_Time, Pump_Central_Frequency,
#                    Pump_Periods, Pump_Phase, Pump_Intensity, Pump_Polarization,
#                    Type_of_Pulse_Probe, Time_Delay_Range, Probe_Central_Frequency,
#                    Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization)
# print("Done Writing Field_pulses_Test File")
#
# gfunc.Execute("dChargeMigrationFT -i CO2_output -o test  -xyz CO2.xyz -stept 1000 -stepw 100 -field Field_pulses_Test -nw 101 -wmax 1.5 -wmin 0.9 -ntw 101 -twmax 2.2 -twmin 0.2")


# a= np.linspace(-100,100,101)
# b= np.exp(-0.0005 * a**2)
#
# plt.plot(a,b)
# plt.show()

for i in range(1,10):
    gfunc.Execute("dChargeMigrationFT -i NMA_output -o sim_" + str(i) + " -xyz NMA.xyz -stept 1000 -stepw 100 -field Field_pulses -nw 101 -wmax 0.7 -wmin 0 -ntw 101 -twmax 0.20 -twmin -0.20  ")