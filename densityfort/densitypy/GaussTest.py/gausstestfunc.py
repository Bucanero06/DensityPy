import numpy as np
import matplotlib.pyplot as plt
import decimal
import time
import fortranformat as ff
import math
from subprocess import Popen, PIPE, CalledProcessError


def Float_Range(start, stop, step):
    while start < stop:
        yield format(start, '.1f')  # float(start)
        start += decimal.Decimal(step)


def ExpDecay(startingvalue,time, decayfactor):
    if abs(decayfactor) > 1:
        decayfactor    = 1
    smoothing_function = 1 / np.sqrt(2 * np.pi) * np.exp(-(abs(time)) ** 2 / 2.)
    valueafterdecay    = (((1 - abs(decayfactor)) ** (abs(time))) * (1 * startingvalue))
    return valueafterdecay


def rFouriertransform(Dipole, tVec, OmegaVec):
    FT_Dipole= np.array([[0] * len(OmegaVec)] * len(OmegaVec), dtype=complex)
    for iOmega in range(0,len(OmegaVec)):
        w = OmegaVec[iOmega]
        for iSim in range(0,len(tVec)):
            t = tVec[iSim]
            zExpFact =  complex(np.cos(w*t),np.sin(w*t))
            FT_Dipole[:,iOmega] = FT_Dipole[:, iOmega] + (zExpFact * Dipole[:, iSim])

    dt = tVec[2] - tVec[1]
    FT_Dipole = FT_Dipole * (dt/(2 * np.pi))
    return FT_Dipole

# >Writes pulses to Field File ,Pump Probe Experiment (0)
def Write_Pulses(Field_File,Type_of_Pulse_Pump,Start_Time,Pump_Central_Frequency,
                           Pump_Periods,Pump_Phase,Pump_Intensity,Pump_Polarization,
                           Type_of_Pulse_Probe,Time_Delay_Range, Probe_Central_Frequency,
                           Probe_Periods, Probe_Phase, Probe_Intensity, Probe_Polarization):

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
        fout.write("XUV;}")#i ruben added "XUV;"

def Execute(command):
    #>Executes to command line
    with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
        for line in p.stdout:
            print(line, end='')  # process line here
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)