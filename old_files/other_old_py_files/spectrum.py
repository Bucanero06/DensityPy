#! /usr/bin/env python3.6
from os import system
import sys
from configparser import ConfigParser

parser = ConfigParser()
parser.read('chargemigration.ini')

project = parser.get('Project settings', 'Project_Name')
sim_folder = parser.get('Charge Migration settings', 'Output Directory')

print("Running SpectrumReconstruction ...")
system('SpectrumReconstruction -i ' + project + '_output -o ' + sim_folder + ' -xyz ' + project + '.xyz -nw 101 '
                                                                                                  '-wmax 0.65 -wmin 0.2 -ntw 101 -twmax 0.20 -twmin -0.2')
print("Running Dipole-Charge.py ...")
system(
    'Dipole-Charge.py -d ' + sim_folder + '/Dipole/DipoleFT_ww -a ' + sim_folder + '/Dipole/DipoleFT_ww_reconstructed -o difference')
print("Done")

###
