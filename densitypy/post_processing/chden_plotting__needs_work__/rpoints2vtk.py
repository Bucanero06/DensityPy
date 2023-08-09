#! /usr/bin/env python3.6
from os import path
import csv
from subprocess import Popen, PIPE, CalledProcessError
import argparse
import glob
import pandas as pd

import time


def Execute(command):
    # >Executes to command line
    with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
        for line in p.stdout:
            print(line, end='')  # process line here
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)


########################################################################################################################
import os


def linecount_wc(fin):
    return int(os.popen('wc -l ' + fin).read().split()[0])


def linecount_1(fin):
    return len(open(fin).readlines())


def Read_Coordinates(filename):
    coords = [] * 3
    x = []
    y = []
    z = []
    with open(filename, 'r')as fin:
        N_points = sum(not line.isspace() for line in fin)
    print("Number of Points = " + str(N_points))

    with open(filename, 'r')as fin:
        stop = 0
        for line in fin:
            stop = stop + 1
            x_, y_, z_, d_ = line.split()
            x.append([float(x_)])
            y.append([float(y_)])
            z.append([float(z_)])
            if (stop == N_points):
                break
    for i in range(0, N_points):
        coords.append((x[i], y[i], z[i]))
    return coords, N_points


def FindDELTA(coords, N_points):
    threshold = 0.0000000001
    DELTA = []
    for iCoord in range(0, 3):
        coord1 = coords[0][iCoord]
        coord1 = (float(str(coord1)[1:-1]))
        delta = 100000
        for iPts in range(0, N_points):
            coord2 = coords[iPts][iCoord]
            coord2 = (float(str(coord2)[1:-1]))
            if (abs(coord1 - coord2) > threshold):
                delta = min(delta, abs(coord1 - coord2))
        DELTA.append(delta)
    return DELTA


def Find_Max_and_Min(coords):
    maxi = []
    mini = []
    for iPol in range(0, 3):
        max = 0
        for i in range(0, len(coords)):
            max_ = float(str(coords[i][iPol])[1:-1])
            if max_ > max:
                max = max_
        min = 0
        for i in range(0, len(coords)):
            min_ = float(str(coords[i][iPol])[1:-1])
            if min_ < min:
                min = min_
        maxi.append(max)
        mini.append(min)
    return maxi, mini


########################################################################################################################

def run(args):
    from natsort import realsorted
    import progressbar

    # Gets Run_Time_Parameters
    inputdirectory = args.inputdirectory

    outputdirectory = args.outputdirectory
    if outputdirectory == None:
        outputdirectory = inputdirectory + "vtk"

        # Create out directory
    if path.exists(str(outputdirectory)):
        Execute("rm -rf " + str(outputdirectory))
        print("Deleted old " + str(outputdirectory) + "directory")
        Execute("mkdir " + str(outputdirectory))
        print("Created empty " + str(outputdirectory) + "directory")
    else:
        Execute("mkdir " + str(outputdirectory))
        print("Created empty " + str(outputdirectory) + "directory")

    # List Files in Correct Order
    list = glob.glob1(str(inputdirectory) + "/", "ChDen*")
    ChDenlisted_folder = realsorted(list)

    ####################################################################################################################
    # > Read file to get GridCoord and Count number of points
    gridcoord, number_of_points = Read_Coordinates(inputdirectory + "/" + ChDenlisted_folder[0])
    # > Get DELTA for each axis
    DELTA = FindDELTA(gridcoord, number_of_points)
    # > Get Max and Min of each axis
    maxi, mini = Find_Max_and_Min(gridcoord)
    xmax = maxi[0]
    ymax = maxi[1]
    zmax = maxi[2]
    xmin = mini[0]
    ymin = mini[1]
    zmin = mini[2]
    # > Find Number of Points per Axis
    X_axis = round((abs(xmax - xmin)) / DELTA[0]) + 1
    Y_axis = round((abs(ymax - ymin)) / DELTA[1]) + 1
    Z_axis = round((abs(zmax - zmin)) / DELTA[2]) + 1
    print(X_axis, Y_axis, Z_axis)
    dimensionsofgrid = str(str(X_axis) + " " + str(Y_axis) + " " + str(Z_axis))
    ####################################################################################################################

    print("Number of Files inside " + inputdirectory + " = " + str(len(ChDenlisted_folder)) +
          "\nNumber of Points per file = " + str(number_of_points))

    # # Separate Values in File using CSV format
    # print("Separating Values in File Using CSV format")
    # for filename in ChDenlisted_folder:
    #     with open(inputdirectory + "/" + str(filename), 'r') as fin:
    #         with open(outputdirectory + "/csv_" + str(filename), 'w') as fout:
    #             for i in fin:
    #                 fout.write(i.replace('   ', ','))
    # print("Done Separating Values in File Using CSV format")

    # Extract Values
    frame = 0
    for filename in progressbar.progressbar(ChDenlisted_folder, redirect_stdout=True):

        print("Creating VTK file from " + str(filename))
        x = []
        y = []
        z = []
        d = []
        csv.register_dialect('skip_space', skipinitialspace=True)
        with open(inputdirectory + "/" + str(filename), 'r') as f:
            reader = csv.reader(f, delimiter=' ', dialect='skip_space')
            for row in reader:
                if row == "\n":
                    next(reader)
                x.append(row[0])
                y.append(row[1])
                z.append(row[2])
                d.append(row[3])

        # Write VTK file
        # with open(outputdirectory + "/" + str(filename) + "_time" + str(frame) + ".vtk", 'w') as fout:
        with open(outputdirectory + "/ChDen" + str(frame) + ".vtk", 'w') as fout:
            fout.write("# vtk DataFile Version 2.0 "
                       "\nTotal electron density"
                       "\nASCII"
                       "\nDATASET STRUCTURED_GRID"
                       "\nDIMENSIONS " + str(dimensionsofgrid) +
                       "\nPOINTS " + str(number_of_points) + " float")
            for i in range(0, number_of_points):
                fout.write("\n{} {} {}".format(x[i], y[i], z[i]))

            fout.write("\nPOINT_DATA " + str(number_of_points) +
                       "\nSCALARS Density float 1"
                       "\nLOOKUP_TABLE my_table")

            for i in range(0, number_of_points):
                fout.write("\n" + d[i])
        frame += 1

        # time.sleep(0.5)


###############################
def main():
    parser = argparse.ArgumentParser(description="Converts xyz points and density output to vtk format.")

    parser.add_argument("-i", help="Name or path to directory containing the ChDen* files with the ChargeDensity. "
                                   "timesteps for the desired pulse sequence",
                        dest="inputdirectory", type=str, required=True, action="store")
    parser.add_argument("-o", help="Output directory for vtk files",
                        dest="outputdirectory", type=str, required=False, action="store")

    # parser.add_argument("-", help="",
    #                     dest="", action='store_true')
    # parser.add_argument("-", help="",
    #                     dest="", action='store_true')
    # parser.add_argument("-", help="",
    #                     dest="", action='store_true')
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
