#! /usr/bin/env python3.7
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import plotly.graph_objects as go
import csv
import fnmatch
import os
from os import path
import moviepy.video.io.ImageSequenceClip
import subprocess
import shutil
import configparser
from configparser import ConfigParser
import numpy as np
import sys
import argparse

 
def run(args):
    ###Variables

    ChDen_timestepsdir = args.input
    images_folder = ChDen_timestepsdir + '/images_folder'
    tempfiles = ChDen_timestepsdir + '/tempfiles'
    ###End of Variables

    ##Creates directories to store temporary information
    if path.exists(tempfiles):
        path_to_dir = tempfiles
        files_in_dir = os.listdir(path_to_dir)  # get list of files in the directory
        for file in files_in_dir:  # loop to delete each file in folder
            os.remove(f'{path_to_dir}/{file}')
    else:
        os.mkdir(tempfiles)
    #
    if path.exists(images_folder):
        path_to_dir = images_folder
        files_in_dir = os.listdir(path_to_dir)  # get list of files in the directory
        for file in files_in_dir:  # loop to delete each file in folder
            os.remove(f'{path_to_dir}/{file}')
    else:
        os.mkdir(images_folder)

    ###End of Creatting Directories

    ##Plot
    ##Reads the last few digits of file name to organize them in the correct order
    def last_8chars(x):
        return (x[-8:])

    file_list = os.listdir(ChDen_timestepsdir)
    ##

    ##Finds all time files
    print('Done extracting Density Data \nNow Plotting')
    for file_name in sorted(file_list, key=last_8chars):
        if fnmatch.fnmatch(file_name, 'ChDen*'):

            ###Changes formating of columns from space to ',' and D to E to avoid bugs
            with open(ChDen_timestepsdir + '/' + file_name, 'r') as fin:
                with open(tempfiles + '/temp' + file_name, 'w') as fout:
                    for i in fin:
                        fout.write(i.replace('  ', ','))

            fin = open(tempfiles + '/temp' + file_name, "rt")
            data = fin.read()
            data = data.replace('D', 'E')
            fin.close()
            fin = open(tempfiles + '/temp' + file_name, "wt")
            fin.write(data)
            fin.close()

            print('Plotting: ' + file_name)
            x = []
            y = []
            z = []
            d = []

            ###Appends rows in time files
            with open(tempfiles + '/temp' + file_name, 'r') as fin:
                xyz = csv.reader(fin, delimiter=',')
                for row in xyz:
                    x.append(row[1])
                    y.append(row[2])
                    z.append(row[3])
                    d.append(row[4])

            ###########VolumePlot
            # ###Plots time files
            # fig = go.Figure()
            # color = 1
            # fig.add_trace(go.Volume(x=x,
            #                         y=y,
            #                         z=z,
            #                         colorscale=[[0, 'blue'], [0.48, 'blue'], [0.49, 'white'], [0.51, 'white'],
            #                                     [0.52, 'red'],[1.0, 'red']],
            #                         value=d,
            #                         isomax=max([float(i) for i in d]),
            #                         isomin=min([float(i) for i in d]),
            #                         opacity=0.001,
            #                         surface_count=1000,
            #                         cmax=color,
            #                         cmin=-color,
            #                         cmid=0
            #                         ))
            #
            # camera = dict(
            #     up=dict(x=0, y=1, z=0),
            #     center=dict(x=0, y=0, z=0),
            #     eye=dict(x=1.0, y=0.5, z=1.0)
            # )
            #
            # fig.update_layout(scene_camera=camera)
            # fig.write_image(images_folder + '/' + file_name + ".png")
            #######################End of Volume Plot

            #########Scatter Plott
            X = np.array(x)
            X = X[:, None]  # col
            np.set_printoptions(threshold=sys.maxsize)
            X = X.astype(float)

            Y = np.array(y)
            Y = Y[:, None]  # col
            np.set_printoptions(threshold=sys.maxsize)
            Y = Y.astype(float)

            Z = np.array(z)
            Z = Z[:, None]  # col
            np.set_printoptions(threshold=sys.maxsize)
            Z = Z.astype(float)

            D = np.array(d)
            D = D[:, None]  # col
            np.set_printoptions(threshold=sys.maxsize)
            D = D.astype(float)

            threshold = 0.005  # sets threshold for Point Value
            X = np.where(abs(D) < threshold, np.nan, X)
            Y = np.where(abs(D) < threshold, np.nan, Y)
            Z = np.where(abs(D) < threshold, np.nan, Z)
            D = np.where(abs(D) < threshold, np.nan, D)

            class MidpointNormalize(mpl.colors.Normalize):
                def __init__(self, vmin, vmax, midpoint=0, clip=False):
                    self.midpoint = midpoint
                    mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

                def __call__(self, value, clip=None):
                    normalized_min = max(0,
                                         1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
                    normalized_max = min(1,
                                         1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
                    normalized_mid = 0.5
                    x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
                    return np.ma.masked_array(np.interp(value, x, y))

            vmin = D.min()
            vmax = D.max()
            norm = MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=0)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(X, Y, Z, c=['blue' if i < 0 else 'red' for i in D], marker='o', norm=norm, alpha=1)
            ax.view_init(91.94805194805241, -88.344155844155416)
 
            ax.set_xlabel('X Label')
            ax.set_ylabel('Y Label')
            ax.set_zlabel('Z Label')

            plt.axis('off')
            plt.savefig(images_folder + '/' + file_name + ".png")
            plt.clf()

        #####End of Scatter Plot

        ###Resseting variables
        x = []
        y = []
        z = []
        d = []

    ###End of plot

    #

    ###Movie Maker
    ##Reads the last few digits of file name to organize them in the correct order
    def last_13chars(x):
        return (x[-13:])

    print('Done Plotting\nNow Making Clip')

    ##
    sorted_folder = sorted(os.listdir(images_folder + '/'), key=last_13chars)

    i = -1
    for filename in sorted_folder:
        i = i + 1
        os.rename(images_folder + '/' + filename, images_folder + '/'
                  + str(i) + filename)
        print('+ Frame ' + str(i))

    fps = 2
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(images_folder, fps=fps)
    clip.write_gif('ChargeMigration_NMA.GIF')
    ##



def main():
    parser = argparse.ArgumentParser(description="Help for density.py")
    parser.add_argument("-i", "--input", help="Directories with time step files",
                        dest="input", type=str,required=True)
    #parser.add_argument("-", help="", dest="", type=)
    #parser.add_argument("-", help="", dest="", type=)
    #parser.add_argument("-", help="", dest="", type=)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()















