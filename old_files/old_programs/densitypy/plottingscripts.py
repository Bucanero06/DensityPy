#! /usr/bin/env python3.6

import argparse
import os
import glob
from natsort import realsorted
from def_functions import Execute,Execute2


# Main
def run(args):
    #Gets Run_Time_Parameters
    inputdirectory     = args.inputdirectory
    gnuplotscript      = args.gnuplotscript
    nvalues            = args.nvalues
    mencodertrue       = args.mencodertrue
    fps                = args.fps
    output             = args.output

    if inputdirectory:
        if os.path.exists("PNG" + str(inputdirectory)):
            Execute("rm -rf PNG" + str(inputdirectory))
            Execute("mkdir PNG" + str(inputdirectory))
        else:
            Execute("mkdir PNG" + str(inputdirectory))

        # Do from sim directory
        list = glob.glob1("ChargeDensity/" + str(inputdirectory) + "/","ChDen*")
        listed_folder= realsorted(list)

        with open("listfiles.txt",'w') as fout:
            for filename in listed_folder:
                fout.write("PNG" + str(inputdirectory) + "/Z0" + filename + ".png\n")
                Execute("Split3DGrid -i ChargeDensity/" + str(inputdirectory) + "/" + filename +
                        " -o PNG" + str(inputdirectory) + "/Z0" + filename + " -n " + str(nvalues))
                Execute("gnuplot -e \"file=\'PNG" + str(inputdirectory) + "/Z0" + filename + "\'\" " + gnuplotscript)
        #Do from sim Directory

    if mencodertrue:
        Execute2("apt-cache policy mencoder", "running mencoder", "mencoder not found or not installed")
        Execute("mencoder mf://@listfiles.txt -mf w=600:h=600:fps=" + str(fps) +
                ":type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o " + output + ".avi")

def main():
    parser = argparse.ArgumentParser(description="Prepares grid for gnuplot, creates png files for each timestep "
                                                 "and compiles them into a video.avi file. Must be ran from inside "
                                                 "the simulation file.")

    # parser.add_argument("-i", "--input", help="Directories with time step files",
    #                     dest="input", type=str, required=True)
    parser.add_argument("-i", help="Name, NOT path of directory inside of ChargeDensity folder containing the ChDen "
                                   "timesteps for the desired pulse sequence",
                        dest="inputdirectory", type=str ,required=False, action="store")
    parser.add_argument("-g", help="Script to be used to make PNG files with GNUPLOT",
                        dest="gnuplotscript", type=str, required=False)
    parser.add_argument("-n", help="number of values in file (e.g., X,Y,Z,n1,n2,n3,n4 ==> -n 4)",
                        dest="nvalues", type= int, required=False)
    parser.add_argument("-mplot", help="whether or not to plot using mencoder.",
                        dest="mencodertrue", action='store_true', required=False)
    parser.add_argument("-fps", help="frames per seconds",
                        dest="fps",  type=int, default="20", required=False)
    parser.add_argument("-o", help="OutPut files",
                        dest="output", type=str, default= "OutPut", required=False)
    # parser.add_argument("-", help="",
    #                     dest="", action='store_true')
    # parser.add_argument("-", help="",
    #                     dest="", action='store_true')
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
###