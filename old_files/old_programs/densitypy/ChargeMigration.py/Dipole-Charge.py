#! /usr/bin/env python3.6
from os import system
import argparse


def run(args):
    # Gets Run_Time_Parameters
    dipole_file = args.dipole
    charge_file = args.charge
    output_file = args.output_file

    system('rm ' + str(output_file))

    system('cat ' + str(
        dipole_file) + ' | awk \'{print $1\" \"$2\" \"$3**2+$4**2\" \"$5**2+$6**2\" \"$7**2+$8**2}\' > temp_dipole')
    system('cat ' + str(charge_file) + ' | awk \'{print $3**2+$4**2\" \"$5**2+$6**2\" \"$7**2+$8**2}\' > temp_charge')
    system(
        'paste temp_dipole temp_charge| awk \'{print $1\" \"$2\" \"$3\" \"$4\" \"$5\" \"$6\" \"$7\" \"$8\" \"$3-$6\" \"$4-$7\" \"$5-$8}\' > temp_difference')

    system("perl -pe \"s/0 0 0 0 0 0   0 0 0/ /g\" temp_difference>" + str(output_file))
    system('rm temp_charge temp_dipole temp_difference')


def main():
    parser = argparse.ArgumentParser(description="subtract charge from dipole")
    parser.add_argument("-d", help="dipoleft_ww input file.",
                        dest="dipole", type=str, action='store')
    parser.add_argument("-a", help="atomicchargeft_ww input file",
                        dest="charge", type=str, action='store')
    parser.add_argument("-o", help="outputfile", default="Dipole_Charge_difference",
                        dest="output_file", type=str, action='store')

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
###
# ')
