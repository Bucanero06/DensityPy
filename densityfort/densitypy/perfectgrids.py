#! /usr/bin/env python3.6
from os import system
import argparse




def make_coor(n):
    phi, theta = np.mgrid[0:pi-0:n, 0:2 * pi:n]
    Coor = namedtuple('Coor', 'r phi theta x y z')
    r = 1
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    return Coor(r, phi, theta, x, y, z)

pts=make_coor(15j)

mlab.figure()
mlab.points3d(pts.x,pts.y,pts.z, scale_factor=0.1)
mlab.mesh(pts.x,pts.y,pts.z)
mlab.savefig('spheretest.png')
mlab.show()


# def run(args):
#     # Gets Run_Time_Parameters
#     dipole_file = args.dipole
#     charge_file = args.charge
#     output_file = args.output_file
#
#     system('rm '+str(output_file))
#
# def main():
#     parser = argparse.ArgumentParser(description="subtract charge from dipole")
#     parser.add_argument("-d", help="dipoleft_ww input file.",
#                         dest="dipole", type=str, action='store')
#     parser.add_argument("-a", help="atomicchargeft_ww input file",
#                         dest="charge", type=str, action='store')
#     parser.add_argument("-o", help="outputfile",default="Dipole_Charge_difference",
#                         dest="output_file", type=str, action='store')
#
#     parser.set_defaults(func=run)
#     args = parser.parse_args()
#     args.func(args)
#
#
# if __name__ == "__main__":
#     main()
###
# ')