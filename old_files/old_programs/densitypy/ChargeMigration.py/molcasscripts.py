#! /usr/bin/env python3.6
# >import molcasscripts as rmolcs

import def_functions as rfunc
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys
from os import remove, path
import numpy as np


# Creates Helpful input file for Molcas with comments to help user begin (0)
def CreateHelpInputFile():
    # >Creates helpful input file guide for pymolcas
    with open('inputhelp.input', 'w') as fout:
        fout.write("&GATEWAY"
                   " \n Title = _name_of_project_"
                   " \n coord = _geometry_file_"
                   " \n basis = _name_of_basis_"
                   " \n group = c1 //this prevents use of symmetry"
                   "\n&SEWARD\n&SCF"
                   "\n&RASSCF"
                   "\n Ras1   = _number_of_RAS1_orbitals_"
                   "\n Ras2   = _number_of_RAS2_orbitals_"
                   "\n RAS3   = _number_of_RAS3_orbitals_"
                   "\n Inactive = _number_of_Inactive_orbitals_"
                   "\n nactel   = _total_ _RAS1_ _RAS3_ "
                   "//# of active electrons, allowed holes, allowed excitations"
                   "\n TDM       //writes DM and TDM to RASSCF.h5"
                   "\n CIRoot\n _highest_root_needed_ _dimention_of_CI_matrix_ _weight_ "
                   "//e.g. 10 10 1 will do a state average calculation of the ground state "
                   "and the lower 4 excited states"
                   "\n&RASSI"
                   "\nNr of JobIph"
                   "\n 1 ALL "
                   "    //e.g. # of job files, # of states in that file ==>> 1 10 or 1 ALL "
                   "\nMEES //writes one-electron properties")
    print("inputhelp.input created")
    sys.exit()


# >Makes points to be used for pymolcas input (0)
def Make_Grid_Coordinates(Molcas_Directory, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax):
    # >Makes points to be used for pymolcas input
    N_points = Nx * Ny * Nz
    print("Number of Points = " + str(N_points))
    with open(Molcas_Directory + '/gridcoord', 'w') as fout:
        for i in range(0, Nx):
            x = xmin + (xmax - xmin) / (Nx - 1) * i
            for j in range(0, Ny):
                y = ymin + (ymax - ymin) / (Ny - 1) * j
                for k in range(0, Ny):
                    z = zmin + (zmax - zmin) / (Nz - 1) * k
                    fout.write("{} {} {}\n".format(x, y, z))


##################################################################################
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


def distance(avec, bvec):
    c = np.square(np.subtract(avec, bvec))
    cvec = np.sqrt(np.sum(c))
    return cvec


def Make_New_Grid(atomic_coordinates, step_size, Boundary, limitedgrid):
    # > Define coordinates matrix
    gridcoord = [] * 3

    maxi, mini = Find_Max_and_Min(atomic_coordinates)
    # # > Find Max Boundries of Grid
    xmax = maxi[0] + Boundary
    ymax = maxi[1] + Boundary
    zmax = maxi[2] + Boundary
    xmin = mini[0] - Boundary
    ymin = mini[1] - Boundary
    zmin = mini[2] - Boundary

    x_ = np.arange(xmin, xmax, step_size)
    y_ = np.arange(ymin, ymax, step_size)
    z_ = np.arange(zmin, zmax, step_size)
    x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

    initial_Npoints = len(x_) * len(y_) * len(z_)
    current_point_index = 0
    for i in range(0, len(x_)):
        for j in range(0, len(y_)):
            for k in range(0, len(z_)):
                X = x[i][j][k]
                Y = y[i][j][k]
                Z = z[i][j][k]
                r_point = [X, Y, Z]
                current_point_index = current_point_index + 1
                percentage_done = round((100 * (current_point_index / initial_Npoints)), 3)
                print('\r' + str(percentage_done) + ' %', end='  ')

                if limitedgrid == True:
                    for icoord in atomic_coordinates:
                        avec_new = []
                        avec_new.append(float(str(icoord[0])[1:-1]))
                        avec_new.append(float(str(icoord[1])[1:-1]))
                        avec_new.append(float(str(icoord[2])[1:-1]))

                        if abs(distance(avec_new, r_point)) <= Boundary:
                            gridcoord.append([X, Y, Z])
                            break
                else:
                    gridcoord.append([X, Y, Z])
    print("Done Making Grid")
    print("\n")
    return gridcoord, initial_Npoints


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


def Make_Better_Grid(Molcas_Directory, geometry, step_size, Boundary, limitedgrid, ):
    print("Making Grid")
    atomic_coordinates = Read_xyz(geometry)
    # > Make_Grid(maxi, mini, Npoints)
    gridcoord, initial_Npoints = Make_New_Grid(atomic_coordinates, step_size, Boundary, limitedgrid)
    N_points = len(gridcoord)
    print("Number of Points = " + str(N_points))
    # > Writes to file
    print("Writing Grid to File")
    with open(Molcas_Directory + '/gridcoord', 'w') as fout:
        for i in range(0, N_points):
            fout.write(
                "{} {} {}\n".format(round(gridcoord[i][0], 2), round(gridcoord[i][1], 2), round(gridcoord[i][2], 2)))
    return N_points


#################################################################################################
# Creates input from users own input file for pymolcas(0)
def CopyInputFileToEdit(pymolcas_input, Project_Name):
    # Creates input from users own input file for pymolcas
    string_list = ['RASSCF', 'RASSI', 'TDM']
    true_values = []

    with open(pymolcas_input, 'r') as fin:
        #
        with open(Project_Name + '.input', 'w') as fout:
            print('Creating ' + Project_Name + '.input from ' + pymolcas_input)
            for line in fin:
                if line.startswith('//') is False:
                    fout.write(line)
                    for string in string_list:
                        if string.upper() in line.upper():
                            true_values.append(string.upper())
    #
    return true_values


# Adds GRID_IT section to copied input file
def AddGridItToManualInputFile(pymolcas_input, Project_Name, Molcas_Directory, orbital_list, N_points):
    # Creates input from users own input file for pymolcas
    Lowest_orbital = min(orbital_list)
    Highest_orbital = max(orbital_list)
    with open(Molcas_Directory + '/gridcoord', 'r') as fin:
        with open(Project_Name + '.input', 'a') as fout:
            fout.write('&Grid_it '
                       '\n select '
                       '\n 1:' + str(Lowest_orbital) + '-' + str(Highest_orbital) +
                       '\n NOLUSCUS '
                       '\n GRID'
                       '\n ')
            fout.write(str(N_points))
            fout.write("\n" + fin.read())


# >Calls pymolcas and writes to scrach folder and log file(1)
def Call_OpenMolcas(Project_Name):
    # >Calls pymolcas and writes to scrach folder and log file
    print("Running OpenMolcas")
    rfunc.ExecutePymolcasWithErrorPrint("pymolcas " + Project_Name + ".input -f", Project_Name)
    # rfunc.Execute("rm " + Project_Name + ".input")


# >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors(1,1,2)
def FetchFromh5_File(Project_Name, Molcas_Directory, pathtofile, filename):
    # >Extracts Information from H5FILEs
    rfunc.ExecuteNoWrite("h5dump -o " + Molcas_Directory + "/" + filename +
                         " -y -w 0 -d " + filename + " " + pathtofile + "/" + Project_Name + ".rasscf.h5")
    print("Extracted " + filename)


def FetchFromh5_Filewidth1(Project_Name, Molcas_Directory, pathtofile, filename):
    # >Extracts Information from H5FILEs
    rfunc.ExecuteNoWrite("h5dump -o " + Molcas_Directory + "/" + filename +
                         " -y -w 1 -d " + filename + " " + pathtofile + "/" + Project_Name + ".rasscf.h5")
    print("Extracted " + filename)


def LoadFromh5File(Project_Name, Molcas_Directory, workingdirectory, true_values, density_matrix,
                   transition_density_matrix,
                   root_energies, ao_multiple_x, ao_multiple_y, ao_multiple_z,
                   mo_energies, mo_vectors, justh5):
    # >Find <ProjectName>.rasscf.h5
    rasscf_h5_filepath = rfunc.Find(Project_Name + ".rasscf.h5", ".", workingdirectory)

    # >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors
    FetchFromh5_File(Project_Name, Molcas_Directory, rasscf_h5_filepath, density_matrix)
    if "TDM" in true_values or justh5:
        FetchFromh5_File(Project_Name, Molcas_Directory, rasscf_h5_filepath, transition_density_matrix)
    else:
        print(f'Found {Project_Name}".rasscf.h5" but TDM keyword was not found in Molcas input file, '
              f'thus no Transition Density Matrix printed')
    FetchFromh5_Filewidth1(Project_Name, Molcas_Directory, rasscf_h5_filepath, root_energies)
    FetchFromh5_File(Project_Name, Molcas_Directory, rasscf_h5_filepath, ao_multiple_x)
    FetchFromh5_File(Project_Name, Molcas_Directory, rasscf_h5_filepath, ao_multiple_y)
    FetchFromh5_File(Project_Name, Molcas_Directory, rasscf_h5_filepath, ao_multiple_z)
    FetchFromh5_Filewidth1(Project_Name, Molcas_Directory, rasscf_h5_filepath, mo_energies)
    FetchFromh5_Filewidth1(Project_Name, Molcas_Directory, rasscf_h5_filepath, mo_vectors)


# >Extracts Density for the grids of each Orbital found oin <Project>.grid (0)
def ExtractGridDensity(orbital_list, output_filename, directorypath, Molcas_Directory):
    Highest_orbital = max(orbital_list)
    stringend = "Title=    0 "
    reset_lines = []
    reset_lines.append(0)
    gridfilepath = rfunc.Find(output_filename + '.grid', ".", directorypath)
    with open(gridfilepath + "/" + output_filename + '.grid', 'r') as fin:
        for num, line in enumerate(fin, 1):
            if stringend in line:
                reset_lines.append(num)

    for iorb in range(len(orbital_list)):
        orbital = orbital_list[iorb]
        string = ["Title=", str(orbital)]

        if (orbital == Highest_orbital):
            pass
        else:
            next_orbital = orbital_list[iorb + 1]
            stringstop = ["Title=", str(next_orbital)]

        # >This commented section is here in case the fix above does not work
        # for orbital in orbital_list:
        #     string = ["Title=", str(orbital)]
        #     stringstop = ["Title=", str(orbital + 1)]
        #     # This means that active space must be continuous, orbital 1,2,3 not 1,3,4.
        #     # Should be okay if using no symmetry

        # Deletes Previous Orbital Density Grids
        Orbital_Grid = Molcas_Directory + '/grid' + str(orbital)
        if path.exists(Orbital_Grid):
            remove(Orbital_Grid)

        # Reads and Writes Densities Orbital by Orbital
        with open(gridfilepath + "/" + output_filename + '.grid', 'r') as fin:
            with open(Molcas_Directory + '/grid' + str(orbital), 'a') as fout:
                for j in reset_lines:
                    copy = False
                    for num, line in enumerate(fin):
                        if num > j:
                            if orbital < Highest_orbital:
                                if all(match in line for match in string):
                                    copy = True
                                elif all(match in line for match in stringstop):
                                    copy = False
                                elif copy:
                                    fout.write(line)
                            elif orbital == Highest_orbital:
                                if all(match in line for match in string):
                                    copy = True
                                elif all(match in line for match in stringend):
                                    copy = False
                                elif copy:
                                    fout.write(line)
        print("grid" + str(orbital))
    print("Extracted Density Grid")


# Gets Dipoles found inside of <Project>.log "RASSI" (2,0)
def CopyXYZDipoles(filein, fileout, linestart, linestop, numberofstates):
    linestop2 = "PROPERTY: MLTPL  2"
    with open(filein, 'r') as fin:
        with open(fileout, 'w') as fout:
            copy = False
            for line in fin:
                if linestart in line:
                    copy = True
                    next(fin)
                    next(fin)
                    next(fin)
                elif "STATE" in line:
                    next(fin)
                elif linestop in line:
                    copy = False
                elif linestop2 in line:
                    copy = False
                elif copy:
                    fout.write(line)
    range_states = range(0, numberofstates, 1)
    with open(fileout + "temp", 'w') as fout:
        fout.write("\n")
        for states in range_states:
            states += 1
            string = " " + str(states) + " "
            value = rfunc.GetDipoleValuesAsArray(fileout, string, "      ")
            fout.write(' '.join([str(f) for f in value]) + "\n")
    rfunc.ExecuteNoWrite("cp " + fileout + "temp " + fileout)
    rfunc.ExecuteNoWrite("rm " + fileout + "temp")


def GetDipolesFromLogFile(Project_Name, outdir, numberofstates, workingdirectory):
    componentlist = list("123")

    logdirectory = rfunc.Find(Project_Name + ".log", ".", workingdirectory)

    for component in componentlist:
        component1 = str(int(float(component) + 1))
        if component == "1":
            dipolename = "X_DIPOLE"
        elif component == "2":
            dipolename = "Y_DIPOLE"
        elif component == "3":
            dipolename = "Z_DIPOLE"
        CopyXYZDipoles(logdirectory + "/" + Project_Name + ".log", outdir + "/" + dipolename,
                       "PROPERTY: MLTPL  1   COMPONENT:   " + component,
                       "PROPERTY: MLTPL  1   COMPONENT:   " + component1, numberofstates)
        print("Extracted " + dipolename)


# >Makes a color matrix refering to Dipole between each state(1)
def Make_MU_HeatMap(directory):
    rows = rfunc.file_len(directory + "/X_DIPOLE")
    cols = rows
    direction_list = list("XYZ")
    for direction in direction_list:
        with open(directory + "/" + direction + '_DIPOLE', 'r') as fin:
            data = []
            next(fin)
            for i in range(1, rows):
                data.append(list(map(float, fin.readline().split()[:cols])))
        df = pd.DataFrame(data)
        # mask = np.zeros_like(data)
        # mask[np.triu_indices_from(mask)] = True # cuts out upper triangular side of the matrix
        f, ax = plt.subplots(figsize=(11, 9))
        cmap = sns.color_palette("coolwarm", as_cmap=True)
        sns.heatmap(
            df,  # The data to plot
            # mask=mask,  # Mask some cells
            cmap=cmap,  # What colors to plot the heatmap as
            annot=True,  # Should the values be plotted in the cells?
            vmax=1,  # The maximum value of the legend. All higher vals will be same color
            vmin=-1,  # The minimum value of the legend. All lower vals will be same color
            center=0,  # The center value of the legend. With divergent cmap, where white is
            square=True,  # Force cells to be square
            linewidths=.5,  # Width of lines that divide cells
            cbar_kws={"shrink": .5}  # Extra kwargs for the legend; in this case, shrink by 50%
        )
        f.savefig(directory + "/" + direction + '_DIPOLE_heatmap')

#
# from io import StringIO
# from pathlib import Path
#
# from matplotlib import pyplot as plt
# import pandas as pd
# import seaborn as sns
# from typing import Iterable
#
# from matplotlib.axes import Axes
# from matplotlib.colors import Colormap
# from matplotlib.figure import Figure
#
#
# class Component:
#     NAMES = 'XYZ'
#     COLORMAP: Colormap = sns.color_palette('coolwarm', as_cmap=True)
#
#     def __init__(self, section: str):
#         *_, component, body = section.split(maxsplit=3)
#         self.index = int(component)
#         self.name = self.NAMES[self.index - 1]
#         self.df = self._parse(body)
#
#     @classmethod
#     def _parse(cls, body: str) -> pd.DataFrame:
#         colsets = [
#             cls._str_to_df(colset)
#             for colset in cls._colset_strings(body)
#         ]
#         without_state = (
#             df.drop('STATE', 1)
#             for df in colsets
#         )
#
#         return pd.concat(
#             (colsets[0].STATE, *without_state),
#             axis=1,
#         )
#
#     @staticmethod
#     def _colset_strings(body: str) -> Iterable[str]:
#         first_col = 'STATE'
#         start = body.find(first_col)
#         while start != -1:
#             end = body.find(first_col, start + len(first_col))
#             yield body[start:end]
#             start = end
#
#     @staticmethod
#     def _str_to_df(body: str) -> pd.DataFrame:
#         with StringIO(body) as f:
#             return pd.read_csv(f, sep=r'\s+')
#
#     def __str__(self):
#         return self.name
#
#     @property
#     def without_state(self) -> pd.DataFrame:
#         return self.df.drop('STATE', axis=1)
#
#     def write_csv(self, path: Path, bare: bool = True):
#         if bare:
#             self.without_state.to_csv(path, index=False, header=False)
#         else:
#             self.df.to_csv(path, index=False)
#
#     def mu_heatmap(self, ax: Axes) -> Axes:
#         ax.set_title(self.name)
#
#         return sns.heatmap(
#             self.without_state,  # The data to plot
#             ax=ax,               # Axes to draw on
#             cmap=self.COLORMAP,  # What colors to plot the heatmap as
#             annot=True,          # Should the values be plotted in the cells?
#             # The maximum value of the legend. All higher vals will be same color
#             vmax=1,
#             # The minimum value of the legend. All lower vals will be same color
#             vmin=-1,
#             # The center value of the legend. With divergent cmap, where white is
#             center=0,
#             square=True,    # Force cells to be square
#             linewidths=.5,  # Width of lines that divide cells
#             # Extra kwargs for the legend; in this case, shrink by 50%
#             cbar_kws={'shrink': .5},
#         )
#
#
# class DipolesProject:
#     def __init__(self, log_path: str):
#         self.components = tuple(self._parse(log_path))
#
#     @staticmethod
#     def _parse(log_path: str):
#         with open(log_path) as f:
#             body = f.read()
#
#         start_i = body.find('PROPERTY: MLTPL  1 ')
#         end_i = body.find('PROPERTY: MLTPL  2 ', start_i)
#         body = body[start_i: end_i]
#         sections = body.split('PROPERTY: MLTPL')[1:]
#
#         for section in sections:
#             yield Component(section)
#
#     def write_csvs(
#         self,
#         directory: str = '.',
#         suffix: str = '_DIPOLE.csv',
#         bare: bool = True,
#     ):
#         directory = Path(directory)
#         for component in self.components:
#             component.write_csv(directory / (component.name + suffix), bare)
#
#     def get_mu_heatmaps(self) -> Figure:
#         fig, ((z, x), (y, unused)) = plt.subplots(
#             2, 2, sharex='col', sharey='row',
#             figsize=(10, 10),
#         )
#         unused.remove()
#
#         for ax, component in zip((x, y, z), self.components):
#             component.mu_heatmap(ax)
#
#         return fig
#
#
# def main():
#     project = DipolesProject('CO2.log')
#     project.write_csvs()
#     fig = project.get_mu_heatmaps()
#     fig.savefig('dipole-heatmap.png')
#
#
# if __name__ == '__main__':
#     main()
