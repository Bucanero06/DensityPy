#! /usr/bin/env python3.10
# >import molcasscripts as rmolcs
# TODO: use other available solutions to these problems when refactoring, preferably from the standard OpenMolcas
#  library
import os
import sys
from os import path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from densitypy.project_utils.def_functions import find, \
    get_dipole_values_as_array, file_len, execute_command, print_molcas_log_errors, copy_file_to, \
    delete_files_or_directories
from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])


# Creates Helpful input file for Molcas with comments to help user begin (0)
def create_help_input_file():
    # >Creates helpful input file guide for pymolcas
    with open('molcas_input_help.input', 'w') as fout:
        fout.write("&GATEWAY"
                   " \n Title = NMA //  _name_of_project_"
                   " \n coord = NMA.xyz //  _geometry_file_"
                   " \n basis = 6-31+g_st__st_.1.molcas //  _name_of_basis_"
                   " \n group = c1 //this prevents use of symmetry"
                   "\n&SEWARD\n&SCF"
                   "\n&RASSCF"
                   "\n Ras1   = 0 //  _number_of_RAS1_orbitals_"
                   "\n Ras2   = 4 //  _number_of_RAS2_orbitals_"
                   "\n RAS3   = 0 //  _number_of_RAS3_orbitals_"
                   "\n Inactive = 17 //  _number_of_Inactive_orbitals_"
                   "\n nactel   = 6 0 0 "
                   "//# of active electrons, allowed holes, allowed excitations.  _total_ _RAS1_ _RAS3_ "
                   "\n TDM       //writes DM and TDM to RASSCF.h5"
                   "\n CIRoot\n 4 4 1 "
                   "//e.g. 10 10 1 will do a state average calculation of the ground state "
                   "and the lower 4 excited states.  _highest_root_needed_ _dimention_of_CI_matrix_ _weight_"
                   "\n&RASSI"
                   "\nNr of JobIph"
                   "\n 1 ALL "
                   "    //e.g. # of job files, # of states in that file ==>> 1 10 or 1 ALL "
                   "\nMEES //writes one-electron properties"
                   )

        logger.info("molcas_input_help.input created")



# >Makes points to be used for pymolcas input (0)
def make_grid_coordinates(molcas_directory, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax):
    # >Makes points to be used for pymolcas input
    n_points = Nx * Ny * Nz
    logger.info("Number of Points = " + str(n_points))
    with open(molcas_directory + '/gridcoord', 'w') as fout:
        for i in range(0, Nx):
            x = xmin + (xmax - xmin) / (Nx - 1) * i
            for j in range(0, Ny):
                y = ymin + (ymax - ymin) / (Ny - 1) * j
                for k in range(0, Ny):
                    z = zmin + (zmax - zmin) / (Nz - 1) * k
                    fout.write("{} {} {}\n".format(x, y, z))


##################################################################################
def read_xyz(xyz_file):
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

    logger.info("filename:         %s" % xyz_file)
    logger.info("title:            %s" % title)
    logger.info("number of atoms:  %d" % n_atoms)

    return atomic_coordinates


def distance(avec, bvec):
    c = np.square(np.subtract(avec, bvec))
    cvec = np.sqrt(np.sum(c))
    return cvec


def make_new_grid(atomic_coordinates, step_size, boundary, limitedgrid, np=None):
    # > Define coordinates matrix
    gridcoord = [] * 3

    maxi, mini = find_max_and_min(atomic_coordinates)
    # # > Find Max Boundries of Grid
    xmax = maxi[0] + boundary
    ymax = maxi[1] + boundary
    zmax = maxi[2] + boundary
    xmin = mini[0] - boundary
    ymin = mini[1] - boundary
    zmin = mini[2] - boundary

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
                logger.info('\r' + str(percentage_done) + ' %', end='  ')

                if limitedgrid == True:
                    for icoord in atomic_coordinates:
                        avec_new = []
                        avec_new.append(float(str(icoord[0])[1:-1]))
                        avec_new.append(float(str(icoord[1])[1:-1]))
                        avec_new.append(float(str(icoord[2])[1:-1]))

                        if abs(distance(avec_new, r_point)) <= boundary:
                            gridcoord.append([X, Y, Z])
                            break
                else:
                    gridcoord.append([X, Y, Z])
    logger.info("Done Making Grid")
    logger.info("\n")
    return gridcoord, initial_Npoints


def find_max_and_min(coords):
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


def make_better_grid(molcas_directory, geometry, step_size, Boundary, limitedgrid, ):
    logger.info("Making Grid")
    atomic_coordinates = read_xyz(geometry)
    # > Make_Grid(maxi, mini, Npoints)
    gridcoord, initial_Npoints = make_new_grid(atomic_coordinates, step_size, Boundary, limitedgrid)
    n_points = len(gridcoord)
    logger.info("Number of Points = " + str(n_points))
    # > Writes to file
    logger.info("Writing Grid to File")
    with open(molcas_directory + '/gridcoord', 'w') as fout:
        for i in range(0, n_points):
            fout.write(
                "{} {} {}\n".format(round(gridcoord[i][0], 2), round(gridcoord[i][1], 2), round(gridcoord[i][2], 2)))
    return n_points


#################################################################################################
# Creates input from users own input file for pymolcas(0)
def copy_input_file_to_edit(pymolcas_input, project_name, molcas_directory):
    # Creates input from users own input file for pymolcas
    string_list = ['RASSCF', 'RASSI', 'TDM']
    true_values = []

    with open(pymolcas_input, 'r') as fin:
        #
        with open(molcas_directory + '/' + project_name + '.input', 'w') as fout:
            logger.info('Creating ' + molcas_directory + '/' + project_name + '.input from ' + pymolcas_input)
            for line in fin:
                if line.startswith('//') is False:
                    fout.write(line)
                    for string in string_list:
                        if string.upper() in line.upper():
                            true_values.append(string.upper())
    #
    return true_values


# Adds GRID_IT section to copied input file
def add_grid_it_to_manual_input_file(pymolcas_input, project_name, molcas_directory, orbital_list, n_points):
    # Creates input from users own input file for pymolcas
    Lowest_orbital = min(orbital_list)
    Highest_orbital = max(orbital_list)
    with open(molcas_directory + '/gridcoord', 'r') as fin:
        # with open(project_name + '.input', 'a') as fout:
        with open(f'{molcas_directory}/{project_name}.input', 'a') as fout:
            logger.info('Adding Grid_it to ' + project_name + '.input')
            fout.write('\n&GRIDIT '
                       '\n select '
                       '\n 1:' + str(Lowest_orbital) + '-' + str(Highest_orbital) +
                       '\n NOLUSCUS '
                       '\n GRID'
                       '\n ')
            fout.write(str(n_points))
            fout.write("\n" + fin.read())


# >Calls pymolcas and writes to scrach folder and log file(1)
def call_open_molcas(project_name, molcas_directory):
    # >Calls pymolcas and writes to scrach folder and log file
    logger.info("Running OpenMolcas")
    #
    current_directory = os.getcwd()
    os.chdir(molcas_directory)

    # rfunc.ExecutePymolcasWithErrorPrint("pymolcas " + project_name + ".input -f", project_name)
    # execute_pymolcas_with_error_print("pymolcas " + project_name + ".input -f", project_name)
    # execute_pymolcas_with_error_print(f'pymolcas {project_name}.input -f', project_name)
    try:
        execute_command(f'pymolcas {project_name}.input -f', _logger=logger)
    except Exception as e:
        logger.error('Error in call_open_molcas')
        print_molcas_log_errors(project_name + ".log", "Timing")
        raise e

    logger.info("Done Running OpenMolcas")
    os.chdir(current_directory)


# >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors(1,1,2)
def fetch_fromh5_file(project_name, molcas_directory, pathtofile, filename):
    # >Extracts Information from H5FILEs
    execute_command("h5dump -o " + molcas_directory + "/" + filename +
                    " -y -w 0 -d " + filename + " " + pathtofile + "/" + project_name + ".rasscf.h5", _logger=logger)
    logger.info("Extracted " + filename)


def fetch_fromh5_filewidth1(project_name, molcas_directory, pathtofile, filename):
    # >Extracts Information from H5FILEs
    execute_command("h5dump -o " + molcas_directory + "/" + filename +
                    " -y -w 1 -d " + filename + " " + pathtofile + "/" + project_name + ".rasscf.h5", _logger=logger)
    logger.info("Extracted " + filename)


def load_fromh5_file(project_name, molcas_directory, workingdirectory, true_values, density_matrix,
                     transition_density_matrix,
                     root_energies, ao_multiple_x, ao_multiple_y, ao_multiple_z,
                     mo_energies, mo_vectors, justh5):
    # >Find <ProjectName>.rasscf.h5
    rasscf_h5_filepath = find(project_name + ".rasscf.h5", ".", workingdirectory)

    # >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors
    fetch_fromh5_file(project_name, molcas_directory, rasscf_h5_filepath, density_matrix)
    if "TDM" in true_values or justh5:
        fetch_fromh5_file(project_name, molcas_directory, rasscf_h5_filepath, transition_density_matrix)
    else:
        logger.info(f'Found {project_name}".rasscf.h5" but TDM keyword was not found in Molcas input file, '
                    f'thus no Transition Density Matrix printed')
    fetch_fromh5_filewidth1(project_name, molcas_directory, rasscf_h5_filepath, root_energies)
    fetch_fromh5_file(project_name, molcas_directory, rasscf_h5_filepath, ao_multiple_x)
    fetch_fromh5_file(project_name, molcas_directory, rasscf_h5_filepath, ao_multiple_y)
    fetch_fromh5_file(project_name, molcas_directory, rasscf_h5_filepath, ao_multiple_z)
    fetch_fromh5_filewidth1(project_name, molcas_directory, rasscf_h5_filepath, mo_energies)
    fetch_fromh5_filewidth1(project_name, molcas_directory, rasscf_h5_filepath, mo_vectors)


# >Extracts Density for the grids of each Orbital found oin <Project>.grid (0)
def extract_grid_density(orbital_list, output_filename, directorypath, molcas_directory):
    Highest_orbital = max(orbital_list)
    stringend = "Title=    0 "
    reset_lines = []
    reset_lines.append(0)
    gridfilepath = find(output_filename + '.grid', ".", directorypath)
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
        orbital_grid = molcas_directory + '/grid' + str(orbital)
        if path.exists(orbital_grid):
            delete_files_or_directories()

        # Reads and Writes Densities Orbital by Orbital
        with open(gridfilepath + "/" + output_filename + '.grid', 'r') as fin:
            with open(molcas_directory + '/grid' + str(orbital), 'a') as fout:
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
        # logger.info("grid" + str(orbital))
        logger.info(f'Extracted Density Grid for Orbital {orbital}')
    logger.info("Successfully Extracted Density Grid")  # todo perhaps... add a check to see if all grids were extracted


# Gets Dipoles found inside of <Project>.log "RASSI" (2,0)
def copy_xyz_dipoles(filein, fileout, linestart, linestop, numberofstates):
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
            value = get_dipole_values_as_array(fileout, string, "      ")
            fout.write(' '.join([str(f) for f in value]) + "\n")

    # execute_command("cp " + fileout + "temp " + fileout)
    # execute_command("rm " + fileout + "temp")

    copy_file_to(f'{fileout}temp', fileout)
    delete_files_or_directories(f'{fileout}temp')


def get_dipoles_from_log_file(project_name, output_dir, number_of_states, working_directory):
    componentlist = list("123")

    logdirectory = find(project_name + ".log", ".", working_directory)

    for component in componentlist:
        component1 = str(int(float(component) + 1))
        if component == "1":
            dipolename = "X_DIPOLE"
        elif component == "2":
            dipolename = "Y_DIPOLE"
        elif component == "3":
            dipolename = "Z_DIPOLE"
        copy_xyz_dipoles(logdirectory + "/" + project_name + ".log", output_dir + "/" + dipolename,
                         "PROPERTY: MLTPL  1   COMPONENT:   " + component,
                         "PROPERTY: MLTPL  1   COMPONENT:   " + component1, number_of_states)
        logger.info("Extracted " + dipolename)


# >Makes a color matrix refering to Dipole between each state(1)
def make_mu_heat_map(directory):
    rows = file_len(directory + "/X_DIPOLE")
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
