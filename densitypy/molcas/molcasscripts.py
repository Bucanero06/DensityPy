#! /usr/bin/env python3.10
# >import molcasscripts as rmolcs
# TODO: use other available solutions to these problems when refactoring, preferably from the standard OpenMolcas
#  library
import csv
import re

import numpy as np
import pandas as pd

from densitypy.project_utils.def_functions import find, \
    execute_command, print_molcas_log_errors, change_directory, \
    extract_datasets_from_h5_to_csv
from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])
MOLECULAR_ORBITAL_DATA = {  # todo not used yet
    'meta_info': {
        # Orbital sym=',i2,' index=',i5,' Energ=',F12.4,' occ=',F4.2,' type=',a1
        'sys': None,
        'energy': None,
        'occ': None,
        'type': None},
    'density': []
}


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
    with open(molcas_directory + '/gridcoord.csv', 'w') as fout:
        fout.write("{},{},{}\n".format("X", "Y", "Z"))
        for i in range(0, Nx):
            x = xmin + (xmax - xmin) / (Nx - 1) * i
            for j in range(0, Ny):
                y = ymin + (ymax - ymin) / (Ny - 1) * j
                for k in range(0, Ny):
                    z = zmin + (zmax - zmin) / (Nz - 1) * k
                    fout.write("{},{},{}\n".format(x, y, z))
    return n_points


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


def make_better_grid(molcas_directory, geometry, step_size, Boundary, limitedgrid, ):
    logger.info("Making Grid")
    atomic_coordinates = read_xyz(geometry)
    # > Make_Grid(maxi, mini, Npoints)
    gridcoord, initial_Npoints = make_new_grid(atomic_coordinates, step_size, Boundary, limitedgrid)
    n_points = len(gridcoord)
    logger.info("Number of Points = " + str(n_points))
    # > Writes to file
    logger.info("Writing Grid to File")
    with open(molcas_directory + '/gridcoord.csv', 'w') as fout:
        fout.write("{},{},{}\n".format("X", "Y", "Z"))
        for i in range(0, n_points):
            fout.write(
                "{},{},{}\n".format(round(gridcoord[i][0], 2), round(gridcoord[i][1], 2), round(gridcoord[i][2], 2)))
    return n_points


def add_grid_it_to_manual_input_file(pymolcas_input, project_name, molcas_directory, orbital_list, n_points):
    Lowest_orbital = min(orbital_list)
    Highest_orbital = max(orbital_list)

    # Construct the output file path
    output_file_path = f'{molcas_directory}/{project_name}.input'

    with open(molcas_directory + '/gridcoord.csv', 'r') as fin, open(output_file_path, 'a') as fout:
        # Using csv.reader to read the CSV file
        csv_reader = csv.reader(fin)

        # Skip the first row (header)
        next(csv_reader)

        # Log the information
        logger.info(f'Adding Grid_it to {project_name}.input')

        # Writing the initial part of the output
        fout.write('\n&GRIDIT '
                   '\n select '
                   f'\n 1:{Lowest_orbital}-{Highest_orbital}'
                   '\n NOLUSCUS '
                   '\n GRID'
                   '\n ')
        fout.write(str(n_points))
        fout.write('\n')

        # Writing each row of the CSV to the output file
        for row in csv_reader:
            fout.write(','.join(row) + '\n')


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


#################################################################################################
# Creates input from users own input file for pymolcas(0)
def copy_and_parse_molcas_input_file_to_edit(pymolcas_input, project_name, molcas_directory):
    """
    Uses the user's input file to create a molcas input file by copying the user's input file and then
    utilizing it to add other necessary sections to the input file.
    """

    KEYWORDS_NEEDED = {
        "TDM": 'Transition Density Matrix',
        "RASSCF": 'Density Matrix',
        "RASSI": 'Density Matrix'
    }
    keywords_needed_found = []

    # Create input from user's own input file for pymolcas
    with open(pymolcas_input, 'r') as fin, open(f'{molcas_directory}/{project_name}.input', 'w') as fout:
        logger.info(f'Creating {molcas_directory}/{project_name}.input from {pymolcas_input}')

        for line in fin:
            if not line.startswith('//'):
                fout.write(line)
                for string in KEYWORDS_NEEDED.keys():
                    if string.upper() in line.upper():
                        keywords_needed_found.append(string.upper())

    # Validate input file
    # TODO: Implement the validation logic
    for keyword, matrix in KEYWORDS_NEEDED.items():
        if keyword not in keywords_needed_found:
            logger.warning(
                f'{keyword} keyword was not found in Molcas input file. As a result, no {matrix} will be created and '
                f'ChargeMigration modules will not be able to run.'
            )
    return keywords_needed_found


# >Calls pymolcas and writes to scrach folder and log file(1)
def call_open_molcas(project_name, molcas_directory):
    # >Calls pymolcas and writes to scratch folder and log file
    logger.info("Running OpenMolcas")

    # Use the context manager to temporarily change the working directory
    with change_directory(molcas_directory):
        try:
            execute_command(f'pymolcas {project_name}.input -f', _logger=logger)
        except Exception as e:
            logger.error('Error in call_open_molcas')
            print_molcas_log_errors(project_name + ".log", "Timing")
            raise e

    logger.info("Done Running OpenMolcas")


# >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors(1,1,2)


# >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors(1,1,2)
def load_project_rasscf_h5_file(project_name, molcas_output_directory, workingdirectory):
    # >Extract Density Matrix, Transition Density Matrix, ENERGIES, AO_Dipoles, MO_energies, MO_vectors
    DATASETS_TO_EXTRACT = ["DENSITY_MATRIX",
                           "TRANSITION_DENSITY_MATRIX", "ROOT_ENERGIES", "AO_MLTPL_X", "AO_MLTPL_Y",
                           "AO_MLTPL_Z", "MO_ENERGIES", "MO_VECTORS"]

    # The keys are the names of the datasets to extract and the values are used for the output files
    dataset_mappings = {dataset: f'{molcas_output_directory}/{dataset}' for dataset in DATASETS_TO_EXTRACT}

    # >Find <ProjectName>.rasscf.h5
    rasscf_h5_dirpath = find(f'{project_name}.rasscf.h5', ".", workingdirectory)
    extract_datasets_from_h5_to_csv(
        h5_filepath=f'{rasscf_h5_dirpath}/{project_name}.rasscf.h5',
        # >Define the dataset mappings using DATASETS_TO_EXTRACT as keys and molcas_output_directory/dataset names as
        # values
        dataset_mapping=dataset_mappings
    )


# >Extracts Density for the grids of each Orbital found oin <Project>.grid (0)
def parse_project_grid_file(output_filename, directorypath):
    gridfilepath = find(output_filename + '.grid', ".", directorypath)

    # Open the file
    with open(gridfilepath + "/" + output_filename + '.grid', 'r') as file:
        lines = file.readlines()

    data = {
        'Natom': None,
        'atoms': [],
        'meta_info': {
            'version': None,
            'N_of_MO': None,
            'N_of_Grids': None,
            'N_of_Points': None,
            'Block_Size': None,
            'N_Blocks': None,
            'Is_cutoff': None,
            'CutOff': None,
            'N_P': None,
            'N_INDEX': None,
        },
        'molecular_orbitals': {}
    }

    atom_flag = False
    store_the_index_of_the_line = None
    for line_index, line in enumerate(lines):  # from line 0 to line N_lines

        if line.startswith('Natom='):
            data['Natom'] = int(line.split('=')[-1].strip())
            atom_flag = True
            continue

        if atom_flag and len(data['atoms']) < data['Natom']:
            parts = line.split()
            atom_info = {
                'element': parts[0][:-1],
                'index': int(parts[0][-1]),
                'x': float(parts[1]),
                'y': float(parts[2]),
                'z': float(parts[3])
            }
            data['atoms'].append(atom_info)
            continue

        if line.startswith('VERSION='):
            data['meta_info']['version'] = line.split('=')[-1].strip()
            continue

        if line.startswith('N_of_MO='):
            data['meta_info']['N_of_MO'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('N_of_Grids='):
            data['meta_info']['N_of_Grids'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('N_of_Points='):
            data['meta_info']['N_of_Points'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('Block_Size='):
            data['meta_info']['Block_Size'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('N_Blocks='):
            data['meta_info']['N_Blocks'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('Is_cutoff='):
            data['meta_info']['Is_cutoff'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('CutOff='):
            data['meta_info']['CutOff'] = float(line.split('=')[-1].strip())
            continue

        if line.startswith('N_P='):
            data['meta_info']['N_P'] = int(line.split('=')[-1].strip())
            continue

        if line.startswith('N_INDEX='):
            data['meta_info']['N_INDEX'] = [int(i) for i in line.split('=')[-1].strip().split()]
            continue

        if line.startswith('GridName='):
            pattern = r"GridName=\s*(\d+)\s*(\d+)\s*([\-\d\.]+)\s*\(([\d\.]+)\)\s*(\d)"
            match = re.search(pattern, line)
            if match:
                # Orbital sym=',i2,' index=',i5,' Energ=',F12.4,' occ=',F4.2,' type=',a1
                orbital_index = int(match.group(2))
                orbital_info = {
                    'sym': int(match.group(1)),
                    'energy': float(match.group(3)),
                    'occ': float(match.group(4)),
                    'type': int(match.group(5))
                }
                data['molecular_orbitals'][orbital_index] = {
                    'meta_info': orbital_info,
                    'density': []
                }
                continue

        if line.startswith('Title='):  # break and continue lines from here
            store_the_index_of_the_line = line_index
            break

    assert store_the_index_of_the_line is not None, "Could not find Title=<MOj> in the grid file"
    density_containing_lines = lines[store_the_index_of_the_line:]
    actual_n_chunks = data['meta_info']['N_Blocks']
    orbitals_to_extract = list(data['molecular_orbitals'].keys()) * actual_n_chunks
    orbital_index = None

    # FIXME DEBUGGING
    # lets check the already outputed file f'{molcas_output_directory}/grid_density.csv' and compare each columun/orbital density
    # with the one extracted here
    # Read the file
    df = pd.read_csv(f'{directorypath}/grid_density.csv')

    print(df.head())
    for line_index, line in enumerate(density_containing_lines):  # use re for matching

        # e.g. Title=   18  or Title= 1 or Title=    194 etc...
        pattern = r'Title=\s*(\d+)'
        match = re.search(pattern, line)

        if match:

            orbital_index = int(match.group(1))


            if orbital_index not in orbitals_to_extract:
                assert orbital_index == 0
                if not orbitals_to_extract:
                    break
                continue
            # data['molecular_orbitals'][orbital_index]['density'].append(float(0))  # fixme appending is expensive

            assert (orbital_index in orbitals_to_extract), f"Orbital {orbital_index} not found in orbitals_to_extract"
            # pop the orbital from the list
            orbitals_to_extract.pop(orbitals_to_extract.index(orbital_index))
            continue

        if not orbital_index == 0:
            data['molecular_orbitals'][orbital_index]['density'].append(float(line))  # fixme appending is expensive

    # print(json.dumps(data, indent=4))
    logger.info(
        "Successfully Extracted Density Grids for each of the orbitals")  # todo perhaps... add a check to see if all grids were extracted

    return data


def write_grid_density_file(grid_file_data_dict, molcas_output_directory):  # todo pythonic data model
    # A Dataframe [N_Orbitals x N_Points]
    grids_meta_info = grid_file_data_dict['meta_info']
    molecular_orbitals_data = grid_file_data_dict['molecular_orbitals']

    # Prepare dataframe
    df = pd.DataFrame(columns=molecular_orbitals_data.keys(),
                      # index=range(grids_meta_info['N_of_Points'])
                      )
    for orbital_index, orbital_data in molecular_orbitals_data.items():
        df[orbital_index] = orbital_data['density']

    # Make sure total number of indexes is the same as the number of points
    # assert len(df) == grids_meta_info['N_of_Points']

    # Write to file
    df.to_csv(f'{molcas_output_directory}/grid_density.csv', index=False)
