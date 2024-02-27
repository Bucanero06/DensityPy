import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import pyvista as pv


def extract_simulation_time_from_filename(filename):
    """
    Extracts the simulation time from a filename with format 'ChDen<time>.csv.gz'.

    Args:
        filename (str): The filename to extract the time from.

    Returns:
        float: The extracted simulation time, or None if the pattern is not found.
    """
    match = re.search(r'ChDen([\-0-9\.]+)\.csv.gz', filename)
    if match:
        return float(match.group(1))
    else:
        return None


def convert_density_csv_to_vtk(csv_file: str, csv_directory: str, vtk_directory: str,
                               x_header: str, y_header: str, z_header: str, density_header: str):
    """
    Converts a CSV file containing density data into a VTK structured grid file.

    Args:
        csv_file (str): The name of the CSV file to convert.
        csv_directory (str): The directory containing the CSV file.
        vtk_directory (str): The directory where the VTK file will be saved.
        x_header (str): The header name for the x-coordinate in the CSV file.
        y_header (str): The header name for the y-coordinate in the CSV file.
        z_header (str): The header name for the z-coordinate in the CSV file.
        density_header (str): The header name for the density data in the CSV file.

    Returns:
        tuple: A tuple containing the simulation time and the name of the generated VTK file.
    """
    simulation_time = extract_simulation_time_from_filename(csv_file)
    df = pd.read_csv(os.path.join(csv_directory, csv_file))
    grid = pv.StructuredGrid()
    grid.dimensions = [len(df[x_header].unique()), len(df[y_header].unique()), len(df[z_header].unique())]
    grid.points = df[[x_header, y_header, z_header]].values
    grid.point_data[density_header] = df[density_header].values
    vtk_filename = os.path.splitext(csv_file)[0] + '.vts'
    grid.save(os.path.join(vtk_directory, vtk_filename))
    return simulation_time, vtk_filename


def convert_all_csvgz_to_vtk(csv_directory, vtk_directory, x_header, y_header, z_header, density_header,
                             max_workers=4):
    """
    Converts all CSV.GZ files in a directory to VTK files in parallel using a thread pool.

    Args:
        csv_directory (str): The directory containing the CSV.GZ files.
        vtk_directory (str): The directory where the VTK files will be saved.
        x_header (str): The header name for the x-coordinate in the CSV files.
        y_header (str): The header name for the y-coordinate in the CSV files.
        z_header (str): The header name for the z-coordinate in the CSV files.
        density_header (str): The header name for the density data in the CSV files.
        max_workers (int): The maximum number of threads to use for parallel conversion.
    """
    if not os.path.exists(vtk_directory):
        os.makedirs(vtk_directory)

    csv_files = [f for f in os.listdir(csv_directory) if f.endswith('.csv.gz')]
    csv_files.sort(key=extract_simulation_time_from_filename)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_csv = {executor.submit(
            convert_density_csv_to_vtk, csv_file, csv_directory, vtk_directory,
            x_header, y_header, z_header, density_header): csv_file for
                         csv_file in csv_files}
        pvd_entries = []

        for future in as_completed(future_to_csv):
            simulation_time, vtk_filename = future.result()
            pvd_entries.append((simulation_time, vtk_filename))
            print(f'Converted {future_to_csv[future]} to {vtk_filename}')

    # Create the PVD file after all conversions are complete
    pvd_filename = os.path.join(vtk_directory, 'timeseries.pvd')
    with open(pvd_filename, 'w') as pvd_file:
        pvd_file.write('<?xml version="1.0"?>\n')
        pvd_file.write('<VTKFile type="Collection" version="0.1">\n')
        pvd_file.write('  <Collection>\n')

        for simulation_time, vtk_filename in sorted(pvd_entries, key=lambda x: x[0]):
            vtk_file_path = os.path.join(vtk_directory, vtk_filename)
            pvd_file.write(f'    <DataSet timestep="{simulation_time}" file="{vtk_file_path}" />\n')

        pvd_file.write('  </Collection>\n')
        pvd_file.write('</VTKFile>\n')


if __name__ == '__main__':
    BASE_PATH = "/home/ruben/PycharmProjects/DensityPy/Studies/ExampleStudy/cm_plot_sim/ChargeDensity"
    csv_directory = f'{BASE_PATH}/ChDenSimPP150'
    vtk_directory = f'{BASE_PATH}/VTK'
    convert_all_csvgz_to_vtk(csv_directory, vtk_directory,
                             x_header='x', y_header='y', z_header='z', density_header='ChargeDensity',
                             max_workers=30)
