import os
import os.path
from multiprocessing import Pool

import imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from matplotlib import colors
from scipy.ndimage import gaussian_filter

from densitypy.molcas.molcasscripts import read_xyz
from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])

'''Pulses(FT) VS Time & Frequency'''


def plot_pulses(study_directory, experiment_directory, time_delays, min_time, max_time, plot_all=False):
    length_of_data_to_match = None
    # Plotting the pulse data
    if plot_all:
        time_delays_to_plot = time_delays
    else:
        zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
        one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
        two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
        max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))

        # Lets also do the same for the minimum time delay and the "thirds" splits as before (remember to get as close as possible to a number
        min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
        one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
        two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

        time_delays_to_plot = [zero_time_delay, one_third_max_time_delay, two_third_max_time_delay, max_time_delay,
                               min_time_delay, one_third_min_time_delay, two_third_min_time_delay]

    from matplotlib import pyplot as plt

    for time_delay in ['XUV'] + time_delays_to_plot:
        file_path = f'{study_directory}/{experiment_directory}/Pulses/pulsePP{time_delay}' if time_delay != 'XUV' else f'{study_directory}/{experiment_directory}/Pulses/pulseXUV'

        try:
            data = pd.read_csv(file_path, delim_whitespace=True, header=None,
                               names=[
                                   'Time',
                                   'Ax',
                                   'Ay',
                                   'Az',
                                   'Real_zw_m1',
                                   'Imag_zw_m1',
                                   'Real_zw_p0',
                                   'Imag_zw_p0',
                                   'Real_zw_p1',
                                   'Imag_zw_p1'
                               ])
        except FileNotFoundError:
            logger.error(f'File {file_path} not found')
            continue
        logger.info(f'Plotting {file_path}')

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'

        # logger.info all columns in pandas
        pd.set_option('display.max_columns', None)
        logger.debug(data.head())
        logger.debug(f'Length of data: {len(data)}')

        plt.figure(figsize=(12, 6))
        for col in data.columns[1:]:
            plt.plot(data['Time'], data[col], label=f'Column {col}')

        plt.xlabel('Time (arbitrary units)')
        plt.ylabel('Value')
        if time_delay != 'XUV':
            plt.title(f'Pulse Characteristics (pulsePP{time_delay})')
        else:
            plt.title(f'Pulse Characteristics (pulseXUV)')

        plt.xlim(min_time, max_time)
        plt.legend()
        plt.grid(True)
        # plt.show()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)
        plt.close('all')


def plot_ft_pulses(study_directory, experiment_directory, time_delays, plot_all=False):
    """
    Plotting the FT pulse data

    Expected file name formats:
        FTpulsePP{delay} and FTpulseXUV

    Columns:
        Freq: The frequency at which the Fourier transform is computed.
        FTAx: The x-component of the Fourier-transformed vector potential.
        FTAy: The y-component of the Fourier-transformed vector potential.
        FTAz: The z-component of the Fourier-transformed vector potential.
        FT_Aminus1_Real and FT_Aminus1_Imag: The real and imaginary parts of the Fourier transform for the mu = -1 component.
        FT_A0_Real and FT_A0_Imag: The real and imaginary parts of the Fourier transform for the mu = 0 component.
        FT_Aplus1_Real and FT_Aplus1_Imag: The real and imaginary parts of the Fourier transform for the mu = +1 component.

    :param study_directory:
    :param experiment_directory:
    :param time_delays:
    :param plot_all:
    :return:

    Output:
        FTpulsePP{delay}.png and FTpulseXUV.png
    """
    # Plotting the pulse data
    if plot_all:
        time_delays_to_plot = time_delays
    else:
        # Lets also do the same for the minimum time delay and the "thirds" splits as before (remember to get as close as possible to a number
        zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
        one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
        two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
        max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))
        min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
        one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
        two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

        time_delays_to_plot = [zero_time_delay, one_third_max_time_delay, two_third_max_time_delay, max_time_delay,
                               min_time_delay, one_third_min_time_delay, two_third_min_time_delay]

    from matplotlib import pyplot as plt
    for time_delay in ['XUV'] + time_delays_to_plot:
        file_path = f'{study_directory}/{experiment_directory}/Pulses/FTpulsePP{time_delay}' if time_delay != 'XUV' else f'{study_directory}/{experiment_directory}/Pulses/FTpulseXUV'
        #
        try:
            data = pd.read_csv(file_path, delim_whitespace=True, header=None,
                               names=[
                                   'Freq',
                                   'FTAx',
                                   'FTAy',
                                   'FTAz',
                                   'FT_Aminus1_Real',
                                   'FT_Aminus1_Imag',
                                   'FT_A0_Real',
                                   'FT_A0_Imag',
                                   'FT_Aplus1_Real',
                                   'FT_Aplus1_Imag'
                               ])
        except FileNotFoundError:
            logger.error(f'File {file_path} not found')
            continue
        logger.info(f'Plotting {file_path}')
        logger.debug(data.head())
        logger.debug(f'Length of data: {len(data)}')

        plt.figure(figsize=(12, 6))
        for col in data.columns[1:]:
            plt.plot(data['Freq'], data[col], label=f'Column {col}')
        # plt.plot(data['Time'], data['Ax'], label=f'Column Ax')
        plt.xlabel('Frequency (arbitrary units)')
        plt.ylabel('Value')
        if time_delay != 'XUV':
            plt.title(f'Pulse Characteristics (FTpulsePP{time_delay})')
        else:
            plt.title(f'Pulse Characteristics (FTpulseXUV)')

        plt.legend()
        plt.grid(True)
        # plt.show()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)
        plt.close('all')


'''Dipoles VS Time'''


def plot_dipoles_v_time(study_directory, experiment_directory, time_delays, min_time, max_time, plot_all=False):
    # Lets plot the Dipolar Reponse vs Time (t)
    # "itime","Time","DipoleX_Re","DipoleX_Im","DipoleY_Re","DipoleY_Im","DipoleZ_Re","DipoleZ_Im"
    """
    Plotting the dipole data vs time

    Expected file name formats:
        DipolePP{delay}.csv and DipoleXUV.csv

    Columns:
        itime: The time step at which the dipole is computed.
        Time: The time at which the dipole is computed.
        DipoleX_Re: The real part of the x-component of the dipole.
        DipoleX_Im: The imaginary part of the x-component of the dipole.
        DipoleY_Re: The real part of the y-component of the dipole.
        DipoleY_Im: The imaginary part of the y-component of the dipole.
        DipoleZ_Re: The real part of the z-component of the dipole.
        DipoleZ_Im: The imaginary part of the z-component of the dipole.

    :param study_directory:
    :param experiment_directory:
    :param time_delays:
    :param min_time:
    :param max_time:
    :param plot_all:
    :return:

    Output:
        DipolePP{delay}.png and DipoleXUV.png
    """

    length_of_data_to_match = None
    if plot_all:
        time_delays_to_plot = time_delays
    else:
        # Time delays to plot = [0, 1/3 max, 2/3 max, max] if not present then use the nearest. Do the same for min
        zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
        one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
        two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
        max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))
        min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
        one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
        two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

        time_delays_to_plot = [zero_time_delay, one_third_max_time_delay, two_third_max_time_delay, max_time_delay,
                               min_time_delay, one_third_min_time_delay, two_third_min_time_delay]

    import matplotlib.pyplot as plt

    for time_delay in ['XUV'] + time_delays_to_plot:
        # Force close all plots
        file_path = f'{study_directory}/{experiment_directory}/Dipole/DipolePP{time_delay}.csv' if time_delay != 'XUV' else f'{study_directory}/{experiment_directory}/Dipole/DipoleXUV.csv'

        try:
            data = pd.read_csv(file_path)
        except FileNotFoundError:
            logger.error(f'File {file_path} not found')
            continue
        logger.info(f'Plotting {file_path}')

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'
            # assert data.columns.tolist() in INDECES_COL_NAMES + FEATURES_COL_NAMES, f'Column names of {file_path} is not the same as the previous file'

        # Extracting relevant data for plotting
        time = data['Time']
        dipole_x_re = data['DipoleX_Re']
        dipole_x_im = data['DipoleX_Im']
        dipole_y_re = data['DipoleY_Re']
        dipole_y_im = data['DipoleY_Im']
        dipole_z_re = data['DipoleZ_Re']
        dipole_z_im = data['DipoleZ_Im']

        # Plotting the real and imaginary parts of the dipole components
        plt.figure(figsize=(16, 12))
        # Add the Main Top Title for the plot (file_path)
        plt.suptitle(file_path)
        plt.subplot(3, 2, 1)
        plt.plot(time, dipole_x_re)
        plt.xlim(min_time, max_time)
        plt.title('Dipole X Real')
        plt.xlabel('Time')
        plt.ylabel('Dipole X Re')

        plt.subplot(3, 2, 2)
        plt.plot(time, dipole_x_im)
        plt.xlim(min_time, max_time)
        plt.title('Dipole X Imaginary')
        plt.xlabel('Time')
        plt.ylabel('Dipole X Im')

        plt.subplot(3, 2, 3)
        plt.plot(time, dipole_y_re)
        plt.xlim(min_time, max_time)
        plt.title('Dipole Y Real')
        plt.xlabel('Time')
        plt.ylabel('Dipole Y Re')

        plt.subplot(3, 2, 4)
        plt.plot(time, dipole_y_im)
        plt.xlim(min_time, max_time)
        plt.title('Dipole Y Imaginary')
        plt.xlabel('Time')
        plt.ylabel('Dipole Y Im')

        plt.subplot(3, 2, 5)
        plt.plot(time, dipole_z_re)
        plt.xlim(min_time, max_time)
        plt.title('Dipole Z Real')
        plt.xlabel('Time')
        plt.ylabel('Dipole Z Re')

        plt.subplot(3, 2, 6)
        plt.plot(time, dipole_z_im)
        plt.xlim(min_time, max_time)
        plt.title('Dipole Z Imaginary')
        plt.xlabel('Time')
        plt.ylabel('Dipole Z Im')

        plt.tight_layout()
        # plt.show()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)
        plt.close('all')


def plot_atomic_dipoles_v_time(study_directory, experiment_directory, time_delays, min_time, max_time,
                               xyz_geometry_path, plot_all=False):
    # Read the XYZ file to get the number of atoms, the name of the atom and the atom names
    # (this is the order in which the atoms are listed in the CSV file)
    with open(xyz_geometry_path, 'r') as f:
        number_of_atoms = int(f.readline())
        name_of_molecule = f.readline().split()[0]
        atom_names = [f.readline().split()[0] for _ in range(number_of_atoms)]

    print(number_of_atoms)
    print(name_of_molecule)
    print(atom_names)

    ATOMIC_COL_NAMES = __generate_atomic_charge_column_names(xyz_geometry_path)

    print(f'{ATOMIC_COL_NAMES  = }')
    assert len(
        ATOMIC_COL_NAMES) == 3 * number_of_atoms, f'Number of columns in {xyz_geometry_path} does not match the number of atoms'

    if plot_all:
        time_delays_to_plot = time_delays
    else:
        # Time delays to plot = [0, 1/3 max, 2/3 max, max] if not present then use the nearest. Do the same for min
        zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
        one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
        two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
        max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))
        min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
        one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
        two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

        time_delays_to_plot = [zero_time_delay, one_third_max_time_delay, two_third_max_time_delay, max_time_delay,
                               min_time_delay, one_third_min_time_delay, two_third_min_time_delay]

    from matplotlib import pyplot as plt

    length_of_data_to_match = None
    for time_delay in ['XUV'] + time_delays_to_plot:
        file_path = f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargePP{time_delay}.csv' if time_delay != 'XUV' else f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeXUV.csv'
        try:
            data = pd.read_csv(file_path)
        except FileNotFoundError:
            logger.error(f'File {file_path} not found')
            continue
        logger.info(f'Plotting {file_path}')
        print(f'{data.head()}')
        print(f'{data.columns  = }')
        print(f'{data.columns.tolist()  = }')
        print(f'{ATOMIC_COL_NAMES  = }')
        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'
            # assert data.columns.tolist() in ATOMIC_COL_NAMES, f'Column names of {file_path} is not the same as the previous file'

        # Extracting relevant data for plotting
        time = data['Time']
        total_charge = data['TotalCharge']
        atom_charges = data.drop(columns=['itime', 'Time', 'TotalCharge'])

        # Subplots for xyz and total charge in one figure , each of the xyz will have all the atoms
        # Plotting charge thus is all real
        plt.figure(figsize=(16, 12))
        # Add the Main Top Title for the plot (file_path)
        plt.suptitle(file_path)

        plt.subplot(2, 2, 1)
        plt.plot(time, total_charge)
        plt.xlim(min_time, max_time)
        plt.title('Total Charge sum(Charge(iPol, iAtom, it))')
        plt.xlabel('Time (t)')
        plt.ylabel('Total Charge')

        plt.subplot(2, 2, 2)
        for col in atom_charges.columns:
            if 'ChargeX' in col:
                plt.plot(time, atom_charges[col], label=f'Column {col}')
        plt.xlim(min_time, max_time)
        plt.title('Atomic Charges X-Axis')
        plt.xlabel('Time (t)')
        plt.ylabel('Atomic Charges')
        plt.legend()

        plt.subplot(2, 2, 3)
        for col in atom_charges.columns:
            if 'ChargeY' in col:
                plt.plot(time, atom_charges[col], label=f'Column {col}')
        plt.xlim(min_time, max_time)
        plt.title('Atomic Charges Y-Axis')
        plt.xlabel('Time (t)')
        plt.ylabel('Atomic Charges')
        plt.legend()

        plt.subplot(2, 2, 4)
        for col in atom_charges.columns:
            if 'ChargeZ' in col:
                plt.plot(time, atom_charges[col], label=f'Column {col}')
        plt.xlim(min_time, max_time)
        plt.title('Atomic Charges Z-Axis')
        plt.xlabel('Time (t)')
        plt.ylabel('Atomic Charges')
        plt.legend()

        plt.tight_layout()
        # plt.show()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)
        plt.close('all')


def difference_between_dipole_and_atomic_charges_v_time(study_directory, experiment_directory, time_delays, min_time,
                                                        max_time,
                                                        xyz_geometry_path):
    import pandas as pd

    # Load the datasets
    atomic_charge_file = f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeXUV.csv'
    dipole_file = f'{study_directory}/{experiment_directory}/Dipole/DipoleXUV.csv'

    atomic_charge_data = pd.read_csv(atomic_charge_file)[::-1]
    dipole_data = pd.read_csv(dipole_file)

    ATOMIC_COL_NAMES = __generate_atomic_charge_column_names(xyz_geometry_path)
    xyz_file_content = read_xyz(xyz_geometry_path)

    # Reconstruct Dipole Response from Atomic Charges
    R_el_bc = pd.read_csv(f'{study_directory}/{experiment_directory}/R_el_bc.csv')
    print(f'{R_el_bc.keys()}')

    atomic_dipole_data = pd.DataFrame()
    for i, atom_name in enumerate(xyz_file_content['atoms'].keys()):
        for j, axis in enumerate(['X', 'Y', 'Z']):
            column_name = ATOMIC_COL_NAMES[3 * i + j]
            # atomic_dipole_data[column_name] = atomic_charge_data[column_name] * R_el_bc[axis + '_Position'][i]
            atomic_dipole_data[column_name] = atomic_charge_data[column_name] * xyz_file_content['atoms'][atom_name][j]

    # Comparing the sum of all dipoles with the total charge
    comparison_with_sum = pd.DataFrame({
        'Time': atomic_charge_data['Time'],
        'Total_Charge': atomic_charge_data['TotalCharge'],
        # dont know why atomic needs to be reversed, is seen across this file
        'Sum_Dipole_across_XYZ_ReIm': dipole_data.iloc[:, 2:].sum(axis=1),
        # 'Sum_Reconstructed_Dipole_across_XYZ': atomic_dipole_data.sum(axis=1)
    })

    print(comparison_with_sum.head())

    # Plotting the comparison
    import matplotlib.pyplot as plt

    plt.figure(figsize=(16, 12))
    plt.suptitle(f'{study_directory}/{experiment_directory} Dipole/Charge XUV Comparison')

    for i, col in enumerate(comparison_with_sum.columns[1:]):
        # Subplots for xyz and total charge in one figure , each of the xyz will have all the atoms
        # Plotting charge thus is all real
        plt.subplot(2, 2, i + 1)
        plt.plot(comparison_with_sum['Time'], comparison_with_sum[col], label=f'Column {col}')
        # subtitle for subplot
        plt.title(f'{col}')

    # Find the Error between the sum of the dipoles and the total charge
    comparison_with_sum['Mean Absolute Error (MAE)'] = comparison_with_sum.apply(
        lambda row: np.abs(-row['Total_Charge'] - row['Sum_Dipole_across_XYZ_ReIm']),
        axis=1)  # fixme this is wrong right? but error big when not negating the charge

    plt.subplot(2, 2, 4)
    plt.plot(comparison_with_sum['Time'], comparison_with_sum['Mean Absolute Error (MAE)'], label=f'MAE')
    plt.title(f'Error')
    plt.legend()

    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)
    # plt.show()
    output_file = f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeXUV_vs_DipoleXUV.png'
    plt.savefig(output_file)
    plt.close('all')


'''Dipoles(FT) VS Time & Frequency'''


def plot_ft_all_dipoles_v_time(study_directory, experiment_directory, dephasing_factor, relaxation_factor,
                               pump_settings, probe_settings, charge_migration_ft_settings, try_color_maps=False):
    """
    Plotting the FT dipole data vs time, vs frequency and vs time and frequency

    Currently is limited to a pump probe experiment with two pulses since, specifically using the central_time_2
    column as the tau time delay axis

    :param study_directory:
    :param experiment_directory:
    :param dephasing_factor:
    :param relaxation_factor:
    :param pump_settings:
    :param probe_settings:
    :param charge_migration_ft_settings:
    :return:
    """

    # "number_of_pulses", "central_time_1", "carrier_frequency_1", "fwhm_1", "carrier_envelope_phase_1",
    #   "intensity_1", "amplitude_1", "period_1", "central_time_...", "carrier_frequency_...", "fwhm_...",
    #   "carrier_envelope_phase_...", "intensity_...", "amplitude_...", "period_...", "iOmega", "OmegaVec", "FTDipoleX_Re",
    #   "FTDipoleX_Im", "FTDipoleY_Re", "FTDipoleY_Im", "FTDipoleZ_Re", "FTDipoleZ_Im"

    FILE_PATH = f'{os.path.join(study_directory, experiment_directory)}/Dipole/DipoleFT_ALL.csv'

    # Pump Settings
    type_of_pulse_pump = pump_settings['typeofpulse']
    start_time = pump_settings['starttime']
    pump_central_frequency = pump_settings['pumpcentralfrequency']
    pump_periods = pump_settings['pumpperiods']
    pump_phase = pump_settings['pumpphase']
    pump_intensity = pump_settings['pumpintensity']
    pump_polarization = pump_settings['pumppolarization']

    # Probe Settings
    type_of_pulse_probe = probe_settings['typeofpulse']
    time_delay_start = probe_settings['timedelaystart']
    time_delay_stop = probe_settings['timedelaystop']
    number_of_pp = probe_settings['numberofpp']
    time_Delay_weight_factor = probe_settings['timedelayweightfactor']
    probe_central_frequency = probe_settings['probecentralfrequency']
    probe_periods = probe_settings['probeperiods']
    probe_phase = probe_settings['probephase']
    probe_intensity = probe_settings['probeintensity']
    probe_polarization = probe_settings['probepolarization']

    # Charge Migration Settings
    ft_time_step = charge_migration_ft_settings['fttimestep']
    ft_width_step = charge_migration_ft_settings['ftwidthstep']

    # Lets plot the Spectra FT Dipolar Reponse vs Time (t)
    data = pd.read_csv(FILE_PATH)
    logger.debug(data.head())
    logger.debug(f'Length of data: {len(data)}')

    # Plotting the real and imaginary parts of the dipole components
    import matplotlib.pyplot as plt
    plt.figure(figsize=(16, 12))
    plt.suptitle(FILE_PATH)
    for i, col_end_name in enumerate(['X_Re', 'X_Im', 'Y_Re', 'Y_Im', 'Z_Re', 'Z_Im']):
        plt.subplot(3, 2, i + 1)
        plt.plot(data['central_time_2'], data[f'FTDipole{col_end_name}'])
        plt.title(f'FTDipole{col_end_name}')
        plt.xlabel('central_time_2')
        plt.ylabel(f'FTDipole{col_end_name}')

    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{study_directory}/{experiment_directory}/Dipole/DipoleFT_ALL_Components_vs_central_time_2.png')
    plt.close('all')

    """
    Average FT Dipole Magnitude: 
    X=central_time_2, Y=OmegaVec, 
    Z=((FTDipoleX_Re**2+FTDipoleX_Im**2+FTDipoleY_Re**2+FTDipoleY_Im**2+FTDipoleZ_Re**2+FTDipoleZ_Im**2)/3)**0.5
    """
    # Calculating the average of the squared magnitudes of the FT dipole components
    average_ft_dipole = np.sqrt((data['FTDipoleX_Re'] ** 2 + data['FTDipoleX_Im'] ** 2 +
                                 data['FTDipoleY_Re'] ** 2 + data['FTDipoleY_Im'] ** 2 +
                                 data['FTDipoleZ_Re'] ** 2 + data['FTDipoleZ_Im'] ** 2) / 3)

    # Creating an interactive plot using Plotly
    _plot_contour_map(x=data['central_time_2'],
                      y=data['OmegaVec'],
                      z=average_ft_dipole,
                      title=f'Contour Map of Average FT Dipole Magnitude\n'
                            f'DipoleFT_ALL.csv\n'
                            f'Pump Settings:  {pump_central_frequency}, {pump_periods}, {pump_phase}, {pump_intensity}, {pump_polarization}\n'
                            f'Probe Settings:  {probe_central_frequency}, {probe_periods}, {probe_phase}, {probe_intensity}, {probe_polarization}\n'
                            f'Dephasing Factor: {dephasing_factor}, Relaxation Factor: {relaxation_factor}'
                            f'FT Time Step: {ft_time_step}, FT Width Step: {ft_width_step}',
                      x_label='Central Time 2 - Probe (tau)',
                      y_label='OmegaVec',
                      output_file=f'{study_directory}/{experiment_directory}/Dipole/DipoleFT_ALL_Contour.html')

    # Plotting
    plt.figure(figsize=(15, 10))

    # Plotting the real parts
    plt.subplot(2, 2, 1)
    plt.plot(data['OmegaVec'], data['FTDipoleX_Re'], label='FTDipoleX_Re')
    plt.plot(data['OmegaVec'], data['FTDipoleY_Re'], label='FTDipoleY_Re')
    plt.plot(data['OmegaVec'], data['FTDipoleZ_Re'], label='FTDipoleZ_Re')
    plt.title('Real Components of FT Dipoles vs OmegaVec')
    plt.xlabel('OmegaVec')
    plt.ylabel('Real Part of FT Dipole')
    plt.legend()

    # Plotting the imaginary parts
    plt.subplot(2, 2, 2)
    plt.plot(data['OmegaVec'], data['FTDipoleX_Im'], label='FTDipoleX_Im')
    plt.plot(data['OmegaVec'], data['FTDipoleY_Im'], label='FTDipoleY_Im')
    plt.plot(data['OmegaVec'], data['FTDipoleZ_Im'], label='FTDipoleZ_Im')
    plt.title('Imaginary Components of FT Dipoles vs OmegaVec')
    plt.xlabel('OmegaVec')
    plt.ylabel('Imaginary Part of FT Dipole')
    plt.legend()

    # Plotting the average FT dipole strength
    plt.subplot(2, 1, 2)
    # plt.plot(omega_vec, average_ft_dipole, label='Average FT Dipole Magnitude', color='black')
    plt.plot(data['OmegaVec'], average_ft_dipole, label='Average FT Dipole Magnitude', color='black')
    plt.title('Average FT Dipole Magnitude vs OmegaVec')
    plt.xlabel('OmegaVec')
    plt.ylabel('Average FT Dipole Magnitude')
    plt.legend()

    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{study_directory}/{experiment_directory}/Dipole/DipoleFT_vs_OmegaVec.png')
    plt.close('all')

    if try_color_maps:
        color_maps = ['viridis', 'seismic', 'inferno', 'twilight_shifted']
    else:
        color_maps = ['seismic']
    # Create a plot per Atom
    for color_map in color_maps:
        """
        Alternative Plots for Understanding the Data (Plotly Contour Maps need to be 2D not 1D so we create mesh)
        """
        # Creating a grid for contour plot
        ct2_grid, omega_grid = np.meshgrid(data['central_time_2'].unique(), data['OmegaVec'].unique(), indexing='ij')

        # Reshaping the data to match the grid
        ft_dipole_x_re_grid = data['FTDipoleX_Re'].values.reshape(ct2_grid.shape)
        ft_dipole_x_im_grid = data['FTDipoleX_Im'].values.reshape(ct2_grid.shape)
        ft_dipole_y_re_grid = data['FTDipoleY_Re'].values.reshape(ct2_grid.shape)
        ft_dipole_y_im_grid = data['FTDipoleY_Im'].values.reshape(ct2_grid.shape)
        ft_dipole_z_re_grid = data['FTDipoleZ_Re'].values.reshape(ct2_grid.shape)
        ft_dipole_z_im_grid = data['FTDipoleZ_Im'].values.reshape(ct2_grid.shape)

        # Calculating the Z value for the contour plot (average FT dipole magnitude)
        average_ft_dipole_grid = np.sqrt((ft_dipole_x_re_grid ** 2 + ft_dipole_x_im_grid ** 2 +
                                          ft_dipole_y_re_grid ** 2 + ft_dipole_y_im_grid ** 2 +
                                          ft_dipole_z_re_grid ** 2 + ft_dipole_z_im_grid ** 2) / 3)

        fig, axes = plt.subplots(3, 2, figsize=(24, 18))

        # Dynamic Range Adjustment using different percentiles
        norm_1_99 = colors.Normalize(
            vmin=np.percentile(average_ft_dipole_grid, 0.01),
            vmax=np.percentile(average_ft_dipole_grid, 99.9))
        plt.colorbar(axes[0, 0].contourf(ct2_grid, omega_grid, average_ft_dipole_grid,
                                         levels=100, cmap=color_map, norm=norm_1_99), ax=axes[0, 0])
        axes[0, 0].set_title('Percentile (0.01% - 99.9%) Adjusted Contour Map')
        axes[0, 0].set_xlabel('Central Time 2')
        axes[0, 0].set_ylabel('OmegaVec')
        #
        norm_5_95 = colors.Normalize(
            vmin=np.percentile(average_ft_dipole_grid, 1),
            vmax=np.percentile(average_ft_dipole_grid, 99))
        plt.colorbar(axes[0, 1].contourf(ct2_grid, omega_grid, average_ft_dipole_grid,
                                         levels=100, cmap=color_map, norm=norm_5_95), ax=axes[0, 1])
        axes[0, 1].set_title('Percentile (1% - 99%) Adjusted Contour Map')
        axes[0, 1].set_xlabel('Central Time 2')
        axes[0, 1].set_ylabel('OmegaVec')

        # Logarithmic Scaling
        log_scaled_data = np.log1p(np.abs(average_ft_dipole_grid))  # Adding 1 to avoid log(0)
        plt.colorbar(axes[1, 0].contourf(ct2_grid, omega_grid, log_scaled_data, levels=100, cmap=color_map),
                     ax=axes[1, 0])
        axes[1, 0].set_title('Logarithmic Scaled Data Contour Map')
        axes[1, 0].set_xlabel('Central Time 2')
        axes[1, 0].set_ylabel('OmegaVec')

        # Derivative Plots
        plt.colorbar(axes[1, 1].contourf(ct2_grid, omega_grid, np.gradient(average_ft_dipole_grid, axis=0),
                                         levels=100, cmap=color_map), ax=axes[1, 1])
        axes[1, 1].set_title('First Derivative of Data')
        axes[1, 1].set_xlabel('Central Time 2')
        axes[1, 1].set_ylabel('OmegaVec')

        # Smoothing (Gaussian Filter) and Filtering
        smoothed_data = gaussian_filter(average_ft_dipole_grid, sigma=1)  # Applying a Gaussian filter
        plt.colorbar(axes[2, 0].contourf(ct2_grid, omega_grid, smoothed_data, levels=100, cmap=color_map),
                     ax=axes[2, 0])
        axes[2, 0].set_title('Smoothed Data Contour Map')
        axes[2, 0].set_xlabel('Central Time 2')
        axes[2, 0].set_ylabel('OmegaVec')

        # Derivative Plots
        plt.colorbar(axes[2, 1].contourf(ct2_grid, omega_grid, np.gradient(smoothed_data, axis=0), levels=100,
                                         cmap=color_map), ax=axes[2, 1])
        axes[2, 1].set_title('First Derivative of Smoothed Data')
        axes[2, 1].set_xlabel('Central Time 2')
        axes[2, 1].set_ylabel('OmegaVec')

        plt.tight_layout()
        # plt.show()
        plt.savefig(f'{study_directory}/{experiment_directory}/Dipole/DipoleFT_ALL_Alternative_{color_map}_Plots.png')
        plt.close('all')


def plot_ft_all_atomic_dipoles_v_time(study_directory, experiment_directory, dephasing_factor, relaxation_factor,
                                      pump_settings, probe_settings, charge_migration_ft_settings,
                                      xyz_geometry_path, try_color_maps=False):
    """
    This function plots the Fourier Transform of the atomic charges/dipoles vs time
    """
    # "number_of_pulses", "central_time_1", "carrier_frequency_1", "fwhm_1", "carrier_envelope_phase_1",
    #     "intensity_1", "amplitude_1", "period_1", "central_time_2", "carrier_frequency_2", "fwhm_2",
    #     "carrier_envelope_phase_2", "intensity_2", "amplitude_...", "period_...", "iOmega", "OmegaVec",
    #     "Atom_O_FTChargeX_Re", "Atom_O_FTChargeX_Im", "Atom_O_FTChargeY_Re", "Atom_O_FTChargeY_Im",
    #     "Atom_O_FTChargeZ_Re", "Atom_O_FTChargeZ_Im", "Atom_N_FTChargeX_Re", "Atom_N_FTChargeX_Im",
    #     "Atom_N_FTChargeY_Re", "Atom_N_FTChargeY_Im", "Atom_N_FTChargeZ_Re", "Atom_N_FTChargeZ_Im",
    #     "Atom_C_FTChargeX_Re", "Atom_C_FTChargeX_Im", "Atom_C_FTCharge..."

    # Read the XYZ file to get the number of atoms, the name of the atom and the atom names
    # (this is the order in which the atoms are listed in the CSV file)
    xyz_file_content = read_xyz(xyz_geometry_path)
    ATOMIC_FT_COL_NAMES = __generate_atomic_FTcharge_column_names(
        xyz_geometry_path)  # Rereads file but is intended for modularity and readability
    R_el_bc = pd.read_csv(f'{study_directory}/{experiment_directory}/R_el_bc.csv')

    FILE_PATH = f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL.csv'

    # Pump Settings
    type_of_pulse_pump = pump_settings['typeofpulse']
    start_time = pump_settings['starttime']
    pump_central_frequency = pump_settings['pumpcentralfrequency']
    pump_periods = pump_settings['pumpperiods']
    pump_phase = pump_settings['pumpphase']
    pump_intensity = pump_settings['pumpintensity']
    pump_polarization = pump_settings['pumppolarization']

    # Probe Settings
    type_of_pulse_probe = probe_settings['typeofpulse']
    time_delay_start = probe_settings['timedelaystart']
    time_delay_stop = probe_settings['timedelaystop']
    number_of_pp = probe_settings['numberofpp']
    time_Delay_weight_factor = probe_settings['timedelayweightfactor']
    probe_central_frequency = probe_settings['probecentralfrequency']
    probe_periods = probe_settings['probeperiods']
    probe_phase = probe_settings['probephase']
    probe_intensity = probe_settings['probeintensity']
    probe_polarization = probe_settings['probepolarization']

    # Charge Migration Settings
    ft_time_step = charge_migration_ft_settings['fttimestep']
    ft_width_step = charge_migration_ft_settings['ftwidthstep']

    # Lets plot the Spectra FT Dipolar Reponse vs Time (t)
    data = pd.read_csv(FILE_PATH)
    logger.debug(data.head())
    logger.debug(f'Length of data: {len(data)}')

    for i in ATOMIC_FT_COL_NAMES:
        if i not in data.columns:
            print(f'{i} not in data.columns')

    # Plotting the real and imaginary parts of the charge components for each atom
    import matplotlib.pyplot as plt

    plt.figure(figsize=(16, 12))
    plt.suptitle(FILE_PATH)
    for i, component in enumerate(['X_Re', 'X_Im', 'Y_Re', 'Y_Im', 'Z_Re', 'Z_Im']):
        columns_to_plot = [col_name for col_name in ATOMIC_FT_COL_NAMES if component in col_name]
        plt.subplot(3, 2, i + 1)
        for col_name in columns_to_plot:
            plt.plot(data['central_time_2'], data[col_name], label=col_name)
        plt.title(f'ChargeFT{component}')
        plt.xlabel('central_time_2')
        plt.ylabel(f'ChargeFT{component}')
        plt.legend()

    plt.tight_layout()
    # plt.show()
    plt.savefig(
        f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL_Components_vs_central_time_2.png')
    plt.close('all')

    # Calculating the average of the squared magnitudes of the FT dipole components
    sum_of_squares = 0
    for i in ATOMIC_FT_COL_NAMES:
        sum_of_squares += (data[i]) ** 2
    average_ft_atomic_dipole = np.sqrt(sum_of_squares / 3)

    # Creating an interactive plot using Plotly
    _plot_contour_map(x=data['central_time_2'],
                      y=data['OmegaVec'],
                      z=average_ft_atomic_dipole,
                      title=f'Contour Map of Average FT Dipole Magnitude\n'
                            f'AtomicChargeFT_ALL.csv\n'
                            f'Pump Settings:  {pump_central_frequency}, {pump_periods}, {pump_phase}, {pump_intensity}, {pump_polarization}\n'
                            f'Probe Settings:  {probe_central_frequency}, {probe_periods}, {probe_phase}, {probe_intensity}, {probe_polarization}\n'
                            f'Dephasing Factor: {dephasing_factor}, Relaxation Factor: {relaxation_factor}'
                            f'FT Time Step: {ft_time_step}, FT Width Step: {ft_width_step}',
                      x_label='Central Time 2 - Probe (tau)',
                      y_label='OmegaVec',
                      output_file=f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL_Contour.html')

    # Plotting
    plt.figure(figsize=(15, 10))

    # Plotting the real parts
    plt.subplot(2, 2, 1)
    for i in ATOMIC_FT_COL_NAMES:
        if i not in data.columns:
            print(f'{i} not in data.columns')
        if '_Re' in i:
            plt.plot(data['OmegaVec'], data[i], label=i)
    plt.title('Real Components of FT Atomic Charges vs OmegaVec')
    plt.xlabel('OmegaVec')
    plt.ylabel('Real Part of FT Atomic Charge')

    # Plotting the imaginary parts
    plt.subplot(2, 2, 2)
    for i in ATOMIC_FT_COL_NAMES:
        if '_Im' in i:
            plt.plot(data['OmegaVec'], data[i], label=i)
    plt.title('Imaginary Components of FT Atomic Charges vs OmegaVec')
    plt.xlabel('OmegaVec')
    plt.ylabel('Imaginary Part of FT Atomic Charge')
    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc='best',
        borderaxespad=0.0,
        title='Legend',
        title_fontsize=12,
        fontsize=12,
        labelspacing=0.1,
        fancybox=True,
        shadow=True,
        ncol=1,
        framealpha=0.5
    )

    # Plotting the average FT dipole strength
    plt.subplot(2, 1, 2)
    # plt.plot(omega_vec, average_ft_dipole, label='Average FT Dipole Magnitude', color='black')
    plt.plot(data['OmegaVec'], average_ft_atomic_dipole, label='Average FT Atomic Dipole Magnitude', color='black')
    plt.title('Average FT Atomic Dipole Magnitude vs OmegaVec')
    plt.xlabel('OmegaVec')
    plt.ylabel('Average FT Atomic Dipole Magnitude')

    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_vs_OmegaVec.png')
    plt.close('all')

    ############
    """
    Alternative Plots for Understanding the Data (Plotly Contour Maps need to be 2D not 1D so we create mesh)
    """
    # Creating a grid for contour plot
    ct2_grid, omega_grid = np.meshgrid(data['central_time_2'].unique(), data['OmegaVec'].unique(), indexing='ij')

    if try_color_maps:
        color_maps = ['viridis', 'seismic', 'inferno', 'twilight_shifted']
    else:
        color_maps = ['seismic']
    # Create a plot per Atom
    for color_map in color_maps:
        for i, _atom_name in enumerate(xyz_file_content['atoms'].keys()):
            fig, axes = plt.subplots(3, 2, figsize=(24, 18))
            # Set the title for the plot
            fig.suptitle(FILE_PATH)
            split_atom_name = _atom_name.split('_')
            if len(split_atom_name) == 1:
                atom_name = _atom_name
                atom_subscript = ''
            elif len(split_atom_name) == 2:
                atom_name = split_atom_name[0]
                atom_subscript = f'.{split_atom_name[1]}'
            else:
                logger.error(f'Atom name {_atom_name} is not valid')

            # Reshaping the data to match the grid
            atomic_ft_dipole_x_re_grid = data[f'Atom_{atom_name}_FTChargeX_Re{atom_subscript}'].values.reshape(
                ct2_grid.shape)
            atomic_ft_dipole_x_im_grid = data[f'Atom_{atom_name}_FTChargeX_Im{atom_subscript}'].values.reshape(
                ct2_grid.shape)
            atomic_ft_dipole_y_re_grid = data[f'Atom_{atom_name}_FTChargeY_Re{atom_subscript}'].values.reshape(
                ct2_grid.shape)
            atomic_ft_dipole_y_im_grid = data[f'Atom_{atom_name}_FTChargeY_Im{atom_subscript}'].values.reshape(
                ct2_grid.shape)
            atomic_ft_dipole_z_re_grid = data[f'Atom_{atom_name}_FTChargeZ_Re{atom_subscript}'].values.reshape(
                ct2_grid.shape)
            atomic_ft_dipole_z_im_grid = data[f'Atom_{atom_name}_FTChargeZ_Im{atom_subscript}'].values.reshape(
                ct2_grid.shape)

            # Calculating the Z value for the contour plot (average FT dipole magnitude)
            average_ft_atomic_dipole_grid = np.sqrt((atomic_ft_dipole_x_re_grid ** 2 + atomic_ft_dipole_x_im_grid ** 2 +
                                                     atomic_ft_dipole_y_re_grid ** 2 + atomic_ft_dipole_y_im_grid ** 2 +
                                                     atomic_ft_dipole_z_re_grid ** 2 + atomic_ft_dipole_z_im_grid ** 2) / 3)

            # Dynamic Range Adjustment using different percentiles
            norm_1_99 = colors.Normalize(
                vmin=np.percentile(average_ft_atomic_dipole_grid, 0.01),
                vmax=np.percentile(average_ft_atomic_dipole_grid, 99.9))
            plt.colorbar(axes[0, 0].contourf(ct2_grid, omega_grid, average_ft_atomic_dipole_grid,
                                             levels=100, cmap=color_map, norm=norm_1_99), ax=axes[0, 0])
            axes[0, 0].set_title(
                f'Percentile (0.01% - 99.9%) Adjusted Contour Map of Atom {_atom_name} FT Charge Averaged')
            axes[0, 0].set_xlabel('Central Time 2')
            axes[0, 0].set_ylabel('OmegaVec')

            #
            norm_5_95 = colors.Normalize(
                vmin=np.percentile(average_ft_atomic_dipole_grid, 5),
                vmax=np.percentile(average_ft_atomic_dipole_grid, 95))
            plt.colorbar(axes[0, 1].contourf(ct2_grid, omega_grid, average_ft_atomic_dipole_grid,
                                             levels=100, cmap=color_map, norm=norm_5_95), ax=axes[0, 1])
            axes[0, 1].set_title(f'Percentile (1% - 99%) Adjusted Contour Map of Atom {_atom_name} FT Charge Averaged')
            axes[0, 1].set_xlabel('Central Time 2')
            axes[0, 1].set_ylabel('OmegaVec')

            # Logarithmic Scaling
            log_scaled_data = np.log1p(np.abs(average_ft_atomic_dipole_grid))
            plt.colorbar(axes[1, 0].contourf(ct2_grid, omega_grid, log_scaled_data, levels=100, cmap=color_map),
                         ax=axes[1, 0])
            axes[1, 0].set_title(f'Logarithmic Scaled Data Contour Map of Atom {_atom_name} FT Charge Averaged')
            axes[1, 0].set_xlabel('Central Time 2')
            axes[1, 0].set_ylabel('OmegaVec')

            # Derivative Plots
            plt.colorbar(axes[1, 1].contourf(ct2_grid, omega_grid, np.gradient(average_ft_atomic_dipole_grid, axis=0),
                                             levels=100, cmap=color_map), ax=axes[1, 1])
            axes[1, 1].set_title(f'First Derivative of Data of Atom {_atom_name} FT Charge Averaged')
            axes[1, 1].set_xlabel('Central Time 2')
            axes[1, 1].set_ylabel('OmegaVec')

            # Smoothing (Gaussian Filter) and Filtering
            smoothed_data = gaussian_filter(average_ft_atomic_dipole_grid, sigma=1)  # Applying a Gaussian filter
            plt.colorbar(axes[2, 0].contourf(ct2_grid, omega_grid, smoothed_data, levels=100, cmap=color_map),
                         ax=axes[2, 0])
            axes[2, 0].set_title(f'Smoothed Data Contour Map of Atom {_atom_name} FT Charge Averaged')
            axes[2, 0].set_xlabel('Central Time 2')

            # Derivative Plots
            plt.colorbar(axes[2, 1].contourf(ct2_grid, omega_grid, np.gradient(smoothed_data, axis=0), levels=100,
                                             cmap=color_map), ax=axes[2, 1])
            axes[2, 1].set_title(f'First Derivative of Smoothed Data of Atom {_atom_name} FT Charge Averaged')
            axes[2, 1].set_xlabel('Central Time 2')
            axes[2, 1].set_ylabel('OmegaVec')

            plt.tight_layout()
            # plt.show()
            plt.savefig(
                f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL_Alternative_{color_map}_Plots_{_atom_name}.png')
            plt.close('all')


'''Dipolar 2D Spectra'''


def plot_2d_spectrum(study_directory, experiment_directory, dephasing_factor, relaxation_factor,
                     pump_settings, probe_settings, charge_migration_ft_settings, match_scales=False):
    OMEGA_TAUOMEGA_FILE_NAMES = [
        'Dipole/DipoleFT_ww.csv',
        'Dipole/DipoleFT_ww_reconstructed.csv',
    ]
    OMEGA_TAUOMEGA_FILE_PATHS = [f'{study_directory}/{experiment_directory}/{file_name}' for file_name in
                                 OMEGA_TAUOMEGA_FILE_NAMES]

    length_of_data_to_match = None
    FEATURES_COL_NAMES = [  # Todo do per polarization images too
        '2DDipoleX_Re',
        '2DDipoleX_Im',
        '2DDipoleY_Re',
        '2DDipoleY_Im',
        '2DDipoleZ_Re',
        '2DDipoleZ_Im'
    ]

    ft_time_step = charge_migration_ft_settings['fttimestep']
    ft_width_step = charge_migration_ft_settings['ftwidthstep']

    _zmax = None
    _zmin = None
    for file_path in OMEGA_TAUOMEGA_FILE_PATHS:
        # Load the CSV file
        logger.info(f'Plotting {file_path}')
        data = pd.read_csv(file_path)

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(
                data), f'Length of {file_path} is not the same as the previous file'

        # Displaying the first few rows of the file to understand its structure
        logger.debug(data.head())
        logger.debug(f'Length of data: {len(data)}')

        # Calculating the average of all polarizations
        data['AveragedDensity'] = (sum(data[col] ** 2 for col in FEATURES_COL_NAMES) / 3) ** 0.5
        # data['AveragedDensity'] = sum(data[col] ** 2 for col in FEATURES_COL_NAMES)

        # Adding the averaged density to the DataFrame
        if 'DipoleFT_ww_reconstructed' in file_path:
            # FIXME Flip accross the frequency axis but keep the same index (idk why this bug exists)
            data['AveragedDensity'] = pd.Series(data['AveragedDensity'].values[::-1],
                                                index=data['AveragedDensity'].index)

        # Plotting the contour map
        if match_scales:
            if _zmax is None:
                _zmax = np.percentile(data['AveragedDensity'], 99.5)
                _zmin = np.percentile(data['AveragedDensity'], 0.01)
        else:
            _zmax = np.percentile(data['AveragedDensity'], 99.9)
            _zmin = np.percentile(data['AveragedDensity'], 0.01)

        _plot_contour_map(
            x=data['TauOmegaVec'] + data['OmegaVec'],
            y=data['OmegaVec'],
            z=data['AveragedDensity'],
            title=f'Contour Map of Average FT Dipole (Pure-Absorption)\n'
                  f'{file_path}\n'
                  f'Pump Settings:  {pump_settings["pumpcentralfrequency"]}, {pump_settings["pumpperiods"]}, {pump_settings["pumpphase"]}, {pump_settings["pumpintensity"]}, {pump_settings["pumppolarization"]}\n'
                  f'Probe Settings:  {probe_settings["probecentralfrequency"]}, {probe_settings["probeperiods"]}, {probe_settings["probephase"]}, {probe_settings["probeintensity"]}, {probe_settings["probepolarization"]}\n'
                  f'FT Time Step: {ft_time_step}, FT Width Step: {ft_width_step}, Dephasing Factor: {dephasing_factor}, Relaxation Factor: {relaxation_factor}',
            x_label='OmegaVec + TauOmegaVec',
            y_label='OmegaVec',
            output_file=f'{file_path}'.replace('.csv', '.html'),
            add_diagonal_line=True,
            contour_kwargs={
                'colorscale': 'Viridis',
                'contours_coloring': 'heatmap',
                'zmin': _zmin,
                'zmax': _zmax
            }
        )


def plot_2d_spectrum_peak_analysis(study_directory, experiment_directory):
    OMEGA_TAUOMEGA_FILE_NAMES = [
        'Dipole/DipoleFT_ww.csv',
        'Dipole/DipoleFT_ww_reconstructed.csv',
    ]
    OMEGA_TAUOMEGA_FILE_PATHS = [f'{study_directory}/{experiment_directory}/{file_name}' for file_name in
                                 OMEGA_TAUOMEGA_FILE_NAMES]

    length_of_data_to_match = None
    FEATURES_COL_NAMES = [
        '2DDipoleX_Re',
        '2DDipoleX_Im',
        '2DDipoleY_Re',
        '2DDipoleY_Im',
        '2DDipoleZ_Re',
        '2DDipoleZ_Im'
    ]

    for file_path in OMEGA_TAUOMEGA_FILE_PATHS:
        # Load the CSV file
        logger.info(f'Plotting {file_path}')
        data = pd.read_csv(file_path)

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(
                data), f'Length of {file_path} is not the same as the previous file'

        # Displaying the first few rows of the file to understand its structure
        logger.debug(data.head())
        logger.debug(f'Length of data: {len(data)}')

        # Calculating the average of all polarizations
        # data['AveragedDensity'] = sum(data[col] ** 2 for col in FEATURES_COL_NAMES)
        data['AveragedDensity'] = (sum(data[col] ** 2 for col in FEATURES_COL_NAMES) / 3) ** 0.5

        if 'DipoleFT_ww_reconstructed' in file_path:
            # FIXME Flip accross the frequency axis but keep the same index (idk why this bug exists)
            data['AveragedDensity'] = pd.Series(data['AveragedDensity'].values[::-1],
                                                index=data['AveragedDensity'].index)

        import matplotlib.pyplot as plt

        # Preparing the data for plotting
        omega_values = data['OmegaVec'].unique()
        tauomega_values = data['TauOmegaVec'].unique()

        # Creating a meshgrid for plotting
        omega_grid, tauomega_grid = np.meshgrid(omega_values, tauomega_values, indexing='ij')

        # Reshaping 'AveragedDensity' to match the shape of the meshgrid
        averaged_density_grid = data['AveragedDensity'].values.reshape(omega_grid.shape)

        ##############
        from scipy.signal import find_peaks
        from sklearn.preprocessing import MinMaxScaler

        # Adjusting peak analysis to be sensitive to both positive and negative values
        # Here, we'll find peaks on the absolute values of the max density along omega
        peaks_positive, _ = find_peaks(averaged_density_grid.max(axis=0))  # positive peaks
        peaks_negative, _ = find_peaks(-averaged_density_grid.max(axis=0))  # negative peaks

        # Dynamic range adjustment for the contour plots
        # Using percentiles to set color scale limits
        percentile_min = np.percentile(averaged_density_grid, 0.01)
        percentile_max = np.percentile(averaged_density_grid, 99.9)

        # Symmetrical color map for first derivative
        first_derivative = np.gradient(averaged_density_grid, axis=0)  # recalculating first derivative
        derivative_min = np.percentile(first_derivative, 0.01)
        derivative_max = np.percentile(first_derivative, 99.9)
        derivative_abs_max = max(abs(derivative_min), abs(derivative_max))

        # Plotting with adjustments
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        plt.suptitle(file_path)

        # Adjusted Peak Analysis Plot
        axes[0, 0].plot(omega_values, averaged_density_grid.max(axis=0))
        axes[0, 0].plot(omega_values[peaks_positive], averaged_density_grid.max(axis=0)[peaks_positive], "x",
                        color='blue')
        axes[0, 0].plot(omega_values[peaks_negative], averaged_density_grid.max(axis=0)[peaks_negative], "x",
                        color='red')
        axes[0, 0].set_title("Peak Analysis" + "\n" + "Blue: Positive Peaks, Red: Negative Peaks")

        # Integrated Density
        scaler = MinMaxScaler()
        normalized_density = scaler.fit_transform(averaged_density_grid)
        integrated_density = np.trapz(normalized_density, axis=0)  # integrating along tauomega
        # The plot for integrated density remains unchanged
        axes[0, 1].plot(omega_values, integrated_density)
        axes[0, 1].set_title("Integrated Density")

        # Adjusted First Derivative
        plt.colorbar(axes[1, 0].contourf(omega_grid + tauomega_grid, omega_grid, first_derivative,
                                         levels=100,
                                         vmin=-derivative_abs_max, vmax=derivative_abs_max,
                                         cmap='viridis'), ax=axes[1, 0])
        axes[1, 0].set_title("First Derivative")

        # Adjusted Normalized Density
        plt.colorbar(axes[1, 1].contourf(omega_grid + tauomega_grid, omega_grid, averaged_density_grid,
                                         levels=100, cmap='viridis',
                                         norm=colors.Normalize(vmin=percentile_min, vmax=percentile_max)
                                         ), ax=axes[0, 1])
        axes[1, 1].set_title("Normalized Density (0.01% and 99.9% Percentiles)")

        plt.tight_layout()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        output_file = output_file.replace('.png', f'_Alternative_Plots.png')
        plt.savefig(output_file)
        # plt.show()
        plt.close('all')


'''__helper_functions__'''


def _plot_contour_map(x, y, z, title, x_label, y_label, output_file, add_diagonal_line=False, contour_kwargs=None):
    if contour_kwargs is None:
        contour_kwargs = {}
    # Creating an interactive plot using Plotly
    import plotly.graph_objects as go
    fig = go.Figure(data=go.Contour(z=z, x=x, y=y,
                                    **contour_kwargs))
    fig.update_layout(title=title, xaxis_title=x_label, yaxis_title=y_label)

    if add_diagonal_line:
        fig.add_shape(
            type="line",
            x0=min(y),
            y0=min(y),
            x1=max(y),
            y1=max(y),
            line=dict(
                color="white",
                width=2,
                dash="dot",
            ),
        )

    fig.write_html(output_file)
    # fig.show()

    return fig


def _calculate_min_max_significant_time(data, target_column, percentile_range=(20, 60)):
    """
    Calculate the min and max time for the significant data points using percentiles.

    Expected data format:
    Time, ..., target_column, Foo1, Foo2, Foo3, ...

    :param data: DataFrame containing the data.
    :param target_column: The name of the column with the target data.
    :param percentile_range: Tuple containing the lower and upper percentiles to consider.
    :return: Tuple with the min and max times corresponding to the significant data range.
    """
    assert target_column in data.columns, f'Column {target_column} not found in {data.columns}'

    # Calculate the lower and upper bounds using percentiles
    lower_bound = data[target_column].quantile(percentile_range[0] / 100)
    upper_bound = data[target_column].quantile(percentile_range[1] / 100)

    # Filter data within the percentile range
    significant_data = data[(data[target_column] >= lower_bound) & (data[target_column] <= upper_bound)]

    # Determine the range for the x-axis
    min_time = significant_data['Time'].min()
    max_time = significant_data['Time'].max()

    return min_time, max_time


def __generate_atomic_charge_column_names(file_path):
    with open(file_path, 'r') as f:
        number_of_atoms = int(f.readline())
        _ = f.readline()  # Skipping the molecule name line
        atom_names = [f.readline().split()[0] for _ in range(number_of_atoms)]

    # Initialize counters for each atom name
    atom_name_counts = {name: 0 for name in set(atom_names)}

    # Initialize the list of column names
    column_names = []

    # Generate column names with suffixes for repeated atom names
    for name in atom_names:
        count = atom_name_counts[name]
        suffix = f".{count}" if count > 0 else ""
        for axis in ['X', 'Y', 'Z']:
            column_names.append(f'Atom_{name}_Charge{axis}{suffix}')
        atom_name_counts[name] += 1

    return column_names


def __generate_atomic_FTcharge_column_names(file_path):
    xyz_file_content = read_xyz(file_path)
    atoms = xyz_file_content['atoms']

    # Initialize counters for each atom name
    atom_name_counts = {name: 0 for name in set([atom_name.split('_')[0] for atom_name in atoms.keys()])}

    # Initialize the list of column names
    column_names = []

    # Generate column names with suffixes for repeated atom names
    for name in atoms.keys():
        _name = name.split('_')[0]
        count = atom_name_counts[_name]
        suffix = f".{count}" if count > 0 else ""
        for axis in ['X', 'Y', 'Z']:
            for component in ['Re', 'Im']:
                column_names.append(f'Atom_{_name}_FTCharge{axis}_{component}{suffix}')
        atom_name_counts[_name] += 1

    return column_names


def __generate_becke_atomic_weghts_column_names(xyz_geometry_path):
    xyz_file_content = read_xyz(xyz_geometry_path)
    atoms = xyz_file_content['atoms']

    # Initialize counters for each atom name
    atom_name_counts = {name: 0 for name in set([atom_name.split('_')[0] for atom_name in atoms.keys()])}

    # Initialize the list of column names
    column_names = []

    # Generate column names with suffixes for repeated atom names
    for name in atoms.keys():
        _name = name.split('_')[0]
        count = atom_name_counts[_name]
        suffix = f".{count}" if count > 0 else ""
        column_names.append(f'Atom_{_name}_ChargeDensity{suffix}')
        atom_name_counts[_name] += 1

    return column_names


'''Becke Plots'''


def plot_becke_weights(study_directory, experiment_directory, xyz_geometry_path, weights_file):
    # "x","y","z","Atom_O_ChargeDensity","Atom_N_ChargeDensity","Atom_C_ChargeDensity","Atom_C_ChargeDensity",
    # "Atom_C_ChargeDensity","Atom_H_ChargeDensity","Atom_H_ChargeDensity","Atom_H_ChargeDensity",
    # "Atom_H_ChargeDensity","Atom_H_ChargeDensity","Atom_H_ChargeDensity","Atom_H_ChargeDensity"
    if weights_file is None:
        # Use the expected default from the ChargeMigration/Becke code
        weights_file = f'{study_directory}/{experiment_directory}/Weights_File_{experiment_directory}.csv'

    # Read the XYZ file to get the number of atoms, the name of the atom and the atom names
    # (this is the order in which the atoms are listed in the CSV file)
    xyz_file_content = read_xyz(xyz_geometry_path)
    number_of_atoms = xyz_file_content['n_atoms']
    name_of_molecule = xyz_file_content['title']
    atoms = xyz_file_content['atoms']
    ATOMIC_FT_COL_NAMES = __generate_becke_atomic_weghts_column_names(xyz_geometry_path)

    print(f'{number_of_atoms = }')
    print(f'{ATOMIC_FT_COL_NAMES = }')

    assert number_of_atoms == len(
        ATOMIC_FT_COL_NAMES), f'{number_of_atoms = } != {len(ATOMIC_FT_COL_NAMES) = } check read_xyz and __generate_becke_atomic_weghts_column_names'

    # Load the CSV file
    logger.info(f'Plotting {weights_file}')
    data = pd.read_csv(weights_file)

    # Displaying the first few rows of the file to understand its structure
    # all columns
    pd.set_option('display.max_columns', None)
    print(data.head())
    logger.debug(f'Length of Weight data: {len(data)} points')

    # Create a 3d plot of the weights with each atom/column as a different color scaled by the density/value of the point
    # Plotting
    import plotly.graph_objects as go

    for col in data.columns[3:]:
        print(f'{col} = {data[col].sum()}')
        if data[col].min() != 0:
            logger.warning(f'{col} has a min value of {data[col].min()} rather than 0')
        if data[col].max() != 1:
            logger.warning(f'{col} has a max value of {data[col].max()} rather than 1')

    # Plot only one atom type for debugging
    fig = go.Figure(data=go.Scatter3d(
        x=data['x'],
        y=data['y'],
        z=data['z'],
        mode='markers',
        marker=dict(
            size=3,
            color=data['Atom_H_ChargeDensity.4'],  # Replace with the correct column name
            colorscale='Viridis',
            opacity=0.5
        )
    ))

    fig.update_layout(title='3D Scatter Plot for Atom_H_ChargeDensity.4')
    fig.show()


'''Correlation Graphs'''


def plot_dipole_correlation_maps(study_directory, experiment_directory):
    import seaborn as sns
    import matplotlib.pyplot as plt

    '''FT WW'''
    atomic_ft_ww = pd.read_csv(f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_ww.csv')
    plt.figure(figsize=(12, 12))
    atomic_ft_ww_real = atomic_ft_ww[[col for col in atomic_ft_ww.columns if '_Re' in col]]
    sns.heatmap(atomic_ft_ww_real.corr(), annot=True, cmap='viridis', vmin=-1, vmax=1, center=0, square=True,
                linewidths=.5, cbar_kws={"shrink": .5}, )
    plt.title('Atomic Charge FT ww Correlation Matrix')
    plt.tight_layout()
    plt.savefig(f'{study_directory}/{experiment_directory}/AtomicChargeFT_ww_CorrelationMatrix_Real.png')

    '''FT ALL'''
    atomic_ft_all = pd.read_csv(f'{study_directory}/{experiment_directory}/AtomicCharge/AtomicChargeFT_ALL.csv')
    plt.figure(figsize=(12, 12))
    atomic_ft_all_real = atomic_ft_all[[col for col in atomic_ft_all.columns if '_Re' in col]]

    sns.heatmap(atomic_ft_all_real.corr(), annot=True, cmap='viridis', vmin=-1, vmax=1, center=0, square=True,
                linewidths=.5, cbar_kws={"shrink": .5}, )
    plt.title('Atomic Charge FT ALL Correlation Matrix')
    plt.tight_layout()
    plt.savefig(f'{study_directory}/{experiment_directory}/AtomicChargeFT_ALL_CorrelationMatrix_Real.png')

    '''Dipole FT WW and Dipole FT WW Reconstructed'''
    OMEGA_TAUOMEGA_FILE_NAMES = [
        'Dipole/DipoleFT_ww.csv',
        'Dipole/DipoleFT_ww_reconstructed.csv',
    ]
    OMEGA_TAUOMEGA_FILE_PATHS = [f'{study_directory}/{experiment_directory}/{file_name}' for file_name in
                                 OMEGA_TAUOMEGA_FILE_NAMES]

    length_of_data_to_match = None
    FEATURES_COL_NAMES = [
        '2DDipoleX_Re',
        '2DDipoleX_Im',
        '2DDipoleY_Re',
        '2DDipoleY_Im',
        '2DDipoleZ_Re',
        '2DDipoleZ_Im'
    ]
    datas = []  # Make sure to vary the names by adding _reconstructed
    suffix_names = []
    for file_path in OMEGA_TAUOMEGA_FILE_PATHS:
        suffix_name = file_path.split('/')[-1].replace('.csv', '')
        suffix_names.append(suffix_name)

        # Load the CSV file
        logger.info(f'Plotting {file_path}')
        data = pd.read_csv(file_path)

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(
                data), f'Length of {file_path} is not the same as the previous file'

        # Displaying the first few rows of the file to understand its structure
        logger.debug(data.head())
        logger.debug(f'Length of data: {len(data)}')

        # Calculating the average of all polarizations
        # data['AveragedDensity'] = sum(data[col] ** 2 for col in FEATURES_COL_NAMES)
        data['AveragedDensity'] = (sum(data[col] ** 2 for col in FEATURES_COL_NAMES) / 3) ** 0.5

        if 'DipoleFT_ww_reconstructed' in file_path:
            # FIXME Flip accross the frequency axis but keep the same index (idk why this bug exists)
            data['AveragedDensity'] = pd.Series(data['AveragedDensity'].values[::-1],
                                                index=data['AveragedDensity'].index)

        datas.append(data)

    plt.figure(
        # Two figures
        figsize=(20, 10)
    )

    for i, data in enumerate(datas):
        # Drop the omega columns
        plot_data = data.drop(columns=['OmegaVec', 'TauOmegaVec', 'AveragedDensity'])
        plt.subplot(1, 2, i + 1)
        plt.title('Correlation Matrix ' + suffix_names[i])
        sns.heatmap(plot_data.corr(),
                    annot=True, cmap='viridis', vmin=-1, vmax=1, center=0,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5},
                    # limit how many decimal places are shown
                    fmt='.1f'

                    )

    plt.tight_layout()
    plt.savefig(f'{study_directory}/{experiment_directory}/2DDipoleCorrelationMatrix.png')


def plot_cm_correlation_maps(study_directory, experiment_directory):
    """This function is still in development and is not yet complete. It should be used as a reference for future work."""
    import seaborn as sns
    import matplotlib.pyplot as plt

    '''FT WW'''
    atomic_cm = pd.read_csv(f'{study_directory}/{experiment_directory}/ChargeDensity/ChDenSimPP-200/ChDen54.3465.csv')
    plt.figure(figsize=(12, 12))
    # Remove the axis columns
    plot_atomic_cm = atomic_cm.drop(columns=['x', 'y', 'z'])
    sns.heatmap(plot_atomic_cm.corr(), annot=True,
                # We are plotting the correlation matrix thus i need a diverging color map but better than viridis
                # cmap='viridis',
                # cmap='seismic',
                # cmap='coolwarm',
                # cmap='twilight_shifted',
                cmap='vlag',
                vmin=-1, vmax=1, center=0, square=True,
                linewidths=.5, cbar_kws={"shrink": .5},
                )
    plt.title('Charge Density Correlation Matrix')
    plt.tight_layout()
    plt.savefig(f'{study_directory}/{experiment_directory}/ChargeDensity_CorrelationMatrix.png')
    # plt.show()

    # exit()
    # sns.pairplot(plot_atomic_cm)
    # plt.savefig(f'{study_directory}/{experiment_directory}/ChargeDensity_PairPlot.png')

    # sns.jointplot(data=plot_atomic_cm)
    # plt.savefig(f'{study_directory}/{experiment_directory}/ChargeDensity_JointPlot.png')


def generate_heatmap(csv_file_name: str, data_directory_path: str):
    atomic_cm = pd.read_csv(os.path.join(data_directory_path, csv_file_name))
    plt.figure(figsize=(12, 12))
    plot_atomic_cm = atomic_cm.drop(columns=['x', 'y', 'z'])
    sns.heatmap(plot_atomic_cm.corr(), annot=True, cmap='vlag', vmin=-1, vmax=1, center=0, square=True,
                linewidths=.5, cbar_kws={"shrink": .5})
    plt.title(f'Charge Density Correlation Matrix ({csv_file_name})')
    plt.tight_layout()

    # Save the heatmap
    heatmap_path = os.path.join(data_directory_path,
                                f'heatmap_{csv_file_name}'.replace('.gz', '').replace('.csv', '.png'))
    print(f'Saving heatmap to {heatmap_path}')
    plt.savefig(heatmap_path)
    plt.close()
    return heatmap_path


def plot_cm_correlation_maps_with_time(study_directory, experiment_directory, cm_pp_dir_name):
    charge_migration_experiment_directory = f'{study_directory}/{experiment_directory}'
    # Directory containing the CSV files
    data_directory = f'{charge_migration_experiment_directory}/ChargeDensity/{cm_pp_dir_name}'

    # List of CSV files sorted by time
    csv_files = sorted(os.listdir(data_directory), key=lambda x: float(x.split('.')[0].split('ChDen')[-1]))
    csv_files = [file for file in csv_files if file.endswith('.csv.gz')]

    assert len(csv_files) > 0, f'No CSV.GZ files found in {data_directory}'

    # Use multiprocessing to generate heatmaps in parallel
    with Pool(processes=32) as pool:
        heatmap_paths = pool.starmap(generate_heatmap, [(csv_file, data_directory) for csv_file in csv_files])

    print('Heatmaps generated!')
    print(f'{heatmap_paths = }')
    # Create a GIF from the heatmaps
    with imageio.get_writer(f'{charge_migration_experiment_directory}/{cm_pp_dir_name}_correlation_heatmaps.gif',
                            mode='I', duration=0.5) as writer:
        print('Creating GIF...')
        for heatmap_path in heatmap_paths:
            image = imageio.imread(heatmap_path)
            writer.append_data(image)
            # Remove the heatmap image to save space
            os.remove(heatmap_path)

    print(
        f'GIF created! Location: {charge_migration_experiment_directory}/ChargeDensity/{cm_pp_dir_name}_correlation_heatmaps.gif')


