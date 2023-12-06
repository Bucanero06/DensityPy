import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])


DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING = {
    'DipoleX_Re': 'blue',
    'DipoleY_Re': 'green',
    'DipoleZ_Re': 'red',
    'DipoleX_Im': 'lightblue',
    'DipoleY_Im': 'lightgreen',
    'DipoleZ_Im': 'pink',
}


def calculate_min_max_significant_time(data, target_column, percentile_range=(20, 60)):
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


def plot_pulses(study_directory, experiment_directory, time_delays,  min_time, max_time, plot_all=False):

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

        logger.info(f'Plotting {file_path}')
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

    plt.clf()


def plot_ft_pulses(study_directory, experiment_directory, time_delays, plot_all=False):
    # FTpulsePP{delay}
    # Freq: The frequency at which the Fourier transform is computed.
    # FTAx: The x-component of the Fourier-transformed vector potential.
    # FTAy: The y-component of the Fourier-transformed vector potential.
    # FTAz: The z-component of the Fourier-transformed vector potential.
    # FT_Aminus1_Real and FT_Aminus1_Imag: The real and imaginary parts of the Fourier transform for the mu = -1 component.
    # FT_A0_Real and FT_A0_Imag: The real and imaginary parts of the Fourier transform for the mu = 0 component.
    # FT_Aplus1_Real and FT_Aplus1_Imag: The real and imaginary parts of the Fourier transform for the mu = +1 component.

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

    # Clear the plot
    plt.clf()


def plot_dipoles_v_time(study_directory, experiment_directory, time_delays, min_time, max_time, plot_all=False):
    # Lets plot the Dipolar Reponse vs Time (t)
    # "itime","Time","DipoleX_Re","DipoleX_Im","DipoleY_Re","DipoleY_Im","DipoleZ_Re","DipoleZ_Im"

    length_of_data_to_match = None
    INDECES_COL_NAMES = ['itime', 'Time']
    FEATURES_COL_NAMES = [
        'DipoleX_Re',
        'DipoleX_Im',
        'DipoleY_Re',
        'DipoleY_Im',
        'DipoleZ_Re',
        'DipoleZ_Im'
    ]

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
        file_path = f'{study_directory}/{experiment_directory}/Dipole/DipolePP{time_delay}.csv' if time_delay != 'XUV' else f'{study_directory}/{experiment_directory}/Dipole/DipoleXUV.csv'

        logger.info(f'Plotting {file_path}')
        try:
            data = pd.read_csv(file_path)
        except FileNotFoundError:
            logger.error(f'File {file_path} not found')
            continue

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'
            assert data.columns.all() in INDECES_COL_NAMES + FEATURES_COL_NAMES, f'Column names of {file_path} is not the same as the previous file'

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

        plt.clf()

    plt.clf()


def plot_atomic_dipoles_v_time(study_directory, experiment_directory, time_delays, plot_all=False):
    # Lets plot the Dipolar Reponse vs Time (t)
    # "itime","Time","DipoleX_Re","DipoleX_Im","DipoleY_Re","DipoleY_Im","DipoleZ_Re","DipoleZ_Im"
    ...





def plot_2d_spectrum(study_directory, experiment_directory):
    OMEGA_TAUOMEGA_FILE_NAMES = [
        'Dipole/DipoleFT_ww.csv',
        'Dipole/DipoleFT_ww_reconstructed.csv',
    ]
    OMEGA_TAUOMEGA_FILE_PATHS = [f'{study_directory}/{experiment_directory}/{file_name}' for file_name in
                                 OMEGA_TAUOMEGA_FILE_NAMES]

    length_of_data_to_match = None
    INDECES_COL_NAMES = ['OmegaVec', 'TauOmegaVec']
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
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'
            assert data.columns.all() in INDECES_COL_NAMES + FEATURES_COL_NAMES, f'Column names of {file_path} is not the same as the previous file'

        # Displaying the first few rows of the file to understand its structure
        logger.debug(data.head())
        logger.debug(f'Length of data: {len(data)}')

        # Calculating the average of all polarizations

        averaged_density = (sum(data[col] ** 2 for col in FEATURES_COL_NAMES) / 3) ** 0.5

        # Adding the averaged density to the DataFrame
        data['AveragedDensity'] = averaged_density

        import matplotlib.pyplot as plt
        import numpy as np

        # Preparing the data for plotting
        omega_values = data['OmegaVec'].unique()
        tauomega_values = data['TauOmegaVec'].unique()

        # Creating a meshgrid for plotting
        omega_grid, tauomega_grid = np.meshgrid(omega_values, tauomega_values, indexing='ij')

        # Reshaping 'AveragedDensity' to match the shape of the meshgrid
        averaged_density_grid = data['AveragedDensity'].values.reshape(omega_grid.shape)

        # Plotting
        plt.figure(figsize=(10, 8))
        plt.contourf(omega_grid, tauomega_grid, averaged_density_grid, levels=100, cmap='viridis')
        plt.colorbar(label='Averaged Density')
        plt.xlabel('OmegaVec')
        plt.ylabel('TauOmegaVec')
        plt.title(f'2D Spectra Plot of Averaged {file_path}')

        # Remove the csv if it has it in the name and replace with png
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)

        ##############
        from scipy.signal import find_peaks
        from sklearn.preprocessing import MinMaxScaler

        # Adjusting peak analysis to be sensitive to both positive and negative values
        # Here, we'll find peaks on the absolute values of the max density along omega
        peaks_positive, _ = find_peaks(averaged_density_grid.max(axis=0))  # positive peaks
        peaks_negative, _ = find_peaks(-averaged_density_grid.max(axis=0))  # negative peaks

        # Dynamic range adjustment for the contour plots
        # Using percentiles to set color scale limits
        percentile_5 = np.percentile(averaged_density_grid, 5)
        percentile_95 = np.percentile(averaged_density_grid, 95)

        # Symmetrical color map for first derivative
        first_derivative = np.gradient(averaged_density_grid, axis=0)  # recalculating first derivative
        derivative_min = np.percentile(first_derivative, 5)
        derivative_max = np.percentile(first_derivative, 95)
        derivative_abs_max = max(abs(derivative_min), abs(derivative_max))

        # Plotting with adjustments
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))

        # Adjusted Peak Analysis Plot
        axes[0, 0].plot(omega_values, averaged_density_grid.max(axis=0))
        axes[0, 0].plot(omega_values[peaks_positive], averaged_density_grid.max(axis=0)[peaks_positive], "x",
                        color='blue')
        axes[0, 0].plot(omega_values[peaks_negative], averaged_density_grid.max(axis=0)[peaks_negative], "x",
                        color='red')
        axes[0, 0].set_title("Adjusted Peak Analysis")

        # Adjusted Normalized Density
        axes[0, 1].contourf(omega_grid, tauomega_grid, averaged_density_grid, levels=100, cmap='viridis',
                            vmin=percentile_5, vmax=percentile_95)
        axes[0, 1].set_title("Adjusted Normalized Density")

        # Adjusted First Derivative
        axes[1, 0].contourf(omega_grid, tauomega_grid, first_derivative, levels=100, cmap='seismic',
                            vmin=-derivative_abs_max, vmax=derivative_abs_max)
        axes[1, 0].set_title("Adjusted First Derivative")

        # Integrated Density
        scaler = MinMaxScaler()
        normalized_density = scaler.fit_transform(averaged_density_grid)
        integrated_density = np.trapz(normalized_density, axis=0)  # integrating along tauomega
        # The plot for integrated density remains unchanged
        axes[1, 1].plot(omega_values, integrated_density)
        axes[1, 1].set_title("Integrated Density")

        plt.tight_layout()
        # plt.show()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)

        from scipy.signal import savgol_filter

        # 1. Rethinking Peak Analysis
        # Applying a more sophisticated peak detection method
        filtered_density_max = savgol_filter(averaged_density_grid.max(axis=0), window_length=51, polyorder=3)
        peaks_refined, _ = find_peaks(filtered_density_max, prominence=0.1)  # using prominence as a criterion

        # 2. Reconsidering Normalization and Baseline Correction
        # Applying a baseline correction using a polynomial fit (for example)
        coefficients = np.polyfit(omega_values, averaged_density_grid.min(axis=0), deg=5)
        baseline_poly = np.polyval(coefficients, omega_values)
        corrected_density_poly = averaged_density_grid - baseline_poly[:, None]

        # Normalization after baseline correction
        normalized_density_poly = scaler.fit_transform(corrected_density_poly)

        # 3. First Derivative Analysis
        # Applying smoothing before differentiation
        smoothed_density = savgol_filter(normalized_density_poly, window_length=51, polyorder=3, axis=0)
        first_derivative_smoothed = np.gradient(smoothed_density, axis=0)

        # 4. Contour Plot Adjustments
        # Using standard deviations for dynamic scaling
        std_dev = np.std(smoothed_density)
        contour_min = np.mean(smoothed_density) - 2 * std_dev
        contour_max = np.mean(smoothed_density) + 2 * std_dev

        # 5. Integration of Density
        # Integration considering the entire range
        integrated_density_smoothed = np.trapz(smoothed_density, axis=0)

        # Plotting with revised strategies
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))

        # Refined Peak Analysis Plot
        axes[0, 0].plot(omega_values, filtered_density_max)
        axes[0, 0].plot(omega_values[peaks_refined], filtered_density_max[peaks_refined], "x")
        axes[0, 0].set_title("Refined Peak Analysis")

        # Normalized Density with Baseline Correction
        axes[0, 1].contourf(omega_grid, tauomega_grid, normalized_density_poly, levels=100, cmap='viridis',
                            vmin=contour_min, vmax=contour_max)
        axes[0, 1].set_title("Normalized Density with Baseline Correction")

        # First Derivative (Smoothed)
        axes[1, 0].contourf(omega_grid, tauomega_grid, first_derivative_smoothed, levels=100, cmap='seismic',
                            vmin=-std_dev, vmax=std_dev)
        axes[1, 0].set_title("First Derivative (Smoothed)")

        # Integrated Density (Smoothed)
        axes[1, 1].plot(omega_values, integrated_density_smoothed)
        axes[1, 1].set_title("Integrated Density (Smoothed)")

        plt.tight_layout()
        # plt.show()
        output_file = file_path.replace('.csv', '.png')
        output_file = output_file if output_file != file_path else f'{file_path}.png'
        plt.savefig(output_file)


def plot_ft_dipoles_v_time(study_directory, experiment_directory, time_delays, plot_all=False):
    # "number_of_pulses", "central_time_1", "carrier_frequency_1", "fwhm_1", "carrier_envelope_phase_1",
    #   "intensity_1", "amplitude_1", "period_1", "central_time_...", "carrier_frequency_...", "fwhm_...",
    #   "carrier_envelope_phase_...", "intensity_...", "amplitude_...", "period_...", "iOmega", "OmegaVec", "FTDipoleX_Re",
    #   "FTDipoleX_Im", "FTDipoleY_Re", "FTDipoleY_Im", "FTDipoleZ_Re", "FTDipoleZ_Im"

    # set pm3d
    # set pm3d interpolate 8,8
    # set view map
    # set hidden3d
    # set ytics font "Arial, 120" offset -1,0,0
    # set xtics font "Arial, 120" offset 1,-7,0
    # set xlabel "Time delay (a.u)"      font "Arial, 156" offset  0,-20.0,0
    # set ylabel "Emission energy (a.u)" font "Arial, 156" offset  -34.0,0,0
    # Average
    # splot [][0.2:0.48] 'sim1/Dipole/DipoleFT_ALL' u 9:17:((($18**2+$19**2+$20**2+$21**2+$22**2+$23**2)/3)**(0.5)) w l notitle

    # Lets plot the Spectra FT Dipolar Reponse vs Time (t)
    # Well use "central_time_2"(PROBE),"OmegaVec","FTDipoleX_Re","FTDipoleX_Im","FTDipoleY_Re","FTDipoleY_Im","FTDipoleZ_Re","FTDipoleZ_Im"
    data = pd.read_csv(f'{study_directory}/{experiment_directory}/Dipole/DipoleFT_ALL.csv')
    logger.debug(data.head())
    logger.debug(f'Length of data: {len(data)}')

    import matplotlib.pyplot as plt

    # Plotting the real and imaginary parts of the dipole components
    plt.figure(figsize=(16, 12))
    # Add the Main Top Title for the plot (file_path)
    plt.suptitle(f'{study_directory}/{experiment_directory}/Dipole/DipoleFT_ALL.csv')
    plt.subplot(3, 2, 1)
    plt.plot(data['central_time_2'], data['FTDipoleX_Re'])
    plt.title('FTDipoleX_Re')
    plt.xlabel('central_time_2')
    plt.ylabel('FTDipoleX_Re')

    plt.subplot(3, 2, 2)
    plt.plot(data['central_time_2'], data['FTDipoleX_Im'])
    plt.title('FTDipoleX_Im')
    plt.xlabel('central_time_2')
    plt.ylabel('FTDipoleX_Im')

    plt.subplot(3, 2, 3)
    plt.plot(data['central_time_2'], data['FTDipoleY_Re'])
    plt.title('FTDipoleY_Re')
    plt.xlabel('central_time_2')
    plt.ylabel('FTDipoleY_Re')

    plt.subplot(3, 2, 4)
    plt.plot(data['central_time_2'], data['FTDipoleY_Im'])
    plt.title('FTDipoleY_Im')
    plt.xlabel('central_time_2')
    plt.ylabel('FTDipoleY_Im')

    plt.subplot(3, 2, 5)
    plt.plot(data['central_time_2'], data['FTDipoleZ_Re'])
    plt.title('FTDipoleZ_Re')
    plt.xlabel('central_time_2')
    plt.ylabel('FTDipoleZ_Re')

    plt.subplot(3, 2, 6)
    plt.plot(data['central_time_2'], data['FTDipoleZ_Im'])
    plt.title('FTDipoleZ_Im')
    plt.xlabel('central_time_2')
    plt.ylabel('FTDipoleZ_Im')

    plt.tight_layout()
    plt.show()
    plt.clf()

    # Calculating the average of the squared magnitudes of the FT dipole components
    average_ft_dipole = np.sqrt((data['FTDipoleX_Re'] ** 2 + data['FTDipoleX_Im'] ** 2 +
                                 data['FTDipoleY_Re'] ** 2 + data['FTDipoleY_Im'] ** 2 +
                                 data['FTDipoleZ_Re'] ** 2 + data['FTDipoleZ_Im'] ** 2) / 3)

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
    plt.show()
    plt.clf()

    # Plot the contour map
    # Plotting X=central_time_2, Y=OmegaVec, Z=((FTDipoleX_Re**2+FTDipoleX_Im**2+FTDipoleY_Re**2+FTDipoleY_Im**2+FTDipoleZ_Re**2+FTDipoleZ_Im**2)/3)**0.5

    # Preparing data for the contour plot
    central_time_2 = data['central_time_2'].unique()
    omega_vec = data['OmegaVec'].unique()

    # Creating a grid for contour plot
    ct2_grid, omega_grid = np.meshgrid(central_time_2, omega_vec, indexing='ij')

    # Calculating the Z value for the contour plot (average FT dipole magnitude)
    # Reshaping the data to match the grid
    ft_dipole_x_re_grid = data['FTDipoleX_Re'].values.reshape(ct2_grid.shape)
    ft_dipole_x_im_grid = data['FTDipoleX_Im'].values.reshape(ct2_grid.shape)
    ft_dipole_y_re_grid = data['FTDipoleY_Re'].values.reshape(ct2_grid.shape)
    ft_dipole_y_im_grid = data['FTDipoleY_Im'].values.reshape(ct2_grid.shape)
    ft_dipole_z_re_grid = data['FTDipoleZ_Re'].values.reshape(ct2_grid.shape)
    ft_dipole_z_im_grid = data['FTDipoleZ_Im'].values.reshape(ct2_grid.shape)

    average_ft_dipole_grid = np.sqrt((ft_dipole_x_re_grid ** 2 + ft_dipole_x_im_grid ** 2 +
                                      ft_dipole_y_re_grid ** 2 + ft_dipole_y_im_grid ** 2 +
                                      ft_dipole_z_re_grid ** 2 + ft_dipole_z_im_grid ** 2) / 3)

    import matplotlib.colors as colors

    # Adjusting color range to be symmetric around 0 and maximizing information display
    # Calculating the maximum absolute value in the data for symmetric color scaling
    max_abs_value = np.max(np.abs(average_ft_dipole_grid))

    # Using a diverging colormap with symmetric range around 0
    norm = colors.TwoSlopeNorm(vmin=-max_abs_value, vcenter=0, vmax=max_abs_value)

    # Plotting with the adjusted color range
    plt.figure(figsize=(12, 8))
    contour = plt.contourf(ct2_grid, omega_grid, average_ft_dipole_grid, levels=100,
                           cmap='seismic',
                           norm=norm)
    plt.colorbar(contour, label='Average FT Dipole Magnitude')
    plt.xlabel('Central Time 2')
    plt.ylabel('OmegaVec')
    plt.title('Adjusted Contour Map of Average FT Dipole Magnitude')
    plt.show()
    plt.clf()

    # Dynamic Range Adjustment using percentiles to set color scale limits
    percentile_5 = np.percentile(average_ft_dipole_grid, 5)
    percentile_95 = np.percentile(average_ft_dipole_grid, 95)

    # Adjusting the normalization to focus on the significant variations
    norm_percentile = colors.Normalize(vmin=percentile_5, vmax=percentile_95)

    # Plotting with percentile-based dynamic range adjustment
    plt.figure(figsize=(12, 8))
    contour_percentile = plt.contourf(ct2_grid, omega_grid, average_ft_dipole_grid, levels=100, cmap='viridis',
                                      norm=norm_percentile)
    plt.colorbar(contour_percentile, label='Average FT Dipole Magnitude')
    plt.xlabel('Central Time 2')
    plt.ylabel('OmegaVec')
    plt.title('Percentile Adjusted Contour Map of Average FT Dipole Magnitude')
    plt.show()
    plt.clf()

    # Applying the requested techniques for feature extraction:

    # 1. Dynamic Range Adjustment using different percentiles
    percentile_1 = np.percentile(average_ft_dipole_grid, 1)
    percentile_99 = np.percentile(average_ft_dipole_grid, 99)
    norm_1_99 = colors.Normalize(vmin=percentile_1, vmax=percentile_99)

    # 3. Smoothing and Filtering
    from scipy.ndimage import gaussian_filter
    smoothed_data = gaussian_filter(average_ft_dipole_grid, sigma=1)  # Applying a Gaussian filter

    # 4. Derivative Plots
    first_derivative = np.gradient(smoothed_data, axis=0)  # First derivative along the central_time_2 axis

    # 5. Interactive Visualization (Skipping as it requires an interactive environment like Jupyter Notebook with interactive widgets)

    # 6. Multiple Plots with Different Scales
    # Applying logarithmic scaling
    log_scaled_data = np.log1p(np.abs(average_ft_dipole_grid))  # Adding 1 to avoid log(0)

    # 7. Statistical Analysis: Clustering (Skipping as it's not typically applied in this type of visualization)

    # Plotting with the applied techniques
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))

    # Plot with Percentile Adjustment (1% to 99%)
    axes[0, 0].contourf(ct2_grid, omega_grid, average_ft_dipole_grid, levels=100, cmap='viridis', norm=norm_1_99)
    axes[0, 0].set_title('Percentile (1% - 99%) Adjusted Contour Map')
    axes[0, 0].set_xlabel('Central Time 2')
    axes[0, 0].set_ylabel('OmegaVec')

    # Plot with Smoothing (Gaussian Filter)
    axes[0, 1].contourf(ct2_grid, omega_grid, smoothed_data, levels=100, cmap='viridis')
    axes[0, 1].set_title('Smoothed Data Contour Map')
    axes[0, 1].set_xlabel('Central Time 2')
    axes[0, 1].set_ylabel('OmegaVec')

    # Plot with First Derivative
    axes[1, 0].contourf(ct2_grid, omega_grid, first_derivative, levels=100, cmap='seismic')
    axes[1, 0].set_title('First Derivative Contour Map')
    axes[1, 0].set_xlabel('Central Time 2')
    axes[1, 0].set_ylabel('OmegaVec')

    # Plot with Logarithmic Scaling
    axes[1, 1].contourf(ct2_grid, omega_grid, log_scaled_data, levels=100, cmap='viridis')
    axes[1, 1].set_title('Logarithmic Scaled Data Contour Map')
    axes[1, 1].set_xlabel('Central Time 2')
    axes[1, 1].set_ylabel('OmegaVec')

    # Skipping Statistical Analysis (Clustering) Plot

    # # Empty subplot as a placeholder for the interactive plot and clustering plot
    # axes[2, 0].axis('off')
    # axes[2, 1].axis('off')

    plt.tight_layout()
    plt.show()
    plt.clf()

    # 7. Statistical Analysis: Clustering
    # Standardizing the data for clustering
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(average_ft_dipole_grid.flatten().reshape(-1, 1))

    # Applying KMeans clustering
    kmeans = KMeans(n_clusters=6, random_state=0).fit(standardized_data)
    clustered_data = kmeans.labels_.reshape(average_ft_dipole_grid.shape)

    # Plotting the result of clustering
    plt.figure(figsize=(12, 8))
    plt.contourf(ct2_grid, omega_grid, clustered_data, cmap='viridis')
    plt.colorbar(label='Cluster ID')
    plt.xlabel('Central Time 2')
    plt.ylabel('OmegaVec')
    plt.title('KMeans Clustering of Average FT Dipole Magnitude')
    plt.show()
    plt.clf()

    import plotly.graph_objects as go

    average_ft_dipole = np.sqrt((data['FTDipoleX_Re'] ** 2 + data['FTDipoleX_Im'] ** 2 +
                                    data['FTDipoleY_Re'] ** 2 + data['FTDipoleY_Im'] ** 2 +
                                    data['FTDipoleZ_Re'] ** 2 + data['FTDipoleZ_Im'] ** 2) / 3)

    # Creating an interactive plot using Plotly
    fig = go.Figure(data=go.Contour(
        z=average_ft_dipole,
        x=data['central_time_2'],  # horizontal axis
        y=data['OmegaVec'],  # vertical axis
    ))

    # Updating layout for better visualization
    fig.update_layout(
        title='Interactive Contour Plot of Average FT Dipole Magnitude',
        xaxis_title='Central Time 2',
        yaxis_title='OmegaVec'
    )

    # Displaying the figure
    fig.show()
