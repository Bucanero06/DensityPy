import pandas as pd

DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING = {
    'DipoleX_Re': 'blue',
    'DipoleY_Re': 'green',
    'DipoleZ_Re': 'red',
    'DipoleX_Im': 'lightblue',
    'DipoleY_Im': 'lightgreen',
    'DipoleZ_Im': 'pink',
}


def plot_pulses(study_directory, experiment_directory, time_delays):
    # pulsePP{delay}
    # Time: The time at which the following measurements are taken, typically in atomic units (a.u.).
    # Ax: The x-component of the vector potential. This is calculated as the real part of the difference between zw_m1 and zw_p1, normalized by sqrt(2).
    # Ay: The y-component of the vector potential. Computed as the imaginary part of the sum of zw_m1 and zw_p1, normalized by sqrt(2).
    # Az: The z-component of the vector potential, calculated as the real part of zw_p0.
    # Real_zw_m1: The real part of the vector potential zw_m1, which represents one of the spherical components of the vector potential at the specified time.
    # Imag_zw_m1: The imaginary part of zw_m1.
    # Real_zw_p0: The real part of zw_p0, another spherical component of the vector potential.
    # Imag_zw_p0: The imaginary part of zw_p0.
    # Real_zw_p1: The real part of zw_p1, the third spherical component of the vector potential.
    # Imag_zw_p1: The imaginary part of zw_p1

    zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
    one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
    two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
    max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))

    # Lets also do the same for the minimum time delay and the "thirds" splits as before (remember to get as close as possible to a number
    min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
    one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
    two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

    length_of_data_to_match = None

    for time_delay in ['XUV'] + time_delays:
        file_path = f'{study_directory}/{experiment_directory}/Pulses/pulsePP{time_delay}' if time_delay != 'XUV' else f'{study_directory}/{experiment_directory}/Pulses/pulseXUV'

        print(f'Plotting {file_path}')
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
            print(f'File {file_path} not found')
            continue

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'

        # print all columns in pandas
        pd.set_option('display.max_columns', None)
        print(data.head())
        print(len(data))

        # Plotting the pulse data (pulsePP0.0)
        from matplotlib import pyplot as plt

        plt.figure(figsize=(12, 6))
        for col in data.columns[1:]:
            plt.plot(data['Time'], data[col], label=f'Column {col}')
        # plt.plot(data['Time'], data['Ax'], label=f'Column Ax')
        plt.xlabel('Time (arbitrary units)')
        plt.ylabel('Value')
        if time_delay != 'XUV':
            plt.title(f'Pulse Characteristics (pulsePP{time_delay})')
        else:
            plt.title(f'Pulse Characteristics (pulseXUV)')
        plt.legend()
        plt.grid(True)
        plt.show()

    plt.clf()



def plot_ft_pulses(study_directory, experiment_directory, time_delays):
    # FTpulsePP{delay}
    # Freq: The frequency at which the Fourier transform is computed.
    # FTAx: The x-component of the Fourier-transformed vector potential.
    # FTAy: The y-component of the Fourier-transformed vector potential.
    # FTAz: The z-component of the Fourier-transformed vector potential.
    # FT_Aminus1_Real and FT_Aminus1_Imag: The real and imaginary parts of the Fourier transform for the mu = -1 component.
    # FT_A0_Real and FT_A0_Imag: The real and imaginary parts of the Fourier transform for the mu = 0 component.
    # FT_Aplus1_Real and FT_Aplus1_Imag: The real and imaginary parts of the Fourier transform for the mu = +1 component.

    zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
    one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
    two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
    max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))

    # Lets also do the same for the minimum time delay and the "thirds" splits as before (remember to get as close as possible to a number
    min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
    one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
    two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

    for time_delay in ['XUV'] + time_delays:
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
            print(f'File {file_path} not found')
            continue

        # print all columns in pandas
        pd.set_option('display.max_columns', None)
        print(data.head())
        print(len(data))

        # Plotting the pulse data (pulsePP0.0)
        from matplotlib import pyplot as plt

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
        plt.show()

    # Clear the plot
    plt.clf()



def plot_dipole_response_vs_time(study_directory, experiment_directory, time_delays):
    # Lets plot the Dipolar Reponse vs Time (t)
    # "itime","Time","DipoleX_Re","DipoleX_Im","DipoleY_Re","DipoleY_Im","DipoleZ_Re","DipoleZ_Im"

    # Use the time_delays to get the files to plot
    # Time delays to plot = [0, 1/3 max, 2/3 max, max] if not present then use the nearest. Do the same for min
    zero_time_delay = min(time_delays, key=lambda x: abs(x - 0))
    one_third_max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays) / 3))
    two_third_max_time_delay = min(time_delays, key=lambda x: abs(x - 2 * max(time_delays) / 3))
    max_time_delay = min(time_delays, key=lambda x: abs(x - max(time_delays)))

    # Lets also do the same for the minimum time delay and the "thirds" splits as before (remember to get as close as possible to a number
    min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays)))
    one_third_min_time_delay = min(time_delays, key=lambda x: abs(x - min(time_delays) / 3))
    two_third_min_time_delay = min(time_delays, key=lambda x: abs(x - 2 * min(time_delays) / 3))

    DIPOLE_PP_FILES_TO_PLOT_PATHS = [
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipoleXUV.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{min_time_delay}.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{one_third_min_time_delay}.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{two_third_min_time_delay}.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{zero_time_delay}.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{one_third_max_time_delay}.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{two_third_max_time_delay}.csv',
        f'{study_directory}/{experiment_directory}/Dipole'
        f'/DipolePP{max_time_delay}.csv',
    ]

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
    import matplotlib.pyplot as plt

    for file_path in DIPOLE_PP_FILES_TO_PLOT_PATHS:
        print(f'Plotting {file_path}')
        try:
            data = pd.read_csv(file_path)
        except FileNotFoundError:
            print(f'File {file_path} not found')
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
        plt.title('Dipole X Real')
        plt.xlabel('Time')
        plt.ylabel('Dipole X Re')

        plt.subplot(3, 2, 2)
        plt.plot(time, dipole_x_im)
        plt.title('Dipole X Imaginary')
        plt.xlabel('Time')
        plt.ylabel('Dipole X Im')

        plt.subplot(3, 2, 3)
        plt.plot(time, dipole_y_re)
        plt.title('Dipole Y Real')
        plt.xlabel('Time')
        plt.ylabel('Dipole Y Re')

        plt.subplot(3, 2, 4)
        plt.plot(time, dipole_y_im)
        plt.title('Dipole Y Imaginary')
        plt.xlabel('Time')
        plt.ylabel('Dipole Y Im')

        plt.subplot(3, 2, 5)
        plt.plot(time, dipole_z_re)
        plt.title('Dipole Z Real')
        plt.xlabel('Time')
        plt.ylabel('Dipole Z Re')

        plt.subplot(3, 2, 6)
        plt.plot(time, dipole_z_im)
        plt.title('Dipole Z Imaginary')
        plt.xlabel('Time')
        plt.ylabel('Dipole Z Im')

        plt.tight_layout()
        plt.show()

        # Lets Get the Dipole Analytical Metrics
        # correlation = data.corr()

    plt.clf()


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
        print(f'Plotting {file_path}')
        data = pd.read_csv(file_path)

        if length_of_data_to_match is None:
            length_of_data_to_match = len(data)
        else:
            assert length_of_data_to_match == len(data), f'Length of {file_path} is not the same as the previous file'
            assert data.columns.all() in INDECES_COL_NAMES + FEATURES_COL_NAMES, f'Column names of {file_path} is not the same as the previous file'

        # Displaying the first few rows of the file to understand its structure
        print(data.head())

        # Lets check length of the data
        print(len(data))

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
        plt.show()

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
        plt.show()

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
        plt.show()
