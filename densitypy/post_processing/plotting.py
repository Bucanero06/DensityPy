import os
import re

from matplotlib import pyplot as plt
from matplotlib.colors import Normalize

DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING = {
    'DipoleX_Re': 'blue',
    'DipoleY_Re': 'green',
    'DipoleZ_Re': 'red',
    'DipoleX_Im': 'lightblue',
    'DipoleY_Im': 'lightgreen',
    'DipoleZ_Im': 'pink',
}
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__)
AGGREGATE_FILES = [  # fixme
    ('all', "DipoleFT_ALL.csv"),
    ('ww', "DipolePP_ww.csv")
]


class PlotDipolesSimulation:  # todo change x axis to be time where appropriate rather than the index of timestep
    def __init__(self, study_directory, experiment_directory, dpi=300):
        self.palette = [
            [0.00, 0.00, 0.00, 0.50],
            [0.10, 0.00, 0.50, 1.00],
            [0.17, 0.00, 1.00, 1.00],
            [0.25, 0.00, 1.00, 0.50],
            [0.34, 0.50, 1.00, 0.00],
            [0.46, 1.00, 1.00, 0.00],
            [0.61, 1.00, 0.50, 0.00],
            [1.00, 0.50, 0.00, 0.00]
        ]
        self.dipoles_dir = f"{study_directory}/{experiment_directory}/Dipole"
        self.atomic_charge_dir = f"{study_directory}/{experiment_directory}/AtomicCharge"
        self.pulse_dir = f"{study_directory}/{experiment_directory}/Pulse"
        self.dpi = dpi
        self.height = 1080
        self.width = 1920

    def load_df(self, filename_or_df):
        if isinstance(filename_or_df, str):
            return pd.read_csv(filename_or_df, sep=',', index_col=0)
        elif isinstance(filename_or_df, pd.DataFrame):
            return filename_or_df
        else:
            raise ValueError('filename_or_df must be a string or a pandas DataFrame')

    def init_plot(self):
        fig, ax1 = plt.subplots(figsize=(self.width / self.dpi, self.height / self.dpi), dpi=self.dpi)
        return fig, ax1

    def plot_dataframe(self, cols, filename_or_df, output_filename):
        logger.info(f'Plotting {output_filename} using columns: {cols}')
        df = self.load_df(filename_or_df)
        fig, ax1 = self.init_plot()

        for col in cols:
            if df[col].std() != 0:  # Don't plot if standard deviation is zero
                ax1.plot(df.index, df[col], color=DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING[col], label=col)

        plt.xlabel('Time (a.u.)')
        plt.ylabel('Dipole (a.u.)')
        plt.legend()
        plt.savefig(output_filename, dpi=self.dpi)
        plt.close()

    def plot_twinx_dipoles(self, cols_re, cols_im, filename_or_df, output_filename):
        logger.info(f'Plotting {output_filename} using columns: {cols_re} and {cols_im}')

        if isinstance(filename_or_df, str):
            df = pd.read_csv(filename_or_df, sep=',', index_col=0)
        elif isinstance(filename_or_df, pd.DataFrame):
            df = filename_or_df
        else:
            raise ValueError('filename_or_df must be a string or a pandas DataFrame')

        plt.figure(figsize=(self.width / self.dpi, self.height / self.dpi), dpi=self.dpi)

        fig, ax1 = plt.subplots(figsize=(self.width / self.dpi, self.height / self.dpi), dpi=self.dpi)
        ax2 = ax1.twinx()

        for col in cols_re:
            if df[col].std() != 0:  # Don't plot if standard deviation is zero
                ax1.plot(df.index, df[col], color=DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING[col], label=col)

        for col in cols_im:
            if df[col].std() != 0:  # Don't plot if standard deviation is zero
                ax2.plot(df.index, df[col], color=DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING[col], label=col)

        ax1.set_xlabel('Time (a.u.)')
        ax1.set_ylabel('Dipole (Real, a.u.)')
        ax2.set_ylabel('Dipole (Imaginary, a.u.)')

        # Align zeros
        ax1.set_ylim(bottom=-max(abs(ax1.get_ylim()[0]), ax1.get_ylim()[1]),
                     top=max(abs(ax1.get_ylim()[0]), ax1.get_ylim()[1]))
        ax2.set_ylim(bottom=-max(abs(ax2.get_ylim()[0]), ax2.get_ylim()[1]),
                     top=max(abs(ax2.get_ylim()[0]), ax2.get_ylim()[1]))

        fig.legend()
        plt.savefig(output_filename, dpi=self.dpi)
        plt.close()

    def plot_interactive_dipoles(self, cols_re, cols_im, filename_or_df, output_filename):
        logger.info(f'Plotting {output_filename} using columns: {cols_re} and {cols_im}')
        if isinstance(filename_or_df, str):
            df = pd.read_csv(f'{self.dipoles_dir}/DipolePP0.0.csv', sep=',', index_col=0)
        elif isinstance(filename_or_df, pd.DataFrame):
            df = filename_or_df
        else:
            raise ValueError('filename_or_df must be a string or a pandas DataFrame')

        fig = make_subplots(specs=[[{"secondary_y": True}]])

        min_re, max_re, min_im, max_im = float('inf'), float('-inf'), float('inf'), float('-inf')
        for col in cols_re:
            if df[col].std() != 0:  # Don't plot if standard deviation is zero
                fig.add_trace(
                    go.Scatter(x=df.index, y=df[col], name=col, line_color=DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING[col]),
                    secondary_y=False)
                min_re, max_re = min(min_re, df[col].min()), max(max_re, df[col].max())

        for col in cols_im:
            if df[col].std() != 0:  # Don't plot if standard deviation is zero
                fig.add_trace(
                    go.Scatter(x=df.index, y=df[col], name=col, line_color=DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING[col]),
                    secondary_y=True)
                min_im, max_im = min(min_im, df[col].min()), max(max_im, df[col].max())

        # Align zeros
        max_range_re = max(abs(min_re), abs(max_re))
        max_range_im = max(abs(min_im), abs(max_im))

        fig.update_layout(height=1080, width=1920, title_text="Dipoles")
        fig.update_xaxes(title_text="Time (a.u.)")
        fig.update_yaxes(title_text="Dipole (Real, a.u.)", range=[-max_range_re, max_range_re], secondary_y=False)
        fig.update_yaxes(title_text="Dipole (Imaginary, a.u.)", range=[-max_range_im, max_range_im], secondary_y=True)

        fig.write_html(output_filename)

    def plot_pp_dipole_file(self, label, file_path):
        logger.info('Plotting all dipoles')
        df = self.load_df(f'{self.dipoles_dir}/{file_path}')

        self.plot_dataframe(['DipoleX_Re', 'DipoleY_Re', 'DipoleZ_Re'], df, f'{self.dipoles_dir}/Dipole_Re_{label}.png')
        self.plot_dataframe(['DipoleX_Im', 'DipoleY_Im', 'DipoleZ_Im'], df, f'{self.dipoles_dir}/Dipole_Im_{label}.png')
        self.plot_dataframe(['DipoleX_Re', 'DipoleX_Im', 'DipoleY_Re', 'DipoleY_Im', 'DipoleZ_Re', 'DipoleZ_Im'], df,
                            f'{self.dipoles_dir}/Dipole_{label}.png')
        self.plot_twinx_dipoles(['DipoleX_Re', 'DipoleY_Re', 'DipoleZ_Re'], ['DipoleX_Im', 'DipoleY_Im', 'DipoleZ_Im'],
                                df, f'{self.dipoles_dir}/Dipole_Re_Im_{label}.png')
        self.plot_interactive_dipoles(['DipoleX_Re', 'DipoleY_Re', 'DipoleZ_Re'],
                                      ['DipoleX_Im', 'DipoleY_Im', 'DipoleZ_Im'], df,
                                      f'{self.dipoles_dir}/Dipole_Re_Im_{label}.html')

    @staticmethod
    def build_sim_tagv(time_number):
        return f'{time_number:.1f}'

    @property
    def get_every_pp_csv_file_names(self):
        pattern = re.compile(r'DipolePP(-?\d+\.?\d*)\.csv')
        return [(float(pattern.match(f).group(1)), f) for f in os.listdir(self.dipoles_dir) if
                pattern.match(f) is not None]

    @property
    def get_every_ftpp_csv_file_names(self):
        pattern = re.compile(r'DipoleFTPP(-?\d+\.?\d*)\.csv')
        return [(float(pattern.match(f).group(1)), f) for f in os.listdir(self.dipoles_dir) if
                pattern.match(f) is not None]

    @property
    def get_every_aggregate_csv_file_names(self):
        return AGGREGATE_FILES

    @staticmethod
    def interpolate_colormap(start_color, end_color, n):
        import numpy as np
        import matplotlib.colors as mcolors

        start_rgb = np.array(mcolors.to_rgb(start_color))
        end_rgb = np.array(mcolors.to_rgb(end_color))

        return [start_rgb + i * (end_rgb - start_rgb) for i in np.linspace(0, 1, n)]

    def plot_ftpp_all(self):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        aggregate_files = self.get_every_aggregate_csv_file_names

        # Filter files with 'all'
        all_file = [filename for keyword, filename in aggregate_files if keyword == 'all'][0]
        # Assuming each file has similar data structure (e.g., same columns)
        # ['number_of_pulses', 'central_time_1  ', 'carrier_frequency', 'fwhm',
        #  'carrier_envelope_phase', 'intensity', 'amplitude', 'period',
        #  'central_time_2  ', 'carrier_frequency.1', 'fwhm.1',
        #  'carrier_envelope_phase.1', 'intensity.1', 'amplitude.1', 'period.1',
        #  'iOmega', 'Omega', 'DipoleX_Re', 'DipoleX_Im', 'DipoleY_Re',
        #  'DipoleY_Im', 'DipoleZ_Re', 'DipoleZ_Im']
        df = self.load_df(f'{self.dipoles_dir}/{all_file}')


        # Compute Z_Value column
        cols = ['DipoleX_Re', # 18
                # 'DipoleX_Im', # 19
                'DipoleY_Re', # 20
                # 'DipoleY_Im', # 21
                'DipoleZ_Re', # 22
                # 'DipoleZ_Im' # 23
                ]

        # splot [][0.2:0.48] 'sim1/Dipole/DipoleFT_ALL' u 9:17:((($18**2+$19**2+$20**2+$21**2+$22**2+$23**2)/3)**(0.5)) w l notitle



        print(df.columns)
        exit()
        X = df["central_time"] # fixme found 2 of this column in the DipoleFT_ALL.csv file
        Y = df["Omega"]
        df['averaged_density'] = (sum(df[col] ** 2 for col in cols) / 3) ** 0.5
        Z = df['averaged_density']

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Create a 3D plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Scatter plot
        ax.scatter(df['central_time'], df['Omega'], df['averaged_density'], c=df['averaged_density'], cmap='viridis')

        # Setting labels
        ax.set_xlabel('Central Time')
        ax.set_ylabel('Omega')
        ax.set_zlabel('Averaged Density')
        ax.set_title('3D Plot of Central Time, Omega and Averaged Density')

        plt.show()



