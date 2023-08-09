from matplotlib import pyplot as plt

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


class PlotSimulation:  # todo change x axis to be time where appropriate rather than the index of timestep
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
            return pd.read_csv(f'{self.dipoles_dir}/DipolePP0.0.csv', sep=',', index_col=0)
        elif isinstance(filename_or_df, pd.DataFrame):
            return filename_or_df
        else:
            raise ValueError('filename_or_df must be a string or a pandas DataFrame')

    def init_plot(self):
        fig, ax1 = plt.subplots(figsize=(self.width / self.dpi, self.height / self.dpi), dpi=self.dpi)
        return fig, ax1

    def plot_dipoles(self, cols, filename_or_df, output_filename):
        logger.info(f'Plotting {output_filename} using columns: {cols}')
        df = self.load_df(filename_or_df)
        fig, ax1 = self.init_plot()

        for col in cols:
            if df[col].std() != 0:  # Don't plot if standard deviation is zero
                ax1.plot(df.index, df[col], color=DEFAULT_DIPOLE_COLUMNS_COLOR_MAPPING[col], label=col)

        plt.xlabel('Time (a.u.)')
        plt.ylabel('Dipole (a.u.)')
        plt.legend()
        plt.savefig(f'{self.dipoles_dir}/{output_filename}', dpi=self.dpi)
        plt.close()

    def plot_twinx_dipoles(self, cols_re, cols_im, filename_or_df, output_filename):
        logger.info(f'Plotting {output_filename} using columns: {cols_re} and {cols_im}')

        if isinstance(filename_or_df, str):
            df = pd.read_csv(f'{self.dipoles_dir}/DipolePP0.0.csv', sep=',', index_col=0)
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
        plt.savefig(f'{self.dipoles_dir}/{output_filename}', dpi=self.dpi)
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

        fig.write_html(f'{self.dipoles_dir}/{output_filename}')

    def plot_all_dipoles(self):
        logger.info('Plotting all dipoles')
        df = self.load_df(f'{self.dipoles_dir}/DipolePP0.0.csv')
        self.plot_dipoles(['DipoleX_Re', 'DipoleY_Re', 'DipoleZ_Re'], df, 'Dipole_Re.png')
        self.plot_dipoles(['DipoleX_Im', 'DipoleY_Im', 'DipoleZ_Im'], df, 'Dipole_Im.png')
        self.plot_dipoles(['DipoleX_Re', 'DipoleX_Im', 'DipoleY_Re', 'DipoleY_Im', 'DipoleZ_Re', 'DipoleZ_Im'],
                          df, 'Dipole.png')
        self.plot_twinx_dipoles(['DipoleX_Re', 'DipoleY_Re', 'DipoleZ_Re'], ['DipoleX_Im', 'DipoleY_Im', 'DipoleZ_Im'],
                                df, 'Dipole_Re_Im.png')
        self.plot_interactive_dipoles(['DipoleX_Re', 'DipoleY_Re', 'DipoleZ_Re'],
                                      ['DipoleX_Im', 'DipoleY_Im', 'DipoleZ_Im'], df, 'Dipole_Re_Im.html')
