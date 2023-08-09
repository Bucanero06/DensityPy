from io import StringIO
from pathlib import Path

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from typing import Iterable, List, Tuple

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure

from densitypy.project_utils.logger import logger


class Component:
    NAMES = 'XYZ'
    COLORMAP: Colormap = sns.color_palette('coolwarm', as_cmap=True)

    def __init__(self, section: str):
        *_, component, body = section.split(maxsplit=3)
        self.index = int(component)
        self.name = self.NAMES[self.index - 1]
        self.df = self._parse_component_data(body)

    @classmethod
    def _parse_component_data(cls, body: str) -> pd.DataFrame:
        colsets = [
            cls._str_to_df(colset)
            for colset in cls._generate_colset_strings(body)
        ]

        df = pd.concat(
            [colsets[0]] + [df.drop('STATE', axis=1) for df in colsets[1:]],
            axis=1,
        )
        df.set_index('STATE', inplace=True)

        return df

    @staticmethod
    def _generate_colset_strings(body: str) -> Iterable[str]:
        first_col = 'STATE'
        start = body.find(first_col)
        while start != -1:
            end = body.find(first_col, start + len(first_col))
            yield body[start:end]
            start = end

    @staticmethod
    def _str_to_df(body: str) -> pd.DataFrame:
        with StringIO(body) as f:
            return pd.read_csv(f, sep=r'\s+')

    def __str__(self):
        return self.name

    @property
    def without_state(self) -> pd.DataFrame:
        if 'STATE' in self.df.columns:
            return self.df.drop('STATE', axis=1)
        return self.df
    def write_csv(self, path: Path, bare: bool = True):
        if bare:
            logger.info(f"Writing {self.name} to {path}")
            self.without_state.to_csv(path, index=True, header=True)
        else:
            self.df.to_csv(path, index=True)


    def mu_heatmap(self, ax: Axes, vmax: float = 1, vmin: float = -1, center: float = 0) -> Axes:
        ax.set_title(self.name)

        return sns.heatmap(
            self.without_state,
            ax=ax,
            cmap=self.COLORMAP,
            annot=True,
            vmax=vmax,
            vmin=vmin,
            center=center,
            square=True,
            linewidths=.5,
            cbar_kws={'shrink': .5},
        )


class DipolesLogParser:
    def __init__(self, log_path: str):
        try:
            self.components: Tuple[Component] = tuple(self._parse_log_file(log_path))
        except FileNotFoundError as fnf_error:
            logger.info(f"No such file or directory: {log_path}")
            raise fnf_error

    @staticmethod
    def _parse_log_file(log_path: str):
        with open(log_path) as f:
            body = f.read()

        start_i = body.find('PROPERTY: MLTPL  1 ')
        end_i = body.find('PROPERTY: MLTPL  2 ', start_i)
        body = body[start_i: end_i]
        sections = body.split('PROPERTY: MLTPL')[1:]

        for section in sections:
            yield Component(section)

    def write_csvs(
        self,
        directory: str = '.',
        suffix: str = '_DIPOLE.csv',
        bare: bool = True,
    ):
        directory = Path(directory)
        for component in self.components:
            component.write_csv(directory / (component.name + suffix), bare)

    def get_mu_heatmaps(self) -> Figure:
        fig, ((z, x), (y, unused)) = plt.subplots(
            2, 2, sharex='col', sharey='row',
            figsize=(10, 10),
        )
        unused.remove()

        for ax, component in zip((x, y, z), self.components):
            component.mu_heatmap(ax)

        return fig





if __name__ == '__main__':
    def parse_dipoles_test():
        logger.warning("This is a test")
        working_dir = '/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/NMA_output'
        project = DipolesLogParser('/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/NMA_output/NMA.log')
        project.write_csvs(directory=working_dir)
        fig = project.get_mu_heatmaps()
        fig.show()
        fig.savefig(f'{working_dir}/dipole-heatmap.png')


    parse_dipoles_test()
