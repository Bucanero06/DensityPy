# Autodocumentation for Fortran using ford

import os

from densitypy.project_utils.command_execution import execute_command
from densitypy.project_utils.file_directory_ops import change_directory_manager
from densitypy.project_utils.string_dict_utils import remove_spaces_and_set_case


def main(project_name, project_description, source_dirs, documentation_dir,exclude_dirs=None,graphs=True, delete_old_files=False):
    if not exclude_dirs:
        exclude_dirs = []

    if delete_old_files and os.path.exists(documentation_dir):
        from densitypy.project_utils.file_directory_ops import delete_files_or_directories
        delete_files_or_directories(documentation_dir)

    if not os.path.exists(documentation_dir):
        os.makedirs(documentation_dir)

    # Create My Project .md file
    new_readme_file_path = os.path.join(documentation_dir,
                                        f'{remove_spaces_and_set_case(project_name, case="u")}.md')

    with open(new_readme_file_path, 'a') as file:
        file.write('---\n')
        file.write(f'project: {project_name}\n')
        file.write('---\n\n')
        file.write(project_description + '\n\n')

    # Run ford
    with change_directory_manager(documentation_dir):
        command_to_run = f"ford {new_readme_file_path}"
        for src_dir in source_dirs:
            command_to_run += f" --src_dir {src_dir}"
        for exclude_dir in exclude_dirs:
            command_to_run += f" --exclude_dir {exclude_dir}"
        if graphs:
            command_to_run += " --graph"
        execute_command(command_to_run)


if __name__ == '__main__':
    main(project_name='Charge Migration',
         project_description=('Fortran code for the calculation of the charge migration in molecules. '
                              'Interfaced with Python\'s DensityPy library.'),
         source_dirs=[
             '../../densityfort/ChargeMigration',
                      # '../../densityfort/ChargeMigrationFT',
                      # '../../densityfort/SpectrumReconstruction'
                      ],
         documentation_dir='/home/ruben/PycharmProjects/DensityPy/Documentation/fortrandocs',
         exclude_dirs=['../../densityfort/ChargeMigration/src/future_reconstruction_using_ci_work',
                       '../../densityfort/ChargeMigrationFT/src/future_reconstruction_using_ci_work'],
         graphs=False,
         delete_old_files=True)
