# Autodocumentation for Fortran using ford

import os

from densitypy.project_utils.command_execution import execute_command
from densitypy.project_utils.file_directory_ops import change_directory_manager
from densitypy.project_utils.string_dict_utils import remove_spaces_and_set_case


def validate_directory_and_get_full_path(path):
    path = path.strip()
    if not os.path.exists(path):
        raise ValueError(f'{path} does not exist.')
    if not os.path.isdir(path):
        raise ValueError(f'{path} is not a directory.')
    return os.path.abspath(path)


def main(project_name, project_description, source_dirs, documentation_dir, exclude_dirs=None, graphs=True,
         delete_old_files=False):
    documentation_dir = validate_directory_and_get_full_path(documentation_dir)
    if not source_dirs:
        source_dirs = ['./']

    if not exclude_dirs:
        exclude_dirs = []
    if isinstance(source_dirs, str):
        source_dirs = source_dirs.split(',')
        source_dirs = [validate_directory_and_get_full_path(source_dir) for source_dir in source_dirs]
    if isinstance(exclude_dirs, str):
        exclude_dirs = exclude_dirs.split(',')
        exclude_dirs = [validate_directory_and_get_full_path(exclude_dir) for exclude_dir in exclude_dirs]

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
    import click

    @click.command()
    @click.option('--project_name', default='My Project', help='Name of the project.')
    @click.option('--project_description', default='My Project Description', help='Description of the project.')
    @click.option('--source_dirs', default=['./'], help='Source directories to be documented.')
    @click.option('--documentation_dir', default='./docs', help='Directory where the documentation will be stored.')
    @click.option('--exclude_dirs', default=[], help='Directories to be excluded from the documentation.')
    @click.option('--graphs', default=True, help='Whether to include graphs in the documentation.')
    @click.option('--delete_old_files', default=False,
                  help='Whether to delete old files in the documentation directory.')
    def main_cli(project_name, project_description, source_dirs, documentation_dir, exclude_dirs, graphs,
                 delete_old_files):
        main(project_name=project_name, project_description=project_description, source_dirs=source_dirs,
             documentation_dir=documentation_dir, exclude_dirs=exclude_dirs, graphs=graphs,
             delete_old_files=delete_old_files)


    main_cli()
