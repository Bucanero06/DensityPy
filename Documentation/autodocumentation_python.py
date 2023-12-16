import logging
import os
import subprocess

from densitypy.project_utils.command_execution import execute_command
from densitypy.project_utils.file_directory_ops import change_directory_manager


def initialize_sphinx(source_dir, project_name, author_name, version='0.0.1'):
    # Run sphinx-quickstart automatically
    with change_directory_manager(source_dir):
        execute_command(f"sphinx-quickstart --quiet "
                        f"--project '{project_name}' "
                        f"--author '{author_name}' "
                        f"--language en "
                        f"-v {version} --release {version} "
                        f"--suffix .rst --master index --makefile --batchfile --sep --dot _"
                        )
        # subprocess.run(['sphinx-quickstart', '--quiet', '--project', str(project_name),  # '--author', str(author_name),
        #                 # '-v', str(version), '--release', str(version),
        #                 '--language', 'en', '--suffix', '.rst',
        #                 '--master', 'index', '--makefile', '--batchfile', '--sep', '--dot', '_'], cwd=source_dir)


def find_python_modules(start_path, ignore_folders=None):
    if ignore_folders is None:
        ignore_folders = []
    exclude_paths = [os.path.join(start_path, p) for p in ignore_folders]

    modules = []
    for root, dirs, files in os.walk(start_path):
        # Exclude the specified directories and their subdirectories
        dirs[:] = [d for d in dirs if os.path.join(root, d) not in exclude_paths]

        for file in files:
            if file.endswith('.py') and not file.startswith('__'):
                module_path = os.path.join(root, file)
                module_name = os.path.relpath(module_path, start_path).replace(os.path.sep, '.').rstrip('.py')
                modules.append(module_name)

    # FIXME HOTFIX for plotdensity file, not sure where the error is coming from but is only this module
    for i, module_path in enumerate(modules):
        if 'chden_plotting_needs_work_to_update.plotdensit' in module_path:
            if 'chden_plotting_needs_work_to_update.plotdensity' in module_path:
                print('plotdensity in modules')
            else:
                # Replace with plotdensity
                logging.warning(f'plotdensit in modules, replacing with plotdensity\n'
                                f'Original module path: {module_path}')
                module_path = module_path.replace('plotdensit', 'plotdensity')
                modules[i] = module_path

    return modules


def update_conf_py(documentation_dir, source_dir):
    conf_py_path = os.path.join(documentation_dir, 'conf.py')
    source_conf_py_path = os.path.join(documentation_dir, 'source', 'conf.py')

    with open(source_conf_py_path, 'r') as file:
        lines = file.readlines()

    with open(conf_py_path, 'w') as f:
        f.write("import os\nimport sys\n")

        # write the lines
        for line in lines:
            f.write(line)

        f.write(f"sys.path.insert(0, os.path.abspath('{source_dir}'))\n")
        f.write("\nhtml_theme = 'sphinx_rtd_theme'\n")
        f.write("\nextensions = [\n")
        f.write("    'sphinx.ext.autodoc',\n"
                "    'sphinx.ext.doctest',\n"
                "    'sphinx.ext.intersphinx',\n"
                "    'sphinx.ext.todo',\n"
                "    'sphinx.ext.coverage',\n"
                "    'sphinx.ext.mathjax',\n"
                "    'sphinx.ext.ifconfig',\n"
                "    'sphinx.ext.viewcode',\n"
                "    'sphinx.ext.githubpages',\n"
                "    'sphinx.ext.napoleon',\n"
                "    'sphinx.ext.autosummary',\n"
                "    'sphinx.ext.autosectionlabel',\n"
                "    'sphinx.ext.autodoc.typehints',\n"
                "    'sphinx.ext.inheritance_diagram',\n"
                "]\n")

    # Just in case copy the conf.py file to the source directory
    subprocess.run(['cp', conf_py_path, source_conf_py_path])


def create_module_rst_files(modules, rst_dir):
    for module in modules:
        with open(os.path.join(rst_dir, f'{module}.rst'), 'w') as f:
            f.write(
                f"{module}\n{'=' * len(module)}\n\n.. automodule:: {module}\n    :members:\n    :undoc-members:\n    :show-inheritance:\n")


def update_index_rst(documentation_dir, modules):
    index_rst_path = os.path.join(documentation_dir, 'index.rst')
    source_index_rst_path = os.path.join(documentation_dir, 'source/index.rst')
    with open(source_index_rst_path, 'r') as file:
        lines = file.readlines()

    with open(index_rst_path, 'w') as file:
        # Write descirption
        in_toctree = False
        for line in lines:
            if '.. toctree::' in line:
                in_toctree = True
                file.write(line)
                continue

            if in_toctree and line.strip() == ':caption: Contents:':
                file.write(line)
                file.write("\n")
                for module in sorted(set(modules)):  # Sort and remove duplicates
                    file.write(f"   {module}\n")
            elif in_toctree and not line.strip():
                # End of the toctree block
                in_toctree = False
                file.write(line)
            else:
                file.write(line)

    # Just in case copy the index.rst file to the source directory
    subprocess.run(['cp', index_rst_path, source_index_rst_path])


def build_sphinx_docs(documentation_dir):
    subprocess.run(['make', 'html'], cwd=documentation_dir)


def main(project_name, author_name, source_dir, documentation_dir, exclude_dirs=None, delete_old_files=False):
    if not exclude_dirs:
        exclude_dirs = []

    if delete_old_files and os.path.exists(documentation_dir):
        from densitypy.project_utils.file_directory_ops import delete_files_or_directories
        delete_files_or_directories(documentation_dir)

    if not os.path.exists(documentation_dir):
        os.makedirs(documentation_dir)

    initialize_sphinx(documentation_dir, project_name, author_name=author_name)
    modules = find_python_modules(source_dir, ignore_folders=exclude_dirs)
    update_conf_py(documentation_dir, source_dir)
    create_module_rst_files(modules, documentation_dir)
    update_index_rst(documentation_dir, modules)
    # Move the source files into the source directory
    execute_command(f"mv {documentation_dir}/*.rst {documentation_dir}/source/")
    build_sphinx_docs(documentation_dir)


if __name__ == '__main__':
    main(
        project_name='DensityPy',
        author_name='Ruben Fernandez Carbon',
        source_dir='/home/ruben/PycharmProjects/DensityPy/densitypy',
        documentation_dir='/home/ruben/PycharmProjects/DensityPy/Documentation/pythondocs',
        exclude_dirs=['venv', 'Studies', 'old_files', 'densityfort',
                      'densitypy/frontend_tbd','frontend_tbd',
                      'Documentation/pythondocs'],
        delete_old_files=True)
