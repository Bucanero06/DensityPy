# DensityPy


Click CLI found in [`cli_run.py`](cli_run.py) and function [`run_densitypy`](densitypy%2Fmain.py) are the entry Way to 
DensityPy. This cli/function is engineered to facilitate computational chemistry simulations with OpenMolcas and the 
ASTRA-ChargeMigration Fortran code.
<p align="center">
  <img src="https://github.com/Bucanero06/Bucanero06/assets/60953006/e9e8a290-9e74-4d45-96ae-4114e423f637" width="400" />
 <img src="https://github.com/Bucanero06/DensityPy/assets/60953006/07edfb34-672b-4151-bb4b-14560f44fac7" height="248" />

</p>

Functionality and Scope:
* [Molcas Integration](densitypy%2Fmolcas%2Fmolcasscripts.py): Incorporates Molcas for quantum chemical calculations, essential for accurate modeling of electronic structures.
* [Charge Migration Analysis](densitypy%2Fcharge_migration%2Fchargemigratonscripts.py): Simulates the charge migration of a system with time under the influence of an external electric field and its environment.
* [Data Visualization and Analysis](densitypy%2Fpost_processing%2Fplotting_module.py): Provides tools for generating plots and graphs, crucial for interpreting complex simulation data.
* [Configurable Workflow](densitypy%2Fproject_utils%2Fconfiguration_parser.py): Leverages a [JSON-based configuration](densitypy%2FDefault_Settings%2Fdefault_config.py) system, allowing for flexible and detailed setup of simulation parameters.
* [Environment Management](densitypy%2Fproject_utils): Implements efficient directory and environment management, crucial for handling extensive simulations.
* [Logging](densitypy%2Fproject_utils%2Flogger.py) and [Debugging](densitypy%2Fproject_utils%2Ffortran_compilation_handler.py): Includes comprehensive logging for monitoring simulation progress and facilitating debugging.
* [Fortran Compatibility](densitypy%2Fproject_utils%2Ffortran_compilation_handler.py): Features integration with Fortran for [performance-critical components](densityfort) of the simulations, often parallelized for optimal performance.
* [AutoDocumentation](Documentation): Using both [Sphinx](https://www.sphinx-doc.org/en/master/) and [FORD](forddocs.readthedocs.io/en/latest/), allowing for easy navigation and understanding of the codebase.


Usage:
```python
run_densitypy(json_config_path="xy_polarized_config.json",
              study_directory="Studies/ExampleStudy",
              molcas_input=False,  # 'molcas_input_help.input', # False just means not running molcas
              run_charge_migration=True,
              run_charge_migration_ft=True,
              run_spectrum_reconstruction=True,
              plot=True,
              #
              field_file_help=False, molcas_input_help=False,
              lus=False, gridit=True, write_charge_migration=None, debug_mode=False,
              justh5=False, justgetdipoles=False, justgetdensity=False, weights_file=None, givenfieldfile=None,
              make_fortran=False, make_fortran_config={'directory': '~/DensityPy/densityfort',
                                                       'make_flags': 'all DEB_FLAG=d'} 
              )
```
or
```bash
python cli_v1.py ~/DensityPy/Studies/ExampleStudy/ configuration_help.json \
  --make --make_directory=densityfort --make_flags='all' \
  --run_charge_migration --run_charge_migration_ft --run_spectrum_reconstruction --plot
```

I am currently having difficulties with my knight's email, please reach out to me here, Github, otherwise either ruben.fernandez.carbon@gmail.org or ruben@carbonyl.org to DM me.

This is a part of a larger codebase, currently still undergoing a code review and refactoring along the manuscript and 
thesis writing process. Threat this codebase as a proof of concept, and not as final; there are still bugs and discrepancies.

## Quick Start
You are welcome to deviate from this quick start guide, but this may be a great place to start prior to diving into the
codebase. Note that DensityPy leverages OpenMolcas (and soon ASTRA) for quantum chemical calculations, essential for 
accurate modeling of electronic structures, 
please [install](https://molcas.gitlab.io/OpenMolcas/sphinx/installation.guide/ig.html) and make sure h5 and 
GridIt(*on its way out*) configurations are enabled.

![Screenshot from 2023-12-17 20-41-16](https://github.com/Bucanero06/DensityPy/assets/60953006/9e90fce9-b7d0-40f5-b09c-63f088a1cce0)
We can begin from a single `.xyz` file, which contains the molecular geometry of the system of
interest. Let our model system be n-methylacetamide (NMA), which has the following geometry:

```
12
NMA
O          0.0000000000          1.2300000000          0.0000000000
N          1.1304300000         -0.7915360000          0.0000000000
C          0.0000000000          0.0000000000          0.0000000000
C          2.4748470000         -0.2483560000          0.0000000000
C         -1.2771790000         -0.7674070000          0.0000000000
H          1.0948330000         -1.8109140000          0.0000000000
H          3.1899510000         -1.0709890000          0.0000000000
H         -2.1320400000         -0.0751540000          0.0000000000
H         -1.2748780000         -1.4665740000         -0.8492120000
H         -1.2748810000         -1.4665680000          0.8492160000
H          2.5809700000          0.4362910000         -0.8414850000
H          2.5809800000          0.4362670000          0.8415030000
```

From here, either bring your own or create a configuration file with the following command:
```bash
python cli_run.py \
# Path to StudyDir - Relative/Absolute Path to Study Directory .
# Configuration Name - We can name the configuration file anything we want, 
#     either it exists or an example will be created by that name.
# We can also include other helpful flags (optional) when starting.
python cli_run.py \
    Studies/ExampleStudy/ nma_configuration.json \
    --molcas_input_help --field_file_help 
``` 
![Screenshot from 2023-12-17 23-29-44](https://github.com/Bucanero06/DensityPy/assets/60953006/639a99f2-e225-4df8-8a1f-56e34044045c)

We already could run the code along the lines of the usage shown above, everything has been prefilled for you. However, 
if we know nothing about a system's orbitals we can begin with some chemical intuitions and a SCF calculation and 
generate any of the following files in order to visualize the orbitals and select an appropriate active space:
* HDF5 files, as generated by some (Open)Molcas modules like SCF, RASSCF or RASSI, if compiled with HDF5 support.
* Molden files.
* Luscus files, generated by the GRID_IT module.
* Grid files (ASCII), generated by the GRID_IT module.
* Cube files (formatted).

To do this we can run the `scforbs` command:
```bash
python cli_run.py Studies/ExampleStudy/ nma_configuration.json --scforbs
```
This will among others generate a quick h5 containing information about each orbital which will then be loaded and 
visualized. 

[Screencast from 12-17-2023 11:47:17 PM.webm](https://github.com/Bucanero06/DensityPy/assets/60953006/fdc60b98-4231-49ee-9d94-5747cb7f711b)

![Screenshot from 2023-12-17 23-41-13](https://github.com/Bucanero06/DensityPy/assets/60953006/0e2a9550-5d0b-44ed-8d5b-ed29bb68deed)

We can then select the orbitals we want to include in our active space, and then run the code again 

The Options are:
```
Options:
  --molcas_input TEXT            Path to the Molcas input file, this enables
                                 molcas execution
  --run_charge_migration         Enable charge migration
  --run_charge_migration_ft      Enable charge migration FT
  --run_spectrum_reconstruction  Enable spectrum reconstruction
  --field_file_help              Enable field file help
  --molcas_input_help            Enable Molcas input help
  --scforbs                      Run SCF Orbital Calculations and plot
  --gridit                       Enable grid flag
  --write_charge_migration       Write charge migration
  --debug_mode                   Enable debug mode
  --justh5                       Enable just h5
  --justgetdipoles               Enable just get dipoles
  --justgetdensity               Enable just get density
  --weights_file TEXT            Path to the weights file
  --givenfieldfile TEXT          Path to the given field file
  --make                         Enable Fortran compilation
  --make_directory TEXT          Path to the Fortran compilation directory
  --make_flags TEXT              Path to the Fortran compilation make flags
  --plot                         Enable plotting
  --help                         Show this message and exit.
```


## Documentation
The documentation is currently being updated, but you can find the latest version [here](https://bucanero06.github.io/DensityPy/).
Running the following command will generate the documentation locally using Sphinx and FORD for 
[python](densitypy%2Fautodocumentation_python.py) and [fortran](densitypy%2Fautodocumentation_fortran.py) respectively:
```bash
bash update_documentation.sh
```
or if needed
```
python -m densitypy.autodocumentation_python 
        -p DensityPy -a {Author} -s {DENSITYPYROOTDIR} -d {OUTPUTDOCUMENTATIONDIR} -e [EXCLUDEDIRS] -r
        
python -m densitypy.autodocumentation_fortran 
        --project_name "Charge Migration" --project_description {...} 
        --source_dirs "densityfort/ChargeMigration,densityfort/ChargeMigrationFT,densityfort/SpectrumReconstruction" 
        --documentation_dir "docs/fortrandocs" 
        --exclude_dirs "densityfort/ChargeMigration/src/future_reconstruction_using_ci_work,densityfort/ChargeMigrationFT/src/future_reconstruction_using_ci_work" 
        --graphs True --remove_old_files
```
The documentation can handle

autodocumentation of docstrings for both modules and functions, handling markdown as well as latex and the generation 
of inheritance diagrams.

