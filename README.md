# DensityPy


Click CLI found in [`cli_run.py`](cli_run.py) and function [`run_densitypy`](densitypy%2Fmain.py) are the entry Way to 
DensityPy. This cli/function facilitates computational chemistry simulations with OpenMolcas and the 
ASTRA-ChargeMigration Fortran code; simulating the charge migration of a molecule from MCSCF calculations and generating plots from the observables of pump-probe experiments. Integration with other frameworks has been experimented with, but is not currently supported out of the box.

<p align="center">
  <img src="https://github.com/Bucanero06/Bucanero06/assets/60953006/e9e8a290-9e74-4d45-96ae-4114e423f637" height="248" />
  <img src="https://github.com/Bucanero06/DensityPy/assets/60953006/46fdf9fc-2a13-4836-98dc-3c81ffd279e5" height="248" />
 <img src="https://github.com/Bucanero06/DensityPy/assets/60953006/190c733a-d935-44e8-9e9e-280d8010714b" height="248" />
</p>

Functionality and Scope:
* [Molcas Integration](densitypy%2Fmolcas%2Fmolcasscripts.py): Incorporates Molcas for quantum chemical calculations, essential for accurate modeling of electronic structures.
* [Charge Migration Analysis](densitypy%2Fcharge_migration%2Fchargemigratonscripts.py): Simulates the charge migration of a system with time under the influence of an external electric field and its environment.
* [Data Visualization and Analysis](densitypy%2Fpost_processing%2Fplotting_module.py): Provides tools for generating plots and graphs, crucial for interpreting complex simulation data.
* [Configurable Workflow](densitypy%2Fproject_utils%2Fconfiguration_parser.py): Leverages a [JSON-based configuration](densitypy%2FDefault_Settings%2Fdefault_config.py) system, allowing for flexible and detailed setup of simulation parameters.
* [Environment Management](densitypy%2Fproject_utils): Implements efficient directory and environment management, crucial for handling extensive simulations.
* [Logging](densitypy%2Fproject_utils%2Flogger.py) and [Debugging](densitypy%2Fproject_utils%2Ffortran_compilation_handler.py): Includes comprehensive logging for monitoring simulation progress and facilitating debugging.
* [Fortran Compatibility](densitypy%2Fproject_utils%2Ffortran_compilation_handler.py): Uses Fortran for [performance-critical components](densityfort) of the simulations, often parallelized for optimal performance.
* [AutoDocumentation](Documentation): Using both [Sphinx](https://www.sphinx-doc.org/en/master/) and [FORD](forddocs.readthedocs.io/en/latest/), allowing for easy navigation and understanding of the codebase.

*all the processes above have plenty of room for improvement*

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
              tidy_up_experiment_dir=False, compress_study_directory=False
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

## Quick Start
You are welcome to deviate from this quick start guide, but this may be a great place to start prior to diving into the
codebase. Note that DensityPy leverages OpenMolcas (and soon ASTRA) for quantum chemical calculations, essential for 
accurate modeling of electronic structures, 
please [install](https://molcas.gitlab.io/OpenMolcas/sphinx/installation.guide/ig.html) and make sure h5 and 
GridIt(*deprecated*) configurations are enabled.

### The Command Line Options are:
```
Options:
  --molcas_input TEXT            Path to the Molcas input file, this enables molcas execution
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
  --tidy_up_experiment_dir       Removes Pulses dir and eletes  PP/FTPP files in AtomicCharge & Dipole
  --compress_study_directory     Compress study directory `tar -czvf {study_directory}.tar.gz {study_directory}`
  --help                         Show this message and exit.
```

### Configuring the Study Directory
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

We already could run the code, everything has been prefilled for you. However, 
if we know nothing about a system's orbitals we can begin with some chemical intuitions and a SCF calculation and 
generate any of the following files in order to visualize the orbitals and select an appropriate active space:
* HDF5 files, as generated by some (Open)Molcas modules like SCF, RASSCF or RASSI, if compiled with HDF5 support.
* Molden files.
* Luscus files, generated by the GRID_IT module.
* Grid files (ASCII), generated by the GRID_IT module.
* Cube files (formatted).

To quickly do this, we can run the `scforbs` command, or for a more automated approach use [`autocas`]()(expensive):
```bash
python cli_run.py Studies/ExampleStudy/ nma_configuration.json --scforbs 
  # --autocas --> automates active-orbital-space selection step in multi-configurational calculations 
```
This will, among others, generate a h5 containing information about each orbital which will then be loaded and 
visualized. 

[Screencast from 12-17-2023 11:47:17 PM.webm](https://github.com/Bucanero06/DensityPy/assets/60953006/fdc60b98-4231-49ee-9d94-5747cb7f711b)

![Screenshot from 2023-12-17 23-41-13](https://github.com/Bucanero06/DensityPy/assets/60953006/0e2a9550-5d0b-44ed-8d5b-ed29bb68deed)

We can then select the orbitals we want to include in our initial active space. For this 
example, to compute the ground and excited states of molecules such as formamide and NMA, let's perform a 
`RAS-SCF` calculation using a `6-31G$^{**}` basis. This basis choice is effective for organic compounds as it 
includes polarization functions on all atoms and adds diffuse functions on heavy atoms. In our calculations, 
we'll focus on critical orbitals, including the highest occupied orbitals, which are typically lone pairs on 
oxygen, and the highest π orbitals. Our active space will also encompass several important unoccupied Rydberg 
and valence orbitals. Although this selection can lead to large and computationally challenging active spaces, 
a more focused approach considering only valence orbitals can simplify the problem. We'll include key unoccupied 
orbitals, like the σ* and π* orbitals, in our active space. This strategy not only makes the calculations more 
manageable but also allows for the inclusion of additional relevant occupied orbitals, such as another lone pair 
on oxygen, thus providing a comprehensive yet efficient and targeted computational run.

```bash
python cli_run.py Studies/ExampleStudy/ nma_configuration.json 
    --molcas_input Studies/ExampleStudy/molcas_input_help.input --gridit 
```

The above command adds the gridcoordinates to the provided file, runs OpenMolcas, and finally extracts the nessesary 
content from the logs and h5 files, preparing them for the charge migration modules.

Until this point, we have been relying on the default configuration file, which is a JSON file that contains the 
parameters for controlling simulations, including molcas. You could have ran the entire pipeline as followed:

```bash
python cli_run.py Studies/ExampleStudy/ nma_configuration.json 
    --molcas_input Studies/ExampleStudy/molcas_input_help.input --gridit   
    --run_charge_migration --run_charge_migration_ft --run_spectrum_reconstruction --plot
```

But this would be inappropriate for any other case, as we would want to change the parameters based on the characteristics
of the states and even conduct multiple geometry optimizations. For reference, this is what the configuration file we 
generated looks like:

```json
{
    "Project settings": {
        "project_name": "NMA",
        "XYZ Molecule Geometry": "NMA.xyz",
        "Number of States": 4,
        "List_of_Active_Orbitals": [
            18,
            19,
            20,
            21
        ],
        "Molcas Output Directory": "NMA_output",
        "Experiment Directory": "example_sim"
    },
    "GRID settings": {
        "Number of Points X axis": 100,
        "Number of Points Y axis": 100,
        "Number of Points Z axis": 100,
        "X MIN": -6,
        "X MAX": 7,
        "Y MIN": -5,
        "Y MAX": 6,
        "Z MIN": -5,
        "Z MAX": 5,
        "Step Size": 0.008,
        "Boundary": 6
    },
    "Charge Migration settings": {
        "Field File": "field_file",
        "Number of Times": 9001,
        "Min Time": -5000,
        "Max Time": 13000,
        "Bath Temperature": 3000,
        "Dephasing Factor": 0.0012,
        "Relaxation Factor": 0.0048
    },
    "Pump Pulses settings": {
        "Type of Pulse": "G",
        "start_time": 0,
        "Pump Central Frequency": 0.41,
        "Pump Periods": 2,
        "Pump Phase": 0,
        "Pump Intensity": 0.12,
        "Pump Polarization": [
            90,
            90
        ]
    },
    "Probe Pulses settings": {
        "Type of Pulse": "G",
        "Time Delay Start": -200,
        "Time Delay Stop": 800,
        "Number Of PP": 120,
        "Time Delay Weight Factor": 0.5,
        "Probe Central Frequency": 0.073, 
        "Probe Periods": 2,
        "Probe Phase": 0,
        "Probe Intensity": 0.12,
        "Probe Polarization": [
            90,
            90
        ]
    },
    "Charge Migration FT settings": {
        "Number of Omegas": 120,
        "Min Omega": 0.18,
        "Max Omega": 0.68,
        "Number of TauOmegas": 120,
        "Min TauOmega": -0.25,
        "Max TauOmega": 0.25,
        "FT TimeStep": 5000,
        "FT WidthStep": 350
    }
}
```
We started with very conservative parameters to quickly iterate over settings and get a feel for the system. For example, 
we can change the number of points and size of the grid, the polarization
of the pump and probe pulses, the number of states, the number of times, omegas, ... as well as other Fourier Transform parameters; 
all which should be considered in order to minimize artifacts and maximize the accuracy of the results. More importantly, 
dephasing & relaxation factors can be changed to account for the environment and the system's dynamics. #TODO search here

```commandline
 Loading Energies
2023-12-20T03:54:06.993796Z [INFO] Charge_Migration: State Energies
2023-12-20T03:54:06.993840Z [INFO] Charge_Migration: 1  -246.525708599498
2023-12-20T03:54:06.993878Z [INFO] Charge_Migration: 2  -246.344501818203
2023-12-20T03:54:06.993913Z [INFO] Charge_Migration: 3  -246.189490935279
2023-12-20T03:54:06.993947Z [INFO] Charge_Migration: 4  -246.047815060530
2023-12-20T03:54:06.993979Z [INFO] Charge_Migration: 
2023-12-20T03:54:06.994015Z [INFO] Charge_Migration: Loading Transition Density Matrix
2023-12-20T03:54:06.994100Z [INFO] Charge_Migration: Number of Active Orbitals: 4
2023-12-20T03:54:06.994459Z [INFO] Charge_Migration: Loading Dipole Matrices for XYZ components
2023-12-20T03:54:06.994520Z [INFO] Charge_Migration: 
2023-12-20T03:54:06.994554Z [INFO] Charge_Migration: NMA_output/X_DIPOLE.csv
2023-12-20T03:54:06.994596Z [INFO] Charge_Migration: 0.347947E+00   0.456999E-06   0.992338E+00   0.274939E+00
2023-12-20T03:54:06.994632Z [INFO] Charge_Migration: 0.456999E-06   0.215518E+00   0.113281E-05  -0.165158E-05
2023-12-20T03:54:06.994667Z [INFO] Charge_Migration: 0.992338E+00   0.113281E-05   0.172030E+01   0.273282E+00
2023-12-20T03:54:06.994703Z [INFO] Charge_Migration: 0.274939E+00  -0.165158E-05   0.273282E+00   0.809354E+00
2023-12-20T03:54:06.994757Z [INFO] Charge_Migration: 
2023-12-20T03:54:06.994790Z [INFO] Charge_Migration: NMA_output/Y_DIPOLE.csv
2023-12-20T03:54:06.994830Z [INFO] Charge_Migration: -0.111155E+01  -0.668804E-06  -0.921473E+00   0.112078E+01
2023-12-20T03:54:06.994866Z [INFO] Charge_Migration: -0.668804E-06  -0.560969E+00  -0.143667E-05   0.850321E-06
2023-12-20T03:54:06.994899Z [INFO] Charge_Migration: -0.921473E+00  -0.143667E-05  -0.205932E+01   0.188877E+00
2023-12-20T03:54:06.994932Z [INFO] Charge_Migration: 0.112078E+01   0.850321E-06   0.188877E+00  -0.203919E+01
2023-12-20T03:54:06.994981Z [INFO] Charge_Migration: 
2023-12-20T03:54:06.995014Z [INFO] Charge_Migration: NMA_output/Z_DIPOLE.csv
2023-12-20T03:54:06.995053Z [INFO] Charge_Migration: 0.585080E-05   0.683533E-01   0.103946E-05   0.630551E-06
2023-12-20T03:54:06.995087Z [INFO] Charge_Migration: 0.683533E-01   0.501887E-05   0.268546E-01  -0.723624E-01
2023-12-20T03:54:06.995121Z [INFO] Charge_Migration: 0.103946E-05   0.268546E-01   0.677975E-05   0.545977E-06
2023-12-20T03:54:06.995156Z [INFO] Charge_Migration: 0.630551E-06  -0.723624E-01   0.545977E-06   0.613159E-05
2023-12-20T03:54:06.995200Z [INFO] Charge_Migration: Loading Geometry
2023-12-20T03:54:06.995274Z [INFO] Charge_Migration: nAtoms          12
2023-12-20T03:54:06.995307Z [INFO] Charge_Migration: Loading Grid
2023-12-20T03:54:34.251353Z [INFO] Charge_Migration: Loading Orbitals
2023-12-20T03:54:46.280446Z [INFO] Charge_Migration: Computing Weights
2023-12-20T03:55:26.524639Z [INFO] Charge_Migration: Writing Weights to File
2023-12-20T03:57:24.469582Z [INFO] Charge_Migration: Computing Barycenters of the Atomic Charges
2023-12-20T03:57:40.215674Z [INFO] Charge_Migration: Computing Becke's Matrix
...
...
```

#TODO EXTEND DOCUMENTATION FOR THESE STEPS
#TODO PARAMSEARCH goes here too
After finishing, the pipeline can be visualized (although changing) through interactive and png plots like: 
#TODO ADD PLOTS --placeholder->

[Screencast from 12-20-2023 05:09:56 PM.webm](https://github.com/Bucanero06/DensityPy/assets/60953006/a1f04802-067d-4d8b-a596-5b06ed591aee)

![Screenshot from 2023-12-20 17-16-45](https://github.com/Bucanero06/DensityPy/assets/60953006/69fee775-9d5f-4ead-a0fc-25bd6dcc4b16)

*caption: Dipole and reconstructed Dipole from Atomic Charges, left and right figures respectively.*

<p align="center">
  <img src="https://github.com/Bucanero06/DensityPy/assets/60953006/9e07f8bf-f922-41a2-bc1b-0d2018bdf595" height="400" />
 <img src="https://github.com/Bucanero06/DensityPy/assets/60953006/e5144cf3-660f-4c27-85b1-bd6c6ec1006f" height="400" />
</p>



`I am currently having difficulties with my knight's email`, please reach out to me here, `Github`, otherwise either 
`ruben.fernandez.carbon@gmail.com` or `ruben@carbonyl.org` to DM me.

This is a part of a larger codebase, currently still undergoing a code review and refactoring along the manuscript and 
thesis writing process. Treat this codebase as a proof of concept, and not as final; there are still bugs and discrepancies.

With modifications, the outputs from OpenMolcas, PySCF, GAMES, and NWChem have been tested and are 
compatible with the Fortran codebase, documentation as it stands is available 
[here](https://bucanero06.github.io/DensityPy) and will be regularly updated as the codebase is refactored.
Right now the python documentation is not showing up as intended on github pages, but it is available locally.


## Documentation
The documentation is currently being updated, but you can find the latest version 
[here.](https://bucanero06.github.io/DensityPy)
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
        --graphs True --remove_old_files
        --exclude_dirs "densityfort/ChargeMigration/src/future_reconstruction_using_ci_work,densityfort/ChargeMigrationFT/src/future_reconstruction_using_ci_work" 
```

![Screenshot from 2023-12-19 22-33-54.png](docs%2Fimages%2FScreenshot%20from%202023-12-19%2022-33-54.png)
![Screenshot from 2023-12-19 22-32-06.png](docs%2Fimages%2FScreenshot%20from%202023-12-19%2022-32-06.png)
The documentation can be generated from docstrings for both modules and functions, handling markdown as well as latex and the generation 
of inheritance diagrams.

![Screenshot from 2023-12-20 13-15-46.png](docs%2Fimages%2FScreenshot%20from%202023-12-20%2013-15-46.png)
![Screenshot from 2023-12-20 13-14-44.png](..%2F..%2FPictures%2FScreenshots%2FScreenshot%20from%202023-12-20%2013-14-44.png)

Notes:
The most convenient way to install a patched PySCF together with Block is to use the build script [PySCF_Block.sh](PySCF_Block.sh)

