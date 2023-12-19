python -m densitypy.autodocumentation_python -p DensityPy -a "Ruben Fernandez Carbon" -s /home/ruben/PycharmProjects/DensityPy -d /home/ruben/PycharmProjects/DensityPy/docs/pythondocs -e venv,Studies,old_files,densityfort,densitypy/frontend_tbd,frontend_tbd,Documentation/pythondocs -r
python -m densitypy.autodocumentation_fortran --project_name "Charge Migration" --project_description "Fortran code for the calculation of the charge migration in molecules. Interfaced with Python's DensityPy library." --source_dirs "densityfort/ChargeMigration,densityfort/ChargeMigrationFT,densityfort/SpectrumReconstruction" --documentation_dir "docs/fortrandocs" --exclude_dirs "densityfort/ChargeMigration/src/future_reconstruction_using_ci_work,densityfort/ChargeMigrationFT/src/future_reconstruction_using_ci_work" --graphs --remove_old_files