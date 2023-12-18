import click


@click.command()
@click.argument("study_directory")
@click.argument("json_config_path")
@click.option("--molcas_input", help="Path to the Molcas input file, this enables molcas execution")
@click.option("--run_charge_migration", is_flag=True, help="Enable charge migration", default=False)
@click.option("--run_charge_migration_ft", is_flag=True, help="Enable charge migration FT", default=False)
@click.option("--run_spectrum_reconstruction", is_flag=True, help="Enable spectrum reconstruction", default=False)
@click.option("--field_file_help", is_flag=True, help="Enable field file help", default=False)
@click.option("--molcas_input_help", is_flag=True, help="Enable Molcas input help", default=False)
@click.option("--scforbs", is_flag=True, help="Run SCF Orbital Calculations and plot", default=False)
@click.option("--gridit", is_flag=True, help="Enable grid flag", default=False)
@click.option("--write_charge_migration", is_flag=True, help="Write charge migration", default=False)
@click.option("--debug_mode", is_flag=True, help="Enable debug mode", default=False)
@click.option("--justh5", is_flag=True, help="Enable just h5", default=False)
@click.option("--justgetdipoles", is_flag=True, help="Enable just get dipoles", default=False)
@click.option("--justgetdensity", is_flag=True, help="Enable just get density", default=False)
@click.option("--weights_file", help="Path to the weights file")
@click.option("--givenfieldfile", help="Path to the given field file")
@click.option("--make", is_flag=True, help="Enable Fortran compilation", default=False)
@click.option("--make_directory", help="Path to the Fortran compilation directory", required=False, default="")
@click.option("--make_flags", help="Path to the Fortran compilation make flags", required=False,
              default="all DEB_FLAG=d")
@click.option("--plot", is_flag=True, help="Enable plotting", default=False)
def cli_run(json_config_path, study_directory, molcas_input, run_charge_migration,
            run_charge_migration_ft, run_spectrum_reconstruction, field_file_help, molcas_input_help, scforbs, gridit,
            write_charge_migration, debug_mode, justh5, justgetdipoles, justgetdensity, weights_file, givenfieldfile,
            make, make_directory, make_flags, plot):
    """
     Entry Way to DensityPy
    This CLI is engineered to facilitate computational chemistry simulations with OpenMolcas and the ASTRA-ChargeMigration Fortran code.

    Functionality and Scope:
        Molcas Integration: Incorporates Molcas for quantum chemical calculations, essential for accurate modeling of electronic structures.
        Charge Migration Analysis: Executes both standard and Fourier-transformed charge migration simulations, offering deep insights into molecular dynamics.
        Data Visualization and Analysis: Provides tools for generating plots and graphs, crucial for interpreting complex simulation data.
        Configurable Workflow: Leverages a JSON-based configuration system, allowing for flexible and detailed setup of simulation parameters.

    Parameters
        json_config_path: String, path to the JSON configuration file.
        study_directory: String, directory path for output and intermediate files.
        molcas_input: String, optional path to the Molcas input file.
        Additional flags and options to control specific features and modes of operation.

    Usage:
        python cli_v1.py configuration_help.json /home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/ --make --make_directory=densityfort --make_flags='all' --run_charge_migration --run_charge_migration_ft --run_spectrum_reconstruction --plot

    """

    from densitypy.main import run_densitypy
    run_densitypy(
        json_config_path=json_config_path,
        study_directory=study_directory,
        molcas_input=molcas_input,  # False just means not running molcas
        run_charge_migration=run_charge_migration,
        run_charge_migration_ft=run_charge_migration_ft,
        run_spectrum_reconstruction=run_spectrum_reconstruction,
        plot=plot,
        #
        field_file_help=field_file_help, molcas_input_help=molcas_input_help,
        scforbs=scforbs, gridit=gridit, write_charge_migration=write_charge_migration, debug_mode=debug_mode,
        justh5=justh5, justgetdipoles=justgetdipoles, justgetdensity=justgetdensity, weights_file=weights_file,
        givenfieldfile=givenfieldfile,
        make_fortran=make,
        make_fortran_config=dict(directory=make_directory, make_flags=make_flags),

    )


if __name__ == "__main__":
    cli_run()
