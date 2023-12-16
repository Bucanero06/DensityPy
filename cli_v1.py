import click


@click.command()
@click.argument("json_config_path")
@click.argument("study_directory")
@click.option("--molcas_input", help="Path to the Molcas input file, this enables molcas execution")
@click.option("--run_charge_migration", is_flag=True, help="Enable charge migration", default=False)
@click.option("--run_charge_migration_ft", is_flag=True, help="Enable charge migration FT", default=False)
@click.option("--run_spectrum_reconstruction", is_flag=True, help="Enable spectrum reconstruction", default=False)
@click.option("--field_file_help", is_flag=True, help="Enable field file help", default=False)
@click.option("--molcas_input_help", is_flag=True, help="Enable Molcas input help", default=False)
@click.option("--lus", is_flag=True, help="Enable LUS", default=False)
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
            run_charge_migration_ft, run_spectrum_reconstruction, field_file_help, molcas_input_help, lus, gridit,
            write_charge_migration, debug_mode, justh5, justgetdipoles, justgetdensity, weights_file, givenfieldfile,
            make, make_directory, make_flags, plot):
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
        lus=lus, gridit=gridit, write_charge_migration=write_charge_migration, debug_mode=debug_mode,
        justh5=justh5, justgetdipoles=justgetdipoles, justgetdensity=justgetdensity, weights_file=weights_file,
        givenfieldfile=givenfieldfile,
        make_fortran=make,
        make_fortran_config=dict(directory=make_directory, make_flags=make_flags),

    )


if __name__ == "__main__":
    cli_run()
