import click


@click.command()
@click.argument("json_config_path")
@click.argument("study_directory")
@click.option("--molcas_input", help="Path to the Molcas input file, this enables molcas execution")
@click.option("--run_charge_migration", is_flag=True, help="Enable charge migration")
@click.option("--run_charge_migration_ft", is_flag=True, help="Enable charge migration FT")
@click.option("--run_spectrum_reconstruction", is_flag=True, help="Enable spectrum reconstruction")
@click.option("--field_file_help", is_flag=True, help="Enable field file help")
@click.option("--molcas_input_help", is_flag=True, help="Enable Molcas input help")
@click.option("--lus", is_flag=True, help="Enable LUS")
@click.option("--gridflag", is_flag=True, help="Enable grid flag")
@click.option("--write_charge_migration", help="Write charge migration")
@click.option("--debug_mode", is_flag=True, help="Enable debug mode")
@click.option("--justh5", is_flag=True, help="Enable just h5")
@click.option("--justgetdipoles", is_flag=True, help="Enable just get dipoles")
@click.option("--justgetdensity", is_flag=True, help="Enable just get density")
@click.option("--weights_file", help="Path to the weights file")
@click.option("--givenfieldfile", help="Path to the given field file")
@click.option("--old_main", is_flag=True, help="Enable old main")
@click.option("--parallel", is_flag=True, help="Enable parallel execution")
@click.option("--save_previous", is_flag=True, help="Enable saving previous results")
@click.option("--make", is_flag=True, help="Enable Fortran compilation", default=False)
@click.option("--make_directory", help="Path to the Fortran compilation directory", required=False, default="")
@click.option("--make_flags", help="Path to the Fortran compilation make flags", required=False, default="")
@click.option("--plot", is_flag=True, help="Enable plotting", default=False)
def cli_run(json_config_path, study_directory, molcas_input, run_charge_migration,
            run_charge_migration_ft, run_spectrum_reconstruction, field_file_help, molcas_input_help, lus, gridflag,
            write_charge_migration, debug_mode, justh5, justgetdipoles, justgetdensity, weights_file, givenfieldfile,
            old_main, parallel, save_previous,
            make, make_directory, make_flags, plot):
    # Modify the run function to include your desired functionality
    from densitypy.main import run_densitypy
    run_densitypy(json_config_path, study_directory, molcas_input, run_charge_migration,
                  run_charge_migration_ft, run_spectrum_reconstruction, field_file_help, molcas_input_help, lus,
                  gridflag, write_charge_migration, debug_mode, justh5, justgetdipoles, justgetdensity, weights_file,
                  givenfieldfile, old_main, parallel, save_previous,
                  make_fortran=make,
                  make_fortran_config=dict(directory=make_directory, make_flags=make_flags),
                  plot=plot

                  )


if __name__ == "__main__":
    cli_run()
