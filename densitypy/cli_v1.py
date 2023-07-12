import click


@click.group()
@click.option("-ini",
              help="Specify configuration file, otherwise chargemigration.ini will be used as default or created",
              type=str, default=None)
@click.pass_context
def main(ctx, ini):
    click.echo("Executing main")

    ctx.obj = {}
    ctx.obj["json_config_path"] = ini

    print("ini", ini)
    print("ctx.obj", ctx.obj)

@main.group()
def molcas():
    pass


@molcas.command()
@click.option("-i", "--input", help="input file for OpenMolcas (pymolcas command). "
                                    "Used to indicate input file for Molcas. Name should differ from "
                                    "<ProjectName>.input to avoid overwriting.",
              type=str, default=None)
@click.option("-i?", "--input-help", help="Helps user by creating a starting input file for OpenMolcas (pymolcas command),"
                                           " which can be edited. To use open molcas with this file use -i "
                                           "followed by filename",
              is_flag=True)
@click.option("-lus", "--luscus", help="If true, Uses Luscus to select active space file "
                                        "($ProjectName.GvOrb)",
              is_flag=True)
@click.option("-h5", "--just-h5", help="This will only extract data found inside <Project>.rasscf.h5 "
                                       "from a previous run of Molcas.(Looks first inside of <Molcas_WORKDIR>)",
              is_flag=True)
@click.option("-getdensities", "--just-get-densities", help="This will only extract data found inside "
                                                             "<Project>.grid from a previous run of Molcas."
                                                             "(Looks first in the current directory)",
              is_flag=True)
@click.option("-getdipoles", "--just-get-dipoles", help="This will only extract the dipoles inside of "
                                                         "<Project>.log from a previous run of Molcas."
                                                         "(Looks first in the current directory)",
              is_flag=True)
@click.option("-gridflag", "--grid-flag", help="Required if grid is needed for post-molcas orbital density calculations",
              is_flag=True)
@click.option("-smartgrid", "--smart-grid", help="Automatically makes a grid for an arbitrary molecule using the size steps "
                                                 "and Boundary(distance from the nuclei to us as cutoff for grid) input. Can be used with optional flag \"-limitedgrid\"."
                                                 "Typically results in a cubic grid with unequal side lengths",
              is_flag=True)
@click.option("-limitedgrid", "--limited-grid",
              help="When used with flag \"-smartgrid\", new grid will only be composed of points \"Boundary\" distance away from each nuclei (reduces # points)."
                   "Typically results in oval grid",
              is_flag=True)
@click.option("-g", "--gridfile", help="input grid for pymolcas. (optional)",
              type=str, default=None)
@click.option("-r", "--just-run-input",
              help="Just run the pymolcas input file",
              is_flag=True)
@click.option("-old", "--old-main",
              help="use old main, debugging, rm",
              is_flag=True)
@click.pass_context
def run(ctx, input, input_help, luscus, just_h5, just_get_densities, just_get_dipoles, grid_flag, smart_grid, limited_grid, gridfile, just_run_input, old_main):
    click.echo("Executing molcas_input")

    ctx.obj["input"] = input
    ctx.obj["input_help"] = input_help
    ctx.obj["lustrue"] = luscus
    ctx.obj["justh5"] = just_h5
    ctx.obj["justgetdensity"] = just_get_densities
    ctx.obj["justgetdipoles"] = just_get_dipoles
    ctx.obj["gridflag"] = grid_flag
    ctx.obj["smartgrid"] = smart_grid
    ctx.obj["limitedgrid"] = limited_grid
    ctx.obj["gridfile"] = gridfile
    ctx.obj["just_run_input"] = just_run_input
    ctx.obj["old_main"] = old_main





@main.group()
def charge_migration():
    pass


@charge_migration.command()
@click.option("-cm", "--chargemigrationtrue", help="If True, calls Charge Migration "
                                                    "using the parameters in chargemigration.ini ",
              is_flag=True)
@click.option("-field", "--givenfieldfile", help="Specifies field file containing the sequences of pulses to be "
                                                 "used by ChargeMigration; otherwise the default is to look for \'Field_pulses\' "
                                                 "or make one using information from chargemigration.ini. Works for both "
                                                 "ChargeMigration or ChargeMigrationFT",
              type=str, default=None)
@click.option("-w", "--weights-file", help="Feed ChargeMigration code a premade File containing Becke's Weights for the Grid, "
                                           "saves time in subsequent calculations after already computed",
              type=str, default=None)
@click.option("-field?", "--fieldfilehelp", help="Creates example template file which can be used with \'-field\' to specify the "
                                                 "sequences of pulses to be used by ChargeMigration.",
              is_flag=True)
@click.option("-sden", "--writeCM", help="Save Charge Density to file when running '-cm'(time-consuming!)",
              nargs=2, type=float)
@click.option("-p", "--parallel", help="Run in parallel",
              is_flag=True)
@click.option("-debug", "--debug", help="Runs ChargeMigration or ChargeMigrationFT in debug mode",
              is_flag=True)
@click.pass_context
def charge_migration_input(ctx, chargemigrationtrue, givenfieldfile, weights_file, fieldfilehelp, writeCM, parallel, debug):
    ctx.obj["chargemigrationtrue"] = chargemigrationtrue
    ctx.obj["givenfieldfile"] = givenfieldfile
    ctx.obj["weights_file"] = weights_file
    ctx.obj["fieldfilehelp"] = fieldfilehelp
    ctx.obj["writeCM"] = writeCM
    ctx.obj["parallel"] = parallel
    ctx.obj["debug"] = debug
    print(ctx.obj)

@main.group()
def charge_migration_ft():
    pass


@charge_migration_ft.command()
@click.option("-ft", help="If True, calls Charge Migration FT"
                          "using the parameters in chargemigration.ini ",
              is_flag=True)
@click.option("-s", help="Calls Reconstruction of Spectrum code and scripts",
              is_flag=True)
@click.option("-save",
              help="Saves FT and difference in spectra files from previous runs and renames them with useful information if "
                   "they are found within the current {SIM} directory in use or the working director in the case of the difference",
              is_flag=True, default=True)
@click.pass_context
def charge_migration_ft_input(ctx, ft, s, save):
    ctx.obj["chargemigrationFTtrue"] = ft
    ctx.obj["SpectrumTrue"] = s
    ctx.obj["save_previous"] = save


@main.group()
def processing_program_x():
    pass


@processing_program_x.command()
@click.option("-ana", help="...", is_flag=True)
@click.option("-pl", help="...", is_flag=True)
@click.option("-ki", help="...", is_flag=True)
@click.option("-study", help="...", type=str, default=None)
@click.pass_context
def processing_program_x_input(ctx, ana, pl, ki, study):
    ctx.obj["run_Analysis"] = ana
    ctx.obj["run_Plot"] = pl
    ctx.obj["run_Kinetics"] = ki
    ctx.obj["study_directory"] = study


if __name__ == "__main__":
    main(obj={})
