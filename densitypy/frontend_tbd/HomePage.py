from h2o_wave import Q, ui, on, data, handle_on, AsyncSite  # noqa F401

from densitypy.charge_migration.chargemigratonscripts import generate_time_delays
from densitypy.project_utils.file_directory_ops import get_folders_in_directory, get_all_files_in_directory
from .util import clear_cards, add_card, load_page_recipe_with_update_models

app = AsyncSite()


async def homepage(q: Q):
    # Clear all cards except the ones needed for this page
    await q.run(clear_cards, q, ignore=['Application_Sidebar'])

    '''Context - Data from Firestore'''

    '''Static Cards'''
    # Add header
    add_card(q, 'ALM_Header', ui.header_card(box='header', title='Dashboard', subtitle='Ruben',
                                             # Color
                                             color='transparent',
                                             icon='AnalyticsView',
                                             icon_color=None,
                                             ))

    # Search the Studies directory for folders (studies)
    studies_names = get_folders_in_directory('/home/ruben/PycharmProjects/DensityPy/Studies/')
    q.user.recent_opened_study_directory = q.user.recent_opened_study_directory or studies_names[0]
    json_files_found_in_study = [json_file for json_file in get_all_files_in_directory(
        f'/home/ruben/PycharmProjects/DensityPy/Studies/{q.user.recent_opened_study_directory}')
                                 if json_file.endswith('.json')]
    q.user.recent_opened_configuration_file = q.user.recent_opened_configuration_file or json_files_found_in_study[0]

    add_card(q, 'ALM_Studies_Directory', ui.form_card(
        box='first_context_1',
        items=[
            ui.dropdown(name='ALM_studies_directory', label='Study Folder', required=True,
                        trigger=True,
                        choices=[
                            ui.choice(name=folder, label=folder) for folder in studies_names
                        ],
                        value=q.user.recent_opened_study_directory
                        ),
            # Configuration File in the study directory selected above,looking for json files

            ui.dropdown(name='ALM_configuration_file', label='Configuration File', required=True,
                        trigger=True,
                        choices=[
                            ui.choice(name=configuration_file, label=configuration_file) for configuration_file in
                            json_files_found_in_study
                        ],
                        value=q.user.recent_opened_configuration_file
                        ),
        ]
    ))

    # If Configuration File is selected, show the configuration file selected
    if q.user.configuration_file:
        from densitypy.project_utils.configuration_parser import parse_configuration_file
        # Lets show the configuration file selected content
        # Split JSON config into sections
        json_config = q.user.configuration_file
        print(f'{json_config = }')
        project_settings = json_config['projectsettings']
        grid_settings = json_config['gridsettings']
        charge_migration_settings = json_config['chargemigrationsettings']
        pump_settings = json_config['pumppulsessettings']
        probe_settings = json_config['probepulsessettings']
        charge_migration_ft_settings = json_config['chargemigrationftsettings']

        # Project Settings
        project_name = project_settings['projectname']
        xyz_geometry = project_settings['xyzmoleculegeometry']
        number_of_states = project_settings['numberofstates']
        list_of_orbitals = project_settings['listofactiveorbitals']
        molcas_output_directory = project_settings['molcasoutputdirectory']
        experiment_directory = project_settings['experimentdirectory']

        # Grid Settings
        nx = grid_settings['numberofpointsxaxis']
        ny = grid_settings['numberofpointsyaxis']
        nz = grid_settings['numberofpointszaxis']
        xmin = grid_settings['xmin']
        xmax = grid_settings['xmax']
        ymin = grid_settings['ymin']
        ymax = grid_settings['ymax']
        zmin = grid_settings['zmin']
        zmax = grid_settings['zmax']
        step_size = grid_settings['stepsize']
        boundary = grid_settings['boundary']

        # Charge Migration Parameters
        field_file = charge_migration_settings['fieldfile']
        number_of_times = charge_migration_settings['numberoftimes']
        min_time = charge_migration_settings['mintime']
        max_time = charge_migration_settings['maxtime']
        bath_temperature = charge_migration_settings['bathtemperature']
        dephasing_factor = charge_migration_settings['dephasingfactor']
        relaxation_factor = charge_migration_settings['relaxationfactor']

        # Pump Settings
        type_of_pulse_pump = pump_settings['typeofpulse']
        start_time = pump_settings['starttime']
        pump_central_frequency = pump_settings['pumpcentralfrequency']
        pump_periods = pump_settings['pumpperiods']
        pump_phase = pump_settings['pumpphase']
        pump_intensity = pump_settings['pumpintensity']
        pump_polarization = pump_settings['pumppolarization']

        # Probe Settings
        type_of_pulse_probe = probe_settings['typeofpulse']
        time_delay_start = probe_settings['timedelaystart']
        time_delay_stop = probe_settings['timedelaystop']
        number_of_pp = probe_settings['numberofpp']
        probe_central_frequency = probe_settings['probecentralfrequency']
        probe_periods = probe_settings['probeperiods']
        probe_phase = probe_settings['probephase']
        probe_intensity = probe_settings['probeintensity']
        probe_polarization = probe_settings['probepolarization']

        # Charge Migration FT Parameters
        number_of_omegas = charge_migration_ft_settings['numberofomegas']
        min_omegas = charge_migration_ft_settings['minomega']
        max_omegas = charge_migration_ft_settings['maxomega']
        number_of_tau_omegas = charge_migration_ft_settings['numberoftauomegas']
        min_tau_omega = charge_migration_ft_settings['mintauomega']
        max_tau_omega = charge_migration_ft_settings['maxtauomega']
        ft_time_step = charge_migration_ft_settings['fttimestep']
        ft_width_step = charge_migration_ft_settings['ftwidthstep']
        Volume = step_size * step_size * step_size

        print("Unit of Volume = " + str(Volume))

        time_delay_range = generate_time_delays(number_of_pp, time_delay_start, time_delay_stop)

        print("number_of_pp = ", number_of_pp)
        print("time_delay_start = ", time_delay_start)


# On Click Events


@on('ALM_studies_directory')
async def update_study_directory(q: Q):
    q.user.recent_opened_study_directory = q.args['ALM_studies_directory']
    # Update the form
    await load_page_recipe_with_update_models(q, homepage)


@on('ALM_configuration_file')
async def update_with_configuration_context(q: Q):
    # Show the configuration file selected
    q.user.recent_opened_configuration_file = q.args['ALM_configuration_file']
    # Load the configuration file
    from densitypy.project_utils.configuration_parser import parse_configuration_file
    q.user.configuration_file = parse_configuration_file(
        f'/home/ruben/PycharmProjects/DensityPy/Studies/{q.user.recent_opened_study_directory}/{q.user.recent_opened_configuration_file}')

    # Update the form
    await load_page_recipe_with_update_models(q, homepage)
