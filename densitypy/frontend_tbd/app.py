import asyncio

from h2o_wave import main, Q, app, ui, on, data, handle_on, run_on  # noqa F401

from densitypy.frontend_tbd.HomePage import homepage
from densitypy.frontend_tbd.util import add_card, clear_cards, load_page_recipe_with_update_models


async def init(q: Q) -> None:
    """
    Q Page Meta (meta_card) Arguments:
        box
            A string indicating how to place this component on the page.
        title
            The title of the page.
        refresh
            Refresh rate in seconds. A value of 0 turns off live-updates. Values != 0 are currently ignored (reserved for future use).
        notification
            Display a desktop notification.
        notification_bar
            Display an in-app notification bar.
        redirect
            Redirect the page to a new URL.
        icon
            Shortcut icon path. Preferably a .png file (.ico files may not work in mobile browsers). Not supported in Safari.
        layouts
            The layouts supported by this page.
        dialog
            Display a dialog on the page.
        side_panel
            Display a side panel on the page.
        theme
            Specify the name of the theme (color scheme) to use on this page. One of 'light', 'neon' or 'h2o-dark'.
        themes
            Themes (color schemes) that define color used in the app.
        tracker
            Configure a tracker for the page (for web analytics).
        scripts
            External Javascript files to load into the page.
        script
            Javascript code to execute on this page.
        stylesheet
            CSS stylesheet to be applied to this page.
        stylesheets
            External CSS files to load into the page.
        commands
            Contextual menu commands for this component.
    """
    # Static Business Website
    # index_file = open('static/html/index.html', 'r').read()

    q.page['meta'] = ui.meta_card(box='',
                                  title='8CoedSports',
                                  layouts=[ui.layout(breakpoint='xs', min_height='100vh', zones=[
                                      ui.zone('main', size='1', direction=ui.ZoneDirection.ROW, zones=[
                                          ui.zone('sidebar', size='208px'),
                                          ui.zone('body', zones=[
                                              ui.zone('header'),
                                              ui.zone('content', zones=[
                                                  # Specify various zones and use the one that is currently needed. Empty zones are ignored.
                                                  ui.zone('first_context', size='0 0 1 4',
                                                          direction=ui.ZoneDirection.ROW,
                                                          zones=[
                                                              ui.zone('first_context_1', size='1 4 0 0'),
                                                              ui.zone('first_context_2', size='1 4 0 0'),
                                                              ui.zone('first_context_3', size='1 4 0 0'),
                                                              ui.zone('first_context_4', size='1 4 0 0'),
                                                          ]),
                                                  ui.zone('second_context', size='0 0 1 4',
                                                          direction=ui.ZoneDirection.ROW,
                                                          zones=[
                                                              ui.zone('second_context_1', size='1 4 0 0'),
                                                              ui.zone('second_context_2', size='1 4 0 0'),
                                                              ui.zone('second_context_3', size='1 4 0 0',
                                                                      direction=ui.ZoneDirection.ROW,
                                                                      zones=[
                                                                          ui.zone('second_context_3_1', size='1 4 0 0'),
                                                                          ui.zone('second_context_3_2', size='1 4 0 0')
                                                                      ]),
                                                          ]),
                                                  ui.zone('details', size='4 4 4 4'),
                                                  ui.zone('horizontal', size='1', direction=ui.ZoneDirection.ROW),
                                                  ui.zone('centered', size='1 1 1 1', align='center'),
                                                  ui.zone('vertical', size='1'),
                                                  ui.zone('grid', direction=ui.ZoneDirection.ROW, wrap='stretch',
                                                          justify='center'),

                                              ]),
                                          ]),
                                      ]),
                                      ui.zone('footer', size='0 1 0 0', direction=ui.ZoneDirection.ROW),
                                  ]),
                                           ],
                                  themes=[
                                      ui.theme(
                                          name='my-awesome-theme',
                                          primary='#8C1B11',  # Header and Sidebaer - Color Light Red
                                          text='#000000',  #
                                          card='#ffffff',
                                          page='#F2F2F2',
                                          # page='#D91A1A',

                                      )
                                  ],
                                  theme='my-awesome-theme'

                                  )
    # Sidebar should be initialized only with non-authenticated pages content only!!!
    add_card(q, 'Application_Sidebar', ui.nav_card(
        box='sidebar', color='primary', title='Demo Admin App',
        subtitle="The team sports you loved as a kid... now all grown up.",
        # Local Image
        image='https://8coedsports.com/wp-content/uploads/2019/12/8coed-logo-color.svg',
        items=[]
    ))
    q.page['footer'] = ui.footer_card(box='footer',
                                      caption='©2023 8CoedSports & Carbonyl LLC. Partnership. All rights reserved.', )

    q.client.initialized = False


async def initialize_client(q: Q):
    q.client.cards = set()
    await init(q)
    q.client.initialized = True
    q.client.token = None  # Initially, no token is set.


async def render_hidden_content(q: Q, context: dict = None):
    """
    Render pages content e.g. homepage or other pages get added here
    """

    # First clear all cards
    await q.run(clear_cards, q, ignore=['Application_Sidebar'])

    # Then add the sidebar
    add_card(q, 'Application_Sidebar', ui.nav_card(
        box='sidebar', color='primary', title='DensityPy',
        subtitle="hi I am so young I don't have a description",
        value=f'#{q.args["#"]}' if q.args['#'] else '#homepage',
        # image='https://wave.h2o.ai/img/h2o-logo.svg', items=[]
        # Loading from local file rather than url
        image='https://avatars.githubusercontent.com/u/60953006?v=4', items=[
            ui.nav_group('Main', items=[
                ui.nav_item(name='#homepage', label='Home', icon='Home'),
            ]),
            ui.nav_group('Docs', items=[
                ui.nav_item(name='#admin_userdocs_link', label='Usage Documentation', icon='TextDocumentShared'),
                ui.nav_item(name='#admin_devdocs_link', label='Codebase Documentation', icon='TextDocumentEdit'),
            ]),
        ],
    ))
    # ui.link(label='Dev Docs', path='https://8coedsports.com/docs/_build/html/index.html', target='_blank')
    if q.args['#'] == 'homepage':
        await load_page_recipe_with_update_models(q, homepage)



@app('/')
async def serve(q: Q):
    """Main application handler."""
    print("Serving")
    if not q.client.initialized:
        print("Initializing")
        await initialize_client(q)

    await render_hidden_content(q),

    await q.page.save()

    await run_on(q)
