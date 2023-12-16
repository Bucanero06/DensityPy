import asyncio
import os
from typing import Optional, List

from h2o_wave import Q
from h2o_wave import ui, on, data  # noqa F401

@on('global_notification_bar.dismissed')
async def on_global_notification_bar_dismissed(q: Q):
    # Delete the notification bar
    q.page['meta'].notification_bar = None
    await q.page.save()


async def update_index(q: Q):
    """
    Updates the indexes of data in the local directories as to have a truth when offering choices to the user.
    :param q:
    :return:
    """
    print("Not implemented yet (update_index)")
    pass

async def load_page_recipe_with_update_models(q: Q, page_callback_function):
    """
    Loads a page after updating the application models.

    :param q: Wave server Q object.
    :param page_callback_function: Callback function for page loading.
    :return: Result of the `page_callback_function` and `update_models` functions.

    """
    await page_callback_function(q)
    await q.page.save()
    return await asyncio.gather(page_callback_function(q), update_index(q))

async def push_notification_bar(q: Q, text: str, type: str = 'info', timeout: int = 5,
                                buttons: Optional[List[str]] = None, position: str = 'top-right',
                                events: Optional[List[str]] = None, name: Optional[str] = None) -> None:
    """
    Notification Bar Arguments:
        text
            The text displayed on the notification bar.
        type
            The icon and color of the notification bar. Defaults to 'info'. One of 'info', 'error', 'warning', 'success', 'danger', 'blocked'. See enum h2o_wave.ui.NotificationBarType.
        timeout
            How long the notification stays visible, in seconds. If set to -1, the notification has to be closed manually. Defaults to 5.
        buttons
            Specify one or more action buttons.
        position
            Specify the location of notification. Defaults to 'top-right'. One of 'top-right', 'bottom-right', 'bottom-center', 'bottom-left', 'top-left', 'top-center'. See enum h2o_wave.ui.NotificationBarPosition.
        events
            The events to capture on this notification bar. One of 'dismissed'.
        name
            An identifying name for this component.
    """

    # Create a notification bar
    q.page['meta'].notification_bar = ui.notification_bar(
        text=text,
        type=type,
        timeout=timeout,
        buttons=buttons,
        position=position,
        events=events,
        name=name,
    )

async def stream_message(q:Q, page_name:str, message:str, delay:float = 0.3):
    stream = ''
    q.page[page_name].data += [stream, False]
    # Show the "Stop generating" button
    q.page[page_name].generating = True
    for w in message.split():
        await asyncio.sleep(delay)
        stream += w + ' '
        q.page[page_name].data[-1] = [stream, False]
        await q.page.save()
    # Hide the "Stop generating" button
    q.page[page_name].generating = False
    await q.page.save()


def remove_card(q: Q, name: str) -> None:
    """
    Remove a specific card from the page based on its name.
    """
    if hasattr(q.client, 'cards') and name in q.client.cards:
        del q.page[name]
        q.client.cards.remove(name)


def add_card(q: Q, name, card) -> None:

    q.client.cards.add(name)
    q.page[name] = card


def clear_cards(q: Q, ignore: Optional[List[str]] = None) -> None:
    """
    Clear cards from the page except those listed in 'ignore'.
    """
    print("Clearing cards")

    if not ignore:
        ignore = []

    # If no cards are present, simply return
    if not hasattr(q.client, 'cards') or not q.client.cards:
        print("No cards")
        return

    # Remove cards not in the ignore list
    for card_name in q.client.cards.copy():
        if card_name not in ignore:
            # del q.page[card_name]
            # q.client.cards.remove(card_name)
            remove_card(q, card_name)



def load_env_file(env_path: str = '.env'):
    """Load environment variables from a file."""
    if os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                if line.startswith('#') or line == '\n':
                    continue
                key, value = line.strip().split('=', 1)
                os.environ[key] = value
    else:
        print(f"File does not exist: {env_path}")
        raise FileNotFoundError(f"File does not exist: {env_path}")


def dynamic_tall_series_stat_card(box, title, main_stat, aux_stat, chart_data, additional_params=None):
    card = ui.tall_series_stat_card(
        box=box,
        title=title,
        value='=${{intl main_val minimum_fraction_digits=2 maximum_fraction_digits=2}}',
        aux_value='={{intl aux_val style="percent" minimum_fraction_digits=1 maximum_fraction_digits=1}}',
        data={
            'main_val': main_stat,
            'aux_val': aux_stat,
        },

        plot_data=chart_data,
        **additional_params  # For any additional customizations
    )
    return card
