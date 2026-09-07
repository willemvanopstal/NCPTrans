from nicegui import ui
import sys
from pathlib import Path

# Import NCPTrans
sys.path.append(str(Path(__file__).parent.parent))
import NCPTrans
from nicegui import ui
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
import NCPTrans

ui.page_title = 'NCPTrans'

#ui.dark_mode().enable()

rows = []

columns = [
    {
        'name': 'key',
        'label': 'Property',
        'field': 'key',
        'sortable': True,
    },
    {
        'name': 'value',
        'label': 'Value',
        'field': 'value',
        'sortable': True,
    },
]

with ui.column().classes(
    'items-center w-full q-pa-xl'
):

    ui.label('NCPTrans').classes('text-h4')

    ui.label(
        'Dutch Continental Shelf Coordinate Transformations'
    ).classes('text-subtitle1')

    #
    # INPUT CARD
    #

    with ui.card().classes('w-full'):

        with ui.row().classes('w-full'):

            source_xy = ui.select(
                ['RD', 'ETRS89'],
                value='RD',
                label='Source XY',
            ).classes('col')

            source_h = ui.select(
                ['NAP', 'LAT', 'ETRS89'],
                value='NAP',
                label='Source Height',
            ).classes('col')

        with ui.row().classes('w-full'):

            x_input = ui.number(
                label='X / Easting / Longitude'
            ).classes('col')

            y_input = ui.number(
                label='Y / Northing / Latitude'
            ).classes('col')

        h_input = ui.number(
            label='Height (optional)'
        ).classes('w-full')

        status = ui.label('Ready')

        def compute():

            try:

                rows.clear()

                t = NCPTrans.Transformer()

                #
                # LOAD XY
                #

                if source_xy.value == 'RD':

                    t.load_rd_xy(
                        float(x_input.value),
                        float(y_input.value),
                    )

                    t.compute_etrs_latlon_from_rd_xy()

                else:

                    t.load_etrs_latlon(
                        float(y_input.value),  # lat
                        float(x_input.value),  # lon
                    )

                    t.compute_rd_xy_from_etrs_latlon()

                #
                # LOAD HEIGHT
                #

                if h_input.value is not None:

                    t.load_height(
                        float(h_input.value),
                        source_h.value,
                    )

                #
                # COMPUTE SEPARATIONS
                #

                t.compute_separations()

                #
                # FILL TABLE
                #

                for key, value in t.__dict__.items():

                    rows.append(
                        {
                            'key': key,
                            'value': str(value),
                        }
                    )
                table.rows[:] = rows
                table.update()

                status.set_text(
                    f'Success ({len(rows)} values)'
                )

            except Exception as ex:

                status.set_text(
                    f'Error: {ex}'
                )

        ui.button(
            'Compute',
            on_click=compute,
        ).props('color=primary')

        status

    #
    # OUTPUT CARD
    #

    with ui.card().classes('w-full'):

        ui.label('Transformer Contents').classes('text-h6')

        table = ui.table(
            columns=columns,
            rows=rows,
            row_key='key',
            pagination=50,
        ).classes('w-full')

ui.run(
    host='0.0.0.0',
    port=8080,
)
