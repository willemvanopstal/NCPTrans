#!/usr/bin/env python3
"""
build_lat_nap_raster.py
========================

Bouwt een LAT-NAP-verschilraster (GeoTIFF) uit de NLGEO2018- en
NLLAT2018-correctiegrids die al in rdnaptrans_grids.sqlite zitten
(dezelfde database die NCPTrans.py / de webapp gebruikt).

Waarde per pixel = NAP-hoogte minus LAT-hoogte op dat punt, exact
volgens dezelfde bilineaire grid-interpolatie als NCPTrans.py zelf
gebruikt in Transformer.compute_separations():

    geoid_lat_separation = geoid_ellipsoid_separation - lat_ellipsoid_separation

Dit is dus 1-op-1 consistent met wat de webapp al toont onder
"NAP - LAT" - alleen nu als raster over de hele kaart, in plaats van
losse puntwaarden.

Het rooster (extent + resolutie) is gelijk aan het NLLAT2018-grid
zoals NCPTrans.py dat operationeel gebruikt (zie constants.py:
NLLAT2018), NIET per se gelijk aan de volledige extent van je eigen
nllat2018.tif-bestand (dat kan een groter, ruwer gepubliceerd gebied
beslaan dan het gebied waarbinnen NCPTrans daadwerkelijk geldige
waarden teruggeeft). Zo sluit dit bestand exact aan op wat de app
zelf ook echt berekent, in plaats van op de losse rasterbestanden.

Vereisten:
    - rdnaptrans_grids.sqlite naast constants.py (dezelfde die de
      webapp gebruikt).
    - numpy en tifffile:
        pip install numpy tifffile --break-system-packages

Gebruik:
    python build_lat_nap_raster.py [output.tif]
    (standaard: latnap2018.tif, in dezelfde map als dit script)

Het resultaat kun je op dezelfde manier in gui/app.py serveren en in
gui/index.html laden als het bestaande NLLAT2018-raster (zelfde
GeoTIFF-tags: EPSG:9067, dus dezelfde proj4-registratie werkt al).
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import numpy as np
import tifffile

REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT))

import constants  # noqa: E402
import NCPTrans  # noqa: E402


def load_table_into_cache(table: str, value_col: str) -> None:
    """Leest een hele grid-tabel in één keer in het geheugen en vult
    Transformer.cache, zodat latere aanroepen van get_val_at_line() /
    get_corrections_at_index() nooit meer de sqlite-database hoeven
    te raken. Zonder deze truc zou dit script voor elk van de
    honderdduizenden roosterpunten losse SQL-queries doen (te traag
    voor een heel grid); dit is puur een snelheidsoptimalisatie en
    verandert niets aan de uitkomst."""
    db_path = Path(constants.RESOURCES) / constants.DB_FILE
    cache = NCPTrans.Transformer.cache.setdefault(table, {})
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        cur.execute(f'SELECT line, {value_col} FROM {table}')
        for line, value in cur.fetchall():
            cache[line] = value


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / 'latnap2018.tif'

    print('Grids inladen in geheugen (eenmalig, scheelt honderdduizenden losse queries)...')
    load_table_into_cache('nlgeo2018', 'geoid')
    load_table_into_cache('nllat2018', 'lat_ellipsoid')

    # Zelfde rooster als NCPTrans.py operationeel gebruikt voor LAT.
    grid = constants.NLLAT2018
    n_rows = int(round(grid['PHI_N']))
    n_cols = int(round(grid['LMBDA_N']))
    phi_min, phi_delta = grid['PHI_MIN'], grid['PHI_DELTA']
    lmbda_min, lmbda_delta = grid['LMBDA_MIN'], grid['LMBDA_DELTA']

    print(
        f'Rooster: {n_rows} x {n_cols} punten '
        f'({phi_min}-{grid["PHI_MAX"]} gr. NB, {lmbda_min}-{grid["LMBDA_MAX"]} gr. OL)'
    )

    t = NCPTrans.Transformer()
    data = np.full((n_rows, n_cols), np.nan, dtype=np.float32)
    skipped = 0

    for row in range(n_rows):
        phi = phi_min + row * phi_delta
        if row % 100 == 0:
            print(f'  rij {row}/{n_rows}...')
        for col in range(n_cols):
            lmbda = lmbda_min + col * lmbda_delta
            try:
                geoid_sep = t.get_geoid_ellipsoid_separation(phi, lmbda, constants.NLGEO2018)
                lat_sep = t.get_lat_ellipsoid_separation(phi, lmbda, constants.NLLAT2018)
            except TypeError:
                # Uiterste rand van het rooster kan door een floating-point
                # afrondingsverschil net buiten de grid-index vallen; dat
                # raakt hooguit een handvol randpixels en wordt hier als
                # nodata behandeld (net als een echte out-of-bounds None).
                skipped += 1
                continue
            if geoid_sep is None or lat_sep is None:
                skipped += 1
                continue
            # De grid-rekenrichting van NCPTrans loopt van zuid naar
            # noord (phi_min -> phi_max); een TIFF-rij 0 hoort bij het
            # noorden. Vandaar de omkering hier.
            tiff_row = n_rows - 1 - row
            data[tiff_row, col] = geoid_sep - lat_sep

    valid = np.count_nonzero(~np.isnan(data))
    print(f'Klaar met rekenen: {valid} geldige punten, {skipped} buiten bereik/overgeslagen.')

    nodata = -32768.0
    out = np.where(np.isnan(data), nodata, data).astype(np.float32)

    print(f'Wegschrijven naar {out_path} ...')
    tifffile.imwrite(
        out_path,
        out,
        photometric='minisblack',
        extratags=[
            # ModelPixelScaleTag (33550): lon-stap, lat-stap, 0
            (33550, 'd', 3, (lmbda_delta, phi_delta, 0.0), True),
            # ModelTiepointTag (33922): pixel (0,0,0) hoort bij model-
            # coordinaat (lmbda_min, phi_max, 0) - linkerbovenhoek.
            (33922, 'd', 6, (0.0, 0.0, 0.0, lmbda_min, grid['PHI_MAX'], 0.0), True),
            # GeoKeyDirectoryTag (34735): zelfde geografische CRS als het
            # bestaande NLLAT-bestand (EPSG:9067 / ETRF2000, ~ETRS89 lon-lat).
            (34735, 'H', 20, (
                1, 1, 1, 4,
                1024, 0, 1, 2,
                1025, 0, 1, 2,
                2048, 0, 1, 9067,
                4096, 0, 1, 7931,
            ), True),
            # GDAL_NODATA (42113): zelfde conventie als het bestaande bestand.
            (42113, 's', 6, '-32768', True),
        ],
    )
    print('Klaar.')


if __name__ == '__main__':
    main()