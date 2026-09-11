#!/usr/bin/env python3
"""
NCPTrans - webapp (v4, geen NiceGUI)
=====================================

Kleine server die uitsluitend de Python-standaardbibliotheek gebruikt
(geen NiceGUI, geen Flask/FastAPI) en:

- de losstaande pagina gui/index.html serveert (vanilla HTML/CSS/JS +
  Leaflet via CDN) op '/',
- een JSON-API biedt op '/api/transform' die de echte berekeningen
  van NCPTrans.py uitvoert (dezelfde grid-gebaseerde RDNAPTRANS-logica
  als voorheen, nu client-side alleen via HTTP aangeroepen),
- de exacte NCP/EEZ-grens server-side ophaalt en doorgeeft op
  '/api/ncp-boundary' (zodat de browser niet zelf tegen CORS van de
  Marine Regions / VLIZ-WFS aanloopt),
- het GeoTIFF-bestand latnap2018.tif (LAT-NAP-verschilraster, gebouwd
  met build_lat_nap_raster.py; indien aanwezig naast NCPTrans.py)
  serveert op '/api/latnap-raster', zodat de browser het zelf met
  georaster/georaster-layer-for-leaflet kan inlezen en tonen,
- daarnaast bij elke transformatie ook de UTM 31N/32N-coordinaten van het
  punt meelevert (standaard UTM-projectie op GRS80, los van RDNAPTRANS),
  plus een ED50-weergave (eenvoudige 3-parameter datumshift) en twee
  "verwisselingsfouten" (RD/Bessel vs ETRS89, en ED50 vs ETRS89) - de
  afstand (in UTM-meters) die je ernaast zou zitten als je een coordinaat
  per ongeluk in het verkeerde datum zou interpreteren. Als invoersysteem
  worden naast RD en ETRS89 ook ED50 en UTM 31N/32N geaccepteerd.

Starten:
    python gui/app.py [poort]      (standaard poort 8080)

Open daarna http://localhost:8080 in de browser.

Vereisten:
    Alleen wat NCPTrans.py zelf al nodig heeft (o.a. numpy, via
    helpers.py) en de sqlite-grid-database (rdnaptrans_grids.sqlite)
    naast constants.py. Zie de hoofd-README van NCPTrans.
"""

from __future__ import annotations

import http.server
import json
import math
import socketserver
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import NCPTrans  # noqa: E402

STATIC_DIR = Path(__file__).resolve().parent
INDEX_FILE = STATIC_DIR / 'index.html'

NCP_EEZ_MRGID = 5668
# We gebruiken hier bewust de "eez_boundaries"-lijnenlaag (de onderhandelde
# maritieme grenzen tussen landen) i.p.v. de volledige "eez"-polygoonlaag.
# Die laatste volgt namelijk ook de kustlijn in detail, wat op de kaart
# rommelig oogt; deze laag geeft alleen de rechte(re) lijnstukken op zee.
NCP_EEZ_WFS_URL = (
    "https://geo.vliz.be/geoserver/MarineRegions/wfs"
    "?service=WFS&version=1.0.0&request=GetFeature"
    f"&typeName=eez_boundaries&cql_filter=mrgid_eez1={NCP_EEZ_MRGID}+OR+mrgid_eez2={NCP_EEZ_MRGID}"
    "&outputFormat=application/json"
)
NCP_EEZ_CACHE_FILE = STATIC_DIR / '.ncp_eez_boundaries_cache.geojson'

# GeoTIFF met het LAT-NAP-verschilraster (zie build_lat_nap_raster.py),
# verwacht naast NCPTrans.py (dus in dezelfde map als constants.py / helpers.py).
LATNAP_TIF_FILE = REPO_ROOT / 'latnap2018.tif'

# -- UTM (standaard Snyder-formules op GRS80, k0=0.9996) -----------------
# Los van RDNAPTRANS: puur voor de UTM 31N/32N-kolommen in de uitvoer.
_UTM_A = 6378137.0
_UTM_F = 1 / 298.257222101  # GRS80
_UTM_E2 = _UTM_F * (2 - _UTM_F)
_UTM_EP2 = _UTM_E2 / (1 - _UTM_E2)
_UTM_K0 = 0.9996


def utm_forward(lat_deg: float, lon_deg: float, zone: int):
    """RD/ETRS89-onafhankelijke standaard UTM-projectie (GRS80, noordelijk halfrond)."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    lon0 = math.radians(zone * 6 - 183)
    e2, ep2, a, k0 = _UTM_E2, _UTM_EP2, _UTM_A, _UTM_K0

    n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    tt = math.tan(lat) ** 2
    cc = ep2 * math.cos(lat) ** 2
    aa = math.cos(lat) * (lon - lon0)

    m = a * (
        (1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * lat
        - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * math.sin(2 * lat)
        + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * math.sin(4 * lat)
        - (35 * e2 ** 3 / 3072) * math.sin(6 * lat)
    )

    easting = k0 * n * (
        aa + (1 - tt + cc) * aa ** 3 / 6
        + (5 - 18 * tt + tt ** 2 + 72 * cc - 58 * ep2) * aa ** 5 / 120
    ) + 500000.0

    northing = k0 * (
        m + n * math.tan(lat) * (
            aa ** 2 / 2
            + (5 - tt + 9 * cc + 4 * cc ** 2) * aa ** 4 / 24
            + (61 - 58 * tt + tt ** 2 + 600 * cc - 330 * ep2) * aa ** 6 / 720
        )
    )
    return easting, northing


def utm_inverse(easting: float, northing: float, zone: int):
    """Inverse van utm_forward (standaard Snyder-formules, GRS80)."""
    a, e2, ep2, k0 = _UTM_A, _UTM_E2, _UTM_EP2, _UTM_K0
    lon0 = math.radians(zone * 6 - 183)
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))

    m = northing / k0
    mu = m / (a * (1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256))

    phi1 = (
        mu
        + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * math.sin(2 * mu)
        + (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32) * math.sin(4 * mu)
        + (151 * e1 ** 3 / 96) * math.sin(6 * mu)
        + (1097 * e1 ** 4 / 512) * math.sin(8 * mu)
    )

    c1 = ep2 * math.cos(phi1) ** 2
    t1 = math.tan(phi1) ** 2
    n1 = a / math.sqrt(1 - e2 * math.sin(phi1) ** 2)
    r1 = a * (1 - e2) / (1 - e2 * math.sin(phi1) ** 2) ** 1.5
    d = (easting - 500000.0) / (n1 * k0)

    lat = phi1 - (n1 * math.tan(phi1) / r1) * (
        d ** 2 / 2
        - (5 + 3 * t1 + 10 * c1 - 4 * c1 ** 2 - 9 * ep2) * d ** 4 / 24
        + (61 + 90 * t1 + 298 * c1 + 45 * t1 ** 2 - 252 * ep2 - 3 * c1 ** 2) * d ** 6 / 720
    )
    lon = lon0 + (
        d - (1 + 2 * t1 + c1) * d ** 3 / 6
        + (5 - 2 * c1 + 28 * t1 - 3 * c1 ** 2 + 8 * ep2 + 24 * t1 ** 2) * d ** 5 / 120
    ) / math.cos(phi1)

    return math.degrees(lat), math.degrees(lon)


# -- ED50 <-> ETRS89 (eenvoudige 3-parameter geocentrische verschuiving) --
# LET OP: dit is een benadering. ED50 (International 1924/Hayford-ellipsoide)
# gebruikt hier de klassieke, over heel West-Europa gemiddelde verschuiving
# (vergelijkbaar met EPSG-transformatie 1311), GEEN Nederland-specifieke
# fit. Nauwkeurigheid ligt in de orde van enkele meters, niet geschikt voor
# nauwkeurigheidseisen waarvoor je normaal ED50 zou raadplegen.
_ED50_A = 6378388.0
_ED50_F = 1 / 297.0
_ED50_E2 = _ED50_F * (2 - _ED50_F)
_ED50_TO_ETRS89_DX = -87.0
_ED50_TO_ETRS89_DY = -98.0
_ED50_TO_ETRS89_DZ = -121.0


def _geodetic_to_geocentric(lat_deg, lon_deg, h, a, e2):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    x = (n + h) * math.cos(lat) * math.cos(lon)
    y = (n + h) * math.cos(lat) * math.sin(lon)
    z = (n * (1 - e2) + h) * math.sin(lat)
    return x, y, z


def _geocentric_to_geodetic(x, y, z, a, e2, iterations=8):
    lon = math.atan2(y, x)
    p = math.sqrt(x * x + y * y)
    lat = math.atan2(z, p * (1 - e2))
    for _ in range(iterations):
        n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - n
        lat = math.atan2(z, p * (1 - e2 * n / (n + h)))
    n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - n
    return math.degrees(lat), math.degrees(lon), h


def ed50_to_etrs89(lat_deg: float, lon_deg: float, h: float = 0.0):
    x, y, z = _geodetic_to_geocentric(lat_deg, lon_deg, h, _ED50_A, _ED50_E2)
    x2 = x + _ED50_TO_ETRS89_DX
    y2 = y + _ED50_TO_ETRS89_DY
    z2 = z + _ED50_TO_ETRS89_DZ
    lat2, lon2, h2 = _geocentric_to_geodetic(x2, y2, z2, _UTM_A, _UTM_E2)
    return lat2, lon2, h2


def etrs89_to_ed50(lat_deg: float, lon_deg: float, h: float = 0.0):
    x, y, z = _geodetic_to_geocentric(lat_deg, lon_deg, h, _UTM_A, _UTM_E2)
    x2 = x - _ED50_TO_ETRS89_DX
    y2 = y - _ED50_TO_ETRS89_DY
    z2 = z - _ED50_TO_ETRS89_DZ
    lat2, lon2, h2 = _geocentric_to_geodetic(x2, y2, z2, _ED50_A, _ED50_E2)
    return lat2, lon2, h2


def run_full_transform(x: float, y: float, h, src_xy: str, src_h):
    """Bereken alle horizontale/verticale representaties van 1 punt.

    Bouwt zelf een NCPTrans.Transformer op (dezelfde stappen als
    NCPTrans.transform()) zodat we ook alle tussenwaarden overhouden
    voor de debug-uitvoer.
    """
    t = NCPTrans.Transformer()

    if src_xy == 'RD':
        t.load_rd_xy(x, y)
    elif src_xy == 'ETRS89':
        t.load_etrs_latlon(y, x)  # load_etrs_latlon(lat, lon)
    else:
        raise ValueError(f'Onbekend src_xy: {src_xy}')

    if h is not None and src_h:
        t.load_height(h, src_h)

    if src_xy == 'RD':
        t.compute_etrs_latlon_from_rd_xy()
    else:
        t.compute_rd_xy_from_etrs_latlon()

    t.compute_separations(include_lat=True)
    return t


def transformer_to_jsonable(t) -> dict:
    out = {}
    for key, value in vars(t).items():
        if key == 'cache':
            continue
        if value is None or isinstance(value, (int, float, str, bool)):
            out[key] = value
        elif isinstance(value, tuple):
            out[key] = [v if isinstance(v, (int, float, str, bool)) or v is None else str(v) for v in value]
        else:
            out[key] = str(value)
    return out


def fetch_ncp_boundary():
    if NCP_EEZ_CACHE_FILE.exists():
        try:
            return json.loads(NCP_EEZ_CACHE_FILE.read_text())
        except Exception:
            pass
    try:
        req = urllib.request.Request(NCP_EEZ_WFS_URL, headers={'User-Agent': 'NCPTrans-gui'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.load(resp)
        try:
            NCP_EEZ_CACHE_FILE.write_text(json.dumps(data))
        except Exception:
            pass
        return data
    except Exception:
        return None


class Handler(http.server.BaseHTTPRequestHandler):
    server_version = 'NCPTrans/4.0'

    # -- helpers ----------------------------------------------------
    def _send_bytes(self, body: bytes, content_type: str, status: int = 200, cache: bool = False) -> None:
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        if cache:
            # Statische, zelden wijzigende bestanden (zoals het raster) mogen
            # door de browser onthouden worden; scheelt een herhaalde download.
            self.send_header('Cache-Control', 'public, max-age=86400')
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, obj, status: int = 200) -> None:
        self._send_bytes(json.dumps(obj).encode('utf-8'), 'application/json; charset=utf-8', status)

    def _send_file(self, path: Path, content_type: str, cache: bool = False) -> None:
        try:
            body = path.read_bytes()
        except FileNotFoundError:
            self.send_error(404, 'Not found')
            return
        self._send_bytes(body, content_type, cache=cache)

    # -- routes -------------------------------------------------------
    def do_GET(self) -> None:  # noqa: N802
        if self.path in ('/', '/index.html'):
            self._send_file(INDEX_FILE, 'text/html; charset=utf-8')
        elif self.path == '/api/ncp-boundary':
            data = fetch_ncp_boundary()
            if data is None:
                self._send_json({'error': 'kon NCP/EEZ-grens niet ophalen (geen internet of WFS onbereikbaar)'}, 502)
            else:
                self._send_json(data)
        elif self.path == '/api/latnap-raster':
            if LATNAP_TIF_FILE.exists():
                self._send_file(LATNAP_TIF_FILE, 'image/tiff', cache=True)
            else:
                self._send_json({'error': f'{LATNAP_TIF_FILE.name} niet gevonden naast NCPTrans.py'}, 404)
        else:
            self.send_error(404, 'Not found')

    def do_POST(self) -> None:  # noqa: N802
        if self.path != '/api/transform':
            self.send_error(404, 'Not found')
            return

        length = int(self.headers.get('Content-Length', 0) or 0)
        try:
            payload = json.loads(self.rfile.read(length) or b'{}')
        except Exception:
            self._send_json({'error': 'ongeldige JSON in request body'}, 400)
            return

        try:
            x = float(payload['x'])
            y = float(payload['y'])
            src_xy = payload['src_xy']
            h = payload.get('h')
            h = float(h) if h is not None else None
            src_h = payload.get('src_h') or None

            # ED50 en UTM zijn geen systemen die NCPTrans zelf kent; die
            # rekenen we hier eerst om naar ETRS89 lon/lat, en laten
            # run_full_transform() daarna gewoon zijn normale werk doen.
            if src_xy == 'ED50':
                lat, lon, _h = ed50_to_etrs89(y, x)
                x, y, src_xy = lon, lat, 'ETRS89'
            elif src_xy in ('UTM31N', 'UTM32N'):
                zone = 31 if src_xy == 'UTM31N' else 32
                lat, lon = utm_inverse(x, y, zone)
                x, y, src_xy = lon, lat, 'ETRS89'

            t = run_full_transform(x, y, h, src_xy, src_h)
        except Exception as exc:
            self._send_json({'error': str(exc)})
            return

        utm31_e = utm31_n = utm32_e = utm32_n = None
        utm_native_zone = None
        ed50_lat = ed50_lon = None
        rd_etrs_mixup_distance = None
        rd_etrs_mixup_bearing = None
        ed50_etrs_mixup_distance = None
        ed50_etrs_mixup_bearing = None

        if t.etrs_lat is not None and t.etrs_lon is not None:
            utm31_e, utm31_n = utm_forward(t.etrs_lat, t.etrs_lon, 31)
            utm32_e, utm32_n = utm_forward(t.etrs_lat, t.etrs_lon, 32)
            utm_native_zone = 31 if t.etrs_lon < 6.0 else 32
            native_e, native_n = (utm31_e, utm31_n) if utm_native_zone == 31 else (utm32_e, utm32_n)

            ed50_lat, ed50_lon, _h = etrs89_to_ed50(t.etrs_lat, t.etrs_lon)

            # Verwisselingsfouten: het punt in een ANDER datum uitgedrukt,
            # maar diens lat/lon-getallen naief (zonder omrekening) door
            # dezelfde UTM-projectie gehaald, vergeleken met de echte
            # ETRS89-UTM-coordinaat - dat simuleert precies de fout die je
            # maakt als je een coordinaat uit het verkeerde datum leest.
            # De richting (bearing, 0-360 gr. vanaf het noorden, met de klok
            # mee) geeft aan waar de ECHTE ETRS89-positie ligt t.o.v. de
            # foutief geinterpreteerde positie.
            if t.bessel_lat_deg is not None and t.bessel_lon_deg is not None:
                bessel_e, bessel_n = utm_forward(t.bessel_lat_deg, t.bessel_lon_deg, utm_native_zone)
                d_e, d_n = native_e - bessel_e, native_n - bessel_n
                rd_etrs_mixup_distance = math.hypot(d_e, d_n)
                rd_etrs_mixup_bearing = math.degrees(math.atan2(d_e, d_n)) % 360

            ed50_e, ed50_n = utm_forward(ed50_lat, ed50_lon, utm_native_zone)
            d_e, d_n = native_e - ed50_e, native_n - ed50_n
            ed50_etrs_mixup_distance = math.hypot(d_e, d_n)
            ed50_etrs_mixup_bearing = math.degrees(math.atan2(d_e, d_n)) % 360

        self._send_json({
            'rd_x': t.rd_x,
            'rd_y': t.rd_y,
            'rd_out_of_bounds': bool(t.rd_correction_out_of_bounds),
            'etrs_lon': t.etrs_lon,
            'etrs_lat': t.etrs_lat,
            'ed50_lat': ed50_lat,
            'ed50_lon': ed50_lon,
            'utm31_easting': utm31_e,
            'utm31_northing': utm31_n,
            'utm32_easting': utm32_e,
            'utm32_northing': utm32_n,
            'utm_native_zone': utm_native_zone,
            'rd_etrs_mixup_distance': rd_etrs_mixup_distance,
            'rd_etrs_mixup_bearing': rd_etrs_mixup_bearing,
            'ed50_etrs_mixup_distance': ed50_etrs_mixup_distance,
            'ed50_etrs_mixup_bearing': ed50_etrs_mixup_bearing,
            'nap_height': t.nap_height,
            'lat_height': t.lat_height,
            'ellipsoidal_height': t.ellipsoidal_height,
            'geoid_ellipsoid_separation': t.geoid_ellipsoid_separation,
            'lat_ellipsoid_separation': t.lat_ellipsoid_separation,
            'geoid_lat_separation': t.geoid_lat_separation,
            'debug': transformer_to_jsonable(t),
        })

    def log_message(self, fmt, *args) -> None:  # rustige console
        sys.stderr.write('%s - %s\n' % (self.address_string(), fmt % args))


def main() -> None:
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
    with socketserver.ThreadingTCPServer(('0.0.0.0', port), Handler) as httpd:
        print(f'NCPTrans draait op http://localhost:{port}  (Ctrl+C om te stoppen)')
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == '__main__':
    main()