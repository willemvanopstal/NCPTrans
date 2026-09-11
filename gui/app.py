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
  Marine Regions / VLIZ-WFS aanloopt).

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
NCP_EEZ_WFS_URL = (
    "https://geo.vliz.be/geoserver/MarineRegions/wfs"
    "?service=WFS&version=1.0.0&request=GetFeature"
    f"&typeName=eez&cql_filter=mrgid={NCP_EEZ_MRGID}&outputFormat=application/json"
)
NCP_EEZ_CACHE_FILE = STATIC_DIR / '.ncp_eez_cache.geojson'


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
    def _send_bytes(self, body: bytes, content_type: str, status: int = 200) -> None:
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, obj, status: int = 200) -> None:
        self._send_bytes(json.dumps(obj).encode('utf-8'), 'application/json; charset=utf-8', status)

    def _send_file(self, path: Path, content_type: str) -> None:
        try:
            body = path.read_bytes()
        except FileNotFoundError:
            self.send_error(404, 'Not found')
            return
        self._send_bytes(body, content_type)

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
            t = run_full_transform(x, y, h, src_xy, src_h)
        except Exception as exc:
            self._send_json({'error': str(exc)})
            return

        self._send_json({
            'rd_x': t.rd_x,
            'rd_y': t.rd_y,
            'etrs_lon': t.etrs_lon,
            'etrs_lat': t.etrs_lat,
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
