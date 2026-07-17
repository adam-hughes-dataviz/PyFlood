import json
import math
import os
import re
import time

from flask import Flask, jsonify, request, send_from_directory
import urllib.request

app = Flask(__name__, static_folder="static")

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)
STATION_CACHE = os.path.join(CACHE_DIR, "stations.json")
FLOW_CACHE_DIR = os.path.join(CACHE_DIR, "flows")
os.makedirs(FLOW_CACHE_DIR, exist_ok=True)

STATION_LIST_URL = (
    "https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info"
    "?station=*&format=json-object"
    "&fields=id,name,river,location,catchment-area,spatial-location"
)

_station_cache_ttl = 60 * 60 * 24 * 7  # 1 week


# ---------- coordinate parsing ----------

def parse_coords(text):
    text = text.strip()
    dms_re = re.compile(r"(\d+)\s*[°:]\s*(\d+)\s*['\u2032:]\s*([\d.]+)\s*[\"\u2033]?\s*([NSEWnsew])")
    matches = dms_re.findall(text)
    if len(matches) >= 2:
        lat = lon = None
        for deg, minutes, sec, direction in matches:
            val = float(deg) + float(minutes) / 60 + float(sec) / 3600
            direction = direction.upper()
            if direction in ("S", "W"):
                val = -val
            if direction in ("N", "S"):
                lat = val
            else:
                lon = val
        if lat is not None and lon is not None:
            return lat, lon

    dec_re = re.compile(r"(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)")
    m = dec_re.search(text)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None


# ---------- NRFA station list ----------

def fetch_stations():
    if os.path.exists(STATION_CACHE):
        age = time.time() - os.path.getmtime(STATION_CACHE)
        if age < _station_cache_ttl:
            with open(STATION_CACHE, "r", encoding="utf-8") as f:
                return json.load(f)

    req = urllib.request.Request(STATION_LIST_URL, headers={"User-Agent": "nrfa-site-tool-web"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))["data"]
    with open(STATION_CACHE, "w", encoding="utf-8") as f:
        json.dump(data, f)
    return data


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return r * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def nearest_station(lat, lon, stations):
    best, best_dist = None, None
    for s in stations:
        if s.get("latitude") is None or s.get("longitude") is None:
            continue
        d = haversine_km(lat, lon, s["latitude"], s["longitude"])
        if best_dist is None or d < best_dist:
            best, best_dist = s, d
    return best, best_dist


# ---------- WGS84 -> OSGB36 ----------

def wgs84_to_osgb36_gridref(lat_deg, lon_deg):
    a1, b1 = 6378137.000, 6356752.3142
    e2_1 = 1 - (b1 ** 2 / a1 ** 2)
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sinlat, coslat = math.sin(lat), math.cos(lat)
    nu1 = a1 / math.sqrt(1 - e2_1 * sinlat ** 2)
    x1 = nu1 * coslat * math.cos(lon)
    y1 = nu1 * coslat * math.sin(lon)
    z1 = (1 - e2_1) * nu1 * sinlat

    tx, ty, tz = -446.448, 125.157, -542.060
    s = -20.4894e-6
    rx = math.radians(-0.1502 / 3600)
    ry = math.radians(-0.2470 / 3600)
    rz = math.radians(-0.8421 / 3600)

    x2 = tx + (1 + s) * x1 + (-rz) * y1 + ry * z1
    y2 = ty + rz * x1 + (1 + s) * y1 + (-rx) * z1
    z2 = tz + (-ry) * x1 + rx * y1 + (1 + s) * z1

    a2, b2 = 6377563.396, 6356256.909
    e2_2 = 1 - (b2 ** 2 / a2 ** 2)
    p = math.sqrt(x2 ** 2 + y2 ** 2)
    philat = math.atan2(z2, p * (1 - e2_2))
    for _ in range(10):
        nu2 = a2 / math.sqrt(1 - e2_2 * math.sin(philat) ** 2)
        philat = math.atan2(z2 + e2_2 * nu2 * math.sin(philat), p)
    lon2 = math.atan2(y2, x2)

    F0 = 0.9996012717
    lat0, lon0 = math.radians(49), math.radians(-2)
    N0, E0 = -100000, 400000
    n = (a2 - b2) / (a2 + b2)

    sinlat2, coslat2 = math.sin(philat), math.cos(philat)
    nu = a2 * F0 / math.sqrt(1 - e2_2 * sinlat2 ** 2)
    rho = a2 * F0 * (1 - e2_2) / (1 - e2_2 * sinlat2 ** 2) ** 1.5
    eta2 = nu / rho - 1

    Ma = (1 + n + 5 / 4 * n ** 2 + 5 / 4 * n ** 3) * (philat - lat0)
    Mb = (3 * n + 3 * n ** 2 + 21 / 8 * n ** 3) * math.sin(philat - lat0) * math.cos(philat + lat0)
    Mc = (15 / 8 * n ** 2 + 15 / 8 * n ** 3) * math.sin(2 * (philat - lat0)) * math.cos(2 * (philat + lat0))
    Md = (35 / 24 * n ** 3) * math.sin(3 * (philat - lat0)) * math.cos(3 * (philat + lat0))
    M = b2 * F0 * (Ma - Mb + Mc - Md)

    tanlat = math.tan(philat)
    tan2lat, tan4lat = tanlat ** 2, tanlat ** 4

    I = M + N0
    II = nu / 2 * sinlat2 * coslat2
    III = nu / 24 * sinlat2 * coslat2 ** 3 * (5 - tan2lat + 9 * eta2)
    IIIA = nu / 720 * sinlat2 * coslat2 ** 5 * (61 - 58 * tan2lat + tan4lat)
    IV = nu * coslat2
    V = nu / 6 * coslat2 ** 3 * (nu / rho - tan2lat)
    VI = nu / 120 * coslat2 ** 5 * (5 - 18 * tan2lat + tan4lat + 14 * eta2 - 58 * tan2lat * eta2)

    dlon = lon2 - lon0
    N = I + II * dlon ** 2 + III * dlon ** 4 + IIIA * dlon ** 6
    E = E0 + IV * dlon + V * dlon ** 3 + VI * dlon ** 5
    return E, N


# ---------- flow data + percentiles ----------

def fetch_flow_series(station_id):
    cache_file = os.path.join(FLOW_CACHE_DIR, f"{station_id}.json")
    if os.path.exists(cache_file):
        age = time.time() - os.path.getmtime(cache_file)
        if age < _station_cache_ttl:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)

    url = (
        f"https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series"
        f"?format=json-object&data-type=gdf&station={station_id}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "nrfa-site-tool-web"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.loads(resp.read().decode("utf-8"))

    values = [v for v in payload.get("data-stream", []) if isinstance(v, (int, float))]
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(values, f)
    return values


def percentile(sorted_vals, p):
    k = (len(sorted_vals) - 1) * (p / 100)
    f, c = math.floor(k), math.ceil(k)
    if f == c:
        return sorted_vals[int(k)]
    return sorted_vals[int(f)] * (c - k) + sorted_vals[int(c)] * (k - f)


# ---------- routes ----------

@app.route("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.route("/api/lookup", methods=["POST"])
def api_lookup():
    body = request.get_json(force=True)
    parsed = parse_coords(body.get("coords", ""))
    if parsed is None:
        return jsonify({"error": "Could not parse those coordinates."}), 400
    lat, lon = parsed

    stations = fetch_stations()
    station, dist_km = nearest_station(lat, lon, stations)
    if station is None or station.get("catchment-area") is None:
        return jsonify({"error": "No usable nearest station found."}), 400

    easting, northing = wgs84_to_osgb36_gridref(lat, lon)

    return jsonify({
        "lat": lat, "lon": lon,
        "easting": round(easting), "northing": round(northing),
        "station": {
            "id": station["id"],
            "name": station["name"],
            "river": station.get("river"),
            "catchment_area": station.get("catchment-area"),
            "distance_km": round(dist_km, 1),
            "url": f"https://nrfa.ceh.ac.uk/data/station/info/{station['id']}",
        },
    })


@app.route("/api/compute", methods=["POST"])
def api_compute():
    body = request.get_json(force=True)
    station_id = body.get("station_id")
    station_area = body.get("station_area")
    site_area = body.get("site_area")
    percentiles = body.get("percentiles", [])

    if not all([station_id, station_area, site_area]):
        return jsonify({"error": "Missing station_id, station_area or site_area."}), 400

    try:
        station_area = float(station_area)
        site_area = float(site_area)
        percentiles = [float(p) for p in percentiles]
    except (ValueError, TypeError):
        return jsonify({"error": "Invalid numeric input."}), 400

    flows = fetch_flow_series(station_id)
    if not flows:
        return jsonify({"error": "No flow data available for this station."}), 400
    flows_sorted = sorted(flows)

    ratio = site_area / station_area
    mean_station = sum(flows) / len(flows)
    mean_site = mean_station * ratio

    rows = []
    for p in percentiles:
        station_val = percentile(flows_sorted, 100 - p)
        rows.append({"label": f"Q{p:g}", "station": round(station_val, 3), "site": round(station_val * ratio, 3)})

    return jsonify({
        "ratio_pct": round(ratio * 100, 1),
        "mean_station": round(mean_station, 3),
        "mean_site": round(mean_site, 3),
        "rows": rows,
    })


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)), debug=True)
