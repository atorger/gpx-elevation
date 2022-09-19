import math
import statistics
import re
import os
import numpy as np
import csv
from shapely import geometry
from itertools import islice
from sortedcontainers import SortedDict
from dataclasses import dataclass
import scipy.signal
import gpxpy
import gpxpy.gpx

@dataclass
class Point:
    x: float
    y: float
    def __hash__(self):
        return hash((self.x, self.y))
    def __getitem__(self, i):
        if i == 0:
            return self.x
        return self.y

@dataclass
class FixPoint:
    x: float
    y: float
    elevation: float
    fixtype: str
    name: str
    dist: float=None
    def __hash__(self):
        return hash((self.x, self.y, self.elevation, self.fixtype))

def find_closest(sorted_dict, key):
    if len(sorted_dict) == 0:
        return None
    keys = list(islice(sorted_dict.irange(minimum=key), 1))
    keys.extend(islice(sorted_dict.irange(maximum=key, reverse=True), 1))
    return min(keys, key=lambda k: abs(key - k))

def apply_lowpass_filter(data, critical_frequency):
    b, a = scipy.signal.butter(3, critical_frequency)
    return scipy.signal.filtfilt(b, a, data)

def sample_curve(curve, step):
    d = SortedDict()
    for point in curve:
        x = point[0]
        y = point[1]
        d[x] = y
    start_x = curve[0][0]
    end_x = curve[-1][0]
    x = start_x
    segments = []
    data = []
    while x <= end_x:
        if x in d:
            if d[x] is None:
                data = []
            else:
                if len(data) == 0:
                    segments.append((x, data))
                data.append(d[x])
        else:
            keys = list(islice(d.irange(minimum=x), 1))
            keys.extend(islice(d.irange(maximum=x, reverse=True), 1))
            y1 = d[keys[0]]
            y2 = d[keys[1]]
            if y1 is None or y2 is None:
                data = []
            else:
                if len(data) == 0:
                    segments.append((x, data))
                mix = (x - keys[0]) / (keys[1] - keys[0])
                y = y1 * (1 - mix) + y2 * mix
                data.append(y)
        x += step
    return segments

def get_clustered_average(values, spans):

    def get_average_within_span(sorted_values, span):
        # get average of the largest group which is within span from eachother
        values = sorted_values
        best_count = None
        for idx in range(len(values)):
            minval = values[idx]
            maxval = minval
            valsum = minval
            i = idx + 1
            while i < len(values) and values[i] <= minval + span:
                valsum += values[i]
                if maxval < values[i]:
                    maxval = values[i]
                i += 1
            count = i - idx
            if best_count is None or count > best_count:
                avgval = valsum / count
                inner_span = maxval - minval
                best_count = count
        return (avgval, best_count, inner_span)

    if len(values) == 0:
        return None, 0, 0
    if len(values) <= 2:
        return statistics.median(values), len(values), max(values) - min(values)
    values = sorted(values)
    for span in spans:
        val, count, inner_span = get_average_within_span(values, span)
        if count >= 3 or count == len(values) or (count == 2 and len(values) == 3 and span > 2):
            return val, count, inner_span
    return statistics.median(values), len(values), max(values) - min(values)

def calculate_elevation_gain(track, get_elevation_lambda):
    it = iter(track)
    prev_ele = get_elevation_lambda(next(it))
    elevation_gain = 0
    for point in it:
        ele = get_elevation_lambda(point)
        if ele > prev_ele:
            elevation_gain += ele - prev_ele
        prev_ele = ele
    return elevation_gain

def dist2d(p1, p2):
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    return math.sqrt(dx * dx + dy * dy)

def get_angle(p1, p2):
    x1 = p1.x
    y1 = p1.y
    x2 = p2.x
    y2 = p2.y
    angle = math.atan2(y2 - y1, x2 - x1) * (180 / math.pi)
    return (angle + 360) % 360

# Check if lines cross.
#
# Lines that lie exactly on top of each other do not cross, so if you need to check if a line
# touches the other you need an additional check.
def line_intersection(n1, n2, p1, p2):
    Nx1 = n1.x
    Ny1 = n1.y
    Nx2 = n2.x
    Ny2 = n2.y
    Px1 = p1.x
    Py1 = p1.y
    Px2 = p2.x
    Py2 = p2.y
    denom = ((Py2 - Py1) * (Nx2 - Nx1)) - ((Px2 - Px1) * (Ny2 - Ny1))
    if denom == 0:
        return None
    a = Ny1 - Py1
    b = Nx1 - Px1
    num1 = ((Px2 - Px1) * a) - ((Py2 - Py1) * b)
    num2 = ((Nx2 - Nx1) * a) - ((Ny2 - Ny1) * b)
    a = num1 / denom
    b = num2 / denom
    if 0 < a < 1 and 0 < b < 1:
        return Point(Nx1 + (a * (Nx2 - Nx1)), Ny1 + (a * (Ny2 - Ny1)))
    return None

def get_center_lat_lon(track):
    lat_min = min([row[0] for row in track])
    lat_max = max([row[0] for row in track])
    lon_min = min([row[1] for row in track])
    lon_max = max([row[1] for row in track])
    # Get Bounds of user selected area
    bound_points_shapely = geometry.MultiPoint([(lon_min, lat_min), (lon_max, lat_max)])
    # get lat/lng from center (below: True centroid - since coords may be multipoint)
    lat = bound_points_shapely.centroid.coords[0][1]
    lon = bound_points_shapely.centroid.coords[0][0]
    return (lat, lon)

def best_utm(lat: float, lon: float):
    """Based on lat and lng, return best utm epsg-code"""
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
        return epsg_code
    epsg_code = '327' + utm_band
    return f'epsg:{epsg_code}'

def write_track_csv(filename, data):
    with open(filename, 'w') as f:
        for item in data:
            f.write(f'{item[0]},{item[1]},{item[2]}\n')

def read_track(filename):
    _, ext = os.path.splitext(filename)

    if ext.lower() == '.gpx':
        data = open(filename).read()

        lat = np.array(re.findall(r'lat="([^"]+)',data),dtype=float)
        lon = np.array(re.findall(r'lon="([^"]+)',data),dtype=float)
        ele = np.array(re.findall(r'<ele>([^\<]+)',data),dtype=float)

        return list(zip(lat,lon,ele))
    if ext.lower() == '.csv':
        with open(filename) as f:
            reader = csv.reader(f)
            return list(map(lambda item: (float(item[0]), float(item[1]), float(item[2])), reader))
    print(f'Unknown/unsupported extension {ext}')
    assert False

def write_gpx(filename, track, transform):
    gpx = gpxpy.gpx.GPX()
    gpx_track = gpxpy.gpx.GPXTrack()
    gpx.tracks.append(gpx_track)
    gpx_segment = gpxpy.gpx.GPXTrackSegment()
    gpx_track.segments.append(gpx_segment)
    for point in track:
        lat, lon = transform(point[0], point[1])
        ele = round(point[2], 1)
        gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(lat, lon, ele))

    with open(filename, 'w') as f:
        f.write(gpx.to_xml())
