import collections
import bisect
import math
from pyproj import Transformer
from dataclasses import dataclass

from twodimsearch import TwoDimSearch
from misc import line_intersection, best_utm, read_track, get_center_lat_lon, dist2d, get_angle, Point

@dataclass
class TrackPoint:
    x: float
    y: float
    elevation: float = 0.0
    angle: float = 0.0
    dist: float = 0.0

class ReferenceTrack:

    def __init__(self, filename):

        print(f'Loading reference track {filename}')
        track = read_track(filename)

        lat, lon = get_center_lat_lon(track)
        self._epsg_name = best_utm(lat, lon)

        print(f'Best UTM for track center (lat={lat}, lon={lon}) is {self._epsg_name}')
        # epsg:4326 = WGS84 (GPS standard coordinate system)
        self._wgs84_to_utm = Transformer.from_crs('epsg:4326', self._epsg_name)

        def make_track_point(item):
            x, y = self._wgs84_to_utm.transform(item[0], item[1])
            return TrackPoint(x, y, elevation=item[2])

        track = list(map(make_track_point, track))

        print('Log reference profile')
        ReferenceTrack._calculate_track_angles_and_dist(track)
        with open('profile-ref.txt', 'w') as f:
            for point in track:
                f.write(f'{point.dist} {point.elevation}\n')

        print('Merging coordinates unreasonably close')
        new_track = []
        for idx, point in enumerate(track):
            if idx == 0:
                new_track.append(point)
            else:
                if dist2d(new_track[-1], track[idx]) < 2.5:
                    x1 = new_track[-1].x
                    x2 = track[idx].x
                    y1 = new_track[-1].y
                    y2 = track[idx].y
                    x = x1 + (x2 - x1) / 2
                    y = y1 + (y2 - y1) / 2
                    new_track.pop()
                    new_track.append(TrackPoint(x, y))
                else:
                    new_track.append(point)
        track = new_track

        print('Fill in gaps to make sure elevation is probed often enough')
        max_seg_len = 20
        new_track = []
        for idx, point in enumerate(track):
            if idx == 0:
                new_track.append(point)
            else:
                dist = dist2d(new_track[-1], track[idx])
                if dist > max_seg_len:
                    count = round(dist / max_seg_len + 1)
                    x1 = new_track[-1].x
                    x2 = track[idx].x
                    y1 = new_track[-1].y
                    y2 = track[idx].y
                    for i in range(1, count):
                        x = x1 + (x2 - x1) * i / count
                        y = y1 + (y2 - y1) * i / count
                        new_track.append(TrackPoint(x, y))
                new_track.append(point)
        track = new_track

        print('Calculating angles and distance in each point')
        ReferenceTrack._calculate_track_angles_and_dist(track)

        self._track = track
        self._track_db = TwoDimSearch()
        for idx, point in enumerate(track):
            self._track_db.insert((point.x, point.y), idx)

    @staticmethod
    def _calculate_track_angles_and_dist(track):
        dist = 0
        for idx, point in enumerate(track):
            if idx == 0:
                angle = get_angle(point, track[idx+1])
            elif idx == len(track) - 1:
                angle = get_angle(track[idx-1], point)
            else:
                a = get_angle(track[idx-1], point)
                b = get_angle(point, track[idx+1])
                diff = ((a-b + 180+360) % 360) - 180
                angle = (360 + b + (diff / 2)) % 360
            point.angle = angle
            point.dist = dist
            if idx != len(track) - 1:
                dist += dist2d(point, track[idx+1])

    def __iter__(self):
        self._iter = iter(self._track)
        return self

    def __next__(self):
        return next(self._iter)

    def __len__(self):
        return len(self._track)

    def distance(self):
        return self._track[-1].dist

    def wgs84_to_utm(self):
        return self._wgs84_to_utm

    def utm_epsg_name(self):
        return self._epsg_name

    def get_point_at_dist(self, dist):
        if dist < 0 or dist > self._track[-1].dist:
            return None
        idx = bisect.bisect_left(self._track, dist, key=lambda x: x.dist)
        if self._track[idx].dist == dist:
            return self._track[idx]
        while self._track[idx].dist > dist:
            idx -= 1
        dist0 = dist - self._track[idx].dist
        dist1 = self._track[idx+1].dist - self._track[idx].dist
        mix = dist0 / dist1
        x = self._track[idx].x * (1 - mix) + self._track[idx+1].x * mix
        y = self._track[idx].y * (1 - mix) + self._track[idx+1].y * mix
        return TrackPoint(x, y, angle=self._track[idx].angle, dist=dist)

    def get_nearest_perpendicular_point(self, point, scan_dist):
        refs = self._track_db.find_all_within((point.x, point.y), 2 * scan_dist)
        if len(refs) == 0:
            return None
        idxs = sorted(list(refs))

        best_match = None
        for idx in idxs:
            if idx == len(self._track)-1:
                continue
            tp1 = self._track[idx]
            tp2 = self._track[idx+1]
            angle = get_angle(tp1, tp2)
            angle = ((angle + 90) % 360) * math.pi / 180
            x1 = point.x + math.cos(angle) * scan_dist
            y1 = point.y + math.sin(angle) * scan_dist
            x2 = point.x - math.cos(angle) * scan_dist
            y2 = point.y - math.sin(angle) * scan_dist
            p1 = Point(x1, y1)
            p2 = Point(x2, y2)
            ip = line_intersection(tp1, tp2, p1, p2)
            if ip is None:
                continue
            ipdist = dist2d(ip, point)
            mix = dist2d(ip, tp1) / dist2d(tp2, tp1)
            if best_match is None or ipdist < best_match[2]:
                track_dist = self._track[idx].dist + mix * (self._track[idx+1].dist - self._track[idx].dist)
                best_match = (ip, track_dist, ipdist)
        if best_match is None:
            return None
        return TrackPoint(best_match[0].x, best_match[0].y, dist=best_match[1])
