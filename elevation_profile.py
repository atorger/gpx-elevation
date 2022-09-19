import os
import math
from functools import cmp_to_key
from sortedcontainers import SortedDict
from dataclasses import dataclass
from itertools import islice
from pyproj import Transformer

from misc import find_closest, get_angle, apply_lowpass_filter, sample_curve, line_intersection, get_clustered_average, dist2d, best_utm, read_track, get_center_lat_lon, FixPoint, Point
from twodimsearch import TwoDimSearch

def join_sampled_segments(seg0, seg0_end_dist, seg1, seg1_start_dist, sampling_step):
    x0 = seg0_end_dist
    y0 = seg0[-1]
    x1 = seg1_start_dist
    y1 = seg1[0]
    x = seg0_end_dist + sampling_step
    while x < seg1_start_dist:
        mix = (x - x0) / (x1 - x0)
        y = y0 * (1-mix) + y1 * mix
        seg0.append(y)
        x += sampling_step
    return seg0 + seg1

def generate_correction_curve(fp_pos_list, track_distance, sampling_step):
    if len(fp_pos_list) == 0:
        return None
    length = math.ceil(track_distance/sampling_step)
    corr_list = []
    idx = 0
    dist = 0
    fp_idx = 0
    last_dist = 0
    last_corr = fp_pos_list[0][1]
    fp_pos_list.append(((length + 2) * sampling_step, fp_pos_list[-1][1]))
    while idx < length:
        dist = idx * sampling_step
        fp_dist = fp_pos_list[fp_idx][0] * sampling_step
        if fp_dist == last_dist:
            mix = 1
        else:
            mix = (fp_dist - dist) / (fp_dist - last_dist)
        if dist > fp_dist:
            last_corr = fp_pos_list[fp_idx][1]
            last_dist = fp_dist
            fp_idx += 1
        else:
            corr = last_corr * mix + fp_pos_list[fp_idx][1] * (1 - mix)
            corr_list.append(corr)
            idx += 1
    return corr_list

def get_minmax_span(profile, idx):
    idx_lo = idx
    idx_hi = idx
    while idx_lo > 0 and profile[idx_lo] == profile[idx]:
        idx_lo -= 1
    while idx_hi < len(profile) and profile[idx_hi] == profile[idx]:
        idx_hi += 1
    return (idx_lo, idx_hi)

def is_maxima(profile, idx):
    idx_lo, idx_hi = get_minmax_span(profile, idx)
    return idx_lo > 0 and profile[idx_lo] < profile[idx] and idx_hi < len(profile) and profile[idx_hi] < profile[idx]

def is_minima(profile, idx):
    idx_lo, idx_hi = get_minmax_span(profile, idx)
    return idx_lo > 0 and profile[idx_lo] > profile[idx] and idx_hi < len(profile) and profile[idx_hi] > profile[idx]

def get_nearest_perpendicular(track, track_db, fp, track_scan_dist):
    # get track points nearby the fixpoint
    refs = track_db.find_all_within((fp.x, fp.y), 2 * track_scan_dist)
    if len(refs) == 0:
        return None
    idxs = sorted(list(refs))

    # find the trackpoint that is nearest the fixpoint, with perpendicular search to track direction
    best_match = None
    del idxs[-1]
    side_scan_dist = 20.0
    for idx in idxs:
        tp1 = track[idx]
        tp2 = track[idx+1]
        angle = get_angle(tp1, tp2)
        angle = ((angle + 90) % 360) * math.pi / 180
        x1 = fp.x + math.cos(angle) * side_scan_dist
        y1 = fp.y + math.sin(angle) * side_scan_dist
        x2 = fp.x - math.cos(angle) * side_scan_dist
        y2 = fp.y - math.sin(angle) * side_scan_dist
        p1 = Point(x1, y1)
        p2 = Point(x2, y2)
        ip = line_intersection(tp1, tp2, p1, p2)
        if ip is None:
            continue
        dist = dist2d(ip, fp)
        if best_match is None or dist < best_match[1]:
            best_match = (idx, dist, ip)
    return best_match

def get_scan_range_for_fix(ref_track, profile_start_dist, profile, sampling_step, fp, track_scan_dist):
    side_scan_dist = 20.0
    point = ref_track.get_nearest_perpendicular_point(fp, side_scan_dist)
    if point is None:
        return (None, None, None)
    ref_dist = point.dist
    if ref_dist < profile_start_dist or ref_dist > profile_start_dist + len(profile) * sampling_step:
        return (None, None, None)

    start_dist = ref_dist - profile_start_dist - track_scan_dist
    end_dist = ref_dist - profile_start_dist + track_scan_dist
    start_idx = math.floor(start_dist / sampling_step)
    end_idx = math.floor(end_dist / sampling_step) + 1
    if start_idx < 0:
        start_idx = 0
    if end_idx > len(profile)-2:
        end_idx = len(profile)-2
    return (start_idx, end_idx, ref_dist - profile_start_dist)

def get_best_match_for_slope_fix(ref_track, start_dist, profile, sampling_step, fp, track_scan_dist):
    start_idx, end_idx, ref_dist = get_scan_range_for_fix(ref_track, start_dist, profile, sampling_step, fp, track_scan_dist)
    if start_idx is None:
        return None

    match_dist = None
    idx = start_idx
    while idx <= end_idx:
        if (profile[idx] < fp.elevation and profile[idx+1] > fp.elevation) or ( profile[idx] > fp.elevation and profile[idx+1] < fp.elevation):
            mix = abs(profile[idx] - fp.elevation) / abs(profile[idx] - profile[idx+1])
            match_dist = (idx + mix) * sampling_step
            break
        idx += 1
    if match_dist is None:
        return None
    corr_dist = match_dist - ref_dist
    return (idx, corr_dist)

def get_best_match_for_minmax_fix(ref_track, start_dist, profile, sampling_step, fp, track_scan_dist):
    start_idx, end_idx, _ = get_scan_range_for_fix(ref_track, start_dist, profile, sampling_step, fp, track_scan_dist)
    if start_idx is None:
        return None

    min_ele_idx = None
    max_ele_idx = None
    idx = start_idx
    while idx <= end_idx:
        if is_minima(profile, idx) and (min_ele_idx is None or profile[idx] < profile[min_ele_idx]):
            min_ele_idx = idx
        if is_maxima(profile, idx) and (max_ele_idx is None or profile[idx] > profile[max_ele_idx]):
            max_ele_idx = idx
        idx += 1

    if (fp.fixtype == 'max' and max_ele_idx is None) or (fp.fixtype == 'min' and min_ele_idx is None):
        # for robustness: if no real maxima or minima was found, then we just go for max/min value
        idx = start_idx
        min_ele_idx = idx
        max_ele_idx = idx
        while idx <= end_idx:
            if profile[idx] < profile[min_ele_idx]:
                min_ele_idx = idx
            if profile[idx] > profile[max_ele_idx]:
                max_ele_idx = idx
            idx += 1

    if fp.fixtype == 'max':
        return max_ele_idx
    elif fp.fixtype == 'min':
        return min_ele_idx
    else:
        raise RuntimeError(f'unknown/unexpected fixtype "{fp.fixtype}"')

@dataclass
class TrackPoint:
    x: float
    y: float
    elevation: float

@dataclass
class ProfileSegment:
    start_dist: float
    profile: list
    lowpass_profile: list=None
    vcorr_profile: list=None
    hvcorr_profile: list=None
    hvcorr_start_dist: float=None
    hvcorr_end_dist: float=None

class ElevationProfile:

    _plot_files = dict()

    def __init__(self, filename):

        # TODO: if we really want to be able to change reference track later, we need to do init with generic
        # projection
        print(f'Loading {filename}')
        track = read_track(filename)

        lat, lon = get_center_lat_lon(track)
        self._epsg_name = best_utm(lat, lon)
        wgs84_to_utm = Transformer.from_crs('epsg:4326', self._epsg_name)

        def make_track_point(item):
            x, y = wgs84_to_utm.transform(item[0], item[1])
            return TrackPoint(x, y, item[2])

        track = list(map(make_track_point, track))

        track_db = TwoDimSearch()
        for point in track:
            track_db.insert((point.x, point.y), point.elevation)

        print('Removing clustered points (likely standing still, bike laying down etc)')
        clean_track = []
        previous_had_multipoints = False
        for point in track:
            points = track_db.find_all_within((point.x, point.y), 2.0)
            if len(points) > 1:
                if not previous_had_multipoints:
                    clean_track.append(point)
                previous_had_multipoints = True
            else:
                previous_had_multipoints = False
                clean_track.append(point)

        print(f'Removed {len(track) - len(clean_track)} points ({len(clean_track)} points remain)')
        track = clean_track

        print('Fill in extra points where points are widely spaced')
        max_seg_len = 15
        new_track = []
        for idx, point in enumerate(track):
            if idx == 0:
                new_track.append(point)
            else:
                dist = dist2d(new_track[-1], track[idx])
                if dist > max_seg_len:
                    count = round(dist / max_seg_len + 1)
                    ele1 = new_track[-1].elevation
                    ele2 = track[idx].elevation
                    x1 = new_track[-1].x
                    x2 = track[idx].x
                    y1 = new_track[-1].y
                    y2 = track[idx].y
                    for i in range(1, count):
                        mix = i / count
                        x = x1 + (x2 - x1) * mix
                        y = y1 + (y2 - y1) * mix
                        ele = ele1 * (1 - mix) + ele2 * mix
                        new_track.append(TrackPoint(x, y, elevation=ele))
                new_track.append(point)
        track = new_track

        self._track = track
        self._track_db = TwoDimSearch()
        for idx, point in enumerate(track):
            self._track_db.insert((point.x, point.y), idx)

    def __len__(self):
        return len(self._profile)

    def _find_local_minmax_points(self, ps):
        minmax_pos = []
        for idx in range(1, len(ps.lowpass_profile)-2):
            if is_maxima(ps.lowpass_profile, idx):
                minmax_pos.append((idx, 'max'))
            elif is_minima(ps.lowpass_profile, idx):
                minmax_pos.append((idx, 'min'))
        prominent_minmax_pos = []
        window_len = math.ceil(200.0 / self._sampling_step)
        last_pos = -window_len / 2
        start_idx = round(ps.hvcorr_start_dist / self._sampling_step)
        for idx, (pos, minmax) in enumerate(minmax_pos):
            if idx < len(minmax_pos) - 1:
                next_pos = minmax_pos[idx+1][0]
            else:
                next_pos = len(minmax_pos) - 1 + window_len / 2
            if pos - window_len > last_pos and pos + window_len < next_pos:
                ele = ps.lowpass_profile[pos]
                prominent_minmax_pos.append((start_idx + pos, ele, minmax))
            last_pos = pos
        return prominent_minmax_pos

    @staticmethod
    def _match_minmax_fixpoints(ps, fp_db, ref_track, sampling_step, lowpass=False):
        fp_minmax_pos_list = []
        if lowpass:
            pr = ps.lowpass_profile
        else:
            pr = ps.profile
        for fpset in fp_db:
            fp = next(iter(fpset))
            if fp.fixtype != 'max' and fp.fixtype != 'min':
                continue
            idx = get_best_match_for_minmax_fix(ref_track, ps.start_dist, pr, sampling_step, fp, 100.0)
            if idx is None:
                continue
            corr = fp.elevation - pr[idx]
            idx += ps.start_dist / sampling_step
            fp_minmax_pos_list.append((idx, corr, fp))
        return fp_minmax_pos_list

    @staticmethod
    def _match_slope_fixpoints(ps, fp_db, ref_track, sampling_step):
        fp_horiz_pos_list = []
        start_idx = ps.start_dist / sampling_step
        for fpset in fp_db:
            fp = next(iter(fpset))
            if fp.fixtype != 'slope':
                continue
            match = get_best_match_for_slope_fix(ref_track, ps.start_dist, ps.vcorr_profile, sampling_step, fp, 100.0)
            if match is None:
                continue
            fp_horiz_pos_list.append((match[0] + start_idx, match[1], fp))
        return fp_horiz_pos_list

    def _apply_hcorr(self):
        for ps in self._profile_segments:
            start_idx = int(ps.start_dist / self._sampling_step)
            if self._hcorr is None:
                ps.hvcorr_start_dist = ps.start_dist
                ps.hvcorr_profile = ps.vcorr_profile
            else:
                ps.hvcorr_start_dist = ps.start_dist - self._hcorr[start_idx]
                ps.hvcorr_profile = []
                idx = start_idx
                while True:
                    if idx < len(self._hcorr)-1:
                        corr = self._hcorr[idx]
                    else:
                        corr = self._hcorr[-1]
                    dist = (idx-start_idx) * self._sampling_step + corr + (ps.hvcorr_start_dist - ps.start_dist)
                    if dist >= 0:
                        idx0 = math.floor(dist / self._sampling_step)
                        if idx0 > len(ps.vcorr_profile) - 1:
                            break
                        mix = (dist % self._sampling_step) / self._sampling_step
                        if mix == 0:
                            ele = ps.vcorr_profile[idx0]
                        else:
                            if idx0 > len(ps.vcorr_profile) - 2:
                                break
                            ele = ps.vcorr_profile[idx0] * (1 - mix) + ps.vcorr_profile[idx0+1] * mix
                        ps.hvcorr_profile.append(ele)
                    idx += 1
            ps.hvcorr_end_dist = ps.hvcorr_start_dist + len(ps.hvcorr_profile) * self._sampling_step

    def _calculate_lowpass_profiles(self):
        for ps in self._profile_segments:
            ps.lowpass_profile = apply_lowpass_filter(ps.hvcorr_profile, self._sampling_step / 300.0)

            # One can make useful analysis with profile gradient, for example horizontal alignment, however
            # got good enough results without it so not used for now
            #x = np.array(ps.lowpass_profile, dtype=float)
            #ps.gradient = list(np.gradient(x))

    def _apply_vcorr(self):
        for ps in self._profile_segments:
            if self._vcorr is None:
                ps.vcorr_profile = ps.profile
                continue
            ps.vcorr_profile = []
            start_idx = round(ps.start_dist / self._sampling_step)
            for idx, ele in enumerate(ps.profile):
                ps.vcorr_profile.append(ele + self._vcorr[start_idx + idx])

    def set_reference_track(self, ref_track, fp_db, side_scan_dist, sampling_step):

        print(f'Making profile for reference track')

        self._sampling_step = sampling_step
        self._ref_track = ref_track
        self._side_scan_dist = side_scan_dist
        self._vcorr = None
        self._hcorr = None

        if ref_track.utm_epsg_name() != self._epsg_name:
            # rare case, need to reproject track and track_db
            reproject = Transformer.from_crs(self._epsg_name, ref_track.utm_epsg_name())

            def make_track_point(item):
                x, y = reproject.transform(item.x, item.y)
                return TrackPoint(x, y, item.elevation)

            self._track = list(map(make_track_point, self._track))
            self._track_db = TwoDimSearch()
            for idx, point in enumerate(self._track):
                self._track_db.insert((point.x, point.y), idx)
            self._epsg_name = ref_track.utm_epsg_name()

        profile = list(map(lambda tp: (tp.dist, self._get_elevation(tp.dist)), ref_track))
        for idx, p in enumerate(profile):
            if idx > 0 and idx < len(profile) and p[1] != None and profile[idx-1][1] is None and profile[idx+1][1] is None:
                print(f'At least two consecutive points with elevation required, ignoring point at {profile[idx][0]}')
                profile[idx] = (profile[idx][0], None)

        profile_segs = sample_curve(profile, self._sampling_step)
        seg_info = list(map(lambda seg: (seg[0], seg[0] + (len(seg[1])-1) * self._sampling_step), profile_segs))
        print(f'The track contains following segment range(s): {seg_info}')

        # bridge small gaps between segments
        bridge_gap = 50.0
        prev_end_dist = None
        joined_profile_segs = []
        for idx, seg in enumerate(profile_segs):
            start_dist = seg[0]
            if prev_end_dist is not None and start_dist - prev_end_dist <= bridge_gap:
                prev = joined_profile_segs.pop()
                joined = join_sampled_segments(prev[1], prev_end_dist, seg[1], start_dist, self._sampling_step)
                joined_profile_segs.append((prev[0], joined))
            else:
                joined_profile_segs.append(seg)
            prev_end_dist = joined_profile_segs[-1][0] + (len(joined_profile_segs[-1][1]) - 1) * self._sampling_step
        if len(joined_profile_segs) < len(profile_segs):
            seg_info = list(map(lambda seg: (seg[0], seg[0] + (len(seg[1])-1) * self._sampling_step), joined_profile_segs))
            print(f'The track contains following segment range(s) after joining nearby segments: {seg_info}')
        profile_segs = joined_profile_segs

        self._profile_segments = []
        fp_minmax_pos_list = []
        for profile_segment in profile_segs:
            ps = ProfileSegment(profile_segment[0], profile_segment[1])

            edge_remove_count = math.floor(200.0 / self._sampling_step)
            if ps.start_dist > 0:
                ps.start_dist += edge_remove_count * self._sampling_step
                ps.profile = ps.profile[edge_remove_count:]
            end_dist = ps.start_dist + len(ps.profile) * self._sampling_step
            if end_dist < profile[-1][0]:
                ps.profile = ps.profile[:-edge_remove_count]

            if len(ps.profile) < 2:
                continue

            ps.profile = apply_lowpass_filter(ps.profile, self._sampling_step / 20.0)
            self._profile_segments.append(ps)

            fp_minmax_pos_list += ElevationProfile._match_minmax_fixpoints(ps, fp_db, self._ref_track, self._sampling_step)

        fp_minmax_pos_list = sorted(fp_minmax_pos_list, key=lambda x: x[0])
        for match in fp_minmax_pos_list:
            print(f'dist={round(ps.start_dist + match[0] * self._sampling_step, 1)} elevation difference={round(match[1], 1)} name={match[2].name}')

        self._vcorr = generate_correction_curve(fp_minmax_pos_list, self._ref_track.distance(), self._sampling_step)
        self._apply_vcorr()

        fp_horiz_pos_list = []
        for ps in self._profile_segments:
            fp_horiz_pos_list += ElevationProfile._match_slope_fixpoints(ps, fp_db, self._ref_track, self._sampling_step)

        fp_horiz_pos_list = sorted(fp_horiz_pos_list, key=lambda x: x[0])
        for match in fp_horiz_pos_list:
            fp = match[2]
            print(f'dist={round(ps.start_dist + match[0] * self._sampling_step,1)} elevation={round(fp.elevation,1)} horizontal distance to fix={round(match[1],1)} name={fp.name}')

        self._hcorr = generate_correction_curve(fp_horiz_pos_list, self._ref_track.distance(), self._sampling_step)
        self._apply_hcorr()
        self._calculate_lowpass_profiles()

    def _get_elevation(self, ref_track_dist):

        tp = self._ref_track.get_point_at_dist(ref_track_dist)
        if tp is None:
            return None
        refs = self._track_db.find_all_within((tp.x, tp.y), 2 * self._side_scan_dist)
        if len(refs) == 0:
            return None

        angle = ((tp.angle + 90) % 360) * math.pi / 180
        x1 = tp.x + math.cos(angle) * self._side_scan_dist
        y1 = tp.y + math.sin(angle) * self._side_scan_dist
        x2 = tp.x - math.cos(angle) * self._side_scan_dist
        y2 = tp.y - math.sin(angle) * self._side_scan_dist
        p1 = Point(x1, y1)
        p2 = Point(x2, y2)

        closest_match = None
        for idx in refs:
            if idx < len(self._track)-1:
                idx1 = idx
                idx2 = idx+1
            else:
                idx1 = idx-1
                idx2 = idx
            tp1 = self._track[idx1]

            if dist2d(tp1, tp) < 0.01: # very close, don't attempt intersection
                return tp1.elevation

            tp2 = self._track[idx2]
            ip = line_intersection(tp1, tp2, p1, p2)
            if ip is not None:
                mix = dist2d(ip, tp1) / dist2d(tp2, tp1)
                ele = tp1.elevation * (1 - mix) + tp2.elevation * mix
                dist = dist2d(ip, tp)
                if closest_match is None or closest_match[0] > dist:
                    closest_match = (dist, ele)
        if closest_match is None:
            return None

        elevation = closest_match[1]

        return elevation

    def _get_datapoint(self, get_sampled_profile, ref_track_dist):
        d = ref_track_dist
        match = False
        for ps in self._profile_segments:
            if d >= ps.hvcorr_start_dist and d < ps.hvcorr_end_dist:
                match = True
                break
        if not match:
            return None

        sampled_profile = get_sampled_profile(ps)
        d -= ps.hvcorr_start_dist
        idx = math.floor(d / self._sampling_step)
        mix = (d - idx * self._sampling_step) / self._sampling_step
        val1 = sampled_profile[idx]
        if mix > 0 and idx < len(sampled_profile) - 1:
            val2 = sampled_profile[idx+1]
            return val1 * (1 - mix) + val2 * mix
        return val1

    def get_elevation(self, ref_track_dist):
        return self._get_datapoint(lambda ps: ps.hvcorr_profile, ref_track_dist)

    def seg_count(self):
        return len(self._profile_segments)

    @staticmethod
    def close_plot_files():
        for f in ElevationProfile._plot_files.values():
            f.close()

    @staticmethod
    def _open_plot_files(plot_dir):
        if len(ElevationProfile._plot_files) > 0:
            return
        names = ['op', 'vc', 'hc', 'vp', 'cp', 'fp']
        for name in names:
            ElevationProfile._plot_files[name] = open(os.path.join(plot_dir, name), 'w')

    def write_plot_data(self, plot_dir):
        ElevationProfile._open_plot_files(plot_dir)
        if self._vcorr is not None:
            f = ElevationProfile._plot_files['vc']
            for idx, val in enumerate(self._vcorr):
                f.write(f'{idx*self._sampling_step}, {val}\n')
            f.write('\n\n')
            f.flush()
        if self._hcorr is not None:
            f = ElevationProfile._plot_files['hc']
            for idx, val in enumerate(self._hcorr):
                f.write(f'{idx*self._sampling_step}, {val}\n')
            f.write('\n\n')
            f.flush()
        for ps in self._profile_segments:
            f = ElevationProfile._plot_files['op']
            for idx, ele in enumerate(ps.profile):
                f.write(f'{ps.start_dist+idx*self._sampling_step}, {ele}\n')
            f.write('\n\n')
            f.flush()
            f = ElevationProfile._plot_files['vp']
            for idx, ele in enumerate(ps.vcorr_profile):
                f.write(f'{ps.start_dist+idx*self._sampling_step}, {ele}\n')
            f.write('\n\n')
            f.flush()
            f = ElevationProfile._plot_files['cp']
            for idx, ele in enumerate(ps.hvcorr_profile):
                f.write(f'{ps.hvcorr_start_dist+idx*self._sampling_step}, {ele}\n')
            f.write('\n\n')
            f.flush()
            f = ElevationProfile._plot_files['fp']
            for idx, ele in enumerate(ps.lowpass_profile):
                f.write(f'{ps.hvcorr_start_dist+idx*self._sampling_step}, {ele}\n')
            f.write('\n\n')
            f.flush()

    @staticmethod
    def _index_local_minmax_points(profiles):
        mmpos = SortedDict()
        for profile in profiles:
            for ps in profile._profile_segments:
                prominent_minmax_pos = ElevationProfile._find_local_minmax_points(profile, ps)
                for (idx, ele, minmax) in prominent_minmax_pos:
                    if idx not in mmpos:
                        mmpos[idx] = [(minmax, ele)]
                    else:
                        mmpos[idx].append((minmax, ele))
        return mmpos

    @staticmethod
    def _make_point_clusters(mmpos, window_len):
        clustered_mmpos = SortedDict()
        while True:
            remaining_mmpos = SortedDict()
            current_minmax = None
            for idx, items in mmpos.items():
                for (minmax, ele) in items:
                    if current_minmax is None or idx >= current_minmax[0] + window_len:
                        current_minmax = (idx, minmax)
                        clustered_mmpos[idx] = [(idx, minmax, ele)]
                    elif minmax != current_minmax[1]:
                        if idx not in remaining_mmpos:
                            remaining_mmpos[idx] = [(minmax, ele)]
                        else:
                            remaining_mmpos[idx].append((minmax, ele))
                    else:
                        clustered_mmpos[current_minmax[0]].append((idx, minmax, ele))
            if len(remaining_mmpos) == 0:
                break
            mmpos = remaining_mmpos
        return clustered_mmpos

    @staticmethod
    def _sort_point_clusters(clustered_mmpos, min_dist_between_fixpoints, small_span, processed_set):
        minmax_values = []
        for idx, items in clustered_mmpos.items():
            closest_idx = find_closest(processed_set, idx)
            if closest_idx is not None and abs(closest_idx - idx) < min_dist_between_fixpoints:
                continue
            minmax = items[0][1]
            eles = [i[2] for i in items]
            val, count, span = get_clustered_average(eles, spans=[small_span])
            if count == 1:
                continue
            idxs = [i[0] for i in items]
            mid_idx, _, _ = get_clustered_average(idxs, spans=[30, 60, 120])
            minmax_values.append((mid_idx, minmax, val, count, span, idx))

        def cmp_minmax(a, b):
            if a[4] <= small_span and b[4] <= small_span:
                if a[1] != b[1]:
                    # prefer max over min
                    if a[1] == 'max':
                        return -1
                    return 1
                if a[3] != b[3]:
                    return b[3] - a[3] # prefer higher count
                if a[2] != b[2]:
                    return b[2] - a[2] # prefer higher elevation
                return a[4] - b[4] # prefer lower span
            return a[4] - b[4] # prefer lower span

        minmax_values.sort(key=cmp_to_key(cmp_minmax))

        return minmax_values

    @staticmethod
    def estimate_fixpoints_and_apply_corrections(profiles, provided_fp_db, plot_dir=None):
        ref_track = profiles[0]._ref_track
        sampling_step = profiles[0]._sampling_step

        vspan = 12.0
        min_dist_between_fixpoints = math.ceil(1500.0 / sampling_step)
        window_length = math.ceil(200.0 / sampling_step)
        processed_set = SortedDict()

        mmpos = ElevationProfile._index_local_minmax_points(profiles)
        clustered_mmpos = ElevationProfile._make_point_clusters(mmpos, window_length)

        sorted_values = ElevationProfile._sort_point_clusters(clustered_mmpos, min_dist_between_fixpoints, vspan, processed_set)
        print(sorted_values)

        provided_fp_positions = SortedDict()
        for fpset in provided_fp_db:
            fp = next(iter(fpset))
            point = ref_track.get_nearest_perpendicular_point(fp, profiles[0]._side_scan_dist)
            if point is not None:
                idx = point.dist / sampling_step
                provided_fp_positions[idx] = True

        spaced_out_values = []
        picked_points = SortedDict()
        for val in sorted_values:
            mid_idx = val[0]

            # use extra distance from provided fixpoints
            closest_idx = find_closest(provided_fp_positions, mid_idx)
            if closest_idx is not None and abs(closest_idx - mid_idx) < 2 * min_dist_between_fixpoints:
                continue

            closest_idx = find_closest(picked_points, mid_idx)
            if closest_idx is not None and abs(closest_idx - mid_idx) < min_dist_between_fixpoints:
                continue
            spaced_out_values.append(val)
            picked_points[mid_idx] = True
        sorted_values = spaced_out_values

        if plot_dir is not None:
            with open(os.path.join(plot_dir, 'derived_minmax'), 'w') as f:
                for val in sorted_values:
                    f.write(f'{val[0]*sampling_step} {val[2]} {val[4]/2}\n')

        if len(sorted_values) == 0:
            print('No derived fixpoints')
            return False

        print(f'Derived {len(sorted_values)} fix points, applying correction')
        fp_db = TwoDimSearch()
        for val in sorted_values:
            dist = val[0] * sampling_step
            fixtype = val[1]
            ele = val[2]
            tp = ref_track.get_point_at_dist(dist)
            fixpoint = FixPoint(tp.x, tp.y, ele, fixtype, 'derived')
            fixpoint.dist = dist
            fp_db.insert(Point(tp.x, tp.y), fixpoint)

        for profile in profiles:
            fp_minmax_pos_list = []
            fp_horiz_pos_list = []
            for ps in profile._profile_segments:
                fp_minmax_pos_list += ElevationProfile._match_minmax_fixpoints(ps, provided_fp_db, ref_track, profile._sampling_step)
                fp_horiz_pos_list += ElevationProfile._match_slope_fixpoints(ps, provided_fp_db, ref_track, profile._sampling_step)

                mmlist = ElevationProfile._match_minmax_fixpoints(ps, fp_db, ref_track, profile._sampling_step, lowpass=True)
                fp_minmax_pos_list += mmlist

                start_idx = ps.start_dist / profile._sampling_step
                for val in mmlist:
                    fp = val[2]
                    dst_idx = fp.dist / profile._sampling_step
                    src_idx = val[0]
                    corr_dist = src_idx * profile._sampling_step - fp.dist
                    fp_horiz_pos_list.append((src_idx, corr_dist, fp))

            fp_minmax_pos_list = sorted(fp_minmax_pos_list, key=lambda x: x[0])
            fp_horiz_pos_list = sorted(fp_horiz_pos_list, key=lambda x: x[0])
            #print(list(map(lambda i: (i[0], i[1]), fp_minmax_pos_list)))
            #print(list(map(lambda i: (i[0], i[1]), fp_horiz_pos_list)))
            profile._vcorr = generate_correction_curve(fp_minmax_pos_list, ref_track.distance(), profile._sampling_step)
            profile._hcorr = generate_correction_curve(fp_horiz_pos_list, ref_track.distance(), profile._sampling_step)
            ElevationProfile._apply_vcorr(profile)
            ElevationProfile._apply_hcorr(profile)
            ElevationProfile._calculate_lowpass_profiles(profile)
        return True
