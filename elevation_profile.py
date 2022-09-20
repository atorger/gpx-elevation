import os
import math
from functools import cmp_to_key
from dataclasses import dataclass
from pyproj import Transformer
from sortedcontainers import SortedDict

from misc import find_closest, apply_lowpass_filter, sample_curve, line_intersection, get_clustered_average, dist2d, best_utm, read_track, get_center_lat_lon, FixPoint, Point
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
    start_idx = max(start_idx, 0)
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
    if fp.fixtype == 'min':
        return min_ele_idx
    raise RuntimeError(f'unknown/unexpected fixtype "{fp.fixtype}"')

@dataclass
class TrackPoint:
    x: float
    y: float
    elevation: float

@dataclass
class ProfileSegment:
    profile: list
    start_dist: float
    end_dist: float
    lowpass_profile: list=None

class ElevationProfile:

    _plot_files = {}

    def __init__(self, filename):

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
        start_idx = round(ps.start_dist / self._sampling_step)
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
    def _match_minmax_fixpoints(ps, fp_db, ref_track, window_len, sampling_step, lowpass=False):
        fp_minmax_pos_list = []
        if lowpass:
            pr = ps.lowpass_profile
        else:
            pr = ps.profile
        for fpset in fp_db:
            fp = next(iter(fpset))
            if fp.fixtype not in ('max', 'min'):
                continue
            idx = get_best_match_for_minmax_fix(ref_track, ps.start_dist, pr, sampling_step, fp, window_len/2)
            if idx is None:
                continue
            corr = fp.elevation - pr[idx]
            idx += ps.start_dist / sampling_step
            fp_minmax_pos_list.append((idx, corr, fp, ps))
        return fp_minmax_pos_list

    @staticmethod
    def _match_slope_fixpoints(ps, fp_db, ref_track, window_len, sampling_step):
        fp_horiz_pos_list = []
        start_idx = ps.start_dist / sampling_step
        for fpset in fp_db:
            fp = next(iter(fpset))
            if fp.fixtype != 'slope':
                continue
            match = get_best_match_for_slope_fix(ref_track, ps.start_dist, ps.profile, sampling_step, fp, window_len/2)
            if match is None:
                continue
            fp_horiz_pos_list.append((match[0] + start_idx, match[1], fp, ps))
        return fp_horiz_pos_list

    def _apply_hcorr(self, hcorr):
        if hcorr is None:
            return
        for ps in self._profile_segments:
            start_idx = int(ps.start_dist / self._sampling_step)
            hvcorr_start_dist = ps.start_dist - hcorr[start_idx]
            hvcorr_profile = []
            idx = start_idx
            while True:
                if idx < len(hcorr)-1:
                    corr = hcorr[idx]
                else:
                    corr = hcorr[-1]
                dist = (idx-start_idx) * self._sampling_step + corr + (hvcorr_start_dist - ps.start_dist)
                if dist >= 0:
                    idx0 = math.floor(dist / self._sampling_step)
                    if idx0 > len(ps.profile) - 1:
                        break
                    mix = (dist % self._sampling_step) / self._sampling_step
                    if mix == 0:
                        ele = ps.profile[idx0]
                    else:
                        if idx0 > len(ps.profile) - 2:
                            break
                        ele = ps.profile[idx0] * (1 - mix) + ps.profile[idx0+1] * mix
                    hvcorr_profile.append(ele)
                idx += 1
            ps.profile = hvcorr_profile
            ps.start_dist = hvcorr_start_dist
            ps.end_dist = ps.start_dist + len(ps.profile) * self._sampling_step

    def _apply_vcorr(self, vcorr):
        if vcorr is None:
            return
        for ps in self._profile_segments:
            vcorr_profile = []
            start_idx = round(ps.start_dist / self._sampling_step)
            for idx, ele in enumerate(ps.profile):
                vidx = min(start_idx + idx, len(vcorr) - 1)
                vcorr_profile.append(ele + vcorr[vidx])
            ps.profile = vcorr_profile

    def _calculate_lowpass_profiles(self):
        for ps in self._profile_segments:
            ps.lowpass_profile = apply_lowpass_filter(ps.profile, self._sampling_step / 300.0)

            # One can make useful analysis with profile gradient, for example horizontal alignment, however
            # got good enough results without it so not used for now
            #x = np.array(ps.lowpass_profile, dtype=float)
            #ps.gradient = list(np.gradient(x))

    def set_reference_track(self, ref_track, side_scan_dist, edge_remove_len, sampling_step):

        print('Making profile for reference track')

        self._sampling_step = sampling_step
        self._ref_track = ref_track
        self._side_scan_dist = side_scan_dist

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
            if 0 < idx < len(profile) and p[1] is not None and profile[idx-1][1] is None and profile[idx+1][1] is None:
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
        for profile_segment in profile_segs:
            ps = ProfileSegment(profile_segment[1], profile_segment[0], -1)

            # Remove edges of segments that don't start or end at the track's start or end.
            # The reason for this is that there is quite often elevation deviations where
            # the track leaves the reference track (=goes onto a different road) and thus
            # causes and end/start of segment.
            if edge_remove_len > 0:
                edge_remove_count = math.floor(edge_remove_len / self._sampling_step)
                if ps.start_dist > 0:
                    ps.start_dist += edge_remove_count * self._sampling_step
                    ps.profile = ps.profile[edge_remove_count:]
                end_dist = ps.start_dist + len(ps.profile) * self._sampling_step
                if end_dist < profile[-1][0]:
                    ps.profile = ps.profile[:-edge_remove_count]

            if len(ps.profile) < 2:
                continue

            # We apply some minimal low pass filtering of the input, as more details than
            # this is most likely just noise that will disturb the matching process.
            ps.profile = apply_lowpass_filter(ps.profile, self._sampling_step / 20.0)
            ps.end_dist = ps.start_dist + len(ps.profile) * self._sampling_step

            self._profile_segments.append(ps)

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

        closest_match = (None, None)
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
                if closest_match[0] is None or closest_match[0] > dist:
                    closest_match = (dist, ele)
        if closest_match is None:
            return None

        elevation = closest_match[1]

        return elevation

    def _get_datapoint(self, get_sampled_profile, ref_track_dist):
        d = ref_track_dist
        match = None
        for ps in self._profile_segments:
            if ps.start_dist <= d < ps.end_dist:
                match = ps
                break
        if match is None:
            return None
        ps = match
        sampled_profile = get_sampled_profile(ps)
        d -= ps.start_dist
        idx = math.floor(d / self._sampling_step)
        mix = (d - idx * self._sampling_step) / self._sampling_step
        val1 = sampled_profile[idx]
        if mix > 0 and idx < len(sampled_profile) - 1:
            val2 = sampled_profile[idx+1]
            return val1 * (1 - mix) + val2 * mix
        return val1

    def get_elevation(self, ref_track_dist):
        return self._get_datapoint(lambda ps: ps.profile, ref_track_dist)

    def seg_count(self):
        return len(self._profile_segments)

    def apply_fixpoint_corrections(self, fp_db, window_length):
        if len(fp_db) == 0:
            return

        fp_minmax_pos_list = []
        for ps in self._profile_segments:
            fp_minmax_pos_list += ElevationProfile._match_minmax_fixpoints(ps, fp_db, self._ref_track, window_length, self._sampling_step)
        fp_minmax_pos_list = sorted(fp_minmax_pos_list, key=lambda x: x[0])
        for match in fp_minmax_pos_list:
            ps = match[3]
            fp = match[2]
            print(f'dist={round(ps.start_dist + match[0] * self._sampling_step, 1)} elevation difference={round(match[1], 1)} name={fp.name}')

        vcorr = generate_correction_curve(fp_minmax_pos_list, self._ref_track.distance(), self._sampling_step)
        self._apply_vcorr(vcorr)

        fp_horiz_pos_list = []
        for ps in self._profile_segments:
            fp_horiz_pos_list += ElevationProfile._match_slope_fixpoints(ps, fp_db, self._ref_track, window_length, self._sampling_step)

        fp_horiz_pos_list = sorted(fp_horiz_pos_list, key=lambda x: x[0])
        for match in fp_horiz_pos_list:
            ps = match[3]
            fp = match[2]
            print(f'dist={round(ps.start_dist + match[0] * self._sampling_step,1)} elevation={round(fp.elevation,1)} horizontal distance to fix={round(match[1],1)} name={fp.name}')

        hcorr = generate_correction_curve(fp_horiz_pos_list, self._ref_track.distance(), self._sampling_step)
        self._apply_hcorr(hcorr)
        self._calculate_lowpass_profiles()

    @staticmethod
    def close_plot_files():
        for f in ElevationProfile._plot_files.values():
            f.close()

    @staticmethod
    def _open_plot_files(plot_dir):
        if len(ElevationProfile._plot_files) > 0:
            return
        names = ['op', 'fp']
        for name in names:
            ElevationProfile._plot_files[name] = open(os.path.join(plot_dir, name), 'w')

    def write_plot_data(self, plot_dir):
        ElevationProfile._open_plot_files(plot_dir)
        for ps in self._profile_segments:
            f = ElevationProfile._plot_files['op']
            for idx, ele in enumerate(ps.profile):
                f.write(f'{ps.start_dist+idx*self._sampling_step}, {ele}\n')
            f.write('\n\n')
            f.flush()
            f = ElevationProfile._plot_files['fp']
            for idx, ele in enumerate(ps.lowpass_profile):
                f.write(f'{ps.start_dist+idx*self._sampling_step}, {ele}\n')
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
            current_minmax = (None, None, None)
            for idx, items in mmpos.items():
                for (minmax, ele) in items:
                    if current_minmax[0] is None or idx >= current_minmax[0] + window_len:
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
    def _sort_point_clusters(clustered_mmpos, min_dist_between_fixpoints, processed_set):
        vertical_span = 12.0
        horizontal_spans = [30, 60, 120]
        minmax_values = []
        for idx, items in clustered_mmpos.items():
            closest_idx = find_closest(processed_set, idx)
            if closest_idx is not None and abs(closest_idx - idx) < min_dist_between_fixpoints:
                continue
            minmax = items[0][1]
            eles = [i[2] for i in items]
            val, count, span = get_clustered_average(eles, [vertical_span])
            if count == 1:
                continue
            idxs = [i[0] for i in items]
            mid_idx, _, _ = get_clustered_average(idxs, horizontal_spans)
            minmax_values.append((mid_idx, minmax, val, count, span, idx))

        def cmp_minmax(a, b):
            if a[4] <= vertical_span and b[4] <= vertical_span:
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
    def estimate_fixpoints_and_apply_corrections(profiles, window_length, plot_dir=None):
        ref_track = profiles[0]._ref_track
        sampling_step = profiles[0]._sampling_step

        min_dist_between_fixpoints = math.ceil(8 * window_length / sampling_step)
        window_length_for_clusters = math.ceil(window_length / sampling_step)
        processed_set = SortedDict()

        mmpos = ElevationProfile._index_local_minmax_points(profiles)
        clustered_mmpos = ElevationProfile._make_point_clusters(mmpos, window_length_for_clusters)

        sorted_values = ElevationProfile._sort_point_clusters(clustered_mmpos, min_dist_between_fixpoints, processed_set)

        spaced_out_values = []
        picked_points = SortedDict()
        for val in sorted_values:
            mid_idx = val[0]

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
                mmlist = ElevationProfile._match_minmax_fixpoints(ps, fp_db, ref_track, window_length, profile._sampling_step, lowpass=True)
                fp_minmax_pos_list += mmlist

                for val in mmlist:
                    fp = val[2]
                    src_idx = val[0]
                    corr_dist = src_idx * profile._sampling_step - fp.dist
                    fp_horiz_pos_list.append((src_idx, corr_dist, fp))

            fp_minmax_pos_list = sorted(fp_minmax_pos_list, key=lambda x: x[0])
            fp_horiz_pos_list = sorted(fp_horiz_pos_list, key=lambda x: x[0])
            #print(list(map(lambda i: (i[0], i[1]), fp_minmax_pos_list)))
            #print(list(map(lambda i: (i[0], i[1]), fp_horiz_pos_list)))
            vcorr = generate_correction_curve(fp_minmax_pos_list, ref_track.distance(), profile._sampling_step)
            hcorr = generate_correction_curve(fp_horiz_pos_list, ref_track.distance(), profile._sampling_step)
            ElevationProfile._apply_vcorr(profile, vcorr)
            ElevationProfile._apply_hcorr(profile, hcorr)
            ElevationProfile._calculate_lowpass_profiles(profile)
        return True
