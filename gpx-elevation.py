#!/usr/bin/env python3

import argparse
import pathlib
import math
import os
import sys
import pyjson5
from pygnuplot import gnuplot
from rdp import rdp
from pyproj import Transformer

from twodimsearch import TwoDimSearch

from reference_track import ReferenceTrack
from elevation_profile import ElevationProfile
from misc import apply_lowpass_filter, calculate_elevation_gain, get_clustered_average, write_gpx, FixPoint, Point

def get_elevation(dist, profiles):
    eles = []
    for p in profiles:
        ele = p.get_elevation(dist)
        if ele is not None:
            eles.append(ele)
    ele, _, _ = get_clustered_average(eles, [1.5, 2, 2.5, 3.0, 5.0])
    return ele

def main():

    parser = argparse.ArgumentParser(description='Derive elevation for a GPX track by averaging elevation from recorded GPX activities')
    parser.add_argument('--track', type=pathlib.Path, help='Reference track to calculate elevation profile for.', required=True)
    parser.add_argument('--activities_dir', type=pathlib.Path, help='Path to directory with recorded activities in GPX or CSV format.', required=True)
    parser.add_argument('--fix_points', type=pathlib.Path, help='Path to JSON file containing provided fix points.')
    parser.add_argument('--plot_dir', type=pathlib.Path, help='Directory to store plot files suitable for gnuplot, and will also activiate plotting during the run.')
    parser.add_argument('--output', type=pathlib.Path, help='Name of output GPX file. Default: output.gpx', default='output.gpx')
    parser.add_argument('--side_scan_dist', type=float, help='Limit of how far to the side of the track in meters to find a point (on another track) and consider it to be a match. Default: 20', default=20.0)
    parser.add_argument('--sampling_step', type=float, help='Span in meters between each point a track is sampled. Default: 5', default=5.0)
    parser.add_argument('--filter_width', type=float, help='Filter width in meters for the lowpass filter applied before final output. Default: 100', default=100.0)
    parser.add_argument('--simplify_max_error', type=float, help='Maximum error in meters that may be introduced when simplifying (reducing number of points) in the final GPX, this applies in 3D. Default: 0.15', default=0.15)

    args = parser.parse_args()

    do_plot = args.plot_dir is not None
    if do_plot:
        gplot = gnuplot.Gnuplot()
        gplot.cmd(f"cd '{args.plot_dir}'")
        screen_width = 1920
        gplot.cmd(f'set terminal qt size {screen_width},{round(screen_width/3)}')
    else:
        gplot = None

    ref_track = ReferenceTrack(args.track)

    fp_db = TwoDimSearch()
    if args.fix_points is not None:
        print('Loading fix points (native coordinates, meter elevation)')
        json = pyjson5.load(open(args.fix_points, 'r'))
        if 'projection' in json:
            src_proj = json['projection']
        else:
            src_proj = 'epsg:4326'
        fp_to_utm = Transformer.from_crs(src_proj, ref_track.utm_epsg_name())

        for fp in json['fix_points']:
            x, y = fp_to_utm.transform(float(fp[0]), float(fp[1]))
            ele = fp[2]
            fixtype = fp[3]
            name = fp[4]
            point = Point(x, y)
            fixpoint = FixPoint(x, y, ele, fixtype, name)
            fp_db.insert(point, fixpoint)

    if do_plot:
        fp_list = []
        for fpset in fp_db:
            fp = next(iter(fpset))
            point = ref_track.get_nearest_perpendicular_point(fp, args.side_scan_dist)
            if point is not None:
                fp_list.append((point.dist, fp.elevation))
        fp_list = sorted(fp_list, key=lambda x: x[0])
        with open(os.path.join(args.plot_dir, 'fixpoints'), 'w') as f:
            for val in fp_list:
                f.write(f'{val[0]} {val[1]} 100\n')

    profiles = []
    seg_base = 0
    seg_count = 0
    for filename in sorted(list(os.listdir(args.activities_dir))):
        f = os.path.join(args.activities_dir, filename)
        if os.path.isfile(f):
            profile = ElevationProfile(f)
            profile.set_reference_track(ref_track, args.side_scan_dist, args.sampling_step)
            seg_count += profile.seg_count()
            if do_plot:
                profile.write_plot_data(args.plot_dir)
                plot_cmd = f"plot for [IDX=0:{seg_count-1}] 'op' index IDX u 1:2 with l lt IDX+1 notitle"
                if len(fp_db) > 0:
                    plot_cmd += ", 'fixpoints' w xerrorbars lc 'red' title 'provided fixpoints'"
                gplot.cmd(plot_cmd)
            profiles.append(profile)
    seg_end = seg_count

    did_apply_estimated = ElevationProfile.estimate_fixpoints_and_apply_corrections(profiles, args.plot_dir)
    if did_apply_estimated:
        seg_base = seg_count
        seg_end = seg_count*2-1

    if len(fp_db) > 0:
        print('Applying provided fix points')
        for profile in profiles:
            profile.apply_fixpoint_corrections(fp_db)

    if do_plot:
        for profile in profiles:
            profile.write_plot_data(args.plot_dir)
        plot_cmd = f"plot \
        for [IDX={seg_base}:{seg_end}] 'fp' index IDX u 1:2 with l lt IDX+1 lw 2 notitle, \
        for [IDX={seg_base}:{seg_end}] 'op' index IDX u 1:2 with l lt IDX+1 notitle"
        if len(fp_db) > 0:
            plot_cmd += ", 'fixpoints' w xerrorbars lc 'red' title 'provided fix points'"
        if did_apply_estimated:
            plot_cmd += ", 'derived_minmax' w yerrorbars title 'derived fix points'"
        gplot.cmd(plot_cmd)

    ElevationProfile.close_plot_files()

    profile = []
    sample_count = math.ceil(ref_track.distance() / args.sampling_step)
    start_dist = None
    end_dist = None
    for idx in range(sample_count+1):
        dist = idx / sample_count * ref_track.distance()
        ele = get_elevation(dist, profiles)
        # we allow missing samples at start and end (which is normal due to horizontal
        # alignment corrections), but not gaps in the middle
        if ele is None:
            if start_dist is not None and end_dist is None:
                end_dist = dist
        elif end_dist is not None:
            print('There is not elevation data to cover the complete track, aborting')
            sys.exit(1)
        else:
            profile.append(ele)
            if start_dist is None:
                start_dist = dist

    lowpass_profile = apply_lowpass_filter(profile, args.sampling_step / args.filter_width)
    if do_plot:
        with open(os.path.join(args.plot_dir, 'profile0'), 'w') as f:
            for idx, ele in enumerate(profile):
                f.write(f'{start_dist+idx*args.sampling_step}, {ele}\n')
        with open(os.path.join(args.plot_dir, 'profile1'), 'w') as f:
            for idx, ele in enumerate(lowpass_profile):
                f.write(f'{start_dist+idx*args.sampling_step}, {ele}\n')
        plot_cmd = f"plot \
        for [IDX={seg_base}:{seg_end}] 'op' index IDX u 1:2 with l lt IDX+1 notitle, \
        'profile1' w l lc rgb '#88333333' lw 4 title 'finished profile', \
        'profile0' w l lc rgb '#333333' title 'raw profile'"
        if len(fp_db) > 0:
            plot_cmd += ", 'fixpoints' w xerrorbars lc 'red' title 'provided fix points'"
        if did_apply_estimated:
            plot_cmd += ", 'derived_minmax' w yerrorbars title 'derived fix points'"
        gplot.cmd(plot_cmd)

    print('Deriving elevation for each point using database')
    rdp_list = []
    distmap = {}
    for tp in ref_track:
        dist = tp.dist - start_dist
        idx = math.floor(dist / args.sampling_step)
        if idx < len(lowpass_profile):
            ele = lowpass_profile[idx]
            mix = (dist - idx * args.sampling_step) / args.sampling_step
            if mix > 0 and idx < len(lowpass_profile) - 1:
                ele2 = lowpass_profile[idx+1]
                ele = ele * (1 - mix) + ele2 * mix
            rdp_list.append([tp.x, tp.y, ele])
            distmap[(tp.x, tp.y, ele)] = tp.dist

    print('Simplify the 3D polyline')
    rdp_result = rdp(rdp_list, epsilon=args.simplify_max_error)
    print(f'Reduced number of points from {len(ref_track)} to {len(rdp_result)}')

    elevation_gain = calculate_elevation_gain(rdp_result, lambda x: x[2])
    print(f'Profile elevation gain: {round(elevation_gain,1)} meters for simplified profile')

    if do_plot:
        with open(os.path.join(args.plot_dir, 'profile2'), 'w') as f:
            for point in rdp_result:
                dist = distmap[(point[0], point[1], point[2])]
                ele = round(point[2], 1)
                f.write(f'{dist}, {ele}\n')
        plot_cmd = f"plot \
        for [IDX={seg_base}:{seg_end}] 'op' index IDX u 1:2 with l lt IDX+1 notitle, \
        'profile1' w l lc rgb '#88333333' lw 4 title 'finished profile', \
        'profile0' w l lc rgb '#333333' title 'raw profile', \
        'profile2' w l lc rgb '#000000' dashtype 2 title 'GPX stored profile'"
        if len(fp_db) > 0:
            plot_cmd += ", 'fixpoints' w xerrorbars lc 'red' title 'provided fix points'"
        if did_apply_estimated:
            plot_cmd += ", 'derived_minmax' w yerrorbars title 'derived fix points'"
        gplot.cmd(plot_cmd)

    print(f'Writing GPX to {args.output}')
    write_gpx(args.output, rdp_result, Transformer.from_crs(ref_track.utm_epsg_name(), 'epsg:4326').transform)

    if do_plot:
        input('Complete! Press Enter to close.')

if __name__ == "__main__":
    main()
