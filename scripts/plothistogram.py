#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import histogramset_pb2

parser = argparse.ArgumentParser(description='Plot histograms of pbdats.')
parser.add_argument('--save-plots', type=bool, nargs='?', const=True, default=False, help='save images of all histograms')
parser.add_argument('--plot-window', type=float, nargs=2, default=[-1.0,-1.0], help='plot integration window.', required=False)
args = parser.parse_args()

print('save plots: ', args.save_plots)
wnd = args.plot_window
# check if values for window exist
plot_window = True
if (wnd[0] == -1) and (wnd[1] == -1):
    plot_window = False
if plot_window:
    if (wnd[0] == -1) or (wnd[1] == -1):
        print('please specify integration window for integration window')
        sys.exit()

DATAEXT = '.pbdat'


def readpbdat(fname):
    try:
        read_data = histogramset_pb2.histograms()
        f = open(fname, "rb")
        read_data.ParseFromString(f.read())
        f.close()

        meastime = read_data.meastime[0]
        patterns = np.array([e.arr for e in read_data.pattern])
        offsets = np.array([e.arr for e in read_data.offsets])
        ccs = np.array([e.arr for e in read_data.cc])
        cc_tags = np.array([e.arr for e in read_data.cc_tags if len(e.arr) > 0], dtype = object)

        return patterns, offsets, ccs, cc_tags, meastime
    except Exception as e:
        print(fname)
        print(e)


fnames = [e for e in os.listdir() if e.endswith(DATAEXT)]
exclusions = ['210824', '210825_0', '210825_10']
for ex in exclusions:
    fnames = [fn for fn in fnames if ex not in fn]
i = 0
for fname in fnames:
    i += 1
    print('{0:d}/{1:d}: {2:s}'.format(i, len(fnames), fname))
    patterns, offsets, ccs, cc_tags, meastime = readpbdat(fname)

    fig, ax = plt.subplots()
    for j in range(0, len(patterns)):
        patlabel = '{'+', '.join([str(e) for e in patterns[j]])+'}'
        ax.plot(offsets[j], ccs[j], label=patlabel)
    ax.legend()
    ax.set_title(fname)

    # plot window
    if plot_window:
        plt.axvspan(wnd[0], wnd[1], color='blue', alpha=0.3)

    if args.save_plots:
        if not os.path.exists('histograms'):
            os.makedirs('histograms')
        plt.savefig('histograms/'+fname.replace('pbdat', 'png'))
        plt.savefig('histograms/'+fname.replace('pbdat', 'pdf'))
    else:
        plt.show()
