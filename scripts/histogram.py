#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import ruamel.yaml as yaml
import sys
import os

import ccset_pb2

parser = argparse.ArgumentParser(description='Plot histograms of pbdats.')
parser.add_argument('--save-plots', type=bool, nargs='?', const=True, default=False, help='save images of all histograms')
parser.add_argument('--plot-window', type=bool, nargs='?', const=True, default=False, help='plot integration window')
parser.add_argument('--plot-window-custom', type=bool, nargs='?', const=True, default=False, help='plot integration window, specify window in code')
parser.add_argument('--h-wnd', type=float, nargs=2, default=[-1.0,-1.0], help='H channel integration window', required=False)
parser.add_argument('--v-wnd', type=float, nargs=2, default=[-1.0,-1.0], help='V channel integration window', required=False)
args=parser.parse_args()

print('save plots: ', args.save_plots)
print('plot integration window: ', args.plot_window)
print('h window: ', args.h_wnd)
print('v window: ', args.v_wnd)

if (args.plot_window):
    if not (args.plot_window_custom):
        if (args.h_wnd[0] == -1) or (args.h_wnd[1] == -1):
            print('please specify integration window for h channel')
            sys.exit()
        if (args.v_wnd[0] == -1) or (args.v_wnd[1] == -1):
            print('please specify integration window for v channel')
            sys.exit()

DATAEXT='.pbdat'
def readtags(fn):
    #print(fn)
    tags=pd.read_csv(fn, sep="\t", header=None, names=["channel", "tag"])
    stags=pd.concat([tags[tags.channel==21],tags[tags.channel==29]])  # tags of signal photon
    ttags=tags[tags.channel==49] # tags of trigger photon

    stags.reset_index(drop=True,inplace=True)
    ttags.reset_index(drop=True,inplace=True)
    starttag=min(ttags.tag[0],stags.tag[0])
    stags.tag-=starttag
    ttags.tag-=starttag
    stags.tag*=cs
    ttags.tag*=cs
    stags.tag=stags.tag.round()
    ttags.tag=ttags.tag.round()
    
    return ttags, stags

def loadcc(fn):
    if os.path.exists(fn):
        data=np.load(fn.replace('.txt.gz','_ccs.npz'))
        ccs_h_1ns = data['ccs_h']
        ccs_v_1ns = data['ccs_v']
        offsets_1ns = data['offsets']
    else:
        print('file does not exist')
    return ccs_h_1ns, ccs_v_1ns, offsets_1ns

def loadh5(fname):
    try:
        read_data = ccset_pb2.ccset_data()
        f=open(fname, "rb")
        read_data.ParseFromString(f.read())
        f.close()
        
        offsets = np.array(read_data.offsets)
        cc_h = np.array(read_data.cc_h)
        cc_v = np.array(read_data.cc_v)
        
        #cc_h_tags=np.array([e.arr for e in read_data.cc_h_tags if len(e.arr)>0], dtype=object)
        #cc_v_tags=np.array([e.arr for e in read_data.cc_v_tags if len(e.arr)>0], dtype=object)
        
        #meastime = read_data.meastime
        
        #num_trigger_tags = read_data.num_trigger_tags
        #num_fpga_tags = read_data.num_fpga_tags
        
        #return offsets, cc_h, cc_v, cc_h_tags, cc_v_tags, num_trigger_tags, num_fpga_tags, meastime
        
        return cc_h, cc_v, offsets
    except Exception as e:
        print(fname)
        print(e)
    
fnames = [e for e in os.listdir() if e.endswith(DATAEXT)]
exclusions = ['210824', '210825_0', '210825_10']
for ex in exclusions:
    fnames = [fn for fn in fnames if not ex in fn]

delays=[e for e in [fn.split('_')[3] for fn in fnames]]
fnames=[y for x,y in sorted(zip(delays,fnames))]
delays.sort()

maxima=[]
maxima_delays=[]
i=0
for delay,fname in zip(delays,fnames):
    offsets=np.arange(2690,2700)
    i+=1
    print('{0:d}/{1:d}: {2:s}'.format(i,len(fnames),fname))
    if DATAEXT=='.pbdat':
        ccs_h_1ns_full, ccs_v_1ns_full, offsets_1ns_full = loadh5(fname)
    elif DATAEXT=='.npz':
        ccs_h_1ns_full, ccs_v_1ns_full, offsets_1ns_full = loadcc(fname)
    pltrng=[0,-1]
    ccs_h_1ns=ccs_h_1ns_full[pltrng[0]:pltrng[1]]
    ccs_v_1ns=ccs_v_1ns_full[pltrng[0]:pltrng[1]]
    offsets_1ns=offsets_1ns_full[pltrng[0]:pltrng[1]]
    fig,ax = plt.subplots(figsize=(12,8), dpi= 100)
    ax.plot(offsets_1ns,ccs_h_1ns, color='blue', label='h')
    ax.plot(offsets_1ns,ccs_v_1ns, color='orange', label='v')
    plt.legend()
    if args.plot_window:
        if args.plot_window_custom:
            fdate = int(fname.split('_')[0])
            ftime = int(fname.split('_')[1])
            
            hwnd = [4411.4, 4411.6]
            vwnd = [4412.0, 4412.2]
            
            if fdate >= 20210825 and ftime > 80000:
                hwnd = [4411.3, 4411.5]
                vwnd = [4411.9, 4412.1]
            
        else:
            hwnd = args.h_wnd
            vwnd = args.v_wnd
        plt.axvspan(hwnd[0], hwnd[1], color='blue', alpha=0.3)
        plt.axvspan(vwnd[0], vwnd[1], color='orange', alpha=0.3)
    ax.set_title(fname)
    maxima.append(max(ccs_h_1ns))
    maxima_delays.append(np.where(ccs_h_1ns_full==max(ccs_h_1ns))[0][0])
    if args.save_plots:
        if not os.path.exists('histograms'):
            os.makedirs('histograms')
        plt.savefig('histograms/'+fname.replace('pbdat','png'))
    else:
        plt.show()
