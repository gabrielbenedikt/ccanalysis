#!/usr/bin/env python3

import pandas as pd
import os
import time

def getnewtagfiles():
	fnames = [e for e in os.listdir() if e.endswith('_tags.txt')]
	analyzed = [e for e in os.listdir() if e.endswith('_tags.h5')]
	for a in analyzed:
		fn = a.replace('.h5','.txt')
		if (fn in fnames):
			fnames.remove(fn)
	return fnames

fnames = getnewtagfiles()
print(fnames)

for fn in fnames:
	print(fn)
	t0=time.time()
	tags=pd.read_csv(fn, sep="\t", header=None, names=["channel", "tag"]).dropna().astype('int64')
	t1=time.time()
	print('time to read: ', t1-t0)
	
	t0=time.time()
	tags.to_hdf(fn.replace('_tags.txt','_tags.h5'), key='tags' ,complib='zlib', complevel=9)
	t1=time.time()
	print('time to write: ', t1-t0)
	
	try:
		print('remove file')
		os.remove(fn)
	except Exception as e:
		print(e)
		
	
