from pixell import enmap,utils
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
import os,sys
import urllib.request
from astropy.io import fits

def mask(imap, ra, dec, width = 0.5, apod_pix = 4):
	ra, dec = np.deg2rad([ra,dec])
	width = np.deg2rad(width)
	box = [[dec-width/2.,ra-width/2.],[dec+width/2.,ra+width/2.]]
	stamp = imap.submap(box)
	taper = enmap.apod(stamp*0+1,apod_pix)
	enmap.write_map('test_taper.fits', taper)
	return taper

mappath = '/scratch/r/rbond/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/'

t = QTable.read('madcows/AdvACT_S18Clusters_v1.0-beta.fits')

boss_mask = enmap.read_map(mappath + "tilec_mask.fits")
boss_map = enmap.read_map(mappath + 'tilec_single_tile_boss_comptony_map_v1.2.0_joint.fits')

ra_temp = t['RADeg']
dec_temp = t['decDeg']
ra, dec = np.array(ra_temp), np.array(dec_temp)

test_taper = mask(boss_map, ra[0], dec[0])

plt.imshow(test_taper)
plt.colorbar()
plt.savefig("test_taper.png")
plt.close()
