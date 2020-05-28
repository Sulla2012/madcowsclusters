from __future__ import print_function
from pixell import enmap,utils, reproject, enplot
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import urllib.request
from astropy.table import QTable
import astropy.units as u

def stack(ras, decs, map1, mask1, width = 20.):
	stack = 0
	divisor = 0
	for i in range(len(ras)):
		tempdec, tempra = np.deg2rad([decs[i], ras[i]])
		tempwid = np.deg2rad(width/60.)
		box = [[tempdec-tempwid,tempra-tempwid],[tempdec+tempwid,tempra+tempwid]]
		maskstamp = mask1.submap(box)
		#Check if the stamp is entirely within borders, i.e. if the mask stamp is entirely ones
		if np.any(maskstamp[0]):
			stamp = reproject.postage_stamp(map1, ras[i], decs[i], width, 0.5)
			if stamp is None: continue
			stack += stamp[0]
			divisor += 1
	stack /= divisor
	print(divisor)
	return stack

mappath = '/scratch/r/rbond/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/'

t = QTable.read('madcows/AdvACT_S18Clusters_v1.0-beta.fits')

boss_mask = enmap.read_map(mappath + "tilec_mask.fits")
boss_map = enmap.read_map(mappath + 'tilec_single_tile_boss_comptony_map_v1.2.0_joint.fits')

ra_temp = t['RADeg']
dec_temp = t['decDeg']
ra, dec = np.array(ra_temp), np.array(dec_temp)
boss_stack = stack(ra, dec, boss_map, boss_mask)
plt.imshow(boss_stack)
plt.colorbar()
plt.savefig("boss_stack_act_y.png")
plt.close()


#plots = enplot.plot(enmap.upgrade(boss_stack,5),grid=False, colorbar=True,color='gray')
#enplot.write("boss_stack_act_y",plots)


