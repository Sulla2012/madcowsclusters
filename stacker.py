from __future__ import print_function
from pixell import enmap,utils, reproject, enplot
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import urllib.request
from astropy.table import QTable
import astropy.units as u

def stack(ras, decs, map1, mask1, map2 = None, mask2 = None, width = 20.):
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
		if map2 is None: continue
		#If a second map is to be checked, repeat above with second map
		maskstamp = mask2.submap(box)
		if np.any(maskstamp[0]):
                        print(i)
                        stamp = reproject.postage_stamp(map2, ras[i], decs[i], width, 0.5)
                        if stamp is None: continue
                        stack += stamp[0]
                        divisor += 1

	stack /= divisor
	print(divisor)
	return stack

#Stack on S18 cluster catalog

"""
if 0:
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

	d56_path = '/scratch/r/rbond/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/'

	d56_mask = enmap.read_map(d56_path + "tilec_mask.fits")
	d56_map = enmap.read_map(d56_path + 'tilec_single_tile_deep56_comptony_map_v1.2.0_joint.fits')


	d56_stack = stack(ra, dec, d56_map, d56_mask)
	plt.imshow(d56_stack)
	plt.colorbar()
	plt.savefig("d56_stack_act_y.png")
	plt.close()

	comb_stack = stack(ra, dec, boss_map, boss_mask, d56_map, d56_mask)
	plt.imshow(comb_stack)
	plt.colorbar()
	plt.savefig("comb_stack_act_y.png")
	plt.close()
"""

"""
#Frequency stamping for Luca
mappath = '/home/r/rbond/sigurdkn/project/actpol/map_coadd/20200228/release/'
	
freqs = ['f090', 'f150', 'f220']

for freq in freqs:
	name = 'act_planck_s08_s18_cmb_{}_night_map.fits'.format(freq)
	print(name)		
	imap = enmap.read_map(mappath + name)
	
	ra, dec = np.deg2rad([20.,1.])
	width = np.deg2rad(0.5)

	box = [[dec-width/2.,ra-width/2.],[dec+width/2.,ra+width/2.]]
	stamp = imap.submap(box)

	plots = enplot.plot(stamp,range=300,mask=0)
	enplot.write(name,plots)

	enmap.write_map('./padded_v1/{}.fits'.format(name), stamp)

"""
