from __future__ import print_function
from pixell import enmap,utils, reproject, enplot
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import urllib.request
from astropy.table import QTable
import astropy.units as u
import yaml
from astropy.io import fits
import matplotlib.patheffects as path_effects

def tileFinder(ra, dec, data):
    #Given an RA and Dec in deg, find the S18d tile containing that RA and Dec
    for i, tile in enumerate(data):
        box = tile['RADecSection']
        if box[0] >= box[1]:
            if (360 >= ra >= box[0] or 0 <= ra <box[1]) and box[2]<=dec<= box[3]:
                return tile['tileName']
        if box[0]<=ra<=box[1] and box[2]<=dec<= box[3]:
            return tile['tileName']
    return None

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
def freqStamp(ra, dec, fmap, name, width = 0.5, write = True):
    
    stamp = reproject.postage_stamp(fmap, ra, dec, width*60, 0.5)
    if write:
        stamp.wcs.wcs.crval = [ra,dec]
        enmap.write_map('./for_tony/{}.fits'.format(name), stamp)
    return stamp

def s18dStamp(ra, dec, data, name, width = 0.5, write = True):
    #Find tile corresponding to RA, Dec
    path = '/scratch/r/rbond/jorlo/S18d_202006/filteredMaps/'
    tileName = tileFinder(ra, dec, data)
    if tileName == None: return None
    tile = enmap.read_map(path+tileName+'/Arnaud_M2e14_z0p4#'+tileName+'_filteredMap.fits')
    
    
    stamp = reproject.postage_stamp(tile, ra, dec, width*60, 0.5)
    if write:
        stamp.wcs.wcs.crval = [ra,dec]
        enmap.write_map('./for_tony/{}.fits'.format(name), stamp)
    return stamp

with open('/scratch/r/rbond/jorlo/S18d_202006/selFn/tileDefinitions.yml') as f:
    
    s18d = yaml.load(f)
"""
mdcw_catalog = fits.open('/home/s/sievers/sdicker/ACTnCOWs/MADCOWSUnion.fits')

ra = mdcw_catalog[1].data['RADeg']
names = mdcw_catalog[1].data['name']
dec = mdcw_catalog[1].data['decDeg']
ra, dec = np.array(ra), np.array(dec)
rich = mdcw_catalog[1].data['Rich']
z = mdcw_catalog[1].data['Photz']

flags = np.where((30<= rich) & (rich <=40))[0]

ra, dec, names, rich = ra[flags], dec[flags], names[flags], rich[flags]

lowest_y = []
central_y = []

for i in range(len(ra)):
    stamp = s18dStamp(ra[i], dec[i], s18d, 'fakename', width = 0.5, write = False)
    if stamp is None: 
        lowest_y.append(0)
        central_y.append(0)
        continue
    lowest_y.append(np.amin(stamp[0][29:31, 29:31]))
    central_y.append(np.mean(stamp[0][29:31, 29:31]))

n=10

central_idx = np.argpartition(central_y, n)
central_idx = central_idx[:n]
lowest_idx = np.argpartition(lowest_y, n)
lowest_idx = lowest_idx[:n]

path = '/home/r/rbond/sigurdkn/project/actpol/map_coadd/20200228/release2/'


maps = ['220', '150', '090']

for i in central_idx:
    for imap in maps:

        cur_map = enmap.read_map(path + 'act_planck_s08_s18_cmb_f{}_night_map.fits'.format(imap))

        #print(str(imap)+names[i].replace(' ', '_'))
        stamp = freqStamp(ra[i], dec[i], cur_map, imap+names[i].replace(' ', '_'))
        
        print(i)
        print(len(names))
        print(len(ra))     
        print(names[i])

        fig = plt.figure()

        plot = plt.imshow(stamp[0], extent = [-30, 30, -30, 30])
        plt.scatter(0,0, marker = '+', color = 'r', alpha = 0.5)
        plt.colorbar(plot, format='%.0e')
        plt.title('Cluster {}, freq: {}GHz'.format(names[i], imap))
        txt = plt.text(-25, 25, r'Central $y_0$ = {:0.3e}\n z: '.format(central_y[i], z[i]), size=11, color='black')
        txt.set_path_effects([path_effects.withStroke(linewidth=5, foreground='w')])
        plt.xlabel('Arcmin')
        plt.ylabel('Arcmin')
        plt.savefig('./plots/low_y/{}_central_{}.pdf'.format(imap, names[i].replace(' ', '_')))
        plt.savefig('./plots/low_y/{}_central_{}.png'.format(imap, names[i].replace(' ', '_'), dpi = 300))
        plt.show()
        plt.close()
"""
#ignore = freqStamp(303.1103611788376,-56.83093699573238, map220, '220_ACT-CL_J2012.4-5649')
#ignore = freqStamp(303.1103611788376,-56.83093699573238, map150, '150_ACT-CL_J2012.4-5649')
#ignore = freqStamp(303.1103611788376,-56.83093699573238, map090, '090_ACT-CL_J2012.4-5649')
maps = ['220', '150', '090']
path = '/home/r/rbond/sigurdkn/project/actpol/map_coadd/20200228/release2/'

dec, ra = -2.38, 33.92
for imap in maps:

    cur_map = enmap.read_map(path + 'act_planck_s08_s18_cmb_f{}_night_map.fits'.format(imap))


    stamp = freqStamp(ra, dec, cur_map, 'ignore')

    fig = plt.figure()

    plot = plt.imshow(stamp[0], extent = [-30, 30, -30, 30])
    plt.scatter(0,0, marker = '+', color = 'r', alpha = 0.5)
    plt.colorbar(plot, format='%.0e')
    plt.title('Test, freq: {}GHz'.format(imap))
    plt.xlabel('Arcmin')
    plt.ylabel('Arcmin')
    plt.savefig('./plots/low_y/{}_central_test.pdf'.format(imap))
    plt.savefig('./plots/low_y/{}_central_test.png'.format(imap, dpi = 300))
    plt.show()
    plt.close()

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

test = enmap.read_map('/home/r/rbond/sigurdkn/project/actpol/map_coadd/20200228/release2/act_planck_s08_s18_cmb_f220_night_map.fits')


