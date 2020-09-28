from __future__ import print_function
from pixell import enmap,utils, reproject, enplot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os,sys
import urllib.request
from scipy import interpolate
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
import csv
import yaml
from scipy.stats import norm
import matplotlib.patheffects as path_effects
import pickle as pk
from matplotlib.ticker import Locator

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

def s18dStamp(ra, dec, data, name, width = 0.5, write = True):
    #Find tile corresponding to RA, Dec

    path = '/project/r/rbond/jorlo/S18d_202006/filteredMaps/'
    tileName = tileFinder(ra, dec, data)
    if tileName == None: return None
    tile = enmap.read_map(path+tileName+'/Arnaud_M2e14_z0p4#'+tileName+'_filteredMap.fits')
    
    
    stamp = reproject.postage_stamp(tile, ra, dec, width*60, 0.5)
    
    if write:
        stamp.wcs.wcs.crval = [ra,dec]
        enmap.write_map('./for_tony/{}.fits'.format(name), stamp)
    return stamp


#########################################################################################
#                          220 Stack over all ACT Clusters                              #
########################################################################################

"""
act_catalog = fits.open('/gpfs/fs0/project/r/rbond/jorlo/cluster_catalogs/DR5_cluster-catalog_v1.0b2.fits')
cur_map = enmap.read_map('/gpfs/fs0/project/r/rbond/jorlo/freq_maps/stitched_Beam220_filteredMap.fits')

ra = act_catalog[1].data['RADeg']
names = act_catalog[1].data['name']
dec = act_catalog[1].data['decDeg']
ra, dec = np.array(ra), np.array(dec)
mass = act_catalog[1].data['M500']

names_cut = {'ACT-CL J0015+0801','ACT-CL J0917+1456','ACT-CL J0936+0336','ACT-CL J1116+1653','ACT-CL J1350+0036','ACT-CL J1355-0114','ACT-CL J2349+0541'}

names_flag = [True]*len(names)

for i, name in enumerate(names):
    if name in names_cut:
        names_flag[i] = False

        
ra = ra[names_flag]
dec = dec[names_flag]
mass = mass[names_flag]
names = names[names_flag]
 

central_t_boot = []

for j in range(200):
    print(j, end='\r')
    flags = np.random.randint(len(ra), size = len(ra))
    
    ra_temp, dec_temp = ra[flags], dec[flags]
    
    stack = 0
    divisor = 0

    for k in range(len(ra_temp)):
        stamp = reproject.postage_stamp(cur_map, ra_temp[k], dec_temp[k], 20., 0.5)
        if stamp is None: continue
        if 0 in stamp[0][19:21, 19:21]: print('skip')
        if stamp[0][20, 20] == 0: print('skip2')
        stack += stamp[0]
        divisor += 1

    stack /= divisor
    central_t_boot.append(np.mean(stack[19:21, 19:21]))
    
print('220boot = {} +/- {}'.format(np.mean(central_t_boot), np.sqrt(np.var(central_t_boot))))



###############################################################################################################
#                         220 Signal stacked on MDCWs positions binned in richness                           #
##############################################################################################################

mdcw_catalog = fits.open('/home/s/sievers/sdicker/ACTnCOWs/MADCOWSUnion.fits')

ra = mdcw_catalog[1].data['RADeg']
names = mdcw_catalog[1].data['name']
dec = mdcw_catalog[1].data['decDeg']
ra, dec = np.array(ra), np.array(dec)
rich = mdcw_catalog[1].data['Rich']
z = mdcw_catalog[1].data['Photz']

names_cut = {'MOO J0015+0801','MOO J0917+1456','MOO J0936+0336','MOO J1116+1653','MOO J1350+0036','MOO J1355-0114','MOO J2349+0541'}

names_flag = [True]*len(names)

for i, name in enumerate(names):
    if name in names_cut:
        names_flag[i] = False

        
ra = ra[names_flag]
dec = dec[names_flag]
rich = rich[names_flag]
names = names[names_flag]
z = z[names_flag]

mdcw_freq_dict = {'090':{}, '150':{}, '220':{}}

boot = True

for i, (freq, cur_dict) in enumerate(mdcw_freq_dict.items()):    
    
    cur_dict['central_t'] = []
    cur_dict['map_var'] = []
    cur_dict['boots'] = []
    
    cur_map = enmap.read_map('/gpfs/fs0/project/r/rbond/jorlo/freq_maps/stitched_Beam{}_filteredMap.fits'.format(str(freq)))
    
    for j in range(0,70, 10):        
        flag = np.where((rich>j) & (rich<=j+10))[0]

        if len(flag) == 0: continue

        ra_temp, dec_temp = ra[flag], dec[flag]
        if boot:
            cur_boot = []
            for k in range(50):
                print(freq, j, k, end = '\r')
                
                flags = np.random.randint(len(ra_temp), size = len(ra_temp))
    
                ra_temp2, dec_temp2 = ra_temp[flags], dec_temp[flags]
    
                stack = 0
                divisor = 0

                for l in range(len(ra_temp)):
                    stamp = reproject.postage_stamp(cur_map, ra_temp2[l], dec_temp2[l], 20., 0.5)
                    if stamp is None: continue
                    if 0 in stamp[0][19:21, 19:21]: continue
        
                    stack += stamp[0]
                    divisor += 1

                stack /= divisor
                cur_boot.append(np.mean(stack[19:21, 19:21]))
            cur_dict['boots'].append(cur_boot)
            cur_dict['central_t'].append(np.mean(cur_boot))
            cur_dict['map_var'].append(np.std(cur_boot))	

        else:
            stack = 0
            divisor = 0

            for k in range(len(ra_temp)):
                stamp = reproject.postage_stamp(cur_map, ra_temp[k], dec_temp[k], 20., 0.5)
                if stamp is None: continue
                if 0 in stamp[0][19:21, 19:21]: continue
                stack += stamp[0]
                divisor += 1

            stack /= divisor

            plot = plt.imshow(stack, extent = [-20, 20, -20, 20])
            plt.scatter(0,0, marker = '+', color = 'r', alpha = 0.5)
            plt.colorbar(plot, format='%.0e')
            plt.title('{} Filtered Maps on all MaDCoWS Positions \n {} < Richness <= {}, Num = {}'.format(freq,j, j+10,len(flag)))

            plt.xlabel('Arcmin')
            plt.ylabel('Arcmin')
            plt.savefig('./plots/PS_freq/mdcw_{}_{}_{}_filtered_mdcw.pdf'.format(j, j+10, freq))
            plt.savefig('./plots/PS_freq/mdcw_{}_{}_{}_filtered_mdcw.png'.format(j, j+10,freq), dpi = 300)
            plt.show()
            plt.close()

            cur_dict['highest_t'].append(np.amin(stack[19:21,19:21]))
            cur_dict['central_t'].append(np.mean(stack[19:21, 19:21]))

            cur_dict['map_var'].append(np.std(stack[0:15, 0:15]))

pk.dump(mdcw_freq_dict, open('ps_dict.p', 'wb'))	

xrange = range(5,75, 10)


for i, (freq, cur_dict) in enumerate(mdcw_freq_dict.items()):
    print(cur_dict['central_t'])
    plt.errorbar(xrange, cur_dict['central_t'], yerr = cur_dict['map_var'], label = freq,fmt='--o', linestyle = 'none')

plt.legend(loc = 3)
plt.title(r'090/150/220GHz MaDCoWS Stacks')
plt.xlabel('Richness Bin (midpoint)')
plt.ylabel(r'Central 2x2 Mean Temperature ($\mu$K)')
plt.ylim(-60,50)
plt.axhline(y=0, color='k')

plt.savefig('./plots/PS_freq/mdcw_rich_central_temp.pdf')
plt.savefig('./plots/PS_freq/mdcw_rich_central_temp.png', dpi = 300)
plt.show()
plt.close()
"""
##########################################################################################################################################
#                                                   ACT 220 Stacks, binned in richness                                                  #
#########################################################################################################################################

act_catalog = fits.open('/gpfs/fs0/project/r/rbond/jorlo/cluster_catalogs/DR5_cluster-catalog_v1.0b2.fits')

ra = act_catalog[1].data['RADeg']
names = act_catalog[1].data['name']
dec = act_catalog[1].data['decDeg']
ra, dec = np.array(ra), np.array(dec)
mass = act_catalog[1].data['M500']

flag = np.where((rich != 999999))[0]
flag2 = np.where((rich == 999999))
ra2, dec2 = ra[flag2], dec[flag2]
ra, dec, rich = ra[flag], dec[flag], rich[flag]

nbins = 7
perc = np.percentile(rich, np.linspace(0,100, nbins+1))  

perc[0] = perc[0]*0.99

act_freq_dict = {'090':{}, '150':{}, '220':{}}

boot = True

for i, (freq, cur_dict) in enumerate(act_freq_dict.items()):

    cur_dict['central_t'] = []
    cur_dict['map_var'] = []
    cur_dict['boots'] = []

    cur_map = enmap.read_map('/gpfs/fs0/project/r/rbond/jorlo/freq_maps/stitched_Beam{}_filteredMap.fits'.format(str(freq)))

    for j in range(nbins):
        flag = np.where((mass>perc[j]) & (mass<=perc[j+1]))[0]


        if len(flag) == 0: continue

        ra_temp, dec_temp = ra[flag], dec[flag]
        if boot:
            cur_boot = []
            for k in range(40):
                print(freq, j, k, end = '\r')

                flags = np.random.randint(len(ra_temp), size = len(ra_temp))

                ra_temp2, dec_temp2 = ra_temp[flags], dec_temp[flags]

                stack = 0
                divisor = 0

                for l in range(len(ra_temp)):
                    stamp = reproject.postage_stamp(cur_map, ra_temp2[l], dec_temp2[l], 20., 0.5)
                    if stamp is None: continue
                    if 0 in stamp[0][19:21, 19:21]: continue

                    stack += stamp[0]
                    divisor += 1

                stack /= divisor
                cur_boot.append(np.mean(stack[19:21, 19:21]))
            cur_dict['boots'].append(cur_boot)
            cur_dict['central_t'].append(np.mean(cur_boot))
            cur_dict['map_var'].append(np.std(cur_boot))

        else:
            stack = 0
            divisor = 0

            for k in range(len(ra_temp)):
                stamp = reproject.postage_stamp(cur_map, ra_temp[k], dec_temp[k], 20., 0.5)
                if stamp is None: continue
                if 0 in stamp[0][19:21, 19:21]: continue
                stack += stamp[0]
                divisor += 1

            stack /= divisor

            plot = plt.imshow(stack, extent = [-20, 20, -20, 20])

            plt.scatter(0,0, marker = '+', color = 'r', alpha = 0.5)
            plt.colorbar(plot, format='%.0e')
            plt.title('{} Filtered Maps on all ACT Positions \n {:0.3e} < Mass (M500) <= {:0.3e}, Num = {}'.format(freq,perc[j], perc[j+1],len(flag)))


            plt.xlabel('Arcmin')
            plt.ylabel('Arcmin')
            plt.savefig('./plots/PS_freq/act_{}_{}_{}.pdf'.format(j, j+10, freq))
            plt.savefig('./plots/PS_freq/act_{}_{}_{}.png'.format(j, j+10,freq), dpi = 300)
            plt.show()
            plt.close()

            cur_dict['highest_t'].append(np.amin(stack[19:21,19:21]))
            cur_dict['central_t'].append(np.mean(stack[19:21, 19:21]))

            cur_dict['map_var'].append(np.std(stack[0:15, 0:15]))


pk.dump(act_freq_dict, open('act_ps_dict.p', 'wb'))

xrange = range(5,75, 10)

for i, (freq, cur_dict) in enumerate(act_freq_dict.items()):
    print(cur_dict['central_t'])
    plt.errorbar(xrange, cur_dict['central_t'], yerr = cur_dict['map_var'], label = freq,fmt='--o', linestyle = 'none')

plt.legend(loc = 3)
plt.title(r'090/150/220GHz ACT Stacks')
plt.xlabel('Mass Bin (M500, midpoint)')
plt.ylabel(r'Central 2x2 Mean Temperature ($\mu$K)')
plt.ylim(-140,20)
plt.axhline(y=0, color='k')

plt.savefig('./plots/PS_freq/act_rich_central_temp.pdf')
plt.savefig('./plots/PS_freq/act_rich_central_temp.png', dpi = 300)
plt.show()
plt.close()

"""
##########################################################################################################################################
#                                                        MDCW 220 Stack                                                                 #
#########################################################################################################################################

cur_map = enmap.read_map('/gpfs/fs0/project/r/rbond/jorlo/freq_maps/stitched_Beam220_filteredMap.fits')
try:
    central_t_boot = pk.load(open('central_t_mdcw.pk', 'rb'))
except:
    central_t_boot = []

mdcw_catalog = fits.open('/home/s/sievers/sdicker/ACTnCOWs/MADCOWSUnion.fits')

hdu = fits.open('/gpfs/fs0/project/r/rbond/jorlo/freq_maps/stitched_Beam220_filteredMap.fits')

ra = mdcw_catalog[1].data['RADeg']
names = mdcw_catalog[1].data['name']
dec = mdcw_catalog[1].data['decDeg']
ra, dec = np.array(ra), np.array(dec)
rich = mdcw_catalog[1].data['Rich']
z = mdcw_catalog[1].data['Photz']

names_cut = {'MOO J0015+0801','MOO J0917+1456','MOO J0936+0336','MOO J1116+1653','MOO J1350+0036','MOO J1355-0114','MOO J2349+0541'}

names_flag = [True]*len(names)

for i, name in enumerate(names):
    if name in names_cut:
        names_flag[i] = False

        
ra = ra[names_flag]
dec = dec[names_flag]
rich = rich[names_flag]
names = names[names_flag]
z = z[names_flag]

for j in range(100):
    print(j,end='\r')
    flags = np.random.randint(len(ra), size = len(ra))
    
    ra_temp, dec_temp = ra[flags], dec[flags]
    
    stack = 0
    divisor = 0

    for k in range(len(ra_temp)):
        stamp = reproject.postage_stamp(cur_map, ra_temp[k], dec_temp[k], 20., 0.5)

        if stamp is None: continue
        if 0 in stamp[0][19:21, 19:21]: 
                continue
        
        stack += stamp[0]
        divisor += 1

    stack /= divisor
    central_t_boot.append(np.mean(stack[19:21, 19:21]))
    
#print('220jk = {} +/- {}'.format(np.mean(central_t_jk), np.var(central_t_jk)))
print('220boot = {} +/- {}'.format(np.mean(central_t_boot), np.std(central_t_boot)))
"""
