from __future__ import print_function
from pixell import enmap,utils, reproject, enplot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os,sys
import urllib.request
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
import csv
import yaml
import pickle as pk

def s18dStack(ras, decs, data, width = 20., weight = False):
        stack = 0
        divisor = 0
        num = 0
        for i in range(len(ras)):
                path = '/scratch/r/rbond/jorlo/S18d_202006/filteredMaps/'
                tileName = tileFinder(ras[i], decs[i], data)
                #print(ras[i], decs[i])

                #print(tileName)
                
                if tileName is None: continue
                    
                tile = enmap.read_map(path+tileName+'/Arnaud_M2e14_z0p4#'+tileName+'_filteredMap.fits')
                stamp = reproject.postage_stamp(tile, ras[i], decs[i], width, 0.5)
                if weight:
                    rms = fits.open('/scratch/r/rbond/jorlo/S18d_202006/selFn/RMSMap_Arnaud_M2e14_z0p4.fits')
                    #print(tileName)
                    for j in range(1,500):
                        if rms[j].header['EXTNAME'] == tileName:
                            break
                    #print(j)
                    
                    rmsmap = rms[j].data
                    #print(ras[i], decs[i])

                    coords = np.deg2rad(np.array((decs[i],ras[i])))
                    ypix,xpix = enmap.sky2pix(tile.shape,tile.wcs,coords)

                    weight = rmsmap[int(ypix),int(xpix)]**2
                    #print(rmsmap[int(ypix),int(xpix)])
                    
                    if weight < 10**-20:
                        continue
                    
                    try: 
                        stack += stamp[0]/weight
                    except:
                        print(tileName)
                        print(tile)
                        print(stamp)
                    num += 1
                    divisor += 1/weight

                else:
                    try: 
                        stack += stamp[0]
                    except:
                        print(tileName)
                        print(tile)
                        print(stamp)

                    divisor += 1
                    num += 1

        
        try:
            stack /= divisor
        except: 
            print("Error: no items in stack")
            return None, None
        print("Number in stack: {}".format(num))
        return stack, num

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
    path = '/scratch/r/rbond/jorlo/S18d_202006/filteredMaps/'
    tileName = tileFinder(ra, dec, data)
    if tileName == None: return None
    tile = enmap.read_map(path+tileName+'/Arnaud_M2e14_z0p4#'+tileName+'_filteredMap.fits')
    
    
    stamp = reproject.postage_stamp(tile, ra, dec, width*60, 0.5)
    if write:
        #tempdec, tempra = np.deg2rad([dec, ra])
        #tempwid = np.deg2rad(width)
        #box = [[tempdec-tempwid,tempra-tempwid],[tempdec+tempwid,tempra+tempwid]]
 
        #stampgeo = tile.submap(box)
        
        box = np.array([[ra-width/2,dec-width/2],[ra+width/2,dec+width/2]]) * utils.degree
        shape, wcs = enmap.geometry(pos=box,res=0.5 * utils.arcmin,proj='car')
        print(shape)
        #print(stampgeo.wcs)
        #print(stamp.wcs)
        stamp.wcs = wcs
        print(stamp.wcs)
        print(stamp[0].shape)
        #plt.imshow(stamp[0])
        #plt.show()
        #Return map
        plot = enplot.plot(stamp,mask=0)
        enplot.show(plot)
        enmap.write_map('./for_tony/{}.fits'.format(name), stamp)
    return stamp


def getFlags(ra, dec, data):
    #Returns indeces of ra, dec that are inside a certain data 
    flags = []
    for i in range(len(ra)):
        if tileFinder(ra[i], dec[i], data) != None:
            flags.append(i)
    return flags


def chunks(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

def binmap(map2d, x, y, rbin):

    '''Bin a 2d map into a 1d symmteric radial profile (i.e. for B_l or radial profile)'''
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2+Y**2)
    d = np.digitize(R.flatten(), rbin)
    rbin1d = np.zeros(len(rbin))
    bin1d = np.zeros(len(rbin))
    binvar = np.zeros(len(rbin))
    binerr = np.zeros(len(rbin))
    for i in range (len(rbin)):
        rbin1d[i] = np.mean(R.flatten()[d==i+1])
        bin1d[i] = np.mean(map2d.flatten()[d==i+1])
        binvar[i] = np.var(map2d.flatten()[d==i+1])
        binerr[i] = np.sqrt(np.var(map2d.flatten()[d==i]) / len(map2d.flatten()[d==i+1]))
    return  rbin1d, bin1d, binvar, binerr

with open('/scratch/r/rbond/jorlo/S18d_202006/selFn/tileDefinitions.yml') as f:
    
    s18d = yaml.load(f)

mdcw_catalog = fits.open('/home/s/sievers/sdicker/ACTnCOWs/MADCOWSUnion.fits')

#Stack on MaDCoWS and ACT centers for known cross matches.
#Stacking on MaDCoWS positions
ra_mdcw = mdcw_catalog[1].data['RADeg']
names_mdcw = mdcw_catalog[1].data['name']
dec_mdcw = mdcw_catalog[1].data['decDeg']
ra_mdcw, dec_mdcw = np.array(ra_mdcw), np.array(dec_mdcw)
rich_mdcw = mdcw_catalog[1].data['Rich']

path = '/scratch/r/rbond/jorlo/S18d_202006/filteredMaps/'

act = fits.open('/home/r/rbond/jorlo/dev/madcowsclusters/AdvACT_S18d_202006_MaDCoWSUnionMatch.fits')
in_act = act[1].data['MADCOWSUnion_name']

names = mdcw_catalog[1].data['name']
flags = [True]*len(ra_mdcw)

for i, name in enumerate(names):
    if name not in in_act:
        flags[i] = False


ra_mdcw = ra_mdcw[flags]
dec_mdcw = dec_mdcw[flags]
rich_mdcw = rich_mdcw[flags]
names_mdcw = names_mdcw[flags]


mdcw_s18d_stack, stack_num = s18dStack(ra_mdcw, dec_mdcw, s18d, weight = True)

t = QTable.read('/scratch/r/rbond/jorlo/AdvACT2.fits')

ra_act = t['RADeg']
dec_act = t['decDeg']
ra_act, dec_act = np.array(ra_act), np.array(dec_act)

path = '/scratch/r/rbond/jorlo/S18d_202006/filteredMaps/'

act = fits.open('/home/r/rbond/jorlo/dev/madcowsclusters/AdvACT_S18d_202006_MaDCoWSUnionMatch.fits')
in_act = act[1].data['Name']

names = t['name']

flags = [True]*len(ra_act)

for i, name in enumerate(names):
    if name not in in_act:
        flags[i] = False


ra_act = ra_act[flags]
dec_act = dec_act[flags]

act_s18d_stack, stack_num = s18dStack(ra_act, dec_act, s18d, weight = True)

xspline = np.linspace(-10, 10, 40)
yspline = xspline

rbin = np.linspace(0, np.sqrt(2)*max(xspline), 20)


act_r_mean, act_bin_data, c, d = binmap(act_s18d_stack,xspline,yspline,rbin)
mdcw_r_mean, mdcw_bin_data, c, d = binmap(mdcw_s18d_stack,xspline,yspline,rbin)

#Bootstrap 

bootstraps = {'act_rs':[], 'act_bin_data':[],'mdcw_rs':[], 'mdcw_bin_data':[]}

for j in range(50):
    flags = np.random.randint(len(ra_act), size = len(ra_act))


    ra_temp_act = ra_act[flags]
    dec_temp_act = dec_act[flags]
    
    ra_temp_mdcw = ra_mdcw[flags]
    dec_temp_mdcw = dec_mdcw[flags]

    stack_act, stack_num_act = s18dStack(ra_temp_act, dec_temp_act, s18d, weight = True)
    act_r_jk, act_bin_data_jk, c, d = binmap(stack_act,xspline,yspline,rbin)
    
    stack_mdcw, stack_num_mdcw = s18dStack(ra_temp_mdcw, dec_temp_mdcw, s18d, weight = True)
    mdcw_r_jk, mdcw_bin_data_jk, c, d = binmap(stack_mdcw,xspline,yspline,rbin)

    bootstraps['act_rs'].append(act_r_jk)
    bootstraps['act_bin_data'].append(act_bin_data_jk)
    bootstraps['mdcw_rs'].append(mdcw_r_jk)
    bootstraps['mdcw_bin_data'].append(mdcw_bin_data_jk)
    print(j)

bootstraps['act_mean'] = act_bin_data
bootstraps['mdcw_mean'] = mdcw_bin_data

pk.dump(jks, open( "/scratch/r/rbond/jorlo/madcowsclusters/weighted_centers_bootstraps.p", "wb" ))




