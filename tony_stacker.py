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
        
        #box = np.array([[ra-width/2,dec-width/2],[ra+width/2,dec+width/2]]) * utils.degree
        #shape, wcs = enmap.geometry(pos=box,res=0.5 * utils.arcmin,proj='car')
        #print(shape)
        #print(stampgeo.wcs)
        #print(stamp.wcs)
        #stamp.wcs = wcs
        #print(stamp.wcs)
        #print(stamp[0].shape)
        #plt.imshow(stamp[0])
        #plt.show()
        #Return map
        stamp.wcs.wcs.crval = [ra,dec]
        #plot = enplot.plot(stamp,mask=0)
        #enplot.show(plot)
        enmap.write_map('./for_tony/{}.fits'.format(name), stamp)
    return stamp

with open('/scratch/r/rbond/jorlo/S18d_202006/selFn/tileDefinitions.yml') as f:
    
    s18d = yaml.load(f)
mdcw_catalog = fits.open('/home/s/sievers/sdicker/ACTnCOWs/MADCOWSUnion.fits')

test= s18dStamp(303.1103611788376,-56.83093699573238, s18d, 'ACT-CL_J2012.4-5649', write = True)


