fromm __future__ import print_function
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
	tempdec, tempra = np.deg2rad([decs[i], ras[i]])
	tempwid = np.deg2rad(width)
	box = [[tempdec-tempwid,tempra-tempwid],[tempdec+tempwid,tempra+tempwid]]
	maskstamp = mask1.submap(box)
	#Check if the stamp is entirely within borders, i.e. if the mask stamp is entirely ones
	if np.any(maskstamp[0]):
		for i in range(len(ras)):
			stamp = reproject.postage_stamp(map, ras[i], decs[i], width, 0.5)
			if stamp is None: continue
			stack += stamp[0]
			divisor += 1
	stack /= divisor
	return stack
