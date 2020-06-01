from pixell import enmap,utils
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import urllib.request

def mask(imap, ra, dec, width = 0.5, apod_pix = 4):
	ra, dec = np.deg2rad([ra,dec])
	width = np.deg2rad(width)
	box = [[dec-width/2.,ra-width/2.],[dec+width/2.,ra+width/2.]]
	stamp = imap.submap(box)
	taper = enmap.apod(stamp*0+1,apod_pix)
	return taper
