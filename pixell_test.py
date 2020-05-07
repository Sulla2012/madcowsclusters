from __future__ import print_function
from pixell import enmap,utils
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import urllib.request
from astropy.table import QTable
import astropy.units as u
from pixell import reproject
from pixell import enplot

def stack(ras, decs, map):
    stack = 0
    for i in range(len(ras)):
        stamp = reproject.postage_stamp(map, ras[i], decs[i], 20., 0.5)
        if stamp is None: continue
        stack += stamp[0]
    stack /= len(ras)
    return stack

print("Starting ACT Maps")

imap = enmap.read_map('/maps/ACTPol_148_D56_PA1_S2_1way_I.fits')

plots = enplot.plot(imap,range=300,mask=0)

enplot.write("ACT_map",plots)

width = np.deg2rad(30./60.)

t = QTable.read('madcows/AdvACT_S18Clusters_v1.0-beta.fits')
ra_temp = t['RADeg']
dec_temp = t['decDeg']
ra, dec = np.array(ra_temp), np.array(dec_temp)
s18stack = stack(ra, dec, imap)
plots = enplot.plot(enmap.upgrade(s18stack,5),grid=False, colorbar=True,color='gray')
enplot.write("s18stack_act",plots)

for i in range(len(ra)):
    #print(str(i)+" "+str(ra[i])+" "+str(dec[i]))
    
    if 37>ra[i]>0 and -10<dec[i]<6:
        tempdec, tempra = np.deg2rad([dec[i], ra[i]])
        box = [[tempdec-width/2.,tempra-width/2.],[tempdec+width/2.,tempra+width/2.]]
        stamp = imap.submap(box)

        plt.imshow(stamp)
        plt.colorbar()
        plt.savefig("madcows/stamps/act_s18_stamp"+str(i)+".png")
        plt.close()

        plots = enplot.plot(enmap.upgrade(stamp,5),grid=False, colorbar=True,color='gray')
        enplot.write("madcows/stamps/act_en_s18_stamp"+str(i),plots)


t_mc = QTable.read('madcows/AdvACT_4sig.fits')

ra, dec = t_mc['RADeg'], t_mc['decDeg']

mcstack = stack(ra, dec, imap)
plots = enplot.plot(enmap.upgrade(mcstack,5),grid=False, colorbar=True,color='gray')
enplot.write("act_mcstack",plots)

for i in range(len(ra)):
    #print(str(i)+" "+str(ra[i])+" "+str(dec[i]))
    
    if 37>ra[i]>0 and -10<dec[i]<6:
        tempdec, tempra = np.deg2rad([dec[i], ra[i]])
        box = [[tempdec-width/2.,tempra-width/2.],[tempdec+width/2.,tempra+width/2.]]
        stamp = imap.submap(box)

        plt.imshow(stamp)
        plt.colorbar()
        plt.savefig("madcows/stamps/act_mc_stamp"+str(i)+".png")
        plt.close()

        plots = enplot.plot(enmap.upgrade(stamp,5),grid=False, colorbar=True,color='gray')
        enplot.write("madcows/stamps/act_en_mc_stamp"+str(i),plots)



#########################################################################################################################################################################################

print("Starting Planck Maps")

imap = enmap.read_map('HFI_SkyMap_143_2048_R2.02_full.fits')

plots = enplot.plot(imap,range=300,mask=0)

enplot.write("planck_map",plots)

t = QTable.read('madcows/AdvACT_S18Clusters_v1.0-beta.fits')
ra_temp = t['RADeg']
dec_temp = t['decDeg']
ra, dec = np.array(ra_temp), np.array(dec_temp)
s18stack = stack(ra, dec, imap)
plots = enplot.plot(enmap.upgrade(s18stack,5),grid=False, colorbar=True,color='gray')
enplot.write("planck_s18stack",plots)



for i in range(len(ra)):
    #print(str(i)+" "+str(ra[i])+" "+str(dec[i]))

    if 37>ra[i]>0 and -10<dec[i]<6:
        tempdec, tempra = np.deg2rad([dec[i], ra[i]])
        box = [[tempdec-width/2.,tempra-width/2.],[tempdec+width/2.,tempra+width/2.]]
        stamp = imap.submap(box)

        plt.imshow(stamp)
        plt.colorbar()
        plt.savefig("madcows/stamps/planck_s18_stamp"+str(i)+".png")
        plt.close()

        plots = enplot.plot(enmap.upgrade(stamp,5),grid=False, colorbar=True,color='gray')
        enplot.write("madcows/stamps/planck_en_s18_stamp"+str(i),plots)


ra, dec = t_mc['RADeg'], t_mc['decDeg']

mcstack = stack(ra, dec, imap)
plots = enplot.plot(enmap.upgrade(mcstack,5),grid=False, colorbar=True,color='gray')
enplot.write("planck_mcstack",plots)


for i in range(len(ra)):
    #print(str(i)+" "+str(ra[i])+" "+str(dec[i]))

    if 37>ra[i]>0 and -10<dec[i]<6:
        tempdec, tempra = np.deg2rad([dec[i], ra[i]])
        box = [[tempdec-width/2.,tempra-width/2.],[tempdec+width/2.,tempra+width/2.]]
        stamp = imap.submap(box)

        plt.imshow(stamp)
        plt.colorbar()
        plt.savefig("madcows/stamps/planck_mc_stamp"+str(i)+".png")
        plt.close()

        plots = enplot.plot(enmap.upgrade(stamp,5),grid=False, colorbar=True,color='gray')
        enplot.write("madcows/stamps/planck_en_mc_stamp"+str(i),plots)

