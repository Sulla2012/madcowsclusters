from astropy.table import QTable
import astropy.units as u
import numpy as np

t = QTable.read('AdvACT_S18Clusters_v1.0-beta.fits')

madcows = t['MaDCoWS']

t_mark = QTable.read('AdvACT_mark.fits')
names_mark_temp = t_mark['name']
names_mark = np.array(names_mark_temp)

i = 0
names_official = []
for a in madcows:
    if a:
        i = i+1
        names_official.append(t['name'][i])

print('Official count: {}'.format(i))
print('Mark count: {}'.format(len(names_mark)))

print(names_official)
print(names_mark)

ra_mark = t_mark['RADeg']
i = 0
ra_official = []
for a in madcows:
    if a:
        ra_official.append(t['RADeg'][i])

print(ra_mark)
print(ra_official)

