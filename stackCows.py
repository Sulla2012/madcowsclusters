"""

Bin MaDCoWS forced photometry by richness, feed catalog through nemoMass, plot

"""

import os, sys
import numpy as np
import astropy.table as atpy
import pylab as plt
from nemo import catalogs
from nemo import plotSettings
import IPython

# Rescaling factor - multiply y by this much to crudely account for miscentering (Jack OS post)
rescaleFactor=1.25
#rescaleFactor=1.0

# Actual source catalog
richTab=atpy.Table().read("MADCOWSUnion.fits")
richTab=richTab[richTab['Rich'].mask == False]
richTab=richTab[richTab['Rich'] > 20]

# Get forced photometry catalog
inFileName="MADCOWSUnion_M500.fits"
if os.path.exists(inFileName) == False:
    os.system("wget https://astro.ukzn.ac.za/~mjh/ACT/releases/MADCOWSUnion_M500.fits --user=act --password=atacamallama")
forcedTab=atpy.Table().read(inFileName)

# Cross match - should have exactly same coords anyway as forced photometry
richTab, forcedTab, sep=catalogs.crossMatch(richTab, forcedTab, radiusArcmin = 0.1)
richTab['fixed_y_c']=forcedTab['fixed_y_c']
richTab['fixed_err_y_c']=forcedTab['fixed_err_y_c']

# Binning
numBins=8
numBootstraps=1000
rich=np.array(richTab['Rich'].data)
log10Rich=np.log10(rich)
binEdges=np.linspace(log10Rich.min(), log10Rich.max(), numBins+1)
binCentres=(binEdges[1:]+binEdges[:-1])/2.
mean_y=np.zeros(len(binCentres))
err_y=np.zeros(len(binCentres))
num_in_bin=np.zeros(len(binCentres))
for i in range(len(binCentres)):
    binMin=binEdges[i]
    binMax=binEdges[i+1]
    mask=np.logical_and(log10Rich >= binMin, log10Rich < binMax)
    yvals=richTab['fixed_y_c'][mask]*rescaleFactor  # Optional rescaling to account for miscentering
    mean_y[i]=np.mean(yvals)
    num_in_bin[i]=mask.sum()
    bsMeans=[]
    for j in range(numBootstraps):
        bsIndices=np.random.randint(0, len(yvals), len(yvals))
        bsMeans.append(np.mean(yvals[bsIndices]))
    err_y[i]=np.std(bsMeans)
    
binTab=atpy.Table()
binTab['name']=np.arange(0, numBins)+1
binTab['log10Rich']=binCentres
binTab['num_in_bin']=num_in_bin
binTab['redshift']=np.mean(richTab['Photz'])
binTab['redshiftErr']=0
binTab['fixed_y_c']=mean_y
binTab['fixed_err_y_c']=err_y
binTab['tileName']='1_10_8'     # Assumed representative! Just for Qs
binTab.write("Binned_Forced_MADCOWSUnion.fits", overwrite = True)

# Calc masses
os.system("nemoMass S18d_202006.yml -c Binned_Forced_MADCOWSUnion.fits -o Binned_Forced_MADCOWSUnion_M500.fits")

# Plot
tab=atpy.Table().read("Binned_Forced_MADCOWSUnion_M500.fits")
massCol='M500'
massSymbol="$M^{\\rm UPP}_{\\rm 500c}$"
outFileName="plot_binned_lambda_%s.png" % (massCol)

plotSettings.update_rcParams()
plt.figure(figsize=(9,6.5))
ax=plt.axes([0.14, 0.11, 0.85, 0.87])

plt.errorbar(np.power(10, tab['log10Rich']), tab[massCol], 
             yerr = [tab['%s_errMinus' % (massCol)], tab['%s_errPlus' % (massCol)]],
             fmt = 'D')

# Gonzalez et al. 2019 scaling relation
alpha=1.65  # +1.45, -0.96
beta=-2.16  # +1.57, -2.38
plotRange=np.linspace(0, 200, 100)
plotM=np.power(10, alpha*np.log10(plotRange)+beta)
plt.plot(plotRange, plotM, 'k--', label = 'Gonzalez et al. (2019)')

plt.xlabel("$\lambda_{15}$")
plt.ylabel("%s (10$^{14}$ $M_{\odot}$)" % (massSymbol))
plt.xlim(10, 110)
plt.legend()
plt.loglog()
plt.savefig(outFileName)
plt.savefig(outFileName.replace(".png", ".pdf"))
plt.close()

# Mean mass for lambda_15 > 20
meanM500=(tab['M500']*tab['num_in_bin']).sum()/tab['num_in_bin'].sum()
meanM500Cal=(tab['M500Cal']*tab['num_in_bin']).sum()/tab['num_in_bin'].sum()
print("rescaleFactor = %.2f applied" % (rescaleFactor))
print("lambda_15 > 20 cut:")
print("    mean M500c UPP = %.3f (10^14 Msun)" % (meanM500))
print("    mean M500c Cal = %.3f (10^14 Msun)" % (meanM500Cal))

