import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.table import Table
import glob
import paths
c = 2.99792e5 # km/s

print("Running fit.py... Imports ok.")

nobs=14
chipmods = {}
chiplams = {}
chisq = {}
vsini = {}
rv = {}
lld = {}
wcoef = {}
modellist = []

for f in sorted(glob.glob(f"{paths.data.as_posix()}/IGRINS_W1049B_K_binned_chipmods_*.fits")):
    modelname = f.split("_")[-1][:12]
    chipmods[modelname] = fits.getdata(f)
    modellist.append(modelname)
for f in sorted(glob.glob(f"{paths.data.as_posix()}/IGRINS_W1049B_K_binned_chiplams_*.fits")):
    modelname = f.split("_")[-1][:12]
    chiplams[modelname] = fits.getdata(f)
for f in sorted(glob.glob(f"{paths.data.as_posix()}/IGRINS_W1049B_K_binned_*.txt")):
    modelname = f.split("_")[-1][:12]
    results = Table.read(f, format='ascii')
    chisq[modelname] = results['chisq']
    vsini[modelname] = results['vsini']
    rv[modelname] = results['rv']
    lld[modelname] = results['lld']
    wcoef[modelname] = results['wcoef']

# find best fitting model
df = pd.DataFrame({
    'model name':[model for model in modellist], 
    'median':[np.median(chisq[model]) for model in modellist], 
    'mean':[np.mean(chisq[model]) for model in modellist], 
    'min':[np.min(chisq[model]) for model in modellist]
})
df = df.sort_values(by=['median'])
print("table:\n", df)
open(paths.output / "fit_table.txt", "w").write(df.to_latex())


# plot best fitting spectrum
model = 't1500g1000f8'
pad = 50
filename = f'IGRINS_W1049B_K_{model}.pickle'
with open(paths.data / filename, 'rb') as f:
    ret = pickle.load(f, encoding="latin1")

fobs = ret['fobs0']
fobs = fobs[:20]
wobs = ret['wobs']
plt.figure(figsize=(15,10))

t=6
for sub in range(4):
    plt.subplot(4,1, 4-sub)
    for jj in range(sub*5, sub*5+5):
        plt.plot(chiplams[model][t,jj, pad:-pad], fobs[t,jj, pad:-pad], linewidth=0.5, color="black", label="observation")
        plt.plot(chiplams[model][t,jj, pad:-pad], chipmods[model][t,jj, pad:-pad], linewidth=0.5, color="r", label="fitted")
        plt.plot(chiplams[model][t,jj,pad:-pad], fobs[t,jj,pad:-pad] - chipmods[model][t,jj,pad:-pad], linewidth=0.5, color="gray")
    plt.ylabel("flux")
    #plt.ylim((0,2))
    if sub==0:
        plt.xlabel("wavelength (micron)")
    if sub==3:
        custom_lines = [Line2D([0], [0], color="black", lw=1),
                        Line2D([0], [0], color="r", lw=1),
                        Line2D([0], [0], color="g", lw=1),]
        plt.legend(custom_lines, ["observation", "best-fitting spectrum", "residual"], loc=2)
plt.title(f"fitted K band spectrum using model Callie {model}, t={t}")
plt.savefig(paths.figures / "fit.pdf", bbox_inches="tight", dpi=300)


# vsini and rv for bestfittin model
print(f'vsini: {np.median(vsini[model]):.2f} ± {np.std(vsini[model]):.2f} km/s, rv: {np.median(rv[model])*c:.2f} ± {np.std(rv[model])*c:.2f} km/s')
colors = [cm.jet_r(x) for x in np.linspace(0, 1, nobs+4)]
plt.figure(figsize=(6,5))
lam_points = np.median(np.median(chiplams[model], axis=0), axis=1) # plot one point per wl
vsini_points = vsini[model].reshape((nobs, 20))
for k in range(nobs):
    plt.plot(lam_points, vsini_points[k], "o", mfc="none", color=colors[k], label=f't={k}')
plt.plot(lam_points, np.median(vsini_points, axis=0), "o", label="median", color="blue")
plt.xlabel("wavelength (micron)")
plt.ylabel("vsini (km/s)")
plt.ylim((10,50))
plt.title("vsini at each observation at each order")
plt.legend(loc=2, bbox_to_anchor=(1,1))
plt.text(chiplams[model].min(), 11, f"vsini: {np.median(vsini[model]):.2f} ± {np.std(vsini[model]):.2f} km/s")
plt.savefig(paths.figures / "fit_vsini.pdf", bbox_inches="tight", dpi=300)
open(paths.output / "fit_vsini.txt", "w").write(f"vsini = {np.median(vsini[model]):.2f} ± {np.std(vsini[model]):.2f} km/s")

plt.figure(figsize=(6,5))
lam_points = np.median(np.median(chiplams[model], axis=0), axis=1) # plot one point per wl
rv_points = rv[model].reshape((nobs, 20)) * c
for k in range(nobs):
    plt.plot(lam_points, rv_points[k], "o", mfc="none", color=colors[k], label=f't={k}')
plt.plot(lam_points, np.median(rv_points, axis=0), "o", label="median", color="blue")
plt.xlabel("wavelength (micron)")
plt.ylabel("rv * c (km/s)")
plt.ylim((10,50))
plt.title("rv at each observation at each order")
plt.legend(loc=2, bbox_to_anchor=(1,1))
plt.text(chiplams[model].min(), 11, f"rv: {np.median(rv[model])*c:.2f} ± {np.std(rv[model])*c:.2f} km/s")
plt.savefig(paths.figures / "fit_rv.pdf", bbox_inches="tight", dpi=300)
open(paths.output / "fit_rv.txt", "w").write(f"rv = {np.median(rv[model])*c:.2f} ± {np.std(rv[model])*c:.2f} km/s")