from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################

from config_run import *

savedir = "igrinsK"
band = "K"

use_eqarea = True
target = "W1049B"
#modelspec = "t1400g1000f8"

#################### Automatic ####################################

if True:
    # Auto consistent options
    contrast = "real"
    noisetype = "real"
    cut = nk - 70

    nobs = nobss[target]

    # set chips to include
    goodchips = goodchips_run[instru][target][band]
    nchip = len(goodchips)

    # set model files to use
    if "t1" in modelspec:
        model_datafile = paths.data / f'{instru}_{target}_{band}_{modelspec}.pickle'
        pmod = f'linbroad_{modelspec}'
        rv = rvs[target]

    line_file = paths.data / f'linelists/{pmod}_edited.clineslsd'
    cont_file = paths.data / f'linelists/{pmod}C.fits'

    # set solver parameters
    period = periods[target]
    inc = incs[target]
    vsini = vsinis[target]
    veq = vsini / np.sin(inc * np.pi / 180)

    # set time and period parameters
    timestamp = timestamps[target]
    phases = timestamp * 2 * np.pi / period # 0 ~ 2*pi in rad
    theta = 360.0 * timestamp / period      # 0 ~ 360 in degree

    kwargs_sim = dict(
        ydeg=ydeg_sim,
        udeg=udeg,
        nc=nc,
        veq=veq,
        inc=inc,
        nt=nobs,
        vsini_max=vsini_max,
        u1=u1,
        theta=theta)

    kwargs_run = kwargs_sim.copy()
    kwargs_run['ydeg'] = ydeg

    kwargs_IC14 = dict(
        phases=phases, 
        inc=inc, 
        vsini=vsini, 
        LLD=LLD, 
        eqarea=use_eqarea, 
        nlat=nlat, 
        nlon=nlon,
        alpha=alpha,
        ftol=ftol
    )

    kwargs_fig = dict(
        goodchips=goodchips,
        noisetype=noisetype,
        contrast=contrast,
        savedir=savedir
    )

assert simulation_on == False
print(f"Using real observation {model_datafile}")

# Load data from pickle fit
mean_spectrum, template, observed, residual, error, wav_nm, wav0_nm = load_data(model_datafile, instru, nobs, goodchips)

# Compute LSD mean profile
intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, vsini, rv, 
                                                     period, timestamps[target], savedir, cut=cut)

bestparamgrid_r, res = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, 
                                         annotate=False, colorbar=False, plot_cells=False, spotfit=False,
                                         create_obs_from_diff=True)


### DIME fit plot
sc_observation_1d = res['sc_observation_1d']
model = res['model_observation']
weight = res['w_observation']
plt.figure(figsize=(15,3))
plt.plot(np.arange(sc_observation_1d.size), sc_observation_1d, '.', color='k', markersize=4)
plt.plot(np.arange(sc_observation_1d.size), model, color='r')
#plt.plot(np.arange(sc_observation_1d.size), weight/weight.max()/15+sc_observation_1d.min(), color='gray')
plt.legend(loc=4,labels=['observation', 'max entropy fit'], bbox_to_anchor=(1,1))
plt.xticks(np.linspace(0+nk/2, sc_observation_1d.size-nk/2, 14), np.arange(1,15))
plt.xlabel("$n_{obs}$")
plt.savefig(paths.figures / "dime_fit.png", bbox_inches="tight", dpi=200, pad_inches=0)


### DIME cell  plot
def plot_map(map_obj):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.grid(False)
    good = (map_obj.projected_area>0)
    for k in range(map_obj.ncell):
        lats = map_obj.corners_latlon[k][0]
        lons = map_obj.corners_latlon[k][1]

        y = np.array([lats[0], lats[1], lats[3], lats[2]]) - np.pi/2
        x = np.array([lons[0], lons[1], lons[3], lons[2]]) - np.pi
        # Plot the polygon
        if good[k]:
            poly = plt.Polygon(np.column_stack((x, y)), facecolor='white', edgecolor='black', linewidth=0.5)
        else:
            poly = plt.Polygon(np.column_stack((x, y)), facecolor='darkgray', edgecolor='black', linewidth=0.5)
        ax.add_patch(poly)
        ax.text(x.mean()-0.1, y.mean()-0.1, f"{k}", size=7.5, alpha=0.7)
        #ax.text(x.mean()-0.1, y.mean()-0.07, f"a:{map_obj.projected_area[k]:.3f}", size=5)

    # Set plot parameters
    ax.set_xticklabels([f'{deg}˚' for deg in [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]], fontsize=8, alpha=0.5)
    ax.set_xticks([])
    ax.set_yticks(np.arange(-72, 73, 18)*np.pi/180, [f'{deg}˚' for deg in np.arange(-72, 73, 18)], fontsize=8, alpha=0.7)

deltaphi = np.pi
inc_plt = (90-inc) *  np.pi / 180
pltmap = ELL_map.map(nlat=10, nlon=18, type='eqarea', inc=inc_plt, deltaphi=deltaphi)
plot_map(pltmap)
plt.savefig(paths.figures / "cells.png", bbox_inches="tight", dpi=200, pad_inches=0)


### DIME Rmatrix plot
plt.figure(figsize=(12,12))
Rmatrix = res['Rmatrix']
Rmat = np.where(Rmatrix==0, np.nan, Rmatrix)
plt.imshow(Rmat, aspect=7, cmap="inferno")
plt.xlabel("$n_{obs} \cdot n_k$")
plt.ylabel("ncell")
plt.xticks(np.arange(0, nk*nobs+1, nk), [f"{n}$ n_k$" for n in np.arange(0, nobs+1)])
plt.colorbar(aspect=25, pad=0.02)
plt.subplots_adjust(bottom=0.5)
from matplotlib.patches import Ellipse
r=9
plt.gca().add_artist(Ellipse((812, 93), 7*r, r, facecolor="none", edgecolor="deepskyblue", linewidth=1))
plt.savefig(paths.figures / "Rmatrix.png", bbox_inches="tight", dpi=200, pad_inches=0)