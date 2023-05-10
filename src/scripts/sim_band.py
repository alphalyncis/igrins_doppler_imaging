from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################
from config_sim import *
modelmap = "1band"
savedir = "sim_band"

##############################################################################
####################      Run!      ##########################################
##############################################################################

assert simulation_on == True

# Load data from fit pickle
mean_spectrum, template, observed, residual, error, wav_nm, wav0_nm = load_data(model_datafile, instru, nobs, goodchips)

# Make mock observed spectra
observed = spectra_from_sim(modelmap, contrast, roll, smoothing, fakemap_nlat, fakemap_nlon, mean_spectrum, error, residual, noisetype, kwargs_sim, savedir, plot_ts=True)

# Compute LSD mean profile
intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, rv, period, savedir)

bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=True)

LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=True)

LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, lr=lr_LSD, niter=niter_LSD, annotate=True)

lin_map = solve_starry_lin(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, annotate=True)
#plt.figure(figsize=(5,3))
#plt.savefig(paths.figures / f"{savedir}/solver4.pdf", bbox_inches="tight", dpi=300)

opt_map = solve_starry_opt(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, lr=lr, niter=niter, annotate=True)

