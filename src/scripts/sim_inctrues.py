from doppler_imaging import *
import numpy as np
import paths
import os

##############################################################################
####################    Configs     ##########################################
##############################################################################

from config_sim import *

target = "W1049A"
contrast = 0.7
roll = 0
noisetype = "random"
#goodchips_sim[instru][band] = [2, 3, 4]

#alpha = 100
#nlat, nlon = 20, 40
#modelmap = "SPOT"

modelmap = "testspots"

tobs = 7
inc_true = 70

maps = []
for inc in [50, 60, 70, 80, 90]:
    savedir = f"sim_inctrues/{inc}"
    if not os.path.exists(paths.figures / savedir):
        os.makedirs(paths.figures / savedir)

    #################### Automatic ####################################

    if True:
        cut = nk - 70
        # Auto consistent options
        if map_type == "eqarea":
            use_eqarea = True

        nobs = nobss[target]

        # set chips to include
        goodchips = goodchips_sim[instru][target][band]
        if use_toy_spec:
            goodchips = [4]
        nchip = len(goodchips)

        # set model files to use
        if "t1" in modelspec:
            if instru == "CRIRES": #TODO: put CRIRES data here
                model_datafiles = {"W1049B": 'fainterspectral-fits_6.pickle', "W1049A":'brighterspectral-fits_6.pickle'}
                model_datafile = model_datafiles[target]
            else:
                model_datafile = paths.data / f'{instru}_{target}_{band}_{modelspec}.pickle'
            pmod = f'linbroad_{modelspec}'
            rv = rvs[target]
            if use_toy_spec:
                print("Using toy spectrum...")
                pmod = f'toy_{modelspec}'
                rv = 0

        elif "lte" in modelspec: #TODO: shall I just delete this option
            if instru == "CRIRES":
                model_datafiles = {"W1049B": 'fainterspectral-fits_6.pickle', "W1049A":'brighterspectral-fits_6.pickle'}
                model_datafile = model_datafiles[target]
            else:
                model_datafile = paths.data / f'{instru}_{target}_{band}_{modelspec}.pickle'
            pmod = 'linbroad_lte015'
            rv = rvs[target]

        line_file = paths.data / f'linelists/{pmod}_edited.clineslsd'
        cont_file = paths.data / f'linelists/{pmod}C.fits'

        # set solver parameters
        inc = inc
        vsini = vsinis[target]
        veq = vsini / np.sin(inc * np.pi / 180)

        # set time and period parameters
        timestamp = np.linspace(0, tobs, nobs)  # observed time points in hours
        phases = timestamp * 2 * np.pi / period # 0 ~ 2*pi in rad     # guessed
        theta = 360.0 * timestamp / period      # 0 ~ 360 in degree   # guessed 

        assert nobs == len(theta)

        kwargs_sim = dict(
            ydeg=ydeg_sim,
            udeg=udeg,
            nc=nc,
            veq=veq,
            inc=inc_true,
            nt=nobs,
            vsini_max=vsini_max,
            u1=u1,
            theta=theta)

        kwargs_run = kwargs_sim.copy()
        kwargs_run['ydeg'] = ydeg
        kwargs_run['inc'] = inc

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


    ##############################################################################
    ####################      Run!      ##########################################
    ##############################################################################

    assert simulation_on == True

    # Load data from fit pickle
    mean_spectrum, template, observed, residual, error, wav_nm, wav0_nm = load_data(model_datafile, instru, nobs, goodchips)

    # Make mock observed spectra
    observed, fakemap = spectra_from_sim(modelmap, contrast, roll, smoothing, mean_spectrum, wav_nm, wav0_nm, error, residual, 
                            noisetype, kwargs_sim, savedir, plot_ts=False, plot_IC14=False, colorbar=False)
    # Compute LSD mean profile
    intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, vsini, rv, 
                                                         period, timestamp, savedir, cut=cut)

    bestparamgrid_r, res = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=False, colorbar=False) #, vmin=40, vmax=160)
    maps.append(bestparamgrid_r)

    #LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=False, colorbar=False)

    #LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, lr=lr_LSD, niter=niter_LSD, annotate=True)

    #lin_map = solve_starry_lin(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, annotate=True)
    #plt.figure(figsize=(5,3))
    #plt.savefig(paths.figures / f"{savedir}/solver4.pdf", bbox_inches="tight", dpi=300)

    #opt_map = solve_starry_opt(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, lr=lr, niter=niter, annotate=True)

