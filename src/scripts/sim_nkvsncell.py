from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################

from config_sim import *

savedir = "sim_nkvsncell"
contrast = 0.15
noisetype = "res+random"

#################### Automatic ####################################

if True:
    # Auto consistent options
    if map_type == "eqarea":
        use_eqarea = True

    nobs = nobss[target]

    # set chips to include
    goodchips = goodchips_sim[instru][band]
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
    period = periods[target]
    inc = incs[target]
    vsini = vsinis[target]
    veq = vsini / np.sin(inc * np.pi / 180)

    # set time and period parameters
    timestamp = np.linspace(0, period, nobs)  # simulate equal time interval obs
    phases = timestamp * 2 * np.pi / period # 0 ~ 2*pi in rad
    theta = 360.0 * timestamp / period      # 0 ~ 360 in degree

    assert nobs == len(theta)

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


##############################################################################
####################      Run!      ##########################################
##############################################################################

assert simulation_on == True
assert savedir == "sim_nkvsncell"

# Load data from fit pickle
mean_spectrum, template, observed, residual, error, wav_nm, wav0_nm = load_data(model_datafile, instru, nobs, goodchips)

# Make mock observed spectra
observed = spectra_from_sim(modelmap, contrast, roll, smoothing, fakemap_nlat, fakemap_nlon, mean_spectrum, wav_nm, wav0_nm, error, residual, noisetype, kwargs_sim, savedir, plot_ts=False)


for nlat, nlon in zip([7, 10, 16, 22], [14, 20, 32, 44]):
    for nk in [31, 51, 101, 151]:
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
        
        # Compute LSD mean profile
        intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, vsini, rv, period, savedir)

        bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=True)

        LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=True)

#LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, lr=lr_LSD, niter=niter_LSD, annotate=True)

#lin_map = solve_starry_lin(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, annotate=True)
#plt.figure(figsize=(5,3))
#plt.savefig(paths.figures / f"{savedir}/solver4.pdf", bbox_inches="tight", dpi=300)

#opt_map = solve_starry_opt(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, lr=lr, niter=niter, annotate=True)

