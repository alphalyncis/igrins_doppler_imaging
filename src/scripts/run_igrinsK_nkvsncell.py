from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################

from config_run import *

savedir = "igrinsK_nkvsncell"
instru = "IGRINS"
band = "H"

#################### Automatic ####################################

if True:
    # Auto consistent options
    contrast = "real"
    noisetype = "real"
    if map_type == "eqarea":
        use_eqarea = True

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

    print(f"Using real observation {model_datafile}")

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


##############################################################################
####################      Run!      ##########################################
##############################################################################

assert simulation_on == False

# Load data from pickle fit
mean_spectrum, template, observed, residual, error, wav_nm, wav0_nm = load_data(model_datafile, instru, nobs, goodchips)

for nlat, nlon in zip([7, 9, 10], [14, 18, 20]):
    for nk in [71, 101]:
        cut = nk - 70
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
        print(f"Running test nk={nk} vs nlatlon={nlat} {nlon}...")
        
        # Compute LSD mean profile
        intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, 
                                                             vsini, rv, period, timestamp, savedir, cut=cut)

        # Solve by IC14
        bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=False, colorbar=False)

        LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=False, colorbar=False)

    print("Run success.")

