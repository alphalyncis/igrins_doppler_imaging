from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################

from config_run import *

savedir = "igrinsHK"
band = "both"
nk = 71
nlat, nlon = 10, 20


#################### Automatic ####################################
if True:
    # Auto consistent options
    contrast = "real"
    noisetype = "real"
    if map_type == "eqarea":
        use_eqarea = True
    
    nobs = nobss[target]

    # set chips to include
    goodchipsK = goodchips_run[instru][target]["K"]
    goodchipsH = goodchips_run[instru][target]["H"]

    # set model files to use
    if "t1" in modelspec:
        model_datafileK = paths.data / f'{instru}_{target}_K_{modelspec}.pickle'
        model_datafileH = paths.data / f'{instru}_{target}_H_{modelspec}.pickle'
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

##############################################################################
####################      Run!      ##########################################
##############################################################################

assert simulation_on == False
assert savedir == "igrinsHK"
print(f"Using real observation {model_datafileK} + {model_datafileH}")

# Load data from pickle fit
print("K band:")
mean_spectrumK, templateK, observedK, residualK, errorK, wav_nmK, wav0_nmK = load_data(model_datafileK, instru, nobs, goodchipsK)
#print("mean_spetrumK:", mean_spectrumK.shape)
#print("observedK:", observedK.shape)

print("\nH band:")
mean_spectrumH, templateH, observedH, residualH, errorH, wav_nmH, wav0_nmH = load_data(model_datafileH, instru, nobs, goodchipsH)
#print("mean_spetrumH:", mean_spectrumH.shape)
#print("observedH:", observedH.shape)

wav_nm = np.concatenate((wav_nmH, wav_nmK), axis=0)
wav0_nm = np.concatenate((wav0_nmH, wav0_nmK), axis=0)
mean_spectrum = np.concatenate((mean_spectrumH, mean_spectrumK), axis=0)
template = np.concatenate((templateH, templateK), axis=1)
observed = np.concatenate((observedH, observedK), axis=1)
residual = np.concatenate((residualH, residualK), axis=1)
error = np.concatenate((errorH, errorK), axis=1)
goodchips = goodchipsH + goodchipsK

print("\nAfter stacking:")
print("mean_spetrum:", mean_spectrum.shape)
print("observed:", observed.shape)
print("wav:", wav_nm.shape)
print("goodchips:", goodchips)

kwargs_fig = dict(
    goodchips=goodchips,
    noisetype=noisetype,
    contrast=contrast,
    savedir=savedir
)

# Compute LSD mean profile
intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, vsini, rv, period, timestamps[target], savedir, cut=1)

# Solve by 5 solvers
#bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=False, colorbar=False)

#LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=False, colorbar=False)

#LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, lr=lr_LSD, niter=2000, annotate=False, colorbar=False)

#lin_map = solve_starry_lin(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, annotate=True)
#plt.figure(figsize=(5,3))
#plt.savefig(paths.figures / f"{savedir}/solver4.pdf", bbox_inches="tight", dpi=300)

#opt_map = solve_starry_opt(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, lr=lr, niter=2000, annotate=True)

print("Run success.")

