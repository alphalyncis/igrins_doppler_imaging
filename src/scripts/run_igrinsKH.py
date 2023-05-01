from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################

if True:
    ####################   Constants    ##########################################

    goodchips_run = {
        "IGRINS":{
            "W1049B":{
                "K": [1, 4, 13], #[0, 1, 2, 3, 4, 5, 13, 14, 15, 16, 18], 
                "H": [1,2,3,4] #[0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19]
            },
            "W1049A":{
                "K": [0, 1, 2, 3, 4, 5, 13, 14, 15, 16, 17, 18], 
                "H": [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19]
            },
            "2M0036":{
                "K": [2,3,4,5,8,10,11,12,13,14,15,18], 
                "H": [1,2,3,4,5,9,10,11,16], #[0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19]
            }
        }
    }

    periods = {"W1049B": 5.28,   "W1049A": 7,      "2M0036": 2.7 }
    incs =    {"W1049B": 80,     "W1049A": 70,     "2M0036": 51  }
    vsinis =  {"W1049B": 29e3,   "W1049A": 21e3,   "2M0036": 32e3}
    rvs =     {"W1049B": 7.4e-5, "W1049A": 9.3e-5, "2M0036": 6e-5}
                    #9e-5 9.3e-5if CRIRES 7.4e-5 5.4e-5 if IGRINS

    timestamps = { # obs times in hour, computed from obs headers (JD-DATE)
        "W1049B": np.array(
            [0.134892  , 0.49418401, 0.85395001, 1.213902  , 1.5732    ,
            1.93294201, 2.30937601, 2.66913001, 3.19404001, 3.61374   ,
            3.987222  , 4.35165001, 4.71316801, 5.07270601]), 
        "W1049A": np.array(
            [0.134892  , 0.49418401, 0.85395001, 1.213902  , 1.5732    ,
            1.93294201, 2.30937601, 2.66913001, 3.19404001, 3.61374   ,
            3.987222  , 4.35165001, 4.71316801, 5.07270601]),
        "2M0036": np.array(
            [0.        , 0.35493599, 0.709224  , 1.063536  , 2.07012   ,
            2.424528  , 2.779728  ])
    }


    ##############################################################################
    ####################    Settings    ##########################################
    ##############################################################################

    simulation_on = False
    savedir = "igrinsHK"
    use_toy_spec = False

    #################### Run settings ####################################

    flux_err = 0.01 if use_toy_spec else 0.025
    instru = "IGRINS"
    target = "W1049B"
    band = "both"
    #solver = "IC14new"
    map_type = "eqarea"
    nobs = 14 if "W1049" in target else 7

    modelspec = "t1500g1000f8"
    LSD = "new"

    ########## IC14 parameters ##########
    nk = 103
    LLD = 1.0
    alpha = 4500
    ftol = 0.01 # tolerance for convergence of maximum-entropy
    nstep = 2000
    nlat, nlon = 10, 20

    ########## Starry parameters ##########
    ydeg_sim = 15
    ydeg = 8
    udeg = 1
    nc = 1
    u1 = 0.3
    vsini_max = 40000.0

    ########## Starry optimization parameters ##########
    lr_LSD = 0.001
    niter_LSD = 5000
    lr = 0.01
    niter = 5000

    #################### Automatic ####################################

    # Auto consistent options
    contrast = "real"
    noisetype = "real"
    if map_type == "eqarea":
        use_eqarea = True

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

    print(f"Using real observation {model_datafileK}")

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

nk = 71
##############################################################################
####################      Run!      ##########################################
##############################################################################

assert simulation_on == False
assert savedir == "igrinsHK"

# Load data from pickle fit
wav=[]
wav0=[]
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
intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, vsini, rv, period, savedir)

# Solve by 5 solvers
bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=True)

LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=True)

LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, lr=lr_LSD, niter=2000, annotate=True)

lin_map = solve_starry_lin(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, annotate=True)
#plt.figure(figsize=(5,3))
#plt.savefig(paths.figures / f"{savedir}/solver4.pdf", bbox_inches="tight", dpi=300)

opt_map = solve_starry_opt(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, lr=lr, niter=2000, annotate=True)

print("Run success.")

