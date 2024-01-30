from doppler_imaging import *
import numpy as np
import paths

##############################################################################
####################    Configs     ##########################################
##############################################################################

from config_run import *

target = "W1049B_0209"
savedir = f"igrinsHK_{target}"
band = "both"
modelspecK = "t1500g1000f8"
modelspecH = "t1500g1000f8"
if "W1049B" in target:
    modelspecH = "t1400g1000f8"

#################### Automatic ####################################
if True:
    if not os.path.exists(paths.figures / savedir):
        os.makedirs(paths.figures / savedir)

    cut = nk - 70

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
    if "t1" in modelspecK:
        model_datafileK = paths.data / f'{instru}_{target}_K_{modelspecK}.pickle'
        model_datafileH = paths.data / f'{instru}_{target}_H_{modelspecH}.pickle'
        pmod = f'linbroad_{modelspecK}'
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
intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, 
                                                    vsini, rv, period, timestamps[target], savedir, cut=cut)

    
# Solve by 5 solvers
bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, 
                                                annotate=False, colorbar=False, spotfit=False)
#plot_IC14_map(bestparamgrid_r, colorbar=False, vmin=85, vmax=110)
mapB_HK = bestparamgrid_r.copy()

print("Run success.")

