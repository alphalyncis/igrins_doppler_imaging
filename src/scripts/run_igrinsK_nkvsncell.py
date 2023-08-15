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
goodchips_run[instru][target][band] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
modelspec = "t1400g1000f8"

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

nks = [53,75,101,125,151,175,201]
std_LPs = []
snr_LPs = []
for nk in nks:
    cut = nk - 50
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
    intrinsic_profiles, obskerns_norm = make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, 
                                                            vsini, rv, period, timestamp, savedir, cut=cut)
    smoothed = savgol_filter(obskerns_norm, 31, 3)
    resid = obskerns_norm - smoothed
    err_pix = np.array([np.abs(resid[:,:,pix] - np.median(resid, axis=2)) for pix in range(nk)]) # error of each pixel in LP by MAD
    err_LP = 1.4826 * np.median(err_pix, axis=0) # error of each LP (at each t and chip)

    signal = 1 - smoothed.min(axis=2).mean(axis=0) # signal = line depth
    #noise = np.r_[obskerns_norm[:,:,:int(cut/2+1)], obskerns_norm[:,:,-int(cut/2+1):]].std(axis=2).mean(axis=0) # noise = std of the flat part
    noise = err_LP.mean(axis=0) # mean error of a chip
    snr_LPs.append(signal/noise) # S/N of LP of each chip
    std_LPs.append(noise)

    # Solve by IC14
    #bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, annotate=False, colorbar=False)

    #LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=False, colorbar=False)

plt.figure(figsize=(5,4))
snr_LPs = np.array(snr_LPs)
for i in range(10):  
    plt.plot(nks, snr_LPs[:, i], marker=".")
plt.xlabel("nk")
plt.ylabel("S/N of line profile")
plt.legend(labels=[f"chip{c}" for c in goodchips[:10]], fontsize=8, bbox_to_anchor=(1,1))

plt.figure(figsize=(5,4))
snr_LPs = np.array(snr_LPs)
for i in range(10,20):  
    plt.plot(nks, snr_LPs[:, i], marker=".")
plt.xlabel("nk")
plt.ylabel("S/N of line profile")
plt.legend(labels=[f"chip{c}" for c in goodchips[10:]], fontsize=8, bbox_to_anchor=(1,1))


plt.figure(figsize=(5,4))
std_LPs = np.array(std_LPs)
for i in range(10):  
    plt.plot(nks, std_LPs[:, i], marker=".")
plt.xlabel("nk")
plt.ylabel("noise level of line profile")
plt.legend(labels=[f"chip{c}" for c in goodchips[:10]], fontsize=8, bbox_to_anchor=(1,1))

plt.figure(figsize=(5,4))
std_LPs = np.array(std_LPs)
for i in range(10,20):  
    plt.plot(nks, std_LPs[:, i], marker=".")
plt.xlabel("nk")
plt.ylabel("noise level of line profile")
plt.legend(labels=[f"chip{c}" for c in goodchips[10:]], fontsize=8, bbox_to_anchor=(1,1))


print("Run success.")

