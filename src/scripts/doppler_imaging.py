# Functions for solving doppler maps using different pipelines.
# Xueqing Chen 24.04.2023
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import starry
import pickle
import os
from astropy.io import fits
import scipy.constants as const
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.ndimage.filters import gaussian_filter
import paths

import pymc3 as pm
import pymc3_ext as pmx
import theano.tensor as tt
from tqdm import tqdm
from matplotlib.colors import Normalize
import emcee

import analysis3 as an # for an.lsq and an.gfit
import ELL_map_class as ELL_map
import dime3 as dime # Doppler Imaging & Maximum Entropy, needed for various funcs
import cartopy.crs as ccrs
from PIL import Image

#TODO: target should name as -> target+night, since can have several nights for one target 
#TODO: run bands separately or together
#TODO: test more parameters for starry solver
#TODO: test sampling rate

################################################################################
####################   Methods for workflow    #################################
################################################################################

def load_data(model_datafile, instru, nobs, goodchips, use_toy_spec=False):
    global npix, npix0, flux_err
    # Load model and data
    with open(model_datafile, "rb") as f:
        data = pickle.load(f, encoding="latin1")
    lams = data["chiplams"][0] # in um
    nchip = len(goodchips)
    npix = lams.shape[1]
    print(f"nobs: {nobs}, nchip: {nchip}, npix: {npix}")

    observed = np.empty((nobs, nchip, npix))
    template = np.empty((nobs, nchip, npix))
    residual = np.empty((nobs, nchip, npix))
    error = np.empty((nobs, nchip, npix))

    if instru == "IGRINS":
        for k in range(nobs):
            for i, jj in enumerate(goodchips):
                observed[k, i] = np.interp(
                    lams[jj], 
                    data["chiplams"][k][jj],
                    data["fobs0"][k][jj] #/ data["chipcors"][k][c+firstchip],
                )
                template[k, i] = np.interp(
                    lams[jj],
                    data["chiplams"][k][jj],
                    data["chipmodnobroad"][k][jj] #/ data["chipcors"][k][c+firstchip]
                )
                residual[k, i] = np.interp(
                    lams[jj], 
                    data["chiplams"][k][jj],
                    data["fobs0"][k][jj] - data["chipmods"][k][jj]
                )
                error[k, i] = np.interp(
                    lams[jj],
                    data["chiplams"][k][jj],
                    remove_spike(data["eobs0"][k][jj])
                )

    elif instru == "CRIRES":
        for k in range(nobs):
            for i, jj in enumerate(goodchips):
                observed[k, i] = np.interp(
                    lams[jj],
                    data["chiplams"][k][jj],
                    data["obs1"][k][jj] / data["chipcors"][k][jj],
                )
                template[k, i] = np.interp(
                    lams[jj],
                    data["chiplams"][k][jj],
                    data["chipmodnobroad"][k][jj] / data["chipcors"][k][jj],
                )
                residual[k, i] = np.interp(
                    lams[jj], 
                    data["chiplams"][k][jj],
                    data["obs1"][k][jj] - data["chipmods"][k][jj]
                )

    pad = 100
    npix0 = len(lams[0])
    npix = len(lams[0][pad:-pad])
    wav_nm = np.zeros((nchip, npix))
    wav0_nm = np.zeros((nchip, npix0))
    mean_spectrum = np.zeros((nchip, npix0))
    observed = observed[:, :, pad:-pad]
    flux_err = eval(f'{error.mean():.3f}') if instru =="IGRINS" else 0.02

    for i, jj in enumerate(goodchips):
        wav_nm[i] = lams[jj][pad:-pad] * 1000 # um to nm
        wav0_nm[i] = lams[jj] * 1000 # um to nm
        if use_toy_spec:
            toy_spec = (
                1.0
                - 0.50 * np.exp(-0.5 * (wav0_nm[i] - 2338) ** 2 / 0.03 ** 2)
                - 0.60 * np.exp(-0.5 * (wav0_nm[i] - 2345) ** 2 / 0.03 ** 2)
                - 0.20 * np.exp(-0.5 * (wav0_nm[i] - 2347) ** 2 / 0.03 ** 2)
            )
            for k in range(nobs):
                template[k, i] = toy_spec
            mean_spectrum[i] = toy_spec
        mean_spectrum[i] = np.mean(template[:, i], axis=0)

    print("mean_spectrum:", mean_spectrum.shape)
    print("template:", template.shape)
    print("observed:", observed.shape)
    print(f"wav: {wav_nm.shape}, wav0: {wav0_nm.shape}")

    return mean_spectrum, template, observed, residual, error, wav_nm, wav0_nm

def spectra_from_sim(modelmap, contrast, roll, smoothing, n_lat, n_lon, mean_spectrum, wav_nm, wav0_nm,
                     error, residual, noisetype, kwargs_sim, savedir, r=33, lat=30, pad=100, plot_ts=False, colorbar=True):
    nobs = error.shape[0]
    # create fakemap
    if modelmap == "1spot":
        spot_brightness = contrast
        print(f"Running spot brightness {spot_brightness*100}% of surrounding")
        fakemap = np.ones((n_lat, n_lon))
        x, y = np.meshgrid(np.linspace(0, 360, n_lon), np.linspace(-90, 90, n_lat), )
        fakemap[np.sqrt((y+lat)**2 + (x-n_lon/2)**2) <= r] = spot_brightness

    if modelmap == "2spot":
        spot_brightness = contrast
        lat2 = 45
        lon2 = 0
        r2 = 20
        print(f"Running spot brightness {spot_brightness*100}% of surrounding")
        fakemap = np.ones((n_lat, n_lon))
        x, y = np.meshgrid(np.linspace(0, 360, n_lon), np.linspace(-90, 90, n_lat), )
        fakemap[np.sqrt((y+lat)**2 + (x-n_lon/2)**2) <= r] = spot_brightness
        fakemap[np.sqrt((y+lat2)**2 + (x-n_lon+lon2)**2) <= r2] = spot_brightness

    elif modelmap == "1band":
        band_width = 15
        band_lat = 30
        amp = 1 - contrast
        print(f"Running wave max amplitude diff {amp*100}%")
        phase = 0.7 #0-1?
        fakemap = np.ones((n_lat, n_lon))
        x, y = np.meshgrid(np.linspace(0, 360, n_lon), np.linspace(-90, 90, n_lat), )
        band_ind = np.s_[(90-band_lat)-band_width:(90-band_lat)+band_width]
        fakemap[band_ind] += amp * np.sin((x[band_ind]/360 - phase) * 2*np.pi)
        #fakemap[band_ind] -= amp

    elif modelmap == "1uniband":
        print(f"Running band with brightness {contrast*100}% of surrounding")
        phase = 0.7 #0-1?
        fakemap = np.ones((n_lat, n_lon))
        x, y = np.meshgrid(np.linspace(0, 360, n_lon), np.linspace(-90, 90, n_lat), )
        band_ind = np.s_[(90-lat)-r:(90-lat)+r]
        fakemap[band_ind] = contrast
        
    elif modelmap == "2band":
        band_width = 10
        band_lat = 45
        band2_lat = 0
        amp = 1 - contrast
        print(f"Running 2 waves max amplitude diff {amp*100}%")
        phase = 0.55 
        phase2 = 0.75
        fakemap = np.ones((n_lat, n_lon))
        x, y = np.meshgrid(np.linspace(0, 360, n_lon), np.linspace(-90, 90, n_lat), )
        band_ind = np.s_[(90-band_lat)-band_width:(90-band_lat)+band_width]
        fakemap[band_ind] += amp * np.sin((x[band_ind]/360 - phase) * 2*np.pi)
        band2_ind = np.s_[(90-band2_lat)-band_width:(90-band2_lat)+band_width]
        fakemap[band2_ind] += amp * np.sin((x[band2_ind]/360 - phase2) * 2*np.pi)

    elif modelmap == "blank":
        fakemap = np.ones((n_lat, n_lon))

    elif modelmap == "gcm":
        fakemap = np.loadtxt(paths.data / 'modelmaps/gcm.txt')
        fakemap /= np.median(fakemap)
        diff = 1 - fakemap
        print(f"Running GCM peak amplitude diff {diff*100}%")
        ampold = diff.max()
        amp = 1 - contrast
        diffnew = diff * amp / ampold
        fakemap = 1 - diffnew
        n_lat, n_lon = fakemap.shape

    fakemap = np.roll(fakemap[::-1, :], shift=int(roll*n_lon), axis=1)

    # Compute simulated flux
    allchips_flux = []
    for i in range(wav_nm.shape[0]):
        sim_map = starry.DopplerMap(lazy=False, wav=wav_nm[i], wav0=wav0_nm[i], **kwargs_sim)
        sim_map.load(maps=[fakemap], smoothing=smoothing)
        sim_map[1] = kwargs_sim["u1"]

        flux_err_add = 0.02
        noise = {
            "none": np.zeros((nobs, npix)),
            "random": np.random.normal(np.zeros((nobs, npix)), flux_err),
            "obserr": error[:, i, pad:-pad],
            "residual": residual[:, i, pad:-pad],
            "res+random": residual[:, i, pad:-pad] + np.random.normal(np.zeros((nobs, npix)), flux_err_add)
        }

        sim_map.spectrum = mean_spectrum[i]
        model_flux = sim_map.flux(kwargs_sim["theta"])
        simulated_flux = model_flux + noise[noisetype]

        allchips_flux.append(simulated_flux)

    allchips_flux = np.array(allchips_flux)

    fig, ax = plt.subplots()
    sim_map.show(ax=ax, projection="moll", colorbar=colorbar)
    plt.savefig(paths.figures / f"{savedir}/fakemap.png", bbox_inches="tight", dpi=100, transparent=True)

    if plot_ts:
        plot_timeseries(sim_map, model_flux, kwargs_sim["theta"], obsflux=simulated_flux, overlap=2)

    observed = np.transpose(allchips_flux, axes=(1,0,2))

    return observed, fakemap

def make_LSD_profile(instru, template, observed, wav_nm, goodchips, pmod, line_file, cont_file, nk, vsini, rv, period, timestamps, savedir, pad=100, cut=30):
    global wav_angs, err_LSD_profiles, dbeta
    print(instru)
    nobs = observed.shape[0]
    nchip = len(goodchips)
    # Read daospec linelist
    lineloc, lineew, _ = dao_getlines(line_file)
    pspec_cont = fits.getdata(cont_file)
    hdr_pspec_cont = fits.getheader(cont_file)
    wspec = hdr_pspec_cont['crval1'] + np.arange(pspec_cont.size)*hdr_pspec_cont['cdelt1']
    factor = 1e11 if "t1" in pmod else 1e5 # don't know why different model needs scaling with a factor
    pspec_cont = pspec_cont/factor
    spline = interpolate.UnivariateSpline(wspec, pspec_cont, s=0.0, k=1.0) #set up interpolation over the continuum measurement
    #plt.figure(figsize=(6,1))
    #plt.plot(wspec, pspec_cont)
    #plt.title("daospec continuum")
    #plt.show()

    wav_angs = np.array(wav_nm) * 10 #convert nm to angstroms

    # Compute LSD velocity grid:
    dbeta = np.diff(wav_angs).mean()/wav_angs.mean()
    print("dbeta", dbeta)
    dx = - dbeta * np.arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5))
    dv = const.c*dx / 1e3 # km/s

    # Compute LSD:
    kerns = np.zeros((nobs, nchip, nk), dtype=float)
    modkerns = np.zeros((nobs, nchip, nk), dtype=float)
    deltaspecs = np.zeros((nobs, nchip, npix), dtype=float)
    for i, jj in enumerate(goodchips): 
        print("chip", jj)
        for kk in range(nobs):
            shift = 1. + rv  # best rv shift for Callie is 9e-5
            deltaspec = make_deltaspec(lineloc*shift, lineew, wav_angs[i], verbose=False, cont=spline(wav_angs[i]))
            m,kerns[kk,i],b,c = dsa(deltaspec, observed[kk,i], nk)
            m,modkerns[kk,i],b,c = dsa(deltaspec, template[kk,i,pad:-pad], nk) 
            deltaspecs[kk,i] = deltaspec

    # Plot lines vs. model
    plt.figure(figsize=(15, 2*nchip))
    t=0
    for i, jj in enumerate(goodchips):
        plt.subplot(nchip,1,i+1)
        plt.plot(wav_angs[i], deltaspecs[t,i], linewidth=0.5, color='C0', label="deltaspec")
        plt.plot(wav_angs[i], template[t,i,pad:-pad], linewidth=0.6, color='C1', label="chipmodnobroad")
        plt.plot(wav_angs[i], observed[t,i], linewidth=0.6, color='k', label="obs")
        plt.text(x=wav_angs[i].min()-10, y=0, s=f"order={jj}")
        if i==0:
            plt.title(f"{pmod} model vs. lines at t={t}")
    plt.legend(loc=4, fontsize=9)
    #plt.show()
    plt.savefig(paths.output / "LSD_deltaspecs.png", transparent=True)
    
    # shift kerns to center
    modkerns, kerns = shift_kerns_to_center(modkerns, kerns, instru, goodchips, dv)
    
    # plot kerns
    #plot_kerns_timeseries(kerns, goodchips, dv, gap=0.02)
    #plot_kerns_timeseries(modkerns, goodchips, dv, gap=0.1)

    err_LSD_profiles = np.median(kerns.mean(1).std(0)) 
    # the error level across different obs of the chip-avged profile, median over nk pixels


    # normalize kerns
    obskerns_norm = cont_normalize_kerns(kerns, instru)
    
    # plot kerns + intrinsic_profile
    plot_kerns_timeseries(modkerns, goodchips, dv, gap=0.1)
    intrinsic_profiles = np.array([modkerns[:,i].mean(0) for i in range(nchip)])
    plot_kerns_timeseries(obskerns_norm, goodchips, dv, gap=0.02, normed=True, intrinsic_profiles=intrinsic_profiles)
    
    ### Plot averaged line shapes
    plot_chipav_kern_timeseries(obskerns_norm, dv, timestamps, savedir, gap=0.02, cut=int(cut/2+1))

    ### Plot deviation map for each chip and mean deviation map
    plot_deviation_map(obskerns_norm, goodchips, dv, vsini, timestamps, savedir, meanby="median", cut=cut)

    return intrinsic_profiles, obskerns_norm

def solve_IC14new(intrinsic_profiles, obskerns_norm, kwargs_IC14, kwargs_fig, clevel=7,
                  ret_both=True, annotate=False, colorbar=False, plot_cells=False, plot_starry=False, 
                  spotfit=False, create_obs_from_diff=True):
    print("*** Using solver IC14new ***")
    nobs, nk = obskerns_norm.shape[0], obskerns_norm.shape[2]

    bestparamgrid, res = solve_DIME(
        obskerns_norm, intrinsic_profiles,
        dbeta, nk, nobs, **kwargs_IC14, plot_cells=plot_cells, spotfit=spotfit,
        create_obs_from_diff=create_obs_from_diff
    )

    bestparamgrid_r = np.roll(
        np.flip(bestparamgrid, axis=1), int(0.5*bestparamgrid.shape[1]), axis=1)
    # TODO: derotate map??? seems like Ic14 maps are flipped and rolled 180 deg

    if plot_starry:
        fig, ax = plt.subplots(figsize=(7,3))
        showmap = starry.Map(ydeg=7)
        showmap.load(bestparamgrid_r)
        showmap.show(ax=ax, projection="moll", colorbar=colorbar)
    
    else:
        plot_IC14_map(bestparamgrid_r, clevel=clevel, sigma=2.) # smoothed contour lines

    map_type = "eqarea" if kwargs_IC14['eqarea'] else "latlon"
    if annotate:
        plt.text(-3.5, -1, f"""
            chip=averaged{kwargs_fig['goodchips']} 
            solver=IC14new {map_type} 
            noise={kwargs_fig['noisetype']} 
            err_level={flux_err} 
            contrast={kwargs_fig['contrast']} 
            limbdark={kwargs_IC14['LLD']}""",
        fontsize=8)
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver1.png", bbox_inches="tight", dpi=50, transparent=True)

    # Plot fit result
    obs_2d = np.reshape(res['sc_observation_1d'], (nobs, nk))
    bestmodel_2d = np.reshape(res['model_observation'], (nobs, nk))
    flatmodel_2d = np.reshape(res['flatmodel'], (nobs, nk))

    plt.figure(figsize=(5, 7))
    for i in range(nobs):
        plt.plot(res['dv'], obs_2d[i] - 0.02*i, color='k', linewidth=1)
        #plt.plot(obs[i] - 0.02*i, '.', color='k', markersize=2)
        plt.plot(res['dv'], bestmodel_2d[i] - 0.02*i, color='r', linewidth=1)
        plt.plot(res['dv'], flatmodel_2d[i] - 0.02*i, '--', color='gray', linewidth=1)
    plt.legend(labels=['obs', 'best-fit map', 'flat map'])
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver1_ts.png", bbox_inches="tight", dpi=100, transparent=True)

    if ret_both:
        return bestparamgrid_r, res
    else:
        return bestparamgrid_r

def solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, annotate=False, colorbar=True):
    print("*** Using solver LSD+starry_lin ***")
        
    mean_profile = np.median(intrinsic_profiles, axis=0) # mean over chips
    observation_norm = np.median(obskerns_norm, axis=1) # mean over chips

    # preprae data
    pad = 100
    model_profile = 1. - np.concatenate((np.zeros(pad), mean_profile, np.zeros(pad)))
    dlam = np.diff(wav_angs[0]).mean() / 10 # angstrom to nm
    lam_ref = wav_angs[0].mean() / 10 # angstrom to nm
    nw = len(model_profile)
    wav0_lsd = np.linspace(start=lam_ref-0.5*dlam*nw, stop=lam_ref+0.5*dlam*nw, num=nw)
    wav_lsd = wav0_lsd[pad:-pad]
    
    map_av = starry.DopplerMap(lazy=False, wav=wav_lsd, wav0=wav0_lsd, **kwargs_run)
    map_av.spectrum = model_profile
    map_av[1] = kwargs_run['u1']

    soln = map_av.solve(
        flux=np.flip(observation_norm, axis=1),
        theta=kwargs_run['theta'],
        normalized=True,
        fix_spectrum=True,
        flux_err=flux_err,
        quiet=os.getenv("CI", "false") == "true",
    )

    image = map_av.render(projection="moll")         
    fig, ax = plt.subplots(figsize=(7,3))
    map_av.show(ax=ax, projection="moll", image=image, colorbar=colorbar)
    if annotate:
        ax.annotate(f"""
            chip={kwargs_fig['goodchips']}
            solver=LSD+starry_lin
            noise={kwargs_fig['noisetype']} 
            err_level={flux_err} 
            contrast={kwargs_fig['noisetype']} 
            limbdark={kwargs_run['u1']}""",
        xy=(-2, -1), fontsize=8)
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver2.png", bbox_inches="tight", dpi=100, transparent=True)

    return map_av

def solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, kwargs_fig, lr=0.001, niter=5000, annotate=False, colorbar=True):
    print("*** Using solver LSD+starry_opt ***")
    print(f"ydeg = {kwargs_run['ydeg']}")

    mean_profile = np.median(intrinsic_profiles, axis=0) # mean over chips
    observation_norm = np.median(obskerns_norm, axis=1) # mean over chips

    # preprae data
    pad = 100
    model_profile = 1. - np.concatenate((np.zeros(pad), mean_profile, np.zeros(pad)))
    dlam = np.diff(wav_angs[0]).mean() / 10 # angstrom to nm
    lam_ref = wav_angs[0].mean() / 10 # angstrom to nm
    nw = len(model_profile)
    wav0_lsd = np.linspace(start=lam_ref-0.5*dlam*nw, stop=lam_ref+0.5*dlam*nw, num=nw)
    wav_lsd = wav0_lsd[pad:-pad]

    with pm.Model() as model:
        # The surface map
        A = starry.DopplerMap(ydeg=kwargs_run['ydeg']).sht_matrix(smoothing=0.075)
        p = pm.Uniform("p", lower=0.0, upper=1.0, shape=(A.shape[1],))
        y = tt.dot(A, p)

        map_pm = starry.DopplerMap(lazy=True, wav=wav_lsd, wav0=wav0_lsd, **kwargs_run)

        map_pm[:, :] = y
        map_pm.spectrum = model_profile
        map_pm[1] = kwargs_run['u1']

        model_flux = map_pm.flux(theta=kwargs_run['theta'], normalize=True)

        # Likelihood term
        pm.Normal(
            f"obs",
            mu=tt.reshape(model_flux, (-1,)),
            sd=flux_err,
            observed=np.flip(observation_norm,axis=1).reshape(-1,)
        )

    # Optimize!
    loss = []
    best_loss = np.inf
    map_soln = model.test_point
    iterator = tqdm(
        pmx.optim.optimize_iterator(pmx.optim.Adam(lr=lr), niter, start=map_soln),
        total=niter,
        disable=os.getenv("CI", "false") == "true",
    )
    with model:
        for obj, point in iterator:
            iterator.set_description(
                "loss: {:.3e} / {:.3e}".format(obj, best_loss)
            )
            loss.append(obj)
            if obj < best_loss:
                best_loss = obj
                map_soln = point

    # Plot the loss
    fig, ax = plt.subplots(1, figsize=(3, 3))
    ax.plot(loss[int(len(loss)/20):])
    ax.set_xlabel("iteration number")
    ax.set_ylabel("loss")
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver3_loss.png", bbox_inches="tight", dpi=100, transparent=True)

    # Plot the MAP map
    map_res = starry.Map(kwargs_run['ydeg'], kwargs_run['udeg'], inc=kwargs_run['inc'], lazy=False)
    map_res[1] = kwargs_run['u1']
    image = map_res.render(projection="rect")

    with model:
        y_map = pmx.eval_in_model(y, point=map_soln)

    map_res[:, :] = y_map   # The spherical harmonic coefficient vector. 
    fig, ax = plt.subplots()
    map_res.show(ax=ax, projection="moll", colorbar=colorbar, cmap="plasma")
    if annotate:
        ax.annotate(f"""
            chip={kwargs_fig['goodchips']}
            solver=LSD+starry_opt(lr={lr}) 
            noise={kwargs_fig['noisetype']} 
            err_level={flux_err} 
            contrast={kwargs_fig['contrast']} 
            limbdark={kwargs_run['u1']}""",
        xy=(-2, -1), fontsize=8)
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver3.png", bbox_inches="tight", dpi=100, transparent=True)
    
    return map_res, image

def solve_starry_lin(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, annotate=False, colorbar=True):
    print("*** Using solver starry_lin***")
    goodchips = kwargs_fig['goodchips']
    nchip = len(goodchips)

    maps = [None for i in range(nchip)]
    images = []
    recimages = []
    successchips = []
    for i, jj in enumerate(goodchips):
        maps[i] = starry.DopplerMap(lazy=False, wav=wav_nm[i], wav0=wav0_nm[i], **kwargs_run)
        maps[i].spectrum = mean_spectrum[i]
        maps[i][1] = kwargs_run['u1']

        try:
            print(f"Solving chip {jj}... [{i+1}/{nchip}]")
            soln = maps[i].solve(
                flux=observed[:,i,:],
                theta=kwargs_run['theta'],
                normalized=True,
                fix_spectrum=True,
                flux_err=flux_err,
                quiet=os.getenv("CI", "false") == "true",
            )
            imag = maps[i].render(projection="moll")
            recimag = maps[i].render(projection="rect")
            images.append(imag)
            recimages.append(recimag)
            successchips.append(jj)

            print("Success!")
        
        except:
            print(f"Solver failed for chip {jj}, moving onto next chip...")
            continue

    # plot map of each chip
    if nchip > 1:
        fig, axs = plt.subplots(len(successchips), 1)
        for i, jj in enumerate(successchips):
            maps[0].show(ax=axs[i], projection="moll", image=images[i], colorbar=False)
            axs[i].annotate(f"chip {jj}", xy=(-1.6, -1))
        plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver4_each.png", bbox_inches="tight", dpi=100, transparent=True)

    # plot chip-averaged map
    images = np.array(images)
    recimages = np.array(recimages)

    fig, ax = plt.subplots()
    maps[0].show(ax=ax, projection="moll", image=np.mean(images, axis=0), colorbar=colorbar)
    if annotate:
        ax.annotate(f"""
            chip=median{goodchips} 
            solver=starry_lin 
            noise={kwargs_fig['noisetype']} 
            err_level={flux_err} 
            contrast={kwargs_fig['contrast']} 
            limbdark={kwargs_run['u1']}""",
        xy=(-2, -1), fontsize=8)
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver4.png", bbox_inches="tight", dpi=100, transparent=True)

    return np.mean(images, axis=0), np.mean(recimages, axis=0)

def solve_starry_opt(mean_spectrum, observed, wav_nm, wav0_nm, kwargs_run, kwargs_fig, lr=0.05, niter=5000, annotate=False, colorbar=True):
    print("*** Using solver starry_opt ***")
    goodchips = kwargs_fig['goodchips']
    nchip = len(goodchips)

    with pm.Model() as model:
        # The surface map
        A = starry.DopplerMap(ydeg=kwargs_run['ydeg']).sht_matrix(smoothing=0.075)
        p = pm.Uniform("p", lower=0.0, upper=1.0, shape=(A.shape[1],))
        y = tt.dot(A, p)

        maps = [None for i in range(nchip)]
        flux_models = [None for i in range(nchip)]
        for i, jj in enumerate(goodchips):
            print(f"Setting chip {jj} ({i+1}/{nchip})...")
            maps[i] = starry.DopplerMap(lazy=True, wav=wav_nm[i], wav0=wav0_nm[i], **kwargs_run)
            maps[i][:, :] = y
            maps[i].spectrum = mean_spectrum[i]
            maps[i][1] = kwargs_run['u1']

            flux_models[i] = maps[i].flux(theta=kwargs_run['theta'], normalize=True)

            # Likelihood term
            pm.Normal(
                f"obs{i}",
                mu=tt.reshape(flux_models[i], (-1,)),
                sd=flux_err,
                observed=observed[:,i,:].reshape(-1,))

    # Optimize!
    loss = []
    best_loss = np.inf
    map_soln = model.test_point
    iterator = tqdm(
        pmx.optim.optimize_iterator(pmx.optim.Adam(lr=lr), niter, start=map_soln),
        total=niter,
        disable=os.getenv("CI", "false") == "true",
    )
    with model:
        for obj, point in iterator:
            iterator.set_description(
                "loss: {:.3e} / {:.3e}".format(obj, best_loss)
            )
            loss.append(obj)
            if obj < best_loss:
                best_loss = obj
                map_soln = point

    # Plot the loss
    fig, ax = plt.subplots(1, figsize=(3, 3))
    ax.plot(loss[int(len(loss)/10):])
    ax.set_xlabel("iteration number")
    ax.set_ylabel("loss")
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver5_loss.png", bbox_inches="tight", dpi=100, transparent=True)

    # Plot the MAP map
    map_res = starry.Map(**kwargs_run)
    map_res[1] = kwargs_run['u1']

    with model:
        y_map = pmx.eval_in_model(y, point=map_soln)

    map_res[:, :] = y_map   # The spherical harmonic coefficient vector. 
    fig, ax = plt.subplots()
    map_res.show(ax=ax, projection="moll", colorbar=colorbar)
    if annotate:
        ax.annotate(f"""
            chip={goodchips}
            solver=starry_opt(lr={lr}) 
            noise={kwargs_fig['noisetype']} 
            err_level={flux_err} 
            contrast={kwargs_fig['contrast']} 
            limbdark={kwargs_run['u1']}""",
        xy=(-2, -1), fontsize=8)
    plt.savefig(paths.figures / f"{kwargs_fig['savedir']}/solver5.png", bbox_inches="tight", dpi=100, transparent=True)

    return map_res


################################################################################
####################   Utils    ################################################
################################################################################

def dao_getlines(f_linelist):
    """
    Read the line locations and equivalent widths from a DAOSPEC output file.

    Example:
      f_linelist = 'model_spec.clines'
      (lineloc, lineew, linespec) = getlines(f_linelist)
    """
    #2009-02-22 10:15 IJC: Initiated

    # Get the line locations and EWs:
    with open(f_linelist, 'r') as f:
        raw = f.readlines()

    #EB - this loop gets the form needed
    dat = np.zeros([len(raw), 2], dtype=float)                                                 
    for i, line in enumerate(raw):                                         
        dat[i,:]= list(map(float, line.split()[0:2]))

    lineloc = dat[:,0]
    lineew = dat[:,1]/1e3
    linespec = [line.split()[-1] for line in raw]
    return (lineloc, lineew, linespec)

def make_deltaspec(loc, ew, win, **kw):
    """
    Create a delta-function line spectrum based on a wavelength grid
    and a list of line locations and equivalent widths.

    :INPUTS:
       loc -- location of lines in the emission frame of reference

       ew  -- equivalent widths of lines, in units of wavelength grid.
               Positive values are emission lines.

       w_in -- wavelength grid in the emission frame, with values
              monotonically increasing (best if it is linearly spaced)

       All inputs should be lists or one-dimensional arrays of scalars

    :OPTIONAL_INPUTS:
       cont=None -- set continuum values in the emission frame;

       nearest=False  -- if True, use full pixels instead of partial

       verbose=False  -- if True, print out various messages

    :OUTPUTS:
      s  -- delta-function line spectrum, with a continuum level of zero
    
    :EXAMPLE: (NEEDS TO BE UPDATED!):
       ::

          w   = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7]
          loc = [2.1, 2.35, 2.62]
          ew  = [0.1, .02, .01]
          s = linespec(loc, ew, w)
          print s  #  --->  [0, 1, 0, 0.1, 0.1, 0, 0.08, 0.02]

    :NOTE:  This may give incorrect results for saturated lines.
    """
    # 2008-12-05 13:31 IJC: Created
    # 2008-12-10 13:30 IJC: Added continuum option, reworked code some.
    # 2008-12-12 12:33 IJC: Removed RV option

    #from pylab import find
    from numpy import where as find

    # Check inputs
    loc = np.array(loc).copy().ravel()
    ew  = np.array(ew ).copy().ravel()
    win = np.array(win).copy().ravel()

    defaults = dict(cont=None, nearest=False, verbose=False)
    for key in defaults:
        #if (not kw.has_key(key)): #EB - This method does not work in python3
        if (not key in kw): #EB update
            kw[key] = defaults[key]
    verbose = bool(kw['verbose'])
    nearest = bool(kw['nearest'])
    contset = kw['cont']!=None

    if contset.all(): #EB update .all() to work in python3
        cont = np.array(kw['cont']).copy()
        if len(cont)!=len(win):
            print( "Wavelength grid and continuum must have the same length!")
            return -1
    else:
        cont = np.ones(win.shape)

    nlines = len(loc)
    if nlines != len(ew):
        if verbose:  print( "len(loc)>>" + str(len(loc)))
        if verbose:  print( "len(ew)>>" + str(len(ew)))
        print( "Line locations and equivalent widths must have same length!")
        return -1

    #Only use lines in the proper wavelength range
    nlineinit = len(loc)
    lind = (loc>=win.min()) * (loc<=win.max())
    loc = loc[lind]
    ew  =  ew[lind]
    nlines = len(loc)

    s = cont.copy()
    d = np.diff(win).mean()

    if verbose:  print( "s>>" + str(s))

    for ii in range(nlines):
        lineloc = loc[ii]
        lineew  = ew[ii]
        index = (win<lineloc).sum() - 1
        if nearest:
            s[index+1] = s[index]-cont[index]*lineew/d
        elif index==len(win):
            s[index] = s[index] - cont[index]*lineew/d
        else:
            s[index] = s[index] - lineew*cont[index]* \
                (win[index+1] - lineloc)/d/d
            s[index+1] = s[index+1] - lineew*cont[index+1] * \
                (lineloc - win[index])/d/d
        
        if verbose:  
            print( "(lineloc, lineew)>>" + str((lineloc, lineew)))
            print( "(index, d)>>" + str((index,d)))

    if verbose:
        print( "(nlineinit, nline)>>" + str((nlineinit, nlines)))
    return s

def rconvolve1d(a, b, mode='valid', extend='nearest'):
    """
    Compute a 1D reverse convolution in the style of Bramich 2008.

    :INPUTS:
        'a' should be longer than 'b' -- i.e., 'b' is the kernel.
        'extend' tells how to extend the boundaries -- either
        'nearest'-neighbor or a number

    :NOTES:
      This is "reversed" from the canonical definition of the convolution.

    :SEE_ALSO:   
      :func:`dsa`
    """
    # 2008-11-14 18:26 IJC: Created
    na = len(a)
    nb = len(b)
    n = max(na, nb)
    
    dx = int(np.floor(nb/2))

    if extend=='nearest':
        X = a[-1]
    else:
        X = extend

    a2 = X + np.zeros(na+nb-1, dtype=float)
    a2[dx:dx+na] = a
    #a2 = concatenate((a, X + zeros(nb-1)))
    #c = zeros(n, dtype='float')
    
    bmat = np.tile(b, (n,1))
    amat = np.zeros((n, nb), dtype='float')
    for ii in range(na):
        amat[ii,:] = a2[range(ii,ii+nb)]

    c = np.sum(amat * bmat, axis=1)
        
    return c

def dsa(r, i, Nk, **kw): #, w=None, verbose=False, noback=False):
    """
    Computational tool for Difference Spectral Analysis (DSA)
    
    :INPUTS:
       R -- reference spectrum.  This should have the highest possible
            signal-to-noise and the highest spectral resolution.
       I -- Current spectrum to be analysed.
       Nk -- number of pixels in the desired convolution kernel

    :OPTIONS:
       w       -- weights of the pixel values in I; typically (sigma)^-2
           (Not HANDELED CORRECTLY?!?!?)
       noback  -- do not fit for a variable background; assume constant.
       tol=1e-10 -- if matrix determinant is less than tol, use
                    pseudoinverse rather than straight matrix
                    inversion
       verbose -- Print output statements and make a plot or two
       retinv -- return a fourth output, the Least Squares inverse
                 matrix (False by default)

    :OUTPUTS:       (M, K, B, C):
       M -- R, convolved to match I
       K -- kernel used in convolution
       B -- background offset
       C -- chisquared of fit. If no weights were specified, weights
            are set to unity for this calculation.

    :OPTIONS:
       I -- inverse matrix

    :NOTES:
        Best results are obtained with proper registration of the spectra.
        Also, beware of edge effects.  As a general rule, anything within
        a kernel width of the edges is suspect.
        Also

    :SEE_ALSO:  
       :func:`dsamulti`    """
#    """Testing Bramich's algorithm for 1D spectra."""
#    Based on the 2D Bramich (2008) DIA algorithm
#    -----
#    2008-11-14 10:56 IJC: Created @ UCLA.
#    2008-11-18 11:12 IJC: Registration now works correctly
#    2008-12-09 16:10 IJC: Somewhat optimized
#    2009-02-26 22:06 IJC: Added retinv, changed optional input format

    defaults = dict(verbose=False, w=None, noback=False, tol=1e-10, \
                        retinv=False)
    for key in defaults:
        if (not (key in kw)):
            kw[key] = defaults[key]
    verbose = bool(kw['verbose'])
    noback = bool(kw['noback'])
    retinv = bool(kw['retinv'])
    w = kw['w']
    if verbose:
        print( "kw>>" + str(kw))

    if noback:
        if verbose: print( "Not fitting for a variable background...")

    tol = 1e-10  # tolerance for singularity

    r = np.array(r, copy=True)
    i = np.array(i, copy=True)
    Nk = int(Nk)  # length of kernel
    dx = int(np.floor(Nk/2))

    if w==None:
        w = np.ones(len(r), dtype=float)

    Nr = len(r)  # length of Referene
    ind = np.arange(Nr-Nk+1, dtype=int)
    wind = w[ind]
        
    if noback:    
        U = np.zeros((Nk,Nk), dtype=float)
        b = np.zeros(Nk, dtype=float)
    else:
        U = np.zeros((Nk+1,Nk+1), dtype=float)
        b = np.zeros(Nk+1, dtype=float)

    # Build the b vector and U matrix
    tempval0 = w[ind+dx] * i[ind+dx]
    for p in range(Nk):
        b[p] = (tempval0 * r[ind+p]).sum()
        tempval2 = wind*r[ind+p]
        for q in range(p, Nk):
            U[p,q] = (tempval2 * r[ind+q]).sum()
            U[q,p] = U[p,q]

    if not noback:
        b[Nk] = (w[ind+dx] * i[ind+dx]).sum()
        for q in range(Nk):
            U[Nk, q] = (wind * r[ind+q]).sum()
            U[q, Nk] = U[Nk, q]

        U[Nk,Nk] = wind.sum()
    
    detU = np.linalg.det(U)
    if verbose: print( "det(U) is:  " + str(detU))

    if detU<tol:
        print( "Singular matrix: det(U) < tol.  Using pseudoinverse...")
        if verbose: 
            print( 'U>>',U)
        invmat = np.linalg.pinv(U)
    else:
        invmat = np.linalg.inv(U)

    a = np.dot(invmat, b)

    if noback:
        K = a
        B0 = 0.0
    else:
        K = a[0:len(a)-1]
        B0 = a[-1]

    m = rconvolve1d(r, K, mode='valid') + B0

    chisq  = ( wind * (i[ind] - m[ind])**2 ).sum()

    if verbose:
        chisq0 = ( wind * (i[ind] - r[ind])**2 ).sum()
        #print "Kernel is:  " + str(K)
        print( "Background: " + str(B0))
        print( "For the (" + str(Nr) + " - " + str(Nk+1) + ") = " + str(Nr-Nk-1) + " DOF:")
        print( "Red. Chisquared (I-R): " + str(chisq0/(Nr-Nk-1)))
        print( "Red. Chisquared (I-M): " + str(chisq/(Nr-Nk-1)))
    
        plt.figure(); plt.subplot(311)
        plt.plot(r, '--'); plt.plot(i, '-x'); plt.plot(m, '-..'); plt.legend('RIM'); 
        plt.subplot(312); plt.plot(r - i, '--'); plt.plot(m - i); plt.legend(['R-I', 'M-I'])
        plt.subplot(313); plt.plot(K, '-o'); plt.grid('on'); plt.legend(['Kernel']); 

    if retinv:
        return (m, K, B0, chisq, invmat)
    else:
        return (m, K, B0, chisq)

def plot_timeseries(map, modelspec, theta, obsflux=None, overlap=8.0, figsize=(5, 11.5)):
    # Plot the "Joy Division" graph
    fig = plt.figure(figsize=figsize)
    ax_img = [
        plt.subplot2grid((map.nt, 8), (t, 0), rowspan=1, colspan=1)
        for t in range(map.nt)
    ]
    ax_f = [plt.subplot2grid((map.nt, 8), (0, 1), rowspan=1, colspan=7)]
    ax_f += [
        plt.subplot2grid(
            (map.nt, 8),
            (t, 1),
            rowspan=1,
            colspan=7,
            sharex=ax_f[0],
            sharey=ax_f[0],
        )
        for t in range(1, map.nt)
    ]

    for t in range(map.nt):
        map.show(theta=theta[t], ax=ax_img[t], res=300)

        for l in ax_img[t].get_lines():
            if l.get_lw() < 1:
                l.set_lw(0.5 * l.get_lw())
                l.set_alpha(0.75)
            else:
                l.set_lw(1.25)
        ax_img[t].set_rasterization_zorder(100)

        # plot the obs data points
        if obsflux is not None:
            ax_f[t].plot(obsflux[t] - modelspec[t, 0], "k.",
                        ms=0.5, alpha=0.75, clip_on=False, zorder=-1)
        # plot the spectrum
        ax_f[t].plot(modelspec[t] - modelspec[t, 0], "C1-", lw=0.5, clip_on=False)

        ax_f[t].set_rasterization_zorder(0)
    fac = (np.max(modelspec) - np.min(modelspec)) / overlap
    ax_f[0].set_ylim(-fac, fac)
 
def plot_map_cells(map_obj):
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.grid(True)
    good = (map_obj.projected_area>0)
    for k in range(map_obj.ncell):
        lats = map_obj.corners_latlon[k][0]
        lons = map_obj.corners_latlon[k][1]

        y = np.array([lats[0], lats[1], lats[3], lats[2]]) - np.pi/2
        x = np.array([lons[0], lons[1], lons[3], lons[2]]) - np.pi
        # Plot the polygon
        if good[k]:
            poly = plt.Polygon(np.column_stack((x, y)), facecolor='gray', edgecolor='black')
            ax.add_patch(poly)
            ax.text(x.mean(), y.mean(), f"{k}", size=5)
        #ax.text(x.mean()-0.1, y.mean()-0.07, f"a:{map_obj.projected_area[k]:.3f}", size=3)

    # Set plot parameters
    ax.set_xticklabels([30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330], fontsize=8)

def solve_DIME(
        obskerns_norm: np.ndarray, 
        intrinsic_profiles: np.ndarray,
        dbeta: float, 
        nk: int, nobs: int, 
        phases: np.ndarray, 
        inc: float, vsini: float, LLD: float, 
        eqarea: bool = True,
        nlat: int = 20, nlon: int = 40,
        alpha: int = 4500, ftol: float = 0.01,
        plot_cells: bool = False,
        plot_unstretched_map: bool = False,
        spotfit: bool = False,
        create_obs_from_diff: bool = True
) -> np.ndarray:
    """
    Copied from IC14orig except kerns used to compute weights should take 
    input from cen_kerns (profiles centered to rv=0).
    ***inc in degrees (90 <-> equator-on).***

    Parameters
    ----------
    obskerns_norm : 3darray, shape=(nobs, nchip, nk)
        The observed line profiles (kerns).

    intrinsic_profiles : 2darray, shape=(nchip, nk)
        The model line profiles (modekerns).

    dbeta: float
        d_lam/lam_ref of the wavelength range that the line profile sits on.

    nk: int
        Size of line profile kernel.

    nobs: int
        Number of observations.

    phases: 1darray, shape=(nobs)
        Phases corresponding to the obs timesteps. In radian (0~2*pi).

    inc: float
        Inclination of star in degrees (common definition, 90 is equator-on)

    Returns
    -------
    bestparamgrid: 2darray, shape=(nlon, nlat)
        Optimized surface map. 
        Cells corresponding to longitude 0~2*pi, latitude 0~pi.

    """
    # Can safely take means over chips now
    mean_profile = np.median(intrinsic_profiles, axis=0) # mean over chips
    observation_2d = np.median(obskerns_norm, axis=1)
    observation_1d = observation_2d.ravel() # mean over chips and ravel to 1d

    # calc error for each obs
    smoothed = savgol_filter(obskerns_norm, 31, 3)
    resid = obskerns_norm - smoothed
    err_pix = np.array([np.abs(resid[:,:,pix] - np.median(resid, axis=2)) for pix in range(nk)]) # error of each pixel in LP by MAD
    err_LP = 1.4826 * np.median(err_pix, axis=0) # error of each LP (at each t and chip)
    err_each_obs = err_LP.mean(axis=1)
    err_observation_norm = np.tile(err_each_obs[:, np.newaxis], (1,nk)).ravel() # look like a step function over different times

    ### Prepare data for DIME
    modIP = 1. - np.concatenate((np.zeros(300), mean_profile, np.zeros(300)))
    modDV = - np.arange(np.floor(-modIP.size/2.+.5), np.floor(modIP.size/2.+.5)) * dbeta * const.c / 1e3
    flineSpline = interpolate.UnivariateSpline(modDV[::-1], modIP[::-1], k=1., s=0.)
    dv = -dbeta * np.arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5)) * const.c / 1e3 # km/s

    ### Reconstruct map

    # initialize Doppler map object
    inc_ = (90 - inc) * np.pi / 180 # IC14 defined 0 <-> equator-on, pi/2 <-> face-on
    if eqarea:
        mmap = ELL_map.map(nlat=nlat, nlon=nlon, type='eqarea', inc=inc_, verbose=True)
    else:
        mmap = ELL_map.map(nlat=nlat, nlon=nlon, inc=inc_) #ELL_map.map returns a class object
    if plot_cells:
        plot_map_cells(mmap)
    ncell = mmap.ncell
    nx = ncell
    dime.setup(observation_1d.size, nk)
    flatguess = 100*np.ones(nx)
    bounds = [(1e-6, 300)]*nx
    allfits = []

    # Compute R matrix
    Rmatrix = np.zeros((ncell, nobs*dv.size), dtype=np.float32) 
    for kk, rot in enumerate(phases):
        speccube = np.zeros((ncell, dv.size), dtype=np.float32) 
        if eqarea:
            this_map = ELL_map.map(nlat=nlat, nlon=nlon, type='eqarea', inc=inc_, deltaphi=-rot)
        else:
            this_map = ELL_map.map(nlat=nlat, nlon=nlon, inc=inc_, deltaphi=-rot)
        this_doppler = 1. + vsini*this_map.visible_rvcorners.mean(1)/const.c/np.cos(inc_) # mean rv of each cell in m/s
        good = (this_map.projected_area>0) * np.isfinite(this_doppler)
        for ii in good.nonzero()[0]:
            speccube[ii,:] = flineSpline(dv + (this_doppler[ii]-1)*const.c/1000.)
        limbdarkening = (1. - LLD) + LLD * this_map.mu
        Rblock = speccube * ((limbdarkening*this_map.projected_area).reshape(this_map.ncell, 1)*np.pi/this_map.projected_area.sum())
        Rmatrix[:,dv.size*kk:dv.size*(kk+1)] = Rblock

    flatmodel = dime.normalize_model(np.dot(flatguess, Rmatrix), nk)
    flatmodel_2d = np.reshape(flatmodel, (nobs, nk))

    # create diff+flat profile
    nchip = obskerns_norm.shape[1]
    uniform_profiles = np.zeros((nchip, nk))
    for c in range(nchip):
        uniform_profiles[c] = obskerns_norm[:,c].mean(axis=0) # time-avged LP for each chip
    mean_dev = np.median(np.array([obskerns_norm[:,c]-uniform_profiles[c] for c in range(nchip)]), axis=0) # mean over chips
    new_observation_2d = mean_dev + flatmodel_2d
    new_observation_1d = new_observation_2d.ravel()

    if len(allfits)==0:  # Properly scale measurement weights:
        #minfm = flatmodel.min()
        #cutoffval = 1. - (1. - minfm) / 22.

        # Mask out non-surface velocity space with weight=0
        width = int(vsini/1e3/np.abs(np.diff(dv).mean())) + 15 # vsini edge plus uncert=3
        central_indices = np.arange(nobs) * nk + int(nk/2)
        mask = np.zeros_like(observation_1d, dtype=bool)
        for central_idx in central_indices:
            mask[central_idx - width:central_idx + width + 1] = True
        w_observation = (mask == True).astype(float) / err_observation_norm**2

        if create_obs_from_diff:
            observation_1d = new_observation_1d

        # Scale the observations to match the model's equivalent width:
        out, eout = an.lsq((observation_1d, np.ones(nobs*nk)), flatmodel, w=w_observation)
        sc_observation_1d = observation_1d * out[0] + out[1]

        # perfect fit is when alpha=0 i.e. no smoothing
        #fitargs = (sc_observation_norm, w_observation, Rmatrix, 0)
        #perfect_fit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, maxiter=1e4, bounds=bounds)
        # the metric to be minimized is (0.5*chisq - alpha*entropy)
        #perfect_model = dime.normalize_model(np.dot(perfect_fit[0], Rmatrix), nk)
        #w_observation /=  w_observation.max() * (sc_observation_norm - perfect_model)[w_observation>0].std()**2

    cc = None
    if spotfit and not eqarea:
        print("Running MCMC spot fitting...")
        nstep = 2500
        # Do a one-spot fit:
        guess = [100,              80,               30,          180,          30          ] 
        #        100?,             spot_brightness%, spotlat=46°, spotlon=132°, spotradius=28°
        limits = [[99.99, 100.01], [0, 200],        [-90, 90],   [0, 360],     [10, 50]] # uniform prior
        #limits = [[100, 0.001],    [90, 20],        [30, 20],    [180, 100],   [30, 20]]
        spotargs0 = (mmap.corners_latlon.mean(2)[:,1].reshape(nlat, nlon), mmap.corners_latlon.mean(2)[:,0].reshape(nlat, nlon) - np.pi/2., Rmatrix)
        spotargs = (dime.profile_spotmap,)+spotargs0 + (sc_observation_1d, w_observation, dict(uniformprior=limits))
        thisfit = an.fmin(an.errfunc, guess, args=spotargs, full_output=True)
        test = dime.profile_spotmap(guess, *spotargs0)
        
        def lnprobfunc(*arg, **kw):
            ret = -0.5 * an.errfunc(*arg, **kw)
            if 'nans_allowed' in kw and (not kw.pop('nans_allowed')) or not (np.isfinite(ret)):
                print( "Whoops -- nan detected, but nans NOT ALLOWED in lnprobfunc!")
                ret = 9e99
            return ret
        # MCMC for one-spot fit:
        ndim = len(guess)
        nwalkers = ndim*30
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobfunc, args=spotargs, threads=6)
        junk = np.abs(thisfit[0])/100.
        junk[0] = 0.01
        p0, testchi = get_emcee_start(thisfit[0], junk, nwalkers, thisfit[1]*2, spotargs, retchisq=True)

        #p1, prob, state = sampler.run_mcmc(p0, nstep) # Burn-in
        #sampler.reset()
        p2, prob, state = sampler.run_mcmc(p0, nstep)
        spotparams = sampler.flatchain[np.nonzero(sampler.lnprobability.ravel()==sampler.lnprobability.ravel().max())[0][0]]

        if False:
            # Do a two-spot fit:
            guess2 = np.concatenate((bestparams, bestparams[1:])) 
            guess2[5] = 2*guess2[0] - guess2[1]
            guess2[6] = -guess2[6]
            guess2[7] = (guess2[3] + np.pi) % (2*np.pi)

            limits2 = limits + limits[1:]
            spotargs2 = (dime.profile_spotmap,)+spotargs0 + (sc_observation_norm, w_observation, dict(uniformprior=limits2))
            thisfit = an.fmin(pc.errfunc, guess2, args=spotargs2, full_output=True)
            test = dime.profile_spotmap(guess2, *spotargs0)

            # MCMC for two-spot fit:
            ndim2 = len(guess2)
            nwalkers2 = ndim2*30
            sampler2 = emcee.EnsembleSampler(nwalkers2, ndim2, pc.lnprobfunc, args=spotargs2, threads=3)
            junk = np.abs(thisfit[0])/100.
            junk[0] = 0.01
            p0, testchi = tools.get_emcee_start(thisfit[0], junk, nwalkers2, thisfit[1]*2, spotargs2, retchisq=True)

            p1, prob, state = sampler2.run_mcmc(p0, nstep) # Burn-in
            sampler2.reset()
            p2, prob, state = sampler2.run_mcmc(p1, nstep) # Burn-in
            spotparams2 = sampler2.flatchain[nonzero(sampler2.lnprobability.ravel()==sampler2.lnprobability.ravel().max())[0][0]]

        cc = sampler.flatchain.copy()
        if inc==0: 
            cc[:,2] = abs(cc[:,2])
        ind  = np.array([1,2,4])
        labs = ['Spot Brightness [%]', 'Spot Latitude [deg]', 'Spot Radius [deg]']
        spotmap = makespot(spotparams[2], spotparams[3], spotparams[4], *spotargs0[0:2])
        print(f"Spot brightness: {spotparams[1]:.0f}\n"
            f"Spot lat: {spotparams[2]:.0f}\n",
              f"spot lon: {spotparams[3]:.0f}\n",
              f"spot radius: {spotparams[4]:.0f}\n", 
              cc)
        plt.figure(figsize=(5,3))
        plt.imshow(spotmap, origin="lower", extent=(0, 360, -90, 90))
        
    ### Solve!
    fitargs = (sc_observation_1d, w_observation, Rmatrix, alpha)
    bfit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, disp=1, maxiter=1e4, bounds=bounds)
    allfits.append(bfit)
    bestparams = bfit[0]
    model_observation = dime.normalize_model(np.dot(bestparams, Rmatrix), nk)
    metric, chisq, entropy = dime.entropy_map_norm_sp(bestparams, *fitargs, retvals=True)
    print("metric:", metric, "chisq:", chisq, "entropy:", entropy)

    # reshape into grid
    if eqarea:
        # reshape into list
        start=0
        bestparamlist = []
        for m in range(this_map.nlat):
            bestparamlist.append(bestparams[start:start+this_map.nlon[m]])
            start = start + this_map.nlon[m]
        # interp into rectangular array
        max_length = max([len(x) for x in bestparamlist])
        stretched_arrays = []
        for array in bestparamlist:
            x_old = np.arange(len(array))
            x_new = np.linspace(0, len(array) - 1, max_length)
            y_new = np.interp(x_new, x_old, array)
            stretched_arrays.append(y_new)

        bestparamgrid = np.vstack(stretched_arrays)

        if plot_unstretched_map:
            # pad into rectangular array
            padded_arrays = []
            for array in bestparamlist:
                left_pad = int((max_length - len(array)) / 2)
                right_pad = max_length - len(array) - left_pad
                padded_array = np.pad(array, (left_pad, right_pad), 'constant')
                padded_arrays.append(padded_array)
                array_2d = np.vstack(padded_arrays)
                plt.imshow(array_2d, cmap='plasma')
                #plt.show()

    else:
        bestparamgrid = np.reshape(bestparams, (-1, nlon))

    res = dict(
        bestparams=bestparams,
        bestparamgrid=bestparamgrid, 
        Q=metric, chisq=chisq, entropy=entropy,
        Rmatrix=Rmatrix,
        model_observation=model_observation,
        sc_observation_1d=sc_observation_1d,
        w_observation=w_observation,
        flatmodel=flatmodel,
        flineSpline=flineSpline,
        mmap=mmap, dv=dv, dbeta=dbeta,
        cc=cc, sampler=sampler
    )
        
    return bestparamgrid, res

def get_emcee_start(bestparams, variations, nwalkers, maxchisq, args, homein=True, retchisq=False, depth=np.inf):
    """Get starting positions for EmCee walkers.

    :INPUTS:
      bestparams : sequence (1D NumPy array)
        Optimal parameters for your fitting function (length N)

      variations : 1D or 2D NumPy array
        If 1D, this should be length N and new trial positions will be
        generated using numpy.random.normal(bestparams,
        variations). Thus all values should be greater than zero!

        If 2D, this should be size (N x N) and we treat it like a
        covariance matrix; new trial positions will be generated using
        numpy.random.multivariate_normal(bestparams, variations). 

      nwalkers : int
        Number of positions to be chosen.

      maxchisq : int
        Maximum "chi-squared" value for a test position to be
        accepted.  In fact, these values are computed with
        :func:`phasecurves.errfunc` as errfunc(test_position, *args)
        and so various priors, penalty factors, etc. can also be
        passed in as keywords.

      args : tuple
        Arguments to be passed to :func:`phasecurves.errfunc` for
        computing 'chi-squared' values.

      homein : bool
        If True, "home-in" on improved fitting solutions. In the
        unlikely event that a randomly computed test position returns
        a better chi-squared than your specified best parameters,
        reset the calculation to start from the new, improved set of
        parameters.

      retchisq : bool
        If True, return the tuple (positions, chisq_at_positions)
        
     :BAD_EXAMPLE:
      ::

        pos0 = tools.get_emcee_start(whitelight_bestfit[0], np.abs(whitelight_bestfit[0])/1000., nwalkers, 10*nobs, mcargs)
        """
    # 2013-05-01 11:18 IJMC: Created
    # 2014-07-24 11:07 IJMC: Fixed typo in warning message.
    
    #get_emcee_start(bestparams, variations, nwalkers, maxchisq, args):

    best_chisq = an.errfunc(bestparams, *args)
    if best_chisq >= maxchisq:
        print("Specified parameter 'maxchisq' is smaller than the chi-squared value for the specified best parameters. Try increasing maxchisq.")
        return -1

    npar = len(bestparams)
    if variations.ndim==2:
        usecovar = True
    else:
        usecovar = False

    pos0 = np.zeros((nwalkers, npar), dtype=float)
    chisq = np.zeros(nwalkers, dtype=float)
    npos = 0
    while npos < nwalkers:
        if usecovar:
            testpos = np.random.multivariate_normal(bestparams, variations)
        else:
            testpos = np.random.normal(bestparams, variations)
        testchi = an.errfunc(testpos, *args)
        if np.isfinite(testchi) and (testchi < best_chisq) and homein and depth>0:
            return get_emcee_start(testpos, variations, nwalkers, maxchisq, args, homein=homein, retchisq=retchisq, depth=depth-1)
        elif testchi < maxchisq:
            pos0[npos] = testpos
            chisq[npos] = testchi
            npos += 1

    if retchisq:
        ret = pos0, chisq
    else:
        ret = pos0
    return ret

def makespot(spotlat, spotlon, spotrad, phi, theta):
    """
    :INPUTS:
      spotlat : scalar
        Latitude of spot center, in degrees, from 0 to 180 (actually from -90 to 90)

      spotlon : scalar
        Longitude of spot center, in degrees, from 0 to 360

      spotrad : scalar
        Radius of spot, in radians. degrees

      phi, theta : 2D NumPy arrays
         output from :func:`makegrid`.  Theta ranges from -pi/2 to +pi/2.

    :EXAMPLE:
      ::

        import maps
        nlat, nlon = 60, 30
        phi, theta = maps.makegrid(nlat, nlon)
        # Make a small spot centered near, but not at, the equator:
        equator_spot = maps.makespot(0, 0, 23, phi, theta)
        # Make a larger spot centered near, but not at, the pole:
        pole_spot = maps.makespot(68, 0, 40, phi, theta)

      ::

        import maps
        nlat, nlon = 60, 30
        map = maps.map(nlat, nlon, i=0., deltaphi=0.)
        phi = map.corners_latlon.mean(2)[:,1].reshape(nlon, nlat)
        theta = map.corners_latlon.mean(2)[:,0].reshape(nlon, nlat) - np.pi/2.
        # Make a small spot centered near, but not at, the equator:
        equator_spot = maps.makespot(0, 0, 23, phi, theta)
        # Make a larger spot centered near, but not at, the pole:
        pole_spot = maps.makespot(68, 0, 40, phi, theta)

    """
    # 2013-08-18 16:01 IJMC: Created

    spotlat *= (np.pi/180)
    spotlon *= (np.pi/180)
    spotrad *= (np.pi/180)

    pi2 = 0.5*np.pi
    xyz = np.array((np.cos(phi) * np.sin(theta + pi2), np.sin(phi) * np.sin(theta + pi2), np.cos(theta + pi2))).reshape(3, phi.size)

    # First rotate around z axis, to align spot with sub-observer meridian
    # Then, rotate around y axis, to align spot with pole.
    zrot = np.array([[np.cos(np.pi-spotlon), -np.sin(np.pi-spotlon), 0], [np.sin(np.pi-spotlon), np.cos(np.pi-spotlon), 0.], [0,0,1]])
    yrot = np.array([[np.cos(spotlat+pi2), 0, np.sin(spotlat+pi2)], [0,1,0], [-np.sin(spotlat+pi2), 0, np.cos(spotlat+pi2)]])
    xyz = np.dot(np.dot(yrot, zrot), xyz)

    # Convert Cartesian to spherical coordinates
    ang = np.arccos(xyz[2])

    # Spot is where (theta - theta_pole) < radius.
    spotmap = ang.T <= spotrad

    return spotmap.reshape(phi.shape)

def remove_spike(data, kern_size=10, lim_denom=5):
    data_pad = np.concatenate([np.ones(kern_size)*np.median(data[:kern_size]), data, np.ones(kern_size)*np.median(data[-kern_size:-1])])
    data_filt = np.copy(data)
    for i, val in enumerate(data):
        i_pad = i + kern_size
        seg = data_pad[i_pad-kern_size:i_pad+kern_size]
        seg = seg[np.abs(seg-np.median(seg))<20]
        lim = np.median(seg)/lim_denom
        if val > np.median(seg) + lim or val < np.median(seg) - lim:
            data_filt[i] = np.median(seg[int(kern_size/5):-int(kern_size/5)])
    return data_filt

def shift_kerns_to_center(modkerns, kerns, instru, goodchips, dv, sim=False):
    '''shift modkerns to center at dv=0 and shift kerns for same amount.'''
    nobs, nchip, nk = modkerns.shape
    cen_modkerns = np.zeros_like(modkerns)
    cen_kerns = np.zeros_like(kerns)
    for i,jj in enumerate(goodchips):
        for k in range(nobs):
            systematic_rv_offset = (modkerns[k,i]==modkerns[k,i].max()).nonzero()[0][0] - (dv==0).nonzero()[0][0] # find the rv offset
            cen_modkerns[k,i] = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, modkerns[k,i]) # shift ip to center at dv=0
            print("chip:", jj , "obs:", k, "offset:", systematic_rv_offset)
            cen_kerns[k,i] = kerns[k,i]
            if not sim: # shift kerns with same amount if not simulation
                if instru != 'CRIRES': # don't shift kerns if crires
                    if k==0:
                        print("kern shifted to same amount")
                    #cen_kerns[k,i] = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, kerns[k,i])
    return cen_modkerns, cen_kerns

def cont_normalize_kerns(cen_kerns, instru):
    '''Continuum-normalize kerns by fitting a line at the flat edges of kern.'''
    nobs, nchip, nk = cen_kerns.shape
    obskerns = 1. - cen_kerns
    obskerns_norm = np.zeros_like(obskerns)
    continuumfit = np.zeros((nobs, nchip, 2))
    side = 15 if instru != "CRIRES" else 7
    for i in range(nchip):
        for n in range(nobs):
            inds = np.concatenate((np.arange(0, side), np.arange(nk-side, nk)))
            continuumfit[n,i] = np.polyfit(inds, obskerns[n,i,inds], 1)
            obskerns_norm[n,i] = obskerns[n,i] / np.polyval(continuumfit[n,i], np.arange(nk))
    return obskerns_norm

def plot_kerns_timeseries(kerns, goodchips, dv, gap=0.03, normed=False, intrinsic_profiles=None):
    '''Plot time series of kerns.'''
    nobs, nchip, nk = kerns.shape
    colors = [cm.gnuplot_r(x) for x in np.linspace(0, 1, nobs+4)]
    plt.figure(figsize=(nchip*3,4))
    for i, jj in enumerate(goodchips):
        plt.subplot(1, nchip, i+1)
        for n in range(nobs):
            if not normed:
                plt.plot(dv, 1 - kerns[n,i] - gap*n, color=colors[n])
            else:
                plt.plot(dv, kerns[n,i] - gap*n, color=colors[n])
        if intrinsic_profiles is not None:
            plt.plot(dv, 1-intrinsic_profiles[i], color='k')
        plt.title(f"chip={jj}")
        plt.xlabel("dv")
    plt.tight_layout()

def plot_chipav_kern_timeseries(obskerns_norm, dv, timestamps, savedir, gap=0.025, cut=17):
    '''Plot time series of chip-averaged kerns.'''
    nobs = obskerns_norm.shape[0]
    colors = [cm.gnuplot_r(x) for x in np.linspace(0, 1, nobs+4)]
    fig, ax = plt.subplots(figsize=(4, 5))
    for n in range(nobs):
        ax.plot(dv[cut:-cut], obskerns_norm[n].mean(axis=0)[cut:-cut] - gap*n, color=colors[n+1])
        #plt.text(dv[cut] + 10, 1 - gap/4 - gap*n, f"{timestamps[n]:.1f}h")
    #plt.plot(dv, 1-intrinsic_profiles.mean(axis=0), color='k', label="intrinsic profile")
    ax.set_xlabel("velocity (km/s)")
    ax.set_xticks([-50, -25, 0, 25, 50])
    ax.set_ylabel("Line intensity")
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ybound())
    ax2.set_yticks([1- gap*n for n in range(nobs)], labels=[f"{t:.1f}h" for t in timestamps], fontsize=9)
    #plt.axvline(x=vsini/1e3, color="k", linestyle="dashed", linewidth=1)
    #plt.axvline(x=-vsini/1e3, color="k", linestyle="dashed", linewidth=1)
    #plt.legend(loc=4, bbox_to_anchor=(1,1))
    plt.savefig(paths.figures / f"{savedir}/tsplot.png", bbox_inches="tight", dpi=300, transparent=True)

def plot_deviation_map(obskerns_norm, goodchips, dv, vsini, timestamps, savedir, meanby="median", cut=30):
    '''Plot deviation map for each chip and mean deviation map'''
    nobs, nchip, nk = obskerns_norm.shape
    uniform_profiles = np.zeros((nchip, nk))

    plt.figure(figsize=(nchip*4,2.5))
    for i, jj in enumerate(goodchips):
        uniform_profiles[i] = obskerns_norm[:,i].mean(axis=0) # is each chip's mean kern over epoches
        plt.subplot(1,nchip+1,i+1)
        plt.imshow(obskerns_norm[:,i]-uniform_profiles[i], 
            extent=(dv.max(), dv.min(), timestamps[-1], 0),
            aspect=int(0.7*vsini/1e3),
            cmap='YlOrBr') # positive diff means dark spot
        plt.xlim(dv.min()+cut, dv.max()-cut),
        plt.xlabel("velocity (km/s)")
        plt.ylabel("Elapsed time (h)")
        plt.colorbar(fraction=0.035)
        plt.title(f"chip={jj}")
    if meanby == "median":
        mean_dev = np.median(np.array([obskerns_norm[:,i]-uniform_profiles[i] for i in range(nchip)]), axis=0) # mean over chips
    elif meanby == "median_each":
        mean_dev = np.median(obskerns_norm, axis=1) - np.median(uniform_profiles,axis=0)
    elif meanby == "mean":
        mean_dev = np.mean(np.array([obskerns_norm[:,i]-uniform_profiles[i] for i in range(nchip)]), axis=0) # mean over chips
    plt.subplot(1, nchip+1,i+2)
    plt.imshow(mean_dev, 
        extent=(dv.max(), dv.min(), timestamps[-1], 0),
        aspect=int(0.7* vsini/1e3),
        cmap='YlOrBr') # positive diff means dark spot
    plt.xlim(dv.min()+cut, dv.max()-cut),
    plt.xlabel("velocity (km/s)")
    plt.ylabel("Elapsed time (h)")
    plt.colorbar(fraction=0.035)
    plt.title(f"{meanby} deviation")
    plt.tight_layout()
    plt.savefig(paths.figures / f"{savedir}/tvplot_full.png", bbox_inches="tight", dpi=100, transparent=True)
    # plot only the mean map
    plt.figure(figsize=(5,3))
    plt.imshow(mean_dev, 
        extent=(dv.max(), dv.min(), timestamps[-1], 0),
        aspect=int(0.7* 29e3/1e3),
        cmap='YlOrBr') # positive diff means dark spot
    plt.xlim(dv.min()+cut, dv.max()-cut),
    plt.xlabel("velocity (km/s)")
    plt.xticks([-50, -25, 0, 25, 50])
    plt.ylabel("Elapsed time (h)")
    plt.vlines(x=vsini/1e3, ymin=0, ymax=timestamps[-1], colors="k", linestyles="dashed", linewidth=1)
    plt.vlines(x=-vsini/1e3, ymin=0, ymax=timestamps[-1], colors="k", linestyles="dashed", linewidth=1)
    plt.colorbar(fraction=0.036)
    #plt.title(f"{meanby} deviation")
    #plt.text(dv.min()+5, 0.5, f"chips={goodchips}", fontsize=8)
    plt.tight_layout()
    plt.savefig(paths.figures / f"{savedir}/tvplot.png", bbox_inches="tight", dpi=100, transparent=True)

def plot_IC14_map(bestparamgrid, colorbar=False, clevel=5, sigma=1):
    '''Plot doppler map from an array.'''
    nlat, nlon = bestparamgrid.shape 
    fig = plt.figure(figsize=(5,3))
    ax = fig.add_subplot(111, projection='mollweide')
    lon = np.linspace(-np.pi, np.pi, nlon)
    lat = np.linspace(-np.pi/2., np.pi/2., nlat)
    Lon,Lat = np.meshgrid(lon,lat)
    im = ax.pcolormesh(Lon, Lat, bestparamgrid, cmap=plt.cm.plasma, shading='gouraud')
    contour = ax.contour(Lon, Lat, gaussian_filter(bestparamgrid, sigma), clevel, colors='white', linewidths=0.5)
    if colorbar:
        fig.colorbar(im)
    ax.set_yticks(np.linspace(-np.pi/2, np.pi/2, 7), labels=[])
    ax.set_xticks(np.linspace(-np.pi, np.pi, 13), labels=[])
    ax.grid('major', color='k', linewidth=0.25)
    for item in ax.spines.values():
        item.set_linewidth(1.2)

def make_gif_map(bestparamgrid, inc, period, savedir, step=15, fps=4):
    fig = plt.figure(figsize=(10,5))
    y, x = bestparamgrid.shape

    # create an interpolation smoothed map
    plt.imshow(bestparamgrid, interpolation='bicubic', cmap='gray', aspect=x/y*0.5)
    plt.axis('off')
    fig.patch.set_visible(False)
    plt.savefig(paths.figures / f"{savedir}/solver1_map.png", bbox_inches='tight', dpi=200, pad_inches=0)
    img = Image.open(paths.figures / f"{savedir}/solver1_map.png")
    img = np.array(img.convert('L'), dtype='float64')
    img /= img.max()

    # make GIF with PIL
    Nlat, Nlon = img.shape
    Lon, Lat = np.meshgrid(np.linspace(-180, 180, Nlon), np.linspace(-90, 90, Nlat))

    from matplotlib.animation import FuncAnimation
    num_frames = int(360/step)
    fig = plt.figure(figsize=(4, 4))
    ax = plt.axes(projection=ccrs.Orthographic(0, 10))
    time_text = ax.text(0.85, 0.05,'', transform=ax.transAxes)
    def update(frame):
        ax = plt.axes(projection=ccrs.Orthographic(frame*step, 90-inc))
        gl = ax.gridlines(xlocs=range(-180, 180, 30), ylocs=range(-90, 90, 30), color='gray', linewidth=0.3)
        im = ax.imshow(img, origin="lower", extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=0, vmax=1)
        time_text.set_text(f'{frame*period/num_frames+0.2:.1f}h')
        for item in ax.spines.values():
            item.set_linewidth(1.5)
        return (ax, gl, im, time_text)

    # Create the animation
    ani = FuncAnimation(plt.gcf(), update, frames=num_frames)

    # Save the animation as a GIF
    output_gif_path = paths.figures / f"{savedir}/solver1_map.gif"
    ani.save(output_gif_path, dpi=100, fps=fps)
    '''
    frames = []
    for view in range(0, 360, step):
        fig = plt.figure(figsize=(4, 4))
        ax = plt.axes(projection=ccrs.Orthographic(view, 90-inc))
        for item in ax.spines.values():
            item.set_linewidth(1.5)
        gl = ax.gridlines(xlocs=range(-180,180,30), ylocs=range(-90,90,30), color='gray', linewidth=0.3)
        ax.imshow(img, origin="lower", extent=(-180, 180, -90, 90),transform=ccrs.PlateCarree(), cmap=plt.cm.plasma, vmin=0, vmax=1)  # Important

        # save fig in bytes format
        temp = paths.figures / f'temp/frame_{view}.png'
        plt.savefig(temp, format='png', dpi=150, bbox_inches='tight', transparent=True)
        #frames.append(imageio.imread(temp))
        frames.append(Image.open(temp))

    gif_file = paths.figures / f"{savedir}/solver1_map.gif"
    frames[0].save(gif_file, format='GIF', append_images=frames[1:], save_all=True, duration=duration, loop=0)
    '''

################################################################################
####################   Tests    ################################################
################################################################################

def test_sim():
    assert simulation_on == True
    global mean_spectrum, template, observed, residual, error, intrinsic_profiles, obskerns_norm
    global lin_map
    mean_spectrum, template, observed, residual, error = load_data(model_datafile)
    observed = spectra_from_sim(modelmap, contrast, roll, smoothing, fakemap_nlat, fakemap_nlon, mean_spectrum, error, residual, plot_ts=True)
    intrinsic_profiles, obskerns_norm = make_LSD_profile(template, observed)
    print(use_toy_spec)
    bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, annotate=True)
    #LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, avgchip=avgchip, annotate=True)
    LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, avgchip=avgchip, lr=lr_LSD, niter=niter_LSD, annotate=True)

    #lin_map = solve_starry_lin(mean_spectrum, observed, annotate=True)
    #opt_map = solve_starry_opt(mean_spectrum, observed, lr=lr, niter=niter, annotate=True)

def test_obs():
    assert simulation_on == False
    global mean_spectrum, template, observed, residual, error, intrinsic_profiles, obskerns_norm
    global bestparamgrid_r, bestparamgrid, LSDlin_map, LSDopt_map, lin_map, opt_map

    mean_spectrum, template, observed, residual, error = load_data(model_datafile)
    intrinsic_profiles, obskerns_norm = make_LSD_profile(template, observed)

    bestparamgrid_r, bestparamgrid = solve_IC14new(intrinsic_profiles, obskerns_norm, annotate=True)

    LSDlin_map = solve_LSD_starry_lin(intrinsic_profiles, obskerns_norm, kwargs_run, avgchip=avgchip, annotate=True)

    LSDopt_map = solve_LSD_starry_opt(intrinsic_profiles, obskerns_norm, kwargs_run, avgchip=avgchip, lr=lr_LSD, niter=2000, annotate=True)

    lin_map = solve_starry_lin(mean_spectrum, observed, annotate=True)
    #plt.figure(figsize=(5,3))
    #plt.savefig(paths.figures / "sim/solver4.png", bbox_inches="tight", dpi=300)

    opt_map = solve_starry_opt(mean_spectrum, observed, lr=lr, niter=2000, annotate=True)

    print("Test ok.")

if __name__ == "__main__":
    pass
    #test_obs()

'''
indices = np.arange(1, observation_norm.size, 5)
new_obs = np.delete(observation_norm, indices)
indices_p = np.arange(1, mean_profile.size, 5)
new_profile = np.delete(mean_profile, indices_p)

obs=np.reshape(res['sc_observation_norm'], (nobs, nk))
model=np.reshape(res['model_observation'], (nobs, nk))
flmodel = np.reshape(flatmodel, (nobs, nk))

plt.figure(figsize=(5, 7))
for i in range(nobs):
    plt.plot(obs[i] - 0.02*i, color='k', linewidth=1)
    #plt.plot(obs[i] - 0.02*i, '.', color='k', markersize=2)
    plt.plot(model[i] - 0.02*i, color='r', linewidth=1)
    plt.plot(flmodel[i] - 0.02*i, '--', color='gray', linewidth=1)
plt.legend(labels=['obs', 'best-fit map', 'flat map'])
'''