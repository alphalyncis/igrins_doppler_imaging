import numpy as np
import paths

##############################################################################
####################   Constants    ##########################################
##############################################################################

OPTIONS = {
    "modelmap": ["1spot", "1band", "2band", "blank"],
    "noisetype": ["residual", "random", "obserr", "none", "res+random"],
    "instru": ["IGRINS", "CRIRES"],
    "target": ["W1049A", "W1049B", "2M0036_1103", "2M0036_1105"],
    "band": ["K", "H", "both"],
    "solver": ["starry_lin", "starry_opt", "IC14new", "IC14orig"],
    "map_type": ["eqarea", "latlon"],
    "modelspec": ["lte015.0-5.0", "t1500g1000f8"],
    "LSD": ["new", "orig", None]
}

goodchips_sim = {
    "IGRINS": {
        "K": [1, 2, 3, 4, 5, 7, 12, 13, 15, 16, 18], #[1,2,4,13],
        "H": [3,4,5,6,18]
    },
    "CRIRES": {
        "K": [0,1,2,3],
        "H": [0,1,2,3]
    }
}

goodchips_sim = {
    "IGRINS": {
        "W1049B":{
            "K": [0, 1, 2, 3, 4, 5, 15, 16, 17, 18, 19], 
            "H": [0, 1, 2, 3, 4, 5, 16, 17, 18, 19] 
        },
        "W1049A":{
            "K": [0, 1, 2, 3, 4, 5, 13, 14, 15, 16, 17, 18], 
            "H": [1, 2, 3, 4, 5, 6, 17, 18, 19]
        },
        "2M0036_1103":{
            "K": [2,4,5,6,12,13,14,15,16,18,19],
            "H": [1,2,3,4,5,6,7,8,9,10,11,13,16,18] 
        },
        "2M0036_1105":{
            "K": [2,4,5,6,11,12,13,16,18], 
        },      
    },
    "CRIRES": {
        "W1049B":{
            "K": [0, 1, 2, 3]
        },
        "W1049A":{
            "K": [0, 1, 2, 3]
        }
    }
}

nobss =   {"W1049B": 14,     "W1049A": 14,     "2M0036_1103": 7   ,  "2M0036_1105": 8}

periods = {"W1049B": 5,   "W1049A": 7,      "2M0036_1103": 2.7 ,  "2M0036_1105": 2.7}
incs =    {"W1049B": 80,     "W1049A": 70,     "2M0036_1103": 51  ,  "2M0036_1105": 51}
vsinis =  {"W1049B": 29e3,   "W1049A": 21e3,   "2M0036_1103": 32e3,  "2M0036_1105": 32e3}
rvs =     {"W1049B": 7.05e-5, "W1049A": 5.4e-5, "2M0036_1103": 6e-5,  "2M0036_1105": 6e-5}
                #9e-5 9.3e-5if CRIRES 7.4e-5 5.4e-5 if IGRINS

##############################################################################
####################    Settings    ##########################################
##############################################################################

simulation_on = True
savedir = "sim"

#################### Simulation settings ############################
use_toy_spec = False
modelmap = "1spot"
noisetype = "res+random"

fakemap_nlat, fakemap_nlon = 181, 361
roll = 0.28 # 0-1
contrast = 0.5
smoothing = 0.1

#################### Run settings ####################################

instru = "IGRINS"
target = "W1049B"
band = "K"
#solver = "IC14new"
map_type = "eqarea"

modelspec = "t1500g1000f8"
LSD = "new"

########## IC14 parameters ##########
nk = 125 if instru != "CRIRES" else 203
cut = nk - 70
LLD = 0.4
alpha = 2000
ftol = 0.01 # tolerance for convergence of maximum-entropy
nstep = 2000
nlat, nlon = 10, 20

########## Starry parameters ##########
ydeg_sim = 16
ydeg = 8
udeg = 1
nc = 1
u1 = LLD
vsini_max = 40000.0

########## Starry optimization parameters ##########
lr_LSD = 0.001
niter_LSD = 5000
lr = 0.01
niter = 5000

#################### Automatic ####################################

if True:
    # Auto consistent options
    flux_err = 0.01 if use_toy_spec else 0.025

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
    period = periods[target]
    inc = incs[target]
    vsini = vsinis[target]
    veq = vsini / np.sin(inc * np.pi / 180)

    # set time and period parameters
    #timestamp = np.linspace(0, period, nobs)  # simulate equal time interval obs
    tobs = 5
    timestamp = np.linspace(0, tobs, nobs)
    phases = timestamp * 2 * np.pi / period # 0 ~ 2*pi in rad; IC14
    theta = 360.0 * timestamp / period      # 0 ~ 360 in degree; starry sim & run

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
