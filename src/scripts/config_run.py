import numpy as np
import paths

##############################################################################
####################   Constants    ##########################################
##############################################################################

OPTIONS = {
    "modelmap": ["1spot", "1band", "2band", "blank"],
    "noisetype": ["residual", "random", "obserr", "none", "res+random"],
    "instru": ["IGRINS", "CRIRES"],
    "target": ["W1049A", "W1049B", "W1049A_0209", "W1049B_0209", "2M0036_1103", "2M0036_1105"],
    "band": ["K", "H", "both"],
    "solver": ["starry_lin", "starry_opt", "IC14new", "IC14orig"],
    "map_type": ["eqarea", "latlon"],
    "modelspec": ["lte015.0-5.0", "t1500g1000f8"],
    "LSD": ["new", "orig", None]
}

goodchips_run = {
    "IGRINS": {
        "W1049B":{
            "K": [0, 1, 2, 3, 4, 5, 15, 16, 17, 18], #[1, 4, 13], #
            "H": [0, 1, 2, 3, 4, 5, 16, 17, 18, 19] 
        },
        "W1049A":{
            "K": [0, 1, 2, 3, 4, 5, 15, 16, 17, 18], 
            "H": [0, 1, 2, 3, 4, 5, 16, 17, 18, 19]
        },
        "W1049B_0209":{
            "K": [2, 3, 4, 5, 16, 17, 18, 19], # snr 20-30: 6, 7, 8, 13, 14, 15, 19
            "H": [0, 1, 2, 3, 4, 5, 16, 17, 18, 19] 
        },
        "W1049A_0209":{
            "K": [0, 1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17], 
            "H": [ 2,  3,  4,  5, 16, 17, 18, 19]
        },
        "2M0036_1103":{
            "K": [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], # [2,3,4,5,8,10,11,12,13,14,15,18], 
            "H": [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]#[0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19]
        },
        "2M0036_1105":{
            "K": [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], #[2,3,4,5,8,10,11,12,13,14,15,18], 
            "H": [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]#[1,2,3,4,5,9,10,11,16], #[0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19]
        },      
    }, # [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],
    "CRIRES": {
        "W1049B":{
            "K": [0, 1, 2, 3]
        },
        "W1049A":{
            "K": [0, 1, 2, 3]
        }
    }
}

nobss =   {"W1049B": 14,     "W1049A": 14,     "W1049B_0209": 14,     "W1049A_0209": 14,     "2M0036_1103": 7   ,  "2M0036_1105": 8}

periods = {"W1049B": 5,      "W1049A": 7,      "W1049B_0209": 5,      "W1049A_0209": 7,      "2M0036_1103": 2.7 ,  "2M0036_1105": 2.7}
incs =    {"W1049B": 80,     "W1049A": 70,     "W1049B_0209": 80,     "W1049A_0209": 70,     "2M0036_1103": 51  ,  "2M0036_1105": 51}
vsinis =  {"W1049B": 29e3,   "W1049A": 21e3,   "W1049B_0209": 29e3,   "W1049A_0209": 21e3,   "2M0036_1103": 33e3,  "2M0036_1105": 33e3}
rvs =     {"W1049B": 7.05e-5,"W1049A": 5.4e-5, "W1049B_0209": 7.05e-5,"W1049A_0209": 5.4e-5, "2M0036_1103": 6.5e-5,"2M0036_1105": 6.5e-5}
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
    "W1049B_0209": np.array(
        [0.15093   , 0.51746401, 0.87756001, 1.23739201, 1.59770401,
        1.97218201, 2.341842  , 2.70773401, 3.07090801, 3.435258  ,
        3.80697001, 4.168056  , 4.52813401, 4.88789401]),
    "W1049A_0209": np.array(
        [0.15093   , 0.51746401, 0.87756001, 1.23739201, 1.59770401,
        1.97218201, 2.341842  , 2.70773401, 3.07090801, 3.435258  ,
        3.80697001, 4.168056  , 4.52813401, 4.88789401]),
    "2M0036_1103": np.array(
        [0.        , 0.35493599, 0.709224  , 1.063536  , 2.07012   ,
        2.424528  , 2.779728  ]),
    "2M0036_1105": np.array(
        [0.        , 0.35536801, 0.727008  , 1.08223201, 1.43724001,
       1.79220001, 2.14716   , 2.50233601]),
}


##############################################################################
####################    Settings    ##########################################
##############################################################################

simulation_on = False
savedir = "igrinsH"
use_toy_spec = False
use_eqarea = True

#################### Run settings ####################################

flux_err = 0.01 if use_toy_spec else 0.025
instru = "IGRINS"
target = "W1049B"
band = "H"
#solver = "IC14new"
map_type = "eqarea"

modelspec = "t1500g1000f8"
LSD = "new"

########## IC14 parameters ##########
nk = 125
cut = nk - 70
LLD = 0.4
alpha = 5000
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
    contrast = "real"
    noisetype = "real"

    if instru == "CRIRES":
        nk = 203
        cut = nk - 70

    nobs = nobss[target]

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
