"""
# Original author: Ian Crossfield (Python 2.7)

Handy Routines for Doppler Imaging and Maximum-Entropy.

To start up with data to fit of size N pixels, and M pixels per
observation, do:

import dime
dime.setup(N, M)
"""

######################################################

# 18-02-2020
# Emma Bubb - change to run on Python 3

######################################################

# EB update the imported toolkits
import numpy as np


nk, j, k, l, jfrac = 0,0,0,0,0

def setup(ndat, nk0):
    global nk, j, k, l, jfrac
    nk = nk0
    j = np.arange(ndat)
    k = (np.floor(1.0*j/nk)*nk).astype(int)
    l = (np.ceil((j+1.0)/nk)*nk - 1).astype(int)
    jfrac = (j % nk) / (nk - 1.0)
    return


# The SciPy way:
def nentropy(x):
  """ Compute Normalized Entropy, Sum(y * ln(y)), where y_i = x_i/sum(x)"""
  # 2013-08-07 21:18 IJMC: Created
  norm_x = x / np.sum(x)
  entropy = -np.sum(norm_x * np.log(norm_x))
  return entropy

def dnentropy_dx(x):
  xsum = 1.0*np.sum(x)
  norm_x = x / xsum
  nx = len(x)
  vec2 = np.log(norm_x) + 1.
  #vec1s_old = (np.tile(x, (nx,1)) + xsum*np.diag(np.ones(nx)))
  vec1s = (-np.tile(x, (nx,1)) + xsum*np.diag(np.ones(nx)))
  grads = -np.dot(vec1s, vec2)/xsum/xsum
  return  grads

def dchisq_dx(x, model):
  nx = len(x)
  dif = (w_observation * (observation_cnorm[0:nk*nobs] - model))
  grads = -2*(dif*Rmatrix).sum(1)
  return grads

def gnorm(unnorm_model, nk):
  """ Compute the normalizing function"""
  return unnorm_model[k] + (unnorm_model[l] - unnorm_model[k]) * jfrac


def dgnorm_dx(Rmatrix, nk):
  return Rmatrix[:,k] + jfrac * (Rmatrix[:,l] - Rmatrix[:,k])

def normalize_model(unnorm_model, nk):
  return unnorm_model / gnorm(unnorm_model, nk)


def dchisq_norm_dx(unnorm_model, nk, data, weights, Rmatrix):
  normalizer = gnorm(unnorm_model, nk)
  dif = (weights * (data - unnorm_model / normalizer))
  Rk = Rmatrix[:,k]
  dfdx = (normalizer * Rmatrix - unnorm_model * (Rk + jfrac * (Rmatrix[:,l] - Rk))) / normalizer / normalizer
  grads = -2*(dif*dfdx).sum(1)
  return grads

def getgrad_norm_sp(x, *args):
  data, weights, R, alpha = args[0:4]
  model = np.dot(x.ravel(), R)
  ds = dnentropy_dx(x)
  dchi = dchisq_norm_dx(model, nk, data, weights, R) 
  return 0.5 * dchi - alpha*ds

def getgrad_sp(x, *args): # No model-normalization
  data, weights, R, alpha = args[0:4]
  model = np.dot(x.ravel(), R)
  ds = dnentropy_dx(x)
  dchi = dchisq_dx(x, model) 
  return 0.5 * dchi - alpha*ds

def entropy_map_sp(map_pixels, *args):
  data, weights, R, alpha = args[0:4]
  if (map_pixels<=0).any():
    map_pixels[map_pixels<=0] = 1e-6 
  #print( 'd', data[0], map_pixels.min() ) # EB updated print statement
  entropy = nentropy(map_pixels) 
  model = np.dot(map_pixels.ravel(), R)
  #dchi = dchisq_dx(x, model) 
  norm_model = model #normalize_model(model, nk, nobs)
  chisq = (weights*(norm_model-data)**2).sum()
  #ds = dnentropy_dx(x)
  #print( 'e', entropy, chisq) # EB updated print statement
  return 0.5*chisq - alpha*entropy #, 0.5*dchi - alpha*ds

def entropy_map_norm_sp(map_pixels, *args, **kw):
  data, weights, R, alpha = args[0:4]
  if (map_pixels<=0).any():
    map_pixels[map_pixels<=0] = 1e-6 #EB: if any pixel values are negative, set to 1e-6 (vv small basically zero)
  #print( 'd', data[0], map_pixels.min() ) # EB updated print statement
  entropy = nentropy(map_pixels) #EB: call function 'nentropy' to calculate the mormalised entropy
  model = np.dot(map_pixels.ravel(), R)
  norm_model = normalize_model(model, nk) #call function 'normalize_model' to normalise the model (basically model/normalising function)
  chisq = (weights*(norm_model-data)**2).sum() # EB: changed method of finding chi squared from calling function to calculating directly

  ret = 0.5*chisq - alpha*entropy # EB add to get rid of error

  for key in kw: #EB update
  #if kw.has_key('retvals') and kw['retvals']==True: #EB - This method does not work in python3
      if key=='retvals': #EB update
        if kw[key]==True:
          ret = 0.5*chisq - alpha*entropy, chisq, entropy
      else:
          ret = 0.5*chisq - alpha*entropy
  return ret #EB: return the entropy 

def profile_spotmap(param, *args, **kw):
    """MOdel line profiles, assuming a simple one-spot model.

    phi, theta, R = args[0:3]
    startemp, spottemp, spotlat, spotlon, spotrad = param[0:5]
     OR
    startemp, temp1, lat1, lon1, rad1, temp2, lat2, lon2, rad2 = param[0:9]
    """
    # 2013-08-19 09:59 IJMC: Created
    # 2013-08-27 10:45 IJMC: Updated to multi-spot-capable
    from maps3 import makespot

    phi, theta, R = args[0:3]
    nparam = len(param)
    nspots = int((nparam-1)/4)
    startemp = param[0]
    map_pixels = np.ones(phi.shape) * param[0]
    for ii in range(nspots):
        spottemp, spotlat, spotlon, spotrad = param[1+ii*4:1+(ii+1)*4]
        boolspot = makespot(spotlat, spotlon, spotrad, phi, theta).astype(np.float32)
        map_pixels -= boolspot * (startemp - spottemp)

    return normalize_model(np.dot(map_pixels.ravel(), R), nk)




