Enable desktop notifications for Gmail.
   OK  No thanks
1 of 1,033
Fw: inclination Monte Carlo code
Inbox

Xueqing Chen <xueqing.chen@ed.ac.uk>
Attachments
5:33 PM (0 minutes ago)
to me



From: Beth Biller <bb@roe.ac.uk>
Sent: 02 November 2023 11:56
To: Xueqing Chen <xueqing.chen@ed.ac.uk>
Subject: inclination Monte Carlo code
 
inclination Monte Carlo code

https://www.dropbox.com/scl/fi/uoltmq9av4979h33s55d0/WISE1049A_calc_inclination.py?rlkey=w75vrg51t45krijm2d606ieaz&dl=0

WISE1049A_calc_inclination.py
Shared with Dropbox
www.dropbox.com


https://www.dropbox.com/scl/fi/uoltmq9av4979h33s55d0/WISE1049A_calc_inclination.py?rlkey=w75vrg51t45krijm2d606ieaz&dl=0

WISE1049A_calc_inclination.py
Shared with Dropbox
www.dropbox.com


************************************************
Dr. Beth Biller 

Personal Chair of Exoplanet Characterisation

Institute for Astronomy
The University of Edinburgh
Royal Observatory
Blackford Hill
Edinburgh EH9 3HJ
U.K. 

Tel: +44 (0) 131 668 8349



My pronouns are she/her.  To find out more about pronouns 
and why I'm sharing mine, click here.
************************************************


The University of Edinburgh is a charitable body, registered in Scotland, with registration number SC005336. Is e buidheann carthannais a th’ ann an Oilthigh Dhùn Èideann, clàraichte an Alba, àireamh clàraidh SC005336.
...

[Message clipped]  View entire message
3
 Attachments
  •  Scanned by Gmail
import numpy as np
import pylab as plt
import IPython
import matplotlib.patches as patches

from astropy.io import fits
from astropy.table import Table

# use vsini results from Bryan et al. 2020
mu, sigma = 26.1, 0.2 # 2.8, mean and standard deviation
vsini_arr = np.random.normal(mu, sigma, 100000)

# generate Gaussian of radius values
mu, sigma = 1.00, 0.1 # 
radius_arr = np.random.normal(mu, sigma, 100000)

# generate Gaussian of period values
mu, sigma = 5.0, 0.2
period_arr = np.random.normal(mu, sigma, 100000)

Rjup = 69911.   # radius of jupiter in km
hourtosec = 3600.  # conversion from hours to seconds

# generate range of equatorial velocities from radius and period
varr = (radius_arr * 2.*3.14159*Rjup) / (period_arr * hourtosec)

# calculate array of sin i
sini_arr = vsini_arr / varr


count, bins, ignored = plt.hist(sini_arr, 30, density=True)

arrgt1 = (sini_arr > 1.)

arrlt1 = (sini_arr < 1.)

np.median(np.arcsin(sini_arr[arrlt1]))*(180./3.14159)
np.mean(np.arcsin(sini_arr[arrlt1]))*(180./3.14159)
np.std(np.arcsin(sini_arr[arrlt1]))*(180./3.14159)

#sini_arr[arrgt1] = 1.0


#count, bins, ignored = plt.hist(inc_arr, 30, normed=True)


import matplotlib.pyplot as pl
hfont = {'fontname':'Helvetica', 'size':'25'}
fig = plt.figure()
ax = fig.add_subplot(141)
ax.hist(varr, 30, density=True, color='blue')
plt.xlabel('equatorial velocity (km/s)', **hfont)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20)

ax = fig.add_subplot(142)
ax.hist(vsini_arr, 30, density=True, color='green')
plt.xlabel('v sin i (km/s)', **hfont)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20)

ax = fig.add_subplot(143)
ax.hist(sini_arr, 30, density=True, color='red')
plt.xlabel('sin i', **hfont)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20)
ax.set_xlim(0.5, 1.2)
ax.add_patch(
    patches.Rectangle(
        (1, 0),   # (x,y)
        (0.2),          # width
        6.,          # height
        color='gray', alpha=0.5, label='nonphysical')
)


plt.plot([1.0, 1.0], [0, 6], 'k')

sini_arr_pinned = sini_arr
arrgt1 = (sini_arr_pinned > 1.)
sini_arr_pinned[arrgt1] = 1.0

inc_arr = np.arcsin(sini_arr)
np.median(inc_arr)*(180./3.14159)
np.mean(inc_arr)*(180./3.14159)
np.std(inc_arr)*(180./3.14159)


ax = fig.add_subplot(144)
ax.hist(inc_arr*(180. / 3.14159), 30, density=True, color='magenta')
plt.xlabel('i (degrees)', **hfont)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20)
#ax.set_xlim(0, 90)
#ax.xaxis.label.set_size(40)