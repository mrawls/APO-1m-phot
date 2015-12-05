from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Meredith Rawls, Dec 2015

Take some TXDUMP-ed files from PHOT, with these extracted fields:
    --> IMAGE, XINIT, YINIT, OTIME, IFILTER, RAPERT, MAG, MERR
(For now, assume there are four apertures all printed out together.)
(Manually replace any INDEF entries with NAN.)

TODO: actual differential photometry (use mag data from other stars!)
TODO: try other apertures!

Make plots of magnitude vs. time! Or maybe even phase!

'''
#napertures = 4 # feature not implemented yet
dir = '../../1m_photometry/KIC03955867/'
txdump_file = 'phot_take2.txt'
period = 33.659962; BJD0 = 54960.866328 # 3955867

# hard wired for now, from the coordinate file you use to do photometry
xinits = [1080.190, 1081.000, 1057.570, 1483.600, 1605.570, 1265.290, 671.490]
yinits = [990.060, 950.250, 534.550, 495.680, 1271.130, 1530.380, 1434.400]

def phasecalc(time, period, zeropoint):
    '''
    Calculates orbital phase of a single observation given period and zeropoint
    '''
    fracP = (time - zeropoint) / period
    phase = fracP % 1
    return(phase)

# adjust 'usecols' as appropriate for your aperture selection
xs, otimes, ifilters, mags, merrs = np.loadtxt(dir+txdump_file, 
    comments='#', usecols=(1,3,4,9,13), unpack=True)

otimeBlist = []
phaseBlist = []
magBlist = []
merrBlist = []
for x, otime, ifilt, mag, merr in zip(xs, otimes, ifilters, mags, merrs):
    if ifilt == 2 and x == xinits[0]: # B images of star 1 only
        otimeBlist.append(otime)
        phaseBlist.append(phasecalc(otime, period, BJD0))
        magBlist.append(mag)
        merrBlist.append(merr)

otimeVlist = []
phaseVlist = []
magVlist = []
merrVlist = []
for x, otime, ifilt, mag, merr in zip(xs, otimes, ifilters, mags, merrs):
    if ifilt == 3 and x == xinits[0]: # V images of star 1 only
        otimeVlist.append(otime)
        phaseVlist.append(phasecalc(otime, period, BJD0))
        magVlist.append(mag)
        merrVlist.append(merr)

otimeRlist = []
phaseRlist = []
magRlist = []
merrRlist = []
for x, otime, ifilt, mag, merr in zip(xs, otimes, ifilters, mags, merrs):
    if ifilt == 4 and x == xinits[0]: # R images of star 1 only
        otimeRlist.append(otime)
        phaseRlist.append(phasecalc(otime, period, BJD0))
        magRlist.append(mag)
        merrRlist.append(merr)

otimeIlist = []
phaseIlist = []
magIlist = []
merrIlist = []
for x, otime, ifilt, mag, merr in zip(xs, otimes, ifilters, mags, merrs):
    if ifilt == 6 and x == xinits[0]: # I images of star 1 only
        otimeIlist.append(otime)
        phaseIlist.append(phasecalc(otime, period, BJD0))
        magIlist.append(mag)
        merrIlist.append(merr)

#plt.axis([56700, 57400, 20, 16])
#plt.plot(otimeBlist, magBlist, ls='None', marker='o', color='b')
#plt.plot(otimeVlist, magVlist, ls='None', marker='o', color='g')
#plt.plot(otimeRlist, magRlist, ls='None', marker='o', color='r')
#plt.plot(otimeIlist, magIlist, ls='None', marker='o', color='k')

plt.axis([0, 1, 20, 16])
plt.plot(phaseBlist, magBlist, ls='None', marker='o', color='b')
plt.plot(phaseVlist, magVlist, ls='None', marker='o', color='g')
plt.plot(phaseRlist, magRlist, ls='None', marker='o', color='r')
plt.plot(phaseIlist, magIlist, ls='None', marker='o', color='k')

plt.show()
        
