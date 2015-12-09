from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Meredith Rawls, Dec 2015

Take some TXDUMP-ed files from PHOT, with these extracted fields:
    --> IMAGE, XINIT, YINIT, OTIME, IFILTER, RAPERT, MAG, MERR
(For now, assume there are four apertures all printed out together.)
(Manually replace any INDEF entries with NAN.)

TODO: write out results to file!

Makes plots of magnitude vs. phase.

'''
dir = '../../1m_observations/KIC03955867/'
txdump_file = 'phot_take8.txt'
photcoord_file = 'photcoords5.txt'
period = 33.659962; BJD0 = 54960.866328 # 3955867
aperturelist = [1,1,3,3,2,4,4] # which aperture to use? set 1, 2, 3, or 4 FOR EACH STAR
compstars_good = [0,1,2,3,4  ] # which comparison stars are OK? at most [0,1,2,3,4,5]
compplot = True # set whether to plot comparison star LCs or not

# read in the coordinate file you used to do photometry
xinits, yinits = np.loadtxt(dir+photcoord_file, usecols=(0,1), unpack=True)

def phasecalc(time, period, zeropoint):
    '''
    Calculates orbital phase of a single observation given period and zeropoint
    '''
    fracP = (time - zeropoint) / period
    phase = fracP % 1
    return phase

def readfilter(filtID = 1, aperturelist = [1,1,1,1,1,1,1]): #filtID, xcoords, obstimes, ifilters, mags, merrs):
    '''
    Returns everything you need to process one filter's worth of info
    Assumes one target star and six comparison stars
    It may not be the most efficient way to do this, but it works, and isn't super slow
    '''
    allotimelist = []; allphaselist = []; allmaglist = []; allmerrlist = []
    for star in range(0, 7): # assuming 1 target star + 6 comparison stars
        # use aperturelist to decide which magnitudes to read in
        # given 4 apertures, the last two values should be a pair of 9,10,11,12 & 13,14,15,16
        if aperturelist[star] == 1: magidx = 9; merridx = 13
        elif aperturelist[star] == 2: magidx = 10; merridx = 14
        elif aperturelist[star] == 3: magidx = 11; merridx = 15
        elif aperturelist[star] == 4: magidx = 12; merridx = 16
        else: print('No aperture index found for star {0}!'.format(star))
        # Read data in from the photometry txdump_file
        xcoords, obstimes, ifilters, mags, merrs = np.loadtxt(dir+txdump_file, 
            comments='#', usecols=(1,3,4,magidx,merridx), unpack=True)
        # loop through the read-in data and save it by filter and star
        otimelist = []; phaselist = []; maglist = []; merrlist = []
        for x, otime, ifilt, mag, merr in zip(xcoords, obstimes, ifilters, mags, merrs):
            if ifilt == filtID and x == xinits[star]: # images for one star in current filter
                otimelist.append(otime) # only need to save this once
                phaselist.append(phasecalc(otime, period, BJD0)) # only need to save this once
                maglist.append(mag) # need to save this each time
                merrlist.append(merr) # need to save this each time
        allmaglist.append(np.array(maglist))
        allmerrlist.append(np.array(merrlist))
    return otimelist, phaselist, allmaglist, allmerrlist 

def compstarcombine(filter, otimes, maglist, merrlist, plot = False):
    '''
    Given a set of comparison star magnitudes and errors, combine them
    --> maglist and merrlist are each a list of N lists, one per comparison star
    --> magcomps and merrcomps are each a single list of comparison magnitude info

    There's an option to plot the comparison star light curves, or not.
    '''
    mag_test = []; merr_test = []; rms = []; rmserr = []
    nstars = len(maglist) 
    print('Filter ID: {0}'.format(filter))
    for star_idx in compstars_good: # use to select decent comparison stars
        print('comp star', star_idx, np.nanmean(maglist[star_idx]), np.nanstd(maglist[star_idx]))
        if plot == True:
            plt.plot(otimes, maglist[star_idx], ls='None', marker='.', label='Comp Star '+str(star_idx))
    for time_idx in range(0, len(maglist[0])): # loop over each image
        magval = []; merrval = []
        for star_idx in compstars_good: # loop over each star you want to use
            magval.append(maglist[star_idx][time_idx]) # build list of all star mags in one image
            merrval.append(merrlist[star_idx][time_idx])

        rms.append( np.sqrt(np.nanmean(np.power(magval,2))) ) # RMS mag for one point in time
        rmserr.append( np.sqrt(np.nanmean(np.power(merrval,2))) )
        #rmserr.append( np.sqrt(np.mean(mag_i - mag_ref)**2) ) # nope

        #mag_test.append(np.nanmean(magval)) # mean mag for one point in time
        #merr_test.append(np.nanmean(merrval))
    magcomps = rms #mag_test
    merrcomps = rmserr #merr_test
    if plot == True:
        plt.errorbar(otimes, magcomps, yerr=merrcomps, ls='None', marker='o', color='k', label='RMS')
        plt.title('filter ID: '+str(filter))
        plt.legend(numpoints=1)
        plt.show()
    return magcomps, merrcomps

def diffmagcalculate(mags, merrs, compmags, comperrs):
    '''
    Given magnitudes (and errors) for a target and a reference, subtract appropriately
    '''
    diffmags = mags - compmags
    differrs = np.sqrt(np.power(merrs,2) + np.power(comperrs,2))
    return diffmags, differrs

# given four apertures, the last two values should be a pair of 9,10,11,12 & 13,14,15,16
#if aperture_idx == 1:
#    magidx = 9; merridx = 13
#elif aperture_idx == 2:
#    magidx = 10; merridx = 14
#elif aperture_idx == 3:
#    magidx = 11; merridx = 15
#elif aperture_idx == 4:
#    magidx = 12; merridx = 16
#else:
#    print('No aperture index found!')

# Read data in from the photometry txdump_file
#xs, otimes, ifilters, mags, merrs = np.loadtxt(dir+txdump_file, 
#    comments='#', usecols=(1,3,4,magidx,merridx), unpack=True)

### B filter ###
Binfo = readfilter(filtID = 2, aperturelist = aperturelist)
otimeBs = Binfo[0]; phaseBs = Binfo[1]
magBs = Binfo[2][0]; merrBs = Binfo[3][0]
magBcomplist = []; merrBcomplist = []
magBcomplist.extend([Binfo[2][1], Binfo[2][2], Binfo[2][3], Binfo[2][4], Binfo[2][5], Binfo[2][6]])
merrBcomplist.extend([Binfo[3][1], Binfo[3][2], Binfo[3][3], Binfo[3][4], Binfo[3][5], Binfo[3][6]])
magBcomps, merrBcomps = compstarcombine('B', otimeBs, magBcomplist, merrBcomplist, plot = compplot)
diffBs, diffBerrs = diffmagcalculate(magBs, merrBs, magBcomps, merrBcomps)

### V filter ###
Vinfo = readfilter(filtID = 3, aperturelist = aperturelist)
otimeVs = Vinfo[0]; phaseVs = Vinfo[1]
magVs = Vinfo[2][0]; merrVs = Vinfo[3][0]
magVcomplist = []; merrVcomplist = []
magVcomplist.extend([Vinfo[2][1], Vinfo[2][2], Vinfo[2][3], Vinfo[2][4], Vinfo[2][5], Vinfo[2][6]])
merrVcomplist.extend([Vinfo[3][1], Vinfo[3][2], Vinfo[3][3], Vinfo[3][4], Vinfo[3][5], Vinfo[3][6]])
magVcomps, merrVcomps = compstarcombine('V', otimeVs, magVcomplist, merrVcomplist, plot = compplot)
diffVs, diffVerrs = diffmagcalculate(magVs, merrVs, magVcomps, merrVcomps)

### R filter ###
Rinfo = readfilter(filtID = 4, aperturelist = aperturelist)
otimeRs = Rinfo[0]; phaseRs = Rinfo[1]
magRs = Rinfo[2][0]; merrRs = Rinfo[3][0]
magRcomplist = []; merrRcomplist = []
magRcomplist.extend([Rinfo[2][1], Rinfo[2][2], Rinfo[2][3], Rinfo[2][4], Rinfo[2][5], Rinfo[2][6]])
merrRcomplist.extend([Rinfo[3][1], Rinfo[3][2], Rinfo[3][3], Rinfo[3][4], Rinfo[3][5], Rinfo[3][6]])
magRcomps, merrRcomps = compstarcombine('R', otimeRs, magRcomplist, merrRcomplist, plot = compplot)
diffRs, diffRerrs = diffmagcalculate(magRs, merrRs, magRcomps, merrRcomps)

### I filter ###
Iinfo = readfilter(filtID = 6, aperturelist = aperturelist)
otimeIs = Iinfo[0]; phaseIs = Iinfo[1]
magIs = Iinfo[2][0]; merrIs = Iinfo[3][0]
magIcomplist = []; merrIcomplist = []
magIcomplist.extend([Iinfo[2][1], Iinfo[2][2], Iinfo[2][3], Iinfo[2][4], Iinfo[2][5], Iinfo[2][6]])
merrIcomplist.extend([Iinfo[3][1], Iinfo[3][2], Iinfo[3][3], Iinfo[3][4], Iinfo[3][5], Iinfo[3][6]])
magIcomps, merrIcomps = compstarcombine('I', otimeIs, magIcomplist, merrIcomplist, plot = compplot)
diffIs, diffIerrs = diffmagcalculate(magIs, merrIs, magIcomps, merrIcomps)

# Plot magnitude vs. orbital phase for all four filters
plt.axis([0, 1, 1, -2.5])
plt.xlabel('Orbital Phase')
plt.ylabel('Differental Mag')
plt.errorbar(phaseBs, diffBs, yerr=diffBerrs, ls='None', marker='o', color='b', label='B')
plt.errorbar(phaseVs, diffVs, yerr=diffVerrs, ls='None', marker='o', color='g', label='V')
plt.errorbar(phaseRs, diffRs, yerr=diffRerrs, ls='None', marker='o', color='r', label='R')
plt.errorbar(phaseIs, diffIs, yerr=diffIerrs, ls='None', marker='o', color='k', label='I')

plt.plot(phaseBs, magBs, ls='None', marker='o', color='b')
plt.plot(phaseVs, magVs, ls='None', marker='o', color='g')
plt.plot(phaseRs, magRs, ls='None', marker='o', color='r')
plt.plot(phaseIs, magIs, ls='None', marker='o', color='k')

plt.plot(phaseBs, magBcomps, ls='None', marker='o', color='c', label='Bcomp')
plt.plot(phaseVs, magVcomps, ls='None', marker='o', color='y', label='Vcomp')
plt.plot(phaseRs, magRcomps, ls='None', marker='o', color='m', label='Rcomp')
plt.plot(phaseIs, magIcomps, ls='None', marker='o', color='0.75', label='Icomp')

#plt.legend()
plt.show()
        
