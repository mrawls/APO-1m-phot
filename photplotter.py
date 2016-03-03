from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Meredith Rawls, 2015

After running PHOT in IRAF, this program does differential photometry!

INPUT:
A TXDUMP-ed text file created from a set of IRAF PHOT outfiles.

You should run TXDUMP with these extracted fields on *.mag.1 (or *.mag.2, or whatever):
    --> IMAGE, XINIT, YINIT, OTIME, IFILTER, RAPERT, MAG, MERR

We assume there are FOUR apertures used, ONE target star, and SIX comparison stars.
If your needs differ, you'll need to manually edit the code. There are comments to help.

!! You must manually replace any INDEF entries in txdump_file with NAN !! (thanks numpy)

OUTPUT:
Intermediate plots of comparison star magnitudes (optional).
Plots of the target star's differential magnitude vs. time and phase.
(The comparison star instrumental magnitudes are also plotted if you span way down! fun!)
Outfile containing differential photometry data.
'''

# EDIT THIS STUFF AS DESIRED!!!
dir =               '../../1m_observations/KIC09291629/'
txdump_file =       'phot_take1.txt'
photcoord_file =    'photcoords1.txt'
outfile =           'BVRI_diffmag_LC1.txt'
#period = 33.659962; BJD0 = 54960.866328 # 3955867
#period = 235.30; BJD0 = 55190.482944 # 9970396
period = 20.686011; BJD0 = 54967.249343 # 9291629
aperturelist = [1,1,2,2,2,2,2] # which aperture to use? choose 1, 2, 3, or 4 FOR EACH STAR
compstars_good = [0,1,2,3    ] # which comparison stars are OK? at most [0,1,2,3,4,5]
compplot =      True    # set whether to plot comparison star LCs or not
phaseoffset =   False   # set whether to offset phase IN FINAL PLOT ONLY by 0.5 (like ELC)
diffmagdim =     0.0    # plot limits used at the very end
diffmagbright = -4.0
# EDIT THIS STUFF AS DESIRED!!!

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
        ##rmserr.append( np.sqrt(np.mean(mag_i - mag_ref)**2) ) # nope, but handy in other cases
        #mag_test.append(np.nanmean(magval)) # mean mag for one point in time
        #merr_test.append(np.nanmean(merrval)) # mean merr for one point in time
    magcomps = rms #mag_test
    merrcomps = rmserr #merr_test
    if plot == True:
        plt.ylabel('BACKWARDS MAGNITUDES ARE BACKWARDS') # seriously
        plt.xlabel('Time')
        plt.errorbar(otimes, magcomps, yerr=merrcomps, ls='None', marker='o', color='k', label='RMS')
        plt.title('filter ID: '+str(filter))
        plt.legend(numpoints=1)
        plt.show()
    return magcomps, merrcomps

def diffmagcalculate(mags, merrs, compmags, comperrs):
    '''
    Given magnitudes (and errors) for a target and a reference, subtract appropriately
    Adding errors in quadrature is a good life choice too
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

# Write results to file
f1 = open(dir+outfile, 'w')
print('# B filter', file=f1)
print('# UTCtime, phase, filter, mag, err, compmag, err, diffmag, err', file=f1)
for time, phase, mag, merr, comp, cerr, diff, derr in zip(otimeBs, phaseBs, magBs, merrBs, magBcomps, merrBcomps, diffBs, diffBerrs):
    print(time, phase, 'B', mag, merr, comp, cerr, diff, derr, file=f1)
print('# V filter', file=f1)
print('# UTCtime, phase, mag, err, compmag, err, diffmag, err', file=f1)
for time, phase, mag, merr, comp, cerr, diff, derr in zip(otimeVs, phaseVs, magVs, merrVs, magVcomps, merrVcomps, diffVs, diffVerrs):
    print(time, phase, 'V', mag, merr, comp, cerr, diff, derr, file=f1)
print('# R filter', file=f1)
print('# UTCtime, phase, mag, err, compmag, err, diffmag, err', file=f1)
for time, phase, mag, merr, comp, cerr, diff, derr in zip(otimeRs, phaseRs, magRs, merrRs, magRcomps, merrRcomps, diffRs, diffRerrs):
    print(time, phase, 'R', mag, merr, comp, cerr, diff, derr, file=f1)
print('# I filter', file=f1)
print('# UTCtime, phase, mag, err, compmag, err, diffmag, err', file=f1)
for time, phase, mag, merr, comp, cerr, diff, derr in zip(otimeIs, phaseIs, magIs, merrIs, magIcomps, merrIcomps, diffIs, diffIerrs):
    print(time, phase, 'I', mag, merr, comp, cerr, diff, derr, file=f1)
print('Data written to {0}'.format(outfile))
f1.close()

# Plot magnitude vs. orbital phase for all four filters
axtop = plt.subplot(2,1,1)
axtop.set_ylim([diffmagdim, diffmagbright])
plt.xlabel('Time (JD-2400000)')
plt.ylabel('Differental Mag')
plt.errorbar(otimeBs, diffBs, yerr=diffBerrs, ls='None', marker='o', color='b', label='B')
plt.errorbar(otimeVs, diffVs, yerr=diffVerrs, ls='None', marker='o', color='g', label='V')
plt.errorbar(otimeRs, diffRs, yerr=diffRerrs, ls='None', marker='o', color='r', label='R')
plt.errorbar(otimeIs, diffIs, yerr=diffIerrs, ls='None', marker='o', color='k', label='I')

plt.errorbar(otimeBs, magBcomps, yerr=merrBcomps, ls='None', marker='o', color='c', label='Bcomp')
plt.errorbar(otimeVs, magVcomps, yerr=merrVcomps, ls='None', marker='o', color='y', label='Vcomp')
plt.errorbar(otimeRs, magRcomps, yerr=merrRcomps, ls='None', marker='o', color='m', label='Rcomp')
plt.errorbar(otimeIs, magIcomps, yerr=merrIcomps, ls='None', marker='o', color='0.75', label='Icomp')

axbot = plt.subplot(2,1,2)
axbot.set_ylim([diffmagdim, diffmagbright])
#plt.axis([0, 1, diffmagdim, diffmagbright])
plt.xlabel('Orbital Phase')
plt.ylabel('Differental Mag')
if phaseoffset == True: # offset phases so primary eclipse appears at 0.5 instead of 0 (1)
    phaseBs = [phase+0.5 if phase < 0.5 else phase-0.5 for phase in phaseBs]
    phaseVs = [phase+0.5 if phase < 0.5 else phase-0.5 for phase in phaseVs]
    phaseRs = [phase+0.5 if phase < 0.5 else phase-0.5 for phase in phaseRs]
    phaseIs = [phase+0.5 if phase < 0.5 else phase-0.5 for phase in phaseIs]
plt.errorbar(phaseBs, diffBs, yerr=diffBerrs, ls='None', marker='o', color='b', label='B')
plt.errorbar(phaseVs, diffVs, yerr=diffVerrs, ls='None', marker='o', color='g', label='V')
plt.errorbar(phaseRs, diffRs, yerr=diffRerrs, ls='None', marker='o', color='r', label='R')
plt.errorbar(phaseIs, diffIs, yerr=diffIerrs, ls='None', marker='o', color='k', label='I')

#plt.plot(phaseBs, magBs, ls='None', marker='o', color='b')
#plt.plot(phaseVs, magVs, ls='None', marker='o', color='g')
#plt.plot(phaseRs, magRs, ls='None', marker='o', color='r')
#plt.plot(phaseIs, magIs, ls='None', marker='o', color='k')

plt.errorbar(phaseBs, magBcomps, yerr=merrBcomps, ls='None', marker='o', color='c', label='Bcomp')
plt.errorbar(phaseVs, magVcomps, yerr=merrVcomps, ls='None', marker='o', color='y', label='Vcomp')
plt.errorbar(phaseRs, magRcomps, yerr=merrRcomps, ls='None', marker='o', color='m', label='Rcomp')
plt.errorbar(phaseIs, magIcomps, yerr=merrIcomps, ls='None', marker='o', color='0.75', label='Icomp')

#plt.legend()
plt.show()
        
