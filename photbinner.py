from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Meredith Rawls, 2015
Takes a file created by photplotter.py and bins it appropriately.
'''
dir =       '../../1m_observations/KIC03955867/'
infile =    'BVRI_diffmag_LC.txt'
target =    '3955867'
period = 33.659962; BJD0 = 54960.866328 # 3955867
time_threshold = 0.1 # time in days by which to bin observations
filtlist = ['B', 'V', 'R', 'I']

def phasecalc(time, period, zeropoint):
    '''
    Calculates orbital phase of a single observation given period and zeropoint
    '''
    fracP = (time - zeropoint) / period
    phase = fracP % 1
    return phase

times, phases, filters, mags, merrs, compmags, cerrs, diffmags, derrs = \
    np.loadtxt(dir+infile, comments='#', unpack=True, dtype={ \
    'names': ('times', 'phases', 'filters', 'mags', 'merrs', 'compmags', 'cerrs', \
    'diffmags', 'derrs'), 'formats': (np.float64, np.float64, '|S2', np.float64, \
    np.float64, np.float64, np.float64, np.float64, np.float64)})

for filt in filtlist:
    # initialize chunk lists
    chunk = 0
    timechunk = [[]]; phasechunk = [[]]
    magchunk = [[]]; merrchunk = [[]]
    compmagchunk = [[]]; cerrchunk = [[]]
    diffmagchunk = [[]]; derrchunk = [[]]
    # load observations for one filter
    for idx, (time, phase, filter, mag, merr, compmag, cerr, diffmag, derr) in enumerate( \
        zip(times, phases, filters, mags, merrs, compmags, cerrs, diffmags, derrs)):
        if filter == filt:
            # calculate time elapsed since previous image
            if idx > 0: 
                dt = (time - times[idx-1])
            else: 
                dt = 0  
            # start a new chunk if dt > some threshold
            if dt >= time_threshold: 
                chunk += 1
                timechunk.append([]); phasechunk.append([])
                magchunk.append([]); merrchunk.append([])
                compmagchunk.append([]); cerrchunk.append([])
                diffmagchunk.append([]); derrchunk.append([])
            # save observations by chunk
            timechunk[chunk].append(time); phasechunk[chunk].append(phase)
            magchunk[chunk].append(mag); merrchunk[chunk].append(merr)
            compmagchunk[chunk].append(compmag); cerrchunk[chunk].append(cerr)
            diffmagchunk[chunk].append(diffmag); derrchunk[chunk].append(derr)
            nchunks = len(timechunk)
    newtime = []; newphase = []
    newmag = []; newmerr = []
    newcomp = []; newcerr = []
    newdiff = []; newderr = []
    for cdx, time in enumerate(timechunk):
        # cdx is chunk index
        # calculate rms for magnitude data points, one per chunk
        newtime.append( (timechunk[cdx][0] + timechunk[cdx][-1]) / 2. )
        newphase.append( phasecalc((timechunk[cdx][0] + timechunk[cdx][-1]) / 2., period, BJD0) )
        newmag.append( np.sqrt(np.nanmean(np.power(magchunk[cdx],2))) ) # RMS mag for one chunk
        newmerr.append( np.sqrt(np.nanmean(np.power(merrchunk[cdx],2))) )
        newcomp.append( np.sqrt(np.nanmean(np.power(compmagchunk[cdx],2))) ) # RMS compmag for one chunk
        newcerr.append( np.sqrt(np.nanmean(np.power(cerrchunk[cdx],2))) )
        newdiff.append( np.sqrt(np.nanmean(np.power(diffmagchunk[cdx],2))) ) # RMS diffmag for one chunk
        newderr.append( np.sqrt(np.nanmean(np.power(derrchunk[cdx],2))) )
        #print(newtime, newphase, newdiff, newderr) #newmag, newmerr, newcomp, newcerr, newdiff, newderr)
        # SAVE THESE VALUES! PLOT THEM! CHECK IF THEY ARE CORRECT??!!
    plt.errorbar(newphase, newdiff, yerr=newderr, ls='None', marker='o')
    plt.show()