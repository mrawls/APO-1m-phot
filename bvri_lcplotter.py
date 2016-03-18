import numpy as np
import matplotlib.pyplot as plt
'''
Makes a multi-pane figure of BVRI light curves for six EBs.
'''
sysnames =  ['08702921',    '09291629', '03955867', '10001167', '05786154', '09970396']
periods =   [19.38446,      20.6864,    33.65685,   120.3903,   197.9182,   235.29852]
zeros =     [54970.2139,    54966.891,  54960.8989, 54957.682,  55162.6140, 55190.5400]
seclocs =   [0.4425,        0.4987,     0.4973,     0.5850,     0.2834,     0.4133]
filenames = []
dir = '../../1m_observations/KIC'
for sysname in sysnames:
    filenames.append(dir+sysname+'/BVRI_diffmag_binned_LC1.txt')

#kdir = '../../RG_ELCmodeling/'
#kfiles = [kdir+'8702921/trial3/modelU.mag', kdir+'9291629/trial5/modelU.mag', 
#          kdir+'3955867/trial5/modelU.mag', kdir+'10001167/trial1/modelU.mag',
#          kdir+'5786154/trial3/modelU.mag', kdir+'9970396/trial7v2/modelU.mag']

kdir = '../../RG_light_curves/KIC_'
ktail = '_LC_mag_var_ecl_Q017.txt'
kfiles = [kdir+'8702921'+ktail, kdir+'9291629'+ktail, kdir+'3955867'+ktail,
         kdir+'10001167'+ktail, kdir+'5786154'+ktail, kdir+'9970396'+ktail]
offsets = [0.10, 0.0, -0.1, -0.15, 0.2, -0.25] # custom magnitude offset from zero

def phasecalc(time, period, zeropoint):
    '''
    Calculates orbital phase of a single observation given period and zeropoint
    '''
    fracP = (time - zeropoint) / period
    phase = fracP % 1
    return phase

# set up figure
fig = plt.figure(1, figsize=(15,8))
fig.text(0.4, 0.02, 'Orbital Phase', ha='center', va='center', size=26)
fig.text(0.07, 0.5, 'Differential Magnitude', ha='center', va='center', size=26, rotation='vertical')

fig.text(0.133, 0.1155, 'KIC 5786154')
fig.text(0.133, 0.385, 'KIC 3955867')
fig.text(0.133, 0.6545, 'KIC 8702921')
fig.text(0.548, 0.1155, 'KIC 9970396')
fig.text(0.548, 0.385, 'KIC 10001167')
fig.text(0.548, 0.6545, 'KIC 9291629')

# manual legend
fig.text(0.7, 0.02, '$B$')
fig.text(0.75, 0.02, '$V$')
fig.text(0.8, 0.02, '$R$')
fig.text(0.85, 0.02, '$I$')
ax1 = fig.add_axes([0.5, 0.01, 0.5, 0.05])
ax1.axis([0.5, 1.0, 0.01, 0.05])
ax1.axis('off')
ax1.errorbar(0.69, 0.025, yerr=0.012, marker='o', ms=10, color='#377eb8')
ax1.errorbar(0.74, 0.025, yerr=0.012, marker='o', ms=10, color='#4daf4a')
ax1.errorbar(0.79, 0.025, yerr=0.012, marker='o', ms=10, color='#e41a1c')
ax1.errorbar(0.84, 0.025, yerr=0.012, marker='o', ms=10, color='#984ea3')

# loop over stars
for idx, (system, file, kfile, period, zero) in enumerate(zip(sysnames, filenames, kfiles, periods, zeros)):
    # create subplot axis
    ax = fig.add_subplot(3, 2, idx+1)
    # plot customization
    plt.subplots_adjust(hspace=0, wspace=0.15)
    if idx != 4 and idx != 5:
        ax.set_xticklabels([])
#    if system == '08702921':
#        ax.set_yticks([-0.2, 0.0, 0.2, 0.4])
#    if system == '03955867':
#        ax.set_yticks([-0.4, -0.2, 0.0, 0.2])
#    if system == '09970396':
#        ax.set_yticks([-0.6, -0.4, -0.2, 0.0, 0.2, 0.4])
    
    # read ELC Kepler model LC data for current star
    print('Reading ELC Kepler model LC data for {0}...'.format(system))
    ktimes, kmags = np.loadtxt(kfile, comments='#', usecols=(0,1), unpack=True)
    # adjust magnitudes
    kmags = kmags - np.median(kmags) + offsets[idx]
    kphases = []
    print('Folding...')
    for ktime in ktimes:
        kphases.append( phasecalc(ktime, period, zero-54833.) )
    newkphases = []
    # adjust kphase to go from 0.2 to 1.2 instead of 0 to 1
    for kphase in kphases:
        if kphase < 0.2:
            newkphase = kphase + 1.0
        else:
            newkphase = kphase
        newkphases.append(newkphase)
    # sort by phase and plot result in the background
    newkphases = np.array(newkphases)
    kmags = kmags[np.argsort(newkphases)]
    newkphases = newkphases[np.argsort(newkphases)]
    plt.plot(newkphases, kmags, color='0.75', marker='.', ls='None', markevery=100)
    
    # plot vertical lines at eclipse locations
    plt.axvline(x=1, color='k', ls=':')
    plt.axvline(x=seclocs[idx], color='k', ls=':')
    
    # read BVRI binned data for current star
    filters, times, phases, mags, merrs = \
        np.loadtxt(file, comments='#', unpack=True, \
        dtype={'names':   ('filters', 'times', 'phases', 'mags', 'merrs'), \
        'formats': ('|S2', np.float64, np.float64, np.float64, np.float64)}) \
    # adjust differential magnitudes
    mags = mags - np.average(mags)
    # set axis limits
    ax.set_xlim([0.2, 1.2])
    ax.set_ylim([np.max(mags)+0.1, np.min(mags)-0.1])
    #ax.text(0.05, -0.6*(np.abs(np.min(mags) - np.max(mags))), system, ha='center', va='center', size=22)
    # loop over observations and make a plot
    for time, phase, mag, merr, filter in zip(times, phases, mags, merrs, filters):
        if filter == 'B':
            plt.errorbar(phase, mag, yerr=merr, ls='None', marker='o', color='#377eb8', ms=10, label='B')
            plt.errorbar(phase+1, mag, yerr=merr, ls='None', marker='o', color='#377eb8', ms=10, label='B')
        elif filter == 'V':
            plt.errorbar(phase, mag, yerr=merr, ls='None', marker='o', color='#4daf4a', ms=10, label='V')
            plt.errorbar(phase+1, mag, yerr=merr, ls='None', marker='o', color='#4daf4a', ms=10, label='V')
        elif filter == 'R':
            plt.errorbar(phase, mag, yerr=merr, ls='None', marker='o', color='#e41a1c', ms=10, label='R')
            plt.errorbar(phase+1, mag, yerr=merr, ls='None', marker='o', color='#e41a1c', ms=10, label='R')
        elif filter == 'I':
            plt.errorbar(phase, mag, yerr=merr, ls='None', marker='o', color='#984ea3', ms=10, label='I')
            plt.errorbar(phase+1, mag, yerr=merr, ls='None', marker='o', color='#984ea3', ms=10, label='I')

    
plt.show()
