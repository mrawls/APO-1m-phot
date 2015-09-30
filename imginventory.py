from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time
from astropy.time import TimeDelta
import os
'''
Read in 1m observation metadata to figure out:
- which stars were imaged when
- whether they were in eclipse, or not
- some kind of image preview or quality flag (NOT IMPLEMENTED YET)
- make a plot of this info

To work, this program needs to be saved somewhere that can see imagedir and reffile.
It assumes there are date-formatted subdirectories in imagedir (e.g. 150212) with FITS
files saved in them. Some of those FITS files might be for targets we don't care about.
'''
imagedir = '/mnt/mrawls/1m_obs/'
reffile = 'RGEB_info_alpha.txt'
outfile = 'imginventory_list.txt'

# Get the paths to the directories in imagedir which are 2014 or 2015 date format
dirs = [x for x in os.listdir(imagedir) if x[0:2] == '14' or x[0:2] == '15']
fulldirs = [imagedir+x+'/' for x in dirs]

# Read in reference data for the targets
refdata = ascii.read(reffile)
KICs = refdata['col1']
Porbs = refdata['col2']
BJD0s = refdata['col3']
RAs = refdata['col7']
Decs = refdata['col8']
# Create astropy Time objects for the zeropoints and orbital periods
Porbs_time = []; BJD0s_time = []
for Porb, BJD0 in zip(Porbs, BJD0s):
    Porbs_time.append(TimeDelta(Porb, format='jd'))			            # duration of one orbit
    BJD0s_time.append(Time(BJD0+2400000.0, format='jd', scale='utc'))   # time of primary eclipse
# Eclipse timing information
pwid = refdata['col4']
swid = refdata['col5']
sep = refdata['col6']

# Find the files that are FITS images
# Save the date, time, RA, Dec, and filter from the header, as well as the filename
# Keep going if a file can't be opened for any reason
dateobs = []; UTobs = []
RAobs = []; Decobs = []
filtnameobs = []; filenamesave = []
for dir in fulldirs:
    filesindir = os.listdir(dir)
    for filename in filesindir:
        # Keep only fits files that are not guiding images
        if filename[-4:] == ('fits' or 'FITS') and 'g.' not in filename:
            fullfile = dir+filename
            try:
                hdu = fits.open(fullfile, ignore_missing_end = True)
                dateobs.append(hdu[0].header['date-obs'])
                UTobs.append(hdu[0].header['UT'])
                RAobs.append(hdu[0].header['RA'])
                Decobs.append(hdu[0].header['Dec'])
                filtnameobs.append(hdu[0].header['filtname'][0:17])
                filenamesave.append(fullfile)
            except:
                print('Error accessing {0}'.format(fullfile))
                continue

# Put RA and Dec values into less annoying formats
print('Done reading image files')
RAs = coord.Angle(RAs, unit=u.hour)
RAobs = coord.Angle(RAobs, unit=u.hour)
Decs = coord.Angle(Decs, unit=u.degree)
Decobs = coord.Angle(Decobs, unit=u.degree)

# Identify which catalog RA and Dec value are closest to the observed ones
# If the closest RA and closest Dec have the same index, assign the appropriate KIC
# Otherwise, assume the file doesn't contain a target of interest
print('Identifying which observations belong to KIC objects...')
KICobs = []
for RA, Dec in zip(RAobs, Decobs):
    idx1 = min(range(len(RAs)), key=lambda i: abs(RAs[i] - RA))
    idx2 = min(range(len(Decs)), key=lambda i: abs(Decs[i] - Dec))
    if idx1 == idx2:
        KICobs.append(KICs[idx1])
    else:
        KICobs.append('None')

# Keep only the good observations that have assigned KICS
# Consolidate the time and date info into a single object
# Keep track of the filters and filenames
KICgoods = []; datetimes = []; RAgoods = []; Decgoods = []; filtgoods = []; filenamegoods = []
for KIC, date, UT, RA, Dec, filtname, file in zip(KICobs, dateobs, UTobs, RAobs, Decobs, filtnameobs, filenamesave):
    if KIC != 'None':
        KICgoods.append(KIC)
        datetimes.append(Time(str(date)+'T'+str(UT), format='isot', scale='utc'))
        RAgoods.append(RA)
        Decgoods.append(Dec)
        filtgoods.append(filtname)
        filenamegoods.append(file)
print('Done')

# Create a new list that contains a list of observation times for each object
# Keep track of the corresponding filters and filenames
# Initialize a parallel list that assumes all observations are not in eclipse, for now
observations = [[] for x in xrange(len(KICs))]
filtrefs = [[] for x in xrange(len(KICs))]
filenamerefs = [[] for x in xrange(len(KICs))]
isEclipse = [[] for x in xrange(len(KICs))]
for idx, (KIC, Porb, BJD0) in enumerate(zip(KICs, Porbs, BJD0s)): # loop over systems
    for KIC_obs, datetime_obs, filt_ref, file_ref in zip(KICgoods, datetimes, filtgoods, filenamegoods): # loop over observations
        if KIC_obs == KIC:
            observations[idx].append(datetime_obs)
            filtrefs[idx].append(filt_ref)
            filenamerefs[idx].append(file_ref)
            isEclipse[idx].append('No ')

# Define start and end of observation windows, with a 10-day buffer for plotting purposes
obs_tstart = min(datetimes) - TimeDelta(10, format='jd')
obs_tend = max(datetimes) + TimeDelta(10, format='jd')

# Calculate eclipse start and end points that fall within the observation window
# (This is blatantly stolen/adapted from 'eclipsefinder.py')
# These eclipse durations do NOT include any extra time for ensuring 100% of ingress/egress
print('Calculating eclipse times...')
print('New BJD0 values immediately before the obs window (KIC, newBJD0, start, end):')
pri_eclipse_mid = [[] for x in xrange(len(KICs))]
pri_eclipse_mid = [[] for x in xrange(len(KICs))]
sec_eclipse_mid = [[] for x in xrange(len(KICs))]
pri_eclipse_start = [[] for x in xrange(len(KICs))]
pri_eclipse_end = [[] for x in xrange(len(KICs))]
sec_eclipse_start = [[] for x in xrange(len(KICs))]
sec_eclipse_end = [[] for x in xrange(len(KICs))]
for j in range(0,len(KICs)):                        # j is the *object* index here
    # Find the most recent bjd0 time right BEFORE the observation window of interest
    newbjd0_float = np.floor((obs_tstart.jd - BJD0s_time[j].jd)/Porbs_time[j].value) * Porbs_time[j].value + BJD0s_time[j].jd
    newbjd0 = Time(newbjd0_float, format='jd', scale='utc')
    print(KICs[j], newbjd0_float, obs_tstart.jd, obs_tend.jd)
    for i in range(0,len(observations[j])):         # i is the *observation* index here
        # Save eclipse midpoints
        pri_eclipse_mid[j].append(newbjd0 + i*Porbs_time[j])
        sec_eclipse_mid[j].append(newbjd0 + i*Porbs_time[j] + sep[j]*Porbs_time[j])
        # Save primary eclipse start & end times
        pri_eclipse_start[j].append(pri_eclipse_mid[j][i] - pwid[j]*Porbs_time[j]/2)
        pri_eclipse_end[j].append(pri_eclipse_mid[j][i] + pwid[j]*Porbs_time[j]/2)
        # Save secondary eclipse start & end times
        sec_eclipse_start[j].append(sec_eclipse_mid[j][i] - swid[j]*Porbs_time[j]/2)
        sec_eclipse_end[j].append(sec_eclipse_mid[j][i] + swid[j]*Porbs_time[j]/2)
print('Done')

# Make a plot as a function of time
# Color-code points by whether they're in eclipse or not (red = primary, blue = secondary)
# Finally, save all the observation info to file
plt.figure(1, figsize=(17,12))
plt.yticks(range(0,len(KICs)), ['%.0f' % a for a in KICs])
plt.axis([obs_tstart.plot_date, obs_tend.plot_date, -1, len(KICs)])
print('Preparing plot and writing to outfile (be patient)...')
outfilelist = open(outfile, 'w')
for idx, KIC in enumerate(KICs): # loop over systems
    print(KIC, '...')
    for jdx, (obs, filt, file) in enumerate(zip(observations[idx], filtrefs[idx], filenamerefs[idx])): # loop over observations
        plt.plot_date(obs.plot_date, idx, marker='o', mec=None, mfc='0.75') # plot all observations
        for start, end in zip(pri_eclipse_start[idx], pri_eclipse_end[idx]): # PRIMARY ECLIPSE
            if start.jd <= obs.jd <= end.jd:
                isEclipse[idx][jdx] = 'Pri'   # change value in isEclipse if appropriate
                plt.plot_date(obs.plot_date, idx, marker='o', mec=None, mfc='r')
        for start, end in zip(sec_eclipse_start[idx], sec_eclipse_end[idx]): # SECONDARY ECLIPSE
            if start.jd <= obs.jd <= end.jd:
                isEclipse[idx][jdx] = 'Sec'   # change value in isEclipse if appropriate
                plt.plot_date(obs.plot_date, idx, marker='o', mec=None, mfc='b')
        print(KIC, obs.jd, isEclipse[idx][jdx], filt, file[-15:], file=outfilelist)
outfilelist.close()
plt.gcf().autofmt_xdate() # for slanty dates on the Time x-axis

# Option to show or save figure
#plt.show()
plt.savefig('imginventory_plot.png')
