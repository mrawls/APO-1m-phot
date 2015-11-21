from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import string
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
from pyraf import iraf
from pyraf.iraf import noao,imred,ccdred
'''
Meredith Rawls, Fall 2015

***Run imginventory.py first!!

NOTE: you likely need to manually edit the "filter" column so there are no spaces!
NOTE: you need to ensure that each image filename begins with yymmdd!
NOTE: if any of the file index numbers are > 999, you might have an interesting time!

Goal of this program: for one target star, do basic image reductions on APO 1m images.
(This means overscan (bias) subtraction and flat field division.)
Save the reduced images in a single directory for that target.

***Run this program once for each target with FlatDivide = False.
***!!LOOK AT THE IMAGES IT MADE AND DECIDE IF THERE ARE SOME SNEAKY GARBAGE ONES OR NOT!!
***Then, run this program once again for each target with FlatDivide = True.

INPUT: - folders named yymmdd in 'workingdir' with images that are named yymmdd.xxx.fits
       - many (not all) folders must contain flat field images, and text files that 
       indicate which images are flats
       - 'imginventory_list.txt' which you can make by running 'imginventory.py' first
       (see important NOTES above)
       - target name corresponding to first column in 'imginventory_list.txt'
       - FlatDivide flag (only set to True once all combined flat field images exist)

OUTPUT: - after first set of runs, combined median flat field images, with one
          assigned to each target image
        - after second set of runs, reduced images organized by target folder
        (you can customize these folder names as desired near the end of the program)

'''
target = '3955867'
#target = '4569590'
#target = '5179609'
#target = '5308778'
#target = '5640750'
#target = '5786154'
#target = '7037405'
#target = '7377422'
#target = '7943602'
#target = '8054233'
#target = '8410637'
#target = '8430105'
#target = '8702921'
#target = '9246715'
#target = '9291629'
#target = '9540226'
#target = '9970396'
#target = '10001167'

workingdir =    '/virgo/mrawls/1mphot/'     # assumes everything is in this directory
inventoryfile = 'imginventory_list2.txt'    # output file from imginventory.py
biassec = '[2110:2200,1:2048]'  # set colstart:colend, rowstart:rowend for overscan region

FlatDivide = False #True
# set to True once you've run the program once for all targets with it set to False
# (this way you won't try to divide by combined flat field images that don't exist yet)

# Read in image inventory file and assign full filename paths for each observation
allobjects, alltimes, alleclipses, allfilters, allfilenames = np.loadtxt(
    workingdir+inventoryfile, usecols=(0,1,2,3,4), 
    dtype={'names': ('obj', 'tim', 'ecl', 'filt', 'file'), 
    'formats': ('|S9', np.float64, '|S3', '|S6', '|S16')}, unpack=True)
allfullfilepaths = []
for filename in allfilenames:
    allfullfilepaths.append(workingdir + str(filename)[0:6] + '/' + str(filename))

# Save only the images corresponding to the present target
times = []
eclipses = []
filters = []
fullfilepaths = []
for object, time, eclipse, filter, fullfilepath in zip(allobjects, alltimes, alleclipses, allfilters, allfullfilepaths):
    if object == target:
        times.append(time)
        eclipses.append(eclipse)
        filters.append(filter)
        fullfilepaths.append(fullfilepath)

## Identify the corresponding flat field images for each target
## Create a median flat with FLATCOMBINE
## If no flat field images exist in the present directory, ### NOT SURE WHAT TO DO ###
flatfiles = [] # will be a parallel list to fullfilepaths containing combined flats
for index, (fullfilepath, filter) in enumerate(zip(fullfilepaths, filters)):   
    if fullfilepath[-9] != '.': # file index is > 999
        dir = fullfilepath[0:-16]
    else:
        dir = fullfilepath[0:-15]
    date = dir[-7:-1]
    #print(fullfilepath, dir, date) # OK   
    ## Identify the appropriate flatid text file based on filter
    ## This if-else section sets the name of the flat to use for each image
    if filter == 'B': flatid = dir + 'flat.b'        
    elif filter == 'V': flatid = dir + 'flat.v' 
    elif filter == 'R': flatid = dir + 'flat.r' 
    elif filter == 'I': flatid = dir + 'flat.i'
    elif filter == 'SDSS-g': flatid = dir + 'flat.sdssg'
    elif filter == 'SDSS-i': flatid = dir + 'flat.sdssi'
    elif filter == 'SDSS-r': flatid = dir + 'flat.sdssr'
    elif filter == 'SDSS-u': flatid = dir + 'flat.sdssu'
    elif filter == 'SDSS-z': flatid = dir + 'flat.sdssz'
    else:
        print('ERROR: flatfield image not assigned for {0}'.format(fullfilepath[-15:]))
        flatid = dir + 'NONE'
    # This is the name of the combined flat we will create (if we haven't already)
    fullflatpath = flatid + '_median.fits'   
    ## This try-except-else-finally section creates and remembers a *combined* flat image 
    try:    # see if the combined flat already exists
        fits.open(fullflatpath)
    except: # if it doesn't exist, create it        
        ## This try-except-else section checks if *individual* flat images exist
        try: # see if the flatid file exists (e.g., 'flat.b')
            open(flatid) # if this succeeds, continue
        except: # no file named flatid exists (there are no flats in this directory)
            print('No flat fields found for {0} in {1}'.format(filter, dir))
            ## We don't need to make a new file, but we DO need to figure out the name of
            ## the combined flat that has been OR will be created in another directory.
            # --> try accessing flatid files in dir with date-1, date+1...
            # --> continue this process until a flatid file is found
            datestring = '20'+date[0:2]+'-'+date[2:4]+'-'+date[4:6] # shoot me
            dts = [-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7, -8, 8, -9, 9]           
            for dt in dts:
                newdate = Time(datestring, in_subfmt='date') + TimeDelta(dt, format='jd')
                newdate = str(newdate)[2:4]+str(newdate)[5:7]+str(newdate)[8:10]
                testflatpath = workingdir+newdate+'/'+flatid[-6:]
                #print(newdate, testflatpath) # OK
                try:
                    open(testflatpath)
                except:
                    continue
                else:
                    fullflatpath = testflatpath + '_median.fits'
                    print('--> Using {0} instead'.format(fullflatpath))
                    break
        else: # if the flatid file DOES exist, we need to make a combined flat
            with open(flatid) as f1: # open flatid to learn which images to combine
                imstart, imend = np.loadtxt(f1, usecols=(0,1))
                imstart = int(imstart); imend = int(imend)
            # Make a list of the flat images (one for raw, one for overscan corrected)
            flatstocombine = []; zflatstocombine = []
            for idx in range(imstart, imend+1):
                if idx < 10: idx = '00'+str(idx)
                elif idx > 10 and idx < 100: idx = '0'+str(idx)
                else: idx = str(idx)
                flatstocombine.append(dir+date+'.'+idx+'.fits')
                zflatstocombine.append(dir+'zflat'+date+'.'+idx+'.fits')
            ## This try-except section applies overscan correction to flats
            try: # see if the overscan correction has already been done
                for zflat in zflatstocombine:
                    test = fits.open(zflat)
                print('Overscan corrected flat fields for {0} already exist in {1}'.format(filter, dir))
            except: # if necessary, do the overscan correction
                zflatsgood = []
                for flat, zflat in zip(flatstocombine, zflatstocombine):
                    try: # if an *individual* flat image is damaged, this will fail
                        iraf.ccdproc('\''+flat+'\'', output='\''+zflat+'\'', ccdtype='none',
                            noproc='no', fixpix='no', oversca='yes', trim='no', zerocor='no', 
                            darkcor='no', flatcor='no', interactive='no', biassec=biassec)
                        zflatsgood.append(zflat)
                        print('Flat field image overscan corrected {0}'.format(flat))
                    except: # skip the problematic flat image and move on
                        print('Error accessing flat field image {0}'.format(flat))
            zflatstring = string.join(zflatsgood, ',')
            print('Flats to combine for {0}: {1}'.format(filter, zflatstring))
            iraf.imcombine(zflatstring, fullflatpath, 
                combine='median', reject = 'avsigclip', mclip = 'yes')                
            print('Flat field combining done for {0} in {1}'.format(filter, dir))                    
    else: # let the user know that the combined flat exists, so we will skip ahead
        print('Combined flat already exists for {0} in {1}'.format(filter, dir))                    
    finally: # at the end of this mess, save the name of the combined flat file
        flatfiles.append(fullflatpath)
print(' ')
print('That was fun! If these numbers are equal and there are no IRAF PANICs')
print('or other errors above, you should have flats for all images.')
print('Target: {0}'.format(target))
print('Number of target images, assigned flats: {0}, {1}'.format(len(fullfilepaths), len(flatfiles)))


# CCDPROC party time: overscan bias subtraction and flatfielding of targets! Finally!
if FlatDivide == True:
## loop over all (fullfilepaths, flatfiles)
## run CCDPROC to divide each image by flat, and save a new image!
    print(' ')
    print('PROCEEDING WITH FLAT DIVISION!')
    for image, flat, eclipse, filter in zip(fullfilepaths, flatfiles, eclipses, filters):
        # picky filenames are picky
        # (they're a function of the number of digits in the target ID and image number)
        if int(target) > 9999999 and image[-9] == '.':
            flattenedimage = workingdir + 'KIC' + target + '/' + image[-15:-5] + eclipse + filter + '.fits'
        elif int(target) <= 9999999 and image[-9] == '.':
            flattenedimage = workingdir + 'KIC0' + target + '/' + image[-15:-5] + eclipse + filter + '.fits'
        elif int(target) > 9999999 and image[-9] != '.':
            flattenedimage = workingdir + 'KIC' + target + '/' + image[-16:-5] + eclipse + filter + '.fits'
        elif int(target) <= 9999999 and image[-9] != '.':
            flattenedimage = workingdir + 'KIC0' + target + '/' + image[-16:-5] + eclipse + filter + '.fits'
        iraf.ccdproc('\''+image+'\'', output='\''+flattenedimage+'\'', ccdtype='none',
            noproc='no', fixpix='no', oversca='yes', trim='no', zerocor='no', 
            darkcor='no', flatcor='yes', interactive='no', biassec=biassec, flat=flat)


# THE HARD PART, for another program
# - identify the target star in the image, and N nearby bright stars too
# - aperture photometry for days
