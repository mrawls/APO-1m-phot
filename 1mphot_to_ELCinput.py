from __future__ import print_function
import numpy as np
'''
Takes output from photbinner.py
Creates four BVRI input files for ELC ** that use Kepler BJDs **
PAY ATTENTION THE TIMESTAMP IS BEING ADJUSTED FROM TRUNCATED UTC JD TO JD-2454833
'''

dir = '../../1m_observations/KIC09291629/'
infile = dir+'BVRI_diffmag_binned_LC1.txt'
filters, UTCtimes, mags, merrs = np.loadtxt(infile, usecols=(0,1,3,4), unpack=True,
    dtype={'names': ('filters', 'UTCtimes', 'mags', 'merrs'), 'formats': ('|S2', np.float64, np.float64, np.float64)})
Keptimes = UTCtimes + 2400000. - 2454833.

Boutfile = dir + 'Bmags.txt'
Voutfile = dir + 'Vmags.txt'
Routfile = dir + 'Rmags.txt'
Ioutfile = dir + 'Imags.txt'

Bout = open(Boutfile, 'w')
Vout = open(Voutfile, 'w')
Rout = open(Routfile, 'w')
Iout = open(Ioutfile, 'w')

for filter, Keptime, mag, merr in zip(filters, Keptimes, mags, merrs):
    if filter == 'B':
        print(Keptime, mag, merr, file=Bout)
    elif filter == 'V':
        print(Keptime, mag, merr, file=Vout)
    elif filter == 'R':
        print(Keptime, mag, merr, file=Rout)
    elif filter == 'I':
        print(Keptime, mag, merr, file=Iout)
    else:
        print('filter not recognized')
        
Bout.close()
Vout.close()
Rout.close()
Iout.close()