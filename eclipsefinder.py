from __future__ import print_function
import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta
import matplotlib.pyplot as plt
'''
Python program to calculate upcoming EB eclipse times and durations.
Useful for figuring out when to observe, or not.
Originally written by Meredith Rawls, March-April 2014

INPUT:
A text file with six columns:
ID - whatever you want to call your EB
P - orbital period
T0 - zeropoint (presumably in BJD)
pwidth - length of primary eclipse in phase units
swidth - length of secondary eclipse in phase units
sep - separation between two eclipses in phase units
RA - 10-character RA value 00:00:00.0
Dec - 10-character Dec value 00:00:00.0

OUTPUT:
A graph of when eclipses will occur.
Eclipse timing data in two text files: one that is human-friendly and one that is 
for use as a control file by the 1m telescope at APO. (You will need to manually edit
the desired exposure times, etc. in this program if you care about the latter.)

Read the damn comments, Jean.
'''

# Function to make a 2D list without stabbing your eye out
def make2dList(rows, cols):
    a=[]
    for row in xrange(rows): a += [[0]*cols]
    return a

infile = '../../1m_observations/RGEB_info_alpha.txt'		# infile you must provide
outfile1 = '../../1m_observations/eclipses_2015_Q3.txt'	# human-friendly outfile
outfile2 = '../../1m_observations/RGEB_1m_2015_test.inp' 	# 1m telescope friendly outfile

# Create a list of possible observation times
obs_tstart = Time('2015-07-01 00:00', format='iso', scale='utc')
obs_tend = Time('2015-09-30 00:00', format='iso', scale='utc')
obs_deltat = TimeDelta(0.5, format='jd')
# arbitrary fraction-of-a-day time resolution... 0.5 is probably sufficient for eclipse
# durations ~1 day (?)
# make sure to set obs_deltat < Porb for all systems or bad things might happen

# Read in data from file
kic, p, t0, pwid, swid, sep, RA, Dec = np.loadtxt(infile, comments='#', usecols=(0,1,2,3,4,5,6,7), 
		dtype={'names': ('kic', 'p', 't0', 'pwid', 'swid', 'sep', 'RA', 'Dec'),
		'formats': (np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,'|S11','|S11')}, unpack=True)
print ('Reading data from ' + infile)

# Create astropy Time objects for the zeropoints and orbital periods
bjd0 = Time(t0+2400000.0, format='jd', scale='utc')	# time of primary eclipse
porb = TimeDelta(p, format='jd')			# duration of one orbit

tobs = Time(np.arange(obs_tstart.jd, obs_tend.jd, obs_deltat.jd), format='jd', scale='utc')

print ('Forecasting eclipses for this time range:')
print (obs_tstart.iso + ' to ' + obs_tend.iso)
print (obs_tstart.jd, obs_tend.jd)
print (' ')

# Access one element of this observation time array in JD format:
#print tobs[0].utc.jd

# What we'd really like is a list of times when eclipses BEGIN and END that fall
# between obs_tstart and obs_tend. Want to work in real time, not phase-time.
# --> astropy Time values DO NOT play nice with numpy arrays!!!! Use lists instead.

# This gives WAY MORE eclipse times than you could ever want, assuming Porb > obs_deltat
newbjd0_float = []
newbjd0 = []
pri_eclipse_mid = make2dList(len(kic), len(tobs))
sec_eclipse_mid = make2dList(len(kic), len(tobs))
pri_eclipse_start = make2dList(len(kic), len(tobs))
pri_eclipse_end = make2dList(len(kic), len(tobs))
sec_eclipse_start = make2dList(len(kic), len(tobs))
sec_eclipse_end = make2dList(len(kic), len(tobs))
for j in range(0,len(kic)):
	# Find the most recent bjd0 time right BEFORE the observation window of interest
	newbjd0_float.append( np.floor((obs_tstart.jd - bjd0[j].jd)/porb[j].value) * porb[j].value + bjd0[j].jd )
	newbjd0.append( Time(newbjd0_float[j], format='jd', scale='utc') )
	i = 0
	for i in range(0,len(tobs)):
		# Save eclipse midpoints
		pri_eclipse_mid[j][i] = newbjd0[j] + i*porb[j]
		sec_eclipse_mid[j][i] = newbjd0[j] + i*porb[j] + sep[j]*porb[j]
		
		# print j, i, pri_eclipse_mid[j][i].iso # IT WORKS!

		# Save primary eclipse start & end times
		pri_eclipse_start[j][i] = pri_eclipse_mid[j][i] - pwid[j]*porb[j]/2
		pri_eclipse_end[j][i] = pri_eclipse_mid[j][i] + pwid[j]*porb[j]/2
		# Save secondary eclipse start & end times
		sec_eclipse_start[j][i] = sec_eclipse_mid[j][i] - swid[j]*porb[j]/2
		sec_eclipse_end[j][i] = sec_eclipse_mid[j][i] + swid[j]*porb[j]/2

# The 1-m telescope 'inp' file needs lots of other info too. Here we go.

## FIRSTLINE INFO ##
Epoch = '\t 2000.0'
Xmin = '\t 1.01'
Xmax = '\t 2.00'
#Year = '\t 2007'
Blank = '\t 24.0' #yeah I don't know either
UTstart = '\t 0.00'
UTend = '\t 24.00'
HAmin = '\t -12.00'
HAmax = '\t 12.00'
MPhase = '\t 1.00'
MDist = '\t 5.00'
Nexp = '\t -1'
DT = '\t 3600' # for file1, make this ~ 1 hour = 3600 sec, and for file2, make this 0
Acquire = '\t -1 \t 0 \t 3 \t 0 \t 0'
# -1 is find guiding star & focus. -2 is just find a guiding star. -3 is just go.
# the trailing 0 3 0 0 are for something with guiding

## SECONDLINE INFO ##
#seq_repeat = '1' #this is calculated below instead
# NOTE exposure times manually correspond with the infile order (sorry)
nU = '0'
tU = '0'
nB = 1
tB = [240,120,120,120,75,75,240,90,240,210,75,10,10,90,5,240,75,75,10]
nV = 1
tV = [120,60,60,60,30,30,120,60,120,120,30,5,5,60,5,120,30,30,5]
nR = 1
tR = [45,30,30,30,15,15,45,30,45,60,15,5,5,30,5,60,15,15,5]
nI = 1
tI = [60,45,45,45,20,20,60,45,60,90,20,5,5,45,5,90,20,20,5]
end_line2 = '0 \t 0.0 \t 0 \t 0.0 \t 0 \t 0.0 \t 0 \t 0.0 \t 0 \t 0.0'

# Plot the eclipses as a function of time
# And write useful info to text files!
f1 = open(outfile1, 'w')
f2 = open(outfile2, 'w')
plt.yticks(range(1,len(kic)+1), ['%.0f' % a for a in kic])
tplotstart = obs_tstart
tplotend = obs_tend
plt.axis([tplotstart.plot_date, tplotend.plot_date, 0, len(kic)+1])
# epic loop for plotting AND writing info to files
for j in range(0,len(kic)):
	# Calculate how many times to repeat the sequence for ~30 min (1800 sec)
	seq_repeat = int(1800 / (tB[j]+20 + tV[j]+20 + tR[j]+20 + tI[j]+20) )
	i = 0
	while (pri_eclipse_start[j][i] < tplotend or sec_eclipse_start[j][i] < tplotend):
	# the while loop below should work, but doesn't, so as it stands (above) there are
	# some extra points that occur before tplotstart. sigh.
	#while ( (pri_eclipse_start[j][i] < tplotend and pri_eclipse_start[j][i] > tplotstart) or (sec_eclipse_start[j][i] < tplotend and sec_eclipse_start[j][i] > tplotstart) ):
		
		# Put stuff on the graph
		# plt.plot([x1,x2],[y1,y2], 'colors-n-stuff') #syntax to draw some lines
		plt.plot_date([pri_eclipse_start[j][i].plot_date, pri_eclipse_end[j][i].plot_date],[j+1, j+1], 'k-', lw=2)#, label='primary' if j==0 else '')
		plt.plot_date([sec_eclipse_start[j][i].plot_date, sec_eclipse_end[j][i].plot_date],[j+1, j+1], 'r-', lw=2)#, label='secondary' if j==0 else '')
		
		# Write just eclipse timing info to a file
		print(int(kic[j]), '\t p \t', pri_eclipse_start[j][i].isot, '\t', 
			pri_eclipse_mid[j][i].isot, '\t', pri_eclipse_end[j][i].isot, file=f1)
		print(int(kic[j]), '\t s \t', sec_eclipse_start[j][i].isot, '\t',
			sec_eclipse_mid[j][i].isot, '\t', sec_eclipse_end[j][i].isot, file=f1)
		
		# Define MS, DS, ME, DE (Month/Day Start/End) windows from eclipse times
		Yearpri = pri_eclipse_start[j][i].iso[0:4]
		MSpri = pri_eclipse_start[j][i].iso[5:7]
		DSpri = pri_eclipse_start[j][i].iso[8:10]
		MEpri = pri_eclipse_end[j][i].iso[5:7]
		DEpri = pri_eclipse_end[j][i].iso[8:10]
		Yearsec = sec_eclipse_start[j][i].iso[0:4]
		MSsec = sec_eclipse_start[j][i].iso[5:7]
		DSsec = sec_eclipse_start[j][i].iso[8:10]
		MEsec = sec_eclipse_end[j][i].iso[5:7]
		DEsec = sec_eclipse_end[j][i].iso[8:10]
		
		# Create new RGEB.inp file for 1m observations
		# LINE 1, primary eclipse
		print(int(kic[j]), '\t \t', RA[j], '\t', Dec[j], Epoch, Xmin, Xmax, Yearpri, '\t', 
			MSpri, '\t', DSpri, '\t', MEpri, '\t', DEpri, '\t', Blank, UTstart, UTend,
			HAmin, HAmax, MPhase, MDist, Nexp, DT, Acquire, file=f2)
		# LINE 2, primary eclipse
		print('\t', seq_repeat, '\t', nU, '\t', tU, '\t', nB, '\t', tB[j], '\t', nV, '\t',
			tV[j], '\t', nR, '\t', tR[j], '\t', nI, '\t', tI[j], '\t', end_line2, file=f2)
		# LINE 1, secondary eclipse
		print(int(kic[j]), '\t \t', RA[j], '\t', Dec[j], Epoch, Xmin, Xmax, Yearsec, '\t', 
			MSsec, '\t', DSsec, '\t', MEsec, '\t', DEsec, '\t', Blank, UTstart, UTend,
			HAmin, HAmax, MPhase, MDist, Nexp, DT, Acquire, file=f2)
		# LINE 2, secondary eclipse
		print('\t', seq_repeat, '\t', nU, '\t', tU, '\t', nB, '\t', tB[j], '\t', nV, '\t',
			tV[j], '\t', nR, '\t', tR[j], '\t', nI, '\t', tI[j], '\t', end_line2, file=f2)
		i = i + 1

# print each system one final time (with no date restrictions, and seq_repeat = 5)
print ('# The following is for file 3! Each system is printed once with no time restrictions', file=f2)
print ('# ', file=f2)
for j in range(0,len(kic)):
	# LINE 1, no date restrictions
	print(int(kic[j]), '\t \t', RA[j], '\t', Dec[j], Epoch, Xmin, Xmax, '2014', '\t', 
		'0', '\t', '0', '\t', '12', '\t', '31', '\t', Blank, UTstart, UTend,
		HAmin, HAmax, MPhase, MDist, '\t 1', DT, Acquire, file=f2)
	# LINE 2, no date restrictions
	print('\t', '5', '\t', nU, '\t', tU, '\t', nB, '\t', tB[j], '\t', nV, '\t',
		tV[j], '\t', nR, '\t', tR[j], '\t', nI, '\t', tI[j], '\t', end_line2, file=f2)

# finish the plot
#plt.xlabel("Time (BJD)")
#plt.ylabel("System")
plt.gcf().autofmt_xdate() #for slanty dates

f1.close()
#plt.legend()

## option to manually plot horizontal lines at dates of interest ##
# windowstart = Time('2015-04-07T09:16:00', format='isot', scale='utc') 
# windowend = Time('2015-04-07T12:35:00', format='isot', scale='utc')
# windowstart2 = Time('2015-04-16T09:16:00', format='isot', scale='utc') 
# windowend2 = Time('2015-04-16T12:24:00', format='isot', scale='utc')
# windowstart3 = Time('2015-04-27T07:16:00', format='isot', scale='utc')
# windowend3 = Time('2015-04-27T12:12:00', format='isot', scale='utc')
# windowstart4 = Time('2015-05-06T07:16:00', format='isot', scale='utc')
# windowend4 = Time('2015-05-06T12:03:00', format='isot', scale='utc')
# plt.axvline(windowstart.plot_date)
# plt.axvline(windowend.plot_date, ls=':')
# plt.axvline(windowstart2.plot_date)
# plt.axvline(windowend2.plot_date, ls=':')
# plt.axvline(windowstart3.plot_date)
# plt.axvline(windowend3.plot_date, ls=':')
# plt.axvline(windowstart4.plot_date)
# plt.axvline(windowend4.plot_date, ls=':')

plt.title('(black primary, red secondary)')
plt.show()
# primary eclipses are plotted in black
# secondary eclipses are plotted in red
