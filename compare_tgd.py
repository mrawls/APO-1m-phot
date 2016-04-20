from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
'''
Compare multiple measurements of Teff, logg, and distance for RG/EBs.
'''
dfile = 'RGEB_distinfo.txt'
tgfile = '../../RGEB_tefflogg.txt'
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']

# distance file column info
KICcol1 = 0
distcol = 7; disterrcol = 8
DR12distcol = 9; DR12disterrcol = 10

# teff & logg file column info
KICcol2 = 0
MOOGtcol = 1;   MOOGterrcol = 2
MOOGgcol = 3;   MOOGgerrcol = 4
MOOGfecol = 5;  MOOGfeerrcol = 6
DR12tcol = 7;   DR12terrcol = 8
DR12gcol = 9;   DR12gerrcol = 10
DR12fecol = 11; DR12feerrcol = 12
Cannontcol = 13; Cannonterrcol = 14
Cannongcol = 15; Cannongerrcol = 16
Cannonfecol = 17; Cannonfeerrcol = 18
KICtcol = 19;   KICterrcol = 20
KICgcol = 21;   KICgerrcol = 22
KICfecol = 23;  KICfeerrcol = 24
MESAtcol = 25;  MESAterrcol = 26
MESAgcol = 27;  MESAgerrcol = 28

# logg column info only
ELCgcol = 29;   ELCgerrcol = 30
Seismicgcol = 31; Seismicgerrcol = 32

# we'll use KICcol1 as the master target list because it's numerically sorted
KICs = np.loadtxt(dfile, comments='#', usecols=(0,), unpack=True)
KICtgs = np.loadtxt(tgfile, comments='#', usecols=(0,), unpack=True)

# create x axis label string list
strKICs = [' ',]
for KIC in KICs:
    strKICs.append(str(int(KIC)))
strKICs.append(' ')
xaxis = np.arange(0,len(KICs))

# load the rest of the data (distances)
dists = np.loadtxt(dfile, comments='#', usecols=(distcol, disterrcol), unpack=True)
DR12ds = np.loadtxt(dfile, comments='#', usecols=(DR12distcol, DR12disterrcol), unpack=True)

# load the rest of the data (teffs, loggs)
MOOGts = np.loadtxt(tgfile, comments='#', usecols=(MOOGtcol, MOOGterrcol), unpack=True)
MOOGgs = np.loadtxt(tgfile, comments='#', usecols=(MOOGgcol, MOOGgerrcol), unpack=True)
MOOGfes = np.loadtxt(tgfile, comments='#', usecols=(MOOGfecol, MOOGfeerrcol), unpack=True)
DR12ts = np.loadtxt(tgfile, comments='#', usecols=(DR12tcol, DR12terrcol), unpack=True)
DR12gs = np.loadtxt(tgfile, comments='#', usecols=(DR12gcol, DR12gerrcol), unpack=True)
DR12fes = np.loadtxt(tgfile, comments='#', usecols=(DR12fecol, DR12feerrcol), unpack=True)
Cannonts = np.loadtxt(tgfile, comments='#', usecols=(Cannontcol, Cannonterrcol), unpack=True)
Cannongs = np.loadtxt(tgfile, comments='#', usecols=(Cannongcol, Cannongerrcol), unpack=True)
Cannonfes = np.loadtxt(tgfile, comments='#', usecols=(Cannonfecol, Cannonfeerrcol), unpack=True)
KICts = np.loadtxt(tgfile, comments='#', usecols=(KICtcol, KICterrcol), unpack=True)
KICgs = np.loadtxt(tgfile, comments='#', usecols=(KICgcol, KICgerrcol), unpack=True)
KICfes = np.loadtxt(tgfile, comments='#', usecols=(KICfecol, KICfeerrcol), unpack=True)
MESAts = np.loadtxt(tgfile, comments='#', usecols=(MESAtcol, MESAterrcol), unpack=True)
MESAgs = np.loadtxt(tgfile, comments='#', usecols=(MESAgcol, MESAgerrcol), unpack=True)
ELCgs = np.loadtxt(tgfile, comments='#', usecols=(ELCgcol, ELCgerrcol), unpack=True)
#Seismicgs = np.loadtxt(tgfile, comments='#', usecols=(Seismicgcol, Seismicgerrcol), unpack=True)

# sort all the teff and logg arrays so they align with xaxis
MOOGts[0] = MOOGts[0][np.argsort(KICtgs)];     MOOGgs[0] = MOOGgs[0][np.argsort(KICtgs)]
MOOGts[1] = MOOGts[1][np.argsort(KICtgs)];     MOOGgs[1] = MOOGgs[1][np.argsort(KICtgs)]
DR12ts[0] = DR12ts[0][np.argsort(KICtgs)];     DR12gs[0] = DR12gs[0][np.argsort(KICtgs)]
DR12ts[1] = DR12ts[1][np.argsort(KICtgs)];     DR12gs[1] = DR12gs[1][np.argsort(KICtgs)]
Cannonts[0] = Cannonts[0][np.argsort(KICtgs)]; Cannongs[0] = Cannongs[0][np.argsort(KICtgs)]
Cannonts[1] = Cannonts[1][np.argsort(KICtgs)]; Cannongs[1] = Cannongs[1][np.argsort(KICtgs)]
KICts[0] = KICts[0][np.argsort(KICtgs)];       KICgs[0] = KICgs[0][np.argsort(KICtgs)]
KICts[1] = KICts[1][np.argsort(KICtgs)];       KICgs[1] = KICgs[1][np.argsort(KICtgs)]
MESAts[0] = MESAts[0][np.argsort(KICtgs)];     MESAgs[0] = MESAgs[0][np.argsort(KICtgs)]
MESAts[1] = MESAts[1][np.argsort(KICtgs)];     MESAgs[1] = MESAgs[1][np.argsort(KICtgs)]
ELCgs[0] = ELCgs[0][np.argsort(KICtgs)]
ELCgs[1] = ELCgs[1][np.argsort(KICtgs)]
#Seismicgs[0] = Seismicgs[0][np.argsort(KICtgs)]
#Seismicgs[1] = Seismicgs[1][np.argsort(KICtgs)]
MOOGfes[0] = MOOGfes[0][np.argsort(KICtgs)];        MOOGfes[1] = MOOGfes[1][np.argsort(KICtgs)]
DR12fes[0] = DR12fes[0][np.argsort(KICtgs)];        DR12fes[1] = DR12fes[1][np.argsort(KICtgs)]
Cannonfes[0] = Cannonfes[0][np.argsort(KICtgs)];    Cannonfes[1] = Cannonfes[1][np.argsort(KICtgs)]
KICfes[0] = KICfes[0][np.argsort(KICtgs)];          KICfes[1] = KICfes[1][np.argsort(KICtgs)]
# WHEW THAT WAS FUN
KICtgs = KICtgs[np.argsort(KICtgs)]

# set up figures
fig = plt.figure(1, figsize=(15,10))
plt.subplots_adjust(hspace=0)

# top panel - TEMPERATURES
ax1 = plt.subplot(4,1,1)
ax1.set_ylabel(r'$T_{\textrm{eff}}$ (K)', size=26)
ax1.set_xlim([-0.5,len(xaxis)-0.5])
ax1.set_ylim([4405, 5200])
ax1.set_xticklabels(strKICs, size=22)
ax1.xaxis.tick_top()
plt.errorbar(xaxis-0.1, DR12ts[0], yerr=DR12ts[1], ls='None', marker='o', ms=10, c=colors[1], label='DR12')
plt.errorbar(xaxis-0.05, Cannonts[0], yerr=Cannonts[1], ls='None', marker='s', ms=10, c=colors[2], label='Cannon')
plt.errorbar(xaxis+0.05, KICts[0], yerr=KICts[1], ls='None', marker='^', ms=10, c=colors[3], label='KIC')
plt.errorbar(xaxis+0.1, MESAts[0], ls='None', marker='h', ms=12, c=colors[4], label='MESA')
plt.errorbar(xaxis, MOOGts[0], yerr=MOOGts[1], ls='None', marker='D', ms=10, c=colors[0], label='MOOG')
plt.axvspan(0.5,1.5, alpha=0.1, color='k')
plt.axvspan(2.5,3.5, alpha=0.1, color='k')
plt.axvspan(4.5,5.5, alpha=0.1, color='k')

#leg1 = ax1.legend(bbox_to_anchor=(1.1,0.68), numpoints=1, loc=1, borderaxespad=0., 
#    frameon=True, handletextpad=0.2, prop={'size':18})
#leg1.get_frame().set_linewidth(0.0)

leg1 = ax1.legend(bbox_to_anchor=(1.1,0.90), numpoints=1, loc=1, borderaxespad=0., 
    frameon=True, handletextpad=0.2, prop={'size':18})
leg1.get_frame().set_linewidth(0.0)

# middle-top panel - LOGGS
ax2 = plt.subplot(4,1,2)
ax2.set_ylabel(r'$\log g$ (cgs)', size=26)
ax2.set_xlim([-0.5,len(xaxis)-0.5])
ax2.set_ylim([1.8, 3.7])
ax2.set_xticklabels([])
plt.errorbar(xaxis-0.1, DR12gs[0], yerr=DR12gs[1], ls='None', marker='o', ms=10, c=colors[1])
plt.errorbar(xaxis-0.05, Cannongs[0], yerr=Cannongs[1], ls='None', marker='s', ms=10, c=colors[2])
plt.errorbar(xaxis+0.05, KICgs[0], yerr=KICgs[1], ls='None', marker='^', ms=10, c=colors[3])
plt.errorbar(xaxis+0.1, MESAgs[0], ls='None', marker='h', ms=12, c=colors[4])
plt.errorbar(xaxis, MOOGgs[0], yerr=MOOGgs[1], ls='None', marker='D', ms=10, c=colors[0])
plt.errorbar(xaxis, ELCgs[0], yerr=ELCgs[1], ls='None', marker='*', ms=20, c=colors[0], label='This work')
#plt.errorbar(xaxis, Seismicgs[0], yerr=Seismicgs[1], ls='None', marker='*', ms=20, c='k', label='Seismic')
plt.axvspan(0.5,1.5, alpha=0.1, color='k')
plt.axvspan(2.5,3.5, alpha=0.1, color='k')
plt.axvspan(4.5,5.5, alpha=0.1, color='k')

#leg2 = ax2.legend(bbox_to_anchor=(1.12,1.01), numpoints=1, loc=1, borderaxespad=0., 
#    frameon=True, handletextpad=0.2, prop={'size':18})
#leg2.get_frame().set_linewidth(0.0)

leg2 = ax2.legend(bbox_to_anchor=(1.12,1.005), numpoints=1, loc=1, borderaxespad=0., 
    frameon=True, handletextpad=0.2, prop={'size':18})
leg2.get_frame().set_linewidth(0.0)

# middle-bottom panel - METALLICITIES
ax3 = plt.subplot(4,1,3)
ax3.set_ylabel(r'[Fe/H]', size=26)
ax3.set_xlim([-0.5,len(xaxis)-0.5])
ax3.set_ylim([-0.9, 0.55])
ax3.set_xticklabels([])
plt.errorbar(xaxis-0.1, DR12fes[0], yerr=DR12fes[1], ls='None', marker='o', ms=10, c=colors[1])
plt.errorbar(xaxis-0.05, Cannonfes[0], yerr=Cannonfes[1], ls='None', marker='s', ms=10, c=colors[2])
plt.errorbar(xaxis+0.05, KICfes[0], yerr=KICfes[1], ls='None', marker='^', ms=10, c=colors[3])
plt.errorbar(xaxis, MOOGfes[0], yerr=MOOGfes[1], ls='None', marker='D', ms=10, c=colors[0])
plt.axvspan(0.5,1.5, alpha=0.1, color='k')
plt.axvspan(2.5,3.5, alpha=0.1, color='k')
plt.axvspan(4.5,5.5, alpha=0.1, color='k')

# bottom panel - DISTANCES
ax4 = plt.subplot(4,1,4)
ax4.set_ylabel(r'$d$ (kpc)', size=26)
ax4.set_xlim([-0.5,len(xaxis)-0.5])
ax4.set_ylim([0.5, 3.9])
ax4.set_xticklabels(strKICs, size=22)
plt.errorbar(xaxis-0.1, DR12ds[0]/1000, yerr=DR12ds[1]/1000, ls='None', marker='o', ms=10, c=colors[1])
plt.errorbar(xaxis, dists[0]/1000, yerr=dists[1]/1000, ls='None', marker='*', ms=20, c=colors[0])
plt.axvspan(0.5,1.5, alpha=0.1, color='k')
plt.axvspan(2.5,3.5, alpha=0.1, color='k')
plt.axvspan(4.5,5.5, alpha=0.1, color='k')



plt.show()