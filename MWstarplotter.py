from __future__ import print_function
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.image as mpimg
import mpl_toolkits.mplot3d.art3d as art3d
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.coordinates import Distance
'''
Reads in RA, Dec, and distance (with error bars) for a set of stars.
Makes plots of where those stars are in the galaxy.
'''
infile = 'RGEB_distinfo.txt'
target_col = 0
RA_col = 3
Dec_col = 4
dist_col = 6
derr_col = 7
FeH_col = 8
usecols = (target_col, RA_col, Dec_col, dist_col, derr_col, FeH_col)

# Read in target information from a text file
targets, RAs, Decs, dists, derrs, FeHs = np.loadtxt(infile, comments='#', usecols=usecols, 
    dtype={'names': ('targets', 'RAs', 'Decs', 'dists', 'derrs', 'FeHs'),
    'formats': (np.int, '|S11','|S11', np.float64, np.float64, np.float64)}, unpack=True)

# Put the RAs, Decs, and distances in a more useful format
RAs = coord.Angle(RAs, unit=u.hour)
RAs_plot = RAs.wrap_at(180*u.degree) # for reasons
Decs = coord.Angle(Decs, unit=u.degree)
dists = dists*u.pc
derrs = derrs*u.pc

# Define a SkyCoord object for each target
starlocs = []
for target, RA, Dec, dist, derr in zip(targets, RAs, Decs, dists, derrs):
    starlocs.append( SkyCoord(ra=RA, dec=Dec, distance=dist) )
    
#print(starlocs) #IT WORKS! (but no way to include distance uncertainty, I don't think)

# Plot the target locations on a sky plane projection? Maybe useful?
#fig = plt.figure(figsize=(8,6))
#ax = fig.add_subplot(111, projection='mollweide')
#ax.scatter(RAs_plot.radian, Decs.radian)
##ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
#ax.grid(True)
#plt.show()

# Make a figure
fig = plt.figure()

# First subplot: galactic (l,b) coordinates in an Aitoff projection
ax = fig.add_subplot(2,1,1, projection='aitoff')
for star in starlocs:
    if star.distance > 0: # only consider the targets with distance info
        #print(star.galactic)
        lrad = star.galactic.l.radian
        #if lrad > np.pi:
        #    lrad = lrad - 2.*np.pi
        brad = star.galactic.b.radian
        ax.scatter(lrad, brad)
ax.grid(True)
ax.set_title('Galactic (l,b)')
#ax.set_xlim(360., 0.)
#ax.set_ylim(-90., 90.)
#ax.set_xlabel('Galactic Longitude')
#ax.set_ylabel('Galactic Latitude')

# Second subplot: cartesian heliocentric coordinates in an X, Y slice
ax2 = fig.add_subplot(2,2,3, aspect='equal')
for star in starlocs:
    if star.distance > 0:
        xcart = star.cartesian.x
        ycart = star.cartesian.y
        zcart = star.cartesian.z
        ax2.scatter(xcart, ycart)
        #star.representation = 'cylindrical'
        #print(xcart, ycart, zcart)
ax2.set_xlim(-2500., 2500.)
ax2.set_ylim(-3000., 3000.)
ax2.set_xlabel('X (pc)')
ax2.set_ylabel('Y (pc)')
plt.plot(0, 0, marker='*', color='y', ms=20)

# Third subplot: cartesian heliocentric coordinates in an X, Z slice
ax3 = fig.add_subplot(2,2,4, aspect='equal')
ax3.set_xlim(-2500., 2500.) # really x
ax3.set_ylim(-3000., 3000.) # actually z
for star in starlocs:
    if star.distance > 0:
        xcart = star.cartesian.x
        ycart = star.cartesian.y
        zcart = star.cartesian.z
        #print(xcart, ycart, zcart)
        ax3.scatter(xcart, zcart)
ax3.set_xlabel('X (pc)')
ax3.set_ylabel('Z (pc)')
plt.plot(0, 0, marker='*', color='y', ms=20)

#plt.show()

# Make a second figure
fig2 = plt.figure()

# Transform stars to galactocentric coordinates (cartesian)
# Set a color scheme as a function of metallicity (different color for every 0.2 dex)
star_galcens = []
colorlist = []
for star, FeH in zip(starlocs, FeHs):
    if star.distance > 0:
        star_galcens.append(star.transform_to(coord.Galactocentric))
        if FeH < -0.8: color='#ffffb2' #yellowest
        elif FeH >= -0.8 and FeH < -0.6: color='#fed976'
        elif FeH >= -0.6 and FeH < -0.4: color='#feb24c'
        elif FeH >= -0.4 and FeH < -0.2: color='#fd8d3c'
        elif FeH >= -0.2 and FeH < 0.0: color='#fc4e2a'
        elif FeH >= 0.0 and FeH < 0.2: color='#e31a1c'
        elif FeH >= 0.2: color='#b10026' #reddest
        colorlist.append(color)

#print(star_galcens)

axnew1 = fig2.add_subplot(1,1,1, projection='3d', aspect='equal')
axnew1.set_axis_off() 
axnew1.grid(False)
axnew1.xaxis.set_ticklabels([])
axnew1.yaxis.set_ticklabels([])
axnew1.zaxis.set_ticklabels([])
axnew1.xaxis.set_ticks([])
axnew1.yaxis.set_ticks([])
axnew1.zaxis.set_ticks([])
for i, star in enumerate(star_galcens):
    #print(star.x, star.y, star.z)
    axnew1.scatter(star.x, star.y, star.z, c=colorlist[i], edgecolors='k', s=150)
axnew1.scatter(0, 0, 0, marker='o', c='k', edgecolors='k', s=50) # galactic center
axnew1.scatter(-8300, 0, 27, marker='*', c='k', edgecolors='k', s=150) # Sun

# Contour-type circles that radiate out from the galactic center for reference
circle1 = plt.Circle((0,0), 2000, color='0.75', fill=False)
circle2 = plt.Circle((0,0), 4000, color='0.75', fill=False)
circle3 = plt.Circle((0,0), 6000, color='0.75', fill=False)
circle4 = plt.Circle((0,0), 8000, color='0.75', fill=False)
circle5 = plt.Circle((0,0), 10000, color='0.75', fill=False)
axnew1.add_patch(circle1)
axnew1.add_patch(circle2)
axnew1.add_patch(circle3)
axnew1.add_patch(circle4)
axnew1.add_patch(circle5)
art3d.pathpatch_2d_to_3d(circle1, z=0, zdir='z')
art3d.pathpatch_2d_to_3d(circle2, z=0, zdir='z')
art3d.pathpatch_2d_to_3d(circle3, z=0, zdir='z')
art3d.pathpatch_2d_to_3d(circle4, z=0, zdir='z')
art3d.pathpatch_2d_to_3d(circle5, z=0, zdir='z')

# Colorbar key
axnew2 = fig2.add_subplot(12,1,10)
cmap = mpl.colors.ListedColormap(['#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c'])
cmap.set_over('#b10026') #reddest, high Fe/H
cmap.set_under('#ffffb2') #yellowest, low Fe/H
bounds = [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb = mpl.colorbar.ColorbarBase(axnew2, cmap=cmap, norm=norm, ticks=bounds, extend='both', 
    boundaries=[-1.0]+bounds+[0.4], spacing='proportional', orientation='horizontal')
cb.set_label('[Fe/H]', size=26)

# manually make a key?
#fig2.text(0.7, 0.7, 'Testing words here', ha='center', va='center', size=26)

# Attempt to plot an image of the Milky Way in the X-Y plane?
#img = mpimg.imread('../../MWimage.png')
#stretch = 1.
#ximg, yimg = np.ogrid[-img.shape[0]/2.*stretch:img.shape[0]/2.*stretch, -img.shape[1]/2.*stretch:img.shape[1]/2.*stretch]
#axnew1.plot_surface(ximg, yimg, 0, rstride=100000, cstride=100000, facecolors=img)
##axnew1.imshow(img)

fig3 = plt.figure()
ax3main = fig3.add_subplot(3,1,2, aspect='equal')
for i, star in enumerate(star_galcens):
    rkpc = np.sqrt(star.x*star.x + star.y*star.y)/1000.
    zkpc = star.z/1000.
    ax3main.scatter(rkpc, zkpc, c=colorlist[i], edgecolor='k', s=150)
#ax3main.scatter(0, 0, marker='o', c='k', edgecolors='k', s=50) # galactic center
ax3main.scatter(8.3, 0.027, marker='*', c='k', edgecolors='k', s=150) # Sun
ax3main.set_xlabel('Galactic radius $R$ (kpc)', size=26)
ax3main.set_ylabel('Height $z$ (kpc)', size=26)
ax3main.set_xlim(6, 9)
#ax3main.set_ylim(-0.1, 0.9)
#plt.xticks( (-8, -6, -4, -2, 0), ('8', '6', '4', '2', '0') )
ax3cb = fig3.add_subplot(15,1,13)
cb3 = mpl.colorbar.ColorbarBase(ax3cb, cmap=cmap, norm=norm, ticks=bounds, extend='both', 
    boundaries=[-1.0]+bounds+[0.4], spacing='proportional', orientation='horizontal')
cb3.set_label('[Fe/H]', size=26)

plt.show()