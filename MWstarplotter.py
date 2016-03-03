from __future__ import print_function
import numpy as np
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
infile = '../../RGEB_distinfo.txt'
target_col = 0
RA_col = 3
Dec_col = 4
dist_col = 6
derr_col = 7
usecols = (target_col, RA_col, Dec_col, dist_col, derr_col)

# Read in target information from a text file
targets, RAs, Decs, dists, derrs = np.loadtxt(infile, comments='#', usecols=usecols, 
    dtype={'names': ('targets', 'RAs', 'Decs', 'dists', 'derrs'),
    'formats': (np.int, '|S11','|S11', np.float64, np.float64)}, unpack=True)

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
star_galcens = []
for star in starlocs:
    if star.distance > 0:
        star_galcens.append(star.transform_to(coord.Galactocentric))

#print(star_galcens)

axnew1 = fig2.add_subplot(1,1,1, projection='3d', aspect='equal')
axnew1.grid(False)
#axnew2 = fig2.add_subplot(2,2,2)
#axnew3 = fig2.add_subplot(2,2,3)
for star in star_galcens:
    #print(star.x, star.y, star.z)
    axnew1.scatter(star.x, star.y, star.z, c='r', edgecolors='r', s=50)
    #axnew1.scatter(star.x, star.y)
    #axnew2.scatter(star.x, star.z)
    #axnew3.scatter(star.y, star.z)
#axnew1.set_xlim(-10000, 1000)
#axnew1.set_ylim(-1000, 5000)
#axnew1.set_zlim(-1000, 1000)
#axnew1.set_xlabel('X (pc)')
#axnew1.set_ylabel('Y (pc)')
#axnew1.set_zlabel('Z (pc)')
axnew1.scatter(0, 0, 0, marker='o', c='k', edgecolors='k', s=50) # galactic center
axnew1.scatter(-8300, 0, 27, marker='*', c='k', edgecolors='k', s=80) # Sun

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

# Attempt to plot an image of the Milky Way in the X-Y plane?
#img = mpimg.imread('../../MWimage.png')
#stretch = 1.
#ximg, yimg = np.ogrid[-img.shape[0]/2.*stretch:img.shape[0]/2.*stretch, -img.shape[1]/2.*stretch:img.shape[1]/2.*stretch]
#axnew1.plot_surface(ximg, yimg, 0, rstride=100000, cstride=100000, facecolors=img)
##axnew1.imshow(img)

plt.show()