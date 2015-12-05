from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import glob
'''
Meredith Rawls, 2015

So you have a bunch of FITS images in a directory.
Let's make a massive PDF with postage stamp preview images for all of them!
There will be pics_per_page images per page.

OSX pro tip: if you get 'IOError: [Errno 24] Too many open files' then you should
try setting 'ulimit -n VALUE' where VALUE is a large number... maybe 4096?
But if you make it too large, your computer will run out of RAM and cry.
(Because you're making it work hard to open a zillion FITS images.)

This program takes a while to run if you have lots of images, which you probably do.
'''
# MODIFY THE FOLLOWING AS DESIRED:
datadir = '../../1m_photometry/KIC09246715/'
outfile = datadir + datadir[-12:-1] + '_fitspreview.pdf' 
imagelist = glob.glob(datadir + '*.fits')
pics_per_page = 20

def plot_maker(imagefilelist):
    '''
    Takes a list of FITS files which must be <= pics_per_page long
    Creates a PDF page of images containing up to pics_per_page images
    '''
    #pdf_page = PdfPages(imagefile+'.pdf') # OPTION: create individual PDF page
    fig = plt.figure()
    pagecols = 4
    pagerows = 5
    for idx, imagefile in enumerate(imagefilelist):
        hdu = fits.open(imagefile, ignore_missing_end = True)
        if 'B' == imagefile[-6]: color = 'b'
        elif 'V' == imagefile[-6]: color = 'g'
        elif 'R' == imagefile[-6]: color = 'r'
        elif 'I' == imagefile[-6]: color = 'k'
        else: color = 'k'
        image_data = hdu[0].data
        ax = fig.add_subplot(pagerows, pagecols, idx+1)
        ax.set_xticklabels(())
        ax.set_yticklabels(())
        plt.axis([48,2096,0,2048])
        plt.imshow(image_data, cmap = 'gray_r', norm = LogNorm())
        # NOTE this normalization isn't always perfect... may need tweaking
        plt.grid(True)
        plt.title(imagefile[-20:], size=10, color=color)
        #plt.colorbar() # uncomment if you want to see the colorbar
        plt.subplots_adjust(wspace=0.05, hspace=0.12)
        hdu.close()
    #pdf_page.savefig() # OPTION: save this individual PDF page
    #pdf_page.close() # OPTION: for individual PDF page only
    pdf_full.savefig() # add the individual page to the giant PDF
    
pdf_full = PdfPages(outfile)
nimages = len(imagelist)
npages = int([np.rint(nimages/pics_per_page) if (np.float(nimages)/pics_per_page)%pics_per_page == 0 else np.rint(nimages/pics_per_page)+1][0])
print('Generating', len(imagelist), 'images on', npages, 'pages...')
for page in range(0, npages): # ADJUST TO UPPER VALUE < npages FOR TESTING!!!
    try:
        if pics_per_page*(page+1) <= len(imagelist):
            plot_maker(imagelist[pics_per_page*page:pics_per_page*(page+1)])
        else: # last page
            plot_maker(imagelist[pics_per_page*page:])
    except:
        print('Error encountered, sad kittens')
        break
    print('Page {0} of {1} done'.format(page+1, npages))
pdf_full.close()
