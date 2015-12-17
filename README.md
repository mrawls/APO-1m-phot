# APO-1m-phot
Tools for dealing with observations using the Apache Point Observatory 1m robotically controlled telescope. Most of these tools are specifically for multiband photometry of eclipsing binaries. You need PyRAF, and patience.

*eclipsefinder.py* - Takes a text infile and makes three things: (1) a plot of when EB eclipses will be, (2) a text file listing these eclipses, and (3) a text file formatted for use by the robotically controlled 1m telescope at APO.

*img_inventory.py* - Looks at a directory tree of APO 1m FITS files, which must be organized by date (default status) and bunzip2-ed (not default status), and figures out WTF they are. You need to feed it a list of target coordinates to search the image headers for, a la `RGEB_info_alpha.txt`. It is slow and spits out a giant text file like `imginventory_list.txt` and a plot. The text file has columns for target name, date of observation, eclipse status (no/pri/sec), filter, and image filename.

*imagereduce.py* - Does overscan correction and flat field division on target images. (First it overscan-corrects and combines the flats, of course.) It works in two stages: (1) flat processing, then (2) flat inventorying + applying processed flats to target images. Uses a slightly tweaked version of `imginventory_list.txt` (necessary tweaks are comments at the top of the code; read them). Puts each target's reduced images in a separate folder. Sometimes you have to run it more than once to finish stage (2) because IRAF gets cranky.

*photplotter.py* - After you align the reduced images and run PHOT in IRAF/PyRAF for your target and some comparison stars, you will have a bunch of `*.mag.*` text files with useful information. This program turns those files into differential photometry. For now, this program assumes there is one target star and six companion stars, each processed with four possible apertures. It makes plots of the comparison star magnitudes (so you can make sure they don't vary), and returns differential magnitudes vs. time (and phase, assuming you're dealing with an eclipsing binary and define period and BJD0 appropriately!).

*photbinner.py* - Once you have differential magnitudes, you may have a bunch on one night and then none for a while and then lots again, etc. This program bins a time series of differential magnitudes with a delta-t threshold you specify.
