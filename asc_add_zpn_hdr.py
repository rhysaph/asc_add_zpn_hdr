#!/usr/bin/env python

# Add the ZPN projection WCS headers to an all sky image

# Example usage:
# ./asc_add_zpn_hdr.py image.fits
#
# This page was helpful
# https://docs.astropy.org/en/stable/wcs/example_create_imaging.html
# 
# Method: 
# Use zenith as reference point, we know the x, y pixel coordinates and we can 
# calculate the RA and Dec of Zenith given our location and date and time.

import sys
#import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
#from astropy.wcs import WCS # changed this to wcs
from astropy import wcs
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy import units as u
from astropy.time import Time
import math

# These are the changeable parameters
# zenith position
# Jamie and Queenette calculated this as  1598.05, 1058.79
xzenith = 1598.05
yzenith = 1058.79
print("Using these zenith coordinates ",xzenith, yzenith)

# Rotation angle (Jamie and Queenette calculated this as -12)
rot_angle = -12
print("Using rotation angle of ",rot_angle)

# print("Not setting ZPN distortion coefficients")

# Plate scale in arcseconds per pixel, we don't think this should be varied
xplate_scale = 0.055
yplate_scale = 0.055

# These probably should not be changed.
lonpole = 180.0
latpole = 0.0 # set to latitude later

# These are the altitude stretch coeficients found by
# 2023-24 students. They need to be divided by the plate
# scale before use.
A1 = 0.05648 / xplate_scale
A2 = 8.227e-06 / xplate_scale
A3 = -1.089e-08 / xplate_scale

print("Using ZPN these coefficients:")
print("A1 = ", A1)
print("A2 = ", A2)
print("A3 = ", A3)

########### Program starts here =====================
imagefile = sys.argv[1]
print("Opening ",imagefile)

hdu = fits.open(imagefile)
hdr = hdu[0].header
data_array = fits.getdata(imagefile,ext=0)
hdu.close()

#print(hdr) 

wcs_orig = wcs.WCS(hdr)
print("Original WCS is \n",wcs_orig)

print("\nCreating new file\n")

# Get some header values
obs_date = hdu[0].header['DATE']
obsgeob = hdu[0].header['OBSGEO-B']
obsgeol = hdu[0].header['OBSGEO-L']

# Create time object
obs_time = Time(obs_date, format='isot', scale='utc')
mjd = obs_time.mjd
print("Modified Julian Date (MJD) is ",mjd)

lat_degrees = float(obsgeob)
lon_degrees = float(obsgeol)

print("latitude is ",lat_degrees)
print("longitude is ",lon_degrees)

#timestr = '2020-10-26T01:30:01'
#lon_degrees =  -2.6021
#lat_degrees = 51.4588

# Change WCS
# Create a new WCS object.  The number of axes must be set
# from the start

#new_wcs = WCS(hdu)
new_wcs = wcs.WCS(naxis=2)

# Set FITS heders to say this is a  ZPN WCS header
new_wcs.wcs.ctype = ["RA---ZPN", "DEC--ZPN"]

# Using values from start of program
new_wcs.wcs.crpix = [xzenith,yzenith]

new_wcs.wcs.cdelt = [xplate_scale, yplate_scale]

# Need to calculate RA of zenith, which is the same as the LST.
#time_obj = Time(timestr, format='isot', scale='utc')

LST = obs_time.sidereal_time('apparent', lon_degrees)
print("LST is ",LST.degree)

# Set zenith coords as reference coords
new_wcs.wcs.crval = [LST.degree,lat_degrees]

# Set CD matrix here for a clockwise rotation of rot_angle degrees
a = math.sin(math.radians(rot_angle))
b = math.cos(math.radians(rot_angle))
new_wcs.wcs.pc=[(b,-a,), (a ,b )]

# Use these values from the start of the file
new_wcs.wcs.lonpole = lonpole
new_wcs.wcs.latpole = lat_degrees


# This sets the ZPN distortion parameters to no correction.
#new_wcs.wcs.set_pv([(2, 1, 1.0)])

# Set 3 orders of ZPN correction
new_wcs.wcs.set_pv([(2, 1, A1),
                    (2, 2, A2),
                    (2, 3, A3)
                ])

print("New WCS is ",new_wcs)

outfile = 'new_zpn_wcs.fits'
print("\nWriting new file", outfile)

new_header = new_wcs.to_header()
new_header.set('MJD-OBS',mjd,'MJD of Exposure Start' )
new_hdu = fits.PrimaryHDU(data_array,header=new_header)

#new_hdu.append('MJD-OBS',mjd)
#new_header[MJDREF] = mjd

new_hdu.writeto(outfile,overwrite=True)


#radec_coords = pixel_to_skycoord(xp, yp, new_wcs, origin=0, mode='all', cls=None)

#hdu2 = fits.open(outfile)
#hdr2 = hdu2[0].header
#hdu2.close()

#print("New header is \n",hdr2) 

#wcs2 = WCS(hdr2)
#print("New WCS is ",wcs2)

