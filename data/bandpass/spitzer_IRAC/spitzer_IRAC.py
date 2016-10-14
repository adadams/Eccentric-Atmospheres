###################################################################
#   SPITZER IRAC DIRECTORY FILE
###################################################################

import numpy as N
from numpy.lib import recfunctions as rfn
import astropy.units as U

#Pull the file data into a dictionary.
bandpass = {'3p6': N.genfromtxt('data/bandpass/spitzer_IRAC/average_3p6micron_response.txt', skip_header=2, names=True, comments='#'),
           '4p5': N.genfromtxt('data/bandpass/spitzer_IRAC/average_4p5micron_response.txt', skip_header=2, names=True, comments='#'),
           '5p8': N.genfromtxt('data/bandpass/spitzer_IRAC/average_5p8micron_response.txt', skip_header=2, names=True, comments='#'),
           '8p0': N.genfromtxt('data/bandpass/spitzer_IRAC/average_8p0micron_response.txt', skip_header=2, names=True, comments='#')}

#Add an extra column in the dictionary data for the weighted response by the energy at each wavelength.
for band in bandpass:
    bandpass[band] = rfn.append_fields(bandpass[band], names='weighted_spectrum', data=bandpass[band]['response'] / U.um.to(U.eV, bandpass[band]['wavelength'], equivalencies=U.spectral()), usemask=False)
