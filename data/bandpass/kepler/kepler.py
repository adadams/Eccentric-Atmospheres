###################################################################
#   KEPLER DIRECTORY FILE
###################################################################

import numpy as N
from numpy.lib import recfunctions as rfn
import astropy.units as U

#Pull the file data into a dictionary.
bandpass = {'KP': N.genfromtxt('data/bandpass/kepler/kepler_response_hires1.txt', skip_header=6, names=True, comments='#')}

#The wavelengths in the file are in nm, so we'll convert them to um.
bandpass['KP']['wavelength'] /= 1000

#Add an extra column in the dictionary data for the weighted response by the energy at each wavelength.
for band in bandpass:
    bandpass[band] = rfn.append_fields(bandpass[band], names='weighted_spectrum', data=bandpass[band]['response'] / U.um.to(U.eV, bandpass[band]['wavelength'], equivalencies=U.spectral()), usemask=False)
