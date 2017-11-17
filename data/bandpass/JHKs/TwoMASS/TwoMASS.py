###################################################################
#   2MASS DIRECTORY FILE
###################################################################

import numpy as N
from numpy.lib import recfunctions as rfn
import astropy.units as U

#Pull the file data into a dictionary.
bandpass = {'J': N.genfromtxt('data/bandpass/JHKs/TwoMASS/J.dat', names=True),
            'H': N.genfromtxt('data/bandpass/JHKs/TwoMASS/H.dat', names=True),
            'Ks': N.genfromtxt('data/bandpass/JHKs/TwoMASS/Ks.dat', names=True)}

#Add an extra column in the dictionary data for the weighted response by the energy at each wavelength.
for band in bandpass:
    bandpass[band] = rfn.append_fields(bandpass[band], names='weighted_spectrum', data=bandpass[band]['response'] / U.um.to(U.eV, bandpass[band]['wavelength'], equivalencies=U.spectral()), usemask=False)
