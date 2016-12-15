import astropy.constants as C
import astropy.units as U
import numpy as N
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.animation as anim
from mpl_toolkits.basemap import Basemap
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

class Planet:
    
    ###########################################################
    # Initial function that takes in a dictionary of established orbital properties.
    ###########################################################   
    def __init__(self, properties):

        self.name = properties['name']
        self.P = properties['orbital period']
        self.a = properties['semimajor axis']
        self.e = properties['eccentricity']
        self.w = properties['argument of periastron']
        self.mp = properties['planet mass']
        self.rp = properties['planet radius']
        self.M = properties['stellar mass']
        self.R = properties['stellar radius']
        self.Teff = properties['stellar temperature']

    ###########################################################
    # Set the spatial and time resolution for the model, and the number of orbits to run. Degrees for longitude and latitude is default.
    ###########################################################   
    def set_resolution(self,
                       longitude_resolution=360,
                       latitude_resolution=180,
                       time_resolution=1000,
                       num_orbits=1,
                       degrees=True):

        self.longitude_resolution = longitude_resolution
        self.latitude_resolution = latitude_resolution

        if degrees:
            self.phi_range = N.linspace(start=0, stop=360, num=longitude_resolution+1) * U.deg
            self.theta_range = N.linspace(start=-90, stop=90, num=latitude_resolution+1) * U.deg

        else:
            self.phi_range = N.linspace(start=0, stop=2*N.pi, num=longitude_resolution+1) * U.rad
            self.theta_range = N.linspace(start=-0.5*N.pi, stop=0.5*N.pi, num=latitude_resolution+1) * U.rad

        self.phis, self.thetas = N.meshgrid(self.phi_range, self.theta_range)

        self.times = N.linspace(start=0.5, stop=num_orbits+0.5, num=time_resolution*num_orbits)*self.P
        self.time_resolution = time_resolution
        self.num_orbits = num_orbits
         

    ###########################################################
    # Given a time and a periastron time in units of the orbital period, this function returns a dictionary of the eccentric and true anomalies.
    # Uses Kepler's equation to find the true anomaly numerically from the eccentric anomaly.
    ###########################################################       
    def anomaly(self, time, periastron_time=None):

        if periastron_time == None: self.t_per = 0.*self.P
        else: self.t_per = periastron_time

        E = 0. * U.rad
        E0 = N.pi * U.rad
        M = ((2*N.pi/self.P) * (time-self.t_per)) % (2*N.pi) * U.rad

        tolerance = 1e-4

        i = 0
        while((N.abs(E-E0)/E0).all > tolerance and i < 99):
            E0 = E
            M0 = E0 - self.e*N.sin(E0)*U.rad
            E = E0 + (M-M0)
            i += 1

        f = 2 * N.angle(N.sqrt(1-self.e)*N.cos(E/2.) + 1j*N.sqrt(1+self.e)*N.sin(E/2.))*U.rad
        return {'ecc': E, 'true': f}

    ###########################################################
    # Take the bandpass ratios and generate a flux ratio timeseries.
    ########################################################### 
    def flux_curve(self, spitzer_dataframe):
        mid = len(spitzer_dataframe.index)/2
        low_wave = N.average(spitzer_dataframe['wavelength'][0:mid], weights=spitzer_dataframe['weighted_spectrum'][0:mid])
        low_resp = N.average(spitzer_dataframe['weighted_spectrum'][0:mid])
        high_wave = N.average(spitzer_dataframe['wavelength'][mid:], weights=spitzer_dataframe['weighted_spectrum'][mid:])
        high_resp = N.average(spitzer_dataframe['weighted_spectrum'][mid:])
        spectral_response = N.array([low_resp, high_resp]) / U.eV
        bandpass_ratio = self.bandpass_ratios(N.array([low_wave, high_wave])*U.um, spectral_response)
        
        #spectral_response = N.array(spitzer_3p6_data['weighted_spectrum']) / U.eV
        #bandpass_ratio = HAT_P_2b.bandpass_ratios(N.array(spitzer_3p6_data['wavelength'])*U.um, spectral_response)
        return bandpass_ratio

    ###########################################################
    # Returns an array of booleans corresponding to whether the planet is in transit or secondary eclipse, as well as the times of transit and eclipse midpoints.
    ########################################################### 
    def calculate_occultation(self, times):
        #times=self.times

        #Use the treatment described on the NASA Exoplanet Archive from Greg Laughlin to calculate time of central transit and duration.
        f = 0.5*N.pi*U.rad - self.w
        E = lambda x: 2 * N.arctan(N.sqrt((1-self.e)/(1+self.e)) * N.tan(0.5*x))
        E_mid = E(f)%(2*N.pi*U.rad)
        Te_t = self.P/(2*N.pi) * (E_mid/U.rad-self.e*N.sin(E_mid))
        df = N.arcsin((self.R+self.rp)/self.a / (1-self.e*N.cos(E_mid)))
        E_end = E(f+df)%(2*N.pi*U.rad)
        half_transit_time = self.P/(2*N.pi) * (E_end/U.rad-self.e*N.sin(E_end)) - Te_t

        transit = N.abs((times - Te_t + 0.5*half_transit_time)%self.P) < half_transit_time

        #Now for the secondary eclipse, which is at the opposite anomaly.
        f = f + N.pi*U.rad
        E_mid = E(f)%(2*N.pi*U.rad)
        Te_e = self.P/(2*N.pi) * (E_mid/U.rad-self.e*N.sin(E_mid))
        df = N.arcsin((self.R+self.rp)/self.a / (1-self.e*N.cos(E_mid)))
        E_end = E(f+df)%(2*N.pi*U.rad)
        half_eclipse_time = self.P/(2*N.pi) * (E_end/U.rad-self.e*N.sin(E_end)) - Te_e

        eclipse = N.abs((times - Te_e + 0.5*half_eclipse_time)%self.P) < half_eclipse_time
        
        return {'transit flag': transit, 'transit': Te_t, 'eclipse flag': eclipse, 'eclipse': Te_e}

    ###########################################################
    # Calculates the pseudosynchronous rotation period.
    ########################################################### 
    def pseudosynchronous_period(self):
        
        omega_ratio = (1. + 7.5*self.e**2 + 45./8.*self.e**4 + 5./16.*self.e**6) / (1. + 3.*self.e**2 + 3./8.*self.e**4) * (1 - self.e**2)**(-1.5)
        rotation_period = self.P/omega_ratio
        
        return rotation_period
    
    ###########################################################
    # Returns the subsolar longitude on the planet at any time, given a rotation period.
    # Zero longitude is defined to be the longitude which is antisolar at the model initialization (first apastron).
    ########################################################### 
    def subsolar_longitude(self, time, rotation_period=None):

        if rotation_period == None: rot_per = self.p_rot
        else: rot_per = rotation_period
        sub_long = (((self.anomaly(time=time)['true']/U.rad).decompose() - 2*N.pi*(time/rot_per).decompose()) % (2*N.pi))*U.rad

        return sub_long
