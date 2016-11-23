###################################################################
#   SURFACE MAP PLOTTING ROUTINE
###################################################################

import numpy as N
import astropy.units as U

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class SurfacePlot:

    def __init__(self, planet_object, map_array, rotation_period):

        self.planet_object = planet_object
        self.rotation_period = rotation_period
        self.phi = planet_object.phis
        self.theta = planet_object.thetas
        self.times = planet_object.times
        self.time_resolution = planet_object.time_resolution

        self.map_array = map_array
        self.map_min = N.nanmin(self.map_array) / self.map_array.unit
        self.map_max = N.nanmax(self.map_array) / self.map_array.unit

        self.long_obs = ((1.5*N.pi + (2*N.pi*(self.times/self.rotation_period).decompose() - ((self.planet_object.w)/U.rad).decompose())) % (2*N.pi))*U.rad

    def draw(self, axis=plt.gca()):
        #Direct focus to the desired axis.
        self.axis = axis
        plt.sca(self.axis)

        #The container for the projection using Basemap.
        #self.projected_map = Basemap(projection='hammer', resolution = 'h', lat_0=0, lon_0=0)
        self.projected_map = Basemap(projection='ortho', resolution = 'h', lat_0=45, lon_0=0)

        #Make sure the array you provide as an argument is unitless.
        lons = N.around((self.phi - self.long_obs[0])/U.deg) % 360
        lats = self.theta / U.deg
        self.globemap = self.projected_map.pcolor(lons, lats, self.map_array[0]/self.map_array.unit, cmap=plt.cm.magma, vmin=self.map_min, vmax=self.map_max, latlon=True)

        #Options to draw lat/lon lines.
        #self.projected_map.drawmeridians(N.arange(-180,180,10))
        #self.projected_map.drawparallels(N.arange(-90,90,30))

        #The temperature legend for the color map.
        cb = self.projected_map.colorbar(self.globemap, "bottom", size="5%", pad="2%")
        cb.set_label('Temperature (K)', size = 14)
        plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), rotation = 90, size = 12)

        #The subsolar longitude on the planet and the line object for that longitude to be drawn on the map. Also draw longitude lines bounding the hemisphere visible to the observer at any given time.
        ss_long = ((self.planet_object.subsolar_longitude(self.times[0], rotation_period=self.rotation_period) - self.long_obs[0])/U.deg).decompose()
        #self.drawn_lines = self.projected_map.drawmeridians([ss_long, long_obs/U.deg-90, long_obs/U.deg+90], linewidth=3, dashes=[3,1])
        self.drawn_lines = self.projected_map.drawmeridians([ss_long], linewidth=3, dashes=[3,1])
            
    def update(self, i):
        print i
        #Remove the existing color map.
        for c in self.globemap.findobj(): c.remove()

        #Redraw the color map for the given time step.
        lons = N.around((self.phi - self.long_obs[i])/U.deg) % 360
        lats = self.theta / U.deg
        self.globemap = self.projected_map.pcolor(lons, lats, self.map_array[i]/self.map_array.unit, cmap=plt.cm.magma, vmin=self.map_min, vmax=self.map_max, latlon=True)
        
        #Remove the existing subsolar longitude line. The data are stored in a dictionary.
        self.ss_keys = (self.drawn_lines).keys()
        for key in self.ss_keys: del self.drawn_lines[key]

        #Redraw the subsolar longitude for the given time step. Also draw longitude lines bounding the hemisphere visible to the observer at any given time.
        ss_long = ((self.planet_object.subsolar_longitude(self.times[i], rotation_period=self.rotation_period) - self.long_obs[i])/U.deg).decompose()
        #long_obs = ((1.5*N.pi + (2*N.pi*(self.times[i]/self.rotation_period).decompose() - ((self.planet_object.w)/U.rad).decompose())) % (2*N.pi))*U.rad
        #self.drawn_lines = self.projected_map.drawmeridians([ss_long, long_obs/U.deg-90, long_obs/U.deg+90], linewidth=3, dashes=[3,1])
        self.drawn_lines = self.projected_map.drawmeridians([ss_long], linewidth=3, dashes=[3,1])
