###################################################################
#   ORBITAL PLOTTING ROUTINE
###################################################################

import numpy as N
import astropy.units as U

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class OrbitPlot:

    def __init__(self, planet_object, rotation_period):
        self.x_points = N.zeros(planet_object.time_resolution)
        self.y_points = N.zeros(planet_object.time_resolution)
        
        self.planet_object = planet_object
        self.line_rot = 2*N.pi * N.array((planet_object.times-planet_object.P/2) / rotation_period)
        self.anomaly = planet_object.anomaly(planet_object.times[0:planet_object.time_resolution])
        self.lum_angle = N.array(self.anomaly['true'].to(U.deg)/U.deg - 90)

    def draw(self, axis=plt.gca()):
        self.axis = axis

        #Draw the dots outlining the orbit using eccentric and true anomalies.
        #for i, time in enumerate(self.planet_object.times[0:self.planet_object.time_resolution]):
        E = self.anomaly['ecc']
        f = self.anomaly['true']
        cos_E = N.cos(E)
        cos_f = N.cos(f)
        sin_f = N.sin(f)

        #Create the coordinate arrays, in units of AU.
        self.orbit_scale = self.planet_object.a / U.AU
        self.x_points = self.orbit_scale * (1 - self.planet_object.e*cos_E) * cos_f
        self.y_points = self.orbit_scale * (1 - self.planet_object.e*cos_E) * sin_f

        #The orbital view. Draw points, subsolar longitude line, and labels.
        orbit_points = self.axis.scatter(self.x_points, self.y_points, color='black', s=2)
        star_point = self.axis.scatter(0, 0, color='black', s=16)

        self.pl_w = self.axis.add_artist(Wedge((self.x_points[0], self.y_points[0]), self.orbit_scale/10, self.lum_angle[0]+180, self.lum_angle[0], fc='w'))
        self.pl_b = self.axis.add_artist(Wedge((self.x_points[0], self.y_points[0]), self.orbit_scale/10, self.lum_angle[0], self.lum_angle[0]+180, fc='k'))

        self.line, = self.axis.plot([self.x_points[0], self.x_points[0]*(1 + 0.1*N.cos(self.line_rot[0]))], [self.y_points[0], self.y_points[0]*(1 - 0.1*N.sin(self.line_rot[0]))], color='black', lw=1)

        #Axis labels.
        plt.setp(plt.getp(self.axis.axes, 'xticklabels'), rotation = 90, size = 12)
        plt.setp(plt.getp(self.axis.axes, 'yticklabels'), size = 12)
        self.axis.set_xlabel(r'$x$ (AU)', size=16)
        self.axis.set_ylabel(r'$y$ (AU)', size=16)
        self.axis.set_aspect('equal')
            
    def update(self, i):
        j = i%self.planet_object.time_resolution

        self.pl_w.remove()
        self.pl_b.remove()
        self.pl_w = self.axis.add_artist(Wedge((self.x_points[j], self.y_points[j]), self.orbit_scale/10, self.lum_angle[j]+180, self.lum_angle[j], fc='w'))
        self.pl_b = self.axis.add_artist(Wedge((self.x_points[j], self.y_points[j]), self.orbit_scale/10, self.lum_angle[j], self.lum_angle[j]+180, fc='k'))

        self.line.set_data([self.x_points[j], self.x_points[j] + 0.2*self.orbit_scale*N.cos(self.line_rot[i])], [self.y_points[j], self.y_points[j] + 0.2*self.orbit_scale*N.sin(self.line_rot[i])])
