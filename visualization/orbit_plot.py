###################################################################
#   ORBITAL PLOTTING ROUTINE
###################################################################

import numpy as N
import astropy.units as U

from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.patches import Arc
from matplotlib.patches import Ellipse
from matplotlib.patches import FancyArrow
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Wedge

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=16)

import colorpy.colormodels as cm
import colorpy.blackbody as bb

from .style.colors import *

import datetime

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class OrbitPlot:

    def __init__(self, planet, time_resolution=None):
        self.planet = planet
        P = planet.P
        e = planet.e
        a = planet.a
        P_PSR = planet.pseudosynchronous_period()
        
        self.data = planet.data
        
        if time_resolution==None:
            self.time_resolution = planet.time_resolution
        else:
            self.time_resolution = time_resolution
            
        self.times = N.linspace(-0.5*P, 0.5*P, num=self.time_resolution+1)
        self.anomaly = planet.anomaly(self.times)

        self.PSR_angle = 2*N.pi*U.rad * N.array(self.times / P_PSR) + planet.w
        self.subsolar_angle = planet.w + self.anomaly['true']
        self.lum_angle = N.array((self.anomaly['true'].to(U.deg)-0*U.deg)/U.deg)

        #Draw the dots outlining the orbit using eccentric and true anomalies.
        E = self.anomaly['ecc']
        f = self.anomaly['true']
        cos_E = N.cos(E)
        cos_f = N.cos(f + planet.w + 180*U.deg)
        sin_f = N.sin(f + planet.w + 180*U.deg)

        #Create the coordinate arrays, in units of AU.
        self.orbit_scale = a / U.AU
        self.x_points = self.orbit_scale * (1 - e*cos_E) * cos_f
        self.y_points = self.orbit_scale * (1 - e*cos_E) * sin_f

    def draw_static(self, axis=plt.gca(), use_data=True, filetypes=['eps', 'pdf', 'png'], save=False):
        self.axis = axis

        #The orbital view. Draw points, subsolar longitude line, and labels.
        P = self.planet.P
        e = self.planet.e
        a = (self.planet.a / U.AU).value

        #If use_data is True, the dots showing the equal time intervals over the orbit will be colored according to where the data in each band exist. Note: this works best when the time sampling of the data is finer than the plotting time resolution.
        if use_data:
            coverage = N.zeros(self.time_resolution+1, dtype=bool)
            partial_phase = True
            t_set = {}
            for i, band in enumerate(['3p6','4p5','8p0']):
                if band in self.data:
                    t_set[band] = N.unique(self.data[band]['t'])
                    data = t_set[band]
                    band_range = N.array([N.any(N.abs(t-data) <= 0.5*P/U.d/self.time_resolution) for t in self.times/U.d])
                    coverage = coverage | band_range
                    #The "stitch" is that the array for checking coverage has its first element overlapping with the last. So if either are covered, both should be covered.
                    coverage_stitch = coverage[0] | coverage[-1]
                    coverage[0] = coverage_stitch
                    coverage[-1] = coverage_stitch
                    #Boolean mask for whether we consider the photometry a partial, rather than full, phase curve. I chose to call any light curve with a time range < 90% of the orbital period as "partial".
                    partial_phase = partial_phase & ((N.max(t_set[band])-N.min(t_set[band]))*U.d/P < 0.8)
                    orbit_range = axis.scatter(self.x_points[band_range], self.y_points[band_range], color=color_datlab[band], s=(6*(i+1)**2), zorder=(3-i)*5)
            orbit_gap = axis.scatter(self.x_points[~coverage], self.y_points[~coverage], color='k', s=2, zorder=100)

        else:
            orbit_points = axis.scatter(self.x_points, self.y_points, color='k', s=2, zorder=1)

        peri_angle = ((self.planet.w + 180*U.deg)%(360*U.deg)).to(U.deg)
        orbit_outline = self.axis.add_artist(Ellipse(xy = (a*e*N.cos(peri_angle+180*U.deg), a*e*N.sin(peri_angle+180*U.deg)),
        #orbit_outline = self.axis.add_artist(Ellipse(xy = (a*e*N.cos(peri_angle+0*U.deg), a*e*N.sin(peri_angle+0*U.deg)),
                                                     width = 2*a,
                                                     height = 2*a*N.sqrt(1-e**2),
                                                     angle = peri_angle.value,
                                                     fill = False, edgecolor = 'k', alpha = 0.3, linestyle = (0, (0.5, 1.5)) ))

        #The star point scales with the radius of the actual star, and is set true to the scale of the plotted orbit semimajor axis.
        star_scale = float(self.planet.R/U.AU)
        #The color is calculated from the blackbody color at the effective temperature.
        star_color = cm.irgb_string_from_xyz(bb.blackbody_color(float(self.planet.Teff/U.K)))
        limb_darkening_span = N.squeeze(N.array([colors.to_rgba_array(cm.irgb_string_from_xyz(bb.blackbody_color(T))) for T in N.linspace(0.8, 1, num=20)*self.planet.Teff/U.K]))

        #star_point = axis.scatter(0, 0, color=star_color, s=star_scale, edgecolors='k')
        for n, shade in enumerate(limb_darkening_span[:,:-1]):
            self.star_point = axis.add_artist(Ellipse(xy=[0,0], width=(star_scale*2)*(1-n/20), height=(star_scale*2)*(1-n/20), angle=0))
            self.star_point.set_facecolor(shade)

        #Draw the line from star to apastron, and the line to periastron.
        if e > 0:
            x_apo = N.cos(self.planet.w) * N.array([star_scale, self.orbit_scale*(1+e)])
            y_apo = N.sin(self.planet.w) * N.array([star_scale, self.orbit_scale*(1+e)])
            apo_line = axis.plot(x_apo, y_apo, color='#444444', linestyle='-', linewidth=1, dashes=[2,4])

            x_peri = N.cos(180*U.deg+self.planet.w) * N.array([star_scale, self.orbit_scale*(1-e)])
            y_peri = N.sin(180*U.deg+self.planet.w) * N.array([star_scale, self.orbit_scale*(1-e)])
            peri_line = axis.plot(x_peri, y_peri, color='#444444', linestyle='-', linewidth=1, dashes=[1,2])

            line_of_nodes = axis.axhline(y=0, linewidth=0.5, color='k', alpha=0.5, linestyle=(0, (0.25, 1.75)), zorder=1)

        #Draw the band representing the geometric transit and eclipse window.
        occultation_band = axis.axvspan(xmin=-star_scale, xmax=star_scale, facecolor='0.2', alpha=0.25)

        #Draw some motion arrows following a slightly larger elliptical arc around the orbit points.
        offset_factor = 0.1
        motion_arcs = [ Arc(xy = (a*e*N.cos(peri_angle+180*U.deg), a*e*N.sin(peri_angle+180*U.deg)),
                                                width = (1+offset_factor)*(2*a),
                                                height = (1+offset_factor)*(2*a*N.sqrt(1-e**2)),
                                                angle = peri_angle.value,
                                                theta1 = i,
                                                theta2 = j ) for (i,j) in [(-20,-10), (160, 170)] ]
        motion_paths = [(motion_arc.get_transform()).transform_path(motion_arc.get_path()) for motion_arc in motion_arcs]
        motion_lines = [axis.add_artist(FancyArrowPatch(path=motion_path, arrowstyle='->', mutation_scale=10, color='#444444')) for motion_path in motion_paths]

        #Remove axis labels and ticks.
        plt.setp(plt.getp(axis.axes, 'xticklabels'), visible=False)
        plt.setp(plt.getp(axis.axes, 'yticklabels'), visible=False)
        axis.set_aspect('equal')

        #Tight boundaries for the figure, using the span of the orbit.
        if use_data:
            if partial_phase:
                orbit_spanx = [N.min(self.x_points[coverage]) - 0.1*(N.max(self.x_points[coverage])-N.min(self.x_points[coverage])), N.max(self.x_points[coverage]) + 0.1*(N.max(self.x_points[coverage])-N.min(self.x_points[coverage]))]
                orbit_spany = [N.min(self.y_points[coverage]) - 0.1*(N.max(self.y_points[coverage])-N.min(self.y_points[coverage])), N.max(self.y_points[coverage]) + 0.1*(N.max(self.y_points[coverage])-N.min(self.y_points[coverage]))]
                axis.set_xlim(orbit_spanx)
                axis.set_ylim(orbit_spany)

            else:
                orbit_spanx = [N.min(self.x_points) - 0.1*(N.max(self.x_points)-N.min(self.x_points)), N.max(self.x_points) + 0.1*(N.max(self.x_points)-N.min(self.x_points))]
                orbit_spany = [N.min(self.y_points) - 0.1*(N.max(self.y_points)-N.min(self.y_points)), N.max(self.y_points) + 0.1*(N.max(self.y_points)-N.min(self.y_points))]
                axis.set_xlim(orbit_spanx)
                axis.set_ylim(orbit_spany)

        #The observer line points from the star to the bottom of the figure, which represents the direction to the observer.
        observer_linelength = self.orbit_scale * (1-e**2) / (1 + e*N.cos(self.planet.w+90*U.deg))
        xmin, xmax = axis.get_xlim()
        ymin, ymax = axis.get_ylim()
        observer_line = axis.arrow(xmin+0.05*(xmax-xmin), ymin+0.16*(ymax-ymin), 0, -0.05*(ymax-ymin), color='#444444', head_width=0.015*(ymax-ymin), linewidth=2, width=0.001*(ymax-ymin))
        earth_label = axis.text(xmin+0.05*(xmax-xmin), ymin+0.08*(ymax-ymin),'$\mathbf{\oplus}$', horizontalalignment='center', verticalalignment='top', fontproperties=fontprops)

        if e > 0:
            PSR_x = self.x_points[1:-1]
            PSR_y = self.y_points[1:-1]
            PSR_lines = [axis.add_artist(FancyArrow(x, y, 0.03*(ymax-ymin)*N.cos(PSR_line), 0.03*(ymax-ymin)*N.sin(PSR_line), color='black', lw=0.25, head_length=0.001, shape='right')) for (x, y, PSR_line) in zip(PSR_x, PSR_y, self.PSR_angle[1:-1])]
            subsolar_lines = [axis.plot([x, x+0.045*(ymax-ymin)*N.cos(subsolar_line)], [y, y+0.045*(ymax-ymin)*N.sin(subsolar_line)], color='#4E89B6', lw=1) for (x, y, subsolar_line) in zip(self.x_points[:-1], self.y_points[:-1], self.subsolar_angle[:-1])]

        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off

        plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left='off',      # ticks along the bottom edge are off
        right='off',         # ticks along the top edge are off
        labelleft='off') # labels along the bottom edge are off

        #Instead of explicit axes we opt for a scale bar, here chosen to be 0.01 AU.
        from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
        scalebar = AnchoredSizeBar(axis.transData,
                                   0.01,
                                   r"0.01 AU",
                                   loc=1,
                                   pad=0.0, borderpad=0.5, sep=5,
                                   frameon=False,
                                   fontproperties=fontprops)
        axis.add_artist(scalebar)

        #If we're plotting the whole orbit, the bounding box for the text showing the orbital period should be able to find room in the lower-right corner. Otherwise, there is a chance something could intersect it, and so we draw a white bounding box around it.
        if partial_phase: textbox_settings = dict(pad=0, fc='white', ec='white')
        else: textbox_settings = dict(pad=0, alpha=0)
        axis.text(0.975, 0.05,'$P = $ {0:.2f}'.format(self.planet.P), horizontalalignment='right',
                       verticalalignment='center',
                       transform=axis.transAxes,
                       fontproperties=fontprops,
                       bbox = textbox_settings)

        if save:
            for filetype in filetypes:
                plt.savefig('{0}_{1}_orbit.{2}'.format(str(datetime.date.today()).replace(" ", ""), (self.planet.name).replace(" ", ""), filetype), bbox_inches='tight', transparent=True)

    def set_rotation_period(self, rotation_period):
        self.rotation_period = rotation_period

    def set_num_orbits(self, num_orbits):
        self.num_orbits = num_orbits

    def draw_animate(self, axis=None):
        if axis != None:
            self.axis = axis
        self.draw_static(self.axis)
        times = N.linspace(0.*self.planet.P, self.num_orbits*self.planet.P, num=(self.num_orbits*self.time_resolution)+1)
        
        rotation_period = self.rotation_period
        self.line_rot = 2*N.pi*U.rad * N.array(times / rotation_period) - 0.5*N.pi*U.rad
        #This extra term is optional and allows the drawn line to point toward the star at the first periastron, to show the pseudosynchronicity of pseudosynchronous rotation.
        #self.line_rot += self.planet.subsolar_longitude(self.planet.times[0], rotation_period) - N.pi*U.rad * (self.planet.P / rotation_period) 
        
        #Draw a planet with one hemisphere illuminated. A line sticking out of the planet shows the rotation.
        self.planet_day = axis.add_artist(Wedge((self.x_points[0], self.y_points[0]), self.orbit_scale/10, self.lum_angle[0], self.lum_angle[0]+180, fc='w'))
        self.planet_night = axis.add_artist(Wedge((self.x_points[0], self.y_points[0]), self.orbit_scale/10, self.lum_angle[0]+180, self.lum_angle[0], fc='k'))
        
        self.line, = axis.plot([self.x_points[0], self.x_points[0] + 0.2*self.orbit_scale*N.cos(self.line_rot[0])], [self.y_points[0], self.y_points[0] + 0.2*self.orbit_scale*N.sin(self.line_rot[0])], color='black', lw=1)
            
    def update(self, i):
        j = i % self.time_resolution
        
        self.planet_day.remove()
        self.planet_night.remove()
        self.planet_day = self.axis.add_artist(Wedge((self.x_points[j], self.y_points[j]), self.orbit_scale/10, self.lum_angle[j], self.lum_angle[j]+180, fc='w'))
        self.planet_night = self.axis.add_artist(Wedge((self.x_points[j], self.y_points[j]), self.orbit_scale/10, self.lum_angle[j]+180, self.lum_angle[j], fc='k'))

        self.line.set_data([self.x_points[j], self.x_points[j] + 0.2*self.orbit_scale*N.cos(self.line_rot[i])], [self.y_points[j], self.y_points[j] + 0.2*self.orbit_scale*N.sin(self.line_rot[i])])
