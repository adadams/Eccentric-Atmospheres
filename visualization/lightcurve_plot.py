###################################################################
#   LIGHT CURVE PLOTTING ROUTINE
###################################################################

import numpy as N
import astropy.units as U

import matplotlib as mpl
from matplotlib import pyplot as plt

from style.colors import *

import datetime

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class LightCurvePlot:

    def __init__(self, planet_object, data, model, parameters):
        self.planet_object = planet_object
        self.data = data
        self.model = model
        self.parameters = parameters

    def draw(self, axis=plt.gca(), save=False):
        for band in self.data:
            low_cut = (self.planet_object.times - (self.planet_object.P*(self.planet_object.num_orbits))) / U.d >= N.min(self.data[band]['t'])
            high_cut = (self.planet_object.times - (self.planet_object.P*(self.planet_object.num_orbits))) / U.d <= N.max(self.data[band]['t'])
            plot_cut = low_cut & high_cut

            plotted_times = self.planet_object.times[plot_cut] - (self.planet_object.P*(self.planet_object.num_orbits))

            for sub_band in self.model:
                subband_label = sub_band.replace('p', '.')
                axis.plot(plotted_times, (self.model[sub_band]['model'].reshape(self.planet_object.time_resolution*self.planet_object.num_orbits))[plot_cut], color=color_modbg[sub_band], label=r'${0} \mu$m Model'.format(subband_label))

            band_label = band.replace('p', '.')
            axis.scatter(self.data[band]['t'], self.data[band]['flux'], color=color_datlab[band], label=r'${0} \mu$m Data'.format(band_label))

            axis.set_xlabel(r'$t$ (Days)')
            axis.set_ylabel(r'$\left(F/F_\star\right)_{3.6 \mu\mathrm{m}}$')
            axis.margins(0.05)
            axis.set_ylim([0.9995, 1.002])
            axis.legend(loc='upper center', fontsize=8, bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True, title=r'$P_{{\mathrm{{rot}}}}\approx$ {0:.1g}, $t_{{\mathrm{{rad}}}} \sim$ {1:.2f}, $T_n=$ {2:.0f}, $A=$ {3:.2f}'.format(*self.parameters[band]))
            
            if save:
                plt.tight_layout
                plt.savefig('{0}_{1}_{2}.pdf'.format(datetime.date.today(), self.planet_object.name, band))
                plt.cla()
