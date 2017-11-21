###################################################################
#   LIGHT CURVE PLOTTING ROUTINE
###################################################################

import numpy as N
import astropy.units as U

import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import interpolate as inter

from .style.colors import *

import datetime

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

class LightCurvePlot:

    def __init__(self, planet, data, model, parameters):
        self.planet = planet
        self.data = data
        self.model = model
        self.parameters = parameters

    def draw(self, axis, s, phase_overlap=0.25, filetypes=['eps', 'pdf', 'png'], combo=False, save=False):
        P = self.planet.P
        time_resolution = self.planet.time_resolution
        num_orbits = self.planet.num_orbits
        
        #We could plot up to 3 full orbits. The main (center) orbit will be from -0.5 to 0.5 times the period, with t=0 representing periastrion passage.
        plotted_times = N.linspace(-1.5, stop=1.5, num=3*time_resolution)*P
        #Time of mid-transit.
        transit_time = self.planet.calculate_occultation(self.planet.times)['transit']

        t_set = {}
        partial_phase = {}
        x = {}
        x_smooth = {}
        low_smoothed = {}
        high_smoothed = {}
        y_model = {}
        y_median = {}
        lowerlim = {}
        upperlim = {}

        for band in sorted(self.data):
            #Pull out upper and lower limits from carefully sampled data, which has multiple sampled points at each time.
            data_upper = []
            data_lower = []
            data_median = []

            #The data may be sampled with multiple fluxes at a given time. We want a list of the unique times in the data for our x-axis.
            t_set[band], t_index = N.unique(self.data[band]['t'], return_index=True)
            #Boolean mask for whether the data point is in transit or eclipse.
            #occulted = self.planet.calculate_occultation(t_set[band]*U.d)['eclipse flag'] | self.planet.calculate_occultation(t_set[band]*U.d)['transit flag']
            occulted = self.data[band]['occultation'][t_index] != b''
            #Boolean mask for whether we consider the photometry a partial, rather than full, phase curve. I chose to call any light curve with a time range < 90% of the orbital period as "partial".
            partial_phase[band] = (N.max(t_set[band])-N.min(t_set[band]))*U.d/P < 0.8
            
            #In order to create a contour representing the spread in the data, we would ideally like at least 2 sampled points per time: one representing the maximum, and other for minimum. The mean data points are generated from the averages, and will be the scatter points which are plotted.
            for t in t_set[band]:
                flux_range = self.data[band]['flux'][self.data[band]['t']==t]
                data_upper.append(N.max(flux_range))
                data_lower.append(N.min(flux_range))
                data_median.append(N.average(flux_range))
                
            #If we want some phase continuity on either side of the main plotted orbit light curve, we ideally want to take parts of the final 3 model orbits. If there are fewer than 3 orbits, tile the single-orbit model on either side. (This will probably result in discontinuities at the borders.)
            if num_orbits < 3:
                y_model[band] = (self.model[band]['model'].reshape(time_resolution*num_orbits))[time_resolution*(num_orbits-1):]
            else:
                y_model[band] = (self.model[band]['model'].reshape(time_resolution*num_orbits))[time_resolution*(num_orbits-3):]

            x[band] = N.array(N.r_[t_set[band] - P/U.d, t_set[band], t_set[band] + P/U.d])

            
            y_median[band] = N.tile(data_median, 3)

            #If we have a partial phase curve, don't draw the data contours as connected between orbits.
            if partial_phase[band]:
                y_low = N.array(data_lower)[~occulted]
                y_high = N.array(data_upper)[~occulted]
                low_smoothed[band] = inter.UnivariateSpline(t_set[band][~occulted], y_low, s=s[band])
                high_smoothed[band] = inter.UnivariateSpline(t_set[band][~occulted], y_high, s=s[band])
                x_smooth[band] = N.linspace(N.min(t_set[band]), N.max(t_set[band]), num=time_resolution)
                lowerlim[band] = N.min([N.min(low_smoothed[band](x_smooth[band])), 1])
                upperlim[band] = N.max([N.max(high_smoothed[band](x_smooth[band])), N.max(self.model[band]['model'])])

            #For full phase curves we will draw connected contours.
            else:
                y_low = N.tile(N.array(data_lower)[~occulted], 3)
                y_high = N.tile(N.array(data_upper)[~occulted], 3)
                low_smoothed[band] = inter.UnivariateSpline(x[band][~N.tile(occulted,3)], y_low, s=s[band])
                high_smoothed[band] = inter.UnivariateSpline(x[band][~N.tile(occulted,3)], y_high, s=s[band])
                x_smooth[band] = N.linspace(-0.5*P/U.d, 0.5*P/U.d, num=time_resolution)
                lowerlim[band] = N.min([N.min(low_smoothed[band](x_smooth[band])), 1])
                upperlim[band] = N.max([N.max(high_smoothed[band](x_smooth[band])), N.max(self.model[band]['model'])])

        combo_lowerlim = N.min(list(lowerlim.values()))
        combo_upperlim = N.max(list(upperlim.values()))

        #Boolean for whether all available phase curves are partial phases.
        all_partial = N.all(list(partial_phase.values()))
        
        for i, band in enumerate(sorted(self.data)):
            band_label = band.replace('p', '.')
            
            axis.plot(plotted_times, y_model[band], linewidth=2, color=color_modbg[band])
            axis.scatter(x[band], y_median[band], color='k', label=r'${0} \mu$m Data'.format(band_label), marker=',', s=8, alpha=0.6, zorder=1)
            axis.fill_between(x_smooth[band], low_smoothed[band](x_smooth[band]), high_smoothed[band](x_smooth[band]), color=color_datlab[band], label=r'${0} \mu$m Data'.format(band_label), alpha=0.4)

            if ~combo & save:

                #If ALL the photometric datasets are partial phase, then we will choose to plot each without phase overlap regions.
                if all_partial:
                    x_min = N.min(N.concatenate(list(t_set.values())))
                    x_max = N.max(N.concatenate(list(t_set.values())))
                    plt.xlim([x_min, x_max])

                #For planets with any full phase photometry, plot shaded regions on either end of the light curve extending beyond the single orbit.
                else:
                    x_min = -(0.5+phase_overlap)*P/U.d
                    x_max = (0.5+phase_overlap)*P/U.d
                    plt.xlim([x_min, x_max])
                    plt.axvspan(x_min, -0.5*P/U.d, facecolor='0.2', alpha=0.4)
                    plt.axvspan(0.5*P/U.d, x_max, facecolor='0.2', alpha=0.4)
                    
                plt.axhline(1, linestyle=(0, (5, 1)))
                
                axis.set_xlabel(r'$t$ (Days)', fontsize=24)
                if i == 0:
                    axis.set_ylabel(r'$F_\lambda/F_{\lambda,\star}$', fontsize=24)
                    plt.setp(axis.get_yticklabels(), fontsize=16)
                else:
                    plt.setp(axis.get_yticklabels(), visible=False)
                plot_lower = combo_lowerlim - 0.05*(combo_upperlim-combo_lowerlim)
                plot_upper = combo_upperlim + 0.05*(combo_upperlim-combo_lowerlim)
                axis.set_ylim([plot_lower, plot_upper])
                
                sig_place = N.floor(N.log10(plot_upper-1))
                y_tickrange = N.linspace(1, (10**sig_place)*N.trunc(10**(-sig_place)*plot_upper), num=3)
                y_ticks = y_tickrange                
                plt.yticks(y_ticks)
                axis.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.4f'))
                
                plt.setp(axis.get_xticklabels(), fontsize=16)
                axis.margins(0.05)
            
                plt.tight_layout
                for filetype in filetypes:
                    plt.savefig('{0}_{1}_{2}.{3}'.format(str(datetime.datetime.today()).replace(" ","_"), self.planet.name.replace(" ",""), band, filetype), bbox_inches='tight', transparent=True)
                plt.cla()

        if combo & save:
            plt.xlim([-(0.5+phase_overlap)*P/U.d, (0.5+phase_overlap)*P/U.d])
            plt.axvspan(-(0.5+phase_overlap)*P/U.d, -0.5*P/U.d, facecolor='0.2', alpha=0.4)
            plt.axvspan(0.5*P/U.d, (0.5+phase_overlap)*P/U.d, facecolor='0.2', alpha=0.4)
            plt.axhline(1, linestyle=(0, (5, 1)))

            plot_lower = combo_lowerlim - 0.05*(combo_upperlim-combo_lowerlim)
            plot_upper = combo_upperlim + 0.05*(combo_upperlim-combo_lowerlim)
            axis.set_ylim([plot_lower, plot_upper])
                
            sig_place = N.floor(N.log10(plot_upper-1))
            y_tickrange = N.linspace(1, (10**sig_place)*N.trunc(10**(-sig_place)*plot_upper), num=3)
            y_ticks = y_tickrange                
            plt.yticks(y_ticks)
            axis.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.4f'))
                
            plt.setp(axis.get_xticklabels(), fontsize=16)
            plt.setp(axis.get_yticklabels(), fontsize=16)
            
            axis.margins(0.05)
            
            plt.tight_layout
            for filetype in filetypes:
                plt.savefig('{0}_{1}_phase.{2}'.format(str(datetime.date.today()).replace(" ", ""), self.planet.name.replace(" ", ""), filetype), bbox_inches='tight', transparent=True)
