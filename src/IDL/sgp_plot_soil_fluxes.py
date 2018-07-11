# -*- coding: utf-8 -*-
'''
Oklahoma (Southern Great Plains) soil flux plots
Reproduce some of the plots in python
'''
import linecache
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy import stats

plt.rcParams.update({'mathtext.default': 'regular'})  # san-serif math

# specify data dir and read data in
data_dir = './input/'
plot_dir = './plots/mdl_plts/'
data_fname = data_dir + 'soil_flux_regular_meas.csv'
data_colnames = linecache.getline(data_fname, 1).replace('\n', '').split(',')

data_ch = np.genfromtxt(data_fname, names=data_colnames,
                        dtype=None, delimiter=',', skip_header=1)

# create aliases
doy_utc = data_ch['doy_utc']
doy_local = data_ch['doy_local']
cos_a = data_ch['cos_a']
co2_a = data_ch['co2_a']
h2o_a = data_ch['h2o_a']
f_cos = data_ch['fcos_corr']
f_co2 = data_ch['fco2']
f_h2o = data_ch['fh2o']
T_air = data_ch['T_air_avg']
T_s = data_ch['soil_temp_avg']
swc = data_ch['swc_avg']

# physical constants
p_atm = 1.01325e5  # Surface pressure (Pa), assumed for the SGP site
p_std = 1.01325e5  # Standard atmosphere pressure (Pa)
R_gas = 8.3144621  # Universal gas constant (J K^-1 mol^-1)
# Air molar concentration (mol m^-3) at the standard condition (273.15 K and 101.325 kPa)
air_conc = p_std / R_gas / 273.15
Mw = 18.016  # molecular weight of water

# data cleaning


'''
Soil COS fluxes and soil variables, a simplified figure
'''

fig1, (f1ax1, f1ax2) = plt.subplots(nrows=2, sharex=True, figsize=(8, 8))
f1ax1.plot([90, 160], [0, 0], '--', color='gray')
f1ax1.plot(doy_local, f_cos, '.-', color='gray',
           markeredgecolor='black', markerfacecolor='black')
f1ax1.annotate(
    'harvest', xy=(144, 15.),
    xytext=(144, 20),
    color='crimson', fontsize=14,
    arrowprops=dict(edgecolor='crimson', arrowstyle='->'),
    horizontalalignment='center', verticalalignment='bottom',)
f1ax1.annotate(
    'growing season', xy=(105., 10),
    color='forestgreen', fontsize=14, horizontalalignment='center',
    verticalalignment='bottom',)
f1ax1.annotate('senescence', xy=(138., 12), color='peru', fontsize=14,
               horizontalalignment='center', verticalalignment='bottom',)
Ts_ln, = f1ax2.plot(doy_local, T_s, '-', color='#cc4c02', label='temp')
f1ax22 = f1ax2.twinx()
swc_ln, = f1ax22.plot(doy_local, swc, '-', color='#41b6c4', label='moisture')
f1leg2 = f1ax2.legend(
    handles=[Ts_ln, swc_ln],
    loc='upper left', labelspacing=0., fontsize=12)


def color_legend_texts(leg):
    '''Color legend texts based on color of corresponding lines'''
    for line, txt in zip(leg.get_lines(), leg.get_texts()):
        txt.set_color(line.get_color())


color_legend_texts(f1leg2)

f1ax2.set_xlim((90, 154))
f1ax2.xaxis.set_minor_locator(AutoMinorLocator())
f1ax2.set_xlabel('Date or day of year (2012)', fontsize=14)
f1ax1.set_ylabel('COS flux (pmol m$^{-2}$ s$^{-1}$)', fontsize=14)
f1ax2.set_ylabel('Soil temperature ($\degree$C)', fontsize=14)
f1ax22.set_ylabel('Soil moisture (m$^3$ m$^{-3}$)', fontsize=14)
f1ax1.set_title('Southern Great Plains, OK, USA\nwheat field', color='navy')
# add description for date
xax_doy_array = f1ax2.xaxis.get_ticklocs()
xax_doy_array = xax_doy_array.astype('int')
xax_date_array = [''] * len(xax_doy_array)
year = 2012
for i in range(len(xax_doy_array)):
    xax_date_array[i] = datetime.date.fromordinal(datetime.date(
        year, 1, 1).toordinal() + xax_doy_array[i]).strftime('%b %d')  # '%m/%d %H:%M'

xax_doy_array = map(str, xax_doy_array)
xax_combined_array = map(''.join, zip(
    xax_date_array, np.repeat('\n', len(xax_doy_array)), xax_doy_array))
# print(xax_date_array)
# print(xax_doy_array)
# print(xax_combined_array)
f1ax2.xaxis.set_ticklabels(xax_combined_array)
fig1.tight_layout()
fig1.savefig(plot_dir + '/fcos_and_soil_var.eps')
fig1.savefig(plot_dir + '/png/fcos_and_soil_var.png')
