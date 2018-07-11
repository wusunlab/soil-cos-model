'''
Extrapolate temperature profile from surface temperature measurements
'''
import linecache
import copy
import numpy as np
from soil_tracers_config import *
import matplotlib.pyplot as plt
met_data_file = './data/soil_temp_and_radiation_2013.csv'

met_data_colnames = linecache.getline(
    met_data_file, 1).replace(
    '\n', '').split(',')
met_data = np.genfromtxt(met_data_file, delimiter=',',
                         skip_header=1, names=met_data_colnames)

doy_local = met_data[:]['doy_local']
temp = met_data[:]['soil_temp_5cm']
temp_10cm = met_data[:]['soil_temp_10cm']
temp_20cm = met_data[:]['soil_temp_20cm']
temp_50cm = met_data[:]['soil_temp_50cm']

# linear interpolation to get even-sampled series for Fourier transform
sampling_interval = 1./144.  # every 10 minutes
doy_local_interp = np.linspace(
    doy_local[0],
    doy_local[-1],
    num=(doy_local[-1] - doy_local[0]) / sampling_interval + 1.)
temp_interp = np.interp(doy_local_interp, doy_local, temp)
# plot to see whether the interpolation makes sense
plt.plot(doy_local, temp)
plt.plot(doy_local_interp, temp_interp, 'r--')
plt.show()
plt.close()

# Fourier transform
N_fft = temp_interp.size
# in day^-1
freq_fft = 1. / (sampling_interval * N_fft) * np.arange(N_fft/2+1)
freq_fft = np.append(freq_fft, 1. / (sampling_interval * N_fft)
                     * np.arange(-N_fft/2+1, 0, 1))
A_k = np.fft.fft(temp_interp)   # amplitudes
plt.plot(freq_fft, np.abs(A_k)**2)
plt.yscale('log')
plt.xlabel('Frequency (d$^{-1}$)')

# assign a reasonable value for the thermal diffusivty at the Stunt Ranch, m2 s-1
therm_diff = 7.5e-7
damping_depth = np.sqrt(therm_diff * 2 / np.abs(freq_fft * 2 * np.pi) * 86400.)

# simulate the temp series at 20 cm
delta_z = 0.5
A_k_prime = copy.copy(A_k)
A_k_prime[1:] = A_k[1:] * np.exp(-delta_z / damping_depth[1:]) * \
    np.exp(np.vectorize(complex)(0, -delta_z / damping_depth[1:]))
temp_prime = np.fft.ifft(A_k_prime)
plt.plot(doy_local_interp, temp_prime, 'r-')
plt.plot(doy_local, temp_50cm, 'b--')
plt.plot(doy_local_interp, temp_interp, 'k--')
