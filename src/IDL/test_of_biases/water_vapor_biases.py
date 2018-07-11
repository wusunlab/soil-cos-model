# test how water vapor biases soil chamber fluxes
import numpy as np
import scipy.constants.constants as const
import matplotlib.pyplot as plt
import scipy.stats as stats

# environmental settings
V_ch = 4076.1 * 1e-6    # in m3
A_ch = 317.8 * 1e-4    # in m2
T_atm = 293.15    # air temperature 20 C
p_atm = const.atm
air_conc = p_atm / const.R / T_atm
flow_rate = 1./60 * air_conc * 1e-3    # flow rate, 1 L/min
C_w_atm = 0.    # no H2O, dry air
C_s_atm = 500. * 1e-12 * air_conc    # COS at 500 pptv
C_c_atm = 400. * 1e-6 * air_conc    # CO2 at 400 ppmv
F_w = (np.arange(4)+1)*5./100 * 1e-3    # H2O fluxes, 0.05, 0.10, 0.15, 0.20 mmol m-2 s-1
F_s = np.linspace(-8,2,101) * 1e-12    # COS fluxes, -8 to +2 pmol m-2 s-1
F_c = np.linspace(0,5,101) * 1e-6    # CO2 fluxes, 0 to 5 umol m-2 s-1
dt = 1.    # time step, 1 sec
duration = 480.    # chamber-close period, 8 min
N_t = int(duration/dt)

i = 3
j = 70
C_w = np.zeros(N_t+1)
C_s = np.zeros(N_t+1)
C_c = np.zeros(N_t+1)
C_w[0] = C_w_atm
C_s[0] = C_s_atm
C_c[0] = C_c_atm
# analytical solution of water concentration
C_w[1:] = C_w_atm + F_w[i] * A_ch / flow_rate * (1. 
	- np.exp(- flow_rate * (np.arange(N_t)+1) / V_ch))
for k in np.arange(N_t)+1:
	C_s[k] = C_s[k-1] + dt * flow_rate / V_ch * (C_s_atm - C_s[k-1] + F_s[j] * A_ch / flow_rate)

chi_s = C_s * air_conc / (air_conc - C_w)
plt.plot(np.arange(N_t+1), C_s/air_conc*1e12, 'b-')
plt.plot(np.arange(N_t+1), chi_s/air_conc*1e12, 'r-')
plt.show()
y_fit = chi_s - C_s_atm
x_fit = 1 - np.exp(-flow_rate * np.arange(N_t+1) / V_ch)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit,y_fit)
F_s_calc = slope * flow_rate / A_ch