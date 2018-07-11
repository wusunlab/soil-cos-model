import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import rc

C_atm = 500.  # COS conc 500 pptv
air_conc = 1.01325e5 / 8.3145 / 298.15
T_soil = 298.15
k = 2e-3  # dep velocity, ~ 1 mm s^-1

def func_D_g(theta_w, theta_sat, T):
    D_g = 1.337e-5 * (theta_sat - theta_w)**2 * (1.-theta_w/theta_sat)**(3/5.3) * (298.15/T)**1.5
    return D_g

def func_D_aq(theta_w, theta_sat, T):
    D_aq = 10**(-1010./T - 1.3246) * 1e-4 * theta_w**2 * (theta_w/theta_sat)**(5.3/3-1.)
    return D_aq

theta_sat1 = 0.50
wfps1 = np.linspace(0,1,201)
theta_w1 = theta_sat1 * wfps1

theta_sat2 = 0.80
wfps2 = np.linspace(0,1,201)
theta_w2 = theta_sat2 * wfps2

D_g1 = func_D_g(theta_w1, theta_sat1, T_soil)
D_aq1 = func_D_aq(theta_w1, theta_sat1, T_soil)

D_g2 = func_D_g(theta_w2, theta_sat2, T_soil)
D_aq2 = func_D_aq(theta_w2, theta_sat2, T_soil)

F_tot1 = np.sqrt(k * D_g1) * C_atm * air_conc + np.sqrt(k * D_aq1) * C_atm * air_conc  # in pmol m^-2 s^-1
F_ratio1 = np.sqrt(D_aq1/D_g1)

F_tot2 = np.sqrt(k * D_g2) * C_atm * air_conc + np.sqrt(k * D_aq2) * C_atm * air_conc  # in pmol m^-2 s^-1
F_ratio2 = np.sqrt(D_aq2/D_g2)

# rc('text', usetex=True)

plt.subplot(2, 2, 1)
plt.plot(wfps1, F_ratio1, 'k-')
plt.title('Total porosity 0.50')
plt.ylabel(r'$\mathsf{F_{aq}/F_{g}}$', fontsize=16)
plt.grid(True)
plt.yscale('log')
plt.subplot(2, 2, 3)
plt.plot(wfps1, F_tot1, 'k-')
plt.ylabel(r'Total flux (pmol m$^\mathsf{-2}$ s$^\mathsf{-1}$)')
plt.grid(True)
plt.xlabel('Water-filled pore space fraction')

plt.subplot(2, 2, 2)
plt.plot(wfps2, F_ratio2, 'k-')
plt.title('Total porosity 0.80')
plt.grid(True)
#plt.ylabel('F_aq/F_g')
plt.yscale('log')
plt.subplot(2, 2, 4)
plt.plot(wfps2, F_tot2, 'k-')
plt.grid(True)
#plt.ylabel('Total flux')
plt.xlabel('Water-filled pore space fraction')

# a single figure
plt.figure(figsize=(5,4))
plt.plot(wfps1, F_ratio1, 'k-')
plt.ylabel(r'$\mathsf{F_{aq}/F_{g}}$', fontsize=16)
plt.xlabel('Water-filled pore space fraction')
plt.grid(True)
plt.yscale('log')
plt.tight_layout()
plt.show()