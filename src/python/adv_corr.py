import copy
import numpy as np
import scipy as sp
from soil_tracers_config import *
from create_grid import *
import scipy.sparse as spm  # scipy sparse matrix
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pres_solver import pres_solver_cyl2D
from flux_solver import *

cyl_grid = EquilGridCyl2D()  # create a uniform cylindrical 2D grid
r_ch = 0.1  # radius of the soil chamber
p_ch = -5.  # pressure difference inside the chamber compared with the ambient
#P_0 = 1.01325e5
pres_field = np.zeros((cyl_grid.N_z+1, cyl_grid.N_r+1))
n_within = np.argmin(np.abs(cyl_grid.r - r_ch))
# a guessed initialization with the top boundary condition satisfied
pres_field[:, 0:n_within+1] = p_ch * np.exp(-cyl_grid.z_grid[:, 0:n_within+1])
pres_field[:, n_within+1:cyl_grid.N_r+1] = p_ch * np.sinh(
    cyl_grid.z_grid[:, n_within+1:cyl_grid.N_r+1]/3.)
pres_field = pres_field + p_ch * 0.05 * np.random.random(
    [pres_field.shape[0], pres_field.shape[1]])

# # initial condition
# pres_field[0,0:n_within+1] = p_ch
# pres_field = pres_field + p_ch * 0.05 * np.random.random([pres_field.shape[0], pres_field.shape[1]])

total_poros = np.zeros(cyl_grid.N_z+1) + 0.75   # 0.45
water_cont = 0.4    # cyl_grid.z * 0.05 + 0.15
air_poros = total_poros - water_cont
soil_temp = 25.     # 25 - 5 * cyl_grid.z
dynam_visc = 1.83e-5
# dynamic viscosity (dry air at 25 C, From Warrick eds Soil physics companion, pp. 299)
int_permeab = 4e-10  # from Lund et al. (1999) GCB, we need our own measurements
min_sat = 0.11  # from Lund et al. (1999) GCB
eff_sat = (water_cont/total_poros - min_sat) / (1 - min_sat)
rel_permeab = (1. - eff_sat)**3
permeab = int_permeab * rel_permeab
#pres_diffusivity = permeab * P_0 / air_poros / dynam_visc
#pres_diffusivity = np.zeros((cyl_grid.n_vert_level, cyl_grid.n_rad_level)) + pres_diffusivity

pres_field_sol, vr_field, vz_field, v_field = pres_solver_cyl2D(
    pres_field, permeab, dynam_visc, soil_temp, p_ch, r_ch, cyl_grid)

# visualize
plt.contourf(cyl_grid.r_grid, cyl_grid.z_grid, pres_field_sol,
             100, cmap=cm.Spectral)  # pressure field
plt.axis([0., 1., 1., 0])
plt.xlabel(r'$r$ (m)', fontsize=14)
plt.ylabel(r'$z$ (m)', fontsize=14)
plt.colorbar(ticks=[-np.arange(11)/10.], label='Pressure (Pa)')
plt.streamplot(cyl_grid.r_grid, cyl_grid.z_grid, vr_field,
               vz_field, color=(1., 1., 1., 0.67))  # streamline
plt.show()

plt.contourf(cyl_grid.r_grid, cyl_grid.z_grid,
             v_field * 1e5, 100, cmap=cm.coolwarm)
plt.axis([0., 1., 1., 0])
plt.xlabel(r'$r$ (m)', fontsize=14)
plt.ylabel(r'$z$ (m)', fontsize=14)
plt.colorbar(
    ticks=[np.arange(13) * 0.5],
    label=r'Advection velocity ($\times$ 10$^{-5}$ m/s)')
plt.quiver(cyl_grid.r_grid, cyl_grid.z_grid,
           vr_field, -vz_field, color=(1., 1., 1., 0.67))
plt.show()


# # interpolate the pressure and velocity fields to the FV cylindrical grid
# cyl_FVGrid = FVGridCyl2D()
# func_interp_pres = sp.interpolate.interp2d(cyl_grid.r_grid, cyl_grid.z_grid, pres_field_sol, kind='linear')
# pres_field_FV = func_interp_pres(cyl_FVGrid.rad_node, cyl_FVGrid.vert_node)
# # rectify the top boundary condition
# pres_field_FV[0,cyl_FVGrid.r_grid[0,:] < r_ch] = p_ch
# pres_field_FV[0,cyl_FVGrid.r_grid[0,:] > r_ch] = 0.
# # evaluate the advection velocities
# func_interp_vr = sp.interpolate.interp2d(cyl_grid.r_grid, cyl_grid.z_grid, vr_field, kind='linear')
# func_interp_vz = sp.interpolate.interp2d(cyl_grid.r_grid, cyl_grid.z_grid, vz_field, kind='linear')
# vr_field_FV = func_interp_vr(cyl_FVGrid.rad_node, cyl_FVGrid.vert_node)
# vz_field_FV = func_interp_vz(cyl_FVGrid.rad_node, cyl_FVGrid.vert_node)
# v_field_FV = np.sqrt(vr_field_FV ** 2 + vz_field_FV **2)


# extract the values for the compressed FV grid with known indices
cyl_FVGrid = FVGridCyl2D()
pres_field_FV = cyl_FVGrid.r_grid * 0.
vr_field_FV = pres_field_FV * 0.
vz_field_FV = pres_field_FV * 0.
for j in range(cyl_FVGrid.n_vert_level):
    pres_field_FV[j, :] = pres_field_sol[cyl_FVGrid.vert_index[j],
                                         cyl_FVGrid.rad_index]
    vr_field_FV[j, :] = vr_field[cyl_FVGrid.vert_index[j], cyl_FVGrid.rad_index]
    vz_field_FV[j, :] = vz_field[cyl_FVGrid.vert_index[j], cyl_FVGrid.rad_index]

v_field_FV = np.sqrt(vr_field_FV ** 2 + vz_field_FV ** 2)


# vr_field_FV = pres_field_FV * 0.
# vz_field_FV = pres_field_FV * 0.

# for i in range(cyl_FVGrid.n_rad_level-1):
#     if i == 0:
#         dr1 = cyl_FVGrid.rad_node[1] - cyl_FVGrid.rad_node[0]
#         dr2 = cyl_FVGrid.rad_node[2] - cyl_FVGrid.rad_node[1]
#         vr_field_FV[:,i] = (- ((dr1+dr2)**2 - dr1**2) * pres_field_FV[:,0] +
#             (dr1+dr2)**2 * pres_field_FV[:,1] - dr1**2 * pres_field_FV[:,2]) / ( dr1*dr2*(dr1+dr2) )
#     else:
#         dr1 = cyl_FVGrid.rad_node[i] - cyl_FVGrid.rad_node[i-1]
#         dr2 = cyl_FVGrid.rad_node[i+1] - cyl_FVGrid.rad_node[i]
#         vr_field_FV[:,i] = (- dr2**2 * pres_field_FV[:,i-1] +
#             (dr2**2 - dr1**2) * pres_field_FV[:,i] + dr1**2 * pres_field_FV[:,i+1]) / ( dr1*dr2*(dr1+dr2) )

# vr_field_FV[:,-1] = 0.  # zero flux at the side boundaries
# permeab_FV = np.interp(cyl_FVGrid.vert_node, cyl_grid.z, permeab)
# vr_field_FV = - np.tile(permeab_FV, (cyl_FVGrid.n_rad_level,1)).transpose() / dynam_visc * vr_field_FV

# for i in range(cyl_FVGrid.n_vert_level-1):
#     if i == 0:
#         dz1 = cyl_FVGrid.vert_node[1] - cyl_FVGrid.vert_node[0]
#         dz2 = cyl_FVGrid.vert_node[2] - cyl_FVGrid.vert_node[1]
#         vz_field_FV[i,:] = (- ((dz1+dz2)**2 - dz1**2) * pres_field_FV[0,:] +
#             (dz1+dz2)**2 * pres_field_FV[1,:] - dz1**2 * pres_field_FV[2,:]) / ( dz1*dz2*(dz1+dz2) )
#     else:
#         dz1 = cyl_FVGrid.vert_node[i] - cyl_FVGrid.vert_node[i-1]
#         dz2 = cyl_FVGrid.vert_node[i+1] - cyl_FVGrid.vert_node[i]
#         vz_field_FV[i,:] = (- dz2**2 * pres_field_FV[i-1,:] +
#             (dz2**2 - dz1**2) * pres_field_FV[i,:] + dz1**2 * pres_field_FV[i+1,:]) / ( dz1*dz2*(dz1+dz2) )

# vz_field_FV[-1,:] = 0.  # zero flux at the bottom boundaries
# vz_field_FV = - np.tile(permeab_FV, (cyl_FVGrid.n_rad_level,1)).transpose() / dynam_visc * vz_field_FV
# v_field_FV = np.sqrt(vr_field_FV ** 2 + vz_field_FV **2)


# visualize
plt.contourf(cyl_FVGrid.r_grid, cyl_FVGrid.z_grid, pres_field_FV,
             100, cmap=cm.Spectral)  # pressure field
plt.axis([0., 1., 1., 0])
plt.xlabel(r'$r$ (m)', fontsize=14)
plt.ylabel(r'$z$ (m)', fontsize=14)
plt.colorbar(ticks=[-np.arange(11)/10.], label='Pressure (Pa)')
plt.streamplot(cyl_FVGrid.r_grid, cyl_FVGrid.z_grid, vr_field_FV,
               vz_field_FV, color=(1., 1., 1., 0.67))  # streamline
plt.show()

plt.contourf(cyl_FVGrid.r_grid, cyl_FVGrid.z_grid,
             v_field_FV * 1e5, 100, cmap=cm.coolwarm)
plt.axis([0., 1., 1., 0])
plt.xlabel(r'$r$ (m)', fontsize=14)
plt.ylabel(r'$z$ (m)', fontsize=14)
plt.colorbar(
    ticks=[np.arange(13) * 0.5],
    label=r'Advection velocity ($\times$ 10$^{-5}$ m/s)')
plt.quiver(cyl_FVGrid.r_grid, cyl_FVGrid.z_grid,
           vr_field_FV, -vz_field_FV, color=(1., 1., 1., 0.67))
plt.show()

# initialize the concentration profile, assuming steady-state
atm_conc = air_conc_std * 4e-4 * 273.15 / (273.15+25.)
conc_field = atm_conc * np.ones((cyl_FVGrid.n_vert_level,
                                 cyl_FVGrid.n_rad_level))
temp_field = 25. * np.ones((cyl_FVGrid.n_vert_level, cyl_FVGrid.n_rad_level))
poros_field = 0.75 * np.ones((cyl_FVGrid.n_vert_level, cyl_FVGrid.n_rad_level))
swc_field = 0.4 * np.ones((cyl_FVGrid.n_vert_level, cyl_FVGrid.n_rad_level))


# co2_ch_conc = ...
# flux_solver_cyl2D('co2', temp_field, poros_field, swc_field, atm_conc=air_conc_std * 4e-4, FVGrid=cyl_FVGrid,
#     steady_state=False, duration=900, chamber=True, chamber_vol=4e-4, chamber_area=3e-2, flow_rate=0.015, save_headspace_conc=True,
#     gas_only=True, vr_field=vr_field_FV, vz_field=vz_field_FV)

# test co2 advection
co2_prod_profile = soil_co2_source(
    temp_field[:, 0],
    swc_field[:, 0],
    cyl_FVGrid.vert_node, vmax=2e-5, decay_depth=0.2)
np.sum(co2_prod_profile[1:] * cyl_FVGrid.vert_cv_size[1:])


# # change the source-sink parameters
# test_source_params = copy.copy(source_params)
# test_source_params['soil_resp_decay_depth'] = 0.05
# test_source_params['soil_co2_source_vmax'] = 4e-5


ch_conc = flux_solver_cyl2D(
    'co2', temp_field, poros_field, swc_field, atm_conc=atm_conc,
    source_profile=co2_prod_profile, adv_vel_U=vr_field_FV,
    adv_vel_V=vz_field_FV, dt=1., FV_grid=cyl_FVGrid, duration=480.,
    flow_rate=1e-3 / 60.)

ch_conc_no_adv = flux_solver_cyl2D(
    'co2', temp_field, poros_field, swc_field, atm_conc=atm_conc,
    source_profile=co2_prod_profile * 5, adv_vel_U=vr_field_FV * 0.,
    adv_vel_V=vz_field_FV * 0., dt=1., FV_grid=cyl_FVGrid, duration=480.,
    flow_rate=1e-3 / 60.)

line1, = plt.plot(ch_conc * 298.15/273.15 / air_conc_std * 1e6, linewidth=2)
plt.xlabel('time since chamber closure (sec)')
plt.ylabel('chamber conc (ppmv)')
line2, = plt.plot(
    ch_conc_no_adv * 298.15/273.15 / air_conc_std * 1e6, linewidth=2)
plt.legend(
    [line1, line2],
    ['1x CO2 production,\n5 Pa under-pressurization',
     '5x CO2 production,\nno advection'],
    loc='lower right')
plt.show()

# test cos advection
