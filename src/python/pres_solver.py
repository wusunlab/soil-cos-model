import numpy as np
from soil_tracers_config import *
from create_grid import *
import scipy.sparse as spm  # scipy sparse matrix
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def pres_solver_cyl2D(pres_field, permeab, dynam_visc, soil_temp, p_ch, r_ch, cyl_grid):
    N_r = cyl_grid.N_r
    N_z = cyl_grid.N_z
    dr = cyl_grid.dr
    dz = cyl_grid.dz
    n_within = np.argmin(np.abs(cyl_grid.r - r_ch))
    pres_vec_0 = pres_field.reshape(pres_field.shape[0] * pres_field.shape[1])
    pres_vec = pres_vec_0
    kappa = permeab / dynam_visc / (soil_temp + 273.15)
    # normalize the kappa coefficients because the values are small
    kappa = kappa / np.median(kappa) * cyl_grid.dr * cyl_grid.dz
    # this accelerates the convergence but does not influence the result
    coef_mat = np.zeros((pres_vec.shape[0], pres_vec.shape[0]))
    # initialize the coefficient matrix
    coef_mat[np.arange(N_r+1),np.arange(N_z+1)] = 1.
    # for grid points at the top boundary
    coef_mat[(N_r+1)*N_z + np.arange(N_r+1), (N_r+1)*(N_z-1) + np.arange(N_r+1)] = 1.
    # for grid points at the bottom boundary
    coef_mat[(N_r+1) * np.arange(1,N_z) + N_r, (N_r+1) * np.arange(1,N_z) + N_r-1 ] = 1.
    # for grid points at the side boundaries

    for j in range(1,N_z):
        # for grid points at the axis
        # coef_mat[(N_r+1)*j, (N_r+1)*j+1] = 10 * kappa[j] / dr**2 / (4 * kappa[j] / dr**2 - 2 * kappa[j] / dz**2) 
        # coef_mat[(N_r+1)*j, (N_r+1)*j+2] = -8 * kappa[j] / dr**2 / (4 * kappa[j] / dr**2 - 2 * kappa[j] / dz**2)
        # coef_mat[(N_r+1)*j, (N_r+1)*j+3] = 2 * kappa[j] / dr**2 / (4 * kappa[j] / dr**2 - 2 * kappa[j] / dz**2)
        # coef_mat[(N_r+1)*j, (N_r+1)*(j+1)] = - (kappa[j] + kappa[j+1])/2. / dz**2 / (
        #   4 * kappa[j] / dr**2 - 2 * kappa[j] / dz**2)
        # coef_mat[(N_r+1)*j, (N_r+1)*(j-1)] = - (kappa[j-1] + kappa[j])/2. / dz**2 / (
        #   4 * kappa[j] / dr**2 - 2 * kappa[j] / dz**2)
        coef_mat[(N_r+1)*j, (N_r+1)*j+1] = 4 * kappa[j] / dr**2 / (4 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2) 
        coef_mat[(N_r+1)*j, (N_r+1)*(j+1)] = (kappa[j] + kappa[j+1])/2. / dz**2 / (
            4 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2)
        coef_mat[(N_r+1)*j, (N_r+1)*(j-1)] = (kappa[j-1] + kappa[j])/2. / dz**2 / (
            4 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2)
    for j in range(1,N_z):
        for i in range(1,N_r):
            coef_mat[(N_r+1)*j+i, (N_r+1)*j+i+1] = (kappa[j] / dr**2 + kappa[j]/2./cyl_grid.r[i]/dr) / (
                2 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2)
            coef_mat[(N_r+1)*j+i, (N_r+1)*j+i-1] = (kappa[j] / dr**2 - kappa[j]/2./cyl_grid.r[i]/dr) / (
                2 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2)
            coef_mat[(N_r+1)*j+i, (N_r+1)*(j+1)+i] = (kappa[j+1] + kappa[j])/2./dz**2 / (
                2 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2)
            coef_mat[(N_r+1)*j+i, (N_r+1)*(j-1)+i] = (kappa[j] + kappa[j-1])/2./dz**2 / (
                2 * kappa[j] / dr**2 + 2 * kappa[j] / dz**2)
    coef_mat_rowsum = np.sum(coef_mat, axis=1)  # test
    coef_mat_spm = spm.coo_matrix(coef_mat)
    idn_mat = spm.coo_matrix(np.eye(pres_vec.shape[0]))
    #
    # use simple BiCGSTAB algorithm to solve the linear system
    # (I-T)P = 0
    # failed attempt
    # resd_0 = -(idn_mat - coef_mat_spm) * pres_vec_0
    # resd_old = resd_0
    # rho_old = 1.
    # alpha = 1.
    # omega_old = 1.
    # p_old = pres_vec * 0.
    # v_old = pres_vec * 0.
    # for loop_counter in range(10000):
    #   rho_new = np.dot(resd_0, resd_old)
    #   beta = (rho_new/rho_old) * (alpha/omega_old)
    #   p_new = resd_old + beta * (p_old - omega_old * v_old)
    #   v_new = (idn_mat-coef_mat_spm) * p_new
    #   alpha = rho_new / np.dot(resd_0, v_new)
    #   s = resd_old - alpha * v_new
    #   t = (idn_mat-coef_mat_spm) * s
    #   omega_new = np.dot(t,s)/np.dot(t,t)
    #   pres_vec_new = pres_vec + alpha * p_new + omega_new * s
    #   conv_fac = np.linalg.norm(pres_vec_new - pres_vec) / np.linalg.norm(pres_vec)
    #   if conv_fac < 1e-6:
    #       break
    #   # update the solution
    #   pres_vec = pres_vec_new
    #   # apply boundary conditions again
    #   pres_vec[0:n_within+1] = p_ch
    #   pres_vec[n_within+1:N_r+1] = 0.
    #   pres_vec[(N_r+1)*N_z + np.arange(N_r+1)] = pres_vec[(N_r+1)*(N_z-1) + np.arange(N_r+1)] 
    #   pres_vec[(N_r+1)*np.arange(N_z+1)+N_r] = pres_vec[(N_r+1)*np.arange(N_z+1)+N_r-1]  
    #   # get the new residual vector
    #   # resd_new = s - omega_new * t
    #   resd_new = -(idn_mat-coef_mat_spm) * pres_vec
    #   # save
    #   rho_old = rho_new
    #   omega_old = omega_new
    #   p_old = p_new
    #   v_old = v_new
    #   resd_old = resd_new
    
    omega = 1.  # over-relaxation factor
    for loop_counter in range(100000):
        # pres_vec_new = np.dot(coef_mat, pres_vec)
        pres_vec_new = coef_mat_spm * pres_vec
        # pres_vec_new = omega * pres_vec_new + (1.-omega) * pres_vec
        pres_vec_new[0:n_within+1] = p_ch
        pres_vec_new[n_within+1:N_r+1] = 0.
        pres_vec_new[(N_r+1)*N_z + np.arange(N_r+1)] = pres_vec_new[(N_r+1)*(N_z-1) + np.arange(N_r+1)] 
        pres_vec_new[(N_r+1)*np.arange(N_z+1)+N_r] = pres_vec_new[(N_r+1)*np.arange(N_z+1)+N_r-1]  
        # if converges, break
        conv_fac = np.linalg.norm(pres_vec_new - pres_vec) / np.linalg.norm(pres_vec)
        # conv_fac2 = np.nanmax(np.abs((pres_vec_new - pres_vec)/pres_vec))
        # conv_fac = np.linalg.norm(-(idn_mat-coef_mat_spm) * pres_vec)
        if conv_fac < 1e-6:
            break
        pres_vec = pres_vec_new
    pres_field_sol = pres_vec.reshape((N_z+1,N_r+1))
    # evaluate the gradients
    vr_field = pres_field_sol * 0.
    vr_field[:,0] = (-1.5 * pres_field_sol[:,0] + 2. * pres_field_sol[:,1] - 0.5 * pres_field_sol[:,2]) / dr
    # vr_field[:,N_r] = (1.5 * pres_field_sol[:,N_r] - 2. * pres_field_sol[:,N_r-1] + 0.5 * pres_field_sol[:,N_r-2]) / dr
    # at the outer boundary, v_r = 0
    vr_field[:,N_r] = 0.
    vr_field[:,np.arange(1,N_r)] = (pres_field_sol[:,np.arange(2,N_r+1)] - pres_field_sol[:,np.arange(0,N_r-1)]) / 2. / dr
    vr_field[0,5:7] = 0.  # at chamber boundary, v_r = 0
    vz_field = pres_field_sol * 0.
    vz_field[0,:] = (-1.5 * pres_field_sol[0,:] + 2. * pres_field_sol[1,:] - 0.5 * pres_field_sol[2,:]) / dz
    # vz_field[N_z,:] = (1.5 * pres_field_sol[N_z,:] - 2. * pres_field_sol[N_z-1,:] + 0.5 * pres_field_sol[N_z-2,:]) / dz
    # at the bottom boundary, v_z = 0
    vz_field[N_z,:] = 0.
    vz_field[np.arange(1,N_z),:] = (pres_field_sol[np.arange(2,N_z+1),:] - pres_field_sol[np.arange(0,N_z-1),:]) / 2. / dz
    # multiply them by the coefficients
    vr_field = - np.tile(permeab, (N_r+1,1)).transpose() / dynam_visc * vr_field
    vz_field = - np.tile(permeab, (N_r+1,1)).transpose() / dynam_visc * vz_field
    v_field = np.sqrt(vr_field ** 2 + vz_field **2)

    return pres_field_sol, vr_field, vz_field, v_field

    # # test
    # plt.contourf(cyl_grid.r_grid, cyl_grid.z_grid, pres_field_sol)
    # plt.axis([0., 1., 1., 0])
    # plt.show()
    # plt.contourf(cyl_grid.r_grid, cyl_grid.z_grid, np.log10(v_field))
    # plt.axis([0., 1., 1., 0])
    # plt.show()
    # plt.streamplot(cyl_grid.r_grid, cyl_grid.z_grid, vr_field, vz_field)
    # plt.axis([0., 1., 1., 0])
    # #plt.colorbar()
    # plt.show()


'''
# an example
cyl_grid = EquilGridCyl2D() # create a uniform cylindrical 2D grid
r_ch = 0.1  # radius of the soil chamber
p_ch = -1.  # pressure difference inside the chamber compared with the ambient
#P_0 = 1.01325e5
pres_field = np.zeros((cyl_grid.N_z+1, cyl_grid.N_r+1))
n_within = np.argmin(np.abs(cyl_grid.r - r_ch))
# a guessed initialization with the top boundary condition satisfied
pres_field[:,0:n_within+1] = p_ch * np.exp(-cyl_grid.z_grid[:,0:n_within+1])
pres_field[:,n_within+1:cyl_grid.N_r+1] = p_ch * np.sinh(cyl_grid.z_grid[:,n_within+1:cyl_grid.N_r+1]/3.)
pres_field = pres_field + p_ch * 0.05 * np.random.random([pres_field.shape[0], pres_field.shape[1]])

# # initial condition
# pres_field[0,0:n_within+1] = p_ch
# pres_field = pres_field + p_ch * 0.05 * np.random.random([pres_field.shape[0], pres_field.shape[1]])

total_poros = np.zeros(cyl_grid.N_z+1) + 0.45
water_cont = cyl_grid.z * 0.05 + 0.15
air_poros = total_poros - water_cont
soil_temp = 25 - 5 * cyl_grid.z
dynam_visc = 1.83e-5
# dynamic viscosity (dry air at 25 C, From Warrick eds Soil physics companion, pp. 299)
int_permeab = 1e-10 # from Lund et al. (1999) GCB, we need our own measurements
min_sat = 0.11 # from Lund et al. (1999) GCB
eff_sat = (water_cont/total_poros - min_sat) / (1 - min_sat)
rel_permeab = (1. - eff_sat)**3
permeab = int_permeab * rel_permeab 
#pres_diffusivity = permeab * P_0 / air_poros / dynam_visc
#pres_diffusivity = np.zeros((cyl_grid.n_vert_level, cyl_grid.n_rad_level)) + pres_diffusivity

pres_field_sol, vr_field, vz_field, v_field = pres_solver_cyl2D(pres_field, permeab, dynam_visc, soil_temp, p_ch, r_ch, cyl_grid)
# visualize
plt.contourf(cyl_grid.r_grid, cyl_grid.z_grid, pres_field_sol, 100, cmap=cm.Spectral) # pressure field
plt.axis([0., 1., 1., 0])
plt.xlabel(r'$r$ (m)',fontsize=14)
plt.ylabel(r'$z$ (m)',fontsize=14)
plt.colorbar(ticks=[-np.arange(11)/10.], label='Pressure (Pa)')
plt.streamplot(cyl_grid.r_grid, cyl_grid.z_grid, vr_field, vz_field, color=(1.,1.,1.,0.67)) # streamline
plt.show()

plt.contourf(cyl_grid.r_grid, cyl_grid.z_grid, v_field * 1e5, 100, cmap=cm.coolwarm)
plt.axis([0., 1., 1., 0])
plt.xlabel(r'$r$ (m)',fontsize=14)
plt.ylabel(r'$z$ (m)',fontsize=14)
plt.colorbar(ticks=[np.arange(13) * 0.5], label=r'Advection velocity ($\times$ 10$^{-5}$ m/s)')
plt.quiver(cyl_grid.r_grid, cyl_grid.z_grid, vr_field, -vz_field, color=(1.,1.,1.,0.67))
plt.show()
'''

