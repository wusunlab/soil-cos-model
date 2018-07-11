# SoilTracers main program
# Wu Sun (2014) UCLA

import numpy as np
import scipy
import matplotlib.pyplot as plt
import create_grid


def eg_co2_profile():
    # an example
    total_poros = 0.45
    air_poros = 0.25
    T_soil = 15
    co2_atm = 400
    co2_atm = co2_atm * 1e-6 * P_STD / (T_soil + 273.15) / R_gas
    dt = 1

    z_grid = create_grid.FVGrid1D()
    decay_depth = 0.4
    source_profile = 5e-6 * np.exp(-z_grid.grid_node/decay_depth)

    co2_flux, co2_profile = flux_solver(z_grid.grid_node, z_grid.cv_size, dt,
        T_soil, total_poros, water_cont, co2_atm, source_profile, f_biomass, lit_wc, litter_layers,
        vmax)

    return 0

def eg_co2_o18_profile():
    # an example
    return 0

def eg_cos_profile():
    # an example
    return 0

def chamber_conc_solver():
    return 0

def flux_solver():
    return 0