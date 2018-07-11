import copy
import numpy as np
from soil_tracers_config import *

from tracer_sources import *
from create_grid import *
import scipy.sparse as spm  # scipy sparse matrix

import os
from gas_diffusivity import *
from solub import *

z_grid = FVGrid1D()
cyl_FVGrid = FVGridCyl2D()

def FV_CN_solver_1D(conc, alpha, beta, gamma, dt, flag_advection=False, 
    adv_vel_plus=0., adv_vel_minus=0., test_mode=False):
    ''' Finite-Volume Crank-Nicolson solver
    '''
    A_coef_mat = np.diag(alpha, 0)
    A_coef_mat_inv = np.diag(1./alpha, 0)
    
    beta_prime = beta[np.r_[1:len(beta),0]]
    # beta_prime = beta
    # beta_prime[0:-1] = beta[1:]
    beta_prime[-1] = 0.
    B_coef_mat = - np.diag(beta) - np.diag(beta_prime)
    for loop_num in range(1,len(beta)):
        B_coef_mat[loop_num, loop_num-1] = beta[loop_num]
        B_coef_mat[loop_num-1, loop_num] = beta[loop_num]

    S_coef_mat = gamma

    if test_mode:
        print (adv_vel_plus)
        print (adv_vel_minus)

    if flag_advection:
        V_plus_mat = np.zeros(( len(adv_vel_plus), len(adv_vel_plus) ))
        V_plus_mat[0,0] = (adv_vel_plus[0] - adv_vel_plus[1]) * 0.5
        V_plus_mat[0,1] = - adv_vel_plus[1] * 0.5
        V_plus_mat[1,0] = (adv_vel_plus[1] + adv_vel_plus[2]) * 0.5
        V_plus_mat[1,1] = 0.5 * adv_vel_plus[1] - 1.5 * adv_vel_plus[2]
        for loop_num in range(2,len(adv_vel_plus)-1):
            V_plus_mat[loop_num,loop_num-2] = -0.5 * adv_vel_plus[loop_num]
            V_plus_mat[loop_num,loop_num-1] = 1.5 * adv_vel_plus[loop_num] + \
                0.5 * adv_vel_plus[loop_num+1]
            V_plus_mat[loop_num,loop_num] = -1.5 * adv_vel_plus[loop_num+1]
        V_plus_mat[-1,-3] = -0.5 * adv_vel_plus[-1]
        V_plus_mat[-1,-2] = 1.5 * adv_vel_plus[-1]

        V_minus_mat = np.zeros(( len(adv_vel_minus), len(adv_vel_minus) ))
        V_minus_mat[0,0] = 1.5 * adv_vel_minus[0]
        V_minus_mat[0,1] = -0.5 * adv_vel_minus[0] - 1.5 * adv_vel_minus[1]
        for loop_num in range(1,len(adv_vel_plus)-2):
            V_minus_mat[loop_num,loop_num] = 1.5 * adv_vel_minus[loop_num]
            V_minus_mat[loop_num,loop_num+1] = -0.5 * adv_vel_minus[loop_num] - \
                1.5 * adv_vel_minus[loop_num+1]
            V_minus_mat[loop_num,loop_num+2] = 0.5 * adv_vel_minus[loop_num+1]
        V_minus_mat[-2,-2] = 1.5 * adv_vel_minus[-2] - 0.5 * adv_vel_minus[-1]
        V_minus_mat[-2,-1] = -0.5 * adv_vel_minus[-2] - 0.5 * adv_vel_minus[-1]
        V_minus_mat[-1,-2] = 0.5 * adv_vel_minus[-1]
        V_minus_mat[-1,-1] = 0.5 * adv_vel_minus[-1]

        B_coef_mat = B_coef_mat + V_plus_mat + V_minus_mat

    D_coef_mat = np.linalg.inv( 2 * A_coef_mat - dt * B_coef_mat )
    conc_updated = np.dot( D_coef_mat, 
        np.dot(2 * A_coef_mat + dt * B_coef_mat, conc) + 2 * dt * S_coef_mat )
    return conc_updated


def FV_RK4_solver(conc, alpha, beta, gamma, dt):
    '''Fourth order Runge-Kutta solver
    Warning: Do not use this solver!
    This method has poor stability and thus requires small time step. 
    '''

    A_coef_mat = np.diag(alpha)
    A_coef_mat_inv = np.diag(1./alpha)
    
    beta_prime = beta[np.r_[1:len(beta),0]]
    # beta_prime = beta
    # beta_prime[0:-1] = beta[1:]
    beta_prime[-1] = 0.
    B_coef_mat = - np.diag(beta) - np.diag(beta_prime)
    for loop_num in range(1,len(beta)):
        B_coef_mat[loop_num, loop_num-1] = beta[loop_num]
        B_coef_mat[loop_num-1, loop_num] = beta[loop_num]

    S_coef_mat = gamma

    K1 = dt * np.dot( A_coef_mat_inv, 
        ( np.dot(B_coef_mat, conc) + S_coef_mat ) )
    K2 = dt * np.dot( A_coef_mat_inv, 
        ( np.dot(B_coef_mat, conc + K1/2.) + S_coef_mat ) )
    K3 = dt * np.dot( A_coef_mat_inv, 
        ( np.dot(B_coef_mat, conc + K2/2.) + S_coef_mat ) )
    K4 = dt * np.dot( A_coef_mat_inv, 
        ( np.dot(B_coef_mat, conc + K3) + S_coef_mat ) )

    #conc_updated = conc + dt * np.dot( A_coef_mat_inv, 
    #   ( np.dot(B_coef_mat, conc) + S_coef_mat ) )
    conc_updated = conc + K1/6. + K2/3. + K3/3. + K4/6.

    return conc_updated

def FV_CN_cyl2D_solver(conc, alpha, beta, gamma, dt):
    return 0.

def flux_solver_1D(species, temp, poros, water_cont, atm_conc=None, air_temp=None, texture='sandy loam', 
    n_lit_layers=0, lit_poros=None, lit_water_cont=None, adv_vel_g=0., adv_vel_aq=0., source_profile=None, 
    f_CA=None, f_biomass=None, custom_src_params=None, 
    init_conc_profile_g=None, init_conc_profile_aq=None, FV_grid=z_grid, flux_solver_options=None):
    '''
    Evaluate flux of a tracer (in 1D)

    Parameters
    ----------
    species: string
        The gas species of interest, currently supports 'co2' and 'cos'
    temp: float
        Soil temperature profile, in degree C
    poros: float
        Soil porosity profile, in m3 m-3
    water_cont: float
        Soil water content profile, in m3 m-3
        When litter is enabled, water_cont in the litter layers represent the volumetric water content in 
        litter pore space, not structural water of the litter.

    (optional)
    atm_conc: float
        Atmospheric concentration (as boundary condition), in mol m-3
        By default, this is set to the mean atmospheric concentration of the species
    air_temp: float
        Air temperature in C. If not specified, air temperature equals to the surface soil temperature.
    texture: string
        Soil texture, used in calculating diffusivities.
    n_lit_layers: integer
        Number of litter layers in current vertical grid settings. If this is set to a valid number N, the top N layers 
        of temperature, porosity, and water content profiles, and other optional profiles (if specified) 
        will be considered as those of litter layers. 
    lit_poros: float
        Litter porosity, in m3 m^-3
    lit_water_cont: float
        Litter structural water content, in mass fraction
    adv_vel_g: float
        Advection velocity of the gaseous phase, in m s-1. Note a positive value means downward transport.
    adv_vel_aq: float
        Advection velocity of the aqueous phase, in m s-1. Note a positive value means downward transport.
    source_profile: float
        Specify the source profile, otherwise the default source functions are called
    f_CA: float
        Ratio of the enhancement of CO2 and COS hydrolysis compared with abiotic reactions (in development)
        This variable scales COS sink strength, and 13C-CO2 kinetic fractionation
    f_biomass: float 
        A scaling factor to scale source-sink terms when they depend on microbial biomass (in development)
    custom_src_params: dictionary
       Custom parameters for microbial sources and sinks 
    init_conc_profile_g, init_conc_profile_aq: float
        Initial concentration profile of gaseous and aqueous phases. If this is not specified, 
        the atmospheric concentration will be used to construct constant profiles. 
    FV_grid: vertical grid, class FVGrid1D
    flux_solver_options: dictionary of flux solver options
        'dt': float
            Time step in sec, by default 1 sec. Use large values if not starting from a steady state.
        'tol': float
            Tolerance, by default 1e-6
        'max_iter': long int
            Maximum number of iterations to calculate for the steady state, by default 1e5.
        'flag_advection': boolean
            By default 'False'. Set it 'True' to include advection effect in calculation
        'flag_chamber': boolean
            If false, the top boundary condition does not change. If true, the top boundary is set 
            as a chamber, and the boundary condition evolves in time. Flow rate needs to be specified 
            at the same time, otherwise an error is raised. Chamber head space concentration will be returned
        'chamber_vol': float
            Chamber volume in m3. The default chamber volume is specified in soil_tracers_config.py
        'chamber_area': float
            Chamber-covered surface area in m2. The default chamber area is specified in soil_tracers_config.py
        'flow_rate': float
            Flow rate through the chamber, in m3 per sec
        'flag_return_conc_profile': boolean
            If True, return both surface flux value and the final concentration profile
        'flag_steady_state': boolean
            Evaluate until a steady state is attained. Set this to 'False' if you want to calculate a transient state.
        'duration': float
            This specifies the duration of calculation for a transient state. This variable is only used when 
            steady_state = False. 
        'flag_solub_eq': boolean
            If True, assume aqueous concentration is always in solubility equilibrium with gaseous concentration.
            If False, perform prognostic calculation of aqueous concentration.
        'flag_ebullition': boolean
            If True, allows aqueous gas transport by bubbling (in development)
        'flag_test_mode': boolean
            This is for the developer to test the code, not useful for a user.

    Examples
    --------
    # no advection
    >>> flux_solver_1D('co2', 25, 0.45, 0.2, flux_solver_options = dict([ ('dt', 60) ])  )
    1.987548610469133e-06
    
    # downward advection (into the soil)
    >>> flux_solver_1D('co2', 25, 0.45, 0.2, adv_vel_g=1e-5, flux_solver_options = dict([ ('dt', 60), ('flag_advection', True), ]) )
    3.9080759580462365e-07
    
    # upward advection (out from the soil)
    >>> flux_solver_1D('co2', 25, 0.45, 0.2, adv_vel_g=-1e-5, flux_solver_options = dict([ ('dt', 60), ('flag_advection', True), ]) )
    1.2365140908145492e-05
    
    # advection with non-constant Darcian velocity (but non-time-varying) 
    >>> flux_solver_1D('co2', 25, 0.5, 0.3, adv_vel_g=-np.arange(z_grid.n_level, 0, -1)*1e-5, flux_solver_options = dict([ ('dt', 60), ('flag_advection', True), ]) )
    3.1299313190814474e-06
    
    # customized source parameters
    >>> test_source_params = copy.copy(src_params)
    >>> test_source_params['soil_resp_swc_dependence'] = True
    >>> flux_solver_1D('co2', 25, 0.45, 0.2, custom_src_params=test_source_params, flux_solver_options = dict([ ('dt', 60.) ])  )
    1.5021164383201848e-05
    
    # for COS
    >>> flux_solver_1D('cos', 15, 0.45, 0.2, flux_solver_options = dict([ ('dt', 5.) ])  )
    -2.1746086267305549e-12
    
    # with chamber
    # COS
    >>> flux, init_conc_profile = flux_solver_1D('cos', 15, 0.45, 0.2, 
    flux_solver_options = dict([ ('dt', 5.), ('flag_return_conc_profile', True) ]) )
    >>> ch_flux, ch_conc_profile, ch_conc = flux_solver_1D('cos', 15, 0.45, 0.2, init_conc_profile_g = init_conc_profile, 
    flux_solver_options = dict([ ('dt', 5.), ('flag_chamber', True), ('flag_steady_state', False), 
    ('duration', 1200), ('flag_return_conc_profile', True), ('flow_rate', 2e-5) ])  )
    # CO2
    >>> flux, init_conc_profile = flux_solver_1D('co2', 25, 0.45, 0.2, 
    flux_solver_options = dict([ ('dt', 60.), ('flag_return_conc_profile', True) ]) )
    >>> ch_flux, ch_conc_profile, ch_conc = flux_solver_1D('co2', 25, 0.45, 0.2, init_conc_profile_g = init_conc_profile, 
    flux_solver_options = dict([ ('dt', 60.), ('flag_chamber', True), ('flag_steady_state', False), 
    ('duration', 1200), ('flag_return_conc_profile', True), ('flow_rate', 2e-5) ])  )
    '''
    # ------ flux solver options ------
    # default flux solver options
    default_flux_solver_options = dict([ ('dt', 1.), ('tol', 1e-6), ('max_iter', 100000L),
        ('flag_advection', False), ('flag_chamber', False), ('chamber_vol', ch_params['volume']),
        ('chamber_area', ch_params['area']), ('flow_rate', np.nan), 
        ('flag_steady_state', True), ('duration', float('Inf')), 
        ('flag_solub_eq', True), ('flag_ebullition', False), 
        ('flag_return_conc_profile', False), 
        ('flag_test_mode', False), ])
    # initialize the flux solver options
    if flux_solver_options==None:
        flux_solver_options = copy.copy(default_flux_solver_options)
    else:
        new_flux_solver_options = copy.copy(default_flux_solver_options)
        new_flux_solver_options.update(flux_solver_options)
        flux_solver_options = copy.copy(new_flux_solver_options)
        del(new_flux_solver_options)
    # if not steady state, set the max iteration number from duration
    if not flux_solver_options['flag_steady_state']:
        flux_solver_options['max_iter'] = np.int(np.ceil( flux_solver_options['duration'] / flux_solver_options['dt'] ))
    # update some runtime parameters
    dt = flux_solver_options['dt']; tol = flux_solver_options['tol']; max_iter = flux_solver_options['max_iter']
    duration = flux_solver_options['duration']

    # if chamber is present, check if necessary parameters are given
    if flux_solver_options['flag_chamber']:
        if not (np.isfinite(flux_solver_options['chamber_vol']) and np.isfinite(flux_solver_options['chamber_area']) 
          and np.isfinite(flux_solver_options['flow_rate']) and np.isfinite(flux_solver_options['duration'])):
            raise ValueError('Chamber-mode error: At least one chamber parameter is not properly given.')
        if flux_solver_options['flag_steady_state']:
            raise RuntimeError('Chamber-mode error: Steady-state option must be turned off.')
        if not ((np.sum(init_conc_profile_g) >= 0.) or (np.sum(init_conc_profile_aq) >= 0.)):
            raise ValueError('Chamber-mode error: Initial concentration profile (gaseous or aqueous) must be given.')
        # create a variable for saving head space concentration
        ch_conc = np.zeros(flux_solver_options['max_iter'] + 1) * np.nan

    # give warning message when using steady-state option for advection mode
    if flux_solver_options['flag_steady_state'] and flux_solver_options['flag_advection']:
        print('Warning: Steady-state calculation is problematic in the advection mode, because zero-flux bottom boundary condition is used. ')
        print('To ensure a good result please use non-steady-state options.')

    # raise error messages for unsupported options
    if (not flux_solver_options['flag_solub_eq']):
        raise RuntimeError('Unsupported option: dynamic solubility equilibrium')
    if flux_solver_options['flag_ebullition']:
        raise RuntimeError('Unsupported option: ebullition')

    # initialize the microbial source-sink parameters
    if custom_src_params == None:
        custom_src_params = copy.copy(src_params)
    # ------ end of flux solver options ------

    # ------ initialize environmental variables ------
    # force the input profiles to be vectors
    temp = np.array(temp, dtype='d') * np.ones(FV_grid.n_level)
    poros = np.array(poros, dtype='d') * np.ones(FV_grid.n_level)
    water_cont = np.array(water_cont, dtype='d') * np.ones(FV_grid.n_level)
    air_poros = poros - water_cont  # air-filled porosity

    if air_temp == None:
        air_temp = temp[0]
    else:
        air_temp = np.array(air_temp, dtype='d')

    # initialize top boundary condition 
    if atm_conc == None:
        atm_conc = mean_conc[species] * air_conc_std * 273.15 / (273.15 + air_temp)
    else:
        atm_conc = np.array(atm_conc, dtype='d')

    # if leaf litter exists, add it
    n_lit_layers = np.int(n_lit_layers)
    if n_lit_layers > 0 and n_lit_layers <= FV_grid.n_level:
        poros[0:n_lit_layers] = lit_poros
        water_cont[0:n_lit_layers] = 0.  # pore-space water content of litter layers is zero
        air_poros = poros - water_cont
    elif not np.equal(n_lit_layers, 0):
        print('Wrong number of litter layers. Will be set 0.')
        n_lit_layers = 0

    # calculate diffusivity profiles
    D_air = diffus_in_air(species, air_temp)
    D_soil_g = diffus_in_soil_air(species, temp, poros, air_poros, soil_texture=texture)
    D_soil_aq = diffus_in_soil_aq(species, temp, poros, air_poros, soil_texture=texture)

    # calculate Bunsen solubility
    if species == 'co2':
        bunsen_solub = solub_co2(temp)
    elif species == 'cos':
        bunsen_solub = solub_cos(temp)

    # initialize concentration profile, gaseous phase
    if init_conc_profile_g == None:
        soil_conc_g = atm_conc * np.ones(FV_grid.n_level)
    else:
        soil_conc_g = init_conc_profile_g * np.ones(FV_grid.n_level)
    # ------ end of the initialization of environmental variables ------

    '''
    Branching
    If solubility equilibrium is turned on (default), evaluate only gaseous concentration profiles
    If solubility equilibrium is turned off, evaluate both gaseous and aqueous concentration profiles
    '''
    if flux_solver_options['flag_solub_eq']: 
        alpha = (air_poros + bunsen_solub * water_cont) * FV_grid.cv_size
        beta = (D_soil_g + D_soil_aq) * np.ones(FV_grid.n_level)
        # beta[0] = (beta[0] + D_air) / 2
        beta[0] = 2 * beta[0] * D_air / (beta[0] + D_air) 
        # use the harmonic mean for soil-air boundary
        beta[1:] = (beta[1:] + beta[0:-1]) / 2
        beta[0] = beta[0] / FV_grid.grid_node[0]
        beta[1:] = beta[1:] / (FV_grid.grid_node[1:] - FV_grid.grid_node[0:-1])
        # for upwind-biased scheme of advective flux evaluation
        adv_vel_g_interface = adv_vel_g * np.ones(FV_grid.n_level)
        adv_vel_g_interface[1:] = (adv_vel_g_interface[1:] + adv_vel_g_interface[0:-1]) / 2. 
        adv_vel_g_plus = np.maximum(adv_vel_g_interface, 0.)
        adv_vel_g_minus = np.minimum(adv_vel_g_interface, 0.)
        
        # gamma = np.zeros(FV_grid.n_level)
        if source_profile == None:
            if species == 'co2':
                gamma = soil_co2_source(temp, water_cont, FV_grid.grid_node, 
                    vmax = custom_src_params['soil_co2_source_vmax'], 
                    decay_depth = custom_src_params['soil_resp_decay_depth'], 
                    depth_func_type = custom_src_params['soil_resp_depth_func_type'], 
                    swc_dependence = custom_src_params['soil_resp_swc_dependence'], 
                    k_w = custom_src_params['soil_resp_water_factor']) * FV_grid.cv_size
                # if leaf litter exists
                if n_lit_layers > 0:  
                    gamma[0:n_lit_layers] = litter_co2_source(temp[0:n_lit_layers], water_cont=lit_water_cont * np.ones(n_lit_layers), 
                        vmax=custom_src_params['litter_co2_source_vmax'], E_resp=custom_src_params['litter_resp_E'], 
                        k_w=custom_src_params['litter_resp_water_factor'], 
                        lwc_dependence=custom_src_params['litter_resp_swc_dependence']) * FV_grid.cv_size[0:n_lit_layers]
                    gamma[n_lit_layers:] = soil_co2_source(temp[n_lit_layers:], water_cont[n_lit_layers:], 
                        FV_grid.grid_node[n_lit_layers:] - FV_grid.cv_face[n_lit_layers-1], 
                        vmax = custom_src_params['soil_co2_source_vmax'], 
                        decay_depth = custom_src_params['soil_resp_decay_depth'], 
                        depth_func_type = custom_src_params['soil_resp_depth_func_type'], 
                        swc_dependence = custom_src_params['soil_resp_swc_dependence'], 
                        k_w = custom_src_params['soil_resp_water_factor']) * FV_grid.cv_size[n_lit_layers:]
                    # note, the depth dependent CO2 source function should have an offset for litter layer thickness
            elif species == 'cos':
                gamma = soil_cos_sink(temp, water_cont, soil_conc_g, 
                    vmax = custom_src_params['soil_cos_sink_vmax'], 
                    K_m = custom_src_params['soil_cos_sink_Km'], 
                    Delta_G_cat = custom_src_params['soil_cos_sink_Delta_G_cat'],
                    Delta_H_eq = custom_src_params['soil_cos_sink_Delta_H_eq'], 
                    T_eq = custom_src_params['soil_cos_sink_T_eq'], 
                    swc_opt = custom_src_params['soil_cos_sink_swc_opt'],
                    swc_dependence = custom_src_params['soil_cos_sink_swc_dependence']) * FV_grid.cv_size
                gamma += soil_cos_source(temp, vmax = custom_src_params['soil_cos_source_vmax'], 
                    q10 = custom_src_params['soil_cos_source_q10']) * FV_grid.cv_size
                # if leaf litter exists
                if n_lit_layers > 0:  
                    gamma[0:n_lit_layers] = litter_cos_sink(temp[0:n_lit_layers], 
                        lit_water_cont * np.ones(n_lit_layers), soil_conc_g[0:n_lit_layers], 
                        vmax=custom_src_params['litter_cos_sink_vmax'], K_m=custom_src_params['litter_cos_sink_Km'], 
                        k_w=custom_src_params['litter_cos_sink_water_factor']) * FV_grid.cv_size[0:n_lit_layers]
        else: 
            gamma += source_profile * FV_grid.cv_size
            
        gamma[0] = gamma[0] + beta[0] * atm_conc
        gamma[0] += 0.5 * adv_vel_g_plus[0] * atm_conc

        # initialize chamber headspace concentration, when chamber option is enabled
        if flux_solver_options['flag_chamber']:
            ch_conc[0] = atm_conc

        for loop_num in range(max_iter):
            # recalculate the source-sink term for COS every time when the conc profile is updated
            if species == 'cos' and source_profile == None:
                gamma = np.zeros(FV_grid.n_level)
                gamma += soil_cos_sink(temp, water_cont, soil_conc_g, 
                    vmax = custom_src_params['soil_cos_sink_vmax'], 
                    K_m = custom_src_params['soil_cos_sink_Km'], 
                    Delta_G_cat = custom_src_params['soil_cos_sink_Delta_G_cat'],
                    Delta_H_eq = custom_src_params['soil_cos_sink_Delta_H_eq'], 
                    T_eq = custom_src_params['soil_cos_sink_T_eq'], 
                    swc_opt = custom_src_params['soil_cos_sink_swc_opt'],
                    swc_dependence = custom_src_params['soil_cos_sink_swc_dependence']) * FV_grid.cv_size
                gamma += soil_cos_source(temp, vmax = custom_src_params['soil_cos_source_vmax'], 
                    q10 = custom_src_params['soil_cos_source_q10']) * FV_grid.cv_size
                # if leaf litter exists
                if n_lit_layers > 0:  
                    gamma[0:n_lit_layers] = litter_cos_sink(temp[0:n_lit_layers], 
                        lit_water_cont * np.ones(n_lit_layers), soil_conc_g[0:n_lit_layers], 
                        vmax=custom_src_params['litter_cos_sink_vmax'], K_m=custom_src_params['litter_cos_sink_Km'], 
                        k_w=custom_src_params['litter_cos_sink_water_factor']) * FV_grid.cv_size[0:n_lit_layers]
                
                gamma[0] = gamma[0] + beta[0] * atm_conc
                gamma[0] += 0.5 * adv_vel_g_plus[0] * atm_conc
            
            # if chamber option is enabled, change the top boundary condition
            if flux_solver_options['flag_chamber']:
                if species == 'cos' and source_profile == None:
                    gamma[0] += beta[0] * (ch_conc[loop_num] - atm_conc) + 0.5 * adv_vel_g_plus[0] * (ch_conc[loop_num] - atm_conc)
                elif species == 'co2' and source_profile == None and loop_num > 0:
                    gamma[0] += beta[0] * (ch_conc[loop_num] - ch_conc[loop_num-1]) + \
                    0.5 * adv_vel_g_plus[0] * (ch_conc[loop_num] - ch_conc[loop_num-1])
                elif source_profile != None:
                    gamma = source_profile * FV_grid.cv_size
                    gamma[0] = gamma[0] + beta[0] * ch_conc[loop_num] + 0.5 * adv_vel_g_plus[0] * ch_conc[loop_num]

            soil_conc_g_updated = FV_CN_solver_1D(soil_conc_g, alpha, beta, gamma, dt, 
                flag_advection = flux_solver_options['flag_advection'], 
                adv_vel_plus = adv_vel_g_plus, adv_vel_minus = adv_vel_g_minus)
            soil_conc_g_updated[-1] = soil_conc_g_updated[-2] # rectify the bottom boundary condition
            
            # if chamber = True, update ch_conc
            # assuming well-mixed chamber air
            if flux_solver_options['flag_chamber']:
                surf_flux = beta[0] * (soil_conc_g[0] - ch_conc[loop_num]) - \
                0.5 * adv_vel_g_plus[0] * (ch_conc[loop_num] + soil_conc_g[0]) - \
                adv_vel_g_minus[0] * (1.5 * soil_conc_g[0] - 0.5 * soil_conc_g[1])
                '''
                # Euler scheme
                ch_conc[loop_num+1] = ch_conc[loop_num] + ( (atm_conc - ch_conc[loop_num]) * \
                flux_solver_options['flow_rate'] + surf_flux * flux_solver_options['chamber_area'] ) * \
                dt / flux_solver_options['chamber_vol']
                '''
                # trapezoidal scheme, which has better stability performance
                ch_conc[loop_num+1] = (ch_conc[loop_num] + ((atm_conc - 0.5 * ch_conc[loop_num]) * flux_solver_options['flow_rate'] + \
                    surf_flux * flux_solver_options['chamber_area'] ) * dt / flux_solver_options['chamber_vol']) / \
                    (1 + 0.5 * flux_solver_options['flow_rate'] * dt / flux_solver_options['chamber_vol'])

                del(surf_flux)

            if loop_num == int(100/dt) and flux_solver_options['flag_test_mode']:
                print 'test'
                print soil_conc_g_updated 
            if (np.linalg.norm(soil_conc_g_updated[0:-1] - soil_conc_g[0:-1]) / np.linalg.norm(soil_conc_g[0:-1]) < tol) and (loop_num > 2):
                soil_conc_g = soil_conc_g_updated
                if flux_solver_options['flag_test_mode']:
                    print loop_num
                break
            soil_conc_g = soil_conc_g_updated

        if flux_solver_options['flag_test_mode']: print soil_conc_g

        # print ('Number of loops executed: ', loop_num)
        if (flux_solver_options['flag_steady_state'] and loop_num+1 >= flux_solver_options['max_iter']):
            print('Warning: Maximum number of iterations reached, still not at steady state.')

        if not flux_solver_options['flag_chamber']:
            surf_flux = beta[0] * (soil_conc_g[0] - atm_conc)
            if flux_solver_options['flag_advection']:
                surf_flux = surf_flux - 0.5 * adv_vel_g_plus[0] * (atm_conc + soil_conc_g[0]) - \
                adv_vel_g_minus[0] * (1.5 * soil_conc_g[0] - 0.5 * soil_conc_g[1])
            if flux_solver_options['flag_return_conc_profile']:
                return surf_flux, soil_conc_g
            else:
                return surf_flux
        else:
            surf_flux = beta[0] * (soil_conc_g[0] - ch_conc[loop_num+1])
            if flux_solver_options['flag_advection']:
                surf_flux = surf_flux - 0.5 * adv_vel_g_plus[0] * (ch_conc[loop_num+1] + soil_conc_g[0]) - \
                adv_vel_g_minus[0] * (1.5 * soil_conc_g[0] - 0.5 * soil_conc_g[1])
            if flux_solver_options['flag_return_conc_profile']:
                return surf_flux, soil_conc_g, ch_conc
            else:
                return surf_flux, ch_conc
        # fin de la branche
    else: 
        # initialize concentration profile, aqueous phase
        if init_conc_profile_aq == None:
            soil_conc_aq = soil_conc_g * bunsen_solub
        else:
            soil_conc_aq = init_conc_profile_aq * np.ones(FV_grid.n_level)

        alpha_g = air_poros * FV_grid.cv_size
        beta_g = D_soil_g * np.ones(FV_grid.n_level) # force to be vector
        beta_g[0] = 2 * beta_g[0] * D_air / (beta_g[0] + D_air) # use the harmonic mean for soil-air boundary
        beta_g[1:] = (beta_g[1:] + beta_g[0:-1]) / 2
        beta_g[0] = beta_g[0] / FV_grid.grid_node[0]
        beta_g[1:] = beta_g[1:] / (FV_grid.grid_node[1:] - FV_grid.grid_node[0:-1])

        alpha_aq = water_cont * FV_grid.cv_size
        beta_aq = D_soil_aq * np.ones(FV_grid.n_level) # force to be vector
        beta_aq[0] = 2 * beta_aq[0] * D_air / (beta_aq[0] + D_air) # in development
        # temporarily, use the harmonic mean for soil-air boundary; a proper treatment could employ the water-film analogous model
        beta_aq[1:] = (beta_aq[1:] + beta_aq[0:-1]) / 2
        beta_aq[0] = beta_aq[0] / FV_grid.grid_node[0]
        beta_aq[1:] = beta_aq[1:] / (FV_grid.grid_node[1:] - FV_grid.grid_node[0:-1])
        
        # for upwind-biased scheme of advective flux evaluation, gaseous phase
        adv_vel_g_interface = adv_vel_g * np.ones(FV_grid.n_level)
        adv_vel_g_interface[1:] = (adv_vel_g_interface[1:] + adv_vel_g_interface[0:-1]) / 2. 
        adv_vel_g_plus = np.maximum(adv_vel_g_interface, 0.)
        adv_vel_g_minus = np.minimum(adv_vel_g_interface, 0.)

        # for upwind-biased scheme of advective flux evaluation, aqueous phase
        adv_vel_aq_interface = adv_vel_aq * np.ones(FV_grid.n_level)
        adv_vel_aq_interface[1:] = (adv_vel_aq_interface[1:] + adv_vel_aq_interface[0:-1]) / 2. 
        adv_vel_aq_plus = np.maximum(adv_vel_aq_interface, 0.)
        adv_vel_aq_minus = np.minimum(adv_vel_aq_interface, 0.)

        # currently assume a linear partitioning of source/sink terms between the two phases
        # i.e. gas-phase source/sink strength is scaled by the air-filled pore space fraction
        # a proper treatment may involve 3D probabilistic distribution of water in soil, which is way too complicated and 
        # is subject to perpetual development
        gamma_g = np.zeros(FV_grid.n_level)
        gamma_aq = np.zeros(FV_grid.n_level)
        if source_profile == None:
            if species == 'co2':
                gamma_g += soil_co2_source(temp, water_cont, FV_grid.grid_node, 
                    vmax = custom_src_params['soil_co2_source_vmax'], 
                    decay_depth = custom_src_params['soil_resp_decay_depth'], 
                    depth_func_type = custom_src_params['soil_resp_depth_func_type'], 
                    swc_dependence = custom_src_params['soil_resp_swc_dependence'], 
                    k_w = custom_src_params['soil_resp_water_factor']) * FV_grid.cv_size
                gamma_aq = gamma_g
                gamma_g *= air_poros / poros
                gamma_aq *= water_cont / poros
            elif species == 'cos':
                gamma_g += soil_cos_sink(temp, water_cont, soil_conc_g, 
                    vmax = custom_src_params['soil_cos_sink_vmax'], 
                    K_m = custom_src_params['soil_cos_sink_Km'], 
                    Delta_G_cat = custom_src_params['soil_cos_sink_Delta_G_cat'],
                    Delta_H_eq = custom_src_params['soil_cos_sink_Delta_H_eq'], 
                    T_eq = custom_src_params['soil_cos_sink_T_eq'], 
                    swc_opt = custom_src_params['soil_cos_sink_swc_opt'],
                    swc_dependence = custom_src_params['soil_cos_sink_swc_dependence']) * FV_grid.cv_size
                gamma_g += soil_cos_source(temp, vmax = custom_src_params['soil_cos_source_vmax'], 
                    q10 = custom_src_params['soil_cos_source_q10']) * FV_grid.cv_size
                gamma_aq = gamma_g
                gamma_g *= air_poros / poros
                gamma_aq *= water_cont / poros
        else: 
            gamma_g += source_profile * FV_grid.cv_size
            gamma_aq = gamma_g
            gamma_g *= air_poros / poros
            gamma_aq *= water_cont / poros
            
        gamma_g[0] = gamma_g[0] + beta_g[0] * atm_conc
        gamma_g[0] += 0.5 * adv_vel_g_plus[0] * atm_conc
        # Solubility factor needs to be considered in the boundary condition of the aqueous phase
        gamma_aq[0] = gamma_aq[0] + beta_aq[0] * bunsen_solub[0] * atm_conc 
        gamma_aq[0] += 0.5 * adv_vel_aq_plus[0] * bunsen_solub[0] * atm_conc

        for loop_num in range(max_iter):
            # recalculate the source-sink term for COS every time when the conc profile is updated
            if species == 'cos' and source_profile == None:
                gamma_g = np.zeros(FV_grid.n_level)
                gamma_g += soil_cos_sink(temp, water_cont, soil_conc_g, 
                    vmax = custom_src_params['soil_cos_sink_vmax'], 
                    K_m = custom_src_params['soil_cos_sink_Km'], 
                    Delta_G_cat = custom_src_params['soil_cos_sink_Delta_G_cat'],
                    Delta_H_eq = custom_src_params['soil_cos_sink_Delta_H_eq'], 
                    T_eq = custom_src_params['soil_cos_sink_T_eq'], 
                    swc_opt = custom_src_params['soil_cos_sink_swc_opt'],
                    swc_dependence = custom_src_params['soil_cos_sink_swc_dependence']) * FV_grid.cv_size
                gamma_g += soil_cos_source(temp, vmax = custom_src_params['soil_cos_source_vmax'], 
                    q10 = custom_src_params['soil_cos_source_q10']) * FV_grid.cv_size
                gamma_aq = gamma_g
                gamma_g *= air_poros / poros
                gamma_aq *= water_cont / poros
                gamma_g[0] = gamma_g[0] + beta_g[0] * atm_conc
                gamma_g[0] += 0.5 * adv_vel_g_plus[0] * atm_conc
                # Solubility factor needs to be considered in the boundary condition of the aqueous phase
                gamma_aq[0] = gamma_aq[0] + beta_aq[0] * bunsen_solub[0] * atm_conc 
                gamma_aq[0] += 0.5 * adv_vel_aq_plus[0] * bunsen_solub[0] * atm_conc
            
            soil_conc_g_updated = FV_CN_solver_1D(soil_conc_g, alpha_g, beta_g, gamma_g, dt, 
                flag_advection = flux_solver_options['flag_advection'], 
                adv_vel_plus = adv_vel_g_plus, adv_vel_minus = adv_vel_g_minus)
            soil_conc_g_updated[-1] = soil_conc_g_updated[-2] # rectify the bottom boundary condition
            soil_conc_aq_updated = FV_CN_solver_1D(soil_conc_aq, alpha_aq, beta_aq, gamma_aq, dt, 
                flag_advection = flux_solver_options['flag_advection'], 
                adv_vel_plus = adv_vel_aq_plus, adv_vel_minus = adv_vel_aq_minus)
            soil_conc_aq_updated[-1] = soil_conc_aq_updated[-2] # rectify the bottom boundary condition
            # air-water gas exchange through solubility equilibrium

            # to be continued
            
            '''
            piston_vel = 0.39 * (1e-2)**2 # use the equation k_660 = 0.39 * U_10^2, assuming U_10 = 0.01 m/s for soil
            # temporary use only, in development
            soil_conc_g += (soil_conc_aq - bunsen_solub * soil_conc_g) * water_cont / air_poros
            soil_conc_aq += (bunsen_solub * soil_conc_g - soil_conc_aq) * air_poros / water_cont
            '''
            # if chamber = True, update atm_conc  # in development
            # recalculate gamma
            if loop_num == int(100/dt) and flux_solver_options['flag_test_mode']:
                print 'test'
                print soil_conc_g_updated 
            if (np.linalg.norm(soil_conc_g_updated[0:-1] - soil_conc_g[0:-1]) / np.linalg.norm(soil_conc_g[0:-1]) < tol) and (loop_num > 2):
                soil_conc_g = soil_conc_g_updated
                if flux_solver_options['flag_test_mode']:
                    print loop_num
                break
            soil_conc_g = soil_conc_g_updated

        if flux_solver_options['flag_test_mode']: print soil_conc_g

        surf_flux = beta[0] * (soil_conc_g[0] - atm_conc)
        if flux_solver_options['flag_advection']:
            surf_flux = surf_flux - 0.5 * adv_vel_g_plus[0] * (atm_conc + soil_conc_g[0]) - \
            adv_vel_g_minus[0] * (1.5 * soil_conc_g[0] - 0.5 * soil_conc_g[1])
        
        if flux_solver_options['flag_return_conc_profile']:
            return surf_flux, soil_conc_g
        else:
            return surf_flux


    

def flux_solver_old(species, temp, poros, water_cont, atm_conc=None, air_temp=None,
    texture='sandy loam', lit_thick=None, lit_water_cont=None, source_profile=None, 
    f_CA=1., f_biomass=1., advection=False, adv_velocity=0., custom_src_params=None,
    FV_grid=z_grid, dt=1., tol=1e-6, max_iter=100000L, steady_state=True, init_conc=None, 
    duration=float('Inf'), chamber=False, chamber_vol=ch_params['volume'], chamber_area=ch_params['area'], 
    flow_rate=float('NaN'), save_headspace_conc=False, save_profile=False, return_conc_only=False,
    test_mode=False):
    '''
    Evaluate flux of a tracer (in 1D)

    Examples
    --------
    # no advection
    >>> flux_solver('co2', 25, 0.45, 0.2, dt=60)
    1.9875484531494665e-06
    # downward advection (into the soil)
    >>> flux_solver('co2', 25, 0.45, 0.2, dt=60, advection=True, adv_velocity=1e-5)
    3.9078287933359341e-07
    # upward advection (out from the soil)
    >>> flux_solver('co2', 25, 0.45, 0.2, dt=60, advection=True, adv_velocity=-1e-5)
    1.2364958134273856e-05
    # advection with non-constant Darcian velocity (but non-time-varying) 
    >>> flux_solver('co2', 25, 0.5, 0.3, dt=60, advection=True, adv_velocity=-np.arange(z_grid.n_level, 0, -1)*1e-5)
    3.1298839561661173e-06
    # customized source parameters
    >>> test_source_params = copy.copy(src_params)
    >>> test_source_params['soil_resp_swc_dependence'] = True
    >>> flux_solver('co2', 25, 0.45, 0.2, dt=60, custom_src_params=test_source_params)
    1.5021163185756575e-05
    # for COS
    >>> flux_solver('cos', 15, 0.45, 0.2, dt=5)
    -2.1744953186421789e-12

    Parameters
    ----------
    species: string
        The gas species of interest, currently supports 'co2' and 'cos'
    temp: float
        Soil temperature profile, in Celsius
    poros: float
        Soil porosity profile, in m3 m-3
    water_cont: float
        Soil water content profile, in m3 m-3

    (optional)
    atm_conc: float
        Atmospheric concentration (as boundary condition), in mol m-3
        By default, this is set to the mean atmospheric concentration of the species
    air_temp: float
        Air temperature in C. If not specified, air temperature equals to the first element of the 
        soil temperature array.
    texture: string
        Soil texture, used in calculating diffusivities.
    lit_thick: float
        Thickness of litter layers in m
    lit_water_cont: float
        Litter water content in g g-1
    source_profile: float
        Specify the source profile, otherwise the default source functions are called
    f_CA: float
        Ratio of the enhancement of CO2 and COS hydrolysis compared with abiotic reactions
        This variable scales COS sink strength, and 13C-CO2 kinetic fractionation
    f_biomass: float
        A scaling factor to scale source-sink terms 
    advection: boolean
        By default 'False'. Set it 'True' to include advection effect in calculation
    adv_velocity: float
        Advection velocity in m s-1. Note a positive value means downward transport.
    custom_src_params: dictionary
       Custom parameters for microbial sources and sinks 
    FV_grid: vertical grid, class FVGrid1D
    dt: float
        Time step in sec, by default 1 sec. Use large values if not starting from a steady state.
    tol: float
        Tolerance, by default 1e-6
    max_iter: long int
        Maximum number of iterations to calculate for the steady state, by default 1e5.
    steady_state: boolean
        Evaluate until a steady state is attained. Set this to 'False' if you want to calculate
        a transient state.
    init_conc: float
        Initial concentration profile. If this is not specified, the atmospheric concentration 
        will be used to construct a constant profile. 
    duration: float
        This specifies the duration of calculation for a transient state. This variable is only
        used when steady_state = False. 
    chamber: boolean
        If false, the top boundary condition does not change. If true, the top boundary is set 
        as a chamber, and the boundary condition evolves in time. Flow rate needs to be specified 
        at the same time, otherwise an error is raised. 
    chamber_vol: float
        Chamber volume in m3. The default chamber volume is specified in soil_tracers_config.py
    chamber_area: float
        Chamber-covered surface area in m2. The default chamber area is specified in 
        soil_tracers_config.py
    flow_rate: float
        Flow rate through the chamber, in m3 per sec
    save_headspace_conc: boolean
        If True, return the headspace concentration as well as the flux
        This is only valid when steady_state = False, and np.isfinite(duration) = True
    save_profile: boolean
        Save the concentration profile at the end of the run
    test_mode: boolean
        This is for the developer to test the code, not useful for a user.
    '''
    
    # initialize the microbial source-sink parameters
    if custom_src_params == None:
        custom_src_params = copy.copy(src_params)
    
    # force the input profiles to be vectors
    temp = np.array(temp, dtype='d') * np.ones(FV_grid.n_level)
    poros = np.array(poros, dtype='d') * np.ones(FV_grid.n_level)
    water_cont = np.array(water_cont, dtype='d') * np.ones(FV_grid.n_level)
    air_poros = poros - water_cont  # air-filled porosity

    if air_temp == None:
        air_temp = temp[0]
    else:
        air_temp = np.array(air_temp, dtype='d')

    D_air = diffus_in_air(species, air_temp)
    D_soil = diffus_in_soil_air(species, temp, poros, air_poros, soil_texture=texture)
    if species == 'co2':
        bunsen_solub = solub_co2(temp)
    elif species == 'cos':
        bunsen_solub = solub_cos(temp)

    if atm_conc == None:
        atm_conc = mean_conc[species] * air_conc_std * 273.15 / (273.15 + air_temp)
    else:
        atm_conc = np.array(atm_conc, dtype='d')

    if init_conc == None:
        soil_conc = atm_conc * np.ones(FV_grid.n_level)
    else:
        soil_conc = init_conc * np.ones(FV_grid.n_level)

    alpha = (air_poros + bunsen_solub * water_cont) * FV_grid.cv_size

    beta = D_soil * np.ones(FV_grid.n_level)
    #beta[0] = (beta[0] + D_air) / 2
    beta[0] = 2 * beta[0] * D_air / (beta[0] + D_air) 
    # use the harmonic mean for soil-air boundary
    beta[1:] = (beta[1:] + beta[0:-1]) / 2
    beta[0] = beta[0] / FV_grid.grid_node[0]
    beta[1:] = beta[1:] / (FV_grid.grid_node[1:] - FV_grid.grid_node[0:-1])
    
    adv_vel_interface = adv_velocity * np.ones(FV_grid.n_level)
    adv_vel_interface[1:] = (adv_vel_interface[1:] + adv_vel_interface[0:-1]) / 2. 
    adv_vel_plus = np.maximum(adv_vel_interface, 0.)
    adv_vel_minus = np.minimum(adv_vel_interface, 0.)
    
    gamma = np.zeros(FV_grid.n_level)
    if source_profile == None:
        if species == 'co2':
            gamma += soil_co2_source(temp, water_cont, FV_grid.grid_node, 
                vmax=custom_src_params['soil_co2_source_vmax'], 
                decay_depth=custom_src_params['soil_resp_decay_depth'], 
                depth_func_type=custom_src_params['soil_resp_depth_func_type'], 
                swc_dependence=custom_src_params['soil_resp_swc_dependence'], 
                k_w=custom_src_params['soil_resp_water_factor']) * FV_grid.cv_size
        elif species == 'cos':
            gamma += soil_cos_sink(temp, water_cont, soil_conc, 
                vmax=custom_src_params['soil_cos_sink_vmax'], 
                K_m=custom_src_params['soil_cos_sink_Km'], 
                Delta_G_cat=custom_src_params['soil_cos_sink_Delta_G_cat'],
                Delta_H_eq=custom_src_params['soil_cos_sink_Delta_H_eq'], 
                T_eq=custom_src_params['soil_cos_sink_T_eq'], 
                swc_opt=custom_src_params['soil_cos_sink_swc_opt'],
                swc_dependence=custom_src_params['soil_cos_sink_swc_dependence']) * FV_grid.cv_size
            gamma += soil_cos_source(temp, vmax=custom_src_params['soil_cos_source_vmax'], 
                q10=custom_src_params['soil_cos_source_q10']) * FV_grid.cv_size
    else: 
        gamma += source_profile * FV_grid.cv_size
        
    gamma[0] = gamma[0] + beta[0] * atm_conc
    gamma[0] += 0.5 * adv_vel_plus[0] * atm_conc    
    
    if not steady_state:
        max_iter = np.int(np.ceil(duration / dt))
    
    '''
    soil_conc_updated = FV_CN_solver_1D(soil_conc, alpha, beta, gamma, dt, 
        flag_advection=advection, adv_vel_plus=adv_vel_plus, adv_vel_minus=adv_vel_minus,
        test_mode=test_mode)
    print soil_conc_updated.shape
    '''
    
    
    for loop_num in range(max_iter):
        # recalculate the source-sink term for COS every time when the conc profile is updated
        if species == 'cos' and source_profile == None:
            gamma = np.zeros(FV_grid.n_level)
            gamma += soil_cos_sink(temp, water_cont, soil_conc, 
                vmax=custom_src_params['soil_cos_sink_vmax'], 
                K_m=custom_src_params['soil_cos_sink_Km'], 
                Delta_G_cat=custom_src_params['soil_cos_sink_Delta_G_cat'],
                Delta_H_eq=custom_src_params['soil_cos_sink_Delta_H_eq'], 
                T_eq=custom_src_params['soil_cos_sink_T_eq'], 
                swc_opt=custom_src_params['soil_cos_sink_swc_opt'],
                swc_dependence=custom_src_params['soil_cos_sink_swc_dependence']) * FV_grid.cv_size
            gamma += soil_cos_source(temp, vmax=custom_src_params['soil_cos_source_vmax'], 
                q10=custom_src_params['soil_cos_source_q10']) * FV_grid.cv_size            
            gamma[0] = gamma[0] + beta[0] * atm_conc
            gamma[0] += 0.5 * adv_vel_plus[0] * atm_conc
        
        soil_conc_updated = FV_CN_solver_1D(soil_conc, alpha, beta, gamma, dt, 
            flag_advection=advection, adv_vel_plus=adv_vel_plus, adv_vel_minus=adv_vel_minus)
        soil_conc_updated[-1] = soil_conc_updated[-2] # rectify the bottom boundary condition
        # if chamber = True, update atm_conc
        # recalculate gamma
        if loop_num == int(100/dt) and test_mode:
            print 'test'
            print soil_conc_updated 
        if ( np.linalg.norm(soil_conc_updated[0:-1] - soil_conc[0:-1]) / 
            np.linalg.norm(soil_conc[0:-1]) < tol) & (loop_num > 2):
            soil_conc = soil_conc_updated
            if test_mode:
                print loop_num
            break
        soil_conc = soil_conc_updated

    if test_mode:
        print soil_conc

    # test mass closure
    # source_func = np.append( [0.], 
    #   soil_co2_source(temp, water_cont, FV_grid.grid_node) )
    # depth_func = np.append([0.], FV_grid.grid_node)
    # print np.trapz(source_func, depth_func)

    surf_flux = beta[0] * (soil_conc[0] - atm_conc)
    if advection:
        surf_flux = surf_flux - 0.5 * adv_vel_plus[0] * (atm_conc + soil_conc[0]) - \
        adv_vel_minus[0] * (1.5 * soil_conc[0] - 0.5 * soil_conc[1])
    
    if save_profile:
        return surf_flux, soil_conc    
    
    if return_conc_only:
        return soil_conc
    else:
        return surf_flux
            
    #class FluxData():
    #return flux


''' to test this routine
from flux_solver import *
flux_solver('co2', 25, 0.45, 0.2, atm_conc=air_conc_std * 4e-4, dt=60)
'''

def flux_solver_cyl2D(species, temp, total_poros, water_cont, air_temp=None,
    litter_thickness=0., litter_water_cont=float('NaN'), source_profile=None, atm_conc=0., 
    ch_conc_init=None, f_CA=1., f_biomass=1., adv_vel_U=0., adv_vel_V=0., dt=1., 
    custom_src_params=None, FV_grid=cyl_FVGrid, tol=1e-6, init_conc=float('NaN'), 
    steady_state=False, duration=900., chamber=True, chamber_rad=ch_params['radius'], chamber_vol=ch_params['volume'], 
    chamber_area=ch_params['area'], flow_rate=0., save_headspace_conc=True, save_profile=False, 
    gas_only=True):
    
    # initialize the microbial source-sink parameters
    if custom_src_params == None:
        custom_src_params = copy.copy(src_params)

    temp = np.array(temp, dtype='d')
    if air_temp == None:
        air_temp = temp[0,0]
    
    total_poros = np.array(total_poros, dtype='d')
    water_cont = np.array(water_cont, dtype='d')
    air_poros = total_poros - water_cont

    if FV_grid.z_grid.shape != temp.shape or FV_grid.z_grid.shape != total_poros.shape or (
        FV_grid.z_grid.shape != water_cont.shape):
        raise ValueError('Input data dimension does not match that of the grid.')

    n_within = np.argmin(np.abs(FV_grid.rad_node - chamber_rad))
    # the number of grid point columns within the chamber
    bdry_conc = FV_grid.rad_node * 0. + atm_conc  # top boundary concentration
    if ch_conc_init == None:
        ch_conc_init = atm_conc
    bdry_conc[0:n_within+1] = ch_conc_init
    n_steps = np.int(duration/dt)
    ch_conc = np.zeros(n_steps)
    ch_conc[0] = ch_conc_init
    ch_flux = np.zeros(n_steps)
    ch_flux_dfs = np.zeros(n_steps)
    ch_flux_adv = np.zeros(n_steps)

    D_air = diffus_in_air(species, air_temp)
    D_soil = diffus_in_soil_air(species, temp, total_poros, air_poros)
    bunsen_solub = solub_co2(temp)

    # by default, use 'gas only' if considering transient changes without dissolution equilibrium
    if gas_only: 
        bunsen_solub = bunsen_solub * 0. 

    # soil_conc = atm_conc * np.ones((FV_grid.n_vert_level, FV_grid.n_rad_level))
    # get the steady-state concentration profile (before the chamber is placed upon the soil column)
    temp_interp = np.interp(z_grid.grid_node, FV_grid.vert_node, temp[:,0])
    poros_interp = np.interp(z_grid.grid_node, FV_grid.vert_node, total_poros[:,0])
    swc_interp = np.interp(z_grid.grid_node, FV_grid.vert_node, water_cont[:,0])
    source_profile_interp = np.interp(z_grid.grid_node, FV_grid.vert_node, source_profile)
    soil_conc = flux_solver('co2', temp_interp, poros_interp, swc_interp, atm_conc=atm_conc, 
        source_profile=source_profile_interp, FV_grid=z_grid, dt=60, return_conc_only=True,
        custom_src_params=custom_src_params)
    soil_conc = np.interp(FV_grid.vert_node, z_grid.grid_node, soil_conc).reshape(25,1)
    soil_conc[0] = atm_conc
    soil_conc = np.tile(soil_conc, (1, len(FV_grid.rad_node)))

    eta = (air_poros + bunsen_solub * water_cont)

    beta = np.copy(D_soil)
    beta[0,:] = 2 * D_soil[0,:] * D_air / (D_soil[0,:] + D_air)
    # use the harmonic mean for soil-air boundary

    # source-sink terms
    if source_profile == None:
        gamma = soil_co2_source(temp, water_cont, FV_grid.z_grid,
            vmax=custom_src_params['soil_co2_source_vmax'], 
            decay_depth=custom_src_params['soil_resp_decay_depth'], 
            depth_func_type=custom_src_params['soil_resp_depth_func_type'], 
            swc_dependence=custom_src_params['soil_resp_swc_dependence'], 
            k_w=custom_src_params['soil_resp_water_factor'])
    else:
        gamma = np.tile(source_profile, (1, len(FV_grid.rad_node)))

    M = FV_grid.z_grid.shape[1] - 1
    N = FV_grid.z_grid.shape[0] - 1

    G_r = np.copy(beta)
    for i in range(0,M):
        G_r[:, i] = (beta[:, i] + beta[:, i+1]) * 0.5 / (FV_grid.rad_node[i+1] - FV_grid.rad_node[i])

    G_r[:, M] = 0.

    G_z = np.copy(beta)
    for j in range(0,N):
        G_z[j, :] = (beta[j, :] + beta[j+1, :]) * 0.5 / (FV_grid.vert_node[j+1] - FV_grid.vert_node[j])

    G_z[N, :] = 0.

    U = np.copy(adv_vel_U)
    U[:, 0:-1] = (adv_vel_U[:, 0:-1] + adv_vel_U[:, 1:]) * 0.5
    U[:, -1] = 0.

    V = np.copy(adv_vel_V)
    V[0:-1, :] = (adv_vel_V[0:-1, :] + adv_vel_V[1:, :]) * 0.5
    V[-1, :] = 0.

    # construct the curvature factors C_L and C_R
    curv_left = np.zeros(len(FV_grid.rad_node))
    curv_right = np.zeros(len(FV_grid.rad_node))
    curv_left[1:] = 2 * FV_grid.rad_cv_face[0:-1] / (FV_grid.rad_cv_face[1:]**2 - FV_grid.rad_cv_face[0:-1]**2)
    curv_left[0] = np.nan
    curv_right[1:] = 2 * FV_grid.rad_cv_face[1:] / (FV_grid.rad_cv_face[1:]**2 - FV_grid.rad_cv_face[0:-1]**2)
    curv_right[0] = np.nan

    # construct GR and GZ, the radial and vertical diffusion operators
    GR_matrix = np.zeros(( (N+1)*(M+1), (N+1)*(M+1) ))
    GZ_matrix = np.zeros(( (N+1)*(M+1), (N+1)*(M+1) ))
    # top boundry
    # GR_matrix[np.arange(0,M+1), :] = 0.
    # GZ_matrix[np.arange(0,M+1), :] = 0.
    # bottom boundry
    # GR_matrix[(M+1)*N + np.arange(0,M+1), :] = 0.
    GZ_matrix[(M+1)*N + np.arange(0,M+1), (M+1)*(N-1) + np.arange(0,M+1)] = 1.
    # lateral boundaries
    GR_matrix[(M+1) * np.arange(0,N+1) + M, (M+1) * np.arange(0,N+1) + M-1] = 1.
    # GZ_matrix[(M+1) * np.arange(0,N+1) + M, :] = 0.
    # axial points
    for j in range(1,N):
        GR_matrix[(M+1)*j, (M+1)*j] = -G_r[j,0] * 2 / FV_grid.rad_cv_face[0]
        GR_matrix[(M+1)*j, (M+1)*j+1] = G_r[j,0] * 2 / FV_grid.rad_cv_face[0]
        GZ_matrix[(M+1)*j, (M+1)*(j-1)] = G_z[j-1,i] / FV_grid.vert_cv_size[j]
        GZ_matrix[(M+1)*j, (M+1)*j] = (-G_z[j-1,i] - G_z[j,i]) / FV_grid.vert_cv_size[j]
        GZ_matrix[(M+1)*j, (M+1)*(j+1)] = G_z[j,i] / FV_grid.vert_cv_size[j]

    # interior points
    for j in range(1,N):
        for i in range(1,M):
            GR_matrix[(M+1)*j+i, (M+1)*j+i-1] = G_r[j,i-1] * curv_left[i]
            GR_matrix[(M+1)*j+i, (M+1)*j+i] = - G_r[j,i-1] * curv_left[i] - G_r[j,i] * curv_right[i]
            GR_matrix[(M+1)*j+i, (M+1)*j+i+1] = G_r[j,i] * curv_right[i]
            
            GZ_matrix[(M+1)*j+i, (M+1)*(j-1)+i] = G_z[j-1,i] / FV_grid.vert_cv_size[j]
            GZ_matrix[(M+1)*j+i, (M+1)*j+i] = (-G_z[j-1,i] - G_z[j,i]) /  FV_grid.vert_cv_size[j]
            GZ_matrix[(M+1)*j+i, (M+1)*(j+1)+i] = G_z[j,i] / FV_grid.vert_cv_size[j]

    # construct U_matrix, the radial advection operator
    U_matrix = np.zeros(( (N+1)*(M+1), (N+1)*(M+1) ))
    V_matrix = np.zeros(( (N+1)*(M+1), (N+1)*(M+1) ))
    # boundary points
    U_matrix[np.arange(0,M+1), :] = 0.
    U_matrix[(M+1)*N + np.arange(0,M+1), :] = 0.
    U_matrix[(M+1) * np.arange(0,N+1) + M, :] = 0.
    # axial points
    for j in range(1,N):
        U_matrix[(M+1)*j, (M+1)*j] = U[j,0] / FV_grid.rad_cv_face[0]
        U_matrix[(M+1)*j, (M+1)*j+1] = U[j,0] / FV_grid.rad_cv_face[0]

    # interior points
    for j in range(1,N):
        for i in range(1,M):
            U_matrix[(M+1)*j+i, (M+1)*j+i-1] = 0.5 * curv_left[i] * U[j,i-1]
            U_matrix[(M+1)*j+i, (M+1)*j+i] = 0.5 * curv_left[i] * U[j,i-1] - 0.5 * curv_right[i] * U[j,i]
            U_matrix[(M+1)*j+i, (M+1)*j+i+1] = - 0.5 * curv_right[i] * U[j,i]
    
    # construct V_matrix, the vertical advection operator
    # boundary points
    V_matrix[np.arange(0,M+1), :] = 0.
    V_matrix[(M+1)*N + np.arange(0,M+1), :] = 0.
    V_matrix[(M+1) * np.arange(0,N+1) + M, :] = 0.
    # all other points
    for i in range(0,M):
        for j in range(1,N):
            V_matrix[(M+1)*j+i, (M+1)*(j-1)+i] = 0.5 * V[j-1, i] / FV_grid.vert_cv_size[j]
            V_matrix[(M+1)*j+i, (M+1)*j+i] = 0.5 * (V[j-1, i] - V[j, i]) / FV_grid.vert_cv_size[j]
            V_matrix[(M+1)*j+i, (M+1)*(j+1)+i] = - 0.5 * V[j, i] / FV_grid.vert_cv_size[j]
    
    # construct a coefficient matrix, dimension (M+1)*(N+1) x (M+1)*(N+1)
    T_matrix = GR_matrix + GZ_matrix + U_matrix + V_matrix
    # other matrices
    A_matrix = np.diag( eta.reshape((M+1)*(N+1),) )
    C_matrix = soil_conc.reshape((M+1)*(N+1),1)
    S_matrix = gamma.reshape((M+1)*(N+1),1)
    
    grid_cross_sec = np.zeros(M+1)
    grid_cross_sec[0] = np.pi * FV_grid.rad_cv_face[0]**2
    grid_cross_sec[1:] = np.pi * (FV_grid.rad_cv_face[1:]**2 - FV_grid.rad_cv_face[0:-1]**2)
    grid_Xsec_ch = grid_cross_sec[0:n_within+1]
    grid_Xsec_ch[-1] += chamber_area - np.sum(grid_Xsec_ch)
    
    surf_flux_dfs = 0.5 * (beta[0, :] + beta[1, :]) * (C_matrix[ M+1 + np.arange(0,M+1)].reshape(20,) - 
        C_matrix[np.arange(0,M+1)].reshape(20,)) / (FV_grid.vert_node[1] - FV_grid.vert_node[0]) 
    surf_flux_adv = - adv_vel_V[0, :] * C_matrix[np.arange(0,M+1)].reshape(20,)

    ch_flux_dfs[0] = np.sum(surf_flux_dfs[0:n_within+1] * grid_Xsec_ch)
    ch_flux_adv[0] = np.sum(surf_flux_adv[0:n_within+1] * grid_Xsec_ch)
    ch_flux[0] = ch_flux_dfs[0] + ch_flux_adv[0]
    # these chamber fluxes are area-integrated
    # divide them by the chamber area to get chamber average flux rates
    
    for counter in range(1, n_steps):
        # call the flux solver to update the concentration field
        C_matrix_update = np.dot( np.linalg.inv(2 * A_matrix - dt * T_matrix), 
            np.dot((2 * A_matrix + dt * T_matrix), C_matrix) + 2 * dt * S_matrix )
        soil_conc_update = C_matrix_update.reshape(N+1, M+1)
        # enforce the boundary conditions again 
        C_matrix_update[np.arange(M+1)] = bdry_conc.reshape(20,1)  # top
        C_matrix_update[(M+1)*np.arange(1,N) + M] = C_matrix_update[(M+1)*np.arange(1,N) + M-1] # side
        C_matrix_update[(M+1)*N + np.arange(M+1)] = C_matrix_update[(M+1)*(N-1) + np.arange(M+1)] # bottom
        # calculate the surface flux rates
        surf_flux_dfs = 0.5 * (beta[0, :] + beta[1, :]) * (C_matrix_update[ M+1 + np.arange(0,M+1)].reshape(20,) - 
            C_matrix_update[np.arange(0,M+1)].reshape(20,)) / (FV_grid.vert_node[1] - FV_grid.vert_node[0]) 
        surf_flux_adv = - adv_vel_V[0, :] * C_matrix_update[np.arange(0,M+1)].reshape(20,)
        ## print surf_flux.shape
        ## print C_matrix_update[np.arange(0,M+1)].shape
        # integrate the surface flux enclosed in the chamber area
        ch_flux_dfs[counter] = np.sum(surf_flux_dfs[0:n_within+1] * grid_Xsec_ch)
        ch_flux_adv[counter] = np.sum(surf_flux_adv[0:n_within+1] * grid_Xsec_ch)
        ch_flux[counter] = ch_flux_dfs[counter] + ch_flux_adv[counter]
        # update the chamber concentration
        ch_conc[counter] = ch_conc[counter-1] + (flow_rate * (atm_conc - ch_conc[counter-1]) + ch_flux[counter] ) * dt / chamber_vol
        bdry_conc[0:n_within+1] = ch_conc[counter]
        C_matrix = C_matrix_update
        # call the 2D flux solver again with changed boundary condition
        # loop this until the chamber period ends

    # continue here
    soil_conc_final = C_matrix.reshape(N+1, M+1)

    return ch_conc

'''
    

    beta = D_soil * np.ones(FV_grid.n_level)
    #beta[0] = (beta[0] + D_air) / 2
    beta[0] = 2 * beta[0] * D_air / (beta[0] + D_air) 
    # use the harmonic mean for soil-air boundary
    beta[1:] = (beta[1:] + beta[0:-1]) / 2
    beta[0] = beta[0] / FV_grid.grid_node[0]
    beta[1:] = beta[1:] / (FV_grid.grid_node[1:] - FV_grid.grid_node[0:-1])

    gamma = soil_co2_source(temp, water_cont, FV_grid.grid_node) * \
        FV_grid.cv_size
    gamma[0] = gamma[0] + beta[0] * atm_conc

    if steady_state:
        max_iter = 1000000
    else:
        max_iter = np.ceil(duration / dt)

    soil_conc_updated = FV_CN_solver_1D(soil_conc, alpha, beta, gamma, dt)
    print soil_conc_updated.shape

    for loop_num in range(max_iter):
        soil_conc_updated = FV_CN_solver_1D(soil_conc, alpha, beta, gamma, dt)
        if loop_num == int(100/dt):
            print 'test'
            print soil_conc_updated 
        if ( np.linalg.norm(soil_conc_updated[0:-1] - soil_conc[0:-1]) / 
            np.linalg.norm(soil_conc[0:-1]) < tol) & (loop_num > 2):
            soil_conc = soil_conc_updated
            print loop_num
            break
        soil_conc = soil_conc_updated

    print soil_conc

    # test mass closure
    # source_func = np.append( [0.], 
    #   soil_co2_source(temp, water_cont, FV_grid.grid_node) )
    # depth_func = np.append([0.], FV_grid.grid_node)
    # print np.trapz(source_func, depth_func)

    surf_flux = beta[0] * (soil_conc[0] - atm_conc)
    return surf_flux
    #class FluxData():
    #return flux
'''
