# -*- coding: utf-8 -*-
'''
Calculate sources and sinks of tracers.
Tracers include: CO2, COS

Revision history
----------------
Wu Sun, 09/18/2014
Wu Sun, 07/30/2015
Wu Sun, 09/30/2015

References
----------
CO2:
Lloyd and Taylor. (1994) Func. Ecol.

COS:
Ogawa et al. (2013). J. Am. Chem. Soc.
Sun et al. (2015). Geosci. Model Dev.
Ogee et al. (2015). Biogeosci. Discuss.
'''

from soil_tracers_config import *
# loading constants and parameters
import numpy as np


def soil_co2_source(temp, water_cont=0., depth=0., vmax=src_params['soil_co2_source_vmax'], 
    decay_depth=src_params['soil_resp_decay_depth'], depth_func_type=src_params['soil_resp_depth_func_type'], 
    swc_dependence=src_params['soil_resp_swc_dependence'], k_w=src_params['soil_resp_water_factor'], 
    kelvin=False):
    """
    Calculate CO2 source in soil.

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    water_cont: float
        Soil volumetric water content. Not necessary if the swc_dependence is set false. 
    depth: float
        Depth of the grid, used to scale respiration rate.
    vmax: float, optional
        CO2 source activity, in mol m^-3 s^-1. By default it is set in the config file. 
    decay_depth: float, optional
        The e-fold decay depth of soil respiration. By default set in the config file.
    swc_dependence: boolean, optional
        If enabled, CO2 production depends exponentially on soil water content.
    k_w: float, optional
        An exponential factor for soil respiration dependence on water content.
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    co2_source_rate: float
        CO2 source rate, in mol m^-3 s^-1.

    Raises
    ------
    ValueError
        If water_cont is negative.

    """
    water_cont = np.array(water_cont, dtype='d')
    if (np.sum(water_cont < 0)):
        raise ValueError('Soil water content cannot be negative. ')
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    depth = np.array(depth, dtype='d')
    
    A = 77.0683
    # pre-exponential factor to normalize the temperature dependence function

    # k_w = 13.60       # water dependence factor
    co2_source_rate = vmax * A * np.exp(-308.56 / (T_k - 227.13) ) * np.exp(- depth / decay_depth )  
    if swc_dependence:
        co2_source_rate *= np.sinh(k_w * water_cont)

    return co2_source_rate

def litter_co2_source(temp, water_cont=0., vmax=src_params['litter_co2_source_vmax'], 
    E_resp=src_params['litter_resp_E'], k_w=src_params['litter_resp_water_factor'], 
    lwc_dependence=src_params['litter_resp_swc_dependence'], kelvin=False):
    """
    Calculate CO2 source in litter.
    Note the current parameterization is derived from Sun et al. (2016) J. Geophys. Res. Biogeosci.

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    water_cont: float
        Litter water content, in g water per g litter.
    vmax: float, optional
        CO2 source activity, in mol m^-3 s^-1. By default it is set in the config file. 
    E_resp: float, optional
        A factor in the Arrhenius equation for temperature dependence
    k_w: float, optional
        An exponential factor for litter respiration dependence on water content.
    lwc_dependence: boolean, optional
        If enabled, CO2 production depends exponentially on litter water content.
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    co2_source_rate: float
        CO2 source rate, in mol m^-3 s^-1.

    Raises
    ------
    ValueError
        If water_cont is negative.

    """
    water_cont = np.array(water_cont, dtype='d')
    if (np.sum(water_cont < 0)):
        raise ValueError('Litter water content cannot be negative. ')
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    
    co2_source_rate = vmax * np.exp(E_resp * (1./298.15 - 1./T_k))
    if lwc_dependence:
        co2_source_rate *= np.sinh(k_w * water_cont)

    return co2_source_rate

def uncat_cos_hydrol(temp, pH, kelvin=False):
    '''
    Uncatalyzed COS hydrolysis rate constant in water as a function of temeprature and pH.
    
    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    k_uncat: float
        COS hydrolysis rate constant, in s^-1.
        Multiply it by the aqueous concentration of COS to get the hydrolysis rate

    Raises
    ------
    None
    '''
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    
    k_uncat = np.exp(24.3 - 10450./T_k) + 10**(-14.+pH) * np.exp(22.8 - 6040./T_k)
    return k_uncat

def soil_cos_sink(temp, water_cont, cos_conc, vmax=src_params['soil_cos_sink_vmax'], 
    K_m=src_params['soil_cos_sink_Km'], Delta_G_cat=src_params['soil_cos_sink_Delta_G_cat'],
    Delta_H_eq=src_params['soil_cos_sink_Delta_H_eq'], T_eq=src_params['soil_cos_sink_T_eq'],
    swc_opt=src_params['soil_cos_sink_swc_opt'], kelvin=False, 
    swc_dependence=src_params['soil_cos_sink_swc_dependence']):
    """
    Calculate COS sink in soil.

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    water_cont: float
        Soil volumetric water content.
    cos_conc: float, optional
        COS aqueous concentration, in mol m^-3.
    vmax: float, optional
        COS sink activity, in mol m^-3 s^-1. By default it is set in the config file. 
    K_m: float, optional
        Michaelis constant for COS hydrolysis via CA enzyme, mol m^-3.
    Delta_G_cat: float, optional
    Delta_H_eq: float, optional
        Enzyme kinetics parameters for the temperature dependence function, J mol^-1
    T_eq: float, optional
        Temperature parameter for the temperature dependence function, K
    swc_opt: float, optional
        Optimal SWC for the SWC dependence function
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.
    swc_dependence: string, optional
        If 'nonlinear', the soil water content dependence has an optimum (default).
        If 'linear', SWC dependence is linear (multiplied by water content).
        If 'none', SWC dependence is off.
    
    Returns
    -------
    cos_sink_rate: float
        COS sink rate, in mol m^-3 s^-1.

    Raises
    ------
    ValueError
        If water_cont is negative.
        If cos_conc is negative.
    """

    water_cont = np.array(water_cont, dtype='d')
    cos_conc = np.array(cos_conc, dtype='d')
    if (np.sum(water_cont < 0)):
        raise ValueError('Soil water content cannot be negative. ')
    if (np.sum(cos_conc < 0)):
        raise ValueError('COS concentration cannot be negative. ')
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    
    # temperature limitation function
    T_opt = T_eq * (1 - R_gas / Delta_H_eq * np.log( (Delta_H_eq - Delta_G_cat) / Delta_G_cat ) * T_eq )
    norm_fac = T_opt * np.exp(-Delta_G_cat / (R_gas*T_opt) ) / (1 + np.exp( -Delta_H_eq / R_gas * (1/T_opt - 1/T_eq) ) ) 
    norm_fac = 1. / norm_fac
    T_limit = norm_fac * T_k * np.exp(-Delta_G_cat / (R_gas*T_k) ) / \
        (1 + np.exp( -Delta_H_eq / R_gas * (1/T_k - 1/T_eq) ) )
    
    # water limitation function
    if swc_dependence == 'nonlinear':
        swc_limit = water_cont / swc_opt**2 * np.exp(- (water_cont/swc_opt)**2/2 ) / ( np.exp(-0.5)/swc_opt )
    elif swc_dependence == 'linear':
        swc_limit = water_cont
    elif swc_dependece == 'none':
        swc_limit = 1.
    
    cos_sink_rate = - vmax * cos_conc / (K_m + cos_conc) * T_limit * swc_limit

    return cos_sink_rate

def soil_cos_sink_Ogee2015(temp, pH, cos_conc, CA_conc=src_params['soil_cos_sink_CA_conc'], 
    K_m=src_params['soil_cos_sink_CA_Km'], k_cat=src_params['soil_cos_sink_CA_kcat'],
    E_a=src_params['soil_cos_sink_CA_Ea'], E_d=src_params['soil_cos_sink_CA_Ed'], 
    S_d=src_params['soil_cos_sink_CA_Sd'], pK_CA=src_params['soil_cos_sink_CA_pKCA'],
    kelvin=False, return_f_CA=False):
    '''
    Calculate COS sink in soil. Parameterization from Ogee et al. (2015). BG.
    Note: multiply the returned value by the soil water content to get the sink rate 
    with respect to the total soil volume (otherwise, it is the uptake rate
    per unit volume of the soil water).  

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    pH: float
        Soil pH.
    cos_conc: float, optional
        COS aqueous concentration, in mol m^-3.
    CA_conc: float, optional
        CA concentration in soil solution, in mol m^-3 s^-1. 
        By default it is set in the config file. 
    K_m: float, optional
        Michaelis constant for COS hydrolysis via CA enzyme, mol m^-3.
    E_a: float, optional
        Enzyme kinetics parameters for the temperature dependence function, J mol^-1
    E_d: float, optional
        Enzyme kinetics parameters for the temperature dependence function, J mol^-1
    S_d: float, optional
        Enzyme kinetics parameters for the temperature dependence function, J mol^-1 K^-1
    pK_CA: float, optional
        Dissociation constant for CA
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.
    return_f_CA: boolean, optional
        If set True, return the f_CA value in addition to cos_sink_rate.

    Returns
    -------
    cos_sink_rate: float
        COS sink rate, in mol m^-3 s^-1.

    Raises
    ------
    None
    
    Examples
    --------
    >>> soil_cos_sink_Ogee2015(25, 7, 2e-8)
    -4.5261580724797873e-09
    
    >>> soil_cos_sink_Ogee2015(25, 7, 2e-8, return_f_CA=True)
    (-4.5261580724797873e-09, 10545.196854938044)
    
    '''
    
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    
    pH_in = 6 + 0.25*pH  # calculate intracellular pH from soil pH
    # expression from Ogee et al. (2015) BG
    
    def __temp_func(T_k, E_a=E_a, E_d=E_d, S_d=S_d):
        return np.exp(-E_a/R_gas/T_k) / (1. + np.exp(-E_d/R_gas/T_k + S_d/R_gas))
    
    def __pH_func(pH, pK_CA=pK_CA):
        return 1./(1.+10**(-pH + pK_CA))
    
    T_limit = __temp_func(T_k) / __temp_func(293.15)
    pH_limit = __pH_func(pH_in) / __pH_func(8.2)
    
    cos_sink_rate = - k_cat * CA_conc * cos_conc / (K_m + cos_conc) * T_limit * pH_limit
    
    # print CA_conc, K_m, k_cat, E_a, E_d, S_d, pK_CA # for test
    
    f_CA = - cos_sink_rate / uncat_cos_hydrol(25, 4.5) / cos_conc
    
    if not return_f_CA:
        return cos_sink_rate
    else:
        return cos_sink_rate, f_CA

def litter_cos_sink(temp, water_cont, cos_conc, vmax=src_params['litter_cos_sink_vmax'],
    K_m=src_params['litter_cos_sink_Km'], k_w=src_params['litter_cos_sink_water_factor'], kelvin=False):
    """
    Calculate COS sink in litter.

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    water_cont: float
        Litter water content, in g water per g litter.
    cos_conc: float, optional
        COS aqueous concentration, in mol m^-3.
    vmax: float, optional
        COS sink activity, in mol m^-3 s^-1. By default it is set in the config file. 
    K_m: float, optional
        Michaelis constant for COS hydrolysis via CA enzyme, mol m^-3.
    k_w: float, optional
        An exponential factor for litter COS sink dependence on water content.
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    cos_sink_rate: float
        COS sink rate, in mol m^-3 s^-1.

    Raises
    ------
    ValueError
        If water_cont is negative.
        If cos_conc is negative.
    """
    
    water_cont = np.array(water_cont, dtype='d')
    cos_conc = np.array(cos_conc, dtype='d')
    if (np.sum(water_cont < 0)):
        raise ValueError('Litter water content cannot be negative. ')
    if (np.sum(cos_conc < 0)):
        raise ValueError('COS concentration cannot be negative. ')
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    
    T_limit = 1.
    # no temperature dependence

    cos_sink_rate = - vmax * cos_conc / (K_m + cos_conc) * T_limit * np.sinh(k_w * water_cont)

    return cos_sink_rate

def soil_cos_source(temp, vmax=src_params['soil_cos_source_vmax'], q10=src_params['soil_cos_source_q10'], 
    kelvin=False, redox_pot_response=False, E_h=None, Eh_ref=src_params['soil_cos_source_Eh_ref'], 
    Eh_div=src_params['soil_cos_source_Eh_div']):
    """
    Calculate COS source in soil.

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    vmax: float, optional
        COS source activity, in mol m^-3 s^-1. By default it is set in the config file. 
    q10: float, optional
        Q10 value for COS source temperature dependence
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.
    redox_pot_response: boolean, optional
        Activate a redox potential dependence term for COS production, see Ogeé et al. (2015). Biogeosci.
    E_h: float, optional
        Redox potential in mV. Must pass this if redox_pot_response is set True
    Eh_ref: float, optional
        Reference redox potential in mV at which soil COS production shoots up
    Eh_div: float, optional
        A divisor in the logistic function of E_h dependence, set arbitrarily. See Ogee et al. (2015). Biogeosci.

    Returns
    -------
    cos_source_rate: float
        COS source rate, in mol m^-3 s^-1.

    Raises
    ------
    ValueError
        - if redox_pot_response==True and E_h==None
    """
    if (redox_pot_response and E_h==None):
        raise ValueError('Must pass the value of E_h if redox potential dependence is enabled. ')

    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin

    T_ref = 298.15
    T_func = np.exp( np.log(q10)/10 * (T_k - T_ref) )  # temperature function
    
    # no moisture dependence
    
    # E_h dependence
    if redox_pot_response:
        Eh_func = 1. / (1. + np.exp(E_h - Eh_ref)/Eh_div)
    else:
        Eh_func = 1.

    cos_source_rate = vmax * T_func * Eh_func

    return cos_source_rate

def litter_cos_source(temp, vmax=src_params['litter_cos_source_vmax'], q10=src_params['litter_cos_source_q10'], 
    kelvin=False):
    """
    Calculate COS source in litter.

    Parameters
    ----------
    temp: float
        Temperature, in Celsius degree by default.
    vmax: float, optional
        COS source activity, in mol m^-3 s^-1. By default it is set in the config file. 
    q10: float, optional
        Q10 value for COS source temperature dependence
    kelvin: boolean, optional
        Temperature input is in Kelvin if enabled.
    
    Returns
    -------
    cos_source_rate: float
        COS source rate, in mol m^-3 s^-1.

    Raises
    ------
    None
    """
    
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    
    T_ref = 298.15
    T_func = np.exp( np.log(q10)/10 * (T_k - T_ref) )
    # temperature function

    # no moisture dependence

    cos_source_rate = vmax * T_func

    return cos_source_rate
