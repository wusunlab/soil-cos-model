# -*- coding: utf-8 -*-
"""
Revision History
----------------
Wu Sun, 09/14/2014

Shape parameter dependency on soil texture is added 
Wu Sun, 05/21/2015

References updated
Wu Sun, 01/31/2016

References
----------
Gaseous diffusivities: H2O, CO2, CH4, CO, SO2, O3, NH3, N2O, NO, NO2, N2, O2
    Massman, W. J. (1998). A review of the molecular diffusivities of H2O, CO2, 
    CH4, CO, O3, SO2, NH3, N2O, NO, and NO2 in air, O2 and N2 near STP. 
    Atmos. Environ., 32(6), 1111-1127. 

COS gaseous diffusivity
    Seibt, U. et al. (2010). A kinetic analysis of leaf uptake of COS and its 
    relation to transpiration, photosynthesis and carbon isotope fractionation.
    Biogeosci., 7, 333–341.

Tortuosity parameter as a function of soil texture 
    Clapp, R. B. and Hornberger, G. M. (1978). Empirical equations for some 
    soil hydraulic properties. Water Resources Res., 14(4), 601-604. 

He, Ne, Kr, Xe, Rn, H2, CH4, CO2 aqueous diffusivities:
    Jähne, B. et al. (1987). Measurement of the diffusion coefficients of 
    sparingly soluble gases in water. J. Geophys. Res., 92(C10), 10767-10776.

CO and NO aqueous diffusivities:
    Wise, D. L. and Houghton, G. (1968). Diffusion coefficients of neon, 
    krypton, xenon, carbon monoxide and nitric oxide in water at 10-6O C. 
    Chem. Eng. Sci., 23, 1211-1216. 

COS aqueous diffusivity
    Ulshöfer, V. S. et al. (1996). Photochemical production and air-sea 
    exchange of sulfide in the eastern Mediterranean Sea. Mar. Chem., 53, 
    25-39. 
"""

import numpy as np
P_STD = 1.01325e5         # standard atmospheric pressure, in Pa

diffus_in_air_stp = dict( [('h2o', 2.178e-5), ('co2', 1.381e-5), 
    ('ch4', 1.952e-5), ('co', 1.807e-5), ('so2', 1.089e-5), ('o3', 1.444e-5), 
    ('nh3', 1.978e-5), ('n2o', 1.436e-5), ('no', 1.802e-5), ('no2', 1.361e-5),
    ('n2', 1.788e-5), ('o2', 1.820e-5), ('cos', 1.381e-5/1.21)] )

diffus_in_water_pre_exp = dict([('he', 818e-9 ), ('ne', 1608e-9), ('kr', 6393e-9),
    ('xe', 9007e-9), ('rn', 15877e-9), ('h2', 3338e-9), ('ch4', 3047e-9), 
    ('co2', 5019e-9), ('cos', 10**-1.3246 * 1e-4), ('co', 0.407e-4),
    ('no', 39.8e-4)])
diffus_in_water_E_act = dict([('he', 11.70e3), ('ne', 14.84e3), ('kr', 20.20e3), 
    ('xe', 21.61e3), ('rn', 23.26e3), ('h2', 16.06e3), ('ch4', 18.36e3), 
    ('co2', 19.51e3), ('cos', np.log(10)*1010*8.3144621), 
    ('co', 5860*4.184), ('no', 8360*4.184)])

def diffus_in_air(gas_name, temp, pressure=P_STD, kelvin=False):
    """
    Calculate the diffusivity of a gas in air (see the list below). 

    Supported gases include H2O, CO2, CH4, CO, SO2, O3, NH3, N2O, NO, NO2, 
    N2, O2 [M98]_ and COS [S10]_. 

    Parameters
    ----------
    gas_name : string
        Chemical name of the gas to be calculated for. Must be in lower case. 
    temp : float or `numpy.ndarray`
        Temperature, in Celsius degree by default.
    pressure : float or `numpy.ndarray`, optional
        Ambient pressure in Pa. Default value is standard atmospheric pressure.
        Note this is not the partial pressure.
    kelvin : boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    D_air : float or `numpy.ndarray`
        Diffusivity (m^2 s^-1) in air. 

    Raises
    ------
    ValueError
        If `gas_name` is not a string.
        If `gas_name` is not a supported gas species.

    See also
    --------
    diffus_in_water : Calculate the diffusivity of a gas in water.
    diffus_in_soil_air : Calculate the diffusivity of a gas in soil air.
    diffus_in_soil_water : Calculate the diffusivity of a gas in soil water.

    References
    ----------
    .. [M98] Massman, W. J. (1998). A review of the molecular diffusivities of 
             H2O, CO2, CH4, CO, O3, SO2, NH3, N2O, NO, and NO2 in air, O2 and 
             N2 near STP. Atmos. Environ., 32(6), 1111-1127. 
    .. [S10] Seibt, U. et al. (2010). A kinetic analysis of leaf uptake of COS 
             and its relation to transpiration, photosynthesis and carbon 
             isotope fractionation. Biogeosci., 7, 333–341.
    """
    if not isinstance(gas_name, str):
        raise ValueError('Not a proper gas name. ')
    if not gas_name in diffus_in_air_stp:
        raise ValueError('This gas is not included in this function. ')
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    pressure = np.array(pressure, dtype='d')
    
    D_air = diffus_in_air_stp[gas_name] * (P_STD / pressure) * \
        (T_k / 273.15)**1.81

    return D_air

def diffus_in_soil_air(gas_name, temp, total_poros, air_poros, 
    pressure=P_STD, shape_param=4.9, soil_texture='', 
    kelvin=False):
    """
    Calculate the diffusivity of a gas in soil air (see the list below). 

    Supported gases include H2O, CO2, CH4, CO, SO2, O3, NH3, N2O, NO, NO2, 
    N2, O2 [M98]_ and COS [S10]_. 

    Parameters
    ----------
    gas_name : string
        Chemical name of the gas to be calculated for. Must be in lower case. 
    temp : float or `numpy.ndarray`
        Temperature, in Celsius degree by default.
    total_poros : float or `numpy.ndarray`
        Total porosity of soil, in m^3 m^-3.
    air_poros : float or `numpy.ndarray`
        Air-filled porosity of soil, in m^3 m^-3.  
    pressure : float or `numpy.ndarray`, optional
        Ambient pressure in Pa. Default value is standard atmospheric pressure.
        Note this is not the partial pressure.
    shape_param : float, optional
        A shape parameter for soil water retention curve that determines the 
        tortuosity of soil gas diffusion. It depends on soil texture [CH78]_. 
        The default value 4.9 is for sandy loam soil.
    soil_texture : string, optional
        Soil texture name, in lower case. When not given properly, sandy loam
        soil is assumed as the default. 
    kelvin : boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    D_soil_air : float or `numpy.ndarray`
        Diffusivity (m^2 s^-1) in soil air. 

    Raises
    ------
    ValueError
        If `gas_name` is not a string.
        If `gas_name` is not a supported gas species.
    
    See also
    --------
    diffus_in_air : Calculate the diffusivity of a gas in air.
    diffus_in_water : Calculate the diffusivity of a gas in water.
    diffus_in_soil_water : Calculate the diffusivity of a gas in soil water.

    References
    ----------
    .. [M98]  Massman, W. J. (1998). A review of the molecular diffusivities of 
              H2O, CO2, CH4, CO, O3, SO2, NH3, N2O, NO, and NO2 in air, O2 and 
              N2 near STP. Atmos. Environ., 32(6), 1111-1127. 
    .. [S10]  Seibt, U. et al. (2010). A kinetic analysis of leaf uptake of COS 
              and its relation to transpiration, photosynthesis and carbon 
              isotope fractionation. Biogeosci., 7, 333–341.
    .. [CH78] Clapp, R. B. and Hornberger, G. M. (1978). Empirical equations 
              for some soil hydraulic properties. Water Resources Res., 14(4), 
              601-604. 
    """
    if not isinstance(gas_name, str):
        raise ValueError('Not a proper gas name. ')
    if not gas_name in diffus_in_air_stp:
        raise ValueError('This gas is not included in this function. ')
    shape_param_list = dict( [('sand', 4.05), ('loamy sand', 4.38), 
        ('sandy loam', 4.9), ('silt loam', 5.30), ('loam', 5.39), 
        ('sandy clay loam', 7.12), ('silty clay loam', 7.75), 
        ('clay loam', 8.52), ('sandy clay', 10.4),
        ('silty clay', 10.4), ('clay', 11.4) ])
    # from Clapp and Hornberg (1978) Water Resources Res.
    if isinstance(soil_texture, str) and (soil_texture in shape_param_list):
        shape_param = shape_param_list[soil_texture]
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    pressure = np.array(pressure, dtype='d')
    total_poros = np.array(total_poros, dtype='d')
    air_poros = np.array(air_poros, dtype='d')
    D_soil_air = diffus_in_air_stp[gas_name] * (P_STD / pressure) * \
        (T_k / 273.15)**1.81 * air_poros**2 * (air_poros/total_poros) ** \
        (3./shape_param)
    return D_soil_air

def diffus_in_water(gas_name, temp, pressure=P_STD, kelvin=False):
    u"""
    Calculate the diffusivity of a gas in water (see the list below). 

    Supported gases include He, Ne, Kr, Xe, Rn, H2, CH4, CO2 [J87]_, 
    CO, NO [WH68]_, and COS [U96]_. 

    Parameters
    ----------
    gas_name : string
        Chemical name of the gas to be calculated for. Must be in lower case. 
    temp : float or `numpy.ndarray`
        Temperature, in Celsius degree by default.
    pressure : float or `numpy.ndarray`, optional
        Ambient pressure in Pa. Default value is standard atmospheric pressure.
        Note this is not the partial pressure.
    kelvin : boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    D_aq : float or `numpy.ndarray`
        Diffusivity (m^2 s^-1) in water. 

    Raises
    ------
    ValueError
        If `gas_name` is not a string.
        If `gas_name` is not a supported gas species.
    
    See also
    --------
    diffus_in_air : Calculate the diffusivity of a gas in air.
    diffus_in_soil_air : Calculate the diffusivity of a gas in soil air.
    diffus_in_soil_water : Calculate the diffusivity of a gas in soil water.

    References
    ----------
    .. [J87]  Jähne, B. et al. (1987). Measurement of the diffusion 
              coefficients of sparingly soluble gases in water. 
              J. Geophys. Res., 92(C10), 10767-10776.
    .. [WH68] Wise, D. L. and Houghton, G. (1968). Diffusion coefficients of 
              neon, krypton, xenon, carbon monoxide and nitric oxide in water 
              at 10-6O C. Chem. Eng. Sci., 23, 1211-1216. 
    .. [U96]  Ulshöfer, V. S. et al. (1996). Photochemical production and 
              air-sea exchange of sulfide in the eastern Mediterranean Sea. 
              Mar. Chem., 53, 25-39. 
    """
    if not isinstance(gas_name, str):
        raise ValueError('Not a proper gas name. ')
    if not gas_name in diffus_in_water_E_act:
        raise ValueError('This gas is not included in this function. ')
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    # pressure = np.array(pressure, dtype='d')
    
    D_aq = diffus_in_water_pre_exp[gas_name] * \
      np.exp( - diffus_in_water_E_act[gas_name] / 8.3144621 / T_k )

    return D_aq

def diffus_in_soil_water(gas_name, temp, total_poros, air_poros, 
    pressure=P_STD, shape_param=4.9, soil_texture='', 
    kelvin=False):
    u"""
    Calculate the diffusivity of a gas in soil water (see the list below). 

    Supported gases include He, Ne, Kr, Xe, Rn, H2, CH4, CO2 [J87]_, 
    CO, NO [WH68]_, and COS [U96]_. 

    Parameters
    ----------
    gas_name : string
        Chemical name of the gas to be calculated for. Must be in lower case. 
    temp : float or `numpy.ndarray`
        Temperature, in Celsius degree by default.
    total_poros : float or `numpy.ndarray`
        Total porosity of soil, in m^3 m^-3.
    air_poros : float or `numpy.ndarray`
        Air-filled porosity of soil, in m^3 m^-3.  
    pressure : float or `numpy.ndarray`, optional
        Ambient pressure in Pa. Default value is standard atmospheric pressure.
        Note this is not the partial pressure.
    shape_param : float, optional
        A shape parameter for soil water retention curve that determines the 
        tortuosity of soil gas diffusion. It depends on soil texture [CH78]_. 
        The default value 4.9 is for sandy loam soil.
    soil_texture : string, optional
        Soil texture name, in lower case. When not given properly, sandy loam
        soil is assumed as the default. 
    kelvin : boolean, optional
        Temperature input is in Kelvin if enabled.

    Returns
    -------
    D_soil_aq : float or `numpy.ndarray`
        Diffusivity (m^2 s^-1) in soil water. 

    Raises
    ------
    ValueError
        If `gas_name` is not a string.
        If `gas_name` is not a supported gas species.

    See also
    --------
    diffus_in_air : Calculate the diffusivity of a gas in air.
    diffus_in_water : Calculate the diffusivity of a gas in water.
    diffus_in_soil_air : Calculate the diffusivity of a gas in soil air.

    References
    ----------
    .. [J87]  Jähne, B. et al. (1987). Measurement of the diffusion 
              coefficients of sparingly soluble gases in water. 
              J. Geophys. Res., 92(C10), 10767-10776.
    .. [WH68] Wise, D. L. and Houghton, G. (1968). Diffusion coefficients of 
              neon, krypton, xenon, carbon monoxide and nitric oxide in water 
              at 10-6O C. Chem. Eng. Sci., 23, 1211-1216. 
    .. [U96]  Ulshöfer, V. S. et al. (1996). Photochemical production and 
              air-sea exchange of sulfide in the eastern Mediterranean Sea. 
              Mar. Chem., 53, 25-39. 
    .. [CH78] Clapp, R. B. and Hornberger, G. M. (1978). Empirical equations 
              for some soil hydraulic properties. Water Resources Res., 14(4), 
              601-604. 
    """
    if not isinstance(gas_name, str):
        raise ValueError('Not a proper gas name. ')
    if not gas_name in diffus_in_water_E_act:
        raise ValueError('This gas is not included in this function. ')
    shape_param_list = dict( [('sand', 4.05), ('loamy sand', 4.38), 
        ('sandy loam', 4.9), ('silt loam', 5.30), ('loam', 5.39), 
        ('sandy clay loam', 7.12), ('silty clay loam', 7.75), 
        ('clay loam', 8.52), ('sandy clay', 10.4),
        ('silty clay', 10.4), ('clay', 11.4) ])
    # from Clapp and Hornberg (1978) Water Resources Res.
    if isinstance(soil_texture, str) and (soil_texture in shape_param_list):
        shape_param = shape_param_list[soil_texture]
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15
    # force temperature to be in Kelvin
    # pressure = np.array(pressure, dtype='d')
    total_poros = np.array(total_poros, dtype='d')
    water_poros = total_poros - np.array(air_poros, dtype='d')
    D_soil_aq = diffus_in_water_pre_exp[gas_name] * \
      np.exp( - diffus_in_water_E_act[gas_name] / 8.3144621 / T_k ) \
      * water_poros**2 * (water_poros/total_poros) ** (shape_param/3.-1.)
    
    return D_soil_aq
