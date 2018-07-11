# -*- coding: utf-8 -*-
''' 
Revision history
----------------
Wu Sun, 09/29/2015

References
----------
Wagner, W. and Pruß, A. (2002). J. Phys. Chem. Ref. Data, 31, 387.

Millero, F. J. and Poisson, A. (1981). Deep Sea Res., 28A(6), 625-629.
''' 

import numpy as np
def water_density(temp, kelvin=False):
    u"""
    Calculate water density (kg m^-3) as a function of temperature.

    Parameters
    ----------
    temp : float or `numpy.ndarray`
        Temperature, in Celsius degree by default.
    kelvin : boolean, optional
        Temperature input is in Kelvin if enabled.
    
    Returns
    -------
    rho_w : float
        Water density (kg m^-3)

    Raises
    ------
    None
    
    References
    ----------
    .. [WP02] Wagner, W. and Pruß, A. (2002). The IAPWS Formulation 1995 for 
              the Thermodynamic Properties of Ordinary Water Substance for 
              General and Scientific Use. J. Phys. Chem. Ref. Data, 31, 387.
    """
    
    T_k = np.array(temp, dtype='d') + (not kelvin) * 273.15     
    # force temperature to be in Kelvin
    
    # critical temperature and density
    T_crit = 647.096
    rho_crit = 322.
    
    # parameters
    theta = 1 - T_k / T_crit
    b = np.array([1.99274064, 1.09965342, -0.510839303, -1.7549349, -45.5170352, -6.74694450e5])
    exponents = np.array([1., 2., 5., 16., 43., 110.]) / 3.
    
    rho_ratio = np.ones_like(theta)  # initialize with 1
    for m in range(b.size):
        rho_ratio += b[m] * theta**exponents[m]
    
    rho_w = rho_crit * rho_ratio
    return rho_w

def seawater_density(temp, salinity, kelvin=False):
    """
    Calculate seawater density (kg m^-3) as a function of temperature and salinity.
	
	Applicable range: 1 atm, 0-40 C, 0.5-43 g kg^-1 salinity.

    Parameters
    ----------
    temp : float or `numpy.ndarray`
        Temperature, in Celsius degree by default.
    salinity : float or `numpy.ndarray`
        Salinity, in per mil mass fraction (g kg^-1 seawater).
    kelvin : boolean, optional
        Temperature input is in Kelvin if enabled.
    
    Returns
    -------
    rho_sw : float or `numpy.ndarray`
        Sea water density (kg m^-3)

    Raises
    ------
    ValueError        
        If salinity is negative or exceeds saturation.
    
    References
    ----------
    .. [MP81] Millero, F. J. and Poisson, A. (1981). International 
              one-atmosphere equation of state of seawater. 
              Deep Sea Res., 28A(6), 625-629.
    """
    if (np.sum(salinity < 0.)):
        raise ValueError('Salinity error, cannot use negative value. ')
    if (np.sum(salinity > 360.)):
        # The saturation concentration of NaCl in water is 360 g/L
        raise ValueError('Salinity error, oversaturation reached. ')

    T_c = np.array(temp, dtype='d') - kelvin * 273.15
    # force temperature to be in Celsius
    salinity = np.array(salinity, dtype='d')
    
    # parameters
    rho_0_coefs = np.array([999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6, 6.536336e-9, ])
    A_coefs = np.array([8.24493e-1, -4.0899e-3, 7.6438e-5, 8.2467e-7, 5.3875e-9, ])
    B_coefs = np.array([-5.72466e-3, 1.0227e-4, -1.6546e-6, ])
    
    # calculate each term
    rho_0 = rho_0_coefs[0] + rho_0_coefs[1] * T_c + rho_0_coefs[2] * T_c**2 + rho_0_coefs[3] * T_c**3 + rho_0_coefs[4] * T_c**4 + rho_0_coefs[5] * T_c**5
    A_term = A_coefs[0] + A_coefs[1] * T_c + A_coefs[2] * T_c**2 + A_coefs[3] * T_c**3 + A_coefs[4] * T_c**4
    B_term = B_coefs[0] + B_coefs[1] * T_c + B_coefs[2] * T_c**2 
    C = 4.8314e-4

    rho_sw = rho_0 + A_term * salinity + B_term * salinity**1.5 + C * salinity
    return rho_sw
