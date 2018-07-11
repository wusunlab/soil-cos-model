# -*- coding: utf-8 -*-
import scipy.constants.constants as sci_const
import copy

'''
Global constants
'''
molec_wt = dict( [('h2o', 18.016), ('co2', 44.01), ('ch4', 16.04), 
    ('co', 28.01), ('cos', 60.07), ('no', 30.01), ('n2o', 44.013), 
    ('nh3', 17.031), ('h2', 2.016), ('cs2', 76.14), ('dms', 62.13)] ) 
# the list of molecular weights can be extended as tracers are added
p_std = sci_const.atm   # standard pressure, in Pa
R_gas = sci_const.R     # Universal gas constant (J K^-1 mol^-1)
air_conc_std = p_std / R_gas / 273.15
# air molar concentration (mol m^-3) at 273.15 K and 101.325 kPa

# mean northern hemisphere concentration in mixing ratio
# for model initialization only
mean_conc = dict( [('co2', 4e-4), ('ch4', 1.9e-6), 
    ('co', 0.), ('cos', 5e-10), ('no', 0.), ('n2o', 3.26e-10), 
    ('nh3', 0.), ('h2', 0.), ('cs2', 0.), ('dms', 0.)] )

'''
Chamber dimension
-----------------
Modify this if you have a different chamber
'''
ch_params = dict([])
ch_params['volume'] = 4.0761e-3  # Li-Cor 8100-104 long-term chamber
ch_params['area'] = 0.03178
ch_params['radius'] = 0.1

'''
Litter default properties
-------------------------
'''
# not used
lit_params = dict([])
lit_params['porosity'] = 0.95

'''
Source and sink parameters for tracers
'''
src_params = dict([])
src_params['soil_co2_source_vmax'] = 1e-5
src_params['soil_resp_decay_depth'] = 0.2     # soil respiration decay depth
src_params['soil_resp_q10'] = 1.5     # Q10 for soil respiration (not used)
src_params['soil_resp_swc_dependence'] = False  # by default, no water dependence
src_params['soil_resp_water_factor'] = 13.60  # exponential water dependence factor for respiration determined 
# from the Stunt Ranch data when soil moisture is low (<0.2)
src_params['soil_resp_depth_func_type'] = 'exp'

src_params['litter_co2_source_vmax'] = 1e-6
src_params['litter_resp_q10'] = 2.0     # Q10 for soil respiration (not used)
src_params['litter_resp_E'] = 1913     # E_resp / R_gas for temperature dependence, Stunt Ranch data
src_params['litter_resp_swc_dependence'] = False     # by default, no water dependence
src_params['litter_resp_water_factor'] = 14.47     # Exponential dependence factor for respiration determined 
# from the Stunt Ranch litter incubation data when litter moisture is low (<0.25)

src_params['soil_cos_sink_vmax'] = 1e-1
src_params['soil_cos_sink_Km'] = 1.9       # Michaelis constant, mol m^-3
src_params['soil_cos_sink_Delta_G_cat'] = 8.41e4        # enzyme kinetics parameter, J mol^-1
src_params['soil_cos_sink_Delta_H_eq'] = 3.59e5         # enzyme kinetics parameter, J mol^-1
src_params['soil_cos_sink_T_eq'] = 288.15       # temperature parameter, K
src_params['soil_cos_sink_swc_opt'] = 0.20        # optimal soil water content for COS uptake
src_params['soil_cos_sink_swc_dependence'] = 'nonlinear'  # sink activity is dependent on soil water content 
# for parameterization in Ogee et al (2015) BG
# CA enzyme kinetics data from Protoschill-Krebs et al. (1996) AE, from pea leaf CA
src_params['soil_cos_sink_CA_conc'] = 1e-4
src_params['soil_cos_sink_CA_kcat'] = 93.  # Protoschill-Krebs et al. (1996). AE
src_params['soil_cos_sink_CA_Km'] = 3.9e-2  # Protoschill-Krebs et al. (1996). AE
src_params['soil_cos_sink_CA_pKCA'] = 7.2  # Rowlett et al. (2002). Arch. Biochem. Biophys.
src_params['soil_cos_sink_CA_Ea'] = 4e4  # activation energy
src_params['soil_cos_sink_CA_Ed'] = 2e5  # deactivation energy
src_params['soil_cos_sink_CA_Sd'] = 660  # deactivation entropy

src_params['litter_cos_sink_vmax'] = 1.68e-3      # from the Stunt Ranch incubation data
src_params['litter_cos_sink_Km'] = 1.9     # Michaelis constant, mol m^-3
src_params['litter_cos_sink_water_factor'] = 11.56063     # exponential dependence factor for litter cos sink 
# from the Stunt Ranch incubation data

src_params['soil_cos_source_vmax'] = 2e-11
src_params['soil_cos_source_q10'] = 1.9    # an average value from Maseyk et al. (2014). PNAS
src_params['soil_cos_source_Eh_ref'] = -100    # in milli-Volt, see Devai and DeLaune (1995). Org. Geochem.
src_params['soil_cos_source_Eh_div'] = 20    # in milli-Volt, the divisor in the logistic function, set arbitrarily

src_params['litter_cos_source_vmax'] = 1.33e-11   # from the Stunt Ranch incubation data
src_params['litter_cos_source_q10'] = 1.9     # set the same as the soil one

'''
Modified parameters for field sites
-----------------------------------
'''

# UCLA Stunt Ranch Reserve, Los Angeles County, CA
src_params_stunt_ranch = copy.copy(src_params)   

# La Selva Biological Station, Costa Rica
src_params_la_selva = copy.copy(src_params)
