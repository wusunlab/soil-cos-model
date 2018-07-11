# calculate thermal diffusivity for the Oklahoma SGP ARM site
# equations from Johansen (1975), Peters-Lidards et al. (1998)
calc_therm_cond <- function(poros, wfps_frac) {
  dry_density = 2700 * (1-poros) # kg m^3
  dry_therm_cond = (0.135 * dry_density + 64.7) / (2700 - 0.947 * dry_density) # W m^-1 K^-1
  qtz_frac = 0.25 # 0.25 for silt loam soil
  qtz_therm_cond = 7.7 # quartz thermal conductivity, W m^-1 K^-1
  other_therm_cond = 2.0 # other mineral thermal conductivity, W m^-1 K^-1
  water_therm_cond = 0.57
  solid_therm_cond = qtz_therm_cond^qtz_frac * other_therm_cond^(1-qtz_frac)
  sat_therm_cond = solid_therm_cond^(1-poros) * water_therm_cond^poros
  
  # if (wfps_frac > 0.1) kstn_no = log10(wfps_frac) + 1. else kstn_no <- 0.
  flag <- (wfps_frac > 0.1)
  kstn_no <- flag * (log10(wfps_frac) + 1.) # Kersten number
  
  soil_therm_cond = kstn_no * (sat_therm_cond - dry_therm_cond) + dry_therm_cond
  return(soil_therm_cond)
}

# reading data
fileName <- './input/soil_flux_regular_meas.csv'
okFlux <- read.csv( fileName, skip = 0, header = TRUE)

porosity <- 0.45
vol_heat_cap <- 2e6 * (1-porosity) + 4.2e6 * okFlux$swc_avg
therm_cond <- calc_therm_cond(porosity, okFlux$swc_avg/porosity)
therm_diff <- therm_cond / vol_heat_cap
therm_decay_depth <- sqrt(2 * therm_diff / (2 * pi / 86400))
avg_therm_decay_depth <- mean(therm_decay_depth, na.rm=TRUE)
