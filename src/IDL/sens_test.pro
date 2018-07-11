; A fast version of the soil flux model
; tolerance 1e-5
; timestep 10 s
; number of layers 25
; for Oklahoma wheat field soil flux
; no litter layers

;function lit_evap, T, rh, lit_wc
;  p = [-0.007658066, 0.001516583, 7.641850059] 
;  return, (p[0] + p[1] * T) * sinh(p[2] * lit_wc)
;end

function calc_cos_solub, T, Kelvin=kelvin
  compile_opt idl2
  if n_params() ne 1 then message,'ERROR: Number of parameters incorrect.'
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  ;if keyword_set(kelvin) then T_c = T - 273.15 else T_c = T
  ;k_H = 0.0008481 * T_c^2 - 0.055736 * T_c + 1.3224 ; k_H is a function of Celsius temp
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15   ; Data from Elliott 1989 Env Sci Tech
  k_H = T_k * exp(4050.32 / T_k - 20.0007) ; k_H is a function of Kelvin temp
  return, k_H       ; dimensionless
end

function calc_co2_solub, T, Kelvin=kelvin
  compile_opt idl2
  if n_params() ne 1 then message,'ERROR: Number of parameters incorrect.'
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  A_1 = - 58.0931
  A_2 = 90.5069
  A_3 = 22.2940
  K_0 = exp( A_1 + A_2 * (100. / T_k) + A_3 * alog(T_k / 100.))   ; in mol L^-1 atm (Weiss 1974 Mar Chem)
  p_std = 1.01325e5           ; Standard atmosphere pressure (N m^-2)
  R_gas = 8.3144621           ; universal gas constant
  air_conc = p_std / R_gas / T_k    ; air molar concentration at temperature T_k
  k_H = 1000 * K_0 / air_conc      ; non-dimensional
  return, k_H       ; dimensionless
end

function calc_co2_diff, T, theta_a, theta_sat, Air=air, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 1 then message,'ERROR: Number of parameters incorrect.'
  mol_diff = 1.618e-5         ; molecular diffusivity of CO2 (m s^-2) ; see Massman 1998 Atmos Env
  T_ref = 298.15            ; reference temperature 298.15 K
  ; Moldrup et al. (2003). Soil Sci.
  n = 1.5                   ; empirical constant 1
  b = 4.9                   ; empirical constant 2, 4.9 for sandy loam, 5.3 for silt loam
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  if keyword_set(air) then D_co2 = mol_diff * (T_k/T_ref)^0.5 $ 
  else begin
    if n_params() ne 3 then message,'ERROR: Number of parameters incorrect.'
    D_co2 = mol_diff * (T_k/T_ref)^n * theta_a^2 * (theta_a / theta_sat)^(3/b)
  endelse
  return, D_co2
end

function calc_cos_diff, T, theta_a, theta_sat, Air=air
  compile_opt idl2
  if n_params() lt 1 then message,'ERROR: Number of parameters incorrect.'
  mol_diff = 1.618e-5/1.21    ; molecular diffusivity of COS (m s^-2)   ; see Seibt 2010 Biogeosci
  T_ref = 298.15            ; reference temperature 298.15 K
  ; Moldrup et al. (2003). Soil Sci.
  n = 1.5                   ; empirical constant 1
  b = 4.9                   ; empirical constant 2
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  if keyword_set(air) then D_cos = mol_diff * (T_k/T_ref)^0.5 $ 
  else begin
    if n_params() ne 3 then message,'ERROR: Number of parameters incorrect.'
    D_cos = mol_diff * (T_k/T_ref)^n * theta_a^2 * (theta_a / theta_sat)^(3/b)
  endelse
  return, D_cos
end

function soil_cos_sink_old, T, water_cont, cos_conc, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 3 then message,'ERROR: Number of parameters incorrect.'
  K_m = 1.9     ; Michaelis constant, mol m^-3
  V_max = vmax 
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ;V_max = 3.73e-9*0   ; for SC1, new litter water content parametrization
  ;V_max = 1.5551091e-08   ; for SC2, new litter water content parametrization
  ;k_w = 10.9  ; water dependence factor 
  ;if water_cont[0] ge 0.18 then V_max = V_max * 2  ; to simulate the burst after rain event
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  T_min = 283.15
  T_max = 296.15
  T_limit = (T_k-T_min) * (T_k-T_max)^2 / 325.4815    ; the envelope of f_cos vs. temp plot  
  bad_T_index = where(T_k le T_min or T_k ge T_max, bad_T_count)
  if bad_T_count ne 0 then T_limit[bad_T_index] = 0.
  ; outside the temperature range, enzyme activity is zero
  ;cos_sink = - V_max * cos_conc / (K_m + cos_conc) * T_limit * exp(water_cont * k_w)
  ; function form changed on 05/16/14
  k_w = 14.13  ; water dependence factor 
  cos_sink = - V_max * cos_conc / (K_m + cos_conc) * T_limit * sinh(k_w * water_cont)
  return, cos_sink
end

function soil_cos_sink, T, swc, cos_conc, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 3 then message,'ERROR: Number of parameters incorrect.'
  K_m = 1.9     ; Michaelis constant, mol m^-3
  V_max = vmax 
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  
  ; constants
  R_gas = 8.3144621
  k_B = 1.3806488e-23
  h = 6.62606957e-34
  
  ; parameters
  Delta_G_cat = 8.41e4 ; 90000 ; 25000
  Delta_H_eq = 3.59e5 ;359000 ; 180000 ; 360000
  T_0 = 273.15
  T_eq = 288.15 ; 288.15 ; 293.15
  swc_opt = 0.14
  T_opt = T_eq*(1-R_gas/Delta_H_eq*alog((Delta_H_eq-Delta_G_cat)/Delta_G_cat)*T_eq)
  norm_fac = k_B * T_opt / h * exp(-Delta_G_cat / (R_gas*T_opt) ) / (1 + exp( -Delta_H_eq / R_gas * (1/T_opt - 1/T_eq) ) ) 
  norm_fac = 1. / norm_fac
  ; normalization factor, max(T_limit) = 1
  
  T_limit = norm_fac * k_B * T_k / h * exp(-Delta_G_cat / (R_gas*T_k) ) / (1 + exp( -Delta_H_eq / R_gas * (1/T_k - 1/T_eq) ) )
  swc_limit = swc/swc_opt^2 * exp(- (swc/swc_opt)^2/2 ) / (exp(-0.5)/swc_opt)
  cos_sink = - V_max * cos_conc / (K_m + cos_conc) * T_limit * swc_limit
  return, cos_sink
end

function litter_cos_sink, T, water_cont, cos_conc, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 3 then message,'ERROR: Number of parameters incorrect.'
  K_m = 1.9     ; Michealis constant, mol m^-3
  V_max = vmax
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ;V_max = 2.1300000e-09   ; fixed
  ;if water_cont[0] ge 0.25 then V_max = V_max * 2  ; to simulate the burst after rain event
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
;  T_min = 283.15
;  T_max = 303.15
;  T_limit = sqrt(T_k-T_min) * (T_k-T_max)^2 / 512.
;  bad_T_index = where(T_k le T_min or T_k ge T_max, bad_T_count)
;  if bad_T_count ne 0 then T_limit[bad_T_index] = 0.
  ; outside the temperature range, enzyme activity is zero
  ; cos_sink = - V_max * cos_conc / (K_m + cos_conc) * T_limit * exp(water_cont/0.15)
  ;cos_sink = - V_max * cos_conc / (K_m + cos_conc) * exp(water_cont/0.15)
  ;k_w = 9.23  ; water dependence factor 
  ;cos_sink = - V_max * cos_conc / (K_m + cos_conc) * exp(water_cont * k_w)
  ; function form changed on 05/16/14
  k_w = 11.56063 ;8.0; 8.50  ; water dependence factor 
  cos_sink = - V_max * cos_conc / (K_m + cos_conc) * sinh(water_cont * k_w)
  return, cos_sink
end

; litter COS source enabled
function litter_cos_source, T, water_cont, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 2 then message,'ERROR: Number of parameters incorrect.'
  V_max = vmax
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ;V_max = 3.4638848e-11   ; fixed
  ;if water_cont[0] ge 0.25 then V_max = 0  ; to simulate the burst after rain event ; no need
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  ;;T_half = 288.15 ;286.15                               ; Half saturation temperature
  ;;T_limit = 1. / (1. + exp(T_half-T_k))    ; Logistic curve
  ;; T_limit = exp(k_T * (T-298.15))    ; not using logistic function, because have no knowledge of the actual relationship
  ; function form changed on 01/04/15
  Q10_prod = 1.9
  T_ref = 298.15
  T_func = exp( alog(Q10_prod)/10 * (T_k - T_ref) )  ; exponential curve
  cos_source = V_max * T_func      ; * (water_cont/0.10)
  return, cos_source
end

; old soil COS source parameterization: logistic function
function soil_cos_source_old, T, water_cont, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 2 then message,'ERROR: Number of parameters incorrect.'
  V_max = vmax
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ;V_max = 1.20e-11*0 ;1.2018677e-11   ; for SC1, new litter water content parametrization
  ;V_max = 5.6949159e-11   ; for SC2, new litter water content parametrization
  ;if water_cont[0] ge 0.25 then V_max = 0  ; to simulate the burst after rain event
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  T_half = 288.15 ;286.15                               ; Half saturation temperature
  ;T_limit = 1. / (1. + exp(T_half-T_k)) + 1.    ; Logistic curve
  T_limit = 1. / (1. + exp(T_half-T_k))    ; Logistic curve
  cos_source = V_max * T_limit      ; * (water_cont/0.10)
  return, cos_source
end

; new soil COS source parameterization: Q10 parameter
function soil_cos_source, T, water_cont, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 2 then message,'ERROR: Number of parameters incorrect.'
  V_max = vmax
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  Q10 = 1.9    ; An average value from Maseyk et al. (2014). PNAS; range 1.7-2.1
  ;T_half = 288.15 ;286.15                               ; Half saturation temperature
  ;T_limit = 1. / (1. + exp(T_half-T_k))    ; Logistic curve
  T_ref = 298.15
  T_func = exp( alog(Q10)/10 * (T_k - T_ref) ) ; exponential curve
  ; T_func = T_func / (1 + exp(-(T_k-288.15))) ; inhibit the production at lower temperatures
  cos_source = V_max * T_func
  return, cos_source
end

function calc_litter_resp, T, water_cont, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 2 then message,'ERROR: Number of parameters incorrect.'
  V_max = vmax
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ;V_max = 6.9259412e-07  ; for SC1, new litter water content parametrization
  ;V_max = 8.0140151e-07; 8.5978683e-07  ; for SC2, new litter water content parametrization
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  k_w = 14.99656
  E_resp = 1.790e3   ; E_resp/R, see results of lit_flux_co2.R
  litter_resp = V_max * exp(E_resp * (1./298.15 - 1./T_k)) * sinh(k_w * water_cont)      ; * (water_cont/0.10)
  return, litter_resp
end

function calc_soil_resp, T, water_cont, depth, vmax, Kelvin=kelvin
  compile_opt idl2
  if n_params() lt 2 then message,'ERROR: Number of parameters incorrect.'
  V_max = vmax
  if ~finite(V_max) then message,'ERROR: Inf or NaN value for V_max.'
  ;V_max = 1.7742266e-06  ; for SC1, new litter water content parametrization
  ;V_max = 9.0908241e-06; 8.9629374e-06  ; for SC2, new litter water content parametrization
  ; T_c: Celsius temperature, T_k: Kelvin temperature
  if keyword_set(kelvin) then T_k = T else T_k = T + 273.15
  ; water dependence function form changed on 05/16/14, exp --> sinh
  
  ; disable sinh water dependence function
  ; k_w = 14.13  ; water dependence factor
  ; swc_limit = sinh(k_w * water_cont)
  
  ; now use linear function
  swc_std = 0.10
  swc_limit = water_cont / swc_std
   
  A = 77.0683 ; pre-exponential factor to normalize the temperature dependence function
  soil_resp = V_max * A * exp(-308.56 / (T_k - 227.13) ) * swc_limit * $
    exp( - depth / 0.2 )  
  ; depth dependence comes from the limitation of oxygen fugacity in soil pore space air
  ; temperature dependence from Lloyd & Taylor 1994 Func Ecol
  return, soil_resp
end

function calc_cos_flux_NSS_FVgrid, T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt, cv_node, cv_size, n_layers, litter_layers, vmax, max_time
  ; NSS means "non-steady state"
  compile_opt idl2 
  tol = 1e-5    ; relative tolerance = 1e-5
  ; vmax input
  ; vmax[0] - soil cos uptake
  ; vmax[1] - litter cos uptake
  ; vmax[2] - soil cos source
  ; vmax[3] - litter cos source
  
  theta_a = theta_sat - theta_w
  co2_k_H = calc_co2_solub(T_soil)
  cos_k_H = calc_cos_solub(T_soil)
  D_co2_soil = calc_co2_diff(T_soil, theta_a, theta_sat)
  D_co2_atm = calc_co2_diff(T_air, /air)
  D_cos_soil = calc_cos_diff(T_soil, theta_a, theta_sat)
  D_cos_atm = calc_cos_diff(T_air, /air)
  
  ; use finite volume grid
  alpha = dblarr(1, n_layers)
  beta = dblarr(1, n_layers)
  gamma = dblarr(1, n_layers)
    
  alpha = (theta_a + cos_k_H * theta_w) * cv_size
    
  beta = D_cos_soil
  beta[0] = 2 * beta[0] * D_cos_atm / (beta[0] + D_cos_atm) ; use the harmonic mean for soil-air boundary
  beta[1:n_layers-1] = (beta[1:n_layers-1] + beta[0:n_layers-2]) / 2
  beta[0] = beta[0] / cv_node[0]
  beta[1:n_layers-1] = beta[1:n_layers-1] / (cv_node[1:n_layers-1] - cv_node[0:n_layers-2])

  for n_time = 0, 100000L do begin
    gamma = ( soil_cos_sink(T_soil, theta_w, cos_soil * cos_k_H, vmax[0]) * f_CA + soil_cos_source(T_soil, theta_w, vmax[2]) ) * cv_size      
    if litter_layers gt 0 then gamma[0:litter_layers-1] = ( litter_cos_sink(T_soil[0:litter_layers-1], lit_wc, cos_soil[0:litter_layers-1] * cos_k_H[0:litter_layers-1], vmax[1]) $
      + litter_cos_source(T_soil[0:litter_layers-1], lit_wc, vmax[3]) ) $
      * cv_size[0:litter_layers-1]
    gamma[0] = gamma[0] + beta[0] * cos_atm
    
    S = solver_tridiag_fvgrid(cos_soil, alpha, beta, gamma, dt)
    fcos_model_old = beta[0] * (cos_soil[0] - cos_atm) * 1e12
    ; surface fcos of the last timestep
    fcos_model_new = beta[0] * (S[0] - cos_atm) * 1e12
    ; surface fcos of this timestep
    if (norm(S[0:n_layers-2] - cos_soil[0:n_layers-2]) / norm(cos_soil[0:n_layers-2]) lt 1e-3 and $
      abs(fcos_model_new - fcos_model_old) / abs(fcos_model_old) lt tol) or n_time gt max_time/dt then break
    ;debug_recorder[n_time] = S[n_layers-1] - cos_soil[n_layers-1]
    cos_soil = S 
  endfor
  print, 'time to run for COS profile = ', n_time * dt, 'seconds'
  fcos_model = beta[0] * (cos_soil[0] - cos_atm) * 1e12
  if litter_layers gt 0 then fcos_model_interface = beta[litter_layers] * (cos_soil[litter_layers] - cos_soil[litter_layers-1]) * 1e12
  ; STOP   ; for debug
  if litter_layers gt 0 then return, [fcos_model, fcos_model_interface] else return, fcos_model
end

function calc_co2_flux_NSS_FVgrid, T_soil, T_air, theta_sat, theta_w, lit_wc, co2_atm, co2_soil, f_biomass, dt, cv_node, cv_size, n_layers, litter_layers, vmax, max_time
  ; NSS means "non-steady state"
  compile_opt idl2 
  tol = 1e-5   ; relative tolerance = 1e-5
  ; vmax input
  ; vmax[0] - soil resp
  ; vmax[1] - litter resp
  
  theta_a = theta_sat - theta_w
  co2_k_H = calc_co2_solub(T_soil)
  cos_k_H = calc_cos_solub(T_soil)
  D_co2_soil = calc_co2_diff(T_soil, theta_a, theta_sat)
  D_co2_atm = calc_co2_diff(T_air, /air)
  D_cos_soil = calc_cos_diff(T_soil, theta_a, theta_sat)
  D_cos_atm = calc_cos_diff(T_air, /air)
  ;lit_wc = replicate(lit_wc, 1, litter_layers)
  
  ; use finite volume grid
  alpha = dblarr(1, n_layers)
  beta = dblarr(1, n_layers)
  gamma = dblarr(1, n_layers)
    
  alpha = (theta_a + co2_k_H * theta_w) * cv_size
    
  beta = D_co2_soil
  beta[0] = 2 * beta[0] * D_co2_atm / (beta[0] + D_co2_atm) ; use the harmonic mean for soil-air boundary
  beta[1:n_layers-1] = (beta[1:n_layers-1] + beta[0:n_layers-2]) / 2
  beta[0] = beta[0] / cv_node[0]
  beta[1:n_layers-1] = beta[1:n_layers-1] / (cv_node[1:n_layers-1] - cv_node[0:n_layers-2])
  ;debug_recorder = dblarr(100001L)
  for n_time = 0L, 100000L do begin
    gamma = calc_soil_resp(T_soil, theta_w, cv_node-(cv_node[litter_layers-1]+cv_node[litter_layers])/2., vmax[0]) * f_biomass $
      * cv_size 
    if litter_layers gt 0 then gamma[0:litter_layers-1] = calc_litter_resp(T_soil[0:litter_layers-1], lit_wc, vmax[1]) * cv_size
    gamma[0] = gamma[0] + beta[0] * co2_atm
    S = solver_tridiag_fvgrid(co2_soil, alpha, beta, gamma, dt)
    fco2_model_old = beta[0] * (co2_soil[0] - co2_atm) * 1e6
    ; surface fco2 of the last timestep
    fco2_model_new = beta[0] * (S[0] - co2_atm) * 1e6
    ; surface fco2 of this timestep
    if (norm(S[0:n_layers-2] - co2_soil[0:n_layers-2]) / norm(co2_soil[0:n_layers-2]) lt 1e-3 and $
      abs(fco2_model_new - fco2_model_old) / abs(fco2_model_old) lt tol) or n_time gt max_time/dt then break
    ;debug_recorder[n_time] = S[n_layers-1] - co2_soil[n_layers-1]
    co2_soil = S 
  endfor
  print, 'time to run for CO2 profile = ', n_time * dt, 'seconds'
  fco2_model = beta[0] * (co2_soil[0] - co2_atm) * 1e6
  if litter_layers gt 0 then fco2_model_interface = beta[litter_layers] * (co2_soil[litter_layers] - co2_soil[litter_layers-1]) * 1e6
  ; STOP   ; for debug
  ; STOP
  ;return, fco2_model
  if litter_layers gt 0 then return, [fco2_model, fco2_model_interface] else return, fco2_model
end

pro sens_test
  compile_opt idl2
  
;;___________________________________________________
;; plot parameters
  !p.background='ffffff'xl
  !p.color='000000'xl
  !x.thick=2
  !y.thick=2
  !p.charsize=2
  !p.thick=2
  !p.font=1
  !x.style=1  
  ;!p.font=0
  ;set_plot, 'ps'
  ;device, helvetica=1
  ;device, isolatin1=1

;;___________________________________________________
;; loading data  

  ; data_dir = './input/'
  plot_dir = './plots/mdl_plts/'
  ; data_file = data_dir + 'soil_flux_regular_meas.csv'
  
  
  ;;___________________________________________________
  ;; set parameters
  ; set physical constants
  p_atm = 1.01325e5           ; Surface pressure (kPa), assumed for the SGP site
  p_std = 1.01325e5           ; Standard atmosphere pressure (N m^-2)
  R_gas = 8.3144621        ; Universal gas constant (J K^-1 mol^-1)
  air_conc = p_std / R_gas / 273.15D    ; Air molar concentration (mol m^-3) at the standard condition (273.15 K and 101.325 kPa)
  Mw = 18.016     ; molecular weight of water 
  
  ; set chamber dimensions
  V_sc = 4076.1D / 1e6                    ; Soil chamber (Licor) volume (m^3), according to website, bowl only  
  A_sc = 317.8D / 1e4                  ; Soil area covered by the chamber (m^2)
  ; sc_collar_ht = 0.     ; Soil chamber collar height (m)
  ; V_sc = V_sc + sc_collar_ht * A_sc    ; Soil chamber volume (m^3)
  
  ; set model configurations
  dt = 10D          ; timestep 10 s
  z_max = 1D        ; lower boundary depth at 1 m
  n_layers = 26     ; number of layers = 26 in fast mode
  litter_layers = 6 ; number of litter layers (0:5 for sc1, 0:10 for sc2)
  cv_node = z_max * exp(0.2 * dindgen(n_layers) - 5)
  cv_size = dindgen(n_layers)
  cv_size[0] = (cv_node[0] + cv_node[1]) / 2
  cv_size[1:n_layers-2] = (cv_node[2:n_layers-1] - cv_node[0:n_layers-3]) / 2
  cv_size[n_layers-1] = cv_node[n_layers-1] - cv_node[n_layers-2]
  cv_face = dindgen(n_layers)
  cv_face[0:n_layers-2] = (cv_node[1:n_layers-1] + cv_node[0:n_layers-2]) / 2
  cv_face[n_layers-1] = cv_node[n_layers-1] + cv_size[n_layers-1] / 2
  
  z_grid = cv_node
  
  test_flag_1 = 0   ; run test for cos uptake vs water filled pore space
  test_flag_2 = 0   ; run test for cos uptake vs vmax
  test_flag_3 = 1   ; run test for cos surface uptake vs vmax for each depth

  if test_flag_1 then begin
    theta_sat = replicate(0.35, n_layers)
    theta_w_test = sin(!dpi * findgen(101)/100 / 2)^2 * theta_sat[0]
    fcos_model = dblarr(101)
    fcos_model_hitemp = dblarr(101)
    fcos_model_lowtemp = dblarr(101)
    fcos_model_medhitemp = dblarr(101)
    fco2_model = dblarr(101)
    vmax_cos_test = [1e-2, 1e-3]
    dt = 1D
    for k = 0, 100 do begin 
      theta_w = replicate(theta_w_test[k], n_layers)
      theta_a = theta_sat - theta_w
      T_soil = replicate(15D, n_layers)
      T_air = 25D
      co2_atm = 400D  ; for ideal conditions, 400 ppmv
      cos_atm = 500D  ; for ideal conditions, 500 pptv
      co2_atm = co2_atm * 1e-6 * air_conc * 273.15D / (T_soil[0]+273.15D)   ; convert to mol m^-3
      cos_atm = cos_atm * 1e-12 * air_conc * 273.15D / (T_soil[0]+273.15D)  ; convert to mol m^-3
      co2_soil = replicate(co2_atm, 1, n_layers)   ; set initial conditions for soil co2 profile
      cos_soil = replicate(cos_atm, 1, n_layers)   ; set initial conditions for soil cos profile
      co2_k_H = calc_co2_solub(T_soil)
      cos_k_H = calc_cos_solub(T_soil)
      D_co2_soil = calc_co2_diff(T_soil, theta_a, theta_sat)
      D_co2_atm = calc_co2_diff(T_air, /air)
      D_cos_soil = calc_cos_diff(T_soil, theta_a, theta_sat)
      D_cos_atm = calc_cos_diff(T_air, /air)
      f_CA = replicate(1D, 1, n_layers)
      
      ; use finite volume grid
      alpha = dblarr(1, n_layers)
      beta = dblarr(1, n_layers)
      gamma = dblarr(1, n_layers)
        
      alpha = (theta_a + cos_k_H * theta_w) * cv_size
        
      beta = D_cos_soil
      beta[0] = 2 * beta[0] * D_cos_atm / (beta[0] + D_cos_atm) ; use the harmonic mean for soil-air boundary
      beta[1:n_layers-1] = (beta[1:n_layers-1] + beta[0:n_layers-2]) / 2
      beta[0] = beta[0] / cv_node[0]
      beta[1:n_layers-1] = beta[1:n_layers-1] / (cv_node[1:n_layers-1] - cv_node[0:n_layers-2])
      
      for n_time = 0, 100000L do begin
        gamma = soil_cos_sink(T_soil, theta_w, cos_soil * cos_k_H, vmax_cos_test[0]) * f_CA $
          * cv_size
        gamma[0] = gamma[0] + beta[0] * cos_atm
        
        S = solver_tridiag_fvgrid(cos_soil, alpha, beta, gamma, dt)
        if norm(S[0:n_layers-2] - cos_soil[0:n_layers-2]) / norm(cos_soil[0:n_layers-2]) lt 1e-6 then break ; relative tolerance = 1e-6
        cos_soil = S 
      endfor
      fcos_model[k] = beta[0] * (cos_soil[0] - cos_atm) * 1e12 
    endfor
    
    print, "Sensitivity test 1/4 done..."
    
    ; at a higher temperature
    for k = 0, 100 do begin 
      theta_w = replicate(theta_w_test[k], n_layers)
      theta_a = theta_sat - theta_w
      T_soil = replicate(22D, n_layers)
      T_air = 25D
      co2_atm = 400D  ; for ideal conditions, 400 ppmv
      cos_atm = 500D  ; for ideal conditions, 500 pptv
      co2_atm = co2_atm * 1e-6 * air_conc * 273.15D / (T_soil[0]+273.15D)   ; convert to mol m^-3
      cos_atm = cos_atm * 1e-12 * air_conc * 273.15D / (T_soil[0]+273.15D)  ; convert to mol m^-3
      co2_soil = replicate(co2_atm, 1, n_layers)   ; set initial conditions for soil co2 profile
      cos_soil = replicate(cos_atm, 1, n_layers)   ; set initial conditions for soil cos profile
      co2_k_H = calc_co2_solub(T_soil)
      cos_k_H = calc_cos_solub(T_soil)
      D_co2_soil = calc_co2_diff(T_soil, theta_a, theta_sat)
      D_co2_atm = calc_co2_diff(T_air, /air)
      D_cos_soil = calc_cos_diff(T_soil, theta_a, theta_sat)
      D_cos_atm = calc_cos_diff(T_air, /air)
      f_CA = replicate(1D, 1, n_layers)
      
      ; use finite volume grid
      alpha = dblarr(1, n_layers)
      beta = dblarr(1, n_layers)
      gamma = dblarr(1, n_layers)
        
      alpha = (theta_a + cos_k_H * theta_w) * cv_size
        
      beta = D_cos_soil
      beta[0] = 2 * beta[0] * D_cos_atm / (beta[0] + D_cos_atm) ; use the harmonic mean for soil-air boundary
      beta[1:n_layers-1] = (beta[1:n_layers-1] + beta[0:n_layers-2]) / 2
      beta[0] = beta[0] / cv_node[0]
      beta[1:n_layers-1] = beta[1:n_layers-1] / (cv_node[1:n_layers-1] - cv_node[0:n_layers-2])
      
      for n_time = 0, 100000L do begin
        gamma = soil_cos_sink(T_soil, theta_w, cos_soil * cos_k_H, vmax_cos_test[0]) * f_CA $
          * cv_size
        gamma[0] = gamma[0] + beta[0] * cos_atm
        
        S = solver_tridiag_fvgrid(cos_soil, alpha, beta, gamma, dt*5)
        if norm(S[0:n_layers-2] - cos_soil[0:n_layers-2]) / norm(cos_soil[0:n_layers-2]) lt 1e-6 then break ; relative tolerance = 1e-6
        cos_soil = S 
      endfor
      fcos_model_hitemp[k] = beta[0] * (cos_soil[0] - cos_atm) * 1e12
    endfor
    
    print, "Sensitivity test 2/4 done..."
    
    ; at a lower temperature
    for k = 0, 100 do begin 
      theta_w = replicate(theta_w_test[k], n_layers)
      theta_a = theta_sat - theta_w
      T_soil = replicate(13D, n_layers)
      T_air = 25D
      co2_atm = 400D  ; for ideal conditions, 400 ppmv
      cos_atm = 500D  ; for ideal conditions, 500 pptv
      co2_atm = co2_atm * 1e-6 * air_conc * 273.15D / (T_soil[0]+273.15D)   ; convert to mol m^-3
      cos_atm = cos_atm * 1e-12 * air_conc * 273.15D / (T_soil[0]+273.15D)  ; convert to mol m^-3
      co2_soil = replicate(co2_atm, 1, n_layers)   ; set initial conditions for soil co2 profile
      cos_soil = replicate(cos_atm, 1, n_layers)   ; set initial conditions for soil cos profile
      co2_k_H = calc_co2_solub(T_soil)
      cos_k_H = calc_cos_solub(T_soil)
      D_co2_soil = calc_co2_diff(T_soil, theta_a, theta_sat)
      D_co2_atm = calc_co2_diff(T_air, /air)
      D_cos_soil = calc_cos_diff(T_soil, theta_a, theta_sat)
      D_cos_atm = calc_cos_diff(T_air, /air)
      f_CA = replicate(1D, 1, n_layers)
      
      ; use finite volume grid
      alpha = dblarr(1, n_layers)
      beta = dblarr(1, n_layers)
      gamma = dblarr(1, n_layers)
        
      alpha = (theta_a + cos_k_H * theta_w) * cv_size
        
      beta = D_cos_soil
      beta[0] = 2 * beta[0] * D_cos_atm / (beta[0] + D_cos_atm) ; use the harmonic mean for soil-air boundary
      beta[1:n_layers-1] = (beta[1:n_layers-1] + beta[0:n_layers-2]) / 2
      beta[0] = beta[0] / cv_node[0]
      beta[1:n_layers-1] = beta[1:n_layers-1] / (cv_node[1:n_layers-1] - cv_node[0:n_layers-2])
      
      for n_time = 0, 100000L do begin
        gamma = soil_cos_sink(T_soil, theta_w, cos_soil * cos_k_H, vmax_cos_test[0]) * f_CA $
          * cv_size
        gamma[0] = gamma[0] + beta[0] * cos_atm
        
        S = solver_tridiag_fvgrid(cos_soil, alpha, beta, gamma, dt*3)
        if norm(S[0:n_layers-2] - cos_soil[0:n_layers-2]) / norm(cos_soil[0:n_layers-2]) lt 1e-6 then break ; relative tolerance = 1e-6
        cos_soil = S 
      endfor
      fcos_model_lowtemp[k] = beta[0] * (cos_soil[0] - cos_atm) * 1e12 
    endfor
    
    print, "Sensitivity test 3/4 done..."
    
    ; at T=20 C
    for k = 0, 100 do begin 
      theta_w = replicate(theta_w_test[k], n_layers)
      theta_a = theta_sat - theta_w
      T_soil = replicate(20D, n_layers)
      T_air = 25D
      co2_atm = 400D  ; for ideal conditions, 400 ppmv
      cos_atm = 500D  ; for ideal conditions, 500 pptv
      co2_atm = co2_atm * 1e-6 * air_conc * 273.15D / (T_soil[0]+273.15D)   ; convert to mol m^-3
      cos_atm = cos_atm * 1e-12 * air_conc * 273.15D / (T_soil[0]+273.15D)  ; convert to mol m^-3
      co2_soil = replicate(co2_atm, 1, n_layers)   ; set initial conditions for soil co2 profile
      cos_soil = replicate(cos_atm, 1, n_layers)   ; set initial conditions for soil cos profile
      co2_k_H = calc_co2_solub(T_soil)
      cos_k_H = calc_cos_solub(T_soil)
      D_co2_soil = calc_co2_diff(T_soil, theta_a, theta_sat)
      D_co2_atm = calc_co2_diff(T_air, /air)
      D_cos_soil = calc_cos_diff(T_soil, theta_a, theta_sat)
      D_cos_atm = calc_cos_diff(T_air, /air)
      f_CA = replicate(1D, 1, n_layers)
      
      ; use finite volume grid
      alpha = dblarr(1, n_layers)
      beta = dblarr(1, n_layers)
      gamma = dblarr(1, n_layers)
        
      alpha = (theta_a + cos_k_H * theta_w) * cv_size
        
      beta = D_cos_soil
      beta[0] = 2 * beta[0] * D_cos_atm / (beta[0] + D_cos_atm) ; use the harmonic mean for soil-air boundary
      beta[1:n_layers-1] = (beta[1:n_layers-1] + beta[0:n_layers-2]) / 2
      beta[0] = beta[0] / cv_node[0]
      beta[1:n_layers-1] = beta[1:n_layers-1] / (cv_node[1:n_layers-1] - cv_node[0:n_layers-2])
      
      for n_time = 0, 100000L do begin
        gamma = soil_cos_sink(T_soil, theta_w, cos_soil * cos_k_H, vmax_cos_test[0]) * f_CA $
          * cv_size
        gamma[0] = gamma[0] + beta[0] * cos_atm
        
        S = solver_tridiag_fvgrid(cos_soil, alpha, beta, gamma, dt*3)
        if norm(S[0:n_layers-2] - cos_soil[0:n_layers-2]) / norm(cos_soil[0:n_layers-2]) lt 1e-6 then break ; relative tolerance = 1e-6
        cos_soil = S 
      endfor
      fcos_model_medhitemp[k] = beta[0] * (cos_soil[0] - cos_atm) * 1e12 
    endfor
    
    print, "Sensitivity test 4/4 done."
    stop
    p_vmax=plot(theta_w_test/theta_sat[0], $
      -soil_cos_sink(15., theta_w_test, cos_atm * calc_cos_solub(15.), vmax_cos_test[0]) / (cos_atm * calc_cos_solub(15.)), $
      'g-', color='005a32'xl, thick=1.5, xrange=[0,1], xtitle='Water-filled pore space fraction', $
      ytitle='Microbial COS turnover rate (s$^{-1}$)', margin=[0.26, 0.10, 0.14, 0.10], dimensions=[720,480], yrange=[0,0.0055])
    p_vmax_hitemp = plot(theta_w_test/theta_sat[0], $
      -soil_cos_sink(22., theta_w_test, cos_atm * calc_cos_solub(22.), vmax_cos_test[0]) / (cos_atm * calc_cos_solub(22.)), $
      'g', color='74c476'xl, thick=1.5, /current, overplot=1)
    p_vmax_lowtemp = plot(theta_w_test/theta_sat[0], $
      -soil_cos_sink(13., theta_w_test, cos_atm * calc_cos_solub(13.), vmax_cos_test[0]) / (cos_atm * calc_cos_solub(13.)), $
      'g', color='238b45'xl, thick=1.5, /current, overplot=1)
    p_vmax_medhitemp = plot(theta_w_test/theta_sat[0], $
      -soil_cos_sink(20., theta_w_test, cos_atm * calc_cos_solub(20.), vmax_cos_test[0]) / (cos_atm * calc_cos_solub(20.)), $
      'g', color='41ab5d'xl, thick=1.5, /current, overplot=1)
    ax = p_vmax.axes
    ;ax[0].hide = 1
    ax[1].hide = 1
    ;ax[2].hide = 1
    ax[3].showtext = 1
    ax[3].color='g'
    p_diff=plot(theta_w_test/theta_sat[0], calc_cos_diff(15., theta_sat[0]-theta_w_test, theta_sat[0]), 'b-', thick=1.5, color='084594'xl, $
      xrange=[0,1], ytitle='COS diffusivity in soil (m$^2$ s$^{-1}$)', margin=[0.26, 0.10, 0.14, 0.10], /current)
    p_diff_hitemp=plot(theta_w_test/theta_sat[0], calc_cos_diff(22., theta_sat[0]-theta_w_test, theta_sat[0]), $
      'b-', color='6baed6'xl, thick=1.5, /current, overplot=1)
    p_diff_lowtemp=plot(theta_w_test/theta_sat[0], calc_cos_diff(13., theta_sat[0]-theta_w_test, theta_sat[0]), $
      'b-', color='2171b5'xl, thick=1.5, /current, overplot=1)
    p_diff_medhitemp=plot(theta_w_test/theta_sat[0], calc_cos_diff(20., theta_sat[0]-theta_w_test, theta_sat[0]), $
      'b-', color='4292c6'xl, thick=1.5, /current, overplot=1)
    
    ax = p_diff.axes
    ax[0].hide = 1
    ax[2].hide = 1
    ax[3].hide = 1
    ax[1].location=[-0.2, 0, 0]
    ax[1].color='b'
    
    
    ;p = plot(cos_soil / (air_conc * 273.15D / (T_soil+273.15D) ) * 1e12, -z_grid, yrange=[-1,0],'k+-', sym_color='red')
    ;p = plot(findgen(n_time) * dt, fcos_model)
    ;STOP
    p_swc=plot(theta_w_test/theta_sat[0], fcos_model, 'k-', color='7a0177'xl, ytitle='Soil COS flux (pmol m$^{-2}$ s$^{-1}$)', $
      thick=2, xrange=[0,1], yrange=[0,-1.6], margin=[0.26, 0.10, 0.14, 0.10], /current)
      
    p_swc_hitemp=plot(theta_w_test/theta_sat[0], fcos_model_hitemp, 'k-', color='fa9fb5'xl, $
      thick=2, /current, overplot=1)
      
    p_swc_lowtemp=plot(theta_w_test/theta_sat[0], fcos_model_lowtemp, 'k-', color='c51b8a'xl, $
      thick=2, /current, overplot=1)
      
    p_swc_medhitemp=plot(theta_w_test/theta_sat[0], fcos_model_medhitemp, 'k-', color='f768a1'xl, $
      thick=2, /current, overplot=1)
      
    ax = p_swc.axes
    ax[0].hide = 1
    ax[2].hide = 1
    ax[3].hide = 1 ; hide right y axis
    ax[1].color = 'purple'
    ;p_vmax=plot(theta_w_test/theta_sat[0], -soil_cos_sink(T_soil[0], theta_w_test, cos_atm * calc_cos_solub(T_soil[0]), vmax_cos_test[0]), 'g-', thick=1.5, $
    ;  ytitle='COS uptake rate (mol m$^{-3}$ s$^{-1}$)', margin=[0.26, 0.10, 0.14, 0.10], /current)
    
    txt1 = text(0.38, 0.56, 'T=15$\deg$C', color='7a0177'xl, font_size=12)
    txt2 = text(0.31, 0.73, 'T=13$\deg$C', color='c51b8a'xl, font_size=12)
    txt3 = text(0.33, 0.35, 'T=20$\deg$C', color='f768a1'xl, font_size=12)
    txt4 = text(0.34, 0.28, 'T=22$\deg$C', color='fa9fb5'xl, font_size=12)
    
    txt_tilted1 = text(0.72, 0.37, 'T=15$\deg$C', color='005a32'xl, font_size=12, orientation=-50)
    txt_tilted2 = text(0.74, 0.47, 'T=13$\deg$C', color='238b45'xl, font_size=12, orientation=-45)
    txt_tilted3 = text(0.72, 0.22, 'T=20$\deg$C', color='41ab5d'xl, font_size=12, orientation=-15)
    txt_tilted4 = text(0.77, 0.14, 'T=22$\deg$C', color='74c476'xl, font_size=12, orientation=-10)
    
    p_swc.save, strjoin([plot_dir, 'sens_test_flux_vs_swc', '.eps'], /single), border=10, resolution=300
    p_swc.save, strjoin([plot_dir, '/png/sens_test_flux_vs_swc', '.png'], /single), border=10, resolution=300
  endif
  
  stop
  
  if test_flag_2 then begin
    litter_layers = 0
    lit_wc = 0.
    theta_w = replicate(0.07, 1, n_layers)
    theta_sat = replicate(0.35D, 1, n_layers)
    theta_a = theta_sat - theta_w
    T_soil = replicate(15., 1, n_layers)
    T_air = 25.
    co2_atm = 400.  ; in ppmv
    cos_atm = 500.  ; in pptv
    co2_atm = co2_atm * 1e-6 * air_conc * p_atm/p_std * 273.15D / (T_soil[0] + 273.15D)  ; convert to mol m^-3
    cos_atm = cos_atm * 1e-12 * air_conc * p_atm/p_std * 273.15D / (T_soil[0] + 273.15D) ; convert to mol m^-3
    co2_soil = replicate(co2_atm, 1, n_layers) * exp(3*z_grid)
    cos_soil = replicate(cos_atm, 1, n_layers) * exp(-0.5*z_grid)  ;* exp(-z_grid/0.4)
    f_CA = replicate(1D, 1, n_layers)
    f_biomass = replicate(1D, 1, n_layers)
    
    max_time = dt*1e6
    fcos_test_flux_vs_vmax = dblarr(101)
    vmax_cos_test = [1e-2, 0., 0., 0.] /2.
    flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt/10., cv_node, cv_size, n_layers, litter_layers, vmax_cos_test, max_time)
    fcos_test_flux_vs_vmax[0] = flux_eval[0]
    for jk = 1, 100 do begin
      vmax_cos_test = [1e-2, 0., 0., 0.] * jk  
      flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt/10., cv_node, cv_size, n_layers, litter_layers, vmax_cos_test, max_time)
      fcos_test_flux_vs_vmax[jk] = flux_eval[0]
    endfor
    stop
    p=plot( [5e-3, 1e-2 * (findgen(100)+1.)], fcos_test_flux_vs_vmax, thick=2, $ 
      yrange=[0,-5], xrange=1e-2 * [0.5,100], xtitle='Soil COS uptake capacity, $V_{SU,max}$ (mol m$^{-3}$ s$^{-1}$)', $
      ytitle='Surface COS flux (pmol m$^{-2}$ s$^{-1}$)', /xlog, dimensions=[480,400], margin=[0.15,0.15,0.05,0.05])
    ;t1 = text(0.22, 0.62, '$\theta_{sat}=0.35, \theta_w=0.07$!C!C$T_S=15\deg$C!C!C[COS]$_{atm}$=500 pptv')
    t1_ln1 = text(0.22, 0.8, '$\theta_{sat}=0.35, \theta_w=0.07$')
    t1_ln2 = text(0.22, 0.7, '$T_S=15\deg$C')
    t1_ln3 = text(0.22, 0.6, '[COS]$_{atm}$=500 pptv')
    p.save, strjoin([plot_dir, 'sens_test_flux_vs_vmax', '.eps'], /single), border=10, resolution=300
    p.save, strjoin([plot_dir, '/png/sens_test_flux_vs_vmax', '.png'], /single), border=10, resolution=300
  endif  
  stop
  
  if test_flag_3 then begin
    litter_layers = 0
    lit_wc = 0.
    theta_w = replicate(0.07, 1, n_layers)
    theta_sat = replicate(0.35D, 1, n_layers)
    theta_a = theta_sat - theta_w
    T_soil = replicate(15., 1, n_layers)
    T_air = 25.
    co2_atm = 400.  ; in ppmv
    cos_atm = 500.  ; in pptv
    co2_atm = co2_atm * 1e-6 * air_conc * p_atm/p_std * 273.15D / (T_soil[0] + 273.15D)  ; convert to mol m^-3
    cos_atm = cos_atm * 1e-12 * air_conc * p_atm/p_std * 273.15D / (T_soil[0] + 273.15D) ; convert to mol m^-3
    co2_soil = replicate(co2_atm, 1, n_layers) * exp(3*z_grid)
    cos_soil = replicate(cos_atm, 1, n_layers) * exp(-0.5*z_grid)  ;* exp(-z_grid/0.4)
    f_CA = replicate(1D, 1, n_layers)
    f_biomass = replicate(1D, 1, n_layers)
    
    max_time = dt*1e6
    fcos_vs_vmax_dz = dblarr(n_layers)
    vmax_cos_test = [1e-2, 0., 0., 0.]
    flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt, cv_node, cv_size, n_layers, litter_layers, vmax_cos_test, max_time)
    fcos_vs_vmax_dz_0 = flux_eval[0]
    for jk = 0, n_layers-1 do begin
      f_CA = replicate(1D, 1, n_layers)
      f_CA[jk] = 10. 
      flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt, cv_node, cv_size, n_layers, litter_layers, vmax_cos_test, max_time)
      fcos_vs_vmax_dz[jk] = flux_eval[0]
    endfor
    
    vmax_cos_test_with_source = [1e-2, 0., 1e-11, 0.]
    fcos_vs_vmax_dz_with_source = dblarr(n_layers)
    flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt, cv_node, cv_size, n_layers, litter_layers, vmax_cos_test_with_source, max_time)
    fcos_vs_vmax_dz_with_source_0 = flux_eval[0]
    for jk = 0, n_layers-1 do begin
      f_CA = replicate(1D, 1, n_layers)
      f_CA[jk] = 10. 
      flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt, cv_node, cv_size, n_layers, litter_layers, vmax_cos_test_with_source, max_time)
      fcos_vs_vmax_dz_with_source[jk] = flux_eval[0]
    endfor
    stop
    
    p = plot((fcos_vs_vmax_dz/fcos_vs_vmax_dz_0-1)*100., z_grid, 'ko-', thick=1.5, sym_thick=1.5, yrange=[1.0,0.005], xrange=[-5,160], /ylog, $
      dimensions=[480,480], margin=[0.15,0.15,0.05,0.05])
    p = plot((fcos_vs_vmax_dz_with_source/fcos_vs_vmax_dz_with_source_0-1)*100., z_grid, 'ro-', thick=1.5, sym_thick=1.5, /current, overplot=1)
    p.xtitle='Percentage increase of surface uptake with!C 10x soil COS uptake capacity at each layer'
    p.ytitle='Depth (m)'
    ;t1 = text(0.35, 0.50, '$\theta_{sat}=0.35, \theta_w=0.07$!C!C$T_S=15\deg$C!C!C[COS]$_{atm}$=500 pptv', font_size=11)
    t1_ln1 = text(0.35, 0.64, '$\theta_{sat}=0.35, \theta_w=0.07$')
    t1_ln2 = text(0.35, 0.57, '$T_S=15\deg$C')
    t1_ln3 = text(0.35, 0.5, '[COS]$_{atm}$=500 pptv')
    leg = legend(position=[0.45, 0.45], font_size=11, transparency=100) 
    leg[0].label = '  $V_{SU, max}=10^{-2}$, $V_{SP, max}=0$'
    leg[1].label = '  $V_{SU, max}=10^{-2}$, $V_{SP, max}=10^{-11}$ (mol m$^{-3}$ s$^{-1}$)'
    
    p.save, strjoin([plot_dir, 'sens_test_flux_vs_vmax_dz', '.eps'], /single), border=10, resolution=300
    p.save, strjoin([plot_dir, '/png/sens_test_flux_vs_vmax_dz', '.png'], /single), border=10, resolution=300
  endif  
  stop
END