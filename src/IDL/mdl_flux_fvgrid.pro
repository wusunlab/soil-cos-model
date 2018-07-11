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
  b = 5.3 ; 4.9             ; empirical constant 2, 4.9 for sandy loam, 5.3 for silt loam
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
  b = 5.3 ; 4.9             ; empirical constant 2, 4.9 for sandy loam, 5.3 for silt loam
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
  T_eq = 283.15 ;288.15 ; 288.15 ; 293.15
  swc_opt = 0.20 ;0.14
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

pro Mdl_flux_FVgrid
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

  data_dir = './input/'
  plot_dir = './plots/mdl_plts/'
  data_file = data_dir + 'soil_flux_regular_meas.csv'
  
  data_query = query_ascii(data_file)
  if data_query then begin 
;    header_name = ['date_local', 'doy_utc doy_local', 'soil_moist_Avg', 'soil_moist_Std', 'soil_temp_C_Avg', $
;      'soil_temp_C_Std', 'COS', 'CO2', 'H2O', 'fCOS_lower', 'fCOS_corr', 'fCOS_upper']
    sgp_flux_data = read_csv(data_file, header=header_name, missing_value=!Values.d_NaN)
    doy_utc = sgp_flux_data.(where(header_name eq 'doy_utc'))
    doy_local = sgp_flux_data.(where(header_name eq 'doy_local'))
    swc = sgp_flux_data.(where(header_name eq 'swc_avg'))
    T_s = sgp_flux_data.(where(header_name eq 'soil_temp_avg'))
    cos_a = sgp_flux_data.(where(header_name eq 'cos_a'))
    co2_a = sgp_flux_data.(where(header_name eq 'co2_a'))
    h2o_a = sgp_flux_data.(where(header_name eq 'h2o_a'))
    fcos = sgp_flux_data.(where(header_name eq 'fcos_corr'))
    fcos_lower = sgp_flux_data.(where(header_name eq 'fcos_lower'))
    fcos_upper = sgp_flux_data.(where(header_name eq 'fcos_upper'))
    fco2 = sgp_flux_data.(where(header_name eq 'fco2'))
    fh2o = sgp_flux_data.(where(header_name eq 'fh2o'))
  endif else begin
    print, 'Data file not found'
    stop
  endelse
  
  ; clean the data
  doy_local[where(doy_local lt 90)] = !Values.f_NaN
  bad_conc_index = where(cos_a lt 300. or co2_a lt 350. or h2o_a lt 6.) 
  cos_a[bad_conc_index] = !Values.d_NaN
  co2_a[bad_conc_index] = !Values.d_NaN
  h2o_a[bad_conc_index] = !Values.d_NaN
  
  finite_index = where(doy_local gt 90. and finite(fcos))
  
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
  
;  ; ========================================================
;  ; new soil temperature parameterization
;  ; interpolate chamber temperatures
  
  ; read gap filled soil temperature data
  T_soil_gap_filled = read_csv(data_dir + 'T_soil_gap_filled.csv')
  doy_local_full = T_soil_gap_filled.field1
  Ts_full = T_soil_gap_filled.field2
  Ts_smooth = smooth(Ts_full, 1441, /NaN, /edge_truncate) ; averaging over 15 days
  ; p=plot(doy_local_full, Ts_full)
  ; p=plot(doy_local_full, Ts_smooth, 'r-', /current, overplot=1)
  Ts_amp = Ts_full - Ts_smooth   ; amplitude of high-frequency variation
  damp_depth = 0.11 
  phase_lag = 1. / (2.*!dpi) / damp_depth * 96.
  T_soil_profile = dblarr(n_elements(Ts_full), n_layers)
  phase_lag_index = round( z_grid * phase_lag )  ; with respect to approx. 5cm soil depth
  T_soil_profile[*, 0] = Ts_full
  for i = 1, 9 do T_soil_profile[*, i] = Ts_full  ; constant profile in approx. upper 5 cm
  for j = 10, n_layers-1 do begin
    T_soil_profile[0:phase_lag_index[j],j] = Ts_amp[0]
    T_soil_profile[phase_lag_index[j]+1:n_elements(Ts_full)-1, j] = Ts_amp[1:n_elements(Ts_full)-1-phase_lag_index[j]]
    T_soil_profile[*, j] = T_soil_profile[*, j] * exp(- (z_grid[j]-z_grid[9]) /damp_depth) + Ts_smooth
  endfor
  ; find the indices
  select_Ts_index = []
  for i = 0, n_elements(finite_index)-1 do begin
    current_index = where( abs(doy_local_full - doy_local[finite_index[i]]) lt 5./1440. )
    if current_index ne -1 then select_Ts_index = [select_Ts_index, current_index]
  endfor 
  T_soil_profile_full = T_soil_profile
  T_soil_profile = T_soil_profile_full[select_Ts_index, *]
  
  stop
;  ; test figure
;  p=plot(doy_local_full, Ts_full) 
;  p=plot(doy_local_full, Ts_smooth, 'r-', /current, overplot=1)
;  ; test figure
;  loadct, 7, rgb_table = temp_ct
;  temp_ct = reverse(temp_ct)
;  p_ts = contour([[Ts_full],[T_soil_profile]], doy_local_full, [0.,z_grid], /fill, rgb_table=temp_ct, c_value=findgen(40)+7, $  ;old c_value = findgen(18)+10
;    title='Soil temperature ($\deg$C)', xtitle="Day of Year", ytitle="depth (m)", xrange=[90,153], yrange=[0.5,0], dimensions=[1200,500], $
;    margin = [0.07, 0.23, 0.07, 0.09], font_size=14)
;  cbar = colorbar(target=p_ts, title='Soil temperature ($\deg$C)', font_size=12)
;  cbar.position=[0.10, 0.085, 0.9, 0.125] ;[0.15,0.05 ,0.9, 0.09]
;  
  
  ; use constant temperature profiles
;  T_soil_profile = dblarr(n_elements(T_s), n_layers)
;  for i = 0, n_elements(T_s)-1 do T_soil_profile[i,*] = T_s[i]
  
  
  ; interpolate soil moisture data to get profiles
  ; read met station data
  met_soil_file = './input/met_station_soil_data.csv'
  met_soil_data = read_csv(met_soil_file, n_table_header=0)
  filter_met_index = where(met_soil_data.field2 gt 0. and met_soil_data.field6 gt 0.)
  doy_local_met = met_soil_data.field1[filter_met_index]
  Ts_met_5 = met_soil_data.field2[filter_met_index]
  Ts_met_15 = met_soil_data.field3[filter_met_index]
  Ts_met_25 = met_soil_data.field4[filter_met_index]
  swc_met_10 = met_soil_data.field5[filter_met_index]
  swc_met_30 = met_soil_data.field6[filter_met_index]
  
  ; interpolate to the chamber data time axis
  swc_30 = interpol(swc_met_30, doy_local_met, doy_local, /NaN) / 100.
  
  ; create soil moisture profile
  swc_profile = dblarr(n_elements(finite_index), n_layers) 
  depth_ind_5cm = where( abs(z_grid-0.05) eq min(abs(z_grid-0.05)))
  depth_ind_30cm = where( abs(z_grid-0.3) eq min(abs(z_grid-0.3)))
  depth_ind_5cm = depth_ind_5cm[0]
  depth_ind_30cm = depth_ind_30cm[0]
  for i = 0, depth_ind_5cm do $
    swc_profile[*, i] = swc[ finite_index ] 
  for i = depth_ind_30cm, n_layers-1 do $ 
    swc_profile[*, i] = swc_30[ finite_index ] 
  ; interpolate between 5 cm and 30 cm using a sigmoid-like function
  for i = depth_ind_5cm+1, depth_ind_30cm-1 do $
    swc_profile[*, i] = swc[ finite_index ] + (swc_30[ finite_index ] - swc[ finite_index ])/(1 + exp( - (z_grid[i]-0.20)/0.03 ) )
    
  stop
  
  vmax_cos = [5e-1, 0., 2e-11, 0.] ; 04/10/2015, test
  vmax_co2 = [2.4530081e-06, 0.]  ; 04/12/2015, test 
  
  flag = intarr(1)
  read, flag, prompt='Calculate SGP cos fluxes? (1 or 0)'
  if flag then begin
  fcos_model = dblarr(n_elements(finite_index))
  fco2_model = dblarr(n_elements(finite_index))
  cos_profile_saved = dblarr(n_elements(finite_index), n_layers) 
  litter_layers = 0  ; disable litter layers
  
  for i = 0, n_elements(finite_index)-1 do begin
  ;for i = 0, 780 do begin
    ; using different time step for before and after harvesting
    if doy_local[finite_index[i]] lt 130. then dt = .5D else dt = 10d
    ; if doy_local[finite_index[i]] lt 130. then vmax_cos = [5e-7, 0., 8e-11, 0.] else vmax_cos = [5e-7, 0., 2e-11, 0.]
    ; if doy_local[finite_index[i]] lt 130. then vmax_cos = [5e-7, 0., 8e-11, 0.] else vmax_cos = [5e-7, 0., 3e-11, 0.] ; T_eq = 5 C
    ; if doy_local[finite_index[i]] lt 130. then vmax_cos = [1e-7, 0., 1e-10, 0.] else vmax_cos = [1e-7, 0., 3e-11, 0.] ; T_eq = 10 C
    if doy_local[finite_index[i]] lt 130. then vmax_cos = [1.2e-1, 0., 1e-10, 0.] $ ; [3e-7, 0., 1.5e-10, 0.]
      else vmax_cos = [1.2e-1, 0., 3.3e-11, 0.] ; T_eq = 10 C
    
    theta_w = dblarr(1, n_layers)
    ; theta_w[0:n_layers-1] = swc[ finite_index[i] ] + (0.30D - swc[ finite_index[i] ]) * (findgen(n_layers)/(n_layers-1.))^0.5
    ; theta_w = theta_w + swc[ finite_index[i] ] ; constant soil moisture profile for test
    ; use new soil moisture profile
    theta_w = swc_profile[i, *]
    theta_sat = replicate(0.50D, 1, n_layers)
    theta_sat[0:5] = 0.60D  ; assuming larger porosity in the A horizon (top 2 cm)
    theta_a = theta_sat - theta_w
    T_soil = T_soil_profile[i, *]
    ; T_air = T_atm[ finite_index[i] ]
    T_air = T_soil[0]
    co2_atm = co2_a[ finite_index[i] ]  ; in ppmv
    cos_atm = cos_a[ finite_index[i] ]  ; in pptv
    if cos_atm lt 400. or ~finite(cos_atm) then cos_atm = 500.
    if co2_atm lt 300. or co2_atm gt 500. or ~finite(co2_atm) then co2_atm = 400.
    co2_atm = co2_atm * 1e-6 * air_conc * p_atm/p_std * 273.15D / (T_soil[0] + 273.15D)  ; convert to mol m^-3
    cos_atm = cos_atm * 1e-12 * air_conc * p_atm/p_std * 273.15D / (T_soil[0] + 273.15D) ; convert to mol m^-3
    if i eq 0 or abs(fcos_model[i-1]) gt 1e3 then begin  ; abs(fcos_model[i-1]) gt 1e3 indicates an overflow
      co2_soil = replicate(co2_atm, 1, n_layers) * exp(3*z_grid)
      cos_soil = replicate(cos_atm, 1, n_layers) * exp(-0.5*z_grid)  ;* exp(-z_grid/0.4)
    endif else begin
      co2_soil[0] = co2_atm
      cos_soil[0] = cos_atm
    endelse
    if total(~finite(cos_soil)) gt 1 then cos_soil = replicate(cos_atm, 1, n_layers) * exp(-0.5*z_grid)
    f_CA = replicate(1D, 1, n_layers)
    ; if doy_local[finite_index[i]] lt 130. then f_CA[0:10] = 100D
    ; f_CA = f_CA * 0.
    ; hr = (doy_local[finite_index[i]] - floor(doy_local[finite_index[i]])) * 24.
    ; if hr gt 8. and  hr lt 16 then f_CA[0:10] = 1.
    f_biomass = replicate(1D, 1, n_layers)
    lit_wc = 0.
    if i eq 0 then max_time = dt*1e5 else max_time = (doy_utc[finite_index[i]] - doy_utc[finite_index[i-1]]) * 86400D
    if cos_atm lt 400. * 1e-12 * air_conc or ~finite(cos_atm) then begin
      fcos_model[i] = !Values.d_NaN
    endif else begin
      flux_eval = calc_cos_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, cos_atm, cos_soil, f_CA, dt, cv_node, cv_size, n_layers, litter_layers, vmax_cos, max_time)
      fcos_model[i] = flux_eval[0]
    endelse
    
;    ; do not evaluate co2
;    if co2_atm lt 300. * 1e-6 * air_conc or co2_atm gt 500. * 1e-6 * air_conc or ~finite(co2_atm) then begin
;      fco2_model[i] = !Values.d_NaN 
;    endif else begin
;      flux_eval = calc_co2_flux_NSS_FVgrid( T_soil, T_air, theta_sat, theta_w, lit_wc, co2_atm, co2_soil, f_biomass, dt, cv_node, cv_size, n_layers, litter_layers, vmax_co2, max_time)
;      fco2_model[i] = flux_eval[0]
;    endelse
    
    ; plot test figures
    ;if (i eq 0) then flag_plot_test_figures=1 else flag_plot_test_figures=0
    flag_plot_test_figures = 0
    if flag_plot_test_figures then begin   ; need revision here
      stop
      p1 = plot(theta_sat, z_grid, 'ko-', title='SGP', $
        ytitle = 'depth (m)', xtitle = 'Porosity (m$^3$ m$^{-3}$)', yrange=[.5,-0.02], xrange=[0.4,0.65], $
        dimension=[320, 480], margin=[0.24,0.15,0.02,0.1], layout=[2,1,1], font_size=12)
      p1 = plot([0.4, 0.65], [0.,0.], color='light gray', /current, overplot=1)
      
      p2 = plot(theta_w, z_grid, 'ko-',  title='SGP', ytitle = 'depth (m)', $
        xtitle = 'Soil moisture (m$^3$ m$^{-3}$)', xrange=[0.,0.45], yrange=[.5, -0.02], $
        /current, layout=[2,1,2], margin=[0.24,0.15,0.02,0.1], font_size=12)
      p2 = plot([0., 0.45], [0.,0.], color='light gray', /current, overplot=1)
      p2_meas = plot(theta_w[[depth_ind_5cm, depth_ind_30cm]], z_grid[[depth_ind_5cm, depth_ind_30cm]], 'rX', sym_size=1.5, sym_thick=1.5, overplot=1, /current)
      p2.xminor=4
      txt_doy = text(0.65, 0.2, 'DOY =' + string(doy_local[finite_index[i]], format='(f8.3)'), target=p2, font_size=10)
      subfig_lbl1 = text(0.07, 0.93, '(a)', target=p1, font_size=12)
      subfig_lbl2 = text(0.57, 0.93, '(b)', target=p2, font_size=12)
      p2.save, strjoin([plot_dir, '/png/model_input_poros_swc_', string(i, format='(i3.3)'), '.png'], /single), resolution=300 ;border=10
      p2.save, strjoin([plot_dir, 'model_input_poros_swc_', string(i, format='(i3.3)'), '.eps'], /single), resolution=300 ;border=10
      
;      p1=plot([lit_wc[0:litter_layers-1], theta_w[litter_layers:n_layers-1]], z_grid, 'ko-', $
;        ytitle = "depth (m)", xtitle = "water content!C (m$^3$ m$^{-3}$ for soil, g g$^{-1}$ for litter)", xrange=[0,0.35], yrange=[1.0,-0.05], $
;        dimension=[320, 480], margin=[0.2,0.15,0.05,0.1], /buffer)
;      ;p1.title = "doy_utc" + string(doy_utc[ch4_index[i]]) + "!C red: litter water content; cross: measured value at 5cm depth"
;      p1.title = "SC1, day of year" + string(doy_local[ch4_index[i]], format='(f9.3)')
;      p1.xtickinterval = 0.10
;      p1=plot(lit_wc, z_grid[0:litter_layers-1], 'ro-', /current, overplot=1)
;      p1=plot([theta_w[litter_layers]], [z_grid[litter_layers]], 'k+', sym_size = 2.0, /current, overplot=1)
;      p1=plot([0, 0.35], [cv_face[litter_layers-1], cv_face[litter_layers-1] ], '--', color='gray', /current, overplot=1)
;      p1=plot([0, 0.35], [0.,0.], color='light gray', /current, overplot=1)
;      p1.save, strjoin([plot_dir, '/png/model_input_swc_ch4_', string(i, format='(i3.3)'), '.png'], /single), border=10, resolution=300
;      p1.save, strjoin([plot_dir, 'model_input_swc_ch4_', string(i, format='(i3.3)'), '.eps'], /single), border=10, resolution=300
;      
;      ; theta_sat profiles
;      p2=plot(theta_sat, z_grid, 'ko-', title = "SC1, day of year" + string(doy_local[ch4_index[i]], format='(f9.3)'), $
;        ytitle = "depth (m)", xtitle = "porosity (m$^3$ m$^{-3}$)", yrange=[1.0,-0.05], $
;        dimension=[320, 480], margin=[0.2,0.15,0.05,0.1], /buffer)
;      p2=plot(theta_sat[0:litter_layers-1], z_grid[0:litter_layers-1], 'ro-', /current, overplot=1)
;      p2=plot([0.2, 1], [cv_face[litter_layers-1], cv_face[litter_layers-1] ], '--', color='gray', /current, overplot=1)
;      p2=plot([0.2, 1], [0.,0.], color='light gray', /current, overplot=1)
;      p2.save, strjoin([plot_dir, '/png/model_input_porosity_ch4_', string(i, format='(i3.3)'), '.png'], /single), border=10, resolution=300
;      p2.save, strjoin([plot_dir, 'model_input_porosity_ch4_', string(i, format='(i3.3)'), '.eps'], /single), border=10, resolution=300
;      
;      ; cos profile
;      p3=plot(cos_soil / (air_conc*p_atm/p_std*273.15/(T_soil+273.15)) * 1e12 , z_grid, 'ko-', title = "SC1, day of year" + string(doy_local[ch4_index[i]], format='(f9.3)'), $
;        ytitle = "depth (m)", xtitle = "COS concentration (pptv)", yrange=[1.0,-0.05], xrange=[0,600], $
;        dimension=[320, 480], margin=[0.2,0.15,0.05,0.1], /buffer)
;      p3.xminor=4
;      p3=plot(cos_soil[0:litter_layers-1] / (air_conc*p_atm/p_std*273.15/(T_soil[0:litter_layers-1]+273.15)) * 1e12, z_grid[0:litter_layers-1], 'ro-', /current, overplot=1)
;      p3=plot([0, 600], [cv_face[litter_layers-1], cv_face[litter_layers-1] ], '--', color='gray', /current, overplot=1)
;      p3=plot([0, 600], [0.,0.], color='light gray', /current, overplot=1)
;      p3=plot([cos_atm, cos_soil[0]] /(air_conc*p_atm/p_std*273.15/(T_air+273.15)) * 1e12, [0., z_grid[0]], 'k-', /current, overplot=1)
;      p3=plot([cos_atm /(air_conc*p_atm/p_std*273.15/(T_air+273.15)) * 1e12 ], [0.], 'bo', /current, overplot=1)  ; atms conc
;      p3.save, strjoin([plot_dir, '/png/model_cos_profile_ch4_', string(i, format='(i3.3)'), '.png'], /single), border=10, resolution=300  
;      p3.save, strjoin([plot_dir, 'model_cos_profile_ch4_', string(i, format='(i3.3)'), '.eps'], /single), border=10, resolution=300  
      STOP
    endif
    cos_profile_saved[i, *] = cos_soil / (air_conc*p_atm/p_std*273.15/(T_soil+273.15)) * 1e12   ; save calculated cos profile (in ppt) for later use
    ;if i mod 20 eq 0 then print, '*', format='(a, $)'
    print, i, fcos_model[i], fcos_model[i]-fcos[finite_index[i]]
    ; print, i, fco2_model[i], fco2_model[i]-fco2[finite_index[i]]
    ;p = plot(cos_soil / air_conc * 1e12, -z_grid, yrange=[-1,0],'k+-', sym_color='red')
    ;STOP
  endfor
  ; stop
  
  ; calculating diagnostics
  mean_obs_fcos = mean(fcos[finite_index])
  mean_mdl_fcos = mean(fcos_model, /NaN)
  rmse_fcos = sqrt(mean((fcos_model - fcos[finite_index])^2, /NaN))
  ; rmse_fco2 = sqrt(mean((fco2_model - fco2[finite_index])^2, /NaN))
  fcos_finite_index = where(finite(fcos[finite_index]) and finite(fcos_model))
  ; fco2_finite_index = where(finite(fco2[finite_index]) and finite(fco2_model))
  r_fcos = correlate(fcos_model[fcos_finite_index], fcos[finite_index[fcos_finite_index]])
  ; r_fco2 = correlate(fco2_model[fco2_finite_index], fco2[finite_index[fco2_finite_index]])
  resid_fcos = fcos[finite_index] - fcos_model
  
  stop
  
  ; bin the data with respect to soil temp and swc
  n_bins1 = 9
  n_bins2 = 8
  bins_array = intarr(n_bins1*n_bins2, 200) - 1
  n_binned = intarr(n_bins1*n_bins2)
  for i=0,8 do begin
    for j=0,7 do begin
      current_bin = where(T_s[finite_index] ge 5.*(i+1) and T_s[finite_index] lt 5.*(i+2) $
        and swc[finite_index] ge 0.05*(j+1) and swc[finite_index] lt 0.05*(j+2))
      if total(current_bin) ne -1 then n_binned[i*n_bins2+j] = n_elements(current_bin) else n_binned[i*n_bins2+j] = 0
      if n_binned[i*n_bins2+j] gt 0 then bins_array[i*n_bins2+j, 0:(n_binned[i*n_bins2+j]-1)] = current_bin
    endfor
  endfor
  fcos_obs_binned = dblarr(n_bins1 * n_bins2) * !Values.d_NaN
  fcos_mdl_binned = dblarr(n_bins1 * n_bins2) * !Values.d_NaN
  for i=0,8 do $
    for j=0,7 do $
      if n_binned[i*n_bins2+j] gt 0 then fcos_obs_binned[i*n_bins2+j] = mean(fcos[finite_index[ bins_array[i*n_bins2+j, 0:(n_binned[i*n_bins2+j]-1)] ]], /NaN)
  
  for i=0,8 do $
    for j=0,7 do $
      if n_binned[i*n_bins2+j] gt 0 then fcos_mdl_binned[i*n_bins2+j] = mean(fcos_model[ bins_array[i*n_bins2+j, 0:(n_binned[i*n_bins2+j]-1)] ], /NaN)
  
  fcos_resid_binned = fcos_obs_binned - fcos_mdl_binned
  n_resid_levels = floor(max(fcos_resid_binned, /NaN)*2.) - floor(min(fcos_resid_binned, /NaN)*2.) + 1
  ; load color tables for the binned plot
  loadct, 33, rgb_table = temp_ct   ; for temperature
  rainbow_ct = bytarr(n_resid_levels, 3) 
  for i=0, n_resid_levels-1 do rainbow_ct[i, *] = temp_ct[i*11+24, *] ; temp_ct[i*16+24, *]
  
  ; a test figure
;  p=plot([!Values.f_NaN], [!Values.f_NaN], xrange=[5,50], yrange=[0.05,0.45], xtitle='Soil temp ($\deg$C)', ytitle='Soil moisture (m$^3$ m$^{-3}$)')
;  for i=0,8 do $
;    for j=0,7 do $
;      if finite(fcos_resid_binned[i*n_bins2+j]) then $
;        poly=polygon([5.*(i+1),5.*(i+2),5.*(i+2),5.*(i+1)], [0.05*(j+1),0.05*(j+1),0.05*(j+2),0.05*(j+2)], /current, /data, /fill_background, $
;          fill_color= transpose(rainbow_ct[ floor(fcos_resid_binned[i*n_bins2+j]*2.) - floor(min(fcos_resid_binned, /NaN)*2.), * ]) )
;  
  
  ; getting flux predictions with other methods for comparison
  ; Berry et al. (2013) function
  ; not applicable to this case because of strong emissions
  ; wfps = swc[finite_index] / 0.45
  ; fcos_b13 = -fco2[finite_index]
  ; Linear SRU-T function, Berkelhammer et al. (2014) modified
  sru = (fcos[finite_index] / fco2[finite_index]) / (cos_a[finite_index]/co2_a[finite_index])
  sru[where(sru lt -3 or sru gt 10)] = !Values.d_NaN
  x_fit = T_s[finite_index]
  y_fit = sru
  xy_finite = where(finite(x_fit) and finite(y_fit))
  x_fit = x_fit[xy_finite]
  y_fit = y_fit[xy_finite]
  fit_sru_t = ladfit(x_fit, y_fit)
  fcos_sru = (fit_sru_t[0] + fit_sru_t[1] * T_s[finite_index]) * (cos_a[finite_index]/co2_a[finite_index]) * fco2[finite_index]
  ; diagnostics for SRU model
  sru_mdl_diag = create_struct('mean_mdl_fcos', !Values.d_NaN, $
    'rmse_fcos', !Values.d_NaN, $
    'fcos_finite_index', where(finite(fcos[finite_index]) and finite(fcos_sru)), $
    'r_fcos', !Values.d_NaN, $
    'resid_fcos', fcos[finite_index] - fcos_sru)
  sru_mdl_diag.mean_mdl_fcos = mean(fcos_sru, /NaN)
  sru_mdl_diag.rmse_fcos = sqrt(mean((fcos_sru - fcos[finite_index])^2, /NaN))
  sru_mdl_diag.r_fcos = correlate(fcos_sru[sru_mdl_diag.fcos_finite_index], fcos[finite_index[sru_mdl_diag.fcos_finite_index]])
  
  stop
  
  ; test figure to show observed results and model results
  ; not saved
  p_obs = plot(doy_local[finite_index], fcos[finite_index], 'ko', sym_size=0.5, xrange=[90,155], yrange=[-10,30], xtitle='Day of Year', $
    ytitle='$f_{COS}$ ($pmol m^{-2} s^{-1}$)', title='SGP')
  p_obs.color = 'gray'
  p_obs.sym_color = 'black'
  ;  title='black: obs, red: model', dimensions=[800,600])
  p_mdl = plot(doy_local[finite_index], fcos_model, 'r-', thick=2, /current, overplot=1)
  p_mdl_sru = plot(doy_local[finite_index], fcos_sru, color='light salmon', thick=1, /current, overplot=1)
  p_zero = plot([90,155], [0,0], '__', color='gray', thick=1.5, /current, overplot=1)
  p_zero.order, /send_to_back
  p_obs.font_size=18  ; for AGU poster
;  ; plot inset figure
;  p_fit = plot([-6.5,4], [-6.5,4], xrange=[-6.5,4], yrange=[-6.5,4], 'r--', thick=1.5, /current, position=[0.4, 0.25, 0.65, 0.5], xtitle='observations', ytitle='model')
;  p_fit = plot(f_cos[ch4_index], fcos_model_ch4, 'ko', /current, overplot=1)
;  p_fit.xtickinterval = 2
;  p_fit.ytickinterval = 2
;  p_fit.sym_size=0.3
;  p_fit.sym_filled=1
;  subfig_lbl1= text(0.85, 0.80, '(a)', font_size=16)
  leg1 = legend(position=[110,25], target=[p_obs, p_mdl], font_size=12, sample_width=0.1, /data, /auto_text_color, shadow=0)
  leg1[0].label = ' observations'
  leg1[1].label = ' simulated flux' 
  leg1.transparency=100
  
  ; figure in the paper
  flag = intarr(1)
  read, flag, prompt='Plot model results and soil variables? (1 or 0) '
  if flag then begin
  p_combined = objarr(4,5)
  p_combined[0,0] = plot(doy_local[finite_index], fcos[finite_index], 'ko-', sym_size=0.5, xrange=[90,155], yrange=[-10,25], xtitle='Day of Year', $
    ytitle='$F_{COS}$ ($pmol m^{-2} s^{-1}$)', title='Southern Great Plains, OK', layout=[1,2,1], dimensions=[800,450], margin=[0.07,0.14,0.4,0.1])
  ; figure aspect ratio: 16 to 9
  p_combined[0,1] = plot(doy_local[finite_index], fcos_model, 'r-', thick=2, /current, overplot=1)
  ; p_combined[0,2] = plot(doy_local[finite_index], fcos_sru, color='cornflower', thick=1, /current, overplot=1)
  ; leg1 = legend(position=[110,23], target=[p_combined[0,0:3]], font_size=10, sample_width=0.05, /data, /auto_text_color, shadow=0)
  leg1 = legend(position=[110,23], target=[p_combined[0,0:2]], font_size=10, sample_width=0.05, /data, /auto_text_color, shadow=0)
  leg1[0].label = ' observations'
  leg1[1].label = ' simulated flux' 
  ;leg1[2].label = ' calculated from SRU-$T$'
  leg1.transparency=100
  
  p_combined[1,0] = plot([-10,25], [-10,25], xrange=[-10,25], yrange=[-10,25], 'r--', $
    thick=1.5, layout=[1,2,1], /current, position=[0.7, 0.64, 0.9, 0.94], xtitle='Observed $F_{COS}$ (pmol m$^{-2}$ s$^{-1}$)', ytitle='Modeled $F_{COS}$ (pmol m$^{-2}$ s$^{-1}$)', font_size=10)
  p_combined[1,1] = plot(fcos[finite_index], fcos_model, 'ko', sym_size=0.3, sym_filled=1, /current, overplot=1)
  t_r2 = text(7.5, -5, '$r^2=$' + string(r_fcos^2, format='(F5.3)'), /data, target=p_combined[1,1],font_size=10)
  ax = p_combined[1,0].axes
  ax[3].showtext=1
  ax[1].showtext=0
  
  
  seg1 = where(doy_local[finite_index] lt 144.)
  seg2 = where(doy_local[finite_index] gt 144.)
  p_combined[2,0] = plot(doy_local[finite_index[seg1]], T_s[finite_index[seg1]], color='orange red', xrange=[90,155], yrange=[0,50], xtitle='Day of Year', $
    ytitle='Soil temperature ($\degC$)', layout=[1,2,2], /current, margin=[0.07,0.18,0.4,0.06])
  p_combined[2,1] = plot(doy_local[finite_index[seg2]], T_s[finite_index[seg2]], color='orange red', layout=[1,2,2], /current, overplot=1)    
  p_combined[2,2] = plot(doy_local[finite_index[seg1]], swc[finite_index[seg1]], color='deep sky blue', xrange=[90,155], yrange=[0,0.5], xtitle='Day of Year', $
    ytitle='Soil moisture ($m^3 m^{-3}$)', layout=[1,2,2], /current, margin=[0.07,0.18,0.4,0.06])
  p_combined[2,3] = plot(doy_local[finite_index[seg2]], swc[finite_index[seg2]], color='deep sky blue', layout=[1,2,2], /current, overplot=1)
  txt_Ts = text(118, 10, '$T_{soil}$', font_size=10, color='orange red', /data, target=p_combined[2,0])
  txt_swc = text(109, 0.4, 'SWC', font_size=10, color='deep sky blue', /data, target=p_combined[2,2])
  ax = p_combined[2,0].axes
  ax[0].hide=1
  ax[2].hide=1
  ax[3].hide=1
  ax = p_combined[2,2].axes
  ax[1].hide=1
  ax[3].showtext=1
  
  p_combined[3,0] = plot([!Values.f_NaN], [!Values.f_NaN], xrange=[5,50], yrange=[0.05,0.45], xtitle='Soil temp ($\deg$C)', ytitle='Soil moisture (m$^3$ m$^{-3}$)', $
    layout=[1,2,2], /current, position=[0.7, 0.21, 0.9, 0.51], font_size=10)
  for i=0,8 do $
    for j=0,7 do $
      if finite(fcos_resid_binned[i*n_bins2+j]) then $
        poly=polygon([5.*(i+1),5.*(i+2),5.*(i+2),5.*(i+1)], [0.05*(j+1),0.05*(j+1),0.05*(j+2),0.05*(j+2)], target=p_combined[3,0], /current, /data, /fill_background, $
          fill_color= transpose(rainbow_ct[ floor(fcos_resid_binned[i*n_bins2+j]*2.) - floor(min(fcos_resid_binned, /NaN)*2.), * ]) )
  p_combined[3,0].xminor=1
  p_combined[3,0].yminor=1
  ax = p_combined[3,0].axes
  ax[3].showtext=1
  ax[1].showtext=0
  ; draw colorbar below the figure
  for i=0, n_resid_levels-1 do $ 
    poly_cb = polygon([0.67-0.012+0.012*i,0.685-0.012+0.012*i,0.685-0.012+0.012*i,0.67-0.012+0.012*i], [0.1,0.1,0.12,0.12], target=p_combined[3,0], /current, /fill_background, $
      fill_color=transpose(rainbow_ct[i, *]) ) 
  for i=0, floor((n_resid_levels)/2) do $
    txt_cb = text(0.66-0.012+0.024*i, 0.07, string(-4+i,format='(i2)'), target=p_combined[3,0], /current, font_size=10)
  ; txt_cb_unit = text(0.885, 0.07, 'pmol m$^{-2}$ s$^{-1}$', target=p_combined[3,0], /current, font_size=10)
  txt_cb_unit_ln1 = text(0.91+0.012, 0.10, 'pmol', target=p_combined[3,0], /current, font_size=10)
  txt_cb_unit_ln2 = text(0.91+0.012, 0.065, 'm$^{-2}$ s$^{-1}$', target=p_combined[3,0], /current, font_size=10)
  txt_cb_lgd = text(0.79+0.012, 0.03, 'Mean residual in the binned domain', font_size=10, alignment=0.5)
  
  subfig_lbl1 = text(0.09,0.89,'(a)', font_size=12)
  subfig_lbl2 = text(0.09,0.41,'(b)', font_size=12)
  subfig_lbl3 = text(0.66,0.92,'(c)', font_size=12)
  subfig_lbl4 = text(0.66,0.47,'(d)', font_size=12)
  
  p_combined[3,0].save, strjoin([plot_dir, 'fcos_model', '.eps'], /single), border=10, resolution=300
  p_combined[3,0].save, strjoin([plot_dir, 'png/fcos_model', '.png'], /single), border=10, resolution=300
  endif 
  stop
  
  ; do not go any further
   
  
  ;p_obs.save, strjoin([plot_dir, 'fcos_model_ch4_fast', '.eps'], /single), border=10, resolution=300
  ;p_obs.save, strjoin([plot_dir, 'png/fcos_model_ch4_fast', '.png'], /single), border=10, resolution=300
  
;  p = plot(doy_local[ch4_index], f_co2[ch4_index], 'ko-', sym_size=0.65, xrange=[90,135], yrange=[0,7], xtitle='Day of Year', $
;    ytitle='$f_{CO2}$ ($\mumol m^{-2} s^{-1}$)', title='Soil Chamber 1')
;  p.ytickinterval=1
;  p.yminor=4
;  p.color = 'gray'
;  p.sym_color = 'black'
;  ;  title='black: obs, red: model', dimensions=[800,600])
;  p = plot(doy_local[ch4_index], fco2_model_ch4, 'r-', thick=2, /current, overplot=1)
;  p = plot(doy_local[ch4_index], fco2_model_ch4_interface, 'b-', thick=1, /current, overplot=1)
;  p.font_size=18  ; for AGU poster
;  ; plot inset figure
;  p_fit = plot([0,5], [0,5], xrange=[0,5], yrange=[0,5], 'r--', thick=1.5, /current, position=[0.35, 0.5, 0.6, 0.75], xtitle='observations', ytitle='model')
;  p_fit = plot(f_co2[ch4_index], fco2_model_ch4, 'ko', /current, overplot=1)
;  p_fit.xtickinterval = 1
;  p_fit.ytickinterval = 1
;  p_fit.xminor=4
;  p_fit.yminor=4
;  p_fit.sym_size=0.3
;  p_fit.sym_filled=1
;  subfig_lbl3= text(0.85, 0.80, '(c)', font_size=16)
;  ; p.save, strjoin([plot_dir, 'fco2_model_ch4_fast', '.eps'], /single), border=10, resolution=300
;  ; p.save, strjoin([plot_dir, 'png/fco2_model_ch4_fast', '.png'], /single), border=10, resolution=300
;  
;  p = plot(Ts_1[ch4_index], fcos_model_ch4 / fco2_model_ch4 * co2_a[ch4_index] / cos_a[ch4_index], $
;    yrange=[-15,5], xtitle='Soil Temperature ($\deg$C) ', ytitle='Modeled soil relative uptake', title='Soil Chamber 1', $
;    'ko', sym_size=0.65, /buffer)
;  p.save, strjoin([plot_dir, 'sru_model_ch4_fast', '.eps'], /single), border=10, resolution=300
;  p.save, strjoin([plot_dir, 'png/sru_model_ch4_fast', '.png'], /single), border=10, resolution=300
;  
  plot_cos_profile_flag = 1
  if plot_cos_profile_flag then begin
    conc_i1=102 ;598; 692
    conc_i2=1174
    p_conc1=plot(cos_profile_saved[conc_i1,*] , z_grid, 'ko-', title = string(replicate(32B,40)) + 'Southern Great Plains, OK', $
      ytitle = 'depth (m)', xtitle = '', yrange=[1.0,-0.05], xrange=[0,500], $
      dimensions=[270,480], margin=[0.24,0.15,0.0,0.1], layout=[2,1,1], font_size=12)
    p_conc1.xminor = 4
    p_conc3=plot([0, 500], [0.,0.], color='light gray', /current, overplot=1)    
    p_conc1atmline = plot([cos_a[finite_index[conc_i1]], cos_profile_saved[conc_i1,0]], [0., z_grid[0]], 'k-', /current, overplot=1)
    p_conc1atm = plot([cos_a[finite_index[conc_i1]]], [0.], 'bo', /current, overplot=1)
    
    p_conc2=plot(cos_profile_saved[conc_i2,*] , z_grid, 'ks-', /current, $
      yrange=[1.0,-0.05], xrange=[400,20000], margin=[0.06,0.15,0.18,0.1], layout=[2,1,2], /xlog, font_size=12)
    p_conc4=plot([400, 400000], [0.,0.], color='light gray', /current, overplot=1)  
    p_conc2atmline = plot([cos_a[finite_index[conc_i2]], cos_profile_saved[conc_i2,0]], [0., z_grid[0]], 'k-', /current, overplot=1)
    p_conc2atm = plot([cos_a[finite_index[conc_i2]]], [0.], 'bs', /current, overplot=1)
    ax = p_conc1.axes
    ax[3].hide=1
    ax = p_conc2.axes
    ax[1].hide=1
    ax[3].showtext=1
    x_txt = text(0.5, 0.05, 'COS concentration (pptv)', target=p_conc2, alignment=0.5)
    
    lgd1_ln1 = text(0.35, 0.6, 'DOY =' + string(doy_local[finite_index[conc_i1]], format='(f9.2)'), layout=[2,1,1], font_size=10)
    lgd1_ln2 = text(0.35, 0.55, 'Modeled COS flux =', layout=[2,1,1], font_size=10)
    lgd1_ln3 = text(0.35, 0.51, string(fcos_model[conc_i1], format='(f7.3)') + ' pmol m$^{-2}$ s$^{-1}$', layout=[2,1,1], font_size=10)
    
    lgd2_ln1 = text(0.58, 0.40, 'DOY =' + string(doy_local[finite_index[conc_i2]], format='(f9.2)'), layout=[2,1,2], font_size=10)
    lgd2_ln2 = text(0.58, 0.35, 'Modeled COS flux =', layout=[2,1,2], font_size=10)
    lgd2_ln3 = text(0.58, 0.31, string(fcos_model[conc_i2], format='(f7.3)') + ' pmol m$^{-2}$ s$^{-1}$', layout=[2,1,2], font_size=10)
    
    subfig_lbl = text(0.07, 0.93, '(a)', font_size=14)
    
;    leg_pconc=legend(target=[p_conc1, p_conc2], position=[0.5,0.3], shadow=0, font_size=9)
;    leg_pconc[0].label = '  DOY =' + string(doy_local[finite_index[conc_i1]], format='(f9.3)')
;    leg_pconc[1].label = '  DOY =' + string(doy_local[finite_index[conc_i2]], format='(f9.3)')
    p_conc1.save, strjoin([plot_dir, 'cos_profile_comp', '.eps'], /single), border=10, resolution=300
    p_conc1.save, strjoin([plot_dir, 'png/cos_profile_comp', '.png'], /single), border=10, resolution=300
  endif
  
  
  ; surface_mean_fcos_ch4 = int_tabulated(doy_utc[ch4_index[fcos_finite_index]], fcos_model_ch4[fcos_finite_index], /double) $
  ;   / (doy_utc[ch4_index[-1]] - doy_utc[ch4_index[0]]) 
  ; soil_mean_fcos_ch4 = int_tabulated(doy_utc[ch4_index[fcos_finite_index]], fcos_model_ch4_interface[fcos_finite_index], /double) $
  ;   / (doy_utc[ch4_index[-1]] - doy_utc[ch4_index[0]]) 
;  surface_mean_fcos_ch4 = mean(fcos_model_ch4, /NaN)
;  soil_mean_fcos_ch4 = mean(fcos_model_ch4_interface, /NaN)
;  litter_mean_fcos_ch4 = surface_mean_fcos_ch4 - soil_mean_fcos_ch4
;  mean_soil_contrib_ch4 = soil_mean_fcos_ch4 / surface_mean_fcos_ch4
;  soil_contrib_ch4 = fcos_model_ch4_interface / fcos_model_ch4
;  ; first5days_index = where( doy_utc[ch4_index[fcos_finite_index]] lt 96.)
;  first5days_index = where( doy_utc[ch4_index] lt 96.)
;  ; soil_contrib_ch4_first5days = int_tabulated(doy_utc[ch4_index[fcos_finite_index[first5days_index]]], fcos_model_ch4_interface[fcos_finite_index[first5days_index]], /double) / $
;  ;   int_tabulated(doy_utc[ch4_index[fcos_finite_index[first5days_index]]], fcos_model_ch4[fcos_finite_index[first5days_index]], /double)
;  soil_contrib_ch4_first5days = mean(fcos_model_ch4_interface[first5days_index], /NaN) / mean(fcos_model_ch4[first5days_index], /NaN)
;  
;  surface_mean_fco2_ch4 = mean(fco2_model_ch4, /NaN)
;  soil_mean_fco2_ch4 = mean(fco2_model_ch4_interface, /NaN)
;  litter_mean_fco2_ch4 = surface_mean_fco2_ch4 - soil_mean_fco2_ch4
;  mean_soil_contrib_fco2_ch4 = soil_mean_fco2_ch4 / surface_mean_fco2_ch4
;  soil_contrib_fco2_ch4 = fco2_model_ch4_interface / fco2_model_ch4
;  soil_contrib_fco2_ch4_first5days = mean(fco2_model_ch4_interface[first5days_index], /NaN) / mean(fco2_model_ch4[first5days_index], /NaN)
;  
;  ; plot litter vs soil contribution to uptake or emissions
;  ; bin the fluxes by day
;  
;  ;  normlz_factor_ch4 = max( abs( [[fcos_model_ch4],[fcos_model_ch4_interface],[fcos_model_ch4-fcos_model_ch4_interface]] ), dimension=2)
;  ;  litter_frac_ch4 = (fcos_model_ch4-fcos_model_ch4_interface)/normlz_factor_ch4 * 100.
;  ;  soil_frac_ch4 = fcos_model_ch4_interface/normlz_factor_ch4 * 100.
;  n_days_ch4 = floor(max(doy_local[ch4_index])) - ceil(min(doy_local[ch4_index])) + 1
;  start_doy_ch4 = ceil(min(doy_local[ch4_index]))
;  fcos_model_ch4_daily_mean = dblarr(n_days_ch4)
;  fcos_model_ch4_interface_daily_mean = dblarr(n_days_ch4)
;  
;  for day_no = 0, n_days_ch4-1 do begin
;    selected_doy_index = where(doy_local[ch4_index] gt start_doy_ch4 + day_no and doy_local[ch4_index] lt start_doy_ch4 + day_no + 1)
;    fcos_model_ch4_daily_mean[day_no] = mean(fcos_model_ch4[selected_doy_index])
;    fcos_model_ch4_interface_daily_mean[day_no] = mean(fcos_model_ch4_interface[selected_doy_index])
;  endfor
;  
;  normlz_fac_ch4_daily_mean = max( abs( [[fcos_model_ch4_daily_mean], [fcos_model_ch4_interface_daily_mean], $
;    [fcos_model_ch4_daily_mean-fcos_model_ch4_interface_daily_mean]] ), dimension=2)
;  litter_fcos_frac_ch4 = (fcos_model_ch4_daily_mean - fcos_model_ch4_interface_daily_mean) / normlz_fac_ch4_daily_mean * 100.
;  soil_fcos_frac_ch4 = fcos_model_ch4_interface_daily_mean / normlz_fac_ch4_daily_mean * 100.
;  same_sign_index_ch4 = where(soil_fcos_frac_ch4 * litter_fcos_frac_ch4 ge 0.)
;  diff_sign_index_ch4 = where(soil_fcos_frac_ch4 * litter_fcos_frac_ch4 lt 0.)
;  
;  p1 = barplot(findgen(41)+91+0.5, litter_fcos_frac_ch4, yrange=[-105,105], fill_color='sky blue', $
;    xtitle='Day of year', ytitle='Uptake (%) $\leftarrow$----   ----$\rightarrow$ Emissions (%)', $
;    title='Soil Chamber 1')
;  p2 = barplot(same_sign_index_ch4+91+0.5, litter_fcos_frac_ch4[same_sign_index_ch4] + soil_fcos_frac_ch4[same_sign_index_ch4], $
;    bottom_values=litter_fcos_frac_ch4[same_sign_index_ch4], fill_color='sandy brown', /overplot)
;  p3 = barplot(diff_sign_index_ch4+91+0.5, soil_fcos_frac_ch4[diff_sign_index_ch4], fill_color='sandy brown', /overplot)
;  
;  s1=symbol(0.2, 0.78, 's', sym_color='black', sym_fill_color='sky blue', /sym_filled, sym_size=1.5)
;  s2=symbol(0.2, 0.70, 's', sym_color='black', sym_fill_color='sandy brown', /sym_filled, sym_size=1.5)
;  t1=text(0.23, 0.75, 'Litter contribution!Cfraction to COS flux')
;  t2=text(0.23, 0.67, 'Soil contribution!Cfraction to COS flux')
;  
;  p1.save, strjoin([plot_dir, 'litter_contrib_ch4', '.eps'], /single), border=10, resolution=300
;  p1.save, strjoin([plot_dir, '/png/litter_contrib_ch4', '.png'], /single), border=10, resolution=300
;  
  ; plot residuals
  ; resd_fcos_ch4 = fcos_model_ch4 - f_cos[ch4_index]
  ; p=plot(doy_local[ch4_index], resd_fcos_ch4, 'ko')
  endif
  stop
  
END