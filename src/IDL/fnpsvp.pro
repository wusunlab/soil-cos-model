;; $Id: fnpsvp.pro,v 1.2 2004/01/15 15:54:29 bhaxo Exp $
;;
;; PURPOSE - Groff-Gratch equation for saturation vapor pressure as
;;           a function of temperature in Celcius (default) or Kelvin
;;
;; INPUTS: temperature, degrees Celsius (celsius) or Kelvin (set keyword)
;; OUTPUT: vapor pressure, Pascals
;;
;;----------------------------------------------------------------------
FUNCTION fnpsvp,t,Kelvin=kelvin

  compile_opt idl2

  IF n_params() NE 1 THEN $
    Message,'ERROR: Number of parameters incorrect.'

  ;; Scalars passed by reference. Assign value to local variable
  IF keyword_set(kelvin) THEN tkelvin = t $
                         ELSE tkelvin = t + 273.16
  u = tkelvin / 373.16
  v = 373.16 / tkelvin

  tmp = -7.90298 * (v - 1.0) + 5.02808 * alog10(v) $
        - 1.3816e-7 * (10.0^(11.344 * (1.0 - u)) - 1.0) $
        + 8.1328e-3 * (10.0^(-3.49149 * (v - 1.0)) - 1.0) $
        + alog10(1013.246)

  esat  = 100.0 * (10.0^tmp)

  return, esat                  ; Saturation vapor pressure, Pascals
;;----------------------------------------------------------------------
END
