# Metadata information for the soil COS flux model project

Wu Sun, 01/05/2016

Modified on 2018-07-11

## General description

This repository contains all the data, figures, code and relevant supporting
information for my authored publication: Sun, W. et al. (2015). Geosci. Model
Dev., 8, 3055–3070. The paper is openly accessible on the journal website. 

## Data files

Data used in this study can be found in the `/data` directory. Description of
each data file is as follows:

- `met_station_soil_data.csv`: Soil temperature and moisture profiles at the
  site.
- `soil_flux_regular_meas.csv`: The part of soil flux data that are relevant to
  my modeling study. The column names of the table should be self-explanatory.
- `soil_moist_temp.csv`: Soil surface moisture and temperature measured at the
  site. *Not used*.
- `stunt_ranch_test_data.csv`: Extracted test dataset from the Stunt Ranch flux
  data. Only the first 15 days of data among the full dataset were used for
  modeling study. For the full dataset, see Sun, W. et al. (2016). J. Geophys.
  Res. Biogeosci.
- `T_soil_gap_filled.csv`: Gap-filled Stunt Ranch soil temperature time series.
- `temp_soil_comm_1.txt`, `Tsoil_SWC_ChTair.xlsx`: Partial record for soil
  temperature and moisture at the Oklahoma site. *Not used* in this modeling
  study.

## Source code

See `/src` folder for the source code. 

- `fit_fcos_vs_temp.R`: Fit a temperature function for COS emission at the
  Oklahoma site. Not used.
- `fnpsvp.pro`: IDL routine inherited from Ulli Seibt, which is used to
  calculate saturation vapor pressure of water.
- `mdl_flux_fvgrid.pro`: Main program (in IDL) for calculating soil COS flux at
  Southern Great Plains, Oklahoma
- `mdl_SR_flux_fvgrid.pro`: Main program (in IDL) for calculating soil COS flux
  at Stunt Ranch, with the test data (15 days, full dataset 40 days).
- `plot_gas2aq_flux_ratio.py`: Plot the ratio of gaseous diffusive flux to
  aqueous diffusive flux under different soil water saturation status. The
  plotted figure is used in the supplement.
- `plot_vert_grid.py`: Plot the vertical grid of the model.
- `sens_test.pro`: Sensitivity tests of soil COS flux against soil moisture and
  Vmax
- `sgp_plot_soil_fluxes.py`: To reproduce some of the figures for the Oklahoma
  site data with Python, for AGU 2015 Fall Meeting presentation. But this
  script was later abandoned and hence not used. 
- `soil_therm_cond.R`: Calculate soil thermal conductivity from soil moisture.
- `solver_tridiag_fvgrid.pro`: Tridiagonal matrix solver for the finite-volume grid

## Figures

Figures are in the `plots` folder. For description of the figures, check the
paper itself or the code that produce the figures.

## Other supporting data

- Supplements (not used): `supplements_for_paper` folder contains data and
  source code that were intended to be released with the paper as the
  supplement. But eventually, they were not included in the supplement. 
- Additional tests for potential biases: `test_of_biases` folder contains
  back-of-the-envelope tests of photochemical production of COS [Whelan and
  Rhew (2015). J. Geophys. Res. Biogeosci.], turnover time of soil COS uptake,
  and potential bias in fluxes that may be caused by water vapor concentration
  changes. 

## Permissions

All source code files in this repository are licensed under the permissive [MIT
License](https://en.wikipedia.org/wiki/MIT_License). However, I (Wu Sun) do not
personally own the *data* used to evaluate the model in this study. Anyone who
wishes to use the data presented in this study must obtain permissions from the
correspondent authors of the following studies.

* Oklahoma ARM site soil COS flux data: Maseyk, K. et al. (2014). Proc. Natl.
  Acad. Sci., 111(25), 9064–9069.
* Stunt Ranch site soil COS flux data: Sun, W. et al. (2016). J. Geophys. Res.
  Biogeosci., 121(2), 438–450.
