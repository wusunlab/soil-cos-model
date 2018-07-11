# Soil carbonyl sulfide flux model

Wu Sun (wu.sun@ucla.edu)

This repository stores legacy code of a model that simulates soil-atmosphere
flux of carbonyl sulfide (COS) — an exotic trace gas that is of interest to
sulfur cycling and ecosystem photosynthesis studies. The model solves a 1D
diffusion–(advection)–reaction equation of COS in the soil column. Details of
the model structure is described in the paper Sun et al. (2015).

There are two implementation of the model, one in IDL (a proprietary high level
language), the other in Python 3. The Python 3 version is the more updated
version and to use it requires `numpy`, `scipy`, and `matplotlib` to be
pre-installed.

**Warning**: The coding style of the Python version does not adhere to the PEP8
standard. Moreover, it does not separate script from functions. Some scripts
lack documentation and have not been sufficiently tested. Use with caution.

## Usage

Since the model is provided as scripts not packages, you can just checkout this
git repository

```shell
git checkout https://github.com/wusunlab/soil-cos-model.git
```

The IDL version is in `src/IDL`, whereas the Python version is in `src/python`.

## License

[MIT License](LICENSE)

## Documentation

There is no documentation (like most academic code), alas. Please consult the
publication [Sun et al. (2015)](https://www.geosci-model-dev.net/8/3055/2015/).

## Issues

Development has been discontinued as of 2016. I would not be able to provide
further maintenance or support. However, a new framework is currently in
development (one soil flux model that rules 'em all).

## Citation

If the model is of use to your research and does not crash your computer,
consider citing our work:

Sun, W., Maseyk, K., Lett, C., and Seibt, U.: A soil diffusion–reaction model
for surface COS flux: COSSM v1, Geosci. Model Dev., 8, 3055-3070,
<https://doi.org/10.5194/gmd-8-3055-2015>, 2015.
