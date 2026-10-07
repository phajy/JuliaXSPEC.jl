# Reference

```@meta
CurrentModule = JuliaXSPEC
```

```@docs
JuliaXSPEC
```

## Spectra

```@docs
BinnedSpectrum
nbins
bin_widths
bin_centres
density
from_density
rebin
```

## Models and parameters

```@docs
Parameter
model_dat_line
AbstractXSPECModel
AdditiveModel
MultiplicativeModel
ConvolutionModel
xspec_type
parameter_names
parameter_defaults
named_parameters
evaluate
```

## Registry and generated files

```@docs
register!
registered_models
lookup_model
c_function_name
model_dat_text
c_wrappers_text
write_xspec_package
```

## Caching on a grid

```@docs
GridAxis
detect_scale
GridInterpolator
grid_size
evaluate_exact
cache_stats
empty_cache!
```

## Kernels and convolution

```@docs
RedshiftKernel
GaussianKernel
cumulative
convolution_matrix
convolve
```

## Analytic spectra

```@docs
gaussian_per_bin
powerlaw_per_bin
```

## Entry point and options

```@docs
juliaxspec_evaluate
verbose
```
