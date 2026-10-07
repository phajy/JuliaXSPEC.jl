"""
    JuliaXSPEC

Write XSPEC spectral models in Julia.

A model is a plain Julia function wrapped in an [`AdditiveModel`](@ref),
[`MultiplicativeModel`](@ref) or [`ConvolutionModel`](@ref) and made known
with [`register!`](@ref). From the registered models JuliaXSPEC generates the
`model.dat` and C wrapper files XSPEC needs, and a single C-callable entry
point ([`juliaxspec_evaluate`](@ref)) forwards every XSPEC call to the right
Julia function.

Two further building blocks are provided for expensive models:

- [`GridInterpolator`](@ref) caches a function on a parameter grid and
  interpolates between the grid points (a table model built on the fly).
- [`convolve`](@ref) blurs a [`BinnedSpectrum`](@ref) with a
  [`RedshiftKernel`](@ref) in `g = E_obs / E_em`, with the bin-integrated
  bookkeeping done explicitly. [`convolve_fft`](@ref) evaluates the same
  operator with a fast Fourier transform.
- [`OGIPTable`](@ref) reads an XSPEC table model, and [`Blurred`](@ref) folds
  one together with a kernel.
"""
module JuliaXSPEC

include("spectrum.jl")
include("parameters.jl")
include("models.jl")
include("registry.jl")
include("options.jl")
include("grid.jl")
include("kernels.jl")
include("convolution.jl")
include("convolution_fft.jl")
include("analytic.jl")
include("tables.jl")
include("composite.jl")
include("entry.jl")

export BinnedSpectrum, nbins, bin_widths, bin_centres, density, from_density, rebin
export Parameter, model_dat_line
export AbstractXSPECModel, AdditiveModel, MultiplicativeModel, ConvolutionModel
export xspec_type, parameter_names, parameter_defaults, named_parameters, evaluate
export register!, registered_models, lookup_model, c_function_name
export model_dat_text, c_wrappers_text, write_xspec_package
export GridAxis, GridInterpolator, grid_size, evaluate_exact, cache_stats, empty_cache!
export disk_loads, cache_memory_used_bytes, cache_limit_bytes, cache_directory
export RedshiftKernel, GaussianKernel, cumulative
export convolution_matrix, convolve, convolve_fft, convolution_method
export gaussian_per_bin, powerlaw_per_bin
export OGIPTable, resolve_table_path, write_ogip_table
export Blurred

end
