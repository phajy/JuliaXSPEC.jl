# Writing a model

A JuliaXSPEC model is a Julia function plus a description of its parameters.
This page covers the three model types, the parameter description, the
registry, and how a model gets into XSPEC.

## Additive models

An additive model computes photons cm⁻² s⁻¹ in each bin for unit
normalisation ([Units, bins and normalisation](units.md)). Its function
receives the bin edges and a `NamedTuple` of parameters:

```julia
using JuliaXSPEC

jlgauss = AdditiveModel(
    "jlgauss",
    [
        Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6, delta = 0.05),
        Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 20.0, soft_max = 10.0, delta = 0.05),
    ];
    description = "Gaussian line",
) do edges, p
    gaussian_per_bin(edges, p.LineE, p.Sigma)
end
```

The function must return `length(edges) - 1` values. `norm` is handled by
XSPEC and is not a parameter of your model.

## Multiplicative models

A multiplicative model returns a dimensionless factor per bin:

```julia
jlabs = MultiplicativeModel("jlabs", [Parameter("nH", 1.0; unit = "1e22", min = 0.0, max = 1e6)]) do edges, p
    E = bin_centres(edges)
    @. exp(-p.nH * 1e22 * cross_section(E))
end
```

## Convolution models

A convolution model receives the spectrum of the components it wraps as a
[`BinnedSpectrum`](@ref) and must return a `BinnedSpectrum` on the same bins:

```julia
jlgconv = ConvolutionModel("jlgconv", [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5, delta = 0.005)]) do spectrum, p
    convolve(spectrum, GaussianKernel(p.SigmaG))
end
```

See [Convolution](convolution.md) for what `convolve` does and for the
band-edge caveats of `con` components.

## Parameters

[`Parameter`](@ref) carries what XSPEC's `model.dat` needs:

```julia
Parameter(name, default; unit = "", min, max, soft_min = min, soft_max = max,
          delta = 1% of the range, frozen = false, kind = :fit)
```

- `min`/`max` are the hard limits: XSPEC never asks for values outside them,
  so your function can rely on that.
- `delta` is the initial fit step; `frozen = true` makes it negative, which
  tells XSPEC to start with the parameter frozen.
- `kind = :scale` or `:switch` produce XSPEC's `*` and `$` parameters.

Names must be identifiers (letters, digits, underscores) because they become
the fields of the `NamedTuple` your function receives.

## Model names

Model names may contain letters, digits and underscores. XSPEC's `initpackage`
wants the C wrapper functions named without underscores, so the C function
for `my_model` is `mymodel`; two names that differ only by underscores are
therefore rejected. The reference models use a `jl` prefix so they are easy
to tell apart from XSPEC's built-in models.

## Registering, and getting into XSPEC

`register!` makes a model known. Do it in a module's `__init__` so that it
also happens when the module is loaded from a precompiled image:

```julia
module MyModels
using JuliaXSPEC
const jlgauss = AdditiveModel(...) do edges, p ... end
function __init__()
    register!(jlgauss)
end
end
```

The reference set in `models/src/JuliaXSPECModels.jl` is the working example;
the simplest way to add models is to put them there and run
`./build-julia.sh` and `./build-xspec.sh` again. From the registry,
[`model_dat_text`](@ref) and [`c_wrappers_text`](@ref) generate the two files
XSPEC needs — you can call them at the REPL to see exactly what will be
written.

## Testing a model without XSPEC

[`evaluate`](@ref) runs a model exactly as the entry point does:

```julia
edges = collect(logrange(0.1, 50.0, 2001))
per_bin = evaluate(jlgauss, edges, [6.4, 0.3])
sum(per_bin)                      # 1.0: one photon per cm² per s

blurred = evaluate(jlgconv, BinnedSpectrum(edges, per_bin), [0.05])
```

The C entry point itself, [`juliaxspec_evaluate`](@ref JuliaXSPEC.juliaxspec_evaluate), can also be called
from Julia with pointers (see `test/runtests.jl`), which is how the test suite
checks the bridge without XSPEC.

## Expensive models

If your model is slow, compute it on a fixed internal energy grid at the
corners of a parameter grid with [`GridInterpolator`](@ref), and `rebin`
onto XSPEC's bins — see [Caching on a grid](caching.md). The `jlgausscached`
reference model shows the pattern in a dozen lines. Give the interpolator
`cache = "a-name"` to keep those corners on disk between XSPEC sessions.

## Table models

[`OGIPTable`](@ref) reads an additive XSPEC table file. The reference model
`jltable` is `xillverD-5.fits`, found via [`resolve_table_path`](@ref) (the
repository root, the current directory, or `JULIAXSPEC_TABLE_DIR`). The first
evaluation reads the whole file, about a gigabyte; later evaluations only
interpolate.

[`Blurred`](@ref) builds an additive model from a table and a kernel, blurring
on the table's energy grid so that photons outside the instrument band can
still redshift into it. `jltableblur` is that model with a Gaussian kernel.
See [Convolution](convolution.md), step 9.
