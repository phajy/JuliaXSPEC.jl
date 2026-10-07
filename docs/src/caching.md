# Caching on a grid

Many physical models are far too slow to evaluate at every step of a fit. The
traditional answer is an XSPEC table model: evaluate the model once on a grid
of parameter values, store the spectra in a FITS file, and let XSPEC
interpolate. [`GridInterpolator`](@ref) does the same thing without the
file: it evaluates the function only at the grid corners a fit actually
visits, keeps those results in memory, and interpolates between them.

## Using it

```julia
const INTERNAL_EDGES = collect(logrange(0.1, 100.0, 4001))

const CACHED_GAUSSIAN = GridInterpolator(
    (LineE = range(5.0, 8.0, length = 61), Sigma = logrange(0.05, 1.0, 14)),
) do LineE, Sigma
    gaussian_per_bin(INTERNAL_EDGES, LineE, Sigma)      # the expensive function
end

spectrum = CACHED_GAUSSIAN(6.43, 0.3)                     # a Vector, interpolated
```

- The `NamedTuple` names the parameters and gives each one's grid points.
- The function takes the parameters as positional arguments, in the same
  order, and returns a vector. Typically that vector is a spectrum on a fixed
  internal energy grid (as above) or a kernel on a fixed ``g`` grid; it must
  not depend on anything that changes between calls, such as XSPEC's current
  bins.
- Calling the interpolator with parameter values returns the multilinear
  interpolation between the ``2^N`` surrounding corners, computing any corner
  that has not been seen before.

In a model, rebin the interpolated internal spectrum onto XSPEC's bins:

```julia
AdditiveModel("jlgausscached", parameters) do edges, p
    internal = BinnedSpectrum(INTERNAL_EDGES, CACHED_GAUSSIAN(p.LineE, p.Sigma))
    rebin(internal, edges).per_bin
end
```

## Linear and logarithmic axes

Each axis is interpolated in its own coordinate. A `range` is interpolated in
the parameter value; a `logrange` in its logarithm; a plain vector is examined
by [`detect_scale`](@ref JuliaXSPEC.detect_scale) — uniform spacing means linear, uniform spacing in
the logarithm means log, and anything else is interpolated piecewise-linearly
between the given values. You can also say so explicitly:

```julia
GridAxis([2.0, 3.0, 5.0, 10.0]; scale = :log)
```

This matters for parameters that span decades (a corona height from 2 to 100
gravitational radii, an ionisation parameter from 1 to 10⁴): with a
logarithmic axis, a value halfway between grid points in ``\log x`` gets equal
weights, which is almost always what the spectra themselves do. XSPEC table
models record the same choice per parameter in their `METHOD` column.

## What interpolation can and cannot do

Multilinear interpolation **blends** the corner vectors with weights that sum
to one. It does not shift features. Halfway between grid points at
`LineE = 6.40` and `6.45` keV, the interpolated spectrum is not a line at
6.425 keV but the average of two lines — a slightly broader, slightly
flat-topped bump. The error shrinks as the grid spacing becomes small compared
to the width of the features, which is the criterion for choosing a grid:
*neighbouring corners must look alike*. For a reflection spectrum blurred by a
relativistic kernel this is easily satisfied; for a narrow Gaussian it is not,
which is exactly what `jlgausscached` versus `jlgauss` demonstrates
([Verification](verification.md)).

Parameter values outside the grid are clamped to its edge. Give the grid the
same limits as the XSPEC parameters (as `jlgausscached` does) and this never
happens during a fit.

## Measuring the error and watching the cache

```julia
exact  = evaluate_exact(CACHED_GAUSSIAN, 6.43, 0.3)      # no grid, no cache
approx = CACHED_GAUSSIAN(6.43, 0.3)
maximum(abs, approx - exact) / maximum(exact)

cache_stats(CACHED_GAUSSIAN)     # (hits = …, misses = …, stored = …, total = 854)
empty_cache!(CACHED_GAUSSIAN)
```

`misses` counts corners computed, `hits` corners reused. In a fit, the first
few evaluations compute the corners around the starting point and the rest
of the fit is almost free, as long as the parameters stay in the same region.

## Keeping corners on disk

Pass a cache name and each newly computed corner is written under
[`cache_directory`](@ref) (`~/.julia/juliaxspec` by default, or
`JULIAXSPEC_CACHE_DIR`):

```julia
GridInterpolator((spin = range(0, 0.998, 11),); cache = "lamppost-v1") do spin
    expensive_line_profile(spin)
end
```

The next session, or a new interpolator with the same name and the same grid,
reads the file instead of calling the function. [`disk_loads`](@ref) counts
those reads. [`empty_cache!`](@ref) forgets the RAM copy only. Changing the
grid under the same name deletes the old files, because they would no longer
be the right corners; changing the name (the `v1`) does the same when you have
changed the function rather than the grid.

All grid caches in the process share a RAM budget of
`JULIAXSPEC_CACHE_LIMIT_GB` gigabytes (default 16; `0` means no limit). Past
the budget, the least recently used corner of a cache is dropped from RAM.
The disk copy, if there is one, stays. [`cache_memory_used_bytes`](@ref)
reports the total.

## Table models use the same interpolator

An [OGIP table](https://heasarc.gsfc.nasa.gov/FTP/caldb/docs/memos/ogip_92_009/ogip_92_009.pdf)
is a grid whose corners were computed in advance and stored as photons per
bin. [`OGIPTable`](@ref) reads one and interpolates it with a
`GridInterpolator`, taking the linear-or-log choice from each parameter's
METHOD flag rather than from the spacing. See [Writing a model](models.md).

A machine-learning emulator trained on the same corners could later replace
the multilinear blend without changing the models that call the interpolator.
