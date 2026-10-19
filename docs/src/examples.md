# Examples

Three scripts under `examples/` exercise the package from Julia, with no
HEASOFT installation. Run them from the repository root:

```sh
julia --project=. examples/additive.jl
julia --project=. examples/convolution.jl
julia --project=. examples/grid.jl
```

Each one prints a short summary and stops with an error if the result is wrong.
The XSPEC session that uses the compiled reference models is on the
[Quick start](quickstart.md) page, and the comparisons against XSPEC's own
models are in [Verification](verification.md).

## An additive model

`examples/additive.jl` builds a Gaussian line. The function returns photons
per bin, and `evaluate` is what the XSPEC entry point calls.

```julia
line = AdditiveModel("example_gauss",
    [Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6, delta = 0.05),
     Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 10.0, delta = 0.05)]) do energy_edges, p
    gaussian_per_bin(energy_edges, p.LineE, p.Sigma)
end

per_bin = evaluate(line, edges, [6.4, 0.1])   # one photon in total
```

## A convolution

`examples/convolution.jl` blurs that kind of line with a Gaussian in
``g = E_{\rm obs}/E_{\rm em}``. The direct matrix and the FFT are asked for
the same spectrum.

```julia
line = BinnedSpectrum(edges, gaussian_per_bin(edges, 6.4, 0.05))
kernel = GaussianKernel(0.05)
direct = convolve(line, kernel; method = :direct)
fast = convolve(line, kernel; method = :fft)
```

## A cached model

`examples/grid.jl` pretends the Gaussian is expensive. It is computed only at
the corners of a `(LineE, Sigma)` grid, on a fixed internal energy grid, and
[`rebin`](@ref) moves the result onto the bins a caller asked for. `LineE` is
linear and `Sigma` is logarithmic.

```julia
cached = GridInterpolator((
    LineE = range(5.0, 8.0, length = 7),
    Sigma = logrange(0.05, 1.0, 5),
)) do LineE, Sigma
    gaussian_per_bin(INTERNAL, LineE, Sigma)
end

on_internal_grid = cached(6.4, 0.3)
delivered = rebin(BinnedSpectrum(INTERNAL, on_internal_grid), xspec_edges)
```
