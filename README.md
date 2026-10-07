# JuliaXSPEC.jl

Write [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) spectral models in Julia.

JuliaXSPEC lets XSPEC call Julia. You write a model as an ordinary Julia
function, JuliaXSPEC generates the `model.dat` and C glue that XSPEC needs,
and a compiled library carries the Julia runtime into your XSPEC session.
On top of that it provides two building blocks for expensive models:

- **Caching on a grid** — `GridInterpolator` evaluates a slow function only
  at the corners of a parameter grid and interpolates between them, like an
  XSPEC table model that is filled in on demand.
- **Convolution** — `convolve` blurs a spectrum with a kernel in the energy
  ratio `g = E_obs / E_em` (relativistic smearing, for example), with the
  bin-integrated bookkeeping that XSPEC uses made explicit.

The package is deliberately readable: each idea has the simplest correct
implementation first, and faster versions are added alongside, never instead.

**Documentation:** [dev](https://phajy.github.io/JuliaXSPEC.jl/dev/)

## Status

Phase 2: the XSPEC bridge, grid caching (in memory and on disk), direct and
FFT convolution, and OGIP table models. The reference models `jlgauss`,
`jlgausscached`, `jlgconv`, `jlgconvfft`, `jltable` and `jltableblur` run
inside XSPEC. `jltable` is checked against `atable{xillverD-5.fits}`, and the
FFT blur against the direct matrix (see the Verification page).

Planned: relativistic blurring via an optional
[Gradus.jl](https://codeberg.org/astro-group/Gradus.jl) extension, verification
against relxill, and examples calling external tools such as
[kerrz](https://git.sr.ht/~fjebaker/kerrz).

## A model in three lines

```julia
using JuliaXSPEC

jlgauss = AdditiveModel("jlgauss", [Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6),
                                    Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 10.0)]) do edges, p
    gaussian_per_bin(edges, p.LineE, p.Sigma)   # photons / cm² / s in each bin
end
register!(jlgauss)
```

## Building and using in XSPEC

Requirements: Julia 1.13 and HEASOFT with `HEADAS` initialised.

```sh
./build-julia.sh     # compile the Julia models into build/lib (several minutes)
./build-xspec.sh     # initpackage + hmake in xspec/
xspec
XSPEC12> lmod juliaxspec xspec
XSPEC12> model jlgconv * gaussian
```

See the documentation for the units conventions, a step-by-step account of
the convolution, and how to write and verify your own models.

## License

MIT. Gradus.jl and kerrz, which optional examples call, are GPL-3.0.
