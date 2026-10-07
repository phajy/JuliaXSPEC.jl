# Quick start

This session assumes both build steps are done ([Installation and build](build.md))
and that you start `xspec` in the repository root.

## The reference models

| Model | Type | Parameters | What it demonstrates |
|-------|------|------------|----------------------|
| `jlgauss` | add | `LineE`, `Sigma` | A model evaluated directly: a Gaussian line, identical to XSPEC's `gaussian` |
| `jlgausscached` | add | `LineE`, `Sigma` | The same line interpolated from a grid of cached evaluations |
| `jlgconv` | con | `SigmaG` | Gaussian blur in ``g``, by the direct matrix |
| `jlgconvfft` | con | `SigmaG` | The same blur, by FFT |
| `jltable` | add | `Gamma`, `A_Fe`, `logXi`, `Dens`, `Incl` | The `xillverD-5` reflection table |
| `jltableblur` | add | `SigmaG` plus the table parameters | That table, blurred on its own energy grid |

## A first look

```text
XSPEC12> lmod juliaxspec xspec
XSPEC12> dummyrsp 0.1 50 2000 log
XSPEC12> model jlgauss
   (accept the defaults with /*)
XSPEC12> setplot energy
XSPEC12> cpd /xs
XSPEC12> plot model
```

The first evaluation starts the Julia runtime. Everything you would do with
`gaussian` works the same way — `newpar`, `freeze`, `fit`, `fakeit` — because
to XSPEC a JuliaXSPEC model is just another local model.

## Direct versus cached evaluation

`jlgausscached` computes the Gaussian only at the corners of a 61 × 14 grid in
(`LineE`, `Sigma`) — linearly spaced in `LineE`, logarithmically in `Sigma` —
and interpolates between them. Compare it with the direct model:

```text
XSPEC12> model jlgauss + jlgausscached
   (set LineE = 6.43 and Sigma = 0.3 for both; norm = 1 and -1 to see the difference)
XSPEC12> plot model
```

Interpolation blends neighbouring spectra rather than shifting them, so the
error is largest for narrow lines and parameter values far from the grid
points. Set `JULIAXSPEC_VERBOSE=1` before starting XSPEC to see when new
corners are computed and when cached ones are reused.

## Convolution

A Gaussian in ``g`` applied to a narrow line at ``E_0`` gives a Gaussian in
energy of width ``\sigma_g E_0``:

```text
XSPEC12> model jlgconv * gaussian
   SigmaG = 0.05, LineE = 6.4, Sigma = 0.01
XSPEC12> plot model
```

should be indistinguishable from `gaussian` with `Sigma = 0.32` (precisely
``\sqrt{(0.05 \times 6.4)^2 + 0.01^2}``). Try `jlgconv * powerlaw` too: the
power law keeps its slope and only changes normalisation, as the
[Convolution](convolution.md) page explains.

Because `jlgconv` is an XSPEC `con` component it can wrap anything, including
table models: `jlgconv * atable{xillverD-5.fits}` blurs a reflection
spectrum.

`jltable` needs `xillverD-5.fits` (see [Installation and build](build.md)).
Its first evaluation reads the file; after that it is ordinary interpolation.
`jlgconvfft * gaussian` should overlie `jlgconv * gaussian`.
`jltableblur` and `jlgconv * atable{xillverD-5.fits}` are the same blur
arranged two ways, and they agree through the middle of the band.

## Scripts

Ready-made versions of these sessions, which also write the model values to
files for comparison, are in `verification/`; see [Verification](verification.md).
