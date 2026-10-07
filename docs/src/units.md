# Units, bins and normalisation

Almost every confusing moment in writing spectral models comes from mixing up
*photons in a bin* with *photons per keV*. This page fixes the conventions
JuliaXSPEC uses, which are XSPEC's conventions, in one place. Everything else
in the package refers back here.

## What XSPEC gives a model, and what it wants back

When XSPEC evaluates a model component it passes

- the **energy bin edges** ``E_1 < E_2 < \dots < E_{n+1}`` in keV (so there
  are ``n`` bins, and bin ``i`` spans ``[E_i, E_{i+1}]``), and
- the **parameter values**, in the order they appear in `model.dat`.

It expects back ``n`` numbers, one per bin: the number of photons
cm⁻² s⁻¹ **in that bin**, for unit normalisation. XSPEC then multiplies by
its own `norm` parameter (additive models only) and folds the result through
the instrument response to predict counts.

In JuliaXSPEC these ``n`` numbers are the `per_bin` field of a
[`BinnedSpectrum`](@ref):

```julia
s = BinnedSpectrum(edges, per_bin)   # edges: n+1 keV values, per_bin: n photons/cm²/s
```

The bins can be anything: XSPEC uses the energy grid of the loaded response,
or whatever `dummyrsp` or `energies` set up, and the grid can change between
calls. A model must therefore never assume a particular binning.

## Per-bin counts versus density

The quantity most people have in mind when they picture a spectrum is the
photon flux **density** ``n(E)``, in photons cm⁻² s⁻¹ keV⁻¹. The two are
related by integrating over the bin:

```math
N_i = \int_{E_i}^{E_{i+1}} n(E)\, dE .
```

If the density is roughly constant across a bin, ``N_i \approx n(E_i)\,\Delta E_i``,
which is what [`density`](@ref) and [`from_density`](@ref) implement:

```julia
density(s)                 # per_bin ./ bin_widths(s)   photons/cm²/s/keV
from_density(edges, n)     # n .* bin_widths(edges)     back to photons per bin
```

Two consequences worth remembering:

- **Per-bin values depend on the binning; densities do not.** Halve the bin
  width and each `per_bin` value halves, while `density` is unchanged. This is
  why `plot model` in XSPEC (which shows density) looks smooth on any grid,
  and why you cannot compare `per_bin` arrays from two different grids.
- **Normalisation is a sum, not an integral.** The total photon flux of a
  component is simply `sum(s.per_bin)`. For XSPEC's `gaussian`, that sum is
  exactly `norm`; for `powerlaw`, `norm` is the density at 1 keV, and the
  per-bin values are the exact integrals ``\int E^{-\Gamma} dE`` over each bin
  (see [`powerlaw_per_bin`](@ref)).

When a model has an analytic density, integrate it over each bin exactly if
you can — [`gaussian_per_bin`](@ref) uses the error function — rather than
evaluating the density at the bin centre and multiplying by the width. The
latter is fine for smooth spectra on fine grids but breaks down for narrow
features.

## The three model types

- An **additive** model (`add`) returns photons per bin for unit `norm`.
- A **multiplicative** model (`mul`) returns a dimensionless factor per bin
  (for example ``e^{-\sigma(E) N_H}``), by which XSPEC multiplies the
  additive components it wraps. Being a ratio, it is the same whether you
  think in counts or density, as long as the factor varies little across a
  bin.
- A **convolution** model (`con`) receives the photons per bin of the
  components it wraps and returns photons per bin on the *same* bins. In
  JuliaXSPEC it receives and returns a `BinnedSpectrum`. What it may do with
  the energies outside the current grid is the subject of the
  [Convolution](convolution.md) page.

## Moving between grids

Because `per_bin` values are counts, moving a spectrum to a different set of
bins is a matter of handing each bin's photons to the destination bins it
overlaps, in proportion to the overlap. That is what [`rebin`](@ref) does.
It conserves the total exactly whenever the destination grid covers the
source grid, and assumes photons are spread uniformly within each source bin
— so rebin *fine to coarse* whenever you have the choice.

```julia
internal = BinnedSpectrum(fine_edges, expensive_model(fine_edges))
xspec_ready = rebin(internal, edges)       # edges: whatever XSPEC asked for
```

## Table models

XSPEC table models (OGIP 92-009 `atable`/`mtable` files) store spectra in
exactly this per-bin form: each row of the `SPECTRA` extension holds the
photons cm⁻² s⁻¹ in each of the table's energy bins, for one combination of
parameter values. Reading one into a `BinnedSpectrum` on the table's own
energy grid therefore involves no unit conversion at all; XSPEC's `atable`
evaluation is the interpolation between parameter rows followed by a
`rebin` onto the current energy grid.
