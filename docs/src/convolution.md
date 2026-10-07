# Convolution

This page builds up, step by step, exactly what JuliaXSPEC's [`convolve`](@ref)
does: from the contract XSPEC imposes on a convolution model, through the
continuous physics, to the discrete operation on bins and the code that
implements it. The aim is that nothing about units, normalisation or edge
effects is left to guesswork.

The running example is relativistic blurring of an emission line or a
reflection spectrum, but the same machinery applies to any kernel that acts
on the energy *ratio*.

## 1. What XSPEC hands a convolution model

An XSPEC convolution component (`con`), used as `jlgconv * gaussian`, is
evaluated after the components it wraps. XSPEC passes it the current energy
bin edges ``E_1 < \dots < E_{n+1}`` and an array of ``n`` numbers: the photons
cm⁻² s⁻¹ **in each bin** produced by the wrapped components, with their
normalisations already applied. The model must overwrite that array with the
blurred photons per bin **on the same bins**. It is not told anything about
photons outside the current energy range, and it cannot return anything
outside it either.

In JuliaXSPEC a [`ConvolutionModel`](@ref)'s function receives this as a
[`BinnedSpectrum`](@ref) and returns one with identical `edges`
(see [Units, bins and normalisation](units.md)).

## 2. Two kinds of smearing

A smearing kernel can act on an energy *difference* or on an energy *ratio*.

- **Additive in energy.** Detector resolution, or XSPEC's `gsmooth`, spreads a
  photon of energy ``E_e`` to ``E_o = E_e + \delta`` with some probability
  ``K(\delta)``. The kernel has units of 1/keV and the operation is an
  ordinary convolution in ``E``.
- **Multiplicative in energy.** Doppler shifts, gravitational redshift and
  cosmological redshift all *multiply* a photon's energy:
  ``E_o = g\,E_e`` with ``g = E_o / E_e``. A photon emitted at 6.4 keV from
  a disc moving towards us may arrive at 7 keV; the same motion takes a 1 keV
  photon to 1.09 keV, not to 1.6 keV. The kernel ``L(g)`` is dimensionless in
  ``g``, and the natural variable is ``\ln E``, where a fixed ``g`` is a fixed
  shift ``\ln g``.

Relativistic blurring is the second kind. Everything below is written for a
kernel ``L(g)`` in the ratio ``g``; if you need the first kind, the
derivation is the same with ``E_o - E_e`` in place of ``E_o / E_e`` and no
Jacobian factor.

## 3. The continuous statement, and what is conserved

Write ``n_{\rm em}(E)`` for the emitted photon flux density (photons cm⁻² s⁻¹
keV⁻¹) and ``L(g)`` for the kernel, normalised so that

```math
\int_0^\infty L(g)\, dg = 1 .
```

``L(g)\,dg`` is the probability that a photon arrives with its energy
multiplied by a factor between ``g`` and ``g + dg``. The observed density is
then

```math
n_{\rm obs}(E_o) = \int_0^\infty L\!\left(\frac{E_o}{E_e}\right) n_{\rm em}(E_e)\,\frac{dE_e}{E_e} .
```

The factor ``1/E_e`` is the Jacobian from ``g`` to ``E_o``: for fixed
``E_e``, ``E_o = g E_e`` so ``dE_o = E_e\,dg`` — the same spread in ``g``
covers a wider range of observed energies when ``E_e`` is larger, and the
density must be correspondingly lower. Substituting ``\tau = \ln E`` turns the
integral into an ordinary convolution in ``\tau`` with the kernel
``L(e^{\tau_o - \tau_e})``, which is the basis of the fast implementation in
step 6.

**Photon number is conserved.** Integrate the observed density over all
energies and swap the order of integration:

```math
\int n_{\rm obs}(E_o)\,dE_o
= \int n_{\rm em}(E_e)\left[\int L\!\left(\frac{E_o}{E_e}\right)\frac{dE_o}{E_e}\right] dE_e
= \int n_{\rm em}(E_e)\left[\int L(g)\,dg\right] dE_e
= \int n_{\rm em}(E_e)\,dE_e .
```

Every emitted photon is observed somewhere; unit area of the kernel is
exactly the statement that no photons are created or lost. This is the
normalisation convention JuliaXSPEC enforces when a
[`RedshiftKernel`](@ref) is constructed.

**Energy flux is not conserved.** Repeating the calculation with an extra
factor ``E_o = g E_e`` gives
``\int E_o\, n_{\rm obs}\,dE_o = \langle g\rangle \int E_e\, n_{\rm em}\,dE_e``,
where ``\langle g \rangle = \int g\,L(g)\,dg``. A kernel whose mean lies below
one (as for gravitational redshift) dims the energy flux even though it
conserves photons. A model that normalises its kernel to unit *energy* flux is
therefore making a different, and for XSPEC's purposes wrong, choice.

## 4. From densities to bins

XSPEC gives us photons per bin, not densities. Define the input and output
bin contents

```math
N^{\rm in}_j = \int_{a_j}^{b_j} n_{\rm em}(E)\,dE, \qquad
N^{\rm out}_i = \int_{c_i}^{d_i} n_{\rm obs}(E)\,dE ,
```

where input bin ``j`` is ``[a_j, b_j]`` and output bin ``i`` is ``[c_i, d_i]``
(the two grids are the same in an XSPEC convolution model, but need not be).
Because the operation is linear, the output is a matrix applied to the input:

```math
N^{\rm out}_i = \sum_j M_{ij}\, N^{\rm in}_j ,
```

and ``M_{ij}`` has a direct meaning: **the fraction of the photons in input
bin ``j`` that are observed in output bin ``i``.**

To write it down, note that a single photon emitted at ``E_e`` lands in
``[c_i, d_i]`` with probability

```math
\int_{c_i}^{d_i} L\!\left(\frac{E_o}{E_e}\right)\frac{dE_o}{E_e}
= \Lambda\!\left(\frac{d_i}{E_e}\right) - \Lambda\!\left(\frac{c_i}{E_e}\right),
\qquad \Lambda(g) = \int_0^g L(g')\,dg' ,
```

where ``\Lambda`` is the kernel's cumulative distribution
([`cumulative`](@ref)), rising from 0 to 1. The photons in input bin ``j``
are spread over emitted energies according to ``n_{\rm em}`` within the bin,
which we do not know; the only information XSPEC gives us is the bin total.
We therefore make the same assumption as `rebin`: that they are spread
**uniformly** across the bin. Averaging over the bin gives

```math
M_{ij} = \frac{1}{b_j - a_j}\int_{a_j}^{b_j}
\left[\Lambda\!\left(\frac{d_i}{E_e}\right) - \Lambda\!\left(\frac{c_i}{E_e}\right)\right] dE_e .
```

Two properties follow immediately:

- ``0 \le M_{ij} \le 1``, since each entry is a probability.
- ``\sum_i M_{ij} = 1`` whenever the output grid covers every energy to which
  the kernel can carry input bin ``j``, because the ``\Lambda`` differences
  then telescope to ``\Lambda(\infty) - \Lambda(0) = 1``. Columns that sum to
  less than one correspond to photons shifted off the grid (step 7).

The uniform-within-bin assumption is the one approximation in this picture.
Its error is second order in the bin width relative to the kernel's width, so
it is negligible on the fine grids typical of X-ray responses, and it can be
made as small as you like by evaluating on a finer internal grid and
rebinning (step 9).

## 5. The direct implementation

[`convolution_matrix`](@ref) is a transcription of the formula above. The
average over the input bin is done with the midpoint rule on `n_sub` equal
sub-bins (4 by default), and for each sample energy only the output bins the
kernel can reach are visited:

```julia
for j in 1:n_in
    a, b = in_edges[j], in_edges[j+1]
    for s in 1:n_sub
        E_em = a + (s - 0.5) * (b - a) / n_sub
        # Only output bins overlapping [g_min E_em, g_max E_em] can receive these photons.
        first_bin = max(searchsortedlast(out_edges, g_min * E_em), 1)
        last_bin = min(searchsortedfirst(out_edges, g_max * E_em) - 1, n_out)
        for i in first_bin:last_bin
            c, d = out_edges[i], out_edges[i+1]
            M[i, j] += (cumulative(kernel, d / E_em) - cumulative(kernel, c / E_em)) / n_sub
        end
    end
end
```

`cumulative` is exact for the tabulated kernel: a [`RedshiftKernel`](@ref) is
linear between its grid points, so ``\Lambda`` is quadratic between them, and
the running sums at the grid points are stored when the kernel is built.
[`convolve`](@ref) then simply forms `M * per_bin`.

The cost is proportional to (input bins) × (output bins within the kernel's
reach) × `n_sub`. For a 2000-bin response and a kernel spanning
``g \in [0.7, 1.3]`` that is a few million cheap operations, well under a
second. For very broad kernels (a maximally spinning black hole spreads
``g`` over roughly 0.1 to 1.5) and fine grids the matrix becomes large, which
motivates the next step. The direct method remains the reference against
which anything faster is checked.

## 6. A faster implementation (planned)

Because the operation is a translation convolution in ``\tau = \ln E``
(step 3), it can be done with a fast Fourier transform: rebin the input onto a
uniform grid in ``\ln E``, convert to density per unit ``\tau``, multiply the
transforms of the spectrum and of ``L(e^{\tau})``, transform back, and rebin
onto the output grid. The cost then grows only as ``m \log m`` with the size
``m`` of the log grid, independent of the kernel width. This will be added in
a later phase alongside the direct method, with tests that the two agree; the
direct method stays as the readable definition.

## 7. Band edges: photons that leave or arrive

A kernel with ``g < 1`` support moves photons *down* in energy; the photons
just above the top of the energy grid would land inside it, and those just
above the bottom of the grid are carried below it. Within an XSPEC convolution
model the first group is simply unavailable (XSPEC never told us about them)
and the second is lost (we cannot return them). The effect is confined to a
band of relative width ``\max(1 - g_{\min},\ g_{\max} - 1)`` at each end of
the grid.

The remedy, as for XSPEC's own `kdblur`, `rdblur` and `relconv`, is XSPEC's
`energies` command, which evaluates the model on a grid extended beyond the
response, for example

```text
XSPEC12> energies extend low 0.01 100 log
XSPEC12> energies extend high 500 100 log
```

Within JuliaXSPEC, composite models that compute the source themselves
(step 9) avoid the problem by evaluating on an internal grid that already
covers the kernel's reach.

The tests check column sums of ``M`` in the interior of the grid: they are 1
to better than ``10^{-6}``.

## 8. Where kernels come from, and their normalisation

A [`RedshiftKernel`](@ref) is just a table of ``L`` on a grid of ``g``; its
constructor rescales the table to unit area, so any non-negative shape can be
supplied.

- **Gaussian in ``g``** ([`GaussianKernel`](@ref)): ``L(g) \propto
  \exp[-(g-1)^2 / 2\sigma_g^2]``, truncated at ``g > 0``. Useful because its
  effect is known in closed form (step 10).
- **A relativistic line profile.** Ray-tracing codes such as Gradus or kerrz
  compute the observed spectrum of a disc emitting a *single* line at rest
  energy ``E_{\rm rest}``. Expressed against ``g = E_o / E_{\rm rest}`` and
  normalised to unit area, that profile *is* ``L(g)``: the probability
  distribution of the energy shift a photon receives from the disc, for a
  given spin, inclination and emissivity. This is what the planned Gradus
  extension provides, cached on a parameter grid with
  [`GridInterpolator`](@ref).

XSPEC's `relconv`, `kdblur` and `rdblur` are convolution models of exactly
this kind. relxill's documentation states that `relconv` conserves photon
number, matching the convention here; when comparing with other codes it is
always worth checking whether photons or energy flux are conserved, because
the two choices differ by the factor ``\langle g \rangle`` of step 3.

## 9. Composite models: blurring the source yourself (planned)

Instead of exposing the kernel as a `con` component and letting XSPEC supply
the input, a Julia model can compute the source spectrum (for example by
interpolating an OGIP reflection table), blur it on the table's own fine
energy grid — which extends well beyond any response — and `rebin` the result
onto XSPEC's bins. This sidesteps the band-edge issue of step 7 and the
uniform-within-bin approximation of step 4 at the same time, at the price of
fixing the source model inside the Julia component. JuliaXSPEC will offer both
routes; the `con` route is what allows direct comparison with `relconv`.

## 10. Checks

Each statement above has a corresponding test in `test/runtests.jl`:

- **Conservation:** interior column sums of ``M`` equal 1.
- **Identity:** a kernel much narrower than a bin returns the input.
- **Pure shift:** a narrow kernel centred on ``g_0`` moves a line from
  ``E_0`` to ``g_0 E_0``, keeping its total.
- **Analytic Gaussian:** a Gaussian kernel of width ``\sigma_g`` applied to a
  narrow line at ``E_0`` gives a Gaussian line of width ``\sigma_g E_0``
  (exactly, since ``n_{\rm obs}(E_o) = L(E_o/E_0)/E_0`` for a delta-function
  input); the test compares with [`gaussian_per_bin`](@ref) and recovers the
  result to better than 0.2%.
- **Power law:** blurring ``E^{-\Gamma}`` returns ``E^{-\Gamma}`` multiplied
  by ``\int L(g)\,g^{\Gamma - 1}\,dg`` — a power law has no features to blur,
  but its normalisation changes, exactly as step 3 predicts.
- **In XSPEC:** `jlgconv * gaussian` with a narrow line is compared against
  XSPEC's own `gaussian` with the predicted width; see
  [Verification](verification.md).
