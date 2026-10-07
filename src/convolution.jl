# Convolution of a binned spectrum with a kernel in g = E_obs / E_em.
#
# The physics, for a photon flux density n(E) (photons cm⁻² s⁻¹ keV⁻¹):
#
#     n_obs(E_o) = ∫ L(E_o / E_e) n_em(E_e) dE_e / E_e,        ∫ L(g) dg = 1.
#
# Every photon emitted at E_e is observed at g E_e with probability L(g) dg;
# the 1/E_e converts dg into dE_o (since E_o = g E_e, dE_o = E_e dg) and makes
# photon number come out conserved. In bins, the operation is a matrix:
# M[i, j] is the fraction of the photons in input bin j that land in output
# bin i, built from the cumulative kernel Λ. This file is the direct, readable
# implementation; see the "Convolution" documentation page for the derivation.

"""
    convolution_matrix(in_edges, out_edges, kernel; n_sub = 4) -> Matrix{Float64}

The matrix `M` with `M[i, j]` = fraction of the photons in input bin `j`
(edges `in_edges[j]`, `in_edges[j+1]`) that are observed in output bin `i`.
Then `per_bin_out = M * per_bin_in`.

For a photon emitted at `E_e`, the fraction observed in `[c, d]` is
`Λ(d / E_e) - Λ(c / E_e)` where `Λ` is the kernel's cumulative distribution
([`cumulative`](@ref)). The photons of an input bin are assumed to be spread
uniformly over the bin, which is approximated by averaging that fraction over
`n_sub` emitted energies at the midpoints of equal sub-bins.

Each column sums to 1 (photons are conserved) whenever the output grid covers
every energy the kernel can shift the input bin to; photons shifted beyond
the output grid are lost.
"""
function convolution_matrix(
    in_edges::AbstractVector{<:Real},
    out_edges::AbstractVector{<:Real},
    kernel::RedshiftKernel;
    n_sub::Integer = 4,
)
    n_sub >= 1 || throw(ArgumentError("n_sub must be at least 1"))
    n_in = length(in_edges) - 1
    n_out = length(out_edges) - 1
    M = zeros(Float64, n_out, n_in)
    g_min, g_max = kernel.g[1], kernel.g[end]

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
    return M
end

"""
    convolve(spectrum, kernel; out_edges = spectrum.edges, n_sub = 4) -> BinnedSpectrum

Blur a [`BinnedSpectrum`](@ref) with a [`RedshiftKernel`](@ref): every photon
at energy `E` is redistributed to `g E` with probability `L(g) dg`. The result
is on `out_edges` (by default the input bins, which is what an XSPEC
convolution model must return).

This is the direct implementation, `convolution_matrix(...) * per_bin`: simple
and exact up to the sub-bin quadrature, with cost proportional to the number
of input bins times the number of output bins within the kernel's reach.
"""
function convolve(
    spectrum::BinnedSpectrum,
    kernel::RedshiftKernel;
    out_edges::AbstractVector{<:Real} = spectrum.edges,
    n_sub::Integer = 4,
)
    M = convolution_matrix(spectrum.edges, out_edges, kernel; n_sub)
    return BinnedSpectrum(collect(float.(out_edges)), M * spectrum.per_bin)
end
