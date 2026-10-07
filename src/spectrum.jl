# The one type that every JuliaXSPEC model speaks: a spectrum on a grid of
# energy bins, holding the number of photons in each bin. This is exactly what
# XSPEC passes to, and expects back from, a model function.

"""
    BinnedSpectrum(edges, per_bin)

A spectrum tabulated on energy bins, in XSPEC's convention.

- `edges`: the `n + 1` bin edges in keV, strictly increasing.
- `per_bin`: `n` values, each the number of photons cm⁻² s⁻¹ **in that bin**.

`per_bin` is a *bin-integrated* quantity, not a density. The photon flux
density (photons cm⁻² s⁻¹ keV⁻¹) is obtained with [`density`](@ref), which
divides by the bin widths. XSPEC's `plot model` shows the density; XSPEC's
model functions, table models, and convolution models all use `per_bin`.

See the "Units, bins and normalisation" page of the documentation.
"""
struct BinnedSpectrum{T<:Real}
    edges::Vector{T}
    per_bin::Vector{T}

    function BinnedSpectrum(edges::Vector{T}, per_bin::Vector{T}) where {T<:Real}
        length(edges) == length(per_bin) + 1 ||
            throw(ArgumentError("need length(edges) == length(per_bin) + 1, got $(length(edges)) and $(length(per_bin))"))
        all(>(0), diff(edges)) ||
            throw(ArgumentError("energy edges must be strictly increasing"))
        return new{T}(edges, per_bin)
    end
end

BinnedSpectrum(edges::AbstractVector, per_bin::AbstractVector) =
    BinnedSpectrum(collect(float.(edges)), collect(float.(per_bin)))

"""
    nbins(spectrum)

Number of energy bins.
"""
nbins(s::BinnedSpectrum) = length(s.per_bin)

"""
    bin_widths(edges)
    bin_widths(spectrum)

Width of each bin in keV.
"""
bin_widths(edges::AbstractVector) = diff(edges)
bin_widths(s::BinnedSpectrum) = bin_widths(s.edges)

"""
    bin_centres(edges)
    bin_centres(spectrum)

Arithmetic midpoint of each bin in keV.
"""
bin_centres(edges::AbstractVector) = (edges[1:end-1] .+ edges[2:end]) ./ 2
bin_centres(s::BinnedSpectrum) = bin_centres(s.edges)

"""
    density(spectrum) -> Vector

Photon flux density in each bin, photons cm⁻² s⁻¹ keV⁻¹: `per_bin ./ bin_widths`.
This is what XSPEC plots with `plot model`.
"""
density(s::BinnedSpectrum) = s.per_bin ./ bin_widths(s)

"""
    from_density(edges, density) -> BinnedSpectrum

Build a [`BinnedSpectrum`](@ref) from a photon flux density (photons cm⁻² s⁻¹
keV⁻¹) tabulated at each bin, multiplying by the bin widths. This assumes the
density is constant across each bin; use finer bins if that matters.
"""
from_density(edges::AbstractVector, density::AbstractVector) =
    BinnedSpectrum(edges, density .* bin_widths(edges))

"""
    rebin(spectrum, new_edges) -> BinnedSpectrum

Move a spectrum onto a different set of energy bins, conserving photons.

Each source bin's photons are assumed to be spread uniformly across the bin,
and the fraction overlapping each destination bin is handed over. Photons
falling outside the destination grid are dropped; the total is conserved
whenever the destination grid covers the source grid.
"""
function rebin(s::BinnedSpectrum, new_edges::AbstractVector)
    out = zeros(Float64, length(new_edges) - 1)
    src = s.edges
    for j in 1:nbins(s)
        lo, hi = src[j], src[j+1]
        width = hi - lo
        width > 0 || continue
        # Destination bins that overlap the source bin [lo, hi).
        first_bin = max(searchsortedlast(new_edges, lo), 1)
        last_bin = min(searchsortedfirst(new_edges, hi) - 1, length(out))
        for i in first_bin:last_bin
            overlap = min(hi, new_edges[i+1]) - max(lo, new_edges[i])
            overlap > 0 || continue
            out[i] += s.per_bin[j] * overlap / width
        end
    end
    return BinnedSpectrum(collect(float.(new_edges)), out)
end
