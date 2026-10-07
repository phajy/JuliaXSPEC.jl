# Simple spectra with closed-form bin integrals. They make good first models
# and, because the answers are known exactly, good tests of everything else.

using SpecialFunctions: erf

"""
    gaussian_per_bin(edges, centre, sigma) -> Vector{Float64}

Photons per bin of a Gaussian emission line with unit total (one photon
cm⁻² s⁻¹ in all), integrated exactly over each bin:

    N_i = ½ [erf((E_{i+1} - centre) / (√2 sigma)) - erf((E_i - centre) / (√2 sigma))]

This is XSPEC's `gaussian` model with `norm = 1`. A `sigma` of zero puts the
whole photon into the bin containing `centre`, as XSPEC does.
"""
function gaussian_per_bin(edges::AbstractVector{<:Real}, centre::Real, sigma::Real)
    sigma >= 0 || throw(ArgumentError("sigma must be non-negative"))
    if sigma == 0
        per_bin = zeros(Float64, length(edges) - 1)
        i = searchsortedlast(edges, centre)
        1 <= i <= length(per_bin) && (per_bin[i] = 1.0)
        return per_bin
    end
    scaled = @. erf((edges - centre) / (sqrt(2) * sigma))
    return 0.5 .* diff(scaled)
end

"""
    powerlaw_per_bin(edges, photon_index) -> Vector{Float64}

Photons per bin of a power law `n(E) = E^(-photon_index)` photons cm⁻² s⁻¹
keV⁻¹ (unit normalisation at 1 keV), integrated exactly over each bin. This is
XSPEC's `powerlaw` model with `norm = 1`.
"""
function powerlaw_per_bin(edges::AbstractVector{<:Real}, photon_index::Real)
    if photon_index == 1
        return diff(log.(edges))
    end
    p = 1 - photon_index
    return diff(edges .^ p) ./ p
end
