"""
    GradusXSPEC

Relativistic kernels from [Gradus.jl](https://codeberg.org/astro-group/Gradus.jl),
registered as JuliaXSPEC models.

The first model is `jllamp`: a lamppost on the spin axis above a thin disc.
It is a convolution, so in XSPEC it smears whatever component it multiplies
(`jllamp * gaussian`, `jllamp * powerlaw`). The line profile is calculated on
a small grid of spin, inclination and height and remembered, on disk, under
the cache name `jllamp-thin-v2`.
"""
module GradusXSPEC

using Gradus
using JuliaXSPEC
using JuliaXSPECModels

# g = E_obs / E_em. Wide enough for a thin-disc lamppost line, and the same
# grid at every grid corner so the cached profiles can be blended.
const G_GRID = collect(range(0.1, 1.5; length = 256))

# Coarse on purpose. The default parameters sit on grid points, so the first
# XSPEC plot calculates one profile. A finer grid can replace this later.
const LAMP_SPIN = [0.2, 0.5, 0.8, 0.998]
const LAMP_INC = [10.0, 30.0, 45.0, 60.0]
const LAMP_HEIGHT = [3.0, 5.0, 10.0, 20.0]

"""
    lamppost_profile(spin, inclination, height) -> Vector{Float64}

Unit-area line profile `L(g)` on `G_GRID` for a lamppost of the given
height, in gravitational radii, above a thin disc. `inclination` is in degrees
from face-on.
"""
function lamppost_profile(spin::Real, inclination::Real, height::Real)
    metric = KerrMetric(M = 1.0, a = Float64(spin))
    disc = ThinDisc(0.0, Inf)
    corona = LampPostModel(h = Float64(height))
    # Gradus's default. For a lamppost this is the number of polar-angle rays
    # from 0.01° to 179.99°, not a Monte Carlo count. 128 leaves the inner
    # disc too coarsely sampled.
    emissivity = emissivity_profile(metric, disc, corona; n_samples = 1000)
    observer = SVector(0.0, 1000.0, deg2rad(Float64(inclination)), 0.0)
    _, flux = lineprofile(
        G_GRID,
        scalar_emissivity(emissivity),
        metric,
        observer,
        disc;
        method = TransferFunctionMethod(),
        maxrₑ = 400.0,
    )
    L = max.(Float64.(flux), 0.0)
    area = 0.0
    for i in 1:(length(G_GRID) - 1)
        area += 0.5 * (L[i] + L[i + 1]) * (G_GRID[i + 1] - G_GRID[i])
    end
    area > 0 || error("lamppost line profile has zero area at spin=$spin, inc=$inclination, h=$height")
    return L ./ area
end

# Lamppost profiles are already scalars. The ring profiles Gradus returns are
# length-1 vectors, and this accepts either so a later ring model can reuse it.
scalar_emissivity(profile) = function (r)
    ε = emissivity_at(profile, r)
    return ε isa Number ? Float64(ε) : Float64(first(ε))
end

const LAMP_PROFILES = GridInterpolator(
    (spin = LAMP_SPIN, inc = LAMP_INC, h = LAMP_HEIGHT);
    cache = "jllamp-thin-v2",
) do spin, inc, h
    lamppost_profile(spin, inc, h)
end

const jllamp = ConvolutionModel(
    "jllamp",
    [
        Parameter("spin", 0.8; min = 0.2, max = 0.998, delta = 0.05),
        Parameter("inc", 30.0; unit = "degrees", min = 10.0, max = 60.0, delta = 5.0),
        Parameter("h", 10.0; unit = "GM/c^2", min = 3.0, max = 20.0, delta = 1.0),
    ];
    description = "Thin-disc lamppost blur. Profiles are cached on a coarse grid.",
) do spectrum, p
    L = max.(LAMP_PROFILES(p.spin, p.inc, p.h), 0.0)
    convolve(spectrum, RedshiftKernel(G_GRID, L); method = :fft)
end

function __init__()
    register!(jllamp)
end

end
