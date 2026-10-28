"""
    GradusXSPEC

Relativistic kernels from [Gradus.jl](https://codeberg.org/astro-group/Gradus.jl),
registered as JuliaXSPEC models.

`jllamp` and `jlring` are convolutions: in XSPEC each smears whatever
component it multiplies. `jllampxill` and `jlringxill` are those kernels
blurring `xillverD-5.fits` on the table's own energy grid. Each line profile
is calculated on a small parameter grid and remembered on disk (`jllamp-thin-v2`,
`jlring-thin-v1`).
"""
module GradusXSPEC

using Gradus
using JuliaXSPEC
using JuliaXSPECModels

# g = E_obs / E_em. The same grid at every cached corner, so neighbouring
# profiles can be blended.
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
    return unit_area(observed_flux(metric, disc, emissivity, inclination),
                     "lamppost at spin=$spin, inc=$inclination, h=$height")
end

# Lamppost profiles are already scalars. A ring profile returns a length-1
# vector even for one radius.
scalar_emissivity(profile) = function (r)
    ε = emissivity_at(profile, r)
    return ε isa Number ? Float64(ε) : Float64(first(ε))
end

function observed_flux(metric, disc, emissivity, inclination)
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
    return flux
end

function unit_area(flux, what)
    L = max.(Float64.(flux), 0.0)
    area = 0.0
    for i in 1:(length(G_GRID) - 1)
        area += 0.5 * (L[i] + L[i + 1]) * (G_GRID[i + 1] - G_GRID[i])
    end
    area > 0 || error("$what has zero area")
    return L ./ area
end

# Gradus takes log2 of a ring's branch emissivities. A slightly negative sample
# otherwise throws. Replace those with the smallest positive value of their type.
function clamp_ring_emissivity!(profile)
    for branch_group in profile.branches
        for branch in branch_group
            for i in eachindex(branch.ε)
                branch.ε[i] > 0 || (branch.ε[i] = eps(typeof(branch.ε[i])))
            end
        end
    end
    return profile
end

# Coarse, and the defaults sit on grid points. Height starts at 2.5: lower than
# that, Gradus can sample the disc inside the horizon.
const RING_SPIN = [0.2, 0.5, 0.8, 0.998]
const RING_INC = [10.0, 30.0, 45.0, 60.0]
const RING_RADIUS = [2.0, 5.0, 10.0, 20.0, 50.0]
const RING_HEIGHT = [2.5, 5.0, 10.0, 20.0]

"""
    ring_profile(spin, inclination, radius, height) -> Vector{Float64}

Unit-area line profile `L(g)` on `G_GRID` for a co-rotating ring of the given
radius and height, in gravitational radii, above a thin disc.
"""
function ring_profile(spin::Real, inclination::Real, radius::Real, height::Real)
    metric = KerrMetric(M = 1.0, a = Float64(spin))
    disc = ThinDisc(0.0, Inf)
    corona = RingCorona(r = Float64(radius), h = Float64(height))
    # 1000 is Gradus's default. The ring solver also requires an even count
    # greater than 160. The previous package used 256.
    emissivity = emissivity_profile(metric, disc, corona; n_samples = 1000)
    clamp_ring_emissivity!(emissivity)
    return unit_area(observed_flux(metric, disc, emissivity, inclination),
                     "ring at spin=$spin, inc=$inclination, r=$radius, h=$height")
end

const LAMP_PROFILES = GridInterpolator(
    (spin = LAMP_SPIN, inc = LAMP_INC, h = LAMP_HEIGHT);
    cache = "jllamp-thin-v2",
) do spin, inc, h
    lamppost_profile(spin, inc, h)
end

const LAMP_PARAMETERS = [
    Parameter("spin", 0.8; min = 0.2, max = 0.998, delta = 0.05),
    Parameter("inc", 30.0; unit = "degrees", min = 10.0, max = 60.0, delta = 5.0),
    Parameter("h", 10.0; unit = "GM/c^2", min = 3.0, max = 20.0, delta = 1.0),
]

"The cached lamppost line profile as a unit-area kernel in g."
function lamppost_kernel(p)
    L = max.(LAMP_PROFILES(p.spin, p.inc, p.h), 0.0)
    return RedshiftKernel(G_GRID, L)
end

const jllamp = ConvolutionModel(
    "jllamp",
    LAMP_PARAMETERS;
    description = "Thin-disc lamppost blur. Profiles are cached on a coarse grid.",
) do spectrum, p
    convolve(spectrum, lamppost_kernel(p); method = :fft)
end

# Kernel parameters first, then the xillver parameters. `inc` is the lamppost
# observer angle; `Incl` is the angle stored in the reflection table.
const jllampxill = Blurred(
    "jllampxill",
    JuliaXSPECModels.xillver_table,
    LAMP_PARAMETERS,
    JuliaXSPECModels.XILLVER_PARAMETERS;
    description = "xillverD-5 blurred by the thin-disc lamppost, on the table energy grid",
    method = :fft,
) do kernel
    lamppost_kernel(kernel)
end

const RING_PROFILES = GridInterpolator(
    (spin = RING_SPIN, inc = RING_INC, r = RING_RADIUS, h = RING_HEIGHT);
    cache = "jlring-thin-v1",
) do spin, inc, r, h
    ring_profile(spin, inc, r, h)
end

const RING_PARAMETERS = [
    Parameter("spin", 0.8; min = 0.2, max = 0.998, delta = 0.05),
    Parameter("inc", 30.0; unit = "degrees", min = 10.0, max = 60.0, delta = 5.0),
    Parameter("r", 5.0; unit = "GM/c^2", min = 2.0, max = 50.0, delta = 1.0),
    Parameter("h", 5.0; unit = "GM/c^2", min = 2.5, max = 20.0, delta = 1.0),
]

"The cached co-rotating ring line profile as a unit-area kernel in g."
function ring_kernel(p)
    L = max.(RING_PROFILES(p.spin, p.inc, p.r, p.h), 0.0)
    return RedshiftKernel(G_GRID, L)
end

const jlring = ConvolutionModel(
    "jlring",
    RING_PARAMETERS;
    description = "Thin-disc ring blur. Profiles are cached on a coarse grid.",
) do spectrum, p
    convolve(spectrum, ring_kernel(p); method = :fft)
end

const jlringxill = Blurred(
    "jlringxill",
    JuliaXSPECModels.xillver_table,
    RING_PARAMETERS,
    JuliaXSPECModels.XILLVER_PARAMETERS;
    description = "xillverD-5 blurred by a thin-disc ring, on the table energy grid",
    method = :fft,
) do kernel
    ring_kernel(kernel)
end

function __init__()
    register!(jllamp)
    register!(jllampxill)
    register!(jlring)
    register!(jlringxill)
end

end
