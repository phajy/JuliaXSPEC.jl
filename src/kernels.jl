# Kernels for convolution in the energy ratio g = E_observed / E_emitted.
# Doppler and gravitational shifts multiply photon energies, so a blurring
# that acts on the ratio g (rather than on an energy difference) is the natural
# description of relativistic smearing. See the "Convolution" documentation page.

"""
    RedshiftKernel(g, L)

A blurring kernel in the energy ratio `g = E_obs / E_em`, tabulated at the
increasing grid points `g` with values `L` and normalised so that
`∫ L(g) dg = 1`.

`L(g) dg` is the fraction of photons emitted at energy `E_em` that are observed
between `g E_em` and `(g + dg) E_em`. Unit area therefore means photon number
is conserved by the convolution. The kernel is zero outside `[g[1], g[end]]`,
and linear between grid points.

The field `Λ` holds the cumulative kernel `Λ(g) = ∫_{g[1]}^{g} L dg'` at the
grid points; see [`cumulative`](@ref).
"""
struct RedshiftKernel
    g::Vector{Float64}
    L::Vector{Float64}
    Λ::Vector{Float64}

    function RedshiftKernel(g::AbstractVector{<:Real}, L::AbstractVector{<:Real})
        length(g) == length(L) || throw(ArgumentError("g and L must have the same length"))
        length(g) >= 2 || throw(ArgumentError("a kernel needs at least two grid points"))
        all(>(0), diff(g)) || throw(ArgumentError("kernel g grid must be strictly increasing"))
        g[1] > 0 || throw(ArgumentError("g must be positive"))
        all(>=(0), L) || throw(ArgumentError("kernel values must be non-negative"))
        area = trapezoid_area(g, L)
        area > 0 || throw(ArgumentError("kernel has zero area"))
        Lnorm = collect(Float64, L) ./ area
        return new(collect(Float64, g), Lnorm, cumulative_trapezoid(g, Lnorm))
    end
end

"Trapezoidal-rule integral of samples `y` on the grid `x`."
function trapezoid_area(x::AbstractVector, y::AbstractVector)
    area = 0.0
    for i in 1:(length(x) - 1)
        area += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return area
end

"Running trapezoidal-rule integral of `y` on `x`, starting from 0 at `x[1]`."
function cumulative_trapezoid(x::AbstractVector, y::AbstractVector)
    Λ = zeros(Float64, length(x))
    for i in 2:length(x)
        Λ[i] = Λ[i-1] + 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    end
    return Λ
end

"""
    (kernel::RedshiftKernel)(g) -> Float64

The kernel value at `g`, interpolating linearly between grid points and zero
outside the grid.
"""
function (k::RedshiftKernel)(g::Real)
    (g < k.g[1] || g > k.g[end]) && return 0.0
    hi = searchsortedfirst(k.g, g)
    hi == 1 && return k.L[1]
    lo = hi - 1
    t = (g - k.g[lo]) / (k.g[hi] - k.g[lo])
    return (1 - t) * k.L[lo] + t * k.L[hi]
end

"""
    cumulative(kernel, g) -> Float64

`Λ(g) = ∫_{g_min}^{g} L(g') dg'`: the fraction of photons observed at a ratio
below `g`. Rises from 0 at the low end of the grid to 1 at the high end.
Because `L` is linear between grid points, `Λ` is quadratic there and is
evaluated exactly.
"""
function cumulative(k::RedshiftKernel, g::Real)
    g <= k.g[1] && return 0.0
    g >= k.g[end] && return 1.0
    hi = searchsortedfirst(k.g, g)
    lo = hi - 1
    x = g - k.g[lo]
    slope = (k.L[hi] - k.L[lo]) / (k.g[hi] - k.g[lo])
    return k.Λ[lo] + k.L[lo] * x + 0.5 * slope * x^2
end

"""
    GaussianKernel(sigma; width = 6, n = 1001) -> RedshiftKernel

A Gaussian in `g` centred on `g = 1` with standard deviation `sigma`,
tabulated on `n` points covering `1 ± width * sigma` (truncated at `g > 0`),
then normalised to unit area.

Applied to a narrow line at `E₀`, it produces a Gaussian line of width
`sigma * E₀`: a convenient analytic test of the convolution machinery.
"""
function GaussianKernel(sigma::Real; width::Real = 6, n::Integer = 1001)
    sigma > 0 || throw(ArgumentError("sigma must be positive"))
    g_lo = max(1 - width * sigma, 1e-3)
    g_hi = 1 + width * sigma
    g = collect(range(g_lo, g_hi, length = n))
    L = @. exp(-0.5 * ((g - 1) / sigma)^2)
    return RedshiftKernel(g, L)
end
