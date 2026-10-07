# Caching an expensive function of a few parameters on a grid, and
# interpolating between the grid points. This is the idea behind XSPEC table
# models, except that here the table is filled in on demand and lives in
# memory: only the grid corners a fit actually visits are ever computed.

"""
    GridAxis(values; scale = :auto)

The grid points of one parameter, together with the coordinate in which to
interpolate along it: `:linear` (in the value itself) or `:log` (in its
logarithm). With `scale = :auto` the spacing of `values` decides: uniform
spacing gives `:linear`, uniform spacing in the logarithm gives `:log`, and
anything else is interpolated linearly between the given values.

Ranges choose their scale by type: `range(a, b, n)` is `:linear` and
`logrange(a, b, n)` is `:log`.
"""
struct GridAxis
    values::Vector{Float64}
    scale::Symbol

    function GridAxis(values::AbstractVector{<:Real}, scale::Symbol)
        length(values) >= 1 || throw(ArgumentError("a grid axis needs at least one point"))
        all(>(0), diff(values)) || throw(ArgumentError("grid axis values must be strictly increasing"))
        scale in (:linear, :log) || throw(ArgumentError("scale must be :linear or :log, got $scale"))
        scale == :log && values[1] <= 0 && throw(ArgumentError("a :log axis needs positive values"))
        return new(collect(Float64, values), scale)
    end
end

GridAxis(values::AbstractVector{<:Real}; scale::Symbol = :auto) =
    GridAxis(values, scale == :auto ? detect_scale(values) : scale)
GridAxis(values::AbstractRange{<:Real}) = GridAxis(values, :linear)
GridAxis(values::Base.LogRange{<:Real}) = GridAxis(values, :log)
GridAxis(axis::GridAxis) = axis

"""
    detect_scale(values) -> :linear or :log

`:log` if the values are uniformly spaced in their logarithm (and not in the
values themselves), `:linear` otherwise.
"""
function detect_scale(values::AbstractVector{<:Real})
    length(values) < 3 && return :linear
    uniform(steps) = maximum(steps) - minimum(steps) <= 1e-6 * abs(sum(steps) / length(steps))
    uniform(diff(values)) && return :linear
    all(>(0), values) && uniform(diff(log.(values))) && return :log
    return :linear
end

Base.length(axis::GridAxis) = length(axis.values)

"The coordinate in which this axis is interpolated."
coordinate(axis::GridAxis, x::Real) = axis.scale == :log ? log(x) : Float64(x)

"""
    bracket(axis, x) -> (lo, hi, t)

Indices of the grid points on either side of `x` and the interpolation weight
`t ∈ [0, 1]` of the upper one, measured in the axis's coordinate. Values outside
the grid are clamped to its ends (so `t` is 0 or 1 there).
"""
function bracket(axis::GridAxis, x::Real)
    v = axis.values
    xc = clamp(Float64(x), v[1], v[end])
    hi = searchsortedfirst(v, xc)
    hi == 1 && return (1, 1, 0.0)
    lo = hi - 1
    t = (coordinate(axis, xc) - coordinate(axis, v[lo])) / (coordinate(axis, v[hi]) - coordinate(axis, v[lo]))
    return (lo, hi, t)
end

"""
    GridInterpolator(f, axes::NamedTuple)

Evaluate `f` only at the corners of a parameter grid, remember those results,
and interpolate multilinearly between them for any other parameter values.

`axes` names the parameters and gives each one's grid points as a range,
`logrange`, vector or [`GridAxis`](@ref). `f` takes the parameter values as
positional arguments in the same order and returns a vector (for example a
spectrum on a fixed internal energy grid, or a kernel on a fixed `g` grid).
The interpolator is called like a function:

```julia
kernel = GridInterpolator((spin = range(0, 0.998, 11), h = logrange(2, 100, 20))) do spin, h
    expensive_line_profile(spin, h)
end
L = kernel(0.5, 7.3)   # interpolates between the 4 surrounding corners, computing any that are new
```

Multilinear interpolation blends the corner vectors; it does not shift
features, so the grid must be fine enough that neighbouring corners look alike.
Parameter values outside the grid are clamped to its edge.
[`evaluate_exact`](@ref) bypasses the grid, and [`cache_stats`](@ref) reports
how often corners were reused.
"""
struct GridInterpolator{F, N}
    f::F
    names::NTuple{N, Symbol}
    axes::NTuple{N, GridAxis}
    corners::Dict{NTuple{N, Int}, Vector{Float64}}
    hits::Base.RefValue{Int}
    misses::Base.RefValue{Int}
end

function GridInterpolator(f, axes::NamedTuple)
    N = length(axes)
    N >= 1 || throw(ArgumentError("GridInterpolator needs at least one axis"))
    grid_axes = ntuple(i -> GridAxis(axes[i]), N)
    return GridInterpolator{typeof(f), N}(
        f, keys(axes), grid_axes, Dict{NTuple{N, Int}, Vector{Float64}}(), Ref(0), Ref(0),
    )
end

"""
    grid_size(interpolator) -> Int

Total number of grid corners (the most that could ever be computed).
"""
grid_size(g::GridInterpolator) = prod(length.(g.axes))

function (g::GridInterpolator{F, N})(x::Vararg{Real, N}) where {F, N}
    brackets = ntuple(i -> bracket(g.axes[i], x[i]), N)
    result = nothing
    # Visit the 2^N corners: on each axis pick the lower (1) or upper (2) grid point.
    for choice in Iterators.product(ntuple(_ -> (1, 2), N)...)
        weight = 1.0
        for i in 1:N
            t = brackets[i][3]
            weight *= choice[i] == 1 ? 1 - t : t
        end
        weight == 0 && continue
        index = ntuple(i -> brackets[i][choice[i]], N)
        value = corner_value(g, index)
        result = result === nothing ? weight .* value : result .+ weight .* value
    end
    return result
end

(g::GridInterpolator{F, N})(x) where {F, N} = throw(ArgumentError("expected $N parameters $(g.names), got 1"))

"The cached value of `f` at one grid corner, computing it on first use."
function corner_value(g::GridInterpolator{F, N}, index::NTuple{N, Int}) where {F, N}
    cached = get(g.corners, index, nothing)
    if cached !== nothing
        g.hits[] += 1
        return cached
    end
    g.misses[] += 1
    point = ntuple(i -> g.axes[i].values[index[i]], N)
    value = Vector{Float64}(g.f(point...))
    g.corners[index] = value
    return value
end

"""
    evaluate_exact(interpolator, x...) -> Vector{Float64}

Call the underlying function directly at `x`, with no grid and no cache. Use
it to measure the interpolation error.
"""
evaluate_exact(g::GridInterpolator{F, N}, x::Vararg{Real, N}) where {F, N} = Vector{Float64}(g.f(x...))

"""
    cache_stats(interpolator) -> (hits, misses, stored, total)

Corner lookups served from the cache, corners computed, corners currently
stored, and the total number of corners in the grid.
"""
cache_stats(g::GridInterpolator) =
    (hits = g.hits[], misses = g.misses[], stored = length(g.corners), total = grid_size(g))

"""
    empty_cache!(interpolator)

Forget all computed corners and reset the counters.
"""
function empty_cache!(g::GridInterpolator)
    empty!(g.corners)
    g.hits[] = 0
    g.misses[] = 0
    return g
end
