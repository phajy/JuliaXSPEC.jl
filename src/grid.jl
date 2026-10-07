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

With `cache = "some-name"` each newly computed corner is also written under
[`cache_directory`](@ref), so a later session reloads it instead of calling
`f` again. Change the name (or the grid) to invalidate those files. All grid
caches in the process share a RAM budget, [`cache_limit_bytes`](@ref); past it,
the least recently used corner of a cache is dropped. Disk copies are kept.
"""
mutable struct GridInterpolator{F, N}
    f::F
    names::NTuple{N, Symbol}
    axes::NTuple{N, GridAxis}
    corners::Dict{NTuple{N, Int}, Vector{Float64}}
    order::Vector{NTuple{N, Int}}          # least recently used first
    nbytes::Dict{NTuple{N, Int}, UInt64}
    cache_name::Union{Nothing, String}
    disk_dir::String
    hits::Base.RefValue{Int}
    misses::Base.RefValue{Int}
    disk_loads::Base.RefValue{Int}
end

# Every live interpolator, so a shared RAM budget can free someone else's
# oldest corner. Weak, so a cache does not stay in memory just because we know
# about it.
const CACHE_LOCK = ReentrantLock()
const CACHE_BYTES = Ref{UInt64}(0)
const LIVE_CACHES = WeakRef[]

function GridInterpolator(f, axes::NamedTuple; cache::Union{Nothing, AbstractString} = nothing)
    N = length(axes)
    N >= 1 || throw(ArgumentError("GridInterpolator needs at least one axis"))
    if cache !== nothing
        occursin(r"^[A-Za-z0-9._-]+$", cache) ||
            throw(ArgumentError("cache name $(repr(cache)) must contain only letters, digits, '.', '_' and '-'"))
    end
    grid_axes = ntuple(i -> GridAxis(axes[i]), N)
    cache_name = cache === nothing ? nothing : String(cache)
    g = GridInterpolator{typeof(f), N}(
        f, keys(axes), grid_axes,
        Dict{NTuple{N, Int}, Vector{Float64}}(),
        NTuple{N, Int}[],
        Dict{NTuple{N, Int}, UInt64}(),
        cache_name,
        cache_name === nothing ? "" : joinpath(cache_directory(), cache_name),
        Ref(0), Ref(0), Ref(0),
    )
    cache_name === nothing || prepare_disk!(g)
    finalizer(release_cache_memory!, g)
    lock(CACHE_LOCK) do
        filter!(w -> w.value !== nothing, LIVE_CACHES)
        push!(LIVE_CACHES, WeakRef(g))
    end
    return g
end

"Drop this cache's contribution to the shared RAM budget. Used as a finalizer."
function release_cache_memory!(g::GridInterpolator)
    lock(CACHE_LOCK) do
        for n in values(g.nbytes)
            CACHE_BYTES[] -= n
        end
        empty!(g.nbytes)
    end
    return nothing
end

"""
    grid_size(interpolator) -> Int

Total number of grid corners (the most that could ever be computed).
"""
grid_size(g::GridInterpolator) = prod(length.(g.axes))

function (g::GridInterpolator{F, N})(x::Vararg{Real}) where {F, N}
    length(x) == N || throw(ArgumentError("expected $N parameter(s) $(g.names), got $(length(x))"))
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

"The cached value of `f` at one grid corner, computing it on first use."
function corner_value(g::GridInterpolator{F, N}, index::NTuple{N, Int}) where {F, N}
    cached = lock(CACHE_LOCK) do
        cached_corner(g, index)
    end
    cached !== nothing && return cached

    point = ntuple(i -> g.axes[i].values[index[i]], N)
    value = Vector{Float64}(g.f(point...))
    return lock(CACHE_LOCK) do
        # Another caller may have filled this corner while f was running.
        existing = get(g.corners, index, nothing)
        if existing !== nothing
            g.hits[] += 1
            touch!(g, index)
            return existing
        end
        g.misses[] += 1
        save_disk!(g, index, value)
        return remember!(g, index, value)
    end
end

"RAM hit, or a disk reload. `nothing` if the corner still has to be computed. Caller holds `CACHE_LOCK`."
function cached_corner(g::GridInterpolator, index)
    value = get(g.corners, index, nothing)
    if value !== nothing
        g.hits[] += 1
        touch!(g, index)
        return value
    end
    loaded = load_disk(g, index)
    loaded === nothing && return nothing
    g.disk_loads[] += 1
    return remember!(g, index, loaded)
end

"Move `index` to the most recently used end of `g.order`. Caller holds `CACHE_LOCK`."
function touch!(g::GridInterpolator, index)
    i = findfirst(==(index), g.order)
    i === nothing || deleteat!(g.order, i)
    push!(g.order, index)
    return nothing
end

"Store `value` in RAM if it fits under the shared budget. Caller holds `CACHE_LOCK`."
function remember!(g::GridInterpolator, index, value::Vector{Float64})
    nb = UInt64(sizeof(value))
    limit = cache_limit_bytes()
    if limit !== nothing
        while CACHE_BYTES[] + nb > limit
            evict_any!() == 0 && return value
        end
    end
    g.corners[index] = value
    g.nbytes[index] = nb
    push!(g.order, index)
    CACHE_BYTES[] += nb
    return value
end

"Free the least recently used corner of some live cache. Caller holds `CACHE_LOCK`."
function evict_any!()
    filter!(w -> w.value !== nothing, LIVE_CACHES)
    for w in LIVE_CACHES
        other = w.value
        other === nothing && continue
        isempty(other.order) && continue
        return evict_oldest!(other)
    end
    return UInt64(0)
end

function evict_oldest!(g::GridInterpolator)
    key = popfirst!(g.order)
    pop!(g.corners, key)
    nb = pop!(g.nbytes, key)
    CACHE_BYTES[] -= nb
    return nb
end

# -- disk -------------------------------------------------------------------

function grid_signature(g::GridInterpolator)
    lines = String["juliaxspec-grid-v1", join(string.(g.names), ",")]
    for axis in g.axes
        push!(lines, string(axis.scale))
        push!(lines, join(string.(axis.values), ","))
    end
    return join(lines, "\n")
end

"Create the cache directory and forget files written for a different grid."
function prepare_disk!(g::GridInterpolator)
    dir = g.disk_dir
    mkpath(dir)
    signature = grid_signature(g)
    sigfile = joinpath(dir, "grid.txt")
    if isfile(sigfile) && read(sigfile, String) != signature
        for file in readdir(dir)
            endswith(file, ".bin") && rm(joinpath(dir, file); force = true)
        end
    end
    write(sigfile, signature)
    return dir
end

function corner_path(g::GridInterpolator, index)
    return joinpath(g.disk_dir, join(string.(index), "_") * ".bin")
end

function load_disk(g::GridInterpolator, index)
    g.cache_name === nothing && return nothing
    path = corner_path(g, index)
    isfile(path) || return nothing
    try
        open(path, "r") do io
            read(io, 4) == b"JXGC" || return nothing
            read(io, UInt32) == UInt32(1) || return nothing
            n = Int(read(io, UInt64))
            value = Vector{Float64}(undef, n)
            read!(io, value)
            return value
        end
    catch
        return nothing
    end
end

function save_disk!(g::GridInterpolator, index, value::Vector{Float64})
    g.cache_name === nothing && return nothing
    path = corner_path(g, index)
    tmp = path * ".tmp"
    try
        open(tmp, "w") do io
            write(io, b"JXGC")
            write(io, UInt32(1))
            write(io, UInt64(length(value)))
            write(io, value)
        end
        mv(tmp, path; force = true)
    catch
        rm(tmp; force = true)   # the disk cache is a convenience; a full disk must not fail the model
    end
    return nothing
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

Forget the corners held in RAM and reset the counters. Files already written
for a `cache` name are left in place, so the next evaluation reloads them.
"""
function empty_cache!(g::GridInterpolator)
    lock(CACHE_LOCK) do
        for n in values(g.nbytes)
            CACHE_BYTES[] -= n
        end
        empty!(g.corners)
        empty!(g.nbytes)
        empty!(g.order)
        g.hits[] = 0
        g.misses[] = 0
        g.disk_loads[] = 0
    end
    return g
end

"""
    disk_loads(interpolator) -> Int

Corners read back from disk rather than computed.
"""
disk_loads(g::GridInterpolator) = g.disk_loads[]

"""
    cache_memory_used_bytes() -> UInt64

RAM currently held by every grid cache in the process.
"""
cache_memory_used_bytes() = CACHE_BYTES[]
