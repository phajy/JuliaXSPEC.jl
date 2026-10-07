# Runtime options. XSPEC's model.dat passes a single-token init string per
# model, which is too little for configuration, so JuliaXSPEC reads
# environment variables instead. Set them before starting xspec.

"""
    verbose() -> Bool

True when `JULIAXSPEC_VERBOSE` is set to `1`, `true` or `yes`. The entry point
then prints one line per model evaluation (model, parameters, bins, timing).
"""
verbose() = lowercase(get(ENV, "JULIAXSPEC_VERBOSE", "0")) in ("1", "true", "yes")

"""
    convolution_method() -> Symbol

`:direct` (the default) or `:fft`, from `JULIAXSPEC_CONVOLVE`.
`matrix` is accepted as a synonym of `direct`, and `fourier` of `fft`.
"""
function convolution_method()
    raw = lowercase(strip(get(ENV, "JULIAXSPEC_CONVOLVE", "direct")))
    if raw in ("direct", "matrix", "")
        return :direct
    elseif raw in ("fft", "fourier")
        return :fft
    else
        throw(ArgumentError("JULIAXSPEC_CONVOLVE must be 'direct' or 'fft' (got $(repr(raw)))"))
    end
end

"""
    cache_limit_bytes() -> Union{UInt64,Nothing}

RAM budget for cached grid corners, from `JULIAXSPEC_CACHE_LIMIT_GB`
(default 16). `nothing` means unlimited (`0` or a negative value).
"""
function cache_limit_bytes()
    raw = get(ENV, "JULIAXSPEC_CACHE_LIMIT_GB", "16")
    gb = tryparse(Float64, raw)
    gb === nothing && (gb = 16.0)
    gb <= 0 && return nothing
    return UInt64(round(gb * (UInt64(1) << 30)))
end

"""
    cache_directory() -> String

Where grid corners are written when a [`GridInterpolator`](@ref) is given a
`cache` name. `JULIAXSPEC_CACHE_DIR` if set, otherwise
`~/.julia/juliaxspec`.
"""
function cache_directory()
    explicit = get(ENV, "JULIAXSPEC_CACHE_DIR", "")
    return isempty(explicit) ? joinpath(first(DEPOT_PATH), "juliaxspec") : explicit
end
