# OGIP / XSPEC additive table models (OGIP memo 92-009).
#
# A table is a grid of spectra, one per corner of a parameter grid, stored as
# photons / cm² / s in each energy bin — the same `per_bin` convention as
# BinnedSpectrum. Reading one is therefore "a GridInterpolator whose corners
# were computed in advance". Each parameter's METHOD flag chooses the
# interpolation coordinate: 0 linear, 1 logarithmic.

using FITSIO

"""
    OGIPTable

An additive XSPEC table model loaded from a FITS file.

`edges` is the table's own energy grid (keV). Calling the table with one
number per parameter interpolates a spectrum on that grid, in photons per bin:

```julia
table = OGIPTable("xillverD-5.fits")
per_bin = table(2.0, 1.0, 2.0, 17.0, 45.0)     # length(edges) - 1 values
```

Interpolation uses each parameter's METHOD flag from the file. Corners are
read from the file, so nothing is recomputed; the [`GridInterpolator`](@ref)
only blends them.
"""
struct OGIPTable{NA, NI, F}
    name::String
    path::String
    parameters::Vector{Parameter}
    edges::Vector{Float64}
    # spectra[energy, i_last, …, i_first]: the last parameter varies fastest,
    # which is the order of the rows in the SPECTRA extension.
    spectra::Array{Float32, NA}
    interpolator::GridInterpolator{F, NI}
end

function OGIPTable(path::AbstractString)
    name, parameters, axes, edges, spectra = read_ogip_table(path)
    names = ntuple(i -> Symbol(parameters[i].name), length(parameters))
    axis_tuple = NamedTuple{names}(ntuple(i -> axes[i], length(axes)))
    interpolator = GridInterpolator(axis_tuple) do values...
        spectrum_at(spectra, values, axes)
    end
    return OGIPTable(String(name), abspath(path), parameters, edges, spectra, interpolator)
end

(table::OGIPTable)(x::Real...) = table.interpolator(x...)

"""
    resolve_table_path(filename) -> String

Find a table file. `JULIAXSPEC_TABLE_DIR` is tried first, then the repository
root (the directory above `src/`), then the current directory.
"""
function resolve_table_path(filename::AbstractString)
    isabspath(filename) && isfile(filename) && return filename
    candidates = String[]
    override = get(ENV, "JULIAXSPEC_TABLE_DIR", "")
    isempty(override) || push!(candidates, joinpath(override, filename))
    push!(candidates, joinpath(@__DIR__, "..", filename))
    push!(candidates, joinpath(pwd(), filename))
    for candidate in candidates
        isfile(candidate) && return abspath(candidate)
    end
    looked = join(candidates, "\n  ")
    error("table file $(repr(filename)) not found. Looked in:\n  $looked")
end

# Row `row` (1-based) of an OGIP SPECTRA extension, last parameter fastest.
function indices_of_row(row::Int, sizes::Tuple)
    idx = Vector{Int}(undef, length(sizes))
    remaining = row - 1
    for k in length(sizes):-1:1
        idx[k] = remaining % sizes[k] + 1
        remaining ÷= sizes[k]
    end
    return Tuple(idx)
end

function row_of_indices(idx::Tuple, sizes::Tuple)
    row = 1
    stride = 1
    for k in length(idx):-1:1
        row += (idx[k] - 1) * stride
        stride *= sizes[k]
    end
    return row
end

"Spectrum at the grid corner whose parameter values are exactly `values`."
function spectrum_at(spectra::AbstractArray, values::Tuple, axes)
    idx = ntuple(length(axes)) do k
        i = searchsortedfirst(axes[k].values, values[k])
        axes[k].values[i] == values[k] ||
            error("internal error: $(values[k]) is not a grid value of parameter $k")
        return i
    end
    # Array layout is (energy, last parameter, ..., first parameter).
    return Float64.(vec(view(spectra, :, reverse(idx)...)))
end

function read_ogip_table(path::AbstractString)
    f = FITS(String(path))
    try
        name = try
            String(strip(read_key(f[1], "MODLNAME")[1]))
        catch
            splitext(basename(path))[1]
        end
        parameters, axes = read_parameter_table(f[2])
        edges = read_energy_edges(f[3])
        spectra = read_spectra(f[4], axes, length(edges) - 1)
        return name, parameters, axes, edges, spectra
    finally
        close(f)
    end
end

function read_parameter_table(hdu)
    names = String.(strip.(read(hdu, "NAME")))
    methods = Vector{Int}(read(hdu, "METHOD"))
    numb = Vector{Int}(read(hdu, "NUMBVALS"))
    initial = Float64.(vec(read(hdu, "INITIAL")))
    delta = Float64.(vec(read(hdu, "DELTA")))
    minimum = Float64.(vec(read(hdu, "MINIMUM")))
    bottom = Float64.(vec(read(hdu, "BOTTOM")))
    top = Float64.(vec(read(hdu, "TOP")))
    maximum = Float64.(vec(read(hdu, "MAXIMUM")))
    raw = read(hdu, "VALUE")

    n = length(names)
    parameters = Parameter[]
    axes = GridAxis[]
    for i in 1:n
        column = parameter_value_column(raw, i, n)
        values = Float64.(column[1:numb[i]])
        scale = if methods[i] == 0
            :linear
        elseif methods[i] == 1
            :log
        else
            error("parameter $(names[i]) has METHOD $(methods[i]); expected 0 (linear) or 1 (log)")
        end
        push!(axes, GridAxis(values, scale))
        # Float32 round-trip in the file can put INITIAL a hair outside the
        # recorded limits. Clamp so the description still satisfies Parameter.
        lo, hi = min(minimum[i], bottom[i]), max(maximum[i], top[i])
        soft_lo, soft_hi = min(bottom[i], top[i]), max(bottom[i], top[i])
        default = clamp(initial[i], soft_lo, soft_hi)
        push!(parameters, Parameter(
            names[i], default;
            min = min(lo, default), max = max(hi, default),
            soft_min = min(soft_lo, default), soft_max = max(soft_hi, default),
            delta = abs(delta[i]) == 0 ? (hi - lo) / 100 : abs(delta[i]),
            frozen = delta[i] < 0,
        ))
    end
    return parameters, axes
end

function parameter_value_column(raw::AbstractMatrix, i::Int, nparams::Int)
    if size(raw, 2) == nparams
        return vec(raw[:, i])
    elseif size(raw, 1) == nparams
        return vec(raw[i, :])
    else
        error("VALUE column has size $(size(raw)); expected one column per parameter ($nparams)")
    end
end

function read_energy_edges(hdu)
    elo = Float64.(vec(read(hdu, "ENERG_LO")))
    ehi = Float64.(vec(read(hdu, "ENERG_HI")))
    length(elo) == length(ehi) || error("ENERG_LO and ENERG_HI have different lengths")
    if length(elo) > 1
        gap = maximum(abs.(elo[2:end] .- ehi[1:end-1]))
        width = minimum(ehi .- elo)
        gap <= 1e-3 * width ||
            error("table energy bins do not meet: largest gap $gap keV, smallest bin $width keV")
    end
    return [elo; ehi[end]]
end

function read_spectra(hdu, axes, n_energy::Int)
    sizes = Tuple(length.(axes))
    nrows = prod(sizes)
    raw = read(hdu, "INTPSPEC")
    matrix = if size(raw) == (n_energy, nrows)
        raw
    elseif size(raw) == (nrows, n_energy)
        permutedims(raw)
    else
        error("INTPSPEC has size $(size(raw)); expected ($n_energy, $nrows) photons per bin")
    end
    check_row_order!(hdu, axes, nrows)
    # (energy, last parameter, ..., first parameter), so one spectrum is contiguous.
    return reshape(matrix, n_energy, reverse(sizes)...)
end

"The SPECTRA rows must list the last parameter varying fastest. Check the ends."
function check_row_order!(hdu, axes, nrows::Int)
    pv = read(hdu, "PARAMVAL")
    n = length(axes)
    row_values = if size(pv, 1) == n
        r -> Float64.(pv[:, r])
    elseif size(pv, 2) == n
        r -> Float64.(pv[r, :])
    else
        error("PARAMVAL has size $(size(pv)); expected $n parameters")
    end
    for row in (1, nrows)
        expected = [axes[k].values[indices_of_row(row, Tuple(length.(axes)))[k]] for k in 1:n]
        actual = row_values(row)
        all(isapprox.(actual, expected; rtol = 1e-5, atol = 1e-5)) ||
            error("SPECTRA row $row has parameters $actual, but OGIP order (last parameter fastest) expects $expected")
    end
    return nothing
end

"""
    write_ogip_table(path, spectra; name, parameters, axes, methods, edges)

Write an additive OGIP table model. `spectra` has shape
`(n_energy, n_last_parameter, …, n_first_parameter)` and holds photons per bin.
`methods[i]` is 0 (linear) or 1 (logarithmic), as in the FITS METHOD column.
`axes[i]` is the grid of parameter `i`. This is the inverse of [`OGIPTable`](@ref)
and is what the tests use to build a small table without a real reflection file.
"""
function write_ogip_table(
    path::AbstractString,
    spectra::AbstractArray{<:Real};
    name::AbstractString,
    parameters::AbstractVector{Parameter},
    axes::AbstractVector,
    methods::AbstractVector{<:Integer},
    edges::AbstractVector{<:Real},
)
    n = length(parameters)
    n == length(axes) == length(methods) ||
        throw(ArgumentError("parameters, axes and methods must have the same length"))
    sizes = Tuple(length.(axes))
    ndims(spectra) == n + 1 || throw(ArgumentError("spectra must have one dimension per parameter, plus energy"))
    size(spectra, 1) == length(edges) - 1 ||
        throw(ArgumentError("spectra has $(size(spectra, 1)) energy bins, edges describe $(length(edges) - 1)"))
    size(spectra)[2:end] == reverse(sizes) ||
        throw(ArgumentError("spectra dimensions $(size(spectra)) do not match the parameter grids $sizes"))

    n_energy = length(edges) - 1
    nrows = prod(sizes)
    maxn = maximum(length.(axes))
    value = zeros(Float32, maxn, n)
    for i in 1:n
        value[1:length(axes[i]), i] .= Float32.(axes[i])
    end
    paramval = zeros(Float32, n, nrows)
    intpspec = zeros(Float32, n_energy, nrows)
    for idx in CartesianIndices(sizes)
        row = row_of_indices(Tuple(idx), sizes)
        for k in 1:n
            paramval[k, row] = axes[k][idx[k]]
        end
        intpspec[:, row] .= spectra[:, reverse(Tuple(idx))...]
    end

    f = FITS(String(path), "w")
    try
        # The first table write creates the primary HDU; the keywords go on it
        # afterwards. REDSHIFT is false so XSPEC does not add its own redshift
        # parameter — JuliaXSPEC models list every parameter themselves.
        write(f,
            ["NAME", "METHOD", "INITIAL", "DELTA", "MINIMUM", "BOTTOM", "TOP", "MAXIMUM", "NUMBVALS", "VALUE"],
            [
                [p.name for p in parameters],
                Int32.(methods),
                Float32.([p.default for p in parameters]),
                Float32.([p.delta for p in parameters]),
                Float32.([p.hard_min for p in parameters]),
                Float32.([p.soft_min for p in parameters]),
                Float32.([p.soft_max for p in parameters]),
                Float32.([p.hard_max for p in parameters]),
                Int32.(length.(axes)),
                value,
            ];
            name = "PARAMETERS",
        )
        write(f, ["ENERG_LO", "ENERG_HI"], [Float32.(edges[1:end-1]), Float32.(edges[2:end])]; name = "ENERGIES")
        write(f, ["PARAMVAL", "INTPSPEC"], [paramval, intpspec];
            name = "SPECTRA", units = Dict("INTPSPEC" => "photons/cm^2/s"))
        write_key(f[1], "MODLNAME", String(name))
        write_key(f[1], "HDUCLAS1", "XSPEC TABLE MODEL")
        write_key(f[1], "ADDMODEL", true)
        write_key(f[1], "REDSHIFT", false)
    finally
        close(f)
    end
    return path
end
