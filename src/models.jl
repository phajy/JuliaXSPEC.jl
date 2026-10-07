# The three kinds of XSPEC model a user can write in Julia. Each is a plain
# struct holding the XSPEC name, the parameter list and the Julia function that
# does the work. Nothing here knows about C or XSPEC's calling convention; that
# is handled once, in entry.jl.

"""
    AbstractXSPECModel

Supertype of [`AdditiveModel`](@ref), [`MultiplicativeModel`](@ref) and
[`ConvolutionModel`](@ref).
"""
abstract type AbstractXSPECModel end

"""
    AdditiveModel(compute, name, parameters; description = "")

An XSPEC additive model (`add`): an emission component.

`compute(edges, p)` receives the `n + 1` energy bin edges in keV and the
parameters as a `NamedTuple` (so `p.LineE`, `p.Sigma`, ...), and must return
`n` values: the photons cm⁻² s⁻¹ **in each bin** for unit normalisation. XSPEC
multiplies the result by its own `norm` parameter, which is therefore not in
`parameters`.

The do-block form reads naturally:

```julia
jlgauss = AdditiveModel("jlgauss", [Parameter("LineE", 6.4; unit = "keV", min = 0, max = 1e6),
                                    Parameter("Sigma", 0.1; unit = "keV", min = 0, max = 10)]) do edges, p
    gaussian_per_bin(edges, p.LineE, p.Sigma)
end
```
"""
struct AdditiveModel{F} <: AbstractXSPECModel
    name::String
    parameters::Vector{Parameter}
    compute::F
    description::String
end

"""
    MultiplicativeModel(compute, name, parameters; description = "")

An XSPEC multiplicative model (`mul`), such as an absorber. `compute(edges, p)`
returns the dimensionless factor applied to each of the `n` bins.
"""
struct MultiplicativeModel{F} <: AbstractXSPECModel
    name::String
    parameters::Vector{Parameter}
    compute::F
    description::String
end

"""
    ConvolutionModel(compute, name, parameters; description = "")

An XSPEC convolution model (`con`), which transforms the spectrum of the
components it wraps, e.g. `jlgconv * gaussian` in XSPEC.

`compute(spectrum, p)` receives the input as a [`BinnedSpectrum`](@ref)
(photons per bin, on XSPEC's current energy bins) and must return a
`BinnedSpectrum` on the same bins. See [`convolve`](@ref) for the usual
implementation.
"""
struct ConvolutionModel{F} <: AbstractXSPECModel
    name::String
    parameters::Vector{Parameter}
    compute::F
    description::String
end

for T in (:AdditiveModel, :MultiplicativeModel, :ConvolutionModel)
    @eval function $T(compute, name::AbstractString, parameters::AbstractVector{Parameter}; description::AbstractString = "")
        occursin(r"^[A-Za-z][A-Za-z0-9_]*$", name) ||
            throw(ArgumentError("model name $(repr(name)) must be letters, digits and underscores, starting with a letter"))
        names = [p.name for p in parameters]
        allunique(names) || throw(ArgumentError("model $name has duplicate parameter names"))
        return $T(String(name), collect(parameters), compute, String(description))
    end
end

"""
    xspec_type(model) -> String

The model type keyword used in `model.dat`: `"add"`, `"mul"` or `"con"`.
"""
xspec_type(::AdditiveModel) = "add"
xspec_type(::MultiplicativeModel) = "mul"
xspec_type(::ConvolutionModel) = "con"

"""
    parameter_names(model) -> Vector{Symbol}
"""
parameter_names(m::AbstractXSPECModel) = [Symbol(p.name) for p in m.parameters]

"""
    parameter_defaults(model) -> Vector{Float64}
"""
parameter_defaults(m::AbstractXSPECModel) = [p.default for p in m.parameters]

"""
    named_parameters(model, values) -> NamedTuple

Pair the parameter values XSPEC passes (in `model.dat` order) with their names,
so model code can write `p.LineE` instead of `values[1]`.
"""
function named_parameters(m::AbstractXSPECModel, values::AbstractVector{<:Real})
    length(values) == length(m.parameters) ||
        throw(ArgumentError("model $(m.name) expects $(length(m.parameters)) parameters, got $(length(values))"))
    return NamedTuple{Tuple(parameter_names(m))}(Tuple(Float64.(values)))
end

"""
    evaluate(model::AdditiveModel, edges, values) -> Vector{Float64}
    evaluate(model::MultiplicativeModel, edges, values) -> Vector{Float64}
    evaluate(model::ConvolutionModel, spectrum, values) -> BinnedSpectrum

Run a model exactly as the XSPEC entry point does, but from Julia: `values`
are the parameter values in `model.dat` order. Handy for testing and plotting.
"""
function evaluate(m::Union{AdditiveModel, MultiplicativeModel}, edges::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    result = m.compute(edges, named_parameters(m, values))
    length(result) == length(edges) - 1 ||
        error("model $(m.name) returned $(length(result)) values for $(length(edges) - 1) bins")
    return Vector{Float64}(result)
end

function evaluate(m::ConvolutionModel, spectrum::BinnedSpectrum, values::AbstractVector{<:Real})
    result = m.compute(spectrum, named_parameters(m, values))
    result isa BinnedSpectrum ||
        error("convolution model $(m.name) must return a BinnedSpectrum, got $(typeof(result))")
    nbins(result) == nbins(spectrum) ||
        error("convolution model $(m.name) changed the number of bins from $(nbins(spectrum)) to $(nbins(result))")
    return result
end
