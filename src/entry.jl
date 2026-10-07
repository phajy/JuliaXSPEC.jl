# The single point where XSPEC's C world meets Julia. Every generated C wrapper
# (see registry.jl) calls `juliaxspec_evaluate` with its model name; this
# function turns raw pointers into Julia arrays, runs the model and writes the
# result back into XSPEC's buffer. It never lets an exception escape to C.

"""
    juliaxspec_evaluate(model_name, energy, nflux, params, spectrum, flux, flux_error, init) -> Cint

C-callable entry point with XSPEC's model-function arguments plus the model
name in front:

- `energy`: `nflux + 1` bin edges in keV.
- `params`: the parameter values, in `model.dat` order (no `norm`).
- `flux`: `nflux` values of photons cm⁻² s⁻¹ per bin. For additive and
  multiplicative models this is output only; for convolution models it holds
  the input spectrum on entry and is overwritten with the result.
- `flux_error`, `spectrum` (the spectrum number) and `init` are accepted for
  completeness and currently unused.

Returns 0 on success. On failure the error is printed, `flux` is set to zeros
(additive), ones (multiplicative) or left unchanged (convolution), and 1 is
returned so the problem is visible rather than silently fitted.
"""
Base.@ccallable function juliaxspec_evaluate(
    model_name::Ptr{Cchar},
    energy::Ptr{Cdouble},
    nflux::Cint,
    params::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    flux_error::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    n = Int(nflux)
    out = unsafe_wrap(Array, flux, n)    # XSPEC's buffer: the only memory we write to
    model = nothing
    try
        model = lookup_model(unsafe_string(model_name))
        edges = copy(unsafe_wrap(Array, energy, n + 1))
        values = copy(unsafe_wrap(Array, params, length(model.parameters)))

        started = time()
        if model isa ConvolutionModel
            input = BinnedSpectrum(edges, copy(out))
            out .= evaluate(model, input, values).per_bin
        else
            out .= evaluate(model, edges, values)
        end

        if verbose()
            settings = join(["$(p.name)=$(v)" for (p, v) in zip(model.parameters, values)], ", ")
            elapsed = round((time() - started) * 1e3; digits = 2)
            lo, hi = round(edges[1]; sigdigits = 4), round(edges[end]; sigdigits = 4)
            println("JuliaXSPEC: $(model.name)($settings) on $n bins, $lo-$hi keV, $elapsed ms")
            flush(stdout)
        end
        return Cint(0)
    catch err
        println(stderr, "JuliaXSPEC: error evaluating model $(model === nothing ? unsafe_string(model_name) : model.name):")
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        flush(stderr)
        if model isa AdditiveModel
            fill!(out, 0.0)
        elseif model isa MultiplicativeModel
            fill!(out, 1.0)
        end
        return Cint(1)
    end
end
