# Exercised while building the library so that the code paths XSPEC will use
# are compiled ahead of time and the first model evaluation in XSPEC is fast.

using JuliaXSPEC, JuliaXSPECModels

edges = collect(logrange(0.1, 50.0, 1001))
for model in registered_models()
    values = parameter_defaults(model)
    if model isa ConvolutionModel
        input = BinnedSpectrum(edges, gaussian_per_bin(edges, 6.4, 0.1))
        evaluate(model, input, values)
    else
        evaluate(model, edges, values)
    end

    # The same evaluations through the C entry point, as XSPEC will call it.
    flux = model isa ConvolutionModel ? gaussian_per_bin(edges, 6.4, 0.1) : zeros(length(edges) - 1)
    flux_error = zeros(length(edges) - 1)
    GC.@preserve edges values flux flux_error begin
        JuliaXSPEC.juliaxspec_evaluate(
            Base.unsafe_convert(Ptr{Cchar}, model.name), pointer(edges), Cint(length(edges) - 1),
            pointer(values), Cint(1), pointer(flux), pointer(flux_error), Base.unsafe_convert(Ptr{Cchar}, ""),
        )
    end
end
