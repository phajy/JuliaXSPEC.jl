# An additive model that owns both the source spectrum and the kernel, so the
# blur happens on the source's energy grid rather than on whatever bins XSPEC
# happens to be using. Photons outside the instrument band can still redshift
# into it. See step 9 of the Convolution documentation page.

"""
    Blurred(kernel_fn, name, table_loader, kernel_parameters, source_parameters;
            description = "", method = :direct) -> AdditiveModel

Build an additive model that interpolates a table, blurs it, and rebins onto
the bins XSPEC asked for.

`table_loader()` returns an [`OGIPTable`](@ref) (called on each evaluation, so
the loader can cache the file). `kernel_fn` receives the kernel parameters as
a `NamedTuple` and returns a [`RedshiftKernel`](@ref). The model's parameters
are the kernel parameters followed by the source parameters.

The blur uses [`convolve`](@ref) with the table's own energy bins as the
*input* and XSPEC's bins as the output, so the whole table contributes.
`method` is `:direct` or `:fft`.

```julia
Blurred("jltableblur", xillver_table, [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5)],
        xillver_parameters) do kernel
    GaussianKernel(kernel.SigmaG)
end
```
"""
function Blurred(
    kernel_fn,
    name::AbstractString,
    table_loader,
    kernel_parameters::AbstractVector{Parameter},
    source_parameters::AbstractVector{Parameter};
    description::AbstractString = "",
    method::Symbol = :direct,
)
    nker = length(kernel_parameters)
    ker_names = ntuple(i -> Symbol(kernel_parameters[i].name), nker)
    parameters = vcat(collect(kernel_parameters), collect(source_parameters))
    return AdditiveModel(name, parameters; description) do edges, p
        vals = values(p)
        kernel = kernel_fn(NamedTuple{ker_names}(vals[1:nker]))
        table = table_loader()
        source = BinnedSpectrum(table.edges, table(vals[nker+1:end]...))
        return convolve(source, kernel; out_edges = edges, method).per_bin
    end
end
