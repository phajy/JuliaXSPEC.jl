# Cache an expensive spectrum on a parameter grid, then move it onto
# whatever bins the caller asked for.
#
# From the repository root:
#   julia --project=. examples/grid.jl

using JuliaXSPEC

# The expensive function is evaluated only at these corners, on this fixed grid.
const INTERNAL = collect(logrange(0.1, 50.0, 501))
calls = Ref(0)
cached = GridInterpolator((
    LineE = range(5.0, 8.0, length = 7),
    Sigma = logrange(0.05, 1.0, 5),
)) do LineE, Sigma
    calls[] += 1
    gaussian_per_bin(INTERNAL, LineE, Sigma)
end

on_internal_grid = cached(6.4, 0.3)          # four corners the first time
again = cached(6.4, 0.3)                     # the same four, from memory
xspec_edges = collect(logrange(1.0, 10.0, 201))
delivered = rebin(BinnedSpectrum(INTERNAL, on_internal_grid), xspec_edges)

println("corners computed: ", calls[])
println("cache: ", cache_stats(cached))
println("photons delivered: ", sum(delivered.per_bin))

calls[] == 4 || error("expected 4 new corners, computed $(calls[])")
cache_stats(cached).hits == 4 || error("the second evaluation should have reused all four corners")
on_internal_grid ≈ again || error("the cached spectrum changed between calls")
abs(sum(delivered.per_bin) - sum(on_internal_grid)) < 1e-12 ||
    error("rebinning onto a covering grid must conserve photons")
println("ok")
