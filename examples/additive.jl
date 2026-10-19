# An additive model, evaluated in Julia with no XSPEC session.
#
# From the repository root:
#   julia --project=. examples/additive.jl

using JuliaXSPEC

edges = collect(logrange(0.1, 50.0, 2001))

line = AdditiveModel(
    "example_gauss",
    [
        Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6, delta = 0.05),
        Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 10.0, delta = 0.05),
    ],
) do energy_edges, p
    gaussian_per_bin(energy_edges, p.LineE, p.Sigma)   # photons / cm² / s in each bin
end

per_bin = evaluate(line, edges, [6.4, 0.1])
println("bins: ", length(per_bin))
println("total photons / cm² / s: ", sum(per_bin))
sum(per_bin) ≈ 1.0 || error("a unit Gaussian must contain one photon in total, got $(sum(per_bin))")
println("ok")
