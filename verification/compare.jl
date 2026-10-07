# Compare the model values written by phase1.xcm, print a summary, and make
# the figures shown on the Verification page of the documentation.
#
# Run from the verification/ directory, after `xspec - phase1.xcm`:
#   julia --project=. -e 'using Pkg; Pkg.develop(path = ".."); Pkg.instantiate()'   # first time
#   julia --project=. compare.jl
#
# Exit status is non-zero if any comparison exceeds its tolerance.

ENV["GKSwstype"] = "100"     # headless plotting
using JuliaXSPEC
using Plots
using Printf

const OUT = joinpath(@__DIR__, "output")
const FIGS = joinpath(@__DIR__, "..", "docs", "src", "assets")
mkpath(FIGS)

"Read a file written by the `dump` proc in phase1.xcm: bin centre, half-width, density."
function read_model(name)
    rows = [parse.(Float64, split(line)) for line in eachline(joinpath(OUT, name * ".dat")) if !isempty(strip(line))]
    E = [r[1] for r in rows]
    dE = [r[2] for r in rows]
    y = [r[3] for r in rows]
    return (E = E, halfwidth = dE, density = y)
end

const failures = Ref(0)

"""
Largest difference between `test` and `reference` as a fraction of the
reference's peak, over the bins selected by `mask`. Prints one line.
"""
function report(label, test, reference; tol, mask = trues(length(test)))
    peak = maximum(reference[mask])
    worst = maximum(abs.(test[mask] .- reference[mask])) / peak
    ok = worst <= tol
    ok || (failures[] += 1)
    @printf("%-58s  max |Δ| / peak = %9.2e   tolerance %7.1e   %s\n", label, worst, tol, ok ? "ok" : "FAILED")
    return worst
end

function overlay(title, E, curves, ylabel; residual = nothing, legend = :topright)
    top = plot(; title, ylabel, legend)
    for (name, y, style) in curves
        plot!(top, E, y; label = name, style...)
    end
    residual === nothing && return top
    bottom = plot(E, residual.y; xlabel = "Energy (keV)", ylabel = residual.label, label = "")
    hline!(bottom, [0.0]; color = :gray, linestyle = :dash, label = "")
    return plot(top, bottom; layout = grid(2, 1, heights = [0.7, 0.3]), link = :x, size = (800, 600))
end

println("Phase 1 verification\n")

# 1. Direct evaluation --------------------------------------------------------
ref = read_model("gaussian_6.4_0.3")
jl = read_model("jlgauss_6.4_0.3")
report("jlgauss vs XSPEC gaussian (LineE 6.4, Sigma 0.3)", jl.density, ref.density; tol = 1e-6)
window = 5.0 .< ref.E .< 7.8
savefig(
    overlay("jlgauss against XSPEC's gaussian", ref.E[window],
        [("XSPEC gaussian", ref.density[window], (linewidth = 3, color = :black)),
         ("jlgauss", jl.density[window], (linewidth = 1.5, color = :orange, linestyle = :dash))],
        "photons / cm² / s / keV";
        residual = (y = (jl.density .- ref.density)[window] ./ maximum(ref.density), label = "(jlgauss − gaussian) / peak")),
    joinpath(FIGS, "verification_jlgauss.png"))

# 2. Grid interpolation --------------------------------------------------------
for (sigma, tol, lo, hi) in (("0.3", 0.02, 5.0, 7.8), ("0.07", 0.5, 6.0, 6.9))
    direct = read_model("jlgauss_6.43_$sigma")
    cached = read_model("jlgausscached_6.43_$sigma")
    report("jlgausscached vs jlgauss (LineE 6.43, Sigma $sigma)", cached.density, direct.density; tol)
    w = lo .< direct.E .< hi
    savefig(
        overlay("Grid interpolation, Sigma = $sigma keV (LineE 6.43 keV, grid step 0.05 keV)", direct.E[w],
            [("direct (jlgauss)", direct.density[w], (linewidth = 3, color = :black)),
             ("interpolated (jlgausscached)", cached.density[w], (linewidth = 1.5, color = :orange, linestyle = :dash))],
            "photons / cm² / s / keV";
            residual = (y = (cached.density .- direct.density)[w] ./ maximum(direct.density), label = "(interpolated − direct) / peak")),
        joinpath(FIGS, "verification_jlgausscached_$sigma.png"))
end

# 3. Convolution of a narrow line ----------------------------------------------
blurred = read_model("jlgconv_0.05_gaussian_6.4_0.01")
analytic = read_model("gaussian_6.4_0.320156")
report("jlgconv(0.05) * gaussian(6.4, 0.01) vs gaussian(6.4, 0.3202)", blurred.density, analytic.density; tol = 5e-3)
window = 5.0 .< analytic.E .< 7.8
fig3 = overlay("jlgconv * gaussian against the analytic result", analytic.E[window],
    [("gaussian(6.4, 0.3202) — expected", analytic.density[window], (linewidth = 3, color = :black)),
     ("jlgconv(SigmaG 0.05) * gaussian(6.4, 0.01)", blurred.density[window], (linewidth = 1.5, color = :orange, linestyle = :dash))],
    "photons / cm² / s / keV";
    residual = (y = (blurred.density .- analytic.density)[window] ./ maximum(analytic.density), label = "(blurred − expected) / peak"))
savefig(fig3, joinpath(FIGS, "verification_jlgconv_line.png"))

# 4. Convolution of a power law --------------------------------------------------
Γ = 2.5
kernel = GaussianKernel(0.05)
factor = JuliaXSPEC.trapezoid_area(kernel.g, kernel.L .* kernel.g .^ (Γ - 1))
pl = read_model("powerlaw_2.5")
blurred_pl = read_model("jlgconv_0.05_powerlaw_2.5")
interior = (pl.E .> 0.1 / kernel.g[1]) .& (pl.E .< 50.0 / kernel.g[end])
report(@sprintf("jlgconv(0.05) * powerlaw(2.5) vs %.5f * powerlaw(2.5), interior", factor),
       blurred_pl.density, factor .* pl.density; tol = 2e-3, mask = interior)
ratio = blurred_pl.density ./ pl.density
fig4 = plot(pl.E, ratio; xscale = :log10, label = "jlgconv * powerlaw / powerlaw", linewidth = 2,
            xlabel = "Energy (keV)", ylabel = "ratio", title = "A blurred power law keeps its slope", size = (800, 420),
            ylims = (0.95, 1.05), bottom_margin = 5Plots.mm, left_margin = 3Plots.mm)
hline!(fig4, [factor]; color = :black, linestyle = :dash, label = @sprintf("∫ L(g) g^(Γ−1) dg = %.5f", factor))
vspan!(fig4, [pl.E[1], 0.1 / kernel.g[1]]; color = :gray, alpha = 0.2, label = "band edges: photons leave the grid")
vspan!(fig4, [50.0 / kernel.g[end], pl.E[end]]; color = :gray, alpha = 0.2, label = "")
savefig(fig4, joinpath(FIGS, "verification_jlgconv_powerlaw.png"))

println("\nFigures written to ", abspath(FIGS))
if failures[] > 0
    println("$(failures[]) comparison(s) exceeded tolerance")
    exit(1)
end
println("All comparisons within tolerance.")
