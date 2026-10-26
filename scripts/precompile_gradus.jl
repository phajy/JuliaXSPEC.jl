# Same warm-up as the core library, plus one lamppost profile at the default
# parameters (a grid point, so this is a single Gradus calculation).

using JuliaXSPEC, JuliaXSPECModels, GradusXSPEC

include(joinpath(@__DIR__, "precompile_models.jl"))
