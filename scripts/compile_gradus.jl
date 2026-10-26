# Compile the reference models and the thin-disc lamppost into one library.
# Run via ./build-julia.sh gradus.
#
# Output, kept apart from the core build:
#   build/gradus/lib/libjuliaxspec_models.{dylib,so}
#   build/gradus/include/julia_init.h
#   build/gradus/xspec/model.dat
#   build/gradus/xspec/juliaxspec_wrappers.c

using PackageCompiler
using JuliaXSPEC
using JuliaXSPECModels
using GradusXSPEC

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "scripts", "stub_tbbmalloc_proxy.jl"))

const DEST = joinpath(ROOT, "build", "gradus")

println("Compiling the reference models and jllamp into build/gradus/")
println("This takes several minutes. The lamppost profile at the default")
println("parameters is calculated once during the build and then remembered.")
create_library(
    joinpath(ROOT, "GradusXSPEC"),
    DEST;
    lib_name = "juliaxspec_models",
    force = true,
    include_transitive_dependencies = true,
    include_lazy_artifacts = true,
    precompile_execution_file = joinpath(ROOT, "scripts", "precompile_gradus.jl"),
)

stub_tbbmalloc_proxy!(DEST)

xspec_dir = joinpath(DEST, "xspec")
mkpath(xspec_dir)
println("Writing $(xspec_dir)/model.dat for: ",
        join([m.name for m in registered_models()], ", "))
write_xspec_package(xspec_dir; julia_init_header = "../include/julia_init.h")
println("Done. Next: ./build-xspec.sh gradus")
