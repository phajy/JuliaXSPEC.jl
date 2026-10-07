# Compile the model set in models/ into a shared library that XSPEC can load,
# and write the XSPEC package files (model.dat and the C wrappers).
#
# Run via ./build-julia.sh, or directly:
#   julia --project=scripts scripts/compile.jl
#
# Output:
#   build/lib/libjuliaxspec_models.{dylib,so}   the Julia runtime plus our models
#   build/include/julia_init.h                  init_julia()/shutdown_julia() for the C wrappers
#   xspec/model.dat, xspec/juliaxspec_wrappers.c

using PackageCompiler
using JuliaXSPEC, JuliaXSPECModels

const ROOT = dirname(@__DIR__)

println("Compiling models/ into build/ (this takes several minutes) ...")
create_library(
    joinpath(ROOT, "models"),
    joinpath(ROOT, "build");
    lib_name = "juliaxspec_models",
    force = true,
    include_transitive_dependencies = true,
    precompile_execution_file = joinpath(ROOT, "scripts", "precompile_models.jl"),
)

println("Writing xspec/model.dat and xspec/juliaxspec_wrappers.c for: ",
        join([m.name for m in registered_models()], ", "))
write_xspec_package(joinpath(ROOT, "xspec"); julia_init_header = "../build/include/julia_init.h")
println("Done. Next: ./build-xspec.sh")
