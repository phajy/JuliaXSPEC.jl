#!/usr/bin/env bash
# Step 1 of 2: compile the Julia model set into build/ and write the XSPEC
# package files into xspec/. Needs Julia only (no HEASOFT).
set -euo pipefail
cd "$(dirname "$0")"

# The models project uses this checkout of JuliaXSPEC; the scripts project
# additionally has PackageCompiler. Both develop/instantiate steps are idempotent.
julia --project=models -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
julia --project=scripts -e 'using Pkg; Pkg.develop([PackageSpec(path = "."), PackageSpec(path = "models")]); Pkg.instantiate()'

julia --project=scripts scripts/compile.jl
