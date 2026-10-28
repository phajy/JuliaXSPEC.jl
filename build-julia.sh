#!/usr/bin/env bash
# Step 1 of 2. Needs Julia only (no HEASOFT).
#
#   ./build-julia.sh           reference models, into build/ and xspec/
#   ./build-julia.sh gradus    reference models plus the thin-disc lamppost,
#                              into build/gradus/
set -euo pipefail
cd "$(dirname "$0")"

variant="${1:-core}"

add_astro_registry() {
    julia --project="$1" -e 'using Pkg; Pkg.Registry.add(url = "https://github.com/astro-group-bristol/AstroRegistry")'
}

case "$variant" in
    core)
        # The models project uses this checkout of JuliaXSPEC; the scripts project
        # additionally has PackageCompiler. Both develop/instantiate steps are idempotent.
        julia --project=models -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
        julia --project=scripts -e 'using Pkg; Pkg.develop([PackageSpec(path = "."), PackageSpec(path = "models")]); Pkg.instantiate()'
        julia --project=scripts scripts/compile.jl
        ;;
    gradus)
        add_astro_registry GradusXSPEC
        julia --project=GradusXSPEC -e 'using Pkg; Pkg.instantiate()'
        add_astro_registry builds/gradus
        julia --project=builds/gradus -e 'using Pkg; Pkg.instantiate()'
        julia --project=builds/gradus scripts/compile_gradus.jl
        ;;
    *)
        echo "Usage: ./build-julia.sh [gradus]" >&2
        echo "  ./build-julia.sh          the reference models" >&2
        echo "  ./build-julia.sh gradus   the reference models plus jllamp" >&2
        exit 1
        ;;
esac
