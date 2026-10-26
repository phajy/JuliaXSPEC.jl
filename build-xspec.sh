#!/usr/bin/env bash
# Step 2 of 2. Needs HEASOFT (HEADAS set, initpackage and hmake on the PATH).
#
#   ./build-xspec.sh           the reference models (after ./build-julia.sh)
#   ./build-xspec.sh gradus    the lamppost library (after ./build-julia.sh gradus)
set -euo pipefail
cd "$(dirname "$0")"
root="$(pwd)"

variant="${1:-core}"
case "$variant" in
    core)
        xspec_dir="$root/xspec"
        lib_dir="$root/build/lib"
        patch_prefix=".."
        ;;
    gradus)
        xspec_dir="$root/build/gradus/xspec"
        lib_dir="$root/build/gradus/lib"
        patch_prefix="../../.."
        ;;
    *)
        echo "Usage: ./build-xspec.sh [gradus]" >&2
        exit 1
        ;;
esac

[[ -n "${HEADAS:-}" ]] || { echo "HEADAS is not set: source \$HEADAS/headas-init.sh first" >&2; exit 1; }
[[ -f "$xspec_dir/model.dat" ]] || {
    if [[ "$variant" == "gradus" ]]; then
        echo "$xspec_dir/model.dat not found: run ./build-julia.sh gradus first" >&2
    else
        echo "$xspec_dir/model.dat not found: run ./build-julia.sh first" >&2
    fi
    exit 1
}

cd "$xspec_dir"
# Remove what a previous initpackage/hmake left behind, then regenerate.
rm -f Makefile pkgIndex.tcl lpack_juliaxspec.* juliaxspecFunctionMap.* ./*.o libjuliaxspec.*
initpackage juliaxspec model.dat .
julia "$patch_prefix/scripts/patch_xspec_makefile.jl" Makefile "$lib_dir"
hmake

echo
echo "Built $(ls libjuliaxspec.*)."
echo "In XSPEC, from any directory:"
echo "  lmod juliaxspec $xspec_dir"
