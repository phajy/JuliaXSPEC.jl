#!/usr/bin/env bash
# Step 2 of 2: build the XSPEC local-model package in xspec/ from the files
# written by build-julia.sh. Needs HEASOFT (HEADAS set, initpackage and hmake
# on the PATH).
set -euo pipefail
cd "$(dirname "$0")/xspec"

[[ -n "${HEADAS:-}" ]] || { echo "HEADAS is not set: source \$HEADAS/headas-init.sh first" >&2; exit 1; }
[[ -f model.dat ]] || { echo "xspec/model.dat not found: run ./build-julia.sh first" >&2; exit 1; }

# Remove what a previous initpackage/hmake left behind, then regenerate.
rm -f Makefile pkgIndex.tcl lpack_juliaxspec.* juliaxspecFunctionMap.* ./*.o libjuliaxspec.*
initpackage juliaxspec model.dat .
julia ../scripts/patch_xspec_makefile.jl Makefile
hmake

echo
echo "Built $(ls libjuliaxspec.*). In XSPEC:  lmod juliaxspec $(pwd)"
