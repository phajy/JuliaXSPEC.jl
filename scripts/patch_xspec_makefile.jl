# XSPEC's initpackage writes a Makefile that knows nothing about our compiled
# Julia library. This adds the link flags for it, right after the -lXS entry
# of HD_SHLIB_LIBS. Safe to run more than once.
#
#   julia scripts/patch_xspec_makefile.jl [xspec/Makefile]

const ROOT = dirname(@__DIR__)
const LIB_DIR = joinpath(ROOT, "build", "lib")
const MARKER = "-ljuliaxspec_models"

makefile = isempty(ARGS) ? joinpath(ROOT, "xspec", "Makefile") : ARGS[1]
isfile(makefile) || error("$makefile not found: run initpackage first (see build-xspec.sh)")

lines = readlines(makefile)
if any(contains(MARKER), lines)
    println("$makefile already links $MARKER")
    exit(0)
end

anchor = findfirst(l -> contains(l, "-lXS") && endswith(rstrip(l), "\\"), lines)
anchor === nothing && error("could not find the HD_SHLIB_LIBS line containing -lXS in $makefile")

flags = "-L$LIB_DIR $MARKER -Wl,-rpath,$LIB_DIR -Wl,-rpath,$(joinpath(LIB_DIR, "julia"))"
insert!(lines, anchor + 1, "\t\t\t  $flags \\")
write(makefile, join(lines, "\n") * "\n")
println("Added $flags to $makefile")
