# Runtime options. XSPEC's model.dat passes a single-token init string per
# model, which is too little for configuration, so JuliaXSPEC reads
# environment variables instead. Set them before starting xspec.

"""
    verbose() -> Bool

True when `JULIAXSPEC_VERBOSE` is set to `1`, `true` or `yes`. The entry point
then prints one line per model evaluation (model, parameters, bins, timing).
"""
verbose() = lowercase(get(ENV, "JULIAXSPEC_VERBOSE", "0")) in ("1", "true", "yes")
