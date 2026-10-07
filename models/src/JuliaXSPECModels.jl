"""
    JuliaXSPECModels

The reference model set that is compiled into the XSPEC library.

1. `jlgauss`       — evaluate a model directly (a Gaussian line).
2. `jlgausscached` — the same model, cached on a parameter grid and interpolated.
3. `jlgconv`       — convolution, by the direct bin-by-bin matrix.
4. `jlgconvfft`    — the same convolution, by FFT.
5. `jltable`       — an OGIP table model (`xillverD-5.fits`).
6. `jltableblur`   — that table blurred on its own energy grid.

Add your own models here (or in a package of your own that depends on
JuliaXSPEC), register them in `__init__`, and rebuild.
"""
module JuliaXSPECModels

using JuliaXSPEC

# ---------------------------------------------------------------------------
# 1. Model evaluation
#
# A model function receives the energy bin edges (keV) and the parameters as a
# NamedTuple, and returns the photons / cm² / s in each bin. XSPEC multiplies
# by `norm` itself. These parameter lines mirror XSPEC's own `gaussian`.

const jlgauss = AdditiveModel(
    "jlgauss",
    [
        Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6, delta = 0.05),
        Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 20.0, soft_max = 10.0, delta = 0.05),
    ];
    description = "Gaussian line, computed directly (same as XSPEC's gaussian)",
) do edges, p
    gaussian_per_bin(edges, p.LineE, p.Sigma)
end

# ---------------------------------------------------------------------------
# 2. Caching on a grid
#
# Pretend the Gaussian were expensive. We compute it only at the corners of a
# (LineE, Sigma) grid, on a fixed fine internal energy grid, and interpolate.
# LineE is spaced linearly, Sigma logarithmically; the interpolator notices and
# interpolates each axis in its own coordinate. Comparing jlgausscached with
# jlgauss in XSPEC shows the interpolation error directly.

const INTERNAL_EDGES = collect(logrange(0.1, 100.0, 4001))   # 0.17% wide bins

const CACHED_GAUSSIAN = GridInterpolator(
    (LineE = range(5.0, 8.0, length = 61), Sigma = logrange(0.05, 1.0, 14)),
) do LineE, Sigma
    gaussian_per_bin(INTERNAL_EDGES, LineE, Sigma)
end

const jlgausscached = AdditiveModel(
    "jlgausscached",
    [
        Parameter("LineE", 6.4; unit = "keV", min = 5.0, max = 8.0, delta = 0.05),
        Parameter("Sigma", 0.3; unit = "keV", min = 0.05, max = 1.0, delta = 0.01),
    ];
    description = "The same Gaussian line, interpolated from a cached parameter grid",
) do edges, p
    on_internal_grid = BinnedSpectrum(INTERNAL_EDGES, CACHED_GAUSSIAN(p.LineE, p.Sigma))
    rebin(on_internal_grid, edges).per_bin
end

# ---------------------------------------------------------------------------
# 3. Convolution
#
# A convolution model receives the spectrum of the components it wraps as a
# BinnedSpectrum and returns the blurred BinnedSpectrum on the same bins.
# A Gaussian kernel in g turns a narrow line at E₀ into a Gaussian line of
# width SigmaG * E₀, so `jlgconv * gaussian` can be checked analytically.

const jlgconv = ConvolutionModel(
    "jlgconv",
    [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5, delta = 0.005)];
    description = "Gaussian blur in g = E_obs / E_em",
) do spectrum, p
    convolve(spectrum, GaussianKernel(p.SigmaG); method = :direct)
end

const jlgconvfft = ConvolutionModel(
    "jlgconvfft",
    [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5, delta = 0.005)];
    description = "Gaussian blur in g, evaluated with an FFT",
) do spectrum, p
    convolve(spectrum, GaussianKernel(p.SigmaG); method = :fft)
end

# ---------------------------------------------------------------------------
# 4. Table model, and a blur done on the table's own energy grid
#
# xillverD-5.fits is not part of the repository (it is about 1 GB). These
# parameter limits are the ones in that file, so model.dat can be generated
# without it; the file is read on the first evaluation. Set
# JULIAXSPEC_TABLE_DIR if it does not live in the repository root.

const XILLVER_PARAMETERS = [
    Parameter("Gamma", 2.0; min = 1.2, max = 3.6, delta = 0.01),
    Parameter("A_Fe", 1.0; min = 0.5, max = 20.0, delta = 0.01),
    Parameter("logXi", 2.0; min = 0.0, max = 4.698969841003418, delta = 0.01),
    Parameter("Dens", 17.0; min = 15.0, max = 19.0, delta = 0.1),
    Parameter("Incl", 45.0; unit = "degrees", min = 18.194873809814453, max = 87.13401794433594, delta = 0.01),
]

const XILLVER = Ref{Union{Nothing, OGIPTable}}(nothing)

function xillver_table()
    if XILLVER[] === nothing
        table = OGIPTable(resolve_table_path("xillverD-5.fits"))
        got = [p.name for p in table.parameters]
        expected = [p.name for p in XILLVER_PARAMETERS]
        got == expected || error("xillverD-5.fits parameters $got do not match $expected")
        XILLVER[] = table
    end
    return XILLVER[]
end

const jltable = AdditiveModel(
    "jltable",
    XILLVER_PARAMETERS;
    description = "xillverD-5 reflection table, interpolated in Julia",
) do edges, p
    table = xillver_table()
    on_table_grid = BinnedSpectrum(table.edges, table(p.Gamma, p.A_Fe, p.logXi, p.Dens, p.Incl))
    return rebin(on_table_grid, edges).per_bin
end

# The blur runs from the table's energy grid (0.07–1000 keV) onto XSPEC's
# bins, so photons outside the instrument band can still redshift into it.
# FFT, because the table has 2999 bins and the direct matrix would be large.
const jltableblur = Blurred(
    "jltableblur",
    xillver_table,
    [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5, delta = 0.005)],
    XILLVER_PARAMETERS;
    description = "xillver reflection blurred by a Gaussian in g, on the table energy grid",
    method = :fft,
) do kernel
    GaussianKernel(kernel.SigmaG)
end

# Registration must happen at load time, not at precompile time, so it lives
# in __init__.
function __init__()
    register!(jlgauss)
    register!(jlgausscached)
    register!(jlgconv)
    register!(jlgconvfft)
    register!(jltable)
    register!(jltableblur)
end

end
