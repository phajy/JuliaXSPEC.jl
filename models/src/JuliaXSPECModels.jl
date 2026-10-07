"""
    JuliaXSPECModels

The reference model set that is compiled into the XSPEC library. It shows the
three things JuliaXSPEC does, each in its simplest form:

1. `jlgauss`       — evaluate a model directly (a Gaussian line).
2. `jlgausscached` — the same model, cached on a parameter grid and interpolated.
3. `jlgconv`       — a convolution model (Gaussian blur in `g = E_obs / E_em`).

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
    convolve(spectrum, GaussianKernel(p.SigmaG))
end

# Registration must happen at load time, not at precompile time, so it lives
# in __init__.
function __init__()
    register!(jlgauss)
    register!(jlgausscached)
    register!(jlgconv)
end

end
