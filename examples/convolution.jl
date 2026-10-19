# Blur a spectrum with a kernel in g = E_obs / E_em.
#
# The direct matrix and the FFT are the same operator. From the repository root:
#   julia --project=. examples/convolution.jl

using JuliaXSPEC

edges = collect(logrange(1.0, 20.0, 513))
line = BinnedSpectrum(edges, gaussian_per_bin(edges, 6.4, 0.05))
kernel = GaussianKernel(0.05)   # unit area in g; a narrow line at E0 becomes width 0.05*E0

direct = convolve(line, kernel; method = :direct)
fast = convolve(line, kernel; method = :fft)

println("direct total: ", sum(direct.per_bin))
println("fft total:    ", sum(fast.per_bin))
rel = sum(abs.(direct.per_bin .- fast.per_bin)) / sum(direct.per_bin)
println("FFT vs direct, relative L1: ", rel)

abs(sum(direct.per_bin) - 1) < 1e-6 || error("direct convolution did not conserve photons")
rel < 1e-3 || error("FFT and direct convolution disagree: relative L1 = $rel")
println("ok")
