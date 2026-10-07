# The same blurring operator as convolution.jl, evaluated as a convolution in
# ln E with an FFT. Step 6 of the Convolution documentation page derives it.
# The direct matrix remains the definition; this is checked against it.

using FFTW

"""
    convolve_fft(spectrum, kernel; out_edges = spectrum.edges, n_bins = nothing) -> BinnedSpectrum

[`convolve`](@ref) via the convolution theorem.

With ``u = \\ln E``, the blurring integral is an ordinary convolution

```math
n_{\\mathrm{obs}}(e^u) = \\int L(e^{u-v})\\, n_{\\mathrm{em}}(e^v)\\, dv ,
```

so on a uniform grid in ``\\ln E`` it is a product of Fourier transforms.
`n_bins` is the number of logarithmic bins. The default is about 0.5% wide,
and at least as many bins as the input spectrum (or `JULIAXSPEC_FFT_NBINS`).
The input is rebinned onto that grid, the
convolution is done there, and the result is rebinned onto `out_edges`.

Tiny negative values from the transform are clipped to zero; on a fine grid
they are numerical noise, far below a bin.
"""
function convolve_fft(
    spectrum::BinnedSpectrum,
    kernel::RedshiftKernel;
    out_edges::AbstractVector{<:Real} = spectrum.edges,
    n_bins::Union{Nothing, Integer} = nothing,
)
    E_min = min(Float64(first(spectrum.edges)), Float64(first(out_edges)))
    E_max = max(Float64(last(spectrum.edges)), Float64(last(out_edges)))
    E_min > 0 || throw(ArgumentError("FFT convolution needs positive energies"))
    n = n_bins === nothing ? default_fft_nbins(E_min, E_max, nbins(spectrum)) : Int(n_bins)
    n >= 16 || throw(ArgumentError("n_bins must be at least 16"))

    log_edges = exp.(range(log(E_min), log(E_max), length = n + 1))
    dτ = log(log_edges[2] / log_edges[1])
    rebinned = rebin(spectrum, log_edges)
    widths = diff(log_edges)
    density_in = rebinned.per_bin ./ widths          # photons / cm² / s / keV

    # Zero-pad by the kernel's width in τ so the circular convolution does not
    # wrap a redshift off one end of the grid onto the other end.
    n_kernel = ceil(Int, (log(kernel.g[end]) - log(kernel.g[1])) / dτ) + 3
    n_fft = nextprod((2, 3, 5), n + n_kernel)

    padded = zeros(n_fft)
    padded[1:n] .= density_in
    kernel_bins = fft_kernel(kernel, dτ, n_fft)

    transformed = rfft(padded)
    transformed .*= rfft(kernel_bins)
    density_out = irfft(transformed, n_fft)

    per_bin = max.(density_out[1:n], 0.0) .* widths
    return rebin(BinnedSpectrum(log_edges, per_bin), out_edges)
end

"Bin count for the logarithmic working grid: at least the input's, and about 0.5% wide."
function default_fft_nbins(E_min::Float64, E_max::Float64, n_input::Int)
    raw = strip(get(ENV, "JULIAXSPEC_FFT_NBINS", ""))
    if !isempty(raw)
        n = tryparse(Int, raw)
        n !== nothing && n >= 16 && return n
    end
    from_spacing = ceil(Int, log(E_max / E_min) / log(1.005))
    return max(1024, n_input, from_spacing)
end

"""
Kernel samples for a circular convolution of spacing `dτ`.

Index 1 is ``τ = 0``. Positive ``τ`` (blueshift, ``g > 1``) follows it;
negative ``τ`` is stored at the end of the array, which is where an FFT puts
negative lags. Each sample is ``L(e^τ)\\, dτ``.
"""
function fft_kernel(kernel::RedshiftKernel, dτ::Float64, n_fft::Int)
    K = zeros(n_fft)
    K[1] = kernel(1.0) * dτ
    half = n_fft ÷ 2
    for k in 1:half
        K[k + 1] = kernel(exp(k * dτ)) * dτ
        if k < n_fft - k
            K[n_fft - k + 1] = kernel(exp(-k * dτ)) * dτ
        end
    end
    return K
end
