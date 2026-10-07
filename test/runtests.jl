using Test
using JuliaXSPEC

# Energy grids used throughout: a log-spaced "response-like" grid and a coarser one.
const EDGES = collect(logrange(0.1, 50.0, 2001))
const COARSE = collect(logrange(0.1, 50.0, 201))

@testset "JuliaXSPEC" begin

@testset "BinnedSpectrum: units and rebinning" begin
    s = BinnedSpectrum(EDGES, gaussian_per_bin(EDGES, 6.4, 0.3))
    @test nbins(s) == 2000
    @test sum(s.per_bin) ≈ 1.0
    @test length(bin_widths(s)) == 2000
    @test bin_centres(s)[1] ≈ (EDGES[1] + EDGES[2]) / 2

    # density is per keV; going back and forth is exact
    @test density(s) ≈ s.per_bin ./ diff(EDGES)
    @test from_density(EDGES, density(s)).per_bin ≈ s.per_bin

    # a flat density has a per-bin content proportional to the bin width
    flat = from_density(EDGES, fill(2.0, 2000))
    @test flat.per_bin ≈ 2.0 .* diff(EDGES)

    @test_throws ArgumentError BinnedSpectrum([1.0, 2.0, 3.0], [1.0])         # length mismatch
    @test_throws ArgumentError BinnedSpectrum([1.0, 3.0, 2.0], [1.0, 1.0])    # not increasing

    @testset "rebin conserves photons" begin
        coarser = rebin(s, COARSE)
        @test sum(coarser.per_bin) ≈ sum(s.per_bin)
        finer = rebin(coarser, EDGES)
        @test sum(finer.per_bin) ≈ sum(s.per_bin)
        # onto bins that are shifted relative to the originals and cover them
        shifted = rebin(s, collect(range(0.05, 60.0, length = 777)))
        @test sum(shifted.per_bin) ≈ sum(s.per_bin)
        # a destination grid that covers only part of the source keeps only that part
        partial = rebin(s, collect(range(6.4, 50.0, length = 100)))
        @test isapprox(sum(partial.per_bin), 0.5; atol = 1e-3)                 # half the Gaussian lies above its centre
        # rebinning onto the same grid is the identity
        @test rebin(s, EDGES).per_bin ≈ s.per_bin
    end
end

@testset "Parameter and model.dat lines" begin
    p = Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6, delta = 0.05)
    @test model_dat_line(p) == "LineE keV 6.4 0.0 0.0 1.0e6 1.0e6 0.05"
    q = Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 20.0, soft_max = 10.0, delta = 0.05)
    @test model_dat_line(q) == "Sigma keV 0.1 0.0 0.0 10.0 20.0 0.05"
    frozen = Parameter("Incl", 30.0; unit = "deg", min = 0.0, max = 90.0, delta = 1.0, frozen = true)
    @test endswith(model_dat_line(frozen), " -1.0")
    unitless = Parameter("spin", 0.9; min = 0.0, max = 0.998)
    @test startswith(model_dat_line(unitless), "spin \" \" 0.9")
    @test model_dat_line(Parameter("redshift", 0.0; min = 0.0, max = 1.0, kind = :scale)) == "*redshift \" \" 0.0"
    @test model_dat_line(Parameter("switch", 2.0; min = 0.0, max = 3.0, kind = :switch)) == "\$switch 2"

    @test_throws ArgumentError Parameter("bad name", 1.0; min = 0.0, max = 2.0)
    @test_throws ArgumentError Parameter("x", 5.0; min = 0.0, max = 2.0)       # default outside limits
    @test_throws ArgumentError Parameter("x", 1.0; min = 0.0, max = 2.0, kind = :weird)
end

@testset "Models" begin
    gauss = AdditiveModel("testgauss", [
        Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6),
        Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 10.0),
    ]) do edges, p
        gaussian_per_bin(edges, p.LineE, p.Sigma)
    end
    @test xspec_type(gauss) == "add"
    @test parameter_names(gauss) == [:LineE, :Sigma]
    @test parameter_defaults(gauss) == [6.4, 0.1]
    @test named_parameters(gauss, [6.0, 0.2]) == (LineE = 6.0, Sigma = 0.2)
    @test evaluate(gauss, EDGES, [6.4, 0.3]) ≈ gaussian_per_bin(EDGES, 6.4, 0.3)
    @test_throws ArgumentError evaluate(gauss, EDGES, [6.4])

    blur = ConvolutionModel("testblur", [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5)]) do s, p
        convolve(s, GaussianKernel(p.SigmaG))
    end
    @test xspec_type(blur) == "con"
    input = BinnedSpectrum(EDGES, gaussian_per_bin(EDGES, 6.4, 0.01))
    output = evaluate(blur, input, [0.05])
    @test output isa BinnedSpectrum
    @test output.edges == input.edges

    # a convolution model that forgets the contract is caught
    broken = ConvolutionModel("broken", Parameter[]) do s, p
        s.per_bin
    end
    @test_throws ErrorException evaluate(broken, input, Float64[])

    @test_throws ArgumentError AdditiveModel((e, p) -> e, "1bad", Parameter[])
    @test_throws ArgumentError AdditiveModel((e, p) -> e, "dup", [
        Parameter("a", 1.0; min = 0.0, max = 2.0), Parameter("a", 1.0; min = 0.0, max = 2.0)])
end

@testset "Registry and generated XSPEC files" begin
    empty!(JuliaXSPEC.REGISTRY)
    a = AdditiveModel((e, p) -> gaussian_per_bin(e, p.LineE, 0.1), "test_line",
                      [Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6, delta = 0.05)];
                      description = "a test line")
    c = ConvolutionModel((s, p) -> s, "testconv", [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5)])
    register!(a)
    register!(c)
    @test [m.name for m in registered_models()] == ["test_line", "testconv"]
    @test lookup_model("testconv") === c
    @test_throws ErrorException lookup_model("nope")
    register!(AdditiveModel((e, p) -> zeros(length(e) - 1), "test_line", Parameter[]))   # replaces, keeps order
    @test [m.name for m in registered_models()] == ["test_line", "testconv"]
    @test isempty(lookup_model("test_line").parameters)
    register!(a)

    @test c_function_name(a) == "testline"
    text = model_dat_text()
    @test text == """
        test_line 1 0. 1.e20 c_testline add 0
        LineE keV 6.4 0.0 0.0 1.0e6 1.0e6 0.05

        testconv 1 0. 1.e20 c_testconv con 0
        SigmaG " " 0.05 0.001 0.001 0.5 0.5 0.005
        """

    wrappers = c_wrappers_text()
    @test occursin("void testline(const double *energy", wrappers)
    @test occursin("juliaxspec_evaluate(\"test_line\", energy", wrappers)
    @test occursin("void testconv(", wrappers)
    @test occursin("a test line", wrappers)

    dir = mktempdir()
    model_dat, cfile = write_xspec_package(dir; julia_init_header = "../build/include/julia_init.h")
    @test read(model_dat, String) == text
    @test read(cfile, String) == wrappers

    # names that collide once underscores are removed are refused
    register!(AdditiveModel((e, p) -> zeros(length(e) - 1), "testline", Parameter[]))
    @test_throws ErrorException write_xspec_package(dir)
    empty!(JuliaXSPEC.REGISTRY)
end

@testset "GridAxis" begin
    @test GridAxis(range(1.0, 2.0, length = 5)).scale == :linear
    @test GridAxis(logrange(1.0, 100.0, 5)).scale == :log
    @test GridAxis([1.0, 2.0, 3.0]).scale == :linear
    @test GridAxis([1.0, 10.0, 100.0]).scale == :log
    @test GridAxis([1.0, 2.0, 10.0]).scale == :linear                 # irregular: plain linear
    @test GridAxis([1.0, 2.0, 10.0]; scale = :log).scale == :log      # or as requested
    @test_throws ArgumentError GridAxis([3.0, 2.0, 1.0])
    @test_throws ArgumentError GridAxis([-1.0, 1.0, 10.0]; scale = :log)

    lin = GridAxis([0.0, 1.0, 2.0])
    @test JuliaXSPEC.bracket(lin, 0.25) == (1, 2, 0.25)
    @test JuliaXSPEC.bracket(lin, 1.0) == (1, 2, 1.0) || JuliaXSPEC.bracket(lin, 1.0) == (2, 3, 0.0)
    @test JuliaXSPEC.bracket(lin, -5.0) == (1, 1, 0.0)                # clamped below
    @test JuliaXSPEC.bracket(lin, 7.0)[3] ≈ 1.0                       # clamped above
    log10axis = GridAxis([1.0, 10.0, 100.0])
    lo, hi, t = JuliaXSPEC.bracket(log10axis, sqrt(10))               # halfway in log
    @test (lo, hi) == (1, 2) && t ≈ 0.5
end

@testset "GridInterpolator" begin
    # Exact for functions linear in each axis coordinate (including the cross term).
    calls = Ref(0)
    g = GridInterpolator((x = range(0.0, 1.0, length = 3), y = logrange(1.0, 100.0, 3))) do x, y
        calls[] += 1
        [1.0 + 2x + 3log(y) + x * log(y), 10x]
    end
    @test grid_size(g) == 9
    v = g(0.3, 7.0)
    @test v ≈ [1.0 + 0.6 + 3log(7.0) + 0.3log(7.0), 3.0]
    @test v ≈ evaluate_exact(g, 0.3, 7.0)
    @test cache_stats(g) == (hits = 0, misses = 4, stored = 4, total = 9)
    g(0.3, 7.0)
    @test cache_stats(g).hits == 4 && cache_stats(g).misses == 4
    g(0.4, 5.0)                                                       # same cell: all hits
    @test cache_stats(g).misses == 4
    g(0.9, 50.0)                                                      # neighbouring cell: shares one corner
    @test cache_stats(g).misses == 7
    @test calls[] == 7 + 1                                            # 7 corners plus the evaluate_exact call

    # Points on a grid line use fewer corners; values outside are clamped.
    g(0.5, 1.0)
    @test g(0.5, 0.01) == g(0.5, 1.0)
    @test g(5.0, 1.0) == g(1.0, 1.0)

    empty_cache!(g)
    @test cache_stats(g) == (hits = 0, misses = 0, stored = 0, total = 9)

    # Single axis, vector input, and the error for the wrong number of arguments.
    h = GridInterpolator((s = [0.0, 0.5, 1.0],)) do s
        [s^2]
    end
    @test h(0.25) ≈ [0.125]                                           # linear between 0 and 0.25
    @test_throws ArgumentError h(0.25, 0.5)
    @test_throws ArgumentError GridInterpolator(s -> [s], NamedTuple())
end

@testset "Analytic spectra" begin
    gp = gaussian_per_bin(EDGES, 6.4, 0.3)
    @test sum(gp) ≈ 1.0
    @test all(>=(0), gp)
    # Agrees with the density at bin centres for narrow bins
    centres = bin_centres(EDGES)
    dens = @. exp(-0.5 * ((centres - 6.4) / 0.3)^2) / (0.3 * sqrt(2π))
    @test isapprox(gp ./ diff(EDGES), dens; rtol = 1e-3)
    delta = gaussian_per_bin(EDGES, 6.4, 0.0)
    @test sum(delta) == 1.0 && count(>(0), delta) == 1
    @test EDGES[findfirst(>(0), delta)] <= 6.4 < EDGES[findfirst(>(0), delta) + 1]
    @test_throws ArgumentError gaussian_per_bin(EDGES, 6.4, -1.0)

    pl = powerlaw_per_bin(EDGES, 2.0)
    @test pl ≈ (1 ./ EDGES[1:end-1] .- 1 ./ EDGES[2:end])
    @test powerlaw_per_bin([1.0, 2.0], 1.0) ≈ [log(2.0)]
end

@testset "RedshiftKernel" begin
    g = collect(range(0.5, 1.5, length = 501))
    k = RedshiftKernel(g, @. exp(-0.5 * ((g - 1) / 0.1)^2))
    @test JuliaXSPEC.trapezoid_area(k.g, k.L) ≈ 1.0
    @test k(1.0) ≈ maximum(k.L)
    @test k(0.4) == 0.0 && k(1.6) == 0.0
    @test cumulative(k, 0.4) == 0.0 && cumulative(k, 1.6) == 1.0
    @test cumulative(k, 1.0) ≈ 0.5 atol = 1e-6
    @test issorted([cumulative(k, x) for x in range(0.5, 1.5, length = 101)])
    # exact at and between grid points
    @test cumulative(k, g[100]) ≈ k.Λ[100]
    mid = (g[100] + g[101]) / 2
    @test cumulative(k, mid) ≈ k.Λ[100] + 0.5 * (k.L[100] + k(mid)) * (mid - g[100])

    gk = GaussianKernel(0.05)
    @test gk.g[1] ≈ 0.7 && gk.g[end] ≈ 1.3
    @test JuliaXSPEC.trapezoid_area(gk.g, gk.L) ≈ 1.0
    @test JuliaXSPEC.trapezoid_area(gk.g, gk.g .* gk.L) ≈ 1.0 atol = 1e-6        # mean g = 1
    @test_throws ArgumentError GaussianKernel(0.0)
    @test_throws ArgumentError RedshiftKernel([1.0, 0.5], [1.0, 1.0])
    @test_throws ArgumentError RedshiftKernel([0.5, 1.0], [0.0, 0.0])
end

@testset "Convolution" begin
    kernel = GaussianKernel(0.05)

    @testset "photons are conserved" begin
        M = convolution_matrix(EDGES, EDGES, kernel)
        colsums = vec(sum(M; dims = 1))
        interior = (EDGES[1:end-1] .> 0.1 / 0.7) .& (EDGES[2:end] .< 50.0 * 0.7)
        @test all(isapprox.(colsums[interior], 1.0; atol = 1e-6))
        @test all(colsums .<= 1.0 + 1e-9)
    end

    @testset "a narrow kernel is the identity" begin
        s = BinnedSpectrum(EDGES, gaussian_per_bin(EDGES, 6.4, 0.3))
        tiny = GaussianKernel(1e-5)
        @test convolve(s, tiny; n_sub = 8).per_bin ≈ s.per_bin rtol = 1e-3
    end

    @testset "a kernel centred at g₀ shifts a line to g₀ E₀" begin
        g = collect(range(0.6, 1.0, length = 2001))
        shift = RedshiftKernel(g, @. exp(-0.5 * ((g - 0.8) / 0.002)^2))
        s = BinnedSpectrum(EDGES, gaussian_per_bin(EDGES, 6.4, 0.05))
        out = convolve(s, shift)
        @test sum(out.per_bin) ≈ 1.0 atol = 1e-6
        centroid = sum(bin_centres(out) .* out.per_bin)
        @test centroid ≈ 0.8 * 6.4 rtol = 1e-3
    end

    @testset "Gaussian in g applied to a narrow line is a Gaussian in energy" begin
        E0, σ_in, σ_g = 6.4, 0.02, 0.05
        s = BinnedSpectrum(EDGES, gaussian_per_bin(EDGES, E0, σ_in))
        out = convolve(s, kernel)
        expected = gaussian_per_bin(EDGES, E0, sqrt((σ_g * E0)^2 + σ_in^2))
        @test out.per_bin ≈ expected rtol = 2e-3
        @test sum(out.per_bin) ≈ 1.0 atol = 1e-6
    end

    @testset "a power law stays a power law, rescaled by ∫ L(g) g^(Γ-1) dg" begin
        Γ = 2.0
        s = BinnedSpectrum(EDGES, powerlaw_per_bin(EDGES, Γ))
        out = convolve(s, kernel)
        factor = JuliaXSPEC.trapezoid_area(kernel.g, kernel.L .* kernel.g .^ (Γ - 1))
        interior = (EDGES[1:end-1] .> 0.1 / 0.7) .& (EDGES[2:end] .< 50.0 * 0.7)
        @test out.per_bin[interior] ≈ factor .* s.per_bin[interior] rtol = 1e-3
    end

    @testset "output grid may differ from the input grid" begin
        s = BinnedSpectrum(EDGES, gaussian_per_bin(EDGES, 6.4, 0.3))
        out = convolve(s, kernel; out_edges = COARSE)
        @test out.edges == COARSE
        @test sum(out.per_bin) ≈ 1.0 atol = 1e-6
        @test out.per_bin ≈ rebin(convolve(s, kernel), COARSE).per_bin rtol = 1e-6
    end

    @test_throws ArgumentError convolution_matrix(EDGES, EDGES, kernel; n_sub = 0)
end

@testset "C entry point" begin
    empty!(JuliaXSPEC.REGISTRY)
    register!(AdditiveModel((e, p) -> gaussian_per_bin(e, p.LineE, p.Sigma), "entryline", [
        Parameter("LineE", 6.4; unit = "keV", min = 0.0, max = 1e6),
        Parameter("Sigma", 0.1; unit = "keV", min = 0.0, max = 10.0)]))
    register!(ConvolutionModel((s, p) -> convolve(s, GaussianKernel(p.SigmaG)), "entryblur",
        [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5)]))
    register!(AdditiveModel((e, p) -> error("boom"), "entryfail", Parameter[]))
    register!(MultiplicativeModel((e, p) -> error("boom"), "entrymulfail", Parameter[]))

    function call(name, values, flux)
        edges = copy(EDGES)
        flux_error = zeros(length(flux))
        GC.@preserve edges values flux flux_error begin
            JuliaXSPEC.juliaxspec_evaluate(
                Base.unsafe_convert(Ptr{Cchar}, name), pointer(edges), Cint(length(edges) - 1),
                pointer(values), Cint(1), pointer(flux), pointer(flux_error), Base.unsafe_convert(Ptr{Cchar}, ""))
        end
    end

    flux = zeros(2000)
    @test call("entryline", [6.4, 0.3], flux) == 0
    @test flux ≈ gaussian_per_bin(EDGES, 6.4, 0.3)

    # convolution: the buffer holds the input and is overwritten with the result
    flux = gaussian_per_bin(EDGES, 6.4, 0.02)
    @test call("entryblur", [0.05], flux) == 0
    @test flux ≈ gaussian_per_bin(EDGES, 6.4, sqrt((0.05 * 6.4)^2 + 0.02^2)) rtol = 2e-3

    # failures are reported, not thrown, and leave a safe result behind
    flux = fill(7.0, 2000)
    stderr_text = mktemp() do path, io
        redirect_stderr(io) do
            @test call("entryfail", Float64[], flux) == 1
        end
        flush(io)
        read(path, String)
    end
    @test occursin("boom", stderr_text)
    @test all(==(0.0), flux)
    flux = fill(7.0, 2000)
    redirect_stderr(devnull) do
        @test call("entrymulfail", Float64[], flux) == 1
    end
    @test all(==(1.0), flux)
    redirect_stderr(devnull) do
        @test call("nonexistent", Float64[], flux) == 1
    end

    # verbose logging goes to stdout
    withenv("JULIAXSPEC_VERBOSE" => "1") do
        out = mktemp() do path, io
            redirect_stdout(io) do
                call("entryline", [6.4, 0.3], zeros(2000))
            end
            flush(io)
            read(path, String)
        end
        @test occursin("entryline(LineE=6.4, Sigma=0.3)", out)
    end
    empty!(JuliaXSPEC.REGISTRY)
end

@testset "FFT convolution" begin
    @test convolution_method() === :direct
    withenv("JULIAXSPEC_CONVOLVE" => "fft") do
        @test convolution_method() === :fft
    end
    withenv("JULIAXSPEC_CONVOLVE" => "nope") do
        @test_throws ArgumentError convolution_method()
    end

    edges = collect(logrange(1.0, 20.0, 257))
    spectrum = BinnedSpectrum(edges, gaussian_per_bin(edges, 6.4, 0.2))
    kernel = GaussianKernel(0.05)
    direct = convolve(spectrum, kernel; method = :direct)
    fast = convolve(spectrum, kernel; method = :fft)
    @test fast.per_bin ≈ convolve_fft(spectrum, kernel).per_bin
    rel = sum(abs.(direct.per_bin .- fast.per_bin)) / sum(direct.per_bin)
    @test rel < 1e-3
    @test sum(fast.per_bin) ≈ 1 atol = 1e-4
    @test all(fast.per_bin .>= 0)
end

@testset "Grid cache on disk and under a memory budget" begin
    dir = mktempdir()
    calls = Ref(0)
    withenv("JULIAXSPEC_CACHE_DIR" => dir) do
        g = GridInterpolator((x = [0.0, 1.0],); cache = "phase2-disk") do x
            calls[] += 1
            [x, 2x]
        end
        @test g(0.0) == [0.0, 0.0]
        @test g(1.0) == [1.0, 2.0]
        @test calls[] == 2
        empty_cache!(g)
        @test g(0.0) == [0.0, 0.0]
        @test calls[] == 2
        @test disk_loads(g) == 1

        g2 = GridInterpolator((x = [0.0, 1.0],); cache = "phase2-disk") do x
            calls[] += 1
            [x, 2x]
        end
        @test g2(1.0) == [1.0, 2.0]
        @test calls[] == 2
        @test disk_loads(g2) == 1

        # A different grid under the same name must not reuse the old files.
        g3 = GridInterpolator((x = [0.0, 2.0],); cache = "phase2-disk") do x
            calls[] += 1
            [10x]
        end
        @test g3(0.0) == [0.0]
        @test calls[] == 3
    end
    @test_throws ArgumentError GridInterpolator((x = [0.0, 1.0],); cache = "has a space") do x
        [x]
    end

    for w in copy(JuliaXSPEC.LIVE_CACHES)
        c = w.value
        c === nothing || empty_cache!(c)
    end
    @test cache_memory_used_bytes() == 0
    calls2 = Ref(0)
    limit_gb = 32 / 2^30          # two corners of two Float64s
    withenv("JULIAXSPEC_CACHE_LIMIT_GB" => string(limit_gb)) do
        h = GridInterpolator((x = [0.0, 1.0, 2.0, 3.0],)) do x
            calls2[] += 1
            [x, x]
        end
        h(0.0)
        h(1.0)
        h(2.0)
        @test cache_stats(h).stored == 2
        @test calls2[] == 3
        h(0.0)                    # the oldest corner was dropped
        @test calls2[] == 4
    end
end

@testset "OGIP table models" begin
    edges = collect(range(1.0, 5.0, length = 5))
    parameters = [
        Parameter("p", 1.5; min = 1.0, max = 2.0, delta = 0.1),
        Parameter("q", 3.0; min = 1.0, max = 10.0, delta = 0.1),
    ]
    # spectra[energy, q, p], values p + 10q at the integer corner labels
    spectra = zeros(Float32, 4, 2, 2)
    for ip in 1:2, iq in 1:2
        spectra[:, iq, ip] .= ip + 10 * iq
    end
    path = tempname() * ".fits"
    write_ogip_table(path, spectra; name = "toy", parameters, axes = [[1.0, 2.0], [1.0, 10.0]], methods = [0, 1], edges)
    table = OGIPTable(path)
    @test table.name == "toy"
    @test table.interpolator.axes[1].scale == :linear
    @test table.interpolator.axes[2].scale == :log
    @test table(1.0, 1.0) ≈ fill(11.0, 4)
    @test table(2.0, 10.0) ≈ fill(22.0, 4)
    # Halfway in p, and halfway in log q, is the mean of the four corners.
    @test table(1.5, sqrt(10)) ≈ fill(16.5, 4)
    @test table(-1.0, 100.0) ≈ table(1.0, 10.0)          # clamped to the grid

    model = Blurred("toyblur", () -> table,
        [Parameter("SigmaG", 0.05; min = 1e-3, max = 0.5, delta = 0.01)],
        table.parameters; method = :direct) do kernel
        GaussianKernel(kernel.SigmaG)
    end
    got = evaluate(model, edges, [0.05, 1.0, 1.0])
    source = BinnedSpectrum(table.edges, table(1.0, 1.0))
    @test got ≈ convolve(source, GaussianKernel(0.05); out_edges = edges, method = :direct).per_bin
end

end
