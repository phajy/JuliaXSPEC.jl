# Replace oneTBB's libtbbmalloc_proxy in a PackageCompiler bundle with an empty
# library (macOS only).
#
# FFTW.jl depends on MKL_jll, which depends on oneTBB_jll, so a compiled
# library initialises oneTBB even though this package never calls it. On macOS
# that initialisation loads libtbbmalloc_proxy, whose only job is to replace
# the process malloc. XSPEC allocated memory before our library was loaded, and
# freeing it through the replacement allocator segfaults. Nothing here uses
# TBB (and MKL does not run on Apple Silicon), so the proxy is replaced with an
# empty library of the same name: oneTBB still opens it, and malloc is left
# alone. Linux does not interpose malloc this way, so there is nothing to do.

function stub_tbbmalloc_proxy!(build_dir::AbstractString = "build")
    if !Sys.isapple()
        println("stub_tbbmalloc_proxy: not macOS; nothing to do")
        return 0
    end
    artifacts = joinpath(build_dir, "share", "julia", "artifacts")
    if !isdir(artifacts)
        println("stub_tbbmalloc_proxy: no artifacts directory at $artifacts")
        return 0
    end

    stub_src = tempname() * ".c"
    write(stub_src, "/* intentionally empty: stub for libtbbmalloc_proxy */\n")

    n = 0
    for (root, _dirs, files) in walkdir(artifacts)
        for file in files
            startswith(file, "libtbbmalloc_proxy") && endswith(file, ".dylib") || continue
            path = joinpath(root, file)
            islink(path) && continue
            stub = path * ".stub.dylib"
            run(`cc -dynamiclib -o $stub $stub_src -install_name @rpath/libtbbmalloc_proxy.2.dylib`)
            rm(path)               # the bundled file is read-only; replace the directory entry
            mv(stub, path)
            chmod(path, 0o555)
            println("stub_tbbmalloc_proxy: replaced $path")
            n += 1
        end
    end
    rm(stub_src; force = true)
    if n == 0
        println("stub_tbbmalloc_proxy: no libtbbmalloc_proxy library found")
    end
    return n
end

if abspath(PROGRAM_FILE) == @__FILE__
    stub_tbbmalloc_proxy!(isempty(ARGS) ? "build" : ARGS[1])
end
