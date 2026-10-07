# Installation and build

## Requirements

- **Julia 1.13** (`juliaup add 1.13` if you use juliaup).
- **HEASOFT** with XSPEC, and the `HEADAS` environment initialised in the
  shell you build from (`source $HEADAS/headas-init.sh`). This provides
  `initpackage` and `hmake`, which XSPEC uses to build local models, and the
  C compiler configuration that goes with them.
- A clone of this repository:

```sh
git clone https://github.com/phajy/JuliaXSPEC.jl
cd JuliaXSPEC.jl
```

Nothing else is needed for the reference models. Later phases will add
optional dependencies (Gradus.jl, reflection tables) for the physics models.

## Step 1: compile the Julia models

```sh
./build-julia.sh
```

This instantiates the `models/` and `scripts/` environments, then runs
`scripts/compile.jl`, which uses
[PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/stable/)
to compile the model set in `models/` — together with JuliaXSPEC and the
Julia runtime — into

```text
build/lib/libjuliaxspec_models.dylib   (or .so on Linux)
build/include/julia_init.h
```

and finally writes, from the registered models,

```text
xspec/model.dat                XSPEC's description of each model and its parameters
xspec/juliaxspec_wrappers.c    one small C function per model
```

Compilation takes several minutes the first time (it builds a full Julia
system image); `build/` is about 300 MB and is not tracked by git. The two
generated files in `xspec/` *are* tracked, so you can read exactly what XSPEC
is given without building.

Repeat this step whenever you change Julia code in `src/` or `models/`.

## Step 2: build the XSPEC local-model package

```sh
./build-xspec.sh
```

With `HEADAS` set, this runs in `xspec/`:

1. `initpackage juliaxspec model.dat .` — XSPEC generates a `Makefile` and the
   C++ glue that registers the models;
2. `julia ../scripts/patch_xspec_makefile.jl` — adds the link flags for
   `libjuliaxspec_models` to that `Makefile` (XSPEC's generated `Makefile`
   does not know about it);
3. `hmake` — compiles the wrappers and links `xspec/libjuliaxspec.dylib`.

Repeat this step whenever `model.dat` changes (new models or parameters).
If only the Julia code behind an existing model changed, step 1 is enough:
XSPEC loads the Julia library at run time.

## Step 3: load in XSPEC

From the repository root:

```text
xspec
XSPEC12> lmod juliaxspec xspec
XSPEC12> model jlgauss
```

`lmod` takes the package name and the directory containing the library. The
first evaluation of the session starts the Julia runtime (a second or two);
after that, calls are as fast as the model itself.

To load the package automatically in every session, add the `lmod` line with
an absolute path to `~/.xspec/xspec.rc`.

## Configuration

JuliaXSPEC reads these environment variables (set them before starting
`xspec`):

- `JULIAXSPEC_VERBOSE=1` — print one line per model evaluation with the
  parameters, number of bins and timing.
- `JULIAXSPEC_CONVOLVE=fft` — make [`convolve`](@ref) use the FFT path unless a
  model asks for a method explicitly. The default is `direct`.
- `JULIAXSPEC_FFT_NBINS` — number of logarithmic bins in an FFT convolution.
- `JULIAXSPEC_TABLE_DIR` — directory to search first for table files such as
  `xillverD-5.fits`. Otherwise the repository root and the current directory
  are tried. The file is not part of the git repository.
- `JULIAXSPEC_CACHE_DIR` — where grid corners are written (default
  `~/.julia/juliaxspec`).
- `JULIAXSPEC_CACHE_LIMIT_GB` — shared RAM budget for cached grid corners
  (default 16; `0` for no limit).

## Troubleshooting

On macOS, FFTW (used by the fast convolution) pulls in Intel's oneTBB
threading library. Its malloc-replacement library crashes XSPEC, which
allocated its memory before our library was loaded. `scripts/compile.jl`
replaces that one library with an empty one after the build; the rest of
oneTBB is unused. Linux does not need this.

- **`lmod` fails to find the library.** Check that `xspec/libjuliaxspec.*`
  exists (step 2 succeeded) and that the path given to `lmod` is the `xspec/`
  directory.
- **The library loads but a model prints `JuliaXSPEC: error evaluating model …`.**
  The Julia exception and stack trace are printed to the terminal; the model
  returns zeros (additive) or leaves its input unchanged (convolution) so that
  the fit does not proceed silently with wrong values.
- **Changes to `models/` have no effect.** Re-run `./build-julia.sh`; the
  compiled library is a snapshot of the Julia code at build time.
- **`jltable` cannot find `xillverD-5.fits`.** Put the file in the repository
  root, or set `JULIAXSPEC_TABLE_DIR` to the directory that contains it, and
  start XSPEC again.
