# Installation and build

## Requirements

- **Julia 1.13** (`juliaup add 1.13` if you use juliaup).
- **HEASOFT** with XSPEC, and the `HEADAS` environment initialised in the
  shell you build from (`source $HEADAS/headas-init.sh`). This provides
  `initpackage` and `hmake`, which XSPEC uses to build local models, and the
  C compiler configuration that goes with them.
There are two ways to get the package, and they are not the same thing.

**The Julia package.** JuliaXSPEC is in the
[AstroRegistry](https://github.com/astro-group-bristol/AstroRegistry), not the
General registry. This installs the functions (`AdditiveModel`, `convolve`,
`GridInterpolator`, and the rest) and nothing else:

```julia
using Pkg
Pkg.Registry.add(url = "https://github.com/astro-group-bristol/AstroRegistry")
Pkg.add("JuliaXSPEC")
```

**The XSPEC library.** `Pkg.add` does not build the shared library XSPEC
loads. That build uses the reference models and the compile scripts in this
repository, so it starts from a clone:

```sh
git clone https://github.com/phajy/JuliaXSPEC.jl
cd JuliaXSPEC.jl
```

Nothing else is needed for the reference models. Gradus and kerrz will be
separate packages, not dependencies of JuliaXSPEC.

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

On Linux, `lmod` can fail with `GLIBCXX_3.4.30 not found`. Julia's C++ library
is newer than `/lib64/libstdc++.so.6`, and XSPEC has already loaded the system
copy. From the repository directory, before starting XSPEC, put the bundled
library first on the path. For the reference models that directory is
`build/lib/julia`. For the Gradus build it is `build/gradus/lib/julia`:

```sh
export LD_LIBRARY_PATH="$PWD/build/gradus/lib/julia:${LD_LIBRARY_PATH}"
```

Setting this inside the XSPEC prompt does not help.

## The thin-disc lamppost

`jllamp` smears whatever it multiplies by the line profile of a lamppost above
a thin disc. From a clone of this repository:

```sh
./build-julia.sh gradus
./build-xspec.sh gradus
```

The first command downloads Gradus and takes several minutes. Load the library
it prints, instead of the reference-model library, and in a new XSPEC session:

```text
model jllamp * gaussian
```

`jllampxill` blurs `xillverD-5.fits` with that same kernel, on the table's own
energy grid. The file goes in the repository folder, or in the folder named by
`JULIAXSPEC_TABLE_DIR`. In XSPEC, `model jllampxill`. Its first three parameters
are the lamppost; the next five are the table, including the table inclination,
which is separate from the lamppost inclination.

The lamppost parameters are spin (0.2 to 0.998), inclination in degrees (10 to
60), and lamp height in gravitational radii (3 to 20). `jlring` is a
co-rotating ring above the same disc, and `jlringxill` blurs the xillver table
with it. The extra parameter is the ring radius (2 to 50), and the ring height
starts at 2.5. In XSPEC, `model jlring * gaussian` or `model jlringxill`.

Each profile is calculated at a handful of parameter values and blended in
between. The first time a combination is needed it can take a few minutes, and
a ring takes longer than a lamppost. It is then remembered. A session can load
only one of these libraries.

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
