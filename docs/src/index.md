```@meta
CurrentModule = SubpixelRegistration
```

# SubpixelRegistration.jl

[![Code](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/JuliaHCI/SubpixelRegistration.jl)
[![Build Status](https://github.com/JuliaHCI/SubpixelRegistration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaHCI/SubpixelRegistration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SubpixelRegistration.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaHCI/SubpixelRegistration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaHCI/SubpixelRegistration.jl)
[![License](https://img.shields.io/github/license/JuliaHCI/SubpixelRegistration.jl?color=yellow)](https://github.com/JuliaHCI/SubpixelRegistration.jl/blob/main/LICENSE)

Image registration with subpixel precision using an upsampled discrete Fourier transform cross-correlation. This uses an efficient matrix-multiplication algorithm for upsampling the cross-correlation following Guizar-Sicairos, Thurman, and Fienup (2008).[^1]

[^1]: Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, ["Efficient subpixel image registration algorithms,"](http://www.opticsinfobase.org/ol/fulltext.cfm?uri=ol-33-2-156&id=148843) Opt. Lett. 33, 156-158 (2008)

## Installation

```julia
julia>] add SubpixelRegistration
```

## Usage

```@example test
using ColorTypes
using ImageIO
using ImageShow
using SubpixelRegistration
using TestImages

image = testimage("cameraman")
```

```@example test
shift = (22.4, -13.32)
source = Float64.(image)
shifted = fourier_shift(source, shift)
Gray.(shifted)
```

```@example test
phase_offset(source, shifted; upsample_factor=1)
```

```@example test
phase_offset(source, shifted; upsample_factor=100)
```

```@example test
registered = register(source, shifted; upsample_factor=100)
Gray.(registered)
```

## Benchmarks

This code has been benchmarked against the [scikit-image](https://github.com/scikit-image/scikit-image) implementation. This benchmark is a measure of the time it takes to measure the offset between two images with various sizes and with various upsample factors. The number of pixels scales with the square of the size, describing the non-linear power law.

**System Information**

```julia
julia> versioninfo()
Julia Version 1.7.0-rc3
Commit 3348de4ea6* (2021-11-15 08:22 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin20.5.0)
  CPU: Intel(R) Core(TM) i5-8259U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, skylake)
Environment:
  JULIA_NUM_THREADS = 1
```

**Benchmark Results**

```@example
using CSV, DataFrames, StatsPlots, SubpixelRegistration # hide
benchdir = joinpath(dirname(pathof(SubpixelRegistration)), "..", "bench") # hide
results = DataFrame(CSV.File(joinpath(benchdir, "benchmark_results.csv"))) # hide
groups = groupby(results, :upsample_factor) # hide
plot(xlabel="image size", ylabel="time (s)") # hide
shapes = [:o :dtriangle :diamond] # hide
for (g, shape) in zip(groups, shapes) # hide
    @df g plot!(:size, [:time_py :time_jl]; c=[1 2], shape, label="") # hide
end # hide
plot!(yscale=:log10) # hide
# create faux-legends # hide
bbox_ = bbox(0, 0, 1, 1, :bottom, :left) # hide
plot!([1 2]; c=[1 2], label=["scikit-image" "SubpixelRegistration.jl"], inset=(1, bbox_), # hide
    bg=:transparent, border=:none, axes=false, sp=2, leg=:topleft, bgcolorlegend=:white) # hide
ups = hcat((string(k.upsample_factor) for k in keys(groups))...) # hide
scatter!([0 0 0]; shape=shapes, c=:black, alpha=0.4, label=ups, inset=(1, bbox_), # hide
    bg=:transparent, border=:none, axes=false, sp=3, ylim=(1, 2), # hide
    legtitle="upsample\nfactor", leg=:bottomright, legendtitlefontsize=9, bgcolorlegend=:white) # hide

```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaHCI/SubpixelRegistration.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaHCI/SubpixelRegistration.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaHCI/SubpixelRegistration.jl/issues).
