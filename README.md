# SubpixelRegistration.jl

[![Build Status](https://github.com/JuliaHCI/SubpixelRegistration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaHCI/SubpixelRegistration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SubpixelRegistration.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaHCI/SubpixelRegistration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaHCI/SubpixelRegistration.jl)
[![License](https://img.shields.io/github/license/JuliaHCI/SubpixelRegistration.jl?color=yellow)](LICENSE)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaHCI.github.io/SubpixelRegistration.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaHCI.github.io/SubpixelRegistration.jl/dev)

Image registration with subpixel precision using an upsampled discrete Fourier transform cross-correlation. This uses an efficient matrix-multiplication algorithm for upsampling the cross-correlation following Guizar-Sicairos, Thurman, and Fienup (2008).[^1]

[^1]: Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, ["Efficient subpixel image registration algorithms,"](http://www.opticsinfobase.org/ol/fulltext.cfm?uri=ol-33-2-156&id=148843) Opt. Lett. 33, 156-158 (2008)

## Installation

```julia
julia>] add SubpixelRegistration
```

## Usage

First, import the package

```julia
julia> using SubpixelRegistration
```

now, let's load two images we wish to register

```julia
julia> source, target = # ...
```

we can use `phase_offset` directly to perform the DFT-upsampling cross-correlation analysis

```julia
julia> phase_offset(source, target; upsample_factor=1)
(shift = (-22.0, 13.0), error = 0.9999999999214572, phasediff = 0.0)

julia> phase_offset(source, target; upsample_factor=4)
(shift = (-22.5, 13.25), error = -11.584699520819848, phasediff = 0.0)
```

the `register` and `coregister` methods directly manipulate images (or stacks of images) using this phase offset method

```julia
julia> register(source, target; upsample_factor=100);

julia> coregister(cube; dims=3, refidx=1, upsample_factor=100);
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaHCI/SubpixelRegistration.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaHCI/SubpixelRegistration.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaHCI/SubpixelRegistration.jl/issues).
