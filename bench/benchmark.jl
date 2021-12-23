using BenchmarkTools
using CSV
using ProgressLogging
using PSFModels: Gaussian
using PythonCall
using Random
using Statistics
using SubpixelRegistration

skr = pyimport("skimage.registration")

rng = Xoshiro(1134)

# image dimension sizes
SIZES = 50:100:1000
UPSAMPLE_FACTORS = (1, 100, 1000)

rows = []
@progress "sizes" for N in SIZES
    center = (N + 1) / 2
    img = Gaussian(x=center, y=center, fwhm=10)[1:N, 1:N]
    shift = 10 .* randn(rng, 2)
    img_shifted = fourier_shift(img, shift)
    @progress "upsample factors" for ups in UPSAMPLE_FACTORS
        bench_jl = @benchmark phase_offset($img, $img_shifted, upsample_factor=$ups)
        bench_py = @benchmark skr.phase_cross_correlation($img, $img_shifted, upsample_factor=$ups, return_error=false)
        bundle = (;
            size=N, upsample_factor=ups, 
            time_jl=median(bench_jl).time/1e9, time_py=median(bench_py).time/1e9
        )
        @info "" bundle...
        push!(rows, bundle)
    end
end

outpath = joinpath(dirname(pathof(SubpixelRegistration)), "..", "bench", "benchmark_results.csv")
CSV.write(outpath, rows)
