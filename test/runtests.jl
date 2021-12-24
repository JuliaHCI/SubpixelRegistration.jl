using PSFModels: Gaussian
using SubpixelRegistration
using Test
using StableRNGs

rng = StableRNG(121)
source = Gaussian(x=100.5, y=100.5, fwhm=10)[1:200, 1:200]

@testset "SubpixelRegistration.jl" begin
    @testset "pixel shift" for _ in 1:100
        shift = (rand(rng, -25:25), rand(rng, -25:25))
        shifted = fourier_shift(source, shift)
        result = phase_offset(source, shifted; upsample_factor=1)
        @test all(result.shift .≈ -1 .* shift)
    end
    @testset "subpixel shift" for _ in 1:100
        shift = 10 .* randn(rng, 2)
        shifted = fourier_shift(source, shift)
        
        f = 10
        result = phase_offset(source, shifted; upsample_factor=f)
        @test all(isapprox.(result.shift, -1 .* shift, atol=inv(f)))
    end
end
