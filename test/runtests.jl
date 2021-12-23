using SubpixelRegistration
using Test
using StableRNGs

rng = StableRNG(121)

@testset "SubpixelRegistration.jl" begin
    reference = reshape(1.0:10000.0, 100, 100)
    @testset "pixel shift" for _ in 1:100
        shift = (rand(rng, -25:25), rand(rng, -25:25))
        shifted = fourier_shift(reference, shift)
        result = phase_offset(reference, shifted; upsample_factor=1)
        @test all(result.shift .≈ -1 .* shift)
    end
    @testset "subpixel shift" for _ in 1:100
        shift = 10 .* randn(rng, 2)
        shifted = fourier_shift(reference, shift)
        
        f = 10
        result = phase_offset(reference, shifted; upsample_factor=f)
        @test all(isapprox.(result.shift, -1 .* shift, atol=inv(f)))
    end
end
