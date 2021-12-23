using SubpixelRegistration
using Test
using StableRNGs

rng = StableRNG(121)

@testset "SubpixelRegistration.jl" begin
    reference = reshape(1.0:10000.0, 100, 100)
    @testset "pixel shift" for _ in 1:100
        shift = (rand(rng, -25:25), rand(rng, -25:25))
        shifted = SubpixelRegistration.fourier_shift(reference, shift)
        _shift = phase_offset(reference, shifted; upsample_factor=1)
        @test all(_shift .≈ -1 .* shift)
    end
    @testset "subpixel shift" for _ in 1:100
        shift = (
            10 * randn(rng),
            10 * randn(rng)
        )
        shifted = SubpixelRegistration.fourier_shift(reference, shift)
        
        f = 10
        _shift = phase_offset(reference, shifted; upsample_factor=f)
        @test all(isapprox.(_shift, -1 .* shift, atol=inv(f)))
    end
end
