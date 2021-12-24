using ImageIO
using SubpixelRegistration
using Test
using TestImages
using StableRNGs

rng = StableRNG(121)

source = Float64.(testimage("cameraman"))

@testset "SubpixelRegistration.jl" begin

    @testset "pixel shift" for _ in 1:100
        shift = (rand(rng, -25:25), rand(rng, -25:25))
        shifted = @inferred fourier_shift(source, shift)
        result = @inferred phase_offset(source, shifted; upsample_factor=1)
        @test all(result.shift .≈ -1 .* shift)

        registered = @inferred register(source, shifted; upsample_factor=1)
        @test registered ≈ fourier_shift(shifted, result.shift, result.phasediff)
    end

    @testset "subpixel shift" for _ in 1:100
        shift = 10 .* randn(rng, 2)
        shifted = fourier_shift(source, shift)

        f = 10
        result = @inferred phase_offset(source, shifted; upsample_factor=f)
        @test all(isapprox.(result.shift, -1 .* shift, atol=inv(f)))

        registered = @inferred register(source, shifted; upsample_factor=f)
        @test registered ≈ fourier_shift(shifted, result.shift, result.phasediff)
    end

    @testset "coregister" begin
        shifts = 10 .* randn(rng, 19, 2)
        cube = cat(source, (fourier_shift(source, shift) for shift in eachrow(shifts))..., dims=3)
        f = 100
        cube_shift = @inferred coregister(cube; dims=3, upsample_factor=f)
        for slice in eachslice(cube_shift, dims=3)
            @test all(isapprox.(slice, source; atol=1e-2))
        end

        # alternate axis of iteration
        cubep = permutedims(cube, (3, 2, 1))
        cubep_shift = @inferred coregister(cubep, dims=1, upsample_factor=f)
        for slice in eachslice(cubep_shift, dims=1)
            @test all(isapprox.(slice', source; atol=1e-2))
        end
    end

end
