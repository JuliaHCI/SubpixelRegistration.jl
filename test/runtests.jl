using ImageIO
using StableRNGs
using SubpixelRegistration
using Test
using TestImages

rng = StableRNG(121)

source = Float64.(testimage("cameraman"))

N_test = 100
upsample_factors = (5, 10, 100)

@testset "SubpixelRegistration.jl (normalize=$normalize)" for normalize in (false, true)

    @testset "pixel shift" for _ = 1:N_test
        shift = (rand(rng, -25:25), rand(rng, -25:25))
        shifted = @inferred fourier_shift(source, shift)
        result = @inferred phase_offset(
            source,
            shifted;
            upsample_factor = 1,
            normalize = normalize,
        )
        @test all(result.shift .≈ -1 .* shift)

        registered =
            @inferred register(source, shifted; upsample_factor = 1, normalize = normalize)
        @test registered ≈ fourier_shift(shifted, result.shift, result.phasediff)
    end

    @testset "subpixel shift (f=$f)" for _ = 1:N_test, f in upsample_factors
        shift = 10 .* randn(rng, 2)
        shifted = fourier_shift(source, shift)

        result = @inferred phase_offset(
            source,
            shifted;
            upsample_factor = f,
            normalize = normalize,
        )
        @test all(isapprox.(result.shift, -1 .* shift, atol = 5 * inv(f)))

        registered =
            @inferred register(source, shifted; upsample_factor = f, normalize = normalize)
        @test registered ≈ fourier_shift(shifted, result.shift, result.phasediff)
    end

    @testset "coregister (f=$f)" for f in upsample_factors
        shifts = 10 .* randn(rng, 19, 2)
        cube = cat(
            source,
            (fourier_shift(source, shift) for shift in eachrow(shifts))...,
            dims = 3,
        )
        cube_shift =
            @inferred coregister(cube; dims = 3, upsample_factor = f, normalize = normalize)
        for slice in eachslice(cube_shift, dims = 3)
            @test all(isapprox.(slice, source; atol = 5 * inv(f)))
        end

        # alternate axis of iteration
        cubep = permutedims(cube, (3, 2, 1))
        cubep_shift = @inferred coregister(
            cubep,
            dims = 1,
            upsample_factor = f,
            normalize = normalize,
        )
        for slice in eachslice(cubep_shift, dims = 1)
            @test all(isapprox.(slice', source; atol = 5 * inv(f)))
        end
    end

end
