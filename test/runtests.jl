using SubpixelRegistration
using FFTW
using Test

## Testing the registration on a small 4d array
@testset "Testing the registration on a small 4d array" begin
    test4d = zeros(40,40,20,2)

    test4d[10:20,10:20,2:10,1] .= 1
    test4d[5:15,15:25,5:13,2] .= 1


    ## No subpixel resolution
    dftResults = stackDftReg(test4d,ufac=1)

    @test length(dftResults) == size(test4d)[4]
    @test dftResults[2]["shift"] == [5,-5,-3]

    ## Factor 2
    dftResults = stackDftReg(test4d,ufac=2)

    @test length(dftResults) == size(test4d)[4]
    @test dftResults[2]["shift"] == [5,-5,-3]

    # More refined (needs dftups)
    dftResults = stackDftReg(test4d,ufac=3)

    @test length(dftResults) == size(test4d)[4]
    @test dftResults[2]["shift"] == [5,-5,-3]

    ## Testing the translation
    back4d = alignFromDict(test4d,dftResults)
    @test all(round.(back4d[:,:,:,2]-back4d[:,:,:,1]).<eps())
end
