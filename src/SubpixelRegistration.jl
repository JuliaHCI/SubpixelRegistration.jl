module SubpixelRegistration

using ComputationalResources, FFTW, Distributed

function __init__()
    # Enable `using` to load additional modules in this folder
    push!(LOAD_PATH, dirname(@__FILE__))
    # Now check for any resources that your package supports
    if haveresource(ArrayFireLibs)
        # User has indicated support for the ArrayFire libraries, so load your relevant code
        @eval using SubpixelRegistrationAF
    end
    # Put additional resource checks here
    # Don't forget to clean up!
    pop!(LOAD_PATH)
end

include("dftReg.jl")

export stackDftReg,
alignFromDict,
subPixShift

end # module
