module SubpixelRegistration

using AbstractFFTs
using FFTW
using Statistics

export phase_offset, register, coregister, coregister!, fourier_shift

"""
    phase_offset(source, target; upsample_factor=1)

Return the shift between `source` and `target` by measuring the maximum in the cross-correlation between the images. This algorithm can achieve `1/upsample_factor` precision by locally upsampling the cross-correlation via a matrix-multiplication DFT.

# Examples

# References

1. 
"""
function phase_offset(source, target; kwargs...)
    plan = plan_fft(source)
    return phase_offset(plan, plan * source, plan * target; kwargs...)
end

function phase_offset(plan, source_freq::AbstractMatrix{<:Complex{T}}, target_freq; upsample_factor=1) where T
    # whole-pixel shift
    # compute cross-correlation via iFFT
    image_product = source_freq .* target_freq'
    # phase normalization
    @. image_product /= max(abs(image_product), 100 * eps(T))
    cross_correlation = plan \ image_product # ifft

    # locate maximums
    maxima, maxidx = findmax(abs, cross_correlation)
    midpoints = map(ax -> (first(ax) + last(ax)) / 2, axes(source_freq))

    shape = size(source_freq)
    shifts = @. ifelse(maxidx.I > midpoints, maxidx.I - shape, maxidx.I) - 1.0

    @debug "found initial shift $shifts" calculate_stats(maxima, source_freq, target_freq)...

    isone(upsample_factor) && return shifts

    # upsample with matrix-multiply DFT
    shifts = @. round(shifts * upsample_factor) / upsample_factor
    upsample_region_size = ceil(upsample_factor * 1.5)
    # center of output array at dftshift + 1
    dftshift = div(upsample_region_size, 2)
    # matmul DFT
    sample_region_offset = @. dftshift - shifts * upsample_factor
    cross_correlation = upsampled_dft(image_product, upsample_region_size, upsample_factor, sample_region_offset)
    maxima, maxidx = findmax(abs, cross_correlation)
    shifts = @. shifts + (maxidx.I - dftshift - 1) / upsample_factor

    @debug "found final shift $shifts" calculate_stats(maxima, source_freq, target_freq)...

    return shifts
end

"""
    SubpixelRegistration.upsampled_dft(data, region_size, upsample_factor, offsets)

Calculate the cross-correlation in a region of size `region_size` via an upsampled DFT. The DFT uses matrix-multiplication to super-sample the input by `upsample_factor`. The frequencies will be shifted and centered around `offsets`.
"""
function upsampled_dft(data::AbstractMatrix{T}, region_size, upsample_factor, offsets) where {T<:Complex}
    im2pi = T(0, 2π)
    shiftrange = 1:region_size
    freqs = fftfreq(size(data, 2), inv(upsample_factor))
    kernel = @. exp(-im2pi * (shiftrange - offsets[2] - 1) * freqs')

    _data = kernel * data'

    freqs = fftfreq(size(data, 1), inv(upsample_factor))
    kernel = @. exp(im2pi * (shiftrange - offsets[1] - 1) * freqs')
    _data = kernel * _data'
    return _data
end

"""
    SubpixelRegistration.calculate_stats(crosscor_maxima, source_freq, target_freq)

Calculate the normalized root-mean-square error (NRMSE) and total phase difference between the two complex arrays, `source_freq` and `target_freq`, with maximum cross-correlation value `crosscor_maxima`
"""
function calculate_stats(crosscor_maxima, source_freq, target_freq)
    source_amp = mean(abs2, source_freq)
    target_amp = mean(abs2, target_freq)
    error = 1 - abs2(crosscor_maxima) / (source_amp * target_amp)
    phasediff = atan(imag(crosscor_maxima), real(crosscor_maxima))
    return (; error, phasediff)
end

"""
    fourier_shift(image, shift)

Shift the given `image` by `shift` along each axis, using the Fourier phase information.
"""
function fourier_shift(image, shift)
    FT = plan_fft(image)
    shifted = fourier_shift!(FT * image, shift)
    return real(FT \ shifted)
end

"""
    fourier_shift!(image_freq::AbstractMatrix{<:Complex}, shift)

Shift the given image, which is already in frequency-space, by `shift` along each axis. Modifies `image_freq` inplace.
"""
function fourier_shift!(image_freq::AbstractMatrix{<:Complex}, shift)
    shape = size(image_freq)
    
    freqs1 = fftfreq(shape[1])'
    freqs2 = fftfreq(shape[2])
    @. image_freq *= exp(-im * 2π * (freqs1 * shift[1] + freqs2 * shift[2]))
    return image_freq
end

function register(reference, target; kwargs...)
    plan = plan_fft(reference)
    target_freq = plan * target
    shift = phase_offset(plan, plan * reference, target_freq; kwargs...)
    shifted = fourier_shift!(target_freq, shift)
    return real(plan \ shifted)
end

function coregister!(cube::AbstractArray; dims, kwargs...)
    @inbounds for idx in axes(cube, dims)[begin + 1:end]
        reference = selectdim(cube, dims, idx - 1)
        target = selectdim(cube, dims, idx)
        # target is a view, update inplace
        target .= register(reference, target; kwargs...)
    end
    return cube
end
coregister(cube::AbstractArray; kwargs...) = coregister!(copy(cube); kwargs...)

end