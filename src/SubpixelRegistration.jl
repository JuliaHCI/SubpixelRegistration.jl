module SubpixelRegistration

using AbstractFFTs
using FFTW
using Statistics

function phase_register(source, target; upsample_factor=1)

    # 1. FFT input data
    source_freq = fft(source)
    target_freq = fft(target)

    # whole-pixel shift
    # compute cross-correlation via iFFT
    image_product = source_freq .* target_freq'
    # phase normalization
    @. image_product /= max(abs(image_product), 100 * eps(eltype(source)))
    cross_correlation = ifft(image_product)

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

function upsampled_dft(data::AbstractArray{T,N}, upsample_region_size, upsample_factor, axis_offsets) where {T<:Complex,N}
    im2pi = T(0, 2π)
    shiftrange = 1:upsample_region_size
    freqs = fftfreq(size(data, 2), inv(upsample_factor))
    kernel = @. exp(-im2pi * (shiftrange - axis_offsets[2] - 1) * freqs')

    _data = kernel * data'

    freqs = fftfreq(size(data, 1), inv(upsample_factor))
    kernel = @. exp(-im2pi * (shiftrange' - axis_offsets[1] - 1) * freqs)
    _data = kernel' * _data'
    return _data
end

function calculate_stats(crosscor_maxima, source_freq, target_freq)
    source_amp = mean(abs2, source_freq)
    target_amp = mean(abs2, target_freq)
    error = abs(1 - crosscor_maxima * crosscor_maxima' / (source_amp * target_amp))
    phasediff = atan(imag(crosscor_maxima), real(crosscor_maxima))
    return (; error, phasediff)
end

end
