module SubpixelRegistration

using AbstractFFTs
using FFTW

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
    shifts = maxidx.I
    shifts = @. ifelse(shifts > midpoints, shifts - shape, shifts) - 1.0

    @info "found shift $shifts" calculate_stats(maxima, source_freq, target_freq)...

    return shifts
end

function calculate_stats(crosscor_maxima, source_freq, target_freq)
    source_amp = sum(real.(source_freq .* source_freq'))
    source_amp /= length(source_freq)
    target_amp = sum(real.(target_freq .* target_freq'))
    target_amp /= length(target_freq)
    error = 1 - crosscor_maxima .* crosscor_maxima' / (source_amp * target_amp)
    phasediff = atan(imag(crosscor_maxima), real(crosscor_maxima))
    return (; error, phasediff)
end

end
