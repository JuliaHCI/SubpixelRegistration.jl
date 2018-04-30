
"
`dftReg{N}(imgRef::AbstractArray{Complex{Float64},N},imgF::AbstractArray{Complex{Float64},N},usfac)`

Main internal function : takes the Fourier transforms of the arrays to register as inputs (`imgRef`/`imgF`) and returns a dictionary containing the residual error and the shift along the dimensions of the arrays, with the level of subpixel precision provided by `usfac`"
function dftReg{T,N}(imgRef::AbstractArray{Complex{T},N},imgF::AbstractArray{Complex{T},N},usfac)
    if usfac==0
        ## Compute error for no pixel shift
        CCmax = sum(imgRef.*conj(imgF))
        rfzero = sum(abs2, imgRef)
        rgzero = sum(abs2, imgF)
        error = 1 - CCmax*conj(CCmax)/(rgzero*rfzero)
        diffphase = atan2(imag(CCmax),real(CCmax))
        output = Dict("error" => error)
    elseif usfac==1
        ## Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the peak
        L = length(imgRef)
        CC = ifft(imgRef.*conj(imgF))
        loc = indmax(abs.(CC))
        CCmax=CC[loc]
        rfzero = sum(abs2, imgRef)/L
        rgzero = sum(abs2, imgF)/L
        error = abs(1 - CCmax*conj(CCmax)/(rgzero*rfzero))
        diffphase = atan2(imag(CCmax),real(CCmax))

        indi = size(imgRef)
        ind2 = tuple([div(x,2) for x in indi]...)
        
        locI = [ind2sub(indi,loc)...]
        
        shift = zeros(size(locI))
        for i in eachindex(locI)
            if locI[i]>ind2[i]
                shift[i]=locI[i]-indi[i]-1
            else shift[i]=locI[i]-1
            end
        end
        
        output = Dict("error"=>error,"shift"=>shift,"diffphase"=>diffphase)
    else
        ## Partial pixel shift
        
        ##First upsample by a factor of 2 to obtain initial estimate
        ##Embed Fourier data in a 2x larger array
        dimR = [size(imgRef)...]
        ranges = [(x+1-div(x,2)):(x+1+div(x-1,2)) for x in dimR]
        dimRL = map(x -> x*2, dimR)
        CC = zeros(Complex{Float64},tuple(dimRL...))
        CC[ranges...] = fftshift(imgRef).*conj(fftshift(imgF))
        
        ##  Compute crosscorrelation and locate the peak 
        CC = ifft(ifftshift(CC))
        loc = indmax(abs.(CC))
        
        indi = size(CC)
        locI = [ind2sub(indi,loc)...]
        CCmax = CC[loc]
        ## Obtain shift in original pixel grid from the position of the crosscorrelation peak 
        
        ind2 = tuple([div(x,2) for x in indi]...)
        
        shift = zeros(size(locI))
        for i in eachindex(locI)
            if locI[i]>ind2[i]
                shift[i]=locI[i]-indi[i]-1
            else shift[i]=locI[i]-1
            end
        end
        shift = shift/2

        ## If upsampling > 2, then refine estimate with matrix multiply DFT
        if usfac > 2
            ### DFT Computation ###
            # Initial shift estimate in upsampled grid
            shift = round.(Integer,shift*usfac)/usfac
            dftShift = div(ceil(usfac*1.5),2) ## center of output array at dftshift+1
            ## Matrix multiplies DFT around the current shift estimate  
            CC = conj(dftups(imgF.*conj(imgRef),ceil(Integer,usfac*1.5),usfac,dftShift-shift*usfac))/(prod(ind2)*usfac^N)
            ## Locate maximum and map back to original pixel grid
            loc = indmax(abs.(CC))
            locI = ind2sub(size(CC),loc)
            CCmax = CC[loc]
            rg00 = dftups(imgRef.*conj(imgRef),1,usfac)[1]/(prod(ind2)*usfac^N)
            rf00 = dftups(imgF.*conj(imgF),1,usfac)[1]/(prod(ind2)*usfac^N)
            locI = map((x) -> x - dftShift - 1,locI)
          
            for i in eachindex(shift)
                shift[i]=shift[i]+locI[i]/usfac
            end
          
        else  
            rg00 = sum(imgRef.*conj(imgRef))/(prod(indi))
            rf00 = sum(imgF.*conj(imgF))/(prod(indi))
        end
        error = 1 - CCmax*conj(CCmax)/(rg00*rf00)
        error = sqrt(abs.(error))
        diffphase = atan2(imag(CCmax),real(CCmax))
        ## If its only one row or column the shift along that dimension has no effect. Set to zero.
        shift[[div(x,2) for x in size(imgRef)].==1]=0
        
        output = Dict("error"=>error,"shift"=>shift,"diffphase"=>diffphase)
       
    end
    output
end

"
`dftups{T,N}(inp::AbstractArray{T,N},no,usfac::Int=1,off=zeros(N))`

Upsampled DFT by matrix multiplication, computes an upsampled DFT in just a small region. `no` is the size of the region in pixels, `offset` the position in the full array, `usfac` the upsampling parameter."
function dftups{T,N}(inp::AbstractArray{T,N},no,usfac::Int=1,offset=zeros(N))
    sz = [size(inp)...]
    permV = 1:N
    for i in permV
        inp = permutedims(inp,[i;deleteat!(collect(permV),i)])
        kern = exp.((-1im*2*pi/(sz[i]*usfac))*((0:(no-1))-offset[i])*(ifftshift(0:(sz[i]-1))-floor(sz[i]/2)).')
        d = size(inp)[2:N]
        inp = kern * reshape(inp, Val{2})
        inp = reshape(inp,(no,d...))
    end
    permutedims(inp,collect(ndims(inp):-1:1))
end


"
`subPixShift(imgft::AbstractArray{Complex{Float64}},shift::Array{Float64,1})`

Shift the image `imgft` (in Fourier space) by the amount provided in the vector `shift`."
function subPixShift{T}(imgft::AbstractArray{Complex{T}},shift::Array{Float64,1},diffphase)
    sz = [size(imgft)...]
    N=0
    for i in eachindex(sz)
        shifti = ifftshift((-div(sz[i],2)):(ceil(Integer,sz[i]/2)-1))*shift[i]/sz[i]
        resh = (repeat([1],inner=[i-1])...,length(shifti))
        N = N .- reshape(shifti,resh)
    end
    
    Greg = imgft .* exp.(2im*pi*N)
    Greg = Greg .* exp(1im*diffphase)
    Greg
end

"
`stackDftReg{T,N,N1}(imgser::AbstractArray{T,N};ref::AbstractArray{T,N1}=reshape(slicedim(imgser,N,1),size(imgser)[1:(N-1)]),ufac::Int=10)`

`dftReg` applied to a full array. Each array along the last dimension of `imgser` is aligned to `ref` (by default the first image of the series, with precision `ufac`. Returns an array of `Dict` containing the translation information."
function stackDftReg{T,N,N1}(imgser::AbstractArray{T,N};ref::AbstractArray{T,N1}=reshape(slicedim(imgser,N,1),size(imgser)[1:(N-1)]),ufac::Int=10)
    if ((N1 != (N - 1)) & (N1 != N))
        error("Reference image has the wrong dimensionality")
    end
    ref = fft(ref)
    imgF = fft(imgser,(1:N1...))
    
    if N1 == (N-1)
        imgF = [reshape(slicedim(imgF,N,i),size(ref)) for i = 1:size(imgF)[N]]
        results = pmap(imgF) do im
            dftReg(ref,im,ufac)
        end
    else
        results = dftReg(ref,imgF,ufac)
    end
    return results
end

"
`alignFromDft{T,N}(img2reg::AbstractArray{T,N},dftRegRes::Array{Any,1})`

Given an array and a `Dict` of translations as returned by `dftReg`, returns the aligned array."
function alignFromDict{T,N}(img2reg::AbstractArray{T,N},dftRegRes::Array{Dict{String,Any},1})
    if (length(dftRegRes) != size(img2reg)[N])
        error("Alignment results and image stack dimensionalities don't match.")
    end

    imRes = similar(img2reg,Float64)
    img2regF = fft(img2reg,(1:(N-1)...))
    strd = stride(imRes,N)
    szF = size(img2reg)[1:(N-1)]
    for i=1:size(img2reg)[N]
        frameft = subPixShift(reshape(slicedim(img2regF,N,i),szF),dftRegRes[i]["shift"],dftRegRes[i]["diffphase"])
        imRes[((i-1)*strd+1):(i*strd)] = real(ifft(frameft))
    end
    
    imRes
end


## Only a single Dict, means we expect one image of the same size
function alignFromDict{T,N}(img2reg::AbstractArray{T,N},dftRegRes::Dict)
    if (length(dftRegRes["shift"]) != N)
        error("Alignment results and image dimensionalities don't match.")
    end    
    frameft = real(ifft(subPixShift(fft(img2reg),dftRegRes["shift"],dftRegRes["diffphase"])))
end
