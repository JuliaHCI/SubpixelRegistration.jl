


## Algorithm modified from the Matlab code accompanying  
## Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
##  "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).

function dftReg{N}(imgRef::AbstractArray{Complex{Float64},N},imgF::AbstractArray{Complex{Float64},N},usfac)
    if usfac==0
        ## Compute error for no pixel shift
        CCmax = sum(imgRef.*conj(imgF))
        rfzero = sumabs2(imgRef)
        rgzero = sumabs2(imgF)
        error = 1 - CCmax*conj(CCmax)/(rgzero*rfzero)
        #diffphase = atan2(imag(CCmax),real(CCmax))
        output = Dict("error" => error)
    elseif usfac==1
        ## Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the peak
        L = length(imgRef)
        CC = ifft(imgRef.*conj(imgF))
        loc = indmax(abs(CC))
        CCmax=CC[loc]
        rfzero = sumabs2(imgRef)/L
        rgzero = sumabs2(imgF)/L
        error = abs(1 - CCmax*conj(CCmax)/(rgzero*rfzero))
        #diffphase = atan2(imag(CCmax),real(CCmax))

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
        
        output = Dict("error"=>error,"shift"=>shift)
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
        loc = indmax(abs(CC))
        
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
            shift = round(Integer,shift*usfac)/usfac
            dftShift = div(ceil(usfac*1.5),2) ## center of output array at dftshift+1
            ## Matrix multiplies DFT around the current shift estimate  
            CC = conj(dftups(imgF.*conj(imgRef),ceil(Integer,usfac*1.5),usfac,dftShift-shift*usfac))/(prod(ind2)*usfac^N)
            ## Locate maximum and map back to original pixel grid
            loc = indmax(abs(CC))
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
        error = sqrt(abs(error))
        #diffphase = atan2(imag(CCmax),real(CCmax))
        ## If its only one row or column the shift along that dimension has no effect. Set to zero.
        shift[[div(x,2) for x in size(imgRef)].==1]=0
        
        output = Dict("error"=>error,"shift"=>shift)
       
    end
    output
end

### Upsampled DFT by matrix multiplication, can compute an upsampled DFT in just a small region
function dftups{T,N}(inp::AbstractArray{T,N},no,usfac::Int=1,off=zeros(N))
    sz = [size(inp)...]
    permV = 1:N
    for i in permV
        inp = permutedims(inp,[i;deleteat!(collect(permV),i)])
        kern = exp((-1im*2*pi/(sz[i]*usfac))*((0:(no-1))-off[i])*(ifftshift(0:(sz[i]-1))-floor(sz[i]/2)).')
        d = size(inp)[2:N]
        inp = kern * inp[:,:]
        inp = reshape(inp,(no,d...))
    end
    permutedims(inp,collect(ndims(inp):-1:1))
end


### Shift the image
function subPixShift(imgft::AbstractArray{Complex{Float64}},shift::Array{Float64,1})
    sz = [size(imgft)...]
    N=0
    for i in eachindex(sz)
        shifti = ifftshift((-div(sz[i],2)):(ceil(Integer,sz[i]/2)-1))*shift[i]/sz[i]
        resh = (repeat([1],inner=[i-1])...,length(shifti))
        N = N .- reshape(shifti,resh)
    end
    
    Greg = imgft .* exp(2im*pi*N)
    #Greg = Greg .* exp(1im*diffphase)
    Greg
end

### dftReg applied to a full array
function stackDftReg{T,N,N1}(imgser::AbstractArray{T,N};ref::AbstractArray{T,N1}=reshape(slicedim(imgser,N,1),size(imgser)[1:(N-1)]),ufac::Int=10)
    if N1 != N - 1
        error("Reference image has the wrong dimensionality")
    end
    ref = fft(ref)
    imgF = fft(imgser,(1:N1...))
    
    
    imgF = [reshape(slicedim(imgF,N,i),size(ref)) for i = 1:size(imgF)[N]]
    pmap(imgF) do im
        dftReg(ref,im,ufac)
    end
end
 
function alignFromDft{T,N}(img2reg::AbstractArray{T,N},dftRegRes::Array{Any,1})
    if length(dftRegRes) != size(img2reg)[N]
        error("Alignment results and image stack dimensionalities don't match.")
    end
    
    imRes = similar(img2reg,Float64)
    img2regF = fft(img2reg,(1:(N-1)...))
    strd = stride(imRes,N)

    szF = size(img2reg)[1:(N-1)]
    for i=1:size(img2reg)[N]
        frameft = subPixShift(reshape(slicedim(img2regF,N,i),szF),dftRegRes[i]["shift"])
        imRes[((i-1)*strd+1):(i*strd)] = real(ifft(frameft))
    end
    imRes
end



