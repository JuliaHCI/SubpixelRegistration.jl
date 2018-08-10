
@afgc function _stackDftReg(resource::ArrayFireLibs,imgRef::AFArray{Complex{T},N1},imgF::AFArray{Complex{T},N},usfac) where {T,N,N1}
    if N1 != (N - 1)
        error("Reference image has the wrong dimensionality")
    end
    nimg = size(imgF,N)
    if usfac==0
        ## Compute error for no pixel shift
        @afgc CCmax = reshape(sum(imgRef.*conj(imgF),1:N1),nimg)
        @afgc rfzero = sum(abs2(imgRef))
        @afgc rgzero = reshape(sum(abs2(imgF),1:N1),nimg)
        @afgc error = 1 - CCmax.*conj(CCmax)./(rgzero*rfzero)
        @afgc diffphase = atan2(imag(CCmax),real(CCmax),false)
        output = error
    elseif usfac==1
        ## Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the peak
        L = length(imgRef)
        if N1 == 2
            @afgc CC = ifft2(imgRef.*conj(imgF))
        elseif N1 == 3
            @afgc CC = ifft3(imgRef.*conj(imgF))
        end
        
        @afgc CC =  abs(CC)
        CCmax = maximum(CC,1)
        
        for i in 2:N1
            CCmax = maximum(CCmax,i)
        end
        
        loc = find(CC .== CCmax) ## Problem if duplicates in a frame
        
        @afgc rfzero = sum(abs2(imgRef))/L
        @afgc rgzero = reshape(sum(abs2(imgF),1:N1),nimg)/L
        @afgc error = abs(1 - CCmax.*conj(CCmax)./(rgzero*rfzero))
        @afgc diffphase = Array(vec(atan2(imag(CCmax),real(CCmax),false)))
        
        indi = size(imgRef)
        ind2 = tuple([div(x,2) for x in indi]...)
        
        locTemp = ind2sub(size(CC),Int.(Array(loc)))
        
        locI = [tuple([locTemp[j][i] for j in 1:N1]...) for i in 1:nimg]
        
        shift = zeros(nimg,length(locI[1])) 
     
        for im in eachindex(locI)
            for i in eachindex(locI[im])
                if locI[im][i]>ind2[i]
                    shift[im,i]=locI[im][i]-indi[i]-1
                else shift[im,i]=locI[im][i]-1
                end
            end
        end
        
        #output = Dict("error"=>error,"shift"=>shift,"diffphase"=>diffphase)
        output= (shift,diffphase,error)
    else
        
        @afgc CC = imgRef.*conj(imgF)
        
        @afgc CC = padAFArray(CC,N1)
        
        #CC[ranges...] = imgRef.*conj(imgF)
        #afgc()
        ##  Compute crosscorrelation and locate the peak
        if N1 == 2
            @afgc CC = ifft2(CC)
        elseif N1 == 3
            @afgc CC = ifft3(CC)
        end

        @afgc CC =  abs(CC)
        CCmax = maximum(CC,1)
        
        for i in 2:N1
           @afgc CCmax = maximum(CCmax,i)
        end
        
        loc = find(CC .== CCmax)
        indi = size(CC)[1:N1]
        ind2 = tuple([div(x,2) for x in indi]...)
        locTemp = ind2sub(size(CC),Int.(Array(loc)))
        
        locI = [tuple([locTemp[j][i] for j in 1:N1]...) for i in 1:nimg]
        
        shift = zeros(nimg,length(locI[1])) 
     
        for im in eachindex(locI)
            for i in eachindex(locI[im])
                if locI[im][i]>ind2[i]
                    shift[im,i]=locI[im][i]-indi[i]-1
                else shift[im,i]=locI[im][i]-1
                end
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
            @afgc CC = conj(dftups(imgF.*conj(imgRef),ceil(Integer,usfac*1.5),usfac,dftShift-shift*usfac,N1))/(prod(ind2)*usfac^N1)
            ## Locate maximum and map back to original pixel grid
            
            @afgc CC =  abs(CC)
            CCmax = maximum(CC,1)
        
            for i in 2:N1
                @afgc CCmax = maximum(CCmax,i)
            end
        
            loc = find(CC .== CCmax)
            indi = size(CC)[1:N1]
            locTemp = ind2sub(size(CC),Int.(Array(loc)))
        
            locI = [tuple([locTemp[j][i] for j in 1:N1]...) for i in 1:nimg]
     
            rg00 = Array(dftups(imgRef.*conj(imgRef),1,usfac,zeros(1,N1),N1))[1]/(prod(ind2)*usfac^N)
            rf00 = Array(dftups(imgF.*conj(imgF),1,usfac,zeros(nimg,N1),N1))[repeat([1],inner=[N1])...,:]./(prod(ind2)*usfac^N)
            
            locI = map((x) -> x .- dftShift .- 1,locI)
            for im in eachindex(locI)
                for i in eachindex(locI[im])
                    shift[im,i]=shift[im,i]+locI[im][i]/usfac
                end
            end
          
        else  
            rg00 = Array(sum(imgRef.*conj(imgRef)))/(prod(indi))
            rf00 = Array(sum(imgF.*conj(imgF),1:N1))/(prod(indi))
        end
        error = Array(1 - vec(CCmax).*vec(conj(CCmax)))./vec(rg00*rf00)
        error = sqrt.(abs.(error))
        diffphase = Array(vec(atan2(imag(CCmax),real(CCmax),false)))
        ## If its only one row or column the shift along that dimension has no effect. Set to zero.
        shift[:,[div(x,2) for x in size(imgRef)].==1]=0
        output = (shift,diffphase,error)
    end
output
end

@afgc function padAFArray(CC,N)
    for r in N:-1:1
        sz = collect(size(CC))
        sz[r]=div(size(CC)[r],2)
        padR = zeros(AFArray{Complex{Float64}},sz...)
        @afgc CC = cat(r,padR,CC)
        @afgc CC = cat(r,CC,padR)
    end
    CC
end

function SubpixelRegistration.stackDftReg(resource::ArrayFireLibs,imgser::AbstractArray{T,N};ref::AbstractArray{T,N1}=reshape(selectdim(imgser,N,1),size(imgser)[1:(N-1)]),ufac::Int=10,chunkSize=100) where {T,N,N1}
    if ((N1 != (N - 1)) & (N1 != N))
        error("Reference image has the wrong dimensionality")
    end
    @afgc ref = fft(AFArray(ref))
    imgF = fft(imgser,(1:N1...,))
    
    if N1 == (N-1)
        chunks = [(i:min(chunkSize+i-1,size(imgF)[N])) for i in 1:chunkSize:size(imgF)[N]]
        shifts = zeros(size(imgF)[N],N1)
        diffphase = zeros(size(imgF)[N])
        error = zeros(size(imgF)[N])
        for ch in chunks
            shifts[ch,:],diffphase[ch],error[ch] = _stackDftReg(ArrayFireLibs(),ref,AFArray(selectdim(imgF,N,ch)),ufac)
        end
        results = (shifts,diffphase,error)
    else
        results = _stackDftReg(ArrayFireLibs(),ref,AFArray(imgF),ufac)
    end
    results = map(1:length(diffphase)) do i
        Dict("shift" => shifts[i,:],"diffphase" => diffphase[i],"error"=>error[i])
    end
    return results
end



"
`dftups{T,N}(inp::AbstractArray{T,N},no,usfac::Int=1,off=zeros(N))`

Upsampled DFT by matrix multiplication, computes an upsampled DFT in just a small region. `no` is the size of the region in pixels, `offset` the position in the full array, `usfac` the upsampling parameter."
@afgc function dftups(inp::AFArray,no,usfac::Int=1,offset=zeros(ndims(inp)),N=ndims(inp))
    sz = [size(inp)...]
    Nf = ndims(inp)
    nim = sz[Nf]
    permV = 1:N
    permFull = 0:3
    if Nf == N
        for i in permV
            pFull = collect(permFull)
            idxs = [pFull[i];deleteat!(pFull,i)]
            @afgc inp = reorder(inp,idxs...)
            leftK = AFArray((0:(no-1)).-offset[:,i]')
            @afgc leftK = reshape(leftK,no)
            rightK = AFArray(ifftshift(0:(sz[i]-1))-floor(sz[i]/2))
            @afgc rightK = rightK'
            kern = leftK*rightK
            @afgc kern = exp.(-1im*2*pi/(sz[i]*usfac)*kern)
            nsz=size(inp)
            d = nsz[2:Nf]
            #kern = reorder(AFArray(repeat(kern,outer=(1,1,nsz[1]))),0,2,1,4)
            @afgc inp = kern * reshape(inp, Val(2))
            @afgc inp = reshape(inp,no,d...)
        end
    else
        for i in permV
            pFull = collect(permFull)
            idxs = [pFull[i];deleteat!(pFull,i)]
            @afgc inp = reorder(inp,idxs...)
            leftK = AFArray((0:(no-1)).-offset[:,i]')
            @afgc leftK = reshape(leftK,no,1,nim)
            rightK = AFArray(repeat(ifftshift(0:(sz[i]-1)),outer=(1,sz[Nf]))-floor(sz[i]/2))
            @afgc rightK = reshape(rightK,1,sz[i],sz[Nf])
            
            kern = leftK*rightK
            @afgc kern = exp.(-1im*2*pi/(sz[i]*usfac)*kern)
            
            nsz=size(inp)
            d = nsz[2:Nf]
            remDim = prod(nsz[2:N])
            
            @afgc inp = kern * reshape(inp, nsz[1],remDim,sz[Nf])
            @afgc inp = reshape(inp,no,d...)
        end
    end
    newIdx = collect(N:-1:1).-1
    extraIdx = setdiff(0:3,newIdx)
    reorder(inp,newIdx...,extraIdx...)
end


"
`subPixShift(imgft::AbstractArray{Complex{Float64}},shift::Array{Float64,1})`

Shift the image `imgft` (in Fourier space) by the amount provided in the vector `shift`."
@afgc function subPixShift(imgft::AFArray{Complex{T},N},shift::Array{Float64,2},diffphase::Array{Float64,1}) where {T,N}
    sz = [size(imgft)[1:(N-1)]...]
    nim = size(imgft)[N]
    Z=0
    resh = [repeat([1],inner=N-1)...,nim]
    for i in eachindex(sz)
        shifti = AFArray(ifftshift((-div(sz[i],2)):(ceil(Integer,sz[i]/2)-1)).*shift[:,i]'/sz[i])
        reshTemp = copy(resh)
        reshTemp[i] = size(shifti)[1]
        Z = Z .- reshape(shifti,reshTemp...)
    end
    
    Greg = imgft .* exp.(2im*pi*Z)
    @afgc Greg = Greg .* reshape(exp.(AFArray(1im*diffphase)),resh...)
    Greg
end

function subPixShift(imgft::AFArray{Complex{T}},shift::Array{Float64,1},diffphase) where {T}
    sz = [size(imgft)...]
    N=0
    for i in eachindex(sz)
        shifti = ifftshift((-div(sz[i],2)):(ceil(Integer,sz[i]/2)-1))*shift[i]/sz[i]
        resh = (repeat([1],inner=[i-1])...,length(shifti))
        N = N .- reshape(shifti,resh...)
    end
    
    Greg = imgft .* exp.(2im*pi*N)
    Greg = Greg .* exp(1im*diffphase)
    Greg
end

function SubpixelRegistration.alignFromDict(resource::ArrayFireLibs,img2reg::AbstractArray{T,N},dftRegRes::Array{Dict{String,Any},1}) where {T,N}
    if (length(dftRegRes) != size(img2reg)[N])
        error("Alignment results and image stack dimensionalities don't match.")
    end

    img2regF = fft(img2reg,(1:(N-1)...,))
    szF = size(img2reg)[1:(N-1)]
    
    imRes = subPixShift(AFArray(img2regF),SubpixelRegistration.shifts_from_res(dftRegRes),SubpixelRegistration.diffphase_from_res(dftRegRes))

    imRes = real(ifft(Array(imRes),1:(N-1)))
    
    imRes
end
