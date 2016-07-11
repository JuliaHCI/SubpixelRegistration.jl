
<a id='subpixelRegistration-documentation-1'></a>

# subpixelRegistration documentation


<a id='User-functions-1'></a>

## User functions

<a id='subpixelRegistration.stackDftReg-Tuple{AbstractArray{T,N}}' href='#subpixelRegistration.stackDftReg-Tuple{AbstractArray{T,N}}'>#</a>
**`subpixelRegistration.stackDftReg`** &mdash; *Method*.



`stackDftReg{T,N,N1}(imgser::AbstractArray{T,N};ref::AbstractArray{T,N1}=reshape(slicedim(imgser,N,1),size(imgser)[1:(N-1)]),ufac::Int=10)`

`dftReg` applied to a full array. Each array along the last dimension of `imgser` is aligned to `ref` (by default the first image of the series, with precision `ufac`. Returns an array of `Dict` containing the translation information.

<a id='subpixelRegistration.alignFromDict-Tuple{AbstractArray{T,N},Array{Any,1}}' href='#subpixelRegistration.alignFromDict-Tuple{AbstractArray{T,N},Array{Any,1}}'>#</a>
**`subpixelRegistration.alignFromDict`** &mdash; *Method*.



`alignFromDft{T,N}(img2reg::AbstractArray{T,N},dftRegRes::Array{Any,1})`

Given an array and a `Dict` of translations as returned by `dftReg`, returns the aligned array.


<a id='Non-exported-functions-1'></a>

## Non exported functions

<a id='subpixelRegistration.subPixShift-Tuple{AbstractArray{Complex{Float64},N},Array{Float64,1}}' href='#subpixelRegistration.subPixShift-Tuple{AbstractArray{Complex{Float64},N},Array{Float64,1}}'>#</a>
**`subpixelRegistration.subPixShift`** &mdash; *Method*.



`subPixShift(imgft::AbstractArray{Complex{Float64}},shift::Array{Float64,1})`

Shift the image `imgft` (in Fourier space) by the amount provided in the vector `shift`.

<a id='subpixelRegistration.dftups' href='#subpixelRegistration.dftups'>#</a>
**`subpixelRegistration.dftups`** &mdash; *Function*.



`dftups{T,N}(inp::AbstractArray{T,N},no,usfac::Int=1,off=zeros(N))`

Upsampled DFT by matrix multiplication, computes an upsampled DFT in just a small region. `no` is the size of the region in pixels, `offset` the position in the full array, `usfac` the upsampling parameter.

<a id='subpixelRegistration.dftReg-Tuple{AbstractArray{Complex{Float64},N},AbstractArray{Complex{Float64},N},Any}' href='#subpixelRegistration.dftReg-Tuple{AbstractArray{Complex{Float64},N},AbstractArray{Complex{Float64},N},Any}'>#</a>
**`subpixelRegistration.dftReg`** &mdash; *Method*.



`dftReg{N}(imgRef::AbstractArray{Complex{Float64},N},imgF::AbstractArray{Complex{Float64},N},usfac)`

Main internal function : takes the Fourier transforms of the arrays to register as inputs (`imgRef`/`imgF`) and returns a dictionary containing the residual error and the shift along the dimensions of the arrays, with the level of subpixel precision provided by `usfac`

