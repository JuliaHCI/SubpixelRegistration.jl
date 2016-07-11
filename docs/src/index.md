# subpixelRegistration documentation

## User functions

```@docs
stackDftReg{T,N}(imgser::AbstractArray{T,N})
```

```@docs
alignFromDict{T,N}(img2reg::AbstractArray{T,N},dftRegRes::Array{Any,1})
```

## Non exported functions

```@docs
subpixelRegistration.subPixShift(imgft::AbstractArray{Complex{Float64}},shift::Array{Float64,1})
```

```@docs
subpixelRegistration.dftups{T,N}(inp::AbstractArray{T,N},no,usfac::Int=1,offset=zeros(N))
```

```@docs
subpixelRegistration.dftReg{N}(imgRef::AbstractArray{Complex{Float64},N},imgF::AbstractArray{Complex{Float64},N},usfac)
```
