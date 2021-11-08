import LinearMaps

struct ZeroMap{T,U,V} <: LinearMap{T}
    range::U
    domain::V
end

ZeroMap{T}(range::U, domain::V) where {T,U,V} = ZeroMap{T,U,V}(range, domain)
LinearMaps.MulStyle(A::ZeroMap) = LinearMaps.FiveArg()

Base.size(A::ZeroMap) = (length(A.range), length(A.domain),)

function LinearAlgebra.mul!(y::AbstractVector, L::ZeroMap, x::AbstractVector,
    α::Number, β::Number)

    y .*= β
end

function LinearAlgebra.mul!(y::AbstractVector, L::ZeroMap, x::AbstractVector)
    y .= 0
end

function Base.:(*)(A::ZeroMap, x::AbstractVector)
    T = eltype(A)
    y = similar(A.range, T)
    fill!(y, zero(T))
    LinearAlgebra.mul!(y,A,x)
end

Base.Matrix{T}(A::ZeroMap) where {T} = zeros(T, length(A.range), length(A.domain))