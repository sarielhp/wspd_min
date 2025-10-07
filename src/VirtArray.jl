#
module VirtArray

struct VArray{T} <: AbstractVector{T}
    vec::Vector{T};
    loc::Vector{Int64};
end

@inline function  Base.getindex(a::VArray{T}, i::Int) where {T}
    return   a.vec[ a.loc[ i ] ];
end

@inline function  Base.size(va::VArray{T}) where {T}
    return  size( va.vec );
end

function  orig_index( va::VArray{T}, i::Int64 )::Int where {T}
    return  va.loc[ i ];
end

function  swap!(va::VArray{T}, i::Int, j::Int) where {T}
   (va.loc[i], va.loc[j]) = (va.loc[j], va.loc[i])
end

function VArray(vec_a::Vector{T}) where {T}
    return VArray{T}(vec_a, [i for i in 1:length(vec_a)])
end

function Base.first(va::VArray{T}) where{T}
    return   va.vec[ va.loc[ 1 ] ]
end

function Base.last(va::VArray{T}) where {T}
    return   va.vec[ va.loc[ end ] ]
end

function Base.length(va::VArray{T}) where {T}
    return   length(va.loc)
end

function Base.eachindex(va::VArray{T}) where{T}
    return   eachindex( va.loc )
end

@inline function Base.iterate(va::VArray{T} ) where{T}
    return   iterate(va, 1)
end

@inline function Base.iterate(va::VArray{T}, i::Int) where {T}
    if i > length(va.loc)
        return   nothing
    else
        return   va.vec[ loc[ i ] ], i + 1
    end
end

export VArray, swap!, orig_index;

end  ## Module

# Example usage:
# va = VArray((10,), i -> i[1]^2)
# println(va[3])  # Output: 9
