module Utilities

using StatsBase
using LinearAlgebra

export deaccumulate, multiaccumulate, makegood, root, ↑

deaccumulate( op, v, init ) = [op( init, v[1] ); op.( v[1:end-1], v[2:end] )]
deaccumulate( op, v, ::Nothing ) = op.( v[1:end-1], v[2:end] )
deaccumulate( op, v; init=nothing ) = deaccumulate( op, v, init )

function multiaccumulate( f, v::AbstractVector{T}, new::AbstractVector{Bool}, init::AbstractVector{T} ) where {T}
    n = length(v)
    out = Vector( undef, n )
    out[1] = f( init[1], v[1] )
    for i = 2:n
        out[i] = f( new[i] ? init[i] : out[i-1], v[i] )
    end
    return out
end

multiaccumulate( f, v::AbstractVector{T}, new::AbstractVector{Bool} ) where {T} =
    multiaccumulate( f, v, new, v )

isgood(args...) = .&([.!isnan.(x) .& .!isinf.(x) for x in args]...)

makegood(args...) = getindex.( args, (isgood(args...),) )

function root( x, n::Int )
    @assert( ( x >= 0 ) || (n % 2 == 1) )
    return sign(x) * abs(x)^(1/n)
end

Base.:|>( i::T, v::AbstractVector{U} ) where {T <: Integer, U} = v[i]
Base.:|>( b::AbstractVector{Bool}, v::AbstractVector{U} ) where {U} = v[b]
Base.:|>( i::T, v::Tuple ) where {T <: Integer} = v[i]
Base.:|>( k::T, d::Dict{T,U} ) where {T, U} = d[k]

↑( f::Function, v::Vector{Vector{T}} ) where {T} = f .↑ v
↑( f::Function, v::Vector{T} ) where {T} = f.(v)

LinearAlgebra.adjoint( v::Vector{T} ) where {T} = reshape( v, 1, : )

# redefine to return one past end of list if not found
function Base.findfirst( b::AbstractVector{Bool} )
    i = 1
    while i <= length(b) && !b[i]
        i += 1
    end
    return i
end

# redefine to return one past end of list if not found
function Base.findlast( b::AbstractVector{Bool} )
    i = 1
    while i <= length(b) && !b[i]
        i += 1
    end
    return i
end

end # module
