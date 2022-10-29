module Utilities

using StatsBase

export deaccumulate, multiaccumulate, momentmap

deaccumulate( op, v, init ) = [op( init, v[1] ); op.( v[1:end-1], v[2:end] )]
deaccumulate( op, v, ::Nothing ) = op.( v[1:end-1], v[2:end] )
deaccumulate( op, v; init=nothing ) = deaccumulate( op, v, init )

function multiaccumulate( f, v::AbstractVector{T}, new::AbstractVector{Bool}, init::AbstractVector{T} ) where {T}
    n = length(v)
    out = zeros( T, n )
    out[1] = f( init[1], v[1] )
    for i = 2:n
        out[i] = f( new[i] ? init[i] : out[i-1], v[i] )
    end
    return out
end

bound( x::AbstractVector{T}, buckets::Int ) where {T} = 
    bounds = quantile.( (x,), (d -> d:d:1-d)(1/buckets) )
    

label( x::AbstractVector{T}, bounds::Vector{T} ) where {T} = 
    searchsortedfirst.( (bounds,), x );

function momentmap( x::AbstractVector{T}, y::AbstractVector{U}, buckets::Int; moments = 1  ) where {T,U}
    good = .!isnan.(x) .& .!isnan.(y) .& .!isinf.(x) .& .!isinf.(y)
    x = x[good]
    y = y[good]
    n = length(x)

    bounds = bound( x, buckets )
    
    labels = label( x, bounds )

    counts = zeros( Int, buckets )
    xmoments = zeros( T, buckets )
    ymoments = zeros( U, buckets, moments )
    for i = 1:n
        l = labels[i]
        counts[l] += 1
        xmoments[l] += x[i]
        for j = 1:moments
            ymoments[l, j] += y[i]^j
        end
    end
    return Dict(
        :bounds => bounds,
        :counts => counts,
        :xmean => xmoments./counts,
        :ymoments => ymoments./counts,
    )
end

end # module
