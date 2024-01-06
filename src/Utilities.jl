module Utilities

using StatsBase

export deaccumulate, multiaccumulate, momentmap, label, bound, multimomentmap, makegood, root

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

bound( x::AbstractVector{T}, buckets::Int ) where {T} = 
    quantile.( (x,), (d -> d:d:1-d)(1/buckets) )

label( x::AbstractVector{T}, bounds::Vector{T} ) where {T} = 
    searchsortedfirst.( (bounds,), x );

label( x::AbstractVector{T}, buckets::Int ) where {T} = label( x, bound( x, buckets ) )

function momentmap( x::AbstractVector{T}, y::AbstractVector{U}, buckets::Int; moments = 1  ) where {T,U}
    (x,y) = makegood(x,y)
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

function root( x, n::Int )
    @assert( ( x >= 0 ) || (n % 2 == 1) )
    return sign(x) * abs(x)^(1/n)
end

function multimomentmap( x::AbstractVector{T}, y::AbstractVector{U}, z::AbstractVector{V}, xbuckets::Int, zbuckets::Int; moments::Int = 1 ) where {T,U,V}
    (x,y,z) = makegood(x,y,z)
    n = length(x)

    zbounds = bound( z, zbuckets )
    zlabels = label( z, zbounds )

    xbounds = bound( x, xbuckets )
    xlabels = label( x, xbounds )

    xmoments = zeros( T, zbuckets, xbuckets, moments + 1 )
    ymoments = zeros( U, zbuckets, xbuckets, moments + 1 )
    for i = 1:n
        lz = zlabels[i]
        lx = xlabels[i]
        for j = 0:moments
            xmoments[lz, lx, j+1] += x[i]^j
            ymoments[lz, lx, j+1] += y[i]^j
        end
    end
    for i in 1:zbuckets
        for j in 1:xbuckets
            for k = 2:moments+1
                xmoments[i,j,k] /= xmoments[i,j,1]
                ymoments[i,j,k] /= ymoments[i,j,1]
            end
            for k = moments+1:-1:3
                sumx = (-xmoments[i,j,2])^(k-1)
                sumy = (-ymoments[i,j,2])^(k-1)
                for l = 1:k-1
                    sumx += binomial( k-1, l ) * (-xmoments[i,j,2])^(k-1-l) * xmoments[i,j,l+1]
                    sumy += binomial( k-1, l ) * (-ymoments[i,j,2])^(k-1-l) * ymoments[i,j,l+1]
                end
                xmoments[i,j,k] = root( sumx, k-1 )
                ymoments[i,j,k] = root( sumy, k-1 )
            end
        end
    end
    return Dict(
        :zbounds => zbounds,
        :xbounds => xbounds,
        :xmoments => xmoments,
        :ymoments => ymoments,
    )
end

end # module
