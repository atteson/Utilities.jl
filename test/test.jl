using Utilities
using Random
using StatsBase

Random.seed!(1)

x = randn( 10_000 );
y = randn( 10_000 );
z = randn( 10_000 );

xbuckets = 10
zbuckets = 5
moments = 4
d = multimomentmap( x, y, z, xbuckets, zbuckets, moments=moments )

xbounds = quantile.( [x], (1:xbuckets-1)./xbuckets );
@assert( maximum( abs.( xbounds - d[:xbounds] ) ) < 1e-8 )

zbounds = quantile.( [z], (1:zbuckets-1)./zbuckets );
@assert( maximum( abs.( zbounds - d[:zbounds] ) ) < 1e-8 )

xlabels = searchsortedfirst.( [xbounds], x );
zlabels = searchsortedfirst.( [zbounds], z );

for i = 1:zbuckets
    for j = 1:xbuckets
        good = (zlabels .== i) .& (xlabels .== j)
        n = sum(good)
        @assert( d[:xmoments][i,j,1] == n, "Error for i=$i, j=$j, k=1" )
        @assert( d[:ymoments][i,j,1] == n, "Error for i=$i, j=$j, k=1" )
        
        xs = x[good]
        ys = y[good]
        
        xmean = mean( xs )
        @assert( abs( d[:xmoments][i,j,2] - xmean ) < 1e-8 )
        ymean = mean( ys )
        @assert( abs( d[:ymoments][i,j,2] - ymean ) < 1e-8 )

        for k = 3:moments+1
            @assert( abs( d[:xmoments][i,j,k] - root( mean((xs .- xmean).^(k-1)), k-1 ) ) < 1e-8, "Error for i=$i, j=$j, k=$k" )
            @assert( abs( d[:ymoments][i,j,k] - root( mean((ys .- ymean).^(k-1)), k-1 ) ) < 1e-8, "Error for i=$i, j=$j, k=$k" )
        end
    end
end
            
