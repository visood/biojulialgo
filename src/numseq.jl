function forwardDiff{T<:Number}(xs::Vector{T})
    ds = xs[2:end] - xs[1:end-1]
    push!(ds, -xs[end])
end

function backwardDiff{T<:Number}(xs::Vector{T})
    ds = [0]
    append!(ds, xs[2:end] - xs[1:end-1])
end


function monotonicity{T<:Number}(xs::Vector{T})
    ds = xs[2:end] - xs[1:end-1]
    apos = any(d -> d > 0, ds)
    aneg = any(d -> d < 0, ds)
    if apos
        aneg ? 0 : 1
    else
        aneg? -1 : 0
    end
end

