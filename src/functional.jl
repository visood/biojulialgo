function takeWhile{T<:Any}(p::Function, xs::Vector{T})
    i = findfirst(x -> p(x), xs)
    i == 0 ? [] : xs[1:i-1]
end

function dropWhile{T <: Any}(p::Function, xs::Vector{T})
    i = findfirst(x -> p(x), xs)
    i == 0 ? T[] : xs[i:end]
end

function Base.insert!{A<:Any, B<:Any}(d::Dict{A, B},
                                      a::A, b::B)
    d[a] = b
    d
end

function Base.insert!{A<:Any, B<:Any}(d::Dict{A, B},
                                      ab::Vector{Tuple{A, B}})
    for (a, b) in ab
        insert!(d, a, b)
    end
    d
end

function insertWith!{A <: Any, B <: Any}(d::Dict{A, B}, a::A, b::B, p::Function)
    v = haskey(d, a) ? p(d[a], b) : b
    d[a] = v
    d
end

