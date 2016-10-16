immutable TakeWhileItr{T<:Any} <: AbstractVector{T}
    pred::Function
    sequence::Vector{T}
end

Base.start{T<:Any}(titr::TakeWhileItr{T}) =
    start(titr.sequence)

Base.next{T<:Any}(titr::TakeWhileItr{T}, pos::Integer) =
    next(titr.sequence, pos)

Base.done{T<:Any}(titr::TakeWhileItr{T}, pos::Integer) =
    (done(titr.sequence, pos) ||
     !titr.pred( first(next(titr.sequence, pos))))

Base.size{T<:Any}(twitr::TakeWhileItr{T}) =
    size(twitr.sequence)

Base.eltype{T<:Any}(::Type{TakeWhileItr{T}}) = T

takeWhile{T<:Any}(pred::Function, xs::Vector{T}) =
    collect(TakeWhileItr(pred, xs))

immutable DropWhileItr{T<:Any} <: AbstractVector{T}
    pred::Function
    sequence::Vector{T}
end

function Base.start{T<:Any}(ditr::DropWhileItr{T})
    titr = TakeWhileItr(ditr.pred, ditr.sequence)
    pos = start(titr)
    while !done(titr, pos)
        x, pos = next(titr, pos)
    end
    pos
end

Base.next{T<:Any}(ditr::DropWhileItr{T}, pos::Integer) =
    next(ditr.sequence, pos)

Base.done{T<:Any}(ditr::DropWhileItr{T}, pos::Integer) =
    done(ditr.sequence, pos)

Base.size{T<:Any}(dwitr::DropWhileItr{T}) =
    size(dwitr.sequence)

Base.eltype{T<:Any}(::Type{DropWhileItr{T}}) = T

dropWhile{T<:Any}(pred::Function, xs::Vector{T}) =
    collect(DropWhileItr(pred, xs))

