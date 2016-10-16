# Iterator that produces the first Fibonacci numbers
immutable FibIt{T<:Integer}
    last2::Tuple{T,T}
    n::Integer
end

fibit(n::Integer) = FibIt((0,1), n)
fibit(n::Integer, T::Type) = FibIt{T}((0,1), n)

Base.start(fi::FibIt) = 1

function Base.next(fi::FibIt, state)
    if state == 1
        return (1,2)
    else
        fi.last2 = fi.last2[2], sum(fi.last2)
        (fi.last2[2], state + 1)
    end
end

Base.done(fi::FibIt, state) = state > fi.n

