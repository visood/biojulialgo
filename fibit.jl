#Iterator produces the first n Fibonacci numbers
immutable FibIt
    max::Integer
end

fibit(n::Integer) = FibIt(n)

Base.start(fi::FibIt) = (0, 1)

function Base.next(fi::FibIt, state)
    state[1], (state[2], sum(state))
end

Base.done(fi::FibIt, state) = state[1] > fi.max

for i = fibit(10)
    print(i, " ")
end
