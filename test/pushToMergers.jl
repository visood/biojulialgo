function pushToMergers1{
    D1<:Located,
    D2<:Located
    } (cs::Vector{D1},
       x::D2,
       popit::Function = (x, y) -> true)
    print("merger to push ", x ,  " into ",
          cs, "\n")
    c = Base.pop!(cs)
    print("left  length ", length(cs), "\n")
    if popit(c, x)
        print("to pop merge\n")
        Base.push!(cs, c)
        print("put it back: ", map(u -> positions(u), cs), ".\n")
        print("length c ", length(c), "\t", "length x ", length(x), "\n")
        print("popmerge ", x, " into ", c, "\n")
        cnew = popmerge!(c, x, y -> popit(x, y))
        print(cnew, "\n")
    else
        print("not popped\n")
        print("concat ", x, " into ", c, "\n")
        cnew = concat!(c, x)
    end
    print("cnew: ", positions(cnew), ".\n")
    print("==============================\n")
    Base.push!(cs, cnew)
end

