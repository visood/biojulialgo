#Iterator to produce kmers from a sequence
immutable KmerItr{T<:AbstractString}
    sequence::T
    k::Integer
end

kmerItr{T<:AbstractString}(seq::T,
                           n::Integer) = KmerItr(seq, n)

Base.start(kt::KmerItr) = 1

function Base.next(kt::KmerItr, state)
    kt.sequence[state:state + kt.k - 1], state + 1
end

Base.done(kt::KmerItr, state) =
    state + kt.k - 1 > length(kt.sequence)


function kmerCounts{T<:AbstractString}(ktr::KmerItr{T})
    wc = Dict{T, Integer}()
    for w in ktr
        wc[w] = 1 + get(wc, w, 0)
    end
    wc
end

kmerCounts{T<:AbstractString}(s::T, k::Integer) =
    kmerCounts(kmerItr(s, k))

function countedAppearances{T}(wc::Dict{T, Integer})
    cw = Dict{Integer, Array{T, 1}}()
    for (w, k) in wc
        cw[k] = push!(get(cw, k, []), w)
    end
    cw
end

function topFrequentKmers{T<:AbstractString}(s::T, k::Integer)
     sort(collect(kmerCounts(kmerItr(s, k))),
          by = wk -> last(wk),
          rev=true)
end

topFrequentKmers{T<:AbstractString}(s::T, k::Integer, n::Integer) =
    topFrequentKmers(s, k)[1:n]

function kmerOccurences{T<:AbstractString}(ktr::KmerItr{T})
    wocs = Dict{T, Vector{Int64}}()
    pos = start(ktr)
    while !done(ktr, pos)
        kmer, posnext = next(ktr, pos)
        wocs[kmer] = push!(get(wocs, kmer, []), pos)
        pos = posnext
    end
    wocs
end

function occurences{T<:AbstractString}(pattern::T, text::T)
    k = length(pattern)
    kitr = KmerItr(text, k)
    occs = []
    pos = start(kitr)
    while !done(kitr, pos)
        kmer, posnext = next(kitr, pos)
        if kmer == pattern
            push!(occs, pos)
        end
        pos = posnext
    end
    occs
end

function clumpAppend(L::Int64, c::Vector{Int64},
                     cs::Vector{Vector{Int64}},
                     x::Int64)
    if x < first(c) + L
        push!(c, x)
    else
        cnew = push!(filter(y -> x < y + L, c), x)
        push!(cs, c)
        c = cnew
    end
    (c, cs)
end

abstract AbstractCluster{T<:Any}

immutable Cluster{T<:Any} <: AbstractCluster{T}
    idpos::Integer
    elementAt::Dict{Integer, T}
    positionsOf::Dict{T, Vector{Integer}}
end

positions{T<:Any}(c<:AbstractCluster{T}, x::T) =
    get(c.positionsOf, x, [])

elements{T<:Any}(c<:AbstractCluster{T}, x::Integer) =
    c.elementAt[x]

startPos(c<:AbstractCluster) =
    minimum(c.keys(elementAt))

endPos(c<:AbstractCluster) =
    maximum(c.keys(elementAt))

range(c<:AbstractCluster) =
    (startPos(c), endPos(c))

size(c<:Cluster) = length(c.elementAt)

idposLeftMost = Base.minimum

cluster{T<:Any}(km::T,
                xs::Vector{Integer},
                idposForCluster::Function) =
                    Cluster{T}(idposForCluster(xs),
                               Dict([(x, km) for x in xs]),
                               Dict(km => xs))

function cluster{T<:Any}(kmxs::Dict{T, Vector{Integer}},
                         idposForCluster::Function = idposLeftMost)
    elementAt = Dict{Integer, T}()
    for t, xs in kmxs
        for x in xs
            elementAt[x] = t
        end
    end
    Cluster{T}(idposForCluster(keys(d)),
               elementAt, kmxs)
end




#Clump = Vector{Int64}

idpos(c::Clump) = c.idpos

distance(c<:AbstractClump, d::AbstractClump) =
    Base.abs(idpos(c) - idpos(d))

distance(c<:AbstractClump, x::Integer) =
    Base.abs(idpos(c) - x)

distance(x::Integer, c<:AbstractClump) =
    distance(c, x)

within(c::AbstractClump, x::Integer) =
    (x >= c.startPos) && (x <= c.startPos + c.range)


clump(x::Int64) = Clump([x])

distance(c::Clump, x::Int64) = abs(x - id(c))
distance(x::Int64, c::Clump) = distance(c, x)
distance(c::Clump, d::Clump) = abs(id(c) - id(d))
distance(x::Int64, y::Int64) = abs(y - x)

inProximacy(L::Int64) = x, y -> distance(x, y) < L
inProximacy(L::Int64, x::Int64) = y -> distance(x, y) < L
inProximacy(L::Int64, c::Clump) = inProximacy(L, id(c))

function pushToClump!(L::Int64,
                      cs::Vector{Clump},
                      x::Int64)
    c = pop!(cs)
    if distance(x, c) < L
        push!(c, x)
    else
        push!(cs, c)
        cnew = push!(dropWhile(inProximacy(L, x), c), x)
        c = cnew
    end
    push!(cs, c)
end


clumps(L::Int64, t::Int64, positions::Vector{Int64}) =
    filter(c -> length(c) >= t,
           foldl( (cs, x) -> pushToClump!(L, cs, x),
                  Vector{Int64}[ [positions[1]]],
                  positions[2:end]
                  )
           )

function clumpyKmers{T<:AbstractString}(L::Int64,
                                        t::Int64,
                                        kitr::KmerItr{T})
    kocs = filter((k,v) -> length(v) >= t && (last(v) - first(v) <= L),
                  kmerOccurences(kitr) )
    collect(keys(kocs))
end

egseq1 = "ababbbababaaab"
xs = [1,2,3,4,7,8,9,16,17,18]
