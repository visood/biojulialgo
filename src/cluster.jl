import DataStructures
using Base.Order
#include("functional.jl")
include("dnaseq.jl")

abstract AbstractCluster{NS<:NucSeqLike, I<:Integer}

leftmost{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}) =
    Base.minimum(positions(c))
rightmost{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}) =
    Base.maximum(positions(c))

contains{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}, ns::NS) =
    in(ns, elements(c))

elematdict{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}) =
    dict(poselems(c))

posofdict{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}) =
    dict(elemposes(c))

function elemat{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}, x::I)
    elemats = poselems(c)
    i = findfirst(pe -> pe[1] == x, elemats)
    i == 0 ? nullable{ns}() : nullable{ns}(elemats[i][2])
end

function posof{NS<:NucSeqLike, I<:Integer}(c::AbstractCluster{NS, I}, ns::NS)
    elemats = poselems(c)
    i = findfirst(pe -> pe[2] == ns, elemats )
    i == 0 ? nullable{i}() : nullable{i}(elemats[i][1])
end

#subsequences are associated with positions. the simplest
#representation of such a subsequence is,
immutable Kmer{NS<:NucSeqLike, I<:Integer} <: AbstractCluster{NS, I}
    seq::NS
    pos::I
end

==(x::Kmer, y::Kmer) =
    x.seq == y.seq && x.pos == y.pos
idpos(km::Kmer) = km.pos
position(km::Kmer) = km.pos
element(km::Kmer) = km.seq
poselem(km::Kmer) = (km.pos, km.seq)
elempos(km::Kmer) = (km.seq, km.pos)
poselems(km::Kmer) =  [poselem(km)]
elemposes(km::Kmer) = [elempos(km)]
positions{NS<:NucSeqLike, I<:Integer}(km::Kmer{NS, I}) = I[km.pos]
elements{NS<:NucSeqLike, I<:Integer}(km::Kmer{NS, I}) = NS[km.seq]
Base.length(km::Kmer) = 1
size(::Kmer) = 1
leftmost(c::Kmer) = c.pos
rightmost(c::Kmer) = c.pos

Base.sort(km::Kmer) = km
Base.sort(kms::Vector{Kmer}) = Base.sort(kms, by = km -> position(km))

#a simple cluster to hold the positions of a sequence
#two consecutive positions of the sequence must be more than
#the sequence's length apart.
#how can we enforce this constraint?

immutable Clumer{NS<:NucSeqLike, I<:Integer} <: AbstractCluster{NS, I}
    idpos::I
    seq::NS
    pos::Vector{I}
end

idpos(c::Clumer) = c.idpos
positions(c::Clumer) = c.pos
Base.length(c::Clumer) = Base.length(c.pos)
size(c::Clumer) = Base.length(c.pos)

Base.sort{NS<:NucSeqLike, I<:Integer}(c::Clumer{NS, I}) =
    Clumer(c.idpos, c.seq, Base.sort(c.pos))

element(c::Clumer) = c.seq
elements{NS<:NucSeqLike, I<:Integer}(c::Clumer{NS, I}) = NS[c.seq]

function canMerge{
    NS<:NucSeqLike,
    I<:Integer,
    }(c1::Clumer{NS, I},
      c2::Clumer{NS, I})
    c1.seq == c2.seq
end

function clumer{
    NS<:NucSeqLike,
    I<:Integer
    }(s::NS, pos::Vector{I},
      idpos = Base.minimum)
    pos = length(pos) == 0 ? I[] : sort!(pos)
    idp = length(pos) == 0 ? typemin(I) : idpos(pos)
    Clumer(idp, s, pos)
end

function clumer{
    NS<:NucSeqLike,
    I<:Integer
    }(s::NS, pos::I)
    Clumer(pos, s, I[pos])
end

clumer(k::Kmer) = clumer(k.seq, k.pos)

function Base.push!{
    NS<:NucSeqLike,
    I<:Integer
    }(c::Clumer{NS, I}, x::I)
    clumer(c.seq, Base.push!(c.pos, x))
end

kmers{NS<:NucSeqLike, I<:Integer}(c::Clumer{NS, I}) =
    map(x -> Kmer(c.seq, x), c.pos)

#cluster should be position => element

immutable Cluster{NS<:NucSeqLike, I<:Integer} <: AbstractCluster{NS, I}
    idpos::I
    kmers::Vector{Kmer{NS, I}}
end

kmers(c::Cluster) = c.kmers
idpos(c::Cluster) = c.idpos
positions{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}) = map(position, c.kmers)
elements{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}) = map(element, c.kmers)
size{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}) = Base.length(c.kmers)
poselems{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}) = map(poselem, c.kmers)
elemposes{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}) = map(elempos, c.kmers)

contains{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}, ns::NS) =
    findfirst(km -> element(km) == ns, c.kmers) != 0

cluster{NS<:NucSeqLike, I<:Integer}(ns::NS, x::I) =
    Cluster{NS, I}(x, Kmer[Kmer(ns, x)])

function cluster{
    NS<:NucSeqLike,
    I<:Integer
    }(kms::Vector{Kmer{NS, I}},
      idposForCluster=Base.minimum)
    kms = Base.sort!(kms, by = position)
    Base.length(kms) > 0 ?
    Cluster(idposForCluster(map(position, kms)), kms) :
    Cluster(typemin(I), kms)
end

function cluster{
    NS<:NucSeqLike,
    I<:Integer
    }(cms::Clumer{NS, I},
      idposForCluster=Base.minimum)
    cluster(map(x -> Kmer{NS, I}(cms.seq, x), cms.pos),
            idposForCluster)
end

function cluster{
    NS<:NucSeqLike,
    I<:Integer
    }(posof::Dict{NS, Vector{I}},
      idposForCluster::Function = Base.minimum)
    cluster(foldl((d, eps) -> append!(d, [kmer(eps[1], p) for p in eps[2]]),
                  Vector{Kmer{NS, I}}(), posof)
            )
end

function cluster{
    NS<:NucSeqLike,
    I<:Integer
    }(elemat::Dict{I, NS},
      idposForCluster::Function = Base.minimum)
    Cluster(idposForCluster(keys(elemat)), collect(elemat))
end

function cluster{
    I<:Integer
    }(posof::Dict{I, Vector{I}},
      idposForCluster::Function = Base.minimum)
    cluster(foldl((d, eps) -> append!(d, [kmer(eps[1], p) for p in eps[2]]),
                  Vector{Kmer{NS, I}}(), posof)
            )
end

canMerge{NS<:NucSeqLike, I<:Integer}(c1::Cluster{NS, I}, c2::Cluster{NS, I}) = true

function concat!{
    NS<:NucSeqLike,
    I<:Integer
    }(c1::Cluster{NS, I}, c2::Cluster{NS, I})
    cluster(append!(c1.kmers, c2.kmers))
end

concat!{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}, k::Kmer{NS, I}) =
    cluster(push!(c.kmers, k))

concat!{NS<:NucSeqLike, I<:Integer}(c::Cluster{NS, I}, k::Clumer{NS, I}) =
    cluster(append!(c.kmers, map(x -> Kmer{NS, I}(k.seq, x), k.pos)))

function concat!{
    NS<:NucSeqLike,
    I<:Integer
    }(c::Clumer{NS, I}, d::Clumer{NS, I})
    if length(c.pos) == 0
        d
    else
        c.seq == d.seq ?
        clumer(c.seq, append!(c.pos, d.pos)) :
        concat!(cluster(c), cluster(d))
    end
   # Clumer(c.idpos, c.seq, append!(c.pos, d.pos))
end

function concat!{
    NS<:NucSeqLike,
    I<:Integer
    }(km1::Kmer{NS, I}, km2::Kmer{NS, I})
    km1.seq == km2.seq ?
    clumer(km1.seq, I[km1.pos, km2.pos]) :
    cluster(Kmer{NS, I}[km1, km2])
end

function concat!{
    NS<:NucSeqLike,
    I<:Integer
    }(c::Clumer{NS, I}, k::Kmer{NS, I})
    if length(c.pos) == 0
        clumer(k)
    else
        c.seq == k.seq ?
        clumer(c.seq, push!(c.pos, k.pos)) :
        cluster(push!(kmers(c), k))
    end
end

#some functional functions
function dropWhile{
    T<:Any
    }(dropit::Function, xs::Vector{T})
    i = Base.findfirst(x -> !dropit(x), xs)
    i == 0 ? T[] : xs[i:end]
end

function dropWhile{
    NS<:NucSeqLike,
    I<:Integer
    }(dropit::Function, cs::Clumer{NS, I})
    poses = dropWhile(dropit, positions(cs))
    clumer(element(cs), poses)
end

function dropWhile{
    NS<:NucSeqLike,
    I<:Integer
    }(dropit::Function, cs::Cluster{NS, I})
    nkmers = dropWhile(xy -> dropit(position(xy)),
                       kmers(cs))
    cluster(nkmers)
end

function dropWhile{
    NS<:NucSeqLike,
    I<:Integer
    }(dropit::Function, cs::Kmer{NS, I})
    cluster(dropit(position(cs))?
            Kmer{NS, I}[] : Kmer{NS, I}[cs])
end

function dropWhile{
    I<:Integer
    }(popit::Function, x::I)
    popit(x) ? I[] : I[x]
end
Located = Union{AbstractCluster, Integer, Vector{Integer}}

inProximity{I<:Integer}(L::I) =
    (x::Located, y::Located) -> distance(x, y) <= L
fartherThan{I<:Integer}(L::I) =
    (x::Located, y::Located) -> distance(x, y) > L

startpos(c::Located) = Base.minimum(positions(c))
endpos(c::Located) = Base.maximum(positions(c))
range(c::Located) = (startpos(c), endpos(c))
size(c::Located) = length(elements(c))
distance(x::Located, y::Located) = Base.abs(idpos(x) - idpos(y))

inProximity{I<:Integer}(L::I, x::Located, y::Located) = distance(x, y) < L
inProximity{I<:Integer}(L::I, x::Located) = y::Located -> distance(x, y) < L

#assuming that end >= start
inside(x::Located, y::Located) =
    (statrpos(y) >= startpos(x)) && (endpos(y) <= endpos(x))
inside(x::Located) = y::Located -> inside(x, y)

overlaps(x::Located, y::Located) =
    ! ( endpos(y) < startpos(x) | startpos(y) > endpos(x))
overlaps(x::Located) = y::Located -> overlaps(x, y)

function mergedGiven{
    D1<:Located,
    D2<:Located
    }(x::D1, y::D2,
      popif::Function = (x, y) -> true)
    canMerge(x, y) ?
    Nullable(popmerge!(x, y, z -> popif(y, z))) :
    Nullable{typeof(x)}()
end

function popmerge!{
    D1<:Located,
    D2<:Located
    }(x::D1, y::D2,
      popit::Function = z::Located -> true,
      mergefun::Function = concat!)
    mergefun(dropWhile(popit, x), y)
end

function mergeif!{
    D1<:Located,
    D2<:Located
    }(x::D1, y::D2,
      pred::Function = (x, y) -> true) pred(x, y) ?
    concat!(x, y) :
    popmerge!(x, y, x -> !pred(x, y))
end

function pushToMergers!{
    D1<:Located,
    D2<:Located
    }(cs::Vector{D1},
      x::D2,
      popit::Function = (x, y) -> true)
    c = Base.pop!(cs)
    if popit(c, x)
        Base.push!(cs, c)
        cnew = popmerge!(c, x, y -> popit(x, y))
    else
        cnew = concat!(c, x)
    end
    Base.push!(cs, cnew)
end

function clumps{
    D1<:Located,
    D2<:Located,
    I<:Integer
    }(L::I, t::I,
       positions::Vector{D1},
       start::Vector{D2} = D2[])
    filter(c -> size(c) > t,
           foldl((cs, x) -> pushToMergers!(cs, x, fartherThan(L)),
                 start, positions)
           )
end

function clumps{
    D<:Located,
    I<:Integer
    }(L::I, t::I,
       positions::Vector{D})
    clumps(L, t, positions[2:end], positions[1:1])
end

function clumps{
    D<:Integer,
    I<:Integer
    }(L::I, t::I,
       positions::Vector{D})
    clumps(L, t, positions[2:end], Vector{D}[positions[1:1]])
end

function clumps{
    NS<:NucSeqLike,
    I<:Integer
    }(L::I, t::I,
      positions::Vector{Kmer{NS, I}})
    clumps(L, t, positions[2:end], Clumer{NS, I}[clumer(positions[1])])
    #clumps(L, t, positions[2:end], Cluster{NS, I}[cluster(positions[1:1])])
end

function clumps{
    NS<:NucSeqLike,
    I<:Integer}(L::I, t::I,
                positions::Vector{Clumer{NS, I}})
    clumps(L, t, positions[2:end], positions[1:1])
    #clumps(L, t, positions[2:end], Cluster{NS, I}[cluster(positions[1])])
end

#function clumps{
    #I<:Integer
    #} (L::I, t::I, positions::Vector{I})
    #clumps(L, t, map(x -> I[x], positions))
#end

#for Integer and [Integer] to be treated as a Located
located{I<:Integer}(x::I) = x
idpos{I<:Integer}(x::I) = Integer(x)
positions{I<:Integer}(x::I) = Integer[x]
startpos{I<:Integer}(x::I) = Integer(x)
endpos{I<:Integer}(x::I) = Integer(x)
elements{I<:Integer}(x::I) = Integer[x]
canMerge{I<:Integer}(x::I, y::I) = true

located{I<:Integer}(xs::Vector{I}) = Base.sort!(xs)
idpos{I<:Integer}(xs::Vector{I}) = xs[1]
positions{I<:Integer}(xs::Vector{I}) = xs
startpos{I<:Integer}(xs::Vector{I}) = Base.minimum(xs)
endpos{I<:Integer}(xs::Vector{I}) = Base.maximum(xs)
elements{I<:Integer}(xs::Vector{I}) = xs
canMerge{I<:Integer}(xs::Vector{I}, ys::Vector{I}) = true

canMerge{I<:Integer}(xs::Vector{I}, y::I) = true
canMerge{I<:Integer}(x::I, ys::Vector{I}) = true
concat!{I<:Integer}(x::I, y::I) = I[x, y]
concat!{I<:Integer}(x::I, ys::Vector{I}) = Base.insert!(ys, 1, x)
concat!{I<:Integer}(xs::Vector{I}, y::I) = Base.push!(xs, y)
concat!{I<:Integer}(xs::Vector{I}, ys::Vector{I}) = Base.append!(xs, ys)

