#a cluster is an associative collection from positions to nucseqs
immutable ClusterEP{NS, I<:Integer} <: AbstractCluster{NS}
    idpos::I
    pos::Dict{NS, Vector{I}}
end

idpos(c::ClusterEP) = c.idpos
positions{NS, I<:Integer}(c::ClusterEP{NS, I}) =
    unique(collect(values(c.pos)))

size{NS, I<:Integer}(c::ClusterEP{NS, I}) =
    sum(map(Base.length, values(c.pos)))

Base.length{NS, I<:Integer}(c::ClusterEP{NS, I}) = size(c)

contains{NS, I<:Integer}(c::ClusterEP{NS, I}, x::I) =
    reduce(|, map(xs -> in(x, xs), values(c)))

elematdict{NS, I<:Integer}(c::ClusterEP{NS, I}) =
    Dict{I, NS}(foldl((d, xy) -> append!(d, map(y -> (y, first(xy)),
                                                last(xy))
                                         ),
                      Vector{Tuple{I, NS}}(),
                      collect(c.pos)
                      )
                )

function elemat{NS, I<:Integer}(c::ClusterEP{NS, I}, x::I)
    d = collect(c.pos)
    i = findfirst(xy -> in(x, xy[2]),  d)
    i > 0 ? Nullable{NS}(d[i][1]) : Nullable{NS}()
end

contains{NS, I<:Integer}(c::ClusterEP{NS, I}, ns::NS) = haskey(c.pos, ns)
posofdict(c::ClusterEP) = c.pos

posof{NS, I<:Integer}(c::ClusterEP{NS, I}, ns::NS) =
    haskey(c.pos, ns) ? Nullable{NS}(c.pos[ns]) : Nullable{NS}()

elements(c::ClusterEP) = keys(c.pos)

clusterEP{NS, I<:Integer}(km::NS, xs::Vector{I},
                        idposForCluster::Function) =
                ClusterEP{NS, I}(idposForCluster(xs),
                            Dict{NS, Vector{I}}(km => xs))

clusterEP{NS, I<:Integer}(km::NS, x::I) =
    ClusterEP{NS, I}(x, Dict{NS, Vector{I}}(km=>[x]))

function clusterEP{NS, I<:Integer}(kmxs::Dict{NS, Vector{I}},
                                 idposForCluster::Function = Base.minimum)
    i = idposForCluster(reduce((xs, ys) -> append!(xs, ys),
                               I[], values(kmxs))
                        )
    ClusterEP{NS, I}(i, kmxs)
end

function clusterEP{NS, I<:Integer}(elemat::Dict{I, NS},
                                 idposForCluster::Function = Base.minimum)
    posof = Dict{NS, Vector{I}}()
    for (x, ns) in elemat
        posof[ns] = haskey(posof, ns) ? push!(posof[ns], x) : I[x]
    end
    ClusterEP(idposForCluster(keys(elemat)), posof, idposForCluster)
end

function clusterEP{I <: Integer}(posof::Dict{I, Vector{I}},
                               idposForCluster::Function = Base.minimum)
    i = idposForCluster(reduce((xs, ys) -> append!(xs, ys),
                               I[], values(posof))
                        )
    ClusterEP(i, posof)
end

canMerge{NS, I<:Integer}(c1::ClusterEP{NS, I}, c2::ClusterEP{NS, I}) =
    true

function combineWith!{K, V}(d1::Dict{K, V}, d2::Dict{K, V}, p::Function)
    for (k, v) in d2
        d1[k] = p(d2[k],
                  haskey(d1, k) ? d1[k] : V[])
    end
end

function merge!{NS, I <: Integer}(c1::ClusterEP{NS, I},
                                  c2::ClusterEP{NS, I},
                                  idposForCluster = Base.minumum)

    if canMerge(c1, c2)
        cluster(foldl((d, xy) -> insert!(d, xy[1],
                                         append!(xy[2],
                                                 haskey(d, xy[1]) ?
                                                 d[xy[1]] :
                                                 I[])
                                         ),
                      c1.pos, collect(c2.pos)
                      ),
                idposForCluster
        )
    else
        c1
    end
end

function merged{NS, I <: Integer}(c1::ClusterEP{NS, I},
                                  c2::ClusterEP{NS, I},
                                  idposForCluster = Base.minumum)
    canMerge(c1, c2) ?
    Nullable(merge!(c1, c2)) :
    Nullable{ClusterEP{NS, I}}()
end

