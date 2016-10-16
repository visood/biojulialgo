#Iterator to produce kmers from a sequence
include("cluster.jl")

#Sense defines a biological property of a Nucleotide Acid Strand
#and should be used only as a property of a Nucleotide sequence,
#KmerItr will always iterate increasing index of a sequence,
#sense of the sequence can be handled when KmerItr is constructed

#A kmer iterator is defined over a sequence,
#and has a direction to iterate through the sequence.
immutable KmerItr{NS<:NucSeqLike}
    nseq::NS
    k::Integer
end

#some constructors for KmerItr
function kmerItr{
    S<:Strand,
    I<:Integer
    }(s::S, k::I)
    isForward(s) ?
    KmerItr(sequence(s), k) :
    KmerItr(reverse(sequence(s)), k)
end

kmerLength(kitr::KmerItr) = kitr.k
sequence(kitr::KmerItr) = kitr.nseq
seqtype{NS<:NucSeqLike}(kitr::KmerItr{NS}) =
    seqtype(NS)

Base.start{NS<:NucSeqLike}(kitr::KmerItr{NS}) =
    Base.start(sequence(kitr))

function Base.next{NS<:NucSeqLike}(kitr::KmerItr{NS}, p)
    seq = sequence(kitr)
    k = kmerLength(kitr)
    (Kmer(seq[p : p + k - 1], p), p + 1)
end

function Base.done{NS<:NucSeqLike}(kitr::KmerItr{NS}, pos)
    seq = sequence(kitr)
    Base.done(seq, pos + kmerLength(kitr) - 1)
end

function kmerCounts{NS<:NucSeqLike}(ktr::KmerItr{NS})
    T = seqtype(NS) #indirect seqtype could cause a problem
    kmerwc = Dict{T, Integer}()
    for kmer in ktr
        kmerwc[element(kmer)] = 1 + get(kmerwc, element(kmer), 0)
    end
    kmerwc
end

function countedAppearances{T}(wc::Dict{T, Integer})
    cw = Dict{Integer, Array{T, 1}}()
    for (w, k) in wc
        cw[k] = push!(get(cw, k, []), w)
    end
    cw
end

topFrequentKmers{NS<:NucSeqLike}(nseq::NS, k::Integer) =
    sort(collect(kmerCounts(KmerItr(nseq, k))),
         by = wk -> last(wk),
         rev = true)

topFrequentKmers{NS<:NucSeqLike}(nseq::NS, k::Integer, n::Integer) =
    topFrequentKmers(nseq, k)[1:n]

function kmerlocs{NS<:NucSeqLike}(ktr::KmerItr{NS})
    T = seqtype(NS) #indirect seqtype could cause a problem
    wocs = Dict{T, Vector{Integer}}()
    for kmer in ktr
        wocs[element(kmer)] = push!(get(wocs, element(kmer), []), position(kmer))
    end
    wocs
end
kmerlocs{NS<:NucSeqLike}(k::Integer, nseq::NS) = kmerlocs(KmerItr(nseq, k))

function kmerlocs{NS<:NucSeqLike}(pattern::NS, nseq::NS)
    k = length(pattern)
    kitr = KmerItr(nseq, k)
    occs = []
    for kmer in kitr
        if element(kmer) == pattern
            push!(occs, position(kmer))
        end
    end
    occs
end

function runningPatternCount{NS<:NucSeqLike}(pattern::NS, nseq::NS)
    k = length(pattern)
    kitr = KmerItr(nseq, k)
    counts = zeros(length(nseq) - k + 1)
    p = start(kitr)
    runcount = 0
    while !done(kitr, p)
        (kmer, p) = next(kitr, p)
        runcount += element(kmer) == pattern ? 1 : 0
        counts[p-1] = runcount
    end
    counts
end

function appxKmerLocs{NS<:NucSeqLike}(pattern::NS, nseq::NS, epsilon = 0)
    k = length(pattern)
    kitr = KmerItr(nseq, k)
    occs = []
    for kmer in kitr
        if (hammingDistance(element(kmer), pattern) <= epsilon)
            push!(occs, position(kmer))
        end
    end
    occs
end

function appxKmerLocs{NS<:NucSeqLike}(ktr::KmerItr{NS}, epsilon = 0)
    T = seqtype(NS)
    aklocs = Dict{T, Vector{Integer}}()
    klocs = kmerlocs(ktr)
    for w1 in keys(klocs)
        wocs = copy(klocs[w1])
        for w2 in keys(klocs)
            d = hammingDistance(w1, w2)
            if d != 0 && d <= epsilon
                append!(wocs, get(klocs, w2,  []) )
            end
        end
        aklocs[w1] = sort(wocs)
    end
    aklocs
end

appxKmerLocs{NS<:NucSeqLike}(k::Integer, nseq::NS, epsilon = 0) =
    appxKmerLocs(KmerItr(nseq, k), epsilon)
