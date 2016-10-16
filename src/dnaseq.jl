using Match
import Base.map

abstract Nucleotide

convert{C<:Nucleotide}(ASCIIString, ns::Vector{C}) =
    Base.join(map(n -> symbol(n), ns))

asString{C<:Nucleotide}(ns::Vector{C}) =
    Base.join(map(n -> symbol(n), ns))

#we can use a Char as a Nucleotide

symbol(c::Char) = c
complement(c::Char) = @match c begin
    'A' => 'T'
    'C' => 'G'
    'T' => 'A'
    'G' => 'C'
    'U' => 'C'
    'N' => 'N'
    '=' => '='
    _   => '_'
end

immutable Nucleochar <: Nucleotide
    symbol::Char
end

symbol(c::Nucleochar) = c.symbol

_A_  = Nucleochar('A')
_C_  = Nucleochar('C')
_G_  = Nucleochar('G')
_T_  = Nucleochar('T')
_U_  = Nucleochar('U')
_N_  = Nucleochar('N')
_eq_ = Nucleochar('=')
_X_  = Nucleochar('_')

complement(b::Nucleochar) = @match b begin
    Nucleochar('A') => Nucleochar('T')
    Nucleochar('C') => Nucleochar('G')
    Nucleochar('G') => Nucleochar('C')
    Nucleochar('T') => Nucleochar('A')
    Nucleochar('U') => Nucleochar('C')
    Nucleochar('N') => Nucleochar('N')
    Nucleochar('=') => Nucleochar('=')
                  _ => Nucleochar('_')
end

symbol(c::Char) = c

nucleochar(x::Char) = @match x begin
    'A' => _A_
    'C' => _C_
    'T' => _T_
    'G' => _G_
    'U' => _U_
    'N' => _N_
    '=' => _eq_
    'a' => _A_
    'c' => _C_
    't' => _T_
    'g' => _G_
    'u' => _U_
    'n' => _N_
    _   => _X_
end

#an exercise in how two bits can be a nucleotide
#two bits can cover only 4 possibilties, so we have only ACGT
#we could use three bits for ACGT + UN=X
immutable NucleoTwoBool <: Nucleotide
    value::Tuple{Bool, Bool}
end

symbol(c::NucleoTwoBool) = @match c begin
    NucleoTwoBool((false, false)) => 'A'
    NucleoTwoBool((false, true))  => 'C'
    NucleoTwoBool((true, false))  => 'G'
    NucleoTwoBool((true, true))   => 'T'
end

complement(n::NucleoTwoBool) =
    NucleoTwoBool((!first(n.value), !last(n.value)))

nucleotwobool(c::Char) = @match c begin
    'A' => NucleoTwoBool((false, false))
    'C' => NucleoTwoBool((false, true))
    'G' => NucleoTwoBool((true, false))
    'T' => NucleoTwoBool((true, true))
    'a' => NucleoTwoBool((false, false))
    'c' => NucleoTwoBool((false, true))
    'g' => NucleoTwoBool((true, false))
    't' => NucleoTwoBool((true, true))
end


#to use all Nucleotide like types together
NucleotideLike = Union{Nucleotide, Char}

#and we can also convert
convert(Char, n::Nucleochar) = symbol(n)
convert(Nucleochar, c::Char) = nucleochar(c)

#Strands of DNA/RNA have a sense,
abstract Sense
immutable Forward <: Sense end
immutable Reverse <: Sense end
immutable SenseLess <: Sense end

isForward(d::Forward) = true
isForward(d::Reverse) = false
isForward(d::SenseLess) = true
isReverse(d::Forward) = false
isReverse(d::Reverse) = true
isReverse(d::SenseLess) = true
anti(d::Forward) = Reverse()
anti(d::Reverse) = Forward()
anti(d::SenseLess) = SenseLess()

#Nucleotides can be put in a sequence
#NucSeq does not have a sense (implicitlt resolved to Forward)
abstract NucSeq

==(x::NucSeq, y::NucSeq) =
    (sequence(x) == sequence(y))

isEmpty{NS<:NucSeq}(nseq::NS) = length(nseq) == 0
Base.start{NS<:NucSeq}(nseq::NS) =
    Base.start(sequence(nseq))

Base.next{NS<:NucSeq}(nseq::NS, pos) =
    Base.next(sequence(nseq), pos)

Base.done{NS<:NucSeq}(nseq::NS, pos) =
    Base.done(sequence(nseq), pos)

Base.getindex{NS<:NucSeq}(ns::NS, posr::UnitRange{Int64}) =
    sequence(ns)[posr]

Base.length{NS<:NucSeq}(ns::NS) = Base.length(sequence(ns))

map{NS<:NucSeq}(func::Function, ns::NS) =
    nucseq(map(func, sequence(ns)))

complement(ns) = map(complement, ns)

#treat an ASCIIString as a NucSeq
isEmpty(ns::ASCIIString) = length(ns) == 0
sequence(ns::ASCIIString) = ns
seqtype(ns::ASCIIString) = ASCIIString
seqtype(::Type(ASCIIString)) = ASCIIString
seqtype(ASCIIString) = ASCIIString
nucseq(ns::ASCIIString) = ns

immutable NucTideSeq{T<:NucleotideLike} <: NucSeq
    seq::Vector{T}
end


nucseq(seq::Vector{Char}) = NucTideSeq(seq)

sequence{NS<:NucTideSeq}(nseq::NS) = nseq.seq
seqtype{T<:NucleotideLike}(::NucTideSeq{T}) = Vector{T}


#nucleotides are put in a sequence on a strand
#in a particular sense
abstract Strand <: NucSeq

# sense of a Strand tells us if the increasing index
# of the stored sequence is sense or antisense.
# the stored sequence will always be iterated
#in its sense sense
isForward{NS<:Strand}(ns::NS) = sense(ns) == Forward
isReverse{NS<:Strand}(ns::NS) = sense(ns) == Reverse

function Base.start{NS<:Strand}(nseq::NS)
    seq = sequence(nseq)
    isForward(nseq) ?
    Base.start(seq) :
    Base.endof(seq)
end

function Base.next{NS<:Strand}(nseq::NS, pos)
    seq = sequence(nseq)
    isForward(nseq) ?
    Base.next(seq, pos) :
    (Base.getindex(seq, pos), Base.prevind(seq, pos))
end

function Base.done{NS<:Strand}(nseq::NS, pos)
    seq = sequence(nseq)
    isForward(nseq) ?
    Base.done(seq, pos) :
    pos < Base.start(seq)
end

Base.getindex{NS<:Strand}(ns::NS, posr::UnitRange{Int64}) =
    sequence(ns)[posr]

Base.reverse{NS<:Strand}(nseq::NS) =
    NS(sequence(nseq), anti(sense(nseq)))

complement{NS<:Strand}(nseq::NS) =
    NS(map(complement, sequence(nseq)),
       anti(sense(nseq)))

isForward{NS<:Strand}(nseq::NS) = isForward(sense(nseq))

Base.reverseind{NS<:Strand}(nseq::NS, i::Integer) =
    Base.reverseind(sequence(nseq), i)

function Base.getindex{NS<:Strand}(nseq::NS, i::Integer)
    j = isForward(nseq) ? i : reverseind(nseq, i)
    getindex(sequence(nseq), j)
end

function getrevindex{NS<:Strand}(nseq::NS, i::Integer)
    j = isForward(nseq) ? reverseind(nseq, i) : i
    getindex(sequence(nseq), j)
end

# a Strand should always return a sequence
# in the sense sense
function raw{NS<:Strand}(nseq::NS)
    seq = sequence(nseq)
    isForward(nseq) ? seq : Base.reverse(seq)
end

function asString{NS<:Strand}(nseq::NS)
    seq = sequence(nseq)
    isForward(nseq) ?
    asString(seq) :
    Base.reverse(asString(seq))
end

immutable EmptyStrand <: Strand
    sense::Sense
end

sequence(ns::EmptyStrand) = []
sense(ns::EmptyStrand) = ns.sense
seqtype(ns::EmptyStrand) = Void
seqtype(::Type{EmptyStrand}) = Void

asString(ns::EmptyStrand) = ""

convert(ASCIIString, seq::EmptyStrand) = ""

complement(ns::EmptyStrand) = EmptyStrand(anti(sense(ns)))

immutable NucCharStrand <: Strand
    sequence::ASCIIString
    sense::Sense
end

sequence(ns::NucCharStrand) = ns.sequence
sense(ns::NucCharStrand) = ns.sense
seqtype(ns::NucCharStrand) = ASCIIString
seqtype(::Type{NucCharStrand}) = ASCIIString

asString(seq::ASCIIString) = seq

convert(ASCIIString, ns::NucCharStrand) = sequence(ns)

immutable NucTideStrand{T <: NucleotideLike} <: Strand
    sequence::Vector{T}
    sense::Sense
end

sequence(ns::NucTideStrand) = ns.sequence
sense(ns::NucTideStrand) = ns.sense
seqtype{T}(ns::NucTideStrand{T}) = Vector{T}
seqtype{T<:Nucleotide}(::Type{NucTideStrand{T}}) = Vector{T}

NucTideStrand{T<:Nucleotide}(seq::Vector{T}) =
    NucTideStrand(seq, Forward())

convert{C<:NucCharStrand}(ASCIIString, nseq::NucTideStrand{C}) =
    convert(ASCIIString, sequence(nseq))

nuctideSeq(seq::AbstractString) =
    NucTideStrand(map(nucleochar, collect(seq)))

nuctwoboolSeq(seq::ASCIIString) =
    NucTideStrand(map(nucleotwobool, collect(seq)), Forward())

NucSeqLike = Union{NucSeq, ASCIIString}



function hammingDistance{NS <: NucSeqLike}(xs::NS, ys::NS)
    length(xs) != length(ys) ?
    Inf32 :
    sum(xy -> xy[1] != xy[2], zip(xs, ys))
end

