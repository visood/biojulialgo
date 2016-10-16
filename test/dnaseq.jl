include("../src/dnaseq.jl")
include("../test/test.jl")


to_test_that("Complement of a Nucleochar is defined correctly") do
    a = nucleochar('A')
    c = nucleochar('C')
    g = nucleochar('G')
    t = nucleochar('T')
    x = nucleochar('P')
    @test (complement(a) == t &&
           complement(c) == g &&
           complement(g) == c &&
           complement(t) == a &&
           complement(x) == x)
end

to_test_that("Complement of a NucleoTwoBool is defined correctly") do
    a = nucleotwobool('A')
    c = nucleotwobool('C')
    g = nucleotwobool('G')
    t = nucleotwobool('T')
    @test (complement(a) == t &&
           complement(c) == g &&
           complement(g) == c &&
           complement(t) == a)
end

to_test_that("ASCIIString can be treated as NucSeq") do
    ns = nucseq("ACGTAAAGT_AG")
    @test (! isEmpty(ns) &&
           sequence(ns) == ns &&
           seqtype(ns) == ASCIIString &&
           length(ns) == 12 &&
           start(ns) == 1 &&
           next(ns, 1) == ('A', 2) &&
           ! done(ns, 1) &&
           ! done(ns, 12) &&
           done(ns, 13) &&
           complement(ns) == "TGCATTTCA_TC")
end

to_test_that("A vector of Chars can be packaged into a NucSeq") do
    ns = nucseq(['A', 'C', 'T', 'G', 'U', '=', 'X', 'N'])
    @test (! isEmpty(ns) &&
           sequence(ns) == ['A', 'C', 'T', 'G', 'U', '=', 'X', 'N'] &&
           seqtype(ns) == Vector{Char} &&
           length(ns) == 8 &&
           start(ns) == 1 &&
           next(ns, 1) == ('A', 2) &&
           ! done(ns, 1) &&
           ! done(ns, 8) &&
           done(ns, 9) &&
           complement(ns) == nucseq(['T', 'G', 'A', 'C', 'C', '=', '_', 'N']))
end

to_test_that("EmptyStrand is empty") do
    e = EmptyStrand(Forward())
    @test (length(e) == 0 &&
           sequence(e) == [] &&
           sense(e) == Forward() &&
           seqtype(e) == Void &&
           asString(e) == "" &&
           start(e) == 1 &&
           done(e, 1) &&
           complement(e) == e )
end

to_test_that("An ASCIIString can be packaged into a NucCharStrand") do
    nsf = NucCharStrand("ACGT", Forward())
    nsb = NucCharStrand("TGCA", Reverse())
    @test (length(nsf) == 4 &&
           length(nsb) == 4 &&
           sequence(nsf) == "ACGT" &&
           sequence(nsb) == "TGCA" &&
           sense(nsf) == Forward() &&
           sense(nsb) == Reverse() &&
           seqtype(nsf) == ASCIIString &&
           seqtype(nsb) == ASCIIString &&
           start(nsf) == 1 &&
           start(nsb) == 4 &&
           next(nsf, 1) == ('A', 2) &&
           next(nsb, 4) == ('A', 3) &&
           ! done(nsf, 1) &&
           ! done(nsb, 4) &&
           done(nsf, 5) &&
           done(nsb, 0) &&
           complement(nsf) == nsb &&
           complement(nsb) == nsf )
end

to_test_that("A vector of Chars can be packaged into a NucTideStrand") do
    nsf = NucTideStrand(['A', 'C', 'G', 'T'], Forward())
    nsb = NucTideStrand(['T', 'G', 'C', 'A'], Reverse())
    @test (length(nsf) == 4 &&
           length(nsb) == 4 &&
           sequence(nsf) == ['A', 'C', 'G', 'T'] &&
           sequence(nsb) == ['T', 'G', 'C', 'A'] &&
           sense(nsf) == Forward() &&
           sense(nsb) == Reverse() &&
           seqtype(nsf) == Vector{Char} &&
           seqtype(nsb) == Vector{Char} &&
           start(nsf) == 1 &&
           start(nsb) == 4 &&
           next(nsf, 1) == ('A', 2) &&
           next(nsb, 4) == ('A', 3) &&
           ! done(nsf, 1) &&
           ! done(nsb, 4) &&
           done(nsf, 5) &&
           done(nsb, 0) &&
           complement(nsf) == nsb &&
           complement(nsb) == nsf )
end


to_test_that("Hamming distance is computed correctly") do
    xs = "AAAAACGTAAAAACGTCCCCACGT"
    ys = "AAAAACGTAAAAACGTCCCCTGCA"
    @test (hammingDistance(xs, ys) == 4 &&
           hammingDistance(xs, ys[2:end]) == Inf32)
end
