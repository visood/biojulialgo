include("../src/cluster.jl")
include("../test/test.jl")
#test cluster, and clumps

#vcfile = open("../data/Vibrio_cholerae.txt")
#vcseq = readline(vcfile)[1:end-1]
#vcnseq = NucCharStrand(vcseq, Forward())

to_test_that("Kmers can be equal") do
    km1 = Kmer("AC", 1)
    km2 = Kmer("AC", 1)
    km3 = Kmer("GT", 2)
    km4 = Kmer("GT", 2)
    km5 = Kmer("AC", 2)
    km6 = Kmer("GT", 1)
    @test (km1 == km1 &&
           km1 == km2 &&
           km1 != km3 &&
           km1 != km4 &&
           km1 != km5 &&
           km1 != km6)
end

to_test_that("NucCharStrand returns a valid result") do
    dcl1 = Dict(NucCharStrand("ACAT", Forward()) => [1],
                NucCharStrand("AATA", Forward()) => [5],
                NucCharStrand("TCGA", Forward()) => [10, 14, 20],
                NucCharStrand("AAGT", Forward()) => [25])
    ns = NucCharStrand("ACAT", Forward())
    @test(sequence(ns) == "ACAT" &&
          length(dcl1) == 4)
end

to_test_that("dropWhile drops till the first element that does not satisfy the drop condition" ) do
    xs = Integer[1,2,3,4,5,6,8,17, 23, 31, 40]
    xds = dropWhile(x -> x < 17, xs)
    ys = Integer[1,2,3,4,5,6,8, 17, 9, 0, 1, 2]
    yds = dropWhile(y -> y < 17, ys)
    @test(xds == Integer[17, 23, 31, 40] &&
          yds == Integer[17, 9, 0, 1, 2])
end

xs = Integer[1,2,3,4,5,7,8,9,11,12,17,23]
L = 5
t = 2
cxs = clumps(L, t, xs)

to_test_that("positions in an integers clump should becloser than distance $(L+1)")  do
    @test(all(map(xs -> maximum(xs) - minimum(xs) <= L, cxs)))
end

to_test_that("an integers clump should be at least size $t") do
    @test all(map(x -> length(x) >= t, cxs))
end

cs = map(x -> cluster("ACTG", x), xs)
ccs = clumps(L, t, cs)

to_test_that("clusters can be clumped") do
    @test length(ccs) > 0
end

to_test_that("positions in a clusters clump should be closer than distance $(L+1)" ) do
    @test all(map(x -> rightmost(x) - leftmost(x) <= L, ccs))
end

kms = map(x -> Kmer("ACTG", x), xs)
ckms = clumps(L, t, kms)
to_test_that("positions in a kmers clump should becloser than distance $(L+1)") do
    @test all(map(x -> rightmost(x) - leftmost(x) <= L, ckms))
end

clms = map(x -> clumer("ACTG", x), xs)
cclms = clumps(L, t, clms)
to_test_that("positions in a clumers clump should be closer than distance $(L+1)") do
    @test all(map(x -> rightmost(x) - leftmost(x) <= L, cclms))
end

xs2 = Integer[6, 10, 13, 14, 15, 16, 18, 21]
kms2 = append!(kms, map(x -> Kmer("TGCA", x), xs2))
ckms2 = clumps(L, t, kms2)
to_test_that("positions in a clump formed from two sets of kmers  with different sequences should be closer than distance $(L+1)") do
    @test all(map(x -> rightmost(x) - leftmost(x) <= L, ckms2))
end


clms = map(x -> clumer("ACTG", x), xs)
cclms = clumps(L, t, clms)
to_test_that("positions in a clumers clump should be closer than distance $(L+1)") do
    @test all(map(x -> rightmost(x) - leftmost(x) <= L, cclms))
end

clms2 = append!(clms, map(x -> clumer("TCGA", x), xs2))
cclms2 = clumps(L, t, clms2)
to_test_that("positions in  clumps formed from two sets of clumers with different sequences should be closer than distance $(L+1)") do
    @test all(map(x -> rightmost(x) - leftmost(x) <= L, cclms2))
end
