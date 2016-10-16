include("../test/test.jl")
include("../src/oric.jl")
include("../src/numseq.jl")

to_test_that("A KmerItr over a Strand respects its sense") do
    fitr = kmerItr(NucCharStrand("ACTG", Forward()), 2)
    ritr = kmerItr(NucCharStrand("ACTG", Reverse()), 2)
    @test (kmerLength(fitr) == 2 &&
           kmerLength(ritr) == 2 &&
           sequence(fitr) == "ACTG" &&
           sequence(ritr) == "GTCA" &&
           seqtype(fitr) == ASCIIString)
end

to_test_that("A KmerItr over a forward strand iterates right") do
    fitr = kmerItr(NucCharStrand("ACGT", Forward()), 2)
    @test (Base.start(fitr) == 1 &&
           !Base.done(fitr, 1) &&
           Base.next(fitr, 1) == (Kmer("AC", 1), 2) &&
           Base.next(fitr, 2) == (Kmer("CG", 2), 3) &&
           Base.next(fitr, 3) == (Kmer("GT", 3), 4) &&
           Base.done(fitr, 4) )
end

to_test_that("A KmerItr over a reverse strand iterates right") do
    ritr = kmerItr(NucCharStrand("ACTG", Reverse()), 2)
    @test (Base.start(ritr) == 1 &&
           !Base.done(ritr, 1) &&
           Base.next(ritr, 1) == (Kmer("GT", 1), 2) &&
           Base.next(ritr, 2) == (Kmer("TC", 2), 3) &&
           Base.next(ritr, 3) == (Kmer("CA", 3), 4) &&
           Base.done(ritr, 4) )
end

seq1 = "ACGTAAAAACGTAAAAACGTCCCCACGT"

to_test_that("KmerItr counts the kmers correctly") do
    seq = "ACGTAAAAACGTAAAAACGTCCCCACGT"
    k = 4
    kitr = KmerItr(seq, k)
    kcounts = kmerCounts(kitr)
    @test(sum(map(x->x[2], kcounts)) == length(sequence(kitr) - k + 1) &&
          kcounts["ACGT"] == 4 &&
          kcounts["AAAA"] == 4 &&
          kcounts["CCCC"] == 1)
end

to_test_that("countedAppearances gives the Kmers that appear a given number of times") do
    kitr = KmerItr(seq1, 4)
    kcounts = kmerCounts(kitr)
    caps = countedAppearances(kcounts)
    @test(in("ACGT", caps[4]) &&
          in("AAAA", caps[4]) &&
          in("CCCC", caps[1]) )
end

to_test_that("topFrequentKmers gives a kmers ordered by their frequency") do
    seq = "ACGTAAAAACGTAAAAACGTCCCCACGT"
    k = 4
    tkmers = topFrequentKmers(seq, k)
    tkmers3 = topFrequentKmers(seq, k, 3)
    counts = map(x -> x[2], tkmers)
    @test (tkmers[1] == Pair("AAAA", 4) &&
           tkmers[2] == Pair("ACGT", 4) &&
           sort(counts, rev = true) == counts &&
           length(tkmers3) == 3 &&
           tkmers3[1] == Pair("AAAA", 4) &&
           tkmers3[2] == Pair("ACGT", 4) &&
           tkmers3[3] == Pair("CGTA", 2) )
end

to_test_that("kmerlocs can locate kmers in a sequence") do
    seq = "ACGTAAAAACGTAAAAACGTCCCCACGT"
    k = 4
    kitr = KmerItr(seq, k)
    klocs = kmerlocs(kitr)
    @test(sum(map(x -> length(x[2]), klocs)) == length(sequence(kitr)) - k + 1 &&
          klocs["ACGT"] == Integer[1, 9, 17, 25])
end

to_test_that("a pattern can be located in a sequence") do
    seq = "ACGTAAAAACGTAAAAACGTCCCCACGT"
    k = 4
    klocs = kmerlocs("ACGT", seq)
    @test( klocs == [1, 9, 17, 25])
end


As = repeat("A", 1000000)
ACTGs = repeat("ACTG", 25000)
Cs = repeat("C", 1000000)
longseq = "$As$ACTGs$Cs"


lskitr = KmerItr(longseq, 4)
lskcounts = kmerCounts(lskitr)
ckmers = countedAppearances(lskcounts)
topkmers = topFrequentKmers(longseq, 4)
lskmocs = kmerlocs(lskitr)

to_test_that("Kmer counting works as expected") do
    @test(lskcounts["AAAA"] == 999998 &&
          lskcounts["ACTG"] == 25000 &&
          lskcounts["CTGA"] == 24999 &&
          lskcounts["TGAC"] == 24999 &&
          lskcounts["GACT"] == 24999 &&
          lskcounts["CTGC"] == 1 &&
          lskcounts["TGCC"] == 1 &&
          lskcounts["GCCC"] == 1 &&
          lskcounts["CCCC"] == 999997)
end


pcACTG = runningPatternCount("ACTG", longseq)

to_test_that("Pattern counts along the sequence are computed correctly") do
    @test(monotonicity(pcACTG) == 1)
end

to_test_that("aeproximate kmers can be located") do
    aseq = "ACTCGCAAGGRTC"
    kitr = KmerItr(aseq, 2)
    kms = map(element, kitr)
    aklocs = appxKmerLocs(2, aseq, 1)
    @test( all(w -> aklocs[w] == appxKmerLocs(w, aseq, 1), kms) )
end


vcfile = open("../data/Vibrio_cholerae.txt")
vcseq = readline(vcfile)[1:end-1]
vcnseq = NucCharStrand(vcseq, Forward())
vckitr = KmerItr(vcnseq, 9)
vckcounts = kmerCounts(vckitr)
ckmers = countedAppearances(vckcounts)
topkmers = topFrequentKmers(vcseq, 9)
vckmocs = kmerlocs(vckitr);

to_test_that("the opened file provided a string") do
    @test(typeof(vcseq) == ASCIIString)
end

