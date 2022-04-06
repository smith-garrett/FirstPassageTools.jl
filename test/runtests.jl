using FirstPassageTools
using Test
using CSV
using DataFrames
using Distributions

T = [-1 1.0; 1 -2]
A = [0 1.0]
p0 = [1.0, 0]

@testset "Testing setup" begin
    Ttest, Atest = setup("goodtest.csv")
    @test isapprox(Ttest, [-3.0 0; 3 -3])
    @test isapprox(Atest, [0 3.0])
    @test_throws AssertionError setup("badtest.csv")
end

@testset "Tests of rescale!()" begin
    # Test with T and A matrices
    Tcopy = deepcopy(T)
    Acopy = deepcopy(A)
    Tcorr = [-2.0 2; 2 -4]
    Acorr = [0 2.0]
    Ttest = rescale!(Tcopy, 2.0)
    Atest = rescale!(Acopy, 2.0)
    @test isapprox(Tcorr, Ttest)
    @test isapprox(Acorr, Atest)

    # Test with full W-matrix
    Wtest = cat(Tcopy, Acopy, dims=1)
    Wtest = hcat(Wtest, zeros(size(Wtest, 1)))
    @test isapprox(3 * Wtest, rescale!(Wtest, 3))
end

@testset "Tests for building an fpdistribution" begin
    # Should pass:
    @test_nowarn fpdistribution(T, A, p0)
    
    # Should fail:
    Tbad = [-1 1.0; 1 -1]
    Abad = [0 1.1]
    p0bad = [1.0, 0.1]
    @test_throws AssertionError fpdistribution(Tbad, A, p0)
    @test_throws AssertionError fpdistribution(T, Abad, p0)
    @test_throws AssertionError fpdistribution(T, A, p0bad)
    @test_throws AssertionError fpdistribution(T, A, [1.0])
end

@testset "Erlang distribution tests" begin
    # goodtest system should be equivalent to an Erlang distribution with 2 transient 
    # states and a scale parameter equal to 1/(total # states)
    erl = Erlang(2, 1/3)
    Terl, Aerl = setup("goodtest.csv")
    p0 = [1.0, 0]
    fperl = fpdistribution(Terl, Aerl, p0)
    times = Array(0.0:5:20)

    @test all(isapprox.(pdf.(erl, times), pdf.(fperl, times)))
    @test all(isapprox.(logpdf.(erl, times), logpdf.(fperl, times)))
    @test all(isapprox.(cdf.(erl, times), cdf.(fperl, times)))
    # Using a looser tolerance for this test because the quantile fn. for fpdistribution
    # relies on numerical root finding
    qtiles = Array(0:0.1:0.9)
    @test all(isapprox.(quantile.(erl, qtiles), quantile.(fperl, qtiles), atol=1e-6))
    @test isapprox(mean(erl), mean(fperl))
    @test isapprox(var(erl), var(fperl))
end

@testset "Multiple absorbing states" begin
    Tm, Am = setup("two_absorbing.csv")
    two_abs = fpdistribution(Tm, Am, [1.0, 0, 0])
    @test all(isapprox.([0.5, 0.5], splittingprobabilities(two_abs)))
    @test_nowarn pdf.(two_abs, [0.0, 5], 1)
    @test_nowarn pdf.(two_abs, [0.0, 5], [1, 2])
end
