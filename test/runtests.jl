using FirstPassageTools
using Test
using CSV
using DataFrames

@testset "Testing setup" begin
    Ttest, Atest = setup("goodtest.csv")
    @test isapprox(Ttest, [-3.0 0; 3 -3])
    @test isapprox(Atest, [0 3.0])
    @test_throws AssertionError setup("badtest.csv")
end

@testset "Tests of rescale!()" begin
    # Test with T and A matrices
    Torig = [-1 1.0; 1 -2]
    Aorig = [0 1.0]
    Tcorr = [-2.0 2; 2 -4]
    Acorr = [0 2.0]
    Ttest = rescale!(Torig, 2.0)
    Atest = rescale!(Aorig, 2.0)
    @test isapprox(Tcorr, Ttest)
    @test isapprox(Acorr, Atest)

    # Test with full W-matrix
    Wtest = cat(Torig, Aorig, dims=1)
    Wtest = hcat(Wtest, zeros(size(Wtest, 1)))
    @test isapprox(3 * Wtest, rescale!(Wtest, 3))
end
