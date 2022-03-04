using FirstPassageTools
using Test

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
