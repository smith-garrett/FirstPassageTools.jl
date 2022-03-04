using FirstPassageTools
using Test

@testset "Tests of rescale!()" begin
    Torig = [-1 1.0; 1 -2]
    Aorig = [0 1.0]
    Tcorr = [-2.0 2; 2 -4]
    Acorr = [0 2.0]
    Ttest, Atest = rescale!(Torig, 2.0)
    @test isapprox(Tcorr, Ttest)
    @test isapprox(Acorr, Atest)
end
