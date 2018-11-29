using Base.Iterators
using Test

include("../src/Gillespie.jl")
include("../src/Model_CRS.jl")
include("../src/Droplets.jl")
include("../src/Data.jl")

@testset "CRS helpers" begin
   
    @test Set(generate_all_binary_seqs(2)) == Set(["A", "B", "AA","AB","BA", "BB"])
    @test Set(generate_all_binary_seqs(3)) == Set(["A", "B", "AA","AB","BA", "BB", "AAA", "AAB","ABA","BAA","BBA","BAB","ABB","BBB"])

    @test all_splits("12345") == (["1","12", "123","1234"], ["2345","345", "45", "5"])
    
    @test countall("10", "1") == 1
    @test countall("10", "0") == 1
    @test countall("111", "1") == 3
    @test countall("101", "1") == 2
    @test countall("BA", "1") == 0
    @test countall("AB", "B") == 1

    @test reduced_mass(1,1) == 0.5
    @test isapprox(reduced_mass(10000,1),1, atol = 1E-4)
    
    @test calculate_reduced_mass("AAA", "AB", 1, 2) == 1.5
    @test isapprox(calculate_reduced_mass("AAA", "BB", 1, 2), 1.71428571429)
    
    @test bimolecular_coef(1,1,1,1,1) == 1.0
    @test isapprox(bimolecular_coef(0.001, 300, 8.5, 1000, 10), 0.00001596871942267131199907024517698) # k = 0.001, T = 300, R = 8.5, V = 1000, m = 10
end


@testset "Reaction Funcs" begin
   
    @test PickReactionID([1.0, 0.0, 0.0]) == 1
    @test PickReactionID([0.0, 0.0, 1.0]) == 3
    @test PickReactionID([0.0, 0.0001, 0.0]) == 2
    ### Include som fancier tests here to make sure relative proportions are correct
    
end
