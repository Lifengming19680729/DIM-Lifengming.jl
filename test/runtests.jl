using Test
using DynamicIntervals

@testset "DynamicIntervals.jl Tests" begin

@testset "Core Construction" begin
    a = DInt(2.0, -0.1, 0.1)
    @test a.b == 2.0
    @test a.δm == -0.1
    @test a.δp == 0.1
    @test lower_bound(a) ≈ 1.9
    @test upper_bound(a) ≈ 2.1
end

@testset "Addition" begin
    a = DInt(1.0, -0.2, 0.3)
    b = DInt(2.0, -0.1, 0.4)
    c = a + b
    @test c.b == 3.0
    @test c.δm ≈ -0.3
    @test c.δp ≈ 0.7
end

@testset "Multiplication" begin
    a = DInt(2.0, -0.1, 0.1)
    b = DInt(3.0, -0.2, 0.2)
    c = a * b
    @test c.b ≈ 6.0
    @test c.δm ≈ (2.0 * -0.2 + 3.0 * -0.1)
    @test c.δp ≈ (2.0 * 0.2 + 3.0 * 0.1)
end

println("All tests passed!")
end
