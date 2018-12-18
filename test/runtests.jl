using MetricsProject
using Test
using Random

@testset "MetricsProject.jl" begin
    # Write your own tests here.
    Random.seed!(1)
    out = Berry1994(50, 10).comparison
    @test out.m_θhatⁱᵛ[1]  ≈ 3.906248115012579
    @test out.m_θhatⁱᵛ[2]  ≈ 1.695701084241484
    @test out.m_θhatⁱᵛ[3]  ≈ -0.8016480548335396
end
