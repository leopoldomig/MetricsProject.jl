using MetricsProject
using Test
using Random

@testset "MetricsProject.jl" begin
    # Tests for Berry1994
    Random.seed!(1)
    out = Berry1994(50, 10).comparison
    @test out.m_θhatⁱᵛ[1]  ≈ 3.906248115012579
    @test out.m_θhatⁱᵛ[2]  ≈ 1.695701084241484
    @test out.m_θhatⁱᵛ[3]  ≈ -0.8016480548335396

    # Tests for sdc (letting the default seed parameter = 0)
    out = sdc(500, 5).comparison
    @test out.θhat_t1ev[5] ≈ 1.566535306145989
    @test out.θhat_normal[5] ≈ 1.331458650920062

end
