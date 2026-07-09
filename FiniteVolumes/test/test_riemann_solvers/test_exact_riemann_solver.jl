const RTOL_TORO = 1e-4
const ATOL_ZERO = 1e-5

@testset "Toro exact Riemann star states" begin

    @testset "Test 1: Sod" begin
        WL = SVector(1.0, 0.0, 1.0)
        WR = SVector(0.125, 0.0, 0.1)
        ps, us, ρs_l, ρs_r = get_star_values(WL, WR, 1.4)
        @test ps   ≈ 0.30313 rtol=RTOL_TORO
        @test us   ≈ 0.92745 rtol=RTOL_TORO
        @test ρs_l ≈ 0.42632 rtol=RTOL_TORO
        @test ρs_r ≈ 0.26557 rtol=RTOL_TORO
    end

    @testset "Test 2: 123 problem" begin
        WL = SVector(1.0, -2.0, 0.4)
        WR = SVector(1.0,  2.0, 0.4)
        ps, us, ρs_l, ρs_r = get_star_values(WL, WR, 1.4)
        @test ps   ≈ 0.00189 rtol=1e-2
        @test us   ≈ 0.0     atol=ATOL_ZERO
        @test ρs_l ≈ 0.02185 rtol=1e-3
        @test ρs_r ≈ 0.02185 rtol=1e-3
    end

    @testset "Test 3: right half of blast wave" begin
        WL = SVector(1.0, 0.0, 1000.0)
        WR = SVector(1.0, 0.0, 0.01)
        ps, us, ρs_l, ρs_r = get_star_values(WL, WR, 1.4)
        @test ps   ≈ 460.894 rtol=RTOL_TORO
        @test us   ≈ 19.5975 rtol=RTOL_TORO
        @test ρs_l ≈ 0.57506 rtol=RTOL_TORO
        @test ρs_r ≈ 5.99924 rtol=RTOL_TORO
    end

    @testset "Test 4: left half of blast wave" begin
        WL = SVector(1.0, 0.0, 0.01)
        WR = SVector(1.0, 0.0, 100.0)
        ps, us, ρs_l, ρs_r = get_star_values(WL, WR, 1.4)
        @test ps   ≈ 46.0950  rtol=RTOL_TORO
        @test us   ≈ -6.19633 rtol=RTOL_TORO
        @test ρs_l ≈ 5.99242  rtol=RTOL_TORO
        @test ρs_r ≈ 0.57511  rtol=RTOL_TORO
    end

    @testset "Test 5: collision of two strong shocks" begin
        WL = SVector(5.99924,  19.5975, 460.894)
        WR = SVector(5.99242, -6.19633, 46.0950)
        ps, us, ρs_l, ρs_r = get_star_values(WL, WR, 1.4)
        @test ps   ≈ 1691.64 rtol=RTOL_TORO
        @test us   ≈ 8.68975 rtol=RTOL_TORO
        @test ρs_l ≈ 14.2823 rtol=RTOL_TORO
        @test ρs_r ≈ 31.0426 rtol=RTOL_TORO
    end
end
