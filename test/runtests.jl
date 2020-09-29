using Test, Genesys

# Parameters
const Δh, nh, ny, ns = 1, 8760, 20, 1

# Initialization
liion = Liion() ; Genesys.preallocate!(liion, nh, ny, ns)
h2tank = H2Tank() ; Genesys.preallocate!(h2tank, nh, ny, ns)
tes = ThermalSto() ; Genesys.preallocate!(tes, nh, ny, ns)
fc = FuelCell() ; Genesys.preallocate!(fc, nh, ny, ns)
elyz = Electrolyzer() ; Genesys.preallocate!(elyz, nh, ny, ns)
heater = Heater() ; Genesys.preallocate!(heater, nh, ny, ns)

### Operation
# Liion test set
@testset "Battery compute_operation_dynamics tests" begin
    @testset "Battery charge/discharge" begin
        x_liion = (Erated=100, soc=0.5, soh=1)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 10., Δh) == (0.5 * (1. - 0.) - 10. / 0.8 * 1. / 100., 1. - 10. / 2. / 2500. / 0.6 / 100., 10.)
        @test Genesys.compute_operation_dynamics(liion, x_liion, -10., Δh) == (0.5 * (1. - 0.) + 0.8 * 10. * 1. / 100., 1. - 10. / 2. / 2500. / 0.6 / 100., -10.)
    end
    @testset "Battery SoC bounds" begin
        x_liion = (Erated=100, soc=0.5, soh=1)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 100., Δh) == (0.5 * (1. - 0.), 1.0, 0.0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, -100., Δh) == (0.5 * (1. - 0.), 1.0, 0.0)
    end
    @testset "Battery power bounds" begin
        x_liion = (Erated=100, soc=0.5, soh=1)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 1000000., Δh) == (0.5 * (1. - 0.), 1.0, 0.0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, -1000000., Δh) == (0.5 * (1. - 0.), 1.0, 0.0)
    end
    @testset "Battery soh" begin
        x_liion = (Erated=100, soc=0.5, soh=0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 100., Δh) == (0.5 * (1. - 0.), 0., 0.0)
        x_liion = (Erated=100, soc=0.2, soh=0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 100., Δh) == (0.2, 0., 0.0)
    end
end

# H2 tank test set
@testset "H2 tank  compute_operation_dynamics tests" begin
    @testset "H2 tank charge/discharge" begin
        x_h2 = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, 10., Δh) == (0.5 * (1. - 0.) - 10 / 1. * 1 / 100, 10)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, -10., Δh) == (0.5 * (1. - 0.) + 1. * 10 * 1 / 100, -10)
    end
    @testset "H2 tank  SoC bounds" begin
        x_h2 = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, 100., Δh) == (0.5 * (1. - 0.), 0.0)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, -100., Δh) == (0.5 * (1. - 0.), 0.0)
    end
    @testset "H2 tank  power bounds" begin
        x_h2 = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, 1000000., Δh) == (0.5 * (1. - 0.), 0.0)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, -1000000., Δh) == (0.5 * (1. - 0.), 0.0)
    end
end

# TES test set
@testset "TES  compute_operation_dynamics tests" begin
    @testset "TES charge/discharge" begin
        x_tes = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(tes, x_tes, 10., Δh) == (0.5 * (1. - 0.01) - 10. / 0.9 * 1. / 100., 10.)
        @test Genesys.compute_operation_dynamics(tes, x_tes, -10., Δh) == (0.5 * (1. - 0.01) + 0.9 * 10. * 1. / 100., -10.)
    end
    @testset "TES  SoC bounds" begin
        x_tes = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(tes, x_tes, 100., Δh) == (0.5 * (1. - 0.01), 0.0)
        @test Genesys.compute_operation_dynamics(tes, x_tes, -100., Δh) == (0.5 * (1. - 0.01), 0.0)
    end
    @testset "TES  power bounds" begin
        x_tes = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(tes, x_tes, 1000000., Δh) == (0.5 * (1. - 0.01), 0.0)
        @test Genesys.compute_operation_dynamics(tes, x_tes, -1000000., Δh) == (0.5 * (1. - 0.01), 0.0)
    end
end
# Heater test set
@testset "Heater compute_operation_dynamics tests" begin
    @testset "Heater power" begin
        x_heater = (powerMax=100,)
        @test Genesys.compute_operation_dynamics(heater, x_heater, -10., Δh) == (-10., 10.)
        @test Genesys.compute_operation_dynamics(heater, x_heater, 10., Δh) == (0., 0.)
    end
    @testset "Heater bounds" begin
        x_heater = (powerMax=100,)
        @test Genesys.compute_operation_dynamics(heater, x_heater, -1000., Δh) == (0., 0.)
    end
end
#  Fuel Cell test set
@testset "Fuel Cell compute_operation_dynamics tests" begin
    @testset "fc power" begin
        x_fc = (powerMax=100., soh=1.)
        @test Genesys.compute_operation_dynamics(fc, x_fc, 50., Δh) == (50., 50. * 0.4 / 0.4, -50. / 0.4, 1. - 1. / 10000.)
        @test Genesys.compute_operation_dynamics(fc, x_fc, -50., Δh) == (0, 0, 0, 1)
    end
    @testset "fc power bounds" begin
        x_fc = (powerMax=100., soh=1.)
        @test Genesys.compute_operation_dynamics(fc, x_fc, 200., Δh)  == (0., 0., 0., 1.)
        @test Genesys.compute_operation_dynamics(fc, x_fc, 2., Δh)  == (0., 0., 0., 1.)
    end
    @testset "fc soh" begin
        x_fc = (powerMax=100., soh=0.)
        @test Genesys.compute_operation_dynamics(fc, x_fc, 50., Δh)  == (0., 0., 0., 0.)
    end
end
#  Electrolyzer test set
@testset "Electrolyzer compute_operation_dynamics tests" begin
    @testset "elyz power" begin
        x_elyz = (powerMax=100., soh=1.)
        @test Genesys.compute_operation_dynamics(elyz, x_elyz, -50., Δh) == (-50., 50. * 0.3, 50. * 0.5, 1. - 1. / 26000.)
        @test Genesys.compute_operation_dynamics(elyz, x_elyz, 50., Δh) == (0, 0, 0, 1)
    end
    @testset "elyz power bounds" begin
        x_elyz = (powerMax=100., soh=1.)
        @test Genesys.compute_operation_dynamics(elyz, x_elyz, -200., Δh)  == (0., 0., 0., 1.)
        @test Genesys.compute_operation_dynamics(elyz, x_elyz, -2., Δh)  == (0., 0., 0., 1.)
    end
    @testset "elyz soh" begin
        x_elyz = (powerMax=100., soh=0.)
        @test Genesys.compute_operation_dynamics(elyz, x_elyz, -50., Δh)  == (0., 0., 0., 0.)
    end
end

### Investment
# Liion test set
@testset "Battery compute_investment_dynamics tests" begin
    x_liion=(Erated=200, soc=0.5, soh=1)
    @test Genesys.compute_investment_dynamics(liion, x_liion, 300) == (300, 0.5, 1) # liion.soc[1,1] = 0 in that case...
    @test Genesys.compute_investment_dynamics(liion, x_liion, 0) == (200, 0.5, 1)
end
# H2 tank test set
@testset "H2tank compute_investment_dynamics tests" begin
    x_tank=(Erated=200, soc=0.5)
    @test Genesys.compute_investment_dynamics(h2tank, x_tank, 300) == (300, 0.5)
    @test Genesys.compute_investment_dynamics(h2tank, x_tank, 0) == (200, 0.5)
end
# TES test set
@testset "TES compute_investment_dynamics tests" begin
    x_tes=(Erated=200, soc=0.5)
    @test Genesys.compute_investment_dynamics(tes, x_tes, 300) == (300, 0.5)
    @test Genesys.compute_investment_dynamics(tes, x_tes, 0) == (200, 0.5)
end
# Elyz Test set
@testset "Elyz compute_investment_dynamics tests" begin
    x_elyz=(powerMax=200, soh=1)
    @test Genesys.compute_investment_dynamics(elyz, x_elyz, 300) == (300, 1)
    @test Genesys.compute_investment_dynamics(elyz, x_elyz, 0) == (200, 1)
end
# FC Test set
@testset "Fc compute_investment_dynamics tests" begin
    x_fc=(powerMax=200, soh=1)
    @test Genesys.compute_investment_dynamics(fc, x_fc, 300) == (300, 1)
    @test Genesys.compute_investment_dynamics(fc, x_fc, 0) == (200, 1)
end
