using Test, Genesys

# Parameters
const Δh, H = 1, 24*365
const nh = length(Δh:Δh:H)
const Δy, Y = 1, 20
const ny = length(Δy:Δy:Y)
const ns = 10


outputGUI = Dict(
"ld" => [],
"pv" => [],
"liion" => (Erated_ini = 0.,α_p_ch = 1.5,α_p_dch = 1.5,η_ch = 0.8,η_dch = 0.8,
η_self = 0.1,α_soc_min = 0.2,α_soc_max = 0.8,α_soc_ini = 0.8,soh_ini=1.,lifetime=10,nCycle = 2500,
dod = 0.6,C_liion = zeros(ny,ns)),
"elyz" => (powerMax_ini = 0.,α_p = 5/100,η_E_H2 = 0.5,η_E_H = 0.3,soh_ini = 1.,lifetime=10,nHoursMax = 26000,C_elyz = zeros(ny,ns)),
"fc" => (powerMax_ini = 0.,α_p = 8/100,η_H2_E = 0.4,η_H2_H = 0.4,soh_ini = 1.,lifetime=5,nHoursMax = 10000,C_fc = zeros(ny,ns)),
"h2tank" => (Erated_ini = 0.,α_p_ch = 1.5,α_p_dch = 1.5,η_ch =0.8,η_dch = 0.8,lifetime=25,
η_self = 0.1,α_soc_min = 0.2,α_soc_max = 0.8,α_soc_ini = 0.8,C_tank = zeros(ny,ns)),
"tes" => (Erated_ini = 0.,α_p_ch = 1.5,α_p_dch = 1.5,η_ch = 0.8,η_dch = 0.8,lifetime=25,
η_self = 0.1,α_soc_min = 0.2,α_soc_max = 0.8,α_soc_ini = 0.8,C_tes = zeros(ny,ns)),
"heater" => (powerMax_ini = 10, η_E_H = 1.,lifetime=25,C_heater = zeros(ny,ns)),
"controller" => [],
"designer" => [],
"grid" => [],
)

# Initialization
liion = Genesys.Liion(outputGUI["liion"], nh, ny, ns)
h2tank = Genesys.H2Tank(outputGUI["h2tank"], nh, ny, ns)
tes = Genesys.ThermalSto(outputGUI["tes"], nh, ny, ns)
fc = Genesys.FuelCell(outputGUI["fc"], nh, ny, ns)
elyz = Genesys.Electrolyzer(outputGUI["elyz"], nh, ny, ns)
heater = Genesys.Heater(outputGUI["heater"], nh, ny, ns)

#-------------------------------------------------------------------------------
# OPERATION
#-------------------------------------------------------------------------------

# Liion test set
@testset "Battery compute_operation_dynamics tests" begin
    @testset "Battery charge/discharge" begin
        x_liion = (Erated=100, soc=0.5, soh=1)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 10., Δh) == (0.5 * (1-0.1) - 10 / 0.8 * 1 / 100, 1-10/2/2500/0.6/100, 10)
        @test Genesys.compute_operation_dynamics(liion, x_liion, -10., Δh) == (0.5 * (1-0.1) + 0.8*10 * 1 / 100, 1-10/2/2500/0.6/100, -10)
    end
    @testset "Battery SoC bounds" begin
        x_liion = (Erated=100, soc=0.5, soh=1)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 100., Δh) == (0.5 *(1-0.1), 1.0, 0.0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, -100., Δh) == (0.5 *(1-0.1), 1.0, 0.0)
    end
    @testset "Battery power bounds" begin
        x_liion = (Erated=100, soc=0.5, soh=1)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 1000000., Δh) == (0.5 * (1-0.1), 1.0, 0.0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, -1000000., Δh) == (0.5 * (1-0.1), 1.0, 0.0)
    end
    @testset "Battery soh" begin
        x_liion = (Erated=100, soc=0.5, soh=0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 100., Δh) == (0.5* (1-0.1), 0, 0.0)
        x_liion = (Erated=100, soc=0.2, soh=0)
        @test Genesys.compute_operation_dynamics(liion, x_liion, 100., Δh) == (0.2, 0, 0.0)
    end
end

# H2 tank test set
@testset "H2 tank  compute_operation_dynamics tests" begin
    @testset "H2 tank charge/discharge" begin
        x_h2 = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, 10., Δh) == (0.5 * (1-0.1) - 10 / 0.8 * 1 / 100, 10)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, -10., Δh) == (0.5 * (1-0.1) + 0.8*10 * 1 / 100, -10)
    end
    @testset "H2 tank  SoC bounds" begin
        x_h2 = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, 100., Δh) == (0.5 *(1-0.1), 0.0)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, -100., Δh) == (0.5 *(1-0.1), 0.0)
    end
    @testset "H2 tank  power bounds" begin
        x_h2 = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, 1000000., Δh) == (0.5 * (1-0.1), 0.0)
        @test Genesys.compute_operation_dynamics(h2tank, x_h2, -1000000., Δh) == (0.5 * (1-0.1), 0.0)
    end
end

# TES test set
@testset "TES  compute_operation_dynamics tests" begin
    @testset "TES charge/discharge" begin
        x_tes = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(tes, x_tes, 10., Δh) == (0.5 * (1-0.1) - 10 / 0.8 * 1 / 100, 10)
        @test Genesys.compute_operation_dynamics(tes, x_tes, -10., Δh) == (0.5 * (1-0.1) + 0.8*10 * 1 / 100, -10)
    end
    @testset "TES  SoC bounds" begin
        x_tes = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(tes, x_tes, 100., Δh) == (0.5 *(1-0.1), 0.0)
        @test Genesys.compute_operation_dynamics(tes, x_tes, -100., Δh) == (0.5 *(1-0.1), 0.0)
    end
    @testset "TES  power bounds" begin
        x_tes = (Erated=100, soc=0.5)
        @test Genesys.compute_operation_dynamics(tes, x_tes, 1000000., Δh) == (0.5 * (1-0.1), 0.0)
        @test Genesys.compute_operation_dynamics(tes, x_tes, -1000000., Δh) == (0.5 * (1-0.1), 0.0)
    end
end
# Heater test set
@testset "Heater compute_operation_dynamics tests" begin
    @testset "Heater power" begin
        x_heater = (powerMax=100,)
        @test Genesys.compute_operation_dynamics(heater, x_heater, -10., Δh) == (-10, 10)
        @test Genesys.compute_operation_dynamics(heater, x_heater, 10., Δh) == (0, 0)
    end
    @testset "Heater bounds" begin
        x_heater = (powerMax=100,)
        @test Genesys.compute_operation_dynamics(heater, x_heater, -1000., Δh) == (0, 0)
    end
end
#  Fuel Cell test set
@testset "Fuel Cell compute_operation_dynamics tests" begin
    @testset "fc power" begin
        x_fc = (powerMax=100., soh=1.)
        @test Genesys.compute_operation_dynamics(fc, x_fc, 50., Δh) == (50, 50*0.4/0.4, -50/0.4, 1-1/10000)
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
        @test Genesys.compute_operation_dynamics(elyz, x_elyz, -50., Δh) == (-50, 50*0.3, 50*0.5, 1-1/26000)
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



#-------------------------------------------------------------------------------
# INVESTMENT
#-------------------------------------------------------------------------------

# Liion test set
@testset "Battery compute_investment_dynamics tests" begin
    x_liion=(Erated=200, soc=0.5, soh=1)
    @test Genesys.compute_investment_dynamics(liion, x_liion, 300) == (300, 0., 1) # liion.soc[1,1] = 0 in that case...
    @test Genesys.compute_investment_dynamics(liion, x_liion, 0) == (200, 0.5, 1)
end
# H2 tank test set
@testset "H2tank compute_investment_dynamics tests" begin
    x_tank=(Erated=200, soc=0.5)
    @test Genesys.compute_investment_dynamics(h2tank, x_tank, 300) == (300, 0.)
    @test Genesys.compute_investment_dynamics(h2tank, x_tank, 0) == (200, 0.5)
end
# TES test set
@testset "TES compute_investment_dynamics tests" begin
    x_tes=(Erated=200, soc=0.5)
    @test Genesys.compute_investment_dynamics(tes, x_tes, 300) == (300, 0.)
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
