# This file includes all the plot functions for the DES

function plot_operation(des::DistributedEnergySystem ; y=2, s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)

    # Parameters
    nh = des.parameters.nh
    Δh = des.parameters.Δh
    hours = range(1, length = nh, step = Δh) / Δh

    isa(des.ld_E, Load) ? ld_E = des.ld_E.power[:,y,s] : ld_E = zeros(nh)
    isa(des.ld_H, Load) ? ld_H = des.ld_H.power[:,y,s] : ld_H = zeros(nh)
    isa(des.pv, Source) ? pv = des.pv.power_E[:,y,s] : pv = zeros(nh)
    isa(des.grid, Grid) ? grid = des.grid.power_E[:,y,s] : grid = zeros(nh)
    if isa(des.liion, Liion)
        liion = des.liion.power_E[:,y,s]
        soc_liion = des.liion.soc[1:end-1,y,s]
    else
        liion, soc_liion = zeros(nh), zeros(nh)
    end
    if isa(des.tes, ThermalSto)
        tes = des.tes.power_H[:,y,s]
        soc_tes = des.tes.soc[1:end-1,y,s]
    else
        tes, soc_tes = zeros(nh), zeros(nh)
    end
    if isa(des.h2tank, H2Tank)
        h2tank = des.h2tank.power_H2[:,y,s]
        soc_h2tank = des.h2tank.soc[1:end-1,y,s]
    else
        h2tank, soc_h2tank = zeros(nh), zeros(nh)
    end
    if isa(des.elyz, Electrolyzer)
        elyz_E = des.elyz.power_E[:,y,s]
        elyz_H = des.elyz.power_H[:,y,s]
    else
        elyz_E, elyz_H = zeros(nh), zeros(nh)
    end
    if isa(des.fc, FuelCell)
        fc_E = des.fc.power_E[:,y,s]
        fc_H = des.fc.power_H[:,y,s]
    else
        fc_E, fc_H = zeros(nh), zeros(nh)
    end
    if isa(des.heater, Heater)
        heater_E = des.heater.power_E[:,y,s]
        heater_H = des.heater.power_H[:,y,s]
    else
        heater_E, heater_H = zeros(nh), zeros(nh)
    end

    # Plots
    figure("POWERS")
    # Node E
    sp=subplot(211)
    plot(hours, ld_E, label="ld_E", color=(0., 0.098, 0.196), linestyle="--")
    plot(hours, ld_E - heater_E, label="ld_E_TOT", color=(0., 0.098, 0.196))
    plot(hours, pv, label="pv", color= (1., 0.784, 0.588))
    plot(hours, liion, label="liion", color=(0., 0.588, 0.314))
    plot(hours, elyz_E + fc_E, label="h2hub", color=(0.39, 1., 0.71))
    plot(hours, grid, label="grid", color=(0.588, 0.588, 0.588))
    ylabel("ELEC. POWER (kW)", weight="bold")
    # Node H
    sp=subplot(212, sharex=sp)
    plot(hours, ld_H, label="ld_H", color=(0., 0.098, 0.196))
    plot(hours, elyz_H + fc_H, label="h2hub", color=(1., 0.588, 0.588))
    plot(hours, tes, label="tes", color=(1., 0., 0.))
    plot(hours, heater_H, label="heater", color=(0.39, 0., 0.))
    ylabel("HEAT. POWER (kW)", weight="bold")

    figure("SoC")
    subplot(311, sharex=sp)
    plot(hours, soc_liion, color=(0., 0.588, 0.314))
    ylabel("BATTERY SOC", weight="bold")

    subplot(312, sharex=sp)
    plot(hours, soc_h2tank, color=(0.39, 1., 0.71))
    ylabel("H2 SOC", weight="bold")

    subplot(313, sharex=sp)
    plot(hours, soc_tes, color=(1., 0., 0.))
    ylabel("TES SOC", weight="bold")
    xlabel("HOURS", weight="bold")
end
function plot_investment(des::DistributedEnergySystem, designer::AbstractDesigner; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.)

    # Years
    years = range(1, length=des.parameters.ny, step=des.parameters.Δy) / des.parameters.Δy

    # Bar plot
    figure("SIZING")
    # PV
    sp=subplot(321)
    bar(years, designer.u.pv[:,s], color="sandybrown")
    ylabel("SOLAR \n PEAK POWER (kWp)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    setp(sp.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # H2 Tank
    sp1=subplot(322, sharex=sp)
    bar(years, designer.u.h2tank[:,s], color="lightgreen")
    ylabel("H2 \n CAPACITY (kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    setp(sp1.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # Liion
    sp2=subplot(323, sharex=sp)
    bar(years, designer.u.liion[:,s], color="steelblue")
    ylabel("BATTERY \n CAPACITY(kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # Electrolyzer
    sp3=subplot(324, sharex=sp)
    bar(years, designer.u.elyz[:,s], color="lightgreen")
    setp(sp3.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    ylabel("ELYZ \n MAX. POWER (kW)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    # TES
    subplot(325, sharex=sp)
    bar(years, designer.u.tes[:,s], color="lightcoral")
    ylabel("TES \n CAPACITY (kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # Fuel cell
    subplot(326, sharex=sp)
    bar(years, designer.u.fc[:,s], color="lightgreen")
    ylabel("FC \n MAX. POWER (kW)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
end
function plot_soh(des::DistributedEnergySystem; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)
    # Parameters
    len = (des.parameters.nh + 1 ) * (des.parameters.ny + 1)
    hours = (1:len) / 8760

    if isa(des.liion, Liion)
        soh_liion = reshape(des.liion.soh[:,:,s],:,1)
    else
        soh_liion = zeros(len)
    end
    if isa(des.elyz, Electrolyzer)
        soh_elyz = reshape(des.elyz.soh[:,:,s],:,1)
    else
        soh_elyz = zeros(len)
    end
    if isa(des.fc, FuelCell)
        soh_fc = reshape(des.fc.soh[:,:,s],:,1)
    else
        soh_fc = zeros(len)
    end

    figure("SoH")
    # Liion
    sp=subplot(311)
    plot(hours, soh_liion, linewidth=3, color = "darkred")
    ylabel("BATTERY SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    # Electrolyzer
    subplot(312, sharex=sp)
    plot(hours, soh_elyz, linewidth=3, color = "darkred")
    ylabel("ELYZ SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    # Fuel Cell
    subplot(313, sharex=sp)
    plot(hours, soh_fc, linewidth=3, color = "darkred")
    ylabel("FC SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    xlabel("YEARS", weight = "black", size = "large"), xlim(1,20), xticks([1,5,10,15,20], weight = "black", size = "large")
end
# Economics
function plot_costs(costs::Costs; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale = 1.5)

    # Horizon
    years = 1:size(costs.capex, 1)

    # Plots
    figure("CASH FLOWS")
    bar(years, costs.capex[:,s] ./ 1000, label="Investment", color="steelblue")
    bar(years, costs.opex[:,s] ./ 1000, label="Income", color="coral")
    ylabel("CASH FLOWS (k€)", weight = "bold"), yticks(weight = "bold")
    xlabel("YEARS", weight = "bold"), xticks(0:5:20, weight = "bold"), xlim(0,21)
    legend(fontsize="xx-large", edgecolor="inherit")
    grid()

    figure("CUMULATIVE NPV")
    bar(years, costs.cumulative_npv[:,s] ./ 1000, color="steelblue")
    ylabel("CUMULATIVE NPV (k€)", weight = "bold"), yticks(weight = "bold")
    xlabel("YEARS", weight = "bold"), xticks(0:5:20, weight = "bold"), xlim(0,21)
    grid()

    figure("NPV")
    hist(costs.npv / 1000, color="darkred")
    ylabel("SCENARIO COUNT", weight = "black", size = "large"), yticks(weight = "black", size = "medium")
    xlabel("NPV (k€)", weight = "black", size = "large"), xticks(weight = "black", size = "medium")
end

# Discarded
function plot_operation_stack(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
     elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
     grid::Grid; y=2, s=1, hmin=1, hmax=8760)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif",font_scale = 2)

    # Parameters
    hours = hmin:hmax

    # Post treatment
    load_E = ld.power_E[hours,y,s] - heater.power_E[hours,y,s]
    load_H = ld.power_H[hours,y,s]

    # Prod
    # Node E
    selfconso = min.(load_E, pv.power_E[hours,y,s])
    liion_in = liion.power_E[hours,y,s] .* (liion.power_E[hours,y,s] .> 0)
    fc_in_E = fc.power_E[hours,y,s]
    grid_in = grid.power_E[hours,y,s]

    # Node H
    coge = fc.power_H[hours,y,s] + elyz.power_H[hours,y,s]
    tes_in = tes.power_H[hours,y,s] .* (tes.power_H[hours,y,s] .> 0)
    heater_in = heater.power_H[hours,y,s]

    # Conso
    # Node E
    liion_out = liion.power_E[hours,y,s] .* (liion.power_E[hours,y,s] .< 0)
    elyz_out = elyz.power_E[hours,y,s]
    curtail = -(selfconso - pv.power_E[hours,y,s] - liion_out - elyz_out)

    # Node H
    tes_out = tes.power_H[hours,y,s] .* (tes.power_H[hours,y,s] .< 0)

    # Plots
    figure("Power")
    # Node E
    sp=subplot(211)
    plot(hours,load_E ,label="Elec. load", color="black", linewidth=1.5, linestyle="--")
    stackplot(hours, selfconso, curtail, liion_in, fc_in_E, grid_in, labels=["PV self-cons.", "Curtail.", "Li-ion", "H2 hub", "Grid"], colors=["sandybrown", "bisque","forestgreen", "lightgreen","lightgray"])
    stackplot(hours, liion_out, elyz_out, colors=["forestgreen", "lightgreen"])
    setp(sp.get_xticklabels(), visible=false), xlim(hours[1],hours[end])
    ylabel("POWER (kW)", weight="black", size="large"), yticks(weight="black")
    legend(loc="lower center",edgecolor="inherit", ncol=6, mode="expand")
    sp.grid()


    # Node H
    sp=subplot(212, sharex=sp)
    plot(hours,load_H ,label="Thermal load", color="black", linewidth=1.5, linestyle="--")
    stackplot(hours, coge, tes_in, heater_in, labels=["H2 Coge.", "TES", "Heater"], colors=["firebrick", "lightcoral", "peachpuff"])
    stackplot(hours, tes_out, colors=["lightcoral"])
    xlabel("HOURS", weight="black", size="large"), xticks(weight="black")
    ylabel("POWER (kW)", weight="black", size="large"), yticks(weight="black")
    legend(loc="lower center",edgecolor="inherit", ncol=4, mode="expand")
    sp.grid()

    f=figure("SoC")
    f.text(0.05, 0.4, "SoC (%)", ha="center", rotation="vertical", weight = "black", size = 20 )
    sp=subplot(311)
    plot(hours,liion.soc[hours,y,s] * 100, label="Li-ion", color="forestgreen", linewidth=2.5)
    yticks(weight="black")
    setp(sp.get_xticklabels(), visible=false), xlim(hours[1],hours[end])
    legend(fontsize="large",edgecolor="inherit")
    sp.grid()

    sp=subplot(312, sharex=sp)
    plot(hours,h2tank.soc[hours,y,s] * 100, label="H2 tank", color="lightgreen", linewidth=2.5)
    setp(sp.get_xticklabels(), visible=false)
    yticks(weight="black")
    legend(fontsize="large",edgecolor="inherit")
    sp.grid()

    sp=subplot(313, sharex=sp)
    plot(hours,tes.soc[hours,y,s] * 100, label="TES", color="lightcoral", linewidth=2.5)
    yticks(weight="black")
    xlabel("HOURS", weight="black", size="large"), xticks(weight="black")
    legend(fontsize="large",edgecolor="inherit")
    sp.grid()
end
