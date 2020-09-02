# This file includes all the plot functions for both DES

# Simple
function plot_operation(ld::Load, pv::Source, liion::Liion, grid::Grid, parameters::NamedTuple; y=2, s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="deep", font="serif")

    # Parameters
    hours = (parameters.Δh:parameters.Δh:parameters.H) / parameters.Δh

    # Plots
    figure("Operation")
    # Node E
    sp=subplot(211)
    plot(hours,ld.power_E[:,y,s],label="ld_E_TOT", color=(0., 0.098, 0.196))
    plot(hours,pv.power_E[:,y,s],label="pv", color= (1., 0.784, 0.588))
    plot(hours,liion.power_E[:,y,s],label="liion", color=(0., 0.588, 0.314))
    plot(hours,grid.power_E[:,y,s],label="grid", color=(0.588, 0.588, 0.588))
    ylabel("POWER (kW)", weight="bold")

    subplot(212, sharex=sp)
    plot(hours,liion.soc[1:end-1,y,s], color=(0., 0.588, 0.314))
    ylabel("BATTERY SOC", weight="bold")
end
function plot_investment(designer::AbstractDesigner, parameters::NamedTuple; s=1)
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif", font_scale=1.5)

    # Years
    years = (parameters.Δy:parameters.Δy:parameters.Y) / parameters.Δy

    # Bar plot
    figure("SIZING")
    sp=subplot(211)
    bar(years, designer.u.u_pv[:,s], color="sandybrown")
    ylabel("SOLAR \n PEAK POWER (kWp)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), #ylim(0,120)
    setp(sp.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    subplot(212)
    bar(years, designer.u.u_liion[:,s], color="steelblue")
    ylabel("BATTERY \n CAPACITY(kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"),
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
end
function plot_soh(liion; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif", font_scale=2.5)
    # Parameters
    nh = size(liion.power_E,1)
    hours = collect(1:length(reshape(liion.soh[:,:,s],:,1))) / nh

    # Liion
    figure("SoH")
    plot(hours,reshape(liion.soh[:,:,s],:,1), linewidth=3, color = "darkred")
    ylabel("BATTERY SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    xlabel("YEARS", weight = "black", size = "large"), xlim(1,20), xticks([1,5,10,15,20], weight = "black", size = "large")
end

# Multi-energy
function plot_operation(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
     elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
     grid::Grid, parameters::NamedTuple; y=2, s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="deep", font="serif")

    # Parameters
    hours = (parameters.Δh:parameters.Δh:parameters.H) / parameters.Δh

    # Plots
    figure("Operation")
    # Node E
    sp=subplot(321)
    plot(hours,ld.power_E[:,y,s],label="ld_E", color=(0., 0.098, 0.196), linestyle="--")
    plot(hours,ld.power_E[:,y,s]-heater.power_E[:,y,s],label="ld_E_TOT", color=(0., 0.098, 0.196))
    plot(hours,pv.power_E[:,y,s],label="pv", color= (1., 0.784, 0.588))
    plot(hours,liion.power_E[:,y,s],label="liion", color=(0., 0.588, 0.314))
    plot(hours,elyz.power_E[:,y,s] + fc.power_E[:,y,s],label="h2hub", color=(0.39, 1., 0.71))
    plot(hours,grid.power_E[:,y,s],label="grid", color=(0.588, 0.588, 0.588))
    ylabel("ELEC. POWER (kW)", weight="bold")

    subplot(322, sharex=sp)
    plot(hours,liion.soc[1:end-1,y,s], color=(0., 0.588, 0.314))
    ylabel("BATTERY SOC", weight="bold")

    subplot(324, sharex=sp)
    plot(hours,h2tank.soc[1:end-1,y,s], color=(0.39, 1., 0.71))
    ylabel("H2 SOC", weight="bold")

    # Node H
    sp=subplot(323, sharex=sp)
    plot(hours,ld.power_H[:,y,s],label="ld_H", color=(0., 0.098, 0.196))
    plot(hours,elyz.power_H[:,y,s] + fc.power_H[:,y,s],label="h2hub", color=(1., 0.588, 0.588))
    plot(hours,tes.power_H[:,y,s],label="tes", color=(1., 0., 0.))
    plot(hours,heater.power_H[:,y,s],label="heater", color=(0.39, 0., 0.))
    ylabel("HEAT. POWER (kW)", weight="bold")

    subplot(326, sharex=sp)
    plot(hours,tes.soc[1:end-1,y,s], color=(1., 0., 0.))
    ylabel("TES SOC", weight="bold")
    xlabel("HOURS", weight="bold")
end
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
function plot_investment_multi_energy(designer::AbstractDesigner, parameters::NamedTuple; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif", font_scale=1)

    # Years
    years = (parameters.Δy:parameters.Δy:parameters.Y) / parameters.Δy

    # Bar plot
    figure("SIZING")
    # PV
    sp=subplot(321)
    bar(years, designer.u.u_pv[:,s], color="sandybrown")
    ylabel("SOLAR \n PEAK POWER (kWp)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    setp(sp.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # H2 Tank
    sp1=subplot(322, sharex=sp)
    bar(years, designer.u.u_tank[:,s], color="lightgreen")
    ylabel("H2 \n CAPACITY (kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    setp(sp1.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # Liion
    sp2=subplot(323, sharex=sp)
    bar(years, designer.u.u_liion[:,s], color="steelblue")
    ylabel("BATTERY \n CAPACITY(kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # Electrolyzer
    sp3=subplot(324, sharex=sp)
    bar(years, designer.u.u_elyz[:,s], color="lightgreen")
    setp(sp3.get_xticklabels(), visible=false), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    ylabel("ELYZ \n MAX. POWER (kW)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    # TES
    subplot(325, sharex=sp)
    bar(years, designer.u.u_tes[:,s], color="lightcoral")
    ylabel("TES \n CAPACITY (kWh)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
    # Fuel cell
    subplot(326, sharex=sp)
    bar(years, designer.u.u_fc[:,s], color="lightgreen")
    ylabel("FC \n MAX. POWER (kW)", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(bottom=0)
    xlabel("YEARS", weight = "black", size = "large"), xlim(0,20), xticks(0:2:20, weight = "black", size = "large")
end
function plot_soh(liion::Liion, elyz::Electrolyzer, fc::FuelCell; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif", font_scale=1.5)
    # Parameters
    nh = size(liion.power_E,1)
    hours = collect(1:length(reshape(liion.soh[:,:,s],:,1))) / nh

    # Liion
    figure("SoH")
    # Liion
    sp=subplot(311)
    plot(hours,reshape(liion.soh[:,:,s],:,1),linewidth=3, color = "darkred")
    ylabel("BATTERY SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    # Electrolyzer
    subplot(312, sharex=sp)
    plot(hours,reshape(elyz.soh[:,:,s],:,1),linewidth=3, color = "darkred")
    ylabel("ELYZ SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    # Fuel Cell
    subplot(313, sharex=sp)
    plot(hours,reshape(fc.soh[:,:,s],:,1),linewidth=3, color = "darkred")
    ylabel("FC SOH", weight = "black", size = "large"), yticks(weight = "black", size = "medium"), ylim(0,1)
    xlabel("YEARS", weight = "black", size = "large"), xlim(1,20), xticks([1,5,10,15,20], weight = "black", size = "large")
end

# Economics
function plot_economics(costs::Costs, parameters::NamedTuple; s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif",font_scale = 2.5)

    # Horizon
    Y = length(costs.capex[:,s])
    γ = 1. ./ (1. + parameters.τ) .^ (1:parameters.Y) # discount factor

    # Plots
    figure("CASH FLOWS")
    bar(1:Y, -γ .* costs.capex[:,s] ./ 1000, label="Investment", color="steelblue")
    bar(1:Y, γ .* costs.opex[:,s] ./ 1000, label="Income", color="coral")
    ylabel("CASH FLOWS (k€)", weight = "bold"), yticks(weight = "bold")
    xlabel("YEARS", weight = "bold"), xticks(0:5:20, weight = "bold"), xlim(0,21)
    legend(fontsize="xx-large", edgecolor="inherit")
    grid()

    figure("CUMULATIVE NPV")
    bar(1:Y, costs.cumulative_npv[:,s] ./ 1000, color="steelblue")
    ylabel("CUMULATIVE NPV (k€)", weight = "bold"), yticks(weight = "bold")
    xlabel("YEARS", weight = "bold"), xticks(0:5:20, weight = "bold"), xlim(0,21)
    grid()
end
function plot_npv(costs)
    Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif", font_scale=1.5)

    # Histogram
    figure("NPV")
    hist(costs.npv / 1000, color="darkred")
    ylabel("SCENARIO COUNT", weight = "black", size = "large"), yticks(weight = "black", size = "medium")
    xlabel("NPV (k€)", weight = "black", size = "large"), xticks(weight = "black", size = "medium")
end
