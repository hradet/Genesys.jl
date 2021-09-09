# This file includes all the plot functions for the MG
function plot_operation(mg::Microgrid ; y=2, s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)

    # Parameters
    nh = mg.parameters.nh
    Δh = mg.parameters.Δh
    hours = range(1, length = nh, step = Δh) / Δh

    # Plots
    # Powers
    f = figure("Powers")
    for (i, type) in enumerate([typeof(Electricity()), typeof(Heat()), typeof(Hydrogen())])
        i == 1 ? subplot(3, 1, i) : subplot(3, 1, i, sharex = f.axes[1])
        # Demands
        for (k, a) in enumerate(mg.demands)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Demand ", k))
            end
        end
        # Generations
        for (k, a) in enumerate(mg.generations)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Generation ", k))
            end
        end
        # Storages
        for (k, a) in enumerate(mg.storages)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Storage ", k))
            end
        end
        # Converters
        for (k, a) in enumerate(mg.converters)
            for c in a.carrier
                if c isa type
                    plot(hours, c.power[:,y,s], label = string("Converter ", k))
                end
            end
        end
        for (k, a) in enumerate(mg.grids)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Grids ", k))
            end
        end
        legend()
    end
    # State of charge
    figure("State-of-charge")
    for (k, a) in enumerate(mg.storages)
        k == 1 ? subplot(length(mg.storages), 1, k) : subplot(length(mg.storages), 1, k, sharex = f.axes[1])
        plot(hours, a.soc[1:end-1, y, s], label = string("Storage ", k))
        legend()
    end
end
function plot_operation(mg::Microgrid, controller::AbstractController; y=2, s=1)
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)

    # Parameters
    nh = mg.parameters.nh
    Δh = mg.parameters.Δh
    hours = range(1, length = nh, step = Δh) / Δh

    # Plots
    # Powers
    f = figure("Powers")
    for (i, type) in enumerate([typeof(Electricity()), typeof(Heat()), typeof(Hydrogen())])
        i == 1 ? subplot(3, 1, i) : subplot(3, 1, i, sharex = f.axes[1])
        # Demands
        for (k, a) in enumerate(mg.demands)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Demand ", k))
            end
        end
        # Generations
        for (k, a) in enumerate(mg.generations)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Generation ", k))
            end
        end
        # Storages
        for (k, a) in enumerate(mg.storages)
            if a.carrier isa type
                plot(hours, controller.decisions.storages[k][:,y,s], label = string("Storage ", k))
            end
        end
        # Converters
        for (k, a) in enumerate(mg.converters)
            for c in a.carrier
                if c isa type
                    plot(hours, controller.decisions.converters[k][:,y,s], label = string("Converter ", k))
                end
            end
        end
        for (k, a) in enumerate(mg.grids)
            if a.carrier isa type
                plot(hours, a.carrier.power[:,y,s], label = string("Grids ", k))
            end
        end
        legend()
    end
end
# Statistics
function plot_metrics(metrics::Metrics; type = "hist")
    # Seaborn configuration
    Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale = 1.5)

    if type == "hist"
        figure("Renewable share")
        hist(reshape(metrics.renewable_share[2:end, :], :) * 100, color="steelblue")
        ylabel("Counts", size = "large"), yticks(size = "medium")
        xlabel("Renewable share (%)", size = "large"), xticks(size = "medium")
        figure("Cost")
        hist(reshape(metrics.eac.total[:, :], :) / 1000, color="steelblue")
        ylabel("Counts", size = "large"), yticks(size = "medium")
        xlabel("Annual cost (k€/y)", size = "large"), xticks(size = "medium")
        if !isa(metrics.lpsp.elec, Nothing)
            figure("LPSP elec")
            hist(reshape(metrics.lpsp.elec[2:end, :], :) * 100, color="steelblue")
            ylabel("Counts", size = "large"), yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.heat, Nothing)
            figure("LPSP heat")
            hist(reshape(metrics.lpsp.heat[2:end, :], :) * 100, color="steelblue")
            ylabel("Counts", size = "large"), yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.hydrogen, Nothing)
            figure("LPSP hydrogen")
            hist(reshape(metrics.lpsp.hydrogen[2:end, :], :) * 100, color="steelblue")
            ylabel("Counts", size = "large"), yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
    elseif type == "violin"
        figure("Renewable share")
        violinplot(reshape(metrics.renewable_share[2:end, :], :) * 100, color="steelblue")
        yticks(size = "medium")
        xlabel("Renewable share (%)", size = "large"), xticks(size = "medium")
        figure("Cost")
        violinplot(reshape(metrics.eac.total[:, :], :) / 1000, color="steelblue")
        yticks(size = "medium")
        xlabel("Annual cost (k€/y)", size = "large"), xticks(size = "medium")
        if !isa(metrics.lpsp.elec, Nothing)
            figure("LPSP elec")
            violinplot(reshape(metrics.lpsp.elec[2:end, :], :) * 100, color="steelblue")
            yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.heat, Nothing)
            figure("LPSP heat")
            violinplot(reshape(metrics.lpsp.heat[2:end, :], :) * 100, color="steelblue")
            yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
        if !isa(metrics.lpsp.hydrogen, Nothing)
            figure("LPSP hydrogen")
            violinplot(reshape(metrics.lpsp.hydrogen[2:end, :], :) * 100, color="steelblue")
            yticks(size = "medium")
            xlabel("LPSP (%)", size = "large"), xticks(size = "medium")
        end
    else
        println("Only 'hist' or 'violin' type accepted")
    end
end

# Discarded
# function plot_operation_stack(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
#      elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
#      grid::Grid; y=2, s=1, hmin=1, hmax=8760)
#     # Seaborn configuration
#     Seaborn.set(context="notebook",style="ticks",palette="muted", font="serif",font_scale = 2)
#
#     # Parameters
#     hours = hmin:hmax
#
#     # Post treatment
#     load_E = ld.power_E[hours,y,s] - heater.power_E[hours,y,s]
#     load_H = ld.power_H[hours,y,s]
#
#     # Prod
#     # Node E
#     selfconso = min.(load_E, pv.power_E[hours,y,s])
#     liion_in = liion.power_E[hours,y,s] .* (liion.power_E[hours,y,s] .> 0)
#     fc_in_E = fc.power_E[hours,y,s]
#     grid_in = grid.power_E[hours,y,s]
#
#     # Node H
#     coge = fc.power_H[hours,y,s] + elyz.power_H[hours,y,s]
#     tes_in = tes.power_H[hours,y,s] .* (tes.power_H[hours,y,s] .> 0)
#     heater_in = heater.power_H[hours,y,s]
#
#     # Conso
#     # Node E
#     liion_out = liion.power_E[hours,y,s] .* (liion.power_E[hours,y,s] .< 0)
#     elyz_out = elyz.power_E[hours,y,s]
#     curtail = -(selfconso - pv.power_E[hours,y,s] - liion_out - elyz_out)
#
#     # Node H
#     tes_out = tes.power_H[hours,y,s] .* (tes.power_H[hours,y,s] .< 0)
#
#     # Plots
#     figure("Power")
#     # Node E
#     sp=subplot(211)
#     plot(hours,load_E ,label="Elec. load", color="black", linewidth=1.5, linestyle="--")
#     stackplot(hours, selfconso, curtail, liion_in, fc_in_E, grid_in, labels=["PV self-cons.", "Curtail.", "Li-ion", "H2 hub", "Grid"], colors=["sandybrown", "bisque","forestgreen", "lightgreen","lightgray"])
#     stackplot(hours, liion_out, elyz_out, colors=["forestgreen", "lightgreen"])
#     setp(sp.get_xticklabels(), visible=false), xlim(hours[1],hours[end])
#     ylabel("POWER (kW)", weight="black", size="large"), yticks(weight="black")
#     legend(loc="lower center",edgecolor="inherit", ncol=6, mode="expand")
#     sp.grid()
#
#
#     # Node H
#     sp=subplot(212, sharex=sp)
#     plot(hours,load_H ,label="Thermal load", color="black", linewidth=1.5, linestyle="--")
#     stackplot(hours, coge, tes_in, heater_in, labels=["H2 Coge.", "TES", "Heater"], colors=["firebrick", "lightcoral", "peachpuff"])
#     stackplot(hours, tes_out, colors=["lightcoral"])
#     xlabel("HOURS", weight="black", size="large"), xticks(weight="black")
#     ylabel("POWER (kW)", weight="black", size="large"), yticks(weight="black")
#     legend(loc="lower center",edgecolor="inherit", ncol=4, mode="expand")
#     sp.grid()
#
#     f=figure("SoC")
#     f.text(0.05, 0.4, "SoC (%)", ha="center", rotation="vertical", weight = "black", size = 20 )
#     sp=subplot(311)
#     plot(hours,liion.soc[hours,y,s] * 100, label="Li-ion", color="forestgreen", linewidth=2.5)
#     yticks(weight="black")
#     setp(sp.get_xticklabels(), visible=false), xlim(hours[1],hours[end])
#     legend(fontsize="large",edgecolor="inherit")
#     sp.grid()
#
#     sp=subplot(312, sharex=sp)
#     plot(hours,h2tank.soc[hours,y,s] * 100, label="H2 tank", color="lightgreen", linewidth=2.5)
#     setp(sp.get_xticklabels(), visible=false)
#     yticks(weight="black")
#     legend(fontsize="large",edgecolor="inherit")
#     sp.grid()
#
#     sp=subplot(313, sharex=sp)
#     plot(hours,tes.soc[hours,y,s] * 100, label="TES", color="lightcoral", linewidth=2.5)
#     yticks(weight="black")
#     xlabel("HOURS", weight="black", size="large"), xticks(weight="black")
#     legend(fontsize="large",edgecolor="inherit")
#     sp.grid()
# end
# function plot_cumulative_distribution(data, data_reduced)
#     # Parameters
#     nh, ny, ns = size(data)
#     # Cumulative distribution of initial data
#     cumu = sort(sum(hcat([data[:,:,s] for s in 1:ns]...), dims=1)[1,:])
#     proba = collect(1:ns*ny) ./ ns ./ ny
#     plot(cumu,proba)
#     # Cumulative distribution of reduced data
#     sum_data_reduced = sum(data_reduced, dims=1)[1,:]
#     cumu = sum_data_reduced[sortperm(sum_data_reduced)]
#     proba = cumsum(results.counts[sortperm(sum_data_reduced)]) / 1000
#     plot(cumu,proba, ds="steps-post")
# end
#
# function plot_operation(mg::Microgrid, designer::MILP ; y=2, s=1)
#     # Seaborn configuration
#     Seaborn.set(context="notebook", style="ticks", palette="muted", font="serif", font_scale=1.5)
#
#     # Parameters
#     nh = mg.parameters.nh
#     Δh = mg.parameters.Δh
#     hours = range(1, length = nh, step = Δh) / Δh
#
#     ld_E = designer.history.ld_E.power[:,1,s]
#     ld_H = designer.history.ld_H.power[:,1,s]
#     pv = value(designer.model[:r_pv]) * designer.history.pv.power[:,1,s]
#     grid = (value.(designer.model[:p_g_in]) + value.(designer.model[:p_g_out]))[:,s]
#     liion = (value.(designer.model[:p_liion_ch]) + value.(designer.model[:p_liion_dch]))[:,s]
#     soc_liion = value.(designer.model[:soc_liion])[1:end-1,s] ./ value.(designer.model[:r_liion])
#     tes = (value.(designer.model[:p_tes_ch]) + value.(designer.model[:p_tes_dch]))[:,s]
#     soc_tes = value.(designer.model[:soc_tes])[1:end-1,s] ./ value.(designer.model[:r_tes])
#     h2tank = (value.(designer.model[:p_h2tank_ch]) + value.(designer.model[:p_h2tank_dch]))[:,s]
#     soc_h2tank = value.(designer.model[:soc_h2tank])[1:end-1,s] ./ value.(designer.model[:r_h2tank])
#     elyz_E = value.(designer.model[:p_elyz_E])[:,s]
#     elyz_H = - mg.elyz.η_E_H .* elyz_E
#     fc_E = value.(designer.model[:p_fc_E])[:,s]
#     fc_H = mg.fc.η_H2_H / mg.fc.η_H2_E .* fc_E
#     heater_E = value.(designer.model[:p_heater_E])[:,s]
#     heater_H = - mg.heater.η_E_H .* heater_E
#
#     # Plots
#     figure("POWERS MILP")
#     # Node E
#     sp=subplot(211)
#     plot(hours, ld_E, label="ld_E", color=(0., 0.098, 0.196), linestyle="--")
#     plot(hours, ld_E - heater_E, label="ld_E_TOT", color=(0., 0.098, 0.196))
#     plot(hours, pv, label="pv", color= (1., 0.784, 0.588))
#     plot(hours, liion, label="liion", color=(0., 0.588, 0.314))
#     plot(hours, elyz_E + fc_E, label="h2hub", color=(0.39, 1., 0.71))
#     plot(hours, grid, label="grid", color=(0.588, 0.588, 0.588))
#     ylabel("ELEC. POWER (kW)", weight="bold")
#     # Node H
#     sp=subplot(212, sharex=sp)
#     plot(hours, ld_H, label="ld_H", color=(0., 0.098, 0.196))
#     plot(hours, elyz_H + fc_H, label="h2hub", color=(1., 0.588, 0.588))
#     plot(hours, tes, label="tes", color=(1., 0., 0.))
#     plot(hours, heater_H, label="heater", color=(0.39, 0., 0.))
#     ylabel("HEAT. POWER (kW)", weight="bold")
#
#     figure("SoC MILP")
#     subplot(311, sharex=sp)
#     plot(hours, soc_liion, color=(0., 0.588, 0.314))
#     ylabel("BATTERY SOC", weight="bold")
#
#     subplot(312, sharex=sp)
#     plot(hours, soc_h2tank, color=(0.39, 1., 0.71))
#     ylabel("H2 SOC", weight="bold")
#
#     subplot(313, sharex=sp)
#     plot(hours, soc_tes, color=(1., 0., 0.))
#     ylabel("TES SOC", weight="bold")
#     xlabel("HOURS", weight="bold")
# end
