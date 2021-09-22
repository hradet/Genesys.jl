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
