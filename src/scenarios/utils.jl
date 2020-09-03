# Find indices of weekend days
isweekend(timestamp::DateTime) = (Dates.dayname.(timestamp) .== "Saturday") .| (Dates.dayname.(timestamp) .== "Sunday")
# Plot statistics
function plot_stats(ω, data)
    Seaborn.set(context="notebook",style="ticks",palette="deep", font="serif", font_scale=1.5)

    # Reshape
    ω_day = reshape(ω, 24, :)
    ω_month = reshape(ω, :, 12)
    data_day = reshape(data, 24, :)
    data_month = reshape(data, :, 12)

    # Stats
    # Mean
    μ_ω_day = mean(ω_day, dims=2)
    μ_data_day = mean(data_day, dims=2)
    μ_ω_month = vec(mean(ω_month, dims=1))
    μ_data_month = vec(mean(data_month, dims=1))
    # Standard deviation
    σ_ω_day = std(ω_day, dims=2)
    σ_data_day = std(data_day, dims=2)
    σ_ω_month = vec(std(ω_month, dims=1))
    σ_data_month = vec(std(data_month, dims=1))


    # Plots
    subplot(121)
    plot(1:24, μ_ω_day, color="gray", label="μ model")
    plot(1:24, μ_data_day, color="darkred", label="μ data")
    plot(1:24, σ_ω_day, color="gray", linestyle="--", label="σ model")
    plot(1:24, σ_data_day, color="darkred", linestyle="--", label="σ data")
    xlabel("Hours", weight = "bold"), xlim(0,24)
    title("Daily indicators", weight = "bold")
    legend()
    subplot(122)
    plot(1:12,μ_ω_month, color="gray", label="μ model")
    plot(1:12,μ_data_month, color="darkred", label="μ data")
    plot(1:12,σ_ω_month, color="gray", linestyle="--", label="σ model")
    plot(1:12,σ_data_month, color="darkred", linestyle="--", label="σ data")
    xlabel("Month", weight = "bold"), xlim(1,12)
    title("Monthly indicators", weight = "bold")
    legend()
end
