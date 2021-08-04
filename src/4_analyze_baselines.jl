include("_common.jl")

using StatsPlots, OnlineStats


calib_start = 1826  # 1826
outflows = Dict()
climate = FULL_CLIMATE
calib_obs = CALIB_OBS[calib_start+1:end]

baseline_path = joinpath(FIG_PATH, "calib_baseline")

for app in APPROACHES
    display_name = replace(app, "_"=>" ")
    sn, n_id = setup_network(joinpath("baselines", "cotter_baseline_IHACRES_$(app).yml"))
    node = sn[n_id]
    run_basin!(sn, climate)

    calib_outflow = node.outflow[calib_start+1:CALIB_LN]
    gw_store = node.gw_store[calib_start+1:CALIB_LN]

    cmd = node.storage[calib_start+1:CALIB_LN]
    log_cmd = log.(cmd)
    max_cmd = maximum(log_cmd)
    quantiles = [0.2, 0.80]

    normal_dist = fit(Normal, log_cmd)
    μ = mean(normal_dist)
    σ = std(normal_dist)
    d_min = minimum(log_cmd)
    d_max = maximum(log_cmd)
    trunc_dist = truncated(Normal(μ, σ), d_min, d_max)
    
    thresholds = quantile(log_cmd, quantiles)
    @info "Log CMD thresholds: [0.1, 0.9 / min / max]" thresholds d_min d_max
    histogram(log_cmd, alpha=0.4, normalize=true,
              title="Distribution of Log CMD values\n($(display_name) calibrated baseline)",
              legend=:topleft,
              label=false,
              xlabel="Log CMD",
              ylabel="Density")
    StatsPlots.plot!(trunc_dist, linewidth=2, label="Fitted normal distribution")
    vline!(thresholds; color="orange", linewidth=1.5, label="State thresholds")
    for (thres, annot) in zip(thresholds, quantiles)
        annotate!(thres, 0.0015, text("q$annot", 8))
    end

    savefig(joinpath(baseline_path, "baseline_cmd_dist_$(app).png"))


    # Save performance over calibration period
    rmse = round(Streamfall.RMSE(calib_obs, calib_outflow), digits=2)
    mkge = round(Streamfall.mKGE(calib_obs, calib_outflow), digits=2)
    nse = round(Streamfall.NSE(calib_obs, calib_outflow), digits=2)
    nnse = round(Streamfall.NNSE(calib_obs, calib_outflow), digits=2)
    nmkge = round(Streamfall.NmKGE(calib_obs, calib_outflow), digits=2)
    mean_nmkge = round(Streamfall.mean_NmKGE(calib_obs, calib_outflow), digits=2)
    @info "$app calibration: "
    report_metrics(calib_obs, calib_outflow)

    plot(CALIB_DATES[calib_start+1:end], calib_obs, label="Historic", title="$app Calibration\n(RMSE: $rmse; NSE: $nse; KGE': $mkge)", legend=:topright)
    plot!(CALIB_DATES[calib_start+1:end], calib_outflow, label="Simulated", alpha=0.7, ylabel="Stream flow [ML]")
    savefig(joinpath(baseline_path, "baseline_calib_perf_$(app).png"))

    qqplot(calib_obs, calib_outflow, title="Calibration Period\n$display_name")
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(joinpath(baseline_path, "baseline_calib_perf_qq_$(display_name).png"))

    # cmd = node.storage[2:CALIB_LN+1]
    # mw = MovingWindow(window, Float64)
    # o = Quantile(quantiles, b=1000)
    # ma = [value(fit!(Mean(weight = ExponentialWeight(.01)), value(fit!(mw, i)))) for i in cmd]
    # thresholds = [value(fit!(o, i)) for i in ma]
    # plot(CALIB_DATES, cmd, label="CMD", title="$app CMD", legend=:topright)
    # plot!(CALIB_DATES, ma, label="MA Simulated", alpha=0.7, ylabel="CMD")
    # plot!(CALIB_DATES, [i[1] for i in thresholds]; color="orange", linewidth=1.5, label="State bounds")
    # plot!(CALIB_DATES, [i[2] for i in thresholds]; color="orange", linewidth=1.5, label="")
    # # plot!(CALIB_DATES, [i[3] for i in thresholds]; color="orange", linewidth=1.5, label="")
    # savefig(joinpath(baseline_path, "baseline_calib_cmd_$(app).png"))

    # logcmd = log.(cmd)
    # replace!(logcmd, -Inf=>0.0, Inf=>0.0)
    # o = Quantile(quantiles, b=1000)
    # ma = [value(fit!(Mean(weight = ExponentialWeight(.01)), value(fit!(mw, i)))) for i in logcmd]
    # thresholds = [value(fit!(o, i)) for i in ma]
    # plot(CALIB_DATES, logcmd, label="Log CMD", title="$app Log CMD", legend=:topright)
    # plot!(CALIB_DATES, ma, label="MA Log CMD", legend=:topright)
    # plot!(CALIB_DATES, [i[1] for i in thresholds]; color="orange", linewidth=1.5, label="State bounds")
    # plot!(CALIB_DATES, [i[2] for i in thresholds]; color="orange", linewidth=1.5, label="")
    # # plot!(CALIB_DATES, [i[3] for i in thresholds]; color="orange", linewidth=1.5, label="")
    # savefig(joinpath(baseline_path, "baseline_calib_logcmd_$(app).png"))

    # gw_store_ts = gw_store

    # log_gw = log.(gw_store_ts)
    # replace!(log_gw, -Inf=>0.0, Inf=>0.0)
    # mw = MovingWindow(window, Float64)
    # o = Quantile(quantiles, b=1000)
    # thresholds = [value(fit!(o, i)) for i in log_gw]
    # plot(CALIB_DATES, log_gw, label="", title="$app GW Store", legend=:topright)
    # plot!(CALIB_DATES, [i[1] for i in thresholds]; color="orange", linewidth=1.5, label="State bounds")
    # plot!(CALIB_DATES, [i[2] for i in thresholds]; color="orange", linewidth=1.5, label="")
    # # plot!(CALIB_DATES, [i[3] for i in thresholds]; color="orange", linewidth=1.5, label="")
    # savefig(joinpath(baseline_path, "baseline_calib_loggw_store_$(app).png"))

    @info "CMD" app minimum(node.storage[2:end]) maximum(node.storage[2:end])
    @info "GW Store" app minimum(node.gw_store[2:end]) maximum(node.gw_store[2:end])
    @info "Quick" app minimum(node.quick_store[2:end]) maximum(node.quick_store[2:end])
    @info "Slow" app minimum(node.slow_store[2:end]) maximum(node.slow_store[2:end])

    # Save performance over validation period
    # valid_cmd = node.storage[CALIB_LN+1:end]
    valid_outflow = node.outflow[CALIB_LN+1:end]
    outflows[app] = valid_outflow

    rmse = round(Streamfall.RMSE(VALID_OBS, valid_outflow), digits=2)
    mkge = round(Streamfall.mKGE(VALID_OBS, valid_outflow), digits=2)
    nse = round(Streamfall.NSE(VALID_OBS, valid_outflow), digits=2)
    nnse = round(Streamfall.NNSE(VALID_OBS, valid_outflow), digits=2)
    nmkge = round(Streamfall.NmKGE(VALID_OBS, valid_outflow), digits=2)
    mean_nmkge = round(Streamfall.mean_NmKGE(VALID_OBS, valid_outflow), digits=2)
    @info "$app validation: "
    report_metrics(VALID_OBS, valid_outflow)

    plot(VALID_DATES, VALID_OBS, label="Historic", title="$display_name Validation\n(RMSE: $rmse; NSE: $nse; KGE': $mkge)", legend=:topright)
    plot!(VALID_DATES, valid_outflow, label="Simulated", alpha=0.7, ylabel="Stream flow [ML]")
    savefig(joinpath(baseline_path, "baseline_valid_perf_$(app).png"))

    qqplot(VALID_OBS, valid_outflow, title="Validation Period\n$display_name")
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(joinpath(baseline_path, "baseline_valid_perf_qq_$(app).png"))

    full_dates = vcat(CALIB_DATES, VALID_DATES)
    plot(full_dates[calib_start:CALIB_LN], CALIB_OBS[calib_start:end], label="Calibration", title="$display_name C&V", legend=:topright)
    plot!(full_dates[CALIB_LN+1:end], VALID_OBS, label="Validation", color="green")
    plot!(full_dates[calib_start:end], node.outflow[calib_start:end], label="Simulated", linestyle=:dashdot, color="red", alpha=0.7, ylabel="Stream flow [ML]")
    savefig(joinpath(baseline_path, "baseline_full_perf_$(app).png"))
end


# Ensemble results
ensemble_outflow = mean.(outflows["NNSE"] .+ outflows["mean_NmKGE"])
rmse = round(Streamfall.RMSE(VALID_OBS, ensemble_outflow), digits=2)
mkge = round(Streamfall.mKGE(VALID_OBS, ensemble_outflow), digits=2)
nse = round(Streamfall.NSE(VALID_OBS, ensemble_outflow), digits=2)

plot(VALID_DATES, VALID_OBS, label="Historic", title="Baseline Ensemble\n(RMSE: $rmse; NSE: $nse; KGE': $mkge)", legend=:topright)
plot!(VALID_DATES, ensemble_outflow, label="Ensemble", alpha=0.8)
savefig(joinpath(FIG_PATH, "baseline_valid_perf_ensemble_NNSE_meanNmKGE.png"))




# keeping the below as it plots out simulation for each bin
# climate = Climate(training, "_rain", "_temp")

# # Set up stream network
# network = YAML.load_file(joinpath(DATA_PATH, "cotter_baseline_calibrated_NNSE.yml"))
# nnse_sn = create_network("Calibrated NNSE", network)

# run_basin!(nnse_sn, climate)

# node = nnse_sn[1]
# nnse_outflow = node.outflow
# nnse_cmd = node.storage[2:end]  # skip the first entry as it is the initial condition

# normal_dist = fit(Normal, nnse_cmd)
# trunc_dist = truncated(Normal(mean(normal_dist), std(normal_dist)), 0.0, maximum(nnse_cmd))

# # save to file
# dist_data = Dict(
#     "model" => ["Baseline_NNSE", "Baseline_NmKGE", "Baseline_RMSE"],
#     "mean" => [mean(normal_dist), ],
#     "std" => [std(normal_dist), ],
#     "max" => [maximum(nnse_cmd), ],
# )

# nnse_bin_bounds = quantile.(trunc_dist, [0.5-0.25, 0.5+0.25, 1.0])
# histogram(nnse_cmd, alpha=0.4, normalize=true,
#           title="Distribution of CMD values\n(NNSE calibrated baseline)",
#           legend=:left,
#           label="CMD",
#           xlabel="Catchment Moisture Deficit",
#           ylabel="Probability")
# StatsPlots.plot!(trunc_dist, linewidth=2, label="Fitted normal distribution")
# vline!(nnse_bin_bounds; color="orange", linewidth=1.5, label="Ensemble bounds")

# savefig(joinpath(FIG_PATH, "NNSE_calibrated_dist.png"))

# ######

# # Set up stream network
# nmkge_network = YAML.load_file(joinpath(DATA_PATH, "cotter_baseline_calibrated_NmKGE.yml"))
# nmkge_sn = create_network("Calibrated NmKGE", nmkge_network)


# run_basin!(nmkge_sn, climate)

# node = nmkge_sn[1]
# nmkge_outflow = node.outflow
# nmkge_cmd = node.storage[2:end]

# normal_dist = fit(Normal, nmkge_cmd)
# trunc_dist = truncated(Normal(mean(normal_dist), std(normal_dist)), 0.0, maximum(nmkge_cmd))

# # Save data
# append!(dist_data["mean"], mean(normal_dist))
# append!(dist_data["std"], std(normal_dist))
# append!(dist_data["max"], maximum(nmkge_cmd))

# nmkge_bin_bounds = quantile.(trunc_dist, [0.5-0.25, 0.5+0.25, 1.0])
# histogram(nmkge_cmd, alpha=0.4, normalize=true,
#           title="Distribution of CMD values\n(NmKGE calibrated baseline)",
#           legend=:left,
#           label="CMD",
#           xlabel="Catchment Moisture Deficit",
#           ylabel="Probability")
# StatsPlots.plot!(trunc_dist, linewidth=2, label="Fitted normal distribution")
# vline!(nmkge_bin_bounds; color="orange", linewidth=1.5, label="Ensemble bounds")

# savefig(joinpath(FIG_PATH, "NmKGE_calibrated_dist.png"))


# ######

# # Set up stream network
# rmse_network = YAML.load_file(joinpath(DATA_PATH, "cotter_baseline_calibrated_RMSE.yml"))
# rmse_sn = create_network("Calibrated RMSE", rmse_network)

# run_basin!(rmse_sn, climate)

# node = rmse_sn[1]
# rmse_outflow = node.outflow
# rmse_cmd = node.storage[2:end]

# normal_dist = fit(Normal, rmse_cmd)
# trunc_dist = truncated(Normal(mean(normal_dist), std(normal_dist)), 0.0, maximum(rmse_cmd))

# # Save data
# append!(dist_data["mean"], mean(normal_dist))
# append!(dist_data["std"], std(normal_dist))
# append!(dist_data["max"], maximum(rmse_cmd))
# CSV.write(joinpath(DATA_PATH, "dist_data.csv"), DataFrame(dist_data))

# nmkge_bin_bounds = quantile.(trunc_dist, [0.5-0.25, 0.5+0.25, 1.0])
# histogram(rmse_cmd, alpha=0.4, normalize=true,
#           title="Distribution of CMD values\n(RMSE calibrated baseline)",
#           legend=:left,
#           label="CMD",
#           xlabel="Catchment Moisture Deficit",
#           ylabel="Probability")
# StatsPlots.plot!(trunc_dist, linewidth=2, label="Fitted normal distribution")
# vline!(nmkge_bin_bounds; color="orange", linewidth=1.5, label="Ensemble bounds")

# savefig(joinpath(FIG_PATH, "RMSE_calibrated_dist.png"))

# #######

# # Performance of baseline and ensemble of baselines
# hist_flow = training_flow["410730"]
# mean_outflow = mean([nnse_outflow, nmkge_outflow], dims=1)[1]
# @info "Performance of baseline and ensemble of baselines"
# @info "NNSE Calibrated" Streamfall.RMSE(hist_flow, nnse_outflow)
# @info "NmKGE Calibrated" Streamfall.RMSE(hist_flow, nmkge_outflow)
# @info "RMSE Calibrated" Streamfall.RMSE(hist_flow, rmse_outflow)
# @info "Mean of Ensemble" Streamfall.RMSE(hist_flow, mean_outflow)


# # Performance within each bin
# nnse_bin_membership = bin_membership(nnse_cmd, nnse_bin_bounds)
# @info "Performance within bins (NNSE)"
# for idx in 1:length(nnse_bin_bounds)+1
#     bin_member = findall(x -> x == idx, nnse_bin_membership)
#     @info "Bin $idx ($(length(bin_member)))"
#     @info "NNSE" Streamfall.RMSE(hist_flow[bin_member], nnse_outflow[bin_member])
#     @info "Ensemble" Streamfall.RMSE(hist_flow[bin_member], mean_outflow[bin_member])
# end

# ######
# println("------------")



# nmkge_bin_membership = bin_membership(nmkge_cmd, nmkge_bin_bounds)
# @info "Performance within bins (NmKGE)"
# for idx in 1:length(nmkge_bin_bounds)+1
#     bin_member = findall(x -> x == idx, nmkge_bin_membership)
#     @info "Bin $idx ($(length(bin_member)))"
#     if length(bin_member) == 0
#         @info "No members"
#         continue
#     end
#     @info "NmKGE" Streamfall.RMSE(hist_flow[bin_member], nmkge_outflow[bin_member])
#     @info "Ensemble" Streamfall.RMSE(hist_flow[bin_member], mean_outflow[bin_member])
# end

# plot(hist_flow, title="Baseline results", label="Observed",
#      xlabel="Time step", ylabel="Streamflow [ML]")
# num_bounds = length(nnse_bin_bounds)+1
# for idx in 1:num_bounds
#     bin_member = findall(x -> x == idx, nnse_bin_membership)
#     if length(bin_member) == 0
#         continue
#     end
#     scatter!(bin_member, nnse_outflow[bin_member],
#              markersize=2, label="Bin $idx",
#              markerstrokewidth=0,
#              alpha=0.6;
#              palette=cgrad(:matter, num_bounds, categorical = true))
# end

# savefig(joinpath(FIG_PATH, "NNSE_baseline_perf.png"))

# plot(hist_flow, title="Baseline results", label="Observed",
#      xlabel="Time step", ylabel="Streamflow [ML]")
# for idx in 1:num_bounds
#     bin_member = findall(x -> x == idx, nmkge_bin_membership)
#     if length(bin_member) == 0
#         continue
#     end
#     scatter!(bin_member, nmkge_outflow[bin_member],
#              markersize=2, label="Bin $idx",
#              markerstrokewidth=0,
#              alpha=0.6;
#              palette=cgrad(:matter, num_bounds, categorical = true))
# end

# savefig(joinpath(FIG_PATH, "NmKGE_baseline_perf.png"))


# v_flow = validation_flow["410730"]
# reset!(nnse_sn)
# v_climate = Climate(validation, "_rain", "_temp")
# run_basin!(nnse_sn, v_climate)


# plot(v_flow, label="Observed", alpha=0.5, title="Validation Period")
# nnse_outflow = nnse_sn[1].outflow
# rmse = round(Streamfall.RMSE(v_flow, nnse_outflow), digits=4)
# plot!(nnse_outflow, label="NNSE calibrated: RMSE = $rmse", linestyle=:dash, alpha=0.9)


# ###
# reset!(nmkge_sn)
# run_basin!(nmkge_sn, v_climate)
# nmkge_outflow = nmkge_sn[1].outflow
# rmse = round(Streamfall.RMSE(v_flow, nmkge_outflow), digits=4)
# plot!(nmkge_outflow, label="NmKGE calibrated: RMSE = $rmse", linestyle=:dashdotdot, alpha=0.8)

# ####
# reset!(rmse_sn)
# run_basin!(rmse_sn, v_climate)
# rmse_outflow = rmse_sn[1].outflow
# rmse = round(Streamfall.RMSE(v_flow, rmse_outflow), digits=4)
# plot!(rmse_outflow, label="RMSE calibrated: RMSE = $rmse", linestyle=:dashdotdot, alpha=0.7)

# # mean_outflow = mean([nnse_outflow, nmkge_outflow], dims=1)[1]
# # rmse = round(Streamfall.RMSE(v_flow, mean_outflow), digits=4)
# # plot!(mean_outflow, label="Ensemble mean: RMSE = $rmse")

# savefig(joinpath(FIG_PATH, "calibrated_performance.png"))




