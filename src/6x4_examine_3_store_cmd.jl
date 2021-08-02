using Streamfall: load_calibration
include("_common.jl")
include("_calib_funcs.jl")

using StatsPlots

ensemble_figs = joinpath(FIG_PATH, "ensemble_results", "cmd") * "/"
mkpath(ensemble_figs)


function read_params(fn)
    calibrated_params = open(fn, "r") do f
        p_ins = Float64[]
        for p in eachline(f)
            p_ins = collect(map(x -> parse(Float64, x), split(p, ",")))
        end

        return p_ins
    end

    return calibrated_params
end


state = :storage
state_str = String(state)
# state_approach = "NNSE"
# state_approach = "RMSE"
# state_approach = "split_NnpKGE"
# state_approach = "split_RMSE"
# state_approach = "split_NmKGE"
# state_approach = "NmKGE"
# state_approach = "split_NNSE"
state_approach = "split_mean_NmKGE"

baseline_approach = "split_NnpKGE"
# baseline_approach = "split_NNSE"
# baseline_approach = "split_NmKGE"
# baseline_approach = "split_mean_NmKGE"
# baseline_approach = "NmKGE"

target_idx = [1,2,5,6,7,8]

baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"

sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$state_approach.yml")
node = sn[n_id]

# Read in parameters
fn = "$(calib_param_path)restricted_log_$(state_str)_$(state_approach)_10_year.txt"
calibrated_params = read_params(fn)

reset!(sn)
calib_ts = 3651 # 1826
_, x0, __ = param_info(node; with_level=false)
mod_params = update_partial(x0, target_idx, calibrated_params)
quantiles = [0.0, 0.1, 0.90]
active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; burn_in=calib_ts, log_transform=true, releases=nothing)

base_sn, base_id = setup_network("baselines/cotter_baseline_IHACRES_$(baseline_approach).yml")
base_node = base_sn[base_id]
Streamfall.run_basin!(base_sn, FULL_CLIMATE)

obs_flow = FULL_DATASET[calib_ts:CALIB_LN, "410730_Q"]
base_flow = base_node.outflow[calib_ts:CALIB_LN]
sim_flow = node.outflow[calib_ts:CALIB_LN]


ensemble_flow = copy(base_flow)
low_store = findall(active_param_set[calib_ts:CALIB_LN] .== 1)
mid_store = findall(active_param_set[calib_ts:CALIB_LN] .== 2)
high_store = findall(active_param_set[calib_ts:CALIB_LN] .== 3)

@assert length(mid_store) > 0

high_base = Streamfall.NmKGE(obs_flow[high_store], base_flow[high_store])
high_state = Streamfall.NmKGE(obs_flow[high_store], sim_flow[high_store])
high_factor = high_state >= (high_base * 1.0) ? 0.99 : 0.01

mid_base = Streamfall.NmKGE(obs_flow[mid_store], base_flow[mid_store])
mid_state = Streamfall.NmKGE(obs_flow[mid_store], sim_flow[mid_store])
mid_factor = mid_state >= (mid_base * 1.0) ? 0.8 : 0.2

low_base = Streamfall.NmKGE(obs_flow[low_store], base_flow[low_store])
low_state = Streamfall.NmKGE(obs_flow[low_store], sim_flow[low_store])
low_factor = low_state >= (low_base * 1.0) ? 1.0 : 0.0

# ensemble_flow[high_store] = ((ensemble_flow[high_store] .* high_factor) .+ (sim_flow[high_store] .* (1.0 - high_factor)))
# ensemble_flow[mid_store] = ((ensemble_flow[mid_store] .* mid_factor) .+ (sim_flow[mid_store] .* (1.0 - mid_factor)))
# ensemble_flow[low_store] = ((ensemble_flow[low_store] .* low_factor) .+ (sim_flow[low_store] .* (1.0 - low_factor)))
ensemble_flow[high_store] = ((sim_flow[high_store] .* high_factor) .+ (ensemble_flow[high_store] .* (1.0 - high_factor)))
ensemble_flow[mid_store] = ((sim_flow[mid_store] .* mid_factor) .+ (ensemble_flow[mid_store] .* (1.0 - mid_factor)))
ensemble_flow[low_store] = ((sim_flow[low_store] .* low_factor) .+ (ensemble_flow[low_store] .* (1.0 - low_factor)))

@info "Calibration metrics for $(baseline_approach) [High/Mid/Low]"
report_metrics(obs_flow[high_store], base_flow[high_store])
report_metrics(obs_flow[mid_store], base_flow[mid_store])
report_metrics(obs_flow[low_store], base_flow[low_store])


@info "Simulation metrics for $(state_approach) [Overall/High/Mid/Low]"
report_metrics(obs_flow, sim_flow)
report_metrics(obs_flow[high_store], sim_flow[high_store])
report_metrics(obs_flow[mid_store], sim_flow[mid_store])
report_metrics(obs_flow[low_store], sim_flow[low_store])

@info "Metrics for Ensemble [High/Mid/Low]"
report_metrics(obs_flow, ensemble_flow)
report_metrics(obs_flow[high_store], ensemble_flow[high_store])
report_metrics(obs_flow[mid_store], ensemble_flow[mid_store])
report_metrics(obs_flow[low_store], ensemble_flow[low_store])


# sanity_check
@info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], sim_flow[mid_store]) Streamfall.RMSE(obs_flow[low_store], sim_flow[low_store])
@info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)


plot_dates = FULL_DATASET.Date[calib_ts:CALIB_LN]
baseline_plot = begin
    plot(plot_dates, obs_flow, title="Baseline\n($baseline_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end

ensemble_plot = begin
    plot(plot_dates, obs_flow, title="Ensemble\n($state_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end
savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)flow_cmd_$(baseline_approach)_$(state_approach)_calib.png")

base_store = base_node.gw_store[calib_ts:CALIB_LN+1]
pos_markers = quantile(base_store, quantiles[2:end])

base_low_store = findall(base_store .< pos_markers[1])
base_mid_store = findall(pos_markers[1] .>= base_store .< pos_markers[2])
base_high_store = findall(base_store .>= pos_markers[2])

score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)
high = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
mid = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
low = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
base_qq = qqplot(obs_flow, base_flow,
                 title="Baseline $(baseline_approach)\n(KGE': $score; Wet: $low Mid: $mid; Dry: $high)",
                 titlefontsize=8,
                 markerstrokewidth=0,
                 markersize=2,
                 alpha=0.5,
                 )

score = round(Streamfall.mKGE(obs_flow, sim_flow), digits=3)
high = round(Streamfall.mKGE(obs_flow[high_store], sim_flow[high_store]), digits=3)
mid = round(Streamfall.mKGE(obs_flow[mid_store], sim_flow[mid_store]), digits=3)
low = round(Streamfall.mKGE(obs_flow[low_store], sim_flow[low_store]), digits=3)
sim_qq = qqplot(obs_flow, sim_flow,
                title="State-based $state_approach\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2,
                alpha=0.5,
                )

sim_qq_col = begin
    # plot(
    #     sim_qq,
    #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
    #     dpi=300
    # )

    low_plot = begin
        scatter(obs_flow[low_store], ensemble_flow[low_store],
                title="State-based Wet\n(KGE': $(low))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    mid_plot = begin
        scatter(obs_flow[mid_store], ensemble_flow[mid_store],
                title="State-based Mid\n(KGE': $(mid))",
                titlefontsize=8,
                markersize=2.0,
                markerstrokewidth=0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="State-based Dry\n(KGE': $(high))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    plot(
        sim_qq,
        low_plot,
        mid_plot,
        high_plot,
        dpi=300
    )
end

savefig(sim_qq_col, "$(ensemble_figs)sim_cmd_$(baseline_approach)_$(state_approach)_qq_calib.png")


score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=3)
high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=3)
mid = round(Streamfall.mKGE(obs_flow[mid_store], ensemble_flow[mid_store]), digits=3)
low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=3)

ensemble_qq = qqplot(obs_flow, ensemble_flow,
                     title="$(baseline_approach) - $(state_approach)\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
                     titlefontsize=8,
                     markerstrokewidth=0,
                     alpha=0.5,
                     markersize=2)

ensemble_col = begin
    # plot(
    #     ensemble_qq,
    #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
    #     dpi=300
    # )
    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="Dry Ensemble\n(KGE': $(high))",
                titlefontsize=8,
                markersize=2.0,
                markerstrokewidth=0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    mid_plot = begin
        scatter(obs_flow[mid_store], ensemble_flow[mid_store],
                title="Mid Ensemble\n(KGE': $(mid))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    low_plot = begin
        scatter(obs_flow[low_store], ensemble_flow[low_store],
                title="Wet Ensemble\n(KGE': $(low))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    plot(
        ensemble_qq,
        low_plot,
        mid_plot,
        high_plot,
        dpi=300
    )
end
savefig(ensemble_col, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_qq_calib.png")


snapshot = begin
    plot(
        base_qq,
        sim_qq,
        ensemble_qq,
        ensemble_plot,
        dpi=300
    )
    plot!(size=(600,400))
end

savefig(snapshot, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_result_calib.png")



high_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
    plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
end

mid_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
    plot(plot_dates[mid_store], obs_flow[mid_store], title="Mid GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[mid_store], base_flow[mid_store], alpha=0.6)
end

low_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
    plot(plot_dates[low_store], obs_flow[low_store], title="Low GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[low_store], base_flow[low_store], alpha=0.6)
end

baseline_qq = begin
    overall_score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)
    high_score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
    mid_score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
    low_score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)

    overall = begin
        scatter(obs_flow, base_flow,
                title="Baseline $(baseline_approach)\n(KGE': $overall_score)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow, obs_flow, alpha=0.6, color="green")
    end

    high_plot = begin
        scatter(obs_flow[high_store], base_flow[high_store],
                title="Dry $(baseline_approach)\n(KGE': $(high_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], alpha=0.6, color="green")
    end

    mid_plot = begin
        scatter(obs_flow[mid_store], base_flow[mid_store],
                title="Mid $(baseline_approach)\n(KGE': $(mid_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[mid_store], obs_flow[mid_store], alpha=0.6, color="green")
    end

    low_plot = begin
        scatter(obs_flow[low_store], base_flow[low_store],
                title="Wet $(baseline_approach)\n(KGE': $(low_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], alpha=0.6, color="green")
    end
    # plot(
    #     qqplot(obs_flow, base_flow, title="Baseline $(baseline_approach)\n(KGE': $overall_score)", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[low_store], base_flow[low_store], title="Wet $(baseline_approach)\n(KGE': $(low_score))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[mid_store], base_flow[mid_store], title="Mid $(baseline_approach)\n(KGE': $(mid_score))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[high_store], base_flow[high_store], title="Dry $(baseline_approach)\n(KGE': $(high_score))", titlefontsize=8, markersize=2.0),
    # )
    plot(
        overall,
        low_plot,
        mid_plot,
        high_plot,
        dpi=300
    )
    plot!(size=(600,400))
end
savefig(baseline_qq, joinpath(ensemble_figs, "baseline_cmd_$(baseline_approach)_$(state_approach)_qq_calib.png"))



### Validation ###

obs_flow = FULL_DATASET[CALIB_LN+1:end, "410730_Q"]
base_flow = base_node.outflow[CALIB_LN+1:end]
sim_flow = node.outflow[CALIB_LN+1:end]


ensemble_flow = copy(base_flow)
low_store = findall(active_param_set[CALIB_LN+1:end] .== 1)
mid_store = findall(active_param_set[CALIB_LN+1:end] .== 2)
high_store = findall(active_param_set[CALIB_LN+1:end] .== 3)

@assert length(mid_store) > 0

ensemble_flow[high_store] = ((sim_flow[high_store] .* high_factor) .+ (ensemble_flow[high_store] .* (1.0 - high_factor)))
ensemble_flow[mid_store] = ((sim_flow[mid_store] .* mid_factor) .+ (ensemble_flow[mid_store] .* (1.0 - mid_factor)))

if length(low_store) > 0
    ensemble_flow[low_store] = ((sim_flow[low_store] .* low_factor) .+ (ensemble_flow[low_store] .* (1.0 - low_factor)))
end

# sanity_check
# @info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], sim_flow[mid_store])
# @info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)
# sanity_check
@info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], sim_flow[mid_store]) Streamfall.RMSE(obs_flow[low_store], sim_flow[low_store])
@info "Bin Metric results" bin_metric(active_param_set[CALIB_LN+1:end], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)

plot_dates = FULL_DATASET.Date[CALIB_LN+1:end]
baseline_plot = begin
    plot(plot_dates, obs_flow, title="Baseline\n($baseline_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end

ensemble_plot = begin
    plot(plot_dates, obs_flow, title="Ensemble\n($state_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end
savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)flow_cmd_$(baseline_approach)_$(state_approach)_valid.png")

base_store = base_node.gw_store[CALIB_LN+1:end]
pos_markers = quantile(base_store, quantiles[2:end])

base_low_store = findall(base_store .< pos_markers[1])
base_mid_store = findall(pos_markers[1] .>= base_store .< pos_markers[2])
base_high_store = findall(base_store .>= pos_markers[2])

score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)

high = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)

if length(mid_store) > 0
    mid = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
else
    mid = NaN
end


if length(low_store) > 0
    low = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)    
else
    low = NaN
end

base_qq = qqplot(obs_flow, base_flow,
                 title="Baseline $(baseline_approach)\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
                 titlefontsize=8,
                 markerstrokewidth=0,
                 markersize=2,
                 alpha=0.5,
                 )

score = round(Streamfall.mKGE(obs_flow, sim_flow), digits=3)

high = round(Streamfall.mKGE(obs_flow[high_store], sim_flow[high_store]), digits=3)

if length(mid_store) > 0
    mid = round(Streamfall.mKGE(obs_flow[mid_store], sim_flow[mid_store]), digits=3)
else
    mid = NaN
end


if length(low_store) > 0
    low = round(Streamfall.mKGE(obs_flow[low_store], sim_flow[low_store]), digits=3)
else
    low = NaN
end

sim_qq = qqplot(obs_flow, sim_flow,
                title="State-based $state_approach\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2,
                alpha=0.5,
                )

sim_qq_col = begin
    # plot(
    #     sim_qq,
    #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
    #     dpi=300
    # )
    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="State-based Dry\n(KGE': $(high))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    mid_plot = begin
        scatter(obs_flow[mid_store], ensemble_flow[mid_store],
                title="State-based Mid\n(KGE': $(mid))",
                titlefontsize=8,
                markersize=2.0,
                markerstrokewidth=0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    low_plot = begin
        scatter(obs_flow[low_store], ensemble_flow[low_store],
                title="State-based Wet\n(KGE': $(low))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    plot(
        sim_qq,
        low_plot,
        mid_plot,
        high_plot,
        dpi=300
    )
end

savefig(sim_qq_col, "$(ensemble_figs)sim_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png")


score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=3)

high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=3)

if length(mid_store) > 0
    mid = round(Streamfall.mKGE(obs_flow[mid_store], ensemble_flow[mid_store]), digits=3)
else
    mid = NaN
end

if length(low_store) > 0
    low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=3)
else
    low = NaN
end

ensemble_qq = qqplot(obs_flow, ensemble_flow,
                     title="$(baseline_approach) - $(state_approach)\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
                     titlefontsize=8,
                     markerstrokewidth=0,
                     alpha=0.5,
                     markersize=2)

ensemble_col = begin
    # plot(
    #     ensemble_qq,
    #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
    #     dpi=300
    # )
    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="Dry Ensemble\n(KGE': $(high))",
                titlefontsize=8,
                markersize=2.0,
                markerstrokewidth=0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    mid_plot = begin
        scatter(obs_flow[mid_store], ensemble_flow[mid_store],
                title="Mid Ensemble\n(KGE': $(mid))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    low_plot = begin
        scatter(obs_flow[low_store], ensemble_flow[low_store],
                title="Wet Ensemble\n(KGE': $(low))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    plot(
        ensemble_qq,
        low_plot,
        mid_plot,
        high_plot,
        dpi=300
    )
end
savefig(ensemble_col, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png")


snapshot = begin
    plot(
        base_qq,
        sim_qq,
        ensemble_qq,
        ensemble_plot,
        dpi=300
    )
    plot!(size=(600,400))
end

savefig(snapshot, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_result_valid.png")



high_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
    plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
end

mid_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
    plot(plot_dates[mid_store], obs_flow[mid_store], title="Mid GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[mid_store], base_flow[mid_store], alpha=0.6)
end

low_store_plot = begin
    if length(low_store) > 0
        score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
    else
        score = NaN
    end
    plot(plot_dates[low_store], obs_flow[low_store], title="Low GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[low_store], base_flow[low_store], alpha=0.6)
end

baseline_qq = begin
    overall_score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)
    high_score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
    mid_score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)

    if length(low_store) > 0
        low_score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
    else
        low_score = NaN
    end

    overall = begin
        scatter(obs_flow, base_flow,
                title="Baseline $(baseline_approach)\n(KGE': $overall_score)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow, obs_flow, alpha=0.6, color="green")
    end

    high_plot = begin
        scatter(obs_flow[high_store], base_flow[high_store],
                title="Dry $(baseline_approach)\n(KGE': $(high_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], alpha=0.6, color="green")
    end

    mid_plot = begin
        scatter(obs_flow[mid_store], base_flow[mid_store],
                title="Mid $(baseline_approach)\n(KGE': $(mid_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[mid_store], obs_flow[mid_store], alpha=0.6, color="green")
    end

    low_plot = begin
        scatter(obs_flow[low_store], base_flow[low_store],
                title="Wet $(baseline_approach)\n(KGE': $(low_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], alpha=0.6, color="green")
    end
    # plot(
    #     qqplot(obs_flow, base_flow, title="Baseline $(baseline_approach)\n(KGE': $overall_score)", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[low_store], base_flow[low_store], title="Wet $(baseline_approach)\n(KGE': $(low_score))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[mid_store], base_flow[mid_store], title="Mid $(baseline_approach)\n(KGE': $(mid_score))", titlefontsize=8, markersize=2.0),
    #     qqplot(obs_flow[high_store], base_flow[high_store], title="Dry $(baseline_approach)\n(KGE': $(high_score))", titlefontsize=8, markersize=2.0),
    # )
    plot(
        overall,
        low_plot,
        mid_plot,
        high_plot,
        dpi=300
    )
    plot!(size=(600,400))
end
savefig(baseline_qq, joinpath(ensemble_figs, "baseline_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png"))
