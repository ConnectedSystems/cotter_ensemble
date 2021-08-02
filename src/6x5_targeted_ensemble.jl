using Streamfall: load_calibration
include("_common.jl")
include("_calib_funcs.jl")

using StatsPlots

ensemble_figs = joinpath(FIG_PATH, "ensemble_results", "targeted") * "/"
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


state = :gw_store
state_str = String(state)
# state_approach = "NNSE"
# state_approach = "RMSE"
state_approach = "split_NnpKGE"
# state_approach = "split_NmKGE"
# state_approach = "split_NNSE"
# state_approach = "split_mean_NmKGE"
target_idx = [5,6,7,8]

baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"

top_sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$state_approach.yml")
top_node = sn[n_id]

# Read in parameters
fn = "$(calib_param_path)target_top-end_log_$(state_str)_$(state_approach).txt"
calibrated_params = read_params(fn)


# Run state-based model for top
reset!(top_sn)
calib_ts = 1826
_, x0, __ = param_info(top_node; with_level=false)
mod_params = update_partial(x0, target_idx, calibrated_params)
quantiles = [0.0, 0.99]
high_active, param_idxs = online_burn_state_node!(top_sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; log_transform=true, releases=nothing)

# Run state-based model for bottom
fn = "$(calib_param_path)calib_params_online_target_low-end_$(state_str)_$(state_approach).txt"
calibrated_params = read_params(fn)
low_sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$state_approach.yml")
low_node = low_sn[n_id]

_, x0, __ = param_info(low_node; with_level=false)
mod_params = update_partial(x0, target_idx, calibrated_params)
quantiles = [0.0, 0.1]
low_active, param_idxs = online_burn_state_node!(low_sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; log_transform=true, releases=nothing)


# Run baseline
baseline_approach = "split_NnpKGE"
# baseline_approach = "split_NNSE"
base_sn, base_id = setup_network("baselines/cotter_baseline_IHACRES_$(baseline_approach).yml")
base_node = base_sn[base_id]
Streamfall.run_basin!(base_sn, FULL_CLIMATE)

obs_flow = FULL_DATASET[calib_ts:CALIB_LN, "410730_Q"]
base_flow = base_node.outflow[calib_ts:CALIB_LN]
low_flow = low_node.outflow[calib_ts:CALIB_LN]
top_flow = top_node.outflow[calib_ts:CALIB_LN]


ensemble_flow = copy(base_flow)
low_store = findall(low_active[calib_ts:CALIB_LN] .== 1)
# mid_store = findall(active_param_set[calib_ts:CALIB_LN] .== 2)
high_store = findall(high_active[calib_ts:CALIB_LN] .== 2)

high_base = Streamfall.NmKGE(obs_flow[high_store], base_flow[high_store])
high_top = Streamfall.NmKGE(obs_flow[high_store], top_flow[high_store])
high_factor = high_top >= (high_base * 1.0) ? 0.99 : 0.01

low_base = Streamfall.NmKGE(obs_flow[low_store], base_flow[low_store])
low_state = Streamfall.NmKGE(obs_flow[low_store], low_flow[low_store])
low_factor = low_state >= (low_base * 1.0) ? 0.99 : 0.01

ensemble_flow[high_store] = ((top_flow[high_store] .* high_factor) .+ (ensemble_flow[high_store] .* (1.0 - high_factor)))
ensemble_flow[low_store] = ((low_flow[low_store] .* low_factor) .+ (ensemble_flow[low_store] .* (1.0 - low_factor)))

# sanity_check
@info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], ensemble_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], ensemble_flow[mid_store])
@info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, ensemble_flow, split_NNSE, nothing)


plot_dates = FULL_DATASET.Date[calib_ts:CALIB_LN]
baseline_plot = begin
    plot(plot_dates, obs_flow, title="Baseline\n($baseline_approach)", titlefontsize=8, legendfontsize=8, linewidth=2, legend=false)
    plot!(plot_dates, base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end

ensemble_plot = begin
    plot(plot_dates, obs_flow, title="Ensemble\n($state_approach)", titlefontsize=8, legendfontsize=8, linewidth=2, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end
savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)flow_$(baseline_approach)_$(state_approach).png")

base_store = base_node.gw_store[calib_ts:CALIB_LN+1]
pos_markers = quantile(base_store, quantiles[1:end])

base_low_store = findall(base_store .< pos_markers[1])
base_high_store = findall(base_store .>= pos_markers[2])

score = round(Streamfall.mKGE(obs_flow, base_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=2)
low = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=2)
base_qq = qqplot(obs_flow, base_flow,
                 title="Baseline $(baseline_approach)\n(KGE': $score; High: $high; Low: $low)",
                 titlefontsize=8,
                 markerstrokewidth=0,
                 markersize=2,
                 alpha=0.5,
                 )

score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=2)
low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=2)
sim_qq = qqplot(obs_flow, ensemble_flow,
                title="State-based $state_approach\n(KGE': $score; High: $high; Mid: $mid; Low: $low)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2,
                alpha=0.5,
                )

sim_qq_col = begin
    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="State-based High-store\n(KGE': $(high))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    low_plot = begin
        scatter(obs_flow[low_store], ensemble_flow[low_store],
                title="State-based Low-store\n(KGE': $(low))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    plot(
        sim_qq,
        high_plot,
        low_plot,
        dpi=300
    )
end

savefig(sim_qq_col, "$(ensemble_figs)sim_$(baseline_approach)_$(state_approach)_qq.png")


score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=2)
low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=2)

ensemble_qq = qqplot(obs_flow, ensemble_flow,
                     title="Ensemble\n(KGE': $score; High: $high; Low: $low)",
                     titlefontsize=8,
                     markerstrokewidth=0,
                     alpha=0.5,
                     markersize=2)

ensemble_col = begin
    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="High-store Ensemble\n(KGE': $(high))",
                titlefontsize=8,
                markersize=2.0,
                markerstrokewidth=0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
    end

    low_plot = begin
        scatter(obs_flow[low_store], ensemble_flow[low_store],
                title="Low-store Ensemble\n(KGE': $(low))",
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
        high_plot,
        dpi=300
    )
end
savefig(ensemble_col, "$(ensemble_figs)ensemble_$(baseline_approach)_$(state_approach)_qq.png")


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

savefig(snapshot, "$(ensemble_figs)ensemble_$(baseline_approach)_$(state_approach)_result.png")



high_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=2)
    plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
end

low_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=2)
    plot(plot_dates[low_store], obs_flow[low_store], title="Low GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[low_store], base_flow[low_store], alpha=0.6)
end

baseline_qq = begin
    overall_score = round(Streamfall.mKGE(obs_flow, base_flow), digits=2)
    high_score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=2)
    mid_score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=2)
    low_score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=2)

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
                title="High-store $(baseline_approach)\n(KGE': $(high_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], alpha=0.6, color="green")
    end

    low_plot = begin
        scatter(obs_flow[low_store], base_flow[low_store],
                title="Low-store $(baseline_approach)\n(KGE': $(low_score))",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2.0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[low_store], obs_flow[low_store], alpha=0.6, color="green")
    end

    plot(
        overall,
        high_plot,
        low_plot,
        dpi=300
    )
    plot!(size=(600,400))
end
savefig(baseline_qq, joinpath(ensemble_figs, "baseline_$(baseline_approach)_$(state_approach)_qq.png"))

### Validation ###

obs_flow = FULL_DATASET[CALIB_LN+1:end, "410730_Q"]
base_flow = base_node.outflow[CALIB_LN+1:end]
ensemble_flow = node.outflow[CALIB_LN+1:end]

ensemble_flow = copy(base_flow)
low_store = findall(active_param_set[CALIB_LN+1:end] .== 1)
high_store = findall(active_param_set[CALIB_LN+1:end] .== 2)

ensemble_flow[high_store] = ((ensemble_flow[high_store] .* high_factor) .+ (ensemble_flow[high_store] .* (1.0 - high_factor)))

if length(low_store) > 0
    ensemble_flow[low_store] = ((ensemble_flow[low_store] .* low_factor) .+ (ensemble_flow[low_store] .* (1.0 - low_factor)))
end

# sanity_check
# @info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], ensemble_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], ensemble_flow[mid_store])
# @info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, ensemble_flow, split_NNSE, nothing)


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
quantiles = [0.1, 0.99]
pos_markers = quantile(base_store, quantiles[1:end])

base_low_store = findall(base_store .< pos_markers[1])
base_high_store = findall(base_store .>= pos_markers[2])

score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)

high = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)

if length(low_store) > 0
    low = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)    
else
    low = NaN
end

base_qq = qqplot(obs_flow, base_flow,
                 title="Baseline $(baseline_approach)\n(KGE': $score; High: $high; Mid: $mid; Low: $low)",
                 titlefontsize=8,
                 markerstrokewidth=0,
                 markersize=2,
                 alpha=0.5,
                 )

score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=3)

high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=3)

if length(low_store) > 0
    low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=3)
else
    low = NaN
end

sim_qq = qqplot(obs_flow, ensemble_flow,
                title="State-based $state_approach\n(KGE': $score; High: $high; Mid: $mid; Low: $low)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2,
                alpha=0.5,
                )

sim_qq_col = begin
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
        high_plot,
        low_plot,
        dpi=300
    )
end

savefig(sim_qq_col, "$(ensemble_figs)sim_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png")


score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=3)

high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=3)

if length(low_store) > 0
    low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=3)
else
    low = NaN
end

ensemble_qq = qqplot(obs_flow, ensemble_flow,
                     title="$(baseline_approach) - $(state_approach)\n(KGE': $score; High: $high; Mid: $mid; Low: $low)",
                     titlefontsize=8,
                     markerstrokewidth=0,
                     alpha=0.5,
                     markersize=2)

ensemble_col = begin
    high_plot = begin
        scatter(obs_flow[high_store], ensemble_flow[high_store],
                title="High Ensemble\n(KGE': $(high))",
                titlefontsize=8,
                markersize=2.0,
                markerstrokewidth=0,
                alpha=0.5,
                legend=false)
        plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
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
        high_plot,
        low_plot,
        dpi=300
    )
end
savefig(ensemble_col, "$(ensemble_figs)ensemble_gw_store_$(baseline_approach)_$(state_approach)_qq_valid.png")


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

savefig(snapshot, "$(ensemble_figs)ensemble_gw_store_$(baseline_approach)_$(state_approach)_result_valid.png")



high_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
    plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false)
    plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
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
                title="High $(baseline_approach)\n(KGE': $(high_score))",
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

    plot(
        overall,
        high_plot,
        mid_plot,
        low_plot,
        dpi=300
    )
    plot!(size=(600,400))
end
savefig(baseline_qq, joinpath(ensemble_figs, "baseline_gw_store_$(baseline_approach)_$(state_approach)_qq_valid.png"))

