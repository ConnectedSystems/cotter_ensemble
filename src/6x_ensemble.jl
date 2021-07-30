using Streamfall: load_calibration
include("_common.jl")
include("_calib_funcs.jl")

using StatsPlots


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
# approach = "split_NnpKGE"
approach = "split_NNSE"
target_idx = [5,6,7,8]

baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"

sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$approach.yml")
node = sn[n_id]

# Read in parameters
# fn = "$(calib_param_path)calib_params_online_burn-in_multi_restricted_alt_params_$(state_str)_$(approach).txt"
# fn = "$(calib_param_path)calib_params_online_burn-in_top-end_$(state_str)_$(approach).txt"
# calib_params_online_burn-in_top-end_gw_store_split_mean_NmKGE
# calib_params_online_burn-in_top-end_gw_store_split_mean_NmKGE
fn = "$(calib_param_path)calib_params_online_burn-in_top-end_$(state_str)_$(approach).txt"
calibrated_params = read_params(fn)

reset!(sn)
calib_ts = 1826
_, x0, __ = param_info(node; with_level=false)
# calibrated_params = [4.1369754778087495,0.07784403998246077,3.879663457396416,0.2748380354619563, 4.1369754778087495,0.07784403998246077,3.879663457396416,0.2748380354619563,5.743645539293278,0.07237237294614068,9.171393274757758,0.398705852449268]
mod_params = update_partial(x0, target_idx, calibrated_params)
quantiles = [0.0, 0.99]
active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; log_transform=false, releases=nothing)

# plot(node.gw_store[2:end], color=active_param_set)


# baseline_approach = "split_NnpKGE"
baseline_approach = "NmKGE"
base_sn, base_id = setup_network("baselines/cotter_baseline_IHACRES_$(baseline_approach).yml")
base_node = base_sn[base_id]
Streamfall.run_basin!(base_sn, FULL_CLIMATE)

# obs_flow = CALIB_OBS[calib_ts:end]
# base_flow = base_node.outflow[calib_ts:CALIB_LN]
# sim_flow = node.outflow[calib_ts:CALIB_LN]
obs_flow = FULL_DATASET[calib_ts:CALIB_LN, "410730_Q"]  # CALIB_OBS[calib_ts:end]
base_flow = base_node.outflow[calib_ts:CALIB_LN]
sim_flow = node.outflow[calib_ts:CALIB_LN]


ensemble_flow = copy(base_flow)
high_flow = findall(active_param_set[calib_ts:CALIB_LN] .== 2)
ensemble_flow[high_flow] = ((ensemble_flow[high_flow] .* 0.1) .+ (sim_flow[high_flow] .* 0.9))

# @info "Calibration metrics for $(baseline_approach) [High/Low]"
# report_metrics(obs_flow[high_flow], base_flow[high_flow])
# report_metrics(obs_flow[Not(high_flow)], base_flow[Not(high_flow)])


@info "Baseline metrics for $(approach) [Overall/High/Low]"
report_metrics(obs_flow, sim_flow)
report_metrics(obs_flow[high_flow], sim_flow[high_flow])
report_metrics(obs_flow[Not(high_flow)], sim_flow[Not(high_flow)])

@info "Metrics for Ensemble [High/Low]"
report_metrics(obs_flow, ensemble_flow)
report_metrics(obs_flow[high_flow], ensemble_flow[high_flow])
report_metrics(obs_flow[Not(high_flow)], ensemble_flow[Not(high_flow)])

# sanity_check
@info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_flow], sim_flow[high_flow]) Streamfall.RMSE(obs_flow[Not(high_flow)], sim_flow[Not(high_flow)])
@info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)

# @info "Pbias?" Streamfall.PBIAS(obs_flow[high_flow], sim_flow[high_flow]) Streamfall.PBIAS(obs_flow[Not(high_flow)], sim_flow[Not(high_flow)])


plot_dates = FULL_DATASET.Date[calib_ts:CALIB_LN]
flow_plot = begin
    plot(plot_dates, obs_flow, title="Streamflow", titlefontsize=8, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.6)
end

base_high_flow_marker = quantile(base_node.gw_store[calib_ts:CALIB_LN], 0.99)
base_high_flow = findall(base_node.gw_store[calib_ts:CALIB_LN] .< base_high_flow_marker)
score = round(Streamfall.mKGE(obs_flow, base_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[base_high_flow], base_flow[base_high_flow]), digits=2)
low = round(Streamfall.mKGE(obs_flow[Not(base_high_flow)], base_flow[Not(base_high_flow)]), digits=2)
base_qq = qqplot(obs_flow, base_flow, 
                 title="Baseline $(baseline_approach)\n(KGE': $score; High Store: $high; Low: Store: $low)",
                 titlefontsize=8,
                 markerstrokewidth=0,
                 markersize=2,
                 )

score = round(Streamfall.mKGE(obs_flow, sim_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_flow], sim_flow[high_flow]), digits=2)
low = round(Streamfall.mKGE(obs_flow[Not(high_flow)], sim_flow[Not(high_flow)]), digits=2)
sim_qq = qqplot(obs_flow, sim_flow,
                title="State-based $approach\n(KGE': $score; High Store: $high; Low: Store: $low)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2,
                )

score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_flow], ensemble_flow[high_flow]), digits=2)
low = round(Streamfall.mKGE(obs_flow[Not(high_flow)], ensemble_flow[Not(high_flow)]), digits=2)
ensemble_qq = qqplot(obs_flow, ensemble_flow, 
                     title="Ensemble\n(KGE': $score; High Store: $high; Low: Store: $low)",
                     titlefontsize=8,
                     markerstrokewidth=0,
                     markersize=2)

savefig(plot(
    base_qq,
    sim_qq,
    ensemble_qq,
    flow_plot
), "$(FIG_PATH)ensemble_$(baseline_approach)_$(approach)_result.png")


high_flow_plot = begin
    score = round(Streamfall.mKGE(obs_flow[high_flow], base_flow[high_flow]), digits=2)
    plot(plot_dates[high_flow], obs_flow[high_flow], title="High GW Store\n(KGE': $score)", legend=false) 
    plot!(plot_dates[high_flow], base_flow[high_flow], alpha=0.6)
end

low_flow_plot = begin
    score = round(Streamfall.mKGE(obs_flow[Not(high_flow)], base_flow[Not(high_flow)]), digits=2)
    plot(plot_dates[Not(high_flow)], obs_flow[Not(high_flow)], title="Low GW Store\n(KGE': $score)", legend=false) 
    plot!(plot_dates[Not(high_flow)], base_flow[Not(high_flow)], alpha=0.6)
end

display(plot(
    qqplot(obs_flow[Not(high_flow)], base_flow[Not(high_flow)], title="Baseline $(baseline_approach)", markersize=2.0),
    high_flow_plot,
    low_flow_plot
))

# display(plot(
#     qqplot(obs_flow[high_flow], base_flow[high_flow], title="Baseline $(baseline_approach)\n[High GW Store]", xlabel="Historic", ylabel="Model"),
#     qqplot(obs_flow[Not(high_flow)], base_flow[Not(high_flow)], title="Baseline $(baseline_approach)\n[Low GW Store]"),
#     # qqplot(obs_flow, sim_flow, title="State-based"),
#     # qqplot(obs_flow, ensemble_flow, title="Ensemble"),
#     # flow_plot
# ))


for app in APPROACHES
    sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$app.yml")
    base_node = sn[n_id]

    # _, x0, _ = param_info(base_node; with_level=false)
    # @info "Parameters:" x0

    run_basin!(sn, FULL_CLIMATE)
    sim = base_node.outflow[calib_ts:CALIB_LN]

    gw_store = base_node.gw_store[calib_ts:CALIB_LN]
    top_end = quantile(gw_store, 0.99)
    low_flow = findall(sim .< top_end)

    @info "Overall/High/low GW $app:" Streamfall.RMSE(obs_flow, sim) Streamfall.RMSE(obs_flow[Not(low_flow)], sim[Not(low_flow)]) Streamfall.RMSE(obs_flow[low_flow], sim[low_flow])
    @info "Overall/High/low GW $app:" Streamfall.mKGE(obs_flow, sim) Streamfall.mKGE(obs_flow[Not(low_flow)], sim[Not(low_flow)]) Streamfall.mKGE(obs_flow[low_flow], sim[low_flow])
    @info "Overall/High/low GW $app:" Streamfall.NSE(obs_flow, sim) Streamfall.NSE(obs_flow[Not(low_flow)], sim[Not(low_flow)]) Streamfall.NSE(obs_flow[low_flow], sim[low_flow])
end





# outflows = Dict()

# target_idx = (5,6,7,8)
# for approach in APPROACHES
#     sn, n_id = setup_network("410730_IHACRES.yml")
#     node = sn[n_id]

#     # Read in parameters
#     infile = "$(calib_param_path)calib_params_online_burn-in_multi_rainfall_$(approach).txt"
#     calibrated_params = read_params(infile)

    # reset!(sn)
    # calib_ts = 1826
    # _, x0, __ = param_info(node; with_level=false)
    # mod_params = update_partial(x0, target_idx, calibrated_params)
    # quantiles = [0.95, 1.0]
    # active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; releases=nothing)

    # obs_flow = CALIB_OBS[calib_ts:end]
    # sim_flow = node.outflow[calib_ts:CALIB_LN]
    # @info "Calibration metrics for $(approach)"
    # report_metrics(obs_flow, sim_flow)


#     display_name = replace(approach, "_" => " ")
#     rmse = round(Streamfall.RMSE(obs_flow, sim_flow), digits=2)
#     mkge = round(Streamfall.mKGE(obs_flow, sim_flow), digits=2)
#     nse = round(Streamfall.NSE(obs_flow, sim_flow), digits=2)
#     plot(CALIB_DATES[calib_ts:end], obs_flow, title="$(display_name) Bin-calibrated Calibration\n(RMSE: $rmse; NSE: $nse; KGE': $mkge)", label="Calibration")
#     plot!(CALIB_DATES[calib_ts:end], sim_flow, label="Simulated", alpha=0.7)
#     splice!(quantiles, 2, 0.8)
#     # pushfirst!(quantiles, 0.6)
#     marks = quantile(sim_flow, quantiles)
#     hline!(marks, label="Thresholds")

#     savefig(joinpath(FIG_PATH, "rainfall_multi_calib_perf_$(approach).png"))

#     @info "Validation Metrics for $(approach)"
#     sim_flow = node.outflow[CALIB_LN+1:end]
#     report_metrics(VALID_OBS, sim_flow)

#     outflows[approach] = node.outflow

#     rmse = round(Streamfall.RMSE(VALID_OBS, sim_flow), digits=2)
#     mkge = round(Streamfall.mKGE(VALID_OBS, sim_flow), digits=2)
#     nse = round(Streamfall.NSE(VALID_OBS, sim_flow), digits=2)
#     plot(VALID_DATES, VALID_OBS, title="$(display_name) Bin-calibrated Validation\n(RMSE: $rmse; NSE: $nse; KGE': $mkge)", label="Validation")
#     plot!(VALID_DATES, sim_flow, label="Simulated", alpha=0.7)
#     savefig(joinpath(FIG_PATH, "rainfall_multi_valid_perf_$(approach).png"))

#     reset!(sn)
# end


# begin
#     # Selected based on goodness-of-fit to calibration data
#     selected_approaches = ["split_NNSE", "split_RMSE", "RMSE", "NNSE"]  # "RMSE", "mean_NmKGE", "split_RMSE", "split_NmKGE" "split_NNSE",

#     ensemble_flow = outflows[selected_approaches[1]]
#     for approach in selected_approaches[2:end]
#         app_flow = outflows[approach]
#         global ensemble_flow = ensemble_flow .+ app_flow
#     end

#     ensemble_flow = ensemble_flow ./ length(selected_approaches)

#     sim_flow = ensemble_flow[CALIB_LN+1:end]
#     rmse = round(Streamfall.RMSE(VALID_OBS, sim_flow), digits=2)
#     mkge = round(Streamfall.NmKGE(VALID_OBS, sim_flow), digits=2)
#     nse = round(Streamfall.NSE(VALID_OBS, sim_flow), digits=2)
#     plot(VALID_DATES, VALID_OBS, title="Ensemble Bin-calibrated Validation\n(RMSE: $rmse; NSE: $nse; KGE': $mkge)", label="Validation")
#     plot!(VALID_DATES, sim_flow, label="Simulated", alpha=0.8)
#     savefig(joinpath(FIG_PATH, "bin_state_ensemble_perf.png"))

#     # @info "Metrics for Ensemble:" ensemble_RMSE ensemble_NSE ensemble_mKGE
#     @info "Ensemble metrics:"
#     report_metrics(VALID_OBS, sim_flow)
# end