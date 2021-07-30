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
state_approach = "split_NnpKGE"
# state_approach = "split_NNSE"
target_idx = [5,6,7,8]

baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"

sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$state_approach.yml")
node = sn[n_id]

# Read in parameters
fn = "$(calib_param_path)restricted_log_$(state_str)_$(state_approach).txt"
calibrated_params = read_params(fn)

reset!(sn)
calib_ts = (1825*2)+1
_, x0, __ = param_info(node; with_level=false)
# calibrated_params = [4.1369754778087495,0.07784403998246077,3.879663457396416,0.2748380354619563, 4.1369754778087495,0.07784403998246077,3.879663457396416,0.2748380354619563,5.743645539293278,0.07237237294614068,9.171393274757758,0.398705852449268]
mod_params = update_partial(x0, target_idx, calibrated_params)
quantiles = [0.0, 0.1, 0.99]
active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; burn_in=calib_ts, log_transform=false, releases=nothing)




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

# prep ensemble

low_store = findall(active_param_set[calib_ts:CALIB_LN] .== 1)
mid_store = findall(active_param_set[calib_ts:CALIB_LN] .== 2)
high_store = findall(active_param_set[calib_ts:CALIB_LN] .== 3)
ensemble_flow[high_store] = ((ensemble_flow[high_store] .* 0.95) .+ (sim_flow[high_store] .* 0.05))
ensemble_flow[mid_store] = ((ensemble_flow[mid_store] .* 0.05) .+ (sim_flow[mid_store] .* 0.95))

@info "Calibration metrics for $(baseline_approach) [High/Low]"
report_metrics(obs_flow[high_store], base_flow[high_store])
report_metrics(obs_flow[Not(high_store)], base_flow[Not(high_store)])


@info "Baseline metrics for $(state_approach) [Overall/High/Low]"
report_metrics(obs_flow, sim_flow)
report_metrics(obs_flow[high_store], sim_flow[high_store])
report_metrics(obs_flow[Not(high_store)], sim_flow[Not(high_store)])

@info "Metrics for Ensemble [High/Low]"
report_metrics(obs_flow, ensemble_flow)
report_metrics(obs_flow[high_store], ensemble_flow[high_store])
report_metrics(obs_flow[Not(high_store)], ensemble_flow[Not(high_store)])

# sanity_check
@info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[Not(high_store)], sim_flow[Not(high_store)])
@info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)

# @info "Pbias?" Streamfall.PBIAS(obs_flow[high_store], sim_flow[high_store]) Streamfall.PBIAS(obs_flow[Not(high_store)], sim_flow[Not(high_store)])


plot_dates = FULL_DATASET.Date[calib_ts:CALIB_LN]
flow_plot = begin
    plot(plot_dates, obs_flow, title="Streamflow", titlefontsize=8, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.6)
end

base_high_store_marker = quantile(base_node.gw_store[calib_ts:CALIB_LN], 0.99)
base_high_store = findall(base_node.gw_store[calib_ts:CALIB_LN] .< base_high_store_marker)
score = round(Streamfall.mKGE(obs_flow, base_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[base_high_store], base_flow[base_high_store]), digits=2)
low = round(Streamfall.mKGE(obs_flow[Not(base_high_store)], base_flow[Not(base_high_store)]), digits=2)
base_qq = qqplot(obs_flow, base_flow, 
                 title="Baseline $(baseline_approach)\n(KGE': $score; High Store: $high; Low: Store: $low)",
                 titlefontsize=8,
                 markerstrokewidth=0,
                 markersize=2,
                 )

score = round(Streamfall.mKGE(obs_flow, sim_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_store], sim_flow[high_store]), digits=2)
low = round(Streamfall.mKGE(obs_flow[Not(high_store)], sim_flow[Not(high_store)]), digits=2)
sim_qq = qqplot(obs_flow, sim_flow,
                title="State-based $state_approach\n(KGE': $score; High Store: $high; Low: Store: $low)",
                titlefontsize=8,
                markerstrokewidth=0,
                markersize=2,
                )

score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=2)
high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=2)
low = round(Streamfall.mKGE(obs_flow[Not(high_store)], ensemble_flow[Not(high_store)]), digits=2)
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
), "$(FIG_PATH)ensemble_$(baseline_approach)_$(state_approach)_result.png")


high_store_plot = begin
    score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=2)
    plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false) 
    plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
end

low_flow_plot = begin
    score = round(Streamfall.mKGE(obs_flow[Not(high_store)], base_flow[Not(high_store)]), digits=2)
    plot(plot_dates[Not(high_store)], obs_flow[Not(high_store)], title="Low GW Store\n(KGE': $score)", legend=false) 
    plot!(plot_dates[Not(high_store)], base_flow[Not(high_store)], alpha=0.6)
end

display(plot(
    qqplot(obs_flow[Not(high_store)], base_flow[Not(high_store)], title="Baseline $(baseline_approach)", markersize=2.0),
    high_store_plot,
    low_flow_plot
))
