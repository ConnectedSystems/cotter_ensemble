using Streamfall: load_calibration
include("_common.jl")
include("_calib_funcs.jl")

using StatsPlots, DataStructures

ensemble_figs = joinpath(FIG_PATH, "ensemble_results", "cmd") * "/"
mkpath(ensemble_figs)

ensemble_res_path = joinpath(DATA_PATH, "ensemble_results") * "/"
mkpath(ensemble_res_path)


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


"""
Converts result Dict to DataFrame and transpose so that each row gives the results for each approach.

Sorts results by KGE'
"""
function collate_results_to_df(res::OrderedDict, approaches::Array)
    tmp_df = DataFrame(res)
    df = DataFrame([[names(tmp_df)]; collect.(eachrow(tmp_df))], [:metric; Symbol.(approaches)])
    sort!(df, Symbol("KGE'"), rev=true)

    return df
end


function calc_metrics!(obs, sim, objs, subset)
    if length(subset) == 0
        return repeat([missing], length(objs))
    end
    
    return [metric(obs[subset], sim[subset]) for metric in objs]
end


baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"
calib_ts = 3651 # 1826

obs_flow = CALIB_OBS[calib_ts:end]

baseline_calib_perf = OrderedDict()
baseline_valid_perf = OrderedDict()

b_objs = vcat([Streamfall.mKGE, Streamfall.NSE], OBJFUNCS)
b_approaches = vcat(["KGE'", "NSE"], APPROACHES)
base_outflows = Dict()
for app in APPROACHES
    # Compare baselines
    base_sn, base_id = setup_network("baselines/cotter_baseline_IHACRES_$(app).yml")
    base_node = base_sn[base_id]
    Streamfall.run_basin!(base_sn, FULL_CLIMATE)

    base_outflows[app] = base_node.outflow

    baseline_calib_perf[app] = [metric(obs_flow, base_node.outflow[calib_ts:CALIB_LN]) for metric in b_objs]
    baseline_valid_perf[app] = [metric(VALID_OBS, base_node.outflow[CALIB_LN+1:end]) for metric in b_objs]
end

# transpose so that each row gives the results for each approach
baseline_calib_df = collate_results_to_df(baseline_calib_perf, b_approaches)
baseline_valid_df = collate_results_to_df(baseline_valid_perf, b_approaches)

state = :storage
state_str = String(state)
target_idx = [1,2,5,6,7,8]
quantiles = [0.0, 0.1, 0.9]
state_calib_overall_perf = OrderedDict()
state_calib_wet_perf = OrderedDict()
state_calib_mid_perf = OrderedDict()
state_calib_dry_perf = OrderedDict()

state_valid_overall_perf = OrderedDict()
state_valid_wet_perf = OrderedDict()
state_valid_mid_perf = OrderedDict()
state_valid_dry_perf = OrderedDict()

baseline_calib_overall_perf = OrderedDict()
baseline_calib_wet_perf = OrderedDict()
baseline_calib_mid_perf = OrderedDict()
baseline_calib_dry_perf = OrderedDict()

baseline_valid_overall_perf = OrderedDict()
baseline_valid_wet_perf = OrderedDict()
baseline_valid_mid_perf = OrderedDict()
baseline_valid_dry_perf = OrderedDict()


state_flows = Dict()
state_active_params = Dict()
for state_approach in APPROACHES
    @info "$(state_approach)"
    # Compare state_based
    # Read in parameters
    calibrated_params = nothing
    try
        fn = "$(calib_param_path)log_$(state_str)_$(state_approach)_10_year.txt"
        calibrated_params = read_params(fn)
    catch err
        if err isa SystemError
            @info "Skipping $(state_approach) as file not available"
            continue
        end
    end

    if length(calibrated_params) != 18
        @info "Skipping $(state_approach) [mismatched params]"
        continue
    end

    sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$state_approach.yml")
    node = sn[n_id]

    _, x0, __ = param_info(node; with_level=false)
    mod_params = update_partial(x0, target_idx, calibrated_params)

    active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; burn_in=calib_ts, log_transform=true, releases=nothing)

    state_flows[state_approach] = node.outflow
    state_active_params[state_approach] = active_param_set

    # Calibration - scores for given approach
    Q_ss = node.outflow[calib_ts:CALIB_LN]
    state_calib_overall_perf["$(state_approach)"] = [metric(obs_flow, Q_ss) for metric in b_objs]

    wet_store = findall(active_param_set[calib_ts:CALIB_LN] .== 1)
    mid_store = findall(active_param_set[calib_ts:CALIB_LN] .== 2)
    dry_store = findall(active_param_set[calib_ts:CALIB_LN] .== 3)

    state_calib_wet_perf["$(state_approach)"] = calc_metrics!(obs_flow, Q_ss, b_objs, wet_store)
    state_calib_mid_perf["$(state_approach)"] = calc_metrics!(obs_flow, Q_ss, b_objs, mid_store)
    state_calib_dry_perf["$(state_approach)"] = calc_metrics!(obs_flow, Q_ss, b_objs, dry_store)

    # Get baseline performance at same temporal locations for each state
    Q_sm = base_outflows[state_approach][calib_ts:CALIB_LN]
    baseline_calib_overall_perf["$(state_approach)"] = [metric(obs_flow, Q_sm) for metric in b_objs]
    baseline_calib_wet_perf["$(state_approach)"] = calc_metrics!(obs_flow, Q_sm, b_objs, wet_store)
    baseline_calib_mid_perf["$(state_approach)"] = calc_metrics!(obs_flow, Q_sm, b_objs, mid_store)
    baseline_calib_dry_perf["$(state_approach)"] = calc_metrics!(obs_flow, Q_sm, b_objs, dry_store)

    # Validation
    Q_ss = node.outflow[CALIB_LN+1:end]
    state_valid_overall_perf["$(state_approach)"] = [metric(VALID_OBS, Q_ss) for metric in b_objs]

    wet_store = findall(active_param_set[CALIB_LN+1:end] .== 1)
    mid_store = findall(active_param_set[CALIB_LN+1:end] .== 2)
    dry_store = findall(active_param_set[CALIB_LN+1:end] .== 3)

    state_valid_wet_perf["$(state_approach)"] = calc_metrics!(VALID_OBS, Q_ss, b_objs, wet_store)
    state_valid_mid_perf["$(state_approach)"] = calc_metrics!(VALID_OBS, Q_ss, b_objs, mid_store)
    state_valid_dry_perf["$(state_approach)"] = calc_metrics!(VALID_OBS, Q_ss, b_objs, dry_store)

    Q_sb = base_outflows[state_approach][CALIB_LN+1:end]
    baseline_valid_overall_perf["$(state_approach)"] = [metric(VALID_OBS, Q_sb) for metric in b_objs]
    baseline_valid_wet_perf["$(state_approach)"] = calc_metrics!(VALID_OBS, Q_sb, b_objs, wet_store)
    baseline_valid_mid_perf["$(state_approach)"] = calc_metrics!(VALID_OBS, Q_sb, b_objs, mid_store)
    baseline_valid_dry_perf["$(state_approach)"] = calc_metrics!(VALID_OBS, Q_sb, b_objs, dry_store)
end


# Calibration sets
state_calib_overall_df = collate_results_to_df(state_calib_overall_perf, b_approaches)

state_calib_wet_df = collate_results_to_df(state_calib_wet_perf, b_approaches)
state_calib_mid_df = collate_results_to_df(state_calib_mid_perf, b_approaches)
state_calib_dry_df = collate_results_to_df(state_calib_dry_perf, b_approaches)

baseline_calib_overall_df = collate_results_to_df(baseline_calib_overall_perf, b_approaches)
baseline_calib_wet_df = collate_results_to_df(baseline_calib_wet_perf, b_approaches)
baseline_calib_mid_df = collate_results_to_df(baseline_calib_mid_perf, b_approaches)
baseline_calib_dry_df = collate_results_to_df(baseline_calib_dry_perf, b_approaches)

# Validation sets
state_valid_overall_df = collate_results_to_df(state_valid_overall_perf, b_approaches)

state_valid_wet_df = collate_results_to_df(state_valid_wet_perf, b_approaches)
state_valid_mid_df = collate_results_to_df(state_valid_mid_perf, b_approaches)
state_valid_dry_df = collate_results_to_df(state_valid_dry_perf, b_approaches)

baseline_valid_overall_df = collate_results_to_df(baseline_valid_overall_perf, b_approaches)
baseline_valid_wet_df = collate_results_to_df(baseline_valid_wet_perf, b_approaches)
baseline_valid_mid_df = collate_results_to_df(baseline_valid_mid_perf, b_approaches)
baseline_valid_dry_df = collate_results_to_df(baseline_valid_dry_perf, b_approaches)


function collate_overview(overall, wet, mid, dry; obj_name, sort_by)
    df_set = [overall[:, ["metric", obj_name]], wet[:, ["metric", obj_name]], mid[:, ["metric", obj_name]], dry[:, ["metric", obj_name]]]

    combined_df = reduce((x,y) -> leftjoin(x,y, on="metric", makeunique=true), df_set)
    rename!(combined_df, ["$(obj_name)_1"=>"$(obj_name) wet", "$(obj_name)_2"=>"$(obj_name) mid", "$(obj_name)_3"=>"$(obj_name) dry"])

    combined_df[:, "Average"] = mean.(skipmissing.(eachrow(combined_df[:, Not(["metric", "KGE'"])])))

    sort!(combined_df, Symbol(sort_by), rev=true)
    return combined_df
end


tgt_metric = "KGE'"
sort_metric = "Average"
df_set = [state_calib_overall_df, state_calib_wet_df, state_calib_mid_df, state_calib_dry_df]
combined_calib_df = collate_overview(df_set...; obj_name=tgt_metric, sort_by=sort_metric)

df_set = [state_valid_overall_df, state_valid_wet_df, state_valid_mid_df, state_valid_dry_df]
combined_valid_df = collate_overview(df_set...; obj_name=tgt_metric, sort_by=sort_metric)


df_set = [baseline_calib_overall_df, baseline_calib_wet_df, baseline_calib_mid_df, baseline_calib_dry_df]
combined_bl_calib_df = collate_overview(df_set...; obj_name=tgt_metric, sort_by=sort_metric)

df_set = [baseline_valid_overall_df, baseline_valid_wet_df, baseline_valid_mid_df, baseline_valid_dry_df]
combined_bl_valid_df = collate_overview(df_set...; obj_name=tgt_metric, sort_by=sort_metric)

# Need to build combinations of baseline/state-based models...


@info "Baseline Calibration" combined_bl_calib_df
@info "State-based Calibration" combined_calib_df


@info "Baseline Validation" combined_bl_valid_df
@info "State-based Validation" combined_valid_df


# Select two to pair up based on performance in the "dry" years

if combined_bl_calib_df[1, "KGE' dry"] > combined_calib_df[1, "KGE' dry"]
    selected_base = combined_bl_calib_df[1, :metric]
    @info "Selected baseline model $(selected_base)"
    base_flow = base_outflows[selected_base][calib_ts:CALIB_LN]

    selected_state_based = combined_calib_df[1, :metric]
else
    selected_base = combined_calib_df[1, :metric]
    @info "Selected state-based model $(selected_base)"
    base_flow = state_flows[selected_base][calib_ts:CALIB_LN]

    selected_state_based = combined_calib_df[5, :metric]
end

@info "Selected model 2: $(selected_state_based)"

bl_flow = base_outflows[selected_base][calib_ts:CALIB_LN]
s_flow = state_flows[selected_state_based][calib_ts:CALIB_LN]
active_params = state_active_params[selected_state_based]

ensemble_flow = copy(base_flow)
low_store = findall(active_params[calib_ts:CALIB_LN] .== 1)
mid_store = findall(active_params[calib_ts:CALIB_LN] .== 2)
high_store = findall(active_params[calib_ts:CALIB_LN] .== 3)

high_base = Streamfall.NmKGE(obs_flow[high_store], base_flow[high_store])
high_state = Streamfall.NmKGE(obs_flow[high_store], s_flow[high_store])
high_factor = high_state >= (high_base * 1.0) ? 0.90 : 0.1

mid_base = Streamfall.NmKGE(obs_flow[mid_store], base_flow[mid_store])
mid_state = Streamfall.NmKGE(obs_flow[mid_store], s_flow[mid_store])
mid_factor = mid_state >= (mid_base * 1.0) ? 0.90 : 0.1

if length(low_store) != 0
    low_base = Streamfall.NmKGE(obs_flow[low_store], base_flow[low_store])
    low_state = Streamfall.NmKGE(obs_flow[low_store], s_flow[low_store])
    low_factor = low_state >= (low_base * 1.0) ? 0.80 : 0.2
else
    low_factor = 0.0
end

ensemble_flow[high_store] = ((s_flow[high_store] .* high_factor) .+ (base_flow[high_store] .* (1.0 - high_factor)))
ensemble_flow[mid_store] = ((s_flow[mid_store] .* mid_factor) .+ (base_flow[mid_store] .* (1.0 - mid_factor)))
ensemble_flow[low_store] = ((s_flow[low_store] .* low_factor) .+ (base_flow[low_store] .* (1.0 - low_factor)))


plot_dates = FULL_DATASET.Date[calib_ts:CALIB_LN]
baseline_plot = begin
    plot(plot_dates, obs_flow, title="Baseline\n($selected_base)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end

ensemble_plot = begin
    plot(plot_dates, obs_flow, title="Ensemble\n($selected_state_based)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end
savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)new_test_$(selected_base)_$(selected_state_based)_calib.png")

@info "Baseline $selected_base" Streamfall.mKGE(obs_flow, bl_flow) Streamfall.mKGE(obs_flow[low_store], bl_flow[low_store]) Streamfall.mKGE(obs_flow[mid_store], bl_flow[mid_store]) Streamfall.mKGE(obs_flow[high_store], bl_flow[high_store])
@info "Model 1 (state-based $selected_base)" Streamfall.mKGE(obs_flow, base_flow) Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]) Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]) Streamfall.mKGE(obs_flow[high_store], base_flow[high_store])
@info "Model 2 (state-based $selected_state_based)" Streamfall.mKGE(obs_flow, s_flow) Streamfall.mKGE(obs_flow[low_store], s_flow[low_store]) Streamfall.mKGE(obs_flow[mid_store], s_flow[mid_store]) Streamfall.mKGE(obs_flow[high_store], s_flow[high_store])
@info "Ensemble" Streamfall.mKGE(obs_flow, ensemble_flow) Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]) Streamfall.mKGE(obs_flow[mid_store], ensemble_flow[mid_store]) Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store])


res = OrderedDict(
    "States" => ["Overall", "Wet", "Usual", "Dry"],
    "Baseline $selected_base" => [Streamfall.mKGE(obs_flow, bl_flow), Streamfall.mKGE(obs_flow[low_store], bl_flow[low_store]), Streamfall.mKGE(obs_flow[mid_store], bl_flow[mid_store]), Streamfall.mKGE(obs_flow[high_store], bl_flow[high_store])],
    "Model 1 (state-based $selected_base)" => [Streamfall.mKGE(obs_flow, base_flow), Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), Streamfall.mKGE(obs_flow[high_store], base_flow[high_store])],
    "Model 2 (state-based $selected_state_based)" => [Streamfall.mKGE(obs_flow, s_flow), Streamfall.mKGE(obs_flow[low_store], s_flow[low_store]), Streamfall.mKGE(obs_flow[mid_store], s_flow[mid_store]), Streamfall.mKGE(obs_flow[high_store], s_flow[high_store])],
    "Ensemble" => [Streamfall.mKGE(obs_flow, ensemble_flow), Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), Streamfall.mKGE(obs_flow[mid_store], ensemble_flow[mid_store]), Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store])],
)
CSV.write(joinpath(ensemble_res_path, "ensemble_calib_results.csv"), DataFrame(res))

# Validation 
v_bl_flow = base_outflows[selected_base][CALIB_LN+1:end]
v_base_flow = state_flows[selected_base][CALIB_LN+1:end]

v_s_flow = state_flows[selected_state_based][CALIB_LN+1:end]
active_params = state_active_params[selected_state_based]

ensemble_flow = copy(v_base_flow)
low_store = findall(active_params[CALIB_LN+1:end] .== 1)
mid_store = findall(active_params[CALIB_LN+1:end] .== 2)
high_store = findall(active_params[CALIB_LN+1:end] .== 3)

# Weighted sum using previously determined weights
ensemble_flow[high_store] = ((v_s_flow[high_store] .* high_factor) .+ (v_base_flow[high_store] .* (1.0 - high_factor)))
ensemble_flow[mid_store] = ((v_s_flow[mid_store] .* mid_factor) .+ (v_base_flow[mid_store] .* (1.0 - mid_factor)))

if length(low_store) != 0
    ensemble_flow[low_store] = ((v_s_flow[low_store] .* low_factor) .+ (v_base_flow[low_store] .* (1.0 - low_factor)))
else
    ensemble_flow[low_store] = v_base_flow[low_store]
end


plot_dates = FULL_DATASET.Date[CALIB_LN+1:end]
baseline_plot = begin
    plot(plot_dates, VALID_OBS, title="Baseline\n($selected_base)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, v_base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end

ensemble_plot = begin
    plot(plot_dates, VALID_OBS, title="Ensemble\n($selected_state_based)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
    plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
    plot!(size=(600,200), dpi=300)
end
savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)new_test_$(selected_base)_$(selected_state_based)_valid.png")

@info "---------------"


if length(low_store) == 0
    bl_res = missing
    base_res = missing
    en2 = missing
    ws_en = missing
else
    bl_res = Streamfall.mKGE(VALID_OBS[low_store], v_bl_flow[low_store])
    base_res = Streamfall.mKGE(VALID_OBS[low_store], v_base_flow[low_store])
    en2 = Streamfall.mKGE(VALID_OBS[low_store], v_s_flow[low_store])
    ws_en = Streamfall.mKGE(VALID_OBS[low_store], ensemble_flow[low_store])
end
    

@info "Baseline $selected_base" Streamfall.mKGE(VALID_OBS, v_bl_flow) bl_res Streamfall.mKGE(VALID_OBS[mid_store], v_bl_flow[mid_store]) Streamfall.mKGE(VALID_OBS[high_store], v_bl_flow[high_store])
@info "Model 1 $selected_base" Streamfall.mKGE(VALID_OBS, v_base_flow) base_res Streamfall.mKGE(VALID_OBS[mid_store], v_base_flow[mid_store]) Streamfall.mKGE(VALID_OBS[high_store], v_base_flow[high_store])
@info "Model 2 $selected_state_based" Streamfall.mKGE(VALID_OBS, v_s_flow) en2 Streamfall.mKGE(VALID_OBS[mid_store], v_s_flow[mid_store]) Streamfall.mKGE(VALID_OBS[high_store], v_s_flow[high_store])
@info "Ensemble" Streamfall.mKGE(VALID_OBS, ensemble_flow) ws_en Streamfall.mKGE(VALID_OBS[mid_store], ensemble_flow[mid_store]) Streamfall.mKGE(VALID_OBS[high_store], ensemble_flow[high_store])


res = OrderedDict(
    "States" => ["Overall", "Wet", "Usual", "Dry"],
    "Baseline $selected_base" => [Streamfall.mKGE(VALID_OBS, v_bl_flow), bl_res, Streamfall.mKGE(VALID_OBS[mid_store], v_bl_flow[mid_store]), Streamfall.mKGE(VALID_OBS[high_store], v_bl_flow[high_store])],
    "Model 1 (state-based $selected_base)" => [Streamfall.mKGE(VALID_OBS, v_base_flow), base_res, Streamfall.mKGE(VALID_OBS[mid_store], v_base_flow[mid_store]), Streamfall.mKGE(VALID_OBS[high_store], v_base_flow[high_store])],
    "Model 2 (state-based $selected_state_based)" => [Streamfall.mKGE(VALID_OBS, v_s_flow), en2, Streamfall.mKGE(VALID_OBS[mid_store], v_s_flow[mid_store]), Streamfall.mKGE(VALID_OBS[high_store], v_s_flow[high_store])],
    "Ensemble" => [Streamfall.mKGE(VALID_OBS, ensemble_flow), ws_en, Streamfall.mKGE(VALID_OBS[mid_store], ensemble_flow[mid_store]), Streamfall.mKGE(VALID_OBS[high_store], ensemble_flow[high_store])],    
)

CSV.write(joinpath(ensemble_res_path, "ensemble_valid_results.csv"), DataFrame(res))


# Write out results

# Collated calib/validations
CSV.write(joinpath(ensemble_res_path, "collated_mKGE_baseline_calib_results.csv"), combined_bl_calib_df)
CSV.write(joinpath(ensemble_res_path, "collated_mKGE_baseline_valid_results.csv"), combined_bl_valid_df)

CSV.write(joinpath(ensemble_res_path, "collated_mKGE_state-based_calib_results.csv"), combined_calib_df)
CSV.write(joinpath(ensemble_res_path, "collated_mKGE_state-based_valid_results.csv"), combined_valid_df)

# All metrics used
CSV.write(joinpath(ensemble_res_path, "overall_baseline_calib_results.csv"), baseline_calib_overall_df)
CSV.write(joinpath(ensemble_res_path, "overall_baseline_valid_results.csv"), baseline_valid_overall_df)

CSV.write(joinpath(ensemble_res_path, "overall_state-based_calib_results.csv"), state_calib_overall_df)
CSV.write(joinpath(ensemble_res_path, "overall_state-based_valid_results.csv"), state_valid_overall_df)

# Individual states
CSV.write(joinpath(ensemble_res_path, "wet_state-based_calib_results.csv"), state_calib_wet_df)
CSV.write(joinpath(ensemble_res_path, "mid_state-based_calib_results.csv"), state_calib_mid_df)
CSV.write(joinpath(ensemble_res_path, "dry_state-based_calib_results.csv"), state_calib_dry_df)

# Individual states for baselines not stored as these depend on the state-based model being compared against



# obs_flow = FULL_DATASET[calib_ts:CALIB_LN, "410730_Q"]
# base_flow = base_node.outflow[calib_ts:CALIB_LN]
# sim_flow = node.outflow[calib_ts:CALIB_LN]


# ensemble_flow = copy(base_flow)
# low_store = findall(active_param_set[calib_ts:CALIB_LN] .== 1)
# mid_store = findall(active_param_set[calib_ts:CALIB_LN] .== 2)
# high_store = findall(active_param_set[calib_ts:CALIB_LN] .== 3)

# @assert length(mid_store) > 0

# high_base = Streamfall.NmKGE(obs_flow[high_store], base_flow[high_store])
# high_state = Streamfall.NmKGE(obs_flow[high_store], sim_flow[high_store])
# high_factor = high_state >= (high_base * 1.0) ? 0.99 : 0.01

# mid_base = Streamfall.NmKGE(obs_flow[mid_store], base_flow[mid_store])
# mid_state = Streamfall.NmKGE(obs_flow[mid_store], sim_flow[mid_store])
# mid_factor = mid_state >= (mid_base * 1.0) ? 0.8 : 0.2

# low_base = Streamfall.NmKGE(obs_flow[low_store], base_flow[low_store])
# low_state = Streamfall.NmKGE(obs_flow[low_store], sim_flow[low_store])
# low_factor = low_state >= (low_base * 1.0) ? 1.0 : 0.0

# # ensemble_flow[high_store] = ((ensemble_flow[high_store] .* high_factor) .+ (sim_flow[high_store] .* (1.0 - high_factor)))
# # ensemble_flow[mid_store] = ((ensemble_flow[mid_store] .* mid_factor) .+ (sim_flow[mid_store] .* (1.0 - mid_factor)))
# # ensemble_flow[low_store] = ((ensemble_flow[low_store] .* low_factor) .+ (sim_flow[low_store] .* (1.0 - low_factor)))
# ensemble_flow[high_store] = ((sim_flow[high_store] .* high_factor) .+ (ensemble_flow[high_store] .* (1.0 - high_factor)))
# ensemble_flow[mid_store] = ((sim_flow[mid_store] .* mid_factor) .+ (ensemble_flow[mid_store] .* (1.0 - mid_factor)))
# ensemble_flow[low_store] = ((sim_flow[low_store] .* low_factor) .+ (ensemble_flow[low_store] .* (1.0 - low_factor)))

# @info "Calibration metrics for $(baseline_approach) [High/Mid/Low]"
# report_metrics(obs_flow[high_store], base_flow[high_store])
# report_metrics(obs_flow[mid_store], base_flow[mid_store])
# report_metrics(obs_flow[low_store], base_flow[low_store])


# @info "Simulation metrics for $(state_approach) [Overall/High/Mid/Low]"
# report_metrics(obs_flow, sim_flow)
# report_metrics(obs_flow[high_store], sim_flow[high_store])
# report_metrics(obs_flow[mid_store], sim_flow[mid_store])
# report_metrics(obs_flow[low_store], sim_flow[low_store])

# @info "Metrics for Ensemble [High/Mid/Low]"
# report_metrics(obs_flow, ensemble_flow)
# report_metrics(obs_flow[high_store], ensemble_flow[high_store])
# report_metrics(obs_flow[mid_store], ensemble_flow[mid_store])
# report_metrics(obs_flow[low_store], ensemble_flow[low_store])


# # sanity_check
# @info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], sim_flow[mid_store]) Streamfall.RMSE(obs_flow[low_store], sim_flow[low_store])
# @info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)


# plot_dates = FULL_DATASET.Date[calib_ts:CALIB_LN]
# baseline_plot = begin
#     plot(plot_dates, obs_flow, title="Baseline\n($baseline_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
#     plot!(plot_dates, base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
#     plot!(size=(600,200), dpi=300)
# end

# ensemble_plot = begin
#     plot(plot_dates, obs_flow, title="Ensemble\n($state_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
#     plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
#     plot!(size=(600,200), dpi=300)
# end
# savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)flow_cmd_$(baseline_approach)_$(state_approach)_calib.png")

# base_store = base_node.gw_store[calib_ts:CALIB_LN+1]
# pos_markers = quantile(base_store, quantiles[2:end])

# base_low_store = findall(base_store .< pos_markers[1])
# base_mid_store = findall(pos_markers[1] .>= base_store .< pos_markers[2])
# base_high_store = findall(base_store .>= pos_markers[2])

# score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)
# high = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
# mid = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
# low = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
# base_qq = qqplot(obs_flow, base_flow,
#                  title="Baseline $(baseline_approach)\n(KGE': $score; Wet: $low Mid: $mid; Dry: $high)",
#                  titlefontsize=8,
#                  markerstrokewidth=0,
#                  markersize=2,
#                  alpha=0.5,
#                  )

# score = round(Streamfall.mKGE(obs_flow, sim_flow), digits=3)
# high = round(Streamfall.mKGE(obs_flow[high_store], sim_flow[high_store]), digits=3)
# mid = round(Streamfall.mKGE(obs_flow[mid_store], sim_flow[mid_store]), digits=3)
# low = round(Streamfall.mKGE(obs_flow[low_store], sim_flow[low_store]), digits=3)
# sim_qq = qqplot(obs_flow, sim_flow,
#                 title="State-based $state_approach\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2,
#                 alpha=0.5,
#                 )

# sim_qq_col = begin
#     # plot(
#     #     sim_qq,
#     #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
#     #     dpi=300
#     # )

#     low_plot = begin
#         scatter(obs_flow[low_store], ensemble_flow[low_store],
#                 title="State-based Wet\n(KGE': $(low))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     mid_plot = begin
#         scatter(obs_flow[mid_store], ensemble_flow[mid_store],
#                 title="State-based Mid\n(KGE': $(mid))",
#                 titlefontsize=8,
#                 markersize=2.0,
#                 markerstrokewidth=0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     high_plot = begin
#         scatter(obs_flow[high_store], ensemble_flow[high_store],
#                 title="State-based Dry\n(KGE': $(high))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     plot(
#         sim_qq,
#         low_plot,
#         mid_plot,
#         high_plot,
#         dpi=300
#     )
# end

# savefig(sim_qq_col, "$(ensemble_figs)sim_cmd_$(baseline_approach)_$(state_approach)_qq_calib.png")


# score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=3)
# high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=3)
# mid = round(Streamfall.mKGE(obs_flow[mid_store], ensemble_flow[mid_store]), digits=3)
# low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=3)

# ensemble_qq = qqplot(obs_flow, ensemble_flow,
#                      title="$(baseline_approach) - $(state_approach)\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
#                      titlefontsize=8,
#                      markerstrokewidth=0,
#                      alpha=0.5,
#                      markersize=2)

# ensemble_col = begin
#     # plot(
#     #     ensemble_qq,
#     #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
#     #     dpi=300
#     # )
#     high_plot = begin
#         scatter(obs_flow[high_store], ensemble_flow[high_store],
#                 title="Dry Ensemble\n(KGE': $(high))",
#                 titlefontsize=8,
#                 markersize=2.0,
#                 markerstrokewidth=0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     mid_plot = begin
#         scatter(obs_flow[mid_store], ensemble_flow[mid_store],
#                 title="Mid Ensemble\n(KGE': $(mid))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     low_plot = begin
#         scatter(obs_flow[low_store], ensemble_flow[low_store],
#                 title="Wet Ensemble\n(KGE': $(low))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     plot(
#         ensemble_qq,
#         low_plot,
#         mid_plot,
#         high_plot,
#         dpi=300
#     )
# end
# savefig(ensemble_col, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_qq_calib.png")


# snapshot = begin
#     plot(
#         base_qq,
#         sim_qq,
#         ensemble_qq,
#         ensemble_plot,
#         dpi=300
#     )
#     plot!(size=(600,400))
# end

# savefig(snapshot, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_result_calib.png")



# high_store_plot = begin
#     score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
#     plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false)
#     plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
# end

# mid_store_plot = begin
#     score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
#     plot(plot_dates[mid_store], obs_flow[mid_store], title="Mid GW Store\n(KGE': $score)", legend=false)
#     plot!(plot_dates[mid_store], base_flow[mid_store], alpha=0.6)
# end

# low_store_plot = begin
#     score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
#     plot(plot_dates[low_store], obs_flow[low_store], title="Low GW Store\n(KGE': $score)", legend=false)
#     plot!(plot_dates[low_store], base_flow[low_store], alpha=0.6)
# end

# baseline_qq = begin
#     overall_score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)
#     high_score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
#     mid_score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
#     low_score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)

#     overall = begin
#         scatter(obs_flow, base_flow,
#                 title="Baseline $(baseline_approach)\n(KGE': $overall_score)",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow, obs_flow, alpha=0.6, color="green")
#     end

#     high_plot = begin
#         scatter(obs_flow[high_store], base_flow[high_store],
#                 title="Dry $(baseline_approach)\n(KGE': $(high_score))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[high_store], obs_flow[high_store], alpha=0.6, color="green")
#     end

#     mid_plot = begin
#         scatter(obs_flow[mid_store], base_flow[mid_store],
#                 title="Mid $(baseline_approach)\n(KGE': $(mid_score))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[mid_store], obs_flow[mid_store], alpha=0.6, color="green")
#     end

#     low_plot = begin
#         scatter(obs_flow[low_store], base_flow[low_store],
#                 title="Wet $(baseline_approach)\n(KGE': $(low_score))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[low_store], obs_flow[low_store], alpha=0.6, color="green")
#     end
#     # plot(
#     #     qqplot(obs_flow, base_flow, title="Baseline $(baseline_approach)\n(KGE': $overall_score)", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[low_store], base_flow[low_store], title="Wet $(baseline_approach)\n(KGE': $(low_score))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[mid_store], base_flow[mid_store], title="Mid $(baseline_approach)\n(KGE': $(mid_score))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[high_store], base_flow[high_store], title="Dry $(baseline_approach)\n(KGE': $(high_score))", titlefontsize=8, markersize=2.0),
#     # )
#     plot(
#         overall,
#         low_plot,
#         mid_plot,
#         high_plot,
#         dpi=300
#     )
#     plot!(size=(600,400))
# end
# savefig(baseline_qq, joinpath(ensemble_figs, "baseline_cmd_$(baseline_approach)_$(state_approach)_qq_calib.png"))



# ### Validation ###

# obs_flow = FULL_DATASET[CALIB_LN+1:end, "410730_Q"]
# base_flow = base_node.outflow[CALIB_LN+1:end]
# sim_flow = node.outflow[CALIB_LN+1:end]


# ensemble_flow = copy(base_flow)
# low_store = findall(active_param_set[CALIB_LN+1:end] .== 1)
# mid_store = findall(active_param_set[CALIB_LN+1:end] .== 2)
# high_store = findall(active_param_set[CALIB_LN+1:end] .== 3)

# @assert length(mid_store) > 0

# ensemble_flow[high_store] = ((sim_flow[high_store] .* high_factor) .+ (ensemble_flow[high_store] .* (1.0 - high_factor)))
# ensemble_flow[mid_store] = ((sim_flow[mid_store] .* mid_factor) .+ (ensemble_flow[mid_store] .* (1.0 - mid_factor)))

# if length(low_store) > 0
#     ensemble_flow[low_store] = ((sim_flow[low_store] .* low_factor) .+ (ensemble_flow[low_store] .* (1.0 - low_factor)))
# end

# # sanity_check
# # @info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], sim_flow[mid_store])
# # @info "Bin Metric results" bin_metric(active_param_set[calib_ts:CALIB_LN], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)
# # sanity_check
# @info "Sanity Check [high/low]" Streamfall.RMSE(obs_flow[high_store], sim_flow[high_store]) Streamfall.RMSE(obs_flow[mid_store], sim_flow[mid_store]) Streamfall.RMSE(obs_flow[low_store], sim_flow[low_store])
# @info "Bin Metric results" bin_metric(active_param_set[CALIB_LN+1:end], param_idxs, obs_flow, sim_flow, split_NNSE, nothing)

# plot_dates = FULL_DATASET.Date[CALIB_LN+1:end]
# baseline_plot = begin
#     plot(plot_dates, obs_flow, title="Baseline\n($baseline_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
#     plot!(plot_dates, base_flow, alpha=0.8, linewidth=1.0, linestyle=:dashdot)
#     plot!(size=(600,200), dpi=300)
# end

# ensemble_plot = begin
#     plot(plot_dates, obs_flow, title="Ensemble\n($state_approach)", titlefontsize=8, legendfontsize=8, linewidth=1.5, legend=false)
#     plot!(plot_dates, ensemble_flow, alpha=0.8, linestyle=:dashdot)
#     plot!(size=(600,200), dpi=300)
# end
# savefig(plot(baseline_plot, ensemble_plot), "$(ensemble_figs)flow_cmd_$(baseline_approach)_$(state_approach)_valid.png")

# base_store = base_node.gw_store[CALIB_LN+1:end]
# pos_markers = quantile(base_store, quantiles[2:end])

# base_low_store = findall(base_store .< pos_markers[1])
# base_mid_store = findall(pos_markers[1] .>= base_store .< pos_markers[2])
# base_high_store = findall(base_store .>= pos_markers[2])

# score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)

# high = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)

# if length(mid_store) > 0
#     mid = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
# else
#     mid = NaN
# end


# if length(low_store) > 0
#     low = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
# else
#     low = NaN
# end

# base_qq = qqplot(obs_flow, base_flow,
#                  title="Baseline $(baseline_approach)\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
#                  titlefontsize=8,
#                  markerstrokewidth=0,
#                  markersize=2,
#                  alpha=0.5,
#                  )

# score = round(Streamfall.mKGE(obs_flow, sim_flow), digits=3)

# high = round(Streamfall.mKGE(obs_flow[high_store], sim_flow[high_store]), digits=3)

# if length(mid_store) > 0
#     mid = round(Streamfall.mKGE(obs_flow[mid_store], sim_flow[mid_store]), digits=3)
# else
#     mid = NaN
# end


# if length(low_store) > 0
#     low = round(Streamfall.mKGE(obs_flow[low_store], sim_flow[low_store]), digits=3)
# else
#     low = NaN
# end

# sim_qq = qqplot(obs_flow, sim_flow,
#                 title="State-based $state_approach\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2,
#                 alpha=0.5,
#                 )

# sim_qq_col = begin
#     # plot(
#     #     sim_qq,
#     #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
#     #     dpi=300
#     # )
#     high_plot = begin
#         scatter(obs_flow[high_store], ensemble_flow[high_store],
#                 title="State-based Dry\n(KGE': $(high))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     mid_plot = begin
#         scatter(obs_flow[mid_store], ensemble_flow[mid_store],
#                 title="State-based Mid\n(KGE': $(mid))",
#                 titlefontsize=8,
#                 markersize=2.0,
#                 markerstrokewidth=0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     low_plot = begin
#         scatter(obs_flow[low_store], ensemble_flow[low_store],
#                 title="State-based Wet\n(KGE': $(low))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     plot(
#         sim_qq,
#         low_plot,
#         mid_plot,
#         high_plot,
#         dpi=300
#     )
# end

# savefig(sim_qq_col, "$(ensemble_figs)sim_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png")


# score = round(Streamfall.mKGE(obs_flow, ensemble_flow), digits=3)

# high = round(Streamfall.mKGE(obs_flow[high_store], ensemble_flow[high_store]), digits=3)

# if length(mid_store) > 0
#     mid = round(Streamfall.mKGE(obs_flow[mid_store], ensemble_flow[mid_store]), digits=3)
# else
#     mid = NaN
# end

# if length(low_store) > 0
#     low = round(Streamfall.mKGE(obs_flow[low_store], ensemble_flow[low_store]), digits=3)
# else
#     low = NaN
# end

# ensemble_qq = qqplot(obs_flow, ensemble_flow,
#                      title="$(baseline_approach) - $(state_approach)\n(KGE': $score; Wet: $low; Mid: $mid; Dry: $high)",
#                      titlefontsize=8,
#                      markerstrokewidth=0,
#                      alpha=0.5,
#                      markersize=2)

# ensemble_col = begin
#     # plot(
#     #     ensemble_qq,
#     #     qqplot(obs_flow[low_store], ensemble_flow[low_store], title="Wet Ensemble\n(KGE': $(low))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[mid_store], ensemble_flow[mid_store], title="Mid Ensemble\n(KGE': $(mid))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[high_store], ensemble_flow[high_store], title="Dry Ensemble\n(KGE': $(high))", titlefontsize=8, markersize=2.0),
#     #     dpi=300
#     # )
#     high_plot = begin
#         scatter(obs_flow[high_store], ensemble_flow[high_store],
#                 title="Dry Ensemble\n(KGE': $(high))",
#                 titlefontsize=8,
#                 markersize=2.0,
#                 markerstrokewidth=0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[high_store], obs_flow[high_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     mid_plot = begin
#         scatter(obs_flow[mid_store], ensemble_flow[mid_store],
#                 title="Mid Ensemble\n(KGE': $(mid))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[mid_store], obs_flow[mid_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     low_plot = begin
#         scatter(obs_flow[low_store], ensemble_flow[low_store],
#                 title="Wet Ensemble\n(KGE': $(low))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[low_store], obs_flow[low_store], label=false, color="green", alpha=0.8)  # 1:1 line
#     end

#     plot(
#         ensemble_qq,
#         low_plot,
#         mid_plot,
#         high_plot,
#         dpi=300
#     )
# end
# savefig(ensemble_col, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png")


# snapshot = begin
#     plot(
#         base_qq,
#         sim_qq,
#         ensemble_qq,
#         ensemble_plot,
#         dpi=300
#     )
#     plot!(size=(600,400))
# end

# savefig(snapshot, "$(ensemble_figs)ensemble_cmd_$(baseline_approach)_$(state_approach)_result_valid.png")



# high_store_plot = begin
#     score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
#     plot(plot_dates[high_store], obs_flow[high_store], title="High GW Store\n(KGE': $score)", legend=false)
#     plot!(plot_dates[high_store], base_flow[high_store], alpha=0.6)
# end

# mid_store_plot = begin
#     score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)
#     plot(plot_dates[mid_store], obs_flow[mid_store], title="Mid GW Store\n(KGE': $score)", legend=false)
#     plot!(plot_dates[mid_store], base_flow[mid_store], alpha=0.6)
# end

# low_store_plot = begin
#     if length(low_store) > 0
#         score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
#     else
#         score = NaN
#     end
#     plot(plot_dates[low_store], obs_flow[low_store], title="Low GW Store\n(KGE': $score)", legend=false)
#     plot!(plot_dates[low_store], base_flow[low_store], alpha=0.6)
# end

# baseline_qq = begin
#     overall_score = round(Streamfall.mKGE(obs_flow, base_flow), digits=3)
#     high_score = round(Streamfall.mKGE(obs_flow[high_store], base_flow[high_store]), digits=3)
#     mid_score = round(Streamfall.mKGE(obs_flow[mid_store], base_flow[mid_store]), digits=3)

#     if length(low_store) > 0
#         low_score = round(Streamfall.mKGE(obs_flow[low_store], base_flow[low_store]), digits=3)
#     else
#         low_score = NaN
#     end

#     overall = begin
#         scatter(obs_flow, base_flow,
#                 title="Baseline $(baseline_approach)\n(KGE': $overall_score)",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow, obs_flow, alpha=0.6, color="green")
#     end

#     high_plot = begin
#         scatter(obs_flow[high_store], base_flow[high_store],
#                 title="Dry $(baseline_approach)\n(KGE': $(high_score))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[high_store], obs_flow[high_store], alpha=0.6, color="green")
#     end

#     mid_plot = begin
#         scatter(obs_flow[mid_store], base_flow[mid_store],
#                 title="Mid $(baseline_approach)\n(KGE': $(mid_score))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[mid_store], obs_flow[mid_store], alpha=0.6, color="green")
#     end

#     low_plot = begin
#         scatter(obs_flow[low_store], base_flow[low_store],
#                 title="Wet $(baseline_approach)\n(KGE': $(low_score))",
#                 titlefontsize=8,
#                 markerstrokewidth=0,
#                 markersize=2.0,
#                 alpha=0.5,
#                 legend=false)
#         plot!(obs_flow[low_store], obs_flow[low_store], alpha=0.6, color="green")
#     end
#     # plot(
#     #     qqplot(obs_flow, base_flow, title="Baseline $(baseline_approach)\n(KGE': $overall_score)", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[low_store], base_flow[low_store], title="Wet $(baseline_approach)\n(KGE': $(low_score))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[mid_store], base_flow[mid_store], title="Mid $(baseline_approach)\n(KGE': $(mid_score))", titlefontsize=8, markersize=2.0),
#     #     qqplot(obs_flow[high_store], base_flow[high_store], title="Dry $(baseline_approach)\n(KGE': $(high_score))", titlefontsize=8, markersize=2.0),
#     # )
#     plot(
#         overall,
#         low_plot,
#         mid_plot,
#         high_plot,
#         dpi=300
#     )
#     plot!(size=(600,400))
# end
# savefig(baseline_qq, joinpath(ensemble_figs, "baseline_cmd_$(baseline_approach)_$(state_approach)_qq_valid.png"))
