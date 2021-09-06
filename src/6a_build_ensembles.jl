# include("6a_comparison.jl")

include("_common.jl")
include("_calib_funcs.jl")

using Streamfall: load_calibration

using StatsPlots, DataStructures
using Infiltrator
using HypothesisTests


ensemble_figs = joinpath(FIG_PATH, "ensemble_results", "cmd") * "/"
mkpath(ensemble_figs)

ensemble_res_path = joinpath(DATA_PATH, "ensemble_results") * "/"
mkpath(ensemble_res_path)


baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"



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


function separate_state_periods(activity)
    low = findall(activity .== 1)
    mid = findall(activity .== 2)
    high = findall(activity .== 3)

    return low, mid, high
end


function separate_state_periods(activity, t_min, t_max=nothing)
    if !isnothing(t_max)
        active = activity[t_min:t_max]
    else
        active = activity[t_min:end]
    end

    return separate_state_periods(active)
end


function build_ensemble_flow(Qb, Qm, state_factors, activity)
    Qe = deepcopy(Qb)

    # Weight towards whichever has the better score in calibration period
    low_factor, mid_factor, high_factor = state_factors

    # Create ensemble outflow for entire simulation period
    low_store, mid_store, high_store = separate_state_periods(activity)

    Qe[high_store] = (Qm[high_store] .* high_factor) .+ (Qb[high_store] .* (1.0 - high_factor))
    Qe[mid_store] = (Qm[mid_store] .* mid_factor) .+ (Qb[mid_store] .* (1.0 - mid_factor))
    Qe[low_store] = (Qm[low_store] .* low_factor) .+ (Qb[low_store] .* (1.0 - low_factor))

    return Qe
end


function combine_ensemble(x, Qb, Qm, active_params)
    if all(isapprox.(x, 1.0)) | all(isapprox.(x, 0.0; atol=1e-6))
        return -9999.9, -9999.9, -9999.9
    end

    # Create ensemble outflow for entire simulation period
    Qe = build_ensemble_flow(Qb, Qm, x, active_params)

    # Subset to calibration period only
    calib_low_store, calib_mid_store, calib_high_store = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
    Qe_calib = Qe[calib_ts+1:CALIB_LN]
    # overall_mKGE = Streamfall.mKGE(Qo_calib, Qe_calib)
    # overall_A2k = KSampleADTest(Qo_calib, Qe_calib).A²k
    low_score = Streamfall.mKGE(Qo_calib[calib_low_store], Qe_calib[calib_low_store])
    mid_score = Streamfall.mKGE(Qo_calib[calib_mid_store], Qe_calib[calib_mid_store])
    high_score = Streamfall.mKGE(Qo_calib[calib_high_store], Qe_calib[calib_high_store])
    
    # low_A2k = KSampleADTest(Qo_calib[calib_low_store], Qe_calib[calib_low_store]).A²k
    # mid_A2k = KSampleADTest(Qo_calib[calib_mid_store], Qe_calib[calib_mid_store]).A²k
    # high_A2k = KSampleADTest(Qo_calib[calib_high_store], Qe_calib[calib_high_store]).A²k
    # mean_A2k = mean([low_A2k, mid_A2k, high_A2k])

    # ensemble_valid = Qe[CALIB_LN+1:end]
    # valid_low_store, valid_mid_store, valid_high_store = separate_state_periods(active_params, CALIB_LN)
    # high_score = Streamfall.mKGE(Qo_valid[valid_high_store], ensemble_valid[valid_high_store])
    # if length(valid_high_store) == 0
    #     return (-9999.9, -9999.9, -9999.9, -9999.9)
    # end

    return (low_score, mid_score, high_score)
end


"""Create base-base ensemble instance"""
function base_base_score(x; model1="split_mean_NmKGE", model2="NmKGE")
    local active_params = state_active_params[model1]
    local Qb = base_Q[model1]
    local Qm = base_Q[model2]

    return combine_ensemble(repeat(x, 3), Qb, Qm, active_params)
end


"""Create base-state ensemble instance"""
function base_state_score(x; model1="split_mean_NmKGE", model2="NmKGE")
    local active_params = state_active_params[model2]
    local Qb = base_Q[model1]
    local Qm = state_Qe[model2]

    return combine_ensemble(x, Qb, Qm, active_params)
end


"""Create state-state ensemble instance"""
function state_state_score(x; model1="split_mean_NmKGE", model2="NmKGE")
    local active_params = state_active_params[model1]
    local Qb = state_Qe[model1]
    local Qm = state_Qe[model2]

    return combine_ensemble(x, Qb, Qm, active_params)
end


function insert_NmKGE_scores(df, approach, Qo, Qm)
    local Qm_calib = Qm[calib_ts+1:CALIB_LN]
    local Qm_valid = Qm[CALIB_LN+1:end]

    if length(valid_low) == 0
        low_score = missing
    else
        low_score = Streamfall.NmKGE(Qo_valid[valid_low], Qm_valid[valid_low])
    end

    if length(valid_high) == 0
        high_score = missing
    else
        high_score = Streamfall.NmKGE(Qo_valid[valid_high], Qm_valid[valid_high])
    end

    push!(df, [approach,
          Streamfall.NmKGE(Qo, Qm),

          Streamfall.NmKGE(Qo_calib, Qm_calib),
          Streamfall.NmKGE(Qo_calib[calib_low], Qm_calib[calib_low]),
          Streamfall.NmKGE(Qo_calib[calib_mid], Qm_calib[calib_mid]),
          Streamfall.NmKGE(Qo_calib[calib_high], Qm_calib[calib_high]),

          Streamfall.NmKGE(Qo_valid, Qm_valid),
          low_score,
          Streamfall.NmKGE(Qo_valid[valid_mid], Qm_valid[valid_mid]),
          high_score,
    ])
end


function insert_A2k_scores(df, approach, Qo, Qm)
    local Qm_calib = Qm[calib_ts+1:CALIB_LN]
    local Qm_valid = Qm[CALIB_LN+1:end]

    AD_low = KSampleADTest(Qo_calib[calib_low], Qm_calib[calib_low]).A²k
    AD_mid = KSampleADTest(Qo_calib[calib_mid], Qm_calib[calib_mid]).A²k
    AD_high = KSampleADTest(Qo_calib[calib_high], Qm_calib[calib_high]).A²k

    if length(valid_low) == 0
        low_score = missing
    else
        low_score = KSampleADTest(Qo_valid[valid_low], Qm_valid[valid_low]).A²k
    end

    if length(valid_high) == 0
        high_score = missing
    else
        high_score = KSampleADTest(Qo_valid[valid_high], Qm_valid[valid_high]).A²k
    end

    push!(df, [approach,
          KSampleADTest(Qo, Qm).A²k,  # Overall A2k score
          mean([AD_low, AD_mid, AD_high]),  # Mean of states in calibration period

          KSampleADTest(Qo_calib, Qm_calib).A²k,
          AD_low,
          AD_mid,
          AD_high,

          KSampleADTest(Qo_valid, Qm_valid).A²k,
          low_score,
          KSampleADTest(Qo_valid[valid_mid], Qm_valid[valid_mid]).A²k,
          high_score,
    ])
end


calib_ts = 3651
state = :storage
state_str = String(state)

target_idx = [1,2,4,5,6,7,8]
quantiles = [0.0, 0.1, 0.9]
# quantiles = [0.0, 0.16, 0.84]

Qo = FULL_DATASET[:, "410730_Q"]
Qo_calib = CALIB_OBS[calib_ts+1:CALIB_LN]
Qo_valid = VALID_OBS

# Recreate results if they do not exist
if (!@isdefined baseline_NmKGE_results) | !isfile(joinpath(ensemble_res_path, "overall_baseline_NmKGE_results.csv"))
    # Record of ensemble flow
    base_Q = Dict()
    state_Qe = Dict()

    # catchment state indicated by each state instance
    state_active_params = Dict()

    # Records of CMD
    base_cmd = Dict()
    state_cmd = Dict()

    for app in APPROACHES
        # baseline
        base_sn, base_id = setup_network("baselines/cotter_baseline_IHACRES_$(app).yml")
        base_node = base_sn[base_id]
        Streamfall.run_basin!(base_sn, FULL_CLIMATE)

        base_Q[app] = base_node.outflow
        base_cmd[app] = base_node.storage[2:end]

        # State-based
        # Read in parameters
        calibrated_state_params = nothing
        try
            fn = "$(calib_param_path)log_$(state_str)_$(app)_10_year.txt"
            calibrated_state_params = read_params(fn)
        catch err
            if err isa SystemError
                @info "Skipping $(app) as file not available"
                continue
            end
        end

        sn, n_id = setup_network("baselines/cotter_baseline_IHACRES_$app.yml")
        node = sn[n_id]

        _, x0, __ = param_info(node; with_level=false)
        mod_params = update_partial(x0, target_idx, calibrated_state_params)

        active_params, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, quantiles; burn_in=calib_ts, log_transform=true, releases=nothing)

        state_Qe[app] = node.outflow
        state_cmd[app] = node.storage
        state_active_params[app] = active_params
    end

    base_base_Q = Dict()
    base_state_Qe = Dict()
    state_state_Qe = Dict()

    base_base_mix = Dict()
    base_state_mix = Dict()
    state_state_mix = Dict()

    # Optimize ensembles and build comparison DFs
    for app in APPROACHES
        for state_approach in APPROACHES
            # Compare baseline-state_based combinations
            local active_params = state_active_params[app]
            local Qb = base_Q[app]
            local Qm = state_Qe[state_approach]

            ensemble_app = (app, state_approach)
            @info "Mixing " ensemble_app

            opt_func = (x) -> base_state_score(x; model1=app, model2=state_approach)
            opt = bbsetup(opt_func; SearchRange=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
                        MaxTime=900, TraceInterval=30,
                        MaxStepsWithoutProgress=50_000,
                        Method=:borg_moea,
                        ϵ=0.05,
                        FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=false))

            res = bboptimize(opt)
            bc_res = best_candidate(res)
            base_state_Qe[ensemble_app] = build_ensemble_flow(Qb, Qm, bc_res, active_params)
            base_state_mix[ensemble_app] = bc_res

            # Compare pairs built with state-based models only
            # Only compare unique combinations
            if app == state_approach
                continue
            end

            # Create base-base ensemble
            local Qb = base_Q[app]
            local Qm = base_Q[state_approach]
            opt_func = (x) -> base_base_score(x; model1=app, model2=state_approach)
            opt = bbsetup(opt_func; SearchRange=[(0.0, 1.0)],
                        MaxTime=900, TraceInterval=30,
                        MaxStepsWithoutProgress=50_000,
                        Method=:borg_moea,
                        ϵ=0.05,
                        FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=false))
            res = bboptimize(opt)
            mixer = best_candidate(res)[1]
            bc_res = [mixer, mixer, mixer]
            base_base_Q[ensemble_app] = build_ensemble_flow(Qb, Qm, bc_res, active_params)
            base_base_mix[ensemble_app] = bc_res

            # Build state-based pairs
            local Qb = state_Qe[app]
            local Qm = state_Qe[state_approach]
            opt_func = (x) -> state_state_score(x; model1=app, model2=state_approach)
            opt = bbsetup(opt_func; SearchRange=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
                        MaxTime=900, TraceInterval=30,
                        MaxStepsWithoutProgress=50_000,
                        Method=:borg_moea,
                        ϵ=0.05,
                        FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=false))

            res = bboptimize(opt)
            bc_res = best_candidate(res)
            state_state_Qe[ensemble_app] = build_ensemble_flow(Qb, Qm, bc_res, active_params)
            state_state_mix[ensemble_app] = bc_res
        end
    end
end


if !isfile(joinpath(ensemble_res_path, "overall_baseline_NmKGE_results.csv"))
    # Collate results using a single identified state series (e.g., by state-based NmKGE)
    baseline_NmKGE_results = DataFrame(Approach=[], Overall_NmKGE=[],
                                       Overall_calib=[], Wet_calib=[], Mid_calib=[], Dry_calib=[],
                                       Overall_valid=[], Wet_valid=[], Mid_valid=[], Dry_valid=[])
    base_base_NmKGE_results = deepcopy(baseline_NmKGE_results)
    base_state_NmKGE_results = deepcopy(baseline_NmKGE_results)
    state_NmKGE_results = deepcopy(baseline_NmKGE_results)
    state_state_NmKGE_results = deepcopy(baseline_NmKGE_results)

    baseline_A2k_results = DataFrame(Approach=[], Overall_A2k=[], Mean_A2k=[],
                                     Overall_calib=[], Wet_calib=[], Mid_calib=[], Dry_calib=[],
                                     Overall_valid=[], Wet_valid=[], Mid_valid=[], Dry_valid=[])
    base_base_A2k_results = deepcopy(baseline_A2k_results)
    base_state_A2k_results = deepcopy(baseline_A2k_results)
    state_A2k_results = deepcopy(baseline_A2k_results)
    state_state_A2k_results = deepcopy(baseline_A2k_results)

    # active_params = state_active_params["NmKGE"]
    # calib_low, calib_mid, calib_high = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
    # valid_low, valid_mid, valid_high = separate_state_periods(active_params, CALIB_LN+1)
    common_cmd = log.(base_cmd["NmKGE"])
    common_cmd_calib = common_cmd[1:calib_ts]
    common_thresholds = quantile(common_cmd_calib, quantiles[2:end])
    calib_low = findall(common_cmd .< common_thresholds[1])
    calib_mid = findall(common_thresholds[1] .<= common_cmd .< common_thresholds[2])
    calib_high = findall(common_thresholds[2] .<= common_cmd)

    active_params = copy(common_cmd)
    active_params[calib_low] .= 1
    active_params[calib_mid] .= 2
    active_params[calib_high] .= 3

    calib_low, calib_mid, calib_high = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
    valid_low, valid_mid, valid_high = separate_state_periods(active_params, CALIB_LN+1)

    for app in APPROACHES
        # Assess baseline instances
        Qs = base_Q[app]
        insert_NmKGE_scores(baseline_NmKGE_results, app, Qo, Qs)
        insert_A2k_scores(baseline_A2k_results, app, Qo, Qs)

        # State-based instances
        Qs = state_Qe[app]
        insert_NmKGE_scores(state_NmKGE_results, app, Qo, Qs)
        insert_A2k_scores(state_A2k_results, app, Qo, Qs)

        for s_app in APPROACHES
            approach = (app, s_app)

            # base-state instances
            Qs = base_state_Qe[approach]
            insert_NmKGE_scores(base_state_NmKGE_results, approach, Qo, Qs)
            insert_A2k_scores(base_state_A2k_results, approach, Qo, Qs)

            if app == s_app
                # skip non-unique combinations
                continue
            end

            # base-base instances
            Qs = base_base_Q[approach]
            insert_NmKGE_scores(base_base_NmKGE_results, approach, Qo, Qs)
            insert_A2k_scores(base_base_A2k_results, approach, Qo, Qs)

            # State-State ensembles
            Qs = state_state_Qe[approach]
            insert_NmKGE_scores(state_state_NmKGE_results, approach, Qo, Qs)
            insert_A2k_scores(state_state_A2k_results, approach, Qo, Qs)
        end
    end

    # Export to file
    CSV.write(joinpath(ensemble_res_path, "overall_baseline_NmKGE_results.csv"), baseline_NmKGE_results)
    CSV.write(joinpath(ensemble_res_path, "overall_state_NmKGE_results.csv"), state_NmKGE_results)
    CSV.write(joinpath(ensemble_res_path, "overall_base-base_NmKGE_results.csv"), base_base_NmKGE_results)
    CSV.write(joinpath(ensemble_res_path, "overall_base-state_NmKGE_results.csv"), base_state_NmKGE_results)
    CSV.write(joinpath(ensemble_res_path, "overall_state-state_NmKGE_results.csv"), state_state_NmKGE_results)

    CSV.write(joinpath(ensemble_res_path, "overall_baseline_A2k_results.csv"), baseline_A2k_results)
    CSV.write(joinpath(ensemble_res_path, "overall_state_A2k_results.csv"), state_A2k_results)
    CSV.write(joinpath(ensemble_res_path, "overall_base-base_A2k_results.csv"), base_base_A2k_results)
    CSV.write(joinpath(ensemble_res_path, "overall_base-state_A2k_results.csv"), base_state_A2k_results)
    CSV.write(joinpath(ensemble_res_path, "overall_state-state_A2k_results.csv"), state_state_A2k_results)

    # Save outflows
    CSV.write(joinpath(ensemble_res_path, "baseline_outflows.csv"), DataFrame(base_Q))
    CSV.write(joinpath(ensemble_res_path, "state_outflows.csv"), DataFrame(state_Qe))
    CSV.write(joinpath(ensemble_res_path, "base_base_outflows.csv"),
                DataFrame([k=>v for (k,v) in zip(["$(b) - $(m)"
                for (b, m) in keys(base_base_Q)], values(base_base_Q))]))
    CSV.write(joinpath(ensemble_res_path, "base_state_outflows.csv"),
                DataFrame([k=>v for (k,v) in zip(["$(b) - $(m)"
                for (b, m) in keys(base_state_Qe)], values(base_state_Qe))]))
    CSV.write(joinpath(ensemble_res_path, "state_state_outflows.csv"),
                DataFrame([k=>v for (k,v) in zip(["$(b) - $(m)"
                for (b, m) in keys(state_state_Qe)], values(state_state_Qe))]))

    # Save instance weightings
    CSV.write(joinpath(ensemble_res_path, "base_base_mix.csv"),
                DataFrame([k=>v for (k,v) in zip(["$(b) - $(m)"
                for (b, m) in keys(base_base_mix)], values(base_base_mix))]))
    CSV.write(joinpath(ensemble_res_path, "base_state_mix.csv"),
                DataFrame([k=>v for (k,v) in zip(["$(b) - $(m)"
                for (b, m) in keys(base_state_mix)], values(base_state_mix))]))
    CSV.write(joinpath(ensemble_res_path, "state_state_mix.csv"),
                DataFrame([k=>v for (k,v) in zip(["$(b) - $(m)"
                for (b, m) in keys(state_state_mix)], values(state_state_mix))]))

    CSV.write(joinpath(ensemble_res_path, "state_active_params.csv"),
                DataFrame(state_active_params))

    CSV.write(joinpath(ensemble_res_path, "base_NmKGE_state.csv"), DataFrame("NmKGE"=>active_params))
end
