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


function combine_ensemble(x, Qb, Qm, active_params)
    # Calibration - scores for given approach
    ensemble_flow = deepcopy(Qb)

    low_factor, mid_factor, high_factor = x

    # Create ensemble outflow for entire simulation period
    low_store, mid_store, high_store = separate_state_periods(active_params)
    ensemble_flow[high_store] = (Qm[high_store] .* high_factor) .+ (Qb[high_store] .* (1.0 - high_factor))
    ensemble_flow[mid_store] = (Qm[mid_store] .* mid_factor) .+ (Qb[mid_store] .* (1.0 - mid_factor))
    ensemble_flow[low_store] = (Qm[low_store] .* low_factor) .+ (Qb[low_store] .* (1.0 - low_factor))

    # Subset to calibration period only
    calib_low_store, calib_mid_store, calib_high_store = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
    ensemble_calib = ensemble_flow[calib_ts+1:CALIB_LN]
    overall = Streamfall.mKGE(Qo_calib, ensemble_calib)
    low_score = Streamfall.mKGE(Qo_calib[calib_low_store], ensemble_calib[calib_low_store])
    mid_score = Streamfall.mKGE(Qo_calib[calib_mid_store], ensemble_calib[calib_mid_store])
    high_score = Streamfall.mKGE(Qo_calib[calib_high_store], ensemble_calib[calib_high_store])

    # ensemble_valid = ensemble_flow[CALIB_LN+1:end]
    # valid_low_store, valid_mid_store, valid_high_store = separate_state_periods(active_params, CALIB_LN)
    # high_score = Streamfall.mKGE(Qo_valid[valid_high_store], ensemble_valid[valid_high_store])
    # if length(valid_high_store) == 0
    #     return (-9999.9, -9999.9, -9999.9, -9999.9)
    # end

    return (overall, low_score, mid_score, high_score)
end


"""
Determine best outflow mix for a base-base ensemble.
"""
function base_base_mix(x; model1="split_mean_NmKGE", model2="NmKGE")
    local active_params = state_active_params[model1]
    local Qb = base_outflows[model1]
    local Qm = base_outflows[model2]

    local ensemble_flow = (Qm .* x[1]) .+ (Qb .* (1.0 - x[1]))

    # ensemble_calib = ensemble_flow[calib_ts+1:CALIB_LN]
    # score = 1.0 - Streamfall.mKGE(Qo_calib, ensemble_calib)
    # return score
    return combine_ensemble(repeat(x, 3), Qb, Qm, active_params)
end


"""
Determine best outflow mix for a base-state ensemble.
"""
function base_state_mix(x; model1="split_mean_NmKGE", model2="NmKGE")
    local active_params = state_active_params[model2]  # state_active_params[model2]  # state_active_params["NmKGE"]
    local Qb = base_outflows[model1]
    local Qm = state_outflows[model2]

    return combine_ensemble(x, Qb, Qm, active_params)
end


"""
Determine best pair of ensembles
"""
function ensemble_mix(x; model1="split_mean_NmKGE", model2="NmKGE")
    local active_params = state_active_params[model1]  # state_active_params[model1]  # state_active_params["NmKGE"]
    local Qb = state_outflows[model1]
    local Qm = state_outflows[model2]

    return combine_ensemble(x, Qb, Qm, active_params)
end


calib_ts = 3651
state = :storage
state_str = String(state)

target_idx = [1,2,5,6,7,8]
quantiles = [0.0, 0.1, 0.9]

Qo = FULL_DATASET[:, "410730_Q"]
Qo_calib = CALIB_OBS[calib_ts+1:CALIB_LN]
Qo_valid = VALID_OBS

if !@isdefined state_active_params
    base_outflows = Dict()
    state_outflows = Dict()
    state_active_params = Dict()

    for app in APPROACHES
        # baseline
        base_sn, base_id = setup_network("baselines/cotter_baseline_IHACRES_$(app).yml")
        base_node = base_sn[base_id]
        Streamfall.run_basin!(base_sn, FULL_CLIMATE)

        base_outflows[app] = base_node.outflow

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

        state_outflows[app] = node.outflow
        state_active_params[app] = active_params
    end
end


function calc_metrics!(primary, secondary, Qb, Qm, active_params, res_score, split_params, AD_test::DataFrame, mKGE_test::DataFrame, ensembles::Dict, mix::Dict)
    calc_metrics!(primary, secondary, Qb, Qm, active_params, split_params, AD_test, mKGE_test, ensembles, mix)

    local x = Array(mKGE_test[nrow(mKGE_test), [:Overall_mKGE_calib, :Wet_calib, :Mid_calib, :Dry_calib]])

    try
        if res_score isa Array
            @assert all(res_score .== x)
        end
    catch e
        if e isa AssertionError
            println("Assertion test failed! Score does not match recalculated scores")
            @infiltrate
        end
    end
end


"""
Calculate all metrics for a given model pair, assigning results to AD_test, mKGE_test and ensembles

mix : Dictionary holding weights used to combine model instances
"""
function calc_metrics!(primary, secondary, Qb, Qm, active_params, split_params, AD_test::DataFrame, mKGE_test::DataFrame, ensembles::Dict, mix::Dict)

    mix[(primary, secondary)] = split_params

    # Calibration - scores for given approach
    ensemble_flow = deepcopy(Qb)

    # Weight towards whichever has the better score in calibration period
    low_factor, mid_factor, high_factor = split_params

    # Create ensemble outflow for entire simulation period
    low_store, mid_store, high_store = separate_state_periods(active_params)
    ensemble_flow[high_store] = (Qm[high_store] .* high_factor) .+ (Qb[high_store] .* (1.0 - high_factor))
    ensemble_flow[mid_store] = (Qm[mid_store] .* mid_factor) .+ (Qb[mid_store] .* (1.0 - mid_factor))
    ensemble_flow[low_store] = (Qm[low_store] .* low_factor) .+ (Qb[low_store] .* (1.0 - low_factor))

    # Subset to calibration period only
    c_params = active_params[calib_ts+1:CALIB_LN]
    calib_low_store, calib_mid_store, calib_high_store = separate_state_periods(active_params, calib_ts+1, CALIB_LN)

    # Get metrics for calibration/validation periods
    ensemble_calib = ensemble_flow[calib_ts+1:CALIB_LN]
    AD_high = KSampleADTest(Qo_calib[calib_high_store], ensemble_calib[calib_high_store]; modified=true).A²k
    AD_mid = KSampleADTest(Qo_calib[calib_mid_store], ensemble_calib[calib_mid_store]; modified=true).A²k
    AD_low = KSampleADTest(Qo_calib[calib_low_store], ensemble_calib[calib_low_store]; modified=true).A²k

    low_KGE = Streamfall.mKGE(Qo_calib[calib_low_store], ensemble_calib[calib_low_store])
    mid_KGE = Streamfall.mKGE(Qo_calib[calib_mid_store], ensemble_calib[calib_mid_store])
    high_KGE = Streamfall.mKGE(Qo_calib[calib_high_store], ensemble_calib[calib_high_store])

    valid_low_store, valid_mid_store, valid_high_store = separate_state_periods(active_params, CALIB_LN+1)
    ensemble_valid = ensemble_flow[CALIB_LN+1:end]

    if length(valid_high_store) == 0
        v_AD_high = missing
    else
        v_AD_high = KSampleADTest(Qo_valid[valid_high_store], ensemble_valid[valid_high_store]; modified=true).A²k
    end

    if length(valid_mid_store) == 0
        v_AD_mid = missing
    else
        v_AD_mid = KSampleADTest(Qo_valid[valid_mid_store], ensemble_valid[valid_mid_store]; modified=true).A²k
    end

    if length(valid_low_store) == 0
        v_AD_low = missing
    else
        v_AD_low = KSampleADTest(Qo_valid[valid_low_store], ensemble_valid[valid_low_store]; modified=true).A²k
    end

    if length(valid_high_store) == 0
        high = missing
    else
        high = Streamfall.mKGE(Qo_valid[valid_high_store], ensemble_valid[valid_high_store])
    end

    if length(valid_mid_store) == 0
        mid = missing
    else
        mid = Streamfall.mKGE(Qo_valid[valid_mid_store], ensemble_valid[valid_mid_store])
    end

    if length(valid_low_store) == 0
        low = missing
    else
        low = Streamfall.mKGE(Qo_valid[valid_low_store], ensemble_valid[valid_low_store])
    end

    push!(AD_test,
        [primary, secondary, mean([AD_low, AD_mid, AD_high]),
         KSampleADTest(Qo_calib, ensemble_calib; modified=true).A²k, AD_low, AD_mid, AD_high,
         KSampleADTest(Qo_valid, ensemble_valid; modified=true).A²k, v_AD_low, v_AD_mid, v_AD_high
        ])

    push!(mKGE_test,
          [primary, secondary,
           Streamfall.mKGE(Qo, ensemble_flow),
           Streamfall.mKGE(Qo_calib, ensemble_calib), Streamfall.mKGE(Qo_valid, ensemble_valid),
           mean(skipmissing([high_KGE, mid_KGE, low_KGE])), low_KGE, mid_KGE, high_KGE,
           mean(skipmissing([low, mid, high])), low, mid, high
        ]
    )

    ensembles[(primary, secondary)] = ensemble_flow
end

if !@isdefined AD_test_results
    AD_test_results = DataFrame(Primary=[], Secondary=[], Mean_A2k=[],
                                Overall_calib=[], Wet_calib=[], Mid_calib=[], Dry_calib=[],
                                Overall_valid=[], Wet_valid=[], Mid_valid=[], Dry_valid=[])
    mKGE_test_results = DataFrame(Primary=[], Secondary=[],
                                Overall_mKGE=[],
                                Overall_mKGE_calib=[], Overall_mKGE_valid=[],
                                mean_mKGE_calib=[], Wet_calib=[], Mid_calib=[], Dry_calib=[],
                                mean_mKGE_valid=[], Wet_valid=[], Mid_valid=[], Dry_valid=[]
                                )
    base_AD_test_results = deepcopy(AD_test_results)
    base_mKGE_test_results = deepcopy(mKGE_test_results)

    bb_AD_test_results = deepcopy(AD_test_results)
    bb_mKGE_test_results = deepcopy(mKGE_test_results)
    base_base_results = DataFrame(Primary=[], Secondary=[], calib_mKGE=[], mix=[])
    bb_ensembles = Dict()
    bb_param_mix = Dict()

    base_state_ensembles = Dict()
    ensembles = Dict()

    base_ensemble_param_mix = Dict()
    ensemble_param_mix = Dict()

    # Optimize ensembles and build comparison DFs
    for app in APPROACHES
        for state_approach in APPROACHES
            # Compare baseline-state_based combinations
            local active_params = state_active_params[app]
            local Qb = base_outflows[app]
            local Qm = state_outflows[state_approach]

            opt_func = (x) -> base_state_mix(x; model1=app, model2=state_approach)
            opt = bbsetup(opt_func; SearchRange=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
                        MaxTime=90, TraceInterval=30,
                        MaxStepsWithoutProgress=100_000,
                        Method=:borg_moea,
                        ϵ=0.01,
                        FitnessScheme=ParetoFitnessScheme{4}(is_minimizing=false))

            res = bboptimize(opt)
            calc_metrics!(app, state_approach, Qb, Qm, active_params, best_fitness(res), best_candidate(res),
                          base_AD_test_results, base_mKGE_test_results, base_state_ensembles, base_ensemble_param_mix)

            # Compare pairs built with state-based models only
            # Only compare unique combinations
            if app == state_approach
                continue
            end

            # Create base-base ensemble
            opt_func = (x) -> base_base_mix(x; model1=app, model2=state_approach)
            opt = bbsetup(opt_func; SearchRange=[(0.0, 1.0)],
                        MaxTime=90, TraceInterval=30,
                        MaxStepsWithoutProgress=100_000,
                        Method=:borg_moea,
                        ϵ=0.01,
                        FitnessScheme=ParetoFitnessScheme{4}(is_minimizing=false))
            res = bboptimize(opt)
            mixer = best_candidate(res)[1]
            push!(base_base_results, [app, state_approach, best_fitness(res), mixer])

            local Qb = base_outflows[app]
            local Qm = base_outflows[state_approach]
            calc_metrics!(app, state_approach, Qb, Qm, active_params, best_fitness(res), [mixer, mixer, mixer],
                          bb_AD_test_results, bb_mKGE_test_results, bb_ensembles, bb_param_mix)

            # Build state-based pairs
            opt_func = (x) -> ensemble_mix(x; model1=app, model2=state_approach)
            opt = bbsetup(opt_func; SearchRange=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
                        MaxTime=90, TraceInterval=30,
                        MaxStepsWithoutProgress=100_000,
                        Method=:borg_moea,
                        ϵ=0.01,
                        FitnessScheme=ParetoFitnessScheme{4}(is_minimizing=false))

            res = bboptimize(opt)

            local Qb = state_outflows[app]
            local Qm = state_outflows[state_approach]
            calc_metrics!(app, state_approach, Qb, Qm, active_params, best_fitness(res), best_candidate(res),
                          AD_test_results, mKGE_test_results, ensembles, ensemble_param_mix)
        end
    end
end


function insert_mKGE_scores(df, approach, Qo, Qm)
    local Qo_calib = Qo[calib_ts+1:CALIB_LN]
    local Qo_valid = Qo[CALIB_LN+1:end]
    local Qm_calib = Qm[calib_ts+1:CALIB_LN]
    local Qm_valid = Qm[CALIB_LN+1:end]

    push!(df, [approach,
        Streamfall.NmKGE(Qo, Qm),

        Streamfall.NmKGE(Qo_calib, Qm_calib),
        Streamfall.NmKGE(Qo_calib[calib_low], Qm_calib[calib_low]),
        Streamfall.NmKGE(Qo_calib[calib_mid], Qm_calib[calib_mid]),
        Streamfall.NmKGE(Qo_calib[calib_high], Qm_calib[calib_high]),

        Streamfall.NmKGE(Qo_valid, Qm_valid),
        Streamfall.NmKGE(Qo_valid[valid_low], Qm_valid[valid_low]),
        Streamfall.NmKGE(Qo_valid[valid_mid], Qm_valid[valid_mid]),
        Streamfall.NmKGE(Qo_valid[valid_high], Qm_valid[valid_high]),
    ])
end

function insert_A2k_scores(df, approach, Qo, Qm)
    local Qm_calib = Qm[calib_ts+1:CALIB_LN]
    local Qm_valid = Qm[CALIB_LN+1:end]

    push!(df, [approach,
            KSampleADTest(Qo, Qm).A²k,

            KSampleADTest(Qo_calib, Qm_calib).A²k,
            KSampleADTest(Qo_calib[calib_low], Qm_calib[calib_low]).A²k,
            KSampleADTest(Qo_calib[calib_mid], Qm_calib[calib_mid]).A²k,
            KSampleADTest(Qo_calib[calib_high], Qm_calib[calib_high]).A²k,

            KSampleADTest(Qo_valid, Qm_valid).A²k,
            KSampleADTest(Qo_valid[valid_low], Qm_valid[valid_low]).A²k,
            KSampleADTest(Qo_valid[valid_mid], Qm_valid[valid_mid]).A²k,
            KSampleADTest(Qo_valid[valid_high], Qm_valid[valid_high]).A²k,
    ])
end


# Collate results using a single identified state series (e.g., by state-based NmKGE)
baseline_results = DataFrame(Approach=[], Overall_mKGE=[],
                                Overall_calib=[], Wet_calib=[], Mid_calib=[], Dry_calib=[],
                                Overall_valid=[], Wet_valid=[], Mid_valid=[], Dry_valid=[])

baseline_A2k_results = DataFrame(Approach=[], Mean_A2k=[],
                                        Overall_calib=[], Wet_calib=[], Mid_calib=[], Dry_calib=[],
                                        Overall_valid=[], Wet_valid=[], Mid_valid=[], Dry_valid=[])

basebase_results = deepcopy(baseline_results)
base_state_results = deepcopy(baseline_results)
state_based_results = deepcopy(baseline_results)
ensemble_results = deepcopy(baseline_results)

basebase_A2k_results = deepcopy(baseline_A2k_results)
base_state_A2k_results = deepcopy(baseline_A2k_results)
state_based_A2k_results = deepcopy(baseline_A2k_results)
ensemble_A2k_results = deepcopy(baseline_A2k_results)

active_params = state_active_params["NmKGE"]
calib_low, calib_mid, calib_high = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
valid_low, valid_mid, valid_high = separate_state_periods(active_params, CALIB_LN+1)
for app in APPROACHES

    # Assess baseline instances
    Qs = base_outflows[app]
    insert_mKGE_scores(baseline_results, app, Qo, Qs)
    insert_A2k_scores(baseline_A2k_results, app, Qo, Qs)

    # State-based instances
    Qs = state_outflows[app]
    insert_mKGE_scores(state_based_results, app, Qo, Qs)
    insert_A2k_scores(state_based_A2k_results, app, Qo, Qs)
end

for app in APPROACHES
    for s_app in APPROACHES

        approach = (app, s_app)
        # base-state instances
        Qs = base_state_ensembles[approach]
        insert_mKGE_scores(base_state_results, approach, Qo, Qs)
        insert_A2k_scores(base_state_A2k_results, approach, Qo, Qs)

        if app == s_app
            continue
        end

        # base-base instances
        Qs = bb_ensembles[(app, s_app)]
        insert_mKGE_scores(basebase_results, approach, Qo, Qs)
        insert_A2k_scores(basebase_A2k_results, approach, Qo, Qs)

        # State-State ensembles
        Qs = ensembles[(app, s_app)]
        insert_mKGE_scores(ensemble_results, approach, Qo, Qs)
        insert_A2k_scores(ensemble_A2k_results, approach, Qo, Qs)
    end
end

# Export to file
# sort!(baseline_results, :Overall_calib, rev=true)
CSV.write(joinpath(ensemble_res_path, "overall_baseline_results.csv"), baseline_results)
CSV.write(joinpath(ensemble_res_path, "overall_state-based_results.csv"), state_based_results)
CSV.write(joinpath(ensemble_res_path, "overall_base-base_results.csv"), basebase_results)
CSV.write(joinpath(ensemble_res_path, "overall_base-state_results.csv"), base_state_results)
CSV.write(joinpath(ensemble_res_path, "overall_state-state_results.csv"), ensemble_results)

CSV.write(joinpath(ensemble_res_path, "overall_baseline_A2k_results.csv"), baseline_A2k_results)
CSV.write(joinpath(ensemble_res_path, "overall_state-based_A2k_results.csv"), state_based_A2k_results)
CSV.write(joinpath(ensemble_res_path, "overall_base-base_A2k_results.csv"), basebase_A2k_results)
CSV.write(joinpath(ensemble_res_path, "overall_base-state_A2k_results.csv"), base_state_A2k_results)
CSV.write(joinpath(ensemble_res_path, "overall_state-state_A2k_results.csv"), ensemble_A2k_results)

calib_bl = baseline_results[:, [:Overall_calib, :Wet_calib, :Mid_calib, :Dry_calib]]
valid_bl = baseline_results[:, [:Overall_valid, :Wet_valid, :Mid_valid, :Dry_valid]]
rename!(calib_bl, Dict(
    "Overall_calib"=>"Overall",
    "Wet_calib"=>"Wet",
    "Mid_calib"=>"Mid",
    "Dry_calib"=>"Dry",
))

rename!(valid_bl, Dict(
    "Overall_valid"=>"Overall",
    "Wet_valid"=>"Wet",
    "Mid_valid"=>"Mid",
    "Dry_valid"=>"Dry",
))

# State-based
calib_en = state_based_results[:, [:Overall_calib, :Wet_calib, :Mid_calib, :Dry_calib]]
valid_en = state_based_results[:, [:Overall_valid, :Wet_valid, :Mid_valid, :Dry_valid]]
rename!(calib_en, Dict(
    "Overall_calib"=>"Overall",
    "Wet_calib"=>"Wet",
    "Mid_calib"=>"Mid",
    "Dry_calib"=>"Dry",
))

rename!(valid_en, Dict(
    "Overall_valid"=>"Overall",
    "Wet_valid"=>"Wet",
    "Mid_valid"=>"Mid",
    "Dry_valid"=>"Dry",
))


baseline_vplot = violin(title="Baseline", ylabel="NKGE'", label=["Calibration", "Validation"], legend=:bottomleft)
for (idx, n) in enumerate(names(calib_bl))
    lab = idx > 1 ? "" : "Calibration"
    violin!([n], calib_bl[:, n], side=:left, color="blue", label=lab)

    lab = idx > 1 ? "" : "Validation"
    violin!([n], valid_bl[:, n], side=:right, color="orange", label=lab)
end

savefig(baseline_vplot, "$(ensemble_figs)baseline_NmKGE_results.png")



state_based_vplot = violin(title="State-based", ylabel="NKGE'", label=["Calibration", "Validation"], legend=:bottomleft)
for (idx, n) in enumerate(names(calib_en))
    lab = idx > 1 ? "" : "Calibration"
    violin!([n], calib_en[:, n], side=:left, color="blue", label=lab)

    lab = idx > 1 ? "" : "Validation"
    violin!([n], valid_en[:, n], side=:right, color="orange", label=lab)
end

savefig(state_based_vplot, "$(ensemble_figs)state-based_NmKGE_results.png")


# Base-Base Ensemble violin plot (display normalized values)
bb_ensemble_vplot = violin(title="Base-Base Ensemble", ylabel="NKGE'", label=["Calibration", "Validation"], legend=:bottomleft)
calib_en = basebase_results[:, [:Overall_calib, :Wet_calib, :Mid_calib, :Dry_calib]]
valid_en = basebase_results[:, [:Overall_valid, :Wet_valid, :Mid_valid, :Dry_valid]]

rename!(calib_en, Dict(
    "Overall_calib"=>"Overall",
    "Wet_calib"=>"Wet",
    "Mid_calib"=>"Mid",
    "Dry_calib"=>"Dry",
))

rename!(valid_en, Dict(
    "Overall_valid"=>"Overall",
    "Wet_valid"=>"Wet",
    "Mid_valid"=>"Mid",
    "Dry_valid"=>"Dry",
))

for (idx, n) in enumerate(names(calib_en))
    lab = idx > 1 ? "" : "Calibration"
    violin!([n], 1 ./ (2 .- calib_en[:, n]), side=:left, color="blue", label=lab)

    lab = idx > 1 ? "" : "Validation"
    tmp = collect(skipmissing(valid_en[:, n]))
    violin!([n], 1 ./ (2 .- tmp), side=:right, color="orange", label=lab)
end

savefig(bb_ensemble_vplot, "$(ensemble_figs)base_base_ensemble_NmKGE_results.png")


# Baseline-State-based Ensemble violin plot (display normalized values)
base_state_ensemble_vplot = violin(title="Base-State Ensemble", ylabel="NKGE'", label=["Calibration", "Validation"], legend=:bottomleft)
calib_en = base_state_results[:, [:Overall_calib, :Wet_calib, :Mid_calib, :Dry_calib]]
valid_en = base_state_results[:, [:Overall_valid, :Wet_valid, :Mid_valid, :Dry_valid]]

rename!(calib_en, Dict(
    "Overall_calib"=>"Overall",
    "Wet_calib"=>"Wet",
    "Mid_calib"=>"Mid",
    "Dry_calib"=>"Dry",
))

rename!(valid_en, Dict(
    "Overall_valid"=>"Overall",
    "Wet_valid"=>"Wet",
    "Mid_valid"=>"Mid",
    "Dry_valid"=>"Dry",
))

for (idx, n) in enumerate(names(calib_en))
    lab = idx > 1 ? "" : "Calibration"
    violin!([n], 1 ./ (2 .- calib_en[:, n]), side=:left, color="blue", label=lab)

    lab = idx > 1 ? "" : "Validation"
    tmp = collect(skipmissing(valid_en[:, n]))
    violin!([n], 1 ./ (2 .- tmp), side=:right, color="orange", label=lab)
end

savefig(base_state_ensemble_vplot, "$(ensemble_figs)base_state_ensemble_NmKGE_results.png")


# State-State ensemble violin plot (display normalized values)
ensemble_vplot = violin(title="Ensemble", ylabel="NKGE'", label=["Calibration", "Validation"], legend=:bottomleft)
calib_en = ensemble_results[:, [:Overall_calib, :Wet_calib, :Mid_calib, :Dry_calib]]
valid_en = ensemble_results[:, [:Overall_valid, :Wet_valid, :Mid_valid, :Dry_valid]]
rename!(calib_en, Dict(
    "Overall_calib"=>"Overall",
    "Wet_calib"=>"Wet",
    "Mid_calib"=>"Mid",
    "Dry_calib"=>"Dry",
))

rename!(valid_en, Dict(
    "Overall_valid"=>"Overall",
    "Wet_valid"=>"Wet",
    "Mid_valid"=>"Mid",
    "Dry_valid"=>"Dry",
))

for (idx, n) in enumerate(names(calib_en))
    lab = idx > 1 ? "" : "Calibration"
    violin!([n], 1 ./ (2 .- calib_en[:, n]), side=:left, color="blue", label=lab)

    lab = idx > 1 ? "" : "Validation"
    tmp = collect(skipmissing(valid_en[:, n]))
    violin!([n], 1 ./ (2 .- tmp), side=:right, color="orange", label=lab)
end

savefig(ensemble_vplot, "$(ensemble_figs)ensemble_NmKGE_results.png")

combined_violin = plot(
    baseline_vplot,
    plot(state_based_vplot, ylabel="", legend=false),
    plot(ensemble_vplot, ylabel="", legend=false),
    layout=(1, 3),
    margin=5Plots.mm,
    size=(1200,400),
    link=:y
)
savefig(combined_violin, "$(ensemble_figs)combined_violin.png")


combined_violin_all = plot(
    baseline_vplot,
    plot(state_based_vplot, ylabel="", legend=false),
    plot(bb_ensemble_vplot, ylabel="", legend=false),
    plot(base_state_ensemble_vplot, ylabel="", legend=false),
    plot(ensemble_vplot, ylabel="", legend=false),
    layout=(1, 5),
    margin=5Plots.mm,
    size=(2000,400),
    link=:y
)
savefig(combined_violin_all, "$(ensemble_figs)combined_violin_all.png")


begin
    constituent_scores = DataFrame(
        :model => [],
        :overall_calib => [],
        :wet_calib => [],
        :usual_calib => [],
        :dry_calib => [],
        :overall_valid => [],
        :wet_valid => [],
        :usual_valid => [],
        :dry_valid => [],
    )

    en_AD_mod = filter(row -> !(ismissing(row.Dry_valid)),  AD_test_results)
    sort!(en_AD_mod, :Mean_A2k)
    sort!(AD_test_results, :Mean_A2k)
    sort!(mKGE_test_results, order(:mean_mKGE_valid, rev=true))

    base_state_AD_mod = filter(row -> !(ismissing(row.Dry_valid)),  base_AD_test_results)
    sort!(base_state_AD_mod, :Mean_A2k)
    sort!(base_mKGE_test_results, order(:mean_mKGE_valid, rev=true))

    bb_AD_mod = filter(row -> !(ismissing(row.Dry_valid)),  bb_AD_test_results)
    sort!(bb_AD_mod, :Mean_A2k)
    sort!(bb_mKGE_test_results, order(:mean_mKGE_valid, rev=true))

    model1 = en_AD_mod[1, :Primary]
    model2 = en_AD_mod[1, :Secondary]
    score = en_AD_mod[1, :Mean_A2k]

    local active_params = state_active_params["NmKGE"]

    # Model 1
    s_calib = state_outflows[model1][calib_ts+1:CALIB_LN]
    s_valid = state_outflows[model1][CALIB_LN+1:end]
    calib_period = Streamfall.mKGE(Qo_calib, s_calib)
    valid_period = Streamfall.mKGE(Qo_valid, s_valid)

    calib_low_store, calib_mid_store, calib_high_store = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
    sc_high = Streamfall.mKGE(Qo_calib[calib_high_store], s_calib[calib_high_store])
    sc_mid = Streamfall.mKGE(Qo_calib[calib_mid_store], s_calib[calib_mid_store])
    sc_low = Streamfall.mKGE(Qo_calib[calib_low_store], s_calib[calib_low_store])

    valid_low_store, valid_mid_store, valid_high_store = separate_state_periods(active_params, CALIB_LN+1)
    sv_high = Streamfall.mKGE(Qo_valid[valid_high_store], s_valid[valid_high_store])
    sv_mid = Streamfall.mKGE(Qo_valid[valid_mid_store], s_valid[valid_mid_store])
    sv_low = Streamfall.mKGE(Qo_valid[valid_low_store], s_valid[valid_low_store])

    push!(
        constituent_scores,
        ["State-based $model1", calib_period, sc_low, sc_mid, sc_high, valid_period, sv_low, sv_mid, sv_high]
    )

    # Model 2
    s_calib = state_outflows[model2][calib_ts+1:CALIB_LN]
    s_valid = state_outflows[model2][CALIB_LN+1:end]
    calib_period = Streamfall.mKGE(Qo_calib, s_calib)
    valid_period = Streamfall.mKGE(Qo_valid, s_valid)

    # Use same temporal locations as model 1
    sc_high = Streamfall.mKGE(Qo_calib[calib_high_store], s_calib[calib_high_store])
    sc_mid = Streamfall.mKGE(Qo_calib[calib_mid_store], s_calib[calib_mid_store])
    sc_low = Streamfall.mKGE(Qo_calib[calib_low_store], s_calib[calib_low_store])

    sv_high = Streamfall.mKGE(Qo_valid[valid_high_store], s_valid[valid_high_store])
    sv_mid = Streamfall.mKGE(Qo_valid[valid_mid_store], s_valid[valid_mid_store])
    sv_low = Streamfall.mKGE(Qo_valid[valid_low_store], s_valid[valid_low_store])
    push!(
        constituent_scores,
        ["State-based $model2", calib_period, sc_low, sc_mid, sc_high, valid_period, sv_low, sv_mid, sv_high]
    )

    # Ensembles
    best_ensemble = (model1, model2)
    flow = ensembles[best_ensemble]
    flow_calib = flow[calib_ts+1:CALIB_LN]
    flow_valid = flow[CALIB_LN+1:end]
    calib_period = Streamfall.mKGE(Qo_calib, flow_calib)
    valid_period = Streamfall.mKGE(Qo_valid, flow_valid)

    @info "Best Performing combination:" model1 model2 score calib_period valid_period

    fc_overall = Streamfall.mKGE(Qo_calib, flow_calib)
    fc_high = Streamfall.mKGE(Qo_calib[calib_high_store], flow_calib[calib_high_store])
    fc_mid = Streamfall.mKGE(Qo_calib[calib_mid_store], flow_calib[calib_mid_store])
    fc_low = Streamfall.mKGE(Qo_calib[calib_low_store], flow_calib[calib_low_store])
    @info "Calibration Low | Mid | High" fc_low fc_mid fc_high

    @info "high/dry store" length(valid_high_store)
    @info "mid store" length(valid_mid_store)
    @info "low/wet store" length(valid_low_store)

    Qbase_valid = base_outflows[model1][CALIB_LN+1:end]
    fv_overall = Streamfall.mKGE(Qo_valid, flow_valid)
    bv_overall = Streamfall.mKGE(Qo_valid, Qbase_valid)

    if length(valid_high_store) > 0
        fv_high = Streamfall.mKGE(Qo_valid[valid_high_store], flow_valid[valid_high_store])
        bv_high = Streamfall.mKGE(Qo_valid[valid_high_store], Qbase_valid[valid_high_store])
    else
        fv_high = missing
        bv_high = missing
    end

    if length(valid_mid_store) > 0
        fv_mid = Streamfall.mKGE(Qo_valid[valid_mid_store], flow_valid[valid_mid_store])
        bv_mid = Streamfall.mKGE(Qo_valid[valid_mid_store], Qbase_valid[valid_mid_store])
    else
        fv_mid = missing
        bv_mid = missing
    end

    if length(valid_low_store) > 0
        fv_low = Streamfall.mKGE(Qo_valid[valid_low_store], flow_valid[valid_low_store])
        bv_low = Streamfall.mKGE(Qo_valid[valid_low_store], Qbase_valid[valid_low_store])
    else
        fv_low = missing
        bv_low = missing
    end

    @info "Validation Low | Mid | High" fv_low fv_mid fv_high
    push!(
        constituent_scores,
        ["Ensemble", fc_overall, fc_low, fc_mid, fc_high, fv_overall, fv_low, fv_mid, fv_high]
    )

    Qbase_calib = base_outflows[model1][calib_ts+1:CALIB_LN]
    bc_overall = Streamfall.mKGE(Qo_calib, Qbase_calib)
    bc_high = Streamfall.mKGE(Qo_calib[calib_high_store], Qbase_calib[calib_high_store])
    bc_mid = Streamfall.mKGE(Qo_calib[calib_mid_store], Qbase_calib[calib_mid_store])
    bc_low = Streamfall.mKGE(Qo_calib[calib_low_store], Qbase_calib[calib_low_store])
    @info "Baseline 1 calib Low | Mid | High" bc_low bc_mid bc_high
    @info "Baseline 1 validation Low | Mid | High" bv_low bv_mid bv_high
    @info "Baseline 1 overall: Calib | Valid" bc_overall bv_overall

    push!(
        constituent_scores,
        ["Baseline $model1", bc_overall, bc_low, bc_mid, bc_high, bv_overall, bv_low, bv_mid, bv_high]
    )

    Qbase_calib = base_outflows[model2][calib_ts+1:CALIB_LN]
    bc_overall = Streamfall.mKGE(Qo_calib, Qbase_calib)
    bc_high = Streamfall.mKGE(Qo_calib[calib_high_store], Qbase_calib[calib_high_store])
    bc_mid = Streamfall.mKGE(Qo_calib[calib_mid_store], Qbase_calib[calib_mid_store])
    bc_low = Streamfall.mKGE(Qo_calib[calib_low_store], Qbase_calib[calib_low_store])

    Qbase_valid = base_outflows[model2][CALIB_LN+1:end]
    bv_overall = Streamfall.mKGE(Qo_valid, Qbase_valid)
    bv_high = Streamfall.mKGE(Qo_valid[valid_high_store], Qbase_valid[valid_high_store])
    bv_mid = Streamfall.mKGE(Qo_valid[valid_mid_store], Qbase_valid[valid_mid_store])
    bv_low = Streamfall.mKGE(Qo_valid[valid_low_store], Qbase_valid[valid_low_store])
    @info "Baseline 2 calib Low | Mid | High" bc_low bc_mid bc_high

    @info "Baseline 2 validation Low | Mid | High" bv_low bv_mid bv_high

    @info "Baseline 2 overall: Calib | Valid" bc_overall bv_overall

    push!(
        constituent_scores,
        ["Baseline $model2", bc_overall, bc_low, bc_mid, bc_high, bv_overall, bv_low, bv_mid, bv_high]
    )
end


CSV.write(joinpath(ensemble_res_path, "state_ensemble_test_results_mKGE.csv"), mKGE_test_results)
CSV.write(joinpath(ensemble_res_path, "state_ensemble_test_results_AD.csv"), AD_test_results)

CSV.write(joinpath(ensemble_res_path, "base_state_test_results_mKGE.csv"), base_mKGE_test_results)
CSV.write(joinpath(ensemble_res_path, "base_state_test_results_AD.csv"), base_AD_test_results)

CSV.write(joinpath(ensemble_res_path, "bb_test_results_mKGE.csv"), bb_mKGE_test_results)


function transpose_df(df)
    df_T = DataFrame([[names(df)]; collect.(eachrow(df))], [:column; Symbol.(axes(df, 1))])
    rename!(df_T, Symbol.(Array(df_T[1, :])))
    df_T = df_T[2:end, :]
    return df_T
end

# transpose dataframe
constituent_scores_T = transpose_df(constituent_scores)
CSV.write(joinpath(ensemble_res_path, "constituents_mKGE.csv"), constituent_scores_T)

DATES = FULL_DATASET[calib_ts+1:end, "Date"]
Qburn = Qo[calib_ts+1:end]
Qbase = base_outflows["NmKGE"][calib_ts+1:end]
Qens = ensembles[best_ensemble][calib_ts+1:end]

# plot 2008 to 2009 (high flows are still an issue, but low flow times are better represented)
subperiod = findall(Date(2008) .<= DATES .< Date(2009))
Streamfall.mKGE(Qburn[subperiod], Qbase[subperiod])
Streamfall.mKGE(Qburn[subperiod], Qens[subperiod])


comp_plot = plot(
    DATES[subperiod], Qburn[subperiod],
    legend=:topleft,
    label="Historic",
    margin=5Plots.mm,
    xlabel="Date",
    ylabel="Streamflow [ML]"
)

plot!(DATES[subperiod], Qbase[subperiod], alpha=0.8, label="Baseline NKGE'", linestyle=:dash)
plot!(DATES[subperiod], Qens[subperiod], alpha=0.8, linewidth=1.5, linestyle=:dashdot, label="Ensemble")

savefig(comp_plot, "$(ensemble_figs)comparison_2008-2009_ensemble.png")


# Plot residuals
baseline_residual_subperiod = Qbase[subperiod] - Qburn[subperiod]
mean_residual_baseline_subperiod = round(mean(baseline_residual_subperiod), digits=2)

end_calib = CALIB_LN-calib_ts
c_baseline_residual = Qbase[1:end_calib] - Qburn[1:end_calib]
c_mean_residual_baseline = round(mean(c_baseline_residual), digits=2)

v_baseline_residual = Qbase[end_calib+1:end] - Qburn[end_calib+1:end]
v_mean_residual_baseline = round(mean(v_baseline_residual), digits=2)

ensemble_residual_subperiod = Qens[subperiod] - Qburn[subperiod]
mean_residual_ensemble_subperiod = round(mean(ensemble_residual_subperiod), digits=2)

c_ensemble_residual = Qens[1:end_calib] - Qburn[1:end_calib]
c_mean_residual_ensemble = round(mean(c_ensemble_residual), digits=2)

v_ensemble_residual = Qens[end_calib+1:end] - Qburn[end_calib+1:end]
v_mean_residual_ensemble = round(mean(v_ensemble_residual), digits=2)


res_plot = plot(
    DATES[subperiod], baseline_residual_subperiod,
    legend=:bottomright,
    label="Baseline NKGE' (C: $(c_mean_residual_baseline); V: $(v_mean_residual_baseline))",
    margin=5Plots.mm,
    rightmargin=15Plots.mm,
    xlabel="Date",
    ylabel="Simulated - Observed [ML]",
    linewidth=1.5,
    color=:orange
)

plot!(
    DATES[subperiod], ensemble_residual_subperiod,
    label="Ensemble (C: $(c_mean_residual_ensemble); V: $(v_mean_residual_ensemble))",
    linestyle=:dash,
    linewidth=1.5,
    alpha=0.8,
    color=:green
)

hline!(repeat([0.0], length(DATES[subperiod])), color=:blue, label="", alpha=0.5)

# Display residuals for subperiod (residuals for entire simulation)
msg = text("$(mean_residual_baseline_subperiod)", :left, :orange, 10)
annotate!(DATES[subperiod[end]+10], baseline_residual_subperiod[end], msg, arrows=true)

msg = text("$(mean_residual_ensemble_subperiod)", :left, :green, 10)
annotate!(DATES[subperiod[end]+10], ensemble_residual_subperiod[end], msg, arrows=true)

savefig(res_plot, "$(ensemble_figs)residuals_2008-2009_ensemble.png")

combined_residual_plot = plot(
    comp_plot,
    res_plot,
    size=(1000,400)
)
savefig(combined_residual_plot, "$(ensemble_figs)combined_residuals_2008-2009_ensemble.png")

# Baseline NmKGE vs Ensemble (dry state)
dry_valid = findall(state_active_params["split_mean_NmKGE"][CALIB_LN+1:end] .== 3)
dry_obs = Qo[CALIB_LN+1:end][dry_valid]
il_x = minimum(dry_obs) # identity line
il_y = maximum(dry_obs)
bl_res_plot = Plots.plot(il_x:il_y, label="", color=:red, xlabel="Observed", ylabel="Baseline NKGE'")
StatsPlots.scatter!(dry_obs, base_outflows["NmKGE"][CALIB_LN+1:end][dry_valid], label="", color=:blue)

en_res_plot = Plots.plot(il_x:il_y, label="", color=:red, xlabel="Observed", ylabel="Ensemble")
StatsPlots.scatter!(dry_obs, ensembles[(best_ensemble)][CALIB_LN+1:end][dry_valid], label="", color=:blue)

combined_residual_plot2 = plot(
    bl_res_plot,
    en_res_plot,
    size=(1000,400)
)
savefig(combined_residual_plot2, "$(ensemble_figs)combined_residuals_bl_NmKGE_ensemble.png")

# Cross-comparison of metrics
Qb = base_outflows[model1]
Qb_calib = Qb[calib_ts+1:CALIB_LN]
Qb_valid = Qb[CALIB_LN+1:end]
Qb2 = base_outflows[model2]
Qb2_calib = Qb2[calib_ts+1:CALIB_LN]
Qb2_valid = Qb2[CALIB_LN+1:end]

Qs = state_outflows[model1]
Qs_calib = Qs[calib_ts+1:CALIB_LN]
Qs_valid = Qs[CALIB_LN+1:end]
Qs2 = state_outflows[model2]
Qs2_calib = Qs2[calib_ts+1:CALIB_LN]
Qs2_valid = Qs2[CALIB_LN+1:end]

Qe = ensembles[best_ensemble]
Qe_calib = Qe[calib_ts+1:CALIB_LN]
Qe_valid = Qe[CALIB_LN+1:end]

prefixes = ["Overall", "Wet", "Mid", "Dry"]
periods = ["overall", "calib", "valid"]
col_names = []
for p in periods
    for pre in prefixes
        if lowercase(p) == lowercase(pre)
            push!(col_names, "$(pre)")
            continue
        end
        push!(col_names, "$(pre)_$(p)")
    end
end
pushfirst!(col_names, "Metric")

metric_comparison = DataFrame([Symbol("$(col)")=>repeat([0.0], length(APPROACHES)*5) for col in col_names])
metric_comparison.Metric = vcat(APPROACHES,
                                "State-based Model1 " .* APPROACHES,
                                "State-based Model2 " .* APPROACHES,
                                "Baseline Model1 " .* APPROACHES,
                                "Baseline Model2 " .* APPROACHES)

active_params = state_active_params["NmKGE"]
overall_low, overall_mid, overall_high = separate_state_periods(active_params, calib_ts+1)
calib_low, calib_mid, calib_high = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
valid_low, valid_mid, valid_high = separate_state_periods(active_params, CALIB_LN+1)


function populate_metric_df!(metric, df, row_id, Qobs, Qsim, low, mid, high, period)
    if period == "overall"
        df[row_id, :Overall] = [metric(Qobs, Qsim)]
        df[row_id, :Wet_overall] = [metric(Qobs[low], Qsim[low])]
        df[row_id, :Mid_overall] = [metric(Qobs[mid], Qsim[mid])]
        df[row_id, :Dry_overall] = [metric(Qobs[high], Qsim[high])]
    elseif period == "calib"
        df[row_id, :Overall_calib] = [metric(Qobs, Qsim)]
        df[row_id, :Wet_calib] = [metric(Qobs[low], Qsim[low])]
        df[row_id, :Mid_calib] = [metric(Qobs[mid], Qsim[mid])]
        df[row_id, :Dry_calib] = [metric(Qobs[high], Qsim[high])]
    elseif period == "valid"
        df[row_id, :Overall_valid] = [metric(Qobs, Qsim)]
        df[row_id, :Wet_valid] = [metric(Qobs[low], Qsim[low])]
        df[row_id, :Mid_valid] = [metric(Qobs[mid], Qsim[mid])]
        df[row_id, :Dry_valid] = [metric(Qobs[high], Qsim[high])]
    else
        println("Unknown period given")
    end
end

# Get metrics for all model instances
for (met_name, metric) in zip(APPROACHES, OBJFUNCS)
    # Ensemble
    row_id = (metric_comparison.Metric .== met_name)
    populate_metric_df!(metric, metric_comparison, row_id, Qo[calib_ts+1:end], Qe[calib_ts+1:end], overall_low, overall_mid, overall_high, "overall")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_calib, Qe_calib, calib_low, calib_mid, calib_high, "calib")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_valid, Qe_valid, valid_low, valid_mid, valid_high, "valid")

    row_id = (metric_comparison.Metric .== "State-based Model1 $met_name")
    populate_metric_df!(metric, metric_comparison, row_id, Qo[calib_ts+1:end], Qs[calib_ts+1:end], overall_low, overall_mid, overall_high, "overall")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_calib, Qs_calib, calib_low, calib_mid, calib_high, "calib")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_valid, Qs_valid, valid_low, valid_mid, valid_high, "valid")

    row_id = (metric_comparison.Metric .== "State-based Model2 $met_name")
    populate_metric_df!(metric, metric_comparison, row_id, Qo[calib_ts+1:end], Qs2[calib_ts+1:end], overall_low, overall_mid, overall_high, "overall")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_calib, Qs2_calib, calib_low, calib_mid, calib_high, "calib")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_valid, Qs2_valid, valid_low, valid_mid, valid_high, "valid")

    row_id = (metric_comparison.Metric .== "Baseline Model1 $met_name")
    populate_metric_df!(metric, metric_comparison, row_id, Qo[calib_ts+1:end], Qb[calib_ts+1:end], overall_low, overall_mid, overall_high, "overall")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_calib, Qb_calib, calib_low, calib_mid, calib_high, "calib")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_valid, Qb_valid, valid_low, valid_mid, valid_high, "valid")

    row_id = (metric_comparison.Metric .== "Baseline Model2 $met_name")
    populate_metric_df!(metric, metric_comparison, row_id, Qo[calib_ts+1:end], Qb2[calib_ts+1:end], overall_low, overall_mid, overall_high, "overall")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_calib, Qb2_calib, calib_low, calib_mid, calib_high, "calib")
    populate_metric_df!(metric, metric_comparison, row_id, Qo_valid, Qb2_valid, valid_low, valid_mid, valid_high, "valid")
end

df = metric_comparison

# metric_comparison |> stack |> @df scatter(:variable, :value, color=1:length(APPROACHES))
df_T = transpose_df(df)

CSV.write(joinpath(ensemble_res_path, "all_metrics_ensemble_results.csv"), df)


replacements = [
    "mean_NmKGE"=>"NμKGE'⁻¹", 
    "NmKGE"=>"NKGE'", 
    "NnpKGE"=>"NKGEnp"
]
for w in replacements
    rename!(df_T, replace.(names(df_T), w))
end
rename!(df_T, replace.(names(df_T), "_"=>" "))

spokes = names(df_T)
ensemble_df = df_T[:, (.!(occursin.("Baseline", spokes)) .& .!(occursin.("State-based", spokes)))]
statebased_df = df_T[:, occursin.("State-based Model1", spokes) .| occursin.("Metric", spokes)]
statebased2_df = df_T[:, occursin.("State-based Model2", spokes) .| occursin.("Metric", spokes)]
rename!(statebased_df, replace.(names(statebased_df), "State-based Model1 "=>""))
rename!(statebased2_df, replace.(names(statebased2_df), "State-based Model2 "=>""))
baseline_df = df_T[:, occursin.("Baseline Model1", spokes) .| occursin.("Metric", spokes)]
baseline2_df = df_T[:, occursin.("Baseline Model2", spokes) .| occursin.("Metric", spokes)]
rename!(baseline_df, replace.(names(baseline_df), "Baseline Model1 "=>""))
rename!(baseline2_df, replace.(names(baseline2_df), "Baseline Model2 "=>""))


using PlotlyJS

instance_radar = []  # Performance of instance: overall, calibration period, validation period, dry states
overall_radar = []  # Performance of instance: Overall, wet periods, usual periods, dry periods
instance_names = ["Ensemble", "State-based Model1", "State-based Model2", "Baseline 1", "Baseline 2"]
for (title, df_t) in zip(instance_names,
                         [ensemble_df, statebased_df, statebased2_df, baseline_df, baseline2_df])
    tmp_df = transpose_df(df_t)[!, [:Metric, :Overall, :Overall_calib, :Overall_valid, :Dry_valid]]

    # `r` is the column to use
    # theta is the name of the spokes
    # `name` is the label text
    r = PlotlyJS.plot([
        scatterpolar(tmp_df, r=Symbol(N), theta=:Metric, mode="lines", name=name=replace(N, "_"=>" "))
        for N in names(tmp_df[!, [:Overall, :Overall_calib, :Overall_valid, :Dry_valid]])
    ], Layout(title="$title Performance"))

    push!(instance_radar, r)

    tmp_df = transpose_df(df_t)[!, [:Metric, :Overall, :Wet_overall, :Mid_overall, :Dry_overall]]
    r = PlotlyJS.plot([
            scatterpolar(tmp_df, r=Symbol(N), theta=:Metric, mode="lines", name=replace(N, "_"=>" "), title="$title Performance")
            for N in ["Overall", "Wet_overall", "Mid_overall", "Dry_overall"]
    ], Layout(title="$title Performance"))

    push!(overall_radar, r)
end

[PlotlyJS.savefig(n, joinpath(ensemble_figs, "instance_radar_$i.png")) 
    for (i, n) in zip(instance_names, instance_radar)]

[PlotlyJS.savefig(n, joinpath(ensemble_figs, "overall_radar_$i.png")) 
    for (i, n) in zip(instance_names, overall_radar)]
