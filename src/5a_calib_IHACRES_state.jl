include("_common.jl")


if nworkers() < 2
    addprocs(min(length(APPROACHES), CPU_CORES), exeflags="--project=..")
end


include_everywhere("_common.jl")
include_everywhere("_calib_funcs.jl")

@everywhere calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"
@everywhere calib_fig_path = joinpath(FIG_PATH, "calib_state") * "/"
mkpath(calib_fig_path)
mkpath(calib_param_path)


@everywhere function run_state_calibration(metric)
    approach, objfunc = metric

    calib_start = 1826  # 5-year burn in
    sn, n_id = setup_network("$(DATA_PATH)410730_IHACRES.yml")
    func = nothing
    if !occursin("RMSE", approach)
        func = (obs, sim) -> abs(1.0 - objfunc(obs[calib_start:end], sim[calib_start:end]))
    else
        func = (obs, sim) -> objfunc(obs[calib_start:end], sim[calib_start:end])
    end

    # parameter 3 is `e` which we ignore as we're using PET data.
    target_idx = (1,2,4,5,6,7,8)
    thresholds = [0.4, 0.8, 1.0]
    popsize = 64 * length(target_idx)^2

    # Calibrate network using the BlackBoxOptim package
    # keyword arguments will be passed to the `bboptimize()` function
    # This will set node parameters to the optimal values found
    metric = (active_param_set, param_idxs, obs, sim) -> func(obs, sim)
    result, optobj = state_based_calibrate(sn, n_id, :gw_store, CALIB_CLIMATE, CALIB_FLOW, metric,
                                           target_idx, thresholds;
                                           PopulationSize=popsize, MaxTime=TRAINING_TIME, TraceInterval=300)

    # Save params
    params = best_candidate(result)
    state_str = String(:gw_store)
    outfile = "$(calib_param_path)calib_params_online_$(state_str)_$approach.txt"
    open(outfile, "w") do f
        print(f, join(params, ","))
    end

    Streamfall.save_calibration!(result, optobj, outfile*".opt")

    node = sn[n_id]
    _, x0, __ = param_info(node; with_level=false)
    mod_params = update_partial(x0, target_idx, params)
    active_param_set, param_idxs = online_state_node!(sn, n_id, FULL_CLIMATE, :gw_store, mod_params, thresholds; releases=nothing)

    sim_flow = node.outflow

    calib_data = FULL_DATASET[calib_start:CALIB_LN, "410730_Q"]
    calib_sim = sim_flow[calib_start:CALIB_LN]
    rmse = round(Streamfall.RMSE(calib_data, calib_sim), digits=2)
    nse = round(Streamfall.NSE(calib_data, calib_sim), digits=2)
    kge = round(Streamfall.KGE(calib_data, calib_sim), digits=2)
    plot(FULL_DATASET.Date[calib_start:CALIB_LN], calib_data,
         title="$approach\n(RMSE: $rmse; NSE: $nse; KGE: $kge)",
         legend=:topright,
         label="Historic",
         xlabel="Date",
         ylabel="Streamflow [ML]")
    plot!(CALIB_DATES[calib_start:end], calib_sim, label="Calibration", alpha=0.7)
    savefig(joinpath(calib_fig_path, "$(state_str)_$(approach)_calibration.png"))

    valid_data = FULL_DATASET[CALIB_LN+1:end, "410730_Q"]
    valid_sim = sim_flow[CALIB_LN+1:end]
    rmse = round(Streamfall.RMSE(valid_data, valid_sim), digits=2)
    nse = round(Streamfall.NSE(valid_data, valid_sim), digits=2)
    kge = round(Streamfall.KGE(valid_data, valid_sim), digits=2)
    plot(FULL_DATASET.Date[CALIB_LN+1:end], valid_data,
         title="$approach\n(RMSE: $rmse; NSE: $nse; KGE: $kge)",
         legend=:topleft,
         label="Historic",
         xlabel="Date",
         ylabel="Streamflow [ML]")
    plot!(VALID_DATES, valid_sim, label="Validation", alpha=0.7)
    savefig(joinpath(calib_fig_path, "$(state_str)_$(approach)_validation.png"))
end

# reset!(sn)
# _, x0, __ = param_info(sn[n_id]; with_level=false)
# mod_params = update_partial(x0, target_idx, params)
# active_param_set, param_idxs = online_gw_store_node!(sn, n_id, CALIB_CLIMATE, mod_params, quantiles; releases=nothing)

# obs_flow = TRAINING_FLOW["410730"][calib_ts:end]
# sim_flow = sn[n_id].outflow[calib_ts:CALIB_LN]
# report_metrics(obs_flow, sim_flow)

# recreated_score = metric(active_param_set, param_idxs, TRAINING_FLOW["410730"], sn[n_id].outflow)
# @assert best_fitness(result) == recreated_score

pmap(run_state_calibration, zip(APPROACHES, OBJFUNCS))
