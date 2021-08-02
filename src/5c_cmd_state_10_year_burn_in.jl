include("_common.jl")
include("_calib_funcs.jl")

if nworkers() < 2
    addprocs(min(length(APPROACHES), CPU_CORES), exeflags="--project=..")
end


include_everywhere("_common.jl")
include_everywhere("_calib_funcs.jl")

@everywhere baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
@everywhere calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"
@everywhere calib_fig_path = joinpath(FIG_PATH, "calib_state") * "/"
mkpath(calib_fig_path)
mkpath(calib_param_path)


@everywhere function run_state_calibration(metric)
    approach, objfunc = metric

    wid = Distributed.myid()
    @info "Worker $wid : Calibrating $approach"

    calib_start = (1825*2)+1  # 10-year burn in
    sn, n_id = setup_network("$(baseline_calib)cotter_baseline_IHACRES_$approach.yml")
    func = nothing
    if !occursin("RMSE", approach)
        func = (obs, sim) -> abs(1.0 - objfunc(obs, sim))
    else
        func = (obs, sim) -> objfunc(obs, sim)
    end

    # Use baseline calibration parameters as the base
    target_idx = (1,2,5,6,7,8)
    thresholds = [0.0, 0.2, 0.90]
    popsize = 64 * length(target_idx)^2
    state = :storage
    state_str = String(state)

    # Calibrate network using the BlackBoxOptim package
    # keyword arguments will be passed to the `bboptimize()` function
    # This will set node parameters to the optimal values found
    metric = (active_param_set, param_idxs, obs, sim) -> bin_metric(active_param_set[calib_start:CALIB_LN], param_idxs, obs[calib_start:CALIB_LN], sim[calib_start:CALIB_LN], func, nothing)
    result, optobj = burn_state_based_calibrate(sn, n_id, state, CALIB_CLIMATE, CALIB_FLOW, metric,
                                                target_idx, thresholds;
                                                burn_in=calib_start,
                                                log_transform=true,
                                                Method=:borg_moea,
                                                FitnessScheme=ParetoFitnessScheme{length(thresholds)}(is_minimizing=true),
                                                PopulationSize=popsize, MaxTime=TRAINING_TIME, TraceInterval=300)

    # Save params
    params = best_candidate(result)
    outfile = "$(calib_param_path)restricted_log_$(state_str)_$(approach)_10_year.txt"
    open(outfile, "w") do f
        print(f, join(params, ","))
    end

    Streamfall.save_calibration!(result, optobj, outfile*".opt")

    best_score = best_fitness(result)

    reset!(sn)
    node = sn[n_id]
    _, x0, __ = param_info(node; with_level=false)
    mod_params = update_partial(x0, target_idx, params)
    active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, thresholds; burn_in=calib_start, log_transform=true, releases=nothing)

    sim_flow = node.outflow
    @assert best_score == metric(active_param_set[1:CALIB_LN], param_idxs, CALIB_OBS[1:CALIB_LN], sim_flow[1:CALIB_LN])

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
    savefig(joinpath(calib_fig_path, "restricted_log_$(state_str)_$(approach)_10_year_calibration.png"))

    qqplot(calib_data, calib_sim, title="4-parameter log Calibration\n($(state_str); $approach)")
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(joinpath(calib_fig_path, "restricted_log_$(state_str)_$(approach)_10_year_qq_calib.png"))

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
    savefig(joinpath(calib_fig_path, "restricted_log_$(state_str)_$(approach)_10_year_validation.png"))

    qqplot(valid_data, valid_sim, title="4-parameter log Validation\n($(state_str); $approach)")
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(joinpath(calib_fig_path, "restricted_log_$(state_str)_$(approach)_10_year_qq_valid.png"))
end


pmap(run_state_calibration, zip(APPROACHES, OBJFUNCS))
