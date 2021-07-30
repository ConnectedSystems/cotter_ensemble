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

    calib_start = 1826  # 5-year burn in
    sn, n_id = setup_network("$(baseline_calib)cotter_baseline_IHACRES_$approach.yml")
    func = nothing
    if !occursin("RMSE", approach)
        func = (obs, sim) -> abs(1.0 - objfunc(obs[calib_start:end], sim[calib_start:end]))
    else
        func = (obs, sim) -> objfunc(obs[calib_start:end], sim[calib_start:end])
    end

    # Use baseline calibration parameters as the base
    # Varying just 2 parameters
    target_idx = (5,6,7,8)
    thresholds = [0.0, 0.99]
    popsize = 64 * length(target_idx)^2
    state = :gw_store
    state_str = String(state)

    # Calibrate network using the BlackBoxOptim package
    # keyword arguments will be passed to the `bboptimize()` function
    # This will set node parameters to the optimal values found
    weightedfitness(f) = f[1] + f[2]*5.0
    metric = (active_param_set, param_idxs, obs, sim) -> bin_metric(active_param_set[calib_start:CALIB_LN], param_idxs, obs, sim, func, nothing)
    result, optobj = burn_state_based_calibrate(sn, n_id, state, CALIB_CLIMATE, CALIB_FLOW, metric,
                                                target_idx, thresholds;
                                                log_transform=false,
                                                Method=:borg_moea,
                                                FitnessScheme=ParetoFitnessScheme{length(thresholds)}(is_minimizing=true, aggregator=weightedfitness),
                                                PopulationSize=popsize, MaxTime=TRAINING_TIME, TraceInterval=300)

    # Save params
    params = best_candidate(result)
    outfile = "$(calib_param_path)weighted_top_end_$(state_str)_$(approach).txt"
    open(outfile, "w") do f
        print(f, join(params, ","))
    end

    Streamfall.save_calibration!(result, optobj, outfile*".opt")

    best_score = best_fitness(result)

    reset!(sn)
    node = sn[n_id]
    _, x0, __ = param_info(node; with_level=false)
    mod_params = update_partial(x0, target_idx, params)
    active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, thresholds; log_transform=false, releases=nothing)

    sim_flow = node.outflow
    @assert best_score == metric(active_param_set[1:CALIB_LN], param_idxs, CALIB_OBS[1:CALIB_LN], sim_flow[1:CALIB_LN])

    calib_sim = sim_flow[calib_start:CALIB_LN]

    calib_data = FULL_DATASET[calib_start:CALIB_LN, "410730_Q"]
    rmse = round(Streamfall.RMSE(calib_data, calib_sim), digits=2)
    nse = round(Streamfall.NSE(calib_data, calib_sim), digits=2)
    kge = round(Streamfall.KGE(calib_data, calib_sim), digits=2)
    calib_flow_plot = plot(FULL_DATASET.Date[calib_start:CALIB_LN], calib_data,
         title="$approach\n(RMSE: $rmse; NSE: $nse; KGE: $kge)",
         titlefontsize=8,
         legend=:topright,
         label="Historic",
         xlabel="Date",
         ylabel="Streamflow [ML]")
    plot!(CALIB_DATES[calib_start:end], calib_sim, label="Calibration", alpha=0.7)
    # savefig(joinpath(calib_fig_path, "weighted_top_end_$(state_str)_$(approach)_calibration.png"))

    calib_flow_qq = qqplot(calib_data, calib_sim, title="Restricted Weighted Top End Calibration\n($(state_str); $approach)", titlefontsize=8)
    xlabel!("Historic")
    ylabel!("Simulated")

    savefig(plot(calib_flow_plot, calib_flow_qq), joinpath(calib_fig_path, "weighted_top_end_$(state_str)_$(approach)_flow_qq_calib.png"))

    valid_data = FULL_DATASET[CALIB_LN+1:end, "410730_Q"]
    valid_sim = sim_flow[CALIB_LN+1:end]
    rmse = round(Streamfall.RMSE(valid_data, valid_sim), digits=2)
    nse = round(Streamfall.NSE(valid_data, valid_sim), digits=2)
    kge = round(Streamfall.KGE(valid_data, valid_sim), digits=2)
    valid_flow_plot = plot(FULL_DATASET.Date[CALIB_LN+1:end], valid_data,
         title="$approach\n(RMSE: $rmse; NSE: $nse; KGE: $kge)",
         titlefontsize=8,
         legend=:topleft,
         label="Historic",
         xlabel="Date",
         ylabel="Streamflow [ML]")
    plot!(VALID_DATES, valid_sim, label="Validation", alpha=0.7)
    # savefig(joinpath(calib_fig_path, "weighted_top_end_$(state_str)_$(approach)_validation.png"))
    # savefig(plot(flow_plot, flow_qq), joinpath(calib_fig_path, "weighted_top_end_$(state_str)_$(approach)_flow_qq_calib.png"))
    # savefig(joinpath(calib_fig_path, "weighted_top_end_$(state_str)_$(approach)_validation.png"))

    valid_flow_qq = qqplot(valid_data, valid_sim, title="Weighted Top End Validation\n($(state_str); $approach)", titlefontsize=8,)
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(plot(valid_flow_plot, valid_flow_qq),
            joinpath(calib_fig_path, "weighted_top_end_$(state_str)_$(approach)_flow_qq_valid.png"))
end


pmap(run_state_calibration, zip(APPROACHES, OBJFUNCS))
