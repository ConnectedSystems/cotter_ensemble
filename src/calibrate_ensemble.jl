calibrate!(sn, CALIB_CLIMATE, CALIB_FLOW; metric=func,
               TraceInterval=300, MaxTime=60*60*6, PopulationSize=1000)



include("_common.jl")
include("_calib_funcs.jl")


# model1 = ...
# model2 = ...

# Two models, one for the range indicated, another for the conditions above
thresholds = [0.0, 0.8]

function state_ensemble_metric(active_param_set, target_state, obs, sim, metric)
    subset = findall(active_param_set .== target_state)

    obs_sub = obs[subset]
    sim_sub = sim[subset]

    if length(sim_sub) == 0
        return 9999.9
    end

    return metric(obs_sub, sim_sub)
end

function test_calibrate(models, thresholds, climate)
end


function state_ensemble_calibrate(models, state_name::Symbol, climate, calib_data, metric, target_idx, thresholds; 
                                  burn_in=1826, log_transform=false, kwargs...)
    # Set defaults as necessary
    defaults = (;
    MaxTime=CALIB_TIME,
    TraceInterval=300.0,
    PopulationSize=125
    )
    kwargs = merge(defaults, kwargs)

    @assert length(thresholds[1:end]) == length(models) || "Given thresholds do not match number of models"

    # setup targets
    for m in models



    # ...
end



# Calibrate each model
# pmap(calibrate!, models)

# Combine, weighting based on catchment state


"""Generic run_node method that handles any single state variable.

Develops an "online" running statistic of states from the first N days, 
which remain static for the rest of the analysis period.
"""
function online_burn_state_node!(sn, v_id, climate, state, params, quantiles; burn_in=1826, releases=nothing, log_transform=false)

    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    if state != :rainfall
        state_var = getfield(node, state)
    else
        tmp_climate = Streamfall.subcatchment_data(node, climate)
        # hard coded, I know :(
        state_var = tmp_climate[:, "410730_P"]
        # TODO: Change over when ready
        # state_var = Streamfall.rainfall_data(node, climate)
    end

    o = Quantile(quantiles, b=1000)
    param_idxs = nothing
    thresholds = nothing
    for ts in (1:timesteps)
        # Update param set based on state
        state_val = state_var[ts]
        if log_transform
            # add constant to avoid log(0)
            tmp = log(1.0 + state_val)
        else
            tmp = state_val
        end

        if ts < burn_in
            # Update thresholds based on state value so far
            fit!(o, tmp)
            thresholds = value(o)
        end

        param_set_id, node_params, param_idxs = find_state_vars(tmp, thresholds, params, num_params, length(thresholds))
        if (ts == 1) || (param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end


# Update parameters, run the model
# in obj function, only consider the conditions indicated by the thresholds
# in running, switch based on state






"""Calibration with online statistics with thresholds based around specified standard deviation offsets.

Generalized to handle any one state variable.
"""
function ensemble_burn_state_based_calibrate(sn, v_id, state_name::Symbol, climate, calib_data, metric, 
                                             target_idx, thresholds; burn_in=1826, log_transform=false, kwargs...)

    # Set defaults as necessary
    defaults = (;
        MaxTime=CALIB_TIME,
        TraceInterval=300.0,
        PopulationSize=125
    )
    kwargs = merge(defaults, kwargs)

    node = sn[v_id]

    # Get node parameters
    _, x0, param_bounds = param_info(node; with_level=false)
    param_bounds = [param_bounds[i] for i in target_idx]  # subset bounds array

    # Create new optimization function
    opt_func = x -> burn_state_obj_func(x, x0, climate, sn, v_id, state_name, 
                                        calib_data, metric, target_idx, thresholds; burn_in=burn_in, log_transform=log_transform)

    # Set up parameters for each state threshold
    n_states = length(thresholds)  # lower, center, upper
    param_bounds = repeat(param_bounds, n_states)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  kwargs...)

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(node.name)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    return res, opt
end









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
    target_idx = (5,6)
    thresholds = [0.0, 0.3, 0.7, 0.9]
    popsize = 64 * length(target_idx)^2
    state = :gw_store
    state_str = String(state)

    # Calibrate network using the BlackBoxOptim package
    # keyword arguments will be passed to the `bboptimize()` function
    # This will set node parameters to the optimal values found
    metric = (active_param_set, param_idxs, obs, sim) -> bin_metric(active_param_set, param_idxs, obs, sim, func, nothing)
    result, optobj = burn_state_based_calibrate(sn, n_id, state, CALIB_CLIMATE, CALIB_FLOW, metric,
                                                target_idx, thresholds;
                                                log_transform=true,
                                                Method=:borg_moea,
                                                FitnessScheme=ParetoFitnessScheme{length(thresholds)}(is_minimizing=true),
                                                PopulationSize=popsize, MaxTime=TRAINING_TIME, TraceInterval=300)

    # Save params
    params = best_candidate(result)
    outfile = "$(calib_param_path)calib_params_online_burn-in_multi_log_alt2_$(state_str)_$(approach).txt"
    open(outfile, "w") do f
        print(f, join(params, ","))
    end

    Streamfall.save_calibration!(result, optobj, outfile*".opt")

    node = sn[n_id]
    _, x0, __ = param_info(node; with_level=false)
    mod_params = update_partial(x0, target_idx, params)
    active_param_set, param_idxs = online_burn_state_node!(sn, n_id, FULL_CLIMATE, state, mod_params, thresholds; log_transform=true, releases=nothing)

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
    savefig(joinpath(calib_fig_path, "burn-in_multi_log_alt2_$(state_str)_$(approach)_calibration.png"))

    qqplot(calib_data, calib_sim, title="Multi-Obj Log Calibration\n($(state_str); $approach)")
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(joinpath(calib_fig_path, "burn-in_multi_log_alt2_$(state_str)_$(approach)_qq_calib.png"))

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
    savefig(joinpath(calib_fig_path, "burn-in_multi_log_$(state_str)_$(approach)_validation.png"))

    qqplot(valid_data, valid_sim, title="Multi-Obj Log Validation\n($(state_str); $approach)")
    xlabel!("Historic")
    ylabel!("Simulated")
    savefig(joinpath(calib_fig_path, "burn-in_multi_log_alt2_$(state_str)_$(approach)_qq_valid.png"))
end


pmap(run_state_calibration, zip(APPROACHES, OBJFUNCS))
