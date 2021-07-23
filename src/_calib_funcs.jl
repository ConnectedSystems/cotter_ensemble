using Statistics
using Distributions
using OnlineStats
import Streamfall: run_node!


function run_node!(sn, v_id, climate, thresholds, params; releases=nothing)
    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    param_idxs = nothing
    for ts in (1:timesteps)
        # Update param set based on CMD state
        cmd = node.storage[ts]

        cmd_param_set_id, node_params, param_idxs = find_state_vars(cmd, thresholds, params, num_params)
        if (ts == 1) || (cmd_param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = cmd_param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end


"""Same as run_node! except determine thresholds on the fly ("online")"""
function run_online_node!(sn, v_id, climate, params; releases=nothing)
    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    o = Series(FitNormal(), Extrema())
    param_idxs = nothing
    for ts in (1:timesteps)
        cmd = node.storage[ts]

        # Update param set based on CMD state
        # Based on mean +/- std
        fit!(o, cmd)
        (μ, σ), (_, d_max) = value(o)
        thresholds = [μ - σ, μ + σ, d_max]

        cmd_param_set_id, node_params, param_idxs = find_state_vars(cmd, thresholds, params, num_params, length(thresholds))
        if (ts == 1) || (cmd_param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = cmd_param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end


"""Runs node online, determining state thresholds around 
quantiles instead of mean.
"""
function run_online_node!(sn, v_id, climate, params, quantiles; releases=nothing)
    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    o = Quantile(quantiles, b=1000)
    param_idxs = nothing
    for ts in (1:timesteps)
        cmd = node.storage[ts]

        # Update param set based on CMD state
        # Quantiles
        fit!(o, cmd)
        thresholds = value(o)

        cmd_param_set_id, node_params, param_idxs = find_state_vars(cmd, thresholds, params, num_params, length(thresholds))
        if (ts == 1) || (cmd_param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = cmd_param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end


"""Same as run_node! except determine thresholds on the fly ("online")

- offsets : vector of offsets from mean to calculate thresholds. Will always include maximum. 
            Will resolve to [μ + (offset * std)..., max]
"""
function rainfall_online_node!(sn, v_id, climate, params, quantiles; releases=nothing)
    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    o = Quantile(quantiles, b=1000)
    mw = MovingWindow(10, Float64)
    param_idxs = nothing
    for ts in (1:timesteps)
        # Update param set based on state
        P, _ = climate_values(node, climate, ts)
        tmp = P == 0.0 ? 0.0 : log(P)
        fit!(mw, tmp)
        w_P = value(fit!(Mean(weight = ExponentialWeight(.01)), value(mw)))
        fit!(o, w_P)
        thresholds = value(o)

        param_set_id, node_params, param_idxs = find_state_vars(w_P, thresholds, params, num_params, length(thresholds))
        if (ts == 1) || (param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end



"""Same as run_node! except determine thresholds on the fly ("online")

- offsets : vector of offsets from mean to calculate thresholds. Will always include maximum. 
            Will resolve to [μ + (offset * std)..., max]
"""
function online_gw_store_node!(sn, v_id, climate, params, quantiles; releases=nothing)
    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    o = Quantile(quantiles, b=1000)
    # o = Series(FitNormal(), Extrema())
    # mw = MovingWindow(365, Float64)
    param_idxs = nothing
    for ts in (1:timesteps)
        # Update param set based on state
        # tmp = node.gw_store[ts]
        # fit!(mw, tmp)
        # gw_store = value(fit!(Mean(weight = ExponentialWeight(.01)), value(mw)))
        # fit!(o, gw_store)
        # thresholds = value(o)
        tmp = node.storage[ts]
        gw_store = tmp == 0.0 ? 0.0 : log(tmp)
        
        fit!(o, gw_store)
        thresholds = value(o)

        # (d_μ, d_σ), (_, d_max) = value(o)
        # thresholds = [d_μ .+ (d_σ .* offsets)..., d_max]
        
        param_set_id, node_params, param_idxs = find_state_vars(gw_store, thresholds, params, num_params, length(thresholds))
        if (ts == 1) || (param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end



"""Custom objective function"""
function obj_func(params, climate, sn, v_id, calib_data, thresholds, metric; extraction=nothing)
    this_node = sn[v_id]

    param_idxs = nothing
    active_param_set = nothing

    try
        active_param_set, param_idxs = run_node!(sn, v_id, climate, thresholds, params; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node_data = this_node.outflow
    h_data = calib_data[this_node.name]

    # Calculate score
    score = metric(active_param_set, param_idxs, h_data, node_data)

    # reset to clear stored values
    Streamfall.reset!(this_node)

    return score
end


"""Online objective function"""
function online_obj_func(params, climate, sn, v_id, calib_data, metric; extraction=nothing)
    this_node = sn[v_id]

    param_idxs = nothing
    active_param_set = nothing

    try
        active_param_set, param_idxs = run_online_node!(sn, v_id, climate, params; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node_data = this_node.outflow
    h_data = calib_data[this_node.name]

    # Calculate score
    score = metric(active_param_set, param_idxs, h_data, node_data)

    # reset to clear stored values
    Streamfall.reset!(this_node)

    return score
end


"""
Update base parameter values with a partial list of new values by index.

# Arguments
-----------
- base_params : the base parameter set to update
- target_idx : the index position of parameters to update
- updated_params : a flattened array of new parameter values in same order as given indices
"""
function update_partial(base_params, target_idx, updated_params)
    # Update base parameter values with given updated values
    n_selected_params = length(target_idx)
    tmp_idxs = collect(Iterators.partition(1:length(updated_params), n_selected_params))
    mod_params = []
    for idxs in tmp_idxs
        bin_params = updated_params[idxs]
        new_params = copy(base_params)
        for (idx, new_value) in zip(target_idx, bin_params)
            new_params[idx] = new_value
        end
        append!(mod_params, new_params)
    end

    return mod_params
end


"""Partial online objective function based on specified quantiles

Targets specific parameters identified through sensitivity analysis.
See script `3b_sensitivity_analysis`

Parameters > 0.05 total order sensitivity are selected.
"""
function partial_online_obj_func(params, base_params, climate, sn, v_id, calib_data, metric, target_idx, quantiles;
                                 extraction=nothing)
    node = sn[v_id]

    param_idxs = nothing
    active_param_set = nothing

    # Update base parameter values with given updated values
    mod_params = update_partial(base_params, target_idx, params)
    try
        active_param_set, param_idxs = run_online_node!(sn, v_id, climate, mod_params, quantiles; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node_data = node.outflow
    h_data = calib_data[node.name]

    # Calculate score
    score = metric(active_param_set, param_idxs, h_data, node_data)

    # reset to clear stored values
    Streamfall.reset!(node)

    return score
end


"""Partial online objective function

Targets specific parameters identified through sensitivity analysis.

Relevant parameters are: d2, f, b, alpha

See script `3b_sensitivity_analysis`

Parameters > 0.1 total order sensitivity are selected.
"""
function partial_online_obj_func(params, base_params, climate, sn, v_id, calib_data, metric, target_idx;
                                 extraction=nothing)
    node = sn[v_id]

    param_idxs = nothing
    active_param_set = nothing

    # Update base parameter values with given updated values
    mod_params = update_partial(base_params, target_idx, params)
    try
        active_param_set, param_idxs = run_online_node!(sn, v_id, climate, mod_params; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node_data = node.outflow
    h_data = calib_data[node.name]

    # Calculate score
    score = metric(active_param_set, param_idxs, h_data, node_data)

    # reset to clear stored values
    Streamfall.reset!(node)

    return score
end


function rainfall_obj_func(x, base_params, climate, sn, v_id, calib_data, metric, target_idx, offsets)
    param_idxs = nothing
    active_param_set = nothing
    extraction = nothing

    # Update base parameter values with given updated values
    mod_params = update_partial(base_params, target_idx, x)

    try
        active_param_set, param_idxs = rainfall_online_node!(sn, v_id, climate, mod_params, offsets; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node = sn[v_id]
    node_data = node.outflow
    h_data = calib_data[node.name]

    # Calculate score
    score = metric(active_param_set, param_idxs, h_data, node_data)

    # reset to clear stored values
    Streamfall.reset!(node)

    return score


end


function gw_store_obj_func(x, base_params, climate, sn, v_id, calib_data, metric, target_idx, offsets)
    param_idxs = nothing
    active_param_set = nothing
    extraction = nothing

    # Update base parameter values with given updated values
    mod_params = update_partial(base_params, target_idx, x)

    try
        active_param_set, param_idxs = online_gw_store_node!(sn, v_id, climate, mod_params, offsets; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node = sn[v_id]
    obs = calib_data[node.name]
    sim = node.outflow

    # Calculate score
    score = metric(active_param_set, param_idxs, obs, sim)

    # reset to clear stored values
    Streamfall.reset!(node)

    return score
end


"""Calibration with online statistics with thresholds based around specified standard deviation offsets.

Generalized to handle any one state variable.
"""
function rainfall_state_based_calibrate(sn, v_id, climate, calib_data, metric, target_idx, offsets; kwargs...)

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
    opt_func = x -> rainfall_obj_func(x, x0, climate, sn, v_id, calib_data, metric, target_idx, offsets)

    # Set up parameters for each CMD state
    n_states = length(offsets)  # lower, center, upper
    param_bounds = repeat(param_bounds, n_states)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  kwargs...)

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(node.name)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    return res, opt
end


"""Generic run_node method that handles any single state variable.

Log-transforms the state value.
"""
function online_state_node!(sn, v_id, climate, state, params, quantiles; releases=nothing)
    node = sn[v_id]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    num_params = length(x0)

    timesteps = sim_length(climate)
    active_param_set = zeros(Int, timesteps)

    o = Quantile(quantiles, b=1000)
    state_var = getfield(node, state)
    param_idxs = nothing
    for ts in (1:timesteps)
        # Update param set based on state
        tmp = state_var[ts]
        if tmp == 0.0
            state_val = 0.0
        else
            state_val = log(tmp)
        end

        # Update param set based on state
        fit!(o, state_val)
        thresholds = value(o)
        param_set_id, node_params, param_idxs = find_state_vars(state_val, thresholds, params, num_params, length(thresholds))
        if (ts == 1) || (param_set_id != active_param_set[ts-1])
            update_params!(node, node_params...)
        end

        # record timestep in which this param was active
        active_param_set[ts] = param_set_id

        Streamfall.run_node!(sn, v_id, climate, ts; extraction=releases)
    end

    return active_param_set, param_idxs
end


function state_obj_func(x, base_params, climate, sn, v_id, state, calib_data, metric, target_idx, offsets)
    param_idxs = nothing
    active_param_set = nothing
    extraction = nothing

    # Update base parameter values with given updated values
    mod_params = update_partial(base_params, target_idx, x)

    try
        active_param_set, param_idxs = online_state_node!(sn, v_id, climate, state, mod_params, offsets; releases=extraction)
    catch err
        if err isa AssertionError
            return 9999.0
        end

        throw(err)
    end

    node = sn[v_id]
    node_data = node.outflow
    h_data = calib_data[node.name]

    # Calculate score
    score = metric(active_param_set, param_idxs, h_data, node_data)

    # reset to clear stored values
    Streamfall.reset!(node)

    return score
end


"""Calibration with online statistics with thresholds based around specified standard deviation offsets.

Generalized to handle any one state variable.
"""
function state_based_calibrate(sn, v_id, state_name::Symbol, climate, calib_data, metric, 
                                       target_idx, thresholds; kwargs...)

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
    opt_func = x -> state_obj_func(x, x0, climate, sn, v_id, state_name, 
                                   calib_data, metric, target_idx, thresholds)

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


"""Calibration with online statistics with thresholds based around specified standard deviation offsets"""
function gw_store_based_calibrate(sn, v_id, climate, calib_data, metric, target_idx, quantiles; kwargs...)

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
    opt_func = x -> gw_store_obj_func(x, x0, climate, sn, v_id, calib_data, metric, target_idx, quantiles)

    # Set up parameters for each CMD state
    n_states = length(quantiles)  # lower, center, upper
    param_bounds = repeat(param_bounds, n_states)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  kwargs...)

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(node.name)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    return res, opt
end


"""
Meta-metric that applies the given metric (`metric`) to each threshold and returns
the combination of these with the metric applied to the entire time series 
(using the provided `agg_func`).
"""
function bin_metric(active_param_set, param_idxs, h_data, n_data, metric, agg_func)
    n_len = length(n_data)
    @assert length(h_data) == n_len
    sub_scores = []
    for param_set in 1:length(param_idxs)
        out_set = findall(active_param_set .== param_set)
        obs = h_data[out_set]
        sim = n_data[out_set]
        n_sim = length(sim)
        push!(sub_scores, metric(obs, sim))  #  * (n_sim / n_len)
    end

    # push!(sub_scores, metric(h_data, n_data))

    return agg_func(sub_scores)
end
