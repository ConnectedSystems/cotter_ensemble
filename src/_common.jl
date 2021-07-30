using DataFrames, Query, Dates, CSV, YAML
using Statistics, Distributions
using Plots, StatsPlots
using Streamfall

using Distributed, BlackBoxOptim


Plots.gr(fmt=:png)

CPU_THREADS = Sys.CPU_THREADS
CPU_CORES = Int(CPU_THREADS / 2)


DATA_PATH = "../data/"
FIG_PATH = "../figures/"

TRAINING_TIME = 60*60*24
CALIB_TIME = 60*60*24  # calibration time in seconds
DATE_FORMAT = "YYYY-mm-dd"


#### common data set up
function setup_network(spec)
    network = YAML.load_file(joinpath(DATA_PATH, spec))
    sn = create_network("Gingera Catchment", network)

    # Override base run function as we're using temperature data
    n_id, node = sn["410730"]

    return sn, n_id
end


function include_everywhere(filepath)
    fullpath = joinpath(@__DIR__, filepath)
    @sync for p in procs()
        @async remotecall_wait(include, p, fullpath)
    end
end


# Load data
FULL_DATASET = DataFrame(CSV.File(joinpath(DATA_PATH, "CAMELS-AUS_410730.csv"), comment="#", dateformat=DATE_FORMAT))
const CALIB_LN = 11134

CALIB_OBS = FULL_DATASET[1:CALIB_LN, "410730_Q"]
VALID_OBS = FULL_DATASET[CALIB_LN+1:end, "410730_Q"]

CALIB_FLOW = Dict(
    "410730" => CALIB_OBS
)

VALID_FLOW = Dict(
    "410730" => VALID_OBS
)

CALIB_CLIMATE = Climate(FULL_DATASET[1:CALIB_LN, ["410730_P", "410730_PET"]], "_P", "_PET")
FULL_CLIMATE = Climate(FULL_DATASET, "_P", "_PET")

CALIB_DATES = FULL_DATASET[1:CALIB_LN, "Date"]
VALID_DATES = FULL_DATASET[CALIB_LN+1:end, "Date"]

APPROACHES = ["NNSE", "NmKGE", "RMSE", "mean_NmKGE", "NnpKGE",
              "split_NNSE", "split_NmKGE", "split_RMSE", "split_mean_NmKGE", "split_NnpKGE"]

split_NNSE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NNSE)
split_NmKGE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NmKGE)
split_RMSE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.RMSE)
split_mean_NmKGE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.mean_NmKGE)
split_NnpKGE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NnpKGE)

OBJFUNCS = [
    Streamfall.NNSE, Streamfall.NmKGE, Streamfall.RMSE, Streamfall.mean_NmKGE,
    Streamfall.NnpKGE,
    split_NNSE, split_NmKGE, split_RMSE, split_mean_NmKGE,
    split_NnpKGE
]


"""
For backwards-compatibility.

Earlier approach defined an additional bin for values which were greater than historic extremes.
Because the approach now uses an adaptive "online" approach, this is no longer necessary as the
extreme always fall into the furthest (right-most) bin.
"""
function find_state_vars(cmd, thresholds, params, num_params)
    n_states = length(thresholds) + 1

    return find_state_vars(cmd, thresholds, params, num_params, n_states)
end


function find_state_vars(val, thresholds, params, num_params, n_states)
    # index of parameters (e.g., if 3 partitions, then create 3 arrays of 8 parameters)
    param_idxs = collect(Iterators.partition(1:length(params), num_params))

    lower = thresholds[1]
    thresholds = thresholds[2:end]
    set_id = n_states
    for (idx, step) in enumerate(thresholds)
        if (lower <= val < step)
            set_id = idx
            break
        end

        lower = step
    end

    param_set = param_idxs[set_id]
    node_params = params[param_set]

    return set_id, node_params, param_idxs
end


function report_metrics(obs, sim)
    rmse = Streamfall.RMSE(obs, sim)
    nse = Streamfall.NSE(obs, sim)
    mKGE = Streamfall.mKGE(obs, sim)
    mean_NmKGE = Streamfall.mean_NmKGE(obs, sim)
    NnpKGE = Streamfall.NnpKGE(obs, sim)
    pbias = Streamfall.PBIAS(obs, sim)
    rsr = Streamfall.RSR(obs, sim)
    npKGE = Streamfall.npKGE(obs, sim)
    
    split_NSE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NSE)
    split_mKGE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.mKGE)
    split_RMSE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.RMSE)
    split_NNSE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NNSE)
    split_NmKGE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NmKGE)
    split_mean_NmKGE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.mean_NmKGE)
    split_NnpKGE = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NnpKGE)

    nnse = Streamfall.NNSE(obs, sim)
    nmkge = Streamfall.NmKGE(obs, sim)

    @info "Metrics:" nnse nmkge mean_NmKGE NnpKGE rmse nse mKGE npKGE pbias rsr split_NNSE split_NmKGE split_RMSE split_mean_NmKGE split_NSE split_mKGE split_NnpKGE
end