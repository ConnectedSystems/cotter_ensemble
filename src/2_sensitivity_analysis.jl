include("_common.jl")

using GlobalSensitivity


function run_sa(X; node, metric)
    n_params = length(X)

    # drop dummy parameter
    X = X[1:n_params-1]
    update_params!(node, X...)

    Streamfall.run_node!(node, CALIB_CLIMATE)

    score = metric(CALIB_OBS[366:CALIB_LN], node.outflow[366:CALIB_LN])
    Streamfall.reset!(node)

    return score
end


N = 8192
approaches = ["NNSE", "NmKGE", "RMSE", "mean_NmKGE", "NnpKGE", "split_NNSE", "split_NmKGE", "split_RMSE", "split_mean_NmKGE", "split_NnpKGE"]

split_NNSE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NNSE)
split_NmKGE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NmKGE)
split_RMSE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.RMSE)
split_mean_NmKGE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.mean_NmKGE)
split_mean_NnpKGE = (obs, sim) -> Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NnpKGE)

metrics = [
    Streamfall.NNSE, Streamfall.NmKGE, Streamfall.RMSE, Streamfall.mean_NmKGE, Streamfall.NnpKGE,
    split_NNSE, split_NmKGE, split_RMSE, split_mean_NmKGE, split_mean_NnpKGE
]

for (app, metric) in zip(approaches, metrics)
    sn, n_id = setup_network(joinpath(DATA_PATH, "cotter_baseline_calibrated_$(app).yml"))
    node = sn[n_id]

    p_names, ini_vals, p_bounds = param_info(node; with_level=false)

    @info "Running $app..."

    # Set up dummy parameter
    push!(p_names, :dummy)
    push!(ini_vals, 0.0)
    push!(p_bounds, (0.0, 1000.0))

    sa_runner = (X) -> run_sa(X; node=node, metric=metric)
    sobol_results = gsa(sa_runner, Sobol(), p_bounds, N=N)

    headers = "$app Parameter, ST, S1, Threshold"
    res_lines = zip(p_names, sobol_results.ST, sobol_results.S1, repeat([0.4], length(p_names)))
    outfile = "$(DATA_PATH)sa_results_w_dummy.csv"

    open(outfile, "a+") do f
        println(f, headers)
        for line in res_lines
            println(f, join(line, ","))
        end

        print(f, "\n")
    end
end
