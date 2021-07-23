include("_common.jl")


if nworkers() < 2
    addprocs(min(length(APPROACHES), CPU_CORES), exeflags="--project=..")
end


include_everywhere("_common.jl")
@everywhere baseline_path = joinpath(DATA_PATH, "baselines") * "/"
mkpath(baseline_path)

@everywhere calib_fig_path = joinpath(FIG_PATH, "calib_baseline") * "/"
mkpath(calib_fig_path)


@everywhere function run_calibration(metric)
    approach, objfunc = metric

    sn, n_id = setup_network("410730_IHACRES.yml")
    node = sn[n_id]

    # Calibrate network using the BlackBoxOptim package
    # keyword arguments will be passed to the `bboptimize()` function
    # This will set node parameters to the optimal values found
    calib_start = 366
    func = nothing
    if !occursin("RMSE", approach)
        func = (obs, sim) -> abs(1.0 - objfunc(obs[calib_start:end], sim[calib_start:end]))
    else
        func = (obs, sim) -> objfunc(obs[calib_start:end], sim[calib_start:end])
    end

    calibrate!(sn, CALIB_CLIMATE, CALIB_FLOW;
               metric=func, TraceInterval=60*60*4, MaxTime=10, PopulationSize=1000)

    # Save calibrated network spec to file
    Streamfall.save_network_spec(sn, "$(baseline_path)cotter_baseline_IHACRES_$approach.yml")

    reset!(node)

    # Run calibrated node
    Streamfall.run_node!(sn, n_id, FULL_CLIMATE)
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
    savefig(joinpath(calib_fig_path, "$(approach)_calibration.png"))

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
    savefig(joinpath(calib_fig_path, "$(approach)_validation.png"))
end

pmap(run_calibration, zip(APPROACHES, OBJFUNCS))
