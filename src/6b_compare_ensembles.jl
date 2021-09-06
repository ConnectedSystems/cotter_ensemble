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


calib_ts = 3651
state = :storage
state_str = String(state)

Qo = FULL_DATASET[:, "410730_Q"]
Qo_calib = CALIB_OBS[calib_ts+1:CALIB_LN]
Qo_valid = VALID_OBS


# Load metric results
baseline_calib = joinpath(DATA_PATH, "baselines") * "/"
calib_param_path = joinpath(DATA_PATH, "calib_params") * "/"

if !@isdefined best_NmKGE_ensemble
    baseline_NmKGE_results = CSV.File(joinpath(ensemble_res_path, "overall_baseline_NmKGE_results.csv")) |> DataFrame
    state_NmKGE_results = CSV.File(joinpath(ensemble_res_path, "overall_state_NmKGE_results.csv")) |> DataFrame
    base_base_NmKGE_results = CSV.File(joinpath(ensemble_res_path, "overall_base-base_NmKGE_results.csv")) |> DataFrame
    base_state_NmKGE_results = CSV.File(joinpath(ensemble_res_path, "overall_base-state_NmKGE_results.csv")) |> DataFrame
    state_state_NmKGE_results = CSV.File(joinpath(ensemble_res_path, "overall_state-state_NmKGE_results.csv")) |> DataFrame

    baseline_A2k_results = CSV.File(joinpath(ensemble_res_path, "overall_baseline_A2k_results.csv")) |> DataFrame
    state_A2k_results = CSV.File(joinpath(ensemble_res_path, "overall_state_A2k_results.csv")) |> DataFrame
    base_base_A2k_results = CSV.File(joinpath(ensemble_res_path, "overall_base-base_A2k_results.csv")) |> DataFrame
    base_state_A2k_results = CSV.File(joinpath(ensemble_res_path, "overall_base-state_A2k_results.csv")) |> DataFrame
    state_state_A2k_results = CSV.File(joinpath(ensemble_res_path, "overall_state-state_A2k_results.csv")) |> DataFrame

    active_params = CSV.File(joinpath(ensemble_res_path, "base_NmKGE_state.csv")) |> DataFrame
    active_params = active_params[:, "NmKGE"]
end


function split_calib_valid_cols(df)
    calib_df = df[:, [:Overall_calib, :Wet_calib, :Mid_calib, :Dry_calib]]
    valid_df = df[:, [:Overall_valid, :Wet_valid, :Mid_valid, :Dry_valid]]
    rename!(calib_df, Dict(
        "Overall_calib"=>"Overall",
        "Wet_calib"=>"Wet",
        "Mid_calib"=>"Usual",
        "Dry_calib"=>"Dry",
    ))

    rename!(valid_df, Dict(
        "Overall_valid"=>"Overall",
        "Wet_valid"=>"Wet",
        "Mid_valid"=>"Usual",
        "Dry_valid"=>"Dry",
    ))

    return calib_df, valid_df
end


function create_violin(calib, valid, title)
    vplot = violin(title=title, ylabel="NKGE'", label=["Calibration", "Validation"], legend=:bottomleft)
    for (idx, n) in enumerate(names(calib))
        lab = idx > 1 ? "" : "Calibration"
        violin!([n], calib[:, n], side=:left, color="lightblue", label=lab)
        dotplot!([n], calib[:, n], side=:left, marker=(:black,stroke(0)), markersize=3, label="")

        # don't plot if no state in validation period
        if !all(ismissing.(valid[:, n]))
            lab = idx > 1 ? "" : "Validation"
            violin!([n], replace(valid[:, n], missing=>0.0), side=:right, color="orange", label=lab)
            dotplot!([n], replace(valid[:, n], missing=>0.0), side=:right, marker=(:black,stroke(0)), markersize=3, label="")
        end
    end

    return vplot
end


function relabel_approaches(dfs; cols=false)
    replacements = [
        "mean_NmKGE"=>"NμKGE'⁻¹",
        "NmKGE"=>"NKGE'",
        "NnpKGE"=>"NKGEnp"
    ]
    
    for df in dfs
        if !cols
            for idx in 1:nrow(df)
                app_name = df[idx, :Approach]
                for w in replacements
                    app_name = replace(app_name, w)
                end
                app_name = replace(app_name, "_"=>" ")
                df[idx, :Approach] = app_name
            end
        else
            for w in replacements
                rename!(df, replace.(names(df), w))
            end
            rename!(df, replace.(names(df), "_"=>" "))
        end
    end
end


df_names = ["Baseline", "State-based", "Base-Base", "Base-State", "State-State"]
NmKGE_dfs = [baseline_NmKGE_results, state_NmKGE_results, base_base_NmKGE_results, base_state_NmKGE_results, state_state_NmKGE_results]
A2k_dfs = [baseline_A2k_results, state_A2k_results, base_base_A2k_results, base_state_A2k_results, state_state_A2k_results]

relabel_approaches(A2k_dfs)
relabel_approaches(NmKGE_dfs)


# Find top performing ensemble
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


function sort_scores(df; col=:Overall_valid, rev=false)
    df = deepcopy(df)

    # Do not consider any model that indicates
    # that drier conditions are not experienced
    filter!(row -> !(ismissing(row.Dry_valid)),  df)

    df[:, :Mean_valid] = mean.(skipmissing.(eachrow(df[:, [:Wet_valid, :Mid_valid, :Dry_valid]])))

    # :Mean_valid, 
    sort!(df, [:Mean_valid, col], rev=rev)
end


function create_readable_name(text)
    if occursin(",", text)
        text = join(eval(Meta.parse(text)), " - ")
    end

    return text
end


function create_NmKGE_label(label, Qo, Qm, calib_pos=nothing, valid_pos=nothing)
    local Qo_calib = Qo[calib_ts+1:CALIB_LN]
    local Qo_valid = Qo[CALIB_LN+1:end]
    local Qm_calib = Qm[calib_ts+1:CALIB_LN]
    local Qm_valid = Qm[CALIB_LN+1:end]

    if !isnothing(calib_pos)
        calib_score = round(Streamfall.NmKGE(Qo_calib[calib_pos], Qm_calib[calib_pos]), digits=2)
        valid_score = round(Streamfall.NmKGE(Qo_valid[valid_pos], Qm_valid[valid_pos]), digits=2)
    else
        calib_score = round(Streamfall.NmKGE(Qo_calib, Qm_calib), digits=2)
        valid_score = round(Streamfall.NmKGE(Qo_valid, Qm_valid), digits=2)
    end

    label = "$(label) (C: $(calib_score); V: $(valid_score))"
    return label
end


sorted_NmKGEs = map((df) -> sort_scores(df, col=:Overall_valid, rev=true), NmKGE_dfs)
sorted_A2ks = map((df) -> sort_scores(df, col=:Overall_valid), A2k_dfs)

best_A2k_ensemble = nothing
best_A2k_score = 9999.9
best_A2k_id = nothing
for (df_id, s_A2k) in enumerate(sorted_A2ks)
    score = mean(skipmissing(s_A2k[1, :Mean_valid]))

    println(df_names[df_id], " - Best A2k score: ", score)

    global best_A2k_score
    global best_A2k_id
    global best_A2k_ensemble
    if score < best_A2k_score
        best_A2k_score = score
        best_A2k_ensemble = s_A2k[1, :Approach]
        best_A2k_id = df_id
    end
end

best_performing_A2k_ensemble = sorted_A2ks[best_A2k_id][1, :]
best_performing_A2k_approach = df_names[best_A2k_id]
println("Best Performing: ", best_performing_A2k_approach, " ", best_performing_A2k_ensemble[:Approach], " Score: ", best_A2k_score)


best_A2k_col = create_readable_name(best_A2k_ensemble)

best_NmKGE_ensemble = nothing
best_NmKGE_score = -1.0
best_NmKGE_id = nothing
for (df_id, s_NmKGE) in enumerate(sorted_NmKGEs)
    score = mean(skipmissing(s_NmKGE[1, :Mean_valid]))

    println(df_names[df_id], " - Best NmKGE: ", score)

    global best_NmKGE_score
    global best_NmKGE_id
    if score > best_NmKGE_score
        best_NmKGE_score = score
        best_NmKGE_ensemble = s_NmKGE[1, :Approach]
        best_NmKGE_id = df_id
    end
end

best_performing_baseline = sorted_NmKGEs[1][1, :]
best_performing_NmKGE_ensemble = sorted_NmKGEs[best_NmKGE_id][1, :]
best_performing_NmKGE_name = df_names[best_NmKGE_id]
best_NmKGE_ensemble = best_performing_NmKGE_ensemble[:Approach]
println("Best Performing: ", best_performing_NmKGE_name, " ", best_NmKGE_ensemble, " Score: ", best_NmKGE_score)

base_Qe = CSV.File(joinpath(ensemble_res_path, "baseline_outflows.csv")) |> DataFrame
state_Qe = CSV.File(joinpath(ensemble_res_path, "state_outflows.csv")) |> DataFrame
base_base_Qe = CSV.File(joinpath(ensemble_res_path, "base_base_outflows.csv")) |> DataFrame
base_state_Qe = CSV.File(joinpath(ensemble_res_path, "base_state_outflows.csv")) |> DataFrame
state_state_Qe = CSV.File(joinpath(ensemble_res_path, "state_state_outflows.csv")) |> DataFrame

outflows = [base_Qe, state_Qe, base_base_Qe, base_state_Qe, state_state_Qe]
relabel_approaches(outflows, cols=true)


best_bl_col = create_readable_name(best_performing_baseline[:Approach])
best_bl_display = "Baseline ($best_bl_col)"

best_NmKGE_col = create_readable_name(best_NmKGE_ensemble)
best_NmKGE_display = "$(best_performing_NmKGE_name) ($(best_NmKGE_col))"

base_NmKGE_Q = base_Qe[:, "NKGE'"]
base_NmKGE_calib = base_NmKGE_Q[calib_ts+1:CALIB_LN]
base_NmKGE_valid = base_NmKGE_Q[CALIB_LN+1:end]

best_bl_Q = base_Qe[:, best_bl_col]
best_bl_calib = best_bl_Q[calib_ts+1:CALIB_LN]
best_bl_valid = best_bl_Q[CALIB_LN+1:end]

best_NmKGE_Q = outflows[best_NmKGE_id][:, best_NmKGE_col]
best_NmKGE_Q_calib = best_NmKGE_Q[calib_ts+1:CALIB_LN]
best_NmKGE_Q_valid = best_NmKGE_Q[CALIB_LN+1:end]

base_NmKGE_label = create_NmKGE_label("Baseline NKGE'", Qo, base_NmKGE_Q)
best_bl_label = create_NmKGE_label(best_bl_display, Qo, best_bl_Q)
best_NmKGE_label = create_NmKGE_label(best_NmKGE_display, Qo, best_NmKGE_Q)

res_qq = qqplot(Qo_valid, base_NmKGE_valid, label=base_NmKGE_label, xlabel="Observed [ML]", ylabel="Modeled Flow [ML]")
qqplot!(Qo_valid, best_bl_valid, label=best_bl_label, xlabel="Observed [ML]", ylabel="Modeled Flow [ML]", alpha=0.5)
qqplot!(Qo_valid, best_NmKGE_Q_valid, alpha=0.4, label=best_NmKGE_label)
savefig("$(ensemble_figs)comparison_baseline_best_ensembles_valid_qqplot.png")

res_scatter = scatter(Qo_valid, base_NmKGE_valid .- Qo_valid, label=base_NmKGE_label, xlabel="Observed [ML]", ylabel="Modeled - Observed [ML]")
scatter!(Qo_valid, best_bl_valid - Qo_valid, alpha=0.5, label=best_bl_label)
scatter!(Qo_valid, best_NmKGE_Q_valid, alpha=0.4, label=best_NmKGE_label)
savefig("$(ensemble_figs)comparison_baseline_best_ensembles_valid_residuals.png")

#######

# Read in state params
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

# Use baseline NmKGE catchment conditions 
calib_low, calib_mid, calib_high = separate_state_periods(active_params, calib_ts+1, CALIB_LN)
valid_low, valid_mid, valid_high = separate_state_periods(active_params, CALIB_LN+1)

# Show each state
base_dry_valid = base_NmKGE_valid[valid_high]
bl_dry_valid = best_bl_valid[valid_high]
NmKGE_dry_valid = best_NmKGE_Q_valid[valid_high]
Qo_dry_valid = Qo_valid[valid_high]
dry_scatter = scatter(Qo_dry_valid, base_dry_valid .- Qo_dry_valid, label=base_NmKGE_label, xlabel="Observed [ML]", ylabel="Modeled [ML]")

base_NmKGE_mid_label = create_NmKGE_label("Baseline NKGE'", Qo, base_NmKGE_Q, calib_mid, valid_mid)
best_bl_mid_label = create_NmKGE_label(best_bl_display, Qo, best_bl_Q, calib_mid, valid_mid)
best_NmKGE_mid_label = create_NmKGE_label(best_NmKGE_display, Qo, best_NmKGE_Q, calib_mid, valid_mid)
base_mid_qq = qqplot(Qo_valid[valid_mid], base_NmKGE_valid[valid_mid], label=base_NmKGE_mid_label, xlabel="Observed [ML]", ylabel="Modeled Flow [ML]")
bl_mid_qq = qqplot!(Qo_valid[valid_mid], best_bl_valid[valid_mid], label=best_bl_mid_label, alpha=0.3)
NmKGE_mid_qq = qqplot!(Qo_valid[valid_mid], best_NmKGE_Q_valid[valid_mid], label=best_NmKGE_mid_label, alpha=0.3)
savefig("$(ensemble_figs)comparison_usual_valid_qq.png")

base_NmKGE_dry_label = create_NmKGE_label("Baseline NKGE'", Qo, base_NmKGE_Q, calib_high, valid_high)
bl_NmKGE_dry_label = create_NmKGE_label(best_bl_display, Qo, best_bl_Q, calib_high, valid_high)
best_NmKGE_dry_label = create_NmKGE_label(best_NmKGE_display, Qo, best_NmKGE_Q, calib_high, valid_high)
base_dry_qq = qqplot(Qo_dry_valid, base_dry_valid, label=base_NmKGE_dry_label, xlabel="Observed [ML]", ylabel="Modeled Flow [ML]", legend=:topleft, log=true)
bl_dry_qq = qqplot!(Qo_dry_valid, bl_dry_valid, label=bl_NmKGE_dry_label, alpha=0.3)
NmKGE_dry_qq = qqplot!(Qo_dry_valid, NmKGE_dry_valid, label=best_NmKGE_dry_label, alpha=0.3)
savefig("$(ensemble_figs)comparison_dry_valid_qq.png")

combined_state_qq = plot(
    plot(base_mid_qq, title="Usual", legend=:bottomright),
    plot(base_dry_qq, title="Dry", legend=:bottomright),
    layout=(1, 2),
    margin=5Plots.mm,
    size=(1000,400),
    markerstrokewidth=0,
    xaxis=:log,
    yaxis=:log
)
savefig(combined_state_qq, "$(ensemble_figs)combined_state_qq.png")
#####

base_NmKGE_label = create_NmKGE_label("Baseline NKGE'", Qo, base_NmKGE_Q)
best_NmKGE_label = create_NmKGE_label(best_NmKGE_display, Qo, best_NmKGE_Q)
res_scatter2 = scatter(Qo_valid, base_NmKGE_valid - Qo_valid, label=base_NmKGE_label, 
                       xlabel="Observed [ML]", ylabel="Modeled - Observed [ML]")
scatter!(Qo_valid, best_bl_valid - Qo_valid, alpha=0.5, label=best_bl_label)
scatter!(Qo_valid, best_NmKGE_Q_valid - Qo_valid, alpha=0.4, label=best_NmKGE_label)

combined_res = plot(
    plot(res_qq, xaxis=:log, yaxis=:log, legend=:bottomright),
    plot(res_scatter2, legend=false),
    layout=(1, 2),
    margin=5Plots.mm,
    size=(1000,400),
    #link=:y,
    markerstrokewidth=0
)
savefig(combined_res, "$(ensemble_figs)combined_residuals.png")

layout_4panel = @layout [a b; c]
combined_4panel = plot(
    plot(res_qq, title="Q-Q Plot", xaxis=:log, yaxis=:log, legend=:bottomright),
    plot(res_scatter2, title="Residuals", legend=false),
    combined_state_qq,
    layout=layout_4panel,
    margin=5Plots.mm,
    size=(1000,860),
    markerstrokewidth=0
)
savefig(combined_4panel, "$(ensemble_figs)combined_residuals_4panel.png")


# QQ of Best performing instance out of each approach
base_NmKGE_label = create_NmKGE_label("Baseline NKGE'", Qo, base_NmKGE_Q)
best_NmKGE_qq = qqplot(Qo_valid, base_Qe[:, "NKGE'"][CALIB_LN+1:end], label=base_NmKGE_label,
                       legend=:bottomright, xaxis=:log, yaxis=:log)
flow_data = [base_Qe, state_Qe, base_base_Qe, base_state_Qe, state_state_Qe]
for (df_name, flow, df) in zip(df_names, flow_data, sorted_NmKGEs)
    best_name = create_readable_name(df[1, :Approach])

    if (df_name == "Baseline") && (best_name == "NKGE'")
        continue
    end

    best_NmKGE_label = create_NmKGE_label(best_name, Qo, flow[:, best_name])
    qqplot!(Qo_valid, flow[CALIB_LN+1:end, best_name], label="$(df_name) [$(best_NmKGE_label)]",
            alpha=0.3, legend=:bottomright, markerstrokewidth=0)
end
savefig(best_NmKGE_qq, "$(ensemble_figs)best_NmKGE_qq.png")

# QQ based on A2k score
# best_A2k_qq = qqplot(Qo_valid, base_Qe[:, "NKGE'"][CALIB_LN+1:end], label="Baseline NKGE'", legend=:topleft)
# for (df_name, flow, df) in zip(df_names, flow_data, sorted_A2ks)
#     best_name = create_readable_name(df[1, :Approach])

#     if (df_name == "Baseline") && (best_name == "NKGE'")
#         continue
#     end

#     ME_calib = round(mean(Qo_calib - flow[calib_ts+1:CALIB_LN, best_name]), digits=2)
#     ME_valid = round(mean(Qo_valid - flow[CALIB_LN+1:end, best_name]), digits=2)
#     qqplot!(Qo_valid, flow[CALIB_LN+1:end, best_name], label="$(df_name) ($(best_name); C: $(ME_calib) | V: $(ME_valid))", 
#             alpha=0.3, legend=:topleft, markerstrokewidth=0)
# end
# savefig(best_A2k_qq, "$(ensemble_figs)best_A2k_qq.png")

# Create violin plots
bl_NmKGE_calib, bl_NmKGE_valid = split_calib_valid_cols(sorted_NmKGEs[1])  # baseline NmKGEs
state_NmKGE_calib, state_NmKGE_valid = split_calib_valid_cols(sorted_NmKGEs[2])  # state-based NmKGEs
base_base_NmKGE_calib, base_base_NmKGE_valid = split_calib_valid_cols(sorted_NmKGEs[3])  # base-base NmKGEs
base_state_NmKGE_calib, base_state_NmKGE_valid = split_calib_valid_cols(sorted_NmKGEs[4])  # base-state NmKGEs
state_state_NmKGE_calib, state_state_NmKGE_valid = split_calib_valid_cols(sorted_NmKGEs[5])  # state-state NmKGEs

baseline_vplot = create_violin(bl_NmKGE_calib, bl_NmKGE_valid, "Baseline")
savefig(baseline_vplot, "$(ensemble_figs)baseline_NmKGE_results.png")

state_based_vplot = create_violin(state_NmKGE_calib, state_NmKGE_valid, "State-based")
savefig(state_based_vplot, "$(ensemble_figs)state-based_NmKGE_results.png")

base_base_vplot = create_violin(base_base_NmKGE_calib, base_base_NmKGE_valid, "Base-Base")
savefig(base_base_vplot, "$(ensemble_figs)base-base_NmKGE_results.png")

base_state_vplot = create_violin(base_state_NmKGE_calib, base_state_NmKGE_valid, "Base-State")
savefig(base_state_vplot, "$(ensemble_figs)base-state_NmKGE_results.png")

state_state_vplot = create_violin(state_state_NmKGE_calib, state_state_NmKGE_valid, "State-State")
savefig(state_state_vplot, "$(ensemble_figs)state-state_NmKGE_results.png")


combined_violin = plot(
    baseline_vplot,
    plot(state_based_vplot, ylabel="", legend=false),
    plot(state_state_vplot, ylabel="", legend=false),
    layout=(1, 3),
    margin=5Plots.mm,
    size=(1200,400),
    link=:y
)
savefig(combined_violin, "$(ensemble_figs)combined_violin.png")


combined_violin_all = plot(
    baseline_vplot,
    plot(state_based_vplot, ylabel="", legend=false),
    plot(base_base_vplot, ylabel="", legend=false),
    plot(base_state_vplot, ylabel="", legend=false),
    plot(state_state_vplot, ylabel="", legend=false),
    layout=(1, 5),
    margin=5Plots.mm,
    size=(2000,400),
    link=:y
)
savefig(combined_violin_all, "$(ensemble_figs)combined_violin_all.png")
