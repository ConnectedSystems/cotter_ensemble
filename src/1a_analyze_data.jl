include("_common.jl")

using CategoricalArrays, StatsPlots


climate_fig_path = joinpath(FIG_PATH, "climate") * "/"
mkpath(climate_fig_path)

# Change plot theme for publication
theme(:default, grid=false)

dataset = copy(FULL_DATASET)
data_years = Dates.year.(dataset.Date)
start_year = minimum(data_years)
end_year = maximum(data_years)
half_decades = cut(data_years, start_year:5:end_year, extend=true)
annual = cut(data_years, start_year:1:end_year, extend=true)

dataset.half_decade_set = half_decades
dataset.annual = annual

half_decade_dataset = groupby(dataset, :half_decade_set)
annual_dataset = groupby(dataset, :annual)


# Plot historic rainfall
thresholds = quantile(FULL_DATASET[:, "410730_P"], [0.6, 0.95])
bar(FULL_DATASET.Date, FULL_DATASET[:, "410730_P"],
	title="Observed Rainfall",
	label="",
	xlabel="Year",
	ylabel="Rainfall [mm]")
hline!(thresholds, label="Bin thresholds")
savefig("$(climate_fig_path)rainfall_record.png")


rain_sym = Symbol("410730_P")
temp_sym = Symbol("410730_max_T")
flow_sym = Symbol("410730_Q")

climate_annual_stats = combine(annual_dataset,
						rain_sym => mean,
						rain_sym => std,
						rain_sym => maximum,
						rain_sym => sum,
						temp_sym => mean,
						temp_sym => std,
						temp_sym => maximum,
						temp_sym => sum)

# Plot annual rain
bar(start_year:end_year, climate_annual_stats[:, "410730_P_sum"],
	title="Observed Annual Rainfall",
	label="",
	xlabel="Year",
	ylabel="Rainfall [mm]")
savefig("$(climate_fig_path)rainfall.png")

climate_half_decade_stats = combine(half_decade_dataset,
							rain_sym => mean,
							rain_sym => std,
							rain_sym => maximum,
							rain_sym => sum)

flow_annual_stats = combine(annual_dataset,
						flow_sym => mean,
						flow_sym => std,
						flow_sym => maximum,
						flow_sym => sum)

flow_half_decade_stats = combine(half_decade_dataset,
							flow_sym => mean,
							flow_sym => std,
							flow_sym => maximum,
							flow_sym => sum)

# flow_half_decade_stats

bar(start_year+1:end_year, climate_annual_stats[:, "410730_P_sum"], label="", 
	title="Total Annual Rainfall",
	xlabel="Year",
	ylabel="Rainfall [mm]")
savefig("$(climate_fig_path)rainfall_sum.png")

p_total = climate_annual_stats[:, "410730_P_sum"]
violin(["Calibration"], p_total[1:30], label="", title="Annual Rainfall", ylabel="Total [mm]")
violin!(["Validation"], p_total[31:end], label="")
dotplot!(["Calibration"], p_total[1:30],
		 marker=(:black, stroke(0)),
		 label="")
dotplot!(["Validation"], p_total[31:end],
		 marker=(:black, stroke(0)),
		 label="")
savefig("$(climate_fig_path)rainfall_comparison.png")


# Average temperatures
plot(start_year+1:end_year, climate_annual_stats[:, "410730_max_T_mean"], color="red",
	 label="", 
	 title="Average Annual Max Temperature",
	 xlabel="Year",
	 ylabel="Average Max Temperature [°C]")
savefig("$(climate_fig_path)average_max_T.png")

T_total = climate_annual_stats[:, "410730_max_T_mean"]
violin(["Calibration"], T_total[1:30], label="", title="Average Max Temperature", ylabel="Average Degrees [°C]")
violin!(["Validation"], T_total[31:end], label="")
dotplot!(["Calibration"], T_total[1:30],
		 marker=(:black, stroke(0)),
		 label="")
dotplot!(["Validation"], T_total[31:end],
		 marker=(:black, stroke(0)),
		 label="")
savefig("$(climate_fig_path)temperature_comparison.png")



# Coefficient of variation
climate_annual_stats.rainfall_variability = climate_annual_stats[:, "410730_P_std"] ./ climate_annual_stats[:, "410730_P_mean"]

p_vari = climate_annual_stats.rainfall_variability
violin(["Calibration"], p_vari[1:30], label="", title="Variability of Rainfall", ylabel="Variability [CV]")
violin!(["Validation"], p_vari[31:end], label="")
dotplot!(["Calibration"], p_vari[1:30],
		 marker=(:black, stroke(0)),
		 label="")
dotplot!(["Validation"], p_vari[31:end],
		 marker=(:black, stroke(0)),
		 label="")
savefig("$(climate_fig_path)rainfall_variability_comparison.png")

bar(start_year:end_year, flow_annual_stats[:, "410730_Q_sum"] ./ 1000,
		xlabel="Year", ylabel="Annual Flow [GL]", label="",
		title="Annual Flow [GL]"
)

savefig("$(climate_fig_path)annual_flow.png")

flow_annual_stats.flow_variability = (flow_annual_stats[:, "410730_Q_std"] ./ flow_annual_stats[:, "410730_Q_mean"])

bar(start_year:end_year, flow_annual_stats.flow_variability,
	xlabel="Year", ylabel="Flow variability [CV]", label="",
	title="Variability of Annual Streamflow"
)

# Qualitatively determined separation between calib start/end
c_end = start_year + 30
v_start = start_year + 31
vline!([v_start], label="")
annotate!((start_year+0.5, 1.6, text("Calibration", font(11), :left)))
annotate!((v_start+0.5, 1.6, text("Validation", font(11), :left)))
savefig("$(climate_fig_path)flow_variability.png")


flow_col = Symbol("flow_variability")
violin(["Calibration"], flow_annual_stats[1:30, flow_col],
		label="", ylabel="Variability of Flow [CV]",
		title="Annual Flow\nVariability")
violin!(["Validation"], flow_annual_stats[31:end, flow_col],
		label="")
dotplot!(["Calibration"], flow_annual_stats[1:30, flow_col],
		 marker=(:black, stroke(0)),
		 label="")
dotplot!(["Validation"], flow_annual_stats[31:end, flow_col],
		 marker=(:black, stroke(0)),
		 label="")
savefig("$(climate_fig_path)flow_variability_comparison.png")



"""
Köppen aridity index

See:
https://doi.org/10.1016/j.palaeo.2013.05.008

AIₖ : arid: < 40; humid: 40–160; perhumid: > 160
"""
AI_k(P, T) = P ./ (T .+ 33.0)


AP = climate_annual_stats[:, "410730_P_sum"]
AT = climate_annual_stats[:, "410730_max_T_mean"]
aridity = AI_k(AP, AT)

bar(start_year:end_year, aridity,
	label="",
	xlabel="Year",
	ylabel="\$AI_{k}\$",
	title="Köppen Aridity Index"
)

vline!([v_start], label="", linewidth=2.0)
annotate!((start_year, 34.0, text("Calibration", font(11), :left)))
annotate!((v_start+0.5, 34.0, text("Validation", font(11), :left)))
savefig("$(climate_fig_path)aridity.png")

violin(["Calibration"], aridity[1:30], label="",
		ylabel="\$AI_{k}\$",
		title="Aridity Characteristics")
violin!(["Validation"], aridity[31:end], label="")
dotplot!(["Calibration"], aridity[1:30],
		 marker=(:black, stroke(0)),
		 label="")
dotplot!(["Validation"], aridity[31:end],
		 marker=(:black, stroke(0)),
		 label="")
savefig("$(climate_fig_path)aridity_comparison.png")

@info "Timeframe:" start_year, end_year
@info minimum(dataset.Date), maximum(dataset.Date)


flow_col = Symbol("410730_Q_sum")
violin(["Calibration"], flow_annual_stats[1:30, flow_col] ./ 1000,
		label="", ylabel="Annual Flow [GL]",
		title="Annual Flow\nCharacteristics")
violin!(["Validation"], flow_annual_stats[31:end, flow_col] ./ 1000,
		label="")
dotplot!(["Calibration"], flow_annual_stats[1:30, flow_col] ./1000,
		 marker=(:black, stroke(0)),
		 label="")
dotplot!(["Validation"], flow_annual_stats[31:end, flow_col] ./1000,
		 marker=(:black, stroke(0)),
		 label="")
savefig("$(climate_fig_path)flow_comparison.png")


describe(flow_annual_stats[1:30, Symbol("410730_Q_sum")] ./ 1000)
# Summary Stats:
# Length:         30        
# Missing Count:  0
# Mean:           48.433578 
# Minimum:        6.393224  
# 1st Quartile:   31.157348 
# Median:         46.762521 
# 3rd Quartile:   62.423117 
# Maximum:        110.813315
# Type:           Float64


describe(flow_annual_stats[31:end, Symbol("410730_Q_sum")] ./ 1000)
# Summary Stats:
# Length:         21
# Missing Count:  0
# Mean:           35.346256
# Minimum:        8.673252
# 1st Quartile:   21.041727
# Median:         35.154885
# 3rd Quartile:   51.299297
# Maximum:        67.310798
# Type:           Float64

# Plot streamflow
plot(CALIB_DATES, FULL_DATASET[1:CALIB_LN, "410730_Q"], label="Calibration", color="darkblue", title="Historic Streamflow")
plot!(VALID_DATES, FULL_DATASET[CALIB_LN+1:end, "410730_Q"], label="Validation", color="lightblue", ylabel="Streamflow [ML]")
savefig("$(climate_fig_path)historic_flow.png")
