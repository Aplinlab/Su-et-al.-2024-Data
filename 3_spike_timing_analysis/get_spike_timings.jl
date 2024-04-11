# Packages
using CSV
using JLD2
using FileIO
using DataFrames
using XLSX
using Dates

# Working directory and dependencies
cd(dirname(@__FILE__))
include("../src/spike_sort_main.jl")
include("../src/spike_sort_plotting.jl")
include("../src/ts_helper.jl")       #! Must remain after other dependencies



df_input = DataFrame(XLSX.readtable("../saved_files_info.xlsx", "Templating"))

l1_spikes = []
m3_spikes = []
m2_spikes = []
m1_spikes = []

pre_window = 7
post_window = 7
fs=30000

for loading_output in eachrow(df_input)
    filename = loading_output["Filename"]
    animal_ID = loading_output["Animal"]
    pos_sheet = loading_output["Position"]
    cell_desc = loading_output["Cell"]
    condition = loading_output["Condition"]
    recovery_time = parse(Float64, loading_output["Recovery Time"])

    # if recovery_time != -1
    #     continue
    # end

    regexParam = r"(?P<position>\d+) (?P<amp>\d\.?\d*) (?P<block_amp>-?\d+) (?P<stim_type>\w{2}) (?P<files_num>\[(?:\d+(?:, )*)+]) (?P<channels>\[(?:\d+(?:, )*)+])"
    parameters = match(regexParam, filename)

    pos_regex = String(parameters["position"])
    if pos_sheet != pos_regex
        error("Positions do not match for $animal_ID $filename")
    end
    pos = parse(Int, pos_sheet)
    block_amp = String(parameters["block_amp"])
    stim_type = String(parameters["stim_type"])

    if block_amp != "0"
        continue
    end

    println("$animal_ID $pos_sheet $cell_desc $condition ($filename)")

    if stim_type == "L1"
        pre_trial_len = round(Int, 1.0 * fs)
        post_trial_len = round(Int,5.0 * fs)
        WOI2 = [500, 2000]
    elseif stim_type == "M3"
        pre_trial_len = round(Int, 0.5 * fs)
        post_trial_len = round(Int,1.5 * fs)
        WOI2 = [50, 400]
    elseif stim_type == "M2"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI2 = [0, 400]
    elseif stim_type == "M1"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI2 = [0, 50]
    else
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI2 = [0, 100]
    end
    WOI1 = stim_type == "L1" ? [-500, 0] : [WOI2[1] - WOI2[2], 0]
    WOI1 .*= (fs/1000)
    WOI2 .*= (fs/1000)
    
    file = load("../1_main_analysis/loading_output/$animal_ID/$filename.jld2")

    signal_matrix, pk_bkgnds, file_lengths, stim_inds, cum_stim_inds = def_template_gen_vars(file)

    pk_inds, _waveforms, _pk_heights = detect_spikes(signal_matrix, pk_bkgnds, stim_inds, [pre_window, post_window], file_lengths = file_lengths)
    pks_plot, _ass_plot = assign_to_trials(cum_stim_inds, pk_inds, zeros(length(pk_inds)), pre_trial_len, post_trial_len)

    pks_plot_zeroed = Vector()
    for pk_trial in pks_plot
        pk_trial_zeroed = pk_trial .- pre_trial_len
        push!(pks_plot_zeroed, pk_trial_zeroed)
    end

    if stim_type == "L1"
        for trial in pks_plot_zeroed
            append!(l1_spikes, spike for spike in trial)
        end
    elseif stim_type == "M3"
        for trial in pks_plot_zeroed
            append!(m3_spikes, spike for spike in trial)
        end
    elseif stim_type == "M2"
        for trial in pks_plot_zeroed
            append!(m2_spikes, spike for spike in trial)
        end
    elseif stim_type == "M1"
        for trial in pks_plot_zeroed
            append!(m1_spikes, spike for spike in trial)
        end
    end
end

max_length = max(length(l1_spikes), length(m3_spikes), length(m2_spikes), length(m1_spikes))

append!(l1_spikes, Vector{Any}(nothing, max_length - length(l1_spikes)))
append!(m3_spikes, Vector{Any}(nothing, max_length - length(m3_spikes)))
append!(m2_spikes, Vector{Any}(nothing, max_length - length(m2_spikes)))
append!(m1_spikes, Vector{Any}(nothing, max_length - length(m1_spikes)))

df_output = DataFrame((laser = l1_spikes, press = m3_spikes, tactile = m1_spikes, proprioceptive = m2_spikes))

date = Dates.today()
CSV.write("spike_timings $date.csv", df_output; transform=(col, val) -> something(val, ""))
println("All files checked.\nOutput saved.")