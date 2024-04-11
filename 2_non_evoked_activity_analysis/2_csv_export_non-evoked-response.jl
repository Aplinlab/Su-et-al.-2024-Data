version = v"4.1.3"
# Exports spikesortingOutput to a CSV after counting spikes.

# * Changelog

# ! Errors
#   Empty assignment workaround may occasionally produce an error (needs more documentation)

#   TODO List

# ? Suggestions & Queries



# Packages
using DSP
using Glob
using LinearAlgebra
using MultivariateStats
using Printf
using PyPlot; pygui(true) # for Display

using Statistics #std
using StatsBase #mode

using CSV
using JLD2
using FileIO
using PyCall
using DataFrames
using XLSX
using Query
using Dates
CheckButtons = matplotlib.widgets.CheckButtons

# Working directory and dependencies
cd(dirname(@__FILE__))
include("../src/spike_sort_main.jl")
include("../src/spike_sort_plotting.jl")
include("../src/ts_helper.jl")       #! Must remain after other dependencies



# * Parameters
# Function Parameters
ignore_unassigned = true
samples_on_x=false  # what units to use for x-axis of raster plots - true for samples or false for milliseconds

detection_threshold = 3.0

################################################################################

df = DataFrame(XLSX.readtable("../saved_files_info.xlsx", "Loading List (pre-stimulus)"))
df_sexes = DataFrame(XLSX.readtable("../saved_files_info.xlsx", "Sexes"))

output_table = DataFrame(
    filename = [],
    animal_ID = [],
    condition = [],
    sex = [],
    position = [],
    stim_type = [],
    stim_amp = [],
    block_amp = [],
    sample_rate = [],

    channel = [],

    spike_count = [],
    duration = [],
    spike_frequency = [],

    voltage_mean = [],
    voltage_std = []
)

global error_log = ""

for loadinginput in eachrow(df)
    filename = loadinginput["rhs_name"]
    animal_ID = loadinginput["animal"]
    condition = animal_ID[1:3]
    position = loadinginput["electrode_pos"]
    stim_type = loadinginput["stim_type"]
    stim_amp = loadinginput["amplitude"]
    block_amp = loadinginput["block_amp"]

    println(filename)

    sex = length(df_sexes[occursin.(df_sexes.Animal, animal_ID), "Sex"]) == 1 ? 
        df_sexes[occursin.(df_sexes.Animal, animal_ID), "Sex"][1] : 
        error("No sex or more than one sex is listed for $animal_ID")

    if isfile("loading_output_non-evoked-response/$animal_ID/$filename.jld2")
        file = load("loading_output_non-evoked-response/$animal_ID/$filename.jld2")
    else
        global error_log *= "$filename - loading output missing\n\n"
        println("$filename - loading output missing (skipped)\n\n")
        continue
    end
    
    # Error checking
    animal_ID_check = file["animalID"]
    condition_check = file["condition"]
    position_check = file["pos"]
    stim_type_check = file["stim_type"]
    stim_amp_check = file["stim_amp"]
    block_amp_check = file["block_amp"]

    if animal_ID != animal_ID_check
        error("Animal IDs do not match")
    elseif condition != condition_check
        error("Conditions do not match")
    elseif position != position_check
        error("Positions do not match")
    elseif stim_type != stim_type_check
        error("Stimulation types do not match")
    elseif stim_amp != stim_amp_check
        error("Stimulation amplitudes do not match")
    elseif block_amp != block_amp_check
        error("Block amplitudes do not match")
    end

    ################################################################################

    signal_matrix = file["signal_matrix"]   # 32x[file_length] matrix
    duration = file["duration"]             # float
    file_length = file["file_length"]       # int
    sample_rate = file["fs"]                # float
    stim_ind = file["stim_ind"]             # int

    for (channel, signal) in enumerate(eachrow(signal_matrix))
        voltage_mean = mean(signal)
        voltage_std = std(signal)

        spike_threshold = voltage_mean + (detection_threshold * voltage_std)
        spikes, spike_properties = findpeaks1d(signal, height=spike_threshold, distance=6)

        spike_count = length(spikes)
        spike_frequency = spike_count / duration

        push!(output_table, (
            filename,
            animal_ID,
            condition,
            sex,
            position,
            stim_type,
            stim_amp,
            block_amp,
            sample_rate,

            channel,

            spike_count,
            duration,
            spike_frequency,

            voltage_mean,
            voltage_std
        ))
        println("Channel $channel counted")
    end
    println("All channels counted\n\n")
end

date = Dates.today()
CSV.write("full_results non-evoked-response $date.csv", output_table)
println("CSV Exporting Complete")

if !isempty(error_log)
    now = Dates.now()
    log_file = open("../error_logs/non-evoked.txt", "a")
    println(log_file, "---CSV Exporting $now---\n\n$error_log")
    close(log_file)
end