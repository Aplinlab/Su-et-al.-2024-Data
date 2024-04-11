version = v"4.3.1"
# Sorts all data files into unit assignments using the previously determined templates. Sorted data
# is saved into /spikesorting_output in the .JLD2 format.

# * Changelog
#   Moved outputTable definition and inner for loop to csvExportBatch
#   Save output up to that point
#   Updated to TSHelperTS4.4.0
#   Updated def_spikesorting_vars() name

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
using PyPlot; pygui(true)

using Statistics #std
using StatsBase #mode

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
# Spikesorting Parameters
fs=30000
onems = round(Int, fs/1000)
pre_window = 7  # peak extraction
post_window = 7 # peak extraction

# Function Parameters
samples_on_x=false  # what units to use for x-axis of raster plots - true for samples or false for milliseconds

################################################################################

df = DataFrame(XLSX.readtable("../saved_files_info.xlsx", "Templating"))
df_sexes = DataFrame(XLSX.readtable("../saved_files_info.xlsx", "Sexes"))

numerals = [Char(x + '0') for x in 0:9]
spikesortedFilesCount = 0

global error_log = ""

for comparison in eachrow(df)
    templateFilename = comparison["Template"]
    comparisonFilename = comparison["Filename"]
    animalID = comparison["Animal"]
    pos = parse(Int, comparison["Position"])
    cell_desc = comparison["Cell"]
    condition = comparison["Condition"]
    recovery_time = parse(Float64, comparison["Recovery Time"])

    println("$animalID $pos $cell_desc $condition ($comparisonFilename)")

    templateFile = load("../spikesorting_templates/$animalID/$templateFilename.jld2")
    templateLoadingFile = load("../loadingOutput/$animalID/$templateFilename.jld2")
    comparisonLoadingFile = load("../loadingOutput/$animalID/$comparisonFilename.jld2")

    wave_signal_matrix, wave_pk_bkgnds, wave_file_lengths, wave_stim_inds, wave_cum_stim_inds, 
        template_waveforms, template_templates, template_assigns, template_thresholds, 
        stim_amp, idc_amp, stim_code, channels, recording_filenames, recording_files_num, template_filenames, template_files_num, good_assigns = def_spikesorting_vars(comparisonLoadingFile, templateLoadingFile, templateFile)
    
    treatment = rstrip(animalID, numerals)
    cell_type = rstrip(cell_desc, numerals)
    trials_count = stim_code == "L1" ? 2 : 50
    sex = length(df_sexes[occursin.(df_sexes.Animal, animalID), "Sex"]) == 1 ? 
        df_sexes[occursin.(df_sexes.Animal, animalID), "Sex"][1] : 
        error("more than one sex is listed for $ref_animalID")
    stim_type = stim_amp == 0.0 ? "spon" : 
        stim_code == "M1" ? "shaker" : 
        stim_code == "M2" ? "motor" : 
        stim_code == "M3" ? "pinch" : 
        stim_code == "L1" ? "laser" : 
        stim_code == "E1" ? "estim" : 
        error("reference stim type is not recognised ($ref_filename)")
    exp_phase = idc_amp != 0 ? "idc" : 
        recovery_time > 0.0 ? "recovery" : 
        recovery_time == -1.0 ? "control" : error("invalid value for recoveryTime ($ref_filename)")
    
    regexTime = r"(?P<hr>\d{2})(?P<min>\d{2})(?P<sec>\d{2}).rhs"
    recordingSearch = match(regexTime, recording_filenames[1])
    recording_hours = parse(Int, recordingSearch["hr"])
    recording_minutes = parse(Int, recordingSearch["min"])
    recording_seconds = parse(Int, recordingSearch["sec"])
    templateSearch = match(regexTime, template_filenames[1])
    template_hours = parse(Int, templateSearch["hr"])
    template_minutes = parse(Int, templateSearch["min"])
    template_seconds = parse(Int, templateSearch["sec"])
    recording_time = recording_hours*3600 + recording_minutes*60 + recording_seconds
    template_time = template_hours*3600 + template_minutes*60 + template_seconds
    delta_t = recording_time - template_time

    if stim_code == "L1"
        pre_trial_len = round(Int, 1.0 * fs)
        post_trial_len = round(Int,5.0 * fs)
        WOI2 = [800, 1500]
    elseif stim_code == "M3"
        pre_trial_len = round(Int, 0.5 * fs)
        post_trial_len = round(Int,1.5 * fs)
        WOI2 = [0, 500]
    elseif stim_code == "M2"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI2 = [0, 400]
    elseif stim_code == "M1"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI2 = [0, 40]
    else
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI2 = [0, 100]
    end
    WOI1 = [WOI2[1] - WOI2[2], 0]
    num_assigns = maximum(template_assigns) + 1

    ################################################################################

    # * template application
    # TODO:    Check modifications by Jason are ok (added x_min and merged two function versions) 

    #' define data file to apply templates
    wave_pk_inds, wave_waveforms, _wave_pk_heights = detect_spikes(
        wave_signal_matrix,
        wave_pk_bkgnds,
        wave_stim_inds,
        [pre_window, post_window],
        file_lengths = wave_file_lengths
    )

    #' Assigning spikes to templates
    wave_assignment = assign_points_to_clusters(wave_waveforms, template_templates, template_thresholds)

    #Plots of assigned waveforms

    # Combine all stim_inds
    pks_plot, ass_plot = assign_to_trials(wave_cum_stim_inds, wave_pk_inds, wave_assignment, pre_trial_len, post_trial_len)

    description = string(animalID*"-"*string(pos)*"-"*cell_desc*" "*condition*"-"*comparisonFilename*"-"*templateFilename)

    trials_actual = size(ass_plot, 1)
    if trials_actual < trials_count
        global error_log *= "$description - only $trials_actual trials\n\n"
        println("$description - only $trials_actual trials (skipped)\n\n")
        continue
    elseif trials_actual > trials_count
        @info "$trials_actual trials ($description)"
    end

    @save "../spikesorted_data/$description.jld2" animalID treatment sex pos channels cell_desc cell_type stim_type stim_code stim_amp idc_amp recovery_time exp_phase WOI1 WOI2 trials_count recording_time template_time delta_t recording_filenames recording_files_num template_filenames template_files_num pks_plot ass_plot pre_trial_len post_trial_len onems pre_window post_window num_assigns good_assigns
    
    ################################################################################

    global spikesortedFilesCount += 1
    println("Spikesorting complete: $condition\n")
end

println("Batch of $spikesortedFilesCount files spikesorted")

if !isempty(error_log)
    now = Dates.now()
    log_file = open("../error_logs/main_analysis.txt", "a")
    println(log_file, "---Spikesorting $now---\n\n$error_log")
    close(log_file)
end