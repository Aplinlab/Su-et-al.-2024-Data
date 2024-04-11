version = v"4.2.2"
# Batch version of loading.
#! THIS VERSION ONLY WORKS IF YOU ALREADY KNOW WHAT PARAMETERS TO FEED INTO LOADING USING THE LOOKUP SPREADSHEET

# * Changelog
#   Completely different for loop
#   Formatting/order changes
#   More error messages
#   Commented out some obsolete/plotting-only variables
#   Define trials_per_file for laser according to files_num
#   See loading for additional changelog

# ! Errors

#   TODO List
#   Move regex out of loop
#   Clean up for upload

# ? Suggestions & Queries



# ! before loading new files, ensure you kill the julia session to remove all global variables
# Packages
using DataFrames
using Query
using XLSX
using PyPlot; pygui(true)
using Statistics # std
using StatsBase # mode
using JLD2
using Dates

# Working directory and dependencies
cd(dirname(@__FILE__))
include("../src/IntanReader.jl")
include("../src/preprocessing_module.jl")        # Preprocessing Tools
include("../src/response_quantification.jl")     # Response quantification Tools
include("../src/RHS2000_MEA_config.jl")



# Filtering data
# ' ## Getting filenames
df = DataFrame(XLSX.readtable("../saved_files_info.xlsx", "Loading List (pre-stimulus)", stop_in_empty_row=false))

datafile_location = "../raw_data/"

logic_high = true               # whether the triggers are peaks or valleys
plot_all = true                 # false will produce a single waveform plot (for use as a figure)

# Stim Artefact
art_pre = 60            # use to control pre-stim blanking period (30 = 1 ms)
art_post = 60           # use to control post-stim blanking period
art_height = 500        # used to threshold artefact blanking (when using artefact blank only)
trig_blank = true       # To control whether the artefact itself is used to blank, or the digital trigger timing (true = trigger)

# Filtering
bpass = 300
bstop = 5000

# Peak Detection
intra_chan_window = 5
inter_chan_window = 6
peak_assignment_err = 5 # 0,3,5
discard_first = 100                 # NOTE Remove first X samples from signal due to recording artefact
elec_stim_thresh = 500 # 180        #Threshold to find peaks
mechan_stim_thresh = 0.75

global error_log = ""

for loadinginput in eachrow(df)

    filename = loadinginput["rhs_name"]
    animalID = loadinginput["animal"]
    digit_pad = loadinginput["digit_pad"]
    pos = loadinginput["electrode_pos"]
    probes = loadinginput["electrodes"]
    stim_type = loadinginput["stim_type"]
    stim_amp = loadinginput["amplitude"]
    block_amp = loadinginput["block_amp"]
    condition = animalID[1:3]

    rhsMEA_config(probes)

    println(filename)  # tells what is in filenames

    active_chs = sort(filter(!iszero, rhsMEA_config(probes)[:]))

    max_channels = maximum(active_chs) # How many channels have recordings

    # ' ## Loading Data
    rhs_name = datafile_location * animalID * "/" * filename * ".rhs"
    read_data(rhs_name)

    # Find triggers
    if stim_type == "M3" # if using piezo stimulation, use piezo stim timing function and don't remove artefacts
        trigthresh = 0.1
        pzChannel = 1
        stim_timing = []
        pz_trig = board_adc_data[pzChannel, :]

        filteredTrigsig = filt(digitalfilter(Lowpass(5000; fs = 30000), Butterworth(4)), pz_trig)
        avtrig = mean(filteredTrigsig[10:20000])
        filteredTrigsig = filteredTrigsig .- avtrig
    
        for sampleNo in 10001:length(pz_trig)
            if filteredTrigsig[sampleNo] > trigthresh && maximum(filteredTrigsig[(sampleNo - 10000):(sampleNo - 1)]) <= trigthresh
                push!(stim_timing, sampleNo-3000)
            end
        end
    else
        stim_timing = logic_high ? board_dig_in_data[2, discard_first:end] : 1 .- board_dig_in_data[2, discard_first:end]
    end

    if stim_type == "M3" # if using piezo stimulation, use piezo stim timing function and don't remove artefacts
        if !isempty(stim_timing)
            stim_ind = stim_timing[1]
        else
            global error_log *= "$filename - no stim triggers detected ($stim_type)\n\n"
            println("$filename - no stim triggers detected (skipped)\n\n")
            continue
        end
    else  # otherwise, use the stim timing from the sync trigger
        stim_inds = []
        for (i, n) in enumerate(stim_timing)
            if n == 1 && i > 1 && stim_timing[i-1] == 0
                push!(stim_inds, i)
            end
        end
        if !isempty(stim_inds)
            stim_ind = stim_inds[1]
        else
            global error_log *= "$filename - no stim triggers detected ($stim_type)\n\n"
            println("$filename - no stim triggers detected (skipped)\n\n")
            continue
        end
    end


    # Labelling Loaded Data
    file_length = length(amplifier_data[1, discard_first:stim_ind-1])
    sig_neg = zeros(max_channels, file_length)
    sig_neg[active_chs,:] = -amplifier_data[active_chs, discard_first:stim_ind-1]
    global fs = frequency_parameters.amplifier_sample_rate
    println("Loading Finished")

    ################################################################################

    signal_matrix = filter_signal(sig_neg)
    duration = file_length / fs

    mkpath("loading_output_non-evoked-response/$animalID")
    @save "loading_output_non-evoked-response/$animalID/$filename.jld2" signal_matrix duration file_length fs stim_ind animalID condition pos stim_amp block_amp stim_type

    println("Loading Saved\n\n\n")

end

if !isempty(error_log)
    now = Dates.now()
    log_file = open("../error_logs/non-evoked.txt", "a")
    println(log_file, "---Loading $now---\n\n$error_log")
    close(log_file)
end