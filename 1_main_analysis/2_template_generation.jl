version = v"4.0.4"
# Applies principal component analysis and k-means clustering to template files in order to produce
# unit assignment templates for spikesorting. Template information, including cluster centres and
# cluster sizes, are saved into /spikesorting_templates in the .JLD2 format.

# * Changelog
#   Updated to TSHelperTS4.4.0
#   Updated def_template_gen_vars() name

# ! Errors
#   Sometimes ISIresult is missing a value

#   TODO List

# ? Suggestions & Queries



# Packages
using DSP
using Glob
using LinearAlgebra
using MultivariateStats
using Printf
using PyPlot; pygui(true)

using Statistics # std
using StatsBase # mode

using CSV
using JLD2
using FileIO
using PyCall
using DataFrames
using XLSX
using Query
using Dates

pe = pyimport("matplotlib.patheffects")
CheckButtons = matplotlib.widgets.CheckButtons

# Working directory and dependencies
cd(dirname(@__FILE__))
include("../src/spike_sort_main.jl")
include("../src/spike_sort_plotting.jl")
include("../src/ts_helper.jl")       #! Must remain after other dependencies

df_t = DataFrame(XLSX.readtable("saved_files_info.xlsx", "Descriptions"))

global error_log = ""

for loadingoutput in eachrow(df_t)
    templateFilename = loadingoutput["Template"]
    animalID = loadingoutput["Animal"]
    stim_type = loadingoutput["Stimtype"]
    println("running next file: ", animalID, " ", templateFilename)

    Random.seed!(42) #set the random seed for kmeans clustering, to allow for consistency

    # # * Parameters
    k_override = 0

    # Spikesorting Parameters
    fs=30000

    if stim_type == "L1"
        pre_trial_len = round(Int, 1.0 * fs)
        post_trial_len = round(Int,5.0 * fs)
        WOI = [800, 1500]
    elseif stim_type == "M3"
        pre_trial_len = round(Int, 0.5 * fs)
        post_trial_len = round(Int,1.5 * fs)
        WOI = [0, 500]
    elseif stim_type == "M2"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI = [0, 400]
    elseif stim_type == "M1"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI = [0, 40]
    else
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI = [0, 100]
    end

    onems = round(Int, fs/1000)
    pre_window = 7  # peak extraction
    post_window = 7 # peak extraction

    # Function Parameters
    samples_on_x=false  # what units to use for x-axis of raster plots - true for samples or false for milliseconds
    showgraphs = false # toggle graph generation for faster debugging
    global ISIresult = [] #set up a global to store the ISI results


    ################################################################################



    if isfile("loading_output/$animalID/$templateFilename.jld2")
        templateFile = load("loading_output/$animalID/$templateFilename.jld2")
    else
        global error_log *= "$filename - loading output missing\n\n"
        println("$filename - loading output missing (skipped)\n\n")
        continue
    end
    signal_matrix, pk_bkgnds, file_lengths, stim_inds, cum_stim_inds = def_template_gen_vars(templateFile)

    # * template construction
    #' ## Peak Detection
    pk_inds, waveforms, _pk_heights = detect_spikes(
        signal_matrix,
        pk_bkgnds,
        stim_inds,
        [pre_window, post_window],
        file_lengths = file_lengths,
        window_of_interest = WOI .* onems
    )
                        
    #' ## Peak Clustering and Template Generation
    cluster_result, templates, PCA_Model, assign_order = cluster_spikes(
        waveforms,
        k_override = k_override,
        apply_pca=true,
        PCA_ratio = 0.99,
        PCA_scale = "off",
        max_clusters = 10,
        num_epocs = 100,
        SHOW_CLUSTER_RESULT = showgraphs
    )
    assigns = cluster_result.assignments
    println("Initial counts: ", counts(assigns), " (k-override = $k_override)")



    # * template application
    #' Get thresholds for clustering
    # use cluster_pratio = 1.0, then reduce until you remove the clearly unwanted waveforms. Conduct the following lines of code iteratively until satisfied
    thresholds = get_cluster_thresholds(waveforms, templates, assigns, cluster_pratio = 1.0, x_min = 10)
    # thresholds = get_cluster_thresholds(waveforms, templates, assigns, cluster_pratio = 0.95, x_min = 10)

    #' define data file to apply templates
    pk_inds_full, waveforms_full, _pk_heights_full = detect_spikes(signal_matrix, pk_bkgnds, stim_inds, [pre_window, post_window], file_lengths = file_lengths)
        
    #' Assigning spikes to templates
    assignment = assign_points_to_clusters(waveforms_full, templates, thresholds)

    #Plots of assigned waveforms
    if showgraphs
        pyplot_waveforms_separately(waveforms_full, assignment, centroids = templates, title = "Assigned from Templates")
    end

    # Combine all stim_inds
    pks_plot, ass_plot = assign_to_trials(cum_stim_inds, pk_inds_full, assignment, pre_trial_len, post_trial_len) 
    # pks_plot, ass_plot = assign_to_trials(cum_stim_inds, pk_inds2, assignment, 0, 300) 

    # pks for plotting, assignment for plotting
    #FIXME:     I've added pre/post trial lengths, but maybe these should go above somewhere, 
    #           e.g. the assignments above?
    #TODO:      Check I did this correctly: # cum_stim_inds = stim_inds2# Should replace for concat_stim_inds
    #           check we need the pre_trial_len and post_ are required

    # separate_events(pks_plot, ass_plot)
    ass_counts = counts(convert.(Int64, vcat(ass_plot...)))
    num_assigns = maximum(vcat(ass_plot...)) + 1
    numTrials = size(ass_plot, 1)
    pksLocs = Array{Any}(undef, maximum(vcat(ass_plot...))+1)
    pksAss = Array{Any}(undef, maximum(vcat(ass_plot...))+1)
    for i = 1:num_assigns
        pksLocs[i]=[]
        pksAss[i]=[]
        for j = 1:numTrials
            push!(pksLocs[i],pks_plot[j][findall(ass_plot[j] .== (i-1))])
            push!(pksAss[i],ass_plot[j][findall(ass_plot[j] .== (i-1))])
        end
        
        # ISI testing
        testSum = []
        for trial in pksLocs[i]
            push!(testSum, !isempty(trial))
        end
        if sum(testSum) > 0
            println("Assignment $i")
            result = isiTests(pksLocs[i], onems)
            println(result, "\n")
            push!(ISIresult, result)
        else println("Assignment $i is empty\n")
            push!(ISIresult, "EMPTY")
        end
    end

    global k_override = 0 # ! this will override every looped retemplate, so should be deprecated
    global templateinit = zeros(length(unique(assigns))) #make a counter to store how many retemplates have occured
    global retemplate = 1
    global templatecounter = []
    global templatecounternew = []

    while true

        if templatecounter == [] 
            templatecounter = templateinit
        else 
            templatecounter = templatecounternew
        end

        if retemplate > length(templatecounter)
            println("All retemplating completed!")
            break
        end

        if templatecounter[retemplate] > 3 || ISIresult[retemplate + 1] != "FAILED"
            retemplate += 1
            println("Templating limit reached for this assignment, moving to next assignment")
            continue

        else
            showgraph = false

            #retemplate -= 1     #! DO NOT REMOVE - this transforms retemplate into the assignment number before unassigned is added
            waveformsRetemplate = fndAssign(waveforms, retemplate, assigns)

            # Retemplate single assignment
            cluster_resultNew, templatesNew, PCA_ModelNew, assign_orderNew = cluster_spikes(waveformsRetemplate,
                        k_override=k_override,
                        apply_pca=true,
                        PCA_ratio=0.99,
                        PCA_scale="off",
                        max_clusters=10,
                        num_epocs=100,
                        SHOW_CLUSTER_RESULT=showgraph)
            assignsNew = cluster_resultNew.assignments
            println("Retemplated counts: ", counts(assignsNew), " (k-override = $k_override)")

            assignsOriginal = assigns
            templatesOriginal = templates
            global ISIresult = []

            assigns, templates, templatecounternew = combine_templates(retemplate, assignsOriginal, templatesOriginal, assignsNew, templatesNew, templatecounter, showgraph)
            thresholds = get_cluster_thresholds(waveforms, templates, assigns, cluster_pratio = 1.0, x_min = 10)
            pk_inds_full, waveforms_full, _pk_heights_full = detect_spikes(signal_matrix, pk_bkgnds, stim_inds, [pre_window, post_window], file_lengths = file_lengths)
            assignment_temp = assign_points_to_clusters(waveforms, templates, thresholds)
            assignment = assign_points_to_clusters(waveforms_full, templates, thresholds)
            if showgraph
                pyplot_waveforms_separately(waveforms_full, assignment, centroids = templates, title = "Assigned from Templates")
            end
            pks_plot, ass_plot = assign_to_trials(cum_stim_inds, pk_inds_full, assignment, pre_trial_len, post_trial_len)
            ass_counts = counts(convert.(Int64, vcat(ass_plot...)))
            num_assigns = maximum(assigns) + 1
            numTrials = size(ass_plot, 1)
            pksLocs = Array{Any}(undef, num_assigns)
            pksAss = Array{Any}(undef, num_assigns)
            for i = 1:num_assigns
                pksLocs[i]=[]
                pksAss[i]=[]
                for j = 1:numTrials
                    push!(pksLocs[i],pks_plot[j][findall(ass_plot[j] .== (i-1))])
                    push!(pksAss[i],ass_plot[j][findall(ass_plot[j] .== (i-1))])
                end
                testSum = []
                for trial in pksLocs[i]
                    push!(testSum, !isempty(trial))
                end
                if sum(testSum) > 0
                    println("Assignment $i")
                    result = isiTests(pksLocs[i], onems)
                    println(result, "\n")
                    push!(ISIresult, result)
                else println("Assignment $i is empty\n")
                    push!(ISIresult, "EMPTY")
                end
            end
            println("Retemplating completed, checking next assignment...")
        end
    end

    #* Save output
    mkpath("../spikesorting_templates/$animalID")
    @save "../spikesorting_templates/$animalID/$templateFilename.jld2" waveforms templates assigns thresholds ISIresult
    println("Template generation complete, moving to next file..")
end
println("Batch templating complete!")

if !isempty(error_log)
    now = Dates.now()
    log_file = open("../error_logs/main_analysis.txt", "a")
    println(log_file, "---Template Generation $now---\n\n$error_log")
    close(log_file)
end