"""
    response_quantification_v1.5.jl
    see loading for changelog

    This file contains functions:
    * hyperlink_figure2(event)
    * response_quantify(values, bins, window)
    * linearise_channels(active_chs)
    * pyplot_fidelity_matrix(responses; reorder=true)
    * pyplot_ratio_diff_grid(response_ratio_Mat, response_diff_Mat)
    * pyplot_trials(trial_data, trial_laser, trial_pks, trial_pksH; resp_ratio, resp_diff)

"""

using PyPlot
using PyCall
np = pyimport("numpy")

"""
    hyperlink_figure2(event)

    Links Figure2(Response Ratio/Diff Figure) with Channel Trials
    Important to update variables so that the link works
"""
function hyperlink_figure2(event)
    global probes, animalID
    global laser_Mat_concat_plot
    global Ch_trials_data_plot, Ch_trials_spikes_plot, Ch_trials_spikes_h_plot
    global response_ratio_Mat_plot, response_diff_Mat_plot #If only Julia is OO, then I could get these values from parent
    i = Int(round(event.ydata))
    j = Int(round(event.xdata))
    grid_ch = rhsMEA_config(probes)
    chNo = grid_ch[i,j]
    println("onclick i: ", i, " j: ", j," Channel:",chNo)
    @warn "Only one linked grid plot can be active"

    pyplot_trials(Ch_trials_data_plot[chNo], laser_Mat_concat_plot[chNo],
        Ch_trials_spikes_plot[chNo], Ch_trials_spikes_h_plot[chNo],
        E_no = chNo,
        resp_ratio = response_ratio_Mat_plot[chNo],
        resp_diff = response_diff_Mat_plot[chNo])
end


"""
    Uses values and bins generated by np.histogram

    Arguments
    `values::Array{T,1}` where values are the values of a Histogram
    `bins::Array{T,1}` the bin edge values of the Histogram
    `prewindow::Vector` prestimulus range in ms [window_start::NegativeReal, window_stop::NegativeReal]
    `postwindow::Vector` poststimulus range in ms [window_start::PositiveReal, window_stop::PositiveReal]

    Alternatively a single argument for the window can be provided:

    `window::PositiveReal` The pre- and post- stimulus range in ms
    
    Example: prewindow = [-100,0]; postwindow = [0, 100] will produce the same output as window = 100. 
    
    **Important:** window1 and window2 must be the same range or the differences and ratios may not return what is expected. A warning will be produced with the difference of window length in ms. 
    
"""
function response_quantify(values, bins, prewindow::Vector, postwindow::Vector)
    zero_stim_ind = findall(x->x == zero(eltype(bins)),bins)[1]
    bin_width = bins[zero_stim_ind+1]
    pre_bin_range1 = round(Int, prewindow[1]/bin_width)
    pre_bin_range2 = round(Int, prewindow[2]/bin_width)
    preStimArea = sum(values[zero_stim_ind-pre_bin_range1:zero_stim_ind-pre_bin_range2])

    post_bin_range1 = round(Int, postwindow[1]/bin_width)
    post_bin_range2 = round(Int, postwindow[2]/bin_width)
    postStimArea = sum( values[zero_stim_ind+post_bin_range1: zero_stim_ind+post_bin_range2])

    prew=pre_bin_range2 - pre_bin_range1
    postw = post_bin_range2 - post_bin_range1
    wind_dif = postw - prew

    if prew != postw
        @warn "Your pre- and post-windows differ by $wind_dif ms. This will affect the results"
    end 

    stim_ratio = postStimArea/preStimArea
    stim_diff = postStimArea - preStimArea
    return stim_ratio, stim_diff
end


function response_quantify(values, bins, window)
    zero_stim_ind = findall(x->x == zero(eltype(bins)),bins)[1]
    bin_width = bins[zero_stim_ind+1]
    bin_range = round(Int, window/bin_width)
    preStimArea = sum( values[zero_stim_ind-bin_range:zero_stim_ind])
    postStimArea = sum( values[zero_stim_ind: zero_stim_ind+bin_range])

    stim_ratio = postStimArea/preStimArea
    stim_diff = postStimArea - preStimArea
    return stim_ratio, stim_diff
end




"""
    function linearise_channels(active_chs)

    Reshapes a X headset probe configuration into a linear configuration with adjacent
    channels of the same shank next to each other
    `active_chs` has 2X columns
"""
function linearise_channels(active_chs)
    num_shanks = round(Int,size(active_chs,2)/2,RoundDown)
    lcc = Vector()
    for i = 1:num_shanks
        append!(lcc, vcat(active_chs[:,(2*i-1):(2*i)]'...))
    end
    return lcc
end

function pyplot_fidelity_matrix(responses; reorder=true)
    global probes, pretrial_time
    NY, NX = size(responses)
    ### Conditional Assignment True:False
    # ch_order = (reorder ? linearise_channels(rhsMEA_config(probes)) : collect(1:max_channels))
    if reorder == true
        ch_order = linearise_channels(rhsMEA_config(probes))
    else
        ch_order = collect(1:max_channels)
    end
    ### Creating Plot
    pyplot_array = PyPlot.matshow(responses[ch_order,:], extent = [0,NX, NY+0.5, 0.5], aspect="auto", cmap="coolwarm")
    PyPlot.colorbar()
    PyPlot.title("Fidelity Matrix for $animalID")
    ### Creating X and Y ticks
    ax = gca()
    # ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100*bin_per_ms))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter((x,pos)->( x/bin_per_ms - pretrial_time*1000 )))
    PyPlot.yticks(1:max_channels, ch_order)
    ### Creating X and Y Labels
    PyPlot.xlabel("Time (Bins are every $(1/bin_per_ms) ms)")
    PyPlot.ylabel("Channel Number")
    ### TESTING Recreating axes every zoom
    # fig = gcf()
    # p = ax.get_window_extent()
    # function verbose(event)
    #     println("Hello $event")
    # end
    #
    # fig.canvas.mpl_connect("draw_event",verbose)
    plt.show()
end

function pyplot_ratio_diff_grid(response_ratio_Mat, response_diff_Mat)
    global resp_ratio_thresh, resp_diff_thresh
    # grid_ch = rhsMEA_config(probes)
    grid_ch = rhsMEA_config(probes)[:, vec(mapslices(col -> any(col .!= 0), sum(rhsMEA_config(probes), dims = 1), dims = 1))] # jp
    grid_size = size(grid_ch)
    resp_ratio_data = reshape(response_ratio_Mat[grid_ch[:]], grid_size)
    resp_diff_data = reshape(response_diff_Mat[grid_ch[:]], grid_size)
    resp_ratio_data2 = resp_ratio_data .- resp_ratio_thresh
    resp_diff_data2 = resp_diff_data .- resp_diff_thresh

    B_ratio = resp_ratio_data2.<=0  # generate bit array where values neg = zero
    resp_ratio_data2 = resp_ratio_data2 .* .!B_ratio

    B_diff = resp_diff_data2.<=0
    resp_diff_data2 = resp_diff_data2 .* .!B_diff
    threshCond_Mat = .!B_ratio .* .!B_diff

    

    fig, ax = plt.subplots(1,3)
    # fig.suptitle([animalID, "Electrode event quantification for $resp_ratio_window ms window"])
    im = ax[1].matshow(threshCond_Mat, extent = [0.5, grid_size[2]+0.5, 16.5, 0.5], cmap = "afmhot")
    ax[1].set_title("ChNos above ratio + difference thresholds")
    im = ax[2].matshow(resp_ratio_data2, extent = [0.5, grid_size[2]+0.5, 16.5, 0.5], cmap = "afmhot")
    ax[2].set_title("Ratios above $resp_ratio_thresh")
    im = ax[3].matshow(resp_diff_data2, extent = [0.5, grid_size[2]+0.5, 16.5, 0.5], cmap = "afmhot")
    ax[3].set_title("Differences above $resp_diff_thresh")
    for i = 1:size(grid_ch, 1)
        for j = 1:size(grid_ch, 2)
            chNo = grid_ch[i, j]
            resp_ratio = round(resp_ratio_data[i, j], digits=1)
            resp_diff = round(resp_diff_data[i, j], digits=1)
            text = ax[1].text(j, i, chNo,
            ha="center", va="center", color="grey", fontsize = 7)
            ax[2].text(j,i, resp_ratio,
            ha="center", va="center", color="grey", fontsize = 7)
            ax[3].text(j,i, resp_diff,
            ha="center", va="center", color="grey", fontsize = 7)
        end
    end
    fig.canvas.mpl_connect("button_press_event", hyperlink_figure2)
end

"""
 pyplot_trials(trial_data, trial_laser, trial_pks, trial_pksH; resp_ratio, resp_diff)

 Plots first and last 5 trials on the left
 Plots a raster of all trials on top right
 Plots a histogram of all trials on bottom right

 Uses global animalID, pretrial_time, trial_t, fs, binwidth
 NOTE: histogram bars may not appear if binwidth is too narrow. Zoom in to check.
"""
function pyplot_trials(trial_data, trial_laser, trial_pks, trial_pksH; E_no, resp_ratio, resp_diff)
    global animalID
    global pretrial_time, trial_t, fs
    global binwidth
    tnt = length(trial_data) #total_num_trials
    FAL5 = [1,2,3,4,5,tnt-4,tnt-3,tnt-2,tnt-1,tnt] #first and last 5
    x_t = ((1:(trial_t*fs)) .* (1000/fs)) .- pretrial_time*1000

    fig = plt.figure(figsize=(16,8))

    # Create First axis first to be used as ref. ax for sharing
    gridSize = plot_all ? (2 * half_height, 2) : (1, 1)
    ax_1 = plt.subplot2grid(gridSize, (0, 0), rowspan=2, colspan=1)
    ax_1.plot(x_t, trial_data[1])
    # ax_1.plot(x_t, 1000*trial_laser[1])
    ax_1.plot((trial_pks[1] ./ (fs/1000)) .- (pretrial_time*1000), trial_pksH[1], "r.", markersize = 6, color = "red")
    if plot_all
        # Plot remaining left trials, sharing axes with ax_1
        for i=2:half_height
            ax_i = plt.subplot2grid((2*half_height,2), (2i-2, 0), rowspan=2, colspan=1, sharex=ax_1, sharey=ax_1)
            ax_i.plot(x_t, trial_data[i])
            # ax_i.plot(x_t, 1000*trial_laser[i])
            x_p = (trial_pks[i] ./ (fs/1000)) .- (pretrial_time*1000)
            ax_i.plot(x_p, trial_pksH[i],"r.", markersize = 6, color = "red")
        end
        ax_1.set_ylabel("Membrane voltage (mV)",position=(0,-4.5))
    else
        ax_1.set_ylabel("Membrane voltage (mV)")
    end
    ax_1.set_title([animalID, "Electrode No: $E_no"])
    PyPlot.xlabel("Time (ms)")

    if plot_all
        # Raster Plot
        ax_2 = plt.subplot2grid((2*half_height,2),(0,1), rowspan = half_height, colspan = 1, sharex=ax_1)
        for i = 1:tnt
            x_p = (trial_pks[i] ./ (fs/1000)) .- (pretrial_time*1000)
            ax_2.plot(x_p, i* ones(length(x_p)), "k.", markersize=6)
        end
        ax_2.set_title([animalID, "Electrode No: $E_no"])
        PyPlot.ylabel("Trial no.")
        # Histogram
        @warn "Histogram may appear to be missing bars due to narrow bins"
        ax_3 = plt.subplot2grid((2*half_height,2),(half_height,1), rowspan = half_height, colspan = 1, sharex=ax_1)
        x_ps = (vcat(trial_pks...) ./ (fs/1000)) .- (pretrial_time*1000)
        bin_edges = -pretrial_time*1000:binwidth:(trial_t-pretrial_time)*1000
        h = ax_3.hist(x_ps, bin_edges)
        hist_ratio = round(resp_ratio, digits = 3)
        hist_diff = round(resp_diff, digits = 3)
        ax_3.annotate("ratio, difference = $hist_ratio, $hist_diff", xy=(0.3,0.9), xycoords="axes fraction")
        PyPlot.xlabel("Time (ms)")
        PyPlot.ylabel("No. of spikes")
    end
end

"""
    function setup_laser_data(lasers_vec, stim_inds_vec)

    Takes in a vector of laser and prepares it into max_channels length vector of trials

"""
function setup_laser_data(lasers_vec, stim_inds_vec)
    #infer activeness from stim_inds
    global max_channels, fs, pre_trial_len, post_trial_len
    laser_trials = [Vector() for ch in 1:max_channels]
    # Should be all equal but active_ch rule means unsure of stim_inds
    for (fno, lasers) in enumerate(lasers_vec)
        stim_inds = stim_inds_vec[fno]
        first_stim_ind = minimum(minimum.(filter(x -> length(x) > 0, stim_inds)))
        bkg_start = max(1, first_stim_ind-fs)
        bkg_mean = mean(lasers[Int(bkg_start):Int(bkg_start+fs)])
        for (ch,stim_ind) in enumerate(stim_inds) #stim_ind is of a channel, stim_inds is entire board
            trials = Vector()
            for si in stim_ind # Do nothing for bad_channels
                data = lasers[(si-pre_trial_len):(si+post_trial_len-1)] .- bkg_mean
                push!(trials,data)
            end
            append!(laser_trials[ch], trials)
        end
    end
    return laser_trials
end
