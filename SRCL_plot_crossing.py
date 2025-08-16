import numpy as np
from scipy.interpolate import LinearNDInterpolator
from matplotlib import pyplot as plt
import pickle as pkl
import sys
import os
import re

os.chdir(os.getcwd() + r"\SRCL_data")

#####################################
# Get the probability distributions #
#####################################

# Construct a probability density in the space of continuous SRCL and E for each pair of indices L and W.

def load_SRCL_data():
    # First we need to check the SRCL_data folder for SRCL data files
    # and take note of which L and W we have data for.
    matched_filenames = []
    data_dict = {}
    max_E_dict = {}
    for filename in os.listdir():
        file_re = re.compile("SRCL_data_L(.*)_W(.*)_c(.*)_(.*).pkl")
        file_re_result = re.match(file_re, filename)
        if file_re_result:
            matched_filenames.append(filename)
            [file_L, file_W] = [int(file_re_result.group(1)), float(file_re_result.group(2))]
            # Initialize empty lists for each pair of L and W for which data is present.
            # Do not enter the data at this point because the list will be reset each loop.
            data_dict[file_L, file_W] = {}
            data_dict[file_L, file_W]["E"] = np.array([])
            data_dict[file_L, file_W]["SRCL"] = np.array([])
            max_E_dict[file_W] = 0
    # Now load the data using the stored filenames.
    # This way we don't need to create any intermediary data holders.
    # We'll get L and W again from the same regular expression.
    for filename in matched_filenames:
        file_re_result = re.match(file_re, filename)
        [file_L, file_W] = [int(file_re_result.group(1)), float(file_re_result.group(2))]
        with open(filename, "rb") as file_in:
            [all_E, all_SRCL] = pkl.load(file_in)
        # If there are multiple data files with the same W and L, their data will be combined data_dict.
        data_dict[file_L, file_W]["E"] = np.concatenate((data_dict[file_L, file_W]["E"], all_E))
        data_dict[file_L, file_W]["SRCL"] = np.concatenate((data_dict[file_L, file_W]["SRCL"], all_SRCL))
        # Find the maximum energy accountirng for symmetry) for each W, comparing over L
        if np.max(np.abs(all_E)) > max_E_dict[file_W]:
            max_E_dict[file_W] = np.max(np.abs(all_E))
    return [data_dict, max_E_dict]

def make_2d_hist(data_dict, nbins_E, nbins_SRCL):
    # Bin the data points by pairs of values (SRCL, E) for each (L, W).
    # This effectively gives a probability density function of the SRCL, E, L, and W.
    SRCL_2d_counts_dict = {}
    for LW_key in data_dict.keys():
        SRCL_2d_counts_dict[LW_key] = {}
        (L,W) = LW_key
        Emax = 6 + W # Bound for the simple cubic lattice rom the Perron–Frobenius theorem
        SRCLmax = 0.55*L # From calculating the SRCL of a plane wave in a periodic system, which is just below 0.5*L - this gives a margin of error for to differences due to finite size or the discrete lattice
        [SRCL_count, E_bin_edges, SRCL_bin_edges] = np.histogram2d(data_dict[LW_key]["E"], data_dict[LW_key]["SRCL"], bins=[nbins_E,nbins_SRCL], range=[[-Emax,Emax],[0,SRCLmax]])
        # histogram2d returns bin edges; take averages of the left and right edges to get the bin centers.
        SRCL_2d_counts_dict[LW_key]["counts"] = SRCL_count
        SRCL_2d_counts_dict[LW_key]["E_ax"] = (E_bin_edges[:-1]+E_bin_edges[1:])/2
        SRCL_2d_counts_dict[LW_key]["SRCL_ax"] = (SRCL_bin_edges[:-1]+SRCL_bin_edges[1:])/2
    return SRCL_2d_counts_dict

def plot_SRCL_prob_dens(LW_key):
    # Use this to inspect the resulting probability density surface visually.
    [EE, SS] = np.meshgrid(SRCL_2d_counts_dict[LW_key]["E_ax"], SRCL_2d_counts_dict[LW_key]["SRCL_ax"], indexing="ij")
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(EE, SS, SRCL_2d_counts_dict[LW_key]["counts"])
    ax.set_xlabel("Energy")
    ax.set_ylabel("SRCL")
    ax.set_zlabel("Counts")
    plt.show()

########################################
# Get the average SRCL for L, W, and E #
########################################

# In principle, the average SRCL is a continuous function of all three of L, W, and E.
# We will interpolate over W and E to get continuous functions for them. We won't need continuous L.
# This will give us surfaces in the space of W and E indexed by L.

def make_average_SRCL_hist(data_dict, nbins_E, nbins_SRCL):
    # Calculate the average SRCL for each triplet (L, W, E).
    SRCL_avg_dict = {}
    for LW_key in data_dict.keys():
        SRCL_avg_dict[LW_key] = {}
        (L,W) = LW_key
        Emax = 6 + W # Bound for the simple cubic lattice from the Perron–Frobenius theorem
        [SRCL_sum, E_bin_edges] = np.histogram(data_dict[LW_key]["E"], bins=nbins_E, range=[-Emax,Emax], weights=data_dict[LW_key]["SRCL"])
        SRCL_count = np.histogram(data_dict[LW_key]["E"], bins=nbins_E, range=[-Emax,Emax])[0]
        SRCL_avg_dict[LW_key]["SRCL_avg"] = SRCL_sum/SRCL_count
        # There will likely be empty bins on the edges with a 0/0 division error, resulting in nan. For now I am ignoring these; I can fill them in with 0 later.
        SRCL_avg_dict[LW_key]["E_ax"] = (E_bin_edges[:-1]+E_bin_edges[1:])/2
    # SRCL_avg_dict has the average SRCL in that bin.
    # Alternatively, it might also be OK to use the midpoint of that bin's edges on the SRCL axis, but this is more precise.
    return SRCL_avg_dict

def plot_SRCL_avg_vs_E(LW_key):
    # Use this to plot to check the average SRCL as a function of E for fixed W
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(SRCL_avg_dict[LW_key]["E_ax"], SRCL_avg_dict[LW_key]["SRCL_avg"])
    ax.set_xlabel("Energy")
    ax.set_ylabel("Average SRCL")
    plt.show()

def make_average_SRCL_interp(SRCL_avg_dict):
    # Now interpolate between W to create a continuous function <SRCL>(W, E) indexed by L
    # To set up for the interpolation function, we need this to be indexed by L only, not by the pair L and W,
    # and to put the data for all W in one array.
    SRCL_avg_interp_data = {}
    for LW_key in SRCL_avg_dict.keys():
        (L,W) = LW_key
        SRCL_avg_interp_data[L] = {}
        SRCL_avg_interp_data[L]["E"] = np.array([])
        SRCL_avg_interp_data[L]["W"] = np.array([])
        SRCL_avg_interp_data[L]["SRCL_avg"] = np.array([])
    # Now enter in the data
    for LW_key in SRCL_avg_dict.keys():
        (L,W) = LW_key
        SRCL_avg_interp_data[L]["E"] = np.concatenate((SRCL_avg_interp_data[L]["E"], SRCL_avg_dict[LW_key]["E_ax"]))
        SRCL_avg_interp_data[L]["W"] = np.concatenate((SRCL_avg_interp_data[L]["W"], W*np.ones(len(SRCL_avg_dict[LW_key]["E_ax"]))))
        SRCL_avg_interp_data[L]["SRCL_avg"] = np.concatenate((SRCL_avg_interp_data[L]["SRCL_avg"], SRCL_avg_dict[LW_key]["SRCL_avg"]))
    # Now construct the interpolations
    SRCL_avg_interps = {}
    for L in SRCL_avg_interp_data.keys():
        SRCL_avg_interps[L] = LinearNDInterpolator(list(zip(SRCL_avg_interp_data[L]["E"], SRCL_avg_interp_data[L]["W"])), SRCL_avg_interp_data[L]["SRCL_avg"])
    return SRCL_avg_interps

def plot_SRCL_avg_vs_W_and_E(L):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    # Plot resolution
    n_plot = 100
    E_lim = np.max([E for E in max_E_dict.values()])
    W_lim = np.max([W for W in max_E_dict.keys()])
    E_plot_ax = np.linspace(-E_lim, E_lim, n_plot)
    W_plot_ax = np.linspace(0, W_lim, n_plot)
    [EE, WW] = np.meshgrid(E_plot_ax, W_plot_ax)
    ax.plot_surface(EE, WW, SRCL_avg_interps[L](EE, WW))
    ax.set_xlabel("Energy")
    ax.set_ylabel("W")
    ax.set_zlabel("Average SRCL")
    plt.show()

###################################################
# Get the cuts across of average SRCL vs W        #
# for a certain fractional X position in the band #
# E = E_min + E_fraction*(E_max-E_min)            #
###################################################

def make_E_max_interp(max_E_dict):
    # Interpolate the maximum energies of the band to estimate the band width for any W
    Ws_sorted = np.sort([W for W in max_E_dict.keys()])
    max_Es_sorted = np.array([E for E in max_E_dict.values()])[np.argsort([W for W in max_E_dict.keys()])]
    E_max_interp = lambda W_arg: np.interp(W_arg, Ws_sorted, max_Es_sorted)
    return E_max_interp

def make_linecut(E_fraction, E_max_interp, n_EW_linecut):
    # A curve in the E-W plane that is (E_fraction)% of the way from -E_max to +E_max
    Ws_sorted = np.sort([W for W in max_E_dict.keys()])
    W_EW_linecut = np.linspace(Ws_sorted[0], Ws_sorted[-1], n_EW_linecut)
    E_EW_linecut = (2*E_fraction-1)*E_max_interp(W_EW_linecut)
    # SRCL_avg along the curve for each L
    SRCL_avg_linecut = {}
    scaled_SRCL_avg_linecut = {}
    for L in SRCL_avg_interps.keys():
        SRCL_avg_linecut[L] = SRCL_avg_interps[L](E_EW_linecut, W_EW_linecut)
        scaled_SRCL_avg_linecut[L] = SRCL_avg_linecut[L]/L
    return [W_EW_linecut, E_EW_linecut, SRCL_avg_linecut, scaled_SRCL_avg_linecut]

def plot_crossings(E_fraction, n_EW_linecut):
    # Plot SRCL/L vs disorder
    fig = plt.figure()
    ax = plt.axes()
    [W_EW_linecut, E_EW_linecut, SRCL_avg_linecut, scaled_SRCL_avg_linecut] = make_linecut(E_fraction, E_max_interp, n_EW_linecut)
    for L in SRCL_avg_interps.keys():
        plt.plot(W_EW_linecut, scaled_SRCL_avg_linecut[L])
    plt.show()

def plot_crossing_error(E_fraction, SRCL_avg_interps, E_max_interp, n_EW_linecut):
    [W_EW_linecut, E_EW_linecut, SRCL_avg_linecut, scaled_SRCL_avg_linecut] = make_linecut(E_fraction, E_max_interp, n_EW_linecut)
    crossing_width = 0
    for L in SRCL_avg_interps.keys():
        for Lp in SRCL_avg_interps.keys():
            crossing_width += (scaled_SRCL_avg_linecut[L] - scaled_SRCL_avg_linecut[Lp])**2
    crossing_width = ( 1/(len(SRCL_avg_interps.keys())*(len(SRCL_avg_interps.keys())-1)) )*crossing_width
    plt.plot(W_EW_linecut, crossing_width)
    plt.show()




nbins_E = 20
nbins_SRCL = 20

[data_dict, max_E_dict] = load_SRCL_data()
SRCL_avg_dict = make_average_SRCL_hist(data_dict, nbins_E, nbins_SRCL)
SRCL_avg_interps = make_average_SRCL_interp(SRCL_avg_dict)
E_max_interp = make_E_max_interp(max_E_dict)

plot_crossings(E_fraction=0.5, n_EW_linecut=100)
plot_crossing_error(0.5, SRCL_avg_interps, E_max_interp, n_EW_linecut=100)
