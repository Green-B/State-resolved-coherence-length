import numpy as np
import pickle as pkl
import sys
import os
import re

nbins_E = 20
nbins_SRCL = 20

os.chdir(os.getcwd() + r"\SRCL_data")

# First we need to check the SRCL_data folder for SRCL data files
# and take note of which L and W we have data for.
matched_filenames = []
data_dict = {}
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

# Now load the data using the stored filenames.
# This way we don't need to create any intermediary data holders.
# We'll get L and W again from the same regular expression.
for filename in matched_filenames:
    file_re_result = re.match(file_re, filename)
    [file_L, file_W] = [int(file_re_result.group(1)), float(file_re_result.group(2))]
    with open(filename, "rb") as file_in:
        [all_E, all_SRCL] = pkl.load(file_in)
    data_dict[file_L, file_W]["E"] = np.concatenate((data_dict[file_L, file_W]["E"], all_E))
    data_dict[file_L, file_W]["SRCL"] = np.concatenate((data_dict[file_L, file_W]["SRCL"], all_SRCL))

# Bin the data points by pairs of values (SRCL, E) for each (L, W).
# This effectively gives a probability density function of the SRCL, E, L, and W.
SRCL_2d_counts_dict = {}
for LW_key in data_dict.keys():
    SRCL_2d_counts_dict[LW_key] = {}
    (L,W) = LW_key
    Emax = 6 + W # Bound for the simple cubic lattice rom the Perron–Frobenius theorem
    SRCLmax = 0.55*L # From calculating the SRCL of a plane wave in a periodic system, which is just below 0.5*L
    [SRCL_count, E_bin_edges, SRCL_bin_edges] = np.histogram2d(data_dict[LW_key]["E"], data_dict[LW_key]["SRCL"], bins=[nbins_E,nbins_SRCL], range=[[-Emax,Emax],[0,SRCLmax]])
    # histogram2d returns bin edges; take averages of the left and right edges to get the bin centers.
    SRCL_2d_counts_dict[LW_key]["counts"] = SRCL_count
    SRCL_2d_counts_dict[LW_key]["E_ax"] = (E_bin_edges[:-1]+E_bin_edges[1:])/2
    SRCL_2d_counts_dict[LW_key]["SRCL_ax"] = (SRCL_bin_edges[:-1]+SRCL_bin_edges[1:])/2

"""
from matplotlib import pyplot as plt

LW_key = (6, 8.0)
[EE, SS] = np.meshgrid(SRCL_2d_counts_dict[LW_key]["E_ax"], SRCL_2d_counts_dict[LW_key]["SRCL_ax"], indexing="ij")
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(EE, SS, SRCL_2d_counts_dict[LW_key]["counts"])
ax.set_xlabel("Energy")
ax.set_ylabel("SRCL")
ax.set_zlabel("Counts")
plt.show()
"""

# Calculate the average SRCL for each triplet (L, W, E).
SRCL_avg_dict = {}
for LW_key in data_dict.keys():
    SRCL_avg_dict[LW_key] = {}
    (L,W) = LW_key
    Emax = 6 + W # Bound for the simple cubic lattice from the Perron–Frobenius theorem
    [SRCL_sum, E_bin_edges] = np.histogram(data_dict[LW_key]["E"], bins=nbins_E, range=[-Emax,Emax], weights=data_dict[LW_key]["SRCL"])
    SRCL_count = np.histogram(data_dict[LW_key]["E"], bins=nbins_E, range=[-Emax,Emax])[0]
    SRCL_avg_dict[LW_key]["SRCL_avg"] = SRCL_sum/SRCL_count
    SRCL_avg_dict[LW_key]["E_ax"] = (E_bin_edges[:-1]+E_bin_edges[1:])/2
    SRCL_avg_dict[LW_key]["SRCL_ax"] = (SRCL_bin_edges[:-1]+SRCL_bin_edges[1:])/2
# SRCL_avg_dict has the average SRCL in that bin.
# Alternatively, it might also be OK to use the midpoint of that bin's edges on the SRCL axis.


