"""
visualize_distances.py
Visualize distances calculated by 'distances_within_between.py'. 
Requires two .npz files as input (short and long reads) and the truth ground labels.
Distances are visualized by violin plots (within cluster distances) and heatmaps (between cluster distances)

Usage
-----
    python visualize_distances.py \\
        --file_1     path/to/short_reads_cluster_distances.npz \\
        --file_2     path/to/long_reads_cluster_distances.npz \\
        --result_png path/to/cluster_comparison.png

Arguments
---------
--file_1: .npz file for short reads from distances_within_between_cluster.py
--file_2: .npz file for long reads from distances_within_between_cluster.py
--result_png: Path to save the output PNG figure.
                Default: calculated_distances.png

Output
------
A single PNG figure with two sections:
    Top    -- Violin plots: within-cluster distances per species,
              short reads (left violin) vs long reads (right violin).
    Bottom -- Heatmaps: pairwise between-cluster distances
"""

import matplotlib.pyplot as plt
import argparse
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec

# calculate number of rows and column for good pleasing grid layout
def nr_rows_cols(nr_plots):
    """
    Computes grid dimensions for the violin plot section based on number of species.
    Parameters:
        nr_plots: Number of species (= number of violin subplots needed)
    returns:
        n_rows: Total number of rows in the figure grid
    rows_violin: Number of rows used by the violin plot section.
    n_cols: Number of columns in the grid (approx. square root of nr_plots)
    """
    n_cols = int(np.ceil(np.sqrt(nr_plots)))
    rows_violin = int(np.ceil(nr_plots/n_cols))
    if n_cols < 4:
        n_rows = rows_violin +4
    else:
        n_rows = rows_violin + 2
    return n_rows,rows_violin,n_cols

# plot violin plots for within cluser distances
def violin_plot_section(fig, dataframe_list, labels, rows_violin, n_cols):
    """
    AddS violin plots for within-cluster distances to the figure (top section).

    One subplot per species; each subplot shows two violins: short reads (left)
    and long reads (right).
    
    Parameters
    fig: The figure to draw into.
    dataframe_list :
        Two DataFrames [short, long], each indexed by species label with a
        column 'dist_to_centr' containing per-read distances to centroid.
    labels: Species labels determining the order of subplots.
    rows_violin: Number of rows in the violin gridspec (from nr_rows_cols()).
    n_cols: Number of columns in the violin gridspec (from nr_rows_cols()).

    """
    gs_violin = gridspec.GridSpec(rows_violin,n_cols,left=0.05,right=0.95,
                                  top=0.95,bottom=0.4,wspace=0.4,hspace=0.2)
    # title for violin plot section
    ax_vio_tit = fig.add_axes((0.05,0.95,1,0.05))
    ax_vio_tit.set_axis_off()
    ax_vio_tit.text(
    0.5, 0.5,
    'Within cluster distances: short vs long reads',
    ha='center', va='center',
    fontsize=18, fontweight='bold'
    )
    cur_row=0
    for ind,lab in enumerate(labels):
        # dataframe_list consists of dataframe for short and long read within cluster distances
        d = [list(df.loc[lab,'dist_to_centr']) for df in dataframe_list]
        cur_col = ind%n_cols
        # add violin plots (long and short) for current cluster
        ax = fig.add_subplot(gs_violin[cur_row,cur_col])
        ax.violinplot(d)
        ax.set_title(lab)
        ax.set_xticks([1,2],labels=['short','long'])
        ax.set_ylabel('distance')
        # if current plot = last plot in this row, increase row number
        if cur_col == n_cols-1 and cur_row < rows_violin-1 and ind < len(labels)-1:
            cur_row += 1
            cur_col = 0

# plot heatmaps for between cluster distances
def heatmap_section(fig, dist_ma1, dist_ma2, n_cols, unique_labels):
    """
    Add heatmaps for between-cluster distances to the figure (bottom section
    Parameters:
        fig: The figure to draw into.
        dist_ma1: Pairwise centroid distance matrix for short reads.
        dist_ma2: Pairwise centroid distance matrix for long reads.
        n_cols: Number of columns from nr_rows_cols(); controls stacked vs side-by-side layout.
        unique_labels : Species labels for heatmap axis tick lab
    """
    all_dists = np.array([dist_ma1,dist_ma2]).flatten()
    p99 = np.percentile(all_dists,99)
    # normalize distances by 99th percentile
    norm_md1 = dist_ma1/p99
    norm_md2 = dist_ma2/p99
    # stack heatmaps vertically if n_cols too small, otherwise side-by-side
    gs_heat = gridspec.GridSpec(2 if n_cols < 4 else 1, 4,left=0.05,right=0.95,top=0.32,bottom=0.05)
    # heatmap section title
    ax_heat_tit = fig.add_axes((0,0.33,1,0.05))
    ax_heat_tit.set_axis_off()
    ax_heat_tit.text(
        0.5,0.5,'Normalized between cluster distances',
        ha='center',va='center',fontsize=18,
        fontweight='bold'
    )
    for ind, ds in enumerate([norm_md1, norm_md2]):
        ax = fig.add_subplot(gs_heat[1 if n_cols < 4 and ind == 1 else 0, ind*2:ind*2+2])
        # reverse color map (dark = long, bright = short distance)
        cmap_rev = plt.get_cmap('cividis').reversed()
        # plot distance matrix
        im = ax.imshow(ds, cmap=cmap_rev,vmin=0,vmax=1)
        # title and label for current subplot
        mode = ['Short', 'Long'][ind]
        ax.set_title(f'{mode} reads', fontsize=13, pad=10)
        ax.set_xticks(range(len(unique_labels)), labels=unique_labels, rotation=45,rotation_mode="anchor",ha="right")
        ax.set_yticks(range(len(unique_labels)), labels=unique_labels)
        ax.set_xlabel('Species')
        ax.set_ylabel('Species')
        cbar = plt.colorbar(im, ax=ax, label='Normalized Distance',shrink=0.8)
        cbar.set_label('Normalized distance',labelpad=8)

# Raise a descriptive error if expected keys are missing from an npz file.
def validate_npz(data, path):
    """
    Validate that a loaded .npz file contains all required keys
    Currently commented out in main() but can be enabled for debugging.
    """
    required_keys = {"distances", "true_label", "between"}
    missing = required_keys - set(data.files)
    if missing:
        raise KeyError(
            f"File '{path}' is missing required keys: {missing}. "
            f"Found: {set(data.files)}"
        )

def main(args):
    file1 = np.load(args.file_1)
    file2 = np.load(args.file_2)
    #validate_npz(file1, args.file_1)
    #validate_npz(file2, args.file_2)
    nr_plots = len(set(file2['true_label']))
    n_rows,rows_violin,n_cols = nr_rows_cols(nr_plots)
    unique_labels=list(set(file2['true_label']))
    fig = plt.figure(figsize=(n_cols*4, n_rows*3))  
    # plot within cluster distances
    dist1 = pd.DataFrame(file1['distances'],index=file1['true_label'],columns=['dist_to_centr'])
    dist2 = pd.DataFrame(file2['distances'],index=file2['true_label'],columns=['dist_to_centr'])
    violin_plot_section(fig, [dist1, dist2], unique_labels, rows_violin, n_cols)
    # plot between cluster distances
    dist_ma1 = file1['between']
    dist_ma2 = file2['between']
    heatmap_section(fig, dist_ma1, dist_ma2, n_cols, unique_labels)
    # save figure
    plt.savefig(args.result_png)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate distances within and between clusters')
    parser.add_argument('--file_1', type=str,default='data/short_reads/distances_short.npz', help='Path to the npz file for short reads generated by distances_within_between_cluster.py')
    parser.add_argument('--file_2', type=str, default='data/long_reads/distances_long.npz', help='Path to the npz file for long reads generated by distances_within_between_cluster.py')
    parser.add_argument('--result_png',type=str,default='calculated_distances.png',help='Path tp npz file to save calculated distances.png')
    args = parser.parse_args()
    main(args)