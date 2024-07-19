import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


hex_colors_red = ["#E16278", "#F7A072"]
custom_cmap_red = mcolors.LinearSegmentedColormap.from_list('HarmonicReds', hex_colors_red)
hex_colors_blue = ["#25719E", "#B3E1F8"]
custom_cmap_blue = mcolors.LinearSegmentedColormap.from_list('HarmonicBlues', hex_colors_blue)
seg_reds = custom_cmap_red(np.linspace(0, 1, 50))
seg_blues = custom_cmap_blue.reversed()(np.linspace(0, 1, 50))
combined_colors = np.vstack((seg_reds, seg_blues))
combined_cmap = mcolors.ListedColormap(combined_colors)

def plot_tsne_spearman(df):
    vcenter = 0
    vmin, vmax = df['diff_spearmanr'].min(), df['diff_spearmanr'].max()
    norm = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
    norm = mcolors.CenteredNorm(vcenter=0)

    fig, ax = plt.subplots(figsize=(20, 14))
    ax = sns.scatterplot(data=df,
                        y='x_1', x='x_0', size=2, linewidth=0, hue_norm=norm,
                        hue='diff_spearmanr', legend=False, palette=combined_cmap.reversed(), 
                        ax=ax)
    ax.set_xlabel('t-SNE 1', size=20)
    ax.set_ylabel('t-SNE 2', size=20)
    ax.tick_params(axis='both', labelsize=20)

    scalarmappaple = cm.ScalarMappable(norm=norm, cmap=combined_cmap.reversed())
    scalarmappaple.set_array(df['diff_spearmanr'])
    cbar = fig.colorbar(scalarmappaple, shrink=0.5, location='left', ax=ax,
                        label='diff Spearman correlation (two-stage - single-stage)')
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_ylabel('Difference in Spearman correlation (two-stage - single-stage)', size=15)
    
    return fig


df = pd.read_csv('data_tsne_spearman_rmse.csv', index_col=0)
fig = plot_tsne_spearman(df)
fig.savefig('fig_cluster_split_spearman.pdf')
