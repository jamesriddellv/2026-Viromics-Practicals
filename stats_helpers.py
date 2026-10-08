"""
stats_helpers.py
----------------
Helper functions and classes for ecological statistics, custom distance matrices,
and visualization utilities without requiring skbio.
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
import seaborn as sns
from sklearn.decomposition import PCA

def get_annotations(votu, df):
  """
  Retrieves annotations for a specific vOTU from a DataFrame.

  Args:
    votu (str): The name of the vOTU to search for.
    df (pd.DataFrame): The DataFrame containing vOTU annotations, expected
                       to have a 'vOTU' column.

  Returns:
    pd.DataFrame: A DataFrame containing all rows where the 'vOTU' column
                  matches the specified vOTU name.
  """
  return df.loc[df['vOTU'] == votu]

# ── 1. Distance Matrix & ANOSIM Helpers ─────────────────────────────────────

class DistanceMatrix:
    """Lightweight replacement for skbio.stats.distance.DistanceMatrix."""

    def __init__(self, data, ids=None):
        data = np.asarray(data, dtype=float)

        # Handle 1D condensed distance matrix from scipy pdist
        if data.ndim == 1:
            self.data = squareform(data)
        elif data.ndim == 2 and data.shape[0] == data.shape[1]:
            self.data = data
        else:
            raise ValueError(
                "Input data must be a square matrix or a condensed 1D distance array."
            )

        self.shape = self.data.shape
        self.ids = (
            list(ids) if ids is not None else [str(i) for i in range(self.shape[0])]
        )
        self._id_to_index = {id_: i for i, id_ in enumerate(self.ids)}

    def __getitem__(self, index):
        if isinstance(index, tuple) and len(index) == 2:
            i, j = index
            i_idx = self._id_to_index.get(i, i) if isinstance(i, str) else i
            j_idx = self._id_to_index.get(j, j) if isinstance(j, str) else j
            return self.data[i_idx, j_idx]
        return self.data[index]

    def redundant(self):
        """Returns the full 2D square matrix (matches skbio behavior)."""
        return self.data

    def condensed(self):
        """Returns 1D condensed array compatible with scipy."""
        return squareform(self.data, checks=False)


def anosim(distance_matrix, grouping, permutations=999, seed=None):
    """Pure NumPy/SciPy replacement for skbio.stats.distance.anosim.

    Parameters
    ----------
    distance_matrix : DistanceMatrix object, 2D numpy array, or 1D condensed vector
    grouping : list-like or pd.Series with group labels for each sample
    permutations : int, number of Monte Carlo permutations for p-value calculation
    seed : int (optional), random seed for reproducibility

    Returns
    -------
    pd.Series containing:
        - method: 'ANOSIM'
        - test statistic: R
        - p-value: p
        - number of samples: N
        - number of permutations: permutations
    """
    if seed is not None:
        np.random.seed(seed)

    # Extract 2D array if DistanceMatrix instance is passed
    dm = (
        distance_matrix.data
        if hasattr(distance_matrix, "data")
        else np.asarray(distance_matrix)
    )
    if dm.ndim == 1:
        dm = squareform(dm)

    n_samples = dm.shape[0]
    grouping = np.asarray(grouping)

    if len(grouping) != n_samples:
        raise ValueError(
            "Grouping length must match the number of samples in distance matrix."
        )

    # Convert lower triangle of distance matrix into rank vector
    tril_indices = np.tril_indices(n_samples, k=-1)
    distances = dm[tril_indices]
    ranks = rankdata(distances)

    # Masks for within-group vs between-group pairwise comparisons
    group_i = grouping[tril_indices[0]]
    group_j = grouping[tril_indices[1]]
    within_mask = group_i == group_j
    between_mask = ~within_mask

    # R-statistic calculation helper
    divisor = (n_samples * (n_samples - 1)) / 4.0

    def calculate_r(within_m, between_m, ranks_arr):
        mean_r_w = np.mean(ranks_arr[within_m])
        mean_r_b = np.mean(ranks_arr[between_m])
        return (mean_r_b - mean_r_w) / divisor

    stat_r = calculate_r(within_mask, between_mask, ranks)

    # Permutation test
    perm_r_stats = np.zeros(permutations)
    for p in range(permutations):
        perm_grouping = np.random.permutation(grouping)
        perm_group_i = perm_grouping[tril_indices[0]]
        perm_group_j = perm_grouping[tril_indices[1]]

        p_within = perm_group_i == perm_group_j
        p_between = ~p_within

        perm_r_stats[p] = calculate_r(p_within, p_between, ranks)

    # Empirical p-value (inclusive of observed statistic)
    p_value = (np.sum(perm_r_stats >= stat_r) + 1.0) / (permutations + 1.0)

    return pd.Series(
        {
            "method": "ANOSIM",
            "test statistic": stat_r,
            "p-value": p_value,
            "number of samples": n_samples,
            "number of permutations": permutations,
        }
    )


# ── 2. Biodiversity Index Helpers ───────────────────────────────────────────

def hill(counts, order=1):
    """Computes the Hill diversity number (effective number of species) for any order q.

    Parameters
    ----------
    counts : array-like
        Abundance or frequency counts of species/OTUs/ASVs.
    order : float or int, default=1
        Diversity order q:
          - q = 0: Species richness
          - q = 1: Exponential Shannon index
          - q = 2: Inverse Simpson index

    Returns
    -------
    float
        Hill number (effective number of species).
    """
    q = order
    counts = np.asarray(counts, dtype=float)

    # Filter out zeros to avoid log(0) and 0^negative power errors
    counts = counts[counts > 0]

    if len(counts) == 0:
        return 0.0

    # Convert counts to relative abundances (proportions p_i)
    p = counts / np.sum(counts)

    # --- Case 1: q = 0 (Richness) ---
    if q == 0:
        return float(len(p))

    # --- Case 2: q -> 1 (Exponential Shannon Diversity) ---
    elif np.isclose(q, 1.0):
        shannon_entropy = -np.sum(p * np.log(p))
        return float(np.exp(shannon_entropy))

    # --- Case 3: General case for any q != 1 (e.g., q = 2) ---
    else:
        return float(np.sum(p**q) ** (1.0 / (1.0 - q)))


# ── 3. Plotting & Visualization Helpers ────────────────────────────────────

def add_confidence_ellipse(x, y, ax, n_std=2.447, **kwargs):
    """Adds a 95% confidence ellipse to a 2D scatterplot (replicates ggplot2's stat_ellipse)."""
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size < 3:
        return  # Need at least 3 points for an ellipse

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    ellipse = Ellipse(
        (0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor="none",
        **kwargs
    )

    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std

    transf = (
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(scale_x, scale_y)
        .translate(np.mean(x), np.mean(y))
    )
    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)


def add_stripplot(ax, plot_df, y_col, marker_map):
    """Draws a stripplot layer onto an existing Matplotlib axis based on timepoints.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to draw the plot on.
    plot_df : pd.DataFrame
        DataFrame containing 'timepoint_str', 'treatment', and y_col.
    y_col : str
        Column name for the y-axis values.
    marker_map : dict
        Mapping of timepoint_str values to matplotlib markers (e.g. {'0': 'o', '7': 's'}).
    """
    for tp, m in marker_map.items():
        subset = plot_df[plot_df["timepoint_str"] == tp]
        if subset.empty:
            continue
        if tp == "7":
            sns.stripplot(
                data=subset,
                x="treatment",
                y=y_col,
                ax=ax,
                marker="s",
                color="black",
                alpha=0.6,
                jitter=False,
                s=9,
            )
            sns.stripplot(
                data=subset,
                x="treatment",
                y=y_col,
                ax=ax,
                marker="x",
                color="white",
                alpha=1.0,
                jitter=False,
                s=6,
                linewidth=1,
            )
        else:
            sns.stripplot(
                data=subset,
                x="treatment",
                y=y_col,
                ax=ax,
                marker=m,
                color="black",
                alpha=0.6,
                jitter=True,
                s=8,
                linewidth=0.5,
            )

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd


def plot_genes(
    df_genes: pd.DataFrame,
    label_col: str = "kegg_hit",
    figsize: tuple = (12, 3.5),
    buffer_bp: int = 500,
    y_lim_bottom: float = -20.0,
    show_plot: bool = True,
):
    """Plots genomic gene tracks as arrow patches with vertical annotation labels.

    Parameters:
    -----------
    df_genes : pd.DataFrame
        DataFrame containing gene metadata with required columns:
        'start_position', 'end_position', 'strandedness', and the specified
        label_col.
    label_col : str, optional
        Column name to use for labeling genes (e.g., 'kegg_hit', 'viral_hit',
        'peptidase_hit', 'cazy_best_hit', 'vogdb_hits',
        'annotation_description'). Default is 'kegg_hit'.
    figsize : tuple, optional
        Dimensions of the matplotlib figure (width, height). Default is (12,
        3.5).
    buffer_bp : int, optional
        Extra sequence padding added to the maximum end position for plot
        limits. Default is 500.
    y_lim_bottom : float, optional
        Lower y-axis limit to accommodate rotated vertical labels. Default is
        -20.0.
    show_plot : bool, optional
        Whether to call plt.show() immediately. Default is True.

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    genome_length = df_genes["end_position"].max() + buffer_bp

    fig, ax = plt.subplots(figsize=figsize)

    for _, row in df_genes.iterrows():
        start, end, strand = (
            row["start_position"],
            row["end_position"],
            row["strandedness"],
        )
        length = end - start
        color = "#4C72B0" if strand in ["+", 1] else "#DD8452"
        head_len = min(length * 0.25, 200)

        dx = length if strand in ["+", 1] else -length
        x_start = start if strand in ["+", 1] else end

        # 1. Draw Arrow
        arrow = patches.FancyArrow(
            x_start,
            0,
            dx,
            0,
            width=0.3,
            head_width=0.5,
            head_length=head_len,
            length_includes_head=True,
            color=color,
            ec="black",
            linewidth=0.5,
        )
        ax.add_patch(arrow)

        # 2. Add Label Below (if column exists and value isn't empty/NaN)
        if label_col in row and pd.notna(row[label_col]):
            label_text = str(row[label_col])
            midpoint = start + (length / 2)  # Center text relative to the gene

            ax.text(
                x=midpoint,
                y=-0.35,  # Placed just below the arrow track
                s=label_text,
                rotation=90,
                ha="center",  # Horizontal center aligned with gene midpoint
                va="top",  # Top of text anchored at y=-0.35 (grows downwards)
                fontsize=6,
                clip_on=True,  # Keeps labels within plot bounds
            )

    # Baseline genome track
    ax.plot([0, genome_length], [0, 0], color="black", linewidth=1, zorder=0)

    # Styling
    ax.set_xlim(0, genome_length)
    ax.set_ylim(y_lim_bottom, 1)
    ax.set_yticks([])
    ax.set_xlabel("Genomic Position (bp)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    plt.tight_layout()

    if show_plot:
        plt.show()

    return fig, ax
