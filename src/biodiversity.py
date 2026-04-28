import numpy as np
import matplotlib.pyplot as plt


def biodiversity_curve(otu_table_df, depth_step=5000, max_depth=None, figsize=(7, 5)):
    """
    Generate and plot rarefaction curves for microbial diversity analysis.

    Parameters:
    -----------
    otu_table_df : pd.DataFrame
        OTU table with Feature ID as column and samples as rows
    depth_step : int, optional
        Step size for sequencing depths (default: 5000)
    max_depth : int, optional
        Maximum sequencing depth (default: None, uses max from data)
    figsize : tuple, optional
        Figure size as (width, height) in inches

    Returns:
    --------
    dict : rarefaction_curves dictionary with sample IDs as keys
    """

    # Process OTU table
    otu_table = otu_table_df.set_index("Feature ID").T
    otu_table = otu_table.round().astype(int)

    def rarefaction_richness(counts, depth):
        if counts.sum() < depth:
            return np.nan

        expanded = np.repeat(np.arange(len(counts)), counts)
        subsample = np.random.choice(expanded, depth, replace=False)
        return len(np.unique(subsample))

    # Set max_depth if not provided
    if max_depth is None:
        max_depth = otu_table.sum(axis=1).max()

    depths = np.arange(100, max_depth, depth_step)

    rarefaction_curves = {}

    for sample in otu_table.index:
        counts = otu_table.loc[sample].values
        rarefaction_curves[sample] = [rarefaction_richness(counts, d) for d in depths]

    # Plot
    plt.figure(figsize=figsize)

    for sample, richness in rarefaction_curves.items():
        plt.plot(depths, richness, label=sample)

    plt.xlabel("Sequencing depth")
    plt.ylabel("Observed OTUs")
    plt.title("Alpha rarefaction curves")
    plt.tight_layout()
    plt.show()

    return rarefaction_curves
