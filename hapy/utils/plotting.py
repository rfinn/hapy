import numpy as np
import matplotlib.pyplot as plt


def raincloud_by_group(
    values_by_group,
    group_labels,
    xlabel,
    title=None,
    figsize=(10, 6),
    jitter=0.05,
    alpha_points=0.15,
    alpha_violin=0.35,
    alpha_box=0.35,
):
    """
    Make a horizontal raincloud plot.

    Parameters
    ----------
    values_by_group : list of 1D arrays
        One array of values per group.
    group_labels : list of str
        Labels for each group.
    xlabel : str
        Label for x-axis.
    title : str, optional
        Plot title.
    """

    fig, ax = plt.subplots(figsize=figsize)

    # drop non-finite values
    clean_values = []
    clean_labels = []
    for vals, lab in zip(values_by_group, group_labels):
        vals = np.asarray(vals)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        clean_values.append(vals)
        clean_labels.append(lab)

    if len(clean_values) == 0:
        print(f"WARNING: no valid values for raincloud plot: {title}")
        plt.close(fig)
        return None, None

    # boxplot
    bp = ax.boxplot(
        clean_values,
        patch_artist=True,
        vert=False,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        widths=0.12,
    )

    for patch in bp["boxes"]:
        patch.set_alpha(alpha_box)

    # violin plot
    vp = ax.violinplot(
        clean_values,
        points=300,
        showmeans=False,
        showextrema=False,
        showmedians=False,
        vert=False,
    )

    for i, body in enumerate(vp["bodies"]):
        yvals = body.get_paths()[0].vertices[:, 1]
        body.get_paths()[0].vertices[:, 1] = np.clip(yvals, i + 1, i + 1.55)
        body.set_alpha(alpha_violin)

    # jittered points
    rng = np.random.default_rng(12345)
    for i, vals in enumerate(clean_values):
        y = np.full(len(vals), i + 0.78)
        y += rng.uniform(low=-jitter, high=jitter, size=len(y))
        ax.scatter(vals, y, s=10, alpha=alpha_points)

    ax.set_yticks(np.arange(1, len(clean_labels) + 1))
    ax.set_yticklabels(clean_labels)
    ax.set_xlabel(xlabel)

    if title is not None:
        ax.set_title(title)

    plt.tight_layout()
    return fig, ax
