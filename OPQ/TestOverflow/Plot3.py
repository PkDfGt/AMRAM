import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def plot_overflow_lambda(prefix, fig_type, folder="Result"):
    Z_values = [3, 4, 5, 6]
    colors = ['blue', 'green', 'orange', 'red']
    light_colors = ['lightblue', 'lightgreen', 'moccasin', 'lightcoral']
    markers = ['o', 's', 'D', '^']

    fsize = 22
    plt.rcParams.update({'font.size': fsize})
    plt.figure(figsize=(7, 5))

    ymax_value = 30
    inf_value = 32
    max_x = 256

    data_dict = {}

    for Z in Z_values:
        filename = f"group1_{prefix}_dbsize_32768_Z_{Z}_Zr_15.csv"
        path = os.path.join(folder, filename)
        data = pd.read_csv(path)

        mask = data['maxBroot'].values <= max_x
        thresholds = data['maxBroot'].values[mask]
        probs = data['overflow_probability'].values[mask]

        lambdas = np.full_like(probs, np.inf, dtype=float)
        nonzero_mask = probs > 0
        lambdas[nonzero_mask] = -np.log2(probs[nonzero_mask])

        data_dict[Z] = {
            'thresholds': thresholds,
            'probs': probs,
            'lambdas': lambdas
        }

    for Z, color, light_color, marker in zip(Z_values, colors, light_colors, markers):
        thresholds = data_dict[Z]['thresholds']
        lambdas = data_dict[Z]['lambdas']

        lambdas_plot = np.where(np.isinf(lambdas), inf_value, lambdas)

        finite_idx = ~np.isinf(lambdas)
        infinite_idx = np.isinf(lambdas)

        if np.any(finite_idx):
            last_finite_idx = np.where(finite_idx)[0][-1]

            plt.plot(
                thresholds[finite_idx],
                lambdas_plot[finite_idx],
                marker=marker,
                linestyle='-',
                label=f"$Z={Z}$",
                color=color,
                markerfacecolor=light_color,
                markeredgecolor=color,
                markeredgewidth=1.5,
                markersize=10,
                zorder=2
            )

            plt.scatter(
                thresholds[last_finite_idx],
                lambdas_plot[last_finite_idx],
                marker=marker,
                color=color,
                facecolors=light_color,
                edgecolors=color,
                linewidths=1.5,
                s=100,
                zorder=5,
                label=None
            )

            if np.any(infinite_idx):
                first_inf_idx = np.where(infinite_idx)[0][0]

                plt.plot(
                    [thresholds[last_finite_idx], thresholds[first_inf_idx]],
                    [lambdas_plot[last_finite_idx], inf_value],
                    linestyle='--',
                    color=color,
                    linewidth=2,
                    zorder=3
                )

                zero_thresholds = thresholds[first_inf_idx:]
                plt.plot(
                    zero_thresholds,
                    [inf_value] * len(zero_thresholds),
                    linestyle='--',
                    color=color,
                    linewidth=2,
                    zorder=3
                )
        else:
            plt.plot(
                thresholds,
                [inf_value] * len(thresholds),
                linestyle='--',
                color=color,
                linewidth=2,
                zorder=3
            )

    plt.xscale('log', base=2)

    all_thresholds = np.concatenate([data_dict[Z]['thresholds'] for Z in Z_values])
    filtered_thresholds = np.unique(all_thresholds)

    xticks = filtered_thresholds[::2]
    xlabels = [f"$2^{{{int(np.log2(t))}}}$" for t in xticks]
    plt.xticks(xticks, xlabels)

    log_min = math.log2(filtered_thresholds.min())
    log_max = math.log2(max_x)
    plt.xlim(2**(log_min - 0.3), 2**(log_max + 0.3))

    ymin = 0
    ymax = inf_value + 1
    y_margin = (ymax - ymin) * 0.1
    plt.ylim(ymin - y_margin, ymax)

    plt.xlabel(r"$|B_{root}|$", fontsize=fsize)
    plt.ylabel(r"$\lambda$", fontsize=fsize)

    yticks = list(range(0, 31, 10)) + [inf_value]
    yticklabels = [str(t) for t in yticks[:-1]] + [r"$\infty$"]

    ax = plt.gca()
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()

    save_path = f"{fig_type}_{prefix}_lambda_final.pdf"
    plt.savefig(save_path, format='pdf', bbox_inches='tight', dpi=300)
    print(f"Saved figure to {save_path}")

    plt.show()
    plt.close()

plot_overflow_lambda("aOPQ", "Over1")
plot_overflow_lambda("nOPQ", "Over1")
