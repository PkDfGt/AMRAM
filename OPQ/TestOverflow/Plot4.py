import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def plot_overflow_dbsize(prefix, fig_type, Z, folder="Result"):
    dbsizes = [2**i for i in range(15, 21)]  # N=2^15,...,2^20
    colors = ['blue', 'green', 'orange', 'red', 'purple', 'brown']
    light_colors = ['lightblue', 'lightgreen', 'moccasin', 'lightcoral', 'plum', 'lightgray']
    markers = ['o', 's', 'D', '^', 'X', 'v']

    fsize = 22
    plt.rcParams.update({'font.size': fsize})
    plt.figure(figsize=(7, 5))

    inf_value = 32  
    max_x = 256  
    data_dict = {}

    max_lambda_overall = 0
    for dbsize in dbsizes:
        dbsize_log2 = int(math.log2(dbsize))
        filename = f"group2_{prefix}_dbsize_{dbsize}_Z_{Z}_Zr_{dbsize_log2}.csv"
        path = os.path.join(folder, filename)
        if not os.path.exists(path):
            print(f"File not found: {path}")
            continue

        data = pd.read_csv(path)
        mask = data['maxBroot'] <= max_x
        thresholds = data['maxBroot'][mask].values
        probs = data['overflow_probability'][mask].values

        lambdas = np.full_like(probs, np.inf, dtype=float)
        nonzero_mask = probs > 0
        lambdas[nonzero_mask] = -np.log2(probs[nonzero_mask])

        if np.any(np.isfinite(lambdas)):
            max_lambda_overall = max(max_lambda_overall, np.nanmax(lambdas[np.isfinite(lambdas)]))

        data_dict[dbsize] = {'thresholds': thresholds, 'lambdas': lambdas}

    for dbsize, color, light_color, marker in zip(dbsizes, colors, light_colors, markers):
        if dbsize not in data_dict:
            continue
        thresholds = data_dict[dbsize]['thresholds']
        lambdas = data_dict[dbsize]['lambdas']

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
                label=f"$N=2^{{{int(math.log2(dbsize))}}}$",
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
                rest_thresholds = thresholds[first_inf_idx:]
                plt.plot(
                    rest_thresholds,
                    [inf_value] * len(rest_thresholds),
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

    all_thresholds = np.concatenate([data_dict[db]['thresholds'] for db in data_dict])
    filtered_thresholds = np.unique(all_thresholds)
    xticks = filtered_thresholds[::2]  # 隔一个
    xlabels = [f"$2^{{{int(np.log2(t))}}}$" for t in xticks]
    plt.xticks(xticks, xlabels)

    log_min = math.log2(filtered_thresholds.min())
    log_max = math.log2(max_x)
    plt.xscale('log', base=2)
    plt.xlim(2**(log_min - 0.3), 2**(log_max + 0.3))

    plt.ylim(0, inf_value + 1)
    yticks = list(range(0, 31, 10)) + [inf_value]
    yticklabels = [str(t) for t in yticks[:-1]] + [r"$\infty$"]

    ax = plt.gca()
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    plt.xlabel(r"$|B_{root}|$", fontsize=fsize)
    plt.ylabel(r"$\lambda$", fontsize=fsize)

    if prefix == "nOPQ" and Z == 3:
        plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()

    save_path = f"{fig_type}_Z{Z}_{prefix}.pdf"
    plt.savefig(save_path, format='pdf', bbox_inches='tight', dpi=300)
    print(f"Saved figure to {save_path}")

    plt.show()
    plt.close()

plot_overflow_dbsize("nOPQ", "Over2", Z=3)
plot_overflow_dbsize("aOPQ", "Over2", Z=3)
plot_overflow_dbsize("nOPQ", "Over2", Z=4)
plot_overflow_dbsize("aOPQ", "Over2", Z=4)
