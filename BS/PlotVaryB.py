import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

fsize = 28
plt.rcParams.update({'font.size': fsize})

normal_df = pd.read_csv("Result/normal_B_N20_vary.csv")
normal_df["time_s"] = normal_df["time_ms"] / 1000
anam_df = pd.read_csv("Result/anamorphic_B_N20_vary.csv")
anam_df["time_s"] = anam_df["time_ms"] / 1000

B_min = 2**4
B_max = 2**14
normal_df = normal_df[(normal_df["BSize"] >= B_min) & (normal_df["BSize"] <= B_max)]
anam_df = anam_df[(anam_df["BSize"] >= B_min) & (anam_df["BSize"] <= B_max)]

Ns = sorted(normal_df["BSize"].unique())
threads = sorted(normal_df["threads"].unique())

normal_cmap = LinearSegmentedColormap.from_list(
    "light_greens", cm.Greens(np.linspace(0.2, 0.6, 256))
)
anam_cmap = LinearSegmentedColormap.from_list(
    "light_blues", cm.Blues(np.linspace(0.2, 0.6, 256))
)

width = 0.208
x = np.arange(len(Ns))
fig, ax = plt.subplots(figsize=(10, 6))

def plot_stacked_bars_with_diff_cmap(ax, df, x_positions, offset, cmap, threads):
    n = len(threads)
    for xi, N in enumerate(Ns):
        sub = df[df["BSize"] == N][["threads", "time_s"]].copy()
        sub_sorted = sub.sort_values("time_s", ascending=True).reset_index(drop=True)

        bottoms = 0
        for i in range(len(sub_sorted)):
            t = int(sub_sorted.loc[i, "threads"])
            if i == 0:
                height = sub_sorted.loc[i, "time_s"]
            else:
                height = sub_sorted.loc[i, "time_s"] - sub_sorted.loc[i - 1, "time_s"]
            color_idx = threads.index(t)
            color = cmap(color_idx / (n - 1))
            ax.bar(x_positions[xi] + offset, height, width=width, bottom=bottoms,
                   color=color, edgecolor='black')
            bottoms += height

plot_stacked_bars_with_diff_cmap(ax, normal_df, x, -width/2, normal_cmap, threads)
plot_stacked_bars_with_diff_cmap(ax, anam_df, x, width/2, anam_cmap, threads)

ax.set_xticks(x)
ax.set_xticklabels([f"$2^{{{int(np.log2(n))}}}$" for n in Ns])
ax.set_xlabel("Block size in bytes")
ax.set_ylabel("Time (s)")
ax.grid(True, linestyle='--', alpha=0.5, axis='y')
ax.set_yscale("log", base = 2)
yticks = [2**i for i in range(2, 15, 4)]  # 2^0,2^2,2^4,2^6,2^8,2^10
ax.set_yticks(yticks)
ax.set_yticklabels([f"$2^{{{i}}}$" for i in range(2, 15, 4)])
ax.set_ylim((1 << 2) * 0.9, (1 << 14) * 1.1)

norm = mcolors.Normalize(vmin=0, vmax=len(threads)-1)

plt.savefig("Result/BS_VaryB_N20_Bar.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.show()
