import pandas as pd
import matplotlib.pyplot as plt
import math

def plot_time_comparison(xvals, labels, yvals_dict, ylabel, filename, logx=True, logy=False):
    colors = ['red', 'green', 'blue']
    light_colors = ['lightcoral', 'lightgreen', 'lightblue']
    markers = ['o', 's', '^']
    fsize = 24

    plt.rcParams.update({'font.size': fsize})
    plt.figure(figsize=(8, 5))

    for i, (label, yvals) in enumerate(yvals_dict.items()):
        plt.plot(
            xvals, yvals,
            marker=markers[i], color=colors[i],
            label=label,
            markerfacecolor=light_colors[i],
            markeredgecolor=colors[i],
            markeredgewidth=1.5,
            markersize=10,
            linestyle='-'
        )

    plt.xlabel("N" if "power2" in filename else "Block size in bytes", fontsize=fsize)
    plt.ylabel(ylabel, fontsize=fsize)

    if logx:
        plt.xscale('log', base=2)
        plt.xticks(xvals, labels)

    if logy:
        plt.yscale('log', base = 2)

    plt.grid(True, linestyle='--', alpha=0.6)
    # plt.legend()
    plt.tight_layout()
    plt.savefig(f"Result/{filename}.pdf", format='pdf', bbox_inches='tight', dpi=300)
    plt.close()


# Read data ===
heap_p2 = pd.read_csv("ResultN/HeapPQ_power2.csv")
normal_p2 = pd.read_csv("ResultN/Normal_power2.csv")
ana_p2 = pd.read_csv("ResultN/Anamorphic_power2.csv")

Ns = heap_p2['N'].tolist()
labels_N = [f"$2^{{{int(math.log2(n))}}}$" for n in Ns]

# Plot
plot_time_comparison(
    Ns, labels_N,
    {
        "Non-oblivious": heap_p2['avg_insert_ms'] * 1000,
        "Normal": normal_p2['avg_insert_ms'] * 1000,
        "Anamorphic": ana_p2['avg_insert_ms'] * 1000
    },
    ylabel="Time ($\mu$s)",
    filename="Extract_Insert_power2",
)

# Plot
plot_time_comparison(
    Ns, labels_N,
    {
        "Non-oblivious": heap_p2['avg_extract_ms'] * 1000,
        "Normal": normal_p2['avg_extract_ms'] * 1000,
        "Anamorphic": ana_p2['avg_extract_ms'] * 1000
    },
    ylabel="Time ($\mu$s)",
    filename="Extract_Extract_power2"
)

# === 读取 bsize 数据 ===
heap_bs = pd.read_csv("Result/N20_HeapPQ_bsize.csv")
normal_bs = pd.read_csv("Result/N20_Normal_bsize.csv")
ana_bs = pd.read_csv("Result/N20_Anamorphic_bsize.csv")

bsizes = heap_bs['bsize'].tolist()
labels_bs = [f"$2^{{{int(math.log2(bsize))}}}$" for bsize in bsizes]
# labels_bs = [str(b) for b in bsizes]

# bsize 插入时间图
plot_time_comparison(
    bsizes, labels_bs,
    {
        "Non-oblivious": heap_bs['avg_insert_ms'] * 1000,
        "Normal": normal_bs['avg_insert_ms'] * 1000,
        "Anamorphic": ana_bs['avg_insert_ms'] * 1000
    },
    ylabel="Time ($\mu$s)",
    filename="Extract_Insert_bsize",
    logy =True
)

# bsize 提取时间图
plot_time_comparison(
    bsizes, labels_bs,
    {
        "Non-oblivious": heap_bs['avg_extract_ms'] * 1000,
        "Normal": normal_bs['avg_extract_ms'] * 1000,
        "Anamorphic": ana_bs['avg_extract_ms'] * 1000
    },
    ylabel="Time ($\mu$s)",
    filename="Extract_Extract_bsize",
    logy =True
)
