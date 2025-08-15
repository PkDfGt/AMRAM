import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fsize = 28
plt.rcParams.update({'font.size': fsize})
plt.rcParams['text.usetex'] = False

width = 0.35
dbsize = 20

normal_df = pd.read_csv("ResultN/access_dbsize_{}.csv".format(dbsize))
anamorphic_df = pd.read_csv("ResultN/aaccess_dbsize_{}.csv".format(dbsize))

for df in [normal_df, anamorphic_df]:
    df['logic_s'] = df['logic_ms'] / 1000.0
    df['enc_dec_s'] = df['enc_dec_ms'] / 1000.0
    df['access_s'] = df['access_ms'] / 1000.0

acsizes = normal_df['acsize'].values
x = np.arange(len(acsizes))
labels = [f"$2^{{{int(np.log2(a))}}}$" for a in acsizes]

fig, ax = plt.subplots(figsize=(12, 7))

ax.bar(
    x - width/2, normal_df['logic_s'],
    width, color='honeydew', edgecolor='green', hatch='/', label='Logic (Normal)', linewidth=2
)
ax.bar(
    x - width/2, normal_df['access_s'],
    width, bottom=normal_df['logic_s'],
    color='honeydew', edgecolor='green', hatch='\\', label='Crypto (Normal)', linewidth=2
)

ax.bar(
    x + width/2, anamorphic_df['logic_s'],
    width, color='lightcyan', edgecolor='blue', hatch='/', label='Logic (Anam)', linewidth=2
)
ax.bar(
    x + width/2, anamorphic_df['access_s'],
    width, bottom=anamorphic_df['logic_s'],
    color='lightcyan', edgecolor='blue', hatch='\\', label='Crypto (Anam)', linewidth=2
)

ax.set_xlabel("#Access")
ax.set_ylabel("Total time (s)")
ax.set_xticks(x)
ax.set_xticklabels(labels)

ax.set_yscale("log", base=2)
yticks = [2**i for i in range(6, 17, 2)]  # 2^0,2^2,2^4,2^6,2^8,2^10
ax.set_yticks(yticks)
ax.set_yticklabels([f"$2^{{{i}}}$" for i in range(6, 17, 2)])
ax.set_ylim((1 << 6) * 0.9, (1 << 16) * 1.1)
ax.grid(True, linestyle='--', alpha=0.5, axis='y')

handles, labels = ax.get_legend_handles_labels()
unique = list(dict(zip(labels, handles)).items())
# ax.legend([h for l, h in unique], [l for l, h in unique], loc='best')

plt.savefig("ResultN/Access_StackedBar_vs_ACSize_dbsize_{}.pdf".format(dbsize), format='pdf', bbox_inches='tight', dpi=300)
plt.show()
