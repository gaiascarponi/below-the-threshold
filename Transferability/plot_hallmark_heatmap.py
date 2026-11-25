# Cross-axis heatmap of the Hallmark sets (NES of the functional GSEA; see 02_functional_gsea_all_axes.R)
# Run from the project root:  python3 Transferability/plot_hallmark_heatmap.py
# Outputs: Transferability/figures/X_GSEA_Hallmark_heatmap.png/.pdf + results/hallmark_NES_matrix.csv
import csv, sys, os
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list
BASE = sys.argv[1] if len(sys.argv) > 1 else "."
OUT  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(BASE, "Transferability")
axes = [("A", "A_mass_trauma", "A\nMass trauma\n(WTC, PTSD)"), ("B", "B_war_trauma", "B\nWar trauma\n(Marines, PTSD)"),
        ("C", "C_poverty_ses", "C\nPoverty\n(HANDLS)"), ("D", "D_early_adversity", "D\nEarly adversity\n(monocytes)")]
nes, padj = {}, {}
for ax, d, _ in axes:
    for r in csv.DictReader(open(os.path.join(BASE, d, "03_results", "tables", f"{ax}02_GSEA_Hallmark.csv"))):
        k = r["pathway"].replace("HALLMARK_", "").replace("_", " ").capitalize()
        k = k.replace("Tnfa signaling via nfkb", "TNFα signaling via NF-κB").replace("Il6 jak stat3", "IL6/JAK/STAT3") \
             .replace("Il2 stat5", "IL2/STAT5").replace("Myc ", "MYC ").replace("E2f ", "E2F ").replace("G2m ", "G2M ") \
             .replace("Mtorc1", "mTORC1").replace("Pi3k akt mtor", "PI3K/AKT/mTOR").replace("Kras ", "KRAS ") \
             .replace("Tgf beta", "TGF-β").replace("Wnt beta catenin", "WNT/β-catenin").replace("Uv response", "UV response") \
             .replace("Dna repair", "DNA repair").replace("P53 ", "p53 ").replace(" dn", " down")
        nes.setdefault(k, {})[ax] = float(r["NES"]); padj.setdefault(k, {})[ax] = float(r["padj"]) if r["padj"] not in ("", "NA") else 1.0
names = sorted(nes)
M = np.array([[nes[n].get(a, np.nan) for a, _, _ in axes] for n in names])
P = np.array([[padj[n].get(a, 1.0) for a, _, _ in axes] for n in names])
# only the sets significant (padj < 0.05) in at least one axis
keep = (P < 0.05).any(axis=1); M, P = M[keep], P[keep]; names = [n for n, k in zip(names, keep) if k]
order = leaves_list(linkage(np.nan_to_num(M), "average", metric="euclidean"))
M, P, names = M[order], P[order], [names[i] for i in order]
with open(os.path.join(OUT, "results", "hallmark_NES_matrix.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["hallmark"] + [f"NES_{a}" for a, _, _ in axes] + [f"padj_{a}" for a, _, _ in axes])
    for n, m, p in zip(names, M, P): w.writerow([n] + [f"{x:.3f}" for x in m] + [f"{x:.3g}" for x in p])
cmap = LinearSegmentedColormap.from_list("div", ["#2E6DA4", "#F2F2F2", "#C0392B"])
vmax = np.nanmax(np.abs(M))
fig, ax = plt.subplots(figsize=(6.4, 0.27 * len(names) + 1.9))
im = ax.imshow(M, cmap=cmap, norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax), aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if P[i, j] < 0.05:
            ax.text(j, i, "*" if P[i, j] >= 0.001 else "**", ha="center", va="center_baseline", fontsize=9,
                    color="white" if abs(M[i, j]) > 0.62 * vmax else "#212121")
ax.set_xticks(range(len(axes))); ax.set_xticklabels([t for _, _, t in axes], fontsize=8)
ax.xaxis.tick_top()
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
ax.set_xticks(np.arange(-.5, len(axes), 1), minor=True); ax.set_yticks(np.arange(-.5, len(names), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.5); ax.tick_params(which="minor", length=0); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03); cb.set_label("NES (red: higher in exposed group)", fontsize=8)
cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
fig.text(0.01, 0.005, "* adjusted p < 0.05   ** adjusted p < 0.001   (Hallmark sets significant in at least one axis)", fontsize=7, color="#555555")
plt.tight_layout(rect=(0, 0.02, 1, 1))
os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
plt.savefig(os.path.join(OUT, "figures", "X_GSEA_Hallmark_heatmap.png"), dpi=300); plt.savefig(os.path.join(OUT, "figures", "X_GSEA_Hallmark_heatmap.pdf"))
print(len(names), "Hallmark sets plotted")
