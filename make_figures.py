import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
os.makedirs("figures", exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 0.9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.dpi": 200, "savefig.bbox": "tight",
})
NAVY, BLUE, RED, GREEN, GREY, ORANGE = "#1F4E79", "#2E75B6", "#C00000", "#2E7D32", "#777777", "#E69F00"
FLOOR, CEIL = 0.267, 1.377

# ---------------- Figure 1: ML information ceiling (SR Pareto front) ----------------
comp = np.array([1, 3, 5, 7, 8, 9])
loss = np.array([3.285586, 1.675298, 1.591946, 1.500809, 1.428088, 1.405166])
rmse = np.sqrt(loss)
fig, ax = plt.subplots(figsize=(6.2, 4.3))
ax.axhspan(0, FLOOR, color=GREEN, alpha=0.10)
ax.axhline(FLOOR, color=GREEN, ls="--", lw=1.3, label=f"natural floor = {FLOOR:.3f} MeV")
ax.axhline(CEIL, color=RED, ls="--", lw=1.3, label=f"ML ceiling (CV) = {CEIL:.2f} MeV")
ax.plot(comp, rmse, "o-", color=NAVY, lw=2, ms=7, label="SR Pareto front (in-sample)")
ax.annotate("V \u2248 diagonal \u00d7 monopole-centroid", xy=(3, rmse[1]), xytext=(4.0, 0.80),
            fontsize=10, color=NAVY,
            arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.1))
ax.set_xlabel("equation complexity"); ax.set_ylabel("RMSE to USDB (MeV)")
ax.set_title("Information ceiling of quantum-number features", fontweight="bold", color=NAVY)
ax.set_ylim(0, 2.0); ax.set_xticks(comp); ax.legend(frameon=False, fontsize=9, loc="upper right")
fig.savefig("figures/fig1_ml_ceiling.png"); plt.close(fig)

# ---------------- Figure 2: zero-parameter baseline scan ----------------
names = ["Kuo\n(G-matrix)", "Brown-Kuo", "Wildenthal\nUSD", "N3LO\n(ab-initio)", "JISP\n(ab-initio)"]
vals = [1.081, 1.144, 1.020, 0.484, 0.466]
cols = [GREY, GREY, GREY, BLUE, NAVY]
fig, ax = plt.subplots(figsize=(6.6, 4.3))
bars = ax.bar(range(len(vals)), vals, color=cols, width=0.62, zorder=3)
ax.axhline(CEIL, color=RED, ls="--", lw=1.3, label=f"ML ceiling = {CEIL:.2f} MeV")
ax.axhline(FLOOR, color=GREEN, ls="--", lw=1.3, label=f"natural floor = {FLOOR:.3f} MeV")
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width()/2, v + 0.03, f"{v:.3f}", ha="center", fontsize=9.5)
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=9.5)
ax.set_ylabel("RMSE to USDB (MeV), zero fit parameters")
ax.set_title("Modern ab-initio is an excellent parameter-free baseline", fontweight="bold", color=NAVY)
ax.set_ylim(0, 1.55); ax.legend(frameon=False, fontsize=9.5, loc="upper right", bbox_to_anchor=(0.99, 0.80))
fig.savefig("figures/fig2_baseline_scan.png"); plt.close(fig)

# ---------------- Figure 3: monopole correction, sd vs fp ----------------
sd_T0 = [-0.368, -0.276, -0.602, -0.567, -0.302, -0.407]
sd_T1 = [-0.032, 0.493, 0.387, 0.164, 0.517, 0.491]
fp_T0_mean, fp_T1_mean = -0.326, -0.059
sd_T0_m, sd_T1_m = np.mean(sd_T0), np.mean(sd_T1)
fig, ax = plt.subplots(figsize=(6.4, 4.3))
rng = np.random.default_rng(0)
ax.scatter(rng.normal(0, 0.05, len(sd_T0)), sd_T0, color=BLUE, alpha=0.65, s=42, zorder=3, label="sd per orbit-pair (JISP)")
ax.scatter(rng.normal(1, 0.05, len(sd_T1)), sd_T1, color=BLUE, alpha=0.65, s=42, zorder=3)
ax.plot([-0.22, 0.22], [sd_T0_m]*2, color=NAVY, lw=3, zorder=4)
ax.plot([0.78, 1.22], [sd_T1_m]*2, color=NAVY, lw=3, zorder=4, label="sd mean")
ax.plot([-0.22, 0.22], [fp_T0_mean]*2, color=ORANGE, lw=3, ls=(0, (4, 2)), zorder=4)
ax.plot([0.78, 1.22], [fp_T1_mean]*2, color=ORANGE, lw=3, ls=(0, (4, 2)), zorder=4, label="fp mean (kb3g)")
ax.axhline(0, color="k", lw=0.7)
ax.annotate("transfers \u2713", xy=(0, fp_T0_mean), xytext=(0.30, -0.52), color=GREEN, fontsize=10, fontweight="bold")
ax.annotate("does not \u2717", xy=(1, fp_T1_mean), xytext=(1.18, 0.18), color=RED, fontsize=10, fontweight="bold")
ax.set_xticks([0, 1]); ax.set_xticklabels(["T = 0  (isoscalar)", "T = 1  (isovector)"])
ax.set_ylabel("monopole correction  USDB \u2212 ab-initio (MeV)")
ax.set_title("Monopole correction: isoscalar part transfers sd \u2192 fp", fontweight="bold", color=NAVY)
ax.set_xlim(-0.45, 1.65); ax.legend(frameon=False, fontsize=9, loc="lower right")
fig.savefig("figures/fig3_monopole_transfer.png"); plt.close(fig)

# ---------------- Figure 4: 18O and 18F level schemes ----------------
def scheme(ax, data, title, ymax):
    # data: dict col -> list of (label, E, color)
    cols = list(data.keys()); xw = 0.6
    for i, col in enumerate(cols):
        for (lab, E, c) in data[col]:
            ax.hlines(E, i - xw/2, i + xw/2, color=c, lw=2.2)
            ax.text(i + xw/2 + 0.04, E, lab, va="center", fontsize=8.5, color=c)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=9.5)
    ax.set_ylim(-0.2, ymax); ax.set_ylabel("excitation energy (MeV)")
    ax.set_title(title, fontweight="bold", color=NAVY); ax.set_xlim(-0.5, len(cols) - 0.25)

T0C, T1C = NAVY, RED
o18 = {
    "exp":  [("0\u207a", 0.0, T1C), ("2\u207a", 1.98, T1C), ("4\u207a", 3.55, T1C)],
    "USDB": [("0\u207a", 0.0, T1C), ("2\u207a", 1.998, T1C), ("4\u207a", 3.527, T1C)],
    "JISP": [("0\u207a", 0.0, T1C), ("2\u207a", 1.732, T1C), ("4\u207a", 3.265, T1C)],
}
f18 = {
    "exp":  [("1\u207a", 0.0, T0C), ("3\u207a", 0.94, T0C), ("0\u207a(T1)", 1.04, T1C), ("5\u207a", 1.12, T0C)],
    "USDB": [("1\u207a", 0.0, T0C), ("3\u207a", 0.945, T0C), ("5\u207a", 1.241, T0C), ("0\u207a(T1)", 1.481, T1C)],
    "JISP": [("1\u207a", 0.0, T0C), ("3\u207a", 0.987, T0C), ("0\u207a(T1)", 1.609, T1C), ("5\u207a", 1.661, T0C)],
}
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
scheme(axes[0], o18, "$^{18}$O   (all T=1)", 4.1)
scheme(axes[1], f18, "$^{18}$F   (T=0 blue, T=1 red)", 2.0)
fig.suptitle("ab-initio reproduces A=18 spectra without fitting", fontweight="bold", color=NAVY, y=1.00)
fig.savefig("figures/fig4_spectra.png"); plt.close(fig)

print("figures written:")
import os
for f in sorted(os.listdir("figures")):
    print("  figures/" + f)
