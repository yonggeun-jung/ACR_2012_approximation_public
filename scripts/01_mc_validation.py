"""Monte Carlo: gravity diagnostic, tables, elasticity figure."""

import csv
import os
import sys
import warnings

import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from distributions import FrechetPareto, LogNormal, ParetoMixture, TruncatedPareto

SIGMA = 4.0
KAPPA = 5.0
N_COUNTRIES = 30
TAU_LOW, TAU_HIGH = 1.1, 3.0
TAU_OLD = 2.0
LIBS = [(1.8, "10%"), (1.5, "25%"), (1.2, "40%")]
SEED = 42

DGPS = [
    ("Frechet/Pareto", FrechetPareto(KAPPA, phi_min=1.0)),
    ("Log-normal", LogNormal(mu=0.0, sigma_ln=1.0)),
    ("Mixture", ParetoMixture(kappa1=5.0, kappa2=3.5, pi=0.9)),
    ("Truncated", TruncatedPareto(kappa=5.0, phi_min=1.0, phi_max=50)),
]

OUTDIR_TAB = os.path.join(os.path.dirname(__file__), "..", "output", "tables")
OUTDIR_FIG = os.path.join(os.path.dirname(__file__), "..", "output", "figures")
os.makedirs(OUTDIR_TAB, exist_ok=True)
os.makedirs(OUTDIR_FIG, exist_ok=True)


class Economy:
    def __init__(self, dist, N, sigma, f_x=1.5, f_d=1.0):
        self.dist = dist
        self.N = N
        self.sigma = sigma
        self.f_x = f_x
        self.f_d = f_d

    def cutoff(self, tau):
        s = self.sigma - 1.0
        pm = getattr(self.dist, "phi_min", 1.0)
        return pm * (self.f_x / self.f_d) ** (1.0 / s) * tau

    def ln_trade(self, tau):
        """Melitz-style: ln Φ(φ*)."""
        ps = self.cutoff(tau)
        Phi = self.dist.Phi(ps, self.sigma)
        return np.log(Phi) if Phi > 0 else -np.inf

    def ln_trade_ek(self, tau):
        """EK-style: ln(1 - F(φ*)) = ln SF(φ*)."""
        ps = self.cutoff(tau)
        s = self.dist.sf(ps)
        return np.log(s) if s > 0 else -np.inf

    def true_eps(self, tau):
        return self.dist.eps_melitz(self.cutoff(tau), self.sigma)

    def true_eps_ek(self, tau):
        return self.dist.eps_ek(self.cutoff(tau))

    def welfare_error(self, tau_old, tau_new, eps_bench, model="melitz"):
        N = self.N
        Phi_d = self.dist.Phi(self.cutoff(1.0), self.sigma)
        n_grid = 500
        t_grid = np.linspace(0, 1, n_grid)
        ln_to, ln_tn = np.log(tau_old), np.log(tau_new)

        ln_lam = np.zeros(n_grid)
        eps_path = np.zeros(n_grid)
        for k, t in enumerate(t_grid):
            tau_t = np.exp((1 - t) * ln_to + t * ln_tn)
            Phi_x = self.dist.Phi(self.cutoff(tau_t), self.sigma)
            if model == "ek":
                eps_path[k] = self.true_eps_ek(tau_t)
            else:
                eps_path[k] = self.true_eps(tau_t)
            ln_lam[k] = np.log(Phi_d / (Phi_d + (N - 1) * Phi_x))

        ln_lam_hat = ln_lam[-1] - ln_lam[0]
        d_ln_lam = np.gradient(ln_lam, t_grid)
        ln_W_exact = np.trapz(-d_ln_lam / eps_path, t_grid)
        ln_W_acr = -(1.0 / eps_bench) * ln_lam_hat
        return ln_W_exact - ln_W_acr, ln_lam_hat


def polynomial_diagnostic(econ, tau_vec, trade_flow="melitz"):
    """Quadratic in ln τ fitted to ln(trade). Use Melitz Φ or EK SF consistently."""
    ln_tau = np.log(tau_vec)
    if trade_flow == "ek":
        ln_X = np.array([econ.ln_trade_ek(t) for t in tau_vec])
    else:
        ln_X = np.array([econ.ln_trade(t) for t in tau_vec])
    valid = np.isfinite(ln_X)
    ln_tau = ln_tau[valid]
    ln_X = ln_X[valid]

    A_lin = np.column_stack([np.ones_like(ln_tau), ln_tau])
    b_lin = np.linalg.lstsq(A_lin, ln_X, rcond=None)[0]
    eps_pool = -b_lin[1]

    ln_tau2 = ln_tau ** 2
    A_quad = np.column_stack([np.ones_like(ln_tau), ln_tau, ln_tau2])
    b_quad = np.linalg.lstsq(A_quad, ln_X, rcond=None)[0]
    beta1, beta2 = b_quad[1], b_quad[2]

    ln_tau_min, ln_tau_max = ln_tau.min(), ln_tau.max()
    eps_at_min = -(beta1 + 2 * beta2 * ln_tau_min)
    eps_at_max = -(beta1 + 2 * beta2 * ln_tau_max)
    eps_lo = min(eps_at_min, eps_at_max)
    eps_hi = max(eps_at_min, eps_at_max)
    Delta_eps = eps_hi - eps_lo

    return {
        "eps_pool": eps_pool,
        "beta2": beta2,
        "eps_min": eps_lo,
        "eps_max": eps_hi,
        "Delta_eps": Delta_eps,
    }


def _setup_pub_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif", "Times", "serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "lines.linewidth": 1.25,
            "lines.solid_capstyle": "round",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
            "axes.unicode_minus": False,
        }
    )


def save_elasticity_figure():
    _setup_pub_style()

    colors = {
        "Frechet/Pareto": "#000000",
        "Log-normal": "#CC0000",
        "Mixture": "#0044AA",
        "Truncated": "#555555",
    }
    line_specs = {
        "Frechet/Pareto": {"linestyle": "-", "marker": "o"},
        "Log-normal": {"linestyle": "--", "marker": "s"},
        "Mixture": {"linestyle": ":", "marker": "^"},
        "Truncated": {"linestyle": "-.", "marker": "D"},
    }

    ln_tau_grid = np.linspace(np.log(TAU_LOW), np.log(TAU_HIGH), 200)
    tau_grid = np.exp(ln_tau_grid)
    markevery = max(1, len(ln_tau_grid) // 12)

    series_ek, series_mel = [], []
    ymin_ek, ymax_ek = np.inf, 0.0
    ymin_mel, ymax_mel = np.inf, 0.0
    for dgp_name, dist in DGPS:
        econ = Economy(dist, N_COUNTRIES, SIGMA)
        ek = np.array([econ.true_eps_ek(t) for t in tau_grid])
        mel = np.array([econ.true_eps(t) for t in tau_grid])
        series_ek.append((dgp_name, ek))
        series_mel.append((dgp_name, mel))
        ymin_ek = min(ymin_ek, float(np.nanmin(ek)))
        ymax_ek = max(ymax_ek, float(np.nanmax(ek)))
        ymin_mel = min(ymin_mel, float(np.nanmin(mel)))
        ymax_mel = max(ymax_mel, float(np.nanmax(mel)))

    def _y_limits(ymin, ymax):
        y_floor = 1.5 if ymin >= 1.5 else ymin * 0.98
        y_top = ymax * 1.06
        return y_floor, y_top

    def _plot_single_fig(series, y_floor, y_top, outfile, legend_dy_axes=0.0):
        fig, ax = plt.subplots(figsize=(5.5, 5.5 / 1.618), constrained_layout=False)
        fig.subplots_adjust(left=0.13, right=0.97, bottom=0.15, top=0.94)
        for dgp_name, eps_vals in series:
            spec = line_specs[dgp_name]
            ax.plot(
                ln_tau_grid,
                eps_vals,
                label=dgp_name,
                color=colors[dgp_name],
                linestyle=spec["linestyle"],
                marker=spec["marker"],
                markevery=markevery,
                markersize=4.0,
                markerfacecolor=colors[dgp_name],
                markeredgecolor=colors[dgp_name],
                markeredgewidth=0.6,
                lw=2.25,
                clip_on=False,
            )
        ax.set_ylim(y_floor, y_top)
        ax.set_xlabel("Trade Costs (log)")
        ax.set_ylabel("Trade Elasticity")
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4, prune=None))
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune=None))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        leg_y = 0.12 / (y_top - y_floor) + legend_dy_axes
        ax.legend(
            loc="lower right",
            bbox_to_anchor=(0.99, leg_y),
            ncol=2,
            frameon=False,
            handlelength=3.0,
            columnspacing=0.9,
            borderaxespad=0.35,
        )
        fig.savefig(outfile, format="pdf", dpi=300)
        plt.close(fig)

    yf_e, yt_e = _y_limits(ymin_ek, ymax_ek)
    yf_m, yt_m = _y_limits(ymin_mel, ymax_mel)
    out_ek = os.path.join(OUTDIR_FIG, "fig_elasticity_profile_ek.pdf")
    out_mel = os.path.join(OUTDIR_FIG, "fig_elasticity_profile_melitz.pdf")
    _plot_single_fig(series_ek, yf_e, yt_e, out_ek, legend_dy_axes=0.15)
    _plot_single_fig(series_mel, yf_m, yt_m, out_mel, legend_dy_axes=0.0)
    plt.rcdefaults()


def main():
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(SEED)
    tau_vec = np.array(
        [rng.uniform(TAU_LOW, TAU_HIGH) for _ in range(N_COUNTRIES * (N_COUNTRIES - 1))]
    )

    tab1 = []
    hdr = (
        f"{'DGP':25s} {'eps_p_m':>8s} {'b2_m':>8s} {'De_m':>8s} "
        f"{'eps_p_e':>8s} {'b2_e':>8s} {'De_e':>8s}"
    )
    print(hdr)
    print("-" * len(hdr))

    for dgp_name, dist in DGPS:
        econ = Economy(dist, N_COUNTRIES, SIGMA)
        diag_m = polynomial_diagnostic(econ, tau_vec, trade_flow="melitz")
        diag_e = polynomial_diagnostic(econ, tau_vec, trade_flow="ek")
        tab1.append(
            {
                "dgp": dgp_name,
                "eps_pool": diag_m["eps_pool"],
                "beta2": diag_m["beta2"],
                "eps_min": diag_m["eps_min"],
                "eps_max": diag_m["eps_max"],
                "Delta_eps": diag_m["Delta_eps"],
                "eps_pool_ek": diag_e["eps_pool"],
                "beta2_ek": diag_e["beta2"],
                "eps_min_ek": diag_e["eps_min"],
                "eps_max_ek": diag_e["eps_max"],
                "Delta_eps_ek": diag_e["Delta_eps"],
            }
        )
        print(
            f"{dgp_name:25s} {diag_m['eps_pool']:8.3f} {diag_m['beta2']:8.4f} "
            f"{diag_m['Delta_eps']:8.4f} {diag_e['eps_pool']:8.3f} {diag_e['beta2']:8.4f} "
            f"{diag_e['Delta_eps']:8.4f}"
        )

    tab2 = []
    print(
        f"\n{'DGP':25s} {'Lib':>4s} {'|ln_lam|':>8s} {'Err_EK%':>8s} "
        f"{'Err_Mel%':>8s} {'Bnd_EK%':>8s} {'Bnd_Mel%':>8s} {'OK_EK':>6s} {'OK_Mel':>7s}"
    )
    print("-" * 105)

    for dgp_name, dist in DGPS:
        econ = Economy(dist, N_COUNTRIES, SIGMA)
        diag_m = polynomial_diagnostic(econ, tau_vec, trade_flow="melitz")
        diag_e = polynomial_diagnostic(econ, tau_vec, trade_flow="ek")

        for tau_new, lib_label in LIBS:
            error_ek, ln_lam = econ.welfare_error(
                TAU_OLD, tau_new, diag_e["eps_pool"], model="ek"
            )
            error_mel, _ = econ.welfare_error(
                TAU_OLD, tau_new, diag_m["eps_pool"], model="melitz"
            )
            true_err_ek = abs(error_ek)
            true_err_mel = abs(error_mel)

            bnd_ek = diag_e["Delta_eps"] / (SIGMA - 1) * abs(ln_lam)

            bnd_mel = diag_m["Delta_eps"] / ((SIGMA - 1) * diag_m["eps_pool"]) * abs(ln_lam)
            holds_ek = true_err_ek <= bnd_ek + 1e-12
            holds_mel = true_err_mel <= bnd_mel + 1e-12
            tab2.append(
                {
                    "dgp": dgp_name,
                    "lib": lib_label,
                    "ln_lam": abs(ln_lam),
                    "error_ek_pct": true_err_ek * 100,
                    "error_melitz_pct": true_err_mel * 100,
                    "bound_ek_pct": bnd_ek * 100,
                    "bound_melitz_pct": bnd_mel * 100,
                    "holds_ek": holds_ek,
                    "holds_melitz": holds_mel,
                }
            )
            print(
                f"{dgp_name:25s} {lib_label:>4s} {abs(ln_lam):8.4f} "
                f"{true_err_ek * 100:7.3f}% {true_err_mel * 100:7.3f}% "
                f"{bnd_ek * 100:7.3f}% {bnd_mel * 100:7.3f}% "
                f"{'Y' if holds_ek else 'N':>6s} {'Y' if holds_mel else 'N':>7s}"
            )
        print()

    all_hold_ek = all(r["holds_ek"] for r in tab2)
    all_hold_mel = all(r["holds_melitz"] for r in tab2)
    print(f"ALL EK BOUNDS HOLD: {'YES' if all_hold_ek else 'NO'}")
    print(f"ALL MELITZ BOUNDS HOLD: {'YES' if all_hold_mel else 'NO'}")

    tab1_fields = [
        "dgp",
        "eps_pool",
        "beta2",
        "eps_min",
        "eps_max",
        "Delta_eps",
        "eps_pool_ek",
        "beta2_ek",
        "eps_min_ek",
        "eps_max_ek",
        "Delta_eps_ek",
    ]
    with open(os.path.join(OUTDIR_TAB, "table_mc_elasticities.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tab1_fields)
        w.writeheader()
        w.writerows(tab1)

    with open(os.path.join(OUTDIR_TAB, "table_mc_bounds.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tab2[0].keys())
        w.writeheader()
        w.writerows(tab2)

    save_elasticity_figure()
    print(f"\n{OUTDIR_FIG}/fig_elasticity_profile_ek.pdf")
    print(f"{OUTDIR_FIG}/fig_elasticity_profile_melitz.pdf")
    print(f"{OUTDIR_TAB}/ table_mc_*.csv")


if __name__ == "__main__":
    main()
