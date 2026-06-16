import torch
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.mixture import GaussianMixture
import math
import numpy as np
import torch.nn.functional as F
import sklearn.metrics as sk
import matplotlib.pyplot as plt

concat = lambda x: np.concatenate(x, axis=0)
to_np = lambda x: x.data.cpu().numpy()


def get_ood_scores_with_single_omics(E_rna, Classifier, omics1_loader, device):
    E_rna.eval()
    Classifier.eval()

    omics1_score = []
    omics1_data = []
    omics1_raw_data = []
    omics1_label = []

    with torch.no_grad():
        for x_raw_1, x_1, y_1 in omics1_loader:
            z_1, _, _ = E_rna(x_1.to(device))

            output = Classifier(z_1)

            smax = to_np(F.softmax(output, dim=1))

            scores = np.max(smax, axis=1)  # np.max(smax, axis=1)

            omics1_data.append(x_1)

            omics1_raw_data.append(x_raw_1)

            omics1_label.append(y_1)

            omics1_score.append(scores[:len(x_1)])

            E_rna.zero_grad()
            Classifier.zero_grad()

    return concat(omics1_score).copy(), torch.cat(omics1_data, dim=0), torch.cat(omics1_label, dim=0), torch.cat(omics1_raw_data, dim=0)


def _count_peaks_kde(scores, grid_size=2048, bandwidth=None, min_peak_ratio=0.2):
    s = np.asarray(scores, dtype=float)
    s = s[np.isfinite(s)]

    if bandwidth is None:
        std = float(np.std(s))
        n = len(s)
        bw = 1.06 * std * (n ** (-1 / 5)) if std > 0 else 0.01
        bandwidth = float(np.clip(bw, 1e-3, 0.2))

    grid = np.linspace(float(s.min()), float(s.max()), grid_size)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(s.reshape(-1, 1))
    dens = np.exp(kde.score_samples(grid.reshape(-1, 1)))

    peaks = np.where((dens[1:-1] > dens[:-2]) & (dens[1:-1] > dens[2:]))[0] + 1

    if len(peaks) == 0:
        mode_idx = int(np.argmax(dens))
        return 1, float(grid[mode_idx]), bandwidth

    peak_heights = dens[peaks]
    order = np.argsort(peak_heights)[::-1]

    main_peak = peak_heights[order[0]]

    if len(order) == 1:
        peak_cnt = 1
    else:
        second_peak = peak_heights[order[1]]
        ratio = second_peak / main_peak

        if ratio < min_peak_ratio:
            peak_cnt = 1
        else:
            peak_cnt = len(peaks)

    mode_idx = peaks[order[0]]
    return int(peak_cnt), float(grid[mode_idx]), bandwidth


def floor_to_one_decimal(eta_0):
    eta = min(max(eta_0, 0.1), 0.9)
    eta = math.floor((eta + 1e-12) * 10) / 10
    return eta


def gmm_threshold_msp(scores,
                      n_init=20,
                      random_state=0,
                      reg_covar=1e-6,
                      mode_offset=0.1,
                      method="equal_likelihood",
                      prior_id=None,
                      prior_ood=None,
                      eta=None):
    x = np.asarray(scores, dtype=float).reshape(-1, 1)
    gmm = GaussianMixture(
        n_components=2, covariance_type="full",
        n_init=n_init, random_state=random_state, reg_covar=reg_covar
    ).fit(x)

    w = gmm.weights_.copy()
    mu = gmm.means_.ravel().copy()
    var = gmm.covariances_.ravel().copy()
    sig = np.sqrt(var)

    id_c = int(np.argmax(mu))
    ood_c = 1 - id_c

    mu_id, mu_ood = mu[id_c], mu[ood_c]
    s_id, s_ood = sig[id_c], sig[ood_c]

    if method == "posterior_0.5":
        w_id, w_ood = w[id_c], w[ood_c]
        k = (w_ood / s_ood) / (w_id / s_id)
    elif method == "equal_likelihood":
        k = (1.0 / s_ood) / (1.0 / s_id)
    elif method == "posterior":
        if prior_id is not None and prior_ood is not None:
            w_id, w_ood = float(prior_id), float(prior_ood)
        else:
            w_id, w_ood = w[id_c], w[ood_c]

        if eta is None:
            eta = w_id / w_ood

        k = (1.0 / float(eta)) * (1.0 / s_ood) / (1.0 / s_id)
    else:
        raise ValueError("method must be one of: posterior_0.5, equal_likelihood, posterior")

    a = (1.0 / (2.0 * s_ood**2)) - (1.0 / (2.0 * s_id**2))
    b = (mu_id / (s_id**2)) - (mu_ood / (s_ood**2))
    c = (mu_ood**2) / (2.0 * s_ood**2) - (mu_id**2) / (2.0 * s_id**2) + np.log(k)

    eps = 1e-12
    if abs(a) < eps:
        tau = -c / b
    else:
        disc = b*b - 4*a*c
        if disc < 0:
            grid = np.linspace(float(x.min()), float(x.max()), 5000)
            left = (1.0/s_id) * np.exp(-(grid-mu_id)**2/(2*s_id**2))
            right = k * (1.0/s_ood) * np.exp(-(grid-mu_ood)**2/(2*s_ood**2))
            tau = float(grid[np.argmin(np.abs(left-right))])
        else:
            r1 = (-b + np.sqrt(disc)) / (2*a)
            r2 = (-b - np.sqrt(disc)) / (2*a)
            lo, hi = min(mu_id, mu_ood), max(mu_id, mu_ood)
            cand = [r for r in (r1, r2) if lo <= r <= hi]
            if len(cand) == 1:
                tau = float(cand[0])
            elif len(cand) == 2:
                mid = 0.5*(mu_id+mu_ood)
                tau = float(min(cand, key=lambda r: abs(r-mid)))
            else:
                tau = float(min((r1, r2), key=lambda r: min(abs(r-lo), abs(r-hi))))

    scores_1d = x.ravel()
    pred_is_id = scores_1d >= tau
    tau = tau - mode_offset

    return tau, {
        "method": method,
        "tau": tau,
        "gmm_weights": w.tolist(),
        "means": mu.tolist(),
        "sigmas": sig.tolist(),
        "id_component": id_c,
        "pred_label": np.where(pred_is_id, "ID", "OOD"),
        "gmm": gmm,
    }


def threshold_msp_simple(scores,
                         mode_offset=0.1,
                         kde_bandwidth=None,
                         gmm_n_init=20,
                         random_state=0):
    s = np.asarray(scores, dtype=float)
    s = s[np.isfinite(s)]
    if len(s) < 30:
        raise ValueError("too little data")

    peak_cnt, mode_val, bw = _count_peaks_kde(s, bandwidth=kde_bandwidth)

    s_min, s_max = float(s.min()), float(s.max())

    if peak_cnt == 1:
        tau = float(np.clip(mode_val - mode_offset, s_min, s_max))
        method = "unimodal: tau = mode - offset"
        info = {"peak_count": peak_cnt, "mode": mode_val, "bandwidth": bw, "mode_offset": mode_offset}
    else:
        tau, info = gmm_threshold_msp(scores, mode_offset=mode_offset, method='equal_likelihood')  # (equal_likelihood) (posterior_0.5) (method="posterior", prior_id=0.5, prior_ood=0.5) (method="posterior", eta=0.3)
        method = "multimodal: 2-GMM posterior=0.5"
        # info = {"peak_count": peak_cnt, "mode": mode_val, "bandwidth": bw, **gmm_info}

    pred = np.where(s >= tau, "ID", "OOD")

    tau = floor_to_one_decimal(tau)

    # print("tau =", tau)

    return tau, pred, {"method": method, **info}


def get_scores_and_density_plot(E, Classifier, omics1_test_id_loader, omics1_test_ood_loader, omics1_train_loader, dataset_name, setting, device):
    if omics1_test_id_loader is None:
        omics1_test_id_loader = omics1_train_loader

    id_scores, _, _, _ = get_ood_scores_with_single_omics(E, Classifier, omics1_test_id_loader, device)
    ood_scores, _, _, _ = get_ood_scores_with_single_omics(E, Classifier, omics1_test_ood_loader, device)

    scores = np.concatenate((id_scores, ood_scores), axis=0)

    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Patch

    plt.rcParams["font.family"] = "DejaVu Serif"
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.titlesize"] = 20
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["xtick.labelsize"] = 20
    plt.rcParams["ytick.labelsize"] = 20
    plt.rcParams["legend.fontsize"] = 15
    plt.rcParams["legend.title_fontsize"] = 16

    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"

    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1.2
    plt.rcParams["grid.alpha"] = 0.0

    plt.figure(figsize=(6, 6))  # (10, 6)

    sns.kdeplot(
        id_scores,
        label="Known",
        fill=True,
        color="blue"
    )

    sns.kdeplot(
        ood_scores,
        label="Unknown",
        fill=True,
        color="red"
    )

    plt.xlabel("Maximum Softmax Probability")
    plt.ylabel("Density")
    plt.grid(False)

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_visible(True)

    plt.savefig(f"./figures/scDiscovery_density_{dataset_name}_{setting}.pdf", dpi=600, bbox_inches="tight")
    plt.show()

    legend_elements = [
        Patch(facecolor="blue", edgecolor="blue", label="Known", alpha=0.6),
        Patch(facecolor="red", edgecolor="red", label="Unknown", alpha=0.6),
    ]

    fig_legend = plt.figure(figsize=(3, 1.2))
    fig_legend.legend(
        handles=legend_elements,
        loc="center",
        ncol=2,
        framealpha=0.1,
        frameon=False
    )

    fig_legend.savefig(
        f"scDiscovery_legend_{dataset_name}_{setting}.pdf",
        dpi=600,
        bbox_inches="tight",
        transparent=True
    )
    plt.close(fig_legend)

    return scores, id_scores, ood_scores


def plot_roc_with_threshold(id_scores, ood_scores, threshold, save_path=None):
    id_scores = np.asarray(id_scores).reshape(-1)
    ood_scores = np.asarray(ood_scores).reshape(-1)

    scores = np.concatenate([id_scores, ood_scores])
    labels = np.concatenate([
        np.ones(len(id_scores)),
        np.zeros(len(ood_scores))
    ])

    fpr_curve, tpr_curve, thresholds = sk.roc_curve(labels, scores, pos_label=1)
    auroc = sk.roc_auc_score(labels, scores)

    id_pred_as_id = id_scores >= threshold
    ood_pred_as_id = ood_scores >= threshold

    tp = np.sum(id_pred_as_id)
    fn = np.sum(~id_pred_as_id)

    fp = np.sum(ood_pred_as_id)
    tn = np.sum(~ood_pred_as_id)

    tpr_at_threshold = tp / (tp + fn)
    fpr_at_threshold = fp / (fp + tn)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_curve, tpr_curve, linewidth=2, label=f"ROC curve, AUROC={auroc:.4f}")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Random")

    plt.scatter(
        fpr_at_threshold,
        tpr_at_threshold,
        s=80,
        color="#800020",
        zorder=10,
        label=f"threshold={threshold:.4f}\nFPR={fpr_at_threshold:.4f}, TPR={tpr_at_threshold:.4f}"
    )

    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    # plt.title("ROC Curve with Threshold")
    plt.legend()
    plt.grid(True)

    if save_path is not None:
        plt.savefig(save_path, dpi=600, bbox_inches="tight")

    plt.show()

    return {
        "auroc": auroc,
        "threshold": threshold,
        "tpr": tpr_at_threshold,
        "fpr": fpr_at_threshold,
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
    }