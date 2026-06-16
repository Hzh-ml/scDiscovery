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
    """
    返回(peak_count, mode_value)

    min_peak_ratio:
        第二峰高度 / 主峰高度 < 该值 → 视为单峰
        推荐:
            0.1 非常宽松（几乎都算单峰）
            0.2–0.3 工程稳健 ⭐推荐
            0.4+ 比较严格
    """
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

    # ⭐ NEW: 峰高度过滤（关键修复）
    peak_heights = dens[peaks]
    order = np.argsort(peak_heights)[::-1]

    main_peak = peak_heights[order[0]]

    if len(order) == 1:
        peak_cnt = 1
    else:
        second_peak = peak_heights[order[1]]
        ratio = second_peak / main_peak

        if ratio < min_peak_ratio:
            peak_cnt = 1   # ⭐ 视为单峰
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
    """
    scores: 1D MSP, 越大越像ID
    method:
      - "posterior_0.5": 用拟合出来的混合权重做 p(ID|s)=p(OOD|s) 边界（更保守）
      - "equal_likelihood": 解 N_id(s)=N_ood(s)（不带混合权重，阈值通常更低）
      - "posterior": 允许你指定 prior_id/prior_ood 或 eta 来调激进程度
          * 若给 prior_id/prior_ood：用它替代拟合权重
          * 若给 eta：用 likelihood ratio > eta（eta 越小越激进抓OOD）
    """
    x = np.asarray(scores, dtype=float).reshape(-1, 1)
    gmm = GaussianMixture(
        n_components=2, covariance_type="full",
        n_init=n_init, random_state=random_state, reg_covar=reg_covar
    ).fit(x)

    w = gmm.weights_.copy()
    mu = gmm.means_.ravel().copy()
    var = gmm.covariances_.ravel().copy()
    sig = np.sqrt(var)

    # MSP越大越ID => 均值更大的component当ID
    id_c = int(np.argmax(mu))
    ood_c = 1 - id_c

    mu_id, mu_ood = mu[id_c], mu[ood_c]
    s_id, s_ood = sig[id_c], sig[ood_c]

    # 决策中使用的“权重”（先验/代价折算）
    if method == "posterior_0.5":
        w_id, w_ood = w[id_c], w[ood_c]
        k = (w_ood / s_ood) / (w_id / s_id)  # 会把混合比例带进去（保守）
    elif method == "equal_likelihood":
        k = (1.0 / s_ood) / (1.0 / s_id)     # 等似然：忽略混合权重
    elif method == "posterior":
        # 你可指定先验；不指定就用拟合权重
        if prior_id is not None and prior_ood is not None:
            w_id, w_ood = float(prior_id), float(prior_ood)
        else:
            w_id, w_ood = w[id_c], w[ood_c]

        # eta 用来直接调 likelihood ratio 阈值：p(s|ood)/p(s|id) > eta
        # 若不给 eta，则默认等后验：eta = w_id/w_ood
        if eta is None:
            eta = w_id / w_ood

        # 从 p_ood/p_id > eta 变形得到： p_id * (w_id) ? 这里等价于把 k 设为 (1/eta)*(1/s_ood)/(1/s_id)
        k = (1.0 / float(eta)) * (1.0 / s_ood) / (1.0 / s_id)
    else:
        raise ValueError("method must be one of: posterior_0.5, equal_likelihood, posterior")

    # 解： (1/s_id) exp(-(t-mu_id)^2/(2 s_id^2))  =  k * (1/s_ood) exp(-(t-mu_ood)^2/(2 s_ood^2))
    # 取log后得到二次方程： a t^2 + b t + c = 0
    a = (1.0 / (2.0 * s_ood**2)) - (1.0 / (2.0 * s_id**2))
    b = (mu_id / (s_id**2)) - (mu_ood / (s_ood**2))
    c = (mu_ood**2) / (2.0 * s_ood**2) - (mu_id**2) / (2.0 * s_id**2) + np.log(k)

    eps = 1e-12
    if abs(a) < eps:
        # 退化为一次方程
        tau = -c / b
    else:
        disc = b*b - 4*a*c
        if disc < 0:
            # 数值兜底：网格找使两边最接近的点
            grid = np.linspace(float(x.min()), float(x.max()), 5000)
            left = (1.0/s_id) * np.exp(-(grid-mu_id)**2/(2*s_id**2))
            right = k * (1.0/s_ood) * np.exp(-(grid-mu_ood)**2/(2*s_ood**2))
            tau = float(grid[np.argmin(np.abs(left-right))])
        else:
            r1 = (-b + np.sqrt(disc)) / (2*a)
            r2 = (-b - np.sqrt(disc)) / (2*a)
            # 通常交点在两个均值之间，优先选这个
            lo, hi = min(mu_id, mu_ood), max(mu_id, mu_ood)
            cand = [r for r in (r1, r2) if lo <= r <= hi]
            if len(cand) == 1:
                tau = float(cand[0])
            elif len(cand) == 2:
                mid = 0.5*(mu_id+mu_ood)
                tau = float(min(cand, key=lambda r: abs(r-mid)))
            else:
                # 都不在区间内就选更靠近区间的
                tau = float(min((r1, r2), key=lambda r: min(abs(r-lo), abs(r-hi))))

    # 最终判别：MSP大->ID
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
    """
    简单逻辑：
      - KDE峰数==1 -> tau = mode - mode_offset
      - 否则 -> 2-GMM posterior=0.5 阈值
    MSP越大越ID：score>=tau 判ID，否则OOD
    """
    s = np.asarray(scores, dtype=float)
    s = s[np.isfinite(s)]
    if len(s) < 30:
        raise ValueError("样本太少，建议至少几十个分数。")

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

    # ===== 全局绘图样式统一 =====
    plt.rcParams["font.family"] = "DejaVu Serif"  # 可改成 "Arial" / "Times New Roman"
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

    # =========================
    # 单独导出 legend
    # =========================
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

    # ID 是正类 1，OOD 是负类 0
    scores = np.concatenate([id_scores, ood_scores])
    labels = np.concatenate([
        np.ones(len(id_scores)),
        np.zeros(len(ood_scores))
    ])

    # ROC 曲线
    fpr_curve, tpr_curve, thresholds = sk.roc_curve(labels, scores, pos_label=1)
    auroc = sk.roc_auc_score(labels, scores)

    # 计算你给定 threshold 下的 TPR / FPR
    id_pred_as_id = id_scores >= threshold
    ood_pred_as_id = ood_scores >= threshold

    tp = np.sum(id_pred_as_id)
    fn = np.sum(~id_pred_as_id)

    fp = np.sum(ood_pred_as_id)
    tn = np.sum(~ood_pred_as_id)

    tpr_at_threshold = tp / (tp + fn)
    fpr_at_threshold = fp / (fp + tn)

    # 画图
    plt.figure(figsize=(8, 6))
    plt.plot(fpr_curve, tpr_curve, linewidth=2, label=f"ROC curve, AUROC={auroc:.4f}")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Random")

    # 标出 threshold 对应的点
    plt.scatter(
        fpr_at_threshold,
        tpr_at_threshold,
        s=80,
        color="#800020",   # 酒红色
        zorder=10,   # 关键：让点压在线上面
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