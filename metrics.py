import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    classification_report,
    f1_score,
    accuracy_score,
    silhouette_score,
    calinski_harabasz_score
)


def evaluate_ncd_discovery_potential(adata, y_true, y_pred, embed_key='X_glue', unknown_label='Unknown'):
    """
    NCD 评估函数：兼顾“未知类识别准确率”与“新类潜在结构质量”

    参数
    ----
    adata         : AnnData (必须包含 embed_key)
    y_true        : str, 真实标签列 (含具体已知类 + 'Unknown')
    y_pred        : str, 预测标签列 (含具体已知类 + 'Unknown')
    embed_key     : str, 嵌入用于评估结构
    unknown_label : str, 标签中代表未知的字符串

    返回
    ----
    dict : 包含 Coarse Metrics (识别) 和 Discovery Metrics (发现潜力)
    """

    # 0. 数据准备
    y_true_arr = adata.obs[y_true].astype(str).values
    y_pred_arr = adata.obs[y_pred].astype(str).values
    emb = adata.obsm[embed_key]

    metrics = {}

    # =========================================================
    # Part 1: Coarse-grained Metrics (粗粒度识别能力)
    # 评估：模型能不能把 Known 分对，能不能把 Unknown 挑出来？
    # =========================================================

    # 1.1 宏观 F1 (Weighted/Macro) - 关注整体分类性能
    # 这里的类别集合是 [Known_A, Known_B, ..., Unknown]
    metrics['Coarse_F1_Macro'] = f1_score(y_true_arr, y_pred_arr, average='macro')
    metrics['Coarse_F1_weighted'] = f1_score(y_true_arr, y_pred_arr, average='weighted')
    metrics['Coarse_ACC'] = accuracy_score(y_true_arr, y_pred_arr)

    # 1.2 专门针对 "Unknown" 的识别能力
    # 把问题转化为二分类：Unknown vs Known
    true_is_unk = (y_true_arr == unknown_label)
    pred_is_unk = (y_pred_arr == unknown_label)

    # Unknown 类的 F1 分数 (最重要的单一分类指标)
    metrics['Detection_F1_Unknown'] = f1_score(true_is_unk, pred_is_unk, pos_label=True)  # average='weighted',

    # 纯度与查全率 (辅助分析)
    # Precision: 预测为 Unknown 的样本里，真的 Unknown 占多少？(越低说明误伤了已知类)
    if pred_is_unk.sum() > 0:
        metrics['Detection_Prec_Unknown'] = (true_is_unk & pred_is_unk).sum() / pred_is_unk.sum()
    else:
        metrics['Detection_Prec_Unknown'] = 0.0

    # Recall: 真实的 Unknown 样本里，被找出来多少？(越低说明把新类误判为已知类)
    if true_is_unk.sum() > 0:
        metrics['Detection_Recall_Unknown'] = (true_is_unk & pred_is_unk).sum() / true_is_unk.sum()
    else:
        metrics['Detection_Recall_Unknown'] = 0.0

    # =========================================================
    # Part 2: Fine-grained Discovery Potential (细粒度发现潜力)
    # 策略 1: 评估【预测为未知】的区域 (User's original, dirty but useful for output check)
    # 这代表了“用户最终拿到手的那堆 Unknown 数据的质量”
    metrics.update(_calculate_cluster_metrics(
        emb, pred_is_unk, prefix='Pred_Unk'
    ))

    # 策略 2: 评估【被正确捕获的真实未知】区域 (Intersection, Pure)
    # 这代表了“模型真正发现的新类样本的质量”，排除了已知类的干扰
    mask_pure_discovery = true_is_unk & pred_is_unk
    metrics.update(_calculate_cluster_metrics(
        emb, mask_pure_discovery, prefix='Pure_Unk'
    ))

    return metrics


def _calculate_cluster_metrics(emb, mask, prefix):
    """辅助函数：对指定 mask 的样本进行 Leiden 聚类并计算指标"""
    import scanpy as sc
    res = {}

    # 样本太少不计算
    if mask.sum() < 20:
        res[f'{prefix}_ASW'] = 0.0
        res[f'{prefix}_N_Clusters'] = 0
        return res

    emb_subset = emb[mask]

    # 虚拟聚类
    try:
        # 创建轻量级临时对象
        temp_adata = sc.AnnData(X=np.zeros((emb_subset.shape[0], 1)))
        temp_adata.obsm['X_emb'] = emb_subset

        # 聚类流程
        sc.pp.neighbors(temp_adata, use_rep='X_emb', n_neighbors=15)
        sc.tl.leiden(temp_adata, resolution=0.5, key_added='sub')

        labels = temp_adata.obs['sub'].values
        n_clusters = len(np.unique(labels))

        res[f'{prefix}_N_Clusters'] = n_clusters

        if n_clusters > 1:
            res[f'{prefix}_ASW'] = silhouette_score(emb_subset, labels)
            res[f'{prefix}_CH'] = calinski_harabasz_score(emb_subset, labels)
        else:
            res[f'{prefix}_ASW'] = 0.0
            res[f'{prefix}_CH'] = 0.0

    except Exception as e:
        print(f"Clustering failed for {prefix}: {e}")
        res[f'{prefix}_ASW'] = 0.0

    return res


def calculate_discovery_asw(X_emb, n_neighbors=15, resolution=0.5):
    """
    对给定的嵌入空间(X_emb)进行无监督聚类(Leiden)，
    并计算轮廓系数(ASW)以评估其内部结构的清晰度。
    """
    # 1. 样本量检查：太少无法计算聚类或轮廓系数
    if X_emb.shape[0] < 10:
        return 0.0

    try:
        # 2. 创建轻量级临时对象，避免污染原数据
        # X 只是占位符，实际运算全靠 obsm
        temp_adata = sc.AnnData(X=np.zeros((X_emb.shape[0], 1)))
        temp_adata.obsm['X_emb'] = X_emb

        # 3. 运行标准的 Scanpy 聚类流程
        # 计算邻居图
        sc.pp.neighbors(temp_adata, use_rep='X_emb', n_neighbors=n_neighbors)
        # 运行 Leiden 聚类
        sc.tl.leiden(temp_adata, resolution=resolution, key_added='leiden_temp')

        # 4. 提取聚类标签
        labels = temp_adata.obs['leiden_temp'].values
        n_clusters = len(np.unique(labels))

        # 5. 计算 ASW (仅当至少发现了 2 个子群时才有意义)
        if n_clusters > 1:
            score = silhouette_score(X_emb, labels)
            return score
        else:
            # 如果只聚成了一类，说明结构单一或无法区分，得分为 0
            return 0.0

    except Exception as e:
        print(f"结构评分计算出错: {e}")
        return 0
