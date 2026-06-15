import torch
import torch.nn as nn
import torch.nn.functional as F
import math


def artanh(x, eps=1e-5):
    x = torch.clamp(x, -1 + eps, 1 - eps)
    return 0.5 * (torch.log1p(x) - torch.log1p(-x))


def project_poincare(x, c=1.0, eps=1e-5):
    sqrt_c = math.sqrt(c)
    max_norm = (1.0 - eps) / sqrt_c

    norm = x.norm(dim=-1, keepdim=True).clamp_min(eps)
    projected = x / norm * max_norm

    return torch.where(norm > max_norm, projected, x)


def expmap0(u, c=1.0, eps=1e-5):
    """
    Map tangent vector u at origin to Poincare ball.
    """
    sqrt_c = math.sqrt(c)
    norm_u = u.norm(dim=-1, keepdim=True).clamp_min(eps)

    z = torch.tanh(sqrt_c * norm_u) * u / (sqrt_c * norm_u)
    z = project_poincare(z, c=c, eps=eps)

    return z


def logmap0(z, c=1.0, eps=1e-5):
    """
    Map Poincare ball point z back to tangent space at origin.
    """
    sqrt_c = math.sqrt(c)
    norm_z = z.norm(dim=-1, keepdim=True).clamp_min(eps)

    u = artanh(sqrt_c * norm_z, eps=eps) * z / (sqrt_c * norm_z)

    return u


class Encoder(nn.Module):
    def __init__(
        self,
        input_size,
        hidden_size=256,
        z_dim=128,
        e_dim=64,
        h_dim=64,
        dropout_rate=0.5,
        curvature=1.0,
    ):
        """
        VAE Encoder with explicit embedding output

        Outputs:
        --------
        embedding : deterministic embedding (same role as old encoder output)
        mu        : mean of q(z|x)
        logvar    : log variance of q(z|x)
        """
        super(Encoder, self).__init__()

        # 第一层
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.bn1 = nn.BatchNorm1d(hidden_size)
        self.dropout1 = nn.Dropout(dropout_rate)

        # 第二层
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.bn2 = nn.BatchNorm1d(hidden_size)
        self.dropout2 = nn.Dropout(dropout_rate)

        # embedding 层（对应你原来的 encoder.fc3）
        self.fc_embed = nn.Linear(hidden_size, z_dim)

        # VAE 的两个 head
        self.fc_mu = nn.Linear(z_dim, z_dim)
        self.fc_logvar = nn.Linear(z_dim, z_dim)

        # # Euclidean projection head
        # self.euclidean_head = nn.Sequential(
        #     nn.Linear(z_dim, z_dim),
        #     nn.ReLU(),
        #     nn.Dropout(dropout_rate),
        #     nn.Linear(z_dim, e_dim),
        # )
        #
        # # Hyperbolic projection head
        # # 先输出切空间向量，再通过 expmap0 映射到双曲空间
        # self.hyperbolic_head = nn.Sequential(
        #     nn.Linear(z_dim, z_dim),
        #     nn.ReLU(),
        #     nn.Dropout(dropout_rate),
        #     nn.Linear(z_dim, h_dim),
        # )
        #
        # self.curvature = curvature

    def forward(self, x):
        # 第一层
        h = self.fc1(x)
        h = self.bn1(h)
        h = F.relu(h)
        h = self.dropout1(h)

        # 第二层
        h = self.fc2(h)
        h = self.bn2(h)
        h = F.relu(h)
        h = self.dropout2(h)

        # deterministic embedding（原 encoder 的输出语义）
        embedding = self.fc_embed(h)

        # 变分参数
        mu = self.fc_mu(embedding)
        logvar = self.fc_logvar(embedding)

        # # 欧氏表示
        # z_e = self.euclidean_head(embedding)
        #
        # # 双曲表示
        # u_h = self.hyperbolic_head(embedding)
        # z_h = expmap0(u_h, c=self.curvature)

        return embedding, mu, logvar  # , z_e, z_h


def reparameterize(mu, logvar):
    std = torch.exp(0.5 * logvar)
    eps = torch.randn_like(std)
    return mu + eps * std


class JointClassificationHead(nn.Module):
    """
    Shared classification head
    """
    def __init__(self, z_dim, n_classes):
        super().__init__()
        self.cls = nn.Linear(z_dim, n_classes)
        self.ood = nn.Linear(z_dim, 1)

    def forward(self, z):
        logits = self.cls(z)
        return logits


class Decoder(nn.Module):
    def __init__(
            self,
            input_size,  # encoder 的 input_size（要重建的维度）
            hidden_size=256,
            z_dim=128,  # encoder 的 output_size
            dropout_rate=0.5
    ):
        """
        对称 Decoder：z -> x_hat
        """
        super(Decoder, self).__init__()

        # 第一层（对应 encoder.fc3 的反向）
        self.fc1 = nn.Linear(z_dim, hidden_size)
        self.bn1 = nn.BatchNorm1d(hidden_size)
        self.dropout1 = nn.Dropout(dropout_rate)

        # 第二层（对应 encoder.fc2 的反向）
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.bn2 = nn.BatchNorm1d(hidden_size)
        self.dropout2 = nn.Dropout(dropout_rate)

        # 第三层（对应 encoder.fc1 的反向）
        self.fc3 = nn.Linear(hidden_size, input_size)

    def forward(self, z):
        # 第一层
        x = self.fc1(z)
        x = self.bn1(x)
        x = F.relu(x)
        x = self.dropout1(x)

        # 第二层
        x = self.fc2(x)
        x = self.bn2(x)
        x = F.relu(x)
        x = self.dropout2(x)

        # 第三层（重建输出，不加激活）
        x_hat = self.fc3(x)
        return x_hat
