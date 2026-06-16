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

        self.fc1 = nn.Linear(input_size, hidden_size)
        self.bn1 = nn.BatchNorm1d(hidden_size)
        self.dropout1 = nn.Dropout(dropout_rate)

        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.bn2 = nn.BatchNorm1d(hidden_size)
        self.dropout2 = nn.Dropout(dropout_rate)

        self.fc_embed = nn.Linear(hidden_size, z_dim)

        self.fc_mu = nn.Linear(z_dim, z_dim)
        self.fc_logvar = nn.Linear(z_dim, z_dim)


    def forward(self, x):
        h = self.fc1(x)
        h = self.bn1(h)
        h = F.relu(h)
        h = self.dropout1(h)

        h = self.fc2(h)
        h = self.bn2(h)
        h = F.relu(h)
        h = self.dropout2(h)

        embedding = self.fc_embed(h)

        mu = self.fc_mu(embedding)
        logvar = self.fc_logvar(embedding)

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
            input_size,
            hidden_size=256,
            z_dim=128,
            dropout_rate=0.5
    ):

        super(Decoder, self).__init__()

        self.fc1 = nn.Linear(z_dim, hidden_size)
        self.bn1 = nn.BatchNorm1d(hidden_size)
        self.dropout1 = nn.Dropout(dropout_rate)

        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.bn2 = nn.BatchNorm1d(hidden_size)
        self.dropout2 = nn.Dropout(dropout_rate)

        self.fc3 = nn.Linear(hidden_size, input_size)

    def forward(self, z):
        x = self.fc1(z)
        x = self.bn1(x)
        x = F.relu(x)
        x = self.dropout1(x)

        x = self.fc2(x)
        x = self.bn2(x)
        x = F.relu(x)
        x = self.dropout2(x)

        x_hat = self.fc3(x)
        return x_hat
