# losses.py
import torch
import torch.nn.functional as F

def ce_loss(logits, labels):
    return F.cross_entropy(logits, labels)

def msp_score(logits):
    probs = F.softmax(logits, dim=1)
    return probs.max(dim=1).values

def oe_loss(ood_logits):
    # Outlier Exposure: force uniform prediction
    loss_oe = 0.5 * -(ood_logits.mean(1) - torch.logsumexp(ood_logits, dim=1)).mean()
    return loss_oe

def align_loss(z1, z2):
    return F.mse_loss(z1, z2)


class RAEReconstructionLoss:
    def __init__(self, lambda_latent=0.0):
        self.lambda_latent = lambda_latent

    def __call__(
        self,
        x_rna,
        x_atac,
        z_rna,
        z_atac,
        D_rna,
        D_atac,
        E_rna,
        E_atac,
    ):
        # -------- decode --------
        x_rna_hat = F.softplus(D_rna(z_rna))
        x_atac_hat = D_atac(z_atac)

        x_atac_from_rna = D_atac(z_rna)
        x_rna_from_atac = F.softplus(D_rna(z_atac))

        # -------- reconstruction (SAFE) --------
        loss_rna = F.l1_loss(x_rna_hat, x_rna)
        loss_atac = F.l1_loss(x_atac_hat, x_atac)

        loss_cross = (
            F.l1_loss(x_atac_from_rna, x_atac)
            + F.l1_loss(x_rna_from_atac, x_rna)
        )

        # -------- latent (OPTIONAL) --------
        if self.lambda_latent > 0:
            z_rna_cycle = E_rna(x_rna_from_atac.detach())[0]
            z_atac_cycle = E_atac(x_atac_from_rna.detach())[0]
            loss_latent = (
                F.mse_loss(z_rna, z_rna_cycle)
                + F.mse_loss(z_atac, z_atac_cycle)
            )
        else:
            loss_latent = 0.0

        total = loss_rna + loss_atac + loss_cross + self.lambda_latent * loss_latent

        return total, {
            "rna": loss_rna.item(),
            "atac": loss_atac.item(),
            "cross": loss_cross.item(),
        }
