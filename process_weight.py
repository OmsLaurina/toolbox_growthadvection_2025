import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.stats import gaussian_kde

plt.close('all')

# Fonction utilitaire pour charger les sorties
def load_output(path):
    return loadmat(path, struct_as_record=False, squeeze_me=True)['final_map']

# Chargement des données
output2D_passive = loadmat('outputs/output2D_region.mat', struct_as_record=False, squeeze_me=True)['output2D']
output_NPZ = {
    'Actif': load_output('outputs/output2D_AB_NPZ_active.mat'),
    'Réactif': load_output('outputs/output2D_AB_NPZ_reactive.mat'),
    'Actif+Réactif': load_output('outputs/output2D_AB_NPZ_activereactive.mat')
}

region = np.array(output2D_passive.Region_flag)
S_passive = np.array(output2D_passive.S.mean)

# Fonction de calcul des R_exp et R_NPZ pour un scénario donné
def compute_R_ratios(output_NPZ_i):
    S = np.array(output2D_passive.S.mean)
    P1 = np.array(output_NPZ_i.P1)
    P2 = np.array(output_NPZ_i.P2)
    R_NPZ_F = np.array(output_NPZ_i.R_P2)[region == 3]

    S_A = S[region == 1][1]
    S_B = S[region == 2][1]
    S_F = S[region == 3]
    alpha = (S_F - S_B) / (S_A - S_B)

    P1_A = P1[region == 1][1]
    P1_B = P1[region == 2][1]
    P2_A = P2[region == 1][1]
    P2_B = P2[region == 2][1]

    P1_exp = P1_B * (1 - alpha) + P1_A * alpha
    P2_exp = P2_B * (1 - alpha) + P2_A * alpha

    R_exp_F = P2_exp / (P1_exp + P2_exp)

    return R_exp_F, R_NPZ_F, alpha

# Préparation des KDE pour chaque scénario
colors = {'exp': 'orange', 'npz': 'steelblue'}

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

for ax, (key, output) in zip(axes, output_NPZ.items()):
    
    if ax == axes[0]:
        x = np.linspace(0., 0.7, 500)
        ax.set_ylim(0, 10)
    if ax == axes[1]:
        x = np.linspace(0.6, 1, 500)
    if ax == axes[2]:
        x = np.linspace(0.6, 1, 500)
        
    R_exp, R_npz, alpha = compute_R_ratios(output)

    kde_exp = gaussian_kde(R_exp)
    kde_npz = gaussian_kde(R_npz)

    pdf_exp = kde_exp(x)
    pdf_npz = kde_npz(x)

    ax.plot(x, pdf_exp, label='R passive', color=colors['exp'], linewidth=2)
    ax.plot(x, pdf_npz, label='R NPZ', color=colors['npz'], linewidth=2)
    ax.fill_between(x, pdf_exp, pdf_npz, color='gray', alpha=0.3)

    area_diff = np.trapz(np.abs(pdf_exp - pdf_npz), x)
    ax.set_title(f'{key} --- ∆ aire = {area_diff:.3f}')
    ax.set_xlabel('R = P2 / (P1 + P2)')
    ax.grid(True)
    if ax == axes[0]:
        ax.set_ylabel('Densité')
        ax.legend(loc='upper left', ncol=2)

plt.tight_layout()
plt.show()
plt.savefig("figures/weightprocess.tif", format="tiff", dpi=600)
