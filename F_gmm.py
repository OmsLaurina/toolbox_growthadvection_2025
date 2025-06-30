import numpy as np
from scipy.io import loadmat, savemat
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd

plt.close('all')

# --- PARAMÈTRES ---
input_mat_file = 'outputs/output2D_AFB.mat'
output_mat_file = 'outputs/output2D_region.mat'
structure_variable_name = 'output2D'
region_variable_name = 'Region'
lon_threshold = 5
boundaries = [38.17, 38.42]  # Fixés à partir du GMM (ou définis a priori)

# --- CHARGEMENT DES DONNÉES ---
data = loadmat(input_mat_file, struct_as_record=False, squeeze_me=True)
output2D = data[structure_variable_name]

salinity = output2D.S  # shape: (lon, lat)
lat = output2D.Latitude       # shape: (lat,)
lon = output2D.Longitude      # shape: (lon,)

# --- AJUSTEMENT GMM AVEC EXCLUSION DES NaNs ---
salinity_flat = salinity.flatten()
valid_mask = ~np.isnan(salinity_flat)
salinity_valid = salinity_flat[valid_mask].reshape(-1, 1)

gmm = GaussianMixture(n_components=2, covariance_type='full', max_iter=1000,
                      reg_covar=1e-6, random_state=42)
gmm.fit(salinity_valid)

means = gmm.means_.flatten()
std_devs = np.sqrt(gmm.covariances_.flatten())

print("Paramètres GMM :")
for i, (mean, std) in enumerate(zip(means, std_devs), 1):
    interval = (mean - std, mean + std)
    print(f"  - Composante {i} : µ = {mean:.3f}, σ = {std:.3f}, intervalle 1σ = [{interval[0]:.3f}, {interval[1]:.3f}]")

# --- CONSTRUCTION DE LA MATRICE DES RÉGIONS ---
region_matrix = np.zeros_like(salinity, dtype=int)

for i in range(salinity.shape[0]):  # lon
    for j in range(salinity.shape[1]):  # lat
        s = salinity[i, j]
        if np.isnan(s):
            region_matrix[i, j] = 0  # Non défini
            continue
        lon_val = lon[j]
        if s >= boundaries[1]:
            region_matrix[i, j] = 1  # Région A
        elif s <= boundaries[0]:
            region_matrix[i, j] = 2  # Région B
        else:
              region_matrix[i, j] = 3 #if lon_val > lon_threshold else 2  # Région F ou B

# --- AJOUT À LA STRUCTURE ET SAUVEGARDE ---
output2D.Region_gmm = region_matrix
savemat(output_mat_file, {structure_variable_name: output2D})

# --- TRAÇAGE DE LA CARTE ---
cmap = ListedColormap(['blue', 'green', 'red'])  # A=1, F=2, B=3
masked_region = np.ma.masked_where(region_matrix == 0, region_matrix)

plt.figure(figsize=(10, 6))
plt.pcolormesh(lon, lat, masked_region.T, cmap=cmap, shading='auto')
plt.colorbar(ticks=[1, 2, 3], label='Région (1=A, 2=B, 3=F)', format='%d')
plt.clim(1, 3)
plt.title('Classification des régions (GMM)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()
plt.show()
plt.savefig("figures/F_gmm.tif", format="tiff", dpi=600)


# Create Gmax from F positions
Gmax1 = np.full_like(region_matrix, np.nan, dtype=float)
Gmax2 = np.full_like(region_matrix, np.nan, dtype=float)

Gmax1[region_matrix == 1] = 3.89
Gmax2[region_matrix == 1] = 0.43

Gmax1[region_matrix == 3] = 3.89 + 1
Gmax2[region_matrix == 3] = 0.43 + 1

Gmax1[region_matrix == 2] = 3.89
Gmax2[region_matrix == 2] = 0.43

df_Gmax1 = pd.DataFrame(Gmax1)
df_Gmax2 = pd.DataFrame(Gmax2)

df_Gmax1.to_csv("outputs/Gmax1_AFB.csv", index=False, header=False)
df_Gmax2.to_csv("outputs/Gmax2_AFB.csv", index=False, header=False)
