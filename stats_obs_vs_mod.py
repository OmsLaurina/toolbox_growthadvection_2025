import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.close('all')


# ---- Data

# Modelised
datastation = pd.read_csv('outputs/adv_biom_AB.csv')
min_lat = 40.6
max_lat = 41.4
min_lon = 5
max_lon = 5.5

datastation = datastation[
    (datastation['Latitude'] >= min_lat) &
    (datastation['Latitude'] <= max_lat) &
    (datastation['Longitude'] >= min_lon) &
    (datastation['Longitude'] <= max_lon)
]


active_vars = ['biomass_E_neg', 'biomass_E_pos', 'biomass_T1', 'biomass_T2', 'Region']
phyto_datastation = datastation[active_vars]


biomass_columnsstation = [col for col in datastation.columns if col.startswith('biomass_')]
TotalBiomassstation = phyto_datastation[biomass_columnsstation].sum(axis=1)

# In situ
datastation_insitu = pd.read_csv('inputs/DATA_insitu/cyto_stations_AFB.csv')

T1 = ['HflrNano_biom', 'RedNano_biom']
T2 = ['RedPico_biom']
E_pos = ['HfNano_biom', 'HsNano_biom', 'HflrPico_biom']
E_neg = ['OraPicoProk_biom']

active_vars_insitu = ['E-', 'E+', 'T1', 'T2', 'Region']

datastation_insitu['T1'] = datastation_insitu[T1].sum(axis=1)
datastation_insitu['T2'] = datastation_insitu[T2].sum(axis=1)
datastation_insitu['E+'] = datastation_insitu[E_pos].sum(axis=1)
datastation_insitu['E-'] = datastation_insitu[E_neg].sum(axis=1)

phyto_datastation_insitu = datastation_insitu[active_vars_insitu]

# biomass_columnsstation_insitu = [col for col in datastation_insitu.columns if col.endswith('_biom')]
# TotalBiomassstation_insitu = phyto_datastation_insitu[biomass_columnsstation_insitu].sum(axis=1)

# ---- Calculations

# Proportion matrix
phyto_datastation[biomass_columnsstation] = phyto_datastation[biomass_columnsstation]#.div(TotalBiomassstation, axis=0)*100

# Mean proportion in each stations
mean_proportions = phyto_datastation.groupby(phyto_datastation['Region'])[biomass_columnsstation].mean()
mean_proportions = mean_proportions.reindex(['A', 'F', 'B'])

# Matrix for sns plot
data_melted = pd.melt(phyto_datastation, id_vars=['Region'], var_name='Group', value_name='Proportion')
data_melted_insitu = pd.melt(phyto_datastation_insitu, id_vars=['Region'], var_name='Group', value_name='Proportion')


# ---- Plot

### Set up

phyto_order = ['biomass_E_neg', 'biomass_E_pos', 'biomass_T1', 'biomass_T2']

phyto_labels = {
    'biomass_E_neg': 'E-',
    'biomass_E_pos': 'E+',
    'biomass_T1': 'T1',
    'biomass_T2': 'T2'
}

phyto_order_insitu = ['E-', 'E+', 'T1', 'T2']


data_plot = data_melted[data_melted['Region'] != '<undefined>'].copy()

data_plot['Region'] = pd.Categorical(data_plot['Region'], categories=['B', 'F', 'A'], ordered=True)

# region_colors = {'A': 'lightblue', 'F': 'lightcoral', 'B': 'lightgreen'}
region_colors_insitu = {'A2': 'lightblue', 'F2': 'lightcoral', 'B2': 'lightgreen'}

data_F_model = data_plot[data_plot['Region'] == 'F']

plt.figure(figsize=(4, 4))
ax = plt.gca()

sns.barplot(
    data=data_melted_insitu, 
    x='Group', 
    y='Proportion', 
    hue='Region', 
    palette=region_colors_insitu, 
    order=phyto_order_insitu,
    dodge=True,
    ax=ax,
    zorder=1,
)

positions = np.arange(len(phyto_order))
width = 0.2

for i, group in enumerate(phyto_order):
        f_model_value = data_F_model[data_F_model['Group'] == group]['Proportion'].values[0]

        f_position = positions[i]
        print(f_position)

        ax.bar(
            x=f_position,
            height=f_model_value,
            width=width,
            color='none',
            edgecolor='black',
            linestyle='--',
            hatch='///',     
            alpha=0.7,
            label="F Modèle" if i == 0 else "",  
            zorder=2
        )
ax.set_xlabel('')
ax.set_ylabel('')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, ['St.B2', 'St.F2', 'St.A2', 'F sim'], title="", loc='upper right', fontsize=10, frameon=True, framealpha=1, facecolor='white', edgecolor='black')

for spine in ["top", "bottom", "right", "left"]:
    ax.spines[spine].set_visible(False)

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xticks(range(len(phyto_order)))


plt.show()

plt.savefig("outputs/barplot_Fsim_vs_Fobs.tif", format="tiff", dpi=600)


plt.figure()
plt.scatter(datastation.Longitude, datastation.Latitude, c=datastation['S'])
