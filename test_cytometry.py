import pyVIA.core as via
import pyVIA.datasets_via as datasets_via
import pandas as pd
import numpy as np
import scanpy as sc     
import umap
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

cytometry_feature_csv_file_path = '/Users/chienyuehliu/Downloads/mcf7_38features.csv'
cytometry_phase_csv_file_path = '/Users/chienyuehliu/Downloads/mcf7_phases.csv'

#load cell cycle cytometry data
df = pd.read_csv(cytometry_feature_csv_file_path) #feature csv
df = df.drop('Unnamed: 0', axis=1)

true_label = pd.read_csv(cytometry_phase_csv_file_path) #phase csv
true_label = list(true_label['phase'].values.flatten())
print('There are ', len(true_label), 'MCF7 cells and ', df.shape[1], 'features')
ad = sc.AnnData(df)
ad.var_names = df.columns


# true_label = [i for i in ad.obs['cell_cycle_phase']]


#normalize features
sc.pp.scale(ad)

sc.tl.pca(ad, svd_solver='arpack')
# Weight the top features (ranked by Mutual Information and Random Forest Classifier)
X_in = ad.X
df_X = pd.DataFrame(X_in)

df_X.columns = [i for i in ad.var_names]

df_X['Area'] = df_X['Area'] * 3
df_X['Dry Mass'] = df_X['Dry Mass'] * 3
df_X['Volume'] = df_X['Volume'] * 20

X_in = df_X.values
ad = sc.AnnData(df_X)

#apply PCA
sc.tl.pca(ad, svd_solver='arpack')
ad.var_names = df_X.columns

f, ax = plt.subplots(figsize = [20,10])

embedding = umap.UMAP().fit_transform(ad.obsm['X_pca'][:, 0:20])
#phate_op = phate.PHATE()
# embedding = phate_op.fit_transform(X_in)


# Plot S, M/G2, G1
cell_dict = {'T1_M1': 'yellow', 'T2_M1': 'yellowgreen', 'T1_M2': 'orange', 'T2_M2': 'darkgreen', 'T1_M3': 'red',
             'T2_M3': 'blue'}
cell_phase_dict = {'T1_M1': 'G1', 'T2_M1': 'G1', 'T1_M2': 'S', 'T2_M2': 'S', 'T1_M3': 'M/G2', 'T2_M3': 'M/G2'}

for key in list(set(true_label)):  # ['T1_M1', 'T2_M1','T1_M2', 'T2_M2','T1_M3', 'T2_M3']:
    loc = np.where(np.asarray(true_label) == key)[0]
    ax.scatter(embedding[loc, 0], embedding[loc, 1], c=cell_dict[key], alpha=.7, label=cell_phase_dict[key])
# plt.legend(markerscale=1.5, fontsize=14)
# plt.show()
true_label = [cell_phase_dict[i] for i in true_label]

#RUN VIA
knn=20
jac_std_global = 0.5
random_seed = 1
root_user = ['G1']
v0 = via.VIA(X_in, true_label, knn=knn,
         too_big_factor=0.3, root_user=root_user, dataset='group', random_seed=random_seed,
          is_coarse=True, preserve_disconnected=True, preserve_disconnected_after_pruning=True,
         pseudotime_threshold_TS=40) 
v0.run_VIA()

via.via_streamplot(via_object=v0,embedding=embedding, scatter_size=200, scatter_alpha=0.1, density_grid=.5, density_stream=1, smooth_transition=1) #you can choose to use either v1 or v0 as the input for the streamplot
plt.show()