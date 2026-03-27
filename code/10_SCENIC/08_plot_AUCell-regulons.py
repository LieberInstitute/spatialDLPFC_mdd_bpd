import sys, pickle
from ctxcore.genesig import openfile

import loompy as lp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


pdir = "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/"
MODULE_PATH = pdir+"processed-data/10_SCENIC/AUCell_regulons/"
RESULTS_PATH = pdir+"processed-data/10_SCENIC/"
PLOT_PATH = pdir+"plots/10_SCENIC/"

MODULE_TYPE = "regulons-top20"


def extract_aucell(cluster, nGenes, modType):
    lf = lp.connect(MODULE_PATH+"spe-n109_seurat-label-"+cluster+"_"+nGenes+"-genes_no-lowUMI_logcounts_"+modType+"_AUCell-fixed-thold-05.loom", mode='r+', validate=False)
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    cell_mtx = pd.DataFrame({"seurat_label": lf.ca.seurat_label, "smoothed_k9_1663": lf.ca.smoothed_k9_1663, #"custom_cluster": cluster,
                             "sex": lf.ca.sex, "condition": lf.ca.condition, "nGene": lf.ca.nGene, "nUMI": lf.ca.nUMI}, index=lf.ca.CellID)
    lf.close()
    return auc_mtx, cell_mtx


auc_mv, cell_mv = extract_aucell("Micro.Vasc", "13162", MODULE_TYPE)
auc_ast, cell_ast = extract_aucell("Astro", "13162", MODULE_TYPE)
auc_l23, cell_l23 = extract_aucell("L2.3", "13162", MODULE_TYPE)
auc_l4, cell_l4 = extract_aucell("L4", "13162", MODULE_TYPE)
auc_inh, cell_inh = extract_aucell("Inhb", "13162", MODULE_TYPE)
auc_l5, cell_l5 = extract_aucell("L5", "13162", MODULE_TYPE)
auc_l6, cell_l6 = extract_aucell("L6", "13162", MODULE_TYPE)
auc_wm, cell_wm = extract_aucell("Oligo", "13162", MODULE_TYPE)


auc_all = pd.concat([auc_mv, auc_ast, auc_l23, auc_l4, auc_inh, auc_l5, auc_l6, auc_wm])
cell_all = pd.concat([cell_mv, cell_ast, cell_l23, cell_l4, cell_inh, cell_l5, cell_l6, cell_wm])

joined_all = pd.concat([auc_all, cell_all], axis=1)
joined_all.to_csv(RESULTS_PATH+"spe-n109_13162-no-lowUMI_"+MODULE_TYPE+"_AUCell.csv")

#precast domains only

joined_precast = joined_all.loc[-joined_all['smoothed_k9_1663'].isin(['Vasc','GABA'])]

#output regulons as merged df
with open(RESULTS_PATH+'spe-n109_13162-no-lowUMI_logcounts.'+MODULE_TYPE+'.dat', 'rb') as file:
    regulons = pickle.load(file)

reg_df = pd.DataFrame({"TF": [item.transcription_factor for item in regulons],
                       #"dir": ["promoter" if "activating" in item.context else "enhancer" for item in regulons],
                       "dir": ["activating" if "activating" in item.context else "repressing" for item in regulons],
                       "set_size": [len(item) for item in regulons],
                      "set_str": ['/' .join(list(item.genes)) for item in regulons]})

reg_df.to_csv(RESULTS_PATH+"spe-n109_13162-no-lowUMI_logcounts_"+MODULE_TYPE+"_merged.csv")



#color palette for plotting

transfer_bright = {'Micro.Vasc': "#911223", 'Astro': "#cfa45c", 'L2.3': "#088F8F", 'L4': "#c2cfcf",
                   'Inhb': "#9377AC", 'L5': "#ddc94e", 'L6': "#E45C5F", 'Oligo': "#D1C4B0"}
smoothed_bright = {'Vasc': "#911223", 'L1': "#cfa45c", 'L2': "#5D9940", 'L3.4': "#5095CD",
                  'GABA': "#9377AC", 'L5': "#ddc94e", 'L6': "#E45C5F", 'WM': "#D1C4B0"}


# plot and save Seurat heatmap

df_scores = joined_all[list(auc_all.select_dtypes('number').columns)+['seurat_label','sex','condition']]
df_means = df_scores.groupby(by=['seurat_label','sex','condition']).mean()

colors_df = pd.DataFrame(index=df_means.index)
colors_df['seurat_label'] = [x[0] for x in df_means.index]
row_colors = colors_df['seurat_label'].map(transfer_bright)

hmp = sns.clustermap(df_means.T, z_score='row', col_colors=row_colors, dendrogram_ratio=.1,
               cmap='bwr', vmin=-3.5, vmax=3.5,
              cbar_pos=(0.02, 0.8, 0.03, 0.15))
hmp.ax_heatmap.set_yticklabels(hmp.ax_heatmap.get_ymajorticklabels(), fontsize = 10)

hmp.figure.savefig(PLOT_PATH+'spe-n109_13162-no-lowUMI_logcounts_'+MODULE_TYPE+'_AUCell-seurat-pc30_heatmap.pdf')



# plot and save PRECAST heatmap

df_scores = joined_all[list(auc_all.select_dtypes('number').columns)+['smoothed_k9_1663','sex','condition']]
df_scores = df_scores.loc[-df_scores['smoothed_k9_1663'].isin(['Vasc','GABA'])]
df_means = df_scores.groupby(by=['smoothed_k9_1663','sex','condition']).mean()

colors_df = pd.DataFrame(index=df_means.index)
colors_df['smoothed_bright'] = [x[0] for x in df_means.index]
row_colors = colors_df['smoothed_bright'].map(smoothed_bright)

hmp = sns.clustermap(df_means.T, z_score='row', col_colors=row_colors, dendrogram_ratio=.1,
               cmap='bwr', vmin=-3.5, vmax=3.5,
              cbar_pos=(0.02, 0.8, 0.03, 0.15))
hmp.ax_heatmap.set_yticklabels(hmp.ax_heatmap.get_ymajorticklabels(), fontsize = 10)

hmp.figure.savefig(PLOT_PATH+'spe-n109_13162-no-lowUMI_logcounts_'+MODULE_TYPE+'_AUCell-smoothed-k9-1663_heatmap.pdf')


# correlate AUCell

def correlate_aucell(cluster, auc_mtx):
    auc_corr = np.corrcoef(auc_mtx.T)
    corrDF = pd.DataFrame(auc_corr, index=auc_mtx.columns, columns= auc_mtx.columns)
    corrDF.to_csv(RESULTS_PATH+"spe-n109_13162-no-lowUMI_"+MODULE_TYPE+"_AUCell-"+cluster+"-correlation.csv")
    return corrDF

corrDF = correlate_aucell("all-spots", auc_all)
#correlate_aucell("Micro.Vasc", auc_mv)
#correlate_aucell("Astro", auc_ast)
#correlate_aucell("L2.3", auc_l23)
#correlate_aucell("L4", auc_l4)
#correlate_aucell("Inhb", auc_inh)
#correlate_aucell("L5", auc_l5)
#correlate_aucell("L6", auc_l6)
#correlate_aucell("Oligo", auc_wm)


hmp = sns.clustermap(corrDF, dendrogram_ratio=.1,
                    cbar_pos=(0.02, 0.8, 0.03, 0.15))
hmp.ax_heatmap.set_yticklabels(hmp.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
hmp.figure.savefig(PLOT_PATH+'spe-n109_13162-no-lowUMI_logcounts_'+MODULE_TYPE+'_AUCell-correlation.pdf')

#sys.exit("Early stopping to check what regulons I have...")

# plot ECDF curves showing differences in seurat label/ domain specificity

fig, axs = plt.subplots(6, 2, figsize=(6, 12))

axs = axs.flatten()

sns.ecdfplot(data=joined_all, x="KLF4(+)", hue="seurat_label", palette=transfer_bright, 
             legend=False, ax=axs[0])
axs[0].set_title("KLF4 (Seurat)")
sns.ecdfplot(data=joined_precast, x="KLF4(+)", hue="smoothed_k9_1663", palette=smoothed_bright, 
             legend=False, ax=axs[1])
axs[1].set_title("KLF4 (PRECAST)")

sns.ecdfplot(data=joined_all, x="CEBPD(+)", hue="seurat_label", palette=transfer_bright, 
             legend=False, ax=axs[2])
axs[2].set_title("CEBPD (Seurat)")
sns.ecdfplot(data=joined_precast, x="CEBPD(+)", hue="smoothed_k9_1663", palette=smoothed_bright, 
             legend=False, ax=axs[3])
axs[3].set_title("CEBPD (PRECAST)")

sns.ecdfplot(data=joined_all, x="SOX2(+)", hue="seurat_label", palette=transfer_bright, 
             legend=False, ax=axs[4])
axs[4].set_title("SOX2 (Seurat)")
sns.ecdfplot(data=joined_precast, x="SOX2(+)", hue="smoothed_k9_1663", palette=smoothed_bright, 
             legend=False, ax=axs[5])
axs[5].set_title("SOX2 (PRECAST)")

sns.ecdfplot(data=joined_all, x="LHX6(+)", hue="seurat_label", palette=transfer_bright, 
             legend=False, ax=axs[6])
axs[6].set_title("LHX6 (Seurat)")
sns.ecdfplot(data=joined_precast, x="LHX6(+)", hue="smoothed_k9_1663", palette=smoothed_bright, 
             legend=False, ax=axs[7])
axs[7].set_title("LHX6 (PRECAST)")

sns.ecdfplot(data=joined_all, x="CUX2(+)", hue="seurat_label", palette=transfer_bright, 
             legend=False, ax=axs[8])
axs[8].set_title("CUX2 (Seurat)")
sns.ecdfplot(data=joined_precast, x="CUX2(+)", hue="smoothed_k9_1663", palette=smoothed_bright, 
             legend=False, ax=axs[9])
axs[9].set_title("CUX2 (PRECAST)")

sns.ecdfplot(data=joined_all, x="THRA(+)", hue="seurat_label", palette=transfer_bright, 
             legend=False, ax=axs[10])
axs[10].set_title("THRA (Seurat)")
sns.ecdfplot(data=joined_precast, x="THRA(+)", hue="smoothed_k9_1663", palette=smoothed_bright, 
             legend=False, ax=axs[11])
axs[11].set_title("THRA (PRECAST)")

plt.tight_layout()
plt.savefig(PLOT_PATH+'spe-n109_13162-no-lowUMI_'+MODULE_TYPE+'_AUCell-layers-ecdf.pdf', 
            bbox_inches='tight')


# plot density curves (helpful for norm)

fig, axs = plt.subplots(5, 5, figsize=(8, 8))

axs_flat = axs.flatten()

for i in range(23):
    gene1 = auc_all.columns[i]
    sns.kdeplot(data=joined_all, x=gene1, ax=axs_flat[i])

plt.tight_layout()
plt.savefig(PLOT_PATH+'spe-n109_13162-no-lowUMI_'+MODULE_TYPE+'_AUCell-density.pdf', 
            bbox_inches='tight')
