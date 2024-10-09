import scanpy as sc
import anndata as ad
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

#####barcodes/features/matrix文件
base_path = "D:/Project/PDAC/scRNA/GSE197177"
groups = ['PDAC1', 'PDAC2', 'PDAC3','Normal']

adatas = []
for group in groups:
    path = os.path.join(base_path, group)
    adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)  # 读取10X数据
    adata.obs['group'] = group  # 给每组添加一个分组标签
    adatas.append(adata)
    
# adata = adatas[0].concatenate(adatas[1:], batch_key='group')
adata = ad.concat(
    adatas, 
    label="group",  # 这将创建一个新的 'group' 列在 adata.obs 中
    keys=groups,  # 使用我们定义的组名作为键
    index_unique='-'  # 这会在重复的索引名后添加组名，例如 'cell1-PDAC1'
)
adata.obs_names_make_unique()

# adata.obs = adata.obs.rename(columns={'group': 'sample'})
print(adata.obs["group"].value_counts())
adata
adata.X.shape
# sc.pl.highest_expr_genes(adata, n_top=20)   




#####Doublet removal 
#!pip install scvi-tools
#！pip install --user scikit-misc
import scvi
adata
#Gene filtered
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False, 
                            flavor="seurat_v3")

#Train Model
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()
#Train SOLO Gene Model
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()
solo.predict()

df=solo.predict()
df['prediction'] = solo.predict(soft=False)
df
df['prediction'].value_counts()
df.groupby('prediction').count()
# df.index = df.index.map(lambda x: x[:-2])
# df
df['dif']=df.doublet - df.singlet
df

import seaborn as sns 
sns.displot(df[df.prediction == 'doublet'], x = 'dif')

# doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
doublets = df[(df.prediction == 'doublet')]
doublets

adata
adata.obs
adata.obs['doublet'] = adata.obs.index.isin(doublets.index)

adata = adata[~adata.obs.doublet]
adata.write('D:/Project/PDAC/scanpy/PDAC1.h5ad')

adata = sc.read('D:/Project/PDAC/scanpy/PDAC1.h5ad')


###Preprocessing
adata.var
# # sc.pl.highest_expr_genes(adata, n_top=20)  
# #去除线粒体基因
# adata.var['mt'] = adata.var.index.str.startswith('MT-')
# adata.var

# #去除核糖体基因

# ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
# ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
# adata.var['ribo']= adata.var_names.isin(ribo_genes[0].values)

# sc.pp.calculate_qc_metrics(
    # adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True
    # )

# adata.var

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var['mt'].value_counts()

# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var['ribo'].value_counts()

# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
adata.var['hb'].value_counts()

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False
)
adata.var
adata.obs
##过滤基因
# adata.var.sort_values('n_cells_by_counts')
# sc.pp.filter_genes(adata, min_cells=3)
adata.var.sort_values('n_cells_by_counts')

# sc.pl.violin(
#     adata,
#     ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
#     jitter=0.4,
#     multi_panel=True,
# )

# sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

sc.pl.violin(
    adata, 
    ['n_genes_by_counts', 'total_counts', 
     'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'], 
    jitter=0.4, 
    multi_panel=True)
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")


upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
upper_lim
adata = adata[adata.obs.n_genes_by_counts < upper_lim]

#去除线粒体/核糖体/血红蛋白
adata.obs['group'].value_counts()
adata = adata[adata.obs.pct_counts_mt < 20]
adata = adata[adata.obs.pct_counts_ribo < 20]
adata = adata[adata.obs.pct_counts_hb < 5]

adata.obs['group'].value_counts()

adata.write('D:/Project/PDAC/scanpy/PDAC2.h5ad')

adata = sc.read('D:/Project/PDAC/scanpy/PDAC2.h5ad')

#####Normalization
adata.X.sum(axis=1)
sc.pp.normalize_total(adata, target_sum=1e4)
# adata.X.sum(axis=1)
sc.pp.log1p(adata)
adata.raw = adata

# 去周期
adata.obs['group'].value_counts()
s_genes = pd.read_csv("D:\Project\scRNA\PDCA2\s_genes.csv", header=None)[0].tolist()
g2m_genes = pd.read_csv("D:\Project\scRNA\PDCA2\g2m_genes.csv", header=None)[0].tolist()
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
adata.obs['phase'].value_counts()
sc.pp.regress_out(adata, keys=['S_score', 'G2M_score'])
adata.obs['group'].value_counts()

###Clustering
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.var

sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 
                          'pct_counts_mt', 'pct_counts_ribo','pct_counts_hb'])
sc.pp.scale(adata,max_value = 10)
sc.tl.pca(adata, svd_solver = "arpack")
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

sc.pp.neighbors(adata, n_pcs=32)
# adata.obsp['connectivities']

sc.tl.umap(adata)
sc.pl.umap(adata)

sc.tl.tsne(adata)
sc.pl.tsne(adata)

# !pip install leidenalg
sc.tl.leiden(adata, resolution=1)
sc.pl.umap(adata, color=['leiden','group'], 
           legend_loc="on data",
           frameon=True,)
adata.obs['leiden'].value_counts()

adata.write('D:/Project/PDAC/scanpy/PDAC3.h5ad')

adata = sc.read('D:/Project/PDAC/scanpy/PDAC3.h5ad')
#Findmakers/Label cell types
sc.tl.rank_genes_groups(adata, 'leiden')

# sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
makers = sc.get.rank_genes_groups_df(adata, None)
makers = makers[(makers.pvals_adj <0.05) & (makers.logfoldchanges > .5)]
makers
makers.to_csv('D:/Project/scRNA/scanpy/CSV/makers.csv', index=False)
# makers = pd.read_csv("D:/Project/scRNA/scanpy/CSV/makers.csv")
# markers_scvi = model.differential_expression(groupby = 'leiden')
# markers_scvi

########Method2
# sc.pl.umap(adata, 
#            color = ['leiden'], 
#            frameon = False, 
#            legend_loc = "on data")

# marker_genes = {
#     'Ductal Cells':["ELF3", "FXYD2", "LCN2", "CEACAM6", 
#                   "FXYD3", "EPCAM", "MUC1", "CEACAM5", 
#                   "KRT19", "AGR2", "TFF2", "TSPAN8", "MMP7"], 
#     'T Cells':["FOS", "JUN", 
#               "IL7R", "IL2RA", "TRAC", "GZMA", 
#               "GZMB", "CD3D", "CD3E", "CCL5", 
#               "SLC30A1", "GZMH"], 
#     'NK Cells':[ "NKG7", "KLRD1"], 
#     'B Cells':["BANK1", "LY9", "TNFRSF13C", "MS4A1"], 
#     'Plasama Cells':["FKBP11","MZB1", "SSR4", "IGHA1"], 
#     'Mast Cells':["CD69", "MS4A2", "TPSB2", "CLU"],
#     'Mac/Mono Cells':["S100A8", "S100A9", "CCL2", "CXCL1", 
# 				"AIF1", "CD74", "HLA-DRA", "CXCL3", 
# 				"CXCL2", "APOE", "CD68", "FCER1G"], 
#     'Endothelial Cells':["PLVAP", "PECAM1"], 
#     'Fibroblasts':["COL3A1", "ACTA2", "LUM", "COL1A1", "DCN", "COL1A2"], 
#     "Endocrine Cells":["PRSS1", "INS", "SST", "CELA2A", 
#                      "LYZ","CPA1", "CELA3B", 'CELA3A', 
#                      "ELF3"]
# }

# sc.pl.dotplot(adata, marker_genes, groupby="leiden", 
#               standard_scale="var", dendrogram=False,
#               show = False)
# plt.savefig(r"D:/Project/scRNA/scanpy/plot/dotplot.tiff", dpi=300)

# adata.obs["cell_type"] = adata.obs["leiden"].map(
#        {
#         "0": "Fibroblasts",
#         "1": "Mac/Mono Cells",
#         "2": "Ductal Cells",
#         "3": "T Cells",
#         "4": "T Cells",
#         "5": "Ductal Cells",
#         "6": "Endocrine Cells",
#         "7": "Ductal Cells",
#         "8": "Ductal Cells",
#         "9": "T Cells",
#         "10": "T Cells",
#         "11": "Ductal Cells",
#         "12": "Plasama Cells",
#         "13": "T Cells",
#         "14": "Endocrine Cells",
#         "15": "Mac/Mono Cells",
#         "16": "Endocrine Cells",
#         "17": "Endocrine Cells",
#         "18": "Fibroblasts",
#         "19": "19",
#         "20": "Mast Cells", 
#         "21":"Ductal Cells", 
#         "22":"Endothelial Cells", 
#         "23":"Fibroblasts" , 
#         "24":"Mac/Mono Cells"            
#     }
# )


# sc.pl.umap(adata, color=['cell_type', 'group'], 
#            legend_loc="on data",
#            frameon=True,)


# adata1 = adata[adata.obs.leiden.isin(["19"])]

# markers_other = {
#     "Macrophages": ["CD68", "CD163", "AIF1", "CSF1R", "MARCO", "CD14"],
#     "Neutrophils": ['S100A9', 'S100A8', 'MMP19', 'ELANE', 'CSTA'],
#     "Complement": ['C1QA', 'C1QB', 'C3'],
#     "Monocytes": ['CD68', 'CD163', 'CD14', 'FCGR3A', 'ITGAM']
# }

# sc.pl.dotplot(adata1, var_names=markers_other, 
#               groupby="leiden",standard_scale="var", 
#               swap_axes=False, dendrogram=False)

# adata.obs["cell_type"] = adata.obs["leiden"].map(
#        {
#         "0": "Ductal Cells",
#         "1": "Ductal Cells",
#         "2": "T Cells",
#         "3": "Endocrine Cells",
#         "4": "Fibroblasts",
#         "5": "Mac/Mono Cells",
#         "6": "Ductal Cells",
#         "7": "T Cells",
#         "8": "Ductal Cells",
#         "9": "T Cells",
#         "10": "Goblet Cells",
#         "11": "T Cells",
#         "12": "Endocrine Cells",
#         "13": "Endocrine Cells",
#         "14": "Mast Cells",
#         "15": "Ductal Cells",
#         "16": "Complement",
#         "17": "Plasama Cells",
#         "18": "Fibroblasts",
#         "19": "B Cells",
#         "20": "Endothelial Cells", 
#         "21":"Fibroblasts", 
#         "22":"Mac/Mono Cells", 
#         "23":"Ductal Cells"             
#     }
# )

########Method2
adata.obs['leiden'].value_counts()
sc.pl.umap(adata, color=["leiden"],legend_loc="on data",
           add_outline=False,      #cluster描边
           legend_fontsize=5,    #字体大小
           legend_fontoutline=1, 
           frameon=True, 
           show = False, 
           title = "leiden")  #字体描边大小
plt.savefig(r"D:/Project/scRNA/scanpy/plot/leiden_plot.tiff", dpi=300)

###免疫细胞鉴定
#Base on: <https://www.cellsignal.com/pathways/immune-cell-markers-human>
#1.分辨髓系/淋巴细胞
#Leukocytes
sc.pl.umap(adata, color=['PTPRC',"CD3E","ITGAM"], frameon=False)
print('PTPRC' in adata.var_names)
makers[makers.names=='PTPRC']
makers[makers.names=='CD3E']
makers[makers.names=='ITGAM']

#鉴定 CD8/CD4
sc.pl.umap(adata, color=['CD8A',"CD4"], frameon=False)
makers[makers.names=='CD8A']
makers[makers.names=='CD4']

#鉴定 B Cells
sc.pl.umap(adata, color=['PTPRC',"PXK","MS4A1"], frameon=False)
sc.pl.umap(adata, color=['PTPRC',"BANK","CD79A"], frameon=False)

#鉴定 Plasma Cells
sc.pl.umap(adata, color=['PTPRC',"MZB1","IGHG1"], frameon=False)
sc.pl.umap(adata, color=['PTPRC',"JCHAIN","SPAG4"], frameon=False)
makers[makers.names=='MZB1']
makers[makers.names=='IGHG1']
makers[makers.names=='JCHAIN']
makers[makers.names=='SPAG4']

#鉴定 NK Cells
sc.pl.umap(adata, color=['PTPRC',"TRDC","NKG7"], frameon=False)
sc.pl.umap(adata, color=['PTPRC',"KLRF1","KLRD1"], frameon=False)
makers[makers.names=='TRDC']
makers[makers.names=='NKG7']
makers[makers.names=='KLRF1']
makers[makers.names=='KLRD1']

#鉴定 NKT Cells
sc.pl.umap(adata, color=['PTPRC',"NCAM1","GATA3"], frameon=False)

#鉴定 巨噬细胞Macrophages
sc.pl.umap(adata, color=['PTPRC',"CD68","NAAA"], frameon=False)
makers[makers.names=='CD68']
makers[makers.names=='NAAA']

#鉴定 经典树突状细胞Conventional Dendritic Cells(cDCs)
sc.pl.umap(adata, color=['PTPRC',"ITGAX","ZBTB46"], frameon=False)
makers[makers.names=='ITGAX']
makers[makers.names=='ZBTB46']

#鉴定 浆细胞样树突状细胞Plasmacytoid Dendritic Cells(pDCs)
sc.pl.umap(adata, color=['PTPRC',"BST2","IL3RA"], frameon=False)
makers[makers.names=='BST2']
makers[makers.names=='IL3RA']

#鉴定 髓系抑制性细胞Myeloid-Derived Suppressor Cells(MDSCs)
sc.pl.umap(adata, color=['PTPRC',"S100A4","CD14"], frameon=False)
makers[makers.names=='S100A4']
makers[makers.names=='CD14']


#鉴定 上皮组织
Ductal_makers = ["ELF3", "FXYD2", "LCN2", "CEACAM6", 
                  "FXYD3", "EPCAM", "MUC1", "CEACAM5", 
                  "KRT19", "AGR2", "TFF2", "TSPAN8", "MMP7"]

sc.pl.umap(adata, color=["ELF3", "FXYD2", "LCN2", "CEACAM6"], frameon=False)
makers[makers.names=='ELF3']
makers[makers.names=='FXYD2']
makers[makers.names=='LCN2']
makers[makers.names=='CEACAM6']



adata.obs["cell_type"] = adata.obs["leiden"].map(
       {
        "0": "",
        "1": "CD4+ T Cells/Macrophages/cDCs/pDCs",
        "2": "",
        "3": "immune",
        "4": "CD8+ T Cells",
        "5": "",
        "6": "",
        "7": "",
        "8": "",
        "9": "immune",
        "10": "CD8+ T Cells",
        "11": "",
        "12": "Plasma Cells",
        "13": "CD4+ T Cells",
        "14": "",
        "15": "",
        "16": "",
        "17": "",
        "18": "",
        "19": "MDSCs",
        "20": "", 
        "21":"", 
        "22":"", 
        "23":"" , 
        "24":"CD4+ T Cells"            
    }
)









sc.pl.umap(adata, color=['PTPRC',"CD3E","CD4"], frameon=False)
sc.pl.umap(adata, color=['PTPRC',"CD3E","CD8A"], frameon=False)
makers[makers.names=='PTPRC']
makers[makers.names=='CD4']
makers[makers.names=='CD8A']

adata.obs["cell_type"] = adata.obs["leiden"].map(
       {
        "0": "Fibroblasts",
        "1": "Mac/Mono Cells",
        "2": "Ductal Cells",
        "3": "T Cells",
        "4": "T Cells",
        "5": "Ductal Cells",
        "6": "Endocrine Cells",
        "7": "Ductal Cells",
        "8": "Ductal Cells",
        "9": "T Cells",
        "10": "T Cells",
        "11": "Ductal Cells",
        "12": "Plasama Cells",
        "13": "T Cells",
        "14": "Endocrine Cells",
        "15": "Mac/Mono Cells",
        "16": "Endocrine Cells",
        "17": "Endocrine Cells",
        "18": "Fibroblasts",
        "19": "19",
        "20": "Mast Cells", 
        "21":"Ductal Cells", 
        "22":"Endothelial Cells", 
        "23":"Fibroblasts" , 
        "24":"Mac/Mono Cells"            
    }
)


sc.pl.umap(adata, color=["cell_type"],legend_loc="on data",
           add_outline=False,      #cluster描边
           legend_fontsize=5,    #字体大小
           legend_fontoutline=1, 
           frameon=True, 
           show = False, 
           title = "Cell Types")  #字体描边大小
plt.savefig(r"D:/Project/scRNA/scanpy/plot/upam_plot.tiff", dpi=300)

sc.pl.tsne(adata, color=["cell_type"],legend_loc="on data",
           add_outline=False,      #cluster描边
           legend_fontsize=5,    #字体大小
           legend_fontoutline=1, 
           frameon=True, 
           show = False, 
           title = "Cell Types")  #字体描边大小
plt.savefig(r"D:/Project/scRNA/scanpy/plot/tsne_plot.tiff", dpi=300)


adata.write('D:/Project/PDAC/scanpy/PDAC4.h5ad')

adata = sc.read('D:/Project/PDAC/scanpy/PDAC4.h5ad')
###
adata.obs['cell_type'].value_counts()
adata.obs['sample'].unique().tolist()


def map_condition(x): 
    if 'PDAC' in x : 
        return 'PDAC'
    else:
        return "Normal"

adata.obs['condition'] = adata.obs['sample'].map(map_condition)
num_tot_cells = adata.obs.groupby(['sample']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
num_tot_cells

cell_type_counts = adata.obs.groupby(['sample','condition','cell_type']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
cell_type_counts
cell_type_counts = cell_type_counts.rename(columns={'_scvi_batch': 'counts'})
cell_type_counts['total_cells'] = cell_type_counts['sample'].map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.counts / cell_type_counts.total_cells
cell_type_counts


import matplotlib.pyplot as plt
import seaborn as sns 
plt.figure(figsize = (10,4))
ax = sns.boxplot(data = cell_type_counts, x = 'cell_type', y = 'frequency', hue = 'condition')
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.show()

adata.write('D:/Project/PDAC/scanpy/PDAC5.h5ad')

adata = sc.read('D:/Project/PDAC/scanpy/PDAC5.h5ad')

#Re Clustering
# t_cells = adata[adata.obs['cell_type'].isin(['T Cells'])].copy()

# sc.pp.normalize_total(t_cells, target_sum=1e4)
# sc.pp.log1p(t_cells)
# sc.pp.pca(t_cells, n_comps=50)
# sc.pp.neighbors(t_cells, n_neighbors=15, n_pcs=50)
# sc.tl.leiden(t_cells, resolution=0.5)  
# t_cells.obs['leiden'].value_counts()
# sc.pl.umap(t_cells, color='leiden', legend_loc='on data')

#DE 
subset = adata[adata.obs["cell_type"].isin(['T Cells','B Cells'])].copy()

import diffxpy.api as de
subset.X = subset.X.toarray()

# import scipy.sparse as sp
# #!pip install diffxpy
# if sp.issparse(subset.X):
#     subset.X = subset.X.toarray()
# # subset.X = subset.X.toarray()

res = de.test.wald(
    data=subset, 
    formula_loc="~1 + cell_type", 
    factor_loc_totest="cell_type"
)
dedf = res.summary()
dedf


# Pseudotemporal ordering
adata.write('D:/Project/PDAC/scanpy/PDAC4.h5ad')

adata = sc.read('D:/Project/PDAC/scanpy/PDAC4.h5ad')
adata

sc.pl.scatter(adata, basis="tsne", color="")

sc.tl.diffmap(adata)
root_ixs = adata.obsm["X_diffmap"][:, 3].argmin()
sc.pl.scatter(
    adata,
    basis="diffmap",
    color=["cell_type"],
    components=[2, 3],
)

adata.uns["iroot"] = root_ixs

sc.tl.dpt(adata)
sc.pl.scatter(
    adata,
    basis="tsne",
    color=["dpt_pseudotime", "palantir_pseudotime"],
    color_map="gnuplot2",
)
sc.pl.violin(
    adata,
    keys=["dpt_pseudotime", "palantir_pseudotime"],
    groupby="clusters",
    rotation=45,
    order=[
        "HSC_1",
        "HSC_2",
        "Precursors",
        "Ery_1",
        "Ery_2",
        "Mono_1",
        "Mono_2",
        "CLP",
        "DCs",
        "Mega",
    ],
)

adata = sc.read('D:/Project/PDAC/scanpy/PDAC4.h5ad')
adata
# RNA velocity
import scvelo as scv 
#pip install scvelo
scv.settings.set_figure_params("scvelo")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

scv.pl.scatter(adata, basis="umap", color="cell_type")

