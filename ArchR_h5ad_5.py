import ArchR_h5ad
import os
import pandas as pd
import numpy as np
import scanpy as sc


arrow_path1 = "./HemeTutorial/ArrowFiles/GB2_ATAC_S1_2.arrow"
arrow_path2 = "./HemeTutorial/ArrowFiles/GB3_ATAC_S3_5.arrow"
arrow_path3 = "./HemeTutorial/ArrowFiles/GB4_ATAC_S6_7.arrow"


adata1 = ArchR_h5ad.read_arrow(arrow_path1, use_matrix="GeneScoreMatrix")
adata2 = ArchR_h5ad.read_arrow(arrow_path2, use_matrix="GeneScoreMatrix")
adata3 = ArchR_h5ad.read_arrow(arrow_path3, use_matrix="GeneScoreMatrix")


adata1.obs['batch']='GB2_ATAC_S1_2'
adata2.obs['batch']='GB3_ATAC_S3_5'
adata3.obs['batch']='GB4_ATAC_S6_7'


adata = sc.AnnData.concatenate(adata1,adata2,adata3,batch_key = 'Batch')

test=pd.DataFrame(data=adata.obs)
test = test.loc[:,['batch','CellNames']]
test['CellNames2'] = ['#'.join(i) for i in test.values]
adata.obs['CellNames2']=pd.Categorical(test.CellNames2)

meta = pd.read_table('./meta.tsv')
adata_sub = adata[adata.obs.CellNames2.isin(meta.Cell), :]

meta['num'] = range(len(meta))
meta.set_index(['Cell'],inplace=True)
test2=pd.DataFrame(data=adata_sub.obs)
meta2 = meta.loc[test2.CellNames2]

adata_sub.obs['annotation2']=pd.Categorical(meta2.annotation2)

x_umap= np.loadtxt("./umap_embedding.tsv")
adata_sub.obsm['X_umap']=x_umap[meta2.num]

adata_sub.write_h5ad('./HRO-5.h5ad')
