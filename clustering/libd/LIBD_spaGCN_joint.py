## Thu Mar 17 07:51:08 2022
### prepare the data 
##===================
import anndata
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt

isam 	= "151673"


for irpt in range(1,11):
	spatial = pd.read_csv(workdir+"/data/simLIBD/newsimLIBD/LIBD_"+str(isam)+"_srtsim_location_rpt"+str(irpt)+"_joint.csv",index_col=0)
	adata = anndata.read_csv(workdir+"/data/simLIBD/newsimLIBD/LIBD_"+str(isam)+"_srtsim_cmat_rpt"+str(irpt)+"_joint.csv",first_column_names=0)
	adata.obs["x_array"]=spatial["row"]
	adata.obs["y_array"]=spatial["col"]
	adata.obs["x_pixel"]=spatial["imagerow"]
	adata.obs["y_pixel"]=spatial["imagecol"]
	adata.obs["label"]= pd.Categorical(spatial["label"])
	adata.var_names=[i.upper() for i in list(adata.var_names)]
	adata.var["genename"]=adata.var.index.astype("str")
	adata.write_h5ad(workdir+"/output/simLIBD/newLIBD/spaGCN/LIBD_"+str(isam)+"_srtsim_anndata_rpt"+str(irpt)+"_joint.h5ad")




## Thu Mar 17 07:51:03 2022
##===================
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt

import sklearn

isam 	= "151673"


for irpt in range(1,11):
	adata=sc.read(workdir+"/output/simLIBD/newLIBD/spaGCN/LIBD_"+str(isam)+"_srtsim_anndata_rpt"+str(irpt)+"_joint.h5ad")
	x_pixel=adata.obs["x_pixel"].tolist()
	y_pixel=adata.obs["y_pixel"].tolist()
	adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
	## for refined result
	x_array=adata.obs["x_array"].tolist()
	y_array=adata.obs["y_array"].tolist()
	adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
	adata.var_names_make_unique()
	spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
	spg.prefilter_specialgenes(adata)
	#Normalize and take log for UMI
	sc.pp.normalize_per_cell(adata)
	sc.pp.log1p(adata)
	#p: percentage of total expression contributed by neighborhoods.
	p = 0.5 
	#l: parameter to control p.
	#Find the l value given p, first use spg.test_l() function to get a rough estimate of the range l falls in
	spg.test_l(adj,[1, 10, 100, 500, 1000,2000,4000])
	#Search l from 100 to 500
	l=spg.find_l(p=p,adj=adj,start=1, end=10,sep=0.1, tol=0.01)
	n_clusters=pd.Categorical(adata.obs["label"]).categories.size
	#Set seed
	pred_score = [0]*10
	refined_pred_score = [0]*10
	for isid in range(0,10):
		r_seed=t_seed=n_seed=isid
		#Search for suitable resolution
		res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
		clf=spg.SpaGCN()
		clf.set_l(l)
		#Set seed
		random.seed(r_seed)
		torch.manual_seed(t_seed)
		np.random.seed(n_seed)
		#Run
		clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
		y_pred, prob=clf.predict()
		adata.obs["pred"]= y_pred
		adata.obs["pred"]=adata.obs["pred"].astype('category')
		refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
		adata.obs["refined_pred"]=refined_pred
		adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
		adata.obs[["pred","refined_pred"]].to_csv(workdir+"/output/simLIBD/newLIBD/spaGCN/LIBD_"+str(isam)+"_srtsim_rpt"+str(irpt)+"_seed"+str(isid)+"_joint_spaGCN.csv")
		

