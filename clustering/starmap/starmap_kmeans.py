

## Thu Apr 21 04:13:52 2022

## joint
##=========================================

import pandas as pd
import stlearn as st 
import sklearn


for irpt in range(1,11):
	count_matrix = pd.read_csv(workdir+"/data/starmap/layer_seg/starmap_layer_srtsim_cmat_rpt"+str(irpt)+"_joint.csv",index_col=0)
	spatial = pd.read_csv(workdir+"/data/starmap/layer_seg/starmap_layer_srtsim_location_rpt"+str(irpt)+"_joint.csv",index_col=0)
	new_spatial = spatial.iloc[:,0:2]
	new_spatial.columns = ['imagerow','imagecol']
	adata = st.create_stlearn(count=count_matrix,spatial=new_spatial,library_id="Sample_test", scale=1,background_color="white")
	numcluster = pd.Categorical(spatial["label"]).categories.size
	adata.obs["label"] = pd.Categorical(spatial["label"])
	st.pp.filter_genes(adata,min_cells=1)
	st.pp.normalize_total(adata)
	st.pp.log1p(adata)
	st.pp.scale(adata)
	st.em.run_pca(adata,n_comps=50,random_state=0)
	for isid in range(0,10):
		st.tl.clustering.kmeans(adata,n_clusters=numcluster,random_state=isid,key_added="kmeans"+str(isid))
	adata.obs.to_csv(workdir+"/output/starmap/layer_seg/stlearn/starmap_layer_srtsim_rpt"+str(irpt)+"_joint_kmeans.csv")




##separate
##=========================================

import pandas as pd
import stlearn as st 
import sklearn

for irpt in range(1,11):
	count_matrix = pd.read_csv(workdir+"/data/starmap/layer_seg/starmap_layer_srtsim_cmat_rpt"+str(irpt)+"_separate.csv",index_col=0)
	spatial = pd.read_csv(workdir+"/data/starmap/layer_seg/starmap_layer_srtsim_location_rpt"+str(irpt)+"_separate.csv",index_col=0)
	new_spatial = spatial.iloc[:,0:2]
	new_spatial.columns = ['imagerow','imagecol']
	adata = st.create_stlearn(count=count_matrix,spatial=new_spatial,library_id="Sample_test", scale=1,background_color="white")
	numcluster = pd.Categorical(spatial["label"]).categories.size
	adata.obs["label"] = pd.Categorical(spatial["label"])
	st.pp.filter_genes(adata,min_cells=1)
	st.pp.normalize_total(adata)
	st.pp.log1p(adata)
	st.pp.scale(adata)
	st.em.run_pca(adata,n_comps=50,random_state=0)
	for isid in range(0,10):
		st.tl.clustering.kmeans(adata,n_clusters=numcluster,random_state=isid,key_added="kmeans"+str(isid))
	adata.obs.to_csv(workdir+"/output/starmap/layer_seg/stlearn/starmap_layer_srtsim_rpt"+str(irpt)+"_separate_kmeans.csv")
 

