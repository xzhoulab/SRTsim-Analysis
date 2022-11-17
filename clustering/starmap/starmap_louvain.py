##=========================================
## joint
##=========================================


import pandas as pd
import stlearn as st 
import sklearn


# nn_option = "ransam"
# inn     = "5"
isam 	= "151673"


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
	for neigh_louvain in [1,3,5,10,15,20,30,40,50,60,70,80,90,100,150,200]:
		adata_tmp = adata.copy()
		st.pp.neighbors(adata_tmp,n_neighbors=neigh_louvain) # 
		for isid in range(0,10):
			st.tl.clustering.louvain(adata_tmp,random_state=isid,key_added="louvain"+str(isid))
		adata_tmp.obs.to_csv(workdir+"/output/starmap/layer_seg/stlearn/louvain_diff_neigh/starmap_layer_srtsim_location_rpt"+str(irpt)+"_joint_louvain"+str(neigh_louvain)+".csv")


##=========================================
## separate
##=========================================


import pandas as pd
import stlearn as st 
import sklearn

# nn_option = "ransam"
# inn     = "5"
isam 	= "151673"


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
	for neigh_louvain in [1,3,5,10,15,20,30,40,50,60,70,80,90,100,150,200]:
		adata_tmp = adata.copy()
		st.pp.neighbors(adata_tmp,n_neighbors=neigh_louvain) # 
		for isid in range(0,10):
			st.tl.clustering.louvain(adata_tmp,random_state=isid,key_added="louvain"+str(isid))
		adata_tmp.obs.to_csv(workdir+"/output/starmap/layer_seg/stlearn/louvain_diff_neigh/starmap_layer_srtsim_location_rpt"+str(irpt)+"_separate_louvain"+str(neigh_louvain)+".csv")

