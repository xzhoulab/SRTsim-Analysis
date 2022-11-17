## separate 
##============================
rm(list=ls())

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)



resol_list <- 1:30/20

for(irpt in 1:10){
	# irpt = 1
	load(paste0(workdir,"/data/starmap/layer_seg/starmap_layer_srtsim_simdata_rpt",irpt,"_separate.rds"))

	# colnames(simCount_mat) <- paste0("cell",1:ncol(simCount_mat))
	scobj <- CreateSeuratObject(counts = simCount_mat, project ="starmap", min.cells = 3, min.features = 100)
	scobj <- SCTransform(scobj, assay = "RNA", verbose = FALSE)
	scobj <- RunPCA(scobj, assay = "SCT", verbose = FALSE)
	scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:30)

	for(isid in 1:10){
		# isid = 1
		scd_loop <- scobj
		for(iresol in resol_list){
			scd_loop <- FindClusters(scd_loop, verbose = FALSE,resolution=iresol,random.seed=isid)
		}
		write.table(scd_loop@meta.data,file=paste0(workdir,"/output/starmap/layer_seg/seurat/starmap_layer_srtsim_simdata_rpt",irpt,"_seed",isid,"_separate_seurat.txt"))
		rm(scd_loop)
	}	
	rm(scobj,simCount_mat)
}




## joint 
##============================
rm(list=ls())

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


resol_list <- 1:30/20

for(irpt in 1:10){
	# irpt = 1
	load(paste0(workdir,"/data/starmap/layer_seg/starmap_layer_srtsim_simdata_rpt",irpt,"_joint.rds"))

	# colnames(simCount_mat) <- paste0("cell",1:ncol(simCount_mat))
	scobj <- CreateSeuratObject(counts = simCount_mat, project ="starmap", min.cells = 3, min.features = 100)
	scobj <- SCTransform(scobj, assay = "RNA", verbose = FALSE)
	scobj <- RunPCA(scobj, assay = "SCT", verbose = FALSE)
	scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:30)

	for(isid in 1:10){
		# isid = 1
		scd_loop <- scobj
		for(iresol in resol_list){
			scd_loop <- FindClusters(scd_loop, verbose = FALSE,resolution=iresol,random.seed=isid)
		}
		write.table(scd_loop@meta.data,file=paste0(workdir,"/output/starmap/layer_seg/seurat/starmap_layer_srtsim_simdata_rpt",irpt,"_seed",isid,"_joint_seurat.txt"))
		rm(scd_loop)
	}	
	rm(scobj,simCount_mat)
}




