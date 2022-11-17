
### joint
##+=======================
rm(list=ls())

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)



isam 	= "151673"



resol_list <- 1:30/20

for(irpt in 1:10){
	# irpt = 1
		if(!file.exists(paste0(workdir,"/output/simLIBD/newLIBD/seurat/LIBD_",isam,"_srtsim_simdata_rpt",irpt,"_seed10_joint_seurat.txt"))){
		load(paste0(workdir,"/data/simLIBD/newsimLIBD/spatialLIBD_",isam,"_unknown_clean_srtsim_simdata_rpt",irpt,"_joint.rds"))

		scobj <- CreateSeuratObject(counts = simCount_mat, project =isam, min.cells = 3, min.features = 100)
		scobj <- SCTransform(scobj, assay = "RNA", verbose = FALSE)
		scobj <- RunPCA(scobj, assay = "SCT", verbose = FALSE)
		scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:30)

		for(isid in 1:10){
			# isid = 1
			scd_loop <- scobj
			for(iresol in resol_list){
				scd_loop <- FindClusters(scd_loop, verbose = FALSE,resolution=iresol,random.seed=isid)
			}
			write.table(scd_loop@meta.data,file=paste0(workdir,"/output/simLIBD/newLIBD/seurat/LIBD_",isam,"_srtsim_simdata_rpt",irpt,"_seed",isid,"_joint_seurat.txt"))
			rm(scd_loop)
		}	
		rm(scobj,simCount_mat)
	}
}





### separate
##+=======================
rm(list=ls())

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

isam 	= "151673"



resol_list <- 1:30/20

for(irpt in 1:10){
	# irpt = 1
		if(!file.exists(paste0(workdir,"/output/simLIBD/newLIBD/seurat/LIBD_",isam,"_srtsim_simdata_rpt",irpt,"_seed10_separate_seurat.txt"))){
		load(paste0(workdir,"/data/simLIBD/newsimLIBD/spatialLIBD_",isam,"_unknown_clean_srtsim_simdata_rpt",irpt,"_separate.rds"))

		scobj <- CreateSeuratObject(counts = simCount_mat, project =isam, min.cells = 3, min.features = 100)
		scobj <- SCTransform(scobj, assay = "RNA", verbose = FALSE)
		scobj <- RunPCA(scobj, assay = "SCT", verbose = FALSE)
		scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:30)

		for(isid in 1:10){
			# isid = 1
			scd_loop <- scobj
			for(iresol in resol_list){
				scd_loop <- FindClusters(scd_loop, verbose = FALSE,resolution=iresol,random.seed=isid)
			}
			write.table(scd_loop@meta.data,file=paste0(workdir,"/output/simLIBD/newLIBD/seurat/LIBD_",isam,"_srtsim_simdata_rpt",irpt,"_seed",isid,"_separate_seurat.txt"))
			rm(scd_loop)
		}	
		rm(scobj,simCount_mat)
	}
}







