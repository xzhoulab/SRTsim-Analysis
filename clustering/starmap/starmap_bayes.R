##==========================
## joint - layer
##==========================

rm(list=ls())
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)

for(irpt in 1:10){
	load(paste0(workdir,"/data/starmap/layer_seg/starmap_layer_srtsim_simdata_rpt",irpt,"_joint.rds"))
	numCluster <- length(table(simLoc_df$label))
	counts <- simCount_mat
	colData <- simLoc_df[,c("x","y")]
	colnames(colData) <- c("row","col")
	sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
	                            colData=colData)

	for(isid in 1:10){
		set.seed(isid)
		bs_obj <- sce
		bs_obj <- spatialPreprocess(bs_obj, platform="Visium", n.HVGs=2000, log.normalize=TRUE)
		bs_obj <- spatialCluster(bs_obj, q=numCluster, platform="Visium",
		                           init.method="mclust", model="t", 
		                           nrep=10000, burn.in=1000,
		                           save.chain=TRUE)
		comp_df <- as.data.frame(colData(bs_obj))

		write.table(comp_df,file=paste0(workdir,"/output/starmap/layer_seg/BayesSpace/starmap_layer_srtsim_simdata_rpt",irpt,"_seed",isid,"_joint_BayesSpace.txt"))
		rm(bs_obj,comp_df)
	}
	rm(simLoc_df,simCount_mat,colData,sce,counts)
}



##==========================
## separate 
##==========================

rm(list=ls())
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)


for(irpt in 1:10){
	load(paste0(workdir,"/data/starmap/layer_seg/starmap_layer_srtsim_simdata_rpt",irpt,"_separate.rds"))
	numCluster <- length(table(simLoc_df$label))
	counts <- simCount_mat
	colData <- simLoc_df[,c("x","y")]
	colnames(colData) <- c("row","col")
	sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
	                            colData=colData)

	for(isid in 1:10){
		set.seed(isid)
		bs_obj <- sce
		bs_obj <- spatialPreprocess(bs_obj, platform="Visium", n.HVGs=2000, log.normalize=TRUE)
		bs_obj <- spatialCluster(bs_obj, q=numCluster, platform="Visium",
		                           init.method="mclust", model="t", 
		                           nrep=10000, burn.in=1000,
		                           save.chain=TRUE)
		comp_df <- as.data.frame(colData(bs_obj))

		write.table(comp_df,file=paste0(workdir,"/output/starmap/layer_seg/BayesSpace/starmap_layer_srtsim_simdata_rpt",irpt,"_seed",isid,"_separate_BayesSpace.txt"))
		rm(bs_obj,comp_df)
	}
	rm(simLoc_df,simCount_mat,colData,sce,counts)
}

