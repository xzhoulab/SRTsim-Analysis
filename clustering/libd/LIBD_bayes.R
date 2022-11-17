rm(list=ls())
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)



library(Matrix)


convert_grid <- function(x){
    uniqx <- unique(x)
    uniqx_order <- uniqx[order(uniqx)]
    return(as.numeric(factor(x,levels=uniqx_order)))
}


ischeme        <- 1
irpt      <- 1
isid      <- 1

scheme_list <- c("joint","separate")

LIBD_layers <- c(paste0("Layer",1:6),"WM")

isam = "151673"

if(!file.exists(paste0(workdir,"/output/simLIBD/segmentation/BayesSpace/LIBD_",isam,"_srtsim_simdata_rpt",irpt,"_seed",isid,"_",scheme_list[ischeme],"_BayesSpace.txt"))){
    

load(paste0(workdir,"/data/simLIBD/newsimLIBD/spatialLIBD_",isam,"_unknown_clean_srtsim_simdata_rpt",irpt,"_",scheme_list[ischeme],".rds"))

numCluster <- length(table(simLoc_df$label))
counts <- simCount_mat

colnames(simLoc_df) <- c("row","col","label")

sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            colData=simLoc_df[,1:2])


set.seed(isid)
bs_obj <- sce
bs_obj <- spatialPreprocess(bs_obj, platform="Visium", n.HVGs=2000, log.normalize=TRUE)
bs_obj <- spatialCluster(bs_obj, q=numCluster, platform="Visium",
                           init.method="mclust", model="t", 
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)
comp_df <- as.data.frame(colData(bs_obj))

write.table(comp_df,file=paste0(workdir,"/output/simLIBD/newLIBD/BayesSpace/LIBD_",isam,"_srtsim_simdata_rpt",irpt,"_seed",isid,"_",scheme_list[ischeme],"_BayesSpace.txt"))

}


