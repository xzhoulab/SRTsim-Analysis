
##=================================================

rm(list=ls())
library(Giotto)

library(Matrix)


my_working_dir = paste0(workdir,"/output/simLIBD/newLIBD/")


ischeme   <- 1
irpt      <- 1
isid     <- 1

scheme_list <- c("joint","separate")

LIBD_layers <- c(paste0("Layer",1:6),"WM")

# isid = 10

isam = "151673"
numGene = 100

myinst  <- createGiottoInstructions(save_plot=F, show_plot=T, save_dir=workdir, python_path="~/anaconda3/bin/python3")

load(paste0(workdir,"/data/simLIBD/newsimLIBD/spatialLIBD_",isam,"_unknown_clean_srtsim_simdata_rpt",irpt,"_",scheme_list[ischeme],".rds"))

numCluster  <- length(table(simLoc_df$label))
obj         <- createGiottoObject(raw_exprs = simCount_mat, spatial_locs = simLoc_df[,1:2],instructions = myinst)

go_obj      <- filterGiotto(gobject = obj, expression_threshold = 1,gene_det_in_min_cells = 50, min_det_genes_per_cell = 100,expression_values = c('raw'),verbose = T)
go_obj      <- normalizeGiotto(gobject = go_obj, scalefactor = 6000, verbose = T)
## add gene & cell statistics
go_obj      <- addStatistics(gobject = go_obj)
## adjust expression matrix for technical or known variables
go_obj      <- adjustGiottoMatrix(gobject = go_obj, expression_values = c('normalized'),batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),return_gobject = TRUE,update_slot = c('custom'))
go_obj      <- calculateHVG(gobject = go_obj, method = 'cov_loess',show_plot=FALSE,save_plot=FALSE, difference_in_cov = 0.1)

## select genes based on HVG and gene statistics, both found in gene metadata
gene_metadata   = fDataDT(go_obj)
featgenes       = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID


go_obj <- runPCA(gobject = go_obj, genes_to_use = featgenes, scale_unit = F, center=T, method="factominer")

# create spatial grid
go_obj <- createSpatialGrid(gobject = go_obj, sdimx_stepsize = 400, sdimy_stepsize = 400, minimum_padding = 0)
go_obj <- createSpatialNetwork(gobject = go_obj, method = 'kNN', k = 5, maximum_distance_knn = 400, name = 'spatial_network')


kmtest <- binSpect(go_obj, spatial_network_name = 'spatial_network')


go_obj <- createSpatialNetwork(gobject = go_obj, minimum_k = 2, name = 'Delaunay_full')
# hmrf_folder <- fs::path(my_working_dir,'/giotto/hmrf/')
# hmrf_folder <- fs::path(my_working_dir,'/giotto/hmrf_all/',isam,irpt,"/")
hmrf_folder = fs::path(my_working_dir,'/giotto/hmrf/')

if(!file.exists(hmrf_folder)){dir.create(hmrf_folder, recursive = T)}


    my_spatial_genes <- kmtest[1:numGene]$genes

    # for(isid in 1:10){
        go_obj_loop <- go_obj
        # do HMRF with different betas
        HMRF_spatial_genes <- doHMRF(gobject = go_obj_loop, expression_values = 'scaled',
                                    spatial_genes = my_spatial_genes,
                                    spatial_network_name = 'Delaunay_full',
                                    k = numCluster,
                                    betas = c(0,5,20), 
                                    seed= isid,
                                    output_folder = paste0(hmrf_folder, "/", "Spatial_genes/SG_top",numGene,"_rpt",irpt,"_seed",isid,"_scaled_",scheme_list[ischeme]))
        ## add HMRF of interest to giotto object
        for(ibetas in seq(0,100,by=5)){
            go_obj_loop <- addHMRF(gobject = go_obj_loop,HMRFoutput = HMRF_spatial_genes,k = numCluster, betas_to_add = c(ibetas),hmrf_name = 'HMRF_2')
        }
        

        comp_df     <- as.data.frame(go_obj_loop@cell_metadata)
        write.table(comp_df,file=paste0(workdir,"/output/simLIBD/newLIBD/giotto/LIBD_",isam,"_srtsim_simdata_numgene",numGene,"_rpt",irpt,"_seed",isid,"_",scheme_list[ischeme],"_hmrf.txt"))



