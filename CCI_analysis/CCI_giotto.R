## refer to the following link to install giotto
## https://rubd.github.io/Giotto_site/articles/installation_issues.html
Giotto_func_simu = function(count_in, location_in, celltype_in,LR_data,my_python_path = "~/anaconda3/bin/python3"){

	library(Giotto)
	library(peakRAM)
	library(tidyverse)
	library(data.table)
	set.seed(seed = 1234)

	celltype_in = data.frame("cell_ID"=colnames(count_in),"cell_types"=celltype_in)

	instrs <- createGiottoInstructions(python_path = my_python_path)
	seqfish_mini <- createGiottoObject(raw_exprs = count_in, spatial_locs = location_in, instructions = instrs)
	seqfish_mini <- filterGiotto(gobject = seqfish_mini, 
                             expression_threshold = 0.5, 
                             gene_det_in_min_cells = 20, 
                             min_det_genes_per_cell = 0)
	seqfish_mini <- normalizeGiotto(gobject = seqfish_mini, scalefactor = 6000, verbose = T)
	seqfish_mini <- addStatistics(gobject = seqfish_mini)
	seqfish_mini <- adjustGiottoMatrix(gobject = seqfish_mini, 
	                                   expression_values = c('normalized'),
	                                   covariate_columns = c('nr_genes', 'total_expr'))
	# celltype_in: data.frame(cellid, cell_type)
	colnames(celltype_in) = c("cell_ID","cell_types")
	seqfish_mini@cell_metadata$cell_types = celltype_in$cell_types[match(seqfish_mini@cell_metadata$cell_ID, celltype_in$cell_ID)]
		# seqfish_mini = annotateGiotto(gobject = seqfish_mini, 
	 #                              annotation_vector = celltype_in, 
	 #                              name = 'cell_types')

	seqfish_mini = createSpatialNetwork(gobject = seqfish_mini, minimum_k = 2, 
                                    maximum_distance_delaunay = 400)

	colnames(LR_data) = c("Ligand","Receptor")
	LR_data=data.table(LR_data)
	LR_data[, ligand_det := ifelse(LR_data$Ligand %in% seqfish_mini@gene_ID, T, F)] # := is function in data.table
	LR_data[, receptor_det := ifelse(LR_data$Receptor %in% seqfish_mini@gene_ID, T, F)]
	LR_data_det = LR_data[ligand_det == T & receptor_det == T]

	select_ligands = LR_data_det$Ligand
	select_receptors = LR_data_det$Receptor


	## get statistical significance of gene pair expression changes based on expression ##
	expr_only_scores = exprCellCellcom(gobject = seqfish_mini,
	                                   cluster_column = 'cell_types',
	                                   random_iter = 500,
	                                   gene_set_1 = select_ligands,
	                                   gene_set_2 = select_receptors,
	                                   verbose = F)

	## get statistical significance of gene pair expression changes upon cell-cell interaction
	spatial_all_scores = spatCellCellcom(seqfish_mini,
	                                     spatial_network_name = 'Delaunay_network',
	                                     cluster_column = 'cell_types',
	                                     random_iter = 500,
	                                     gene_set_1 = select_ligands,
	                                     gene_set_2 = select_receptors,
	                                     adjust_method = 'fdr',
	                                     do_parallel = T,
	                                     cores = 4,
	                                     verbose = "none")
	## * plot communication scores ####
	## select top LR ##
	selected_spat = spatial_all_scores[p.adj <= 0.05 ]
	data.table::setorder(selected_spat, -PI)
	top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)
	top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)
	return(list("selected_spat"=selected_spat,
				"top_LR_ints"=top_LR_ints,
				"top_LR_cell_ints"=top_LR_cell_ints,
				"spatial_all_scores"=spatial_all_scores,
				"expr_only_scores"=expr_only_scores))

}

giotto_result_collect_spatial = function(result_giotto,LR_dat,method, rep, fold){
    input_LR_pairs = unique(paste0(unlist(LR_dat[,1]),"-",unlist(LR_dat[,2])))
    LR_AB = paste0(result_giotto$selected_spat$LR_comb, "_",result_giotto$selected_spat$LR_cell_comb)
    TRUE_LR_AB = paste0(LR_dat[,1],"-",LR_dat[,2],"_",LR_dat[,3],"--",LR_dat[,4])
    num_LR_all = length(input_LR_pairs)
    num_LR_detected = length(unique(result_giotto$selected_spat$LR_comb))
    num_LR_right = length(intersect(result_giotto$selected_spat$LR_comb, input_LR_pairs))
    num_LRAB_all = length(TRUE_LR_AB)
    num_LRAB_detected = length(unique(LR_AB))
    num_LRAB_right = length(intersect(LR_AB, TRUE_LR_AB))
    res_collect = data.frame(method, rep, fold, num_LR_all,num_LR_detected, num_LR_right, num_LRAB_all, num_LRAB_detected, num_LRAB_right )
    return(res_collect)
}


giotto_result_collect_no_spatial = function(result_giotto,LR_dat,method, rep, fold){

    result_giotto$selected_spat = result_giotto$expr_only_scores[p.adj <= 0.05]
    data.table::setorder(result_giotto$selected_spat, -PI)
    top_LR_ints = unique(result_giotto$selected_spat[order(-abs(PI))]$LR_comb)
    top_LR_cell_ints = unique(result_giotto$selected_spat[order(-abs(PI))]$LR_cell_comb)
    input_LR_pairs = unique(paste0(unlist(LR_dat[,1]),"-",unlist(LR_dat[,2])))
    LR_AB = paste0(result_giotto$selected_spat$LR_comb, "_",result_giotto$selected_spat$LR_cell_comb)
    TRUE_LR_AB = paste0(LR_dat[,1],"-",LR_dat[,2],"_",LR_dat[,3],"--",LR_dat[,4])
    num_LR_all = length(input_LR_pairs)
    num_LR_detected = length(unique(result_giotto$selected_spat$LR_comb))
    num_LR_right = length(intersect(result_giotto$selected_spat$LR_comb, input_LR_pairs))
    num_LRAB_all = length(TRUE_LR_AB)
    num_LRAB_detected = length(unique(LR_AB))
    num_LRAB_right = length(intersect(LR_AB, TRUE_LR_AB))
    res_collect = data.frame(method, rep, fold, num_LR_all,num_LR_detected, num_LR_right, num_LRAB_all, num_LRAB_detected, num_LRAB_right )
    return(res_collect)
}


print("Running Giotto")
start_time <- Sys.time()
result_giotto = Giotto_func_simu(count_in=example_CCI_free@simCounts, location_in=as.data.frame(example_CCI_free@simcolData[,1:2]),celltype_in= example_CCI_free@simcolData$celltype,LR_dat[,1:2])
end_time <- Sys.time()
T=end_time-start_time

result_giotto_spatial = giotto_result_collect_spatial(result_giotto,LR_dat,method="giotto_w_spatial", rep=isid, fold=foldchange)
result_giotto_nospatial = giotto_result_collect_no_spatial(result_giotto,LR_dat,method="giotto_wo_spatial",rep=isid, fold=foldchange)

