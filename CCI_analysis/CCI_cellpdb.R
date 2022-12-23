

cellphonedb_func_simu = function(count_in, celltype_in,region_in, path_data,path_cellphonedb, addspatial=TRUE){
	

	dir.create(file.path(path_data), showWarnings = TRUE)
	setwd(path_data)
	
	# make count file
	write.table(count_in, file = paste0(path_data,"/count.txt"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
	print("check0")
	# make meta file
	celltype_in = data.frame("Cell"=colnames(count_in),"cell_type"=celltype_in)
	write.table(celltype_in, file = paste0(path_data,"/meta.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
	print("check1")
	# make microenvironment file
	region = data.frame("cell_type"=celltype_in$cell_type,"microenvironment"=paste0(region_in))
	region=unique(region)
	write.table(region, file = paste0(path_data,"/microenvironment.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
	
	print("check2")
	# load package
	system("source activate cpdb3.8")
	system(paste0("source ",path_cellphonedb))

	print("check3")
	# prepare cellphoneDB input 
	A="cellphonedb method statistical_analysis "
	metafile=paste0(" ",path_data,"/meta.txt ")
	countfile=paste0(" ",path_data,"/count.txt ")
	microenvs_file = paste0(" ",path_data,"/microenvironment.txt ")
	outfolder = paste0(path_data,"/out")

	# run cellphoneDB
	if(addspatial==TRUE){
		print("check4")
		system(paste0(A, metafile, countfile, " --microenvs ", microenvs_file,  " --counts-data hgnc_symbol  ", "--output-path ",outfolder))
	}else if(addspatial==FALSE){
		print("check5")
		system(paste0(A, metafile, countfile, " --counts-data hgnc_symbol  ", "--output-path ",outfolder))
	}
	print("check6")
	res_pval = read.table(paste0(path_data,"/out/pvalues.txt"),header = TRUE, fill = TRUE)
	cellphonedb_res = na.omit(reshape::melt(res_pval[,c(2,12:dim(res_pval)[2])]))
	cellphonedb_res_sig = unique(cellphonedb_res[which(cellphonedb_res$value<0.05),c(1:2)])
	compare_cellphonedb_res_sig = paste0(cellphonedb_res_sig[,1],"_",cellphonedb_res_sig[,2])
	compare_true = paste0(LR_dat[,1],"_",LR_dat[,2],"_","celltype_",LR_dat[,3],".","celltype_",LR_dat[,4])
	common_pair=length(intersect(compare_cellphonedb_res_sig, compare_true))
	return(list("cellphonedb_result"=cellphonedb_res,"cellphonedb_res_sig"=cellphonedb_res_sig,"common_pair"=common_pair))

}




library(tidyverse)
library(data.table)

cellphonedb_result_collect = function(result_cellphonedb,LR_dat,method, rep, fold){
    LR_cellphoneDB = unique(result_cellphonedb$cellphonedb_res_sig$interacting_pair)
    LR_AB_cellphoneDB = paste0(result_cellphonedb$cellphonedb_res_sig[,1],"_",result_cellphonedb$cellphonedb_res_sig[,2])
    TRUE_LR_AB = paste0(LR_dat[,1],"_",LR_dat[,2],"_",LR_dat[,3],".",LR_dat[,4])
    length(TRUE_LR_AB)
    length(intersect(LR_AB_cellphoneDB, TRUE_LR_AB))
    num_LR_all = dim(LR_dat)[1]
    num_LR_detected = length(unique(LR_cellphoneDB))
    num_LR_right = length(intersect(LR_cellphoneDB, paste0(LR_dat[,1],"_",LR_dat[,2])))
    num_LRAB_all = length(TRUE_LR_AB)
    num_LRAB_detected = length(unique(LR_AB_cellphoneDB))
    num_LRAB_right = length(intersect(LR_AB_cellphoneDB, TRUE_LR_AB))
    res_collect = data.frame(method, rep, fold, num_LR_all,num_LR_detected, num_LR_right, num_LRAB_all, num_LRAB_detected, num_LRAB_right )
    return(res_collect)
}

