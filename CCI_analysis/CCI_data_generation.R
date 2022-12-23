rm(list=ls())
library(dplyr)
library(SRTsim)

isid = 1
foldchange = 3 # foldchange, 1, 5, 10
set.seed(isid)

load(paste0("~/starmap_1020_0410_seurat_filter_layer.rds"))
info2   <- info %>% select(c(x,y,label)) %>% 
				filter(label!="Smc") %>% as.data.frame()

region_celltype_df <- matrix(0.1,4,4)
diag(region_celltype_df) <- 0.7
rownames(region_celltype_df) <- paste0("Region",1:4)
colnames(region_celltype_df) <- paste0("Celltype",1:4)

nGene = 2000
nLoc = 5000

##========================
## simulate spatial spot
##========================
simLoc <- simNewLocs(newN= nLoc,lay_out="random",preLoc = info2)

## assign the region labels
simLoc %<>% mutate(region_label = case_when(
		x <= 0.5*median(x) ~ "Region1",
		x > 0.5*median(x) & x <= median(x)~ "Region2",
		x > median(x) & x <= 1.5*median(x)~ "Region3",
		TRUE ~ "Region4"
	)) %>% as.data.frame()


## load the Ligand-Receptor data
load("~/LR_dat.RData")



##========================
## Data Generation
##========================

## reference free
example_CCI_free = srtsim_cci_free(
					zero_prop_in = 0,
					disper_in = Inf,
					mu_in = 1, 
					numGene = nGene,
					location_in  = simLoc[,c("x","y","region_label")],
					# region_label = region_label,
					region_cell_map = region_celltype_df,
					sim_seed = isid,
					fc = foldchange,
					LR_in = LR_dat)



## reference based
simSRT  <- createSRT(count_in=sp_count[,rownames(info2)],loc_in =info2)
simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")

example_CCI_ref = srtsim_cci_ref(
					EstParam = simSRT1@EstParam,
					numGene = nGene,
					location_in  = simLoc[,c("x","y","region_label")],
					# region_label = region_label,
					region_cell_map = region_celltype_df,
					sim_seed = isid,
					fc = foldchange,
					LR_in = LR_dat)



