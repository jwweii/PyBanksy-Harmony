library(tidyverse)
library(Seurat)
library(qs)
library(googlesheets4)
library(viridis, lib.loc = "/diskmnt/Projects/Users/Evan.p/tools/miniconda3/envs/seurat5_env/lib/R/library")

# Parameter - !!CHANGE THIS!!
project_name = 'mCRC_5K16'
# out
out_dir ='/diskmnt/Projects/MetNet_analysis_2/Colorectal/Xenium/Banksy/20250125_outs_5K18'


# Xenium mouse sheet processing------------------------------------------
# !!CHANGE THIS!! to the sample tile needs to filter
# Note CasetteID is used as a batch parameter in the later Banksy + Harmony workflow
gs4_deauth()
sheet =read_sheet('https://docs.google.com/spreadsheets/d/1kKb-DC3jobabcp29VNkkrN024XcHBDavkiFEBxKsBTg/edit?usp=sharing', sheet = 'Xenium')
sheet_mouse = sheet %>% 
	filter(
		Banksy_Include =='Yes',
		file.exists(`Output file path`),
		)
sheet_use = sheet_mouse[,c('Section ID', 'Casette ID', 'Output file path')] %>%  # Casette ID as batch
	setNames(c('SectionID','CasetteID','FilePath'))


## Create object ------------------------------------------
# functions: Timer, Slim seurat, Shared genes, getSize
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/2_SlimMergedSeurat/script/src/function_SlimSeurat_v1.R')
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/1_LeanSeurat/script/src/function_Timer.R')
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/1_LeanSeurat/script/src/function_sharedGenes.R')
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/1_LeanSeurat/script/src/function_getSize.R')


# Load samples
# Oneliner to create object 
Timer()
sheet_use_formated = sheet_use %>% dplyr::rename(SampleID = SectionID, FolderPath = FilePath)
obj_merged = MakeMergeXeniumSlim(sheet_use_formated)
obj_merged = AddSampleLevelMeta(
	obj_merged, 
	meta_df = sheet_use_formated, 
	meta_sample_column = 'SampleID', 
	obj_sample_column = 'orig.ident')
Timer()

genes_to_remove <- read.delim("/diskmnt/Projects/MetNet_analysis_2/Colorectal/Xenium/src/banksy/python/5K_gene_to_remove.tsv", 
                              header = FALSE, 
                              stringsAsFactors = FALSE)

genes_to_remove <- genes_to_remove[[1]]
genes_to_remove <- gsub("_", "-", genes_to_remove)
print(genes_to_remove)

obj_merged <- subset(obj_merged, features = setdiff(rownames(obj_merged), genes_to_remove))
obj_merged


getSize(obj_merged)
qsave(obj_merged, str_glue('{out_dir}/1_xenium_slim_{project_name}.qs', nthreads = 30))


# Log Noramlize Data
obj_merged <- subset(obj_merged, subset = nCount_Xenium >= 20)
obj_merged <- subset(obj_merged, subset = nFeature_Xenium >= 10)
obj_merged_data = NormalizeData(obj_merged)
getSize(obj_merged_data) # Obj Size increased to 9.9GB. On Disk 2.5 GB
qsave(obj_merged_data, str_glue('{out_dir}/1_xenium_slim_{project_name}_normalized.qs', nthreads = 30))



# Plot with custom function
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Plot/ImagePlotFaster_v1.R')
pdf(str_glue('{out_dir}/2_nCount_spatial.pdf'), width = 30, height = 30)
	print(ImagePlotFasterAllFOV(obj_merged, features = 'nCount_Xenium'))
dev.off()