yl <- get(load(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/0_prepGeneData/pipeline_geneList.Rdata")))
td <- get(load(file.path("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/0_prepGeneData/pipeline_geneList.Rdata")))

cat(paste0("length interesect\t=\t", length(intersect(yl, td)), "\n"))
cat(paste0("length setdiff YL\t=\t", length(setdiff(yl, td)), "/", length(yl), "\n"))
cat(paste0("length setdiff TD\t=\t", length(setdiff(td,yl)), "/", length(td), "\n"))


yl <- get(load(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/0_prepGeneData/pipeline_regionList.Rdata")))
td <- get(load(file.path("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/0_prepGeneData/pipeline_regionList.Rdata")))

cat(paste0("length YL TADs\t=\t", length(yl), "\n"))
cat(paste0("length TD TADs\t=\t", length(td), "\n"))

# I WILL NEED TO CHECK WHICH DATA USED FOR MEAN CORR !!!

