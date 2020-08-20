options(scipen=100)

# Rscript run_topdom.R

require(TopDom)
require(foreach)
require(doMC)
registerDoMC(40)

cell_line <- "ENCSR489OCU_NCI-H460"

binSize <- 40000

all_chromos <- paste0("chr", c(1:22))


outFolder <- file.path(paste0(cell_line, "_", binSize/1000, "kb"), "FINAL_DOMAINS")
dir.create(outFolder, recursive = TRUE)

# output should look like for yuanlong
# head ../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr21_YL_40kb_final_domains.txt 
# chr21   9640001 10640000
# chr21   10680001        10960000


foo <- foreach(chromo = all_chromos) %dopar% {
  
  ##### 1) run TopDom  
  topDom_file <- file.path(cell_line, paste0("mat_", chromo, "_", binSize/1000, "kb_ob_TopDom.txt"))
  stopifnot(file.exists(topDom_file))
  topDom_out <- TopDom(topDom_file, window.size=5)
  
  
  ##### 2) prepare output
  topDom_tads_dt <- topDom_out[["bed"]]
  topDom_tads_dt <- topDom_tads_dt[as.character(topDom_tads_dt$name) == "domain",]
  domainsDT <- topDom_tads_dt[,c("chrom", "chromStart","chromEnd")]
  stopifnot(is.numeric(domainsDT$chromStart))
  stopifnot(is.numeric(domainsDT$chromEnd))
  
  stopifnot( (domainsDT$chromStart/binSize)%%1 == 0)
  stopifnot( (domainsDT$chromEnd/binSize)%%1 == 0)
  stopifnot( ( (domainsDT$chromEnd-domainsDT$chromStart)/binSize)%%1 == 0)
  
  # to 1-based start positions:
  domainsDT$chromStart <- domainsDT$chromStart + 1 
  # ensure ordering
  domainsDT <- domainsDT[order(domainsDT$chromStart, domainsDT$chromEnd),]
  stopifnot(domainsDT$end > domainsDT$start)
  stopifnot(diff(domainsDT$end) > 0)
  stopifnot(diff(domainsDT$start) > 0)
  
  outFile <- file.path(outFolder, paste0(cell_line, "_", chromo, "_TopDom_", binSize/1000, "kb_final_domains.txt"))
  write.table(domainsDT, file=outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
  
  cat(paste0("... ", nrow(domainsDT) , " domains written in: ", outFile, "\n"))
}