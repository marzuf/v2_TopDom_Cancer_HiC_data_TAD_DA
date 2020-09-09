options(scipen=100)

# Rscript run_topdom.R

require(TopDom)
require(foreach)
require(doMC)
registerDoMC(40)

cell_line <- "ENCSR489OCU_NCI-H460"
cell_line <- "ENCSR444WCZ_A549"
cell_line <- "LG1"

binSize <- 40000

all_chromos <- paste0("chr", c(1:22))



all_cell_lines <- c(
  "Barutcu_MCF-10A",
  "Barutcu_MCF-7",
  "GM12878",
  "K562",
  "Rao_HCT-116_2017",
  "LG2",
  "LI",
  "PA2",
  "PA3",
  "GSE109229_BT474",
  "GSE118514_22Rv1",
  "GSE118588_Panc_beta",
  "ENCSR079VIJ_G401",
  "ENCSR312KHQ_SK-MEL-5",
  "ENCSR346DCU_LNCaP",
  "ENCSR401TBQ_Caki2",
  "ENCSR444WCZ_A549",
  "ENCSR489OCU_NCI-H460",
  "ENCSR549MGQ_T47D",
  "ENCSR862OGI_RPMI-7951",
  "GSE105194_cerebellum",
  "GSE105194_spinal_cord",
  "GSE105318_DLD1",
  "GSE105381_HepG2",
  "GSE118514_RWPE1",
  "Panc1_rep12",
  "ENCSR504OTV_transverse_colon",
  "GSE109229_SKBR3",
  "GSE99051_786_O",
  "HMEC"
)

all_cell_lines=c("LG2", "GSE105194_cerebellum")

for(cell_line in all_cell_lines) {
  

  outFolder <- file.path(paste0(cell_line, "_", binSize/1000, "kb"), "FINAL_DOMAINS")
  dir.create(outFolder, recursive = TRUE)
  
  # output should look like for yuanlong
  # head ../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr21_YL_40kb_final_domains.txt 
  # chr21   9640001 10640000
  # chr21   10680001        10960000
  
  
  foo <- foreach(chromo = all_chromos) %dopar% {
    
    ##### 1) run TopDom  
    topDom_file <- file.path(cell_line, paste0("mat_", chromo, "_", binSize/1000, "kb_ob_TopDom.txt"))
    # stopifnot(file.exists(topDom_file))
    if(!file.exists(topDom_file)) return(NULL)
    topDom_out <- TopDom(topDom_file, window.size=5)
    
    
    ##### 2) prepare output
    topDom_tads_dt <- topDom_out[["bed"]]
    topDom_tads_dt <- topDom_tads_dt[as.character(topDom_tads_dt$name) == "domain",]
    domainsDT <- topDom_tads_dt[,c("chrom", "chromStart","chromEnd")]
    stopifnot(is.numeric(domainsDT$chromStart))
    stopifnot(is.numeric(domainsDT$chromEnd))
    # !!! found a dataset with 434  chr2  242840000 242840000
    domainsDT <- domainsDT[domainsDT$chromEnd > domainsDT$chromStart,]
    
    stopifnot( (domainsDT$chromStart/binSize)%%1 == 0)
    stopifnot( (domainsDT$chromEnd/binSize)%%1 == 0)
    stopifnot( ( (domainsDT$chromEnd-domainsDT$chromStart)/binSize)%%1 == 0)
    
    # to 1-based start positions:
    domainsDT$chromStart <- domainsDT$chromStart + 1 
    # ensure ordering
    domainsDT <- domainsDT[order(domainsDT$chromStart, domainsDT$chromEnd),]
    stopifnot(domainsDT$end > domainsDT$start)
    stopifnot(diff(domainsDT$chromStart) > 0)
    stopifnot(diff(domainsDT$chromEnd) > 0)
    
    outFile <- file.path(outFolder, paste0(cell_line, "_", chromo, "_TopDom_", binSize/1000, "kb_final_domains.txt"))
    write.table(domainsDT, file=outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    
    cat(paste0("... ", nrow(domainsDT) , " domains written in: ", outFile, "\n"))
  }
}