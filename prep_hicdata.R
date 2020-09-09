options(scipen=100)

# Rscript prep_hicdata.R

require(Matrix)
require(doMC)
registerDoMC(40)
require(foreach)


all_chromos <- paste0("chr", c(1:22))

setDir <- "/media/electron"
setDir <- ""

binSize <- 40000

cell_line <- "mega_ENCSR489OCU_NCI-H460"
cell_line <- "mega_ENCSR444WCZ_A549"
cell_line <- "Compendium_LG1"

all_cell_lines <- c(
  "AWS_Barutcu_MCF-10A",
  "AWS_Barutcu_MCF-7",
  "AWS_GM12878",
  "AWS_K562",
  "AWS_Rao_HCT-116_2017",
  "Compendium_LG2",
  "Compendium_LI",
  "Compendium_PA2",
  "Compendium_PA3",
  "GSE109229_BT474",
  "GSE118514_22Rv1",
  "GSE118588_Panc_beta",
  "mega_ENCSR079VIJ_G401",
  "mega_ENCSR312KHQ_SK-MEL-5",
  "mega_ENCSR346DCU_LNCaP",
  "mega_ENCSR401TBQ_Caki2",
  "mega_ENCSR444WCZ_A549",
  "mega_ENCSR489OCU_NCI-H460",
  "mega_ENCSR549MGQ_T47D",
  "mega_ENCSR862OGI_RPMI-7951",
  "mega_GSE105194_cerebellum",
  "mega_GSE105194_spinal_cord",
  "mega_GSE105318_DLD1",
  "mega_GSE105381_HepG2",
  "mega_GSE118514_RWPE1",
  "mega_Panc1_rep12",
  "ENCSR504OTV_transverse_colon",
  "GSE109229_SKBR3",
  "GSE99051_786_O",
  "AWS_HMEC"
)
all_cell_lines=c("Compendium_LG2", "mega_GSE105194_cerebellum")


for(cell_line in all_cell_lines) {
  foo <- foreach(chromo = all_chromos) %dopar% {
    # for(chromo in all_chromos){
    
    
    
    base_file <- paste0("mat_", chromo, "_", binSize/1000, "kb_ob.txt.gz")
    
    inFile <- file.path(setDir, 
                        "/mnt/ndata/Yuanlong/2.Results/1.Juicer",
                        cell_line, 
                        "contact_mat", base_file)
    stopifnot(file.exists(inFile))
    
    out_cell_line <- gsub("^mega_", "", cell_line)
    out_cell_line <- gsub("^Compendium_", "", out_cell_line)
    out_cell_line <- gsub("^AWS_", "", out_cell_line)
    outFolder <- file.path(gsub("mega_", "", out_cell_line))
    dir.create(outFolder, recursive = TRUE)
    
    outBase <- gsub(".txt.gz", "_TopDom.txt", base_file)
    
    test_dt <- data.frame(binA = c(0, 40000, 120000),
                          binB = c(0, 80000, 200000),
                          count=c(55, 66, 77), stringsAsFactors = FALSE)
    test_dt$binA <- test_dt$binA/binSize
    test_dt$binB <- test_dt$binB/binSize
    stopifnot(test_dt$binA %%1  == 0)
    stopifnot(test_dt$binB %%1  == 0)
    
    # max_dim <- max(test_dt$binB+1)
    # stopifnot(max_dim >= (test_dt$binA+1))
    # 
    # in_mat <- sparseMatrix(i=test_dt$binA+1, j = test_dt$binB+1, x = test_dt$count)
    # in_mat
    # 4 x 6 sparse Matrix of class "dgCMatrix"
    
    # 
    # [1,] 55 .  . . .  .
    # [2,]  . . 66 . .  .
    # [3,]  . .  . . .  .
    # [4,]  . .  . . . 77
    # as.matrix(in_mat)
    # [,1] [,2] [,3] [,4] [,5] [,6]
    # [1,]   55    0    0    0    0    0
    # [2,]    0    0   66    0    0    0
    # [3,]    0    0    0    0    0    0
    # [4,]    0    0    0    0    0   77
    # as.data.frame(as.matrix(in_mat))
    # V1 V2 V3 V4 V5 V6
    # 1 55  0  0  0  0  0
    # 2  0  0 66  0  0  0
    # 3  0 66  0  0  0  0
    # 4  0  0  0  0  0 77
    # 5  0  0  0  0  0  0
    # 6  0  0  0 77  0  0
    
    # in_mat <- sparseMatrix(i=test_dt$binA+1, j = test_dt$binB+1, x = test_dt$count, dims=c(max_dim, max_dim))
    # in_mat
    
    
    in_dt <- read.table(inFile, col.names = c("binA", "binB", "count"), header=FALSE, stringsAsFactors = FALSE)
    stopifnot(in_dt$binB >= in_dt$binA)
    in_dt <- na.omit(in_dt)
    init_dt <- in_dt
    in_dt$binA <- in_dt$binA/binSize
    in_dt$binB <- in_dt$binB/binSize
    stopifnot(in_dt$binA %%1  == 0)
    stopifnot(in_dt$binB %%1  == 0)
    max_dim <- max(in_dt$binB+1)
    stopifnot(max_dim >= (in_dt$binA+1))
    in_mat <- sparseMatrix(i=in_dt$binA+1, j = in_dt$binB+1, x = in_dt$count, dims=c(max_dim, max_dim))
    in_dt <- as.data.frame(as.matrix(in_mat))
    stopifnot(dim(in_dt)[1] == dim(in_dt)[2]) 
    
    pos_dt <- data.frame(
      chromo = chromo,
      start= seq(from = 0, by = binSize, length.out=nrow(in_dt)),
      end= seq(from = binSize, by = binSize, length.out=nrow(in_dt)),
      stringsAsFactors = FALSE
    )
    
    out_dt <- cbind(pos_dt, in_dt)
    
    stopifnot(nrow(out_dt) == ncol(out_dt) - 3)
    
    # if the max of init_dt$binB is 0,
    # then in out_dt, the row will be chrX 0 40000
    stopifnot(max(init_dt$binB) == max(out_dt$start))
    
    cat(paste0("dim out_dt \t=\t", dim(out_dt), "\n"))
    
    outFile <- file.path(outFolder, outBase)
    dir.create(dirname(outFile), recursive = TRUE)
    write.table(out_dt, file=outFile, col.names = FALSE, row.names=FALSE, append=F, sep="\t", quote=FALSE)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
  
  
}

