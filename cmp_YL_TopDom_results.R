require(ggpubr)
require(ggsci)

outFolder <- "CMP_YL_TOPDOM_RESULTS"
dir.create(outFolder, recursive = TRUE)
myHeight <- 5
myWidth <- 7
plotType <- "svg"

plotTypeP <- "png"
myHeightP <- myWidthP <- 400
plotCex <- 1.2


# Rscript cmp_YL_TopDom_results.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
  

td_folder <- "v2_TopDom_Cancer_HiC_data_TAD_DA"
yl_folder <- "v2_Yuanlong_Cancer_HiC_data_TAD_DA"

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"

setDir <- "/media/electron"
setDir <- ""
############################################################################### COMPARE TAD P-VALUES

td_pvals <- get(load(file.path(setDir,
  "/mnt/etemp/marie", td_folder, 
  "PIPELINE/OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata")))

td_pvals <- p.adjust(td_pvals, method="BH")

yl_pvals <- get(load(file.path(setDir,
  "/mnt/etemp/marie", yl_folder, 
  "PIPELINE/OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata")))

yl_pvals <- p.adjust(yl_pvals, method="BH")

pvals_dt <- data.frame(
  pval_log10 = c(-log10(td_pvals), -log10(yl_pvals)),
  tads_data = c(rep("TopDom", length(td_pvals)), rep("YL", length(yl_pvals))),
  stringsAsFactors =FALSE
)
legTitle <- "TAD caller"
plotTit <- "Comparison TAD adj. p-vals."

my_cols <- c(pal_jama()(5)[c(3, 2)])

signifThresh <- 0.01

tmp1 <- by(pvals_dt, pvals_dt$tads_data, function(x) sum(x$pval_log10 >= -log10(signifThresh)))
tmp2 <- by(pvals_dt, pvals_dt$tads_data, function(x) length(x$pval_log10))

mySub <- paste0("# TADs p-val <= ", signifThresh, ": ", paste0(names(tmp1), "= ", as.numeric(tmp1), "/", as.numeric(tmp2), collapse="; "))

p_pvals <- ggdensity(pvals_dt,
                  title="TAD pvals [-log10]",
                  x = "pval_log10",
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = "TAD p-val [-log10]", 
                  # add = "median",                  # Add median line. 
                  rug = FALSE,                      # Add marginal rug
                  color = "tads_data", 
                  fill = "tads_data",
                  palette = "jco"
) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(values=my_cols)+
  ggtitle(plotTit, subtitle = mySub)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density")  +
  theme(  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
          plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"))
  

outFile <- file.path(outFolder, paste0("cmp_YL_TopDom_TADpvals_density.", plotType))
ggsave(p_pvals, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))


check_dt <- get(load(file.path(setDir,
                               "/mnt/etemp/marie", yl_folder, 
                               "CREATE_FINAL_TABLE//all_result_dt.Rdata")))
stopifnot(sum(check_dt$adjPvalComb[check_dt$hicds == hicds & check_dt$exprds == exprds] <= signifThresh) == tmp1["YL"])



# GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata

gt_td_dt <- get(load(file.path("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))
gt_yl_dt <- get(load(file.path("..", yl_folder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))

gt_dt <- merge(gt_td_dt, gt_yl_dt, all.x=F, all.y=F,by=c("hicds", "exprds", "entrezID"), suffixes=c("_td", "_yl"))

nDS <- length(unique(file.path(gt_dt$hicds, gt_dt$exprds)))

all_plot_vars <- c("gene_rank", "tad_rank", "tad_adjCombPval", "adj.P.Val")

for(plot_var in all_plot_vars) {
  my_x <- gt_dt[,paste0(plot_var, "_yl")]
  my_y <- gt_dt[,paste0(plot_var, "_td")]
  
  outFile <- file.path(outFolder, paste0(plot_var, "_TopDom_vs_YL_densplot.", plotTypeP)) 
  do.call(plotTypeP, list(outFile, height=myHeightP, width=myWidthP))
  densplot(x=my_x, y=my_y,
           main=paste0(plot_var),
           xlab=paste0("YL TADs"), ylab=paste0("TopDom TADs"),
       cex.main=plotCex,cex.lab=plotCex,cex.axis=plotCex,
       cex=0.7, pch=16)
  mtext(side=3, text=paste0("# DS = ", nDS))
  addCorr(x=my_x,y=my_y,bty="n")
  
  foo <- dev.off()
  cat(paste0("...written:", outFile,"\n"))
    
}




