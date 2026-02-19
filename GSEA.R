rm(list = ls())

library(openxlsx)
library(dplyr)

GSEA <- read.xlsx("GSEA.xlsx", rowNames = F)
GSEA <- dplyr::arrange(GSEA, desc(NES))

Pathway <- c("_nf", "pi3k|akt", "mapk", "jak", "tgf", "wnt", "notch", "hedgehog", "hippo", "_il", "tp53")
GSEA_output2 <- GSEA[0, 0]
for(i in Pathway){
  GSEA_output3 <- grep(i, GSEA$ID, ignore.case = T, value = T)
  GSEA_output2 <- rbind(GSEA_output2, GSEA[GSEA$ID %in% GSEA_output3,])
}
GSEA_output2_adjp <- GSEA_output2[GSEA_output2$p.adjust < 0.05&GSEA_output2$qvalue < 0.25,]
GSEA_output2_p <- GSEA_output2[GSEA_output2$pvalue < 0.05,]
if(nrow(GSEA_output2_adjp) > 4){
  write.xlsx(GSEA_output2_adjp, "GSEA_output_adjp.xlsx")
}else{
    write.xlsx(GSEA_output2_p, "GSEA_output_pvalue.xlsx")} 




