library(tidyverse)

## calculate fold change from count data ##
con <- ALDH3A2_NormCount %>%
  select(c("hgnc_symbol", "XH9", "XH21", "XH29")) %>%
  mutate(con =rowMeans(.[, -1])) %>%
  select(c("hgnc_symbol", "con"))
ALDH3A2_FC_count <- ALDH3A2_NormCount %>% merge(., con, 
                                                by = "hgnc_symbol")
ALDH3A2_FC_count <- ALDH3A2_FC_count[!duplicated(ALDH3A2_FC_count$hgnc_symbol),]
row.names(ALDH3A2_FC_count)<- ALDH3A2_FC_count$hgnc_symbol
ALDH3A2_FC_count <- ALDH3A2_FC_count[,3:14]-ALDH3A2_FC_count$con 
 
## create heatmap in EMT geneset ##
ALDH3A2_EMT <- filter(ALDH3A2_FC_count, row.names(ALDH3A2_FC_count) %in% HALLMARK_EMT$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)

library(ComplexHeatmap)
Heatmap(ALDH3A2_EMT[,c(1:3, 7:12)], name = "log2 FC",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(ALDH3A2_EMT[,c(1:3, 7:12)],
                                               rep(c("control", "TGFb","TGFb+siALDH3A2"), 
                                                   times = c(3, 3, 3))),
        row_dend_reorder = T, column_dend_reorder = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = 4))

ALDH3A2_Fib <- filter(ALDH3A2_FC_count, row.names(ALDH3A2_FC_count) %in% FibAct$X1)
Heatmap(ALDH3A2_Fib[,c(1:3, 7:12)], name = "log2 FC",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(ALDH3A2_Fib[,c(1:3, 7:12)],
                                               rep(c("control", "TGFb","TGFb+siALDH3A2"), 
                                                   times = c(3, 3, 3))),
        row_dend_reorder = T, column_dend_reorder = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = 12))

ALDH3A2_secretion <- filter(ALDH3A2_FC_count, row.names(ALDH3A2_FC_count) %in% GO_NEGATIVE_REGULATION_OF_PEPTIDE_SECRETION_$GO_NEGATIVE_REGULATION_OF_PEPTIDE_SECRETION)
Heatmap(ALDH3A2_secretion[,c(1:3, 7:12)], name = "log2 FC",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(ALDH3A2_secretion[,c(1:3, 7:12)],
                                               rep(c("control", "TGFb","TGFb+siALDH3A2"), 
                                                   times = c(3, 3, 3))),
        row_dend_reorder = T, column_dend_reorder = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = 12))
