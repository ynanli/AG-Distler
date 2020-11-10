library(tidyverse)

SE58095_GENISOS_sample_list <- GSE58095_GENISOS_sample_list %>%
  mutate(col = paste(.$group, .$`Diffuse vs limited`, sep = "_"))
GENISOS <- GSE58095_neqc_GENISOS_normalized %>% `colnames<-`(c("Gene", SE58095_GENISOS_sample_list$col))

## retrieve ALDH gene expression ##
GENISOS_ALDH <- GENISOS[grep("ALDH", GENISOS$Gene), ]
GENISOS_ALDH <- GENISOS_ALDH[order(rowMeans(GENISOS_ALDH[-1]), decreasing = T), ]
GENISOS_ALDH <- GENISOS_ALDH[!duplicated(GENISOS_ALDH$Gene),]

#generate heapmap for ALDH CPM
GENISOS_CPM <- select(GENISOS_ALDH, c( ends_with(c("NA", "Lim", "Dif"))))
row.names(GENISOS_CPM) <- GENISOS_ALDH$Gene

library(ComplexHeatmap)
Heatmap(GENISOS_CPM, name = "log2 CPM",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(GENISOS_CPM,
                                               rep(c("control", "lcSSc", "dcSSc"), 
                                                   times = c(36, 19, 47))),
        row_dend_reorder = T, column_dend_reorder = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = 10))

## calculate fold change for heatmap ##
con <- GENISOS_ALDH %>% select(c("Gene", ends_with("_NA"))) %>% 
  mutate(con =rowMeans(.[, -1])) %>%
  select(c("Gene", "con"))
GENISOS_FC <- merge(GENISOS_ALDH, con, by="Gene")
GENISOS_FC <- GENISOS_FC[-1]-GENISOS_FC$con
row.names(GENISOS_FC) <- GENISOS_ALDH$Gene
GENISOS_FC <- select(GENISOS_FC, c( ends_with(c("NA", "Lim", "Dif"))))

Heatmap(GENISOS_FC, name = "log2 FC",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(GENISOS_FC,
                                               rep(c("control", "lcSSc", "dcSSc"), 
                                                   times = c(36, 19, 47))),
        row_dend_reorder = T, column_dend_reorder = T, show_column_dend = F,
        row_names_gp = gpar(fontsize = 10))
