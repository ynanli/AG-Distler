library(tidyverse)


PRESS <- PREE_SSc_Total_normailized_data

## retrieve ALDH fold change ##
con <- PRESS %>% select(c("X1", ends_with("_N"))) %>% 
  mutate(con =rowMeans(.[, -1])) %>%
  select(c("X1", "con"))
PRESS <- merge(PRESS_annotation, PRESS, by.x = "Ensembl_ID", by.y = "X1") %>%
  .[grep("ALDH", .$Gene_Name), ] %>%
  merge(., con, by.x = "Ensembl_ID", by.y="X1") %>% distinct()

#heatmap for logCPM
PRESS_CPM <- select(PRESS, c("Gene_Name", ends_with(c("_N", "_P"))))
row.names(PRESS_CPM) <- PRESS_CPM$Gene_Name
PRESS_CPM <- PRESS_CPM[-1]

library(ComplexHeatmap)
Heatmap(PRESS_CPM, name = "log2 CPM",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(PRESS_CPM,
                                               rep(c("control", "SSc"), 
                                                   times = c(33,58))),
        row_dend_reorder = T, column_dend_reorder = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = 10))

## heatmap for log fold change ##
PRESS_FC <- PRESS[, 3:94]-PRESS$con
row.names(PRESS_FC) <- PRESS$Gene_Name
PRESS_FC <- select(PRESS_FC, c(ends_with(c("_N", "_P"))))

Heatmap(PRESS_FC, name = "log2 FC",rect_gp = gpar(col = "white",lwd=0.5),
        cluster_columns = cluster_within_group(PRESS_FC,
                                               rep(c("control", "SSc"), 
                                                   times = c(33,58))),
        row_dend_reorder = T, column_dend_reorder = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = 10))
