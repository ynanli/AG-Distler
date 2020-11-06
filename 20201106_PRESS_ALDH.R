library(tidyverse)
library(edgeR)

PRESS <- PREE_SSc_Total_normailized_data

## retrieve ALDH fold change ##
con <- PRESS %>% select(c("X1", ends_with("_N"))) %>% 
  mutate(con =rowMeans(.[, -1])) %>%
  select(c("X1", "con"))
PRESS <- merge(PRESS_annotation, PRESS, by.x = "Ensembl_ID", by.y = "X1") %>%
  .[grep("ALDH", .$Gene_Name), ] %>%
  merge(., con, by.x = "Ensembl_ID", by.y="X1") %>% unique()
