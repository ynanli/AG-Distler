## load packages ##
library(biomaRt)
library(tidyverse)

## retrieve HGNC symbol ##
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- ALDH3A_rb_normalized_data$X1
genes <- getBM(attributes = c('ensembl_gene_id',
                       "hgnc_symbol"),
        filters = 'ensembl_gene_id', 
        values = gene_ids,
        mart = mart)
ALDH3A2_NormCount<-ALDH3A_rb_normalized_data %>% merge(genes,
                          by.x = "X1", by.y = "ensembl_gene_id") %>%
  filter(hgnc_symbol != "")
write.csv(ALDH3A2_NormCount, "ALDH3A2_NormCount.csv", row.names = F)
