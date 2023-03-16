# this code is to practice RNAseq with GSE183947
# setwd("/Users/rezamoosavi_1/Desktop/RNAseq")

#now we need to load the libraries

library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)

# reading the data
dat <- read.csv(file = "/Users/rezamoosavi_1/Desktop/RNAseq/GSE183947_fpkm.csv")
head(dat)
dim(dat)

# importing metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
gse

#exporting the first element of gse list as metadata
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# The easier way to generate a subset of metadata
metadata.subset <- select(metadata, c(1,10,11, 17))
head(metadata.subset)

# now we need to generate a subset from metadata and trim it with pipe operator from dplyr
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

metadata.modified

# now we need to change our expression data to long format
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "samples", value = "FPKM", -gene)

# we merge our metadata.modified to our dat.long dataframe
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# explore data by extracting expression for BRCA genes
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM), median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM) %>%
  head()

# now we start to visualize our data with barplots for BRCA1
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) + geom_col()
  
# generation of density plot for BRCA1
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) + geom_density(alpha = 0.3)


# generation of a boxplot for BRCA1
dat.long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM)) + geom_boxplot()

# generation of correlation plot for BRCA1 and BRCA2
dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = "FPKM") %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# generation of a heatmap
genes.of.interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")
plot1 <- dat.long %>% 
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

ggsave(plot1, fileName("heatmap1.pdf"), width = 10, height = 8)



