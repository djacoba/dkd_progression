library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(devtools)
library(WeightedTreemaps)


####Run ORA using Reactome on Glom DEGs####
glom.degs <- read_delim("entrez_glom.txt")


reactome.glom <- enrichPathway(
  glom.degs$`NCBI gene ID`,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE)

reactome.glom.df <- reactome.glom@result

####Extract parent pathway name from Reactome####

pathway.relations <- read.delim("ReactomePathwaysRelation.txt", header = F)
colnames(pathway.relations) <- c("parent", "child")

pathway.list <- read.delim("ReactomePathways.txt", header = F)
colnames(pathway.list) <- c("Pathway identifier", "Pathway name", "Species")
pathway.list %>% filter(Species == "Homo sapiens") -> pathway.list
pathway.list$Species <- NULL

left_join(reactome.glom.df, pathway.relations, by = c("ID" = "child")) -> reactome.glom.df2
left_join(reactome.glom.df2, pathway.list, by = c("parent" = "Pathway identifier")) -> reactome.glom.df2

reactome.glom.df2 %>% dplyr::filter('Pathway name' != 'Disease') -> reactome.glom.df2

reactome.glom.df2 %>% mutate(`Pathway name` = case_when(is.na(`Pathway name`) ~ `Description`,
                                                        TRUE ~ `Pathway name`)) -> reactome.glom.df2

reactome.glom.df2 %>% dplyr::filter(p.adjust < 0.05) -> reactome.glom.df2

####Plot as weighted treemap####

tm <- voronoiTreemap(
  data = reactome.glom.df2,
  levels = c("Pathway name", "Description"),
  cell_size = "Count",
  shape = "rounded_rect",
  seed = 123
)

drawTreemap(tm, label_size = 9,
            color_type = "categorical", color_level = 1,
            label_level = c(2),
            border_color = grey(0.4), label_color = grey(0.1))


dev.copy(png,"glomtreemap.png", width=1200, height=1000)
dev.off()


####Run ORA using Reactome on Tub DEGs####

tub.degs <- read_delim("entrez_tub.txt")

reactome.tub <- enrichPathway(
  tub.degs$`NCBI gene ID`,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE)

reactome.tub.df <- reactome.tub@result

left_join(reactome.tub.df, pathway.relations, by = c("ID" = "child")) -> reactome.tub.df2
left_join(reactome.tub.df2, pathway.list, by = c("parent" = "Pathway identifier")) -> reactome.tub.df2

reactome.tub.df2 %>% dplyr::filter('Pathway name' != 'Disease') -> reactome.tub.df2

reactome.tub.df2 %>% mutate(`Pathway name` = case_when(is.na(`Pathway name`) ~ `Description`,
                                                       TRUE ~ `Pathway name`)) -> reactome.tub.df2

reactome.tub.df2 %>% dplyr::filter(p.adjust < 0.05) -> reactome.tub.df2


tm2 <- voronoiTreemap(
  data = reactome.tub.df2,
  levels = c("Pathway name", "Description"),
  cell_size = "Count",
  shape = "rounded_rect",
  seed = 123
)

drawTreemap(tm2, label_size = 9,
            color_type = "categorical", color_level = 1,
            label_level = c(2),
            border_color = grey(0.4), label_color = grey(0.1))


dev.copy(png,"tubtreemap.png", width=1200, height=1000)
dev.off()