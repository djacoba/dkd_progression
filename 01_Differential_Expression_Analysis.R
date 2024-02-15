library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(gridExtra)

####Load data####

load("./KI-DN_clin_hist.RData")

load("./KI_DN_clinBL_transcriptomics.RData")

####Data preparation of Glom data####

glom.meta <- dn.meta[dn.meta$Compartment=="Glom",]
glom.meta <- glom.meta[order(glom.meta$Exiqon.ID),]
glom.cases <- which(dn.meta$Compartment=="Glom")
glom.counts <- dn.counts[,c(glom.cases,79)]

glom.cases <- which(dn.meta$Status=="Glom-DN")
glom.counts.ckd <- dn.counts[,c(glom.cases,79)]
glom.meta.ckd <- dn.meta[dn.meta$Status=="Glom-DN",]

glom.meta <- merge(glom.meta, hist, by.x="Individual_ID", by.y="ID")
glom.meta <- merge(glom.meta, clinvar, by.x="Individual_ID", by.y="Pat.ID")
glom.meta.ckd <- glom.meta[glom.meta$Status=="Glom-DN",]
glom.meta.ckd <- glom.meta.ckd[order(glom.meta.ckd$Exiqon.ID),]

glom.counts.ckd.filt <- glom.counts.ckd
glom.counts.ckd.filt <- glom.counts.ckd.filt[apply(glom.counts.ckd.filt[,c(1:19)], 1, function(row) all(row !=0 )), ]

glom.tpm.ckd.filt <- dn.tpm[,c(glom.cases,79)]
glom.tpm.ckd.filt <- glom.tpm.ckd.filt[apply(glom.tpm.ckd.filt[,c(1:19)], 1, function(row) all(row >=1 )), ]

tpm.pass <- rownames(glom.tpm.ckd.filt)
glom.counts.ckd.filt <- glom.counts.ckd.filt[rownames(glom.counts.ckd.filt) %in% tpm.pass, ]

tpm.cts.pass <- rownames(glom.counts.ckd.filt)
glom.tpm.ckd.filt <- glom.tpm.ckd.filt[rownames(glom.tpm.ckd.filt) %in% tpm.cts.pass, ]

classification <- readxl::read_xlsx("patient_classification2.xlsx")

left_join(glom.meta.ckd, classification) -> glom.meta.ckd

glom.meta.ckd %>% mutate(UACR_baseline = UACR_BL) -> glom.meta.ckd
glom.meta.ckd[7, 145] = 104.02

dn.tpm %>% dplyr::select(symbol) -> gene_symbol
rownames_to_column(gene_symbol, "ensembl_id") -> gene_symbol 

####Data preparation of Tub data####

tub.meta <- dn.meta[dn.meta$Compartment=="Tub",]
tub.meta <- tub.meta[order(tub.meta$Exiqon.ID),]
tub.cases <- which(dn.meta$Compartment=="Tub")
tub.counts <- dn.counts[,c(tub.cases,79)]

tub.cases <- which(dn.meta$Status=="Tub-DN")
tub.counts.ckd <- dn.counts[,c(tub.cases,79)]
tub.meta.ckd <- dn.meta[dn.meta$Status=="Tub-DN",]

tub.meta <- merge(tub.meta, hist, by.x="Individual_ID", by.y="ID")
tub.meta <- merge(tub.meta, clinvar, by.x="Individual_ID", by.y="Pat.ID")
tub.meta.ckd <- tub.meta[tub.meta$Status=="Tub-DN",]
tub.meta.ckd <- tub.meta.ckd[order(tub.meta.ckd$Exiqon.ID),]

tub.counts.ckd.filt <- tub.counts.ckd
tub.counts.ckd.filt <- tub.counts.ckd.filt[apply(tub.counts.ckd.filt[,c(1:19)], 1, function(row) all(row !=0 )), ]

tub.tpm.ckd.filt <- dn.tpm[,c(tub.cases,79)]
tub.tpm.ckd.filt <- tub.tpm.ckd.filt[apply(tub.tpm.ckd.filt[,c(1:19)], 1, function(row) all(row >=1 )), ]

tpm.pass.tub <- rownames(tub.tpm.ckd.filt)
tub.counts.ckd.filt <- tub.counts.ckd.filt[rownames(tub.counts.ckd.filt) %in% tpm.pass.tub, ]

tpm.cts.pass.tub <- rownames(tub.counts.ckd.filt)
tub.tpm.ckd.filt <- tub.tpm.ckd.filt[rownames(tub.tpm.ckd.filt) %in% tpm.cts.pass.tub, ]

classification <- readxl::read_xlsx("patient_classification2.xlsx")

left_join(tub.meta.ckd, classification) -> tub.meta.ckd

tub.meta.ckd %>% mutate(UACR_baseline = UACR_BL) -> tub.meta.ckd
tub.meta.ckd[5, 145] = 104.02


####DEA for Glom data based on eGFR slope#### 


dds.glom.egfr <- DESeqDataSetFromMatrix(countData=glom.counts.ckd.filt[,c(1:12,15:19)],
                                        colData=glom.meta.ckd[c(1:12,15:19),],
                                        design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+kidney_function_decline)

dds.glom.egfr$kidney_function_decline <- relevel(dds.glom.egfr$kidney_function_decline, ref = "non_rapid")
dds.glom.egfr <- DESeq(dds.glom.egfr)
dds.glom.egfr <- results(dds.glom.egfr)
summary(dds.glom.egfr)

dds.glom.egfr.df <- as.data.frame(dds.glom.egfr)


p1 <- ggplot(data=dds.glom.egfr.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.glom.egfr$padj < 0.05 & dds.glom.egfr$log2FoldChange > 0, na.rm=TRUE)
sum(dds.glom.egfr$padj < 0.05 & dds.glom.egfr$log2FoldChange < 0, na.rm=TRUE)


dds.glom.egfr.df <- tibble::rownames_to_column(dds.glom.egfr.df, "ensembl_id")
dds.glom.egfr.df <- left_join(dds.glom.egfr.df, gene_symbol)

writexl::write_xlsx(dds.glom.egfr.df, "redo_glom_egfr.xlsx")

keyvals <- ifelse(
  dds.glom.egfr.df$log2FoldChange < 0 & dds.glom.egfr.df$padj < 0.05, 'blue',
  ifelse(dds.glom.egfr.df$log2FoldChange > 0 & dds.glom.egfr.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

c <- EnhancedVolcano(dds.glom.egfr.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-3,6),
                     ylim = c(0,5),
                     title = "Glomeruli",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

####DEA for Tub data based on eGFR slope####

dds.tub.egfr <- DESeqDataSetFromMatrix(countData=tub.counts.ckd.filt[,c(1,3:5,7:19)],
                                       colData=tub.meta.ckd[c(1,3:5,7:19),],
                                       design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+kidney_function_decline)

dds.tub.egfr$kidney_function_decline <- relevel(dds.tub.egfr$kidney_function_decline, ref = "non_rapid")
dds.tub.egfr <- DESeq(dds.tub.egfr)
dds.tub.egfr <- results(dds.tub.egfr)
summary(dds.tub.egfr)

dds.tub.egfr.df <- as.data.frame(dds.tub.egfr)


p1 <- ggplot(data=dds.tub.egfr.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.tub.egfr$padj < 0.05 & dds.tub.egfr$log2FoldChange > 0, na.rm=TRUE)
sum(dds.tub.egfr$padj < 0.05 & dds.tub.egfr$log2FoldChange < 0, na.rm=TRUE)


dds.tub.egfr.df <- tibble::rownames_to_column(dds.tub.egfr.df, "ensembl_id")
dds.tub.egfr.df <- left_join(dds.tub.egfr.df, gene_symbol)

writexl::write_xlsx(dds.tub.egfr.df, "redo_tub_egfr.xlsx")

keyvals <- ifelse(
  dds.tub.egfr.df$log2FoldChange < 0 & dds.tub.egfr.df$padj < 0.05, 'blue',
  ifelse(dds.tub.egfr.df$log2FoldChange > 0 & dds.tub.egfr.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'

d <- EnhancedVolcano(dds.tub.egfr.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-3,6),
                     ylim = c(0,5),
                     title = "Tubulointerstitium",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

####DEA for Glom data based on CKD advancement####

dds.glom.gfrstage <- DESeqDataSetFromMatrix(countData=glom.counts.ckd.filt[,c(1:12,15:19)],
                                            colData=glom.meta.ckd[c(1:12,15:19),],
                                            design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+GFR_stage_progression)

dds.glom.gfrstage$GFR_stage_progression <- relevel(dds.glom.gfrstage$GFR_stage_progression, ref = "non_rapid")
dds.glom.gfrstage <- DESeq(dds.glom.gfrstage)
dds.glom.gfrstage <- results(dds.glom.gfrstage)
summary(dds.glom.gfrstage)

dds.glom.gfrstage.df <- as.data.frame(dds.glom.gfrstage)


p1 <- ggplot(data=dds.glom.gfrstage.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.glom.gfrstage$padj < 0.05 & dds.glom.gfrstage$log2FoldChange > 0, na.rm=TRUE)
sum(dds.glom.gfrstage$padj < 0.05 & dds.glom.gfrstage$log2FoldChange < 0, na.rm=TRUE)


dds.glom.gfrstage.df <- tibble::rownames_to_column(dds.glom.gfrstage.df, "ensembl_id")
dds.glom.gfrstage.df <- left_join(dds.glom.gfrstage.df, gene_symbol)

writexl::write_xlsx(dds.glom.gfrstage.df, "redo_glom_gfrstage.xlsx")

keyvals <- ifelse(
  dds.glom.gfrstage.df$log2FoldChange < 0 & dds.glom.gfrstage.df$padj < 0.05, 'blue',
  ifelse(dds.glom.gfrstage.df$log2FoldChange > 0 & dds.glom.gfrstage.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

e <- EnhancedVolcano(dds.glom.gfrstage.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-2.5,3.5),
                     ylim = c(0,7),
                     title = "Glomeruli",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

####DEA for Tub data based on CKD advancement####

dds.tub.gfrstage <- DESeqDataSetFromMatrix(countData=tub.counts.ckd.filt[,c(1,3:5,7:19)],
                                           colData=tub.meta.ckd[c(1,3:5,7:19),],
                                           design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+GFR_stage_progression)

dds.tub.gfrstage$GFR_stage_progression <- relevel(dds.tub.gfrstage$GFR_stage_progression, ref = "non_rapid")
dds.tub.gfrstage <- DESeq(dds.tub.gfrstage)
dds.tub.gfrstage <- results(dds.tub.gfrstage)
summary(dds.tub.gfrstage)

dds.tub.gfrstage.df <- as.data.frame(dds.tub.gfrstage)


p1 <- ggplot(data=dds.tub.gfrstage.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.tub.gfrstage$padj < 0.05 & dds.tub.gfrstage$log2FoldChange > 0, na.rm=TRUE)
sum(dds.tub.gfrstage$padj < 0.05 & dds.tub.gfrstage$log2FoldChange < 0, na.rm=TRUE)


dds.tub.gfrstage.df <- tibble::rownames_to_column(dds.tub.gfrstage.df, "ensembl_id")
dds.tub.gfrstage.df <- left_join(dds.tub.gfrstage.df, gene_symbol)

writexl::write_xlsx(dds.tub.gfrstage.df, "redo_tub_gfrstage.xlsx")

keyvals <- ifelse(
  dds.tub.gfrstage.df$log2FoldChange < 0 & dds.tub.gfrstage.df$padj < 0.05, 'blue',
  ifelse(dds.tub.gfrstage.df$log2FoldChange > 0 & dds.tub.gfrstage.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'

f <- EnhancedVolcano(dds.tub.gfrstage.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-2.5,3.5),
                     ylim = c(0,7),
                     title = "Tubulointerstitium",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)


####DEA for Glom data based on UACR change####

dds.glom.uacr <- DESeqDataSetFromMatrix(countData=glom.counts.ckd.filt[,c(1,2,4:7,9:12,15,19)],
                                        colData=glom.meta.ckd[c(1,2,4:7,9:12,15,19),],
                                        design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+UACR_change)

dds.glom.uacr$UACR_change <- relevel(dds.glom.uacr$UACR_change, ref = "decrease")
dds.glom.uacr <- DESeq(dds.glom.uacr)
dds.glom.uacr <- results(dds.glom.uacr)
summary(dds.glom.uacr)

dds.glom.uacr.df <- as.data.frame(dds.glom.uacr)


p1 <- ggplot(data=dds.glom.uacr.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.glom.uacr$padj < 0.05 & dds.glom.uacr$log2FoldChange > 0, na.rm=TRUE)
sum(dds.glom.uacr$padj < 0.05 & dds.glom.uacr$log2FoldChange < 0, na.rm=TRUE)


dds.glom.uacr.df <- tibble::rownames_to_column(dds.glom.uacr.df, "ensembl_id")
dds.glom.uacr.df <- left_join(dds.glom.uacr.df, gene_symbol)

writexl::write_xlsx(dds.glom.uacr.df, "redo_glom_uacr2.xlsx")


keyvals <- ifelse(
  dds.glom.uacr.df$log2FoldChange < 0 & dds.glom.uacr.df$padj < 0.05, 'blue',
  ifelse(dds.glom.uacr.df$log2FoldChange > 0 & dds.glom.uacr.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

g <- EnhancedVolcano(dds.glom.uacr.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-8.5,5.5),
                     ylim = c(0,21),
                     title = "Glomeruli",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)


####DEA for Tub data based on UACR change####

dds.tub.uacr <- DESeqDataSetFromMatrix(countData=tub.counts.ckd.filt[,c(4,5,7:13,15,16,18)],
                                       colData=tub.meta.ckd[c(4,5,7:13,15,16,18),],
                                       design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+UACR_change)

dds.tub.uacr$UACR_change <- relevel(dds.tub.uacr$UACR_change, ref = "decrease")
dds.tub.uacr <- DESeq(dds.tub.uacr)
dds.tub.uacr <- results(dds.tub.uacr)
summary(dds.tub.uacr)

dds.tub.uacr.df <- as.data.frame(dds.tub.uacr)


p1 <- ggplot(data=dds.tub.uacr.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.tub.uacr$padj < 0.05 & dds.tub.uacr$log2FoldChange > 0, na.rm=TRUE)
sum(dds.tub.uacr$padj < 0.05 & dds.tub.uacr$log2FoldChange < 0, na.rm=TRUE)


dds.tub.uacr.df <- tibble::rownames_to_column(dds.tub.uacr.df, "ensembl_id")
dds.tub.uacr.df <- left_join(dds.tub.uacr.df, gene_symbol)

writexl::write_xlsx(dds.tub.uacr.df, "redo_tub_uacr2.xlsx")

keyvals <- ifelse(
  dds.tub.uacr.df$log2FoldChange < 0 & dds.tub.uacr.df$padj < 0.05, 'blue',
  ifelse(dds.tub.uacr.df$log2FoldChange > 0 & dds.tub.uacr.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'

h <- EnhancedVolcano(dds.tub.uacr.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-8.5,5.5),
                     ylim = c(0,21),
                     title = "Tubulointerstitium",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

####DEA for Glom data based on NRA####

dds.glom.alb <- DESeqDataSetFromMatrix(countData=glom.counts.ckd.filt[,c(1:7,9:12,15:19)],
                                       colData=glom.meta.ckd[c(1:7,9:12,15:19),],
                                       design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+NRA_last_followup)

dds.glom.alb$NRA_last_followup <- relevel(dds.glom.alb$NRA_last_followup, ref = "non_rapid")
dds.glom.alb <- DESeq(dds.glom.alb)
dds.glom.alb <- results(dds.glom.alb)
summary(dds.glom.alb)

dds.glom.alb.df <- as.data.frame(dds.glom.alb)


p1 <- ggplot(data=dds.glom.alb.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.glom.alb$padj < 0.05 & dds.glom.alb$log2FoldChange > 0, na.rm=TRUE)
sum(dds.glom.alb$padj < 0.05 & dds.glom.alb$log2FoldChange < 0, na.rm=TRUE)


dds.glom.alb.df <- tibble::rownames_to_column(dds.glom.alb.df, "ensembl_id")
dds.glom.alb.df <- left_join(dds.glom.alb.df, gene_symbol)

writexl::write_xlsx(dds.glom.alb.df, "redo_glom_alb.xlsx")


keyvals <- ifelse(
  dds.glom.alb.df$log2FoldChange < 0 & dds.glom.alb.df$padj < 0.05, 'blue',
  ifelse(dds.glom.alb.df$log2FoldChange > 0 & dds.glom.alb.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

i <- EnhancedVolcano(dds.glom.alb.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-4,5),
                     ylim = c(0,4.5),
                     title = "Glomeruli",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

####DEA for Tub data based on NRA####

dds.tub.alb <- DESeqDataSetFromMatrix(countData=tub.counts.ckd.filt[,c(3:5,7:19)],
                                      colData=tub.meta.ckd[c(3:5,7:19),],
                                      design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+NRA_last_followup)

dds.tub.alb$NRA_last_followup <- relevel(dds.tub.alb$NRA_last_followup, ref = "non_rapid")
dds.tub.alb <- DESeq(dds.tub.alb)
dds.tub.alb <- results(dds.tub.alb)
summary(dds.tub.alb)

dds.tub.alb.df <- as.data.frame(dds.tub.alb)


p1 <- ggplot(data=dds.tub.alb.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.tub.alb$padj < 0.05 & dds.tub.alb$log2FoldChange > 0, na.rm=TRUE)
sum(dds.tub.alb$padj < 0.05 & dds.tub.alb$log2FoldChange < 0, na.rm=TRUE)


dds.tub.alb.df <- tibble::rownames_to_column(dds.tub.alb.df, "ensembl_id")
dds.tub.alb.df <- left_join(dds.tub.alb.df, gene_symbol)

writexl::write_xlsx(dds.tub.alb.df, "redo_tub_alb.xlsx")

keyvals <- ifelse(
  dds.tub.alb.df$log2FoldChange < 0 & dds.tub.alb.df$padj < 0.05, 'blue',
  ifelse(dds.tub.alb.df$log2FoldChange > 0 & dds.tub.alb.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

j <- EnhancedVolcano(dds.tub.alb.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-4,5),
                     ylim = c(0,4.5),
                     title = "Tubulointerstitium",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

####DEA for Glom data based on outcome####

dds.glom.out <- DESeqDataSetFromMatrix(countData=glom.counts.ckd.filt[,c(1:13,15:19)],
                                       colData=glom.meta.ckd[c(1:13,15:19),],
                                       design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+composite_outcome)

dds.glom.out$composite_outcome <- relevel(dds.glom.out$composite_outcome, ref = "non_rapid")
dds.glom.out <- DESeq(dds.glom.out)
dds.glom.out <- results(dds.glom.out)
summary(dds.glom.out)

dds.glom.out.df <- as.data.frame(dds.glom.out)


p1 <- ggplot(data=dds.glom.out.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.glom.out$padj < 0.05 & dds.glom.out$log2FoldChange > 0, na.rm=TRUE)
sum(dds.glom.out$padj < 0.05 & dds.glom.out$log2FoldChange < 0, na.rm=TRUE)


dds.glom.out.df <- tibble::rownames_to_column(dds.glom.out.df, "ensembl_id")
dds.glom.out.df <- left_join(dds.glom.out.df, gene_symbol)

writexl::write_xlsx(dds.glom.out.df, "redo_glom_out2.xlsx")


keyvals <- ifelse(
  dds.glom.out.df$log2FoldChange < 0 & dds.glom.out.df$padj < 0.05, 'blue',
  ifelse(dds.glom.out.df$log2FoldChange > 0 & dds.glom.out.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

k <- EnhancedVolcano(dds.glom.out.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-5,5),
                     ylim = c(0,3),
                     title = "Glomeruli",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)



####DEA for Tub data based on outcome####

dds.tub.out <- DESeqDataSetFromMatrix(countData=tub.counts.ckd.filt[,c(1,3:19)],
                                      colData=tub.meta.ckd[c(1,3:19),],
                                      design=~Center+Gender+Age_Bx+BMI+eGFR_BL+UACR_baseline+composite_outcome)

dds.tub.out$composite_outcome <- relevel(dds.tub.out$composite_outcome, ref = "non_rapid")
dds.tub.out <- DESeq(dds.tub.out)
dds.tub.out <- results(dds.tub.out)
summary(dds.tub.out)

dds.tub.out.df <- as.data.frame(dds.tub.out)


p1 <- ggplot(data=dds.tub.out.df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() #+xlim(-3,4) #+ theme_minimal()
p2 <- p1 + geom_hline(yintercept=-log10(0.05), col="red")

sum(dds.tub.out$padj < 0.05 & dds.tub.out$log2FoldChange > 0, na.rm=TRUE)
sum(dds.tub.out$padj < 0.05 & dds.tub.out$log2FoldChange < 0, na.rm=TRUE)


dds.tub.out.df <- tibble::rownames_to_column(dds.tub.out.df, "ensembl_id")
dds.tub.out.df <- left_join(dds.tub.out.df, gene_symbol)

writexl::write_xlsx(dds.tub.out.df, "redo_tub_out2.xlsx")

keyvals <- ifelse(
  dds.tub.out.df$log2FoldChange < 0 & dds.tub.out.df$padj < 0.05, 'blue',
  ifelse(dds.tub.out.df$log2FoldChange > 0 & dds.tub.out.df$padj < 0.05, 'red',
         'gray'))

names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'downregulated'
#names(keyvals)[keyvals == 'gray'] <- 'not significant'

l <- EnhancedVolcano(dds.tub.out.df,
                     lab = "",
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-5,5),
                     ylim = c(0,3),
                     title = "Tubulointerstitium",
                     subtitle = "",
                     caption = "",
                     legendPosition = 'none',
                     subtitleLabSize = 16,
                     legendLabSize = 16,
                     legendIconSize = 6,
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 6.0,
                     colCustom = keyvals,
                     boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)



####Save figures####

grid.arrange(c, d, ncol = 2)

dev.copy(png,"volcano_egfr.png", width=1200, height=600)
dev.off()

grid.arrange(e, f, ncol = 2)

dev.copy(png,"volcano_gfr2.png", width=1200, height=600)
dev.off()


grid.arrange(g, h, ncol = 2)

dev.copy(png,"volcano_uacr.png", width=1200, height=600)
dev.off()


grid.arrange(i, j, ncol = 2)

dev.copy(png,"volcano_alb.png", width=1200, height=600)
dev.off()

grid.arrange(k, l, ncol = 2)

dev.copy(png,"volcano_out2.png", width=1200, height=600)
dev.off()