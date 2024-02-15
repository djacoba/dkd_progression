library(tidyverse)
library(reshape2)
library(ggsankey)

####Load data####

load("./KI-DN_clin_hist.RData")

load("./KI_DN_clinBL_transcriptomics.RData")


####Read eGFR values####

egfr <- readxl::read_xlsx("egfr_uacr_values.xlsx")

####Time-series plot for eGFR####

wide.egfr <- t(egfr)[c(2,3,4,5,6,7),]
colnames(wide.egfr) <- c("GU-001", "GU-002", "GU-003", "GU-004", "GU-005", "GU-006", "GU-007", "GU-008", "GU-009", "GU-010", "KI-020", "KI-021", "KI-022", "KI-023", "KI-024", "KI-025", "KI-027", "KI-028", "KI-029")
df.wide.egfr <- as.data.frame(wide.egfr)
df.wide.egfr[] <- sapply(df.wide.egfr, as.numeric)

df.wide2.egfr <- data.frame(x=seq_along(df.wide.egfr[,1]), df.wide.egfr)
df.wide2.egfr <- df.wide2.egfr %>% mutate(x = x-1)

df.long.egfr <- melt(df.wide2.egfr, id.vars = "x")

ggplot(df.long.egfr, aes(x = x, y = value, color = variable)) +
  geom_point(size = 5) + geom_line(size = 1) + theme_classic() + theme(legend.position = "none") +
  scale_color_viridis_d() +
  labs(x = "Year", y = "eGFR (ml/min/year)") +
  theme(text=element_text(size=20))

dev.copy(png,"REDOegfr.png", width=600, height=600)
dev.off()


####Time-series plot for UACR####

UACR_ts <- data.frame(clinvar$Pat.ID, clinvar$UACR_BL, clinvar$UACR_FU1y, clinvar$UACR_FU2y, clinvar$UACR_FU3y, clinvar$UACR_FU4y, clinvar$UACR_FU5y)
UACR_ts$clinvar.UACR_FU4y <- as.numeric(UACR_ts$clinvar.UACR_FU4y)

wide <- t(UACR_ts)[c(2,3,4,5,6,7),]
colnames(wide) <- c("GU-001", "GU-002", "GU-003", "GU-004", "GU-005", "GU-006", "GU-007", "GU-008", "GU-009", "GU-010", "KI-020", "KI-021", "KI-022", "KI-023", "KI-024", "KI-025", "KI-027", "KI-028", "KI-029")
df.wide <- as.data.frame(wide)
as.data.frame(sapply(df.wide[, c(1:19)], as.numeric))
df.wide[] <- sapply(df.wide, as.numeric)

df.wide2 <- data.frame(x=seq_along(df.wide[,1]), df.wide)
sapply(df.wide2, class)
df.wide2 %>% mutate(x = x-1)
df.wide2 <- df.wide2 %>% mutate(x = x-1)

df.long <- melt(df.wide2, id.vars = "x")

ggplot(df.long, aes(x = x, y = log(value), color = variable)) +
  geom_point(size=6) + geom_line(size=1) + theme_classic() + theme(legend.position = "none") +
  scale_color_viridis_d() +
  labs(x = "Year", y = "log(UACR) (mg/mmol)") +
  theme(text=element_text(size=20))

dev.copy(png,"REDOuacr.png", width=600, height=600)
dev.off()


####Histogram for eGFR slope####

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

df2 <- glom.meta.ckd %>% dplyr::select(Individual_ID, eGFR_FU1y, eGFR_FU2y, eGFR_FU3y, eGFR_FU4y, eGFR_FU5y)

df2.long <- df2 %>% 
  pivot_longer(!Individual_ID, names_to="Year", values_to = "eGFR")

df2.long <- df2.long %>% 
  mutate(Year = case_when(Year == "eGFR_FU1y" ~ 1,
                          Year == "eGFR_FU2y" ~ 2,
                          Year == "eGFR_FU3y" ~ 3,
                          Year == "eGFR_FU4y" ~ 4,
                          Year == "eGFR_FU5y" ~ 5))

baseline <- glom.meta.ckd %>% dplyr::select(Individual_ID, Gender, eGFR_BL, Diabetes_type)
baseline$Diabetes_type <- as.factor(baseline$Diabetes_type)

df2.long <- merge(df2.long, baseline, by.x="Individual_ID")


fm52 <- lmerTest::lmer(eGFR ~ Year + eGFR_BL + (0 + Year|Individual_ID), df2.long)
coef.fm52 <- as.data.frame(coef(fm52)$Individual_ID)

glom.meta.ckd <- merge(glom.meta.ckd, coef.fm52, by.x="Individual_ID", by.y="row.names")
glom.meta.ckd <- dplyr::rename(glom.meta.ckd, eGFR_slope = Year)
glom.meta.ckd$eGFR_BL.y <- NULL
glom.meta.ckd$`(Intercept)` <- NULL
glom.meta.ckd <- glom.meta.ckd[order(glom.meta.ckd$Exiqon.ID),]

glom.meta.ckd[7, 84] = 104.02


ggplot(data = glom.meta.ckd, aes(x = eGFR_slope)) +
  geom_histogram(binwidth = 3.25, aes(fill = eGFR_slope < -5)) +
  geom_vline(aes(xintercept=-5),
             color="white", linetype="dashed", size=1)+
  scale_fill_manual(values = c("blue", "red")) +
  theme_classic() + labs(x = "eGFR slope (ml/min/year)", y = "Number of patients") +
  theme(text=element_text(size=20)) + theme(legend.position = "none")

dev.copy(png,"histegfr.png", width=600, height=600)
dev.off()

####Histogram for UACR change####

uacr_change <- readxl::read_xlsx("uacr_change_plot.xlsx")

uacr_change$group = ifelse(uacr_change$`UACR change (%)` > 30, 'inc', ifelse(uacr_change$`UACR change (%)` < -30, 'dec', 'exc'))

ggplot(data = uacr_change, aes(x = `UACR change (%)`, fill = group)) +
  geom_histogram(bins=60, binwidth = 25) +
  scale_fill_manual(values = c('inc' = 'red',
                               'dec' = 'blue',
                               'exc' = 'gray')) +
  geom_vline(aes(xintercept=30), color="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept = -30), color="black", linetype="dashed", size=1) +
  theme_classic() + labs(x = "UACR change (%)", y = "Number of patients") +
  scale_x_continuous(limits = c(-250,400), breaks=c(-100, -30, 30, seq (100,350,by=50))) +
  theme(text=element_text(size=20)) + theme(legend.position = "none")

dev.copy(png,"histuacr2.png", width=600, height=600)
dev.off()

####Sankey plot for CKD advancement####

df <- glom.meta.ckd

#Binning eGFR into CKD stages

df$eGFRbin_BL <- cut(df$eGFR_BL, breaks = c(0,14,29,59,89,114), labels = c("G5", "G4", "G3", "G2", "G1"))
df$eGFRbin_FU1y <- cut(df$eGFR_FU1y, breaks = c(0,14,29,59,89,114), labels = c("G5", "G4", "G3", "G2", "G1"))
df$eGFRbin_FU2y <- cut(df$eGFR_FU2y, breaks = c(0,14,29,59,89,114), labels = c("G5", "G4", "G3", "G2", "G1"))
df$eGFRbin_FU3y <- cut(df$eGFR_FU3y, breaks = c(0,14,29,59,89,114), labels = c("G5", "G4", "G3", "G2", "G1"))
df$eGFRbin_FU4y <- cut(df$eGFR_FU4y, breaks = c(0,14,29,59,89,114), labels = c("G5", "G4", "G3", "G2", "G1"))
df$eGFRbin_FU5y <- cut(df$eGFR_FU5y, breaks = c(0,14,29,59,89,114), labels = c("G5", "G4", "G3", "G2", "G1"))

allu.test <- df %>%
  make_long(eGFRbin_BL, eGFRbin_FU1y, eGFRbin_FU2y, eGFRbin_FU3y, eGFRbin_FU4y, eGFRbin_FU5y)

allu.test$node <- factor(allu.test$node, levels=c('G5', 'G4', 'G3', 'G2', 'G1'))
allu.test$next_node <- factor(allu.test$next_node, levels=c('G5', 'G4', 'G3', 'G2', 'G1'))

p1 <- ggplot(allu.test, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = c("green4", "green2", "yellow3", "orange2", "red"), na.value="white") +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("CKD stages over time") + theme(text=element_text(size=30))
p1 <- p1 + scale_x_discrete("Year", labels = c("eGFRbin_BL" = "0","eGFRbin_FU1y" = "1", "eGFRbin_FU2y" = "2","eGFRbin_FU3y" = "3","eGFRbin_FU4y" = "4", "eGFRbin_FU5y" = "5"))

dev.copy(png,"REDOckdsankey.png", width=1200, height=600)
dev.off()

####Sankey plot for NRA####

df$UACR_FU4y <- as.numeric(df$UACR_FU4y)

#Binning UACR into NRA or non-NRA

df$UACRbin_BL <- cut(df$UACR_BL, breaks = c(0,220,909), labels = c("non-NRA", "NRA"))
df$UACRbin_FU1y <- cut(df$UACR_FU1y, breaks = c(0,220,909), labels = c("non-NRA", "NRA"))
df$UACRbin_FU2y <- cut(df$UACR_FU2y, breaks = c(0,220,909), labels = c("non-NRA", "NRA"))
df$UACRbin_FU3y <- cut(df$UACR_FU3y, breaks = c(0,220,909), labels = c("non-NRA", "NRA"))
df$UACRbin_FU4y <- cut(df$UACR_FU4y, breaks = c(0,220,909), labels = c("non-NRA", "NRA"))
df$UACRbin_FU5y <- cut(df$UACR_FU5y, breaks = c(0,220,909), labels = c("non-NRA", "NRA"))

allu.uacr <- df %>%
  make_long(UACRbin_BL, UACRbin_FU1y, UACRbin_FU2y, UACRbin_FU3y, UACRbin_FU4y, UACRbin_FU5y)

allu.uacr$node <- factor(allu.uacr$node, levels=c('non-NRA', 'NRA'))
allu.uacr$next_node <- factor(allu.uacr$next_node, levels=c('non-NRA', 'NRA'))

p2 <- ggplot(allu.uacr, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = c("green3", "red"), na.value="white") +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Albuminuria over time") + theme(text=element_text(size=30))
p2 <- p2 + scale_x_discrete("Year", labels = c("UACRbin_BL" = "0","UACRbin_FU1y" = "1", "UACRbin_FU2y" = "2","UACRbin_FU3y" = "3","UACRbin_FU4y" = "4", "UACRbin_FU5y" = "5"))

dev.copy(png,"REDOnrasankey.png", width=1200, height=600)
dev.off()
