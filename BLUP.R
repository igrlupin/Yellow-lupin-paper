library(xlsx)
library(openxlsx)
library(tidyr)
library(metan)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(agricolae)
library(rstatix)
library(lmerTest)
library(car)
library(moments)
library(MVN)
library(emmeans)



df <- read.csv('TableS7.csv', sep = '\t')
df <- read.csv('TableS7_updated.csv', sep = '\t')
df$Name[df$Name == "Parys "] <- "Parys"

df1 <- df[,1:15]


BE_full <- df1 %>% pivot_longer(cols = starts_with("BE"),names_to = "Plant", 
                            values_to = "BE", values_drop_na = TRUE)

BE <- BE_full %>% filter(Variant == 'non-vernalized')
  
vBE <- BE_full %>% filter(Variant == 'vernalized')

BE_model <- gamem_met(BE, env = Year, gen = Name, rep = Plant, random = 'gen', resp = everything())
BE_data <- get_model_data(BE_model)
print(BE_data)
BE_data_d <- get_model_data(BE_model, "details")
print(BE_data_d)
BE_blupmod <- get_model_data(BE_model, what = "blupg")

vBE_model <- gamem_met(vBE, env = Year, gen = Name, rep = Plant, random = 'gen', resp = everything())
vBE_data <- get_model_data(vBE_model)
print(vBE_data)
vBE_blupmod <- get_model_data(vBE_model, what = "blupg")

colnames(vBE_data)[2] <- "vBE"

df2 <- df[,c(1:5,16:25)]


SF_full <- df2 %>% pivot_longer(cols = starts_with("SF"),names_to = "Plant", 
                                values_to = "SF", values_drop_na = TRUE)

SF <- SF_full %>% filter(Variant == 'non-vernalized')

vSF <- SF_full %>% filter(Variant == 'vernalized')

SF_model <- gamem_met(SF, env = Year, gen = Name, rep = Plant, random = 'gen', resp = everything())
SF_data <- get_model_data(SF_model)
print(SF_data)
SF_blupmod <- get_model_data(SF_model, what = "blupg")

vSF_model <- gamem_met(vSF, env = Year, gen = Name, rep = Plant, random = 'gen', resp = everything())
vSF_data <- get_model_data(vSF_model)
print(vSF_data)
vSF_blupmod <- get_model_data(vSF_model, what = "blupg")

colnames(vSF_data)[2] <- "vSF"

df3 <- df[,c(1:5,26:35)]


EF_full <- df3 %>% pivot_longer(cols = starts_with("EF"),names_to = "Plant", 
                                values_to = "EF", values_drop_na = TRUE)
EF <- EF_full %>% filter(Variant == 'non-vernalized')

vEF <- EF_full %>% filter(Variant == 'vernalized')

EF_model <- gamem_met(EF, env = Year, gen = Name, rep = Plant, random = 'gen', resp = everything())
EF_data <- get_model_data(EF_model)
print(EF_data)
EF_blupmod <- get_model_data(EF_model, what = "blupg")

vEF_model <- gamem_met(vEF, env = Year, gen = Name, rep = Plant, random = 'gen', resp = everything())
vEF_data <- get_model_data(vEF_model)
print(vEF_data)
vEF_blupmod <- get_model_data(vEF_model, what = "blupg")

colnames(vEF_data)[2] <- "vEF"

df_stab9 <- cbind(BE_data,vBE_data[,2], SF_data[,2], vSF_data[,2], EF_data[,2], vEF_data[,2])

df_stab9[,2:7] <- round(df_stab9[,2:7],4)
df_stab9[c(2,4),2:7] <- df_stab9[c(2,4),2:7]*100
df_stab9[,2:7] <- round(df_stab9[,2:7],2)

write.xlsx(df_stab9, '~/konas13@gmail.com/Dokumenty/Łubin/Publikacje/FT_zolty_Poznań/TableS9_update.xlsx', 
           sheetName = "stab9", col.names = TRUE, row.names = FALSE, append = FALSE)


#########################
# Korelacja Geno_Feno ###
#########################

dane_geno <- read.csv('geno_nowe.csv', sep = '\t', dec = ',')
dane_geno <- dane_geno[order(factor(dane_geno$Marker.name, levels = levels(EF_blupmod$GEN))),]

all(dane_geno$Name == BE_blupmod$GEN)

ss <- vBE_blupmod$BE - BE_blupmod$BE

pval.vec <- c()
pval.vec.ver <- c()
cor.vec <- c()
cor.vec.ver <- c()
for (i in 2:22){
  df.cor <- data.frame(ver=BE_blupmod$BE, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec <- c(pval.vec, unname(cor.dat$p.value))
  cor.vec <- c(cor.vec, unname(cor.dat$estimate))
  df.cor <- data.frame(ver=ss, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec.ver <- c(pval.vec.ver, unname(cor.dat$p.value))
  cor.vec.ver <- c(cor.vec.ver, unname(cor.dat$estimate))
}

df.cor.res.BE <- data.frame(BE=cor.vec,BE.pval=pval.vec,vBE=cor.vec.ver,vBE.pval=pval.vec.ver)

ss <- vSF_blupmod$SF - SF_blupmod$SF

pval.vec <- c()
pval.vec.ver <- c()
cor.vec <- c()
cor.vec.ver <- c()
for (i in 2:22){
  df.cor <- data.frame(ver=SF_blupmod$SF, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec <- c(pval.vec, unname(cor.dat$p.value))
  cor.vec <- c(cor.vec, unname(cor.dat$estimate))
  df.cor <- data.frame(ver=ss, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec.ver <- c(pval.vec.ver, unname(cor.dat$p.value))
  cor.vec.ver <- c(cor.vec.ver, unname(cor.dat$estimate))
}

df.cor.res.SF <- data.frame(SF=cor.vec,SF.pval=pval.vec,vSF=cor.vec.ver,vSF.pval=pval.vec.ver)

ss <- vEF_blupmod$EF - EF_blupmod$EF

pval.vec <- c()
pval.vec.ver <- c()
cor.vec <- c()
cor.vec.ver <- c()
for (i in 2:22){
  df.cor <- data.frame(ver=EF_blupmod$EF, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec <- c(pval.vec, unname(cor.dat$p.value))
  cor.vec <- c(cor.vec, unname(cor.dat$estimate))
  df.cor <- data.frame(ver=ss, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec.ver <- c(pval.vec.ver, unname(cor.dat$p.value))
  cor.vec.ver <- c(cor.vec.ver, unname(cor.dat$estimate))
}

df.cor.res.EF <- data.frame(EF=cor.vec,EF.pval=pval.vec,vEF=cor.vec.ver,vEF.pval=pval.vec.ver)

df.cor.res <- cbind(df.cor.res.BE,df.cor.res.SF, df.cor.res.EF)
rownames(df.cor.res) <- colnames(dane_geno[,2:22])
df.cor.res.xlsx <- t(df.cor.res)
mycolumn <- seq(1,12,2)
df.cor.res.heatmap <- df.cor.res[,mycolumn]
rownames(df.cor.res.heatmap) <- colnames(dane_geno[,2:22])
df.cor.res.heatmap <- df.cor.res.heatmap %>% filter(rownames(df.cor.res.heatmap)!='FTa1a_F12_R12')

svg("fig_heatmapa.svg", width = 8, height = 6, pointsize = 12)
Heatmap(as.matrix(df.cor.res.heatmap), cluster_columns=FALSE,
        row_names_side = "left",
        row_dend_sid = "left",
        row_names_gp=gpar(cex=1.4),
        row_km = 2,
        heatmap_legend_param = list(title = "cor"),
        column_names_gp = gpar(cex = 1.6))
dev.off()


df.cor.res <- df.cor.res %>% filter(rownames(df.cor.res)!='FTa1a_F12_R12')
df.cor.res$marker.name <- rownames(df.cor.res)

wb <- createWorkbook()
addWorksheet(wb, "stab_n")
posStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
writeData(wb, "stab_n", df.cor.res)
conditionalFormatting(wb, "stab_n",
                      cols = 1,
                      rows = 2:21, rule = "B2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 2,
                      rows = 2:21, rule = "B2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 3,
                      rows = 2:21, rule = "D2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 4,
                      rows = 2:21, rule = "D2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 5,
                      rows = 2:21, rule = "F2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 6,
                      rows = 2:21, rule = "F2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 7,
                      rows = 2:21, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 8,
                      rows = 2:21, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 9,
                      rows = 2:21, rule = "J2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 10,
                      rows = 2:21, rule = "J2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 11,
                      rows = 2:21, rule = "L2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 12,
                      rows = 2:21, rule = "L2<0.05", style = posStyle
)

headerStyle <- createStyle(
  fontSize = 14, fontColour = "black", halign = "center",
  fgFill = "white", border = "TopBottom", borderColour = "black",
  textDecoration = "bold"
)

addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:13, gridExpand = TRUE)

saveWorkbook(wb, "stab_n_update.xlsx", TRUE)

#############################
####### korelacja średnie ###
#############################


BE_full$Year <- as.factor(BE_full$Year)
BE_full$Name <- as.factor(BE_full$Name)
BE_full$Variant <- as.factor(BE_full$Variant)

BE_mean <- BE_full %>% group_by(Name, Variant) %>% summarise(mean(BE))
BE_mean1 <- BE_mean[BE_mean$Variant=='non-vernalized',]
BE_mean2 <- BE_mean[BE_mean$Variant=='vernalized',]
all(dane_geno$Name == BE_mean1$Name)

pval.vec <- c()
pval.vec.ver <- c()
cor.vec <- c()
cor.vec.ver <- c()
for (i in 4:22){
  df.cor <- data.frame(ver=BE_mean1$`mean(BE)`, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec <- c(pval.vec, unname(cor.dat$p.value))
  cor.vec <- c(cor.vec, unname(cor.dat$estimate))
  df.cor <- data.frame(ver=(BE_mean2$`mean(BE)`- BE_mean1$`mean(BE)`), allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec.ver <- c(pval.vec.ver, unname(cor.dat$p.value))
  cor.vec.ver <- c(cor.vec.ver, unname(cor.dat$estimate))
}

df.cor.res.BE <- data.frame(BE=cor.vec,BE.pval=pval.vec,vBE=cor.vec.ver,vBE.pval=pval.vec.ver)

SF_full$Year <- as.factor(SF_full$Year)
SF_full$Name <- as.factor(SF_full$Name)
SF_full$Variant <- as.factor(SF_full$Variant)

SF_mean <- SF_full %>% group_by(Name, Variant) %>% summarise(mean(SF))
SF_mean1 <- SF_mean[SF_mean$Variant=='non-vernalized',]
SF_mean2 <- SF_mean[SF_mean$Variant=='vernalized',]
all(dane_geno$Name == SF_mean1$Name)

pval.vec <- c()
pval.vec.ver <- c()
cor.vec <- c()
cor.vec.ver <- c()
for (i in 4:22){
  df.cor <- data.frame(ver=SF_mean1$`mean(SF)`, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec <- c(pval.vec, unname(cor.dat$p.value))
  cor.vec <- c(cor.vec, unname(cor.dat$estimate))
  df.cor <- data.frame(ver=(SF_mean2$`mean(SF)`- SF_mean1$`mean(SF)`), allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec.ver <- c(pval.vec.ver, unname(cor.dat$p.value))
  cor.vec.ver <- c(cor.vec.ver, unname(cor.dat$estimate))
}

df.cor.res.SF <- data.frame(SF=cor.vec,SF.pval=pval.vec,vSF=cor.vec.ver,vSF.pval=pval.vec.ver)

EF_full$Year <- as.factor(EF_full$Year)
EF_full$Name <- as.factor(EF_full$Name)
EF_full$Variant <- as.factor(EF_full$Variant)

EF_mean <- EF_full %>% group_by(Name, Variant) %>% summarise(mean(EF))
EF_mean1 <- EF_mean[EF_mean$Variant=='non-vernalized',]
EF_mean2 <- EF_mean[EF_mean$Variant=='vernalized',]
all(dane_geno$Name == EF_mean1$Name)

pval.vec <- c()
pval.vec.ver <- c()
cor.vec <- c()
cor.vec.ver <- c()
for (i in 4:22){
  df.cor <- data.frame(ver=EF_mean1$`mean(EF)`, allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec <- c(pval.vec, unname(cor.dat$p.value))
  cor.vec <- c(cor.vec, unname(cor.dat$estimate))
  df.cor <- data.frame(ver=(EF_mean2$`mean(EF)`- EF_mean1$`mean(EF)`), allele = dane_geno[,i])
  df.cor <- df.cor %>% drop_na(allele) %>% drop_na(ver)
  df.cor <- df.cor[order(df.cor$allele),]
  cor.dat <- cor.test(df.cor$ver,df.cor$allele,method = 'spearman')
  pval.vec.ver <- c(pval.vec.ver, unname(cor.dat$p.value))
  cor.vec.ver <- c(cor.vec.ver, unname(cor.dat$estimate))
}

df.cor.res.EF <- data.frame(EF=cor.vec,EF.pval=pval.vec,vEF=cor.vec.ver,vEF.pval=pval.vec.ver)

df.cor.res <- cbind(df.cor.res.BE,df.cor.res.SF,df.cor.res.EF)
rownames(df.cor.res) <- colnames(dane_geno[,4:22])
df.cor.res.xlsx <- t(df.cor.res)
mycolumn <- seq(1,12,2)
df.cor.res.heatmap <- df.cor.res[,mycolumn]
rownames(df.cor.res.heatmap) <- colnames(dane_geno[,4:22])

Heatmap(df.cor.res.heatmap, cluster_columns=FALSE,
        row_names_side = "left",
        row_dend_sid = "left",
        row_names_gp=gpar(cex=0.6))

df.cor.res$marker.name <- rownames(df.cor.res)

wb <- createWorkbook()
addWorksheet(wb, "stab_n")
posStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
writeData(wb, "stab_n", df.cor.res)
conditionalFormatting(wb, "stab_n",
                      cols = 1,
                      rows = 2:20, rule = "B2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 2,
                      rows = 2:20, rule = "B2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 3,
                      rows = 2:20, rule = "D2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 4,
                      rows = 2:20, rule = "D2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 5,
                      rows = 2:20, rule = "F2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 6,
                      rows = 2:20, rule = "F2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 7,
                      rows = 2:20, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 8,
                      rows = 2:20, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 9,
                      rows = 2:20, rule = "J2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 10,
                      rows = 2:20, rule = "J2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 11,
                      rows = 2:20, rule = "L2<0.05", style = posStyle
)
conditionalFormatting(wb, "stab_n",
                      cols = 12,
                      rows = 2:20, rule = "L2<0.05", style = posStyle
)

headerStyle <- createStyle(
  fontSize = 14, fontColour = "black", halign = "center",
  fgFill = "white", border = "TopBottom", borderColour = "black",
  textDecoration = "bold"
)

addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:13, gridExpand = TRUE)

saveWorkbook(wb, "stab_nn.xlsx", TRUE)


################################
### Wilxon test     ############
################################

pwc <- BE_full %>% group_by(Year,Name)
pwc2016 <- BE_full %>% filter(Year=='2016') %>% filter(Name!='Talar (WTD 2504)') %>% group_by(Name, Variant) %>% summarise(mean(BE), sd(BE)) %>% 
  wilcox_test(BE~Variant, p.adjust.method = "bonferroni")
pwc2017 <- BE_full %>% filter(Year=='2017') %>% group_by(Name) %>% wilcox_test(BE~Variant, p.adjust.method = "bonferroni")
pwc2019 <- BE_full %>% filter(Year=='2019') %>% group_by(Name, Variant) %>% summarise(mean(BE))# %>% wilcox_test(BE~Variant, p.adjust.method = "bonferroni")

ss <- BE_full %>% filter(Year=='2016') %>% filter(Variant=='non-vernalized') %>% pull(Name)
sss <- BE_full %>% filter(Year=='2016') %>% filter(Variant=='vernalized') %>% pull(Name)

length(unique(pwc2016$Name)) == length(unique(pwc2017$Name))
setdiff(unique(pwc2017$Name), unique(pwc2016$Name))

pwc2016a <- BE_full %>% filter(Name=='Talar (WTD 2504)') %>% filter(Year=='2016') 

pwc2016a %>% group_by(Year, Variant) %>% summarise(mean(BE), sd(BE))

pwct1 <- pwc2016a %>% filter(Variant=='non-vernalized')
pwct2 <- pwc2016a %>% filter(Variant=='vernalized')

wilcox.test(pwct1$BE,pwct2$BE)

pwc2016w <- pwc2016 %>% wilcox_test(BE~Variant, p.adjust.method = "bonferroni")

unique(pwc2019$Name) == unique(pwc1$Name)

summary(pwc2019$BE)

pwc1 <- BE_full %>% group_by(Name) %>% pairwise_t_test(BE~Variant, p.adjust.method = "bonferroni")

pwc2 <- SF_full %>% group_by(Name) %>% pairwise_t_test(SF~Variant, p.adjust.method = "bonferroni")

pwc3 <- EF_full %>% group_by(Name) %>% pairwise_t_test(EF~Variant, p.adjust.method = "bonferroni")

wb <- createWorkbook()
addWorksheet(wb, "BE")
writeData(wb, "BE", as.data.frame(pwc1))
addWorksheet(wb, "SF")
writeData(wb, "SF", as.data.frame(pwc2))
addWorksheet(wb, "EF")
writeData(wb, "EF", as.data.frame(pwc3))


addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:8, gridExpand = TRUE)
addStyle(wb, sheet = 2, headerStyle, rows = 1, cols = 1:8, gridExpand = TRUE)
addStyle(wb, sheet = 3, headerStyle, rows = 1, cols = 1:8, gridExpand = TRUE)

conditionalFormatting(wb, "BE",
                      cols = 1,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 2,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 3,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 4,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 5,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 6,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 7,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "BE",
                      cols = 8,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)

conditionalFormatting(wb, "SF",
                      cols = 1,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 2,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 3,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 4,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 5,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 6,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 7,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "SF",
                      cols = 8,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)

conditionalFormatting(wb, "EF",
                      cols = 1,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 2,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 3,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 4,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 5,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 6,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 7,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)
conditionalFormatting(wb, "EF",
                      cols = 8,
                      rows = 2:112, rule = "H2<0.05", style = posStyle
)

saveWorkbook(wb, "stab_8n.xlsx", TRUE)

## modele

BE_full$Year <- as.factor(BE_full$Year)

model_mieszany_full_BE <- lme4::lmer(BE ~ Variant*Name + (1|Year) + (1|Year:Plant) + (1|Plant), data = BE_full)

coompare1_BE <- emmeans(model_mieszany_full_BE, list(pairwise ~ Variant | Name), adjust = "tukey")
BE_compare <- summary(coompare1_BE)
BE_compare_df <- data.frame(Name = BE_compare$`pairwise differences of Variant | Name`$Name, pval = BE_compare$`pairwise differences of Variant | Name`$p.value)
BE_blupul <- BE_model$BE$BLUPgen
vBE_blupul <- vBE_model$BE$BLUPgen
BE_blupul <- BE_blupul[order(factor(BE_blupul$GEN, levels = levels(EF_blupmod$GEN))),]
vBE_blupul <- vBE_blupul[order(factor(vBE_blupul$GEN, levels = levels(EF_blupmod$GEN))),]
BE_compare_df <- BE_compare_df[order(factor(BE_compare_df_BE$Name, levels = levels(EF_blupmod$GEN))),]

BE_blup_pval <- cbind(BE_blupul[,c(2,5,6,7)],vBE_blupul[,c(5,6,7)],BE_compare_df[,2])
colnames(BE_blup_pval) <- c("Name", "BE","BE.LL","BE.UL","vBE","vBE.LL","vBE.UL","pval")
BE_blup_pval[,2:7] <- round(BE_blup_pval[,2:7],2)

SF_full$Year <- as.factor(SF_full$Year)

model_mieszany_full_SF <- lme4::lmer(SF ~ Variant*Name + (1|Year) + (1|Year:Plant) + (1|Plant), data = SF_full)

coompare1_SF <- emmeans(model_mieszany_full_SF, list(pairwise ~ Variant | Name), adjust = "tukey")
SF_compare <- summary(coompare1_SF)
SF_compare_df <- data.frame(Name = SF_compare$`pairwise differences of Variant | Name`$Name, pval = SF_compare$`pairwise differences of Variant | Name`$p.value)
SF_blupul <- SF_model$SF$BLUPgen
vSF_blupul <- vSF_model$SF$BLUPgen
SF_blupul <- SF_blupul[order(factor(SF_blupul$GEN, levels = levels(EF_blupmod$GEN))),]
vSF_blupul <- vSF_blupul[order(factor(vSF_blupul$GEN, levels = levels(EF_blupmod$GEN))),]
SF_compare_df <- SF_compare_df[order(factor(SF_compare_df$Name, levels = levels(EF_blupmod$GEN))),]

SF_blup_pval <- cbind(SF_blupul[,c(2,5,6,7)],vSF_blupul[,c(5,6,7)],SF_compare_df[,2])
colnames(SF_blup_pval) <- c("Name", "SF","SF.LL","SF.UL","vSF","vSF.LL","vSF.UL","pval")
SF_blup_pval[,2:7] <- round(SF_blup_pval[,2:7],2)

EF_full$Year <- as.factor(EF_full$Year)

model_mieszany_full_EF <- lme4::lmer(EF ~ Variant*Name + (1|Year) + (1|Year:Plant) + (1|Plant), data = EF_full)

coompare1_EF <- emmeans(model_mieszany_full_EF, list(pairwise ~ Variant | Name), adjust = "tukey")
EF_compare <- summary(coompare1_EF)
EF_compare_df <- data.frame(Name = EF_compare$`pairwise differences of Variant | Name`$Name, pval = EF_compare$`pairwise differences of Variant | Name`$p.value)
EF_blupul <- EF_model$EF$BLUPgen
vEF_blupul <- vEF_model$EF$BLUPgen
EF_blupul <- EF_blupul[order(factor(EF_blupul$GEN, levels = levels(EF_blupmod$GEN))),]
vEF_blupul <- vEF_blupul[order(factor(vEF_blupul$GEN, levels = levels(EF_blupmod$GEN))),]
EF_compare_df <- EF_compare_df[order(factor(EF_compare_df$Name, levels = levels(EF_blupmod$GEN))),]

EF_blup_pval <- cbind(EF_blupul[,c(2,5,6,7)],vEF_blupul[,c(5,6,7)],EF_compare_df[,2])
colnames(EF_blup_pval) <- c("Name", "EF","EF.LL","EF.UL","vEF","vEF.LL","vEF.UL","pval")
EF_blup_pval[,2:7] <- round(EF_blup_pval[,2:7],2)

wb <- createWorkbook()
addWorksheet(wb, "BE")
writeData(wb, "BE", as.data.frame(BE_blup_pval))
addWorksheet(wb, "SF")
writeData(wb, "SF", as.data.frame(SF_blup_pval))
addWorksheet(wb, "EF")
writeData(wb, "EF", as.data.frame(EF_blup_pval))
saveWorkbook(wb, "TableS8_update.xlsx", TRUE)

BE_srednie <- BE_full %>% group_by(Name, Variant) %>% summarise(mean.BE=mean(BE), sd.BE=sd(BE))
SF_srednie <- SF_full %>% group_by(Name, Variant) %>% summarise(mean.SF=mean(SF), sd.BE=sd(SF))
EF_srednie <- EF_full %>% group_by(Name, Variant) %>% summarise(mean.BE=mean(EF), sd.BE=sd(EF))

wb <- createWorkbook()
addWorksheet(wb, "BE")
writeData(wb, "BE", as.data.frame(BE_srednie))
addWorksheet(wb, "SF")
writeData(wb, "SF", as.data.frame(SF_srednie))
addWorksheet(wb, "EF")
writeData(wb, "EF", as.data.frame(EF_srednie))
saveWorkbook(wb, "TableSX_srednie.xlsx", TRUE)

car::qqPlot(residuals(model_mieszany_full))
stz <- residuals(model_mieszany_full)

por1 <- BE_compare_df[BE_compare_df$pval>0.05,1]
por1 <- as.vector(por1)
por2 <- pwc1 %>% filter(fdr>0.05) %>% pull(Name)


setdiff(por1,por2)
setdiff(por2,por1)
intersect(por1,por2)

gen <- "Sulfa"

my_df <- BE_full %>% filter(Name==gen)
ggplot(my_df, aes(x=Name, y=BE, fill=Variant)) + 
  geom_boxplot()

my_df %>% group_by(Variant) %>% summarise(srednia=mean(BE),odchyl=sd(BE))
BE_compare_df %>% filter(Name==gen)
pwc1 %>% filter(Name==gen)

my_df %>% pairwise_t_test(BE~Variant, p.adjust.method = "bonferroni")
pwc1 <- BE_full %>% group_by(Name) %>%  t_test(BE~Variant, p.adjust.method = "bonferroni")

pwc1 <- as.data.frame(pwc1)
pwc1$fdr <- p.adjust(pwc1$p, method = 'fdr')

BE_blupmod %>% filter(GEN==gen)
vBE_blupmod %>% filter(GEN==gen)

dane_geno %>% filter(Name%in%por1)
