library(openxlsx)
library(tidyr)
library(dplyr)
library(mratios)


#################
###   8h  #######
################


dane <- read.xlsx('cq_reference genes.xlsx',sheet = '8H')
wb <- loadWorkbook("Supplementary_Tables_S1_S17.xlsx")

gen_list <- unique(dane$Target)
line_list <- unique(dane$Line)

pval <- c()
for (k in gen_list) {
  df1 <- dane %>% filter(Vern=="NWER") %>%
    filter(Target==k)
  df2 <- dane %>% filter(Vern=="WER") %>%
    filter(Target==k) 
  
  m1 <- mean(df1$Mean.Efficiency.Corrected.Cq)
  m2 <- mean(df2$Mean.Efficiency.Corrected.Cq)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Mean.Efficiency.Corrected.Cq,df2$Mean.Efficiency.Corrected.Cq)
  print(paste(k,k,sep = ':'))
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Mean.Efficiency.Corrected.Cq, df2$Mean.Efficiency.Corrected.Cq)
  
  if (vt$p.value < 0.0001){
    pval <- c(pval,"-")
    next
  }
  if (test$p.value< 0.001){
    gwiazdki <- "***"
  } else if (test$p.value < 0.01){
    gwiazdki <- "**"
  } else if (test$p.value < 0.05) {
    gwiazdki <- "*"
  } else {
    gwiazdki <- "NS"
  }
  if (test$p.value<0.05){
    pvag <- paste(round(test$p.value,2),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pvag <- paste(round(test$p.value,2),'NS',sep = " ") 
    pval <- c(pval,pvag)
  }
  for (n in line_list){
    df1 <- dane %>% filter(Vern=="NWER") %>% filter(Line==n) %>% 
      filter(Target==k)
    df2 <- dane %>% filter(Vern=="WER") %>% filter(Line==n) %>% 
      filter(Target==k) 
    
    m1 <- mean(df1$Mean.Efficiency.Corrected.Cq)
    m2 <- mean(df2$Mean.Efficiency.Corrected.Cq)
    
    ratio <- m1/m2
    test <- ttestratio(df1$Mean.Efficiency.Corrected.Cq,df2$Mean.Efficiency.Corrected.Cq)
    print(paste(n,k,sep = ':'))
    print(ratio)
    print(test$p.value)
    vt <- var.test(df1$Mean.Efficiency.Corrected.Cq, df2$Mean.Efficiency.Corrected.Cq)
    
    if (vt$p.value < 0.0001){
      pval <- c(pval,"-")
      next
    }
    if (test$p.value< 0.001){
      gwiazdki <- "***"
    } else if (test$p.value < 0.01){
      gwiazdki <- "**"
    } else if (test$p.value < 0.05) {
      gwiazdki <- "*"
    } else {
      gwiazdki <- "NS"
    }
    if (test$p.value<0.05){
      pvag <- paste(round(test$p.value,2),gwiazdki,sep = " ")
      pval <- c(pval,pvag)
    } else {
      pvag <- paste(round(test$p.value,2),'NS',sep = " ") 
      pval <- c(pval,pvag)
    }
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 17, dfwrite, startCol = 9, startRow = 16, colNames = F, rowNames = F)

#################
###   16h #######
################


dane <- read.xlsx('cq_reference genes.xlsx',sheet = '16H')

gen_list <- unique(dane$Target)
line_list <- unique(dane$Line)

pval <- c()
for (k in gen_list) {
  df1 <- dane %>% filter(Vern=="NWER") %>%
    filter(Target==k)
  df2 <- dane %>% filter(Vern=="WER") %>%
    filter(Target==k) 
  
  m1 <- mean(df1$Mean.Efficiency.Corrected.Cq)
  m2 <- mean(df2$Mean.Efficiency.Corrected.Cq)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Mean.Efficiency.Corrected.Cq,df2$Mean.Efficiency.Corrected.Cq)
  print(paste(k,k,sep = ':'))
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Mean.Efficiency.Corrected.Cq, df2$Mean.Efficiency.Corrected.Cq)
  
  if (vt$p.value < 0.0001){
    pval <- c(pval,"-")
    next
  }
  if (test$p.value< 0.001){
    gwiazdki <- "***"
  } else if (test$p.value < 0.01){
    gwiazdki <- "**"
  } else if (test$p.value < 0.05) {
    gwiazdki <- "*"
  } else {
    gwiazdki <- "NS"
  }
  if (test$p.value<0.05){
    pvag <- paste(round(test$p.value,2),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pvag <- paste(round(test$p.value,2),'NS',sep = " ") 
    pval <- c(pval,pvag)
  }
  for (n in line_list){
    df1 <- dane %>% filter(Vern=="NWER") %>% filter(Line==n) %>% 
      filter(Target==k)
    df2 <- dane %>% filter(Vern=="WER") %>% filter(Line==n) %>% 
      filter(Target==k) 
    
    m1 <- mean(df1$Mean.Efficiency.Corrected.Cq)
    m2 <- mean(df2$Mean.Efficiency.Corrected.Cq)
    
    ratio <- m1/m2
    test <- ttestratio(df1$Mean.Efficiency.Corrected.Cq,df2$Mean.Efficiency.Corrected.Cq)
    print(paste(n,k,sep = ':'))
    print(ratio)
    print(test$p.value)
    vt <- var.test(df1$Mean.Efficiency.Corrected.Cq, df2$Mean.Efficiency.Corrected.Cq)
    
    if (vt$p.value < 0.0001){
      pval <- c(pval,"-")
      next
    }
    if (test$p.value< 0.001){
      gwiazdki <- "***"
    } else if (test$p.value < 0.01){
      gwiazdki <- "**"
    } else if (test$p.value < 0.05) {
      gwiazdki <- "*"
    } else {
      gwiazdki <- "NS"
    }
    if (test$p.value<0.05){
      pvag <- paste(round(test$p.value,2),gwiazdki,sep = " ")
      pval <- c(pval,pvag)
    } else {
      pvag <- paste(round(test$p.value,2),'NS',sep = " ") 
      pval <- c(pval,pvag)
    }
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 17, dfwrite, startCol = 9, startRow = 6, colNames = F, rowNames = F)

saveWorkbook(wb,"Supplementary_Tables_S1_S17_1.xlsx",overwrite = T)
