library(openxlsx)
library(tidyr)
library(dplyr)
library(mratios)

dane <- read.xlsx('TableS13_expression_all_data.xlsx')
wb <- loadWorkbook("TableS14_expression_means.xlsx")
wb <- loadWorkbook("~/Pobrane/TableS14_expression_means_n2.xlsx")

term_list <- unique(dane$Term)

#################
### col1 #######
################

################
#### Gen 1 #####
################

k=0*26
genes = "FTa1a"

po_list <- c(term_list[1:5], tail(term_list,1))

pval <- c()
for (n in 1:5) {
  term = po_list[n+1]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:11]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
  if (vt$p.value < ){
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c(term_list[1:5], tail(term_list,1))

pval <- c()
for (n in 1:5) {
  term = po_list[n+1]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:11]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c(term_list[1:5], tail(term_list,1))

pval <- c()
for (n in 1:5) {
  term = po_list[n+1]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:11]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c(term_list[1:5], tail(term_list,1))

pval <- c()
for (n in 1:5) {
  term = po_list[n+1]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:11]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 20, startRow = 22+k, colNames = F, rowNames = F)


#################
### col2 #######
################

################
#### Gen 1 #####
################

k=0*26
genes = "FTa1a"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+1]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25, startRow = 22+k, colNames = F, rowNames = F)

#################
### col3 #######
################

################
#### Gen 1 #####
################

k=0*26
l=1*5
genes = "FTa1a"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

#################
### col4 #######
################

################
#### Gen 1 #####
################

k=0*26
l=2*5
genes = "FTa1a"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- term_list[1:5]

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- term_list[6:10]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- term_list[16:20]

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- term_list[12:15]

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

#################
### col5 #######
################

################
#### Gen 1 #####
################

k=0*26
l=3*5
genes = "FTa1a"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

#################
### col6 #######
################

################
#### Gen 1 #####
################

k=0*26
l=4*5
genes = "FTa1a"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:5) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

#################
### col7 #######
################

################
#### Gen 1 #####
################

k=0*26
l=5*5
genes = "FTa1a"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

#################
### col8 #######
################

################
#### Gen 1 #####
################

k=0*26
l=6*5
genes = "FTa1a"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c(term_list[16:20],term_list[1:5])

pval <- c()
for (n in 1:4) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+5]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[6:11])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15],term_list[1:5])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 14+k, colNames = F, rowNames = F)

po_list <- c(term_list[12:15], term_list[16:20])

pval <- c()
for (n in 1:3) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+4]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 19+k, colNames = F, rowNames = F)

#################
### col9 #######
################

################
#### Gen 1 #####
################

k=0*26
l=7*5
genes = "FTa1a"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  5","Pa  6","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  4","Wo  5","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  5","Pa  6","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  4","Wo  5","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  5","Pa  6","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  4","Wo  5","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  5","Pa  6","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  4","Wo  5","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

#################
### col10 #######
################

################
#### Gen 1 #####
################

k=0*26
l=8*5
genes = "FTa1a"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  3","PR  4","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="8-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

#################
### col11 #######
################

################
#### Gen 1 #####
################

k=0*26
l=9*5
genes = "FTa1a"

po_list <- c("Po  5","Po 6","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  4","Pa  5","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c("Po  5","Po 6","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  4","Pa  5","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c("Po  5","Po 6","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  4","Pa  5","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c("Po  5","Po 6","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  4","Pa  5","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="non-vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

#################
### col12 #######
################

################
#### Gen 1 #####
################

k=0*26
l=10*5
genes = "FTa1a"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 2 #####
################

k=1*26
genes = "FTa1b"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 3 #####
################

k=2*26
genes = "FTc1"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

################
#### Gen 4 #####
################

k=3*26
genes = "FTc2"

po_list <- c("Po  4","Po  5","Po  1","Po  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 2+k, colNames = F, rowNames = F)

po_list <- c("Pa  3","Pa  4","Pa  1","Pa  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 9+k, colNames = F, rowNames = F)

po_list <- c("Wo  3","Wo  4","Wo  1","Wo  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 16+k, colNames = F, rowNames = F)

po_list <- c("PR  2","PR  3","PR  1","PR  1")

pval <- c()
for (n in 1:2) {
  term = po_list[n]
  df1 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes)
  term = po_list[n+2]
  df2 <- dane %>% filter(Photoperiod=="16-hour") %>% filter(Term==term) %>% 
    filter(Vernalization.variant=="vernalized") %>% filter(Gene==genes) 
  
  m1 <- mean(df1$Normalized.expression)
  m2 <- mean(df2$Normalized.expression)
  
  ratio <- m1/m2
  test <- ttestratio(df1$Normalized.expression,df2$Normalized.expression)
  print(term)
  print(ratio)
  print(test$p.value)
  vt <- var.test(df1$Normalized.expression, df2$Normalized.expression)
  
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
    pvag <- paste(round(test$p.value,3),gwiazdki,sep = " ")
    pval <- c(pval,pvag)
  } else {
    pval <- c(pval,"NS")
  }
}

dfwrite <- data.frame(pval=pval)

writeData(wb, sheet = 1, dfwrite, startCol = 25+l, startRow = 22+k, colNames = F, rowNames = F)

saveWorkbook(wb,"TableS14_expression_means_n.xlsx",overwrite = T)

