"""This is a common script used for creating heatmaps for the complete data as well as multidrug resistant data for which different data files are used"""

#Importing libraries
library(readr)
library(plyr)
library("gplots")
library(reshape2)

set.seed(100)

#Data Files for all records
amr_data <- read_csv(file="data/final/same_zone.csv")

#Filtering organism-wise data
filter_data <- function(df, organism) {
  df <- df[df$organism_name == organism, ]

  if (organism == "Enterococcus") {
    antibiotics <- c("Vancomycin", "Nitrofurantoin", "Linezolid", "Teicoplanin", "Ampicillin", "Gentamicin HL", "Ciprofloxacin", "Fosfomycin")
  } else if (organism == "Escherichia coli ") {
    antibiotics <- c("Amikacin", "Piperacillin-tazobactam", "Nitrofurantoin", "Imipenem", "Meropenem", "Ciprofloxacin", "Trimethoprim-sulfamethoxazole", "Cefotaxime", "Ertapenem", "Levofloxacin", "Fosfomycin", "Cefazolin", "Colistin")
  } else if (organism == "Klebsiella") {
    antibiotics <- c("Amikacin", "Piperacillin-tazobactam", "Nitrofurantoin", "Imipenem", "Meropenem", "Ciprofloxacin", "Trimethoprim-sulfamethoxazole", "Ertapenem", "Cefotaxime", "Levofloxacin", "Fosfomycin", "Cefazolin", "Colistin")
  } else if (organism == "Pseudomonas aeruginosa") {
    antibiotics <- c("Amikacin", "Meropenem", "Piperacillin-tazobactam", "Ciprofloxacin", "Imipenem", "Cefepime", "Ceftazidime", "Gentamicin", "Levofloxacin", "Tobramycin", "Colistin")
  }

  list(df = df, antibiotics = antibiotics)
}

organism = "Enterococcus"
amr_data <- filter_data(amr_data, organism)

summary(amr_data)

# Splitting data by 'zone'.
amr_data_zone <- dlply(.data = amr_data, .variables = "zone")
# Counting the number of rows in each zone and creating a dataframe.
zone_count <- sapply(amr_data_zone,nrow)
zone_df <- data.frame(zone=names(zone_count), count=zone_count) 
# Function to summarize AMR data for each antibiotic.
getSummaryAMR <- function(xx) {
  colNames <- antibiotics
  summary_amr <- list()
  print(xx)
  for(i in 1:length(colNames)){
    print(colNames[i])
    ytab <- table(xx[,colNames[i]])
    if(length(ytab) == 0){
      df_ <- data.frame(Var1=c('I','R','S'),Freq=c(0,0,0))
    }else{
      df_ <- data.frame(ytab)
    }
    sum_count <- sum(df_$Freq)
    if(sum_count > 0) {
      df_ <- transform(df_, Freq = Freq/sum_count) 
    }
    print(df_)
    df_$antibiotic <- colNames[i]
    summary_amr[[i]] <- df_
  }
  summary_amr_df <- ldply(summary_amr)
  return(summary_amr_df)
}

# Applying the summary function to each zone.
amr_summary <- lapply(amr_data_zone, getSummaryAMR)

names(amr_summary) <- names(amr_data_zone)
# Combining the summaries into a single dataframe.
amr_summary_df <- ldply(amr_summary,.id = 'zone')
# Filtering out rows with frequency 0.
amr_summary_df <- amr_summary_df[amr_summary_df$Freq != 0,]

# Preparing data for heatmap.
amr_summary_df_patter <- dlply(amr_summary_df,'Var1')

xx <- dcast(amr_summary_df_patter[[2]], formula = zone~antibiotic,value.var = "Freq")
xx_mat <-as.matrix(xx[,-1])
rownames(xx_mat) <- xx$zone
xx_mat[is.na(xx_mat)] <- 0
xx_mat <- xx_mat *100

# Creating the heatmap.
col <- rev(heat.colors(999))

hm <- heatmap.2(xx_mat, trace="none",cexCol=0.9, cexRow=0.9, cellnote=round(xx_mat), notecol="black", density.info="none", margins=c(12,10),  col=col, scale = "column")

