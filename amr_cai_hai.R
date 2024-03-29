setwd ("~/Documents/AMR_model/amr_urine")

library(corrplot) 
library(ggcorrplot)
library(GGally)
library(readr)
library(tseries)
set.seed(123)  # For reproducibility


# st <- c("Enterococcus", "Escherichiacoli", "Klebsiella", "Pseudomonasaeruginosa", "Acinetobacter", "Staphylococcus")

st <- c("Escherichiacoli")

########### data paths ############
main_path <- "data/time_series_processed/inf/"
save_path <- "results/CAI_HAI_lead_lag/"

check_stationary <- function(series) {
  # print(" ********* Inside check_stationary")
  series_p <- adf.test(series)$p.value
  return(series_p <= 0.05)
}

my_diff <- function(x){x = c(diff(x))}
temp_ccf<- NULL
temp_ccf_pvalue <- NULL

make_stationary_and_calc_ccf <- function(series1, series2, anitibiotic_name) {
  # print("Call to make_stationary")
  stat_series1 <- series1
  stat_series2 <- series2
  count<- 0
  
  if(length(stat_series1) != length(stat_series2)) {
    print("------ Diff Num rows detected -------------")
    print(length(stat_series1))
    print(length(stat_series2))
    sink()
    stop("Different number of rows detected")
  }

  while(check_stationary(stat_series1) == FALSE || check_stationary(stat_series2) == FALSE) {
    stat_series1<- sapply(data.frame(stat_series1),my_diff)
    stat_series2<- sapply(data.frame(stat_series2),my_diff)
    # print("---- Pass inside make_stationary")
    count<- count+1
    cat(paste0(antibiotic_name, ",", count,"\n"))
  }
  
  ccf_values <- ccf(as.vector(stat_series1), 
                    as.vector(stat_series2), 
                    main = paste0(anitibiotic_name,"_cai", " ", anitibiotic_name, "_hai"))
  
  ps <- 2 *(1 - pnorm(abs(ccf_values$acf), mean = 0, sd = 1/sqrt(ccf_values$n.used)))
  
  # Construct a data frame containing lags, correlations, and p-values
  ccf_results <- data.frame(
    lag = ccf_values$lag,
    correlation = ccf_values$acf,
    p_value = ps
  )
  file_name_ccf <- paste0(anitibiotic_name,"_ecoli_cai_hai.csv")
  write.csv(ccf_results, file_name_ccf)
  
}

for(nam in st) {
  sink(paste0(save_path, nam, ".csv"))
  ####### data loading #########
  dat_cai <- read.csv(paste0(main_path, paste0("month_year_", nam, "_CAI.csv")))
  dat_hai <- read.csv(paste0(main_path, paste0("month_year_", nam, "_HAI.csv")))
  dat_cai$date<- NULL
  dat_hai$date<- NULL
  
  pdf(paste0("./results/CAI_HAI_lead_lag/", nam, ".pdf"), height = 6, width = 8)
  par(mar =c(10,4,5,2))
  cai_colname <- colnames(dat_cai)
  hai_colname <- colnames(dat_hai)
  for(i in 1:length(cai_colname)){
    for(j in 1:length(hai_colname)){
      if(cai_colname[i] == hai_colname[j]) {
        #print(paste0("Processing column - ", cai_colname[i]))
        cai_series <- dat_cai[,i]
        hai_series <- dat_hai[,j]
        antibiotic_name <- cai_colname[i]
        make_stationary_and_calc_ccf(cai_series, hai_series, antibiotic_name)
        break
      }
    }
    # break
  }
  sink()
  dev.off()
}

