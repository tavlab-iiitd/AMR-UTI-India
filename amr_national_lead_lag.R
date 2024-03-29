# Set the working directory to the specified path
setwd ("~/Documents/AMR_model/amr_urine")

set.seed(123)  # For reproducibility

# Load necessary libraries for plotting, data manipulation, and statistical analysis
library(corrplot) 
library(ggcorrplot)
library(GGally)
library(readr)
library(tseries)

#st <- c("Enterococcus", "Escherichiacoli", "Klebsiella", "Pseudomonasaeruginosa")
st <- c("Escherichiacoli")

# Define paths for data input and figure output
########### data paths ############
main_path <- "data/for_r/national/"
save_path <- "figures/national_lead_lag/2024_new_"

# Function to check if a time series is stationary using the Augmented Dickey-Fuller test
check_stationary <- function(series) {
  # print(" ********* Inside check_stationary")
  series_p <- adf.test(series)$p.value
  return(series_p <= 0.05)
}

# Function to differentiate a series (used to help achieve stationarity)
my_diff <- function(x){x = c(diff(x))}

# Initialize variables to store cross-correlation function (CCF) results
temp_ccf<- NULL
temp_ccf_pvalue <- NULL

# Function to make series stationary and calculate CCF between two antibiotic time series
make_stationary_and_calc_ccf <- function(series1, series2, anitibiotic_name1, anitibiotic_name2) {
  # print("Call to make_stationary")
  stat_series1 <- series1
  stat_series2 <- series2
  count<- 0
  
  # Check if the two series have the same length, exit if not
  if(length(stat_series1) != length(stat_series2)) {
    print("------ Diff Num rows detected -------------")
    print(length(stat_series1))
    print(length(stat_series2))
    exit(-1)
  }

  # Differentiate series until both are stationary
  while(check_stationary(stat_series1) == FALSE || check_stationary(stat_series2) == FALSE) {
    stat_series1<- sapply(data.frame(stat_series1),my_diff)
    stat_series2<- sapply(data.frame(stat_series2),my_diff)
    # print("---- Pass inside make_stationary")
    count<- count+1
  }
  print(paste0(anitibiotic_name1, " x ", anitibiotic_name2, "--------------", count))
  
  # Calculate cross-correlation function
  ccf_values <- ccf(as.vector(stat_series1), 
                    as.vector(stat_series2), 
                    main = paste0(anitibiotic_name1, " ", anitibiotic_name2))
  # Calculate p-values for the correlations
  ps <- 2 *(1 - pnorm(abs(ccf_values$acf), mean = 0, sd = 1/sqrt(ccf_values$n.used)))
  
  # Construct a data frame containing lags, correlations, and p-values
  ccf_results <- data.frame(
    lag = ccf_values$lag,
    correlation = ccf_values$acf,
    p_value = ps
  )
  file_name_ccf <- paste0(antibiotic_name1, paste0("_", antibiotic_name2, ".csv"))
  write.csv(ccf_results, file_name_ccf)
  
}
for(nam in st) {
  ####### data loading #########
  dat <- read.csv(paste0(main_path, paste0("month_year_", nam, ".csv")))
  dat$date<- NULL
  
  pdf(paste0(save_path, nam, ".pdf"), height = 6, width = 8)
  par(mar =c(10,4,5,2))
  dat_colnames <- colnames(dat)
  # Loop through pairs of antibiotics to process
  for(i in 1:(length(dat_colnames) - 1)){
    for(j in (i + 1):length(dat_colnames)){
      #print(paste0("Processing columns - ", dat_colnames[i], " x ", dat_colnames[j]))
      series1 <- dat[,i]
      series2 <- dat[,j]
      antibiotic_name1 <- dat_colnames[i]
      antibiotic_name2 <- dat_colnames[j]
      make_stationary_and_calc_ccf(series1, series2, antibiotic_name1, antibiotic_name2)
    }
  }
  dev.off()
}