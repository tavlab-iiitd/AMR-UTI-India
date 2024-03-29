# Set the working directory and load necessary libraries for data manipulation, statistical modeling, and plotting
setwd ("~/Documents/AMR_model/amr_urine")

library(dplyr)
library(lubridate)
library(nlme) 
library(lattice)
library(lmtest)
library(sandwich)
library(stats)
set.seed(123)  # For reproducibility

# Load antibiotic resistance data
antibiotic_data <- read.csv("data/final/all_with_age.csv")


antibiotic = "Levofloxacin"

# Subset data and Convert the date to a month-year format
#Add month year in excel from collection date
subset_data <- antibiotic_data[, c("amr_id", "month_year", "organism_name", "hospital_name", antibiotic)]
org_data <- subset(subset_data, organism_name == 'Klebsiella')

# Calculate the resistance percentage for each hospital, each month
monthly_data <- org_data %>%
  group_by(hospital_name, month_year) %>%
  summarize(total_samples = n(),
            resistant_samples = sum( .data[[antibiotic]] == 'R', na.rm = TRUE),
            resistance_percentage = (resistant_samples / total_samples) * 100)

# Calculate national-level data
national_data <- org_data %>%
  group_by(month_year) %>%
  summarize(hospital_name = "National",  # Assign a label for national data
            total_samples = n(),
            resistant_samples = sum(.data[[antibiotic]] == 'R', na.rm = TRUE),
            resistance_percentage = (resistant_samples / total_samples) * 100)

# Append national data to the existing monthly_data
monthly_data_new <- rbind(monthly_data, national_data)

# Convert month_year string to a date time 
monthly_data_new$date <- parse_date_time(monthly_data_new$month_year, "%m-%Y")

data_decomp <- monthly_data_new[, c("resistance_percentage", "date", "hospital_name")]

# Decomposition
decomposed_data <- data_decomp %>%
  group_by(hospital_name) %>%
  do({
    # Ensure the data is ordered by time before decomposition
    .data <- arrange(., date)
    
    # Perform time series decomposition
    ts_data <- ts(.data$resistance_percentage, frequency = 12)
    decomposed <- tryCatch({
      decompose(ts_data, type = "additive")
    }, error = function(e) {
      return(NULL)
    })
    # If decomposition is successful, create a data frame with time and trend
    if (!is.null(decomposed)) {
      data.frame(date = .data$date, trend = as.numeric(decomposed$trend))
    } else {
      # If decomposition fails, return a data frame with NA values
      data.frame(date = .data$date, trend = NA)
    }
  }) %>% ungroup()

decomposed_data <- decomposed_data %>% filter(!is.na(trend))

# Define the start date
start_date <- as.Date("2017-07-01")

# Calculate the number of months since the start date for each date in the dataset
decomposed_data <- decomposed_data %>%
  mutate(seq_num = interval(start_date, date) / months(1))

# Ensure seq_num is an integer
decomposed_data$seq_num <- as.integer(decomposed_data$seq_num)
decomposed_data$hospital_name <- as.factor(decomposed_data$hospital_name)
decomposed_data$hospital_name <- relevel(decomposed_data$hospital_name, ref = "National")

filtered_data <- decomposed_data %>% filter(hospital_name == "National")

#Code peice to set anonyms for hospitals has been removed from public code file.

decomposed_data$hosp_anonym <- hospital_anonym[decomposed_data$hospital_name]

decomposed_data <- decomposed_data %>% 
  mutate(hospital_code = recode(hospital_name, !!!hospital_anonym))

#Mixed effects model
filtered_data_notnan <- decomposed_data %>% filter(hospital_name != "National")

model4 <- lme(trend ~ seq_num, random = ~ 1 + seq_num | hospital_code, data = filtered_data_notnan, method = "REML")

# Open the file for writing
file_name <- paste0("mixed_", antibiotic, "_ecoli.txt")
sink(file_name)

# Print the summary to the file
summary(model4)
random.effects(model4)

# Close the file
close(sink())
