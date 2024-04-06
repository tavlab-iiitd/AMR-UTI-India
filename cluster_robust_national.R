setwd ("~/Documents/AMR_model/amr_urine")

library(dplyr)
library(lubridate)
library(lmtest)
library(sandwich)
library(stats)
set.seed(123)  

antibiotic_data <- read.csv("data/final/all_with_age.csv")

# Initialize an empty dataframe to store results
all_results <- data.frame(Organism = character(),
                          Antibiotic = character(),
                          Estimate = numeric(),
                          StdError = numeric(),
                          TValue = numeric(),
                          PValue = numeric(),
                          DegreeOfFreedom = integer(),
                          RSquared = numeric(),
                          stringsAsFactors = FALSE)

# Function to read and preprocess dataset (simplified for integration)
read_and_filter_data <- function(status, organism) {
  df <- antibiotic_data %>% 
    mutate(gender = tools::toTitleCase(gender),
           state = tools::toTitleCase(state),
           zone = tools::toTitleCase(zone))
  
  if (status == "covid") {
    df <- df %>% filter(covid == 1)
  } else if (status == "noncovid") {
    df <- df %>% filter(covid == 0)
  }
  # Filtering by organism
  df <- df %>% filter(organism_name == organism)
  return(df)
}

# Define organisms and their respective antibiotics
organisms_antibiotics <- list(
  "Escherichia coli " = c("Amikacin" , "Piperacillin.tazobactam" , "Nitrofurantoin" , "Imipenem" , "Meropenem" , "Ciprofloxacin" , "Trimethoprim.sulfamethoxazole" , "Cefotaxime" , "Ertapenem" , "Levofloxacin" , "Fosfomycin" , "Cefazolin"),
  "Klebsiella" = c("Amikacin" , "Piperacillin.tazobactam" , "Nitrofurantoin" , "Imipenem" , "Meropenem" , "Ciprofloxacin" , "Trimethoprim.sulfamethoxazole" , "Ertapenem" , "Cefotaxime" , "Levofloxacin" , "Fosfomycin" , "Cefazolin"),
  "Pseudomonas aeruginosa" = c("Amikacin" , "Meropenem" , "Piperacillin.tazobactam" , "Ciprofloxacin" , "Imipenem" , "Cefepime" , "Ceftazidime" , "Gentamicin" , "Levofloxacin" , "Tobramycin"),
  "Enterococcus" = c("Vancomycin" , "Nitrofurantoin" , "Linezolid" , "Teicoplanin" , "Ampicillin" , "Gentamicin.HL" , "Ciprofloxacin" , "Fosfomycin")

)


# Loop over each organism and antibiotic
for (organism in names(organisms_antibiotics)) {
  for (antibiotic in organisms_antibiotics[[organism]]) {
    print(paste("Processing.....................", antibiotic))
    # Load and filter data for the current organism
    df <- read_and_filter_data("all", organism) 
    # Subset data and Convert the date to a month-year format
    #subset_data <- df[, c("amr_id", "month_year", "hospital_name", antibiotic)]
    # Calculate the resistance percentage for each hospital, each month
    monthly_data <- df %>%
      group_by(hospital_name, month_year) %>%
      summarize(total_samples = n(),
                resistant_samples = sum( .data[[antibiotic]] == 'R', na.rm = TRUE),
                resistance_percentage = (resistant_samples / total_samples) * 100)
  
    # Calculate national-level data
    national_data <- df %>%
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
    
    hospital_anonym <- c("Post Graduate Institute of Medical Education and Research" = "A",
                         "All India Institute Of Medical Sciences" = "B",
                         "Jawaharlal Institute of Postgraduate Medical Education and Research" = "C",
                         "Christian Medical College" = "D",
                         "P. D. Hinduja National Hospital and Medical Research Centre" = "E",
                         "Sir Ganga Ram Hospital" = "F",
                         "Armed Forces Medical College, Pune" = "G",
                         "Apollo Hospital" = "H",
                         "Mahatma Gandhi Institute of Medical Sciences" = "I",
                         "Tata Medical Center" = "J",
                         "Assam Medical College and Hospital" = "K",
                         "King George's Medical University" = "L",
                         "Nizam's Institute of Medical Sciences" = "M",
                         "All India Institute of Medical Sciences, Jodhpur" = "N",
                         "Kasturba Medical College,Manipal" = "O",
                         "Sher-i-kashmir Institute of Medical Sciences, Srinagar" = "P",
                         "All India Institute of Medical Sciences, Bhopal" = "Q",
                         "Regional Institute of Medical Sciences Hospital, Imphal" = "R",
                         "Institute of Post Graduate Medical Education &amp; Research" = "S",
                         "Lokmanya Tilak Municipal Medical College & General Hospital" = "T",
                         "JPN Apex Trauma Center, AIIMS" = "U",
                         "National" = "National"
    )
    
    decomposed_data <- decomposed_data %>% 
      mutate(hospital_code = recode(hospital_name, !!!hospital_anonym))
    
    model <- lm(trend ~ seq_num * factor(hospital_code), data = decomposed_data)
    summary_model <- summary(model)
    print(summary_model)
    robust_se <- coeftest(model, vcov = vcovHAC(model, type = "HC", cluster = "hospital_code", group = decomposed_data$hospital_code))
    print(robust_se)
    seq_num_info <- robust_se["seq_num", ]
    model_estimates <- data.frame(
      Organism = organism, 
      Antibiotic = antibiotic, 
      Estimate = seq_num_info["Estimate"],
      StdError = seq_num_info["Std. Error"],
      TValue = seq_num_info["t value"],
      PValue = seq_num_info["Pr(>|t|)"],  
      DegreeOfFreedom = summary_model$df[2],  
      RSquared = summary_model$r.squared
    )
    # Append to the results dataframe
    all_results <- rbind(all_results, model_estimates)
    write.csv(all_results, "national_robust_se.csv")
    }
}
