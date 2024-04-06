setwd ("~/Documents/AMR_model/amr_urine")

set.seed(100)
my_plot.decomposed.ts = function(x, title="", ...) {
  xx <- x$x
  if (is.null(xx)) 
    xx <- with(x, if (type == "additive") 
      random + trend + seasonal
      else random * trend * seasonal)
  plot(cbind(observed = xx, trend = x$trend, seasonal = x$seasonal, random = x$random), 
       main=title, ...)
}


do_for_ab <- function(out_file_name, ab, df, org_name, file_category) {
  A <- ts(df[ab],frequency = 12)
  dc <- decompose(A,type = "additive")
  trend <- as.numeric(na.omit(dc$trend))
  df <- data.frame("Month"=1:length(trend),
                   "Value"=trend)

  lm <- lm(Value~Month,data=df)
  
  # Construct a filename for the model
  model_filename <- paste0("code/lm/simple_lm_age/models/", file_category, "_", org_name, "_", ab, ".rds")
  
  # Save the model object
  saveRDS(lm, model_filename)
  
  #print(paste("Saved model for", org_name, "-", ab, "as", model_filename))
  
  return(summary(lm))
}

get_all_files_for_type <- function(dat_type) {
  return(list.files(path = paste0("data/for_r_2024/",dat_type), pattern="*.csv"))
}

get_trend_p_value <- function(s) {
  p <- pf(s$fstatistic[1],            
          s$fstatistic[2],
          s$fstatistic[3],
          lower.tail = FALSE)
  pValue <- getElement(p, "value")
  trend <- s$coefficients[2]
  rsquared <- s$r.squared
  stdError <- s$coefficients[2, 2] # Standard error of the trend coefficient
  tValue <- s$coefficients[2, 3] # t-value of the trend coefficient
  degreeOfFreedom <- s$df[2] # Residual degrees of freedom
  return(list(pValue=pValue, trend=trend, rsquared=rsquared, stdError=stdError, tValue=tValue,  degreeOfFreedom=degreeOfFreedom))
  }

handle_file_inner <- function(csv_file, output_map, output_cat, categories) {
  print(paste0("--- Handling   ", csv_file))
  full_csv_path <- paste0("data/for_r_2024/", output_cat, "/", csv_file)
  
  for(category in categories) {
    if(is.null(output_map[[paste0(output_cat, "_", category, "_ab")]]))
      output_map[[paste0(output_cat, "_", category, "_ab")]] <- c()
    
    if(is.null(output_map[[paste0(output_cat, "_", category, "_p_values")]]))
      output_map[[paste0(output_cat, "_", category, "_p_values")]] <- c()
    
    if(is.null( output_map[[paste0(output_cat, "_", category, "_trend")]]))
      output_map[[paste0(output_cat, "_", category, "_trend")]] <- c()
    
    if(is.null( output_map[[paste0(output_cat, "_", category, "_rsquared")]]))
      output_map[[paste0(output_cat, "_", category, "_rsquared")]] <- c()
    
    if(is.null( output_map[[paste0(output_cat, "_", category, "_stdError")]]))
      output_map[[paste0(output_cat, "_", category, "_stdError")]] <- c()
    if(is.null( output_map[[paste0(output_cat, "_", category, "_tValue")]]))
      output_map[[paste0(output_cat, "_", category, "_tValue")]] <- c()
    if(is.null( output_map[[paste0(output_cat, "_", category, "_degreeOfFreedom")]]))
      output_map[[paste0(output_cat, "_", category, "_degreeOfFreedom")]] <- c()
  }
  
  org_name <- NULL
  file_category <- NULL
  # Attempt to directly match categories in the filename
  for(category in categories) {
    if(grepl(paste0("_", category), csv_file)) {
      # If category is preceded by an underscore
      parts <- strsplit(csv_file, paste0("_", category))[1]
      org_name <- unlist(parts)[1]
      file_category <- category
      break
    } else if(grepl(category, csv_file)) {
      # Direct concatenation of organism name and category
      category_index = regexpr(category, csv_file)
      if(category_index[1] != -1) {
        org_name <- substr(csv_file, 1, category_index - 1)
        file_category <- category
        break
      }
    }
  }
  
  if(is.null(org_name)) {
    print("\n\n\n\n----------------Error----------------\n\n\n")
    return(output_map)
  }
  
  dat <- read.csv(full_csv_path)
  print(paste("Rows:", nrow(dat), "Columns:", ncol(dat)))
  
  abs <- colnames(dat)[2:ncol(dat)]
  print(abs)
  
  for(ab in abs) {
    out_file_name <- paste0(output_cat,"/", org_name, "_", ab, "_", file_category)
    
    ab_list_key <- paste0(output_cat, "_", file_category, "_ab") 
    output_map[[ab_list_key]] <- c(output_map[[ab_list_key]], paste0(org_name, "_", ab))
    
    s <- do_for_ab(out_file_name, ab, dat, org_name, file_category)
    s_summary_vals <- unlist(get_trend_p_value(s))
    
    p_value_list_key <- paste0(output_cat, "_", file_category, "_p_values")
    trend_list_key <- paste0(output_cat, "_", file_category, "_trend")
    rsquared_list_key <- paste0(output_cat, "_", file_category, "_rsquared")  
    stderror_list_key <- paste0(output_cat, "_", file_category, "_stdError")
    tvalue_list_key <- paste0(output_cat, "_", file_category, "_tValue")
    dof_list_key <- paste0(output_cat, "_", file_category, "_degreeOfFreedom")  
    
    output_map[[p_value_list_key]] <- c(output_map[[p_value_list_key]], s_summary_vals[1])
    output_map[[trend_list_key]] <- c(output_map[[trend_list_key]], s_summary_vals[2])
    output_map[[rsquared_list_key]] <- c(output_map[[rsquared_list_key]], s_summary_vals[3])  
    output_map[[stderror_list_key]] <- c(output_map[[stderror_list_key]], s_summary_vals[4])
    output_map[[tvalue_list_key]] <- c(output_map[[tvalue_list_key]], s_summary_vals[5])
    output_map[[dof_list_key]] <- c(output_map[[dof_list_key]], s_summary_vals[6])  
  }
  
  return(output_map)
}

read_dir_inner <- function(category, category_types) {
  all_files <- get_all_files_for_type(category)
  output_map <- list()
  for(csv_file in all_files) {
    #if(grepl("Infant" ,csv_file)) {
      output_map <- handle_file_inner(csv_file, output_map, category, category_types)
      #if (is.null(output_map[["age_Infants_ab"]])) {
        #print("Hello")
      #}
    #}
  }
  
  for(category_type in category_types) {
    Antibiotic <- output_map[[paste0(category, "_", category_type, "_ab")]]
    p_values <- output_map[[paste0(category, "_", category_type, "_p_values")]]
    trend <- output_map[[paste0(category, "_", category_type, "_trend")]]
    rsquared <- output_map[[paste0(category, "_", category_type, "_rsquared")]]
    stdError<- output_map[[paste0(category, "_", category_type, "_stdError")]]
    tValue <- output_map[[paste0(category, "_", category_type, "_tValue")]]
    degreeOfFreedom <- output_map[[paste0(category, "_", category_type, "_degreeOfFreedom")]]
    df_category <- data.frame(Antibiotic, p_values, trend, rsquared, stdError, tValue, degreeOfFreedom)
    write.csv(x=df_category, file=paste0("code/lm/simple_lm_age/", category_type, ".csv"))
  }
}

read_dir_inner("age_ped", c("Infants", "Pre-schooler", "Middle childhood", "Toddlers", "Teenagers"))
read_dir_inner("age", c("Neonate", "Pediatric", "Middle", "Reproductive", "Elderly"))