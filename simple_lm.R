set.seed(123)  # For reproducibility

library(nlme)
# Load data from a CSV file into a DataFrame 'dat'
dat <- read.csv("data.csv")
# Convert the antibiotic column of 'dat' into a time series object 'A' with monthly frequency
A <- ts(dat$Imipenem,frequency = 12)
# Decompose the time series 'A' into trend, seasonal, and irregular components using an additive model
dc <- decompose(A,type = "additive")
dc$figure
plot(dc)
# Extract the trend component from the decomposition, omitting NA values, and convert it to numeric
trend <- as.numeric(na.omit(dc$trend))
# Create a new DataFrame 'df' with two columns: 'Month' (time period) and 'Value' (trend values)
df <- data.frame("Month"=1:length(trend),
                    "Value"=trend)
# Plot the trend values against time and add the linear regression line
fit_int <- lme(value~ month, data = df, random = ~ 1|zone)
plot(df$Month,df$Value,type="l")
abline(lm)
# Display a summary of the linear regression model, including coefficients and statistical measures
summary(lm)

summ_freq<- summary(fit_int)

summ_freq$coefficients




