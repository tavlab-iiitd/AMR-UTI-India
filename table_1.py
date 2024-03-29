#Import required libraries
import pandas as pd
import numpy as np
# urine_lib is a custom library, containing functions for parsing dates,reading data, and custom aggregation functions
from urine_lib import parse_date, read_df_for_sample_type, custom_aggregate_number, custom_aggregate

# Define a list of organisms to be analyzed
organisms = ["Escherichia coli ", "Enterococcus", "Pseudomonas aeruginosa", "Klebsiella"]

# Loop through each organism in the list
for organism in organisms:
    status = "all"
    df, antibiotics = read_df_for_sample_type(organism, status)
    print("---------------------------------------------" + organism)
    for antibiotic in antibiotics:
        result = df[antibiotic].agg(custom_aggregate_number)
        result2 = df[antibiotic].agg(custom_aggregate)
        print(antibiotic + ", " + str(result) + "(" + str(result2)+ ") " )