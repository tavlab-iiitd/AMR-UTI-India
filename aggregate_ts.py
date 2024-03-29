"""
This code is used for generating time-series data for monthly antibiotic resistance percentages.
Here national time series are being generated. Same code is used for generating state-wise, zone-wise, age-wise and other time series.
"""

# Importing necessary libraries
import pandas as pd
import numpy as np
from urine_lib import parse_date

# Function to read data for a specific organism and filter it along with corresponding antibiotics
def read_df_for_sample_type(organism):
    df = pd.read_csv("../data/final/uti_final.csv", keep_default_na=True)

    df["gender"] = df["gender"].str.title()
    df["state"] = df["state"].str.title()
    df["zone"] = df["zone"].str.title()
    # Additional datasets for state and zone
    df_state = pd.read_csv("../data/final/same_state.csv")
    df_zone = pd.read_csv("../data/final/same_zone.csv")
    # Filtering data based on the specific organism and assigning relevant antibiotics
    if (organism == "Enterococcus faecium") or (organism == "Enterococcus faecalis"):
        df = df[df["organism_complete"] == organism]
        df_state = df_state[df_state["organism_complete"] == organism]
        df_zone = df_zone[df_zone["organism_complete"] == organism]

        antibiotics = ["Vancomycin" , "Nitrofurantoin" , "Linezolid" , "Teicoplanin" , "Ampicillin" , "Gentamicin HL" , "Ciprofloxacin" , "Fosfomycin"]

    elif organism == "Escherichia coli":
        df = df[df["organism_name"] == organism]
        df_state = df_state[df_state["organism_name"] == organism]
        df_zone = df_zone[df_zone["organism_name"] == organism]
        antibiotics = ["Amikacin" , "Piperacillin-tazobactam" , "Nitrofurantoin" , "Imipenem" , "Meropenem" , "Ciprofloxacin" , "Trimethoprim-sulfamethoxazole" , "Cefotaxime" , "Ertapenem" , "Levofloxacin" , "Fosfomycin" , "Cefazolin" , "Colistin"]

    elif organism == "Klebsiella":
        df = df[df["organism_name"] == organism]
        df_state = df_state[df_state["organism_name"] == organism]
        df_zone = df_zone[df_zone["organism_name"] == organism]
        antibiotics = ["Amikacin" , "Piperacillin-tazobactam" , "Nitrofurantoin" , "Imipenem" , "Meropenem" , "Ciprofloxacin" , "Trimethoprim-sulfamethoxazole" , "Ertapenem" , "Cefotaxime" , "Levofloxacin" , "Fosfomycin" , "Cefazolin" , "Colistin"]

    elif organism == "Pseudomonas aeruginosa":
        df = df[df["organism_name"] == organism]
        df_state = df_state[df_state["organism_name"] == organism]
        df_zone = df_zone[df_zone["organism_name"] == organism]
        antibiotics = ["Amikacin" , "Meropenem" , "Piperacillin-tazobactam" , "Ciprofloxacin" , "Imipenem" , "Cefepime" , "Ceftazidime" , "Gentamicin" , "Levofloxacin" , "Tobramycin" , "Colistin"]

    else:
        assert 0

    assert antibiotics is not None

    return df, df_state, df_zone, antibiotics

# Function to aggregate data, calculating the resistance percentage
def custom_aggregate(series, susceptibility="R"):
    unique_counts = round(series.value_counts(normalize=True) * 100, 2)
    return unique_counts[susceptibility] if susceptibility in unique_counts else np.nan

organisms =["Escherichia coli"]#, "Pseudomonas aeruginosa", "Klebsiella", "Enterococcus"]

for org in organisms:
    print(".......................Processing........................." + org)
    df, df_state, df_zone, antibiotics = read_df_for_sample_type(org)
    # Parsing and formatting the collection date
    df["collection_date"] = df.apply(lambda x: parse_date(x["collection_date"]).strftime('%m-%Y'), axis=1)
    month_year_val = {}
    for antibiotic in antibiotics:
        result = {}
        for month_year in df["collection_date"].unique():
            df_month_year = df[df["collection_date"] == month_year][antibiotic]
            if(len(df_month_year) == 0):
                result[month_year] = np.nan
                print("For ", antibiotic, "+", month_year, "no data is available")
            else:
                result[month_year] = df_month_year.agg(custom_aggregate)
        month_year_val[antibiotic] = result
        # Converting the aggregated data into a DataFrame and saving to CSV
        new_year = pd.DataFrame.from_dict(month_year_val)
        new_year.index.name = 'date'
    print(new_year)

    new_year.to_csv("../data/month_year_" + org + ".csv")