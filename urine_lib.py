"""This code comprises functions essential for the functionality of other codes"""

#Load required libraries
import pandas as pd
from datetime import datetime, timedelta


# This function loads a dataset, filters it based on the presence of certain organisms and the COVID status of the samples, and then returns the filtered DataFrame along with a list of relevant antibiotics for the specified organism.
def read_df_for_sample_type(organism, status):
    df = pd.read_csv("data.csv", keep_default_na=True)
    df["gender"] = df["gender"].str.title()
    df["state"] = df["state"].str.title()
    df["zone"] = df["zone"].str.title()

    if status == "covid":
        df = df[df["covid"] == 1]

    elif status == "noncovid":
        df = df[df["covid"] == 0]

    elif status == "all":
        pass

    if organism == "Enterococcus":
        df = df[df["organism_name"] == organism]
        antibiotics = ["Vancomycin" , "Nitrofurantoin" , "Linezolid" , "Teicoplanin" , "Ampicillin" , "Gentamicin HL" , "Ciprofloxacin" , "Fosfomycin"]

    elif organism == "Escherichia coli ":
        df = df[df["organism_name"] == organism]
        antibiotics = ["Amikacin" , "Piperacillin-tazobactam" , "Nitrofurantoin" , "Imipenem" , "Meropenem" , "Ciprofloxacin" , "Trimethoprim-sulfamethoxazole" , "Cefotaxime" , "Ertapenem" , "Levofloxacin" , "Fosfomycin" , "Cefazolin" , "Colistin"]

    elif organism == "Klebsiella":
        df = df[df["organism_name"] == organism]
        antibiotics = ["Amikacin" , "Piperacillin-tazobactam" , "Nitrofurantoin" , "Imipenem" , "Meropenem" , "Ciprofloxacin" , "Trimethoprim-sulfamethoxazole" , "Ertapenem" , "Cefotaxime" , "Levofloxacin" , "Fosfomycin" , "Cefazolin" , "Colistin"]

    elif organism == "Pseudomonas aeruginosa":
        df = df[df["organism_name"] == organism]
        antibiotics = ["Amikacin" , "Meropenem" , "Piperacillin-tazobactam" , "Ciprofloxacin" , "Imipenem" , "Cefepime" , "Ceftazidime" , "Gentamicin" , "Levofloxacin" , "Tobramycin" , "Colistin"]

    else:
        assert 0


    assert antibiotics is not None
    return df, antibiotics

def parse_date(date_str):
    #This function is used to parse dates from strings using multiple possible date formats

    dob_obj = None
    date_format_strs = ["%d/%m/%y", "%Y-%m-%d", "%b-%y", "%m-%y", "%m-%Y"]  ##, "%d/%m/-%y", "%Y/%m/%d", "%Y/%m/%d", "%d-%m-%Y", "%Y-%m-%d", "%d-%m-%Y %H:%M:%S"]
    for date_format_str in date_format_strs:
        try:
            dob_obj = datetime.strptime(date_str, date_format_str)
            break
        except Exception as e:
            pass

    if dob_obj != None and dob_obj.year > 2024:
        dob_obj = dob_obj.replace(year=(dob_obj.year-100))

    if(dob_obj == None):
        print("Bad date", date_str)

    assert(dob_obj != None)

    return dob_obj

def process_for_zone(row):
    #This function determines the geographical zone of India based on the state.

    if row["state"] in ["Himachal Pradesh", "Punjab", "Uttarakhand", "Uttar Pradesh", "Haryana", "Ladakh", "Chandigrah",
                        "Chandigarh", "Delhi", "DELHI", "Jammu And Kashmir", 'HIMACHAL PRADESH', 'HARYANA', 'PUNJAB']:
        return "North"
    if row["state"] in ["Andhra Pradesh", "Karnataka", "Kerala", "Telangana", "Tamil Nadu", "Andaman And Nicobar",
                        "Andaman and Nicobar", "Lakshadweep", 'LAKSHWADEEP', "Lakshwadeep", "Puducherry", 'Pondicherry',
                        'ANDHRA PRADESH', 'TAMIL NADU']:
        return "South"
    if row["state"] in ["Bihar", "Odisha", "Jharkhand", "West Bengal", 'WEST BENGAL', 'BIHAR' 'JHARKHAND']:
        return "East"
    if row["state"] in ["Rajasthan", 'RAJASTHAN', "Gujarat", 'Gujrat', "Goa", "Maharashtra", "Dadra And Nagar Haveli",
                        "Daman And Diu"]:
        return "West"
    if row["state"] in ["Madhya Pradesh", "madhya pradesh", "Chhattisgarh"]:
        return "Central"
    if row["state"] in ["Assam", "Sikkim", "Nagaland", "Meghalaya", "Manipur", "Mizoram", "Tripura",
                        "Arunachal Pradesh", 'ASSAM']:
        return "North-east"
    assert(0)


def process_for_hosp_state(row):
    #This function is used to return the hospital's state from IAMRSN network sites

    ## Site names excluded from public code
    assert(0)

def process_for_hosp_zone(row):
    # This function determines the hospital's geographical zone based on the hospital's state.

    if row["hosp_state"] in ["Himachal Pradesh", "Punjab", "Uttarakhand", "Uttar Pradesh", "Haryana", "Ladakh", "Chandigrah",
                        "Chandigarh", "Delhi", "Jammu and Kashmir"]:
        return "North"
    if row["hosp_state"] in ["Andhra Pradesh", "Karnataka", "Kerala", "Telangana", "Tamil Nadu", "Andaman And Nicobar",
                        "Andaman and Nicobar", "Lakshadweep", "Puducherry", 'Pondicherry']:
        return "South"
    if row["hosp_state"] in ["Bihar", "Odisha", "Jharkhand", "West Bengal"]:
        return "East"
    if row["hosp_state"] in ["Rajasthan", "Gujarat", 'Gujrat', "Goa", "Maharashtra", "Dadra and Nagar Haveli",
                        "Daman and Diu"]:
        return "West"
    if row["hosp_state"] in ["Madhya Pradesh", "madhya pradesh", "Chhattisgarh"]:
        return "Central"
    if row["hosp_state"] in ["Assam", "Sikkim", "Nagaland", "Meghalaya", "Manipur", "Mizoram", "Tripura",
                        "Arunachal Pradesh"]:
        return "North-east"

    assert(0)

def custom_aggregate(series, susceptibility="R"):
    # Calculates the percentage of a specified susceptibility for an antibiotic.
    unique_counts = round(series.value_counts(normalize=True) * 100, 2)
    return unique_counts[susceptibility] if susceptibility in unique_counts else np.nan

def custom_aggregate_number(series, susceptibility="R"):
    #This function returns the count of the specified susceptibility type and the total count in a formatted string.    
    unique_counts = series.value_counts(normalize=False)
    total_count = pd.Series.sum(unique_counts)
    return str(unique_counts[susceptibility])+ "/" + str(total_count) if susceptibility in unique_counts else np.nan