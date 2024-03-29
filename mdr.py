"""This script essentially processes data to determine multi-drug resistance in various organisms.
It does so by reading organism-specific data, analyzing antibiotic resistance, and then calculating the proportion of MDR cases.
Finally, it outputs the results to organism specific CSV files."""

#Importing necessary libraries
from urine_lib import read_df_for_sample_type

# List of organisms for analysis.
organisms = ["Escherichia coli ", "Enterococcus", "Pseudomonas aeruginosa", "Klebsiella" ]

# Mapping of organisms to their respective antibiotic groups and antibiotics
antibiotic_group_data = {
    "Escherichia coli ": {
        "cephems": ["Cefazolin", "Cefotaxim"],
        "Beta lactam inhibitor": ["Piperacillin-tazobactam"],
        "Carbapenems": ["Ertapenem", "Imipenem", "Meropenem"],
        "Lipopeptides": ["Colistin"],
        "Aminoglycosides":["Amikacin"],
        "Quinolones/Fluorquinolones": ["Ciprofloxacin", "Levofloxacin"],
        "Folate Pathway Antagonists": ["Trimethoprim-sulfamethoxazole"],
        "Fosfomycins": ["Fosfomycin"],
        "Nitrofurans": ["Nitrofurantoin"]
    },
    "Klebsiella": {
        "cephems": ["Cefazolin", "Cefotaxim", "Cefepime"],
        "Beta lactam inhibitor": ["Piperacillin-tazobactam"],
        "Carbapenems": ["Ertapenem", "Imipenem", "Meropenem"],
        "Lipopeptides": ["Colistin"],
        "Aminoglycosides":["Amikacin"],
        "Quinolones/Fluorquinolones": ["Ciprofloxacin", "Levofloxacin"],
        "Folate Pathway Antagonists": ["Trimethoprim-sulfamethoxazole"],
        "Fosfomycins": ["Fosfomycin"],
        "Nitrofurans": ["Nitrofurantoin"]
    },
    "Pseudomonas aeruginosa": {
        "cephems": ["Ceftazidime", "Cefepime"],
        "Beta lactam inhibitor": ["Piperacillin-tazobactam"],
        "Carbapenems": ["Imipenem", "Meropenem"],
        "Lipopeptides": ["Colistin"],
        "Aminoglycosides": ["Amikacin", "Gentamicin", "Tobramycin"],
        "Quinolones/Fluorquinolones": ["Ciprofloxacin", "Levofloxacin"]
    },
    "Enterococcus": {
        "Penicillins": ["Ampicillin"],
        "Glycopeptides": ["Vancomycin"],
        "Lipoglycopeptides": ["Teicoplanin"],
        "Fluorquinolones": ["Ciprofloxacin"],
        "Nitrofurans": ["Nitrofurantoin"],
        "Aminoglycosides": ["Gentamicin HL"],
        "Fosfomycins": ["Fosfomycin"],
        "Oxazolidinones": ["Linezolid"]
    }
}

for organism in organisms:
    assert organism in antibiotic_group_data

    status = "all"
    # Reading data for the given organism and status.
    df, antibiotics = read_df_for_sample_type(organism, status)

    # Dictionary to track unique amr_id values.
    group_dict = {}

    # Function to process each row in the dataframe.
    def process_groups(row):
        # Ensure amr_id is unique.
        assert row["amr_id"] not in group_dict
        ## Check amr id is unique
        group_dict[row["amr_id"]] = True
        # print("\n\n-----------------------------------------------------------------")
        mdr_count = 0
        # Counting the number of antibiotics for which the organism is resistant.
        for group in antibiotic_group_data[organism]:
            for antibiotic in antibiotics:
                if antibiotic in antibiotic_group_data[organism][group] and row[antibiotic] == "R":
                    mdr_count +=1
                    break

        # Returns True if the count of resistances is greater than 2, indicating MDR.
        return mdr_count > 2

    # Applying the process_groups function to each row and creating a new 'mdr' column.
    df["mdr"] = df.apply(lambda x: process_groups(x), axis=1)

    # Printing the proportion of MDR cases for each organism.
    print(organism, len(df[df["mdr"] == True])/len(df))

    # Saving the processed data to a CSV file.
    df.to_csv("../data/final/mdr/" + organism + "_mdr.csv")


