# Importing necessary libraries
import pandas as pd
import os
# Setting paths for the raw data and destination of processed data
root_path = "../data/nm_ped_new/ts_raw/"
dest_path = "../data/nm_ped_new/ts_processed/"

# Identifying directories and files within the root path
root_dirs = {}
for dir_name in os.listdir(root_path):
    if os.path.isdir(os.path.join(root_path, dir_name)) and "mdr" not in dir_name:
        root_dirs[dir_name] = []

for root_dir in root_dirs:
    root_dir_full = os.path.join(root_path, root_dir)
    for file_name in os.listdir(root_dir_full):
        if os.path.isfile(os.path.join(root_dir_full, file_name)) and ".csv" in file_name:
            root_dirs[root_dir].append(file_name)

print(root_dirs)

# Function to generate output file names based on the input parameters
def get_outfile_name(is_combined: bool, file_p: str, antibiotic: str = None):
    dir_p = os.path.dirname(file_p)
    dir_p = dir_p.replace(root_path, dest_path)
    dir_p = os.path.dirname(dir_p) + "/ts_decomp/" +   os.path.basename(dir_p)
    file_p = os.path.basename(file_p).split(".csv")[0]
    file_p = file_p.split("month_year_")[1]
    if is_combined:
        return dir_p + "/" + file_p + "_combined.csv"
    else:
        return dir_p + "/" + file_p + "_" + antibiotic + ".csv"

# Function to check the date range in the DataFrame
def check_date_range(df):
    # Check if DataFrame is empty
    if df.empty:
        return False

    # Get the first and last dates from the 'Date' column
    first_date = pd.to_datetime(df['date'].iloc[0])
    last_date = pd.to_datetime(df['date'].iloc[-1])

    return ((last_date.to_period('M') - first_date.to_period('M')).n + 1) - len(df)

# Function to insert missing rows into the DataFrame
def insert_missing_rows(df):
    # Check if DataFrame is empty
    if df.empty:
        return df

    # Convert 'date' column to datetime format
    df['date'] = pd.to_datetime(df['date'], format='(%B,%Y)')

    # Get the first and last dates from the 'date' column
    first_date = df['date'].iloc[0]
    last_date = df['date'].iloc[-1]

    # Calculate the number of months between the first and last dates
    num_months = (last_date.year - first_date.year) * 12 + (last_date.month - first_date.month) + 1

    # Generate a new DataFrame with the desired date range
    date_range = pd.date_range(start=first_date, end=last_date, freq='MS')
    new_df = pd.DataFrame({'date': date_range})

    # Merge the new DataFrame with the original DataFrame
    merged_df = pd.merge(new_df, df, on='date', how='left')

    return merged_df

# Main function to process each file
def do_for_file(file_p: str):
    print(" ----- Reading file", file_p)
    dat = pd.read_csv(file_p)
    try:
        dat['date'] = pd.to_datetime(dat['date'])
    except Exception as e:
        print("Error",dat['date'], e)
        exit(-1)

    dat = dat.sort_values('date')

    gap_dat = check_date_range(dat)
    if gap_dat != 0:
        print("******************** There are gaps in the time series", file_p, gap_dat)
        dat = insert_missing_rows(dat)
        # print(new_dat)
        # new_file_p = file_p.replace("month_year", "new_month_year")
        # new_dat.to_csv(new_file_p)
        # return


    Abs = dat.columns[1:]  # First column is date
    # Series with max 4 missing values in a stretch in between will be processed else discarded
    def fix_na(s: pd.Series, start_id: int, end_id: int):
        max_allowed_gap = 4
        if (end_id - start_id + 1) <= max_allowed_gap:
            avg_val = (s.get(start_id - 1) + s.get(end_id + 1)) / 2
            for i in range(start_id, end_id + 1):
                s.at[i] = avg_val
            return True
        return False

   #Only data with > 10 entries and less than 20% missing values will be processed
    def is_good_for_na(s: pd.Series):
        threshold = 0.2
        min_number_of_enteries = 10
        return (len(s) > min_number_of_enteries) and (s.isna().sum() <= threshold * len(s))

    def has_any_missing(s: pd.Series):
        return s.isna().sum() > 0

    na_df = dat.ffill().isna() | dat.bfill().isna()

    missing_abs = []
    for ab in Abs:
        df_ab = dat[["date", ab]]
        # print(df_ab)
        if has_any_missing(df_ab[ab]):
            missing_abs.append(ab)

    for ab in missing_abs:
        print("------ Processing missing ab", ab)
        df_ab = dat[ab]
        if not is_good_for_na(df_ab):
            print("-------------------------- Skipping missing ab", ab)
            continue


        df_ab = dat[["date", ab]][na_df[ab] == False]
        df_ab_no_index = df_ab.reset_index()
        df_ab_dates = df_ab_no_index["date"]
        df_ab = df_ab_no_index[ab]

        nas = df_ab.isna()
        found_na = None
        is_still_good_enough = True
        for i in range(len(df_ab)):
            if nas.get(i, ab) == True and found_na == None:
                ## We have found start of nan strech
                found_na = i
            elif nas.get(i) == False and found_na != None:
                ## We have found end of nan strech
                is_still_good_enough = is_still_good_enough and fix_na(df_ab, found_na, i - 1)
                found_na = None
                if is_still_good_enough == False:
                    break

        if is_still_good_enough == False:
            print("-------------------------- Skipping ab having very large gaps", ab)
            continue
        df_concat = pd.concat([df_ab_dates, df_ab], axis=1)
        df_concat['date'] = df_concat['date'].dt.strftime('%m-%Y')
        out_file_name = get_outfile_name(False, file_p, ab)
        print(out_file_name)
        df_concat.to_csv(out_file_name, index=False)

    if(len(missing_abs) != len(Abs)):
        data = dat.drop(missing_abs, axis=1)
        data['date'] = data['date'].dt.strftime('%m-%Y')
        out_file_name = get_outfile_name(True, file_p)
        print("-----------------------------------")
        print(out_file_name)
        data.to_csv(out_file_name, index=False)
    print("\n\n")

# Loop through each directory and file, and process them
for root_dir in root_dirs:
    root_dir_full = os.path.join(root_path, root_dir)
    for file_name in root_dirs[root_dir]:
        full_file_path = os.path.join(root_dir_full, file_name)
        do_for_file(full_file_path)

