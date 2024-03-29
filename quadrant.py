import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

df = pd.read_csv('cluster_robust/quad/robust_se_Imipenem_ecoli.csv')

org = "Ecoli"
ab= "Imipenem"

# Separate the national average from the hospital data
national_average = df[df['hospital'] == 'National'].iloc[0]
hospitals = df[df['hospital'] != 'National']

# Plotting
plt.figure(figsize=(10, 8))

# Plot hospital data with labels
for index, row in hospitals.iterrows():
    plt.scatter(row['slope'], row['intercept'], color = "blue", s=12)
    # Adjust these values to change label offset
    x_offset = -0.04
    y_offset = 0.5
    plt.text(row['slope'] + x_offset, row['intercept'] + y_offset, row['hospital'], fontsize=10)


# Add lines for national average slope and intercept
plt.axhline(y=national_average['intercept'], color='red', linestyle='--', label='National Avg Intercept')
plt.axvline(x=national_average['slope'], color='green', linestyle='--', label='National Avg Slope')

# Labeling
plt.title('Slope and Intercept of Hospitals wrt National for '+ ab+ " in " + org)
plt.xlabel('Slope')
plt.ylabel('Intercept')
plt.legend()

# Save the figure
plt.savefig('cluster_robust/fig/' + org + "_" + ab + ".png", dpi=300, bbox_inches='tight')


#Correlation between slope and intercept
slope = hospitals['slope']
intercept = hospitals['intercept']
correlation = pearsonr(slope, intercept)

print("Correlation between slope and intercept is", correlation)
