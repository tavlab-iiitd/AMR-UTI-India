#Import required libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv("../data/final/all_with_age.csv")

df["state"] = df["state"].str.title()
print(df["state"].unique())

# Get the state counts
state_counts = df['state'].value_counts()

# Create a bar chart
plt.figure(figsize=(8, 8))
state_counts.plot(kind='bar', color='mediumslateblue', width = 0.8)
# Add value labels above each bar
for i, v in enumerate(state_counts.values):
  plt.text(i, v + 0.1, str(v), ha='center', va='bottom', fontsize=6)  # Adjust v + 0.1 for label position

plt.title('State wise number of isolates', fontsize = 14)
plt.xlabel('State', fontsize = 12)
plt.ylabel('Number of isolates', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.tight_layout()
plt.savefig('state.pdf')

df["gender"] = df["gender"].str.title()

organism_counts = df.groupby('gender')['organism_name'].value_counts().unstack(fill_value=0)

# Color dictionary
ab_color_dict = {
    "Pseudomonas aeruginosa": "dodgerblue",
    "Klebsiella": "lightcoral",
    "Escherichia coli ": "sandybrown",  # Ensure this key exactly matches your DataFrame's value
    "Enterococcus": "mediumslateblue"
}

# Check that all organism names have an entry in the color dictionary, otherwise default to a specified color
default_color = "gray"
organism_colors = [ab_color_dict.get(x, default_color) for x in organism_counts.columns]

# Plotting with the defined colors for each segment
plt.figure(figsize=(15, 10))
# Width of 20 inches and height of 10 inches
ax = organism_counts.plot(kind='bar', stacked=True, color=organism_colors)

plt.title('Distribution of Organisms Across Genders', fontsize=12)
plt.xlabel('Gender', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.xticks(fontsize=12)  # Rotate labels for better readability
plt.yticks(fontsize=12)

# Adjust legend to top right
plt.legend(loc='upper right', fontsize=8)

# Add data labels within the bars, adjusted for stacked bars
for rect in ax.patches:  # iterate through all bars
    # Calculate width and height of each bar
    width, height = rect.get_width(), rect.get_height()
    x, y = rect.get_xy()
    label_text = f'{int(height)}'  # The label is the height of the bar

    # Only add labels to non-zero bars
    if height > 0:
        ax.text(x + width / 2, y + height / 2, label_text,
                ha='center', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('gender_org.pdf')

organism_counts = df.groupby('infection_type')['organism_name'].value_counts().unstack(fill_value=0)

# Color dictionary
ab_color_dict = {
    "Pseudomonas aeruginosa": "dodgerblue",
    "Klebsiella": "lightcoral",
    "Escherichia coli ": "sandybrown",  # Ensure this key exactly matches your DataFrame's value
    "Enterococcus": "mediumslateblue"
}

# Check that all organism names have an entry in the color dictionary, otherwise default to a specified color
default_color = "gray"
organism_colors = [ab_color_dict.get(x, default_color) for x in organism_counts.columns]

# Plotting with the defined colors for each segment
plt.figure(figsize=(15, 10))
# Width of 20 inches and height of 10 inches
ax = organism_counts.plot(kind='bar', stacked=True, color=organism_colors)

plt.title('Distribution of Organisms Across Infection type', fontsize=12)
plt.xlabel('Infection type', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.xticks(fontsize=12)  # Rotate labels for better readability
plt.yticks(fontsize=12)

# Adjust legend to top right
plt.legend(loc='upper left', fontsize=8)

# Add data labels within the bars, adjusted for stacked bars
for rect in ax.patches:  # iterate through all bars
    # Calculate width and height of each bar
    width, height = rect.get_width(), rect.get_height()
    x, y = rect.get_xy()
    label_text = f'{int(height)}'  # The label is the height of the bar

    # Only add labels to non-zero bars
    if height > 0:
        ax.text(x + width / 2, y + height / 2, label_text,
                ha='center', va='center', fontsize=8)
plt.tight_layout()

plt.savefig('inf_org.pdf')

organism_counts = df.groupby('year')['organism_name'].value_counts().unstack(fill_value=0)

# Color dictionary
ab_color_dict = {
    "Pseudomonas aeruginosa": "dodgerblue",
    "Klebsiella": "lightcoral",
    "Escherichia coli ": "sandybrown",  # Ensure this key exactly matches your DataFrame's value
    "Enterococcus": "mediumslateblue"
}

# Check that all organism names have an entry in the color dictionary, otherwise default to a specified color
default_color = "gray"
organism_colors = [ab_color_dict.get(x, default_color) for x in organism_counts.columns]

# Plotting with the defined colors for each segment
plt.figure(figsize=(15, 10))
# Width of 20 inches and height of 10 inches
ax = organism_counts.plot(kind='bar', stacked=True, color=organism_colors)

plt.title('Year-wise distribution of Organisms', fontsize=12)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.xticks(fontsize=12)  # Rotate labels for better readability
plt.yticks(fontsize=12)
plt.tight_layout()

# Adjust legend to top right
plt.legend(loc='upper left', fontsize=8)

# Add data labels within the bars, adjusted for stacked bars
for rect in ax.patches:  # iterate through all bars
    # Calculate width and height of each bar
    width, height = rect.get_width(), rect.get_height()
    x, y = rect.get_xy()
    label_text = f'{int(height)}'  # The label is the height of the bar

    # Only add labels to non-zero bars
    if height > 0:
        ax.text(x + width / 2, y + height / 2, label_text,
                ha='center', va='center', fontsize=8)
plt.savefig('year_org.pdf')

