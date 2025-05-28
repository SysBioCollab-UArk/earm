import pandas as pd
import numpy as np

# Load your EARM dataset
input_path = "EC-RP_IMS-RP_IC-RP_data_for_models.csv"

df = pd.read_csv(input_path)

# Rename time column
df = df.rename(columns={"# Time": "time"})

# Define the observables and matching normalized/variance columns
observables = {
    "IC-RP": ("norm_IC-RP", "nrm_var_IC-RP"),
    "IMS-RP": ("IMS-RP step", "VAR"),  # this might need to be adjusted!
    "EC-RP": ("norm_EC-RP", "nrm_var_EC-RP")
}

# Create an empty list to hold the reformatted rows
reformatted_rows = []

# Go through each observable and reshape the data
for obs_name, (avg_col, var_col) in observables.items():
    for _, row in df.iterrows():
        average = row[avg_col]
        stderr = np.sqrt(row[var_col]) if pd.notnull(row[var_col]) else np.nan

        reformatted_rows.append({
            "observable": obs_name,
            "time": row["time"],
            "time_units": "sec",
            "average": average,
            "stderr": stderr,
            "amount_units": "a.u.",
            "expt_id": "EARM",
            "alt_expt_id": "ModelData2024"
        })

# Convert to DataFrame
output_df = pd.DataFrame(reformatted_rows)

# Export to CSV
output_df.to_csv("EARM_reformatted_output.csv", index=False)

print
