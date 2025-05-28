import pandas as pd
import numpy as np
import re
import os


def clean_column_name(name):
    name = re.sub(r'[^0-9a-zA-Z_]', '', name)  # remove special characters
    name = re.sub(r'\s+', '_', name)  # convert spaces to underscores
    return name


# Load the EARM dataset
input_file = "EC-RP_IMS-RP_IC-RP_data_for_models.csv"
df = pd.read_csv(input_file)
df.columns = [clean_column_name(col) for col in df.columns]  # Clean column names

# Define the observables and matching average/variance columns
observables = {
    "mBid": ("ICRP", "nrm_var_ICRP"),
    "aSmac": ("IMSRP", "VAR"),
    "cPARP": ("ECRP", "nrm_var_ECRP")
}

# Build a long-format DataFrame
records = []
for obs_name, (avg_col, var_col) in observables.items():
    # Calculate normalized average
    avg = df[avg_col]
    norm_avg = (avg - avg.min()) / (avg.max() - avg.min())
    # Calculate standard error
    var = df[var_col]
    stderr = np.where(var_col == "VAR", 0.01, np.sqrt(var / 50))  # ~50 samples per data point (Albeck et al., 2008)
    stderr = np.where(pd.isnull(var), np.nan, stderr)

    obs_df = pd.DataFrame({
        "observable": obs_name,
        "time": df["Time"],
        "time_units": "seconds",
        "average": norm_avg,
        "stderr": stderr,
        "amount_units": "fraction",
        "expt_id": "A",
        "alt_expt_id": "Albeck2008"
    })
    records.append(obs_df)

# Concatenate all observables into one DataFrame
output_df = pd.concat(records, ignore_index=True)

# Export to CSV
output_file = os.path.splitext(os.path.basename(input_file))[0] + "_PyDREAM.csv"
output_df.to_csv(output_file, index=False)
