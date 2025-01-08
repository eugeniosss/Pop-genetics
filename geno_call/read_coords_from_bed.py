import sys
import json

# Read the input file
input_file = snakemake.input[0]
output_log = snakemake.output[0]

# Initialize the dictionary
DIC = {}

with open(input_file, "r") as file:
    for idx, line in enumerate(file, start=1):
        values = line.strip().split()
        DIC[idx] = values

# Write the dictionary to the log file
with open(output_log, "w") as log_file:
    json.dump(DIC, log_file, indent=4)

# Store the dictionary in a global variable
globals()["global_dictionary"] = DIC
