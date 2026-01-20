import pandas as pd
import json

# Generate clean CSV, FASTA, and JSON files of the dataset

# === File paths ===
mechanism_path = "processed_annotatedcopy_peptide_data_updated.csv"
peptipedia_path = "Peptipedia_Antibacterial_Features.csv"

# === Load CSV files ===
mechanism_df = pd.read_csv(mechanism_path)
peptipedia_df = pd.read_csv(peptipedia_path, low_memory=False)

# === Normalize column names ===
mechanism_df.columns = [col.strip().lower().replace(" ", "_") for col in mechanism_df.columns]
peptipedia_df.columns = [col.strip().lower().replace(" ", "_") for col in peptipedia_df.columns]

# === Filter out entries with no mechanism ===
filtered_df = mechanism_df[mechanism_df["mechanism"].notna() & (mechanism_df["mechanism"].str.strip() != "")]

# === Deduplicate sequences: keep all PMC IDs and mechanisms ===
grouped = filtered_df.groupby("sequence").agg({
    "mechanism": lambda x: ";".join(sorted(set(
        m.strip() for i in x.dropna() for m in str(i).split(",")
    ))),
    "pmc_id": lambda x: ";".join(sorted(set(str(i).strip() for i in x.dropna()))),
    "gram": lambda x: ";".join(sorted(set(x.dropna().astype(str))))
}).reset_index()

# === Mechanism mapping for substring matching ===
mechanism_mapping = [
    ("toroidal pore", ("pore", "toroidal_pore")),
    ("barrel-stave", ("pore", "barrel_stave")),
    ("barrel_stave", ("pore", "barrel_stave")),
    ("ion channel formation", ("pore", "ion_channel")),
    ("carpet", ("carpet", None)),
    ("membrane micelle formation", ("carpet", "membrane_micelle_formation")),
    ("pore", ("pore", None)),
    ("non-lytic", ("non-lytic", None)),
    ("membrane permeability", ("membrane_disruption", None)),
    ("membrane disruption", ("membrane_disruption", None)),
    ("membrane", ("membrane_disruption", None)),
    ("fusion of vesicles", ("membrane_disruption", "fusion of vesicles")),
    ("electrostatic interaction", ("membrane_disruption", "electrostatic_interaction")),
    ("electrostatic", ("membrane_disruption", "electrostatic_interaction")),
    ("biofilm destruction", ("biofilm_destruction", None)),
    ("antibiofilm", ("biofilm_destruction", None)),
    ("coaggregation of ribosomal protein", ("intracellular_targeting", "coaggregation_of_ribosomal_protein")),
    ("bacterial membrane external protrusion", ("membrane_disruption", "bacterial_membrane_external_protrusion")),
    ("immunomodulatory", ("immunomodulatory", None)),
    ("intracellular targeting", ("intracellular_targeting", None)),
    ("inhibition of outer membrane protein synthesis", ("non-lytic", "inhibition_outer_membrane_protein")),
]

# === Mechanism mapping function with safe sorting ===
def map_mechanism_to_labels(mechanism_string):
    mechanism_string = mechanism_string.lower()
    matches = set()
    for key, label in mechanism_mapping:  # <=== FIXED HERE
        if key in mechanism_string:
            matches.add(label)

    if not matches:
        return [("unspecified", "")]

    sorted_matches = sorted([(m[0], m[1] if m[1] is not None else "") for m in matches])
    return sorted_matches

# === Assign mechanism labels ===
grouped[["primary_mechanism", "subtype"]] = grouped["mechanism"].apply(
    lambda x: pd.Series(map_mechanism_to_labels(x)[0])
)

# === Create ID mapping from Peptipedia ===
sequence_to_id = {}
if "id" in peptipedia_df.columns:
    for _, row in peptipedia_df.iterrows():
        sequence = row["sequence"]
        seq_id = row["id"]
        sequence_to_id[sequence] = seq_id

# === Assign IDs to sequences ===
id_prefix = "ABP"
id_counter = 1
assigned_ids = []
for seq in grouped["sequence"]:
    if seq in sequence_to_id:
        assigned_ids.append(sequence_to_id[seq])
    else:
        while True:
            candidate_id = f"{id_prefix}{id_counter:04d}"
            id_counter += 1
            if candidate_id not in sequence_to_id.values():
                break
        sequence_to_id[seq] = candidate_id
        assigned_ids.append(candidate_id)
grouped.insert(0, "id", assigned_ids)

# === Save as CSV ===
grouped.to_csv("filtered_antibacterial_mechanisms.csv", index=False)

# === Write FASTA file ===
with open("filtered_antibacterial_mechanisms.fasta", "w") as fasta_file:
    for _, row in grouped.iterrows():
        header = f">{row['id']} | mechanism={row['primary_mechanism']} | subtype={row['subtype']} | gram={row['gram']} | pmc_id={row['pmc_id']}"
        sequence = row["sequence"]
        fasta_file.write(f"{header}\n{sequence}\n")

# === Write JSON file ===
json_output = {}
for _, row in grouped.iterrows():
    row_data = row.to_dict()
    seq_id = row_data.pop("id")
    json_output[seq_id] = row_data

with open("filtered_antibacterial_mechanisms.json", "w") as json_file:
    json.dump(json_output, json_file, indent=2)

print("CSV, FASTA, and JSON files created successfully.")