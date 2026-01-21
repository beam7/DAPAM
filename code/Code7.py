from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re

#Adding Uniprot references for more usability

# === Load DAPAM AMP FASTA file ===
abp_records = []
for record in SeqIO.parse("/Users/beatricemihalache/Desktop/RADS/dapam_abp_sequences.fasta", "fasta"):
    # Example header: >ABP0001 | mechanism=membrane_disruption | subtype= | gram=- | pmc_id=PMC6267271
    match = re.match(r"^(\S+)\s+\|\s+mechanism=(.*?)\s+\|\s+subtype=(.*?)\s+\|\s+gram=(.*?)\s+\|\s+pmc_id=(.*)", record.description)
    if match:
        abp_id = match.group(1)
        mechanism = match.group(2)
        subtype = match.group(3)
        gram = match.group(4)
        pmc_id = match.group(5)
        abp_records.append({
            "abp_id": abp_id,
            "sequence": str(record.seq),
            "mechanism": mechanism,
            "subtype": subtype,
            "gram": gram,
            "pmc_id": pmc_id
        })

abp_df = pd.DataFrame(abp_records)

# === Load Swiss-Prot UniProt proteins ===
uniprot_records = list(SeqIO.parse("/Users/beatricemihalache/Desktop/RADS/uniprot_sprot.fasta", "fasta"))

# === Match each peptide to UniProt ===
matches = []
for i, abp in abp_df.iterrows():
    seq = abp['sequence']
    found_match = False

    for uniprot in uniprot_records:
        full_seq = str(uniprot.seq)
        header_parts = uniprot.id.split("|")

        if len(header_parts) >= 3:
            accession = header_parts[1]
            entry_name = header_parts[2]
        else:
            accession = None
            entry_name = None

        # Exact match
        if seq == full_seq:
            matches.append((
                abp['abp_id'], seq, accession, entry_name, "exact_match", 0
            ))
            found_match = True
            break

        # Subsequence match (not equal in length)
        elif seq in full_seq:
            position = full_seq.index(seq)
            matches.append((
                abp['abp_id'], seq, accession, entry_name, "subsequence", position
            ))
            found_match = True
            break

    if not found_match:
        matches.append((
            abp['abp_id'], seq, None, None, "no_match", None
        ))

# === Output final DataFrame with all info ===
match_df = pd.DataFrame(matches, columns=[
    "abp_id", "sequence", "uniprot_accession", "uniprot_entry",
    "match_type", "match_position"
])

# Merge with original metadata
final_df = pd.merge(abp_df, match_df, on=["abp_id", "sequence"], how="left")

# Save
final_df.to_csv("/Users/beatricemihalache/Desktop/RADS/dapam_with_uniprot_matches.csv", index=False)

fasta_records = []

for i, row in final_df.iterrows():
    header_parts = [
        f"mechanism={row['mechanism']}",
        f"subtype={row['subtype']}",
        f"gram={row['gram']}",
        f"pmc_id={row['pmc_id']}",
        f"uniprot={row['uniprot_accession']}" if pd.notnull(row["uniprot_accession"]) else "uniprot=None",
        f"entry={row['uniprot_entry']}" if pd.notnull(row["uniprot_entry"]) else "entry=None",
        f"match_type={row['match_type']}"
    ]
    full_description = " | ".join(header_parts)
    record = SeqRecord(
        Seq(row["sequence"]),
        id=row["abp_id"],
        description=full_description
    )
    fasta_records.append(record)

SeqIO.write(fasta_records, "/Users/beatricemihalache/Desktop/RADS/dapam_with_uniprot.fasta", "fasta")
final_df.to_json("/Users/beatricemihalache/Desktop/RADS/dapam_with_uniprot.json", orient="records", indent=2)
