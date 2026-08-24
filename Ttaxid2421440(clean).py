import csv

# ---- EDIT THESE PATHS ----
TYPING_CSV = "C:/Users/Dell/Downloads/(taxid2421440) Ater Typing tooresults (2).csv"
META_CSV   = "C:/Users/Dell/Downloads/(Ttaxid2421440)sequences.csv"
FASTA_IN   = "C:/Users/Dell/Downloads/(taxid2421440) Ater Typing tool_sequences.fasta"

TYPING_OUT = "C:/Users/Dell/Downloads/Asia WNV project/(taxid2421440) typing_results_Lineage2.csv"
META_OUT   = "C:/Users/Dell/Downloads/Asia WNV project/(taxid2421440)_sequences_Lineage2.csv"
FASTA_OUT  = "C:/Users/Dell/Downloads/Asia WNV project/(taxid2421440)_sequences_Lineage2.fasta"

TARGET_LINEAGE = "Lineage 2"
# ---------------------------

# Step 1: get accessions that are Lineage 2 from the typing CSV
lineage2_names_full = set()   # e.g. "PX625775.1" (with version)
lineage2_names_base = set()   # e.g. "PX625775"    (no version)

with open(TYPING_CSV, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    rows = [row for row in reader if row["type"].strip() == TARGET_LINEAGE]
    fieldnames = reader.fieldnames

with open(TYPING_OUT, "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

for row in rows:
    name = row["name"].strip()
    lineage2_names_full.add(name)
    lineage2_names_base.add(name.split(".")[0])

print(f"Lineage 2 records found in typing results: {len(rows)}")

# Step 2: filter metadata CSV by Accession
# Handles Accession column whether it includes the version (e.g. "PQ493387.1")
# or not (e.g. "PQ493387") by checking both forms.
with open(META_CSV, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    meta_fieldnames = reader.fieldnames
    meta_rows = []
    for row in reader:
        acc = row["Accession"].strip()
        acc_base = acc.split(".")[0]
        if acc in lineage2_names_full or acc_base in lineage2_names_base:
            meta_rows.append(row)

with open(META_OUT, "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=meta_fieldnames)
    writer.writeheader()
    writer.writerows(meta_rows)

print(f"Lineage 2 records found in metadata CSV: {len(meta_rows)}")

# Step 3: filter FASTA sequences by header (matches full accession with version)
def parse_fasta(path):
    header = None
    seq_lines = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)

written = 0
with open(FASTA_OUT, "w", encoding='utf-8') as out:
    for header, seq in parse_fasta(FASTA_IN):
        acc = header.split()[0].strip()
        if acc in lineage2_names_full:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i+70] + "\n")
            written += 1

print(f"Lineage 2 sequences written to FASTA: {written}")