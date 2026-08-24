import csv
import os

OUT_DIR = "C:/Users/Dell/Downloads" # where the merged files will be saved

CSV_OUT   = os.path.join(OUT_DIR, "Merged_sequences_Lineage2.csv")
FASTA_OUT = os.path.join(OUT_DIR, "Merged_sequences_Lineage2.fasta")

REQUIRED_COL = "Accession"

# ---- EDIT: put the full path to each file here ----
CSV_FILES = [
    "C:/Users/Dell/Downloads/Asia WNV project/Asia_sequences_Lineage2.csv",
    "C:/Users/Dell/Downloads/Asia WNV project/Africa_sequences_Lineage2.csv",
    "C:/Users/Dell/Downloads/Asia WNV project/Europe _sequences_Lineage2.csv",
    "C:/Users/Dell/Downloads/Asia WNV project/(taxid2421440)_sequences_Lineage2.csv",
]

FASTA_FILES = [
    "C:/Users/Dell/Downloads/Asia WNV project/Africa_sequences_results_Lineage2.fasta",
    "C:/Users/Dell/Downloads/Asia WNV project/Asia_sequences_results_Lineage2.fasta",
    "C:/Users/Dell/Downloads/Asia WNV project/Europe_sequences_results_Lineage2.fasta",
    "C:/Users/Dell/Downloads/Asia WNV project/(taxid2421440)_results_sequences_Lineage2.fasta",
]
# -----------------------------------------------------

# ---- Step 1: Merge CSVs ----
all_rows = []
all_fieldnames = []
seen_accessions = set()

for path in CSV_FILES:
    if not os.path.exists(path):
        print(f"skipped (not found): {path}")
        continue
    with open(path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        if REQUIRED_COL not in fieldnames:
            print(f"skipped (no Accession column): {path}")
            continue
        for col in fieldnames:
            if col not in all_fieldnames:
                all_fieldnames.append(col)
        count = 0
        for row in reader:
            acc = row.get("Accession", "").strip()
            if not acc:
                continue
            acc_base = acc.split(".")[0]
            if acc_base in seen_accessions:
                continue
            seen_accessions.add(acc_base)
            all_rows.append(row)
            count += 1
        print(f"used: {os.path.basename(path)} -> {count} rows")

with open(CSV_OUT, "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=all_fieldnames)
    writer.writeheader()
    for row in all_rows:
        writer.writerow({col: row.get(col, "") for col in all_fieldnames})

print(f"Total merged CSV rows: {len(all_rows)} -> {CSV_OUT}")

# ---- Step 2: Merge FASTA files ----
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

seen_fasta_accessions = set()
written = 0

with open(FASTA_OUT, "w", encoding='utf-8') as out:
    for path in FASTA_FILES:
        if not os.path.exists(path):
            print(f"skipped (not found): {path}")
            continue
        count = 0
        for header, seq in parse_fasta(path):
            acc = header.split()[0].strip()
            acc_base = acc.split(".")[0]
            if acc_base in seen_fasta_accessions:
                continue
            seen_fasta_accessions.add(acc_base)
            out.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i+70] + "\n")
            written += 1
            count += 1
        print(f"used: {os.path.basename(path)} -> {count} sequences")

print(f"Total merged FASTA sequences: {written} -> {FASTA_OUT}")