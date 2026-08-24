# This script filters a FASTA file of West Nile virus (WNV) sequences down
# to only those sequences present in Nextclade's CSV output AND that have
# Envelope (env) gene coverage - i.e. exactly the same sequence set used in
# the WNV_Envelope_Mutation_Summary workbook. Sequences with no envelope
# gene coverage (or entirely absent from the Nextclade output) are dropped.

import csv

# ---- EDIT: input/output paths ----
NEXTCLADE_CSV = "C:/Users/Dell/Downloads/nextclade (4).csv"
INPUT_FASTA   = "C:/Users/Dell/Downloads/Asia WNV project/Merged_sequences_Lineage2.fasta"
OUTPUT_FASTA  = "C:/Users/Dell/Downloads/Asia WNV project/Merged_sequences_Lineage2_env_only.fasta"
NEXTCLADE_DELIMITER = ";"  # Nextclade's default CSV export delimiter
TARGET_GENE = "env"        # Envelope protein gene name as used by Nextclade
# -----------------------------------


def base_accession(acc: str) -> str:
    """Strip version suffix (e.g. 'PX619833.1' -> 'PX619833')."""
    return acc.strip().split(".")[0]


def parse_cds_coverage(cds_coverage_str: str) -> dict:
    """Parse Nextclade's 'cdsCoverage' field, e.g. 'capsid:1,env:0.98,NS1:1'
    into {'capsid': 1.0, 'env': 0.98, 'NS1': 1.0}."""
    coverage = {}
    for part in cds_coverage_str.split(","):
        part = part.strip()
        if not part or ":" not in part:
            continue
        gene, value = part.split(":", 1)
        gene = gene.strip()
        value = value.strip()
        try:
            coverage[gene] = float(value) if value else 0.0
        except ValueError:
            coverage[gene] = 0.0
    return coverage


def parse_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file."""
    header = None
    seq_lines = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)


# ---- Step 1: Determine which accessions have envelope gene coverage ----
env_accessions = set()
no_env_accessions = set()

with open(NEXTCLADE_CSV, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter=NEXTCLADE_DELIMITER)
    for row in reader:
        acc = base_accession(row.get("seqName", ""))
        coverage = parse_cds_coverage(row.get("cdsCoverage", ""))
        if coverage.get(TARGET_GENE, 0.0) > 0.0:
            env_accessions.add(acc)
        else:
            no_env_accessions.add(acc)

print(f"Accessions with envelope gene coverage: {len(env_accessions)}")
print(f"Accessions without envelope gene coverage: {len(no_env_accessions)}")

# ---- Step 2: Filter the FASTA file ----
written = 0
skipped_no_env = 0
skipped_not_in_nextclade = 0

with open(OUTPUT_FASTA, "w", encoding="utf-8") as out:
    for header, seq in parse_fasta(INPUT_FASTA):
        acc = base_accession(header.split()[0])

        if acc in env_accessions:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i + 70] + "\n")
            written += 1
        elif acc in no_env_accessions:
            skipped_no_env += 1
        else:
            skipped_not_in_nextclade += 1

print(f"Sequences written:                    {written}")
print(f"Skipped (no envelope gene coverage):  {skipped_no_env}")
print(f"Skipped (not found in Nextclade CSV): {skipped_not_in_nextclade}")
print(f"Filtered FASTA written to '{OUTPUT_FASTA}'")