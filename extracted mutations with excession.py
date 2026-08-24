# This script merges West Nile virus (WNV) Nextclade output with the NCBI
# sequence metadata CSV, and extracts per-sequence nucleotide mutations plus
# Envelope (env) protein amino-acid mutations into a single analysis-ready CSV.
#
# It is adapted from an earlier Genome Detective XML-parsing script used for
# Rift Valley fever virus (RVFV) M-segment mutation extraction. Since WNV
# genomes are translated as a single polyprotein (rather than RVFV's M
# segment with Gn/Gc glycoproteins), this version:
#   - reads Nextclade's CSV export instead of Genome Detective's XML export
#   - merges by sequence ID (accession, version suffix stripped from both
#     sides, since the two files don't consistently include it)
#   - extracts ONLY Envelope ("env") amino-acid substitutions from
#     Nextclade's aaSubstitutions column (which lists ALL genes as
#     "GENE:mutation" pairs) and counts them
#   - extracts nucleotide substitutions and counts them (plus total
#     nucleotide mutations = substitutions + deletions + insertions)
#   - in the frequency sheet, also lists WHICH accessions carry each
#     envelope mutation (overall, and split by Historical/Novel range)

import csv
import os
import re
from collections import Counter, defaultdict

# ---- EDIT: input/output paths ----
NEXTCLADE_CSV = "C:/Users/Dell/Downloads/nextclade (4).csv"
METADATA_CSV  = "C:/Users/Dell/Downloads/Asia WNV project/Merged_sequences_Lineage2_classified.csv"
OUTPUT_CSV    = "C:/Users/Dell/Downloads/WNV_Envelope_Mutation_Summary.csv"
FREQUENCY_CSV = "C:/Users/Dell/Downloads/WNV_Envelope_Mutation_Frequency.csv"
NEXTCLADE_DELIMITER = ";"  # Nextclade's default CSV export delimiter
TARGET_GENE = "env"        # Envelope protein gene name as used by Nextclade

# Exact category strings used in the metadata's Range_Classification column
HISTORICAL_LABEL = "Historical/Native Range (Endemic)"
NOVEL_LABEL = "Novel/Introduced Range (Recently Emerged)"

# Separator used to join the list of accessions in the frequency sheet's
# Accessions columns (comma keeps it human-readable while still being a
# single CSV cell; change to e.g. "; " if you prefer)
ACCESSION_JOIN_SEP = ", "
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

# ---- Step 1: Load metadata, keyed on base accession ----
metadata = {}
with open(METADATA_CSV, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        acc = base_accession(row.get("Accession", ""))
        if acc:
            metadata[acc] = row

print(f"Loaded metadata for {len(metadata)} accessions from {METADATA_CSV}")

# ---- Step 2: Parse Nextclade output, merge by sequence ID, extract mutations ----
records = []
# (Range_Classification, Accession, [env substitution strings]) per sequence
mutation_tracking = []
matched, unmatched = 0, 0
skipped_no_env = 0

with open(NEXTCLADE_CSV, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter=NEXTCLADE_DELIMITER)
    for row in reader:
        seq_name = row.get("seqName", "").strip()
        acc = base_accession(seq_name)

        # ---- Skip sequences that don't have the envelope gene ----
        # Nextclade's cdsCoverage field lists coverage per gene; if 'env'
        # is absent (or has 0 coverage), the gene wasn't aligned/present
        # in this sequence, so it's excluded entirely.
        cds_coverage = parse_cds_coverage(row.get("cdsCoverage", ""))
        if cds_coverage.get(TARGET_GENE, 0.0) <= 0.0:
            skipped_no_env += 1
            continue

        meta = metadata.get(acc, {})
        if meta:
            matched += 1
        else:
            unmatched += 1

        # ---- Nucleotide mutations ----
        nt_substitutions = row.get("substitutions", "").strip()
        total_nt_substitutions = int(row.get("totalSubstitutions", "0").strip() or 0)
        total_nt_deletions = int(row.get("totalDeletions", "0").strip() or 0)
        total_nt_insertions = int(row.get("totalInsertions", "0").strip() or 0)
        total_nucleotide_mutations = total_nt_substitutions + total_nt_deletions + total_nt_insertions

        # ---- Amino-acid substitutions: filter to Envelope (env) only ----
        # Nextclade's aaSubstitutions column lists every gene as
        # "GENE:mutation", comma-separated, e.g. "env:E55D,NS1:Y34F,..."
        aa_substitutions = row.get("aaSubstitutions", "").strip()
        aa_deletions = row.get("aaDeletions", "").strip()

        env_substitutions = []
        if aa_substitutions:
            for mut in aa_substitutions.split(","):
                mut = mut.strip()
                if not mut or ":" not in mut:
                    continue
                gene, change = mut.split(":", 1)
                if gene.strip() == TARGET_GENE:
                    env_substitutions.append(change.strip())

        env_deletions = []
        if aa_deletions:
            for mut in aa_deletions.split(","):
                mut = mut.strip()
                if not mut or ":" not in mut:
                    continue
                gene, change = mut.split(":", 1)
                if gene.strip() == TARGET_GENE:
                    env_deletions.append(change.strip())

        total_env_substitutions = len(env_substitutions)
        total_env_deletions = len(env_deletions)
        total_env_mutations = total_env_substitutions + total_env_deletions

        records.append({
            "Accession": acc,
            "Sequence_ID": seq_name,
            "Country": meta.get("Country", ""),
            "Geo_Location": meta.get("Geo_Location", ""),
            "Host": meta.get("Host", ""),
            "Collection_Date": meta.get("Collection_Date", ""),
            "Isolate": meta.get("Isolate", ""),
            "Length": meta.get("Length", ""),
            "Nuc_Completeness": meta.get("Nuc_Completeness", ""),
            "Range_Classification": meta.get("Range_Classification", ""),
            "Nextclade_Clade": row.get("clade", "").strip(),
            "QC_Overall_Status": row.get("qc.overallStatus", "").strip(),
            "Nucleotide_Substitutions": nt_substitutions,
            "Total_Nt_Substitutions": total_nt_substitutions,
            "Total_Nt_Deletions": total_nt_deletions,
            "Total_Nt_Insertions": total_nt_insertions,
            "Total_Nucleotide_Mutations": total_nucleotide_mutations,
            "Envelope_Substitutions": ", ".join(env_substitutions),
            "Total_Envelope_Substitutions": total_env_substitutions,
            "Envelope_Deletions": ", ".join(env_deletions),
            "Total_Envelope_Deletions": total_env_deletions,
            "Total_Envelope_Mutations": total_env_mutations,
        })
        mutation_tracking.append((meta.get("Range_Classification", ""), acc, env_substitutions))

print(f"Nextclade rows processed: {len(records)}")
print(f"  skipped (no envelope gene): {skipped_no_env}")
print(f"  matched to metadata:   {matched}")
print(f"  no metadata match:     {unmatched}")

# ---- Step 3: Write merged output ----
fieldnames = list(records[0].keys()) if records else []
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(records)

print(f"Envelope mutation summary written to '{OUTPUT_CSV}' ({len(records)} rows)")

# ---- Step 4: Build the envelope mutation-frequency sheet ----
# For every distinct envelope mutation seen across all sequences, count how
# many sequences carry it: overall, within the Historical/Native Range
# (Endemic) group, and within the Novel/Introduced Range (Recently Emerged)
# group. Frequency = count / number of sequences in that group. Mutations
# are also flagged as Historical-only, Novel-only, Shared, or (rarely)
# found only in sequences with an unclassified/unknown range.
#
# NEW: alongside each count, we also record WHICH accessions carry the
# mutation (overall, and split into Historical/Novel accession lists) so
# the frequency sheet is traceable back to individual sequences.

total_all = len(mutation_tracking)
total_historical = sum(1 for rc, _, _ in mutation_tracking if rc == HISTORICAL_LABEL)
total_novel = sum(1 for rc, _, _ in mutation_tracking if rc == NOVEL_LABEL)

total_count = Counter()
historical_count = Counter()
novel_count = Counter()

# mutation -> list of accessions (overall / historical / novel)
total_accessions = defaultdict(list)
historical_accessions = defaultdict(list)
novel_accessions = defaultdict(list)

for range_classification, acc, muts in mutation_tracking:
    for mut in set(muts):  # count/list each mutation once per sequence
        total_count[mut] += 1
        total_accessions[mut].append(acc)
        if range_classification == HISTORICAL_LABEL:
            historical_count[mut] += 1
            historical_accessions[mut].append(acc)
        elif range_classification == NOVEL_LABEL:
            novel_count[mut] += 1
            novel_accessions[mut].append(acc)

def mutation_sort_key(mut: str):
    match = re.search(r"(\d+)", mut)
    return int(match.group(1)) if match else 0

frequency_rows = []
for mut in sorted(total_count.keys(), key=mutation_sort_key):
    t_count = total_count[mut]
    h_count = historical_count.get(mut, 0)
    n_count = novel_count.get(mut, 0)

    if h_count > 0 and n_count > 0:
        presence = "Shared (Historical & Novel)"
    elif h_count > 0:
        presence = "Historical/Native Range only"
    elif n_count > 0:
        presence = "Novel/Introduced Range only"
    else:
        presence = "Unclassified range only"

    frequency_rows.append({
        "Mutation": mut,
        "Total_Count": t_count,
        "Total_Frequency": round(t_count / total_all, 4) if total_all else 0,
        "Accessions": ACCESSION_JOIN_SEP.join(sorted(total_accessions[mut])),
        "Historical_Count": h_count,
        "Historical_Frequency": round(h_count / total_historical, 4) if total_historical else 0,
        "Historical_Accessions": ACCESSION_JOIN_SEP.join(sorted(historical_accessions[mut])),
        "Novel_Count": n_count,
        "Novel_Frequency": round(n_count / total_novel, 4) if total_novel else 0,
        "Novel_Accessions": ACCESSION_JOIN_SEP.join(sorted(novel_accessions[mut])),
        "Presence": presence,
    })

freq_fieldnames = [
    "Mutation", "Total_Count", "Total_Frequency", "Accessions",
    "Historical_Count", "Historical_Frequency", "Historical_Accessions",
    "Novel_Count", "Novel_Frequency", "Novel_Accessions", "Presence",
]
with open(FREQUENCY_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=freq_fieldnames)
    writer.writeheader()
    writer.writerows(frequency_rows)

shared_n = sum(1 for r in frequency_rows if r["Presence"] == "Shared (Historical & Novel)")
hist_only_n = sum(1 for r in frequency_rows if r["Presence"] == "Historical/Native Range only")
novel_only_n = sum(1 for r in frequency_rows if r["Presence"] == "Novel/Introduced Range only")

print(f"Envelope mutation frequency sheet written to '{FREQUENCY_CSV}' ({len(frequency_rows)} distinct mutations)")
print(f"  Sequences -> Historical: {total_historical}, Novel: {total_novel}, Total: {total_all}")
print(f"  Shared mutations: {shared_n}")
print(f"  Historical-only mutations: {hist_only_n}")
print(f"  Novel-only mutations: {novel_only_n}")