# This script merges West Nile virus (WNV) Nextclade output with the NCBI
# sequence metadata CSV, extracts Envelope (env) amino-acid mutations, and
# writes a single Excel workbook (.xlsx) with five tabs:
#
#   1. Dataset Summary     - overview counts (sequences, ranges, mutations)
#   2. Mutation Frequencies - every envelope mutation, with counts/frequency
#                             overall and within each range group
#   3. Historical-Native Range - mutations found ONLY in the
#                             Historical/Native Range (Endemic) group
#   4. Novel-Introduced Range  - mutations found ONLY in the
#                             Novel/Introduced Range (Recently Emerged) group
#   5. Shared Mutations    - mutations found in BOTH groups
#
# It is adapted from an earlier Genome Detective XML-parsing script used for
# Rift Valley fever virus (RVFV) M-segment mutation extraction, updated for
# WNV's single-polyprotein genome (reads Nextclade's CSV export, merges by
# accession with version suffixes stripped from both sides since the two
# source files don't consistently include them, and keeps only the
# envelope gene's amino-acid substitutions).

import csv
import re
from collections import Counter, defaultdict

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

# ---- EDIT: input/output paths ----
NEXTCLADE_CSV = "C:/Users/Dell/Downloads/nextclade (4).csv"
METADATA_CSV  = "C:/Users/Dell/Downloads/Asia WNV project/Merged_sequences_Lineage2_classified.csv"
OUTPUT_XLSX   = "C:/Users/Dell/Downloads/Asia WNV project/WNV_Envelope_Mutation_Summary.xlsx"
NEXTCLADE_DELIMITER = ";"  # Nextclade's default CSV export delimiter
TARGET_GENE = "env"        # Envelope protein gene name as used by Nextclade

# Exact category strings used in the metadata's Range_Classification column
HISTORICAL_LABEL = "Historical/Native Range (Endemic)"
NOVEL_LABEL = "Novel/Introduced Range (Recently Emerged)"
# -----------------------------------

FONT_NAME = "Arial"
HEADER_FILL = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
HEADER_FONT = Font(name=FONT_NAME, bold=True, color="FFFFFF")
TITLE_FONT = Font(name=FONT_NAME, bold=True, size=14)
LABEL_FONT = Font(name=FONT_NAME, bold=True)
BODY_FONT = Font(name=FONT_NAME)


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


def mutation_sort_key(mut: str):
    match = re.search(r"(\d+)", mut)
    return int(match.group(1)) if match else 0


def style_header_row(ws, row_idx, ncols):
    for col in range(1, ncols + 1):
        cell = ws.cell(row=row_idx, column=col)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(horizontal="center", vertical="center")


def autosize_columns(ws, widths):
    for i, width in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = width


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
mutation_tracking = []  # (Range_Classification, [env substitution strings]) per sequence
matched, unmatched = 0, 0
skipped_no_env = 0

with open(NEXTCLADE_CSV, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter=NEXTCLADE_DELIMITER)
    for row in reader:
        seq_name = row.get("seqName", "").strip()
        acc = base_accession(seq_name)

        # Skip sequences that don't have the envelope gene: Nextclade's
        # cdsCoverage field lists coverage per gene; if 'env' is absent
        # (or has 0 coverage), the gene wasn't aligned/present.
        cds_coverage = parse_cds_coverage(row.get("cdsCoverage", ""))
        if cds_coverage.get(TARGET_GENE, 0.0) <= 0.0:
            skipped_no_env += 1
            continue

        meta = metadata.get(acc, {})
        if meta:
            matched += 1
        else:
            unmatched += 1

        # Nucleotide mutations
        nt_substitutions = row.get("substitutions", "").strip()
        total_nt_substitutions = int(row.get("totalSubstitutions", "0").strip() or 0)
        total_nt_deletions = int(row.get("totalDeletions", "0").strip() or 0)
        total_nt_insertions = int(row.get("totalInsertions", "0").strip() or 0)
        total_nucleotide_mutations = total_nt_substitutions + total_nt_deletions + total_nt_insertions

        # Amino-acid substitutions: filter to Envelope (env) only.
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

        range_classification = meta.get("Range_Classification", "")

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
            "Range_Classification": range_classification,
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
        mutation_tracking.append((range_classification, env_substitutions))

print(f"Nextclade rows processed: {len(records)}")
print(f"  skipped (no envelope gene): {skipped_no_env}")
print(f"  matched to metadata:   {matched}")
print(f"  no metadata match:     {unmatched}")

# ---- Step 3: Compute mutation frequency tables ----
total_all = len(mutation_tracking)
total_historical = sum(1 for rc, _ in mutation_tracking if rc == HISTORICAL_LABEL)
total_novel = sum(1 for rc, _ in mutation_tracking if rc == NOVEL_LABEL)

total_count = Counter()
historical_count = Counter()
novel_count = Counter()

for range_classification, muts in mutation_tracking:
    for mut in set(muts):  # count each mutation once per sequence
        total_count[mut] += 1
        if range_classification == HISTORICAL_LABEL:
            historical_count[mut] += 1
        elif range_classification == NOVEL_LABEL:
            novel_count[mut] += 1

frequency_rows = []
unique_endemic_rows = []
unique_epidemic_rows = []
shared_rows = []

for mut in sorted(total_count.keys(), key=mutation_sort_key):
    t_count = total_count[mut]
    h_count = historical_count.get(mut, 0)
    n_count = novel_count.get(mut, 0)
    t_freq = round(t_count / total_all, 4) if total_all else 0
    h_freq = round(h_count / total_historical, 4) if total_historical else 0
    n_freq = round(n_count / total_novel, 4) if total_novel else 0

    if h_count > 0 and n_count > 0:
        presence = "Shared (Historical & Novel)"
        shared_rows.append([mut, h_count, h_freq, n_count, n_freq, t_count, t_freq])
    elif h_count > 0:
        presence = "Unique to Historical/Native Range"
        unique_endemic_rows.append([mut, h_count, h_freq])
    elif n_count > 0:
        presence = "Unique to Novel/Introduced Range"
        unique_epidemic_rows.append([mut, n_count, n_freq])
    else:
        presence = "Unclassified range only"

    frequency_rows.append([mut, t_count, t_freq, h_count, h_freq, n_count, n_freq, presence])

print(f"Distinct envelope mutations: {len(frequency_rows)}")
print(f"  Unique Endemic (Historical-only): {len(unique_endemic_rows)}")
print(f"  Unique Epidemic (Novel-only):     {len(unique_epidemic_rows)}")
print(f"  Shared:                           {len(shared_rows)}")

# ---- Step 4: Build the workbook ----
wb = Workbook()

# --- Tab 1: Dataset Summary ---
ws = wb.active
ws.title = "Dataset Summary"
ws["A1"] = "WNV Envelope Mutation Analysis - Dataset Summary"
ws["A1"].font = TITLE_FONT
ws.merge_cells("A1:B1")

summary_items = [
    ("Nextclade sequences processed", len(records) + skipped_no_env),
    ("Sequences skipped (no envelope gene)", skipped_no_env),
    ("Sequences included in analysis", total_all),
    ("  - Historical/Native Range (Endemic)", total_historical),
    ("  - Novel/Introduced Range (Recently Emerged)", total_novel),
    ("  - Unclassified / Unknown range", total_all - total_historical - total_novel),
    ("Sequences matched to metadata", matched),
    ("Sequences with no metadata match", unmatched),
    ("", ""),
    ("Distinct envelope mutations (total)", len(frequency_rows)),
    ("Unique to Historical/Native Range", len(unique_endemic_rows)),
    ("Unique to Novel/Introduced Range", len(unique_epidemic_rows)),
    ("Shared mutations (Historical & Novel)", len(shared_rows)),
]

row_idx = 3
for label, value in summary_items:
    ws.cell(row=row_idx, column=1, value=label).font = LABEL_FONT if label and not label.startswith("  ") else BODY_FONT
    ws.cell(row=row_idx, column=2, value=value).font = BODY_FONT
    row_idx += 1

autosize_columns(ws, [42, 14])

# --- Tab 2: Mutation Frequencies ---
ws = wb.create_sheet("Mutation Frequencies")
headers = [
    "Mutation", "Total_Count", "Total_Frequency",
    "Historical_Count", "Historical_Frequency",
    "Novel_Count", "Novel_Frequency", "Presence",
]
ws.append(headers)
style_header_row(ws, 1, len(headers))
for r in frequency_rows:
    ws.append(r)
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = BODY_FONT
    row[2].number_format = "0%"
    row[4].number_format = "0%"
    row[6].number_format = "0%"
ws.freeze_panes = "A2"
autosize_columns(ws, [12, 12, 14, 16, 18, 12, 14, 26])

# --- Tab 3: Historical/Native Range (mutations unique to this group) ---
ws = wb.create_sheet("Historical-Native Range")
headers = ["Mutation", "Count", "Frequency (within Historical/Native Range group)"]
ws.append(headers)
style_header_row(ws, 1, len(headers))
for r in unique_endemic_rows:
    ws.append(r)
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = BODY_FONT
    row[2].number_format = "0%"
ws.freeze_panes = "A2"
autosize_columns(ws, [12, 10, 36])

# --- Tab 4: Novel/Introduced Range (mutations unique to this group) ---
ws = wb.create_sheet("Novel-Introduced Range")
headers = ["Mutation", "Count", "Frequency (within Novel/Introduced Range group)"]
ws.append(headers)
style_header_row(ws, 1, len(headers))
for r in unique_epidemic_rows:
    ws.append(r)
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = BODY_FONT
    row[2].number_format = "0%"
ws.freeze_panes = "A2"
autosize_columns(ws, [12, 10, 36])

# --- Tab 5: Shared Mutations ---
ws = wb.create_sheet("Shared Mutations")
headers = [
    "Mutation", "Historical_Count", "Historical_Frequency",
    "Novel_Count", "Novel_Frequency", "Total_Count", "Total_Frequency",
]
ws.append(headers)
style_header_row(ws, 1, len(headers))
for r in shared_rows:
    ws.append(r)
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = BODY_FONT
    row[2].number_format = "0%"
    row[4].number_format = "0%"
    row[6].number_format = "0%"
ws.freeze_panes = "A2"
autosize_columns(ws, [12, 14, 16, 14, 16, 12, 14])

wb.save(OUTPUT_XLSX)
print(f"Workbook written to '{OUTPUT_XLSX}'")