# This script extracts the Envelope (env) gene's coding nucleotide sequence
# from each WNV genome, using Nextclade's PER-QUERY-SEQUENCE genome
# annotation GFF (--output-annotation-gff / the "genome annotation of query
# sequences" export). Unlike the reference genome annotation, this file
# gives exact env CDS coordinates in EACH sequence's own coordinate frame -
# i.e. it already accounts for that specific sequence's indels/truncations
# relative to the reference - so slicing directly from the raw (unaligned)
# genome FASTA using these coordinates gives a correct, in-frame env CDS
# for every sequence without needing a whole-genome alignment step at all.
#
# Output is an UNALIGNED nucleotide FASTA of just the env CDS per sequence.
# This still needs a multiple sequence alignment (e.g. MAFFT or Nextclade/
# Nextalign against just the env reference) and a phylogenetic tree before
# it can be used with HyPhy FEL/MEME - but extracting the correct per-
# sequence CDS first, then aligning only ~1500nt instead of the whole
# ~11kb genome, is simpler and less error-prone.
#
# Includes the same QC check as before: any sequence whose extracted env
# CDS contains a premature stop codon is flagged and excluded, since that
# usually indicates a frameshift/assembly artifact that would break
# downstream codon-based analysis.

import re

# paths
GFF_PATH       = "C:/Users/Dell/Downloads/nextclade.gff"
RAW_FASTA      = "C:/Users/Dell/Downloads/Asia WNV project/Merged_sequences_Lineage2.fasta"
OUTPUT_FASTA   = "C:/Users/Dell/Downloads/WNV_env_CDS_unaligned.fasta"
QC_REPORT_TXT  = "C:/Users/Dell/Downloads/Asia WNV project/WNV_env_CDS_QC_report.txt"
TARGET_GENE = "env"


CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

COMPLEMENT = str.maketrans("ACGTacgtNn-", "TGCAtgcaNn-")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def parse_fasta(path):
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


def find_premature_stop(cds_nt: str):
    """Translate (ignoring gap/ambiguous codons) and check for a stop codon
    before the final codon. Returns 1-based codon position of the first
    premature stop, or None if the sequence looks clean."""
    codons = [cds_nt[i:i + 3] for i in range(0, len(cds_nt) - len(cds_nt) % 3, 3)]
    n_codons = len(codons)
    for i, codon in enumerate(codons):
        if "-" in codon or "N" in codon.upper():
            continue
        aa = CODON_TABLE.get(codon.upper())
        if aa == "*" and i < n_codons - 1:
            return i + 1
    return None


# Parse per-sequence env CDS coordinates from the GFF
# Also capture the GFF 'phase' field (column 8): for CDS fragments that
# don't start on a natural codon boundary (common for partial/truncated
# env-only submissions in this dataset), phase indicates how many bases
# to skip from the start of the feature to reach the first complete codon.
# Ignoring this would silently shift the reading frame for those fragments.
env_coords = {}  # accession (with version) -> list of (start, end, strand, phase)

with open(GFF_PATH, encoding="utf-8") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = cols
        if ftype != "CDS":
            continue
        attr_dict = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
        if attr_dict.get("gene", "") == TARGET_GENE:
            phase_int = int(phase) if phase in ("0", "1", "2") else 0
            is_partial = "truncated-left" in attr_dict or "truncated-right" in attr_dict
            env_coords.setdefault(seqid, []).append(
                (int(start), int(end), strand, phase_int, is_partial)
            )

print(f"Sequences with '{TARGET_GENE}' CDS annotated in GFF: {len(env_coords)}")

#  Extract, QC, and write the env CDS for each sequence 
raw_sequences = {}
for header, seq in parse_fasta(RAW_FASTA):
    raw_sequences[header.split()[0]] = seq

written = []
excluded = []
partial_flagged = []

for seqid, segments in env_coords.items():
    if seqid not in raw_sequences:
        excluded.append((seqid, "accession not found in raw genome FASTA"))
        continue

    genome_seq = raw_sequences[seqid]
    segments_sorted = sorted(segments, key=lambda s: s[0])
    strand = segments_sorted[0][2]
    is_partial = any(seg[4] for seg in segments_sorted)

    pieces = []
    ok = True
    for i, (start, end, seg_strand, phase, _seg_partial) in enumerate(segments_sorted):
        if end > len(genome_seq):
            excluded.append((seqid, f"CDS end {end} exceeds sequence length {len(genome_seq)}"))
            ok = False
            break
        piece = genome_seq[start - 1:end]
        if i == 0 and phase:
            # Skip 'phase' bases so this fragment starts on a codon boundary
            piece = piece[phase:]
        pieces.append(piece)
    if not ok:
        continue

    cds_nt = "".join(pieces)
    if strand == "-":
        cds_nt = revcomp(cds_nt)

    if len(cds_nt) < 30:
        excluded.append((seqid, f"env CDS too short after extraction ({len(cds_nt)} nt)"))
        continue

    # Trim any trailing partial codon rather than rejecting outright -
    # legitimate for partial/truncated fragment submissions.
    remainder = len(cds_nt) % 3
    if remainder:
        cds_nt = cds_nt[:-remainder]

    stop_pos = find_premature_stop(cds_nt)
    if stop_pos is not None:
        excluded.append((seqid, f"premature stop codon at codon {stop_pos} (likely frameshift/assembly artifact)"))
        continue

    written.append((seqid, cds_nt))
    if is_partial:
        partial_flagged.append(seqid)

with open(OUTPUT_FASTA, "w", encoding="utf-8") as out:
    for seqid, cds_nt in written:
        out.write(f">{seqid}\n")
        for i in range(0, len(cds_nt), 70):
            out.write(cds_nt[i:i + 70] + "\n")

with open(QC_REPORT_TXT, "w", encoding="utf-8") as out:
    out.write("Envelope CDS extraction QC report\n")
    out.write(f"Sequences written:  {len(written)}\n")
    out.write(f"  of which partial/truncated fragments: {len(partial_flagged)}\n")
    out.write(f"Sequences excluded: {len(excluded)}\n\n")
    if partial_flagged:
        out.write("Partial/truncated fragments included (verify these suit your analysis):\n")
        for seqid in partial_flagged:
            out.write(f"  {seqid}\n")
        out.write("\n")
    out.write("Excluded sequences:\n")
    for seqid, reason in excluded:
        out.write(f"{seqid}\t{reason}\n")

print(f"Sequences written to unaligned env CDS FASTA: {len(written)}")
print(f"  of which partial/truncated fragments: {len(partial_flagged)}")
print(f"Sequences excluded (see QC report): {len(excluded)}")
print(f"Unaligned envelope CDS FASTA written to '{OUTPUT_FASTA}'")
print(f"QC report written to '{QC_REPORT_TXT}'")
print()
print("NEXT STEP: this FASTA is unaligned. Align it (e.g. MAFFT, or run it")
print("back through Nextclade/Nextalign against just the env reference) and")
print("build a phylogenetic tree before using it with HyPhy FEL/MEME.")