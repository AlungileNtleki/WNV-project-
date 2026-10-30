#!/usr/bin/env python3
"""
West Nile Virus Envelope Protein Translation Script
Translates DNA sequences to protein and identifies V188I and V488I mutations

Requirements:
- Python 3.6+
- BioPython (install with: pip install biopython)

Usage:
python3 translate_wnv_sequences.py <input_fasta> <output_fasta>

Example:
python3 translate_wnv_sequences.py WNV_env_CDS_aligned.fasta protein_output.fasta
"""

import sys
from Bio import SeqIO

# ============================================================================
#  DEFINE THE CODON TRANSLATION TABLE
# ============================================================================

CODON_TABLE = {
    'ATG': 'M', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TTT': 'F', 'TTC': 'F',
    'TGT': 'C', 'TGC': 'C',
    'TGG': 'W',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y',
    'AAT': 'N', 'AAC': 'N',
    'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TAA': '*', 'TAG': '*', 'TGA': '*',
}


# ============================================================================
#  DEFINE THE TRANSLATION FUNCTION
# ============================================================================

def translate_dna_to_protein(dna_sequence):
    """Translate a DNA sequence to protein"""
    dna_clean = dna_sequence.replace('-', '').replace(' ', '').upper()
    protein = []

    for i in range(0, len(dna_clean) - 2, 3):
        codon = dna_clean[i:i + 3]
        if 'N' in codon or len(codon) < 3:
            protein.append('X')
        else:
            protein.append(CODON_TABLE.get(codon, 'X'))

    return ''.join(protein)


# ============================================================================
# MAIN SCRIPT
# ============================================================================

def main():
    """Main function to process sequences"""

    if len(sys.argv) < 2:
        print("=" * 80)
        print("West Nile Virus Envelope Protein Translation Script")
        print("=" * 80)
        print("\nUsage:")
        print("  python translate_wnv_sequences.py <input_fasta> [output_fasta]")
        print("\nArguments:")
        print("  input_fasta   : Path to your aligned DNA FASTA file")
        print("  output_fasta  : (Optional) Path to save protein sequences")
        print("\nExample:")
        print("  python translate_wnv_sequences.py WNV_env_CDS_aligned.fasta protein.fasta")
        print("=" * 80)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    try:
        sequences = {}
        print("\n" + "=" * 80)
        print("READING DNA SEQUENCES")
        print("=" * 80)
        print(f"\nReading from: {input_file}\n")

        for record in SeqIO.parse(input_file, "fasta"):
            dna_seq = str(record.seq)
            protein_seq = translate_dna_to_protein(dna_seq)
            sequences[record.id] = {'dna': dna_seq, 'protein': protein_seq}
            print(f"OK {record.id:20} | DNA: {len(dna_seq.replace('-', '')):5} bp -> Protein: {len(protein_seq):3} aa")

        print(f"\nOK Total sequences translated: {len(sequences)}\n")

    except FileNotFoundError:
        print(f"\nERROR: File not found: {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)

    # ========================================================================
    # EXTRACT POSITIONS 188 AND 488
    # ========================================================================

    print("=" * 80)
    print("EXTRACTING POSITIONS 188 & 488")
    print("=" * 80)
    print()

    mutations = {188: {}, 488: {}}

    print(f"{'Sequence ID':<25} {'Position 188':<15} {'Position 488':<15} {'Notes':<20}")
    print("-" * 75)

    for seq_id, data in sequences.items():
        protein = data['protein']

        aa188 = protein[187] if len(protein) > 187 else 'N/A'
        aa488 = protein[487] if len(protein) > 487 else 'N/A'

        notes = []
        if aa188 == 'I':
            notes.append("V188I")
        if aa488 == 'I':
            notes.append("V488I")

        note_str = ", ".join(notes) if notes else "Wild-type"

        print(f"{seq_id:<25} {aa188:<15} {aa488:<15} {note_str:<20}")

        if aa188 not in mutations[188]:
            mutations[188][aa188] = []
        mutations[188][aa188].append(seq_id)

        if aa488 not in mutations[488]:
            mutations[488][aa488] = []
        mutations[488][aa488].append(seq_id)

    # ========================================================================
    # SUMMARY STATISTICS
    # ========================================================================

    print("\n" + "=" * 80)
    print("SUMMARY: AMINO ACID FREQUENCIES")
    print("=" * 80)

    print("\nPosition 188:")
    print("-" * 40)
    for aa in sorted(mutations[188].keys(), key=lambda x: len(mutations[188][x]), reverse=True):
        count = len(mutations[188][aa])
        percent = (count / len(sequences)) * 100
        print(f"  {aa:3} : {count:4} sequences ({percent:5.1f}%)")
        if count <= 5 and aa != 'N/A':
            examples = mutations[188][aa][:3]
            print(f"       Examples: {', '.join(examples)}")

    print("\nPosition 488:")
    print("-" * 40)
    for aa in sorted(mutations[488].keys(), key=lambda x: len(mutations[488][x]), reverse=True):
        count = len(mutations[488][aa])
        percent = (count / len(sequences)) * 100
        print(f"  {aa:3} : {count:4} sequences ({percent:5.1f}%)")
        if count <= 5 and aa != 'N/A':
            examples = mutations[488][aa][:3]
            print(f"       Examples: {', '.join(examples)}")

    # ========================================================================
    # IDENTIFY V188I AND V488I
    # ========================================================================

    print("\n" + "=" * 80)
    print("MUTATION IDENTIFICATION: V188I & V488I")
    print("=" * 80)

    v188i_strains = mutations[188].get('I', [])
    v488i_strains = mutations[488].get('I', [])

    print(f"\nV188I (Valine -> Isoleucine at position 188):")
    if v188i_strains:
        print(f"  FOUND in {len(v188i_strains)} strain(s):")
        for strain in v188i_strains:
            print(f"    - {strain}")
    else:
        print("  NOT found in any sequences")

    print(f"\nV488I (Valine -> Isoleucine at position 488):")
    if v488i_strains:
        print(f"  FOUND in {len(v488i_strains)} strain(s):")
        for strain in v488i_strains:
            print(f"    - {strain}")
    else:
        print("  NOT found in any sequences")

    # ========================================================================
    # SAVE PROTEIN SEQUENCES TO FILE
    # ========================================================================

    if output_file:
        print("\n" + "=" * 80)
        print("SAVING PROTEIN SEQUENCES")
        print("=" * 80)

        try:
            with open(output_file, 'w') as f:
                for seq_id, data in sequences.items():
                    protein = data['protein']
                    f.write(f">{seq_id}\n")
                    for i in range(0, len(protein), 60):
                        f.write(protein[i:i + 60] + '\n')

            print(f"\nOK Protein sequences saved to: {output_file}\n")

        except Exception as e:
            print(f"\nERROR writing output file: {e}\n")

    # ========================================================================
    # SAVE RESULTS TO CSV
    # ========================================================================

    csv_file = input_file.replace('.fasta', '_translation_results.csv').replace('.fa', '_translation_results.csv')

    try:
        with open(csv_file, 'w') as f:
            f.write("Sequence_ID,DNA_Length,Protein_Length,Position_188,Position_488,V188I,V488I\n")

            for seq_id, data in sequences.items():
                protein = data['protein']
                dna_len = len(data['dna'].replace('-', ''))
                aa188 = protein[187] if len(protein) > 187 else 'N/A'
                aa488 = protein[487] if len(protein) > 487 else 'N/A'
                v188i = "YES" if aa188 == 'I' else "NO"
                v488i = "YES" if aa488 == 'I' else "NO"

                f.write(f"{seq_id},{dna_len},{len(protein)},{aa188},{aa488},{v188i},{v488i}\n")

        print(f"OK Results also saved to CSV: {csv_file}\n")

    except Exception as e:
        print(f"\nWARNING: Could not save CSV: {e}\n")

    print("=" * 80)
    print("OK TRANSLATION COMPLETE!")
    print("=" * 80)

    print("\nNext steps:")
    print("1. Upload protein FASTA to ColabFold")
    print("2. Use protein sequence to generate 3D structure")
    print("3. Map positions 188 and 488 onto structure\n")


if __name__ == "__main__":
    main()