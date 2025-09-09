# GC Content and AT/GC Ratio Calculator
import streamlit as st
import pandas as pd

# ORF Finder Functions
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence"""
    complement = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return "".join(complement.get(base, base) for base in seq[::-1])

def find_orfs(sequence, min_length=100):
    sequence = sequence.upper().replace("\n","").replace(" ","")
    orfs = []

    def scan(seq, strand):
        seq_len = len(seq)
        for frame in range(3):
            i = frame
            while i < seq_len - 2:
                codon = seq[i:i+3]
                if codon == START_CODON:
                    for j in range(i+3, seq_len-2, 3):
                        stop_codon = seq[j:j+3]
                        if stop_codon in STOP_CODONS:
                            orf_seq = seq[i:j+3]
                            if len(orf_seq) >= min_length:
                                orfs.append((strand, frame+1, i, j+3, orf_seq))
                            break
                i += 1  # shift by 1 to catch overlapping ORFs

    # Scan forward strand
    scan(sequence, '+')

    # Scan reverse strand
    rev_seq = reverse_complement(sequence)
    scan(rev_seq, '-')

    return orfs

# ------------------- GC CONTENT -------------------

def calculate_gc_content(sequence):
    """Takes a DNA sequence and calculates GC content and AT/GC ratio."""

    # Step 1: Clean sequence properly
    sequence = sequence.strip()   # remove extra spaces
    sequence = sequence.upper()   # convert to uppercase

    if sequence.startswith(">"):  # if FASTA header is present
        sequence = "\n".join(sequence.split("\n")[1:])  # drop the header line

    # Now remove spaces and newlines safely
    sequence = sequence.replace(" ", "").replace("\n", "")

    # Step 2: Count nucleotides
    A = sequence.count("A")
    T = sequence.count("T")
    G = sequence.count("G")
    C = sequence.count("C")
    total = A + T + G + C

    # Step 3: Avoid division errors
    if total == 0:
        st.error("⚠️ Error: Empty or invalid sequence!")
        return

    # Step 4: Calculations
    gc_content = ((G + C) / total) * 100
    at_content = ((A + T) / total) * 100
    at_gc_ratio = (A + T) / (G + C) if (G + C) != 0 else "Infinity"

    # Step 5: Display results
    st.success("✅ Results:")
    st.write(f"**Sequence length:** {total}")
    st.write(f"**A:** {A}, **T:** {T}, **G:** {G}, **C:** {C}")
    st.write(f"**GC Content:** {gc_content:.2f}%")
    st.write(f"**AT Content:** {at_content:.2f}%")
    st.write(f"**AT/GC Ratio:** {at_gc_ratio}")
    
# ------------------- Protein Translation -------------------
# Codon table
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_dna(sequence):
    """Convert DNA to protein, works with raw or FASTA input"""
    if not sequence:
        return ""

    sequence = sequence.strip()
    
    # Remove FASTA header if present
    if sequence.startswith(">"):
        lines = sequence.split("\n")
        sequence = "".join(lines[1:])  # join all lines after header

    # Clean sequence: remove spaces, newlines, convert to uppercase
    sequence = sequence.replace(" ", "").replace("\n", "").upper()

    protein = ""
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]
        amino_acid = CODON_TABLE.get(codon, 'X')
        protein += amino_acid

    return protein

# ------------------- MAIN APP -------------------
st.title("🔬 Bioinformatics Toolkit")

tool = st.sidebar.selectbox(
    "Choose a Tool",
    ["ORF Finder", "GC Content Calculator", "DNA → Protein Translator"]
)

if tool == "ORF Finder":
    st.subheader("🔎 ORF Finder Tool")
    seq = st.text_area("Paste DNA sequence (raw or FASTA):", key="orf")
    if st.button("Find ORFs"):
        if seq:
            orfs = find_orfs(seq)
            if orfs:
                st.success(f"Found {len(orfs)} ORFs")

                # Build results table
                data = []
                for strand, frame, start, end, orf_seq in orfs:
                    data.append({
                        "Strand": strand,
                        "Frame": frame,
                        "Start": start,
                        "End": end,
                        "Length (bp)": len(orf_seq),
                        "Sequence": orf_seq
                    })

                df = pd.DataFrame(data)
                st.dataframe(df)   # 👈 interactive table output
            else:
                st.warning("No ORFs found.")
        else:
            st.warning("Please paste a DNA sequence first.")

elif tool == "GC Content Calculator":
    st.subheader("🧬 GC Content Calculator")
    seq = st.text_area("Paste DNA sequence (raw or FASTA):", key="gc")
    if st.button("Calculate GC Content"):
        if seq:
            calculate_gc_content(seq)
        else:
            st.warning("Please paste a DNA sequence first.")
            
elif tool == "DNA → Protein Translator":
    st.subheader("🧬 DNA → Protein Translation")
    seq_input = st.text_area("Paste DNA sequence (raw or FASTA):", key="protein")
    if st.button("Translate DNA"):
        if seq_input:
            protein_seq = translate_dna(seq_input)
            st.success("✅ Translation Complete")
            st.text_area("Protein Sequence", protein_seq, height=200)
        else:
            st.warning("Paste a DNA sequence first.")