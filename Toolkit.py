# GC Content and AT/GC Ratio Calculator
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import Entrez, SeqIO
from io import StringIO

Entrez.email = "dhruv130903@gmail.com"

def fetch_nucleotide_sequence(accession_id):
    """Fetch nucleotide sequence from NCBI using accession ID"""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        if not fasta_data.strip():
            return None

        # Parse FASTA to get SeqRecord
        fasta_io = StringIO(fasta_data)
        seq_record = SeqIO.read(fasta_io, "fasta")
        return seq_record

    except Exception as e:
        st.error(f"❌ Error: {str(e)}")
        return None

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

# ------------------- Protein Translation -------------------
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

def clean_sequence(seq):
    """Remove FASTA headers, spaces, numbers, and non-ATGC characters"""
    seq = seq.upper()
    seq = re.sub(r">.*\n", "", seq)     # remove FASTA header
    seq = re.sub(r"[^ATGC]", "", seq)   # keep only ATGC
    return seq

def codon_usage(sequence):
    """Calculate codon usage from DNA sequence"""
    sequence = clean_sequence(sequence)
    
    codon_count = {codon:0 for codon in CODON_TABLE.keys()}
    
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]
        if codon in codon_count:
            codon_count[codon] += 1
    
    # Convert to DataFrame
    df = pd.DataFrame([
        {"Codon": codon,
         "Amino Acid": CODON_TABLE[codon],
         "Count": count}
        for codon, count in codon_count.items()
    ])
    
    total_codons = df["Count"].sum()
    if total_codons > 0:
        df["Frequency (%)"] = (df["Count"] / total_codons * 100).round(2)
    else:
        df["Frequency (%)"] = 0
    
    return df.sort_values(by="Codon")

# ------------------- Protein Property Calculator -------------------
def clean_fasta_sequence(fasta_str: str) -> str:
    lines = fasta_str.strip().splitlines()
    seq = "".join(line for line in lines if not line.startswith(">"))
    return seq.upper().replace(" ", "")

def protein_property_calculator(fasta_input: str) -> dict:
    sequence = clean_fasta_sequence(fasta_input)
    analysed_seq = ProteinAnalysis(sequence)
    return {
        "Length": len(sequence),
        "Molecular Weight (Da)": round(analysed_seq.molecular_weight(), 2),
        "Theoretical pI": round(analysed_seq.isoelectric_point(), 2),
        "Amino Acid Composition (%)": {aa: round(val*100, 2) for aa,val in analysed_seq.get_amino_acids_percent().items()},
        "Aromaticity": round(analysed_seq.aromaticity(), 3),
        "Instability Index": round(analysed_seq.instability_index(), 2),
        "GRAVY (Hydropathicity)": round(analysed_seq.gravy(), 3),
        "Extinction Coefficient (1 mg/ml, 280 nm)": analysed_seq.molar_extinction_coefficient()
    }

# ------------------- MAIN APP -------------------
st.title("🔬 Bioinformatics Toolkit")

tool = st.sidebar.radio(
    "Choose a Tool",
    ["Nucleotide Sequence Retrieval", "ORF Finder", "GC Content Calculator", "DNA → Protein Translator", "Codon Usage Analyzer", "Protein Property Calculator"]
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
            
elif tool == "Codon Usage Analyzer":
    st.subheader("📊 Codon Usage Analyzer")
    seq = st.text_area("Paste DNA sequence (raw or FASTA):", key="codon")
    if st.button("Analyze Codon Usage"):
        if seq:
            df = codon_usage(seq)   # ✅ now directly use the dataframe

            if not df.empty:
                st.success("✅ Codon usage calculated successfully!")
                st.dataframe(df)

                # Plot codon usage frequencies
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.bar(df["Codon"], df["Frequency (%)"])
                ax.set_xlabel("Codon")
                ax.set_ylabel("Frequency (%)")
                ax.set_title("Codon Usage Frequency")
                plt.xticks(rotation=90)
                st.pyplot(fig)
            else:
                st.warning("No codons found. Check your sequence.")
        else:
            st.warning("Please paste a DNA sequence first.")
    
elif tool == "Protein Property Calculator":
    st.subheader("🧬 Protein Property Calculator")
    fasta_input = st.text_area("Paste Protein Sequence (FASTA or plain):", key="protein_prop")
    if st.button("Calculate Properties"):
        if fasta_input:
            try:
                results = protein_property_calculator(fasta_input)
                st.success("✅ Properties calculated successfully!")
                st.write(f"**Sequence Length:** {results['Length']}")
                st.write(f"**Molecular Weight:** {results['Molecular Weight (Da)']} Da")
                st.write(f"**Theoretical pI:** {results['Theoretical pI']}")
                st.write(f"**Aromaticity:** {results['Aromaticity']}")
                st.write(f"**Instability Index:** {results['Instability Index']}")
                st.write(f"**GRAVY (Hydropathicity):** {results['GRAVY (Hydropathicity)']}")
                st.write(f"**Extinction Coefficient:** {results['Extinction Coefficient (1 mg/ml, 280 nm)']}")
                st.subheader("Amino Acid Composition (%)")
                st.table(results["Amino Acid Composition (%)"])
            except Exception as e:
                st.error(f"❌ Error: {e}")
                
elif tool == "Nucleotide Sequence Retrieval":
    st.subheader("📥 Nucleotide Sequence Retrieval (NCBI)")
    accession_id = st.text_input("Enter NCBI Accession ID (e.g., NM_007294):")
    if st.button("Fetch Sequence"):
        if accession_id:
            seq_record = fetch_nucleotide_sequence(accession_id)
            if seq_record:
                st.success(f"✅ Sequence retrieved successfully: {seq_record.id}")
                st.write(f"**Description:** {seq_record.description}")
                st.write(f"**Length:** {len(seq_record.seq)} bp")
                st.text_area("Sequence (FASTA format)", 
                             f">{seq_record.description}\n{seq_record.seq}", height=200)
            else:
                st.warning("⚠️ No sequence found for this ID.")
        else:
            st.warning("⚠️ Please enter a valid accession ID.")
