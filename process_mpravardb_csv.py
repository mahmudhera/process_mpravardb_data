"""
In this script, we read a csv file downloaded from MPRAVarDb and 
extract the exact 200bp sequences for each row such that
- the logfc is a valid number
- the ref and alt alleles are single nucleotides
All columns are kept in the output, but
- description is dropped
- mpra_study field is made smaller (only the first author's name)
- the ref and alt sequences are added to the output, with the snp in the center and 300 bp on each side
final sequence length is 601 bp
"""


import pandas as pd
import argparse
from Bio import SeqIO

UPSTREAM_BP = 300
DOWNSTREAM_BP = 300


# Extract the exact 200bp sequences for each row, keeping SNP in the center, extra 200 bp on each side
def extract_sequence(row, hg19_index, hg38_index):
    chromosome = row['chr']
    ref_nucleotide = row['ref']
    alt_nucleotide = row['alt']
    
    # Convert chromosome to the correct format
    chromosome = 'chr' + str(int(chromosome))

    # check if chromosome exists in the index
    if chromosome not in hg19_index and row['genome'] == 'hg19':
        raise ValueError(f"Chromosome {chromosome} not found in HG19 index")
    if chromosome not in hg38_index and row['genome'] == 'hg38':
        raise ValueError(f"Chromosome {chromosome} not found in HG38 index")
    
    # get the entire reference sequence
    if row['genome'] == 'hg19':
        ref_seq = hg19_index[chromosome]
    elif row['genome'] == 'hg38':
        ref_seq = hg38_index[chromosome]
    else:
        raise ValueError(f"Unsupported genome: {row['genome']}")

    # verify that the position is valid
    if row['pos'] - UPSTREAM_BP - 1 < 0 or row['pos'] + DOWNSTREAM_BP > len(ref_seq):
        raise ValueError(f"Position {row['pos']} with upstream {UPSTREAM_BP} and downstream {DOWNSTREAM_BP} exceeds chromosome bounds for chromosome {chromosome}")

    # verify that ref nucleotide at that position is valid
    found_ref_nucleotide = ref_seq.seq[row['pos'] - 1]
    found_ref_nucleotide_before = ref_seq.seq[row['pos'] - 2]
    found_ref_nucleotide_after = ref_seq.seq[row['pos'] ]
    if found_ref_nucleotide != ref_nucleotide:
        print(f"Found ref nucleotide before: {found_ref_nucleotide_before}, at pos: {found_ref_nucleotide}, after: {found_ref_nucleotide_after}")
        print(f"Row info: chromosome: {chromosome}, pos: {row['pos']}, ref_nucleotide: {ref_nucleotide}, alt_nucleotide: {alt_nucleotide}, genome: {row['genome']}")
        raise ValueError(f"Ref nucleotide {ref_nucleotide} does not match reference genome at position {row['pos']} on chromosome {chromosome}")

    # extract the sequence around the SNP
    start = row['pos'] - UPSTREAM_BP - 1
    end = row['pos'] + DOWNSTREAM_BP
    ref_seq = str(ref_seq.seq[start:end])

    # create the alt sequence by replacing the ref nucleotide with the alt nucleotide
    alt_seq = ref_seq[:UPSTREAM_BP+1] + alt_nucleotide + ref_seq[UPSTREAM_BP+2:]
    return ref_seq, alt_seq

def parse_args():
    parser = argparse.ArgumentParser(description="Process MPRAVarDb CSV file")
    parser.add_argument("input_csv", help="Path to the input CSV file")
    parser.add_argument("output_csv", help="Path to save the output CSV file")
    parser.add_argument("hg19_path", help="Path to the HG19 FASTA file")
    parser.add_argument("hg38_path", help="Path to the HG38 FASTA file")
    return parser.parse_args()


def process_csv(input_csv, output_csv, hg19_path, hg38_path):
    # Load the CSV_csv into a DataFrame
    df = pd.read_csv(input_csv)
    
    # Drop the 'description' column
    df.drop(columns=['Description'], inplace=True)

    # Filter out rows where logfc is not a valid number
    df = df[pd.to_numeric(df['log2FC'], errors='coerce').notnull()]

    # Filter out rows where ref or alt alleles are not single nucleotides
    df = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)]
    
    # Make the 'mpra_study' field smaller (only the first author's name)
    df['MPRA_study'] = df['MPRA_study'].str.split('(').str[-1].str.split(' ').str[0]
    
    # Load the HG19 and HG38 FASTA files 
    hg19_index = SeqIO.to_dict(SeqIO.parse(hg19_path, "fasta"))
    hg38_index = SeqIO.to_dict(SeqIO.parse(hg38_path, "fasta"))
    df['ref_seq'], df['alt_seq'] = zip(*df.apply(lambda row: extract_sequence(row, hg19_index, hg38_index), axis=1))
    
    # Save the modified DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)


def main():    
    args = parse_args()
    process_csv(args.input_csv, args.output_csv, args.hg19_path, args.hg38_path)


if __name__ == "__main__":
    main()