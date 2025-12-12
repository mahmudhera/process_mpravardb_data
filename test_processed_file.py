"""
test the processed file. Check that:
- ref_seq's 301st base matches ref
- alt_seq's 301st base matches alt
- length of ref_seq and alt_seq is 601 bp
- except for the 301st base, ref_seq and alt_seq are identical
"""

import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Test processed MPRAVarDb CSV file")
    parser.add_argument("--input_csv", required=True, help="Path to the processed MPRAVarDb CSV file")
    return parser.parse_args()

def main():
    args = parse_args()
    df = pd.read_csv(args.input_csv)

    for index, row in df.iterrows():
        ref_seq = row['ref_seq']
        alt_seq = row['alt_seq']
        ref_nucleotide = row['ref']
        alt_nucleotide = row['alt']

        # Check that the length of ref_seq and alt_seq is 601 bp
        assert len(ref_seq) == 601, f"Row {index}: ref_seq length is {len(ref_seq)}, expected 601"
        assert len(alt_seq) == 601, f"Row {index}: alt_seq length is {len(alt_seq)}, expected 601"

        # Check that the 301st base of ref_seq matches ref
        assert ref_seq[300] == ref_nucleotide, f"Row {index}: ref_seq[300] is {ref_seq[300]}, expected {ref_nucleotide}"

        # Check that the 301st base of alt_seq matches alt
        assert alt_seq[300] == alt_nucleotide, f"Row {index}: alt_seq[300] is {alt_seq[300]}, expected {alt_nucleotide}"

        # Check that except for the 301st base, ref_seq and alt_seq are identical
        assert ref_seq[:300] + ref_seq[301:] == alt_seq[:300] + alt_seq[301:], f"Row {index}: ref_seq and alt_seq differ outside the 301st base"

    print("All tests passed!")

if __name__ == "__main__":
    main()