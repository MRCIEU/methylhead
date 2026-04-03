
import argparse
from pathlib import Path
import itertools

from methylseqextractor import (
  MethylSeqDataset,
  ClonalFlipCounter,
  Utils
)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("bam", help="Path to the sorted BAM file")
    parser.add_argument("fasta", help="Path to the reference FASTA file")
    parser.add_argument("csv", help="Path to the output CSV file")
    parser.add_argument(
        "--size", 
        type=int, 
        default=100, 
        help="Size of the genomic windows in base pairs (default: 100)"
    )
    parser.add_argument(
        "--chunk", 
        type=int, 
        default=5000, 
        help="Number of rows to buffer before writing (default: 5000)"
    )

    args = parser.parse_args()

    bam_file = args.bam
    fasta_file = args.fasta
    csv_file = args.csv

    if not Path(bam_file).exists():
        print(f"Error: Input file '{bam_file}' does not exist.")
        return

    if not Path(fasta_file).exists():
        print(f"Error: Reference file '{fasta_file}' does not exist.")
        return
    
    print(f"Processing {bam_file}...")
    
    dataset = MethylSeqDataset(bam_file, fasta_file)
    flips = ClonalFlipCounter(dataset, size=args.size, min_depth=10)
    chromosomes = ["chr"+str(i) for i in range(1,23)] + ["chrX"]
    iterators = [flips.calculate(chrom) for chrom in chromosomes]
    iterator = itertools.chain.from_iterable(iterators)
    Utils.to_csv(iterator, csv_file, chunk_size=args.chunk)
    
    print(f"Done! Data saved to {csv_file}")

if __name__ == "__main__":
    main()

