#!/usr/bin/python

import pysam
import argparse

def generate_flanking_sequences(reference_genome_file, vcf_file, flanking_length, output_file):
    # Load the reference genome sequence
    reference = pysam.FastaFile(reference_genome_file)

    # Open the VCF file
    vcf = pysam.VariantFile(vcf_file)

    with open(output_file, 'w') as output:
        for record in vcf:
            chrom = record.contig
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]  # Assuming there's only one alternate allele

            # Determine the start and end positions for flanking sequences
            start = max(0, pos - flanking_length)
            end = min(reference.get_reference_length(chrom), pos + flanking_length)

            # Extract the flanking sequences
            left_flank = reference.fetch(chrom, start, pos)
            right_flank = reference.fetch(chrom, pos + len(ref), end)

            # Create the output string with the custom format
            variant_name = f"Chr{chrom}_{pos}_{ref}_{alt}"
            output_line = f"{variant_name}\t{chrom}\t{left_flank}[{ref}/{alt}]{right_flank}\n"

            output.write(output_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate flanking sequences for variants in a VCF file.")
    parser.add_argument("-r", "--reference_genome", required=True, help="Path to the reference genome file (FASTA format).")
    parser.add_argument("-v", "--vcf_file", required=True, help="Path to the VCF file.")
    parser.add_argument("-l", "--flanking_length", type=int, required=True, help="Length of flanking sequences.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file.")

    args = parser.parse_args()

    generate_flanking_sequences(args.reference_genome, args.vcf_file, args.flanking_length, args.output_file)
