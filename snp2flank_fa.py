import pysam
import argparse

def extract_flanking_sequence(reference_file, vcf_file, output_file, flank_length):
    # vcf
    vcf = pysam.VariantFile(vcf_file)
    
    # ref genome
    reference = pysam.FastaFile(reference_file)
    
    with open(output_file, 'w') as output:
        for record in vcf:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]  # 

            # flanking seq.
            start = max(1, pos - flank_length)
            end = pos + len(ref) - 1 + flank_length
            sequence = reference.fetch(chrom, start - 1, end)  # pysam中的索引从0开始

            # output
            output.write(f">chrom={chrom}_pos={pos}_ref={ref}_alt={alt}\n")
            output.write(sequence + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract flanking sequences from VCF file and reference genome.')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('-v', '--vcf', required=True, help='VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-l', '--flank_length', type=int, default=50, help='Flanking sequence length (default: 50)')
    args = parser.parse_args()

    extract_flanking_sequence(args.reference, args.vcf, args.output, args.flank_length)
