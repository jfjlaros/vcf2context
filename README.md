# VCF to context converter
Convert a VCF file to two FASTA files, one containing the reference allele the
other one containing the alternative allele. Both variants are flanked by  the
reference sequence.

## Usage

    vcf2context input.vcf reference.fa ref.fa alt.fa

## Notes
Please note that only the reference sequence is used to make the flanking
regions. This means that if two variants are close together, you will not see
the neighbouring variant in the flanking region.
