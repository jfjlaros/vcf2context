#!/usr/bin/env python

import argparse

import vcf
from Bio import SeqIO, SeqRecord, Seq

from . import usage, version


def vcf2context(input_handle, reference_handle, output_handles, flank=9):
    """
    Convert a VCF file to two FASTA files, one containing the reference allele
    the other one containing the alternative allele. Both variants are flanked
    by  the reference sequence.

    :arg stream input_handle: Input file in VCF format.
    :arg stream reference_handle: Reference file in FASTA format
    :arg list(stream) output_handles: Output files
    :arg int flank: Length of the flanking region.
    """
    reference = dict(map(lambda x: (x.id, x.seq),
        SeqIO.parse(reference_handle, 'fasta')))

    for record in vcf.Reader(input_handle):
        if record.is_snp:
            chrom = reference[record.CHROM]
            left_flank = chrom[record.start - flank:record.start]
            right_flank = chrom[record.start + 1: record.start + flank + 1]
            ref = '{}{}{}'.format(left_flank, record.REF, right_flank)
            alt = '{}{}{}'.format(left_flank, record.ALT[0], right_flank)
            snp_id = '{}:{}'.format(record.CHROM, record.start)
            SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(ref), snp_id, '', ''),
                output_handles[0], 'fasta')
            SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(alt), snp_id, '', ''),
                output_handles[1], 'fasta')


def main():
    """
    Command line argument parsing.
    """
    parser = argparse.ArgumentParser(description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='input file in VCF format')
    parser.add_argument('reference_handle', metavar='REFERENCE',
        type=argparse.FileType('r'), help='reference file in FASTA format')
    parser.add_argument('output_handles', metavar='OUTPUT',
        type=argparse.FileType('w'), nargs=2, help='output file')
    parser.add_argument('-f', dest='flank', type=int, default=9,
        help='length of the flanking region (%(type)s default=%(default)s)')
    parser.add_argument('-v', action='version', version=version(parser.prog))

    try:
        arguments = parser.parse_args()
    except IOError as error:
        parser.error(error)

    try:
        vcf2context(**dict((k, v) for k, v in vars(arguments).items()
            if k not in ('func', 'subcommand')))
    except ValueError as error:
        parser.error(error)


if __name__ == '__main__':
    main()
