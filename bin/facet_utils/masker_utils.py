#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

from pathlib import Path
from shutil import copyfileobj
from os import listdir, remove

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .defaults import ProgDefaults as dv


# driver function for masker module
def masker_driver(blast_outdir, args):
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]

    # index the genome sequence to access information w/o taking up large chunks of memory
    genome_seqs = SeqIO.index(args.genome, 'fasta')
    outfmt_specs = list(args.outfmt.keys())
    if args.verbose:
        print()
    blast_outfiles = ["%s/%s" % (blast_outdir, i) for i in listdir(blast_outdir)]
    for aln_file in blast_outfiles:
        # gets the tig_id from the end of the
        tig_id = "".join(Path(aln_file).suffixes)[1:]
        if args.verbose:
            print("Adding BLAST alignment ranges to coverage list for %s...." % tig_id)
        curr_seq = genome_seqs[tig_id].seq
        curr_description = genome_seqs[tig_id].description

        tig_len = len(curr_seq)
        contig_coverage = [0] * tig_len  # generates a list of zeroes as long as the tig

        with open(aln_file, 'r') as curr_file:
            while True:
                line = curr_file.readline()
                if not line:
                    break

                line = line.rstrip().split('\t')

                line[sstart] = dv.BLASTN_TYPES["sstart"](line[sstart])
                line[send] = dv.BLASTN_TYPES["send"](line[send])

                # switches minus hits so they're [min,max]
                if line[sstart] > line[send]:
                    line[sstart], line[send] = min(line[send], line[sstart]), max(line[send], line[sstart])

                # adds ranges from combined hits to the tic_coverage list #
                contig_coverage[line[sstart] - 1] += 1
                if line[send] < tig_len:  # CHANGED OPERATOR FROM != TO <
                    contig_coverage[line[send]] -= 1
        remove(aln_file)
        if args.verbose:
            print("Converting operational list to coverage list....")
        # converts operational list to actual coverage list #
        for basepair in range(1, tig_len):
            contig_coverage[basepair] += contig_coverage[basepair - 1]
        # masks each contig in the genome and writes a file
        if args.verbose:
            print("Masking %s...." % tig_id)
        curr_seq = Seq(list_get_masked_genome_seq(curr_seq, contig_coverage, args))
        curr_description += " %s_masked_at_%s_cov" % (dv.PROG_NAME, args.cov_depth)
        SeqIO.write(SeqRecord(curr_seq, description=curr_description, id=tig_id), aln_file, 'fasta')
        if args.verbose:
            print("%s masked!\n" % tig_id)

    # globs all the fasta files together into a master masked file
    with open("%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower()), 'w') as final_file:
        for out_file in blast_outfiles:
            with open(out_file, 'r') as of:
                copyfileobj(of, final_file)

    if args.verbose:
        print("Masked genome written to \'%s_%smasked.fasta\'" % (Path(args.genome).stem, dv.PROG_NAME.lower()))


# generates a masked chromosome sequence from the coverage list
def list_get_masked_genome_seq(genome_seq, cov_list, args):
    masked_seq = ""
    for base in range(0, len(genome_seq)):
        # if the base in the coverage string has >= coverage depth specified by user [default: 2]
        if cov_list[base] >= args.cov_depth:
            # masks that base in the provided genome using character specified by user [default: n]
            masked_seq += args.mask_char
        else:
            masked_seq += genome_seq[base]
    return masked_seq


# splits genome file into individual contigs and dumps them to the dump_path
def split_contigs(dump_path, genome_file):
    genome_index = SeqIO.index(str(genome_file), 'fasta')
    for record in genome_index:
        SeqIO.write(genome_index[record], "%s/%s.fasta" % (dump_path, genome_index[record].id), 'fasta')
    genome_index.close()
