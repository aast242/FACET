#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from shutil import rmtree

from . import blast_utils
from .defaults import ProgDefaults as dv


def masker_list_driver(args):
    # creates a temporary BLAST database #
    if Path("%s_masker_tempdb" % dv.PROG_NAME.lower()).exists():
        print("FATAL: \'%s_masker_tempdb\' already exists! Delete it if you want to run the masker."
              % dv.PROG_NAME.lower())
        exit()
    if args.verbose:
        print("\nCreating a temporary BLAST database....")

    dbfilepath = blast_utils.make_blast_db("%s_masker_tempdb" % dv.PROG_NAME.lower(), args.genome, ".", args)

    blast_aln = blast_utils.masker_blast(dbfilepath, args.genome, args)

    genome_dict = SeqIO.to_dict(SeqIO.parse(Path(args.genome).resolve(), 'fasta'))

    if args.verbose:
        print()

    for contig in genome_dict.keys():
        tig_len = len(genome_dict[contig])
        contig_coverage = [0] * tig_len  # generates a list of zeroes as long as the tig

        if args.verbose:
            print("Adding BLAST alignment ranges to coverage list for %s...." % contig)

        blast_aln.setdefault(contig, [])  # prevents errors if no blast hits on contig

        # adds ranges from combined hits to the tic_coverage list #
        # creates an operational list that needs to be converted  #
        for hit_records in blast_aln[contig]:
            contig_coverage[hit_records[dv.SSTART_INDEX] - 1] += 1
            if hit_records[dv.SEND_INDEX] < tig_len:  # CHANGED OPERATOR FROM != TO <
                contig_coverage[hit_records[dv.SEND_INDEX]] -= 1

        if args.verbose:
            print("Converting operational list to coverage list....")
        # converts operational list to actual coverage list #
        for basepair in range(1, tig_len):
            contig_coverage[basepair] += contig_coverage[basepair - 1]
        # masks each contig in the genome
        if args.verbose:
            print("Masking %s...." % contig)
        genome_dict[contig].seq = Seq(list_get_masked_genome_seq(genome_dict[contig].seq, contig_coverage, args))
        if args.verbose:
            print("%s masked!\n" % contig)

    if args.verbose:
        print("Temporary BLAST database removed!")
    rmtree(Path(dbfilepath).parent.absolute())

    if args.verbose:
        print("Masked genome file written to \'%s\'!" %
              Path("%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower())).absolute()
              )
    SeqIO.write(genome_dict.values(), "%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower()), "fasta")


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


# -------------------------------------- #
# ALL THE FUNCTIONS BELOW ARE DEPRECATED #
# -------------------------------------- #
# TODO: use sed to change the character 'sed s/./N/5' (5 is the same as str[4])

def masker_str_driver(args):
    # creates a temporary BLAST database #
    if args.verbose:
        print("\nCreating a temporary BLAST database....")

    dbfilepath = blast_utils.make_blast_db("%s_masker_tempdb" % dv.PROG_NAME.lower(), args.genome, ".", args)

    blast_aln = blast_utils.masker_blast(dbfilepath, args.genome, args)

    genome_dict = SeqIO.to_dict(SeqIO.parse(Path(args.genome).resolve(), 'fasta'))

    if args.verbose:
        print("Adding BLAST alignment ranges to coverage string....")
    for contig in genome_dict.keys():
        tig_len = len(genome_dict[contig])
        cov_string = str(dv.MASKER_MID_CHAR) * tig_len

        # adds ranges from combined hits to the tic_coverage list #
        # creates an operational list that needs to be converted  #
        for hit_records in blast_aln[contig]:
            # if the character at the start doesn't equal the ceiling character
            if cov_string[hit_records[dv.SSTART_INDEX] - 1] != dv.MASKER_CEILING_CHAR:
                # increment the character at the start by one
                cov_string = cov_string[:hit_records[dv.SSTART_INDEX] - 1] + \
                             increment_ascii(cov_string[hit_records[dv.SSTART_INDEX] - 1]) + \
                             cov_string[hit_records[dv.SSTART_INDEX]:]
                # if the alignment doesn't stop at the end of the contig (stops don't need to be added for the end)
                if hit_records[dv.SEND_INDEX] < tig_len:
                    # AND the final character is not equal to the floor
                    if cov_string[hit_records[dv.SEND_INDEX]] != dv.MASKER_FLOOR_CHAR:
                        # decrement the final character
                        cov_string = cov_string[:hit_records[dv.SEND_INDEX]] + \
                                     decrement_ascii(cov_string[hit_records[dv.SEND_INDEX]]) + \
                                     cov_string[hit_records[dv.SEND_INDEX] + 1:]
                    # if the final character was equal to the floor, set the starting character back to what it was
                    # before we incremented it bc the alignment doesn't need to be added
                    else:
                        cov_string = cov_string[:hit_records[dv.SSTART_INDEX] - 1] + \
                                     decrement_ascii(cov_string[hit_records[dv.SSTART_INDEX] - 1]) + \
                                     cov_string[hit_records[dv.SSTART_INDEX]:]
        # masks each contig in the genome
        if args.verbose:
            print("Masking %s...." % contig)
        genome_dict[contig].seq = get_masked_genome_seq(
            genome_dict[contig].seq, convert_operational_str(cov_string), args)
    rmtree(Path(dbfilepath).parent.absolute())
    SeqIO.write(genome_dict, "%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower()), "fasta")


def get_masked_genome_seq(genome_seq, cov_str, args):
    if genome_seq != cov_str:
        print("FATAL: contig length and coverage string length do not match!")
        exit()
    if args.verbose:
        print("Masking any base in \'%s\' with >= %s coverage depth with \'%s\'...." %
              (Path(args.genome).resolve(), args.cov_depth, args.mask_char))
    for base in range(0, len(genome_seq)):
        # if the base in the coverage string has >= coverage depth specified by user [default: 2]
        if ascii_to_int(cov_str[base]) - dv.MASKER_MID_INT >= args.cov_depth:
            # masks that base in the provided genome using character specified by user [default: n]
            genome_seq = genome_seq[:base] + args.mask_char + genome_seq[base + 1:]

    return genome_seq


def convert_operational_str(cov_string):
    ceiling_debt = 0
    ascii_offset = dv.MASKER_MID_INT
    for cov_index in range(1, len(cov_string)):
        curr_ascii = ascii_to_int(cov_string[cov_index]) - ascii_offset
        prev_ascii = ascii_to_int(cov_string[cov_index - 1]) - ascii_offset

        # If you're subtracting from the current ascii and there is unpaid debt
        if curr_ascii < 0 and ceiling_debt > 0:
            # See if the debt has been paid
            payment = ceiling_debt + curr_ascii
            # If the subtracting didn't pay off the debt
            if payment > 0:
                # new ascii is still the masker ceiling
                new_ascii = dv.MASKER_CEILING_INT - ascii_offset
            # If subtracting DID pay off the debt
            else:
                # new_ascii is the positive amount of remainder
                new_ascii = -payment
        # any other situation (no debt/subtracting or some debt but adding)
        else:
            # calculate the new ascii
            new_ascii = curr_ascii + prev_ascii
            # if the new ascii breaks the ceiling
            if new_ascii > dv.MASKER_CEILING_INT - ascii_offset:
                # add anything over 45 to the ceiling debt
                ceiling_debt += new_ascii - (dv.MASKER_CEILING_INT - ascii_offset)
                # set new ascii equal to the ceiling (45 in our case)
                new_ascii = dv.MASKER_CEILING_INT - ascii_offset

        # Set the coverage character to the newly calculated value
        cov_string = cov_string[:cov_index] + int_to_ascii(new_ascii+ascii_offset) + cov_string[cov_index + 1:]
    return cov_string


def increment_ascii(ascii_char):
    return int_to_ascii(ascii_to_int(ascii_char) + 1)


def decrement_ascii(ascii_char):
    return int_to_ascii(ascii_to_int(ascii_char) - 1)


def ascii_to_int(ascii_str):
    return ord(ascii_str)


def int_to_ascii(user_int):
    return str(chr(user_int))
