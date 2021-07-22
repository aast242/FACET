#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import csv
import os
from sys import exc_info
from pathlib import Path
from shutil import copyfileobj
from pprint import pprint

from Bio import SeqIO

from .defaults import ProgDefaults as dv
from . import btop_utils


# driver function for the variant call format module
def vc_driver(args):
    final_vc_outdir = "%s/%s_%s_VCF" % (dv.OUTDIR_NAME, Path(args.genome).stem.replace(".", "_"),
                                        Path(args.vcf_file).stem.replace(".", "_"))

    true_snps_fp = "%s/%s_%s_%s" % (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                    Path(args.vcf_file).stem.replace(".", "_"), dv.VCF_TRUE_SNP_FILE)
    false_snps_fp = "%s/%s_%s_%s"% (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                     Path(args.vcf_file).stem.replace(".", "_"), dv.VCF_FALSE_SNP_FILE)
    summary_file_fp = "%s/%s_%s_%s" % (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                       Path(args.vcf_file).stem.replace(".", "_"), dv.VCF_SUMM_FILE)
    final_rplot_fp = "%s/%s_%s_%s" % (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                      Path(args.vcf_file).stem.replace(".", "_"), dv.RPLOT_FINAL_OUT)
    # get the vc driver to where the tempdir is filled with BLAST outfiles and get an os.listdir at that point
    # so that you have the outfiles to iterate through saved and can safely output the vcf_header and variants.vcf files
    genome_indx = SeqIO.index(args.genome, 'fasta')

    # gets organised BLAST hits that have not been cleaned
    blast_outfiles = ["%s/%s" % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR)]

    # get VCF file that contains only SNPs
    if args.verbose:
        print("\nParsing VCF file to find bases with SNP calls....")
    flush_vcf_file(Path(args.vcf_file).resolve())

    # split SNP_only vcf up into contigs
    try:
        os.mkdir(dv.VCF_TEMPDIR)
    except Exception:
        print("FATAL: problem making VCF tempdir: %s" % exc_info()[0])
    if args.verbose:
        print("Flushing SNPs to files based on contig ID")

    flush_var_to_tigs(dv.VCF_VARIANT_TEMPFILE, dv.VCF_TEMPDIR, dv.VCF_FMT_DICT["chrom"])

    true_cnt = 0
    true_snp = []
    false_cnt = 0
    false_snp = []
    cycle_count = 0
    for of in blast_outfiles:
        curr_tigid = Path(of).suffix.replace(".", "")

        if args.verbose:
            print("\nDetermining which %s SNPs are called in repetitive regions...." % curr_tigid)
        curr_vcf = "%s/%s.vcf" % (dv.VCF_TEMPDIR, curr_tigid)
        if not Path(curr_vcf).resolve().exists():
            if args.verbose:
                print("\'%s\' does not exist! Skipping..." % curr_vcf)
            continue
        curr_lst = gen_cov_lst(of, genome_indx, args)
        # ^ might seem like a lot to load cov_list into memory, but the largest human chromosome only took ~2gb of mem

        # IF YOU ARE GENERATING AN RPLOT
        if args.rplot:
            write_rplot_file(curr_vcf, curr_lst, curr_tigid, args)

        with open(curr_vcf, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                cycle_count += 1
                line = line.rstrip().split("\t")
                line[dv.VCF_FMT_DICT["pos"]] = int(line[dv.VCF_FMT_DICT["pos"]])

                # if the region is considered repetitive, add the variant to the false_snp list
                if curr_lst[line[dv.VCF_FMT_DICT["pos"]] - 1] >= args.cov_depth:
                    false_snp.append(line)
                    false_cnt += 1
                # if the region is not considered repetitive, add variant to true_snp list
                else:
                    true_snp.append(line)
                    true_cnt += 1

                # if the cycle count is bigger than the flush length, flush variants to files
                if cycle_count >= dv.FLUSH_LEN:
                    if false_snp:
                        with open("%s/false_snp.vcf" % dv.VCF_TEMPDIR, 'a', newline='') as csvfile:
                            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                    lineterminator="\n")
                            for records in false_snp:
                                filewriter.writerow(records)
                        false_snp = []
                    if true_snp:
                        with open("%s/true_snp.vcf" % dv.VCF_TEMPDIR, 'a', newline='') as csvfile:
                            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                    lineterminator="\n")
                            for records in true_snp:
                                filewriter.writerow(records)
                        true_snp = []
                    cycle_count = 0

        # clear curr_list each iteration to preserve memory
        del curr_lst

    # dump alns at the end if cycle count is greater than zero
    if cycle_count > 0:
        if false_snp:
            with open("%s/false_snp.vcf" % dv.VCF_TEMPDIR, 'a', newline='') as csvfile:
                filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                        lineterminator="\n")
                for records in false_snp:
                    filewriter.writerow(records)
            del false_snp
        if true_snp:
            with open("%s/true_snp.vcf" % dv.VCF_TEMPDIR, 'a', newline='') as csvfile:
                filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                        lineterminator="\n")
                for records in true_snp:
                    filewriter.writerow(records)
            del true_snp

    total_cnt = false_cnt + true_cnt

    try:
        os.mkdir(final_vc_outdir)
    except FileExistsError:
        pass
    except Exception:
        print("Problem creating final VCF output directory! Exiting....")
        raise

    # write summary file
    with open(summary_file_fp, 'w') as f:
        if args.verbose:
            print("\n~~~~~~~~~~~~~~~~~~")
            print("# of total SNPs: %s" % total_cnt)
            print("# of SNPs in repetitive regions (coverage >= %s): %s (%s/%s, %s%%)" %
                  (args.cov_depth, false_cnt, false_cnt, total_cnt, "{:.2f}".format(float(false_cnt / total_cnt) * 100))
                  )
            print("# of SNPs in non-repetitive regions (coverage < %s): %s (%s/%s, %s%%)" %
                  (args.cov_depth, true_cnt, true_cnt, total_cnt, "{:.2f}".format(float(true_cnt / total_cnt) * 100))
                  )
            print("~~~~~~~~~~~~~~~~~~")
        f.write("Genome file: %s" % Path(args.genome).resolve())
        f.write('\n')
        f.write("Initial VCF file: %s" % Path(args.vcf_file).resolve())
        f.write('\n')
        f.write("# of total SNPs: %s" % total_cnt)
        f.write('\n')
        f.write("# of SNPs in repetitive regions (coverage >= %s): %s (%s/%s, %s%%)" %
                (args.cov_depth, false_cnt, false_cnt, total_cnt, "{:.2f}".format(float(false_cnt / total_cnt) * 100))
                )
        f.write('\n')
        f.write("# of SNPs in non-repetitive regions (coverage < %s): %s (%s/%s, %s%%)" %
                (args.cov_depth, true_cnt, true_cnt, total_cnt, "{:.2f}".format(float(true_cnt / total_cnt) * 100))
                )
        f.write('\n')

    # glob header and false SNPs together into master_false
    if args.verbose:
        print("\nWriting SNPs in repetitive regions to \'%s\'...." % false_snps_fp)
    with open(false_snps_fp, 'wb') as tmp:
        with open(dv.VCF_HEADER_TEMPFILE, 'rb') as of:
            copyfileobj(of, tmp)
        with open("%s/false_snp.vcf" % dv.VCF_TEMPDIR, 'rb') as of:
            copyfileobj(of, tmp)

    # glob header and true SNPs together into master_true
    if args.verbose:
        print("Writing SNPs in non-repetitive regions to \'%s\'...." % true_snps_fp)
    with open(true_snps_fp, 'wb') as tmp:
        with open(dv.VCF_HEADER_TEMPFILE, 'rb') as of:
            copyfileobj(of, tmp)
        with open("%s/true_snp.vcf" % dv.VCF_TEMPDIR, 'rb') as of:
            copyfileobj(of, tmp)

    # copy over RPLOT file
    if args.rplot:
        if args.verbose:
            print("Copying temp rplot file to \'%s\'...." % final_rplot_fp)
        with open(final_rplot_fp, 'wb') as tmp:
            with open(dv.RPLOT_TEMP_FILE, 'rb') as rpt:
                copyfileobj(rpt, tmp)

# parses a vcf file and dumps 2 files to tempdir: vcf_header, variants.vcf
def flush_vcf_file(vcf_filepath):
    header = []
    w_alt = []

    cycle_count = 0
    # read file line by line
    with open(vcf_filepath, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break

            # if there is a leading pound, the line is a header
            if line[0] == "#":
                cycle_count += 1
                line = line.rstrip().split("\t")
                header.append(line)

            # no leading pound means the line contains variant information
            else:
                line = line.rstrip().split("\t")

                # if there are variants called, add it to the with alt list
                if line[dv.VCF_FMT_DICT["alt"]] != ".":
                    cycle_count += 1
                    w_alt.append(line)

            # if cycle count is larger than the flush length, flush out to files
            if cycle_count >= dv.FLUSH_LEN:
                if header:  # if the header is not empty, write to file
                    with open(dv.VCF_HEADER_TEMPFILE, 'a', newline='') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter='\t', escapechar='', quotechar='',
                                                quoting=csv.QUOTE_NONE, lineterminator="\n")
                        for records in header:
                                filewriter.writerow(records)

                    header = []

                if w_alt:  # if with alt list is not empty, flush to file
                    with open(dv.VCF_VARIANT_TEMPFILE, 'a', newline='') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                lineterminator="\n")
                        for records in w_alt:
                            filewriter.writerow(records)
                    w_alt = []
                cycle_count = 0

    if cycle_count > 0:
        if header:  # if the header is not empty, flush to file
            with open(dv.VCF_HEADER_TEMPFILE, 'a', newline='') as csvfile:
                filewriter = csv.writer(csvfile, delimiter='\t',escapechar='', quotechar='', quoting=csv.QUOTE_NONE,
                                        lineterminator="\n")
                for records in header:
                    filewriter.writerow(records)
            del header
        if w_alt:  # if with alt list is not empty, flush to file
            with open(dv.VCF_VARIANT_TEMPFILE, 'a', newline='') as csvfile:
                filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                        lineterminator="\n")
                for records in w_alt:
                    filewriter.writerow(records)
            del w_alt


# flush variant file to tig files
def flush_var_to_tigs(var_file, flushdir, flushindx):
    tig_dict = {}
    cycle_count = 0

    with open(var_file, 'r') as f:  # iterate through the file
        while True:
            # read each line
            line = f.readline()
            # if the line is the last line of the file, exit the loop
            if not line:
                break
            cycle_count += 1

            # parse the line
            line = line.rstrip().split("\t")

            # if the tig id is not in the dictionary, add it
            tig_dict.setdefault(line[flushindx], [])

            # append the line to the tig id list in the dictionary
            tig_dict[line[flushindx]].append(line)

            # flush each contig to a file if there have been enough cycles
            if cycle_count >= dv.FLUSH_LEN:
                for tig_id in tig_dict.keys():
                    if tig_dict[tig_id]:
                        temp_filename = "%s/%s.vcf" % (flushdir, tig_id)
                        with open(temp_filename, 'a') as csvfile:
                            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                    lineterminator="\n")
                            for records in tig_dict[tig_id]:
                                filewriter.writerow(records)
                            # sets the tig list to empty after it's flushed
                        tig_dict[tig_id] = []
                # resets cycle count
                cycle_count = 0

    # if there are still hits in the tig dictionary
    if cycle_count > 0:
        for tig_id in tig_dict.keys():
            if tig_dict[tig_id]:
                temp_filename = "%s/%s.vcf" % (flushdir, tig_id)
                with open(temp_filename, 'a') as csvfile:
                    filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                            lineterminator="\n")
                    for records in tig_dict[tig_id]:
                        filewriter.writerow(records)
                    # sets the tig list to empty after it's flushed
                tig_dict[tig_id] = []


# generates a list as long as a contig that contains coverage information at each base
def gen_cov_lst(selfblast_outfile, genome_indx, args):
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]

    chr_name = Path(selfblast_outfile).suffix.replace(".", "")
    chr_len = len(genome_indx[chr_name])
    chr_lst = [0] * chr_len
    with open(Path(selfblast_outfile).resolve(), 'r') as curr_file:
        while True:
            line = curr_file.readline()
            if not line:
                break

            line = line.rstrip().split('\t')

            # # removes FACET IDs from the beginning of the line
            # if "%s_" % dv.PROG_NAME in line[0]:
            #     del line[0]

            line[sstart] = dv.BLASTN_TYPES["sstart"](line[sstart])
            line[send] = dv.BLASTN_TYPES["send"](line[send])

            # switches minus hits so they're [min,max]
            if line[sstart] > line[send]:
                line[sstart], line[send] = min(line[send], line[sstart]), max(line[send], line[sstart])

            # adds ranges from combined hits to the tic_coverage list #
            chr_lst[line[sstart] - 1] += 1
            if line[send] < chr_len:  # CHANGED OPERATOR FROM != TO <
                chr_lst[line[send]] -= 1

    # converts operational list to actual coverage list #
    for basepair in range(1, chr_len):
        chr_lst[basepair] += chr_lst[basepair - 1]

    return chr_lst


# gets fraction_RIP, VCF_coverage, and BLAST_coverage for a window
def get_window_stats(window_start_indx, vcf_list, blast_cov_list, args):
    num_ripd = 0
    non_transition_muts = 0
    avg_VCF_cov = 0
    avg_BLAST_cov = 0

    # iterates through the vcf list starting at window_start_index and going until it has iterated window_size times
    for i in range(0, args.rp_window_size):
        curr_indx = window_start_indx + i

        # checks ref base against alt base and sees if RIP is possible
        if btop_utils.check_if_ripd(vcf_list[curr_indx][dv.VCF_FMT_DICT["ref"]],
                                    vcf_list[curr_indx][dv.VCF_FMT_DICT["alt"]]):
            num_ripd += 1
        else:
            non_transition_muts += 1

        # gets the first sample value (only one we care about for our purposes)
        search_term_indx = 0
        search_list = vcf_list[curr_indx][dv.VCF_FMT_DICT["format"]].split(":")

        for k in range(0, len(search_list)):
            if search_list[k] == dv.RPLOT_SEARCH_TERM:
                search_term_indx = k
                break

        avg_VCF_cov += int(vcf_list[curr_indx][len(dv.VCF_FMT_DICT)].split(":")[search_term_indx].split(",")[1])

        # gets blast coverage for indicated base
        avg_BLAST_cov += blast_cov_list[int(vcf_list[curr_indx][dv.VCF_FMT_DICT["pos"]]) - 1]

    if non_transition_muts == 0:
        non_transition_muts = 1

    # computes values for each stat
    num_ripd = float(num_ripd / non_transition_muts)
    avg_VCF_cov = float(avg_VCF_cov / args.rp_window_size)
    avg_BLAST_cov = float(avg_BLAST_cov / args.rp_window_size)

    return num_ripd, avg_VCF_cov, avg_BLAST_cov


# takes variants.vcf and flushes out entries that have <= cutoff in AD[0]
def flush_filtered_vcf(varfile):
    cycle_count = 0
    filtered_vcf = []
    outfile_path = "%s/%s_filtered.vcf" % (dv.VCF_TEMPDIR, Path(varfile).stem.replace(".", "_"))
    with open(varfile, 'r') as vtf:
        while True:
            # read each line
            line = vtf.readline()
            # if the line is the last line of the file, exit the loop
            if not line:
                break

            cycle_count += 1
            line = line.rstrip().split("\t")

            fmt_lst = line[dv.VCF_FMT_DICT["format"]].split(":")

            search_indx = 0

            # finds where the search term is in the format string
            for i in range(0, len(fmt_lst)):
                if fmt_lst[i] == dv.RPLOT_SEARCH_TERM:
                    search_indx = i
                    break

            # if seachstring value is <= cutoff, add to filter list
            if int(line[len(dv.VCF_FMT_DICT)].split(":")[search_indx].split(",")[0]) <= dv.RPLOT_REF_AGREE_CTOFF:
                filtered_vcf.append(line)

            if cycle_count >= dv.FLUSH_LEN:
                with open(outfile_path, 'a') as csvfile:
                    filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                            lineterminator="\n")
                    for records in filtered_vcf:
                        filewriter.writerow(records)
                cycle_count = 0
                filtered_vcf = []
    if cycle_count > 0:
        with open(outfile_path, 'a') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            for records in filtered_vcf:
                filewriter.writerow(records)
    return outfile_path


def write_rplot_file(varfile, blast_cov_list, chrid, args):
    # adds header if the file doesn't exist #
    if not Path(dv.RPLOT_TEMP_FILE).exists():
        with open(dv.RPLOT_TEMP_FILE, 'a') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            for records in [dv.RPLOT_HEADERS]:
                filewriter.writerow(records)

    filtered_vcf = flush_filtered_vcf(varfile)

    cycle_count = 0
    with open(filtered_vcf, 'r') as f:

        # initialize first window
        current_window = [next(f) for i in range(args.rp_window_size)]
        for i in range(0, len(current_window)):
            current_window[i] = current_window[i].rstrip().split('\t')
        rip_val, vcf_cov, blast_cov = get_window_stats(0, current_window, blast_cov_list, args)
        rplot_stats = [[chrid, current_window[0][dv.VCF_FMT_DICT["pos"]],
                        round(rip_val, 2), vcf_cov, blast_cov]]

        while True:
            line = f.readline()

            if not line:
                break

            cycle_count += 1

            line = line.rstrip().split('\t')

            # add new line and remove first line
            current_window.append(line)
            del current_window[0]

            # if there has been an appropriate step, add in new stats
            if cycle_count % args.rp_window_step == 0:
                rip_val, vcf_cov, blast_cov = get_window_stats(0, current_window, blast_cov_list, args)
                rplot_stats.append([chrid, current_window[0][dv.VCF_FMT_DICT["pos"]],
                                    round(rip_val, 2), vcf_cov, blast_cov])

                if cycle_count >= dv.FLUSH_LEN:
                    with open(dv.RPLOT_TEMP_FILE, 'a') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                lineterminator="\n")
                        for records in rplot_stats:
                            filewriter.writerow(records)
                    rplot_stats = []
                    cycle_count = 0

    if cycle_count > 0:
        with open(dv.RPLOT_TEMP_FILE, 'a') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            for records in rplot_stats:
                filewriter.writerow(records)


def check_vc_overwrite(args):
    final_vc_outdir = "%s/%s_%s_VCF" % (dv.OUTDIR_NAME, Path(args.genome).stem.replace(".", "_"),
                                        Path(args.vcf_file).stem.replace(".", "_"))

    true_snps = "%s/%s_%s_%s" % (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                 Path(args.vcf_file).stem.replace(".", "_"), dv.VCF_TRUE_SNP_FILE)
    false_snps = "%s/%s_%s_%s"% (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                  Path(args.vcf_file).stem.replace(".", "_"), dv.VCF_FALSE_SNP_FILE)
    summary_file = "%s/%s_%s_%s" % (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                    Path(args.vcf_file).stem.replace(".", "_"), dv.VCF_SUMM_FILE)
    if args.rplot:

        final_rplot = "%s/%s_%s_%s" %  (final_vc_outdir, Path(args.genome).stem.replace(".", "_"),
                                        Path(args.vcf_file).stem.replace(".", "_"), dv.RPLOT_FINAL_OUT)

        if Path(final_rplot).exists() and args.force is False:
            print("FATAL: \'%s\' already exists! use --force to proceed" % final_rplot)
            exit()
        elif Path(final_rplot).exists() and args.force is True:
            os.remove(final_rplot)

    if Path(true_snps).exists() and args.force is False:
        print("FATAL: \'%s\' already exists! use --force to proceed" % true_snps)
        exit()
    elif Path(true_snps).exists() and args.force is True:
        os.remove(true_snps)

    if Path(false_snps).exists() and args.force is False:
        print("FATAL: \'%s\' already exists! use --force to proceed" % false_snps)
        exit()
    elif Path(false_snps).exists() and args.force is True:
        os.remove(false_snps)

    if Path(summary_file).exists() and args.force is False:
        print("FATAL: \'%s\' already exists! use --force to proceed" % summary_file)
        exit()
    elif Path(summary_file).exists() and args.force is True:
        os.remove(summary_file)