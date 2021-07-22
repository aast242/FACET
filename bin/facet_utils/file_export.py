#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import subprocess
import os
import csv
from pathlib import Path
from mmap import mmap
from shutil import copyfileobj

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import btop_utils
from .defaults import ProgDefaults as dv


# checks each of the files that would be created and sees if the file already exists           #
# If it does, refuses to overwrite unless --force is specified                                 #
# If force is specified, deletes files that would be overwritten so that flushes work properly #
def check_file_overwrite(datafilepath, filepath, args):
    ok_check = True
    if args.writegff:
        if Path(datafilepath + "/" + filepath + ".gff").exists():
            ok_check = False
            if args.force:
                os.remove(Path(datafilepath + "/" + filepath + ".gff"))
    if args.writesam:
        if args.verboseoutput:
            if Path(datafilepath + "/" + filepath + "_verbose.sam").exists():
                ok_check = False
                if args.force:
                    os.remove(Path(datafilepath + "/" + filepath + "_verbose.sam"))
        else:
            if Path(datafilepath + "/" + filepath + ".sam").exists():
                ok_check = False
                if args.force:
                    os.remove(Path(datafilepath + "/" + filepath + ".sam"))
    if args.writebam:
        if args.verboseoutput:
            if Path(datafilepath + "/" + filepath + "_verbose.bam").exists() or \
                    Path(datafilepath + "/" + filepath + "_verbose.bai").exists():
                ok_check = False
                if args.force:
                    os.remove(Path(datafilepath + "/" + filepath + "_verbose.bam"))
                    os.remove(Path(datafilepath + "/" + filepath + "_verbose.bai"))
        else:
            if Path(datafilepath + "/" + filepath + ".bam").exists() or \
                    Path(datafilepath + "/" + filepath + ".bai").exists():
                ok_check = False
                if args.force:
                    os.remove(Path(datafilepath + "/" + filepath + ".bam"))
                    os.remove(Path(datafilepath + "/" + filepath + ".bai"))
    if args.nocsv:
        if args.verboseoutput:
            if Path(datafilepath + "/" + filepath + "_verbose.csv").exists():
                ok_check = False
                if args.force:
                    os.remove(Path(datafilepath + "/" + filepath + "_verbose.csv"))
        else:
            if Path(datafilepath + "/" + filepath + ".csv").exists():
                ok_check = False
                if args.force:
                    os.remove(Path(datafilepath + "/" + filepath + ".csv"))
    if args.writetic:
        if Path(datafilepath + "/" + filepath + ".tic").exists():
            ok_check = False
            if args.force:
                os.remove(Path(datafilepath + "/" + filepath + ".tic"))
    if args.writefasta > 0:
        if Path(datafilepath + "/" + filepath + "_" + str(args.writefasta) + ".fasta").exists():
            ok_check = False
            if args.force:
                os.remove(Path(datafilepath + "/" + filepath + "_" + str(args.writefasta) + ".fasta"))
    # exit()
    return ok_check


# Writes a gff file containing the locations of query elements in the reference #
def write_gff_new(aln_files, filepath, data_outdir, queryfile, args):
    sseqid = args.outfmt["sseqid"]
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]
    qseqid = args.outfmt["qseqid"]
    qstart = args.outfmt["qstart"]
    qend = args.outfmt["qend"]
    sstrand = args.outfmt["sstrand"]

    if args.verbose:
        print("Writing hits to \'" + data_outdir + "/" + filepath + ".gff\'")

    # gets lengths of query sequences and stores them in a dictionary #
    querydict = SeqIO.to_dict(SeqIO.parse(queryfile, 'fasta'))
    for i in querydict.keys():
        querydict[i] = len(querydict[i].seq)

    gfflist = [["##gff-version 3"]]
    cycle_count = 0
    outfmt_specs = list(args.outfmt.keys())
    uid_num = 0
    # goes through each tig in the temp directory
    for tig in aln_files:
        with open(tig, "r") as f:
            while True:
                # reads line-by-line to reduce memory load
                line = f.readline()
                if not line:
                    break

                # increments cycle count
                cycle_count += 1
                uid_num += 1

                line = line.rstrip().split('\t')

                # casts elements to the correct type
                for outfmt_feat in outfmt_specs:
                    line[args.outfmt[outfmt_feat]] = \
                        dv.BLASTN_TYPES[outfmt_feat](line[args.outfmt[outfmt_feat]])

                # gets length of the query element
                ele_len = querydict[line[qseqid]]

                # if hit contains both the beginning and the end, make it cyan
                if line[qend] >= ele_len - args.buffer and line[qstart] <= args.buffer:
                    color_str = dv.BOTH_ENDS_COLOR
                # if hit contains only the end, make it yellow
                elif line[qend] >= ele_len - args.buffer:
                    color_str = dv.ONLY_END_COLOR
                # if hit contains only the beginning, make it orange
                elif line[qstart] <= args.buffer:
                    color_str = dv.ONLY_BEGINNING_COLOR
                # if hit contains neither, make it grey
                else:
                    color_str = dv.NEITHER_END_COLOR

                if line[sstrand] == 'plus':
                    strand = '+'
                else:
                    strand = '-'

                gfflist.append([line[sseqid], "%s" % dv.PROG_NAME, line[qseqid], line[sstart], line[send],
                                '.', strand, '.', 'color=%s; ID=%s_%s;Name=%s;Query=start: %s, end: %s;'
                                                  'Subject=start: %s, end: %s;Length=%s'
                                % (color_str, dv.PROG_NAME, uid_num, line[qseqid], line[qstart], line[qend],
                                   line[sstart], line[send], line[send] - line[sstart])])
                # if the number of items in the list reaches the flush length, flush everything out to the GFF file
                if cycle_count >= dv.FLUSH_LEN:
                    with open(data_outdir + "/" + filepath + ".gff", 'a', newline='') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                lineterminator="\n")
                        for records in gfflist:
                            filewriter.writerow(records)
                    gfflist = []
                    cycle_count = 0

    # writes any remaining hits in the GFF list that haven't been flushed
    if cycle_count > 0:
        with open(data_outdir + "/" + filepath + ".gff", 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            for hits in gfflist:
                filewriter.writerow(hits)
        gfflist = None
        del gfflist


# Writes a sam/bam file containing the locations of query elements in the reference #
# This file can include sequence/alignment info if -q is specified                  #
def write_sambam_new(contigs, aln_files, data_outdir, filepath, args):

    sseqid = args.outfmt["sseqid"]
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]
    qseqid = args.outfmt["qseqid"]
    qstart = args.outfmt["qstart"]
    qend = args.outfmt["qend"]
    sstrand = args.outfmt["sstrand"]
    if "btop" in args.outfmt:
        btop = args.outfmt["btop"]

    if args.writebam and not args.writesam:
        samstring = data_outdir + "/" + "%s_temp_verbose.sam" % dv.PROG_NAME
    else:
        if args.verboseoutput:
            samstring = data_outdir + "/" + filepath + "_verbose.sam"
        else:
            samstring = data_outdir + "/" + filepath + ".sam"

    samlist = [["@HD", "VN:1.6"]]
    # TODO: the SeqIO.parse here might not be great for longer sequences
    for record in SeqIO.parse(contigs, "fasta"):
        samlist.append(["@SQ", "SN:" + record.id, "LN:" + str(len(record.seq))])
    if args.verboseoutput:
        query_dict = SeqIO.to_dict(SeqIO.parse(args.query, 'fasta'))
        if args.verbose:
            print("Parsing BTOP sequences and converting to CIGAR sequences....")
    if not args.writebam or args.writesam:
        if args.verbose:
            print("Writing hits to \'" + samstring + "\'")

    cycle_count = 0
    uid_count = 0
    outfmt_specs = list(args.outfmt.keys())
    # goes through each tig in the temp directory
    for tig in aln_files:
        with open(tig, "r") as f:
            while True:
                # reads line-by-line to reduce memory load
                line = f.readline()
                if not line:
                    break

                # increments cycle count
                cycle_count += 1
                uid_count += 1

                line = line.rstrip().split('\t')

                # casts elements to the correct type
                for outfmt_feat in outfmt_specs:
                    line[args.outfmt[outfmt_feat]] = \
                        dv.BLASTN_TYPES[outfmt_feat](line[args.outfmt[outfmt_feat]])

                if args.verboseoutput:
                    qseq = str(query_dict[line[qseqid]].seq)[line[qstart] - 1:line[qend]]
                if line[sstrand] == 'plus':
                    flag = 0
                    if args.verboseoutput:
                        cigar_seq = btop_utils.btop_to_cigar(btop_utils.parse_btop(line[btop]))
                else:
                    flag = 16
                    if args.verboseoutput:
                        # BTOP sequence needs to be reverse complimented, +/-
                        cigar_seq = btop_utils.btop_to_cigar(btop_utils.parse_btop(
                            btop_utils.reverse_compliment(line[btop])))
                        qseq = str(Seq(qseq).reverse_complement())
                if args.verboseoutput:
                    qseq = qseq.replace("-", "")
                    samlist.append(["%s_%s" % (dv.PROG_NAME, uid_count), flag, line[sseqid], line[sstart],
                                    255, cigar_seq, "*", 0, 0, qseq, "*",
                                    "CO:Z:query(" + line[qseqid] + "), qloc(" + str(line[qstart]) +
                                    "," + str(line[qend]) + ")," + " sloc(" + str(line[sstart]) + ","
                                    + str(line[send]) + "), len(" + str(line[send] - line[sstart] + 1) + ")"])
                else:
                    samlist.append(["%s_%s" % (dv.PROG_NAME, uid_count), flag, line[sseqid], line[sstart], 255,
                                    str(line[send] - line[sstart] + 1) + "M", "*", 0, 0, "*", "*",
                                    "CO:Z:query(" + line[qseqid] + "), qloc(" + str(line[qstart]) +
                                    "," + str(line[qend]) + ")," + " sloc(" + str(line[sstart]) + ","
                                    + str(line[send]) + "), len(" + str(line[send] - line[sstart] + 1) + ")"])
                if cycle_count >= dv.FLUSH_LEN:
                    with open(samstring, 'a', newline='') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                lineterminator="\n")
                        for records in samlist:
                            filewriter.writerow(records)
                    cycle_count = 0
                    samlist = []

    # flush any remaining records into the sam file
    if len(samlist) > 0:
        with open(samstring, 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            for hits in samlist:
                filewriter.writerow(hits)

    # convert written sam file into a bam file and index it
    if args.writebam:
        if args.verboseoutput:
            file_end = "_verbose.bam"
        else:
            file_end = ".bam"
        if args.verbose:
            print("Writing hits to \'" + data_outdir + "/" + filepath
                  + file_end + "\'")
        subprocess.call(["samtools", "view", "-bS", samstring],
                        stdout=open(data_outdir + "/" + filepath + file_end, 'w'),
                        stderr=open(os.devnull, 'w'))
        if args.verbose:
            print("Indexing bam file....")
        subprocess.call(["samtools", "index", data_outdir + "/" + filepath + file_end,
                         data_outdir + "/" + filepath + file_end[:-1] + "i"])
        if not args.writesam:
            os.remove(samstring)
    samlist = None
    del samlist


# Writes a tab delimited csv file containing information on the location and % ident of query elements in the ref #
# TODO: could honestly just glob files together, right? No need to open and parse
def write_csv_new(aln_files, data_outdir, filepath, args):

    if args.verboseoutput:
        csvname = data_outdir + "/" + filepath + "_verbose.csv"
    else:
        csvname = data_outdir + "/" + filepath + ".csv"

    if args.verbose:
        print("Writing hits to \'" + csvname + "\'")

    hits = []
    cycle_count = 0
    # goes through each tig in the temp directory
    for tig in aln_files:
        with open(tig, "r") as f:
            while True:
                # reads line-by-line to reduce memory load
                line = f.readline()
                if not line:
                    break

                cycle_count += 1
                line = line.rstrip().split('\t')
                hits.append(line)

                # if the number of cycles is > the flush value, flush to file
                if cycle_count >= dv.FLUSH_LEN:
                    with open(csvname, 'a', newline='') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='|',
                                                quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
                        for alns in hits:
                            filewriter.writerow(alns)
                    cycle_count = 0
                    hits = []
    # if there are hits still in the stream, flush to file
    if cycle_count > 0:
        with open(csvname, 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='|',
                                    quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
            for alns in hits:
                filewriter.writerow(alns)


# Writes a TIg Coverage file that contains information about coverage depth on the contig.                   #
# Really only useful in Genome x Genome comparisons that aren't cleaned to try to find the location of novel #
# repeat elements and inform evolutionary history of existing insertions.                                    #
# This file format is in development and is currently slightly stable                                        #
def write_tic_new(contigs, aln_files, data_outdir, filepath, args):
    sstart = args.outfmt["sstart"]  # 2
    send = args.outfmt["send"] # 3

    # the following code produces a dictionary called tig_lengths such that {"tig01": tiglength} #
    tig_lengths = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta'))
    tempdict = {}
    for tignames in tig_lengths.keys():
        tempdict[tignames] = len(tig_lengths[tignames].seq)
    tig_lengths = tempdict
    tempdict = None
    del tempdict

    if args.verbose:
        print()

    cycle_count = 0
    final_tl = []
    # goes through each tig in the temp directory
    for tig in aln_files:
        with open(tig, "r") as f:

            # adds the tig name to the beginning of each tig
            tig_names = tig.split(".")
            del tig_names[0]
            tig_names = ".".join(tig_names)
            final_tl.append(">" + tig_names)
            cycle_count += 1

            # generates a list of zeroes as long as the tig
            tic_coverage = [0] * tig_lengths[tig_names]

            if args.verbose:
                print("Writing %s to tic file...." % str(tig_names))

            # reads line-by-line to reduce memory load
            while True:
                line = f.readline()
                if not line:
                    break

                line = line.rstrip().split('\t')

                # casts necessary parts to ints
                line[sstart] = int(line[sstart])
                line[send] = int(line[send])

                # adds ranges from combined hits to the tic_coverage list #
                # creates an operational list that needs to be converted  #
                tic_coverage[line[sstart] - 1] += 1
                if line[send] - 1 < tig_lengths[tig_names] - 1:  # CHANGED OPERATOR FROM != TO <
                    tic_coverage[line[send]] -= 1

        # converts operational list to actual coverage list #
        for basepair in range(1, tig_lengths[tig_names]):
            tic_coverage[basepair] += tic_coverage[basepair - 1]

        # tic_coverage = tic_coverage[5191538:5197197]

        # adds the beginning of every contig to the tic file #
        final_tl.append("1:%s" % (tic_coverage[0]))
        cycle_count += 1

        for bases in range(1, len(tic_coverage) - 1):
            cur_base = tic_coverage[bases]
            pre_base = tic_coverage[bases - 1]
            # if exiting a unique region, adds non-unique base to tic file #
            if pre_base == 1 and cur_base != 1:
                final_tl.append("%s:x" % (bases + 1))
            # if entering a unique region, adds starting unique base to tic file #
            elif cur_base == 1 and pre_base != 1:
                final_tl.append("%s:e" % (bases + 1))
            # if there is a large change in coverage, adds the base to the tic file #
            elif cur_base - pre_base >= dv.TIC_DELTA or cur_base - pre_base <= -dv.TIC_DELTA:
                final_tl.append("%s:%s" % (bases + 1, cur_base))
            cycle_count += 1
            if cycle_count >= dv.FLUSH_LEN:
                with open(data_outdir + "/" + filepath + ".tic", 'a', newline='') as ticfile:
                    for entry in final_tl:
                        ticfile.write(entry)
                        ticfile.write("\n")
                cycle_count = 0
                final_tl = []
        # adds the end of every contig to the tic tile #
        final_tl.append("%s:%s" % (len(tic_coverage), tic_coverage[len(tic_coverage) - 1]))

        # clears the tic_coverage variable every iteration to free up memory #
        tic_coverage = None

    if cycle_count > 0:
        with open(data_outdir + "/" + filepath + ".tic", 'a', newline='') as ticfile:
            filewriter = csv.writer(ticfile, delimiter='\n', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            filewriter.writerow(final_tl)

    if args.verbose:
        print("TIC file written to \'" + data_outdir + "/" + filepath + ".tic\'")


# Writes a FASTA file containing the sequence of all identified occurances of the query in the reference longer than #
# the specified integer. If the input is a multiple FASTA file with sequences of varying length, the program still   #
# uses the same cutoff for all elements. As such, this is most useful with elements of similar length or a single    #
# FASTA input. All sequences are written to a single FASTA file.                                                     #
def write_qhit_fasta_new(contigs, aln_files, data_outdir, filepath, args):
    sseqid = args.outfmt["sseqid"]  # 1
    sstart = args.outfmt["sstart"]  # 2
    send = args.outfmt["send"]  # 3
    qseqid = args.outfmt["qseqid"]
    sstrand = args.outfmt["sstrand"]  # 7

    if args.verbose:
        print("Writing fasta file containing alignments longer than %s bp...." % args.writefasta)
    tig_seqdic = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta'))
    fasta_seqs = []
    file_chunk = 0
    cycle_count = 0
    const_cc = 0
    # goes through each tig in the temp directory
    for tig in aln_files:
        with open(tig, "r") as f:
            while True:
                # reads line-by-line to reduce memory load
                line = f.readline()
                if not line:
                    break

                line = line.rstrip().split('\t')

                # casts necessary parts to ints
                line[sstart] = int(line[sstart])
                line[send] = int(line[send])

                if line[send] - line[sstart] >= args.writefasta:
                    cycle_count += 1
                    const_cc += 1
                    if line[sstrand] == "plus":
                        fasta_seqs.append(SeqRecord(Seq(str(tig_seqdic[line[sseqid]].seq)[line[sstart]:line[send]]),
                                                    id=line[qseqid] + "_%s" % dv.PROG_NAME + str(const_cc),
                                                    description="%s(%s:%s) %s" %
                                                                (line[sseqid], line[sstart], line[send],
                                                                 line[send] - line[sstart]),
                                                    annotations={"molecule_type": "DNA"}))
                    else:
                        fasta_seqs.append(SeqRecord(
                            Seq(str(tig_seqdic[line[sseqid]].seq)[line[sstart]:line[send]]).reverse_complement(),
                            id=line[qseqid] + "_%s" % dv.PROG_NAME + str(const_cc),
                            description="%s(%s:%s) %s" % (line[sseqid], line[sstart], line[send],
                                                          line[send] - line[sstart]),
                            annotations={"molecule_type": "DNA"}))
                if cycle_count >= dv.FLUSH_LEN:
                    file_chunk += 1
                    SeqIO.write(fasta_seqs, "%s/%s_%s__%s.fasta" %
                                (dv.PROG_TEMP_DIR, filepath, str(args.writefasta), file_chunk), 'fasta')
                    cycle_count = 0
                    fasta_seqs = []
    if cycle_count > 0:
        file_chunk += 1
        SeqIO.write(fasta_seqs, "%s/%s_%s__%schunk.fasta" %
                    (dv.PROG_TEMP_DIR, filepath, str(args.writefasta), file_chunk), 'fasta')

    fasta_fps = ["%s/%s" % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR) if "chunk.fasta" in i]

    # glob all FASTA chunks together in the outdir
    with open(data_outdir + "/" + filepath + "_" + str(args.writefasta) + ".fasta", 'wb') as fnl_fna:
        for out_file in fasta_fps:
            with open(out_file, 'rb') as of:
                copyfileobj(of, fnl_fna)

    # remove temp FASTA chunks
    for fp in fasta_fps:
        os.remove(fp)

    if args.verbose:
        print("FASTA file written to %s" % (data_outdir + "/" + filepath + "_" + str(args.writefasta) + ".fasta"))


# exports data generated by FACET
def data_export_driver_new(contigs, aln_files, data_outdir, filepath, query, args):
    if args.writegff is False and args.writesam is False and args.writebam is False and args.nocsv is False \
            and args.writetic is False and args.writefasta == 0:  # and args.derip is False
        print("FATAL: No files to export!")
        exit()

    # sets combined hits equal to all the files in the aln_dir
    aln_files = ["%s/%s" % (aln_files, i) for i in os.listdir(aln_files)]

    # counts number of BLAST hits that remain after cleaning
    if args.verbose:
        hit_count = 0
        for tigs in aln_files:
            # iterate through each file and count lines
            with open(tigs, "r+") as f:
                f_buffer = mmap(f.fileno(), 0)
                readline = f_buffer.readline
                while readline():
                    hit_count += 1
        print("%s BLAST hits remain after cleaning" % hit_count)
        print("\nBeginning data export...")

    # writes tic file #
    if args.writetic:
        write_tic_new(contigs, aln_files, data_outdir, filepath, args)

    # writes multiple fasta containing hit sequences (all sequences are on the plus strand) #
    if args.writefasta > 0:
        write_qhit_fasta_new(contigs, aln_files, data_outdir, filepath, args)

    # writes allhits file #
    if args.nocsv:
        write_csv_new(aln_files, data_outdir, filepath, args)

    # writes GFF file #
    if args.writegff:
        write_gff_new(aln_files, filepath, data_outdir, str(Path(query).resolve().absolute()), args)

    # writes BAM and SAM files #
    if args.writesam or args.writebam:
        write_sambam_new(contigs, aln_files, data_outdir, filepath, args)

    if args.verbose:
        print("Data export finished!")
