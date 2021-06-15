#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import subprocess
import os
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import btop_utils
from .defaults import ProgDefaults as dv


# checks each of the files that would be created and sees if the file already exists #
# If it does, refuses to overwrite unless --force is specified                       #
def check_file_overwrite(datafilepath, filepath, args):
    ok_check = True
    if args.writegff:
        if Path(datafilepath + "/" + filepath + ".gff").exists():
            ok_check = False
    if args.writesam:
        if args.verboseoutput:
            if Path(datafilepath + "/" + filepath + "_verbose.sam").exists():
                ok_check = False
        else:
            if Path(datafilepath + "/" + filepath + ".sam").exists():
                ok_check = False
    if args.writebam:
        if args.verboseoutput:
            if Path(datafilepath + "/" + filepath + "_verbose.bam").exists() or \
                    Path(datafilepath + "/" + filepath + "_verbose.bai").exists():
                ok_check = False
        else:
            if Path(datafilepath + "/" + filepath + ".bam").exists() or \
                    Path(datafilepath + "/" + filepath + ".bai").exists():
                ok_check = False
    if args.nocsv:
        if args.verboseoutput:
            if Path(datafilepath + "/" + filepath + "_verbose.csv").exists():
                ok_check = False
        else:
            if Path(datafilepath + "/" + filepath + ".csv").exists():
                ok_check = False
    """
    if args.derip:
        if Path(datafilepath + "/" + queryFilename + "_derip_" + str(
                int(args.deripswitch * 100)) + "perc.fasta").exists() or \
                Path(datafilepath + "/" + queryFilename + "_new_consensus.fasta").exists():
            ok_check = False
    """
    if args.writetic:
        if Path(datafilepath + "/" + filepath + ".tic").exists():
            ok_check = False
    if args.writefasta > 0:
        if Path(datafilepath + "/" + filepath + "_" + str(args.writefasta) + ".fasta").exists():
            ok_check = False
    return ok_check


# Writes a gff file containing the locations of query elements in the reference #
def write_gff(flatcombinedhits, filepath, datafilepath, queryfile, args):
    sseqid = dv.SSEQID_INDEX + 1  # 1
    sstart = dv.SSTART_INDEX + 1  # 2
    send = dv.SEND_INDEX + 1  # 3
    qseqid = dv.QSEQID_INDEX + 1  # 4
    qstart = dv.QSTART_INDEX + 1  # 5
    qend = dv.QEND_INDEX + 1  # 6
    sstrand = dv.SSTRAND_INDEX + 1  # 7
    if args.subparser_id == "outfile":
        sseqid = dv.FMTSIX_SSEQID + 1
        sstart = dv.FMTSIX_SSTART + 1
        send = dv.FMTSIX_SEND + 1
        qseqid = dv.FMTSIX_QSEQID + 1
        qstart = dv.FMTSIX_QSTART + 1
        qend = dv.FMTSIX_QEND + 1
        sstrand = len(flatcombinedhits[0]) - 1  # should always be last element because it is added at the end
    if args.verbose:
        print("Writing " + str(len(flatcombinedhits)) + " hits to \'" + datafilepath + "/" + filepath + ".gff\'")
    querydict = SeqIO.to_dict(SeqIO.parse(queryfile, 'fasta'))
    for i in querydict.keys():
        querydict[i] = len(querydict[i].seq)

    gfflist = [["##gff-version 3"]]
    for hits in flatcombinedhits:

        ele_len = querydict[hits[qseqid]]

        # if hit contains both the beginning and the end, make it cyan
        if hits[qend] >= ele_len - args.buffer and hits[qstart] <= args.buffer:
            color_str = dv.BOTH_ENDS_COLOR
        # if hit contains only the end, make it yellow
        elif hits[qend] >= ele_len - args.buffer:
            color_str = dv.ONLY_END_COLOR
        # if hit contains only the beginning, make it orange
        elif hits[qstart] <= args.buffer:
            color_str = dv.ONLY_BEGINNING_COLOR
        # if hit contains neither, make it grey
        else:
            color_str = dv.NEITHER_END_COLOR

        if hits[sstrand] == 'plus':
            strand = '+'
        else:
            strand = '-'

        gfflist.append([hits[sseqid], "%s" % dv.PROG_NAME, hits[qseqid], hits[sstart], hits[send], '.', strand, '.',
                        'color=%s; ID=%s;Name=%s;Query=start: %s, end: %s;Subject=start: %s, end: %s;Length=%s' %
                        (color_str, hits[0], hits[qseqid], hits[qstart], hits[qend], hits[sstart], hits[send],
                         hits[send] - hits[sstart])])

    with open(datafilepath + "/" + filepath + ".gff", 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for hits in gfflist:
            filewriter.writerow(hits)
    gfflist = None
    del gfflist


# Writes a sam/bam file containing the locations of query elements in the reference #
# This file can include sequence/alignment info if -q is specified                  #
def write_sambam(contigs, flatcombinedhits, datafilepath, filepath, args):
    samlist = [["@HD", "VN:1.6"]]
    for record in SeqIO.parse(contigs, "fasta"):
        samlist.append(["@SQ", "SN:" + record.id, "LN:" + str(len(record.seq))])
    if args.verboseoutput:
        if args.verbose:
            print("Parsing BTOP sequences and converting to CIGAR sequences....")
    query_dict = SeqIO.to_dict(SeqIO.parse(args.query, 'fasta'))
    for hits in flatcombinedhits:
        qseq = str(query_dict[hits[dv.QSEQID_INDEX + 1]].seq)[hits[dv.QSTART_INDEX + 1] - 1:hits[dv.QEND_INDEX + 1]]
        if hits[dv.SSTRAND_INDEX + 1] == 'plus':
            flag = 0
            if args.verboseoutput:
                cigar_seq = btop_utils.btop_to_cigar(btop_utils.parse_btop(hits[dv.BTOP_INDEX + 1]))

        else:
            flag = 16
            if args.verboseoutput:
                # BTOP sequence needs to be reverse complimented, +/-
                cigar_seq = btop_utils.btop_to_cigar(btop_utils.parse_btop(btop_utils.reverse_compliment(hits[9])))
                qseq = str(Seq(qseq).reverse_complement())
        if args.verboseoutput:
            qseq = qseq.replace("-", "")
            samlist.append([hits[0], flag, hits[dv.SSEQID_INDEX + 1], hits[dv.SSTART_INDEX + 1], 255, cigar_seq, "*",
                            0, 0, qseq, "*",
                            "CO:Z:query(" + hits[dv.QSEQID_INDEX + 1] + "), qloc(" + str(hits[dv.QSTART_INDEX + 1]) +
                            "," + str(hits[dv.QEND_INDEX + 1]) + ")," + " sloc(" + str(hits[dv.SSTART_INDEX + 1]) + ","
                            + str(hits[dv.SEND_INDEX + 1]) + "), len(" + str(hits[dv.SEND_INDEX + 1] -
                                                                             hits[dv.SSTART_INDEX + 1] + 1) + ")"])
        else:
            samlist.append([hits[0], flag, hits[dv.SSEQID_INDEX + 1], hits[dv.SSTART_INDEX + 1], 255,
                            str(hits[dv.SEND_INDEX + 1] - hits[dv.SSTART_INDEX + 1] + 1) + "M", "*", 0, 0, "*", "*",
                            "CO:Z:query(" + hits[dv.QSEQID_INDEX + 1] + "), qloc(" + str(hits[dv.QSTART_INDEX + 1]) +
                            "," + str(hits[dv.QEND_INDEX + 1]) + ")," + " sloc(" + str(hits[dv.SSTART_INDEX + 1]) + ","
                            + str(hits[dv.SEND_INDEX + 1]) + "), len(" + str(hits[dv.SEND_INDEX + 1] -
                                                                             hits[dv.SSTART_INDEX + 1] + 1) + ")"])

    if args.writebam and not args.writesam:
        samstring = datafilepath + "/" + "FACET_temp_verbose.sam"
    else:
        if args.verboseoutput:
            if args.verbose:
                print("Writing " + str(len(flatcombinedhits)) + " hits to \'" + datafilepath + "/" + filepath
                      + "_verbose.sam\'")
            samstring = datafilepath + "/" + filepath + "_verbose.sam"
        else:
            if args.verbose:
                print("Writing " + str(len(flatcombinedhits)) + " hits to \'" + datafilepath + "/" + filepath
                      + ".sam\'")
            samstring = datafilepath + "/" + filepath + ".sam"
    with open(samstring, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for hits in samlist:
            filewriter.writerow(hits)
    if args.writebam:
        if args.verboseoutput:
            file_end = "_verbose.bam"
        else:
            file_end = ".bam"
        if args.verbose:
            print("Writing " + str(len(flatcombinedhits)) + " hits to \'" + datafilepath + "/" + filepath
                  + file_end + "\'")
        subprocess.call(["samtools", "view", "-bS", samstring],
                        stdout=open(datafilepath + "/" + filepath + file_end, 'w'),
                        stderr=open(os.devnull, 'w'))
        if args.verbose:
            print("Indexing bam file....")
        subprocess.call(["samtools", "index", datafilepath + "/" + filepath + file_end,
                         datafilepath + "/" + filepath + file_end[:-1] + "i"])
        if not args.writesam:
            os.remove(samstring)
    samlist = None
    del samlist


# Writes a tab delimited csv file containing information on the location and % ident of query elements in the ref #
def write_csv(flatcombinedhits, datafilepath, filepath, args):
    """
    if args.derip:
        slice_num = 1

    if args.verboseoutput:
        slice_num = 0  # USED TO BE 3
    """
    if args.verboseoutput:
        csvname = datafilepath + "/" + filepath + "_verbose.csv"
    else:
        csvname = datafilepath + "/" + filepath + ".csv"
    if args.verbose:
        print("Writing " + str(len(flatcombinedhits)) + " hits to \'" + csvname + "\'")
    with open(csvname, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for hits in flatcombinedhits:
            """
            if args.verboseoutput:  # or args.derip
                filewriter.writerow(hits[:-slice_num])
            else:
            """
            filewriter.writerow(hits)


# Writes a TIg Coverage file that contains information about coverage depth on the contig.                   #
# Really only useful in Genome x Genome comparisons that aren't cleaned to try to find the location of novel #
# repeat elements and inform evolutionary history of existing insertions.                                    #
# This file format is in development and is currently slightly stable                                        #
def write_tic(contigs, combinedhits, datafilepath, filepath, args):
    sstart = dv.SSTART_INDEX + 1  # 2
    send = dv.SEND_INDEX + 1  # 3
    if args.subparser_id == "outfile":
        sstart = dv.FMTSIX_SSTART + 1
        send = dv.FMTSIX_SEND + 1

    # the following code produces a dictionary called tig_lengths such that {"tig01": tiglength} #
    tig_lengths = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta'))
    tempdict = {}
    for tignames in tig_lengths.keys():
        tempdict[tignames] = len(tig_lengths[tignames].seq)
    tig_lengths = tempdict
    tempdict = None
    del tempdict

    if args.verbose:
        print("Finding tic file coverage....")
    for tig_names in combinedhits.keys():

        if args.verbose:
            print("writing %s to tic file...." % str(tig_names))
        """
        # adding start (a) and stop (o) locations #
        a_list = []
        o_list = []
        for hit_records in combinedhits[tig_names]:
            a_list.append(hit_records[2] - 1)
            o_list.append(hit_records[3] - 1)
        a_list = [ele for ele in Counter(a_list).items()]
        o_list = [ele for ele in Counter(o_list).items()]
        """

        tic_coverage = [0] * tig_lengths[tig_names]  # generates a list of zeroes as long as the tig

        # adds ranges from combined hits to the tic_coverage list #
        # creates an operational list that needs to be converted  #
        for hit_records in combinedhits[tig_names]:
            tic_coverage[hit_records[sstart] - 1] += 1
            if hit_records[send] - 1 < tig_lengths[tig_names] - 1:  # CHANGED OPERATOR FROM != TO <
                tic_coverage[hit_records[send]] -= 1

        # converts operational list to actual coverage list #
        for basepair in range(1, tig_lengths[tig_names]):
            tic_coverage[basepair] += tic_coverage[basepair - 1]
        """
        tic_coverage = [[basepair, 0, 0] for basepair in tic_coverage]  # casts all elements to lists

        for ao_sites in range(0, max(len(a_list), len(o_list))):  # adds a/o sites
            try:
                tic_coverage[a_list[ao_sites][0]][1] = a_list[ao_sites][1]
            except IndexError:
                pass
            try:
                tic_coverage[o_list[ao_sites][0]][2] = o_list[ao_sites][1]
            except IndexError:
                pass
        """
        # tic_coverage = tic_coverage[5191538:5197197]

        # adds the beginning of every contig to the tic file #
        final_tl = ["1:%s" % (tic_coverage[0])]
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
        # adds the end of every contig to the tic tile #
        final_tl.append("%s:%s" % (len(tic_coverage), tic_coverage[len(tic_coverage) - 1]))

        # clears the tic_coverage variable every iteration to free up memory #
        tic_coverage = None

        with open(datafilepath + "/FACET_temp_tic.tic", 'a', newline='') as ticfile:
            filewriter = csv.writer(ticfile, delimiter='\n', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow([">" + tig_names] + final_tl)
        # clears the final_tl variable every iteration to free up memory #
        final_tl = None
    os.rename(datafilepath + "/FACET_temp_tic.tic", datafilepath + "/" + filepath + ".tic")  # renames temp tic file
    # clears variables and deletes them to free up memory #
    tig_lengths = None
    del tig_lengths
    del tic_coverage
    del final_tl
    if args.verbose:
        print("TIC file written to \'" + datafilepath + "/" + filepath + ".tic\'")


# Writes a FASTA file containing the sequence of all identified occurances of the query in the reference longer than #
# the specified integer. If the input is a multiple FASTA file with sequences of varying length, the program still   #
# uses the same cutoff for all elements. As such, this is most useful with elements of similar length or a single    #
# FASTA input. All sequences are written to a single FASTA file.                                                     #
def write_qhit_fasta(contigs, flatcombinedhits, datafilepath, filepath, args):
    sseqid = dv.SSEQID_INDEX + 1  # 1
    sstart = dv.SSTART_INDEX + 1  # 2
    send = dv.SEND_INDEX + 1  # 3
    sstrand = dv.SSTRAND_INDEX + 1  # 7
    if args.subparser_id == "outfile":
        sseqid = dv.FMTSIX_SSEQID + 1
        sstart = dv.FMTSIX_SSTART + 1
        send = dv.FMTSIX_SEND + 1
        sstrand = len(flatcombinedhits[0]) - 1  # should always be last element because it is added at the end

    if args.verbose:
        print("Writing fasta file containing hit sequences....")
    tig_seqdic = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta'))
    fasta_seqs = []
    count = 0
    for hit in flatcombinedhits:
        if hit[send] - hit[sstart] >= args.writefasta:
            count += 1
            if hit[sstrand] == "plus":
                fasta_seqs.append(SeqRecord(Seq(str(tig_seqdic[hit[sseqid]].seq)[hit[sstart]:hit[send]]),
                                            id=hit[4] + "_" + str(count),
                                            description="%s(%s:%s),%s" %
                                                        (hit[sseqid], hit[sstart], hit[send], hit[send] - hit[sstart]),
                                            annotations={"molecule_type": "DNA"}))
            else:
                fasta_seqs.append(SeqRecord(Seq(str(tig_seqdic[hit[sseqid]].seq)[hit[sstart]:hit[send]]).reverse_complement(),
                                            id=hit[4] + "_" + str(count),
                                            description="%s(%s:%s),%s" %
                                                        (hit[sseqid], hit[sstart], hit[send], hit[send] - hit[sstart]),
                                            annotations={"molecule_type": "DNA"}))
    SeqIO.write(fasta_seqs, datafilepath + "/" + filepath + "_" + str(args.writefasta) + ".fasta", "fasta")
    if args.verbose:
        print("Wrote %s hits to %s" % (count, datafilepath + "/" + filepath + "_" +
                                       str(args.writefasta) + ".fasta"))


# exports data generated by FACET
def data_export_driver(contigs, combinedhits, datafilepath, filepath, query, args):
    if args.writegff is False and args.writesam is False and args.writebam is False and args.nocsv is False \
            and args.writetic is False and args.writefasta == 0:  # and args.derip is False
        if args.verbose:
            print("FATAL: No files to export!")
        exit()

    if args.verbose:
        hit_count = 0
        for tigs in combinedhits.keys():
            hit_count += len(combinedhits[tigs])
        print("%s BLAST hits remain after cleaning" % hit_count)
        print("\nBeginning data export...")

    # writes tic file #
    if args.writetic:
        write_tic(contigs, combinedhits, datafilepath, filepath, args)

    # flattens the combined hits dictionary #
    flatcombinedhits = []
    for tigs in combinedhits.keys():
        flatcombinedhits.extend(combinedhits[tigs])
    combinedhits = None
    del combinedhits

    # writes multiple fasta containing hit sequences (all sequences are on the plus strand) #
    if args.writefasta > 0:
        write_qhit_fasta(contigs, flatcombinedhits, datafilepath, filepath, args)

    # writes allhits file #
    if args.nocsv:
        write_csv(flatcombinedhits, datafilepath, filepath, args)

    # writes GFF file #
    if args.writegff:
        write_gff(flatcombinedhits, filepath, datafilepath, str(Path(query).resolve().absolute()), args)

    # writes BAM and SAM files #
    if args.writesam or args.writebam:
        write_sambam(contigs, flatcombinedhits, datafilepath, filepath, args)

    """
    # de-RIP process #
    if args.derip:
        derip_query(query, flatcombinedhits, datafilepath, queryFilename)
    """

    if args.verbose:
        print("Data export finished!")
