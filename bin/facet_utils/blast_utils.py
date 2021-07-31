#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import subprocess
import os
import csv
import multiprocessing as mp
from pathlib import Path
from shutil import copyfileobj, rmtree

from Bio import SeqIO

from .defaults import ProgDefaults as dv
from . import btop_utils
from . import masker_utils


# takes a list of lists and retains elements that have no complete overlap with other elements  #
# Input looks like: [[9,20], [2,3],[5,7],[2,8],[1,9],[1,10],[8,15]], returns [[1,10],[9,20]]    #
# index1 and index2 allow you to have lists that have more than two elements but you still want #
# to compare them as if they did. Input ranges must be [min, max]                               #
def clean_list_new(cluttered, index1, index2, pident_index, pident_check, args):
    # sorts cluttered list by first part of each sublist
    cluttered = sorted(cluttered, key=lambda x: x[index1])
    # reduces function calls #
    clutteredremove = cluttered.remove

    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within other ranges remain

        exitcon = False  # sets exit condition to false each time the loop starts

        # loop checks to see if either list is containing the other #
        i = 1
        while i < len(cluttered):
            if len(cluttered) <= 1:
                break

            curr_contained = False
            prev_contained = False

            # checks to see if cluttered[i] is contained within a buffered cluttered[i - 1]
            if cluttered[i][index1] >= cluttered[i - 1][index1] - args.buffer and \
                    cluttered[i][index2] <= cluttered[i - 1][index2] + args.buffer:
                curr_contained = True

            # checks to see if cluttered[i - 1] is contained within a buffered cluttered[i]
            if cluttered[i - 1][index1] >= cluttered[i][index1] - args.buffer and \
                    cluttered[i - 1][index2] <= cluttered[i][index2] + args.buffer:
                prev_contained = True

            # checks to see if hits are approx. the same length #
            if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                    (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer and \
                    (prev_contained or curr_contained):

                # if we're looking at pident
                if pident_check:
                    # if cluttered[i - 1] has a lower pident, remove it #
                    if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                        clutteredremove(cluttered[i - 1])

                    # if the pidents are the same, keep the longer hit #
                    elif cluttered[i][pident_index] - cluttered[i - 1][pident_index] == 0:
                        if (cluttered[i][index2] - cluttered[i][index1]) > \
                                (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                            clutteredremove(cluttered[i - 1])
                        else:
                            clutteredremove(cluttered[i])

                    # if cluttered[i] has a lower pident, remove it
                    else:
                        clutteredremove(cluttered[i])

                # if we're looking at length
                else:
                    # if cluttered[i] is longer, remove cluttered[i - 1]
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])

                    # If both alignments are the same length, check pident
                    elif (cluttered[i][index2] - cluttered[i][index1]) == \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        else:
                            clutteredremove(cluttered[i])
                    # if cluttered[i-1] is longer, remove cluttered[i]
                    else:
                        clutteredremove(cluttered[i])
                # keeps the outer loop going
                exitcon = True

                # decrements iterator because an element was removed
                i -= 1
            # if the hits are not ~ the same length, remove hit based on beginning checks
            else:
                # checks to see if cluttered[i] is contained within a buffered cluttered[i - 1]
                if curr_contained:
                    clutteredremove(cluttered[i])
                    # keeps the loop going
                    exitcon = True

                    # decrements iterator because an element was removed
                    i -= 1

                # checks to see if cluttered[i - 1] is contained within a buffered cluttered[i]
                elif prev_contained:
                    clutteredremove(cluttered[i - 1])
                    # keeps the loop going
                    exitcon = True

                    # decrements iterator because an element was removed
                    i -= 1
            # increments iterator each loop
            i += 1

    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within two surrounding ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        i = 1
        while i < len(cluttered) - 1:

            # This situation literally cannot happen with 2 or fewer alignments, so just break out of loop
            if len(cluttered) <= 2:
                break

            # if the range in question is overlapping both the one before and after it
            # AND the range in question is contained by the range before and after
            if max(cluttered[i][index1], cluttered[i - 1][index1]) <= \
                    min(cluttered[i][index2], cluttered[i - 1][index2]) and \
                    max(cluttered[i][index1], cluttered[i + 1][index1]) <= \
                    min(cluttered[i][index2], cluttered[i + 1][index2]) and \
                    cluttered[i - 1][index2] >= cluttered[i + 1][index1]:
                # keeps the loop going
                exitcon = True

                # removes the hit contained within two other hits
                clutteredremove(cluttered[i])

                # decrements the iterator if an element is removed
                i -= 1

            # increments the iterator each loop cycle
            i += 1

    return cluttered


# takes the name of the new directory and the file used to create the database and makes a blast database #
# also returns the filepath of the newly created blast database                                           #
def make_blast_db(dbname, dbfile, out_dir, args):
    outfilename = out_dir + "/" + dbname + "/" + str(
        Path(Path(dbfile).resolve().name).with_suffix('')) + "_db.fasta"
    try:
        os.mkdir(out_dir + "/" + dbname)
    except FileExistsError:
        if args.verbose:
            print("Directory \'" + dbname + '\' already exists; using it')
    # checks to see if a blast database of the same name has already been created at the target #
    # If it has, checks to see if permission to overwrite has been granted                      #
    if Path(outfilename + ".nhr").exists() and Path(outfilename + ".nin").exists() and \
            Path(outfilename + ".nsq").exists() and args.overwritedb is False:
        if args.verbose:
            print("Using existing database, as permission to overwrite has not been granted.\n"
                  "Use \'--overwritedb\' if you want to overwrite the existing BLAST database")
    else:
        p1 = subprocess.run(["makeblastdb", "-in", dbfile, "-dbtype", "nucl", "-out", outfilename],
                            stdout=subprocess.PIPE, universal_newlines=True, env={'PATH': os.getenv('PATH')})
        if p1.returncode != 0:
            print("FATAL: Unable to create BLAST database!")
            exit()
        elif args.verbose:
            print("\nBLAST database created at %s" % outfilename)
    return outfilename


# gets infile blast dictionary and cleans it up #
def infile_driver(blast_infile, args):
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]
    pident = args.outfmt["pident"]

    if args.verboseoutput:
        try:
            args.outfmt["btop"]
        except KeyError:
            print("FATAL: outfmt does not contain btop and \'--verboseoutput\' specifed!")
            print("       provided outfmt does not contain btop information and cannot be"
                  " used to make verbose outputs!")
            exit()

    verify_outfmt(blast_infile, args)

    # flush alignments to tig files
    if args.verbose:
        print("\nFlushing alignments from master file based on contig ID....")
    flush_outfile_to_tigs(blast_infile, dv.PROG_TEMP_DIR, "sseqid", True, blast_infile, args)

    # master becomes a list of paths to the files in the temp_dir
    master_file = ['%s/%s' % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR)]

    # Parses each tig
    if args.verbose:
        print("Parsing each contig file....")
    multi_pool = mp.Pool(processes=dv.AVAIL_CORES)
    for tmpfile in master_file:
        # if we aren't cleaning, no need to break up by query
        if args.noclean:
            multi_pool.apply_async(write_blast_outfile, args=(parse_blast_outfile(tmpfile, args), tmpfile, ))
        # if we are cleaning, break each tig up by query
        else:
            multi_pool.apply_async(individual_tig_parser, args=(tmpfile, args, sstart, send, pident, ))

    multi_pool.close()
    multi_pool.join()
    del multi_pool

    # concatenates repeat data if noclean and nocat have not been called
    if args.nocat is False and args.noclean is False:
        if args.verbose:
            print("\nConcatenating repeat data....")
        for contig_files in master_file:
            grimy = parse_blast_outfile(contig_files, args)

            grimy = clean_list_new(grimy, sstart, send, pident, True, args)

            # write the cleaned up file for each contig
            with open(contig_files, 'w') as k:
                filewriter = csv.writer(k, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                        lineterminator="\n")
                for records in grimy:
                    filewriter.writerow(records)

    if args.verbose:
        print("\nSorting each contig's alignments by sstart....")


    # sort each contig's alignments by sstart
    multi_pool = mp.Pool(processes=dv.AVAIL_CORES)
    for outfile_name in ["%s/%s" % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR)]:
        multi_pool.apply_async(sstart_sort_outfile, args=(outfile_name, args,))

    multi_pool.close()
    multi_pool.join()
    del multi_pool

    if args.verbose:
        print()
    # all hits are stored in dv.PROG_TEMP_DIR, just operate on that dir
    return dv.PROG_TEMP_DIR


def individual_tig_parser(temp_tig_file, args, sstart, send, pident):
    curr_id_dir = "%s/%s_outdir" % (dv.PROG_TEMP_DIR, Path(temp_tig_file).stem)
    # make a dir to store files exported based on qseqid
    Path(curr_id_dir).mkdir(parents=True, exist_ok=True)

    # flush alignments to query files and delete the master
    flush_outfile_to_tigs(temp_tig_file, curr_id_dir, "qseqid", False, temp_tig_file, args)

    qfiles = ['%s/%s' % (curr_id_dir, i) for i in os.listdir(curr_id_dir)]
    # open file containing a single contig's matches with a single query
    for qfile in qfiles:
        grimy = parse_blast_outfile(qfile, args)
        # if we are cleaning, clean up the grimy list; else, do no cleaning
        if args.noclean is False:
            # pident is set to -1 here because the larger hit should be picked, not the one with the highest pident
            grimy = clean_list_new(grimy, sstart, send, pident, False, args)

        # write the cleaned up file
        with open(qfile, 'w') as k:
            filewriter = csv.writer(k, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                    lineterminator="\n")
            for records in grimy:
                filewriter.writerow(records)

    # combines all of the queryid files into an outfile for the contig
    with open(temp_tig_file, 'wb') as tmp:
        for out_file in qfiles:
            with open(out_file, 'rb') as of:
                copyfileobj(of, tmp)
    rmtree(curr_id_dir)


# makes sure the first entry of the blast_list is in the correct format #
def verify_outfmt(blast_list, args):
    bl = open(blast_list, 'r')
    line = bl.readline().rstrip()
    bl.close()
    try:
        check_list = line.split('\t')
    except Exception:
        print("\nFATAL: infile could not be accessed or is empty")
        exit()

    if "blastHit_" in check_list[0] or "shortHit_" in check_list[0]:
        check_list.pop(0)

    # casts each part of the blast list to the defined type in the default values dictionary
    outfmt_specs = list(args.outfmt.keys())
    for j in range(0, len(outfmt_specs)):
        try:
            check_list[j] = dv.BLASTN_TYPES[outfmt_specs[j]](check_list[j])
        except Exception:
            print("FATAL: first row of input does not match provided outfmt:")
            print("\"%s\"" % generate_outfmt_str(args))
            exit()

    if args.verboseoutput:
        try:
            btop_utils.btop_to_cigar(btop_utils.parse_btop(check_list[args.outfmt["btop"]]))
        except Exception:
            print("FATAL: BTOP cannot be converted to CIGAR")
            exit()


# parse outfmt strings to allow users to use custom outfmts #
def outfmt_parser(args):
    # if the user is just using the facet outfmt, set the outfmt to the base one
    if args.outfmt.lower() == 'facet':
        args.outfmt = dv.BASE_INDX
        # if btop is not required, remove it from the dictionary
        if args.subparser_id not in dv.MASKER_ALIAS and args.subparser_id not in dv.VC_ALIAS:
            if not args.verboseoutput:
                del args.outfmt["btop"]
        else:
            del args.outfmt["btop"]

    # if the user is using outfmt 6, set the outfmt indices using outfmt six
    elif args.outfmt == '6' or args.outfmt.lower() == 'six':
        args.outfmt = dv.OUT_SIX_INDX

    # if the user is using a custom outfmt, make sure it is formatted correctly and contains necessary data
    else:
        outfmt_list = args.outfmt.split()

        # makes sure outfmt begins with a six, like blastn wants
        if outfmt_list[0] == '6' or outfmt_list[0].lower() == 'six':
            custom_outfmt = {}

            # runs through each of the provided outfmt specs
            for fmt_spec in range(1, len(outfmt_list)):
                # checks that the spec is in the list of acceptable specs
                if outfmt_list[fmt_spec].lower() in dv.BLASTN_TYPES.keys():
                    # if it is, adds it to the dictionary with its value as its index
                    custom_outfmt.setdefault(outfmt_list[fmt_spec].lower(), fmt_spec - 1)

                # if it isnt, return a fatal error
                else:
                    print("FATAL: Invalid format specifier (%s)" % outfmt_list[fmt_spec])
                    exit()

            # defines required keys based on parser type and options
            if args.subparser_id in dv.MASKER_ALIAS:
                required_keys = dv.MASKER_OUTFMT.split()[1:]
            else:
                if args.verboseoutput:
                    required_keys = dv.VERBOSE_OUTFMT_STR.split()[1:]
                else:
                    required_keys = dv.BASE_OUTFMT_STR.split()[1:]
                required_keys.remove("sstrand")  # the sstrand is not required because it will be added if not present

            # runs through each of the required keys to make sure they are present in the provided outfmt
            for spec in required_keys:
                try:
                    custom_outfmt[spec]
                except KeyError:
                    print("FATAL: Required outformat specifier is not present (%s)" % spec)
                    print("The following outformat specifiers are required: ")
                    print(" ".join(required_keys))
                    exit()

            # if all checks are passed, sets the outfmt arg equal to the custom outfmt
            args.outfmt = custom_outfmt
        # if outfmt isn't specified correctly, fatal error
        else:
            print("FATAL: invalid custom outformat! Custom outfmt should look something like this: ")
            print("\"6 sseqid qseqid sstart sstop qstart qstop\"")
            exit()


# generate an outfmt string from the keys in an outfmt dictionary
def generate_outfmt_str(args):
    outfmt_str = "6"
    for specs in args.outfmt.keys():
        outfmt_str += " %s" % specs
    return outfmt_str


def blast_to_outfile(db, query, outpath, args):
    if args.verbose:
        print("\nRunning BLASTn with %s as db and %s as query...." % (Path(db).resolve(), Path(query).resolve()))
    # blast run contains only the necessary information
    blastcommand = ['blastn', '-db', str(Path(db).resolve()), "-query", query]
    if args.num_threads > 0:
        blastcommand.extend(["-num_threads", str(args.num_threads)])

    outfmt_string = generate_outfmt_str(args)

    blastcommand.extend(["-evalue", args.blaste, "-outfmt", outfmt_string, "-out", "%s/%s.out" %
                         (outpath, Path(query).stem)])
    blast = subprocess.run(blastcommand, stdout=subprocess.PIPE, universal_newlines=True)
    blastcommand[len(blastcommand) - 1] = "\"" + blastcommand[len(blastcommand) - 1] + "\""
    if args.verbose:
        print("BLAST was run with these options: ")
        print(" ".join(blastcommand))


# takes a blastn outfile and flushes its contents to separate files based on flush key [sseqid or qseqid] #
def flush_outfile_to_tigs(outfile, flush_dir, flush_key, add_file_end, name_file, args):
    tig_dict = {}
    cycle_count = 0

    outfile_name = Path(name_file).stem.replace(".", "")  # removes . from outfile name

    with open(outfile, 'r') as f:
        # iterate through the file
        while True:
            # read each line
            line = f.readline()
            # if the line is the last line of the file, exit the loop
            if not line:
                break
            else:
                cycle_count += 1
                # parse the line
                line = line.rstrip().split("\t")

                # if the line is empty, don't try to parse it
                if len(line) == 0:
                    cycle_count -= 1
                    continue

                # add the contig to the tig dict if it isn't already there
                tig_dict.setdefault(line[args.outfmt[flush_key]], [])

                # add the line to the tig dict
                tig_dict[line[args.outfmt[flush_key]]].append(line)

                # write each tig to a separate file
                if cycle_count >= dv.FLUSH_LEN:
                    for tig_id in tig_dict.keys():
                        if add_file_end:
                            temp_file_name = "%s/%s_%s_%s_%stemp.%s" % \
                                             (flush_dir, outfile_name, tig_id.replace(".", ""),
                                              dv.TIME_STR, dv.PROG_NAME, tig_id)
                        else:
                            temp_file_name = "%s/%s_%s_%s_%stemp" % \
                                             (flush_dir, outfile_name, tig_id.replace(".", ""), dv.TIME_STR,
                                              dv.PROG_NAME)
                        with open(temp_file_name, 'a') as csvfile:
                            filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                                    lineterminator="\n")
                            for records in tig_dict[tig_id]:
                                filewriter.writerow(records)
                        # sets the tig list to empty after it's flushed
                        tig_dict[tig_id] = []
                    # resets cycle count
                    cycle_count = 0
        # if there are any hits remaining in the tig_dict, write them out to files
        # file look like FACET_temp_20210628120709/chr1_20210628120709_FACETtemp
        if cycle_count > 0:
            for tig_id in tig_dict.keys():
                if add_file_end:
                    temp_file_name = "%s/%s_%s_%s_%stemp.%s" % \
                                     (flush_dir, outfile_name, tig_id.replace(".", ""),
                                      dv.TIME_STR, dv.PROG_NAME, tig_id)
                else:
                    temp_file_name = "%s/%s_%s_%s_%stemp" % \
                                     (flush_dir, outfile_name, tig_id.replace(".", ""), dv.TIME_STR, dv.PROG_NAME)
                with open(temp_file_name, 'a') as csvfile:
                    filewriter = csv.writer(csvfile, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                            lineterminator="\n")
                    for records in tig_dict[tig_id]:
                        filewriter.writerow(records)


def blast_driver(dbfilepath, subject, query, args):
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]
    pident = args.outfmt["pident"]

    delete_master = True
    # if using outfile module, master file is just the blast outfile specified
    if args.subparser_id in dv.OUTFILE_ALIAS:
        master_file = args.blast_outfile
        delete_master = False

    # elif using variant calling module and outfile has been specified
    elif args.subparser_id in dv.VC_ALIAS and args.outfile != "":
        master_file = args.outfile
        delete_master = False

    # elif using masking module and outfile has been specified
    elif args.subparser_id in dv.MASKER_ALIAS and args.outfile != "":
        master_file = args.outfile
        delete_master = False

    # elif (running vc/masker without outfile OR running db/dbf) AND using large method
    elif args.large:
        delete_master = False
        # makes a temporary directory to store fasta files and blast outfiles #
        temp_fasta_fp = "%s/%s_individual_contigs" % (dv.PROG_TEMP_DIR, Path(subject).stem.replace(".", "_"))
        if Path(temp_fasta_fp).exists():
            print("FATAL: \'%s\' already exists! How did this happen?!" % temp_fasta_fp)
            exit()
        else:
            os.mkdir(temp_fasta_fp)

        # creates a fasta file for each contig present in the genome file #
        if args.verbose:
            print("\nGetting individual contigs from genome file")
        masker_utils.split_contigs(temp_fasta_fp, Path(subject).resolve())
        contig_list = ["%s/%s" % (temp_fasta_fp, i) for i in os.listdir(temp_fasta_fp)]

        # blasts each contig individually against the genome database #
        avail_blast_cores = dv.AVAIL_CORES
        try:
            avail_blast_cores = int(avail_blast_cores/args.num_threads)
        except AttributeError:
            pass

        multi_pool = mp.Pool(processes=avail_blast_cores)

        for tig in contig_list:
            multi_pool.apply_async(execute_blast, args=(dbfilepath, tig, "%s/%s" % (dv.PROG_TEMP_DIR, Path(tig).stem),
                                                        "IMPLEMENT VVQ", args, ))

        multi_pool.close()
        multi_pool.join()
        del multi_pool

        master_file = "%s/%s_%s_%s_master.out" % \
                      (dv.PROG_TEMP_DIR, Path(subject).stem, Path(query).stem, dv.TIME_STR)
        # removes the directory containing the individual contig files
        rmtree(temp_fasta_fp)

    # else (running vc/masker without outfile OR running db/dbf) AND not using large method
    else:
        # run blastn and get the filepath to the master
        master_file = "%s/%s_%s_%s_master.out" % \
                      (dv.PROG_TEMP_DIR, Path(subject).stem, Path(query).stem, dv.TIME_STR)
        master_file = execute_blast(dbfilepath, query, master_file, "IMPLEMENT VVQ", args)

    # flush alignments to tig files
    if args.large:
        indv_outfiles = ["%s/%s" % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR)]
        multi_pool = mp.Pool(processes=dv.AVAIL_CORES)
        for i in indv_outfiles:
            multi_pool.apply_async(flush_outfile_to_tigs, args=(i, dv.PROG_TEMP_DIR, "sseqid", True, master_file, args))

        multi_pool.close()
        multi_pool.join()
        del multi_pool
    else:
        flush_outfile_to_tigs(master_file, dv.PROG_TEMP_DIR, "sseqid", True, master_file, args)

    # remove the original files
    if args.large:
        for i in indv_outfiles:
            os.remove(i)

    # delete the master if it isn't the main operating file of the program
    if delete_master:
        os.remove(master_file)

    # master becomes a list of paths to the files in the temp_dir
    master_file = ['%s/%s' % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR)]

    # TODO: add multiprocessing here
    for tmpfile in master_file:
        if args.noclean:
            pass
        else:
            curr_id_dir = "%s/%s_outdir" % (dv.PROG_TEMP_DIR, Path(tmpfile).stem)
            # make a dir to store files exported based on qseqid
            Path(curr_id_dir).mkdir(parents=True, exist_ok=True)

            # flush alignments to query files and delete the master
            flush_outfile_to_tigs(tmpfile, curr_id_dir, "qseqid", False, tmpfile, args)

            qfiles = ['%s/%s' % (curr_id_dir, i) for i in os.listdir(curr_id_dir)]
            # open file containing a single contig's matches with a single query
            for qfile in qfiles:
                grimy = parse_blast_outfile(qfile, args)
                # if we are cleaning, clean up the grimy list; else, do no cleaning
                if args.noclean is False:
                    # pident is set to -1 here because the larger hit should be picked,
                    # not the one with the highest pident
                    grimy = clean_list_new(grimy, sstart, send, pident, False, args)

                # write the cleaned up file
                with open(qfile, 'w') as k:
                    filewriter = csv.writer(k, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                            lineterminator="\n")
                    for records in grimy:
                        filewriter.writerow(records)

            # combines all of the queryid files into an outfile for the contig
            with open(tmpfile, 'wb') as tmp:
                for out_file in qfiles:
                    with open(out_file, 'rb') as of:
                        copyfileobj(of, tmp)
            rmtree(curr_id_dir)

    # concatenates repeat data if noclean and nocat have not been called
    if args.nocat is False and args.noclean is False:
        if args.verbose:
            print("\nConcatenating repeat data....")

        for contig_files in master_file:
            grimy = parse_blast_outfile(contig_files, args)

            grimy = clean_list_new(grimy, sstart, send, pident, True, args)

            # write the cleaned up file for each contig
            with open(contig_files, 'w') as k:
                filewriter = csv.writer(k, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                        lineterminator="\n")
                for records in grimy:
                    filewriter.writerow(records)

    if args.verbose:
        print("\nSorting each contig's alignments by sstart....")

    multi_pool = mp.Pool(processes=dv.AVAIL_CORES)
    for outfile_name in ["%s/%s" % (dv.PROG_TEMP_DIR, i) for i in os.listdir(dv.PROG_TEMP_DIR)]:
        multi_pool.apply_async(sstart_sort_outfile, args=(outfile_name, args,))

    multi_pool.close()
    multi_pool.join()
    del multi_pool

    if args.verbose:
        print()

    # all hits are stored in dv.PROG_TEMP_DIR, just operate on that dir
    return dv.PROG_TEMP_DIR


# Runs blast with the given database and query and outputs data sorted by tig #
# if shortq == true, the command runs a shortblast instead of a normal blast  #
# shortblast has an evalue of 1e-5, normal blast has an evalue of 1e-15       #
# TODO: implement inverse (run one blast w -sub x -query y; one with -sub y query x); need to normalize output format
def execute_blast(db, query, outfilename,  inverse, args):
    # if running db, masker, or vcf module, use -db flag in BLASTn
    if args.subparser_id in dv.DB_ALIAS or args.subparser_id in dv.MASKER_ALIAS or args.subparser_id in dv.VC_ALIAS:
        blastcommand = ['blastn', '-db', db, "-query", query]
    # if running any other module, use -subject flag in BLASTn
    else:
        blastcommand = ['blastn', '-subject', db, '-query', query]

    # remove the dust filter if indicated
    try:
        if args.dust == "no":
            blastcommand.extend(["-dust", "no"])
    except AttributeError:
        pass

    # try to use a higher number of threads
    try:
        if args.num_threads > 0:
            blastcommand.extend(["-num_threads", str(args.num_threads)])
    except AttributeError:
        pass

    # generate the outfmt string to be used
    outfmt_string = generate_outfmt_str(args)

    # add in the evalue parameter
    blastcommand.extend(["-evalue", args.blaste])

    # allows user to run BLAST in different modes
    if args.subparser_id in dv.DB_ALIAS or args.subparser_id in dv.FREE_ALIAS:
        # Megablast is the default operation for blastn, so don't add anything if task is default, else add task
        if args.task != "megablast":
            blastcommand.extend(["-task", args.task])

    # add in outfmt string
    blastcommand.extend(["-outfmt", outfmt_string])

    # make blast pipe the results to a file in the temp directory
    blastcommand.extend(["-out", outfilename])

    if args.verbose:
        printcommand = blastcommand[:-2]  # removes the -out command from the printed command
        printcommand[len(printcommand) - 1] = "\"" + printcommand[len(printcommand) - 1] + "\""
        print("\nRunning BLAST with these options:\n%s\n" % (" ".join(printcommand)), end="")

    blast = subprocess.run(blastcommand, stdout=subprocess.PIPE, universal_newlines=True)

    if blast.returncode != 0:
        print("FATAL: BLASTn process did not have a 0 return code")
        exit()

    # returns the path to the outfile
    return outfilename


# TODO: this function probably should read outfiles line by line, but currently does not
def parse_blast_outfile(outfile_path, args):
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]
    pident = args.outfmt["pident"]

    if args.notigcov:
        # index the genome so you can access information about it
        genome_indx = SeqIO.to_dict(SeqIO.parse(str(Path(args.subject).resolve()), 'fasta'))

    with open(Path(outfile_path).resolve(), 'r') as k:
        grimy = k.readlines()
        # TODO: add case for qfile being empty, i don't think it will happen but jic
        # minus hits are changed to min,max; so sstrand needs to be present to preserve this information
        # this adds the sstrand qualifier, even if it was not present in the outfmt string
        if "sstrand" not in args.outfmt:
            args.outfmt.setdefault("sstrand", len(args.outfmt))

        # gets the length of the contig if you're removing hits that cover the contig
        if args.notigcov:
            curr_sseq_len = grimy[0].rstrip().split('\t')
            # if "FACET_" in curr_sseq_len[0]:
            #     curr_sseq_len.pop(0)
            curr_sseq_len = curr_sseq_len[args.outfmt["sseqid"]]
            curr_sseq_len = len(genome_indx[curr_sseq_len])

        # iterates through the grimy list
        aln = 0
        while aln < len(grimy):
            grimy[aln] = grimy[aln].rstrip().split('\t')
            # # removes FACET IDs from outfiles
            # if "FACET_" in grimy[aln][0]:
            #     grimy[aln].pop(0)

            # casts each part of the blast list to the defined type in the default values dictionary
            grimy[aln][sstart] = dv.BLASTN_TYPES["sstart"](grimy[aln][sstart])
            grimy[aln][send] = dv.BLASTN_TYPES["send"](grimy[aln][send])
            grimy[aln][pident] = dv.BLASTN_TYPES["pident"](grimy[aln][pident])

            # switches minus hits so they're [min,max]
            if grimy[aln][sstart] > grimy[aln][send]:
                # appends strand to outformats that don't already have the sstrand
                if len(grimy[aln]) == args.outfmt["sstrand"]:
                    grimy[aln].append('minus')
                grimy[aln][sstart], grimy[aln][send] = \
                    min(grimy[aln][send], grimy[aln][sstart]), \
                    max(grimy[aln][send], grimy[aln][sstart])
            else:
                # appends strand to outformats that don't already have the sstrand
                if len(grimy[aln]) == args.outfmt["sstrand"]:
                    grimy[aln].append('plus')

            if args.notigcov:
                if grimy[aln][sstart] < dv.TIG_SPAN_BUFFER and \
                        curr_sseq_len - dv.TIG_SPAN_BUFFER < grimy[aln][send]:
                    del grimy[aln]
                    aln -= 1
            aln += 1
    return grimy


def write_blast_outfile(grimy_list, outfile_name):
    with open(outfile_name, 'w') as k:
        filewriter = csv.writer(k, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                lineterminator="\n")
        for records in grimy_list:
            filewriter.writerow(records)


# sorts outfiles in outfile_dir by sstart
def sstart_sort_outfile(outfile_path, args):
    sstart = args.outfmt["sstart"]
    send = args.outfmt["send"]
    if args.noclean:
        if "sstrand" not in args.outfmt:
            args.outfmt.setdefault("sstrand", len(args.outfmt))

    if args.verbose:
        print("Sorting %s...." % outfile_path.split("_%stemp." % dv.PROG_NAME)[-1])
    # reads in the entire outfile #
    with open(outfile_path, 'r') as k:
        unsorted = k.readlines()
    # makes it into a list #
    for k in range(0, len(unsorted)):
        unsorted[k] = unsorted[k].rstrip().split('\t')
        unsorted[k][sstart] = dv.BLASTN_TYPES["sstart"](unsorted[k][sstart])
        unsorted[k][send] = dv.BLASTN_TYPES["send"](unsorted[k][send])
        # switches minus hits so they're [min,max]
        if unsorted[k][sstart] > unsorted[k][send]:
            # appends strand to outformats that don't already have the sstrand
            if len(unsorted[k]) == args.outfmt["sstrand"]:
                unsorted[k].append('minus')
            unsorted[k][sstart], unsorted[k][send] = min(unsorted[k][send], unsorted[k][sstart]),\
                                                        max(unsorted[k][send], unsorted[k][sstart])
        else:
            # appends strand to outformats that don't already have the sstrand
            if len(unsorted[k]) == args.outfmt["sstrand"]:
                unsorted[k].append('plus')
    # sorts it #
    unsorted = sorted(unsorted, key=lambda x: x[sstart])
    # writes it back out #
    with open(outfile_path, 'w', newline='') as k:
        filewriter = csv.writer(k, delimiter='\t', quotechar='\"', quoting=csv.QUOTE_MINIMAL,
                                lineterminator="\n")
        for records in unsorted:
            filewriter.writerow(records)
