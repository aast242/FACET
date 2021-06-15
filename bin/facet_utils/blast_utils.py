#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import subprocess
import os
from pathlib import Path

from Bio import SeqIO

from .defaults import ProgDefaults as dv
from . import btop_utils


# takes a list of lists and retains elements that have no complete overlap with other elements  #
# Input looks like: [[9,20], [2,3],[5,7],[2,8],[1,9],[1,10],[8,15]], returns [[1,10],[9,20]]    #
# index1 and index2 allow you to have lists that have more than two elements but you still want #
# to compare them as if they did. Input ranges must be [min, max]                               #
def clean_list(cluttered, index1, index2, pident_index, args):
    # sorts cluttered list by first part of each sublist
    cluttered = sorted(cluttered, key=lambda x: x[index1])
    # reduces function calls #
    clutteredremove = cluttered.remove
    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within other ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        # for loop checks to see if either list is containing the other #
        for i in range(1, len(cluttered)):
            if cluttered[i][index1] >= cluttered[i - 1][index1] - args.buffer and \
                    cluttered[i][index2] <= cluttered[i - 1][index2] + args.buffer:
                if pident_index >= 0:  # checks percent identity to retain matches that match the best
                    # checks to see if hits are approx. the same length #
                    if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer:
                        # removes the hit with lower pident #
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        else:
                            clutteredremove(cluttered[i])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break
                else:  # if we don't check by pident, remove the shorter element
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break
            elif cluttered[i - 1][index1] >= cluttered[i][index1] - args.buffer and \
                    cluttered[i - 1][index2] <= cluttered[i][index2] + args.buffer:
                if pident_index >= 0:  # checks percent identity to retain matches that match the best
                    # checks to see if hits are approx. the same length #
                    if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer:
                        # removes the hit with lower pident #
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        else:
                            clutteredremove(cluttered[i])
                    else:
                        clutteredremove(cluttered[i - 1])
                    exitcon = True
                    break
                else:  # if we don't check by pident, remove the shorter element
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break

    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within surrounding ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        for i in range(1, len(cluttered) - 1):
            # if the range in question is overlapping both the one before and after it
            if max(cluttered[i][index1], cluttered[i - 1][index1]) <= min(cluttered[i][index2],
                                                                          cluttered[i - 1][index2]):
                if max(cluttered[i][index1], cluttered[i + 1][index1]) <= min(cluttered[i][index2],
                                                                              cluttered[i + 1][index2]):
                    # if the range in question is contained by the range before and after
                    if cluttered[i - 1][index2] >= cluttered[i + 1][index1]:
                        exitcon = True
                        clutteredremove(cluttered[i])
                        break
    return cluttered


# takes the name of the new directory and the file used to create the database and makes a blast database #
# also returns the filepath of the newly created blast database                                           #
def make_blast_db(dbname, dbfile, facet_dir, args):
    outfilename = facet_dir + "/" + dbname + "/" + str(
        Path(Path(dbfile).resolve().name).with_suffix('')) + "_db.fasta"
    try:
        os.mkdir(facet_dir + "/" + dbname)
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


# Runs blast with the given database and query and outputs data sorted by tig #
# if shortq == true, the command runs a shortblast instead of a normal blast  #
# shortblast has an evalue of 1e-1, normal blast has an evalue of 1e-15       #
# TODO: implement inverse (run one blast w -sub x -query y; one with -sub y query x); need to normalize output format
def run_blast(db, query, shortq, inverse, args):
    # blast run contains only the necessary information
    if args.makedb:
        blastcommand = ['blastn', '-db', db, "-query", query]
    else:
        blastcommand = ['blastn', '-subject', db, '-query', query]
    try:
        if args.num_threads > 0:
            blastcommand.extend(["-num_threads", str(args.num_threads)])
    except AttributeError:
        pass
    outfmt_string = "6 sseqid sstart send qseqid qstart qend sstrand pident"
    if args.verboseoutput:  # or args.derip
        outfmt_string += " btop"
    """
    if args.verboseoutput:
        outfmt_string += " qseq"
    """
    if shortq:
        blastcommand.extend(["-evalue", args.shorte, "-task", "blastn-short",
                             "-outfmt", outfmt_string])
    else:
        blastcommand.extend(["-evalue", args.blaste,
                             "-outfmt", outfmt_string])
    blast = subprocess.run(blastcommand, stdout=subprocess.PIPE, universal_newlines=True)
    blastcommand[len(blastcommand) - 1] = "\"" + blastcommand[len(blastcommand) - 1] + "\""
    if args.verbose:
        print("BLAST was run with these options: ")
        print(" ".join(blastcommand))
    blast = blast.stdout.strip().split('\n')  # splits each line of the output
    blast = list(filter(None, blast))  # removes empty lines
    if args.verbose:
        print("BLAST found %s hits" % len(blast))
    return blast


# makes a blast dictionary from the blast list #
def build_blast_dict(blast, args):
    blast_dict = {}

    if len(blast) > 0:
        for i in range(0, len(blast)):  # splits each line into its components and casts values to int/float
            blast[i] = blast[i].split('\t')
            blast[i][1] = int(blast[i][1])
            blast[i][2] = int(blast[i][2])
            blast[i][4] = int(blast[i][4])
            blast[i][5] = int(blast[i][5])
            blast[i][7] = float(blast[i][7])
            if blast[i][1] > blast[i][2]:  # switches minus hits so that they're [less, greater]
                blast[i][1], blast[i][2] = min(blast[i][2], blast[i][1]), max(blast[i][2], blast[i][1])

        if args.verbose:
            print("Building BLAST dictionary.... ")

        # creates a dictionary of the BLAST output where the first level of keys are the repeat IDs and the second
        # level of keys are the tig IDs. Tig ID dicts contain a list of alignments for that element and that tig
        for align in blast:
            ele_key = align[3]
            blast_dict.setdefault(ele_key, {})
            tig_key = align[0]
            blast_dict[ele_key].setdefault(tig_key, []).append(align)

        # Sorts each Tig ID list by alignment start site on the tig (setup for filtering algorithm)
        for ele_id in blast_dict.keys():
            for tig_id in blast_dict[ele_id].keys():
                blast_dict[ele_id][tig_id] = sorted(blast_dict[ele_id][tig_id], key=lambda x: x[1])

    return blast_dict


# Checks rmDict against constDict and removes hits from rmlist that are covered by hits in constlist #
def compCombDict(constDict, rmDict, args):
    for tig in rmDict.keys():
        try:
            constDict[tig]
        except KeyError:
            # If the blast report doesn't have any hits in the tig, no need to compare #
            continue
        for rmranges in reversed(rmDict[tig]):
            for constranges in constDict[tig]:
                if rmranges[1] >= constranges[1] - args.buffer and rmranges[2] <= constranges[2] + args.buffer:
                    rmDict[tig].remove(rmranges)
                    break
    # removes empty tigs from the rmDict #
    for tigs in reversed(list(rmDict.keys())):
        if len(rmDict[tigs]) == 0:
            del rmDict[tigs]
    return rmDict


# Combines blast hits from different sources (i.e. BLASTn, and BLASTn-short). The combinedDictionary input    #
# is a dictionary that contains each of the blast hit dictionaries within it. The identifiers for each type   #
# are the keys used to call the dictionary. For a run with BLASTn and BLASTn-short, it would look like this:  #
# {blast: BLASTDICT, short: SHORTDICT} #
def combineData(combinedDictionary):
    tigKeysList = []
    combinedKeys = list(combinedDictionary.keys())
    # gets keys for every tig that had hits #
    for blasttype in combinedKeys:
        tigKeysList.extend(list(combinedDictionary[blasttype].keys()))
    tigKeysList = sorted(list(set(tigKeysList)))
    # initializes the final output dictionary and the list of zeroes to keep track of id numbers for each hit #
    finCombDict = {}
    idNums = [0] * len(combinedKeys)

    # iterates through each tig in the tigKeysList #
    for tigs in tigKeysList:
        templst = []
        # goes through each blasttype for each of the tigKeys #
        for blasttype in range(0, len(combinedKeys)):
            # tries to iterate through each of the hits in the tig #
            try:
                for hits in combinedDictionary[combinedKeys[blasttype]][tigs]:
                    # iterates the hitID number #
                    idNums[blasttype] += 1
                    # adds unique identifier to each hit #
                    hits.insert(0, combinedKeys[blasttype] + "Hit_" + str(idNums[blasttype]))
                # extends the templst with the modified tig #
                templst.extend(combinedDictionary[combinedKeys[blasttype]][tigs])
            # if there aren't any hits in the tig, passes #
            except KeyError:
                pass
        # adds a tig key with the combined hits list to the final dictionary #
        templst = sorted(templst, key=lambda x: x[2])
        finCombDict[tigs] = templst
    return finCombDict


# Fairly self-explanitory; goes in after cleaning and combines repeats into a single list by tig id #
def combineRepeatsByTig(dataList, args):
    if args.outfmt_six:
        sseqid = 1
        sstart = 8
    else:
        sseqid = dv.SSEQID_INDEX
        sstart = dv.SSTART_INDEX


    flattenedlist = []
    # appends every hit post-clean to a list #
    for repeats in dataList.keys():
        for tigs in dataList[repeats].keys():
            for hits in dataList[repeats][tigs]:
                flattenedlist.append(hits)
    # sorts 'flat' list by tig number #
    flattenedlist = sorted(flattenedlist, key=lambda x: x[sseqid])
    fin = []
    temp = []
    # breaks up the 'flat' list by tig #
    for i in range(0, len(flattenedlist)):
        if i != len(flattenedlist) - 1:
            if flattenedlist[i][sseqid] != flattenedlist[i + 1][sseqid]:
                temp.append(flattenedlist[i])
                fin.append(temp)
                temp = []
            else:
                temp.append(flattenedlist[i])
        else:
            if flattenedlist[i][sseqid] != flattenedlist[i - 1][sseqid]:
                fin.append([flattenedlist[i]])
            else:
                temp.append(flattenedlist[i])
                fin.append(temp)
    flattenedlist = fin
    # sorts tigs by start site and creates a dictionary from the tigs #
    tigDict = {}
    for tigs in range(0, len(flattenedlist)):
        flattenedlist[tigs] = sorted(flattenedlist[tigs], key=lambda x: x[sstart])
        tigDict[flattenedlist[tigs][0][sseqid]] = flattenedlist[tigs]
    return tigDict


# runs blastn and cleans up the output
def blastn_driver(dbfilepath, subj, query, args):
    blastdata = build_blast_dict(run_blast(dbfilepath, query, False, "IMPLEMENT VVQ", args), args)
    if len(blastdata) < 1:
        if not args.runshort or args.shortonly:
            print("FATAL: No BLAST hits")
            exit()
        else:
            print("No BLASTn hits, trying BLASTn-short....")
    else:
        if args.verbose:
            print("BLASTn completed!")
        if args.notigcov:
            if args.verbose:
                print("Removing hits covering the entire contig....")
            for seqrecord in SeqIO.parse(subj, "fasta"):
                try:
                    blastdata[seqrecord.id][seqrecord.id]
                except KeyError:
                    # If the blast report doesn't have any hits in the tig, no need to look #
                    continue
                for blasthits in blastdata[seqrecord.id][seqrecord.id]:
                    if int(blasthits[1]) < dv.TIG_SPAN_BUFFER and len(seqrecord.seq) - dv.TIG_SPAN_BUFFER < \
                            int(blasthits[2]) < len(seqrecord.seq) + dv.TIG_SPAN_BUFFER:
                        blastdata[seqrecord.id][seqrecord.id].remove(blasthits)
                    break
        if not args.noclean:
            if args.verbose:
                print("Cleaning BLASTn....")
            for repeats in blastdata.keys():
                for tigs in blastdata[repeats].keys():
                    # cleans each tig list in each repeat #
                    blastdata[repeats][tigs] = clean_list(blastdata[repeats][tigs], 1, 2, -1, args)
            if args.verbose:
                print("BLASTn cleaned!")
    return blastdata


# runs blastn-short and cleans up the output
def short_driver(dbfilepath, subj, query, blastn_results, args):
    shortblastdata = build_blast_dict(run_blast(dbfilepath, query, True, "IMPLEMENT VVQ", args), args)
    if len(shortblastdata) < 1 and len(blastn_results) < 1:
        print("FATAL: No BLAST hits")
        exit()
    elif len(shortblastdata) < 1 <= len(blastn_results):
        print("No BLASTn-short hits")
    else:
        if args.verbose:
            print("BLASTn-short completed!")
    if args.notigcov:
        if args.verbose:
            print("Removing hits covering the entire contig....")
        for seqrecord in SeqIO.parse(subj, "fasta"):
            try:
                shortblastdata[seqrecord.id][seqrecord.id]
            except KeyError:
                # If the blast report doesn't have any hits in the tig, no need to look #
                continue
            for blasthits in shortblastdata[seqrecord.id][seqrecord.id]:
                if int(blasthits[1]) < dv.TIG_SPAN_BUFFER and len(seqrecord.seq) - dv.TIG_SPAN_BUFFER < \
                        int(blasthits[2]) < len(seqrecord.seq) + dv.TIG_SPAN_BUFFER:
                    shortblastdata[seqrecord.id][seqrecord.id].remove(blasthits)
                    break
    if not args.noclean:
        if args.verbose:
            print("Cleaning BLASTn-short....")
        for repeats in shortblastdata.keys():
            for tigs in shortblastdata[repeats].keys():
                # cleans each tig list in each repeat #
                shortblastdata[repeats][tigs] = clean_list(shortblastdata[repeats][tigs], 1, 2, -1, args)
        if args.verbose:
            print("BLASTn-short cleaned!")
    return shortblastdata


# gets infile blast dictionary and cleans it up #
def infile_driver(blast_list, args):
    if args.outfmt_six and args.verboseoutput:
        print("FATAL: both \'--outfmt_six\' and \'--verboseoutput\' specifed!")
        print("       outfmt 6 does not contain btop information and cannot be used to make verbose outputs!")
        exit()
    verify_outfmt(blast_list, args)
    infile_data, aln_fmt_dict = infile_to_blast_dict(blast_list, args)
    if args.verbose:
        print("User infile verified!")
    if not args.noclean:
        if args.verbose:
            print("Cleaning user infile....")
        for repeats in infile_data.keys():
            for tigs in infile_data[repeats].keys():
                # cleans each tig list in each repeat #
                infile_data[repeats][tigs] = clean_list(infile_data[repeats][tigs], aln_fmt_dict["sstart_index"],
                                                        aln_fmt_dict["send_index"], -1, args)
        if args.verbose:
            print("User infile cleaned!")

    infile_data = combineRepeatsByTig(infile_data, args)
    if args.nocat is False and args.noclean is False:
        if args.verbose:
            print("\nConcatenating repeat data....")
        for tigs in infile_data.keys():
            infile_data[tigs] = clean_list(infile_data[tigs], aln_fmt_dict["sstart_index"],
                                           aln_fmt_dict["send_index"], aln_fmt_dict["pident_index"], args)
    infile_data = combineData({'blast': infile_data})
    blastdata = None
    del blastdata

    return infile_data


# combines results from blastn and blastn-short
def cat_megablast_short(mega, short, args):
    if args.runshort or args.shortonly:
        if args.nocat is False and args.noclean is False:
            if args.verbose:
                print("\nConcatenating repeat data....")
            # combines the repeats by tig for each dataset #
            shortblastdata = combineRepeatsByTig(short, args)
            blastdata = combineRepeatsByTig(mega, args)
            # removes repeats from shortblastdata overlapping with repeats from blastdata #
            shortblastdata = compCombDict(blastdata, shortblastdata, args)

            # blast list is cleaned by tig after being cleaned by repeat #
            for tigs in blastdata.keys():
                blastdata[tigs] = clean_list(blastdata[tigs], 1, 2, 7, args)

            # blastn-short list is cleaned by tig after being cleaned by repeat #
            for tigs in shortblastdata.keys():
                shortblastdata[tigs] = clean_list(shortblastdata[tigs], 1, 2, 7, args)

            # combines the data into a readable list #
            combinedhits = combineData({"blast": blastdata, "short": shortblastdata})
        else:
            if args.verbose:
                print("\nNot concatenating repeat data!")
            # combines the repeats by tig for each dataset #
            shortblastdata = combineRepeatsByTig(short, args)
            blastdata = combineRepeatsByTig(mega, args)

            # combines the data from each dataset #
            combinedhits = combineData({'blast': blastdata, 'short': shortblastdata})
        blastdata = None
        shortblastdata = None
        del blastdata
        del shortblastdata
        if args.verbose:
            print("All hits combined!")
    else:
        if args.verbose:
            print("\nNot running BLASTn-short!")
        blastdata = combineRepeatsByTig(mega, args)
        if args.nocat is False and args.noclean is False:
            if args.verbose:
                print("\nConcatenating repeat data....")
            for tigs in blastdata.keys():
                blastdata[tigs] = clean_list(blastdata[tigs], 1, 2, 7, args)
        combinedhits = combineData({'blast': blastdata})
        blastdata = None
        del blastdata
    return combinedhits


# makes sure the first entry of the blast_list is in the correct format #
# TODO: check subject and query files to make sure qseqid and sseqid are present in the files
def verify_outfmt(blast_list, args):
    try:
        check_list = blast_list[0].split('\t')
    except Exception:
        print("\nFATAL: infile could not be accessed or is empty")
        exit()

    if args.outfmt_six:
        if args.verbose:
            print("\nVerifying input is in outfmt 6....")
        try:
            str(check_list[0])
            str(check_list[1])
            float(check_list[2])
            int(check_list[3])
            int(check_list[4])
            int(check_list[5])
            int(check_list[6])
            int(check_list[7])
            int(check_list[8])
            int(check_list[9])
            float(check_list[10])
            int(check_list[11])
        except Exception:
            print("FATAL: first row of input does not match outfmt 6")
            exit()
    else:
        # print information about what the program is going to do #
        if args.verbose:
            print("\nNo outfmt specified for input file!")
            if not args.verboseoutput:
                print("Assuming default %s csv outfmt (\"%s\")" % (dv.PROG_NAME, dv.BASE_OUTFMT_STR))
            else:
                print("Assuming default %s outfmt (\"%s\")" % (dv.PROG_NAME, dv.VERBOSE_OUTFMT_STR))
                print("Verifying input has BTOP entries verbose output files....")

        if "blastHit_" in check_list[0] or "shortHit_" in check_list[0]:
            check_list.pop(0)

        # verify data types of entries 1-8 #
        try:
            str(check_list[0])
            int(check_list[1])
            int(check_list[2])
            str(check_list[3])
            int(check_list[4])
            int(check_list[5])
            str(check_list[6])
            float(check_list[7])
        except Exception:
            print("FATAL: provided file outfmt does not match %s's default outfmt" % dv.PROG_NAME)
            print("If your input file was generated using \'-outfmt 6\' in blastn, use the flag --outfmt_six")
            quit()

        # checks to see if BTOP can be converted to CIGAR and if the sequence is okay#
        if args.verboseoutput:
            try:
                btop_utils.btop_to_cigar(btop_utils.parse_btop(check_list[dv.BTOP_INDEX]))
            except Exception:
                print("FATAL: BTOP cannot be converted to CIGAR")
                quit()


# gets an input file to the point where it can be turned into a blast dictionary #
def infile_to_blast_dict(blast_list, args):
    # TODO: implement evalue filtering
    blast_dict = {}

    if args.outfmt_six:
        aln_fmt_dict = {
            "qseqid_index": 0,
            "sseqid_index": 1,
            "pident_index": 2,
            "sstart_index": 8,
            "send_index": 9
        }

        for i in range(0, len(blast_list)):  # splits each line into its components and casts values

            blast_list[i] = blast_list[i].split('\t')

            blast_list[i][0] = str(blast_list[i][0])  # qseqid
            blast_list[i][1] = str(blast_list[i][1])  # sseqid
            blast_list[i][2] = float(blast_list[i][2])  # pident
            blast_list[i][3] = int(blast_list[i][3])  # length
            blast_list[i][4] = int(blast_list[i][4])  # mismatch
            blast_list[i][5] = int(blast_list[i][5])  # gapopen
            blast_list[i][6] = int(blast_list[i][6])  # qstart
            blast_list[i][7] = int(blast_list[i][7])  # qend
            blast_list[i][8] = int(blast_list[i][8])  # sstart
            blast_list[i][9] = int(blast_list[i][9])  # send
            blast_list[i][10] = float(blast_list[i][10])  # evalue
            blast_list[i][11] = float(blast_list[i][11])  # bitscore

            # switches minus hits so they're [min,max]
            if blast_list[i][aln_fmt_dict["sstart_index"]] > blast_list[i][aln_fmt_dict["send_index"]]:
                blast_list[i].append('minus')
                blast_list[i][aln_fmt_dict["sstart_index"]], blast_list[i][aln_fmt_dict["send_index"]] = \
                    min(blast_list[i][aln_fmt_dict["send_index"]], blast_list[i][aln_fmt_dict["sstart_index"]]), \
                    max(blast_list[i][aln_fmt_dict["send_index"]], blast_list[i][aln_fmt_dict["sstart_index"]])
            else:
                blast_list[i].append('plus')

    else:
        aln_fmt_dict = {
            "qseqid_index": dv.QSEQID_INDEX,
            "sseqid_index": dv.SSEQID_INDEX,
            "pident_index": dv.PIDENT_INDEX,
            "sstart_index": dv.SSTART_INDEX,
            "send_index": dv.SEND_INDEX
        }
        for i in range(0, len(blast_list)):  # splits each line into its components and casts values

            blast_list[i] = blast_list[i].split('\t')

            # if ids from a previous FACET run are present, remove them
            if "blastHit_" in blast_list[i][0] or "shortHit_" in blast_list[i][0]:
                blast_list[i].pop(0)

            blast_list[i][0] = str(blast_list[i][0])  # sseqid
            blast_list[i][1] = int(blast_list[i][1])  # sstart
            blast_list[i][2] = int(blast_list[i][2])  # send
            blast_list[i][3] = str(blast_list[i][3])  # qseqid
            blast_list[i][4] = int(blast_list[i][4])  # qstart
            blast_list[i][5] = int(blast_list[i][5])  # qend
            blast_list[i][6] = str(blast_list[i][6])  # sstrand
            blast_list[i][7] = float(blast_list[i][7])  # pident

            if args.verboseoutput:
                blast_list[i][8] = str(blast_list[i][8])    # btop

            # switches minus hits so they're [min,max]
            if blast_list[i][aln_fmt_dict["sstart_index"]] > blast_list[i][aln_fmt_dict["send_index"]]:
                blast_list[i][aln_fmt_dict["sstart_index"]], blast_list[i][aln_fmt_dict["send_index"]] = \
                    min(blast_list[i][aln_fmt_dict["send_index"]], blast_list[i][aln_fmt_dict["sstart_index"]]), \
                    max(blast_list[i][aln_fmt_dict["send_index"]], blast_list[i][aln_fmt_dict["sstart_index"]])

    if args.verbose:
        print("Creating BLAST dictionary.... ")

    # creates a dictionary of the BLAST output where the first level of keys are the repeat IDs and the second
    # level of keys are the tig IDs. Tig ID dicts contain a list of alignments for that element and that tig
    for align in blast_list:
        ele_key = align[aln_fmt_dict["qseqid_index"]]
        blast_dict.setdefault(ele_key, {})
        tig_key = align[aln_fmt_dict["sseqid_index"]]
        blast_dict[ele_key].setdefault(tig_key, []).append(align)

    # Sorts each Tig ID list by alignment start site on the tig (setup for filtering algorithm)
    for ele_id in blast_dict.keys():
        for tig_id in blast_dict[ele_id].keys():
            blast_dict[ele_id][tig_id] = sorted(blast_dict[ele_id][tig_id],
                                                key=lambda x: x[aln_fmt_dict["sstart_index"]])

    return blast_dict, aln_fmt_dict


def masker_blast(db, query, args):
    if args.verbose:
        print("Running a self-blast on %s...." % Path(query).resolve())
    # blast run contains only the necessary information
    blastcommand = ['blastn', '-db', str(Path(db).resolve()), "-query", query]
    if args.num_threads > 0:
        blastcommand.extend(["-num_threads", str(args.num_threads)])

    outfmt_string = "6 sseqid sstart send qseqid qstart qend sstrand pident"

    blastcommand.extend(["-evalue", args.blaste, "-outfmt", outfmt_string])
    blast = subprocess.run(blastcommand, stdout=subprocess.PIPE, universal_newlines=True)
    blastcommand[len(blastcommand) - 1] = "\"" + blastcommand[len(blastcommand) - 1] + "\""
    if args.verbose:
        print("BLAST was run with these options: ")
        print(" ".join(blastcommand))
    blast = blast.stdout.strip().split('\n')  # splits each line of the output
    blast = list(filter(None, blast))  # removes empty lines
    if args.verbose:
        print("BLAST found %s hits" % len(blast))

    blast_dict = {}

    if len(blast) > 0:
        for i in range(0, len(blast)):  # splits each line into its components and casts values to int/float
            blast[i] = blast[i].split('\t')
            blast[i][1] = int(blast[i][1])
            blast[i][2] = int(blast[i][2])
            blast[i][4] = int(blast[i][4])
            blast[i][5] = int(blast[i][5])
            blast[i][7] = float(blast[i][7])
            if blast[i][1] > blast[i][2]:  # switches minus hits so that they're [less, greater]
                blast[i][1], blast[i][2] = min(blast[i][2], blast[i][1]), max(blast[i][2], blast[i][1])

        if args.verbose:
            print()
            print("Building BLAST dictionary.... ")

        for align in blast:
            tig_key = align[0]
            blast_dict.setdefault(tig_key, []).append(align)

        # Sorts each Tig ID list by alignment start site on the tig (setup for filtering algorithm)
        for tig_id in blast_dict.keys():
            blast_dict[tig_id] = sorted(blast_dict[tig_id], key=lambda x: x[1])

    return blast_dict